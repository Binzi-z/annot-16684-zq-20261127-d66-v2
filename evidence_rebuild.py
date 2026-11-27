#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
evidence_rebuild.py — 合成电子证据条的只读检查与重建工具（Python 3.12 标准库，Windows 可运行）。

场景：一批合成证据按“数据块 + 校验块”分条带保存。校验块由输入清单给出的
有限域线性关系生成（每条方程  sum_i c_i * block_i = 0，系数与块取值均在同一
有限域内）。抽检发现：有些块缺失，有些块虽存在但与预存 SHA-256 摘要不明珠不符。

本工具：
  * 只读检查：输入块只以 'rb' 打开，从不改写、从不删除；
  * 每条带判定：哪些缺块可唯一恢复、哪些存在多解、哪些现有块相互矛盾；
  * 输出可独立复核的证据：矩阵秩 / 零化度 / 每块实际与预期摘要 / 见证摘要；
  * 多解时给出两个不同的合法见证（均通过方程复核）；
  * 小条带（解数 <= --enumerate-max）穷举每个符号位置上的解集合；
  * 较大输入下按“擦除模式”缓存行最简形，统计实际命中次数与峰值条目数；
  * rebuild 子命令把唯一可恢复的块写入一个“新目录”，原目录绝不写入。

清单（manifest）JSON 格式：
{
  "format": "evidence-pack/1",
  "field":  {"kind": "gf256", "poly": 283, "generator": 3}   // 或 {"kind":"prime","p":251}
  "stripes": [
    {
      "id": "s0",
      "block_size": 16,                      // 可选；缺省时取摘要有效块的多数长度
      "blocks": ["s0/d0.bin", "s0/d1.bin", "s0/p0.bin"],   // 相对清单所在目录
      "equations": [[1,1,1], [1,2,1]],        // m 行 × n 列，sum c_i*block_i = 0
      "sha256": {"s0/d0.bin": "<hex>", ...}   // 预存摘要
    }
  ]
}

用法：
  python evidence_rebuild.py check   MANIFEST [--report R.json] [--enumerate-max N] [--cache-size N]
  python evidence_rebuild.py rebuild MANIFEST --out OUTDIR [--report R.json] [--no-copy-ok]
  python evidence_rebuild.py verify  MANIFEST --report R.json     # 独立复核报告内容
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from collections import OrderedDict
from datetime import datetime, timezone

DEFAULT_POLY = 0x11B   # AES 不可约多项式 x^8+x^4+x^3+x+1
DEFAULT_GEN = 3        # 0x11B 下的本原元


# ---------------------------------------------------------------- 有限域

class GF256:
    """GF(2^8)，符号即字节。"""
    kind = "gf256"

    def __init__(self, poly: int = DEFAULT_POLY, gen: int = DEFAULT_GEN):
        self.poly = poly
        self.gen = gen
        self.q = 256
        exp = [0] * 512
        log = [0] * 256
        x = 1
        for i in range(255):
            exp[i] = x
            log[x] = i
            x = self._mul_raw(x, gen)
        if len(set(exp[:255])) != 255 or x != 1:
            raise ValueError("generator %d is not primitive for poly 0x%03x" % (gen, poly))
        for i in range(255, 512):
            exp[i] = exp[i - 255]
        self.exp = exp
        self.log = log

    def _mul_raw(self, a: int, b: int) -> int:
        r = 0
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
            if a & 0x100:
                a ^= self.poly
        return r & 0xFF

    def add(self, a, b): return a ^ b
    def sub(self, a, b): return a ^ b
    def neg(self, a): return a

    def mul(self, a, b):
        if a == 0 or b == 0:
            return 0
        return self.exp[self.log[a] + self.log[b]]

    def inv(self, a):
        if a == 0:
            raise ZeroDivisionError("gf256: inv(0)")
        return self.exp[255 - self.log[a]]

    def describe(self):
        return {"kind": "gf256", "poly": self.poly, "generator": self.gen}


class GFp:
    """素域 GF(p)，p <= 256，符号以单字节存储且必须 < p。"""
    kind = "prime"

    def __init__(self, p: int):
        if p < 2 or p > 256 or not _is_prime(p):
            raise ValueError("p must be a prime in [2, 256], got %r" % (p,))
        self.p = p
        self.q = p

    def add(self, a, b): return (a + b) % self.p
    def sub(self, a, b): return (a - b) % self.p
    def neg(self, a): return (-a) % self.p
    def mul(self, a, b): return (a * b) % self.p

    def inv(self, a):
        if a % self.p == 0:
            raise ZeroDivisionError("GF(%d): inv(0)" % self.p)
        return pow(a, self.p - 2, self.p)

    def describe(self):
        return {"kind": "prime", "p": self.p}


def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def build_field(spec: dict):
    kind = (spec or {}).get("kind", "gf256")
    if kind == "gf256":
        return GF256(spec.get("poly", DEFAULT_POLY), spec.get("generator", DEFAULT_GEN))
    if kind == "prime":
        return GFp(int(spec["p"]))
    raise ValueError("unknown field kind: %r" % (kind,))


# ---------------------------------------------------------------- 矩阵（行最简形 + 变换矩阵）

def rref_aug(field, A):
    """对 m×n 矩阵 A 做 Gauss-Jordan，返回 (R, T, pivots)，满足 T·A = R（R 为行最简形）。"""
    m = len(A)
    n = len(A[0]) if m else 0
    M = [list(row) + [1 if i == j else 0 for j in range(m)] for i, row in enumerate(A)]
    pivots = []
    row = 0
    for col in range(n):
        piv = None
        for r in range(row, m):
            if M[r][col] != 0:
                piv = r
                break
        if piv is None:
            continue
        M[row], M[piv] = M[piv], M[row]
        inv = field.inv(M[row][col])
        M[row] = [field.mul(v, inv) for v in M[row]]
        for r in range(m):
            if r != row and M[r][col] != 0:
                f = M[r][col]
                M[r] = [field.sub(a, field.mul(f, b)) for a, b in zip(M[r], M[row])]
        pivots.append(col)
        row += 1
        if row == m:
            break
    R = [r[:n] for r in M]
    T = [r[n:] for r in M]
    return R, T, pivots


class ReducedSystem:
    """未知列子矩阵 A（m×u）的行最简形、零空间基等可复用结果。"""

    def __init__(self, field, A, unknown_cols):
        self.unknown_cols = list(unknown_cols)
        self.u = len(unknown_cols)
        self.m = len(A)
        self.R, self.T, self.pivots = rref_aug(field, A)
        self.rank = len(self.pivots)
        self.free = [c for c in range(self.u) if c not in self.pivots]
        self.nullity = len(self.free)
        self.null_basis = []
        for f in self.free:
            v = [0] * self.u
            v[f] = 1
            for i, pc in enumerate(self.pivots):
                v[pc] = field.neg(self.R[i][f])
            self.null_basis.append(v)


class SystemCache:
    """按擦除模式（方程指纹 + 未知列集合）缓存 ReducedSystem，统计命中与峰值。"""

    def __init__(self, maxsize=None):
        self.maxsize = maxsize
        self.store = OrderedDict()
        self.hits = 0
        self.misses = 0
        self.peak = 0
        self.uses = {}

    def get_or_build(self, key, build):
        self.uses[key] = self.uses.get(key, 0) + 1
        if key in self.store:
            self.hits += 1
            self.store.move_to_end(key)
            return self.store[key]
        self.misses += 1
        val = build()
        self.store[key] = val
        if self.maxsize is not None:
            while len(self.store) > self.maxsize:
                self.store.popitem(last=False)
        if len(self.store) > self.peak:
            self.peak = len(self.store)
        return val

    def stats(self):
        return {
            "lookups": self.hits + self.misses,
            "hits": self.hits,
            "misses": self.misses,
            "peak_entries": self.peak,
            "max_size": self.maxsize,
            "patterns": [
                {
                    "equations_fingerprint": fp[:12],
                    "unknown_columns": list(cols),
                    "uses": uses,
                }
                for (fp, cols), uses in sorted(self.uses.items(), key=lambda kv: (kv[0][0], kv[0][1]))
            ],
        }


# ---------------------------------------------------------------- 工具函数

def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def safe_join(base: str, rel: str):
    p = os.path.abspath(os.path.join(base, rel))
    try:
        if os.path.commonpath([p, base]) != base:
            return None
    except ValueError:
        return None
    return p


def erasure_runs(cols):
    """把未知列下标合并成连续区间（用于刻画突发连续损坏）。"""
    out = []
    for c in cols:
        if out and c == out[-1][1] + 1:
            out[-1][1] = c
        else:
            out.append([c, c])
    return out


def assignment_ok(field, eqs, full_blocks) -> bool:
    """复核一组完整块赋值是否满足全部线性方程。"""
    L = len(full_blocks[0])
    for eq in eqs:
        acc = bytearray(L)
        for c, data in zip(eq, full_blocks):
            if c == 0:
                continue
            if c == 1:
                for pos in range(L):
                    acc[pos] = field.add(acc[pos], data[pos])
            else:
                for pos in range(L):
                    acc[pos] = field.add(acc[pos], field.mul(c, data[pos]))
        if any(acc):
            return False
    return True


def enumerate_values(field, part, null_basis, L):
    """穷举每个未知列、每个符号位置上的所有可能取值（仅在小条带上调用）。"""
    u = len(part)
    d = len(null_basis)
    per_pos = [[set() for _ in range(L)] for _ in range(u)]
    combo = [0] * d
    for _ in range(field.q ** d):
        for col in range(u):
            delta = 0
            for i in range(d):
                if combo[i]:
                    delta = field.add(delta, field.mul(combo[i], null_basis[i][col]))
            prow = part[col]
            sets = per_pos[col]
            for pos in range(L):
                sets[pos].add(field.add(prow[pos], delta))
        for i in range(d):  # 进位计数
            combo[i] += 1
            if combo[i] == field.q:
                combo[i] = 0
            else:
                break
    return per_pos


# ---------------------------------------------------------------- 条带处理

def process_stripe(field, stripe, base, cache, enumerate_max):
    """返回 (结果 dict, 可重建写入列表 [(rel, bytes)])。结果 dict 全部为 JSON 可序列化。"""
    blocks = stripe["blocks"]
    n = len(blocks)
    eqs = stripe.get("equations", [])
    digests = stripe.get("sha256", {})
    for eq in eqs:
        if len(eq) != n:
            raise ValueError("stripe %s: equation length %d != block count %d"
                             % (stripe.get("id"), len(eq), n))
        for c in eq:
            if not (0 <= c < field.q):
                raise ValueError("stripe %s: coefficient %d out of field" % (stripe.get("id"), c))
    m = len(eqs)
    eq_fp = hashlib.sha256(json.dumps(eqs, separators=(",", ":")).encode()).hexdigest()

    # 1) 只读读取 + 摘要核对
    infos = []
    for rel in blocks:
        info = {"path": rel, "sha256_expected": digests.get(rel)}
        apath = safe_join(base, rel)
        if apath is None:
            info["state"] = "bad-path"
            info["data"] = None
        elif not os.path.isfile(apath):
            info["state"] = "missing"
            info["data"] = None
        else:
            with open(apath, "rb") as f:      # 只读，绝不改写
                data = f.read()
            actual = sha256_hex(data)
            info["sha256_actual"] = actual
            if info["sha256_expected"] is not None and actual != info["sha256_expected"]:
                info["state"] = "corrupt"     # 与预存摘要不符：按擦除处理
                info["data"] = None
            elif field.kind == "prime" and any(b >= field.p for b in data):
                info["state"] = "symbol-out-of-range"
                info["data"] = None
            else:
                info["state"] = "ok"
                info["data"] = data
        infos.append(info)

    # 2) 长度一致性（短块同样适用）
    expected = stripe.get("block_size")
    ok_lens = [len(i["data"]) for i in infos if i["state"] == "ok"]
    if expected is None and ok_lens:
        expected = max(set(ok_lens), key=ok_lens.count)
    for i in infos:
        if i["state"] == "ok" and expected is not None and len(i["data"]) != expected:
            i["state"] = "length-mismatch"
            i["data"] = None

    unknown = [idx for idx, i in enumerate(infos) if i["state"] != "ok"]
    known = [idx for idx, i in enumerate(infos) if i["state"] == "ok"]

    # 3) 未知列子矩阵 —— 走擦除模式缓存
    A = [[eq[j] for j in unknown] for eq in eqs]
    key = (eq_fp, tuple(unknown))
    rsys = cache.get_or_build(key, lambda: ReducedSystem(field, A, unknown))

    res = {
        "id": stripe.get("id"),
        "blocks": [{k: v for k, v in i.items() if k != "data"} for i in infos],
        "unknown_blocks": [blocks[j] for j in unknown],
        "erasure_runs": erasure_runs(unknown),
        "matrix": {
            "equations": m,
            "blocks": n,
            "unknowns": len(unknown),
            "rank": rsys.rank,
            "nullity": rsys.nullity,
            "pivot_columns": rsys.pivots,
            "free_columns": rsys.free,
        },
    }

    L = expected
    if L is None:
        res["status"] = "indeterminate"
        res["note"] = "no digest-valid block and no block_size; rank evidence only"
        return res, []

    # 4) 右端项 b = -(已知列贡献)，再施加行变换 T
    b = [[0] * L for _ in range(m)]
    for r, eq in enumerate(eqs):
        brow = b[r]
        for j in known:
            c = eq[j]
            if c == 0:
                continue
            data = infos[j]["data"]
            if c == 1:
                for pos in range(L):
                    brow[pos] = field.add(brow[pos], data[pos])
            else:
                for pos in range(L):
                    brow[pos] = field.add(brow[pos], field.mul(c, data[pos]))
    for r in range(m):
        brow = b[r]
        for pos in range(L):
            brow[pos] = field.neg(brow[pos])

    bp = [[0] * L for _ in range(m)]
    for i in range(m):
        Ti = rsys.T[i]
        out = bp[i]
        for j, t in enumerate(Ti):
            if t == 0:
                continue
            src = b[j]
            if t == 1:
                for pos in range(L):
                    out[pos] = field.add(out[pos], src[pos])
            else:
                for pos in range(L):
                    out[pos] = field.add(out[pos], field.mul(t, src[pos]))

    # 5) 相容性：秩之外的行若右端非零，则现有块相互矛盾
    bad_pos = [pos for pos in range(L)
               if any(bp[i][pos] != 0 for i in range(rsys.rank, m))]
    res["matrix"]["consistent"] = not bad_pos

    if bad_pos:
        involved_eqs = sorted({
            j
            for i in range(rsys.rank, m)
            if any(bp[i][pos] for pos in bad_pos)
            for j, t in enumerate(rsys.T[i]) if t != 0
        })
        involved_blocks = sorted({
            blocks[j]
            for r in involved_eqs
            for j, c in enumerate(eqs[r])
            if c != 0 and j in known
        })
        res["status"] = "contradictory"
        res["contradiction"] = {
            "inconsistent_position_count": len(bad_pos),
            "first_positions": bad_pos[:16],
            "equations_involved": involved_eqs,
            "known_blocks_involved": involved_blocks,
        }
        return res, []

    # 6) 特解（自由变量取 0）
    part = [[0] * L for _ in unknown]
    for i, pc in enumerate(rsys.pivots):
        part[pc] = list(bp[i])

    def full_with(assign):
        full = [None] * n
        for j in known:
            full[j] = infos[j]["data"]
        for col, j in enumerate(unknown):
            full[j] = assign[col]
        return full

    if not unknown:
        res["status"] = "ok"
        return res, []

    if rsys.nullity == 0:
        # 唯一可恢复
        res["status"] = "unique"
        recovered = []
        writes = []
        for col, j in enumerate(unknown):
            data = bytes(part[col])
            digest = sha256_hex(data)
            recovered.append({
                "path": blocks[j],
                "was": infos[j]["state"],
                "sha256": digest,
                "matches_manifest": digest == digests.get(blocks[j]),
            })
            writes.append((blocks[j], data))
        res["recovered"] = recovered
        res["assignment_verified"] = assignment_ok(
            field, eqs, full_with([bytes(part[c]) for c in range(len(unknown))]))
        return res, writes

    # 7) 多解：给两个不同的合法见证
    res["status"] = "multiple"
    d = rsys.nullity
    count = field.q ** d
    res["solution_space"] = {"dimension": d, "solutions_per_position": count}

    witnesses = []
    for combo in ([0] * d, [1] + [0] * (d - 1)):
        assign = []
        for col in range(len(unknown)):
            buf = bytearray(part[col])
            delta = 0
            for i, cf in enumerate(combo):
                if cf:
                    delta = field.add(delta, field.mul(cf, rsys.null_basis[i][col]))
            if delta:
                for pos in range(L):
                    buf[pos] = field.add(buf[pos], delta)
            assign.append(bytes(buf))
        witnesses.append(assign)

    wrep = []
    for widx, assign in enumerate(witnesses, 1):
        verified = assignment_ok(field, eqs, full_with(assign))
        wrep.append({
            "index": widx,
            "verified": verified,
            "blocks": {
                blocks[unknown[c]]: {
                    "sha256": sha256_hex(assign[c]),
                    **({"hex": assign[c].hex()} if L <= 64 else {}),
                }
                for c in range(len(unknown))
            },
        })
    res["witnesses"] = wrep
    res["witnesses_differ"] = witnesses[0] != witnesses[1]

    # 8) 小条带穷举解集合
    if count <= enumerate_max:
        per_pos = enumerate_values(field, part, rsys.null_basis, L)
        pv = {}
        for col, j in enumerate(unknown):
            sets = per_pos[col]
            if L <= 16:
                pv[blocks[j]] = {
                    "per_position": ["".join("%02x" % v for v in sorted(s)) for s in sets]
                }
            else:
                union = sorted(set().union(*sets))
                pv[blocks[j]] = {"union": "".join("%02x" % v for v in union)}
        res["enumeration"] = {
            "solutions_per_position": count,
            "possible_values": pv,
        }
    else:
        res["enumeration"] = {
            "solutions_per_position": count,
            "skipped": "exceeds --enumerate-max (%d)" % enumerate_max,
        }
    return res, []


# ---------------------------------------------------------------- 运行 / 报告

def load_manifest(path):
    ap = os.path.abspath(path)
    with open(ap, "rb") as f:
        raw = f.read()
    manifest = json.loads(raw.decode("utf-8"))
    return manifest, raw, os.path.dirname(ap)


def run_all(manifest, base, enumerate_max, cache_size):
    field = build_field(manifest.get("field"))
    cache = SystemCache(cache_size or None)
    results = []
    writes = []
    for stripe in manifest["stripes"]:
        res, w = process_stripe(field, stripe, base, cache, enumerate_max)
        results.append(res)
        writes.extend(w)
    return field, results, writes, cache


def build_report(manifest_path, manifest_raw, base, field, results, cache, rebuild_info):
    return {
        "tool": "evidence-rebuild/1.0",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "manifest_path": manifest_path,
        "manifest_sha256": sha256_hex(manifest_raw),
        "input_dir": base,
        "field": field.describe(),
        "stripes": results,
        "cache": cache.stats(),
        "rebuild": rebuild_info,
    }


def print_summary(results, cache):
    for r in results:
        m = r["matrix"]
        print("stripe %-18s status=%-14s unknowns=%d rank=%d nullity=%d runs=%s"
              % (r["id"], r["status"], m["unknowns"], m["rank"], m["nullity"],
                 r.get("erasure_runs")))
    s = cache.stats()
    print("cache: lookups=%d hits=%d misses=%d peak_entries=%d"
          % (s["lookups"], s["hits"], s["misses"], s["peak_entries"]))


def dirs_overlap(a, b):
    try:
        return os.path.commonpath([a, b]) in (a, b)
    except ValueError:
        return False


# ---------------------------------------------------------------- 子命令

def cmd_check(args):
    manifest, raw, base = load_manifest(args.manifest)
    field, results, _writes, cache = run_all(manifest, base, args.enumerate_max, args.cache_size)
    print_summary(results, cache)
    if args.report:
        report = build_report(os.path.abspath(args.manifest), raw, base, field, results, cache, None)
        with open(args.report, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=1, ensure_ascii=False)
        print("report written: %s" % args.report)
    return 0


def cmd_rebuild(args):
    manifest, raw, base = load_manifest(args.manifest)
    out_dir = os.path.abspath(args.out)
    if dirs_overlap(out_dir, base):
        print("error: --out directory must be a NEW directory outside the input tree", file=sys.stderr)
        return 2
    field, results, writes, cache = run_all(manifest, base, args.enumerate_max, args.cache_size)

    written = []
    skipped = []
    # 唯一可恢复的缺块/腐坏块
    for rel, data in writes:
        target = safe_join(out_dir, rel)
        os.makedirs(os.path.dirname(target), exist_ok=True)
        with open(target, "wb") as f:
            f.write(data)
        written.append(rel)
    # 默认把摘要有效的现存块一并复制，使新目录成为完整可用副本
    if not args.no_copy_ok:
        for res in results:
            for binfo in res["blocks"]:
                if binfo["state"] != "ok":
                    continue
                rel = binfo["path"]
                if rel in written:
                    continue
                src = safe_join(base, rel)
                with open(src, "rb") as f:
                    data = f.read()
                target = safe_join(out_dir, rel)
                os.makedirs(os.path.dirname(target), exist_ok=True)
                with open(target, "wb") as f:
                    f.write(data)
                written.append(rel)
    # 多解 / 矛盾条带中的未知块不写入
    for res in results:
        if res["status"] in ("multiple", "contradictory", "indeterminate"):
            skipped.extend(res["unknown_blocks"])

    print_summary(results, cache)
    print("rebuild: %d blocks written to %s ; %d ambiguous blocks skipped"
          % (len(written), out_dir, len(skipped)))
    if args.report:
        report = build_report(os.path.abspath(args.manifest), raw, base, field, results, cache,
                              {"out_dir": out_dir, "written": written, "skipped": skipped})
        with open(args.report, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=1, ensure_ascii=False)
        print("report written: %s" % args.report)
    return 0


COMPARE_KEYS = ["status", "matrix", "blocks", "unknown_blocks", "erasure_runs",
                "recovered", "assignment_verified", "witnesses", "witnesses_differ",
                "solution_space", "enumeration", "contradiction"]


def cmd_verify(args):
    """独立复核：从清单与原始输入重算全部证据，与报告逐项比对。"""
    with open(args.report, "r", encoding="utf-8") as f:
        rep = json.load(f)
    manifest, raw, base = load_manifest(args.manifest)
    problems = []
    if sha256_hex(raw) != rep.get("manifest_sha256"):
        problems.append("manifest sha256 mismatch with report")
    field, results, _w, cache = run_all(manifest, base, 256, 0)
    old = {s["id"]: s for s in rep.get("stripes", [])}
    for res in results:
        o = old.get(res["id"])
        if o is None:
            problems.append("stripe %s missing from report" % res["id"])
            continue
        for k in COMPARE_KEYS:
            if res.get(k) != o.get(k):
                problems.append("stripe %s: field %r differs" % (res["id"], k))
    cs = cache.stats()
    for k in ("lookups", "hits", "misses", "peak_entries"):
        if rep.get("cache", {}).get(k) != cs[k]:
            problems.append("cache stat %r differs: report=%r recomputed=%r"
                            % (k, rep.get("cache", {}).get(k), cs[k]))
    if problems:
        for p in problems:
            print("verify FAIL: %s" % p)
        return 1
    print("verify OK: %d stripes re-checked; ranks, digests, witnesses and cache stats all match"
          % len(results))
    return 0


def _build_parser():
    ap = argparse.ArgumentParser(
        prog="evidence_rebuild",
        description="read-only check & rebuild for finite-field striped evidence packs")
    sub = ap.add_subparsers(dest="cmd", required=True)

    def add_common(p):
        p.add_argument("manifest")
        p.add_argument("--report")
        p.add_argument("--enumerate-max", type=int, default=256,
                       help="enumerate solution set when solutions-per-position <= N (default 256)")
        p.add_argument("--cache-size", type=int, default=0,
                       help="LRU cap on erasure-pattern cache; 0 = unbounded (default)")

    add_common(sub.add_parser("check", help="read-only inspection"))
    pr = sub.add_parser("rebuild", help="rebuild uniquely-recoverable blocks into a NEW directory")
    add_common(pr)
    pr.add_argument("--out", required=True, help="new output directory (must be outside input tree)")
    pr.add_argument("--no-copy-ok", action="store_true",
                    help="do not copy digest-valid blocks into the output directory")
    pv = sub.add_parser("verify", help="independently recompute and compare against a report")
    pv.add_argument("manifest")
    pv.add_argument("--report", required=True)
    return ap


def cli(argv=None):
    ap = _build_parser()
    args = ap.parse_args(argv)
    if args.cmd == "check":
        return cmd_check(args)
    if args.cmd == "rebuild":
        return cmd_rebuild(args)
    if args.cmd == "verify":
        return cmd_verify(args)
    return 2


if __name__ == "__main__":
    sys.exit(cli())
