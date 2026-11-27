#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_sample.py — 生成一个“受损”的合成证据包并运行 evidence_rebuild.py 做端到端自检。

覆盖场景：
  s-clean            全部块完好
  s-unique-missing   缺 1 个数据块            -> 唯一恢复
  s-unique-corrupt   1 个数据块与预存摘要不符 -> 按擦除处理，唯一恢复且摘要吻合
  s-unique-parity    缺 1 个校验块            -> 唯一恢复
  s-multi            缺 3 块（2 条方程）      -> 多解，给两个不同见证 + 穷举（短块 L=4）
  s-contra           摘要有效但违反线性关系   -> 现有块相互矛盾
  s-burst-00..23     较大输入：短块(L=8)、连续突发损坏（两种擦除模式各重复 12 次）
                     -> 同一擦除模式命中缓存，统计实际命中次数与峰值
"""
import hashlib
import json
import os
import random
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SAMPLE = os.path.join(HERE, "sample")
OUT = os.path.join(HERE, "rebuilt")
REPORT = os.path.join(HERE, "report.json")

sys.path.insert(0, HERE)
from evidence_rebuild import GF256  # noqa: E402

F = GF256()


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def make_pack():
    if os.path.exists(SAMPLE):
        shutil.rmtree(SAMPLE)
    os.makedirs(SAMPLE)
    rng = random.Random(20261128)
    stripes = []

    def add_stripe(sid, k, eqs, L, damage=None, manifest_tamper=None):
        m = len(eqs)
        n = k + m
        names = (["%s/d%d.bin" % (sid, i) for i in range(k)]
                 + ["%s/p%d.bin" % (sid, j) for j in range(m)])
        datas = [rng.randbytes(L) for _ in range(k)]
        # 由方程求校验块：第 r 条方程中校验列 k+r 系数为 1，其余校验列为 0
        par = []
        for r, eq in enumerate(eqs):
            assert eq[k + r] == 1 and all(eq[k + j] == 0 for j in range(m) if j != r)
            acc = bytearray(L)
            for i in range(k):
                c = eq[i]
                if c:
                    for pos in range(L):
                        acc[pos] = F.add(acc[pos], F.mul(c, datas[i][pos]))
            par.append(bytes(F.neg(v) for v in acc))
        allb = list(datas) + par
        if manifest_tamper:
            # 篡改块内容，但清单记录的是篡改后的摘要 -> 摘要有效却违反线性关系
            allb = manifest_tamper(allb)
        os.makedirs(os.path.join(SAMPLE, sid))
        for rel, data in zip(names, allb):
            with open(os.path.join(SAMPLE, rel), "wb") as f:
                f.write(data)
        digests = {rel: sha(data) for rel, data in zip(names, allb)}
        stripes.append({
            "id": sid,
            "block_size": L,
            "blocks": names,
            "equations": eqs,
            "sha256": digests,
        })
        if damage:
            damage(names)

    def delete(*idx):
        def dmg(names):
            for i in idx:
                os.remove(os.path.join(SAMPLE, names[i]))
        return dmg

    def corrupt_first_byte(idx):
        def dmg(names):
            p = os.path.join(SAMPLE, names[idx])
            with open(p, "rb") as f:
                b = bytearray(f.read())
            b[0] ^= 0xFF
            with open(p, "wb") as f:
                f.write(b)
        return dmg

    H5 = [[1, 1, 1, 1, 0],
          [1, 2, 3, 0, 1]]
    add_stripe("s-clean", 3, H5, 16)
    add_stripe("s-unique-missing", 3, H5, 16, damage=delete(1))
    add_stripe("s-unique-corrupt", 3, H5, 16, damage=corrupt_first_byte(0))
    add_stripe("s-unique-parity", 3, H5, 16, damage=delete(4))
    add_stripe("s-multi", 3, H5, 4, damage=delete(0, 1, 3))

    def tamper(allb):
        b = bytearray(allb[2])
        b[0] ^= 0x01
        allb[2] = bytes(b)
        return allb
    add_stripe("s-contra", 3, H5, 16, manifest_tamper=tamper)

    HB = [[1, 1, 1, 1, 1, 1, 1, 0],
          [1, 2, 3, 4, 5, 6, 0, 1]]
    for t in range(24):
        drop = (2, 3) if t < 12 else (4, 5)   # 两种连续突发损坏模式
        add_stripe("s-burst-%02d" % t, 6, HB, 8, damage=delete(*drop))

    manifest = {
        "format": "evidence-pack/1",
        "field": {"kind": "gf256", "poly": 0x11B, "generator": 3},
        "stripes": stripes,
    }
    mpath = os.path.join(SAMPLE, "manifest.json")
    with open(mpath, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=1)
    return mpath


def snapshot_input():
    """记录受损后输入目录中每个文件的摘要，用于事后证明工具未改写原块。"""
    snap = {}
    for root, _dirs, files in os.walk(SAMPLE):
        for fn in files:
            p = os.path.join(root, fn)
            with open(p, "rb") as f:
                snap[os.path.relpath(p, SAMPLE)] = sha(f.read())
    return snap


def run_tool(mpath):
    if os.path.exists(OUT):
        shutil.rmtree(OUT)
    tool = os.path.join(HERE, "evidence_rebuild.py")
    subprocess.run([sys.executable, tool, "rebuild", mpath,
                    "--out", OUT, "--report", REPORT], check=True)
    print("-" * 70)
    subprocess.run([sys.executable, tool, "verify", mpath,
                    "--report", REPORT], check=True)


def check_expectations(mpath, snap):
    with open(REPORT, "r", encoding="utf-8") as f:
        rep = json.load(f)
    by_id = {s["id"]: s for s in rep["stripes"]}
    fails = []

    def expect(cond, msg):
        if not cond:
            fails.append(msg)

    expect(by_id["s-clean"]["status"] == "ok", "s-clean should be ok")
    expect(by_id["s-unique-missing"]["status"] == "unique", "s-unique-missing should be unique")
    expect(by_id["s-unique-corrupt"]["status"] == "unique", "s-unique-corrupt should be unique")
    expect(by_id["s-unique-parity"]["status"] == "unique", "s-unique-parity should be unique")
    expect(by_id["s-multi"]["status"] == "multiple", "s-multi should be multiple")
    expect(by_id["s-contra"]["status"] == "contradictory", "s-contra should be contradictory")
    for t in range(24):
        expect(by_id["s-burst-%02d" % t]["status"] == "unique",
               "s-burst-%02d should be unique" % t)

    # 多解：两个不同且都通过复核的见证
    sm = by_id["s-multi"]
    expect(sm["matrix"]["nullity"] == 1, "s-multi nullity should be 1")
    expect(sm["solution_space"]["solutions_per_position"] == 256, "s-multi should have 256 solutions/pos")
    expect(sm["witnesses_differ"] is True, "s-multi witnesses must differ")
    expect(all(w["verified"] for w in sm["witnesses"]), "s-multi witnesses must verify")
    w1 = {b: v["sha256"] for b, v in sm["witnesses"][0]["blocks"].items()}
    w2 = {b: v["sha256"] for b, v in sm["witnesses"][1]["blocks"].items()}
    expect(w1 != w2, "s-multi witness digests must differ")
    expect("enumeration" in sm and "possible_values" in sm["enumeration"],
           "s-multi small stripe should be enumerated")

    # 矛盾：指认涉及的现有块
    sc = by_id["s-contra"]
    expect(sc["contradiction"]["inconsistent_position_count"] >= 1, "s-contra needs bad positions")
    expect("s-contra/d2.bin" in sc["contradiction"]["known_blocks_involved"],
           "s-contra should implicate d2")

    # 唯一恢复的块：重建后摘要与预存摘要一致
    with open(mpath, "r", encoding="utf-8") as f:
        manifest = json.load(f)
    man = {s["id"]: s for s in manifest["stripes"]}
    for sid in ("s-unique-missing", "s-unique-corrupt", "s-unique-parity"):
        for rec in by_id[sid]["recovered"]:
            expect(rec["matches_manifest"] is True,
                   "%s recovered %s must match manifest digest" % (sid, rec["path"]))
            rp = os.path.join(OUT, rec["path"])
            expect(os.path.isfile(rp), "%s must exist in rebuild output" % rp)
            if os.path.isfile(rp):
                with open(rp, "rb") as f:
                    expect(sha(f.read()) == man[sid]["sha256"][rec["path"]],
                           "%s rebuilt %s digest mismatch" % (sid, rec["path"]))
    for t in range(24):
        sid = "s-burst-%02d" % t
        for rec in by_id[sid]["recovered"]:
            expect(rec["matches_manifest"] is True,
                   "%s recovered %s must match manifest digest" % (sid, rec["path"]))

    # 原块从不改写：输入目录中每个文件（含故意腐坏的块）的摘要与运行前快照一致
    for rel, dg in snap.items():
        p = os.path.join(SAMPLE, rel)
        expect(os.path.isfile(p), "input file %s disappeared!" % rel)
        if os.path.isfile(p):
            with open(p, "rb") as f:
                expect(sha(f.read()) == dg, "input block %s was modified!" % rel)

    # 缓存：同一擦除模式的实际命中次数与峰值
    c = rep["cache"]
    expect(c["lookups"] == 30, "cache lookups should be 30, got %d" % c["lookups"])
    expect(c["misses"] == 7, "cache misses should be 7, got %d" % c["misses"])
    expect(c["hits"] == 23, "cache hits should be 23, got %d" % c["hits"])
    expect(c["peak_entries"] == 7, "cache peak should be 7, got %d" % c["peak_entries"])

    print("-" * 70)
    if fails:
        for m_ in fails:
            print("CHECK FAIL:", m_)
        return 1
    print("all expectation checks passed")
    print("cache evidence: hits=%d misses=%d peak_entries=%d"
          % (c["hits"], c["misses"], c["peak_entries"]))
    print("s-multi witness1 digests:", json.dumps(w1, indent=1))
    print("s-multi witness2 digests:", json.dumps(w2, indent=1))
    print("s-contra contradiction:", json.dumps(sc["contradiction"], indent=1))
    return 0


def main():
    mpath = make_pack()
    print("sample pack created at: %s" % SAMPLE)
    snap = snapshot_input()
    run_tool(mpath)
    return check_expectations(mpath, snap)


if __name__ == "__main__":
    sys.exit(main())
