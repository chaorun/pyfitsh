#!/usr/bin/env python3
"""
run_match.py -- grmatch 管道测试脚本 (Cython 版本)。

功能: 从 testdata/ 读取原始数据, 调用 Cython 的 Matcher,
       输出匹配结果并与 C dump 基线对比。
构建: python3 setup.py build_ext --inplace
用法: python3 run_match.py
"""
import sys
import time
import struct
import os
from grmatch_cy import Matcher
PARAMS = {
    "ref_file": "testdata/temp_ref.txt",  # 参考星表
    "inp_file": "testdata/temp_obs.txt",  # 目标星表
    "ref_cols": "2,3",   # col2=X, col3=Y (1-indexed)
    "inp_cols": "2,3",   # 同上
    "ref_ord":  "-4",    # col4 升序 (负号)
    "inp_ord":  "-4",    # 同上
    "max_dist": 2.0,     # --max-distance
    "order":    2,       # -a / --order
    "tri_mode": "auto",  # --triangulation
    "tri_unitarity": 0.01,  # unitarity=0.01
}
def parse_colspec(spec):
    """'2,3' -> [1,2] (0-indexed)"""
    return [int(s) - 1 for s in spec.split(",")]
def parse_ordering(spec):
    """'-4' -> idx=3, ascending=True"""
    v = int(spec)
    if v > 0:
        return v - 1, False
    return -v - 1, True
def read_points(filepath, colspec, ordspec):
    """读取数据文件, 返回 (xlist, ylist, ordlist)。"""
    cols = parse_colspec(colspec)
    ord_col, ascending = parse_ordering(ordspec)
    xs, ys, ords = [], [], []
    with open(filepath) as f:
        for raw in f:
            line = raw.split("#")[0].strip()
            if not line:
                continue
            tokens = line.split()
            try:
                x = float(tokens[cols[0]])
                y = float(tokens[cols[1]])
                w = float(tokens[ord_col])
                if ascending:
                    w = -w  # 负号 = 升序
                xs.append(x)
                ys.append(y)
                ords.append(w)
            except (IndexError, ValueError):
                continue
    return xs, ys, ords
def compare_baseline(nhit, vfits_dx, vfits_dy, stats):
    """对比 testdata/dump/ 中的 C grmatch 输出 (如已生成)。"""
    baseline = "testdata/dump/do_pointmatch_after"
    if not os.path.exists(f"{baseline}.vfits.bin"):
        print("\n  (No baseline dump -- skipped.)")
        return
    with open(f"{baseline}.hits.bin", "rb") as f:
        bl_nhit = struct.unpack("<i", f.read(4))[0]
    with open(f"{baseline}.vfits.bin", "rb") as f:
        f.read(4)  # skip unused header int
        bl_nvar = struct.unpack("<i", f.read(4))[0]
        bl_dx = list(struct.unpack(f"<{bl_nvar}d", f.read(bl_nvar * 8)))
        bl_dy = list(struct.unpack(f"<{bl_nvar}d", f.read(bl_nvar * 8)))
    print()
    print("=" * 60)
    print("COMPARISON WITH C grmatch (BIN dump)")
    print("=" * 60)
    vf_ok = True
    for i in range(min(len(vfits_dx), bl_nvar)):
        if abs(vfits_dx[i] - bl_dx[i]) > 1e-14:
            print(f"  VFITS dx[{i}] MISMATCH: {vfits_dx[i]:.15g} vs {bl_dx[i]:.15g}")
            vf_ok = False
        if abs(vfits_dy[i] - bl_dy[i]) > 1e-14:
            print(f"  VFITS dy[{i}] MISMATCH")
            vf_ok = False
    match_ok = (nhit == bl_nhit)
    print(f"  VFITS: {'ALL MATCH' if vf_ok else 'MISMATCH'}")
    print(f"  Matched: {nhit} (baseline: {bl_nhit}) {'OK' if match_ok else 'MISMATCH'}")
    print("=" * 60)
def main():
    """主流程: 读文件 → 匹配 → 输出 → 对比基线。"""
    print("=== grmatch Cython pipeline ===", file=sys.stderr)
    print(f"reading {PARAMS['ref_file']} ...", file=sys.stderr)
    ref_x, ref_y, ref_o = read_points(PARAMS["ref_file"],
                                      PARAMS["ref_cols"],
                                      PARAMS["ref_ord"])
    print(f"  -> {len(ref_x)} valid", file=sys.stderr)
    print(f"reading {PARAMS['inp_file']} ...", file=sys.stderr)
    inp_x, inp_y, inp_o = read_points(PARAMS["inp_file"],
                                      PARAMS["inp_cols"],
                                      PARAMS["inp_ord"])
    print(f"  -> {len(inp_x)} valid", file=sys.stderr)
    m = Matcher(order=PARAMS["order"],    # 多项式阶数
                maxdist=PARAMS["max_dist"],  # 最大距离
                unitarity=PARAMS["tri_unitarity"])  # unitarity 阈值
    m.set_reference(ref_x, ref_y, ref_o)
    m.set_input(inp_x, inp_y, inp_o)
    print(f"running do_pointmatch (order={PARAMS['order']}) ...",
          file=sys.stderr)
    t0 = time.time()
    r = m.run()
    dt = time.time() - t0
    print(f"  -> {r.nhit} matches in {dt:.3f}s", file=sys.stderr)
    nref, ninp = len(ref_x), len(inp_x)
    rat = r.nhit / min(nref, ninp) * 100
    print()
    print("=" * 60)
    print("POINT MATCHING RESULTS")
    print("=" * 60)
    print(f"  Matched:  {r.nhit} / {nref} (ref) / {ninp} (inp) = {rat:.1f}%")
    ns_val = r.stats['nsigma']
    ws_val = r.stats['wsigma']
    print(f"  Residual: {ns_val:.6f} (native)  {ws_val:.6f} (weighted)")
    print(f"  Unitarity: {r.stats['unitarity']:.8f}")
    print(f"  Hull cov:  {r.stats['hull_coverage']:.4f}")
    print(f"  Tri lvl:   {r.stats['tri_level']}  iters: {r.stats['nmiter']}")
    dx_list = ', '.join(f'{v:.10g}' for v in r.vfits_dx)
    dy_list = ', '.join(f'{v:.10g}' for v in r.vfits_dy)
    print(f"  VFITS dx:  [{dx_list}]")
    print(f"  VFITS dy:  [{dy_list}]")
    print("=" * 60)
    compare_baseline(r.nhit, r.vfits_dx, r.vfits_dy, r.stats)
    return 0
if __name__ == "__main__":
    sys.exit(main())
