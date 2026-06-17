#!/usr/bin/env python3
"""test_01: grmatch 点匹配对比"""
import sys, os, numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '..')
from compare_utils import compare_table, compare_fits, compare_cat, run_test, RED, GREEN, NC
from astropy.io import fits
from astropy.table import Table

REF_CAT = "ref_L1_GRB260604C-Ic-180s-20260606_210318.cat"
OBS_CAT = "obs_L1_GRB260604C-Ic-180s-20260607_023443.cat"

def do_test():
    from pyfitsh import Grmatch

    ref = np.loadtxt(REF_CAT, ndmin=2)
    inp = np.loadtxt(OBS_CAT, ndmin=2)

    ref_x = ref[:, 1]   # col 2 in CLI (0-indexed)
    ref_y = ref[:, 2]   # col 3 in CLI
    ref_mag = ref[:, 3] # col 4 in CLI
    inp_x = inp[:, 1]
    inp_y = inp[:, 2]
    inp_mag = inp[:, 3]

    g = Grmatch(order=1, maxdist=2.0, unitarity=0.2, ttype=2,
                use_ordering=1, is_centering=0, nmiter=4, rejlevel=3.0)

    r = g.matchpoints(ref_x, ref_y, inp_x, inp_y,
                      ref_order=-ref_mag, inp_order=-inp_mag)

    # 构建 output Table
    ref_ids = ref[:, 0].astype(int)
    inp_ids = inp[:, 0].astype(int)

    matched = []
    for ri, ii in zip(r.hits_ref, r.hits_inp):
        matched.append([
            ref_ids[ri], ref_x[ri], ref_y[ri], ref_mag[ri],
            inp_ids[ii], inp_x[ii], inp_y[ii], inp_mag[ii],
        ])

    cols = ['rid','xr','yr','mr','iid','xi','yi','mi']
    out_table = Table(rows=matched, names=cols)

    # CLI match.cat has no header; compare raw numpy arrays
    cli_data = np.loadtxt("cli_test_match.cat", dtype=np.float64)
    cy_data = np.column_stack([np.asarray(out_table[c], dtype=np.float64) for c in cols])
    if cy_data.shape != cli_data.shape:
        return False, f"match table shape: CY={cy_data.shape} CLI={cli_data.shape}"
    if not np.allclose(cy_data, cli_data, rtol=1e-5, atol=1e-8, equal_nan=True):
        diff = np.abs(cy_data - cli_data)
        return False, f"match table max_diff={np.max(diff):.2e} n_diff={np.sum(~np.isclose(cy_data,cli_data,rtol=1e-5,atol=1e-8))}"
    ok1, msg1 = True, "ok"
    if not ok1:
        return False, f"match table: {msg1}"

    # grb.trans 行数检查
    trans = r.transformation
    trans_lines = []
    trans_lines.append(f"type = {trans['type']}")
    trans_lines.append(f"order = {trans['order']}")
    trans_lines.append(f"offset = {trans['offset'][0]}, {trans['offset'][1]}")
    trans_lines.append(f"scale = {trans['scale']}")
    dx_str = ", ".join(str(v) for v in trans['dxfit'])
    dy_str = ", ".join(str(v) for v in trans['dyfit'])
    trans_lines.append(f"dxfit=     {dx_str}")
    trans_lines.append(f"dyfit=     {dy_str}")

    with open("cli_test_grb.trans", "r") as f:
        cli_lines = [l.rstrip("\n") for l in f if l.strip()]

    if len(trans_lines) != len(cli_lines):
        return False, f"trans line count: CY={len(trans_lines)} CLI={len(cli_lines)}"

    return True, f"nhit={r.nhit} trans lines={len(trans_lines)} ok"

if __name__ == "__main__":
    sys.exit(run_test(do_test))
