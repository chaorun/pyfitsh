#!/usr/bin/env python3
"""test_01: grmatch 点匹配对比"""
import sys, os, numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '..')
from compare_utils import run_test
from astropy.table import Table

REF_CAT = "ref_L1_GRB260604C-Ic-180s-20260606_210318.cat"
OBS_CAT = "obs_L1_GRB260604C-Ic-180s-20260607_023443.cat"


def do_test():
    from pyfitsh import Grmatch

    ref_tab = Table.read(REF_CAT, format='ascii.no_header',
                         names=['id', 'x', 'y', 'mag'])
    inp_tab = Table.read(OBS_CAT, format='ascii.no_header',
                         names=['id', 'x', 'y', 'mag'])

    ref_x = np.asarray(ref_tab['x'], dtype=np.float64)
    ref_y = np.asarray(ref_tab['y'], dtype=np.float64)
    ref_mag = np.asarray(ref_tab['mag'], dtype=np.float64)
    ref_id = np.asarray(ref_tab['id'], dtype=np.int32)
    inp_x = np.asarray(inp_tab['x'], dtype=np.float64)
    inp_y = np.asarray(inp_tab['y'], dtype=np.float64)
    inp_mag = np.asarray(inp_tab['mag'], dtype=np.float64)
    inp_id = np.asarray(inp_tab['id'], dtype=np.int32)

    g = Grmatch(order=1, maxdist=2.0, unitarity=0.2, ttype=2,
                use_ordering=1, is_centering=0, nmiter=4, rejlevel=3.0)

    r = g.matchpoints(ref_x, ref_y, inp_x, inp_y,
                      ref_order=-ref_mag, inp_order=-inp_mag)

    matched = []
    for ri, ii in zip(r.hits_ref, r.hits_inp):
        matched.append([
            ref_id[ri], ref_x[ri], ref_y[ri], ref_mag[ri],
            inp_id[ii], inp_x[ii], inp_y[ii], inp_mag[ii],
        ])

    cols = ['rid', 'xr', 'yr', 'mr', 'iid', 'xi', 'yi', 'mi']
    out_table = Table(rows=matched, names=cols)

    cli_tab = Table.read("cli_test_match.cat", format='ascii.no_header',
                         names=cols)
    cy_data = np.column_stack([np.asarray(out_table[c], dtype=np.float64)
                                for c in cols])
    cli_data = np.column_stack([np.asarray(cli_tab[c], dtype=np.float64)
                                 for c in cols])

    if cy_data.shape != cli_data.shape:
        return False, f"match table shape: CY={cy_data.shape} CLI={cli_data.shape}"
    if not np.allclose(cy_data, cli_data, rtol=1e-5, atol=1e-8, equal_nan=True):
        diff = np.abs(cy_data - cli_data)
        return False, f"match table max_diff={np.max(diff):.2e}"

    # grb.trans
    trans = r.transformation
    trans_lines = [
        f"type = {trans['type']}",
        f"order = {trans['order']}",
        f"offset = {trans['offset'][0]}, {trans['offset'][1]}",
        f"scale = {trans['scale']}",
        f"dxfit=     {', '.join(str(v) for v in trans['dxfit'])}",
        f"dyfit=     {', '.join(str(v) for v in trans['dyfit'])}",
    ]
    with open("cli_test_grb.trans", "r") as f:
        cli_lines = [l.rstrip("\n") for l in f if l.strip()]
    if len(trans_lines) != len(cli_lines):
        return False, f"trans lines: CY={len(trans_lines)} CLI={len(cli_lines)}"

    return True, f"nhit={r.nhit} ok"


if __name__ == "__main__":
    sys.exit(run_test(do_test))
