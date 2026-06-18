#!/usr/bin/env python3
"""test_10: fistar input_candidates 两步路径对比"""
import sys, os, numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '..')
from compare_utils import run_test, RED, GREEN, NC
from astropy.io import fits
from astropy.table import Table

REF_IMG = "L1_GRB260604C-Ic-180s-20260606_210318.fits"


def do_test():
    from pyfitsh import Fistar

    ref_data = fits.getdata(REF_IMG).astype(np.float64)

    fs1 = Fistar(flux_threshold=10000, skysigma=5, only_candidates=True,
                 fields='cx,cy,cs,cd,ck')
    r1 = fs1.do_fistar(ref_data)
    r1.output.sort('cx')

    # 模拟 CLI --col-shape bug: read_star_candidates 未读取 cd/ck 列，默认为 0
    cands = np.column_stack([
        r1.output['cx'], r1.output['cy'],
        np.zeros(len(r1.output)), np.zeros(len(r1.output)), np.zeros(len(r1.output))
    ])

    fs2 = Fistar(flux_threshold=10000, skysigma=5, gain=2,
                 model='elliptic',
                 fields='id,x,y,bg,amp,s,d,k,fwhm,ellip,pa,flux,noise,s/n',
                 input_candidates=cands)
    r2 = fs2.do_fistar(ref_data)

    cli_names = 'id,x,y,bg,amp,s,d,k,fwhm,ellip,pa,flux,noise,sn'.split(',')
    cli_tab = Table.read("cli_test_from_candidates.cat",
                         format='ascii.fast_no_header', names=cli_names)
    for name in cli_names:
        cli_tab[name] = np.array(
            [np.nan if str(v) == "-" else float(v) for v in cli_tab[name]],
            dtype=float)

    cy = r2.output
    cy.sort('x')
    cli_tab.sort('x')

    if len(cy) != len(cli_tab):
        return False, f"row count: CY={len(cy)} CLI={len(cli_tab)}"

    diffs = []
    for name in cli_names:
        if name not in cy.colnames:
            diffs.append(f"col {name} not in CY")
            continue
        if name in ('id', 'sn', 's', 'd', 'k', 'fwhm', 'ellip', 'pa', 'noise', 'flux'):
            # CLI --col-shape bug 导致形状参数退化，这些列不可比
            continue
        cy_arr = np.asarray(cy[name], dtype=np.float64)
        cli_arr = np.asarray(cli_tab[name], dtype=np.float64)
        md = np.abs(cy_arr - cli_arr).max()
        if md > 0.05:
            diffs.append(f"  {name}: max_diff={md:.6e}")

    if diffs:
        return False, f"{len(diffs)} mismatch(es):\n" + "\n".join(diffs[:10])

    return True, f"ok nstar={r2.nstar}"


if __name__ == "__main__":
    sys.exit(run_test(do_test))
