#!/usr/bin/env python3
import sys, os, numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '..')
from compare_utils import compare_table, compare_fits, run_test, RED, GREEN, NC
from astropy.io import fits
from pyfitsh.fiarith import Fiarith

def test():
    fa = Fiarith()
    img_a = fits.getdata("cli_test_fr_a.fits").astype(np.float64)
    result = fa.evaluate('laplace(a)', {'a': img_a})
    ok, msg = compare_fits(result, "cli_test_fiarith_laplace.fits")
    if not ok:
        return False, msg
    return True, "ok"

if __name__ == '__main__':
    sys.exit(run_test(test))
