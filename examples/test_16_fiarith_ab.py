#!/usr/bin/env python3
import sys, os, numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '..')
from compare_utils import compare_table, compare_fits, run_test, RED, GREEN, NC
from astropy.io import fits
from pyfitsh.fiarith import Fiarith

def test():
    a_data = fits.getdata("cli_test_fr_a.fits").astype(np.float64)
    b_data = fits.getdata("cli_test_fr_b.fits").astype(np.float64)
    fa = Fiarith()
    result = fa.evaluate('a-b', {'a': a_data, 'b': b_data})
    ok, msg = compare_fits(result, "cli_test_fiarith_a_minus_b.fits")
    if not ok:
        return False, msg
    return True, "ok"

if __name__ == '__main__':
    sys.exit(run_test(test))
