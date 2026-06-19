#!/usr/bin/env python3
"""test_08: ficalib trim (FITS coord 100:100:1900:1900 → numpy [100:1901,100:1901])"""
import sys, os, numpy as np
TESTDIR=os.path.dirname(os.path.realpath(__file__)); ROOT=os.path.dirname(os.path.dirname(TESTDIR))
sys.path.insert(0,ROOT); os.chdir(TESTDIR)
from pyfitsh import Ficalib; from astropy.io import fits
P="ficalib_test"
def test():
    # CLI --image 100:100:1900:1900 uses FITS 1-indexed coords
    # Equivalent numpy slice: [100:1901, 100:1901] (0-indexed inclusive)
    sci=fits.getdata(f"{P}_sci_a.fits").astype(np.float64)[100:1901,100:1901]
    r=Ficalib().calibrate(sci)
    cli=fits.getdata(f"{P}_trim.fits").astype(np.float64)
    return np.allclose(r['data'],cli,rtol=1e-5,atol=1e-8)
if __name__=="__main__": print(f"{'PASS' if test() else 'FAIL'}")
