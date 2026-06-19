#!/usr/bin/env python3
"""test_04: ficalib bias + dark + flat (no exptime correction)"""
import sys, os, numpy as np
TESTDIR=os.path.dirname(os.path.realpath(__file__)); ROOT=os.path.dirname(os.path.dirname(TESTDIR))
sys.path.insert(0,ROOT); os.chdir(TESTDIR)
from pyfitsh import Ficalib; from astropy.io import fits
P="ficalib_test"
def test():
    sci=fits.getdata(f"{P}_sci_a.fits").astype(np.float64)
    bias=fits.getdata(f"{P}_bias.fits").astype(np.float64)
    dark=fits.getdata(f"{P}_dark.fits").astype(np.float64)
    flat=fits.getdata(f"{P}_flat.fits").astype(np.float64)
    r=Ficalib().calibrate(sci,bias=bias,dark=dark,flat=flat)
    cli=fits.getdata(f"{P}_bdf_noexp.fits").astype(np.float64)
    return np.allclose(r['data'],cli,rtol=1e-5,atol=1e-8)
if __name__=="__main__": print(f"{'PASS' if test() else 'FAIL'}")
