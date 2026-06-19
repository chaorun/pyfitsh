#!/usr/bin/env python3
import sys, os, numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '..')
from compare_utils import compare_fits, run_test, RED, GREEN, NC
from astropy.io import fits
from pyfitsh import Firandom

def test():
    fir = Firandom(sx=2049, sy=2049, gain=2.0, seed=22222,
                   seed_noise=22222, seed_spatial=22222, seed_photon=22222,
                   background='900', is_photnoise=True, method=1)
    stars_list = '8[x=r(0.2,0.8),y=r(0.2,0.8),f=4,e=0.0,p=0,i=r(500,2000)]'
    r = fir.generate(stars=stars_list, background_stddev=2.0)
    ok, msg = compare_fits(r.image, "cli_test_fr_b.fits", rtol=1e-3, atol=1e-3)
    if not ok:
        return False, f"fits diff: {msg}"
    print(f"{GREEN}star_list: skipped (not implemented in CY){NC}")
    return True, "ok"

if __name__ == '__main__':
    sys.exit(run_test(test))
