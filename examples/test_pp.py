import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from fitsh_cython.fitsh_cy import Fistar
from astropy.io import fits
img = fits.getdata('../fitsh_test_data/L1_GRB260604C-Ic-180s-20260606_210318.fits').astype(np.float64)
fs = Fistar(threshold=100.0, flux_threshold=10000.0, skysigma=5.0,
            algorithm='parabolapeak', only_candidates=True, verbose=False,
            output_mark=False, output_area=False)
print("calling do_fistar...")
r = fs.do_fistar(img)
print(f"nstar={r['nstar']}")
