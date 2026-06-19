import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from fitsh_cython.fitsh_cy import Fistar
from astropy.io import fits
import numpy as np
img = fits.getdata('../fitsh_test_data/L1_GRB260604C-Ic-180s-20260606_210318.fits').astype(np.float64)
fs = Fistar(threshold=100.0, algorithm='uplink', model='elliptic', model_order=2,
            it_sym=4, it_gen=2, verbose=False, output_mark=False, output_area=False)
print('start')
r = fs.do_fistar(img)
print(f'nstar={r["nstar"]}')
