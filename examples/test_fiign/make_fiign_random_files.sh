#!/bin/bash
# make_fiign_random_files.sh — 生成 fiign 测试用的随机 FITS 文件

P="fiign_test"
SX=256
SY=256
RPY="python3 -c \"import numpy as np; import sys"

set -ex
cd "$(dirname "$0")"

FIIGN="/Users/chaorun/Code/Githubs/fitsh-0.9.4/origincode/src/fiign"

# --- 01: test image with positive values only ---
python3 -c "
import numpy as np
from astropy.io import fits
np.random.seed(42)
img = np.random.uniform(10, 1000, (${SY}, ${SX})).astype(np.float32)
fits.writeto('${P}_01_input.fits', img, overwrite=True)
"

# --- 02: image with negative and zero values ---
python3 -c "
import numpy as np
from astropy.io import fits
np.random.seed(42)
img = np.random.uniform(-50, 200, (${SY}, ${SX})).astype(np.float32)
img[10:20, 10:20] = 0.0
img[30:40, 30:40] = -10.0
fits.writeto('${P}_02_mixed.fits', img, overwrite=True)
"

# --- 03: image with some very bright pixels (saturation test) ---
python3 -c "
import numpy as np
from astropy.io import fits
np.random.seed(42)
img = np.random.uniform(100, 1000, (${SY}, ${SX})).astype(np.float32)
img[50, 50] = 60000.0
img[51, 51] = 62000.0
img[52, 50] = 58000.0
fits.writeto('${P}_03_sat.fits', img, overwrite=True)
"

# --- 04: image with cosmic-ray-like spikes ---
python3 -c "
import numpy as np
from astropy.io import fits
np.random.seed(42)
img = np.random.normal(100, 5, (${SY}, ${SX})).astype(np.float32)
img[100, 100] = 500.0   # very bright cosmic
img[100, 101] = 3.0     # very dark cosmic
img[150, 150] = 800.0
img[200, 200] = 1.0
fits.writeto('${P}_04_cosmics.fits', img, overwrite=True)
"

# --- 05: image with mask block test ---  
python3 -c "
import numpy as np
from astropy.io import fits
np.random.seed(42)
img = np.random.uniform(100, 1000, (${SY}, ${SX})).astype(np.float32)
fits.writeto('${P}_05_block.fits', img, overwrite=True)
"

# --- 06: mask input (external mask file) ---
python3 -c "
import numpy as np
from astropy.io import fits
np.random.seed(42)
img = np.random.uniform(100, 1000, (${SY}, ${SX})).astype(np.float32)
fits.writeto('${P}_06_maskin.fits', img, overwrite=True)
mask = np.zeros((${SY}, ${SX}), dtype=np.uint8)
mask[60:80, 60:80] = 1   # MASK_FAULT
fits.writeto('${P}_06_mask.fits', mask, overwrite=True)
"

# --- 07: integer limit test (bitpix=8) ---
python3 -c "
import numpy as np
from astropy.io import fits
np.random.seed(42)
img = np.random.uniform(-200, 300, (${SY}, ${SX})).astype(np.float32)
img[10, 10] = -200.0
img[20, 20] = 300.0
fits.writeto('${P}_07_intlimit.fits', img, overwrite=True)
"

echo "All test data generated."
ls -la ${P}_*.fits
