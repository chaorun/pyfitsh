#!/usr/bin/env python3
"""ficonv usage examples — matching ficonv --long-help keywords.

Requires: fitsh_cython built (python3 setup.py build_ext --inplace)
"""
import numpy as np
from fitsh_cy import Fitrans

# ---- mock data: 128x128 images ----
rng = np.random.default_rng(42)
ref = np.ones((128, 128), dtype=np.float64) * 1000
# img = ref * 1.5 + background gradient + poisson noise
yy, xx = np.mgrid[:128, :128]
bg = 200 + 0.3 * xx + 0.2 * yy
img = (ref * 1.5 + bg + rng.normal(0, 10, (128, 128))).astype(np.float64)

# ---- mask: mark edges ----
ref_mask = np.zeros((128, 128), dtype=np.uint8)
ref_mask[:2, :] = 1; ref_mask[-2:, :] = 1; ref_mask[:, :2] = 1; ref_mask[:, -2:] = 1
img_mask = ref_mask.copy()

# ===================================================================
# 1. Basic convolution (normal mode)
# ===================================================================
c = Fitrans.Convolution(
    kernel='i/0;b/1;g=3,1.0,2/0',   # -k: identity + background + gaussian
    iterations=3,                      # -n: 3 rejection iterations
    rejection_level=3.0,               # -s: 3-sigma rejection
    verbose=True,                      # --verbose
)
cnv, mask, kl = c.fit(ref, img, output_kernel_list=True)
print("1. Basic fit:", cnv.shape, "kernels:", kl[:120])

# ===================================================================
# 2. Masked mode — fit only foreground pixels
# ===================================================================
c2 = Fitrans.Convolution(
    kernel='i/0;b/1',
    masked=True,                       # -m: masked fitting
    divide=16,                         # -d: block divide
)
cnv2, mask2 = c2.fit(ref, img, ref_mask=ref_mask, img_mask=img_mask)
print("2. Masked:", cnv2.shape, "masked pix:", mask2.sum())

# ===================================================================
# 3. Weighted mode — weights from image noise (implies masked)
# ===================================================================
c3 = Fitrans.Convolution(
    kernel='i/0;b/1;g=3,1.0',
    weighted=True,                     # -w: weighted fitting (implies -m)
    gain=1.5,                          # -g: detector gain
)
cnv3, mask3 = c3.fit(ref, img, ref_mask=ref_mask, img_mask=img_mask)
print("3. Weighted:", cnv3.shape)

# ===================================================================
# 4. Background-iterative mode
# ===================================================================
c4 = Fitrans.Convolution(
    kernel='i/0;b/1',
    background_iterative=True,         # -b: background-iterative fitting
    iterations=5,
)
cnv4, mask4 = c4.fit(ref, img)
print("4. Bg-iterative:", cnv4.shape)

# ===================================================================
# 5. Subtract output — useful for difference imaging
# ===================================================================
c5 = Fitrans.Convolution(kernel='i/0;b/1;g=3,1.0,2/0')
cnv5, mask5, sub5 = c5.fit(ref, img, output_subtracted=True)
print(f'5. Subtracted RMS: {sub5.std():.2f}')

# ===================================================================
# 6. Add-to — add template to convolved result
# ===================================================================
add_img = np.ones_like(ref) * 500
c6 = Fitrans.Convolution(kernel='i/0;b/1')
cnv6, mask6, added6 = c6.fit(ref, img, add_to=add_img)
print(f'6. Add-to center: {added6[64, 64]:.1f}')

# ===================================================================
# 7. Stamps — fit only in specified regions
# ===================================================================
c7 = Fitrans.Convolution(kernel='i/0;b/1')
cnv7, mask7 = c7.fit(ref, img,
    stamps=[(32, 32, 32, 32), (80, 80, 32, 32)])  # -t: stamp regions
print("7. Stamps:", cnv7.shape)

# ===================================================================
# 8. Unity kernels — normalize identity kernel flux to 1.0
# ===================================================================
c8 = Fitrans.Convolution(
    kernel='i/0;b/1',
    unity_kernels=True,                # -u: normalize to unit flux
)
cnv8, mask8, kl8 = c8.fit(ref, img, output_kernel_list=True)
print("8. Unity kernels:", kl8[:120])

print("\n=== All examples complete ===")
