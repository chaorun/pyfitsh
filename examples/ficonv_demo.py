#!/usr/bin/env python3
"""ficonv 完整用法示例：三种核型（Gauss / 离散 / 空间变化）及 kernel_dict 跳过拟合"""

import sys, os, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from pyfitsh.ficonv.ficonv import Ficonv
from pyfitsh.utils import decode_maskinfo
from astropy.io import fits

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'fitsh_test_data')

ref = fits.getdata(os.path.join(DATA, 'L1_GRB260604C-Ic-180s-20260606_210318.fits')).astype(np.float64)

with fits.open(os.path.join(DATA, 'v1.fits')) as hdul:
    img = hdul[0].data.astype(np.float64)
    img_mask = decode_maskinfo(hdul[0])  # CLI 对比必须

# ========================================================================
# 核型 1: 多高斯核（等价 -k "i/0;b/1;g=3,1.0,2/0;g=5,2.0,2/0;g=7,4.0,2/0"）
# ========================================================================
print("=== 核型 1: 多高斯 ===")
fc = Ficonv(kernel='i/0;b/1;g=3,1.0,2/0;g=5,2.0,2/0;g=7,4.0,2/0',
            iterations=3, rejection_level=3.0)
t0 = time.time()
conv, msk, sub, kd = fc.fit(ref, img, img_mask=img_mask,
                             output_subtracted=True, output_kernel_list=True)
print(f"  耗时 {time.time()-t0:.1f}s, {len(kd['kernels'])} kernels, ref_conv shape={conv.shape}")

# kernel_dict 跳过拟合
conv2, msk2 = fc.fit(ref, img, kernel_dict=kd)
print(f"  kdict skip: fit vs skip max_diff={np.abs(conv-conv2).max():.2e}")

# ========================================================================
# 核型 2: 离散核（等价 -k "i/0;b/0;d=2/0"）
# ========================================================================
print("\n=== 核型 2: 离散核 ===")
fc2 = Ficonv(kernel='i/0;b/0;d=2/0', iterations=2, rejection_level=3.0)
t0 = time.time()
conv_d, _, sub_d, kd_d = fc2.fit(ref, img, img_mask=img_mask,
                                   output_subtracted=True, output_kernel_list=True)
print(f"  耗时 {time.time()-t0:.1f}s, {len(kd_d['kernels'])} kernels")
conv_d2, _ = fc2.fit(ref, img, kernel_dict=kd_d)
print(f"  kdict skip: max_diff={np.abs(conv_d-conv_d2).max():.2e}")

# ========================================================================
# 核型 3: 空间变化核（等价 -k "i/1;b/1;g=5,2.0,2/1"）
# ========================================================================
print("\n=== 核型 3: 空间变化核 ===")
fc3 = Ficonv(kernel='i/1;b/1;g=5,2.0,2/1', iterations=2, rejection_level=3.0)
t0 = time.time()
conv_s, _, sub_s, kd_s = fc3.fit(ref, img, img_mask=img_mask,
                                   output_subtracted=True, output_kernel_list=True)
print(f"  耗时 {time.time()-t0:.1f}s, {len(kd_s['kernels'])} kernels, 每核 {len(kd_s['kernels'][0]['coeff'])} 系数")
conv_s2, _ = fc3.fit(ref, img, kernel_dict=kd_s)
print(f"  kdict skip: max_diff={np.abs(conv_s-conv_s2).max():.2e}")

# ========================================================================
# 检查
# ========================================================================
errors = []
for name, a, b in [('多高斯', conv, conv2), ('离散核', conv_d, conv_d2), ('空间变化', conv_s, conv_s2)]:
    d = np.abs(a - b).max()
    ok = d < 1e-10
    errors += [not ok]
    print(f"\n{name}: fit vs kdict = {d:.2e}  {'PASS' if ok else 'FAIL'}")

if any(errors):
    raise SystemExit(1)
print("\nALL PASS")
