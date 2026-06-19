#!/usr/bin/env python3
"""
CLI fiphot 全流程 vs Cython 全流程对比测试

CLI 命令：
  # 第一步：标准测光（带 --sky-fit）
  fiphot --input L1_GRB260604C-Ic-180s-20260606_210318.fits \
    --input-list fiphot_positions.dat --col-id 1 --col-xy 2,3 \
    --apertures 10:18:12 --gain 2 \
    --sky-fit median,iterations=2,sigma=3 \
    --format IXY,BbFfMm --output ref_rawcheck.cat \
    --output-raw-photometry ref_rawcheck.raw

  # 第二步：减影测光（无 --sky-fit，使用默认 bgmode）
  fiphot -s sub_spatial1.fits \
    --input-raw-photometry ref_rawcheck.raw \
    --input-kernel kernel_spatial1.txt \
    --format IXY,BbFfMm --output fiphot_sub_spatial1.cat

  # 第三步：ficonv 核拟合（产生 kernel_spatial1.txt 和 sub_spatial1.fits）
  ficonv -r L1_GRB260604C-Ic-180s-20260606_210318.fits -i v1.fits \
    -k "i/1;b/1;g=5,2.0,2/1" --output-kernel-list kernel_spatial1.txt \
    -o ref_conv_spatial1.fits --output-subtracted sub_spatial1.fits -n 2 -s 3

关键发现：
  CLI 两步使用不同 bgmode：第一步有 --sky-fit，第二步无 --sky-fit。
  Cython 必须用两个 Fiphot 实例分别设置，否则 bgmode 不一致导致偏差。
"""
import sys, numpy as np
sys.path.insert(0, '/Users/chaorun/Code/Githubs/fitsh-0.9.4/testgrmatch')
from fitsh_cython.ficonv.ficonv import Ficonv
from fitsh_cython.fitsh_cy import Fiphot
from fitsh_cython.utils import decode_maskinfo
from astropy.io import fits

DATA = '/Users/chaorun/Code/Githubs/fitsh-0.9.4/testgrmatch/fitsh_test_data'

ref = fits.getdata(DATA+'/L1_GRB260604C-Ic-180s-20260606_210318.fits').astype(np.float64)
img = fits.getdata(DATA+'/v1.fits').astype(np.float64)
with fits.open(DATA+'/v1.fits') as hdul:
    img_mask = decode_maskinfo(hdul[0])
img_mask[np.isnan(img)] = 1

stars = np.loadtxt(DATA+'/fiphot_positions.dat')

# 第一步：photometry（CLI 有 --sky-fit）
fp1 = Fiphot(apertures='10:18:12', gain='2', mag_flux=(10.0, 10000.0),
             sky_fit='median,iterations=2,sigma=3')

# 第二步：photometry_from_raw（CLI 无 --sky-fit）
fp2 = Fiphot(apertures='10:18:12', gain='2', mag_flux=(10.0, 10000.0))

r = fp1.photometry(ref, stars=stars, col_xy=(1,2), col_id=0)

fc = Ficonv(kernel='i/1;b/1;g=5,2.0,2/1', iterations=2, rejection_level=3.0)
_, _, _, kd = fc.fit(ref, img, img_mask=img_mask, output_subtracted=True, output_kernel_list=True)

sub = fits.getdata(DATA+'/sub_spatial1.fits').astype(np.float64)
with fits.open(DATA+'/sub_spatial1.fits') as hdul:
    sub_mask = decode_maskinfo(hdul[0])
r2 = fp2.photometry_from_raw(r.dict['output_raw_photometry_3d_data'], sub, mask=sub_mask, kernel_dict=kd)

# 按 ID 对齐对比
cy_by_id = {}
for i, sid in enumerate(r.output['id']):
    cy_by_id[int(sid)] = (r2.dict['flux'][i, 0], r2.dict['bgmedian'][i, 0])

cli_by_id = {}
with open(DATA+'/fiphot_sub_spatial1.cat') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        p = line.split()
        try:
            cli_by_id[int(float(p[0]))] = (float(p[5]), float(p[3]))
        except:
            pass

common = sorted(set(cy_by_id) & set(cli_by_id))
diffs = []
for s in common:
    cf, cb = cli_by_id[s]
    yf, yb = cy_by_id[s]
    diffs.append((abs(cf - yf), cf, yf, cb, yb, s))
diffs.sort(reverse=True)

print(f'{"星号":>4s} | {"CLI flux":>13s} | {"CY flux":>13s} | {"diff":>8s} | {"ratio":>8s}')
print('-' * 70)
for d, cf, yf, cb, yb, s in diffs[:10]:
    r = yf / cf if cf != 0 else 0
    print(f'  {s:3d} | {cf:13.2f} | {yf:13.2f} | {d:8.2f} | {r:8.4f}')

dr = np.array([d[0] for d in diffs])
print('-' * 70)
print(f'{len(common)}颗星, max={dr.max():.2f}, mean={dr.mean():.4f}, median={np.median(dr):.4f}')
print()
print('PASS' if dr.max() < 100 else 'FAIL')
