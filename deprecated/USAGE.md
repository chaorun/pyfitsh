# fitsh_cython 使用指南 v0.02

## 安装

```bash
pip install cython numpy astropy
cd fitsh_cython
python3 setup.py build_ext --inplace
```

## 导入

```python
from fitsh_cy import (Grmatch, MatchResult, Grtrans, Ficonv,
                      PolyFitter, Fitrans, Fistar,
                      parse_trans_file)
```

---

## 1. 星表匹配 (Grmatch → grmatch)

```python
import numpy as np
from fitsh_cy import Grmatch

# 加载星表: (id, x, y, mag) 格式
ref = np.loadtxt('ref.cat')
obs = np.loadtxt('obs.cat')

m = Grmatch(order=2, maxdist=2.0, unitarity=0.01, verbose=True)
m.set_reference(ref[:,1], ref[:,2], -ref[:,3])  # X, Y, 排序值(负星等=亮星优先)
m.set_input(obs[:,1], obs[:,2], -obs[:,3])
result = m.run()

print(result.nhit)           # 匹配对数
print(result.vfits_dx)       # 变换系数 (X)
print(result.transform)      # dict, 可直接传入 Grtrans
```

---

## 2. 坐标变换 (Grtrans → grtrans)

```python
# 从匹配结果
t = Grtrans(result.transform)
new_x, new_y = t.apply(old_x, old_y)

# 从 .trans 文件
t = Grtrans(parse_trans_file('grb.trans'))

# 与 C .trans 格式互转
t = Grtrans.from_trans_string(trans_str)
s = t.to_trans_string()

# 恒等变换
t_id = Grtrans(order=1, dxfit=[0,1,0], dyfit=[0,0,1])
```

### WCS 拟合

```python
w = Grtrans.WCSFitter(ra0=189.0, dec0=62.0, proj='tan', order=3)
w.fit(ra_matched, dec_matched, x_matched, y_matched)
xi, eta = w.project(ra, dec)
x, y = w.proj_to_pix(xi, eta)
```

### 多项式拟合

```python
f = PolyFitter(order=2, niter=3, rejlevel=3.0)
f.fit(x, y, values)
print(f.coeff, f.residual, f.nrej)
```

---

## 3. 卷积核拟合 (Ficonv → ficonv)

```python
from fitsh_cy import Ficonv

# 基本用法
fc = Ficonv(kernel='i/0; b/0; g=3,1.0,0/0')
cnv, mask = fc.fit(ref_data, img_data)

# 带迭代剔点
fc = Ficonv(kernel='i/0; b/1; g=5,2.0,2/0', iterations=3, rejection_level=3.0)
cnv, mask = fc.fit(ref, img)

# 带 mask + 附加输出
cnv, mask, subtracted = fc.fit(ref, img, output_subtracted=True)
cnv, mask, kernel_list = fc.fit(ref, img, output_kernel_list=True)

# stamp 限定区域
cnv, mask = fc.fit(ref, img, stamps=[(x, y, sx, sy), ...])

# 也可用 Fitrans.Convolution(...), 接口相同
```

### kernel 语法 (与 CLI 一致)

```
i/<spatial>                identity kernel (通量项)
b/<spatial>                background offset
g=<hsize>,<sigma>,<order>/<spatial>    Gaussian kernel
d=<hsize>/<spatial>        discrete kernel
```

多核用 `;` 分隔。示例:
- `'i/0;b/1'` — identity + 线性背景
- `'g=3,1.0,0/0'` — 单个 hsize=3, sigma=1.0 的高斯核
- `'g=3,1.0,0/0;g=5,2.0,0/0;g=7,4.0,0/0'` — 三尺度高斯核

---

## 4. 图像变换 (Fitrans → fitrans)

### 几何变换 (ImageTransformation)

```python
f = Fitrans.ImageTransformation(method='lanczos3', transformation=result.transform)
out, mask = f.apply(img)
out, mask = f.apply(img, out_shape=(nsy, nsx))
```

方法: `bilinear`, `integrate`, `bicubic`, `spline_integrate`, `lanczos3`, `lanczos4`.

### 噪声 / 放缩 / 扩展 / 平滑

```python
# 噪声估计
out, mask = Fitrans.Noise().apply(img)

# 放大
out, mask = Fitrans.Zoom(3).apply(img)            # 带子像素插值
out, mask = Fitrans.Zoom(2, raw=True).apply(img)   # 无插值

# 缩并
out, mask = Fitrans.Shrink(2).apply(img)           # 均值
out, mask = Fitrans.Shrink(2, median=True).apply(img)  # 中值

# 重复扩展 (--repetitive-xy)
out, mask = Fitrans.Repetitive(2, 3).apply(img)
out, mask = Fitrans.Repetitive(2, 3, offset_x=5, offset_y=10).apply(img)

# 交错扩展 (--interleave-xy)
out, mask = Fitrans.Interleave(4, 6).apply(img)
out, mask = Fitrans.Interleave(4, 6, median=True).apply(img)

# 大尺度平滑
out, mask = Fitrans.Smooth('spline', order=2).apply(img)
out, mask = Fitrans.Smooth('polynomial', order=3, prefilter='median', niter=3).apply(img)
```

### 解码 MASKINFO

```python
from astropy.io import fits
mask = Fitrans.decode_maskinfo(fits.open('image.fits')[0])
# 返回 uint8 ndarray, 等价于 C: fits_mask_read_from_header()
```

---

## 5. 星检测与测光 (Fistar → fistar)

### 基本用法

```python
r = Fistar(
    threshold=50.0,
    flux_threshold=0.0,
    algorithm='uplink',
    model='gauss',
    sort='x',
    fields='x,y,flux,magnitude,fwhm',
    output_mark=True,
    output_area=True,
    psf='native,order=2',
    psf_output=True,
).search(img_data, mask=mask)

r['nstar']         # 检测到的星数
r['table']         # astropy.table.Table (星表)
r['output_mark']   # ndarray (标记图像)
r['output_area']   # ndarray (星区域图像)
r['psf']           # fits.ImageHDU (PSF 数据立方)
```

### 输入候选星 / 位置

```python
# 候选星 (替代 -C)
cands = np.array([[1024.0, 512.0], [800.0, 300.0]], dtype=np.float64)
r = Fistar(input_candidates=cands).search(img)

# 指定位置测光 (替代 -P)
r = Fistar(input_positions=[(1024.5, 512.3), (800.1, 300.7)]).search(img)
# r['table'] 含: input_id, input_x, input_y, mag, merr, flux, ferr
```

### 参数速查

| 参数 | 默认 | CLI |
|------|------|-----|
| `threshold` | 100.0 | `-t` |
| `algorithm` | `'uplink'` | `--algorithm` |
| `model` | `'elliptic'` | `--model` (gauss/elliptic) |
| `model_order` | 2 | `--model order=N` |
| `psf` | None | `--psf` (native/integral/circle) |
| `gain` | 1.0 | `-g` |
| `sort` | `'x'` | `-s` (x/y/peak/fwhm/amp/flux/noise/sn) |
| `mag_flux` | `(10.0, 10000.0)` | `--mag-flux` |
| `input_candidates` | None | `-C` |
| `input_positions` | None | `-P` |
| `candidate_radius` | 2.0 | `-R` |
| `collective_fit` | None | `--collective-fit` |

---

## 典型处理流程

```python
import numpy as np
from astropy.io import fits
from fitsh_cy import *

# 1. 星表匹配 → 获取变换
ref = np.loadtxt('ref.cat')
obs = np.loadtxt('obs.cat')
g = Grmatch(order=2)
g.set_reference(ref[:,1], ref[:,2], -ref[:,3])
g.set_input(obs[:,1], obs[:,2], -obs[:,3])
r = g.run()
print(f'Matched {r.nhit} pairs')

# 2. 星检测
img = fits.getdata('obs.fits').astype(np.float64)
sr = Fistar(threshold=100.0, algorithm='uplink').search(img)
print(f'Found {sr["nstar"]} stars')

# 3. 图像变换 (几何校正)
t = Fitrans.ImageTransformation(transformation=r.transform, method='lanczos3')
corrected, _ = t.apply(img)

# 4. 核拟合 → 卷积匹配
fc = Ficonv(kernel='i/0; b/0; g=5,2.0,2/0', iterations=3)
cnv, mask = fc.fit(corrected, img)

# 5. 平滑
smooth, _ = Fitrans.Smooth('spline', order=2).apply(cnv)
```
