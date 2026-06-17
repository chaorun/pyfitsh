# fitsh_cython — Python 原生天文图像处理扩展 v0.02

基于 [fitsh](https://fitsh.net/) 的 C 算法库，编译为 Python 原生 `.so` 扩展，
覆盖全部 5 个 CLI 程序：grmatch, grtrans, ficonv, fitrans, fistar。

已通过原始 CLI 二进制交叉验证（39 项测试全部通过）。

## 构建

```bash
cd fitsh_cython
pip install cython numpy astropy
python3 setup.py build_ext --inplace
```

## 快速导入

```python
from fitsh_cy import (Grmatch, MatchResult, Grtrans, Ficonv,
                      PolyFitter, Fitrans, Fistar,
                      parse_trans_file)
```

## 模块一览

| 模块 | 原始 CLI | 功能 |
|------|---------|------|
| `Grmatch` | grmatch | 星表匹配 + 三角剖分 |
| `Grtrans` | grtrans | 多项式坐标变换 + WCS 拟合 |
| `Ficonv` | ficonv | 核拟合 + 图像卷积 |
| `Fitrans` | fitrans | 图像变换 (Zoom/Shrink/Repetitive/Interleave/Noise/Smooth/几何变换) |
| `Fistar` | fistar | 星检测 + PSF 拟合 + 测光 |

---

## Grmatch — 星表匹配

```python
m = Grmatch(order=1, maxdist=2, unitarity=0.2, verbose=True)
m.matchxy(ref_x, ref_y, obj_x, obj_y)
# 带星等排序
m.matchxy(ref_x, ref_y, obj_x, obj_y,
          ref_ord=-ref_mag, inp_ord=-obj_mag)
r = m.run()

r.nhit          # 匹配对数
r.hits_ref      # ref 表匹配行索引 (0-based)
r.hits_inp      # obj 表匹配行索引
r.vfits_dx      # X 方向变换多项式系数
r.vfits_dy      # Y 方向变换多项式系数
r.transform     # dict (可直接传入 Grtrans)
r.stats         # {'nsigma': 0.59, 'unitarity': 0.0001, ...}
```

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `order` | 2 | 多项式阶数 |
| `maxdist` | 2.0 | 最大匹配距离 (像素) |
| `unitarity` | 0.01 | 三角剖分单位性阈值 |
| `ttype` | 3 | 三角剖分模式: 0=delaunay, 1=expanded, 2=huber, 3=auto |
| `parity` | 0 | 奇偶约束: 0=忽略, 1=保持, -1=翻转 |
| `use_ordering` | 1 | 按排序值选点 |
| `nmiter` | 0 | 自洽迭代次数 |
| `rejlevel` | 3.0 | 迭代剔点阈值 (sigma) |
| `verbose` | False | 打印匹配进度与统计 |

---

## Grtrans — 坐标变换

```python
# 从匹配结果构造
t = Grtrans(r.transform)
new_x, new_y = t.apply(pixel_x, pixel_y)

# 手动构造
t = Grtrans(order=1, dxfit=[0,1,0], dyfit=[0,0,1], invert=True)

# 从 .trans 文件 / 字符串
t = Grtrans(parse_trans_file('grb.trans'))
t = Grtrans.from_trans_string(trans_str)
trans_str = t.to_trans_string()
```

### WCSFitter

```python
w = Grtrans.WCSFitter(ra0=189.0, dec0=62.0, projection='tan', order=3)
w.fit(ra_list, dec_list, x_pix, y_pix)
xi, eta = w.project(ra, dec)    # (RA,Dec) -> 投影平面
x, y = w.proj_to_pix(xi, eta)   # 投影平面 -> 像素
```

### PolyFitter

```python
f = PolyFitter(order=2, niter=3, rejlevel=3.0)
f.fit(x, y, values, weights=None)
print(f.coeff, f.residual, f.nrej)
```

---

## Ficonv — 核拟合 + 卷积

```python
fc = Ficonv(kernel='i/0; b/0; g=3,1.0,0/0',
            iterations=3, rejection_level=3.0,
            divide=32, gain=1.0, unity_kernels=False)

cnv, mask = fc.fit(ref, img)
cnv, mask = fc.fit(ref, img, ref_mask=rm, img_mask=im)
cnv, mask, subtracted = fc.fit(ref, img, output_subtracted=True)
cnv, mask, kernel_list = fc.fit(ref, img, output_kernel_list=True)
# stamp 区域拟合
cnv, mask = fc.fit(ref, img, stamps=[(x, y, sx, sy), ...])
```

kernel 格式 (同 ficonv CLI)：

```
i/<spatial order>                    identity kernel (通量项)
b/<spatial order>                    background offset
d=<hsize>/<spatial order>            discrete kernel
g=<hsize>,<sigma>,<order>/<spatial>  Gaussian kernel
```

示例: `'i/0;b/1'`, `'i/0;b/1;g=5,2.0,2/1'`, `'g=3,1.0,0/0;g=5,2.0,0/0;g=7,4.0,0/0'`

> 也可以通过 `Fitrans.Convolution(...)` 使用，接口完全一致。

---

## Fitrans — 图像变换

### decode_maskinfo

```python
from astropy.io import fits
mask = Fitrans.decode_maskinfo(fits.open('image.fits')[0])
# 等价于 C: fits_mask_read_from_header(), 返回 uint8 ndarray
```

### ImageTransformation — 几何变换

```python
f = Fitrans.ImageTransformation(method='lanczos3', transformation=r.transform)
out, mask = f.apply(img)
out, mask = f.apply(img, out_shape=(4096, 4096))
```

method: `bilinear`, `integrate`, `bicubic`, `spline_integrate`, `lanczos3`, `lanczos4`.

### Noise / Zoom / Shrink / Repetitive / Interleave / Smooth

```python
# Noise
out, mask = Fitrans.Noise().apply(img)

# Zoom (整数倍放大)
out, mask = Fitrans.Zoom(3).apply(img)           # 带双二次子像素插值
out, mask = Fitrans.Zoom(2, raw=True).apply(img)  # 无插值

# Shrink (整数倍缩并，流量守恒)
out, mask = Fitrans.Shrink(2).apply(img)                 # 均值
out, mask = Fitrans.Shrink(2, median=True).apply(img)    # 中值

# Repetitive (重复扩展 → --repetitive-xy)
out, mask = Fitrans.Repetitive(2, 3).apply(img)
out, mask = Fitrans.Repetitive(2, 3, offset_x=5).apply(img)

# Interleave (交错扩展 → --interleave-xy)
out, mask = Fitrans.Interleave(4, 6).apply(img)
out, mask = Fitrans.Interleave(4, 6, median=True).apply(img)

# Smooth (大尺度平滑)
out, mask = Fitrans.Smooth('spline', order=2).apply(img)
out, mask = Fitrans.Smooth('polynomial', order=3, prefilter='median', niter=3).apply(img)
```

---

## Fistar — 星检测与测光

```python
r = Fistar(
    threshold=100.0,
    flux_threshold=0.0,
    model='gauss',
    sort='x',
    fields='x,y,flux,magnitude,fwhm',
    output_mark=True,
    output_area=True,
    psf='native,order=2',
    psf_output=True,
).search(img_data, mask=mask)

r['table']         # astropy.table.Table
r['output_mark']   # ndarray (标记图像)
r['output_area']   # ndarray (星区域图像)
r['psf']           # fits.ImageHDU (PSF 数据立方)
r['nstar']         # int
```

| 参数 | 默认 | CLI 对应 |
|------|------|---------|
| `threshold` | 100.0 | `-t` |
| `algorithm` | `'uplink'` | `--algorithm` |
| `model` | `'elliptic'` | `--model` (gauss/elliptic) |
| `model_order` | 2 | `--model order=N` |
| `psf` | None | `--psf` (native/integral/circle + 子参数) |
| `gain` | 1.0 | `-g` |
| `sort` | `'x'` | `-s` (x/y/peak/fwhm/amp/flux/noise/sn) |
| `mag_flux` | `(10.0, 10000.0)` | `--mag-flux` |
| `input_candidates` | None | `-C` (numpy (N,2) 或 (N,3)) |
| `input_positions` | None | `-P` (list of (x,y) tuples) |
| `candidate_radius` | 2.0 | `-R` |

---

## 测试

```bash
cd fitsh_cython && python3 test_comprehensive.py
# 预期: RESULTS: 39/39 passed, 0 failed
```

---

## 模块结构

```
fitsh_cython/
  setup.py                构建配置
  fitsh_cy.pyx            Cython 接口 (~1800行)
  fitsh_cy.pxd            类型声明
  __init__.py             包导出
  _core.c / _trans_core.c / _wcs_core.c    gmatch/grtrans 核心
  _fitrans_core.c / _fitrans_ops.c         fitrans 核心 + 算子
  _polyfit_core.c                          多项式拟合核心
  ficonv_pipeline.c       卷积管线 (来自 ficonv.c)
  fistar_pipeline.c       星检测管线 (来自 fistar.c)
  fits_stubs.c            FITS 文件 IO stub
  kernel-base.c / kernel-io.c  核函数 + 解析
  star-*.c / psf-*.c      星检测 / PSF 算法 (逐字一致于原版)
  math/                   数学库 (完全一致于原版)
  test_comprehensive.py   完整测试 (39 项)
  deprecated/             已弃用文件 (备查)
```
