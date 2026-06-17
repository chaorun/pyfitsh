# pyfitsh 用法

## Fiarith — 逐像素算术表达式求值器

```python
from pyfitsh.fiarith import Fiarith
import numpy as np

fa = Fiarith()
result = fa.evaluate(expr, operands, mask=None)
```

**参数说明**：

| Python 参数 | CLI 参数 | 说明 |
|-------------|----------|------|
| `expr` | `<expression>` | fiarith 算术表达式字符串 |
| `operands` | 表达式内嵌 `'file.fits'` | dict，key=变量名 value=numpy 2D float64 数组 |
| `mask` | `--input-mask` | 可选的 mask 数组，uint8 |

**表达式语法**：

| 类别 | 函数 | 说明 |
|------|------|------|
| 算术 | `+ - * /` | 加/减/乘/除（除零自动得零） |
| 通用函数 | `sin cos tan asin acos atan atan2 arg` | 三角函数 |
| 通用函数 | `abs sq sqrt exp log` | 绝对值/平方/平方根/指数/对数 |
| 图像函数 | `min(...) max(...) mean(...)` | 逐像素极值/均值 |
| 图像函数 | `norm(x) sign(x) theta(x)` | 均值/符号/阶跃 |
| 图像函数 | `laplace(x) scatter(x)` | 拉普拉斯变换/噪声估计 |
| 图像函数 | `corr(a,b)` | 两图像互相关 |
| 图像函数 | `smooth(img,sigma,size)` | 高斯平滑 |

**示例**：

```python
# 简单运算
result = fa.evaluate('a+b', {'a': img_a, 'b': img_b})
# 图像函数
result = fa.evaluate('smooth(a, 2.0, 5)', {'a': img_a})
# 子表达式 (UDF): body 仅支持通用函数 (sin/cos/exp/sqrt...)，不支持图像函数
result = fa.evaluate('[a](sqrt(a))', {'a': img_a})
```

原版限制：UDF `[...](body)` 不支持图像函数（sq/smooth/laplace 等），与 CLI 一致。

---

## Fistar — 星检测与拟合

```python
from pyfitsh import Fistar
from astropy.io import fits
import numpy as np

img = fits.getdata('image.fits').astype(np.float64)
fs = Fistar(threshold=100.0, algorithm='uplink', model='elliptic',
            model_order=2, it_sym=4, it_gen=2, gain=1.0,
            sort='x', fields='x,y,bg,amp,fwhm,ellip,flux,s/n,magnitude',
            verbose=False, output_mark=False, output_area=False)
r = fs.do_fistar(img)
tab = r.output
print(f'nstar={r.nstar}')
```

| Python 参数 | CLI 参数 | 说明 |
|-------------|----------|------|
| `threshold` | `-t` | 检测阈值 (default 100.0) |
| `flux_threshold` | `-f` | 通量阈值，>0 时替代 threshold 用于 UPLINK |
| `skysigma` | `-d` | 天空背景噪声 |
| `gain` | `-g` | 探测器增益 (e-/ADU) |
| `algorithm='uplink'\|'parabolapeak'` | `--algorithm` | 候选检测算法 |
| `model='elliptic'\|'gauss'\|'deviated'` | `--model` | 星点模型 |
| `only_candidates=True` | `--only-candidates` | 仅候选检测，跳过拟合 |
| `psf='native,order=2,grid=4'` | `--psf` | PSF 拟合参数 |
| `sort='x'\|'y'\|'peak'\|'flux'\|'fwhm'` | `-s` | 排序方式 |
| `section=(x1,x2,y1,y2)` | `--section` | 图像区域限制 |
| `output_mark=True` | `--mark-output` | 输出标记图像 |
| `output_area=True` | `--output-area` | 输出区域图像 |
| `it_sym` | `--it-sym` | 对称拟合迭代次数 |
| `it_gen` | `--it-gen` | 通用拟合迭代次数 |
| `model_order` | `--model-order` | 模型阶数 (0/1/2) |
| `verbose=True` | `--verbose` | 详细输出 |

**高斯模型 + PSF 拟合示例**（对应 CLI `-d 5 -f 10000 --model gauss --psf native,order=2`）：

```python
fs = Fistar(threshold=100.0, flux_threshold=10000.0, skysigma=5.0,
            algorithm='uplink', model='gauss', model_order=2,
            it_sym=4, it_gen=2, gain=2.0,
            sort='x', psf='native,order=2', psf_output=True,
            verbose=True)
r = fs.do_fistar(img_data, mask=mask_data)
# r.output (Table), r.nstar (int), r.psf (astropy HDU)
```

> **已知问题**：gauss 模型下 `l`、`sigma`、`delta`、`kappa`、`fwhm`、`pa`、`flux`、`noise`、`sn` 列与 CLI 输出有差异。根因疑为 `star-model.c` 中 `shape.gl` 赋值逻辑——该文件与 origincode 一致，但 CLI 在 PSF 拟合后更新了 `shape.gl`，Cython 端未同步。`x`、`y`、`bg`、`amp`、`s`、`d`、`k` 完全一致。（挂起）

---

## Fiphot — 孔径测光

### 基本用法

```python
from pyfitsh import Fiphot
import numpy as np

stars = np.loadtxt('positions.dat')
fp = Fiphot(apertures='10:18:12', gain='2', mag_flux=(10.0, 10000.0))
r = fp.photometry(img, stars=stars, col_xy=(1, 2), col_id=0)
```

多孔径：`apertures='8:18:12,10:18:12,12:20:15'`（逗号分隔）。

| Python 参数 | CLI 参数 | 说明 |
|-------------|----------|------|
| `apertures` | `-a` | 孔径规格 `rad:r_in:r_out` |
| `gain` | `-g` | 增益 (e-/ADU) |
| `mag_flux` | `--mag-flux` | (mag, flux) 零点校准 |
| `zero` | `--zero` | 星等零点值 |
| `dark` | `--dark` | 暗场值 |
| `flat` | `--flat` | 平场值 |
| `sort` | `-s` | 排序方式 |
| `output_raw=True` | `--output-raw` | 输出 raw 测光数据 |
| `cols_xy` | (N/A) | 输入星表中 x,y 列索引 |
| `cols_id` | (N/A) | 输入星表中 ID 列索引 |
| `mask` | `-K` | 蒙版 |

### photometry_from_raw

对应 CLI `--input-raw-photometry`。第一个参数应传 3D raw 数组 `(nstar, nap, 12)`，而非 output Table：

```python
r = fp.photometry(img, stars=stars, col_xy=(1,2), col_id=0)
raw_3d = r.dict['output_raw_photometry_3d_data']  # shape (nstar, nap, 12)
r2 = fp.photometry_from_raw(raw_3d, subtracted_img)
# 可选：kernel 减影模式
r2 = fp.photometry_from_raw(raw_3d, subtracted_img,
                            kernel_spec='gauss:3:3', normalize_kernel=True)
```

### subtracted_photometry / magfit

```python
r = fp.subtracted_photometry(img, subtracted_img, stars=stars)
r = fp.magfit(mag_arr, flux_arr)
```

---

## Grmatch — 星点匹配

```python
from pyfitsh import Grmatch

ref = np.loadtxt('ref.cat')
inp = np.loadtxt('obs.cat')
g = Grmatch(order=2, maxdist=2.0, unitarity=0.2, ttype=2,
            use_ordering=1, is_centering=1, nmiter=4, rejlevel=3.0)
r = g.matchpoints(ref[:,1], ref[:,2], inp[:,1], inp[:,2],
                  ref_order=-ref[:,3], inp_order=-inp[:,3])
print(f'nhit={r.nhit}, trans={r.transformation}')
```

| Python 参数 | CLI 参数 | 说明 |
|-------------|----------|------|
| `order` | `-a` | 多项式变换阶数 |
| `maxdist` | `--max-distance` | 匹配点最大接受距离 |
| `unitarity` | `--unitarity` | 对称匹配阈值 |
| `ttype` | `--triangulation` | 三角剖分类型 |
| `use_ordering` | `--use-ordering` | 是否使用星等排序 |
| `is_centering` | `--centering` | 是否对坐标做中心化 |
| `nmiter` | `--nm-iter` | 三角匹配迭代次数 |
| `rejlevel` | `-r` | 拒绝 sigma 级别 |
| `min_matches` | `--min-matches` | 最少匹配点数 |

**坐标输入**：

```python
# 基本用法（无星等信息）
r = g.matchpoints(ref_x, ref_y, inp_x, inp_y)

# 带星等排序
r = g.matchpoints(ref_x, ref_y, inp_x, inp_y,
                  ref_order=ref_mag, inp_order=inp_mag)

# 带 ID 匹配
r = g.matchpoints(ref_x, ref_y, inp_x, inp_y,
                  ref_id=ref_id_arr, inp_id=inp_id_arr)

# 带 PSF 权重
r = g.matchpoints(ref_x, ref_y, inp_x, inp_y,
                  ref_weight=ref_w, inp_weight=inp_w)
```

**返回**：
- `r.nhit` — 匹配星点数
- `r.transformation` — 变换参数字典
- `r.matched` — 匹配后的坐标数组
- `r.excluded_ref` / `r.excluded_inp` — 被排除的点
- `r.dict` — 向后兼容字典

---

## Grtrans — 坐标变换拟合 / 评估

```python
from pyfitsh import Grtrans

# 拟合变换
g = Grtrans(order=2, nmiter=5, rejlevel=3.0)
r = g.fit_transform(x_data, y_data, fit1_data, fit2_data,
                    weight=w_data)
print(r.transformation)

# 应用变换
trans = r.transformation
result = g.evaluate(x_data, y_data, transformation=trans)
```

| Python 参数 | CLI 参数 | 说明 |
|-------------|----------|------|
| `order` | `-a` | 多项式阶数 |
| `nmiter` | `-n` | 拒绝迭代次数 |
| `rejlevel` | `-r` | 拒绝 sigma 级别 |
| `transformation` | `-T` | 输入变换（用于评估模式） |

**返回**：
- `r.transformation` — 变换参数字典，direct/inverse
- `r.output` — 变换后的坐标点

**变换字典结构**（与 fitrans 共用）：

```python
trans = {
    'type': 'polynomial',     # 变换类型
    'order': 2,               # 多项式阶数
    'offset': (0.0, 0.0),     # 坐标偏移
    'scale': 1.0,             # 坐标缩放
    'dxfit': [...],           # X 方向多项式系数
    'dyfit': [...],           # Y 方向多项式系数
}
```

---

## Fitrans — 图像变换

```python
from pyfitsh import Fitrans

trans = {'type':'polynomial', 'order':1, 'offset':(0.0,0.0), 'scale':1.0,
         'dxfit':[-148.653, 0.99898, 0.04297],
         'dyfit':[140.218, -0.04291, 0.99912]}
ft = Fitrans.ImageTransformation(transformation=trans, inverse=True, method='bilinear')
output, mask = ft.apply(img)
```

| Python 参数 | CLI 参数 | 说明 |
|-------------|----------|------|
| `transformation` | `-t` / `-T` | 变换参数字典或文件路径 |
| `inverse=True` | `--reverse` | 使用逆变换 |
| `method='bilinear'` | `-l` / `-c` / `-m` / `-k` | 插值方法 |
| `size=(sx,sy)` | `-s` | 输出图像尺寸 |
| `offset=(x,y)` | `-f` | 输出图像零点坐标 |

插值方法对照：

| Python | CLI | 说明 |
|--------|-----|------|
| `'simple'` | `-m` | 简单线性插值（不保通量） |
| `'bilinear'` | `-l` | 线性插值（保通量） |
| `'bicubic'` | `-c` | 双三次样条插值 |
| `'biquad'` | `-k` | 双二次插值（保通量） |

---

## Ficonv — 卷积/减影

```python
from pyfitsh.ficonv.ficonv import Ficonv
from pyfitsh.utils import decode_maskinfo
from astropy.io import fits

ref = fits.getdata('ref.fits').astype(np.float64)
img = fits.getdata('v1.fits').astype(np.float64)

with fits.open('v1.fits') as hdul:
    img_mask = decode_maskinfo(hdul[0])

fc = Ficonv(kernel='i/0;b/1;g=3,1.0,2/0', iterations=3, rejection_level=3.0)
convolved, mask, subtracted, kernel_dict = fc.fit(ref, img, img_mask=img_mask,
                                                   output_subtracted=True,
                                                   output_kernel_list=True)
```

| Python 参数 | CLI 参数 | 说明 |
|-------------|----------|------|
| `kernel` | `-k` | kernel 规格字符串 |
| `iterations` | `-n` | 拒绝迭代次数 |
| `rejection_level` | `-s` | 拒绝 sigma |
| `masked` | `-m` | 蒙版拟合模式 |
| `weighted` | `-w` | 加权拟合（隐含 -m） |
| `background_iterative` | `-b` | 背景迭代拟合 |
| `divide` | `-d` | 分块因子 (default 32) |
| `gain` | `-g` | 增益 (e-/ADU) |
| `unity_kernels` | `-u` | 归一化 identity kernel |

**用 kernel_dict 跳过拟合**（等价 CLI `--input-kernel-list`）：

```python
convolved, mask = fc.fit(ref, img, kernel_dict=kernel_dict)
```

`kernel_dict` 结构：

```python
{
    'global': {'type': 1, 'offset': [1024.0, 1024.0], 'scale': 1024.0},
    'kernels': [
        {'index': 0, 'type': 1, 'order': 1, 'coeff': [253.3, -5.98, -23.42]},
        {'index': 1, 'type': 2, 'order': 0, 'coeff': [0.939]},
        {'index': 2, 'type': 3, 'order': 0, 'coeff': [2.132], 'hsize': 3, 'sigma': 1.0},
    ]
}
```

类型码：1=background, 2=identity, 3=gaussian。

**Mask 处理注意事项**：
- CLI 从 FITS header 自动读取 MASKINFO；Cython 版需显式调用 `decode_maskinfo(hdu)`
- ref 图像通常无 MASKINFO，传入全零 mask
- NaN 像素需在 Python 层标记为 mask：`mask[np.isnan(data)] = 1`

---

## Firandom — 人工图像生成

```python
from pyfitsh import Firandom
import numpy as np

fir = Firandom(sx=512, sy=512, gain=2.0, seed=12345,
               background='1000', is_photnoise=True, zoom=1)
result = fir.generate(
    stars='20[x=r(0.1,0.9),y=r(0.1,0.9),f=6,e=0.05,p=r(0,180),i=r(50000,200000)]',
    background_stddev=5.0)
img = result.image
```

| Python 参数 | CLI 参数 | 说明 |
|-------------|----------|------|
| `sx, sy` | `--size` | 图像尺寸 |
| `gain` | `--gain` | 增益 (e-/ADU) |
| `background` | `--sky` | 背景表达式 |
| `background_stddev` | `--sky-noise` | 背景额外噪声 |
| `is_photnoise` | `--photon-noise` | 光子噪声仿真 |
| `seed` | `--seed` | 通用随机种子（fallback for seed_noise/seed_spatial/seed_photon） |
| `seed_noise` | `--seed-noise` | 背景噪声种子 |
| `seed_spatial` | `--seed-spatial` | 空间坐标种子 |
| `seed_photon` | `--seed-photon` | 光子噪声种子 |
| `dontquantize` | (默认) | 默认 True，匹配 CLI |
| `method=1` | `--integral` | 积分绘制模式 |
| `zoom` | `--zoom` | 过采样因子 |

**Star 输入格式**：

1. firandom 表达式字符串（对应 CLI `--list`）：

```python
fir.generate(stars='5[x=r(0.2,0.8),y=r(0.2,0.8),f=4,e=0.1,p=r(0,180),i=50000]')
```

2. Python list-of-dicts：

```python
stars = [
    {'x': 128, 'y': 128, 'flux': 10000, 'fwhm': 3.0, 'ellip': 0.1, 'pa': 45},
    {'x': 256, 'y': 256, 'flux': 5000,  's': 1.5, 'd': 0.0, 'k': 0.0},
]
fir.generate(stars=stars)
```

每个 dict 支持：`x, y, flux, s/d/k`（sigma/delta/kappa）或 `fwhm/ellip/pa`。

**独立绘制星点**：

```python
img = np.zeros((512, 512), dtype=np.float64)
fir.draw_stars(img, stars)
```

---

## 通用注意事项

- 所有模块返回 `types.SimpleNamespace`，用 `r.output` 访问主结果
- 图像数据均为 float64 numpy 数组，C-order (row-major)
- mask 为 uint8 numpy 数组，0 = 有效像素
- `decode_maskinfo(hdu)` 读取 FITS 头中 MASKINFO 关键字生成 mask
- UDF 子表达式 `[...](body)` 仅支持通用函数，不支持图像函数（与 CLI 一致）
