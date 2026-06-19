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
| `threshold` | `-t` | 检测阈值 (ADU, default 100.0) |
| `flux_threshold` | `-f` | 通量阈值，>0 时替代 threshold |
| `critical_prominence` | `-p` | UPLINK 算法的临界相对突出度 |
| `skysigma` | `-d` | 天空背景噪声 |
| `gain` | `-g` | 探测器增益 (e-/ADU) |
| `algorithm='uplink'\|'parabolapeak'` | `--algorithm` | 候选检测算法 |
| `model='elliptic'\|'gauss'\|'deviated'` | `--model` | 星点模型 |
| `model_order` | `--model-order` | 模型阶数 (deviated 时 2-4) |
| `it_sym` | `--iterations symmetric=` | 对称拟合迭代次数 (default 4) |
| `it_gen` | `--iterations general=` | 通用拟合迭代次数 (default 2) |
| `only_candidates=True` | `--only-candidates` | 仅候选检测，跳过模型拟合 |
| `psf='native,order=2,grid=4'` | `--psf` | PSF 参数：type,order,grid,halfsize,symmetrize,biquad,kappa,circlewidth,circleorder |
| `psf_output=True` | `--output-psf` | 输出 PSF 3D FITS 文件 |
| `sort='x'\|'y'\|'peak'\|'fwhm'\|'amp'\|'flux'\|'noise'\|'s/n'` | `-s` | 排序方式 |
| `fields='id,x,y,bg,amp,...'` | `--format` | 输出列（逗号分隔） |
| `section=(x1,x2,y1,y2)` | `--section` | 图像区域限制 |
| `input_candidates` | `-C\|--input-candidates` | 输入候选文件/数组 (N,2) 或 (N,5) |
| `candidate_radius` | `--candidate-radius` | 候选半径 (default 2.0) |
| `input_positions` | `--input-positions` | 输入位置匹配列表 |
| `mag_flux=(mag, flux)` | `--mag-flux` | 星等-通量校准 (default 10.0, 10000.0) |
| `output_mark=True` | `--mark-output` | 输出标记图像 |
| `output_area=True` | `--output-area` | 输出区域图像 |
| `mark_symbol='dot'` | `--mark-symbol` | 标记符号 |
| `mark_size` | `--mark-size` | 标记大小 |
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

> **已修复**：PSF type 常量偏移 bug（`_psf_type` 赋值 `native→0` 应为 `native→1`，`integral→1` 应为 `integral→2`，`circle→2` 应为 `circle→3`），该 bug 导致 `pbg`/`pamp` 列为 nan。（2026-06-17）
> 
> **已修复**：`fistar_core.c` input_candidates 路径三处修复：① 删除 s/d/k 默认值覆写 ② 增加 `cleanup_candlist` 调用 ③ 结果填充循环中 `wc` 指针加边界校验。修复后不再崩溃。（2026-06-17）
> 
> **已知 CLI bug**：`--input-candidates --col-shape 3,4,5` 时 `read_star_candidates` 未读取 cd/ck 列（默认为 0），导致 `make_star_candidates` 中 sxx=syy=cs, sxy=0，形状参数丢失，`peak`/`amp`/`flux` 为 0。这是 origincode CLI 的 bug，非 Cython 端问题。Cython 端通过 `input_candidates` 传入正确 cs/cd/ck 值可获得准确结果。测试脚本中为匹配 CLI 行为，将 cd/ck 列置零。（2026-06-18）
> 
> **已知差异**：test_04（ficonv from fitted kernel）中卷积图像有 2/4,194,304 像素差异，max_diff=0.107。原因为 `kernel_dict` 在 Python 层 round-trip 后产生浮点精度损失，非算法问题。（2026-06-18）

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
| `apertures` | `-a` | 孔径规格 `rad:r_in:r_out`，逗号分隔 |
| `gain` | `-g` | 增益多项式 (e-/ADU)，支持空间变化 |
| `gain_vmin` | `--gain-vmin` | 最小增益值 |
| `mag_flux` | `--mag-flux` | (mag, flux) 零点校准 |
| `sky_fit` | `--sky-fit` | 天空拟合: mean\|median\|mode,iterations=N,sigma=S |
| `correlation_length` | `--correlation-length` | 背景相关长度 (default 1.0) |
| `zoom` | `-z` | 坐标/孔径缩放因子 |
| `serial` | `--serial` | 序列标识符 |
| `mask_ignore` | `--aperture-mask-ignore` | 孔径内忽略的 mask 类型: saturated,hot,outer |
| `nan='-'` | `--nan` | 坏测光标记字符串 |
| `output_raw=True` | `--output-raw-photometry` | 输出 raw 测光数据 |
| `cols_xy` | `--col-xy` | 输入星表中 x,y 列索引 (1-indexed) |
| `cols_id` | `--col-id` | 输入星表中 ID 列索引 |
| `cols_ap` | `--col-ap` | 输入星表中每星自定义孔径列 |
| `cols_mag` | `--col-mag` | 参考星等列 |
| `cols_col` | `--col-col` | 测光颜色列 |
| `cols_err` | `--col-err` | 星等误差列 |
| `mask` | `-M` | 蒙版图像 |

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
| `zoom` | `-z` | 整数倍放大（双二次插值，保通量） |
| `shrink` | `-r` | 整数倍缩小 |
| `median_shrink` | `-d` | 缩小使用中值平均 |
| `smooth` | `-a` | 图像平滑参数 |
| `noise` | `-n` | 噪声估计图 |
| `layer` | `-y` | 数据立方体切层 |
| `explode` | `-x` | 数据立方体展开 |
| `bitpix` | `-b` | 输出 FITS bitpix |

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
| `kernel` | `-k` | kernel 规格字符串，分号分隔 |
| `iterations` | `-n` | 拒绝迭代次数 |
| `rejection_level` | `-s` | 拒绝 sigma |
| `masked` | `-m` | 蒙版拟合模式 |
| `weighted` | `-w` | 加权拟合（隐含 -m） |
| `background_iterative` | `-b` | 背景迭代拟合 |
| `divide` | `-d` | 分块因子 (default 32) |
| `gain` | `-g` | 增益 (e-/ADU) |
| `unity_kernels` | `-u` | 归一化 identity kernel |
| `verbose` | `--verbose` | 详细输出 |

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
| `background` | `--sky` | 背景表达式（支持 x,y 函数） |
| `background_stddev` | `--sky-noise` | 背景额外噪声 |
| `is_photnoise` | `--photon-noise` | 光子噪声仿真 |
| `seed` | `--seed` | 通用随机种子 |
| `seed_noise` | `--seed-noise` | 背景噪声种子 |
| `seed_spatial` | `--seed-spatial` | 空间坐标种子 |
| `seed_photon` | `--seed-photon` | 光子噪声种子 |
| `zoom` | `--zoom` | 过采样因子 |
| `method=1` | `--integral` | 积分绘制模式 (default) |
| `method=0` | `--monte-carlo` | Monte-Carlo 绘制模式 |
| `dontquantize=True` | -- | 默认不量化 (匹配 CLI `--no-quantize`) |
| `nsuppress=10000.0` | `--noise-suppression` | 积分边界抑制级别 |

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

## Ficalib — 图像校准（bias/flat/dark）

```python
from pyfitsh import Ficalib
import numpy as np
from astropy.io import fits

fc = Ficalib()
r = fc.calibrate(sci, bias=bias_data, flat=flat_data,
                 dark=dark_data, flat_mean=1.0, dark_time=0.0)
calibrated = r['data']
```

| Python 参数 | CLI 参数 | 说明 |
|-------------|----------|------|
| `img_data` | `-i` | 科学图像 float64 2D |
| `mask` | — | 可选 mask uint8 2D |
| `bias` | `-B` | bias 帧 float64 2D |
| `dark` | `-D` | dark 帧 float64 2D |
| `flat` | `-F` | flat 帧 float64 2D |
| `flat_mean` | `--flat-mean` | flat 归一化因子 (default 1.0) |
| `dark_time` | `--dark-time` | dark 曝光时间比例 (0 = 使用 1.0) |

**返回**：`{'data': ndarray float64, 'mask': ndarray uint8}`

---

## Fiign — 像素 mask 操作

```python
from pyfitsh import Fiign
import numpy as np

fg = Fiign()
r = fg.apply(img, mask=mask,
             ignore_nonpos=True, saturation=50000,
             ignore_cosmics=True, replace_cosmics=True,
             th_low=5.0, th_high=10.0,
             apply_mask=True, mask_value=-999)
# r['data'], r['mask'], r['history']
```

| Python 参数 | CLI 参数 | 说明 |
|-------------|----------|------|
| `img` | `-i` | 输入图像 float64 2D |
| `mask` | `-M` | 可选输入 mask uint8 2D（多个外部 mask 需 Python 层 `\|=` 预合并） |
| `saturation` | `-s` | 饱和阈值 (>0 启用) |
| `saturation_img` | `-S` | 逐像素饱和阈值 float64 2D |
| `leak_method` | `--lu`/`--lr`/`--an` | blooming 方向：0=none, 1=垂直, 2=水平, 3=任意 |
| `ignore_nonpos` | `-n` | 标记非正像素 (≤0) 为 MASK_FAULT |
| `ignore_neg` | `-g` | 标记负像素为 MASK_FAULT |
| `ignore_zero` | `-z` | 标记零值为 MASK_FAULT |
| `ignore_cosmics` | `-c` | 启用宇宙线检测 |
| `replace_cosmics` | `-r` | 用局部均值替换宇宙线像素 |
| `th_low` | `--threshold-low` | 宇宙线低阈值 (default 10.0) |
| `th_high` | `--threshold-high` | 宇宙线高阈值 (default 50.0) |
| `sky_sigma` | `-d` | 天空背景噪声 (>0 时限制 sig≤2×skysigma) |
| `expand_hsize` | `-x` | mask 膨胀半径 (default 0) |
| `apply_mask` | `-a` | 将被 mask 像素替换为 `mask_value` |
| `mask_value` | `-m` | 被 mask 像素的替换值 (default 0.0) |
| `bitpix` | `-b` | 整数 bitpix (8/16/32)，负数=浮点 (default -32) |
| `convert_list` | `--convert` | mask 转换规格列表，如 `['fault:any:none:cosmic']` |
| `mask_block_list` | `-q` | 几何形状列表，如 `['block:hot:100,100:300,300', 'circle:cosmic:500,500:50']` |

**mask block 格式**：
- `block:<mask>:<x1>,<y1>:<x2>,<y2>` — 矩形
- `circle:<mask>:<xc>,<yc>:<radius>` — 圆形
- `pixel:<mask>:<x>,<y>` — 单像素
- `line:<mask>:<x1>,<y1>:<x2>,<y2>[:<width>]` — 线段

**mask 名称**：`none`, `clear`, `fault`, `hot`, `cosmic`, `outer`, `oversaturated`, `leaked`, `bloomed`, `saturated`, `interpolated`, `all`, `bad`

**返回**：`{'data': ndarray float64, 'mask': ndarray uint8, 'history': str}`

---

## 通用注意事项

- Fiign/Ficalib 返回普通 `dict`，其余模块返回 `types.SimpleNamespace`
- 图像数据均为 float64 numpy 数组，C-order (row-major)
- mask 为 uint8 numpy 数组，0 = 有效像素
- `decode_maskinfo(hdu)` 读取 FITS 头中 MASKINFO 关键字生成 mask
- UDF 子表达式 `[...](body)` 仅支持通用函数，不支持图像函数（与 CLI 一致）

---

## Lfit — 通用曲线拟合

完整复刻 origincode `lfit` 命令行工具（除文件 I/O 外），参数名对标 `lfit --long-help`。

### 基本用法

```python
from pyfitsh.lfit import Lfit, lfit_fit
import numpy as np

data = np.loadtxt('data.dat')

# 类式 API
lf = Lfit(method='clls', rejection_niter=3, rejection_level=3.0)
r = lf.fit(data, variables='a,b', columns='x:1,y:2',
           function='a+b*x', dependent='y')

# 函数式 API
r = lfit_fit(data, variables='a,b', columns='x:1,y:2',
             function='a+b*x', dependent='y')

# 结果访问 (属性模式 + dict 模式均可)
r.params          # array([1.99992, 3.0])
r['param_dict']   # {'a': 1.99992, 'b': 3.0}
r.chi2            # 6.9e-06
r['nrow']         # 77
r.ok              # True
r.keys()          # 所有可用 key
```

### 回归分析方法 (9种)

| Python `method=` | CLI 参数 | 说明 | 对应 `--long-help` |
|------------------|----------|------|--------------------|
| `'clls'` | `-L, --clls, --linear` | 经典线性最小二乘 | 模型函数须可微且线性 |
| `'nllm'` | `-N, --nllm, --nonlinear` | 非线性 Levenberg-Marquardt | 模型函数须可微，需初始值 `a=1,b=1` |
| `'lmnd'` | `-U, --lmnd` | Levenberg-Marquardt (数值导数) | 同 NLLM，但偏导数用数值近似，需指定 `differences` |
| `'dhsx'` | `-D, --dhsx, --downhill` | Downhill simplex | 需初始值+不确定度 `a=1:0.5` 或 `parameters='fisher'` |
| `'mcmc'` | `-M, --mcmc` | Markov Chain Monte-Carlo | 需初始值+不确定度，不要求可微 |
| `'mchi'` | `-K, --mchi, --chi2` | Chi2 网格搜索 | 变量用网格语法 `a=[0:1:3]` |
| `'emce'` | `-E, --emce` | Monte-Carlo 误差估计 | 需指定主方法 `parameters='clls'` 等 |
| `'xmmc'` | `-X, --xmmc` | 扩展 MCMC | Fisher 协方差转移分布，含 DHSX 初始化 |
| `'fima'` | `-A, --fima` | Fisher 信息矩阵分析 | 初始值应为"最佳拟合"值 |

### 完整参数对照表 (严格对标 `lfit --long-help`)

#### General options

| CLI 参数 | Python | 说明 |
|----------|--------|------|
| `-h, --help` | — | 信息输出，不适用 |
| `--long-help, --help-long` | — | 信息输出，不适用 |
| `--wiki-help, --help-wiki` | — | 信息输出，不适用 |
| `--version, --version-short` | — | 信息输出，不适用 |
| `--functions, --list-functions` | — | 信息输出，不适用 |
| `--wiki-functions` | — | 信息输出，不适用 |
| `--examples` | — | 信息输出，不适用 |

#### Common options for regression analysis

| CLI 参数 | Python | 类型 | 默认值 | 说明 |
|----------|--------|------|--------|------|
| `-v, --variable, --variables <list>` | `variables` | str | 必填 | 拟合变量定义。格式：`'a,b'` / `'a=1,b=1,c=0.1'`（含初始值）/ `'a=1:0.5,b=1:0.5'`（含步长）/ `'a,b,c:=0.5'`（含约束） |
| `-c, --column, --columns <cols>` | `columns` | str | 必填 | 输入列定义。格式：`'x:1,y:2,e:3'`，列号 1-indexed |
| `-f, --function <func>` | `function` | str | 必填 | 模型函数表达式。可含内建函数、宏、拟合变量、列变量 |
| `-y, --dependent <expr>` | `dependent` | str/None | None | 因变量表达式。None 时为求值模式（无拟合） |
| `-o, --output <file>` | — | — | — | 文件 I/O，不适用（结果通过返回值获取） |

#### Common options for function evaluation

| CLI 参数 | Python | 说明 |
|----------|--------|------|
| `-f, --function <func>[,...]` | `function` | 同上，求值模式时可逗号分隔多个表达式 |
| `-o, --output <file>` | — | 文件 I/O，不适用 |

注：`dependent=None` 时自动进入求值模式。

#### Regression analysis methods

| CLI 参数 | Python `method=` | 说明 |
|----------|------------------|------|
| `-L, --clls, --linear` | `'clls'` | 经典线性最小二乘。模型函数须可微且线性 |
| `-N, --nllm, --nonlinear` | `'nllm'` | 非线性 Levenberg-Marquardt。模型函数须可微，需初始值 |
| `-U, --lmnd` | `'lmnd'` | Levenberg-Marquardt 数值导数。偏导数数值近似，需 `differences` |
| `-D, --dhsx, --downhill` | `'dhsx'` | Downhill simplex。需初始值+步长 `a=1:0.5`，或 `parameters='fisher'` |
| `-M, --mcmc` | `'mcmc'` | Markov Chain Monte-Carlo。需初始值+步长，不要求可微 |
| `-K, --mchi, --chi2` | `'mchi'` | Chi2 网格搜索。变量用网格语法 `a=[0:1:3]`，固定变量 `c=0.5` |
| `-E, --emce` | `'emce'` | Monte-Carlo 误差估计。需指定主方法 `parameters='clls'` 等 |
| `-X, --xmmc` | `'xmmc'` | 扩展 MCMC。Fisher 协方差转移分布，含 DHSX 初始化 |
| `-A, --fima` | `'fima'` | Fisher 信息矩阵分析。初始值应为最佳拟合值 |

#### Fine-tuning of regression analysis methods

| CLI 参数 | Python | 类型 | 默认值 | 说明 |
|----------|--------|------|--------|------|
| `-e, --error <expr>` | `error` | str/None | None | 误差/sigma 表达式。0 或负值 → 权重为零 → 丢弃该行 |
| `-w, --weight <expr>` | `weight` | str/None | None | 权重表达式（= 1/sigma）。与 `error` 互斥 |
| `-P, --parameters <params>` | `parameters` | str/None | None | 方法微调参数（逗号分隔字符串，见下表） |
| `-q, --difference <diffs>` | `differences` | str/None | None | LMND 数值偏导数差分步长，如 `'a=0.001,b=0.001'` |
| `-k, --separate <vars>` | `separate` | str/None | None | 非线性方法中线性子集变量名（逗号分隔） |
| `--perturbations <noise>` | `perturbations` | str/None | None | EMCE 合成数据集附加白噪声级别，如 `'0.01'` |

#### `-P, --parameters` 支持的全部关键字

| 关键字 | 适用方法 | 说明 |
|--------|----------|------|
| `default`, `defaults` | 所有 | 使用默认微调参数 |
| `clls`, `linear` | EMCE | EMCE 主方法：经典线性最小二乘 |
| `nllm`, `nonlinear` | EMCE | EMCE 主方法：非线性 LM |
| `lmnd` | EMCE | EMCE 主方法：LM 数值导数 |
| `dhsx`, `downhill` | EMCE | EMCE 主方法：downhill simplex |
| `mc`, `montecarlo` | EMCE / FIMA | EMCE 主方法：原始 MC 扩散；FIMA：生成 mock 高斯分布 |
| `fisher` | DHSX / EMCE | 从 Fisher 协方差导出初始 simplex 大小 |
| `skip` | EMCE / XMMC | 跳过初始最小化 |
| `lambda=<value>` | NLLM / LMND | LM 算法 lambda 初始值 (default 0.001) |
| `multiply=<value>` | NLLM / LMND | LM 算法 lambda 乘数 (default 10.0) |
| `iterations=<max>` | NLLM / LMND / XMMC | LM 最大迭代数 (default 10)；XMMC 附加迭代数 |
| `accepted` | MCMC / XMMC | 计数已接受的转移（默认） |
| `nonaccepted` | MCMC / XMMC | 计数全部转移（含未接受） |
| `gibbs` | MCMC | 使用 Gibbs 采样器 |
| `adaptive` | XMMC | 自适应 XMMC（每次接受后重算 Fisher 协方差） |
| `window=<size>` | XMMC | 自相关长度窗口大小 (default 16) |

#### Additional parameters for Monte-Carlo analysis

| CLI 参数 | Python | 类型 | 默认值 | 说明 |
|----------|--------|------|--------|------|
| `-s, --seed <seed>` | `seed` | int | 0 | 随机种子。0=固定可复现，>0=指定种子，<0=自动（/dev/urandom+系统时间） |
| `-i, --mcmc-iterations <n>` | `mc_iterations` | int | 1000 | MCMC/EMCE/XMMC/FIMA 的迭代次数 |

#### Clipping outlier data points

| CLI 参数 | Python | 类型 | 默认值 | 说明 |
|----------|--------|------|--------|------|
| `-r, --sigma, --rejection-level <level>` | `rejection_level` | float | 3.0 | 拒绝阈值（sigma 单位） |
| `-n, --iterations <n>` | `rejection_niter` | int | 0 | 最大拒绝迭代次数。0=不做 outlier clipping |
| `--weighted-sigma` | `weighted_sigma=1` | int | 0 | 计算标准差时按权重加权 |
| `--no-weighted-sigma` | `weighted_sigma=0` | int | 0 | 不按权重加权（默认） |

仅 CLLS/NLLM/LMND 方法支持 outlier clipping。

#### Multiple data blocks

| CLI 参数 | Python | 说明 |
|----------|--------|------|
| `-i<key> <file>` | — | 文件 I/O，不适用。Python 层传入单个 numpy 数组 |
| `-c<key> <cols>` | — | 多数据块列定义，不适用。当前仅支持单数据块 |
| `-f<key> <func>` | — | 多数据块模型函数，不适用 |
| `-y<key> <expr>` | — | 多数据块因变量，不适用 |
| `-e<key> <err>` | — | 多数据块误差，不适用 |
| `-w<key> <wgt>` | — | 多数据块权重，不适用 |

#### Constraints

| CLI 参数 | Python | 类型 | 默认值 | 说明 |
|----------|--------|------|--------|------|
| `-t, --constraint, --constraints <expr>` | `constraints` | str/None | None | 线性约束表达式（逗号分隔或多次指定） |
| `-v <name>:=<value>` | 在 `variables` 中使用 `:=` | — | — | 约束语法糖，等价于 `-v name=value -t name=value` |

#### User-defined functions

| CLI 参数 | Python | 类型 | 默认值 | 说明 |
|----------|--------|------|--------|------|
| `-x, --define, --macro <def>` | `macros` | list/None | None | 用户宏定义列表。每个元素格式 `'name(params)=expr'`。宏定义中不可引用外部变量 |

#### Dynamically loaded extensions

| CLI 参数 | Python | 说明 |
|----------|--------|------|
| `-d, --dynamic <lib:array>` | — | 动态加载扩展库，不适用。Python 层使用内建函数 |

#### More on outputs

| CLI 参数 | Python | 类型 | 默认值 | 说明 |
|----------|--------|------|--------|------|
| `--errors, --error-line, --error-columns` | `errdump` | int | 0 | 非零时在结果中包含拟合变量不确定度 |
| `-F, --format <fmt>` | `format` | str/None | None | 变量输出 printf 格式，如 `'a=%16.10g,b=%16.10g'` (default `%12.6g`) |
| `-C, --correlation-format <fmt>` | `correlation_format` | str/None | None | 相关系数矩阵格式 (default `%6.3f`) |
| `-g, --derived-variable[s] <expr>` | `derived_variables` | str/None | None | 派生变量定义，如 `'sum=a+b+c'` |
| `-z, --columns-output <cols>` | `columns_output` | str/None | None | 求值模式中替换输出的列索引 |
| `--delta` | `is_dump_delta` | int | 0 | 在输出中写入残差列 |
| `--delta-comment` | `is_dump_delta` | int | 0 | 同 `--delta`，但残差以注释形式写入 |
| `--residual` | `resdump` | int | 0 | 在输出中写入最终拟合残差值 |
| `-u, --output-fitted <file>` | — | — | — | 文件 I/O，不适用 |
| `-j, --output-rejected <file>` | — | — | — | 文件 I/O，不适用 |
| `-a, --output-all <file>` | — | — | — | 文件 I/O，不适用 |
| `-p, --output-expression <file>` | — | — | — | 文件 I/O，不适用 |
| `-l, --output-variables <file>` | — | — | — | 文件 I/O，不适用 |

#### 内部参数

| Python | 类型 | 默认值 | 说明 |
|--------|------|--------|------|
| `force_nonlinear` | int | 0 | 强制非线性求解。通常自动从 method 推导，无需手动设置 |

### 求值模式 (对应 `--long-help` "Common options for function evaluation")

当 `dependent=None` 时，lfit 进入求值模式（无拟合），仅对输入数据计算表达式值：

```python
r = lfit_fit(data, variables='a,b', columns='x:1,y:2', function='2+3*x')
```

### 完整示例

```python
from pyfitsh.lfit import Lfit, lfit_fit
import numpy as np

# 1. 线性拟合
data = np.loadtxt('lfit_test_linear.dat')
r = lfit_fit(data, 'a,b', 'x:1,y:2', 'a+b*x', 'y')
print(r.params)  # [1.99992, 3.0]

# 2. 二次拟合 + sigma rejection
data = np.loadtxt('lfit_test_outlier.dat')
r = lfit_fit(data, 'a,b,c', 'x:1,y:2', 'a+b*x+c*x*x', 'y',
             rejection_niter=3, rejection_level=3.0)
print(f'nused={r.nused}/{r.nrow}')  # 75/77

# 3. 带误差权重的拟合
data = np.loadtxt('lfit_test_weighted.dat')
r = lfit_fit(data, 'a,b,c', 'x:1,y:2,e:3', 'a+b*x+c*x*x', 'y', error='e')

# 4. 非线性 LM + 自定义宏
r = lfit_fit(data, 'a=1,b=1,c=0.1', 'x:1,y:2', 'a+b*x+c*x*x', 'y',
             method='nllm', parameters='iterations=50')

# 5. MCMC
r = lfit_fit(data, 'a=1:0.2,b=1:0.2,c=0.1:0.1', 'x:1,y:2',
             'a+b*x+c*x*x', 'y',
             method='mcmc', parameters='accepted', seed=12345, mc_iterations=20000)

# 6. EMCE 误差估计 (CLLS 主方法)
r = lfit_fit(data, 'a,b,c', 'x:1,y:2', 'a+b*x+c*x*x', 'y',
             method='emce', parameters='clls', seed=12345, mc_iterations=2000)

# 7. Fisher 信息矩阵
r = lfit_fit(data, 'a=1,b=2,c=0.5', 'x:1,y:2', 'a+b*x+c*x*x', 'y',
             method='fima', derived_variables='sum=a+b+c')

# 8. 星表配对拟合
pairs = np.loadtxt('lfit_test_pairs.dat')
r = lfit_fit(pairs, 'a,b,c', 'xr:1,yr:2,xo:3,yo:4,mr:5,mo:6',
             'a+b*xr+c*yr', 'xo-xr')
```

### CLI ↔ Python 复杂对照示例

**示例 A：LMND + 差分步长 + 自定义宏 + 约束固定变量**

```bash
# CLI
lfit -U -v a=1,b=1,c:=0.5 -c x:1,y:2 \
  -x "poly(a,b,c,t)=a+b*t+c*t*t" \
  -f "poly(a,b,c,x)" -y "y" \
  -q a=0.001,b=0.001 \
  -P iterations=50 \
  data.dat -o result.out
```

```python
# Python
r = lfit_fit(data,
    variables='a=1,b=1,c:=0.5',
    columns='x:1,y:2',
    function='poly(a,b,c,x)',
    dependent='y',
    method='lmnd',
    macros=['poly(a,b,c,t)=a+b*t+c*t*t'],
    differences='a=0.001,b=0.001',
    parameters='iterations=50')
print(r.param_dict)  # {'a': 1.0, 'b': 2.0, 'c': 0.5}
```

**示例 B：XMMC 加权拟合 + 自适应 + 窗口 + 派生变量**

```bash
# CLI
lfit -X -v a=1:0.2,b=1:0.2,c=0.1:0.1 \
  -c x:1,y:2,e:3 -f "a+b*x+c*x*x" -y "y" -e "e" \
  -P adaptive,window=20 \
  -g "sum=a+b+c" \
  -s 12345 -i 10000 \
  weighted.dat -o result.out
```

```python
# Python (类式 API)
lf = Lfit(method='xmmc',
          parameters='adaptive,window=20',
          seed=12345, mc_iterations=10000,
          derived_variables='sum=a+b+c')
r = lf.fit(data,
    variables='a=1:0.2,b=1:0.2,c=0.1:0.1',
    columns='x:1,y:2,e:3',
    function='a+b*x+c*x*x',
    dependent='y', error='e')
print(r['params'], r['nrow'])
```

**示例 C：EMCE 误差估计 + DHSX 主方法 + 扰动 + sigma rejection**

```bash
# CLI
lfit -E -v a=1:0.5,b=1:0.5,c=0.1:0.1 \
  -c x:1,y:2 -f "a+b*x+c*x*x" -y "y" \
  -P dhsx --perturbations 0.01 \
  -r 3 -n 3 \
  -s 12345 -i 2000 --errors \
  outlier.dat -o result.out
```

```python
# Python
r = lfit_fit(data,
    variables='a=1:0.5,b=1:0.5,c=0.1:0.1',
    columns='x:1,y:2',
    function='a+b*x+c*x*x',
    dependent='y',
    method='emce',
    parameters='dhsx',
    perturbations='0.01',
    rejection_niter=3, rejection_level=3.0,
    seed=12345, mc_iterations=2000,
    errdump=1)
print(r.params, r.errors)
```

### 返回值 (LfitResult)

| 属性/key | 类型 | 说明 |
|----------|------|------|
| `params` | ndarray | 拟合参数值 |
| `errors` | ndarray | 参数不确定度 (Fisher/EMCE) |
| `chi2` | float | 拟合残差 chi2 |
| `nrow` | int | 总数据行数 |
| `nused` | int | 实际使用行数 (rejected 后) |
| `nrejected` | int | 被拒绝行数 (= nrow - nused) |
| `residual_sigma` | float | 残差标准差 |
| `acceptance` | float | MCMC/XMMC 接受率 |
| `ok` | bool | 拟合是否成功 |
| `param_dict` | dict | {变量名: 值} |
| `variables` | list | 变量名列表 |
| `used_mask` | list | 每行是否被使用 |
| `chain` | ndarray/None | 采样链 (chain_count, nvar+1)，最后一列为 chi2。MCMC/XMMC/EMCE/FIMA/MCHI 自动返回 |
| `cov_matrix` | ndarray/None | 协方差矩阵 (nvar, nvar)。FIMA/XMMC 自动返回 |
| `corr_matrix` | ndarray/None | 相关矩阵 (nvar, nvar)。从 cov_matrix 计算 |
| `eval_data` | ndarray/None | eval 模式输出 (nrow, ncol)。dependent=None 时自动返回 |
| `error_code` | int | 错误码 (0=成功) |
| `error_msg` | str | 错误信息 |

支持属性访问 (`r.params`) 和 dict 访问 (`r['params']`)，以及 `r.keys()`、`r.items()` 遍历。

### 测试覆盖

对标 origincode 3 个测试脚本：

| 脚本 | 测试数 | 对比方式 | 结果 |
|------|--------|---------|------|
| test_lfit.sh | 19 | params 逐值对比 | 19/19 PASS |
| test_lfit_montecarlo.sh | 19 | chain 逐行对比 + FIMA params/errors/corr | 19/19 PASS |
| longtest_lfit.sh | 9 | chain 逐行对比 + FIMA | 7/9 PASS, 2 DIFF (CLI 精度截断) |

### 性能 (-O3 编译)

pairs MCMC 20000 iterations 对比：

| 版本 | 时间 | 相对 C 原版 |
|------|------|------------|
| C 原版 (-O2) | 221s / 236s | 1.00x |
| Cython -O3 | 211s / 221s | 0.95x (快 5%) |
