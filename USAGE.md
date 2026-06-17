# pyfitsh 用法

## Fistar — 星检测与拟合

```python
from pyfitsh.fitsh_cy import Fistar
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

参数匹配 CLI `--long-help`：

| Python 参数 | CLI 参数 | 说明 |
|-------------|----------|------|
| `threshold` | `-t` | 检测阈值 (default 100.0) |
| `flux_threshold` | `-f` | 通量阈值，>0 时替代 threshold 用于 UPLINK 算法 |
| `skysigma` | `-d` | 天空背景噪声 (sky sigma) |
| `gain` | `-g` | 探测器增益 (e-/ADU) |
| `algorithm='uplink'\|'parabolapeak'` | `--algorithm` | 候选检测算法 |
| `model='elliptic'\|'gauss'\|'deviated'` | `--model` | 星点模型 |
| `only_candidates=True` | `--only-candidates` | 仅候选检测，跳过拟合 |
| `psf='native,order=2,grid=4'` | `--psf` | PSF 拟合参数 |
| `sort='x'\|'y'\|'peak'\|'flux'\|'fwhm'` | `-s` | 排序方式 |
| `section=(x1,x2,y1,y2)` | `--section` | 图像区域限制 |
| `output_mark=True` | `--mark-output` | 输出标记图像 |
| `output_area=True` | `--output-area` | 输出区域图像 |

高斯模型 + PSF 拟合示例（对应 CLI 的 `-d 5 -f 10000 --model gauss --psf native,order=2`）：

```python
fs = Fistar(threshold=100.0, flux_threshold=10000.0, skysigma=5.0,
            algorithm='uplink', model='gauss', model_order=2,
            it_sym=4, it_gen=2, gain=2.0,
            sort='x', psf='native,order=2', psf_output=True,
            verbose=True)
r = fs.do_fistar(img_data, mask=mask_data)
# r.output (Table), r.nstar (int), r.psf (astropy HDU)
```

**已知问题**：parabolapeak 算法拟合时 bg/amp/flux 与 CLI 有 ~1% 差异；`input_candidates` 二次引用候选取值异常。

## Fiphot — 孔径测光

### 基本用法

```python
from pyfitsh.fitsh_cy import Fiphot
import numpy as np

stars = np.loadtxt('positions.dat')
fp = Fiphot(apertures='10:18:12', gain='2', mag_flux=(10.0, 10000.0))
r = fp.photometry(img, stars=stars, col_xy=(1, 2), col_id=0)
# r.output (Table), r.dict (backward compat dict)
```

多孔径：`apertures='8:18:12,10:18:12,12:20:15'`（逗号分隔）。

### photometry_from_raw

对应 CLI `--input-raw-photometry <file> -s <subtracted_image>`，从 pre-computed raw 测光数据在减影图像上执行 subtraced photometry：

```python
# 1. 标准测光获取 raw 数据
r = fp.photometry(img, stars=stars, col_xy=(1,2), col_id=0)

# 2. 从 raw 做 subtracted photometry
r2 = fp.photometry_from_raw(r.output, subtracted_img)
# r2.output (Table), r2.dict

# 可选：传入参考星等
r2 = fp.photometry_from_raw(r['raw'], subtracted_img,
                            ref_mag=ref_mag_arr, ref_col=ref_col_arr)

# 可选：kernel 减影模式
r2 = fp.photometry_from_raw(r['raw'], subtracted_img,
                            kernel_spec='gauss:3:3', normalize_kernel=True)
```

已验证与 CLI `--input-raw-photometry` 一致，49 颗星对比 flux 差异 < 0.6%（来源：CLI ASCII raw 文件精度截断，非实现错误）。

### subtracted_photometry

对应 CLI 减影测光模式：`subtracted_photometry()`。

### magfit

对应 CLI `--magfit`：`magfit()`。

## Grmatch — 星点匹配

```python
from pyfitsh.fitsh_cy import Grmatch

ref = np.loadtxt('ref.cat')
inp = np.loadtxt('obs.cat')
g = Grmatch(order=2, maxdist=2.0, unitarity=0.2, ttype=2,
            use_ordering=1, is_centering=1, nmiter=4, rejlevel=3.0)
r = g.matchpoints(ref[:,1], ref[:,2], inp[:,1], inp[:,2],
                  ref_order=-ref[:,3], inp_order=-inp[:,3])
print(f'nhit={r.nhit}, trans={r.transformation}')
```

参数匹配 CLI 的 `grmatch --long-help`。

## Fitrans — 图像变换

```python
from pyfitsh.fitsh_cy import Fitrans

trans = {'type':'polynomial', 'order':1, 'offset':(0.0,0.0), 'scale':1.0,
         'dxfit':[-148.653, 0.99898, 0.04297],
         'dyfit':[140.218, -0.04291, 0.99912]}
ft = Fitrans.ImageTransformation(transformation=trans, inverse=True, method='bilinear')
output, mask = ft.apply(img)
```

`method='bilinear'` 匹配 CLI 默认 `steplike+interpolate`。

## Ficonv — 卷积/减影

完整示例见 `examples/ficonv_demo.py`（含三种核型：多高斯、离散核、空间变化核）。

**基本用法**：

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

**用 kernel_dict 跳过拟合（等价 CLI `--input-kernel-list`）**：

```python
# kernel_dict 从上一次拟合获得（或保存/加载）
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
        # ...
    ]
}
```
类型码：1=background, 2=identity, 3=gaussian。

**与 CLI 对比时必须注意 mask 处理**：
- CLI 的 `fits_mask_read_from_header()` 自动从 FITS header 读取 MASKINFO 并生成 mask。Cython 版需显式调用 `decode_maskinfo(hdu)` 解析
- 若 FITS 文件包含 MASKINFO（如边缘坏区标记），不传入 mask 将导致拟合在坏像素上进行，结果与 CLI 完全不同
- ref 图像通常无 MASKINFO，传入全零 mask 即可
- 若图像存在 NaN 像素，需在 Python 层标记为 mask（`mask[np.isnan(data)] = 1`），匹配 CLI 的 `fits_mask_mark_nans()` 行为

参数匹配 CLI `--long-help`：

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

## Firandom — 人工图像生成

```python
from pyfitsh.fitsh_cy import Firandom
import numpy as np

fir = Firandom(sx=512, sy=512, gain=2.0, seed=12345,
               background='1000', is_photnoise=True, zoom=1)
result = fir.generate(
    stars='20[x=r(0.1,0.9),y=r(0.1,0.9),f=6,e=0.05,p=r(0,180),i=r(50000,200000)]',
    background_stddev=5.0)
img = result.image  # numpy 2D float64 array, 逐像素匹配 CLI
```

**种子参数**：
- `seed` — 通用种子，当 `seed_noise`/`seed_spatial`/`seed_photon` 未指定时的默认值
- `seed_noise` — 背景噪声种子（对应 CLI `--seed-noise`）
- `seed_spatial` — 空间坐标种子（对应 CLI `--seed-spatial`）
- `seed_photon` — 光子噪声种子（对应 CLI `--seed-photon`）

**Star 输入格式**（两种）：

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

**绘制星点到已有图像**：
```python
img = np.zeros((512, 512), dtype=np.float64)
fir.draw_stars(img, stars)
```

参数对照 CLI：

| Python 参数 | CLI 参数 | 说明 |
|-------------|----------|------|
| `sx, sy` | `--size` | 图像尺寸 |
| `gain` | `--gain` | 增益 (e-/ADU) |
| `background` | `--sky` | 背景表达式 |
| `background_stddev` | `--sky-noise` | 背景额外噪声 |
| `is_photnoise` | `--photon-noise` | 光子噪声仿真 |
| `dontquantize` | (默认) | 默认 True，匹配 CLI 不量化行为 |
| `method=1` | `--integral` | 积分绘制模式 |
| `seed` | `--seed` | 通用随机种子 |

## 通用注意事项

- Fistar / Fiphot 返回 `types.SimpleNamespace`，使用 `r.output` 访问主结果 Table，`r.nstar` 访问星数，`r.dict` 向后兼容字典
- 图像数据均为 float64 numpy 数组，C-order (row-major)
- mask 为 uint8 numpy 数组，0 = 有效像素
- `decode_maskinfo(hdu)` 读取 FITS 头中 MASKINFO 关键字生成 mask
