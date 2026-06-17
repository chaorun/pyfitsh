# fiphot 输出规范

所有 `fiphot.Fiphot` 公开测光方法均返回 `types.SimpleNamespace` 对象，

提供 `.output`（`astropy.table.Table`）和 `.dict`（向后兼容 `dict`）两种访问方式。

---

## 1. `photometry()`

对应 CLI:
```
fiphot --input <image> --input-list <positions> --apertures <apertures> ...
       --format IXY,BbFfMm --output <cat> --output-raw-photometry <raw>
```

| 属性 | 类型 | 说明 |
|------|------|------|
| `.output` | `astropy.table.Table` | 主测光结果，每行一个孔径 |
| `.output_raw_photometry` | `astropy.table.Table` | raw photometry 数据，CLI 固定格式 |
| `.dict['output_raw_photometry_3d_data']` | `numpy.ndarray` (nstar, nap, 12) | 内部 raw 数组 |

### `.output` 列

| 列名 | 类型 | CLI tag | 说明 |
|------|------|---------|------|
| `id` | str/int | I | 星标识 |
| `x` | float64 | X | 质心 X |
| `y` | float64 | Y | 质心 Y |
| `flux` | float64 | F | 流量 |
| `flux_err` | float64 | f | 流量误差 |
| `bg` | float64 | B | 背景水平 |
| `bg_err` | float64 | b | 背景散射 |
| `mag` | float64 | M | 星等 |
| `mag_err` | float64 | m | 星等误差 |

### `.output_raw_photometry` 列

| 列名 | 类型 | 说明 |
|------|------|------|
| `id` | str/int | 星标识 |
| `x` | float64 | 质心 X |
| `y` | float64 | 质心 Y |
| `nap` | int | 孔径数量 |
| `ref_mag` | float64 | 参考星等（无则为 `np.nan`） |
| `ref_col` | float64 | 参考色指数（无则为 `np.nan`） |
| `r0` | float64 | 孔径半径 |
| `ra` | float64 | 背景内半径 |
| `da` | float64 | 背景宽度 |
| `flux` | float64 | 流量 |
| `flux_err` | float64 | 流量误差 |
| `flag` | int | 标志位（hex 格式） |

### `.dict` 键

`flux`, `fluxerr`, `magnitude`, `magerr`, `bgarea`, `bgflux`, `bgmedian`, `bgsigma`,
`cntr_x`, `cntr_y`, `cntr_width`, `cntr_w_d`, `cntr_w_k`,
`cntr_x_err`, `cntr_y_err`, `cntr_w_err`,
`flag`, `rtot`, `rbad`, `rign`, `atot`, `abad`,
`id`（可选）, `optimal_r0`/`optimal_ra`/`optimal_da`（当 `calc_optimal=True`）

---

## 2. `photometry_from_raw()`

对应 CLI:
```
fiphot -s <sub_image> --input-raw-photometry <raw> --input-kernel <kernel> ...
```

| 属性 | 类型 | 说明 |
|------|------|------|
| `.output` | `astropy.table.Table` | 减影测光结果 |
| `.dict` | `dict` | 向后兼容 |

Table 列与 `photometry().output` 相同。

### 输入参数 `raw_data`

支持两种形式：
- `numpy.ndarray` (nstar, nap, 12) — 3D raw 数组
- `astropy.table.Table` — `photometry().output_raw_photometry` 的返回值

### `.dict` 键

与 `photometry().dict` 相同（无 `output_raw_photometry_3d_data`）。

---

## 3. `subtracted_photometry()`

对应 CLI 减影测光模式。

| 属性 | 类型 | 说明 |
|------|------|------|
| `.output` | `astropy.table.Table` | 减影测光结果 |
| `.dict` | `dict` | 向后兼容 |

Table 列与 `photometry().output` 相同。

---

## 4. `magfit()`

对应 CLI:
```
fiphot --magfit orders=<c0>:<c1>:...,iterations=<d>,sigma=<s>
```

| 属性 | 类型 | 说明 |
|------|------|------|
| `.output` | `astropy.table.Table` | 星等拟合结果 |
| `.dict` | `dict` | 向后兼容 |

### `.output` 列

| 列名 | 类型 | 说明 |
|------|------|------|
| `mag` | float64 | 校准星等 |

### `.dict` 键

`mag` (numpy array), `ninit` (int), `nrejs` (int), `nstar` (int), `naperture` (int)

---

## 使用示例

```python
from fitsh_cython.fitsh_cy import Fiphot

fp = Fiphot(apertures='10:18:12', gain='2', mag_flux=(10.0, 10000.0))

# photometry
r = fp.photometry(img, stars=stars, col_xy=(1,2), col_id=0)
print(r.output)                    # Table
print(r.output['id'], r.output['flux'])

# photometry_from_raw（支持 Table 输入）
r2 = fp.photometry_from_raw(r.output_raw_photometry, sub_img, kernel_dict=kd)
# 或
r2 = fp.photometry_from_raw(r.dict['output_raw_photometry_3d_data'], sub_img, kernel_dict=kd)

# 向后兼容 dict 访问
print(r.dict['flux'], r.dict['bgmedian'])
```
