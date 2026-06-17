# fiphot CLI vs Cython 全流程对比测试

## 测试命令

### CLI 三步流程

```bash
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

# 第三步：核拟合
ficonv -r L1_GRB260604C-Ic-180s-20260606_210318.fits -i v1.fits \
  -k "i/1;b/1;g=5,2.0,2/1" --output-kernel-list kernel_spatial1.txt \
  -o ref_conv_spatial1.fits --output-subtracted sub_spatial1.fits -n 2 -s 3
```

### Cython 等价代码

见 `test_fiphot_full.py`。

## 关键发现：bgmode 不匹配

**原因**：CLI 两步使用不同的 bgmode（第一步有 `--sky-fit`，第二步无 `--sky-fit`，使用默认值）。

Cython 中 `photometry()` 和 `photometry_from_raw()` 共享同一个 `Fiphot` 实例的 `sky_fit` 参数。必须用两个实例分别设置：
- 第一步：`Fiphot(sky_fit='median,iterations=2,sigma=3')`
- 第二步：`Fiphot()`（默认 bgmode）

**之前的错误**：一直用同一个 `sky_fit` 参数导致 `photometry_from_raw` 中 `bgm.rejniter=2` 而 CLI 第二步 `bgm.rejniter=0`，背景估计差异放大为 ~241 flux 偏差（0.36%）。

## 最终结果

| 指标 | 值 |
|------|-----|
| 对比星数 | 45 |
| 最大差异 | 4.00 |
| 平均差异 | 0.35 |
| 中位差异 | 0.02 |

差异来源：CLI ficonv（C 编译器）与 Cython Ficonv（Python/Cython）的独立拟合过程产生的浮点精度差。
