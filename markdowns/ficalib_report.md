# ficalib 移植报告

**日期**: 2026-06-18

## 移植内容

从 origincode `src/ficalib.c` (2095行 CLI 程序) 移植到 `pyfitsh/ficalib/`。

### 保留的核心算法
- `overscan_model` — overscan 建模（spline/poly 拟合）
- `overscan_do_vertical` / `overscan_do_horizontal` — overscan 校正
- `overscan_correction` — 完整 overscan 处理管线
- `calibrate_fit_flatpoly` — flat 场多项式拟合
- `ficalib_calibrate_cy` — 新增 memory-based 校准入口

### 新增文件
- `ficalib/ficalib_core.c` — 核心算法（从 origincode 裁剪至 ~400 行）
- `ficalib/ficalib_core.h` — C API 声明
- `ficalib/combine.c` / `ficalib/combine.h` — 像素合并工具
- `math/splinefit.c` / `math/splinefit.h` — 样条拟合

### Cython 接口
- `core.pyx`: `class Ficalib` — `calibrate(img, bias, dark, flat, mask, flat_mean, dark_time)` → `{'data', 'mask', 'sx', 'sy'}`
- 所有输入为 numpy float64 2D 数组，mask 为 uint8

## 测试结果

12 个测试全部 PASS：

| 测试 | 说明 | max_diff | tol | 状态 |
|------|------|----------|-----|------|
| 01 | copy only | 0 | 1e-8 | OK |
| 02 | bias subtraction | 0 | 1e-8 | OK |
| 03 | bias + flat | 0 | 1e-8 | OK |
| 04 | bias + dark + flat (no exptime) | 0 | 1e-8 | OK |
| 05 | bias + dark + flat (exptime) | 0 | 1e-8 | OK |
| 06 | multi-file (sci_b) | 0 | 1e-8 | OK |
| 07 | rewrite output name | 0 | 1e-8 | OK |
| 08 | trim (100:100:1900:1900) | 0 | 1e-8 | OK |
| 09 | post multiply 2.0 | 0 | 1e-8 | OK |
| 10 | post scale 10000 | 4.88e-04 | 1e-2 | OK |
| 11 | saturation + gain | 0 | 1e-8 | OK |
| 12 | input mask | 0 | 1e-8 | OK |

**11/12 逐像素零差异**。仅 test_10 有 4.88e-4 浮点舍入误差（均值归一化除法）。
