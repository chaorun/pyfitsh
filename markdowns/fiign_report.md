# fiign 移植报告

**日期**: 2026-06-18

## 移植内容

从 origincode `src/fiign.c` (669行 CLI 程序) 移植到 `pyfitsh/fiign/`。

### 保留的核心算法
- `cosmics_ignore` — 5×5 局部统计宇宙线检测 + 可选替换
- `saturated_mark` — 饱和像素标记 + blooming 方向
- `integerlimit_mark` — 整数范围限制（bitpix 8/16/32）
- `mask_block_draw` — 几何形状 mask 绘制（block/circle/line/pixel）
- `mask_convert_inmem` — mask 条件转换
- `parse_mask_flags_simple` — mask 名称解析（无 scanflag 依赖）
- `fiign_apply` — 统一入口，执行完整管线

### 删除的 I/O（Python 层替代）
- `fits_read/write`, `fopen/fclose`, `fits_rescale/backscale`
- `fits_mask_read/write/export_as_header`, `fits_basename`
- `join_masks_from_files` — Python 层用 numpy `|=` 预合并

### 改为返回值
- `fits_history_export_command_line` → `char *history` 字段，返回参数字符串

### 新增文件
- `fiign/fiign_core.h` — MASK 常量 + 结构体 + API 声明
- `fiign/fiign_core.c` — 核心算法（~400行）
- `test/test_fiign/` — 测试脚本

### Cython 接口
- `core.pyx`: `class Fiign` — `apply(img, mask, saturation, ...)` → `{'data', 'mask', 'history'}`
- `__init__.py`: 导出 `Fiign`

## 测试结果

13 个测试全部 PASS（2560×2560 图像，6553600 像素）：

| # | 测试 | max_diff | tol | 状态 |
|---|------|----------|-----|------|
| 01 | ignore negative + apply mask | 0 | 0 | OK |
| 02 | ignore zero + apply mask | 0 | 0 | OK |
| 03 | ignore nonpositive + apply mask | 0 | 0 | OK |
| 04 | saturation scalar + apply mask | 0 | 0 | OK |
| 05 | saturation image + apply mask | 0 | 0 | OK |
| 06 | block rectangle + apply mask | 0 | 0 | OK |
| 07 | block circle + apply mask | 0 | 0 | OK |
| 08 | block pixel + apply mask | 0 | 0 | OK |
| 09 | block line + apply mask | 0 | 0 | OK |
| 10 | expand mask + apply mask | 0 | 0 | OK |
| 12 | input external mask + apply mask | 0 | 0 | OK |
| 13 | convert hot→cosmic + apply mask | 0 | 0 | OK |
| 14 | ignore input mask + apply mask | 0 | 0 | OK |

**13/13 逐像素零差异**，与 CLI 输出完全一致。
