# Cython 版 ficonv 代码复用分类报告

---

## 一、直接复用（源码完全一致，diff 无输出）

### .c 源文件（17 个）

| 文件 | 提供的功能 |
|------|-----------|
| `kernel-io.c` | kernel 解析、`create_kernels_from_kernelarg`、`kernel_init_images`（定义了但 Cython 路径未调用） |
| `tensor.c` | 张量分配 |
| `fbase.c` | 基础数据类型 |
| `statistics.c` | 统计函数（median, truncated_mean 等） |
| `projection.c` | WCS 投影 |
| `math/tpoint.c` | 点集合操作 |
| `math/cpmatch.c` | 坐标匹配 |
| `math/delaunay.c` | Delaunay 三角剖分 |
| `math/trimatch.c` | 三角形匹配 |
| `math/poly.c` | `eval_2d_poly`、`eval_2d_monoms`、`calc_2d_unitarity` |
| `math/polyfit.c` | `fit_2d_poly` |
| `math/spmatrix.c` | 稀疏矩阵 |
| `math/convexhull.c` | 凸包 |
| `math/fit/lmfit.c` | `solve_gauss`、矩阵/向量分配 |
| `math/spline/bicubic.c` | 双三次样条 |
| `math/spline/biquad.c` | 双二次样条 |
| `math/spline/biquad-isc.c` | 双二次子像素积分 |
| `math/spline/spline.c` | 自然样条 |

### .h 头文件（2 个）

| 文件 | 说明 |
|------|------|
| `longhelp.h` | 长帮助宏 |
| `tensor.h` | 张量声明 |

---

## 二、直接复用的源码，但有局部修改

| 文件 | 修改内容 | 影响 |
|------|---------|------|
| `kernel-base.c` | + profiling 计数器（24行）、+ `fitrans_profile_dump()`、+ 3 处计时包裹、+ `FITRANS_OMP` pragma（未启用） | 无算法改动，轻微运行时开销 |

---

## 三、修改过的头文件

| 文件 | 修改内容 |
|------|---------|
| `fitsh.h` | `#include "../config.h"` → `#include "config.h"`；`logmsg()` 函数声明改为 `static inline` 空实现 |
| `kernel.h` | 顶部加了 `#include <fits/fits.h>` |
| `config.h` | 原版无此文件（原版由 configure 生成 `config.h.in`），Cython 版自定义了 FITSH_VERSION 等宏 |

---

## 四、重新实现（Cython 自建 .c/.h，非原版文件的拷贝）

### 核心文件

| Cython 文件 | 对应原版实现 | 策略 |
|-------------|-------------|------|
| `_core.c` + `_core.h` | `grmatch_lib.c` | **部分拷贝**：`do_pointmatch()` 和其 static 辅助函数逐行复制，外加 flat wrapper（`do_pointmatch_flat`、`coord_match_flat`、`id_match_flat`）将 iline 结构体展开为平面数组。**跳过**了 `read_match_data_points`、`normalize_columns`、`get_number_list`、`do_coordmatch`、`do_idmatch` 等 IO/UI 函数 |
| `_trans_core.c` + `_trans_core.h` | `transform.c` | **重新实现**：`trans_apply_flat()` flat wrapper，内部 Jacobi + 正向/逆向求值 |
| `_wcs_core.c` + `_wcs_core.h` | `wcsfit.c`（原版） | **重新实现**：`wcs_fit_flat()` + 辅助函数 |
| `_polyfit_core.c` + `_polyfit_core.h` | `math/polyfit.c`（部分） | **重新实现**：`poly_fit_flat()`、`poly_compose_affine_flat()` 平面数组版 |

### fitrans 相关

| Cython 文件 | 对应原版 | 策略 |
|-------------|---------|------|
| `_fitrans_ops.c` + `_fitrans_ops.h` | `fitrans.c` 中 noise/zoom/shrink/smooth 函数 | **提取重写**：4 个函数改为 `_flat` 平面数组版，去掉 FITS 头依赖 |
| `_fitrans_core.c` + `_fitrans_core.h` | `fitrans.c` 中图像变换/插值部分 | **重新实现**：6 种插值方法的平面数组版 |

### ficonv 相关（核心差异区）

| Cython 文件 | 对应原版 | 策略 |
|-------------|---------|------|
| `_ficonv_core.c` + `_ficonv_core.h` | `ficonv.c` 的全部 main 流程 | **大幅缩短**：直接调 `fit_kernel_poly_coefficients_block` + `convolve_with_kernel_set`，跳过了中间层 |
| `_ficonv_kernel.c` | `kernel-io.c` 的 `create_kernels_from_kernelarg` | **重新实现但未编译**（不在 setup.py sources 中），属于死代码 |

### Python 层

| 文件 | 策略 |
|------|------|
| `fitsh_cy.pyx` | **全新**：Cython 包装类（Grmatch、Grtrans、WCSFitter、PolyFitter、Fitrans、_ConvolutionOp 等） |
| `fitsh_cy.pxd` | **全新**：C 函数/类型声明 |
| `setup.py` | **全新**：构建配置 |
| `__init__.py` | **全新**：包导出 |

---

## 五、被完全跳过的原版函数

### 来自 `ficonv.c`

| 函数 | 作用 | 跳过后果 |
|------|------|---------|
| `fit_kernels()` | 迭代拟合包装（mask 合并、扩张、niter 循环剔除） | 无迭代、无 mask 预处理 |
| `fit_kernels_native()` | 前景 mask + 加权拟合 | 无背景迭代模式 |
| `create_level_mask()` | 多级 mask | 缺失 |
| `create_foreground_mask()` | 前景 mask | 缺失 |
| `create_background_mask()` | 背景 mask | 缺失 |
| `create_image_weight()` | 图像加权 | 缺失 |
| `make_subtracted_image()` | 减法残差图像 | 缺失 |
| `create_weights()` | 权重数组构建 | 缺失 |
| `mark_stamps()` | stamp 标记 | 缺失 |
| `stamp_parse_argument()` | stamp 命令行解析 | 缺失 |
| `stamp_read_file()` | stamp 文件读取 | 缺失 |
| `ficonv_dump_*()` | dump 中间图像 | 缺失 |
| `fprint_ficonv_usage/long_help()` | 帮助文本 | 不需要 |
| `kernel_init_images()` | **kernel image 初始化** | **P0 缺陷：`k->image` 始终 NULL** |

### 来自 `fitrans.c`

| 函数 | 作用 | 跳过后果 |
|------|------|---------|
| `affine_img_transformation_spline_fit()` | 预样条拟合 | 缺失 |
| `affine_img_transformation_spline_eval()` | 预样条求值 | 缺失 |

### 来自 `grmatch_lib.c`

| 函数 | 作用 | 跳过原因 |
|------|------|---------|
| `read_match_data_points()` | 文件读取 | Python 侧替代 |
| `normalize_columns()` | 列归一化 | Python 侧替代 |
| `get_number_list()` | 数字列表解析 | Python 侧替代 |
| `colinfo_reset()` | 列信息重置 | 不需要 |
| `do_coordmatch()` | 坐标匹配主函数 | 由 `coord_match_flat` 替代 |
| `do_idmatch()` | ID 匹配主函数 | 由 `id_match_flat` 替代 |

### 整个被跳过的源文件

| 原版文件 | 作用 | 跳过原因 |
|---------|------|---------|
| `ui.c` | UI 工具函数 | 不需要 |
| `io/iof.c` | 文件 I/O | 不需要 |
| `io/scanarg.c` | 命令行解析 | 不需要 |
| `io/tokenize.c` | 字符串分割 | 不需要（内核已含） |
| `io/format.c` | 格式化 | 不需要 |
| `common.c` / `common.h` | 公共工具 | 不需要 |
| `history.c` / `history.h` | FITS 历史 | 不需要 |
| `fitsmask.c` / `fitsmask.h` | FITS mask 操作 | 不需要 |
| `longhelp.c` | 长帮助格式化 | 不需要 |

---

## 六、流程被缩短或替换的关键路径

| 原版流程 | Cython 流程 | 缩短/替换程度 |
|---------|------------|-------------|
| `main()` → `scanarg` → FITS I/O → mask 合并 → `kernel_init_images`×2 → `fit_kernels()` → `fit_kernels_native()` → mask 扩张 → 前景mask → `fit_kernel_poly_coefficients_block(bdc=32)` × niter 轮 → `convolve_with_kernel_set` → mask/dump/output | `ficonv_fit_flat()` → `create_kernels_from_kernelarg` → 直接 `fit_kernel_poly_coefficients_block(bdc=1)` 一次 → `convolve_with_kernel_set` | 中间层全部省略，总共跳过了 ~400 行流程代码 |
| `convolve_point()` 访问 `k->image` | 同函数 | `k->image` 未初始化 → 访问 NULL |
| xlist = NULL (空列表) | xl.kernels = kl.kernels (自别名) | 语义完全不同的错误构造 |
| bdc=32 分块 | bdc=1 整图单块 | 矩阵组装策略改变 |

---

## 小结

- **直接复用**：17 个 .c（全部 math/ + tensor + fbase + statistics + projection + kernel-io）+ 2 个 .h（longhelp.h, tensor.h）
- **复用但有修改**：1 个 .c（kernel-base.c，profiling + OMP）+ 3 个 .h（fitsh.h, kernel.h, config.h）
- **重新实现**：7 个 .c + 7 个 .h（_core, _trans_core, _wcs_core, _polyfit_core, _fitrans_core, _fitrans_ops, _ficonv_core）
- **被跳过**：~20 个函数 + 8 个完整源文件
- **死代码**：`_ficonv_kernel.c`（未编译进 .so）
