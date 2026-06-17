# 原版 ficonv vs Cython (`fitsh_cython/`) 逐项差异报告

---

## 一、整体流程对照

| # | 原版 `ficonv.c:main()` 步骤 | Cython `_ficonv_core.c:ficonv_fit_flat()` | 差异与影响 |
|---|------|------|------|
| 1 | 命令行解析 (`scanarg`) | 无（Python 侧传参） | 设计差异，非 bug |
| 2 | FITS 文件 I/O + rescale | 无（Python 侧 numpy） | 设计差异 |
| 3 | mask 读取+合并 (`fits_mask_and`) | **无 mask 合并** | **结果可能不同** |
| 4 | stamp 解析/文件读取 | **无** | 缺少 stamp 支持 |
| 5 | `kernel_init_images(klist)` | **`kernel_init_images` 未调用** | **致命：高斯 kernel image 未初始化，`convolve_point()` 访问 NULL 指针** |
| 6 | `kernel_init_images(xlist)` | 未调用 (xlist 为别名) | 同上 |
| 7 | 预样条变换 (`affine_img_transformation_spline_fit/eval`) | **无** | 缺少 `-p` 功能 |
| 8 | `fit_kernels()` → `fit_kernels_native()` | 直接调 `fit_kernel_poly_coefficients_block` | 跳过了 mask 处理、迭代拟合 |
| 9 | 迭代拟合循环 (niter, 像素剔除) | **niter 参数传入但未使用** | **结果可能不同** |
| 10 | 前景 mask 构建 (`create_foreground_mask`) | **无** | 影响 background-iterative 模式 |
| 11 | 加权拟合 (`create_image_weight`) | **无** | 缺少 `-w` 功能 |
| 12 | mask 扩张 (`fits_mask_expand_false`) | **无** | kernel 边缘像素处理不同 |
| 13 | `convolve_with_kernel_set` → 卷积 | 同样调用 | 函数调用一致，但内核状态不同 |

---

## 二、源码文件级别逐项差异

### 2.1 `ficonv.c` (原版) vs `_ficonv_core.c` (Cython)

| 差异项 | 原版 | Cython | 后果 |
|--------|------|--------|------|
| **块大小 bdc** | `bdc=32` (32×32 子块) | `1, 1` (全图一整个块) | **性能：4M 像素的单次矩阵组装 vs 1024 个 4K 像素的小块。计算总量相同但缓存行为不同** |
| **xlist 构造** | `xlist->nkernel=0, xlist->kernels=NULL` (空) | `xl.kernels=kl.kernels, xl.nkernel=kl.nkernel` (别名 klist) | **危险性：subblock 内 realloc 自拷贝操作** |
| **`kernel_init_images` 调用** | main 中调用 2 次 | **未调用** | **`k->image` 始终为 NULL，`convolve_point` 行 210 访问 `k->image[hsize+i][hsize+j]` 应 segfault** |
| **迭代拟合循环** | 完整实现（niter 轮 + 每轮像素剔除） | 未实现（声明了参数但函数体内没用） | 结果差异 |
| **背景前景 mask** | `create_foreground_mask` + `create_background_mask` + `create_weights` | 全无 | 功能缺失 |
| **mask 合并** | `fits_mask_and(mask, ..., mask_ref); fits_mask_and(mask, ..., mask_img); fits_mask_and(mask, ..., inmask)` | 无（直接用 Python 传入的 mask） | 输入 mask 语义不同 |
| **mask 扩张** | `fits_mask_expand_false(mask, sx, sy, hsize, ...)` | 无 | kernel 边缘像素未被 mask |
| **dump 功能** | 有 (`is_dump`, `ficonv_dump_image`) | 无 | - |
| **output kernel list** | 有 (`kernel_info_write`) | 无 | - |
| **unity kernels** | 有 (`-u` 归一化) | 无 | - |
| **output subtracted** | `make_subtracted_image` 含 xlist 卷积修正 | 无 | - |
| **添加图像 (-a)** | 有 (`addimg`) | 无 | - |

### 2.2 `kernel-base.c` (原版) vs `fitsh_cython/kernel-base.c`

| 差异项 | 原版 (origincode) | fitsh_cython | 影响 |
|--------|------------------|-------------|------|
| `#include <time.h>` | 无 | 有 | 无影响 |
| profiling 计数器+`fitrans_profile_dump()` | 无 | 有（24 行） | **轻微运行时开销** |
| `kernel_image_calc_gaussian` 计时 | 无 | `clock_t _pg0` + 2 行累加 | **轻微开销** |
| `fit_kernel_poly_coefficients_subblock` 入口 | 变量声明后直接 `sx=...` | 添加 `clock_t _pt0` + profile block (5 行) | **轻微开销 + 打乱了 `monoms`, `kernels`, `ox,oy,scale` 等变量的声明位置** |
| `solve_gauss` 调用处 | 直接调 | 包裹 `{ clock_t _ps0; solve_gauss(...); 累加 }` | **轻微开销** |
| subblock 末尾 | 直接 `return(0)` | `_prof_subblock_time += ...` 再 `return(0)` | 轻微开销 |
| **FITRANS_OMP pragma** | 无 | `#ifdef FITRANS_OMP` / `#pragma omp parallel for collapse(2) private(j)`（卷积函数内） | **若无 `-DFITRANS_OMP` 编译则无效；若定义了，并行会引入线程开销和结果不确定性** |

### 2.3 `kernel-io.c` — 完全相同

两个目录的 `kernel-io.c` diff 无输出，完全一致。

### 2.4 `_ficonv_kernel.c` — 完全未编译

`setup.py` 的 sources 列表中**没有** `_ficonv_kernel.c`，它是一个独立重新实现的 kernel 解析器，但**从未被编译进 `.so`**。实际使用的是 `kernel-io.c` 的原版 `create_kernels_from_kernelarg`。

---

## 三、编译参数对比

| 参数 | C 二进制 (Makefile) | Cython (.so) |
|------|-------------------|-------------|
| 优化级别 | `-O2 -g` | `-O3 -DNDEBUG -ffast-math` |
| 标准 | `c99` | (Cython 自动) |
| 告警 | `-Wall -Wextra -Wno-unused-parameter` | 无 |
| OpenMP | 无 | 无（FITRANS_OMP pragma 未生效） |
| include 路径 | 本地 | 额外加了 conda `/include` 和 openblas `/include`（`fits/fits.h` 可能被覆盖） |

---

## 四、可能导致结果或速度变化的关键差异排序

| 优先级 | 差异 | 位置 | 可能后果 |
|--------|------|------|------|
| **P0** | `kernel_init_images()` 未调用 | `_ficonv_core.c:28` vs `ficonv.c:main` | **k->image 为 NULL。高斯核应在 `convolve_point` 行 210 崩溃。若未崩溃说明内存地址 0 被映射了或访问被编译器优化掉，但结果必然错误。** |
| **P0** | xlist 别名 klist | `_ficonv_core.c:29` | **subblock 内 realloc 自拷贝，重新分配数组时 xlist->kernels 可能变成悬空指针** |
| **P1** | 块大小 `bi=1, bj=1` vs `bi=bdc, bj=bdc` | `_ficonv_core.c:56` vs `ficonv.c:fit_kernels_native` | 单块 4M 像素的矩阵组装 vs 1024 个 4K 小块。总计算量相同但内存局部性差异显著。C 二进制中 `fit_kernels=3.97s`，全图单块理论上不应慢 100×，但叠加 P0 问题后会走向异常路径。 |
| **P1** | 无 mask 合并/扩张 | `_ficonv_core.c` | 更多像素参与拟合（C 二进制会 mask 掉 kernel 边缘和 bad 像素），结果不同 |
| **P2** | 无迭代拟合循环 | `ficonv.c:fit_kernels` | 不做 outlier 剔除，系数可能受异常像素污染 |
| **P2** | profiling 开销 | `kernel-base.c` (fitsh_cython) | 每像素的 `clock()` 调用在 profile block 外，开销极小。但 solve_gauss 处的 `clock()` 是单次调用，无影响 |
| **P2** | ficonv.c 子函数缺失 | `_ficonv_core.c` | 缺少 foreground mask / pre-spline / subtract / add / stamps 等功能，但非速度影响 |

---

## 五、总结

Cython 版 `_ficonv_core.c` 是原版 `ficonv.c` 的**大幅度简化重写**，并非薄封装。它直接跳过了 `fit_kernels()` / `fit_kernels_native()` 的完整流程，改为直接调用底层 `fit_kernel_poly_coefficients_block`。两个最关键的差异：

1. **`kernel_init_images()` 缺失** — 高斯 kernel 的 `k->image` 从未被分配，`convolve_point` 应 crash 而非超时
2. **块大小 1 vs 32** — 改变矩阵组装的分块策略

若当前 `.so` 版本确实"跑不出来"且"非常慢"但不崩溃，说明存在其他未追踪的初始化路径或编译优化改变了行为。建议在 `.so` 版本中先加 `assert(k->image != NULL)` 确认 image 是否真的为 NULL，再做下一步判断。
