# fiphot 进展

## 已验证

| 功能 | 状态 | 详情 |
|------|------|------|
| `photometry()` | ✅ | 与 CLI 完全一致 |
| `photometry_from_raw()` 无 kernel | ✅ | 49 星 flux diff < 0.6%，bg diff < 0.84。差异来源：CLI ASCII raw 精度截断 |
| ficonv `fit()` | ✅ | convolved/subtracted 与 CLI 一致，kernel list coeff 与 CLI 一致 |

## 待实现：`photometry_from_raw()` 带 kernel

对应 CLI 命令：
```
fiphot -s sub_spatial1.fits --input-raw-photometry ref_rawcheck.raw --input-kernel kernel_spatial1.txt ...
```

### 第一层（源头）：C 层 kernel list 返回字符串缺少 basis

| 项 | 详情 |
|------|------|
| **位置** | `ficonv/ficonv_core.c:826-840` |
| **问题** | 构建 `klist_out` 返回字符串时，只输出了 `kernel=`, `type=`, `order=`, `hsize=`, `sigma=`, `coeff=`，**缺少 `b=%d,%d`（basis）** |
| **影响** | `Fitrans_ConvolutionOp.fit()` 在 `ficonv_internals.pxi:203` 调用 `parse_kernel_list(raw)` 解析后，kernel dict 中无 `bx`/`by` 字段 |
| **对比** | 写入文件的路径 `kernel-io.c:107` 用 `fprintf(fw,"basis = %d,%d\n",...)` 正确输出 basis，但 C 层 `klist_out` 用另一套格式化代码，缺失 basis |

### 第二层（解析）：`parse_kernel_list` 缺少 bx/by 解析

| 项 | 详情 |
|------|------|
| **位置** | `ficonv/ficonv_internals.pxi:2-45` |
| **当前解析字段** | `kernel`, `type`, `order`, `hsize`, `sigma`, `coeff` |
| **缺失** | 无 `b=` 或 `basis=` 解析分支 |
| **结果** | 即使源头修好，解析层也不认识 basis 字段 |

### 第三层（传递）：`Fitrans_ConvolutionOp.fit()` 和 `Fiphot.photometry_from_raw()` 传递链

**3a. Fitrans_ConvolutionOp.fit() prefit 传递**

| 项 | 详情 |
|------|------|
| **位置** | `ficonv/ficonv_internals.pxi:151-180` |
| **问题** | 从 `kernel_dict['kernels']` 提取数据传给 C 层时，只传了 `types`, `orders`, `ncoeffs`, `coeffs`，**没传 `hsize`, `sigma`, `bx`, `by`** |
| **影响** | 即使用户从外部传入带 basis 的 kernel_dict，C 层也无法收到 |

**3b. Fiphot.photometry_from_raw() 缺失 kernel_dict 路径**

| 项 | 详情 |
|------|------|
| **位置** | `fitsh_cy.pyx:1382-1384`（签名）+ `fiphot_core.c:1807-1826`（C 函数） |
| **当前支持** | 仅 `kernel_spec` 字符串参数（走 `create_kernels_from_kernelarg`，只建空框架不设系数） |
| **缺失** | 无 `kernel_dict` 参数，无法将 ficonv 拟合好的 kernel 数据传入 C 层 |
| **缺失的 C 函数** | `kernel_info_read_dicts` — 此函数**当前不存在**于任何 `.c`/`.h` 文件中 |

### 第四层（副作用）：`kernel_init_images` 对 gaussian kernel 的误处理

| 项 | 详情 |
|------|------|
| **位置** | `algorithms/kernel-io.c:419-455` |
| **问题** | 当 `bx=0 && by=0`（默认值/缺失 basis 时的值），所有 6 个 gaussian kernel 都满足 `bx%2==0 && by%2==0` 条件 |
| **误触发** | Line 437-438：第一个满足条件的 gaussian kernel 被错误当作"偶数 basis"的首个核 |
| **误操作** | Line 440-445：后续 gaussian kernel 中心像素被减 1.0（`k->image[k->hsize][k->hsize]-=1.0`） |
| **结论** | 即使 bx/by 只影响这个逻辑，但若所有核的 bx/by 都是 0，会产生不正确的中心像素修正 |

## 函数命名清理

| 变更 | 状态 |
|------|------|
| `_fiphot_*` → `fiphot_*` | ✅ |
| `_prof_*` → `prof_*` | ✅ |
| `_min4`/`_max4`/... `_integrate` → 去下划线 | ✅ |
| `_flat` → `_cy` | ✅ |
| `_parse_kernel_list` → `parse_kernel_list` | ✅ |
| `ConvolutionOp` → `Fitrans_ConvolutionOp`（ficonv 侧） | ✅ |
| 迁移 `parse_kernel_list` + `Fitrans_ConvolutionOp` 到 `ficonv_internals.pxi` | ✅ |

## 调试中临时代码（最终需清理）

| 文件 | 临时代码 | 用途 |
|------|----------|------|
| `fiphot/fiphot_core.c` | `dump_helper.h` include + 输入/输出 dump 调用 | CLI vs Cython 对比验证 |
| `fiphot/dump_helper.h` | dump 辅助函数 | **最终删除** |
| `tmp/origincode/src/fiphot.c` | `dump_helper.h` include + dump 调用 | CLI 侧对比验证 |
| `tmp/origincode/src/dump_helper.h` | dump 辅助函数副本 | **最终删除** |
