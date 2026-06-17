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

### 缺少的组件

1. **`kernel_info_read_dicts`** — C 函数，接受扁平化的核数据数组，构建 `kernellist` 结构
2. **basis 字段传递** — `kernel_dict` 中缺少 `bx`/`by`（空间变化核方向信息）
   - `_parse_kernel_list`（`fitrans_internals.pxi`）需添加 `b=` 解析
   - `ficonv_core.c` kernel list 输出需添加 `b=%d,%d`
3. **Cython 层编排** — pyx 中 `photometry_from_raw` 需支持 `kernel_dict=None` 参数，将 dict 扁平化后传给 C 层

### 技术约束

- 禁止 Cython 侧读写文件（`kernel_info_read` 不可用）
- kernel dict 由 ficonv 输出，需完整保留与 CLI txt 文件一致的信息（含 bx/by/coeff/image）
- `create_kernels_from_kernelarg` 只建核框架，不填充系数，不适用于此场景

## 函数命名清理

| 变更 | 状态 |
|------|------|
| `_fiphot_*` → `fiphot_*` | ✅ |
| `_prof_*` → `prof_*` | ✅ |
| `_min4`/`_max4`/... `_integrate` → 去下划线 | ✅ |
| `_flat` → `_cy` | ✅ |
