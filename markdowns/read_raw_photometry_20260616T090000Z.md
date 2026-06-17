# read_raw_photometry: CLI 原版 vs Cython 版本 — 逐字段对比

## 对比对象

| | CLI 原版 | Cython 版本 |
|---|---|---|
| **文件** | `tmp/origincode/src/fiphot-io.c` | `fitsh_cython/fiphot/fiphot_core.c` |
| **函数** | `read_raw_photometry()` | `fiphot_photometry_from_raw_cython()` |
| **行号** | 275-374 | 1809-1840 |
| **输入** | 文本流 (FILE*)，逐行 ASCII | numpy float64 数组 (nstar × nap × 12) |

---

[旧版逐字段对比内容与 read_raw_photometry.md 一致]
