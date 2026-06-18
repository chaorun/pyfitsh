# fistar input_candidates 路径问题报告

**日期**: 2026-06-17 ~ 2026-06-18

## 1. 问题概述

`Fistar(input_candidates=cands).do_fistar()` 调用崩溃，访问 `stars[i].cand` 时 `cand` 为野指针。

## 2. 崩溃分析

### 2.1 崩溃位置

`fistar_search_cy +6736`（ARM64汇编）:
```asm
ldr x17, [x16, #0x128]    ; x17 = stars[i].cand (offset 296)
cbz x17, skip              ; NULL check passes
ldp w16, w17, [x17]        ; CRASH: read cand->ix, cand->iy
```

`stars[i].cand` 非 NULL 但指向无效地址 `0x260001120e3918`。

### 2.2 崩溃路径

```
input_candidates (cand_xy != NULL)
  → make_star_candidates() → cands (line 398)
  → fit_star_single_model() → stars (line 498)
  → sort stars (line 640+)
  → result filling loop → CRASH (line 700+)
```

### 2.3 根因

`fistar_core.c:392-394` 中 s/d/k 被默认值覆写（bug），且缺少边界校验。

## 3. 修复方案（已实施）

### 3.1 删除 s/d/k 覆写 (`fistar_core.c:392-394`)

```diff
- icands[ci].s = 2.0;   /* default shape */
- icands[ci].d = 0.0;
- icands[ci].k = 0.0;
```

三元表达式已提供默认值（`cand_ncol >= 3 ? ... : 2.0`），无需覆写。

### 3.2 增加 `cleanup_candlist` (`fistar_core.c:398` 后)

```diff
  make_star_candidates(icands, in_ncand, cand_radius, &img, mask, &cands, &ncand);
+ cleanup_candlist(&cands, &ncand);
  free(icands);
```

origincode 同样路径中无此调用，但添加无害（无 `marked` 候选时是 no-op）。

### 3.3 结果填充循环 wc 边界校验

```diff
- candidate *wc = ws->cand;
+ candidate *wc = ws->cand;
+ int is_wc = (wc != NULL && wc >= cands && wc < cands + ncand);
- result->ix[i] = wc ? wc->ix + 1 : 0;
+ result->ix[i] = is_wc ? wc->ix + 1 : 0;
```

## 4. 修复后状态

- **崩溃已消除**：不再 crash
- **输出行数不匹配**：CY nstar=35，CLI nstar=49

## 5. fit_star_single_model 输入/输出完整对比

为定位 35 vs 49 差异，对 CLI 和 CY 的 `fit_star_single_model` 所有参数进行了完整二进制 dump 对比。

### 5.1 对比方法

CLI 和 CY 分别在 `fit_star_single_model` 调用前后，将全部参数（含结构体）写入二进制文件：
- `dump_cli_fit_io.bin` (37769684 bytes)
- `dump_cy_fit_io.bin` (37765428 bytes)

文件结构：
```
[img header: 144B] [img data: sx*sy*8B] [mask: sx*sy*1B] [ncand: 4B]
[cands: ncand*120B] [sfp: 12B] [model/order: 8B] [nstar: 4B] [stars: nstar*304B]
```

### 5.2 对比结果

| 参数 | CLI | CY | 一致? |
|------|-----|-----|-------|
| **img header** | -- | -- | 19/144B 差 (FITS元数据) |
| sx | 2048 | 2048 | ✓ |
| sy | 2048 | 2048 | ✓ |
| bit | -32 | 0 | FITS元数据，不影响 |
| data* | 0x06fc0000 | 0x0111235a | 运行时堆地址，不影响 |
| dim | 2 | 0 | FITS元数据，不影响 |
| naxis[0] | 2048 | 0 | FITS元数据，不影响 |
| naxis[1] | 2048 | 0 | FITS元数据，不影响 |
| vdata* | 0x06fc0000 | NULL | FITS文件指针，不影响 |
| allocdata* | 0x06fc0000 | NULL | FITS释放标记，不影响 |
| curr_bscale | 1 | 0 | FITS缩放，不影响 |
| read_bscale | 1 | 0 | FITS缩放，不影响 |
| **img data** | 33,554,432B | 33,554,432B | ✓ **完全逐字节一致** |
| **mask** | 4,194,304B | 4,194,304B | ✓ **完全逐字节一致** |
| **ncand** | 49 | 49 | ✓ |
| **cands** | -- | -- | ✗ **3,063/5,880B 差** |
| cand[0].ix | 842 | 73 | ✗ |
| cand[0].iy | 40 | 1875 | ✗ |
| cand[0].cx | 842.482 | 73.532 | ✗ |
| cand[0].cy | 40.226 | 1875.288 | ✗ |
| cand[0].peak | **0.0** | 0.028 | ✗ |
| cand[0].amp | **0.0** | 0.028 | ✗ |
| cand[0].flux | **0.0** | 1.0 | ✗ |
| cand[0].sxx | 0.148 | 0.185 | ✗ |
| cand[0].sxy | **0.0** | -0.020 | ✗ |
| cand[0].nipoint | 11 | 11 | ✓ |
| **sfp** | -- | -- | ✓ **完全一致** |
| fit_flags | 15 | 15 | ✓ |
| iter_symmetric | 4 | 4 | ✓ |
| iter_general | 2 | 2 | ✓ |
| **model/order** | 2/2 | 2/2 | ✓ **完全一致** |
| **nstar** | **49** | **35** | ✗ |
| **stars** | 49×304B | 35×304B | ✗ |
| star[0].loc.gcx | 842.482 | -55051.7 | ✗ |
| star[0].loc.gcy | 40.226 | 28045.6 | ✗ |
| star[0].loc.gamp | **0.0** | 2897.1 | ✗ |
| star[0].loc.gbg | **0.0** | 2189.4 | ✗ |
| star[0].flux | **0.0** | 0.56 | ✗ |

### 5.3 关键发现

1. **img_data 和 mask 完全逐字节一致**：排除图像数据差异。
2. **sfp 和 model/order 完全一致**：排除拟合参数差异。
3. **cands 输入本质不同**：
   - CLI 候选的 **peak/amp/flux 全为 0**，sxy=0
   - CY 候选有实际值（peak≈0.028, amp≈0.028, flux=1.0, sxy≈-0.02）
   - CLI 候选坐标 (842, 40)，CY 候选坐标 (73, 1875)——**顺序完全不同**
4. **CLI 的 star[0] 所有值均为 0**：说明 CLI 的 `fit_star_single_model` 产生了未拟合的 star（通过 `convert_candidates` fallback 补入），而 CY 的 star[0] 有正常拟合值。

## 6. 当前状态

- **修复后 CY nstar=35 (已拟合的星)**
- **CLI nstar=49 (含 14 个未拟合的空 star，通过 fallback 补入)**
- **根因**: CLI 和 CY 的 `make_star_candidates` 产生了不同的候选数据（peak/amp/flux 值不同，坐标顺序不同）
- **待排查**: `read_star_candidates` (CLI) vs 直接数组构建 (CY) 的 icands 差异
