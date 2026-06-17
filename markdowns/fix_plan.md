# fitsh_cython 修复方案（逐行精确对照）

## 概览：6 个差异

| # | 文件 | 差异类型 | 需修复 |
|---|------|---------|--------|
| 1 | `fitsh_cython/fiphot/fiphot_core.c` | 末尾新增 | **已完成** |
| 2 | `fitsh_cython/fiphot/fiphot_core.h` | 新增声明 | 否 — 正确修改 |
| 3 | `fitsh_cython/fiphot/fiphot_cdef.pxd` | 新增声明 | 否 — 正确修改 |
| 4 | `fitsh_cython/fitsh_cy.pyx` | unzip 恢复后新增方法 | **已完成** |
| 5 | `fitsh_cython/fitsh_cy.c` | 被删除 → 重新生成 | **已完成** |
| 6 | `tmp/origincode/src/fiphot.c` | 内联 dump 代码 | 已注释 |

## 修复步骤汇总

### 步骤 A（已完成）：`read_raw_photometry_cy` 替换旧函数
### 步骤 B（已完成）：`fiphot_photometry_from_raw_cy` 编排函数 + pyx 方法
### 步骤 C（已完成）：重新编译
### 步骤 D（已完成）：清理 debug 代码

## 验证结果（20260616T104000Z）

CLI `--input-raw-photometry` vs Cython `photometry_from_raw`：
- 49 颗星，flux 差异 < 0.6%（最大 0.45%，源自 CLI ASCII raw 文件精度截断）
- bg 差异 < 0.84

## 变更历史

### 20260616T102000Z — 重构 `read_raw_photometry_cy`
### 20260616T103500Z — 添加编排函数 + pyx 方法
### 20260616T104000Z — 编译验证通过
### 20260616T105000Z — 下划线前缀清理（3 文件，17 函数 + 10 变量）
### 20260616T110000Z — `_flat` → `_cy` 重命名（22 文件），结论存档
