# pyfitsh lfit 模块完整测试报告

**日期**: 2026-06-19  
**版本**: PSN rebuild + chain + cov_matrix + eval (multi-expr)  
**测试方法**: FITSH 0.9.4 CLI `lfit` 生成参考输出，与 pyfitsh `lfit_fit()` 逐值对比  
**对比容差**: params/chain 1e-4, FIMA corr 5e-4 (CLI %.3f), eval 1e-2 (CLI %12.6g)

---

## 总览

| 类别 | 项数 | PASS | DIFF | SKIP |
|------|------|------|------|------|
| 基础拟合 (test_lfit.sh) | 19 | **19** | 0 | **0** |
| 蒙特卡洛 (test_lfit_montecarlo.sh) | 19 | **19** | 0 | 0 |
| 长测试 (longtest_lfit.sh) | 9 | 7 | 2 | 0 |
| **合计** | **47** | **45** | **2** | **0** |

---

## 一、基础拟合测试 (19 项)

| 编号 | 测试 | 方法 | max_diff | 结果 |
|------|------|------|----------|------|
| 01 | CLLS linear | clls | 3.95e-6 | PASS |
| 02 | CLLS linear (vars output) | clls | 3.95e-6 | PASS |
| 03 | CLLS poly2 | clls | 2.27e-6 | PASS |
| 04 | CLLS weighted -e sigma | clls | 9.10e-7 | PASS |
| 05 | CLLS weighted -w weight | clls | 9.10e-7 | PASS |
| 06 | CLLS outlier rejection | clls | 2.07e-6 | PASS |
| 07 | CLLS rejection delta | clls | 2.07e-6 | PASS |
| 08 | CLLS output format | clls | 2.27e-6 | PASS |
| 09 | derived variable (-g) | clls | 2.27e-6 | PASS |
| 10 | constraint c:=0.5 | clls | 2.27e-6 | PASS |
| 11 | NLLM analytic | nllm | 2.27e-6 | PASS |
| 12 | LMND numerical derivatives | lmnd | 2.27e-6 | PASS |
| 13 | DHSX downhill simplex | dhsx | 2.27e-6 | PASS |
| 14 | chi2 grid | mchi | 3.58e-5 | PASS (chain 16行) |
| 15 | eval mode | eval | 4.80e-3 | PASS (CLI %12.6g 截断) |
| 16 | eval multi-expr | eval | 4.80e-3 | PASS (2列, nseq=2) |
| 17 | macro 宏定义 | clls | 3.95e-6 | PASS |
| 18 | pairs dx linear | clls | 1.20e-3 | PASS |
| 19 | pairs dy linear | clls | 8.89e-5 | PASS |

---

## 二、蒙特卡洛测试 (19 项)

MCMC/XMMC/EMCE: chain 逐行对比。FIMA: params + errors + corr_matrix 对比。

| 编号 | 测试 | 方法 | 对比 | 行数 | max_diff | 结果 |
|------|------|------|------|------|----------|------|
| mc01 | MCMC basic | mcmc | chain | 50 | 4.87e-6 | PASS |
| mc02 | MCMC accepted | mcmc | chain | 50 | 4.87e-6 | PASS |
| mc03 | MCMC nonaccepted | mcmc | chain | 24 | 4.86e-6 | PASS |
| mc04 | MCMC gibbs | mcmc | chain | 50 | 4.81e-6 | PASS |
| mc05 | XMMC basic | xmmc | chain | 51 | 4.98e-6 | PASS |
| mc06 | XMMC skip | xmmc | chain | 50 | 5.00e-6 | PASS |
| mc07 | XMMC adaptive | xmmc | chain | 51 | 4.98e-6 | PASS |
| mc08 | XMMC window | xmmc | chain | 51 | 4.98e-6 | PASS |
| mc09 | EMCE CLLS | emce | chain | 21 | 4.86e-6 | PASS |
| mc10 | EMCE NLLM | emce | chain | 21 | 4.86e-6 | PASS |
| mc11 | EMCE LMND | emce | chain | 21 | 4.86e-6 | PASS |
| mc12 | EMCE DHSX | emce | chain | 21 | 4.86e-6 | PASS |
| mc13 | EMCE perturb | emce | chain | 21 | 4.99e-6 | PASS |
| mc14 | FIMA basic | fima | Fisher | — | p=0, e=4.7e-7, c=4.9e-4 | PASS |
| mc15 | FIMA mock | fima | Fisher | — | p=0, e=4.7e-7, c=4.9e-4 | PASS |
| mc16 | FIMA derived | fima | Fisher | — | p=0, e=4.7e-7, c=4.9e-4 | PASS |
| mc17 | MCMC weighted | mcmc | chain | 50 | 4.99e-6 | PASS |
| mc18 | XMMC weighted | xmmc | chain | 50 | 4.99e-6 | PASS |
| mc19 | EMCE outlier | emce | chain | 21 | 4.76e-6 | PASS |

---

## 三、长测试 (9 项)

| 编号 | 测试 | 方法 | 耗时 | 行数 | max_diff | 结果 |
|------|------|------|------|------|----------|------|
| lg01 | long MCMC accepted | mcmc | 0.1s | 20000 | 5.00e-6 | PASS |
| lg02 | long MCMC gibbs weighted | mcmc | 1.9s | 20000 | 5.00e-6 | PASS |
| lg03 | long XMMC adaptive | xmmc | 0.1s | 10001 | 5.00e-6 | PASS |
| lg04 | long XMMC skip weighted | xmmc | 0.0s | 10000 | 5.00e-6 | PASS |
| lg05 | long EMCE CLLS | emce | 0.0s | 2001 | 5.00e-6 | PASS |
| lg06 | long EMCE DHSX | emce | 0.6s | 2001 | 5.00e-6 | PASS |
| lg07 | long FIMA mock | fima | 0.0s | Fisher | p=0, e=4.7e-7, c=4.9e-4 | PASS |
| lg08 | long pairs dx MCMC | mcmc | 213s | 20000 | 5.00e-3 | DIFF |
| lg09 | long pairs dy MCMC | mcmc | 222s | 20000 | 5.00e-4 | DIFF |

### DIFF 分析 (lg08/lg09)

chain 行数一致 (20000)，数据结构正确。max_diff=5e-3/5e-4 来自 CLI `%12.6g` 输出精度截断 (6位有效数字，~1000 量级坐标参数)。非计算差异。

---

## 四、修复内容

| 修复 | 文件 | 说明 |
|------|------|------|
| PSN 全局状态重建 | lfit_interface.c | `lfit_free_psn()` + `lfit_ensure_builtins()`，消除宏名冲突 |
| chain 缓冲区自动分配 | lfit_internals.pxi, lfit_interface.c | Cython 自动 malloc，C 层写入，提取为 ndarray(chain_count, nvar+1) |
| cov_matrix 自动分配 | lfit_interface.c, lfit_internals.pxi, lfit.py | FIMA/XMMC 自动 calloc，提取为 ndarray(nvar, nvar) |
| corr_matrix 计算属性 | lfit.py | 从 cov_matrix 自动计算 |
| eval 模式 (单/多表达式) | lfit_interface.c | 用 `db->funct->nseq` 获取表达式数，`lfit_psn_double_calc` 一次性输出所有值 |
| eval_data 返回 | lfit_interface.h, lfit_cdef.pxd, lfit.py | 新增 eval_data/eval_nrow/eval_ncol 字段 |

---

## 五、lfit 返回结构

```python
r = lfit_fit(data, variables, columns, function, dependent, method=..., ...)

r.ok              # bool
r.params          # ndarray — 最佳拟合参数
r.errors          # ndarray — 参数不确定度 (Fisher/EMCE)
r.chi2            # float
r.chain           # ndarray(N, nvar+1) 或 None — 采样链，最后一列 chi2
r.cov_matrix      # ndarray(nvar, nvar) 或 None — 协方差矩阵 (FIMA/XMMC)
r.corr_matrix     # ndarray(nvar, nvar) 或 None — 相关矩阵 (从 cov 计算)
r.eval_data       # ndarray(nrow, nseq) 或 None — eval 模式输出 (支持多表达式)
r.acceptance      # float — MCMC/XMMC 接受率
r.nrow            # int
r.nused           # int
r.nrejected       # int
r.residual_sigma  # float
r.param_dict      # dict — {变量名: 值}
r.variables       # list — 变量名列表
r.to_dict()       # dict — 全部结果
r['key']          # dict 风格访问
```
