# lfit 移植方案

**日期**: 2026-06-18
**状态**: 方案审批中，尚未执行

## 架构决策

lfit 与 fiign/ficombine 模式一致——创建独立文件，编译进 pyfitsh .so，Cython 包装。
区别：lfit 核心约 4000 行，采用"最小修改 origincode + 新写 interface"策略。

## 文件结构

```
pyfitsh/lfit/
├── lfit.c              ← origincode 代码，最小修改
├── lfit.h              ← 扩展：暴露给 interface 的函数声明
├── lfit_interface.c    ← 我们新写：fitbuf、lfit_apply()、Python 入口
├── lfit-builtin.c      ← 内置函数注册表（无修改）
├── lfit-builtin.h      ← 内置函数声明
├── lfit-info.c         ← 帮助信息（后续可能删除）
├── lfit-info.h
├── lfit-stubs.c        ← stub：I/O/CLI 函数空实现（后续可能删除）
├── iof.h               ← stub：I/O 函数声明（后续可能删除）
├── scanarg.h           ← stub：CLI 解析声明（后续可能删除）
├── longhelp.h          ← 长帮助宏/类型
├── xfunct.c            ← 天文扩展函数（暂不编译，需 astro.h）
└── xfunct.h            ← 天文扩展声明
```

## lfit.c 需要的最小修改

9 个 `fit_*()` 函数已经是 non-static（main() 调用它们），lfit_interface 可以直接引用：

| 函数 | 行号 | 方法 |
|------|------|------|
| `fit_linear_or_nonlinear()` | 1898 | CLLS / NLLM / LMND |
| `fit_markov_chain_monte_carlo()` | 2856 | MCMC |
| `fit_map_chi2()` | 3314 | MCHI |
| `fit_error_monte_carlo_estimation()` | 3820 | EMCE |
| `fit_extended_markov_chain_mc()` | 4143 | XMMC |
| `fit_fisher_matrix_analysis()` | 4628 | FIMA |
| `fit_downhill_simplex()` | 4876 | DHSX |

唯一需要暴露的辅助函数（~3-5 个 `static` → 去掉 `static`）：
- PSN 表达式编译相关
- PSN symbol 创建/注册
- 数据行解析辅助函数

## lfit_interface.c 设计

### 核心思路：用 fmemopen() 避免修改 fitout

lfit.c 中 `fitout` 使用 FILE*，不修改其定义。lfit_interface.c 使用 POSIX `fmemopen()` 创建内存 FILE*：

```c
FILE *fw = fmemopen(buf, size, "w");
// 传给 fitout，fit_*() 函数的所有 fprintf(off->f_write, ...) 自动写入 buf
// 结束时 fclose(fw)，从 buf 提取结果
```

`fmemopen()` 是 POSIX 标准，macOS/Linux 均支持。这样 lfit.c 零修改。

### 需要实现的函数

| 函数 | 说明 | 估计行数 |
|------|------|----------|
| `lfit_init_data()` | 从 Python 数组构建 lfitdata + datablock | ~150 |
| `lfit_apply()` | PSN 编译 + lfit_init + fit_*() + 提取结果 | ~150 |
| `lfit_result_free()` | 释放 lfit_result | ~20 |

### fitout 处理

```c
// lfit_interface.c 中
typedef struct {
    char *data_buffer;   size_t data_len;
    char *used_buffer;   size_t used_len;
    char *rejd_buffer;   size_t rejd_len;
    char *all_buffer;    size_t all_len;
    char *expr_buffer;   size_t expr_len;
    char *vval_buffer;   size_t vval_len;
} fitbuf;

fitout *create_fitout_from_buffers(fitbuf *fb);  // 用 fmemopen 创建 FILE*
void    extract_fitout_results(fitout *off, fitbuf *fb, lfit_result *r);
```

### API 设计

```c
// lfit_apply 入口
lfit_result *lfit_apply(
    // 数据：每个 datablock
    int ndb,                        // datablock 数量
    double **db_x, int *db_nrow, int *db_ncol,  // 独立变量 (扁平数组)
    double **db_y,                  // 因变量
    double **db_err,                // 误差
    int *db_errtype,                // 0=error, 1=weight

    // 模型和参数
    char **db_model,                // 模型函数表达式
    char **db_dep,                  // 因变量表达式（NULL = 评估模式）
    char **db_err_expr,             // 误差/权重表达式
    char **db_col_names,            // 列名
    char **db_keys,                 // datablock key

    // 拟合变量
    char *variables,                // var=init[:err],...

    // 拟合设置
    int fit_method,                 // FIT_METHOD_*
    double sigma, int niter,        // 外点拒绝
    int mc_iterations,              // MC 迭代次数

    // 各方法微调参数
    double lm_lambda, double lm_lambda_mpy, int lm_max_iter,
    int use_numeric_derivs,
    int dhsx_use_fisher,
    // ... 等等
);
```

### 返回结构体

```c
typedef struct {
    int    nvar;                    // 变量数量
    char   **var_names;             // 变量名列表
    double *var_values;             // 最佳拟合值
    double *var_errors;             // 误差估计
    double *covariance;             // nvar×nvar 协方差矩阵（扁平）
    double residual;                // 拟合残差
    double red_chi2;                // 约化 chi^2
    
    int    nused;                   // 参与拟合的行数
    int    nrejd;                   // 被拒绝的行数
    double *used_mask;              // 参与/拒绝标记
    
    char   *history;                // 命令历史字符串
    char   *expr_text;              // 最佳拟合表达式文本
    char   *vval_text;              // 变量=值文本
} lfit_result;
```

## Cython 接口

```python
class Lfit:
    def fit(self, x_data, y_data, model, variables,
            errors=None, method='clls',
            sigma=3.0, niter=5,
            lm_lambda=0.001, lm_max_iter=100,
            mc_iterations=10000,
            ...):
        """→ {'params': {name: value}, 'errors': {name: err},
              'cov': ndarray, 'residual': float, 'history': str}"""

    def evaluate(self, expressions, variables, x_data):
        """仅表达式求值模式 (FIT_METHOD_NONE)"""
```

## 编译

setup.py 需要：
1. `lfit.c` + `lfit-builtin.c` + `lfit_interface.c` → 编译
2. 移除 `lfit-stubs.c` 和相关 stub 头文件
3. 保留 `-D_FITSH_SOURCE` 编译标志

## 测试

利用 `fitsh_test_data/test_lfit/` 目录中的测试数据：
1. 运行 CLI lfit 生成参考输出
2. 编写 Python 测试脚本
3. 逐项对比参数值、误差、协方差矩阵

## 已知限制

- `lfit-info.c` / `lfit-stubs.c` 在最终版本中移除（stub 仅用于过渡编译验证）
- xfunct 天文扩展暂时不启用（需 astro/astro.h 依赖链）
- 动态加载扩展模块不支持（macOS 有 dlopen 但 pyfitsh 不引入）
- bitpix 相关的 FITS 量化不适用（Python 层处理）
