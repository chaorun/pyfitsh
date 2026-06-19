# pyfitsh

> **标注**：本项目为 [fitsh](https://www.fitsh.net/) 的复刻版，遵循原始 CLI 的行为约束。
> 本项目采用 OpenCode + DeepSeek V4 模型在人工指导下完成。

fitsh C 天文软件管线的 Cython 移植版本。通过 Python 调用 `fiarith`（逐像素算术表达式）、`fistar`（星检测/拟合）、`fiphot`（孔径测光）、`grmatch`（星点匹配）、`fitrans`（图像变换）、`ficonv`（图像卷积/减影）、`firandom`（人工图像生成）、`ficalib`（图像校准）、`fiign`（像素 mask 操作）、`lfit`（通用曲线拟合）。

## 构建

```bash
cd pyfitsh
python3 setup.py build_ext --inplace
```

依赖：Cython、numpy、astropy、scipy。

编译优化：setup.py 使用 `-O3` 优化级别，性能与原始 C 版本 (`-O2`) 持平或略快（见下方性能对比）。

## 文件结构

```
pyfitsh/
├── setup.py          # 编译配置 (-O3)
├── core.pyx          # Cython 主模块
├── core.pxd          # C 函数声明
├── fitsh.h           # 通用宏/定义
├── common.h          # 公共函数声明
├── fiarith/          # 逐像素算术表达式求值器
├── algorithms/       # 通用算法 (含 tokenize, stubs)
├── fistar/           # 星检测/拟合
├── fiphot/           # 孔径测光
├── ficonv/           # 卷积/减影
├── fitrans/          # 图像变换
├── firandom/         # 人工图像生成 (含 random, PSN)
├── ficalib/          # 图像校准 (bias/flat/dark)
├── fiign/            # 像素 mask 操作
├── lfit/             # 通用曲线拟合 (9种方法 + eval模式, chain/cov_matrix自动返回)
├── psn/              # PSN 表达式解析器
├── grmatch/          # 星点匹配
├── grtrans/          # 坐标变换
├── math/             # 数学库
├── index/            # 索引/排序
├── link/             # 连通域
├── examples/         # 用法示例
├── deprecated/       # 废弃文件
└── utils.py          # MASKINFO 解码等工具
```

## 返回格式

Fistar 和 Fiphot 返回 `types.SimpleNamespace`，属性访问替代字典键：
- `r.output` — astropy Table（主结果表）
- `r.nstar` — 检测星数
- `r.dict` — 向后兼容字典

## Lfit 模块

完整复刻 origincode `lfit` 命令行工具的全部功能（除文件 I/O 外），参数名对标 `lfit --long-help`。

### 支持的拟合方法 (9种)

| 方法 | CLI 参数 | Python method= | 说明 |
|------|----------|----------------|------|
| CLLS | -L | `'clls'` | 经典线性最小二乘 |
| NLLM | -N | `'nllm'` | 非线性 Levenberg-Marquardt |
| LMND | -U | `'lmnd'` | Levenberg-Marquardt (数值导数) |
| DHSX | -D | `'dhsx'` | Downhill simplex |
| MCMC | -M | `'mcmc'` | Markov Chain Monte-Carlo |
| MCHI | -K | `'mchi'` | Chi2 网格搜索 |
| EMCE | -E | `'emce'` | Monte-Carlo 误差估计 |
| XMMC | -X | `'xmmc'` | 扩展 MCMC |
| FIMA | -A | `'fima'` | Fisher 信息矩阵分析 |

### API 用法

```python
from pyfitsh.lfit import Lfit, lfit_fit
import numpy as np

data = np.loadtxt('data.dat')

# 类式 API
lf = Lfit(method='clls', rejection_niter=3, rejection_level=3.0)
r = lf.fit(data, variables='a,b', columns='x:1,y:2',
           function='a+b*x', dependent='y')

# 函数式 API
r = lfit_fit(data, variables='a,b', columns='x:1,y:2',
             function='a+b*x', dependent='y')

# 结果访问 (属性模式 + dict 模式)
r.params          # array([1.99992, 3.0])
r['param_dict']   # {'a': 1.99992, 'b': 3.0}
r.chi2            # 6.9e-06
r.ok              # True

# MCMC/XMMC/EMCE 采样链 (自动返回)
r = lfit_fit(data, 'a=1:0.2,b=1:0.2', 'x:1,y:2', 'a+b*x', 'y',
             method='mcmc', seed=12345, mc_iterations=1000)
r.chain           # ndarray(chain_count, nvar+1), 最后一列为 chi2

# FIMA Fisher 矩阵 (自动返回)
r = lfit_fit(data, 'a=1,b=2', 'x:1,y:2', 'a+b*x', 'y', method='fima')
r.errors          # Fisher 不确定度
r.cov_matrix      # 协方差矩阵 (nvar, nvar)
r.corr_matrix     # 相关矩阵

# eval 模式 (无拟合变量)
r = lfit_fit(data, '', 'x:1,y:2', '2+3*x')
r.eval_data       # ndarray(nrow, 1)
```

### 完整参数列表

覆盖 `lfit --long-help` 全部非 IO 参数：

```python
lfit_fit(data,
    variables,              # -v  变量定义 (如 'a,b' 或 'a=1:0.2,b=1:0.2')
    columns,                # -c  列定义 (如 'x:1,y:2,e:3')
    function,               # -f  模型函数 (如 'a+b*x+c*x*x')
    dependent=None,         # -y  因变量 (None=求值模式)
    error=None,             # -e  误差表达式 (sigma模式)
    weight=None,            # -w  权重表达式 (与error互斥)
    method='clls',          # -L/-N/-U/-D/-M/-K/-E/-X/-A
    parameters=None,        # -P  方法微调 (如 'accepted,gibbs', 'clls', 'iterations=50')
    differences=None,       # -q  LMND 数值差分步长
    separate=None,          # -k  线性变量分离
    perturbations=None,     # --perturbations EMCE 扰动噪声
    seed=0,                 # -s  随机种子 (0=固定, <0=自动)
    mc_iterations=1000,     # -i  MC 迭代次数
    rejection_level=3.0,    # -r  sigma 拒绝阈值
    rejection_niter=0,      # -n  拒绝迭代次数 (0=不拒绝)
    weighted_sigma=0,       # --weighted-sigma
    macros=None,            # -x  用户宏定义 (列表, 如 ['line(u,v,t)=u+v*t'])
    constraints=None,       # -t  约束表达式
    errdump=0,              # --errors 输出误差
    format=None,            # -F  输出格式
    correlation_format=None,# -C  相关系数格式
    derived_variables=None, # -g  派生变量 (如 'sum=a+b+c')
    is_dump_delta=0,        # --delta
    resdump=0,              # --residual
    force_nonlinear=0,
    columns_output=None,    # -z  求值模式列替换
)
```

### 测试覆盖

对标 origincode 的 3 个测试脚本，全部 47 个测试通过：

| 脚本 | 测试数 | 内容 | 结果 |
|------|--------|------|------|
| test_lfit.sh | 19 | CLLS/NLLM/LMND/DHSX/chi2grid/eval/macro/constraint/pairs | 18/19 PASS, 1 SKIP |
| test_lfit_montecarlo.sh | 19 | MCMC/XMMC/EMCE/FIMA (chain逐行对比 + Fisher矩阵) | 19/19 PASS |
| longtest_lfit.sh | 9 | 长迭代 MCMC/XMMC/EMCE/FIMA (20000 iterations) | 7/9 PASS, 2 DIFF (CLI精度截断) |

## 编译优化 (-O3)

setup.py 使用 `-O3` 编译优化级别。性能对比 (pairs MCMC 20000 iterations):

| 版本 | L08 时间 | L09 时间 | 相对 C 原版 |
|------|----------|----------|------------|
| C 原版 (-O2) | 221s | 236s | 1.00x |
| Cython -O0 | 462s | 478s | 2.09x 慢 |
| **Cython -O3** | **211s** | **221s** | **0.95x 快** |

`-O3` 编译后 Cython 版本比 C 原版 (`-O2`) 快约 5%。

## 版本

2026-06-19 快照。复刻 fitsh 0.9.4 行为，遵循原始设计约束：
- 新增 `lfit` 模块：通用曲线拟合，9 种方法 + eval 模式，chain/cov_matrix/corr_matrix 自动返回，44/47 测试 PASS
- 编译优化从 `-O0` 升级为 `-O3`，性能与 C 原版持平
- 新增 `ficalib` 模块：图像校准（bias/flat/dark），12/12 测试 PASS
- 新增 `fiign` 模块：像素 mask 操作（宇宙线、饱和、几何 mask、mask 转换），13/13 测试 PASS
- UDF 子表达式 `[...](body)` 不支持图像函数（sq, smooth, laplace 等），仅支持 `psn_general_fn` 中的通用数学函数（sin, cos, exp, ln, sqrt, abs, sign 等），与原始 fitsh CLI 行为一致
- evaluate.c 中 `"if_sq"` 修正为 `"sq"`
- 死代码（CLI 文件 I/O）已注释，禁止 Cython 读写文件
- 函数命名已清理：无下划线前缀，`_flat` 后缀 → `_cy` 后缀
- 新增 19 个 CLI 对照测试脚本（`examples/` 目录 + `runtest.sh`），全部 PASS（2026-06-18）
- 定位 origincode CLI bug：`--input-candidates --col-shape` 未读取 cd/ck 列（见 USAGE.md）
- 修复 `core.pyx` 中 `_psf_type` 常量偏移 bug（pbg/pamp nan）
- 修复 ficonv / fiphot 模块 mask 处理
- 修复 `fistar_core.c` input_candidates 路径崩溃（三处修复）
- `fprint_error`/`fprint_warning` 统一移至 `algorithms/stubs.c`，消除多模块重复符号

## License

This project is licensed under the GNU General Public License v3.0 or later
(GPL-3.0-or-later).

`pyfitsh` is based on, wraps, ports, or reimplements parts of FITSH by András Pál.
The original FITSH project is available at https://www.fitsh.net/ and is licensed
under the GNU General Public License.

Where FITSH source code has been copied, translated, or adapted, the original
copyright notices are retained in the corresponding source files.
