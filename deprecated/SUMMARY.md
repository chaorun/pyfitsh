# grmatch Cython 化总结

## 目标

将 C 语言编写的 `grmatch` 天文星表匹配程序的 `do_pointmatch()` 核心函数暴露给 Python 调用，用 Cython 编译为自包含的 `.so` 扩展模块。

## 架构演进

### 第一阶段：代码分离（C 层面）

原始 `grmatch.c`（2111 行）包含 `main()` 和所有库函数，无法直接复用。

```
grmatch.c  ──拆分──→  grmatch.c (仅 main, 583行)
                      grmatch_lib.c (全部库函数, 1439行)
                      grmatch_lib.h (类型 + 函数声明, 110行)
```

- `grmatch_lib.h` 声明 `matchpointtune`、`iline`、`cphit` 等类型
- `gratch_lib.c` 定义 `do_pointmatch()`、`read_match_data_points()` 等
- Makefile 新增 `grmatch_lib.c` 编译

### 第二阶段：dump 机制与 ctypes 验证

`do_pointmatch()` 内置 `dump_pointmatch_stage()`：在函数入口/出口将全部内存数据写入 `dump/*.bin`（二进制精度）+ `dump/*.txt`（`%g` 文本）。

**grmatch_simple.c** — 从 `dump/` 加载 `.bin`，重建结构体，重新调用 `do_pointmatch()`，输出到 `dump_out/`。用于验证 dump 的完整可重放性。

**BIN vs TXT 输入对比实验**：
- TXT 的 `%g` 截断（6 位有效数字）导致 78% 配对不同
- 但在 500×500 像场上实际坐标偏差仅 **0.0005 像素**，无实质影响

### 第三阶段：ctypes + dylib（Python 桥接）

```
grmatch_bindings.py  ── ctypes 结构体定义 + CDLL 加载
grmatch_reader.py    ── 列配置 + 文件读取 + iline 构建
```

- `Matchpointtune` (136B), `Iline` (72B), `Cphit` (16B) 等全部手工对齐
- 通过 `ctypes.CDLL("libgrmatch.dylib")` 调用 C 函数
- 验证：ctypes 结果与 C grmatch_simple **完全一致**（12 个 vfits 系数精确至 1e-14）

### 第四阶段：Cython（编译为 Python 原生扩展）

**目标**：消除 `ctypes` 的 FFI 开销、手工结构体对齐、`.dylib` 分离部署。

**三步策略**：
1. 从 `grmatch_lib.c` 提取 `do_pointmatch()` + 3 个 static 辅助函数 → `_core.c`
2. 新增 `do_pointmatch_flat()` 包装：接受平坦 `double[]` 数组代替 `iline` 结构体
3. Cython `.pyx` 包装为 `Matcher` 类

**依赖裁剪**：`do_pointmatch` 只需 10 个 `.c` + 15 个 `.h`（纯数学核心），不需要 `io/`、`ui.c`、`longhelp.c` 等。

**自包含化**：所有依赖复制到 `grmatch_cython/` 目录，零上级目录引用。

## 最终文件结构

```
grmatch_cython/          ← 独立、可移植
├── run_match.py         # python3 run_match.py 一键测试
├── setup.py             # python3 setup.py build_ext --inplace 编译
├── _core.c / _core.h    # do_pointmatch 提取 + 平坦数组包装
├── grmatch_cy.pxd        # Cython 类型声明
├── grmatch_cy.pyx        # Matcher 类 + MatchResult 类
├── grmatch_cy.*.so       # 编译产物（169KB）
├── *.h                  # grmatch_lib.h, transform.h, longhelp.h, fitsh.h, config.h
├── math/                # trimatch, cpmatch, poly, polyfit, delaunay, spmatrix,
│   │                    #   convexhull, tpoint (各 .c + .h)
│   └── fit/lmfit.*      # Levenberg-Marquardt 求解器
└── testdata/            # temp_ref.txt, temp_obs.txt, dump/
```

## ctypes vs Cython 对比

| | ctypes | Cython |
|---|---|---|
| 导入方式 | `ctypes.CDLL("libgrmatch.dylib")` | `from grmatch_cy import Matcher` |
| 类型安全 | 运行时 FFI，错配难排查 | 编译时 `.pxd` 检查 |
| 接口 | C 结构体指针，手工 `byref` | `m.set_reference(x, y, ord)` |
| 构建 | `make dylib` | `python3 setup.py build_ext` |
| 部署 | `.dylib` + `.py` 分离 | 单文件 `.so` |
| 结构体 | 手工 Python 类对齐 | `cdef extern from ".h"` 自动 |

## 关键踩坑记录

1. **Cython `__cinit__` vs `__init__`**：cdef class 的 C 级构造用 `__cinit__`，Python 属性必须在 `__init__` 中赋值，且需声明为 `cdef public`

2. **Cython 生成文件名冲突**：`.pyx` 名与手写 `.c` 同名时 Cython 拒绝覆盖 → 改手写 C 文件为 `_core.c`

3. **setup.py 路径**：`cythonize()` 从 CWD 解析相对路径 → 统一使用 `os.path.join(HERE, ...)` 绝对路径

4. **`%g` 截断传播**：输入坐标 6 位有效数字截断 → 三角形匹配算法重分配配对（78% 变化）→ 但低阶变换系数几乎不变，实际坐标偏差 0.0005px

5. **自包含化陷阱**：看似只需 `do_pointmatch`，但编译器需要解析所有顶层 `#include`，间接依赖链拉入了 `longhelp.h → fitsh.h → config.h` 等纯类型头文件

## 调用示例

```python
from grmatch_cy import Matcher

# 准备数据 (与 C read_match_data_points 等价)
ref_x, ref_y, ref_o = read_points("testdata/temp_ref.txt", "2,3", "-4")
inp_x, inp_y, inp_o = read_points("testdata/temp_obs.txt", "2,3", "-4")

# 运行匹配
m = Matcher(order=2, maxdist=2.0, unitarity=0.01)
m.set_reference(ref_x, ref_y, ref_o)
m.set_input(inp_x, inp_y, inp_o)
r = m.run()

# 结果
print(r.nhit)        # 498
print(r.vfits_dx)    # [-2.27037, 1.00003, 0.000152, ...]
print(r.stats)       # {'nsigma': 0.59055, 'unitarity': 0.000111, ...}
```

## 验证结果

所有路径收敛于相同变换系数：

| 路径 | 匹配数 | vfits | 状态 |
|------|--------|-------|------|
| C grmatch 命令行 | 498 | baseline | — |
| dump → C grmatch_simple | 498 | 精确一致 | OK |
| dump → Python ctypes | 498 | 精确一致 | OK |
| raw files → Python ctypes | 498 | 精确一致 | OK |
| **raw files → Cython** | **498** | **精确一致** | **OK** |
