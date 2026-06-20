# gropt 移植计划

**日期**: 2026-06-20  
**工具**: gropt — 几何光学计算（geometrical optics）  
**版本**: FITSH 0.9.4, gropt 0.1.0  
**源码**: `old_c_code/gropt.c` (859行) + `old_c_code/optcalc.c` (674行)

---

## 一、独立性分析

gropt 是 fitsh 中最独立的模块——不依赖 fits、fistar、fiphot、fiarith 等任何其他模块。

核心只有两个文件：
- `optcalc.c` (674行) — 纯数学光线追踪库
- `gropt.c` (859行) — 主程序（参数解析 + 功能调用）

### 编译依赖（完整清单）

```bash
gcc -I. gropt.c optcalc.c math/matrixvector.c math/polygon.c \
    io/tokenize.c longhelp.c ui.c io/scanarg.c -lm -o gropt
```

| 文件 | 行数 | 用途 | 移植需要？ |
|------|------|------|-----------|
| optcalc.c/h | 674+96 | 核心光学计算 | **是** — 复制到 pyfitsh/gropt/ |
| math/matrixvector.c/h | 229+36 | 3D 向量矩阵 | **是** — 复制到 pyfitsh/math/ |
| math/polygon.c/h | 241+? | 多边形运算(PSF) | 已有 (pyfitsh/math/) |
| io/tokenize.c/h | — | 字符串分词 | 已有 (pyfitsh/algorithms/) |
| list.h | 150 | 链表宏(PSF) | **是** — 复制到 pyfitsh/gropt/ |
| longhelp.c/h | — | CLI 帮助 | **不需要** |
| ui.c | — | CLI 版本打印 | **不需要** |
| io/scanarg.c/h | — | CLI 参数解析 | **不需要** |
| fitsh.h | — | 版本宏 | **不需要** |

---

## 二、功能模块分析

### longhelp 参数一览

| CLI 参数 | Python 参数 | 类型 | 说明 |
|----------|-------------|------|------|
| -i, --input | optics_text (str) | 输入 | 光学系统描述文本 |
| -l, --wavelength | wavelength (float) | 公共参数 | 波长（微米） |
| -f, --focus | focus (float) | 公共参数 | 焦平面位置（毫米） |
| -x, --scale | pixel_scale (float) | 公共参数 | 像素尺度（毫米） |
| -a, --angle | angle (float/tuple) | 公共参数 | 入射角 |
| -s, --spot-aperture | aperture_radius, nrings | spot | 孔径和环数 |
| -o, --output-spot | — | spot | 返回 ndarray |
| -t, --transfer | — | transfer | 返回 dict |
| -z, --psf-half-size | psf_half_size (int) | PSF | PSF 半尺寸 |
| -p, --output-psf | — | PSF | 返回 ndarray |
| -d, --output-scad | — | export | 返回 str |
| -e, --output-eps | — | export | 返回 str |
| --speed-test | — | 测试 | 暂不移植 |

### 功能模块与代码位置

| 模块 | gropt.c 行号 | C 代码量 | 输出类型 |
|------|-------------|---------|----------|
| 光斑图 spot | 394-442 | ~50行 | 数值 → ndarray(n,2) |
| 传输矩阵 transfer | 444-459 | ~15行 | 数值 → dict |
| OpenSCAD 导出 | 462-517 | ~55行 | 文本 → str |
| EPS 导出 | 519-668 | ~150行 | 文本 → str |
| PSF 计算 | 670-856 | ~186行 | 数值 → ndarray(2H+1, 2H+1) |

---

## 三、I/O 消除方案

### 原则：全部 I/O 替换为 Python 类型

| 方向 | 原始 CLI | 移植后 |
|------|----------|--------|
| 输入 | `FILE*` → fgets 读 .opt 文件 | Python `str` → C 层字符串解析 |
| 数值输出 | `fprintf` 到文件 | `double` 数组 → numpy `ndarray` |
| 文本输出 | `fprintf` 到文件 | `sprintf` 到字符串缓冲区 → Python `str` |

### optcalc.c 的处理

optcalc.c 共 22 个函数：

**纯计算函数（18个，完整保留不动）：**
- optcalc_raytrace_reset / free / append(static)
- optcalc_get_refraction_index — Sellmeier 折射率
- optcalc_glass_refraction_precompute
- optcalc_snell_descartes — Snell 折射定律
- optcalc_reflection
- optcalc_surface_aspheric_eval / diff — 非球面求值/微分
- optcalc_surface_aspheric_ray_trace — 非球面交点+法线
- optcalc_refract_quadratic_ray_trace
- optcalc_refract_aspheric_ray_trace
- optcalc_surface_reset / is_quadratic
- optcalc_lens_ray_trace — 单透镜追踪
- optcalc_ray_trace — 完整系统追踪
- optcalc_transfer_matrix_compose
- optcalc_compute_transfer_matrix
- optcalc_reset

**字符串解析辅助函数（3个，保留，被改写版 read_optics 复用）：**
- optcalc_read_curvature — 解析 "+1/100" → double（纯 sscanf，无 I/O）
- optcalc_read_aspheric_data — 解析 "0.001:0.0002"（纯 tokenize+sscanf）
- optcalc_read_radii — 解析 "25:30"（纯 sscanf）

**FILE I/O 函数（1个，就地改写）：**
- `optcalc_read_optics(FILE *fr, ...)` → 注释掉原版，改写为 `optcalc_read_optics_from_string(const char *text, ...)`
  - 将 `fgets` 循环改为逐行解析字符串（用 `\n` 分割）
  - 三个辅助函数完全复用，不需要改
  - optcalc.h 中的声明同步更新

### gropt.c 功能搬迁到 gropt_interface.c

gropt.c 的 `main()` 中有 5 个功能块，全部搬到 gropt_interface.c 中作为库函数：

**数值输出函数**（fprintf → 写入预分配 double 数组）：
```c
// 光斑图：环形采样 + ray_trace，结果写入 out_xy[]
int gropt_spot_diagram(struct optics *opt,
    double wavelength, double aperture_radius, int nrings,
    double angle_nx, double angle_ny, double zstart, double pixel_scale,
    double *out_xy, int *out_npoints);

// 传输矩阵：compute_transfer_matrix → 两个 double
int gropt_transfer(struct optics *opt, double wavelength,
    double *out_focal_plane, double *out_eff_focus);

// PSF：三角形自适应细分 + ray_trace + polygon 面积
int gropt_psf(struct optics *opt,
    double wavelength, double aperture_radius, int nrings,
    double angle_nx, double angle_ny, double zstart, double pixel_scale,
    int psf_hsize, double *out_psf);

// 单条光线追踪路径
int gropt_raytrace(struct optics *opt, double wavelength,
    double x0, double y0, double z0,
    double nx, double ny, double nz,
    double *out_points, int *out_npoints);
```

**文本输出函数**（fprintf(fw,...) → sprintf 追加到动态字符串缓冲区）：
```c
// OpenSCAD 3D 模型文本
int gropt_export_scad(struct optics *opt, char **out_str, size_t *out_len);

// EPS 平面图文本（含可选光线追踪）
int gropt_export_eps(struct optics *opt,
    double aperture_radius, int nrings,
    double angle_nx, double angle_ny, double zstart,
    double wavelength,
    char **out_str, size_t *out_len);
```

文本输出的实现方式——动态字符串缓冲区：
```c
typedef struct { char *data; size_t len; size_t cap; } strbuf;
static void strbuf_appendf(strbuf *sb, const char *fmt, ...);
// 所有 fprintf(fw, ...) 改为 strbuf_appendf(&sb, ...)
// 最终 *out_str = sb.data; *out_len = sb.len;
```

---

## 四、目录结构

```
pyfitsh/gropt/
├── optcalc.c            # 从 old_c_code 复制，注释掉 read_optics(FILE*)，
│                        # 新增 read_optics_from_string(const char*)
├── optcalc.h            # 同步更新声明
├── gropt_interface.c    # 新写：spot/transfer/psf/raytrace/scad/eps 库函数
├── gropt_interface.h    # 接口头文件
├── gropt_cdef.pxd       # Cython 声明（结构体 + 函数）
├── gropt.pxi            # Cython 实现（调用 C → 返回 ndarray/str）
├── gropt.py             # Python API（Optics 类 + 便捷函数）
└── list.h               # 链表宏（PSF 功能需要）

pyfitsh/math/
├── matrixvector.c       # 新增（从 old_c_code/math/ 复制）
├── matrixvector.h       # 新增
├── polygon.c            # 已有
└── polygon.h            # 已有
```

---

## 五、Python API 设计

```python
from pyfitsh.gropt import Optics

# 从字符串加载光学系统
opt = Optics.from_string("""
glass BK7 1.03961212:0.00600069867:0.231792344:0.0200179144:1.01046945:103.560653
lens  50.0  5.0  25.0  +1/100  -1/200  BK7
focal 120.0
""")

# 光斑图 → ndarray(n, 2)
spots = opt.spot(wavelength=0.6, aperture_radius=10.0, nrings=5,
                 angle=0.01, pixel_scale=0.015)

# 传输矩阵 → dict
tm = opt.transfer(wavelength=0.6)
# {'focal_plane': 120.5, 'effective_focus': 500.0}

# PSF → ndarray(2H+1, 2H+1)
psf = opt.psf(wavelength=0.6, half_size=3, pixel_scale=0.015,
              aperture_radius=10.0, angle=0.01)

# 光线追踪路径 → ndarray(npoints, 3)
trace = opt.raytrace(wavelength=0.6, x0=0, y0=5, z0=-0.1,
                     nx=0, ny=0, nz=1)

# 导出文本 → str
scad_text = opt.to_openscad()
eps_text = opt.to_eps(aperture_radius=10, nrings=3, wavelength=0.6)

# Python 层写文件（不在 C 层）
with open("system.scad", "w") as f:
    f.write(scad_text)
```

---

## 六、实施步骤

### 阶段 1：基础搭建
1. 创建 `pyfitsh/gropt/` 目录
2. 复制 `math/matrixvector.c/h` → `pyfitsh/math/`
3. 复制 `optcalc.c/h` → `pyfitsh/gropt/`
4. 复制 `list.h` → `pyfitsh/gropt/`
5. 在 optcalc.c 中：注释掉 `optcalc_read_optics(FILE*)`，新增 `optcalc_read_optics_from_string(const char*)`
6. 在 optcalc.h 中：同步更新声明，去掉 `#include <stdio.h>`（不再需要 FILE*）
7. 编译验证

### 阶段 2：核心功能移植
1. 写 `gropt_interface.c/h` — 实现 spot/transfer/psf/raytrace + strbuf 动态字符串
2. 写 `gropt_cdef.pxd` — Cython 声明
3. 写 `gropt.pxi` — Cython 实现
4. 更新 `setup.py`
5. 编译验证

### 阶段 3：Python API + 测试
1. 写 `gropt.py` — Optics 类 + from_string
2. 寻找或创建 .opt 测试数据
3. 用 CLI gropt 生成参考输出
4. 与 pyfitsh API 对比验证

### 未来升级（方案 B）
Python 层完全解析 .opt 文件，构建 Python 对象（Glass/Lens/Optics 类），Cython 层从 Python 对象填充 C 结构体。优点：
- Python 层可直接访问各透镜参数
- 支持动态修改光学系统（增删透镜、改参数）
- OpenSCAD/EPS 导出可纯 Python 实现（不需要 C 层 strbuf）
- 更好的 Python 交互体验（如 `opt.lenses[0].thickness = 5.0`）

本次不实现，留作日后升级。

---

## 七、工作量估计

| 内容 | 预计行数 | 难度 |
|------|----------|------|
| optcalc.c 改写 read_optics → from_string | ~60 | 低 |
| gropt_interface.c (spot/transfer/psf/raytrace) | ~250 | 中 |
| gropt_interface.c (scad/eps + strbuf) | ~250 | 中 |
| gropt_cdef.pxd | ~100 | 低 |
| gropt.pxi | ~150 | 中 |
| gropt.py | ~120 | 低 |
| setup.py 更新 | ~20 | 低 |
| 测试 | ~100 | 中 |
| **总计** | **~1050** | **中等** |
