# read_raw_photometry: CLI 原版 vs Cython 版本 — 逐字段对比

## 对比对象

| | CLI 原版 | Cython 版本 |
|---|---|---|
| **文件** | `tmp/origincode/src/fiphot-io.c` | `fitsh_cython/fiphot/fiphot_core.c` |
| **函数** | `read_raw_photometry()` | `read_raw_photometry_cy()` |
| **行号** | 275-374 | 1735-1805 |
| **输入** | 文本流 (FILE*)，逐行 ASCII | numpy float64 数组 (nstar × nap × 12) |

---

## 逐字段对比

### 1. x

| | CLI | Cython |
|---|---|---|
| 代码 | `sscanf(cmd[1],"%lg",&x)` → `ps[np].x=x` | `p->x = raw_data[i*nap*12 + 0*12 + 1]` |
| 数据来源 | 文本行第 2 列 | numpy 列 1 |
| 差异 | 相同 | |

### 2. y

| | CLI | Cython |
|---|---|---|
| 代码 | `sscanf(cmd[2],"%lg",&y)` → `ps[np].y=y` | `p->y = raw_data[i*nap*12 + 0*12 + 2]` |
| 数据来源 | 文本行第 3 列 | numpy 列 2 |
| 差异 | 相同 | |

### 3. n / ninap

| | CLI | Cython |
|---|---|---|
| 代码 | `sscanf(cmd[3],"%d",&n)` → `ps[np].n = ps[np].ninap = n` | `p->n = nap; p->ninap = nap` |
| 数据来源 | 文本行第 4 列（每行可变） | numpy shape[1]（全矩阵统一） |
| **差异** | **每 star 的 n 独立** | **所有 star 统一 nap** |

### 4. use_ref

| | CLI | Cython |
|---|---|---|
| 代码 | `if(sscanf(cmd[4],"%lg",&ref_mag)==1) ps[np].use_ref=1; else ps[np].use_ref=0;` | 由 `ref_mag_arr[i] != 0.0` 决定 |
| 数据来源 | 文本第 5 列（"-" 表示无 ref） | 外部传入数组 |
| **差异** | **内建解析** | **Python 侧传入** |

### 5. ref_mag

| | CLI | Cython |
|---|---|---|
| 代码 | `ps[np].ref_mag = ref_mag;` (有值) / `ps[np].ref_mag = 0.0;` (无值) | `p->ref_mag = ref_mag_arr[i]` 或 `0.0` |
| **差异** | **内建解析** | **外部传入** |

### 6. ref_col

| | CLI | Cython |
|---|---|---|
| 代码 | `sscanf(cmd[5],"%lg",&ref_col)` → `ps[np].ref_col = ref_col;` 或 `= 0.0` | `p->ref_col = ref_col_arr ? ref_col_arr[i] : 0.0` |
| **差异** | **内建解析** | **外部传入** |

### 7. ref_err

| | CLI | Cython |
|---|---|---|
| 代码 | **未设置**（realloc 不含初始化，垃圾数据） | `p->ref_err = ref_err_arr ? ref_err_arr[i] : 0.0` |
| 差异 | CLI 为未定义行为 | Cython 显式零初始化 |

### 8. inaps

| | CLI | Cython |
|---|---|---|
| 代码 | `inaps = malloc(sizeof(apgeom)*n); for(i..){ aperture_definition_parse(apdef,&ap); inaps[i]=ap; } ps[np].inaps=inaps;` | `p->inaps = malloc(sizeof(apgeom)*nap); for(j..){ memset(&p->inaps[j],0,sizeof(apgeom)); p->inaps[j].apgeom_type=APGEOM_TYPE_CIRCULAR; p->inaps[j].r0=... p->inaps[j].ra=... p->inaps[j].da=... }` |
| 所有权 | **每 star 独立 malloc** | **每 star 独立 malloc**（修复后一致） |
| 构建方式 | `aperture_definition_parse` (sscanf "r0:ra:da") | 直接赋值 r0/ra/da |
| **差异** | **支持 CIRCULAR / POLYGON 两种类型** | **仅 CIRCULAR** |

### 9. rfflux

| | CLI | Cython |
|---|---|---|
| 代码 | `rfflux = malloc(sizeof(photflux)*n); for(i..){ sscanf flux,fluxerr,flag; 其余字段=0; } ps[np].rfflux=rfflux;` | `p->rfflux = malloc(sizeof(photflux)*nap); for(j..){ flux/err/flag 从 numpy 列 10/11/4; 其余字段=0; }` |
| 所有权 | **每 star 独立 malloc** | **每 star 独立 malloc** |
| **差异** | **flag 通过 sscanf 解析 hex 并检查 r<3 设置 MASK_BAD** | **flag 直接 (int) 截断，无条件赋值** |

### 10. fluxes

| | CLI | Cython |
|---|---|---|
| 代码 | `ps[np].fluxes = NULL` | `p->fluxes = NULL` |
| 差异 | 相同 | |

### 11. id

| | CLI | Cython |
|---|---|---|
| 代码 | `ps[np].id = strdup(id)` (或 NULL) | `p->id = NULL` |
| **差异** | **保留星标识符字符串** | **始终为 NULL** |

---

## 关键差异汇总

| # | 字段 | 严重程度 | 说明 |
|---|------|---------|------|
| 1 | `use_ref/ref_mag/ref_col` | **已解决** | Python 侧传入 ref 数组，C 层动态设置 |
| 2 | `id` | 低 | CLI 保留星标识符用于输出；Cython 丢失。不影响测光 |
| 3 | `n/ninap` | 中 | CLI 每星独立 n；Cython 统一 nap。当前用例（nap=1）无影响 |
| 4 | `rfflux.flag` 解析 | 中 | CLI 对 hex 解析失败设置 MASK_BAD + 清零；Cython 无保护 |
| 5 | `inaps` 构建 | 低 | CLI 支持 POLYGON；Cython 仅 CIRCULAR |
| 6 | `inaps` 内存所有权 | **已修复** | 每星独立，无 double-free |

---

## 验证结果（20260616T104000Z）

CLI `--input-raw-photometry` vs Cython `photometry_from_raw` 对比：
- 49 颗星，flux 差异 < 0.6%（最大 0.59%），bg 差异 < 0.84
- 差异来源：CLI ASCII raw 文件精度截断
