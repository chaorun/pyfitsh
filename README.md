# pyfitsh

fitsh C 天文软件管线 的 Cython 移植版本。通过 Python 调用 `fistar`（星检测/拟合）、`fiphot`（孔径测光）、`grmatch`（星点匹配）、`fitrans`（图像变换）、`ficonv`（图像卷积/减影）、`firandom`（人工图像生成）。

## 已知问题

- **ficonv 图像相减未修复**：`ref_conv.fits`（卷积结果）匹配 CLI，但 `--output-subtracted`（减影结果）与 CLI 不一致（max_diff ~46400），mask 像素的零值处理逻辑有差异。
- **fistar parabolapeak 拟合**：bg/amp/flux 与 CLI 有 ~1% 差异（编译优化所致）。uplink (ALG_LNK) + gauss 模型已通过验证（ncand/nstar 与 CLI 一致）。`input_candidates` 二次引用候选时拟合结果异常（35 vs 49 待修）。

## fiphot 已验证功能

- **`photometry()`** — 标准多孔径测光，与 CLI 一致（含 `raw` 输出）
- **`photometry_from_raw()`** — 从 pre-computed raw 数据做 subtracted photometry，与 CLI `--input-raw-photometry` 对比：49 颗星 flux 差异 < 0.6%，bg 差异 < 0.84。差异来源为 CLI ASCII raw 文件精度截断而非实现错误
- **`subtracted_photometry()`** — 从 star 列表直接做 subtracted photometry
- **`magfit()`** — 星等校准

## firandom 已验证功能

- **`generate()`** — 人工图像生成，背景 + 星点 + 光子噪声，与 CLI 逐像素一致（262144/262144 像素完全匹配）
- 支持 firandom 表达式语法（`--list` 参数格式）
- 支持 Python list-of-dicts 格式输入星点
- PSN 表达式解析器已集成（`psn/` 目录）
- **已知限制**：PSN 常量 > 32767 会被 16-bit 截断（与 CLI 行为一致，原版同样存在）
- CLI 的 `--seed` 参数存在已知 bug：无条件覆盖 `--seed-noise`/`--seed-spatial`/`--seed-photon`（Python 端已正确实现为 fallback 逻辑）

## 构建

```bash
cd pyfitsh
python3 setup.py build_ext --inplace
```

依赖：Cython、numpy、astropy、scipy。

## 文件结构

```
pyfitsh/
├── setup.py          # 编译配置（-O3 优化）
├── fitsh_cy.pyx      # Cython 主模块
├── fitsh_cy.pxd      # C 函数声明
├── fitsh.h           # 通用宏/定义
├── common.h          # 公共函数声明
├── image.h, mask.h   # 轻量 FITS 类型替代
├── algorithms/       # 通用算法 (含 tokenize)
├── fistar/           # 星检测/拟合
├── fiphot/           # 孔径测光
├── ficonv/           # 卷积/减影
├── fitrans/          # 图像变换
├── firandom/         # 人工图像生成 (含 random, PSN)
├── psn/              # PSN 表达式解析器
├── grmatch/          # 星点匹配
├── grtrans/          # 坐标变换
├── math/             # 数学库
├── index/            # 索引/排序
├── link/             # 连通域
├── examples/         # 用法示例
├── utils.py          # MASKINFO 解码等工具
```

## 返回格式

Fistar 和 Fiphot 返回 `types.SimpleNamespace`，属性访问替代字典键：
- `r.output` — astropy Table（主结果表）
- `r.nstar` — 检测星数
- `r.dict` — 向后兼容字典

## 版本

2026-06-17 快照。编译选项 `-O3`，与 CLI 的 `Makefile -O3` 对齐。
- `io/` 目录已移除，tokenize 移至 `algorithms/`
- 死代码（CLI 文件 I/O）已 `//` 注释，禁止 Cython 读写文件
- `firandom` 模块集成完成，逐像素验证通过（262144/262144）
- PSN 表达式库集成（`psn/` 目录）
- 函数命名已清理：无下划线前缀（`_fiphot_` → `fiphot_`），`_flat` 后缀 → `_cy` 后缀
