# pyfitsh

fitsh C 天文软件管线的 Cython 移植版本。通过 Python 调用 `fiarith`（逐像素算术表达式）、`fistar`（星检测/拟合）、`fiphot`（孔径测光）、`grmatch`（星点匹配）、`fitrans`（图像变换）、`ficonv`（图像卷积/减影）、`firandom`（人工图像生成）。

## 构建

```bash
cd pyfitsh
python3 setup.py build_ext --inplace
```

依赖：Cython、numpy、astropy、scipy。

## 文件结构

```
pyfitsh/
├── setup.py          # 编译配置
├── core.pyx          # Cython 主模块
├── core.pxd          # C 函数声明
├── fitsh.h           # 通用宏/定义
├── common.h          # 公共函数声明
├── fiarith/          # 逐像素算术表达式求值器
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
├── deprecated/       # 废弃文件
└── utils.py          # MASKINFO 解码等工具
```

## 返回格式

Fistar 和 Fiphot 返回 `types.SimpleNamespace`，属性访问替代字典键：
- `r.output` — astropy Table（主结果表）
- `r.nstar` — 检测星数
- `r.dict` — 向后兼容字典

## 版本

2026-06-17 快照。复刻 fitsh 0.9.4 行为，遵循原始设计约束：
- UDF 子表达式 `[...](body)` 不支持图像函数（sq, smooth, laplace 等），仅支持 `psn_general_fn` 中的通用数学函数（sin, cos, exp, ln, sqrt, abs, sign 等），与原始 fitsh CLI 行为一致
- evaluate.c 中 `"if_sq"` 修正为 `"sq"`
- 死代码（CLI 文件 I/O）已注释，禁止 Cython 读写文件
- 函数命名已清理：无下划线前缀，`_flat` 后缀 → `_cy` 后缀
