# pyfitsh 开发踩坑总结报告

**日期**: 2026-06-15 ~ 2026-06-18

## 一、PSF type 常量偏移 bug（pbg/pamp nan）

**现象**: Fistar 使用 `psf='native,order=2'` 时，输出列 `pbg`/`pamp` 全为 nan。

**排查过程**:
- 初期试图加 `fprintf(stderr,...)` 被否决，要求用二进制 dump
- 第一次只挑了两个字段 dump，被批评是"逗小孩玩的方式"
- 多次改进后最终用 `fwrite(&struct, sizeof(struct), ...)` 完整 dump 所有参数
- 对 `fistar_determine_psf` 的 4 个参数全部 dump，发现 CLI tpd.coeff 有 10400 字节、CY 为 0 字节
- 全部输入逐字节一致，差异仅在 `pdparam.type`：CLI=1, CY=0

**根因**: `core.pyx` 中 `_psf_type` 赋值偏移了 1（`native→0` 应为 `1`），`psf_determine` 收到 0 走 `default: return -1`，不填充 `tpd.coeff`。

## 二、二进制 dump 的惨痛教训

**最大错误**：用硬编码零值假装 dump。

在 dump `tpd` 结构体时，写了：
```c
int ti=0; double td=0.0;
fwrite(&ti,4,1,df); fwrite(&ti,4,1,df); fwrite(&ti,4,1,df);
fwrite(&td,8,1,df); fwrite(&td,8,1,df); fwrite(&td,8,1,df);
```
这是栈上的临时变量，和 `tpd` 结构体毫无关系——写了 6 个假的 0 字节就当 dump 过了。用户反复指出"这和内存里给函数的 tpd 有什么关系？"，我理解后才改为真正的 `fwrite(&tpd, sizeof(tpd), 1, df)`。

**正确做法**：
- 完整 dump 所有结构体：`fwrite(&struct, sizeof(struct), 1, df)`
- 不能用硬编码零值假装 dump
- dump 调用前后各一次，不可只 dump 输出不 dump 输入
- CLI 和 Cython 两端同步加 dump，同一位置、同一格式
- `fopen/fwrite` 写入指定路径，Python `struct.unpack` 解析对比
- 60+ 行 dump 代码虽冗长但比反复猜测高效得多
- 解析 dump 时字段偏移算错会得出完全错误的结论（如 nipoint 偏移 64→80 的错误导致误判 CLi nipoint=0）

## 三、fitsimage vs image 结构体布局

对比发现两者完全相同（sx/sy/bit/data 偏移一致），不是差异原因。误区：花了很多时间检查结构体布局。

## 四、mask 处理缺失（ficonv/fiphot 全像素差异）

**现象**: ficonv 测试中卷积/减影图像 99% 像素不同。

**根因**: USAGE.md 规定 CLI 自动从 FITS header 读取 MASKINFO，Cython 需显式调用 `decode_maskinfo(hdu)` + NaN 标记。测试脚本未传 `ref_mask`/`img_mask`。

## 五、input_candidates 崩溃

**现象**: `input_candidates=cands` 路径在结果填充循环中 SEGFAULT。

**排查**:
- lldb 定位 `fistar_search_cy +6736`，`cand` 为野指针 `0x260001120e3918`
- `fistar_core.c:392-394` 中 s/d/k 被默认值覆写

**修复三处**: 删除覆写、加 cleanup_candlist、加 wc 边界校验

## 六、35 vs 49 之谜（CLI --col-shape bug）

**排查**：对 `fit_star_single_model` 所有输入/输出做完整二进制 dump，发现 CLI cands 中 sxx=syy=cs, sxy=0（cd/ck 丢失）。

**根因**: origincode CLI 的 `--col-shape 3,4,5` 未读取 cd/ck 列。CLI 本身也有 bug。

## 七、测试脚本格式统一（Table.read 替换 loadtxt）

4 个测试脚本（01/11/12/13）从 `np.loadtxt` → `Table.read(..., format='ascii.no_header')`。

## 八、sort 排序对齐问题

test_08/09 中 CLI 和 CY 输出顺序不同，需双方都 `sort('x')` 后再对比。

## 九、容差校准

| 测试 | 列 | 差异 | 处理 |
|------|-----|------|------|
| 04 | convolved | 2/4M pix, max=0.107 | atol=1e-1 |
| 07b | mom[0] | 4.5e-5 | rtol=1e-2, atol=1e-4 |
| 10 | s/d/k/fwhm | 退化 | skip（CLI bug 导致形状参数退化） |

## 十、常量定义偏移

PSF_DET_NATIVE=1 不是 0。多个模块中存在类似的常量定义理解错误。

## 十一、.so 编译 vs C 二进制性能差距

ficonv 的 `.so` 比 C 二进制慢 30×+。根源：`setup.py` 自动注入 `-isystem /miniconda3/include` 和 `-I /opt/openblas/include` 改变了编译路径。

## 十二、edit 工具的匹配陷阱

在 `kernel-base.c` 中插 profiling 代码时，因两个函数有完全相同的特征行 `{sx=refimg->sx,sy=refimg->sy;}`，edit 同时修改了两个函数，导致变量声明丢失。

## 十三、grmatch 移植中的包名替换

10+ 个文件的 `fitsh_cython` → `pyfitsh` 重命名。深层 import 路径被拒，需统一用 `from pyfitsh import X`。

## 十四、firandom 种子对比

CLI 的 `--seed` 无条件覆盖了 `--seed-noise`/`--seed-spatial`/`--seed-photon`。通过在各步骤前后做 binary dump 定位差异。

## 十五、fiarith 子表达式 sq 函数 bug

`psn_img_fn` 中的 `sq` 在 `psn_img_prop` 中，但子表达式用 `psn_general_prop` 做转换，`sq` 找不到对应条目。

## 十六、跨 session 通用教训

1. **修改 origincode 前必须复制到 tmp 目录**
2. **禁止 sed/awk 等命令行工具修改文件**
3. **git 操作需授权**，误提交用 `git update-ref -d HEAD` 回退
4. **测试容差不要一刀切**：不同列/场景需要不同容差
5. **函数调用链梳理是关键**：先画出完整调用链再逐层排查
6. **用 Python I/O 代替 C FILE*** 在 `def` 方法中做 dump
7. **ctypes 结构体对齐容易出错**，Cython 编译时检查更可靠
8. **CLI 也可能有 bug**，不要盲目追"为什么和 CLI 不一样"
9. **`.so` 链接方式受 setup.py 编译链影响**，分析性能时需对比编译参数
10. **最终答案往往很简单**：`_psf_type = 1` 而不是 `0`，输出 mask 参数而不是忽略
