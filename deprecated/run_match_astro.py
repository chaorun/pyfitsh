#!/usr/bin/env python3
"""
run_match_astro.py -- grmatch + astropy 匹配示例。

依赖: astropy numpy grmatch_cy
构建: python3 setup.py build_ext --inplace
用法: python3 run_match_astro.py
"""
import time
import numpy as np
from astropy.table import Table, hstack
from grmatch_cy import Matcher
# 用 astropy 读入参考星表和目标星表
ref = Table.read("testdata/temp_ref.txt", format="ascii.fast_no_header")
obj = Table.read("testdata/temp_obs.txt", format="ascii.fast_no_header")
ref.rename_columns(ref.colnames, ["id","x","y","ord"])
obj.rename_columns(obj.colnames, ["id","x","y","ord"])
nref, ninp = len(ref), len(obj)
print(f"ref: {nref} rows | obj: {ninp} rows")
# 提取坐标和排序值 (col2=x, col3=y, col4=ord 取负=升序)
ref_x = ref["x"].astype(np.float64).tolist()
ref_y = ref["y"].astype(np.float64).tolist()
obj_x = obj["x"].astype(np.float64).tolist()
obj_y = obj["y"].astype(np.float64).tolist()
ref_ord = (-ref["ord"].astype(np.float64)).tolist()
obj_ord = (-obj["ord"].astype(np.float64)).tolist()
# 运行点匹配 (order=2 maxdist=2 unitarity=0.01)
m = Matcher(order=2, maxdist=2.0, unitarity=0.01)
m.set_reference(ref_x, ref_y, ref_ord)
m.set_input(obj_x, obj_y, obj_ord)
t0 = time.time()
r = m.run()
dt = time.time() - t0
ns = r.stats['nsigma']
un = r.stats['unitarity']
hc = r.stats['hull_coverage']
print(f"{r.nhit} pairs matched in {dt:.3f}s")
print(f"  nsigma={ns:.6f}  unitarity={un:.8f}  hull={hc:.4f}")
# 构建匹配表: hstack([ref[idx_ref], obj[idx_inp]])
matched = hstack([ref[np.array(r.hits_ref)], obj[np.array(r.hits_inp)]])
print(f"\n匹配表 (前 5 行):")
print(matched[:5])
# 输出变换系数 (dxfit dyfit 多项式系数)
dx_vfits = r.vfits_dx
dy_vfits = r.vfits_dy
dx_str = ' + '.join(f'{v:.10g}*T{k}' for k, v in enumerate(dx_vfits))
dy_str = ' + '.join(f'{v:.10g}*T{k}' for k, v in enumerate(dy_vfits))
print(f"\n变换系数 (order={r.order}, nvar={r.nvar}):")
print(f"  dx = {dx_str}")
print(f"  dy = {dy_str}")
def eval_2d_poly(x, y, coeff, order):
    """计算 2D 多项式, 与 C 版本 eval_2d_poly 等价。"""
    if order == 0:
        return coeff[0]
    elif order == 1:
        return coeff[0] + x * coeff[1] + y * coeff[2]
    elif order == 2:
        return (coeff[0] + x * coeff[1] + y * coeff[2]
                + x * x * coeff[3] / 2.0  # x^2/2! 项
                + x * y * coeff[4]        # x*y 项
                + y * y * coeff[5] / 2.0) # y^2/2! 项
    # 通用金字塔展开 (order>=3)
    idx = 0
    result = 0.0
    for k in range(order + 1):
        for j in range(k + 1):
            term = 1.0
            for p in range(1, k - j + 1):
                term *= x / p  # x^i / i!
            for p in range(1, j + 1):
                term *= y / p  # y^j / j!
            result += coeff[idx] * term
            idx += 1
    return result
# 验证变换: 对匹配的 ref 坐标应用多项式, 与 obj 实际坐标对比
print("\n变换验证 (前 3 对):")
for i in range(3):
    ri = r.hits_ref[i]     # ref 表中匹配行的索引
    oi = r.hits_inp[i]     # obj 表中匹配行的索引
    px = eval_2d_poly(ref_x[ri], ref_y[ri], dx_vfits, r.order)
    py = eval_2d_poly(ref_x[ri], ref_y[ri], dy_vfits, r.order)
    dx = px - obj_x[oi]    # x 方向残差
    dy = py - obj_y[oi]    # y 方向残差
    err = np.sqrt(dx * dx + dy * dy)  # 二维距离
    print(f"  ref[{ri}] ({ref_x[ri]:.2f},{ref_y[ri]:.2f})"
          f" -> pred ({px:.2f},{py:.2f})"
          f" vs obj[{oi}] ({obj_x[oi]:.2f},{obj_y[oi]:.2f})"
          f"  err={err:.4f}")
