"""公共对比工具：表格逐行逐列 diff，FITS 图像像素 diff"""
import sys, os, re, math
import numpy as np
from astropy.io import fits
from astropy.table import Table

RED = '\033[91m'
GREEN = '\033[92m'
NC = '\033[0m'

ROOT = os.path.dirname(os.path.abspath(__file__))


def _is_close(a, b, rtol=1e-5, atol=1e-8):
    if np.isnan(a) and np.isnan(b):
        return True
    if np.isnan(a) or np.isnan(b):
        return False
    return math.isclose(a, b, rel_tol=rtol, abs_tol=atol)


def compare_table(cy_table, cli_path, skip_cols=None, rtol=1e-5, atol=1e-8):
    """对比 Cython Table 和 CLI .cat 文件，逐行逐列"""
    cli_tab = Table.read(cli_path, format='ascii')

    if len(cy_table) != len(cli_tab):
        return False, f"row count: CY={len(cy_table)} CLI={len(cli_tab)}"

    skip = set(skip_cols or [])

    cy_cols = [c for c in cy_table.colnames if c not in skip]
    cli_cols = [c for c in cli_tab.colnames if c not in skip]

    if set(cy_cols) != set(cli_cols):
        only_cy = set(cy_cols) - set(cli_cols)
        only_cli = set(cli_cols) - set(cy_cols)
        msg = f"column mismatch: only CY={only_cy}, only CLI={only_cli}"
        return False, msg

    diffs = []
    for col in sorted(cy_cols):
        cy_arr = np.asarray(cy_table[col], dtype=np.float64)
        cli_arr = np.asarray(cli_tab[col], dtype=np.float64)

        if cy_arr.dtype.kind not in 'fc':
            if not np.array_equal(cy_arr, cli_arr):
                for i in range(len(cy_arr)):
                    if cy_arr[i] != cli_arr[i]:
                        diffs.append(f"  {col}[{i}]: CY={cy_arr[i]} CLI={cli_arr[i]}")
                        if len(diffs) >= 30:
                            break
            continue

        for i in range(len(cy_arr)):
            if not _is_close(cy_arr[i], cli_arr[i], rtol=rtol, atol=atol):
                diffs.append(f"  {col}[{i}]: CY={cy_arr[i]:.8g} CLI={cli_arr[i]:.8g}")
                if len(diffs) >= 30:
                    break
        if len(diffs) >= 30:
            break

    if diffs:
        return False, f"{len(diffs)} mismatch(es) first:\n" + "\n".join(diffs[:30])
    return True, "ok"


def compare_fits(cy_fits, cli_path, rtol=1e-5, atol=1e-8):
    """对比 FITS / numpy 数组与 CLI FITS 文件，计算像素 diff"""
    with fits.open(cli_path) as hdul:
        cli_data = hdul[0].data.astype(np.float64)

    if isinstance(cy_fits, np.ndarray):
        cy_data = cy_fits.astype(np.float64)
    elif hasattr(cy_fits, 'data'):
        cy_data = np.asarray(cy_fits.data, dtype=np.float64)
    else:
        cy_data = np.asarray(cy_fits, dtype=np.float64)

    if cy_data.shape != cli_data.shape:
        return False, f"shape: CY={cy_data.shape} CLI={cli_data.shape}"

    if np.allclose(cy_data, cli_data, rtol=rtol, atol=atol, equal_nan=True):
        return True, "ok"

    diff = np.abs(cy_data - cli_data)
    max_diff = np.nanmax(diff)
    mean_diff = np.nanmean(diff)
    n_diff = np.sum(~np.isclose(cy_data, cli_data, rtol=rtol, atol=atol, equal_nan=True))

    cy_nan = np.sum(np.isnan(cy_data))
    cli_nan = np.sum(np.isnan(cli_data))

    detail = f"max_diff={max_diff:.6g} mean_diff={mean_diff:.6g} n_diff={n_diff}/{cy_data.size}"
    if cy_nan != cli_nan:
        detail += f" nan: CY={cy_nan} CLI={cli_nan}"
    return False, detail


def compare_cat(cy_table, cli_path, **kw):
    """兼容名"""
    return compare_table(cy_table, cli_path, **kw)


def run_test(test_fn):
    name = os.path.basename(sys.argv[0]).replace('.py','')
    print(f"{'='*60}")
    print(f"  {name}")
    print(f"{'='*60}")
    try:
        ok, msg = test_fn()
    except Exception as e:
        import traceback
        traceback.print_exc()
        ok, msg = False, str(e)

    if ok:
        print(f"\n{GREEN}PASS{NC} {msg}")
        return 0
    else:
        print(f"\n{RED}FAIL{NC} {msg}")
        return 1
