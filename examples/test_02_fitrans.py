#!/usr/bin/env python3
"""test_02: fitrans 逆变换对比"""
import sys, os, numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '..')
from compare_utils import compare_table, compare_fits, compare_cat, run_test, RED, GREEN, NC
from astropy.io import fits

OBS_IMG = "L1_GRB260604C-Ic-180s-20260607_023443.fits"

def parse_trans_file(filepath):
    """解析 cli_test_grb.trans 文本格式为 transformation 字典"""
    trans = {'type': 'polynomial', 'order': 1, 'offset': [0.0, 0.0], 'scale': 1.0}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            key, _, value = line.partition('=')
            key = key.strip()
            value = value.strip()
            if key == 'type':
                trans['type'] = value
            elif key == 'order':
                trans['order'] = int(value)
            elif key == 'offset':
                parts = [float(x) for x in value.split(',')]
                trans['offset'] = parts
            elif key == 'scale':
                trans['scale'] = float(value)
            elif key == 'dxfit':
                trans['dxfit'] = [float(x) for x in value.split(',')]
            elif key == 'dyfit':
                trans['dyfit'] = [float(x) for x in value.split(',')]
    return trans

def do_test():
    from pyfitsh import Fitrans

    img = fits.getdata(OBS_IMG).astype(np.float64)
    trans = parse_trans_file("cli_test_grb.trans")

    ft = Fitrans.ImageTransformation(transformation=trans, inverse=True,
                                     method='bilinear')
    output, mask = ft.apply(img)

    ok, msg = compare_fits(output, "cli_test_v1.fits")
    if not ok:
        return False, msg
    return True, "ok"

if __name__ == "__main__":
    sys.exit(run_test(do_test))
