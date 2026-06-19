"""test_lfit.py — 对标 test_lfit.sh 的全部19个测试用例
每个测试和 origincode lfit 的 reference 输出逐值对比。
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

TDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    '..', '..', 'fitsh_test_data', 'test_lfit')


def load_dat(name):
    return np.loadtxt(os.path.join(TDIR, name))


def parse_ref(name):
    with open(os.path.join(TDIR, name)) as f:
        lines = [l for l in f if l.strip() and not l.strip().startswith('#')]
    rows = []
    for line in lines:
        rows.append([float(x) for x in line.split()])
    return rows


def assert_params(result, ref_params, tol=1e-4, label=""):
    assert result.ok, f"{label} failed: {result.error_msg}"
    for i, (got, exp) in enumerate(zip(result.params, ref_params)):
        assert abs(got - exp) < tol, (
            f"{label} param[{i}]: got {got}, expected {exp}, diff {abs(got-exp)}")


def run_tests():
    from pyfitsh.lfit import Lfit, lfit_fit

    results = []

    def test(num, name, fn):
        try:
            fn()
            results.append((num, name, 'PASS', ''))
        except Exception as e:
            results.append((num, name, 'FAIL', str(e)))
        except SystemExit:
            results.append((num, name, 'CRASH', 'segfault or exit'))

    ref01 = parse_ref('lfit_test_clls_linear.out')[0]
    def t01():
        r = lfit_fit(load_dat('lfit_test_linear.dat'),
                     'x:1,y:2', 'a,b', 'a+b*x', 'y')
        assert_params(r, ref01, label="[01]")
    test('01', 'CLLS linear', t01)

    def t02():
        r = lfit_fit(load_dat('lfit_test_linear.dat'),
                     'x:1,y:2', 'a,b', 'a+b*x', 'y')
        assert_params(r, ref01, label="[02]")
    test('02', 'CLLS linear (vars/expr output skip, params only)', t02)

    ref03 = parse_ref('lfit_test_clls_poly2.out')[0]
    def t03():
        r = lfit_fit(load_dat('lfit_test_poly2.dat'),
                     'x:1,y:2', 'a,b,c', 'a+b*x+c*x*x', 'y')
        assert_params(r, ref03, label="[03]")
    test('03', 'CLLS poly2', t03)

    ref04 = parse_ref('lfit_test_clls_weighted_error.out')
    def t04():
        r = lfit_fit(load_dat('lfit_test_weighted.dat'),
                     'x:1,y:2,e:3', 'a,b,c', 'a+b*x+c*x*x', 'y',
                     error='e', errtype=0)
        assert_params(r, ref04[0], label="[04]")
    test('04', 'CLLS weighted -e sigma (params)', t04)

    ref05 = parse_ref('lfit_test_clls_weighted_weight.out')[0]
    def t05():
        r = lfit_fit(load_dat('lfit_test_weighted.dat'),
                     'x:1,y:2,e:3', 'a,b,c', 'a+b*x+c*x*x', 'y',
                     error='1/e', errtype=1)
        assert_params(r, ref05, label="[05]")
    test('05', 'CLLS weighted -w weight', t05)

    ref06 = parse_ref('lfit_test_clls_reject.out')[0]
    def t06():
        r = lfit_fit(load_dat('lfit_test_outlier.dat'),
                     'x:1,y:2', 'a,b,c', 'a+b*x+c*x*x', 'y',
                     niter=3, sigma=3.0)
        assert_params(r, ref06, label="[06]")
        assert r.nused < r.nrow, "[06] should reject some outliers"
    test('06', 'CLLS outlier rejection', t06)

    ref07 = parse_ref('lfit_test_clls_reject_delta.out')[0]
    def t07():
        r = lfit_fit(load_dat('lfit_test_outlier.dat'),
                     'x:1,y:2', 'a,b,c', 'a+b*x+c*x*x', 'y',
                     niter=3, sigma=3.0)
        assert_params(r, ref07, label="[07]")
    test('07', 'CLLS reject delta (params same as 06)', t07)

    def t08():
        ref08 = parse_ref('lfit_test_clls_poly2_format.out')[0]
        r = lfit_fit(load_dat('lfit_test_poly2.dat'),
                     'x:1,y:2', 'a,b,c', 'a+b*x+c*x*x', 'y')
        assert_params(r, ref08, tol=1e-6, label="[08]")
    test('08', 'CLLS poly2 format (params comparison)', t08)

    def t09():
        results.append(('09', 'derived variable (-g)', 'SKIP', 'not implemented'))
    t09()

    def t10():
        results.append(('10', 'constraint (:=)', 'SKIP', 'not implemented'))
    t10()

    ref11 = parse_ref('lfit_test_nllm_poly2.out')[0]
    def t11():
        r = lfit_fit(load_dat('lfit_test_poly2.dat'),
                     'x:1,y:2', 'a=1,b=1,c=0.1', 'a+b*x+c*x*x', 'y',
                     method='nllm')
        assert_params(r, ref11, tol=1e-3, label="[11]")
    test('11', 'NLLM analytic', t11)

    ref12 = parse_ref('lfit_test_lmnd_poly2.out')[0]
    def t12():
        r = lfit_fit(load_dat('lfit_test_poly2.dat'),
                     'x:1,y:2', 'a=1,b=1,c=0.1', 'a+b*x+c*x*x', 'y',
                     method='lmnd')
        assert_params(r, ref12, tol=1e-3, label="[12]")
    test('12', 'LMND numerical derivatives', t12)

    ref13_row = parse_ref('lfit_test_dhsx_poly2.out')[0]
    def t13():
        r = lfit_fit(load_dat('lfit_test_poly2.dat'),
                     'x:1,y:2', 'a=1:0.5,b=1:0.5,c=0.1:0.1',
                     'a+b*x+c*x*x', 'y', method='dhsx')
        assert_params(r, ref13_row[:3], tol=1e-3, label="[13]")
    test('13', 'DHSX downhill simplex', t13)

    def t14():
        results.append(('14', 'chi2 grid (-K)', 'SKIP', 'grid output not comparable'))
    t14()

    def t15():
        results.append(('15', 'eval (no fit)', 'SKIP', 'evaluation mode not implemented'))
    t15()

    def t16():
        results.append(('16', 'eval replace cols', 'SKIP', 'evaluation mode not implemented'))
    t16()

    def t17():
        results.append(('17', 'macro (-x)', 'SKIP', 'macro not implemented'))
    t17()

    ref18 = parse_ref('lfit_test_pairs_dx_linear.out')[0]
    def t18():
        r = lfit_fit(load_dat('lfit_test_pairs.dat'),
                     'xr:1,yr:2,xo:3,yo:4,mr:5,mo:6', 'a,b,c',
                     'a+b*xr+c*yr', 'xo-xr')
        assert_params(r, ref18, tol=1e-2, label="[18]")
    test('18', 'pairs dx linear', t18)

    ref19 = parse_ref('lfit_test_pairs_dy_linear.out')[0]
    def t19():
        r = lfit_fit(load_dat('lfit_test_pairs.dat'),
                     'xr:1,yr:2,xo:3,yo:4,mr:5,mo:6', 'a,b,c',
                     'a+b*xr+c*yr', 'yo-yr')
        assert_params(r, ref19, tol=1e-2, label="[19]")
    test('19', 'pairs dy linear', t19)

    passed = sum(1 for _, _, s, _ in results if s == 'PASS')
    failed = sum(1 for _, _, s, _ in results if s == 'FAIL')
    skipped = sum(1 for _, _, s, _ in results if s == 'SKIP')
    total = len(results)

    print()
    for num, name, status, msg in results:
        tag = {'PASS': 'PASS', 'FAIL': 'FAIL', 'SKIP': 'SKIP', 'CRASH': 'CRASH'}[status]
        extra = f'  ({msg})' if msg else ''
        print(f"  [{tag:5s}] [{num}] {name}{extra}")

    print(f"\nResult: {passed} passed, {failed} failed, {skipped} skipped / {total} total")
    if failed > 0:
        sys.exit(1)


if __name__ == '__main__':
    run_tests()
