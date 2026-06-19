#!/usr/bin/env python3
"""run_all_lfit_tests.py — 对标3个shell脚本的全部47个测试"""
import subprocess, sys, os, time

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
os.environ['PYTHONPATH'] = os.getcwd()

TDIR = 'fitsh_test_data/test_lfit'

def run_one(num, name, code, timeout=30):
    t0 = time.time()
    try:
        p = subprocess.run(
            [sys.executable, '-c', code],
            capture_output=True, text=True, timeout=timeout)
        dt = time.time() - t0
        out = p.stdout.strip()
        err = p.stderr.strip()
        if p.returncode == 0 and out:
            return (num, name, 'PASS', out, dt)
        elif p.returncode == -11:
            return (num, name, 'CRASH', 'SIGSEGV', dt)
        else:
            msg = err[-200:] if err else out[-200:] if out else f'rc={p.returncode}'
            return (num, name, 'FAIL', msg, dt)
    except subprocess.TimeoutExpired:
        return (num, name, 'TIMEOUT', f'>{timeout}s', time.time()-t0)

def mk(num, name, call, ref_params=None, tol=1e-3, timeout=30):
    check = ""
    if ref_params is not None:
        ref = repr(ref_params)
        check = f"""
ref={ref}
for i,(g,e) in enumerate(zip(r.params,ref)):
    assert abs(g-e)<{tol}, f'param[{{i}}]: got {{g}}, expected {{e}}'
print(f'params={{[round(p,6) for p in r.params]}} OK')
"""
    else:
        check = "print(f'rc={r.return_code} ok={r.ok} nrow={r.nrow}')"

    code = f"""
import numpy as np, sys; sys.path.insert(0,'.')
from pyfitsh.lfit import lfit_fit
TDIR='{TDIR}'
def ld(n): return np.loadtxt(f'{{TDIR}}/{{n}}')
{call}
assert r.ok, f'FAILED: {{r.error_msg}}'
{check}
"""
    return run_one(num, name, code, timeout)


results = []
print("=" * 60)
print("test_lfit.sh (19 tests)")
print("=" * 60)

results.append(mk('T01','CLLS linear',
    "r=lfit_fit(ld('lfit_test_linear.dat'),'a,b','x:1,y:2','a+b*x','y')",
    [1.99992, 3]))

results.append(mk('T02','CLLS linear (same as 01)',
    "r=lfit_fit(ld('lfit_test_linear.dat'),'a,b','x:1,y:2','a+b*x','y')",
    [1.99992, 3]))

results.append(mk('T03','CLLS poly2',
    "r=lfit_fit(ld('lfit_test_poly2.dat'),'a,b,c','x:1,y:2','a+b*x+c*x*x','y')",
    [1, 2, 0.499998]))

results.append(mk('T04','CLLS weighted -e',
    "r=lfit_fit(ld('lfit_test_weighted.dat'),'a,b,c','x:1,y:2,e:3','a+b*x+c*x*x','y',error='e')",
    [1, 2, 0.5]))

results.append(mk('T05','CLLS weighted -w',
    "r=lfit_fit(ld('lfit_test_weighted.dat'),'a,b,c','x:1,y:2,e:3','a+b*x+c*x*x','y',weight='1/e')",
    [1, 2, 0.5]))

results.append(mk('T06','CLLS reject',
    "r=lfit_fit(ld('lfit_test_outlier.dat'),'a,b,c','x:1,y:2','a+b*x+c*x*x','y',rejection_niter=3,rejection_level=3.0)",
    [1, 2, 0.499998]))

results.append(mk('T07','CLLS reject delta',
    "r=lfit_fit(ld('lfit_test_outlier.dat'),'a,b,c','x:1,y:2','a+b*x+c*x*x','y',rejection_niter=3,rejection_level=3.0,is_dump_delta=1)",
    [1, 2, 0.499998]))

results.append(mk('T08','CLLS poly2 format',
    "r=lfit_fit(ld('lfit_test_poly2.dat'),'a,b,c','x:1,y:2','a+b*x+c*x*x','y')",
    [1, 2, 0.5], tol=1e-5))

results.append(mk('T09','derived variable',
    "r=lfit_fit(ld('lfit_test_poly2.dat'),'a,b,c','x:1,y:2','a+b*x+c*x*x','y',derived_variables='sum=a+b+c')",
    [1, 2, 0.499998]))

results.append(mk('T10','constraint c:=0.5',
    "r=lfit_fit(ld('lfit_test_poly2.dat'),'a,b,c:=0.5','x:1,y:2','a+b*x+c*x*x','y')",
    [1, 2, 0.5]))

results.append(mk('T11','NLLM analytic',
    "r=lfit_fit(ld('lfit_test_poly2.dat'),'a=1,b=1,c=0.1','x:1,y:2','a+b*x+c*x*x','y',method='nllm',parameters='iterations=50')",
    [1, 2, 0.499998]))

results.append(mk('T12','LMND numerical',
    "r=lfit_fit(ld('lfit_test_poly2.dat'),'a=1,b=1,c=0.1','x:1,y:2','a+b*x+c*x*x','y',method='lmnd',differences='a=0.001,b=0.001,c=0.001',parameters='iterations=50')",
    [1, 2, 0.499998]))

results.append(mk('T13','DHSX simplex',
    "r=lfit_fit(ld('lfit_test_poly2.dat'),'a=1:0.5,b=1:0.5,c=0.1:0.1','x:1,y:2','a+b*x+c*x*x','y',method='dhsx')",
    [1, 2, 0.499998]))

results.append(mk('T14','chi2 grid',
    "r=lfit_fit(ld('lfit_test_poly2.dat'),'a=[0:1:3],b=[1:1:4],c=0.5','x:1,y:2','a+b*x+c*x*x','y',method='mchi')",
    None))

results.append(mk('T15','eval linear',
    "r=lfit_fit(ld('lfit_test_linear.dat'),'a,b','x:1,y:2','2+3*x')",
    None))

results.append(mk('T16','eval replace cols',
    "r=lfit_fit(ld('lfit_test_linear.dat'),'a,b','x:1,y:2','2+3*x,y-(2+3*x)',columns_output='1,2')",
    None))

results.append(mk('T17','macro',
    "r=lfit_fit(ld('lfit_test_linear.dat'),'a,b','x:1,y:2','line(a,b,x)','y',macros=['line(u,v,t)=u+v*t'])",
    [1.99992, 3]))

results.append(mk('T18','pairs dx',
    "r=lfit_fit(ld('lfit_test_pairs.dat'),'a,b,c','xr:1,yr:2,xo:3,yo:4,mr:5,mo:6','a+b*xr+c*yr','xo-xr')",
    [1132.72, -1.03116, -0.122157], tol=0.01))

results.append(mk('T19','pairs dy',
    "r=lfit_fit(ld('lfit_test_pairs.dat'),'a,b,c','xr:1,yr:2,xo:3,yo:4,mr:5,mo:6','a+b*xr+c*yr','yo-yr')",
    [305.122, -0.0547023, -0.282015], tol=0.01))

for num, name, status, msg, dt in results:
    print(f"  [{status:7s}] [{num}] {name} ({dt:.1f}s) {msg[:80] if status!='PASS' else ''}")

print()
print("=" * 60)
print("test_lfit_montecarlo.sh (19 tests)")
print("=" * 60)

mc_tests = [
    ('M01','MCMC basic',        "method='mcmc',parameters=None",                      50),
    ('M02','MCMC accepted',     "method='mcmc',parameters='accepted'",                 50),
    ('M03','MCMC nonaccepted',  "method='mcmc',parameters='nonaccepted'",              50),
    ('M04','MCMC gibbs',        "method='mcmc',parameters='gibbs'",                    50),
    ('M05','XMMC basic',        "method='xmmc'",                                       50),
    ('M06','XMMC skip',         "method='xmmc',parameters='skip'",                     50),
    ('M07','XMMC adaptive',     "method='xmmc',parameters='adaptive'",                 50),
    ('M08','XMMC window',       "method='xmmc',parameters='window=10'",                50),
]
for num, name, extra, iters in mc_tests:
    call = f"r=lfit_fit(ld('lfit_test_poly2.dat'),'a=1:0.2,b=1:0.2,c=0.1:0.1','x:1,y:2','a+b*x+c*x*x','y',{extra},seed=12345,mc_iterations={iters})"
    results.append(mk(num, name, call, None, timeout=60))

emce_tests = [
    ('M09','EMCE clls',  "'a,b,c'",       "parameters='clls'",  None, 20),
    ('M10','EMCE nllm',  "'a=1,b=1,c=0.1'","parameters='nllm'", None, 20),
    ('M11','EMCE lmnd',  "'a=1,b=1,c=0.1'","parameters='lmnd',differences='a=0.001,b=0.001,c=0.001'", None, 20),
    ('M12','EMCE dhsx',  "'a=1:0.5,b=1:0.5,c=0.1:0.1'","parameters='dhsx'", None, 20),
    ('M13','EMCE perturb',"'a,b,c'",       "parameters='clls',perturbations='0.01'", None, 20),
]
for num, name, vstr, extra, ref, iters in emce_tests:
    call = f"r=lfit_fit(ld('lfit_test_poly2.dat'),{vstr},'x:1,y:2','a+b*x+c*x*x','y',method='emce',{extra},seed=12345,mc_iterations={iters})"
    results.append(mk(num, name, call, None, timeout=60))

results.append(mk('M14','FIMA basic',
    "r=lfit_fit(ld('lfit_test_poly2.dat'),'a=1,b=2,c=0.5','x:1,y:2','a+b*x+c*x*x','y',method='fima')",
    None))

results.append(mk('M15','FIMA mock',
    "r=lfit_fit(ld('lfit_test_poly2.dat'),'a=1,b=2,c=0.5','x:1,y:2','a+b*x+c*x*x','y',method='fima',parameters='montecarlo',seed=12345,mc_iterations=50)",
    None, timeout=60))

results.append(mk('M16','FIMA derived',
    "r=lfit_fit(ld('lfit_test_poly2.dat'),'a=1,b=2,c=0.5','x:1,y:2','a+b*x+c*x*x','y',method='fima',derived_variables='sum=a+b+c')",
    None))

results.append(mk('M17','MCMC weighted',
    "r=lfit_fit(ld('lfit_test_weighted.dat'),'a=1:0.2,b=1:0.2,c=0.1:0.1','x:1,y:2,e:3','a+b*x+c*x*x','y',error='e',method='mcmc',parameters='accepted',seed=12345,mc_iterations=50)",
    None, timeout=60))

results.append(mk('M18','XMMC weighted',
    "r=lfit_fit(ld('lfit_test_weighted.dat'),'a=1:0.2,b=1:0.2,c=0.1:0.1','x:1,y:2,e:3','a+b*x+c*x*x','y',error='e',method='xmmc',parameters='skip',seed=12345,mc_iterations=50)",
    None, timeout=60))

results.append(mk('M19','EMCE outlier',
    "r=lfit_fit(ld('lfit_test_outlier.dat'),'a,b,c','x:1,y:2','a+b*x+c*x*x','y',method='emce',parameters='clls',seed=12345,mc_iterations=20)",
    None, timeout=60))

for num, name, status, msg, dt in results[19:]:
    print(f"  [{status:7s}] [{num}] {name} ({dt:.1f}s) {msg[:80] if status!='PASS' else ''}")

print()
print("=" * 60)
print("longtest_lfit.sh (9 tests)")
print("=" * 60)

MCMC_N = 20000
XMMC_N = 10000
EMCE_N = 2000
FIMA_N = 20000

long_tests = [
    ('L01','long MCMC accepted',
     f"r=lfit_fit(ld('lfit_test_poly2.dat'),'a=1:0.2,b=1:0.2,c=0.1:0.1','x:1,y:2','a+b*x+c*x*x','y',method='mcmc',parameters='accepted',seed=12345,mc_iterations={MCMC_N})", 300),
    ('L02','long MCMC gibbs weighted',
     f"r=lfit_fit(ld('lfit_test_weighted.dat'),'a=1:0.2,b=1:0.2,c=0.1:0.1','x:1,y:2,e:3','a+b*x+c*x*x','y',error='e',method='mcmc',parameters='gibbs',seed=12345,mc_iterations={MCMC_N})", 300),
    ('L03','long XMMC adaptive',
     f"r=lfit_fit(ld('lfit_test_poly2.dat'),'a=1:0.2,b=1:0.2,c=0.1:0.1','x:1,y:2','a+b*x+c*x*x','y',method='xmmc',parameters='adaptive,window=20',seed=12345,mc_iterations={XMMC_N})", 300),
    ('L04','long XMMC skip weighted',
     f"r=lfit_fit(ld('lfit_test_weighted.dat'),'a=1:0.2,b=1:0.2,c=0.1:0.1','x:1,y:2,e:3','a+b*x+c*x*x','y',error='e',method='xmmc',parameters='skip',seed=12345,mc_iterations={XMMC_N})", 300),
    ('L05','long EMCE clls',
     f"r=lfit_fit(ld('lfit_test_poly2.dat'),'a,b,c','x:1,y:2','a+b*x+c*x*x','y',method='emce',parameters='clls',perturbations='0.01',seed=12345,mc_iterations={EMCE_N})", 300),
    ('L06','long EMCE dhsx',
     f"r=lfit_fit(ld('lfit_test_poly2.dat'),'a=1:0.5,b=1:0.5,c=0.1:0.1','x:1,y:2','a+b*x+c*x*x','y',method='emce',parameters='dhsx',perturbations='0.01',seed=12345,mc_iterations={EMCE_N})", 300),
    ('L07','long FIMA mock',
     f"r=lfit_fit(ld('lfit_test_poly2.dat'),'a=1,b=2,c=0.5','x:1,y:2','a+b*x+c*x*x','y',method='fima',parameters='montecarlo',derived_variables='sum=a+b+c',seed=12345,mc_iterations={FIMA_N})", 300),
    ('L08','long pairs dx MCMC',
     f"r=lfit_fit(ld('lfit_test_pairs.dat'),'a=0:1,b=0:0.01,c=0:0.01','xr:1,yr:2,xo:3,yo:4,mr:5,mo:6','a+b*xr+c*yr','xo-xr',method='mcmc',parameters='accepted',seed=12345,mc_iterations={MCMC_N})", 300),
    ('L09','long pairs dy MCMC',
     f"r=lfit_fit(ld('lfit_test_pairs.dat'),'a=0:1,b=0:0.01,c=0:0.01','xr:1,yr:2,xo:3,yo:4,mr:5,mo:6','a+b*xr+c*yr','yo-yr',method='mcmc',parameters='accepted',seed=12345,mc_iterations={MCMC_N})", 300),
]

for num, name, call, timeout in long_tests:
    r = mk(num, name, call, None, timeout=timeout)
    results.append(r)
    status, msg, dt = r[2], r[3], r[4]
    print(f"  [{status:7s}] [{num}] {name} ({dt:.1f}s) {msg[:80] if status!='PASS' else ''}")
    sys.stdout.flush()

print()
print("=" * 60)
passed = sum(1 for r in results if r[2]=='PASS')
failed = sum(1 for r in results if r[2]=='FAIL')
crashed = sum(1 for r in results if r[2]=='CRASH')
timedout = sum(1 for r in results if r[2]=='TIMEOUT')
total = len(results)
print(f"TOTAL: {passed} PASS, {failed} FAIL, {crashed} CRASH, {timedout} TIMEOUT / {total}")
total_time = sum(r[4] for r in results)
print(f"Total time: {total_time:.0f}s ({total_time/60:.1f}min)")
if failed + crashed > 0:
    print("\nFAILED/CRASHED tests:")
    for r in results:
        if r[2] in ('FAIL','CRASH'):
            print(f"  [{r[0]}] {r[1]}: {r[3][:120]}")
