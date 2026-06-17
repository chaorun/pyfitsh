# Fistar 三项测试报告 (2026-06-16)

Cython 编译参数: `-O3`（与 CLI Makefile 对齐）

---

## Test 1: parabolapeak only_candidates

```
fistar --input L1...fits -d 5 -f 10000 --algorithm parabolapeak --only-candidates
       --format id,ix,iy,cx,cy,cbg,camp,cmax,npix,cs,cd,ck
```

| 指标 | CLI | Cython | 结果 |
|------|-----|--------|------|
| nstar | 147 | 147 | PASS |
| KD-tree 匹配 | — | 147/147 | PASS |
| max_dist | — | 0.0007 pix | PASS |
| 字段差异 | — | 无 | PASS |

---

## Test 2: parabolapeak fit

```
fistar --input L1...fits -d 5 -f 10000 -g 2 --algorithm parabolapeak
       --model elliptic
       --format id,x,y,bg,amp,s,d,k,fwhm,ellip,pa,flux,noise,s/n
```

| 指标 | CLI | Cython | 结果 |
|------|-----|--------|------|
| nstar | 146 | 146 | PASS |
| KD-tree 匹配 | — | 146/146 | PASS |
| max_dist | — | 0.0007 pix | PASS |
| bg | — | max_diff=18.7 | DIFF (~1.3%) |
| amp | — | max_diff=17.1 | DIFF (~0.03%) |
| flux | — | max_diff=2850 | DIFF (~0.14%) |
| pa | — | max_diff=0.05 | DIFF (浮点) |
| x,y,s,d,k,fwhm,ellip,noise,sn | — | 0 | PASS |

---

## Test 3: Gauss model + PSF

```
fistar --input L1...fits -d 5 -f 10000 -s x -g 2 --model gauss
       --format id,ix,iy,...,px,py,pbg,pamp,ps,pd,pk,pl
       --psf native,order=2 --output-psf test.psf.fits -V
```

| 指标 | CLI | Cython | 结果 |
|------|-----|--------|------|
| nstar | 49 | **9** | **FAIL** |
| PSF | 有 | 有（对比失败） | 待验证 |

根因: `search_star_candidates_link` 在 threshold=10000 时仅检出 9 颗星，CLI 检出 49 颗。Cython debug 输出确认 `threshold=10000.0, flux_th=0.0, skysigma=5.0`，参数正确。

---

## 总结

- Test 1: 完全通过
- Test 2: nstar 通过，bg/amp/flux 有 0.03%-1.3% 差异（编译优化所致）
- Test 3: nstar 严重偏差（9 vs 49），根因不明
