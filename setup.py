"""setup.py -- 编译 fitsh_cy Cython 扩展 (自包含)。"""
import os
import numpy as np
import setuptools._distutils.sysconfig as distutilsc
# Force exact same compile/link flags as working dylib build
distutilsc.get_config_vars()['CFLAGS'] = '-O0 -g -fPIC'
distutilsc.get_config_vars()['OPT'] = '-O0 -g'
distutilsc.get_config_vars()['BASECFLAGS'] = ''
distutilsc.get_config_vars()['CC'] = 'gcc'
distutilsc.get_config_vars()['LDSHARED'] = 'gcc -shared -undefined dynamic_lookup'
distutilsc.get_config_vars()['BLDSHARED'] = 'gcc -shared -undefined dynamic_lookup'
# Keep only Python include, strip everything else distutils adds
import sysconfig
py_include = sysconfig.get_path('include')
distutilsc.get_config_vars()['CFLAGS'] += f' -isystem {py_include}'
from setuptools import setup, Extension
from Cython.Build import cythonize
HERE = os.path.dirname(os.path.abspath(__file__))
MATH = os.path.join(HERE, "math")
FIT  = os.path.join(MATH, "fit")
SPL  = os.path.join(MATH, "spline")
IO_DIR = os.path.join(HERE, "io")
LINK = os.path.join(HERE, "link")
INDEX = os.path.join(HERE, "index")
EXPINT = os.path.join(MATH, "expint")
INTERSEC = os.path.join(MATH, "intersec")
ELLIPTIC = os.path.join(MATH, "elliptic")
FIPHOT = os.path.join(HERE, "fiphot")
FICONV = os.path.join(HERE, "ficonv")
FISTAR = os.path.join(HERE, "fistar")
FICALIB = os.path.join(HERE, "ficalib")
FIIGN   = os.path.join(HERE, "fiign")
FITRANS = os.path.join(HERE, "fitrans")
GRMATCH = os.path.join(HERE, "grmatch")
GRTRANS = os.path.join(HERE, "grtrans")
FIRANDOM = os.path.join(HERE, "firandom")
FIARITH  = os.path.join(HERE, "fiarith")
PSN_DIR  = os.path.join(HERE, "psn")
DFT_DIR  = os.path.join(MATH, "dft")
# 所有 C 源码文件 (核心 + 数学库 + 样条 + 投影 + 图像变换 + kernel + 星检测)
sources = [
    "grmatch/grmatch_core.c",
    "algorithms/trans_core.c",
    "algorithms/wcs_core.c",
    "algorithms/polyfit_core.c",
    "fitrans/fitrans_core.c",
    "fitrans/fitrans_ops.c",
    "grtrans/grtrans_core.c",
    "ficonv/ficonv_core.c",
    "fistar/fistar_core.c",
    "fiphot/fiphot_core.c",
    "algorithms/stubs.c",
    "algorithms/str.c",
    "algorithms/kernel-base.c",
    "algorithms/kernel-io.c",
    "algorithms/mask.c",
    "algorithms/projection.c",
    "algorithms/tensor.c",
    "algorithms/fbase.c",
    "algorithms/statistics.c",
    # star search (all 4 algorithms)
    "fistar/star-base.c",
    "fistar/star-cand-lnk.c",
    "fistar/star-cand-pp.c",
    "fistar/star-cand-trb.c",
    "fistar/star-cand-biq.c",
    "fistar/star-draw.c",
    "fistar/star-model.c",
    "fistar/background.c",
    # PSF
    "fistar/psf-base.c",
    "fistar/psf-determine.c",
    # aperture photometry (fiphot)
    "fiphot/apphot.c",
    "fiphot/fiphot-io.c",
    "fiphot/magnitude.c",
    "fiphot/weight-gen.c",
    "fiphot/weight-star.c",
    # link algorithm deps (all needed at link time)
    "link/floodfill.c",
    "link/linkblock.c",
    "link/linkpoint.c",
    # math extras
    "math/expint/expint.c",
    "math/intersec/intersec.c",
    "math/intersec/intersec-cri.c",
    "math/elliptic/elliptic.c",
    "math/elliptic/ntiq.c",
    "math/fit/downhill.c",
    "index/sort.c",
    "index/multiindex.c",
    # existing math
    "math/tpoint.c",
    "math/cpmatch.c",
    "math/delaunay.c",
    "math/trimatch.c",
    "math/poly.c",
    "math/polyfit.c",
    "math/spmatrix.c",
    "math/convexhull.c",
    "math/polygon.c",
    "math/fit/lmfit.c",
    "math/spline/bicubic.c",
    "math/spline/biquad.c",
    "math/spline/biquad-isc.c",
    "math/spline/spline.c",
    # image helpers (no fits lib dependency)
    "algorithms/image.c",
    # io
    "algorithms/tokenize.c",
    # firandom
    "firandom/firandom_core.c",
    "firandom/firandom-eval.c",
    "firandom/random.c",
    # PSN expression parser
    "psn/psn.c",
    "psn/psn-general.c",
    "psn/psn-general-ds.c",
    # fiarith expression evaluator
    "fiarith/evaluate.c",
    "fiarith/limitpix.c",
    "fiarith/maskjoin.c",
    # pbfft (needed by fiarith corr())
    "math/dft/pbfft.c",
]
ext = Extension(
    name="pyfitsh.core",
    sources=[os.path.join(HERE, "core.pyx")]
          + [os.path.join(HERE, s) for s in sources]
           + [os.path.join(HERE, "ficalib", "ficalib_core.c"),
              os.path.join(HERE, "ficalib", "combine.c"),
              os.path.join(HERE, "math", "splinefit.c"),
              os.path.join(HERE, "fiign", "fiign_core.c")],
    include_dirs=[os.path.join(HERE,"algorithms"), HERE, MATH, FIT, SPL, IO_DIR, LINK, INDEX, EXPINT, INTERSEC, ELLIPTIC, FIPHOT, FICONV, FISTAR, FITRANS, GRMATCH, GRTRANS, FIRANDOM, FIARITH, PSN_DIR, DFT_DIR, FICALIB, FIIGN],
    extra_compile_args=["-O3"],
    language="c",
)

# fiarith Cython extension (separate module, same C sources)
fiarith_ext = Extension(
    name="pyfitsh.fiarith.fiarith_core",
    sources=[os.path.join(HERE, "fiarith/fiarith_core.pyx")]
          + [os.path.join(HERE, s) for s in sources],
    include_dirs=[np.get_include(),
                  os.path.join(HERE,"algorithms"), HERE, MATH, FIT, SPL, IO_DIR, LINK, INDEX, EXPINT, INTERSEC, ELLIPTIC, FIPHOT, FICONV, FISTAR, FITRANS, GRMATCH, GRTRANS, FIRANDOM, FIARITH, PSN_DIR, DFT_DIR],
    extra_compile_args=["-O3"],
    language="c",
)

setup(
    name="pyfitsh",
    version="0.01",
    license="GPL-3.0-or-later",
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
    ],
    ext_modules=cythonize(
        [ext, fiarith_ext],
        compiler_directives={"language_level": "3"},
    ),
    zip_safe=False,
)
