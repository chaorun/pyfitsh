"""Ficonv — kernel fitting + image convolution (Python wrapper)."""


class Ficonv:
    """Kernel fitting + image convolution (ficonv).

    Parameters (matching ficonv --long-help keywords):
    kernel              : -k/--kernel       kernel specification string
    iterations          : -n/--iterations   number of rejection iterations
    rejection_level     : -s/--rejection-level  rejection sigma
    masked              : -m/--masked       masked fitting mode
    weighted            : -w/--weighted     weighted fitting (implies masked)
    background_iterative: -b/--bg-iterative background-iterative fitting
    divide              : -d/--divide       block divide factor
    gain                : -g/--gain         detector gain (e-/ADU)
    verbose             : --verbose         print progress
    unity_kernels       : -u/--unity-kernels normalize identity kernel to flux=1

    Usage:
        fc = Ficonv(kernel='gauss 1.0 2', iterations=3)
        convolved, mask = fc.fit(ref_data, img_data[, ref_mask, img_mask])
    """

    def __init__(self, kernel=None, iterations=0, rejection_level=3.0,
                 masked=False, weighted=False, background_iterative=False,
                 divide=32, gain=1.0, verbose=False, unity_kernels=False):
        from pyfitsh.fitsh_cy import Ficonv_ConvolutionOp
        self._op = Ficonv_ConvolutionOp(kernel, iterations, rejection_level,
                                   masked=masked, weighted=weighted,
                                   background_iterative=background_iterative,
                                   divide=divide, gain=gain,
                                   verbose=verbose, unity_kernels=unity_kernels)

    def fit(self, ref_data, img_data, ref_mask=None, img_mask=None,
            stamps=None, add_to=None, maskinfo=None,
            output_subtracted=False, output_kernel_list=False,
            kernel_dict=None):
        return self._op.fit(ref_data, img_data, ref_mask, img_mask,
                            stamps, add_to, maskinfo,
                            output_subtracted, output_kernel_list,
                            kernel_dict)
