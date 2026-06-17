# utils.pxi — misc helpers
from pyfitsh.utils import parse_maskinfo  # noqa: F401

def parse_trans_file(path):
    """Parse a .trans file and return a dict suitable for Grtrans.

    Usage:
        gt = Grtrans(parse_trans_file('grb.trans'))
    """
    with open(path) as f:
        return Grtrans.from_trans_string(f.read()).to_dict()
