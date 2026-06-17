"""Shared utilities: mask info decode/encode, trans file parsing."""


def decode_maskinfo(hdu):
    """Read MASKINFO from FITS HDU header into uint8 mask array.

    Equivalent to the C function mask_read_from_header().
    Returns numpy array shape (NAXIS2, NAXIS1), zero where valid.
    """
    import re
    import numpy as np
    FITS_MASK_MAX = 0x7F

    header = hdu.header
    sx = header['NAXIS1']
    sy = header['NAXIS2']

    lines = []
    for card in header.cards:
        if card.keyword == 'MASKINFO':
            img = card.image
            start = img.index("'") + 1
            end = img.rindex("'")
            lines.append(img[start:end])

    if not lines:
        return np.zeros((sy, sx), dtype=np.uint8)

    mask = np.zeros((sy, sx), dtype=np.uint8)
    return parse_maskinfo(lines, sx, sy)


def parse_maskinfo(maskinfo_lines, sx, sy):
    """Parse MASKINFO-format strings into a uint8 mask array.

    Parameters
    ----------
    maskinfo_lines : str or list of str
        MASKINFO header strings. Format: space-separated entries
        ``x,y:lx,ly``, ``x,y``, or control codes (0/+/−N).
    sx, sy : int
        Image dimensions.

    Returns
    -------
    mask : ndarray of uint8, shape (sy, sx). 0 = valid pixel.
    """
    import numpy as np
    if isinstance(maskinfo_lines, str):
        maskinfo_lines = [maskinfo_lines]

    MASK_MAX = 0x7F
    MASK_DEF = 0x01

    mask = np.zeros((sy, sx), dtype=np.uint8)
    xprev = yprev = 0
    use_diff = False
    data = MASK_DEF

    for line in maskinfo_lines:
        if not line:
            continue
        if not isinstance(line, bytes):
            line = line.encode('ascii', errors='replace')

        parts = line.split(b' ')
        for part in parts:
            if not part:
                continue
            n = part.count(b',') + 1
            nums = []
            for token in part.replace(b':', b',').split(b','):
                try:
                    nums.append(int(token))
                except (ValueError, OverflowError):
                    nums = []
                    break
            if not nums:
                continue

            x = y = lx = ly = 0
            if len(nums) == 1:
                if nums[0] > 0:
                    use_diff = True
                elif nums[0] < 0:
                    data = (-nums[0]) & MASK_MAX
                else:
                    use_diff = False
                x = y = lx = ly = 0
            elif len(nums) == 2:
                x, y = nums[0], nums[1]
                lx = ly = 1
            elif len(nums) == 3:
                x, y = nums[0], nums[1]
                if nums[2] > 1:
                    lx, ly = nums[2], 1
                elif nums[2] < -1:
                    lx, ly = 1, -nums[2]
                else:
                    lx = ly = 1
            elif len(nums) >= 4:
                x, y, lx, ly = nums[0], nums[1], nums[2], nums[3]

            if lx > 0 and ly > 0:
                if use_diff:
                    x += xprev
                    y += yprev
                if x < 0:
                    lx += x
                    x = 0
                if y < 0:
                    ly += y
                    y = 0
                if x + lx > sx:
                    lx = sx - x
                if y + ly > sy:
                    ly = sy - y
                xprev, yprev = x, y
                if ly > 0 and lx > 0:
                    mask[y:y + ly, x:x + lx] |= data

    return mask


def encode_maskinfo(mask):
    """Encode a uint8 mask array into MASKINFO header strings.

    Equivalent to the C function fits_mask_export_as_header().
    Returns a list of strings, each fit for a MASKINFO header card.

    Parameters
    ----------
    mask : ndarray of uint8, shape (sy, sx). Non-zero = masked.

    Returns
    -------
    lines : list of str
        Each string is a MASKINFO value (the part between the quotes).
        Cards should be written as ``MASKINFO = '...'``.
    """
    import numpy as np
    sy, sx = mask.shape
    FITS_MASK_MAX = 0x7F

    # Collect horizontal runs grouped by value
    blocks = {}  # value -> list of (x1, x2, y)
    for y in range(sy):
        x = 0
        while x < sx:
            if mask[y, x] != 0:
                v = int(mask[y, x])
                x1 = x
                while x < sx and mask[y, x] == v:
                    x += 1
                x2 = x
                blocks.setdefault(v, []).append((x1, x2, y))
            else:
                x += 1

    lines = []
    for value, runs in sorted(blocks.items()):
        # Merge consecutive rows with same x-range into rectangles
        rects = []
        for x1, x2, y in runs:
            merged = False
            for ri, (rx1, rx2, ry1, ry2) in enumerate(rects):
                if rx1 == x1 and rx2 == x2 and ry2 + 1 == y:
                    rects[ri] = (rx1, rx2, ry1, y)
                    merged = True
                    break
            if not merged:
                rects.append((x1, x2, y, y))

        parts = [f'-{value & FITS_MASK_MAX}']
        for x1, x2, y1, y2 in rects:
            lx = x2 - x1
            ly = y2 - y1 + 1
            if lx == 1 and ly == 1:
                parts.append(f'{x1},{y1}')
            elif ly == 1:
                parts.append(f'{x1},{y1}:{lx}')
            else:
                parts.append(f'{x1},{y1}:{lx},{ly}')
        lines.append(' '.join(parts))

    return lines if lines else ['']


def load_kernel_spatial_txt_to_dict(filepath):
    """Parse CLI kernel_spatial1.txt file into kernel_dict format.

    Parses the [global] and [kernel] blocks produced by ficonv --output-kernel-list.
    Returns a dict compatible with kernel_dict parameter of photometry_from_raw().
    """
    TYPE_MAP = {'background': 1, 'identity': 2, 'gaussian': 3, 'ddelta': 4}

    with open(filepath, 'r') as f:
        lines = f.readlines()

    result = {'global': {}, 'kernels': []}
    in_global = False
    in_kernel = False
    current_kernel = {}

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            continue

        if stripped == '[global]':
            in_global = True
            in_kernel = False
            continue
        elif stripped == '[end]':
            if in_kernel and current_kernel:
                result['kernels'].append(current_kernel)
                current_kernel = {}
            in_global = False
            in_kernel = False
            continue
        elif stripped.startswith('[kernel]'):
            # extract index from "[kernel] # 0"
            in_kernel = True
            in_global = False
            current_kernel = {}
            try:
                idx_str = stripped.split('#')[-1].strip()
                current_kernel['index'] = int(idx_str)
            except (ValueError, IndexError):
                current_kernel['index'] = len(result['kernels'])
            continue

        if '=' not in stripped:
            continue

        key, _, value = stripped.partition('=')
        key = key.strip()
        value = value.strip()

        if in_global:
            if key == 'type':
                result['global']['type'] = int(value)
            elif key == 'offset':
                parts = value.split(',')
                result['global']['offset'] = [float(parts[0]), float(parts[1])]
            elif key == 'scale':
                result['global']['scale'] = float(value)

        elif in_kernel:
            if key == 'type':
                current_kernel['type'] = TYPE_MAP.get(value, 0)
            elif key == 'order':
                current_kernel['order'] = int(value)
            elif key == 'hsize':
                current_kernel['hsize'] = int(value)
            elif key == 'sigma':
                current_kernel['sigma'] = float(value)
            elif key == 'basis':
                parts = value.split(',')
                current_kernel['bx'] = int(parts[0])
                current_kernel['by'] = int(parts[1])
            elif key == 'coeff':
                current_kernel['coeff'] = [float(x) for x in value.split(',')]

    return result
