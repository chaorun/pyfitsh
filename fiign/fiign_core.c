/*****************************************************************************/
/* fiign_core.c — memory-based pixel mask operations (no FITS I/O)           */
/* Extracted from origincode: fiign.c, common.c, fitsmask.c, maskdraw.c      */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "fitsh.h"
#include "image.h"
#include "mask.h"
#include "tokenize.h"

#include "fiign_core.h"

/*****************************************************************************/

#define MASKDRAW_PIXEL   1
#define MASKDRAW_BLOCK   2
#define MASKDRAW_LINE    3
#define MASKDRAW_CIRCLE  4

/*****************************************************************************/

static double **data_alloc(int sx, int sy)
{
    double **ret;
    double *r0;
    int i;

    ret = (double **)malloc(sizeof(double *) * (size_t)sy +
                            sizeof(double) * (size_t)sx * (size_t)sy);
    if (ret == NULL) return NULL;
    r0 = (double *)(ret + (size_t)sy);
    for (i = 0; i < sy; i++)
        ret[i] = r0 + (size_t)i * (size_t)sx;
    return ret;
}

static void data_free(double **d)
{
    if (d != NULL) free(d);
}

/*****************************************************************************/

int parse_mask_flags_simple(char *list)
{
    char *tmp, *tok;
    int flag = 0;

    if (list == NULL || *list == '\0')
        return -1;

    tmp = strdup(list);
    tok = strtok(tmp, ",");
    while (tok != NULL) {
        while (isspace((unsigned char)*tok)) tok++;
        if (strcmp(tok, "none") == 0 || strcmp(tok, "clear") == 0)
            flag |= MASK_OK;
        else if (strcmp(tok, "fault") == 0)
            flag |= MASK_FAULT;
        else if (strcmp(tok, "hot") == 0)
            flag |= MASK_HOT;
        else if (strcmp(tok, "cosmic") == 0)
            flag |= MASK_COSMIC;
        else if (strcmp(tok, "outer") == 0)
            flag |= MASK_OUTER;
        else if (strcmp(tok, "oversaturated") == 0)
            flag |= MASK_OVERSATURATED;
        else if (strcmp(tok, "leaked") == 0 || strcmp(tok, "bloomed") == 0)
            flag |= MASK_LEAKED;
        else if (strcmp(tok, "saturated") == 0)
            flag |= MASK_SATURATED;
        else if (strcmp(tok, "interpolated") == 0)
            flag |= MASK_INTERPOLATED;
        else if (strcmp(tok, "all") == 0 || strcmp(tok, "bad") == 0)
            flag |= MASK_ALL;
        tok = strtok(NULL, ",");
    }
    free(tmp);
    return (flag > 0) ? (flag & MASK_ALL) : -1;
}

/*****************************************************************************/

static void maskdraw_draw_block(char **mask, int sx, int sy,
    int maskval, int x0, int y0, int bx, int by)
{
    int k, l;
    for (k = 0; k < by; k++) {
        if (y0 + k < 0 || y0 + k >= sy) continue;
        for (l = 0; l < bx; l++) {
            if (x0 + l < 0 || x0 + l >= sx) continue;
            mask[y0 + k][x0 + l] |= maskval;
        }
    }
}

static void maskdraw_draw_line(char **mask, int sx, int sy,
    int maskval, int x1, int y1, int x2, int y2, int width)
{
    int len, lx, ly, i, x, y;
    lx = (x2 > x1 ? x2 - x1 : x1 - x2);
    ly = (y2 > y1 ? y2 - y1 : y1 - y2);
    len = (lx > ly ? lx : ly);
    if (width < 0) width = 0;
    if (len <= 0) {
        maskdraw_draw_block(mask, sx, sy, maskval, x1 - width, y1 - width,
                            2 * width + 1, 2 * width + 1);
        return;
    }
    for (i = 0; i <= len; i++) {
        x = ((2 * x1 + 1) + 2 * (x2 - x1) * i / len) / 2;
        y = ((2 * y1 + 1) + 2 * (y2 - y1) * i / len) / 2;
        maskdraw_draw_block(mask, sx, sy, maskval, x - width, y - width,
                            2 * width + 1, 2 * width + 1);
    }
}

static void maskdraw_draw_circle(char **mask, int sx, int sy,
    int maskval, int x1, int y1, int radius)
{
    int k, l, dx, dy2, r2;
    r2 = radius * (radius + 1);
    for (k = y1 - radius; k <= y1 + radius; k++) {
        if (!(0 <= k && k < sy)) continue;
        dy2 = (k - y1) * (k - y1);
        for (l = x1 - radius; l <= x1 + radius; l++) {
            if (!(0 <= l && l < sx)) continue;
            dx = l - x1;
            if (dy2 + dx * dx <= r2)
                mask[k][l] |= maskval;
        }
    }
}

static int mask_block_draw_one(char **mask, int sx, int sy, char *maskblock)
{
    char *mtmp, *cmd[8], **cc;
    int n, is_invalid, maskval, x0, y0, bx, by, br, w, type;

    if (maskblock == NULL) return 1;
    if (mask == NULL || sx <= 0 || sy <= 0) return -1;

    mtmp = strdup(maskblock);
    n = tokenize_char(mtmp, cmd, ':', 7);
    is_invalid = 0;
    x0 = y0 = 0;
    bx = by = -1;
    maskval = 0;

    if (strcmp(cmd[0], "block") == 0)
        type = MASKDRAW_BLOCK, cc = cmd + 1, n--;
    else if (strcmp(cmd[0], "pixel") == 0)
        type = MASKDRAW_PIXEL, cc = cmd + 1, n--;
    else if (strcmp(cmd[0], "line") == 0)
        type = MASKDRAW_LINE, cc = cmd + 1, n--;
    else if (strcmp(cmd[0], "circle") == 0)
        type = MASKDRAW_CIRCLE, cc = cmd + 1, n--;
    else
        type = MASKDRAW_BLOCK, cc = cmd;

    if (type == MASKDRAW_BLOCK) {
        if (n < 2 || n > 3) is_invalid = 1;
        else if ((maskval = parse_mask_flags_simple(cc[0])) < 0) is_invalid = 2;
        else if (sscanf(cc[1], "%d,%d", &x0, &y0) < 2) is_invalid = 3;
        else if (n >= 3 && sscanf(cc[2], "%d,%d", &bx, &by) < 2) is_invalid = 4;
        if (bx >= 0 && by >= 0) { bx = bx - x0 + 1; by = by - y0 + 1; }
        else { bx = 1; by = 1; }
        if (!is_invalid)
            maskdraw_draw_block(mask, sx, sy, maskval, x0, y0, bx, by);
    } else if (type == MASKDRAW_PIXEL) {
        if (n != 2) is_invalid = 1;
        else if ((maskval = parse_mask_flags_simple(cc[0])) < 0) is_invalid = 2;
        else if (sscanf(cc[1], "%d,%d", &x0, &y0) < 2) is_invalid = 3;
        if (!is_invalid)
            maskdraw_draw_block(mask, sx, sy, maskval, x0, y0, 1, 1);
    } else if (type == MASKDRAW_LINE) {
        w = 0;
        if (n < 3 || n > 4) is_invalid = 1;
        else if ((maskval = parse_mask_flags_simple(cc[0])) < 0) is_invalid = 2;
        else if (sscanf(cc[1], "%d,%d", &x0, &y0) < 2) is_invalid = 3;
        else if (sscanf(cc[2], "%d,%d", &bx, &by) < 2) is_invalid = 3;
        else if (n == 4 && (sscanf(cc[3], "%d", &w) < 1 || w < 0)) is_invalid = 3;
        if (!is_invalid)
            maskdraw_draw_line(mask, sx, sy, maskval, x0, y0, bx, by, w);
    } else if (type == MASKDRAW_CIRCLE) {
        w = 0;
        if (n != 3) is_invalid = 1;
        else if ((maskval = parse_mask_flags_simple(cc[0])) < 0) is_invalid = 2;
        else if (sscanf(cc[1], "%d,%d", &x0, &y0) < 2) is_invalid = 3;
        else if (sscanf(cc[2], "%d", &br) < 1) is_invalid = 3;
        if (!is_invalid)
            maskdraw_draw_circle(mask, sx, sy, maskval, x0, y0, br);
    } else {
        is_invalid = 5;
    }
    free(mtmp);
    return is_invalid;
}

/*****************************************************************************/

static int saturated_mark(double **data, char **mask, int sx, int sy,
    double *sat_img, double param, int method)
{
    int i, j, c;

    if (mask == NULL) return -1;

    if (sat_img != NULL) {
        for (i = 0; i < sy; i++) {
            for (j = 0; j < sx; j++) {
                if (data[i][j] >= sat_img[i * sx + j] * param)
                    mask[i][j] |= MASK_OVERSATURATED;
            }
        }
    } else {
        for (i = 0; i < sy; i++) {
            for (j = 0; j < sx; j++) {
                if (data[i][j] >= param)
                    mask[i][j] |= MASK_OVERSATURATED;
            }
        }
    }

    for (i = 0; i < sy; i++) {
        for (j = 0; j < sx; j++) {
            c = mask[i][j];
            if (!(c & MASK_OVERSATURATED)) continue;
            if (method & 1) {
                if (i > 0)    mask[i - 1][j] |= MASK_LEAKED;
                if (i < sy - 1) mask[i + 1][j] |= MASK_LEAKED;
            }
            if (method & 2) {
                if (j > 0)    mask[i][j - 1] |= MASK_LEAKED;
                if (j < sx - 1) mask[i][j + 1] |= MASK_LEAKED;
            }
        }
    }
    for (i = 0; i < sy; i++) {
        for (j = 0; j < sx; j++) {
            c = mask[i][j];
            if ((c & MASK_OVERSATURATED) && (c & MASK_LEAKED))
                mask[i][j] &= ~MASK_LEAKED;
        }
    }
    return 0;
}

/*****************************************************************************/

static int integerlimit_mark(double **data, char **mask, int sx, int sy,
    int bitpix, int is_corr, int mvlo, int mvhi)
{
    double llo, lhi;
    int k, l;

    if (bitpix < 0) return 0;
    else if (bitpix == 8)   { llo = -128.0;          lhi = 127.0; }
    else if (bitpix == 16)  { llo = -32768.0;        lhi = 32767.0; }
    else if (bitpix == 32)  { llo = -2147483648.0;   lhi = 2147483647.0; }
    else return -1;

    for (k = 0; k < sy; k++) {
        for (l = 0; l < sx; l++) {
            if (data[k][l] < llo) {
                if (mask != NULL) mask[k][l] |= mvlo;
                if (is_corr) data[k][l] = llo;
            } else if (data[k][l] > lhi) {
                if (mask != NULL) mask[k][l] |= mvhi;
                if (is_corr) data[k][l] = lhi;
            } else if (is_corr) {
                data[k][l] = floor(data[k][l]);
            }
        }
    }
    return 0;
}

/*****************************************************************************/

static int cosmics_ignore(double **data, char **mask, int sx, int sy,
    double th_low, double th_high, int is_repl, double skysigma)
{
    int i, j, k, l, hsize, fsize, ii, jj;
    double arr[25], s, s2, sig, w;
    int is_low, is_high;

    if (data == NULL) return 0;
    is_low = is_high = 0;
    if (th_low  > 0.0) is_low  = 1;
    if (th_high > 0.0) is_high = 1;

    hsize = 1;
    fsize = 2 * hsize + 1;

    for (i = hsize; i < sy - hsize; i++) {
        for (j = hsize; j < sx - hsize; j++) {
            l = 0;
            for (k = 0, ii = jj = -hsize; k < fsize * fsize; k++) {
                if (!(mask[i + ii][j + jj] & MASK_FAULT) && !(ii == 0 && jj == 0)) {
                    w = data[i + ii][j + jj];
                    arr[l] = w; l++;
                }
                if (jj == hsize) ii++, jj = -hsize;
                else jj++;
            }
            if (l < fsize * fsize - 1) continue;

            s = s2 = 0.0;
            for (k = 0; k < l; k++) s += arr[k], s2 += arr[k] * arr[k];
            s /= l; s2 /= l; sig = sqrt(s2 - s * s);

            if (skysigma > 0.0 && sig > 2 * skysigma) continue;
            if (skysigma > 0.0) sig = skysigma;

            if (is_low && data[i][j] < s - th_low * sig) {
                mask[i][j] |= MASK_COSMIC;
                if (is_repl) {
                    data[i][j] = s;
                    mask[i][j] |= MASK_INTERPOLATED;
                }
            }
            if (is_high && data[i][j] > s + th_high * sig) {
                mask[i][j] |= MASK_COSMIC;
                if (is_repl) {
                    data[i][j] = s;
                    mask[i][j] |= MASK_INTERPOLATED;
                }
            }
        }
    }
    return 0;
}

/*****************************************************************************/

static void mask_convert_inmem(char **mask, int sx, int sy,
    char **convert_list)
{
    int i, j, k, nmcl;
    int match, value, reset, set, is_any;
    char *cmd[4], *cc;

    if (convert_list == NULL) return;

    for (k = 0; convert_list[k] != NULL; k++) {
        cc = strdup(convert_list[k]);
        int n = tokenize_char(cc, cmd, ':', 4);
        match = value = reset = set = -1;
        is_any = 0;
        if (n == 4) {
            match = parse_mask_flags_simple(cmd[0]);
            if (strcmp(cmd[1], "any") == 0) is_any = 1, value = 0;
            else { is_any = 0; value = parse_mask_flags_simple(cmd[1]); }
            reset = parse_mask_flags_simple(cmd[2]);
            set   = parse_mask_flags_simple(cmd[3]);
        }
        free(cc);
        if (match <= 0 || value < 0 || reset < 0 || set < 0 || (!is_any && value == 0))
            continue;

        for (i = 0; i < sy; i++) {
            for (j = 0; j < sx; j++) {
                int m = mask[i][j];
                if ((m & match) == (value & match))
                    mask[i][j] = (m & (~reset)) | set;
            }
        }
    }
}

/*****************************************************************************/

static char *history_build(double saturation, int leak_method,
    int ignore_nonpos, int ignore_neg, int ignore_zero,
    int ignore_cosmics, int replace_cosmics,
    double th_low, double th_high, double sky_sigma,
    int expand_hsize, int apply_mask, double mask_value,
    int bitpix, char **convert_list, char **mask_block_list)
{
    char buf[4096];
    int pos = 0;

    pos += snprintf(buf + pos, sizeof(buf) - pos,
        "fiign v%s: ", FIIGN_VERSION);

    if (saturation > 0.0)
        pos += snprintf(buf + pos, sizeof(buf) - pos, "saturation=%.1f ", saturation);
    if (leak_method)
        pos += snprintf(buf + pos, sizeof(buf) - pos, "leak=%d ", leak_method);
    if (ignore_nonpos)
        pos += snprintf(buf + pos, sizeof(buf) - pos, "nonpos=1 ");
    if (ignore_neg)
        pos += snprintf(buf + pos, sizeof(buf) - pos, "neg=1 ");
    if (ignore_zero)
        pos += snprintf(buf + pos, sizeof(buf) - pos, "zero=1 ");
    if (ignore_cosmics)
        pos += snprintf(buf + pos, sizeof(buf) - pos,
            "cosmics=1 th_low=%.1f th_high=%.1f ", th_low, th_high);
    if (replace_cosmics)
        pos += snprintf(buf + pos, sizeof(buf) - pos, "replace=1 ");
    if (sky_sigma > 0.0)
        pos += snprintf(buf + pos, sizeof(buf) - pos, "skysigma=%.2f ", sky_sigma);
    if (expand_hsize > 0)
        pos += snprintf(buf + pos, sizeof(buf) - pos, "expand=%d ", expand_hsize);
    if (apply_mask)
        pos += snprintf(buf + pos, sizeof(buf) - pos, "apply_mask=1 mask_value=%.6g ", mask_value);
    if (bitpix > 0)
        pos += snprintf(buf + pos, sizeof(buf) - pos, "bitpix=%d ", bitpix);
    if (convert_list != NULL) {
        pos += snprintf(buf + pos, sizeof(buf) - pos, "convert=[");
        for (int i = 0; convert_list[i] != NULL; i++)
            pos += snprintf(buf + pos, sizeof(buf) - pos, "%s%s", i > 0 ? "," : "", convert_list[i]);
        pos += snprintf(buf + pos, sizeof(buf) - pos, "] ");
    }
    if (mask_block_list != NULL) {
        pos += snprintf(buf + pos, sizeof(buf) - pos, "blocks=[");
        for (int i = 0; mask_block_list[i] != NULL; i++)
            pos += snprintf(buf + pos, sizeof(buf) - pos, "%s%s", i > 0 ? "," : "", mask_block_list[i]);
        pos += snprintf(buf + pos, sizeof(buf) - pos, "] ");
    }
    buf[pos] = '\0';
    return strdup(buf);
}

/*****************************************************************************/

fiign_result *fiign_apply(
    double **img, char **mask_in, int sx, int sy,
    double saturation, double *sat_img,
    int leak_method,
    int ignore_nonpos, int ignore_neg, int ignore_zero,
    int ignore_cosmics, int replace_cosmics,
    double th_low, double th_high, double sky_sigma,
    int expand_hsize,
    int apply_mask, double mask_value,
    int bitpix,
    char **convert_list,
    char **mask_block_list)
{
    int i, j;
    fiign_result *r;
    char **mask;

    r = (fiign_result *)calloc(1, sizeof(fiign_result));
    if (r == NULL) return NULL;

    r->sx = sx;
    r->sy = sy;

    r->data = data_alloc(sx, sy);
    r->mask = mask_create_empty(sx, sy);
    if (r->data == NULL || r->mask == NULL) {
        fiign_result_free(r);
        return NULL;
    }

    for (i = 0; i < sy; i++) {
        memcpy(r->data[i], img[i], sx * sizeof(double));
        if (mask_in != NULL)
            memcpy(r->mask[i], mask_in[i], sx * sizeof(char));
        else
            memset(r->mask[i], 0, sx * sizeof(char));
    }

    mask = r->mask;

    /* 1. mark nans */
    mask_mark_nans(NULL, mask, MASK_NAN);
    {
        double **d = r->data;
        for (i = 0; i < sy; i++) {
            for (j = 0; j < sx; j++) {
                if (!isfinite(d[i][j]))
                    mask[i][j] |= MASK_NAN;
            }
        }
    }

    /* 2. mark saturated pixels */
    if (saturation > 0.0 || sat_img != NULL) {
        saturated_mark(r->data, mask, sx, sy, sat_img,
            sat_img != NULL ? 1.0 : saturation, leak_method);
    }

    /* 3. mask block drawing */
    if (mask_block_list != NULL) {
        for (i = 0; mask_block_list[i] != NULL; i++)
            mask_block_draw_one(mask, sx, sy, mask_block_list[i]);
    }

    /* 4. ignore nonpos / neg / zero */
    if (ignore_nonpos || ignore_neg || ignore_zero) {
        int cm;
        for (i = 0; i < sy; i++) {
            for (j = 0; j < sx; j++) {
                if (r->data[i][j] == 0.0) cm = 1;
                else if (r->data[i][j] < 0.0) cm = 2;
                else cm = 0;
                if (ignore_nonpos && cm != 0) mask[i][j] |= MASK_FAULT;
                else if (ignore_neg && cm == 2) mask[i][j] |= MASK_FAULT;
                else if (ignore_zero && cm == 1) mask[i][j] |= MASK_FAULT;
            }
        }
    }

    /* 5. ignore cosmics */
    if (ignore_cosmics) {
        cosmics_ignore(r->data, mask, sx, sy,
            th_low, th_high, replace_cosmics, sky_sigma);
    }

    /* 6. mask convert */
    if (convert_list != NULL) {
        mask_convert_inmem(mask, sx, sy, convert_list);
    }

    /* 7. mask expand */
    if (expand_hsize > 0) {
        mask_expand_logic(mask, sx, sy, expand_hsize, -1, 0, 0);
    }

    /* 8. apply mask to data */
    if (apply_mask) {
        for (i = 0; i < sy; i++) {
            for (j = 0; j < sx; j++) {
                if (mask[i][j])
                    r->data[i][j] = mask_value;
            }
        }
    }

    /* 9. mark integer limited pixels */
    if (bitpix > 0) {
        integerlimit_mark(r->data, mask, sx, sy, bitpix, 1,
            MASK_OVERSATURATED, MASK_OVERSATURATED);
    }

    /* 10. build history */
    r->history = history_build(saturation, leak_method,
        ignore_nonpos, ignore_neg, ignore_zero,
        ignore_cosmics, replace_cosmics,
        th_low, th_high, sky_sigma,
        expand_hsize, apply_mask, mask_value,
        bitpix, convert_list, mask_block_list);

    return r;
}

void fiign_result_free(fiign_result *r)
{
    if (r == NULL) return;
    if (r->data != NULL) data_free(r->data);
    if (r->mask != NULL) mask_free(r->mask);
    if (r->history != NULL) free(r->history);
    free(r);
}

/*****************************************************************************/
