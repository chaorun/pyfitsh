/*****************************************************************************/
/* ficombine_core.c — memory-based image combination (no FITS I/O)           */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "fitsh.h"
#include "image.h"
#include "mask.h"
#include "tensor.h"
#include "statistics.h"
#include "combine.h"

#include "ficombine_core.h"

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

int combine_parse_mode_simple(char *modstr, combine_parameters *cp)
{
    char *tmp, *tok;
    double sigma_val = 0.0;

    if (modstr == NULL || *modstr == '\0')
        return 0;

    combine_parameters_reset(cp);
    cp->mode = COM_MODE_MEAN;

    tmp = strdup(modstr);
    tok = strtok(tmp, ",");
    while (tok != NULL) {
        while (isspace((unsigned char)*tok)) tok++;

        if (strcmp(tok, "mean") == 0)
            cp->mode = COM_MODE_MEAN;
        else if (strcmp(tok, "median") == 0)
            cp->mode = COM_MODE_MEDIAN;
        else if (strcmp(tok, "min") == 0 || strcmp(tok, "minimum") == 0)
            cp->mode = COM_MODE_MIN;
        else if (strcmp(tok, "max") == 0 || strcmp(tok, "maximum") == 0)
            cp->mode = COM_MODE_MAX;
        else if (strcmp(tok, "rejmed") == 0 || strcmp(tok, "rejection") == 0)
            cp->mode = COM_MODE_REJ_MEDIAN;
        else if (strcmp(tok, "rejmean") == 0)
            cp->mode = COM_MODE_REJ_MEAN;
        else if (strcmp(tok, "truncated") == 0)
            cp->mode = COM_MODE_TRUNC_MEAN;
        else if (strcmp(tok, "winsorized") == 0)
            cp->mode = COM_MODE_WINS_MEAN;
        else if (strcmp(tok, "sum") == 0)
            cp->mode = COM_MODE_SUM;
        else if (strcmp(tok, "squaresum") == 0)
            cp->mode = COM_MODE_SQSUM;
        else if (strcmp(tok, "scatter") == 0 || strcmp(tok, "stddev") == 0)
            cp->mode = COM_MODE_SCT;
        else if (strcmp(tok, "or") == 0)
            cp->logicalmethod = 0;
        else if (strcmp(tok, "and") == 0)
            cp->logicalmethod = 1;
        else if (strcmp(tok, "ignorenegative") == 0)
            cp->ignore_flag |= COM_IGNORE_NEGATIVE;
        else if (strcmp(tok, "ignorezero") == 0)
            cp->ignore_flag |= COM_IGNORE_ZERO;
        else if (strcmp(tok, "ignorepositive") == 0)
            cp->ignore_flag |= COM_IGNORE_POSITIVE;
        else if (strncmp(tok, "iterations=", 11) == 0)
            cp->niter = atoi(tok + 11);
        else if (strncmp(tok, "lower=", 6) == 0)
            cp->lower = atof(tok + 6);
        else if (strncmp(tok, "upper=", 6) == 0)
            cp->upper = atof(tok + 6);
        else if (strncmp(tok, "sigma=", 6) == 0) {
            sigma_val = atof(tok + 6);
            cp->lower = cp->upper = sigma_val;
        }
        else if (strncmp(tok, "lowest=", 7) == 0)
            cp->lowest = atoi(tok + 7);
        else if (strncmp(tok, "highest=", 8) == 0)
            cp->highest = atoi(tok + 8);
        else if (strncmp(tok, "discard=", 8) == 0) {
            int d = atoi(tok + 8);
            cp->lowest = cp->highest = d;
        }

        tok = strtok(NULL, ",");
    }
    free(tmp);
    return 0;
}

/*****************************************************************************/

ficombine_result *ficombine_apply(
    double **images, char **masks,
    int nimg, int sx, int sy,
    combine_parameters *cp,
    int apply_mask)
{
    int i, j;
    ficombine_result *r;
    char *outmask_line;
    double **line_ptrs;
    char **mask_line_ptrs;
    double *out_line;

    r = (ficombine_result *)calloc(1, sizeof(ficombine_result));
    if (r == NULL) return NULL;

    r->sx = sx;
    r->sy = sy;
    r->data = data_alloc(sx, sy);
    r->mask = mask_create_empty(sx, sy);
    if (r->data == NULL || r->mask == NULL) {
        ficombine_result_free(r);
        return NULL;
    }

    out_line = (double *)malloc(sizeof(double) * sx);
    outmask_line = (char *)malloc(sizeof(char) * sx);
    if (out_line == NULL || outmask_line == NULL) {
        ficombine_result_free(r);
        return NULL;
    }

    line_ptrs = (double **)malloc(sizeof(double *) * (size_t)nimg);
    mask_line_ptrs = (char **)malloc(sizeof(char *) * (size_t)nimg);

    for (i = 0; i < sy; i++) {
        for (j = 0; j < nimg; j++) {
            line_ptrs[j] = images[j] + (size_t)i * (size_t)sx;
            if (masks != NULL)
                mask_line_ptrs[j] = masks[j] + (size_t)i * (size_t)sx;
            else
                mask_line_ptrs[j] = NULL;
        }
        combine_lines(line_ptrs, nimg, sx, out_line, cp,
                      (masks != NULL) ? mask_line_ptrs : NULL,
                      outmask_line);
        memcpy(r->data[i], out_line, sizeof(double) * sx);
        memcpy(r->mask[i], outmask_line, sizeof(char) * sx);
    }

    free(out_line);
    free(outmask_line);
    free(line_ptrs);
    free(mask_line_ptrs);

    if (apply_mask) {
        for (i = 0; i < sy; i++) {
            for (j = 0; j < sx; j++) {
                if (r->mask[i][j])
                    r->data[i][j] = 0.0;
            }
        }
    }

    r->history = strdup("ficombine v" FICOMBINE_VERSION);
    return r;
}

void ficombine_result_free(ficombine_result *r)
{
    if (r == NULL) return;
    if (r->data != NULL) data_free(r->data);
    if (r->mask != NULL) mask_free(r->mask);
    if (r->history != NULL) free(r->history);
    free(r);
}

/*****************************************************************************/
