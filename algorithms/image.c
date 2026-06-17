/*****************************************************************************/
/* image.c — image_duplicate, image_free, image_draw_line (2D only)          */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "image.h"

/*****************************************************************************/

int image_duplicate(image *ret, image *src, int flag)
{
    int i;

    ret->sx = src->sx;
    ret->sy = src->sy;
    ret->bit = src->bit;

    if (src->data && src->sx > 0 && src->sy > 0) {
        ret->data = malloc(src->sy * sizeof(double *));
        if (!ret->data) return 1;
        for (i = 0; i < src->sy; i++) {
            ret->data[i] = malloc(src->sx * sizeof(double));
            if (!ret->data[i]) return 1;
            if (flag)
                memcpy(ret->data[i], src->data[i], src->sx * sizeof(double));
            else
                memset(ret->data[i], 0, src->sx * sizeof(double));
        }
    } else {
        ret->data = NULL;
        ret->sx = ret->sy = 0;
    }
    return 0;
}

/*****************************************************************************/

void image_free(image *img)
{
    int i;
    if (img->data) {
        for (i = 0; i < img->sy; i++)
            free(img->data[i]);
        free(img->data);
        img->data = NULL;
        img->sx = img->sy = 0;
    }
}

/*****************************************************************************/

int image_draw_line(image *img, int x, int y, int dx, int dy, double col, int style)
{
    int sx, sy, s, n;

    if (!img || !img->data || img->sx <= 0 || img->sy <= 0)
        return 1;

    style = style & 0xffff;
    style = style | (style << 16);

    if (style & 1) {
        if (x >= 0 && y >= 0 && x < img->sx && y < img->sy)
            img->data[y][x] = col;
    }

    dx -= x;
    dy -= y;
    if (dx == 0 && dy == 0)
        return 0;

    if (dx < 0) { dx = -dx; sx = -1; } else sx = +1;
    if (dy < 0) { dy = -dy; sy = -1; } else sy = +1;

    if (dx >= dy) {
        for (s = dx / 2, n = dx; n; n--) {
            style = (style << 1) | ((style >> 31) & 1);
            x += sx;
            s += dy;
            if (s >= dx) { s -= dx; y += sy; }
            if (style & 1) {
                if (x >= 0 && y >= 0 && x < img->sx && y < img->sy)
                    img->data[y][x] = col;
            }
        }
    } else {
        for (s = dy / 2, n = dy; n; n--) {
            style = (style << 1) | ((style >> 31) & 1);
            y += sy;
            s += dx;
            if (s >= dy) { s -= dy; x += sx; }
            if (style & 1) {
                if (x >= 0 && y >= 0 && x < img->sx && y < img->sy)
                    img->data[y][x] = col;
            }
        }
    }
    return 0;
}
