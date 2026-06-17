#ifndef DUMP_HELPER_H
#define DUMP_HELPER_H

#include <stdio.h>
#include <string.h>

static void _dump_do_sub_header(FILE *df, int phase) {
    char magic[8] = "FPDUMP01";
    int ver = 1;
    fwrite(magic, 1, 8, df);
    fwrite(&ver, 4, 1, df);
    fwrite(&phase, 4, 1, df);
}

static void _dump_do_sub_input(
    FILE *df,
    int sx, int sy, double **data, char **mask,
    photstar *ps, int np,
    apgeom *inaps, int numap,
    spatialgain *sg, double correlation_length,
    xphotpar *xpp,
    kernellist *klist, int normalize_kernel)
{
    int i, j, n, iv, nv, has, imgsz;

    _dump_do_sub_header(df, 0);

    fwrite(&sx, 4, 1, df);
    fwrite(&sy, 4, 1, df);
    for (i = 0; i < sy; i++) fwrite(data[i], 8, sx, df);
    for (i = 0; i < sy; i++) fwrite(mask[i], 1, sx, df);

    fwrite(&np, 4, 1, df);
    for (i = 0; i < np; i++) {
        fwrite(&ps[i].x, 8, 1, df);
        fwrite(&ps[i].y, 8, 1, df);
        n = ps[i].n;
        fwrite(&n, 4, 1, df);
        fwrite(&ps[i].ninap, 4, 1, df);
        fwrite(&ps[i].use_ref, 4, 1, df);
        fwrite(&ps[i].ref_mag, 8, 1, df);
        fwrite(&ps[i].ref_col, 8, 1, df);
        fwrite(&ps[i].ref_err, 8, 1, df);

        iv = ps[i].id ? (int)strlen(ps[i].id) : 0;
        fwrite(&iv, 4, 1, df);
        if (iv) fwrite(ps[i].id, 1, iv, df);

        for (j = 0; j < n; j++) {
            fwrite(&ps[i].inaps[j].r0, 8, 1, df);
            fwrite(&ps[i].inaps[j].ra, 8, 1, df);
            fwrite(&ps[i].inaps[j].da, 8, 1, df);
        }

        has = ps[i].rfflux != NULL ? 1 : 0;
        fwrite(&has, 4, 1, df);
        if (has) {
            for (j = 0; j < n; j++) {
                fwrite(&ps[i].rfflux[j].flux, 8, 1, df);
                fwrite(&ps[i].rfflux[j].fluxerr, 8, 1, df);
                fwrite(&ps[i].rfflux[j].flag, 4, 1, df);
                fwrite(&ps[i].rfflux[j].bgarea, 8, 1, df);
                fwrite(&ps[i].rfflux[j].bgflux, 8, 1, df);
                fwrite(&ps[i].rfflux[j].bgmedian, 8, 1, df);
                fwrite(&ps[i].rfflux[j].bgsigma, 8, 1, df);
                fwrite(&ps[i].rfflux[j].mag, 8, 1, df);
                fwrite(&ps[i].rfflux[j].magerr, 8, 1, df);
                fwrite(&ps[i].rfflux[j].rtot, 4, 1, df);
                fwrite(&ps[i].rfflux[j].rbad, 4, 1, df);
                fwrite(&ps[i].rfflux[j].rign, 4, 1, df);
                fwrite(&ps[i].rfflux[j].atot, 4, 1, df);
                fwrite(&ps[i].rfflux[j].abad, 4, 1, df);
            }
        }
    }

    iv = inaps != NULL ? 1 : 0;
    fwrite(&iv, 4, 1, df);
    if (iv) {
        fwrite(&numap, 4, 1, df);
        for (j = 0; j < numap; j++) {
            fwrite(&inaps[j].r0, 8, 1, df);
            fwrite(&inaps[j].ra, 8, 1, df);
            fwrite(&inaps[j].da, 8, 1, df);
        }
    }

    fwrite(&sg->order, 4, 1, df);
    fwrite(&sg->vmin, 8, 1, df);
    nv = (sg->order + 1) * (sg->order + 2) / 2;
    fwrite(&nv, 4, 1, df);
    fwrite(sg->coeff, 8, nv, df);

    fwrite(&correlation_length, 8, 1, df);
    fwrite(&xpp->disjoint_radius, 8, 1, df);
    fwrite(&xpp->use_biquad, 4, 1, df);
    fwrite(&xpp->bgm.type, 4, 1, df);
    fwrite(&xpp->bgm.scatter, 4, 1, df);
    fwrite(&xpp->bgm.rejniter, 4, 1, df);
    fwrite(&xpp->bgm.rejlower, 8, 1, df);
    fwrite(&xpp->bgm.rejupper, 8, 1, df);
    fwrite(&xpp->maskignore, 4, 1, df);
    fwrite(&xpp->use_sky, 4, 1, df);
    fwrite(&xpp->sky, 8, 1, df);

    iv = klist != NULL ? 1 : 0;
    fwrite(&iv, 4, 1, df);
    if (iv) {
        fwrite(&klist->nkernel, 4, 1, df);
        fwrite(&klist->ox, 8, 1, df);
        fwrite(&klist->oy, 8, 1, df);
        fwrite(&klist->scale, 8, 1, df);
        fwrite(&klist->type, 4, 1, df);
        for (j = 0; j < klist->nkernel; j++) {
            fwrite(&klist->kernels[j].type, 4, 1, df);
            fwrite(&klist->kernels[j].sigma, 8, 1, df);
            fwrite(&klist->kernels[j].bx, 4, 1, df);
            fwrite(&klist->kernels[j].by, 4, 1, df);
            fwrite(&klist->kernels[j].hsize, 4, 1, df);
            fwrite(&klist->kernels[j].order, 4, 1, df);
            fwrite(&klist->kernels[j].flag, 4, 1, df);
            fwrite(&klist->kernels[j].target, 4, 1, df);
            fwrite(&klist->kernels[j].offset, 8, 1, df);
            nv = (klist->kernels[j].order + 1) * (klist->kernels[j].order + 2) / 2;
            fwrite(&nv, 4, 1, df);
            fwrite(klist->kernels[j].coeff, 8, nv, df);
            has = klist->kernels[j].image != NULL ? 1 : 0;
            fwrite(&has, 4, 1, df);
            if (has) {
                imgsz = 2 * klist->kernels[j].hsize + 1;
                fwrite(&imgsz, 4, 1, df);
                for (iv = 0; iv < imgsz; iv++)
                    fwrite(klist->kernels[j].image[iv], 8, imgsz, df);
            }
        }
    }

    fwrite(&normalize_kernel, 4, 1, df);
}

static void _dump_do_sub_output(FILE *df, photstar *ps, int np)
{
    int i, j, n;
    _dump_do_sub_header(df, 1);

    for (i = 0; i < np; i++) {
        n = ps[i].n;
        fwrite(&n, 4, 1, df);
        for (j = 0; j < n; j++) {
            fwrite(&ps[i].fluxes[j].flux, 8, 1, df);
            fwrite(&ps[i].fluxes[j].fluxerr, 8, 1, df);
            fwrite(&ps[i].fluxes[j].flag, 4, 1, df);
            fwrite(&ps[i].fluxes[j].bgarea, 8, 1, df);
            fwrite(&ps[i].fluxes[j].bgflux, 8, 1, df);
            fwrite(&ps[i].fluxes[j].bgmedian, 8, 1, df);
            fwrite(&ps[i].fluxes[j].bgsigma, 8, 1, df);
            fwrite(&ps[i].fluxes[j].mag, 8, 1, df);
            fwrite(&ps[i].fluxes[j].magerr, 8, 1, df);
            fwrite(&ps[i].fluxes[j].rtot, 4, 1, df);
            fwrite(&ps[i].fluxes[j].rbad, 4, 1, df);
            fwrite(&ps[i].fluxes[j].rign, 4, 1, df);
            fwrite(&ps[i].fluxes[j].atot, 4, 1, df);
            fwrite(&ps[i].fluxes[j].abad, 4, 1, df);
            fwrite(&ps[i].fluxes[j].cntr_x, 8, 1, df);
            fwrite(&ps[i].fluxes[j].cntr_y, 8, 1, df);
            fwrite(&ps[i].fluxes[j].cntr_width, 8, 1, df);
            fwrite(&ps[i].fluxes[j].cntr_w_d, 8, 1, df);
            fwrite(&ps[i].fluxes[j].cntr_w_k, 8, 1, df);
            fwrite(&ps[i].fluxes[j].cntr_x_err, 8, 1, df);
            fwrite(&ps[i].fluxes[j].cntr_y_err, 8, 1, df);
            fwrite(&ps[i].fluxes[j].cntr_w_err, 8, 1, df);
        }
    }
}

#endif
