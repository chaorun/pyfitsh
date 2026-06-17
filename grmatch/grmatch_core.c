/*****************************************************************************/
/* grmatch_cy.c — minimal do_pointmatch extract for Cython                  */
/* Only includes what do_pointmatch plus its static helpers need.           */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

#include "grmatch_core.h"
#include "math/poly.h"
#include "math/polyfit.h"
#include "math/tpoint.h"
#include "math/trimatch.h"
#include "math/convexhull.h"

/*****************************************************************************/

int   is_verbose = 0;
int   is_comment = 0;
char *progbasename = "grmatch_cy";

/*****************************************************************************/
/* Static dump helpers (from grmatch_lib.c)                                  */
/*****************************************************************************/

static FILE *dump_open_file(char *stage, char *ext)
{
    char fname[256];
    sprintf(fname, "dump/do_pointmatch_%s.%s", stage, ext);
    return fopen(fname, ext[0] == 'b' ? "wb" : "w");
}

static void dump_write_string(FILE *fb, char *s)
{
    int n;
    if (s == NULL)
    {
        n = -1;
        fwrite(&n, sizeof(n), 1, fb);
    }
    else
    {
        n = (int)strlen(s);
        fwrite(&n, sizeof(n), 1, fb);
        fwrite(s, 1, n, fb);
    }
}

static void dump_pointmatch_stage(char *stage, iline *refls, int nref,
    iline *inpls, int ninp, matchpointtune *mptp, int order,
    double **vfits, matchpointstat *mps, cphit *hits, int nhit)
{
    FILE *ft, *fb;
    int i, nvar;
    nvar = (order + 1) * (order + 2) / 2;

    ft = dump_open_file(stage, "txt");
    fb = dump_open_file(stage, "bin");
    if (ft != NULL) { fprintf(ft, "stage=%s\n", stage); fclose(ft); }
    if (fb != NULL) { fwrite(stage, 1, strlen(stage), fb); fclose(fb); }

    ft = dump_open_file(stage, "mptp.txt");
    fb = dump_open_file(stage, "mptp.bin");
    if (ft != NULL)
    {
        fprintf(ft, "nmiter=%d friter=%d rejlevel=%g\n", mptp->nmiter, mptp->friter, mptp->rejlevel);
        fprintf(ft, "maxdist=%g unitarity=%g parity=%d ttype=%d use_ordering=%d maxnum_ref=%d maxnum_inp=%d\n",
            mptp->maxdist, mptp->unitarity, mptp->parity, mptp->ttype, mptp->use_ordering, mptp->maxnum_ref, mptp->maxnum_inp);
        fprintf(ft, "wcat=%d w_magnitude=%d wpower=%g is_centering=%d refcx=%g refcy=%g inpcx=%g inpcy=%g maxcenterdist=%g hintorder=%d\n",
            mptp->wcat, mptp->w_magnitude, mptp->wpower, mptp->is_centering, mptp->refcx, mptp->refcy, mptp->inpcx, mptp->inpcy, mptp->maxcenterdist, mptp->hintorder);
        fclose(ft);
    }
    if (fb != NULL) { fwrite(mptp, sizeof(*mptp), 1, fb); fclose(fb); }

    ft = dump_open_file(stage, "refls.txt");
    fb = dump_open_file(stage, "refls.bin");
    if (ft != NULL)
    {
        fprintf(ft, "count=%d\n", nref);
        for (i = 0; i < nref; i++)
            fprintf(ft, "[%d] x=%g y=%g ordseq=%g weight=%g id=%s line=%s\n",
                i, refls[i].p.d.x, refls[i].p.d.y, refls[i].ordseq, refls[i].weight,
                refls[i].id ? refls[i].id : "(null)", refls[i].line ? refls[i].line : "(null)");
        fclose(ft);
    }
    if (fb != NULL)
    {
        fwrite(&nref, sizeof(nref), 1, fb);
        for (i = 0; i < nref; i++)
        {
            fwrite(&refls[i].p.d.x, sizeof(double), 1, fb);
            fwrite(&refls[i].p.d.y, sizeof(double), 1, fb);
            fwrite(&refls[i].ordseq, sizeof(double), 1, fb);
            fwrite(&refls[i].weight, sizeof(double), 1, fb);
            dump_write_string(fb, refls[i].id);
            dump_write_string(fb, refls[i].line);
        }
        fclose(fb);
    }

    ft = dump_open_file(stage, "inpls.txt");
    fb = dump_open_file(stage, "inpls.bin");
    if (ft != NULL)
    {
        fprintf(ft, "count=%d\n", ninp);
        for (i = 0; i < ninp; i++)
            fprintf(ft, "[%d] x=%g y=%g ordseq=%g weight=%g id=%s line=%s\n",
                i, inpls[i].p.d.x, inpls[i].p.d.y, inpls[i].ordseq, inpls[i].weight,
                inpls[i].id ? inpls[i].id : "(null)", inpls[i].line ? inpls[i].line : "(null)");
        fclose(ft);
    }
    if (fb != NULL)
    {
        fwrite(&ninp, sizeof(ninp), 1, fb);
        for (i = 0; i < ninp; i++)
        {
            fwrite(&inpls[i].p.d.x, sizeof(double), 1, fb);
            fwrite(&inpls[i].p.d.y, sizeof(double), 1, fb);
            fwrite(&inpls[i].ordseq, sizeof(double), 1, fb);
            fwrite(&inpls[i].weight, sizeof(double), 1, fb);
            dump_write_string(fb, inpls[i].id);
            dump_write_string(fb, inpls[i].line);
        }
        fclose(fb);
    }

    ft = dump_open_file(stage, "hits.txt");
    fb = dump_open_file(stage, "hits.bin");
    if (ft != NULL)
    {
        fprintf(ft, "count=%d\n", nhit);
        if (hits != NULL)
            for (i = 0; i < nhit; i++)
                fprintf(ft, "[%d] distance=%g idx0=%d idx1=%d\n", i, hits[i].distance, hits[i].idx[0], hits[i].idx[1]);
        fclose(ft);
    }
    if (fb != NULL)
    {
        fwrite(&nhit, sizeof(nhit), 1, fb);
        if (hits != NULL)
            for (i = 0; i < nhit; i++)
            {
                fwrite(&hits[i].distance, sizeof(double), 1, fb);
                fwrite(hits[i].idx, sizeof(int), 2, fb);
            }
        fclose(fb);
    }

    ft = dump_open_file(stage, "vfits.txt");
    fb = dump_open_file(stage, "vfits.bin");
    if (ft != NULL)
    {
        fprintf(ft, "order=%d nvar=%d\n", order, nvar);
        if (vfits != NULL)
            for (i = 0; i < nvar; i++)
                fprintf(ft, "%d %g %g\n", i, vfits[0][i], vfits[1][i]);
        fclose(ft);
    }
    if (fb != NULL)
    {
        fwrite(&order, sizeof(order), 1, fb);
        fwrite(&nvar, sizeof(nvar), 1, fb);
        if (vfits != NULL)
        {
            fwrite(vfits[0], sizeof(double), nvar, fb);
            fwrite(vfits[1], sizeof(double), nvar, fb);
        }
        fclose(fb);
    }

    ft = dump_open_file(stage, "mps.txt");
    fb = dump_open_file(stage, "mps.bin");
    if (ft != NULL)
    {
        if (mps != NULL)
            fprintf(ft, "wsigma=%g nsigma=%g unitarity=%g time_total=%g time_trimatch=%g time_symmatch=%g nmiter=%d tri_level=%d hull_coverage=%g\n",
                mps->wsigma, mps->nsigma, mps->unitarity, mps->time_total, mps->time_trimatch, mps->time_symmatch, mps->nmiter, mps->tri_level, mps->hull_coverage);
        fclose(ft);
    }
    if (fb != NULL)
    {
        if (mps != NULL) fwrite(mps, sizeof(*mps), 1, fb);
        fclose(fb);
    }
}

static int match_compare_ordering(const void *vi1, const void *vi2)
{
    if (((iline *)vi1)->ordseq < ((iline *)vi2)->ordseq) return 1;
    else return -1;
}

/*****************************************************************************/
/* do_pointmatch (from grmatch_lib.c, identical to original)                 */
/*****************************************************************************/

int do_pointmatch(iline *refls, int nref, iline *inpls, int ninp,
    matchpointtune *mptp, cphit **rhits, int *rnhit, int order,
    double **vfits, matchpointstat *mps)
{
    tpoint *refps, *rftps, *inpps, *wp;
    point *pfits;
    tpointarr arrref, arrinp;
    cphit *hits;
    int nhit, i, j, k, nmin, ntriref, ntriinp;
    int i_iter, nmiter, tri_level;
    double rejlevel, x, y, w;
    double *xfit, *yfit;
    trimatchpar tmp;
    trimatchlog tml;
    iline *wi;
    double time_trimatch, time_symmatch, time_total;
    clock_t t0, t1;

    dump_pointmatch_stage("before", refls, nref, inpls, ninp, mptp, order, vfits, mps, NULL, 0);

    t0 = clock();

    nmiter = mptp->nmiter;
    if (nmiter <= 0) nmiter = 0;
    rejlevel = mptp->rejlevel;
    if (rejlevel <= 0.0) nmiter = 0;

    if (mptp->use_ordering)
    {
        qsort(inpls, ninp, sizeof(iline), match_compare_ordering);
        qsort(refls, nref, sizeof(iline), match_compare_ordering);
    }

    xfit = vfits[0];
    yfit = vfits[1];

    refps = (tpoint *)malloc(sizeof(tpoint) * nref);
    inpps = (tpoint *)malloc(sizeof(tpoint) * ninp);
    rftps = (tpoint *)malloc(sizeof(tpoint) * nref);

    hits = NULL;

    time_trimatch = 0.0;
    time_symmatch = 0.0;

    tri_level = 0;

    for (i_iter = 0; i_iter <= nmiter; i_iter++)
    {
        if (i_iter <= 0)
        {
            for (i = 0; i < nref; i++)
            {
                wp = &refps[i];
                wp->id = i;
                wp->xcoord = refls[i].p.d.x;
                wp->ycoord = refls[i].p.d.y;
            }
        }
        else
        {
            for (i = 0; i < nref; i++)
            {
                wp = &refps[i];
                wp->id = i;
                x = refls[i].p.d.x;
                y = refls[i].p.d.y;
                wp->xcoord = eval_2d_poly(x, y, order, xfit, 0, 0, 1);
                wp->ycoord = eval_2d_poly(x, y, order, yfit, 0, 0, 1);
            }
        }

        for (i = 0; i < ninp; i++)
        {
            wp = &inpps[i];
            wp->id = i;
            wp->xcoord = inpls[i].p.d.x;
            wp->ycoord = inpls[i].p.d.y;
        }

        hits = NULL;

        tmp.level = mptp->ttype;
        tmp.maxdist = -1.0;
        tmp.unitarity = mptp->unitarity;
        tmp.parity = mptp->parity;

        if (ninp > nref) nmin = nref;
        else nmin = ninp;

        if (mptp->use_ordering && mptp->maxnum_inp > 0 && mptp->maxnum_ref > 0)
        {
            if (nref >= mptp->maxnum_ref && ninp >= mptp->maxnum_inp)
            {
                ntriref = mptp->maxnum_ref;
                ntriinp = mptp->maxnum_inp;
            }
            else
            {
                double nimr, nrmi;
                ntriref = nref;
                ntriinp = ninp;
                nimr = (double)ntriinp * (double)mptp->maxnum_ref;
                nrmi = (double)ntriref * (double)mptp->maxnum_inp;
                if (nimr <= nrmi)
                {
                    ntriref = (int)(nimr / (double)mptp->maxnum_inp);
                }
                else
                {
                    ntriinp = (int)(nrmi / (double)mptp->maxnum_ref);
                }
                if (ntriref > mptp->maxnum_ref) ntriref = mptp->maxnum_ref;
                if (ntriinp > mptp->maxnum_inp) ntriinp = mptp->maxnum_inp;
                if (ntriref > nref) ntriref = nref;
                if (ntriinp > ninp) ntriinp = ninp;
            }
        }
        else if (mptp->use_ordering)
        {
            ntriref = nmin;
            ntriinp = nmin;
        }
        else
        {
            ntriref = nref;
            ntriinp = ninp;
        }

        if (hits != NULL)
        {
            free(hits);
            hits = NULL;
        }

        t0 = clock();
        trimatch(refps, ntriref, inpps, ntriinp, order, &tmp, &hits, &nhit, xfit, yfit, &tml);
        t1 = clock();

        time_trimatch += (double)(t1 - t0) / (double)CLOCKS_PER_SEC;
        if (tml.level_used > tri_level) tri_level = tml.level_used;

        for (i = 0; i < nref; i++)
        {
            x = refps[i].xcoord;
            y = refps[i].ycoord;
            rftps[i].xcoord = eval_2d_poly(x, y, order, xfit, 0, 0, 1);
            rftps[i].ycoord = eval_2d_poly(x, y, order, yfit, 0, 0, 1);
            rftps[i].id = refps[i].id;
        }
        if (hits != NULL) free(hits);

        arrref.points = rftps, arrref.length = nref;
        arrinp.points = inpps, arrinp.length = ninp;

        t0 = clock();
        hits = cpmatch_symmetric(&arrref, &arrinp, &nhit, 0.0, mptp->maxdist);
        t1 = clock();
        time_symmatch += (double)(t1 - t0) / (double)CLOCKS_PER_SEC;

        pfits = (point *)malloc(sizeof(point) * nhit);
        for (i = 0; i < nhit; i++)
        {
            j = hits[i].idx[0];
            pfits[i].x = refls[j].p.d.x;
            pfits[i].y = refls[j].p.d.y;
            pfits[i].weight = 1.0;
        }
        if (mptp->wcat >= 0)
        {
            for (i = 0; i < nhit; i++)
            {
                if (!mptp->wcat) wi = &refls[hits[i].idx[0]];
                else wi = &inpls[hits[i].idx[1]];
                w = wi->weight;
                if (mptp->w_magnitude) { w = exp(-0.4 * M_LN10 * (w - 20.0)); }
                if (mptp->wpower != 1.0) { w = pow(w, mptp->wpower); }
                pfits[i].weight = w;
            }
        }
        for (i = 0; i < nhit; i++) { pfits[i].value = inpls[hits[i].idx[1]].p.d.x; }
        fit_2d_poly(pfits, nhit, order, xfit, 0, 0, 1);
        for (i = 0; i < nhit; i++) { pfits[i].value = inpls[hits[i].idx[1]].p.d.y; }
        fit_2d_poly(pfits, nhit, order, yfit, 0, 0, 1);

        if (mptp->is_centering && mptp->maxcenterdist > 0.0)
        {
            x = eval_2d_poly(mptp->refcx, mptp->refcy, order, xfit, 0, 0, 1) - mptp->inpcx;
            y = eval_2d_poly(mptp->refcx, mptp->refcy, order, yfit, 0, 0, 1) - mptp->inpcy;
            if (x * x + y * y > mptp->maxcenterdist * mptp->maxcenterdist)
                break;
        }
    }

    t1 = clock();
    time_total = (double)(t1 - t0) / (double)CLOCKS_PER_SEC;

    if (mps != NULL)
    {
        mps->time_total = time_total;
        mps->time_trimatch = time_trimatch;
        mps->time_symmatch = time_symmatch;
        mps->nmiter = nmiter + 1;
        mps->tri_level = tri_level;
    }

    if (mps != NULL && nhit > 0)
    {
        double ws, ns, wdd, ndd, nx, ny, dx, dy, dd;
        point *hullpoints;
        int nhull;
        double area_total, area_matched;

        ws = wdd = 0.0;
        ns = ndd = 0.0;
        for (i = 0; i < nhit; i++)
        {
            j = hits[i].idx[0];
            k = hits[i].idx[1];

            if (mptp->wcat < 0) wi = NULL;
            else if (!mptp->wcat) wi = &refls[j];
            else wi = &inpls[k];
            if (wi != NULL)
            {
                w = wi->weight;
                if (mptp->w_magnitude) { w = exp(-0.4 * M_LN10 * (w - 20.0)); }
                if (mptp->wpower != 1.0) { w = pow(w, mptp->wpower); }
            }
            else
                w = 1.0;

            x = refls[j].p.d.x;
            y = refls[j].p.d.y;
            nx = eval_2d_poly(x, y, order, xfit, 0, 0, 1);
            ny = eval_2d_poly(x, y, order, yfit, 0, 0, 1);
            dx = nx - inpls[k].p.d.x;
            dy = ny - inpls[k].p.d.y;
            dd = dx * dx + dy * dy;

            ws += w; wdd += dd * w;
            ns += 1.0; ndd += dd;
        }
        if (ws > 0.0) mps->wsigma = sqrt(wdd / ws);
        else mps->wsigma = 0.0;
        if (ns > 0.0) mps->nsigma = sqrt(ndd / ns);
        else mps->nsigma = 0.0;

        mps->unitarity = calc_2d_unitarity(xfit, yfit, order);

        hullpoints = (point *)malloc(sizeof(point) * ninp);

        for (i = 0; i < ninp; i++)
        {
            hullpoints[i].x = inpls[i].p.d.x;
            hullpoints[i].y = inpls[i].p.d.y;
        }
        nhull = convexhull_compute(hullpoints, ninp);
        if (3 <= nhull)
            area_total = convexhull_area(hullpoints, nhull);
        else
            area_total = 0.0;

        for (i = 0; i < nhit; i++)
        {
            hullpoints[i].x = inpls[hits[i].idx[1]].p.d.x;
            hullpoints[i].y = inpls[hits[i].idx[1]].p.d.y;
        }
        nhull = convexhull_compute(hullpoints, nhit);
        if (3 <= nhull)
            area_matched = convexhull_area(hullpoints, nhull);
        else
            area_matched = 0.0;

        if (0 < area_total)
            mps->hull_coverage = area_matched / area_total;
        else
            mps->hull_coverage = 0.0;

        free(hullpoints);
    }
    else if (mps != NULL)
    {
        mps->wsigma = 0.0;
        mps->nsigma = 0.0;
        mps->unitarity = -1.0;
        mps->hull_coverage = 0.0;
    }

    dump_pointmatch_stage("after", refls, nref, inpls, ninp, mptp, order, vfits, mps, hits, nhit);

    free(rftps);
    free(inpps);
    free(refps);

    if (rhits != NULL) *rhits = hits;
    if (rnhit != NULL) *rnhit = nhit;

    return 0;
}

/*****************************************************************************/
/* Wrapper: flat arrays → iline structs → do_pointmatch → flat arrays       */
/* This is the function exposed to Cython / Python.                         */
/*****************************************************************************/

int do_pointmatch_cy(
    /* reference inputs */
    double *ref_x, double *ref_y, double *ref_ord, int nref,
    /* input inputs */
    double *inp_x, double *inp_y, double *inp_ord, int ninp,
    /* matchpointtune fields */
    int ttype, double maxdist, double unitarity, int parity,
    int use_ordering, int maxnum_ref, int maxnum_inp,
    int nmiter, int friter, double rejlevel,
    int wcat, int w_magnitude, double wpower,
    int is_centering, double refcx, double refcy,
    double inpcx, double inpcy, double maxcenterdist,
    /* polynomial */
    int order,
    /* offset/scale (for output transformation centering) */
    double ox, double oy, double scale,
    /* hint transformation (optional, NULL = unused) */
    double *hint_dxfit, double *hint_dyfit, int hint_order,
    /* outputs  — caller must free: hits_idx0, hits_idx1, vfits_dx, vfits_dy */
    int **hits_idx0, int **hits_idx1, int *rnhit,
    double **vfits_dx, double **vfits_dy,
    /* statistics */
    matchpointstat *mps)
{
    iline *refls, *inpls;
    matchpointtune mptp;
    cphit *hits;
    int nhit, nvar, i;
    double **vfits, *xfit, *yfit;

    /* build refls */
    refls = (iline *)malloc((size_t)nref * sizeof(iline));
    for (i = 0; i < nref; i++)
    {
        memset(&refls[i], 0, sizeof(iline));
        refls[i].p.d.x = ref_x[i];
        refls[i].p.d.y = ref_y[i];
        refls[i].ordseq = ref_ord ? ref_ord[i] : 1.0;
        refls[i].weight = 0.0;
        refls[i].line = NULL;
        refls[i].id = NULL;
    }

    /* build inpls */
    inpls = (iline *)malloc((size_t)ninp * sizeof(iline));
    for (i = 0; i < ninp; i++)
    {
        memset(&inpls[i], 0, sizeof(iline));
        inpls[i].p.d.x = inp_x[i];
        inpls[i].p.d.y = inp_y[i];
        inpls[i].ordseq = inp_ord ? inp_ord[i] : 1.0;
        inpls[i].weight = 0.0;
        inpls[i].line = NULL;
        inpls[i].id = NULL;
    }

    /* store original indices in id field (cast to char*), used to remap
     * after do_pointmatch's internal qsort reorders the arrays */
    for (i = 0; i < nref; i++)
        refls[i].id = (char *)(intptr_t)i;
    for (i = 0; i < ninp; i++)
        inpls[i].id = (char *)(intptr_t)i;

    /* build matchpointtune */
    memset(&mptp, 0, sizeof(mptp));
    mptp.ttype = ttype;
    mptp.maxdist = maxdist;
    mptp.unitarity = unitarity;
    mptp.parity = parity;
    mptp.use_ordering = use_ordering;
    mptp.maxnum_ref = maxnum_ref;
    mptp.maxnum_inp = maxnum_inp;
    mptp.nmiter = nmiter;
    mptp.friter = friter;
    mptp.rejlevel = rejlevel;
    mptp.wcat = wcat;
    mptp.w_magnitude = w_magnitude;
    mptp.wpower = wpower;
    mptp.is_centering = is_centering;
    mptp.refcx = refcx; mptp.refcy = refcy;
    mptp.inpcx = inpcx; mptp.inpcy = inpcy;
    mptp.maxcenterdist = maxcenterdist;
    mptp.htf = NULL;
    mptp.hintorder = -1;

    /* build hint transformation from hint coefficients or offset/scale */
    {
        transformation trf_hint;
        memset(&trf_hint, 0, sizeof(trf_hint));
        if (hint_dxfit != NULL && hint_dyfit != NULL && hint_order > 0) {
            int nv = (hint_order+1)*(hint_order+2)/2;
            trf_hint.type = 1;
            trf_hint.order = hint_order;
            trf_hint.ox = (ox != 0.0 || oy != 0.0 || scale != 1.0) ? ox : 0.0;
            trf_hint.oy = (ox != 0.0 || oy != 0.0 || scale != 1.0) ? oy : 0.0;
            trf_hint.scale = (ox != 0.0 || oy != 0.0 || scale != 1.0) ? scale : 1.0;
            trf_hint.nval = 2;
            trf_hint.vfits = (double **)malloc(sizeof(double *) * 2);
            trf_hint.vfits[0] = (double *)malloc(sizeof(double) * nv);
            trf_hint.vfits[1] = (double *)malloc(sizeof(double) * nv);
            memcpy(trf_hint.vfits[0], hint_dxfit, sizeof(double) * nv);
            memcpy(trf_hint.vfits[1], hint_dyfit, sizeof(double) * nv);
            /* need to heap-allocate so it survives the block */
            mptp.htf = (transformation *)malloc(sizeof(transformation));
            memcpy(mptp.htf, &trf_hint, sizeof(transformation));
            mptp.hintorder = hint_order;
        } else if (ox != 0.0 || oy != 0.0 || scale != 1.0) {
            trf_hint.type = 1;
            trf_hint.order = 1;
            trf_hint.ox = 0.0; trf_hint.oy = 0.0; trf_hint.scale = 1.0;
            trf_hint.nval = 2;
            trf_hint.vfits = (double **)malloc(sizeof(double *) * 2);
            trf_hint.vfits[0] = (double *)malloc(sizeof(double) * 3);
            trf_hint.vfits[1] = (double *)malloc(sizeof(double) * 3);
            trf_hint.vfits[0][0] = -ox / scale;
            trf_hint.vfits[0][1] = 1.0 / scale;
            trf_hint.vfits[0][2] = 0.0;
            trf_hint.vfits[1][0] = -oy / scale;
            trf_hint.vfits[1][1] = 0.0;
            trf_hint.vfits[1][2] = 1.0 / scale;
            mptp.htf = (transformation *)malloc(sizeof(transformation));
            memcpy(mptp.htf, &trf_hint, sizeof(transformation));
            mptp.hintorder = 1;
        }
    }

    /* allocate vfits */
    nvar = (order + 1) * (order + 2) / 2;
    xfit = (double *)calloc((size_t)nvar, sizeof(double));
    yfit = (double *)calloc((size_t)nvar, sizeof(double));
    vfits = (double **)malloc(sizeof(double *) * 2);
    vfits[0] = xfit;
    vfits[1] = yfit;

    /* call do_pointmatch */
    do_pointmatch(refls, nref, inpls, ninp, &mptp, &hits, &nhit, order, vfits, mps);

    /* free hint transformation if we built one */
    if (mptp.htf != NULL) {
        if (mptp.htf->vfits) {
            if (mptp.htf->vfits[0]) free(mptp.htf->vfits[0]);
            if (mptp.htf->vfits[1]) free(mptp.htf->vfits[1]);
            free(mptp.htf->vfits);
        }
        free(mptp.htf);
    }

    /* extract hits */
    *rnhit = nhit;
    if (nhit > 0 && hits != NULL)
    {
        *hits_idx0 = (int *)malloc((size_t)nhit * sizeof(int));
        *hits_idx1 = (int *)malloc((size_t)nhit * sizeof(int));
        for (i = 0; i < nhit; i++)
        {
            (*hits_idx0)[i] = (int)(intptr_t)refls[hits[i].idx[0]].id;
            (*hits_idx1)[i] = (int)(intptr_t)inpls[hits[i].idx[1]].id;
        }
    }
    else
    {
        *hits_idx0 = NULL;
        *hits_idx1 = NULL;
    }

    /* extract vfits */
    *vfits_dx = xfit;
    *vfits_dy = yfit;
    free(vfits);  /* only free the pointer array, keep the data */

    /* cleanup */
    free(refls);
    free(inpls);
    if (hits) free(hits);

    return 0;
}

/*****************************************************************************/
/* coord_match_cy — flat wrapper for do_coordmatch                          */
/*****************************************************************************/

static int match_compare_firstcoord(const void *vi1, const void *vi2)
{
    if (((iline *)vi1)->p.n.x[0] < ((iline *)vi2)->p.n.x[0]) return -1;
    else return 1;
}

static double get_distance(double *x1, double *x2, int dim)
{
    double r = 0.0;
    for (; dim > 0; dim--, x1++, x2++)
        r += ((*x1) - (*x2)) * ((*x1) - (*x2));
    return r;
}

static int coordmatch_search_nearest(iline *ils, int nil, double *x0, int dim)
{
    int best, min, max, mid, i, j, brk, outl, outr, pl, pr;
    double xx, cc, dd, cdist, edist, mindist, maxdist;
    double dist[5], xc0[5], *xc;

    if (ils == NULL || nil <= 0 || x0 == NULL || dim <= 0 || dim > 5) return 0;
    min = mid = 0; max = nil;
    while (max > min) {
        mid = (min + max) / 2;
        cdist = x0[0] - ils[mid].p.n.x[0];
        if (fabs(cdist) < 1e-10) break;
        else if (cdist > 0.0) min = mid + 1;
        else max = mid;
    }
    if (dim == 1) {
        best = mid; mindist = fabs(x0[0] - ils[mid].p.n.x[0]);
        if (mid > 0) { cdist = fabs(x0[0] - ils[mid-1].p.n.x[0]); if (cdist < mindist) { mindist = cdist; best = mid - 1; } }
        if (mid < nil - 1) { cdist = fabs(x0[0] - ils[mid+1].p.n.x[0]); if (cdist < mindist) { mindist = cdist; best = mid + 1; } }
        return best;
    }
    mindist = maxdist = 0.0;
    for (i = 0; i < dim; i++) {
        xc0[i] = ils[mid].p.n.x[i]; dist[i] = fabs(x0[i] - xc0[i]);
        mindist += dist[i] * dist[i]; maxdist += dist[i];
    }
    best = mid;
    for (outl = outr = 1, pl = best - 1, pr = best + 1; outl || outr; pl--, pr++) {
        if (outl) {
            if (pl >= 0) {
                xc = &ils[pl].p.n.x[0]; xx = fabs(x0[0] - xc[0]);
                if (xx < maxdist) {
                    brk = 0; cdist = xx;
                    for (i = 1; i < dim; i++) { cc = fabs(xc[i] - x0[i]); cdist += cc; if (cc <= dist[i]) { dist[i] = cc; brk = 1; } }
                    if (brk) { edist = xx * xx; for (j = 1; j < dim; j++) { dd = xc[j] - x0[j]; edist += dd * dd; }
                        if (edist <= mindist) { best = pl; mindist = edist; if (cdist < maxdist) maxdist = cdist; } }
                }
            } else outl = 0;
        }
        if (outr) {
            if (pr < nil) {
                xc = &ils[pr].p.n.x[0]; xx = fabs(x0[0] - xc[0]);
                if (xx < maxdist) {
                    brk = 0; cdist = xx;
                    for (i = 1; i < dim; i++) { cc = fabs(xc[i] - x0[i]); cdist += cc; if (cc <= dist[i]) { dist[i] = cc; brk = 1; } }
                    if (brk) { edist = xx * xx; for (j = 1; j < dim; j++) { dd = xc[j] - x0[j]; edist += dd * dd; }
                        if (edist <= mindist) { best = pr; mindist = edist; if (cdist < maxdist) maxdist = cdist; } }
                }
            } else outr = 0;
        }
    }
    return best;
}

int coord_match_cy(
    double *ref_x, double *ref_y, int nref,
    double *inp_x, double *inp_y, int ninp,
    double maxdist,
    int **hits_idx0, int **hits_idx1, int *rnhit)
{
    iline *refls, *inpls;
    cphit *hits;
    int i, j, nhit, nmin, *refhs, *inphs;
    double mxd2;

    refls = (iline *)malloc((size_t)nref * sizeof(iline));
    inpls = (iline *)malloc((size_t)ninp * sizeof(iline));
    for (i = 0; i < nref; i++) {
        memset(&refls[i], 0, sizeof(iline)); refls[i].p.n.x[0] = ref_x[i]; refls[i].p.n.x[1] = ref_y[i];
        refls[i].line = (char *)(intptr_t)i;
    }
    for (i = 0; i < ninp; i++) {
        memset(&inpls[i], 0, sizeof(iline)); inpls[i].p.n.x[0] = inp_x[i]; inpls[i].p.n.x[1] = inp_y[i];
        inpls[i].line = (char *)(intptr_t)i;
    }
    qsort(refls, nref, sizeof(iline), match_compare_firstcoord);
    qsort(inpls, ninp, sizeof(iline), match_compare_firstcoord);

    nmin = (nref > ninp ? ninp : nref);
    hits = (cphit *)malloc(sizeof(cphit) * nmin);
    refhs = (int *)malloc(sizeof(int) * nref);
    inphs = (int *)malloc(sizeof(int) * ninp);

    for (i = 0; i < nref; i++) refhs[i] = coordmatch_search_nearest(inpls, ninp, refls[i].p.n.x, 2);
    for (i = 0; i < ninp; i++) inphs[i] = coordmatch_search_nearest(refls, nref, inpls[i].p.n.x, 2);

    nhit = 0; mxd2 = maxdist * maxdist;
    for (i = 0; i < nref; i++) {
        j = refhs[i];
        if (inphs[j] == i && (maxdist < 0 || get_distance(refls[i].p.n.x, inpls[j].p.n.x, 2) <= mxd2)) {
            hits[nhit].idx[0] = i; hits[nhit].idx[1] = j; nhit++;
        }
    }
    free(inphs); free(refhs);

    *rnhit = nhit;
    if (nhit > 0) {
        *hits_idx0 = (int *)malloc(sizeof(int) * nhit);
        *hits_idx1 = (int *)malloc(sizeof(int) * nhit);
        for (i = 0; i < nhit; i++) {
            (*hits_idx0)[i] = (int)(intptr_t)refls[hits[i].idx[0]].line;
            (*hits_idx1)[i] = (int)(intptr_t)inpls[hits[i].idx[1]].line;
        }
    } else { *hits_idx0 = NULL; *hits_idx1 = NULL; }
    free(hits); free(refls); free(inpls);
    return 0;
}

/*****************************************************************************/
/* id_match_cy — flat wrapper for do_idmatch                                */
/*****************************************************************************/

static int stridcmp(char *id1, char *id2)
{
    if (id1 == NULL && id2 == NULL) return 0;
    else if (id1 == NULL) return 1;
    else if (id2 == NULL) return -1;
    else return strcmp(id1, id2);
}

static int id_compare(const void *vi1, const void *vi2)
{
    return stridcmp(((iline *)vi1)->id, ((iline *)vi2)->id);
}

static int search_id_limiters(iline *inpls, int n, char *id, int *rleft, int *rright)
{
    int min, mid, max;
    if (n <= 0) return 1;
    if (rleft == NULL || rright == NULL) return -1;
    min = 0; max = n;
    if (!(0 <= stridcmp(inpls[n-1].id, id))) { *rleft = 1; *rright = 0; return 1; }
    while (max > min + 1) { mid = (min + max) / 2;
        if (0 <= stridcmp(inpls[mid-1].id, id)) max = mid; else min = mid; }
    *rleft = min;
    min = 0; max = n;
    if (!(stridcmp(inpls[0].id, id) <= 0)) { *rleft = 1; *rright = 0; return 1; }
    while (max > min + 1) { mid = (min + max) / 2;
        if (stridcmp(inpls[mid].id, id) <= 0) min = mid; else max = mid; }
    *rright = min;
    if (*rleft > *rright) return 1; else return 0;
}

#define AMBIG_NONE  0
#define AMBIG_FIRST 1
#define AMBIG_ANY   2
#define AMBIG_FULL  3

int id_match_cy(
    char **ref_ids, int nref, char **inp_ids, int ninp, int ambig,
    int **hits_idx0, int **hits_idx1, int *rnhit)
{
    iline *refls, *inpls;
    cphit *hits;
    int i, j, k, nhit, ahit, hr, hi, rl, rr, il, ir;
    char *id;

    refls = (iline *)malloc((size_t)nref * sizeof(iline));
    inpls = (iline *)malloc((size_t)ninp * sizeof(iline));
    for (i = 0; i < nref; i++) { memset(&refls[i], 0, sizeof(iline)); refls[i].id = ref_ids[i]; }
    for (i = 0; i < ninp; i++) { memset(&inpls[i], 0, sizeof(iline)); inpls[i].id = inp_ids[i]; }

    /* store original indices before sorting */
    {
        char *tmp_buf = (char *)malloc((size_t)nref * sizeof(char *));
        for (i = 0; i < nref; i++) ((char **)tmp_buf)[i] = (char *)(intptr_t)i;
        for (i = 0; i < nref; i++) refls[i].line = ((char **)tmp_buf)[i]; /* reuse line ptr */
        for (i = 0; i < ninp; i++) ((char **)tmp_buf)[i] = (char *)(intptr_t)i;
        for (i = 0; i < ninp; i++) inpls[i].line = ((char **)tmp_buf)[i];
        free(tmp_buf);
    }

    qsort(inpls, ninp, sizeof(iline), id_compare);
    qsort(refls, nref, sizeof(iline), id_compare);

    ahit = nref; hits = (cphit *)malloc(sizeof(cphit) * ahit); nhit = 0;

    for (i = 0; i < nref; ) {
        id = refls[i].id; if (id == NULL) { i++; continue; }
        hr = search_id_limiters(refls, nref, id, &rl, &rr);
        if (hr) { i++; continue; }
        hi = search_id_limiters(inpls, ninp, id, &il, &ir);
        if (hi) { i += rr - rl + 1; continue; }
        hr = rr - rl + 1; hi = ir - il + 1;
        switch (ambig) {
            case AMBIG_FIRST: k = 1; break;
            case AMBIG_ANY: k = (hr < hi ? hr : hi); break;
            case AMBIG_FULL: k = hr * hi; break;
            default: k = (hr == 1 && hi == 1) ? 1 : 0; break;
        }
        if (k <= 0) { i += rr - rl + 1; continue; }
        if (nhit + k > ahit) { ahit = nhit + k; hits = (cphit *)realloc(hits, sizeof(cphit) * ahit); }
        switch (ambig) {
            case AMBIG_FIRST: hits[nhit].idx[0] = rl; hits[nhit].idx[1] = il; nhit++; break;
            case AMBIG_ANY: for (j = 0; j < hr && j < hi; j++) { hits[nhit].idx[0] = rl + j; hits[nhit].idx[1] = il + j; nhit++; } break;
            case AMBIG_FULL: for (j = 0; j < hr; j++) { for (k = 0; k < hi; k++) { hits[nhit].idx[0] = rl + j; hits[nhit].idx[1] = il + k; nhit++; } } break;
            default: if (hr == 1 && hi == 1) { hits[nhit].idx[0] = rl; hits[nhit].idx[1] = il; nhit++; } break;
        }
        i += rr - rl + 1;
    }
    *rnhit = nhit;
    if (nhit > 0) {
        *hits_idx0 = (int *)malloc(sizeof(int) * nhit);
        *hits_idx1 = (int *)malloc(sizeof(int) * nhit);
        for (i = 0; i < nhit; i++) {
            (*hits_idx0)[i] = (int)(intptr_t)refls[hits[i].idx[0]].line;
            (*hits_idx1)[i] = (int)(intptr_t)inpls[hits[i].idx[1]].line;
        }
    } else { *hits_idx0 = NULL; *hits_idx1 = NULL; }
    free(hits); free(refls); free(inpls);
    return 0;
}
