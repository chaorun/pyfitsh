#ifndef GRMATCH_LIB_H
#define GRMATCH_LIB_H

#include <stdio.h>
#include <stdarg.h>

#include "transform.h"
#include "math/cpmatch.h"

#define FITSH_GRMATCH_VERSION "0.9e0"

#define MAX_COORDMATCH_DIM 5

#define MATCH_POINTS      1
#define MATCH_IDS         2
#define MATCH_COORDS      3

#define MATCH_READ_COORDS 0x01
#define MATCH_READ_ID     0x02

#define AMBIG_NONE  0
#define AMBIG_FIRST 1
#define AMBIG_ANY   2
#define AMBIG_FULL  3

typedef struct {
    int  *colcoords;
    int   ncolcoord;
    int   colord;
    int   colwgh;
    int  *colids;
    int   ncolid;
    int   neg_ordering;
    int   separator;
} colinfo;

typedef struct {
    double x, y;
} ilinepoint2d;

typedef struct {
    double x[MAX_COORDMATCH_DIM];
} ilinepointnd;

typedef union {
    ilinepoint2d d;
    ilinepointnd n;
} ilinepoint;

typedef struct {
    ilinepoint p;
    double     ordseq, weight;
    char      *line, *id;
} iline;

typedef struct {
    int   nmiter, friter;
    double rejlevel;
    double maxdist, unitarity;
    int   parity;
    int   ttype, use_ordering, maxnum_ref, maxnum_inp;
    int   wcat, w_magnitude;
    double wpower;
    int   is_centering;
    double refcx, refcy;
    double inpcx, inpcy;
    double maxcenterdist;
    transformation *htf;
    int   hintorder;
} matchpointtune;

typedef struct {
    double wsigma;
    double nsigma;
    double unitarity;
    double time_total;
    double time_trimatch;
    double time_symmatch;
    int   nmiter;
    int   tri_level;
    double hull_coverage;
} matchpointstat;

extern int   is_verbose, is_comment;
extern char *progbasename;

// int  fprint_error(char *expr, ...);
// int  fprint_warning(char *expr, ...);

int  colinfo_reset(colinfo *col);
int  normalize_columns(colinfo *col, char *strcolcoord, char *strcolid);
int  read_match_data_points(FILE *fr, colinfo *col, iline **rils, int *rnil, int mtype);

int  do_pointmatch(iline *refls, int nref, iline *inpls, int ninp,
                   matchpointtune *mptp, cphit **rhits, int *rnhit,
                   int order, double **vfits, matchpointstat *mps);
int  do_coordmatch(iline *refls, int nref, iline *inpls, int ninp,
                   cphit **rhits, int *rnhit, int dim, double maxdist);
int  do_idmatch(iline *refls, int nref, iline *inpls, int ninp,
                cphit **rhits, int *rnhit, int ambig);

int  fprint_grmatch_usage(FILE *fw);
int  fprint_pointmatch_stat(FILE *fw, int nhit, int nref, int ninp, matchpointstat *mps);
int  fprint_general_stat(FILE *fw, int nhit, int nref, int ninp);
int  fwrite_matched_or_excluded_lines(FILE *fw, iline *refls, int nref,
                                      cphit *hits, int nhit, int type, int is_excluded);

#endif
