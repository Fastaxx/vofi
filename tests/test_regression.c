/**
 * @file      test_regression.c
 * @brief     Compare the current build against tests/golden.dat.
 *
 * The 2D and 3D paths must be untouched by the 4D work. Tolerance is tight
 * on purpose: the dimension refactor is meant to be a no-op, and the only
 * differences it can legitimately produce are in the finite stand-ins for
 * "no box constraint along this direction" inside the minimisers.
 *
 *   ./test_regression tests/golden.dat
 */

#include "sweep.h"

#define TOL 1.0e-12

typedef struct {
  FILE  *fp;
  int    ncell, nbad;
  double worst;
  const char *cname;
} ctx_t;

static int read_row(FILE *fp, double *v)
{
  int  i;
  char line[4096];

  for (;;) {
    if (!fgets(line, sizeof(line), fp)) return 0;
    if (line[0] == '#' || line[0] == '\n') continue;
    break;
  }
  {
    char *p = line, *e;
    for (i = 0; i < VOFI_NREC; i++) {
      v[i] = strtod(p, &e);
      if (p == e) return 0;
      p = e;
    }
  }
  return 1;
}

static void check(const double *rec, void *vctx)
{
  ctx_t *c = (ctx_t *) vctx;
  double ref[VOFI_NREC], d;
  int    i;

  if (!read_row(c->fp, ref)) {
    fprintf(stderr, "golden file exhausted at cell %d of case %s\n",
            c->ncell, c->cname);
    exit(2);
  }
  for (i = 0; i < VOFI_NREC; i++) {
    d = fabs(rec[i] - ref[i]);
    if (d > c->worst) c->worst = d;
    if (!(d <= TOL)) {
      if (c->nbad < 10)
        fprintf(stderr, "MISMATCH %s cell %d field %d: got %.17g want %.17g\n",
                c->cname, c->ncell, i, rec[i], ref[i]);
      c->nbad++;
    }
  }
  c->ncell++;
}

int main(int argc, char **argv)
{
  ctx_t c;
  int   i, nbad = 0;

  if (argc < 2) {
    fprintf(stderr, "usage: %s golden.dat\n", argv[0]);
    return 2;
  }
  c.fp = fopen(argv[1], "r");
  if (!c.fp) {
    fprintf(stderr, "cannot open %s\n", argv[1]);
    return 2;
  }
  c.worst = 0.;
  for (i = 0; i < VOFI_NCASE_23 + VOFI_NCASE_4; i++) {
    const vofi_case *cs = (i < VOFI_NCASE_23) ? &vofi_cases_23[i]
                                              : &vofi_cases_4[i-VOFI_NCASE_23];
    c.ncell = c.nbad = 0;
    c.cname = cs->name;
    vofi_run_case(cs, check, &c);
    printf("%-14s %5d cells  %s\n", c.cname, c.ncell,
           c.nbad ? "FAIL" : "ok");
    nbad += c.nbad;
  }
  fclose(c.fp);
  printf("worst absolute deviation: %.3e  (tol %.1e)\n", c.worst, TOL);
  return nbad ? 1 : 0;
}
