/**
 * @file      sweep.h
 * @brief     The deterministic cell sweep shared by the golden-data
 *            generator and the regression test.
 */

#ifndef VOFI_TEST_SWEEP_H
#define VOFI_TEST_SWEEP_H

#include "cases.h"
#include "vofi.h"

/* run one case, calling sink() once per cell with VOFI_NREC doubles */
static void vofi_run_case(const vofi_case *c,
                          void (*sink)(const double *, void *), void *ctx)
{
  int    nex[2] = {1,1}, npt[6] = {0,0,0,0,0,0}, nvis[2] = {0,0};
  int    i, j, k, l, n = c->n, nd = c->ndim;
  double x0[4], h0[4], xex[8], xgam[4], rec[VOFI_NREC];
  double par[9];

  memcpy(par, c->par, sizeof(par));
  for (i = 0; i < 4; i++)
    h0[i] = c->hcell[i];

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < (nd >= 3 ? n : 1); k++)
        for (l = 0; l < (nd >= 4 ? n : 1); l++) {
          x0[0] = c->xlo[0] + i*c->hcell[0];
          x0[1] = c->xlo[1] + j*c->hcell[1];
          x0[2] = (nd >= 3) ? c->xlo[2] + k*c->hcell[2] : 0.;
          x0[3] = (nd >= 4) ? c->xlo[3] + l*c->hcell[3] : 0.;
          memset(xex,  0, sizeof(xex));
          memset(xgam, 0, sizeof(xgam));
          rec[0] = vofi_get_cc_gam(c->func, par, x0, h0, xex, xgam,
                                   nex, npt, nvis, nd);
          rec[1] = xex[0];  rec[2] = xex[1];  rec[3] = xex[2];
          rec[4] = xex[3];  rec[5] = (nd == 4) ? xex[4] : 0.;
          rec[6] = xgam[0]; rec[7] = xgam[1]; rec[8] = xgam[2];
          rec[9] = (nd == 4) ? xgam[3] : 0.;
          rec[10] = (double) vofi_get_cell_type(c->func, par, x0, h0, nd);
          sink(rec, ctx);
        }
}

#endif
