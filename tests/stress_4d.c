/**
 * @file      stress_4d.c
 * @brief     Randomised self-consistency sweep for the 4D kernel.
 *
 * For each random shape and cell, the fraction, the centroid and the
 * interface measure of ONE cell are compared with the sum over its sixteen
 * children. Quadrature error cancels to a large extent in that comparison,
 * so what it is really looking for is a MISSED TOPOLOGY CHANGE: a
 * subdivision the parent failed to place and a child, seeing a smaller
 * piece of the interface, did not need. Those show up as O(1) discrepancies,
 * far above the tolerance.
 *
 * The shapes are chosen to cover what the targeted suite does not: a slab
 * gives two interface sheets crossing one height line, a 4D torus gives a
 * cross-section that changes genus through the cell, and the cell aspect
 * ratios are drawn over three decades.
 *
 *   ./stress_4d [ncase] [seed] [-v]
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vofi.h"

/* a small reproducible generator, so a failure can always be replayed */
static unsigned long rng_s = 1;
static double rnd(void)
{
  rng_s = rng_s*6364136223846793005UL + 1442695040888963407UL;
  return (double)((rng_s >> 11) & 0x1FFFFFFFFFFFFFUL)/9007199254740992.0;
}
static double rndr(double a, double b) { return a + (b - a)*rnd(); }

/* ------------------------------------------------------------------ */
static double f_ball(const double x[], void *p)
{
  const double *q = (const double *) p;
  double d = 0., t;
  int i;
  for (i = 0; i < 4; i++) { t = x[i] - q[i]; d += t*t; }
  return d - q[4]*q[4];
}

static double f_ellip(const double x[], void *p)
{
  const double *q = (const double *) p;
  double d = 0., t;
  int i;
  for (i = 0; i < 4; i++) { t = (x[i] - q[i])/q[4+i]; d += t*t; }
  return d - 1.;
}

static double f_plane(const double x[], void *p)
{
  const double *q = (const double *) p;
  return q[0]*x[0] + q[1]*x[1] + q[2]*x[2] + q[3]*x[3] - q[4];
}

/* two parallel sheets: a height line along any direction crosses twice, so
   the cell carries two cut sectors at once                              */
static double f_slab(const double x[], void *p)
{
  const double *q = (const double *) p;
  double s = q[0]*x[0] + q[1]*x[1] + q[2]*x[2] + q[3]*x[3];
  return fabs(s - q[4]) - q[5];
}

/* a 4D torus: the cross-section opens a hole part-way through the cell */
static double f_torus(const double x[], void *p)
{
  const double *q = (const double *) p;
  double a = x[0]-q[0], b = x[1]-q[1], c = x[2]-q[2], d = x[3]-q[3];
  double r = sqrt(a*a + b*b) - q[4];
  return r*r + c*c + d*d - q[5]*q[5];
}

/* a sphere whose radius varies with the fourth coordinate */
static double f_grow(const double x[], void *p)
{
  const double *q = (const double *) p;
  double d = 0., t, R;
  int i;
  for (i = 0; i < 3; i++) { t = x[i] - q[i]; d += t*t; }
  R = q[3] + q[4]*x[3];
  return d - R*R;
}

typedef double (*fptr)(const double [], void *);

/* ------------------------------------------------------------------ */
int main(int argc, char **argv)
{
  int  nex[2] = {1,1}, npt[6] = {0,0,0,0,0,0}, nvis[2] = {0,0};
  int  ncase = 400, verbose = 0, i, v, ish, nbad = 0, ncut = 0, a;
  double x0[4], h[4], hs[4], xs[4], par[8];
  double xex[8], xgam[4], xe2[8], xg2[4];
  double cc, sub, msub, cs[4], ms[4], w, worstv = 0., worstm = 0.;
  double worstc = 0., d;
  fptr f;

  if (argc > 1) ncase = atoi(argv[1]);
  if (argc > 2) rng_s = (unsigned long) atoi(argv[2]);
  if (argc > 3) verbose = 1;

  for (i = 0; i < ncase; i++) {
    double hmax = 0., pt[4], nv[4], nn, sc;

    /* an anisotropic cell, aspect ratios over three decades */
    for (a = 0; a < 4; a++) {
      h[a]  = pow(10., rndr(-1.5,0.));
      x0[a] = rndr(-0.6,0.6);
      if (h[a] > hmax) hmax = h[a];
    }
    /* a point inside the cell and a direction: every shape below is built
       so that its interface passes through that point, and so that its
       radius of curvature is several cell diagonals. An interface the
       grid does NOT resolve is outside what the method claims -- a cell
       and its children then legitimately disagree -- so those are not
       generated here.                                                   */
    nn = 0.;
    for (a = 0; a < 4; a++) {
      pt[a] = x0[a] + rndr(0.15,0.85)*h[a];
      nv[a] = rndr(-1.,1.);
      nn += nv[a]*nv[a];
    }
    nn = sqrt(nn) + 1.e-30;
    for (a = 0; a < 4; a++) nv[a] /= nn;

    ish = (int)(rnd()*5.);
    memset(par,0,sizeof(par));
    switch (ish) {
    case 0:                                                     /* ball */
      f = f_ball;
      par[4] = rndr(3.,12.)*hmax;
      for (a = 0; a < 4; a++) par[a] = pt[a] - par[4]*nv[a];
      break;
    case 1:                                                 /* ellipsoid */
      f = f_ellip;
      for (a = 0; a < 4; a++) par[4+a] = rndr(4.,12.)*hmax;
      for (a = 0; a < 4; a++)
        par[a] = pt[a] - 0.7*par[4+a]*nv[a];
      sc = 0.;                       /* rescale so it passes through pt */
      for (a = 0; a < 4; a++) {
        double t = (pt[a] - par[a])/par[4+a];
        sc += t*t;
      }
      sc = sqrt(sc);
      for (a = 0; a < 4; a++) par[4+a] *= sc;
      break;
    case 2:                                                    /* plane */
      f = f_plane;
      par[4] = 0.;
      for (a = 0; a < 4; a++) {
        par[a] = nv[a];
        par[4] += nv[a]*pt[a];
      }
      break;
    case 3:                                                     /* slab */
      f = f_slab;
      par[5] = rndr(1.5,6.)*hmax;
      par[4] = 0.;
      for (a = 0; a < 4; a++) {
        par[a] = nv[a];
        par[4] += nv[a]*pt[a];
      }
      par[4] += (rnd() < 0.5 ? par[5] : -par[5]);
      break;
    default:                                        /* torus or growing */
      if (rnd() < 0.5) {
        double q0, q1;
        f = f_torus;
        for (a = 0; a < 4; a++) par[a] = pt[a] - rndr(3.,9.)*hmax*nv[a];
        par[4] = rndr(4.,12.)*hmax;
        q0 = sqrt((pt[0]-par[0])*(pt[0]-par[0]) +
                  (pt[1]-par[1])*(pt[1]-par[1])) - par[4];
        q1 = q0*q0 + (pt[2]-par[2])*(pt[2]-par[2])
                   + (pt[3]-par[3])*(pt[3]-par[3]);
        par[5] = sqrt(q1);            /* the tube passes through pt */
        if (par[5] < 3.*hmax) par[5] = 3.*hmax;
      }
      else {
        double d3 = 0.;
        f = f_grow;
        for (a = 0; a < 3; a++) {
          par[a] = pt[a] - rndr(3.,10.)*hmax*nv[a];
          d3 += (pt[a]-par[a])*(pt[a]-par[a]);
        }
        par[4] = rndr(-0.6,0.6);
        par[3] = sqrt(d3) - par[4]*pt[3];     /* R(t) passes through pt */
        if (par[3] < 3.*hmax) par[3] = 3.*hmax;
      }
      break;
    }

    memset(xex,0,sizeof(xex));
    memset(xgam,0,sizeof(xgam));
    cc = vofi_get_cc_gam(f,par,x0,h,xex,xgam,nex,npt,nvis,4);
    if (cc <= 0. || cc >= 1.)
      continue;                       /* only cut cells say anything */
    ncut++;

    for (a = 0; a < 4; a++) hs[a] = 0.5*h[a];
    sub = msub = 0.;
    for (a = 0; a < 4; a++) cs[a] = ms[a] = 0.;
    for (v = 0; v < 16; v++) {
      for (a = 0; a < 4; a++)
        xs[a] = x0[a] + ((v >> a) & 1)*hs[a];
      memset(xe2,0,sizeof(xe2));
      memset(xg2,0,sizeof(xg2));
      w = vofi_get_cc_gam(f,par,xs,hs,xe2,xg2,nex,npt,nvis,4);
      sub  += w;
      msub += xe2[4];
      for (a = 0; a < 4; a++) {
        cs[a] += w*xe2[a];
        ms[a] += xe2[4]*xg2[a];
      }
    }
    sub /= 16.;

    d = fabs(cc - sub);
    if (d > worstv) worstv = d;
    if (d > 1.e-5) nbad++;
    if (msub > 0.) {
      double dm = fabs(xex[4] - msub)/msub;
      if (dm > worstm) worstm = dm;
      if (dm > 1.e-3) nbad++;
    }
    if (sub > 0.)
      for (a = 0; a < 4; a++) {
        double dc = fabs(xex[a] - cs[a]/(16.*sub))/h[a];
        if (dc > worstc) worstc = dc;
        if (dc > 1.e-4) nbad++;
      }
    if (verbose && (d > 1.e-5 || (msub > 0. && fabs(xex[4]-msub)/msub > 1.e-3)))
      printf("  case %d shape %d cc %.12g sub %.12g  meas %.10g sub %.10g\n"
             "    x0 %.6g %.6g %.6g %.6g  h %.6g %.6g %.6g %.6g\n",
             i, ish, cc, sub, xex[4], msub,
             x0[0],x0[1],x0[2],x0[3], h[0],h[1],h[2],h[3]),
      printf("    par %.9g %.9g %.9g %.9g %.9g %.9g\n",
             par[0],par[1],par[2],par[3],par[4],par[5]);
  }

  printf("stress_4d: %d cases, %d cut, %d suspect\n", ncase, ncut, nbad);
  printf("  worst |cc - children|          %.3e  (tol 1e-5)\n", worstv);
  printf("  worst relative measure gap     %.3e  (tol 1e-3)\n", worstm);
  printf("  worst centroid gap / cell edge %.3e  (tol 1e-4)\n", worstc);

  return nbad ? 1 : 0;
}
