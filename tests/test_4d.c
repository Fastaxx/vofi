/**
 * @file      test_4d.c
 * @brief     Validation of the 4D kernel.
 *
 * The tests are ordered by how much they can prove:
 *
 *  1. HALF-SPACE. A hyperplane cut of a 4D box has a closed form, and the
 *     whole method is exact for it -- the heights are linear, the sector
 *     boundaries are exactly the eight edge crossings, and Gauss-Legendre
 *     integrates the piecewise-cubic V(u) exactly. Anything above round-off
 *     here is a structural defect, not an accuracy shortfall.
 *  2. HYPERBALL. No closed form per cell, so the grid sum is checked
 *     against pi^2 R^4 / 2.
 *  3. SPACE-TIME SLAB. A 3D ball whose radius varies with the fourth
 *     coordinate: the motivating case, and the only one whose cross-section
 *     is not symmetric in u.
 *  4. REFINEMENT CONSISTENCY. One cell against its sixteen children, for
 *     both the fraction and the centroid.
 *  5. DEGENERATE CONFIGURATIONS. Full, empty, tangent, vertex-touching,
 *     flat cells, and each of the four axis-aligned orientations.
 *  6. CELL TYPE. vofi_get_cell_type must not contradict vofi_get_cc.
 *  7. 3-FACE CAP. The one configuration that needs a genuine minimisation
 *     inside a cubic 3-face, checked against the exact spherical cap.
 *  8. ARRAY-SIZE CONTRACT. Exactly-sized caller arrays, so that a run
 *     under AddressSanitizer proves the library reads no further.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vofi.h"

#define MAXD(a,b) ((a) > (b) ? (a) : (b))

static int nfail = 0, ntest = 0;

static void ck(const char *what, double got, double want, double tol)
{
  double d = fabs(got - want);
  ntest++;
  if (!(d <= tol)) {
    nfail++;
    printf("  FAIL %-38s got %.16g want %.16g  |d|=%.3e > %.1e\n",
           what, got, want, d, tol);
  }
  else
    printf("  ok   %-38s %.16g  (|d|=%.2e)\n", what, got, d);
}

/* ------------------------------------------------------------------ */
/* implicit functions                                                  */

static double f_plane(const double x[], void *p)
{
  const double *q = (const double *) p;
  return q[0]*x[0] + q[1]*x[1] + q[2]*x[2] + q[3]*x[3] - q[4];
}

static double f_ball(const double x[], void *p)
{
  const double *q = (const double *) p;
  double d = 0., t;
  int i;
  for (i = 0; i < 4; i++) { t = x[i] - q[i]; d += t*t; }
  return d - q[4]*q[4];
}

/* a 3D ball whose radius grows with the fourth coordinate */
static double f_slab(const double x[], void *p)
{
  const double *q = (const double *) p;
  double d = 0., t, R;
  int i;
  for (i = 0; i < 3; i++) { t = x[i] - q[i]; d += t*t; }
  R = q[3] + q[4]*x[3];
  return d - R*R;
}

/* ------------------------------------------------------------------ */
/* exact fraction of the box [x0, x0+h] lying in the half-space n.x <= d
   (inclusion-exclusion over the corners; every n[i]*h[i] must be non-zero) */
static double plane_frac(const double n[4], double d, const double x0[4],
                         const double h[4])
{
  double a[4], c = d, prod = 1., s = 0., t;
  int i, S, bits;

  for (i = 0; i < 4; i++) {
    a[i] = n[i]*h[i];
    c   -= n[i]*x0[i];
  }
  for (i = 0; i < 4; i++)                     /* y_i -> 1 - y_i if a_i < 0 */
    if (a[i] < 0.) { c -= a[i]; a[i] = -a[i]; }
  for (i = 0; i < 4; i++) prod *= a[i];

  for (S = 0; S < 16; S++) {
    t = c;
    bits = 0;
    for (i = 0; i < 4; i++)
      if (S & (1 << i)) { t -= a[i]; bits++; }
    if (t > 0.)
      s += ((bits & 1) ? -1. : 1.)*t*t*t*t;
  }
  return s/(24.*prod);
}

/* 3-volume of the interface a hyperplane cuts in the box. The coarea
   identity for a linear f gives H^3 = |n| * d(Volume)/d(offset), and the
   volume is the closed form above, so this is exact too.               */
static double plane_meas(const double n[4], double d, const double x0[4],
                         const double h[4])
{
  double a[4], c = d, prod = 1., s = 0., t, nn = 0., vol = 1.;
  int i, S, bits;

  for (i = 0; i < 4; i++) {
    a[i] = n[i]*h[i];
    c   -= n[i]*x0[i];
    nn  += n[i]*n[i];
    vol *= h[i];
  }
  for (i = 0; i < 4; i++)
    if (a[i] < 0.) { c -= a[i]; a[i] = -a[i]; }
  for (i = 0; i < 4; i++) prod *= a[i];

  for (S = 0; S < 16; S++) {
    t = c;
    bits = 0;
    for (i = 0; i < 4; i++)
      if (S & (1 << i)) { t -= a[i]; bits++; }
    if (t > 0.)
      s += ((bits & 1) ? -1. : 1.)*t*t*t;
  }
  return sqrt(nn)*vol*s/(6.*prod);
}

/* ------------------------------------------------------------------ */
static double cc4(double (*f)(const double [], void *), void *par,
                  const double x0[], const double h[], double xex[])
{
  int nex[2] = {1,0}, npt[6] = {0,0,0,0,0,0}, nvis[2] = {0,0};
  double loc[8];
  double *px = xex ? xex : loc;

  memset(px, 0, 8*sizeof(double));

  return vofi_get_cc(f, par, x0, h, px, nex, npt, nvis, 4);
}

/* the same with the interface measure and its centroid */
static double cg4(double (*f)(const double [], void *), void *par,
                  const double x0[], const double h[], double xex[],
                  double xgam[])
{
  int nex[2] = {1,1}, npt[6] = {0,0,0,0,0,0}, nvis[2] = {0,0};

  memset(xex,  0, 8*sizeof(double));
  memset(xgam, 0, 4*sizeof(double));

  return vofi_get_cc_gam(f, par, x0, h, xex, xgam, nex, npt, nvis, 4);
}

/* ------------------------------------------------------------------ */
static void test_halfspace(void)
{
  /* fixed normals with every component away from zero, so the closed form
     stays well conditioned; offsets sweep the whole range of n.x          */
  static const double nn[6][4] = {
    { 0.6, 0.8, 0.5, 0.4}, {-0.7, 0.3, 0.9,-0.5}, { 1.0,-0.6, 0.4, 0.8},
    { 0.35,0.95,-0.75,0.55}, {-0.4,-0.9, 0.6,-0.3}, { 0.5, 0.5, 0.5, 0.5}
  };
  double x0[4] = {0.13,-0.27, 0.41, 0.05};
  double h[4]  = {0.7, 1.3, 0.9, 1.1};
  double par[5], cc, ex, worst = 0., d;
  int i, k, m, cut = 0;
  char lbl[64];

  printf("1. half-space cuts (the method is exact here)\n");
  for (i = 0; i < 6; i++)
    for (k = 0; k < 15; k++) {
      double lo = 0., hi = 0., off, a, b;
      for (m = 0; m < 4; m++) {
        a = nn[i][m]*x0[m];
        b = nn[i][m]*(x0[m] + h[m]);
        lo += (a < b) ? a : b;
        hi += (a < b) ? b : a;
      }
      off = lo + (hi - lo)*(k + 0.5)/15.;
      memcpy(par, nn[i], 4*sizeof(double));
      par[4] = off;
      cc = cc4(f_plane, par, x0, h, NULL);
      ex = plane_frac(nn[i], off, x0, h);
      d = fabs(cc - ex);
      if (d > worst) worst = d;
      if (ex > 1.e-6 && ex < 1.-1.e-6) cut++;
    }
  snprintf(lbl,sizeof(lbl), "worst error over 90 cuts (%d cut)", cut);
  ck(lbl, worst, 0., 1.e-12);
}

/* ------------------------------------------------------------------ */
static void test_axis_aligned(void)
{
  double x0[4] = {-0.2, 0.3, 1.1, -0.7};
  double h[4]  = {0.5, 0.8, 0.4, 0.6};
  double par[5], xex[8], cc, ex, worst = 0., wcen = 0., d;
  int a, b, k;

  printf("2. axis-aligned interfaces (all four orientations)\n");
  for (a = 0; a < 4; a++)
    for (k = 1; k < 8; k++) {
      memset(par, 0, sizeof(par));
      par[a] = 1.;
      par[4] = x0[a] + h[a]*k/8.;
      cc = cc4(f_plane, par, x0, h, xex);
      ex = k/8.;
      if (fabs(cc - ex) > worst) worst = fabs(cc - ex);
      /* the wet region is a box: its centroid is known exactly */
      for (b = 0; b < 4; b++) {
        d = fabs(xex[b] - (x0[b] + ((b == a) ? 0.5*ex : 0.5)*h[b]));
        if (d > wcen) wcen = d;
      }
    }
  ck("worst error, 4 axes x 7 offsets", worst, 0., 1.e-13);
  ck("worst centroid error (exact box)", wcen, 0., 1.e-13);
}

/* ------------------------------------------------------------------ */
static double ball_sum(int N, double L, double R, int *ncut)
{
  double x0[4], h[4], par[5] = {0.,0.,0.,0., 0.}, tot = 0., cc;
  int i, j, k, l;

  par[4] = R;
  for (i = 0; i < 4; i++) h[i] = L/N;
  *ncut = 0;
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      for (k = 0; k < N; k++)
        for (l = 0; l < N; l++) {
          x0[0] = -0.5*L + i*h[0]; x0[1] = -0.5*L + j*h[1];
          x0[2] = -0.5*L + k*h[2]; x0[3] = -0.5*L + l*h[3];
          cc = cc4(f_ball, par, x0, h, NULL);
          if (cc > 0. && cc < 1.) (*ncut)++;
          tot += cc*h[0]*h[1]*h[2]*h[3];
        }
  return tot;
}

static void test_hyperball(void)
{
  double R = 0.8, exact = 0.5*M_PI*M_PI*R*R*R*R, v4, v6;
  int c4, c6;
  char lbl[64];

  printf("3. hyperball: the grid sum must be pi^2 R^4 / 2\n");
  v4 = ball_sum(4, 2.0, R, &c4);
  v6 = ball_sum(6, 2.0, R, &c6);
  snprintf(lbl,sizeof(lbl), "4^4 grid (%d cut cells)", c4);
  ck(lbl, v4/exact, 1., 1.e-6);
  snprintf(lbl,sizeof(lbl), "6^4 grid (%d cut cells)", c6);
  ck(lbl, v6/exact, 1., 1.e-7);
}

/* ------------------------------------------------------------------ */
static void test_slab(void)
{
  /* |x|^2 = R(t)^2 with R(t) = R0 + g t, t in [0,1]: the 4-volume is the
     time integral of the ball volume, pi ((R0+g)^4 - R0^4)/(3 g)        */
  double par[5] = {0.,0.,0., 0.55, 0.15};
  double R0 = 0.55, g = 0.15, exact;
  double x0[4], h[4], tot = 0., cc;
  int i, j, k, l, N = 5, M = 4, ncut = 0;
  char lbl[64];

  printf("4. space-time slab: a 3D ball growing along the 4th axis\n");
  exact = M_PI*(pow(R0+g,4.) - pow(R0,4.))/(3.*g);
  for (i = 0; i < 3; i++) h[i] = 2.0/N;
  h[3] = 1.0/M;
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      for (k = 0; k < N; k++)
        for (l = 0; l < M; l++) {
          x0[0] = -1. + i*h[0]; x0[1] = -1. + j*h[1];
          x0[2] = -1. + k*h[2]; x0[3] = l*h[3];
          cc = cc4(f_slab, par, x0, h, NULL);
          if (cc > 0. && cc < 1.) ncut++;
          tot += cc*h[0]*h[1]*h[2]*h[3];
        }
  snprintf(lbl,sizeof(lbl), "5^3 x 4 grid (%d cut cells)", ncut);
  ck(lbl, tot/exact, 1., 1.e-6);
}

/* ------------------------------------------------------------------ */
static void test_refinement(void)
{
  double x0[4] = {0.15, -0.35, -0.05, 0.10};
  double h[4]  = {0.6, 0.5, 0.7, 0.55};
  double hs[4], xs[4], xex[8], xexs[8];
  double cc, sub, csub[4], w, wtot, par[5];
  int i, v, ic;
  char lbl[64];

  printf("5. refinement consistency: one cell against its 16 children\n");
  for (i = 0; i < 4; i++) hs[i] = 0.5*h[i];

  for (ic = 0; ic < 2; ic++) {
    double (*fn)(const double [], void *) = ic ? f_plane : f_ball;
    if (ic) {
      par[0] = 0.62; par[1] = -0.41; par[2] = 0.77; par[3] = 0.35;
      par[4] = 0.62*0.45 - 0.41*(-0.10) + 0.77*0.30 + 0.35*0.37;
    }
    else {
      par[0] = 0.31; par[1] = -0.12; par[2] = 0.07; par[3] = 0.22;
      par[4] = 0.83;
    }
    cc = cc4(fn, par, x0, h, xex);
    sub = wtot = 0.;
    for (i = 0; i < 4; i++) csub[i] = 0.;
    for (v = 0; v < 16; v++) {
      for (i = 0; i < 4; i++)
        xs[i] = x0[i] + ((v >> i) & 1)*hs[i];
      w = cc4(fn, par, xs, hs, xexs);
      sub += w;
      if (w > 0.)
        for (i = 0; i < 4; i++) csub[i] += w*xexs[i];
      wtot += w;
    }
    sub /= 16.;
    for (i = 0; i < 4; i++)
      csub[i] = (wtot > 0.) ? csub[i]/wtot : 0.;
    snprintf(lbl,sizeof(lbl), "%s fraction", ic ? "plane" : "ball ");
    ck(lbl, cc, sub, ic ? 1.e-13 : 2.e-6);
    for (i = 0; i < 4; i++) {
      snprintf(lbl,sizeof(lbl), "%s centroid[%d]", ic ? "plane" : "ball ", i);
      ck(lbl, xex[i], csub[i], ic ? 1.e-12 : 2.e-6);
    }
  }
}

/* ------------------------------------------------------------------ */
static void test_degenerate(void)
{
  double x0[4] = {0.,0.,0.,0.}, h[4] = {1.,1.,1.,1.};
  double par[5], xex[8], cc;

  printf("6. degenerate configurations\n");

  par[0]=par[1]=par[2]=par[3]=0.5; par[4]=100.;
  ck("ball swallowing the cell (full)", cc4(f_ball,par,x0,h,NULL), 1., 0.);
  par[4]=0.01;
  ck("unresolved ball (empty, as in 2D/3D)",
     cc4(f_ball,par,x0,h,NULL), 0., 0.);

  memset(par,0,sizeof(par)); par[0]=1.; par[4]=0.;
  ck("hyperplane on the lower 3-face", cc4(f_plane,par,x0,h,NULL), 0., 0.);
  par[4]=1.;
  ck("hyperplane on the upper 3-face", cc4(f_plane,par,x0,h,NULL), 1., 0.);

  par[0]=par[1]=par[2]=0.5; par[3]=-0.5; par[4]=0.5;
  ck("ball tangent to a 3-face from outside",
     cc4(f_ball,par,x0,h,NULL), 0., 1.e-12);
  par[0]=par[1]=par[2]=par[3]=-0.5; par[4]=1.0;
  ck("ball touching a single vertex", cc4(f_ball,par,x0,h,NULL), 0., 1.e-9);

  par[0]=par[1]=par[2]=par[3]=0.5; par[4]=1.0;
  ck("diagonal hyperplane through centre",
     cc4(f_plane,par,x0,h,xex), 0.5, 1.e-14);
  /* the centroid of {x in [0,1]^4 : sum x <= 2} is (23/60,...) exactly */
  ck("  its centroid (exact, 23/60)", xex[0], 23./60., 1.e-13);
  ck("  its centroid, x0 against x3", xex[0], xex[3], 1.e-13);

  /* a very flat cell: the interface has to stay resolved in the three
     long directions, so tilt the hyperplane rather than shrink a ball  */
  h[3] = 1.e-3;
  par[0]=0.6; par[1]=0.5; par[2]=0.4; par[3]=0.7;
  par[4]=0.6*0.5 + 0.5*0.5 + 0.4*0.5 + 0.7*0.5e-3;
  cc = cc4(f_plane,par,x0,h,xex);
  ck("flat cell (h3/h0 = 1e-3), hyperplane",
     cc, plane_frac(par,par[4],x0,h), 1.e-12);
  par[0]=par[1]=par[2]=0.2; par[3]=0.5e-3; par[4]=0.4;
  cc = cc4(f_ball,par,x0,h,NULL);
  ck("flat cell, ball cutting it",
     (cc > 0. && cc < 1.) ? 1. : 0., 1., 0.);
  h[3] = 1.;
}

/* ------------------------------------------------------------------ */
static void test_cell_type(void)
{
  double x0[4], h[4] = {0.4,0.4,0.4,0.4}, par[5] = {0.,0.,0.,0., 0.7};
  double cc;
  int i, j, k, l, N = 5, bad = 0, ityp;

  printf("7. cell type against the fraction\n");
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      for (k = 0; k < N; k++)
        for (l = 0; l < N; l++) {
          x0[0] = -1. + i*h[0]; x0[1] = -1. + j*h[1];
          x0[2] = -1. + k*h[2]; x0[3] = -1. + l*h[3];
          ityp = vofi_get_cell_type(f_ball, par, x0, h, 4);
          cc = cc4(f_ball, par, x0, h, NULL);
          if (ityp == 1 && cc != 1.) bad++;
          if (ityp == 0 && cc != 0.) bad++;
        }
  ck("full/empty verdicts contradicting cc", (double) bad, 0., 0.);
}

/* ------------------------------------------------------------------ */
/* A cap poking into the cell through the INTERIOR of a cubic 3-face: no
   vertex changes sign, and the cap touches none of that 3-face's own
   square faces either. This is the one configuration that needs a genuine
   3D minimisation inside a 3-face -- vofi_check_cell_consistency plus
   vofi_get_cell_min -- and it is missed by anything that only inspects the
   2-faces. The reference is the 4D spherical cap
       V = int_d^R (4/3) pi (R^2 - z^2)^{3/2} dz,
   quadrature of the exact integrand to 15 digits.                      */
static void test_face_cap(void)
{
  double x0[4] = {0.,0.,0.,0.}, h[4] = {1.,1.,1.,1.};
  double par[5] = {0.5,0.5,0.5,-0.2, 0.35};
  double exact = 0.00667223772686801, cc;
  int a, ityp;
  char lbl[64];

  printf("8. cap through the interior of a 3-face\n");
  for (a = 0; a < 4; a++) {         /* the same cap on each of the 4 axes */
    memset(par, 0, sizeof(par));
    par[0] = par[1] = par[2] = par[3] = 0.5;
    par[a] = -0.2;
    par[4] = 0.35;
    ityp = vofi_get_cell_type(f_ball, par, x0, h, 4);
    cc = cc4(f_ball, par, x0, h, NULL);
    snprintf(lbl,sizeof(lbl), "cap through the axis-%d face", a);
    ck(lbl, cc, exact, 5.e-6);
    snprintf(lbl,sizeof(lbl), "  reported cut, axis %d", a);
    ck(lbl, (double) ityp, -1., 0.);
  }
}

/* ------------------------------------------------------------------ */
/* the array-size contract of vofi.h. Run under AddressSanitizer this is
   the check that the library never reads past the caller's arrays: the
   arrays here are sized EXACTLY as documented, not padded.            */
static void test_api(void)
{
  int nex[2] = {1,1}, npt4[4] = {0,0,0,0}, npt6[6] = {0,0,0,0,0,0};
  int nvis[2] = {0,0};
  double x2[2] = {0.1,0.2}, h2[2] = {0.4,0.5};
  double x3[3] = {0.1,0.2,0.3}, h3[3] = {0.4,0.5,0.6};
  double x4[4] = {0.1,0.2,0.3,0.4}, h4[4] = {0.4,0.5,0.6,0.7};
  double xex4[4], xex5[5], xgam3[3], xgam4[4];
  double par[5] = {0.6,0.5,0.4,0.7, 0.}, cc;
  int i;

  printf("9. array-size contract (exactly sized caller arrays)\n");

  par[4] = 0.6*0.3 + 0.5*0.45;                              /* 2D */
  cc = vofi_get_cc_gam(f_plane,par,x2,h2,xex4,xgam3,nex,npt4,nvis,2);
  ck("2D: 2-element xin/h0, 4-element xex", (cc > 0. && cc < 1.) ? 1. : 0.,
     1., 0.);

  par[4] = 0.6*0.3 + 0.5*0.45 + 0.4*0.6;                    /* 3D */
  cc = vofi_get_cc_gam(f_plane,par,x3,h3,xex4,xgam3,nex,npt4,nvis,3);
  ck("3D: 3-element xin/h0, 4-element xex", (cc > 0. && cc < 1.) ? 1. : 0.,
     1., 0.);

  par[4] = 0.6*0.3 + 0.5*0.45 + 0.4*0.6 + 0.7*0.75;         /* 4D */
  for (i = 0; i < 5; i++) xex5[i] = -1.;
  for (i = 0; i < 4; i++) xgam4[i] = -1.;
  cc = vofi_get_cc_gam(f_plane,par,x4,h4,xex5,xgam4,nex,npt6,nvis,4);
  ck("4D: 4-element xin/h0, 5-element xex", (cc > 0. && cc < 1.) ? 1. : 0.,
     1., 0.);
  ck("4D: interface measure in xex[4]", xex5[4],
     plane_meas(par,par[4],x4,h4), 1.e-12);
  ck("4D: xgam inside the cell", (xgam4[3] > x4[3] &&
     xgam4[3] < x4[3] + h4[3]) ? 1. : 0., 1., 0.);
}

/* ------------------------------------------------------------------ */
/* The interface measure: the 3-volume of the hypersurface, and its
   centroid. For a hyperplane the integrand |grad f| / |df/dp| is constant
   and the projected domain is integrated exactly, so the measure is exact
   too. For a curved interface the integrand is unbounded where the
   interface turns over in p, which caps the convergence rate.          */
static void test_measure(void)
{
  static const double nn[4][4] = {
    { 0.6, 0.8, 0.5, 0.4}, {-0.7, 0.3, 0.9,-0.5},
    { 0.35,0.95,-0.75,0.55}, { 0.5, 0.5, 0.5, 0.5}
  };
  double x0[4] = {0.13,-0.27, 0.41, 0.05};
  double h[4]  = {0.7, 1.3, 0.9, 1.1};
  double xex[8], xgam[4], par[5], ex, worst = 0., d;
  double xc[4] = {0.,0.,0.,0.}, hc[4] = {1.,1.,1.,1.};
  double x1[4], hb[4], tot, R = 0.8;
  int i, k, m, N, j, l;

  printf("10. interface measure and its centroid\n");

  for (i = 0; i < 4; i++)
    for (k = 0; k < 9; k++) {
      double lo = 0., hi = 0., off, a, b;
      for (m = 0; m < 4; m++) {
        a = nn[i][m]*x0[m];
        b = nn[i][m]*(x0[m] + h[m]);
        lo += (a < b) ? a : b;
        hi += (a < b) ? b : a;
      }
      off = lo + (hi - lo)*(k + 0.5)/9.;
      memcpy(par, nn[i], 4*sizeof(double));
      par[4] = off;
      cg4(f_plane, par, x0, h, xex, xgam);
      ex = plane_meas(nn[i], off, x0, h);
      d = fabs(xex[4] - ex);
      if (ex > 0.) d = d/ex;
      if (d > worst) worst = d;
    }
  ck("hyperplane, worst relative error", worst, 0., 1.e-11);

  /* the diagonal section of the unit hypercube: 3-volume 4/3, centred */
  par[0]=par[1]=par[2]=par[3]=0.5; par[4]=1.0;
  cg4(f_plane, par, xc, hc, xex, xgam);
  ck("diagonal section, measure (exact 4/3)", xex[4], 4./3., 1.e-12);
  worst = 0.;
  for (i = 0; i < 4; i++)
    worst = MAXD(worst, fabs(xgam[i] - 0.5));
  ck("diagonal section, centroid at centre", worst, 0., 1.e-12);

  /* the 3-sphere: the grid sum must be 2 pi^2 R^3 */
  N = 6;
  for (i = 0; i < 4; i++) hb[i] = 2.0/N;
  memset(par, 0, sizeof(par));
  par[4] = R;
  tot = 0.;
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      for (k = 0; k < N; k++)
        for (l = 0; l < N; l++) {
          x1[0] = -1. + i*hb[0]; x1[1] = -1. + j*hb[1];
          x1[2] = -1. + k*hb[2]; x1[3] = -1. + l*hb[3];
          cg4(f_ball, par, x1, hb, xex, xgam);
          tot += xex[4];
        }
  ck("3-sphere, grid sum / 2 pi^2 R^3",
     tot/(2.*M_PI*M_PI*R*R*R), 1., 1.e-8);

  /* a grid sum can hide cancelling errors: check ONE cut cell against its
     sixteen children, for the measure and for the interface centroid    */
  memset(par, 0, sizeof(par));
  par[0] = 0.31; par[1] = -0.12; par[2] = 0.07; par[3] = 0.22; par[4] = 0.83;
  x0[0] = 0.15; x0[1] = -0.35; x0[2] = -0.05; x0[3] = 0.10;
  h[0] = 0.6; h[1] = 0.5; h[2] = 0.7; h[3] = 0.55;
  cg4(f_ball, par, x0, h, xex, xgam);
  {
    double hs[4], xs[4], xe2[8], xg2[4], msub = 0., csub[4] = {0.,0.,0.,0.};
    for (i = 0; i < 4; i++) hs[i] = 0.5*h[i];
    for (k = 0; k < 16; k++) {
      for (i = 0; i < 4; i++)
        xs[i] = x0[i] + ((k >> i) & 1)*hs[i];
      cg4(f_ball, par, xs, hs, xe2, xg2);
      msub += xe2[4];
      for (i = 0; i < 4; i++)
        csub[i] += xe2[4]*xg2[i];
    }
    ck("ball measure, cell vs 16 children", xex[4], msub, 1.e-7);
    worst = 0.;
    for (i = 0; i < 4; i++)
      worst = MAXD(worst,fabs(xgam[i] - csub[i]/msub));
    ck("ball interface centroid, same", worst, 0., 1.e-7);
  }
}

/* ------------------------------------------------------------------ */
int main(void)
{
  test_halfspace();
  test_axis_aligned();
  test_hyperball();
  test_slab();
  test_refinement();
  test_degenerate();
  test_cell_type();
  test_face_cap();
  test_api();
  test_measure();

  printf("\n%d checks, %d failed\n", ntest, nfail);

  return nfail ? 1 : 0;
}
