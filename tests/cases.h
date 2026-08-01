/**
 * @file      cases.h
 * @brief     Shared implicit functions and cell sweeps used by BOTH the
 *            golden-data generator and the regression test, so the two can
 *            never drift apart.
 *
 * The sweeps are deterministic: same order, same cells, same calls. Adding a
 * shape or changing a sweep INVALIDATES the golden file and requires
 * regenerating it -- which must only ever be done from a build known to be
 * good.
 *
 * Every implicit function reads exactly par[8] coordinates of x[]. The
 * library hands the callback an array of NDIM reals with the components above
 * ndim0 zeroed, so reading further would be legal but would couple the test
 * to that padding; reading exactly ndim0 keeps the callback honest and lets
 * AddressSanitizer catch a library that hands over a short array.
 */

#ifndef VOFI_TEST_CASES_H
#define VOFI_TEST_CASES_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------------- *
 * implicit functions: the reference phase is f < 0                        *
 * par[8] is always the dimension; the rest of the layout is per-shape     *
 * ---------------------------------------------------------------------- */

/* par: c[0..3], R  --  a ball/sphere/circle in any dimension  */
static double f_ball(const double x[], void *par)
{
  const double *p = (const double *) par;
  int nd = (int) p[8], i;
  double d = 0., t;
  for (i = 0; i < nd; i++) {
    t = x[i] - p[i];
    d += t*t;
  }
  return d - p[4]*p[4];
}

/* par: n[0..3], d  --  half-space n.x < d  */
static double f_plane(const double x[], void *par)
{
  const double *p = (const double *) par;
  int nd = (int) p[8], i;
  double s = 0.;
  for (i = 0; i < nd; i++)
    s += p[i]*x[i];
  return s - p[4];
}

/* par: c[0..2], a[0..2] (semi-axes), unused, nd  --  ellipsoid  */
static double f_ellipsoid(const double x[], void *par)
{
  const double *p = (const double *) par;
  int nd = (int) p[8], i;
  double d = 0., t;
  for (i = 0; i < nd; i++) {
    t = (x[i] - p[i])/p[4+i];
    d += t*t;
  }
  return d - 1.;
}

/* par: amp, k, off, dim, nd -- a wavy interface normal to axis (int)p[3]:
   f = x[dim] - off - amp*prod_{i != dim} cos(k*x[i])  */
static double f_wave(const double x[], void *par)
{
  const double *p = (const double *) par;
  int dim = (int) p[3], nd = (int) p[8], i;
  double w = 1.;
  for (i = 0; i < nd; i++)
    if (i != dim) w *= cos(p[1]*x[i]);
  return x[dim] - p[2] - p[0]*w;
}

/* par: c[0..2], R, r -- torus about the z axis (3D only) */
static double f_torus(const double x[], void *par)
{
  const double *p = (const double *) par;
  double dx = x[0]-p[0], dy = x[1]-p[1], dz = x[2]-p[2];
  double q = sqrt(dx*dx + dy*dy) - p[3];
  return q*q + dz*dz - p[4]*p[4];
}

/* par: c[0..2], R0, growth -- a 3D ball whose radius varies with the 4th
   coordinate: the space-time slab that motivates the 4D kernel (4D only) */
static double f_slab(const double x[], void *par)
{
  const double *p = (const double *) par;
  double d = 0., t, R;
  int i;
  for (i = 0; i < 3; i++) {
    t = x[i] - p[i];
    d += t*t;
  }
  R = p[3] + p[4]*x[3];
  return d - R*R;
}

/* ---------------------------------------------------------------------- *
 * the sweep                                                               *
 * ---------------------------------------------------------------------- */

typedef double (*vofi_test_func)(const double [], void *);

typedef struct {
  const char   *name;
  vofi_test_func func;
  double        par[9];   /* par[8] MUST be the dimension */
  int           ndim;
  int           n;        /* cells per direction                */
  double        xlo[4];   /* sweep origin                       */
  double        hcell[4]; /* cell size                          */
} vofi_case;

/* 2D and 3D sweeps -- these lock the PRE-EXISTING behaviour and must not
   change. New shapes go at the END of the list. */
static const vofi_case vofi_cases_23[] = {
  {"circle-2d",    f_ball,      {0.5,0.5,0.,0., 0.31,0.,0.,0., 2.}, 2, 8,
   {0.,0.,0.,0.}, {0.125,0.125,0.,0.}},
  {"circle-2d-c",  f_ball,      {0.5,0.5,0.,0., 0.12,0.,0.,0., 2.}, 2, 8,
   {0.,0.,0.,0.}, {0.125,0.125,0.,0.}},
  {"plane-2d",     f_plane,     {0.6,0.8,0.,0., 0.57,0.,0.,0., 2.}, 2, 8,
   {0.,0.,0.,0.}, {0.125,0.125,0.,0.}},
  {"ellipse-2d",   f_ellipsoid, {0.5,0.5,0.,0., 0.4,0.22,1.,1., 2.}, 2, 8,
   {0.,0.,0.,0.}, {0.125,0.125,0.,0.}},
  {"wave-2d",      f_wave,      {0.13,7.0,0.5,1., 0.,0.,0.,0., 2.}, 2, 8,
   {0.,0.,0.,0.}, {0.125,0.125,0.,0.}},
  {"sphere-3d",    f_ball,      {0.5,0.5,0.5,0., 0.33,0.,0.,0., 3.}, 3, 5,
   {0.,0.,0.,0.}, {0.2,0.2,0.2,0.}},
  {"sphere-3d-c",  f_ball,      {0.5,0.5,0.5,0., 0.09,0.,0.,0., 3.}, 3, 5,
   {0.,0.,0.,0.}, {0.2,0.2,0.2,0.}},
  {"plane-3d",     f_plane,     {0.48,0.62,0.62,0., 0.86,0.,0.,0., 3.}, 3, 5,
   {0.,0.,0.,0.}, {0.2,0.2,0.2,0.}},
  {"ellipsoid-3d", f_ellipsoid, {0.5,0.5,0.5,0., 0.42,0.3,0.2,1., 3.}, 3, 5,
   {0.,0.,0.,0.}, {0.2,0.2,0.2,0.}},
  {"torus-3d",     f_torus,     {0.5,0.5,0.5, 0.3,0.13, 0.,0.,0., 3.}, 3, 5,
   {0.,0.,0.,0.}, {0.2,0.2,0.2,0.}},
  {"wave-3d",      f_wave,      {0.1,6.0,0.5,2., 0.,0.,0.,0., 3.}, 3, 5,
   {0.,0.,0.,0.}, {0.2,0.2,0.2,0.}}
};

#define VOFI_NCASE_23 ((int)(sizeof(vofi_cases_23)/sizeof(vofi_cases_23[0])))

/* 4D sweeps. Small grids on purpose: a 4D cut cell costs ~2e4 function
   evaluations, so the point is coverage of the code paths, not statistics
   -- stress_4d does the statistics. */
static const vofi_case vofi_cases_4[] = {
  {"ball-4d",      f_ball,      {0.5,0.5,0.5,0.5, 0.62,0.,0.,0., 4.}, 4, 3,
   {0.,0.,0.,0.}, {0.34,0.34,0.34,0.34}},
  {"plane-4d",     f_plane,     {0.43,0.71,0.52,0.66, 1.05,0.,0.,0., 4.}, 4, 3,
   {0.,0.,0.,0.}, {0.34,0.34,0.34,0.34}},
  {"slab-4d",      f_slab,      {0.5,0.5,0.5, 0.45,0.35, 0.,0.,0., 4.}, 4, 3,
   {0.,0.,0.,0.}, {0.34,0.34,0.34,0.34}},
  {"ellipsoid-4d", f_ellipsoid, {0.5,0.5,0.5,0.5, 0.7,0.5,0.62,0.44, 4.}, 4, 3,
   {0.,0.,0.,0.}, {0.34,0.34,0.34,0.34}}
};

#define VOFI_NCASE_4 ((int)(sizeof(vofi_cases_4)/sizeof(vofi_cases_4[0])))

/* values recorded per cell: cc, xex[0..4] (centroid + interface measure;
   in 1D/2D/3D xex[4] is unused and recorded as zero), xgam[0..3]
   (interface centroid), cell type */
#define VOFI_NREC 11

#endif
