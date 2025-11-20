#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "vofi.h"

#define M_PI 3.14159265358979323846

/* mesh nodes (33 nodes -> 32 cells per direction) */
static const double x_nodes[33] = { -0.96875, -0.90625, -0.84375, -0.78125, -0.71875, -0.65625, -0.59375, -0.53125, -0.46875, -0.40625, -0.34375, -0.28125, -0.21875, -0.15625, -0.09375, -0.03125, 0.03125, 0.09375, 0.15625, 0.21875, 0.28125, 0.34375, 0.40625, 0.46875, 0.53125, 0.59375, 0.65625, 0.71875, 0.78125, 0.84375, 0.90625, 0.96875, 1.03125};
static const double y_nodes[33] = {-0.96875, -0.90625, -0.84375, -0.78125, -0.71875, -0.65625, -0.59375, -0.53125, -0.46875, -0.40625, -0.34375, -0.28125, -0.21875, -0.15625, -0.09375, -0.03125, 0.03125, 0.09375, 0.15625, 0.21875, 0.28125, 0.34375, 0.40625, 0.46875, 0.53125, 0.59375, 0.65625, 0.71875, 0.78125, 0.84375, 0.90625, 0.96875, 1.03125};

/* t nodes (5 nodes -> 4 time cells) */
static const double t_nodes[4] = {0.05859375, 0.06221, 0.076402425, 0.37223125};


/* time-dependent body: negative inside body */
double impl_time_body(const double xs[], void *par)
{
  /* use time from xs[2] (vofi supplies full 3D point) */
  double t = xs[2];
  double x = xs[0];
  double y = xs[1];

  /* domain extents (from nodes): Lx = Ly = x_nodes[32] - x_nodes[0] */
  double Lx = x_nodes[32] - x_nodes[0];
  double Ly = y_nodes[32] - y_nodes[0];
  double Lmin = (Lx < Ly) ? Lx : Ly;

  double radius0 = 0.15 * Lmin;
  double ampl = 0.05 * Lmin;
  double freq = 2.0 * M_PI;
  double center_x = 0.0, center_y = 0.0;

  double r = sqrt((x - center_x)*(x - center_x) + (y - center_y)*(y - center_y));
  double iso = r - (radius0 + ampl * sin(freq * t));

  /* vofi convention: return negative inside the body */
  return -iso;
}

int main(void)
{
  const int ncx = 32, ncy = 32, nt = 3;
  const double hx = x_nodes[1] - x_nodes[0];
  const double hy = y_nodes[1] - y_nodes[0];

  printf("3D integration over %dx%dx%d cells\n", ncx, ncy, nt);
  printf("hx=%g hy=%g\n", hx, hy);

  double total_volume = 0.0;

  /* arrays for vofi_get_cc call */
  double xin[3];
  double h0[3];
  double xex[4];
  int nex[2] = {1,1};
  int npt[1] = {0};
  int nvis[2] = {0,0};

  for (int kt = 0; kt < nt; ++kt) {
    double t0 = t_nodes[kt];
    double dt = t_nodes[kt+1] - t_nodes[kt];
    h0[2] = dt;

    for (int i = 0; i < ncx; ++i) {
      h0[0] = hx;
      xin[0] = x_nodes[i];
      for (int j = 0; j < ncy; ++j) {
        h0[1] = hy;
        xin[1] = y_nodes[j];
        xin[2] = t0;

        /* call vofi_get_cc in 3D mode: par is unused (impl reads xs[2]) */
        double cc = vofi_get_cc((integrand)impl_time_body, NULL, xin, h0,
                                xex, nex, npt, nvis, 3);
        if (isnan(cc) || isinf(cc)) {
          fprintf(stderr, "Non-finite cc at cell (%d,%d,%d): cc=%g\n", i, j, kt, cc);
          continue;
        }
        total_volume += cc * (hx * hy * dt);
      }
    }
    printf("time slab %d [t=%g .. %g] integrated partial volume so far = %.12e\n",
           kt, t0, t0+dt, total_volume);
  }

  printf("Total integrated 3D volume = %.12e\n", total_volume);
  return 0;
}