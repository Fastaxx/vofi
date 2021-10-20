/**
 * @file      checkconsistency.c
 * @authors   Andrea Chierici, Leonardo Chirco, Vincent Le Chenadec, 
 *            Ruben Scardovelli, Philip Yecko and Stéphane Zaleski 
 * @date      April 15, 2021
 * @brief     Routines to check consistency with a minimum along  a cell side 
 *            (2D), a cell face (3D), a line in the secondary direction (3D)
 * @version   Vofi 2.0
 * @copyright GPLv3 license
 **/

#include "vofi_stddecl.h"

/* -------------------------------------------------------------------------- *
 * GPL Licence                                                                *
 *                                                                            *
 *     This file is part of VOFI.                                             *
 *                                                                            *
 *     VOFI is free software: you can redistribute it and/or modify           *
 *     it under the terms of the GNU General Public License as published by   *
 *     the Free Software Foundation, either version 3 of the License, or      *
 *     (at your option) any later version.                                    *
 *                                                                            *
 *     VOFI is distributed in the hope that it will be useful,                *
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the           *
 *     GNU General Public License for more details.                           *
 *                                                                            *
 *     You should have received a copy of the GNU General Public License      *
 *     along with VOFI. If not, see <http://www.gnu.org/licenses/>.           *
 * -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * Routines to check consistency with a minimum along  a cell side (2D),      *
 *  a cell face (3D), a line in the secondary direction (3D)                  *
 * INPUT:  pointer to the implicit function, arrays of function parameters    * 
 *         par, coordinates of minimum vertex x0, side direction dir (2D) or  *
 *         dir1, dir2 (3D).function value on endpoints fse (2D) or on         *
 *         vertices fv (3D), cell edges h0 (3D) or side length h0 (2D)        *
 * OUTPUT: integer consi (2D) or dir_data ipsc (3D) (see vofi_stddecl.h)      *
 * -------------------------------------------------------------------------- */
vofi_int vofi_check_side_consistency(integrand impl_func,vofi_void_cptr par,
                                     vofi_creal x0[],vofi_creal dir[],
                                     vofi_creal fse[],vofi_creal h0)
{
  vofi_int i,consi;
  vofi_real xs[NDIM],dh,f0,f1,fs,ft;

  fs = fse[0]+fse[1];
  if (fs > 0.) 
    consi = 1;
  else if (fs < 0.) 
    consi = -1;
  else                                  /* - */
    consi = 0;

  /* - */
  /* - */
  if (consi != 0) {
    dh = MAX(EPS_M*h0,EPS_ROOT);
    f0 = fabs(fse[0]);
    f1 = fabs(fse[1]);
    if (f0 <= f1)
      ft = f0;
    else {
      ft = f1;
      dh = h0 - dh;
    }
    for (i=0; i<NDIM; i++)
      xs[i] = x0[i] + dh*dir[i];
    fs = consi*impl_func(xs,par);
    if (fs >= ft)           /* - */
      consi = 0;
  }
  
  return consi;
}

/* -------------------------------------------------------------------------- */
dir_data vofi_check_face_consistency(integrand impl_func,vofi_void_cptr par,
                         vofi_creal x0[],vofi_creal h0[],vofi_creal dir1[],
                         vofi_creal dir2[],vofi_creal fv[])
{
  vofi_int i,is1,is2,consi;
  vofi_real xx[NDIM],x1[NDIM],x2[NDIM],fl[NVER],f0,f1,f2,dh1,dh2,h1,h2;
  dir_data ipsc;
  
  h2 = h1 = f0 = 0.;
  ipsc.ind2 = ipsc.ind1 = ipsc.swt2 = ipsc.swt1 = 0;
  for (i=0;i<NDIM;i++) { 
    h1 += dir1[i]*h0[i];
    h2 += dir2[i]*h0[i];
  }
  for (i=0;i<NVER;i++)
    f0 += fv[i];
  
  if (f0 > 0.) 
    ipsc.consi = 1;
  else if (f0 < 0.) 
    ipsc.consi = -1;
  else                                    /* - */
    ipsc.consi = 0;

  /* - */
  /* - */
  if (ipsc.consi != 0) {
    dh1 = MAX(EPS_M*h1,EPS_ROOT);
    dh2 = MAX(EPS_M*h2,EPS_ROOT);
    f0 = fabs(f0);
    for (i=0;i<NVER;i++)
      fl[i] = fabs(fv[i]);
    
    if (fl[0] < f0) {
      f0 = fl[0];
      is1 = is2 = 1;
    }
    if (fl[1] < f0) {
      f0 = fl[1];
      ipsc.ind1 = 1;  
      is1 = -1; is2 = 1; 
    }  
    if (fl[2] < f0) {
      f0 = fl[2];
      ipsc.ind1 = 0; ipsc.ind2 = 1;  
      is1 = 1; is2 = -1; 
    }  
    if (fl[3] < f0) {
      f0 = fl[3];
      ipsc.ind1 = ipsc.ind2 = 1;
      is1 = is2 = -1;
    }
  
    consi = 0;
    for (i=0;i<NDIM;i++) {
      xx[i] = x0[i] + h1*ipsc.ind1*dir1[i] + h2*ipsc.ind2*dir2[i];
      x1[i] = xx[i] + dh1*is1*dir1[i];
      x2[i] = xx[i] + dh2*is2*dir2[i];
    }
    
    f1 = ipsc.consi*impl_func(x1,par);
    if (f1 < f0) {
      consi = ipsc.consi;
      ipsc.swt1 = 1;
    }
    
    f2 = ipsc.consi*impl_func(x2,par);
    if (f2 < f0) {
      consi = ipsc.consi;
      ipsc.swt2 = 1;
    }

    ipsc.consi = consi;
  }

  return ipsc;
}

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * Routine to check consistency with a region present beyond the interface    *
 * intersection with a cell tertiary side                                     *
 * EXTRA INPUT: which endpoint of the line n, min_data pointer                *
 * -------------------------------------------------------------------------- */
vofi_int vofi_check_line_consistency(integrand impl_func,vofi_void_cptr par,
                                     vofi_creal x0[],vofi_creal dir[],
                                     vofi_creal h0,vofi_cint n,min_data *xfs)
{
  vofi_int i,consi;
  vofi_real xs[NDIM],dh,f0,f1;

  consi = 0;
  dh = MAX(EPS_M*h0,EPS_ROOT);
  /* - */
  for (i=0; i<NDIM; i++)
    xs[i] = x0[i];
  for (i=0; i<NDIM; i++)
    xs[i] = x0[i] + (1.-2.*n)*h0*dir[i];
  f1 = impl_func(xs,par);
  for (i=0; i<NDIM; i++)
    xs[i] = x0[i] + (1.-2.*n)*dh*dir[i];
  f0 = impl_func(xs,par);
  if (f0*f1 <=0.) {
    consi = 1;
  for (i=0; i<NDIM; i++)
    xfs->xval[i] = xs[i];
  xfs->fval = f0;
  xfs->sval = dh;
  xfs->isc[0] = 1;
  xfs->isc[1] = 1;
  }
  
  return consi;
}

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * Routine to check consistency and compute the intersection along the secon- * 
 * dary side beyond the intersection on a tertiary side on the edge plane     *
 * -------------------------------------------------------------------------- */
void vofi_check_edge_consistency(integrand impl_func,vofi_void_cptr par,
                     vofi_real fse[],vofi_creal x0[],vofi_real base[],
                     vofi_creal dir[],vofi_creal h0,vofi_cint nsub)
{
  vofi_int i,f2neg;
  vofi_real s0[4],xs[NDIM],dh,dhl;

  dh = MAX(EPS_M*h0,EPS_ROOT);
  /* - */
  if (fabs(fse[0]) < fabs(fse[1])) { 
    for (i=0; i<NDIM; i++)
      xs[i] = x0[i] + dh*dir[i];
    fse[0] = impl_func(xs,par);
    if (fse[0]*fse[1] > 0.)
      base[nsub] = 0.;
    else {
      if (fse[0] < 0.)
        f2neg = 1;
      else
        f2neg = -1;
      s0[0] = h0 - dh;
      s0[1] = 0.;
      s0[2] = fse[0];
      s0[3] = (fse[1]-fse[0])/s0[0];
      dhl = vofi_get_segment_zero(impl_func,par,xs,dir,s0,f2neg);
      if (f2neg < 0)
        dhl = s0[0] - dhl;
      base[nsub] = dhl + dh;
    }
  }
  else {
    for (i=0; i<NDIM; i++)
      xs[i] = x0[i] + (h0-dh)*dir[i];
    fse[1] = impl_func(xs,par);
    if (fse[0]*fse[1] > 0.)
      base[nsub] = h0;
    else {
      if (fse[0] < 0.)
        f2neg = 1;
      else
        f2neg = -1;
      s0[0] = h0 - dh;
      s0[1] = s0[0];
      s0[2] = fse[1];
      s0[3] = (fse[1]-fse[0])/s0[0];
      dhl = vofi_get_segment_zero(impl_func,par,xs,dir,s0,f2neg);
      if (f2neg < 0)
        dhl = s0[0] - dhl;
      base[nsub] = dhl;
    }
  }
  
  return;
}
