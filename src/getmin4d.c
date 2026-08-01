/**
 * @file      getmin4d.c
 * @brief     Consistency check and minimum search inside a cubic 3-face of
 *            a 4D cell: the three-direction extension of
 *            vofi_check_face_consistency and vofi_get_face_min.
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
 * The boundary of a 4D cell is made of eight CUBIC 3-faces, so an interface  *
 * that reaches the cell boundary without changing the sign at any vertex is  *
 * found by minimising f inside a 3-face -- one dimension more than the face  *
 * minimisation the 3D kernel needs. Both routines below are the term-by-term *
 * extension of their 2D counterparts in checkconsistency.c and getmin.c:     *
 * same starting-vertex rule, same conjugate-gradient outer loop, same        *
 * Brent line search, same early exit the moment f goes negative (the sign    *
 * change is all the caller asked for; the true minimum is not needed).       *
 * INPUT:  pointer to the implicit function, arrays of function parameters    *
 *         par, corner x0, cell edges h0, the three free directions dir1,     *
 *         dir2, dir3, function values fve at the eight vertices of the       *
 *         3-face, indexed ind1 + 2*ind2 + 4*ind3                             *
 * OUTPUT: dir_data ipsc / integer flag sign_change and min_data *xfs_pt      *
 * FUNCTIONS:                                                                 *
 * dir_data vofi_check_cell_consistency: is a minimum inside possible at all  *
 * vofi_int vofi_get_cell_min          : 3D CG + Brent, stop on a sign change *
 * -------------------------------------------------------------------------- */
dir_data vofi_check_cell_consistency(integrand impl_func,vofi_void_cptr par,
                         vofi_creal x0[],vofi_creal h0[],vofi_creal dir1[],
                         vofi_creal dir2[],vofi_creal dir3[],vofi_creal fve[])
{
  vofi_int i,k,v,imin,consi,ind[NDIM3],is[NDIM3];
  vofi_real xx[NDIM],xp[NDIM],hh[NDIM3],dh,f0,fmin,fk;
  vofi_creal *dirs[NDIM3];
  dir_data ipsc;

  dirs[0] = dir1; dirs[1] = dir2; dirs[2] = dir3;
  ipsc.ind1 = ipsc.ind2 = ipsc.ind3 = 0;
  ipsc.swt1 = ipsc.swt2 = ipsc.swt3 = 0;
  for (k=0;k<NDIM3;k++) {
    hh[k] = 0.;
    for (i=0;i<NDIM;i++)
      hh[k] += dirs[k][i]*h0[i];
  }

  f0 = 0.;
  for (v=0;v<NVERC;v++)
    f0 += fve[v];

  if (f0 > 0.)
    ipsc.consi = 1;
  else if (f0 < 0.)
    ipsc.consi = -1;
  else                                    /* - */
    ipsc.consi = 0;

  if (ipsc.consi != 0) {
    /* start from the vertex where |f| is smallest and walk inwards along
       each free direction: a minimum inside is possible only if f drops
       below its vertex value along at least one of them                  */
    imin = 0;
    fmin = fabs(fve[0]);
    for (v=1;v<NVERC;v++)
      if (fabs(fve[v]) < fmin) {
        fmin = fabs(fve[v]);
        imin = v;
      }
    for (k=0;k<NDIM3;k++) {
      ind[k] = (imin >> k) & 1;
      is[k]  = (ind[k] == 0) ? 1 : -1;
    }
    ipsc.ind1 = ind[0]; ipsc.ind2 = ind[1]; ipsc.ind3 = ind[2];

    for (i=0;i<NDIM;i++)
      xx[i] = x0[i] + hh[0]*ind[0]*dir1[i] + hh[1]*ind[1]*dir2[i]
                    + hh[2]*ind[2]*dir3[i];

    consi = 0;
    for (k=0;k<NDIM3;k++) {
      dh = MAX(EPS_M*hh[k],EPS_ROOT);
      for (i=0;i<NDIM;i++)
        xp[i] = xx[i] + dh*is[k]*dirs[k][i];
      fk = ipsc.consi*impl_func(xp,par);
      if (fk < fmin) {
        consi = ipsc.consi;
        if (k == 0)      ipsc.swt1 = 1;
        else if (k == 1) ipsc.swt2 = 1;
        else             ipsc.swt3 = 1;
      }
    }
    ipsc.consi = consi;
  }

  return ipsc;
}

/* -------------------------------------------------------------------------- */
/* central differences of f along the three free directions, in the sign
   convention that makes f positive at the starting point                    */
static void vofi_cell_derivs(integrand impl_func,vofi_void_cptr par,
                             vofi_creal xs0[],vofi_creal *dirs[],
                             vofi_creal dh,vofi_cint f2pos,vofi_creal fs0,
                             vofi_real df[],vofi_real d2f[])
{
  vofi_int i,k,neg;
  vofi_real xf[NDIM],xb[NDIM],ff,fb;

  neg = 0;
  for (k=0;k<NDIM3;k++) {
    for (i=0;i<NDIM;i++) {
      xf[i] = xs0[i] + dh*dirs[k][i];
      xb[i] = xs0[i] - dh*dirs[k][i];
    }
    ff = f2pos*impl_func(xf,par);
    fb = f2pos*impl_func(xb,par);
    df[k]  = -0.5*(ff - fb)/dh;
    d2f[k] = (ff + fb - 2.*fs0)/(dh*dh);
    if (d2f[k] <= 0.)
      neg = 1;
  }
  if (neg)                                  /* - */
    for (k=0;k<NDIM3;k++)
      d2f[k] = 1.;

  return;
}

/* -------------------------------------------------------------------------- */
vofi_int vofi_get_cell_min(integrand impl_func,vofi_void_cptr par,
                           vofi_creal x0[],vofi_creal h0[],vofi_creal dir1[],
                           vofi_creal dir2[],vofi_creal dir3[],
                           vofi_creal fve[],min_data *xfs_pt,dir_data ipsc)
{
  vofi_int i,k,kk,not_conv,iter,iss,f2pos,sign_change = 0;
  vofi_real xs0[NDIM],xs1[NDIM];
  vofi_real res[NDIM],hes[NDIM],rs0[NDIM],hs0[NDIM];
  vofi_real pcrs[NDIM],nmdr[NDIM],cndr[NDIM],ss[NDIM],fse[NSE];
  vofi_real df[NDIM3],d2f[NDIM3];
  vofi_int  ind[NDIM3],swt[NDIM3];
  vofi_creal *dirs[NDIM3];
  vofi_real eps2,fs0,mcd,ss0,ss1,beta;
  vofi_real del0,delnew,delold,delmid,d1,d2,a1,a2;
  vofi_creal dh = 1.e-04;         /* - */

  dirs[0] = dir1; dirs[1] = dir2; dirs[2] = dir3;
  ind[0] = ipsc.ind1; ind[1] = ipsc.ind2; ind[2] = ipsc.ind3;
  swt[0] = ipsc.swt1; swt[1] = ipsc.swt2; swt[2] = ipsc.swt3;

  eps2 = EPS_E*EPS_E;
  for (i=0;i<NDIM;i++) {                           /* - */
    xs0[i] = x0[i] + h0[i]*(ind[0]*dir1[i] + ind[1]*dir2[i] + ind[2]*dir3[i]);
    rs0[i] = 0.;
    hs0[i] = 1. - dir1[i] - dir2[i] - dir3[i];
  }
  fse[0] = fve[ipsc.ind1 + 2*ipsc.ind2 + 4*ipsc.ind3];
  f2pos = ipsc.consi;
  fs0 = f2pos*fse[0];
  vofi_cell_derivs(impl_func,par,xs0,dirs,dh,f2pos,fs0,df,d2f);

  mcd = del0 = 0.;
  for (i=0;i<NDIM;i++) {
    res[i] = rs0[i];
    hes[i] = hs0[i];
    for (k=0;k<NDIM3;k++) {
      res[i] += df[k]*swt[k]*dirs[k][i];
      hes[i] += d2f[k]*dirs[k][i];
    }
    pcrs[i] = res[i]/hes[i];                       /* - */
    mcd += pcrs[i]*pcrs[i];
    del0 += res[i]*pcrs[i];                        /* - */
  }

  /* - */
  mcd = sqrt(mcd + EPS_NOT0);
  for (i=0;i<NDIM;i++) {
    cndr[i] = pcrs[i];                             /* - */
    nmdr[i] = cndr[i]/mcd;                         /* - */
    d1 = SGN0P(nmdr[i]);
    d2 = fabs(nmdr[i]) + EPS_NOT0;
    if (d2 < EPS_ROOT)
      ss[i] = SS_FREE;                             /* - */
    else {
      a1 = (x0[i] - xs0[i])/(d1*d2);
      a2 = (x0[i] + h0[i] - xs0[i])/(d1*d2);
      ss[i] = MAX(a1,a2);
    }
  }
  ss0 = ss[0];
  for (i=1;i<NDIM;i++)
    ss0 = MIN(ss0,ss[i]);
  for (i=0;i<NDIM;i++)
    xs1[i] = xs0[i] + ss0*nmdr[i];
  fse[1] = impl_func(xs1,par);

  delnew = del0;
  not_conv = 1;
  iter = kk = 0;
  while (not_conv && iter < MAX_ITER_MINI) {               /* - */
    sign_change = vofi_get_segment_min(impl_func,par,xs0,nmdr,fse,xfs_pt,
                                       ss0,f2pos);
    for (i=0;i<NDIM;i++)
      xs0[i] = xfs_pt->xval[i];
    fse[0] = xfs_pt->fval;
    fs0 = f2pos*fse[0];
    if (sign_change)
      not_conv = 0;                                   /* - */
    else {
      ss0 = xfs_pt->sval;
      vofi_cell_derivs(impl_func,par,xs0,dirs,dh,f2pos,fs0,df,d2f);
      delold = delnew;
      delmid = delnew = 0.;
      for (i=0;i<NDIM;i++) {
        res[i] = rs0[i];
        hes[i] = hs0[i];
        for (k=0;k<NDIM3;k++) {
          res[i] += df[k]*dirs[k][i];
          hes[i] += d2f[k]*dirs[k][i];
        }
        delmid += res[i]*pcrs[i];
        pcrs[i] = res[i]/hes[i];
        delnew += res[i]*pcrs[i];
      }

      beta = (delnew-delmid)/delold;
      kk++;

      if (kk == NDIM3 || beta <= 0.) {
        beta = 0.;
        kk = 0;
      }
      mcd = 0.;
      for (i=0;i<NDIM;i++) {
        cndr[i] = pcrs[i] + beta*cndr[i];
        mcd += cndr[i]*cndr[i];
      }
      mcd = sqrt(mcd + EPS_NOT0);
      for (i=0;i<NDIM;i++) {
        nmdr[i] = cndr[i]/mcd;                   /* - */
        d1 = SGN0P(nmdr[i]);
        d2 = fabs(nmdr[i]) + EPS_NOT0;
        if (d2 < EPS_ROOT)
          ss[i] = SS_FREE;
        else {
          a1 = (x0[i] - xs0[i])/(d1*d2);
          a2 = (x0[i] + h0[i] - xs0[i])/(d1*d2);
          ss[i] = MAX(a1,a2);
        }
      }
      ss1 = ss[0];
      for (i=1;i<NDIM;i++)
        ss1 = MIN(ss1,ss[i]);
      ss0 = MIN(1.2*ss0,ss1);

      /* - */
      if (delnew < eps2*del0 || ss0 < EPS_ROOT)
        not_conv = 0;
      else {
        for (i=0;i<NDIM;i++)
          xs1[i] = xs0[i] + ss0*nmdr[i];
        fse[1] = impl_func(xs1,par);
        iss = 0;                                      /* - */
        while (f2pos*fse[1] < fs0 && iss < 3 && ss0 < ss1) {
          ss0 = MIN(3.*ss0,ss1);
          if (iss == 2)
            ss0 = ss1;
          for (i=0;i<NDIM;i++)
            xs1[i] = xs0[i] + ss0*nmdr[i];
          fse[1] = impl_func(xs1,par);
          iss++;
        }
      }
      iter++;
    }
  }

  return sign_change;
}
