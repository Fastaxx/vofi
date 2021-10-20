/**
 * @file      getintersections.c
 * @authors   Andrea Chierici, Leonardo Chirco, Vincent Le Chenadec, 
 *            Ruben Scardovelli, Philip Yecko and Stéphane Zaleski 
 * @date      April 15, 2021
 * @brief     Routines to compute the interface intersections with a 
 *            cell side and the external limits of a cap-like intersection
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
 * compute the interface intersections with a cell side (1D) or the external  *
 * limits of a cap-like intersection (2D), these will be internal/external    *
 * limits of integration                                                      *
 * INPUT:  pointer to the implicit function, arrays of function parameters    * 
 *         par, function value at endpoints fse, starting point x0, min_data  *
 *         xfsp (see vofi_stddecl.h), side subdivision base, side direction   * 
 *         dir (directions sdir, tdir), side length hl (cell edges h0),       * 
 *         number of input subdivision nsub, type of intersection isc         *
 * OUTPUT: number of new intersections inters, updated side subdivision base  *
 * FUNCTIONS:                                                                 *
 * vofi_int vofi_get_side_intersections: 1D problem                           *
 * vofi_int vofi_get_ext_intersections : 2D problem                           *
 * -------------------------------------------------------------------------- */
vofi_int vofi_get_side_intersections(integrand impl_func,vofi_void_cptr par,
                       vofi_real fse[],vofi_creal x0[],min_data xfsp,
                       vofi_real base[],vofi_creal dir[],vofi_creal hl,
                       vofi_int nsub,vofi_cint isc)
{
  vofi_int f2neg,inters;
  vofi_real s0[4],dhl;    

  inters = 0;
  /* - */
  if (isc < 0) {
    if (fse[0] < 0.)
      f2neg = 1;
    else
      f2neg = -1;
    s0[0] = hl;
    if ( fabs(fse[0]) < fabs(fse[1]) ) {
      s0[1] = 0.;
      s0[2] = fse[0];
    }
    else {
      s0[1] = hl;
      s0[2] = fse[1];
    }
    s0[3] = (fse[1] - fse[0])/hl;
    dhl = vofi_get_segment_zero(impl_func,par,x0,dir,s0,f2neg);
    if (f2neg < 0)
      dhl = hl - dhl;
    base[nsub] = dhl;
    inters++;
  }
  /* - */
  else {
    if ((fse[0]+fse[1]) <= 0.)
      f2neg = 1;
    else
      f2neg = -1;
    s0[0] = xfsp.sval;
    if ( fabs(fse[0]) < fabs(xfsp.fval) ) {
      s0[1] = 0.;
      s0[2] = fse[0];
    }
    else {
      s0[1] = xfsp.sval;
      s0[2] = xfsp.fval;
    }
    s0[3] = (xfsp.fval - fse[0])/xfsp.sval;
    dhl = vofi_get_segment_zero(impl_func,par,x0,dir,s0,f2neg);
    if (fse[0] > 0.0 || xfsp.fval < 0.0)
      dhl = xfsp.sval - dhl;
    base[nsub] = dhl;
    inters++;
    f2neg = - f2neg;
    s0[0] = hl - xfsp.sval;
    if ( fabs(xfsp.fval) < fabs(fse[0]) ) {
      s0[1] = 0.;
      s0[2] = xfsp.fval;
    }
    else {
      s0[1] = s0[0];
      s0[2] = fse[1];
    }
    s0[3] = (fse[1]-xfsp.fval)/s0[0];
    dhl = vofi_get_segment_zero(impl_func,par,xfsp.xval,dir,s0,f2neg);
    if (xfsp.fval > 0.0 || fse[1] < 0.0)
      dhl = s0[0] - dhl;
    base[nsub+1] = xfsp.sval + dhl;
    inters++;
  }
  
  return inters;
}    

/* -------------------------------------------------------------------------- */
vofi_int vofi_get_ext_intersections(integrand impl_func,vofi_void_cptr par,
              vofi_creal x0[],vofi_creal h0[],min_data xfsp,vofi_real base[],
              vofi_creal sdir[],vofi_creal tdir[],vofi_cint nsub)
{
  vofi_int i,k,iter,js,jt,not_conv,ipt,ist,f2neg,inters;
  vofi_real pt0[NDIM],pt1[NDIM],pt2[NDIM],ptt[NDIM],mp0[NDIM],mp1[NDIM];
  vofi_real ss[NDIM],ext_dir[NDIM],int_dir[NDIM],fse[NSE];
  vofi_real s0[4],ss0,ds0,fpt0,sss,sst,ssx,ssy,tol2,normdir,d1,d2,a1,a2;
  vofi_creal tol = EPS_M; 
 
  tol2 = 2.*tol;
  js = jt = 0;
  inters = 0;
  for (i=0;i<NDIM;i++)
    pt0[i] = xfsp.xval[i];
  if (xfsp.fval < 0.)
    f2neg = 1;
  else
    f2neg = -1;
  fse[0] = f2neg*xfsp.fval;

  /* - */
  for (i=0;i<NDIM;i++) {
    int_dir[i] = sdir[i];
    js += i* (vofi_int) sdir[i];
    jt += i* (vofi_int) tdir[i];
  }

  /* - */
  for (i=0;i<NDIM;i++)
    pt2[i] = pt1[i] = pt0[i];
  pt2[js] = x0[js] + h0[js];            
  s0[0] = pt2[js] - pt0[js];         
  fse[1] = f2neg*impl_func(pt2,par);
  if (fse[1] > 0.) {
    if ( fabs(fse[0]) < fabs(fse[1]) ) {
      s0[1] = 0.;
      s0[2] = fse[0];
    }
    else {
      s0[1] = s0[0];
      s0[2] = fse[1];
    }
    s0[3] = (fse[1]-fse[0])/s0[0];
    ds0 = vofi_get_segment_zero(impl_func,par,pt0,int_dir,s0,f2neg);
    pt2[js] = pt0[js] + ds0;
  } 

  /* - */
  pt1[js] = x0[js];                  
  s0[0] = pt0[js] - pt1[js];    
  int_dir[js] = -1.;
  fse[1] = f2neg*impl_func(pt1,par);
  if (fse[1] > 0.) {
    if ( fabs(fse[0]) < fabs(fse[1]) ) {
      s0[1] = 0.;
      s0[2] = fse[0];
    }
    else {
      s0[1] = s0[0];
      s0[2] = fse[1];
    }
    s0[3] = (fse[1]-fse[0])/s0[0];
    ds0 = vofi_get_segment_zero(impl_func,par,pt0,int_dir,s0,f2neg);
    pt1[js] = pt0[js] - ds0;
  }

  /* - */
  for (i=0;i<NDIM;i++)
    pt0[i] = 0.5*(pt1[i] + pt2[i]);                  
  fpt0 = f2neg*impl_func(pt0,par);
  ss0 = pt2[js]-pt1[js];

  /* - */
  for (k=-1;k<=1;k=k+2) {
    iter = 0;
    not_conv = 1;
    /* - */
    for (i=0;i<NDIM;i++)                     
      ext_dir[i] = tdir[i];
    ext_dir[jt] = k;
    if (k < 0)
      sst = pt0[jt] - x0[jt];   
    else
      sst = x0[jt] + h0[jt] - pt0[jt];   
    for (i=0;i<NDIM;i++) {
      mp1[i] = pt0[i];
      pt1[i] = mp1[i] + sst*ext_dir[i];
    }
    fse[0] = fpt0;
    fse[1] = f2neg*impl_func(pt1,par);
    sss = ss0;
    /* - */
    while (not_conv  && iter < MAX_ITER_MINI) {    
      s0[0] = sst;
      if (fse[1] > 0.) {
        if ( fabs(fse[0]) < fabs(fse[1]) ) {
          s0[1] = 0.;
          s0[2] = fse[0];
        }
        else {
          s0[1] = s0[0];
          s0[2] = fse[1];
        }
        s0[3] = (fse[1]-fse[0])/s0[0];
        ds0 = vofi_get_segment_zero(impl_func,par,mp1,ext_dir,s0,f2neg);
        sst = ds0;
      }
      /* - */
      for (i=0;i<NDIM;i++) {
        pt1[i] = mp1[i] + sst*ext_dir[i];         
        mp0[i] = mp1[i];
        mp1[i] = pt2[i] = ptt[i] = pt1[i];
      }

      /* - */
      ipt = ist = 0;             
      ptt[js] += tol;            
      fse[0] = f2neg*impl_func(ptt,par);
      if (fse[0] < 0.) {
        ipt = 1;
        ssx = x0[js] + h0[js] - ptt[js];
        int_dir[js] = 1.;
      }
      else {
        ptt[js] -= tol2;
        fse[0] = f2neg*impl_func(ptt,par);
        if (fse[0] < 0.) {
          ipt = 1;
          ssx = ptt[js] - x0[js];
          int_dir[js] = -1.;
        }
      }
      /* - */
      if (ipt) {             
        sss = MIN(1.2*sss,ssx); 
        for (i=0;i<NDIM;i++) 
          pt2[i] = ptt[i] + sss*int_dir[i];
        fse[1] = f2neg*impl_func(pt2,par);
        while (fse[1] < 0. && ist < 3 && sss < ssx) {
          sss = MIN(3.*sss,ssx);
          if (ist == 2)
            sss = ssx;
          for (i=0;i<NDIM;i++) 
            pt2[i] = ptt[i] + sss*int_dir[i];
          fse[1] = f2neg*impl_func(pt2,par);
          ist++;
        }
        /* - */
        if (fse[0]*fse[1] < 0.) {         
          s0[0] = sss;
          if ( fabs(fse[0]) < fabs(fse[1]) ) {
            s0[1] = 0.;
            s0[2] = fse[0];
          }
          else {
            s0[1] = s0[0];
            s0[2] = fse[1];
          }
          s0[3] = (fse[1]-fse[0])/s0[0];
          ds0 = vofi_get_segment_zero(impl_func,par,ptt,int_dir,s0,f2neg);
          for (i=0;i<NDIM;i++) 
            pt2[i] = ptt[i] + ds0*int_dir[i];
        }
        /* - */
        for (i=0;i<NDIM;i++)                        
          mp1[i] = 0.5*(pt1[i] + pt2[i]);
        fse[0] = f2neg*impl_func(mp1,par);
        sss = fabs(pt1[js]-pt2[js]);
      }
      /* - */
      normdir = 0.;                             
      for (i=0;i<NDIM;i++) {
        ext_dir[i] =  mp1[i] - mp0[i];                      
        normdir += ext_dir[i]*ext_dir[i];
      }  
      normdir = sqrt(normdir) + EPS_NOT0; 
      for (i=0;i<NDIM;i++) {              
        ext_dir[i] = ext_dir[i]/normdir;
        d1 = SGN0P(ext_dir[i]);
        d2 = fabs(ext_dir[i]) + EPS_NOT0;
        a1 = (x0[i] - mp1[i])/(d1*d2);
        a2 = (x0[i] + h0[i] - mp1[i])/(d1*d2);
        ss[i] = MAX(a1,a2);
      }
      ssy = MIN(ss[0],ss[1]);      
      ssy = MIN(ssy,ss[2]);             
      sst = MIN(1.2*sst,ssy);
      /* - */
      if (!ipt || sss < tol2 || sst < EPS_ROOT) {       
        not_conv = 0;                       
      }
      else {
        for (i=0;i<NDIM;i++) 
          pt1[i] = mp1[i] + sst*ext_dir[i];
        fse[1] = f2neg*impl_func(pt1,par);
        ist = 0;
        /* - */
        while (fse[1] < 0. && ist < 3 && sst < ssy) { 
          sst = MIN(3.*sst,ssy);
          if (ist == 2)
            sst = ssy;
          for (i=0;i<NDIM;i++) 
            pt1[i] = mp1[i] + sst*ext_dir[i];
          fse[1] = f2neg*impl_func(pt1,par);
          ist++;
        }
      }  
      iter++;
    }
    base[nsub+inters] = mp1[jt] - x0[jt];
    inters++;
  }

  return inters;
}
