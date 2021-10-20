/**
 * @file      getmin.c
 * @authors   Andrea Chierici, Leonardo Chirco, Vincent Le Chenadec, 
 *            Ruben Scardovelli, Philip Yecko and Stéphane Zaleski 
 * @date      April 15, 2021
 * @brief     Routines to compute the function minimum either along a segment 
 *            or in a cell face, the search is stopped if a sign change is 
 *            detected 
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
 * compute the function minimum either along a segment or in a cell face,     *
 * the search is immediately stopped if a sign change is detected             *
 * INPUT:  pointer to the implicit function, arrays of function parameters    *
 *         par, starting point x0, segment direction dir in 2D (plane         *
 *         directions dir1 and dir2 in 3D), function value fse at endpoints   * 
 *         (function value fve at 4 vertices), segment length s0 (cell edges  *
 *         h0), sign ifsign to have function positive at endpoints            *
 *         (dir_data ipsc)                                                    *
 * OUTPUT: integer flag sign_change (1/0: detected or not a sign change),     *
 *         min_data *xfs_pt (see vofi_stddecl.h)                              *
 * FUNCTIONS:                                                                 *
 * vofi_int vofi_get_segment_min: 1D problem with Brent's method              *
 * vofi_int vofi_get_face_min: 2D with CG + Brent's method                    *
 * -------------------------------------------------------------------------- */
vofi_int vofi_get_segment_min(integrand impl_func,vofi_void_cptr par,
                          vofi_creal x0[],vofi_creal dir[],vofi_creal fse[],
                          min_data *xfs_pt,vofi_creal s0,vofi_cint ifsign)
{
  vofi_int i,j,iter,not_conv,igold,iseca,sign_change;
  vofi_real xs[NDIM],fs,fu,ft,fv,fa,fb,p,q,r;
  vofi_real sa,sb,ss,su,st,sv,se,sd,sc,sm,fm,sp,fp,sz,fz;
  vofi_real GRIS,tol,t2;
 
  /* - */
  GRIS = 0.5*(3.0 - sqrt(5.0));
  sign_change = 0;
  igold = 1;
  not_conv = 1;
  
  /* - */
  fa = ifsign*fse[0];
  sa = 0.;
  fb = ifsign*fse[1];
  sb = s0;
  if (fa <= fb) {
    st = 0.;
    sv = s0;
    ft = fa;
    fv = fb;
  }
  else {
    st = s0;
    sv = 0.;
    ft = fb;
    fv = fa;
  }

  ss = sa + GRIS*(sb - sa);
  for (i=0; i<NDIM; i++)
    xs[i] = x0[i] + ss*dir[i];
  fs = ifsign*impl_func(xs,par);
  if (fs > ft) {
    SHFT4(fu,ft,fs,fu);
    SHFT4(su,st,ss,su);
  }
  se = st - sv;
  sd = ss - st;
  if (fs < 0.)                            /* - */
    not_conv = 0;
 
  iter = 0;   
  while (not_conv  && iter < MAX_ITER_MINI) { 
    sc = 0.5*(sa + sb);
    tol = EPS_M*fabs(ss) + EPS_LOC;
    t2 = 2.0*tol;

    /* - */
    if ( fabs(ss - sc) <= t2 - 0.5*(sb - sa) )      
      not_conv = 0;
    /* - */
    else {
      iter++;
      r = 0.0;
      q = r;
      p = q;
      if (fabs(se) > tol) {
        r = (ss - st)*(fs - fv);
        q = (ss - sv)*(fs - ft);
        p = (ss - sv)*q - (ss - st)*r;
        q = 2.0*(q - r);
        if (q > 0.0)
          p = - p;
        else
          q = fabs(q);
        r  = se;
        se = sd;
      }
      /* - */
      /* - */
      if ( fabs(p) < fabs(0.5*q*r) && p > q*(sa - ss) && p < q*(sb - ss) ) {
        sd = p/q;
        su = ss + sd;
        if ((su - sa) < t2 || (sb - su) < t2) {
          if (ss < sc)
            sd = tol;
          else
            sd = - tol;
        }
      }  
      /* - */
      else {
        if (ss < sc)
          se = sb - ss;
        else
          se = sa - ss;
        sd = GRIS*se; 
        igold++;
      }
      /* - */
      if (fabs(sd) >= tol)
        su = ss + sd;
      else if (sd > 0.0)
        su = ss + tol;
      else
        su = ss - tol;
      
      /* - */
      for (i=0; i<NDIM; i++)
         xs[i] = x0[i] + su*dir[i];
      fu = ifsign*impl_func(xs,par);
      if (fu < 0.)                                   /* - */
        not_conv = 0;
      /* - */
      if (fu <= fs) { 
        if (su < ss) 
          { CPSF(sb,ss,fb,fs); }
        else 
          { CPSF(sa,ss,fa,fs); }

        SHFT4(fv,ft,fs,fu);
        SHFT4(sv,st,ss,su);
      }
      /* - */
      else {
        if (su < ss)
          { CPSF(sa,su,fa,fu); }
        else 
          { CPSF(sb,su,fb,fu); }
        
        if (fu <= ft || st == ss) {
          CPSF(sv,st,fv,ft);
          CPSF(st,su,ft,fu);
        }
        else if (fu <= fv || sv == ss || sv == st)
          { CPSF(sv,su,fv,fu); }
      }
      /* - */
      if (igold == 2 && not_conv == 1) {
        igold = 0;
        CPSF(sm,ss,fm,fs);
        if (st < sm && fabs(st-sa)>t2)
          { CPSF(sm,st,fm,ft); }
        if (sv < sm && fabs(sv-sa)>t2)
          { CPSF(sm,sv,fm,fv); }          

        CPSF(sp,ss,fp,fs);
        if (st > sp && fabs(st-sb)>t2)
          { CPSF(sp,st,fp,ft); }
        if (sv > sp && fabs(sv-sb)>t2)
          { CPSF(sp,sv,fp,fv); }

        p = (sa-sm)*(fp*sb - fb*sp) + (sp-sb)*(fm*sa - fa*sm);
        q = (sa-sm)*(fp-fb) + (sp-sb)*(fm-fa);
        if (q < 0.0) {
          p = - p;
          q = - q;
        }

        /* - */
        sm = MIN(sa,sm);
        sp = MAX(sb,sp);
        if ( p > q*sm  && p < q*sp ) {
          su = p/q;
          for (i=0; i<NDIM; i++)
            xs[i] = x0[i] + su*dir[i];
          fu = ifsign*impl_func(xs,par);
          iseca = 0;
          if (fu < fs) {    
            if (fu < 0.)                             /* - */
              not_conv = 0;
            tol = EPS_M*fabs(su) + EPS_LOC;  
            for (j=-1;j<=1;j=j+2) {  
              sz = su +j*tol;
              for (i=0; i<3; i++)
                xs[i] = x0[i] + sz*dir[i];
              fz = ifsign*impl_func(xs,par);
              if (fz > fu)
                iseca++;
            }
            if (iseca == 2) {
                  CPSF(ss,su,fs,fu);
              sb = sz;
              sa = sz - 2.*tol;
            }
          }
        }
      } /* - */
    }   /* - */
  }     /* - */
  
  for (i=0; i<NDIM; i++)
    xfs_pt->xval[i] = xs[i];
  xfs_pt->fval = ifsign*fs;
  xfs_pt->sval = ss;
  if (fs < 0.) 
    sign_change = 1;
  
  return sign_change;
}

/* -------------------------------------------------------------------------- */
vofi_int vofi_get_face_min(integrand impl_func,vofi_void_cptr par,
                           vofi_creal x0[],vofi_creal h0[],vofi_creal dir1[],
			   vofi_creal dir2[],vofi_creal fve[],min_data *xfs_pt,
			   dir_data ipsc)
{
  vofi_int i,not_conv,iter,k,ipt,iss,f2pos,sign_change;
  vofi_real xs0[NDIM],xs1[NDIM],x1f[NDIM],x1b[NDIM],x2f[NDIM],x2b[NDIM];
  vofi_real res[NDIM],hes[NDIM],rs0[NDIM],hs0[NDIM];
  vofi_real pcrs[NDIM],nmdr[NDIM],cndr[NDIM],ss[NDIM],fse[NSE];
  vofi_real eps2,fs0,f1f,f1b,f2f,f2b,df1,df2,d2f1,d2f2,mcd,ss0,ss1,beta;
  vofi_real del0,delnew,delold,delmid,d1,d2,a1,a2;
  vofi_creal dh = 1.e-04;         /* - */

  eps2 = EPS_E*EPS_E;
  for (i=0;i<NDIM;i++) {                           /* - */
    xs0[i] =  x0[i] + h0[i]*(ipsc.ind1*dir1[i] + ipsc.ind2*dir2[i]);
    x1f[i] = xs0[i] + dh*dir1[i];
    x1b[i] = xs0[i] - dh*dir1[i];
    x2f[i] = xs0[i] + dh*dir2[i];
    x2b[i] = xs0[i] - dh*dir2[i];
    rs0[i] = 0.;
    hs0[i] = 1. - dir1[i] - dir2[i];
  }
  fse[0] = fve[ipsc.ind1+2*ipsc.ind2];
  f2pos = ipsc.consi;
  fs0 = f2pos*fse[0];
  f1f = f2pos*impl_func(x1f,par);
  f1b = f2pos*impl_func(x1b,par);
  f2f = f2pos*impl_func(x2f,par);
  f2b = f2pos*impl_func(x2b,par);
  
  /* - */
  df1 = -0.5*(f1f-f1b)/dh;
  df2 = -0.5*(f2f-f2b)/dh;
  d2f1 = (f1f+f1b-2.*fs0)/(dh*dh);
  d2f2 = (f2f+f2b-2.*fs0)/(dh*dh);
  if (d2f1 <= 0. || d2f2 <= 0.)             /* - */
    d2f1 = d2f2 = 1.;

  mcd = del0 = 0.;                
  for (i=0;i<NDIM;i++) {
    res[i] = rs0[i] + df1*ipsc.swt1*dir1[i] + df2*ipsc.swt2*dir2[i];
    hes[i] = hs0[i] + d2f1*dir1[i] + d2f2*dir2[i];
    pcrs[i] = res[i]/hes[i];                       /* - */
    mcd += pcrs[i]*pcrs[i];
    del0 += res[i]*pcrs[i];                               /* - */
  }

  /* - */
  mcd = sqrt(mcd + EPS_NOT0); 
  for (i=0;i<NDIM;i++) {
    cndr[i] = pcrs[i];                               /* - */
    nmdr[i] = cndr[i]/mcd;                           /* - */
    d1 = SGN0P(nmdr[i]);
    d2 = fabs(nmdr[i]) + EPS_NOT0;
    if (d2 < EPS_ROOT)
      ss[i] = 1000.*h0[i];
    else {
      a1 = (x0[i] - xs0[i])/(d1*d2);
      a2 = (x0[i] + h0[i] - xs0[i])/(d1*d2);
      ss[i] = MAX(a1,a2);
    }
  }
  ss0 = MIN(ss[0],ss[1]);
  ss0 = MIN(ss0,ss[2]);
  for (i=0;i<NDIM;i++) 
    xs1[i] = xs0[i] + ss0*nmdr[i];
  fse[1] = impl_func(xs1,par);
  
  delnew = del0;
  not_conv = 1;
  iter = k = ipt = 0;
  while (not_conv  && iter < MAX_ITER_MINI) {              /* - */
    sign_change = vofi_get_segment_min(impl_func,par,xs0,nmdr,fse,xfs_pt,
                                       ss0,f2pos);
    for (i=0;i<NDIM;i++)                    
      xs0[i] = xfs_pt->xval[i];
    fse[0] = xfs_pt->fval;                   
    fs0 = f2pos*fse[0];                   
    if (sign_change)
      not_conv = 0;                                   /* - */
    else {        
      for (i=0;i<NDIM;i++) {                          /* - */
        x1f[i] = xs0[i] + dh*dir1[i];
        x1b[i] = xs0[i] - dh*dir1[i];
        x2f[i] = xs0[i] + dh*dir2[i];
        x2b[i] = xs0[i] - dh*dir2[i];
      }
      ss0 = xfs_pt->sval;

      f1f = f2pos*impl_func(x1f,par);
      f1b = f2pos*impl_func(x1b,par);
      f2f = f2pos*impl_func(x2f,par);
      f2b = f2pos*impl_func(x2b,par);
      df1 = -0.5*(f1f-f1b)/dh;
      df2 = -0.5*(f2f-f2b)/dh;
      d2f1 = (f1f+f1b-2.*fs0)/(dh*dh);
      d2f2 = (f2f+f2b-2.*fs0)/(dh*dh);
      if (d2f1 <= 0. || d2f2 <= 0.) 
        d2f1 = d2f2 = 1.;
      delold = delnew;
      delmid = delnew = 0.;
      for (i=0;i<=2;i++) {
        res[i] = rs0[i] + df1*dir1[i] + df2*dir2[i];
        delmid += res[i]*pcrs[i];
        hes[i] = hs0[i] + d2f1*dir1[i] + d2f2*dir2[i];
        pcrs[i] = res[i]/hes[i]; 
        delnew += res[i]*pcrs[i];               
      }

      beta = (delnew-delmid)/delold;      
      k++;          

      if (k == 2 || beta <= 0.) {
        beta = 0.;
        k = 0;
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
        a1 = (x0[i] - xs0[i])/(d1*d2);
        a2 = (x0[i] + h0[i] - xs0[i])/(d1*d2);
        ss[i] = MAX(a1,a2);
      }
      ss1 = MIN(ss[0],ss[1]);
      ss1 = MIN(ss1,ss[2]);             
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
      iter ++;
    }
  }
  
  return sign_change;
} 
