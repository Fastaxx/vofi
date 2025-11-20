/**
 * @file      getzero.c
 * @authors   Andrea Chierici, Leonardo Chirco, Vincent Le Chenadec, 
 *            Ruben Scardovelli, Philip Yecko and Stéphane Zaleski 
 * @date      April 15, 2021
 * @brief     Root finding routine along an oriented segment
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
 * Root finding routine along an oriented segment                             *
 * INPUT:  pointer to the implicit function, pointer to function parameters   *
 *         par, starting point x0, segment direction dir, s0, sign to have    *
 *         a negative function value in x0                                    *
 * OUTPUT: length of segment section where the function is negative           *
 * -------------------------------------------------------------------------- *
 * DESCRIPTION of array s0:                                                   *
 * s0[0]: segment length,       s0[1]: starting position ss,                  *
 * s0[2]: function value at ss, s0[3]: function derivative at ss along dir    *
 * -------------------------------------------------------------------------- */
vofi_real vofi_get_segment_zero(integrand impl_func,vofi_void_cptr par,
                           vofi_creal x0[],vofi_creal dir[],vofi_real s0[],
                           vofi_cint f_sign)
{
  vofi_int not_conv,i,iter,gensec;
  vofi_real xs[NDIM],sl,sr,ss,fs,dss,dsold,ds2,fps,sz,s1,s2,f1,f2;
  vofi_real fl=-EPS_SEGM,fr=EPS_SEGM;
  vofi_real sv[NDIM],fv[NDIM];
  
  sl = 0.; sr = s0[0];
  not_conv = 1;
  dsold = dss = s0[0];
  ss = s0[1];
  fs = f_sign*s0[2];                                            /* - */
  fps= f_sign*s0[3];
  sv[0] = sv[1] = sv[2] = ss;
  fv[0] = fs; fv[1] = fps; fv[2] = 0.;
  gensec = 0;
  if (fs < 0.0) {
    sl = ss; fl = fs;
  }
  else if (fs > 0.0) {
    sr = ss; fr = fs;
  }
  else
    not_conv = 0;                                                 /* - */
  iter = 0;
  while (not_conv  && iter < MAX_ITER_ROOT) {              /* - */

    if ( ((ss-sr)*fps-fs)*((ss-sl)*fps-fs) > 0.0 ||
         fabs(2.*fs) > fabs(dsold*fps) ) {                 /* - */
      dsold = dss;
      dss = 0.5*(sr - sl);
      ss = sl + dss;
      gensec = 0;
    }
    else {                                                    /* - */
      dsold = dss;
      dss = fs/fps;
      ss = ss - dss;
    }
    iter++;
    if (fabs(dss) < EPS_ROOT) {                     /* - */
      not_conv = 0; s0[3] = f_sign*fps;
    }
    
    if (not_conv) {                          /* - */
      for (i=0; i<NDIM; i++)
        xs[i] = x0[i] + ss*dir[i];
      fs = f_sign*impl_func(xs,par);                       /* - */
      fps = (fs-fv[0])/(ss-sv[2]);
      if (fs < 0.0) {                                    /* - */
        sl = ss; fl = fs;
      }
      else if (fs > 0.0) {
        sr = ss; fr = fs;
      }
      else {
       not_conv = 0; s0[3] = f_sign*fps;
      }
      sv[0] = sv[1]; sv[1] = sv[2]; sv[2] = ss;
      ds2 = sv[2]-sv[0];
      if (gensec > 0 && fabs(ds2) > EPS_ROOT)  
        fv[2] = (fps - fv[1])/ds2;          /* - */
      else
        fv[2] = 0.;
      fv[0] = fs; fv[1] = fps;
      gensec = 1;
      fps = fv[1] + fv[2]*(sv[2]-sv[1]);     /* - */
    }
  }

  if (!not_conv)                                               /* - */
    sz = f_sign*ss + 0.5*(1-f_sign)*s0[0];     /* - */
  else {                                              
    s1 = 0.;
    f1 = f_sign*impl_func(x0,par);
    s2 = s0[0];
    for (i=0; i<NDIM; i++)
      xs[i] = x0[i] + s2*dir[i];
    f2 = f_sign*impl_func(xs,par);   
    if (f1*f2 <= 0.) {
      if (sl > 0.) {
        s1 = sl; f1 = fl;
      }
      if (sr < s2) {
        s2 = sr; f2 = fr;
      }
      ss = s1 - f1*(s2-s1)/(f2-f1);
      sz = f_sign*ss + 0.5*(1-f_sign)*s0[0];
      s0[3] = f_sign*fps;     
    }
    else {
      printf("Warning: getzero bracket check failed (f1*f2 > 0).\n");
      printf("f1: %17.10e f2: %17.10e -- returning fallback.\n",f1,f2);
      /* fallback: linear interpolation between endpoints */
      if (f2 - f1 != 0.0) {
        ss = -f1*(s2 - s1)/(f2 - f1);
        sz = f_sign*ss + 0.5*(1-f_sign)*s0[0];
      } else {
        /* if values are equal (rare), return midpoint */
        sz = 0.5*s0[0];
      }
    }
  }
  
  return sz;
}
