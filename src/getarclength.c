/**
 * @file      getarclength.c
 * @authors   Andrea Chierici, Leonardo Chirco, Vincent Le Chenadec, 
 *            Ruben Scardovelli, Philip Yecko and Stéphane Zaleski 
 * @date      April 15, 2021
 * @brief     Routine to compute the length of the interface line inside 
 *            a cell (2D) 
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
 * compute the length of the interface line inside a cell (2D problem)        *
 * INPUT:  pointer to the implicit function, arrays of function parameters    * 
 *         par, coordinates of minimum vertex x0, cell edges h0, coordinate   * 
 *         directions pdir and sdir, len_data  xhhp (see vofi_stddecl.h),     *
 *         printing flag ipf                                                  *
 * OUTPUT: length of the interface line arc                                   *
 * FUNCTIONS:                                                                 *
 * vofi_real vofi_interface_length                                            *
 * -------------------------------------------------------------------------- */
vofi_real vofi_interface_length(integrand impl_func,vofi_void_cptr par,
                         vofi_creal x0[],vofi_creal h0[],vofi_creal pdir[],
                         vofi_creal sdir[],len_data xhhp[],vofi_cint ipf)
{
  FILE *fp;
  vofi_int i,j,k,npt,f_sign,it0,nseg,j0,j1,j2,j3;
  vofi_real hp,a1,a2,b1,b2,ratio,xm,hm,hpm,arc,d1,d2,hsum;
  vofi_real dx1,dx2,dx12,dc1,dc2,x0b,h0b,hpb,xc,hc,xl,hl,xr,hr;
  vofi_real x20[NDIM],x21[NDIM],s0[4];

  if (ipf == 1)
    fp = fopen("arcline.dat","a");    
  arc = hp = 0.;
  nseg = 0;
  if (xhhp[1].np0 > 0)
    nseg = 2;
  else if (xhhp[0].np0 > 0)
    nseg = 1;
  for (i=0;i<NDIM;i++)  
    hp += pdir[i]*h0[i];    
  s0[0] = hp;
  for (it0=0; it0<nseg; it0++) {
    npt = xhhp[it0].np0;
    f_sign = xhhp[it0].f_sign;
    j0 = 0; j1 = 1; j2 = 2; j3 = 3;
    /* - */
    for (j=0;j<NSE;j++) {
      dx1  = xhhp[it0].xt0[j1] - xhhp[it0].xt0[j2];
      dx2  = xhhp[it0].xt0[j2] - xhhp[it0].xt0[j3];
      dx12 = xhhp[it0].xt0[j1] - xhhp[it0].xt0[j3];
      dc1  = xhhp[it0].xt0[j0] - xhhp[it0].xt0[j1];

      dc2  = xhhp[it0].xt0[j0] - xhhp[it0].xt0[j2];
      a1  = (xhhp[it0].ht0[j1] - xhhp[it0].ht0[j2])/dx1;
      a2  = (xhhp[it0].ht0[j2] - xhhp[it0].ht0[j3])/dx2;
      b1  = (xhhp[it0].htp[j1] - xhhp[it0].htp[j2])/dx1;
      b2  = (xhhp[it0].htp[j2] - xhhp[it0].htp[j3])/dx2;
      s0[1] = xhhp[it0].ht0[j1] + a1*dc1 + (a1-a2)*dc1*dc2/dx12;
      s0[3] = xhhp[it0].htp[j1] + b1*dc1 + (b1-b2)*dc1*dc2/dx12;
      if (f_sign < 0)
        s0[1] = hp - s0[1];
      ratio = s0[1]/hp;
      if (ratio < NEAR_EDGE_RATIO)
        s0[1] = 0.;
      else if (ratio > (1.-NEAR_EDGE_RATIO))
        s0[1] = hp;   

      for (i=0;i<NDIM;i++) {
        x20[i] = x0[i] + sdir[i]*xhhp[it0].xt0[j0];
        x21[i] = x20[i] + pdir[i]*s0[1];
      }
      s0[2] = impl_func(x21,par);
      xhhp[it0].ht0[j0] = vofi_get_segment_zero(impl_func,par,
                                           x20,pdir,s0,f_sign);
      xhhp[it0].htp[j0] = s0[3];
      j0 = npt+1; j1 = npt; j2 = npt-1; j3 = npt-2;
    }
    
    /* - */
    x0b = xhhp[it0].xt0[0];
    h0b = xhhp[it0].ht0[0];
    hpb = xhhp[it0].htp[0];
    xm = 0.5*(x0b + xhhp[it0].xt0[1]);
    xl = x0b;
    hl = h0b;
    xr = xhhp[it0].xt0[1];
    hr = xhhp[it0].ht0[1];
    k = 0;
    if (ipf == 1) {
      fprintf(fp," ZONE I=1, J=%d, F=POINT \n",2*npt+3);
      tecplot_arcline(x0,pdir,sdir,xl,hl,hp,f_sign,fp);
    }      
    for (j=0;j<=npt;j++) {
      dx1  = x0b - xhhp[it0].xt0[k+1];
      dx2  = xhhp[it0].xt0[k+1] - xhhp[it0].xt0[k+2];
      dx12 = x0b - xhhp[it0].xt0[k+2];
      dc1  = xm  - x0b;
      dc2  = xm  - xhhp[it0].xt0[k+1];
      a1 = (h0b - xhhp[it0].ht0[k+1])/dx1;
      b1 = (hpb - xhhp[it0].htp[k+1])/dx1;
      a2 = (xhhp[it0].ht0[k+1] - xhhp[it0].ht0[k+2])/dx2;
      b2 = (xhhp[it0].htp[k+1] - xhhp[it0].htp[k+2])/dx2;
      s0[1] = h0b + a1*dc1 + (a1-a2)*dc1*dc2/dx12;
      s0[3] = hpb + b1*dc1 + (b1-b2)*dc1*dc2/dx12;

      if (f_sign < 0)
        s0[1] = hp - s0[1];
      ratio = s0[1]/hp;
      if (ratio < NEAR_EDGE_RATIO)
        s0[1] = 0.;
      else if (ratio > (1.-NEAR_EDGE_RATIO))
        s0[1] = hp;   
      for (i=0;i<NDIM;i++) {
        x20[i] = x0[i] + sdir[i]*xm;
        x21[i] = x20[i] + pdir[i]*s0[1];
      }
      s0[2] = impl_func(x21,par);
      hm = vofi_get_segment_zero(impl_func,par,x20,pdir,s0,f_sign);
      hpm = s0[3]; 
      xc = xm;
      hsum = hl + hr;
      hc = 0.5*hsum + sqrt(3.)/3.*(2.*hm - hsum);
      d1 = (xl-xc)*(xl-xc) + (hl-hc)*(hl-hc);
      d2 = (xr-xc)*(xr-xc) + (hr-hc)*(hr-hc);
      arc += sqrt(d1) + sqrt(d2);
      if (ipf == 1) {
        tecplot_arcline(x0,pdir,sdir,xc,hc,hp,f_sign,fp);   
        tecplot_arcline(x0,pdir,sdir,xr,hr,hp,f_sign,fp);
      }      
    
      x0b = xm;
      h0b = hm;
      hpb = hpm;
      k = MIN(j,npt-1);
      xl = xhhp[it0].xt0[k+1];
      hl = xhhp[it0].ht0[k+1];
      xr = xhhp[it0].xt0[k+2];
      hr = xhhp[it0].ht0[k+2];
      xm = 0.5*(xl+xr);
    }
  }
  if (ipf == 1)
    fclose(fp);
  
  return arc;
}
