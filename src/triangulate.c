/**
 * @file      triangulate.c
 * @authors   Andrea Chierici, Leonardo Chirco, Vincent Le Chenadec, 
 *            Ruben Scardovelli, Philip Yecko and Stéphane Zaleski 
 * @date      April 15, 2021
 * @brief     Routines to approximate the interface with triangles in
 *            a cut cell, compute their area and print them, if required,
 *            in Tecplot data file triangles.dat  
 * @version   Vofi 2.0
 * @copyright GPLv3 license
 **/

#include "vofi_stddecl.h"
#include "vofi_GL_nodes.h"

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
 * compute vertices of all triangles in (t,s,p) coordinates, their area and   *
 * print them in (x,y,z) coordinates in Tecplot data file triangles.dat       *
 * INPUT:  pointer to the implicit function, arrays of function parameters    *
 *         par, coordinates of minimum vertex x0, cell edges h0, tertiary     *
 *         coordinates xt, ordered coordinate directions pdir, sdir and tdir, *
 *         len_data xhpn and xhpo (see vofi_stddecl.h), plane index k, number *
 *         of planes nexpt, printing flag ipf                                 *
 * OUTPUT: triangles area                                                     *
 * FUNCTIONS:                                                                 *
 * vofi_real vofi_interface_surface: compute triangles from height points     *
 * vofi_real vofi_triarea: calculate triangle area                            *
 * -------------------------------------------------------------------------- */
vofi_real vofi_interface_surface(integrand impl_func,vofi_void_cptr par,
                         vofi_creal x0[],vofi_creal h0[],vofi_creal xt[],
                         vofi_creal pdir[],vofi_creal sdir[],vofi_creal tdir[],
                         len_data xhpn[],len_data xhpo[],vofi_cint k,
                         vofi_cint nexpt,vofi_cint ipf)
{
  vofi_int i,j,nsec,npn,npo,f_sign,it0,km;
  vofi_int djl,djc,djr,npa,nmin,nmax;
  vofi_real xa[NDIM],xb[NDIM],xc[NDIM],x1[NDIM],x2[NDIM];
  vofi_real hp,surfer,s0[4],sc,tc,hc,t1,t2,dxl,dxr,hsum,ratio;
  vofi_real *pts1,*pte1,*pts2,*pte2,*psa,*psb,*pth1,*pth2,*pha,*phb;
    
  surfer = hp = 0.;
  nsec = 0;
  km = k - 1;
  if (k == 2 || k > nexpt)
    km--;
  if (xhpn[1].np0 > 0)
    nsec = 2;
  else if (xhpn[0].np0 > 0)
    nsec = 1;
  for (i=0;i<NDIM;i++)  
    hp += pdir[i]*h0[i];    
  s0[0] = hp;
  for (it0=0; it0<nsec; it0++) {
    npn = xhpn[it0].np0;
    npo = xhpo[it0].np0;
    f_sign = xhpn[it0].f_sign;
    djl = djr = djc = -1;
    nmin = MIN(npn,npo);
    nmax = MAX(npn,npo);
    if (npn >= npo) {
      pts1 = &xhpn[it0].xt0[0]; pte1 = &xhpn[it0].xt0[npn-1];
      pts2 = &xhpo[it0].xt0[0]; pte2 = &xhpo[it0].xt0[npo-1];
      pth1 = &xhpn[it0].ht0[0]; pth2 = &xhpo[it0].ht0[0];
      t1 = xt[k]; t2 = xt[km];
    }
    else {
      pts1 = &xhpo[it0].xt0[0]; pte1 = &xhpo[it0].xt0[npo-1];
      pts2 = &xhpn[it0].xt0[0]; pte2 = &xhpn[it0].xt0[npn-1];
      pth1 = &xhpo[it0].ht0[0]; pth2 = &xhpn[it0].ht0[0];
      t2 = xt[k]; t1 = xt[km];
    }
    if (nmin == 1)
      djl = nmax-2;
    else {
      djc = nmin-2;
      npa = nmax; psa = pts1; psb = pte1;
      dxl = (*pts2) - (*psa);
      dxr = (*psb) - (*pte2);
      while (npa > nmin) {
        if (dxr >= dxl) {
          djr++; psb--; dxr = (*psb) - (*pte2);
        }
        else {
          djl++; psa++; dxl = (*pts2) - (*psa);
        }
        npa--;
      }
    }

    /* - */
    for (j=0;j<=djl;j++) {
      xa[0] = t1; xa[1] = *pts1; xa[2] = *pth1;
      xb[0] = t2; xb[1] = *pts2; xb[2] = *pth2;      
      pts1++; pth1++;
      xc[0] = t1; xc[1] = *pts1; xc[2] = *pth1;
      surfer += vofi_triarea(xa,xb,xc); 
      if (ipf)
        tecplot_triangle(x0,pdir,sdir,tdir,xa,xb,xc,hp,f_sign);
    }
    
    /* - */
    tc = 0.5*(t1+t2);
    s0[3] = 0.5*(xhpo[it0].htp[1] + xhpn[it0].htp[1]);
    for (j=0;j<=djc;j++) {
      psa = pts1; psb = pts2; pha = pth1; phb = pth2;
      pts1++; pts2++; pth1++; pth2++; 
      sc = 0.25*((*psa) + (*psb) + (*pts1) + (*pts2));
      hsum = (*pha) + (*phb) + (*pth1) + (*pth2);
      s0[1] = 0.25*hsum;
      if (f_sign < 0)
        s0[1] = hp - s0[1];
      ratio = s0[1]/hp;
      if (ratio < NEAR_EDGE_RATIO)
        s0[1] = 0.;
      else if (ratio > (1.-NEAR_EDGE_RATIO))
        s0[1] = hp;   
      for (i=0;i<NDIM;i++) {
        x1[i] = x0[i] + tc*tdir[i] + sc*sdir[i];
        x2[i] = x1[i] + s0[1]*pdir[i];
      }
      s0[2] = impl_func(x2,par);
      hc = vofi_get_segment_zero(impl_func,par,x1,pdir,s0,f_sign);
      /* - */
      xa[0] = t1; xa[1] = *psa;  xa[2] = *pha;
      xb[0] = tc; xb[1] = sc;    xb[2] = hc;
      xc[0] = t1; xc[1] = *pts1; xc[2] = *pth1;
      surfer += vofi_triarea(xa,xb,xc);
      if (ipf)
        tecplot_triangle(x0,pdir,sdir,tdir,xa,xb,xc,hp,f_sign);
      xc[0] = t2; xc[1] = *psb;  xc[2] = *phb;
      surfer += vofi_triarea(xa,xb,xc); 
      if (ipf)
        tecplot_triangle(x0,pdir,sdir,tdir,xa,xb,xc,hp,f_sign);
      xa[0] = t2; xa[1] = *pts2; xa[2] = *pth2;
      surfer += vofi_triarea(xa,xb,xc); 
      if (ipf)
        tecplot_triangle(x0,pdir,sdir,tdir,xa,xb,xc,hp,f_sign);
      xc[0] = t1; xc[1] = *pts1;  xc[2] =*pth1;
      surfer += vofi_triarea(xa,xb,xc); 
      if (ipf)
        tecplot_triangle(x0,pdir,sdir,tdir,xa,xb,xc,hp,f_sign);
    }

    /* - */
    for (j=0;j<=djr;j++) {
      xa[0] = t1; xa[1] = *pts1; xa[2] = *pth1;
      xb[0] = t2; xb[1] = *pts2; xb[2] = *pth2;
      pts1++; pth1++;
      xc[0] = t1; xc[1] = *pts1; xc[2] = *pth1;
      surfer += vofi_triarea(xa,xb,xc); 
      if (ipf)
        tecplot_triangle(x0,pdir,sdir,tdir,xa,xb,xc,hp,f_sign);
    }
  }

  return surfer;
}

/* -------------------------------------------------------------------------- */
vofi_real vofi_triarea(vofi_creal x1[],vofi_creal x2[],vofi_creal x3[])
{
  vofi_int i;
  vofi_real u[NDIM],v[NDIM],area;

  for (i=0;i<NDIM;i++) {
    u[i] = x2[i] - x1[i]; v[i] = x3[i] - x1[i];
  }
  area = 0.5*sqrt( Sq(u[1]*v[2] - u[2]*v[1]) + Sq(u[2]*v[0] - u[0]*v[2]) +
                   Sq(u[0]*v[1] - u[1]*v[0]) );

  return area;
}

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * compute left and right endpoints that substitute first and last height     *
 * points on an integration plane                                             *
 * INPUT:  pointer to the implicit function, arrays of function parameters    *
 *         par, starting point along tertiary direction x0, cell edges h0,    *
 *         ordered coordinate directions pdir and sdir, len_data xhhp (see    *
 *         vofi_stddecl.h)                                                    *
 * OUTPUT: updated heights data in xhhp                                       *
 * FUNCTIONS:                                                                 *
 * void vofi_end_points                                                       *
 * -------------------------------------------------------------------------- */
void vofi_end_points(integrand impl_func,vofi_void_cptr par,
                     vofi_creal x0[],vofi_creal h0[],vofi_creal pdir[],
                     vofi_creal sdir[],len_data xhhp[])
{
  vofi_int i,j,npt,f_sign,it0,nseg,j0,j1,j2,j3;
  vofi_real hp,a1,a2,b1,b2,ratio;
  vofi_real dx1,dx2,dx12,dc1,dc2;
  vofi_real x20[NDIM],x21[NDIM],s0[4];

  nseg = 0;
  hp = 0.;
  if (xhhp[1].np0 > 0)
    nseg = 2;
  else if (xhhp[0].np0 > 0)
    nseg = 1;
  for (i=0;i<NDIM;i++)  
    hp += pdir[i]*h0[i];    
  s0[0] = hp;
  for (it0=0; it0<nseg; it0++) {
    npt = xhhp[it0].np0;
    if (npt > 1) {
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
        xhhp[it0].ht0[j0] = vofi_get_segment_zero(impl_func,par,x20,pdir,
        					  s0,f_sign);
        xhhp[it0].htp[j0] = s0[3];
        j0 = npt+1; j1 = npt; j2 = npt-1; j3 = npt-2;
      }
      for (i=2;i<=npt-1;i++) {
        xhhp[it0].xt0[i-1] = xhhp[it0].xt0[i];	
        xhhp[it0].ht0[i-1] = xhhp[it0].ht0[i];
      }
      xhhp[it0].xt0[npt-1] = xhhp[it0].xt0[npt+1];	
      xhhp[it0].ht0[npt-1] = xhhp[it0].ht0[npt+1];
    }
    else
      xhhp[it0].ht0[0] = xhhp[it0].ht0[1];
  }
  
  return;
}

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * compute height points on edge plane                                        *
 * INPUT:  pointer to the implicit function, arrays of function parameters    *
 *         par, starting point along tertiary direction x0, cell edges h0,    *
 *         secondary side subdivisions base, ordered coordinate directions    *
 *         pdir and sdir, len_data xhp (see vofi_stddecl.h), integration      *
 *         points of next plane npt, number of subdivisions nsub, sectors     * 
 *         type nsect, heights direction ndire                                *
 * OUTPUT: heights data on edge plane in xhp                                  *
 * FUNCTIONS:                                                                 *
 * void vofi_edge_points                                                      *
 * -------------------------------------------------------------------------- */
void vofi_edge_points(integrand impl_func,vofi_void_cptr par,vofi_creal x0[],
                     vofi_creal h0[],vofi_creal base[],vofi_creal pdir[],
                     vofi_creal sdir[],len_data xhp[],vofi_cint npt[],
                     vofi_cint nsub,vofi_int nsect[],vofi_int ndire[])
{
  vofi_int i,j,k,ns,npts,f_sign,it0;
  vofi_real x1[NDIM],x20[NDIM],x21[NDIM],s0[4],fse[NSE];
  vofi_real hp,hs,mdpt,ds,tmp,ratio;
  vofi_real a1,a2,b1,b2,dxm1,dxp1,dxm2,dxp2;
  vofi_creal *ptx;

  hp = hs = 0.;
  it0 = 0;
  for (i=0;i<NDIM;i++) { 
    x1[i] = x0[i] + pdir[i]*h0[i];
    hp += pdir[i]*h0[i];
    hs += sdir[i]*h0[i];
  }
  for (ns=1;ns<=nsub;ns++) {                     /* - */
    if (nsect[ns-1] < 0)     {                 /* - */
      ds = base[ns] - base[ns-1];
      mdpt = 0.5*(base[ns] + base[ns-1]);
      if (ds < 2.*EPS_LOC)
        npts = 1;
      else
        npts = npt[it0];
      xhp[it0].np0 = npts;
      f_sign = ndire[ns-1];
      xhp[it0].f_sign = f_sign;
      j = MAX(0,npts - 3);
      ptx = (double *) csipt[j];
      
      a1 = a2 = b1 = b2 = 0.;
      xhp[it0].ht0[0] = xhp[it0].htp[0] = 0.;
      xhp[it0].xt0[0] = base[ns-1];
      xhp[it0].xt0[npts+1] = base[ns];      
      xhp[it0].xt0[1] = mdpt + 0.5*ds*(*ptx);
      for (i=0;i<NDIM;i++) {
        tmp = sdir[i]*xhp[it0].xt0[1];
        x20[i] = x0[i] + tmp;
        x21[i] = x1[i] + tmp;
      }    
      fse[0] = impl_func(x20,par);
      fse[1] = impl_func(x21,par);
      s0[0] = hp;                                  /* - */
      if ( fabs(fse[0]) < fabs(fse[1]) ) {
        s0[1] = 0.; s0[2] = fse[0];
      }
      else { 
        s0[1] = hp; s0[2] = fse[1];
      }
      s0[3] = (fse[1]-fse[0])/hp;
      for (k=1;k<=npts;k++) {            /* - */
        xhp[it0].ht0[k] = vofi_get_segment_zero(impl_func,par,x20,pdir,s0,
						f_sign);
        xhp[it0].htp[k] = s0[3];
        ptx++;
        if (k < npts) {
          a2 = a1;
          b2 = b1;
          xhp[it0].xt0[k+1] = mdpt + 0.5*ds*(*ptx);
          s0[1] = xhp[it0].ht0[k];
          s0[3] = xhp[it0].htp[k];
          if (k > 1) {
            dxm1 = (xhp[it0].xt0[k] - xhp[it0].xt0[k-1]);
            dxp1 = (xhp[it0].xt0[k+1] - xhp[it0].xt0[k]);
            a1 = (xhp[it0].ht0[k] - xhp[it0].ht0[k-1])/dxm1;
            s0[1] += a1*dxp1;
            b1 = (xhp[it0].htp[k] - xhp[it0].htp[k-1])/dxm1;
            s0[3] += b1*dxp1;
            if (k > 2) {
              dxm2 = (xhp[it0].xt0[k] - xhp[it0].xt0[k-2]);
              dxp2 = (xhp[it0].xt0[k+1] - xhp[it0].xt0[k-1]);
              s0[1] += (a1-a2)*dxp1*dxp2/dxm2;
              s0[3] += (b1-b2)*dxp1*dxp2/dxm2;
            }
          }
          if (f_sign < 0)
            s0[1] = hp - s0[1];
          ratio = s0[1]/hp;
          if (ratio < NEAR_EDGE_RATIO)
            s0[1] = 0.;
          else if (ratio > (1.-NEAR_EDGE_RATIO))
            s0[1] = hp;   
          for (i=0;i<NDIM;i++) {
            x20[i] = x0[i] + sdir[i]*xhp[it0].xt0[k+1];
            x21[i] = x20[i] + pdir[i]*s0[1];
          }
          s0[2] = impl_func(x21,par);
        }
      }
      it0++;
    }
  }

  return;
}
