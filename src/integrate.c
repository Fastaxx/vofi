/**
 * @file      integrate.c
 * @authors   Andrea Chierici, Vincent Le Chenadec, Ruben Scardovelli, 
 *            Philip Yecko and Stéphane Zaleski 
 * @date      February 15, 2021
 * @brief     Routines to compute the cut area in 2D (1D G-L integration)
 *            and the cut volume in 3D (2D G-L integration) 
 * @version   Vofi 2.0
 * @copyright GPLv3 license
  **/

#include "vofi_stddecl.h"
#include "vofi_GL_nodes.h"
#include "vofi_GL_weights.h"

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
 * DESCRIPTION OF 2D ROUTINE:                                                 *
 * compute  in each full/cut rectangle the area occupied by the reference     *
 * phase and, if required, the centroid coordinates                           *
 * INPUT:  pointer to the implicit function, arrays of function parameters    * 
 *         par, coordinates of minimum vertex x0, cell edges h0, secondary    *
 *         side subdivisions base, ordered coordinate directions pdir and     *
 *         sdir, len_data xhp (see vofi_stddecl.h), centroid coordinates,     *
 *         centroid flag ncen, user's number of points npt, number of         *
 *         subdivisions nsub, tentative number of integration points nptmp,   *
 *         sectors type nsect, heights direction ndire                        * 
 * OUTPUT: normalized area, centroid coordinates                              *
 * FUNCTIONS:                                                                 * 
 * vofi_real vofi_get_area: 2D problem, 1D integration                        *
 * -------------------------------------------------------------------------- */
vofi_real vofi_get_area(integrand impl_func,vofi_void_cptr par,vofi_creal x0[],
                        vofi_creal h0[],vofi_creal base[],vofi_creal pdir[],
                        vofi_creal sdir[],len_data xhp[],vofi_real centroid[],
                        vofi_cint ncen,vofi_cint npt[],vofi_cint nsub,
                        vofi_cint nptmp,vofi_int nsect[],vofi_int ndire[])
{
  vofi_int i,j,k,ns,npts,f_sign,it0;
  vofi_real x1[NDIM],x20[NDIM],x21[NDIM],s0[4],fse[NSE];
  vofi_real hp,hs,hm,mdpt,area,ds,tmp;
  vofi_real xp,xs,al,quada,quadp,quads,ratio;
  vofi_real a1,a2,b1,b2,dxm1,dxp1,dxm2,dxp2;
  vofi_creal *ptw, *ptx;

  area = hp = hs = xp = xs = 0.;
  it0 = 0;
  for (i=0;i<NDIM;i++) { 
    x1[i] = x0[i] + pdir[i]*h0[i];
    hp += pdir[i]*h0[i];
    hs += sdir[i]*h0[i];
  }
  hm = MAX(h0[0],h0[1]);
  hm = MAX(hm,h0[2]);
  for (ns=1;ns<=nsub;ns++) {                     /* - */
    ds = base[ns] - base[ns-1];
    mdpt = 0.5*(base[ns] + base[ns-1]);
    if (nsect[ns-1] > 0) {                                 /* - */
      al = ds*hp;
      area += al;
      if (ncen > 0) {
        xp += 0.5*hp*al;
        xs += mdpt*al;
      }
    }
    else if (nsect[ns-1] < 0) {      /* - */
      npts = (vofi_int) (18.*ds/hm); /* - */
      npts = MIN(20,npts+3);                               /* - */
      npts = MIN(nptmp,npts);                            /* - */
      if (npt[1] >= 3 && npt[1] <= 20)
        npts = MIN(npt[1],npts);                /* - */
      if (npt[0] >= 3 && npt[0] <= 20)
        npts = MAX(npt[0],npts);       /* - */
      xhp[it0].np0 = npts;
      f_sign = ndire[ns-1];
      xhp[it0].f_sign = f_sign;
      j = npts - 3;
      ptx = (double *) csipt[j];
      ptw = (double *) wgtpt[j];
      
      quada = quadp = quads = a1 = a2 = b1 = b2 = 0.;
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
        quada += (*ptw)*xhp[it0].ht0[k];
        if (ncen > 0) {
          quadp += (*ptw)*0.5*xhp[it0].ht0[k]*xhp[it0].ht0[k];
          quads += (*ptw)*xhp[it0].ht0[k]*xhp[it0].xt0[k];
        }
        ptx++;
        ptw++;
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
      quada = 0.5*ds*quada;
      area += quada; 
      if (ncen > 0) {
        quadp = 0.5*ds*quadp/quada;
        quads = 0.5*ds*quads/quada;
        if (f_sign < 0)
          quadp = hp - quadp;
        xp += quadp*quada;
        xs += quads*quada;
      }
      it0++;
    }
  }
  centroid[0] = xp;
  centroid[1] = xs;

  return area;
}

/* -------------------------------------------------------------------------- *
 * DESCRIPTION OF 3D ROUTINE:                                                 *
 * compute  in each full/cut cuboid the volume occupied by the reference      *
 * phase and, if required, the centroid coordinates                           *
 * INPUT:  pointer to the implicit function, arrays of function parameters    * 
 *         par, coordinates of minimum vertex x0, cell edges h0, tertiary     *
 *         side subdivisions base_ext, ordered coordinate directions pdir,    *
 *         sdir and tdir, centroid coordinates and interface area centroid,   *
 *         centroid and surface flags nex, user's number of points npt,       *
 *         number of external subdivisions nsub_ext, tentative number of      * 
 *         integration points nptmp, printing flags nvis                      *
 * OUTPUT: normalized volume, centroid coordinates and interface area         *
 * FUNCTIONS:                                                                 * 
 * vofi_real vofi_get_volume: 3D problem, 2D integration                      *
 * -------------------------------------------------------------------------- */
vofi_real vofi_get_volume(integrand impl_func,vofi_void_cptr par,vofi_creal x0[],
                   vofi_creal h0[],vofi_creal base_ext[],vofi_creal pdir[],
                   vofi_creal sdir[],vofi_creal tdir[],vofi_real centroid[],
                   vofi_cint nex[],vofi_cint npt[],vofi_cint nsub_ext,
                   vofi_cint nptmp,vofi_cint nvis[])
{
  vofi_int i,j,nt,k,nexpt,sect_hexa,nsub_int,nintmp;
  vofi_int nsect[NSEG],ndire[NSEG],nptin[2];
  vofi_real x1[NDIM],cent_2D[NDIM],base_int[NSEG],xmidt[NGLM+2];
  vofi_real volume,vol,hp,hs,ht,xp,xs,xt,hm;
  vofi_real dt,mdpt,xit,area,quadv,quadp,quads,quadt,surfer;
  vofi_creal *ptw_ext,*ptx_ext;
  len_data xhpn[2],xhpo[2];
  min_data xfs;
  
  volume = surfer = hp = hs = ht = xp = xs = xt = 0.;
  for (i=0;i<NDIM;i++) { 
    hp += pdir[i]*h0[i];
    hs += sdir[i]*h0[i];
    ht += tdir[i]*h0[i];
  }
  hm = MAX(h0[0],h0[1]);
  hm = MAX(hm,h0[2]);
 
  for (nt=1;nt<=nsub_ext;nt++) {                    /* - */
    dt = base_ext[nt] - base_ext[nt-1];        
    mdpt = 0.5*(base_ext[nt] + base_ext[nt-1]);
    for (i=0;i<NDIM;i++)  
      x1[i] = x0[i] + tdir[i]*mdpt;
    for (i=0;i<NDIM;i++)
      xfs.isc[i] = 0;
    sect_hexa = vofi_check_plane(impl_func,par,x1,h0,&xfs,base_int,pdir,sdir,
                                 nsect,ndire);
    
    if (!sect_hexa) {                                   /* - */
      if (nsect[0] == 1) {
        vol = dt*hs*hp;
        volume += vol;
        if (nex[0] > 0) {
          xp += 0.5*hp*vol;
          xs += 0.5*hs*vol;
          xt += mdpt*vol;
        }
      }
      else
        ;
    }
    else {                     /* - */
      nexpt = (vofi_int) (18.*dt/hm); /* - */
      nexpt = MIN(20,nexpt+3);          /* - */
      if (npt[3] >= 3 && npt[3] <= 20)
        nexpt = MIN(npt[3],nexpt);                        /* - */
      if (npt[2] >= 3 && npt[2] <= 20)
        nexpt = MAX(npt[2],nexpt);               /* - */
      j = nexpt - 3;
      ptx_ext = (double *) csipt[j];
      ptw_ext = (double *) wgtpt[j];

      quadv = quadp = quads = quadt = 0.;
      xmidt[0] = base_ext[nt-1];
      xmidt[nexpt+1] = base_ext[nt];
      for (k=1;k<=nexpt;k++) {
        xit = mdpt + 0.5*dt*(*ptx_ext);
        xmidt[k] = xit;
        for (i=0;i<NDIM;i++) 
          x1[i] = x0[i] + tdir[i]*xit;
        nsub_int = vofi_get_limits_inner_2D(impl_func,par,x1,h0,&xfs,base_int,
                                            pdir,sdir,nsect,ndire,sect_hexa);
        xhpn[0].np0 = xhpn[1].np0 = 0;
        area = vofi_get_area(impl_func,par,x1,h0,base_int,pdir,sdir,xhpn,
                             cent_2D,nex[0],npt,nsub_int,nptmp,nsect,ndire);
        if (nvis[0] > 0)
          tecplot_heights(x1,h0,pdir,sdir,xhpn);
        if (nex[1] > 0) {      
          vofi_end_points(impl_func,par,x1,h0,pdir,sdir,xhpn);
          if (k == 1) {
            for (i=0;i<NDIM;i++) 
              x1[i] = x0[i] + tdir[i]*xmidt[0];
            nintmp = vofi_get_limits_edge_2D(impl_func,par,x1,h0,&xfs,
                                             base_int,pdir,sdir,nsub_int);
            nptin[0] = xhpn[0].np0; nptin[1] = xhpn[1].np0; 
            xhpo[0].np0 = xhpo[1].np0 = 0;
            vofi_edge_points(impl_func,par,x1,h0,base_int,pdir,sdir,xhpo,
                             nptin,nintmp,nsect,ndire);
            vofi_end_points(impl_func,par,x1,h0,pdir,sdir,xhpo);
          }
          else if (k > 1 && k < nexpt) {
            surfer += vofi_interface_surface(impl_func,par,x0,h0,xmidt,pdir,
                                             sdir,tdir,xhpn,xhpo,k,nexpt,nvis[1]);
            xhpo[0] = xhpn[0]; xhpo[1] = xhpn[1];
          }
          else {
            for (i=0;i<NDIM;i++) 
              x1[i] = x0[i] + tdir[i]*xmidt[nexpt+1];
            nintmp = vofi_get_limits_edge_2D(impl_func,par,x1,h0,&xfs,
                                             base_int,pdir,sdir,nsub_int);
            nptin[0] = xhpo[0].np0; nptin[1] = xhpo[1].np0; 
            xhpn[0].np0 = xhpn[1].np0 = 0;
            vofi_edge_points(impl_func,par,x1,h0,base_int,pdir,sdir,xhpn,
                             nptin,nintmp,nsect,ndire);
            vofi_end_points(impl_func,par,x1,h0,pdir,sdir,xhpn);
            surfer += vofi_interface_surface(impl_func,par,x0,h0,xmidt,pdir,  
                                     sdir,tdir,xhpn,xhpo,k+1,nexpt,nvis[1]);
          }
        }
        quadv += (*ptw_ext)*area;
        quadp += (*ptw_ext)*cent_2D[0];
        quads += (*ptw_ext)*cent_2D[1];
        quadt += (*ptw_ext)*area*xit;
        ptx_ext++;
        ptw_ext++;
      }
      quadv = 0.5*dt*quadv;
      volume += quadv;
      if (nex[0] > 0) {
        xp += 0.5*dt*quadp;
        xs += 0.5*dt*quads;
        xt += 0.5*dt*quadt;
      }
    }
  }
  centroid[0] = xp;
  centroid[1] = xs;
  centroid[2] = xt;
  centroid[3] = surfer;
  
  return volume;
}
