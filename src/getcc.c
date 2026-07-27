/**
 * @file      getcc.c
 * @authors   Andrea Chierici, Leonardo Chirco, Vincent Le Chenadec, 
 *            Ruben Scardovelli, Philip Yecko and Stéphane Zaleski 
 * @date      April 15, 2021
 * @brief     Driver to compute the area/volume fraction in a given cell
 *            and, if required, the centroid, the interface length/area 
 *            and to print data in Tecplot files
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
 * driver to compute the area/volume fraction in a given cell and, if         *
 * required, the centroid, the interface length/area and to print data in     *
 * Tecplot files                                                              *
 * INPUT:  pointer to the implicit function, arrays of function parameters    *
 *         par, coordinates of minimum vertex xin, cell edges h0, integer     *
 *         flags to compute centroid and interface length/area nex, user's    *
 *         number of points npt, printing flags nvis, space dimensions ndim0  *
 * OUTPUT: area/volume fraction cc,  centroid coordinates and interface       *
 *         length/area xex                                                    *
 * -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- *
 * vofi_get_cc: the historical entry point, unchanged.                        *
 * vofi_get_cc_gam: the same, plus the INTERFACE CENTROID in xgam (3 reals,   *
 * both ndim0 == 2 and ndim0 == 3). The interface length/area vofi already    *
 * returns is a sum of chords (2D) or of triangles (3D), so the centroid      *
 * handed back here is the centroid of that same polyline/polyhedral surface  *
 * -- the identical quadrature, not a separate approximation. Where the cell  *
 * carries no interface, xgam is the cell centre. Pass xgam = NULL (or        *
 * nex[1] == 0) to skip the extra work.                                       *
 * -------------------------------------------------------------------------- */
vofi_real vofi_get_cc(integrand impl_func,vofi_void_cptr par,vofi_creal xin[],
                      vofi_creal h0[],vofi_real xex[],vofi_cint nex[],
                      vofi_cint npt[],vofi_cint nvis[],vofi_cint ndim0)
{
  return vofi_get_cc_gam(impl_func,par,xin,h0,xex,NULL,nex,npt,nvis,ndim0);
}

vofi_real vofi_get_cc_gam(integrand impl_func,vofi_void_cptr par,
                      vofi_creal xin[],vofi_creal h0[],vofi_real xex[],
                      vofi_real xgam[],vofi_cint nex[],
                      vofi_cint npt[],vofi_cint nvis[],vofi_cint ndim0)
{
  vofi_int  i,icc,nsub;
  vofi_int nsect[NSEG],ndire[NSEG];
  vofi_real f03D[NSE][NSE][NSE],f02D[NSE][NSE],base[NSEG];
  vofi_real centroid[NDIM+1+NDIM],x0[NDIM],area,volume,cc;
  vofi_real pdir[NDIM]={0.,0.,0.},sdir[NDIM]={0.,0.,0.},tdir[NDIM]={0.,0.,0.};
  min_data  xfsp[5]={{{0.,0.,0.},0.,0.,{0,0,0},0},{{0.,0.,0.},0.,0.,{0,0,0},0},
                     {{0.,0.,0.},0.,0.,{0,0,0},0},{{0.,0.,0.},0.,0.,{0,0,0},0},
                     {{0.,0.,0.},0.,0.,{0,0,0},0}}; 
  len_data xhp[2];
  
  xhp[0].np0 = xhp[1].np0 = 0;
  if (xgam != NULL) {              /* the cell centre, unless an interface */
    for (i=0;i<NDIM;i++)           /* is found and asked for below         */
      xgam[i] = 0.;
    for (i=0;i<ndim0;i++)
      xgam[i] = xin[i] + 0.5*h0[i];
  }
  for (i=0;i<=NDIM;i++)
    xex[i] = 0.0;
  if (ndim0 == 2) {                                               /* - */
    x0[0] = xin[0]; x0[1] = xin[1]; x0[2] = 0.;
    icc = vofi_order_dirs_2D(impl_func,par,x0,h0,pdir,sdir,f02D,&xfsp[0]);
    if (icc >= 0) {
      cc = (vofi_real) icc;
      if (icc > 0 && nex[0] > 0) {
        for (i=0;i<NSE;i++)
          xex[i] = x0[i] + 0.5*h0[i];
      }
      return cc;
    }
    nsub = vofi_get_limits_2D(impl_func,par,x0,h0,f02D,xfsp[0],base,
			      pdir,sdir,nsect,ndire);
    area = vofi_get_area(impl_func,par,x0,h0,base,pdir,sdir,xhp,centroid,
			 nex[0],npt,nsub,xfsp[0].ipt,nsect,ndire);
    cc = area/(h0[0]*h0[1]);
    if (nvis[0] > 0)
      tecplot_heights(x0,h0,pdir,sdir,xhp);
    if (nex[0] > 0 && area > 0.) {
      centroid[0] = centroid[0]/area;
      centroid[1] = centroid[1]/area;
      centroid[2] = 0.;
      for (i=0;i<NSE;i++)
        xex[i] = x0[i] + centroid[0]*pdir[i] + centroid[1]*sdir[i];
    }
    if (nex[1] > 0) {
      vofi_real scent[NDIM],*psc;
      scent[0] = scent[1] = scent[2] = 0.;
      psc = (xgam != NULL) ? scent : NULL;
      xex[3] = vofi_interface_length(impl_func,par,x0,h0,pdir,sdir,xhp,psc,
                                     nvis[1]);
      if (psc != NULL && xex[3] > 0.)
        for (i=0;i<NSE;i++)
          xgam[i] = x0[i] + (scent[0]/xex[3])*pdir[i] +
                            (scent[1]/xex[3])*sdir[i];
    }
  }
  else if (ndim0 == 3) {                                          /* - */
    x0[0] = xin[0]; x0[1] = xin[1]; x0[2] = xin[2];
    icc = vofi_order_dirs_3D(impl_func,par,x0,h0,pdir,sdir,tdir,f03D,xfsp);
    if (icc >= 0) {
      cc = (vofi_real) icc;
      if (icc > 0 && nex[0] > 0) {
        for (i=0;i<NDIM;i++)
          xex[i] = x0[i] + 0.5*h0[i];
      }
      return cc;
    }
    nsub = vofi_get_limits_3D(impl_func,par,x0,h0,f03D,xfsp,base,pdir,sdir,
                              tdir);
    volume = vofi_get_volume(impl_func,par,x0,h0,base,pdir,sdir,tdir,centroid,
                             nex,npt,nsub,xfsp[4].ipt,nvis);
    cc = volume/(h0[0]*h0[1]*h0[2]);

    if (nex[0] > 0 && volume > 0.) {
      centroid[0] = centroid[0]/volume;
      centroid[1] = centroid[1]/volume;
      centroid[2] = centroid[2]/volume;
      for (i=0;i<NDIM;i++)
        xex[i] = x0[i] + centroid[0]*pdir[i] + centroid[1]*sdir[i] + 
                 centroid[2]*tdir[i];
    }
    if (nex[1] > 0) {
      xex[3] = centroid[3];
      if (xgam != NULL && centroid[3] > 0.)
        for (i=0;i<NDIM;i++)
          xgam[i] = x0[i] + (centroid[4]/centroid[3])*pdir[i] +
                            (centroid[5]/centroid[3])*sdir[i] +
                            (centroid[6]/centroid[3])*tdir[i];
    }
  }
  else {                                                          /* - */
    printf(" EXIT: wrong value of variable ndim0! \n");
    exit(1);
  }

  return cc;
}
