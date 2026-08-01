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
 *         (1, 2, 3 or 4)                                                     *
 * OUTPUT: length/area/volume fraction cc, centroid coordinates and           *
 *         interface count/length/area xex                                    *
 * -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- *
 * vofi_get_cc: the historical entry point, unchanged.                        *
 * vofi_get_cc_gam: the same, plus the INTERFACE CENTROID in xgam (3 reals,   *
 * ndim0 == 1, 2 and 3). The interface measure vofi already returns is a      *
 * count of points (1D), a sum of chords (2D) or of triangles (3D), so the    *
 * centroid handed back here is the centroid of that same point set /         *
 * polyline / polyhedral surface -- the identical quadrature, not a separate  *
 * approximation. Where the cell carries no interface, xgam is the cell       *
 * centre. Pass xgam = NULL (or nex[1] == 0) to skip the extra work.          *
 *                                                                           *
 * ndim0 == 4 orders all four directions once, sweeps the cell along udir and *
 * hands every hyperplane u = const to the 3D kernel IN THAT SAME FRAME, so   *
 * the nesting u -> t -> s -> root is a genuine 4D height-function            *
 * integration and the interface is one graph x_p = H(s,t,u) over the cell.   *
 * The measure returned in xex[4] is the 3-volume of that graph, integrated   *
 * on the same nodes -- the analogue of the polyline in 2D and the            *
 * triangulation in 3D, though computed as an integral rather than            *
 * reconstructed, so it is not a sum over facets.                             *
 *                                                                           *
 * ndim0 == 1 is the degenerate case and is worth stating plainly: on a       *
 * segment the wet region is an interval, so cc is the root position, the     *
 * centroid xex[0] is the interval midpoint, and the interface is a single    *
 * POINT -- measure xex[3] = 1 (a count, dimensionless, NOT a length) with    *
 * xgam[0] the point itself. Nothing is integrated; the whole content is the  *
 * root, which is why it is solved exactly rather than by a secant.           *
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
  vofi_int  i,icc,nsub,nxex,ngam;
  vofi_int nsect[NSEG],ndire[NSEG];
  vofi_real f04D[NVERH],f03D[NSE][NSE][NSE],f02D[NSE][NSE];
  vofi_real base[NSEG];
  vofi_real centroid[2*NDIM+1],x0[NDIM],h0l[NDIM],area,volume,cc;
  vofi_real pdir[NDIM]={0.,0.,0.,0.},sdir[NDIM]={0.,0.,0.,0.};
  vofi_real tdir[NDIM]={0.,0.,0.,0.},udir[NDIM]={0.,0.,0.,0.};
  min_data  xfsp[NXFS];
  len_data xhp[2];

  if (ndim0 < 1 || ndim0 > NDIM) {
    printf(" EXIT: wrong value of variable ndim0! \n");
    exit(1);
  }
  memset(xfsp,0,sizeof(xfsp));
  memset(centroid,0,sizeof(centroid));
  xhp[0].np0 = xhp[1].np0 = 0;

  /* Work on LOCAL copies padded with zeros: the caller only ever supplies
     ndim0 components, while every direction-driven routine inside the
     library loops over all NDIM of them.                                 */
  for (i=0;i<NDIM;i++) {
    x0[i]  = (i < ndim0) ? xin[i] : 0.;
    h0l[i] = (i < ndim0) ? h0[i]  : 0.;
  }
  /* xex holds ndim0 centroid components followed by the interface measure,
     so it is 4 reals in 1D/2D/3D (the historical layout) and 5 in 4D;
     xgam holds ndim0 components, historically 3.                         */
  nxex = (ndim0 == 4) ? 5 : 4;
  ngam = (ndim0 == 4) ? 4 : 3;
  if (xgam != NULL) {              /* the cell centre, unless an interface */
    for (i=0;i<ngam;i++)           /* is found and asked for below         */
      xgam[i] = 0.;
    for (i=0;i<ndim0;i++)
      xgam[i] = x0[i] + 0.5*h0l[i];
  }
  for (i=0;i<nxex;i++)
    xex[i] = 0.0;
  if (ndim0 == 2) {                                               /* - */
    icc = vofi_order_dirs_2D(impl_func,par,x0,h0l,pdir,sdir,f02D,&xfsp[0]);
    if (icc >= 0) {
      cc = (vofi_real) icc;
      if (icc > 0 && nex[0] > 0) {
        for (i=0;i<NSE;i++)
          xex[i] = x0[i] + 0.5*h0l[i];
      }
      return cc;
    }
    nsub = vofi_get_limits_2D(impl_func,par,x0,h0l,f02D,xfsp[0],base,
			      pdir,sdir,nsect,ndire);
    area = vofi_get_area(impl_func,par,x0,h0l,base,pdir,sdir,xhp,centroid,
			 nex[0],npt,nsub,xfsp[0].ipt,nsect,ndire,NULL);
    cc = area/(h0l[0]*h0l[1]);
    if (nvis[0] > 0)
      tecplot_heights(x0,h0l,pdir,sdir,xhp);
    if (nex[0] > 0 && area > 0.) {
      centroid[0] = centroid[0]/area;
      centroid[1] = centroid[1]/area;
      centroid[2] = 0.;
      for (i=0;i<NSE;i++)
        xex[i] = x0[i] + centroid[0]*pdir[i] + centroid[1]*sdir[i];
    }
    if (nex[1] > 0) {
      vofi_real scent[NDIM],*psc;
      for (i=0;i<NDIM;i++)
        scent[i] = 0.;
      psc = (xgam != NULL) ? scent : NULL;
      xex[3] = vofi_interface_length(impl_func,par,x0,h0l,pdir,sdir,xhp,psc,
                                     nvis[1]);
      if (psc != NULL && xex[3] > 0.)
        for (i=0;i<NSE;i++)
          xgam[i] = x0[i] + (scent[0]/xex[3])*pdir[i] +
                            (scent[1]/xex[3])*sdir[i];
    }
  }
  else if (ndim0 == 3) {                                          /* - */
    icc = vofi_order_dirs_3D(impl_func,par,x0,h0l,pdir,sdir,tdir,f03D,xfsp);
    if (icc >= 0) {
      cc = (vofi_real) icc;
      if (icc > 0 && nex[0] > 0) {
        for (i=0;i<ndim0;i++)
          xex[i] = x0[i] + 0.5*h0l[i];
      }
      return cc;
    }
    nsub = vofi_get_limits_3D(impl_func,par,x0,h0l,f03D,xfsp,base,pdir,sdir,
                              tdir);
    volume = vofi_get_volume(impl_func,par,x0,h0l,base,pdir,sdir,tdir,centroid,
                             nex,npt,nsub,xfsp[4].ipt,nvis,NULL);
    cc = volume/(h0l[0]*h0l[1]*h0l[2]);

    if (nex[0] > 0 && volume > 0.) {
      centroid[0] = centroid[0]/volume;
      centroid[1] = centroid[1]/volume;
      centroid[2] = centroid[2]/volume;
      for (i=0;i<ndim0;i++)
        xex[i] = x0[i] + centroid[0]*pdir[i] + centroid[1]*sdir[i] +
                 centroid[2]*tdir[i];
    }
    if (nex[1] > 0) {
      xex[3] = centroid[3];
      if (xgam != NULL && centroid[3] > 0.)
        for (i=0;i<ndim0;i++)
          xgam[i] = x0[i] + (centroid[4]/centroid[3])*pdir[i] +
                            (centroid[5]/centroid[3])*sdir[i] +
                            (centroid[6]/centroid[3])*tdir[i];
    }
  }
  else if (ndim0 == 4) {                                          /* - */
    icc = vofi_order_dirs_4D(impl_func,par,x0,h0l,pdir,sdir,tdir,udir,f04D,
                             xfsp);
    if (icc >= 0) {
      cc = (vofi_real) icc;
      if (icc > 0 && nex[0] > 0) {
        for (i=0;i<ndim0;i++)
          xex[i] = x0[i] + 0.5*h0l[i];
      }
      return cc;
    }
    nsub = vofi_get_limits_4D(impl_func,par,x0,h0l,f04D,xfsp,base,pdir,sdir,
                              tdir,udir);
    volume = vofi_get_hypervolume(impl_func,par,x0,h0l,base,pdir,sdir,tdir,
                                  udir,centroid,nex,npt,nsub,xfsp[IXC4].ipt);
    cc = volume/(h0l[0]*h0l[1]*h0l[2]*h0l[3]);

    if (nex[0] > 0 && volume > 0.) {
      for (i=0;i<ndim0;i++)
        centroid[i] = centroid[i]/volume;
      for (i=0;i<ndim0;i++)
        xex[i] = x0[i] + centroid[0]*pdir[i] + centroid[1]*sdir[i] +
                 centroid[2]*tdir[i] + centroid[3]*udir[i];
    }
    if (nex[1] > 0) {
      xex[4] = centroid[4];
      if (xgam != NULL && centroid[4] > 0.)
        for (i=0;i<ndim0;i++)
          xgam[i] = x0[i] + (centroid[5]/centroid[4])*pdir[i] +
                            (centroid[6]/centroid[4])*sdir[i] +
                            (centroid[7]/centroid[4])*tdir[i] +
                            (centroid[8]/centroid[4])*udir[i];
    }
  }
  else if (ndim0 == 1) {                                          /* - */
    vofi_real dir[NDIM]={1.,0.,0.,0.},xs[NDIM]={0.,0.,0.,0.},s0[4];
    vofi_real f0,f1,hh,sz;
    vofi_int fsign;

    hh = h0l[0];
    xs[0] = x0[0];      f0 = impl_func(xs,par);
    xs[0] = x0[0] + hh; f1 = impl_func(xs,par);

    /* Full and empty are decided by the endpoint SIGNS alone, which also
       disposes of every exactly-zero endpoint: (0,+) and (+,0) are empty,
       (0,-) and (-,0) are full. A pair of roots interior to the segment
       is therefore reported full or empty -- a sub-grid feature this
       returns no information about, the same convention as
       CartesianGeometry.jl's 1D kernel. */
    if (f0 <= 0. && f1 <= 0.) {
      if (nex[0] > 0)
        xex[0] = x0[0] + 0.5*hh;
      return 1.;
    }
    if (f0 >= 0. && f1 >= 0.) {
      if (nex[0] > 0)
        xex[0] = x0[0] + 0.5*hh;
      return 0.;
    }

    /* One strict sign change: the root, and hence the wet length, is
       resolved by the SAME safeguarded solver the 2D/3D height functions
       use -- not by the secant. In 1D the wet volume is the root
       position, so a linearly-interpolated crossing is a first-order
       error in V itself; in 2D/3D it only perturbs an aperture. */
    fsign = (f0 < 0.) ? 1 : -1;
    s0[0] = hh;
    s0[1] = hh*f0/(f0 - f1);            /* secant crossing: the seed */
    xs[0] = x0[0] + s0[1];
    s0[2] = impl_func(xs,par);
    s0[3] = (f1 - f0)/hh;               /* slope estimate for Newton */
    sz = vofi_get_segment_zero(impl_func,par,x0,dir,s0,fsign);

    cc = sz/hh;                         /* sz IS the wet length */
    if (nex[0] > 0)
      xex[0] = (fsign > 0) ? x0[0] + 0.5*sz : x0[0] + hh - 0.5*sz;
    if (nex[1] > 0) {
      /* the interface is a POINT: its measure is the count 1, and its
         centroid is the point itself */
      xex[3] = 1.;
      if (xgam != NULL)
        xgam[0] = (fsign > 0) ? x0[0] + sz : x0[0] + hh - sz;
    }
  }
  else {                                                          /* - */
    printf(" EXIT: wrong value of variable ndim0! \n");
    exit(1);
  }

  return cc;
}
