/**
 * @file      getlimits.c
 * @authors   Andrea Chierici, Leonardo Chirco, Vincent Le Chenadec, 
 *            Ruben Scardovelli, Philip Yecko and Stéphane Zaleski 
 * @date      April 15, 2021
 * @brief     Routines to compute the limits of integration along the 
 *            secondary side in 2D, in 3D along the tertiary side for the 
 *            external integration and along the secondary side for the 
 *            internal one  
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
 * subdivide the side along the secondary/tertiary (2D/3D) direction          *
 * to define rectangles/right cuboids (2D/3D) with or without the interface   *
 * INPUT:  pointer to the implicit function, arrays of function parameters    * 
 *         par, coordinates of minimum vertex x0, cell edges h0, function     *
 *         value at vertices f0, min_data xfsp (see vofi_stddecl.h),          *
 *         ordered coordinate directions pdir and sdir (and tdir in 3D)       *
 * OUTPUT: number of subdivisions nsub, array of side subdivisions base,      * 
 *         sector type nsect, direction of height ndire                       *
 * FUNCTIONS:                                                                 *
 * vofi_int vofi_get_limits_2D: limits of integration in 2D                   *
 * vofi_int vofi_get_limits_3D: external limits of integration in 3D          *
 * vofi_int vofi_check_plane: check midplane of each sector                   *
 * vofi_int vofi_get_limits_inner_2D: limits of integration in each plane of  *
 * a cut sector (the fast alternative requires the number of subdivisions     *
 * from vofi_int vofi_check_plane)                                            *
 * vofi_int vofi_get_limits_edge_2D: limits of integration in edge planes of  *
 * a cut cuboid (only for triangulation, it requires the number of            *
 * subdivisions from the closest plane, the case where they are different has *
 * not been implemented: it should not occur!)                                *
 * void vofi_reorder: order limits in ascending order                         *
 * vofi_int vofi_rm_segs: remove zero-length segments                         *
 * void vofi_sector_new: in each sector nsect=1/0/-1: full/empty/cut; if      *
 * full/empty: ndire=0, if cut: ndire=1 and height bottom-up, or ndir=-1 and  *
 * height top-down; no extra function call is required                        *
 * void vofi_sector_old: old style: compute sign function in midsection, on   *
 * top and bottom of secondary sides; two extra function calls per sector     *
 * -------------------------------------------------------------------------- */
vofi_int vofi_get_limits_2D(integrand impl_func,vofi_void_cptr par,
                  vofi_creal x0[],vofi_creal h0[],vofi_real f0[][NSE],
                  min_data xfsp,vofi_real baser[],vofi_creal pdir[],
                  vofi_creal sdir[],vofi_int nsect[],vofi_int ndire[])
{
  vofi_int i,k,nsub,inters,down2up,atleast1,ncheck,iside,isect;
  vofi_int basei[NSEG],nbt[NSE]={0,0},sign_sect[NSE][NDIM]={{0,0,0},{0,0,0}};
  vofi_real x1[NDIM],fse[NSE],hs,fsum;
  
  baser[0] = hs = 0.;
  basei[0] = nsub = down2up = 1;
  atleast1 = 0;
  iside = isect = 0;
  
  for (i=0;i<NDIM;i++) 
    hs += sdir[i]*h0[i];
  for (k=0;k<NSE;k++) {           /* - */
    fse[0] = f0[k][0];
    fse[1] = f0[k][1];
    fsum = fse[0] + fse[1];
    if (xfsp.isc[k+1] == 0) {
      atleast1 = 1;
      nbt[k] = 1;
      iside = k;
      isect = 1 - k;
      if (fsum > 0.)
        sign_sect[k][0] = 1;
      else if (fsum < 0.)
        sign_sect[k][0] = -1;
      if ((k == 0 && sign_sect[k][0] == 1) ||
	  (k == 1 && sign_sect[k][0] == -1))
        down2up = -1;
    }
    else {
      for (i=0;i<NDIM;i++)
        x1[i] = x0[i] + k*pdir[i]*h0[i];
      inters = vofi_get_side_intersections(impl_func,par,fse,x1,xfsp,baser,
                                           sdir,hs,nsub,xfsp.isc[k+1]);
      nsub += inters;
      for (i=1;i<=inters;i++)
	basei[nsub-i] = 1;
      nbt[k] = inters + 1;
      if (inters == 1) {
        if (fse[0] < 0.) {
          sign_sect[k][0] = -1; sign_sect[k][1] = 1;
        }
        else {
          sign_sect[k][0] = 1; sign_sect[k][1] = -1;
        }
      }
      else {
        if (fsum > 0.) {
          sign_sect[k][0] = sign_sect[k][2] = 1; sign_sect[k][1] = -1;
        }
        else {
          sign_sect[k][0] = sign_sect[k][2] = -1; sign_sect[k][1] = 1;
        }
      }
    }
  }

  ncheck = atleast1*MAX(nbt[0],nbt[1]);
  baser[nsub] = hs;
  basei[nsub] = 1;
  
  /* - */
  vofi_reorder(baser,basei,nsub);

  /* - */
  nsub = vofi_rm_segs(baser,basei,nsub);

  /* - */
  if (ncheck == nsub) 
    vofi_sector_new(sign_sect,nsect,ndire,nsub,iside,isect,down2up);
  /* - */
  else 
    vofi_sector_old(impl_func,par,x0,h0,baser,pdir,sdir,nsect,ndire,nsub);
  
  return nsub;    
}

/* -------------------------------------------------------------------------- */
vofi_int vofi_get_limits_3D(integrand impl_func,vofi_void_cptr par,
              vofi_creal x0[],vofi_creal h0[],vofi_real f0[][NSE][NSE],
              min_data xfsp[],vofi_real baser[],vofi_creal pdir[],
              vofi_creal sdir[],vofi_creal tdir[])
{
  vofi_int  basei[NSEG],i,l0,l,m,n,nsub,inters,consi;
  vofi_real xp[NDIM],xs[NDIM],xt[NDIM],fse[NSE],ht,hs;
  min_data xfsl={{0.,0.,0.},0.,0.,{0,0,0},0};
  
  baser[0] = 0.;
  basei[0] = nsub = 1;
  hs = ht = 0.;
  for (i=0;i<NDIM;i++) {
    hs += sdir[i]*h0[i];
    ht += tdir[i]*h0[i];
  }
  
  /* - */
  for (m=0;m<NSE;m++) {
    for (l=0;l<NDIM;l++)
      xp[l] = x0[l] + m*pdir[l]*h0[l];
    /* - */
    for (n=0;n<NSE;n++) {
      l0 = 2*m + n;
      if (xfsp[l0].isc[1] != 0) {
        for (l=0;l<NDIM;l++)
          xs[l] = xp[l] + n*sdir[l]*h0[l];
        fse[0] = f0[m][n][0];
        fse[1] = f0[m][n][1];
        inters = vofi_get_side_intersections(impl_func,par,fse,xs,xfsp[l0],
                               baser,tdir,ht,nsub,xfsp[l0].isc[1]);
        nsub += inters;
	for (i=1;i<=inters;i++)
	  basei[nsub-i] = 1;
        for (l=0;l<NDIM;l++)
          xt[l] = xs[l] + baser[nsub-1]*tdir[l];	
        consi = vofi_check_line_consistency(impl_func,par,xt,sdir,hs,n,&xfsl);
        if (inters > 1 && consi == 0) {
          for (l=0;l<NDIM;l++)
            xt[l] = xs[l] + baser[nsub-2]*tdir[l];
          consi = vofi_check_line_consistency(impl_func,par,xt,sdir,hs,n,&xfsl);
        }
        if (consi > 0) {
          inters = vofi_get_ext_intersections(impl_func,par,xp,h0,xfsl,baser,
                                              sdir,tdir,nsub);
          nsub += inters;
	  for (i=1;i<=inters;i++)
	    basei[nsub-i] = 0;
        }
      }
    }
    if (xfsp[4].isc[m+1] != 0) {
      inters = vofi_get_ext_intersections(impl_func,par,xp,h0,xfsp[4],baser,
                                          sdir,tdir,nsub);
      nsub += inters;
      for (i=1;i<=inters;i++)
	basei[nsub-i] = 0;
    }
  }

  baser[nsub] = ht;
  basei[nsub] = 1;

  /* - */
  vofi_reorder(baser,basei,nsub);
  
  /* - */
  nsub = vofi_rm_segs(baser,basei,nsub);
  
  return nsub;    
}

/* -------------------------------------------------------------------------- */
vofi_int vofi_check_plane(integrand impl_func,vofi_void_cptr par,
                    vofi_creal x0[],vofi_creal h0[],min_data *xfs_pt,
                    vofi_real baser[],vofi_creal pdir[],vofi_creal sdir[],
                    vofi_int nsect[],vofi_int ndire[])
{
  vofi_int i,k,consi,f2pos,sign_change,nsub;
  vofi_int inters,nointer,down2up,atleast1,ncheck,iside,isect;
  vofi_int basei[NSEG],nbt[NSE]={0,0},sign_sect[NSE][NDIM]={{0,0,0},{0,0,0}};
  vofi_real x1[NDIM],x2[NDIM],fse[NSE],hs,fsum;
  min_data xfsl={{0.,0.,0.},0.,0.,{0,0,0},0};

  baser[0] = hs = 0.;
  basei[0] = nsub = down2up = 1;
  atleast1 = 0;
  iside = isect = 0;
  for (i=0;i<NDIM;i++) {
    hs += sdir[i]*h0[i];
    xfs_pt->isc[i] = 0;
  }
  
  for (k=0;k<NSE;k++) {           /* - */
    for (i=0;i<NDIM;i++) {
      x1[i] = x0[i] + k*pdir[i]*h0[i];
      x2[i] = x1[i] + sdir[i]*hs;
    }
    fse[0] = impl_func(x1,par);
    fse[1] = impl_func(x2,par);
    fsum = fse[0] + fse[1];

    if (fse[0]*fse[1] < 0.) {
      xfs_pt->isc[0]   = xfsl.isc[0]   = 1;
      xfs_pt->isc[k+1] = xfsl.isc[k+1] = -1;
      inters = vofi_get_side_intersections(impl_func,par,fse,x1,xfsl,baser,sdir,
                                           hs,nsub,xfsl.isc[k+1]);
      nsub += inters;
      for (i=1;i<=inters;i++)
	basei[nsub-i] = 1;
      nbt[k] = inters + 1;
      if (fse[0] < 0.) {
        sign_sect[k][0] = -1; sign_sect[k][1] = 1;
      }
      else {
        sign_sect[k][0] = 1; sign_sect[k][1] = -1;
      }
    }
    else {
      nointer = 0;
      consi = vofi_check_side_consistency(impl_func,par,x1,sdir,fse,hs); 
      if (consi == 0) {
	nointer = 1;
      }
      if (consi != 0) {
	f2pos = consi;
        sign_change = vofi_get_segment_min(impl_func,par,x1,sdir,fse,&xfsl,
                                           hs,f2pos);
        if (sign_change == 0) {
          nointer = 1;
        }
        if (sign_change != 0) {
          xfs_pt->isc[0]   = xfsl.isc[0]   = 1;
          xfs_pt->isc[k+1] = xfsl.isc[k+1] = 1;
          inters = vofi_get_side_intersections(impl_func,par,fse,x1,xfsl,baser,
                                 sdir,hs,nsub,xfsl.isc[k+1]);
          nsub += inters;
	  for (i=1;i<=inters;i++)
	    basei[nsub-i] = 1;
          nbt[k] = inters + 1;
          xfs_pt->sval = 0.5*(baser[nsub-1]+baser[nsub-2]);
          if (fsum > 0.) {
            sign_sect[k][0] = sign_sect[k][2] = 1; sign_sect[k][1] = -1;
          }
          else {
            sign_sect[k][0] = sign_sect[k][2] = -1; sign_sect[k][1] = 1;
          }
        }
      }

      /* - */
      if (nointer == 1) {
        atleast1 = 1;
        nbt[k] = 1;
        iside = k;
        isect = 1 - k;
        if (fsum > 0.)
          sign_sect[k][0] = 1;
        else if (fsum < 0.)
          sign_sect[k][0] = -1;
        if ((k == 0 && sign_sect[k][0] == 1) ||
            (k == 1 && sign_sect[k][0] == -1))
          down2up = -1;
      }
    }
  }

  ncheck = atleast1*MAX(nbt[0],nbt[1]);
  baser[nsub] = hs;
  basei[nsub] = 1;
  
  /* - */
  vofi_reorder(baser,basei,nsub);

  /* - */
  nsub = vofi_rm_segs(baser,basei,nsub);

  /* - */
  if (ncheck == nsub) 
    vofi_sector_new(sign_sect,nsect,ndire,nsub,iside,isect,down2up);
  /* - */
  else 
    vofi_sector_old(impl_func,par,x0,h0,baser,pdir,sdir,nsect,ndire,nsub);
  
  if (nsub == 1 && ndire[0] == 0)
    nsub = 0;
  
  return nsub;
}

/* -------------------------------------------------------------------------- */
vofi_int vofi_get_limits_inner_2D(integrand impl_func,vofi_void_cptr par,
                         vofi_creal x0[],vofi_creal h0[],min_data *xfs_pt,
                         vofi_real baser[],vofi_creal pdir[],vofi_creal sdir[],
                         vofi_int nsect[],vofi_int ndire[], vofi_cint nsub_int)
{
  vofi_int basei[NSEG],i,k,f2pos,sign_change,nsub,inters;
  vofi_real x1[NDIM],x2[NDIM],fse[NSE],hs,fs,fsum;
  min_data xfsl;

  xfsl = *xfs_pt;
  baser[0] = hs = 0.;
  basei[0] = nsub = 1;
  for (i=0;i<NDIM;i++) 
    hs += sdir[i]*h0[i];

  for (k=0;k<NSE;k++) {           /* - */
    sign_change = 1;
    for (i=0;i<NDIM;i++) {
      x1[i] = x0[i] + k*pdir[i]*h0[i];
      x2[i] = x1[i] + sdir[i]*hs;
    }
    fse[0] = impl_func(x1,par);
    fse[1] = impl_func(x2,par);
    if (fse[0]*fse[1] < 0.) {
      inters = vofi_get_side_intersections(impl_func,par,fse,x1,xfsl,baser,sdir,
                                           hs,nsub,-1);
      nsub += inters;
      for (i=1;i<=inters;i++)
	basei[nsub-i] = 1;
    }
    else if (xfsl.isc[k+1] == 1) {
      for (i=0;i<NDIM;i++) 
        xfsl.xval[i] = x1[i] + xfsl.sval*sdir[i];
      fs = impl_func(xfsl.xval,par);
      xfsl.fval = fs;
      fsum = fse[0] + fse[1];
      if (fs*fsum >= 0.) {
        if (fsum > 0.)
          f2pos = 1;
        else if (fsum < 0.)
          f2pos = -1;
        else
          f2pos = 0;
        if (f2pos != 0)
          sign_change = vofi_get_segment_min(impl_func,par,x1,sdir,fse,&xfsl,
                                             hs,f2pos);
      }
      if (sign_change != 0) {
        inters = vofi_get_side_intersections(impl_func,par,fse,x1,xfsl,baser,
                                             sdir,hs,nsub,1);
        nsub += inters;
        xfs_pt->sval = 0.5*(baser[nsub-1]+baser[nsub-2]);
	for (i=1;i<=inters;i++)
	  basei[nsub-i] = 1;
      }
    }
    else
      ;
  }
  
  baser[nsub] = hs;
  basei[nsub] = 1;
  /* - */
  vofi_reorder(baser,basei,nsub);

  /* - */
  nsub = vofi_rm_segs(baser,basei,nsub);
  
  /* - */
  if (nsub == nsub_int)
    ;
  /* - */
  else 
    vofi_sector_old(impl_func,par,x0,h0,baser,pdir,sdir,nsect,ndire,nsub);

  return nsub;    
}

/* -------------------------------------------------------------------------- */
vofi_int vofi_get_limits_edge_2D(integrand impl_func,vofi_void_cptr par,
                         vofi_creal x0[],vofi_creal h0[],min_data *xfs_pt,
                         vofi_real baser[],vofi_creal pdir[],vofi_creal sdir[],
                         vofi_cint nsub_int)
{
  vofi_int basei[NSEG],i,k,f2pos,sign_change,nsub,inters;
  vofi_real x1[NDIM],x2[NDIM],fse[NSE],hs,fs,fsum;
  min_data xfsl;

  xfsl = *xfs_pt;
  baser[0] = hs = 0.;
  basei[0] = nsub = 1;
  for (i=0;i<NDIM;i++) 
    hs += sdir[i]*h0[i];

  for (k=0;k<NSE;k++) {           /* - */
    sign_change = 1;
    for (i=0;i<NDIM;i++) {
      x1[i] = x0[i] + k*pdir[i]*h0[i];
      x2[i] = x1[i] + hs*sdir[i];
    }
    fse[0] = impl_func(x1,par);
    fse[1] = impl_func(x2,par);
    if (xfsl.isc[k+1] == -1) {
      if (fse[0]*fse[1] < 0.) {
        inters = vofi_get_side_intersections(impl_func,par,fse,x1,xfsl,baser,
                                             sdir,hs,nsub,-1);
      }
      else
        vofi_check_edge_consistency(impl_func,par,fse,x1,baser,sdir,hs,nsub);
      basei[nsub] = 1;
      nsub++;
    }
    else if (xfsl.isc[k+1] == 1) {
      for (i=0;i<NDIM;i++) 
        xfsl.xval[i] = x1[i] + xfsl.sval*sdir[i];
      fs = impl_func(xfsl.xval,par);
      xfsl.fval = fs;
      fsum = fse[0] + fse[1];
      if (fs*fsum >= 0.) {
        if (fsum > 0.)
          f2pos = 1;
        else if (fsum < 0.)
          f2pos = -1;
        else
          f2pos = 0;
        if (f2pos != 0)
          sign_change = vofi_get_segment_min(impl_func,par,x1,sdir,fse,&xfsl,
                                             hs,f2pos);
      }
      if (sign_change == 0) {
        baser[nsub] = xfsl.sval;
        baser[nsub+1] = xfsl.sval;
        basei[nsub] = 1;
        basei[nsub+1] = 1;
        nsub += 2;
      }
      else {
        inters = vofi_get_side_intersections(impl_func,par,fse,x1,xfsl,baser,
                                             sdir,hs,nsub,1);
        nsub += inters;
	for (i=1;i<=inters;i++)
	  basei[nsub-i] = 1;
      }
    }
    else
      ;
  }
  baser[nsub] = hs;
  basei[nsub] = 1;
  
  /* - */
  vofi_reorder(baser,basei,nsub);
  
  if (nsub == nsub_int)
    ;
  else {
    printf(" EXIT: in vofi_get_limits_edge_2D: never occured! \n"); 
    exit(1);
  }
  
  return nsub;    
}

/* -------------------------------------------------------------------------- */
void vofi_reorder(vofi_real baser[],vofi_int basei[],vofi_int nsub)
{
  vofi_int i,j,ks;
  vofi_real ls;
  
  for (j=2;j<nsub;j++) {            
    ls = baser[j];
    ks = basei[j];
    i = j-1;
    while (i > 0 && baser[i] > ls) {
      baser[i+1] = baser[i];
      basei[i+1] = basei[i];
      i--;
    }
    baser[i+1] = ls;
    basei[i+1] = ks;
  }

  return;
}

/* -------------------------------------------------------------------------- */
vofi_int vofi_rm_segs(vofi_real baser[],vofi_int basei[],vofi_int nsub)
{
  vofi_int i,j,ks;
  vofi_real ds,bs,be,eps;
  vofi_real eps2[2]={2.*EPS_LOC,EPS_SEGM};

  i = 0;
  bs = baser[0];
  be = baser[nsub];
  ds = baser[1]-baser[0];
  ks = basei[1]*basei[0];
  eps = eps2[ks];
  while (i<nsub) {                           
    if (ds < eps) {
      if (basei[i] == 1)
	;
      else if (basei[i+1] == 1)
	baser[i] = baser[i+1];
      for (j=i+1;j<nsub;j++) {
        baser[j] = baser[j+1];
	basei[j] = basei[j+1];
      }
      nsub--;
    }
    ds = baser[i+1]-baser[i];
    ks = basei[i+1]*basei[i];
    eps = eps2[ks];
    
    if (ds >= eps)     
    i++;
  }
  baser[0] = bs;
  baser[nsub] = be;

  return nsub;
}
/* -------------------------------------------------------------------------- */
void vofi_sector_new(vofi_int sign_sect[][NDIM],vofi_int nsect[],
		     vofi_int ndire[],vofi_cint nsub,vofi_cint iside,
		     vofi_cint isect,vofi_cint down2up)
{
  vofi_int i;
  
  for (i=1; i<=nsub; i++) {
    if (sign_sect[iside][0] * sign_sect[isect][i-1] > 0) {
      if (sign_sect[iside][0] < 0)
        nsect[i-1] = 1;
      else
        nsect[i-1] = 0;
      ndire[i-1] = 0;
    }
    else {
      nsect[i-1] = -1;
      ndire[i-1] = down2up;
    }
  }
  
  return;
}

/* -------------------------------------------------------------------------- */
void vofi_sector_old(integrand impl_func,vofi_void_cptr par,vofi_creal x0[],
                     vofi_creal h0[],vofi_creal base[],vofi_creal pdir[],
                     vofi_creal sdir[],vofi_int nsect[],vofi_int ndire[],
                     vofi_cint nsub)
{
  vofi_int i,j;
  vofi_real x1[NDIM],x2[NDIM],fse[NSE],mdpt;
  
  for (i=1; i<=nsub; i++) {
    mdpt = 0.5*(base[i] + base[i-1]);
    for (j=0;j<NDIM;j++) {
      x1[j] = x0[j] + sdir[j]*mdpt;
      x2[j] = x1[j] + pdir[j]*h0[j];
    }    
    fse[0] = impl_func(x1,par);
    fse[1] = impl_func(x2,par);
    if (fse[0]*fse[1] <= 0.) {
      nsect[i-1] = -1;
      if (fse[0] < 0.0 || fse[1] > 0.0)
        ndire[i-1] = 1;
      else 
        ndire[i-1] = -1;
    }
    else {
      if (fse[0] < 0.0) 
        nsect[i-1] = 1;
      else
        nsect[i-1] = 0;
      ndire[i-1] = 0;
    }
  }

  return;
}
