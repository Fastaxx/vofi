/**
 * @file      orderdirs.c
 * @authors   Andrea Chierici, Leonardo Chirco, Vincent Le Chenadec, 
 *            Ruben Scardovelli, Philip Yecko and Stéphane Zaleski 
 * @date      April 15, 2021
 * @brief     Routines to compute the cell type; if the cell is cut then 
 *            order coordinate directions and compute number of integration
 *            points along the secondary direction
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
 * a) calculate f on the cell vertices and the threshold value fth,           *
 * b) with fth, consistency and minimum checks determine if full/empty/cut    *
 * c) if cut, order coordinate directions and compute number of integration   *
 *    points along the secondary direction                                    *
 * INPUT:  pointer to the implicit function, arrays of function parameters    * 
 *         par, coordinates of minimum vertex x0, cell edges h0               *
 * OUTPUT: arrays of ordered coordinate directions pdir, sdir, (tdir),        *
 *         function value on vertices f0, min_data (see vofi_stddecl.h)       *
 * FUNCTIONS:                                                                 *
 * vofi_int vofi_order_dirs_2D: 2D problem                                    *
 * vofi_int vofi_order_dirs_3D: 3D problem                                    *
 * void vofi_xyz2pst          : from [x][y][z] to [p][s][t] coordinates       *
 * -------------------------------------------------------------------------- */
vofi_int vofi_order_dirs_2D(integrand impl_func,vofi_void_cptr par,
                         vofi_creal x0[],vofi_creal h0[],vofi_real pdir[],
                         vofi_real sdir[],vofi_real f0[][NSE],min_data *xfs_pt)
{
  vofi_int n0[NSE][NSE],nc[NDIM][NDIM];
  vofi_int i,j,np0,nm0,icc,check_dir,nmax0,jp,js,npt,iwgt,nix,niy;
  vofi_real fc[NDIM][NDIM],x1[NDIM],hh[NDIM],fgrad[NSE];
  vofi_real f0mod,fgradmod,fgradsq;
  vofi_real tmp,hm,fth,have,Kappa;
  vofi_real fx,fy,fxx,fyy,fxy;
  vofi_creal MIN_GRAD=1.0e-04;
  vofi_creal rwgt[5]={100.,50.,10.,2.,1.};
  vofi_creal a0=2.30477, a1=28.5312, a2=-46.2729, a3=56.9179;
  
  /* - */
  np0 = nm0 = 0;
  icc = check_dir = -1;
  nmax0 = 4;  
  for (i=0;i<NDIM;i++) 
    hh[i] = 0.5*h0[i];

  /* - */
  x1[2] = x0[2];
  for (i=0;i<NSE;i++)
    for (j=0;j<NSE;j++) {
      x1[0] = x0[0] + i*h0[0];
      x1[1] = x0[1] + j*h0[1];
      f0[i][j] = impl_func(x1,par);
      if (f0[i][j] > 0.) 
        np0++;
      else if (f0[i][j] < 0.) 
        nm0++;
  }
  
  /* - */
  fgrad[0] = 0.5*((f0[1][1]+f0[1][0]) - (f0[0][1]+f0[0][0]))/h0[0];
  fgrad[1] = 0.5*((f0[1][1]+f0[0][1]) - (f0[1][0]+f0[0][0]))/h0[1];
  fgradsq = Sq2(fgrad);
  fgradmod = MAX(sqrt(fgradsq),MIN_GRAD);
  hm = MAX(hh[0],hh[1]);
  fth = fgradmod*hm;
  
  if (np0*nm0 == 0) {    
    
    /* - */
    /* - */
    np0 = nm0 = 0;
    for (i=0;i<NSE;i++)
      for (j=0;j<NSE;j++) {
        f0mod = fabs(f0[i][j]); 
        if (f0mod > fth) {
          n0[i][j] = 0;
          if (f0[i][j] < 0.)
            nm0++;
          else
            np0++;
        }
        else
          n0[i][j] = 1;
    }  

    /* - */
    if (nm0 == nmax0) {                         
      icc = 1;
      return icc;
    }
    else if (np0 == nmax0) {                                         
      icc = 0;
      return icc;
    }

    /* - */
    check_dir = vofi_check_boundary_line(impl_func,par,x0,h0,f0,xfs_pt,n0);
    
    /* - */
    if (check_dir < 0) {
      if (nm0 > 0) {
        icc = 1;
      }
      else {
        icc = 0;
      }
      return icc;
    }
  }
  
  /* - */
  have = 0.5*(h0[0] + h0[1]);
  for (i=0;i<NSE;i++)
    for (j=0;j<NSE;j++) 
      fc[2*i][2*j] = f0[i][j];
  
  for (i=0;i<NDIM;i++)
    x1[i] = x0[i] + hh[i];

  fc[1][1] = impl_func(x1,par);
  
  for (i=0;i<=2;i+=2) {
    x1[0] = x0[0] + i*hh[0];
    fc[i][1] = impl_func(x1,par);
  }
  x1[0] -= hh[0];
  for (j=0;j<=2;j+=2) {
    x1[1] = x0[1] + j*hh[1];
    fc[1][j] = impl_func(x1,par);
  }

  /* - */
  for (i=0;i<NDIM;i++)
    for (j=0;j<NDIM;j++) {
      if (fc[i][j] > 0.)       nc[i][j] =  1;
      else if (fc[i][j] < 0.)  nc[i][j] = -1;
      else     	               nc[i][j] =  0;
  }

  /* - */
  /* - */
  nix = niy = 0;
  for (i=0; i<=NSE;i=i+2)
    if (nc[i][2]*nc[i][0] < 0) niy++;
  for (j=0; j<=NSE;j=j+2)
    if (nc[2][j]*nc[0][j] < 0) nix++;

  /* compute p-s directions based on intersections or a more */
  /* accurate gradient near the interface */
  if (niy > nix) {
    jp =1; js = 0;
  }
  else if (nix > niy) {
    jp = 0; js = 1;
  }
  else {
    fgrad[0] = fgrad[1] = 0.;
    for (i=0;i<NSE;i++)
      for (j=0;j<NSE;j++) {
	fx = 0.5*((fc[i+1][j+1]+fc[i+1][j])-(fc[i][j+1]+fc[i][j]))/hh[0];
	fy = 0.5*((fc[i+1][j+1]+fc[i][j+1])-(fc[i+1][j]+fc[i][j]))/hh[1];
	tmp = sqrt(fx*fx+fy*fy);
	fx = fx/tmp;
	fy = fy/tmp;
	iwgt = abs(nc[i+1][j+1]+nc[i+1][j]+nc[i][j+1]+nc[i][j]);
        tmp = rwgt[iwgt];
	fgrad[0] += fx*tmp;
	fgrad[1] += fy*tmp;
      }
    for (i=0;i<NSE;i++)
      fgrad[i] = fabs(fgrad[i]);  
    if (fgrad[0] >= fgrad[1]) {
      jp = 0; js = 1; 
    }
    else {
      jp = 1; js = 0;
    }
  }
  if (check_dir >= 0 && check_dir != jp) {
    js = jp; jp = check_dir;
  }
  pdir[jp] = 1.; sdir[js] = 1.;
  
  /* - */
  if (jp == 1) {
    tmp = f0[0][1]; f0[0][1] = f0[1][0]; f0[1][0] = tmp;
  }

  /* - */
  /* - */
  if (check_dir < 0) 
    vofi_check_secondary_side(impl_func,par,x0,h0,pdir,sdir,f0,xfs_pt,fth);

  /* - */
  /* - */
  fx  = (fc[2][1] - fc[0][1])*have/(2.*hh[0]);
  fy  = (fc[1][2] - fc[1][0])*have/(2.*hh[1]);
  fxx = (fc[2][1] + fc[0][1] - 2.*fc[1][1])*have*have/
        (hh[0]*hh[0]);
  fyy = (fc[1][2] + fc[1][0] - 2.*fc[1][1])*have*have/
        (hh[1]*hh[1]);
  fxy = (fc[2][2] - fc[2][0] - fc[0][2] + fc[0][0])*
        have*have/(4.*hh[0]*hh[1]);
  tmp = fx*fx + fy*fy;
  tmp = sqrt(tmp*tmp*tmp);
  Kappa = fabs(fxx*fy*fy - 2.*fx*fy*fxy + fx*fx*fyy)/tmp;
  tmp = a0 + Kappa*(a1 + Kappa*(a2 + a3*Kappa));
  npt = (vofi_int) ceil(tmp);
  xfs_pt->ipt = MAX(4,MIN(npt,NGLM));
  
  return icc;
}

/* -------------------------------------------------------------------------- */
vofi_int vofi_order_dirs_3D(integrand impl_func,vofi_void_cptr par,
              vofi_creal x0[],vofi_creal h0[],vofi_real pdir[],
              vofi_real sdir[],vofi_real tdir[],vofi_real f0[][NSE][NSE],
              min_data xfsp[])
{
  vofi_int n0[NSE][NSE][NSE],i,j,k,i0,j0,k0,ii,jj,kk;
  vofi_int np0,nm0,icc,check_dir,nmax0,jp,js,jt,npt;
  vofi_real fc[NDIM][NDIM][NDIM],fd[NDIM][NDIM],x1[NDIM],hh[NDIM],fgrad[NDIM];
  vofi_real f0mod,fgradmod,fgradsq;
  vofi_real sumf[NDIM],curv[NDIM],tmp,hm,fth,have,Kappa;
  vofi_real fx,fy,fz,fxx,fyy,fxy;
  vofi_creal MIN_GRAD=1.0e-04;
  vofi_creal a0=2.34607, a1=16.5515, a2=-5.53054, a3=54.0866;
  
  /* - */
  np0 = nm0 = 0;
  icc = check_dir = -1;
  nmax0 = 8;
  for (i=0;i<NDIM;i++) 
    hh[i] = 0.5*h0[i];

  /* - */
  for (i=0;i<NSE;i++)
    for (j=0;j<NSE;j++) 
      for (k=0;k<NSE;k++) {
        x1[0] = x0[0] + i*h0[0];
        x1[1] = x0[1] + j*h0[1];
        x1[2] = x0[2] + k*h0[2];
        f0[i][j][k] = impl_func(x1,par);
        if (f0[i][j][k] > 0.)  
          np0++;
        else if (f0[i][j][k] < 0.)  
          nm0++;
  }
  
  /* - */
  fgrad[0] = 0.25*( (f0[1][1][1]+f0[1][0][1]+f0[1][1][0]+f0[1][0][0]) -
             (f0[0][1][1]+f0[0][0][1]+f0[0][1][0]+f0[0][0][0]) )/h0[0];
  fgrad[1] = 0.25*( (f0[1][1][1]+f0[0][1][1]+f0[1][1][0]+f0[0][1][0]) -
             (f0[1][0][1]+f0[0][0][1]+f0[1][0][0]+f0[0][0][0]) )/h0[1];
  fgrad[2] = 0.25*( (f0[1][1][1]+f0[1][0][1]+f0[0][1][1]+f0[0][0][1]) -
             (f0[1][1][0]+f0[1][0][0]+f0[0][1][0]+f0[0][0][0]) )/h0[2]; 
  fgradsq = Sq3(fgrad);
  fgradmod = MAX(sqrt(fgradsq),MIN_GRAD);
  hm = MAX(hh[0],hh[1]);
  hm = MAX(hh[2],hm);
  fth = sqrt(2.)*fgradmod*hm;
  
  if (np0*nm0 == 0) {    
    
    /* - */
    /* - */
    np0 = nm0 = 0;
    for (i=0;i<NSE;i++)
      for (j=0;j<NSE;j++) 
        for (k=0;k<NSE;k++) {
          f0mod = fabs(f0[i][j][k]); 
          if (f0mod > fth) {
            n0[i][j][k] = 0;
            if (f0[i][j][k] < 0.)
              nm0++;
            else
              np0++;
          }
          else
            n0[i][j][k] = 1;
    }  

    /* - */
    if (nm0 == nmax0) {                         
      icc = 1;
      return icc;
    }
    else if (np0 == nmax0) {                                         
      icc = 0;
      return icc;
    }

    /* - */
    check_dir = vofi_check_boundary_surface(impl_func,par,x0,h0,f0,xfsp,n0);
    
    /* - */
    if (check_dir < 0) {
      if (nm0 > 0) {
        icc = 1;
      }
      else {
        icc = 0;
      }
      return icc;
    }
  }
  
  /* - */
  /* - */
  for (i=0;i<NSE;i++)
    for (j=0;j<NSE;j++) 
      for (k=0;k<NSE;k++) {
        fc[2*i][2*j][2*k] = f0[i][j][k];
  }
  x1[2] = x0[2] + hh[2];
  for (i=0;i<=NSE;i+=2)
    for (j=0;j<=NSE;j+=2) {
      x1[0] = x0[0] + i*hh[0];
      x1[1] = x0[1] + j*hh[1];
      fc[i][j][1] = impl_func(x1,par);
    }
  
  x1[0] = x0[0] + hh[0];
  x1[1] = x0[1] + hh[1];
  for (k=0;k<NDIM;k++) {
    x1[2] = x0[2] + k*hh[2];
    fc[1][1][k] = impl_func(x1,par);
  }
  
  x1[1] = x0[1] + hh[1];
  for (i=0;i<=NSE;i+=2) {
    x1[0] = x0[0] + i*hh[0];
    for (k=0;k<NDIM;k++) {
      x1[2] = x0[2] + k*hh[2];
      fc[i][1][k] = impl_func(x1,par);
    }
  }

  x1[0] = x0[0] + hh[0];
  for (j=0;j<=NSE;j+=2) {
    x1[1] = x0[1] + j*hh[1];
    for (k=0;k<NDIM;k++) {
      x1[2] = x0[2] + k*hh[2];
      fc[1][j][k] = impl_func(x1,par);
    }
  }
  /* - */
  /* - */
  fgrad[0] = fgrad[1] = fgrad[2] = 0.;
  for (i=0;i<NSE;i++)
    for (j=0;j<NSE;j++) 
      for (k=0;k<NSE;k++) {
        fx = 0.25*((fc[i+1][j+1][k+1]+fc[i+1][j][k+1]+fc[i+1][j+1][k]+
                    fc[i+1][j][k]) - (fc[i][j+1][k+1]+fc[i][j][k+1]+
                    fc[i][j+1][k]+fc[i][j][k]))/hh[0];
        fy = 0.25*((fc[i+1][j+1][k+1]+fc[i][j+1][k+1]+fc[i+1][j+1][k]+
                    fc[i][j+1][k]) - (fc[i+1][j][k+1]+fc[i][j][k+1]+  
                    fc[i+1][j][k]+fc[i][j][k]))/hh[1];
        fz = 0.25*((fc[i+1][j+1][k+1]+fc[i+1][j][k+1]+fc[i][j+1][k+1]+
                    fc[i][j][k+1]) - (fc[i+1][j+1][k]+fc[i+1][j][k]+
                    fc[i][j+1][k]+fc[i][j][k]))/hh[2];

        tmp = fc[i+1][j+1][k+1]+fc[i+1][j][k+1]+fc[i][j+1][k+1]+fc[i][j][k+1]+
              fc[i+1][j+1][k]+fc[i+1][j][k]+fc[i][j+1][k]+fc[i][j][k];
        tmp = MAX(fabs(tmp),MIN_GRAD);
        fgrad[0] += fx/tmp;
        fgrad[1] += fy/tmp;
        fgrad[2] += fz/tmp;
  }
  
  for (i=0;i<NDIM;i++)
    fgrad[i] = fabs(fgrad[i]);  
  /* - */
  jt = 2;
  if (fgrad[0] >= fgrad[1]) {
    jp = 0; js = 1; 
  }
  else {
    jp = 1; js = 0;
  }
  if (fgrad[2] > fgrad[jp]) {
    jt = js; js = jp; jp = 2;
  }
  else if (fgrad[2] > fgrad[js]) {
    jt = js; js = 2;
  }
  
  if (check_dir == 0 && xfsp[jp].isc[0] != 1) { 
    if (xfsp[js].isc[0] == 1) {
      tmp = js; js = jp; jp = tmp;
    }
    else {
      tmp = jt; jt = js; js = jp; jp = tmp;
    }
  }
  pdir[jp] = 1.; sdir[js] = 1.; tdir[jt] = 1.;
  
  /* - */
  vofi_xyz2pst(f0,jp,js,jt);

  /* - */
  /* - */
  if (check_dir >= 0)
    xfsp[4] = xfsp[jp];
  else
    vofi_check_secter_face(impl_func,par,x0,h0,pdir,sdir,tdir,f0,&xfsp[4],fth);

  /* - */
  vofi_check_tertiary_side(impl_func,par,x0,h0,pdir,sdir,tdir,f0,xfsp,fth);
  
  /* - */
  /* - */
  have = 0.5*(h0[jp] + h0[js]);
  for (k=0;k<NDIM;k++) {
    i0 = k*tdir[0]; j0 = k*tdir[1]; k0 = k*tdir[2];
    sumf[k] = 0.;
    for (i=0;i<NDIM;i++)
      for (j=0;j<NDIM;j++) {
	ii = i0 + i* (vofi_int) sdir[0] + j* (vofi_int) pdir[0];
	jj = j0 + i* (vofi_int) sdir[1] + j* (vofi_int) pdir[1];
	kk = k0 + i* (vofi_int) sdir[2] + j* (vofi_int) pdir[2];
	fd[i][j] = fc[ii][jj][kk];
	sumf[k] += fd[i][j];
      }
    sumf[k] = MAX(fabs(sumf[k]),EPS_NOT0);
    sumf[k] = 1./sumf[k];
    fx  = (fd[2][1] - fd[0][1])*have/(h0[js]);
    fy  = (fd[1][2] - fd[1][0])*have/(h0[jp]);
    fxx = (fd[2][1] + fd[0][1] - 2.*fd[1][1])*have*have/
          (hh[js]*hh[js]);
    fyy = (fd[1][2] + fd[1][0] - 2.*fd[1][1])*have*have/
          (hh[jp]*hh[jp]);
    fxy = (fd[2][2] - fd[2][0] - fd[0][2] + fd[0][0])*
          have*have/(h0[js]*h0[jp]);
    tmp = fx*fx + fy*fy;
    tmp = sqrt(tmp*tmp*tmp);
    curv[k] = fabs(fxx*fy*fy - 2.*fx*fy*fxy + fx*fx*fyy)/tmp;
  }
  Kappa = tmp = 0.;
  for (k=0;k<NDIM;k++) {
    Kappa += sumf[k]*curv[k];
    tmp += sumf[k];
  }
  Kappa = Kappa/tmp;
  tmp = a0 + Kappa*(a1 + Kappa*(a2 + a3*Kappa));
  npt = (vofi_int) ceil(tmp);
  xfsp[4].ipt = MAX(4,MIN(npt,NGLM));
  
  return icc;
}

/* -------------------------------------------------------------------------- */
void vofi_xyz2pst(vofi_real g0[][NSE][NSE],vofi_cint jp,vofi_cint js,
                  vofi_cint jt)
{
  vofi_real tmp;
  
  if (jt == 2 && js == 0) {      
      tmp = g0[0][1][0]; g0[0][1][0] = g0[1][0][0]; g0[1][0][0]=tmp;
      tmp = g0[0][1][1]; g0[0][1][1] = g0[1][0][1]; g0[1][0][1]=tmp;  
  }
  else if (jp == 2) {
    tmp = g0[0][0][1]; g0[0][0][1] = g0[1][0][0]; g0[1][0][0]=tmp;
    tmp = g0[0][1][1]; g0[0][1][1] = g0[1][1][0]; g0[1][1][0]=tmp;
    if (jt == 1) {
      tmp = g0[0][0][1]; g0[0][0][1] = g0[0][1][0]; g0[0][1][0]=tmp;
      tmp = g0[1][0][1]; g0[1][0][1] = g0[1][1][0]; g0[1][1][0]=tmp;
    }
  }
  else if (js == 2) {
    tmp = g0[0][0][1]; g0[0][0][1] = g0[0][1][0]; g0[0][1][0]=tmp;
    tmp = g0[1][0][1]; g0[1][0][1] = g0[1][1][0]; g0[1][1][0]=tmp;
    if (jp == 1) {
      tmp = g0[0][0][1]; g0[0][0][1] = g0[1][0][0]; g0[1][0][0]=tmp;
      tmp = g0[0][1][1]; g0[0][1][1] = g0[1][1][0]; g0[1][1][0]=tmp;
    }
  }
  
  return;
}
