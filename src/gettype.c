/**
 * @file      gettype.c
 * @authors   Andrea Chierici, Leonardo Chirco, Vincent Le Chenadec, 
 *            Ruben Scardovelli, Philip Yecko and Stéphane Zaleski 
 * @date      April 15, 2021
 * @brief     Routines to compute if the given cell is full, empty or cut
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
 * Routines to compute if the given cell is full/empty/cut                    *
 * INPUT:  pointer to the implicit function, arrays of function parameters    *
 *         par, coordinates of minimum vertex x0, cell edges h0, space        * 
 *         dimension ndim0                                                    * 
 * OUTPUT: icc (1/0/-1: full/empty/cut cell)                                  *
 * FUNCTIONS:                                                                 *
 * vofi_int vofi_get_cell_type: driver for the two cases                      *
 * vofi_int vofi_cell_type_2D : compute the cell type in 2D                   *
 * vofi_int vofi_cell_type_3D : compute the cell type in 3D                   *
 * -------------------------------------------------------------------------- */
vofi_int vofi_get_cell_type(integrand impl_func,vofi_void_cptr par,
			     vofi_creal xin[],vofi_creal h0[],vofi_cint ndim0) 
{
  vofi_int  icc;
  vofi_real x0[NDIM];
  
  /* - */
  if (ndim0 == 2) {
    x0[0] = xin[0]; x0[1] = xin[1]; x0[2] = 0.;
    icc = vofi_cell_type_2D(impl_func,par,x0,h0);
  }
  /* - */
  else if (ndim0 == 3) {
    x0[0] = xin[0]; x0[1] = xin[1]; x0[2] = xin[2];
    icc = vofi_cell_type_3D(impl_func,par,x0,h0);
  }
  /* - */
  else {                                                           
    printf(" EXIT: wrong value of variable ndim0! \n");
    exit(1);
  }

  return icc;
}

/* -------------------------------------------------------------------------- */
vofi_int vofi_cell_type_2D(integrand impl_func,vofi_void_cptr par,
		            vofi_creal x0[],vofi_creal h0[])
{
  vofi_int n0[NSE][NSE],i,j,np0,nm0,icc,check_dir,nmax0;
  vofi_real f0[NSE][NSE],x1[NDIM],fgrad[NSE];
  vofi_real f0mod,fgradmod,fgradsq,hm,fth;
  vofi_creal MIN_GRAD=1.0e-04;
  min_data  xfsp = {{0.,0.,0.},0.,0.,{0,0,0},0};

  /* - */
  np0 = nm0 = 0;
  icc = check_dir = -1;
  nmax0 = 4;  

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
  hm = 0.5*MAX(h0[0],h0[1]);
  fth = fgradmod*hm;
  
  if (np0*nm0 == 0) {    
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
    if (nm0 == nmax0)
      icc = 1;
    else if (np0 == nmax0) 
      icc = 0;

    /* - */
    if (icc < 0)
      check_dir = vofi_check_boundary_line(impl_func,par,x0,h0,f0,&xfsp,n0);
    
    /* - */
    if (check_dir < 0) {
      if (nm0 > 0) 
        icc = 1;
      else 
        icc = 0;
    }
  }
  
  return icc;
}

/* -------------------------------------------------------------------------- */
vofi_int vofi_cell_type_3D(integrand impl_func,vofi_void_cptr par,
                           vofi_creal x0[],vofi_creal h0[])
{
  vofi_int n0[NSE][NSE][NSE],i,j,k;
  vofi_int np0,nm0,icc,check_dir,nmax0;
  vofi_real f0[NSE][NSE][NSE],x1[NDIM],fgrad[NDIM];
  vofi_real f0mod,fgradmod,fgradsq,hm,fth;
  vofi_creal MIN_GRAD=1.0e-04;
  min_data  xfsp[3]={{{0.,0.,0.},0.,0.,{0,0,0},0},
	    {{0.,0.,0.},0.,0.,{0,0,0},0},{{0.,0.,0.},0.,0.,{0,0,0},0}};
  
  /* - */
  np0 = nm0 = 0;
  icc = check_dir = -1;
  nmax0 = 8;  

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
  hm = MAX(h0[0],h0[1]);
  hm = 0.5*MAX(h0[2],hm);
  fth = fgradmod*hm/sqrt(2.);
  
  if (np0*nm0 == 0) {    
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
    if (nm0 == nmax0)                          
      icc = 1;
    else if (np0 == nmax0) 
      icc = 0;

    /* - */
    if (icc < 0)
      check_dir = vofi_check_boundary_surface(impl_func,par,x0,h0,f0,xfsp,n0);
    
    /* - */
    if (check_dir < 0) {
      if (nm0 > 0) 
        icc = 1;
      else 
        icc = 0;
    }
  }
  
  return icc;
}
