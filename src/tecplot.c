/**
 * @file      tecplot.c
 * @authors   Andrea Chierici, Leonardo Chirco, Vincent Le Chenadec, 
 *            Ruben Scardovelli, Philip Yecko and Stéphane Zaleski 
 * @date      April 15, 2021
 * @brief     Routines to print the heights in 2D in Tecplot data file 
 *            heights.dat, points in Tecplot data file arcline.dat, 
 *            triangles in Tecplot data file triangles.dat 
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
 * print data to files *.dat in Tecplot format                                *
 * FUNCTIONS:                                                                 *
 * void tecplot_heights : append endpoints of heights segment to file         *
 *                        "heights.dat"                                       *
 * void tecplot_arcline : append arc points to file "arcline.dat"             *
 * void tecplot_triangle: append triangle data to file "triangles.dat"        *
 * -------------------------------------------------------------------------- */
void tecplot_heights(vofi_creal x0[],vofi_creal h0[],vofi_creal pdir[],
                     vofi_creal sdir[],len_data xhp[])
{
  FILE *fp;
  vofi_int i,j,npt,f_sign,nseg,it0;
  vofi_real xb[NDIM],xt[NDIM],hp,xx,yy,zz;

  nseg = 0;
  if (xhp[1].np0 > 0)
    nseg = 2;
  else if (xhp[0].np0 > 0)
    nseg = 1;
  hp = 0.;
  for (i=0;i<NDIM;i++) 
    hp += pdir[i]*h0[i];
  fp = fopen("heights.dat","a");
  for (it0=0; it0<nseg; it0++) {
    npt = xhp[it0].np0;
    f_sign = xhp[it0].f_sign;
    fprintf(fp," ZONE I=%d, J=2, F=POINT \n",npt);
    for(j=1;j<=npt;j++) {
      for (i=0;i<NDIM;i++) {
	xb[i] = x0[i] + sdir[i]*xhp[it0].xt0[j];
	xt[i] = xb[i] + pdir[i]*hp;
      }
      if (f_sign > 0) {
        xx = xb[0] + pdir[0]*xhp[it0].ht0[j];
        yy = xb[1] + pdir[1]*xhp[it0].ht0[j];
        zz = xb[2] + pdir[2]*xhp[it0].ht0[j];
      }
      else {
        xx = xt[0] - pdir[0]*xhp[it0].ht0[j];
        yy = xt[1] - pdir[1]*xhp[it0].ht0[j];
        zz = xt[2] - pdir[2]*xhp[it0].ht0[j];
      }
      fprintf(fp," %13.6E    %13.6E    %13.6E \n",xx,yy,zz);
    }	
    for(j=1;j<=npt;j++) {
      for (i=0;i<NDIM;i++) {
        xb[i] = x0[i] + sdir[i]*xhp[it0].xt0[j];
        xt[i] = xb[i] + pdir[i]*hp;
      }
      if (f_sign > 0) {
        xx = xb[0];
        yy = xb[1];
        zz = xb[2];
      }
      else {
        xx = xt[0];
        yy = xt[1];
        zz = xt[2];
      }
      fprintf(fp," %13.6E    %13.6E    %13.6E \n",xx,yy,zz);
    }
  }
  fclose(fp);
  
  return;
}

/* -------------------------------------------------------------------------- */
void tecplot_arcline(vofi_creal x0[],vofi_creal pdir[],vofi_creal sdir[],
                     vofi_creal xs,vofi_creal xp,vofi_creal hp,
                     vofi_cint f_sign,FILE *fp)
{
  vofi_int i;
  vofi_real xx[NDIM],hh;

  hh = xp;
  if (f_sign < 0)
    hh = hp - xp;
  for (i=0;i<NDIM;i++) 
    xx[i] = x0[i] + sdir[i]*xs + pdir[i]*hh;

  fprintf(fp," %13.6E    %13.6E    %13.6E \n",xx[0],xx[1],0.);

  return;
}

/* -------------------------------------------------------------------------- */
void tecplot_triangle(vofi_creal x0[],vofi_creal pdir[],vofi_creal sdir[],
                      vofi_creal tdir[],vofi_creal xa[],vofi_creal xb[],
                      vofi_creal xc[],vofi_creal hp, vofi_cint f_sign)
{
  FILE *fp;
  vofi_int i;
  vofi_real x1[NDIM],ht;

  fp = fopen("triangles.dat","a");
  fprintf(fp," ZONE N = 3, E = 1, F=FEPOINT, ET=TRIANGLE \n");
  ht = xa[2];
  if (f_sign < 0)
    ht = hp - xa[2];
  for (i=0;i<NDIM;i++) 
    x1[i] = x0[i] + tdir[i]*xa[0] + sdir[i]*xa[1] + pdir[i]*ht;
  fprintf(fp," %13.6E    %13.6E    %13.6E \n",x1[0],x1[1],x1[2]);

  ht = xb[2];
  if (f_sign < 0)
    ht = hp - xb[2];
  for (i=0;i<NDIM;i++) 
    x1[i] = x0[i] + tdir[i]*xb[0] + sdir[i]*xb[1] + pdir[i]*ht;
  fprintf(fp," %13.6E    %13.6E    %13.6E \n",x1[0],x1[1],x1[2]);

  ht = xc[2];
  if (f_sign < 0)
    ht = hp - xc[2];
  for (i=0;i<NDIM;i++) 
    x1[i] = x0[i] + tdir[i]*xc[0] + sdir[i]*xc[1] + pdir[i]*ht;
  fprintf(fp," %13.6E    %13.6E    %13.6E \n",x1[0],x1[1],x1[2]);
  fprintf(fp,"         1         2         3 \n");
  fclose(fp);

  return;
}
