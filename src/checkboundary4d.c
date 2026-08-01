/**
 * @file      checkboundary4d.c
 * @brief     4D analogues of checkboundary.c: minima inside the eight cubic
 *            3-faces of a hypercube, inside the two 3-faces normal to pdir,
 *            and on the eight cell edges parallel to udir. Plus the
 *            direction-generic form of the 3D face check, needed because a
 *            hyperplane of a 4D cell inherits its frame instead of choosing
 *            it.
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
 * One level up from checkboundary.c, with the same three jobs and the same   *
 * order:                                                                     *
 *   vofi_check_boundary_hyper : is there an interface hiding in the boundary *
 *                               of the cell? (2D looks on the sides, 3D on   *
 *                               the faces, 4D inside the cubic 3-faces)      *
 *   vofi_check_secter_cell    : minima in the two 3-faces normal to pdir     *
 *                               (the 3D vofi_check_secter_face)              *
 *   vofi_check_quaternary_side: sign changes on the eight edges parallel to  *
 *                               udir (the 3D vofi_check_tertiary_side)       *
 * and one extra:                                                             *
 *   vofi_check_boundary_cell  : vofi_check_boundary_surface with the three   *
 *                               directions passed in rather than hard-wired  *
 *                               to x,y,z, so that a hyperplane of a 4D cell  *
 *                               can be checked in the frame the 4D cell      *
 *                               ordered.                                     *
 * -------------------------------------------------------------------------- */

/* vertex of the 4D cell in (p,s,t,u) order */
#define VIDX4(m,n,q,r) ((m) + 2*(n) + 4*(q) + 8*(r))

vofi_int vofi_check_boundary_cell(integrand impl_func,vofi_void_cptr par,
                    vofi_creal x0[],vofi_creal h0[],vofi_real f0[][NSE][NSE],
                    min_data xfs[],vofi_int n0[][NSE][NSE],vofi_creal dir1[],
                    vofi_creal dir2[],vofi_creal dir3[])
{
  vofi_int i,j,k,m,check_dir;
  vofi_int nd[NDIM3][NSE]={{1,1},{1,1},{1,1}},sign_change;
  vofi_real fve[NVER],x1[NDIM],hh[NDIM3];
  vofi_creal *dirs[NDIM3];
  dir_data ipsc;
  min_data xfsl;

  dirs[0] = dir1; dirs[1] = dir2; dirs[2] = dir3;
  check_dir = -1;
  for (k=0;k<NDIM3;k++) {
    hh[k] = 0.;
    for (i=0;i<NDIM;i++)
      hh[k] += dirs[k][i]*h0[i];
  }

  for (i=0;i<NSE;i++)
    for (j=0;j<NSE;j++)
      for (k=0;k<NSE;k++) {
        if (n0[i][j][k] <= 0)
          continue;
        if (nd[0][i] > 0) {                    /* the face dir1 = i */
          nd[0][i] = 0;
          fve[0] = f0[i][0][0]; fve[1] = f0[i][1][0];
          fve[2] = f0[i][0][1]; fve[3] = f0[i][1][1];
          for (m=0;m<NDIM;m++)
            x1[m] = x0[m] + i*hh[0]*dir1[m];
          ipsc = vofi_check_face_consistency(impl_func,par,x1,h0,dir2,dir3,fve);
          if (ipsc.consi != 0) {
            memset(&xfsl,0,sizeof(xfsl));
            sign_change = vofi_get_face_min(impl_func,par,x1,h0,dir2,dir3,
                                            fve,&xfsl,ipsc);
            if (sign_change != 0) {
              xfs[0] = xfsl;
              xfs[0].isc[0] = 1; xfs[0].isc[i+1] = 1;
              check_dir = 0;
            }
          }
        }
        if (nd[1][j] > 0) {                    /* the face dir2 = j */
          nd[1][j] = 0;
          fve[0] = f0[0][j][0]; fve[1] = f0[1][j][0];
          fve[2] = f0[0][j][1]; fve[3] = f0[1][j][1];
          for (m=0;m<NDIM;m++)
            x1[m] = x0[m] + j*hh[1]*dir2[m];
          ipsc = vofi_check_face_consistency(impl_func,par,x1,h0,dir1,dir3,fve);
          if (ipsc.consi != 0) {
            memset(&xfsl,0,sizeof(xfsl));
            sign_change = vofi_get_face_min(impl_func,par,x1,h0,dir1,dir3,
                                            fve,&xfsl,ipsc);
            if (sign_change != 0) {
              xfs[1] = xfsl;
              xfs[1].isc[0] = 1; xfs[1].isc[j+1] = 1;
              check_dir = 1;
            }
          }
        }
        if (nd[2][k] > 0) {                    /* the face dir3 = k */
          nd[2][k] = 0;
          fve[0] = f0[0][0][k]; fve[1] = f0[1][0][k];
          fve[2] = f0[0][1][k]; fve[3] = f0[1][1][k];
          for (m=0;m<NDIM;m++)
            x1[m] = x0[m] + k*hh[2]*dir3[m];
          ipsc = vofi_check_face_consistency(impl_func,par,x1,h0,dir1,dir2,fve);
          if (ipsc.consi != 0) {
            memset(&xfsl,0,sizeof(xfsl));
            sign_change = vofi_get_face_min(impl_func,par,x1,h0,dir1,dir2,
                                            fve,&xfsl,ipsc);
            if (sign_change != 0) {
              xfs[2] = xfsl;
              xfs[2].isc[0] = 1; xfs[2].isc[k+1] = 1;
              check_dir = 2;
            }
          }
        }
        n0[i][j][k] = 0;
      }

  return check_dir;
}

/* -------------------------------------------------------------------------- *
 * f0 is in GLOBAL axis order here (bit a = index along axis a): this runs    *
 * before the directions are ordered. The per-axis results go to xfs[0..3],   *
 * which vofi_order_dirs_4D consults to make the axis carrying the minimum    *
 * the primary direction -- exactly as vofi_order_dirs_3D does with the       *
 * output of vofi_check_boundary_surface.                                     *
 * -------------------------------------------------------------------------- */
vofi_int vofi_check_boundary_hyper(integrand impl_func,vofi_void_cptr par,
                          vofi_creal x0[],vofi_creal h0[],vofi_creal f0[],
                          min_data xfs[],vofi_int n0[])
{
  vofi_int a,b,k,v,w,sd,check_dir,sign_change,done[NDIM][NSE];
  vofi_int jf[NDIM3];
  vofi_real fve[NVERC],x1[NDIM],edir[NDIM][NDIM];
  dir_data ipsc;
  min_data xfsl;

  check_dir = -1;
  memset(done,0,sizeof(done));
  memset(edir,0,sizeof(edir));
  for (a=0;a<NDIM;a++)
    edir[a][a] = 1.;

  for (v=0;v<NVERH;v++) {
    if (n0[v] <= 0)
      continue;
    for (a=0;a<NDIM;a++) {                  /* the four 3-faces through v */
      sd = (v >> a) & 1;
      if (done[a][sd] || h0[a] <= 0.)
        continue;
      done[a][sd] = 1;

      k = 0;                                /* free axes, ascending order */
      for (b=0;b<NDIM;b++)
        if (b != a)
          jf[k++] = b;
      for (w=0;w<NVERC;w++)
        fve[w] = f0[ (sd << a) | (((w   ) & 1) << jf[0])
                              | (((w >> 1) & 1) << jf[1])
                              | (((w >> 2) & 1) << jf[2]) ];
      for (b=0;b<NDIM;b++)
        x1[b] = x0[b];
      x1[a] = x0[a] + sd*h0[a];

      ipsc = vofi_check_cell_consistency(impl_func,par,x1,h0,edir[jf[0]],
                                         edir[jf[1]],edir[jf[2]],fve);
      if (ipsc.consi != 0) {
        memset(&xfsl,0,sizeof(xfsl));
        sign_change = vofi_get_cell_min(impl_func,par,x1,h0,edir[jf[0]],
                                        edir[jf[1]],edir[jf[2]],fve,&xfsl,ipsc);
        if (sign_change != 0) {
          xfs[a] = xfsl;
          xfs[a].isc[0] = 1; xfs[a].isc[sd+1] = 1;
          check_dir = a;
        }
      }
    }
    n0[v] = 0;
  }

  return check_dir;
}

/* -------------------------------------------------------------------------- *
 * minima inside the two 3-faces normal to pdir: the 4D vofi_check_secter_face *
 * -------------------------------------------------------------------------- */
void vofi_check_secter_cell(integrand impl_func,vofi_void_cptr par,
                vofi_creal x0[],vofi_creal h0[],vofi_creal pdir[],
                vofi_creal sdir[],vofi_creal tdir[],vofi_creal udir[],
                vofi_creal f0[],min_data *xfs_pt,vofi_creal fth)
{
  vofi_int i,m,n,q,r,w,np0,nm0,small,sign_change;
  vofi_real x1[NDIM],fve[NVERC];
  dir_data ipsc;
  min_data xfsl;

  for (i=0;i<=NDIM;i++)
    xfs_pt->isc[i] = 0;

  for (m=0;m<NSE;m++) {
    np0 = nm0 = small = 0;
    for (n=0;n<NSE;n++)
      for (q=0;q<NSE;q++)
        for (r=0;r<NSE;r++) {
          w = n + 2*q + 4*r;
          fve[w] = f0[VIDX4(m,n,q,r)];
          if (fve[w] > 0.)
            np0++;
          else if (fve[w] < 0.)
            nm0++;
          if (fabs(fve[w]) <= fth)
            small++;
        }

    if (nm0*np0 > 0)             /* the 3-face is cut through its vertices */
      ;
    else if (small == 0)         /* every vertex is too far from the zero  */
      ;
    else {
      for (i=0;i<NDIM;i++)
        x1[i] = x0[i] + m*pdir[i]*h0[i];
      ipsc = vofi_check_cell_consistency(impl_func,par,x1,h0,sdir,tdir,udir,
                                         fve);
      if (ipsc.consi != 0) {
        memset(&xfsl,0,sizeof(xfsl));
        sign_change = vofi_get_cell_min(impl_func,par,x1,h0,sdir,tdir,udir,
                                        fve,&xfsl,ipsc);
        if (sign_change != 0) {
          *xfs_pt = xfsl;
          xfs_pt->isc[0] = 1; xfs_pt->isc[m+1] = 1;
        }
      }
    }
  }

  return;
}

/* -------------------------------------------------------------------------- *
 * sign changes on the eight edges parallel to udir: the 4D                   *
 * vofi_check_tertiary_side. Edge l0 = 4m + 2n + q for (p,s,t) = (m,n,q).     *
 * -------------------------------------------------------------------------- */
void vofi_check_quaternary_side(integrand impl_func,vofi_void_cptr par,
                vofi_creal x0[],vofi_creal h0[],vofi_creal pdir[],
                vofi_creal sdir[],vofi_creal tdir[],vofi_creal udir[],
                vofi_creal f0[],min_data xfs[],vofi_creal fth)
{
  vofi_int i,l0,m,n,q,consi,sign_change,f2pos;
  vofi_real x1[NDIM],fse[NSE],hu;
  min_data xfsl;

  hu = 0.;
  for (i=0;i<NDIM;i++)
    hu += udir[i]*h0[i];
  for (l0=0;l0<NVERC;l0++)
    for (i=0;i<=NDIM;i++)
      xfs[l0].isc[i] = 0;

  for (m=0;m<NSE;m++)
    for (n=0;n<NSE;n++)
      for (q=0;q<NSE;q++) {
        l0 = 4*m + 2*n + q;
        fse[0] = f0[VIDX4(m,n,q,0)];
        fse[1] = f0[VIDX4(m,n,q,1)];
        if (fse[0]*fse[1] < 0.) {
          xfs[l0].isc[0] = 1; xfs[l0].isc[1] = -1;
        }
        else {
          if (fabs(fse[0]) > fth && fabs(fse[1]) > fth)
            ;
          else {
            for (i=0;i<NDIM;i++)
              x1[i] = x0[i] + m*pdir[i]*h0[i] + n*sdir[i]*h0[i]
                            + q*tdir[i]*h0[i];
            consi = vofi_check_side_consistency(impl_func,par,x1,udir,fse,hu);
            if (consi != 0) {
              f2pos = consi;
              memset(&xfsl,0,sizeof(xfsl));
              sign_change = vofi_get_segment_min(impl_func,par,x1,udir,fse,
                                                 &xfsl,hu,f2pos);
              if (sign_change != 0) {
                for (i=0;i<NDIM;i++)
                  xfs[l0].xval[i] = xfsl.xval[i];
                xfs[l0].sval = xfsl.sval; xfs[l0].fval = xfsl.fval;
                xfs[l0].isc[0] = 1; xfs[l0].isc[1] = 1;
              }
            }
          }
        }
      }

  return;
}
