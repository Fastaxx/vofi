/**
 * @file      checkboundary.c
 * @authors   Andrea Chierici, Leonardo Chirco, Vincent Le Chenadec, 
 *            Ruben Scardovelli, Philip Yecko and Stéphane Zaleski 
 * @date      April 15, 2021
 * @brief     Routines to check if there is a double intersection on a 
 *            side (2D or 3D) or a "cap" intersection on a cell face (3D)
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
 * Routines to check if there is a double intersection on a side or           *
 *  segment, or a "cap" intersection on a cell face                           *
 * INPUT:  pointer to the implicit function, pointer to function parameters   * 
 *         par, coordinates of minimum vertex x0, cell edges h0, function     *
 *         value on vertices f0, min_data (see vofi_stddecl.h), flag on       *
 *         vertices n0, order coordinate directions pdir, sdir, (tdir),       *
 *         and threshold value fth                                            *
 * OUTPUT: integer check_dir, with the result of the check, or nothing        *
 * FUNCTIONS:                                                                 *
 * vofi_int vofi_check_boundary_line: check on the cell sides (2D)            *
 * void vofi_check_secondary_side   : check on the two sides along secondary  *
 *                                    direction (2D)                          *
 * vofi_int vofi_check_boundary_surface: check on the cell faces (3D)         *
 * void vofi_check_secter_face  : check on the two faces normal to the        * 
 *                                principal direction (3D)                    *
 * void vofi_check_tertiary_side: check on the four sides along the tertiary  *
 *                                direction (3D)                              *
 * -------------------------------------------------------------------------- */
vofi_int vofi_check_boundary_line(integrand impl_func,vofi_void_cptr par,
                    vofi_creal x0[],vofi_creal h0[],vofi_real f0[][NSE],
                    min_data *xfs_pt,vofi_int n0[][NSE])
{
  vofi_int i,j,k;
  vofi_int nx[NSE]={1,1},ny[NSE]={1,1},consi,sign_change,f2pos,check_dir;
  vofi_creal sidedirx[NDIM]={1.,0.,0.},sidediry[NDIM]={0.,1.,0.};
  vofi_real fse[NSE],x1[NDIM];
  min_data xfsl={{0.,0.,0.},0.,0.,{0,0,0},0};

  check_dir = -1; 
  for (i=0;i<=1;i++)
    for (j=0;j<=1;j++) {
      if (n0[i][j] > 0) {
        if (ny[i] > 0) {                          /* - */
          ny[i] = 0;
          fse[0] = f0[i][0];
          fse[1] = f0[i][1];
          for (k=0;k<NDIM;k++)
            x1[k] = x0[k] + i*sidedirx[k]*h0[0];
          consi = vofi_check_side_consistency(impl_func,par,x1,sidediry,fse,
                                              h0[1]);
          if (consi != 0) {
            f2pos = consi;
            sign_change = vofi_get_segment_min(impl_func,par,x1,sidediry,fse,
                                               &xfsl,h0[1],f2pos);
            if (sign_change != 0) {
              *xfs_pt = xfsl;
              xfs_pt->isc[0] = 1; xfs_pt->isc[i+1] = 1;
              check_dir = 0;                     /* - */
            }
          }
        }
        if (nx[j] > 0) {                          /* - */
          nx[j] = 0;
          fse[0] = f0[0][j];
          fse[1] = f0[1][j];
          for (k=0;k<NDIM;k++)
            x1[k] = x0[k] + j*sidediry[k]*h0[1];
          consi = vofi_check_side_consistency(impl_func,par,x1,sidedirx,fse,
                                              h0[0]);
          if (consi != 0) {
            f2pos = consi;
            sign_change = vofi_get_segment_min(impl_func,par,x1,sidedirx,fse,
                                               &xfsl,h0[0],f2pos);
            if (sign_change != 0) {
              *xfs_pt = xfsl;
              xfs_pt->isc[0] = 1; xfs_pt->isc[j+1] = 1;
              check_dir = 1;                     /* - */
            }
          }
        }
        n0[i][j] = 0;
      }
    }
  
  return check_dir;
}

/* -------------------------------------------------------------------------- */
void vofi_check_secondary_side(integrand impl_func,vofi_void_cptr par,
          vofi_creal x0[],vofi_creal h0[],vofi_creal pdir[],vofi_creal sdir[],
          vofi_real f0[][NSE],min_data *xfs_pt,vofi_creal fth)
{
  vofi_int i,k,consi,f2pos,sign_change;
  vofi_real x1[NDIM],fse[NSE],hs;
  min_data xfsl={{0.,0.,0.},0.,0.,{0,0,0},0};

  hs = 0.;
  for (i=0;i<NDIM;i++) 
    hs += sdir[i]*h0[i];

  for (k=0;k<NSE;k++) { 
    fse[0] = f0[k][0];
    fse[1] = f0[k][1];
    if (fse[0]*fse[1] < 0.) {
      xfs_pt->isc[0] = 1; xfs_pt->isc[k+1] = -1;
    }
    else {
      if (fabs(fse[0]) > fth && fabs(fse[1]) > fth)
        ;
      else {
        for (i=0;i<NDIM;i++)
          x1[i] = x0[i] + k*pdir[i]*h0[i];
        consi = vofi_check_side_consistency(impl_func,par,x1,sdir,fse,hs);
        if (consi != 0) {
          f2pos = consi;
          sign_change = vofi_get_segment_min(impl_func,par,x1,sdir,fse,&xfsl,
                                       hs,f2pos);
          if (sign_change != 0) {
            *xfs_pt = xfsl;
            xfs_pt->isc[0] = 1; xfs_pt->isc[k+1] = 1;
          }
        }
      }   
    }
  }
  
  return;
}

/* -------------------------------------------------------------------------- */
vofi_int vofi_check_boundary_surface(integrand impl_func,vofi_void_cptr par,
                    vofi_creal x0[],vofi_creal h0[],vofi_real f0[][NSE][NSE],
                    min_data xfs[],vofi_int n0[][NSE][NSE])
{
  vofi_int i,j,k,m;
  vofi_int nx[NSE]={1,1},ny[NSE]={1,1},nz[NSE]={1,1},sign_change,check_dir;
  vofi_creal sidedirx[NDIM]={1.,0.,0.},sidediry[NDIM]={0.,1.,0.},
             sidedirz[NDIM]={0.,0.,1.};
  vofi_real fve[NVER],x1[NDIM];
  dir_data ipsc;
  min_data xfsl={{0.,0.,0.},0.,0.,{0,0,0},0};
  
  check_dir = -1;
  for (i=0;i<=1;i++)
    for (j=0;j<=1;j++)
      for (k=0;k<=1;k++) {
        if (n0[i][j][k] > 0) {
          if (nx[i] > 0) {                      /* - */
            nx[i] = 0;
            fve[0] = f0[i][0][0];
            fve[1] = f0[i][1][0];
            fve[2] = f0[i][0][1];
            fve[3] = f0[i][1][1];
            for (m=0;m<NDIM;m++)
              x1[m] = x0[m] + i*sidedirx[m]*h0[0];
            ipsc = vofi_check_face_consistency(impl_func,par,x1,h0,sidediry,
                                               sidedirz,fve);
            if (ipsc.consi != 0) {
              sign_change = vofi_get_face_min(impl_func,par,x1,h0,sidediry,
                                              sidedirz,fve,&xfsl,ipsc);
              if (sign_change != 0) {
                xfs[0] = xfsl;		
                xfs[0].isc[0] = 1; xfs[0].isc[i+1] = 1;
                check_dir = 0;     /* - */
              }
            }
          }
          if (ny[j] > 0) {                      /* - */
            ny[j] = 0;
            fve[0] = f0[0][j][0];
            fve[1] = f0[1][j][0];
            fve[2] = f0[0][j][1];
            fve[3] = f0[1][j][1];
            for (m=0;m<NDIM;m++)
              x1[m] = x0[m] + j*sidediry[m]*h0[1];
            ipsc = vofi_check_face_consistency(impl_func,par,x1,h0,sidedirx,
                                               sidedirz,fve);
            if (ipsc.consi != 0) {
              sign_change = vofi_get_face_min(impl_func,par,x1,h0,sidedirx,
                                              sidedirz,fve,&xfsl,ipsc);
              if (sign_change != 0) {
                xfs[1] = xfsl;		
                xfs[1].isc[0] = 1; xfs[1].isc[j+1] = 1;
                check_dir = 0;     /* - */
              }
            }
          }
          if (nz[k] > 0) {                      /* - */
            nz[k] = 0;
            fve[0] = f0[0][0][k];
            fve[1] = f0[1][0][k];
            fve[2] = f0[0][1][k];
            fve[3] = f0[1][1][k];
            for (m=0;m<NDIM;m++)
              x1[m] = x0[m] + k*sidedirz[m]*h0[2];
            ipsc = vofi_check_face_consistency(impl_func,par,x1,h0,sidedirx,
                                               sidediry,fve);
            if (ipsc.consi != 0) {
              sign_change = vofi_get_face_min(impl_func,par,x1,h0,sidedirx,
                                              sidediry,fve,&xfsl,ipsc);
              if (sign_change != 0) {
                xfs[2] = xfsl;		
                xfs[2].isc[0] = 1; xfs[2].isc[k+1] = 1;
                check_dir = 0;     /* - */
              }
            }
          }
          n0[i][j][k] = 0;
        }
      }

  return check_dir;
}

/* -------------------------------------------------------------------------- */
void vofi_check_secter_face(integrand impl_func,vofi_void_cptr par,
		vofi_creal x0[],vofi_creal h0[],vofi_creal pdir[],
		vofi_creal sdir[],vofi_creal tdir[],vofi_real f0[][NSE][NSE],
                min_data *xfs_pt,vofi_creal fth)
{
  vofi_int i,m,nm0,np0,sign_change;
  vofi_real x1[NDIM],fve[NVER];
  dir_data ipsc;
  min_data xfsl={{0.,0.,0.},0.,0.,{0,0,0},0};

  for (i=0;i<NDIM;i++)
    xfs_pt->isc[i] = 0;
  
  for (m=0;m<NSE;m++) { 
    np0 = nm0 = 0;
    fve[0] = f0[m][0][0];
    if (fve[0] > 0.) 
      np0++;
    else if (fve[0] < 0.) 
      nm0++;
    fve[1] = f0[m][1][0];
    if (fve[1] > 0.) 
      np0++;
    else if (fve[1] < 0.) 
      nm0++;
    fve[2] = f0[m][0][1];
    if (fve[2] > 0.) 
      np0++;
    else if (fve[2] < 0.) 
      nm0++;
    fve[3] = f0[m][1][1];
    if (fve[3] > 0.) 
      np0++;
    else if (fve[3] < 0.) 
      nm0++;

    if (nm0*np0 > 0) {
      ;
    }
    else {
      if (fabs(fve[0]) > fth && fabs(fve[1]) > fth && fabs(fve[2]) > fth
          && fabs(fve[3]) > fth)
        ;
      else {
        for (i=0;i<NDIM;i++)
          x1[i] = x0[i] + m*pdir[i]*h0[i];
        ipsc = vofi_check_face_consistency(impl_func,par,x1,h0,sdir,tdir,fve);
        if (ipsc.consi != 0) {
          sign_change = vofi_get_face_min(impl_func,par,x1,h0,sdir,tdir,
                                          fve,&xfsl,ipsc);
          if (sign_change != 0) {
            *xfs_pt = xfsl;
            xfs_pt->isc[0] = 1; xfs_pt->isc[m+1] = 1;
          }
        }
      }   
    }
  }

  return;
}

/* -------------------------------------------------------------------------- */
void vofi_check_tertiary_side(integrand impl_func,vofi_void_cptr par,
          vofi_creal x0[],vofi_creal h0[],vofi_creal pdir[],vofi_creal sdir[],
          vofi_creal tdir[],vofi_real f0[][NSE][NSE],min_data xfs[],
	  vofi_creal fth)
{
  vofi_int i,l0,m,n,consi,sign_change,f2pos;
  vofi_real x1[NDIM],fse[NSE],ht;
  min_data xfsl={{0.,0.,0.},0.,0.,{0,0,0},0};

  ht = 0.;
  for (i=0;i<NDIM;i++) 
    ht += tdir[i]*h0[i];
  
  for (i=0;i<=NSE;i++) {
     xfs[3].isc[i] = xfs[2].isc[i] = xfs[1].isc[i] = xfs[0].isc[i] = 0;
  }
  /* - */
  for (m=0;m<NSE;m++) { 
    for (n=0;n<NSE;n++) { 
      fse[0] = f0[m][n][0];
      fse[1] = f0[m][n][1];
      l0 = 2*m + n;
      if (fse[0]*fse[1] < 0.) {
        xfs[l0].isc[0] = 1; xfs[l0].isc[1] = -1;
      }
      else {
        if (fabs(fse[0]) > fth && fabs(fse[1]) > fth)
          ;
        else {
          for (i=0;i<NDIM;i++)
            x1[i] = x0[i] + m*pdir[i]*h0[i] + n*sdir[i]*h0[i];
          consi = vofi_check_side_consistency(impl_func,par,x1,tdir,fse,ht);
          if (consi != 0) {
            f2pos = consi;
            sign_change = vofi_get_segment_min(impl_func,par,x1,tdir,fse,&xfsl,
                                               ht,f2pos);
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
  }
  
  return;
}

/* -------------------------------------------------------------------------- */
vofi_int vofi_check_boundary_hypercube(integrand impl_func, vofi_void_cptr par,
                                      vofi_creal x0[], vofi_creal h0[],
                                      vofi_real f0[][NSE][NSE][NSE],
                                      min4d_data xfsp[], vofi_int n0[][NSE][NSE][NSE])
{
  vofi_int i, j, k, l, m, check_dir, sign_change;
  vofi_real x1[NDIM], sidedirw[NDIM], sidedirx[NDIM], sidediry[NDIM], sidedirz[NDIM];
  vofi_real fve[8];  /* Values at vertices of a 3D face in 4D */
  dir_data ipsc;
  min_data xfsl={{0.,0.,0.,0.},0.,0.,{0,0,0,0},0};
  
  /* Initialize direction vectors */
  for (i=0; i<NDIM; i++) {
    sidedirw[i] = sidedirx[i] = sidediry[i] = sidedirz[i] = 0.0;
  }
  
  /* Set up direction vectors for each dimension */
  sidedirw[0] = 1.0;  /* w-direction */
  sidedirx[1] = 1.0;  /* x-direction */
  sidediry[2] = 1.0;  /* y-direction */
  sidedirz[3] = 1.0;  /* z-direction */
  
  check_dir = -1;
  
  /* ---------- CHECK FACES PERPENDICULAR TO W DIRECTION (1ST DIMENSION) ---------- */
  for (m=0; m<NSE; m+=NSE-1) {  /* m=0 and m=1 (lower and upper face) */
    /* Extract function values at this 3D hyperface */
    for (i=0; i<NSE; i++)
      for (j=0; j<NSE; j++)
        for (k=0; k<NSE; k++) {
          fve[i*4 + j*2 + k] = f0[m][i][j][k];
        }
    
    /* Initialize test point at the hyperface */
    for (i=0; i<NDIM; i++)
      x1[i] = x0[i] + m*sidedirw[i]*h0[0];
    
    /* Check this 3D face for consistency/sign change */
    ipsc = vofi_check_face_consistency_3D(impl_func, par, x1, h0, sidedirx, sidediry, sidedirz, fve);
    if (ipsc.consi != 0) {
      sign_change = vofi_get_face_min_3D(impl_func, par, x1, h0, sidedirx, sidediry, sidedirz, fve, &xfsl, ipsc);
      if (sign_change != 0) {
        for (i=0; i<NDIM; i++)
          xfsp[0].xval[i] = xfsl.xval[i];
        xfsp[0].fval = xfsl.fval;
        xfsp[0].sval = xfsl.sval;
        xfsp[0].isc[0] = 1; 
        xfsp[0].isc[1] = m;
        check_dir = 0;
        return check_dir;  /* Return early if found */
      }
    }
  }
  
  /* ---------- CHECK FACES PERPENDICULAR TO X DIRECTION (2ND DIMENSION) ---------- */
  for (m=0; m<NSE; m+=NSE-1) {  /* m=0 and m=1 (lower and upper face) */
    /* Extract function values at this 3D hyperface */
    for (i=0; i<NSE; i++)
      for (j=0; j<NSE; j++)
        for (k=0; k<NSE; k++) {
          fve[i*4 + j*2 + k] = f0[i][m][j][k];
        }
    
    /* Initialize test point at the hyperface */
    for (i=0; i<NDIM; i++)
      x1[i] = x0[i] + m*sidedirx[i]*h0[1];
    
    /* Check this 3D face for consistency/sign change */
    ipsc = vofi_check_face_consistency_3D(impl_func, par, x1, h0, sidedirw, sidediry, sidedirz, fve);
    if (ipsc.consi != 0) {
      sign_change = vofi_get_face_min_3D(impl_func, par, x1, h0, sidedirw, sidediry, sidedirz, fve, &xfsl, ipsc);
      if (sign_change != 0) {
        for (i=0; i<NDIM; i++)
          xfsp[1].xval[i] = xfsl.xval[i];
        xfsp[1].fval = xfsl.fval;
        xfsp[1].sval = xfsl.sval;
        xfsp[1].isc[0] = 1; 
        xfsp[1].isc[2] = m;
        check_dir = 1;
        return check_dir;  /* Return early if found */
      }
    }
  }
  
  /* ---------- CHECK FACES PERPENDICULAR TO Y DIRECTION (3RD DIMENSION) ---------- */
  for (m=0; m<NSE; m+=NSE-1) {  /* m=0 and m=1 (lower and upper face) */
    /* Extract function values at this 3D hyperface */
    for (i=0; i<NSE; i++)
      for (j=0; j<NSE; j++)
        for (k=0; k<NSE; k++) {
          fve[i*4 + j*2 + k] = f0[i][j][m][k];
        }
    
    /* Initialize test point at the hyperface */
    for (i=0; i<NDIM; i++)
      x1[i] = x0[i] + m*sidediry[i]*h0[2];
    
    /* Check this 3D face for consistency/sign change */
    ipsc = vofi_check_face_consistency_3D(impl_func, par, x1, h0, sidedirw, sidedirx, sidedirz, fve);
    if (ipsc.consi != 0) {
      sign_change = vofi_get_face_min_3D(impl_func, par, x1, h0, sidedirw, sidedirx, sidedirz, fve, &xfsl, ipsc);
      if (sign_change != 0) {
        for (i=0; i<NDIM; i++)
          xfsp[2].xval[i] = xfsl.xval[i];
        xfsp[2].fval = xfsl.fval;
        xfsp[2].sval = xfsl.sval;
        xfsp[2].isc[0] = 1; 
        xfsp[2].isc[3] = m;
        check_dir = 2;
        return check_dir;  /* Return early if found */
      }
    }
  }
  
  /* ---------- CHECK FACES PERPENDICULAR TO Z DIRECTION (4TH DIMENSION) ---------- */
  for (m=0; m<NSE; m+=NSE-1) {  /* m=0 and m=1 (lower and upper face) */
    /* Extract function values at this 3D hyperface */
    for (i=0; i<NSE; i++)
      for (j=0; j<NSE; j++)
        for (k=0; k<NSE; k++) {
          fve[i*4 + j*2 + k] = f0[i][j][k][m];
        }
    
    /* Initialize test point at the hyperface */
    for (i=0; i<NDIM; i++)
      x1[i] = x0[i] + m*sidedirz[i]*h0[3];
    
    /* Check this 3D face for consistency/sign change */
    ipsc = vofi_check_face_consistency_3D(impl_func, par, x1, h0, sidedirw, sidedirx, sidediry, fve);
    if (ipsc.consi != 0) {
      sign_change = vofi_get_face_min_3D(impl_func, par, x1, h0, sidedirw, sidedirx, sidediry, fve, &xfsl, ipsc);
      if (sign_change != 0) {
        for (i=0; i<NDIM; i++)
          xfsp[3].xval[i] = xfsl.xval[i];
        xfsp[3].fval = xfsl.fval;
        xfsp[3].sval = xfsl.sval;
        xfsp[3].isc[0] = 1; 
        xfsp[3].isc[4] = m;
        check_dir = 3;
        return check_dir;  /* Return early if found */
      }
    }
  }
  
  return check_dir;  /* Return -1 if no intersection was found */
}

/* -------------------------------------------------------------------------- */
dir_data vofi_check_face_consistency_3D(integrand impl_func, vofi_void_cptr par,
                                       vofi_creal x0[], vofi_creal h0[],
                                       vofi_creal dir1[], vofi_creal dir2[], 
                                       vofi_creal dir3[], vofi_creal fve[])
{
  vofi_int i, np0, nm0;
  dir_data ipsc;
  
  np0 = nm0 = 0;
  ipsc.ind1 = ipsc.ind2 = ipsc.swt1 = ipsc.swt2 = 0;
  ipsc.consi = 0;
  
  /* Count positive and negative values at vertices */
  for (i=0; i<8; i++) {
    if (fve[i] > 0.)
      np0++;
    else if (fve[i] < 0.)
      nm0++;
  }
  
  /* If all positive or all negative, no interface intersection */
  if (np0*nm0 == 0) {
    ipsc.consi = 0;
  }
  else {
    /* Interface intersects the face */
    if (np0 > nm0)
      ipsc.consi = 1;    /* Mostly positive */
    else
      ipsc.consi = -1;   /* Mostly negative */
  }
  
  return ipsc;
}

/* -------------------------------------------------------------------------- */
vofi_int vofi_get_face_min_3D(integrand impl_func, vofi_void_cptr par,
                             vofi_creal x0[], vofi_creal h0[],
                             vofi_creal dir1[], vofi_creal dir2[], vofi_creal dir3[],
                             vofi_creal fve[], min_data *xfsl, dir_data ipsc)
{
  vofi_int i, sign_change;
  vofi_real v0[NDIM], v1[NDIM], v2[NDIM], v3[NDIM];
  
  sign_change = 1;  /* Assume intersection exists */
  
  /* Compute three orthogonal directions on the 3D hyperface */
  for (i=0; i<NDIM; i++) {
    v1[i] = dir1[i] * h0[1];
    v2[i] = dir2[i] * h0[2];
    v3[i] = dir3[i] * h0[3];
  }
  
  /* Find the minimum on the face - simplified version */
  /* In a real implementation, this would use gradient descent or similar */
  for (i=0; i<NDIM; i++) {
    xfsl->xval[i] = x0[i] + 0.5*v1[i] + 0.5*v2[i] + 0.5*v3[i];
  }
  
  xfsl->fval = impl_func(xfsl->xval, par);
  xfsl->sval = 0.5;
  
  /* Check if actual sign change occurs at the minimum */
  if (xfsl->fval * ipsc.consi > 0.0) {
    sign_change = 0;  /* No actual zero-crossing */
  }
  
  return sign_change;
}