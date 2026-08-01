/**
 * @file      getlimits4d.c
 * @brief     Limits of integration in 4D: along udir for the external
 *            integration, and the analysis of one hyperplane u = const for
 *            the internal one.
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
 * vofi_get_limits_4D is vofi_get_limits_3D one level up, statement for       *
 * statement: walk the eight cell edges parallel to udir collecting the       *
 * crossings, and wherever the interface heads into the cell from one of      *
 * them, or a minimum was found inside a 3-face normal to pdir, add the       *
 * extrema of u on that interface. The first kind are internal limits         *
 * (basei = 1), the second external (basei = 0), and the same reorder and     *
 * zero-length removal close the routine.                                     *
 *                                                                            *
 * vofi_check_hyperplane is the 4D counterpart of vofi_check_plane: it takes  *
 * one hyperplane u = const, decides whether the 3D cross-section there is    *
 * full, empty or cut, and if cut returns its subdivision along tdir. The     *
 * frame is NOT re-ordered for the hyperplane -- it inherits (pdir,sdir,tdir) *
 * from the cell, which is why the direction-generic vofi_check_boundary_cell *
 * is used instead of vofi_check_boundary_surface.                            *
 * INPUT:  pointer to the implicit function, arrays of function parameters    *
 *         par, coordinates of minimum vertex x0, cell edges h0, the 16       *
 *         vertex values f0 in (p,s,t,u) order, min_data table xfsp, ordered  *
 *         directions                                                         *
 * OUTPUT: number of subdivisions and the array of subdivisions               *
 * FUNCTIONS:                                                                 *
 * vofi_int vofi_get_limits_4D  : external limits, along udir                 *
 * vofi_int vofi_check_hyperplane: cross-section type and its tdir limits     *
 * -------------------------------------------------------------------------- */

#define VIDX4(m,n,q,r) ((m) + 2*(n) + 4*(q) + 8*(r))

vofi_int vofi_get_limits_4D(integrand impl_func,vofi_void_cptr par,
                            vofi_creal x0[],vofi_creal h0[],vofi_creal f0[],
                            min_data xfsp[],vofi_real base[],vofi_creal pdir[],
                            vofi_creal sdir[],vofi_creal tdir[],
                            vofi_creal udir[])
{
  vofi_int basei[NSEG],i,l,l0,m,n,q,nsub,inters,consi;
  vofi_real xp[NDIM],xs[NDIM],xu[NDIM],fse[NSE],hs,ht,hu;
  min_data xfsl;

  /* Capacity: each of the eight u edges can contribute two crossings and
     two external limits, and each of the two 3-faces normal to pdir two
     more, so nsub cannot exceed 8*4 + 2*2 + 2 = 38 against NSEG = 48. The
     guards below are therefore unreachable today and are kept only so
     that raising the edge count can never silently truncate.          */
  base[0] = 0.;
  basei[0] = nsub = 1;
  hs = ht = hu = 0.;
  for (i=0;i<NDIM;i++) {
    hs += sdir[i]*h0[i];
    ht += tdir[i]*h0[i];
    hu += udir[i]*h0[i];
  }

  for (m=0;m<NSE;m++) {
    for (l=0;l<NDIM;l++)
      xp[l] = x0[l] + m*pdir[l]*h0[l];

    for (n=0;n<NSE;n++)
      for (q=0;q<NSE;q++) {
        l0 = 4*m + 2*n + q;
        if (xfsp[l0].isc[1] == 0 || nsub > NSEG-6)
          continue;
        for (l=0;l<NDIM;l++)
          xs[l] = xp[l] + n*sdir[l]*h0[l] + q*tdir[l]*h0[l];
        fse[0] = f0[VIDX4(m,n,q,0)];
        fse[1] = f0[VIDX4(m,n,q,1)];
        inters = vofi_get_side_intersections(impl_func,par,fse,xs,xfsp[l0],
                                 base,udir,hu,nsub,xfsp[l0].isc[1]);
        nsub += inters;
        for (i=1;i<=inters;i++)
          basei[nsub-i] = 1;

        /* does the reference phase continue into the 3-face from that
           crossing? If it does, the interface has an extremum in u
           somewhere inside and it has to bound a sector.               */
        memset(&xfsl,0,sizeof(xfsl));
        for (l=0;l<NDIM;l++)
          xu[l] = xs[l] + base[nsub-1]*udir[l];
        consi = vofi_check_line_consistency(impl_func,par,xu,sdir,hs,n,&xfsl);
        if (consi == 0)
          consi = vofi_check_line_consistency(impl_func,par,xu,tdir,ht,q,&xfsl);
        if (inters > 1 && consi == 0) {
          for (l=0;l<NDIM;l++)
            xu[l] = xs[l] + base[nsub-2]*udir[l];
          consi = vofi_check_line_consistency(impl_func,par,xu,sdir,hs,n,&xfsl);
          if (consi == 0)
            consi = vofi_check_line_consistency(impl_func,par,xu,tdir,ht,q,
                                                &xfsl);
        }
        if (consi > 0) {
          inters = vofi_get_uext_intersections(impl_func,par,xp,h0,xfsl,base,
                                               sdir,tdir,udir,nsub);
          nsub += inters;
          for (i=1;i<=inters;i++)
            basei[nsub-i] = 0;
        }
      }

    if (xfsp[IXC4].isc[m+1] != 0 && nsub <= NSEG-4) {
      inters = vofi_get_uext_intersections(impl_func,par,xp,h0,xfsp[IXC4],base,
                                           sdir,tdir,udir,nsub);
      nsub += inters;
      for (i=1;i<=inters;i++)
        basei[nsub-i] = 0;
    }
  }

  base[nsub] = hu;
  basei[nsub] = 1;

  /* - */
  vofi_reorder(base,basei,nsub);

  /* - */
  nsub = vofi_rm_segs(base,basei,nsub);

  return nsub;
}

/* -------------------------------------------------------------------------- */
vofi_int vofi_check_hyperplane(integrand impl_func,vofi_void_cptr par,
                          vofi_creal x1[],vofi_creal h0[],vofi_creal pdir[],
                          vofi_creal sdir[],vofi_creal tdir[],
                          vofi_real base_t[],min_data xfs3[],vofi_int *ityp)
{
  vofi_int i,m,n,q,np0,nm0,nneg,check_dir,n0[NSE][NSE][NSE];
  vofi_real f0s[NSE][NSE][NSE],xv[NDIM],fgrad[NDIM3];
  vofi_real hm,fth,fgradmod,fgradsq,fp,fmi;
  vofi_creal MIN_GRAD=1.0e-04;
  vofi_creal *dirs[NDIM3];

  dirs[0] = pdir; dirs[1] = sdir; dirs[2] = tdir;
  np0 = nm0 = 0;
  for (m=0;m<NSE;m++)
    for (n=0;n<NSE;n++)
      for (q=0;q<NSE;q++) {
        for (i=0;i<NDIM;i++)
          xv[i] = x1[i] + (m*pdir[i] + n*sdir[i] + q*tdir[i])*h0[i];
        f0s[m][n][q] = impl_func(xv,par);
        if (f0s[m][n][q] > 0.)
          np0++;
        else if (f0s[m][n][q] < 0.)
          nm0++;
      }

  /* the 3D threshold, in the frame the cell ordered */
  fgradsq = 0.;
  hm = 0.;
  for (i=0;i<NDIM3;i++) {
    fp = fmi = 0.;
    for (m=0;m<NSE;m++)
      for (n=0;n<NSE;n++)
        for (q=0;q<NSE;q++) {
          if ((i == 0 ? m : (i == 1 ? n : q)) == 1)
            fp += f0s[m][n][q];
          else
            fmi += f0s[m][n][q];
        }
    fgrad[i] = 0.;
    for (m=0;m<NDIM;m++)
      if (dirs[i][m] > 0.5 && h0[m] > 0.) {
        fgrad[i] = 0.25*(fp - fmi)/h0[m];
        hm = MAX(hm,0.5*h0[m]);
      }
    fgradsq += fgrad[i]*fgrad[i];
  }
  fgradmod = MAX(sqrt(fgradsq),MIN_GRAD);
  fth = sqrt(2.)*fgradmod*hm;

  *ityp = -1;
  if (np0*nm0 == 0) {
    nneg = nm0;
    np0 = nm0 = 0;
    for (m=0;m<NSE;m++)
      for (n=0;n<NSE;n++)
        for (q=0;q<NSE;q++) {
          if (fabs(f0s[m][n][q]) > fth) {
            n0[m][n][q] = 0;
            if (f0s[m][n][q] < 0.)
              nm0++;
            else
              np0++;
          }
          else
            n0[m][n][q] = 1;
        }

    if (nm0 == NVERC) {
      *ityp = 1;
      return 0;
    }
    else if (np0 == NVERC) {
      *ityp = 0;
      return 0;
    }

    check_dir = vofi_check_boundary_cell(impl_func,par,x1,h0,f0s,xfs3,n0,
                                         pdir,sdir,tdir);
      /* nm0 and np0 were overwritten by the recount above: if NO vertex is
         confident they are both zero, and the old test picked empty
         whatever the sign really was -- which reported a cell lying
         entirely inside the reference phase as empty whenever fth
         exceeded |f| at every vertex (reachable on strongly anisotropic
         cells). Decide on the RAW sign instead; this branch is only
         entered when there is no sign change, so it is unambiguous. */
    if (check_dir < 0) {
      *ityp = (nneg > 0) ? 1 : 0;
      return 0;
    }
  }

  vofi_check_secter_face(impl_func,par,x1,h0,pdir,sdir,tdir,f0s,&xfs3[4],fth);
  vofi_check_tertiary_side(impl_func,par,x1,h0,pdir,sdir,tdir,f0s,xfs3,fth);

  return vofi_get_limits_3D(impl_func,par,x1,h0,f0s,xfs3,base_t,pdir,sdir,
                            tdir);
}
