/**
 * @file      orderdirs4d.c
 * @brief     4D cell type and, if the cell is cut, the ordered coordinate
 *            directions pdir, sdir, tdir, udir and the number of integration
 *            points along the secondary direction.
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
 * a) evaluate f at the 16 vertices and the threshold value fth,              *
 * b) with fth, consistency and minimum checks, decide full / empty / cut,    *
 * c) if cut, order ALL FOUR coordinate directions and compute the tentative  *
 *    number of integration points.                                           *
 *                                                                            *
 * The ordering is by |grad|, largest first: pdir carries the height, udir is *
 * the outermost sweep. It is chosen ONCE for the whole cell and every        *
 * hyperplane inherits it, which is what 2D does with (p,s) across its        *
 * sectors and 3D with (p,s,t) across its planes. A single frame is also what *
 * makes the interface a single graph x_p = H(s,t,u) over the cell, so the    *
 * interface measure can be integrated on the very nodes the volume uses.     *
 * As in 3D, a minimum found inside a 3-face overrides the gradient and       *
 * becomes the primary direction.                                             *
 * INPUT:  pointer to the implicit function, arrays of function parameters    *
 *         par, coordinates of minimum vertex x0, cell edges h0               *
 * OUTPUT: icc (1 full / 0 empty / -1 cut), pdir, sdir, tdir, udir, the 16    *
 *         vertex values f0 REORDERED into (p,s,t,u), and the min_data table  *
 * FUNCTIONS:                                                                 *
 * vofi_int vofi_cell_type_4D : cell type only                                *
 * vofi_int vofi_order_dirs_4D: cell type, ordered directions, point count    *
 * -------------------------------------------------------------------------- */

/* f at the 16 vertices in GLOBAL axis order, bit a = index along axis a */
static void vofi_hyper_vertices(integrand impl_func,vofi_void_cptr par,
                                vofi_creal x0[],vofi_creal h0[],vofi_real f0[],
                                vofi_int *np0,vofi_int *nm0)
{
  vofi_int a,v;
  vofi_real x1[NDIM];

  *np0 = *nm0 = 0;
  for (v=0;v<NVERH;v++) {
    for (a=0;a<NDIM;a++)
      x1[a] = x0[a] + ((v >> a) & 1)*h0[a];
    f0[v] = impl_func(x1,par);
    if (f0[v] > 0.)
      (*np0)++;
    else if (f0[v] < 0.)
      (*nm0)++;
  }

  return;
}

/* -------------------------------------------------------------------------- */
/* the shared full/empty/cut decision; f0 and fgrad come back filled, and
   xfs[0..3] carry, per axis, any minimum found inside a 3-face            */
static vofi_int vofi_hyper_type(integrand impl_func,vofi_void_cptr par,
                                vofi_creal x0[],vofi_creal h0[],vofi_real f0[],
                                vofi_real fgrad[],min_data xfs[],
                                vofi_int *check_dir)
{
  vofi_int v,np0,nm0,nneg,icc,n0[NVERH];
  vofi_real fth;

  icc = -1;
  *check_dir = -1;
  vofi_hyper_vertices(impl_func,par,x0,h0,f0,&np0,&nm0);
  fth = vofi_hyper_fth(f0,h0,fgrad);

  if (np0*nm0 == 0) {
    /* no sign change on the vertices: only a vertex whose |f| exceeds fth
       can rule out an interface, the others have to be probed           */
    nneg = nm0;
    np0 = nm0 = 0;
    for (v=0;v<NVERH;v++) {
      if (fabs(f0[v]) > fth) {
        n0[v] = 0;
        if (f0[v] < 0.)
          nm0++;
        else
          np0++;
      }
      else
        n0[v] = 1;
    }

    if (nm0 == NVERH)
      return 1;
    else if (np0 == NVERH)
      return 0;

    *check_dir = vofi_check_boundary_hyper(impl_func,par,x0,h0,f0,xfs,n0);

      /* nm0 and np0 were overwritten by the recount above: if NO vertex is
         confident they are both zero, and the old test picked empty
         whatever the sign really was -- which reported a cell lying
         entirely inside the reference phase as empty whenever fth
         exceeded |f| at every vertex (reachable on strongly anisotropic
         cells). Decide on the RAW sign instead; this branch is only
         entered when there is no sign change, so it is unambiguous. */
    if (*check_dir < 0)
      icc = (nneg > 0) ? 1 : 0;
  }

  return icc;
}

/* -------------------------------------------------------------------------- */
vofi_int vofi_cell_type_4D(integrand impl_func,vofi_void_cptr par,
                           vofi_creal x0[],vofi_creal h0[])
{
  vofi_int check_dir;
  vofi_real f0[NVERH],fgrad[NDIM];
  min_data xfsl[NDIM];

  memset(xfsl,0,sizeof(xfsl));

  return vofi_hyper_type(impl_func,par,x0,h0,f0,fgrad,xfsl,&check_dir);
}

/* -------------------------------------------------------------------------- */
/* curvature of the interface in the (ja,jb) plane through xc, on a 3-point
   stencil of half-cell spacing: the 2D estimate of vofi_order_dirs_2D. The
   weight is the reciprocal of |sum f| over the stencil, so that sections
   passing close to the interface dominate the average, as in 3D.           */
static vofi_real vofi_kappa_plane(integrand impl_func,vofi_void_cptr par,
                                  vofi_creal xc[],vofi_creal hh[],vofi_cint ja,
                                  vofi_cint jb,vofi_real *wgt)
{
  vofi_int i,j,a;
  vofi_real fd[NSTC][NSTC],x1[NDIM],sumf,have,tmp;
  vofi_real fx,fy,fxx,fyy,fxy;

  sumf = 0.;
  for (i=0;i<NSTC;i++)
    for (j=0;j<NSTC;j++) {
      for (a=0;a<NDIM;a++)
        x1[a] = xc[a];
      x1[ja] += (i-1)*hh[ja];
      x1[jb] += (j-1)*hh[jb];
      fd[i][j] = impl_func(x1,par);
      sumf += fd[i][j];
    }

  have = hh[ja] + hh[jb];                    /* = 0.5*(h0[ja] + h0[jb]) */
  fx  = (fd[2][1] - fd[0][1])*have/(2.*hh[ja]);
  fy  = (fd[1][2] - fd[1][0])*have/(2.*hh[jb]);
  fxx = (fd[2][1] + fd[0][1] - 2.*fd[1][1])*have*have/(hh[ja]*hh[ja]);
  fyy = (fd[1][2] + fd[1][0] - 2.*fd[1][1])*have*have/(hh[jb]*hh[jb]);
  fxy = (fd[2][2] - fd[2][0] - fd[0][2] + fd[0][0])*have*have/
        (4.*hh[ja]*hh[jb]);
  tmp = fx*fx + fy*fy;
  tmp = sqrt(tmp*tmp*tmp) + EPS_NOT0;
  *wgt = 1./MAX(fabs(sumf),EPS_NOT0);

  return fabs(fxx*fy*fy - 2.*fx*fy*fxy + fx*fx*fyy)/tmp;
}

/* -------------------------------------------------------------------------- */
vofi_int vofi_order_dirs_4D(integrand impl_func,vofi_void_cptr par,
                            vofi_creal x0[],vofi_creal h0[],vofi_real pdir[],
                            vofi_real sdir[],vofi_real tdir[],vofi_real udir[],
                            vofi_real f0[],min_data xfsp[])
{
  vofi_int a,k,m,n,q,r,icc,check_dir,jp,js,jt,ju,ord[NDIM],npt,itmp;
  vofi_real fgrad[NDIM],hh[NDIM],xc[NDIM],fp[NVERH];
  vofi_real fth,Kappa,wsum,w,tmp;
  vofi_creal a0=2.34607, a1=16.5515, a2=-5.53054, a3=54.0866;

  icc = vofi_hyper_type(impl_func,par,x0,h0,f0,fgrad,xfsp,&check_dir);
  if (icc >= 0)
    return icc;

  /* order the axes by |grad|, largest first */
  for (a=0;a<NDIM;a++) {
    fgrad[a] = fabs(fgrad[a]);
    ord[a] = a;
  }
  for (a=1;a<NDIM;a++) {
    k = a;
    while (k > 0 && fgrad[ord[k-1]] < fgrad[ord[k]]) {
      itmp = ord[k-1]; ord[k-1] = ord[k]; ord[k] = itmp;
      k--;
    }
  }
  jp = ord[0]; js = ord[1]; jt = ord[2]; ju = ord[3];

  /* a minimum inside a 3-face overrides the gradient: that 3-face has to be
     normal to pdir for vofi_check_secter_cell to see it again, exactly as
     in vofi_order_dirs_3D                                                 */
  if (check_dir >= 0 && xfsp[jp].isc[0] != 1) {
    for (a=0;a<NDIM;a++)
      if (ord[a] == check_dir) {
        for (k=a;k>0;k--)
          ord[k] = ord[k-1];
        ord[0] = check_dir;
        break;
      }
    jp = ord[0]; js = ord[1]; jt = ord[2]; ju = ord[3];
  }

  for (a=0;a<NDIM;a++)
    pdir[a] = sdir[a] = tdir[a] = udir[a] = 0.;
  pdir[jp] = sdir[js] = tdir[jt] = udir[ju] = 1.;

  /* reorder the vertex values from global axes into (p,s,t,u) */
  fth = vofi_hyper_fth(f0,h0,fgrad);
  for (m=0;m<NSE;m++)
    for (n=0;n<NSE;n++)
      for (q=0;q<NSE;q++)
        for (r=0;r<NSE;r++)
          fp[m + 2*n + 4*q + 8*r] =
            f0[(m << jp) | (n << js) | (q << jt) | (r << ju)];
  for (k=0;k<NVERH;k++)
    f0[k] = fp[k];

  /* the minima that set the limits of integration along udir */
  if (check_dir >= 0 && xfsp[jp].isc[0] == 1)
    xfsp[IXC4] = xfsp[jp];
  else
    vofi_check_secter_cell(impl_func,par,x0,h0,pdir,sdir,tdir,udir,f0,
                           &xfsp[IXC4],fth);
  vofi_check_quaternary_side(impl_func,par,x0,h0,pdir,sdir,tdir,udir,f0,
                             xfsp,fth);

  /* tentative number of integration points, from the interface curvature in
     the (p,s) plane averaged over the three levels of tdir                */
  for (a=0;a<NDIM;a++) {
    hh[a] = 0.5*h0[a];
    xc[a] = x0[a] + hh[a];
  }
  Kappa = wsum = 0.;
  for (k=0;k<NSTC;k++) {
    xc[jt] = x0[jt] + k*hh[jt];
    Kappa += vofi_kappa_plane(impl_func,par,xc,hh,jp,js,&w)*w;
    wsum += w;
  }
  Kappa = Kappa/MAX(wsum,EPS_NOT0);
  tmp = a0 + Kappa*(a1 + Kappa*(a2 + a3*Kappa));
  npt = (vofi_int) ceil(tmp);
  xfsp[IXC4].ipt = MAX(4,MIN(npt,NGLM));

  return icc;
}

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * threshold below which the function value at a vertex is too small to rule  *
 * out an interface elsewhere in the cell. 2D uses |grad|*hmax and 3D         *
 * sqrt(2)*|grad|*hmax, i.e. the distance from a vertex to the centre of the  *
 * cell measured in half-edges; in 4D that factor is sqrt(3). The gradient is *
 * the average of the eight opposite-vertex differences along each axis and   *
 * is returned in fgrad, which the caller also uses to order the directions.  *
 * -------------------------------------------------------------------------- */
vofi_real vofi_hyper_fth(vofi_creal f0[],vofi_creal h0[],vofi_real fgrad[])
{
  vofi_int a,v;
  vofi_real fp,fm,fgradsq,fgradmod,hm;
  vofi_creal MIN_GRAD=1.0e-04;

  fgradsq = 0.;
  hm = 0.;
  for (a=0;a<NDIM;a++) {
    fp = fm = 0.;
    for (v=0;v<NVERH;v++) {
      if (v & (1 << a))
        fp += f0[v];
      else
        fm += f0[v];
    }
    fgrad[a] = (h0[a] > 0.) ? 0.125*(fp - fm)/h0[a] : 0.;
    fgradsq += fgrad[a]*fgrad[a];
    hm = MAX(hm,0.5*h0[a]);
  }
  fgradmod = MAX(sqrt(fgradsq),MIN_GRAD);

  return sqrt(3.)*fgradmod*hm;
}
