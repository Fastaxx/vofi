/**
 * @file      getintersections4d.c
 * @brief     External limits of integration in 4D: the extrema along udir of
 *            the interface inside a 3-face normal to pdir.
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
 * This is vofi_get_ext_intersections one dimension up, and it works the same  *
 * way. In 3D the reference phase inside a face normal to pdir is a 2D region; *
 * the routine tracks the MIDPOINT OF ITS CHORD along the secondary direction  *
 * while marching along the tertiary one, and the tertiary extremum is where   *
 * the chord shrinks to nothing. In 4D the region inside a 3-face is a volume, *
 * so the chord midpoint becomes the centre of the two chords through the      *
 * point, along sdir and along tdir, and the march runs along udir. The stop   *
 * condition is unchanged: the section has pinched off.                        *
 *                                                                            *
 * Marching on the centre rather than on the interface itself is what makes    *
 * this robust near the extremum, where the interface turns over and any       *
 * height function in u becomes singular: the centre stays in the interior     *
 * and its motion, not a derivative, sets the next direction.                  *
 *                                                                            *
 * INPUT:  pointer to the implicit function, arrays of function parameters    *
 *         par, corner x0 of the 3-face and cell edges h0, min_data xfsp       *
 *         carrying a point where the reference phase is present, current      *
 *         subdivision count nsub, the ordered directions sdir, tdir, udir     *
 * OUTPUT: number of new limits (2), updated subdivision array base            *
 * FUNCTIONS:                                                                 *
 * vofi_int vofi_get_uext_intersections                                        *
 * -------------------------------------------------------------------------- */

/* distance from pt to the root along dir, given the raw function values fa at
   pt and fb at pt + len*dir; the sign convention is the one every other
   caller of vofi_get_segment_zero uses                                       */
static vofi_real vofi_edge_root(integrand impl_func,vofi_void_cptr par,
                                vofi_creal pt[],vofi_creal dir[],
                                vofi_creal len,vofi_creal fa,vofi_creal fb)
{
  vofi_int f2neg;
  vofi_real s0[4],ds;

  f2neg = (fa < 0.) ? 1 : -1;
  s0[0] = len;
  if (fabs(fa) < fabs(fb)) {
    s0[1] = 0.;  s0[2] = fa;
  }
  else {
    s0[1] = len; s0[2] = fb;
  }
  s0[3] = (fb - fa)/len;
  ds = vofi_get_segment_zero(impl_func,par,pt,dir,s0,f2neg);
  if (f2neg < 0)
    ds = len - ds;

  return ds;
}

/* -------------------------------------------------------------------------- */
/* the chord of the reference phase through pt along axis j, clipped to the
   cell; pt must be strictly inside (f2neg*f(pt) < 0). Returns its length.  */
static vofi_real vofi_chord_axis(integrand impl_func,vofi_void_cptr par,
                                 vofi_creal x0[],vofi_creal h0[],
                                 vofi_creal pt[],vofi_cint j,vofi_cint f2neg,
                                 vofi_real *lo,vofi_real *hi)
{
  vofi_int i,k,sgn;
  vofi_real xs[NDIM],dir[NDIM],fa,fb,len,bnd[NSE];

  fa = impl_func(pt,par);
  for (k=0;k<NSE;k++) {
    sgn = (k == 0) ? -1 : 1;
    len = (k == 0) ? pt[j] - x0[j] : x0[j] + h0[j] - pt[j];
    bnd[k] = pt[j] + sgn*len;
    if (len > EPS_ROOT) {
      for (i=0;i<NDIM;i++) {
        xs[i] = pt[i];
        dir[i] = 0.;
      }
      dir[j] = (vofi_real) sgn;
      xs[j] = bnd[k];
      fb = impl_func(xs,par);
      if (f2neg*fb > 0.)
        bnd[k] = pt[j] + sgn*vofi_edge_root(impl_func,par,pt,dir,len,fa,fb);
    }
  }
  *lo = bnd[0];
  *hi = bnd[1];

  return bnd[1] - bnd[0];
}

/* -------------------------------------------------------------------------- */
/* move the s and t coordinates of pc onto the centre of the reference phase,
   starting from a point that may sit exactly ON the interface: probe a
   tolerance either way first, as the 3D routine does. Returns 0 if no
   interior point was found along either axis, otherwise the larger of the
   two chord lengths in *csize.                                            */
static vofi_int vofi_recentre(integrand impl_func,vofi_void_cptr par,
                              vofi_creal x0[],vofi_creal h0[],vofi_real pc[],
                              vofi_cint jax[],vofi_cint f2neg,
                              vofi_real *csize)
{
  vofi_int i,j,jj,k,inside,ipt;
  vofi_real ptt[NDIM],lo,hi,chord;
  vofi_creal tol = EPS_M;

  ipt = 0;
  *csize = 0.;
  for (j=0;j<NSE;j++) {
    jj = jax[j];
    inside = 0;
    for (k=0;k<NSE && !inside;k++) {
      for (i=0;i<NDIM;i++)
        ptt[i] = pc[i];
      ptt[jj] += (k == 0) ? tol : -tol;
      if (ptt[jj] < x0[jj] || ptt[jj] > x0[jj] + h0[jj])
        continue;
      if (f2neg*impl_func(ptt,par) < 0.)
        inside = 1;
    }
    if (!inside)
      continue;
    ipt = 1;
    chord = vofi_chord_axis(impl_func,par,x0,h0,ptt,jj,f2neg,&lo,&hi);
    pc[jj] = 0.5*(lo + hi);
    *csize = MAX(*csize,chord);
  }

  return ipt;
}

/* -------------------------------------------------------------------------- */
vofi_int vofi_get_uext_intersections(integrand impl_func,vofi_void_cptr par,
              vofi_creal x0[],vofi_creal h0[],min_data xfsp,vofi_real base[],
              vofi_creal sdir[],vofi_creal tdir[],vofi_creal udir[],
              vofi_cint nsub)
{
  vofi_int i,k,a,iter,not_conv,ist,ipt,inters,f2neg,ju,jax[NSE];
  vofi_real pt0[NDIM],pt1[NDIM],mp0[NDIM],mp1[NDIM],ext_dir[NDIM],ss[NDIM];
  vofi_real lo,hi,csize,sst,ssy,fa,fb,normdir,d1,d2,a1,a2;
  vofi_creal tol2 = 2.*EPS_M;

  ju = jax[0] = jax[1] = 0;
  for (a=0;a<NDIM;a++) {
    if (udir[a] > 0.5) ju = a;
    if (sdir[a] > 0.5) jax[0] = a;
    if (tdir[a] > 0.5) jax[1] = a;
  }
  inters = 0;
  f2neg = (xfsp.fval < 0.) ? 1 : -1;

  for (i=0;i<NDIM;i++)
    pt0[i] = xfsp.xval[i];
  /* the seed can sit a rounding step outside; clip it onto the 3-face */
  for (k=0;k<NSE;k++) {
    a = jax[k];
    pt0[a] = MAX(pt0[a],x0[a]);
    pt0[a] = MIN(pt0[a],x0[a] + h0[a]);
  }
  pt0[ju] = MAX(pt0[ju],x0[ju]);
  pt0[ju] = MIN(pt0[ju],x0[ju] + h0[ju]);

  /* centre the seed inside the section it belongs to */
  for (k=0;k<NSE;k++) {
    a = jax[k];
    if (f2neg*impl_func(pt0,par) < 0.) {
      vofi_chord_axis(impl_func,par,x0,h0,pt0,a,f2neg,&lo,&hi);
      pt0[a] = 0.5*(lo + hi);
    }
  }

  for (k=-1;k<=1;k=k+2) {                       /* the two u directions */
    for (i=0;i<NDIM;i++) {
      ext_dir[i] = 0.;
      mp1[i] = pt0[i];
    }
    ext_dir[ju] = (vofi_real) k;
    sst = (k < 0) ? pt0[ju] - x0[ju] : x0[ju] + h0[ju] - pt0[ju];

    iter = 0;
    not_conv = 1;
    while (not_conv && iter < MAX_ITER_MINI) {
      /* march to the interface along the current direction, growing the
         step while the far end is still inside the reference phase      */
      fa = impl_func(mp1,par);
      ist = 0;
      for (i=0;i<NDIM;i++)
        pt1[i] = mp1[i] + sst*ext_dir[i];
      fb = impl_func(pt1,par);
      while (f2neg*fb < 0. && ist < 3 && sst > EPS_ROOT) {
        ssy = sst;
        for (i=0;i<NDIM;i++) {          /* how far can the box be crossed */
          d1 = SGN0P(ext_dir[i]);
          d2 = fabs(ext_dir[i]) + EPS_NOT0;
          if (d2 < EPS_ROOT)
            ss[i] = SS_FREE;
          else {
            a1 = (x0[i] - mp1[i])/(d1*d2);
            a2 = (x0[i] + h0[i] - mp1[i])/(d1*d2);
            ss[i] = MAX(a1,a2);
          }
        }
        ssy = ss[0];
        for (i=1;i<NDIM;i++)
          ssy = MIN(ssy,ss[i]);
        if (sst >= ssy)
          break;
        sst = MIN(3.*sst,ssy);
        for (i=0;i<NDIM;i++)
          pt1[i] = mp1[i] + sst*ext_dir[i];
        fb = impl_func(pt1,par);
        ist++;
      }
      if (f2neg*fb > 0. && sst > EPS_ROOT) {
        sst = vofi_edge_root(impl_func,par,mp1,ext_dir,sst,fa,fb);
        for (i=0;i<NDIM;i++)
          pt1[i] = mp1[i] + sst*ext_dir[i];
      }

      /* re-centre the section through the new point */
      for (i=0;i<NDIM;i++) {
        mp0[i] = mp1[i];
        mp1[i] = pt1[i];
      }
      ipt = vofi_recentre(impl_func,par,x0,h0,mp1,jax,f2neg,&csize);

      if (!ipt) {                       /* the section has pinched off */
        for (i=0;i<NDIM;i++)
          mp1[i] = pt1[i];
        not_conv = 0;
      }
      else {
        /* the next direction is the motion of the centre */
        normdir = 0.;
        for (i=0;i<NDIM;i++) {
          ext_dir[i] = mp1[i] - mp0[i];
          normdir += ext_dir[i]*ext_dir[i];
        }
        normdir = sqrt(normdir) + EPS_NOT0;
        for (i=0;i<NDIM;i++) {
          ext_dir[i] = ext_dir[i]/normdir;
          d1 = SGN0P(ext_dir[i]);
          d2 = fabs(ext_dir[i]) + EPS_NOT0;
          if (d2 < EPS_ROOT)
            ss[i] = SS_FREE;
          else {
            a1 = (x0[i] - mp1[i])/(d1*d2);
            a2 = (x0[i] + h0[i] - mp1[i])/(d1*d2);
            ss[i] = MAX(a1,a2);
          }
        }
        ssy = ss[0];
        for (i=1;i<NDIM;i++)
          ssy = MIN(ssy,ss[i]);
        sst = MIN(1.2*sst,ssy);
        if (csize < tol2 || sst < EPS_ROOT || normdir < EPS_ROOT)
          not_conv = 0;
      }
      iter++;
    }
    base[nsub+inters] = mp1[ju] - x0[ju];
    inters++;
  }

  return inters;
}
