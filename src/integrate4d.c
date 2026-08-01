/**
 * @file      integrate4d.c
 * @brief     Cut hypervolume of a 4D cell (3D Gauss-Legendre integration),
 *            its centroid, and the 3-volume of the interface with its
 *            centroid.
 * @version   Vofi 2.0
 * @copyright GPLv3 license
 **/

#include "vofi_stddecl.h"
#include "vofi_GL_nodes.h"
#include "vofi_GL_weights.h"

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
 * The outermost of the four nested quadratures. In each sector of the u axis *
 * the cross-section is uniformly full, empty or cut, so V(u) is smooth there *
 * and one Gauss-Legendre rule integrates it; each node is a 3D problem       *
 * solved by vofi_get_volume IN THE FRAME THE CELL ORDERED, so the nesting is *
 *   u (here) -> t (vofi_get_volume) -> s (vofi_get_area) -> root along p.    *
 *                                                                            *
 * vofi_meas_add is called from the innermost loop, once per height point.    *
 * Because the whole cell shares one frame the interface is the graph         *
 * x_p = H(s,t,u), and its 3-volume is                                        *
 *      int sqrt(1 + |grad H|^2) ds dt du  =  int |grad f| / |df/dp|,         *
 * so the measure is a weighted sum over the very nodes the hypervolume       *
 * already visits. The integrand is unbounded where the interface turns over  *
 * in p; that happens at the EDGE of each sector, which the limits routines   *
 * put on a sector boundary and where Gauss-Legendre never samples, so the    *
 * sum stays finite but converges more slowly there than the volume does.     *
 * INPUT:  pointer to the implicit function, arrays of function parameters    *
 *         par, coordinates of minimum vertex x0, cell edges h0, sweep        *
 *         subdivisions base, the four ordered directions, flags nex, user's  *
 *         number of points npt, number of sectors nsub, tentative number of  *
 *         integration points nptmp                                           *
 * OUTPUT: cut hypervolume; centroid[0..3] the volume-weighted first moment   *
 *         in the (p,s,t,u) frame, centroid[4] the interface measure and      *
 *         centroid[5..8] its measure-weighted first moment in the same frame *
 * FUNCTIONS:                                                                 *
 * void vofi_meas_add            : one height point's contribution            *
 * vofi_real vofi_get_hypervolume: 4D problem, 3D integration                 *
 * -------------------------------------------------------------------------- */
void vofi_meas_add(meas_acc *macc,vofi_creal x20[],vofi_creal pdir[],
                   vofi_creal w,vofi_creal hp,vofi_creal ht,vofi_creal xt,
                   vofi_cint f_sign)
{
  vofi_int a;
  vofi_real pl,xpt[NDIM],xf[NDIM],xb[NDIM],grad[NDIM];
  vofi_real gn,gp,r,loc[NDIM];
  vofi_creal RMAX = 1./EPS_M;      /* a true tangency contributes no area */

  /* the height is measured from the far face when f_sign < 0 */
  pl = (f_sign > 0) ? ht : hp - ht;
  for (a=0;a<NDIM;a++)
    xpt[a] = x20[a] + pdir[a]*pl;

  gn = 0.;
  for (a=0;a<NDIM;a++) {
    memcpy(xf,xpt,sizeof(xf));
    memcpy(xb,xpt,sizeof(xb));
    xf[a] += macc->dh[a];
    xb[a] -= macc->dh[a];
    grad[a] = 0.5*(macc->func(xf,macc->par) - macc->func(xb,macc->par))/
              macc->dh[a];
    gn += grad[a]*grad[a];
  }
  gn = sqrt(gn);
  gp = fabs(grad[macc->jp]);
  r = (gp > gn/RMAX) ? gn/gp : RMAX;

  loc[0] = pl; loc[1] = xt; loc[2] = macc->tloc; loc[3] = macc->uloc;
  gn = w*macc->wout*r;
  macc->meas += gn;
  for (a=0;a<NDIM;a++)
    macc->mom[a] += gn*loc[a];

  return;
}

/* -------------------------------------------------------------------------- */
vofi_real vofi_get_hypervolume(integrand impl_func,vofi_void_cptr par,
                        vofi_creal x0[],vofi_creal h0[],vofi_creal base[],
                        vofi_creal pdir[],vofi_creal sdir[],vofi_creal tdir[],
                        vofi_creal udir[],vofi_real centroid[],vofi_cint nex[],
                        vofi_cint npt[],vofi_cint nsub,vofi_cint nptmp)
{
  vofi_int a,i,j,k,nu,nsub_t,ityp,nexpt;
  vofi_int nex3[2],nvis3[2]={0,0};
  vofi_real x1[NDIM],base_t[NSEG],cent3[2*NDIM+1],mom[NDIM],quadm[NDIM];
  vofi_real hp,hs,ht,hu,hm,volume,vol,v3,du,mdpt,xu,quadv;
  vofi_creal *ptw,*ptx;
  min_data xfs3[NXFS];
  meas_acc macc,*pmacc;

  hp = hs = ht = hu = 0.;
  for (a=0;a<NDIM;a++) {
    hp += pdir[a]*h0[a];
    hs += sdir[a]*h0[a];
    ht += tdir[a]*h0[a];
    hu += udir[a]*h0[a];
  }
  hm = h0[0];
  for (a=1;a<NDIM;a++)
    hm = MAX(hm,h0[a]);

  nex3[0] = nex[0];
  nex3[1] = 0;                      /* no 3D triangulation inside a slice */
  pmacc = NULL;
  if (nex[1] > 0) {
    memset(&macc,0,sizeof(macc));
    macc.func = impl_func;
    macc.par  = (void *) par;
    macc.jp   = 0;
    for (a=0;a<NDIM;a++) {
      if (pdir[a] > 0.5)
        macc.jp = a;
      macc.dh[a] = MAX(1.e-04*h0[a],EPS_ROOT);
    }
    pmacc = &macc;
  }

  volume = 0.;
  for (a=0;a<NDIM;a++)
    mom[a] = 0.;

  for (nu=1;nu<=nsub;nu++) {                            /* - */
    du = base[nu] - base[nu-1];
    mdpt = 0.5*(base[nu] + base[nu-1]);
    for (a=0;a<NDIM;a++)
      x1[a] = x0[a] + udir[a]*mdpt;
    memset(xfs3,0,sizeof(xfs3));
    nsub_t = vofi_check_hyperplane(impl_func,par,x1,h0,pdir,sdir,tdir,base_t,
                                   xfs3,&ityp);

    if (ityp == 0)                                      /* - */
      ;
    else if (ityp > 0) {                                /* - */
      vol = du*hp*hs*ht;
      volume += vol;
      if (nex[0] > 0) {
        mom[0] += vol*0.5*hp;
        mom[1] += vol*0.5*hs;
        mom[2] += vol*0.5*ht;
        mom[3] += vol*mdpt;
      }
    }
    else {                                              /* - */
      nexpt = (vofi_int) (18.*du/hm);
      nexpt = MIN(NGLM,nexpt+3);
      if (npt[5] >= 3 && npt[5] <= NGLM)
        nexpt = MIN(npt[5],nexpt);
      if (npt[4] >= 3 && npt[4] <= NGLM)
        nexpt = MAX(npt[4],nexpt);
      nexpt = MAX(3,MIN(NGLM,nexpt));
      j = nexpt - 3;
      ptx = (double *) csipt[j];
      ptw = (double *) wgtpt[j];

      quadv = 0.;
      for (a=0;a<NDIM;a++)
        quadm[a] = 0.;
      for (k=1;k<=nexpt;k++) {
        xu = mdpt + 0.5*du*(*ptx);
        for (a=0;a<NDIM;a++)
          x1[a] = x0[a] + udir[a]*xu;
        memset(xfs3,0,sizeof(xfs3));
        memset(cent3,0,sizeof(cent3));
        nsub_t = vofi_check_hyperplane(impl_func,par,x1,h0,pdir,sdir,tdir,
                                       base_t,xfs3,&ityp);
        if (ityp == 0)
          v3 = 0.;
        else if (ityp > 0) {
          v3 = hp*hs*ht;
          cent3[0] = v3*0.5*hp;
          cent3[1] = v3*0.5*hs;
          cent3[2] = v3*0.5*ht;
        }
        else {
          if (pmacc != NULL) {
            pmacc->wout = 0.5*du*(*ptw);
            pmacc->uloc = xu;
          }
          v3 = vofi_get_volume(impl_func,par,x1,h0,base_t,pdir,sdir,tdir,
                               cent3,nex3,npt,nsub_t,nptmp,nvis3,pmacc);
        }
        quadv += (*ptw)*v3;
        if (nex[0] > 0) {
          for (i=0;i<NDIM3;i++)
            quadm[i] += (*ptw)*cent3[i];
          quadm[3] += (*ptw)*v3*xu;
        }
        ptx++;
        ptw++;
      }
      volume += 0.5*du*quadv;
      if (nex[0] > 0)
        for (a=0;a<NDIM;a++)
          mom[a] += 0.5*du*quadm[a];
    }
  }

  for (a=0;a<NDIM;a++)
    centroid[a] = mom[a];
  centroid[NDIM] = (pmacc != NULL) ? macc.meas : 0.;
  for (a=0;a<NDIM;a++)
    centroid[NDIM+1+a] = (pmacc != NULL) ? macc.mom[a] : 0.;

  return volume;
}
