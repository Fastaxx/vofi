/****************************************************************************
 * Copyright (C) 2021 by Andrea Chierici(a), Vincent Le Chenadec(b),        *
 * Ruben Scardovelli(a,f), Philip Yecko(c) and Stéphane Zaleski(d,e)        *
 * (a) DIN–Lab. di Montecuccolino, Università di Bologna,                   *
 *     Via dei Colli 16, 40136 Bologna, Italy                               *
 * (b) Campus de Marne-la-Vallée, Université Gustave Eiffel,                *
 *     77454, Marne-la-Vallée Cedex 2, France                               *
 * (c) Physics Department, Cooper Union, New York, NY, USA                  *
 * (d) Sorbonne Universités, UPMC Univ Paris 06, UMR 7190,                  *
 *     Institut Jean Le Rond d’Alembert, F-75005, Paris, France             *
 * (e) CNRS, UMR 7190, Institut Jean Le Rond d’Alembert, F-75005,           *
 *     Paris, France                                                        *
 * (f) e-mail: ruben.scardovelli@unibo.it                                   *
 *                                                                          *
 * Vofi library: version 2.0                                                *
 * You should have received a copy of the CPC license along with Vofi.      *
 * If not, see http://cpc.cs.qub.ac.uk/licence/licence.html.                *
 ****************************************************************************/
/**
 * @file      vofi.h
 * @authors   Andrea Chierici, Vincent Le Chenadec, Ruben Scardovelli, 
 *            Philip Yecko and Stéphane Zaleski 
 * @date      February 15, 2021
 * @brief     Header file for the Vofi library
 * @version   Vofi 2.0
 * @copyright CPC license
 **/
/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * header file for the Vofi library                                           *
 * -------------------------------------------------------------------------- */

#ifndef VOFI_H
#define VOFI_H

typedef       double  vofi_real;
typedef const double  vofi_creal;
typedef       int  vofi_int;
typedef const int  vofi_cint;
typedef int * const vofi_int_cpt;
typedef void * const vofi_void_cptr;
typedef double (*integrand) (vofi_creal [],vofi_void_cptr);

#ifdef __cplusplus
extern "C" {
#endif

/* ndim0 = 1, 2, 3 or 4. The caller's arrays must be sized for ndim0:
 *
 *   xin, h0   ndim0 reals (only ndim0 components are ever read)
 *   xex       4 reals for ndim0 <= 3 -- xex[0..2] centroid, xex[3] the
 *             interface measure -- and 5 for ndim0 == 4, where the centroid
 *             takes xex[0..3] and the measure moves to xex[4]
 *   npt       4 ints for ndim0 <= 3; 6 for ndim0 == 4, the extra pair
 *             npt[4], npt[5] being the min and max number of
 *             Gauss-Legendre points along the 4D sweep direction. Passing
 *             a shorter array in 4D reads out of bounds.
 *   nex, nvis 2 ints
 *
 * 4D NOTES. Everything is computed the same way as in 2D and 3D: the four
 * directions are ordered once for the cell, and one more level of nested
 * height-function quadrature is added, u -> t -> s -> root along p. The
 * hypervolume fraction, the 4-component centroid and the interface measure
 * are all exact to round-off for a hyperplane cut of a hypercube.
 *
 * xex[4] is the 3-VOLUME of the interface and xgam its centroid. Unlike the
 * 2D polyline and the 3D triangulation, it is not a sum over facets but the
 * integral of sqrt(1 + |grad H|^2) over the same nodes the hypervolume uses,
 * H being the height x_p = H(s,t,u); the two agree in the limit and this one
 * needs no tetrahedralisation. nvis is ignored in 4D.                    */
vofi_real vofi_get_cc(integrand,vofi_void_cptr,vofi_creal [],
                      vofi_creal [],vofi_real [],vofi_cint [],
                      vofi_cint [],vofi_cint [],vofi_cint);

/* as vofi_get_cc, with the INTERFACE CENTROID returned in the extra array
   (3 reals, or 4 when ndim0 == 4; only when nex[1] > 0 -- the cell centre
   otherwise). It is the centroid of the very same point set / polyline /
   triangulated surface / height graph whose measure vofi_get_cc returns in
   xex[3] (xex[4] in 4D) -- the identical quadrature, not a separate
   approximation. */
vofi_real vofi_get_cc_gam(integrand,vofi_void_cptr,vofi_creal [],
                          vofi_creal [],vofi_real [],vofi_real [],
                          vofi_cint [],vofi_cint [],vofi_cint [],vofi_cint);

/* 1 full / 0 empty / -1 cut, for ndim0 = 2, 3 or 4 */
 vofi_int vofi_get_cell_type(integrand,vofi_void_cptr,vofi_creal [],
			     vofi_creal [],vofi_cint);

#ifdef __cplusplus
}
#endif

#endif
