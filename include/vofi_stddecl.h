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
 * You should have received a copy of the CPC license along with Vofi       *
 * If not, see http://cpc.cs.qub.ac.uk/licence/licence.html.                *
 ****************************************************************************/
/**
 * @file      vofi_stddecl.h
 * @authors   Andrea Chierici, Vincent Le Chenadec, Ruben Scardovelli, 
 *            Philip Yecko and Stéphane Zaleski 
 * @date      February 15, 2021
 * @brief     Header file containing the functions prototype of the library
 * @version   Vofi 2.0
 * @copyright CPC license
 **/

#ifndef VOFI_STDDECL_H
#define VOFI_STDDECL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef __cplusplus
#define Extern extern "C"
#else
#define Extern extern
#endif

#define PREFIX(s) s

#if NOUNDERSCORE
#define SUFFIX(s) s
#else
#define SUFFIX(s) s##_
#endif

#define EXPORT(s) EXPORT_(PREFIX(s))
#define EXPORT_(s) SUFFIX(s)

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SGN0P(a) ((a<0) ? -1 : 1)
#define Sq(a) ((a)*(a))
#define Sq2(a) (a[0]*a[0] + a[1]*a[1])
#define Sq3(a) (a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
#define Sqd3(a,b) ((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))
#define SHFT4(a,b,c,d)  (a)=(b); (b)=(c); (c)=(d)
#define CPSF(s,t,f,g) (s)=(t); (f)=(g)

#define EPS_M    1.5e-07
#define EPS_LOC  1.5e-07
#define EPS_E    5.0e-07
#define EPS_SEGM 1.0e-12
#define EPS_ROOT 1.0e-14
#define EPS_NOT0 1.0e-90
#define NEAR_EDGE_RATIO 2.0e-02
#define MAX_ITER_ROOT 15
#define MAX_ITER_MINI 50  
/* NDIM is the STORAGE dimension: every coordinate and direction array in the
   library is this long, and the components above ndim0 are zero. It is 4 so
   that one set of direction-vector primitives serves 1D, 2D, 3D and 4D. Two
   constants that used to ride on NDIM but never meant "dimension" are now
   spelled out: NSTC is the width of the 3-point central stencils, NSCT the
   largest number of sectors a side can be cut into.                        */
#define NDIM   4
#define NDIM3  3     /* a genuine 3-vector (local triangle coordinates)     */
#define NSTC   3     /* points per axis in the central-difference stencils  */
#define NSCT   3     /* max sectors along a side handled by vofi_sector_new */
#define NVER   4     /* vertices of a 2D face                               */
#define NVERC  8     /* vertices of a 3D cell (a 3-face of a 4D cell)       */
#define NVERH 16     /* vertices of a 4D cell                               */
#define NSE    2
#define NSEG  48
#define NGLM  20
/* stand-in for "this direction is not constrained by the box": used where a
   search direction has a vanishing component along an axis, so that the
   axis' box constraint must not enter the MIN that bounds the step        */
#define SS_FREE 1.0e+30

typedef double vofi_real;
typedef const double vofi_creal;
typedef const int vofi_cint;
typedef int  vofi_int;
typedef int * const vofi_int_cpt;
typedef void * const vofi_void_cptr;
typedef double (*integrand) (vofi_creal [],vofi_void_cptr);

/* min_data structure: 
   xval: coordinates of the minimum or where the function sign has changed, 
   fval: function value at xval, 
   sval: distance from the starting point, applicable only on a side,
   isc:  in 2D (on a side, at the end of direction-ordering function) 
         [0]) (1/0) --> (yes/no) change of sign, 
         [1]) (-1/1/0) --> (1/2/0) changes of sign on the lower secon. side, 
         [2]) (-1/1/0) --> (1/2/0) changes of sign on the upper secon. side,

         in 3D (on a tertiary side) 
         [0]) (1/0) --> (yes/no) change of sign, 
         [1]) (-1/1/0) --> (1/2/0) changes of sign on the side, 
         [2]) not used

         in 3D (on a face, anywhere) 
         [0]) (1/0) --> (yes/no) change of sign, 
         [1]) (0/1) without/with sign change on the lower sec./ter. face, 
         [2]) (0/1) without/with sign change on the upper sec./ter. face, 
   ipt:  tentative number of integration points                               */
typedef struct {
  vofi_real xval[NDIM];
  vofi_real fval;
  vofi_real sval;
  vofi_int isc[NDIM+1];
  vofi_int ipt;
} min_data;

/* dir_data structure:
   ind1,ind2,(ind3): indices in {0,1} to locate the vertex in a face (xy,yz,xz)
                     or, in 4D, in one of the eight cubic 3-faces,
   swt1,swt2,(swt3): switches to turn on/off the components of the initial
                     gradient,
   consi: if = 0 no sign change is possible, otherwise sign to have f>0       */
typedef struct {
  int ind1; int ind2; int ind3;
  int swt1; int swt2; int swt3;
  int consi;
} dir_data;

/* len_data structure:
   np0: actual number of internal nodes,
   f_sign: (+1/-1) --> local height stems from (lower/upper) boundary, 
   NGLM +2: maximum number of internal nodes + boundary nodes,
   xt0: nodes coordinate along the secondary direction (1D),
   ht0: local height at xt0,
   htp: local derivative along the primary direction (1D)                     */
typedef struct {
  vofi_int np0;
  vofi_int f_sign;
  vofi_real xt0[NGLM+2];
  vofi_real ht0[NGLM+2];
  vofi_real htp[NGLM+2];
} len_data; 

/* meas_acc: the 4D interface measure, accumulated where the heights are
   actually computed. With one frame for the whole cell the interface is the
   graph x_p = H(s,t,u), so its 3-volume is the integral of
   sqrt(1 + |grad H|^2) = |grad f| / |df/dp| over the same (s,t,u) nodes the
   hypervolume quadrature already visits -- no second traversal, no separate
   surface reconstruction. wout carries the product of the outer (t and u)
   quadrature weights; tloc and uloc the current coordinates in the local
   frame. NULL in 1D, 2D and 3D, which keep their geometric measures.     */
typedef struct {
  integrand  func;
  void      *par;             /* not vofi_void_cptr: this member is assigned */
  vofi_int   jp;              /* global axis carrying pdir                */
  vofi_real  wout;            /* outer quadrature weight                  */
  vofi_real  tloc, uloc;      /* current t and u in the local frame       */
  vofi_real  dh[NDIM];        /* central-difference step per global axis  */
  vofi_real  meas;            /* accumulated 3-volume of the interface    */
  vofi_real  mom[NDIM];       /* accumulated first moment, local frame    */
} meas_acc;

/*------------ function prototypes ------------*/

/* function to compute the root along an oriented segment */
vofi_real vofi_get_segment_zero(integrand,vofi_void_cptr,vofi_creal [],
                                vofi_creal [],vofi_real [],vofi_cint);

/* functions to check consistency with a minimum along a cell side, */
/* a line, a cell face or and edge intersection                     */
vofi_int vofi_check_side_consistency(integrand,vofi_void_cptr,vofi_creal [],
                                     vofi_creal [],vofi_creal [],vofi_creal);
vofi_int vofi_check_line_consistency(integrand,vofi_void_cptr,vofi_creal [],
                                     vofi_creal [],vofi_creal ,vofi_cint,
                                     min_data *);
dir_data vofi_check_face_consistency(integrand,vofi_void_cptr,vofi_creal [],
                                     vofi_creal [],vofi_creal [],
                                     vofi_creal [],vofi_creal []);
void vofi_check_edge_consistency(integrand,vofi_void_cptr,vofi_real [],
				 vofi_creal [],vofi_real [],vofi_creal [],
                                 vofi_creal,vofi_cint);

/* functions to check if there is a double intersection along a */
/* cell side, a segment or a "cap" intersection on a cell face  */
vofi_int vofi_check_boundary_line(integrand,vofi_void_cptr,vofi_creal [],
                                  vofi_creal [],vofi_real [][NSE],
                                  min_data *,vofi_int [][NSE]);
void vofi_check_secondary_side(integrand,vofi_void_cptr,vofi_creal [],
                               vofi_creal [],vofi_creal [],vofi_creal [],
                               vofi_real [][NSE],min_data *,vofi_creal);
vofi_int vofi_check_boundary_surface(integrand,vofi_void_cptr,vofi_creal [],
                                     vofi_creal [],vofi_real [][NSE][NSE],
                                     min_data[],vofi_int [][NSE][NSE]);
void vofi_check_secter_face(integrand,vofi_void_cptr,vofi_creal [],
                            vofi_creal [],vofi_creal [],vofi_creal [],
                            vofi_creal [],vofi_real [][NSE][NSE],min_data *,
                            vofi_creal);
void vofi_check_tertiary_side(integrand,vofi_void_cptr,vofi_creal [],
                              vofi_creal [],vofi_creal [],vofi_creal [],
                              vofi_creal [],vofi_real [][NSE][NSE],min_data[],
                              vofi_creal);

/* function to compute the arclength of the interface in a cell (2D),  */
/* and, if scent != NULL, the length-weighted first moment of the same */
/* polyline in the local (pdir,sdir) frame                             */
double vofi_interface_length(integrand,vofi_void_cptr,vofi_creal [],
                             vofi_creal [],vofi_creal [],vofi_creal [],
                             len_data [],vofi_real [],vofi_cint);

/* functions to compute the interface intersections with a cell side */
/* or the external limits of a cap-like intersections                */
vofi_int vofi_get_side_intersections(integrand,vofi_void_cptr,vofi_real [],
                                     vofi_creal [],min_data,vofi_real [],
                                     vofi_creal [],vofi_creal,vofi_int,
                                     vofi_cint);
vofi_int vofi_get_ext_intersections(integrand,vofi_void_cptr,vofi_creal [],
                                    vofi_creal [],min_data,vofi_real [],
                                    vofi_creal [],vofi_creal [],vofi_cint);

/* functions to compute the limits of integration along the */
/* secondary and tertiary directions                        */
vofi_int vofi_get_limits_2D(integrand,vofi_void_cptr,vofi_creal [],vofi_creal [],
                            vofi_real [][NSE],min_data,vofi_real [],
                            vofi_creal [],vofi_creal [],vofi_int[],vofi_int[]);
vofi_int vofi_get_limits_3D(integrand,vofi_void_cptr,vofi_creal [],vofi_creal [],
                            vofi_real [][NSE][NSE],min_data [],vofi_real [],
                            vofi_creal [],vofi_creal [],vofi_creal []);
vofi_int vofi_check_plane(integrand,vofi_void_cptr,vofi_creal [],vofi_creal [],
                          min_data *,vofi_real [],vofi_creal [],vofi_creal [],
                          vofi_int [],vofi_int []);
vofi_int vofi_get_limits_inner_2D(integrand,vofi_void_cptr,vofi_creal [],
                                  vofi_creal [],min_data *,vofi_real [],
                                  vofi_creal [],vofi_creal [],vofi_int [],
                                  vofi_int [],vofi_cint);
vofi_int vofi_get_limits_edge_2D(integrand,vofi_void_cptr,vofi_creal [],
                                 vofi_creal [],min_data *,vofi_real [],
                                 vofi_creal [],vofi_creal []);
void vofi_reorder(vofi_real [],vofi_int [],vofi_int);
vofi_int vofi_rm_segs(vofi_real [],vofi_int [],vofi_int);
void vofi_sector_new(vofi_int [][NSCT],vofi_int [],vofi_int [],vofi_cint,
                     vofi_cint,vofi_cint,vofi_cint);
void vofi_sector_old(integrand,vofi_void_cptr,vofi_creal [],vofi_creal [],
                     vofi_creal [],vofi_creal [],vofi_creal [],vofi_int [],
                     vofi_int [],vofi_cint);

/* threshold on |f| at the 16 vertices of a 4D cell, and the gradient it is
   built from (see vofi_hyper_fth in orderdirs4d.c)                        */
vofi_real vofi_hyper_fth(vofi_creal [],vofi_creal [],vofi_real []);

/* Size of the min_data table handed around by the direction-ordering
   routines. 3D uses [0..3] for the four tertiary sides and [4] for the
   two faces normal to pdir. 4D uses [0..7] for the eight cell edges
   parallel to udir -- edge le carries the {0,1} indices of the three
   non-u axes in the bits of le, least significant bit = lowest axis --
   and [IXC4] for the tentative point count plus the seed point left by
   the boundary check.                                                     */
#define NXFS   9
#define IXC4   8

/* consistency check and minimum search inside a cubic 3-face of a 4D cell:
   the 3-direction extension of vofi_check_face_consistency/vofi_get_face_min */
dir_data vofi_check_cell_consistency(integrand,vofi_void_cptr,vofi_creal [],
                                     vofi_creal [],vofi_creal [],vofi_creal [],
                                     vofi_creal [],vofi_creal []);
vofi_int vofi_get_cell_min(integrand,vofi_void_cptr,vofi_creal [],vofi_creal [],
                           vofi_creal [],vofi_creal [],vofi_creal [],
                           vofi_creal [],min_data *,dir_data);

/* look for an interface that lives entirely inside the boundary of the 4D
   cell, i.e. inside one of its eight cubic 3-faces                       */
vofi_int vofi_check_boundary_hyper(integrand,vofi_void_cptr,vofi_creal [],
                                   vofi_creal [],vofi_creal [],min_data *,
                                   vofi_int []);

/* cell type, ordered directions, external limits and quadrature (4D) */
vofi_int vofi_cell_type_4D(integrand,vofi_void_cptr,vofi_creal [],vofi_creal []);
vofi_int vofi_order_dirs_4D(integrand,vofi_void_cptr,vofi_creal [],vofi_creal [],
                            vofi_real [],vofi_real [],vofi_real [],vofi_real [],
                            vofi_real [],min_data []);
vofi_int vofi_get_limits_4D(integrand,vofi_void_cptr,vofi_creal [],vofi_creal [],
                            vofi_creal [],min_data [],vofi_real [],
                            vofi_creal [],vofi_creal [],vofi_creal [],
                            vofi_creal []);
vofi_int vofi_check_hyperplane(integrand,vofi_void_cptr,vofi_creal [],
                               vofi_creal [],vofi_creal [],vofi_creal [],
                               vofi_creal [],vofi_real [],min_data [],
                               vofi_int *);
vofi_int vofi_get_uext_intersections(integrand,vofi_void_cptr,vofi_creal [],
                               vofi_creal [],min_data,vofi_real [],
                               vofi_creal [],vofi_creal [],vofi_creal [],
                               vofi_cint);
vofi_real vofi_get_hypervolume(integrand,vofi_void_cptr,vofi_creal [],
                               vofi_creal [],vofi_creal [],vofi_creal [],
                               vofi_creal [],vofi_creal [],vofi_creal [],
                               vofi_real [],vofi_cint [],vofi_cint [],
                               vofi_cint,vofi_cint);
void vofi_meas_add(meas_acc *,vofi_creal [],vofi_creal [],vofi_creal,
                   vofi_creal,vofi_creal,vofi_creal,vofi_cint);

/* the direction-generic form of vofi_check_boundary_surface, and the 4D
   counterparts of vofi_check_secter_face / vofi_check_tertiary_side      */
vofi_int vofi_check_boundary_cell(integrand,vofi_void_cptr,vofi_creal [],
                               vofi_creal [],vofi_real [][NSE][NSE],min_data [],
                               vofi_int [][NSE][NSE],vofi_creal [],
                               vofi_creal [],vofi_creal []);
void vofi_check_secter_cell(integrand,vofi_void_cptr,vofi_creal [],
                               vofi_creal [],vofi_creal [],vofi_creal [],
                               vofi_creal [],vofi_creal [],vofi_creal [],
                               min_data *,vofi_creal);
void vofi_check_quaternary_side(integrand,vofi_void_cptr,vofi_creal [],
                               vofi_creal [],vofi_creal [],vofi_creal [],
                               vofi_creal [],vofi_creal [],vofi_creal [],
                               min_data [],vofi_creal);

/* functions to compute the function minimum along a segment or */
/* in a cell face                                               */
vofi_int vofi_get_segment_min(integrand,vofi_void_cptr,vofi_creal [],
                              vofi_creal [],vofi_creal [],min_data *,
                              vofi_creal,vofi_cint);
vofi_int vofi_get_face_min(integrand,vofi_void_cptr,vofi_creal [],
                           vofi_creal [],vofi_creal [],vofi_creal [],
                           vofi_creal [],min_data *,dir_data);

/* functions to compute the cell type */
vofi_int vofi_cell_type_2D(integrand,vofi_void_cptr,vofi_creal [],
                           vofi_creal []);
vofi_int vofi_cell_type_3D(integrand,vofi_void_cptr,vofi_creal [],
                           vofi_creal []);

/* functions to integrate in 2D (1D Gauss-Legendre integration) */
/* and in 3D (2D Gauss-Legendre integration)                    */
vofi_real vofi_get_area(integrand,vofi_void_cptr,vofi_creal [],vofi_creal [],
                        vofi_creal [],vofi_creal [],vofi_creal [],
                        len_data [],vofi_real [],vofi_cint,vofi_cint [],
                        vofi_cint,vofi_cint,vofi_int [],vofi_int [],
                        meas_acc *);
double vofi_get_volume(integrand,vofi_void_cptr,vofi_creal [],vofi_creal [],
                       vofi_creal [],vofi_creal [],vofi_creal [],
                       vofi_creal [],vofi_real [],vofi_cint [],vofi_cint [],
                       vofi_cint,vofi_cint,vofi_cint [],meas_acc *);

/* functions to compute the cell type and, if cut, to order coordinate */ 
/* directions and to compute tentative number of integration points    */
vofi_int vofi_order_dirs_2D(integrand,vofi_void_cptr,vofi_creal [],
                            vofi_creal [],vofi_real [],vofi_real [],
                            vofi_real [][NSE],min_data *);
vofi_int vofi_order_dirs_3D(integrand,vofi_void_cptr,vofi_creal [],
                            vofi_creal [],vofi_real [],vofi_real [],
                            vofi_real [],vofi_real [][NSE][NSE],
                            min_data []);
void vofi_xyz2pst(vofi_real [][NSE][NSE],vofi_cint,vofi_cint,vofi_cint);

/* functions to print data in Tecplot format (ASCII data format: *.dat) */
void tecplot_heights(vofi_creal [],vofi_creal[],vofi_creal [],
                     vofi_creal [],len_data []);
void tecplot_arcline(vofi_creal [],vofi_creal [],vofi_creal [],
                     vofi_creal,vofi_creal,vofi_creal,vofi_cint,FILE *);
void tecplot_triangle(vofi_creal [],vofi_creal [],vofi_creal [],
	              vofi_creal [],vofi_creal [],vofi_creal [],
		      vofi_creal [],vofi_creal, vofi_cint);

/* functions to triangulate the interface by adding end points along     */
/* the secondary direction on integration planes, and all points on edge */
/* planes, and to compute the triangles area */
/* public entry points (see vofi.h): vofi_get_cc is the historical one,   */
/* vofi_get_cc_gam adds the interface centroid                            */
double vofi_get_cc(integrand,vofi_void_cptr,vofi_creal [],vofi_creal [],
                   double [],vofi_cint [],vofi_cint [],vofi_cint [],vofi_cint);
double vofi_get_cc_gam(integrand,vofi_void_cptr,vofi_creal [],vofi_creal [],
                       double [],double [],vofi_cint [],vofi_cint [],
                       vofi_cint [],vofi_cint);
vofi_int vofi_get_cell_type(integrand,vofi_void_cptr,vofi_creal [],
                            vofi_creal [],vofi_cint);

/* the last argument accumulates the area-weighted interface centroid in  */
/* the (pdir,sdir,tdir) frame; pass NULL not to compute it               */
double vofi_interface_surface(integrand,vofi_void_cptr,vofi_creal [],
                              vofi_creal [],vofi_creal [],vofi_creal [],
                              vofi_creal [],vofi_creal [],len_data [],
                              len_data [],vofi_cint,vofi_cint,vofi_cint,
                              double []);
void vofi_end_points(integrand,vofi_void_cptr,vofi_creal [],vofi_creal [],
                     vofi_creal [],vofi_creal [],len_data []);
void vofi_edge_points(integrand,vofi_void_cptr,vofi_creal [],
                      vofi_creal [],vofi_creal [],vofi_creal [],
                      vofi_creal [],len_data [],vofi_cint [],
                      vofi_cint ,vofi_int [],vofi_int []);
vofi_real vofi_triarea(vofi_creal [],vofi_creal [],vofi_creal []);

#endif
