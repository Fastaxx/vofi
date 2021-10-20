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
#define NDIM   3
#define NVER   4
#define NSE    2
#define NSEG  10
#define NGLM  20

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
  vofi_int isc[NDIM]; 
  vofi_int ipt;
} min_data;

/* dir_data structure:
   ind1,ind2: indices in {0,1} to locate the vertex in a face (xy,yz,xz), 
   swt1,swt2: switches to turn on/off the components of the initial gradient, 
   consi: if = 0 no sign change is possible, otherwise sign to have f>0       */
typedef struct { 
  int ind1; int ind2;               
  int swt1; int swt2; 
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

/* function to compute the arclength of the interface in a cell (2D) */
double vofi_interface_length(integrand,vofi_void_cptr,vofi_creal [],
                             vofi_creal [],vofi_creal [],vofi_creal [],
                             len_data [],vofi_cint);

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
                                 vofi_creal [],vofi_creal [],vofi_cint);
void vofi_reorder(vofi_real [],vofi_int [],vofi_int);
vofi_int vofi_rm_segs(vofi_real [],vofi_int [],vofi_int);
void vofi_sector_new(vofi_int [][NDIM],vofi_int [],vofi_int [],vofi_cint,
                     vofi_cint,vofi_cint,vofi_cint);
void vofi_sector_old(integrand,vofi_void_cptr,vofi_creal [],vofi_creal [],
                     vofi_creal [],vofi_creal [],vofi_creal [],vofi_int [],
                     vofi_int [],vofi_cint);

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
                        vofi_cint,vofi_cint,vofi_int [],vofi_int []);
double vofi_get_volume(integrand,vofi_void_cptr,vofi_creal [],vofi_creal [],
                       vofi_creal [],vofi_creal [],vofi_creal [],
                       vofi_creal [],vofi_real [],vofi_cint [],vofi_cint [],
                       vofi_cint,vofi_cint,vofi_cint []);

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
double vofi_interface_surface(integrand,vofi_void_cptr,vofi_creal [],
                              vofi_creal [],vofi_creal [],vofi_creal [],
                              vofi_creal [],vofi_creal [],len_data [],
                              len_data [],vofi_cint,vofi_cint,vofi_cint);
void vofi_end_points(integrand,vofi_void_cptr,vofi_creal [],vofi_creal [],
                     vofi_creal [],vofi_creal [],len_data []);
void vofi_edge_points(integrand,vofi_void_cptr,vofi_creal [],
                      vofi_creal [],vofi_creal [],vofi_creal [],
                      vofi_creal [],len_data [],vofi_cint [],
                      vofi_cint ,vofi_int [],vofi_int []);
vofi_real vofi_triarea(vofi_creal [],vofi_creal [],vofi_creal []);

#endif
