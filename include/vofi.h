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

vofi_real vofi_get_cc(integrand,vofi_void_cptr,vofi_creal [],
                      vofi_creal [],vofi_real [],vofi_cint [],
                      vofi_cint [],vofi_cint [],vofi_cint);

 vofi_int vofi_get_cell_type(integrand,vofi_void_cptr,vofi_creal [],
			     vofi_creal [],vofi_cint);

#ifdef __cplusplus
}
#endif

#endif
