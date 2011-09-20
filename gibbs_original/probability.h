/****************************************************************/
/* Gibbs - A program for detecting subtle sequence signals      */
/*                                                              */
/* Please acknowledge the program authors on any publication of */
/* scientific results based in part on use of the program and   */
/* cite the following articles in which the program was         */
/* described.                                                   */
/* For data involving protein sequences,                        */
/* Detecting subtle sequence signals: A Gibbs sampling          */
/* strategy for multiple alignment. Lawrence, C. Altschul,      */
/* S. Boguski, M. Liu, J. Neuwald, A. and Wootton, J.           */
/* Science, 262:208-214, 1993.                                  */
/* and                                                          */
/* Bayesian models for multiple local sequence alignment        */
/* and Gibbs sampling strategies, Liu, JS. Neuwald, AF. and     */
/* Lawrence, CE. J. Amer Stat. Assoc. 90:1156-1170, 1995.       */
/* For data involving nucleotide sequences,                     */
/* Gibbs Recursive Sampler: finding transcription factor        */
/* binding sites. W. Thompson, E. C. Rouchka and                */
/* C. E. Lawrence, Nucleic Acids Research, 2003,                */
/* Vol. 31, No. 13 3580-3585.                                   */
/*                                                              */
/* Copyright (C) 2006   Health Research Inc.                    */
/* HEALTH RESEARCH INCORPORATED (HRI),                          */
/* ONE UNIVERSITY PLACE, RENSSELAER, NY 12144-3455.             */
/* Email:  gibbsamp@wadsworth.org                               */
/*                                                              */
/****************************************************************/
/*                                                              */
/* Changes Copyright (C) 2007   Brown University                */
/* Brown University                                             */
/* Providence, RI 02912                                         */
/* Email:  gibbs@brown.edu                                      */
/*                                                              */
/* For the Centroid sampler, please site,                       */
/* Thompson, W.A., Newberg, L., Conlan, S.P., McCue, L.A. and   */
/* Lawrence, C.E. (2007) The Gibbs Centroid Sampler             */
/* Nucl. Acids Res., doi:10.1093/nar/gkm265                     */
/*                                                              */
/* For the Phylogenetic Gibbs Sampler, please site,             */
/* Newberg, L., Thompson, W.A., Conlan, S.P., Smith, T.M.,      */
/* McCue, L.A. and Lawrence, C.E. (2007) A phylogenetic Gibbs   */
/* sampler that yields centroid solutions for cis regulatory    */
/* site prediction., Bioinformatics,                            */
/* doi:10.1093/bioinformatics/btm241.                           */
/*                                                              */
/****************************************************************/
/*                                                              */
/* This program is free software; you can redistribute it       */
/* and/or modify it under the terms of the GNU General Public   */
/* License as published by the Free Software Foundation;        */
/* either version 2 of the License, or (at your option)         */
/* any later version.                                           */
/*                                                              */
/* This program is distributed in the hope that it will be      */
/* useful, but WITHOUT ANY WARRANTY; without even the implied   */
/* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR      */
/* PURPOSE. See the GNU General Public License for more         */
/* details.                                                     */
/*                                                              */
/* You should have received a copy of the GNU General Public    */
/* License along with this program; if not, write to the        */
/* Free Software Foundation, Inc., 51 Franklin Street,          */
/* Fifth Floor, Boston, MA  02110-1301, USA.                    */
/****************************************************************/
/**************************************************************************/
/* FILE : probability.h                                                   */
/*                                                                        */
/* Author :       Eric C. Rouchka                                         */
/* Last Update :  May 29, 1996                                            */
/**************************************************************************/

#ifndef PROB_H
#define PROB_H

#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "stdio.h"
#include "math.h"
#include "sampling.h"
#include "width.h"
#include "recur.h"
#include "motif.h"
#include "counts.h"
#include "fels.h"
#include "impsamp.h"
#include "weights.h"
#include "freqbkgnd.h"
#include "species_sample.h"

#define RCH2INT(n,R)	 ((int)(R)[(n)] - 97)

double ln_gamma(double a);
void update_prob(ProbStruct *P, Model B, short bUpdateBkgrnd);
void   init_prob(ProbStruct *P, IPtype IP);
double in_motif_prob(Model B, register ProbStruct P, int i, int t, 
		     short rev_comp, short bSamplingRatio, 
		     double *dMotifProb, double *dBGProb);
long double LnBeta(double a, double b);
double      FindMapProb(Model B, int t);
void set_posterior_prob(IPtype IP, Counts *C);
void update_posterior_prob(Model B);
void CheckMotifProb(Model B, int n, int t, PoSition **Pos, 
                    ProbStruct *P, Mlist M, int iterations, double ***dFinalProb);
void free_prob(ProbStruct *P, IPtype IP);
double CalcMapProb( Model B, short bCalcPalin );
double CalcMotifMap( Model B, int t, short bCalcPalin );
double CalcBetaMap( Model B, short bCalcPalin );
double CalcBkgndMap( Model B, short bCalcPalin );
double CalcMotifFragMap( Model B, int t, short bCalcPalin );
double CalcNullMap( Model B );
double CalcSitePerSeqMap( Model B );
double CalcBetaPartialSum(Model B, int nSeq, int t, int nPrevSite, int nPrevType, int curPos, int nSites );
int PossibleBkgndPos( Model B, int pos, int offset );

#endif
