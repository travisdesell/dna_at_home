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
/* FILE : sampling.h                                                      */
/*                                                                        */
/* Author :       Eric C. Rouchka                                         */
/* Last Update :  July  8, 1996                                           */
/**************************************************************************/

#ifndef SAMPLING_H
#define SAMPLING_H

/*-------------------------------------------------------------------------*/
/*----------------------------- INCLUDE FILES -----------------------------*/
/*-------------------------------------------------------------------------*/

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/times.h>
#include <limits.h>

#include "column.h"
#include "common.h"
#include "counts.h"
#include "expectmax.h"
#include "fragadj.h"
#include "fragment.h"
#include "initseq.h"
#include "maximize.h"
#include "motif.h"
#include "motifsamp.h"
#include "nearopt.h"
#include "parameters.h"
#include "printdata.h"
#include "probability.h"
#include "sitesamp.h"
#include "subopt.h"
#include "mem_mgmt.h"
#ifdef _GUI_
#include "GUI_common.h"
#endif
#include "width.h"
#include "recur.h"
#include "motifpair.h"
#include "rsample.h"
#include "weights.h"
#include "centroid.h"


/*-------------------------------------------------------------------------*/
/*-------------------------- FUNCTION PROTOTYPES --------------------------*/
/*-------------------------------------------------------------------------*/

void       find_best_sites(register Model B);
MaxResults motif_sampler(Model B, PoSition **Pos,Mlist M);
MaxResults setMaxData(Ftype F, IPtype IP, double dProb, int iteration,
                      PoSition **Pos, int *last_increase, double *dMaxProb,
                      double *dLocMax, MaxResults maxData, double **dProbArray,
		      Model B);
void       reset_values(int *iterations, int *last_increase, double *dLocMax,
                  ProbStruct *P, IPtype IP);
double **  setElementProb(Model B, PoSition **Pos, ProbStruct P);
void       CalcPosProb( Model B, ProbStruct P, PoSition **Pos, int n, double **posProb );
void       setMotifNum(IPtype IP, MaxResults maxData);
void       init_values(double ***dFinalProb, PoSition **Pos, int n, int t);
void       check_map(Model B, double *dCurrProb, int iterations,PoSition **Pos,
                     int *last_increase, double *dMaxProbability, 
                     double *dLocMax, MaxResults *maxData, ProbStruct P);
void       PrintDoubleAsBytes( double *d );
void       print_Pos(PoSition **Pos, Model B, Mlist M);


/*-------------------------------------------------------------------------*/

#endif

