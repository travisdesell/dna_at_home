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
#ifndef RSITESAMP_H
#define RSITESAMP_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/times.h>
#include <sys/types.h>
#include <limits.h>
#include <time.h>
#include <unistd.h>
/* #include <float.h> */
#include "common.h"
#include "sampling.h"
#include "probability.h"
#include "mem_mgmt.h"
#include "motifsamp.h"
#include "printdata.h"
#include "recur.h"
#include "prior.h"
#include "footprint.h"
#include "nearopt.h"
#include "width.h"
#include "trans.h"
#include "position.h"
#include "seqmap.h"
#include "parameters.h"
#include "exchange.h"
#include "model.h"
#include "weights.h"

#ifdef _MPI_INC_
#include "hier.h"
#endif

#ifdef _CYGWIN_
#define isnand isnan
#endif

#define TICKS_PER_MILLISEC  (((double) sysconf(_SC_CLK_TCK)) / 1000.0)
#define EXCHNG_TIME         500
#define WAIT_PERIOD         5 

#ifndef localtime_r
extern struct tm *localtime_r(const time_t *, struct tm *);
#endif


MaxResults rsite_sampler(Model B, PoSition **Pos, Mlist M, int seed_run );
MaxResults nearopt_rsite_samp(Model B, PoSition **Pos,  Mlist M, int **good, 
                             ProbStruct *P, int ****occurence,
			      double dSubOptProb );
double rsample( Model B, PoSition **Pos, Mlist M, int ****occurence, 
		ProbStruct *P, int iter, int seed_run, short bAccumulateCounts,
		double *dMap, double bestMap );
short PossStartPos( PoSition **Pos, int n, int motif_type, Model B);
void PrintProbModel( Model B, ProbStruct *P, int t );
double ExDistrib( double p, double temp );
void DumpProbModels( Model B, ProbStruct *P, FILE *fpt );

#endif


