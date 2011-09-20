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
#ifndef RECUR_H
#define RECUR_H

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "common.h"
#include "fragadj.h"
#include "counts.h"
#include "probability.h"
#include "rsample.h"
#include "width.h"
#include "prior.h"

#define GRAPH_SIZE 8192

int SequenceStartPos( Model B, int nSeq );
int SequenceEndPos( Model B, int nSeq );
int SequenceLength( Model B, int nSeq );
int SequenceFromPosition( Model B, int nPos );
int MinLenMotif( Model B );
int MaxLenMotif( Model B );
int MaxSequenceLen( Model B );
int TotalNumMotifs( Model B );
int MotifWidth( Model B, int t );
int SpeciesInc( Model B, int seq );
int SeqInc( Model B, int seq );
int OffsetInMASS( Model B, int seq );

int IsPhyloSeq( Model B, int seq );
int PosInFirstMASSSeq( Model B, int seqPos );
int IsPhyloSeqBySpecies( Model B, int seq );
int NumWidths( Model B );
int MaxPossibleWidth( Model B );
int MinPossibleWidth( Model B );
void InitRProbStruct( Model B );
void InitRProbArray( Model B, int nSeq );
void FreeRProb( Model B );
void FreeRProbArray( Model B );
void CalcRProbArray( Model B,  PoSition **Pos, int nSeq, ProbStruct *P, int iter );
void DumpProbArray( Model B, int nSeq, int kmax, int iter, int seed );
void CountAlignments( Model B, PoSition **Pos );
void FreeAlignmentCounts( Model B );
double BiCoef( int n, int k);
double LnBiCoef( int n, int k);
void CalcAlignProbFromFootprint( Model B );
double BkgndProbFromComp( Model B, int nSeq, int nPos, int t, short bRevComp );
void SaveAlignmentCounts( Model B ); 
void RestoreAlignmentCounts( Model B );   
short PossMotifStartPos(PoSition **Pos, int n, int t, Model B); 
short PossibleStartPos(PoSition **Pos, int n, int t, int nOffset, int nSeq, Model B);
char *GibbsTempnam(const char *dir, const char *pfx);
#endif

