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
/* FILE : parameters.h                                                    */
/*                                                                        */
/* Author :       Eric C. Rouchka                                         */
/* Last Update :  April 29, 1996                                          */
/**************************************************************************/

#ifndef PARAMETERS_H
#define PARAMETERS_H
#include "common.h"
#ifdef _GUI_
#include "GUI_common.h"
#endif
#include "fragadj.h" 
#include "bayes.h" 
#include <unistd.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <unistd.h>

#define POST_PLATEAU_SAMPLES  100

short stripargs(int argc, char **argv, Model B, short bGui );	  /* BT 4/29/98 */
/* void stripargs(int argc, char **argv, IPtype IP); */
short    ParseReals(char *str, double *values);             /* BT 7/21/97 */
int      ParseIntegers(char *str, int *values, char *msg);  /* BT 7/21/97 */

void confirm_data(IPtype IP);
void set_filename(IPtype IP);
void set_sequence_type(IPtype IP);
void set_numMotifs(IPtype IP);
void get_inputs(Model B);
void set_motifLen(IPtype IP);
void set_numSeeds(IPtype IP);
void set_plateauPer(IPtype IP);
void set_pseudoWt(IPtype IP);
void set_MaxIterations(IPtype IP);
void set_SeedVal(IPtype IP);
void set_outputFile(IPtype IP);
void set_pseudoSiteWt(IPtype IP);
void set_priorFile(IPtype IP);

void print_options( Model B );			/* BT 1/27/97 */
void printargs(int argc, char **argv, IPtype IP, FILE *tmpFpt );   /* BT 1/9/98 */

void PrintTempOut(FILE *fpt, char *fmt, ...);
void PrintSeqDescriptions( Model M );

#endif

