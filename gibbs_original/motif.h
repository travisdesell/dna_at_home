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
/***************************************************************************/
/* FILE : motif.h                                                          */
/*                                                                         */
/* Author :      Eric C. Rouchka                                           */
/* Last Update : July  3, 1996                                             */ 
/***************************************************************************/

#ifndef MOTIF_H
#define MOTIF_H


#include "common.h"
#include "printdata.h"
#include "recur.h"


/*-------------------------------------------------------------------------*/
/*                              FUNCTION PROTOTYPES                        */
/*-------------------------------------------------------------------------*/
int **copy_motif_num (IPtype IP);
void  zero_motifs    (IPtype IP);
void  reset_motif_num(int **nNumMotifs, IPtype IP); /* BT 4/4/97 */ 
short PossibleOverlapsPrevMotif( Model B, int curPos, int nMotifType,  PoSition **Pos );
short OverlapsPrevMotif(int CurPos,int nNumMotif,PoSition **Pos);
short PossibleMotifPos(int i, PoSition **Pos, IPtype IP, int t, Model B );     /* BT 7/5/97 */
Mlist initializeMotifs(IPtype IP);
void  set_in_motif(PoSition **Pos,int i, Model B,int t,short RevComp);
void  not_in_motif(PoSition **Pos, int i,Model B, int t);
int   overlap(int i, PoSition **Pos,int *nMotifLen,int *newpos,int t, int t1);
Mlist set_motif_info(IPtype IP, int **StartPos, Stringstruct *Seq);

void  add_motif(IPtype IP, Stype Seq, int pos, Mlist M, int t, short RevComp);
void  delete_motif(Model B, int pos, Mlist M, int t);		/* BT 3/31/97 */
void  print_motifs(Model B, Mlist M, int t);
void  free_motifs(Model B, Mlist M);
int AnyOverlap(Model B, int n, int t, PoSition **Pos, int *newpos, int *typ);

#endif

