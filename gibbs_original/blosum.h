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
/* $Header: /home/Bill/cvs/gibbs/src/blosum.h,v 1.4 2007/05/23 18:19:55 Bill Exp $                                                             */
/*                                                                        */
/* $Author: Bill $                                                              */
/* $Date: 2007/05/23 18:19:55 $                                                                */
/* $Description:  Substituion matrix for amino acids                      */
/* The detail about the matrix can be found in  Henikofss &Henikoff       */
/* Amino acid substitution matrices from protein blocks in Proc.Natl.Acad */
/* Sci. USA Vol.89, pp.10915-10919, Nov. 1992                             */
/**************************************************************************/

#ifndef BLOSUM_62
#define BLOSUM_62

/* the order like below is because "initseq.c" uses this order */
static short related[21][21] = { /* BLOSUM matrix for probability of */
                          /* mutation from one amino acid sequence to   */
                          /* another                                    */
/* A  V  C  D  E  F  G  H  I  W  K  L  M  N  Y  P  Q  R  S  T */
  {4, 0, 0,-2,-1,-2, 0,-2,-1,-3,-1,-1,-1,-2,-2,-1,-1,-1, 1, 0,-9},
  {0, 4,-1,-3,-2,-1,-3,-3, 3,-3,-2, 1, 1,-3,-1,-2,-2,-3,-2, 0,-9},
  {0,-1, 9,-3,-4,-2,-3,-3,-1,-2,-3,-1,-1,-3,-2,-3,-3,-3,-1,-1,-9},
 {-2,-3,-3, 6, 2,-3,-1,-1,-3,-4,-1,-4,-3, 1,-3,-1, 0,-2, 0,-1,-9},
 {-1,-2,-4, 2, 5,-3,-2, 0,-3,-3, 1,-3,-2, 0,-2,-1, 2, 0, 0,-1,-9},
 {-2,-1,-2,-3,-3, 6,-3,-1, 0, 1,-3, 0, 0,-3, 3,-4,-3,-3,-2,-2,-9},
 { 0,-3,-3,-1,-2,-3, 6,-2,-4,-2,-2,-4,-3, 0,-3,-2,-2,-2, 0,-2,-9},
 {-2,-3,-3,-1, 0,-1,-2, 8,-3,-2,-1,-3,-2, 1, 2,-2, 0, 0,-1,-2,-9},
 {-1, 3,-1,-3,-3, 0,-4,-3, 4,-3,-3, 2, 1,-3,-1,-3,-3,-3,-2,-1,-9},
 {-3,-3,-2,-4,-3, 1,-2,-2,-3,11,-3,-2,-1,-4, 2,-4,-2,-3,-3,-2,-9},
 {-1,-2,-3,-1, 1,-3,-2,-1,-3,-3, 5,-2,-1, 0,-2,-1, 1, 2, 0,-1,-9},
 {-1, 1,-1,-4,-3, 0,-4,-3, 2,-2,-2, 4, 2,-3,-1,-3,-2,-2,-2,-1,-9},
 {-1, 1,-1,-3,-2, 0,-3,-2, 1,-1,-1, 2, 5,-2,-1,-2, 0,-1,-1,-1,-9},
 {-2,-3,-3, 1, 0,-3, 0, 1,-3,-4, 0,-3,-2, 6,-2,-2, 0, 0, 1, 0,-9},
 {-2,-1,-2,-3,-2, 3,-3, 2,-1, 2,-2,-1,-1,-2, 7,-3,-1,-2,-2,-2,-9},
 {-1,-2,-3,-1,-1,-4,-2,-2,-3,-4,-1,-3,-2,-2,-3, 7,-1,-2,-1,-1,-9},
 {-1,-2,-3, 0, 2,-3,-2, 0,-3,-2, 1,-2, 0, 0,-1,-1, 5, 1, 0,-1,-9},
 {-1,-3,-3,-2, 0,-3,-2, 0,-3,-3, 2,-2,-1, 0,-2,-2, 1, 5,-1,-1,-9},
 { 1,-2,-1, 0, 0,-2, 0,-1,-2,-3, 0,-2,-1, 1,-2,-1, 0,-1, 4, 1,-9},
 { 0, 0,-1,-1,-1,-2,-2,-2,-1,-2,-1,-1,-1, 0,-2,-1,-1,-1, 1, 5,-9},
 {-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9}};


#endif
