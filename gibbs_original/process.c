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
/***********************************************************************/
/* $Id: process.c,v 1.5 2007/05/23 18:19:57 Bill Exp $        */
/*                                                                     */
/* Author :       Eric C. Rouchka, April 29, 1996                      */
/*                Jun Zhu, October 20, 1996                            */
/*                                                                     */
/* Description :  This file contains functions that are used to        */
/*                preprocess the input strings                         */
/***********************************************************************/

#include "blosum.h"
#include "process.h"

/********************** Xnu_Sequence *********************/
/*                                                       */
/* DESCRIPTION : This function goes through the sequence */
/*               and X's out particular sequences that   */
/*               are closely related using a blosum      */
/*               matrix                                  */
/*=======================================================*/

char *Xnu_Sequence(Model B)
{
  long i, j, low, high, len, v;
  double *pr, lambda, K, H;
  double *freq;
  char *retval;

  if((B->IP->nAlphaLen == 20) && (!B->IP->is_defined[cl_x])) {
     NEW(freq, B->IP->nAlphaLen+1, double);
     for(i = 0; i < B->IP->nAlphaLen; i++) 
        freq[i] = (double)B->First->fCounts[0][BG][i] / 
                  (double)B->First->nTotBack; 
     freq[B->IP->nAlphaLen] = 0.0;
     low = -9;
     high = 11;
     len = high - low + 1;
     NEW(pr, B->IP->nAlphaLen+1, double);
     for(i = 0; i <= B->IP->nAlphaLen; i++)
        pr[i] = 0.0;
     for(i = 0; i <= B->IP->nAlphaLen; i++) {
        for(j = 0; j <= B->IP->nAlphaLen; j++) {
           v = related[i][j] - low;
           pr[v] += freq[i] * freq[j];
        }
     }
     if(!karlin(low, high, pr, &lambda, &K, &H)) {
        fprintf(stderr, "\nusing blast amino acid frequencies\n");
     }
     else
        H = ExpectedInformation(lambda, freq);
     retval = ProcessSeq(B, lambda, K, H);
     free(freq);
     free(pr);

  }
  else retval = (*B->Seq->R);
  return retval;
}

char *ProcessSeq(Model B, double lambda,double K,double H)

   /*=============================================================*/
   /* FUNCTION NAME : ProcessSeq                                  */
   /*                                                             */
   /* DESCRIPTION :                                               */
   /*=============================================================*/

{
        long    i,k,t, off,sum,beg,end,top,noff;
        long    topcut,fallcut;
        char    *hit;
        char    *processed;
        double  s0;
        long ascend=1;
        long descend=1;
        long ncut=4;
        long mcut=0;
        double pcut=0.001;
        long scut=0;

        B->IP->nNumProcessed = 0;
        NEW(processed, B->IP->nSeqLen, char);
        NEW(hit, B->IP->nSeqLen, char);
        for (i=0; i<B->IP->nSeqLen; i++) hit[i]=FALSE;
        noff = B->IP->nSeqLen-1;
        if(ncut>0) noff=ncut;
        if(scut!=0) topcut = scut;
        else {
                s0 = - log( pcut*H / (noff*K) ) / lambda;
                if (s0>0) topcut = floor(s0 + log(s0)/lambda + 0.5);
                else topcut = 0;
        }
        fallcut = (long)log(K/0.001)/lambda;
        for (off=mcut; off<=noff; off++) {
                sum=top=0; beg=off; end=0;
                for(i=off; i<B->IP->nSeqLen; i++) {
                   sum += related[(int)(*B->Seq->R)[i] - 97][(int)(*B->Seq->R)[i-off] - 97];
                        if (sum>top) { top=sum; end=i; }
                        if (top>=topcut && top-sum>fallcut) {
                                for (k=beg; k<=end; k++) {
                                        if (ascend) hit[k] = TRUE;
                                        if (descend) hit[k-off] = TRUE;
                                }
                                sum=top=0; beg=end=i+1;
                        } else if (top-sum>fallcut) {
                                sum=top=0;
                                beg=end=i+1;
                        }
                        if (sum<0) { beg=end=i+1; sum=top=0; }
                }
        }
        for (i=0; i < B->IP->nSeqLen; i++) {
            if(hit[i]) {
               for(t = 0; t < B->IP->nNumMotifTypes; t++)
                  B->First->fCounts[t][BG][(int)(*B->Seq->R)[i] - 97] = 
                      B->First->fCounts[t][BG][(int)(*B->Seq->R)[i]-97]-1;
               processed[i] = 'X';
               B->IP->nNumProcessed++;
            }
            else
               processed[i] = (*B->Seq->R)[i];
        }

        free(hit);
        return(processed);
}
 
