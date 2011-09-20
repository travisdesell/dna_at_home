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
#include "fels.h"

#define MSG_LEN 256

void CalcFels2( Model B, PhyloTree tree, int seqPos, int rev_comp, int seqLen );

/* On return the current substitution matrix will be the one calculated for the motif at this position */
/* CalcFelsBySeq does not have to update it                                                            */

void CalcFels( Model B, ProbStruct P, int seqPos, int motifPos, int t, int rev_comp, 
	       double *motifProb, double *bkgndProb )
{
  IPtype     IP;
  PhyloType  PH;
  int        seq;
  int        seqLen;
  int        seqStart;
  PhyloTree  tree;
  double     *prob;
  double     *bgProb;
  char       *msg;

  IP = B->IP;
  PH = B->Phylo;

  seq = SequenceFromPosition( B, seqPos );
  tree = PH->phyloTree[PH->phyloIndex[seq]];

  if( ! IsPhyloSeq( B, seq ) )
    {
      NEW( msg, MSG_LEN, char );      
      sprintf( msg, "CalcFels: sequence %d is not a phylogentic sequence", seq + 1);
      p_internal_error( msg );
    }

  if( seq != PH->phyloTreeStartSeq[PH->phyloIndex[seq]] )
    {
      NEW( msg, MSG_LEN, char );      
      sprintf( msg, "CalcFels: sequence %d is not the start of a MASS", seq + 1);
      p_internal_error( msg );
    }

  seqLen = SequenceLength( B, seq );
  seqStart = SequenceStartPos( B, seq );

  if( B->IP->is_defined[cl_B] )
    bgProb = B->BP->dBkgndProb[seq][seqPos - seqStart];
  else
    bgProb = P.dvInBGProb[t];
  CalcSubsMatrix( B, tree, bgProb );
  *bkgndProb = CalcFelsProb( B, tree, bgProb, seqPos, FALSE );

  prob = P.dvInMotifProb[t][motifPos];
  CalcSubsMatrix( B, tree, prob );
  *motifProb = CalcFelsProb( B, tree, prob, seqPos, rev_comp );
}


void CalcFels2( Model B, PhyloTree tree, int seqPos, int rev_comp, int seqLen )
{
  int        a;
  int        b;
  int        index;
  double     sum;
  double     sum2;
  PhyloTree  lChild;
  PhyloTree  rChild;
  double     *likelihood;
  double     *subs;

  lChild = tree->leftChild;
  rChild = tree->rightChild;

  if( lChild == NULL && rChild == NULL )
    {
      likelihood = tree->likelihood;

      for( b = 0; b < B->IP->nAlphaLen; b++ ) 
	likelihood[b] = 0.0;
      
      index = (! rev_comp) ? RCH2INT(seqPos + tree->offset*seqLen, *(B)->Seq->R ) :
	                 nComp[RCH2INT(seqPos + tree->offset*seqLen, *(B)->Seq->R )]; 
      if( index < 0 || index > 3 )
	{
	  for( b = 0; b < B->IP->nAlphaLen; b++ ) 
	    likelihood[b] = 0.25;
	}
	/*	p_internal_error( "CalcFels2: invalid character index" ); */
      else
	likelihood[index] = 1.0;
    }
  else
    {
      CalcFels2( B, lChild, seqPos, rev_comp, seqLen );
      CalcFels2( B, rChild, seqPos, rev_comp, seqLen );
      
      for( a = 0; a < B->IP->nAlphaLen; a++ )
	{
	  likelihood = lChild->likelihood;
	  subs = lChild->subs[a];
	  for( sum = 0, b = 0; b < B->IP->nAlphaLen; b++ )     
	    sum += likelihood[b] * subs[b];

	  likelihood = rChild->likelihood;
	  subs = rChild->subs[a];
	  for( sum2 = 0, b = 0; b < B->IP->nAlphaLen; b++ )     
	    sum2 += likelihood[b] * subs[b];
	  
	      /* sum +=  tree->leftChild->likelihood[b] * tree->leftChild->subs[b][a];   */ 
	      /* this is wrong -- check with tree dist far enough to be independent*/
	  
	  tree->likelihood[a] = sum * sum2;
	}
    }
}


double CalcFelsProb( Model B, PhyloTree tree, double *theta, int seqPos, int rev_comp )
{
  IPtype     IP;
  int        seq;
  int        seqLen;
  int        i;
  double     prob = 0.0;

  IP = B->IP;

  seq = SequenceFromPosition( B, seqPos );
  seqLen = SequenceLength( B, seq );
  
  CalcFels2( B, tree, seqPos, rev_comp, seqLen );
      
  for( i = 0; i < IP->nAlphaLen; i++ ) 
    {
      prob += theta[i] * tree->likelihood[i];
    }

  return( prob );  
}


