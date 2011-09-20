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
#include "species_sample.h"

#define MSG_LEN 256

double CalcFelsByAlpha( Model B, PhyloTree tree, double *theta, int seqPos, int rev_comp, 
			int seqInMASS, int index );
void CalcFelsByAlpha2( Model B, PhyloTree tree, int seqPos, int rev_comp, int seqLen, 
		       int seqInMASS, int index );


void CalcFelsBySeq( Model B, ProbStruct P, int seqPos, int motifPos, int t, int rev_comp, 
		    double *seqProb, double *bkProb )
{
  IPtype     IP;
  PhyloType  PH;
  PhyloTree  tree;
  double     *prob;
  char       *msg;
  int        seq;
  int        i;
  int        posInFirstMASSSeq;
  int        seqInMASS;
  int        seqStart;
  double     denom = 0;
  double     felsProb;
  double     bgFelsProb;

  IP = B->IP;
  PH = B->Phylo;

  seq = SequenceFromPosition( B, seqPos );
  seqStart = SequenceStartPos( B, seq );
  seqInMASS = OffsetInMASS( B, seq );
  tree = PH->phyloTree[PH->phyloIndex[seq]];

  if( ! IsPhyloSeq( B, seq ) )
    {
      NEW( msg, MSG_LEN, char );      
      sprintf( msg, "CalcFelsBySeq: sequence %d is not a phylogentic sequence", seq + 1);
      p_internal_error( msg );
    }

  posInFirstMASSSeq = PosInFirstMASSSeq( B, seqPos );

  CalcFels( B, P, posInFirstMASSSeq, motifPos, t, rev_comp, &felsProb, &bgFelsProb );

  prob = P.dvInMotifProb[t][motifPos];
  CalcSubsMatrix( B, tree, prob );
  for( i = 0; i < IP->nAlphaLen; i++ )
    denom += CalcFelsByAlpha( B, tree, prob, posInFirstMASSSeq, rev_comp, seqInMASS, i );
  *seqProb = felsProb / denom;

  if( IP->is_defined[cl_B] )
    prob = B->BP->dBkgndProb[seq][seqPos - seqStart];
  else
    prob = P.dvInBGProb[t];
  CalcSubsMatrix( B, tree, prob );
  for( denom = 0, i = 0; i < IP->nAlphaLen; i++ )
    denom += CalcFelsByAlpha( B, tree, prob, posInFirstMASSSeq, rev_comp, seqInMASS, i );
  *bkProb = bgFelsProb / denom;
}


double CalcFelsByAlpha( Model B, PhyloTree tree, double *theta, int seqPos, int rev_comp, 
			int seqInMASS, int index )
{
  IPtype     IP;
  int        seq;
  int        seqLen;
  int        i;
  double     prob = 0.0;

  IP = B->IP;

  seq = SequenceFromPosition( B, seqPos );
  seqLen = SequenceLength( B, seq );
  
  CalcFelsByAlpha2( B, tree, seqPos, rev_comp, seqLen, seqInMASS, index );
      
  for( i = 0; i < IP->nAlphaLen; i++ ) 
    {
      prob += theta[i] * tree->likelihood[i];
    }

  return( prob );  
}



void CalcFelsByAlpha2( Model B, PhyloTree tree, int seqPos, int rev_comp, int seqLen, 
		       int seqInMASS, int index )
{
  int        a;
  int        b;
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

      if( tree->offset != seqInMASS )
	index = (! rev_comp) ? RCH2INT(seqPos + tree->offset*seqLen, *(B)->Seq->R ) :
	                       nComp[RCH2INT(seqPos + tree->offset*seqLen, *(B)->Seq->R )];

      if( index < 0 || index > 3 )
	{
	  for( b = 0; b < B->IP->nAlphaLen; b++ ) 
	    likelihood[b] = 0.25;
	}
      else
	likelihood[index] = 1.0;
    }
  else
    {
      CalcFelsByAlpha2( B, lChild, seqPos, rev_comp, seqLen, seqInMASS, index );
      CalcFelsByAlpha2( B, rChild, seqPos, rev_comp, seqLen, seqInMASS, index );

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

	  tree->likelihood[a] = sum * sum2;
	}
    }
}


double CalcFelsSeqProb( Model B, PhyloTree tree, double *theta, int seqPos, int rev_comp )
{
  IPtype     IP;
  PhyloType  PH;
  char       *msg;
  int        seq;
  int        seqStart;
  int        seqInMASS;
  int        posInFirstMASSSeq;
  int        i;
  double     felsProb;
  double     denom = 0;
  double     seqProb;

  IP = B->IP;
  PH = B->Phylo;

  seq = SequenceFromPosition( B, seqPos );
  seqStart = SequenceStartPos( B, seq );
  seqInMASS = OffsetInMASS( B, seq );
 
  if( ! IsPhyloSeq( B, seq ) )
    {
      NEW( msg, MSG_LEN, char );      
      sprintf( msg, "CalcFelsBySeq: sequence %d is not a phylogentic sequence", seq + 1);
      p_internal_error( msg );
    }

  posInFirstMASSSeq = PosInFirstMASSSeq( B, seqPos );

  felsProb = CalcFelsProb( B,  tree, theta, posInFirstMASSSeq, rev_comp );

  for( i = 0; i < IP->nAlphaLen; i++ )
    denom += CalcFelsByAlpha( B, tree, theta, posInFirstMASSSeq, rev_comp, seqInMASS, i );
  seqProb = felsProb / denom;

  return seqProb;
}
