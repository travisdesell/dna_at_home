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
#include "recur.h"

void CalcRProbArraySpacing( Model B,  PoSition **Pos, int nSeq, ProbStruct *P, int iter );
double CalcGapSum( Model B, int k, int j, int t, int dir, double dGapPr,
		   int nStart, 
		   double **dPrevSum, double **dPrevNorm,
		   int bCutOff );
double GapProduct( Model B, int nStart, int nEnd, double *dProb, double *dGap, double *fp,
		   double *pSum, double *nSum );
void CalcRProbArraySpacingNoOPt( Model B,  PoSition **Pos, int nSeq, ProbStruct *P, int iter );
void CountAlignments2( Model B,  PoSition **Pos, int nSeq  );
void DumpAlignCnts( Model B, int nSeq, double **algn );

#define DEFAULT_SITE_PROB 1.0;
#define DEFAULT_NONSITE_PROB 1.0;
#define DUMP_CUTOFF 1.0E-06


/*******************************************************************/
/*                          Uitlity routines                       */
/*                                                                 */
/*******************************************************************/


int SequenceStartPos( Model B, int nSeq )
{
  if(nSeq == 0)
    return 0;
  else
    return (*B->Seq->nvEndLocs)[nSeq - 1];
}


int SequenceEndPos( Model B, int nSeq )
{
  return (*B->Seq->nvEndLocs)[nSeq] - 1;
}


int SequenceLength( Model B, int nSeq )
{
  /*  return (SequenceEndPos( B, nSeq ) - SequenceStartPos( B, nSeq ) + 1); */

  return B->Seq->SeqLen[nSeq];  
}


int SequenceFromPosition( Model B, int nPos )
{
  return (B->RP->nSeqPos[nPos]);
}


int MinLenMotif( Model B )
{
  int t;
  int index = -1;
  int minLen = INT_MAX;

  for( t = 0; t < B->IP->nNumMotifTypes; t++ )
    {
      if( MotifWidth( B, t ) < minLen )
	{
	  index = t;
	  minLen =  MotifWidth( B, t );
	}
    }

  return index;
}


int MaxLenMotif( Model B )
{
  int t;
  int index = -1;
  int maxLen = 0;

  for( t = 0; t < B->IP->nNumMotifTypes; t++ )
    {
      if(  MotifWidth( B, t ) > maxLen )
	{
	  index = t;
	  maxLen =  MotifWidth( B, t );
	}
    }

  return index;
}


int MaxSequenceLen( Model B )
{
  int len;
  int maxLen = -1;
  int seq;

   for(seq = 0; seq < B->IP->nNumSequences; seq++) 
     {
       len = SequenceLength( B, seq );
       if( len > maxLen )
	 maxLen = len;
     }

   return maxLen;
}


int TotalNumMotifs( Model B )
{
  int t;
  int totalMotifs = 0;

  for( t = 0; t < B->IP->nNumMotifTypes; t++)     
    totalMotifs += NUMMOTIFS(B->IP->nNumMotifs[t]);

  return totalMotifs;
}


int MotifWidth( Model B, int t )
{
  if(B->IP->is_defined[cl_F]) 
    return( B->IP->nMotifLen[t] );
  else
    return( B->F->FragWidth[t] );
}


int SpeciesInc( Model B, int seq )
{
  if( B->IP->is_defined[cl_D] )
    {
      if( seq > B->Phylo->maxPhyloSeq )
	return 1;
      else
	return B->Phylo->phyloSpecies[0];
    }
  else
    {
      if( B->Phylo->treeCount == 0 )
	return 1;
      else if( seq > B->Phylo->maxPhyloSeq )
	return 1;
      else
	return B->Phylo->phyloSpecies[B->Phylo->phyloIndex[seq]];
    }
}


int SeqInc( Model B, int seq )
{
  if( B->IP->is_defined[cl_D] )
    {
      if( seq > B->Phylo->maxPhyloSeq )
	return 1;
      else
	return B->Phylo->phyloSpecies[0];
    }
  else
    {
      if( B->Phylo->treeCount == 0 )
	return 1;
      else if( seq > B->Phylo->maxPhyloSeq )
	return 1;
      else if( B->Phylo->phyloSpeciesSample[B->Phylo->phyloIndex[seq]] )
	return 1;
      else
	return B->Phylo->phyloSpecies[B->Phylo->phyloIndex[seq]];
    }
}


int OffsetInMASS( Model B, int seq )
{
  if( B->Phylo->treeCount == 0 )
    return -1;
  else if( seq > B->Phylo->maxPhyloSeq )
    return -1;
  else
    return (seq - B->Phylo->phyloTreeStartSeq[B->Phylo->phyloIndex[seq]]);
}


int PosInFirstMASSSeq( Model B, int seqPos )
{
  int  seq = 0;
  int  offset;

  seq = SequenceFromPosition( B, seqPos );
  if( ! IsPhyloSeq( B,seq )  )
    return seqPos;
  else 
    {
      offset = seqPos - SequenceStartPos( B, seq );	  
      return SequenceStartPos( B, B->Phylo->phyloTreeStartSeq[B->Phylo->phyloIndex[seq]] ) + offset;
    }
  }


int IsPhyloSeq( Model B, int seq )
{
  IPtype     IP;
  PhyloType  PH;

  IP = B->IP;
  PH = B->Phylo;

  return( (IP->is_defined[cl_D] || PH->treeCount > 0) && seq <= PH->maxPhyloSeq );
}


int IsPhyloSeqBySpecies( Model B, int seq )
{
  IPtype     IP;
  PhyloType  PH;

  IP = B->IP;
  PH = B->Phylo;

  return( (IP->is_defined[cl_D] || PH->treeCount > 0) && 
	  (seq <= PH->maxPhyloSeq && B->Phylo->phyloSpeciesSample[B->Phylo->phyloIndex[seq]]) );
}


int NumWidths( Model B )
{
  return (MaxPossibleWidth( B ) - MinPossibleWidth( B ) + 1);
}


int MaxPossibleWidth( Model B )
{
  int t;
  int maxWidth = -1;

  if( B->IP->is_defined[cl_F] )
    {
      for( t = 0; t < B->IP->nNumMotifTypes; t++ )
	{
	  if( B->IP->nMotifLen[t] > maxWidth )
	    maxWidth = B->IP->nMotifLen[t]; 
	}
    }
  else
    {
      for( t = 0; t < B->IP->nNumMotifTypes; t++ )
	{
	  if( B->IP->nMaxFragWidth[t] > maxWidth )
	    maxWidth = B->IP->nMaxFragWidth[t]; 
	}
    }

  return( maxWidth );
}


int MinPossibleWidth( Model B )
{
  int t;
  int minWidth = INT_MAX;

  for( t = 0; t < B->IP->nNumMotifTypes; t++ )
    {
      if( B->IP->nMotifLen[t] < minWidth )
	minWidth =  B->IP->nMotifLen[t]; 
    }

  return( minWidth );
}


/***************************************************************/
/*  InitRProb - initialize Probability array for recursion     */
/*                                                             */
/*  Initialize dProb[0...nSeqLen-1][0...nNumMotifTypes - 1]    */
/*                  [0...nMaxBlocks][Forward...Reverse]        */
/*                                                             */
/*  Note the endpoints - this wastes a couple of bytes, but    */
/*  makes keeping track of things easier.                      */
/*                                                             */
/***************************************************************/

void InitRProbStruct( Model B )
{
  IPtype IP;
  RPType RP;
  ALType AP;
  int    i;
  int    j;
  int    k;
  int    t;
  int    totalMotifs;
  int    nCount;
  double dTrans;
  int    nSeq;
  int    n;
    
  NEW( B->RP, 1, RProbStruct );

  IP = B->IP;
  RP = B->RP;

  RP->nMaxBlocks = IP->nMaxBlocks;  
  RP->nMinBlocks = IP->nMinBlocks;  

  RP->dExpBlkWt = DEF_BLOCK_WT;
  /*  B->RP->dKSampleMap = -DBL_MAX; */ /* BT 02/11/03 */
  RP->dKSampleMap = IP->dKSampleMap;  

  RP->nSeq = 0;

  RP->nMaxSeqLen = MaxSequenceLen( B );

  RP->dProb = NULL;

  dTrans =  IP->nNumMotifTypes *  IP->nNumMotifTypes;
  dTrans = 1.0;     /* BT 9/6/2000 - test */
  NEWP3( RP->dPTrans, IP->nNumMotifTypes, double );
  NEWP3( RP->dInitialPTrans, IP->nNumMotifTypes, double );
  for( i = 0; i < IP->nNumMotifTypes; i++ )
    {
      NEWPP( RP->dPTrans[i], 2, double );
      NEWP( RP->dPTrans[i][FORWARD], IP->nNumMotifTypes, double );
      NEWP( RP->dPTrans[i][REVERSE], IP->nNumMotifTypes, double );
      NEWPP( RP->dInitialPTrans[i], 2, double );
      NEWP( RP->dInitialPTrans[i][FORWARD], IP->nNumMotifTypes, double );
      NEWP( RP->dInitialPTrans[i][REVERSE], IP->nNumMotifTypes, double );
      for( j = 0; j < IP->nNumMotifTypes; j++ )
	{
	  NEW( RP->dPTrans[i][FORWARD][j], 2, double );
	  NEW( RP->dPTrans[i][REVERSE][j], 2, double );
	  NEW( RP->dInitialPTrans[i][FORWARD][j], 2, double );
	  NEW( RP->dInitialPTrans[i][REVERSE][j], 2, double );
	  RP->dPTrans[i][FORWARD][j][FORWARD] = 1.0 / dTrans;    
	  RP->dInitialPTrans[i][FORWARD][j][FORWARD] = 1.0 / dTrans;    
	  if(IP->RevComplement)
	    {
	      RP->dPTrans[i][FORWARD][j][REVERSE] = 1.0 / dTrans; 
	      RP->dPTrans[i][REVERSE][j][FORWARD] = 1.0 / dTrans; 
	      RP->dPTrans[i][REVERSE][j][REVERSE] = 1.0 / dTrans; 
	      RP->dInitialPTrans[i][FORWARD][j][REVERSE] = 1.0 / dTrans; 
	      RP->dInitialPTrans[i][REVERSE][j][FORWARD] = 1.0 / dTrans; 
	      RP->dInitialPTrans[i][REVERSE][j][REVERSE] = 1.0 / dTrans; 
	    }
	}
    }

  NEWP( RP->dPriorTransCnts, IP->nNumMotifTypes, double );
  NEWP( RP->dTransCnts,  IP->nNumMotifTypes, double );
  for( i = 0; i < IP->nNumMotifTypes; i++ )
    {
      NEW( RP->dPriorTransCnts[i], IP->nNumMotifTypes, double );
      NEW( RP->dTransCnts[i], IP->nNumMotifTypes, double );
    }

  /* Gap Penalties - dGap[d][t] - prob of motif starting a distance d after end of motif */
  /* type t.                                                                             */

  RP->nMaxGap = RP->nMaxSeqLen - 1; 
  NEWPP( RP->dPGap, IP->nNumMotifTypes, double );
  for( i = 0; i < IP->nNumMotifTypes; i++ )
    {
      NEWP(  RP->dPGap[i], IP->nNumMotifTypes, double );
      for( j = 0; j < IP->nNumMotifTypes; j++ )
	{
	  NEW(  RP->dPGap[i][j],  RP->nMaxGap + 1, double );
	  for( k = 0; k <=  RP->nMaxGap; k++ )
	    RP->dPGap[i][j][k] = 1.0;   /* BT 09/12/06 */  /* Set defualt penalty to 1 */
	  /*	    RP->dPGap[i][j][k] = 1.0 / (RP->nMaxGap + 1);  */
	}
    }
  
  NEWP(  RP->nMaxGapWidth, IP->nNumMotifTypes, int );
  for( j = 0; j < IP->nNumMotifTypes; j++ )
    {
      NEW( RP->nMaxGapWidth[j], IP->nNumMotifTypes, int );
      for( k = 0; k < IP->nNumMotifTypes; k++ )
	RP->nMaxGapWidth[j][k] = 0.0;
    }

  NEWP( RP->dPEndGap, RP->nMaxGap + 1, double );
  for( i = 0; i <= RP->nMaxGap; i++ )
    {
      NEW(  RP->dPEndGap[i], IP->nNumMotifTypes, double );
      for( j = 0; j < IP->nNumMotifTypes; j++ )
	{
	  RP->dPEndGap[i][j] = 1.0;         
	}
    }

  NEWP( RP->dPSpacing, IP->nNumSequences, double );
  for( i = 0; i < IP->nNumSequences; i++ )
    {
      NEW( RP->dPSpacing[i], SequenceLength( B, i ), double );
    }

  RP->prSum = NULL;

  if( IP->is_defined[cl_E] )
    {
      if( ! IP->nNumMotifs[0][FORWARD] )
	{
	  /*	  totalMotifs = IP->nNumMotifTypes * IP->nNumSequences; */
	  totalMotifs = IP->nNumSequences;
	  RP->dPriorSites = totalMotifs; 
	  nCount = 0;
	  while( nCount < RP->dPriorSites )
	    {
	      for( t = 0; t <  IP->nNumMotifTypes && nCount < RP->dPriorSites; t++ )
		{
		  IP->nNumMotifs[t][FORWARD]++;
		  nCount++;
		}
	      
	      if( IP->RevComplement &&  nCount < RP->dPriorSites ) 
		{
		  for( t = 0; t <  IP->nNumMotifTypes && nCount < RP->dPriorSites; t++ )
		    {
		      IP->nNumMotifs[t][REVERSE]++;
		      nCount++;
		    }
		}
	    }	  
	}
      else
	{
	  totalMotifs = 0;
	  for( t = 0; t <  IP->nNumMotifTypes; t++ )
	    {
	      totalMotifs += IP->nNumMotifs[t][FORWARD];
	      totalMotifs += IP->nNumMotifs[t][REVERSE];
	    }
	  RP->dPriorSites = totalMotifs; 
	}
      
      RP->dPriorSeq = IP->nNumSequences;
    }
  else if( IP->site_samp )
    {
      totalMotifs = IP->nNumMotifTypes * IP->nNumSequences;
      RP->dPriorSites = totalMotifs; 
      RP->dPriorSeq = IP->nNumSequences;
    }
  else
    {
      for( totalMotifs = 0, t = 0; t < IP->nNumMotifTypes; t++ )
	totalMotifs +=  NUMMOTIFS( IP->nNumMotifs[t] ); 
    
      RP->dPriorSites = totalMotifs; 
      RP->dPriorSeq = IP->nNumSequences;
    }

  RP->dSiteProbCons = 0.5;
  RP->dSiteProbNonCons = 0.5;
  RP->dSiteProbWt = 1;

  NEW( RP->dExpBlkCnt, RP->nMaxBlocks + 1, double );
  for( i = 0; i <= RP->nMaxBlocks; i++ )                      /* 03/11/03 */
    RP->dExpBlkCnt[i] = 1.0 / (RP->nMaxBlocks + 1);    /* changed default from  DEF_EXPBLK_CNT */

  NEWP( RP->dProbCons, B->IP->nNumSequences, double );
  for( i = 0; i < B->IP->nNumSequences; i++ )
    {
      NEW(RP->dProbCons[i], SequenceLength( B, i ), double );
      for( j = 0; j < SequenceLength( B, i ); j++ )
	RP->dProbCons[i][j] = DBL_MAX;
    }

  NEWP( RP->dProbSite, B->IP->nNumSequences, double );
  NEWP( RP->dProbNonSite, B->IP->nNumSequences, double );
  for( i = 0; i < B->IP->nNumSequences; i++ )
    {
      NEW(RP->dProbSite[i], SequenceLength( B, i ), double );
      NEW(RP->dProbNonSite[i], SequenceLength( B, i ), double );
    }

  NEW( RP->nSites, B->IP->nNumSequences, int );

  NEW( RP->nSeqPos, B->IP->nSeqLen, int );
  nSeq = 0;
  for( i = 0; i < B->IP->nSeqLen; i++ )
    {
      if( i > SequenceEndPos( B, nSeq ) )
	  nSeq++;
      RP->nSeqPos[i] = nSeq;
    }

  NEW( RP->nInMotif, B->IP->nSeqLen, int );

  NEWPP( RP->sitePos,  B->IP->nNumSequences, SiteStruct );
  for( i = 0; i < B->IP->nNumSequences; i++ )
    {
      NEWP(RP->sitePos[i], SequenceLength( B, i ), SiteStruct );
      for( j = 0; j < SequenceLength( B, i ); j++ )
	{
	  NEW( RP->sitePos[i][j], IP->nNumMotifTypes, SiteStruct );
	}
    }  
  
  NEW( RP->alphaCons, B->IP->nNumMotifTypes, double );
  NEW( RP->alphaNonCons, B->IP->nNumMotifTypes, double );
  NEW( RP->betaCons, B->IP->nNumMotifTypes, double );
  NEW( RP->betaNonCons, B->IP->nNumMotifTypes, double );

  NEW( RP->dEndCnts, IP->nNumMotifTypes, double );
  NEW( RP->dPriorEndCnts, IP->nNumMotifTypes, double );
  NEW( RP->dEndSiteProb, IP->nNumMotifTypes, double );
  NEW( RP->dBeginCnts, IP->nNumMotifTypes, double );
  NEW( RP->dPriorBeginCnts, IP->nNumMotifTypes, double );
  NEW( RP->dBeginSiteProb, IP->nNumMotifTypes, double );
  NEW( RP->dInitialEndSiteProb, IP->nNumMotifTypes, double );
  NEW( RP->dInitialBeginSiteProb, IP->nNumMotifTypes, double );
  for( j = 0; j < IP->nNumMotifTypes; j++ )
    {
/*       RP->dInitialEndSiteProb[j] = 1.0 / IP->nNumMotifTypes; */
/*       RP->dInitialBeginSiteProb[j] = 1.0 / IP->nNumMotifTypes; */
/*       RP->dEndSiteProb[j] = 1.0 / IP->nNumMotifTypes; */
/*       RP->dBeginSiteProb[j] = 1.0 / IP->nNumMotifTypes; */
      RP->dInitialEndSiteProb[j] = 1.0;
      RP->dInitialBeginSiteProb[j] = 1.0;
      RP->dEndSiteProb[j] = 1.0;
      RP->dBeginSiteProb[j] = 1.0;
     }

  NEWP( RP->dInitialPosMatrix, RP->nMaxBlocks, double );
  NEWP( RP->dPriorPosMatrixCnts, RP->nMaxBlocks, double );
  NEWP( RP->dPosMatrixCnts, RP->nMaxBlocks, double );
  NEWP( RP->dPosMatrix, RP->nMaxBlocks, double );
  for( i = 0; i <  RP->nMaxBlocks; i++ )
    {
      NEW( RP->dInitialPosMatrix[i], IP->nNumMotifTypes, double );
      NEW( RP->dPriorPosMatrixCnts[i], IP->nNumMotifTypes, double );
      NEW( RP->dPosMatrixCnts[i], IP->nNumMotifTypes, double );
      NEW( RP->dPosMatrix[i], IP->nNumMotifTypes, double );
    }

  NEW( B->AP, 1, AlignStruct );

  AP = B->AP;
  AP->nMaxBlocks = IP->nMaxBlocks;    /* Temporary */

  NEWP( AP->dAlignCount, AP->nMaxBlocks + 1, double );
  for( i = 0; i <= AP->nMaxBlocks; i++ )
    {
      NEW( AP->dAlignCount[i], IP->nNumSequences, double );
      for( j = 0; j < IP->nNumSequences; j++ )
	AP->dAlignCount[i][j] = -DBL_MAX;
    }
  for( j = 0; j < IP->nNumSequences; j++ )
    AP->dAlignCount[0][j] = 1;

  NEWP( AP->dAlignCountSave, AP->nMaxBlocks + 1, double );
  for( i = 0; i <= AP->nMaxBlocks; i++ )
    {
      NEW( AP->dAlignCountSave[i], IP->nNumSequences, double );
    }

  NEW( AP->dAlignWt, AP->nMaxBlocks + 1, double );
  for( i = 0; i <= AP->nMaxBlocks; i++ )
    {
      AP->dAlignWt[i] = 1.0;
    }

  RP->nonCutoffSites = B->IP->nSeqLen;
  /* RP->dMinSiteMap = DBL_MAX;
     RP->nMinSiteMotifs = INT_MAX; */
  RP->dMinSiteMap = -DBL_MAX; /* BT 04/09/04 */
  RP->nMinSiteMotifs = 0; 

  RP->nUseFixedBlockProbs = TRUE;  /* BT 05/02/04 */

  NEWP( RP->nPriorMotifSites, IP->nNumMotifTypes, int );
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      NEW( RP->nPriorMotifSites[t], 2, int );      
    }

  NEW( RP->probK, RP->nMaxBlocks + 1, double );
  NEW( RP->priorSitesPerSeq, RP->nMaxBlocks + 1, double );
  for( k = 0; k <= RP->nMaxBlocks; k++ )
    RP->priorSitesPerSeq[k] = RP->dExpBlkWt * RP->dPriorSeq * RP->dExpBlkCnt[k];

  if( IP->is_defined[cl_bayes] )
    {
      NEWP( RP->kTotalCounts, IP->nNumSequences, int );
      NEWPP( RP->totalPosCounts, IP->nNumSequences, int );
      for( n = 0; n < IP->nNumSequences; n++ )
	{
	  NEW( RP->kTotalCounts[n], RP->nMaxBlocks + 1, int );

	  NEWP( RP->totalPosCounts[n], SequenceLength( B, n), int );
	  for( i = 0; i < SequenceLength( B, n); i++ )
	    NEW( RP->totalPosCounts[n][i], NumWidths( B ), int);
	}

      if( ! IP->is_defined[cl_no_cred] )
	{
	  NEWP3( RP->samples, IP->nNumSequences, SampleStruct );
	  for( n = 0; n < IP->nNumSequences; n++ )
	    {
	      NEWPP( RP->samples[n], IP->nSeeds, SampleStruct );
	      for( i = 0; i < IP->nSeeds; i++ )
		{
		  NEWP( RP->samples[n][i], IP->bayesSamplePeriod, SampleStruct );
		  for( j = 0; j < IP->bayesSamplePeriod; j++ )
		    NEW( RP->samples[n][i][j], RP->nMaxBlocks, SampleStruct );
		}
	    }
	}
    }
}


void InitRProbArray( Model B, int nSeq )
{
  IPtype IP;
  RPType RP;
  int    j;
  int    k;
  int    d;
    
  IP = B->IP;
  RP = B->RP;

  RP->nSeq = nSeq;
  RP->nSeqLen = SequenceLength( B, nSeq );

  NEWP3( RP->dProb, RP->nMaxBlocks + 1, double );
  for( k = 0; k <= RP->nMaxBlocks; k++ )
    {
      NEWPP( RP->dProb[k], IP->nNumMotifTypes, double );
      for( j = 0; j < IP->nNumMotifTypes; j++ )
	{
	  NEWP( RP->dProb[k][j], 2, double );
	  for( d = FORWARD; d <= REVERSE; d++ )
	    {
	      NEW( RP->dProb[k][j][d],  RP->nSeqLen, double );
	    }
	}
    }

  NEWP( RP->dPFootprint, IP->nNumMotifTypes, double );
  for( j = 0; j < IP->nNumMotifTypes; j++ )
    {
      NEW( RP->dPFootprint[j],  RP->nSeqLen, double );
    }

  NEW( RP->prSum, RP->nMaxBlocks + 1, double );
}


/***************************************************************/
/*  FreeRProb - free Probability arrays                        */
/*                                                             */
/***************************************************************/

void FreeRProb( Model B )
{
  IPtype IP;
  RPType RP;
  int    i;
  int    d;
  int    j;
  int    k;

  if( B->RP == NULL )
    return;
    
  IP = B->IP;
  RP = B->RP;

  if( RP->dProb != NULL )
    {
      for( k = 0; k <= RP->nMaxBlocks; k++ )
	{
	  for( j = 0; j < IP->nNumMotifTypes; j++ )
	    {
	      for( d = FORWARD; d <= REVERSE; d++ )
		free( RP->dProb[k][j][d] );
	      free( RP->dProb[k][j] );
	    }
	  free( RP->dProb[k] );
	}
      free( RP->dProb );
      RP->dProb = NULL;  
    }

  for( i = 0; i < IP->nNumMotifTypes; i++ )
    {
      for( j = 0; j < IP->nNumMotifTypes; j++ )
	{
 	  free( RP->dPTrans[i][FORWARD][j] );
 	  free( RP->dPTrans[i][REVERSE][j] );
 	  free( RP->dInitialPTrans[i][FORWARD][j] );
 	  free( RP->dInitialPTrans[i][REVERSE][j] );
	}
      free( RP->dPTrans[i][FORWARD] );
      free( RP->dPTrans[i][REVERSE] );
      free( RP->dInitialPTrans[i][FORWARD] );
      free( RP->dInitialPTrans[i][REVERSE] );
      
      free( RP->dPTrans[i] );
      free( RP->dInitialPTrans[i] );
      free( RP->dTransCnts[i] );
      free( RP->dPriorTransCnts[i] );
    }
  free( RP->dPTrans );
  RP->dPTrans = NULL;
  free( RP->dInitialPTrans );
  RP->dInitialPTrans = NULL;
  free( RP->dTransCnts );
  RP->dTransCnts = NULL;
  free( RP->dPriorTransCnts );
  RP->dPriorTransCnts = NULL;

  free( RP->nSeqPos );
  free( RP->nInMotif );

  FREEPP( RP->dPGap,  IP->nNumMotifTypes, IP->nNumMotifTypes );
  FREEP( RP->dPEndGap, RP->nMaxGap + 1 );
  FREEP( RP->dPSpacing, IP->nNumSequences );
  FREEP( RP->nMaxGapWidth, IP->nNumMotifTypes );
  free( RP->dExpBlkCnt );

  if( RP->prSum != NULL )
    free( RP->prSum );

  for( i = 0; i < B->IP->nNumSequences; i++ )
    {
      FREEP( RP->sitePos[i], SequenceLength( B, i ) );
    }
  free( RP->sitePos );

  FREEP( RP->dProbCons, B->IP->nNumSequences );
  FREEP( RP->dProbSite, B->IP->nNumSequences );
  FREEP( RP->dProbNonSite, B->IP->nNumSequences );
  free( RP->dEndSiteProb );
  free( RP->dInitialEndSiteProb );
  free( RP->dEndCnts );
  free( RP->dPriorEndCnts );
  free( RP->dBeginSiteProb );
  free( RP->dInitialBeginSiteProb );
  free( RP->dBeginCnts );
  free( RP->dPriorBeginCnts );

  if( RP->nSites != NULL )
    free( RP->nSites );

  if( RP->alphaCons != NULL )
    free( RP->alphaCons );
  if( RP->alphaNonCons != NULL )
    free( RP->alphaNonCons );
  if( RP->betaCons != NULL )
    free( RP->betaCons );
  if( RP->betaNonCons != NULL )
    free( RP->betaNonCons );

  if( RP->dInitialPosMatrix != NULL )
    FREEP( RP->dInitialPosMatrix, RP->nMaxBlocks );
  if( RP->dPriorPosMatrixCnts != NULL )
    FREEP( RP->dPriorPosMatrixCnts, RP->nMaxBlocks );
  if( RP->dPosMatrixCnts != NULL )
    FREEP( RP->dPosMatrixCnts, RP->nMaxBlocks );
  if( RP->dPosMatrix != NULL )
    FREEP( RP->dPosMatrix, RP->nMaxBlocks );

  FREEP( RP->nPriorMotifSites, IP->nNumMotifTypes );
  free( RP->priorSitesPerSeq );
  free( RP->probK );

  if( IP->is_defined[cl_bayes] )
    {
      for( i = 0; i < IP->nNumSequences; i++ )
	{
	  FREEP( RP->totalPosCounts[i], SequenceLength( B, i ) );
	  if( ! IP->is_defined[cl_no_cred] )
	    FREEPP( RP->samples[i], IP->nSeeds, IP->bayesSamplePeriod );
	}
      if( ! IP->is_defined[cl_no_cred] )
	free( RP->samples );
      free( RP->totalPosCounts );
      FREEP( RP->kTotalCounts, IP->nNumSequences );
    }

  free( B->RP );
  B->RP = NULL;
}


void FreeRProbArray( Model B )
{
  IPtype IP;
  RPType RP;
  int    j;
  int    k;
  int    d;

  if( B->RP == NULL )
    return;
    
  IP = B->IP;
  RP = B->RP;

  for( k = 0; k <= RP->nMaxBlocks; k++ )
    {
      for( j = 0; j < IP->nNumMotifTypes; j++ )
	{
	  for( d = FORWARD; d <= REVERSE; d++ )
	    free( RP->dProb[k][j][d] );
	  free( RP->dProb[k][j] );
	}
      free( RP->dProb[k] );
    }
  free( RP->dProb );

  free( RP->prSum );

  FREEP( RP->dPFootprint, IP->nNumMotifTypes );

  RP->prSum = NULL;
  RP->dProb = NULL;  
}


void CalcRProbArray( Model B,  PoSition **Pos, int nSeq, ProbStruct *P, int iter )
{
  if( ! B->IP->is_defined[cl_z] )
    CalcRProbArraySpacingNoOPt( B, Pos, nSeq, P, iter );
  else
    CalcRProbArraySpacing( B, Pos, nSeq, P, iter );
}


void CalcRProbArraySpacing( Model B,  PoSition **Pos, int nSeq, ProbStruct *P, int iter  )
{
  IPtype     IP;
  RPType     RP;
  ALType     AP;
  int        n;
  int        nLen;
  int        nEndPos;
  int        t;
  int        dir;
  int        j;
  int        t1;
  double     sum;
  int        k;
  int        minMotif;
  int        nStart;
  double     pr1;
  double     **dCutOff;
  double     dGapPr;
  double     dMotifProb;
  double     dBGProb;
  double     fpSiteProb;
  double     fpBGProb;
  double     fpSum;
  int        nGap;
  int        nMinMotifWidth;
  int        nMotifWidth;
  double     sum2;
  int        nPartialSumIndex;
  double     **dPrevSum;
  double     **dPrevNorm;
  int        bPossiblePos;
  int        bCutOff;
  double     dBGProbSum = 1.0;
#ifdef _DEBUG_
  double     prf[GRAPH_SIZE];   /* DEBUG */
  double     prr[GRAPH_SIZE];
  double     fpr[GRAPH_SIZE];
#endif
 
  IP = B->IP;
  RP = B->RP;
  AP = B->AP;

  NEWP( dCutOff, IP->nNumMotifTypes, double );
  NEWP( dPrevSum, IP->nNumMotifTypes, double );
  NEWP( dPrevNorm, IP->nNumMotifTypes, double );
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      NEW( dPrevSum[t], 2, double );
      NEW( dPrevNorm[t], 2, double );
      NEW( dCutOff[t], 2, double );
      dCutOff[t][FORWARD] = -DBL_MAX;
      dCutOff[t][REVERSE] = -DBL_MAX;
      dBGProbSum -= IP->dposterior_prob[t];
    }
         
  /* Calculate first for 1 block */
  nLen = SequenceLength( B, nSeq );
  nEndPos = nLen - 1;
  nStart = SequenceStartPos( B, nSeq );
    
  if( B->IP->is_defined[cl_T] )
    {
      UpdateFootprintParameters( B, Pos, nSeq );
    }

  for( t = 0; t < B->IP->nNumMotifTypes; t++ ) 
    {
      if( B->IP->is_defined[cl_T] )
	{
	  fpSum = 0.0;
	  for( n = 0; n < nLen; n++ )
	    {
	      RP->dPFootprint[t][n] = SiteProbFromFootprint( B, n, t, nSeq, 
							     &fpSiteProb, &fpBGProb );
	      fpSum += RP->dPFootprint[t][n];
	    }

	  /* for( n = 0; n < nLen; n++ )   
	     RP->dPFootprint[t][n] /= fpSum;  */  /* TEST */
	  
	}
      else if( RP->nUseSpacingProb )
	{
	  fpSum = RP->sitePos[nSeq][nLen-1][t].dPartialSum;
	  for( n = 0; n < nLen; n++ )   
	    {
	      RP->dPFootprint[t][n] = RP->sitePos[nSeq][n][t].dSpacingProb;
	    }
	  
	  for( n = 0; n < nLen; n++ )   
	    RP->dPFootprint[t][n] /= fpSum;
	}      
      else
	{
	  fpSum = 1.0;
	  dBGProb = IP->dposterior_prob[t] / dBGProbSum;
	  for( n = 0; n < nLen; n++ )   
	    {
	      if( RP->bUseFlatGap )  
		RP->dPFootprint[t][n] = RP->dFlatGapProb;
	      else 
		RP->dPFootprint[t][n] = 1.0 / nLen;  /* TEST !!!  10/12/04 */
	    }
	}
    }
  
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      nMotifWidth = MotifWidth(B, t);
      for( n = 0; n <= nEndPos - nMotifWidth + 1; n++ )
	{
	  if( PossibleStartPos( Pos, n + nStart, t, nLen, nSeq, B ) )
	    {
	      in_motif_prob( B, *P, n + nStart, t, FALSE, TRUE, 
			      &dMotifProb, &dBGProb );
	      RP->dProb[1][t][FORWARD][n] = dMotifProb / dBGProb;

	      if(IP->RevComplement)
		{
		  in_motif_prob( B, *P, n + nStart, t, TRUE, TRUE, 
				 &dMotifProb, &dBGProb );
		  RP->dProb[1][t][REVERSE][n] = dMotifProb / dBGProb;  /* BT 04/16/03 */
		}

	      if( FALSE && RP->bUsePosMatrix )
		{
		  RP->dProb[1][t][FORWARD][n] *= RP->dPosMatrix[0][t];
		  RP->dProb[1][t][REVERSE][n] *= RP->dPosMatrix[0][t];
		}

	      /*	      RP->dProb[1][t][FORWARD][n] *= (fp[t][n] / fpSum[t]);
			      RP->dProb[1][t][REVERSE][n] *= (fp[t][n] / fpSum[t]); */
	      
	        RP->prSum[1] += (RP->dProb[1][t][FORWARD][n] + RP->dProb[1][t][REVERSE][n])
		  * RP->dEndSiteProb[t] * RP->dBeginSiteProb[t]; /* BT 4/16/03 */
		/* RP->prSum[1] += (RP->dProb[1][t][FORWARD][n] + RP->dProb[1][t][REVERSE][n])
		 * RP->dBeginSiteProb[t]; BT 06/15/04 */

	      if( RP->dProb[1][t][FORWARD][n] > dCutOff[t][FORWARD] )
		dCutOff[t][FORWARD] =  RP->dProb[1][t][FORWARD][n];
	      if( RP->dProb[1][t][REVERSE][n] > dCutOff[t][REVERSE] )
		dCutOff[t][REVERSE] =  RP->dProb[1][t][REVERSE][n];
	    }
	}
    }
    
#ifdef _DEBUG_
  memset( prf, 0, sizeof( prf ) );
  memset( prr, 0, sizeof( prr ) );
  memset( fpr, 0, sizeof( fpr ) );
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      for( j = 0; j < min( nLen, GRAPH_SIZE ) ; j++ )      /* DEBUG */
	{
	  prf[j] = max( prf[j], RP->dProb[1][t][FORWARD][j] );
	  prr[j] = max( prr[j],  RP->dProb[1][t][REVERSE][j]) ;
	  fpr[j] = max( fpr[j], RP->dPFootprint[t][j] );
	}
    }
  if( nSeq == 0 )     /* DEBUG */
    nSeq = nSeq;      
#endif

  RP->prTotalSum = RP->prSum[1];

  minMotif = MinLenMotif( B );

  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      dCutOff[t][FORWARD] *= IP->dRCutoff;
      dCutOff[t][REVERSE] *= IP->dRCutoff;
    }
  
  /* Now calcultate for other block sizes */
  nMinMotifWidth = MotifWidth( B, minMotif);
  for( k = 2; k <= RP->nMaxBlocks; k++ )
    {
      if( RP->dProbCons[nSeq][0] > 1.0 )
	{
	  if( nLen - (k-1) * nMinMotifWidth > 0 )
	    dGapPr = 1.0 / (nLen - (k-1) * nMinMotifWidth + 1.0);
	  else
	    dGapPr = 0;
	}
      else  /* TEST */  /* 04/06/04 */ /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
	dGapPr = 1.0;
      
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  nMotifWidth = MotifWidth( B, t );
	  for( dir = FORWARD; 
	       dir <= REVERSE * ( IP->RevComplement ? 1 : 0); 
	       dir++ )
	    {
	      for( t1 = 0; t1 < IP->nNumMotifTypes; t1++ )
		{
		  dPrevSum[t1][FORWARD] = 0;
		  dPrevSum[t1][REVERSE] = 0;
		  dPrevNorm[t1][FORWARD] = 0;
		  dPrevNorm[t1][REVERSE] = 0;
		}
	      bCutOff = FALSE;
	      nPartialSumIndex = 0;
	      for(nGap = 0, j = (k-1) *  nMinMotifWidth; 
		  j <= nEndPos - nMotifWidth + 1; 
		  j++, nGap++ )
		{	
		  pr1 = RP->dProb[1][t][dir][j];
		  sum = 0.0;
		  if( /* RP->bUseFlatGap || */ pr1 > dCutOff[t][dir] ) 
		    {
		      bPossiblePos =  PossibleStartPos( Pos, j + nStart, t, nLen, nSeq, B );
		      if( RP->bUseFlatGap || bPossiblePos )
			{
			  sum = CalcGapSum( B, k, j, t, dir, dGapPr, 
					    (k-1) *  nMinMotifWidth, dPrevSum, dPrevNorm,
					    bCutOff); 
			  bCutOff = FALSE;
			}
		      
		      if( bPossiblePos )
			{
			  sum2 = sum * pr1;
			  RP->dProb[k][t][dir][j] = sum2; 
#ifdef _DEBUG_
			  if( ! finite( RP->dProb[k][t][dir][j] ) )
			    {
			      sum2 = sum2;
			    }			  
#endif			  
			  RP->prSum[k] += sum2 * RP->dEndSiteProb[t] * 
			     RP->dBeginSiteProb[t];  /* BT 04/16/03 */
			  /* RP->prSum[k] += sum2 * RP->dEndSiteProb[t];   BT 10/20/03 */
			    }
		    }
		  else
		    bCutOff = TRUE;
		}
	    }
	}
      
#ifdef _DEBUG_
      memset( prf, 0, sizeof( prf ) );
      memset( prr, 0, sizeof( prr ) );
      memset( fpr, 0, sizeof( fpr ) );
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  for( j = 0; j < min( nLen, GRAPH_SIZE ) ; j++ )      /* DEBUG */
	    {
	      prf[j] = max( prf[j], RP->dProb[k][t][FORWARD][j] );
	      prr[j] = max( prr[j], RP->dProb[k][t][REVERSE][j] );
	      fpr[j] = max( fpr[j], RP->dPFootprint[t][j] );
	    }
	}
      if( nSeq == 0 )     /* DEBUG */
	nSeq = nSeq;      
#endif

      RP->prTotalSum += RP->prSum[k];
    }

  FREEP( dCutOff, IP->nNumMotifTypes );
  FREEP( dPrevSum, IP->nNumMotifTypes );
  FREEP( dPrevNorm, IP->nNumMotifTypes );  
  }


double CalcGapSum( Model B, int k, int j, int t, int dir, double dGapPr,
		   int nStart, 
		   double **dPrevSum, double **dPrevNorm,
		   int bCutOff )
{
  double     sum2;
  double     *pdPGap;
  double     *fpt1;
  int        nPosDiff;
  int        nPosDiff2;
  double     dGapWt;
  RPType     RP;
  IPtype     IP;
  int        t1;
  int        d;
  int        v;
  double     *pr_kt1;
  int        nMotifWidtht1;
  int        nMaxGapWidth;
  double     dPr1 =0;
  double     sum = 0.0;
  double     partSum;
  double     normSum;
  double     pSum;
  double     nSum;
  int        nSub;
  double     ***prk;
  double     **prkt;
  double     **dPTrans;
  double     *dPTranst;
  double     *dPrevSumt1;
  double     *dPrevNormt1;
  int        foR;

  RP = B->RP;
  IP = B->IP;

  prk = RP->dProb[k-1];
  dPTrans = RP->dPTrans[t][dir];    
  foR = REVERSE * ( IP->RevComplement ? 1 : 0);

  for( t1 = 0; t1 < IP->nNumMotifTypes; t1++ )
    {
      nMaxGapWidth = RP->nMaxGapWidth[t1][t];
      nMotifWidtht1 = MotifWidth( B, t1 );
      nPosDiff = j - nMotifWidtht1 - 1;  /* BT 07/30/2002 */
      if( nPosDiff >= 0 )  /* BT 1/29/2001 */
	{
      /*      nPosDiff2 = max( 0, j - nMaxGapWidth - nMotifWidtht1 ); */
	  nPosDiff2 = max( max( 0, nStart - nMotifWidtht1 ), j - nMaxGapWidth - nMotifWidtht1 );
	  nSub = max( 0, nPosDiff2 - 1);
	  pdPGap = RP->dPGap[t1][t];
	  fpt1 =  RP->dPFootprint[t1];
	  prkt = prk[t1]; 
	  dPTranst = dPTrans[t1]; 
	  dGapWt = dGapPr; 
	  dPrevSumt1 =  dPrevSum[t1];
	  dPrevNormt1 = dPrevNorm[t1];
	  if( RP->bUsePosMatrix )
	    dPr1= RP->dPosMatrix[k-2][t1];
	  for( d = FORWARD; d <= foR; d++ )
	    {	  
	      pr_kt1 = prkt[d];     /* p(R|k-1,t1,d)  */
	      if( ! RP->bUsePosMatrix )
		{
		  dPr1 = dPTranst[d];   /* dPTrans[t][dir][t1][d] probability of t following t1 */
		  /* dPr1 *= 4; */ /* TEST */ /* BT 01/28/03 */
		}

	      if( FALSE && ! bCutOff && RP->bUseFlatGap )  /* TEST */  /* BT 07/03/2001 */
		{
		  partSum = pr_kt1[nPosDiff] * fpt1[nPosDiff];
		  normSum = fpt1[nPosDiff];
		  if( j - nMaxGapWidth - nMotifWidtht1 > 0 && j - nMaxGapWidth > nStart )
		    {
		      partSum -= pr_kt1[nSub] * fpt1[nSub];
		      normSum -= fpt1[nSub];
		    }
		  partSum *= dPr1;		  
		  partSum += dPrevSumt1[d];
		  normSum += dPrevNormt1[d];
		  partSum = max( 0, partSum );  /* BT 07/02/2001 */
		  normSum = max( 0, normSum ); 
		  if( normSum > 0.0 )
		    sum += partSum / normSum;
		  dPrevSumt1[d] = partSum;
		  dPrevNormt1[d] = normSum;
		}	  	  
	      else if( ! RP->nUseGap )
		{
		  for( sum2 = 0.0, v = 0; 
		       v < j -  nMotifWidtht1 + 1; v++ ) 
		    {
		     /* sum2 += pr_kt1[v] * fpt1[v]; */ /* BT 2/15/2001 */
		      sum2 += pr_kt1[v];  /* BT 12/14/2000 */
		    }
		  sum += sum2 * dGapWt * dPr1;
		}
	      else
		{
		  sum2 = GapProduct( B, nPosDiff2, nPosDiff, pr_kt1, 
				     pdPGap, fpt1, &pSum, &nSum ); 
		  sum += sum2 * dPr1;
		  dPrevSum[t1][d] = pSum * dPr1 / pdPGap[0];
		  dPrevNorm[t1][d] = nSum / pdPGap[0]; 
		}
#ifdef _DEBUG_
	      if( ! finite( sum ) || sum < 0.0 )
		sum = sum;
#endif
	    }
	}
    }

  return( sum );
}


/* Unroll the loop to calculate the gap probabilities in CalcRProbArraySpacing */
/* that loop is where the routine spends most of its time                      */
/* Parameters:                                                                 */
/*             nStart - start of gap                                           */
/*             nEnd - length of gap + 1, gap runs from nStart to nEnd          */
/*             dProb - address of starrting point in prob array                */
/*             dGap - gap probability array                                    */
double GapProduct( Model B, int nStart, int nEnd, double *dProb, double *dGap, double *fp,
		   double *pSum, double *nSum )
{
  double sum = 0.0;
  double normsum = 0.0;
  int    v;
  int    nLimit;
  double *dGapProb;
  int    i;
  
  nLimit = nStart + 8 * ((nEnd - nStart + 1) / 8);
  v = nStart;
  i = 0;
  if( B->RP->bUseFlatGap )
    {
      dGapProb = &fp[nStart];

      while( v < nLimit )
	{
	  normsum += dGapProb[i] +
	    dGapProb[i + 1] +
	    dGapProb[i + 2] +
	    dGapProb[i + 3] +
	    dGapProb[i + 4] +
	    dGapProb[i + 5] +
	    dGapProb[i + 6] +
	    dGapProb[i + 7];
	  
	  i += 8;
	  v += 8;
	}
      
      if( v <= nEnd )
	{
	  for( ; v <= nEnd; v++, i++ )
	    {
	      normsum += dGapProb[i];
	    }
	}
    }
  else 
    {
      NEW( dGapProb, nEnd - nStart + 1, double );

      while( v < nLimit )
	{
	  dGapProb[i] = dGap[nEnd - v] * fp[v];
	  dGapProb[i + 1] = dGap[nEnd - v - 1] * fp[v + 1];
	  dGapProb[i + 2] = dGap[nEnd - v - 2] * fp[v + 2];
	  dGapProb[i + 3] = dGap[nEnd - v - 3] * fp[v + 3];
	  dGapProb[i + 4] = dGap[nEnd - v - 4] * fp[v + 4];
	  dGapProb[i + 5] = dGap[nEnd - v - 5] * fp[v + 5];
	  dGapProb[i + 6] = dGap[nEnd - v - 6] * fp[v + 6];
	  dGapProb[i + 7] = dGap[nEnd - v - 7] * fp[v + 7]; 
	  
	  normsum += dGapProb[i] +
	    dGapProb[i + 1] +
	    dGapProb[i + 2] +
	    dGapProb[i + 3] +
	    dGapProb[i + 4] +
	    dGapProb[i + 5] +
	    dGapProb[i + 6] +
	    dGapProb[i + 7];
	  
	  i += 8;
	  v += 8;
	}
      
      if( v <= nEnd )
	{
	  for( ; v <= nEnd; v++, i++ )
	    {
	      dGapProb[i] = dGap[nEnd - v] * fp[v]; 
	      normsum += dGapProb[i];
	    }
	}
    }

  v = nStart;
  i = 0;
  while( v < nLimit )
    {
      sum += dProb[v] * dGapProb[i] +
	dProb[v + 1] * dGapProb[i + 1] +
	dProb[v + 2] * dGapProb[i + 2] +
	dProb[v + 3] * dGapProb[i + 3] +
	dProb[v + 4] * dGapProb[i + 4] +
	dProb[v + 5] * dGapProb[i + 5] +
	dProb[v + 6] * dGapProb[i + 6] +
	dProb[v + 7] * dGapProb[i + 7];

      v += 8;
      i += 8;
    }
  
  if( v <= nEnd )
    {
      for( ; v <= nEnd; v++, i++ )
	{
	  sum += dProb[v] * dGapProb[i];
	}
    }

  if( dGapProb != &fp[nStart] )
    free( dGapProb );

  *pSum = sum * B->RP->dFlatGapProb;
  *nSum = normsum * B->RP->dFlatGapProb;

#ifdef _DEBUG_
  if( sum > 0 )
    sum = sum;
#endif

  if( normsum > 0 )
    sum /= normsum;
    
  return( sum );
}


void CalcRProbArraySpacingNoOPt( Model B,  PoSition **Pos, int nSeq, ProbStruct *P, int iter  )
{
  IPtype     IP;
  RPType     RP;
  ALType     AP;
  int        n;
  int        nLen;
  int        nEndPos;
  int        t;
  int        dir;
  int        j;
  int        t1;
  double     sum;
  int        k;
  int        d;
  int        minMotif;
  int        nStart;
  double     **dCutOff;
  double     dGapPr;
  double     dMotifProb;
  double     dBGProb;
  double     fpSiteProb;
  double     fpBGProb;
  double     fpSum;
  int        nMotifWidth;
  double     sum2;
  double     dBGProbSum = 1.0;
  int        motifWidth;
  int        minMotifWidth;
  int        nPos;
  int        nMotifWidtht1;
  int        nMax;
#ifdef _DEBUG_
  double     prf[GRAPH_SIZE];   /* DEBUG */
  double     prr[GRAPH_SIZE];
  double     fpr[GRAPH_SIZE];
#endif

  IP = B->IP;
  RP = B->RP;
  AP = B->AP;

  NEWP( dCutOff, IP->nNumMotifTypes, double );
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      NEW( dCutOff[t], 2, double );
      dCutOff[t][FORWARD] = -DBL_MAX;
      dCutOff[t][REVERSE] = -DBL_MAX;
      dBGProbSum -= IP->dposterior_prob[t];
    }
         
  /* Calculate first for 1 block */
  nLen = SequenceLength( B, nSeq );
  nEndPos = nLen - 1;
  nStart = SequenceStartPos( B, nSeq );
    
  if( B->IP->is_defined[cl_T] )
    {
      UpdateFootprintParameters( B, Pos, nSeq );
    }

  for( t = 0; t < B->IP->nNumMotifTypes; t++ ) 
    {
      if( B->IP->is_defined[cl_T] )
	{
	  fpSum = 0.0;
	  for( n = 0; n < nLen; n++ )
	    {
	      RP->dPFootprint[t][n] = SiteProbFromFootprint( B, n, t, nSeq, 
							     &fpSiteProb, &fpBGProb );
	      fpSum += RP->dPFootprint[t][n];
	    }
	}
      else if( RP->nUseSpacingProb )
	{
	  fpSum = RP->sitePos[nSeq][nLen-1][t].dPartialSum;
	  for( n = 0; n < nLen; n++ )   
	    {
	      RP->dPFootprint[t][n] = RP->sitePos[nSeq][n][t].dSpacingProb;
	    }
	  
	  for( n = 0; n < nLen; n++ )   
	    RP->dPFootprint[t][n] /= fpSum;
	}      
      else
	{
	  fpSum = 1.0;
	  dBGProb = IP->dposterior_prob[t] / dBGProbSum;
	  for( n = 0; n < nLen; n++ )   
	    {
	      if( RP->bUseFlatGap )  
		RP->dPFootprint[t][n] = RP->dFlatGapProb;
	      else 
		RP->dPFootprint[t][n] = 1.0;		
	      /* else if( RP->nUseGap )
		RP->dPFootprint[t][n] = 1.0;		
	      else 
		      RP->dPFootprint[t][n] = 1.0 / nLen;  */ /* TEST 07/10/06 */ /* TEST !!!  10/12/04 */
	    }
	}
    }
  
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      nMotifWidth = MotifWidth(B, t);
      for( n = 0; n <= nEndPos - nMotifWidth + 1; n++ )
	{
	  if( PossibleStartPos( Pos, n + nStart, t, nLen, nSeq, B ) )
	    {
	      in_motif_prob( B, *P, n + nStart, t, FALSE, TRUE, 
			      &dMotifProb, &dBGProb );
	      RP->dProb[1][t][FORWARD][n] = dMotifProb / dBGProb;

	      if(IP->RevComplement)
		{
		  in_motif_prob( B, *P, n + nStart, t, TRUE, TRUE, 
				 &dMotifProb, &dBGProb );
		  RP->dProb[1][t][REVERSE][n] = dMotifProb / dBGProb;
		}

	      if( RP->bUsePosMatrix )
		{
		  RP->dProb[1][t][FORWARD][n] *= RP->dPosMatrix[0][t];
		  RP->dProb[1][t][REVERSE][n] *= RP->dPosMatrix[0][t];
		}

	      RP->prSum[1] += (RP->dProb[1][t][FORWARD][n] + RP->dProb[1][t][REVERSE][n])
		* RP->dEndSiteProb[t]; 
	      
	      if( RP->dProb[1][t][FORWARD][n] > dCutOff[t][FORWARD] )
		dCutOff[t][FORWARD] =  RP->dProb[1][t][FORWARD][n];
	      if( RP->dProb[1][t][REVERSE][n] > dCutOff[t][REVERSE] )
		dCutOff[t][REVERSE] =  RP->dProb[1][t][REVERSE][n];
	    }
	}
    }
    
#ifdef _DEBUG_
  memset( prf, 0, sizeof( prf ) );
  memset( prr, 0, sizeof( prr ) );
  memset( fpr, 0, sizeof( fpr ) );
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      for( j = 0; j < min( nLen, GRAPH_SIZE ) ; j++ )      /* DEBUG */
	{
	  prf[j] = max( prf[j], RP->dProb[1][t][FORWARD][j] );
	  prr[j] = max( prr[j],  RP->dProb[1][t][REVERSE][j]) ;
	  fpr[j] = max( fpr[j], RP->dPFootprint[t][j] );
	}
    }  
  if( nSeq == 0 )     /* DEBUG */
    nSeq = nSeq;      
#endif
  
  RP->prTotalSum = RP->prSum[1];

  minMotif = MinLenMotif( B );
  minMotifWidth = MotifWidth( B, minMotif );

  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      dCutOff[t][FORWARD] *= IP->dRCutoff;
      dCutOff[t][REVERSE] *= IP->dRCutoff;
    }

  if( RP->bUseFlatGap )
    nMax = RP->nFlatGapLen;
  else
    nMax = nLen;
  
  /* Now calcultate for other block sizes */
  for( k = 2; k <= RP->nMaxBlocks; k++ )
    {
      if( RP->bUseFlatGap )  /* BT 4/27/05 */
	dGapPr = RP->dFlatGapProb;
      else if( RP->nUseGap )
	dGapPr = 1.0;
      else
	/* dGapPr = 1.0 / (nLen - (k-1) * minMotifWidth + 1.0);  */
	dGapPr = 1.0; 

      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  motifWidth =  MotifWidth( B, t );
	  for( dir = FORWARD; 
	       dir <= REVERSE * ( IP->RevComplement ? 1 : 0); 
	       dir++ )
	    {
	      for(j = (k-1) * minMotifWidth;  
		  j <= nEndPos - motifWidth + 1; 
		  j++ )
		{
		  if( PossibleStartPos( Pos, j + nStart, t, nLen, nSeq, B ) )
		    {
		      sum = 0.0;
		      for( t1 = 0; t1 <  IP->nNumMotifTypes; t1++ )
			{
			  sum2 = 0.0;
			  nMotifWidtht1 = MotifWidth( B, t1 );
			  for( d = FORWARD; 
			       d <= REVERSE * (IP->RevComplement ? 1 : 0); 
			       d++ )
			    {
			      for( n = 0; n <= nMax; n++ )
				{
				  nPos = j - nMotifWidtht1 - n;
				  if( nPos >= 0 )
				    {
				      if( RP->nUseGap && ! RP->bUseFlatGap )
					sum2 += RP->dProb[k-1][t1][d][nPos] * RP->dPGap[t][t1][n];
				      else
					sum2 += RP->dProb[k-1][t1][d][nPos];
				    }
				  else
				    break;
				}
			    }
			  sum2 *= RP->dPTrans[t][0][t1][0]; /* we ignore the direction in the trans matrix */
			  sum += sum2;
			}

		      sum *= dGapPr * RP->dProb[1][t][dir][j];
		      RP->dProb[k][t][dir][j] = sum;
		      RP->prSum[k] += RP->dProb[k][t][dir][j] * RP->dEndSiteProb[t]; 
		    }
		}	      
	    }
	}
      
      
#ifdef _DEBUG_
      memset( prf, 0, sizeof( prf ) );
      memset( prr, 0, sizeof( prr ) );
      memset( fpr, 0, sizeof( fpr ) );
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  for( j = 0; j < min( nLen, GRAPH_SIZE ) ; j++ )      /* DEBUG */
	    {
	      prf[j] = max( prf[j], RP->dProb[k][t][FORWARD][j] );
	      prr[j] = max( prr[j], RP->dProb[k][t][REVERSE][j] );
	      fpr[j] = max( fpr[j], RP->dPFootprint[t][j] );
	    }
	}
      if( nSeq == 0 )     /* DEBUG */
	nSeq = nSeq;      
#endif

      RP->prTotalSum += RP->prSum[k];
    }

  FREEP( dCutOff, IP->nNumMotifTypes );
}


void DumpProbArray( Model B, int nSeq, int kmax, int iter, int seed )
{
  IPtype IP;
  RPType RP;
  ALType AP;
  int    i;
  int    j;
  int    k;
  double sum;
  char   fName[1024];
  FILE*  fpt;

  if( B->RP == NULL )
    return;
    
  IP = B->IP;
  RP = B->RP;
  AP = B->AP;

  strcpy( fName, "/cascv/ccmb/data/thompson/gibbs/output/dump/dump_" );
  sprintf( &fName[strlen(fName)], "%d_%d_%d", seed, nSeq, iter );

  fpt = fopen( fName, "w" );

  fprintf( fpt, "seed iter seq pos motif sites probf probr align\n" );
  
  for( k = 1; k <= kmax; k++ )
    {
      sum = 0;
      for( j = 0; j < IP->nNumMotifTypes; j++ )
	{
	  for( i = 0; i < RP->nSeqLen; i++ )
	    {
	      fprintf( fpt, "%d %d %d %d %d %d %g %g %g\n",  
		       seed, iter, nSeq, i, j, k, RP->dProb[k][j][FORWARD][i], 
		       RP->dProb[k][j][REVERSE][i], AP->dAlignCount[k][nSeq] );
	      sum +=  RP->dProb[k][j][FORWARD][i] + RP->dProb[k][j][REVERSE][i];
	    }
	}
    }    

  fclose( fpt );
}


void FreeAlignmentCounts( Model B )
{
  if( B->AP != NULL )
    {
      FREEP( B->AP->dAlignCount, B->AP->nMaxBlocks + 1 );
      FREEP( B->AP->dAlignCountSave, B->AP->nMaxBlocks + 1 );
      free( B->AP->dAlignWt );
      
      free( B->AP );
      B->AP = NULL;
    }
}


double BiCoef( int n, int k)
     /* Binomial Coefficient from Numerical Recipes p 215 */
{
  if( n == 0 )
    return 1.0;
  
  return floor( 0.5 + exp( ln_factrl( n ) - ln_factrl( k ) - ln_factrl( n - k )));
}


double LnBiCoef( int n, int k)
{
  if( n == 0 )
    return 0;
  
  return (ln_factrl( n ) - ln_factrl( k ) - ln_factrl( n - k ));
}


void CalcAlignProbFromFootprint( Model B )
{
  IPtype IP;
  RPType RP;
  int    seq;
  int    len;
  int    n;
  double sum = 0; 

#ifdef _DEBUG_
  double     prf[GRAPH_SIZE];   /* DEBUG */
  double     prr[GRAPH_SIZE];   /* DEBUG */
  double     prn[GRAPH_SIZE];   /* DEBUG */
#endif
 
  if( B->RP == NULL )
    return;
    
  IP = B->IP;
  RP = B->RP;
   
  for(seq = 0; seq < IP->nNumSequences; seq++) 
     {
       len = SequenceLength( B, seq );
       if( RP->dProbCons[seq][0] > 1.0 )
	 {
	   for( n = 0; n < len; n++ )
	     {
	       RP->dProbSite[seq][n] = DEFAULT_SITE_PROB;
	       RP->dProbNonSite[seq][n] = DEFAULT_NONSITE_PROB;
	       sum += RP->dProbSite[seq][n];
	     }
	 }
       else
	 {
	   for( n = 0; n < len; n++ )
	     {
	       RP->dProbSite[seq][n] = RP->dProbCons[seq][n] * RP->dSiteProbCons +
		 (1.0 - RP->dProbCons[seq][n]) * RP->dSiteProbNonCons; 
	       RP->dProbNonSite[seq][n] = DEFAULT_NONSITE_PROB;
	       sum += RP->dProbCons[seq][n];
	     }
	 }
              
#ifdef _DEBUG_
      memset( prf, 0, sizeof( prf ) );
      memset( prr, 0, sizeof( prr ) );
      memset( prn, 0, sizeof( prn ) );
      for( n = 0; n < min( len, GRAPH_SIZE ) ; n++ )      /* DEBUG */
	{
	  prf[n] = RP->dProbSite[seq][n];
	  prn[n] = RP->dProbNonSite[seq][n];
	  prr[n] = RP->dProbSite[seq][n] / RP->dProbNonSite[seq][n];
	}
      n = n;
#endif
       
     }
}


double BkgndProbFromComp( Model B, int nSeq, int nPos, int t, short bRevComp )
{
  IPtype  IP;
  BkgType BP;
  int     i;
  double  bgProb = 1.0;
  int     nStart;
  int     nEndPos;
  int     nLen;
  int     n;
  int     index;

  IP = B->IP;
  BP = B->BP;
  
  nLen = SequenceLength( B, nSeq );
  nEndPos = nLen - 1;
  nStart = SequenceStartPos( B, nSeq );

  for( i = 0; i < IP->nMotifLen[t]; i++ )
    {
      if( ! bRevComp )
	{
	  n = nStart + nPos + i;
	  index = CH2INT(n,B);
	  bgProb *= BP->dBkgndProb[nSeq][nPos +i][index];
	}
      else
	{
	  n = nStart + nPos + (IP->nMotifLen[t] - i - 1);
	  index = CH2INT(n, B);
	  bgProb *= BP->dBkgndProb[nSeq][n-nStart][index];
	}
    }
  return bgProb;
}


void CountAlignments( Model B,  PoSition **Pos  )    
{
  int  nSeq;

  for( nSeq = 0; nSeq < B->IP->nNumSequences; nSeq++ )
    {
#ifdef _DEBUG_
      printf("Counting alignments for sequence %d\n", nSeq);
#endif
      CountAlignments2( B,  Pos, nSeq  );
    }
}


void CountAlignments2( Model B,  PoSition **Pos, int nSeq  )
{
  IPtype     IP;
  RPType     RP;
  ALType     AP;
  int        n;
  int        nLen;
  int        nEndPos;
  int        t;
  int        j;
  int        t1;
  int        v;
  int        k;
  int        minMotif;
  int        maxMotif;
  int        nStart;
  double     ***algn;
  double     dGapPr;
  double     fpSiteProb;
  double     fpBGProb;
  double     sum;
  char*      residual;
  int        nMotifWidth;
  int        nMinMotifWidth;
  int        nMaxMotifWidth;
  int        nMaxGapWidth;
  double     *algnkt;
  double     *algnk1t1;
  double     *algn1t;
  int        nPossSite;
  int        p;
  int        pStart;
  
  residual=(*B->Seq->R);
 
  IP = B->IP;
  RP = B->RP;
  AP = B->AP;

  /* Calculate first for 1 block */
  nLen = SequenceLength( B, nSeq );
  nEndPos = nLen - 1;
  nStart = SequenceStartPos( B, nSeq );

  NEWPP( algn, RP->nMaxBlocks + 1, double );
  for( k = 0; k <= RP->nMaxBlocks; k++ )
    {
      NEWP( algn[k], IP->nNumMotifTypes, double );  
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	NEW( algn[k][t], nLen, double );
    }

  if( IP->is_defined[cl_T] )    /* BT 07/21/2000 */
    UpdateFootprintParameters( B, Pos, nSeq );

  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      algn1t = algn[1][t];
      nMotifWidth = MotifWidth( B, t );
      sum = 0.0;
      for( n = 0; n <= nEndPos - nMotifWidth + 1; n++ )
	{
	  RP->sitePos[nSeq][n][t].dSpacingProb = 0;
	  nPossSite = PossMotifStartPos( Pos, n + nStart, t, B );
	  if( nPossSite && IsPhyloSeq( B, nSeq ) )
	    {
	      pStart = SequenceStartPos( B, B->Phylo->phyloTreeStartSeq[B->Phylo->phyloIndex[nSeq]] );
	      for( p = 0; p < SeqInc( B, nSeq ); p++ )
		nPossSite &= PossMotifStartPos( Pos, n + pStart + p * nLen, t, B );
	    } 
	  if( nPossSite 
	      || residual[n + nStart] == 'x'        /* BT 03/21/2000 */
	      || residual[n + nStart] == 'X' )      /* X out prob - revisit */
	    {
	      algn1t[n]++;

	      if( IP->is_defined[cl_T] )    /* BT 11/17/99 */
		{
		  RP->sitePos[nSeq][n][t].dSpacingProb = 
		    SiteProbFromFootprint( B, n, t, nSeq, 
					   &fpSiteProb, &fpBGProb );
		  sum += RP->sitePos[nSeq][n][t].dSpacingProb;
		}
	      else if( RP->nUseSpacingProb )
		{
		  RP->sitePos[nSeq][n][t].dSpacingProb = RP->dPSpacing[nSeq][n];
		}

	      if(IP->RevComplement)
		{
		  algn1t[n]++;
		}
	    }
	}

      if( IP->is_defined[cl_T] )
	{
	  for( n = 0; n <= nEndPos - nMotifWidth + 1; n++ )
	    RP->sitePos[nSeq][n][t].dSpacingProb /= sum;	  
	}
      
      if( RP->nUseSpacingProb || IP->is_defined[cl_T] )
	{
	  RP->sitePos[nSeq][0][t].dPartialSum = RP->sitePos[nSeq][0][t].dSpacingProb;
	  for( n = 1; n < nLen; n++ )
	    {
	      RP->sitePos[nSeq][n][t].dPartialSum = 
		RP->sitePos[nSeq][n-1][t].dPartialSum + 
		RP->sitePos[nSeq][n][t].dSpacingProb;
	    }
	}
    }
    
  minMotif = MinLenMotif( B );
  nMinMotifWidth = MotifWidth( B, minMotif ); 
  maxMotif = MaxLenMotif( B );
  nMaxMotifWidth = MotifWidth( B, maxMotif ); 

  nMaxGapWidth = 0;
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      for( t1 = 0; t1 < IP->nNumMotifTypes; t1++ )
	{
	  if( RP->nMaxGapWidth[t1][t] > nMaxGapWidth )
	    nMaxGapWidth =  RP->nMaxGapWidth[t1][t];
	}
    }

  for( AP->dAlignCount[1][nSeq] = 0, j = 0; j < nLen; j++ )
    {
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{	      
	  AP->dAlignCount[1][nSeq] += algn[1][t][j];        
	}
    }
  
  /* Now calcultate for other block sizes */
  for( k = 2; k <= RP->nMaxBlocks; k++ )
    {
      if( RP->dProbCons[nSeq][0] > 1.0 )
	{
	  if( nLen - (k-1) *  nMinMotifWidth > 0 )
	    dGapPr = 1.0 / (nLen - (k-1) *  nMinMotifWidth + 1.0);
	  else
	    dGapPr = 0;
	}
      else
	dGapPr = 1.0;

      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  algnkt = algn[k][t];
	  algn1t = algn[1][t];
	  for(j = (k-1) *  nMinMotifWidth; 
	      j <= nEndPos - nMinMotifWidth + 1; 
	      j++ )
	    {
	      if( algn1t[j] )
		{
		  if( nMaxGapWidth == 0 )
		    {
		      if( dGapPr > 0 )
			{
			  for( t1 = 0; t1 < IP->nNumMotifTypes; t1++ )
			    {
			      algnk1t1 = algn[k-1][t1];
			      nMotifWidth = MotifWidth( B, t1 );
			      for( v = 0;
				   v <= j - nMotifWidth; v++ )  /* BT 1/26/2001 */
				{
				  algnkt[j] += algnk1t1[v];
#ifdef _DEBUG_
				  if( ! finite( algnkt[j] ) )
				    p_error("**** non finite result ****");
#endif	      	      
				}
			    }
			}
		    }
		  else
		    {
		      for( t1 = 0; t1 < IP->nNumMotifTypes; t1++ )
			{
			  algnk1t1 = algn[k-1][t1];
			  nMotifWidth = MotifWidth( B, t1 );
			  for( v = max( 0, j -  nMaxGapWidth - nMotifWidth);
			       v <= j - nMotifWidth; v++ )  /* BT 1/26/2001 */
			    {
			      algnkt[j] += algnk1t1[v];
			    }
			}
		    }
		  algnkt[j] *= algn1t[j];
		}
	    }
	}

      for( AP->dAlignCount[k][nSeq] = 0, j = 0; j < nLen; j++ )
	{
	  for( t = 0; t < IP->nNumMotifTypes; t++ )
	    {	      
	      AP->dAlignCount[k][nSeq] += algn[k][t][j];        	    
#ifdef _DEBUG_
	      if( ! finite( AP->dAlignCount[k][nSeq] ) )
		p_error("**** non finite result ****");
#endif	      	      
	    }
	}
    }

  FREEPP( algn, RP->nMaxBlocks + 1, IP->nNumMotifTypes );  
}


void CountAlignments2x( Model B,  PoSition **Pos, int nSeq  )
{
  IPtype     IP;
  RPType     RP;
  ALType     AP;
  int        n;
  int        nLen;
  int        nEndPos;
  int        t;
  int        dir;
  int        j;
  int        t1;
  int        v;
  int        k;
  int        minMotif;
  int        maxMotif;
  int        nStart;
  double     **algn;
  double     dGapPr;
  double     fpSiteProb;
  double     fpBGProb;
  double     sum;
  char*      residual;
  int        nMotifWidth;
  int        nMinMotifWidth;
  int        nMaxMotifWidth;
  int        nMaxGapWidth;
  int        nOffset;
  double     *algnk;
  double     *algnk1;
  int        nPossSite;
#ifdef _DEBUG_
  double     prf[GRAPH_SIZE];   /* DEBUG */
  double     prr[GRAPH_SIZE];
#endif
  
  residual=(*B->Seq->R);
 
  IP = B->IP;
  RP = B->RP;
  AP = B->AP;

  /* Calculate first for 1 block */
  nLen = SequenceLength( B, nSeq );
  nEndPos = nLen - 1;
  nStart = SequenceStartPos( B, nSeq );

  NEWP( algn, RP->nMaxBlocks + 1, double );
  for( k = 0; k <= RP->nMaxBlocks; k++ )
    NEW( algn[k], nLen, double );

  if( IP->is_defined[cl_T] )    /* BT 07/21/2000 */
    UpdateFootprintParameters( B, Pos, nSeq );
      
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      nMotifWidth = MotifWidth( B, t );
      sum = 0.0;
      for( n = 0; n <= nEndPos - nMotifWidth + 1; n++ )
	{
	  RP->sitePos[nSeq][n][t].dSpacingProb = 0;
	  nPossSite = PossMotifStartPos( Pos, n + nStart, t, B );
	  if( nPossSite && IsPhyloSeq( B, nSeq ) )
	    {
	      if( nSeq % 2 == 0 )
		nPossSite &= PossMotifStartPos( Pos, n + nStart + nLen, t, B );
	      else
		nPossSite &= PossMotifStartPos( Pos, n + nStart - nLen, t, B );
	    } 
	  if( nPossSite 
	      || residual[n + nStart] == 'x'        /* BT 03/21/2000 */
	      || residual[n + nStart] == 'X' )      /* X out prob - revisit */
	    {
	      algn[1][n]++;

	      if( IP->is_defined[cl_T] )    /* BT 11/17/99 */
		{
		  RP->sitePos[nSeq][n][t].dSpacingProb = 
		    SiteProbFromFootprint( B, n, t, nSeq, 
					   &fpSiteProb, &fpBGProb );
		  sum += RP->sitePos[nSeq][n][t].dSpacingProb;
		}
	      else if( RP->nUseSpacingProb )
		{
		  RP->sitePos[nSeq][n][t].dSpacingProb = RP->dPSpacing[nSeq][n];
		}

	      if(IP->RevComplement)
		{
		  algn[1][n]++;
		}
	    }
	}

      if( IP->is_defined[cl_T] )
	{
	  for( n = 0; n <= nEndPos - nMotifWidth + 1; n++ )
	    RP->sitePos[nSeq][n][t].dSpacingProb /= sum;	  
	}
      
      if( RP->nUseSpacingProb || IP->is_defined[cl_T] )
	{
	  RP->sitePos[nSeq][0][t].dPartialSum = RP->sitePos[nSeq][0][t].dSpacingProb;
	  for( n = 1; n < nLen; n++ )
	    {
	      RP->sitePos[nSeq][n][t].dPartialSum = 
		RP->sitePos[nSeq][n-1][t].dPartialSum + 
		RP->sitePos[nSeq][n][t].dSpacingProb;
	    }
	}

#ifdef _DEBUG_
      memset( prf, 0, sizeof( prf ) );
      memset( prr, 0, sizeof( prr ) );
      for( n = 0; n < min( nLen, GRAPH_SIZE ) ; n++ )      /* DEBUG */
	{
	  prf[n] = RP->sitePos[nSeq][n][t].dSpacingProb;
	  prr[n] = RP->dProbCons[nSeq][n];
	}
      
      if( nSeq == 16 )     /* DEBUG */
	nSeq = nSeq;            
#endif      
    }
    
  minMotif = MinLenMotif( B );
  nMinMotifWidth = MotifWidth( B, minMotif ); 
  maxMotif = MaxLenMotif( B );
  nMaxMotifWidth = MotifWidth( B, maxMotif ); 

  nMaxGapWidth = 0;
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      for( t1 = 0; t1 < IP->nNumMotifTypes; t1++ )
	{
	  if( RP->nMaxGapWidth[t1][t] > nMaxGapWidth )
	    nMaxGapWidth =  RP->nMaxGapWidth[t1][t];
	}
    }

  /*   if( ! RP->nUseSpacingProb  && ! IP->is_defined[cl_T]  )  
       { */   /* BT 2/14/2001 */
  for( AP->dAlignCount[1][nSeq] = 0, j = 0; j < nLen; j++ )
    AP->dAlignCount[1][nSeq] += algn[1][j];
  
  /* Now calcultate for other block sizes */
  for( k = 2; k <= RP->nMaxBlocks; k++ )
    {
      algnk = algn[k];
      algnk1 = algn[k-1];
      if( RP->dProbCons[nSeq][0] > 1.0 )
	{
	  if( nLen - (k-1) *  nMinMotifWidth > 0 )
	    dGapPr = 1.0 / (nLen - (k-1) *  nMinMotifWidth + 1.0);
	  else
	    dGapPr = 0;
	}
      else
	dGapPr = 1.0;
      
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  nMotifWidth =  MotifWidth( B, t );
	  /* nOffset = (k-2) *  nMotifWidth; */
	  nOffset = 0;  /* BT 2/13/2001 */
	  for( dir = FORWARD; 
	       dir <= REVERSE * ( IP->RevComplement ? 1 : 0); 
	       dir++ )
	    {
	      for(j = (k-1) *  nMinMotifWidth; 
		  j <= nEndPos -  nMotifWidth + 1; 
		  j++ )
		{
		  /* nPossSite = PossMotifStartPos( Pos, j + nStart, t, B ); 
		  if( nPossSite && IP->is_defined[cl_D] )
		    {
		      if( nSeq % 2 == 0 )
			nPossSite &= PossMotifStartPos( Pos, j + nStart + nLen, t, B );
		      else
			nPossSite &= PossMotifStartPos( Pos, j + nStart - nLen, t, B );
			} */

#ifdef _DEBUG_
		  if( nSeq == 30 && j == 1477 )
		    j = j;
#endif
		  
		  nPossSite = PossFragStartPos( Pos, j + nStart, t, B);
		  if( IsPhyloSeq( B, nSeq ) )
		    {
		      if( nSeq % 2 == 0 )
			nPossSite &= PossFragStartPos( Pos, j + nStart + nLen, t, B ); 
		      else
			nPossSite &= PossFragStartPos( Pos, j + nStart - nLen, t, B );
		    } 

		  if( nPossSite )
		    {
		      if( nMaxGapWidth == 0 )
			{
			  if( dGapPr > 0 )
			    {
			      for( v = 0; 
				   v <= j - nMotifWidth; v++ )  /* BT 1/26/2001 */
				/* v <= j - (k-1) * nMotifWidth; v++ ) */
				algnk[j] += algnk1[nOffset + v];
			    }
			}
		      else
			{
			  for( v = max( 0, j - nMaxGapWidth - 
					nMinMotifWidth);
			       /* nMotifWidth); */ /* BT 09/03/2002 */
			       v <= j - nMotifWidth; v++ )  /* BT 1/26/2001 */
			    /* v <= j - (k-1) *  nMotifWidth; v++ ) */
			    {
			      algnk[j] += algnk1[nOffset + v];				
			    }
			}
		    }
		}
	    }
	}
	  
      for( AP->dAlignCount[k][nSeq] = 0, j = 0; j < nLen; j++ )
	AP->dAlignCount[k][nSeq] += algn[k][j];  
      
    }
  /* } */

#ifdef _DEBUG_
  if( nSeq == 30 || nSeq == 31 )
    DumpAlignCnts( B, nSeq, algn );
#endif

  FREEP( algn, RP->nMaxBlocks + 1 );  
}


/* BT 04/06/2009 */
short PossMotifStartPos(PoSition **Pos, int n, int t, Model B)
{
   return Pos[t][n].nPossStartPos;
}


short PossMotifStartPos2(PoSition **Pos, int n, int t, Model B)
{
   int j;
   PoSition *Pos_t;  /* Pos_t=Pos[t] */

   Pos_t=Pos[t];
   if(!B->IP->is_defined[cl_F]) 
   {	
     for(j = 0; 
	 j < B->F->FragWidth[t] - 
	   (int) ((1.0 - B->IP->glOverlapParam) * B->IP->nMotifLen[t]) + 1; 
	 j++)
     {
       if(!Pos_t[n + j].nPossStartPos)
	 return FALSE;  
     }
   }
   else{ /* not fragmentation */
     if( !Pos_t[n].nPossStartPos)
       return FALSE;
   }

   return TRUE;
}


void SaveAlignmentCounts( Model B )    
{
  int     nSeq;
  ALType  AP;
  int     k;
  int     n;
  int     t;
  RPType  RP;
  IPtype  IP;

  AP = B->AP;
  RP = B->RP;
  IP = B->IP;

  for( nSeq = 0; nSeq < B->IP->nNumSequences; nSeq++ )
    {
      for( k = 0; k <= B->RP->nMaxBlocks; k++ )
	{
	  AP->dAlignCountSave[k][nSeq] = AP->dAlignCount[k][nSeq];
	}
      for( n = 0; n < SequenceLength( B, nSeq ); n++ )
	{
	  for( t = 0; t < IP->nNumMotifTypes; t++ )
	    {
	      RP->sitePos[nSeq][n][t].dSpacingProbSave = 
		RP->sitePos[nSeq][n][t].dSpacingProb;
	      RP->sitePos[nSeq][n][t].dPartialSumSave = 
		RP->sitePos[nSeq][n][t].dPartialSum;
	    }
	}
    }
}


void RestoreAlignmentCounts( Model B )    
{
  int     nSeq;
  ALType  AP;
  int     k;
  int     n;
  int     t;
  RPType  RP;
  IPtype  IP;

  AP = B->AP;
  RP = B->RP;
  IP = B->IP;

  for( nSeq = 0; nSeq < B->IP->nNumSequences; nSeq++ )
    {
      for( k = 0; k <= B->RP->nMaxBlocks; k++ )
	{
	  AP->dAlignCount[k][nSeq] = AP->dAlignCountSave[k][nSeq];
	}
      for( n = 0; n < SequenceLength( B, nSeq ); n++ )
	{
	  for( t = 0; t < IP->nNumMotifTypes; t++ )
	    {
	      RP->sitePos[nSeq][n][t].dSpacingProb = 
		RP->sitePos[nSeq][n][t].dSpacingProbSave;
	      RP->sitePos[nSeq][n][t].dPartialSum = 
		RP->sitePos[nSeq][n][t].dPartialSumSave;
	    }
	}
    }
}


/* PossibleStartPos - is this a possible start of a site */
/* nOffset is the length of the current sequence - only used for homologous pair */
short PossibleStartPos(PoSition **Pos, int n, int t, int nOffset, int seq, Model B)
{
  int p;
  
  if( ! IsPhyloSeq( B, seq ) )
    return PossFragStartPos( Pos, n, t, B);
  else
    {      
      for( p = 0; p < SeqInc( B, seq ); p++ )
	{
	  if( ! PossFragStartPos( Pos, n + p * nOffset, t, B) )
	    return FALSE;
	}
      return TRUE;
    }      

  /*  if( ! B->IP->is_defined[cl_D] )
    return PossMotifStartPos( Pos, n, t, B);
  else
  return (PossMotifStartPos( Pos, n, t, B) && PossMotifStartPos( Pos, n + nOffset, t, B)); */
}


void DumpAlignCnts( Model B, int nSeq, double **algn )
{
  IPtype IP;
  RPType RP;
  ALType AP;
  int    j;
  int    k;
  char   fName[1024];
  FILE*  fpt;

  if( B->RP == NULL )
    return;
    
  IP = B->IP;
  RP = B->RP;
  AP = B->AP;

  strcpy( fName, "/local/bioproj/people/thompson/Dump/dump_algn" );
  sprintf( &fName[strlen(fName)], "%d", nSeq );

  fpt = fopen( fName, "w" );

  fprintf( fpt, "k\tpos\tcount\n" );
  
  for( k = 0; k <= RP->nMaxBlocks; k++ )
    {
      for( j = 0; j < SequenceLength( B, nSeq ); j++ )	
	fprintf( fpt, "%d %d %g\n",  k, j, algn[k][j] );
    }    

  fclose( fpt );
}


char *GibbsTempnam(const char *dir, const char *pfx)
{
  return tempnam( dir, pfx );
}
