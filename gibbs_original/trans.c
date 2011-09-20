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
#include "trans.h"


void InitializeTransMatrix( Model B )
     /* Note - we are not using FORWARD/REVERSE in the trans matrix. They are just carried along */
     /* for possible future updates */
{
  RPType RP;
  IPtype IP;
  int    t;
  int    t1;
  
  RP = B->RP;
  IP = B->IP;
  
  /* nNumMotifs = 0;
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      nNumMotifs += NUMMOTIFS( IP->nNumMotifs[t] );
    }  
  nNumMotifs -= IP->nNumSequences;
  nNumMotifs = max( 0, nNumMotifs ); */ /* BT 01/21/03 */

  RP->dTransTotalPriorCnts = 0.0;
  RP->dTotalPriorEndCnts = 0.0;
  
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      for( t1 = 0; t1 < IP->nNumMotifTypes; t1++ )
	{
	  RP->dPTrans[t][FORWARD][t1][FORWARD] = RP->dInitialPTrans[t][FORWARD][t1][FORWARD];
	  RP->dPTrans[t][FORWARD][t1][REVERSE] = RP->dInitialPTrans[t][FORWARD][t1][REVERSE];
	  RP->dPTrans[t][REVERSE][t1][FORWARD] = RP->dInitialPTrans[t][REVERSE][t1][FORWARD];
	  RP->dPTrans[t][REVERSE][t1][REVERSE] = RP->dInitialPTrans[t][REVERSE][t1][REVERSE]; 
	  
	  /*	  RP->dPriorTransCnts[t][t1] = RP->dPTrans[t][FORWARD][t1][FORWARD] * nNumMotifs; */
	  RP->dPriorTransCnts[t][t1] = RP->dPTrans[t][FORWARD][t1][FORWARD] *
	    NUMMOTIFS( IP->nNumMotifs[t] );    /* BT 01/21/03 */
	  RP->dPriorTransCnts[t][t1] = max( RP->dPriorTransCnts[t][t1], MIN_PSEUDOCNT );
	  RP->dTransTotalPriorCnts += RP->dPriorTransCnts[t][t1];
	}
      RP->dEndSiteProb[t] = RP->dInitialEndSiteProb[t];
      RP->dPriorEndCnts[t] = max( RP->dEndSiteProb[t] * RP->dPriorSeq, MIN_PSEUDOCNT ); /* BT 1/21/03 */
      RP->dTotalPriorEndCnts += RP->dPriorEndCnts[t];

      RP->dBeginSiteProb[t] = RP->dInitialBeginSiteProb[t];
      RP->dPriorBeginCnts[t] = max( RP->dBeginSiteProb[t] * RP->dPriorSeq, MIN_PSEUDOCNT ); /* BT 1/21/03 */
      RP->dTotalPriorBeginCnts += RP->dPriorBeginCnts[t];
    }
}


void UpdateTransMatrix( Model B )
{
  RPType RP;
  IPtype IP;
  int    t;
  int    t1;
  int    nNumMotifs;
  int    nSeq;
  int    nPrevType;
  int    nLen;
  int    j;
  double tr;
  double sum;
  double weight;
  double dTotalTransCnts = 0;
  double denom = 0;  
  
  RP = B->RP;
  IP = B->IP;

  nNumMotifs = 0;
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      nNumMotifs += NUMMOTIFS( IP->nNumMotifs[t] );
    }
  
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      for( t1 = 0; t1 < IP->nNumMotifTypes; t1++ )
	{
	  RP->dTransCnts[t1][t] = 0;
	}
      RP->dEndCnts[t] = 0;
      RP->dBeginCnts[t] = 0;
    }
  RP->dPosTotalCnts = 0;
  
  for( nSeq = 0; nSeq < IP->nNumSequences; nSeq++ )
    {
      nPrevType = -1;
      nLen = SequenceLength( B, nSeq );
      for( j = 0; j < nLen; j++ )
	{
	  for( t = 0; t < B->IP->nNumMotifTypes; t++ ) 
	    {
	      if( RP->sitePos[nSeq][j][t].nMotifStart )
		{
		  if( nPrevType != -1 )
		    {
		      RP->dTransCnts[t][nPrevType]++;
		      dTotalTransCnts++;
		    }
		  else
		    RP->dBeginCnts[t]++;
		  nPrevType = t;
		}
	    }
	}
      if( nPrevType != -1 )
	{
	  RP->dEndCnts[nPrevType]++;      
	  RP->dPosTotalCnts++;
	}
    }
  RP->dTransTotalCnts = dTotalTransCnts;

  weight = RP->dTransWt / (1.0 - RP->dTransWt); 
  
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      for( t1 = 0; t1 < IP->nNumMotifTypes; t1++ )  /* BT 4/19/2002 */
	{
	  denom += RP->dTransCnts[t][t1] + weight * RP->dPriorTransCnts[t][t1];
	}

      sum = 0;
      for( t1 = 0; t1 < IP->nNumMotifTypes; t1++ )
	{
	  /*	  denom = (dTotalTransCnts + weight * RP->dTransTotalPriorCnts ); 
	  if( denom > 0 )
	    tr = (RP->dTransCnts[t][t1] + weight * RP->dPriorTransCnts[t][t1]) /
	      (dTotalTransCnts + weight * RP->dTransTotalPriorCnts );
	  else
	  tr = 0; */

	  if( denom > 0 )  /* BT 4/19/2002 */
	    tr = (RP->dTransCnts[t][t1] + weight * RP->dPriorTransCnts[t][t1]) / denom;
	  else
	    tr = 0;
	  
	  RP->dPTrans[t][FORWARD][t1][FORWARD] = tr;
	  RP->dPTrans[t][FORWARD][t1][REVERSE] = tr; 
	  RP->dPTrans[t][REVERSE][t1][FORWARD] = tr;
	  RP->dPTrans[t][REVERSE][t1][REVERSE] = tr;

	  if( IP->RevComplement )
	    {
	      /* sum += 4 * tr; */
	      sum += tr;   /* BT 01/28/03 */
	    }
	  else
	    {
	      sum += tr;
	    }	  
	}

      /* Normalize */ /* BT 12/03/2002 */
        if( sum > 0 ) 
	{
	  for( t1 = 0; t1 < IP->nNumMotifTypes; t1++ )
	    {
	      RP->dPTrans[t][FORWARD][t1][FORWARD] /= sum;
	      RP->dPTrans[t][FORWARD][t1][REVERSE] /= sum;
	      RP->dPTrans[t][REVERSE][t1][FORWARD] /= sum;
	      RP->dPTrans[t][REVERSE][t1][REVERSE] /= sum;
	    } 
	} 
    }  

  /* Normalize */ /* BT 12/03/2002 */
  /*  for( t1 = 0; t1 < IP->nNumMotifTypes; t1++ )
    { 
      sum = 0;
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  sum += RP->dInitialPTrans[t][FORWARD][t1][FORWARD];
	}

      if( IP->RevComplement ) 
	sum *= 4;
	  
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  RP->dPTrans[t][FORWARD][t1][FORWARD] /= sum;
	  RP->dPTrans[t][FORWARD][t1][REVERSE] /= sum;
	  RP->dPTrans[t][REVERSE][t1][FORWARD] /= sum;
	  RP->dPTrans[t][REVERSE][t1][REVERSE] /= sum;
	}
	} */
  
  if( RP->bSymTrans )
    {
      for( t = 0; t < IP->nNumMotifTypes - 1; t++ )
	{
	  for( t1 = t+1; t1 < IP->nNumMotifTypes; t1++ )
	    {
	      tr = (RP->dPTrans[t][FORWARD][t1][FORWARD] + 
		    RP->dPTrans[t1][FORWARD][t][FORWARD]) / 2.0;
	      RP->dPTrans[t][FORWARD][t1][FORWARD] = tr;
	      RP->dPTrans[t1][FORWARD][t][FORWARD] = tr;
	      RP->dPTrans[t][FORWARD][t1][REVERSE] = tr;
	      RP->dPTrans[t1][FORWARD][t][REVERSE] = tr;
	      RP->dPTrans[t][REVERSE][t1][FORWARD] = tr;
	      RP->dPTrans[t1][REVERSE][t][FORWARD] = tr;
	      RP->dPTrans[t][REVERSE][t1][REVERSE] = tr;
	      RP->dPTrans[t1][REVERSE][t][REVERSE] = tr;
	    } 
	}	  
    }

  sum = 0;
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      tr = (RP->dEndCnts[t] + weight * RP->dPriorEndCnts[t]) /
        	(IP->nNumSequences + weight * RP->dTotalPriorEndCnts);
      RP->dEndSiteProb[t] = tr;
      sum += tr;
    }

  if( sum > 0 )
    {
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  RP->dEndSiteProb[t] /= sum;
	}
    }

  sum = 0;
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      tr = (RP->dBeginCnts[t] + weight * RP->dPriorBeginCnts[t]) /
        	(IP->nNumSequences + weight * RP->dTotalPriorBeginCnts);
      RP->dBeginSiteProb[t] = tr;
      sum += tr;
    }

  if( sum > 0 )
    {
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  RP->dBeginSiteProb[t] /= sum;
	}
    }
  
#ifdef _DEBUG_
  DumpTrans( B, stdout );   /* DEBUG */
#endif
}

void DumpTrans( Model B, FILE *fpt )
{
  RPType RP;
  IPtype IP;
  int    t;
  int    t1;
  int    nNumMotifs;
  
  RP = B->RP;
  IP = B->IP;

  nNumMotifs = 0;
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      nNumMotifs += NUMMOTIFS( IP->nNumMotifs[t] );
    }
  
  fprintf( fpt, "Transition Counts %d Motifs\n", nNumMotifs );
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      for( t1 = 0; t1 < IP->nNumMotifTypes; t1++ )
	{
	  fprintf( fpt, "%5.1f ", RP->dTransCnts[t][t1] );
	}
      fprintf( fpt, "\n" );
    }
  fprintf( fpt, "\n" );  

  fprintf( fpt, "Transition Matrix\n" );
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      for( t1 = 0; t1 < IP->nNumMotifTypes; t1++ )
	{
	  fprintf( fpt, "%6.3f ", RP->dPTrans[t][FORWARD][t1][FORWARD] );
	}
      fprintf( fpt, "\n" );
    }
  fprintf( fpt, "\n" );  

  fprintf( fpt, "End Counts %5.1f seqs\n",  RP->dPriorSeq );
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      fprintf( fpt, "%5.1f ", RP->dEndCnts[t] );
    }
  fprintf( fpt, "\n" );  

  fprintf( fpt, "End Probs\n" );
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      fprintf( fpt, "%6.3f ", RP->dEndSiteProb[t] );
    }
  fprintf( fpt, "\n" );  

  fprintf( fpt, "Begin Counts %5.1f seqs\n",  RP->dPriorSeq );
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      fprintf( fpt, "%5.1f ", RP->dBeginCnts[t] );
    }
  fprintf( fpt, "\n" );  

  fprintf( fpt, "Begin Probs\n" );
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      fprintf( fpt, "%6.3f ", RP->dBeginSiteProb[t] );
    }
  fprintf( fpt, "\n" );  

  fflush( fpt );
}
