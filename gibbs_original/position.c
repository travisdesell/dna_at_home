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
#include "position.h"

void DumpPosMatrix( Model B );


void InitializePosMatrix( Model B )
{
  RPType RP;
  IPtype IP;
  int    t;
  int    i;
  
  RP = B->RP;
  IP = B->IP;
  
  RP->dPosTotalPriorCnts = 0.0;
  
  for( i = 0; i <  RP->nModuleSites; i++ )
    {
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  RP->dPosMatrix[i][t] = RP->dInitialPosMatrix[i][t];
	  RP->dPriorPosMatrixCnts[i][t] = RP->dPosMatrix[i][t] *  RP->dPriorSeq;
	  RP->dPosTotalPriorCnts += RP->dPriorPosMatrixCnts[i][t];
	}
    }
}


void UpdatePosMatrix( Model B )
{
  RPType RP;
  IPtype IP;
  int    t;
  int    i;
  int    nSeq;
  int    nPos;
  int    nLen;
  int    j;
  double tr;
  double sum;
  double weight;
  double dTotalPosCnts = 0;
  double denom;
  char   *msg;
  
  RP = B->RP;
  IP = B->IP;

  for( i = 0; i < RP->nModuleSites; i++ )
    {
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  RP->dPosMatrixCnts[i][t] = 0;
	}
    }
  
  for( nSeq = 0; nSeq < IP->nNumSequences; nSeq++ )
    {
      nPos = 0;
      nLen = SequenceLength( B, nSeq );
      for( j = 0; j < nLen; j++ )
	{
	  for( t = 0; t < B->IP->nNumMotifTypes; t++ ) 
	    {
	      if( RP->sitePos[nSeq][j][t].nMotifStart )
		{
		  if( nPos >= RP->nModuleSites )
		    {
		      NEW( msg, 512, char );
		      sprintf( msg, "UpdatePosMatrix: module should have %d sites. %d found.",
			       RP->nModuleSites, nPos + 1 );
		      p_internal_error( msg );
		    }
		  RP->dPosMatrixCnts[nPos][t]++;
		  dTotalPosCnts++;
		  nPos++;
		}
	    }
	}
    }

  weight = RP->dPosWt / (1.0 - RP->dPosWt); 
  denom = (dTotalPosCnts + weight * RP->dPosTotalPriorCnts );
  
  for( i = 0; i < RP->nModuleSites; i++ )
    {
      sum = 0;
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  if( denom > 0 )
	    tr = (RP->dPosMatrixCnts[i][t] + weight * RP->dPriorPosMatrixCnts[i][t]) / denom;
	  else
	    tr = 0;
	  
	  RP->dPosMatrix[i][t] = tr;
	  sum += tr;
	}

      /* Normalize */
      if( sum > 0 )
	{
	  for( t = 0; t < IP->nNumMotifTypes; t++ )
	    {
	      RP->dPosMatrix[i][t] /= sum;
	    } 
	}
    }  
  
#ifdef _DEBUG_
  DumpPosMatrix( B );   /* DEBUG */
#endif  
}

void DumpPosMatrix( Model B )
{
  RPType RP;
  IPtype IP;
  int    t;
  int    i;
  
  RP = B->RP;
  IP = B->IP;
  
  fprintf( stdout, "Position Matrix Counts\n");
  for( i = 0; i < RP->nModuleSites; i++ )
    {
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  fprintf( stdout, "%5.1f ", RP->dPosMatrixCnts[i][t] );
	}
      fprintf( stdout, "\n" );
    }
  fprintf( stdout, "\n" );  

  fprintf( stdout, "Position Matrix\n" );
  for( i = 0; i < RP->nModuleSites; i++ )
    {
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  fprintf( stdout, "%6.3f ", RP->dPosMatrix[i][t] );
	}
      fprintf( stdout, "\n" );
    }
  fprintf( stdout, "\n" );  

  fflush( stdout );
}
