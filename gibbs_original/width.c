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
/* $Id: width.c,v 1.7 2007/11/19 17:46:40 Bill Exp $          */
/*                                                                        */
/* Author :       Bill Thompson 5/9/97                                    */
/*                                                                        */
/* Description :  This file contains the functions which are used to find */
/*                the best motif width. These routines are based on notes */
/*                Jun Liu 3/5/97                                          */
/*                                                                        */
/**************************************************************************/

#include "width.h"

#define LEFT    0
#define RIGHT   1

void CalcWidthProb( Model B, PoSition **Pos, Mlist M, long double ***dWidthProb, short sample );
short PossibleWidth( int w, int t, Model B, PoSition **Pos, Mlist M, 
		     int dir, int oldLen );
short PossibleMotifWidth( int w, int t, Model B, PoSition **Pos, MlistEl *pMotif, 
			  int dir, int oldLen );
short PossibleOverlap( Model B, int curPos, int w, int dir, Mlist M,  
		       MlistEl *pCurMotif, int oldLen );
short CheckPossibleOverlap( Model B, int curPos, int w, int dir, Mlist M,  
			    MlistEl *pCurMotif, int oldLen );
short CheckPossibleOverlapRight( Model B, int curPos, int w, int dir, Mlist M,  
				MlistEl *pCurMotif, MlistEl *pMotif, int oldLen,
				int t );
short CheckPossibleMotifRight( Model B, int curPos, int w, int dir, Mlist M,  
				MlistEl *pCurMotif, MlistEl *pMotif, int oldLen,
				int t );
short CheckPossibleOverlapLeft( Model B, int curPos, int w, int dir, Mlist M,  
			       MlistEl *pCurMotif, MlistEl *pMotif, int oldLen,
			       int t );
short CheckPossibleMotifLeft( Model B, int curPos, int w, int dir, Mlist M,  
				MlistEl *pCurMotif, MlistEl *pMotif, int oldLen,
				int t );
void ResetLengths( int w, int t, Model B, PoSition **Pos, 
		   Mlist M, int dir, int oldLen, short sample );
void CopyMotifs( Model B, Mlist M, Mlist M2 );
void CopyTheMotif( MlistEl *pMotif, int t, Model B, Mlist M, Mlist M2 );
void ReleaseMotifs( Model B, Mlist M2 );
void ResetGoodPosition( int t, Model B, PoSition **Pos, Mlist M );
void DumpMotifs( Model B, Mlist M );
int GetMaxMotifWidth( Model B );
void DeleteMotif(Model B,  MlistEl *pCurMotif, Mlist M, int t); 		/* BT 11/26/97 */

#ifndef __LONGDOUBLE128
extern long double expl (long double __x) ; 
extern long double powl (long double __x, long double __y) ; 
#endif

void ResizeMotif( Model B, PoSition **Pos, Mlist M, short sample )
{
  int         i;
  int         t;
  IPtype      IP;
  Ctype       C;
  long double ***dWidthProb;
#ifdef _DEBUG_		  
  double      dPrevProb;
#endif

  IP= B->IP;
  C = B->C;

  NEWPP( dWidthProb, IP->nNumMotifTypes, long double );
  for( t= 0; t < IP->nNumMotifTypes; t++ )
    {
      NEWP( dWidthProb[t], 2, long double );
      for( i = 0; i < 2; i++ )
	{
	  NEW( dWidthProb[t][i], IP->nMaxMotifLen[t] + 1, long double );
	}
    }
  
#ifdef _DEBUG_		  
  dPrevProb = CalcMapProb(B, B->IP->is_defined[cl_R]);
#endif

  CalcWidthProb( B, Pos, M, dWidthProb, sample ); 

#ifdef _DEBUG_		  
  fprintf( stdout, "Prev Map: %g New Map: %g\n", 
	   dPrevProb, CalcMapProb(B, B->IP->is_defined[cl_R]) );
#endif

  for( t= 0; t < IP->nNumMotifTypes; t++ )
    {
      for( i = 0; i < 2; i++ )
	free( dWidthProb[t][i] );
      free( dWidthProb[t] );
    }
  free( dWidthProb );
}


void CalcWidthProb( Model B, PoSition **Pos, Mlist M, long double ***dWidthProb, short sample )
{
  IPtype       IP;
  Ctype        C;
  int          t;
  int          i;
  int          j;
  int          w; 
  int          col;
  MlistEl      *pMotif;
  long double  dProb;
  int          oldLen;
  long double  dTot;
  long double  rnd;
  short        bFoundit;
  int          maxLen;
  int          minLen;
  int          dir;
  double       dProb2;
  int          len; 
  double       pw;
  double       ***dMap;
  Mlist        M2;
  int          nMotifCnt[2];
  int          ***nMotifCntSave;
#ifdef _DEBUG_
  long double  dPr[256];
  long double  dPrl[256];
  long double  dPrr[256];
#endif

  CheckCounts( B );
  IP= B->IP;
  C = B->C;

  NEWP( M2, IP->nNumMotifTypes, Mlist_struct );
  for( i = 0; i < IP->nNumMotifTypes; i++ )
    {
      NEW( M2[i], 1, Mlist_struct );
      M2[i]->Motifs = NULL;
      M2[i]->nNumMotifs = 0;
    }

  NEWPP( dMap, IP->nNumMotifTypes, double );
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      NEWP( dMap[t], 2, double );
      for( dir = LEFT; dir <= RIGHT; dir++ )
        NEW( dMap[t][dir], IP->nMaxMotifLen[t] + 1, double );
    }

  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      oldLen = IP->nMotifLen[t];
      nMotifCnt[FORWARD] = IP->nNumMotifs[t][FORWARD];
      nMotifCnt[REVERSE] = IP->nNumMotifs[t][REVERSE];
      bFoundit = FALSE;
      minLen = IP->nMinMotifLen[t];
      maxLen = IP->nMaxMotifLen[t];

      dir = RIGHT;

      NEWPP( nMotifCntSave, 2, int );
      for( i = 0; i < 2; i++ )
	{
	  NEWP( nMotifCntSave[i], 2, int );
	  for( j = 0; j < 2; j++ )
	    {
	      NEW( nMotifCntSave[i][j], maxLen + 1, int );
	    }
	}
  
      pw = log( IP->nMaxMotifLen[t] - IP->nMinMotifLen[t] + 1 );    /* p(w) */

      pMotif = M[t]->Motifs;
      while( pMotif != NULL )               /* remove old motif counts */
	{
	  adjust_counts( B, DELETE, pMotif->pos, t, pMotif->RevComp );
	  pMotif = pMotif->next;
	}
      
      for( w = minLen; w <= maxLen; w++ )
	{
	  for( dir = LEFT; dir <= RIGHT; dir++ )
	    {
	      IP->nMotifLen[t] = oldLen;
	      IP->nNumMotifs[t][FORWARD] = nMotifCnt[FORWARD];
	      IP->nNumMotifs[t][REVERSE]= nMotifCnt[REVERSE];

	      CopyMotifs(  B, M, M2 );
	      if( PossibleWidth( w, t, B, Pos, M2, dir, oldLen ) )
		{
		  nMotifCntSave[dir][FORWARD][w] = IP->nNumMotifs[t][FORWARD];
		  nMotifCntSave[dir][REVERSE][w] = IP->nNumMotifs[t][REVERSE];

		  IP->nMotifLen[t] = w;
		  col = w;
		  pMotif = M2[t]->Motifs;
		  while( pMotif != NULL )
		    {
		      if( (dir == RIGHT && !pMotif->RevComp) ||
			  (dir == LEFT && pMotif->RevComp))
			adjust_counts( B, ADD, pMotif->pos, t, pMotif->RevComp );
		      else
			{
			  if( IP->is_defined[cl_R] )
			    adjust_counts( B, ADD, pMotif->pos - ((w - oldLen) / 2), t, 
					   pMotif->RevComp );
			  else			    
			    adjust_counts( B, ADD, pMotif->pos - (w - oldLen), t, 
					   pMotif->RevComp );
			}
		      pMotif = pMotif->next;
		    }

#ifdef _DEBUG_		  
		  if( ! CheckCounts( B ) )
		    w = w;   /*****/ /* DEBUG */
#endif

		  set_indicator_vector(Pos, B->Seq, B );
		  SetPossibleSites( B, t );
		  B->C->dtot_sites[t] = IP->nPossSites[t];
		  update_posterior_prob(B);

		  if( IP->is_defined[cl_E] )
		    CountAlignments( B,  Pos );
		  
		  /*		  nSiteSamp = IP->site_samp; */
		  /*		  B->IP->site_samp = 1; */
		  dMap[t][dir][w] = CalcMapProb(B, B->IP->is_defined[cl_R]);
		  /*		  B->IP->site_samp = nSiteSamp; */
		  
		  dProb2 = dMap[t][dir][w] - pw;
		  if( IP->site_samp )
		    {
		      for( i = 0; i < IP->nNumSequences; i++ )
			{
			  len = SequenceLength( B, i );
			  dProb2 -= log(len - w);
			}  
		    }
		  /* else
		    {
		      dProb3 = ln_factrl(  NUMMOTIFS( IP->nNumMotifs[t] ) );
		      for( i = 0; i <  NUMMOTIFS( IP->nNumMotifs[t] ); i++ )
			{
			  if( IP->nPossSites[t] - i * w <= 0 )
			    break;
			  dProb3 -= log( IP->nPossSites[t] - i * w );
			}
		      dProb2 += dProb3; 
		      } */


		  dWidthProb[t][dir][w] = expl( dProb2 );

		  pMotif = M2[t]->Motifs;
		  while( pMotif != NULL )
		    {
		      if( (dir == RIGHT && !pMotif->RevComp) ||
			  (dir == LEFT && pMotif->RevComp))
			adjust_counts( B, DELETE, pMotif->pos, t, pMotif->RevComp );
		      else
			{
			  if( IP->is_defined[cl_R] )
			    adjust_counts( B, DELETE, pMotif->pos -((w - oldLen) / 2), t, 
					   pMotif->RevComp );
			  else			    
			    adjust_counts( B, DELETE, pMotif->pos - (w - oldLen), t, 
					   pMotif->RevComp );
			}
		      pMotif = pMotif->next;
		    }
		}
	      ReleaseMotifs( B, M2 );
	    }
	}

      IP->nMotifLen[t] = oldLen;
      IP->nNumMotifs[t][FORWARD] = nMotifCnt[FORWARD];
      IP->nNumMotifs[t][REVERSE] = nMotifCnt[REVERSE];

      for( dTot = 0.0, w = minLen; w <= maxLen; w++ )
	{
	  dTot += dWidthProb[t][LEFT][w] + dWidthProb[t][RIGHT][w]; 
#ifdef _DEBUG_		  
	  dPrl[w] = dWidthProb[t][LEFT][w];    /* DEBUG */
	  dPrr[w] = dWidthProb[t][RIGHT][w];    /* DEBUG */
	  dPr[w] = dWidthProb[t][LEFT][w] + dWidthProb[t][RIGHT][w];    /* DEBUG */
#endif
	}

      if( dTot > 0.0 )
	{ 
	  dProb = 0;
	  rnd = drand();
	  bFoundit = FALSE;
	  for( w = minLen; w <= maxLen && !bFoundit; w++ )
	    /*	  for( bFoundit = FALSE, dir = LEFT; dir <= RIGHT && !bFoundit; dir++ ) */
	    {
	      /*	      for( w = minLen; w <= maxLen && !bFoundit; w++ ) */
	      for( dir = LEFT; dir <= RIGHT && !bFoundit; dir++ ) 
		{
		  dProb += dWidthProb[t][dir][w] / dTot;
		  if( dProb > rnd )
		    {
		      ResetLengths( w, t, B, Pos, M, dir, oldLen, sample );
		      bFoundit = TRUE;
		      break;
		    }
	        }
	    }
	  if( ! bFoundit )
	    ResetLengths( oldLen, t, B, Pos, M, RIGHT, oldLen, sample );
	}
      else
	ResetLengths( oldLen, t, B, Pos, M, RIGHT, oldLen, sample );

      if( IP->is_defined[cl_E] )
	CountAlignments( B,  Pos );

#ifdef _DEBUG_		  
      if( ! CheckCounts( B ) )
	w = w;   /*****/ /* DEBUG */
#endif
      
      FREEPP( nMotifCntSave, 2, 2 );
    }
  
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    SetPossibleSites( B, t );

  update_posterior_prob( B );

#ifdef _DEBUG_
  dProb2 = CalcMapProb(B, B->IP->is_defined[cl_R]);
#endif
  
  for( t= 0; t < IP->nNumMotifTypes; t++ )
    {
      for( dir = LEFT; dir <= RIGHT; dir++ )
	free( dMap[t][dir] );
      free( dMap[t] );
    }
  free( dMap );

  FREEP( M2, IP->nNumMotifTypes );

#ifdef _DEBUG_		  
  CheckCounts (B ); /* DEBUG */
#endif
}


short PossibleWidth( int w, int t, Model B, PoSition **Pos, Mlist M, int dir, int oldLen )
{
  IPtype  IP;
  MlistEl *pMotif;
  MlistEl *pCurMotif;

  IP = B->IP;
  pMotif = M[t]->Motifs;
  while( pMotif != NULL )
    {
      if( (dir == RIGHT && !pMotif->RevComp) ||
	  (dir == LEFT && pMotif->RevComp))
	{
	  if( pMotif->pos + w >= IP->nSeqLen || 
	      /* Pos[t][pMotif->pos].nSeq != Pos[t][pMotif->pos + w - 1].nSeq  ||
		 PossibleOverlap( B, pMotif->pos, w, dir, M,  pMotif, oldLen ) */
	      Pos[t][pMotif->pos].nSeq != Pos[t][pMotif->pos + w - 1].nSeq ||
	      pMotif->pos + w - 1 > pMotif->right )
	    {
	      return FALSE;
	      if( IP->site_samp )
		return FALSE;
	      B->IP->nNumMotifs[t][pMotif->RevComp ? REVERSE : FORWARD]--;
	      pCurMotif = pMotif;
	      pMotif = pMotif->next;
	      DeleteMotif( B, pCurMotif, M, t);
	    }
	  else
	    {
	      if( CheckPossibleOverlap( B, pMotif->pos, w, dir, M,  pMotif, oldLen ) )
		return FALSE;
	      pMotif = pMotif->next;
	    }
	}
      else
	{
	  if( IP->is_defined[cl_R] )
	    {
	      if( pMotif->pos - ((w - oldLen) / 2) < 0 || 
		  pMotif->pos - ((w - oldLen) / 2) >= IP->nSeqLen || 
		  Pos[t][pMotif->pos].nSeq != Pos[t][pMotif->pos - ((w - oldLen) / 2)].nSeq ||
		  Pos[t][pMotif->pos].nSeq != 
		    Pos[t][pMotif->pos - ((w - oldLen) / 2) + w - 1].nSeq ||
		  pMotif->pos - ((w - oldLen) / 2) + w - 1 > pMotif->right ||
		  pMotif->pos - ((w - oldLen) / 2) < pMotif->left )   /* BT 10/17/97 */
		{
		  return FALSE;
		  if( IP->site_samp )
		    return FALSE;
		  B->IP->nNumMotifs[t][pMotif->RevComp ? REVERSE : FORWARD]--;
		  pCurMotif = pMotif;
		  pMotif = pMotif->next;
		  DeleteMotif( B, pCurMotif, M, t);
		}
	      else
		{
		  if( CheckPossibleOverlap( B, pMotif->pos, w, dir, M,  pMotif, oldLen ) )
		    return FALSE;
		  pMotif = pMotif->next;
		}
	    }
	  else
	    {
	      if( pMotif->pos - (w - oldLen) < 0 || 
		  pMotif->pos - (w - oldLen) >= IP->nSeqLen || 
		  Pos[t][pMotif->pos].nSeq != Pos[t][pMotif->pos - (w - oldLen)].nSeq ||
		  pMotif->pos - (w - oldLen) < pMotif->left )   /* BT 10/17/97 */
		{
		  return FALSE;
		  if( IP->site_samp )
		    return FALSE;
		  B->IP->nNumMotifs[t][pMotif->RevComp ? REVERSE : FORWARD]--;
		  pCurMotif = pMotif;
		  pMotif = pMotif->next;
		  DeleteMotif( B, pCurMotif, M, t);
		}
	      else
		{
		  if( CheckPossibleOverlap( B, pMotif->pos, w, dir, M,  pMotif, oldLen ) )
		    return FALSE;
		  pMotif = pMotif->next;
		}
	    }
	}
    }

  return TRUE;
}


short CheckPossibleOverlap( Model B, int curPos, int w, int dir, Mlist M,  
			    MlistEl *pCurMotif, int oldLen )
{
  IPtype   IP;
  int      t;
  MlistEl  *pMotif;
  MlistEl  *pNext;
  short    overlap = FALSE;

  IP = B->IP;
  for( t = 0; t < IP->nNumMotifTypes && ! overlap; t++ )
    {
      pMotif = M[t]->Motifs;
      while( pMotif != NULL && ! overlap )
	{
	  pNext = pMotif->next;
	  if( pMotif->pos != pCurMotif->pos )
	    {
	      if( pCurMotif->seq_num == pMotif->seq_num )
		{
		  if( dir == RIGHT )
		    {
		      if( t == pCurMotif->Mtype )
			overlap = CheckPossibleOverlapRight( B, curPos, w, dir, M,  
							     pCurMotif, pMotif, oldLen, t );
		      else
			overlap = CheckPossibleMotifRight( B, curPos, w, dir, M,  
							   pCurMotif, pMotif, oldLen, t );
		    }
		  else
		    {
		      if( t == pCurMotif->Mtype )
			overlap = CheckPossibleOverlapLeft( B, curPos, w, dir, M,  
							    pCurMotif, pMotif, oldLen, t );
		      else
			overlap = CheckPossibleMotifLeft( B, curPos, w, dir, M,  
							  pCurMotif, pMotif, oldLen, t );
		    }
		}
	    }
	  pMotif = pNext;
	}
    }

  return overlap;
}


short CheckPossibleOverlapRight( Model B, int curPos, int w, int dir, Mlist M,  
				MlistEl *pCurMotif, MlistEl *pMotif, int oldLen,
				int t )
{
  if( pCurMotif->pos < pMotif->pos )
    {
      if( !pCurMotif->RevComp )
	{
	  if( !pMotif->RevComp )
	    {
	      if( pCurMotif->pos + w - 1 >= pMotif->pos )
		{
		  if( B->IP->site_samp )
		    return TRUE;
		  DeleteMotif( B, pMotif, M, t);
		  B->IP->nNumMotifs[t][FORWARD]--;
		}
	    }
	  else
	    {
	      if( pCurMotif->pos + w - 1 >= pMotif->pos - (w - oldLen) )
		{
		  if( B->IP->site_samp )
		    return TRUE;
		  DeleteMotif( B, pMotif, M, t);
		  B->IP->nNumMotifs[t][REVERSE]--;
		}
	    }
	}
      else
	{
	  if( pMotif->RevComp )
	    {
	      if( pCurMotif->pos - (w - oldLen) + w - 1 >= pMotif->pos - (w - oldLen) )
		{
		  if( B->IP->site_samp )
		    return TRUE;
		  DeleteMotif( B, pMotif, M, t);
		  B->IP->nNumMotifs[t][REVERSE]--;
		}
	    }
	}
    }
  else
    {
      if( !pCurMotif->RevComp )
	{
	  if( !pMotif->RevComp )
	    {
	      if( pMotif->pos + w - 1 >= pCurMotif->pos )
		{
		  if( B->IP->site_samp )
		    return TRUE;
		  DeleteMotif( B, pMotif, M, t);
		  B->IP->nNumMotifs[t][FORWARD]--;
		}
	    }
	}
      else
	{
	  if( !pMotif->RevComp )
	    {
	      if( pMotif->pos + w - 1 >= pCurMotif->pos - (w - oldLen) )
		{
		  if( B->IP->site_samp )
		    return TRUE;
		  DeleteMotif( B, pMotif, M, t);
		  B->IP->nNumMotifs[t][FORWARD]--;
		}
	    }
	  else
	    {
	      if( pMotif->pos - (w - oldLen) + w - 1 >= pCurMotif->pos - (w - oldLen) )
		{
		  if( B->IP->site_samp )
		    return TRUE;
		  DeleteMotif( B, pMotif, M, t);
		  B->IP->nNumMotifs[t][REVERSE]--;
		}
	    }
	}			  
    }

  return FALSE;
}

short CheckPossibleMotifRight( Model B, int curPos, int w, int dir, Mlist M,  
			       MlistEl *pCurMotif, MlistEl *pMotif, int oldLen,
			       int t )
{
  /* Check for overlap with other motif types */
  if( pCurMotif->pos < pMotif->pos )
    {
      if( !pCurMotif->RevComp )
	{
	  if( pCurMotif->pos + w - 1 >= pMotif->pos )
	    return TRUE;
	}
    }
  else
    {
      if( pCurMotif->RevComp )
	{
	  if( pMotif->pos + B->IP->nMotifLen[pMotif->Mtype] - 1 >= pCurMotif->pos - (w - oldLen) )
	    return TRUE;
	}
    }

  return FALSE;
}


short CheckPossibleOverlapLeft( Model B, int curPos, int w, int dir, Mlist M,  
				MlistEl *pCurMotif, MlistEl *pMotif, int oldLen,
				int t )
{
  if( pCurMotif->pos < pMotif->pos )
    {
      if( !pCurMotif->RevComp )
	{
	  if( !pMotif->RevComp )
	    {
	      if( B->IP->is_defined[cl_R] )
		{
		  if( pCurMotif->pos + oldLen - 1 >= pMotif->pos - ((w - oldLen) / 2) )
		    {
		      if( B->IP->site_samp )
			return TRUE;
		      DeleteMotif( B, pMotif, M, t);
		      B->IP->nNumMotifs[t][FORWARD]--;
		    }
		}
	      else
		{
		  if( pCurMotif->pos + oldLen - 1 >= pMotif->pos - (w - oldLen) )
		    {
		      if( B->IP->site_samp )
			return TRUE;
		      DeleteMotif( B, pMotif, M, t);
		      B->IP->nNumMotifs[t][FORWARD]--;
		    }
		}
	    }
	}
      else
	{
	  if( ! pMotif->RevComp )
	    {
	      if( pCurMotif->pos + w - 1 >= pMotif->pos - (w - oldLen) )
		{
		  if( B->IP->site_samp )
		    return TRUE;
		  DeleteMotif( B, pMotif, M, t);
		  B->IP->nNumMotifs[t][FORWARD]--;
		}
	    }
	  else
	    {
	      if( pCurMotif->pos + w - 1 >= pMotif->pos )
		{
		  if( B->IP->site_samp )
		    return TRUE;
		  DeleteMotif( B, pMotif, M, t);
		  B->IP->nNumMotifs[t][REVERSE]--;
		}
	    }
	}
    }
  else
    {
      if( !pCurMotif->RevComp )
	{
	  if( !pMotif->RevComp )
	    {
	      if( B->IP->is_defined[cl_R] )
		{
		  if( pMotif->pos + w - 1 >= pCurMotif->pos - ((w - oldLen) / 2) )
		    {
		      if( B->IP->site_samp )
			return TRUE;
		      DeleteMotif( B, pMotif, M, t);
		      B->IP->nNumMotifs[t][FORWARD]--;
		    }
		}
	      else
		{
		  if( pMotif->pos + w - 1 >= pCurMotif->pos - (w - oldLen) )
		    {
		      if( B->IP->site_samp )
			return TRUE;
		      DeleteMotif( B, pMotif, M, t);
		      B->IP->nNumMotifs[t][FORWARD]--;
		    }
		}
	    }
	  else
	    {
	      if( pMotif->pos + (w - oldLen) + w - 1 >=  pCurMotif->pos - (w - oldLen) )
		{
		  if( B->IP->site_samp )
		    return TRUE;
		  DeleteMotif( B, pMotif, M, t);
		  B->IP->nNumMotifs[t][REVERSE]--;
		}
	    }
	}
      else
	{
	  if( pMotif->RevComp )
	    {
	      if( pMotif->pos + w - 1 >= pCurMotif->pos )
		{
		  if( B->IP->site_samp )
		    return TRUE;
		  DeleteMotif( B, pMotif, M, t);
		  B->IP->nNumMotifs[t][REVERSE]--;
		}
	    }
	}			  
    }		
      
  return FALSE;
}


short CheckPossibleMotifLeft( Model B, int curPos, int w, int dir, Mlist M,  
				MlistEl *pCurMotif, MlistEl *pMotif, int oldLen,
				int t )
{
  /* Check for overlap with other motif types */
  if( pCurMotif->pos < pMotif->pos )
    {
      if( pCurMotif->RevComp )
	{
	  if( pCurMotif->pos + w - 1 >= pMotif->pos )
	    return TRUE;
	}
    }
  else
    {
      if( !pCurMotif->RevComp )
	{
	  if( B->IP->is_defined[cl_R] )
	    {	      
	      if( pMotif->pos + B->IP->nMotifLen[pMotif->Mtype] - 1 >= 
		  pCurMotif->pos - ((w - oldLen) / 2) )
		return TRUE;
	    }
	  else
	    {
	      if( pMotif->pos + B->IP->nMotifLen[pMotif->Mtype] - 1 >= 
		  pCurMotif->pos - (w - oldLen) )
		return TRUE;
	    }
	}
    }		
      
  return FALSE;
}


short PossibleMotifWidth( int w, int t, Model B, PoSition **Pos, MlistEl *pMotif, int dir, int oldLen )
{
  IPtype  IP;

  IP = B->IP;

  if( dir == RIGHT )
    {
      if( pMotif->pos + w >= IP->nSeqLen || Pos[t][pMotif->pos].nSeq != Pos[t][pMotif->pos + w].nSeq|| 
	  PossibleOverlapsPrevMotif( B, pMotif->pos, t,  Pos ) )
	return FALSE;
    }
  else
    {
      if( pMotif->pos - (w - oldLen) < 0 || 
	  Pos[t][pMotif->pos].nSeq != Pos[t][pMotif->pos - (w - oldLen)].nSeq || 
	  pMotif->pos - (w - oldLen) >= IP->nSeqLen || 
	  PossibleOverlapsPrevMotif( B, pMotif->pos - (w - oldLen), t,  Pos ) )
	return FALSE;
    }

  return TRUE;
}


void ResetLengths( int w, int t, Model B, PoSition **Pos, 
		   Mlist M, int dir, int oldLen, short sample )
{
  IPtype   IP;
  MlistEl  *pMotif;
  MlistEl  *pOldMotif;
  MlistEl  *next;
  int      pos;
  int      t1;

  IP = B->IP;

  pMotif = M[t]->Motifs;
  while( pMotif != NULL )
    {
      not_in_motif( Pos, pMotif->pos, B, t );
      pMotif = pMotif->next;
    }

   if( ! PossibleWidth( w, t, B, Pos, M, dir, oldLen ) )
       p_internal_error( "Invalid Length");  

  IP->nMotifLen[t] = w;
  M[t]->nMotifLen = w;

  set_indicator_vector(Pos, B->Seq, B );
  for( t1 = 0; t1 < IP->nNumMotifTypes; t1++ )
    {
      if( t1 != t )
	{            
	  pMotif = M[t1]->Motifs;
	  while( pMotif != NULL )
	    {
	      set_in_motif( Pos, pMotif->pos, B, t1, pMotif->RevComp );
	      pMotif = pMotif->next;
	    }
	}
    }
  
  IP->nNumMotifs[t][FORWARD] = 0;
  IP->nNumMotifs[t][REVERSE]= 0;
      
  pMotif = M[t]->Motifs;
  pOldMotif = M[t]->Motifs;
  M[t]->Motifs = NULL;
  M[t]->nNumMotifs = 0;
  while( pMotif != NULL )
    {
      if( (dir == RIGHT && !pMotif->RevComp) ||
	  (dir == LEFT && pMotif->RevComp))
	pos = pMotif->pos;
      else
	{
	  if( IP->is_defined[cl_R] )
	    pos = pMotif->pos - ((w - oldLen) / 2);
	  else
	    pos = pMotif->pos - (w - oldLen);
	}
      adjust_counts( B, ADD, pos, t, pMotif->RevComp );
      set_in_motif( Pos, pos, B, t, pMotif->RevComp );
      add_motif( B->IP, B->Seq, pos, M, t,  pMotif->RevComp );
      if( pMotif->RevComp )
	IP->nNumMotifs[t][REVERSE]++;
      else
	IP->nNumMotifs[t][FORWARD]++;
#ifdef _DEBUG_		  
      CheckCounts( B );   /* Debug */
#endif      
      pMotif = pMotif->next;
    } 
  
  pMotif = pOldMotif;
  while( pMotif != NULL )
    {
      next = pMotif->next;
      free( pMotif );
      pMotif = next;
    }

#ifdef _DEBUG_		  
  CheckCounts( B );   /* Debug */
#endif

  if( sample )
    IP->nWidthCnts[t][dir][w]++;

  if( ! B->IP->is_defined[cl_Z] )
    fprintf( stdout, "\nMotif %d width set to %d dir = %d sites = %d\n", 
	     t, w, dir, NUMMOTIFS( IP->nNumMotifs[t] ) );
}


/*
void CopyMotifs( int t, Model B, PoSition **Pos, Mlist M, Mlist M2 )
{
  IPtype   IP;
  MlistEl  *pMotif;

  M2[t]->Motifs = NULL;
  M2[t]->nNumMotifs = 0;
  
  pMotif = M[t]->Motifs;
  while( pMotif != NULL )
    {
      add_motif( B->Seq, pMotif->pos, M2, t,  pMotif->RevComp );      
      pMotif = pMotif->next;
    }
}
*/

void CopyMotifs( Model B, Mlist M, Mlist M2 )
{
  int t;

  for( t = 0; t < B->IP->nNumMotifTypes; t++ )
    {
      M2[t]->Motifs = NULL;
      M2[t]->nNumMotifs = 0;
      CopyTheMotif( M[t]->Motifs, t, B, M, M2 );
    }
}


void CopyTheMotif( MlistEl *pMotif, int t, Model B, Mlist M, Mlist M2 )
{
  if( pMotif == NULL )
    return;

  CopyTheMotif( pMotif->next, t, B, M, M2 );
  add_motif( B->IP, B->Seq, pMotif->pos, M2, t,  pMotif->RevComp );      
}


void ReleaseMotifs( Model B, Mlist M2 )
{
  MlistEl  *pMotif;
  MlistEl  *next;
  int      t;

  for( t = 0; t < B->IP->nNumMotifTypes; t++ )
    {
      pMotif = M2[t]->Motifs;
      while( pMotif != NULL )
	{
	  next = pMotif->next;
	  free( pMotif );
	  pMotif = next;
	}
      
      M2[t]->Motifs = NULL;
      M2[t]->nNumMotifs = 0;
    }
}


void ResetGoodPosition( int t, Model B, PoSition **Pos, Mlist M )
{
  IPtype   IP;
  MlistEl  *pMotif;
  MlistEl  *pOldMotif;
  MlistEl  *next;
  int      goodPos[18] = {17, 17, 76, 63, 50, 7, 42, 39, 9, 14, 29, 41, 48, 71, 17, 53, 5, 78};
  int      i;

  IP = B->IP;
  pOldMotif = M[t]->Motifs;
  M[t]->Motifs = NULL;
  M[t]->nNumMotifs = 0;
  for( i = 0; i < 18; i++ )
    {
      set_in_motif( Pos, goodPos[i], B, t, FALSE );
      add_motif( B->IP, B->Seq, goodPos[i], M, t,  FALSE );      
    }

  pMotif = pOldMotif;
  while( pMotif != NULL )
    {
      next = pMotif->next;
      free( pMotif );
      pMotif = next;
    }
}


void DumpCounts( Model B, int t, FILE *fpt )
{
  IPtype   IP;
  int      i;
  int      j;
  double   tot;
  double   tot2;
  double   sum1 = 0;
  double   sum2 = 0;

  fprintf( fpt, "\n Motif %d\n", t);
  IP = B->IP;
  for( i = 0; i < IP->nAlphaLen; i++ )
    {
      for( tot = 0, tot2 = 0, j = 0; j <= IP->nMotifLen[t]; j++ )
	{
	  fprintf( fpt, "%6.2f ", B->C->wCounts[t][j][i] );
	  if( j > 0 )
	    tot += B->C->wCounts[t][j][i];
	  else 
	    tot2 += B->C->wCounts[t][j][i];
	}
      fprintf( fpt, "%f\n", tot );
      sum1 += tot;
      sum2 += tot2;
      if( tot2 == 0 )
	tot2 = tot2;
    }
  fprintf( fpt, "sum1 = %f\n", sum1 );
  fprintf( fpt, "sum2 = %f\n", sum2 );
  fprintf( fpt, "Motifs: %3d %3d\n", B->IP->nNumMotifs[t][0],  B->IP->nNumMotifs[t][1] );
}


void DumpMotifs( Model B, Mlist M )
{
  IPtype   IP;
  int      t;
  MlistEl  *pMotif;

  IP = B->IP;
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      pMotif = M[t]->Motifs;
      while( pMotif != NULL )
	{
	  fprintf( stdout, "%4d %4d, %4d %4d %4d %4d\n", pMotif->Mtype, pMotif->seq_num, pMotif->pos, pMotif->left, pMotif->right, pMotif->RevComp );
	  pMotif = pMotif->next;
	}
    }
}


void DumpPositions( int t, Model B, PoSition **Pos )
{
  IPtype   IP;
  PoSition *pos_t;
  int      i;

  pos_t = Pos[t];
  IP = B->IP;
  fprintf( stdout, "\nt = %d\n", t );
  
  for( i = 0; i < IP->nSeqLen; i++ )
    {
      if( pos_t[i].nMotifStartPos )
	  fprintf( stdout, "start = %d ", i );
      if( pos_t[i].nInMotif )
	  fprintf( stdout, "in = %d ", i );
    }
  fprintf( stdout, "\n" );
}


int GetMaxMotifWidth( Model B )
{
  int      i;
  int      len;
  int      t;
  int      minLen = INT_MAX;
  IPtype   IP;

  IP = B->IP;

  for( i = 0; i < IP->nNumSequences; i++ )
    {
      if( i == 0 )
	len = (*B->Seq->nvEndLocs)[i] - 1;
      else
	len = (*B->Seq->nvEndLocs)[i] - (*B->Seq->nvEndLocs)[i - 1] - 1;
      if( len < minLen )
	minLen = len;
    }  

  for( t = 0; t < IP->nNumMotifTypes; t++ )
    if( IP->nMaxMotifLen[t] < minLen )
      minLen = IP->nMaxMotifLen[t];

  return minLen;
}


double factrl( int n )
/* Return n! as a double. Adapted from Numerical Recipes p214 */
{
  static int ntop = 4;
  static double a[33] = {1.0, 1.0, 2.0, 6.0, 24.0};
  int j;

  if( n < 0 )
    p_internal_error( "Negative factorial in factrl." );

  if( n > 32 )
    return exp( ln_gamma( n + 1.0 ) );

  while( ntop < n )
    {
      j = ntop++;
      a[ntop] = a[j] * ntop;
    }

  return a[n];
}


double ln_factrl( int n )
{
  static double a[101];

  if( n < 0 )
    p_internal_error( "Negative factorial in ln_factrl." );

  if( n <= 1 )
    return 0.0;

  if( n <= 32 )
    return (a[n] ? a[n] : (a[n] = log( factrl( n ) )));
  else if( n <= 100 )
    return (a[n] ? a[n] : (a[n] = ln_gamma( n + 1.0 )));
  else
    return ln_gamma( n + 1.0 );
}


short CheckCounts( Model B )
{
  IPtype   IP;
  int      i;
  int      j;
  double   tot;
  int      t;

  IP = B->IP;
  for( t = 0; t < IP->nNumMotifTypes; t++ ) 
    {
      for( j = 1; j <= IP->nMotifLen[t]; j++ )
	{
	  for( tot = 0, i = 0; i < IP->nAlphaLen; i++ )
	    tot += B->C->fCounts[t][j][i];

	  if( tot != B->IP->nNumMotifs[t][0] + B->IP->nNumMotifs[t][1] )
	    {
	      p_internal_error("Width Sampler - Invalid Count" );
	      return FALSE;
	    }
	}
    }
  return TRUE;
}


void DeleteMotif(Model B,  MlistEl *pCurMotif, Mlist M, int t) 		/* BT 11/26/97 */
{
  delete_motif( B, pCurMotif->pos, M, t);
}



