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
/*  $Id: fragment.c,v 1.8 2009/04/23 18:43:53 Bill Exp $         */
/*                                                                        */
/*  Author:  Eric Eric C. Rouchka, 1996.7                                 */
/*                Jun Zhu, 1996.8                                         */
/*  Fragmentation Routines -- Programs originally written by A. Neuwald   */
/*	modified by Judy Kilday.					  */
/*									  */
/*                               					  */
/**************************************************************************/

#include "fragment.h"

typedef struct
{
  int    count;
  double ratio;
} RatioStruct;

#define RATIO_ALLOC_SIZE 256
#define MAX_SAMPLES 1000
#define LMAX_SAMPLES log( (double) MAX_SAMPLES )

char FragmentSample( Model B, int typ, Mlist M );
char FragmentEnum( Model B, int typ, Mlist M  );
int ChooseColumns( Model B, int typ, int width, int offset, int *lgIndex );
void SetFragmentation( Model B, int typ, int offset, int *c, int *lgIndex );
int PossibleFragmentation( Model B, int typ, int mid, int offset, 
			   int firstPal, int lastPal, int firstRep, int lastRep, int *c, int *lgIndex );
int FragmentOverlap( Model B, int t, Mlist M );
int FirstFragPos( Model B, int typ );
int LastFragPos( Model B, int typ );
int FirstPalPos( Model B, int typ, int mid );
int LastPalPos( Model B, int typ, int mid );
int FirstRepPos( Model B, int typ, int mid );
int LastRepPos( Model B, int typ, int mid );
int intcompare( const void *p1, const void *p2 );


int     Frag( Model B, int typ, Mlist M  );
double  CalcWeight( Model B, int t, int lemon, int new_col, int *Mark );
double  CalcRatio( Model B, int t, int lemon, int new_col, int *Mark );
int	LemonFModel( Model B, int typ , int *Mark );
double	RatioFModel(int new,int new_col, int old, int old_col, Model B, int typ);
double	mv_col_weight(int c, int w, int n, int lemon, Model B, int typ);
int     CheckRepeat(int lemon, int pos, Model B, int typ, int *Mark);
int	MvColumnFModel(int lemon, int pos, Model B, int typ, int *Mark);
int	GetColPos(Model B, int typ, int colnum);
char	move_column_gibbs(Model B, int lemon, int typ, int *Mark, Mlist M );
int     FragOverlap( Model B, int t, int lemon, int new_col, int *Mark,  Mlist M );
void    FragmentFromMiddle( Model B, int typ, int *Mark, Mlist M );
int     FragmentFromMiddleLeftHalf( Model B, int typ, int *Mark, Mlist M );
int     FragmentFromMiddleRightHalf( Model B, int typ, int *Mark, Mlist M );
double  bico(int N,int k);
double	lnfact(int n);
double  lnbico(register int N, register int k);

#ifndef __LONGDOUBLE128
extern long double expl (long double __x) ; 
extern long double powl (long double __x, long double __y) ; 
#endif

/*========================================================================*/
/*  Allows for fragementation of a motif with a fixed number of columns	  */
/*  and allowing gaps between residues. The maximum possible width is 	  */
/*  set with the -M option                                                */
/*========================================================================*/


char Fragment( Model B, int typ, Mlist M  )
{
  if( B->IP->is_defined[cl_frag] )
    return Frag( B, typ, M );
  else if( B->IP->is_defined[cl_R] || B->IP->is_defined[cl_I] || B->IP->is_defined[cl_J] )
    return FragmentEnum( B, typ, M );
  else
    return FragmentSample( B, typ, M );
}


char FragmentSample( Model B, int typ, Mlist M  )
{
  IPtype      IP;
  Ftype       F;
  Ctype       C;
  int         i;
  int         j;
  int         *lgIndex;
  long double *lg;
  long double sum;
  long double rld;
  long double denom;
  double      r;
  double      width;
  double      **countSave;
  int         *fragSave;
  double      fragWidthSave;
  double      motifMap;
  double      bkgndMap;
  double      fragMap;
  double      currentMap;
  PhyloTree   *phyloSave;
  int         bkgndSave;
  int         phyloCalcSave;
  double      prob;
  double      cSum;
  double      pSum;
  
  IP = B->IP;
  F = B->F;
  C = B->C;    
  
  width = LastFragPos( B, typ ) - FirstFragPos( B, typ ) + 1;
  if( width == IP->nMotifLen[typ] )
    return FALSE;

  NEWP( countSave, IP->nMotifLen[typ] + 1, double );
  for( i = 0; i <= IP->nMotifLen[typ]; i++ )
    NEW( countSave[i], IP->nAlphaLen, double );

  for( i = 0; i <= IP->nMotifLen[typ]; i++ )
    {
      for( j = 0; j < IP->nAlphaLen; j++ )
	countSave[i][j] = C->wCounts[typ][i][j];
    }

  NEW( fragSave, F->nMaxLen[typ], int );
  for( i = 0; i < F->nMaxLen[typ]; i++ )
    fragSave[i] = F->nColMask[typ][i];
  fragWidthSave = F->FragWidth[typ];

  phyloSave = B->Phylo->phyloTree;
  phyloCalcSave = B->Phylo->bCalcPhylo;
  bkgndSave = IP->is_defined[cl_B];
  B->Phylo->phyloTree = NULL;
  B->Phylo->bCalcPhylo = FALSE;
  IP->is_defined[cl_B] = 0;

  motifMap = CalcMotifMap( B, typ, B->IP->is_defined[cl_R] );
  bkgndMap = CalcBkgndMap( B, B->IP->is_defined[cl_R] );
  fragMap = CalcMotifFragMap( B, typ,  B->IP->is_defined[cl_R] );	  
  currentMap = motifMap + bkgndMap + fragMap;  

  NEW( lgIndex, IP->nMotifLen[typ], int );
  NEW( lg, F->nMaxLen[typ], long double );

  for( i = 0; i < F->nMaxLen[typ]; i++ )
    {
      if( F->nColMask[typ][i] == COL_ON || F->nColMask[typ][i] == COL_OFF )
	{
	  for( cSum = 0, pSum = 0, j = 0; j < IP->nAlphaLen; j++ )
	    {
	      lg[i] += lgamma( F->nvFragCnts[typ][i][j] + C->dPseudoCounts[typ][BG][j] ) -
		lgamma( C->dPseudoCounts[typ][BG][j] );
	      cSum += F->nvFragCnts[typ][i][j];
	      pSum += C->dPseudoCounts[typ][BG][j];
	    }
	  lg[i] += lgamma( pSum ) - lgamma( cSum + pSum );
	  lg[i] = expl( lg[i] );
	}
    }

  for( i = 0; i < IP->nMotifLen[typ]; i++ )
    {
      lgIndex[i] = -1;
      for( denom = 0.0, j = 0; j < F->nMaxLen[typ]; j++ )
	denom += lg[j];
      sum = 0;
      rld = ((long double) drand()) * denom;
      for( j = 0; j < F->nMaxLen[typ]; j++ )
	{
	  sum += lg[j];
	  if( sum > rld )
	    {
	      lgIndex[i] = j;
	      lg[j] = 0;
	      break;
	    }
	}
    }
  qsort( lgIndex, IP->nMotifLen[typ], sizeof( int ), intcompare );

  for( i = 0; i < F->nMaxLen[typ]; i++ )
    {
      if( F->nColMask[typ][i] == COL_ON )
	{
	  F->nColMask[typ][i] = COL_OFF;
	  for( j = 0; j < IP->nAlphaLen; j++ )
	    {
	      C->wCounts[typ][BG][j] += F->nvFragCnts[typ][i][j];
	    }
	}
    }
  
  for( i = 0; i < IP->nMotifLen[typ]; i++ )
    {
      F->nColMask[typ][lgIndex[i]] = COL_ON;
      for( j = 0; j < B->IP->nAlphaLen; j++ )
	{
	  C->wCounts[typ][i+1][j] = F->nvFragCnts[typ][lgIndex[i]][j];
	  C->wCounts[typ][BG][j] -= F->nvFragCnts[typ][lgIndex[i]][j];
	}      
    }

  F->FragWidth[typ] = LastCol( F, typ ) - FirstCol( F, typ ) + 1;
  
  motifMap = CalcMotifMap( B, typ, B->IP->is_defined[cl_R] );
  bkgndMap = CalcBkgndMap( B, B->IP->is_defined[cl_R] );
  fragMap = CalcMotifFragMap( B, typ,  B->IP->is_defined[cl_R] ); 
  r = motifMap + bkgndMap + fragMap; 
  
  prob = exp( r - currentMap );
  if( FragmentOverlap( B, typ, M ) || prob < drand() )
    {
      for( i = 0; i < F->nMaxLen[typ]; i++ )
	F->nColMask[typ][i] = fragSave[i];
      F->FragWidth[typ] = fragWidthSave;
    }

  for( i = 0; i <= B->IP->nMotifLen[typ]; i++ )
    {
      for( j = 0; j < B->IP->nAlphaLen; j++ )
	C->wCounts[typ][i][j] = countSave[i][j];
    }
  FREEP( countSave, IP->nMotifLen[typ] + 1 );

  free( fragSave );
  B->Phylo->phyloTree = phyloSave;
  B->Phylo->bCalcPhylo = phyloCalcSave;
  IP->is_defined[cl_B] = bkgndSave;

  free( lg );
  free( lgIndex );

  return TRUE;
}


/* This routine generates all possible fragmentations using algorithm 7.2.1.3T */
/* from Knuth Art of Computer Programming Vol 4 preprint                       */
char FragmentEnum( Model B, int typ, Mlist M  )
{
  int         *c;
  int         j;
  int         retVal = FALSE;
  IPtype      IP;
  Ftype       F;
  Ctype       C;
  int         mid;
  int         firstPal;
  int         lastPal;
  int         firstRep;
  int         lastRep;
  int         n = 0;
  int         count = 0;
  double      denom = 0.0;
  double      sum;
  int         selected = -1;
  int         i;
  double      motifMap;
  double      bkgndMap;
  double      fragMap;
  double      currentMap;
  double      r;
  PhyloTree   *phyloSave;
  int         phyloCalcSave;
  int         bkgndSave;
  RatioStruct *ratio;
  double      **countSave;
  int         width;
  int         offset;
  int         x;
  int         *lgIndex;
  int         *fragSave;
  int         fragWidthSave;
  int         lgwidth;
  
  IP = B->IP;
  F = B->F;
  C = B->C;

  offset = FirstFragPos( B, typ );
  lgwidth = LastFragPos( B, typ ) - offset + 1;
  
  if( lgwidth == IP->nMotifLen[typ] )
    return FALSE;

  NEW( lgIndex, lgwidth, int );
  width = ChooseColumns( B, typ, lgwidth, offset, lgIndex );

  NEWP( countSave, IP->nMotifLen[typ] + 1, double );
  for( i = 0; i <= IP->nMotifLen[typ]; i++ )
    NEW( countSave[i], IP->nAlphaLen, double );

  for( i = 0; i <= IP->nMotifLen[typ]; i++ )
    {
      for( j = 0; j < IP->nAlphaLen; j++ )
	countSave[i][j] = C->wCounts[typ][i][j];
    }

  NEW( fragSave, F->nMaxLen[typ], int );
  for( i = 0; i < F->nMaxLen[typ]; i++ )
    fragSave[i] = F->nColMask[typ][i];
  fragWidthSave = F->FragWidth[typ];

  phyloSave = B->Phylo->phyloTree;
  phyloCalcSave = B->Phylo->bCalcPhylo;
  bkgndSave = IP->is_defined[cl_B];
  B->Phylo->phyloTree = NULL;
  IP->is_defined[cl_B] = 0;

  motifMap = CalcMotifMap( B, typ, B->IP->is_defined[cl_R] );
  bkgndMap = CalcBkgndMap( B, B->IP->is_defined[cl_R] );
  fragMap = CalcMotifFragMap( B, typ,  B->IP->is_defined[cl_R] );	  
  currentMap = motifMap + bkgndMap + fragMap;  

  mid = IP->nMotifLen[typ] / 2 + IP->nMotifLen[typ] % 2;
  firstPal = FirstPalPos( B, typ, mid );
  lastPal = LastPalPos( B, typ, mid );
  firstRep = FirstRepPos( B, typ, mid );
  lastRep = LastRepPos( B, typ, mid );

  NEW( c, IP->nMotifLen[typ] + 3, int );
  NEW( ratio, RATIO_ALLOC_SIZE, RatioStruct );

  for( j = 1; j <= IP->nMotifLen[typ]; j++ )
    c[j] = j - 1;

  c[IP->nMotifLen[typ]+1] = width;
  c[IP->nMotifLen[typ]+2] = 0;

  denom = 0.0;
  j = IP->nMotifLen[typ];
  while( TRUE )
    {
      if( PossibleFragmentation( B, typ, mid, offset, firstPal, lastPal, firstRep, lastRep, c, lgIndex ) )
	{
	  SetFragmentation( B, typ, offset, c, lgIndex );
	  if( ! FragmentOverlap( B, typ, M ) )
	    {
	      motifMap = CalcMotifMap( B, typ, B->IP->is_defined[cl_R] );
	      bkgndMap = CalcBkgndMap( B, B->IP->is_defined[cl_R] );
	      fragMap = CalcMotifFragMap( B, typ,  B->IP->is_defined[cl_R] ); 
	      r = motifMap + bkgndMap + fragMap; 

	      ratio[n].ratio = exp( r - currentMap );
	      ratio[n].count = count;
	      denom += ratio[n].ratio;
	      n++;
	    }
	} 
      count++;

      if( j > 0 )
	{
	  x = j;
	  c[j] = x;
	  j--;
	}
      else
	{
	  if( c[1] + 1 < c[2] )
	    c[1] = c[1] + 1;
	  else
	    {
	      j = 2;
	      c[j-1] = j - 2;
	      x = c[j] + 1;
	    
	      while( x == c[j+1] )
		{
		  j++;
		  c[j-1] = j - 2;
		  x = c[j] + 1;
	      }
	    
	      if( j > IP->nMotifLen[typ] )
		break;
	      
	      c[j] = x;
	      j--;
	    }
	} 
      
      if( n % RATIO_ALLOC_SIZE == 0 )
	ratio = realloc( ratio, ((n / RATIO_ALLOC_SIZE) + 1) * RATIO_ALLOC_SIZE * sizeof( RatioStruct ) );
    }

  if( denom < drand() )
    {
      for( i = 0; i < F->nMaxLen[typ]; i++ )
	F->nColMask[typ][i] = fragSave[i];
      F->FragWidth[typ] = fragWidthSave;

      for( i = 0; i <= B->IP->nMotifLen[typ]; i++ )
	{
	  for( j = 0; j < B->IP->nAlphaLen; j++ )
	    C->wCounts[typ][i][j] = countSave[i][j];
	}
      FREEP( countSave, IP->nMotifLen[typ] + 1 );
      free( fragSave );
      free( c );  
      free( ratio );
      B->Phylo->phyloTree = phyloSave;
      IP->is_defined[cl_B] = bkgndSave;
      
      free( lgIndex );
      
      return FALSE;
    }

  r = drand() * denom;
  for( sum = 0.0, j = 0; j < n; j++ )
    {
      sum += ratio[j].ratio;
      if( sum >= r )
	{
	  selected = ratio[j].count;
	  retVal = (ratio[j].ratio != 1.0);
	  break;
	}
    }
  
  for( j = 1; j <= IP->nMotifLen[typ]; j++ )
    c[j] = j - 1;

  c[IP->nMotifLen[typ]+1] = F->nMaxLen[typ];
  c[IP->nMotifLen[typ]+2] = 0;

  count = 0;
  j = IP->nMotifLen[typ];
  while( TRUE )
    {
      if( count == selected )
	{
	  SetFragmentation( B, typ, offset, c, lgIndex );
	  break;
	}
      count++;

      if( j > 0 )
	{
	  x = j;
	  c[j] = x;
	  j--;
	}
      else
	{
	  if( c[1] + 1 < c[2] )
	    c[1] = c[1] + 1;
	  else
	    {
	      j = 2;
	      c[j-1] = j - 2;
	      x = c[j] + 1;
	      
	      while( x == c[j+1] )
		{
		  j++;
		  c[j-1] = j - 2;
		  x = c[j] + 1;
		}
	      
	      if( j > IP->nMotifLen[typ] )
		p_internal_error( "Fragment: Invalid fragmentation" );
	      
	      c[j] = x;
	      j--;
	    }
	} 
    }
  
  for( i = 0; i <= B->IP->nMotifLen[typ]; i++ )
    {
      for( j = 0; j < B->IP->nAlphaLen; j++ )
	C->wCounts[typ][i][j] = countSave[i][j];
    }
  FREEP( countSave, IP->nMotifLen[typ] + 1 );

  free( fragSave );
  free( c );  
  free( ratio );
  B->Phylo->phyloTree = phyloSave;
  B->Phylo->bCalcPhylo = phyloCalcSave;
  IP->is_defined[cl_B] = bkgndSave;

  free( lgIndex );
  
  return retVal;
}


int ChooseColumns( Model B, int typ, int width, int offset, int *lgIndex )
{
  IPtype      IP;
  Ftype       F;
  Ctype       C;
  int         i;
  int         j;
  long double denom;
  long double sum;
  long double *lg;
  long double r;
  int         normalize = FALSE;
  int         newWidth;
  
  IP = B->IP;
  F = B->F;
  C = B->C;
  
  NEW( lg, width, long double );

  for( i = 0; i < width; i++ )
    {
      lgIndex[i] = i;
      if( F->nColMask[typ][i+offset] == COL_ON || F->nColMask[typ][i+offset] == COL_OFF )
	{
	  for( j = 0; j < IP->nAlphaLen; j++ )
	    lg[i] += lgamma( F->nvFragCnts[typ][i+offset][j] + C->dPseudoCounts[typ][0][j] );
	  lg[i] = expl( lg[i] );
	}
    }

  if( width > IP->nMotifLen[typ] && 
      lnbico( width, IP->nMotifLen[typ] ) > LMAX_SAMPLES ) 
    normalize = TRUE;

  newWidth = min( width, IP->nMotifLen[typ] + 4 );
  if( normalize )
    {
      for( i = 0; i < newWidth; i++ )
	{
	  for( denom = 0.0, j = 0; j < width; j++ )
	    denom += lg[j];
	  sum = 0;
	  r = ((long double) drand()) * denom;
	  for( j = 0; j < width; j++ )
	    {
	      sum += lg[j];
	      if( sum > r )
		{
		  lgIndex[i] = j;
		  lg[j] = 0;
		  break;
		}
	    }
	}
      width = IP->nMotifLen[typ] + 4;
    }
  qsort( lgIndex, newWidth, sizeof( int ), intcompare );
  
  free( lg );

  return newWidth;
}


void SetFragmentation( Model B, int typ, int offset, int *c, int *lgIndex )
{
  IPtype IP;
  Ctype  C;
  Ftype  F;
  int    i;
  int    j;

  IP = B->IP;
  C = B->C;
  F = B->F;

  for( i = 0; i < F->nMaxLen[typ]; i++ )
    {
      if( F->nColMask[typ][i] == COL_ON )
	{
	  F->nColMask[typ][i] = COL_OFF;
	  for( j = 0; j < B->IP->nAlphaLen; j++ )
	    {
	      C->wCounts[typ][BG][j] += F->nvFragCnts[typ][i][j];
	    }
	}      
    }
  
  for( i = 1; i <= IP->nMotifLen[typ]; i++ )
    {
      F->nColMask[typ][lgIndex[c[i]]+offset] = COL_ON;
      for( j = 0; j < B->IP->nAlphaLen; j++ )
	{
	  C->wCounts[typ][i][j] = F->nvFragCnts[typ][lgIndex[c[i]]+offset][j];
	  C->wCounts[typ][BG][j] -= F->nvFragCnts[typ][lgIndex[c[i]]+offset][j];
	}
    }

  F->FragWidth[typ] = lgIndex[c[IP->nMotifLen[typ]]] - lgIndex[c[1]] + 1;
}


int PossibleFragmentation( Model B, int typ, int mid, int offset, 
			   int firstPal, int lastPal, int firstRep, int lastRep, int *c, int *lgIndex )
{
  IPtype IP;
  Ftype  F;
  int    i;
  int    sum;
  int    diff;

  IP = B->IP;
  F = B->F;

  sum = lgIndex[c[firstPal]] + lgIndex[c[lastPal]];
  diff = lgIndex[c[lastRep]] - lgIndex[c[firstRep]];

  for( i = 1; i <= IP->nMotifLen[typ]; i++ )
    {
      if( F->nColMask[typ][lgIndex[c[i]]+offset] == COL_BLOCKED || 
	  F->nColMask[typ][lgIndex[c[i]]+offset] == COL_DONT_USE )
	return FALSE;
      
      if( i <= mid )
	{
	  if( IP->is_defined[cl_J] && i < mid )
	    {	 
	      if( lgIndex[c[i+1]] - lgIndex[c[i]] != 1 )
		return FALSE;
	      if( lgIndex[c[mid+i+1]] - lgIndex[c[mid+i]] != 1 )
		return FALSE;
	    }
	      
	  if( IP->AltModel->Palandromic[typ][i-1] ) 
	    {
	      if( lgIndex[c[IP->nMotifLen[typ] - i + 1]] + lgIndex[c[i]] != sum )
		return FALSE;
	    }
	  else if( IP->AltModel->Repeat[typ][i-1] )
	    {
	      if( lgIndex[c[mid + i]] - lgIndex[c[i]] != diff )
		return FALSE;
	    }
	}
    }

  return TRUE;
}


/* Make sure that fragmenting dooesn't cause an overlap with any motif.  */
/* This is called after fragmentation is set */

int FragmentOverlap( Model B, int t, Mlist M )
{
  MlistEl      *curr;
  MlistEl      *curr2;
  int          **old;
  int          **new;
  int          f1;
  int          f2;
  int          shift;
  int          revshift;
  int          l1;
  int          l2;
  int          seq;
  int          t2;
  int          pos;
  int          end;
  int          pos2;
  int          end2;

  old = B->F->nOldColMask;
  new = B->F->nColMask;
  f1 = FirstCol(B->F, t);
  l1 = LastCol(B->F, t);
  B->F->nColMask = old;
  f2 = FirstCol(B->F, t);
  l2 = LastCol(B->F, t);
  shift = f1 - f2;
  revshift = l1 - l2;
  
  B->F->nColMask = new;			/* BT 2/7/97 */

  curr = M[t]->Motifs;
  while( curr != NULL )
    {
      if( ! curr->RevComp ) 
	pos = curr->pos + shift; 
      else               
	pos = curr->pos - revshift;
      
      if( pos < 0 )
	return TRUE;
      
      seq = SequenceFromPosition( B, curr->pos );
      if( seq != SequenceFromPosition( B, pos ) )
	return TRUE;

      end = pos + B->F->FragWidth[t] - 1;
      if( seq != SequenceFromPosition( B, end ) )
	return TRUE;
      
      for(t2 = 0; t2 < B->IP->nNumMotifTypes; t2++) 
	{
	  curr2 = M[t2]->Motifs;
	  while( curr2 != NULL )
	    {
	      if( curr != curr2 )
		{
		  pos2 = curr2->pos + B->F->shift[t2];
		  end2 = pos2 + B->F->FragWidth[t2] - 1;
		  if( pos < pos2 && end >= pos2 )
		    return TRUE;
		  if( pos > pos2 && pos <= end2 )
		    return TRUE;
		}
	      curr2 = curr2->next;
	    }
	}
      curr = curr->next;
    }

  return FALSE;
}


int FirstFragPos( Model B, int typ )
{
  Ftype  F;
  int    i;
  
  F = B->F;

  for( i = 0; i < F->nMaxLen[typ]; i++ )
    { 
      if( F->nColMask[typ][i] == COL_ON || F->nColMask[typ][i] == COL_OFF ) 
	return i;
    }

  return -1;
}


int LastFragPos( Model B, int typ )
{
  Ftype  F;
  int    i;
  
  F = B->F;

  for( i = F->nMaxLen[typ] - 1; i >= 0; i-- )
    { 
      if( F->nColMask[typ][i] == COL_ON || F->nColMask[typ][i] == COL_OFF ) 
	return i;
    }

  return -1;
}


int FirstPalPos( Model B, int typ, int mid )
{
  IPtype IP;
  int    i;

  IP = B->IP;
  
  for( i = 1; i <= mid; i++ )    
    {      
      if( IP->AltModel->Palandromic[typ][i-1] ) 
	return i;
    }
  
  return 0;
}


int LastPalPos( Model B, int typ, int mid )
{
  IPtype IP;
  int    i;

  IP = B->IP;
  
  for( i = IP->nMotifLen[typ]; i > mid; i-- )    
    {      
      if( IP->AltModel->Palandromic[typ][i-1] ) 
	return i;
    }
  
  return 0;
}


int FirstRepPos( Model B, int typ, int mid )
{
  IPtype IP;
  int    i;

  IP = B->IP;
  
  for( i = 1; i <= mid; i++ )    
    {      
      if( IP->AltModel->Repeat[typ][i-1] ) 
	return i;
    }
  
  return 0;
}


int LastRepPos( Model B, int typ, int mid )
{
  IPtype IP;
  int    i;

  IP = B->IP;
  
  for( i = mid + 1; i <= IP->nMotifLen[typ]; i++ )    
    {      
      if( IP->AltModel->Repeat[typ][i-1] ) 
	return i;
    }
  
  return 0;
}


int intcompare(const void *p1, const void *p2)
{
  int i = *((int *)p1);
  int j = *((int *)p2);
  
  if (i > j)
    return (1);
  if (i < j)
    return (-1);
  return (0);
}


/* =========================================================================== */
/* Old style fragmentation - kept to maintain compatbility with older versions */
/* BT 04/26/04                                                                 */
/* =========================================================================== */


int Frag( Model B, int typ, Mlist M  )
{
  int lemon, i, fcols, retval=FALSE;
  int *Mark;

  NEW(Mark, B->F->nMaxLen[typ], int);
  
  if( B->IP->is_defined[cl_J] && 
      (! B->IP->is_defined[cl_R] || B->IP->nMotifLen[typ] % 2 == 0) )
    {
      FragmentFromMiddle( B, typ, Mark, M );
    }
  else
    {
      fcols = max(B->IP->nMotifLen[typ]/2,2);
      fcols = 1;
      for(i = 0; i < fcols; i++) 
	{ 
	  lemon = LemonFModel( B, typ, Mark );	/*choose an On column to eliminate*/
	  if( lemon < 0 )
	    break;
	  retval = move_column_gibbs(B, lemon, typ, Mark, M);
	  if(!retval)
	    break;
	}
    }
  free(Mark);
  return retval;
}


/* Sample a column in model proportional to how BAD it is. */
int	LemonFModel( Model B, int typ, int *Mark )
{
  double  r,rand,total;
  int	  i,j,k;
  int	  colnum=0,start;
  double  *tmp_val;
  int     mid;
  int     lemon_comp;
  
  NEW(tmp_val,B->F->nMaxLen[typ],double);
  start = FirstCol(B->F,typ);

  mid = B->IP->nMotifLen[typ] / 2 + B->IP->nMotifLen[typ] % 2;
		
  for (total=0.0,i=0; i<B->F->nMaxLen[typ]; i++)
    {
      if (B->F->nColMask[typ][i] == COL_ON) 
	{
	  if( !Mark[i] ) 
	    {
	      if( (B->IP->AltModel->Palandromic[typ][colnum] ||  
		   B->IP->AltModel->Repeat[typ][colnum]) &&
		  colnum >= mid )	    
		tmp_val[i] = 0.0;	  
	      else if( B->IP->AltModel->Repeat[typ][colnum] )
		{
		  if(  CheckRepeat( i, -1, B, typ, Mark) )
		    {
		      r = CalcRatio( B, typ , i, -1, Mark ); 
		      tmp_val[i] = r; 
		      total += r;		  
		    }
		  else
		    tmp_val[i] = 0.0;
		}
	      else
		{
		  r = CalcRatio( B, typ , i, -1, Mark ); 
		  /* r = RatioFModel( start, (int)0, i, colnum, M, typ ); */
		  lemon_comp = LastCol( B->F, typ );
		  for( j = 0; j < colnum; j++ )
		    {
		      lemon_comp = PrevCol( B->F, typ, lemon_comp);
		    }

		  if( (!B->IP->AltModel->Palandromic[typ][colnum]) ||
		      ((B->IP->AltModel->Palandromic[typ][colnum]) && 
		       B->F->nColMask[typ][lemon_comp] == COL_ON) )  
		    {
		      tmp_val[i] = r; 
		      total += r;
		    }
		  else
		    tmp_val[i] = 0.0;	  		    
		}
	    }
	  colnum++;
	}
      else 
	tmp_val[i] = 0.0;
    }

  rand = ((double) Random()/(double) INT32_MAX)*total;
  for (k=0,i=0; i< B->F->nMaxLen[typ]; i++,k++)
    {
      if ((B->F->nColMask[typ][i]==COL_ON) && !Mark[i])
	{
	  rand -= tmp_val[k];
	  if(rand <= 0.0) {Mark[k] = TRUE; free(tmp_val); return k;}
	}
    }

  p_warning(" Fragmentation: LemonFModel( )... this should not happen.");
  free(tmp_val); 
  return( -1 );
}


double	RatioFModel(int col1, int ps_col1, int col2, int ps_col2, Model B, int typ )
/* Return the ratio of new to old observed.  */
{
  double   sum,sumlgamma;
  int	   b;
  double   *mobserved,*observed;
  int      col1_pal, ps_col1_pal;
  int      col2_pal, ps_col2_pal;
  double   s1;
  double   dCnt;
  double   ps_sum;
  double   psum;
  int      i;

/* fix pseudo counts to correlate with the correct ratio values */

  mobserved = B->F->nvFragCnts[typ][col2];
  psum = 0.0;
  dCnt = 0.0;
  for(sumlgamma=0.0,b=0;b<B->IP->nAlphaLen; b++) 
    {
      if (B->C->dPseudoCounts[typ][ps_col2][b]!=0.0)
	{
	  s1 = mobserved[b]+
	    B->C->dPseudoCounts[typ][ps_col2][b];
	  dCnt += mobserved[b];
	  ps_sum = B->C->dPseudoCounts[typ][ps_col2][b];
	  if( B->IP->AltModel->Palandromic[typ][ps_col2] || B->IP->AltModel->Repeat[typ][ps_col2] )
	    {
	      ps_col2_pal = B->IP->nMotifLen[typ] - ps_col2 - 1;
	      
	      col2_pal = LastCol( B->F, typ ) - (col2 - FirstCol( B->F, typ ));

	      if( B->F->nColMask[typ][col2_pal] != COL_BLOCKED )		
		{
		  if( B->IP->AltModel->Palandromic[typ][ps_col2] )
		    {
		      s1 += B->F->nvFragCnts[typ][col2_pal][nComp[b]]+
			B->C->dPseudoCounts[typ][ps_col2_pal][nComp[b]];
		      dCnt += B->F->nvFragCnts[typ][col2_pal][nComp[b]];
		      ps_sum += B->C->dPseudoCounts[typ][ps_col2_pal][nComp[b]];
		    }
		  else if( B->IP->AltModel->Repeat[typ][ps_col2] )
		    {
		      s1 += B->F->nvFragCnts[typ][col2_pal][b]+
			B->C->dPseudoCounts[typ][ps_col2_pal][b];
		      dCnt += B->F->nvFragCnts[typ][col2_pal][b];
		      ps_sum += B->C->dPseudoCounts[typ][ps_col2_pal][b];
		    }
		}
	      else
		return 0.0;
	    }
	  psum += ps_sum;
	  sumlgamma += ln_gamma( s1 ) - ln_gamma( ps_sum );
	}
    }
  sumlgamma += ln_gamma( psum ) - ln_gamma( dCnt + psum );
  
  observed = B->F->nvFragCnts[typ][col1];
  psum = 0.0;
  dCnt = 0.0;
  for(sum=0.0,b=0;b<B->IP->nAlphaLen; b++) 
    {
      if(B->C->dPseudoCounts[typ][ps_col1][b]!=0.0) 
	{
	  s1 = observed[b]+
	    B->C->dPseudoCounts[typ][ps_col1][b];
	  dCnt += mobserved[b];
	  ps_sum = B->C->dPseudoCounts[typ][ps_col1][b];
	  if( B->IP->AltModel->Palandromic[typ][ps_col1] || B->IP->AltModel->Repeat[typ][ps_col1] )
	    {
	      ps_col1_pal = B->IP->nMotifLen[typ] - ps_col1 - 1;

	      col1_pal = LastCol( B->F, typ );
	      for( i = 0; i < col1; i++ )
		{
		  col1_pal = PrevCol( B->F, typ, col1_pal);
		}

	      if( B->F->nColMask[typ][col1_pal] != COL_BLOCKED )		
		{
		  if( B->IP->AltModel->Palandromic[typ][ps_col1] )
		    {
		      s1 +=  B->F->nvFragCnts[typ][col1_pal][nComp[b]]+
			B->C->dPseudoCounts[typ][ps_col1_pal][nComp[b]];
		      dCnt += B->F->nvFragCnts[typ][col1_pal][nComp[b]];
		      ps_sum += B->C->dPseudoCounts[typ][ps_col1_pal][nComp[b]];
		    }
		  else if( B->IP->AltModel->Repeat[typ][ps_col1] )
		    {
		      s1 +=  B->F->nvFragCnts[typ][col1_pal][b]+
			B->C->dPseudoCounts[typ][ps_col1_pal][b];
		      dCnt += B->F->nvFragCnts[typ][col1_pal][b];
		      ps_sum += B->C->dPseudoCounts[typ][ps_col1_pal][b];
		    }
		}
	      else
		return 0.0;
	    }
	  psum += ps_sum;
	  sum += ln_gamma( s1 ) - ln_gamma( ps_sum );
	}
    }
  sum += ln_gamma( psum ) - ln_gamma( dCnt + psum );
  
  return exp(sum-sumlgamma);
}



double	mv_col_weight(int c,int w, int n, int lemon, Model B, int typ)
/********************************************************************

   If a fmodel has c on columns and a length of w.  Then there are 
   w-2 choose c-2  configurations for that model.

          n     1 lemon       w_old
          |     |  |           |
        ........*??*???????????*.............
                 |--- (w-2) --|
               .....
              -10234
      n = -inf..0 or +2..+inf.	 (i.e., n != 1)

       lemon = 1..w and n = -f..end.

 ********************************************************************/
{
	int start=1, end=w, w_new = 0,w_old,first_col;

	w_old = w;
	first_col = FirstCol(B->F,typ);
	lemon = lemon - first_col + 1;
	
	if(lemon == 1) {     /** case 1: **/
	   if(n < 1){          /** case 1a: increase width by |n| **/
		w_new = w_old - n + 1;  /** note: n is <= 0 **/
	   } else if(n > 1){   /** case 1b: **/
		for(start=1;B->F->nColMask[typ][start+first_col] != COL_ON; start++) ;
		if(n < w_old){
		    start = min(start,n);
		    w_new = w_old - start + 1;
		} else if(n > w_old){
		    w_new = n - start + 1;   
		} else p_error("Fragment1: input error in mv_col_weight( )");
	   } else p_error("Fragment2: input error in mv_col_weight( )");
	} else if(lemon == w_old){	/** case 2: **/
	   for(end = w_old-1;B->F->nColMask[typ][end+first_col-1] == COL_OFF; end--) ;
	   if(n < 1){
		w_new = end - n + 1;  /** note: n is <= 0 **/
	   } else if(n > 1){
		if(n < end) w_new = end;
		else if(n > end) w_new = n; 
		else p_error("Fragment3: input error in mv_col_weight( )");
	   } else p_error("Fragment4: input error in mv_col_weight( )");
	} else { 		/** case 3: 1 < lemon < w_old **/
	   if(n < 1){
		w_new = w_old - n + 1;  /** note: n is <= 0 **/
	   } else if(n > 1){
		if(n < w_old) w_new = w_old; 
		else if(n > w_old) w_new = n;
		else p_error("Fragment5: input error in mv_col_weight( )");
	   } else p_error("Fragment6: input error in mv_col_weight( )");
	}

	return (bico(w_old-2,c-2)/bico(w_new-2,c-2));
}


int CheckRepeat(int lemon, int pos, Model B, int typ, int *Mark)
{
  int	oldstart=0;
  int	i,end;
  int   lemon_pos = 0;
  int   origin;
  int   lemon_comp;
  int   pos_comp;
  int   pos2 = -1;
  int   colnum = 0;
  int   zero_count = 0;
  int   mid;
  int   lemon_rep = 0;
  int   midPos;
  int   pos_rep;
  int   start;
  int   j;
  int   *fragVec;
  int   repeatOK = TRUE;
  int   repCnt;

  if( pos >= B->F->nMaxLen[typ] )
    return FALSE;

  if( lemon >= B->F->nMaxLen[typ] )
    return FALSE;

  NEW( fragVec, B->F->nMaxLen[typ], int );

  for( i = 0; i < B->F->nMaxLen[typ]; i++ )
    fragVec[i] = B->F->nColMask[typ][i];

  origin = FirstCol(B->F,typ); 
  for (i=origin; i<=lemon; i++)
    {
      if (B->F->nColMask[typ][i]==COL_OFF) 
	zero_count++;
    }
  lemon_pos = lemon-origin-zero_count; 
  lemon_comp = LastCol( B->F, typ );
  for( i = 0; i < lemon_pos; i++ )
    {
      lemon_comp = PrevCol( B->F, typ, lemon_comp);
    }

  mid = B->IP->nMotifLen[typ] / 2 + B->IP->nMotifLen[typ] % 2 - 1;
  midPos = FirstCol( B->F, typ );
  for( i = 0; i < mid; i++ )
    {
      midPos = NextCol( B->F, typ, midPos );
    }

  lemon_rep = midPos;
  for( i = 0; i <= lemon_pos; i++ )
    {
      lemon_rep = NextCol( B->F, typ, lemon_rep );
    }
  if( lemon_rep == LastCol( B->F, typ ) )
    lemon_rep = midPos + 1;
  else if( pos < lemon )
    lemon_rep++;

  for (end=origin,i=origin; i<B->F->nMaxLen[typ]; i++)
    if (B->F->nColMask[typ][i]!=COL_BLOCKED) end=i;

  for( i=0; i <= pos; i++)
    {
      if  (B->F->nColMask[typ][i] == COL_ON) colnum++;
      if  ((B->F->nColMask[typ][i]==COL_OFF)) 
	{
	  if ( i < lemon )
	    pos2 = colnum;
	  else
	    pos2 = colnum-1;
	}
    }
  if( pos2 == -1 )
    pos2 = colnum - 1;
  pos_comp = LastCol( B->F, typ ) - (pos - FirstCol( B->F, typ ));
  if( B->IP->is_defined[cl_J] )
    pos_rep = LastCol( B->F, typ ) + 1;
  else
    {
      pos_rep = midPos;
      for( i = 0; i <= pos; i++ )
	{
	  pos_rep = NextCol( B->F, typ, pos_rep );
	}
    }    
	  
  oldstart = FirstCol(B->F,typ);
  fragVec[lemon] = COL_OFF;
  if( pos != - 1 )
    {
      if( pos >= B->F->nMaxLen[typ] )
	repeatOK = FALSE;
      else
	fragVec[pos] = COL_ON;
    }
      
  if( B->IP->is_defined[cl_J] )
    {
      if( B->IP->AltModel->Repeat[typ][lemon_pos] )
	{
	  if( lemon_rep >= B->F->nMaxLen[typ] )
	    repeatOK = FALSE;
	  else
	    fragVec[lemon_rep] = COL_OFF;
	}
      
      if( pos != -1 && B->IP->AltModel->Repeat[typ][pos2] )
	{
	  if( pos_rep >= B->F->nMaxLen[typ] )
	    repeatOK = FALSE;
	  else if( fragVec[pos_rep] != COL_OFF )
	    repeatOK = FALSE;
	  else
	    fragVec[pos_rep] = COL_ON;
	}
    }
  else
    {
      if( B->IP->AltModel->Repeat[typ][lemon_pos] )
	{
	  i = 0;
	  start = 0;
	  while (B->F->nColMask[typ][i++] == COL_BLOCKED)  
	    start++;

	  midPos = FirstCol( B->F, typ );
	  for( i = 0; i < mid; i++ )
	    {
	      midPos = NextCol( B->F, typ, midPos );
	    }

	  for( j = 1, i = start; i <= midPos && j + midPos < B->F->nMaxLen[typ]; i++, j++ )
	    {
	      if( j + midPos >= B->F->nMaxLen[typ] )
		repeatOK = FALSE;
	      else if( fragVec[i] == COL_BLOCKED )
		repeatOK = FALSE;
	      else if( B->F->nColMask[typ][midPos + j] == COL_BLOCKED )
		repeatOK = FALSE;
	      else
		fragVec[midPos + j] = fragVec[i];
	    }

	  for( i = midPos + j; i < B->F->nMaxLen[typ]; i++ )
	    {
	      if( B->F->nColMask[typ][i] != COL_BLOCKED )
		fragVec[i]= COL_OFF;
	    }
	}      
    }

  if( pos != -1 && repeatOK )
    {
      for( repCnt = 0, i = 0; i < B->F->nMaxLen[typ]; i++ )
	{
	  if( fragVec[i] == COL_ON )
	    repCnt++;
	}
      if( repCnt != B->IP->nMotifLen[typ] )
	repeatOK = FALSE;
    }
  
  free( fragVec );

  return repeatOK;
}



int	MvColumnFModel(int lemon, int pos, Model B, int typ, int *Mark)
/*  Remove column at position lemon and add observed to pos in M */
{
  int	newstart=0,oldstart=0;
  int	i,end;
  int   lemon_pos = 0;
  int   origin;
  int   lemon_comp;
  int   pos_comp;
  int   pos2 = -1;
  int   colnum = 0;
  int   zero_count = 0;
  int   mid;
  int   lemon_rep = 0;
  int   midPos;
  int   pos_rep;
  int   start;
  int   j;

  origin = FirstCol(B->F,typ); 
  for (i=origin; i<=lemon; i++)
    {
      if (B->F->nColMask[typ][i]==COL_OFF) 
	zero_count++;
    }
  lemon_pos = lemon-origin-zero_count; 
  lemon_comp = LastCol( B->F, typ );
  for( i = 0; i < lemon_pos; i++ )
    {
      lemon_comp = PrevCol( B->F, typ, lemon_comp);
    }

  mid = B->IP->nMotifLen[typ] / 2 + B->IP->nMotifLen[typ] % 2 - 1;
  midPos = FirstCol( B->F, typ );
  for( i = 0; i < mid; i++ )
    {
      midPos = NextCol( B->F, typ, midPos );
    }

  lemon_rep = midPos;
  for( i = 0; i <= lemon_pos; i++ )
    {
      lemon_rep = NextCol( B->F, typ, lemon_rep );
    }
  if( lemon_rep == LastCol( B->F, typ ) )
    lemon_rep = midPos + 1;
  else if( pos < lemon )
    lemon_rep++;

  for (end=origin,i=origin; i<B->F->nMaxLen[typ]; i++)
    if (B->F->nColMask[typ][i]!=COL_BLOCKED) end=i;

  for( i=0; i <= pos; i++)
    {
      if  (B->F->nColMask[typ][i] == COL_ON) colnum++;
      if  ((B->F->nColMask[typ][i]==COL_OFF)) 
	{
	  if ( i < lemon )
	    pos2 = colnum;
	  else
	    pos2 = colnum-1;
	}
    }
  if( pos2 == -1 )
    pos2 = colnum - 1;
  pos_comp = LastCol( B->F, typ ) - (pos - FirstCol( B->F, typ ));
  if( B->IP->is_defined[cl_J] )
    pos_rep = LastCol( B->F, typ ) + 1;
  else
    {
      pos_rep = midPos;
      for( i = 0; i <= pos; i++ )
	{
	  pos_rep = NextCol( B->F, typ, pos_rep );
	}
    }    
	  
  oldstart = FirstCol(B->F,typ);
  B->F->nColMask[typ][lemon] = COL_OFF;
  if( pos != - 1 )
    B->F->nColMask[typ][pos] = COL_ON;
      
  if( B->IP->AltModel->Palandromic[typ][lemon_pos] )
    {
      B->F->nColMask[typ][lemon_comp] = COL_OFF;
    }

  if( pos != -1 && B->IP->AltModel->Palandromic[typ][pos2] )
    {
      B->F->nColMask[typ][pos_comp] = COL_ON;            
    }

  if( B->IP->is_defined[cl_J] )
    {
      if( B->IP->AltModel->Repeat[typ][lemon_pos] )
	{
	  B->F->nColMask[typ][lemon_rep] = COL_OFF;
	}
      
      if( pos != -1 && B->IP->AltModel->Repeat[typ][pos2] )
	{
	  B->F->nColMask[typ][pos_rep] = COL_ON;            
	}
    }
  else
    {
      if( B->IP->AltModel->Repeat[typ][lemon_pos] )
	{
	  i = 0;
	  start = 0;
	  while (B->F->nColMask[typ][i++] == COL_BLOCKED)  
	    start++;

	  midPos = FirstCol( B->F, typ );
	  for( i = 0; i < mid; i++ )
	    {
	      midPos = NextCol( B->F, typ, midPos );
	    } 

	  for( j = 1, i = start; i <= midPos && j + midPos < B->F->nMaxLen[typ]; i++, j++ )
	    {
	      B->F->nColMask[typ][midPos + j] = B->F->nColMask[typ][i];
	    }
	  for( i = midPos + j; i < B->F->nMaxLen[typ]; i++ )
	    {
	      if( B->F->nColMask[typ][i] != COL_BLOCKED )
		B->F->nColMask[typ][i]= COL_OFF;
	    }
	}      
    }
#ifdef _DEBUG_
  CheckCounts( B );
  for( j = 0; j <  B->IP->nMotifLen[typ]; j++ )
    {
      for( i = 0; i < B->IP->nAlphaLen; i++ )
	{
	  if( B->C->fCounts[typ][j+1][i] < -EPS )
	    p_error( "Negative frequency count." );	   
	}
    }
#endif

  newstart = FirstCol(B->F,typ);
  for (i=newstart,end=newstart;i<B->F->nMaxLen[typ];i++)
    if (B->F->nColMask[typ][i] == COL_ON) end = i;

  B->F->FragWidth[typ] = end - newstart + 1;

  return  (B->F->shift[typ] = newstart - oldstart);
}


double CalcWeight( Model B, int t, int lemon, int new_col, int *Mark )
{
  double c1;
  double c2;
  int    *fragMaskSave;
  int    widthSave;
  int    shiftSave;
  int    i;

  NEW( fragMaskSave, B->F->nMaxLen[t], int );

  for( i = 0; i < B->F->nMaxLen[t]; i++ )
    {
      fragMaskSave[i] = B->F->nColMask[t][i];
    }
  widthSave = B->F->FragWidth[t];
  shiftSave = B->F->shift[t];
 
  c1 = CalcMotifFragMap( B, t,  B->IP->is_defined[cl_R] );

  MvColumnFModel( lemon, new_col, B, t, Mark);

  c2 = CalcMotifFragMap( B, t,  B->IP->is_defined[cl_R] );

  /*
  printf( "new = %d old = %d  Frag c1 = %f c2 = %f diff = %f\n", 
	  new_col, lemon, c1, c2, c2 - c1 );
  */

  for( i = 0; i < B->F->nMaxLen[t]; i++ )
    {
      B->F->nColMask[t][i] = fragMaskSave[i];
    }
  free( fragMaskSave );
  B->F->FragWidth[t] = widthSave;
  B->F->shift[t] = shiftSave;

  return( exp(c2 - c1) );
}


double CalcRatio( Model B, int t, int lemon, int new_col, int *Mark )
{
  double b1;
  double b2;
  double c1;
  double c2;  
  int    *fragMaskSave;
  int    widthSave;
  int    shiftSave;
  int    i;
  int    j;
  int    cnt;
  double **countSave;
  double **pcountSave;
  int    phyloSave;

  phyloSave = B->Phylo->bCalcPhylo;
  B->Phylo->bCalcPhylo = FALSE;

  NEW( fragMaskSave, B->F->nMaxLen[t], int );

  NEWP( countSave, B->IP->nMotifLen[t] + 1, double );
  for( i = 0; i <= B->IP->nMotifLen[t]; i++ )
    NEW( countSave[i], B->IP->nAlphaLen, double );

  NEWP( pcountSave, B->IP->nMotifLen[t] + 1, double );
  for( i = 0; i <= B->IP->nMotifLen[t]; i++ )
    NEW( pcountSave[i], B->IP->nAlphaLen, double );

  for( i = 0; i < B->F->nMaxLen[t]; i++ )
    {
      fragMaskSave[i] = B->F->nColMask[t][i];
    }
  widthSave = B->F->FragWidth[t];
  shiftSave = B->F->shift[t];

  for( i = 0; i <= B->IP->nMotifLen[t]; i++ )
    {
      for( j = 0; j < B->IP->nAlphaLen; j++ )
	{
	  countSave[i][j] = B->C->wCounts[t][i][j];
	  pcountSave[i][j] = B->C->dPseudoCounts[t][i][j];	 
	} 
    }
  
  c1 = CalcMotifMap( B, t, B->IP->is_defined[cl_R] );
  b1 = CalcBkgndMap( B, B->IP->is_defined[cl_R] );

  MvColumnFModel( lemon, new_col, B, t, Mark);

  cnt = 1;
  for( i = 0; i < B->F->nMaxLen[t]; i++ )
    {
      if( i == lemon )
	{
	  for( j = 0; j < B->IP->nAlphaLen; j++ )
	    {
	      B->C->wCounts[t][BG][j] += B->F->nvFragCnts[t][i][j];
	    }
	}
      if( i == new_col )
	{
	  for( j = 0; j < B->IP->nAlphaLen; j++ )
	    {
	      B->C->wCounts[t][BG][j] -= B->F->nvFragCnts[t][i][j];
	    }
	}
      if( B->F->nColMask[t][i] == COL_ON )
	{
	  for( j = 0; j < B->IP->nAlphaLen; j++ )
	    {
	      B->C->wCounts[t][cnt][j] = B->F->nvFragCnts[t][i][j];
	    }
	  cnt++;
	}
    }

  c2 = CalcMotifMap( B, t, B->IP->is_defined[cl_R] );
  b2 = CalcBkgndMap( B, B->IP->is_defined[cl_R] );

  /*
  printf( "new = %d old = %d  Map c1 = %f b1 = %f c2 = %f b2 = %f diff = %f\n", 
	  new_col, lemon, c1, b1, c2, b2, (c2 + b2) - (c1 + b1) );
  */

  for( i = 0; i < B->F->nMaxLen[t]; i++ )
    {
      B->F->nColMask[t][i] = fragMaskSave[i];
    }
  free( fragMaskSave );
  B->F->FragWidth[t] = widthSave;
  B->F->shift[t] = shiftSave;

  for( i = 0; i <= B->IP->nMotifLen[t]; i++ )
    {
      for( j = 0; j < B->IP->nAlphaLen; j++ )
	{
	  B->C->wCounts[t][i][j] = countSave[i][j];
	  B->C->dPseudoCounts[t][i][j] = pcountSave[i][j];	  
	}
    }
  FREEP( countSave, B->IP->nMotifLen[t] + 1 );
  FREEP( pcountSave, B->IP->nMotifLen[t] + 1 );

  B->Phylo->bCalcPhylo = phyloSave;

  return( exp((c2 + b2) - (c1 + b1)) );
}



char	move_column_gibbs(Model B, int lemon, int typ, int *Mark, Mlist M )
/*******************************************************************
  Sample a column to remove and then sample a column to replace it.
  if lemon = new site then don't bother to move ??? 
********************************************************************/
{
  int	  i=0,d,leng,col_ref,end;
  double  *ratio,total,rand_no;
  int	  colnum=0,start=0,zero_count=0;
  int	  ncol,origin,old_col_pos,new_col_pos;
  double  weight;
  int     mid;
  int     pos_comp;
  int     j;
  int     frag_mid;
  int     lemon_comp;
  short   nonPal;
  int     f;
  int     l;
  int     *colpos;
  int     *tMask;
  double  r;
  double  sum;
  int     repCnt;
  int     isPal = FALSE;
	
  NEW(ratio,B->F->nMaxLen[typ],double);
  NEW(colpos,B->F->nMaxLen[typ],int);

  ncol = B->IP->nMotifLen[typ];
  leng = B->F->FragWidth[typ];
  while (B->F->nColMask[typ][i++] == COL_BLOCKED)  start++;
  origin = FirstCol(B->F,typ); 
  for (i=origin;i<=lemon;i++)
    if (B->F->nColMask[typ][i]==COL_OFF) zero_count++;
  old_col_pos = lemon-origin-zero_count; 
  
  for (end=origin,i=origin;i<B->F->nMaxLen[typ];i++)
    if (B->F->nColMask[typ][i]!=COL_BLOCKED) end=i;

  mid = B->IP->nMotifLen[typ] / 2 + B->IP->nMotifLen[typ] % 2 - 1; /* BT 11/01/2001 */
  frag_mid = B->F->nMaxLen[typ] / 2 + B->F->nMaxLen[typ] % 2;	  

  for( i = 0; i < ncol; i++ )
    {
      if( B->IP->AltModel->Palandromic[typ][i] )
	{
	  isPal = TRUE;
	  break;
	}
    }

  if( isPal && old_col_pos > mid - (B->IP->nMotifLen[typ] % 2) )
    {
      free(ratio);
      free(colpos);
      
      return FALSE;
    }

  for(total=1.0,i=start;i<=end; i++)
    {
      if( B->F->nColMask[typ][i] == COL_ON ) 
	colnum++;
      if( (B->F->nColMask[typ][i]==COL_OFF) && 
	  ((! B->IP->is_defined[cl_J]) || (i < lemon))) 
	{
	  if ( i < lemon )
	    new_col_pos = colnum;
	  else
	    new_col_pos = colnum-1;

	  colpos[i] = new_col_pos;

	  pos_comp = LastCol( B->F, typ ) - (i - FirstCol( B->F, typ ));

	  if( B->F->nColMask[typ][i]==COL_DONT_USE )
	      ratio[i] = 0.0;
	  else if( B->IP->AltModel->Palandromic[typ][new_col_pos] && 
	      B->F->nColMask[typ][pos_comp] != COL_OFF )
	    ratio[i] = 0.0;	  
	  else if( B->IP->AltModel->Palandromic[typ][new_col_pos] && 
		   (new_col_pos > mid - (B->IP->nMotifLen[typ] % 2) ) ) /* BT 11/01/2001 */
	    ratio[i] = 0.0;
	  else if( B->IP->AltModel->Palandromic[typ][old_col_pos] && 
		   (old_col_pos > mid - (B->IP->nMotifLen[typ] % 2) ) ) /* BT 11/01/2001 */
	    ratio[i] = 0.0;              /* BT 11/01/2001 */
	  else if( B->IP->AltModel->Palandromic[typ][old_col_pos] && 
		   ! B->IP->AltModel->Palandromic[typ][new_col_pos] )
	    ratio[i] = 0.0;		   
	  else if( ! B->IP->AltModel->Palandromic[typ][old_col_pos] && 
		   B->IP->AltModel->Palandromic[typ][new_col_pos] )
	    ratio[i] = 0.0;
	  else if( B->IP->AltModel->Repeat[typ][new_col_pos] && 
		   (new_col_pos > mid - (B->IP->nMotifLen[typ] % 2) ) ) /* BT 11/01/2001 */
	    ratio[i] = 0.0;
	  else if( B->IP->AltModel->Repeat[typ][old_col_pos] && 
		   (old_col_pos > mid - (B->IP->nMotifLen[typ] % 2) )  )/* BT 11/01/2001 */
	    ratio[i] = 0.0;              /* BT 11/01/2001 */
	  else if( B->IP->AltModel->Repeat[typ][old_col_pos] && 
		   ! B->IP->AltModel->Repeat[typ][new_col_pos] )
	    ratio[i] = 0.0;		   
	  else if( ! B->IP->AltModel->Repeat[typ][old_col_pos] && 
		   B->IP->AltModel->Repeat[typ][new_col_pos] )
	    ratio[i] = 0.0;
	  else if( (B->IP->AltModel->Repeat[typ][old_col_pos] || 
		    B->IP->AltModel->Repeat[typ][new_col_pos] ) &&
		   ! CheckRepeat( lemon, i, B, typ, Mark ) )
	    ratio[i] = 0.0;
	  else if( Mark[i] )
	    ratio[i] = 0.0;
	  else if( (B->IP->AltModel->Palandromic[typ][new_col_pos] ||  
		    B->IP->AltModel->Repeat[typ][new_col_pos])  &&
		   pos_comp >= B->F->nMaxLen[typ] ) 
		   ratio[i] = 0.0;   /* BT 09/10/03 */	    
	  else
	    {		
	      nonPal = FALSE;
	      if( new_col_pos == old_col_pos  )
		{
		  NEW(tMask,B->F->nMaxLen[typ],int);
		  for( j = 0; j < B->F->nMaxLen[typ]; j++ )
		    {
		      tMask[j] = B->F->nColMask[typ][j];
		    }
		  tMask[i] = COL_ON;
		  tMask[lemon] = COL_OFF;
		  if( B->IP->AltModel->Palandromic[typ][new_col_pos] || B->IP->AltModel->Repeat[typ][new_col_pos] )
		    {
		      tMask[pos_comp] = COL_ON;
		      lemon_comp = LastCol( B->F, typ );
		      for( j = 0; j < old_col_pos; j++ )
			{
			  lemon_comp = PrevCol( B->F, typ, lemon_comp);
			}
		      tMask[lemon_comp] = COL_OFF;
		      f = 0;
		      while( tMask[f] != COL_ON )
			f++;
		      l = B->F->nMaxLen[typ] - 1;
		      while( tMask[l] != COL_ON )
			l--;
		      while( f <= l )
			{
			  if( tMask[f] != tMask[l] )
			    {
			      ratio[i] = 0.0;
			      nonPal = TRUE;			
			      break;
			    }
			  f++;
			  l--;
			}		  
		    }
		  
		  free( tMask );
		}
	      
	      if( ! nonPal && ! FragOverlap( B, typ, lemon, i, Mark, M))
		{
		  /* ratio[i] = RatioFModel(i,new_col_pos,lemon,old_col_pos, M, typ); */
		  ratio[i] = CalcRatio( B, typ, lemon, i, Mark );
		  col_ref = i-origin + 1;
		  /* weight = mv_col_weight(ncol, leng, col_ref, lemon, M, typ); */
		  weight = CalcWeight( B, typ, lemon, i, Mark );

		  ratio[i] *= weight;
		  total += ratio[i];
		}
	    }
	}
      else 
	ratio[i] = 0.0;
    }
  
  r = drand();
  rand_no = r*total;
  
  for(sum = 0, i=start;i<=end; i++)
    {
      /*      if(((rand_no -= ratio[i]) <= 0.0) && (ratio[i] != 0.0) && 
	      (! FragOverlap( B, typ, lemon, i, Mark, M) ) ) 
	      { */ /* BT 1018/2001 */
      sum += ratio[i];
      if( sum > rand_no && (ratio[i] != 0.0) && ! FragOverlap( B, typ, lemon, i, Mark, M) ) 
	{
	  d = MvColumnFModel(lemon,i,B,typ,Mark);	  
	  B->F->shift[typ] = d;
	  if( ! B->IP->is_defined[cl_Z] )
	    {
	      fprintf(stdout," (typ=%d; new=%d; old=%d; w=%d shift = %d) %g %g\n",
		      typ, i, lemon, B->F->FragWidth[typ], B->F->shift[typ], 
		      ratio[i], log( ratio[i] ));
	      printFragVec( B );
	    }

#ifdef _DEBUG_
	  if( old_col_pos == mid )
	    i = i;
#endif

	  for( repCnt = 0, j = 0; j < B->F->nMaxLen[typ]; j++ )
	    {
	      if( B->F->nColMask[typ][j] == COL_ON )
		repCnt++;
	    }
	  if( repCnt != B->IP->nMotifLen[typ] )
	    p_error( "Fragment: invalid number of columns" );

	  free(ratio);
	  free(colpos);
	  return TRUE;
	} 
    }

  free(ratio);
  free(colpos);

  return FALSE;
}


/* Make sure that fragmenting dooesn't cause an overlap with any motif. 
   If it does, return false and reject new_col as a possible column to enter the motif 
   BT 5/31/2001 */

int FragOverlap( Model B, int t, int lemon, int new_col, int *Mark,  Mlist M )
{
  int     *fragMaskSave;
  int     widthSave;
  int     shiftSave;
  int     i;
  int     shift;
  int     t2;
  int     pos;
  int     pos2;
  int     end;
  int     end2;
  int     overlap;
  MlistEl *curr; 
  MlistEl *curr2; 
  int     f1;
  int     f2;
  int     l1;
  int     l2;
  int     revshift;


  NEW( fragMaskSave, B->F->nMaxLen[t], int );

  for( i = 0; i < B->F->nMaxLen[t]; i++ )
    {
      fragMaskSave[i] = B->F->nColMask[t][i];
    }
  widthSave = B->F->FragWidth[t];
  shiftSave = B->F->shift[t];
 
  f1 = FirstCol(B->F, t);
  l1 = LastCol(B->F, t);

  shift = MvColumnFModel( lemon, new_col, B, t, Mark);
  f2 = FirstCol(B->F, t);
  l2 = LastCol(B->F, t);

  shift = f1 - f2;
  revshift = l1 - l2;

  overlap = FALSE;
  curr = M[t]->Motifs;
  while( (curr != NULL) && (! overlap) )
    {
      if( ! curr->RevComp ) 
	pos = curr->pos + shift; 
      else               
	pos = curr->pos - revshift;
      
      if( pos < 0 )
	overlap = TRUE;
      
      if( SequenceFromPosition( B, curr->pos ) != SequenceFromPosition( B, pos ) )
	overlap = TRUE;

      end = pos + B->F->FragWidth[t] - 1;
      if( SequenceFromPosition( B, curr->pos ) != SequenceFromPosition( B, end ) )
	overlap = TRUE;

      for(t2 = 0; (t2 < B->IP->nNumMotifTypes) && (! overlap); t2++) 
	{
	  curr2 = M[t2]->Motifs;
	  while( (curr2 != NULL) && (! overlap) )
	    {
	      if( curr != curr2 )
		{
		  pos2 = curr2->pos + B->F->shift[t2];
		  end2 = pos2 + B->F->FragWidth[t2] - 1;
		  if( pos < pos2 && end >= pos2 )
		    overlap = TRUE;
		  if( pos > pos2 && pos <= end2 )
		    overlap = TRUE;
		}
	      curr2 = curr2->next;
	    }
	}
      curr = curr->next;
    }

  for( i = 0; i < B->F->nMaxLen[t]; i++ )
    {
      B->F->nColMask[t][i] = fragMaskSave[i];
    }
  free( fragMaskSave );

  B->F->FragWidth[t] = widthSave;
  B->F->shift[t] = shiftSave;

  return( overlap );
}


void FragmentFromMiddle( Model B, int typ, int *Mark, Mlist M )
{
  char *temp;
  int  count;
  int  i;

  /*  if( ! FragmentFromMiddleLeftHalf( B, typ, Mark, M ) )
      FragmentFromMiddleRightHalf( B, typ, Mark, M ); */
  FragmentFromMiddleLeftHalf( B, typ, Mark, M );
  FragmentFromMiddleRightHalf( B, typ, Mark, M ); 

  for( count = 0, i = 0; i < B->F->nMaxLen[typ]; i++ )
    {
      if( B->F->nColMask[typ][i] == COL_ON )
	count++;
    }
  if( count != B->IP->nMotifLen[typ] )
    {
      NEW( temp, 1024, char );
      sprintf( temp,"FragmentFromMiddle: invalid column count motif %d count %d", typ, count );
      p_error( temp );
    }  
}


int FragmentFromMiddleLeftHalf( Model B, int typ, int *Mark, Mlist M )
{
  int    first;
  int    last;
  int    mid;
  int    i;
  int    leftshift;
  int    rightshift;
  int    lemon;
  double ratio;
  double weight = 1.0;
  int    done;
  int    shifted;
  int    in;
  double r;
  int    midshift;
  int    lemonPos;
  int    j;
  
 
  first = FirstCol( B->F, typ );
  last = LastCol( B->F, typ );
  mid = first + (last - first) / 2;

  leftshift = 0;
  i = first - 1;
  while( i >= 0 && B->F->nColMask[typ][i] == COL_OFF )
    {
      leftshift++;
      i--;
    }

  i = mid;
  midshift = mid;
  while( i >= 0 && B->F->nColMask[typ][i] == COL_OFF )
    {
      midshift--;
      i--;
    }

  done = FALSE;
  shifted = FALSE;
  for( i = 1; i <= leftshift && ! done  ; i++ )
    {
      lemon = midshift - i + 1;
      in = first - i;

      lemonPos = 0;
      for( j = 0; j <= lemon; j++ )
	{
	  if( B->F->nColMask[typ][i] == COL_ON) 
	    lemonPos++;
	}

      if( lemon < 0 || lemon >=  B->F->FragWidth[typ] || 
	  B->F->nColMask[typ][lemon] != COL_ON )
	{
	  done = TRUE;
#ifdef _DEBUG_
	  printf( "FragmentFromMiddleLeftHalf - 1: invalid old column: %d new column: %d typ: %d\n", 
		  lemon, in, typ );
#endif
	}

      if( in < 0 || in >=  B->F->FragWidth[typ] || 
	  B->F->nColMask[typ][in] != COL_OFF )
	{
	  done = TRUE;
#ifdef _DEBUG_	  
	  printf( "FragmentFromMiddleLeftHalf - 1: invalid new column: %d old column: %d typ: %d\n", 
		  lemon, in, typ );
#endif
	}
      
      if( ! done && B->F->nColMask[typ][lemon] != COL_DONT_USE &&
	  B->F->nColMask[typ][in] != COL_DONT_USE &&
	  ( ! B->IP->AltModel->Repeat[typ][lemonPos] ||
	    CheckRepeat( lemon, in, B, typ, Mark ) ) &&
	  ! FragOverlap( B, typ, lemon, in, Mark, M ) ) 
	{
	  ratio = CalcRatio( B, typ, lemon, in, Mark ); 
	  if( B->IP->nAlphaLen != 4 )   /* BT 04/12/04 */
	    weight = CalcWeight( B, typ, lemon, in, Mark );
	  ratio *= weight;
	}
    } 

  if( ! shifted )
    {
      done = FALSE;
      rightshift = mid - PrevCol( B->F, typ, mid + 1);
      for( i = 1; i <= rightshift && ! done; i++ )
	{
	  lemon = first + i - 1;
	  in = mid - rightshift + i;

	  lemonPos = 0;
	  for( j = 0; j <= lemon; j++ )
	    {
	      if( B->F->nColMask[typ][i] == COL_ON) 
		lemonPos++;
	    }
	  
	  if( lemon < 0 || lemon >=  B->F->FragWidth[typ] || 
	      B->F->nColMask[typ][lemon] != COL_ON )
	    {
	      done = TRUE;
#ifdef _DEBUG_	      
	      printf( "FragmentFromMiddleLeftHalf - 3: invalid old column: %d new column: %d typ: %d\n", 
		      lemon, in, typ );
#endif
	    }
	  
	  if( in < 0 || in >=  B->F->FragWidth[typ] || 
	      B->F->nColMask[typ][in] != COL_OFF )
	    {
#ifdef _DEBUG_
	      printf( "FragmentFromMiddleLeftHalf - 4: invalid new column: %d old column: %d typ: %d\n", 
		      lemon, in, typ );
#endif
	    }

	  if( B->IP->AltModel->Repeat[typ][in] && 
	      LastCol( B->F, typ ) + 1 >= B->F->nMaxLen[typ] )
	    done = TRUE;

	  if( ! done && B->F->nColMask[typ][lemon] != COL_DONT_USE &&
	      B->F->nColMask[typ][in] != COL_DONT_USE &&
	      ( ! B->IP->AltModel->Repeat[typ][lemonPos] ||
		CheckRepeat( lemon, in, B, typ, Mark ) ) &&
	      ! FragOverlap( B, typ, lemon, in, Mark, M ) ) 
	    {
	      ratio = CalcRatio( B, typ, lemon, in, Mark );      
	      if( B->IP->nAlphaLen != 4 )   /* BT 04/12/04 */
		weight = CalcWeight( B, typ, lemon, in, Mark );
	      ratio *= weight;
	      if( ratio >= 1.0 )
		{
		  MvColumnFModel(lemon,in,B,typ,Mark);
		  if( ! B->IP->is_defined[cl_Z] )
		    fprintf(stdout,"ll (typ=%d; new=%d; old=%d; w=%d; shift=%d) %g %g\n",
			    typ, in, lemon, B->F->FragWidth[typ],  B->F->shift[typ],
			    ratio, log( ratio ));
		  shifted = TRUE; 
		}
	      else
		{
		  r = drand();
		  if( ratio > r )
		    {
		      MvColumnFModel(lemon,in,B,typ,Mark);
		      if( ! B->IP->is_defined[cl_Z] )
			fprintf(stdout,"4l (typ=%d; new=%d; old=%d; w=%d; shift=%d) %g %g\n",
				typ, in, lemon, B->F->FragWidth[typ],  B->F->shift[typ],
				ratio, log( ratio ));
		      shifted = TRUE; 
		    }
		  else
		    done = TRUE;
		}
	    }
	  else
	    done = TRUE;      
	}
    }

 return( shifted );
}


int FragmentFromMiddleRightHalf( Model B, int typ, int *Mark, Mlist M )
{
  int    first;
  int    last;
  int    mid;
  int    i;
  int    leftshift;
  int    rightshift;
  int    lemon;
  double ratio;
  double weight = 1.0;
  int    done;
  int    shifted;
  int    in;
  double r;
  int    midshift = -1;
  int    lemonPos;
  int    j;

  first = FirstCol( B->F, typ );
  last = LastCol( B->F, typ );
  mid = first + (last - first) / 2;

  leftshift = 0;
  i = last;
  while( i >= 0 && B->F->nColMask[typ][i] == COL_ON )
    {
      i--;
    }
  if( i >= mid )
    {
      i++;
      midshift = i;
      while( i > mid && B->F->nColMask[typ][i] == COL_OFF )
	{
	  i--;
	  leftshift++;
	}
    }

  done = FALSE;
  shifted = FALSE;
  for( i = 1; i <= leftshift && ! done; i++ )
    {
      in = midshift - i;
      lemon = last + i;

      lemonPos = 0;
      for( j = 0; j <= lemon; j++ )
	{
	  if( B->F->nColMask[typ][i] == COL_ON) 
	    lemonPos++;
	}

      if( lemon < 0 || lemon >=  B->F->FragWidth[typ] || 
	  B->F->nColMask[typ][lemon] != COL_ON )
	{
	  done = TRUE;
#ifdef _DEBUG_	  
	  printf( "FragmentFromMiddleRightHalf - 1: invalid old column: %d new column: %d typ: %d\n", 
		  lemon, in, typ );
#endif
	}

      if( in < 0 || in >=  B->F->FragWidth[typ] || 
	  B->F->nColMask[typ][in] != COL_OFF )
	{
	  done = TRUE;
#ifdef _DEBUG_
	  printf( "FragmentFromMiddleRightHalf - 1: invalid new column: %d old column: %d typ: %d\n", 
		  lemon, in, typ );
#endif	  
	}
      
      if( ! done && ! B->IP->AltModel->Palandromic[typ][lemon] && ! B->IP->AltModel->Repeat[typ][lemon] )
	{
	  if( B->F->nColMask[typ][lemon] != COL_DONT_USE &&
	      B->F->nColMask[typ][in] != COL_DONT_USE &&
	      ( ! B->IP->AltModel->Repeat[typ][lemonPos] ||
		CheckRepeat( lemon, in, B, typ, Mark ) ) &&
	      ! FragOverlap( B, typ, lemon, in, Mark, M ) ) 
	    {
	      ratio = CalcRatio( B, typ, lemon, in, Mark );      
	      if( B->IP->nAlphaLen != 4 )   /* BT 04/12/04 */
		weight = CalcWeight( B, typ, lemon, in, Mark );
	      ratio *= weight;
	      if( ratio >= 1.0 )
		{
		  MvColumnFModel(lemon,in,B,typ,Mark);
		  if( ! B->IP->is_defined[cl_Z] )
		    fprintf(stdout,"1r (typ=%d; new=%d; old=%d; w=%d; shift=%d) %g %g\n",
			    typ, in, lemon, B->F->FragWidth[typ],  B->F->shift[typ],
			    ratio, log( ratio ));
		  shifted = TRUE;
		}
	      else
		{
		  r = drand();
		  if( ratio > r )
		    {
		      MvColumnFModel(lemon,in,B,typ,Mark);
		      if( ! B->IP->is_defined[cl_Z] )
			fprintf(stdout,"2r (typ=%d; new=%d; old=%d; w=%d; shift=%d) %g %g\n",
				typ, in, lemon, B->F->FragWidth[typ],  B->F->shift[typ],
				ratio, log( ratio ));
		      shifted = TRUE;
		    }
		  else
		    done = TRUE;
		}
	    }
	  else
	    done = TRUE;      
	}
    }

  if( ! shifted )
    {
      done = FALSE;
      rightshift = 0;
      i = last + 1;
      while(  i <  B->F->FragWidth[typ] && B->F->nColMask[typ][i] == COL_OFF )
	{
	  rightshift++;
	  i++;
	}

      for( i = 1; i <= rightshift && ! done; i++ )
	{
	  in = last + i;
	  lemon = mid + i;

	  lemonPos = 0;
	  for( j = 0; j <= lemon; j++ )
	    {
	      if( B->F->nColMask[typ][i] == COL_ON) 
		lemonPos++;
	    }

	  if( lemon < 0 || lemon >=  B->F->FragWidth[typ] || 
	      B->F->nColMask[typ][lemon] != COL_ON )
	    {
	      done = TRUE;
#ifdef _DEBUG_
	      printf( "FragmentFromMiddleRightHalf - 1: invalid old column: %d new column: %d typ: %d\n", 
		      lemon, in, typ );
#endif
	    }

	  if( in < 0 || in >=  B->F->FragWidth[typ] || 
	      B->F->nColMask[typ][in] != COL_OFF )
	    {
	      done = TRUE;
#ifdef _DEBUG_	      
	      printf( "FragmentFromMiddleRightHalf - 1: invalid new column: %d old column: %d typ: %d\n", 
		      lemon, in, typ );
#endif
	    }
	  
	  if( ! done && ! B->IP->AltModel->Palandromic[typ][lemon] && ! B->IP->AltModel->Repeat[typ][lemon] )
	    {
	      if( B->F->nColMask[typ][lemon] != COL_DONT_USE &&
		  B->F->nColMask[typ][in] != COL_DONT_USE &&
		  ( ! B->IP->AltModel->Repeat[typ][lemonPos] ||
		    CheckRepeat( lemon, in, B, typ, Mark ) ) &&
		  ! FragOverlap( B, typ, lemon, in, Mark, M ) ) 
		{
		  ratio = CalcRatio( B, typ, lemon, in, Mark );      
		  if( B->IP->nAlphaLen != 4 )   /* BT 04/12/04 */
		    weight = CalcWeight( B, typ, lemon, in, Mark );
		  ratio *= weight;
		  if( ratio >= 1.0 )
		    {
		      MvColumnFModel(lemon,in,B,typ,Mark);
		      if( ! B->IP->is_defined[cl_Z] )
			fprintf(stdout,"3r (typ=%d; new=%d; old=%d; w=%d; shift=%d) %g %g\n",
				typ, in, lemon, B->F->FragWidth[typ],  B->F->shift[typ],
				ratio, log( ratio ));
		      shifted = TRUE;
		    }
		  else
		    {
		      r = drand();
		      if( ratio > r )
			{
			  MvColumnFModel(lemon,in,B,typ,Mark);
			  if( ! B->IP->is_defined[cl_Z] )
			    fprintf(stdout,"4r (typ=%d; new=%d; old=%d; w=%d; shift=%d) %g %g\n",
				    typ, in, lemon, B->F->FragWidth[typ],  B->F->shift[typ],
				    ratio, log( ratio ));
			  shifted = TRUE;
			}
		      else
			done = TRUE;
		    }
		}
	      else
		done = TRUE;      
	    }
	}
    }
  return( shifted );      
}


double  bico(int N,int k)
{ return floor(0.5+exp(lnfact(N)-lnfact(k)-lnfact(N-k))); }

double  lnbico(register int N, register int k)
{ return lnfact(N)-lnfact(k)-lnfact(N-k); }

double	lnfact(int n)
/* static variables are guaranteed to be initialized to zero */
{
	static double lnft[101];

	if (n <= 1) return 0.0;
	if (n <= 100) return lnft[n] ? lnft[n] : (lnft[n]=lgamma(n+1.0));
	else return lgamma(n+1.0);
}
