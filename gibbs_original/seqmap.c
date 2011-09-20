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
#include "seqmap.h"

double CalcSeqMotifMap( Model B, int t, double ***dCounts );
double CalcSeqBetaMap( Model B, int seq );
double CalcSeqBkgndMap( Model B, int seq, double ***dCounts );
double CalcSeqSitePerSeqMap( Model B, int seq );
double CalcSeqNullMap( Model B, int seq, double ***dCounts );
void AddCounts( Model B, int begin_loc, int t, short RevComp, double ***dCounts );


double CalcSeqMapProb( Model B, int seq, double *motifMap )
{
  double ***dCounts;
  IPtype IP;
  RPType RP;
  int    t;
  int    j;
  int    n;
  int    nSeqStart;
  int    nLen;
  double result;
  double motifPart;
  double betaPart;
  double bgPart;
  double fragPart = 0.0;
  double seqPart;
  double nullPart;
  char   processed;
  char   *Res_ptr=B->Seq->R[0];  /* for quick access sequence residue */
  int    index;
  int    *seqMotifs;
  double fp;
  int    p;

  IP = B->IP;
  RP = B->RP;

  if( IsPhyloSeq( B, seq ) && seq % SpeciesInc( B, seq ) != 0 )
    p_error( "CalcSeqMapProb: Only the first sequences in a group should be used with -D option" );

  nSeqStart = SequenceStartPos( B, seq );
  nLen = SequenceLength( B, seq );

  NEWPP( dCounts, IP->nNumMotifTypes, double );
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      NEWP( dCounts[t], IP->nMotifLen[t] + 1, double );
      for( j = 0; j <= IP->nMotifLen[t]; j++ ) 
	{
	  NEW( dCounts[t][j], IP->nAlphaLen, double );
	}
    }
  NEW( seqMotifs, IP->nNumMotifTypes, int );

  for( n = 0; n < nLen; n++ )
    {
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  processed = Res_ptr[nSeqStart + n];
	  if( (IP->nAlphaLen == 20 && (processed != 'x') && (processed != 'X') &&
	       (processed != 'u')) ||
	      (IP->nAlphaLen == 4 && (processed != 'x') && (processed != 'X') && 
	       (processed != 'n') && (processed != 'N')) ) 
	    {
	      index = (int)(processed) - 97;
	      dCounts[t][BG][index] += 1.0;
	    }

	  if( IsPhyloSeq( B, seq ) )
	    {
	      for( p = 1; p < SpeciesInc( B, seq ); p++ )
		{
		  processed = Res_ptr[nSeqStart + n + p * nLen];
		  if( (IP->nAlphaLen == 20 && (processed != 'x') && (processed != 'X') &&
		       (processed != 'u')) ||
		      (IP->nAlphaLen == 4 && (processed != 'x') && (processed != 'X') && 
		       (processed != 'n') && (processed != 'N')) ) 
		    {
		      index = (int)(processed) - 97;
		      dCounts[t][BG][index] += 1.0;
		    }
		}
	    }
	}
    }

  for( n = 0; n < nLen; n++ )
    {
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  if( RP->sitePos[seq][n][t].nMotifStart )
	    {
	      seqMotifs[t]++;
	      AddCounts( B, nSeqStart + n, t, RP->sitePos[seq][n][t].nRevComp, dCounts );
	      if( IsPhyloSeq( B, seq ) )
		{
		  for( p = 1; p < SpeciesInc( B, seq ); p++ )
		    {
		      seqMotifs[t]++;
		      AddCounts( B, nSeqStart + n + p * nLen, t, RP->sitePos[seq+p][n][t].nRevComp, dCounts );		
		    }
		}
	    }
	}
    }

  nullPart = CalcSeqNullMap( B, seq, dCounts );

  /*************************/
  /* Motif element portion */
  /*************************/ 

  for( motifPart = 0, t=0; t < IP->nNumMotifTypes; t++)
    {    
      if( seqMotifs[t] ) 
	motifMap[t] = CalcSeqMotifMap( B, t, dCounts );
      else
	motifMap[t] = 0;
      motifPart += motifMap[t];
    }

  /*-----------*/
  /* Beta dist.*/
  /*-----------*/
  
  betaPart = CalcSeqBetaMap( B, seq ); 
  
  /*-------------------*/
  /* background portion*/
  /*-------------------*/

  bgPart = CalcSeqBkgndMap( B, seq, dCounts );

  /*-----------------------*/
  /* fragmentation portion */
  /*-----------------------*/

  if( ! B->IP->is_defined[cl_F] )
    {
      for( t=0; t < IP->nNumMotifTypes; t++)
	{
	  if( seqMotifs[t] )
	    {
	      fp = CalcMotifFragMap( B, t, TRUE );
	      fragPart += fp;
	      motifMap[t] += fp;
	    }
	}
    }

  seqPart =  CalcSeqSitePerSeqMap( B, seq );

  result = motifPart + betaPart + bgPart + fragPart + seqPart - nullPart; 
  /* result = motifPart + bgPart + fragPart + betaPart - nullPart; */

#ifdef _DEBUG_
  if( ! finite( result ) )
    {
      result = result;
    }
#endif

  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      FREEP( dCounts[t], IP->nMotifLen[t]+1 );
    }
  free( dCounts );
  free( seqMotifs );

  return result;
}


void AddCounts( Model B, int begin_loc, int t, short RevComp, double ***dCounts )
{
  int    i, j, index, mid;
  int    n, J, compindex, orig_index;
  int    last = -1, curr, first = -1;
  int    curr_fwd, curr_back, last_fwd = -1, last_back;
  char   processed;
  char   *Res_ptr=B->Seq->R[0];  /* for quick access sequence residue */
  int    compindex2;
  int    frag_width = 0;
   
  if(!B->IP->is_defined[cl_F]) 
    {
      first = FirstCol(B->F, t);
      last = last_fwd = -1;
      last_back = B->F->nMaxLen[t];
      frag_width = LastCol(B->F, t) - first + 1;
    }
  
  J = 0;
  for( i = 0; i < B->IP->nNumSequences; i++ )  /* BT 3/24/97 - allow overlap */
    {
      if( begin_loc < (*B->Seq->nvEndLocs)[i] )
	{
	  J = min( B->IP->nMotifLen[t], (*B->Seq->nvEndLocs)[i] - begin_loc );
	  break;
	}
    }
  
  if( J == 0 )
    return;		/* BT 4/18/97 */
  
  mid = J / 2 + J % 2 - 1;
  
  /*--------------------------*/
  /* DEALING WITH AMINO ACIDS */
  /*--------------------------*/
  
  if(B->IP->nAlphaLen == 20) 
    {
      for(j = 0; j < J; j++) 
	{
	  if(!B->IP->is_defined[cl_F]) 
	    {     /* Fragmentation is used */
	      curr = NextCol(B->F, t, last);
	      if(curr == -1) 
		p_internal_error("NO NEXT COLUMN");
	      last = curr;   curr -= first;
	      if(curr < 0) 
		p_internal_error("INVALID LOCATION");
	    }
	  else 
	    curr = j;
	  
	  processed = Res_ptr[begin_loc + curr];
	  if((processed != 'x') && (processed != 'X') && (processed != 'u'))
	    {
	      index = (int)(processed) - 97;
	      dCounts[t][j+1][index] += 1.0;
	      for(i = 0; i < B->IP->nNumMotifTypes; i++) 
		{
		  dCounts[i][BG][index] -= 1.0;
		}
	    }
	}
    }
  
  /*--------------------------*/
  /* DEALING WITH NUCLEOTIDES */
  /*--------------------------*/
  
  else 
    {    
      for(j = 0; j < J; j++) 
	{
	  if(!B->IP->is_defined[cl_F]) 
	    {/* Fragmentation */
	      curr_fwd = NextCol(B->F, t, last_fwd);
	      /* curr_back = PrevCol(B->F, t, last_back); */ /* 9/8/99 */
	      last_fwd = curr_fwd;   
	      /* last_back = curr_back; */ /* 9/8/99 */
	      curr_fwd -= first;     
	      /* curr_back -= first; */  /* 9/8/99 */
	      
	      curr_back = frag_width - curr_fwd - 1;   /* BT 08/11/99 */
	      
	      if((curr_fwd < 0) || (curr_fwd > B->F->nMaxLen[t])) 
		{
		  printf("first %d curr_fwd = %d\n", first, curr_fwd);
		  p_internal_error("INVALID FORWARD COLUMN");
		}
            if((curr_back < 0) || (curr_back > B->F->nMaxLen[t])) 
	      p_internal_error("INVALID BACKWARD COLUMN");
	    }
	  else 
	    {      /* No Fragmentation */
	      curr_fwd = j;
	      curr_back = J - j - 1;
	    }
	  if(!RevComp) {curr = curr_fwd;}
	  else         {curr = curr_back;}
	  
	  
	  processed = Res_ptr[begin_loc + curr_fwd];
	  orig_index = (int)processed - 97;
	  if((processed == 'X') || (processed == 'n') || 
            (processed == 'x') || (processed == 'N'))  /* BT 4/9/97 */
	    {
	      orig_index=4;
	    }
	  
	  if(!RevComp) 
	    { /* FORWARD MOTIF */
	      n = begin_loc + curr_fwd; 
	      index = orig_index;
	      compindex = CH2INTCOMP(n, B);
	      compindex2 = CH2INTCOMP(begin_loc + curr_back, B);	     
	    }
	  else 
	    { /* REVERSE COMPLEMENT */
	      n = begin_loc + curr_back;
	      index = CH2INTCOMP(n, B);
	      compindex = orig_index;
	      compindex2 = CH2INTCOMP(begin_loc + curr_fwd, B);	     
	    }
	  
	  dCounts[t][j+1][index] += 1.0;
	  for(i = 0; i < B->IP->nNumMotifTypes; i++) 
	    {
	      dCounts[i][BG][orig_index] -= 1.0;
	    }
	}      
    }
}


double CalcSeqMotifMap( Model B, int t, double ***dCounts )
{
  int        i;
  int        j;
  IPtype     IP;
  int        nNumMotifs;
  double     result;
  double     dPseudo; 
  double     value;
  double     **fCounts;       /* for quick access */
  double     **dPseudoCounts; /* for quick access */
  Ctype      C;
  double     dCnt;
  double     dTotalCnt;
  double     dPseudoCnt;
  int        mid;

  IP = B->IP;
  C = B->C;
  
  result = 0.0;
  value = 0.0;
  fCounts = dCounts[t];
  dPseudoCounts = C->dPseudoCounts[t];
  
  nNumMotifs = NUMMOTIFS(IP->nNumMotifs[t]); 
  
  /*************************/
  /* Motif element portion */
  /*************************/ 
  mid = IP->nMotifLen[t] / 2 + IP->nMotifLen[t] % 2;	  
  
  for(j = 1; j <= IP->nMotifLen[t]; j++) 
    {  
      dPseudo = 0.0;
      if( IP->AltModel->Palandromic[t][j-1] )
	{
	  if( j <= mid )
	    {
	      dTotalCnt = 0;
	      for(i = 0; i < IP->nAlphaLen; i++) 
		{     
		  dCnt = fCounts[j][i] + fCounts[IP->nMotifLen[t] - j + 1][nComp[i]];

		  dPseudoCnt = 
		    dPseudoCounts[j][i] + 
		    dPseudoCounts[IP->nMotifLen[t] - j + 1][nComp[i]];
		  value += ln_gamma(dCnt + dPseudoCnt ); 
		  value -= ln_gamma( dPseudoCnt );
		  dPseudo += dPseudoCnt;  

		  dTotalCnt += dCnt;
		} 
	      value += ln_gamma( dPseudo );
	      value -= ln_gamma( dTotalCnt + dPseudo );
	    }
	}
      else if( IP->AltModel->Repeat[t][j-1] )
	{
	  if( j <= mid )
	    {
	      dTotalCnt = 0;
	      for(i = 0; i < IP->nAlphaLen; i++) 
		{     
		  dCnt = fCounts[j][i] + fCounts[mid + j][i];

		  dPseudoCnt = 
		    dPseudoCounts[j][i] + 
		    dPseudoCounts[mid + j][i];
		  value += ln_gamma(dCnt + dPseudoCnt ); 
		  value -= ln_gamma( dPseudoCnt );
		  dPseudo += dPseudoCnt;  

		  dTotalCnt += dCnt;
		} 
	      value += ln_gamma( dPseudo );
	      value -= ln_gamma( dTotalCnt + dPseudo );
	    }
	}
      else if( IP->AltModel->Collapsed[t][j-1])  /* BT 12/28/99 */
      	{
	  dTotalCnt = 0;
	  for(i = 0; i < IP->nAlphaLen; i += 2) 
	    {     
	      dCnt = fCounts[j][i] + fCounts[j][nComp[i]];
	      
	      dPseudoCnt = 
		dPseudoCounts[j][i] + 
		dPseudoCounts[j][nComp[i]];
	      value += ln_gamma(dCnt + dPseudoCnt ); 
	      value -= ln_gamma( dPseudoCnt );
	      value -= ln_gamma( dCnt + 2 );
	      value += ln_gamma( (double) fCounts[j][i] + 1.0 ) + 
		ln_gamma( (double) fCounts[j][i+1] + 1.0 );
	      dPseudo += dPseudoCnt;
	      
	      dTotalCnt += dCnt;
	    } 
	  value += ln_gamma( dPseudo );
	  value -= ln_gamma( dTotalCnt + dPseudo );
	}
      else
	{
	  dTotalCnt = 0;
	  for(i = 0; i < IP->nAlphaLen; i++) 
	    {     
	      value += ln_gamma( fCounts[j][i] + dPseudoCounts[j][i] );
	      value -= ln_gamma( dPseudoCounts[j][i] );
	      dPseudo += dPseudoCounts[j][i];
	      dTotalCnt += fCounts[j][i];
	    } 
	  value += ln_gamma( dPseudo );
	  value -= ln_gamma( dTotalCnt + dPseudo );
	}
    } 

  result = value;
        
  return result;
}


double CalcSeqBetaMap( Model B, int seq )
{
  int        i;
  int        j;
  int        t;
  IPtype     IP;
  int        nonsites, nNumMotifs;
  double     alpha;
  double     dBetaPart;
  Ctype      C;
  RPType     RP;
  int        nSites;
  int        nLen = 0;
  double     weight;
  double     dTotPseudo;
  double     dBgPseudo;
  int        nPrevSite = -1;
  int        nPrevType = -1;
  int        nPrevDir = 0;
  double     dPartialSum;
  int        p;

  IP = B->IP;
  C = B->C;
  RP = B->RP;

  /*-----------*/
  /* Beta dist.*/
  /*-----------*/
  
  dBetaPart = 0.0;
  if( IP->is_defined[cl_E] )  
    {
      for( i = seq; i <= seq + SpeciesInc( B, seq ) - 1; i++ )
	{
	  if( RP->nUseSpacingProb || IP->is_defined[cl_T] )   /* BT 11/17/99 */
	    {
	      nPrevType = -1;
	      nLen = SequenceLength( B, i );
	      nSites = 0;
	      for( j = nLen - 1; j >= 0; j-- )
		{
		  for( t = 0; t < B->IP->nNumMotifTypes; t++ ) 
		    {
		      if( RP->sitePos[i][j][t].nMotifStart )
			{
			  dBetaPart += log( RP->sitePos[i][j][t].dSpacingProb );
			  dPartialSum = CalcBetaPartialSum(B, i, t, nPrevSite, nPrevType, j, nSites );
			  if( dPartialSum > 0 )
			    dBetaPart -= log( dPartialSum ); 
			  
			  if( nPrevType != -1 )
			    {
			      if( RP->sitePos[i][j][t].nRevComp )  /* BT 2/15/2001 */
				dBetaPart += log( RP->dPTrans[nPrevType][nPrevDir][t][REVERSE] );
			      else
				dBetaPart += log( RP->dPTrans[nPrevType][nPrevDir][t][FORWARD] ); 
			    }
			  
			  nSites++;
			  nPrevSite = j;
			  nPrevType = t;
			  if( RP->sitePos[i][j][t].nRevComp ) 
			    nPrevDir = REVERSE;
			  else
			    nPrevDir = FORWARD;			    
			}
		    }
		}
	    }      
	  else
	    {
	      nSites = min( B->RP->nSites[i], B->RP->nMaxBlocks );
	      if( B->AP->dAlignCount[nSites][i] != -DBL_MAX )
		dBetaPart += -log( B->AP->dAlignCount[nSites][i] ); 
	    } 
	}
    }
  else
    {  
      if( IP->is_defined[cl_T] )    /* BT 11/10/99 */
	{
	  weight = IP->dPseudoSiteWt / (1.0 - IP->dPseudoSiteWt);
	  for( t=0; t < IP->nNumMotifTypes; t++)
	    {
	      nNumMotifs = RP->nSites[seq];
	      nonsites = nLen - nNumMotifs;
	      if( IsPhyloSeq( B, seq ) )
		{
		  for( p = 1; p < SpeciesInc( B, seq ); p++ )
		    {
		      nNumMotifs += RP->nSites[seq+p];
		      nonsites += nLen - nNumMotifs;
		    }
		}
	      dTotPseudo = RP->nonCutoffSites * weight;
	      dBgPseudo = dTotPseudo - 
		(B->First->dmodel_sites[t][FORWARD] + B->First->dmodel_sites[t][REVERSE])  * weight;
	      alpha = dTotPseudo - dBgPseudo;

	      dBetaPart += LnBeta((double)nNumMotifs + alpha,             
				  (double)nonsites + dBgPseudo) -  
		LnBeta((double)alpha, dBgPseudo );       
	    }	  
	}
      else
	{
	  for( t=0; t < IP->nNumMotifTypes; t++)
	    {
	      nNumMotifs = RP->nSites[seq];
	      nonsites = nLen - nNumMotifs;
	      if( IsPhyloSeq( B, seq ) )
		{
		  for( p = 1; p < SpeciesInc( B, seq ); p++ )
		    {
		      nNumMotifs += RP->nSites[seq+p];
		      nonsites += nLen - nNumMotifs;
		    }
		}
	      
	      alpha = C->dtot_pseudo[t] - C->dbg_pseudo[t];
	      
	      dBetaPart += LnBeta((double)nNumMotifs + alpha,             
				  (double)nonsites + C->dbg_pseudo[t]) -  
		LnBeta((double)alpha, C->dbg_pseudo[t]);       
	    }
	}
    }

  return dBetaPart;
}


double CalcSeqBkgndMap( Model B, int seq, double ***dCounts )
{
  double     **fCounts;       /* for quick access */
  double     **dPseudoCounts; /* for quick access */
  double     bgValue; 
  double     bgPseudo; 
  double     total; 
  int        i;
  IPtype     IP;
  Ctype      C;
  RPType     RP;
  int        nLen;
  int        t;
  int        j;
  BkgType    BP;
  int        start;
  double     *bgCounts;
  double     *dPseudo;
  int        p;

  IP = B->IP;
  C = B->C;
  RP = B->RP;
  BP = B->BP;

  fCounts = dCounts[BG];
  dPseudoCounts = C->dPseudoCounts[BG];

  NEW( dPseudo, IP->nAlphaLen, double );  /* BT 08/06/2001 */
  for( i=0; i < IP->nAlphaLen; i++ ) 
    {
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  dPseudo[i] += C->dPseudoCounts[t][BG][i];
	}
    }
  
  bgValue = 0.0;
  bgPseudo = 0.0;
  total = 0.0;

  if( IP->is_defined[cl_B] )
    {
      NEW( bgCounts, IP->nAlphaLen, double );

      nLen = SequenceLength( B, seq );
      start = SequenceStartPos( B, seq );
      j = 0;
      while( j < nLen )
	{
	  if( PossibleBkgndPos( B, start + j, nLen ) )
	    {
	      for( i = 0; i < IP->nAlphaLen; i++)
		{
		  bgCounts[i] += BP->dBkgndProb[seq][j][i];		  
		  if( IsPhyloSeq( B, seq ) )
		    for( p = 1; p < SpeciesInc( B, seq ); p++ )
		      bgCounts[i] += BP->dBkgndProb[seq+p][j][i];
		}
	    }
	  j++;
	}
      
      for(i=0; i < IP->nAlphaLen; i++) 
	{                    
	  bgValue += ln_gamma(bgCounts[i] +  
			      dPseudo[i]);
	  bgValue -= ln_gamma( dPseudo[i] );      	      
	  total += bgCounts[i];
	  bgPseudo += dPseudo[i];   /* BT 8/7/98 */
	}
      free( bgCounts );
    }
  else
    {    
      for(i=0; i < IP->nAlphaLen; i++) 
	{                    
	  bgValue += ln_gamma((double)fCounts[BG][i] +  
			      dPseudo[i]);
	  bgValue -= ln_gamma( dPseudo[i] );      
	  total += fCounts[BG][i];
	  bgPseudo += dPseudo[i];   /* BT 8/7/98 */
	}
    }

  bgValue -= ln_gamma(total + bgPseudo);    /* BT 8/7/98 */
  bgValue += ln_gamma( bgPseudo ) ;

  free( dPseudo );

  return( bgValue );
}


double CalcSeqSitePerSeqMap( Model B, int seq )
{
  RPType RP;
  IPtype IP;
  double *cntSitesPerSeq;
  int    i;
  int    k;
  double siteMap = 0.0;
  double cntSum = 0.0;
  double pseudoSum = 0.0;

  if( B->IP->is_defined[cl_E] && B->RP->nUseFixedBlockProbs )
    {
      RP = B->RP;
      IP = B->IP;
      
      NEW( cntSitesPerSeq,  RP->nMaxBlocks + 1, double );
      
      for( i = seq; i <= seq + SpeciesInc( B, seq ) - 1; i++ )
	{
	  if( IP->is_defined[cl_V] )
	    {
	      if( i >= IP->nVerifySeq )
		cntSitesPerSeq[RP->nSites[i]]++;
	    }
	  else
	    cntSitesPerSeq[RP->nSites[i]]++;	    
	}
      
      for( k = 0; k <= RP->nMaxBlocks; k++ )
	{
	  siteMap += ln_gamma( cntSitesPerSeq[k] + RP->priorSitesPerSeq[k] );
	  siteMap -= ln_gamma( RP->priorSitesPerSeq[k] );
	  cntSum += cntSitesPerSeq[k];
	  pseudoSum += RP->priorSitesPerSeq[k];
	}
      siteMap += ln_gamma( pseudoSum );
      siteMap -= ln_gamma( cntSum + pseudoSum );
      
      free( cntSitesPerSeq );
    }

  return( siteMap );
}


double CalcSeqNullMap( Model B, int seq, double ***dCounts )
{
  double     **fCounts;       /* for quick access */
  double     **dPseudoCounts; /* for quick access */
  double     bgValue; 
  double     bgPseudo; 
  double     total; 
  int        i;
  IPtype     IP;
  Ctype      C;
  RPType     RP;
  int        nLen;
  int        t;
  int        j;
  BkgType    BP;
  int        start;
  double     *bgCounts;
  double     *dPseudo;
  int        nSitesSave;
  int        nSitesSave2 = 0;
  int        p;

  IP = B->IP;
  C = B->C;
  RP = B->RP;
  BP = B->BP;

  fCounts = dCounts[BG];
  dPseudoCounts=C->dPseudoCounts[BG];

  NEW( dPseudo, IP->nAlphaLen, double );  /* BT 08/06/2001 */
  for( i=0; i < IP->nAlphaLen; i++ ) 
    {
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  dPseudo[i] += C->dPseudoCounts[t][BG][i];
	}
    }
  
  bgValue = 0.0;
  bgPseudo = 0.0;
  total = 0.0;

  if( IP->is_defined[cl_B] )
    {
      NEW( bgCounts, IP->nAlphaLen, double );

      nLen = SequenceLength( B, seq );
      start = SequenceStartPos( B, seq );
      j = 0;
      while( j < nLen )
	{
	  if( PossibleBkgndPos( B, start + j, nLen ) )
	    {
	      for( i = 0; i < IP->nAlphaLen; i++)
		{
		  bgCounts[i] += BP->dBkgndProb[seq][j][i];
		  if( IsPhyloSeq( B, seq ) )
		    {
		      for( p = 1; p < SpeciesInc( B, seq ); p++ )
			bgCounts[i] += BP->dBkgndProb[seq+p][j][i];
		    }
		}
	    }
	  j++;
	}
      
      for(i=0; i < IP->nAlphaLen; i++) 
	{                    
	  bgValue += ln_gamma(bgCounts[i] +  
			      dPseudo[i]);
	  bgValue -= ln_gamma( dPseudo[i] );      
	  total += bgCounts[i];
	  bgPseudo += dPseudo[i];   /* BT 8/7/98 */
	}
      free( bgCounts );
    }
  else
    {    
      for(i=0; i < IP->nAlphaLen; i++) 
	{                    
	  bgValue += ln_gamma((double)fCounts[BG][i] +  
			      dPseudo[i]);
	  bgValue -= ln_gamma( dPseudo[i] );      
	  total += fCounts[BG][i];
	  bgPseudo += dPseudo[i];   /* BT 8/7/98 */
	}
    }

  bgValue -= ln_gamma(total + bgPseudo);    /* BT 8/7/98 */
  bgValue += ln_gamma( bgPseudo ) ;

  nSitesSave = RP->nSites[seq];
  RP->nSites[seq] = 0;
  if( IsPhyloSeq( B, seq ) )
    {
      for( p = 1; p < SpeciesInc( B, seq ); p++ )
	{
	  nSitesSave2 = RP->nSites[seq+p];
	  RP->nSites[seq+p] = 0;
	}
    }
  
  bgValue += CalcSeqSitePerSeqMap( B, seq );

  RP->nSites[seq] = nSitesSave; 
  if( IsPhyloSeq( B, seq ) )
    for( p = 1; p < SpeciesInc( B, seq ); p++ )
      RP->nSites[seq+p] = nSitesSave2; 

  free( dPseudo );

  return( bgValue );
}
