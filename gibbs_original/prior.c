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
#include "prior.h"

static char *delims = ",\t ";
static char cIncludeFile[MAX_LENGTH];
static int bkgndPseudo = FALSE;
 
void ReadPriorFile2( Model B, Ctype C, FILE* fpt );
void ReadPriors( Model B, Ctype C, int t, FILE* fpt );
void ReadPseudo( Model B, Ctype C, int t, FILE* fpt );
void ReadGapPenalty( Model B, FILE* fpt );
void ReadTransProbability( Model B, double dTransWt, int nUpDate, int nSym, FILE* fpt );
void ReadEndPenalty( Model B, int t, FILE* fpt );
void ReadFootPrintData( Model B, int nSeq, FILE* fpt );
void ReadSequenceBlocks( Model B, FILE* fpt );
void ReadGroups( Model B, FILE* fpt );
void ReadSiteProb( Model B, FILE* fpt );
void ReadFlatGapProb( Model B, FILE* fpt );
void ReadBlockProb( Model B, double weight, FILE* fpt );
void ReadSitePseudo( Model B, double weight, FILE* fpt );
void ReadInitPos( Model B, int seq, FILE* fpt );
void ReadBackProb( Model B, FILE* fpt );
void ReadSpacingProb( Model B, FILE* fpt );
void ReadSimpleBkgnd( Model B, FILE* fpt );
void ReadSeqWeights( Model B, FILE* fpt );
void ReadWeights( Model B, FILE* fpt );
void ReadWeightsFromFile( Model B, FILE* wtFile );
void ReadEndSiteProb( Model B, FILE* fpt );
void ReadBeginSiteProb( Model B, FILE* fpt );
void ReadPosMatrix( Model B, double dTransWt, int nUpDate, FILE* fpt );
void ReadMask( Model B, int t, FILE* fpt );
void ReadComment( Model B, FILE* fpt );
void ReadGibbsName( Model B, FILE* fpt );
void GetBkgndPseudo( Model B, Ctype C, int t, FILE* fpt );
void ReadMinSiteMap( Model B, FILE* fpt );
void ReadKSampleMap( Model B, FILE* fpt );
void ReadAlignWeight( Model B, FILE* fpt );
void ReadTree( Model B, int treeStartSeq, int treeEndSeq, int speciesSample, FILE* fpt );
void ReadBkgndFreqFile( Model B, char *seqStr, FILE *fpt );
void RunOut( FILE *fpt, char *line, int maxLength );
double Poisson( int n, double mean );
void ReCalcBkgndPseudoCnts( Model B, Ctype C );
void CheckNmers( Model B );


void ReadPriorFile( Model B, Ctype C )
{
  IPtype IP;

  IP = B->IP;

  if( IP->Datafiles->prior_fpt == NULL )
    return;

  ReadPriorFile2( B, C, IP->Datafiles->prior_fpt );
  
  /* if backgrnd comp is used, we will recalc after reading file */
  if( ! IP->is_defined[cl_B] && ! bkgndPseudo ) 
    ReCalcBkgndPseudoCnts( B, C );   

  if( (B->IP->is_defined[cl_D] || B->Phylo->treeCount > 0) && B->IP->is_defined[cl_freq_background] )
    CheckNmers( B );
}


void ReadPriorFile2( Model B, Ctype C, FILE* fpt )
{
  IPtype IP;
  char   *line;
  char   *msg;
  char   *token;
  char   *ptr = NULL;
  int    t;
  FILE*  fpt_save;
  double wt;
  char   *ptr2 = NULL;
  int    nUpDate;
#ifndef _NO_PHYLO_
  int    treeStartSeq;
  int    treeEndSeq;
#endif
  int    sym;
  int    nDontAdjustPrior;

  IP = B->IP;

  if( fpt == NULL )
    return;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );
  if( GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    ptr = GetAToken( line, token );

  while( ptr != NULL )
    {
      if( strncasecmp( token, ">PRIOR", 6 ) == 0 )
	{
	  ptr = GetAToken( ptr, token );
	  t = max( atoi( token ) - 1, 0 );
	  if( t < IP->nNumMotifTypes )	    
	    ReadPriors( B, C, t, fpt );
	  else
	    RunOut( fpt, line, MAX_LENGTH );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">PSEUDO", 7 ) == 0 )
	{
	  ptr = GetAToken( ptr, token );
	  t = max( atoi( token ) - 1, 0 );
	  if( t < IP->nNumMotifTypes )	    
	    ReadPseudo( B, C, t, fpt );
	  else
	    RunOut( fpt, line, MAX_LENGTH );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">SEQ", 4 ) == 0 )
	{
	  ReadSequenceBlocks( B, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">CONS", 4 ) == 0 )
	{
	  ptr = GetAToken( ptr, token );
	  t = max( atoi( token ) - 1, 0 );
	  if( t < IP->nNumSequences )	    
	    ReadFootPrintData( B, t, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">BLOCK", 6 ) == 0 ) /* BT 0526/99 */
	{
	  ptr = GetAToken( ptr, token );
	  if( ptr != NULL )
	    wt = strtod( token, &ptr2);
	  else
	    wt = DEF_BLOCK_WT;
	  ptr = GetAToken( ptr, token );
	  if( ptr != NULL )
	    nDontAdjustPrior = atoi( token );
	  else
	    nDontAdjustPrior = FALSE;
	  ReadBlockProb( B, wt, fpt );
	  if( nDontAdjustPrior )
	    {
	      B->RP->nScaleZeroBlocksProb = TRUE;
	      B->RP->nDontAdjustPrior = TRUE;
	    }
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
       else if( strncasecmp( token, ">SITEPSEUDO", 11 ) == 0 ) 
	{
	  ptr = GetAToken( ptr, token );
	  if( ptr != NULL )
	    wt = strtod( token, &ptr2);
	  else
	    wt = DEF_SPSEUDO_WT;
	  ReadSitePseudo( B, wt, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">SBLOCK", 7 ) == 0 ) /* use to raise prob of 0 blocks */
	{
	  ReadBlockProb( B, DEF_BLOCK_WT, fpt );
 	  B->RP->nScaleZeroBlocksProb = TRUE;
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">SITE", 5 ) == 0 )
	{
	  ReadSiteProb( B, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">INIT", 5 ) == 0 )
	{
	  ptr = GetAToken( ptr, token );
	  t = max( atoi( token ) - 1, 0 );
	  if( t < IP->nNumSequences )	    
	    ReadInitPos( B, t, fpt );
	  else
	    RunOut( fpt, line, MAX_LENGTH );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">BACK", 5 ) == 0 )
	{
	  ReadBackProb( B, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">SPACING", 8 ) == 0 )
	{
	  ReadSpacingProb( B, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">BKGND", 6 ) == 0 )
	{
	  ReadSimpleBkgnd( B, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">COMMENT", 8 ) == 0 )
	{
	  ReadComment( B, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">GIBBS", 8 ) == 0 )
	{
	  ReadGibbsName( B, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">BGPSEUDO", 5 ) == 0 )
	{
	  ptr = GetAToken( ptr, token );
	  t = max( atoi( token ) - 1, 0 );
	  if( t < IP->nNumMotifTypes )	    
	    GetBkgndPseudo( B, C, t, fpt );
	  else
	    RunOut( fpt, line, MAX_LENGTH );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
#ifndef _NOTRANS_
      else if( strncasecmp( token, ">FLAT", 5 ) == 0 )
	{
	  ReadFlatGapProb( B, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">GROUP", 6 ) == 0 )
	{
	  ReadGroups( B, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strcasecmp( token, ">ENDGAP" ) == 0 )
	{
	  ptr = GetAToken( ptr, token );
	  t = max( atoi( token ) - 1, 0 );
	  if( t < IP->nNumMotifTypes )	    
	    ReadEndPenalty( B, t, fpt );
	  else
	    RunOut( fpt, line, MAX_LENGTH );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">TREEWEIGHT", 7 ) == 0 )
	{
	  ReadWeights( B, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">WEIGHT", 5 ) == 0 )
	{
	  ReadSeqWeights( B, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">END", 4 ) == 0 )
	{
	  ReadEndSiteProb( B, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">BEGIN", 6 ) == 0 )
	{
	  ReadBeginSiteProb( B, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strcasecmp( token, ">GAP" ) == 0 )
	{
	  ReadGapPenalty( B, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">TRANS", 6 ) == 0 )
	{
	  ptr = GetAToken( ptr, token );  /* BT 09/14/05 - changed order of options */
	  if( ptr != NULL )
	    wt = strtod( token, &ptr2);
	  else
	    wt = DEF_TRANS_WT;
	  ptr = GetAToken( ptr, token );
	  if( ptr != NULL )
	    nUpDate = strtod( token, &ptr2);
	  else
	    nUpDate = TRUE;             /* BT 09/14/05 - changed default to true */
	  ptr = GetAToken( ptr, token );
	  if( ptr != NULL )
	    sym = strtod( token, &ptr2);
	  else
	    sym = FALSE;
	  ReadTransProbability( B, wt, nUpDate, sym, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">POS", 4 ) == 0 )
	{
	  ptr = GetAToken( ptr, token );
	  if( ptr != NULL )
	    nUpDate = strtod( token, &ptr2);
	  else
	    nUpDate = FALSE;
	  ptr = GetAToken( ptr, token );
	  if( ptr != NULL )
	    wt = strtod( token, &ptr2);
	  else
	    wt = DEF_POS_WT;
	  ReadPosMatrix( B, wt, nUpDate, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
#endif
#ifndef _NO_PHYLO_
      else if( strncasecmp( token, ">PHYLO", 6 ) == 0 ||strncasecmp( token, ">TREE", 5 ) == 0  )
	{
	  ptr = GetAToken( ptr, token );
	  if( ptr != NULL )
	    treeStartSeq =  max( atoi( token ) - 1, 0 );
	  else
	    treeStartSeq = 0;
	  ptr = GetAToken( ptr, token );
	  if( ptr != NULL )
	    treeEndSeq =  max( atoi( token ) - 1, 0 );
	  else
	    treeEndSeq = B->IP->nNumSequences - 1;
	  ptr = GetAToken( ptr, token );
	  if( ptr != NULL )
	    ReadTree( B, treeStartSeq, treeEndSeq, 1, fpt );
	  else
	    ReadTree( B, treeStartSeq, treeEndSeq, 0, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
#endif      
      else if( strncasecmp( token, ">MASK", 5 ) == 0 )
	{
	  ptr = GetAToken( ptr, token );
	  t = max( atoi( token ) - 1, 0 );
	  if( t < IP->nNumMotifTypes )	    
	    ReadMask( B, t, fpt );
	  else
	    RunOut( fpt, line, MAX_LENGTH );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">MINSITEMAP", 10 ) == 0 )
	{
	  ReadMinSiteMap( B, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">KSAMPLEMAP", 11 ) == 0 )
	{
	  ReadKSampleMap( B, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">ALIGN", 6 ) == 0 )
	{
	  ReadAlignWeight( B, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strncasecmp( token, ">NMER", 5 ) == 0 )
	{
	  ReadBkgndFreqFile( B, ptr, fpt );
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else if( strcasecmp( token, ">" ) == 0 )
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
      else if( token[0] == '!'  )
	GetALine( fpt, line, MAX_LENGTH, TRUE );
      else if( strncasecmp( token, "#INCLUDE", 8 ) == 0 )
	{
	  ptr = GetAToken( ptr, token );
	  fpt_save = fpt;
	  strncpy( cIncludeFile, token, sizeof( cIncludeFile ) );
	  fpt = fopen(token, "r");
	  if( fpt == NULL )
	    {
	      sprintf( line, "Unable to open %s\n", token );
	      p_error( "Unable to open include file." );
	    }
	  ReadPriorFile2( B, C, fpt );
	  fclose( fpt );
	  fpt = fpt_save;
	  GetALine( fpt, line, MAX_LENGTH, TRUE );
	}
      else
	{
	  NEW( msg, MAX_LENGTH, char );
	  sprintf( msg, "Error in priors file: %s", line );
	  p_error( msg );
	}
      ptr = GetAToken( line, token );           
    } 

  free(token );
  free( line );
}


void ReadPriors( Model B, Ctype C, int t, FILE* fpt )
{
  IPtype IP;
  char   *line;
  char   *token;
  char   *ptr;
  char   *ptr2;
  double pseudocnt;
  double weight[21];
  int    nNumMotifs; 
  char   aminoacids[20] = "acdefghiklmnpqrstbjo";
  int    nEnd;
  int    i;
  int    j;
  int    index;

  IP = B->IP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  nNumMotifs = NUMMOTIFS(IP->nNumMotifs[t]);    

  if( IP->is_defined[cl_d] )        
    nEnd = IP->nMaxMotifLen[t];
  else
    nEnd = IP->nMotifLen[t];
 
  for( j = 1; j <= nEnd; j++ )
    {
      if( ! GetALine( fpt, line, MAX_LENGTH, TRUE ) )
	break;
      
      ptr = line;
      ptr = GetAToken( ptr, token );
      if( strncmp( token, ">", 1 ) == 0 )
	break;
      for( i = 0; i < IP->nAlphaLen; i++ )
	{
	  weight[i] = strtod( token, &ptr2 );
	  ptr = GetAToken( ptr, token );
	}
      
      if( ptr == NULL )
	weight[IP->nAlphaLen] = IP->dPseudoCntWt;
      else
	{
	  weight[IP->nAlphaLen] = strtod(token, &ptr2);
	  if( weight[IP->nAlphaLen] == 0.0 )
	    weight[IP->nAlphaLen] = IP->dPseudoCntWt;
	}

      for(i = 0; i < IP->nAlphaLen; i++) 
	{
	  pseudocnt = weight[i]/100.0;
	  if(IP->nAlphaLen == 20)  
	    index = (int)aminoacids[i] - 97;
	  else
	    index = i;  
	  C->dPseudoCounts[t][j][index] = pseudocnt * weight[IP->nAlphaLen] *
	                                  (double)nNumMotifs;
	  if( C->dPseudoCounts[t][j][index]< MIN_PSEUDOCNT)
	    C->dPseudoCounts[t][j][index] = MIN_PSEUDOCNT;
       }
    }
  
  if( GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    {
      ptr = line;
      ptr = GetAToken( ptr, token );
      if( strcasecmp( token, ">" ) != 0 )
	{
	  sprintf( line, "Error in Priors: %s\n", token );
	  p_error( line );
	}
    }

  free(token );
  free( line );  
}


void ReadPseudo( Model B, Ctype C, int t, FILE* fpt )
{
  IPtype IP;
  char   *line;
  char   *token;
  char   *ptr;
  char   *ptr2;
  double pseudocnt;
  double weight[21];
  int    nNumMotifs; 
  char   aminoacids[20] = "acdefghiklmnpqrstbjo";
  int    nEnd;
  int    i;
  int    j;
  int    index;
  double pseudoSum;

  IP = B->IP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  nNumMotifs = NUMMOTIFS(IP->nNumMotifs[t]);    

  if( IP->is_defined[cl_d] )        
    nEnd = IP->nMaxMotifLen[t];
  else
    nEnd = IP->nMotifLen[t];
 
  for( j = 1; j <= nEnd; j++ )
    {
      if( ! GetALine( fpt, line, MAX_LENGTH, TRUE ) )
	break;
     
      ptr = line;
      ptr = GetAToken( ptr, token );
      if( strncmp( token, ">", 1 ) == 0 )
	{	  
	  sprintf( line, "Error in Pseudo %d - should be %d lines\n", 
		   t+1, nEnd );
	  p_error( line );
	}

      pseudoSum = 0;
      for( i = 0; i < IP->nAlphaLen; i++ )
	{
	  weight[i] = strtod( token, &ptr2 );
	  pseudoSum += weight[i];
	  ptr = GetAToken( ptr, token );
	}

      if( pseudoSum == 0 )
	p_error( "Pseudocounts sum to 0" );	  
      
      if( ptr == NULL )
	{
	  if( IP->is_defined[cl_Y] )	   
	    weight[IP->nAlphaLen] = (1.0 / pseudoSum) * (0.2 * IP->nNumSequences);
	  else
	    weight[IP->nAlphaLen] = IP->dPseudoCntWt;
	}	  
      else
	{
	  weight[IP->nAlphaLen] = strtod( token, &ptr2);
	  if( weight[IP->nAlphaLen] == 0.0 )
	    weight[IP->nAlphaLen] = IP->dPseudoCntWt;
	}

      for(i = 0; i < IP->nAlphaLen; i++) 
	{
	  pseudocnt = weight[i];
	  if(IP->nAlphaLen == 20)  
	    index = (int)aminoacids[i] - 97;
	  else
	    index = i;  
	  C->dPseudoCounts[t][j][index] = pseudocnt * weight[IP->nAlphaLen];
	  if( C->dPseudoCounts[t][j][index]< MIN_PSEUDOCNT)
	    C->dPseudoCounts[t][j][index] = MIN_PSEUDOCNT;
       }
    }
  
  if( GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    {
      ptr = line;
      ptr = GetAToken( ptr, token );
      if( strcasecmp( token, ">" ) != 0 )
	{
	  sprintf( line, "Error in Pseudo: %s\n", token );
	  p_error( line );
	}
    }

  set_posterior_prob( B->IP, B->C);

  free(token );
  free( line );  
}

/* BT 04/27/05 -- Use same gap for all motif types */
void ReadGapPenalty( Model B, FILE* fpt )
{
  IPtype IP;
  RPType RP;
  int    j;
  char   *ptr;
  char   *ptr2;
  char   *line;
  char   *token;
  int    len = 0;
  int    t;
  int    t1;
  double prob;
  double sum = 0.0;

  IP = B->IP;
  RP = B->RP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  while( GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    {
      ptr = line;
      while( (ptr = GetAToken( ptr, token )) )
	{
	  if( strncmp( token, ">", 1 ) == 0 )
	    break;
	  prob = strtod( token, &ptr2 );
	  if( prob < 0 )
	    p_error( "ReadGapPenalty: Gap Probabilities weight must be > 0." );
	  sum += prob;
	  for( t1 = 0; t1 < IP->nNumMotifTypes; t1++ )
	    {
	      for( t = 0; t < IP->nNumMotifTypes; t++ )
		RP->dPGap[t1][t][len] = prob;
	    }
	  len++;
	}

      if( len >  RP->nMaxGap )
	p_error( "ReadGapPenalty: Gap probabilities are longer than the maximum sequence length." );

      if( strncmp( token, ">", 1 ) == 0 )
	break;
    }

  for( j = 0; j < len; j++ )
    {
      for( t1 = 0; t1 < IP->nNumMotifTypes; t1++ )
	{
	  for( t = 0; t < IP->nNumMotifTypes; t++ )
	    RP->dPGap[t1][t][j] /= sum;
	}
    }

  for( j = len + 1; j < RP->nMaxGap; j++ )
    {
      for( t1 = 0; t1 < IP->nNumMotifTypes; t1++ )
	{
	  for( t = 0; t < IP->nNumMotifTypes; t++ )
	    RP->dPGap[t1][t][j] = 0;
	}
    }

  RP->nUseGap = TRUE;

  free( line );
  free( token );
}

/* RP->dPTrans[t][dir][t1][d] - probability of site of type t following a site of type t1 */
/* it's this order so we can speed up array access in CalcGapSum */
/* we're not using direction in transition matrix right now */
/* It's assumed that the original matrix is symmetric */
void ReadTransProbability( Model B, double dTransWt, int nUpDate, int nSym, FILE* fpt )
{
  IPtype IP;
  RPType RP;
  int    t;
  int    t1;
  char   *ptr;
  char   *line;
  char   *token;
  double tr;
  double sum;
  double totalSum = 0;

  IP = B->IP;
  RP = B->RP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      if( ! GetALine( fpt, line, MAX_LENGTH, TRUE ) )
	break;
	  
      ptr = line;
      ptr = GetAToken( ptr, token );
      if( strncmp( token, ">", 1 ) == 0 )
	break;
      sum = 0;
      for( t1 = 0; t1 < IP->nNumMotifTypes; t1++ )
	{
	  if( sscanf( token, "%lf", &tr ) != 1 )
	    p_error( "ReadTransProbability: Error in transition probability table." );
	  /*	  tr = strtod( token, &ptr2 ); */
	  RP->dInitialPTrans[t][FORWARD][t1][FORWARD] = tr;
	  RP->dInitialPTrans[t][FORWARD][t1][REVERSE] = tr;
	  RP->dInitialPTrans[t][REVERSE][t1][FORWARD] = tr;
	  RP->dInitialPTrans[t][REVERSE][t1][REVERSE] = tr;
	  
	  if( IP->RevComplement )  /* BT 2/7/2001 */
	    {
	      /* sum += 4 * tr; */
	      sum += tr;    /* Bt 01/28/03 */
	      totalSum += tr;
	    }
	  else
	    {
	      sum += tr;
	      totalSum += tr;
	    }
	  
	  ptr = GetAToken( ptr, token );		  
	}

      /* Normalize */  /* BT 12/03/2002 */
      if( sum > 0 ) 
	{
	  for( t1 = 0; t1 < IP->nNumMotifTypes; t1++ )
	    {
	      RP->dInitialPTrans[t][FORWARD][t1][FORWARD] /= sum;
	      RP->dInitialPTrans[t][FORWARD][t1][REVERSE] /= sum;
	      RP->dInitialPTrans[t][REVERSE][t1][FORWARD] /= sum;
	      RP->dInitialPTrans[t][REVERSE][t1][REVERSE] /= sum;
	    } 
	} 
    }
  
  /* Normalize */ /* BT 12/03/2002 */
  /* for( t1 = 0; t1 < IP->nNumMotifTypes; t1++ )
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
	  RP->dInitialPTrans[t][FORWARD][t1][FORWARD] /= sum;
	  RP->dInitialPTrans[t][FORWARD][t1][REVERSE] /= sum;
	  RP->dInitialPTrans[t][REVERSE][t1][FORWARD] /= sum;
	  RP->dInitialPTrans[t][REVERSE][t1][REVERSE] /= sum;
	}
	} */ 
  
  if( dTransWt <= 0.0 || dTransWt >= 1.0 )
    p_error( "ReadTransProbability: Trans weight must be between 0 and 1." );
    
  RP->dTransWt = dTransWt;
  RP->bUseTrans = TRUE;
  RP->bUpdateTrans = nUpDate;
  RP->bSymTrans = nSym;

  if( RP->bSymTrans )
    {
      for( t = 0; t < IP->nNumMotifTypes - 1; t++ )
	{
	  for( t1 = t+1; t1 < IP->nNumMotifTypes; t1++ )
	    {
	      tr = (RP->dInitialPTrans[t][FORWARD][t1][FORWARD] + 
		    RP->dInitialPTrans[t1][FORWARD][t][FORWARD]) / 2.0;
	      RP->dInitialPTrans[t][FORWARD][t1][FORWARD] = tr;
	      RP->dInitialPTrans[t1][FORWARD][t][FORWARD] = tr;
	      RP->dInitialPTrans[t][FORWARD][t1][REVERSE] = tr;
	      RP->dInitialPTrans[t1][FORWARD][t][REVERSE] = tr;
	      RP->dInitialPTrans[t][REVERSE][t1][FORWARD] = tr;
	      RP->dInitialPTrans[t1][REVERSE][t][FORWARD] = tr;
	      RP->dInitialPTrans[t][REVERSE][t1][REVERSE] = tr;
	      RP->dInitialPTrans[t1][REVERSE][t][REVERSE] = tr;
	    } 
	}	  
    }

  free( line );
  free( token );
}


void ReadEndPenalty( Model B, int t, FILE* fpt )
{
  IPtype IP;
  RPType RP;
  int    j;
  char   *ptr;
  char   *ptr2;
  char   *line;
  char   *token;

  IP = B->IP;
  RP = B->RP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  for( j = 0; j <= RP->nMaxGap; j++ )
    {
      if( ! GetALine( fpt, line, MAX_LENGTH, TRUE ) )
	break;
      
      ptr = line;
      ptr = GetAToken( ptr, token );
      if( strncmp( token, ">", 1 ) == 0 )
	break;
      RP->dPEndGap[j][t] = strtod( token, &ptr2 );
    }
   
  free( line );
  free( token );
}


void ReadFootPrintData( Model B, int nSeq, FILE* fpt )
{
  IPtype IP;
  RPType RP;
  char   *ptr = 0;
  char   *ptr2;
  char   *line;
  char   *token;
  int    len;
  int    cnt;

  IP = B->IP;
  RP = B->RP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  if( GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    ptr = GetAToken( line, token );
  else
    strcpy( token, ">" );

  len = SequenceLength( B, nSeq );
  cnt = 0;
  while( cnt < len && strcmp( token, ">" ) != 0 )
    {
      RP->dProbCons[nSeq][cnt] = strtod( token, &ptr2 );

      if( RP->dProbCons[nSeq][cnt] < FOOTPRINT_CUTOFF )   /* BT 11/10/99 */
	RP->nonCutoffSites--;

      cnt++;
      ptr = GetAToken( ptr, token );
      if( ptr == NULL )
	{
	  if( GetALine( fpt, line, MAX_LENGTH, TRUE ) )
	    ptr = GetAToken( line, token );
	  else
	    strcpy( token, ">" );
	}
    }

  IP->is_defined[cl_T] = TRUE;

  free( line );
  free( token );
}


void ReadSiteProb( Model B, FILE* fpt )
{
  IPtype IP;
  RPType RP;
  char   *ptr = 0;
  char   *ptr2;
  char   *line;
  char   *token;
  double alpha;
  double beta;
  double wt;

  IP = B->IP;
  RP = B->RP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  if( GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    ptr = GetAToken( line, token );

  if( strncmp( token, ">", 1 ) == 0 )
    return;

  alpha = strtod( token, &ptr2 );
  ptr = GetAToken( ptr, token );
  beta = strtod( token, &ptr2 );
  if( ptr != NULL )
    {
      ptr = GetAToken( ptr, token );
      wt = strtod( token, &ptr2 );      
      if( wt < 0 )
	p_error( "ReadSiteProb: Site Probabilities weight must be > 0." );
      RP->dSiteProbWt = wt;
    }

  if( alpha < 0 || beta < 0 )
    p_error( "ReadSiteProb: Negative Site Probability value." );

  if( alpha > 1.0 || beta > 1.0 )
    p_error( "ReadSiteProb: Site Probabilities must be >= 0 and <= 1.0" );

  RP->dSiteProbCons = alpha;
  RP->dSiteProbNonCons = beta;
  IP->is_defined[cl_T] = TRUE;

  free( line );
  free( token );
}


void ReadFlatGapProb( Model B, FILE* fpt )
{
  RPType RP;
  IPtype IP;
  char   *ptr;
  char   *ptr2;
  char   *line;
  char   *token;
  double dGapProb;
  int    nLen;
  int    i;
  int    j;
  int    k;

  RP = B->RP;
  IP = B->IP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  if( ! GetALine( fpt, line, MAX_LENGTH, TRUE ) )  /* BT 04/14/05 */
    return;

  ptr = GetAToken( line, token );

  if( strncmp( token, ">", 1 ) == 0 )
    return;

  nLen = atoi( token );
  if( nLen < 0 )
    p_error( "Negative gap length." );
  ptr = GetAToken( ptr, token );
  if( ptr == NULL )
    {
      dGapProb = 1.0 / (double) nLen; 
    }
  else
    {
      dGapProb = strtod( token, &ptr2 );
      if( dGapProb < 0 )
	p_error( "Negtive gap Probability value." );
    }

  nLen = min( MaxSequenceLen( B ) - 1, nLen );

  RP->bUseFlatGap = TRUE;
  RP->nFlatGapLen = nLen;
  RP->dFlatGapProb = dGapProb;
  RP->nUseGap = TRUE;

  for( i = 0; i < IP->nNumMotifTypes; i++ )
    {
      for( j = 0; j < IP->nNumMotifTypes; j++ )
	{
	  RP->nMaxGapWidth[i][j] = nLen;
	  for( k = 0; k <= nLen; k++ )
	    RP->dPGap[i][j][k] = dGapProb;
	  for( k = nLen + 1; k <= RP->nMaxGap; k++ ) 
	    RP->dPGap[i][j][k] = 0;
	}
    }

  /*  B->IP->is_defined[cl_nopt] = FALSE;
      B->IP->is_defined[cl_m] = TRUE; */

  free( line );
  free( token );
}


void ReadBlockProb( Model B, double weight, FILE* fpt )
{
  IPtype IP;
  RPType RP;
  char   *ptr = 0;
  char   *ptr2;
  char   *line;
  char   *token;
  int    k;
  double alpha;
  double  sum = 0;

  IP = B->IP;
  RP = B->RP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  if( GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    ptr = GetAToken( line, token );

  if( strncmp( token, ">", 1 ) == 0 )
    return;

  for( k = 0; k <= RP->nMaxBlocks; k++ )
    {
      if( strncmp( token, ">", 1 ) == 0 || strncmp( token, "", 1 ) == 0 )
	break;
      alpha = strtod( token, &ptr2 );
      ptr = GetAToken( ptr, token );
      if( alpha < 0 )
	p_error( "ReadBlockProb: Probability values must be  >= 0." );
      RP->dExpBlkCnt[k] = alpha;
    }

  for( k = 0; k <= RP->nMaxBlocks; k++ )
    sum += RP->dExpBlkCnt[k];

  if( sum > 0 )    /* Normalize */
    {
      for( k = 0; k <= RP->nMaxBlocks; k++ )
	{
	  RP->dExpBlkCnt[k] /= sum;
	}
    }

  RP->nUseFixedBlockProbs = TRUE;
  RP->dExpBlkWt = weight;

  for( k = 0; k <= RP->nMaxBlocks; k++ )
    {
      RP->priorSitesPerSeq[k] = RP->dExpBlkWt * RP->dPriorSeq * RP->dExpBlkCnt[k];
      if( RP->priorSitesPerSeq[k] < MIN_PSEUDOCNT )
	RP->priorSitesPerSeq[k] = MIN_PSEUDOCNT;
    }

  free( line );
  free( token );
}


void ReadSitePseudo( Model B, double weight, FILE* fpt )
{
  IPtype IP;
  RPType RP;
  char   *ptr = 0;
  char   *ptr2;
  char   *line;
  char   *token;
  int    k;
  double alpha;

  IP = B->IP;
  RP = B->RP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  if( GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    ptr = GetAToken( line, token );

  if( strncmp( token, ">", 1 ) == 0 )
    return;

  for( k = 0; k <= RP->nMaxBlocks; k++ )
    {
      if( strncmp( token, ">", 1 ) == 0 || strncmp( token, "", 1 ) == 0 )
	break;
      alpha = strtod( token, &ptr2 );
      ptr = GetAToken( ptr, token );
      if( alpha <= 0 )
	p_error( "ReadSitePseudo: Pseudocounts must be  > 0." );
      RP->priorSitesPerSeq[k] = alpha;
      if( RP->priorSitesPerSeq[k] < MIN_PSEUDOCNT )
	RP->priorSitesPerSeq[k] = MIN_PSEUDOCNT;
    }

  RP->nUseFixedBlockProbs = TRUE;
  RP->dExpBlkWt = weight;

  free( line );
  free( token );
}


void ReadSequenceBlocks( Model B, FILE* fpt )
{
  IPtype IP;
  RPType RP;
  char   *ptr = 0;
  char   *ptr2;
  char   *line;
  char   *token;
  double *alpha;
  double beta;
  int    i;
  int    nCount;
  int    t;
  double alphasum = 0.0;
  int    k;

  IP = B->IP;
  RP = B->RP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );
  NEW( alpha, IP->nNumMotifTypes, double );

  if( GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    ptr = GetAToken( line, token );

  if( strncmp( token, ">", 1 ) == 0 )
    return;

  if( B->IP->site_samp )
    p_error( "ReadSequenceBlocks: The >SEQ option is not valid for the Site Sampler." );
         
  beta = strtod( token, &ptr2 );
  if( beta < 0 )
    p_error( "ReadSequenceBlocks: Negtive value for the number of sequences." );
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      ptr = GetAToken( ptr, token );
      if( ptr == NULL )
	{
	  /*	  sprintf( line, "ReadSequenceBlocks: Number of sites for motif %d must be > 0.", t );
		  p_error( line ); */
	  alpha[t] = 0;
	}
      alpha[t] = strtod( token, &ptr2 );
      if( alpha[t] < 0 )
	{
	  sprintf( line, "ReadSequenceBlocks: Number of sites for motif %d must be >= 0.", t );
	  p_error( line );
	}      
      if( B->IP->is_defined[cl_D] && (int) alpha[t] % B->Phylo->phyloSpecies[0] != 0 )
	{
	  for( k = 1; k < B->Phylo->phyloSpecies[0]; k++ )
	    alpha[t]++;
	}
      alphasum += alpha[t];
    }

  RP->dPriorSites = alphasum;
  RP->dPriorSeq = beta;

  if( alphasum >= 0 && beta >= 0 )
    {
      for( i = 0; i <= RP->nMaxBlocks; i++ )
	if( alphasum > 0 && beta > 0 )   /* added 03/13/03 */
	  RP->dExpBlkCnt[i] = NegBinomialPrior( alphasum, beta, i );

      /* distribute total motif counts */
      for( t = 0; t <  IP->nNumMotifTypes; t++ )
	{
	  IP->nNumMotifs[t][FORWARD] = 0;
	  IP->nNumMotifs[t][REVERSE] = 0;

	  nCount = 0;
	  while( nCount < alpha[t] )
	    {
	      IP->nNumMotifs[t][FORWARD]++;
	      nCount++;
	      
	      if( IP->RevComplement &&  nCount < alpha[t] ) 
		{
		  IP->nNumMotifs[t][REVERSE]++;
		  nCount++;
		}
	    }
	  if( IP->is_defined[cl_D] )
	    {
	      if( IP->nNumMotifs[t][FORWARD] % 2 == 1 )
		{
		  IP->nNumMotifs[t][FORWARD]++;
		  IP->nNumMotifs[t][REVERSE]++;
		}
	    }
	}
    }
  
  if( GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    ptr = GetAToken( line, token );

  for( k = 0; k <= RP->nMaxBlocks; k++ )
    {
      RP->priorSitesPerSeq[k] = RP->dExpBlkWt * RP->dPriorSeq * RP->dExpBlkCnt[k];
      if( RP->priorSitesPerSeq[k] < MIN_PSEUDOCNT )
	RP->priorSitesPerSeq[k] = MIN_PSEUDOCNT;
    }
  
  free( line );
  free( token );
  free( alpha );
}


void ReadGroups( Model B, FILE* fpt )
{
  IPtype IP;
  RPType RP;
  char   *ptr;
  char   *line;
  char   *token;
  int    *groups;
  int    cnt = 0;
  int    i;
  int    j;
  int    gr;

  IP = B->IP;
  RP = B->RP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );
  NEW( groups, IP->nNumSequences, int );

  if( GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    {
      ptr = line;
      while( TRUE )
	{
	  ptr = GetAToken( ptr, token );
	  if( ptr == NULL )
	    {
	      if( ! GetALine( fpt, line, MAX_LENGTH, TRUE ) )		
		strcpy( token, ">" );
	      else
		ptr = GetAToken( line, token );
	    }
	  if( strncmp( token, ">", 1 ) == 0 )
	    break;            
	  
	  gr = atoi( token ) - 1;
	  if( gr < 0 )
	    p_error( "ReadGroups: Negative sequence number." );
	  if( gr >= IP->nNumSequences )
	    break;
	    /*	    p_error( "ReadGroups: Invalid sequence number." ); */
	  groups[cnt] = gr;
	  cnt++;
	  if( cnt >= IP->nNumSequences ) 
	    break;
	  /*	    p_error( "ReadGroups: Too many groups." ); */
	}
    }

  if( cnt > 0 )
    {
      NEW( RP->groups, IP->nNumSequences, int );
      RP->nGroupCnt = cnt;
      for( i = 0; i < cnt; i++ )
	{
	  for( j = groups[i]; j < IP->nNumSequences; j++ )
	    {
	      RP->groups[j] = i;
	    }
	}
    }

  free( line );
  free( token );
  free( groups );
}


void ReadInitPos( Model B, int seq, FILE* fpt )
{
  IPtype IP;
  RPType RP;
  int    j;
  char   *ptr;
  char   *ptr2;
  char   *line;
  char   *token;
  int    startPos;
  int    t;
  int    loc;
  int    pos;  

  IP = B->IP;
  RP = B->RP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  if( B->InitPos == NULL )
    {
      NEWP(B->InitPos, B->IP->nNumMotifTypes, PoSition );
      for( t=0; t<B->IP->nNumMotifTypes; t++)
	{
	  NEW( B->InitPos[t], B->IP->nSeqLen, PoSition );
	  IP->nNumMotifs[t][FORWARD] = 0;
	  IP->nNumMotifs[t][REVERSE] = 0;
	}
    }

  startPos = SequenceStartPos( B, seq );

  while( GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    {
      ptr = line;
      ptr = GetAToken( ptr, token );
      if( strncmp( token, ">", 1 ) == 0 )
	break;
      pos = strtod( token, &ptr2 );
      ptr = GetAToken( ptr, token );
      t = strtod( token, &ptr2 );
      ptr = GetAToken( ptr, token );

      if( t >= IP->nNumMotifTypes )
	{	  
	  sprintf( token, "Invalid INIT line - motif types must be <= 0 < %d: %s", IP->nNumMotifTypes, line );
	  p_error( token );
	}

      loc = startPos + pos - 1;
      if( loc < 0 || pos > SequenceEndPos( B, seq ) )
	{
	  sprintf( token, "Invalid INIT line - position for sequence %d  must be <= 0 <= %d: %s", 
		   seq, SequenceEndPos( B, seq ), line );
	  p_error( token );
	}

      B->InitPos[t][loc].nMotifStartPos = TRUE;
      for(j = 0; j < B->IP->nMotifLen[t]; j++) 
	B->InitPos[t][j + loc].nInMotif = TRUE;
      if( strcmp( token, "F" ) == 0 )
	{
	  B->InitPos[t][loc].RevComp = FORWARD;
	  IP->nNumMotifs[t][FORWARD]++;
	}
      else
	{
	  B->InitPos[t][loc].RevComp = REVERSE;
	  IP->nNumMotifs[t][REVERSE]++;
	}
   }
  
  free( line );
  free( token );
}


void ReadBackProb( Model B, FILE* fpt )
{
  IPtype IP;
  RPType RP;
  char   *ptr;
  char   *line;
  char   *token;
  char   *tmpstr;
  int    i;
  int    j;

  IP = B->IP;
  RP = B->RP;

  NEW( tmpstr, MAX_LENGTH, char );
  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  if( GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    ptr = GetAToken( line, token );

  if( strncmp( token, ">", 1 ) == 0 )
    return;

  IP->Datafiles->bkgnd_fpt = fopen(token, "r");

  if(IP->Datafiles->bkgnd_fpt == NULL)	
    {                
      sprintf(tmpstr, "Cannot open background composition file %s", token );  
      fprintf( stderr, "In ReadBackProb\n");
      p_error(tmpstr);
      return;
    }
  
  IP->is_defined[cl_B] = TRUE;
  NEW( IP->Datafiles->BkgndFileName, strlen(token) + 1, char);   /* BT 1/27/97 */
  strcpy( IP->Datafiles->BkgndFileName, token );
  
  free( line );
  free( token );
  free( tmpstr );

  if( B->BP == NULL )
    {
      NEW( B->BP, 1, BkgndStruct );
      NEWPP( B->BP->dBkgndProb, IP->nNumSequences, double );   
      for( i = 0; i < IP->nNumSequences; i++ )
	{
	  NEWP( B->BP->dBkgndProb[i], SequenceLength( B, i ), double );
	  for( j = 0; j < SequenceLength( B, i ); j++ )
	    NEW( B->BP->dBkgndProb[i][j], IP->nAlphaLen, double );
	}
    }
}


void ReadSpacingProb( Model B, FILE* fpt )
{
  IPtype IP;
  RPType RP;
  char   *line;
  char   *token;
  char   *tmpstr;
  FILE   *spFile;
  char   *ptr;

  IP = B->IP;
  RP = B->RP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  if( GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    ptr = GetAToken( line, token );

  if( strncmp( token, ">", 1 ) == 0 )
    return;

  spFile = fopen(token, "r");

  if(spFile == NULL)	
    {                
      NEW( tmpstr, MAX_LENGTH, char );
      sprintf(tmpstr, "Cannot open spacing probability file %s", token );  
      fprintf( stderr, "In ReadSpacingProb\n");
      p_error(tmpstr);
      return;
    }

  ReadSpacingFile( B, spFile );

  fclose( spFile );

  free( line );
  free( token );
}


void ReadSpacingFile( Model B, FILE *spFile )
{
  IPtype IP;
  RPType RP;
  char   *ptr;
  char   *ptr2;
  char   *line;
  char   *token;
  char   *tmpstr;
  int    i;
  int    j;
  int    j2;
  double *dSpacingProb;
  int    cnt;
  int    spaceCnt;

  IP = B->IP;
  RP = B->RP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  if(spFile == NULL)	
    {                
      NEW( tmpstr, MAX_LENGTH, char );
      sprintf(tmpstr, "Cannot open spacing probability file %s", token );  
      fprintf( stderr, "In ReadSpacingFile\n");
      p_error(tmpstr);
      return;
    }

  spaceCnt = 0;
  while(  GetALine( spFile, line, MAX_LENGTH, TRUE ) )   /* BT 10/5/2000 */
    {
      ptr = line;
      while( (ptr = GetAToken( ptr, token )) )
	{
	  spaceCnt++;
	}
    }

  fseek( spFile, 0, SEEK_SET );

  /*   NEW( dSpacingProb, RP->nMaxSeqLen, double ); */
  NEW( dSpacingProb, spaceCnt, double );

  cnt = 0;
  while( (cnt <  spaceCnt) && GetALine( spFile, line, MAX_LENGTH, TRUE ) )
    {
      ptr = line;
      /*       while( (cnt <  RP->nMaxSeqLen) && (ptr = GetAToken( ptr, token )) ) */
      while( (ptr = GetAToken( ptr, token )) )
	{
	  dSpacingProb[cnt] = strtod( token, &ptr2 );
	  if( *ptr2 != '\0' )
	    {
	      NEW( tmpstr, MAX_LENGTH, char );
	      sprintf( tmpstr, "ReadSpacingFile: non numerical character (%c) in spacing file.", ptr2[0] ); 
	      p_error( tmpstr );
	    }
	  cnt++;
	}
    }

  if( cnt == 0 )
    p_error( "ReadSpacingFile: Spacing file is empty" );

  for( i = 0; i < IP->nNumSequences; i++ )
    {
      for( j2 = 0, j = SequenceLength( B, i ) - 1; (j >= 0) && (cnt - j2 - 1 >= 0); j--, j2++ )
	{
	  RP->dPSpacing[i][j] = dSpacingProb[cnt - j2 - 1]; 
	}
    }
  
  RP->nUseSpacingProb = TRUE;
  
  free( line );
  free( token );
  free( dSpacingProb );
}


void ReadSpacingProb2( Model B, FILE* fpt )
{
  IPtype IP;
  RPType RP;
  char   *ptr;
  char   *ptr2;
  char   *line;
  char   *token;
  char   *tmpstr;
  int    i;
  int    j;
  int    j2;
  FILE   *spFile;
  double *dSpacingProb;
  int    cnt;

  IP = B->IP;
  RP = B->RP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  if( GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    ptr = GetAToken( line, token );

  if( strncmp( token, ">", 1 ) == 0 )
    return;

  spFile = fopen(token, "r");

  if(spFile == NULL)	
    {                
      NEW( tmpstr, MAX_LENGTH, char );
      sprintf(tmpstr, "Cannot open spacing probability file %s", token );  
      fprintf( stderr, "In ReadSpacingProb\n");
      p_error(tmpstr);
      return;
    }

  NEW( dSpacingProb, RP->nMaxSeqLen, double );

  cnt = 0;
  while( (cnt <  RP->nMaxSeqLen) && GetALine( spFile, line, MAX_LENGTH, TRUE ) )
    {
      ptr = line;
      while( (cnt <  RP->nMaxSeqLen) && (ptr = GetAToken( ptr, token )) )
	{
	  dSpacingProb[cnt] = strtod( token, &ptr2 );
	  cnt++;
	}
    }

  fclose( spFile );

  for( i = 0; i < IP->nNumSequences; i++ )
    {
      for( j2 = 0, j = SequenceLength( B, i ) - 1; (j >= 0) && (cnt - j2 - 1 >= 0); j--, j2++ )
	{
	  RP->dPSpacing[i][j] = dSpacingProb[cnt - j2 - 1]; 
	}
    }
  
  RP->nUseSpacingProb = TRUE;
  
  free( line );
  free( token );
  free( dSpacingProb );
}


void ReadBkgndComposition( Model B, Ctype C )
{
  IPtype    IP;
  int       t;
  int       i;
  int       j;
  int       n;
  int       totalMotifs;
  double    dPseudo;
  char      *ptr;
  char      *token;
  char      *line;
  char      *ptr2;
  double    *dProb;
  int       seq;
  int       seq2;
  int       pos;
  double    BGDenom; 
  int       len = 0;
  double    sum;
  char      *msg;
  PhyloType PH;
  BkgType   BP;
  double    theta;

  IP = B->IP;
  PH = B->Phylo;
  BP = B->BP;

  if(IP->Datafiles->bkgnd_fpt == NULL)	
    return;        /* BT 12/5/2000 */
  
  for( totalMotifs = 0, t = 0; t < IP->nNumMotifTypes; t++) 
    {
      totalMotifs += NUMMOTIFS(IP->nNumMotifs[t]) * IP->nMotifLen[t];
    }
  
  for( dPseudo = 0.0, n = 0; n < IP->nAlphaLen; n++ )  
    dPseudo += C->dPseudoCounts[0][BG][n];

  if( B->WT == NULL )
    {
      if(IP->nAlphaLen == 20)      
	BGDenom = (double)(C->nTotBack - IP->nNumProcessed - 
			   totalMotifs) + dPseudo;
      else
	BGDenom = (double)(C->nTotBack) + dPseudo; 
    }
  else
    BGDenom = C->dBackCnt + dPseudo;
  
  for( i = 0; i < IP->nNumSequences; i++ )
    {
      for( j = 0; j < SequenceLength( B, i ); j++ )
	{
	  for( n = 0; n < IP->nAlphaLen; n++ )
	    B->BP->dBkgndProb[i][j][n] = ((double)C->fCounts[0][BG][n] + 
					  C->dPseudoCounts[0][BG][n]) / BGDenom;
	}
    }

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );
  NEW( dProb, IP->nAlphaLen, double );
  
  while( GetALine( IP->Datafiles->bkgnd_fpt, line, MAX_LENGTH, TRUE ) )
    {
      ptr = line;
      ptr = GetAToken( ptr, token );
      seq = strtod( token, &ptr2 ) - 1;
      if( seq >= IP->nNumSequences || seq < 0 )
	{
	  NEW( msg, 256, char );
	  sprintf( msg, "ReadBkgndComposition: Invalid sequence number in background file: %d", seq + 1 );
	  p_error( msg );
	}
      ptr = GetAToken( ptr, token );
      for( i = 0; i < IP->nAlphaLen; i++ )
	{
	  ptr = GetAToken( ptr, token );
	  dProb[i] = strtod( token, &ptr2 );
	}
      ptr = GetAToken( ptr, token );
      pos = strtod( token, &ptr2 ) - 1;
      len = SequenceLength( B, seq );
      if( pos < len )
	{	
	  sum = 0;
	  for( i = 0; i < IP->nAlphaLen; i++ )
	    {
	      B->BP->dBkgndProb[seq][pos][i] = dProb[i];
	      sum += dProb[i];
	    }
	  if( fabs( sum - 1.0 ) > 0.001 )
	    {
	      NEW( msg, 256, char );
	      sprintf( msg, "ReadBkgndComposition: Seq: %d Position: %d probabilities do not sum to 1", 
		       seq + 1, pos );
	      p_error( msg );
	    }
	}
      else
	{
	  NEW( msg, 256, char );
	  sprintf( msg, "ReadBkgndComposition: Invalid position in background file: seq: %d pos: %d", 
		   seq + 1, pos );
	  p_error( msg );
	}
    }

  /* recalculate background for phylo calculations */
  if( PH->phyloTree )
    {
      for( seq = 0; seq < IP->nNumSequences; seq += SpeciesInc( B, seq ) )
	{
	  if( IsPhyloSeq( B, seq ) )
	    {
	      len = SequenceLength( B, seq );
	      
	      sum = 0.0;   /* normalize the weights */
	      for( seq2 = seq; seq2 < seq + SpeciesInc( B, seq ); seq2++ )
		sum += GetSeqWeight( B, seq2, 0 );		  
	      
	      for( j = 0; j < len; j++ )
		{
		  for( n = 0; n < IP->nAlphaLen; n++ )
		    {
		      theta = 0.0;
		      for( seq2 = seq; seq2 < seq + SpeciesInc( B, seq ); seq2++ )
			{
			  theta += GetSeqWeight( B, seq2, 0 ) * BP->dBkgndProb[seq2][j][n];
			}		      
		      theta /= sum;
		      
		      /* set all background probs for this tree to the wighted average */
		      for( seq2 = seq; seq2 < seq + SpeciesInc( B, seq ); seq2++ )
			{
			  BP->dBkgndProb[seq2][j][n] = theta;
			}
		    }
		}
	    }
	}
    }
  
  free( token );
  free( line );
  free( dProb );
  
  if( ! bkgndPseudo )
    ReCalcBkgndPseudoCnts( B, C ); 
}


void ReCalcBkgndPseudoCnts( Model B, Ctype C )
{
  IPtype     IP;
  int        t;
  int        i;
  int        j;
  double     fraction;
  int        newlen;
  int        nNumMotifs;
  int        seq;
  int        nLen;
  int        start;
  BkgType    BP;
  double     *bgCounts;

  IP = B->IP;
  BP = B->BP;

  C->dTotSumMotifPseudo = 0.0;
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      for( j = 0; j < IP->nMotifLen[t]; j++ )
	{
	  for( i = 0; i < IP->nAlphaLen; i++ )
	    C->dTotSumMotifPseudo += C->dPseudoCounts[t][j][i];
	}
    }
  
  C->dSumBGPseudo = 0.0;                    
 
  if( IP->is_defined[cl_B] )
    {
      NEW( bgCounts, IP->nAlphaLen, double );

      newlen = 0;
      for( seq = 0; seq < IP->nNumSequences; seq++ )	 
	{
	  nLen = SequenceLength( B, seq );
	  start = SequenceStartPos( B, seq );
	  for( j = 0; j < nLen; j++ )
	    {
	      for( i = 0; i < IP->nAlphaLen; i++)
		bgCounts[i] += BP->dBkgndProb[seq][j][i];
	      newlen++;
	    }
	}
      
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  nNumMotifs = NUMMOTIFS(IP->nNumMotifs[t]); 
	  for(i = 0; i < IP->nAlphaLen; i++) 
	    {
	      fraction = bgCounts[i] / (double)newlen;
	      C->dPseudoCounts[t][BG][i] = 
		fraction * IP->dPseudoCntWt * (double)nNumMotifs;
	      if(C->dPseudoCounts[t][BG][i] < MIN_PSEUDOCNT)
		C->dPseudoCounts[t][BG][i] = MIN_PSEUDOCNT;
	      if(t == 0)
		C->dSumBGPseudo += C->dPseudoCounts[t][BG][i]; 
	    }
	}

      free( bgCounts );
    }
  else
    {
      if(IP->nAlphaLen == 20)   
	newlen = IP->nSeqLen - IP->nNumProcessed;
      else
	newlen = C->nTotBack;

      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  nNumMotifs = NUMMOTIFS(IP->nNumMotifs[t]);
	  for(i = 0; i < IP->nAlphaLen; i++) 
	    {
	      fraction = (double)C->fCounts[t][BG][i] / (double)newlen;
	      C->dPseudoCounts[t][BG][i] = 
		fraction * IP->dPseudoCntWt * (double)nNumMotifs;
	      if(C->dPseudoCounts[t][BG][i] < MIN_PSEUDOCNT)
		C->dPseudoCounts[t][BG][i] = MIN_PSEUDOCNT;
	      if(t == 0)
		C->dSumBGPseudo += C->dPseudoCounts[t][BG][i]; 
	    }
	}
    }
}


void ReadSimpleBkgnd( Model B, FILE* fpt )
{
  IPtype IP;
  char   *line;
  char   *token;
  char   *ptr;
  char   *ptr2;
  double prob[21];
  int    i;
  int    j;
  int    pos;
  int    seq;
  double sum = 0.0;

  IP = B->IP;

  if( B->BP == NULL )
    {
      NEW( B->BP, 1, BkgndStruct );
      NEWPP( B->BP->dBkgndProb, IP->nNumSequences, double );   
      for( i = 0; i < IP->nNumSequences; i++ )
	{
	  NEWP( B->BP->dBkgndProb[i], SequenceLength( B, i ), double );
	  for( j = 0; j < SequenceLength( B, i ); j++ )
	    NEW( B->BP->dBkgndProb[i][j], IP->nAlphaLen, double );
	}
    }
  
  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  if( ! GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    {
      free(token );
      free( line );  
      return;
    }
      
  ptr = line;
  ptr = GetAToken( ptr, token );
  if( strncmp( token, ">", 1 ) == 0 )
    {
      free(token );
      free( line );  
      return;
    }

  for( i = 0; i < IP->nAlphaLen; i++ )
    {
      prob[i] = strtod( token, &ptr2 );
      if( prob[i] < MIN_BKGND_PROB )
	prob[i] = MIN_BKGND_PROB;
      sum += prob[i];
      ptr = GetAToken( ptr, token );
    }
  
  for( i = 0; i < IP->nAlphaLen; i++ )
    prob[i] /= sum; 
  
  for( seq = 0; seq < IP->nNumSequences; seq++ )
    {
      for( pos = 0; pos < SequenceLength( B, seq ); pos++ )
	{
	  for(i = 0; i < IP->nAlphaLen; i++) 
	    B->BP->dBkgndProb[seq][pos][i] = prob[i];
	}
    }

  IP->is_defined[cl_B] = TRUE;

  free(token );
  free( line );  
}


void ReadSeqWeights( Model B, FILE* fpt )
{
  IPtype     IP;
  char       *line;
  char       *token;
  char       *ptr;
  char       *ptr2;
  int        nSeq;
  char       *tmpstr;
  double     wt;

  IP = B->IP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  NEW( B->WT, 1, SeqWeights );
  NEW( B->WT->weight, IP->nNumSequences, double );

  for( nSeq = 0; nSeq < IP->nNumSequences; nSeq++ )
    B->WT->weight[nSeq] = 1.0;      
  
  while( GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    {
      ptr = line;
      ptr = GetAToken( ptr, token );
      if( strncmp( token, ">", 1 ) == 0 )
	break;
      nSeq = atoi( token ) - 1;
      if( nSeq < 0 || nSeq >= IP->nNumSequences )
	{                
	  NEW( tmpstr, MAX_LENGTH, char );
	  sprintf(tmpstr, "ReadSeqWeights: Invalid sequence number %s", token );  
	  p_error(tmpstr);
	}
	      
      ptr = GetAToken( ptr, token );
      wt = strtod( token, &ptr2 );      
      if( wt < 0 || strncmp( token, ">", 1 ) == 0  )
	{                
	  NEW( tmpstr, MAX_LENGTH, char );
	  sprintf(tmpstr, "ReadSeqWeights: Invalid sequence weight %s", token );  
	  p_error(tmpstr);
	}
      B->WT->weight[nSeq] = wt;
    }  

  set_counts( B );
  
  free(token );
  free( line );  
}


void ReadWeightFile( Model B )
{
  char  *tmpstr;

  if( B->IP->Datafiles->weight_fpt == NULL )
    {                
      NEW( tmpstr, MAX_LENGTH, char );
      sprintf(tmpstr, "Invalid weight file");  
      fprintf( stderr, "In ReadWeightFile\n");
      p_error(tmpstr);
      return;
    }

  ReadWeightsFromFile( B, B->IP->Datafiles->weight_fpt );    
}


void ReadWeights( Model B, FILE* fpt )
{
  char       *line;
  char       *token;
  char       *ptr;
  FILE       *wtFile;
  char       *tmpstr;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  if( GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    ptr = GetAToken( line, token );

  if( strncmp( token, ">", 1 ) == 0 )
    return;

  wtFile = fopen(token, "r");

  if(wtFile == NULL)	
    {                
      NEW( tmpstr, MAX_LENGTH, char );
      sprintf(tmpstr, "Cannot open weight file %s", token );  
      fprintf( stderr, "In ReadWeights\n");
      p_error(tmpstr);
      return;
    }

  ReadWeightsFromFile( B, wtFile );
  fclose( wtFile );

  free(token );
  free( line );  
}


void ReadWeightsFromFile( Model B, FILE* wtFile )
{
  IPtype     IP;
  char       *line;
  char       *token;
  char       *ptr;
  char       *ptr2;
  SeqWeights *wts;
  int        nSeqTypes;
  int        cnt;
  int        wtId;
  char       desc[3];
  int        seq;

  IP = B->IP;
  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  NEW( wts, 1, SeqWeights );

  GetALine( wtFile, line, MAX_LENGTH, TRUE );  
  ptr = GetAToken( line, token );
  nSeqTypes = atoi( token );
  cnt = 0;
  while( (ptr = GetAToken( ptr, token )) && (cnt < nSeqTypes) )
    {
      wts->seqIdStr[2*cnt] = toupper( token[0] );
      wts->seqIdStr[2*cnt + 1] = toupper( token[1] );
      cnt++;
    }
  
  wts->nSeqTypes = nSeqTypes;
  wts->nTableSize = 1 << nSeqTypes;
  NEW( wts->weightTable, wts->nTableSize, WeightStruct );
  NEW( wts->seqCodes, IP->nNumSequences, int );

  while( GetALine( wtFile, line, MAX_LENGTH, TRUE ) )
    {
      ptr = line;
      ptr = GetAToken( ptr, token );
      wtId = atoi( token );
      wts->weightTable[wtId].code = wtId;
      ptr = GetAToken( ptr, token );
      strcpy( wts->weightTable[wtId].id, token );

      NEW( wts->weightTable[wtId].weight, nSeqTypes, double );
            
      cnt = 0;
      while( (ptr = GetAToken( ptr, token )) && (cnt < nSeqTypes) )
	{
	  wts->weightTable[wtId].weight[cnt] = strtod( token, &ptr2 );
	  cnt++;
	}
    }

  desc[2] = '\0';
  ptr = wts->seqIdStr;
  for( seq = 0; seq < B->IP->nNumSequences; seq++ )
    {
      desc[0] = B->IP->fastaHeader[seq][0];
      desc[1] = B->IP->fastaHeader[seq][1];

      ptr2 = strstr( wts->seqIdStr, desc );
      if( ptr2 != NULL )
	wts->seqCodes[seq] = ((ptr2 - ptr) / 2);
    }
  
  B->WT = wts;

  free(token );
  free( line );  
}


void ReadEndSiteProb( Model B, FILE* fpt )
{
  IPtype IP;
  RPType RP;
  char   *ptr = 0;
  char   *line;
  char   *token;
  int    t;
  double tr;
  double sum = 0;

  IP = B->IP;
  RP = B->RP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  if( GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    ptr = GetAToken( line, token );

  if( strncmp( token, ">", 1 ) == 0 )
    return;

  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      if( sscanf( token, "%lf", &tr ) != 1 )
	p_error( "ReadEndSite: Error in end site probability." );
      RP->dInitialEndSiteProb[t] = tr;
      sum += tr;
      ptr = GetAToken( ptr, token );		  
    }

  if( sum > 0 )
    {
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  RP->dInitialEndSiteProb[t] /= sum;
	}
    }
  else
    p_error( "ReadEndSite:Sum of probabilities = 0." );

  free( line );
  free( token );
}


void ReadBeginSiteProb( Model B, FILE* fpt )
{
  IPtype IP;
  RPType RP;
  char   *ptr = 0;
  char   *line;
  char   *token;
  int    t;
  double tr;
  double sum = 0;

  IP = B->IP;
  RP = B->RP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  if( GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    ptr = GetAToken( line, token );

  if( strncmp( token, ">", 1 ) == 0 )
    return;

  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      if( sscanf( token, "%lf", &tr ) != 1 )
	p_error( "ReadBeginSite: Error in starting site probability." );
      RP->dInitialBeginSiteProb[t] = tr;
      sum += tr;
      ptr = GetAToken( ptr, token );		  
    }

  if( sum > 0 )
    {
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  RP->dInitialBeginSiteProb[t] /= sum;
	}
    }
  else
    p_error( "ReadBeginSite:Sum of probabilities = 0." );

  free( line );
  free( token );
}


/* position matrix POS[i][j] - prob that position i in module is motif type j
   rows - number of elements in module
   cols - number of motif types */

void ReadPosMatrix( Model B, double dTransWt, int nUpDate, FILE* fpt )
{
  IPtype IP;
  RPType RP;
  char   *ptr = 0;
  char   *line;
  char   *token;
  int    t;
  int    i;
  double sum;
  double tr;

  IP = B->IP;
  RP = B->RP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  RP->nModuleSites = 0;
  for( i = 0; i < RP->nMaxBlocks; i++ )
    {  
      if( GetALine( fpt, line, MAX_LENGTH, TRUE ) )
	ptr = GetAToken( line, token );
      
      if( strncmp( token, ">", 1 ) == 0 )
	return;

      sum = 0;
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  if( sscanf( token, "%lf", &tr ) != 1 )
	    p_error( "ReadPosMatrix: Error in position probability." );
	  if( tr < 0 )
	    p_error( "ReadPosMatrix: Position probability must be >= 0." );
	  RP->dInitialPosMatrix[i][t] = tr;
	  sum += tr;
	  ptr = GetAToken( ptr, token );		  
	}
      RP->nModuleSites++;

      if( sum > 0 )
	{
	  for( t = 0; t < IP->nNumMotifTypes; t++ )
	    {
	      RP->dInitialPosMatrix[i][t] /= sum;
	    }
	}
    }
  
  RP->bUsePosMatrix = TRUE;
  RP->dPosWt = dTransWt;
  RP->bUpdatePosMatrix = nUpDate;

  free( line );
  free( token );
}

void ReadMask( Model B, int t, FILE* fpt )
{
  IPtype IP;
  RPType RP;
  int    j;
  char   *ptr;
  char   *line;
  char   *token;
  int    code;
  char   *msg;
  int    ones;

  IP = B->IP;
  RP = B->RP;

  if( IP->is_defined[cl_F] )
    {
      p_error( "ReadMask: Remove -F option if specifying a fragmentation mask." );
    }

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  if( ! GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    return;
      
  ptr = line;
  ptr = GetAToken( ptr, token );
  if( strncmp( token, ">", 1 ) == 0 )
    return;

  if( strlen( token ) != B->F->nMaxLen[t] )
    {
      NEW( msg, MAX_LENGTH, char );
      sprintf( msg, "ReadMask: Invalid mask length. Should be %d", B->F->nMaxLen[t] );
      p_error( msg );
    }	

  for( ones = 0, j = 0; j < B->F->nMaxLen[t]; j++ )
    {
      if( token[j] == '\0' )
	{
	  NEW( msg, MAX_LENGTH, char );
	  sprintf( msg, "ReadMask: Invalid mask length. Should be %d", B->F->nMaxLen[t] );
	  p_error( msg );
	}	
      code = token[j] - '0';
      if( code < 0 || code > 3 )
	p_error( "ReadMask: Mask codes must be >= 0 and <= 3." );
      if( code == 1 )
	ones++;
      B->F->fragInitMask[t][j] = code;
    }

  if( ones != B->IP->nMotifLen[t] )
    p_error( "ReadMask: The number of ON columns (1's) must equal the number of conserved positions." );
    
  initMask( B );

  free( line );
  free( token );
}


void ReadComment( Model B, FILE* fpt )
{
  IPtype IP;
  char   *line;
  char   *buffer;
  int    bufSize;
  int    done;

  IP = B->IP;

  NEW( line, MAX_LENGTH, char );
  NEW( buffer, MAX_LENGTH, char );

  bufSize = MAX_LENGTH;

  if( ! GetALine( fpt, line, MAX_LENGTH, FALSE ) )
    {
      free( line );  
      free( buffer );
      return;
    }
      
  done = FALSE;
  while( ! done )
    {
      if( strncmp( line, ">", 1 ) == 0 )
	{
	  fseek( fpt, -(strlen( line ) + 1), SEEK_CUR );
	  done = TRUE;
	}
      else
	{
	  if( strlen( line ) + 1 + strlen( buffer ) >= bufSize )
	    {
	      bufSize += MAX_LENGTH;
	      buffer = realloc( buffer, bufSize );
	    }
	  strcat( buffer, line );
	  strcat( buffer, "\n" );      
	  
	  if( ! GetALine( fpt, line, MAX_LENGTH, FALSE ) )
	    {
	      done = TRUE;
	    }
	}
    }

  IP->comment = buffer;

  free( line );  
}


void ReadGibbsName( Model B, FILE* fpt )
{
  IPtype IP;
  char   *line;
  char   *buffer;
  int    bufSize;

  IP = B->IP;

  NEW( line, MAX_LENGTH, char );
  NEW( buffer, MAX_LENGTH, char );

  bufSize = MAX_LENGTH;

  if( ! GetALine( fpt, line, MAX_LENGTH, FALSE ) )
    {
      free( line );  
      free( buffer );
      return;
    }
      
  if( strncmp( line, ">", 1 ) == 0 )
    fseek( fpt, -(strlen( line ) + 1), SEEK_CUR );
  else
    {
      if( strlen( line ) + 1 + bufSize >= bufSize )
	{
	  bufSize += MAX_LENGTH;
	  buffer = realloc( buffer, bufSize );
	}
      strcpy( buffer, line );	  
    }
  
  IP->programName = buffer;

  free( line );  
}


void GetBkgndPseudo( Model B, Ctype C, int t, FILE* fpt )
{
  IPtype IP;
  char   *line;
  char   *token;
  char   *ptr;
  char   *ptr2;
  double prob;
  int    i;

  IP = B->IP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  if( ! GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    {
      free(token );
      free( line );  
      return;
    }
      
  ptr = line;
  ptr = GetAToken( ptr, token );
  if( strncmp( token, ">", 1 ) == 0 )
    {
      free(token );
      free( line );  
      return;
    }

  for( i = 0; i < IP->nAlphaLen; i++ )
    {
      prob = strtod( token, &ptr2 );
      if( prob < MIN_BKGND_PROB )
	prob = MIN_BKGND_PROB;

      C->dPseudoCounts[t][BG][i] = prob;
      ptr = GetAToken( ptr, token );
    }
  
  bkgndPseudo = TRUE;

  free(token );
  free( line );  
}


void ReadMinSiteMap( Model B, FILE* fpt )
{
  RPType RP;
  char   *ptr;
  char   *ptr2;
  char   *line;
  char   *token;

  RP = B->RP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  if( ! GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    return;
      
  ptr = line;
  ptr = GetAToken( ptr, token );
  if( strncmp( token, ">", 1 ) == 0 )
    return;

  RP->dMinSiteMap = strtod( token, &ptr2 );
  ptr = GetAToken( ptr, token );
  if( strncmp( token, ">", 1 ) == 0 )
    return;
  RP->nMinSiteMotifs = atoi( token );

  free( line );
  free( token );
}


void ReadKSampleMap( Model B, FILE* fpt )
{
  RPType RP;
  char   *ptr;
  char   *ptr2;
  char   *line;
  char   *token;

  RP = B->RP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  if( ! GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    return;
      
  ptr = line;
  ptr = GetAToken( ptr, token );
  if( strncmp( token, ">", 1 ) == 0 )
    return;

  RP->dKSampleMap = strtod( token, &ptr2 );

  free( line );
  free( token );
}


void ReadAlignWeight( Model B, FILE* fpt )
{
  IPtype IP;
  RPType RP;
  ALType AP;
  char   *ptr = 0;
  char   *ptr2;
  char   *line;
  char   *token;
  int    k;
  double wt;

  IP = B->IP;
  RP = B->RP;
  AP = B->AP;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  if( GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    ptr = GetAToken( line, token );

  if( strncmp( token, ">", 1 ) == 0 )
    return;

  for( k = 0; k <= RP->nMaxBlocks; k++ )
    {
      if( strncmp( token, ">", 1 ) == 0 || strncmp( token, "", 1 ) == 0 )
	break;
      wt = strtod( token, &ptr2 );
      if( wt < 0 )
	p_error( "ReadAlignWeight: Weight values must be >= 0." );
      AP->dAlignWt[k] = wt;
      ptr = GetAToken( ptr, token );
    }

  free( line );
  free( token );
}


void ReadTree( Model B, int treeStartSeq, int treeEndSeq, int speciesSample, FILE* fpt )
{
  IPtype    IP;
  PhyloType PH;
  char      *line;
  char      *buffer;
  int       bufSize;
  int       done;
  PhyloTree tree;
  int       seq;
  int       length;

  IP = B->IP;
  PH = B->Phylo;

  if( IP->is_defined[cl_D] )
    p_error( "ReadTree: Don't use -D with phylogenetic trees." );

  if( ! IP->is_defined[cl_E] )
    p_error( "ReadTree: -E option must be used with phylognenetic trees." );

  if( PH->treeCount == 0 )
    {
      NEW( PH->phyloTree, PHYLO_ALLOC_COUNT, PhyloTree );
      NEW( PH->phyloSpecies, PHYLO_ALLOC_COUNT, int );
      NEW( PH->phyloSeq, PHYLO_ALLOC_COUNT, int );
      NEW( PH->phyloIndex, IP->nNumSequences, int );
      NEW( PH->phyloSpeciesSample, PHYLO_ALLOC_COUNT, int );
      NEW( PH->phyloTreeStartSeq, PHYLO_ALLOC_COUNT, int );
      for( seq = 0; seq < IP->nNumSequences; seq++ )
	PH->phyloIndex[seq] = -1;
      PH->phyloAllocCnt += PHYLO_ALLOC_COUNT;
    }
  else if( PH->treeCount % PHYLO_ALLOC_COUNT == 0 )
    {
      PH->phyloTree = (PhyloTree *) realloc( PH->phyloTree, 
					     (PH->treeCount + PHYLO_ALLOC_COUNT) * sizeof( PhyloTree * ) );
      PH->phyloSpecies = (int *) realloc( PH->phyloSpecies,
					  (PH->treeCount + PHYLO_ALLOC_COUNT) * sizeof( int * ) );
      PH->phyloSeq = (int *) realloc( PH->phyloSeq,
				      (PH->treeCount + PHYLO_ALLOC_COUNT) * sizeof( int * ) );
      PH->phyloSpeciesSample = (int *) realloc( PH->phyloSpeciesSample,
						(PH->treeCount + PHYLO_ALLOC_COUNT) * sizeof( int * ) );
      PH->phyloTreeStartSeq = (int *) realloc( PH->phyloTreeStartSeq,
					       (PH->treeCount + PHYLO_ALLOC_COUNT) * sizeof( int * ) );
      PH->phyloAllocCnt += PHYLO_ALLOC_COUNT;
    }

  PH->treeCount++;

  NEW( line, MAX_LENGTH, char );
  NEW( buffer, MAX_LENGTH, char );

  bufSize = MAX_LENGTH;

  if( ! GetALine( fpt, line, MAX_LENGTH, FALSE ) )
    {
      free( line );  
      free( buffer );
      return;
    }
      
  done = FALSE;
  while( ! done )
    {
      if( strncmp( line, ">", 1 ) == 0 )
	{
	  fseek( fpt, -(strlen( line ) + 1), SEEK_CUR );
	  done = TRUE;
	}
      else
	{
	  if( strlen( line ) + 1 + strlen( buffer ) >= bufSize )
	    {
	      bufSize += MAX_LENGTH;
	      buffer = realloc( buffer, bufSize );
	    }
	  strcat( buffer, line );
	  strcat( buffer, "\n" );      
	  
	  if( ! GetALine( fpt, line, MAX_LENGTH, FALSE ) )
	    {
	      done = TRUE;
	    }
	}
    }

  if( ! IP->is_defined[cl_bayes] )
    PH->bCalcPhylo = TRUE;

  PH->phyloSpecies[PH->treeCount-1] = ReadPhylTree( B, buffer, &tree );
  PH->phyloTree[PH->treeCount-1] = tree;
  PH->phyloSeq[PH->treeCount-1] = treeEndSeq - treeStartSeq + 1;
  PH->phyloSpeciesSample[PH->treeCount-1] = speciesSample;
  PH->phyloTreeStartSeq[PH->treeCount-1] = treeStartSeq;
  if( speciesSample )
    PH->phyloSpeciesCount += treeEndSeq - treeStartSeq + 1;
  else
    PH->phyloSpeciesCount++;
  
  for( seq = treeStartSeq; seq <= treeEndSeq; seq++ )
    {
      if( PH->phyloIndex[seq] != -1 )
	{
	  sprintf( buffer, "ReadTree: The sequence range for tree %d overlaps the range for another tree.", PH->treeCount );      
	  p_error( buffer );
	}

      PH->phyloIndex[seq] = PH->treeCount - 1;
    }

  length = SequenceLength( B, treeStartSeq );
  for( seq = treeStartSeq + 1; seq <= treeEndSeq; seq++ )
    {
      if( SequenceLength( B, seq ) != length )
	{
	  sprintf( buffer, "ReadTree: %d %d - All sequences in a MASS must be aligned and have the same length.", 
		   treeStartSeq + 1, treeEndSeq + 1 );      
	  p_error( buffer );
	}
    }

  if( B->Phylo->phyloSeq[PH->treeCount-1] % B->Phylo->phyloSpecies[PH->treeCount-1] != 0 )
    {
      sprintf( buffer, "ReadTree: The sequence count is not an even multiple of the species in tree %d", PH->treeCount );      
      p_error( buffer );
    }
  
  if( treeEndSeq > PH->maxPhyloSeq )
    PH->maxPhyloSeq = treeEndSeq;
    
  free( line );  
  free( buffer );  
}


void ReadBkgndFreqFile( Model B, char *seqStr, FILE *fpt )
{
  IPtype    IP;
  int       *seqList = NULL;
  char      *token;
  char      *fileName;
  char      *ptr;
  int       seq;
  double    *nMerArray;
  int       size;
  int       i;

  IP = B->IP;

  if( ! IP->is_defined[cl_E] )
    p_error( "ReadBkgndFreqFile: -E option must be used with nMer background model." );

  NEW( fileName, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );
  NEW( seqList, IP->nNumSequences, int );

  if( ! GetALine( fpt, fileName, MAX_LENGTH, FALSE ) )
    p_error( "ReadBkgndFreqFile: missing file name." );

  if( ! FileExists( fileName ) )
    {
      sprintf( token, "ReadBkgndFreqFile: File %s could not be found.", fileName );
      p_error( token );
    }
	     
  if( strlen( seqStr ) == 0 )
    {
      for( i = 0; i < IP->nNumSequences; i++ )
	seqList[i] = 1;	
    }
  else
    {
      ptr = seqStr;
      while( (ptr = GetAToken( ptr, token )) != NULL )
	{
	  seq = atoi( token );
	  if( seq < 1 || seq > IP->nNumSequences )
	    p_error( "ReadBkgndFreqFile: invalid sequence number." );
	  seqList[seq-1] = 1;
	}
    }

  if( ! IP->freqSeqs )
    {
      NEWP( IP->freqSeqs, IP->nNumSequences, double );
      NEW( IP->nMerSize, IP->nNumSequences, int );
    }

  ReadFreqBackground( B, fileName, &nMerArray, &size );

  for( i = 0; i < IP->nNumSequences; i++ )
    {
      if( seqList[i] )
	{
	  IP->freqSeqs[i] = nMerArray;
	  IP->nMerSize[i] = size;
	}
    }

  IP->is_defined[cl_freq_background] = TRUE;

  free( token );
  free( seqList );
  free( fileName );
}


short GetALine( FILE *fpt, char *line, int maxLength, int skipBlankLines )
{
  char ch;
  int  pos = 0;

  line[0] = '\0';

  ch = fgetc( fpt );
  if( skipBlankLines )
    {
      while( (ch != (char) EOF) && ((ch == '\n') || (ch == '\r')) )     /* skip blank lines */ /* BT 04/21/05 */
	ch = fgetc( fpt );
    }
  
  if( ch == (char) EOF ) /* BT 12/14/99 */
    return FALSE;

  while( (ch !=(char)  EOF) && (ch != '\n') && (pos < maxLength ) ) /* BT 12/14/99 */
    {
      if( ch != '\r' )    /* BT 04/21/05 */
	line[pos++] = ch;
      ch = fgetc( fpt );
    }

  line[pos] = '\0';
    
  return TRUE;
}



char *GetAToken( char *line, char *token )
{
  int lPos = 0;
  int tPos = 0;

  token[0] = '\0';

  if( line == NULL )
    return NULL;

  if( line[lPos] == '\0' )
    return NULL;

  while( line[lPos] != '\0' && strchr( delims, line[lPos] ) != NULL )
    lPos++;

  if( line[lPos] == '\0' )
    return NULL;

  token[tPos++] = line[lPos++];
  while( line[lPos] != '\0' && 
	 strchr( delims, line[lPos] ) == NULL  &&
	 line[lPos] != '!' )
    token[tPos++] = line[lPos++];

  token[tPos] = '\0';

  return &line[lPos];
}


void RunOut( FILE *fpt, char *line, int maxLength )
{
  char   *ptr;
  char   *token;

  NEW( line, MAX_LENGTH, char );
  NEW( token, MAX_LENGTH, char );

  if( ! GetALine( fpt, line, MAX_LENGTH, TRUE ) )
    {
      free( line );
      free( token );
      return;
    }
      
  ptr = line;
  ptr = GetAToken( ptr, token );
  while( strncmp( token, ">", 1 ) != 0 )
    {
      if( ! GetALine( fpt, line, MAX_LENGTH, TRUE ) )
	break;
      ptr = line;
      ptr = GetAToken( ptr, token );
    }

  free( line );
  free( token );
}


double Poisson( int n, double mean )
{
  return exp( n * log( mean ) - mean - ln_gamma( n + 1.0 ) );
}


double NegBinomialPrior( double alpha, double beta, int n )
{
  return( BiCoef( alpha + n - 1, alpha - 1 ) * 
	  pow( beta / (beta + 1), alpha ) * 
	  pow(1 / (beta + 1), n) );
}


void CheckNmers( Model B )
{
  int  i;
  char *msg;
 
  for( i = 0; i < B->Phylo->maxPhyloSeq; i++ )
    {
      if( ! B->IP->freqSeqs[i] )
	{
	  NEW( msg, 256, char );
	  sprintf( msg, "You must supply a frequency file for sequences 1 through %d",  B->Phylo->maxPhyloSeq );
	  p_error( msg );
	}
    }
}
