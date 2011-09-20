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
#include "impsamp.h"

typedef struct
{
  int  pos;
  int  rev;
  int  seq;
  int  seqPos;
  int  seqPair;
  int  seqLen;
} MotifStruct;

void ImpGetBkgndCountsPos( Model B, int seq, int pos, double *bgCounts );
int GetMotifList( Model B, int t, MotifStruct **motifs );
void GetPhyloCounts( Model B, int t, int num, MotifStruct **motifs, double **counts );


double ImpSampleMotifMap( Model B, int t )
{
  IPtype      IP;
  Ctype       C;
  PhyloType   PH;
  char        *R;
  int         i;
  int         j;
  double      prob = 0.0;
  double      felsProb;
  double      dirichletProb;
  double      ratio;
  double      sum;
  double      sumSq;
  double      avg;
  double      var;
  int         c;
  int         p;
  double      ***fgTheta;
  MotifStruct **motifs;
  double      **motifCounts;
  double      pr;
  int         phyloMotifs;
  int         impSamples = IMP_SAMPLES;

  IP = B->IP;
  C = B->C;
  PH = B->Phylo;
  R = (*(B)->Seq->R);

  if( PH->phyloSpeciesSample )
    impSamples /= 10;
    
  NEWP( motifs, IP->nMotifLen[t], MotifStruct );
  for( i = 0; i < IP->nMotifLen[t]; i++ )
    NEW( motifs[i], NUMMOTIFS( IP->nNumMotifs[t] ), MotifStruct );  
  phyloMotifs = GetMotifList( B, t, motifs );

  if( phyloMotifs == 0 )
    {
      FREEP( motifs, IP->nMotifLen[t] );
      return 0.0;
    } 

  NEWP( motifCounts, IP->nMotifLen[t] + 1, double );
  for( i = 0; i <= IP->nMotifLen[t]; i++ )
    NEW( motifCounts[i], IP->nAlphaLen, double );
  GetPhyloCounts( B, t, phyloMotifs, motifs, motifCounts );

  NEWPP( fgTheta, impSamples, double );
  for( c = 0; c < impSamples; c++ )
    {
      NEWP( fgTheta[c], IP->nMotifLen[t], double );
      for( j = 0;  j < IP->nMotifLen[t]; j++ )
	NEW( fgTheta[c][j], IP->nAlphaLen, double );
      SampleAModel( B, t, fgTheta[c], motifCounts );
    }

  for( j = 0; j < IP->nMotifLen[t]; j++ )
    {
      sum = 0.0;
      sumSq = 0.0;
      for( c = 0; c < impSamples; c++ )
	{
	  felsProb = 0.0;
	  for( i = 0; i < phyloMotifs; i++ )
	    {	      
	      CalcSubsMatrix( B, PH->phyloTree[PH->phyloIndex[motifs[j][i].seq]], fgTheta[c][j] );
	      if( PH->phyloSpeciesSample[PH->phyloIndex[motifs[j][i].seq]] )
		pr = CalcFelsSeqProb( B, PH->phyloTree[PH->phyloIndex[motifs[j][i].seq]], 
				      fgTheta[c][j], motifs[j][i].pos, motifs[j][i].rev );
	      else
		pr = CalcFelsProb( B, PH->phyloTree[PH->phyloIndex[motifs[j][i].seq]], 
				   fgTheta[c][j], motifs[j][i].pos, motifs[j][i].rev );
	      felsProb += log( pr );
#ifdef _DEBUG_
	      if( ! finite( felsProb ) )
		  p_internal_error("**** non finite result ****");
#endif
	    }

	  dirichletProb = 0.0;
	  for( p = 0; p < IP->nAlphaLen; p++ )
	    {
	      if( fgTheta[c][j][p] > 0 )
		dirichletProb += log( fgTheta[c][j][p] ) * motifCounts[j+1][p];
	    }

	  ratio = exp( felsProb - dirichletProb );

#ifdef _DEBUG_
      if( ! finite( ratio ) )
	p_error( "**** Non finite value ****" );
#endif 
	  sum += ratio;
	  sumSq += ratio * ratio;
	}
      
      avg = sum / impSamples;
      var = (sumSq / impSamples) - avg * avg;
            
      prob += log( avg );
    }

  FREEP( motifs, IP->nMotifLen[t] );
  FREEP( motifCounts, IP->nMotifLen[t] + 1 );
  FREEPP( fgTheta, impSamples, IP->nMotifLen[t] );
  
  return( prob );
}


double ImpCalcBkgndMap( Model B )
{
  IPtype       IP;
  PhyloType    PH;
  MotifStruct  **motifs;
  int          t;
  int          i;
  int          j;
  double       bkgndProb;
  int          seqPos;
  int          seq;
  int          phyloMotifs;
  
  IP = B->IP;
  PH = B->Phylo;

  bkgndProb = PH->phyloNull;

  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      NEWP( motifs, IP->nMotifLen[t], MotifStruct );
      for( i = 0; i < IP->nMotifLen[t]; i++ )
	NEW( motifs[i], NUMMOTIFS( IP->nNumMotifs[t] ), MotifStruct );
      
      phyloMotifs = GetMotifList( B, t, motifs );
      for( j = 0; j < IP->nMotifLen[t]; j++ )
	{
	  for( i = 0; i < phyloMotifs; i++ )
	    {
	      seq = motifs[j][i].seqPair;
	      seqPos = motifs[j][i].seqPos;
	      bkgndProb -= PH->phyloBkgnd[seq][seqPos];
#ifdef _DEBUG_
	      if( ! finite( bkgndProb ) )
		p_internal_error( "**** non-finite value ****" );
#endif
	    }
	}

      FREEP( motifs, IP->nMotifLen[t] );
    }

  return( bkgndProb );
}


double ImpSampleNullMap( Model B )
{
  double       *bgTheta;
  IPtype       IP;
  RPType       RP;
  Ctype        C;
  Ftype        F;
  PhyloType    PH;
  int          seq;
  int          nLen;
  int          pos;
  int          i;
  double       prob;
  double       felsProb;
  double       dirichletProb;
  double       ratio;
  double       sum;
  double       sumSq;
  double       avg;
  double       var;
  double       *pseudoCounts;
  double       *bgCounts;
  double       *counts;
  int          c;
  int          t;
  int          count;
  int          start;
  double       scale = 0;
  double       ***bgCountsPos;
  
  IP = B->IP;
  RP = B->RP;
  C = B->C;
  F = B->F;
  PH = B->Phylo;
  
  PH->phyloNull = 0.0;

  NEWP( PH->phyloBkgnd, PH->phyloSpeciesCount, double );
  NEWPP( bgCountsPos, PH->phyloSpeciesCount, double );
  for( count = 0, seq = 0; seq <= PH->maxPhyloSeq; seq += SeqInc( B, seq ), count++ )
    {
      nLen = SequenceLength( B, seq );

      NEW( PH->phyloBkgnd[count], nLen, double );
      NEWP( bgCountsPos[count], nLen, double );

      for( pos = 0; pos < nLen; pos++ )
	{
	  NEW( bgCountsPos[count][pos], IP->nAlphaLen, double );
	  ImpGetBkgndCountsPos( B, seq, pos, bgCountsPos[count][pos] );	
	}
    }
  
  NEW( bgCounts, IP->nAlphaLen, double );
  NEW( counts, IP->nAlphaLen, double );
  NEW( pseudoCounts, IP->nAlphaLen, double );
  
  NEW( bgTheta, IP->nAlphaLen, double );
   
  ImpGetBkgndCounts( B, bgCounts );
  for( i=0; i < IP->nAlphaLen; i++ ) 
    {
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  pseudoCounts[i] += C->dPseudoCounts[t][BG][i];
	}
      counts[i] = bgCounts[i] + pseudoCounts[i];
    }

  sum = 0.0;
  sumSq = 0.0;
  for( c = 0; c < IMP_SAMPLES_BKGND; c++ )
    {
      rdirichlet( B, IP->nAlphaLen, counts, bgTheta );

      felsProb = 0.0;
      for( count = 0, seq = 0; seq <= PH->maxPhyloSeq; seq += SeqInc( B, seq ), count++ )
	{
	  CalcSubsMatrix( B, PH->phyloTree[PH->phyloIndex[seq]], bgTheta );
	  nLen = SequenceLength( B, seq );
	  start = SequenceStartPos( B, seq );
	  for( pos = 0; pos < nLen; pos++ )
	    {
	      if( PH->phyloSpeciesSample[PH->phyloIndex[seq]] || PossibleAlignedBkgndPos( B, pos + start, nLen, seq ) )
		{
		  if( PH->phyloSpeciesSample[PH->phyloIndex[seq]] )
		    prob = CalcFelsSeqProb( B, PH->phyloTree[PH->phyloIndex[seq]], bgTheta, 
					    start + pos, FALSE );
		  else
		    prob = CalcFelsProb( B, PH->phyloTree[PH->phyloIndex[seq]], bgTheta, 
					 start + pos, FALSE );
		  felsProb += log( prob );
#ifdef _DEBUG_
		  if( ! finite( felsProb ) )
		    p_internal_error( "**** non-finite value ****" );
#endif

		  dirichletProb = 1.0;
		  for( i = 0; i < IP->nAlphaLen; i++ )
		    dirichletProb *= pow( bgTheta[i], bgCountsPos[count][pos][i] );
	  
		  PH->phyloBkgnd[count][pos] += prob / dirichletProb;
		}
	    }
	}
		  
      dirichletProb = 0.0;
      for( i=0; i < IP->nAlphaLen; i++ ) 
	{
	  dirichletProb += log( bgTheta[i] ) * bgCounts[i];
	}

      if( c == 0 )
	scale = felsProb - dirichletProb;

      ratio = exp( felsProb - dirichletProb - scale ); 

#ifdef _DEBUG_
      if( ! finite( ratio ) )
	p_internal_error( "**** non-finite value ****" );
#endif
      sum += ratio;
      sumSq += ratio * ratio;
    }
	      
  avg = sum / IMP_SAMPLES_BKGND;
  var = sqrt( fabs( ((sumSq / IMP_SAMPLES_BKGND) - (avg * avg)) ) / IMP_SAMPLES_BKGND );
  
  for( count = 0, seq = 0; seq <= PH->maxPhyloSeq; seq += SeqInc( B, seq ), count++ )
    {
      nLen = SequenceLength( B, seq );
      for( pos = 0; pos < nLen; pos++ )
	{
	  PH->phyloBkgnd[count][pos] = log( PH->phyloBkgnd[count][pos] ) - log( IMP_SAMPLES_BKGND );
	}
    }

  free( counts );
  free( bgCounts );
  free( pseudoCounts );
  free( bgTheta );

  for( count = 0, seq = 0; seq < PH->maxPhyloSeq; seq += SeqInc( B, seq ), count++ )
    {
      nLen = SequenceLength( B, seq );
      for( pos = 0; pos < nLen; pos++ )
	free( bgCountsPos[count][pos] );
      free(  bgCountsPos[count] );
    }
  free( bgCountsPos );

  return log( avg ) + scale;
}


/* return an array of K random variates from a Dirichlet distribution of order K-1. */
/* uses the method of sampling from a gamma distribution                            */
/* K - number of samples to return                                                  */
/* alpha - dirichlet parameters (counts + pseudocounts)                             */
/* theta - returned dirichlet samples                                                */
void rdirichlet( Model B, size_t K, double *alpha, double *theta )
{
  double *x;
  int    i;
  double sum = 0.0;

  NEW( x, K, double );
  for( i = 0; i < K; i++ )
    {
      if( alpha[i] < 0 )
	p_internal_error( "rdirichlet: Invalid alpha parameter" );
      x[i] = rgamma( alpha[i] );
      if( x[i] < MIN_DIRICHLET_PROB )
	x[i] = MIN_DIRICHLET_PROB;  /* BT 09/12/05 -- to prevent problems with 0 probs */
      sum += x[i];
    }
  
  for( i = 0; i < K; i++ )
    {
      theta[i] = x[i] / sum;
#ifdef _DEBUG_
      if( ! finite( theta[i] ) )
	p_internal_error( "rdirichlet: **** Non finite value ****" );		
      if( theta[i] <= 0)
	p_internal_error("rdirichlet: **** Theta == 0 ****");		
#endif
    }
  
  free( x );
}
  
  
/* return a model from a dirichlet distribution using counts + pseudocounts for each motif position */
/* take palindromes etc. into account */
void SampleAModel( Model B, int t, double **theta, double **motifCounts )
{
  int            n;
  int            j;
  int            mid;
  double         *counts;
  Ctype          C;
  IPtype         IP;
  double         *conc;
  int            c;
  
  C  = B->C;
  IP = B->IP;

  NEW( counts, IP->nAlphaLen, double );
  
  if( (! IP->is_defined[cl_R]) && (! IP->is_defined[cl_c]) && (! IP->is_defined[cl_a]) )
    {
      for(j = 0; j < IP->nMotifLen[t]; j++) 
	{
	  for(n = 0; n < IP->nAlphaLen; n++)
	    {
	      counts[n] = motifCounts[j + 1][n] + C->dPseudoCounts[t][j+1][n];
	    }
	  rdirichlet( B, IP->nAlphaLen, counts, theta[j] );
	}
    }
  else
    {
      mid = IP->nMotifLen[t] / 2 + IP->nMotifLen[t] % 2 - 1;
	   
      for( j = 0; j < IP->nMotifLen[t]; j++ )
	{
	  if( IP->AltModel->Palandromic[t][j] )
	    {
	      if( j <= mid )
		{
		  for( n = 0; n < IP->nAlphaLen; n++ )
		    {
		      counts[n] = (motifCounts[j + 1][n] + C->dPseudoCounts[t][j+1][n]) +
			(motifCounts[IP->nMotifLen[t]-j][nComp[n]] + 
			 C->dPseudoCounts[t][IP->nMotifLen[t]-j][nComp[n]]);
		    }
		  rdirichlet( B, IP->nAlphaLen, counts, theta[j] );
		}
	      else
		{
		  for( n = 0; n < IP->nAlphaLen; n++ )
		    theta[j][n] = theta[IP->nMotifLen[t] - j - 1][nComp[n]];
		}
	    }
	  else if( IP->AltModel->Repeat[t][j] )
	    {
	      if( j <= mid )
		{
		  for( n = 0; n < IP->nAlphaLen; n++ )
		    {
		      counts[n] = (motifCounts[j+1][n] + C->dPseudoCounts[t][j+1][n]) +
			(motifCounts[j + mid + 1][n] + C->dPseudoCounts[t][j+mid+1][n]);
		    }

		  rdirichlet( B, IP->nAlphaLen, counts, theta[j] );
		}
	      else
		{
		  for( n = 0; n < IP->nAlphaLen; n++ )
		    theta[j][n] = theta[j - mid - 1][n];
		}
	    }
	  else if( IP->AltModel->Collapsed[t][j] )  
	    {
	      for( n = 0; n < IP->nAlphaLen; n++ )
		{
		  counts[n] = (motifCounts[j+1][n] + C->dPseudoCounts[t][j+1][n]) +
		    (motifCounts[j+1][nComp[n]] +  C->dPseudoCounts[t][j+1][nComp[n]]);
		}
	      rdirichlet( B, IP->nAlphaLen, counts, theta[j] );
	    }
	  else
	    {
	      for( n = 0; n < IP->nAlphaLen; n++ )
		{
		  counts[n] = motifCounts[j+1][n] + C->dPseudoCounts[t][j+1][n];
		}
	      rdirichlet( B, IP->nAlphaLen, counts, theta[j] );
	    }
	  
	  if( IP->is_defined[cl_a] ) 
	    {
	      for(c = 1; c <= IP->AltModel->Concen->NumConcen[t]; c++) 
		{ 
		  NEW( conc, IP->nAlphaLen, double );
		  for( j = 0; j < IP->nMotifLen[t]; j++ )
		    {
		      if( IP->AltModel->Concen->Concentrated[t][j] == c ) 
			{			   
			  for( n = 0; n < IP->nAlphaLen; n++ )
			    {
			      conc[n] += counts[n];  
			    }
			}
		    }
		  for( j = 0; j < IP->nMotifLen[t]; j++ )
		    {
		      if( IP->AltModel->Concen->Concentrated[t][j] == c ) 
			{
			  for( n = 0; n < IP->nAlphaLen; n++ )
			    {
			      counts[n] = conc[n];
			    }
			}
		    }
		  free( conc );
		}
	      rdirichlet( B, IP->nAlphaLen, counts, theta[j] );	      
	    }
	}
     }
    
  free( counts );
}


void ImpGetBkgndCounts( Model B, double *bgCounts )
{
   IPtype     IP;
   BkgType    BP;
   PhyloType  PH;
   char       *R;
   int        seq;
   int        nLen;
   int        start;
   int        j;
   int        p;
   int        i;
   double     wt;
        
   IP = B->IP;
   BP = B->BP;
   R = (*(B)->Seq->R);
   PH = B->Phylo;
   
   for( i = 0; i < IP->nAlphaLen; i++)
     bgCounts[i] = 0.0;
   
   for( seq = 0; seq <= PH->maxPhyloSeq; seq += SeqInc( B, seq ) )	 
     {
       nLen = SequenceLength( B, seq );
       start = SequenceStartPos( B, seq );
       for( j = 0; j < nLen; j++ )
	 {
	   if( PossibleAlignedBkgndPos( B, start + j, nLen, seq ) )
	     {
	       for( p = 0; p < SeqInc( B, seq); p++ )
		 {
		   wt = GetSeqWeight( B, seq + p, 0 );
		   bgCounts[RCH2INT(j + start + p * nLen, R)] += wt;
		 }
	     }
	 }             
     }
}


void ImpGetBkgndCountsPos( Model B, int seq, int pos, double *bgCounts )
{
   IPtype     IP;
   BkgType    BP;
   PhyloType  PH;
   char       *R;
   int        nLen;
   int        start;
   int        p;
   int        i;
   int        nPos;
        
   IP = B->IP;
   BP = B->BP;
   R = (*(B)->Seq->R);
   PH = B->Phylo;
   
   for( i = 0; i < IP->nAlphaLen; i++)
     bgCounts[i] = 0.0;
   
   nLen = SequenceLength( B, seq );
   start = SequenceStartPos( B, seq );
   nPos = start + pos;

   if( PossibleAlignedBkgndPos( B, nPos, nLen, seq ) )
     {
       for( p = 0; p < SeqInc( B, seq ); p++ )
	 bgCounts[RCH2INT(nPos + p * nLen, R)] += GetSeqWeight( B, seq + p, 0 );
     }
}


/* returns the number of unique sites = maxPhyloSeq/phyoSpecies */
int GetMotifList( Model B, int t, MotifStruct **motifs )
{
  IPtype     IP;
  RPType     RP;
  Ftype      F;
  PhyloType  PH;
  int        seq;
  int        pos;
  int        sitePos;
  int        count = 0;
  int        nLen;
  int        j;
  int        col;
  int        pair;
  int        last = -1;
  int        rev;
  int        first = -1;

  IP = B->IP;
  RP = B->RP;
  F = B->F;
  PH = B->Phylo;

  if(! IP->is_defined[cl_F]) 
    {
      first = FirstCol( F, t );
      last = LastCol( F, t );
    }

  for( pair = 0, seq = 0; seq <= PH->maxPhyloSeq; seq += SeqInc( B, seq ), pair++ )
    {
      nLen = SequenceLength( B, seq );
      for( pos = 0; pos < nLen; pos++ )
	{
	  if( RP->sitePos[seq][pos][t].nMotifStart )
	    {
	      sitePos = SequenceStartPos( B, seq ) + pos;
	      rev = RP->sitePos[seq][pos][t].nRevComp;
	      for( j = 0; j < IP->nMotifLen[t]; j++ )
		{
		  if( ! IP->is_defined[cl_F] )
		    {
		      if( ! rev )
			col = F->fragPos[t][j];
		      else
			col = last - first - F->fragPos[t][j];
		    }
		  else
		    {
		      if( ! rev )
			col = j;
		      else
			col = IP->nMotifLen[t] - j - 1;
		    }
		  motifs[j][count].pos = sitePos + col;
		  motifs[j][count].rev = rev;
		  motifs[j][count].seq = seq;
		  motifs[j][count].seqPos = pos + col;
		  motifs[j][count].seqPair = pair;
		  motifs[j][count].seqLen = nLen;
		}
	      count++;
	    }
	}
    }
  return count;
}


int PossibleAlignedBkgndPos( Model B, int pos, int offset, int seq )
{
  int  n;
  int  seqPos;
  char ch;

  if( B->RP->nInMotif[pos] )
    return FALSE;

  if( B->Phylo->phyloTree ) 
    {    
      seqPos = PosInFirstMASSSeq( B, pos );
      for( n = 0; n < SpeciesInc( B, B->Phylo->phyloTreeStartSeq[B->Phylo->phyloIndex[seq]] ); n++ )
	{
	  ch = (*B->Seq->R)[seqPos + n * offset];
	  if( ch == 'n' || ch == 'x' )
	    return FALSE;
	}
    }
  
  return TRUE;
}


/* get the motif counts in the phylogentic sequences */
/* this routine assumes that counts is zeroed */
/* counts is 1 longer than the motif length to make it compatable with C->wCounts */
void GetPhyloCounts( Model B, int t, int num, MotifStruct **motifs, double **counts )
{
  PhyloType  PH;
  IPtype     IP;
  int        index;
  int        i;
  int        j;
  int        k;
  int        rev;
  int        pos;
  int        offset;
  int        seq;

  IP = B->IP;
  PH = B->Phylo;

  for( i = 0; i < num; i++ )
    {
      for( j = 0; j < IP->nMotifLen[t]; j++ )
	{
	  rev = motifs[j][i].rev;
	  pos = motifs[j][i].pos;
	  offset = motifs[j][i].seqLen;
	  seq = motifs[j][i].seq;
	  for( k = 0; k < SeqInc( B, seq ); k++ )
	    {
	      index = (! rev) ?	
		RCH2INT( pos + k * offset, *(B)->Seq->R ) :
		nComp[RCH2INT(pos + k * offset, *(B)->Seq->R )]; 

	      counts[j+1][index] += GetSeqWeight( B, seq + k, t );
	    }
	}
    }    
}

