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
#include "freqbkgnd.h"

#define LINE_SIZE   2048
#define TOKEN_SIZE  1024

int ReadFreqLine( FILE *fpt, char *line, int lineSize, char *token, int tokenSize,
		  char *nMer, int *count, double *prob );
int NmerIndex( Model B, int seq, int *nMer );


void ReadFreqBackground( Model B, char *fileName, double **nMerArray, int *size )
{
  IPtype   IP;
  FILE     *fpt;
  char     *line;
  char     *token;
  char     *nMer;
  int      count=0;
  double   prob;
  int      arraySize;
  int      index;

  IP = B->IP;

  NEW( line, LINE_SIZE, char );
  NEW( token, TOKEN_SIZE, char );
  NEW( nMer, TOKEN_SIZE, char );

  fpt = fopen( fileName, "r" );
  if( fpt == NULL )
    {
      sprintf( line, "ReadFreqBackground: Can't open frequency file: %s\n", fileName );
      p_error( line );
    }

  ReadFreqLine( fpt, line, LINE_SIZE, token, TOKEN_SIZE,
		nMer, &count, &prob );
  *size = strlen( nMer );
  if( count < 1 || *size > MAX_NMER || *size < 1 )
    {
      sprintf( line, "ReadFreqBackground: error getting n-mer size: %s\n", nMer );
      p_error( line );
    }

  if( *size > IP->nMotifLen[MinLenMotif( B )] )
    {
      sprintf( line, "ReadFreqBackground: n-mer size must be <= %d\n", MinLenMotif( B ) );
      p_error( line );
    }

  arraySize = pow( 4, *size );
  NEW( *nMerArray, arraySize, double );

  fseek( fpt, 0, SEEK_SET );
  while( ReadFreqLine( fpt, line, LINE_SIZE, token, TOKEN_SIZE,
		       nMer, &count, &prob ) )
    {
      index = NmerToPosition( nMer, *size );
      if( index >= arraySize || index < 0 || prob < 0 || prob > 1 )
	{
	  sprintf( line, "ReadFreqBackground: invalid n-mer: %s\n", nMer );
	  p_error( line );	  
	}
      (*nMerArray)[index] = prob;
    }
  
  
  free( line );
  free( token );
  free( nMer );
  fclose( fpt );
}


int ReadFreqLine( FILE *fpt, char *line, int lineSize, char *token, int tokenSize,
		  char *nMer, int *count, double *prob ) 
{
  char *ptr;
  char *ptr2;

  if( ! GetALine( fpt, line, lineSize, FALSE ) )
    return FALSE;

  ptr = line;
  ptr = GetAToken( ptr, token );
  strncpy( nMer, token, tokenSize );
  ptr = GetAToken( ptr, token );
  *count = atoi( token );
  ptr = GetAToken( ptr, token );
  *prob = strtod( token, &ptr2 );

  return TRUE;
}


int NmerToPosition( char *nMer, int nMerSize )
{
  int   i;
  int   pos = 0;
  int   n = 0;
  char  *msg;

  for( i = 0; i < nMerSize; i++ )
    {
      switch( toupper( nMer[i] ) )
	{
	case 'A':
	  n = 0;
	  break;

	case 'T':
	  n = 1;
	  break;

	case 'C':
	  n = 2;
	  break;

	case 'G':
	  n = 3;
	  break;
	  
	default:
	  NEW( msg, LINE_SIZE, char );
	  sprintf( msg, "NmerToPosition: invalid nucleotide %c in %s", nMer[i], nMer );
	  p_error( msg );
	  break;
	}

      pos += n * (1 << 2*i);
    }

  return pos;
}


void FreeNMerArrays( Model B )
{
  IPtype   IP;
  int      i;
  int      j;

  IP = B->IP;

  for( i = 0; i < IP->nNumSequences; i++ )
    {
      if( IP->freqSeqs[i] )
	{
	  for( j = i + 1; j < IP->nNumSequences; j++ )
	    {
	      if( IP->freqSeqs[i] == IP->freqSeqs[j] )
		IP->freqSeqs[j] = NULL;
	    }
	  free( IP->freqSeqs[i] );
	}
    }
}


int IsNmerBkgnd( Model B, int seq )
{
  return (B->IP->is_defined[cl_freq_background] && 
	  (B->IP->freqSeqs[seq] != NULL) );
}


double NmerProb( Model B, int seq, int t, int *nMer )
{
  IPtype   IP;
  int      index;
  int      len;
  int      size;
  int      *nMer2;
  int      i;
  int      j;
  int      k;
  double   prob;
  double   sum;

  IP = B->IP;

  if( IP->nMerSize[seq] == IP->nMotifLen[t] )
    {
      index = NmerIndex( B, seq, nMer );
      return IP->freqSeqs[seq][index];
    }
  else
    {
      len = IP->nMotifLen[t];
      size = IP->nMerSize[seq];
      NEW( nMer2, size, int );
      for( i = 0; i < size; i++ )
	nMer2[i] = nMer[i];
      index = NmerIndex(B, seq, nMer2 );
      prob = IP->freqSeqs[seq][index];

      for( i = size; i < len; i++ )
	{
	  for( k = 0, j = i - size + 1; j <= i; j++, k++ )
	    nMer2[k] = nMer[j];
	  index = NmerIndex(B, seq, nMer2 );
	  prob *= IP->freqSeqs[seq][index];

	  for( sum = 0, j = 0; j < 4; j++ )
	    {
	      nMer2[size-1] = j;
	      index = NmerIndex(B, seq, nMer2 );
	      sum += IP->freqSeqs[seq][index];
	    }
	  prob /= sum;
	}

      free( nMer2 );

      return prob;
    }
}


int NmerIndex( Model B, int seq, int *nMer )
{
  IPtype   IP;
  int      nMerSize;
  int      index = 0;
  int      i;

  IP = B->IP;
  nMerSize = IP->nMerSize[seq];

  for( i = 0; i < nMerSize; i++ )
    {
      index += nMer[i] * (1 << 2*i);      
    }
  
  return index;
}



