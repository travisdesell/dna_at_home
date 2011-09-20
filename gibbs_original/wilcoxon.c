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
/* $Id: wilcoxon.c,v 1.9 2007/05/23 18:19:58 Bill Exp $                   */
/*                                                                        */
/* Author :       Bill Thompson 6/30/97                                   */
/*                                                                        */
/* Description :  This file contains the functions which are used to      */
/*                perform the Wilcoxon signed rank test.                  */
/*                                                                        */
/**************************************************************************/

#include "wilcoxon.h"

void ShuffleSequence( IPtype IP, int* Seq, int start, int end );
short CompareSequences( Model B );
void WriteWilcoxonSequences( Model B, FILE *fpt );

void Wilcoxon( Model B )
{
  IPtype IP;
  int    i;
  FILE   *fpt;
  int    nSeqCount;

  IP= B->IP;
  sRandom(B, IP->lSeedVal);  

  get_strings(B);           /* Read the sequences */

  IP->Datafiles->WilcoxonFilename = GibbsTempnam( NULL, "wil" ); 

  fpt = fopen( IP->Datafiles->WilcoxonFilename, "w" );
  if( fpt == NULL )
    p_error( "Unable to open Wilcoxon output file." );

  chmod( IP->Datafiles->WilcoxonFilename, S_IRUSR | S_IWUSR );

  if( IP->is_defined[cl_D] || B->Phylo->bCalcPhylo )
    p_error( "Wilcoxon: With a tree or -D, include sequences in fasta file and use -wilcox" );
  else
    WriteWilcoxonSequences( B, fpt );
  
  fclose( fpt );

  IP->nOrigSeqCount = IP->nNumSequences;
  nSeqCount = IP->nNumSequences;

  fclose( B->IP->Datafiles->fpt );	    
  free(B->IP->Datafiles->filename);

  IP->Datafiles->filename = (char *)malloc(strlen(IP->Datafiles->WilcoxonFilename ) * 
					   sizeof(char) + 1);
  strcpy(IP->Datafiles->filename, IP->Datafiles->WilcoxonFilename);
  IP->Datafiles->fpt = fopen(IP->Datafiles->filename, "r"); 
             
  free(*B->Seq->R);
  free(*B->Seq->ProcessedSTR);
  free(B->Seq->ProcessedSTR);
  free(*B->Seq->nvEndLocs);    /* BT 8/27/97 */
  free(B->Seq->nvEndLocs);    /* BT 8/27/97 */
  free(B->Seq->R);
  for(i=0; i < B->IP->nSeqAllocCnt; i++ )    /* BT 8/27/97 */
    free(B->Seq->Orig[i]);
  free(B->Seq->Orig);
  free(B->Seq->SeqLen);
  free(B->Seq);
  
  NEW(B->Seq, 1, Stringstruct);
  NEWP(B->Seq->R, 1, char);
  NEWP(B->Seq->nvEndLocs, 1, int);
  NEWP(B->Seq->ProcessedSTR, 1, char);

  get_strings( B );  

  if( ! CompareSequences( B ) )
    p_error( "Wilcoxon control file error." );
}


void ShuffleSequence( IPtype IP, int* Seq, int start, int end )
{
  char  alpha[20] = "actgdefhiklmnpqrsvwy";
  int   j;
  int   in;
  int   out;
  short tmp;

  for( j = start; j < end; j++ )
    {
      in = j;
      if( (Seq[in] == 'x' || Seq[in] == 'X' || Seq[in] == 'u') ||
	  ((IP->nAlphaLen == 4) &&
	   (Seq[in] == 'n' || Seq[in] == 'N'))) 
	Seq[in] = alpha[RandomInterval( 0, IP->nAlphaLen - 1)];
      else
	{
	  out = RandomInterval( in, end );
	  if( (Seq[out] == 'x' || Seq[out] == 'X' || Seq[out] == 'u') ||
	      ((IP->nAlphaLen == 4) &&
	       (Seq[out] == 'n' || Seq[out] == 'N'))) 
	    Seq[out] = alpha[RandomInterval( 0, IP->nAlphaLen - 1)];
	  if( in != out )
	    {
	      tmp = Seq[in];
	      Seq[in] = Seq[out];
	      Seq[out] = tmp;
	    }
	}
    }
}


void ShuffleSequence2( IPtype IP, short* Seq, int start, int end )
{
  char  alpha[20] = "actgdefhiklmnpqrsvwy";
  int   i;
  int   j;
  int   in;
  int   out;
  short tmp;

  for(i = 0; i < 7; i++ )
    {
      for( j = start; j <= end; j++ )
	{
	  in = j;
	  if( (Seq[in] == 'x' || Seq[in] == 'X' || Seq[in] == 'u') ||
	      ((IP->nAlphaLen == 4) &&
	       (Seq[in] == 'n' || Seq[in] == 'N'))) 
	    Seq[in] = alpha[RandomInterval( 0, IP->nAlphaLen - 1)];
	  else
	    {
	      out = RandomInterval( start, end );
	      if( (Seq[out] == 'x' || Seq[out] == 'X' || Seq[out] == 'u') ||
		  ((IP->nAlphaLen == 4) &&
		   (Seq[out] == 'n' || Seq[out] == 'N'))) 
		Seq[out] = alpha[RandomInterval( 0, IP->nAlphaLen - 1)];
	      if( in != out )
		{
		  tmp = Seq[in];
		  Seq[in] = Seq[out];
		  Seq[out] = tmp;
		}
	    }
	}
    }  
}


void WriteWilcoxon( Model B, MaxResults maxData, int t )
{
  IPtype IP;
  char   *tempname;
  FILE*  fpt;
  int    i;
  int    seq;
  int    loc;
  char   *outfile;
  char   *comStr;
  int    ch;
  char   *msg;
#ifdef _CYGWIN_
  char   *comStr2;
#endif
  
  IP = B->IP;

  NEW( comStr, 3 * FILENAME_MAX + 32, char ); 

  tempname = GibbsTempnam( NULL, "wil" ); 
  outfile = GibbsTempnam( NULL, "wil" );

  fpt = fopen( tempname, "w" );
  if( fpt == NULL )
    p_error( "Unable to open Wilcoxon output file." );
      
  for( i = 0; i < maxData.nNumMotifs[t]; i++ ) 
    {
      loc = maxData.nMotifLoc[i][t];
      seq = SequenceFromPosition( B, loc );
      
      if( strstr( IP->fastaHeader[seq], "Wilcoxon" ) != IP->fastaHeader[seq] )
	fprintf( fpt, "0.0 %.4f\n", maxData.dvMotifProb[i][t] );
      else
	fprintf( fpt, "%.4f 0.0\n", maxData.dvMotifProb[i][t] );
    }
  fclose( fpt );
  
  remove( outfile );

  for( i = strlen( IP->argv[0] ) - 1; i >= 0; i-- )
    {
      if( IP->argv[0][i] == '/' )
	break;
    }
  strncpy( comStr,IP->argv[0] , i + 1 );
  if( strstr( IP->argv[0], ".x86" ) == NULL )
    strcat( comStr, "wilcox.test" );
  else
    strcat( comStr, "wilcox.test.x86" );

  if( ! FileExists( comStr ) )
    {
      NEW( msg, 3 * FILENAME_MAX + 32, char );
      sprintf( msg, "WriteWilcoxon: can't find %s.", comStr );
      p_warning( msg );
      free( msg );

      remove( outfile );
      remove( tempname );  
      free( comStr );
      free( outfile );
      free( tempname );    

      return;
    }

#ifdef _CYGWIN_
  /* deal with spaces in path name for Cygwin */
  NEW( comStr2, strlen( comStr ) + 32, char );
  strcpy( comStr2, "'" );
  strcat( comStr2, comStr );
  strcat( comStr2, "'" );
  strcpy( comStr, comStr2 );
  free( comStr2 );
#endif
    
  strcat( comStr, " " );
  strcat( comStr, tempname );
  strcat( comStr, " -P -L >" );
  strcat( comStr, outfile );
  
  system( comStr );

  fpt = fopen( outfile, "r" );
  if( fpt == NULL )
    p_error( "Unable to open results file for Wilcoxon test." );
  
  fprintf( IP->Datafiles->out_fpt, 
	   "\n\n======================= Wilcoxon Signed-rank Test for Motif %c ========================\n",
	   (char)('a' + t));
  
  while( (ch = getc( fpt )) != EOF  )
    putc( ch, IP->Datafiles->out_fpt );
  fclose( fpt );

  remove( outfile );
  remove( tempname );
  
  free( comStr );
  free( outfile );
  free( tempname );
}


void WriteWilcoxonSPlus( Model B, MaxResults maxData, int t )
{
  IPtype IP;
  char   *tempname;
  FILE*  fpt;
  int    i;
  int    j;
  int    seq;
  int    loc;
  char   *splusfile;
  char   *splusOutfile;
  char   *splusVar;
  char   *comStr;
  int    ch;
  
  IP = B->IP;

  NEW( comStr, 256, char ); 

  tempname = GibbsTempnam( NULL, "wil" ); 
  splusfile = GibbsTempnam( NULL, "wil" ); 
  splusOutfile = GibbsTempnam( NULL, "wil" );

  splusVar = GibbsTempnam( NULL, "wil" );
  for( i = 0; i < strlen( splusVar ); i++ )    /* BT 8/29/97 */
    {
      if( splusVar[i] == '_' )
	splusVar[i] = '.';
    }

  fpt = fopen( tempname, "w" );
  if( fpt == NULL )
    p_error( "Unable to open Wilcoxon output file." );
      
  for( i = 0; i < maxData.nNumMotifs[t]; i++ ) 
    {
      loc = maxData.nMotifLoc[i][t];
      j = 0;
      while( loc > (*B->Seq->nvEndLocs)[j])
	j++;
      
      seq = j;
      
      if( seq < IP->nOrigSeqCount )
	fprintf( fpt, "0.0 %.4f\n", maxData.dvMotifProb[i][t] );
      else
	fprintf( fpt, "%.4f 0.0\n", maxData.dvMotifProb[i][t] );
    }
  fclose( fpt );
  
  remove( splusOutfile );
  
  fpt = fopen( splusfile, "w" );
  if( fpt == NULL )
    p_error( "Unable to open SPlus script file." );
  
  sprintf( comStr, "%s_matrix( scan(\"%s\"), 2, %d )\n", 
	   &splusVar[5], tempname, maxData.nNumMotifs[t] );
  
  /*      fprintf( fpt, "xxx_matrix( scan(\"%s\"), 2, %d )\n", tempname, maxData.nNumMotifs[t] ); */
  fprintf( fpt, "%s", comStr );
  fprintf( fpt, "sink(\"%s\")\n", splusOutfile );
  fprintf( fpt, "wilcox.test(%s[1,],%s[2,], alt = \"less\",paired=T)\n",
	   &splusVar[5], &splusVar[5] );
  fprintf( fpt, "sink()\n" );
  fprintf( fpt, "remove(\"%s\")\n",  &splusVar[5] ); 
  fclose( fpt );
  
  strcpy( comStr, "Splus <" );
  strcat( comStr, splusfile );
  
  system( comStr );
  fpt = fopen( splusOutfile, "r" );
  if( fpt == NULL )
    p_error( "Unable to open SPlus results file." );
  
  fprintf( IP->Datafiles->out_fpt, 
	   "\n\n======================= Wilcoxon Signed-rank Test for Motif %c ========================\n",
	   (char)('a' + t));
  
  while( (ch = getc( fpt )) != EOF  )
    putc( ch, IP->Datafiles->out_fpt );
  fclose( fpt );

  remove( splusOutfile );
  remove( splusfile );
  remove( tempname );
  
  free( comStr );
  free( splusVar );
  free( splusOutfile );
  free( splusfile );
  free( tempname );
}


short CompareSequences( Model B )
{
  IPtype IP;
  char   msg[256];
  int    i;

  IP = B->IP;

  if(  IP->nNumSequences % 2 != 0 )
    {
      p_error( "The Wilcoxon control file must contain the same number of sequences as the data file." );
      return FALSE;
    }

  for( i = 0; i < IP->nNumSequences / 2; i++ )
    {
      if( SequenceLength( B, i ) != SequenceLength( B, i + (IP->nNumSequences / 2) ) )
	{
	  sprintf( msg, "Wilcoxon Control file: Sequence %d has a different length than the sequence in the data file.", i + 1 );
	  p_error( msg );
	  return FALSE;
	}
    }

  return TRUE;      
}


void WriteWilcoxonSequences( Model B, FILE *fpt )
{
  int    i;
  int    j;
  int    k;
  IPtype IP;
  int    start;
  int    end;
  int*   newSeq = NULL;
  char   tmpch;
  
  IP = B->IP;

  if( IP->Datafiles->control_fpt == NULL )
    {
      NEW( newSeq, IP->nSeqLen, int );
      for( i = 0; i < IP->nSeqLen; i++ )
	newSeq[i] = (*B->Seq->R)[i];

      for( i = 0; i < IP->nNumSequences; i++ )
	{
	  ShuffleSequence( IP, newSeq, SequenceStartPos( B, i ), SequenceEndPos( B, i ) );
	}
    }

  for( i = 0; i < IP->nNumSequences; i++ )
    {
      start = SequenceStartPos( B, i );
      end = SequenceEndPos( B, i );
      fprintf( fpt, ">%s", B->IP->fastaHeader[i] );
      for( k = 0, j = start; j <= end; j++, k++ )
	{
	  if( k % 60 == 0 )
	    fprintf( fpt, "\n" );
	  fprintf( fpt, "%c", curr_ch( (*B->Seq->R)[j], IP->nAlphaLen, TRUE ) );
	}
      fprintf( fpt, "\n\n" );
    }

  if( IP->Datafiles->control_fpt != NULL )
    {
      fprintf( fpt, "\n" );
      while(fscanf(B->IP->Datafiles->control_fpt, "%c", &tmpch) != EOF)
	{
	  fprintf( fpt, "%c", tmpch );
	  if( tmpch == '>' )
	    {
	      fprintf( fpt, "Wilcoxon " );
	    }
	}
      fprintf( fpt, "\n\n" );
      fclose( B->IP->Datafiles->control_fpt );
    }
  else
    {
      for( i = 0; i < IP->nNumSequences; i++ )
	{
	  start = SequenceStartPos( B, i );
	  end = SequenceEndPos( B, i );
	  fprintf( fpt, ">Wilcoxon Shuffled Sequence %d", i + 1 );
	  for( k = 0, j = start; j <= end; j++, k++ )
	    {
	      if( k % 60 == 0 )
		fprintf( fpt, "\n" );
	      fprintf( fpt, "%c", curr_ch( newSeq[j], IP->nAlphaLen, TRUE ) );
	    }
	  fprintf( fpt, "\n\n" );
	}
    }

  if( IP->Datafiles->control_fpt == NULL )
    free( newSeq );
 }
