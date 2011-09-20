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
#include "bayes.h"


void PrintSampleResults( Model B, int ***counts, int seed_run );
#ifdef _MPI_
void SendSampleResults( Model B, int ***posCounts, int seed_run);
#endif
void CopyCounts( Model B, int ***counts, int ***occurence, int seed_run, int iter );
void CopyTotalCounts( Model B, int ***counts, int **kCounts );


MaxResults bayes_sampler(Model B, PoSition **Pos, Mlist M, int seed_run)
{
   int             iter, last_increase;
   int             n, t;
   int             i;
   int             k;
   double          dLocMax = -DBL_MAX, dCurrProb, dMaxProbability = 0.0;
   double          **dProb, **dFinalProb;
   ProbStruct      P;
   MaxResults      maxData;
   IPtype          IP;
   RPType          RP;
   int             ***occurence;
   double          **dProbArray;
   double          dMap;
   short           bCount;
   int             nFragWidth;
   struct tms      timeBuffer;
   clock_t         currTime;
   struct tm       *localTime;
   struct tm       tm_buf;
   time_t          tloc;
   int             ***posCounts;
   int             nWidths;
#ifdef _MPI_
   double         prevTemp;
#endif

   IP = B->IP;
   RP = B->RP;

   if( ! IP->is_defined[cl_E] )
     p_error( "bayes: -E must be users with -bayes." );

   RP->nIncludeAlign = TRUE;
   RP->nScaleZeroBlocksProb = TRUE;
   RP->nDontAdjustPrior = TRUE;
   RP->nUseFixedBlockProbs = TRUE;

   if( IP->is_defined[cl_X] )
     {
       if( ! B->AN->bExchange )
	 B->AN->currTemp = B->AN->dMaxTemp;
       currTime = (times( &timeBuffer ) / IP->ticks) * 1000;
       B->AN->lastExchangeTime = currTime;
       time( &tloc );
       localTime = localtime_r( &tloc, &tm_buf );
       if( ! B->IP->is_defined[cl_Z] )
	 {
	   printf( "Current temperature: %7.2f at %d:%d:%d\n", B->AN->currTemp, 
		   localTime->tm_hour, localTime->tm_min, localTime->tm_sec );
	   PrintTempOut(  IP->Datafiles->mpiTemp_fpt, "currTime = %d lastExchange = %d exchangePeriod = %d\n", 
			  currTime,  B->AN->lastExchangeTime, B->AN->exchangePeriod );   /* DEBUG */
	 }
     }

   init_maxdata(&maxData);
   NEWP(dFinalProb, B->IP->nSeqLen, double);
   for(n = 0; n < B->IP->nSeqLen; n++)
     NEW(dFinalProb[n], 2, double); 
   
   NEWP(dProb, B->IP->nSeqLen, double);
   for(n = 0; n < B->IP->nSeqLen; n++)
     NEW(dProb[n], 2, double);
   
   NEWPP(occurence, IP->nSeqLen, int);
   for(n = 0; n < IP->nSeqLen; n++) 
     {
       NEWP(occurence[n], IP->nNumMotifTypes, int);
       for(t = 0; t < IP->nNumMotifTypes; t++) 
	 NEW(occurence[n][t], 3, int);
     }

   NEWP( RP->kCounts, IP->nNumSequences, int );
   for( n = 0; n < IP->nNumSequences; n++ )
     NEW( RP->kCounts[n], RP->nMaxBlocks + 1, int );
   
   reset_values(&iter, &last_increase, &dLocMax, &P, B->IP); 
   update_prob(&P, B, (! B->IP->is_defined[cl_b]));

   if( B->RP->bUsePosMatrix )
     InitializePosMatrix( B );
   else if( B->RP->bUseTrans )
     InitializeTransMatrix( B );

   if( B->InitPos != NULL && seed_run == 1 && B->IP->is_defined[cl_opt] )
     {
       if( ! B->IP->is_defined[cl_Z] )
	 {
	   time( &tloc );
	   localTime = localtime_r( &tloc, &tm_buf );
	   
	   fprintf( stdout, "Max set %f at %d motifs: %d %d:%d:%d\n", 
		  B->RP->dInitProb, iter, TotalNumMotifs( B ), 
		  localTime->tm_hour, localTime->tm_min, localTime->tm_sec);
	 }
       dProbArray = setElementProb(B, Pos, P);
       free_maxdata(&maxData, B->IP);
       maxData = setMaxData(B->F, B->IP, B->RP->dInitProb, 0, Pos, &last_increase,
			    &dMaxProbability, &dLocMax, maxData, dProbArray, B);
       FREEP(dProbArray, findMaxNumMotif(B->IP));
       } 

   if( ! B->IP->is_defined[cl_Z] )
     {
       fprintf(stdout, "Burn-in\n"); 
       fflush( stdout );    /* BT 9/19/97 */
     }

   IP->bayesSampleIter = -1;
   for( iter = 0; iter < IP->burnInPeriod; iter++ )
     {
       if( ! IP->is_defined[cl_Z] )
	 {
	   fprintf(stdout, "\r%d", iter); 
	   fflush( stdout );    /* BT 9/19/97 */
	 }
       dCurrProb = rsample( B, Pos, M, &occurence, &P, iter, 
			    seed_run, FALSE, &dMap, dLocMax );	
       
       if( dCurrProb > dLocMax && iter > WAIT_PERIOD && B->IP->is_defined[cl_opt] )
	 {
	   if( (! B->IP->is_defined[cl_Z]) )
	     {
	       time( &tloc );
	       localTime = localtime_r( &tloc, &tm_buf );
	       
	       fprintf( stdout, "Max set %f at %d motifs: %d %d:%d:%d\n", 
			dCurrProb, iter, TotalNumMotifs( B ),
			localTime->tm_hour, localTime->tm_min, localTime->tm_sec);
	     }
	   dProbArray = setElementProb(B, Pos, P);
	   free_maxdata(&maxData, B->IP);
	   maxData = setMaxData(B->F, B->IP, dCurrProb, iter, Pos, &last_increase,
				&dMaxProbability, &dLocMax, maxData, dProbArray, B);
	   FREEP(dProbArray, findMaxNumMotif(B->IP));
	 } 

       if(iter % IP->nAdjustPeriod == 0) 
	 {
	   if(!B->IP->is_defined[cl_F]) 
	     {
	       bCount = FALSE;
	       for(t = 0; t < B->IP->nNumMotifTypes; t++)
		 {
		   StoreMask(B); 
		   setFragVec(B, M, Pos);   /* BT 4/26/04 */
		   setFragCnts(B, M);
		   nFragWidth = B->F->FragWidth[t];
		   Fragment(B, t, M);
		   adjust_fragment(B, M, Pos, &P);
		   if( nFragWidth != B->F->FragWidth[t] )
		     bCount = TRUE;
		 }
	       /* adjust_fragment(B, M, Pos, &P); */
	       if( bCount )
		 CountAlignments( B, Pos ); 
	     }
	   else
	     {
	       if( B->IP->is_defined[cl_d] )           
		 ResizeMotif( B, Pos, M, FALSE );            
	       ColumnShift(B, M, Pos, &P);        
	     }

	   if( RP->bUsePosMatrix && RP->bUpdatePosMatrix )
	     UpdatePosMatrix( B );
	   else if( RP->bUseTrans && RP->bUpdateTrans )
	     UpdateTransMatrix( B );
	 }

       if(iter > WAIT_PERIOD && B->IP->is_defined[cl_opt] ) 
	 {
	   check_map(B, &dCurrProb, iter, Pos, &last_increase,
		     &dMaxProbability, &dLocMax, &maxData, P); 	   
	 }
     }
   
   for( n = 0; n < IP->nNumSequences; n++ )
     {
       for( k = 0; k <= RP->nMaxBlocks; k++ )
	 RP->kCounts[n][k] = 0;
     }

   nWidths = NumWidths( B );
   NEWPP( posCounts, IP->nSeqLen, int );
   for(n = 0; n < IP->nSeqLen; n++) 
     {
       NEWP( posCounts[n], IP->nNumMotifTypes, int);
       for(t = 0; t < IP->nNumMotifTypes; t++) 
	 NEW( posCounts[n][t], nWidths, int );
     }
   
   if( ! IP->is_defined[cl_Z] )
     {
       fprintf(stdout, "Sampling\n"); 
       fflush( stdout );    /* BT 9/19/97 */
     }

   for( iter = 0; iter < IP->bayesSamplePeriod; iter++ )
     {
       IP->bayesSampleIter = iter;

       if( ! B->IP->is_defined[cl_Z] )
	 {
	   fprintf(stdout, "\r%d", iter); 
	   fflush( stdout );    /* BT 9/19/97 */
	 }

       for(n = 0; n < IP->nSeqLen; n++) 
	 {
	   for(t = 0; t < IP->nNumMotifTypes; t++) 
	     {
	       for( i = 0; i < 3; i++ )
		 occurence[n][t][i] = 0;
	     }
	 }

       dCurrProb = rsample( B, Pos, M, &occurence, &P, iter, 
			    seed_run, IP->is_defined[cl_Q], &dMap, dLocMax );	

       CopyCounts( B, posCounts, occurence, seed_run, iter );
       
       if( dCurrProb > dLocMax && iter > WAIT_PERIOD && B->IP->is_defined[cl_opt] )
	 {
	   if( (! B->IP->is_defined[cl_Z]) )
	     {
	       time( &tloc );
	       localTime = localtime_r( &tloc, &tm_buf );
	       
	       fprintf( stdout, "Max set %f at %d motifs: %d %d:%d:%d\n", 
			dCurrProb, iter, TotalNumMotifs( B ),
			localTime->tm_hour, localTime->tm_min, localTime->tm_sec);
	     }
	   dProbArray = setElementProb(B, Pos, P);
	   free_maxdata(&maxData, B->IP);
	   maxData = setMaxData(B->F, B->IP, dCurrProb, iter, Pos, &last_increase,
				&dMaxProbability, &dLocMax, maxData, dProbArray, B);
	   FREEP(dProbArray, findMaxNumMotif(B->IP));
	 } 

       if(iter % IP->nAdjustPeriod == 0) 
	 {
	   if(!B->IP->is_defined[cl_F]) 
	     {
	       bCount = FALSE;
	       for(t = 0; t < B->IP->nNumMotifTypes; t++)
		 {
		   StoreMask(B); 
		   setFragVec(B, M, Pos);  
		   setFragCnts(B, M);
		   nFragWidth = B->F->FragWidth[t];
		   Fragment(B, t, M);
		   adjust_fragment(B, M, Pos, &P);
		   if( nFragWidth != B->F->FragWidth[t] )
		     bCount = TRUE;
		 }
	       if( bCount )
		 CountAlignments( B, Pos ); 
	     }
	   else
	     {
	       if( B->IP->is_defined[cl_d] )           
		 ResizeMotif( B, Pos, M, FALSE );            
	       ColumnShift(B, M, Pos, &P);        
	     }

	   if( RP->bUsePosMatrix && RP->bUpdatePosMatrix )
	     UpdatePosMatrix( B );
	   else if( RP->bUseTrans && RP->bUpdateTrans )
	     UpdateTransMatrix( B );
	 }

       if(iter > WAIT_PERIOD && B->IP->is_defined[cl_opt] ) 
	 {
	   check_map(B, &dCurrProb, iter, Pos, &last_increase,
		     &dMaxProbability, &dLocMax, &maxData, P); 
	   
#ifdef _MPI_
	   if( IP->is_defined[cl_X] && 
	       B->AN->bExchange )
	     {
	       currTime = (times( &timeBuffer ) / IP->ticks) * 1000;
	       if( ! B->IP->is_defined[cl_Z] )
		 PrintTempOut(  IP->Datafiles->mpiTemp_fpt, 
				"currTime = %d lastExchange = %d exchangePeriod = %d diff = %d\n", 
				currTime,  B->AN->lastExchangeTime, B->AN->exchangePeriod,
				currTime - B->AN->lastExchangeTime );   /* DEBUG */
	       if( (B->AN->exchangePeriod == 0 && (iter % B->AN->nExchangeIterations == 0)) ||
		   (B->AN->exchangePeriod > 0 && (currTime - B->AN->lastExchangeTime > B->AN->exchangePeriod)) )
		 {
		   prevTemp = B->AN->currTemp;
		   ExchangeTemps( B, dCurrProb, seed_run, iter );
		   currTime = (times( &timeBuffer ) / IP->ticks) * 1000;
		   B->AN->lastExchangeTime = currTime;
		 }   
	     }

	   if( IP->is_defined[cl_hm] && iter % IP->hier_iter == 0 )
	     ExchangeK( B, seed_run, iter );
#endif	 
	 }
     }

   free_prob(&P, B->IP); 
   maxData.frequency = NULL;
   FREEP(dFinalProb, B->IP->nSeqLen);
   FREEP(dProb, B->IP->nSeqLen);
   if(!B->IP->is_defined[cl_F]) 
     {  
       if(  B->F->nvFragCnts != NULL )
	 {
	   for( t = 0; t < B->IP->nNumMotifTypes; t++ ) 	/* BT 4/2/97 */
	     FREEP( B->F->nvFragCnts[t], B->F->nMaxLen[t] );
	   free( B->F->nvFragCnts );
	   B->F->nvFragCnts = NULL;
	 }
     }   
   
   free(dFinalProb);
   free(dProb);
   
#ifdef _MPI_
       SendSampleResults( B, posCounts, seed_run );
#else
   if( IP->is_defined[cl_bayes_counts] )
       PrintSampleResults( B, posCounts, seed_run );
#endif
      
   CopyTotalCounts( B, posCounts, RP->kCounts );
		
   FREEPP(occurence, IP->nSeqLen, IP->nNumMotifTypes);
   FREEPP(posCounts, IP->nSeqLen, IP->nNumMotifTypes);
   FREEP( RP->kCounts, IP->nNumSequences );

   return(maxData);
}


#ifdef _MPI_
void SendSampleResults( Model B, int ***posCounts, int seed_run)
{
  IPtype      IP;
  RPType      RP;
  double      dTempInfo[4];
  int         *counts;
  MPI_Status  status;
  int         len;
  int         w;
  int         j;
  int         t;
  int         n;
  int         seq;
  int         nWidths;
  int         k;
  
  IP = B->IP;
  RP = B->RP;

  dTempInfo[0] = seed_run;
  Gibbs_MPI_Send( B, dTempInfo, 4, MPI_DOUBLE, 0, G_MPI_BAYES, MPI_COMM_WORLD );
  Gibbs_MPI_Recv( B, dTempInfo, 4, MPI_DOUBLE, 0, G_MPI_BAYES, MPI_COMM_WORLD, &status);
  len = (int) dTempInfo[0];

  NEW( counts, len, int );

  nWidths = NumWidths( B );

  for( j = 0, n = 0; n < IP->nSeqLen; n++) 
    {
      for(t = 0; t < IP->nNumMotifTypes; t++) 
	{
	  for( w = 0; w < nWidths; w++ )
	    {
	      counts[j] = posCounts[n][t][w];
	      j++;
	    }
	}
    }

  for( seq = 0; seq < IP->nNumSequences; seq++ )
    {
      for( k = 0; k <= RP->nMaxBlocks; k++ )
	{
	  counts[j] = RP->kCounts[seq][k];
	  j++;
	}	
    }

  Gibbs_MPI_Send( B, counts, len, MPI_INT, 0, G_MPI_BAYES_DATA, MPI_COMM_WORLD );

  counts[0] = RP->acceptCount;
  counts[1] = RP->bayesSampleCount;
  Gibbs_MPI_Send( B, counts, 2, MPI_INT, 0, G_MPI_BAYES_MH, MPI_COMM_WORLD );
  RP->acceptCount = 0;
  RP->bayesSampleCount = 0;

  free( counts );

  if( ! IP->is_defined[cl_no_cred] )
    {
      len =  3 * IP->bayesSamplePeriod * RP->nMaxBlocks;
      NEW( counts, len, int );
      
      IP->Datafiles->mpiSample_fpt = fopen( IP->Datafiles->mpiSampleFileName, "a" );

      for( seq = 0; seq < IP->nNumSequences; seq++ )
	{
	  for( j = 0, n = 0; n < IP->bayesSamplePeriod; n++ )
	    {
	      for( k = 0; k < RP->nMaxBlocks; k++ )
		{
		  counts[j] = RP->samples[seq][seed_run-1][n][k].pos;
		  counts[j+1] = RP->samples[seq][seed_run-1][n][k].width;
		  counts[j+2] = RP->samples[seq][seed_run-1][n][k].t;
		  j+= 3;
		}
	    }
	  fwrite( counts, sizeof( int ), len, IP->Datafiles->mpiSample_fpt );
	  fflush( IP->Datafiles->mpiSample_fpt );
	}
      
      fsync( fileno( IP->Datafiles->mpiSample_fpt ) );
      fclose( IP->Datafiles->mpiSample_fpt );
      
      free( counts );
    }
}
#endif


void PrintSampleResults( Model B, int ***counts, int seed_run )
{
  IPtype  IP;
  RPType  RP;
  int     n;
  int     t;
  int     seq;
  int     start;
  int     pos;
  int     minWidth;
  int     maxWidth;
  int     nWidths;
  int     w;
  FILE    *fpt;
  int     k;
  int     total;
   
  IP = B->IP;
  RP = B->RP;

  minWidth = MinPossibleWidth( B );
  maxWidth = MaxPossibleWidth( B );
  nWidths = NumWidths( B );

  fpt = IP->Datafiles->out_fpt;
  
  fprintf( fpt, "+++++++++++++++++++++++++ Samples seed %d +++++++++++++++++++++++++\n",
	   seed_run );
  for(n = 0; n < IP->nSeqLen; n++) 
    {
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  for( total = 0, w = 0; w < nWidths; w++ )
	    total += counts[n][t][w];
	  
	  if( total > 0 )
	    {
	      seq = SequenceFromPosition( B, n );
	      start = SequenceStartPos( B, seq );
	      pos = n - start;
	      fprintf( fpt, "%d\t%d\t%d\t%d\t%d\t%g", 
		       seed_run, t, seq + 1, pos + 1, 
		       total, 
		       ((double) total/IP->bayesSamplePeriod ) ); 
	      for( w = 0; w < nWidths; w++ )
		{
		  if( counts[n][t][w] )
		    fprintf( fpt, "\t%d\t%d", w + minWidth, counts[n][t][w] );
		}
	      fprintf( fpt, "\n" );
	    }
	}
    }

  fprintf( fpt, "\n+++++++++++++++++++++++++ Sites seed %d +++++++++++++++++++++++++\n",
	   seed_run );
  for( seq = 0; seq < IP->nNumSequences; seq++ )
    {
      fprintf( fpt, "%d\t%d", seed_run, seq + 1 );
      for( k = 0; k <= RP->nMaxBlocks; k++ )
	{
	  fprintf( fpt, "\t%d\t%d\t%g", k, RP->kCounts[seq][k], ((double)  RP->kCounts[seq][k])/IP->bayesSamplePeriod );
	}
      fprintf( fpt, "\n" );
    }
}


void CopyCounts( Model B, int ***counts, int ***occurence, int seed_run, int iter )
{
  IPtype  IP;
  RPType  RP;
  int     n;
  int     t;
  int     width;
  int     seq;
  int     prevSeq = -1;
  int     k = 0;
  int     pos;
  
  IP = B->IP;
  RP = B->RP;

  for(n = 0; n < IP->nSeqLen; n++) 
    {
      seq = SequenceFromPosition( B, n );
      if( seq != prevSeq )
	k = 0;
      prevSeq = seq;
      
      for(t = 0; t < IP->nNumMotifTypes; t++) 
	{
	  width = MotifWidth( B, t ) - MinPossibleWidth( B );
	  counts[n][t][width] += occurence[n][t][0];
	  if( ! IP->is_defined[cl_no_cred] && occurence[n][t][0] )
	    {
	      pos = n - SequenceStartPos( B, seq );
	      RP->samples[seq][seed_run-1][iter][k].pos = pos;
	      RP->samples[seq][seed_run-1][iter][k].width = MotifWidth( B, t );
	      RP->samples[seq][seed_run-1][iter][k].t = t;
	      k++;
	    }
	}
    }
}


void CopyTotalCounts( Model B, int ***counts, int **kCounts )
{
  IPtype  IP;
  RPType  RP;
  int     n;
  int     t;
  int     k;
  int     nWidths;
  int     w;
  int     seq;
  int     nLen;
  int     start;
  
  IP = B->IP;
  RP = B->RP;

  nWidths = NumWidths( B );

  for( seq = 0; seq < IP->nNumSequences; seq++ )
    {
      for( k = 0; k <= RP->nMaxBlocks; k++ )	
	RP->kTotalCounts[seq][k] += kCounts[seq][k];

      nLen = SequenceLength( B, seq );
      start = SequenceStartPos( B, seq );
      for(n = 0; n < nLen; n++) 
	{
	  for(t = 0; t < IP->nNumMotifTypes; t++) 
	    {
	      for( w = 0; w < nWidths; w++ )
		RP->totalPosCounts[seq][n][w] += counts[n+start][t][w];
	    }
	}
    }
}

