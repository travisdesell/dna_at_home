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
#include "gibbsmpi.h"

#define _MPI_DEBUG_ 1

int TempCompare( const void *p1, const void *p2 );   /* reverse sort */
void SwapTemp( MPITemp *p1, MPITemp *p2 );
void MPIPrintBayesCounts( Model M, FILE *fpt, int *nCounts );
void MPICopyTotalCounts( Model M, int *nCounts, int accumLen );

void Gibbs_MPI( int argc, char **argv )
{
  int         i;
  int         n;
  Model       M;
  PoSition    **Pos;
  MaxResults  maxData, optmax;
  SimTime     S;
  int         j, k, t, m;
  FILE        *fpt;
  int         nPos;
  
  int         my_rank;       /* rank of process      */
  int         p;             /* number of processes  */
  int         source;        /* rank of sender       */
  double      maxtmp;        /* random number        */
  char        message[G_MPI_MESSAGE_SIZE];  /* storage for message  */
  char        name[100];     /* name of computer     */
  int         name_length;   /* computer name length */
  MPI_Status  status;        /* return status for    */
  /* receive              */
  int         doit_processor; /* rank of the chosen processor */
  int         *nCounts;       /* to hold accumulated count data */
  int         accumLen;
  int         cnt;
  int         sum;
  struct tms  tBuffer1;
  struct tms  tBuffer2;
  int         newlen;
  int         nSeedSent;
  int         nSeedRecv;
  double      dSeedVal;
  int         **nNumMotifs;   /* BT 06/22/2000 */
  double      dTempInfo[4];
  long        lRandInfo[2];
  int         nBusyProcessors;
  MPITemp     **mpiTemp;
  int         mpiTempReceived;
  int         nMPICountsReceived;
  double      ratio;
  struct tm   *currTime;
  struct tm   tm_buf;
  time_t      tloc;
  double      *dProcessorTemp;
  long        nRandSeed;
  int         ****bayesCounts = 0;
  int         bayesRecv = 0;
  int         ***bayeskCounts = 0;
  int         nWidths = 1;
  int         bayesSize = 0;
  int         w;
  char        *sampleFiles = 0;
  FILE        *sample_fpt;
  int         filesize = 0;
  int         bytesread;

  /* Find out process rank  */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  /* get my name */
  MPI_Get_processor_name(name,&name_length);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  
  if( strcmp( argv[1], "-PBernoulli" ) == 0 || strcmp( argv[1], "-PMBernoulli" ) == 0 )
    {
      GUI_argc=argc-1;
      for(i=1;i<=GUI_argc-1;i++)
	{
	  GUI_argv[i]=argv[i+1];
	}
      GUI_argv[0]=argv[0];
    }
  else
    {
      GUI_argc=argc;
      for(i=0;i<=GUI_argc-1;i++)
	{
	  GUI_argv[i]=argv[i];
	}
    }
  
  /* call the main function of Bernoulli */
  M = alloc_Model();
  M->IP->nMPI = TRUE;
  M->IP->nRank = my_rank;
  M->IP->nMPIProcesses = p;
  
  times( &tBuffer1);     /* 12/15/99 */
  
  if( my_rank == 0 )
    {
      BeginTime(&S);
      stripargs(GUI_argc, GUI_argv, M, FALSE); 
      fflush( M->IP->Datafiles->out_fpt ); /* BT 12/15/99 */
      sprintf( message, "Initialize" );
      MPI_Bcast(message, G_MPI_MESSAGE_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
  else
    {
      MPI_Bcast(message, G_MPI_MESSAGE_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
      stripargs(GUI_argc, GUI_argv, M, FALSE); 
    }
  
  NEW( M->IP->Datafiles->mpiTempFileName, FILENAME_MAX, char);   /* BT 1/27/97 */
  sprintf( M->IP->Datafiles->mpiTempFileName, "temp.process.%d.%d.out", my_rank, getpid() );
  remove( M->IP->Datafiles->mpiTempFileName );
#ifndef _RPI_
  if( ! M->IP->is_defined[cl_Z] )
    {
      MPE_IO_Stdout_to_file( M->IP->Datafiles->mpiTempFileName, O_TRUNC | O_RDWR | O_CREAT ); 
      chmod( M->IP->Datafiles->mpiTempFileName, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH );   /* BT 1/10/2000 */
      M->IP->Datafiles->mpiTemp_fpt = stdout;
    }
  else
    {
      M->IP->Datafiles->mpiTemp_fpt = fopen( M->IP->Datafiles->mpiTempFileName, "w" );
    }
#else
  M->IP->Datafiles->mpiTemp_fpt = fopen(M->IP->Datafiles->mpiTempFileName, "w" );
#endif

  if( M->IP->is_defined[cl_bayes] )
    {
      NEW( M->IP->Datafiles->mpiSampleFileName, FILENAME_MAX, char);   /* BT 1/27/97 */
      sprintf( M->IP->Datafiles->mpiSampleFileName, "sample.%d.%d.out", my_rank, getpid() );
      remove( M->IP->Datafiles->mpiSampleFileName );

      if( my_rank == 0 )
	{
	  NEW( sampleFiles, M->IP->nMPIProcesses * FILENAME_MAX, char );
	  for( source = 1; source < M->IP->nMPIProcesses; source++ )
	    {
	      Gibbs_MPI_Recv( M, &sampleFiles[source*FILENAME_MAX], FILENAME_MAX, MPI_CHAR, source, 
			      G_MPI_BAYES_SAMPLE, MPI_COMM_WORLD, &status );
	      PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, "Sample file: %d %s\n",
			    source, &sampleFiles[source*FILENAME_MAX] );
	    }
	      
	}
      else
	Gibbs_MPI_Send( M, M->IP->Datafiles->mpiSampleFileName, FILENAME_MAX, MPI_CHAR, 0, G_MPI_BAYES_SAMPLE, 
			MPI_COMM_WORLD );
    }

  PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, "Processor %d id is %s\n", my_rank,name);  /* BT 1/10/2000 */
  
  get_inputs(M);
  M->IP->is_defined[cl_u] = FALSE;    /* override suboptimal output */
  M->IP->is_defined[cl_t] = FALSE;    /* 12/15/99 */

  /* M->IP->lSeedVal = (long)time(NULL) + 1000 * my_rank; */  /* Get a new seed just in case */
  if( my_rank == 0 )
    {
      M->IP->lSeedVal0 = M->IP->lSeedVal;
      for( source = 1; source < M->IP->nMPIProcesses; source++ )
	{	      
	  nRandSeed = M->IP->lSeedVal + 1009 * source;
	  lRandInfo[0] = nRandSeed;  /* BT 02/25/03 */
	  lRandInfo[1] = M->IP->lSeedVal;
	  Gibbs_MPI_Send( M, lRandInfo, 2, MPI_LONG, source, G_MPI_RAND_SEED, MPI_COMM_WORLD );
	}
    }
  else
    {
      Gibbs_MPI_Recv( M, lRandInfo, 2, MPI_LONG, 0, G_MPI_RAND_SEED, MPI_COMM_WORLD,
		      &status );
      M->IP->lSeedVal = lRandInfo[0];
      M->IP->lSeedVal0 = lRandInfo[1];	  
    }

  PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, "seed = %d\n", M->IP->lSeedVal );    
  
  if( M->IP->is_defined[cl_l] ) 
    Wilcoxon( M );
  else
    get_strings(M);           /* Read the sequences */
  
  InitRProbStruct( M );
  
  alloc_Counts(M);
  set_counts(M);                        /* Set the initial counts     */
  set_posterior_prob(M->IP, M->First);  /* Set initial posterior prob */
  *M->Seq->ProcessedSTR=Xnu_Sequence(M);/* preprocess the sequence    */

  set_pseudo_counts(M, M->IP, M->First);   /* Set initial pseudo counts  */

  if(M->IP->nAlphaLen == 20)   
    newlen = M->IP->nSeqLen - M->IP->nNumProcessed;
  else
    newlen = M->First->nTotBack;
  PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, "newlen %d IP.nSeqLen %d\n", newlen, M->IP->nSeqLen);   

  if( M->IP->is_defined[cl_B] )
    ReadBkgndComposition( M, M->First );
  
  if( M->IP->is_defined[cl_U] )
    ReadSpacingFile( M, M->IP->Datafiles->spacing_fpt );
  
  if(M->IP->Datafiles->weight_fpt != NULL)      
    ReadWeightFile( M );
  
  CalcAlignProbFromFootprint( M );
  
  /* allocate space for alignment */
  /*  Pos[MotifType][SeqLen] */
  NEWP(Pos,M->IP->nNumMotifTypes,PoSition);
  for(t=0;t<M->IP->nNumMotifTypes;t++)
    NEW(Pos[t], M->IP->nSeqLen, PoSition);
  
  nNumMotifs = copy_motif_num(M->IP);   /* Save the input motif count */   /* BT 06/22/2000 */
  
  fpt = M->IP->Datafiles->out_fpt;
  init_maxdata(&maxData);
  init_maxdata(&optmax);
  
  NEW( dProcessorTemp, M->IP->nMPIProcesses, double );
  NEWP( mpiTemp, M->IP->nMPIProcesses, MPITemp );
  for( i = 0; i < M->IP->nMPIProcesses; i++ )
    {
      NEW( mpiTemp[i], 1, MPITemp );
    }
  
  if( M->IP->is_defined[cl_bayes] && my_rank == 0 )
    {
      nWidths = NumWidths( M );
      NEWP3( bayesCounts, M->IP->nSeeds, int );
      NEWPP( bayeskCounts, M->IP->nSeeds, int );
      for( i = 0; i < M->IP->nSeeds; i++ )
	{
	  bayesSize = 0;
	  NEWPP( bayesCounts[i], M->IP->nSeqLen, int);
	  for( n = 0; n < M->IP->nSeqLen; n++) 
	    {
	      NEWP(bayesCounts[i][n], M->IP->nNumMotifTypes, int);
	      for(t = 0; t < M->IP->nNumMotifTypes; t++) 
		{
		  NEW( bayesCounts[i][n][t], nWidths, int );
		  bayesSize += nWidths;
		}
	    }

	  NEWP( bayeskCounts[i], M->IP->nNumSequences, int );
	  for( n = 0; n < M->IP->nNumSequences; n++ )
	    NEW( bayeskCounts[i][n], M->RP->nMaxBlocks + 1, int );
	}
    }

  if( my_rank > 0 )
    maxData = suboptSampler(M, Pos, &S);    /* Find Suboptimal Solution   */
  else
    {
      nSeedSent = 0;
      nBusyProcessors = 0;
      for( source = 1; source < M->IP->nMPIProcesses && nSeedSent < M->IP->nSeeds; source++ )
	{
	  if( M->IP->is_defined[cl_X] && M->AN->bExchange )
	    {
	      dTempInfo[0] = 1.0 / ((1.0 / M->AN->dMinTemp) - 
				    ((((double) source) - 1.0) / (((double) M->IP->nMPIProcesses) - 2.0)) *
				    ((1.0 / M->AN->dMinTemp) - (1.0 / M->AN->dMaxTemp)));
	    } 
	  dTempInfo[2] = nSeedSent;
	  Gibbs_MPI_Send( M, dTempInfo, 4, MPI_DOUBLE, source, G_MPI_SEED, MPI_COMM_WORLD );
	  nSeedSent++;
	  nBusyProcessors++;
	}
      
      maxtmp = -DBL_MAX;
      doit_processor = -1;
      nSeedRecv = 0;
      mpiTempReceived = 0;
      while( nSeedRecv < M->IP->nSeeds )
	{
	  Gibbs_MPI_Recv( M, dTempInfo, 4, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
			  &status );
	  source = status.MPI_SOURCE;
	  if( status.MPI_TAG == G_MPI_DATA )
	    {
	      dSeedVal = dTempInfo[1];
	      PrintTempOut(M->IP->Datafiles->mpiTemp_fpt, 
			   "Seed %d received with value %f %d motifs from %d\n", 
			   nSeedRecv, dSeedVal, (int) dTempInfo[3], source);	  
	      nSeedRecv++;
	      nBusyProcessors--;
	      if( dSeedVal > maxtmp )
		{
		  maxtmp = dSeedVal;
		  doit_processor = source;
		}
	      
	      if( nSeedSent < M->IP->nSeeds )
		{
		  dTempInfo[2] = nSeedSent;
		  Gibbs_MPI_Send( M, dTempInfo, 4, MPI_DOUBLE, source, G_MPI_SEED, MPI_COMM_WORLD );
		  nSeedSent++;
		  nBusyProcessors++;
		}
	      PrintTempOut(M->IP->Datafiles->mpiTemp_fpt, "%d processors are busy\n", 
			   nBusyProcessors );
	    }

	  if( status.MPI_TAG == G_MPI_BAYES )
	    {
	      accumLen = bayesSize + (M->IP->nNumSequences * (M->RP->nMaxBlocks + 1));
	      NEW( nCounts, accumLen, int );

	      dTempInfo[0] = accumLen;
	      Gibbs_MPI_Send( M, dTempInfo, 4, MPI_DOUBLE, source, G_MPI_BAYES, MPI_COMM_WORLD );
	      Gibbs_MPI_Recv( M, nCounts, accumLen, MPI_INT, source, G_MPI_BAYES_DATA,
			      MPI_COMM_WORLD, &status);

	      nWidths = NumWidths( M );
	      for( j = 0, n = 0; n < M->IP->nSeqLen; n++) 
		{
		  for(t = 0; t < M->IP->nNumMotifTypes; t++) 
		    {
		      for( w = 0; w < nWidths; w++ )
			{
			  bayesCounts[bayesRecv][n][t][w] = nCounts[j];
			  j++;
			}
		    }
		}

	      for( n = 0; n < M->IP->nNumSequences; n++ )
		{
		  for( k = 0; k <= M->RP->nMaxBlocks; k++ )
		    {
		      bayeskCounts[bayesRecv][n][k] = nCounts[j];
		      j++;
		    }
		}

	      bayesRecv++;
	      
	      Gibbs_MPI_Recv( M, nCounts, 2, MPI_INT, source, G_MPI_BAYES_MH,
			      MPI_COMM_WORLD, &status);
	      M->RP->acceptCount += nCounts[0];
	      M->RP->bayesSampleCount += nCounts[1];

	      free( nCounts );
	    }

	  if( M->IP->is_defined[cl_X] && 
	      (status.MPI_TAG == G_MPI_TEMP || mpiTempReceived == nBusyProcessors) ) 
	    {
	      if( status.MPI_TAG == G_MPI_TEMP )
		{
		  time( &tloc );
		  currTime = localtime_r( &tloc, &tm_buf );
		  PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, 
				"%d has sent temp = %f prob = %f seed = %d at %d:%d:%d\n",
				source, dTempInfo[0], dTempInfo[1], (int) dTempInfo[2],
				currTime->tm_hour, currTime->tm_min, currTime->tm_sec );
		  mpiTemp[mpiTempReceived]->temp = dTempInfo[0];
		  mpiTemp[mpiTempReceived]->currProb = dTempInfo[1];
		  mpiTemp[mpiTempReceived]->seed = dTempInfo[2];
		  mpiTemp[mpiTempReceived]->source = source;
		  mpiTempReceived++;
		}
	      
	      if( mpiTempReceived == nBusyProcessors )
		{
		  PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, 
				"Sorting...\n" );		      
		  qsort( mpiTemp, nBusyProcessors, sizeof( MPITemp * ), TempCompare );
		  
		  for( j = 0; j < EXCHNG_COUNT; j++ )
		    {
		      for( i = nBusyProcessors - 1; i > 0; i-- )
			{
			  ratio = exp( (mpiTemp[i]->currProb - mpiTemp[i-1]->currProb) * 
				       ((1.0 / mpiTemp[i-1]->temp) - (1.0 / mpiTemp[i]->temp)) );
			  PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, 
					"comparing %d temp = %f prob = %f to %d temp = %f prob = %f ratio = %f\n",
					mpiTemp[i-1]->source, mpiTemp[i-1]->temp, mpiTemp[i-1]->currProb,
					mpiTemp[i]->source, mpiTemp[i]->temp, mpiTemp[i]->currProb, 
					ratio );
			  if( ratio > drand() )
			    {
			      PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, "swapping %d and %d\n",
					    mpiTemp[i-1]->source, mpiTemp[i]->source );
			      SwapTemp( mpiTemp[i], mpiTemp[i-1] );
			    }
			}		
		    }

		  for( i = 0; i < nBusyProcessors;  i++ )
		    {
		      time( &tloc );
		      currTime = localtime_r( &tloc, &tm_buf );
		      PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, 
				    "Sending temp = %f prob = %f seed = %d to %d at  %d:%d:%d\n",
				    mpiTemp[i]->temp, mpiTemp[i]->currProb, (int) mpiTemp[i]->seed, 
				    mpiTemp[i]->source,
				    currTime->tm_hour, currTime->tm_min, currTime->tm_sec );
		      dProcessorTemp[mpiTemp[i]->source] = mpiTemp[i]->temp;
		      Gibbs_MPI_Send( M, mpiTemp[i], 4, MPI_DOUBLE, mpiTemp[i]->source, 
				G_MPI_TEMP, MPI_COMM_WORLD );
		    }
		  mpiTempReceived = 0;
		}
	    } 
	}
      
      /* tell everyone to quit */
      for( source = 1; source < p ; source++ )
	{
	  dTempInfo[2] = nSeedSent;
	  Gibbs_MPI_Send( M, dTempInfo, 4, MPI_DOUBLE, source, G_MPI_DONE, MPI_COMM_WORLD );
	}
      
      PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, "The max is %f from processor %d\n", 
		    maxtmp, doit_processor );
    }
  
  if( my_rank == 0 )
    {	  
      MPI_Bcast( &doit_processor, 1, MPI_INT, 0, MPI_COMM_WORLD);

      if( M->IP->is_defined[cl_bayes] )
	{
	  accumLen = M->IP->nSeeds * (bayesSize + (M->IP->nNumSequences * (M->RP->nMaxBlocks + 1)));

	  NEW( nCounts, accumLen, int );
	  
	  nWidths = NumWidths( M );
	  for( j = 0, i = 0; i < M->IP->nSeeds; i++ )
	    {
	      for(n = 0; n < M->IP->nSeqLen; n++) 
		{
		  for( t = 0; t < M->IP->nNumMotifTypes; t++ )
		    {
		      for( w = 0; w < nWidths; w++ )
			{
			  nCounts[j] = bayesCounts[i][n][t][w];
			  j++;
			}
		    }
		}	      
	      
	      for( n = 0; n < M->IP->nNumSequences; n++ )
		{
		  for( k = 0; k <= M->RP->nMaxBlocks; k++ )
		    {
		      nCounts[j] = bayeskCounts[i][n][k];
		      j++;
		    }
		}
	    }

	  Gibbs_MPI_Send( M, &accumLen, 1, MPI_INT, doit_processor, G_MPI_BAYES_DATA, MPI_COMM_WORLD );
	  for(i = 0; i < accumLen; i += G_MPI_TRANSFER_SIZE)
	    {	      
	      Gibbs_MPI_Send( M, &nCounts[i], min(accumLen - i, G_MPI_TRANSFER_SIZE), MPI_INT, 
			      doit_processor, G_MPI_BAYES_DATA, MPI_COMM_WORLD );
	    }

	  for( source = 1; source < M->IP->nMPIProcesses; source++ )
	    Gibbs_MPI_Send( M, &sampleFiles[source*FILENAME_MAX], FILENAME_MAX, MPI_CHAR, 
			    doit_processor, G_MPI_BAYES_SAMPLE, 
			    MPI_COMM_WORLD );

	  nCounts[0] = M->RP->acceptCount;
	  nCounts[1] = M->RP->bayesSampleCount;
	  Gibbs_MPI_Send( M, nCounts, 2, MPI_INT, doit_processor, 
			       G_MPI_BAYES_MH, MPI_COMM_WORLD );

	  free( nCounts );
	}

      Gibbs_MPI_Recv( M, dTempInfo, 4, MPI_DOUBLE, MPI_ANY_SOURCE, G_MPI_DATA, MPI_COMM_WORLD,
		      &status );
      dSeedVal = dTempInfo[1];
    }           
  else 
    {
      MPI_Bcast( &doit_processor, 1, MPI_INT, 0, MPI_COMM_WORLD);
      
      if( my_rank == doit_processor )  /* have I been selected ? */ 
	{
	  /***************** near optimal processing ************/
	  PrintTempOut(M->IP->Datafiles->mpiTemp_fpt, "I, %s have been chosen \n", name);
	  
	  /* BT 06/22/2000 */
	  reset_motif_num(nNumMotifs, M->IP);   /* restore the motif count so it appears corrctly in output */
	  
	  printargs(GUI_argc, GUI_argv, M->IP, M->IP->Datafiles->mpiTemp_fpt );
	  print_options( M );
	  PrintSeqDescriptions( M ); 

	  if( M->IP->is_defined[cl_bayes] )
	    {
	      Gibbs_MPI_Recv( M, &accumLen, 1, MPI_INT, 0, G_MPI_BAYES_DATA,
			      MPI_COMM_WORLD, &status);
	      PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, "accumlen =  %d\n",
			    accumLen );
	  
	      NEW( nCounts, accumLen, int );
	      for(i = 0; i < accumLen; i += G_MPI_TRANSFER_SIZE)
		{	      
		  Gibbs_MPI_Recv( M, &nCounts[i], min(accumLen - i, G_MPI_TRANSFER_SIZE), MPI_INT, 
				  0, G_MPI_BAYES_DATA, MPI_COMM_WORLD, &status );
		}

	      MPICopyTotalCounts( M, nCounts, accumLen );

	      if( M->IP->is_defined[cl_bayes_counts] )
		MPIPrintBayesCounts( M, fpt, nCounts );
	      
	      free( nCounts );

	      if( ! M->IP->is_defined[cl_no_cred] )
		{
		  accumLen = 3 * M->IP->bayesSamplePeriod * M->RP->nMaxBlocks;
		  NEW( nCounts, accumLen, int );
		  NEW( sampleFiles, FILENAME_MAX, char );
		  
		  for( m = 0, source = 1; source < M->IP->nMPIProcesses; source++ )
		    {
		      Gibbs_MPI_Recv( M, sampleFiles, FILENAME_MAX, MPI_CHAR, 0, 
				      G_MPI_BAYES_SAMPLE, MPI_COMM_WORLD, &status );
		      
		      PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, "Sample file: %d %s seed: %d\n",
				    source, sampleFiles, m );
		      
		      if( ! FileExists( sampleFiles ) )
			{
			  PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, "Sample file: %d %s seed: %d was not found.\n",
					source, sampleFiles, m );
			  continue;
			}

		      filesize = FileSize( sampleFiles ) / sizeof( int );
		      bytesread = 0;		      
		      sample_fpt = fopen( sampleFiles, "r" );
		      while( bytesread < filesize )
			{
			  for( i = 0; i < M->IP->nNumSequences; i++ )
			    {
			      sum = fread( nCounts, sizeof( int ), accumLen,  sample_fpt );
			      bytesread += sum;
			      
			      PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, 
					    "Read %d ints from Sample file: %d %s seed: %d seq: %d\n",
					    sum, source, sampleFiles, m, i );
		      
			      for( j = 0, n = 0; n < M->IP->bayesSamplePeriod; n++ )
				{
				  for( k = 0; k < M->RP->nMaxBlocks; k++ )
				    {
				      M->RP->samples[i][m][n][k].pos = nCounts[j];
				      M->RP->samples[i][m][n][k].width = nCounts[j+1];
				      M->RP->samples[i][m][n][k].t = nCounts[j+2];
				      j += 3;
				    }
				}
			    }
			  m++;
			}		      
		      PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, "Read %d ints from Sample file: %d %s seed: %d\n",
				    bytesread, source, sampleFiles, m );
		      fclose( sample_fpt );
		      remove( sampleFiles ); 
		    }

		  free(sampleFiles );
		  
		  Gibbs_MPI_Recv( M, nCounts, 2, MPI_INT, 0, 
				  G_MPI_BAYES_MH, MPI_COMM_WORLD, &status );
		  M->RP->acceptCount = nCounts[0];
		  M->RP->bayesSampleCount = nCounts[1];

		  free( nCounts );
		}
	    }
	  
	  if( ! M->IP->is_defined[cl_nopt] )
	    {
	      PrintTempOut(fpt, "\n\n\n\n\n",
			   "======================================================================\n");
	      PrintTempOut(fpt, 
			   "======================== NEAR OPTIMAL RESULTS ========================\n");
	      PrintTempOut(fpt, 
			   "======================================================================\n",
			   "\n");
	    }
	  
	  if( M->IP->is_defined[cl_opt] )
	    optmax = NearOptSampler(M, maxData, Pos);   /* Find Near optimal Solution */
	  
	  if(M->IP->is_defined[cl_opt] && ! M->IP->is_defined[cl_m] )
	    {             /* Maximize the MAP value */
	      PrintTempOut(fpt, "\n\n\n\n\n",
			   "======================================================================\n");
	      PrintTempOut(fpt, 
			   "======================== MAP MAXIMIZATION RESULTS ====================\n");
	      PrintTempOut(fpt,
			   "======================================================================\n",
			   "\n");
	      optmax = maximize_map(optmax, M, Pos);   /* BT 2/21/97 */
	    }
	      
	  if(  M->IP->is_defined[cl_bayes] )
	    {
	      fprintf(fpt, "\n\n\n\n\n======================================================================\n" );
	      fprintf(fpt, 		   
		      "========================== CENTROID RESULTS ==========================\n" );
	      fprintf(fpt, 		   
		      "======================================================================\n\n" );
	      Centroid( M );
	    }

	  times( &tBuffer2);     /* 12/15/99 */
	  PrintTempOut( fpt, "Elapsed time: %f secs\n", 
			(double) (tBuffer2.tms_utime - tBuffer1.tms_utime) / M->IP->ticks );	
	  
	  dTempInfo[0] = 0;
	  dTempInfo[1] = optmax.dProbability;
	  Gibbs_MPI_Send( M, dTempInfo, 4, MPI_DOUBLE, 0, G_MPI_DATA, MPI_COMM_WORLD );
	}
    }
  
  if( M->IP->is_defined[cl_Q] )    
    {
      for( accumLen = 0, i = 0; i < M->IP->nNumSequences; i++ )
	accumLen += SequenceLength( M, i );
      accumLen *= M->IP->nNumMotifTypes;
      NEW( nCounts, accumLen, int );
      
      if( my_rank == 0 )
	{	      
	  nMPICountsReceived = 0;
	  while( nMPICountsReceived < M->IP->nMPIProcesses )
	    {
	     Gibbs_MPI_Recv( M, nCounts, accumLen, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
			     MPI_COMM_WORLD, &status);		  
	      if( status.MPI_TAG == G_MPI_COUNT_DATA )
		{
		  PrintTempOut(M->IP->Datafiles->mpiTemp_fpt, 
			       "received %d optimal counts from %d final temp = %f\n", 
			       accumLen, status.MPI_SOURCE, dProcessorTemp[status.MPI_SOURCE]);
		  for( cnt = 0, i = 0; i < M->IP->nNumSequences; i++ )
		    {
		      for( j = 0; j < SequenceLength( M, i ); j++ )
			{
			  for( t = 0;  t < M->IP->nNumMotifTypes; t++ )
			    {
			      M->IP->nAlignCnts[1][t][i][j] += nCounts[cnt];
			      if( dProcessorTemp[status.MPI_SOURCE] == 1.0 )
				M->IP->nAlignCnts[2][t][i][j] += nCounts[cnt];
			      cnt++;
			    }
			}
		    }
		}
	      else if( status.MPI_TAG == G_MPI_COUNT_DATA0 )
		{
		  PrintTempOut(M->IP->Datafiles->mpiTemp_fpt, "received %d near optimal counts from %d\n", 
			       accumLen, status.MPI_SOURCE );
		  for( cnt = 0, i = 0; i < M->IP->nNumSequences; i++ )
		    {
		      for( j = 0; j < SequenceLength( M, i ); j++ )
			{
			  for( t = 0;  t < M->IP->nNumMotifTypes; t++ )
			    {
			      M->IP->nAlignCnts[0][t][i][j] += nCounts[cnt++];
			    }
			}
		    }
		}
	      nMPICountsReceived++;
	    }
	  
	  PrintTempOut( M->IP->Datafiles->occur_fpt, "seq pos pos2 motif count\n" );
	  PrintTempOut( M->IP->Datafiles->near_occur_fpt, "seq pos pos2 motif count\n" );
	  for( i = 0; i < M->IP->nNumSequences; i++ )
	    {
	      for( j = 0; j < SequenceLength( M, i ); j++ )
		{
		  nPos = SequenceStartPos( M, i ) + j;
		  for( t = 0;  t < M->IP->nNumMotifTypes; t++ )
		    {
		      PrintTempOut( M->IP->Datafiles->occur_fpt, 
				    "%5d %5d %5d %5d %5d\n",
				    i, j, nPos, t, M->IP->nAlignCnts[1][t][i][j] );
		      
		      /* Changed from M->IP->nAlignCnts[2][t][i][j] so we can print near opt counts */
		      PrintTempOut( M->IP->Datafiles->near_occur_fpt, 
				    "%5d %5d %5d %5d %5d\n",
				    i, j, nPos, t,  M->IP->nAlignCnts[0][t][i][j] );  /* BT 08/23/05 */
		    }
		}
	    }
	}
      else
	{
	  PrintTempOut(M->IP->Datafiles->mpiTemp_fpt, "%d is copying %d counts\n", 
		       my_rank, accumLen);	

	  for( cnt = 0, i = 0; i < M->IP->nNumSequences; i++ )
	    {
	      for( j = 0; j < SequenceLength( M, i ); j++ )
		{
		  for( t = 0;  t < M->IP->nNumMotifTypes; t++ )
		    {
		      for( sum = 0, k = 1; k <= M->IP->nSeeds; k++ )
			{
			  sum += M->IP->nAlignCnts[k][t][i][j];
			}
		      nCounts[cnt++] += sum;
		    }
		}
	    }
	  
	  PrintTempOut(M->IP->Datafiles->mpiTemp_fpt, "%d is sending %d counts\n", my_rank, cnt);
	  Gibbs_MPI_Send( M, nCounts, accumLen, MPI_INT, 0, G_MPI_COUNT_DATA, 
			  MPI_COMM_WORLD );
	  
	  if( my_rank == doit_processor )
	    {
	      PrintTempOut(M->IP->Datafiles->mpiTemp_fpt, "%d is copying %d near optimal counts\n", 
			   my_rank, accumLen);	
	      for( cnt = 0, i = 0; i < M->IP->nNumSequences; i++ )
		{
		  for( j = 0; j < SequenceLength( M, i ); j++ )
		    {
		      for( t = 0;  t < M->IP->nNumMotifTypes; t++ )
			{
			  nCounts[cnt++] = M->IP->nAlignCnts[0][t][i][j];
			}
		    }
		}
	      
	      PrintTempOut(M->IP->Datafiles->mpiTemp_fpt, "%d is sending %d near optimal counts\n", 
			   my_rank, cnt);
	      Gibbs_MPI_Send( M, nCounts, accumLen, MPI_INT, 0, G_MPI_COUNT_DATA0, 
			      MPI_COMM_WORLD );
	    }
	}
      free( nCounts );	      
    }

  EndTime(M, &S, M->IP->Datafiles->mpiTemp_fpt);
  
  free_maxdata(&maxData, M->IP);
  free_maxdata(&optmax, M->IP);
  
  FREEP(Pos,M->IP->nNumMotifTypes);  
  free(Pos);
  
  free( dProcessorTemp );
  FREEP( mpiTemp, M->IP->nMPIProcesses );
  
  if( M->IP->is_defined[cl_bayes] && my_rank == 0 )
    {
      for( i = 0; i < M->IP->nSeeds; i++ )
	FREEPP(bayesCounts[i], M->IP->nSeqLen, M->IP->nNumMotifTypes);      
      FREEPP( bayeskCounts, M->IP->nSeeds, M->IP->nNumSequences );
      free( sampleFiles );
    }
  
  fclose( M->IP->Datafiles->mpiTemp_fpt );
  FreeData(M); 
}


int TempCompare( const void *p1, const void *p2 )
{
  MPITemp *t1;
  MPITemp *t2;
 
  t1 = *(MPITemp **) p1;
  t2 = *(MPITemp **) p2;

  if( t1->temp > t2->temp )
    return( 1 );
  if( t1->temp < t2->temp )
    return( -1 );
  return( 0 );
}


void SwapTemp( MPITemp *p1, MPITemp *p2 )
{
  MPITemp tp;

  tp.currProb = p1->currProb;
  tp.seed = p1->seed;
  tp.source = p1->source;
  p1->currProb = p2->currProb;
  p1->seed = p2->seed;
  p1->source = p2->source;
  p2->currProb = tp.currProb;
  p2->seed = tp.seed;
  p2->source = tp.source;
}


int Gibbs_MPI_Send( Model B, void *buf, int count, MPI_Datatype datatype, int dest, 
		    int tag, MPI_Comm comm )
{
  int   ret;
#ifdef _MPI_DEBUG_
  int  dataSize;

  MPI_Type_size( datatype, &dataSize );
  PrintTempOut(B->IP->Datafiles->mpiTemp_fpt, "Sending %d elements of size %d to %d tag: %d\n", 
	       count, dataSize, dest, tag);
#endif

  ret = MPI_Send( buf, count, datatype, dest, tag, comm );

#ifdef _MPI_DEBUG_
  PrintTempOut(B->IP->Datafiles->mpiTemp_fpt, "ret = %d\n", ret );
#endif

  return ret;

}


int Gibbs_MPI_Recv( Model B, void *buf, int count, MPI_Datatype datatype, int source, 
		    int tag, MPI_Comm comm, MPI_Status *status )
{
  int  ret;
#ifdef _MPI_DEBUG_
  int  dataSize;

  MPI_Type_size( datatype, &dataSize );
#endif

  ret = MPI_Recv( buf, count, datatype, source, tag, comm, status );

#ifdef _MPI_DEBUG_
  PrintTempOut(B->IP->Datafiles->mpiTemp_fpt, "Received %d elements of size %d from %d tag: %d ret: %d\n", 
	       count, dataSize, source, tag, ret);
#endif

  return ret;
}


void MPIPrintBayesCounts( Model M, FILE *fpt, int *nCounts )
{
  int   i;
  int   j;
  int   n;
  int   total;
  int   t;
  int   w;
  int   minWidth;
  int   seq;
  int   start;
  int   nPos;
  int   nWidths;
  int   k;
  
  nWidths = NumWidths( M );
  minWidth = MinPossibleWidth( M );

  for( j = 0, i = 0; i < M->IP->nSeeds; i++ )
    {
      fprintf( fpt, "+++++++++++++++++++++++++ Samples seed %d +++++++++++++++++++++++++\n",
	       i + 1 );
		      
      for(n = 0; n < M->IP->nSeqLen; n++) 
	{
	  seq = SequenceFromPosition( M, n );
	  start = SequenceStartPos( M, seq );
	  nPos = n - start;

	  for( t = 0; t < M->IP->nNumMotifTypes; t++ )
	    {	      
	      for( total = 0, w = 0; w < nWidths; w++ )
		total += nCounts[j+w];
	      
	      if( total > 0 )
		{
		  fprintf( fpt, "%d\t%d\t%d\t%d\t%d\t%g", 
			   i+1, t, seq + 1, nPos + 1, 
			   total,
			   ((double) total)/ M->IP->bayesSamplePeriod ); 
		  
		  for( w = 0; w < nWidths; w++ )
		    {
		      if( nCounts[j+w] )
			fprintf( fpt, "\t%d\t%d", w + minWidth, nCounts[j+w] );
		    }
		  fprintf( fpt, "\n" );
		}
	      j += nWidths;
	    }
	}
		      
      fprintf( fpt, "\n+++++++++++++++++++++++++ Sites seed %d +++++++++++++++++++++++++\n",
	       i + 1 );
      for( seq = 0; seq < M->IP->nNumSequences; seq++ )
	{
	  fprintf( fpt, "%d\t%d", i+1, seq + 1 );
	  for( k = 0; k <= M->RP->nMaxBlocks; k++ )
	    {
	      fprintf( fpt, "\t%d\t%d\t%g", k, nCounts[j], ((double) nCounts[j])/ M->IP->bayesSamplePeriod );
	      j++;
	    }
	  fprintf( fpt,  "\n" );
	}
    }  
}


void MPICopyTotalCounts( Model M, int *nCounts, int accumLen )
{
  IPtype  IP;
  RPType  RP;
  int     i;
  int     j;
  int     n;
  int     k;
  int     t;
  int     w;
  int     seq;
  int     nWidths;
  int     nLen;
  char    *msg;

  IP = M->IP;
  RP = M->RP;

  nWidths = NumWidths( M );

  for( seq = 0; seq < IP->nNumSequences; seq++ )
    {
      for( k = 0; k <= RP->nMaxBlocks; k++ )	
	RP->kTotalCounts[seq][k] = 0;

      nLen = SequenceLength( M, seq );
      for(n = 0; n < nLen; n++) 
	{
	  for( w = 0; w < nWidths; w++ )
	    RP->totalPosCounts[seq][n][w] = 0;
	}
    }

  for( j = 0, i = 0; i < IP->nSeeds; i++ )
    {
      for( seq = 0; seq < IP->nNumSequences; seq++ )
	{
	  nLen = SequenceLength( M, seq );
	  for(n = 0; n < nLen; n++) 
	    {
	      for(t = 0; t < IP->nNumMotifTypes; t++) 
		{
		  for( w = 0; w < nWidths; w++ )
		    {
		      if( j >= accumLen )
			{
			  NEW( msg, 1024, char );
			  sprintf( msg, "MPICopyTotalCounts: length error accumLen = %d seq = %d n = %d w= %d j = %d",
				   accumLen, seq, n, w, j );
			  p_internal_error( msg );
			}
		      RP->totalPosCounts[seq][n][w] += nCounts[j];
		      j++;
		    }
		}
	    }
	}

      for( seq = 0; seq < M->IP->nNumSequences; seq++ )
	{
	  for( k = 0; k <= M->RP->nMaxBlocks; k++ )
	    {
	      if( j >= accumLen )
		{
		  NEW( msg, 1024, char );
		  sprintf( msg, "MPICopyTotalCounts: length error accumLen = %d seq = %d k = %d j = %d",
			   accumLen, seq, k, j );
		  p_internal_error( msg );
		}
	      RP->kTotalCounts[seq][k] += nCounts[j];
	      j++;
	    }
	}
    }
}
