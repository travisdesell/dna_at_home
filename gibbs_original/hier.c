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
#include "hier.h"

#define CMD_ALLOC_SIZE  20
#define CMD_SIZE        1024
#define TOKEN_SIZE      256
#define K_BUFFER_SIZE   1024
#define MSG_BUFFER_SIZE 1024
#define TOL             1.0e-03
#define MAX_ITER        100
#define MIN_PSEUDO      1.0e-06

char** ReadCommandFile( char *commandFile, int *nCmds );
int ParseCommandLine( char *commandLine, char **argv );
void PrintReceivedData( Model M, double *buffer, MPI_Status *status, int seed, int iter );
void PrintMPIMsg( Model B, char *fmt, ...);
void SendMessageToAllChildren( Model M, double *dTempInfo, int dataLength, int tag, char *tagName, int nProceses  );
int MinimizeK( Model M, double *kPseudo, double **kCnts, int kmax, int kReceived );
int MinimizeK2( Model M, double *kPseudo, double **kCnts, int kmax, int kReceived );
void SendCounts( Model M, int seed, int iter );
void CalcProbK( Model B, int kmax, int kReceived, double *kCnts, double *kPseudo, double *kProb );
void SendProb( Model B, int kmax, double *kProb, int *done, int nProcesses );


void Gibbs_Hierarchy( int argc, char **argv )
{
  int         i;
  Model       M;
  PoSition    **Pos;
  MaxResults  maxData, optmax;
  SimTime     S;
  int         k, t;
  FILE        *fpt;
  
  int         my_rank;       /* rank of process      */
  int         p;             /* number of processes  */
  int         source;        /* rank of sender       */
  char        name[100];     /* name of computer     */
  int         name_length;   /* computer name length */
  MPI_Status  status;        /* return status for    */
                             /* receive              */
  struct tms  tBuffer1;
  struct tms  tBuffer2;
  int         newlen;
  int         **nNumMotifs;   /* BT 06/22/2000 */
  long        lRandInfo[2];
  int         nBusyProcessors;
  long        nRandSeed;
  char        **commandList = NULL;
  int         nCmds = 0;
  char        *cmdBuffer;
  int         *done;
  double      *kCnts = NULL;
  double      *kPseudo = NULL;
  double      *kBuffer;
  double      *kProb = NULL;
  int         kmax = 0;
  int         kReceived;
  int         kmaxPrev = -1;
  double      dTempInfo[4];
  int         totalProcessors;
  int         nSeeds;
  int         nSeedsPrev = -1;
  
  /* Find out process rank  */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  /* get my name */
  MPI_Get_processor_name(name,&name_length);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
 
  GUI_argc=argc;
  for(i=0;i<=GUI_argc-1;i++)
    GUI_argv[i]=argv[i];

  /* call the main function of Bernoulli */
  M = alloc_Model();
  M->IP->nMPI = TRUE;
  M->IP->nRank = my_rank;
  M->IP->nMPIProcesses = p;

  times( &tBuffer1); 

  if( my_rank == 0 )
    {
      BeginTime(&S);
      commandList = ReadCommandFile( GUI_argv[2], &nCmds );

      if( nCmds > M->IP->nMPIProcesses - 1 )
	p_error( "Hierarchical model: the number of command lines exceeds the number of processors" );
 
      for( i = 0; i < nCmds; i++ )
	Gibbs_MPI_Send( M, commandList[i], strlen( commandList[i] ) + 1, MPI_CHAR, i+1, 
			G_MPI_COMMAND_LINE, MPI_COMM_WORLD );
      
      FREEP( commandList, nCmds );

      for( i = 3; i < argc; i += 2 )
	{
	  if( strcmp( GUI_argv[i] , "-s" ) == 0 )
	    M->IP->lSeedVal = atol(argv[i + 1]);
	}

      if(! M->IP->lSeedVal )  
	M->IP->lSeedVal = (long)time(NULL);
      NEW(M->IP->Datafiles, 1, files);
    }
  else
    {
      NEW( cmdBuffer, CMD_SIZE, char );
      Gibbs_MPI_Recv( M, cmdBuffer, CMD_SIZE, MPI_CHAR, 0, G_MPI_COMMAND_LINE, MPI_COMM_WORLD, &status );

      GUI_argc = ParseCommandLine( cmdBuffer, GUI_argv );
      free( cmdBuffer );

      stripargs(GUI_argc, GUI_argv, M, FALSE); 

      if( M->IP->is_defined[cl_X] )
	p_error( "Hierarchical model: -X cannot be used with hierarchical models" );
    }
	       
  M->IP->is_defined[cl_hm] = TRUE;

  NEW( M->IP->Datafiles->mpiTempFileName, FILENAME_MAX, char);
  sprintf( M->IP->Datafiles->mpiTempFileName, "temp.process.%d.%d.out", my_rank, getpid() );
  remove( M->IP->Datafiles->mpiTempFileName );
#ifndef _RPI_
  if( ! M->IP->is_defined[cl_Z] )
    {
      MPE_IO_Stdout_to_file( M->IP->Datafiles->mpiTempFileName, O_TRUNC | O_RDWR | O_CREAT ); 
      chmod( M->IP->Datafiles->mpiTempFileName, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH );   
      M->IP->Datafiles->mpiTemp_fpt = stdout;
    }
  else
    {
      M->IP->Datafiles->mpiTemp_fpt = fopen( M->IP->Datafiles->mpiTempFileName, "w" );
    }
#else
  M->IP->Datafiles->mpiTemp_fpt = fopen(M->IP->Datafiles->mpiTempFileName, "w" );
#endif
  
  PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, "Processor %d id is %s\n", my_rank, name);  

  M->IP->is_defined[cl_u] = FALSE;    /* override suboptimal output */
  M->IP->is_defined[cl_t] = FALSE;    /* 12/15/99 */

  if( my_rank == 0 )
    {
      M->IP->lSeedVal0 = M->IP->lSeedVal;
      for( source = 0; source < nCmds; source++ )
	{	      
	  nRandSeed = M->IP->lSeedVal + 1009 * source;
	  lRandInfo[0] = nRandSeed;  /* BT 02/25/03 */
	  lRandInfo[1] = M->IP->lSeedVal;
	  Gibbs_MPI_Send( M, lRandInfo, 2, MPI_LONG, source + 1, G_MPI_RAND_SEED, MPI_COMM_WORLD );
	}
      PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, "seed = %d\n", M->IP->lSeedVal );

      NEW( kBuffer, K_BUFFER_SIZE, double );
      NEW( done, nCmds + 1, int );

      /* get initial seed */
      for( i = 0; i < nCmds; i++ )
	{
	  Gibbs_MPI_Recv( M, &nSeeds, 1, MPI_INT, MPI_ANY_SOURCE, G_MPI_SEED, MPI_COMM_WORLD,
			  &status );
	  if( i > 0 && nSeeds != nSeedsPrev )
	    p_error( "Hierarchical model: number of seeds must be the same for all processes" );
	  nSeedsPrev = nSeeds;
	}

      /* get initial counts                                */
      /* we're assuming everyone has the same priors for k */
      
      for( i = 0; i < nCmds; i++ )
	{
	 Gibbs_MPI_Recv( M, kBuffer, K_BUFFER_SIZE, MPI_DOUBLE, MPI_ANY_SOURCE, G_MPI_DATA, MPI_COMM_WORLD,
			 &status );
	  PrintReceivedData( M, kBuffer, &status, 0, 0 );
	  
	  kmax = kBuffer[0];
	  if( i > 0 && kmax != kmaxPrev )
	    p_error( "Hierarchical model: kmax must be the same for all processes" );
	  kmaxPrev = kmax;

	  for( k = 0; k <= kmax; k++ )
	    {
	      if( kBuffer[k+1] < 0 )
		p_error( "Hierarchical model: counts per sequence must >= 0." );	      
	      if( kBuffer[1+kmax+k+1] <= 0 )
		p_error( "Hierarchical model: pseudocounts per sequence must > 0." );	      
	    }
	}
            
      /*=================================================================*/
      /* max likelihood method of estimatinng pseudoCounts from          */
      /* http://research.microsoft.com/~minka/papers/dirichlet/          */
      /*=================================================================*/

      totalProcessors = nCmds;
      while( totalProcessors > 0 )
	{
	  nBusyProcessors = totalProcessors;
	  while( nBusyProcessors > 0 )
	    {
	      NEW( kCnts, kmax+1, double );
	      NEW( kPseudo, kmax+1, double );
	      NEW( kProb, kmax+1, double );
	      
	      kReceived = 0;
	      while( kReceived < nBusyProcessors )
		{
		  Gibbs_MPI_Recv( M, kBuffer, K_BUFFER_SIZE, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
				  &status );
		  source = status.MPI_SOURCE;

		  PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, 
				"%d subopt has finished -- %d processors left\n",
				source, totalProcessors, nBusyProcessors );

		  if( status.MPI_TAG == G_MPI_DATA )
		    {
		      PrintReceivedData( M, kBuffer, &status, 0, 0 );
		      
		      for( k = 0; k <= kmax; k++ )
			{
			  if( kBuffer[k+1] < 0 )
			    p_error( "Hierarchical model: sites per sequence must be >= 0." );
			  if( kBuffer[kmax+1+k+1] <= 0 )
			    p_error( "Hierarchical model: pseudocounts per sequence must be > 0." );
			  kCnts[k] += kBuffer[k+1];
			  kPseudo[k] += kBuffer[kmax+1+k+1];
			}
		      kReceived++;
		    }
		  else if( status.MPI_TAG == G_MPI_SUBOPT_DONE )
		    {
		      totalProcessors--;
		      nBusyProcessors--;
		      done[source] = TRUE;
		      PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, 
				    "%d subopt has finished -- %d processors left\n",
				    source, totalProcessors );
		      if( totalProcessors == 0 )
			break;
		    }
		  else if( status.MPI_TAG == G_MPI_SEED_DONE )
		    {
		      nBusyProcessors--;
		      done[source] = TRUE;
		      PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, "%d subopt seed %d has finished\n",
				    source, (int) kBuffer[0] );
		    }
		}
	      
	      if( kReceived )
		{
		  CalcProbK( M, kmax, kReceived, kCnts, kPseudo, kProb );
		  SendProb( M, kmax, kProb, done, nCmds );
		}
	      else if( totalProcessors ) /* all seeds not done tell everyone to do another */
		{
		  for( source = 1; source <= nCmds; source++ )
		    done[source] = FALSE;
		  SendMessageToAllChildren( M, dTempInfo, 4, G_MPI_CONTINUE, "continue", nCmds );
		}
	      else
		break;
	    }
	  free( kCnts );
	  free( kPseudo );
	  free( kProb );
	}

      free( kBuffer );
      free( done );
      
      /* tell everyone to finish up */
      SendMessageToAllChildren( M, dTempInfo, 4, G_MPI_FINISH, "finish", nCmds );

      /* wait until they are done */
      for( source = 1; source <= nCmds; source++ )
	{
	  PrintMPIMsg( M, "waiting -- %d\n", source );
	  Gibbs_MPI_Recv( M, dTempInfo, 4, MPI_DOUBLE, MPI_ANY_SOURCE, G_MPI_ALL_DONE, MPI_COMM_WORLD,
			  &status );
	  PrintMPIMsg( M, "%d is all done.\n", status.MPI_SOURCE );
	}
      
     /* tell everyone to quit */
      SendMessageToAllChildren( M, dTempInfo, 4, G_MPI_QUIT, "quit", nCmds );

      PrintMPIMsg( M, "Quiting...\n" );
      EndTime(M, &S, M->IP->Datafiles->mpiTemp_fpt);
      
      fclose( M->IP->Datafiles->mpiTemp_fpt );
      FreeData(M);
    }
  else
    {
      Gibbs_MPI_Recv( M, lRandInfo, 2, MPI_LONG, 0, G_MPI_RAND_SEED, MPI_COMM_WORLD,
		      &status );
      M->IP->lSeedVal = lRandInfo[0];
      M->IP->lSeedVal0 = lRandInfo[1];	  
      PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, "seed = %d\n", M->IP->lSeedVal );
      M->IP->is_defined[cl_s] = TRUE;

      get_inputs(M);
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

      /* send master the number of seeds  and counts*/
      Gibbs_MPI_Send( M, &(M->IP->nSeeds), 1, MPI_INT, 0, G_MPI_SEED, MPI_COMM_WORLD );
      SendCounts( M , 0, 0);
	      
      /* allocate space for alignment */
      /*  Pos[MotifType][SeqLen] */
      NEWP(Pos,M->IP->nNumMotifTypes,PoSition);
      for(t=0;t<M->IP->nNumMotifTypes;t++)
	NEW(Pos[t], M->IP->nSeqLen, PoSition);

      nNumMotifs = copy_motif_num(M->IP);   /* Save the input motif count */   /* BT 06/22/2000 */
     
      fpt = M->IP->Datafiles->out_fpt;
      init_maxdata(&maxData);
      init_maxdata(&optmax);
    
      maxData = suboptSampler(M, Pos, &S);    /* Find Suboptimal Solution   */
      reset_motif_num(nNumMotifs, M->IP);   /* restore the motif count so it appears corrctly in output */

      printargs(GUI_argc, GUI_argv, M->IP, M->IP->Datafiles->mpiTemp_fpt );
      print_options( M );
      PrintSeqDescriptions( M ); 

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
      
      optmax = NearOptSampler(M, maxData, Pos);   /* Find Near optimal Solution */
      
      if(!M->IP->is_defined[cl_m] )
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
      
      free_maxdata(&maxData, M->IP);
      free_maxdata(&optmax, M->IP);
      
      FREEP(Pos,M->IP->nNumMotifTypes);  
      free(Pos);

      times( &tBuffer2);     /* 12/15/99 */
      PrintTempOut( fpt, "Elapsed time: %f secs\n", 
		    (double) (tBuffer2.tms_utime - tBuffer1.tms_utime) / M->IP->ticks );	

      /* send master a message saying we're done */
      PrintMPIMsg( M, "Sending all done message \n");
      Gibbs_MPI_Send( M, dTempInfo, 4, MPI_DOUBLE, 0, G_MPI_ALL_DONE, MPI_COMM_WORLD );
      PrintMPIMsg( M, "Waiting for quit message \n");
      Gibbs_MPI_Recv( M, dTempInfo, 4, MPI_DOUBLE, 0, G_MPI_QUIT, MPI_COMM_WORLD,
		      &status );
      PrintMPIMsg( M, "Received quit message \n");
      PrintMPIMsg( M, "Sleeping...\n" );
      sleep( 1 );
      PrintMPIMsg( M, "Done sleeping\n" ); 

      fclose( M->IP->Datafiles->mpiTemp_fpt );
      FreeData(M); 
    }
}


char** ReadCommandFile( char *commandFile, int *nCmds )
{
  FILE *fp;
  int  nAlloc;
  char **commandList;
  char *buffer;

  fp = fopen( commandFile, "r" );
  if( fp == NULL )
    {
      NEW( buffer,  FILENAME_MAX + 64, char ); 
      sprintf( buffer, "Command file %s could not be opened.", commandFile );      
      p_error( buffer );
    }

  NEWP( commandList, CMD_ALLOC_SIZE, char );
  nAlloc = CMD_ALLOC_SIZE;

  NEW( buffer, CMD_SIZE, char );
 
  *nCmds = 0;
  while( GetALine( fp, buffer, CMD_SIZE, TRUE ) )  
    {
      NEW( commandList[*nCmds], strlen( buffer ) + 1, char );
      strcpy( commandList[*nCmds], buffer ); 
      
     (*nCmds)++;
      if( (*nCmds % CMD_ALLOC_SIZE) == 0 )
	{
	  commandList = (char **) realloc( commandList, (nAlloc + CMD_ALLOC_SIZE) * sizeof(char *) );
	  nAlloc += CMD_ALLOC_SIZE;
	}
    }
  fclose( fp );

  commandList = (char **) realloc( commandList, (*nCmds) * sizeof(char *) );

  return commandList;
}


int ParseCommandLine( char *commandLine, char **argv )
{
  char token[TOKEN_SIZE];
  int  argc = 1;
  char *ptr;

  ptr = commandLine;
  while( (ptr = GetAToken( ptr, token )) )
    {
      NEW( argv[argc], strlen(token ) + 1, char );
      strcpy( argv[argc], token );
      argc++;
    }
  
  return argc;
}


void ExchangeK( Model B, int seed, int iter )
{
  int         k;
  RPType      RP;
  MPI_Status  status;
  double      *kBuffer;

  RP = B->RP;

  NEW( kBuffer, RP->nMaxBlocks + 1, double );
    
  SendCounts( B, seed, iter );

  Gibbs_MPI_Recv( B, kBuffer, RP->nMaxBlocks  + 1, 
		  MPI_DOUBLE, 0, G_MPI_DATA, MPI_COMM_WORLD, &status ); 
  PrintReceivedData( B, kBuffer, &status, seed, iter );

  for( k = 0; k <= RP->nMaxBlocks; k++ )
    RP->probK[k] = kBuffer[k];
  RP->useProbK = TRUE; 

  free( kBuffer );
}


void PrintReceivedData( Model M, double *buffer, MPI_Status *status, int seed, int iter )
{
  int         received;
  int         i;
  struct tm   *currTime;
  time_t      tloc;

  MPI_Get_count( status, MPI_DOUBLE, &received );
  time( &tloc );
  currTime = localtime( &tloc );
  PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, "%d has sent %d numbers at %d:%d:%d seed = %d iter = %d\n", 
		 status->MPI_SOURCE, received,
		 currTime->tm_hour, currTime->tm_min, currTime->tm_sec,
		 seed, iter );
  for( i = 0; i < received; i++ )
    PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, "%g ", buffer[i] );
  PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, "\n" );
}


void PrintMPIMsg( Model B, char *fmt, ...)
{
  struct tm   *currTime;
  time_t      tloc;
  va_list     argp;
  
  time( &tloc );
  currTime = localtime( &tloc );
  PrintTempOut( B->IP->Datafiles->mpiTemp_fpt, 
		"%d:%d:%d --", 
		currTime->tm_hour, currTime->tm_min, currTime->tm_sec );

  va_start(argp, fmt);
  vfprintf(B->IP->Datafiles->mpiTemp_fpt, fmt, argp);
  va_end(argp);
  fflush( B->IP->Datafiles->mpiTemp_fpt ); 
}


void SendMessageToAllChildren( Model M, double *dTempInfo, int dataLength, int tag, char *tagName, int nProcesses )
{
  int source;

  for( source = 1; source <= nProcesses; source++ )
    {
      PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, "sending %s to %d\n",
		    tagName, source );
     Gibbs_MPI_Send( M, dTempInfo, dataLength, MPI_DOUBLE, source, tag, MPI_COMM_WORLD );
    }
}


/* mimimize the P(k|pseudo) over kPseudo                                               */
/* Uses iterative method from http://research.microsoft.com/~minka/papers/dirichlet/   */

int MinimizeK( Model M, double *kPseudo, double **kCnts, int kmax, int kReceived )
{
  double *prevPseudo;
  double maxPrev;
  int    i;
  int    k;
  double pseudoSum;
  double num;
  double denom;
  double ksum;
  double diff;
  double iter = 0;
  int    code = 0;

  NEW( prevPseudo, kmax + 1, double );

  do
    {
      for( maxPrev = 1.0, k = 0; k <= kmax; k++ )
	{
	  prevPseudo[k] = kPseudo[k];
	  maxPrev = max( maxPrev, prevPseudo[k] );
	}
		      
      for( pseudoSum = 0.0, k = 0; k <= kmax; k++ )
	pseudoSum += kPseudo[k];		    
		      
      for( denom = 0.0, i = 0; i < kReceived; i++ )
	{
	  for( ksum = 0, k = 0; k <= kmax; k++ )
	    ksum += kCnts[i][k];
	  denom += digamma( ksum + pseudoSum );
	}
      denom -= kReceived * digamma( pseudoSum );
		      
      for( diff = -1.0, k = 0; k <= kmax; k++ )
	{
	  for( ksum = 0, i = 0; i < kReceived; i++ )
	    ksum += digamma( kCnts[i][k] + kPseudo[k] );
	  num = ksum - kReceived * digamma( kPseudo[k] );;
	  kPseudo[k] *= (num / denom );
	  kPseudo[k] = max( MIN_PSEUDO, kPseudo[k] );
	  diff = max( diff, fabs( kPseudo[k] - prevPseudo[k] ) );
	  PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, 
			"k: %d pseudo: %g delta %g num %g denom %g n/d %g\n", 
			k, kPseudo[k], kPseudo[k]-prevPseudo[k],
			num, denom, num/denom);
	}
		      
      diff /= maxPrev;
      iter++;
      PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, "iter: %d diff: %g\n", 
		    iter, diff );

    } while( iter < MAX_ITER && diff > TOL );

  free( prevPseudo );

  if( iter >= MAX_ITER )
    code = 1;

  return code;
}


/* this version assumes that all k are drawn from same distribution */

int MinimizeK2( Model M, double *kPseudo, double **kCnts, int kmax, int kReceived )
{
  double *prevPseudo;
  double *kCntsSum;
  double maxPrev;
  int    i;
  int    k;
  double pseudoSum;
  double num;
  double denom;
  double ksum;
  double diff;
  int    iter = 0;
  int    code = 0;

  NEW( prevPseudo, kmax + 1, double );
  NEW( kCntsSum, kmax + 1, double );

  for( k = 0; k <= kmax; k++ )
    {
      for( i = 0; i < kReceived; i++ )
	{
	  kCntsSum[k] += kCnts[i][k];
	}
    }

  do
    {
      for( maxPrev = 1.0, k = 0; k <= kmax; k++ )
	{
	  prevPseudo[k] = kPseudo[k];
	  maxPrev = max( maxPrev, prevPseudo[k] );
	}
		      
      for( ksum = 0, pseudoSum = 0.0, k = 0; k <= kmax; k++ )
	{
	  ksum += kCntsSum[k];
	  pseudoSum += kPseudo[k];
	}

      denom = digamma( ksum + pseudoSum ) - digamma( pseudoSum );
		      
      for( diff = -1.0, k = 0; k <= kmax; k++ )
	{
	  num = digamma( kCntsSum[k] + kPseudo[k] ) - digamma( kPseudo[k] );
	  kPseudo[k] *= (num / denom );
	  kPseudo[k] = max( MIN_PSEUDO, kPseudo[k] );
	  diff = max( diff, fabs( kPseudo[k] - prevPseudo[k] ) );
	  PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, 
			"k: %d pseudo: %g delta %g num %g denom %g n/d %g\n", 
			k, kPseudo[k], kPseudo[k]-prevPseudo[k],
			num, denom, num/denom);
	}
		      
      diff /= maxPrev;
      iter++;
      PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, "iter: %d diff: %g\n", 
		    iter, diff );

    } while( iter < MAX_ITER && diff > TOL );

  free( prevPseudo );
  free( kCntsSum );

  if( iter >= MAX_ITER )
    code = 1;

  return code;
}


void SendCounts( Model M, int seed, int iter )
 {
   double      *kBuffer;
   double      *cntSitesPerSeq;
   int         k;
   int         i;
   struct tm   *currTime;
   time_t      tloc;
   
   /* send master the counts & pseudocounts */
   NEW( kBuffer, 2 * (M->RP->nMaxBlocks + 1) + 1, double );
   NEW( cntSitesPerSeq,  M->RP->nMaxBlocks + 1, double );
   
   for( i = 0; i < M->IP->nNumSequences; i++ )
     cntSitesPerSeq[M->RP->nSites[i]] += GetSeqWeight( M, i, 0 );
   
   kBuffer[0] = M->RP->nMaxBlocks;
   for( k = 0; k <= M->RP->nMaxBlocks; k++ )
     {	  
       kBuffer[k+1] = cntSitesPerSeq[k];
       kBuffer[1 + M->RP->nMaxBlocks + k + 1] = M->RP->priorSitesPerSeq[k];
     }
   
   time( &tloc );
   currTime = localtime( &tloc );
   PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, 
		 "\n%d is sending %d numbers at %d:%d:%d seed = %d iter = %d\n", 
		 M->IP->nRank, 2 * (M->RP->nMaxBlocks + 1) + 1,
		 currTime->tm_hour, currTime->tm_min, currTime->tm_sec, seed, iter );
   for( k = 0; k < 2 * (M->RP->nMaxBlocks + 1) + 1; k++ )
     PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, "%g ", kBuffer[k] );
   PrintTempOut( M->IP->Datafiles->mpiTemp_fpt, "\n" );
   
   Gibbs_MPI_Send( M, kBuffer, 2 * (M->RP->nMaxBlocks + 1) + 1, MPI_DOUBLE, 0, G_MPI_DATA, MPI_COMM_WORLD );
   
   free( cntSitesPerSeq );
   free( kBuffer );
 }


void CalcProbK( Model B, int kmax, int kReceived, double *kCnts, double *kPseudo, double *kProb )
 {
   int    k;
   double denom = 0;

   for( k = 0; k <= kmax; k++ )
     kPseudo[k] /= kReceived;

   for( k = 0; k <= kmax; k++ )
     {
       denom += (kCnts[k] + kPseudo[k]);
     }

   for( k = 0; k <= kmax; k++ )
     kProb[k] = (kCnts[k] + kPseudo[k])/ denom;
 }


void SendProb( Model B, int kmax, double *kProb, int *done, int nProcesses )
{
  int         source;
  struct tm   *currTime;
  time_t      tloc;
  int         k;

  for( source = 1; source <= nProcesses; source++ )
    {
      if( ! done[source] )
	{
	  time( &tloc );
	  currTime = localtime( &tloc );
	  PrintTempOut( B->IP->Datafiles->mpiTemp_fpt, 
			"sending %d probs to %d at %d:%d:%d\n", 
			kmax + 1, source,
			currTime->tm_hour, currTime->tm_min, currTime->tm_sec );
	  for( k = 0; k <= kmax; k++ )
	    PrintTempOut( B->IP->Datafiles->mpiTemp_fpt, "%g ", kProb[k] );
	  PrintTempOut( B->IP->Datafiles->mpiTemp_fpt, "\n" );
	  
	  Gibbs_MPI_Send( B, kProb, kmax + 1, MPI_DOUBLE, source, G_MPI_DATA, MPI_COMM_WORLD );
	}
    }
}
