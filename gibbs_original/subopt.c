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
/*************************************************************************/
/* $Id: subopt.c,v 1.19 2009/04/23 18:43:55 Bill Exp $           */  
/*                                                                       */
/* Author : Eric C. Rouchka   July 17, 1996                              */
/*          Jun Zhu, October 20, 1996                                    */
/*	    Bill Thompson     2/7/97					 */ 
/*                                                                       */
/*************************************************************************/

#include "subopt.h"

void GetWidthCounts( Model B );
#ifdef _MPI_
void SendSuboptMsg( Model B, int tag, int seed );
#endif

/*===================   suboptSampler   =================================*/
/*                                                                       */
/* DESCRIPTION  : Finds a suboptimal motif alignment that will be better */
/*=======================================================================*/

MaxResults suboptSampler(Model B, PoSition **Pos, SimTime *S)
 
{
   int        **nNumMotifs, 
              **startPos,   /* startPos[MotifType][motifnum] */
              i, j, seed_run = 0, t, k; 
   int        maxnum;  /* maximum number of motif in one type motif */
   register   Mlist M; 
   MaxResults maxData, locMax;
   RPType     RP;
   IPtype     IP;      /* for quick access of input data    */
   PoSition   *Pos_t;  /* Pos_t=Pos[t] */
   int	      n;				/* BT 2/7/97 */ 
   int        nPos;
   int        sum;
   double     dCurrProb;
#ifdef _MPI_
   MPI_Status status;
   double     dTempInfo[4];
#endif

   IP=B->IP;
   RP = B->RP;

   init_maxdata(&maxData);
   init_maxdata(&locMax);
   BeginTime(S);

#ifdef _MPI_
   if( ! IP->is_defined[cl_hm] )
     {
       /* get the go ahead to start */
       Gibbs_MPI_Recv( B, &dTempInfo, 4, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, 
		       &status );
       if( status.MPI_TAG == G_MPI_DONE )
	 return maxData;
       
       if( IP->is_defined[cl_X] && B->AN->bExchange )
	 {
	   B->AN->currTemp = dTempInfo[0];
	   B->AN->dMinTemp = B->AN->currTemp;
	   B->AN->dMaxTemp = B->AN->currTemp;
	   PrintTempOut(IP->Datafiles->mpiTemp_fpt, "currTemp = %f maxTemp = %f minTemp = %f\n",
			B->AN->currTemp, B->AN->dMaxTemp, B->AN->dMinTemp );
	 }
     }
#endif
   
   /* save the inital number of motifs */
   for(t = 0; t < IP->nNumMotifTypes; t++)		/* BT 2/7/97 */
     {
       RP->nPriorMotifSites[t][FORWARD] = IP->nNumMotifs[t][FORWARD];
       RP->nPriorMotifSites[t][REVERSE] = IP->nNumMotifs[t][REVERSE];
     }

   /* copy the original counts into current counts */
   copy_counts(B);  
   nNumMotifs = copy_motif_num(IP);
   zero_motifs(IP); /* set motif count to zero */

   /* If not running Wilcox test initialize random number gerator. */
   /* Wilcox already did this */

   if( ! IP->is_defined[cl_l] )      /* BT 7/23/97 */ 
       sRandom(B, IP->lSeedVal);  

   /* calculate the probability of NULL model */
   /* IP->dnull_map = CalcMapProb(B, IP->is_defined[cl_R]); */
   IP->dnull_map = CalcNullMap(B);  /* 2/16/2001 */

   reset_motif_num(nNumMotifs,IP);

   /* mark all possible motif sites */
   if( ! IP->inCentroidAlign )
     set_indicator_vector(Pos, B->Seq, B);
   maxnum = findMaxNumMotif(IP);
   
   if( IP->is_defined[cl_E] ||
       RP->nUseSpacingProb || IP->is_defined[cl_T]  )  /* BT 11/20/2000 */
     {
       CountAlignments( B, Pos ); 
       SaveAlignmentCounts( B );
     }
       
   /* allocate space for start positions of motifs */
   NEWP(startPos,IP->nNumMotifTypes, int);
   for(i = 0; i <IP->nNumMotifTypes ; i++)
     NEW(startPos[i],maxnum, int);

   /* Loop through for a number of different seeds */
   while(seed_run++ < IP->nSeeds) 
     {  
       if( IP->is_defined[cl_X] )
	 {
	   /* if( ! B->AN->bExchange )
	      B->AN->currTemp = B->AN->dMaxTemp; */
	   if( ! IP->is_defined[cl_Z] )	    
	     printf( "Current temperature: %7.2f\n", B->AN->currTemp );
	 }
       IP->nSeedRun = seed_run;
        
       if(!IP->is_defined[cl_F])
	 initMask(B);

       if( IP->is_defined[cl_X] )
	 init_maxdata( &(B->AN->results[seed_run]) );
       
       /* initialize the Pos */
       for(t=0;t<IP->nNumMotifTypes;t++)
	 {
	   Pos_t=Pos[t];
	   for(i = 0; i < IP->nSeqLen; i++) 
	     {
	       Pos_t[i].nInMotif=FALSE;
	       Pos_t[i].nMotifStartPos=FALSE;
	     }
	 }

       for( t= 0; t < IP->nNumMotifTypes; t++ )          /* BT 5/23/97 */
	 {
	   IP->nMotifLen[t] = IP->nInputMotifLen[t];        /* restore original lengths */
	   SetPossibleSites( B, t );                       /* restore original site count */
	 }
       if( IP->is_defined[cl_d] && ! IP->inCentroidAlign ) 
	 set_indicator_vector(Pos, B->Seq, B);
       
       if( ! IP->site_samp )                      /* BT 7/16/97 */
	 set_posterior_prob(IP, B->C);         /* Set initial posterior prob */
       
       /* random select a set of motifs */
       if( ! IP->is_defined[cl_V] ) /* BT 04/16/03 */
	 {
	   if( B->InitPos == NULL || seed_run > 1 )
	     {
	       if( IP->is_defined[cl_A] )
		 {
		   set_posterior_prob(IP, B->C);     
		   for( t = 0; t < IP->nNumMotifTypes; t++ )
		     {
		       IP->nNumMotifs[t][FORWARD] = 0;
		       IP->nNumMotifs[t][REVERSE] = 0;
		     }
		 }
	       set_random_sequences(B, Pos, startPos); /* set alignment &*/
	     }
	   else
	     {
	       SetInitSequences( B, Pos, startPos );  /* BT 8/5/98 */
	       
	       set_posterior_prob(IP, B->First);     

	       if( RP->bUsePosMatrix )
		 InitializePosMatrix( B );
	       else if( RP->bUseTrans )
		 InitializeTransMatrix( B );
	     }
	 }
       else
	 {
	   if( B->InitPos == NULL )
	     {
	       for( t = 0; t < IP->nNumMotifTypes; t++ )
		 {
		   IP->nNumMotifs[t][FORWARD] = 0;
		   IP->nNumMotifs[t][REVERSE] = 0;
		 }
	     }
	   else
	     {
	       SetInitSequences( B, Pos, startPos ); 
	       set_counts( B );  /* BT 10/23/03 */
	       set_posterior_prob(IP, B->First);     
	     }
	   set_posterior_prob(IP, B->C);     
	 }

       reset_counts(B, startPos, Pos);             /* its counts     */

       if(  B->InitPos != NULL && seed_run == 1 )
	 {
	   RP->dInitProb = CalcMapProb( B, IP->is_defined[cl_R] );
	   if( ! IP->inCentroidAlign )
	     {
	       fprintf( IP->Datafiles->out_fpt, "seed = %d Initial MAP =  %.5f sites = %d \n",  
			IP->nSeedRun, 
			RP->dInitProb,
			TotalNumMotifs( B ));
	       for( t = 0; t < IP->nNumMotifTypes; t++ )
		 {
		   fprintf(IP->Datafiles->out_fpt, 
			    "Motif %d Map =  %.5f Frag =  %.5f\n", t, 
			    CalcMotifMap(B, t, IP->is_defined[cl_R]),
			    CalcMotifFragMap(B, t, IP->is_defined[cl_R]) );
		 }
	     }


#ifdef _MPI_
	   PrintTempOut( IP->Datafiles->mpiTemp_fpt,
			 "rank = %d process = %d seed = %d Initial MAP =  %.5f sites = %d\n",  
			 IP->nRank,
			 IP->nMPIProcesses,
			 IP->nSeedRun, 
			  RP->dInitProb,
			 TotalNumMotifs( B ) );

	   for( t = 0; t < IP->nNumMotifTypes; t++ )
	     {
	       PrintTempOut( IP->Datafiles->mpiTemp_fpt,
			     "%d Motif =  %.5f Frag =  %.5f\n", t, 
			     CalcMotifMap(B, t, IP->is_defined[cl_R]),
			     CalcMotifFragMap(B, t, IP->is_defined[cl_R]) );
	     }
	   PrintTempOut( IP->Datafiles->mpiTemp_fpt,
			 "Bkgnd =  %.5f Beta =  %.5f Null = %.5f\n",  
			 CalcBkgndMap(B, IP->is_defined[cl_R]),
			 CalcBetaMap( B, IP->is_defined[cl_R]),
			 IP->dnull_map );	  
#else
	   fprintf( stdout, "seed = %d Initial MAP =  %.5f sites = %d\n",  
		    IP->nSeedRun, 
		    B->RP->dInitProb,
		    TotalNumMotifs( B ) );

	   for( t = 0; t < IP->nNumMotifTypes; t++ )
	     {
	       fprintf( stdout, 
			"%d Motif =  %.5f Frag =  %.5f\n", t, 
			CalcMotifMap(B, t, IP->is_defined[cl_R]),
			CalcMotifFragMap(B, t, IP->is_defined[cl_R]) );
	     }
	   fprintf( stdout, 
		    "Bkgnd =  %.5f Beta =  %.5f Seq = %.5f Null = %.5f\n",  
		    CalcBkgndMap(B, IP->is_defined[cl_R]),
		    CalcBetaMap( B, IP->is_defined[cl_R]),
		    CalcSitePerSeqMap( B ), 
		    IP->dnull_map );	  
#endif

	   for( t = 0; t < IP->nNumMotifTypes; t++ )
	     {
	       fprintf( stdout, 
			"-----------------------------------------------------------\n" );
	       fprintf( stdout, "                          MOTIF %c\n\n", (char)(97 + t));
	       DumpMotifPositions( t, B, Pos, stdout );
	       fprintf( stdout, "%d sites\n", NUMMOTIFS( IP->nNumMotifs[t]) );
	     }
	   fflush( stdout );
	 }
       
       /* create a motif element list */
       M = set_motif_info(IP, startPos, B->Seq);

#ifdef _DEBUG_
       /* ==================================================================== */
       for( t = 0; t < IP->nNumMotifTypes; t++ )
	 {
	   fprintf( stdout, 
		    "-----------------------------------------------------------\n" );
	   fprintf( stdout, "                          MOTIF %c\n\n", (char)(97 + t));
	   DumpMotifPositions( t, B, Pos, stdout );
	   fprintf( stdout, "%d sites\n", NUMMOTIFS( IP->nNumMotifs[t]) );
	   DumpCounts( B, t, stdout );
	   fflush( stdout );
	 }
       /*  ==================================================================== */
#endif

       if( ! IP->is_defined[cl_Z] )
	 {
	   put_prior(B);
	   fprintf(stdout, "\r** %d **\n", seed_run);
	   fflush( stdout );    /* BT 9/19/97 */
	 }

       if(seed_run > 1) 
	 {
	   RestoreAlignmentCounts( B );
	   if(IP->site_samp) 
	     locMax = site_sampler(B, Pos, M); 
	   else if( IP->is_defined[cl_bayes] )
	     locMax = bayes_sampler(B, Pos, M, seed_run); 
	   else if( IP->is_defined[cl_E] )
	     locMax = rsite_sampler(B, Pos, M, seed_run); 
	   else
	     locMax = motif_sampler(B,Pos,M);

	   if( IP->is_defined[cl_X] )
	     CopyMaxResults( &(B->AN->results[seed_run]), &locMax, B );
	   
	   dCurrProb = locMax.dProbability;

	   if(locMax.dProbability > maxData.dProbability) 
	     {
	       /* Check to see if current run is max */
	       free_maxdata(&maxData, IP);                  
	       maxData = locMax;     
	       maxData.nSuboptSeed = seed_run;
	       if( (! IP->is_defined[cl_Z]) && IP->is_defined[cl_opt] )
		 print_maxData(IP->nNumMotifTypes, maxData);
	       }
	   else 
	     free_maxdata(&locMax, IP); 
	 }
       else 
	 {  /* Set Maximum first time through */
	     if(IP->site_samp) 
	       maxData = site_sampler(B, Pos, M); 
	     else if( IP->is_defined[cl_bayes] )
	       maxData = bayes_sampler(B, Pos, M, seed_run); 
	     else if( IP->is_defined[cl_E] )
	       maxData = rsite_sampler(B, Pos, M, seed_run); 
	     else              
	       maxData = motif_sampler(B,Pos,M);
	     
	     dCurrProb = maxData.dProbability;
	     maxData.nSuboptSeed = seed_run;

	     if( IP->is_defined[cl_X] )
	       CopyMaxResults( &(B->AN->results[seed_run]), &maxData, B );
	     if( (! IP->is_defined[cl_Z]) && IP->is_defined[cl_opt] )
	       print_maxData(IP->nNumMotifTypes, maxData);
	 }
              
       reset_motif_num(nNumMotifs, IP);
       copy_counts(B);                        /* Counts w/o motifs */
       free_motifs(B, M); 

#ifdef _MPI_
       if( IP->is_defined[cl_hm] )
	 {
	   if( seed_run == IP->nSeeds )
	     SendSuboptMsg( B, G_MPI_SUBOPT_DONE, seed_run );
	   else
	     SendSuboptMsg( B, G_MPI_SEED_DONE, seed_run );
	 }
       else
	 {
	   dTempInfo[1] = dCurrProb;
	   dTempInfo[3] = 0;
	   if( IP->is_defined[cl_opt] )  
	     {
	       for(t = 0; t < IP->nNumMotifTypes; t++) 
		 {
		   dTempInfo[3] += (double) maxData.nNumMotifs[t];
		 }
	     }
	   
	   Gibbs_MPI_Send( B, dTempInfo, 4, MPI_DOUBLE, 0, G_MPI_DATA, MPI_COMM_WORLD );
	 }

       Gibbs_MPI_Recv( B, dTempInfo, 4, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, 
		       &status );

       if( status.MPI_TAG == G_MPI_DONE )
	 {
	   PrintTempOut( IP->Datafiles->mpiTemp_fpt, "Received MPI_DONE signal\n" );
	   break;
	 }
       else if( status.MPI_TAG == G_MPI_FINISH )
	 {
	   PrintTempOut( IP->Datafiles->mpiTemp_fpt, "Received subopt end signal\n" );
	   break;
	 }
#endif       
     }/* end of while loop, at this time we have MaxData */
   FREEP(startPos,IP->nNumMotifTypes );
   
   if(!IP->is_defined[cl_F] && maxData.F) 
     {
       for(t = 0; t < IP->nNumMotifTypes; t++)		/* BT 2/7/97 */
	 { 
	   for(n = 0; n < B->F->nMaxLen[t]; n++)
	     B->F->nColMask[t][n] = maxData.F->nColMask[t][n];
	   B->F->FragWidth[t] = maxData.F->FragWidth[t];
	   for( n = 0; n < IP->nMotifLen[t]; n++ )
	     B->F->fragPos[t][n] = maxData.F->fragPos[t][n];
	 } 
     }
   
   for(t = 0; t < IP->nNumMotifTypes; t++)
     {            /* Add in motifs   */
       if( maxData.nMotifLen )
	 IP->nMotifLen[t] = maxData.nMotifLen[t];
       /* Reset Pos in case width changed */
       if( ! IP->inCentroidAlign )
	 set_indicator_vector(Pos, B->Seq, B);       /* BT 9/12/97 */  
       if(  maxData.nNumMotifs )
	 {
	   for(i = 0; i < maxData.nNumMotifs[t]; i++)        /* from the maximum */
	     adjust_counts(B, ADD, maxData.nMotifLoc[i][t], /* alignment       */
			   t, maxData.RevComp[i][t]);
	 }
     }
   
   if(  maxData.nNumMotifs )
     {
       setMotifNum(IP, maxData);      /* BT 5/30/97 */
     }
   
   if( IP->is_defined[cl_u] )
     print_info(B, maxData, TRUE, SUBOPT);		/* BT 3/19/97 */
   
   if( IP->is_defined[cl_d] && IP->is_defined[cl_q] )
     GetWidthCounts( B );

   FREEP(nNumMotifs, IP->nNumMotifTypes);
   free(nNumMotifs);

   if( IP->is_defined[cl_Q] )     /* BT 3/27/98 */
     {
       fprintf( IP->Datafiles->occur_fpt, "seq pos pos2 motif count\n" );
       for( i = 0; i < IP->nNumSequences; i++ )
	 {
	   for( j = 0; j < SequenceLength( B, i ); j++ )
	     {
	       nPos =  SequenceStartPos( B, i ) + j;
	       for(  t = 0;  t < IP->nNumMotifTypes; t++ )
		 {
		   for( sum = 0, k = 1; k <= IP->nSeeds; k++ )
		     {
		       sum += IP->nAlignCnts[k][t][i][j];
		     }
		   fprintf( IP->Datafiles->occur_fpt, "%5d %5d %5d %5d %5d\n",
			    i, j, nPos, t, sum );
		 }
	     }
	   fflush(  IP->Datafiles->occur_fpt );
	 }
     }

   if( ! IP->is_defined[cl_Z] && ! IP->is_defined[cl_nopt] )
     printf( "Max subopt MAP found on seed %d\n", maxData.nSuboptSeed );
   
   return maxData;
}


void GetWidthCounts( Model B )
{
  IPtype  IP;
  int     t;
  int     i;
  double  dWtCnt[256];
  FILE    *fpt;

  IP=B->IP;
  
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      for( i = 0; i <= IP->nMaxMotifLen[t]; i++ )
	dWtCnt[i] = IP->nWidthCnts[t][0][i] + IP->nWidthCnts[t][1][i];

      fpt = IP->Datafiles->out_fpt;
      if( fpt != NULL )
	{
	  fprintf( fpt, "\n\nWidth Sample Counts for motif %c\n", 'a' + t );
	  fprintf( fpt, "Width     Count\n" ); 
	  for( i =0; i <= IP->nMaxMotifLen[t]; i++ )
	    fprintf( fpt, "%4d  %8.2f\n",
		     i, IP->nWidthCnts[t][0][i] + IP->nWidthCnts[t][1][i] );
	  fprintf( fpt, "\n" );
	}
      t = t; 
    }
}

#ifdef _MPI_
void SendSuboptMsg( Model B, int tag, int seed )
{
  RPType  RP;
  double  *buffer;

  RP = B->RP;
  
  NEW( buffer,  RP->nMaxBlocks + 1, double );  
  buffer[0] = seed;
  Gibbs_MPI_Send( B, buffer, RP->nMaxBlocks + 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD );
  free( buffer );
}
#endif



