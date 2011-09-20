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
/* $Id: sampling.c,v 1.13 2009/04/23 18:43:55 Bill Exp $          */
/*                                                                        */
/* Author :       Eric C. Rouchka,      July 15, 1996                     */
/*                Jun Zhu, October,20, 1996                               */
/* Description :  This file contains the functions which are used to find */
/*                the best motif alignment including :                    */
/*                                                                        */
/*                find_best_sites                                         */
/**************************************************************************/

#include "sampling.h" 

void PrintPos( Model B, PoSition **Pos, int t );


/**********************  find_best_sites **********************************/
/*                                                                        */
/* This function finds the best motif alignment                           */
/**************************************************************************/

void find_best_sites(register Model B)

{
   PoSition   **Pos;
   MaxResults maxData, optmax;
   SimTime    S;
   int        i, t;
   FILE       *fpt;
   struct tms tBuffer1;
   struct tms tBuffer2;
   
   times( &tBuffer1);     /* DEBUG */

   if( ! B->IP->is_defined[cl_Z] )
     printf( "\nSeed = %ld\n\n", B->IP->lSeedVal );

   /* allocate space for alignment */
   /*  Pos[MotifType][SeqLen] */
   NEWP(Pos,B->IP->nNumMotifTypes,PoSition);
   for(t=0;t<B->IP->nNumMotifTypes;t++)
     NEW(Pos[t], B->IP->nSeqLen, PoSition);
     
   fpt = B->IP->Datafiles->out_fpt;
   init_maxdata(&maxData);
   init_maxdata(&optmax);
   fprintf(fpt, "\n\n\n\n\n");
   if( B->IP->is_defined[cl_u])
     { 
       fprintf(fpt,
	       "=====================================================================\n");
       fprintf(fpt, 
	       "=========================SUBOPTIMAL SAMPLING RESULTS=================\n");
       fprintf(fpt, 
	       "=====================================================================\\n\n");
     }	/* BT 3/19/97 */
   maxData = suboptSampler(B, Pos, &S);    /* Find Suboptimal Solution   */
   
   if( ! B->IP->is_defined[cl_V] )
     {
       if( ! B->IP->is_defined[cl_nopt] )
	 {
	   fprintf(fpt, "\n\n\n\n\n======================================================================\n");
	   fprintf(fpt, 
		   "======================== NEAR OPTIMAL RESULTS ========================\n");
	   fprintf(fpt, 
		   "======================================================================\n\n");
	 }
	   
       if( B->IP->is_defined[cl_opt] )
	 optmax = NearOptSampler(B, maxData, Pos);/* Find Near optimal Solution */
       
       if(B->IP->is_defined[cl_opt] && ! B->IP->is_defined[cl_m] )
	 {             /* Maximize the MAP value */
	   fprintf(fpt, "\n\n\n\n\n======================================================================\n");
	   fprintf(fpt, 
		   "======================== MAP MAXIMIZATION RESULTS ====================\n");
	   fprintf(fpt,
		   "======================================================================\n\n");
	   optmax = maximize_map(optmax, B, Pos);   /* BT 2/21/97 */
	 }

       if(B->IP->is_defined[cl_e]) 
	 {    /* Run Expection Maximization */
	   fprintf(fpt, "\n\n\n\n\n=======================================================================\n");
	   fprintf(fpt, 
		   "======================== EXPECTATION MAXIMIZATION RESULTS =============\n");
	   fprintf(fpt, 
		   "=======================================================================\n\n");
	   EM(B, Pos, optmax); 
	 }

       if(  B->IP->is_defined[cl_bayes] )
	 {
	   fprintf(fpt, "\n\n\n\n\n======================================================================\n" );
	   fprintf(fpt, 		   
		   "========================== CENTROID RESULTS ==========================\n" );
	   fprintf(fpt, 		   
		   "======================================================================\n\n" );
	   Centroid( B );
	 }
     }
   EndTime(B, &S, stdout);

   times( &tBuffer2);     /* DEBUG */
   fprintf( fpt, "Elapsed time: %f secs\n", 
	    (double) (tBuffer2.tms_utime - tBuffer1.tms_utime) / B->IP->ticks );

   free_maxdata(&maxData, B->IP);
   free_maxdata(&optmax, B->IP);
   
   FREEP(Pos,B->IP->nNumMotifTypes);   /* BT 8/27/97 */
   free(Pos);
   
   for(t = 0; t < B->IP->nNumMotifTypes; t++) 
     {
       B->IP->DOF[t] = B->IP->nMotifLen[t] * (B->IP->nAlphaLen - 1);
       for(i = 0; i < B->IP->nMotifLen[t]; i++) 
	 {
	   if(B->IP->AltModel->Collapsed[t][i])
	     B->IP->DOF[t] -= (B->IP->nAlphaLen) / 2.0 - 1.0;
	   if((B->IP->AltModel->Palandromic[t][i]) && 
	      (i < (B->IP->nMotifLen[t] / 2 + B->IP->nMotifLen[t] % 2)))
	     B->IP->DOF[t] -= B->IP->nAlphaLen - 1.0;
	   if((B->IP->AltModel->Repeat[t][i]) && 
	      (i < (B->IP->nMotifLen[t] / 2 + B->IP->nMotifLen[t] % 2)))
	     B->IP->DOF[t] -= B->IP->nAlphaLen - 1.0;
	 }

       if( ! B->IP->is_defined[cl_Z] )   
	 printf("DOF[%d] = %d\n", t, B->IP->DOF[t]);
     }   
}


/************************  motif_sampler **********************************/
/*                                                                        */
/*                                                                        */
/* RETURN VALUE:                                                          */
/*       Maximal alignment for the given seed                             */
/*                                                                        */
/* DESCRIPION : This function finds the maximal alignment for a given seed*/
/**************************************************************************/
MaxResults motif_sampler(register Model B, PoSition **Pos, Mlist M)
{
   ProbStruct P;
   MaxResults maxData;
   int n, t,  iter, last_increase=0, ***occurence;
   double dLocMax, dCurrProb, dMaxProbability=0.0;
   IPtype IP; /* for quick access of input data */
#ifdef _GUI_
   XEvent event;                        /* BT 6/17/97 */
#endif
   int    sampleCnt = 0;
   short  done = FALSE;
   short  samplePost = FALSE;
   
   IP=B->IP;
   init_maxdata(&maxData);
   /* occurrence[SeqLen][MotifType][ / /] */
   NEWPP(occurence, IP->nSeqLen, int);
   for(n = 0; n < IP->nSeqLen; n++) {
      NEWP(occurence[n], IP->nNumMotifTypes, int);
      for(t = 0; t < IP->nNumMotifTypes; t++) 
         NEW(occurence[n][t], 3, int);
   }

   /* initiate the values */
   iter= 0;
   last_increase = 0;
   dLocMax = -DBL_MAX;
   P.update = TRUE;
   for(t = 0; t < IP->nNumMotifTypes; t++) 
     IP->col_shift[t] = 0;
   init_prob(&P,IP); /* allocate space for probability */
       
   while( ! done )
     {
       iter++;
       if( ! B->IP->is_defined[cl_Z] )
	 {
	   fprintf(stdout, "\r%d", iter);
	   fflush( stdout );    /* BT 9/19/97 */
	 }
       
       /* BT 6/17/97  */
#ifdef _GUI_
       if( GL_GUI )
	 {
	   XmUpdateDisplay(toplevel);     /* BT 6/17/97 */
	   if( XCheckMaskEvent( XtDisplay( stopButton ), 
				ButtonPressMask | ButtonReleaseMask, &event ) )
	     {
	       if( event.xany.window == XtWindow( stopButton ) ) 
		 XtDispatchEvent( &event );
	     }
	 }
#endif       
       
       if( ! IP->is_defined[cl_D] )
	 sample_motif(B, Pos,M,iter, &occurence, &P);
       else
	 sample_motif_pair(B, Pos,M,iter, &occurence, &P);

       if(iter % IP->nAdjustPeriod == 0) 
	 {
	   if(!IP->is_defined[cl_F]) 
	     {  /* Run Fragmentation */
	       /* setFragVec(B, M, Pos);
	       setFragCnts(B, M);
	       StoreMask(B); */
	       for(t = 0; t < IP->nNumMotifTypes; t++)
		 {
		   StoreMask(B); 
		   setFragVec(B, M, Pos);   /* BT 4/26/04 */
		   setFragCnts(B, M);
		   Fragment(B, t, M);
		   adjust_fragment(B, M, Pos, &P);
		 }
	       /*		 Fragment(B, t, M);
				 adjust_fragment(B, M, Pos, &P); */
	     }
	   else  /* Do not fragment, but shift columns */
	     {
	       if( IP->is_defined[cl_d] )           /* BT 6/9/97 */
		 {
		   ResizeMotif( B, Pos, M, samplePost );      /* BT 5/12/97 */
		   for( t = 0; t < IP->nNumMotifTypes; t++ ) /* BT 11/24/97 */
	             B->C->dtot_sites[t] = IP->nPossSites[t];
		   if( samplePost )
		     sampleCnt++;
		 }
	       
	       ColumnShift(B, M, Pos, &P);  
	     }
	 } 

       if(iter > 3 &&  ! samplePost ) 
	 check_map(B, &dCurrProb, iter, Pos, &last_increase,
		   &dMaxProbability, &dLocMax, &maxData, P);  

#ifdef _DEBUG_
       if( last_increase == iter )
	 {
	   for( t = 0; t < IP->nNumMotifTypes; t++ )
	     {
	       fprintf( stdout, 
			"-----------------------------------------------------------\n" );
	       fprintf(stdout, "                          MOTIF %c\n\n", (char)(97 + t));
	       DumpMotifPositions( t, B, Pos, stdout );
	       fprintf( stdout, "Motif MAP = %g Frag Map = %g\n", 
			CalcMotifMap(B, t, B->IP->is_defined[cl_R]), 
			CalcMotifFragMap( B, t, B->IP->is_defined[cl_R] ));
	     }
	   fprintf( stdout, "iter = %d MAP = %g Bkgnd = %g Beta = %g seq = %g\n", 
		    iter, CalcMapProb(B, B->IP->is_defined[cl_R]),
		    CalcBkgndMap( B,  B->IP->is_defined[cl_R] ),
		    CalcBetaMap( B,  B->IP->is_defined[cl_R] ),
		    CalcSitePerSeqMap( B ) );
	 }
#endif

       if( iter >= IP->nMaxIterations ) 
	 done = TRUE;
       else if( iter - last_increase >= IP->nPlateauPeriods )
	 {
	   /* if( IP->is_defined[cl_X] &&
	       ! B->AN->bExchange &&
	       B->AN->currTemp > B->AN->dMinTemp )
	     {
	       B->AN->currTemp -= B->AN->dTempStep;
	       last_increase = iter;
	       if( ! B->IP->is_defined[cl_Z] )	    
		 printf( "Current temperature: %7.2f\n", B->AN->currTemp );
		 } 
		 else */ if( ! IP->is_defined[cl_d] || ! IP->is_defined[cl_q] ) 
	     done = TRUE;
	   else
	     {
	       if( sampleCnt >= IP->nSampleCnt )
		 done = TRUE;
	       samplePost = TRUE;
	     }
	 }
     }
   free_prob(&P, IP);
   /*   if(!IP->is_defined[cl_F])
	FREEPP(B->F->nvFragCnts, IP->nNumMotifTypes, B->F->nMaxLen[t]); */
   FREEPP(occurence, IP->nSeqLen, IP->nNumMotifTypes);
   maxData.frequency = NULL;
   return(maxData);
}

void print_Pos(PoSition **Pos, Model B, Mlist M)
{
   int t,j,n;
   PoSition *Pos_t;
   MlistEl  *curr;

   /* ???????? */
   for(t=0; t<B->IP->nNumMotifTypes; t++){
     j=0;
     Pos_t=Pos[t];
     for(n = 0; n < B->IP->nSeqLen; n++){ 
       if(Pos_t[n].nMotifStartPos) {
	 j++;
	 printf("%5d ",n);
	 if(j%10==0) printf("\n");
       }
     }
     printf("\n");
     printf("motif=%d\n",j);
     /* print start from the Motif list */
     j=0;
     curr=M[t]->Motifs;
     while(curr!=NULL){
       j++;
       printf("%5d ",curr->pos);
       if(j%10==0) printf("\n");
       curr=curr->next;
     }
     printf("\n");
     printf("motif from list =%d %d \n",j, M[t]->nNumMotifs);     
   }
}


/*************************   check_map ***************************************/
/*                                                                           */
/*                                                                           */
/* OUTPUT PARAMETERS :                                                       */
/*       dCurrProb      : MAP Probability value                              */
/*       C              : Current motif and background counts                */
/*       IP             : Parameters input at the command line (and updated) */
/*       iterations     : current iteration number                           */
/*       Pos            : Position structure                                 */
/*       last_increase  : Time period when last increase in the MAP occurred */
/*       dMaxProbability: Maximum MAP value                                  */
/*       dLocMax        : Maximum for the current seed run                   */
/*       maxData        : structure containing values of max MAP alignment   */
/*       dProbArray     : Probability for each element in the alignment      */
/*                                                                           */
/* DESCRIPION : This function finds the MAP value for the current alignment  */
/*              and resets the maximal information if it is the maximum      */
/*###########################################################################*/

void check_map(Model B, double *dCurrProb, int iterations, 
               register PoSition **Pos, int *last_increase,
               double *dMaxProbability, double *dLocMax, MaxResults *maxData, 
               ProbStruct P) 

{
   int      t, nNumMotifs;
   double   *dMotifMAPProb;
   register double **dProbArray;

#ifdef _DEBUG_
   if( B->IP->nSeedRun == 2 && iterations > 106 )
     iterations = iterations;
#endif

   (*dCurrProb) = CalcMapProb(B, B->IP->is_defined[cl_R]);

   if( *dCurrProb > *dLocMax ) 
     {
       if( ! B->IP->is_defined[cl_Z] )
	 {
	   if(B->IP->nNumMotifTypes > 1) 
	     {
	       NEW(dMotifMAPProb, B->IP->nNumMotifTypes, double);

	       for(t = 0; t < B->IP->nNumMotifTypes; t++)
		 {
		   dMotifMAPProb[t] = CalcMotifMap(B, t, B->IP->is_defined[cl_R]);
		 }

	       for(t = 0; t < B->IP->nNumMotifTypes; t++) {
		 if(!B->IP->is_defined[cl_F]) print_shift(B->F->shift[t]);
		 else                         print_shift(B->IP->col_shift[t]);
		 
		 nNumMotifs = NUMMOTIFS(B->IP->nNumMotifs[t]);
		 if( ! B->IP->is_defined[cl_Z] )
		   printf("motif %c cycle %d AP %.1f (%d sites)\n", 
			  (char)((int)'A'+t),iterations, dMotifMAPProb[t], 
			  nNumMotifs);
	       }
	       printf("Total Map : %g Prev: %g Diff: %g Motifs: %d\n", 
		      *dCurrProb, *dLocMax, *dCurrProb - *dLocMax, TotalNumMotifs( B ) );
	       free(dMotifMAPProb);   /* BT 8/27/97 */
	     }
	   else
	     {
	       if(!B->IP->is_defined[cl_F]) print_shift(B->F->shift[0]);
	       else                         print_shift(B->IP->col_shift[0]);
	       printf("cycle %d MAP %g (%d sites) Prev: %g Diff: %g Motifs: %d\n", iterations,
		      *dCurrProb, NUMMOTIFS(B->IP->nNumMotifs[0]),
		      *dLocMax, *dCurrProb - *dLocMax, TotalNumMotifs( B ) );
	     }
	 }
       dProbArray = setElementProb(B, Pos, P);
       free_maxdata(maxData, B->IP); 
       *maxData = setMaxData(B->F, B->IP, *dCurrProb, iterations, Pos, 
			     last_increase,dMaxProbability, dLocMax, 
			     *maxData, dProbArray, B);
       FREEP(dProbArray, findMaxNumMotif(B->IP));
     }
}


double **setElementProb(Model B, register PoSition **Pos, ProbStruct P)

   /*************************************************************************/ 
   /* FUNCTION NAME : setElementProb                                        */
   /*                                                                       */
   /* DESCRIPTION : When the maximum alignment is found, this will reset    */
   /*               the probability for each element given the whole        */
   /*               alignment                                               */
   /*************************************************************************/ 
  
{
   int             n, t, maxnum;
   register double **Prob;
   PoSition        *Pos_t;  /* Pos_t=Pos[t] */
   int             *count;
   double          **posProb;

   maxnum = findMaxNumMotif(B->IP);
   NEWP(Prob, maxnum, double);
   for(n = 0; n < maxnum; n++)
      NEW(Prob[n], B->IP->nNumMotifTypes, double);

   NEWP( posProb, B->IP->nNumMotifTypes, double );   
   for(t=0; t<B->IP->nNumMotifTypes; t++)
     NEW( posProb[t], 2, double );

   NEW(count, B->IP->nNumMotifTypes, int);

   update_prob( &P, B, TRUE );
   
   for(n = 0; n < B->IP->nSeqLen; n++)
     { 
       for(t=0; t<B->IP->nNumMotifTypes; t++)
	 {
	   Pos_t=Pos[t];
	   if(Pos_t[n].nMotifStartPos) 
	     {	       
	       CalcPosProb( B, P, Pos, n, posProb );
	       if( ! Pos_t[n].RevComp )
		 Prob[count[t]][t] = posProb[t][FORWARD];
	       else
		 Prob[count[t]][t] = posProb[t][REVERSE];
	       count[t]++;
	     }
	 }
     }

   FREEP( posProb, B->IP->nNumMotifTypes );
   free( count );
   
   return(Prob);
}


void CalcPosProb( Model B, ProbStruct P, PoSition **Pos, int n, double **posProb )
   /*************************************************************************/ 
   /* FUNCTION NAME : CalcPosProb                                           */
   /*                                                                       */
   /* DESCRIPTION : Calculates the probility of a site for each motif type  */
   /*               starting at pos n.                                      */
   /*************************************************************************/ 
{
   int             j, t;
   double          dMotifProb;
   double          dBGProb;
   double          **odds;
   double          denom;
   int             nSeq;
   int             nPos;
   int             bPairs;
   int             treeSave;

   NEWP(odds, B->IP->nNumMotifTypes, double);
   for( t = 0; t < B->IP->nNumMotifTypes; t++ )
     NEW( odds[t], 2, double );
   
   for( t = 0; t < B->IP->nNumMotifTypes; t++ )
     {
       posProb[t][FORWARD] = 0;
       posProb[t][REVERSE] = 0;
     }

   bPairs = B->IP->is_defined[cl_D]; /* turn off pairwise calculation of probs */
   B->IP->is_defined[cl_D] = FALSE;  /* we only want to calculate based on position n*/
   treeSave = B->Phylo->treeCount;
   B->Phylo->treeCount = 0;

   if( B->RP->nUseSpacingProb || B->IP->is_defined[cl_T] )   /* 9/27/2000 */
     {
       denom = 1.0;
       nSeq = SequenceFromPosition( B, n );
       nPos = n - SequenceStartPos( B, nSeq );
       for(t=0; t < B->IP->nNumMotifTypes; t++)
	 denom -= B->RP->sitePos[nSeq][nPos][t].dSpacingProb;

       for(t=0; t < B->IP->nNumMotifTypes; t++)
	 {
	   if( PossFragStartPos(Pos, n, t, B) )
	     {
	       for( j = FORWARD; j <= REVERSE; j++ )
		 {
		   in_motif_prob(B, P, n, t, j, FALSE, &dMotifProb, &dBGProb);
		   odds[t][j] = (dMotifProb / dBGProb) * B->RP->sitePos[nSeq][nPos][t].dSpacingProb;
		   denom += odds[t][j];
		   if( ! B->IP->RevComplement )
		     break;
		 }
	     }	   
	 }
       
       for(t=0; t < B->IP->nNumMotifTypes; t++)
	 {
	   for( j = FORWARD; j <= REVERSE; j++ )
	     posProb[t][j] = odds[t][j] / denom;
	 }
     }
   else
     {       
       nSeq = SequenceFromPosition( B, n );

       denom = 1.0;
       for(t=0; t < B->IP->nNumMotifTypes; t++)
	 denom -= B->IP->dposterior_prob[t];

       for(t=0; t < B->IP->nNumMotifTypes; t++)
	 {
	   if( GetSeqWeight( B, nSeq, t ) == 0 )
	     {
	       odds[t][FORWARD] = 0.0;
	       odds[t][REVERSE] = 0.0;
	     }	     
	   else if( PossFragStartPos(Pos, n, t, B) )
	     {
	       for( j = FORWARD; j <= REVERSE; j++ )
		 {
		   in_motif_prob(B, P, n, t, j, FALSE, &dMotifProb, &dBGProb);
		   odds[t][j] = (dMotifProb / dBGProb) * B->IP->dposterior_prob[t] ;
		   denom += odds[t][j];
		   if( ! B->IP->RevComplement )
		     break;
		 }
	     }       
	 }

       for(t=0; t < B->IP->nNumMotifTypes; t++) /* BT 12/13/2000 */
	 {
	   for( j = FORWARD; j <= REVERSE; j++ )
	     posProb[t][j] = odds[t][j] / denom;
	 }	   
     }

   FREEP( odds, B->IP->nNumMotifTypes );
   B->IP->is_defined[cl_D] = bPairs;
   B->Phylo->treeCount = treeSave;
}


/***************************************************************/
/* FUNCTION NAME : setMaxData                                  */
/*                                                             */
/* DESCRIPTION : This function sets the data for the maximum   */
/*               alignment that has been found                 */
/***************************************************************/

MaxResults setMaxData(Ftype F, IPtype IP, double dProb, int iteration,
                      register PoSition **Pos, int *last_increase, 
                      double *dMaxProb,
                      double *dLocMax, MaxResults maxData, 
                      register double **dProbArray,
		      Model B )
{
   MaxResults temp;
   int        n, i, j, t, maxlen;
   PoSition   *Pos_t;  /* Pos_t=Pos[t] */
   int        phyloSave;

   init_maxdata(&temp);
   maxlen = findMaxNumMotif(IP);

   temp.nMaxLen = maxlen;

   NEW(temp.nNumMotifs, IP->nNumMotifTypes, int);
   NEW(temp.nMotifLen, IP->nNumMotifTypes, int);     /* BT 5/28/97 */
   NEW(temp.dMap, IP->nNumMotifTypes, double);       /* BT 3/4/98 */
   NEW(temp.dFragMap, IP->nNumMotifTypes, double);       /* BT 3/4/98 */
   NEWP(temp.nMotifLoc, maxlen, int);
   NEWP(temp.dvMotifProb, maxlen, double);
   NEWP(temp.RevComp, maxlen, short);
   NEW(temp.F, 1, FragStruct);
   temp.F->nMaxLen = NULL;
   temp.F->nvFragCnts = NULL;

   if(!IP->is_defined[cl_F]) {
      NEW(temp.F->FragWidth, IP->nNumMotifTypes, int);
      NEW(temp.F->nMaxLen, IP->nNumMotifTypes, int);
      NEWP(temp.F->nColMask, IP->nNumMotifTypes, int);
      NEWP(temp.F->fragPos, IP->nNumMotifTypes, int);
      for(t = 0; t < IP->nNumMotifTypes; t++)  {
         temp.F->FragWidth[t] = F->FragWidth[t];
         NEW(temp.F->nColMask[t], F->nMaxLen[t] , int);
         NEW(temp.F->fragPos[t], IP->nMotifLen[t] , int);
	 temp.F->nMaxLen[t] = F->nMaxLen[t];
         for(i = 0; i < F->nMaxLen[t]; i++) 	   
	   temp.F->nColMask[t][i] = F->nColMask[t][i];
         for(i = 0; i < IP->nMotifLen[t]; i++) 	   
	   temp.F->fragPos[t][i] = F->fragPos[t][i];	  	     
      }
   }
   for (i = 0; i < maxlen; i++) {
      NEW(temp.nMotifLoc[i], IP->nNumMotifTypes, int);
      NEW(temp.dvMotifProb[i], IP->nNumMotifTypes, double);
      NEW(temp.RevComp[i], IP->nNumMotifTypes, short);
   }
   *dLocMax = *dMaxProb = temp.dProbability = dProb;
   *last_increase = temp.nIterationNum = iteration;
   temp.nseed = IP->lSeedVal;
   for(t = 0; t < IP->nNumMotifTypes; t++) {
      temp.nNumMotifs[t] = NUMMOTIFS(IP->nNumMotifs[t]);
      temp.nMotifLen[t] = IP->nMotifLen[t];      /* BT 5/23/97 */
   }
   
   for(t=0;t < IP->nNumMotifTypes; t++){
     Pos_t=Pos[t];
     j=0;
     for(n = 0; n < IP->nSeqLen; n++) {
       if(Pos_t[n].nMotifStartPos){
	 temp.nMotifLoc[j][t] = n;
	 temp.RevComp[j][t]= Pos_t[n].RevComp;
	 j++;
	 /* skip the motif length */
	 n+=(IP->nMotifLen[t]-1);
       }
     }
   }
 
   for(t = 0; t < IP->nNumMotifTypes; t++) 
      for(i = 0; i < temp.nNumMotifs[t]; i++)
         temp.dvMotifProb[i][t] = dProbArray[i][t];

   for( t = 0; t < IP->nNumMotifTypes; t++)  /* BT 3/4/98 */
     {
       temp.dMap[t] = CalcMotifMap(B, t, B->IP->is_defined[cl_R] );
       temp.dFragMap[t] = CalcMotifFragMap(B, t, B->IP->is_defined[cl_R] );
     }

   temp.dTotalMap = CalcMapProb( B, B->IP->is_defined[cl_R] );
   temp.dBkgndMap = CalcBkgndMap( B, B->IP->is_defined[cl_R] );
   temp.dBetaMap = CalcBetaMap( B, B->IP->is_defined[cl_R] );
   temp.dSeqMap = CalcSitePerSeqMap( B );
   temp.dTotalNonPalinMap = CalcMapProb( B, FALSE );
   if( B->Phylo->phyloTree )
     {
       phyloSave = B->Phylo->bCalcPhylo;
       B->Phylo->bCalcPhylo = FALSE;
       temp.dTotalNonPhyloMap = CalcMapProb( B, B->IP->is_defined[cl_R] ) + B->Phylo->phyloNull;
       B->Phylo->bCalcPhylo = phyloSave;
     }

   return(temp);
}


void reset_values(int *iterations, int *last_increase, double *dLocMax,
                  ProbStruct *P, IPtype IP)
{
  int t;

  *iterations = 0;
  *last_increase = 0;
  *dLocMax = -DBL_MAX;
  P->update = TRUE;
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    IP->col_shift[t] = 0;
  init_prob(P, IP);
}

/***************** void init_values ********************************/

void init_values(double ***dFinalProb, PoSition **Pos, int n, int t)
{
   (*dFinalProb)[t][0] = 0.0;
   (*dFinalProb)[t][1] = 0.0;
   Pos[t][n].nInMotif = FALSE;
}


void setMotifNum(IPtype IP, MaxResults maxData)
{
   int i, t;

   for(t = 0; t < IP->nNumMotifTypes; t++)  {
      IP->nNumMotifs[t][FORWARD] = IP->nNumMotifs[t][REVERSE] = 0;
      for(i = 0; i < maxData.nNumMotifs[t]; i++) 
        IP->nNumMotifs[t][maxData.RevComp[i][t]]++;
   }
}

void PrintPos( Model B, PoSition **Pos, int t )
{
  int i;

  for( i = 0; i < B->IP->nSeqLen; i++ )
    {
      if( Pos[t][i].nMotifStartPos )
	fprintf( stdout, "%5d %2d", i, Pos[t][i].nSeq );
    }
  fprintf( stdout, "\n" );
}


void PrintDoubleAsBytes( double *d )
{
  int            i;
  unsigned char  *p = (unsigned char *) d;
 
  for( i = 0; i < 8; i++ )
    {
      printf( "%x ", p[i] );
    }
}
