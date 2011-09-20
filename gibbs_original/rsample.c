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
#include "rsample.h"

#define WAIT_PERIOD_SITE  100
#define PLOT_SEQ          4
#define TARGET_SEED       8
#define MIN_ITER          27
#define MAX_ITER          32
#define MAX_TRYS          100000
#define MOTIF_NUM         0
#define SITE_WAIT_PERIOD  5
#define MIN_SITE_PROB     1e-10
#define MIN_SITE_MOTIFS   160
#define MOTIF_SAMPLE_PERIOD 100

double SampleSites( Model B, PoSition **Pos, Mlist M, int ****occurence, 
		    ProbStruct *P, int iter, int seed_run, short bAccumulateCounts,
		    double *dMap );
double SampleSitesFromAllK( Model B, PoSition **Pos, Mlist M, int ****occurence, 
			    ProbStruct *P, int iter, int seed_run, short bAccumulateCounts,
			    double *dMap );
int SampleMaxBlocks( Model B, int seq, int iter, int seed_run, short bAccumulateCounts );
int SampleMaxBlocksK( Model B, int seq, int iter, int seed_run, short bAccumulateCounts );
int SampleMaxBlocksHier( Model B, int seq, int iter, int seed_run, short bAccumulateCounts );
double CalcAlignmentProb( Model B, PoSition **Pos, Mlist M, ProbStruct *P, double *prob );
void UpdateBackgroundModel( Model B, PoSition **Pos, ProbStruct *P );
void RemoveAllSites( Model B, PoSition **Pos, Mlist M, int ****occurence, 
		     ProbStruct *P );
void DumpBackground( Model B );
void AdjustSiteSeqCounts( Model B, PoSition **Pos, Mlist M );
void RemoveSites( Model B, int seq,  PoSition **Pos, Mlist M );

#ifndef __LONGDOUBLE128
extern long double expl (long double __x) ; 
extern long double powl (long double __x, long double __y) ; 
#endif


/*************    MaxResult rsite_sampler *********************************/
/*                                                                       */
/*    Pos[MotifType][SeqLen]                                             */
/*=======================================================================*/

MaxResults rsite_sampler(Model B, PoSition **Pos, Mlist M, int seed_run)
 
{
   int             iter, last_increase;
   int             n, t, i;
   double          dLocMax = -DBL_MAX, dCurrProb, dMaxProbability = 0.0;
   double          **dProb, **dFinalProb;
   ProbStruct      P;
   MaxResults      maxData;
#ifdef _GUI_
   XEvent          event;                        /* BT 6/17/97 */
#endif
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
#ifdef _DEBUG_   
   double          map1, map2;
   double          diffTotal = 0.0;
   double          b1, m1, f1, bt1, sq1;
   double          b2 ,m2, f2, bt2, sq2;
   double          prf[GRAPH_SIZE];   /* DEBUG */
   int             nGrCnt = 0;
#endif   
#ifdef _MPI_
   double         prevTemp;
#endif

#ifdef _DEBUG_   
   memset( prf, 0, sizeof( prf ) );
#endif

   IP = B->IP;
   RP = B->RP;

   if( IP->is_defined[cl_V] )
     {
       RP->nIncludeAlign = TRUE;
       UpdateTransMatrix( B );   /* BT 1/22/03 */
     }
   else
     {
       RP->nIncludeAlign = FALSE; 
     }

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
   
   reset_values(&iter, &last_increase, &dLocMax, &P, B->IP); 
   update_prob(&P, B, (! B->IP->is_defined[cl_b]));

   if( B->RP->bUsePosMatrix )
     InitializePosMatrix( B );
   else if( B->RP->bUseTrans )
     InitializeTransMatrix( B );

   if( B->InitPos != NULL && seed_run == 1 )
     {
       if( (! B->IP->is_defined[cl_Z]) && B->IP->is_defined[cl_opt] )
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

   /* ======================================================================= */
   if( IP->is_defined[cl_L] && ! B->Phylo->phyloTree )
     {
       IP->is_defined[cl_E] = 0;
       while( TRUE )
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
	       if(!B->IP->is_defined[cl_F]) 
		 {
		   /* setFragVec(B, M, Pos);
		      setFragCnts(B, M); 
		      StoreMask(B); */
		   bCount = FALSE;
#ifdef _DEBUG_   
		   map1 =  CalcMapProb(B, B->IP->is_defined[cl_R]);
		   b1 =CalcBkgndMap(B, B->IP->is_defined[cl_R]);
		   m1 =CalcMotifMap(B, 0, B->IP->is_defined[cl_R]);
		   f1 =CalcMotifFragMap(B, 0, B->IP->is_defined[cl_R]);
		   bt1 = CalcBetaMap( B, B->IP->is_defined[cl_R]);
		   sq1 = CalcSitePerSeqMap( B );
		   printf("bk = %f f = %f m = %f bt = %f seq = %f sum = %f\n", 
			  b1, f1, m1, bt1, sq1, b1 + f1 + m1 + bt1 + sq1 );
#endif
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
#ifdef _DEBUG_   
		   map2 =  CalcMapProb(B, B->IP->is_defined[cl_R]);
		   b2 =CalcBkgndMap(B, B->IP->is_defined[cl_R]);
		   m2 =CalcMotifMap(B, 0, B->IP->is_defined[cl_R]);
		   f2 =CalcMotifFragMap(B, 0, B->IP->is_defined[cl_R]);
		   bt2 = CalcBetaMap( B, B->IP->is_defined[cl_R]);
		   sq2 = CalcSitePerSeqMap( B );
		   printf("bk = %f f = %f m = %f bt = %f seq = %f sum = %f\n", 
			  b2, f2, m2, bt2, sq2, b2 + f2 + m2 + bt2 + sq2 );
		   
		   diffTotal += map2 - map1;
		   printf( "before = %f after = %f diff = %f total = %f\n", 
			   map1, map2, map2 - map1, diffTotal );
#endif  
		 }
	       else
		 {
		   if( B->IP->is_defined[cl_d] )           
		     ResizeMotif( B, Pos, M, FALSE );            
		   ColumnShift(B, M, Pos, &P);        
		 }
	     }
	   
	   if(iter > WAIT_PERIOD ) 
	     check_map(B, &dCurrProb, iter, Pos, &last_increase,
		       &dMaxProbability, &dLocMax, &maxData, P);
	   
	   if( dLocMax > 0 )
	     break;
	   
	   if( iter >= IP->nMaxIterations ) 
	     break;
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
		 break;
	     }
	 }
       
       IP->is_defined[cl_E] = 1;
       AdjustSiteSeqCounts( B, Pos, M ); 
       dCurrProb = CalcMapProb(B, B->IP->is_defined[cl_R]);
       dProbArray = setElementProb(B, Pos, P);
       free_maxdata(&maxData, B->IP);
       maxData = setMaxData(B->F, B->IP, dCurrProb, iter, Pos, &last_increase,
			    &dMaxProbability, &dLocMax, maxData, dProbArray, B);
       FREEP(dProbArray, findMaxNumMotif(B->IP));
       
#ifdef _DEBUG_
       for( t = 0; t < B->IP->nNumMotifTypes; t++ )
	 {
	   fprintf( stdout, 
		    "-----------------------------------------------------------\n" );
	   fprintf(stdout, "                          MOTIF %c\n\n", (char)(97 + t));
	   DumpMotifPositions( t, B, Pos, stdout );
	   fprintf( stdout, "Motif MAP = %g\n", 
		    CalcMotifMap(B, t, B->IP->is_defined[cl_R]));
	 }
       fprintf( stdout, "iter = %d MAP = %g\n", 
		0, CalcMapProb(B, B->IP->is_defined[cl_R]));
#endif
       
     }
   /* =============================================== */

   while( TRUE )
     {
       if( iter > WAIT_PERIOD || IP->is_defined[cl_V] )
	 {
	   if( iter >= IP->nMaxIterations ) 
	     break;
	   else if( iter - last_increase >= IP->nPlateauPeriods )
	     break;
	 }
       
       iter++;
       if( ! B->IP->is_defined[cl_Z] )
	 {
	   fprintf(stdout, "\r%d", iter); 
	   fflush( stdout );    /* BT 9/19/97 */
	 }
       
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

       dCurrProb = rsample( B, Pos, M, &occurence, &P, iter, 
				seed_run, FALSE, &dMap, dLocMax );	
       
#ifdef _DEBUG_
       if( nGrCnt < GRAPH_SIZE )
	 {
	   prf[nGrCnt] = dCurrProb;
	   nGrCnt++;
	 }
#endif

       if( dCurrProb > dLocMax && iter > WAIT_PERIOD )
	 {
	   if( (! B->IP->is_defined[cl_Z]) && B->IP->is_defined[cl_opt] )
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
	       /* setFragVec(B, M, Pos);
		  setFragCnts(B, M); 
		  StoreMask(B); */
	       bCount = FALSE;
#ifdef _DEBUG_   
	       map1 =  CalcMapProb(B, B->IP->is_defined[cl_R]);
	       b1 =CalcBkgndMap(B, B->IP->is_defined[cl_R]);
	       m1 =CalcMotifMap(B, 0, B->IP->is_defined[cl_R]);
	       f1 =CalcMotifFragMap(B, 0, B->IP->is_defined[cl_R]);
	       bt1 = CalcBetaMap( B, B->IP->is_defined[cl_R]);
 	       sq1 = CalcSitePerSeqMap( B );
	       printf("bk = %f f = %f m = %f bt = %f seq = %f sum = %f\n", 
		      b1, f1, m1, bt1, sq1, b1 + f1 + m1 + bt1 + sq1 );

	       for( t = 0; t < B->IP->nNumMotifTypes; t++ )
		 {
		   fprintf( stdout, 
			    "-----------------------------------------------------------\n" );
		   fprintf(stdout, "                          MOTIF %c\n\n", (char)(97 + t));
		   DumpMotifPositions( t, B, Pos, stdout );
		   fprintf( stdout, "Motif MAP = %g\n", 
			    CalcMotifMap(B, t, B->IP->is_defined[cl_R]));
		 }
	       fprintf( stdout, "iter = %d MAP = %g\n", 
			0, CalcMapProb(B, B->IP->is_defined[cl_R]));
#endif
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
#ifdef _DEBUG_   
	       map2 =  CalcMapProb(B, B->IP->is_defined[cl_R]);
	       b2 =CalcBkgndMap(B, B->IP->is_defined[cl_R]);
	       m2 =CalcMotifMap(B, 0, B->IP->is_defined[cl_R]);
	       f2 =CalcMotifFragMap(B, 0, B->IP->is_defined[cl_R]);
	       bt2 = CalcBetaMap( B, B->IP->is_defined[cl_R]);
 	       sq2 = CalcSitePerSeqMap( B );
	       printf("bk = %f f = %f m = %f bt = %f seq = %f sum = %f\n", 
		      b2, f2, m2, bt2, sq2, b2 + f2 + m2 + bt2 + sq2 );

	       diffTotal += map2 - map1;
	       printf( "before = %f after = %f diff = %f total = %f\n", 
		       map1, map2, map2 - map1, diffTotal );
#endif  
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

       if(iter > WAIT_PERIOD ) 
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
   
   if( IP->is_defined[cl_Q] )
     {
       update_prob(&P, B, TRUE);
              
       if( ! B->IP->is_defined[cl_Z] )
	 fprintf( stdout, "\nSampling\n" );
       for(i = 1; i <= B->IP->nPostPlateauIter; i++) 
	 {
	   if( ! B->IP->is_defined[cl_Z] && i % 5 == 0)
	     {
	       fprintf(stdout, "\r%d", i);
	       fflush( stdout );  
	     }
	   
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

	   dCurrProb = rsample( B, Pos, M, &occurence, &P, i, seed_run, TRUE, &dMap, dLocMax );
	   if( dCurrProb > dLocMax )
	     {
	       if( (! B->IP->is_defined[cl_Z]) && B->IP->is_defined[cl_opt] )
		 {
		   time( &tloc );
		   localTime = localtime_r( &tloc, &tm_buf );
		   
		   fprintf( stdout, "Max set %f at %d motifs: %d %d:%d:%d\n", 
			    dCurrProb, iter + i, TotalNumMotifs( B ),
			    localTime->tm_hour, localTime->tm_min, localTime->tm_sec);
		 }
	       dProbArray = setElementProb(B, Pos, P);
	       free_maxdata(&maxData, B->IP);
	       maxData = setMaxData(B->F, B->IP, dCurrProb, iter + i, Pos, &last_increase,
				    &dMaxProbability, &dLocMax, maxData, dProbArray, B);
	       FREEP(dProbArray, findMaxNumMotif(B->IP));
	     } 

#ifdef _MPI_
	   if( IP->is_defined[cl_X] && 
	       B->AN->bExchange )
	     {
	       currTime = (times( &timeBuffer ) / IP->ticks) * 1000;
	       if( ! B->IP->is_defined[cl_Z] )
		 PrintTempOut( IP->Datafiles->mpiTemp_fpt, 
			       "currTime = %d lastExchange = %d exchangePeriod = %d diff = %d\n", 
			       currTime,  B->AN->lastExchangeTime, B->AN->exchangePeriod,
			       currTime - B->AN->lastExchangeTime );   /* DEBUG */
	       if( (B->AN->exchangePeriod == 0 && ((iter + i) % B->AN->nExchangeIterations == 0)) ||
		   (B->AN->exchangePeriod > 0 && (currTime - B->AN->lastExchangeTime > B->AN->exchangePeriod)) )
		 {
		   ExchangeTemps( B, dCurrProb, seed_run, iter + i );
		   currTime = (times( &timeBuffer ) / IP->ticks) * 1000;
		   B->AN->lastExchangeTime = currTime;
		 }    /* BT 1/29/2001 */
	     }
	   
	   if( IP->is_defined[cl_hm] && (iter + i) % IP->hier_iter == 0 )
	     ExchangeK( B, seed_run, iter + i );
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
   FREEPP(occurence, IP->nSeqLen, IP->nNumMotifTypes);

   return(maxData);
}


MaxResults nearopt_rsite_samp(Model B, PoSition **Pos,  Mlist M, int **good, 
			      ProbStruct *P, int ****occurence,
			      double dSubOptProb )
   /*=======================================================================*/
   /* FUNCTION NAME : nearopt_site_samp                                     */
   /*                                                                       */
   /*=======================================================================*/
{
   int            i=0, n, t, last_increase = 0;  /* BT 1/27/97 */
   double         dCurrProb, dLocMax = -DBL_MAX, dMaxProbability=0.0;
   double         **dProb, **dFinalProb, **dProbArray;
   MaxResults     currMax;
#ifdef _GUI_
   XEvent         event;                        /* BT 6/17/97 */
#endif
   int            **posSave;			/* BT 3/26/97 */
   double         dMap;
   struct tm      *localTime;
   struct tm      tm_buf;
   time_t         tloc;
#ifdef _DEBUG_   
   double          prf[GRAPH_SIZE];   /* DEBUG */
   int             nGrCnt = 0;
#endif   

#ifdef _DEBUG_   
   memset( prf, 0, sizeof( prf ) );
#endif

   B->IP->is_defined[cl_b] = FALSE;
   if( B->IP->is_defined[cl_X] )
     B->AN->currTemp = 1.0;

   init_maxdata(&currMax);
   NEWP(dFinalProb, B->IP->nSeqLen, double);
   for(n = 0; n < B->IP->nSeqLen; n++)
      NEW(dFinalProb[n], 2, double);
   
   NEWP(dProb, B->IP->nSeqLen, double);
   for(n = 0; n < B->IP->nSeqLen; n++)
      NEW(dProb[n], 2, double);
   
   if( B->IP->is_defined[cl_V] )
     {
       B->RP->nIncludeAlign = TRUE;
     }
   /* else
     {
       B->RP->nIncludeAlign = FALSE;
       } */  /* BT 7/10/2002 */

   dCurrProb = dSubOptProb;   

   dProbArray = setElementProb(B, Pos, *P);
   currMax = setMaxData(B->F, B->IP, dCurrProb, i, Pos, &last_increase,
			&dMaxProbability, &dLocMax, currMax, dProbArray, B);
   FREEP(dProbArray, findMaxNumMotif(B->IP));

#ifdef _DEBUG_    
    CalcMapProb( B, B->IP->is_defined[cl_R] );
#endif    

   NEWP( posSave, B->IP->nNumMotifTypes, int);                /* BT 3/26/97 */
   for(i = 0; i < B->IP->nNumMotifTypes; i++)
      NEW( posSave[i], B->IP->nSeqLen, int);

   for(t = 0; t < B->IP->nNumMotifTypes; t++)
     {         
       for(n = 0; n < B->IP->nSeqLen; n++)
	 {
	   posSave[t][n] = Pos[t][n].nPossStartPos; 
	   if(!good[n][t])
	     {
	       Pos[t][n].nPossStartPos = FALSE;
	     }
	 }
     }

   if( B->RP->bUsePosMatrix )
     {
       InitializePosMatrix( B );
       if( B->RP->bUpdatePosMatrix )
	 UpdatePosMatrix( B );
     }
   else if( B->RP->bUseTrans )
     {
       InitializeTransMatrix( B );
       if( B->RP->bUpdateTrans )
	 UpdateTransMatrix( B );
     }

   last_increase = 0;
   for(n = 0; n < B->IP->nSeqLen; n++) 
     {
       for(t = 0; t < B->IP->nNumMotifTypes; t++) 
	 {
	   for( i = 0; i < 3; i++ )
	     (*occurence)[n][t][i] = 0;
	 }
     }
   
   if( ! B->IP->is_defined[cl_Z] )
     {
       fprintf(stdout, "Sampling\n" ); 
       fflush( stdout );   
     }
   for(i = 1; i <= B->IP->nMaxIterations; i++) 
     {
       if( (! B->IP->is_defined[cl_Z]) && i % 5 == 0)
	 {
	   fprintf(stdout, "\r%d", i);
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
       
       dCurrProb = rsample( B, Pos, M, occurence, P, i, 0, B->IP->is_defined[cl_Q], &dMap, dLocMax );  /* BT 07/05/05 */
       
#ifdef _DEBUG_
       if( nGrCnt < GRAPH_SIZE )
	 {
	   prf[nGrCnt] = dCurrProb;
	   nGrCnt++;
	 }
#endif
       if(dCurrProb > dLocMax)      /* BT 7/30/97 */
	 {
	   if( (! B->IP->is_defined[cl_Z]) && B->IP->is_defined[cl_opt] )
	     {
	       time( &tloc );
	       localTime = localtime_r( &tloc, &tm_buf );
	       
	       fprintf( stdout, "Max set %f at %d motifs: %d %d:%d:%d\n", 
			dCurrProb, i, TotalNumMotifs( B ),
			localTime->tm_hour, localTime->tm_min, localTime->tm_sec);
	     }
	   dProbArray = setElementProb(B, Pos, *P);
	   free_maxdata(&currMax, B->IP);
	   currMax = setMaxData(B->F, B->IP, dCurrProb, i, Pos, &last_increase,
				&dMaxProbability, &dLocMax, currMax, dProbArray, B);
	   FREEP(dProbArray, findMaxNumMotif(B->IP));
	 }
     }
   FREEP(good, B->IP->nSeqLen);
   FREEP(dProb, B->IP->nSeqLen);
   
   /* Restore possible site starting positions for next phase */
   for(t = 0; t < B->IP->nNumMotifTypes; t++)		/* BT 3/26/97 */
   {         
     for(n = 0; n < B->IP->nSeqLen; n++)
       Pos[t][n].nPossStartPos = posSave[t][n]; 
   }
   FREEP( posSave, B->IP->nNumMotifTypes);
   
   FREEP(dFinalProb, B->IP->nSeqLen);

   return (currMax);
}


double rsample( Model B, PoSition **Pos, Mlist M, int ****occurence, 
		ProbStruct *P, int iter, int seed_run, short bAccumulateCounts,
		double *dMap, double bestMap )
{
  if( B->IP->is_defined[cl_K] && (bestMap >= B->RP->dKSampleMap || B->IP->is_defined[cl_V]) ) 
    return SampleSitesFromAllK( B, Pos, M, occurence,
				P, iter, seed_run, bAccumulateCounts,
				dMap );
  else
    return SampleSites( B, Pos, M, occurence,
			P, iter, seed_run, bAccumulateCounts,
			dMap ); 
}


double SampleSites( Model B, PoSition **Pos, Mlist M, int ****occurence, 
		ProbStruct *P, int iter, int seed_run, short bAccumulateCounts,
		double *dMap )
{
  int      seq, first, last;
  int      t;
  int      n;
  int      kmax;
  short    RevComp;                      /* BT 8/15/97 */
  IPtype   IP;
  RPType   RP;
  int      k;
  int      lastpos;
  int      motif_type = 0;
  int      index;
  double   rand_num;
  double   dtot;
  int      type = -1;
  double   denom;
  int      cnt = 0;
  int      tcnt;
  double   psum = 0;
  double   dCurrProb = 0;
  double   dProb;
  double   prob;
  int      alignSum = 0;
  short    bMaxSet = FALSE;
  int      nWidth;
  int      nOffset;
  PoSition **savePos = 0;
  int      nSitesSelected;
  double   *motifMap = 0;
  int      p;
  int      j;

  static double psumMax = -DBL_MAX;
  static double dProbMax = -DBL_MAX;
  static int    maxIter; 
  static int    maxIter3; 

#ifdef _DEBUG_
  double     prMax;
  int        nMaxPos;
  double     prf[GRAPH_SIZE];   /* DEBUG */
  double     prr[GRAPH_SIZE];
  double     fpr[GRAPH_SIZE];
  int        nPrevMotifPos;
#endif
   
  IP = B->IP;
  RP = B->RP;

  B->AP->dTotalAlgnCnt = 0;
  for(seq = 0; seq < B->IP->nNumSequences; seq += SpeciesInc( B, seq ) ) 
    {
      B->AP->dTotalAlgnCnt += B->AP->dAlignCount[RP->nSites[seq]][seq];
    }

  if( IP->is_defined[cl_b] )
    UpdateBackgroundModel( B, Pos, P );

  if( IP->is_defined[cl_V] )
    {
      update_posterior_prob(B);   /* 4/18/03 */
      NEWP( savePos, B->IP->nNumMotifTypes, PoSition);
      NEW( motifMap, B->IP->nNumMotifTypes, double);
      for( t=0; t < B->IP->nNumMotifTypes; t++ )
	{
	  NEW( savePos[t], B->IP->nSeqLen, PoSition );
	  for( n = 0; n < B->IP->nSeqLen; n++ )
	    {
	      savePos[t][n] = Pos[t][n];
	    }
	}      
    }

  for(seq = 0; seq < B->IP->nNumSequences; seq += SeqInc( B, seq ) ) 
    {      
#ifdef _DEBUG_
      nPrevMotifPos = -1;
#endif

      if( IP->is_defined[cl_V] && seq >= IP->nVerifySeq )   /* BT 04/30/03 */
	break;

      first = SequenceStartPos( B, seq );
      nOffset = SequenceLength( B, seq );
      
      if( IP->is_defined[cl_sample_model] )
	{
	  SampleModel( B, P);
	  RemoveSites( B, seq, Pos, M );
	}
      else
	{
	  RemoveSites( B, seq, Pos, M );
	  update_prob(P, B, (! B->IP->is_defined[cl_b]));	
	}

      B->AP->dTotalAlgnCnt = 0;
      for(n = 0; n < B->IP->nNumSequences; n++) 
	{
	  if( n != seq && B->AP->dAlignCount[RP->nSites[n]][n] != -DBL_MAX )
	    B->AP->dTotalAlgnCnt += B->AP->dAlignCount[RP->nSites[n]][n];
	}
      if( B->AP->dTotalAlgnCnt == 0 )
	B->AP->dTotalAlgnCnt = 1;

      InitRProbArray( B, seq );     /* BT 12/29/97 */
      CalcRProbArray( B,  Pos, seq, P, iter );
  
      kmax = SampleMaxBlocks( B, seq, iter, seed_run, bAccumulateCounts );
      
      lastpos = SequenceEndPos( B, seq );
      for( k = kmax; k > 0; k-- )
	{
#ifdef _DEBUG_
	  for( t = 0; t < IP->nNumMotifTypes; t++ )
	    {
	      memset( prf, 0, sizeof( prf ) );
	      memset( prr, 0, sizeof( prr ) );
	      memset( fpr, 0, sizeof( fpr ) );
	      prMax = -DBL_MAX;
	      for( n = first; n < lastpos; n++ ) 
		{		  
		  if( n - first >= GRAPH_SIZE )
		    break;

		  fpr[n-first] = RP->dPFootprint[t][n-first];
		  		  
		  if( PossStartPos( Pos, n, t, B ) )
		    {
		      if( k == kmax || type == -1 )
			{
			  if(  RP->bUsePosMatrix )
			    {
			      prf[n-first] = RP->dProb[k][t][FORWARD][n - first] * 
				RP->dPFootprint[t][n-first] * RP->dPosMatrix[k-1][t];
			      prr[n-first] = RP->dProb[k][t][REVERSE][n - first] * 
				RP->dPFootprint[t][n-first] * RP->dPosMatrix[k-1][t];
			    }
			  else
			    {
			      prf[n-first] = RP->dProb[k][t][FORWARD][n - first] * 
				RP->dPFootprint[t][n-first] * RP->dEndSiteProb[t];
			      prr[n-first] = RP->dProb[k][t][REVERSE][n - first] * 
				RP->dPFootprint[t][n-first] * RP->dEndSiteProb[t];
			    }
			}
		      else
			{
			  nWidth = MotifWidth( B, t );
			  if(   RP->bUsePosMatrix )
			    {
			      prf[n-first] =
				RP->dPGap[t][motif_type][lastpos-n-nWidth + 1] *
				RP->dPosMatrix[k-1][t] *
				RP->dProb[k][t][FORWARD][n - first] * RP->dPFootprint[t][n-first];
			      prr[n-first] = 
				RP->dPGap[t][motif_type][lastpos-n-nWidth + 1] *
				RP->dPosMatrix[k-1][t] *
				RP->dProb[k][t][REVERSE][n - first] * RP->dPFootprint[t][n-first];
			    }
			  else
			    {
			      prf[n-first] =
				RP->dPGap[t][motif_type][lastpos-n-nWidth + 1] *
				RP->dPTrans[motif_type][type][t][FORWARD] *
				RP->dProb[k][t][FORWARD][n - first] * RP->dPFootprint[t][n-first];
			      prr[n-first] = 
				RP->dPGap[t][motif_type][lastpos-n-nWidth + 1] *
				RP->dPTrans[motif_type][type][t][REVERSE] *
				RP->dProb[k][t][REVERSE][n - first] * RP->dPFootprint[t][n-first];
			    }
			}

		      if( prf[n-first] > prMax )
			{
			  prMax = prf[n-first];
			  nMaxPos = n-first;
			}
		      if( prr[n-first] > prMax )
			{
			  prMax = prr[n-first];
			  nMaxPos = n-first;
			}
		    }
		}
	      
	      if( IP->is_defined[cl_X] && seed_run > 0 )
		{
		  for( denom = 0, n = 0; n < GRAPH_SIZE; n++ )
		    {
		      prf[n] = ExDistrib( prf[n], B->AN->currTemp );
		      prr[n] = ExDistrib( prr[n], B->AN->currTemp );
		      denom += prf[n] + prr[n];
		    }
		  
		  for( n = 0; n < GRAPH_SIZE; n++ )
		    {
		      prf[n] /= denom;
		      prr[n] /= denom;
		    }
		}

	   /*	      fprintf( stdout, "t = %d seq = %d k = %d prMax = %f nMaxPos = %d\n", 
		       t, seq, k, prMax, nMaxPos ); */	    
	      
	      if( seq == 0 )
		seq = seq;
	    }
#endif

	  for( denom = 0.0, t = 0; t < IP->nNumMotifTypes; t++ )
	    {
	      if( k == kmax || type == -1 )
		{
		  for( n = lastpos; n >= first; n-- )
		    {
		      if( PossStartPos( Pos, n, t, B ) )
			{
			  if(   RP->bUsePosMatrix )
			    {
			      prob = (RP->dProb[k][t][FORWARD][n - first] + 
				      RP->dProb[k][t][REVERSE][n - first]) * 
				RP->dPFootprint[t][n-first] * RP->dPosMatrix[k-1][t];
			    }
			  else
			    {
			      prob = (RP->dProb[k][t][FORWARD][n - first] + 
				      RP->dProb[k][t][REVERSE][n - first]) * 
				RP->dPFootprint[t][n-first] * RP->dEndSiteProb[t] ;
			    }
			  if( IP->is_defined[cl_X] && seed_run > 0 )
			    prob =  ExDistrib( prob, B->AN->currTemp );  
			  denom += prob;
			}
		    }
		}
	      else
		{
		  for( n = lastpos; n >= first; n-- )
		    {
		      if( PossStartPos( Pos, n, t, B ) )
			{
			  nWidth = MotifWidth( B, t );
			  if(   RP->bUsePosMatrix )
			    {
			      prob= RP->dPGap[t][motif_type][lastpos - n - nWidth + 1] *
				 RP->dPosMatrix[k-1][t] * 
				(RP->dProb[k][t][FORWARD][n - first] +
				 RP->dProb[k][t][REVERSE][n - first] ) *
				RP->dPFootprint[t][n-first];
			    }
			  else
			    {
			      prob= RP->dPGap[t][motif_type][lastpos - n - nWidth + 1] *
				(RP->dProb[k][t][FORWARD][n - first] *
				 RP->dPTrans[motif_type][type][t][FORWARD] +
				 RP->dProb[k][t][REVERSE][n - first] *
				 RP->dPTrans[motif_type][type][t][REVERSE]) * 
				RP->dPFootprint[t][n-first];
			    }
			  if( IP->is_defined[cl_X] && seed_run > 0 )
			    prob =  ExDistrib( prob, B->AN->currTemp );  
			  denom += prob;
			}
		    }
		}
	    }

	  index = -1;
	  rand_num = drand()  * denom;
	  dtot = 0.0;
	  for( n = lastpos; n >= first && index == -1; n-- )
	    {
	      cnt = cnt %  B->IP->nNumMotifTypes;   
	      for( tcnt = cnt; tcnt < cnt + B->IP->nNumMotifTypes; tcnt++ )
		{
		  t = tcnt %  B->IP->nNumMotifTypes;
		  nWidth = MotifWidth( B, t );
		  if( PossStartPos( Pos, n, t, B ) )
		    {
		      if( k == kmax || type == -1 )
			{
			  if( RP->bUsePosMatrix )
			    {
			      prob = RP->dProb[k][t][FORWARD][n - first] * 
				RP->dPFootprint[t][n-first] * RP->dPosMatrix[k-1][t];
			    }
			  else
			    {
			      prob = RP->dProb[k][t][FORWARD][n - first] * 
				RP->dPFootprint[t][n-first] * RP->dEndSiteProb[t];
			    }
			   if( IP->is_defined[cl_X] && seed_run > 0 )
			    prob =  ExDistrib( prob, B->AN->currTemp );  
			}
		      else
			{
			  if(   RP->bUsePosMatrix )
			    {
			      prob = RP->dPGap[t][motif_type][lastpos - n - nWidth + 1] *
				RP->dPosMatrix[k-1][t] *
				RP->dProb[k][t][FORWARD][n - first] * RP->dPFootprint[t][n-first];
			    }
			  else
			    {			     
			      prob = RP->dPGap[t][motif_type][lastpos - n - nWidth + 1] *
				RP->dPTrans[motif_type][type][t][FORWARD] *
				RP->dProb[k][t][FORWARD][n - first] * RP->dPFootprint[t][n-first];
			    }
			  if( IP->is_defined[cl_X] && seed_run > 0 )
			    prob =  ExDistrib( prob, B->AN->currTemp );
			}

		      dtot += prob;
		      if( dtot >= rand_num )
			{
			  index = n;
			  type = FORWARD;
			  motif_type = t;
			  break;
			}
			  
		      if( k == kmax || type == -1 )
			{
			  if(   RP->bUsePosMatrix )
			    {
			      prob = RP->dProb[k][t][REVERSE][n - first] * 
				RP->dPFootprint[t][n-first] * RP->dPosMatrix[k-1][t];
			    }
			  else			    
			    {
			      prob = RP->dProb[k][t][REVERSE][n - first] * 
				RP->dPFootprint[t][n-first] * RP->dEndSiteProb[t];
			    }
			  if( IP->is_defined[cl_X] && seed_run > 0 )
			    prob =  ExDistrib( prob, B->AN->currTemp ); 
			}
		      else
			{
			  if(   RP->bUsePosMatrix )
			    {
			      prob = RP->dPGap[t][motif_type][lastpos - n - nWidth + 1] *
				RP->dPosMatrix[k-1][t] *
				RP->dProb[k][t][REVERSE][n - first] * RP->dPFootprint[t][n-first];
			    }
			  else
			    {
			      prob = RP->dPGap[t][motif_type][lastpos - n - nWidth + 1] *
				RP->dPTrans[motif_type][type][t][REVERSE] *
				RP->dProb[k][t][REVERSE][n - first] * RP->dPFootprint[t][n-first];
			    }
			  if( IP->is_defined[cl_X] && seed_run > 0 )
			    prob =  ExDistrib( prob, B->AN->currTemp );  
			}

		      dtot += prob;
		      if( dtot >= rand_num )
			{
			  index = n;
			  type = REVERSE;
			  motif_type = t;
			  break;
			}
		    }
		}
	    }
	  
	  if( index != -1 && dtot > 0.0)
	    {
	      if( IP->is_defined[cl_Q] && bAccumulateCounts )     /* BT 3/27/98 */
		{
		  IP->nAlignCnts[seed_run][motif_type][seq][index - first]++;
		  if( IP->is_defined[cl_X] )		    
		    fprintf( IP->Datafiles->sites_fpt, "%d %d %d %d %d %d %d %f\n",
			     seed_run, iter, motif_type, seq, index - first, nOffset, 
			     MotifWidth( B, motif_type ), B->AN->currTemp );
		  else
		    fprintf( IP->Datafiles->sites_fpt, "%d %d %d %d %d %d %d\n",
			     seed_run, iter, motif_type, seq, index - first, nOffset,
			     MotifWidth( B, motif_type ) );
		  fflush( IP->Datafiles->sites_fpt );
		}
	      
	      if( k == kmax )
		alignSum += B->AP->dAlignCount[kmax][seq];
	      
	      (*occurence)[index][motif_type][0]++;                /* BT 8/26/98 */
	      if(type == FORWARD) (*occurence)[index][motif_type][1]++;
	      else                (*occurence)[index][motif_type][2]++;
	      
	      adjust_counts(B, ADD,index, motif_type, type);
	      B->IP->nNumMotifs[motif_type][type]++;
	      add_motif(B->IP, B->Seq, index, M, motif_type, type);
	      set_in_motif(Pos, index, B, motif_type, type);
	      RP->nSites[seq]++;
	      
	      if( IsPhyloSeq( B, seq ) && ! IsPhyloSeqBySpecies( B, seq ) )
		{
		  for( p = 1; p < SpeciesInc( B, seq ); p++ )
		    {
		      (*occurence)[index+ p * nOffset][motif_type][0]++;                /* BT 8/26/98 */
		      if(type == FORWARD) (*occurence)[index+ p * nOffset][motif_type][1]++;
		      else                (*occurence)[index+ p * nOffset][motif_type][2]++;
		  
		      adjust_counts(B, ADD,index + p * nOffset, motif_type, type);
		      B->IP->nNumMotifs[motif_type][type]++;
		      add_motif(B->IP, B->Seq, index+ p * nOffset, M, motif_type, type);
		      set_in_motif(Pos, index+ p * nOffset, B, motif_type, type);
		      RP->nSites[seq+p]++;		  
		      
		      if( IP->is_defined[cl_Q] && bAccumulateCounts )     /* BT 3/27/98 */
			{
			  IP->nAlignCnts[seed_run][motif_type][seq+p][index - first]++;
			  if( IP->is_defined[cl_X] )		    
			    fprintf( IP->Datafiles->sites_fpt, "%d %d %d %d %d %d %d %f\n",
				     seed_run, iter, motif_type, seq + p, index - first, nOffset, 
				     MotifWidth( B, motif_type ), B->AN->currTemp );
			  else
			    fprintf( IP->Datafiles->sites_fpt, "%d %d %d %d %d %d %d\n",
				     seed_run, iter, motif_type, seq + p, index - first,  
				     nOffset , MotifWidth( B, motif_type ) );
			  fflush( IP->Datafiles->sites_fpt );
			}
		    }
		}
	      
	      if( (iter > WAIT_PERIOD) || bAccumulateCounts )
		update_posterior_prob(B);   /* BT 9/31/98 */
	      
	      lastpos = index - 1;
	    }
	  cnt++;       
	}
      
      psum += RP->prSum[kmax];
	 	 
      FreeRProbArray( B );

      if( IP->is_defined[cl_V] )
	{
	  first = SequenceStartPos( B, seq );
	  nOffset = SequenceLength( B, seq );

	  for( t = 0; t < IP->nNumMotifTypes; t++ )
	    {
	      fprintf( IP->Datafiles->out_fpt, 
		       "-----------------------------------------------------------\n" );
	      fprintf(IP->Datafiles->out_fpt, "                          MOTIF %c\n\n", (char)(97 + t));
	      DumpMotifPositions( t, B, Pos, IP->Datafiles->out_fpt );
	      fprintf( IP->Datafiles->out_fpt, "%d sites\n",NUMMOTIFS(B->IP->nNumMotifs[t]) );
	    }
	  fprintf( IP->Datafiles->out_fpt, "MAP = %g\n", CalcMapProb(B, B->IP->is_defined[cl_R]) );
	  fprintf( IP->Datafiles->out_fpt, "MAP for seq %d = %g\n", seq, CalcSeqMapProb( B, seq, motifMap ) );
	      
	  nSitesSelected = 0;
	  for (t = 0; t < B->IP->nNumMotifTypes; t++) 
	    {
	      last = SequenceEndPos( B, seq ) - IP->nMotifLen[t] + 1;
	      
	      for( n = first; n <= last; n++ )    /* BT 8/15/97 */
		{		  
		  if( Pos[t][n].nMotifStartPos )
		    {
		      nSitesSelected++;
		      RevComp = Pos[t][n].RevComp;        
		      adjust_counts(B,DELETE,n,t,RevComp);
		      B->IP->nNumMotifs[t][RevComp]--;
		      not_in_motif(Pos,n,B,t);    
		      delete_motif(B, n, M, t);
		      RP->nSites[seq]--;	      
		      if( IsPhyloSeq( B, seq ) )
			{
			  for( p = 1; p < SpeciesInc( B, seq ); p++ )
			    {
			      nSitesSelected++;
			      adjust_counts(B,DELETE,n+ p * nOffset,t,RevComp);
			      B->IP->nNumMotifs[t][RevComp]--;
			      not_in_motif(Pos,n+ p * nOffset,B,t);    
			      delete_motif(B, n+ p * nOffset, M, t);
			      RP->nSites[seq+p]--;
			    }
			}
		    }
		  Pos[t][n] = savePos[t][n];
		}

	      for( n = first; n <= last; n++ )    /* BT 8/15/97 */
		{		  
		  if( Pos[t][n].nMotifStartPos )
		    {
		      RevComp = Pos[t][n].RevComp;        
		      adjust_counts( B, ADD, n, t, RevComp );
		      B->IP->nNumMotifs[t][RevComp]++;
		      set_in_motif( Pos, n, B, t, RevComp );    
		      add_motif( B->IP, B->Seq, n, M, t, RevComp );
		      RP->nSites[seq]++;
		      if( IsPhyloSeq( B, seq ) )
			{
			  for( p = 1; p < SpeciesInc( B, seq ); p++ )
			    {
			      adjust_counts( B, ADD, n + p * nOffset, t, RevComp );
			      B->IP->nNumMotifs[t][RevComp]++;
			      set_in_motif( Pos, n + p * nOffset, B, t, RevComp );    
			      add_motif( B->IP, B->Seq, n + p * nOffset, M, t, RevComp );
			      RP->nSites[seq+p]++;
			    }
			}
		    }
		}
	    }
	  if( IsPhyloSeq( B, seq ) )
	    {	  
	      for( p = 1; p < SpeciesInc( B, seq ); p++ )	      
		fprintf( IP->Datafiles->out_fpt, 
			 "seed: %d iter: %d Seqs: %d, %d sites selected: %d total sites: %d\n\n", 
			 seed_run, iter, seq, seq + p, nSitesSelected, TotalNumMotifs( B ) );
	    }
	  else
	    {
	      fprintf( IP->Datafiles->out_fpt, 
		       "seed: %d iter: %d Seq: %d sites selected: %d total sites: %d\n\n", 
		       seed_run, iter, seq, nSitesSelected, TotalNumMotifs( B ) );
	    }	  
	}
    }

  if( ! IP->is_defined[cl_sample_model] )
    update_prob(P, B, (! IP->is_defined[cl_b]));

  if( IP->is_defined[cl_b] )
    UpdateBackgroundModel( B, Pos, P );

  /*  psum2 = log( psum2 ); */   /* DEBUG */
  psum = log( psum ); 

  if( psum > psumMax && ((iter > WAIT_PERIOD)  || bAccumulateCounts ))
    {
      psumMax = psum;
      maxIter = iter;
    }
      
  if( IP->is_defined[cl_opt] )
    dCurrProb = CalcMapProb(B, B->IP->is_defined[cl_R]);
  *dMap = dCurrProb;
  
  if(  /* iter > WAIT_PERIOD  && */ TotalNumMotifs( B ) > RP->nMinSiteMotifs && dCurrProb > RP->dMinSiteMap )
    {
      if( ! RP->nIncludeAlign && ! B->IP->is_defined[cl_Z] )
	fprintf( stdout, "nIncludeAlign set to TRUE at iter %d\n", iter );
      RP->nIncludeAlign = TRUE;
    }

  dProb = dCurrProb;
  
  if( dProb > dProbMax && ((iter > WAIT_PERIOD)  || bAccumulateCounts) )
    {
      dProbMax = dProb;
      maxIter3 = iter;
      bMaxSet = TRUE;
    }
 
#ifdef _DEBUG_ 
  if( iter <= WAIT_PERIOD || bMaxSet )  
    { 
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  fprintf( stdout, 
		   "-----------------------------------------------------------\n" );
	  fprintf(stdout, "                          MOTIF %c\n\n", (char)(97 + t));
	  DumpMotifPositions( t, B, Pos, stdout );
	  fprintf( stdout, "%d sites\n",NUMMOTIFS(B->IP->nNumMotifs[t]) );
	}
      fflush( stdout );
  }  
#endif 
  
  B->AP->dTotalAlgnCnt = alignSum;

  if( IP->is_defined[cl_V] )
    {
      FREEP( savePos, B->IP->nNumMotifTypes ); 
      free( savePos );
      free( motifMap );
    }

  if( IP->is_defined[cl_Q] && bAccumulateCounts )     /* BT 7/26/05 */
    {
      for( t = 0; t <  B->IP->nNumMotifTypes; t++ )
	{
	  for( j = 0; j < IP->nMotifLen[t]; j++ )
	    {
	      fprintf( IP->Datafiles->prob_fpt, "%d %d %d %d", 
		       seed_run, iter, t, j );
	      for( n = 0; n < IP->nAlphaLen; n++ )
		fprintf( IP->Datafiles->prob_fpt, " %g", P->dvInMotifProb[t][j][n] );
	      if( IP->is_defined[cl_X] )
		fprintf( IP->Datafiles->prob_fpt, " %g\n", B->AN->currTemp );
	      else
		fprintf( IP->Datafiles->prob_fpt, " 0\n" );
	    }

	}
    }

  return dProb;  
}


double SampleSitesFromAllK( Model B, PoSition **Pos, Mlist M, int ****occurence, 
			    ProbStruct *P, int iter, int seed_run, short bAccumulateCounts,
			    double *dMap )
{
  int         seq, first, last;
  int         t;
  int         n;
  int         i;
  int         kmax;
  short       RevComp;                      /* BT 8/15/97 */
  IPtype      IP;
  RPType      RP;
  int         k;
  int         lastpos;
  int         motif_type = 0;
  int         index;
  long double rand_num;
  long double dtot = 0.0;
  int         type = -1;
  long double denom;
  int         cnt = 0;
  int         tcnt;
  double      dCurrProb;
  double      dProb;
  double      prob;
  int         alignSum = 0;
  int         nWidth;
  int         nOffset;
  PoSition    **savePos = 0;
  int         nSitesSelected;
  double      *dMapK;
  long double *dMapProb;
  double      *dSeqMap;
  double      *dSeqMapProb;
  double      cntDenom;
  SiteStruct  ***kSitePos;
  double      *priorK;
  double      *cntSitesPerSeq;
  double      **motifMap;
  int         p;
  double      dProbMax = -DBL_MAX;
  int         bMaxSet = FALSE;
  int         j;

#ifdef _DEBUG_
  double     prMax;
  int        nMaxPos;
  double     prf[GRAPH_SIZE];   /* DEBUG */
  double     prr[GRAPH_SIZE];
  double     fpr[GRAPH_SIZE];
#endif
   
  IP = B->IP;
  RP = B->RP;

  NEW( dMapK, RP->nMaxBlocks + 1, double );
  NEW( dMapProb, RP->nMaxBlocks + 1, long double );
  NEW( dSeqMap, RP->nMaxBlocks + 1, double );
  NEW( dSeqMapProb, RP->nMaxBlocks + 1, double );

  if( IP->is_defined[cl_b] )
    UpdateBackgroundModel( B, Pos, P );

  if( IP->is_defined[cl_V] )
    {
      update_posterior_prob(B);   /* 4/18/03 */
      NEWP( savePos, IP->nNumMotifTypes, PoSition);
      for( t=0; t < IP->nNumMotifTypes; t++ )
	{
	  NEW( savePos[t], IP->nSeqLen, PoSition );
	  for( n = 0; n < IP->nSeqLen; n++ )
	    {
	      savePos[t][n] = Pos[t][n];
	    }
	}      
    }
      
  for(seq = 0; seq < B->IP->nNumSequences; seq += SpeciesInc( B, seq ) ) 
    {      
      if( IP->is_defined[cl_V] && seq >= IP->nVerifySeq )   /* BT 04/30/03 */
	break;

      NEWP( motifMap, RP->nMaxBlocks + 1, double );
      NEW( motifMap[0], IP->nNumMotifTypes, double );

      NEWPP( kSitePos, RP->nMaxBlocks + 1, SiteStruct );
      for( k = 1; k <= RP->nMaxBlocks; k++ )
	{
	  NEW( motifMap[k], IP->nNumMotifTypes, double );

	  NEWP( kSitePos[k], SequenceLength( B, seq ), SiteStruct );
	  for( n = 0; n < SequenceLength( B, seq ); n++ )
	    {
	      NEW( kSitePos[k][n], IP->nNumMotifTypes, SiteStruct );
	    }
	}  

      first = SequenceStartPos( B, seq );
      nOffset = SequenceLength( B, seq );
      
      if( IP->is_defined[cl_sample_model] )
	SampleModel( B, P);

      for (t = 0; t < B->IP->nNumMotifTypes; t++) 
	{
	  last = SequenceEndPos( B, seq ) - IP->nMotifLen[t] + 1;
	      
	  for( n = first; n <= last; n++ )    /* BT 8/15/97 */
	    {
	      if( Pos[t][n].nMotifStartPos )
		{
		  RevComp = Pos[t][n].RevComp;        
		  adjust_counts(B,DELETE,n,t,RevComp);
		  B->IP->nNumMotifs[t][RevComp]--;
		  not_in_motif(Pos,n,B,t);    
		  delete_motif(B, n, M, t);
		  if( IsPhyloSeq( B, seq ) )
		    {
		      for( p = 1; p < SpeciesInc( B, seq ); p++ )
			{
			  adjust_counts(B,DELETE,n+ p * nOffset,t,RevComp);
			  B->IP->nNumMotifs[t][RevComp]--;
			  not_in_motif(Pos,n+ p * nOffset,B,t);    
			  delete_motif(B, n+p * nOffset, M, t);
			}
		    }
		}
	    }
	}

	update_prob(P, B, (! B->IP->is_defined[cl_b])); 

      RP->nSites[seq] = 0;
      if( IsPhyloSeq( B, seq ) )
	{
	  for( p = 1; p < SpeciesInc( B, seq ); p++ )
	    RP->nSites[seq+p] = 0;
	}

      NEW( cntSitesPerSeq,  RP->nMaxBlocks + 1, double );
      NEW( priorK,  RP->nMaxBlocks + 1, double );

      for( i = 0; i < IP->nNumSequences; i++ )
	{
	  if( IP->is_defined[cl_V] )
	    {
	      if( i >= IP->nVerifySeq )
		cntSitesPerSeq[RP->nSites[i]]++;
	    }
	  else
	    cntSitesPerSeq[RP->nSites[i]]++;	    
	}
      
      for( cntDenom = 0, k = 0; k <= RP->nMaxBlocks; k++ )
	{
	  cntDenom += (cntSitesPerSeq[k] + RP->priorSitesPerSeq[k]);
	}

      for( k = 0; k <= RP->nMaxBlocks; k++ )
	{
	  priorK[k] = (cntSitesPerSeq[k] + RP->priorSitesPerSeq[k]) / cntDenom;
	}

      if( RP->nMinBlocks == 0 )
	{
	  dMapK[0] = CalcMapProb( B, B->IP->is_defined[cl_R]);
	  dMapProb[0] = expl( dMapK[0] ) * priorK[0];
	  dSeqMap[0] = CalcSeqMapProb( B, seq, motifMap[0] );
	  dSeqMapProb[0] = exp( dSeqMap[0] ) * priorK[0];
	}

      InitRProbArray( B, seq );     /* BT 12/29/97 */
      CalcRProbArray( B,  Pos, seq, P, iter );
  
      if( IP->is_defined[cl_V] )
	SampleMaxBlocks( B, seq, iter, seed_run, bAccumulateCounts );

      for( kmax = RP->nMinBlocks; kmax <= RP->nMaxBlocks; kmax++ )
	{
	  lastpos = SequenceEndPos( B, seq );      
	  for( k = kmax; k > 0; k-- )
	    {
#ifdef _DEBUG_
	      for( t = 0; t < IP->nNumMotifTypes; t++ )
		{
		  memset( prf, 0, sizeof( prf ) );
		  memset( prr, 0, sizeof( prr ) );
		  memset( fpr, 0, sizeof( fpr ) );
		  prMax = -DBL_MAX;
		  for( n = first; n < lastpos; n++ ) 
		    {		  
		      if( n - first >= GRAPH_SIZE )
			break;
		      
		      fpr[n-first] = RP->dPFootprint[t][n-first];
		      
		      if( PossStartPos( Pos, n, t, B ) )
			{
			  if( k == kmax || type == -1 )
			    {
			      if(  RP->bUsePosMatrix )
				{
				  prf[n-first] = RP->dProb[k][t][FORWARD][n - first] * 
				    RP->dPFootprint[t][n-first] * RP->dPosMatrix[k-1][t];
				  prr[n-first] = RP->dProb[k][t][REVERSE][n - first] * 
				    RP->dPFootprint[t][n-first] * RP->dPosMatrix[k-1][t];
				}
			      else
				{
				  prf[n-first] = RP->dProb[k][t][FORWARD][n - first] * 
				    RP->dPFootprint[t][n-first] * RP->dEndSiteProb[t];
				  prr[n-first] = RP->dProb[k][t][REVERSE][n - first] * 
				    RP->dPFootprint[t][n-first] * RP->dEndSiteProb[t];
				}
			    }
			  else
			    {
			      nWidth = MotifWidth( B, t );
			      if( RP->bUsePosMatrix )
				{
				  prf[n-first] =
				    RP->dPGap[t][motif_type][lastpos-n-nWidth + 1] *
				    RP->dPosMatrix[k-1][t] *
				    RP->dProb[k][t][FORWARD][n - first] * RP->dPFootprint[t][n-first];
				  prr[n-first] = 
				    RP->dPGap[t][motif_type][lastpos-n-nWidth + 1] *
				    RP->dPosMatrix[k-1][t] *
				    RP->dProb[k][t][REVERSE][n - first] * RP->dPFootprint[t][n-first];
				}
			      else
				{
				  prf[n-first] =
				    RP->dPGap[t][motif_type][lastpos-n-nWidth + 1] *
				    RP->dPTrans[motif_type][type][t][FORWARD] *
				    RP->dProb[k][t][FORWARD][n - first] * RP->dPFootprint[t][n-first];
				  prr[n-first] = 
				    RP->dPGap[t][motif_type][lastpos-n-nWidth + 1] *
				    RP->dPTrans[motif_type][type][t][REVERSE] *
				    RP->dProb[k][t][REVERSE][n - first] * RP->dPFootprint[t][n-first];
				  if( k == 1 )
				    {
				      prf[n-first] *= RP->dBeginSiteProb[t];
				      prr[n-first] *= RP->dBeginSiteProb[t];
				    }   /* BT 04/16/03 */
				}
			    }
			  
			  if( prf[n-first] > prMax )
			    {
			      prMax = prf[n-first];
			      nMaxPos = n-first;
			    }
			  if( prr[n-first] > prMax )
			    {
			      prMax = prr[n-first];
			      nMaxPos = n-first;
			    }
			}
		    }
		  
		  if( IP->is_defined[cl_X] && seed_run > 0 )
		    {
		      for( denom = 0, n = 0; n < GRAPH_SIZE; n++ )
			{
			  prf[n] = ExDistrib( prf[n], B->AN->currTemp );
			  prr[n] = ExDistrib( prr[n], B->AN->currTemp );
			  denom += prf[n] + prr[n];
			}
		      
		      for( n = 0; n < GRAPH_SIZE; n++ )
			{
			  prf[n] /= denom;
			  prr[n] /= denom;
			}
		    }
		  if( seq == 0 )
		    seq = seq;
		}
#endif

	      for( denom = 0.0, t = 0; t < IP->nNumMotifTypes; t++ )
		{
		  if( k == kmax || type == -1 )
		    {
		      for( n = lastpos; n >= first; n-- )
			{
			  if( PossStartPos( Pos, n, t, B ) )
			    {
			      if( RP->bUsePosMatrix )
				{
				  prob = (RP->dProb[k][t][FORWARD][n - first] + 
					  RP->dProb[k][t][REVERSE][n - first]) * 
				    RP->dPFootprint[t][n-first] * RP->dPosMatrix[k-1][t];
				}
			      else
				{
				  prob = (RP->dProb[k][t][FORWARD][n - first] + 
					  RP->dProb[k][t][REVERSE][n - first]) * 
				    RP->dPFootprint[t][n-first] * RP->dEndSiteProb[t] ;
				}
			      if( IP->is_defined[cl_X] && seed_run > 0 )
				prob =  ExDistrib( prob, B->AN->currTemp );  
			      denom += prob;
			    }
			}
		    }
		  else
		    {
		      for( n = lastpos; n >= first; n-- )
			{
			  if( PossStartPos( Pos, n, t, B ) )
			    {
			      nWidth = MotifWidth( B, t );
			      if( RP->bUsePosMatrix )
				{
				  prob= RP->dPGap[t][motif_type][lastpos - n - nWidth + 1] *
				    RP->dPosMatrix[k-1][t] * 
				    (RP->dProb[k][t][FORWARD][n - first] +
				     RP->dProb[k][t][REVERSE][n - first] ) *
				    RP->dPFootprint[t][n-first];
				}
			      else
				{
				  prob= RP->dPGap[t][motif_type][lastpos - n - nWidth + 1] *
				    (RP->dProb[k][t][FORWARD][n - first] *
				     RP->dPTrans[motif_type][type][t][FORWARD] +
				     RP->dProb[k][t][REVERSE][n - first] *
				     RP->dPTrans[motif_type][type][t][REVERSE]) * 
				    RP->dPFootprint[t][n-first];
				  if( k == 1 )
				    prob *= RP->dBeginSiteProb[t];
				}
			      if( IP->is_defined[cl_X] && seed_run > 0 )
				prob =  ExDistrib( prob, B->AN->currTemp );  
			      denom += prob;
			    }
			}
		    }
		}

	      index = -1;
	      rand_num = drand()  * denom;
	      dtot = 0.0;
	      for( n = lastpos; n >= first && index == -1; n-- )
		{
		  cnt = cnt %  B->IP->nNumMotifTypes;   
		  for( tcnt = cnt; tcnt < cnt + B->IP->nNumMotifTypes; tcnt++ )
		    {
		      t = tcnt %  B->IP->nNumMotifTypes;
		      nWidth = MotifWidth( B, t );
		      if( PossStartPos( Pos, n, t, B ) )
			{
			  if( k == kmax || type == -1 )
			    {
			      if( RP->bUsePosMatrix )
				{
				  prob = RP->dProb[k][t][FORWARD][n - first] * 
				    RP->dPFootprint[t][n-first] * RP->dPosMatrix[k-1][t];
				}
			      else
				{
				  prob = RP->dProb[k][t][FORWARD][n - first] * 
				    RP->dPFootprint[t][n-first] * RP->dEndSiteProb[t];
				}
			      if( IP->is_defined[cl_X] && seed_run > 0 )
				prob =  ExDistrib( prob, B->AN->currTemp );  
			    }
			  else
			    {
			      if( RP->bUsePosMatrix )
				{
				  prob = RP->dPGap[t][motif_type][lastpos - n - nWidth + 1] *
				    RP->dPosMatrix[k-1][t] *
				    RP->dProb[k][t][FORWARD][n - first] * RP->dPFootprint[t][n-first];
				}
			      else
				{			     
				  prob = RP->dPGap[t][motif_type][lastpos - n - nWidth + 1] *
				    RP->dPTrans[motif_type][type][t][FORWARD] *
				    RP->dProb[k][t][FORWARD][n - first] * RP->dPFootprint[t][n-first];
				  if( k == 1 )
				    prob *= RP->dBeginSiteProb[t];
				}
			      if( IP->is_defined[cl_X] && seed_run > 0 )
				prob =  ExDistrib( prob, B->AN->currTemp );
			    }
			  
			  dtot += prob;
			  if( dtot >= rand_num )
			    {
			      index = n;
			      type = FORWARD;
			      motif_type = t;
			      break;
			    }
			  
			  if( k == kmax || type == -1 )
			    {
			      if(   RP->bUsePosMatrix )
				{
				  prob = RP->dProb[k][t][REVERSE][n - first] * 
				    RP->dPFootprint[t][n-first] * RP->dPosMatrix[k-1][t];
				}
			      else			    
				{
				  prob = RP->dProb[k][t][REVERSE][n - first] * 
				    RP->dPFootprint[t][n-first] * RP->dEndSiteProb[t];
				}
			      if( IP->is_defined[cl_X] && seed_run > 0 )
				prob =  ExDistrib( prob, B->AN->currTemp ); 
			    }
			  else
			    {
			      if( RP->bUsePosMatrix )
				{
				  prob = RP->dPGap[t][motif_type][lastpos - n - nWidth + 1] *
				    RP->dPosMatrix[k-1][t] *
				    RP->dProb[k][t][REVERSE][n - first] * RP->dPFootprint[t][n-first];
				}
			      else
				{
				  prob = RP->dPGap[t][motif_type][lastpos - n - nWidth + 1] *
				    RP->dPTrans[motif_type][type][t][REVERSE] *
				    RP->dProb[k][t][REVERSE][n - first] * RP->dPFootprint[t][n-first];
				  if( k == 1 )
				    prob *= RP->dBeginSiteProb[t];
				}
			      if( IP->is_defined[cl_X] && seed_run > 0 )
				prob =  ExDistrib( prob, B->AN->currTemp );  
			    }
			  
			  dtot += prob;
			  if( dtot >= rand_num )
			    {
			      index = n;
			      type = REVERSE;
			      motif_type = t;
			      break;
			    }
			}
		    }
		}
	  
	      if( index != -1 && dtot > 0.0)
		{
		  if( k == kmax )
		    alignSum += B->AP->dAlignCount[kmax][seq];
		  		  
		  adjust_counts(B, ADD,index, motif_type, type);
		  B->IP->nNumMotifs[motif_type][type]++;
		  add_motif(B->IP, B->Seq, index, M, motif_type, type);
		  set_in_motif(Pos, index, B, motif_type, type);
		  RP->nSites[seq]++;

		  kSitePos[kmax][index-first][motif_type].nMotifStart = TRUE;
		  kSitePos[kmax][index-first][motif_type].nRevComp = type;
		  		  
		  if( IsPhyloSeq( B, seq ) )
		    {
		      for( p = 1; p < SpeciesInc( B, seq ); p++ )
			{
			  adjust_counts(B, ADD,index + p * nOffset, motif_type, type);
			  B->IP->nNumMotifs[motif_type][type]++;
			  add_motif(B->IP, B->Seq, index+ p * nOffset, M, motif_type, type);
			  set_in_motif(Pos, index+ p * nOffset, B, motif_type, type);
			  RP->nSites[seq+p]++;		  
			}
		    }
		  
		  if( (iter > WAIT_PERIOD) || bAccumulateCounts )
		    update_posterior_prob(B);   /* BT 9/31/98 */
		  
		  lastpos = index - 1;
		} 	
	      cnt++;       
	    }

	  first = SequenceStartPos( B, seq );
	  nOffset = SequenceLength( B, seq );
	  
	  if( dtot > 0.0 || kmax == 0 )
	    {
	      dMapK[kmax]  = CalcMapProb(B, B->IP->is_defined[cl_R]);
	      dMapProb[kmax] = expl( dMapK[kmax] ) * priorK[kmax];
	      dSeqMap[kmax] = CalcSeqMapProb( B, seq, motifMap[kmax] );
	      dSeqMapProb[kmax] = exp( dSeqMap[kmax] ) * priorK[kmax];
	    }
	  else
	    {
	      dMapK[kmax]  = 0;
	      dMapProb[kmax] = 0;
	      dSeqMap[kmax] = 0;
	      dSeqMapProb[kmax] = 0;
	    }
	  
	  /*
#ifdef _DEBUG_
	  update_posterior_prob(B);  
	  for( motifSum = 0, fragSum = 0, t = 0; t < IP->nNumMotifTypes; t++ )
	    {
	      fprintf( stdout, 
		       "-----------------------------------------------------------\n" );
	      fprintf( stdout, "                          MOTIF %c\n\n", (char)(97 + t));
	      DumpMotifPositions( t, B, Pos, P, stdout ); 
	      fprintf( stdout, "%d sites\n",NUMMOTIFS(B->IP->nNumMotifs[t]) );
	      fprintf( stdout, "motif map = %g frag map = %g\n", 
		       CalcMotifMap(B, t, B->IP->is_defined[cl_R]),
		       CalcMotifFragMap(B, t, B->IP->is_defined[cl_R]) );
	      motifSum += CalcMotifMap(B, t, B->IP->is_defined[cl_R]);
	      fragSum += CalcMotifFragMap(B, t, B->IP->is_defined[cl_R]);
	    }
	  fprintf( stdout, "MAP = %g\n", CalcMapProb(B, B->IP->is_defined[cl_R]) );
	  fprintf( stdout, "m = %g f = %g bkgnd = %g al = %g s = %g null = %g\n",
		   motifSum, 
		   fragSum, 
		   CalcBkgndMap(B, B->IP->is_defined[cl_R]),
		   CalcBetaMap( B, B->IP->is_defined[cl_R]),
		   CalcSitePerSeqMap( B ),
		   IP->dnull_map );

	  fprintf( stdout, "k = %d Probability = %Lg\n", kmax, dMapProb[kmax] );
#endif	  
	  */
	  
	  nSitesSelected = 0;
	  for (t = 0; t < B->IP->nNumMotifTypes; t++) 
	    {
	      last = SequenceEndPos( B, seq ) - IP->nMotifLen[t] + 1;
	      
	      for( n = first; n <= last; n++ )    /* BT 8/15/97 */
		{		  
		  if( Pos[t][n].nMotifStartPos )
		    {
		      nSitesSelected++;
		      RevComp = Pos[t][n].RevComp;        
		      adjust_counts(B,DELETE,n,t,RevComp);
		      B->IP->nNumMotifs[t][RevComp]--;
		      not_in_motif(Pos,n,B,t);    
		      delete_motif(B, n, M, t);
		      RP->nSites[seq]--;	      
		      if( IsPhyloSeq( B, seq ) )
			{
			  for( p = 1; p < SpeciesInc( B, seq ); p++ )
			    {
			      nSitesSelected++;
			      adjust_counts(B,DELETE,n+ p * nOffset,t,RevComp);
			      B->IP->nNumMotifs[t][RevComp]--;
			      not_in_motif(Pos,n+ p * nOffset,B,t);    
			      delete_motif(B, n+ p * nOffset, M, t);
			      RP->nSites[seq+p]--;
			    }
			}
		    }
		}	     
	    }
	  update_prob(P, B, (! B->IP->is_defined[cl_b]));
	}

      for( denom = 0.0, kmax = 0; kmax <= RP->nMaxBlocks; kmax++ )
	{
	  denom += dMapProb[kmax];
	}

      k = -1;
      rand_num = drand() * denom;
      dtot = 0.0;
      for( kmax = RP->nMinBlocks; kmax <= RP->nMaxBlocks; kmax++ )
	{
	  dtot += dMapProb[kmax];
	  if( dtot > rand_num )
	    {
	      k = kmax;
	      break;
	    }
	}
      /*
#ifdef _DEBUG_
      if( k > 0 )
	fprintf( stdout, "seq %d - %d sites selected\n", seq, k );
      else
	fprintf( stdout, "seq %d - 0 sites selected\n", seq );
      for( kmax = RP->nMinBlocks; kmax <= RP->nMaxBlocks; kmax++ )
	{
	  printf( "k = %d MAP = %g prob = %Lg\n", kmax, dMapK[kmax],  dMapProb[kmax] );
	}
      printf( "\n" );
#endif
      */

      if( k > 0 )
	{
	  update_posterior_prob(B);   /* BT 04/06/04 */
	  for (t = 0; t < B->IP->nNumMotifTypes; t++) 
	    {
	      last = SequenceEndPos( B, seq ) - IP->nMotifLen[t] + 1;
	      
	      for( n = first; n <= last; n++ )    /* BT 8/15/97 */
		{		  
		  if( kSitePos[k][n-first][t].nMotifStart )
		    {		      
		      RevComp = kSitePos[k][n-first][t].nRevComp;        

		      (*occurence)[n][t][0]++;               
		      if(RevComp == FORWARD) 
			(*occurence)[n][t][1]++;
		      else                
			(*occurence)[n][t][2]++;

		      adjust_counts( B, ADD, n, t, RevComp );
		      B->IP->nNumMotifs[t][RevComp]++;
		      set_in_motif( Pos, n, B, t, RevComp );    
		      add_motif( B->IP, B->Seq, n, M, t, RevComp );
		      RP->nSites[seq]++;
		      
		      if( IP->is_defined[cl_Q] && bAccumulateCounts )     /* BT 3/27/98 */
			{
			  IP->nAlignCnts[seed_run][t][seq][n - first]++;
			  if( IP->is_defined[cl_X] )		    
			    fprintf( IP->Datafiles->sites_fpt, "%d %d %d %d %d %d %d %f\n",
				     seed_run, iter, t, seq, n - first, nOffset, 
				     MotifWidth( B, t ), B->AN->currTemp );
			  else
			    fprintf( IP->Datafiles->sites_fpt, "%d %d %d %d %d %d %d\n",
				     seed_run, iter, t, seq, n - first, 
				     nOffset, MotifWidth( B, t ) );
			  fflush( IP->Datafiles->sites_fpt );
			}

		      if( IsPhyloSeq( B, seq ) )
			{
			  for( p = 1; p < SpeciesInc( B, seq ); p++ )
			    {
			      (*occurence)[n + p * nOffset][t][0]++;               
			      if(RevComp == FORWARD) 
				(*occurence)[n + p * nOffset][t][1]++;
			      else                
				(*occurence)[n + p * nOffset][t][2]++;
			      
			      adjust_counts( B, ADD, n + p * nOffset, t, RevComp );
			      B->IP->nNumMotifs[t][RevComp]++;
			      set_in_motif( Pos, n + p * nOffset, B, t, RevComp );    
			      add_motif( B->IP, B->Seq, n + p * nOffset, M, t, RevComp );
			      RP->nSites[seq+p]++;
			      
			      if( IP->is_defined[cl_Q] && bAccumulateCounts )     /* BT 3/27/98 */
				{
				  IP->nAlignCnts[seed_run][t][seq+p][n - first]++;
				  if( IP->is_defined[cl_X] )		    
				    fprintf( IP->Datafiles->sites_fpt, "%d %d %d %d %d %d %d %f\n",
					     seed_run, iter, t, seq + p, n - first, nOffset, 
					     MotifWidth( B, t ), B->AN->currTemp );
				  else
				    fprintf( IP->Datafiles->sites_fpt, "%d %d %d %d %d %d %d\n",
					     seed_run, iter, t, seq + p, n - first, 
					     nOffset, MotifWidth( B, t ) );
				  fflush( IP->Datafiles->sites_fpt );
				}
			    }
			}
		    }
		}
	    }
	}
      
      if( IP->is_defined[cl_V] && iter <= IP->nMaxIterations && seed_run > 0 )
	{	  
	  if( (iter > WAIT_PERIOD) || bAccumulateCounts )
	    update_posterior_prob(B);   /* 4/18/03 */

	  first = SequenceStartPos( B, seq );
	  nOffset = SequenceLength( B, seq );

	  for( t = 0; t < IP->nNumMotifTypes; t++ )
	    {
	      fprintf( IP->Datafiles->out_fpt, 
		       "-----------------------------------------------------------\n" );
	      fprintf(IP->Datafiles->out_fpt, "                          MOTIF %c\n\n", (char)(97 + t));
	      DumpMotifPositions( t, B, Pos, IP->Datafiles->out_fpt );
	      fprintf( IP->Datafiles->out_fpt, "%d sites\n",NUMMOTIFS(B->IP->nNumMotifs[t]) );
	    }
	  
	  fprintf( IP->Datafiles->out_fpt, "MAP: " );
	  for( kmax = 0; kmax <= RP->nMaxBlocks; kmax++ )
	    {
	      fprintf( IP->Datafiles->out_fpt, "%g ", dMapK[kmax] );
	    }
	  fprintf( IP->Datafiles->out_fpt, "\n" );

	  fprintf( IP->Datafiles->out_fpt, "dMapProb: " );
	  for( kmax = 0; kmax <= RP->nMaxBlocks; kmax++ )
	    {
	      fprintf( IP->Datafiles->out_fpt, "%Lg ", dMapProb[kmax] );
	    }
	  fprintf( IP->Datafiles->out_fpt, "\n" );

	  fprintf( IP->Datafiles->out_fpt, "SeqMap: " );
	  for( kmax = 0; kmax <= RP->nMaxBlocks; kmax++ )
	    {
	      fprintf( IP->Datafiles->out_fpt, "%g ", dSeqMap[kmax] );
	    }
	  fprintf( IP->Datafiles->out_fpt, "\n" );

	  fprintf( IP->Datafiles->out_fpt, "dSeqMapProb: " );
	  for( kmax = 0; kmax <= RP->nMaxBlocks; kmax++ )
	    {
	      fprintf( IP->Datafiles->out_fpt, "%g ", dSeqMapProb[kmax] );
	    }
	  fprintf( IP->Datafiles->out_fpt, "\n" );

	  fprintf( IP->Datafiles->out_fpt, "P(R|k): 1 " );
	  for( kmax = 1; kmax <= RP->nMaxBlocks; kmax++ )
	    {
	      fprintf( IP->Datafiles->out_fpt, "%g ", RP->prSum[kmax] );
	    }
	  fprintf( IP->Datafiles->out_fpt, "\n" );

	  fprintf( IP->Datafiles->out_fpt, "AP: " );
	  for( kmax = 0; kmax <= RP->nMaxBlocks; kmax++ )
	    {
	      fprintf( IP->Datafiles->out_fpt, "%g ", B->AP->dAlignCount[kmax][seq] );
	    }
	  fprintf( IP->Datafiles->out_fpt, "\n" );

	  fprintf( IP->Datafiles->out_fpt, "P(k): " );
	  for( kmax = 0; kmax <= RP->nMaxBlocks; kmax++ )
	    {
	      fprintf( IP->Datafiles->out_fpt, "%g ", priorK[kmax] );
	    }
	  fprintf( IP->Datafiles->out_fpt, "\n" );
	  	  
	  fprintf( IP->Datafiles->out_fpt, "MAP = %g\n", CalcMapProb(B, B->IP->is_defined[cl_R]) );
	  fprintf( IP->Datafiles->out_fpt, "MAP for seq %d = %g\n", seq, log( dSeqMapProb[k] ) );

	  DumpTrans( B, IP->Datafiles->out_fpt );

	  nSitesSelected = 0;
	  for (t = 0; t < IP->nNumMotifTypes; t++) 
	    {
	      last = SequenceEndPos( B, seq ) - IP->nMotifLen[t] + 1;
	      
	      for( n = first; n <= last; n++ )    /* BT 8/15/97 */
		{		  
		  if( Pos[t][n].nMotifStartPos )
		    {
		      nSitesSelected++;
		      RevComp = Pos[t][n].RevComp;        
		      adjust_counts(B,DELETE,n,t,RevComp);
		      B->IP->nNumMotifs[t][RevComp]--;
		      not_in_motif(Pos,n,B,t);    
		      delete_motif(B, n, M, t);
		      RP->nSites[seq]--;	      
		      if( IsPhyloSeq( B, seq ) )
			{
			  for( p = 1; p < SpeciesInc( B, seq ); p++ )
			    {
			      nSitesSelected++;
			      adjust_counts(B,DELETE,n+ p * nOffset,t,RevComp);
			      B->IP->nNumMotifs[t][RevComp]--;
			      not_in_motif(Pos,n + p * nOffset,B,t);    
			      delete_motif(B, n + p * nOffset, M, t);
			      RP->nSites[seq+p]--;
			    }
			}
		    }
		  Pos[t][n] = savePos[t][n];
		}
	    }

	  for (t = 0; t < IP->nNumMotifTypes; t++) 
	    {
	      last = SequenceEndPos( B, seq ) - IP->nMotifLen[t] + 1;
	      for( n = first; n <= last; n++ )    /* BT 8/15/97 */
		{		  
		  if( Pos[t][n].nMotifStartPos )
		    {
		      RevComp = Pos[t][n].RevComp;        
		      adjust_counts( B, ADD, n, t, RevComp );
		      B->IP->nNumMotifs[t][RevComp]++;
		      set_in_motif( Pos, n, B, t, RevComp );    
		      add_motif( B->IP, B->Seq, n, M, t, RevComp );
		      RP->nSites[seq]++;
		      if( IsPhyloSeq( B, seq ) )
			{
			  for( p = 1; p < SpeciesInc( B, seq ); p++ )
			    {
			      adjust_counts( B, ADD, n + p * nOffset, t, RevComp );
			      B->IP->nNumMotifs[t][RevComp]++;
			      set_in_motif( Pos, n + p * nOffset, B, t, RevComp );    
			      add_motif( B->IP, B->Seq, n + p * nOffset, M, t, RevComp );
			      RP->nSites[seq+p]++;
			    }
			}
		    }
		}
	    }

	  if( IsPhyloSeq( B, seq ) )
	    {	  
	      for( p = 1; p < SpeciesInc( B, seq ); p++ )
		fprintf( IP->Datafiles->out_fpt, 
			 "seed: %d iter: %d Seqs: %d, %d sites selected: %d total sites: %d\n\n", 
			 seed_run, iter, seq, seq + p, nSitesSelected, TotalNumMotifs( B ) );
	    }
	  else
	    {
	      fprintf( IP->Datafiles->out_fpt, 
		       "seed: %d iter: %d Seq: %d sites selected: %d total sites: %d\n\n", 
		       seed_run, iter, seq, nSitesSelected, TotalNumMotifs( B ) );
	    }	  
	}

      FreeRProbArray( B );
      
      for( k = 1; k <= RP->nMaxBlocks; k++ )	
	{
	  FREEP( kSitePos[k], SequenceLength( B, seq ) );
	  free( motifMap[k] );
	}
      free( kSitePos );
      free( motifMap );

      free( cntSitesPerSeq );
      free( priorK );
    }

  update_prob(P, B, (! B->IP->is_defined[cl_b]));

  if( IP->is_defined[cl_b] )
    UpdateBackgroundModel( B, Pos, P );
      
  dCurrProb = CalcMapProb(B, B->IP->is_defined[cl_R]);
  *dMap = dCurrProb;
  
  if( /* iter > WAIT_PERIOD  && */ TotalNumMotifs( B ) > RP->nMinSiteMotifs && dCurrProb > RP->dMinSiteMap )
    {
      if( ! RP->nIncludeAlign && ! B->IP->is_defined[cl_Z] )
	fprintf( stdout, "nIncludeAlign set to TRUE at iter %d\n", iter );
      RP->nIncludeAlign = TRUE;
    }

  dProb = dCurrProb;
  
  if( dProb > dProbMax && ((iter > WAIT_PERIOD)  || bAccumulateCounts) )
    {
      dProbMax = dProb;
      bMaxSet = TRUE;
    }
 
#ifdef _DEBUG_ 
  if( iter <= WAIT_PERIOD || bMaxSet )  
    { 
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  fprintf( stdout, 
		   "-----------------------------------------------------------\n" );
	  fprintf(stdout, "                          MOTIF %c\n\n", (char)(97 + t));
	  DumpMotifPositions( t, B, Pos, stdout );
	  fprintf( stdout, "%d sites\n",NUMMOTIFS(B->IP->nNumMotifs[t]) );
	}
    }

  fprintf( stdout, " %d %d %g %d\n", 
	     seed_run, iter, dCurrProb, TotalNumMotifs( B ) );  
#endif 
  
  B->AP->dTotalAlgnCnt = alignSum;

 
  if( IP->is_defined[cl_Q] && bAccumulateCounts )     /* BT 7/26/05 */
    {
      for( t = 0; t <  B->IP->nNumMotifTypes; t++ )
	{
	  for( j = 0; j < IP->nMotifLen[t]; j++ )
	    {
	      fprintf( IP->Datafiles->prob_fpt, "%d %d %d %d", 
		       seed_run, iter, t, j );
	      for( n = 0; n < IP->nAlphaLen; n++ )
		fprintf( IP->Datafiles->prob_fpt, " %g", P->dvInMotifProb[t][j][n] );
	      if( IP->is_defined[cl_X] )
		fprintf( IP->Datafiles->prob_fpt, " %g\n", B->AN->currTemp );
	      else
		fprintf( IP->Datafiles->prob_fpt, " 0\n" );
	    }

	}
    }

 if( IP->is_defined[cl_V] )
    {
      FREEP( savePos, IP->nNumMotifTypes ); 
      free( savePos );
    }

  free( dMapK );
  free( dMapProb );
  free( dSeqMap );
  free( dSeqMapProb );

  return dProb;  
}


short PossStartPos( PoSition **Pos, int n, int motif_type, Model B)
{
  int      t;
  int      newpos;
  PoSition *Pos_t;  /* Pos_t=Pos[t]  */
  int      width;

  if( ! PossFragStartPos( Pos, n, motif_type, B) )
    return FALSE;

  width = MotifWidth( B, motif_type );

  for( t=0; t < B->IP->nNumMotifTypes; t++)
    {
      Pos_t = Pos[t];
      for( newpos=n; newpos < n + width; newpos++)
	{
	  if( Pos_t[newpos].nInMotif )
	    return FALSE;
	  if( Pos_t[newpos].nMotifStartPos )
	    return FALSE;
	}
    }

  return TRUE;
}

int SampleMaxBlocks( Model B, int seq, int iter, int seed_run, short bAccumulateCounts )
{
  if( B->IP->is_defined[cl_hm] && B->RP->useProbK )
    return SampleMaxBlocksHier( B, seq, iter, seed_run, bAccumulateCounts );
  else
    return SampleMaxBlocksK( B, seq, iter, seed_run, bAccumulateCounts );
}


int SampleMaxBlocksHier( Model B, int seq, int iter, int seed_run, short bAccumulateCounts )
{
  IPtype  IP;
  RPType  RP;
  ALType  AP;
  int     k;
  double  rnd;
  double  dtot = 0.0;
  double  dtot2 = 0.0;
  double* prob;
  int     blocks = 0;
  double  prob0;

  IP = B->IP;
  RP = B->RP;
  AP = B->AP;

  NEW( prob,  RP->nMaxBlocks + 1, double );

  prob0 = RP->probK[0];
  prob[0] = 1;
  dtot = prob[0];

  for( k = 1; k <= RP->nMaxBlocks; k++)
    {
      prob[k] = RP->prSum[k] * RP->probK[k] / prob0;
      prob[k] *= (AP->dAlignWt[k] / AP->dAlignCount[k][seq]);  
      dtot += prob[k]; 
    }  
  
  rnd = drand() * dtot;  /* BT 4/20/98 */
  dtot2 = 0.0;
  blocks = 0;
  for( k = RP->nMaxBlocks; k > 0; k--)
    {
      dtot2 += prob[k];
      if( dtot2 >= rnd )
	{
	  blocks = k;
	  break;
	}
    }

  if( IP->is_defined[cl_Q] && bAccumulateCounts )
    {
      fprintf( IP->Datafiles->dist_fpt, "%d %d %d %d",
	       seed_run, iter, seq, blocks );
      for( k = 0; k <= RP->nMaxBlocks; k++)
	{
	  fprintf( IP->Datafiles->dist_fpt, " %g", prob[k] / dtot );
	}
      fprintf( IP->Datafiles->dist_fpt, "\n" );      
    }
  
  free( prob );

  return blocks;
}


int SampleMaxBlocksK( Model B, int seq, int iter, int seed_run, short bAccumulateCounts )
{
  IPtype  IP;
  RPType  RP;
  ALType  AP;
  int     t;
  int     k;
  int     i;
  double  rnd;
  double  dtot = 0.0;
  double  dtot2 = 0.0;
  double* prob;
  int     blocks = 0;
  int     nNumMotifs;
  double  prob0;
  double  pr;
  double  *cntSitesPerSeq = 0;
  double  denom = 0;
  double  weight;

  IP = B->IP;
  RP = B->RP;
  AP = B->AP;

  weight = IP->dPseudoSiteWt / (1.0 - IP->dPseudoSiteWt);  /* BT 11/23/99 */

  NEW( prob,  RP->nMaxBlocks + 1, double );

  if( RP->nUseFixedBlockProbs )
    {
      NEW( cntSitesPerSeq,  RP->nMaxBlocks + 1, double );
      
      for( i = 0; i < IP->nNumSequences; i++ )
	{
	  if( IP->is_defined[cl_V] )
	    {
	      if( i >= IP->nVerifySeq )
		cntSitesPerSeq[RP->nSites[i]]++;
	    }
	  else
	    cntSitesPerSeq[RP->nSites[i]] += GetSeqWeight( B, i, 0 );	    
	}
      for( denom = 0, k = 0; k <= RP->nMaxBlocks; k++ )
	{
	  denom += (cntSitesPerSeq[k] + RP->priorSitesPerSeq[k]);
	}
    }
  
  for( nNumMotifs = 0, t = 0; t < IP->nNumMotifTypes; t++ )
    nNumMotifs += NUMMOTIFS(IP->nNumMotifs[t]);
  
  if( RP->nUseFixedBlockProbs )
    {
      if( ! RP->nScaleZeroBlocksProb ) 
	{
	  prob0 =  (cntSitesPerSeq[0] + RP->priorSitesPerSeq[0]) / denom;
	  prob[0] = 1;
	}
      else
	{
	  prob0 =  max( RP->dExpBlkCnt[0], MIN_SITE_PROB ); 
	  /* if( IP->is_defined[cl_bayes] && iter < IP->burnInPeriod )
	    prob[0] = RP->dExpBlkCnt[0];
	    else */ if( RP->nDontAdjustPrior )
	    {
	      prob[0] = 1;
	    }
	  else
	    prob[0] = RP->dExpBlkCnt[0];
	}
    }
  else
    {
      prob0 = NegBinomialPrior( weight * RP->dPriorSites + nNumMotifs, 
				weight * RP->dPriorSeq + IP->nNumSequences - 1, 0 ); 
      prob[0] = 1;
    }
  dtot = prob[0];
  
  for( k = 1; k <= RP->nMaxBlocks; k++)
    {
      if( RP->nUseFixedBlockProbs )
	{
	  /* if( IP->is_defined[cl_bayes] && iter < IP->burnInPeriod )
	    prob[k] = RP->dExpBlkCnt[k];
	    else */ if( RP->nScaleZeroBlocksProb ) 
	    {
	      if( RP->nDontAdjustPrior )
		{
		  pr = (RP->prSum[k] * RP->dExpBlkCnt[k]) / AP->dAlignCount[k][seq];
		  prob[k] = pr / prob0;		  
		}
	      else
		prob[k] = RP->dExpBlkCnt[k]; 	 
	    }
	  else
	    {
	      if( AP->dAlignCount[k][seq] > 0 )
		{
		  pr = (cntSitesPerSeq[k] + RP->priorSitesPerSeq[k]) / denom;
		  pr *= RP->prSum[k]; 
		  prob[k] = pr / prob0;
		  if( RP->nIncludeAlign )
		    {
		      prob[k] *= (AP->dAlignWt[k] / AP->dAlignCount[k][seq]);  
		    }
		}
	    }
	}
      else
	{
	  pr = NegBinomialPrior( weight * RP->dPriorSites + nNumMotifs, 
				 weight * RP->dPriorSeq + IP->nNumSequences - 1, k );
	  prob[k] = RP->prSum[k] * pr / prob0;
	  if( RP->nIncludeAlign )
	    {
	      if( AP->dAlignCount[k][seq] > 0 )
		prob[k] /= AP->dAlignCount[k][seq];    /* TEST */
	    }
	}
      
      dtot += prob[k]; 
    }  

  rnd = drand() * dtot;  /* BT 4/20/98 */
  dtot2 = 0.0;
  blocks = 0;
  for( k = RP->nMaxBlocks; k > 0; k--)
    {
      dtot2 += prob[k];
      if( dtot2 >= rnd )
	{
	  blocks = k;
	  break;
	}
    }

  if( IP->is_defined[cl_bayes] && seed_run > 0 )
    RP->kCounts[seq][blocks]++;
  if( IP->is_defined[cl_V] && iter <= IP->nMaxIterations && seed_run > 0 )
    {
      fprintf( IP->Datafiles->out_fpt, "seed: %d iter: %d Seq: %d P(k|R): ", 
	       seed_run, iter, seq );
      for( k = 0; k <= RP->nMaxBlocks; k++)
	{
	  fprintf( IP->Datafiles->out_fpt, "%g ", prob[k] / dtot );
	}
      fprintf( IP->Datafiles->out_fpt, "\n" );
    }

  if( IP->is_defined[cl_Q] && bAccumulateCounts )
    {
      fprintf( IP->Datafiles->dist_fpt, "%d %d %d %d",
	       seed_run, iter, seq, blocks );
      for( k = 0; k <= RP->nMaxBlocks; k++)
	{
	  fprintf( IP->Datafiles->dist_fpt, " %g", prob[k] / dtot );
	}
       for( k = 0; k <= RP->nMaxBlocks; k++)
	{
	  fprintf( IP->Datafiles->dist_fpt, " %g", RP->prSum[k] );
	}
      fprintf( IP->Datafiles->dist_fpt, "\n" );      
    }
  
  free( prob );

  if( RP->nUseFixedBlockProbs )
    {
      free( cntSitesPerSeq );
    }

  return blocks;
}


double CalcAlignmentProb( Model B, PoSition **Pos, Mlist M, ProbStruct *P, double *prob )
{
  int        i;
  int        t;
  double     result;
  double     dMotifProb;
  double     dBGProb;
  int        seq;
  
  result = CalcMapProb( B, B->IP->is_defined[cl_R] ); 
      
  for( *prob = 0, t=0; t<B->IP->nNumMotifTypes; t++)
    {
      for( i = 0; i < B->IP->nSeqLen; i++)
	{ 
	   if(Pos[t][i].nMotifStartPos) 
	     {
	       seq = SequenceFromPosition( B, i );
	       if( IsPhyloSeq( B, seq ) )
		 {
		   if( seq % SpeciesInc( B, seq ) == 0 )
		     in_motif_prob(B, *P, i, t, Pos[t][i].RevComp, FALSE, 
				   &dMotifProb, &dBGProb);
		 }
	       else
		 in_motif_prob(B, *P, i, t, Pos[t][i].RevComp, FALSE, 
			       &dMotifProb, &dBGProb);

	       *prob += log( dMotifProb / dBGProb );  /* BT 2/22/99 */
	       /*  *prob += log( (dMotifProb / dBGProb) * B->IP->dposterior_prob[t] ); */
	       
	       i += (B->IP->nMotifLen[t]-1);
	     }
	}
    }
  
  return result;
}


void UpdateBackgroundModel( Model B, PoSition **Pos, ProbStruct *P )
{
  Ctype      C;
  IPtype     IP;
  RPType     RP;
  int        **bgCounts;
  int        **bgPos;
  int        t;
  int        n;
  int        k;
  int        i;
  int        index;
  int        nSites;
  int        nSiteCount;
  short      bOK;
  int        nTrys = 0;
  double     dPseudo =0;
  int        totalMotifs;
  int        nMotifType;
  double     dBGDenom;
  char       *cSeq = B->Seq->R[0];
   
  C = B->C;
  IP = B->IP;
  RP = B->RP;

  NEWP( bgCounts, IP->nAlphaLen, int );
  for( n = 0; n < IP->nAlphaLen; n++ )
    {
      NEW( bgCounts[n],  IP->nNumMotifTypes, int );
    } 

  NEWP( bgPos, IP->nSeqLen, int );
  for( n = 0; n < IP->nSeqLen; n++ )
    {
      NEW( bgPos[n],  IP->nNumMotifTypes, int );
    } 

  nSites = 0;
  nSiteCount = 0;
  nMotifType = 0;
  while( nSites < 3 * RP->dPriorSites )
    {
      n = RandomInterval( 0, IP->nSeqLen - 1 );
      if( PossStartPos( Pos, n, nMotifType, B ) )
	{
	  bOK = TRUE;
	  for( t = 0; t < IP->nNumMotifTypes && bOK; t++ )
	    {
	      for( i = 0; i < IP->nMotifLen[nMotifType]; i++ )
		{
		  if( bgPos[n+i][t] )
		    {
		      bOK = FALSE;
		      break;
		    }		  
		}
	    }

	  if( bOK )
	    {
	      for( i = 0; i < IP->nMotifLen[nMotifType]; i++ )
		{
		  bgPos[n+i][nMotifType] = TRUE;
		  index = (int)(cSeq[n+i] - 97);
		  for( t = 0; t < IP->nNumMotifTypes; t++ )
		    bgCounts[index][t]++;
		  nSiteCount++;
		}
	      nSites++;
	      nMotifType = (nMotifType + 1) % IP->nNumMotifTypes;
	    }
	}

      nTrys++;
      if( nTrys > MAX_TRYS )
	p_error("UpdateBackgroundModel: Unable to sample required number of background sites." );
    }

  for( totalMotifs = 0, t = 0; t < IP->nNumMotifTypes; t++) 
    {
      totalMotifs += NUMMOTIFS(IP->nNumMotifs[t]) * IP->nMotifLen[t];
      
      for( dPseudo = 0.0, k = 0; k < IP->nAlphaLen; k++ )  
	dPseudo += C->dPseudoCounts[t][BG][k];
    }
  
  dBGDenom = nSiteCount + dPseudo;
  
  for(n = 0; n < IP->nAlphaLen; n++) 
    {
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  P->dvInBGProb[t][n] = ((double)bgCounts[n][t] +         /* BT 12/21/98 */
				 C->dPseudoCounts[t][BG][n]) / dBGDenom;  
	}
    }

  FREEP( bgCounts, IP->nAlphaLen );
  FREEP( bgPos, IP->nSeqLen );
}


void RemoveAllSites( Model B, PoSition **Pos, Mlist M, int ****occurence, 
		     ProbStruct *P )
{
  int        first;
  int        last;
  int        seq;
  int        n;
  int        t;
  int        i;
  short      RevComp;
 
  for(seq = 0; seq < B->IP->nNumSequences; seq++) 
    {
      first = SequenceStartPos( B, seq );
      
      for (t = 0; t < B->IP->nNumMotifTypes; t++) 
	{
	  last = (*B->Seq->nvEndLocs)[seq] - (*B->IP).nMotifLen[t] + 1;
	      
	  for( n = first; n < last; n++ )    /* BT 8/15/97 */
	    {
	      if( Pos[t][n].nMotifStartPos )
		{
		  RevComp = Pos[t][n].RevComp;        
		  adjust_counts(B,DELETE,n,t,RevComp);
		  B->IP->nNumMotifs[t][RevComp]--;
		  not_in_motif(Pos,n,B,t);    
		  delete_motif(B, n, M, t);
		  update_prob(P, B, (! B->IP->is_defined[cl_b]));
		}
	    }
	}
    }  

   for( i = 0; i < B->IP->nNumSequences; i++ )
     B->RP->nSites[i] = 0;
}


void PrintProbModel( Model B, ProbStruct *P, int t )
{
  IPtype     IP;
  int        n;
  int        i;

  IP = B->IP;

  for( i = 0; i < IP->nMotifLen[t]; i++ )
    {
      for(n = 0; n < IP->nAlphaLen; n++) 
	{
	  printf( " %5.3f", P->dvInMotifProb[t][i][n] ); 
	}
      printf( "\n" );
    }
}


double ExDistrib( double p, double temp )
{
  if( p > 0 )
    /*    return( exp( log( p ) / temp ) ); */
    return( pow( p, 1 / temp ) );
  else
    return 0;
}


void DumpProbModels( Model B, ProbStruct *P, FILE *fpt )
{
  IPtype   IP;  
  int      t;
  int      j;
  int      n;

  IP = B->IP;

  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      fprintf( fpt, "Motif type: %d\n", t );
      for( j = 0; j < IP->nMotifLen[t]; j++ )
	{
	  for( n = 0; n < IP->nAlphaLen; n++ )
	    fprintf( fpt, "%8.4f ", P->dvInMotifProb[t][j][n] );
	  fprintf( fpt, "\n" );
	}
      fprintf( fpt, "\n" );      
    }
}


void DumpBackground( Model B )
{
  IPtype IP;
  RPType RP;
  ALType AP;
  int    i;
  int    j;
  int    n;
  char   fName[1024];
  FILE*  fpt;

  static int iter = 0;

  if( B->RP == NULL )
    return;
    
  IP = B->IP;
  RP = B->RP;
  AP = B->AP;

  strcpy( fName, "/local/bioproj/people/thompson/Dump/dump_" );
  sprintf( &fName[strlen(fName)], "%d", iter );

  fpt = fopen( fName, "w" );

  for( i = 0; i < IP->nNumSequences; i++ )
    {
      for( j = 0; j < SequenceLength( B, i ); j++ )
	{
	  fprintf( fpt, "%d %d ", i, j );
	  for( n = 0; n < IP->nAlphaLen; n++ )
	    fprintf( fpt, "%8.5f ", B->BP->dBkgndProb[i][j][n] );	  
	  fprintf( fpt, "\n" );
	}      
    }

  iter++;
  fclose( fpt );
}


void AdjustSiteSeqCounts( Model B, PoSition **Pos, Mlist M )
{
  int      n;
  int      t;
  int      seq;
  int      first;
  int      last;
  int      RevComp;
  RPType   RP;
  IPtype   IP;
  PoSition **startPos;

  RP = B->RP;
  IP = B->IP;
    
  NEWP( startPos, B->IP->nNumMotifTypes, PoSition );
  for( t=0; t<B->IP->nNumMotifTypes; t++)
    {
      NEW( startPos[t], IP->nSeqLen, PoSition );
      for( n = 0; n < IP->nSeqLen; n++ )
	{
	  startPos[t][n] = Pos[t][n];
	}
    }

  for(seq = 0; seq < B->IP->nNumSequences; seq++) 
    RP->nSites[seq] = 0;

  for(seq = 0; seq < B->IP->nNumSequences; seq++) 
    {
      first = SequenceStartPos( B, seq );
      
      for (t = 0; t < B->IP->nNumMotifTypes; t++) 
	{
	  last = SequenceEndPos( B, seq );
	      
	  for( n = first; n <= last; n++ )    /* BT 8/15/97 */
	    {
	      if( Pos[t][n].nMotifStartPos )
		{
		  RevComp = Pos[t][n].RevComp;        
		  adjust_counts(B,DELETE,n,t,RevComp);
		  B->IP->nNumMotifs[t][RevComp]--;
		  not_in_motif(Pos,n,B,t);    
		  delete_motif(B, n, M, t);
		}
	    }
	}
    }

  for(seq = 0; seq < B->IP->nNumSequences; seq++) 
    {
      first = SequenceStartPos( B, seq );
      last = SequenceEndPos( B, seq );
      
      for (t = 0; t < B->IP->nNumMotifTypes; t++) 
	{	      
	  for( n = first; n <= last; n++ )    /* BT 8/15/97 */
	    {
	      if( RP->nSites[seq] < RP->nMaxBlocks )
		{
		  if( startPos[t][n].nMotifStartPos )
		    {
		      RevComp = Pos[t][n].RevComp;        
		      adjust_counts(B,ADD,n,t,RevComp);
		      IP->nNumMotifs[t][RevComp]++;
		      set_in_motif( Pos, n, B, t, RevComp );    
		      add_motif( IP, B->Seq, n, M, t, RevComp );
		      RP->nSites[seq]++;
		    }
		}
	    }
	}
    }
  
  FREEP( startPos, IP->nNumMotifTypes );  
}


void RemoveSites( Model B, int seq,  PoSition **Pos, Mlist M )
{
  int t;
  int first;
  int last;
  int n;
  int RevComp;
  int p;
  int nOffset;

  first = SequenceStartPos( B, seq );
  nOffset = SequenceLength( B, seq );

  for (t = 0; t < B->IP->nNumMotifTypes; t++) 
    {
      last = SequenceEndPos( B, seq ) - B->IP->nMotifLen[t] + 1;
      
      for( n = first; n <= last; n++ )    /* BT 8/15/97 */
	{
	  if( Pos[t][n].nMotifStartPos )
	    {
	      RevComp = Pos[t][n].RevComp;        
	      adjust_counts(B,DELETE,n,t,RevComp);
	      B->IP->nNumMotifs[t][RevComp]--;
	      not_in_motif(Pos,n,B,t);    
	      delete_motif(B, n, M, t);
	      if( IsPhyloSeq( B, seq ) && ! IsPhyloSeqBySpecies( B, seq ) )
		{
		  for( p = 1; p < SeqInc( B, seq ); p++ )
		    {
		      adjust_counts(B,DELETE,n+ p * nOffset,t,RevComp);
		      B->IP->nNumMotifs[t][RevComp]--;
		      not_in_motif(Pos, n + p * nOffset,B,t);    
		      delete_motif(B, n + p * nOffset, M, t);
		    }
		}
	    }
	}
    }

  B->RP->nSites[seq] = 0;
  if( IsPhyloSeq( B, seq ) )
    for( p = 1; p < SeqInc( B, seq ); p++ )
      B->RP->nSites[seq+p] = 0; 
}
