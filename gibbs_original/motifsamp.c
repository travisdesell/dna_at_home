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
/***************************************************************************/
/* $Id: motifsamp.c,v 1.7 2007/05/23 18:19:57 Bill Exp $                                                                  */
/* AUTHOR    : Eric C. Rouchka,  July 2, 1996                              */
/*             Jun Zhu,   October 20, 1996                                 */
/*                                                                         */
/***************************************************************************/

#include "motifsamp.h"


void PrintFootprint( Model B, int n, int t );


/*******************  sample_motif  ***************************************/
/* PARAMETERS :                                                           */
/*                                                                        */
/* DESCRIPTION :                                                          */
/*    Sampling motifs according to posterior probability distribution     */
/*    of elements in the model.                                           */
/*========================================================================*/


void sample_motif(Model B, PoSition **Pos, Mlist M, int iterations, 
		  int ****occurence, ProbStruct *P)

{
  double tprob, rand_num;
  int i,t, n;
  int count;
  short FOUND=FALSE;
  short RevComp;
  double sum;
  double **dFinalProb;
  double **dLenProb;      /* BT 9/3/97 */
  int    nMotifType = -1; 
  int    cnt = 0;        /* BT 11/5/97 */
  int    tcnt;           /* BT 11/5/97 */
  double dMotifProb;
  double dBGProb;

#ifdef _DEBUG_
  double pf[GRAPH_SIZE];
  double pr[GRAPH_SIZE];

  for( n = 0; n < GRAPH_SIZE; n++ )
    {
      pf[n] = 0;
      pr[n] = 0;
    }

#endif
  
  /* dFinalProb[Motiftype][FORWARD/REVERSE] */
  /* dFinalProb is the probability that the segment derived from */
  /* the model motif */
   NEWP(dFinalProb, B->IP->nNumMotifTypes, double);
   for(t = 0; t < B->IP->nNumMotifTypes; t++)
     NEW(dFinalProb[t], 2, double);
   
   n=0;
   /* we should avoid to always sample from type 1 to N */
   count=0;
   while(n<B->IP->nSeqLen)		/* BT 4/14/97 */
     {
       tprob = 0.0;
       for(t = 0; t < B->IP->nNumMotifTypes; t++)
	 {
	   dFinalProb[t][FORWARD] = 0.0;
	   dFinalProb[t][REVERSE] = 0.0;
	   
	   if(PossFragStartPos(Pos, n, t, B) &&
	      !OverlapsPrevMotif(n,B->IP->nNumMotifTypes,Pos))
	     {
	       if( (! B->IP->is_defined[cl_T]) || InFootprintRegion( B, n, t ) )
		 {
		   CheckMotifProb(B, n, t, Pos, P, M, iterations,&dFinalProb);  /* BT 11/09/99 */
		   if( B->IP->is_defined[cl_X] && B->IP->nSeedRun > 0 )
		     {
		       dFinalProb[t][FORWARD] =  ExDistrib( dFinalProb[t][FORWARD], B->AN->currTemp );
		       if( B->IP->RevComplement )
			 dFinalProb[t][REVERSE] =  ExDistrib( dFinalProb[t][REVERSE], B->AN->currTemp );
		     }
		   tprob += dFinalProb[t][FORWARD]; /* Find the total prob that */
		   /* a position begins here   */
		   if( B->IP->RevComplement )
		     tprob += dFinalProb[t][REVERSE];
		 }
	     }
	 } 

#ifdef _DEBUG_
       if( n < GRAPH_SIZE )
	 {
	   pf[n] = dFinalProb[0][FORWARD];
	   pr[n] = dFinalProb[0][REVERSE];
	 }
       CheckCounts( B );
#endif

       FOUND = FALSE;
       if( tprob > 0.0 )
	 {
	   sum = 0.0;
	   RevComp = FALSE;
	   rand_num = ((double) 1.0 + tprob) * drand();    /* BT 7/11/97 */

	   cnt = cnt %  B->IP->nNumMotifTypes;    /* BT 11/5/97 */
	   for( tcnt = cnt; tcnt < cnt + B->IP->nNumMotifTypes && !FOUND; tcnt++ )
	     {
	       t = tcnt %  B->IP->nNumMotifTypes;
	       sum += dFinalProb[t][FORWARD]; 
	       if( sum >= rand_num )
		 {
		   FOUND = TRUE;
		   nMotifType = t;  
		 } 
	       else
		 {
		   sum += dFinalProb[t][REVERSE];
		   if( sum >= rand_num )
		     {
		       FOUND = TRUE;
		       nMotifType = t;
		       RevComp = TRUE;   
		     }  
		 }
	     } 
	   
	   if( FOUND )
	     {
	       if( B->IP->is_defined[cl_g] )    /* BT 9/10/97 */
		 {
		   NEWP( dLenProb, 2, double );
		   for( i = 0; i < 2; i++ )
		     NEW( dLenProb[i], B->IP->nMotifLen[nMotifType], double );
		   
		   for( sum = 0.0, i = 0; i < B->IP->nMotifLen[nMotifType]; i++ )  
		     {
		       if( PossibleMotifPos( n + i, Pos, B->IP, nMotifType, B ) )
			 {
			   dLenProb[FORWARD][i]=in_motif_prob(B,*P,n+i,nMotifType,FALSE, 
							      TRUE,  &dMotifProb, &dBGProb);
			   if( B->IP->RevComplement )
			     dLenProb[REVERSE][i]=in_motif_prob(B,*P,n+i,nMotifType,TRUE, 
								TRUE,  &dMotifProb, &dBGProb);
			   sum += dLenProb[FORWARD][i] + dLenProb[REVERSE][i];
			 }
		     }
		 
		   rand_num = drand();   /* BT 7/11/97 */
		   for( tprob = 0.0, i = 0; i < B->IP->nMotifLen[nMotifType]; i++ )
		     {
		       tprob += dLenProb[FORWARD][i] / sum;
		       if( tprob >= rand_num )
			 {
			   n += i;
			   RevComp = FALSE;
			   break;
			 }
		       tprob += dLenProb[REVERSE][i] / sum;  
		     if( tprob >= rand_num )
		       {
			 n += i;
			 RevComp = TRUE;
			 break;
		       }
		     }
		   
		   FREEP(dLenProb, 2);   /* BT 9/3/97 */
		 }

	       (*occurence)[n][nMotifType][0]++;
	       if(RevComp) { (*occurence)[n][nMotifType][2]++; }
	       else        { (*occurence)[n][nMotifType][1]++; }
	       P->update = TRUE;
	       adjust_counts(B,ADD,n,nMotifType,RevComp);
	       B->IP->nNumMotifs[nMotifType][RevComp]++;
	       set_in_motif(Pos, n, B, nMotifType, RevComp);
	       add_motif(B->IP, B->Seq,n, M, nMotifType, RevComp);

	       if(iterations > 3) 
		 update_posterior_prob(B);
	  
	       /* skip the length of this motif */
	       n+=(B->IP->nMotifLen[nMotifType]-1);   /* -1 here, because n will */
	                                              /* increase by 1 at the end*/
                                         	      /* of for loop */
	       cnt++;                                 /* BT 11/5/97 */
	     }
	 }
       n++;
     }
     
   FREEP(dFinalProb, B->IP->nNumMotifTypes);
}


MaxResults nearopt_motif_samp(Model B, PoSition **Pos, Mlist M, int **good,  
                              ProbStruct *P, int ****occurence) 
{
   int i, n, t, last_increase;
   double dCurrProb, dLocMax, dMaxProbability=0.0;
   double **dProbArray;
   MaxResults currMax;
   int **posSave;			/* BT 3/26/97 */
   int nMotifs;                         /* BT 5/21/97 */
#ifdef _GUI_
   XEvent event;                        /* BT 6/17/97 */
#endif

   NEWP( posSave, B->IP->nNumMotifTypes, int);                /* BT 3/26/97 */
   for(i = 0; i < B->IP->nNumMotifTypes; i++)
      NEW( posSave[i], B->IP->nSeqLen, int);

   init_maxdata(&currMax);
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

   FREEP(good, B->IP->nSeqLen);

   /* Calculate initial probability */          /* BT 7/30/97 */
   for(nMotifs = 0,  t=0; t < B->IP->nNumMotifTypes; t++)
     {
       nMotifs += B->IP->nNumMotifs[t][FORWARD];     
       if( B->IP->RevComplement )
	 nMotifs += B->IP->nNumMotifs[t][REVERSE];
     }
   dCurrProb = CalcMapProb(B, B->IP->is_defined[cl_R]);
   dProbArray = setElementProb(B, Pos, *P);
   currMax = setMaxData(B->F, B->IP, dCurrProb, i, Pos, &last_increase,
			&dMaxProbability, &dLocMax, currMax, dProbArray, B);
   FREEP(dProbArray, findMaxNumMotif(B->IP)); 

   for(i = 1; i <= B->IP->nMaxIterations; i++) {
     if( (! B->IP->is_defined[cl_Z]) && i % 5 == 0 )
	{
	  fprintf(stdout, "\r%d", i);
 	  fflush( stdout );      /* BT 9/19/97 */
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

      if( ! B->IP->is_defined[cl_D] )
	sample_motif(B, Pos, M, i, occurence, P);
      else
	sample_motif_pair(B, Pos, M, i, occurence, P);

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
		i, CalcMapProb(B, B->IP->is_defined[cl_R]));
#endif

      for(nMotifs = 0, t=0; t < B->IP->nNumMotifTypes; t++)
	{
	  nMotifs += B->IP->nNumMotifs[t][FORWARD];      /* BT 5/21/97 */
	  if( B->IP->RevComplement )
	    nMotifs += B->IP->nNumMotifs[t][REVERSE];
	}
      dCurrProb = CalcMapProb(B, B->IP->is_defined[cl_R]);
      /*      if( ((dCurrProb > dLocMax) && (nMotifs != 0)) || (i == 1)) */      /* BT 5/21/97 */
      if( (dCurrProb > dLocMax) && (nMotifs != 0))     /* BT 5/21/97 */
	{ 
	  if( ! B->IP->is_defined[cl_Z] )
	    printf("Max set %f at %d\n", dCurrProb, i);
          dProbArray = setElementProb(B, Pos, *P);
          free_maxdata(&currMax, B->IP);
          currMax = setMaxData(B->F, B->IP, dCurrProb, i, Pos, &last_increase,
			       &dMaxProbability, &dLocMax, currMax, dProbArray, B);
          FREEP(dProbArray, findMaxNumMotif(B->IP));
	}
   }

   /* Restore possible site starting positions for next phase */
   for(t = 0; t < B->IP->nNumMotifTypes; t++)		/* BT 3/26/97 */
   {         
     for(n = 0; n < B->IP->nSeqLen; n++)
       Pos[t][n].nPossStartPos = posSave[t][n]; 
   }
   FREEP( posSave, B->IP->nNumMotifTypes);

   return currMax;
}


void DumpMotifPositions( int t, Model B, PoSition **Pos, FILE *fpt )
{
  IPtype      IP;
  PoSition    *pos_t;
  int         i;
  int         last_seq = -1;
  int         site_num = -1;
  double      prob;
  double      **dProbArray;
  int         count;
  ProbStruct  P;

  pos_t = Pos[t];
  IP = B->IP;
  fprintf( fpt, "\n" );

  init_prob( &P, IP );

  dProbArray = setElementProb(B, Pos, P);

  count = 0;
  for( i = 0; i < IP->nSeqLen; i ++ )
    {
      if( pos_t[i].nMotifStartPos )
	{
	  /*	  prob = in_motif_prob(B, *P, i, t,  pos_t[i].RevComp, FALSE, 
			       &dMotifProb, &dBGProb); */
	  prob = dProbArray[count][t];
	  print_it( B, t, i, fpt, prob, 
		    pos_t[i].RevComp, &last_seq, &site_num, NULL, 0);
	  count++;
	}
    }

  free_prob( &P, IP );
  FREEP(dProbArray, findMaxNumMotif(B->IP));
}


short InFootprintRegion( Model B, int n, int t )
{
  int nSeq;
  int nPos;
  int nOffset;

  if( ! B->IP->is_defined[cl_T] )
    return TRUE;

  nSeq = SequenceFromPosition( B, n );
  nPos = n - SequenceStartPos( B, nSeq );

  nOffset = nPos + MotifWidth( B, t ) / 2;

  if( nOffset > SequenceLength( B, nSeq ) )
    return FALSE;
  else if (B->RP->dProbCons[nSeq][nOffset] >= FOOTPRINT_CUTOFF)
    return TRUE;
  else 
    return FALSE;    
}


void PrintFootprint( Model B, int n, int t )
{
  int nSeq;
  int nPos;
  int nOffset;

  nSeq = SequenceFromPosition( B, n );
  nPos = n - SequenceStartPos( B, nSeq );
  nOffset = nPos + B->IP->nMotifLen[t] / 2;

  fprintf( stdout, "seq = %d, pos = %d motif = %d n = %d offset = %d fp = %f\n",
	   nSeq, nPos, t, n, nOffset, 
	   B->RP->dProbCons[nSeq][nOffset] );
}


