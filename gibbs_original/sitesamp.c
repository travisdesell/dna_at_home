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
/* $Id: sitesamp.c,v 1.6 2007/05/23 18:19:58 Bill Exp $         */
/* AUTHOR      : Eric C. Rouchka  July 11, 1996                          */
/*               Jun Zhu October 20, 1996                                */
/*               Bill Thompson  1/27/97                                  */
/*                                                                       */
/* DESCRIPTION : This file contains the site sampler code that will go   */
/*               through each sequence and sample one and only one motif */
/*               of each type with the advantage being that none of the  */
/*               motif elements of one type will overlap another, and    */
/*               thus an alignment that is more likely to be maximal will*/
/*                result                                                 */
/*************************************************************************/

#include "sitesamp.h"


MaxResults nearopt_site_samp(Model B, PoSition **Pos,  Mlist M, int **good, 
                             ProbStruct *P, int ****occurence)

   /*=======================================================================*/
   /* FUNCTION NAME : nearopt_site_samp                                     */
   /*                                                                       */
   /*=======================================================================*/
   
{
   int     i = 0, seq, first, last, n, t, index, last_increase, j;  /* BT 1/27/97 */
   double  dCurrProb, dLocMax, dMaxProbability=0.0;
   double  **dProb, **dFinalProb, **dProbArray;
   double  dtot, rand_num;
   unsigned short type = FORWARD, found;
   MaxResults currMax;
#ifdef _GUI_
   XEvent event;                        /* BT 6/17/97 */
#endif
   int    nMotifs;                      /* BT 7/30/97 */
   short  RevComp;                      /* BT8/15/97 */

   init_maxdata(&currMax);
   NEWP(dFinalProb, B->IP->nSeqLen, double);
   for(n = 0; n < B->IP->nSeqLen; n++)
      NEW(dFinalProb[n], 2, double);
   
   NEWP(dProb, B->IP->nSeqLen, double);
   for(n = 0; n < B->IP->nSeqLen; n++)
      NEW(dProb[n], 2, double);

   /* Calculate initial probability */          /* BT 7/30/97 */
   for(nMotifs = 0, t=0; t < B->IP->nNumMotifTypes; t++)
     {
       nMotifs += B->IP->nNumMotifs[t][FORWARD];      /* BT 5/21/97 */
       if( B->IP->RevComplement )
	 nMotifs += B->IP->nNumMotifs[t][REVERSE];
     }
   dCurrProb = CalcMapProb(B, B->IP->is_defined[cl_R]);
   dProbArray = setElementProb(B, Pos, *P);
   currMax = setMaxData(B->F, B->IP, dCurrProb, i, Pos, &last_increase,
			&dMaxProbability, &dLocMax, currMax, dProbArray, B);
   FREEP(dProbArray, findMaxNumMotif(B->IP));

   for(i = 1; i <= B->IP->nMaxIterations; i++) {
     if( ! B->IP->is_defined[cl_Z] && i % 5 == 0)
	{
	  fprintf(stdout, "\r%d", i);
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
     
      for(seq = 0; seq < B->IP->nNumSequences; seq++) 
	{
	  first = SequenceStartPos( B, seq );

	  for(t = 0; t < B->IP->nNumMotifTypes; t++) 
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
		      update_prob(P, B, (! B->IP->is_defined[cl_b]));
		      break;
		    }
		}

	      dtot = 0.0;
	      
             for(n = first; n <= last; n++) {
               init_values(&dFinalProb, Pos, n, t);
               if(PossFragStartPos(Pos, n, t, B) &&
                   good[n][t] &&
                   !PossibleOverlapsPrevMotif( B, n, t, Pos ) )   /* BT 4/4/97 */
		 CheckMotifProb(B, n, t, Pos, P, M, i, &dFinalProb);
	       dProb[n][0] = dFinalProb[t][0];
               dProb[n][1] = dFinalProb[t][1];
               dtot += dFinalProb[t][0] + dFinalProb[t][1];
            }

            for(n = first; n <= last; n++) 
	      {
		/* if(dtot > 0.0000001)  */
		if(dtot > 0.0)  /* BT 7/21/05 */
		  {
		    dProb[n][0] /= dtot;
		    if(B->IP->RevComplement)
		      dProb[n][1] /= dtot;
		  }
	      }
            rand_num = (double)Random()/(double)INT32_MAX;
            dtot = 0.0;
            index = first-1;
            found = FALSE;
            while(!found) {
               index++;
               if(index > last) 
		 p_internal_error("Site Sampler - cannot find site \n"); 	       
               dtot += dProb[index][0];
               if(dtot >= rand_num)
                    { found = TRUE; type = FORWARD; }
               else if(B->IP->RevComplement){
                  dtot += dProb[index][1];
                  if(dtot >= rand_num)
                    { found = TRUE; type = REVERSE; }
               }
            }
            (*occurence)[index][t][0]++;
            if(type == FORWARD) (*occurence)[index][t][1]++;
            else                (*occurence)[index][t][2]++;
            adjust_counts(B, ADD, index,t, type);
	    B->IP->nNumMotifs[t][type]++;
            add_motif(B->IP, B->Seq, index, M, t, type);
	    Pos[t][index].nMotifStartPos=TRUE;
	    Pos[t][index].RevComp=type;
	    for(j=0;j<B->IP->nMotifLen[t];j++)    /* BT 1/27/97 */
	      Pos[t][index+j].nInMotif=TRUE;
         }
      }
      dCurrProb  = CalcMapProb(B, B->IP->is_defined[cl_R]);
      if((dCurrProb > dLocMax) /*|| (i == 1) */)      /* BT 7/30/97 */
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
   FREEP(good, B->IP->nSeqLen);
   FREEP(dProb, B->IP->nSeqLen);
   FREEP(dFinalProb, B->IP->nSeqLen);
   return (currMax);
}


/*************    MaxResult site_sampler *********************************/
/*                                                                       */
/*    Pos[MotifType][SeqLen]                                             */
/*=======================================================================*/

MaxResults site_sampler(Model B, PoSition **Pos, Mlist M)
 
{
   int             iter, last_increase, seq, first, last;
   int             n, t, index;
   double          dLocMax, dCurrProb, dMaxProbability = 0.0;
   double          **dProb, **dFinalProb, dtot, rand_num;
   ProbStruct      P;
   MaxResults      maxData;
   unsigned short  type = 0, found;
#ifdef _GUI_
   XEvent          event;                        /* BT 6/17/97 */
#endif
   short           RevComp;                      /* BT 8/15/97 */

   init_maxdata(&maxData);
   NEWP(dFinalProb, B->IP->nSeqLen, double);
   for(n = 0; n < B->IP->nSeqLen; n++)
      NEW(dFinalProb[n], 2, double);
   
   NEWP(dProb, B->IP->nSeqLen, double);
   for(n = 0; n < B->IP->nSeqLen; n++)
      NEW(dProb[n], 2, double);

   reset_values(&iter, &last_increase, &dLocMax, &P, B->IP); 
   while(((iter - last_increase) < B->IP->nPlateauPeriods) &&
         (iter < B->IP->nMaxIterations)) 
     {
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
      
      for(seq = 0; seq < B->IP->nNumSequences; seq++) 
	{
	  first = SequenceStartPos( B, seq );

	  for (t = 0; t < B->IP->nNumMotifTypes; t++)   /* BT 12/07/2000 */
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
		      break;
		    }
		}
	    }
	  update_prob(&P, B, (! B->IP->is_defined[cl_b]));
	      
	  for (t = 0; t < B->IP->nNumMotifTypes; t++) 
	    {
	      last = SequenceEndPos( B, seq ) - B->IP->nMotifLen[t] + 1;
	      dtot = 0.0;
	      for(n = first; n <= last; n++) 
		{
		  init_values(&dFinalProb, Pos, n, t);
		  if(PossFragStartPos(Pos, n, t, B) &&
		     !PossibleOverlapsPrevMotif( B, n, t, Pos ) )   /* BT 4/4/97 */
		    CheckMotifProb(B, n, t, Pos, &P, M, iter, &dFinalProb);
		  dProb[n][0] = dFinalProb[t][0];
		  dProb[n][1] = dFinalProb[t][1];
		  dtot += dFinalProb[t][0] + dFinalProb[t][1];
		}

	      for(n = first; n <= last; n++) 
		{
		  dProb[n][0] /= dtot;
		  if(B->IP->RevComplement) 
		    dProb[n][1] /= dtot;
		}
	      rand_num = (double)Random()/(double)INT32_MAX;
	      dtot = 0.0;
	      index = first-1;
	      found = FALSE;
	      while(!found) 
		{
		  index++;
		  if(index > last) {
		    print_Pos( Pos, B, M );                     /* BT 6/2/97 */
		    p_internal_error("Site sampler -- cannot find site");
		  }
		  dtot += dProb[index][0];
		  if(dtot >= rand_num)
                  { found = TRUE; type = FORWARD; }
		  else if(B->IP->RevComplement){
		    dtot += dProb[index][1];
		    if(dtot >= rand_num)
		      { found = TRUE; type = REVERSE; }
		  }
		}    
	      adjust_counts(B, ADD,index,t, type);
	      B->IP->nNumMotifs[t][type]++;
	      add_motif(B->IP, B->Seq, index, M, t, type);
	      set_in_motif(Pos, index, B, t, type);	      
	    }

	}
      if(iter % B->IP->nAdjustPeriod == 0) 
	{
	  if(!B->IP->is_defined[cl_F]) 
	    {
	      /*   setFragVec(B, M, Pos);
	      setFragCnts(B, M);
	      StoreMask(B);
	      for(t = 0; t < B->IP->nNumMotifTypes; t++)
		Fragment(B, t, M);
		adjust_fragment(B, M, Pos, &P); */
	       for(t = 0; t < B->IP->nNumMotifTypes; t++)
		 {
		   StoreMask(B); 
		   setFragVec(B, M, Pos);   /* BT 4/26/04 */
		   setFragCnts(B, M);
		   Fragment(B, t, M);
		   adjust_fragment(B, M, Pos, &P);
		 }
	    }
	  else
	    {
	      if( B->IP->is_defined[cl_d] )           /* BT 6/9/97 */
		ResizeMotif( B, Pos, M, FALSE );             /* BT 5/12/97 */
	      ColumnShift(B, M, Pos, &P);        
	     }
	}
      if(iter > 3)
	check_map(B, &dCurrProb, iter, Pos, &last_increase, &dMaxProbability, 
		  &dLocMax, &maxData, P);
     }
   free_prob(&P, B->IP); 
   maxData.frequency = NULL;
   FREEP(dFinalProb, B->IP->nSeqLen);
   FREEP(dProb, B->IP->nSeqLen);
   if(!B->IP->is_defined[cl_F]) 
     {  
       if( B->F->nvFragCnts )
	 {
	   for( t = 0; t < B->IP->nNumMotifTypes; t++ ) 	/* BT 4/2/97 */
	     FREEP( B->F->nvFragCnts[t], B->F->nMaxLen[t] );
	   free( B->F->nvFragCnts );
	 }
       B->F->nvFragCnts = NULL;
     }   
   
   /*      FREEPP(B->F->nvFragCnts, B->IP->nNumMotifTypes, B->F->nMaxLen[t]); */
   free(dFinalProb);
   free(dProb);
   return(maxData);
}





