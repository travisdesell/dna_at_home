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
#include "motifpair.h"


void sample_motif_pair(Model B, PoSition **Pos, Mlist M, int iterations, 
		  	int ****occurence, ProbStruct *P)
{
  double tprob;
  double rand_num;
  int i,t, n;
  int count;
  short FOUND=FALSE;
  short RevComp;
  double sum;
  double **dFinalProb;
  double **dLenProb;      /* BT 9/3/97 */
  int    nMotifType = 0; 
  int    cnt = 0;        /* BT 11/5/97 */
  int    tcnt;           /* BT 11/5/97 */
  double dMotifProb;
  double dBGProb;
  int    seq;
  int    first;
  int    lastpos;
  int    nOffset;
#ifdef _DEBUG_
  double prf[GRAPH_SIZE];
  double prr[GRAPH_SIZE];
#endif
  
  /* dFinalProb is the probability that the segment derived from */
  /* the model motif */

  NEWP(dFinalProb, B->IP->nNumMotifTypes, double);
  for(t = 0; t < B->IP->nNumMotifTypes; t++)
    NEW(dFinalProb[t], 2, double);
       
   n=0;
   count=0;
   for(seq = 0; seq < B->IP->nNumSequences; seq += 2) 
     {
       first = SequenceStartPos( B, seq );       
       nOffset = SequenceLength( B, seq );
       lastpos = first + nOffset - 1;
       
       for( n = first; n < lastpos; n++ )
	 {
	   tprob = 0.0;
	   for(t = 0; t < B->IP->nNumMotifTypes; t++)
	     {
	       dFinalProb[t][FORWARD] = 0.0;
	       dFinalProb[t][REVERSE] = 0.0;
	   
	       if(PossFragStartPos(Pos, n, t, B) &&
		  PossFragStartPos(Pos, n + nOffset, t, B) &&
		  !OverlapsPrevMotif(n, B->IP->nNumMotifTypes, Pos) &&
		  !OverlapsPrevMotif(n + nOffset, B->IP->nNumMotifTypes, Pos))
		 {
		   CheckMotifProb(B, n, t, Pos, P, M, iterations,&dFinalProb);
		   if( B->IP->is_defined[cl_X] && B->IP->nSeedRun > 0 )
		     {
		       dFinalProb[t][FORWARD] =  ExDistrib( dFinalProb[t][FORWARD], B->AN->currTemp );
		       if( B->IP->RevComplement )
			 dFinalProb[t][REVERSE] =  ExDistrib( dFinalProb[t][REVERSE], B->AN->currTemp );
		     }
		   tprob += dFinalProb[t][FORWARD];
		   if( B->IP->RevComplement )
		     tprob += dFinalProb[t][REVERSE];
		 }
	     }

#ifdef _DEBUG_	   
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
		       sum +=dFinalProb[t][REVERSE];
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
		   
		   (*occurence)[n + nOffset][nMotifType][0]++;
		   if(RevComp) { (*occurence)[n + nOffset][nMotifType][2]++; }
		   else        { (*occurence)[n + nOffset][nMotifType][1]++; }
		   P->update = TRUE;
		   adjust_counts(B,ADD,n + nOffset,nMotifType,RevComp);
		   B->IP->nNumMotifs[nMotifType][RevComp]++;
		   set_in_motif(Pos, n + nOffset, B, nMotifType, RevComp);
		   add_motif(B->IP, B->Seq,n + nOffset, M, nMotifType, RevComp);
		   
		   if(iterations > 3) 
		     update_posterior_prob(B);
		   
		   /* skip the length of this motif */
		   n+=(B->IP->nMotifLen[nMotifType]-1);
		   cnt++;                                  
		 }
	     }
	 }

#ifdef _DEBUG_	   
       if( seq == 0 )
	 {
	   for(t = 0; t < B->IP->nNumMotifTypes; t++)
	     {
	       for( n = first; n <= min(lastpos - first, GRAPH_SIZE - 1); n++ )
		 {
		   if(PossFragStartPos(Pos, n, t, B) &&
		      PossFragStartPos(Pos, n + nOffset, t, B) )
		     {
		       prf[n - first] = in_motif_prob(B,*P,n,t,FALSE, 
						      TRUE,  &dMotifProb, &dBGProb);		   
		       prr[n - first] = in_motif_prob(B,*P,n,t,TRUE, 
						      TRUE,  &dMotifProb, &dBGProb);
		     }
		   else
		     {
		       prf[n - first] = 0;
		       prr[n - first] = 0;
		     }
		 }
	       seq = seq;
	     }
	 }
#endif
       
     }
       
   FREEP(dFinalProb, B->IP->nNumMotifTypes);
}
