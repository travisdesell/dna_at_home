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
/************************************************************************/
/* $Id: expectmax.c,v 1.6 2007/05/23 18:19:56 Bill Exp $       */
/*                                                                      */
/* AUTHOR      : Eric C. Rouchka,  August 5, 1996                       */
/*               Jun Zhu, October 20, 1996                              */
/*               Bill Thompson  2/14/97                                 */  
/*                                                                      */
/* DESCRIPTION : Contains the functions to perform the expectation      */
/*               maximization once a near optimal solution has been     */
/*               found.                                                 */
/*                                                                      */
/************************************************************************/

/*###########################################################*/
/* NOTE : STILL DOES NOT WORK CORRECTLY WITH MULTIPLE MOTIFS */
/* OR REVERSE COMPLEMENTS OR ANYTHING ELSE YET               */
/* NEED TO WORK ON IMPLEMENTING OTHER MODELS AS WELL         */
/*###########################################################*/

/* NEED TO WORK IN PALANDROMES, etc */

#include "expectmax.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define TOT_BG 0

/*----------------------------------------------------------------------*/
/*-----------------------   LOCAL FUNCTION PROTOTYPES ------------------*/
/*----------------------------------------------------------------------*/
void set_values(Model B, MaxResults optmax);
void print_header(FILE *fpt);
void init_expected(Model B, double ****expected, double ****StartMotifProb, 
                   PoSition **Pos, ProbStruct P);
double log_epsilon(double *total, double tot_back, double epsilon);
double log_cost(Model B, PoSition **Pos, double ***StartMotifProb);
double log_likelihood(IPtype IP, double ***expected, double *total, 
                      double tot_back);
void init_total(double **total, double ***expected, IPtype IP);
void print_EM_results(Model B, PoSition **Pos, double ***StartMotifProb);
void init_stats(EMStats *S);
void print_EM_stats(FILE *fpt, EMStats S);
void calculate_EM_stats(EMStats *S, Model B, double ***expected,
                        double *total, double tot_back, double epsilon,
                        PoSition **Pos, double ***StartMotifProb);

/*-----------------------------------------------------------------------*/
/*--------------------------- FUNCTION DECLARATIONS ---------------------*/
/*-----------------------------------------------------------------------*/

void set_values(Model B, MaxResults optmax)

/*=======================================================================*/
/* FUNCTION NAME : set_values                                            */
/*                                                                       */
/* DESCRIPTION : Sets the number of motifs in both the forward and       */
/*               reverse complement directions according to the optimal  */
/*               alignment results                                       */
/*=======================================================================*/

{
   int t, n;

   for(t = 0; t < B->IP->nNumMotifTypes; t++) {
     B->IP->nMotifLen[t] = optmax.nMotifLen[t];   /* BT 5/23/97 */
      B->C->nTotBack += optmax.nNumMotifs[t] * B->IP->nMotifLen[t];
      B->IP->nNumMotifs[t][FORWARD] = 0;
      B->IP->nNumMotifs[t][REVERSE] = 0;
      for(n = 0; n < B->IP->nSeqLen; n++) {
        if(optmax.frequency[n][t].frequency > B->IP->dNearOptCutoff ) {
          B->C->nTotBack -= B->IP->nMotifLen[t];
          if(optmax.frequency[n][t].fwd_freq > optmax.frequency[n][t].rev_freq)
             B->IP->nNumMotifs[t][FORWARD]++;
          else
             B->IP->nNumMotifs[t][REVERSE]++;
         }
      }
   }
}


void print_header(FILE *fpt)
{
   fprintf(fpt, "\t\tEXPECTATION/MAXIMIZATION PARAMETERS\n");
   fprintf(fpt, "\t\t___________________________________\n\n");
   fprintf(fpt, "   logL    Tot pos.  Cost     Total     Change    Pct. ");
   fprintf(fpt, "logL Cng,    %%\n");
   fprintf(fpt, "---------- -------- ------ ---------- ---------- ----- ");
   fprintf(fpt, "---------- -----\n");
}

void init_expected(Model B, double ****expected, double ****StartMotifProb, 
                   PoSition **Pos, ProbStruct P)


/*===========================================================================*/
/* FUNCTION NAME : init_expected                                             */
/*                                                                           */
/* DESCRIPTION : This is the expected step of the Expectation/Maximization   */
/*               algorithm.  It calculates the counts of each residue in each*/
/*               PoSition of the motif based on the probability of each      */
/*               position being a motif start position                       */
/*===========================================================================*/

{
   int n, t, index, j, i, start, end, palindex = 0, mid;
   double prob, ratio, tot;
   double fwd=0.0, rev=0.0;
   int nLen;			/* BT 4/21/97 */
   double dMotifProb;
   double dBGProb;
   int    repindex = 0;
   
   ratio = (double)B->IP->nNumMotifs[0][FORWARD] / 
          (double)(B->IP->nNumMotifs[0][FORWARD] + B->IP->nNumMotifs[0][REVERSE]);
  
   for(n = 0; n < B->IP->nAlphaLen; n++) 
      for(t = 0; t < B->IP->nNumMotifTypes; t++) 
         for(j = 0; j <= B->IP->nMotifLen[t]; j++) 
            (*expected)[n][t][j] = 0.0;

   for(t = 0; t < B->IP->nNumMotifTypes; t++)			/* BT 2/24/97 */
   {
      for( n = 0; n < B->IP->nSeqLen; n++)
      {
         (*StartMotifProb)[n][t][FORWARD] = (*StartMotifProb)[n][t][REVERSE] = 0.0;
      }
    }

/****************** 
   for(i = 0; j < B->IP->nNumSequences; i++) {
      if(i == 0) start = 0;
      else       start = (*B->Seq->nvEndLocs)[i-1];
      end = (*B->Seq->nvEndLocs)[i] - 1;
      for (tot = 0.0, n = start; n <= end; n++) {
         for(t = 0; t < B->IP->nNumMotifTypes; t++) {
            if(PossFragStartPos(Pos, n, t, B)) {
               (*StartMotifProb)[n][t][FORWARD] = in_motif_prob(B, P, n, t, FALSE);
               if(B->IP->RevComplement)
                  (*StartMotifProb)[n][t][REVERSE] = in_motif_prob(B, P, n, t,TRUE);
               prob = (*StartMotifProb)[n][t][FORWARD];
               for (j = 1; j <= B->IP->nMotifLen[t]; j++) {
                  index = CH2INT(n + j - 1, B);
                  prob = (*StartMotifProb)[n][t][FORWARD];
                  if(B->IP->RevComplement) {
                     if(j == 1) {
                        fwd += (*StartMotifProb)[n][t][FORWARD];
                        rev += (*StartMotifProb)[n][t][REVERSE];
                     }
                     prob *= ratio;   tot += prob;
                     (*expected)[index][t][j] += prob;
                     (*expected)[index][TOT_BG][BG] += ratio - prob;
                     index = CH2INTCOMP(n + B->IP->nMotifLen[t] - j, B);
                     tot += prob;
                     (*expected)[index][t][j] += prob;
                     (*expected)[index][TOT_BG][BG] += 1.0 - ratio - prob;
                  }
                  else {
                     tot += prob;
                     (*expected)[index][t][j] += prob;
                     (*expected)[index][TOT_BG][BG] += 1.0 - prob;
                  }
               }
            }
         }
      }
      if(B->IP->site_samp)
         for(n = start; n <= end; n++) 
            for(t = 0; t < B->IP->nNumMotifTypes; t++) 
               for(j = 0; j < B->IP->nAlphaLen; j ++) {
                  (*expected)[j][t][n] /= tot;
                  (*expected)[j][TOT_BG][BG] = 1.0 - (*expected)[j][t][n];
               }
   }
   if(B->IP->RevComplement) {
      B->IP->nNumMotifs[0][FORWARD] = (int)fwd;
      B->IP->nNumMotifs[0][REVERSE] = (int)rev;
   }
   for(n = 0; n < B->IP->nAlphaLen; n++) 
      for(t = 0; t < B->IP->nNumMotifTypes; t++) 
         (*expected)[n][TOT_BG][BG] /= B->IP->nMotifLen[t];
******************/

 
                      
/****************OLD **********/
   for(i = 0; i < B->IP->nNumSequences; i++) {
      if(i == 0) start = 0;
      else start = (*B->Seq->nvEndLocs)[i - 1];
      end = (*B->Seq->nvEndLocs)[i] - 1;
      for(t = 0; t < B->IP->nNumMotifTypes; t++) {
         for(tot = 0.0, n = start; n <= end; n++) {
            if(PossFragStartPos(Pos, n, t, B)) {
               (*StartMotifProb)[n][t][FORWARD] = 
		 in_motif_prob(B, P, n, t, FALSE, FALSE,  &dMotifProb, &dBGProb);
               tot += (*StartMotifProb)[n][t][FORWARD];
               if(B->IP->RevComplement) {
                  (*StartMotifProb)[n][t][REVERSE] = 
		    in_motif_prob(B, P, n,t,TRUE, FALSE,  &dMotifProb, &dBGProb);
                  tot += (*StartMotifProb)[n][t][REVERSE];
               }
            }
         }
         if((tot > 0.0) && B->IP->site_samp ) 
            for(n = start; n <= end; n++) 
               if(PossFragStartPos(Pos, n, t, B)) {
                  (*StartMotifProb)[n][t][FORWARD] /= tot;
                  if(B->IP->RevComplement)
                     (*StartMotifProb)[n][t][REVERSE] /= tot;
               }
      }
   }        
         
   for(n = 0; n < B->IP->nSeqLen; n++) {
      for(t = 0; t < B->IP->nNumMotifTypes; t++) 
      {
          nLen = min( (*B->Seq->nvEndLocs)[Pos[t][n].nSeq] - n,  B->IP->nMotifLen[t] ); /* BT 4/21/97 */
          mid = nLen / 2 + nLen % 2 - 1;
          if(PossFragStartPos(Pos, n, t, B)) {
             prob = (*StartMotifProb)[n][t][FORWARD];
             for(j = 1; j <= nLen; j++) 
             {
                index = (int)((*B->Seq->R)[n + j - 1]) - 97;
                if( index < B->IP->nAlphaLen )		/* BT 4/18/97 */
                {
                   if(B->IP->nAlphaLen == 4) {
		     if(j > mid) {
		       palindex = CH2INTCOMP(n+j-1, B);
		       repindex = CH2INT(j-mid-1, B);
		     }
		     else {
		       palindex = index;
		       repindex = index;
		     }
                   }
                   prob = (*StartMotifProb)[n][t][FORWARD];
                   if(B->IP->RevComplement) {
                       if(j == 1) {
                          fwd += (*StartMotifProb)[n][t][FORWARD];
                          rev += (*StartMotifProb)[n][t][REVERSE];
                       }
                       prob *= ratio;
                       if(B->IP->AltModel->Palandromic[t][j])
                          (*expected)[palindex][t][j] += prob;
                       else if(B->IP->AltModel->Repeat[t][j])
                          (*expected)[repindex][t][j] += prob;
                       else
                          (*expected)[index][t][j] += prob;
                       (*expected)[index][TOT_BG][BG] += ratio - prob;
                       index = CH2INTCOMP(n + nLen - j, B);
                       if( j > mid) {
			 palindex = CH2INTCOMP(n + nLen - j, B);
			 repindex = CH2INT(n + j - mid, B);
		       }
                       else {
			 palindex = index;
			 repindex = index;
		       }
                       prob = (*StartMotifProb)[n][t][REVERSE] * (1.0 - ratio);
                       if(B->IP->AltModel->Palandromic[t][j])
                          (*expected)[palindex][t][j] += prob;
                       if(B->IP->AltModel->Repeat[t][j])
			 (*expected)[repindex][t][j] += prob;
                       else
                          (*expected)[index][t][j] += prob;
                       (*expected)[index][TOT_BG][BG] += 1.0 - ratio - prob;
                       if( index >= B->IP->nAlphaLen )
                          p_error( "Invalid Index" );		/* BT 4/21/97 --- Debug */
                   }
                   else {
                      if(B->IP->AltModel->Palandromic[t][j])
                         (*expected)[palindex][t][j] += prob;
                      else
                         (*expected)[index][t][j] += prob;
                      (*expected)[index][TOT_BG][BG] += 1.0 - prob;
                   }
                }
             }
          }
      }
   }
   if(B->IP->RevComplement && B->IP->site_samp) {	/* BT 4/23/97 */
      B->IP->nNumMotifs[0][FORWARD] = (int)fwd;
      B->IP->nNumMotifs[0][REVERSE] = (int)rev;
   }
   for(n = 0; n < B->IP->nAlphaLen; n++) 
      for(t = 0; t < B->IP->nNumMotifTypes; t++) 
         (*expected)[n][TOT_BG][BG] /= B->IP->nMotifLen[t];
/*******************/
}


double log_likelihood(IPtype IP, double ***expected, double *total, 
                      double tot_back)

/*===========================================================================*/
/* FUNCTION NAME : log_likelihood                                            */
/*                                                                           */
/* DESCRIPTION : calculates the log likelihood value for the alignment given */
/*             the expected number of each residue, the total number of      */
/*             residues in the motifs, the expected counts of residues in the*/
/*             background and the total number of background residues        */
/*===========================================================================*/

{
   int j, t=0, n;
   double logl, bg_part;
   double tot_logl;

   for(tot_logl = 0.0, t = 0; t < IP->nNumMotifTypes; t++) {
      for(logl = 0.0, j = 1; j <= IP->nMotifLen[t]; j++) 
         for(n = 0; n < IP->nAlphaLen; n++) 
         {  
            if( expected[n][t][j] != 0.0 )			/* BT 2/14/97 */
               logl +=  (expected[n][t][j] / total[t]) * 
                        log(expected[n][t][j]/total[t]);
         }
      logl *= total[t];
      for(bg_part = 0.0, n = 0; n < IP->nAlphaLen; n++)
      { 
         if( expected[n][TOT_BG][BG] != 0.0 )			/* BT 2/14/97 */
            bg_part += (expected[n][TOT_BG][BG] / tot_back) *
                       log(expected[n][TOT_BG][BG] / tot_back);
     } 
      bg_part *= tot_back;
      logl += bg_part;
      tot_logl += logl;
   }
   return tot_logl;
}

double log_epsilon(double *total, double tot_back, double epsilon)

/*===========================================================================*/
/* FUNCTION NAME : log_epsilon                                               */
/*                                                                           */
/*DESCRIPTION : calculates the log likelihood value that describes how likely*/
/*             the counts are given the total number of motifs and the number*/
/*             of expected motifs                                            */
/*===========================================================================*/

{
   int t = 0;
   double likelihood;

   likelihood = total[t] * log(epsilon);
   likelihood += tot_back * log(1.0 - epsilon);
   return likelihood;
}

double log_cost(Model B, PoSition **Pos, double ***StartMotifProb)

/*===========================================================================*/
/* FUNCTION NAME : log_cost                                                  */
/*                                                                           */
/* DESCRIPTION : calcs the log cost describing how likely the alignment is   */
/*             given the probability of each position being a motif starting */
/*             position                                                      */
/*===========================================================================*/

{
   double cost=0.0, prob;
   int n, t=0;
   double tot_cost;

   for(tot_cost = 0.0, t = 0; t < B->IP->nNumMotifTypes; t++) {
      for(cost = 0.0, n = 0; n < B->IP->nSeqLen; n++) {
         if(PossFragStartPos(Pos, n, t, B)) {
            prob = StartMotifProb[n][t][FORWARD];
            if(B->IP->RevComplement) {
               prob /= 2.0;
               prob += StartMotifProb[n][t][REVERSE]/2.0;
            }
            if( prob != 0.0 ) 			/* BT 2/14/97 */
               cost += prob * log(prob);
            if(cost > 0.0)
               printf("prob = %f cost %f\n", prob, cost);
         }
      }
      tot_cost += cost;
   }
   return cost;
}

void init_total(double **total, double ***expected, IPtype IP)
{
   int t, n;

   for(t = 0; t < IP->nNumMotifTypes; t++) {
      (*total)[t] = 0.0;
      for(n = 0; n < IP->nAlphaLen; n++)  
         (*total)[t] += expected[n][t][1];
      if(IP->RevComplement)
         (*total)[t] *= 2.0;
   }
   
}

void print_EM_results(Model B, PoSition **Pos, double ***StartMotifProb)

   /*=====================================================================*/
   /* FUNCTION NAME : print_EM_results                                    */
   /*=====================================================================*/

{
   double     **dProbArray, dMax, dLocMax, EMmap;
   int        n, t=0, last_seq, site_num, i;
   int        last_inc, num = 0; 
   FILE       *fpt;
   MaxResults M;
   int        nMotifNumSave;		/* BT 4/25/97 */

   init_maxdata(&M);
   copy_counts(B);
   NEWP(dProbArray, B->IP->nSeqLen, double);
   for(i = 0; i < B->IP->nSeqLen; i++) {
      NEW(dProbArray[i], B->IP->nNumMotifTypes, double);
      for(t = 0; t < B->IP->nNumMotifTypes; t++) {
         dProbArray[i][t] = 0.0;
         Pos[t][i].nMotifStartPos = FALSE;
         Pos[t][i].nInMotif = FALSE;			/* BT 2/19/97 */
      }
   }
   t = 0;
   fpt = B->IP->Datafiles->out_fpt;
   last_seq = -1; site_num = 1;
   
   B->IP->nNumMotifs[t][FORWARD] = 0;
   B->IP->nNumMotifs[t][REVERSE] = 0;

   for(n = 0; n < B->IP->nSeqLen; n++) {
      if(PossFragStartPos(Pos, n, t, B)) {
         if(StartMotifProb[n][t][FORWARD] > 0.05 /*B->IP->dCutoff*/) {
               if(B->IP->RevComplement && (StartMotifProb[n][t][REVERSE] > 
                                           StartMotifProb[n][t][FORWARD])) {
                  dProbArray[num][t] = StartMotifProb[n][t][REVERSE];

                  adjust_counts(B,ADD,n,t, TRUE);
		  B->IP->nNumMotifs[t][REVERSE]++;
                  set_in_motif(Pos, n, B, t, TRUE);
		  n+=B->IP->nMotifLen[t]-1;
               }
               else {
		 dProbArray[num][t] = StartMotifProb[n][t][FORWARD];
		 adjust_counts(B,ADD,n,t, FALSE);
		 B->IP->nNumMotifs[t][FORWARD]++;
		 set_in_motif(Pos, n, B, t, FALSE);
		 n+=B->IP->nMotifLen[t]-1;
               }
                 
               num++;
         }
         else if(B->IP->RevComplement &&
                 (StartMotifProb[n][t][REVERSE] > 0.05 /* B->IP->dCutoff*/)){
               dProbArray[num][t] = StartMotifProb[n][t][REVERSE];
               adjust_counts(B,ADD,n,t,TRUE);
	       B->IP->nNumMotifs[t][REVERSE]++;
               set_in_motif(Pos, n, B, t, TRUE);
	       n+=B->IP->nMotifLen[t]-1;
               num++;
         }
      }
   }
   EMmap = CalcMapProb(B, B->IP->is_defined[cl_R]);   
   M = setMaxData(B->F, B->IP, EMmap, 1, Pos, &last_inc, &dMax,
                  &dLocMax, M, dProbArray, B);
   M.frequency = NULL;

   nMotifNumSave  = B->IP->nNumMotifTypes;		/* BT 4/25/97 */
   B->IP->nNumMotifTypes = 1;				/* Only the 1st motif type counts for EM */
   print_info(B, M, TRUE, OTHER );
   B->IP->nNumMotifTypes = nMotifNumSave;

   free_maxdata(&M, B->IP);
   for(i = 0; i < B->IP->nSeqLen; i++) 
      free(dProbArray[i]);
   free(dProbArray);
}

void init_stats(EMStats *S)
{
   S->logl = S->likelihood = S->cost = S->new_tot = S->abs_tot =
   S->pct_tot = S->abs_logl = S->pct_logl = S->last_tot = S->last_logl = 0.0;
}

void calculate_EM_stats(EMStats *S, Model B, double ***expected,
                        double *total, double tot_back, double epsilon,
                        PoSition **Pos, double ***StartMotifProb)

   /*=================================================================*/
   /* FUNCTION NAME : calculate_EM_stats                              */
   /* DESCRIPTION   : Calculates the Expectation/Maximization stats   */
   /*                 including the log likelihood and cost ratios    */
   /*=================================================================*/

{
   S->logl = log_likelihood(B->IP, expected, total, tot_back);
   S->likelihood = log_epsilon(total, tot_back, epsilon);
   S->cost = log_cost(B, Pos, StartMotifProb);
   S->new_tot = S->logl - S->cost;
   if(!B->IP->site_samp)
      S->new_tot += S->likelihood;
   S->abs_tot = S->new_tot - S->last_tot;
   if(fabs(S->last_tot) > 0.001)
      S->pct_tot = S->abs_tot / S->last_tot * 100.0;
   S->abs_logl = S->logl - S->last_logl;
   if(fabs(S->last_logl) > 0.001)
      S->pct_logl = S->abs_logl / S->last_logl * 100.0;
   S->last_tot = S->new_tot;
   S->last_logl = S->logl;
}

void print_EM_stats(FILE *fpt, EMStats S)

   /*==================================================================*/
   /* FUNCTION NAME : print_EM_stats                                   */
   /* DESCRIPTION   : Prints the statistics for the Expect./Max.       */
   /*                 algorithm                                        */
   /*==================================================================*/

{
   fprintf(fpt,"%10.2f %8.2f %6.2f %10.2f %10.2f"
             , S.logl, S.likelihood, S.cost, S.new_tot, S.abs_tot);
   if(fabs(S.pct_tot) > 0.0001)
      fprintf(fpt, "%6.1f ", S.pct_tot);
   else
      fprintf(fpt, "   0.0 ");
   fprintf(fpt, "%10.2f ", S.abs_logl);
   if(fabs(S.pct_logl) > 0.001)
      fprintf(fpt, "%6.1f\n", S.pct_logl);
   else
      fprintf(fpt, "   0.0\n");

}
void EM(Model B, PoSition **Pos, MaxResults optmax)
{
   /*==================================================================*/
   /* FUNCTION NAME : EM                                               */
   /* DESCRIPTION : Performs Expectation/Maximization for the data     */
   /*               given the alignment found in near optimal sampling */
   /*==================================================================*/

   int t, j, n, q, i;
   EMStats S;
   ProbStruct P;
   double ***StartMotifProb, *total, ***expected;
   double prob, tot_back, epsilon, dmodel_pseudo;
   double old_avg, new_avg, diff, high, low;
   double recent[6];					/* BT 2/14/97 */
   
   init_stats(&S);
   /* set_values(B, optmax);  */			/* BT 2/24/97 */
   print_header(B->IP->Datafiles->out_fpt);

   NEW(total, B->IP->nNumMotifTypes, double);
   NEWPP(expected, B->IP->nAlphaLen + 1, double);	/* BT 4/18/97 */
   for(n = 0; n < B->IP->nAlphaLen + 1; n++) 		/* BT 4/18/97 */
   {
      NEWP(expected[n], B->IP->nNumMotifTypes, double);
      for(t = 0; t < B->IP->nNumMotifTypes; t++) 
         NEW(expected[n][t], (*B->IP).nMotifLen[t]+1, double);  
   }

   NEWPP(StartMotifProb, B->IP->nSeqLen, double);
   for(n = 0; n < B->IP->nSeqLen; n++) {
      NEWP(StartMotifProb[n], B->IP->nNumMotifTypes, double);
      for(t = 0; t < B->IP->nNumMotifTypes; t++)
         NEW(StartMotifProb[n][t], 2, double);
   }

   init_prob(&P, B->IP);
   update_prob(&P, B, TRUE);
   update_posterior_prob(B); 

   init_expected(B, &expected, &StartMotifProb, Pos, P);
   init_total(&total, expected, B->IP);
   t = 0;
   tot_back = B->IP->nSeqLen - total[t] * (*B->IP).nMotifLen[t];
   epsilon = (*B->IP).dposterior_prob[t];
   calculate_EM_stats(&S, B, expected, total, tot_back, epsilon, Pos,
                      StartMotifProb);
   print_EM_stats(B->IP->Datafiles->out_fpt, S);
   dmodel_pseudo = (*B->C).dmodel_pseudo[t][FORWARD];
   if(B->IP->RevComplement)
      dmodel_pseudo += (*B->C).dmodel_pseudo[t][REVERSE];
   S.pct_logl = 2.0;
   q = 0;
   old_avg = new_avg = S.new_tot;
   diff = 1.0;
/**/
   while( ((fabs(diff) > 0.1) || q < 50) &&  (q < B->IP->nMaxIterations) ) { /* BT 2/19/97 */
/**/
/**   while(q < 10) { **/
      recent[q % 6] = S.new_tot;
      old_avg = new_avg;
      q++;
      for(tot_back = 0.0, n = 0; n < B->IP->nAlphaLen; n++)
         tot_back += expected[n][0][BG];
      t = 0;
      epsilon = (total[t] + dmodel_pseudo) /
                ((*B->C).dtot_sites[t] + (*B->C).dtot_pseudo[t] - 1.0);
      (*B->IP).dposterior_prob[t] = epsilon;
      for(n = 0; n < B->IP->nAlphaLen; n++) {
         for(t = 0; t < B->IP->nNumMotifTypes; t++) {
            P.dvInBGProb[t][n] = expected[n][TOT_BG][BG] / tot_back;  /* BT 12/21/99 */
            for(j = 1; j <= (*B->IP).nMotifLen[t]; j++) {
               prob = expected[n][t][j];
               P.dvInMotifProb[t][j-1][n] = prob / total[t];   /* BT 12/21/99 */
            }
         }
      }
      init_expected(B, &expected, &StartMotifProb, Pos, P);
      for(tot_back = 0.0, n = 0; n < B->IP->nAlphaLen; n++)
         tot_back += expected[n][0][BG];

      init_total(&total, expected, B->IP);
      calculate_EM_stats(&S, B, expected, total, tot_back, epsilon, Pos, 
                         StartMotifProb);
      print_EM_stats(B->IP->Datafiles->out_fpt, S);      
      new_avg = (old_avg * (double)q + S.new_tot) / (double)(q + 1.0);
      if(q > 6) {
         high = low = recent[0];
         for(i = 1; i < 6; i++) {
            if(recent[i] > high)
               high = recent[i];
            else if(recent[i] < low)
               low = recent[i];
         }
         diff = high - low;
      }
   }
   free_prob(&P, B->IP);
   free(total);
   for(n = 0; n < B->IP->nAlphaLen + 1; n++) 		/* BT 4/18/97 */
   {
      for(t = 0; t < B->IP->nNumMotifTypes; t++)
         free(expected[n][t]);
      free(expected[n]);
   }
   free(expected);
   B->IP->Logl = S.logl;
   print_EM_results(B, Pos, StartMotifProb);
   for(n = 0; n < B->IP->nSeqLen; n++) {
      for(t = 0; t < B->IP->nNumMotifTypes; t++) 
          free(StartMotifProb[n][t]);
      free(StartMotifProb[n]);
   }
   free(StartMotifProb);
}

