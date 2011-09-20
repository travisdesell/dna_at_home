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
/* $Id: printdata.c,v 1.10 2007/05/23 18:19:57 Bill Exp $         */
/*                                                                        */
/* Author :       Eric C. Rouchka July 17, 1996                           */
/*                Jun Zhu October 20, 1996                                */
/*                Bill Thompson  2/3/97                                   */
/* Description :  This file contains the functions that are used in       */
/*                printing the alignment data                             */
/**************************************************************************/

/*------------------------------------------------------------------------*/
/*                           INCLUDE FILES                                */
/*------------------------------------------------------------------------*/
#include "printdata.h"


/*------------------------------------------------------------------------*/
/*                        LOCAL FUNCTION PROTOTYPES                       */
/*------------------------------------------------------------------------*/
void print_strings(Model B, MaxResults maxData, int t, int bPrintOpt, int bPrintType );
void Calculate_information(Model B, IPtype IP, Ctype C, int t, int nNumMotifs,
                           FILE *fpt, double *observed, FILE *prout_fpt);       /* BT 9/26/97 */
void print_probability(Model B, IPtype IP, Ctype C, int t, int nNumMotifs,
		       FILE *fpt );                /* BT 1/5/97 */
void print_freq_results(Model B, FILE *fpt, MaxResults maxData, int t, FILE *prout_fpt );
void PrintSeqData( Model B, MaxResults maxData  );
void PrintSeqFreqData( Model B, MaxResults maxData  );
int  MaxPrintLenMotif( Model B );

/*------------------------------------------------------------------------*/
/*                          FUNCTION DECLARATIONS                         */
/*------------------------------------------------------------------------*/

void print_shift(int col_shift)

   /*=========================================================*/
   /* FUNCTION NAME : print_shift                             */
   /* DESCRIPTION   : Prints the number of columns that shift */
   /*=========================================================*/
{
   int i;
   
   printf("[");
   if(col_shift < 0) for(i = col_shift; i < 0; i++)   printf("-");
   else              for(i = 0; i < col_shift; i++)   printf("+");
   printf("] ");
}

void print_maxData(int nNumMotifTypes, MaxResults maxData)

   /*=========================================================*/
   /* FUNCTION NAME : print_maxData                           */
   /* DESCRIPTION   :                                         */
   /*=========================================================*/

{
   int t;

   printf("\nMAX :: %f (Seed = %ld, Iteration = %d   ",
          maxData.dProbability, maxData.nseed, maxData.nIterationNum);
   if( maxData.nNumMotifs )
     {
       for(t = 0; t < nNumMotifTypes; t++) 
	 printf("Motif %c = %d ", (char)((int)'A' + t), 		
              maxData.nNumMotifs[t]);
     }
   printf(")\n");
}


void put_prior(Model B)

   /*===========================================================*/
   /* FUNCTION NAME : put_prior                                 */
   /* DESCRIPTION   :                                           */
   /*===========================================================*/
 
{
   double variance, stdev, model, bg;
   int t, nNumMotifs, sites;
  
   for(t = 0; t < B->IP->nNumMotifTypes; t++) {
      nNumMotifs = NUMMOTIFS(B->IP->nNumMotifs[t]);
      model = B->C->dmodel_pseudo[t][FORWARD];
      bg = B->C->dmodel_pseudo[t][FORWARD];
      if(B->IP->RevComplement) {
         model += B->C->dmodel_pseudo[t][REVERSE];
         bg += B->C->dmodel_pseudo[t][REVERSE];
      }
      bg = B->C->dbg_pseudo[t] - model * B->IP->nMotifLen[t];
      sites = B->IP->nPossSites[t] - nNumMotifs * B->IP->nMotifLen[t]; 
      variance = (model + (double)nNumMotifs);
      variance *= (bg + sites - (double)nNumMotifs);
      variance /= (double)(B->C->dTot[t] + sites)*(B->C->dTot[t] + sites) *
                  (B->C->dTot[t] + sites + 1.0);
      stdev = sqrt(variance) * 2.0 * sites;
      printf("motif %c: %d (+/- %.2f) out of %d", (char) t + 'A', nNumMotifs, stdev, sites);
      printf("   a = %g; b = %g; p = %g\n", model, bg, model/B->C->dtot_pseudo[t]);
   }
}

void print_info(Model B, MaxResults maxData, int bPrintOpt, int bPrintType )

   /*==================================================================*/
   /* FUNCTION NAME : print_info                                       */
   /*==================================================================*/
{
   int    t;
   FILE   *fpt;
   int    oldLen;
   double dMap;
   double freq_map;

   fpt = B->IP->Datafiles->out_fpt;

   if( B->IP->nNumMotifTypes > 1 && bPrintOpt )
     {
       fprintf( fpt, "\n\n" );
       fprintf( fpt, "=============================================================\n" );
       fprintf( fpt, "======                Results by Sequence               =====\n" );
       fprintf( fpt, "=============================================================\n" );
       fprintf( fpt, "\n\n" );
       PrintSeqData( B, maxData ); 

       fprintf(fpt, "\n\n");
       fprintf(fpt, "\n\n");
       fprintf(fpt, "Column 1 :  Sequence Number, Site Number\n");
       fprintf(fpt, "Column 2 :  Motif type\n");
       fprintf(fpt, "Column 3 :  Left End Location\n");
       fprintf(fpt, "Column 4 :  Motif Element\n");
       fprintf(fpt, "Column 5 :  Right End Location\n");
       fprintf(fpt, "Column 6 :  Probability of Element\n");
       fprintf(fpt, "Column 7 :  Forward Motif (F) or Reverse Complement (R) \n");
       fprintf(fpt, "Column 8 :  Sequence Description from Fast A input\n\n");

#ifdef _MPI_
       fsync( fileno( fpt ) );
#endif
     }

   for(t = 0; t < B->IP->nNumMotifTypes; t++) 
     {
       oldLen = B->IP->nMotifLen[t];
       B->IP->nMotifLen[t] = maxData.nMotifLen[t];    /* set motif length to best length */ /* BT 5/28/97 */

       fprintf( fpt, 
		"-------------------------------------------------------------------------\n" );  /* BT 4/9/97 */
       fprintf(fpt, "                          MOTIF %c\n\n", (char)(97 + t));
      
       /* if(B->IP->nNumMotifs[t][FORWARD] != 0 || B->IP->nNumMotifs[t][REVERSE] != 0) */  /* BT 5/28/ 97 */
       if( maxData.nNumMotifs[t] ) 
	   print_strings(B, maxData, t, bPrintOpt, bPrintType ); 
       else
	 fprintf(fpt, "               No Motifs Detected\n\n");
       B->IP->nMotifLen[t] = oldLen;    /* restore old length */ /* BT 5/28/97 */
#ifdef _MPI_
       fsync( fileno( fpt ) );
#endif
     }

   if( bPrintOpt && ! B->IP->inCentroidAlign )
     {
       dMap = maxData.dTotalMap;
       /*  fprintf(fpt, "Difference of Logs of Maps = %.5f\n", dMap );  */
       
       fprintf(fpt, "\n" );         
       fprintf(fpt, "Log Background portion of Map = %.5f\n", maxData.dBkgndMap );         
       fprintf(fpt, "Log Alignment portion of Map = %.5f\n", maxData.dBetaMap );         
       fprintf(fpt, "Log Site/seq portion of Map = %.5f\n", maxData.dSeqMap );         
       fprintf(fpt, "Log Null Map = %.5f\n", B->IP->dnull_map );      
       if( B->IP->is_defined[cl_R] || B->IP->is_defined[cl_c] )
	 fprintf(fpt, "Non Palindromic Map = %.5f\n", maxData.dTotalNonPalinMap );
       if( B->Phylo->phyloTree )	
	 fprintf(fpt, "Non Phylogenetic Map = %.5f\n", maxData.dTotalNonPhyloMap );
       fprintf(fpt, "Log Map = %.5f\n", dMap );
     }

   if( B->IP->is_defined[cl_m] && maxData.frequency != NULL  && ! B->IP->inCentroidAlign )
     {
       copy_counts(B);
       freq_map = set_freq_map(B, maxData);
       if( ! B->IP->is_defined[cl_y] )
	 fprintf(fpt, "\nFrequency Map = %f\n", freq_map);
       fprintf(fpt, "Nearopt Map   = %f\n\n", maxData.dTotalMap);   
     }

   if( bPrintOpt && ! B->IP->inCentroidAlign  )
     fprintf( fpt, 
	      "\n\nlog MAP = sum of motif and fragmentation parts of MAP + background + alignment + sites/seq - Null\n\n" );

#ifdef _MPI_
   fsync( fileno( fpt ) );
#endif

   if( B->IP->nNumMotifTypes > 1 && maxData.frequency != NULL && ! B->IP->is_defined[cl_y] )
     {
       fprintf( fpt, "\n\n" );
       fprintf( fpt, "=============================================================\n" );
       fprintf( fpt, "======                Results by Sequence               =====\n" );
       fprintf( fpt, "====== ELEMENTS OCCURRING GREATER THAN %3d%% OF THE TIME =====\n", 
		(int)(100 * B->IP->dNearOptCutoff ));
       fprintf( fpt, "=============================================================\n" );
       fprintf( fpt, "\n\n" );
       PrintSeqFreqData( B, maxData );     

       fprintf(fpt, "\n\n");
       fprintf(fpt, "Column 1 :  Sequence Number, Site Number\n");
       fprintf(fpt, "Column 2 :  Motif type\n");
       fprintf(fpt, "Column 3 :  Left End Location\n");
       fprintf(fpt, "Column 4 :  Motif Element\n");
       fprintf(fpt, "Column 5 :  Right End Location\n");
       fprintf(fpt, "Column 6 :  Probability of Element\n");
       fprintf(fpt, "Column 7 :  Forward Motif (F) or Reverse Complement (R) \n");
       fprintf(fpt, "Column 8 :  Sequence Description from Fast A input\n\n");
     }
#ifdef _MPI_
   fsync( fileno( fpt ) );
#endif
}

void print_counts(Model B, Ctype C, IPtype IP, int t, int nNumMotifs, FILE *prout_fpt)

   /*====================================================================*/
   /* FUNCTION NAME : print_counts                                       */
   /*                                                                    */
   /* DESCRIPTION : This function prints out the residue frequencies for */
   /*               each of the amino acids                              */
   /*====================================================================*/

{
   int i, j, pct;
   FILE *fpt;
   double *observed, total;

   NEW(observed, IP->nAlphaLen + 1, double);
   observed[0] = 0.0;
   for(j = 0; j < IP->nAlphaLen; j++) 
     {
       observed[j+1] = (double)C->fCounts[t][BG][j];
       for(i = 1; i <= IP->nMotifLen[t]; i++) 
         observed[j+1] += (double)C->fCounts[t][i][j]; 
       observed[0] += observed[j+1];
     }
   fpt = IP->Datafiles->out_fpt;
   fprintf(fpt, "\n\nMotif model (residue frequency x 100)\n");
   fprintf(fpt, "____________________________________________\n");
   fprintf(fpt, "Pos. #  ");
   if(IP->nAlphaLen == 2)
     fprintf(fpt, "a-t  c-g Info");
   else 
     for(j = 0; j < IP->nAlphaLen; j++)
       fprintf(fpt, "%4c", curr_ch((char)(j + 97), IP->nAlphaLen, FALSE));
   fprintf(fpt, "  Info\n_____________");
   for(j = 0; j < IP->nAlphaLen; j++) fprintf(fpt, "____"); 
   fprintf(fpt, "\n"); 
   Calculate_information(B, IP, C, t, nNumMotifs, fpt, observed, prout_fpt);   /* BT 9/26/97 */
   fprintf(fpt, "\nnonsite ");
   for(j = 0; j < IP->nAlphaLen; j++) 
     {
       if((observed[0] == 0) || (observed[j + 1] == 0))
         pct = 0;
       else
         pct = (int)(100.0 * observed[j+1] / observed[0]);
       if(pct > 0) fprintf(fpt, "%4d", pct);
       else        fprintf(fpt, "  . ");
     }
   fprintf(fpt, "\nsite    ");
   for(total = 0, j = 0; j < IP->nAlphaLen; j++) 
     total += (observed[j+1] - (double)C->fCounts[t][BG][j]);

   for(j = 0; j < IP->nAlphaLen; j++) 
     {
       if(total > 0.0)
         pct = (int)(100.0 * (observed[j+1] - (double)C->fCounts[t][BG][j]) /
		     total);
       else
         pct = 0;
       if(pct > 0)  fprintf(fpt, "%4d", pct);
       else         fprintf(fpt, "  . ");       
   } 

   print_probability(B, IP, C, t, nNumMotifs, fpt);   /* BT 1/5/97 */

   fprintf(fpt, "\n\n%d columns\nNum Motifs: %d\n",IP->nMotifLen[t],nNumMotifs);
   free(observed);
}


void Calculate_information(Model B, IPtype IP, Ctype C, int t, int nNumMotifs,
                           FILE *fpt, double *observed, FILE *prout_fpt)
 
   /*====================================================================*/
   /* FUNCTION NAME : Calculate_information                              */
   /* DESCRIPTION   : Calculates the information per parameter info      */
   /*                 in log base 2 for ease in telling information bits */
   /*====================================================================*/

{
   int     i, j, k, pct;
   double  info, p, q, ii;
   int     row = 0; 
   int     first = 0;
   short   bGap;			/* BT 2/5/97 */ 
   double  dPseudo;

   if( ! IP->is_defined[cl_F] )			/* BT 2/12/97 */
	first = FirstCol( B->F, t ); 

   if( prout_fpt != NULL )                       
     fprintf(prout_fpt, ">PSEUDO %d\n", t+1 ); 

   for(i = 1; i <= IP->nMotifLen[t]; i++) 
     {
       bGap = FALSE; 						/* BT 2/5/97 */
       if( ! IP->is_defined[cl_F] )
	 { 
	   while( B->F->nColMask[t][first + row] != COL_ON )		/* BT 2/3/97 */
	     { 
	       row++;
	       bGap = TRUE;
	     }
	 } 
       if( bGap )
         fprintf( fpt, "\n" );   				/* BT 2/5/97 */
      fprintf(fpt, "%4d |  ", row + 1);
      row++;
/*      fprintf(fpt, "%4d |  ", i);  */

      for(dPseudo = 0.0, j = 0; j < IP->nAlphaLen; j++) 
	{     
          dPseudo += C->dPseudoCounts[t][i][j];
	} 

      
      for(info = 0.0, j = 0; j < IP->nAlphaLen; j++)  
	{
	  pct = (int)(100.0 * (double)C->fCounts[t][i][j] /(double)nNumMotifs);
	  /*	  p = ((double)C->fCounts[t][i][j]+C->dPseudoCounts[t][i][j]) / 
	    ((double)nNumMotifs + C->dSumBGPseudo);  */
	  p = ((double)C->fCounts[t][i][j]+C->dPseudoCounts[t][i][j]) /    /* BT 12/22/97 */
	    ((double)nNumMotifs + dPseudo); 
	  
	  if(p > 0.0) 
	    {
	      q = observed[j + 1] / observed[0]; 
	      ii = p*log(p/q)/log(2.0);
	      info += ii;
	    }
	  if(pct > 0)   fprintf(fpt, "%4d", pct);
	  else          fprintf(fpt, "  . ");

	  if( prout_fpt != NULL )
	    fprintf(prout_fpt, "%4.1f ", (double) pct / 10.0);
	}
      fprintf(fpt, "  %3.1f\n", info);
      if( prout_fpt != NULL )
	fprintf(prout_fpt, "\n");
     }
   
   if( prout_fpt != NULL )
     {
       fprintf(prout_fpt, ">\n\n! Counts -- remove this line for dscan\n!");

       if( ! B->IP->is_defined[cl_F] )
	 {
	   for(i = 0; i < B->F->FragWidth[t]; i++)
	     {
	       if(B->F->nColMask[t][i + first] == COL_ON)  
		 {
		   fprintf(prout_fpt, "*");
		 }
	       else
		 {
		   fprintf(prout_fpt, ".");
		 }
	     }
	   fprintf(prout_fpt, "\n");
	   
	   for( k = 1, i = 0; i < B->F->FragWidth[t]; i++ ) 
	     {
	       fprintf(prout_fpt, "! ");
	       if(B->F->nColMask[t][i + first] == COL_ON)  
		 {		   
		   for( j = 0; j < IP->nAlphaLen; j++ )  
		     {
		       fprintf( prout_fpt, "%4d", (int) C->fCounts[t][k][j] );	   
		     }
		   k++;
		 }
	       else
		 {
		   for( j = 0; j < IP->nAlphaLen; j++ )  
		     {
		       fprintf( prout_fpt, "%4d", 0 );	   
		     }
		 }
	       fprintf( prout_fpt, "\n" );
	     }
	   fprintf(prout_fpt, "\n");       
	 }
       else
	 {
	   for( i = 0; i < B->IP->nMotifLen[t]; i++)
	     {
	       fprintf(prout_fpt, "*");
	     }	   
	   fprintf(prout_fpt, "\n");
	   	   
	   for(i = 1; i <= IP->nMotifLen[t]; i++) 
	     {
	       fprintf(prout_fpt, "! ");
	       for( j = 0; j < IP->nAlphaLen; j++ )  
		 {
		   fprintf( prout_fpt, "%4d", (int) C->fCounts[t][i][j] );	   
		 }
	       fprintf( prout_fpt, "\n" );
	     }
	   fprintf(prout_fpt, "\n");       
	 }
     }
}


void print_probability(Model B, IPtype IP, Ctype C, int t, int nNumMotifs,
		       FILE *fpt )                /* BT 1/5/97 */
   /*====================================================================*/
   /* FUNCTION NAME : print_probability                                  */
   /* DESCRIPTION   : Prints the probability model                       */ 
   /*====================================================================*/
{
   int     i, j;
   double  p;
   int     row = 0; 
   int     first = 0;
   short   bGap;			/* BT 2/5/97 */ 
   double  dPseudo;
   ProbStruct P;

   P.update = TRUE;
   init_prob(&P, IP);
   update_prob(&P, B, TRUE);
   
   fprintf(fpt, "\n\nMotif probability model\n");
   fprintf(fpt, "____________________________________________\n");
   fprintf(fpt, "Pos. #  ");
   for(j = 0; j < IP->nAlphaLen; j++)
     fprintf(fpt, "  %c   ", curr_ch((char)(j + 97), IP->nAlphaLen, FALSE));
   fprintf(fpt, "\n");
   fprintf(fpt, "____________________________________________\n");

   if( ! IP->is_defined[cl_F] )			/* BT 2/12/97 */
	first = FirstCol( B->F, t ); 

   for(i = 1; i <= IP->nMotifLen[t]; i++) 
     {
       bGap = FALSE; 						/* BT 2/5/97 */
       if( ! IP->is_defined[cl_F] )
	 { 
	   while( B->F->nColMask[t][first + row] != COL_ON )		/* BT 2/3/97 */
	     { 
	       row++;
	       bGap = TRUE;
	     }
	 } 
       if( bGap )
         fprintf( fpt, "\n" );   				/* BT 2/5/97 */
      fprintf(fpt, "%4d |  ", row + 1);
      row++;

      for(dPseudo = 0.0, j = 0; j < IP->nAlphaLen; j++) 
      {     
          dPseudo += C->dPseudoCounts[t][i][j];
      } 
      
      for(j = 0; j < IP->nAlphaLen; j++)  
	{
	  p = ((double)C->fCounts[t][i][j]+C->dPseudoCounts[t][i][j]) /    /* BT 12/22/97 */
	    ((double)nNumMotifs + dPseudo); 
	  fprintf(fpt, "%4.3f ", p);
	}
      fprintf(fpt, "\n");
     }
   fprintf(fpt, "\n");

   fprintf(fpt, "\n\nBackground probability model\n");
   fprintf(fpt, "        " );
   for(j = 0; j < IP->nAlphaLen; j++)  
     fprintf(fpt, "%4.3f ", P.dvInBGProb[t][j]);   /* BT 12/21/99 */
   fprintf(fpt, "\n\n");

   free_prob(&P, IP);
}


void print_strings(Model B, MaxResults maxData, int t, int bPrintOpt, int bPrintType )

   /*==============================================================*/
   /* FUNCTION NAME : print_strings                                */
   /* DESCRIPTION : This function prints out the maximal alignment */
   /*               of motif elements that the program found       */
   /****************************************************************/

{
   int    i, last_seq=-1, site_num=1;
   FILE  *fpt;
   FILE  *sn_fpt = NULL;    /* BT 9/17/97 */
   FILE  *prout_fpt = NULL;    /* BT 9/26/97 */
   int   j;
   int   first;

   if(B->IP->is_defined[cl_o])  fpt = B->IP->Datafiles->out_fpt;
   else                         fpt = stdout;
      
   fprintf(fpt, "\n\n");
   if( maxData.nNumMotifs[t] > 0 )    /* BT 5/7/97 */
     {
       if( maxData.frequency != NULL && 
	   (B->IP->is_defined[cl_N] || B->IP->is_defined[cl_O]) )   /* BT 9/17/97 */
	 {
	   sn_fpt = B->IP->Datafiles->sn_fpt;
	   prout_fpt = B->IP->Datafiles->prout_fpt;

	   if( B->IP->is_defined[cl_O] )
	     {
	       fprintf( prout_fpt, "! Pseudocounts are generated so that with the default pseudocount weight of 0.1\n" );
	       fprintf( prout_fpt, "! there will be 1 pseudocount per position. To increase or decrease,\n" );
	       fprintf( prout_fpt, "! use the -w option.\n\n" );
	     }

	   if(B->IP->is_defined[cl_F]) 
	     {
	       if( B->IP->is_defined[cl_O] )
		 {
		   fprintf( prout_fpt, "! To use this file as input to dscan, delete this comment,\n" );
		   fprintf( prout_fpt, "! delete the ! from the line of *'s and the counts\n" );
		   fprintf( prout_fpt, "! and delete the pseudocounts.\n\n" );
		 }

	       if( B->IP->is_defined[cl_N] )
		 {
		   for(j = 0; j < B->IP->nMotifLen[t]; j++)
		     {
		       fprintf(sn_fpt, "*");
		     }
		 }
	     }
	   else 
	     {
	       if( B->IP->is_defined[cl_O] )
		 {
		   fprintf( prout_fpt, "! To use this file as input to dscan, delete this comment\n" );
		   fprintf( prout_fpt, "! delete the ! from the line of *'s and the counts\n" );
		   fprintf( prout_fpt, "! and delete the >PSEUDO and the mask.\n\n" );
		 }

	       first = FirstCol(B->F, t);
	       if( B->IP->is_defined[cl_O] )
		 {
		   fprintf( prout_fpt, ">MASK %d\n", t + 1 );
		   for( i = 0; i < B->F->nMaxLen[t]; i++ )
		     {
		       if( B->F->nColMask[t][i] == 2 || B->F->nColMask[t][i] == 3 )		   
			 fprintf( prout_fpt, "0" );
		       else
			 fprintf( prout_fpt, "%d", B->F->nColMask[t][i] );
		     }
		   fprintf( prout_fpt, "\n>\n\n" );
		 }

	       if( B->IP->is_defined[cl_N] )
		 {
		   for(i = 0; i < B->F->FragWidth[t]; i++)
		     {
		       if(B->F->nColMask[t][i + first] == COL_ON)  
			 {
			   fprintf(sn_fpt, "*");
			 }
		       else
			 {
			   fprintf(sn_fpt, ".");
			 }
		     }
		 }
	     }
	   if( B->IP->is_defined[cl_N] )
	     fprintf(sn_fpt, "\n");
	   if( B->IP->is_defined[cl_O] )
	     fprintf(prout_fpt, "\n");	     
	 }
       
       if( bPrintOpt )
	 {
	   print_counts(B, B->C, B->IP, t, maxData.nNumMotifs[t], NULL);

	   if( B->IP->is_defined[cl_N] && 
	       ((B->IP->is_defined[cl_m] && bPrintType == NEAROPT) || 
		(! B->IP->is_defined[cl_m] && bPrintType == MAPMAX)) )  
	     sn_fpt = B->IP->Datafiles->sn_fpt;
	   else
	     sn_fpt = NULL;
	     
	   for(i = 0; i < maxData.nNumMotifs[t]; i++) 
	     print_it(B, t, maxData.nMotifLoc[i][t], fpt, maxData.dvMotifProb[i][t], 
		      maxData.RevComp[i][t], &last_seq, &site_num, sn_fpt, 0);
	   
	   print_used_col(B, fpt, t);
	   print_col_desc(fpt); 
	   SetPossibleSites( B, t );

	   if( ! B->IP->inCentroidAlign )
	     {
	       fprintf( fpt, "Log Motif portion of MAP for motif %c = %.5f\n", (char)(97 + t), 
			maxData.dMap[t] );  /* BT 3/4/98 */
	       fprintf( fpt, "Log Fragmentation portion of MAP for motif %c = %.5f\n\n", (char)(97 + t), 
			maxData.dFragMap[t] );  /* BT 3/4/98 */
	     }
	 }
       
       if( B->IP->is_defined[cl_wilcox] )
	 WriteWilcoxon( B, maxData, t ); 

       if(maxData.frequency != NULL)
	 {
	   if( !B->IP->is_defined[cl_y] )
	     print_freq_results(B, fpt, maxData, t, prout_fpt);  
	 }
     }
   else
     fprintf( fpt, "No Motifs Detected.\n" );
}

void print_freq_results(Model B, FILE *fpt, MaxResults maxData, int t, FILE *prout_fpt )
{
   int        i, j, nNumMotifs, flag, first = -1, last;
   int        curr_fwd, curr_back, last_fwd = -1, last_back = -1;
   int        curr, index, last_seq, site_num;
   FILE       *sn_fpt = NULL;         /* BT 9/26/97 */
   short      foundIt = FALSE;
   int        frag_width = 0;

   /* BT 4/7/97 */
   fprintf( fpt, "\n\n=============================================================\n" );
   fprintf( fpt, "====== ELEMENTS OCCURRING GREATER THAN %3d%% OF THE TIME =====\n", 
	    (int)(100 * B->IP->dNearOptCutoff ));
   fprintf( fpt, "======                    Motif %c                       =====\n", 
	    (char)(97 + t));      /* BT 7/16/97 */
   fprintf( fpt, "=============================================================\n" );

   fprintf(fpt, "\n\nListing of those elements occurring greater than ");
   fprintf(fpt, "%d%% of the time\n", (int)(100 * B->IP->dNearOptCutoff));
   fprintf(fpt, "in near optimal sampling using %d iterations\n\n", 
                 B->IP->nMaxIterations);

   for(j = 0; j < B->IP->nAlphaLen; j++) 
     {
     for(i = 1; i <= B->IP->nMotifLen[t]; i++) 
	{
	  B->C->fCounts[t][BG][j] += B->C->fCounts[t][i][j];
	  B->C->fCounts[t][i][j] = 0;
	}
     }

      /*-- Update the counts from the near optimal sampling --*/
   nNumMotifs = 0;

   for(i = 0; i < B->IP->nSeqLen; i++) 
     {
       if(maxData.frequency[i][t].frequency > B->IP->dNearOptCutoff ) 
	 {
	   if(!B->IP->is_defined[cl_F]) 
	     {
	       first = FirstCol(B->F, t);
	       last = last_fwd = -1;
	       last_back = B->F->nMaxLen[t];
	       frag_width = LastCol(B->F, t) - first + 1;
	     }

	   for(j = 0; j < B->IP->nMotifLen[t]; j++) 
	     {
	       if(!B->IP->is_defined[cl_F]) 
		 {
		   curr_fwd = NextCol(B->F, t, last_fwd);
		   curr_back = PrevCol(B->F, t, last_back);
		   last_fwd = curr_fwd; last_back = curr_back;
		   curr_fwd -= first; curr_back -= first;

		   curr_back = frag_width - curr_fwd - 1;  /* BT 08/11/99 */

		   if((curr_fwd < 0) || (curr_fwd > B->F->nMaxLen[t])) 
		     p_internal_error("INVALID FORWARD COLUMN");
		   if((curr_back < 0) || (curr_back > B->F->nMaxLen[t])) 
		     p_internal_error("INVALID BACKWARD COLUMN");
		 }
	       else 
		 {
		   curr_fwd = j;
		   curr_back = B->IP->nMotifLen[t] - j - 1;
		 }
             
	       if(maxData.frequency[i][t].fwd_freq > 
		  maxData.frequency[i][t].rev_freq) 
		 {
		   curr = curr_fwd;
		   index = CH2INT(i+curr, B);
		 }
	       else 
		 {
		   curr = curr_back;
		   index = CH2INTCOMP(i + curr_back, B);
		 }
	       B->C->fCounts[t][BG][index]-=1.0;
	       B->C->fCounts[t][j+1][index]+=1.0;
	     }
	   nNumMotifs++;
	 }
     }

   if(nNumMotifs != 0)
     {
       print_counts(B, B->C, B->IP, t, nNumMotifs, prout_fpt);
     }
   else
      fprintf( fpt, "No Motifs Detected.\n");
   last_seq = -1; site_num = 1;
   for(i = 0; i < B->IP->nSeqLen; i++) 
     {
       if(maxData.frequency[i][t].frequency > B->IP->dNearOptCutoff )  
	 {
	   if(maxData.frequency[i][t].fwd_freq > 
	      maxData.frequency[i][t].rev_freq)
	     flag = FALSE;
	   else
	     flag = TRUE;

	   print_it(B, t, i, fpt, maxData.frequency[i][t].frequency,
                    flag, &last_seq, &site_num, sn_fpt, 0);
	   foundIt = TRUE;
	 }
     }
   if( foundIt )
     print_used_col(B, fpt, t);
   fprintf(fpt, "\n");
}


void print_used_col(Model B, FILE *fpt, int t)
{
   int i, j, first;

   fprintf(fpt, "                 ");    /* BT 9/12/97 */
   for(j = 0; j <= NEIGHBOR_RESIDUES; j++) {fprintf(fpt, " ");}

   if(B->IP->is_defined[cl_F]) 
      for(j = 0; j < B->IP->nMotifLen[t]; j++)
         fprintf(fpt, "*");
   else {
      first = FirstCol(B->F, t);
      for(i = 0; i < B->F->FragWidth[t]; i++)
         if(B->F->nColMask[t][i + first] == COL_ON)  fprintf(fpt, "*");
         else   /* COL NOT USED */                   fprintf(fpt, " ");
   }
   fprintf(fpt, "\n");
}

void print_col_desc(FILE *fpt)
{
   fprintf(fpt, "\n\n");
   fprintf(fpt, "Column 1 :  Sequence Number, Site Number\n");
   fprintf(fpt, "Column 2 :  Left End Location\n");
   fprintf(fpt, "Column 4 :  Motif Element\n");
   fprintf(fpt, "Column 6 :  Right End Location\n");
   fprintf(fpt, "Column 7 :  Probability of Element\n");
   fprintf(fpt, "Column 8 :  Forward Motif (F) or Reverse Complement (R) \n");
   fprintf(fpt, "Column 9 :  Sequence Description from Fast A input\n\n");
}

void print_it(Model B, int t, int location, FILE *fpt, double probability, 
              short RevComp, int *last_seq, int *site_num, FILE* snFpt,
	      int nMaxWidth)

   /*=====================================================================*/
   /* FUNCTION NAME : print_it                                            */
   /* DESCRIPTION   :                                                     */
   /*=====================================================================*/

{
   int  j, flag, seq_start, seq_num=0;
   int  beg, fin, beg_sp, fin_sp;
   int  length;
   char **desc;
   int  nLen;		/* BT 3/24/97 */
   int  first;
   int  totalLength;
   int  k;

   if(!B->IP->is_defined[cl_F])  length = B->F->FragWidth[t];
   else                          length = B->IP->nMotifLen[t]; 

   if(!B->IP->is_defined[cl_F])       
     {
       first = FirstCol(B->F, t);	     
       totalLength =  B->F->FragWidth[t] + (B->F->FragWidth[t] - B->IP->nMotifLen[t]);
     }
   else
     totalLength = length;

   nLen = length;		/* BT 3/24/97 */

   NEWP( desc, B->IP->nNumSequences, char );     /* BT 5/21/97 */
   for( j = 0; j < B->IP->nNumSequences; j++ )
     {
       NEW( desc[j], 64, char );
       strncpy( desc[j], B->IP->fastaHeader[j], 50 );
     }

   j = 0;
   flag = FALSE;
   while((j < B->IP->nNumSequences) && !flag)
     {
       if(location < (*B->Seq->nvEndLocs)[j]) 
	 {
	   flag = TRUE;
	   seq_num = j;
	   length = min( length, (*B->Seq->nvEndLocs)[j] - location );	/* BT 3/24/97 */
	 }
       j++;
     }
   if(seq_num == *last_seq)  (*site_num)++;
   else                      (*site_num) = 1;
   (*last_seq) = seq_num;
   beg_sp = fin_sp = 0;
   if(seq_num != 0)
     seq_start = location - ((*B->Seq->nvEndLocs)[seq_num - 1]);
   else
     seq_start = location; 

   beg = location - NEIGHBOR_RESIDUES;
   if(seq_num != 0) 
     {
      if(beg <= (*B->Seq->nvEndLocs)[seq_num-1]) 
	{
	  beg_sp = (*B->Seq->nvEndLocs)[seq_num - 1] - beg;
	  beg = (*B->Seq->nvEndLocs)[seq_num -1];
	}
     }
   else 
     {
       if(beg < 0) 
	 {
	   beg_sp = 0 - beg;
         beg = 0;
	 }
     }

   fin = location + nLen + NEIGHBOR_RESIDUES - 1;   /* BT 3/24/97 */
   if(fin >= (*B->Seq->nvEndLocs)[seq_num]) 
     {               /* BT 6/2/97 */
       fin_sp = fin - (*B->Seq->nvEndLocs)[seq_num] + 1;
       fin = (*B->Seq->nvEndLocs)[seq_num] - 1;
     }
   
   if( nMaxWidth == 0 )  /* BT 5/4/98 */
     fprintf(fpt, "%4d,%3d ", seq_num+1,*site_num);
   else
     fprintf(fpt, "%4d,%3d,%3d ", seq_num+1,*site_num, t);

   if(RevComp) 
     {
      fprintf(fpt, " %6d", seq_start + length);       /* BT 9/12/97 */
      fprintf(fpt, " ");
      for(j = 0; j < fin_sp; j++)
	fprintf(fpt, " ");
      for(j = fin; j >= location + length; j--)
	fprintf(fpt, "%c", curr_ch(complement((*B->Seq->R)[j]), 
				   B->IP->nAlphaLen, FALSE));

      fprintf(fpt, " ");
      
      if( length < nMaxWidth )           /* BT 5/4/98 */
	for( j = 0; j < nMaxWidth - length; j++ )
	  fprintf(fpt, " ");
      
      for(k = 0, j = length - 1; j >= 0; j--, k++)
	{      
	  fprintf(fpt, "%c", curr_ch(complement((*B->Seq->R)[location + j]),
				     B->IP->nAlphaLen, TRUE));
	  /*	  if( ! B->IP->is_defined[cl_F] &&  
	      B->F->nColMask[t][k + first] != COL_ON)
	    fprintf(fpt, " " ); */
	  if( snFpt != NULL )                                    /* BT 9/17/97 */  
	    fprintf(snFpt, "%c", curr_ch(complement((*B->Seq->R)[location + j]),      
					 B->IP->nAlphaLen, TRUE));
	}

      fprintf(fpt, " ");
      for(j = location - 1; j >= beg; j--)
	fprintf(fpt, "%c", curr_ch(complement((*B->Seq->R)[j]),
				   B->IP->nAlphaLen, FALSE));
      for(j = 0; j <= beg_sp; j++)
	fprintf(fpt, " ");
      fprintf(fpt, " %6d", seq_start + 1);            /* BT 10/1/97 */
     }
   else 
     {
       fprintf(fpt, " %6d ", seq_start + 1);           /* BT 9/12/97 */
       for(j = 0; j < beg_sp; j++)                    
         fprintf(fpt, " ");                         
       for(j = beg; j < location; j++)   /* preceding residues */
         fprintf(fpt, "%c", curr_ch((*B->Seq->R)[j], B->IP->nAlphaLen, FALSE));
       fprintf(fpt, " ");
       
       for(j = 0; j < length; j++)                 /* motif residues */
	 {                                         /* BT 9/17/97 */	   
	   /*	   if( ! B->IP->is_defined[cl_F] &&  
	       B->F->nColMask[t][length - j - 1 + first] != COL_ON)
	     fprintf(fpt, " " ); */

	   fprintf(fpt, "%c", curr_ch((*B->Seq->R)[location + j],
				      B->IP->nAlphaLen, TRUE));

	   if( snFpt != NULL )                      /* BT 9/17/97 */  
	     fprintf(snFpt, "%c", curr_ch((*B->Seq->R)[location + j], 
					  B->IP->nAlphaLen, TRUE));
	 }
       
       if( length < nMaxWidth )           /* BT 5/4/98 */
	 for( j = 0; j < nMaxWidth - length; j++ )
	   fprintf(fpt, " ");
       
       fprintf(fpt, " ");
       
       for(j = location + length; j <= fin; j++) /*tail*/
         fprintf(fpt, "%c", curr_ch((*B->Seq->R)[j], B->IP->nAlphaLen, FALSE));
       for(j = 0; j <= fin_sp; j++)
         fprintf(fpt, " ");
       fprintf(fpt, " %6d", seq_start +length);    /* 10/1/97 */
     }
   fprintf(fpt, " %6.2f ", probability);
   if(RevComp)  fprintf(fpt, "R ");
   else         fprintf(fpt, "F ");
   fprintf(fpt, "%s\n", desc[seq_num]);
   
   if( snFpt != NULL )
     fprintf(snFpt, "\n" );
   
   FREEP( desc, B->IP->nNumSequences );    /* BT 5/21/97 */
}

void print_curr(MCstruct *curr)
{
   int i, len, diff;
   
   len = max(strlen(curr->desc), strlen(curr->cmp_desc));
   printf("DESCRIPTION");
   for(i = 0; i < len - 11; i++)
      printf(" ");
   printf("     Logl    lrt   pval  DOF\n___________");
   diff = strlen(curr->desc);
   for(i = 0; i < len - 11; i++)
      printf("_");
   printf("_____________________________\n%s", curr->cmp_desc);
   if(len > 11) {
      diff = len - strlen(curr->cmp_desc);
      for(i = 0; i <= diff; i++)
         printf(" ");
   }
   else {
      diff = 11 - strlen(curr->cmp_desc);
      if(diff > 0) {
         for(i = 0; i <= diff; i++)
            printf(" ");
      }
   }
   printf("%10.1f\n%s ", curr->cmp_logl, curr->desc);
   if(len > 11) {
      diff = len - strlen(curr->desc);
      for(i = 0; i < diff; i++)
         printf(" ");
   }
   else {
      diff = 11 - strlen(curr->desc);
      if(diff > 0) {
         for(i = 0; i < diff; i++)
            printf(" ");
      }
   }

   printf("%10.1f %6.1f %4.2f %3.0f\n\n", curr->logl, curr->lrt, curr->pval, 
           curr->dof);



}

void print_all_desc(MCstruct *First)
{
   int i, len, diff;
   int currsmall, currlarge;
   int small=1000, large=0;
   MCstruct *curr;
  
   curr = First;

   while (curr != NULL) {
      currsmall = min(strlen(curr->desc), strlen(curr->cmp_desc));
      currlarge = max(strlen(curr->desc), strlen(curr->cmp_desc));
      if(currsmall < small)
         small = currsmall;
      if(currlarge > large)
         large = currlarge;
      curr = curr->next;
   }
   len = large;
   printf("DESCRIPTION");
   for(i = 0; i < len - 11; i++)
      printf(" ");
   printf("     Logl    lrt   pval  DOF\n___________");
   for(i = 0; i < len - 11; i++)
      printf("_");
   printf("_____________________________\n");
   curr = First;
   while (curr != NULL) { 
      diff = strlen(curr->desc);
      printf("%s", curr->cmp_desc);
  
      if(len > 11) {
         diff = len - strlen(curr->cmp_desc);
         for(i = 0; i <= diff; i++)
            printf(" ");
      }
      else {
         diff = 11 - strlen(curr->cmp_desc);
         if(diff > 0) {
            for(i = 0; i <= diff; i++)
               printf(" ");
         }
      }
      printf("%10.1f\n%s ", curr->cmp_logl, curr->desc);
      if(len > 11) {
         diff = len - strlen(curr->desc);
         for(i = 0; i < diff; i++)
            printf(" ");
      }
      else {
         diff = 11 - strlen(curr->desc);
         if(diff > 0) {
            for(i = 0; i < diff; i++)
               printf(" ");
         }
      }

      printf("%10.1f %6.1f %4.2f %3.0f\n\n", curr->logl, curr->lrt, curr->pval, 
              curr->dof);
      curr = curr->next;
   }
} 
     
   
char curr_ch(char c, int alphalen, int caps)

   /*==============================================================*/
   /* FUNCTION NAME : curr_ch                                      */
   /* DESCRIPTION   : Takes the current character in the sequence  */
   /*                 and converts it to the appropriate character */
   /*                 capitalizing it if it is a part of a motif   */
   /*==============================================================*/

{
  switch(c) 
    {
    case 'b' :
      if((alphalen == 4) || (alphalen == 2))
	c = 't';
      else 
	c = 'v';
      break;
    case 'd' :
      if((alphalen == 4) || (alphalen == 2))
	c = 'g'; break;
    case 'j':
      c = 'w'; break;
    case 'o':
      c = 'y'; break;
    case 'u':
      c = 'x'; break; 
    default:
      break;
    }
  c = (char)((int)c - caps * 32);
  return(c);       
}


void PrintSeqData( Model B, MaxResults maxData  )
{
  int     seq;
  int     first;
  int     last;
  IPtype  IP;
  int     i;
  int     j;
  int     t;
  int     last_seq = -1;  /* BT 10/23/2001 */
  int     site_num = 0;   /* BT 10/23/2001 */
  FILE    *fpt;
  FILE    *sn_fpt = NULL;    /* BT 9/17/97 */
  int     nMaxWidth;
  short   bDone;
  int     nMotifCnt = 0;

  IP = B->IP;
  
  if(IP->is_defined[cl_o])  
    fpt = IP->Datafiles->out_fpt;
  else
    fpt = stdout;
  
  nMaxWidth = MaxPrintLenMotif( B );

  for(seq = 0; seq < B->IP->nNumSequences; seq++) 
    {
      if(seq == 0)
	first = 0;
      else
	first = (*B->Seq->nvEndLocs)[seq - 1];
      last = (*B->Seq->nvEndLocs)[seq];

      for( i = first; i < last; i++ )
	{
	  bDone = FALSE;
	  for (t = 0; t < B->IP->nNumMotifTypes && !bDone; t++) 
	    {
	      for(j = 0; j < maxData.nNumMotifs[t]; j++) 
		{
		  if( maxData.nMotifLoc[j][t] == i )
		    {
		      print_it(B, t, maxData.nMotifLoc[j][t], fpt, 
			       maxData.dvMotifProb[j][t], 
			       maxData.RevComp[j][t], &last_seq, &site_num, sn_fpt,
			       IP->nMotifLen[nMaxWidth]);
		      bDone = TRUE;
		      nMotifCnt++;
		      break;
		    }
		  else if( maxData.nMotifLoc[j][t] > i )
		    break;
		}
	    }
	}
    }
  fprintf( fpt, "%d motifs\n", nMotifCnt );
}


void PrintSeqFreqData( Model B, MaxResults maxData  )
{
  int     seq;
  int     first;
  int     last;
  IPtype  IP;
  int     i;
  int     t;
  int     last_seq = 1;
  int     site_num = 1;
  FILE    *fpt;
  FILE    *sn_fpt = NULL;    /* BT 9/17/97 */
  int     nMaxWidth;
  short   bDone;
  int     nMotifCnt = 0;
  short   flag;

  IP = B->IP;
  
  if(IP->is_defined[cl_o])  
    fpt = IP->Datafiles->out_fpt;
  else
    fpt = stdout;
  
  nMaxWidth = MaxPrintLenMotif( B );

  for(seq = 0; seq < B->IP->nNumSequences; seq++) 
    {
      if(seq == 0)
	first = 0;
      else
	first = (*B->Seq->nvEndLocs)[seq - 1];
      last = (*B->Seq->nvEndLocs)[seq];

      for( i = first; i < last; i++ )
	{
	  bDone = FALSE;
	  for (t = 0; t < B->IP->nNumMotifTypes && !bDone; t++) 
	    {
	      if(maxData.frequency[i][t].frequency > B->IP->dNearOptCutoff )  
		{
		  if(maxData.frequency[i][t].fwd_freq > 
		     maxData.frequency[i][t].rev_freq)
		    flag = FALSE;
		  else
		    flag = TRUE;
		  
		  print_it(B, t, i, fpt, 
			    maxData.frequency[i][t].frequency, 
			   flag, &last_seq, &site_num, sn_fpt,
			   IP->nMotifLen[nMaxWidth]);
		  bDone = TRUE;
		  nMotifCnt++;
		  break;
		}
	    }
	}
    }
  fprintf( fpt, "%d motifs\n", nMotifCnt );
}


int MaxPrintLenMotif( Model B )
{
  int t;
  int index = -1;
  int maxLen = 0;
  int length;

  for( t = 0; t < B->IP->nNumMotifTypes; t++ )
    {
      if(!B->IP->is_defined[cl_F])  
	length = B->F->FragWidth[t];
      else                          
	length = B->IP->nMotifLen[t]; 
      if( length > maxLen )
	{
	  index = t;
	  maxLen = B->IP->nMotifLen[t];
	}
    }

  return index;
}


