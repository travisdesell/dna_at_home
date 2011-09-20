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
/* $Id: probability.c,v 1.17 2009/04/23 18:43:54 Bill Exp $       */
/*                                                                        */
/* Author:        (Eric C. Rouchka  July 1, 1996)                         */
/*                (Jun Zhu, October 20, 1996)                             */
/*		  (Bill Thompson 3/3/97)                                  */
/* Description :  This file contains functions that will calculate various*/ 
/*                probabilities, including the probability that the       */
/*                current position is a motif element and the probability */
/*                that the current alignment is the correct one           */ 
/**************************************************************************/

/*------------------------------------------------------------------------*/
/*                                INCLUDE FILES                           */
/*------------------------------------------------------------------------*/
#include "probability.h"

/*------------------------------------------------------------------------*/
/*                           FUNCTION DECLARATIONS                        */
/*------------------------------------------------------------------------*/

const char *dnaAlpha = "abcdnxX";
const char *compDna =  "badcnee";
char *GL_P__;

#define COMP(ch) ((GL_P__ = strchr( dnaAlpha, (char)(ch))) ? \
			 compDna[GL_P__-dnaAlpha] :    \
			 (printf("CH = %c", (char) (ch)),      \
			  p_internal_error("INVALID NUCLEOTIDE ENCOUNTERED"), 0))

#define RCH2INTCOMP(n,R) ((int)(COMP((R)[(n)])) - 97)	/* BT 2/21/97 */


/*------------------------------------------------------------------------*/

double ln_gamma(double xx)
{
  double value;

  value = lgamma( xx );
#ifdef _DEBUG_
  if( ! finite( value ) )
    p_error( "**** Non finite value ****" );
#endif 
  return value;
} 


/* BT 3/5/97 */
/* Adopted from Numerical Recipes */

/* double ln_gamma(double xx)
{
   double x,y,tmp,ser;
   static double cof[6]={76.18009172947146,-86.50532032941677, 24.01409824083091,
                         -1.231739572450155, 0.1208650973866179e-2,-0.5395239384953e-5};
   int j;

   y=x=xx; 
   tmp=x+5.5; 
   tmp -= (x+0.5)*log(tmp); 
   ser=1.000000000190015; 
   for (j=0;j<=5;j++) 
      ser += cof[j]/++y; 
      
   return -tmp+log(2.5066282746310005*ser/x);
   }  */


double xln_gamma(double a)
{
  /* NOTE: This ln_gamma is deprecated. It overflows for small values (< 1.e-16) */
  /* gammln is called instead */

   /*=============================================================*/
   /* FUNCTION NAME : ln_gamma                                    */
   /*                                                             */
   /* DESCRIPTION : This function returns the natural log of the  */
   /*               gamma function.  It is adapted from a version */
   /*               found in numerical recipies                   */
   /*=============================================================*/
  
   /*------- LOCAL VARIABLES -------*/ 
   double x, tmp;
   int j;
   double cof[6] = {76.18009173, -86.50532033, 24.01409822, 
                    -1.231739516, 0.00120858003, -.00000536382};
   double stp = 2.50662827465;
   double ser;
   
   x = a - 1.0;
   tmp = x + 5.5;
   tmp = (x + 0.5) * log(tmp) - tmp;
   ser = 1.0;
   for(j = 0; j < 6; j++) {
      x = x + 1.0;
      ser = ser + cof[j] / x;
   }

   return (tmp + log(stp * ser));
}


/******************   FindMapProb  *********************************/
/*                                                                 */
/* INPUT PARAMETERS :                                              */
/*          C  : Background and motif residue counts               */
/*          IP : Input Parameters from the command line            */ 
/*          t  : Motif Type finding MAP for                        */
/*                                                                 */
/* DESCRIPTION :                                                   */
/*     This function returns the natural log of the probability    */
/*     that the current motif alignment is the one we desire       */
/*=================================================================*/

double FindMapProb(Model B, int t)		/* BT 3/3/97 */
{ 
   int     i, j, nonsites, nNumMotifs;
   double  result, alpha, total=0.0;
   double  dPseudo; 
   double  bgValue = 0.0; 
   double  dBetaPart = 0.0;
   double  value=0.0;
   double  **fCounts;       /* for quick access */
   double  **dPseudoCounts; /* for quick access */
   double  bgPseudo = 0.0;
   Ctype   C;
   IPtype  IP;

   C = B->C;
   IP = B->IP;

   fCounts=C->fCounts[t];
   dPseudoCounts=C->dPseudoCounts[t];

   nNumMotifs = NUMMOTIFS(IP->nNumMotifs[t]); 
  
   /*************************/
   /* Motif element portion */
   /*************************/ 
   for(j = 1; j <= IP->nMotifLen[t]; j++) 
   {  
      dPseudo = 0.0;
      for(i = 0; i < IP->nAlphaLen; i++) 
      {     
          value += ln_gamma((double)fCounts[j][i]+  
                           dPseudoCounts[j][i]);     
          value -= ln_gamma( dPseudoCounts[j][i] );
          dPseudo += dPseudoCounts[j][i];
      } 
      value += ln_gamma( dPseudo );
      value -= ln_gamma( (double) nNumMotifs + dPseudo );
   } 

   /*-------------------*/
   /* background portion*/
   /*-------------------*/
   for(i=0; i < IP->nAlphaLen; i++) 
   {                    
      bgValue += ln_gamma((double)fCounts[BG][i] +  
                       dPseudoCounts[BG][i]);
      bgValue -= ln_gamma( dPseudoCounts[BG][i] );      
      total += fCounts[BG][i];
      bgPseudo += dPseudoCounts[BG][i];   /* BT 8/7/98 */
   }
   /*   bgValue -= ln_gamma(total + C->dSumBGPseudo);
	bValue += ln_gamma( C->dSumBGPseudo ) ;  */
   bgValue -= ln_gamma(total + bgPseudo);    /* BT 8/7/98 */
   bgValue += ln_gamma( bgPseudo ) ;

   result = value + bgValue; 

   nonsites = IP->nPossSites[t] - nNumMotifs;
   alpha = C->dtot_pseudo[t] - C->dbg_pseudo[t];

   if( IP->is_defined[cl_E] )  
     {
       for( i = 0; i < IP->nNumSequences; i++ )	 
	 dBetaPart -= log( B->AP->dAlignCount[B->RP->nSites[i]][i] );
     }
   else
     {  
       /*-----------*/
       /* Beta dist.*/
       /*-----------*/
       dBetaPart = LnBeta((double)nNumMotifs + alpha,             
			  (double)nonsites + C->dbg_pseudo[t]) -  
	           LnBeta((double)alpha, C->dbg_pseudo[t]);       
     }

   result += dBetaPart;
   
   return result;
} 




long double LnBeta(double a, double b)
{
   /*=================================================================*/
   /* FUNCTION NAME : LnBeta                                          */
   /*                                                                 */
   /* INPUT PARAMETERS : a, b: values to find beta                    */
   /*                                                                 */
   /* DESCRIPTION : This function returns the natural log of the beta */
   /*               function of two values (it is a ratio of gammas)  */
   /*=================================================================*/
   if((a == 0.0) || (b == 0.0))
      return 0.0;
   else 
      return (ln_gamma(a) + ln_gamma(b) - ln_gamma((a+b)));
}

/****************  set_posterior_prob  ***********************************/
/*                                                                       */
/*                                                                       */
/* OUTPUT PARAMETERS :                                                   */
/*    IP : Input parameters from the command line tweaked and calculated */
/*    C  : Counts for each of the motif and background types             */
/*                                                                       */
/* DESCRIPTION :                                                         */
/*    This function sets the posterior probabilities based on a weight   */
/*    and the number of observed sites                                   */
/*=======================================================================*/

void set_posterior_prob(IPtype IP, Ctype C)
{

   double weight;
   int t, i;
   int nNumMotifs;
   double dmodel_sites, dmodel_pseudo;

   /* the description of weight is in Protein Science (1995) 4:p.1631 */
   weight = IP->dPseudoSiteWt / (1.0 - IP->dPseudoSiteWt);

   for(t = 0; t < IP->nNumMotifTypes; t++) {
       nNumMotifs = NUMMOTIFS(IP->nNumMotifs[t]);
       C->dtot_sites[t] = IP->nPossSites[t];
       C->dtot_pseudo[t] = C->dtot_sites[t] * weight; 
       C->dTot[t] = C->dtot_sites[t];
       C->dbg_pseudo[t] = C->dtot_pseudo[t];
       dmodel_sites = 0.0; dmodel_pseudo = 0.0;
       for(i = 0; i < IP->RevComplement + 1; i++) {
          C->dmodel_sites[t][i] = (double)IP->nNumMotifs[t][i];
          C->dmodel_pseudo[t][i] = C->dmodel_sites[t][i] * weight;
          C->dTot[t] += C->dmodel_sites[t][i];
          C->dbg_pseudo[t] -= C->dmodel_pseudo[t][i];
          dmodel_sites += C->dmodel_sites[t][i];
          dmodel_pseudo += C->dmodel_pseudo[t][i];
       }   
       /* the description of this posterior in  */
       /* Protein Science (1995) 4: p.1631.     */
       IP->dposterior_prob[t] = (dmodel_sites + dmodel_pseudo) /
	 (C->dtot_sites[t] + C->dtot_pseudo[t]);  
       /* IP->dposterior_prob[t] = (dmodel_sites + dmodel_pseudo) /
	  (C->dtot_sites[t] + C->dtot_pseudo[t] - 1); */  /* BT 8/11/97 */
   }
}

void update_posterior_prob(Model B)

   /*=================================================================*/
   /* FUNCTION NAME : update_posterior_prob                           */
   /*                                                                 */
   /* DESCRIPTION : updates the posterior probability once it has been*/
   /*               set                                               */
   /*=================================================================*/

{
   int t, nNumMotifs;
   double dmodel_pseudo;

    for(t = 0; t < B->IP->nNumMotifTypes; t++) {
       nNumMotifs = NUMMOTIFS(B->IP->nNumMotifs[t]);
       B->C->dmodel_sites[t][FORWARD] = (double)B->IP->nNumMotifs[t][FORWARD];
       dmodel_pseudo = B->C->dmodel_pseudo[t][FORWARD];
       if(B->IP->RevComplement) {
         B->C->dmodel_sites[t][REVERSE] = (double)B->IP->nNumMotifs[t][REVERSE];
         dmodel_pseudo += B->C->dmodel_pseudo[t][REVERSE];
       }
       B->IP->dposterior_prob[t] = ((double)nNumMotifs + dmodel_pseudo) / 
                           (B->C->dtot_sites[t] + B->C->dtot_pseudo[t]);
       /*       B->IP->dposterior_prob[t] = ((double)nNumMotifs + dmodel_pseudo) / 
		(B->C->dtot_sites[t] + B->C->dtot_pseudo[t] - 1.0); */
    } 
}


void init_prob(ProbStruct *P, IPtype IP)
{
   /*==================================================*/
   /* FUNCTION NAME : init_prob                        */
   /*                                                  */
   /* DESCRIPTION : allocated the structures that will */
   /*               store the probabilities            */
   /*==================================================*/

   int t, j, maxlen;
  
   maxlen = findMaxMotifLen(IP);
   NEWPP(P->dvInMotifProb,IP->nNumMotifTypes, double);   /* BT 12/21/99 */
   NEWP(P->dvInBGProb, IP->nNumMotifTypes, double);
   for( t = 0; t < IP->nNumMotifTypes; t++ )
     {
       NEW(P->dvInBGProb[t], IP->nAlphaLen, double); /* BT 12/21/99 */
       NEWP(P->dvInMotifProb[t], maxlen, double);
       for(j = 0; j < maxlen; j++) 
	 NEW( P->dvInMotifProb[t][j], IP->nAlphaLen, double );
     }
}

void free_prob(ProbStruct *P, IPtype IP)
{
  int maxlen;
  
  maxlen = findMaxMotifLen(IP);
  FREEPP(P->dvInMotifProb, IP->nNumMotifTypes, maxlen );
  FREEP(P->dvInBGProb, IP->nNumMotifTypes);  
}


void update_prob(ProbStruct *P, Model B, short bUpdateBkgrnd)
{
   /*======================================================*/
   /* FUNCTION NAME : update_prob                          */
   /*                                                      */
   /* DESCRIPTION : This function updates the probability  */
   /*               that each residue is in the background */
   /*               or in the motif                        */
   /*======================================================*/

   int            n, j, t, nNumMotifs;
   double         BGDenom, MotifDenom;
   double         dPseudo;     /* BT 5/14/97 */
   int            k;
   int            totalMotifs;
   int            mid;
   double         **fq;
   double         *denom;
   register Ctype C;
   IPtype         IP;
   double         *conc;
   int            c;
   
   C  = B->C;
   IP = B->IP;
   
   for( totalMotifs = 0, t = 0; t < IP->nNumMotifTypes; t++) 
     {
       totalMotifs += NUMMOTIFS(IP->nNumMotifs[t]) * IP->nMotifLen[t];
     }
   
   for(t = 0; t < IP->nNumMotifTypes; t++) 
     {
       nNumMotifs= NUMMOTIFS(IP->nNumMotifs[t]);
       for( dPseudo = 0.0,k = 0; k < IP->nAlphaLen; k++ )  
	 dPseudo += C->dPseudoCounts[t][BG][k];
       
       for(n = 0; n < IP->nAlphaLen; n++) 
	 {
	   if( bUpdateBkgrnd )
	     {
	       if( B->WT == NULL )
		 {
		   if(IP->nAlphaLen == 20)       /* BT 8/7/98 */
		     BGDenom = (double)(C->nTotBack - IP->nNumProcessed - 
					totalMotifs) + dPseudo;
		   else
		     BGDenom = (double)(C->nTotBack) + dPseudo; 
		 }
	       else		 
		 BGDenom = C->dBackCnt + dPseudo; 	       
	       P->dvInBGProb[t][n] = ((double)C->wCounts[t][BG][n] + 
				      C->dPseudoCounts[t][BG][n]) / BGDenom;  
	     }
	 }
     }

   if( (! IP->is_defined[cl_R]) && (! IP->is_defined[cl_c]) && (! IP->is_defined[cl_a]) )
     {
       for(t = 0; t < IP->nNumMotifTypes; t++) 
	 {
	   nNumMotifs= NUMMOTIFS(IP->nNumMotifs[t]);
	   for(j = 0; j < IP->nMotifLen[t]; j++) 
	     {
	       for( MotifDenom = 0, dPseudo = 0.0,n = 0; n < IP->nAlphaLen; n++ ) 
		 {
		   dPseudo += C->dPseudoCounts[t][j+1][n];
		   MotifDenom += C->wCounts[t][j+1][n] + C->dPseudoCounts[t][j+1][n];
		 }
	       
	       /* MotifDenom = (double) nNumMotifs + dPseudo; */
	       
	       for(n = 0; n < IP->nAlphaLen; n++)
		 {
		   if(C->dPseudoCounts[t][j+1][n]>  0.00000000000001) 
		     {
		       P->dvInMotifProb[t][j][n] = (C->wCounts[t][j+1][n] + 
						    C->dPseudoCounts[t][j+1][n])/MotifDenom; 
		     }
		   else 
		     {
		       p_internal_error("Update Prob: 0 pseudocount\n");
		     }
		 }
	     }
	 }
     }
   else
     {
       for( t = 0; t < IP->nNumMotifTypes; t++ )
	 {
	   mid = IP->nMotifLen[t] / 2 + IP->nMotifLen[t] % 2 - 1;
	   
	   NEWP( fq, IP->nMotifLen[t], double );
	   for( j = 0; j < IP->nMotifLen[t]; j++ )
	     NEW( fq[j], IP->nAlphaLen, double );
	   NEW( denom, IP->nMotifLen[t], double );
	   
	   for( j = 0; j < IP->nMotifLen[t]; j++ )
	     {
	       for( n = 0; n < IP->nAlphaLen; n++ )
		 {
		   if( IP->AltModel->Palandromic[t][j] )
		     {
		       if( j <= mid )
			 {
			   fq[j][n] = 
			     (((double)C->wCounts[t][j+1][n] + C->dPseudoCounts[t][j+1][n]) +
			      ((double)C->wCounts[t][IP->nMotifLen[t]-j][nComp[n]] + 
			       C->dPseudoCounts[t][IP->nMotifLen[t]-j][nComp[n]]));
			 }
		       else
			 {
			   fq[j][n] = fq[IP->nMotifLen[t] - j - 1][nComp[n]];
			 }			     
		     }
		   else if( IP->AltModel->Repeat[t][j] )
		     {
		       if( j <= mid )
			 {
			   fq[j][n] = 
			     (((double)C->wCounts[t][j+1][n] + C->dPseudoCounts[t][j+1][n]) +
			      ((double)C->wCounts[t][j + mid + 1 + 1][n] + 
			       C->dPseudoCounts[t][j+mid+1 + 1][n]));   /* BT 04/21/04 */
			 }
		       else
			 {
			   fq[j][n] = fq[j - mid - 1][n];
			 }			     
		     }
		   else if( IP->AltModel->Collapsed[t][j] )  /* BT 12/28/99 */
		     {
		       fq[j][n] = 
			 (((double)C->wCounts[t][j+1][n] + C->dPseudoCounts[t][j+1][n]) +
			  ((double)C->wCounts[t][j+1][nComp[n]] + 
			   C->dPseudoCounts[t][j+1][nComp[n]]));
		     }
		   else
		     {
		       fq[j][n] = 
			 ((double)C->wCounts[t][j+1][n] + C->dPseudoCounts[t][j+1][n]);
		     }
		   denom[j] += fq[j][n];		   
		 }
	     }
	 
	   if( IP->is_defined[cl_a] )  /* BT 12/08/03 */
	     {
	       for(c = 1; c <= IP->AltModel->Concen->NumConcen[t]; c++) 
		 { 
		   NEW( conc, IP->nAlphaLen, double );
		   for( j = 0; j < IP->nMotifLen[t]; j++ )
		     {
		       if( IP->AltModel->Concen->Concentrated[t][j] == c ) 
			 {			   
			   for( n = 0; n < IP->nAlphaLen; n++ )
			     {
			       conc[n] += fq[j][n];			       
			     }
			 }
		     }
		   for( j = 0; j < IP->nMotifLen[t]; j++ )
		     {
		       if( IP->AltModel->Concen->Concentrated[t][j] == c ) 
			 {
			   denom[j] = 0;
			   for( n = 0; n < IP->nAlphaLen; n++ )
			     {
			       fq[j][n] = conc[n];
			       denom[j] += conc[n];
			     }
			 }
		     }
		   free( conc );
		 }
	     }

	   for( j = 0; j < IP->nMotifLen[t]; j++ )
	     {
	       for( n = 0; n < IP->nAlphaLen; n++ )
		 P->dvInMotifProb[t][j][n] = fq[j][n] / denom[j];  /* BT 12/21/99 */
	     }
	       
	   FREEP( fq, IP->nMotifLen[t] );
	   free( denom );
	 }
     }

   P->update = FALSE;
}
      
/* BT 12/17/97 */
double in_motif_prob(Model B, register ProbStruct P, int i, int t, 
		     short rev_comp, short bSamplingRatio, 
		     double *dMotifProb, double *dBGProb)
{
   /*===========================================================*/
   /* FUNCTION NAME : in_motif_prob                             */
   /*                                                           */
   /* DESCRIPTION : This function returns the probability that  */
   /*               the current position is in the alignment    */ 
   /*               Any position can be a palindromic position  */
   /*               or a collapsed alphabet position or both    */
   /*               while also being a reverse complementary    */
   /*               motif                                       */
   /*===========================================================*/

  register int     j, index, n, J, mid, compindex;
  register double  dOdds=1.0, dFinalProb, motifProb=1.0, bgProb=1.0;
  int              curr_fwd, curr_back, last_fwd, last_back, last, first, curr;
  double           dPosBGProb;	/* BT 4/19/97 */
  int              nSeq;
  int              nPos;
  int              nOffset = 0;
  int              index2 = 0;
  int              bgIndex;
  int              bgIndex2 = 0;
  BkgType          BP;
  IPtype           IP;
  Ftype            F;  
  PhyloType        PH;
  int              frag_width = 0;
  double           *dInMotifProb;
  double           **dvInMotifProb;
  double           *dInBGProb;
  char             *processedSTR;
  int              *palandromic;
  int              *collapsed;
  short            *concentrated;
  char             *R;
  int              *fragPos = 0;
  short            *is_defined;
  int              nDefault;
  int              nSeqStart = 0;
  int              nAlphaLen;
  double           fgProb;
  double           bkProb;
  int              nSpecies;
  int              p;
  double           wt;
  int              *nMer = NULL;
  int              *nMerPos = NULL;
  short            *orig = NULL;
  double           sum = 0;
  int              bPhyloSeq;
  int              bNmerBkgnd;

  BP = B->BP;
  IP = B->IP;
  F = B->F;
  R = *(B)->Seq->R;
  PH = B->Phylo;

  nAlphaLen = IP->nAlphaLen;
  
  processedSTR = (*B->Seq->ProcessedSTR); /* 12/22/99 */
  palandromic = IP->AltModel->Palandromic[t];
  concentrated = IP->AltModel->Concen->Concentrated[t];
  collapsed = IP->AltModel->Collapsed[t];
  dvInMotifProb = P.dvInMotifProb[t];
  is_defined = IP->is_defined;

  nDefault = (! is_defined[cl_R] ) && ( ! is_defined[cl_a] ) && 
             (! is_defined[cl_c] ) && ( ! is_defined[cl_D] );
  
  dInBGProb = P.dvInBGProb[t]; /* BT 12/21/99 */
  
  nSeq = SequenceFromPosition( B, i );
  bPhyloSeq = IsPhyloSeq( B, nSeq );
  bNmerBkgnd = IsNmerBkgnd( B, nSeq );
  nSpecies = SpeciesInc( B, nSeq );
  wt = GetSeqWeight( B, nSeq, t );
  if( is_defined[cl_B] || is_defined[cl_D] || PH->treeCount > 0 || bNmerBkgnd )
    {
      if( (is_defined[cl_D] || PH->treeCount > 1) && 
	  nSeq == B->Phylo->phyloTreeStartSeq[B->Phylo->phyloIndex[nSeq]] )
	nOffset = SequenceLength( B, nSeq ); 
      nSeqStart = SequenceStartPos( B, nSeq );
      nPos = i - nSeqStart;     
    }

  if( IP->is_defined[cl_freq_background] )
    {
      orig = B->Seq->Orig[nSeq];
      NEW( nMer, IP->nMotifLen[t], int );
      NEW( nMerPos, IP->nMotifLen[t], int );
    }

  J = min( IP->nMotifLen[t], SequenceEndPos( B, nSeq ) - i + 1 );

  mid = J / 2 + J % 2 - 1;		/* BT 4/23/97 */
  
  if(!is_defined[cl_F]) 
    {
      first = FirstCol(F, t);
      last = last_fwd = -1;
      last_fwd = first - 1;  /* BT 12/22/99 */
      last_back = F->nMaxLen[t];
      frag_width = LastCol(F, t) - first + 1;
      fragPos = F->fragPos[t];
    }
    
  for(j = 0; j < J; j++) 
    {
      if(!is_defined[cl_F]) 
	{
	  curr_fwd = fragPos[j];	  
	  curr_back = frag_width - curr_fwd - 1;  /* BT 08/11/99 */	  
	  curr = curr_fwd;
	}   
      else 
	{
	  curr = curr_fwd = j; 
	  curr_back = J - j - 1;
	}
      
      if(!rev_comp)  
	{
	  n = i+curr_fwd;
	  index = RCH2INT(n,R); /* BT 12/22/99 */
	  if( bPhyloSeq ) 
	    {
	      index2 = RCH2INT(n + nOffset, R); /* BT 12/22/99 */
	      bgIndex2 = index2;
	    }

	  bgIndex = index;
	  
	  if(IP->nAlphaLen == 4)
	    {
	      if( IP->RevComplement )
		compindex = nComp[index];
	      else
		compindex = index;
	    }
	}
      else 
	{
	  n = i + curr_back;
	  compindex = RCH2INT(n, R); /* BT 12/22/99 */
	  index = nComp[compindex];
	  if( bPhyloSeq ) 
	    {
	      index2 = RCH2INT(n + nOffset, R); /* BT 12/22/99 */
	      bgIndex2 = index2;
	      index2 = nComp[index2];
	    }

	   bgIndex = compindex; 
	}
      
      dInMotifProb = dvInMotifProb[j];
      if( bNmerBkgnd )
	{
	  if( ! rev_comp )
	    {
	      nMer[j] = orig[n - nSeqStart];
	      nMerPos[j] = n - nSeqStart;
	    }
	  else
	    {
	      nMer[J-j-1] = orig[n - nSeqStart];
	      nMerPos[J-j-1] = n - nSeqStart;
	    }	    
	}

      if( is_defined[cl_B] )
	dInBGProb = BP->dBkgndProb[nSeq][n - nSeqStart];     

      if( bPhyloSeq && ! bNmerBkgnd )
	{
	  if( PH->phyloTree )
	    {
	      if( PH->phyloSpeciesSample[PH->phyloIndex[nSeq]] )
		CalcFelsBySeq( B, P, n, j, t, rev_comp, &fgProb, &bkProb );
	      else
		CalcFels( B, P, n, j, t, rev_comp, &fgProb, &bkProb );

	      
	      motifProb *= fgProb;
	      bgProb *= bkProb;
	    }		
	  else
	    {
	      for( p = 0; p < nSpecies; p++ )
		{
		  if(!rev_comp)  
		    {
		      index = RCH2INT( n + p * nOffset, R );
		      bgIndex = index;
		    }
		  else
		    {		  
		      compindex = RCH2INT(n + p * nOffset, R); 
		      index = nComp[compindex];
		      bgIndex = compindex; 
		    }
		  motifProb *= dInMotifProb[index];
		  if( is_defined[cl_B] )
		    bgProb *= BP->dBkgndProb[nSeq+p][n - nSeqStart][bgIndex];
		  else
		    bgProb *= dInBGProb[bgIndex];  
		} 
	    }
	}
      else
	{
	  motifProb *= dInMotifProb[index];  /* BT 12/21/99 */
	  bgProb *= dInBGProb[bgIndex];  /* BT 12/21/98 */
	}
     }

  if( bNmerBkgnd )
    {
      if( bPhyloSeq )
	{
	  for( bgProb = 0, p = 0; p < nSpecies; p++ )
	    {
	      orig = B->Seq->Orig[nSeq + p];
	      for( j = 0; j < IP->nMotifLen[t]; j++ )
		nMer[j] = orig[nMerPos[j]];
	      bgProb += GetSeqWeight( B, nSeq + p, t ) * NmerProb( B, nSeq + p, t, nMer );
	      sum += GetSeqWeight( B, nSeq + p, t ); 
	    }
	  bgProb /= sum;
	}
      else
	bgProb = NmerProb( B, nSeq, t, nMer );
    }
      
  dOdds = motifProb / bgProb;
  
  if( bSamplingRatio ) 	/* BT 4/16/97 */
    {
      dPosBGProb = 1.0;
      for( j = 0; j < IP->nNumMotifTypes; j++ )	/* BT 4/14/97 */
	{
	  dPosBGProb -= IP->dposterior_prob[j];
	}
      
      if( dPosBGProb <= 0.0 )
	dFinalProb = dOdds;
      else if( bPhyloSeq )  
	dFinalProb = pow( IP->dposterior_prob[t] / dPosBGProb, SpeciesInc( B, nSeq ) ) * dOdds;  /* BT 12/12/2000 */
      else
	dFinalProb = (IP->dposterior_prob[t] / dPosBGProb) * dOdds; 
    }
  else
    dFinalProb = dOdds * IP->dposterior_prob[t] /
      ((1.0 - IP->dposterior_prob[t]) +  dOdds*IP->dposterior_prob[t]); 
  
  *dMotifProb = motifProb;
  *dBGProb = bgProb;
  
  if( IP->is_defined[cl_freq_background] )
    {
      free( nMer );
      free( nMerPos );
    }
  
  return(dFinalProb);
}


/*********** void CheckMoitfProb ******************************/

void CheckMotifProb(Model B, int n, int motif_type, PoSition **Pos,
                    ProbStruct *P, Mlist M, int iterations, double ***dFinalProb)

{
   int      newpos,t;
   short    RevComp;
   PoSition *Pos_t;  /* Pos_t=Pos[t]  */
   double   dMotifProb;
   double   dBGProb;
   int      nLen;
   int      nSeq;
   int      nOffset;
   int      k;

   if(B->IP->is_defined[cl_F])  /* BT 9/11/98 */
     nLen = B->IP->nMotifLen[motif_type];
   else
     nLen = B->F->FragWidth[motif_type];

   nSeq = SequenceFromPosition( B, n );
   nOffset = SequenceLength( B, nSeq );
   
   /* seq overlap, delete old motif data */
   for(t=0;t<B->IP->nNumMotifTypes;t++)
     {
       Pos_t=Pos[t];
       for(newpos=n;newpos<n+nLen;newpos++)
	 {
	   if( Pos_t[newpos].bEndOfSeq )
	     break;				/* BT 4/18/97 */ 
	   if(Pos_t[newpos].nMotifStartPos)
	     {
	       RevComp = Pos_t[newpos].RevComp;        
	       adjust_counts(B,DELETE,newpos,t,RevComp);
	       B->IP->nNumMotifs[t][RevComp]--;
	       not_in_motif(Pos,newpos,B,t);    
	       P->update = TRUE;
	       delete_motif(B, newpos, M, t);			/* BT 3/31/97 */
	       if( IsPhyloSeq( B, nSeq ) )
		 {	
		   for( k = 1; k < SpeciesInc( B, nSeq ); k++ )
		     {
		       RevComp = Pos_t[newpos + k * nOffset].RevComp;        
		       adjust_counts(B,DELETE,newpos+ k * nOffset,t,RevComp);
		       B->IP->nNumMotifs[t][RevComp]--;
		       not_in_motif(Pos,newpos+ k * nOffset,B,t);    
		       delete_motif(B, newpos + k * nOffset, M, t);
		     }
		 }
	       newpos+=(B->IP->nMotifLen[t]-1);
	     }
	 }
     }
   
   if(P->update) 
     update_prob(P, B, (! B->IP->is_defined[cl_b]));
   
   (*dFinalProb)[motif_type][0]=in_motif_prob(B,*P,n,motif_type,FALSE, 
					      TRUE, &dMotifProb, &dBGProb); 
   if(B->IP->RevComplement)
     {
       (*dFinalProb)[motif_type][1]=in_motif_prob(B,*P,n,motif_type,
						  TRUE, TRUE, &dMotifProb, &dBGProb);
     }
}


double CalcMapProb( Model B, short bCalcPalin )
{
  int        i;
  int        j;
  int        t;
  IPtype     IP;
  double     result;
  double     bgValue; 
  double     dBetaPart, value;
  Ctype      C;
  RPType     RP;
  double     fragValue = 0.0;
  double     siteValue = 0.0;
  int        nStart;
  int        bRev;

  IP = B->IP;
  C = B->C;
  RP = B->RP;

  if( B->WT != NULL )  
    {
      for(t = 0; t < B->IP->nNumMotifTypes; t++)   
	{                                          
	  for( i = 0; i < IP->nNumSequences; i++ )
	    {
	      nStart = SequenceStartPos( B, i );
	      for( j = 0; j < SequenceLength( B, i ); j++ )
		{
		  if(  RP->sitePos[i][j][t].nMotifStart )
		    {
		      bRev = RP->sitePos[i][j][t].nRevComp;
		      adjust_counts(B, DELETE, nStart + j, t, bRev);
		      adjust_counts(B, ADD, nStart + j, t, bRev);
		    }
		}
	    }
	}
    }
  
  result = 0.0;
  value = 0.0;

  /*************************/
  /* Motif element portion */
  /*************************/ 

  for( t=0; t < IP->nNumMotifTypes; t++)
    {    
      value += CalcMotifMap( B, t, bCalcPalin );
    }
  
  /*-----------*/
  /* Beta dist.*/
  /*-----------*/
  
  dBetaPart = CalcBetaMap( B, bCalcPalin ); 
  
  /*-------------------*/
  /* background portion*/
  /*-------------------*/

  bgValue = CalcBkgndMap( B, bCalcPalin );

  if( ! B->IP->is_defined[cl_F] )
    {
      for( t=0; t < IP->nNumMotifTypes; t++)       
	fragValue += CalcMotifFragMap( B, t, bCalcPalin );
    }

  siteValue = CalcSitePerSeqMap( B );

  result = value + bgValue + dBetaPart + fragValue + siteValue - IP->dnull_map;  

#ifdef _DEBUG_
  if( ! finite( result ) )
      p_internal_error("**** non finite result ****");
#endif

  return result;
}


double CalcBetaMap( Model B, short bCalcPalin )
{
  int        i;
  int        j;
  int        t;
  IPtype     IP;
  int        nonsites, nNumMotifs;
  double     alpha;
  double     dBetaPart;
  Ctype      C;
  RPType     RP;
  int        nSites;
  int        nLen;
  double     weight;
  double     dTotPseudo;
  double     dBgPseudo;
  int        nPrevSite = -1;
  int        nPrevType;
  int        nPrevDir = 0;
  double     dPartialSum;
  double     dTransPart;
  double     dScale;
  double     dEndPart;
  int        p;

  IP = B->IP;
  C = B->C;
  RP = B->RP;

  /*-----------*/
  /* Beta dist.*/
  /*-----------*/
  
  dBetaPart = 0.0;
  if( IP->is_defined[cl_E] )  
    {
      if( RP->nUseSpacingProb || IP->is_defined[cl_T] )   /* BT 11/17/99 */
	{
	  for( i = 0; i < IP->nNumSequences; i++ )	 
	    {
	      nPrevType = -1;
	      nLen = SequenceLength( B, i );
	      nSites = 0;
	      for( j = nLen - 1; j >= 0; j-- )
		{
		  for( t = 0; t < B->IP->nNumMotifTypes; t++ ) 
		    {
		      if( RP->sitePos[i][j][t].nMotifStart )
			{
			  dBetaPart += log( RP->sitePos[i][j][t].dSpacingProb );
			  dPartialSum = CalcBetaPartialSum(B, i, t, nPrevSite, nPrevType, j, nSites );
			  if( dPartialSum > 0 )
			    dBetaPart -= log( dPartialSum ); 

			  if( nPrevType != -1 )
			    {
			      if( RP->sitePos[i][j][t].nRevComp )  /* BT 2/15/2001 */
				dBetaPart += log( RP->dPTrans[nPrevType][nPrevDir][t][REVERSE] );
			      else
				dBetaPart += log( RP->dPTrans[nPrevType][nPrevDir][t][FORWARD] ); 
			    }

			  nSites++;
			  nPrevSite = j;
			  nPrevType = t;
			  if( RP->sitePos[i][j][t].nRevComp ) 
			    nPrevDir = REVERSE;
			  else
			    nPrevDir = FORWARD;			    
			}
		    }
		}
	    }      
	}
      else
	{
	  dTransPart = 0.0;
	  dEndPart = 0.0;
	  if( IP->RevComplement )
	    dScale = 4.0;
	  else
	    dScale = 1.0;

	  /* if( RP->bUseTrans && RP->bUpdateTrans )
	     UpdateTransMatrix( B ); */

	  for( i = 0; i < IP->nNumSequences; i += SeqInc( B, i ) )	 
	    {
	      nSites = min( B->RP->nSites[i], B->RP->nMaxBlocks );
	      if( B->AP->dAlignCount[nSites][i] != -DBL_MAX )
		dBetaPart -= log( B->AP->dAlignCount[nSites][i] ); 

	      if( IsPhyloSeq( B, i ) && ! B->Phylo->phyloTree  )
		for( p = 1; p < SeqInc( B, i ); p++ )
		  dBetaPart -= log( B->AP->dAlignCount[nSites][i+p] ); 

#ifdef _DEBUG_
	      if( ! finite( dBetaPart ) )
		{
		  p_error("**** non finite result ****");
		  dBetaPart = dBetaPart;
		}
#endif	      	      
	      
	      /* if( RP->bUseTrans && RP->bUpdateTrans )
		{
		  nPrevType = -1;
		  nLen = SequenceLength( B, i );
		  nSites = 0;
		  for( j = nLen - 1; j >= 0; j-- )
		    {
		      for( t = 0; t < B->IP->nNumMotifTypes; t++ ) 
			{
			  if( RP->sitePos[i][j][t].nMotifStart )
			    {
			      if( nSites > 0 )
				{
				  dTransPart += log( dScale * 
						     RP->dPTrans[nPrevType][FORWARD][t][FORWARD] );
				}
			      else
				{
				  dTransPart += log( RP->dEndSiteProb[t] );
				}
			      nSites++;
			      nPrevSite = j;
			      nPrevType = t;
			      if( RP->sitePos[i][j][t].nRevComp ) 
				nPrevDir = REVERSE;
			      else
				nPrevDir = FORWARD;
			    }
			} 
		    }  
		    } */
	    }

	  /* if( RP->bUseTrans && RP->bUpdateTrans )
	    {
	      UpdateTransMatrix( B );
	      weight = RP->dTransWt / (1.0 - RP->dTransWt);
	      
	      for( t = 0; t < B->IP->nNumMotifTypes; t++ ) 
		{
		  dTransRowSum = 0;
		  dPriorTransRowSum = 0;
		  for( t1 = 0; t1 < B->IP->nNumMotifTypes; t1++ ) 
		    {
		      dTransPart += ln_gamma( RP->dTransCnts[t][t1] +  weight * RP->dPriorTransCnts[t][t1] );
		      dTransPart -= ln_gamma( weight * RP->dPriorTransCnts[t][t1] );
		      dTransRowSum += RP->dTransCnts[t][t1];
		      dPriorTransRowSum += RP->dPriorTransCnts[t][t1];
		    }
		  dEndPart += ln_gamma( RP->dEndCnts[t] + weight * RP->dPriorEndCnts[t] );
		  dEndPart -= ln_gamma( weight * RP->dPriorEndCnts[t] );
		  dTransPart -= ln_gamma( dTransRowSum + weight * dPriorTransRowSum );
		  dTransPart += ln_gamma( weight * dPriorTransRowSum );
		} 
	      
	      dEndPart -= ln_gamma( RP->dPosTotalCnts + weight * RP->dTotalPriorEndCnts );
	      dEndPart += ln_gamma( weight * RP->dTotalPriorEndCnts );
	      } */
	  
	  /* dBetaPart += (dTransPart + dEndPart);   */
	}
    }
  else
    {  
      if( IP->is_defined[cl_T] )    /* BT 11/10/99 */
	{
	  weight = IP->dPseudoSiteWt / (1.0 - IP->dPseudoSiteWt);
	  for( t=0; t < IP->nNumMotifTypes; t++)
	    {
	      nNumMotifs = NUMMOTIFS(IP->nNumMotifs[t]); 	  
	      nonsites = RP->nonCutoffSites - nNumMotifs;
	      dTotPseudo = RP->nonCutoffSites * weight;
	      dBgPseudo = dTotPseudo - 
		(B->First->dmodel_sites[t][FORWARD] + B->First->dmodel_sites[t][REVERSE])  * weight;
	      alpha = dTotPseudo - dBgPseudo;

	      dBetaPart += LnBeta((double)nNumMotifs + alpha,             
				  (double)nonsites + dBgPseudo) -  
		LnBeta((double)alpha, dBgPseudo );       
	    }	  
	}
      else
	{
	  for( t=0; t < IP->nNumMotifTypes; t++)
	    {
	      nNumMotifs = NUMMOTIFS(IP->nNumMotifs[t]); 
	      
	      nonsites = IP->nPossSites[t] - nNumMotifs;
	      alpha = C->dtot_pseudo[t] - C->dbg_pseudo[t];
	      
	      dBetaPart += LnBeta((double)nNumMotifs + alpha,             
				  (double)nonsites + C->dbg_pseudo[t]) -  
		LnBeta((double)alpha, C->dbg_pseudo[t]);      
	    }
	}
    }

  return dBetaPart;
}


double CalcBetaPartialSum(Model B, int nSeq, int t, int nPrevSite, int nPrevType, int curPos, int nSites )
{
  RPType     RP;
  IPtype     IP;
  double     dSum = 0;
  int        n;
  double     *pdPGap;
  int        nWidth;

  RP = B->RP;
  IP = B->IP;

  pdPGap = RP->dPGap[t][nPrevType];
  nWidth = MotifWidth( B, t );

  if( RP->nUseGap )
    {
      if( nSites > 0 )
	{
	  for( n = 0; n < nPrevSite - nWidth + 1; n++ )
	    {
	      dSum += RP->sitePos[nSeq][n][t].dSpacingProb *  
		pdPGap[nPrevSite - n - nWidth];
	    }
	}
    }
  else
    {
      /*      if( curPos > 0 )
	{
	  if( nSites == 0 )
	    {
	      nLen = SequenceLength( B, nSeq );
	      dSum =  RP->sitePos[nSeq][nLen - 1][t].dPartialSum;
	    }
	  else
	    dSum = RP->sitePos[nSeq][curPos - 1][t].dPartialSum;
	    } */

      if( nSites > 0 )
	{
	  for( n = 0; n < nPrevSite - nWidth + 1; n++ )
	    {
	      dSum += RP->sitePos[nSeq][n][t].dSpacingProb;
	    }        /* BT 2/15/2001 */
	}
    }

  return( dSum );
}
 

double CalcBkgndMap( Model B, short bCalcPalin )
{
  double     **fCounts;       /* for quick access */
  double     **dPseudoCounts; /* for quick access */
  double     bgValue; 
  double     bgPseudo; 
  double     total; 
  int        i;
  IPtype     IP;
  Ctype      C;
  RPType     RP;
  int        nLen;
  int        t;
  int        j;
  BkgType    BP;
  int        seq;
  int        start;
  double     *bgCounts;
  double     *dPseudo;
  double     impResult = 0.0;
  int        p;

  IP = B->IP;
  C = B->C;
  RP = B->RP;
  BP = B->BP;

  fCounts=C->wCounts[BG];
  dPseudoCounts=C->dPseudoCounts[BG];

  NEW( dPseudo, IP->nAlphaLen, double );  /* BT 08/06/2001 */
  for( i=0; i < IP->nAlphaLen; i++ ) 
    {
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  dPseudo[i] += C->dPseudoCounts[t][BG][i];
	}
    }
  
  bgValue = 0.0;
  bgPseudo = 0.0;
  total = 0.0;

  if( IP->is_defined[cl_B] )
    {
      NEW( bgCounts, IP->nAlphaLen, double );

      for( seq = 0; seq < IP->nNumSequences; seq += SpeciesInc( B, seq ) )	 
	{
	  nLen = SequenceLength( B, seq );
	  start = SequenceStartPos( B, seq );
	  j = 0;
	  while( j < nLen )
	    {
	      if( PossibleBkgndPos( B, start + j, nLen ) ) 
		{		  
		  for( i = 0; i < IP->nAlphaLen; i++)
		    {
		      bgCounts[i] += BP->dBkgndProb[seq][j][i]  * GetSeqWeight( B, seq, 0 );
		      if( IsPhyloSeq( B, seq ) )
			for( p = 1; p < SpeciesInc( B, seq ); p++ )
			  bgCounts[i] += BP->dBkgndProb[seq+p][j][i] * GetSeqWeight( B, seq+p, 0 );;		  
		    }
		}
	      j++;
	    }
	}
      
      for(i=0; i < IP->nAlphaLen; i++) 
	{                    
	  bgValue += ln_gamma(bgCounts[i] +  
			      dPseudo[i]);
	  bgValue -= ln_gamma( dPseudo[i] );      	      
	  total += bgCounts[i];
	  bgPseudo += dPseudo[i];   /* BT 8/7/98 */
	}
      free( bgCounts );
    }
  else
    {    
      for(i=0; i < IP->nAlphaLen; i++) 
	{                    
	  bgValue += ln_gamma((double)fCounts[BG][i] +  
			      dPseudo[i]);
	  bgValue -= ln_gamma( dPseudo[i] );      
	  total += fCounts[BG][i];
	  bgPseudo += dPseudo[i];   /* BT 8/7/98 */
	}
    }

  bgValue -= ln_gamma(total + bgPseudo);    /* BT 8/7/98 */
  bgValue += ln_gamma( bgPseudo ) ;

#ifdef _DEBUG_
  if( ! finite( bgValue ) )
    p_error( "**** Non finite value ****" );
#endif

  free( dPseudo );

  if( B->Phylo->bCalcPhylo  ) 
    impResult = ImpCalcBkgndMap( B );

#ifdef _DEBUG_
  if( ! finite( bgValue + impResult ) )
    p_error( "**** Non finite value ****" );
#endif


  return( bgValue + impResult );
}


double CalcMotifMap( Model B, int t, short bCalcPalin )
{
  int        i;
  int        j;
  IPtype     IP;
  int        nNumMotifs;
  double     result, total;
  double     dPseudo; 
  double     value;
  double     **fCounts;       /* for quick access */
  double     **dPseudoCounts; /* for quick access */
  Ctype      C;
  double     dCnt;
  double     dTotalCnt;
  double     dPseudoCnt;
  int        mid;
  double     *conc;
  double     *pseudoconc;
  int        *concProcessedCols;
  int        c;
  double     impResult = 0.0;

  IP = B->IP;
  C = B->C;

  result = 0.0;
  value = 0.0;
  fCounts=C->wCounts[t];
  dPseudoCounts=C->dPseudoCounts[t];
  
  nNumMotifs = NUMMOTIFS(IP->nNumMotifs[t]); 

  /* keep track of columns processed for concentrated positions */
  NEW( concProcessedCols, IP->nMotifLen[t] + 1, int ); 

  if( IP->is_defined[cl_a] )
    {
      for(c = 1; c <= IP->AltModel->Concen->NumConcen[t]; c++) 
	{ 
	  NEW( conc, IP->nAlphaLen, double );
	  NEW( pseudoconc, IP->nAlphaLen, double );
	  for( j = 1; j <= IP->nMotifLen[t]; j++ )
	    {
	      if( IP->AltModel->Concen->Concentrated[t][j-1] == c ) 
		{			   
		  concProcessedCols[j] = 1;
		  for( i = 0; i < IP->nAlphaLen; i++ )
		    {
		      conc[i] += fCounts[j][i];
		      pseudoconc[i] += dPseudoCounts[j][i];
		    }
		}
	    }
	  
	  total = 0;
	  dPseudo = 0;
	  for( i=0; i < IP->nAlphaLen; i++ ) 
	    {                    
	      value += ln_gamma(conc[i] +  pseudoconc[i] );
	      value -= ln_gamma( pseudoconc[i] );      
	      total += conc[i];
	      dPseudo += pseudoconc[i];
	    }
      
	  value -= ln_gamma( total + dPseudo ); 
	  value += ln_gamma( dPseudo ) ;
      
	  free( conc );
	  free( pseudoconc );
	}
    }
  
  /*************************/
  /* Motif element portion */
  /*************************/ 
  mid = IP->nMotifLen[t] / 2 + IP->nMotifLen[t] % 2;	  
  
  for(j = 1; j <= IP->nMotifLen[t]; j++) 
    {  
      if( ! concProcessedCols[j] )
	{
	  dPseudo = 0.0;
	  if( IP->AltModel->Palandromic[t][j-1] && bCalcPalin )
	    {
	      if( j <= mid )
		{
		  dTotalCnt = 0;
		  for(i = 0; i < IP->nAlphaLen; i++) 
		    {     
		      dCnt = fCounts[j][i] + fCounts[IP->nMotifLen[t] - j + 1][nComp[i]];
		      
		      dPseudoCnt = dPseudoCounts[j][i] + 
			dPseudoCounts[IP->nMotifLen[t] - j + 1][nComp[i]];
		      value += ln_gamma(dCnt + dPseudoCnt ); 
		      value -= ln_gamma( dPseudoCnt );
		      dPseudo += dPseudoCnt;  
		      
		      dTotalCnt += dCnt;
		    } 
		  value += ln_gamma( dPseudo );
		  value -= ln_gamma( dTotalCnt + dPseudo );
		}
	    }
	  else if( IP->AltModel->Repeat[t][j-1] && bCalcPalin )
	    {
	      if( j <= mid )
		{
		  dTotalCnt = 0;
		  for(i = 0; i < IP->nAlphaLen; i++) 
		    {     
		      dCnt = fCounts[j][i] + fCounts[mid + j][i];
		      
		      dPseudoCnt = 
			dPseudoCounts[j][i] + dPseudoCounts[mid + j][i];
		      value += ln_gamma(dCnt + dPseudoCnt ); 
		      value -= ln_gamma( dPseudoCnt );
		      dPseudo += dPseudoCnt;  
		      
		      dTotalCnt += dCnt;
		    } 
		  value += ln_gamma( dPseudo );
		  value -= ln_gamma( dTotalCnt + dPseudo );
		}
	    }
	  else if( IP->AltModel->Collapsed[t][j-1] && bCalcPalin )  /* BT 12/28/99 */
	    {
	      dTotalCnt = 0;
	      for(i = 0; i < IP->nAlphaLen; i += 2) 
		{     
		  dCnt = fCounts[j][i] + fCounts[j][nComp[i]];
		  
		  dPseudoCnt = 
		    dPseudoCounts[j][i] + 
		    dPseudoCounts[j][nComp[i]];
		  value += ln_gamma(dCnt + dPseudoCnt ); 
		  value -= ln_gamma( dPseudoCnt );
		  value -= ln_gamma( dCnt + 2 );
		  value += ln_gamma( (double) fCounts[j][i] + 1.0 ) + 
		    ln_gamma( (double) fCounts[j][i+1] + 1.0 );
		  dPseudo += dPseudoCnt;  
		  
		  dTotalCnt += dCnt;
		} 
	      value += ln_gamma( dPseudo );
	      value -= ln_gamma( dTotalCnt + dPseudo );
	    }
	  else
	    {
	      dTotalCnt = 0;
	      for(i = 0; i < IP->nAlphaLen; i++) 
		{     
		  value += ln_gamma((double)fCounts[j][i]+  
				    dPseudoCounts[j][i]);     
		  value -= ln_gamma( dPseudoCounts[j][i] );
		  dPseudo += dPseudoCounts[j][i];
		  dTotalCnt += (double)fCounts[j][i];
		} 
	      value += ln_gamma( dPseudo );
	      value -= ln_gamma( dTotalCnt + dPseudo );
	    }
	}
#ifdef _DEBUG_
      if( ! finite( value ) )
	p_error( "**** Non finite value ****" );
#endif 
    }
  result = value;

  free( concProcessedCols );

  if( B->Phylo->phyloTree && B->Phylo->bCalcPhylo  ) 
    impResult = ImpSampleMotifMap( B, t );
  
  return result + impResult;
}


double CalcMotifFrag2( Model B, int t, short bCalcPalin )
{
  int        i;
  int        j;
  IPtype     IP;
  int        nNumMotifs;
  double     result;
  double     bgValue; 
  double     value;
  double     **fCounts;       /* for quick access */
  double     **dPseudoCounts; /* for quick access */
  Ctype      C;
  ProbStruct P;

  IP = B->IP;
  C = B->C;
  fCounts=C->fCounts[t];
  dPseudoCounts=C->dPseudoCounts[t];
    
  nNumMotifs = NUMMOTIFS(IP->nNumMotifs[t]); 

  value = CalcMotifMap( B, t, bCalcPalin );
  
  init_prob( &P, IP);
  update_prob( &P, B, (! IP->is_defined[cl_b]) );

  bgValue = 0.0;
  for(j = 1; j <= IP->nMotifLen[t]; j++) 
    {  
      for(i=0; i < IP->nAlphaLen; i++) 
	{        
	  bgValue += (double) fCounts[j][i] * log( P.dvInBGProb[t][i] );  /* BT 12/21/98 */
	}
    }
  
  result = value - bgValue;
        
  free_prob( &P, IP );

  return result;
}


double CalcMotifFragMap( Model B, int t, short bCalcPalin )
{
  int    first;
  int    last;
  int    span = 0;
  int    width = 0;
  int    mid;
  int    cnt;
  int    i;
  int    prev;
  IPtype IP;
  int    fragWidth;

  IP = B->IP;

  if( IP->is_defined[cl_F] )
    return 0.0;

  if( NUMMOTIFS(IP->nNumMotifs[t]) == 0 )
    return 0.0;

  first = FirstCol( B->F, t );
  if( first == -1 )
    return 0.0;

  last = LastCol( B->F, t );
  if( last == -1 )
    return 0.0;

  /* mid = IP->nMotifLen[t] / 2 + IP->nMotifLen[t] % 2 - 1;       */
  fragWidth = MotifWidth( B, t ); 
  mid = fragWidth / 2 + fragWidth % 2 - 1;  

  i = first;
  prev = i;
  cnt = 0;
  while( i != -1 )
    {
      if( (IP->AltModel->Palandromic[t][cnt] || IP->AltModel->Repeat[t][cnt]) && bCalcPalin )
	{
	  if( cnt <= mid )
	    {
	      span += (i - prev);
	      width++;
	    }
	  else
	    break;
	}
      else
	{
	  span += (i - prev);
	  width++;
	}
      prev = i;
      i = NextCol( B->F, t, i);
      cnt++;
    }

  /*  span = last - first - 1; */
  if(  bCalcPalin && (IP->is_defined[cl_R] || IP->is_defined[cl_I]) )
    {
      /* span = mid; */
      /* width--; */
      span = fragWidth / 2 - 1;
      width = IP->nMotifLen[t] / 2 - 1;  /* 03/21/03 */
    }
  else
    {
      span--;
      width -= 2;
    }
  
  return (-LnBiCoef( span, width ));
}


double CalcSitePerSeqMap( Model B )
{
  RPType RP;
  IPtype IP;
  double *cntSitesPerSeq;
  int    seq;
  int    k;
  double siteMap = 0.0;
  double cntSum = 0.0;
  double pseudoSum = 0.0;

  if( B->IP->is_defined[cl_E] )  /* removed && B->RP->nUseFixedBlockProbs 03/11/03 */
    {
      RP = B->RP;
      IP = B->IP;
      
      NEW( cntSitesPerSeq,  RP->nMaxBlocks + 1, double );
      
      for( seq = 0; seq < IP->nNumSequences; seq++ )
	{
	  if( IP->is_defined[cl_V] )
	    {
	      if( seq >= IP->nVerifySeq )
		cntSitesPerSeq[RP->nSites[seq]]++;
	    }
	  else
	    cntSitesPerSeq[RP->nSites[seq]] += GetSeqWeight( B, seq, 0 );	    
	}
      
      for( k = 0; k <= RP->nMaxBlocks; k++ )
	{
	  siteMap += ln_gamma( cntSitesPerSeq[k] + RP->priorSitesPerSeq[k] );
	  siteMap -= ln_gamma( RP->priorSitesPerSeq[k] );
	  cntSum += cntSitesPerSeq[k];
	  pseudoSum += RP->priorSitesPerSeq[k];
	}
      siteMap += ln_gamma( pseudoSum );
      siteMap -= ln_gamma( cntSum + pseudoSum );
      
#ifdef _DEBUG_
      if( ! finite( siteMap ) )
	{
	  siteMap = siteMap;
	}
#endif

      free( cntSitesPerSeq );
    }

  return( siteMap );
}


double CalcNullMap( Model B )
{
  double     **fCounts;       /* for quick access */
  double     **dPseudoCounts; /* for quick access */
  double     bgValue; 
  double     bgPseudo; 
  double     total; 
  int        i;
  IPtype     IP;
  Ctype      C;
  RPType     RP;
  int        nLen;
  int        t;
  int        j;
  BkgType    BP;
  int        seq;
  double     *bgCounts;
  double     *dPseudo;
  int        p;
  int        start;
  double     impResult = 0.0;
  double     seqMap;

  IP = B->IP;
  C = B->C;
  RP = B->RP;
  BP = B->BP;

  fCounts=C->wCounts[BG];
  dPseudoCounts=C->dPseudoCounts[BG];

  NEW( dPseudo, IP->nAlphaLen, double );  /* BT 08/06/2001 */
  for( i=0; i < IP->nAlphaLen; i++ ) 
    {
      for( t = 0; t < IP->nNumMotifTypes; t++ )
	{
	  dPseudo[i] += C->dPseudoCounts[t][BG][i];
	}
    }
  
  bgValue = 0.0;
  bgPseudo = 0.0;
  total = 0.0;

  if( IP->is_defined[cl_B] )
    {
      NEW( bgCounts, IP->nAlphaLen, double );

      for( seq = 0; seq < IP->nNumSequences; seq += SpeciesInc( B, seq ) )	 
	{
	  nLen = SequenceLength( B, seq );
	  start = SequenceStartPos( B, seq );
	  j = 0;
	  while( j < nLen )
	    {
	      if( PossibleBkgndPos( B, start + j, nLen ) ) 
		{
		  for( i = 0; i < IP->nAlphaLen; i++)
		    {
		      bgCounts[i] += BP->dBkgndProb[seq][j][i] * GetSeqWeight( B, seq, 0 );
		      if( IsPhyloSeq( B, seq ) )
			for( p = 1; p < SpeciesInc( B, seq ); p++ )
			  bgCounts[i] += BP->dBkgndProb[seq+p][j][i] * GetSeqWeight( B, seq+p, 0 );
		    }
		}
	      j++;
	    }
	}
      
      for(i=0; i < IP->nAlphaLen; i++) 
	{                    
	  bgValue += ln_gamma(bgCounts[i] +  
			      dPseudo[i]);
	  bgValue -= ln_gamma( dPseudo[i] );      
	  total += bgCounts[i];
	  bgPseudo += dPseudo[i];   /* BT 8/7/98 */
	}
      free( bgCounts );
    }
  else
    {    
      for(i=0; i < IP->nAlphaLen; i++) 
	{                    
	  bgValue += ln_gamma((double)fCounts[BG][i] +  
			      dPseudo[i]);
	  bgValue -= ln_gamma( dPseudo[i] );      
	  total += fCounts[BG][i];
	  bgPseudo += dPseudo[i];   /* BT 8/7/98 */
	}
    }

  bgValue -= ln_gamma(total + bgPseudo);    /* BT 8/7/98 */
  bgValue += ln_gamma( bgPseudo ) ;

  seqMap = CalcSitePerSeqMap( B );
  bgValue += seqMap;
  B->RP->sitesPerSeqNullMap = seqMap;

  free( dPseudo );
  
  if( B->Phylo->bCalcPhylo  ) 
    {
      impResult = ImpSampleNullMap( B );
      B->Phylo->phyloNull = impResult;
    }
  
  return( bgValue + impResult );
}


int PossibleBkgndPos( Model B, int pos, int offset )
{
  if( B->RP->nInMotif[pos] )
    return FALSE;

  if( B->IP->is_defined[cl_V] )   /* TEST */  /* BT 06/15/04 */
    return TRUE;

  /*  if( B->Phylo->phyloTree) 
    {    
      for( n = 0; n < B->Phylo->phyloSpecies; n++ )
	{
	  ch = (*B->Seq->R)[pos + n * offset];
	  if( ch == 'n' || ch == 'x' )
	    return FALSE;
	}
	} */
  
  return TRUE;
}


/*******************************************/
