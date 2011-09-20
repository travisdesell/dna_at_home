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
#include "footprint.h"


void UpdateFootprintParameters( Model B, PoSition **Pos, int nSeq )
{
   int      t;
   Ctype    C;
   IPtype   IP;
   int      n;
   PoSition *Pos_t;
   RPType   RP;
   double   motifCnt;
   double   modelSitesCons;
   double   modelPseudoCons;
   double   modelSitesNonCons;
   double   modelPseudoNonCons;
   double   bgSitesCons;
   double   bgPseudoCons;
   double   bgSitesNonCons;
   double   bgPseudoNonCons;
   double   lenConserved;
   int      nPos;
   double   dModelPseudo;
   double   dBGPseudo;
   double   beta0;
   double   beta1;
   double   frac;
   int      seq;
   double   weight;
   int      nLen;
   int      nStart;
#ifdef _DEBUG_
   char     temp[128];
#endif
   
   C = B->C;
   IP = B->IP;
   RP = B->RP;   

   /*   weight =  IP->dPseudoSiteWt / (1 -  IP->dPseudoSiteWt); */
   weight = RP->dSiteProbWt;
   dModelPseudo = RP->dPriorSites * weight;

   lenConserved = 0;
   modelSitesCons = 0;       

   for( seq = 0; seq < IP->nNumSequences; seq++ )
     {
       if( seq != nSeq )
	 {
	   motifCnt = 0;
	   nLen = SequenceLength( B, seq );
	   nStart =  SequenceStartPos( B, seq );
	   for( n = 0; n < nLen; n++ )
	     {
	       lenConserved += RP->dProbCons[seq][n];
#ifdef _DEBUG_
	       if( ! finite( RP->dProbCons[seq][n] ) || ! finite( lenConserved ) )
		 {
		   sprintf( temp, "Infinite conservation probability - seq: %d pos: %d\n", 
			    seq, n );
		   p_error( temp );
		 }
#endif
	       nPos = nStart + n;

	       for(t = 0; t < IP->nNumMotifTypes; t++) 
		 {
		   Pos_t=Pos[t];
		   if( Pos_t[nPos].nMotifStartPos )
		     motifCnt += RP->dProbCons[seq][n];
		 }
	     }
		   
	   modelSitesCons += motifCnt;
	 }
     }

   modelSitesNonCons = TotalNumMotifs( B ) - modelSitesCons;

   modelPseudoCons = RP->dSiteProbCons * dModelPseudo; 
   modelPseudoNonCons = (1 - RP->dSiteProbCons) * dModelPseudo; 
		     
   frac = lenConserved / IP->nSeqLen;
   beta1 =  RP->dSiteProbNonCons;
   beta0 = 1. - RP->dSiteProbNonCons;
   dBGPseudo = weight * (IP->nSeqLen - RP->dPriorSites);

   bgSitesCons = lenConserved;
   bgPseudoCons = (beta1 / (beta1 + beta0)) * dBGPseudo;
   bgSitesNonCons = IP->nSeqLen - lenConserved;
   bgPseudoNonCons = (beta0 / (beta1 + beta0)) * dBGPseudo;
   
   for(t = 0; t < IP->nNumMotifTypes; t++) 
     {
       RP->alphaCons[t] = modelSitesCons + modelPseudoCons;
       RP->alphaNonCons[t] = modelSitesNonCons + modelPseudoNonCons;
       RP->betaCons[t] = bgSitesCons + bgPseudoCons + modelSitesCons;
       RP->betaNonCons[t] = bgSitesNonCons + bgPseudoNonCons + modelSitesNonCons;
     }
}


double SiteProbFromFootprint( Model B, int nPos, int t, int nSeq, 
			      double *siteProb, double *bgProb )
{
  RPType RP;
  double denomCons;
  double denomNonCons;
  
  RP = B->RP;
  
  denomCons = RP->alphaCons[t] + RP->betaCons[t] + 1.0;
  denomNonCons = RP->alphaNonCons[t] + RP->betaNonCons[t] + 1.0;
  
  /*  nPos += MotifWidth( B, t ) / 2;
      nPos = min( nPos , SequenceLength( B, nSeq ) - 1 ); */    /* TEST - 09/29/2000 */

  *siteProb = RP->dProbCons[nSeq][nPos] * ((RP->alphaCons[t] + 1) / denomCons) +  
            (1.0 - RP->dProbCons[nSeq][nPos]) * ((RP->alphaNonCons[t] + 1) / denomNonCons);
  *bgProb = RP->dProbCons[nSeq][nPos] * ((RP->betaCons[t] + 1) / denomCons) +  
            (1.0 - RP->dProbCons[nSeq][nPos]) * ((RP->betaNonCons[t] + 1) / denomNonCons);

  return( *siteProb / *bgProb );
}
