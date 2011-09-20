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
/******************************************************************************/
/* FILE NAME :                                                                */
/* AUTHOR    : Eric C. Rouchka                                                */
/* LAST UPDATE : July 22, 1996                                                */
/******************************************************************************/

#include "calcstats.h"

double Gser(double A, double X);
double GCF(double A, double X);


/* Taken from Numerical Recipies, pp 160 - 165 */

double Pval(double dof, double chisqr)

   /* Calculates the incomplete gamma function P(A, X) */
{
   double Gammp;
   double A, X;

   printf("dof = %f chisqr = %f\n", dof, chisqr);
   A = dof / 2.0; X = chisqr / 2.0;
   if((X < 0.0) || (A <= 0.0))
      perror("ERROR in finding Pvalue");
   if(X < A + 1.0)                        /* Use the Series Representation */
      Gammp = Gser(A, X);
   else                                   /* Use comeplement of the continued */
      Gammp = (1.0 - GCF(A, X));          /* fraction representation          */

   return 1.0 - Gammp;
}

double Gser(double A, double X)

   /* Returns the incomplete gamma function P(a, x) evaluated */
   /* by its series representation as gamser                  */
{
   double gln;
   double AP, sum, del;
   int n;

   gln = ln_gamma(A);
   if(X <= 0.0) {
      if(X < 0.0)
         perror("ERROR in Gser::X < 0.0");
      return 0.0;
   }
   AP = A;
   del = sum = 1.0 / A;
   for(n = 1; n <= ITMAX; n++) {
      AP++;
      del *= (X/AP);
      sum += del;
      if(fabs(del) < fabs(sum)*EPS) 
         return sum * exp(-X + A * log(X) - gln);
   }
   fprintf(stderr, "WARNING :: in Gser : A to large, itmax too small\n");
   return (sum * exp(-X + A * log(X) - gln));
}

double GCF(double A, double X)
{
   double gln, Gold, G=0;
   double A0, A1, B0, B1, AN, ANA, FAC, ANF;
   int n;

   gln = ln_gamma(A);
   Gold = B0 = 0.0;           
   A0 = B1 = FAC = 1.0;
   A1 = X;
   for(n = 1; n <= ITMAX; n++) {
      AN = (double)n;
      ANA = AN - A;
      A0 = (A1 + A0 * ANA) * FAC;
      B0 = (B1 + B0 * ANA) * FAC;
      ANF = AN * FAC;
      A1 = X * A0 + ANF * A1;
      B1 = X * B0 + ANF * B1;
      if (A1 != 0.0) {
         FAC = 1.0 / A1;
         G = B1 * FAC;
         if(fabs((G - Gold) / G) < EPS) 
            return (exp(-X + A * log(X) - gln) * G);
         Gold = G;
      }
   }
   fprintf(stderr, "WARNING :: in Gser : A too large, ITMAX too small\n");
   return(exp(-X + A * log(X) -gln) * G);
}
