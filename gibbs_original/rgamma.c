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
/* 
 * File:   rgamma.c
 * Author: thompson
 *
 * Created on June 23, 2004, 2:43 PM
 */

#include "rgamma.h"

static unsigned int iz;
static unsigned int jz;
static unsigned int jsr = 349576937;
static unsigned int kn[128];
static int hz1;
static double wn[128];
static double fn[128];

double nfix();

#define SHR3 (jz=jsr, jsr ^= (jsr << 13), jsr ^= (jsr >> 17), jsr ^= (jsr << 5), jz+jsr)
#define UNI  (SHR3 * 0.2328306e-9)
#define RNOR (hz1=SHR3, iz=hz1 & 127, (abs(hz1) < kn[iz]) ? hz1 * wn[iz] : nfix())

/*
  rgamma - Generate random samples from a gamma distribution
  based on Marsaglia, G. and Tsang W. W. "A simple Method for Generating Gamma Variable", 
  ACM Transactions on Mathematical Software, Vol 26, No 3, Sept. 2000, p 363-372

 */


double rgamma( double a) {
    double d;   
    double c;
    double x;
    double v;
    double u;

    if( a < 0.0 )
      p_internal_error( "rgamma: negative shape parameter" );
    
    if( a < 1.0 )
        return pow( UNI, 1.0 / a ) * rgamma( 1.0 + a );
    
    d = a - 1.0/3.0;
    c = 1.0 / sqrt( 9.0 * d );
    while( TRUE ) {
        do {
            x = RNOR;
            v= 1.0 + c * x;
        } while( v <= 0.0 );
        
        v = v * v * v;
        u = UNI;
        if( u < 1.0 - 0.0331 * (x * x) * (x * x) )
            return( d * v );
        if( log( u ) < 0.5* x * x + d * (1.0 - v + log( v ) ) )
            return( d * v );
    }
}

double nfix() {
    const double r = 3.713086;
    static double x;
    static double y;
    
    while( TRUE ) {
        x = hz1 * wn[iz];
        if( iz == 0 ) {
            if( UNI * 0.002669629 < exp( - 0.5 * x * x) )
                return x;
            do {
                x = -log( UNI ) * 0.2693178;
                y = -log( UNI );
            } while( y + y < x * x );
            return( (hz1 > 0) ? r + x : -r - x );
        }
        if( fn[iz] + UNI * (fn[iz-1] - fn[iz]) < exp( -0.5 * x * x ) )
            return x;
        hz1 = SHR3;
        iz = hz1 & 127;
        if( abs( hz1 ) < kn[iz] )
            return( hz1 * wn[iz] );
    }
}

void zigset( unsigned int jsrseed ) {
    const double m1 = 2147483648.0;
    double dn = 3.442619855899;
    double tn = dn;
    double vn = 9.91256303526217e-3;
    double q;
    int i;
    
    jsr = jsrseed;
    q = vn / exp( -0.5 * dn * dn );
    kn[0] = (dn / q) * m1;
    kn[1] = 0;
    wn[0] = q / m1;
    wn[127] = dn / m1;
    fn[0] = 1.0;
    fn[127] = exp( - 0.5 * dn * dn );
    for( i = 126; i >= 1; i-- ) {
        dn = sqrt( -2.0 * log( vn/dn + exp( -0.5 * dn * dn) ) );
        kn[i+1] = (dn / tn) * m1;
        tn = dn;
        fn[i] = exp( -0.5 * dn * dn);
        wn[i] = dn / m1;
    }
}
