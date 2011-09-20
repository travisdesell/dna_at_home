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
/*************************************
 A C implementation of the digamma-function based
 on the Chebyshev expansion proposed in appendex E of 
 http://arXiv.org/abs/math.CA/0403344

Richard J. Mathar, 2005-05-15
**************************************/

/* downloaded from http://www.strw.leidenuniv.nl/~mathar/progs/ */


#include "digammal.h"

 /* constants: pi, Euler's constant, and log(2) */
#ifndef M_PIl
#define M_PIl 3.1415926535897932384626433832795029L
#endif
#ifndef M_GAMMAl
#define M_GAMMAl 0.5772156649015328606065120900824024L
#endif
#ifndef M_LN2l
#define M_LN2l 0.6931471805599453094172321214581766L
#endif

#ifdef _SOLARIS_
extern long double tanl(long double x);
#endif


double digamma( double x )
{
  return (double) digammal( (long double) x );
}


long double digammal(long double x)
{
  /* force into the interval 1..3 */
  if( x < 0.0L )
    return digammal(1.0L-x)+M_PIl/tanl(M_PIl*(1.0L-x)) ;	/* reflection formula */
  else if( x < 1.0L )
    return digammal(1.0L+x)-1.0L/x ;
  else if ( x == 1.0L)
    return -M_GAMMAl ;
  else if ( x == 2.0L)
    return 1.0L-M_GAMMAl ;
  else if ( x == 3.0L)
    return 1.5L-M_GAMMAl ;
  else if ( x > 3.0L)
    /* duplication formula */
    return 0.5L*(digammal(x/2.0L)+digammal((x+1.0L)/2.0L))+M_LN2l ;
  else
    {
      static long double Kncoe[] = { .30459198558715155634315638246624251L,
				     .72037977439182833573548891941219706L, -.12454959243861367729528855995001087L,
				     .27769457331927827002810119567456810e-1L, -.67762371439822456447373550186163070e-2L,
				     .17238755142247705209823876688592170e-2L, -.44817699064252933515310345718960928e-3L,
				     .11793660000155572716272710617753373e-3L, -.31253894280980134452125172274246963e-4L,
				     .83173997012173283398932708991137488e-5L, -.22191427643780045431149221890172210e-5L,
				     .59302266729329346291029599913617915e-6L, -.15863051191470655433559920279603632e-6L,
				     .42459203983193603241777510648681429e-7L, -.11369129616951114238848106591780146e-7L,
				     .304502217295931698401459168423403510e-8L, -.81568455080753152802915013641723686e-9L,
				     .21852324749975455125936715817306383e-9L, -.58546491441689515680751900276454407e-10L,
				     .15686348450871204869813586459513648e-10L, -.42029496273143231373796179302482033e-11L,
				     .11261435719264907097227520956710754e-11L, -.30174353636860279765375177200637590e-12L,
				     .80850955256389526647406571868193768e-13L, -.21663779809421233144009565199997351e-13L,
				     .58047634271339391495076374966835526e-14L, -.15553767189204733561108869588173845e-14L,
				     .41676108598040807753707828039353330e-15L, -.11167065064221317094734023242188463e-15L } ;
      
      register long double Tn_1 = 1.0L ;	/* T_{n-1}(x), started at n=1 */
      register long double Tn = x-2.0L ;	/* T_{n}(x) , started at n=1 */
      register long double resul = Kncoe[0] + Kncoe[1]*Tn ;
      int                  n;
      
      x -= 2.0L ;
      
      for(n = 2 ; n < sizeof(Kncoe)/sizeof(long double) ;n++)
	{
	  const long double Tn1 = 2.0L * x * Tn - Tn_1 ;	/* Chebyshev recursion, Eq. 22.7.4 Abramowitz-Stegun */
	  resul += Kncoe[n]*Tn1 ;
	  Tn_1 = Tn ;
	  Tn = Tn1 ;
	}
      return resul ;
    }
}

