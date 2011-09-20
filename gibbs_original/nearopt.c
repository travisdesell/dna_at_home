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
/* $Id: nearopt.c,v 1.11 2009/04/23 18:43:54 Bill Exp $           */
/*                                                                        */
/* AUTHOR    : Eric C. Rouchka, July 8, 1996                              */
/*             Jun Zhu, October 20, 1996                                  */
/**************************************************************************/

#include "nearopt.h"


int compare( const void *v1, const void *v2 );


/*------------------------------------------------------------------------*/
/*----------------------------- FUNCTION DECLARATIONS --------------------*/
/*------------------------------------------------------------------------*/
void FindGoodSites(Model B, int ****occurence, MaxResults maxData, 
                   PoSition **Pos, int ***good, ProbStruct P, int seed_run)
{
   int    i, t, j=0;
   int    k;
   int    ngood=0, nmap=0, nbad=0;
   double prob = 0.0;
   short  foundIt;	/* BT 4/4/97 */ 
   int    ndiff;
   int    nSeq;
   int    nOffset = 0;
   double ***posProb;
   double denom;
   int    first;
   Mlist  M;  /* BT 11/12/99 */
   double dCutoff;
   int    p;
   
   NEWP((*good), B->IP->nSeqLen, int);                /* ALLOCATE MEMORY */
   for(i = 0; i < B->IP->nSeqLen; i++)
     NEW((*good)[i], B->IP->nNumMotifTypes, int);

   NEWPP( posProb, B->IP->nSeqLen, double);
   for( i = 0; i  < B->IP->nSeqLen; i++ ) 
     {
       NEWP( posProb[i], B->IP->nNumMotifTypes, double);
       for(t = 0; t < B->IP->nNumMotifTypes; t++) 
	 {
	   NEW( posProb[i][t], 2, double);
	 }		
     }

   /* BT 11/12/99 */
   copy_counts(B);
   for( i = 0; i < B->IP->nNumSequences; i++ )
     B->RP->nSites[i] = 0;
   
   M = initializeMotifs(B->IP);                    /* now add in the */
   for(t = 0; t < B->IP->nNumMotifTypes; t++)      /* motifs in the  */
     {                                             /* max alignment  */
       B->IP->nNumMotifs[t][FORWARD] = 0;
       B->IP->nNumMotifs[t][REVERSE] = 0;
       for(i = 0; i < maxData.nNumMotifs[t]; i++) 
	 { 
	   adjust_counts(B, ADD, maxData.nMotifLoc[i][t], t, maxData.RevComp[i][t]);
	   add_motif(B->IP, B->Seq, maxData.nMotifLoc[i][t],M,t, 
		     maxData.RevComp[i][t]);
	   set_in_motif(Pos,maxData.nMotifLoc[i][t],B,t,
			maxData.RevComp[i][t]);
	   B->IP->nNumMotifs[t][maxData.RevComp[i][t]]++;
	   B->RP->nSites[SequenceFromPosition( B,  maxData.nMotifLoc[i][t] )]++;
	 }
     }
   update_prob(&P, B, TRUE); 
   update_posterior_prob( B );  /* BT 04/06/04 */
   
   if(  B->IP->is_defined[cl_E] )
     {
       if( B->IP->is_defined[cl_Q] && seed_run == 0 )
	 {
	   /* average counts over all motif types because there's no way to guarantee */
	   /* that the motif types are the same between seeds                         */
	   for( k = 0; k < B->IP->nSeeds; k++ )
	     {
	       for( nSeq = 0; nSeq < B->IP->nNumSequences; nSeq++ )
		 {
		   first = SequenceStartPos( B, nSeq );
		   for( t = 0; t < B->IP->nNumMotifTypes; t++) 
		     {
		       for( i = 0; i < SequenceLength( B, nSeq ); i++ )
			 {
			   posProb[i+first][0][FORWARD] += 
			     ((double) B->IP->nAlignCnts[k][t][nSeq][i]);
			 }
		     }
		 }
	     }

	   denom = B->IP->nPostPlateauIter * B->IP->nSeeds;
	   for( i = 0; i < B->IP->nSeqLen; i++ ) 
	     {
	       posProb[i][0][FORWARD] /= denom;
	       for( t = 1; t < B->IP->nNumMotifTypes; t++) 
		 {
		   posProb[i][t][FORWARD] =  posProb[i][0][FORWARD];
		 }
	     }
	 }
       else
	 {
	   for( i = 0; i <  B->IP->nSeqLen; i++ )
	     CalcPosProb( B, P, Pos, i, posProb[i] );
	 }
     }
   else if( B->IP->is_defined[cl_D] )
     {
       for( nSeq = 0; nSeq < B->IP->nNumSequences; nSeq += SpeciesInc( B, nSeq ) )
	 {
	   first = SequenceStartPos( B, nSeq );
	   nOffset = SequenceLength( B, nSeq );
	   for( i = first; i < nOffset + first - 1; i++ )
	     {
	       for( t = 0; t < B->IP->nNumMotifTypes; t++ )
		 { 
		   if(PossFragStartPos(Pos, i, t, B) &&
		      PossFragStartPos(Pos, i + nOffset, t, B) &&
		      !OverlapsPrevMotif(i, B->IP->nNumMotifTypes, Pos) &&
		      !OverlapsPrevMotif(i + nOffset, B->IP->nNumMotifTypes, Pos))
		     {
		       CalcPosProb( B, P, Pos, i, posProb[i] );
		       CalcPosProb( B, P, Pos, i + nOffset, posProb[i+nOffset] );
		     }
		 }	       
	     }
	 }
     }
   else
     {
       for( i = 0; i <  B->IP->nSeqLen; i++ )
	 CalcPosProb( B, P, Pos, i, posProb[i] );
     }

   
#ifdef _DEBUG_ 
   for( t = 0; t < B->IP->nNumMotifTypes; t++ )
     {
       fprintf( stdout, 
		"%d Motif =  %.5f Frag =  %.5f\n", t, 
		CalcMotifMap(B, t, B->IP->is_defined[cl_R]),
		CalcMotifFragMap(B, t, B->IP->is_defined[cl_R]) );
     }
   fprintf( stdout, 
	    "Bkgnd =  %.5f Beta =  %.5f Seq = %.5f Null = %.5f\n",  
	    CalcBkgndMap(B, B->IP->is_defined[cl_R]),
	    CalcBetaMap( B, B->IP->is_defined[cl_R]),
	    CalcSitePerSeqMap( B ), 
	    B->IP->dnull_map );	  
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
	    0, CalcMapProb(B, B->IP->is_defined[cl_R]));
#endif 

   dCutoff = B->IP->dCutoff;

   for(t = 0; t < B->IP->nNumMotifTypes; t++) 
     {
       i = 0;
       while( i < B->IP->nSeqLen )    /* BT 9/9/98 */
	 {
	   nSeq = SequenceFromPosition( B, i );
	   nOffset = SequenceLength( B, nSeq );     
	   
	   Pos[t][i].nMotifStartPos = FALSE;
	   Pos[t][i].nInMotif=FALSE;
	   (*good)[i][t]=TRUE;
	   if(maxData.nNumMotifs[t] != 0) 
	     {
	       /* BT 4/4/97 */
	       for( foundIt = FALSE, j = 0; j < maxData.nNumMotifs[t]; j++ ) 
		 {
		   if(i == maxData.nMotifLoc[j][t]) 
		     {
		       foundIt = TRUE; 
		       (*good)[i][t] = TRUE;
		       if(B->IP->is_defined[cl_t])
			 {
			   prob = posProb[i][t][FORWARD] + posProb[i][t][REVERSE];
			   nSeq = SequenceFromPosition( B, i );
			   fprintf(B->IP->Datafiles->out_fpt, "%3d     %d   %d    %d   %g\n", t, i, 
				   nSeq, i - SequenceStartPos( B, nSeq ), prob);
			 }

		       ngood++;
		       nmap++;
		       if(!B->IP->is_defined[cl_F]) 
			 {	/* BT 9/2/98 */
			   ndiff = B->F->FragWidth[t] - 
			     (int) ((1.0 - B->IP->glOverlapParam) * 
				    B->IP->nMotifLen[t]) + 1; 
			   for(k = 1; k < ndiff; k++ )
			     {
			       (*good)[i+k][t] = TRUE;
			       if(B->IP->is_defined[cl_t])
				 {
				   prob = posProb[i + k][t][FORWARD] + posProb[i + k][t][REVERSE];
				   nSeq = SequenceFromPosition( B, i + k );
				   fprintf(B->IP->Datafiles->out_fpt, "%3d     %d   %d    %d   %g\n", t, i + k, 
					   nSeq, i + k - SequenceStartPos( B, nSeq ), prob);
				 }

			       ngood++;
			       nmap++;			   
			     }
			   i += ndiff-1;
			 }       
		       break;
		     }
		 }
	       
	       if( !foundIt )
		 { 
		   if( PossFragStartPos(Pos, i, t, B) ) 
		     {
		       (*good)[i][t] = TRUE;
		       if( ! IsPhyloSeq( B, nSeq ) ||
			   nSeq == B->Phylo->phyloTreeStartSeq[B->Phylo->phyloIndex[nSeq]] )
			 {
			   /* prob=in_motif_prob(B, P, i, t, FALSE, FALSE, &dMotifProb, &dBGProb);
			   if(B->IP->RevComplement)
			     prob += in_motif_prob(B, P, i, t, TRUE, FALSE,  &dMotifProb, &dBGProb);
			     */

			   prob = posProb[i][t][FORWARD] + posProb[i][t][REVERSE];
			 
			   /* if (prob > B->IP->dCutoff || B->IP->is_defined[cl_E] ) */
			   if (prob >= dCutoff ) /* BT 11/12/99 */
			     {
			       ngood++; (*good)[i][t] =TRUE; 
			       if(B->IP->is_defined[cl_t])
				 {
				   nSeq = SequenceFromPosition( B, i );
				   fprintf(B->IP->Datafiles->out_fpt, "%3d     %d   %d    %d   %g\n", t, i, 
					   nSeq, i - SequenceStartPos( B, nSeq ), prob);
				 }

			       if( IsPhyloSeq( B, nSeq) )
				 {
				   for( p = 1; p < SpeciesInc( B, nSeq ); p++ )
				     {
				       ngood++; 
				       (*good)[i + p * nOffset][t] =TRUE; 
				       if(B->IP->is_defined[cl_t])
					 {
					   prob = posProb[i + p * nOffset][t][FORWARD] + 
					     posProb[i + p * nOffset][t][REVERSE];
					   nSeq = SequenceFromPosition( B, i + p * nOffset );
					   fprintf(B->IP->Datafiles->out_fpt, "%3d     %d   %d    %d   %g\n", t, 
						   i + p * nOffset, 
						   nSeq, i + p * nOffset - SequenceStartPos( B, nSeq ), prob);
					 }
				     }
				 }
			     }
			   else 
			     {
			       (*good)[i][t] = FALSE; 
			       nbad++;
			       if( IsPhyloSeq( B, nSeq ) )
				 {
				   for( p = 1; p < SpeciesInc( B, nSeq ); p++ )
				     {
				       nbad++; 
				       (*good)[i + p * nOffset][t] = FALSE; 
				     }
				 }
			     }
			 }
		     }
		   else 
		     {
		       (*good)[i][t] = FALSE; 
		       nbad++;
		     }
		 }
	     }
	   else if(PossFragStartPos(Pos, i, t, B)) 
	     {
	       if( ! IsPhyloSeq( B, nSeq ) || 
		   nSeq == B->Phylo->phyloTreeStartSeq[B->Phylo->phyloIndex[nSeq]] )
		 {
		   /*prob=in_motif_prob(B, P, i, t, FALSE, FALSE, &dMotifProb, &dBGProb);
		   if(B->IP->RevComplement)
		     prob += in_motif_prob(B, P, i, t, TRUE, FALSE,  &dMotifProb, &dBGProb);
		     */

		   prob = posProb[i][t][FORWARD] + posProb[i][t][REVERSE];

		   /* if (prob > B->IP->dCutoff || B->IP->is_defined[cl_E] ) */
		   if (prob >= dCutoff) /* BT 11/12/99 */ 
		     {
		       ngood++; (*good)[i][t] =TRUE; 
		       nSeq = SequenceFromPosition( B, i );
		       if(B->IP->is_defined[cl_t])
			 {
			   prob = posProb[i][t][FORWARD] + posProb[i][t][REVERSE];
			   fprintf(B->IP->Datafiles->out_fpt, "%3d     %d   %d    %d   %g\n", t, i, 
				   nSeq, i - SequenceStartPos( B, nSeq ), prob);
			 }
		       if( IsPhyloSeq( B, nSeq ) )
			 {
			   for( p = 1; p < SpeciesInc( B, nSeq ); p++ )
			     {
			       ngood++; 
			       (*good)[i + p * nOffset][t] =TRUE; 
			       if(B->IP->is_defined[cl_t])
				 {
				   nSeq = SequenceFromPosition( B, i + p * nOffset);
				   prob = posProb[i][t][FORWARD] + posProb[i][t][REVERSE];
				   fprintf(B->IP->Datafiles->out_fpt, "%3d     %d   %d    %d   %g\n", 
					   t, i + p * nOffset, 
					   nSeq, i + p * nOffset - SequenceStartPos( B, nSeq ), prob);
				 }
			     }
			 }
		     }
		   else 
		     {
		       (*good)[i][t] = FALSE; 
		       nSeq = SequenceFromPosition( B, i );
		       nbad++;
		       if( IsPhyloSeq( B, nSeq ) )
			 {
			   for( p = 0; p < SpeciesInc( B, nSeq ); p++ )
			     {
			       nbad++; 
			       (*good)[i + p * nOffset][t] =FALSE; 
			     }
			 }
		     }
		 }	       
	     }
	   else {
	     (*good)[i][t] = FALSE; 
	     nbad++;
	   }
	   i++;
	 }
     }
#ifndef _MPI_
   if( ! B->IP->is_defined[cl_Z] &&  ! B->IP->is_defined[cl_nopt] )
     fprintf(B->IP->Datafiles->out_fpt, "MAP = %d maybe = %d discard = %d\n",
	     nmap, ngood, nbad);
#endif

   FREEPP( posProb, B->IP->nSeqLen, B->IP->nNumMotifTypes );
   free_motifs( B, M );

   NEWPP((*occurence), B->IP->nSeqLen, int);
   for( i = 0; i  < B->IP->nSeqLen; i++ )     /* BT 9/9/98 */
     {
       NEWP((*occurence)[i], B->IP->nNumMotifTypes, int);
       for(t = 0; t < B->IP->nNumMotifTypes; t++) 
	 {
	   NEW((*occurence)[i][t], 3, int);
	 }		
     }
}

  
MaxResults NearOptSampler(Model B, MaxResults maxData, PoSition **Pos)

   /*======================================================================*/
   /* FUNCTION NAME : NearOptSampler                                       */
   /* DESCRIPTION   : This function finds the near optimal alignment for   */
   /*                 the current run                                      */
   /*======================================================================*/
{
  int n, i, t, j;
  int **good, ***occurence;
  ProbStruct P;
  MaxResults currMax;
  Mlist M;
  int    nPos;
  
  init_maxdata(&currMax);
  P.update = TRUE;
  init_prob(&P, B->IP);
  update_prob(&P, B, TRUE); 
  if(B->IP->is_defined[cl_t]) 
    {
      fprintf(B->IP->Datafiles->out_fpt, "Near Optimal Sampling sites\n");
      fprintf(B->IP->Datafiles->out_fpt, "___________________________\n");
      fprintf(B->IP->Datafiles->out_fpt, "Motif   Site\n");
    }
  
  if( B->IP->is_defined[cl_E] && (! B->IP->is_defined[cl_F] || B->IP->is_defined[cl_d] ) )
    CountAlignments( B, Pos ); 
  
  B->IP->nSeedRun = 0;

  FindGoodSites(B, &occurence, maxData, Pos, &good, P, 0);
  
  copy_counts(B);
  for( i = 0; i < B->IP->nNumSequences; i++ )
    B->RP->nSites[i] = 0;
  
  M = initializeMotifs(B->IP);                    /* now add in the */
  for(t = 0; t < B->IP->nNumMotifTypes; t++)      /* motifs in the  */
    {                                             /* max alignment  */
      B->IP->nNumMotifs[t][FORWARD] = 0;
      B->IP->nNumMotifs[t][REVERSE] = 0;
      for(i = 0; i < maxData.nNumMotifs[t]; i++) 
	{ 
	  adjust_counts(B, ADD, maxData.nMotifLoc[i][t], t, maxData.RevComp[i][t]);
	  add_motif(B->IP, B->Seq, maxData.nMotifLoc[i][t],M,t, 
		    maxData.RevComp[i][t]);
	  set_in_motif(Pos,maxData.nMotifLoc[i][t],B,t,
		       maxData.RevComp[i][t]);
	  B->IP->nNumMotifs[t][maxData.RevComp[i][t]]++;
	  B->RP->nSites[SequenceFromPosition( B,  maxData.nMotifLoc[i][t] )]++;
	}
    }
  update_prob(&P, B, TRUE);
  currMax.dvMotifProb = NULL;
    
  if(B->IP->site_samp)
    currMax = nearopt_site_samp(B, Pos, M, good, &P, &occurence);
  else if( B->IP->is_defined[cl_E] )
    currMax = nearopt_rsite_samp(B, Pos, M, good, &P, &occurence, maxData.dProbability);
  else
    currMax = nearopt_motif_samp(B, Pos, M, good, &P, &occurence);
  
  free_motifs(B, M);
  NEWP(currMax.frequency, B->IP->nSeqLen, Frequency);    
  for(n = 0; n < B->IP->nSeqLen; n++) 
    {
      NEW(currMax.frequency[n], B->IP->nNumMotifTypes, Frequency);
      for(t = 0; t < B->IP->nNumMotifTypes; t++)
	{
	  currMax.frequency[n][t].frequency = (double)occurence[n][t][0] / 
	    (double)B->IP->nMaxIterations;
	  currMax.frequency[n][t].fwd_freq = (double)occurence[n][t][1];
	  currMax.frequency[n][t].rev_freq = (double)occurence[n][t][2];
	}       
    }
  currMax.nSuboptSeed = maxData.nSuboptSeed;
  free_prob(&P, B->IP);
  FREEPP(occurence, B->IP->nSeqLen, B->IP->nNumMotifTypes);
  copy_counts(B);
  for(t = 0; t < B->IP->nNumMotifTypes; t++)
    {
      B->IP->nMotifLen[t] = currMax.nMotifLen[t];     /* restore best motif len */  /* BT 5/28/97 */
      for(i = 0; i < currMax.nNumMotifs[t]; i++)
	{
	  adjust_counts(B, ADD, currMax.nMotifLoc[i][t],
			t, currMax.RevComp[i][t]);
	}
      setMotifNum(B->IP, maxData);
    }   
  
  if( (! B->IP->is_defined[cl_nopt] || ! B->IP->is_defined[cl_y]) && ! B->IP->inCentroidAlign )
    print_info(B, currMax, ! B->IP->is_defined[cl_nopt], NEAROPT );
    
  if( B->IP->is_defined[cl_Q] )     /* BT 3/27/98 */
    {
      fprintf( B->IP->Datafiles->near_occur_fpt, "seq pos pos2 motif count\n" );
      for( i = 0; i < B->IP->nNumSequences; i++ )
	{
	  for( j = 0; j < SequenceLength( B, i ); j++ )
	    {
	      nPos =  SequenceStartPos( B, i ) + j;
	      for(  t = 0;  t < B->IP->nNumMotifTypes; t++ )
		{
		  fprintf( B->IP->Datafiles->near_occur_fpt, "%5d %5d %5d %5d %5d\n",
			   i, j, nPos, t,  B->IP->nAlignCnts[0][t][i][j] );
		}
	    }
	  fflush( B->IP->Datafiles->near_occur_fpt );
	}
    }
  
  return(currMax);
}


int compare( const void *v1, const void *v2 )
{
  if( *((double *) v1) > *((double *) v2) )
    return( -1 );
  else if( *((double *) v1) < *((double *) v2) )
    return( 1 );
  else
    return( 0 );
}
