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
/* $Id: maximize.c,v 1.10 2009/04/23 18:43:53 Bill Exp $         */
/*                                                                       */
/* AUTHOR :      Eric C. Rouchka July  8, 1996                           */
/*               Jun Zhu, October, 20, 1996                              */
/*                                                                       */
/* DESCRIPTION : Contains the code used to maximize the map value after  */
/*               near optimal sampling has taken place.  This is         */
/*               accomplished by checking each individual site and seeing*/
/*               whether or not its inclusion into the alignment improves*/
/*               the map.                                                */
/*************************************************************************/

#include "maximize.h"

/*-----------------------------------------------------------------------*/
/*--------------------- LOCAL FUNCTION PROTOTYPES -----------------------*/
/*-----------------------------------------------------------------------*/
MaxResults adjust_max_map  (Model B, PoSition **Pos, double *adj_map); 
void       print_adj_tail  (FILE *fpt, double freq_map, double nearopt_map,
                            double adj_map);
void CopyBayesianProb( Model B, MaxResults max,  MaxResults nearoptmax, PoSition **Pos );

/*------------------------------------------------------------------------*/
/*-------------------------FUNCTION DECLARATIONS -------------------------*/
/*------------------------------------------------------------------------*/

double set_freq_map(Model B, MaxResults optmax)

   /*========================================================================*/
   /* FUNCTION NAME : set_freq_map                                           */
   /*                                                                        */
   /* DESCRIPTION : calculates the map value based on the frequency of times */
   /*               that various motif elements got sampled in               */
   /*========================================================================*/

{
   double freq_map = 0.0;
   int t, i;
   int nSeq;

   for( i = 0; i < B->IP->nNumSequences; i++ )
     B->RP->nSites[i] = 0;
   
   for(t = 0; t < B->IP->nNumMotifTypes; t++) 
     {
       B->IP->nMotifLen[t] = optmax.nMotifLen[t];   /* BT 5/23/97 */
       B->C->nTotBack = B->IP->nSeqLen;
       B->C->dBackCnt = B->First->dBackCnt;
       
       B->IP->nNumMotifs[t][FORWARD] = 0;
       B->IP->nNumMotifs[t][REVERSE] = 0;
       for(i = 0; i < B->IP->nSeqLen; i++) 
	 {	   
	   /* Include those elements which occur with a frequency greater */
	   /* than a predefined cutoff value into the frequency alignment */
	   
	   nSeq = SequenceFromPosition( B,  i );
	   if(optmax.frequency[i][t].frequency > B->IP->dNearOptCutoff &&
	      ((B->IP->is_defined[cl_E] && 
		B->RP->nSites[nSeq] < B->RP->nMaxBlocks) || 
	       (! B->IP->is_defined[cl_E])) ) 
	     {
	       if(optmax.frequency[i][t].fwd_freq > 
		  optmax.frequency[i][t].rev_freq)
		 {
		   adjust_counts(B, ADD, i, t, FALSE);
		   B->IP->nNumMotifs[t][FORWARD]++;
		 }
	       else 
		 {
		   adjust_counts(B, ADD, i,t, TRUE);
		   B->IP->nNumMotifs[t][REVERSE]++;
		 }
	       B->RP->nSites[nSeq]++;
	     }
	 }
     }
   freq_map = CalcMapProb(B, B->IP->is_defined[cl_R]);
   return(freq_map);
}

double set_max_map(Model B, MaxResults optmax, PoSition **Pos,
                   double nearopt_map, double freq_map)

   /*======================================================================*/
   /* FUNCTION NAME : set_max_map                                          */
   /*                                                                      */
   /* DESCRIPTION : Sets the alignment that yields the highest map value   */
   /*               found thus far                                         */
   /*======================================================================*/

{
  int      i, n, t;
   double   adj_map;
   PoSition *Pos_t;  /* Pos_t=Pos[t] */
   int      nSeq;

   copy_counts(B);

   /* Clear all of the motif starting positions  */
   for(t = 0; t < B->IP->nNumMotifTypes; t++)
     { 
       Pos_t=Pos[t];
       for(n = 0; n < B->IP->nSeqLen; n++)       
	 {
	   Pos_t[n].nMotifStartPos = FALSE;
	   Pos_t[n].nInMotif = FALSE;	
	 } 
     }

   if(nearopt_map >= freq_map) {                   /* near optimal is max  */
     for( i = 0; i < B->IP->nNumSequences; i++ )
       B->RP->nSites[i] = 0;

      for(t = 0; t < B->IP->nNumMotifTypes; t++) {/* so set its alignment */
	B->IP->nMotifLen[t] = optmax.nMotifLen[t];   /* BT 5/23/97 */
	for(i = 0; i < optmax.nNumMotifs[t]; i++) {/* as maximal alignment */
	  adjust_counts(B, ADD, optmax.nMotifLoc[i][t], 
			t, optmax.RevComp[i][t]);
	  set_in_motif(Pos,optmax.nMotifLoc[i][t],B,t,
		       optmax.RevComp[i][t]);
	  B->RP->nSites[SequenceFromPosition( B,  optmax.nMotifLoc[i][t] )]++;
	}
	setMotifNum(B->IP, optmax);
      }
      adj_map = nearopt_map;
   }
   else 
     { /* frequency map is max so set max alignment */
       for( i = 0; i < B->IP->nNumSequences; i++ )
	 B->RP->nSites[i] = 0;
       
       for(t = 0; t < B->IP->nNumMotifTypes; t++) 
	 {   
	   B->IP->nMotifLen[t] = optmax.nMotifLen[t];   /* BT 5/23/97 */
	   B->IP->nNumMotifs[t][FORWARD] = 0;
	   B->IP->nNumMotifs[t][REVERSE] = 0;
	   for(i = 0; i < B->IP->nSeqLen; i++) 
	     {
	       nSeq = SequenceFromPosition( B,  i );
	       /* BT 04/14/05 */ /* make sure we don't overflow max blocks */
	       if( optmax.frequency[i][t].frequency > B->IP->dNearOptCutoff &&
		   ((B->IP->is_defined[cl_E] && 
		     B->RP->nSites[nSeq] < B->RP->nMaxBlocks) || 
		    (! B->IP->is_defined[cl_E])) ) 
		 {
		   if(optmax.frequency[i][t].fwd_freq >
		      optmax.frequency[i][t].rev_freq) 
		     {
		       adjust_counts(B,ADD,i,t,FALSE);
		       B->IP->nNumMotifs[t][FORWARD]++;
		       set_in_motif(Pos,i,B,t,FALSE);
		     }
		   else 
		     {
		       adjust_counts(B,ADD,i,t,TRUE);
		       B->IP->nNumMotifs[t][REVERSE]++;
		       set_in_motif(Pos,i,B,t,TRUE);
		     }
		   B->RP->nSites[nSeq]++;
		   /* i+=B->IP->nMotifLen[t]; */
		 }
	     }
	 }
       adj_map = freq_map;
     }
   return adj_map;
}


MaxResults adjust_max_map(Model B, PoSition **Pos,  double *adj_map)

   /*=======================================================================*/
   /* FUNCTION NAME : adjust_max_map                                        */
   /*                                                                       */
   /* DESCRIPTION : Improves the maximal alignment by including motif       */
   /*               motif elements whose inclusion (whether it is forward   */
   /*               or reverse complement) improves upon the maximum        */
   /*               calculated map value                                    */
   /*=======================================================================*/
 
{
  double     map_fwd, map_rev, map_same;
  int        n, t, newpos; 
  MaxResults M; /** Have to set these -- will then be returned **/
  int        last_inc;
  double     dMax, dLocMax;
  double     **dProbArray;
  int        maxnum;
  ProbStruct P;		/* BT 2/21/97 */
  int        start, end;
  int        nSeq;
  int        typ;
  int        p;
  int        len;

  init_maxdata(&M);
  map_fwd = map_rev = map_same = (*adj_map);
  for( nSeq = 0; nSeq < B->IP->nNumSequences; nSeq += SeqInc( B, nSeq ) )
    {
      len = SequenceLength( B, nSeq ); 
      start = SequenceStartPos( B, nSeq );
      end = SequenceEndPos( B, nSeq );
      for(n = start; n <= end; n++) 
	{      
	  for(t = 0; t < B->IP->nNumMotifTypes; t++) 
	    {
	      if(((B->IP->is_defined[cl_E] && 
		   B->RP->nSites[nSeq] < B->RP->nMaxBlocks) || 
		  (! B->IP->is_defined[cl_E])) &&
		 PossibleStartPos( Pos, n, t, len, nSeq, B) &&
		 !OverlapsPrevMotif(n,B->IP->nNumMotifTypes,Pos) &&  
		 !Pos[t][n].nMotifStartPos) 
		{
		  newpos = -1;	/* BT 3/12/97 */
		  if( AnyOverlap(B, n, t, Pos, &newpos, &typ) )
		    {
		      for( p = 0; p < SeqInc( B, nSeq ); p++ )
			{
			  adjust_counts(B,DELETE,newpos + p * len, typ, Pos[typ][newpos + p * len].RevComp);
			  B->IP->nNumMotifs[typ][Pos[typ][newpos + p * len].RevComp]--;
			  B->RP->nSites[nSeq + p]--;
			  not_in_motif(Pos, newpos + p * len,B,typ);
			}
		    } 
		  
		  /* Calculate the map values for the current position */
		  /* being excluded (map_same), included as a forward  */
		  /* motif (map_fwd), or being included as a reverse   */
		  /* complement motif (map_rev).                       */
		  
		  map_same = CalcMapProb(B, B->IP->is_defined[cl_R]); 
		  
		  if( PossibleStartPos( Pos, n, t, len, nSeq, B) )
		    {		  
		      for( p = 0; p < SeqInc( B, nSeq ); p++ )
			{
			  adjust_counts(B, ADD, n + p * len,t,FALSE);
			  B->RP->nSites[nSeq+p]++;	      
			  B->IP->nNumMotifs[t][FORWARD]++;		/* BT 3/12/97 */
			}
		      map_fwd = CalcMapProb(B, B->IP->is_defined[cl_R]);	      
		      for( p = 0; p < SeqInc( B, nSeq ); p++ )
			{
			  adjust_counts(B, DELETE, n + p * len,t,FALSE);
			  B->RP->nSites[nSeq+p]--;
			  B->IP->nNumMotifs[t][FORWARD]--;		/* BT 3/12/97 */
			}
		      
		      if(B->IP->RevComplement) 
			{
			  for( p = 0; p < SeqInc( B, nSeq ); p++ )
			    {
			      adjust_counts(B, ADD, n + p * len, t,TRUE);
			      B->RP->nSites[nSeq+p]++;
			      B->IP->nNumMotifs[t][REVERSE]++;		/* BT 3/12/97 */
			    }
			  map_rev = CalcMapProb(B, B->IP->is_defined[cl_R]);
			  for( p = 0; p < SeqInc( B, nSeq ); p++ )
			    {
			      adjust_counts(B, DELETE, n + p * len, t, TRUE);
			      B->RP->nSites[nSeq+p]--;
			      B->IP->nNumMotifs[t][REVERSE]--;		/* BT 3/12/97 */
			    }
			}
		      
		      /* See if the forward map improves the maximal */
		      /* map calculated so far                       */
		      
		      if((map_fwd >= map_same)  && (map_fwd >= map_rev) &&
			 (map_fwd >= (*adj_map)))
			{
			  for( p = 0; p < SeqInc( B, nSeq ); p++ )
			    {
			      adjust_counts(B, ADD, n + p * len,t, FALSE);
			      B->RP->nSites[nSeq+p]++;
			      B->IP->nNumMotifs[t][FORWARD]++;
			      set_in_motif(Pos,n + p * len,B,t,FALSE);
			    }
			  (*adj_map) = map_fwd;
			  
			  /* as this segment in alignment, we need to check */
			  /* segments do not overlap with this segment */
			  n += (MotifWidth( B, t ) - 1);
			}		
		      
		      /* See if the reverse complement map improves */
		      /* the maximal map calculated so far          */
		      
		      else if(B->IP->RevComplement && 
			      (map_rev >= map_same) && (map_rev >= map_fwd) && 
			      (map_rev >= (*adj_map))) 
			{
			  for( p = 0; p < SeqInc( B, nSeq ); p++ )
			    {
			      adjust_counts(B, ADD, n + p * len,t, TRUE);
			      B->RP->nSites[nSeq+p]++;
			      B->IP->nNumMotifs[t][REVERSE]++;
			      set_in_motif(Pos,n + p * len,B,t,TRUE);
			    }
			  (*adj_map) = map_rev;
			  /* as this segment in alignment, we need to check */
			  /* segments do not overlap with this segment */
			  n += (MotifWidth( B, t ) - 1);
			}     		  
		      /* if we removed a motif and things did not improve when we */
		      /* added n, put the old motif back */
		      /* BT 3/12/97 */ 
		      else if( newpos >= 0 )
			{
			  for( p = 0; p < SeqInc( B, nSeq ); p++ )
			    {
			      adjust_counts(B, ADD,newpos + p * len,typ,Pos[typ][newpos + p * len].RevComp);
			      B->IP->nNumMotifs[typ][Pos[typ][newpos+p*len].RevComp]++;
			      set_in_motif(Pos,newpos + p*len,B,typ, Pos[typ][newpos+p*len].RevComp);
			      B->RP->nSites[nSeq+p]++;
			    }
			}           		
		      else if(Pos[t][n].nMotifStartPos) 
			/* as this segment in alignment, we need to check */
			/* segments do not overlap with this segment */
			n += (MotifWidth( B, t ) - 1);
		    }
		}
	    }
	}
    }
   
   maxnum = findMaxNumMotif(B->IP);
      
   init_prob( &P, B->IP );			/* BT 2/21/97 */
   update_prob( &P, B, TRUE  );
   update_posterior_prob(B);		/* BT 3/17/97 */

   dProbArray = setElementProb( B, Pos, P );

   M = setMaxData(B->F, B->IP, (*adj_map), 1, Pos, &last_inc, &dMax,
                  &dLocMax, M, dProbArray, B);
   FREEP(dProbArray, maxnum); 
   M.frequency = NULL; 

   free_prob( &P, B->IP );
    
   return M;
}


void print_adj_tail(FILE *fpt, double freq_map, double nearopt_map,
                    double adj_map)

   /*===================================================================*/
   /* FUNCTION NAME : print_adj_tail                                    */
   /*                                                                   */
   /* DESCRIPTION : Prints out the map values for the various stages of */
   /*               the Bernoulli sampler                               */
   /*===================================================================*/

{
   fprintf(fpt, "Frequency Map = %f\n", freq_map);
   fprintf(fpt, "Nearopt Map   = %f\n", nearopt_map);
   fprintf(fpt, "Maximal Map   = %f\n", adj_map);
}


void CopyBayesianProb( Model B, MaxResults max,  MaxResults nearoptmax, PoSition **Pos )
{
  int      t;
  int      n;
  int      *count;
  PoSition *Pos_t; 

   NEW(count, B->IP->nNumMotifTypes, int);

   for(n = 0; n < B->IP->nSeqLen; n++)
     { 
       for(t=0; t<B->IP->nNumMotifTypes; t++)
	 {
	   Pos_t=Pos[t];
	   if(Pos_t[n].nMotifStartPos) 
	     {	       
	       if( ! Pos_t[n].RevComp )
		 max.dvMotifProb[count[t]][t] = nearoptmax.frequency[n][t].frequency;
	       else
		 max.dvMotifProb[count[t]][t] = nearoptmax.frequency[n][t].frequency;
	       count[t]++;
	     }
	 }
     }

   free( count );
  
}


MaxResults maximize_map(MaxResults optmax, Model B, PoSition **Pos)  /* BT 2/21/97 */

   /*=================================================================*/
   /* FUNCTION NAME : maximize_map                                    */
   /*                                                                 */
   /* DESCRIPTION : maximizes the map by sampling in individual       */
   /*               residues that increase the map value              */
   /*=================================================================*/

{
   double nearopt_map, adj_map, freq_map;
   MaxResults M;

   copy_counts(B);
   nearopt_map = optmax.dProbability;   
   freq_map    = set_freq_map(B, optmax);
   adj_map     = set_max_map(B, optmax,  Pos, nearopt_map, freq_map);
   M = adjust_max_map(B, Pos, &adj_map); 
   M.nSuboptSeed = optmax.nSuboptSeed;
   CopyBayesianProb( B, M, optmax, Pos );
   print_info(B, M, TRUE, MAPMAX );
   if( ! B->IP->inCentroidAlign )
     print_adj_tail(B->IP->Datafiles->out_fpt, freq_map, nearopt_map, adj_map);

   M.frequency = optmax.frequency; /* BT 2/21/97 */
   optmax.frequency = NULL;

   free_maxdata( &optmax, B->IP );
   
   return M;					/* BT 2/21/97 */
}

