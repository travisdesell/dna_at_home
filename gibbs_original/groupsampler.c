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
#include "groupsampler.h"

void RemoveGroup( Model B, int grp, PoSition **Pos, Mlist M,  ProbStruct *P );
int FirstGroupSeq( Model B, int grp );
void AddGroup( Model B, int grp, int n, int t, int type, PoSition **Pos, Mlist M, 
	       ProbStruct *P );
int MinGroupSeqLength( Model B, int grp );
void InitGroups( Model B, PoSition **Pos, Mlist M, ProbStruct *P );
double BayesFactor( Model B,  int *dir, int *n, int t, 
		    PoSition **Pos, Mlist M, ProbStruct *P,
		    double *p1Model, double *p2Model );
void PrintGroupInfo( Model B, double prob, double bf, PoSition **Pos, ProbStruct *P,
		     double p1Model, double p2Model );
int PossGroupFragStartPos(PoSition **Pos, int grp, int n, int t, Model B );
double KLDistance( Model B,  int *dir, int *n, int t, 
		   PoSition **Pos, Mlist M, ProbStruct *P,
		   int grp1, int grp2,
		   double *kldist12, double *kldist21 );


#define LOG2( x )  (log((x)) / log( 2.0 ))


void GroupSiteSampler( Model B )
{
  RPType      RP;
  IPtype      IP;
  int         nGrp;
  int         first;
  int         last;
  int         t;
  double      bfTotal = 0.0;
  double      **prob;
  double      totalProb;
  double      dtot;
  double      probTotal = 0.0;
  int         probCnt = 0;
  PoSition    *Pos_t;
  PoSition    **Pos;
  Mlist       M; 
  int         **nNumMotifs;
  int         maxnum;
  int         i;
  int         seq;
  int         n;
  ProbStruct  P;
  double      maxProb = -DBL_MAX;
  double      maxbf = -DBL_MAX;  
  int         *ps;
  int         *dir;
  int         dr;
  double      p1Total = 0.0;
  double      p2Total = 0.0;
  int         found;
  int         index;
  double      sum;
  int         *nGrpLen;

#ifdef _DEBUG_
  double      bf;
  double      p1Model;
  double      p2Model;
#endif

  RP = B->RP;
  IP = B->IP;

  NEWP( Pos, IP->nNumMotifTypes, PoSition);
  for(t=0; t < IP->nNumMotifTypes;t++)
    NEW(Pos[t], B->IP->nSeqLen, PoSition);

  NEW( ps, RP->nGroupCnt, int );
  NEW( dir, RP->nGroupCnt, int );
  
  NEW( nGrpLen, RP->nGroupCnt, int );
  for( nGrp = 0; nGrp < RP->nGroupCnt; nGrp++ )
    nGrpLen[nGrp] = MinGroupSeqLength( B, nGrp );

  copy_counts(B);  
  nNumMotifs = copy_motif_num(IP);
  zero_motifs(IP); /* set motif count to zero */

   /* calculate the probability of NULL model */
  IP->dnull_map = CalcMapProb(B, B->IP->is_defined[cl_R]);
  reset_motif_num(nNumMotifs,IP);
  
   /* mark all possible motif sites */
  set_indicator_vector(Pos, B->Seq, B);
  maxnum = findMaxNumMotif(IP);
  
  CountAlignments( B, Pos ); 
  SaveAlignmentCounts( B );
  
  for(t=0;t<IP->nNumMotifTypes;t++)
    {
      Pos_t=Pos[t];
      for(i = 0; i < IP->nSeqLen; i++) 
	{
	  Pos_t[i].nInMotif=FALSE;
	  Pos_t[i].nMotifStartPos=FALSE;
	}
    }

  M = initializeMotifs( IP );
  init_prob( &P, IP ); 

  InitGroups( B, Pos, M, &P );
  /*  set_posterior_prob(B->IP, B->C); */

  for( nGrp = 0; nGrp < RP->nGroupCnt; nGrp++ )
    {
      RemoveGroup( B, nGrp, Pos, M, &P );
      seq = FirstGroupSeq( B, nGrp );      
      first = SequenceStartPos( B, seq );

      NEWP( prob, nGrpLen[nGrp], double );
      for( i = 0; i < RP->nGroupCnt; i++ )
	NEW( prob[i], 2, double );	    
    	      
      for (t = 0; t < IP->nNumMotifTypes; t++) 
	{
	  last = SequenceEndPos( B, seq ) - IP->nMotifLen[t] + 1;
	  totalProb = 0.0;
	  
	  for( n = first; n <= last; n++ )
	    {
	      if(PossGroupFragStartPos(Pos, nGrp, n - first, t, B) )
		{
		  for( dr = FORWARD; 
		       dr <= REVERSE * ( IP->RevComplement ? 1 : 0); 
		       dr++ )
		    {
		      AddGroup( B, nGrp, n - first, t, dr, Pos, M, &P );
		      dir[nGrp] = dr;
		      ps[nGrp] = n - first;
		      prob[n - first][dr]  = CalcMapProb(B, B->IP->is_defined[cl_R]);
		      totalProb += prob[n - first][dr];
		      RemoveGroup( B, nGrp, Pos, M, &P );
		    }
		}
	    }
	  
	  dtot = totalProb * drand();
	  found = FALSE;
	  index = 0;
	  sum = 0.0;
	  while( (! found) && (index < nGrpLen[nGrp])  );
	  {
	    sum += prob[index][FORWARD];
	    if( sum >= dtot )
	      {
		dir[nGrp] = FORWARD;
		AddGroup( B, nGrp, index, t, FORWARD, Pos, M, &P );
		ps[nGrp] = index;
		found = TRUE;
	      }
	    else if( IP->RevComplement )
	      {
		sum += prob[index][REVERSE];
		if( sum >= dtot )
		  {
		    dir[nGrp] = REVERSE;
		    AddGroup( B, nGrp, index, t, REVERSE, Pos, M, &P );
		    ps[nGrp] = index;
		    found = TRUE;
		  }					    
	      }
	    index++;
	  }
	  
#ifdef _DEBUG_					  
	  if( found )
	    {
	      probTotal += prob[ps[nGrp]][dir[nGrp]];
	      probCnt++;
	      bf = BayesFactor( B,  dir, ps, t, Pos, M, &P,
				&p1Model, &p2Model);
	      bfTotal += bf;
	      p1Total += p1Model;
	      p2Total += p2Model;			  
	      PrintGroupInfo( B,  prob[ps[nGrp]][dir[nGrp]], bf, Pos, &P, p1Model, p2Model ); 
	    }
	  else	    
	    p_warning( "GroupSiteSampler: Could not find site." );
#endif	  
	}
      FREEP( prob, nGrpLen[nGrp] );
    }

  fprintf( B->IP->Datafiles->out_fpt, "Count = %d\n", probCnt );
  fprintf( B->IP->Datafiles->out_fpt, "max prob = %g\navg prob = %g\n", 
	  maxProb + IP->dnull_map, 
	  (probTotal / (double) probCnt) + IP->dnull_map );

  fprintf( B->IP->Datafiles->out_fpt, "BF at max = %g\navg BF = %g\nBF total = %g\n", 
	  log( maxbf ), log( bfTotal / (double) probCnt ), log( bfTotal ) );

  fprintf( B->IP->Datafiles->out_fpt, "p1Total = %g\np2Total = %g\n", 
	   log( p1Total ), log( p2Total ) );
  fprintf( B->IP->Datafiles->out_fpt, "p1Avg = %g\np2Avg = %g\n",
	   log( p1Total / (double) probCnt ), log( p2Total / (double) probCnt ) );
  fprintf( B->IP->Datafiles->out_fpt, "BF ratio  = %g\n", 
	   log( p1Total ) - log ( p2Total ) );
  
  
  free_prob( &P, IP);

  FREEP(Pos, IP->nNumMotifTypes);  
  free(Pos);
  free( ps );
  free( dir );
  free( nGrpLen );
}


void GroupSample( Model B )
{
  RPType      RP;
  IPtype      IP;
  int         nGrp;
  int         first;
  int         last;
  int         t;
  int         nGrp2 = 0;
  int         first2;
  int         last2;
  double      bf;
  double      bfTotal = 0.0;
  double      prob;
  double      probTotal = 0.0;
  int         probCnt = 0;
  PoSition    *Pos_t;
  PoSition    **Pos;
  Mlist       M; 
  int         **nNumMotifs;
  int         maxnum;
  int         i;
  int         seq;
  int         n;
  int         seq2;
  int         n2;
  ProbStruct  P;
  double      maxProb = -DBL_MAX;
  double      maxbf = -DBL_MAX;  
  int         *ps;
  int         *dir;
  int         dr;
  int         dr2;
  double      p1Model;
  double      p2Model;
  double      p1Total = 0.0;
  double      p2Total = 0.0;
  MaxResults  maxData;
  double      dLocMax = -DBL_MAX; 
  int         last_increase = 0;
  int         iter = 0;
  double      **dProbArray;
  double      kldist = 0;
  double      kldist12;
  double      kldist21;

  RP = B->RP;
  IP = B->IP;

  init_maxdata(&maxData);  /* BT 11/29/99 */

  NEWP( Pos, IP->nNumMotifTypes, PoSition);
  for(t=0; t < IP->nNumMotifTypes;t++)
    NEW(Pos[t], B->IP->nSeqLen, PoSition);

  NEW( ps, RP->nGroupCnt, int );
  NEW( dir, RP->nGroupCnt, int );

  copy_counts(B);  
  nNumMotifs = copy_motif_num(IP);
  zero_motifs(IP); /* set motif count to zero */

   /* calculate the probability of NULL model */
  IP->dnull_map = CalcMapProb(B, B->IP->is_defined[cl_R]);
  reset_motif_num(nNumMotifs,IP);
  
   /* mark all possible motif sites */
  set_indicator_vector(Pos, B->Seq, B);
  maxnum = findMaxNumMotif(IP);
  
  CountAlignments( B, Pos ); 
  SaveAlignmentCounts( B );
  
  for(t=0;t<IP->nNumMotifTypes;t++)
    {
      Pos_t=Pos[t];
      for(i = 0; i < IP->nSeqLen; i++) 
	{
	  Pos_t[i].nInMotif=FALSE;
	  Pos_t[i].nMotifStartPos=FALSE;
	}
    }

  M = initializeMotifs( IP );
  init_prob( &P, IP ); 

  InitGroups( B, Pos, M, &P );
  /*  set_posterior_prob(B->IP, B->C); */

  for( nGrp = 0; nGrp < RP->nGroupCnt; nGrp++ )
    {
      RemoveGroup( B, nGrp, Pos, M, &P );
      seq = FirstGroupSeq( B, nGrp );      
      first = SequenceStartPos( B, seq );
    	      
      for (t = 0; t < IP->nNumMotifTypes; t++) 
	{
	  last = SequenceEndPos( B, seq ) - IP->nMotifLen[t] + 1;
	  
	  for( n = first; n <= last; n++ )
	    {
	      iter++;
	      if(PossGroupFragStartPos(Pos, nGrp, n - first, t, B) )
		{
		  for( dr = FORWARD; 
		       dr <= REVERSE * ( IP->RevComplement ? 1 : 0); 
		       dr++ )
		    {
		      AddGroup( B, nGrp, n - first, t, dr, Pos, M, &P );
		      dir[nGrp] = dr;
		      ps[nGrp] = n - first;
		      if( RP->nGroupCnt == 1 )
			{
			  prob = CalcMapProb(B, B->IP->is_defined[cl_R]);
			  probTotal += prob;
			  probCnt++;
			  bf = BayesFactor( B,  dir, ps, t, Pos, M, &P,
					    &p1Model, &p2Model);
			  bfTotal += bf;
			  p1Total += p1Model;
			  p2Total += p2Model;			  
#ifdef _DEBUG_					  
			  PrintGroupInfo( B, prob, bf, Pos, &P, 
					  p1Model, p2Model ); 
#endif
			  if( prob > maxProb )
			    {
			      kldist = KLDistance( B,  dir, ps, t, 
						    Pos, M, &P,
						    nGrp, nGrp2,
						    &kldist12, &kldist21 );
			      maxProb = prob;
			      maxbf = bf;    /* BT 11/29/99 */
			      dProbArray = setElementProb(B, Pos, P);
			      maxData = setMaxData(B->F, B->IP, 
						   prob, iter,
						   Pos, 
						   &last_increase,
						   &maxProb, 
						   &dLocMax, 
						   maxData, dProbArray, 
						   B);
			      FREEP(dProbArray, findMaxNumMotif(B->IP));
			    }			  
			}
		      else
			{
			  for( nGrp2 = 0; nGrp2 < RP->nGroupCnt; nGrp2++ )
			    {
			      if( nGrp2 != nGrp )
				{			
				  RemoveGroup( B, nGrp2, Pos, M, &P );
				  seq2 = FirstGroupSeq( B, nGrp2 );      
				  first2 = SequenceStartPos( B, seq2 );
				  
				  last2 = SequenceEndPos( B, seq2 ) - IP->nMotifLen[t] + 1;
				  
				  for( n2 = first2; n2 <= last2; n2++ )
				    {
				      if(PossGroupFragStartPos(Pos, nGrp2, n2 - first2, t, B))
					{
					  for( dr2 = FORWARD; 
					       dr2 <= REVERSE * ( IP->RevComplement ? 1 : 0); 
					       dr2++ )
					    {
					      AddGroup( B, nGrp2, n2 - first2, t, dr2,
							Pos, M, &P );
					      
					      prob = CalcMapProb(B, B->IP->is_defined[cl_R]);
					      probTotal += prob;
					      probCnt++;
					      
					      dir[nGrp2] = dr2;
					      ps[nGrp2] = n2 - first2;
					      bf = BayesFactor( B,  dir, ps, t, Pos, M, &P,
								&p1Model, &p2Model);
					      bfTotal += bf;
					      p1Total += p1Model;
					      p2Total += p2Model;
#ifdef _DEBUG_					  
					      PrintGroupInfo( B, prob, bf, Pos, &P, 
							      p1Model, p2Model ); 
#endif
					      if( prob > maxProb )
						{
						  kldist = KLDistance( B,  dir, ps, t, 
								     Pos, M, &P,
								     nGrp, nGrp2,
								     &kldist12, &kldist21 );
						  maxProb = prob;
						  maxbf = bf;    /* BT 11/29/99 */
						  dProbArray = setElementProb(B, Pos, P);
						  maxData = setMaxData(B->F, B->IP, 
								       prob, iter,
								       Pos, 
								       &last_increase,
								       &maxProb, 
								       &dLocMax, 
								       maxData, dProbArray, 
								       B);
						  FREEP(dProbArray, findMaxNumMotif(B->IP));
						}
					      RemoveGroup( B, nGrp2, Pos, M, &P );
					    }
					}
				    }
				}
			    }		
			}
		      RemoveGroup( B, nGrp, Pos, M, &P );
		    }
		}
	    }
	}
    }

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

  print_info( B, maxData, TRUE, OTHER );  /* BT 11/29/99 */

  fprintf( B->IP->Datafiles->out_fpt, "Count = %d\n", probCnt );
  fprintf( B->IP->Datafiles->out_fpt, "max prob = %g\navg prob = %g\n", 
	  maxProb + IP->dnull_map, 
	  (probTotal / (double) probCnt) + IP->dnull_map );

  fprintf( B->IP->Datafiles->out_fpt, "BF at max = %g\navg BF = %g\nBF total = %g\n", 
	  log( maxbf ), log( bfTotal / (double) probCnt ), log( bfTotal ) );

  fprintf( B->IP->Datafiles->out_fpt, "KLDist at max = %g klDist12 = %g klDist21 = %g\n", 
	   kldist, kldist12, kldist21 );

  fprintf( B->IP->Datafiles->out_fpt, "p1Total = %g\np2Total = %g\n", 
	   log( p1Total ), log( p2Total ) );
  fprintf( B->IP->Datafiles->out_fpt, "p1Avg = %g\np2Avg = %g\n",
	   log( p1Total / (double) probCnt ), log( p2Total / (double) probCnt ) );
  fprintf( B->IP->Datafiles->out_fpt, "BF ratio  = %g\n", 
	   log( p1Total ) - log ( p2Total ) );
  
  
  free_prob( &P, IP);

  FREEP(Pos, IP->nNumMotifTypes);  
  free(Pos);
  free( ps );
  free( dir );
}


void RemoveGroup( Model B, int grp, PoSition **Pos, Mlist M, ProbStruct *P )
{
  int     seq;
  int     first;
  int     last;
  int     t;
  int     RevComp;
  int     n;
  
  for(seq = 0; seq < B->IP->nNumSequences; seq++) 
    {
      if( B->RP->groups[seq] == grp )
	{
	  first = SequenceStartPos( B, seq );
	  
	  for (t = 0; t < B->IP->nNumMotifTypes; t++) 
	    {
	      last = SequenceEndPos( B, seq ) - B->IP->nMotifLen[t] + 1;
	      
	      for( n = first; n <= last; n++ )
		{
		  if( Pos[t][n].nMotifStartPos )
		    {
		      RevComp = Pos[t][n].RevComp;        
		      adjust_counts(B,DELETE,n,t,RevComp);
		      B->IP->nNumMotifs[t][RevComp]--;
		      not_in_motif(Pos,n,B,t);    
		      delete_motif(B, n, M, t);
		      update_prob(P, B, (! B->IP->is_defined[cl_b]));
		      break;
		    }
		}
	    }
	}
    }
}


void AddGroup( Model B, int grp, int n, int t, int type, PoSition **Pos, Mlist M, 
	       ProbStruct *P )
{
  int  seq;
  int  index;

  for(seq = 0; seq < B->IP->nNumSequences; seq++) 
    {
      if( B->RP->groups[seq] == grp )
	{	  
	  index = n + SequenceStartPos( B, seq );

	  adjust_counts(B, ADD,index,t, type);
	  B->IP->nNumMotifs[t][type]++;
	  add_motif(B->IP, B->Seq, index, M, t, type);
	  set_in_motif(Pos, index, B, t, type);
	  update_prob( P, B, (! B->IP->is_defined[cl_b]));
	}
    }
}


int FirstGroupSeq( Model B, int grp )
{
  int     seq;
  
  for(seq = 0; seq < B->IP->nNumSequences; seq++) 
    {
      if( B->RP->groups[seq] == grp )
	{
	  return( seq );
	}
    }

  return(-1 );
}


int MinGroupSeqLength( Model B, int grp )
{
  int     seq;
  int     nLen;
  int     minLength = INT_MAX;
  
  for(seq = 0; seq < B->IP->nNumSequences; seq++) 
    {
      if( B->RP->groups[seq] == grp )
	{
	  nLen = SequenceLength( B, seq );	  
	  if( minLength > nLen )
	    minLength = nLen;
	}
    }

  return( minLength );
}


void InitGroups( Model B, PoSition **Pos, Mlist M, ProbStruct *P )
{
  int    seq;
  int    nGrp;
  int    t;
  int    first;
  int    last;
  int    n;
  short  bDone;

    	      
  for (t = 0; t < B->IP->nNumMotifTypes; t++) 
    {
      B->IP->nNumMotifs[t][FORWARD] = 0;
      B->IP->nNumMotifs[t][REVERSE] = 0;
    }

  for( nGrp = 0; nGrp < B->RP->nGroupCnt; nGrp++ )
    {
      bDone = FALSE;
      seq = FirstGroupSeq( B, nGrp );      
      first = SequenceStartPos( B, seq );
    	      
      for (t = 0; t < B->IP->nNumMotifTypes && ! bDone; t++) 
	{
	  last = SequenceEndPos( B, seq ) - B->IP->nMotifLen[t] + 1;
	  
	  for( n = first; n <= last && ! bDone; n++ )
	    {
	      if(PossGroupFragStartPos(Pos, nGrp, n - first, t, B) )
		{
		  AddGroup( B, nGrp, n - first, t, FORWARD, Pos, M, P );
		  bDone = TRUE;
		}
	    }
	}
    }  
}


int PossGroupFragStartPos(PoSition **Pos, int grp, int n, int t, Model B )
{
  int  seq;
  int  first;

  for(seq = 0; seq < B->IP->nNumSequences; seq++) 
    {
      if( B->RP->groups[seq] == grp )
	{
	  first = SequenceStartPos( B, seq );

	  if( ! PossFragStartPos(Pos,  n + first, t, B) )
	    return FALSE;
	}
    }
  return TRUE;
}


double BayesFactor( Model B,  int *dir, int *n, int t, 
		    PoSition **Pos, Mlist M, ProbStruct *P,
		    double *p1Model, double *p2Model )
{
  int    i;
  int    j;
  double prob1Model;
  double prob2Model;
  
  prob1Model = CalcMotifMap( B, t,  B->IP->is_defined[cl_R] );    
  
  prob2Model = 0.0;
  for( i = 0; i < B->RP->nGroupCnt; i++ )
    {
      for( j = 0; j < B->RP->nGroupCnt; j++ )
	{	
	  if( i != j )
	      RemoveGroup( B, j, Pos, M, P );
	}

      prob2Model += CalcMotifMap( B, t,  B->IP->is_defined[cl_R] );

      for( j = 0; j < B->RP->nGroupCnt; j++ )
	{	
	  if( i != j )
	    AddGroup( B, j, n[j], t, dir[j], Pos, M, P );
	}
    }
  
  *p1Model = exp( prob1Model );
  *p2Model = exp( prob2Model );

  return exp( prob1Model - prob2Model );
}


void PrintGroupInfo( Model B, double prob, double bf, PoSition **Pos, ProbStruct *P,
		     double p1Model, double p2Model )
{
  int  t;

  for( t = 0; t < B->IP->nNumMotifTypes; t++ )
    {
      fprintf( stdout, 
	       "-----------------------------------------------------------\n" );
      fprintf(stdout, "                          MOTIF %c\n\n", (char)(97 + t));
      DumpMotifPositions( t, B, Pos, stdout );
    }
  fprintf( stdout, "prob = %g Bayes Factor = %g\n", prob, bf );
  fprintf( stdout, "p1Model = %g p2Model = %g\n", p1Model, p2Model );
#ifdef _DEBUG_
  if( ! finite( prob ) || ! finite( bf ) || 
      ! finite( p1Model ) || ! finite( p2Model ) )
    {
      prob = prob;
    }
#endif
}


double KLDistance( Model B,  int *dir, int *n, int t, 
		   PoSition **Pos, Mlist M, ProbStruct *P,
		   int grp1, int grp2,
		   double *kldist12, double *kldist21 )
{
  int    i;
  int    j;
  IPtype IP;
  double **prob;

  IP = B->IP;

  NEWP( prob, IP->nMotifLen[t], double );
  for( i = 0; i < IP->nMotifLen[t]; i++ )
    {
      NEW( prob[i], IP->nAlphaLen, double );
    }

  for( i = 0; i < B->RP->nGroupCnt; i++ )
    {
      if( i != grp1 )
	RemoveGroup( B, i, Pos, M, P );
    }

  for( i = 0; i < IP->nMotifLen[t]; i++ )
    {
      for(j = 0; j < IP->nAlphaLen; j++) 
	{
	  prob[i][j] = P->dvInMotifProb[t][i][j];
	}
    }

  RemoveGroup( B, grp1, Pos, M, P );
  AddGroup( B, grp2, n[grp2], t, dir[grp2], Pos, M, P );
  
  *kldist12 = 0.0;
  *kldist21 = 0.0;
  for( i = 0; i < IP->nMotifLen[t]; i++ )
    {
      for(j = 0; j < IP->nAlphaLen; j++) 
	{
	  *kldist12 += prob[i][j] * LOG2( prob[i][j] / P->dvInMotifProb[t][i][j] );
	  *kldist21 += P->dvInMotifProb[t][i][j] * LOG2( P->dvInMotifProb[t][i][j] / prob[i][j] );
	}
    }

  for( i = 0; i < B->RP->nGroupCnt; i++ )
    {
      if( i != grp2 )
	AddGroup( B, i, n[i], t, dir[i], Pos, M, P );
    }

  FREEP( prob, IP->nMotifLen[t] );

  return( (*kldist12 + *kldist21) / 2.0 );
}
     

