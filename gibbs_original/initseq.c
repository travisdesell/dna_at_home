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
/* $Id: initseq.c,v 1.13 2008/09/16 19:50:20 Bill Exp $           */
/*                                                                        */
/* Author: Eric C. Rouchka, 1996.7                                        */
/*         Jun Zhu, 1996.8                                                */
/*         Bill Thompson 1/24/97                                          */
/* Description :  This file contains the functions that initialize the    */
/*                sequences.  They include:                               */
/*                    set_random_sequences                                */
/*                    set_indicator_vector                                */
/*                    get_strings                                         */ 
/**************************************************************************/

#include "initseq.h"

void CheckSequenceLengths( Model B );
void ResetMotifCounts( Model B, PoSition **Pos, int **startPos );


/************************  set_random_sequences   ********************/
/*                                                                   */
/*                                                                   */
/* DESCRIPTION : This function randomly sets the initial alignment of*/
/*               motif elements based on the user inputs.  If an     */
/*               initial alignment cannot be found, an error message */
/*               is produced and execution is terminated             */
/*********************************************************************/

void set_random_sequences(Model B, PoSition **Pos,
                          int **startPos)
{
   int      i, j, t;
   int      loc, valid;
   long     it_num;
   char     *error_msg;
   int      nNumMotifs;
   PoSition *Pos_t;       /* Pos_t=Pos[t] */
   int      nSeq;
   int      nOffset;
   int      nPos;
   int      nMotifLen;
   int      k;

   NEW(error_msg, 100, char); 
   strcpy( error_msg, 
	   "cannot allocate random motifs\n shorten the number of motifs or motif length \0") ;	/* BT 1/24/97 */
 
   for( i = 0; i < B->IP->nNumSequences; i++ )
     B->RP->nSites[i] = 0;

   for(t = 0; t < B->IP->nNumMotifTypes; t++) 
     {
       nMotifLen = MotifWidth( B, t );
       Pos_t=Pos[t];
       nNumMotifs = NUMMOTIFS(B->IP->nNumMotifs[t]);       

       i = 0;
       while( i < nNumMotifs )
	 {
	   valid = FALSE;
	   it_num = 0;
	   while(!valid) 
	     {
	       /* randomly select a location and make sure it's a valid */
	       /* position.                                             */ 
	       if(B->IP->site_samp )
		 {
		   if( i >= B->IP->nNumSequences )
		     p_internal_error( "Site Sampler: Number of random sites exceeds number of sequences");
		   loc = RandomInterval( SequenceStartPos( B, i ),    /* BT 12/7/2000 */
					 SequenceEndPos( B, i ) - nMotifLen );
		 }
	       else
		 {
		   nSeq = (int) (drand() * (double) B->IP->nNumSequences);
		   loc = RandomInterval( SequenceStartPos( B, nSeq ),   
					 SequenceEndPos( B, nSeq ) - nMotifLen );
		 }
	       
	       valid = TRUE;  
	       if( loc < 0 )
		 valid = FALSE;  
	       else if(!Pos_t[loc].nPossStartPos)
		 valid = FALSE;
	       else if( ! PossMotifStartPos( Pos, loc, t, B) )
		 valid = FALSE;  /* BT 10/29/2001 */
	       else if( loc > (*B->Seq->nvEndLocs)[Pos_t[loc].nSeq] - nMotifLen )   /* BT 4/25/97 */
		 valid = FALSE;
	       else if( B->IP->is_defined[cl_E] && 
			B->RP->nSites[Pos_t[loc].nSeq] >= B->RP->nMaxBlocks )
		 valid = FALSE;
	       else 
		 {
		   if( B->IP->is_defined[cl_T] )
		     {
		       if( B->IP->is_defined[cl_E] )  /* BT 11/09/99 */
			 {
			   nPos = loc - SequenceStartPos( B, Pos_t[loc].nSeq );
			   valid = (B->RP->dProbCons[Pos_t[loc].nSeq][nPos] > drand());
			 }
		       else
			 valid = InFootprintRegion( B, loc, t );   /* BT 11/09/99 */
		     }
		   if( valid && IsPhyloSeq( B, Pos_t[loc].nSeq ) )
		     {
		       nOffset = loc - SequenceStartPos( B, Pos_t[loc].nSeq );
		       nSeq = B->Phylo->phyloTreeStartSeq[B->Phylo->phyloIndex[Pos_t[loc].nSeq]];
		       loc = SequenceStartPos( B, nSeq ) + nOffset;
		       if( nSeq + SpeciesInc( B, nSeq ) > B->Phylo->maxPhyloSeq + 1 )
			 valid = FALSE;
		       else
			 {
			   nOffset = SequenceLength( B, nSeq ); 
			   for( k = 0; k < SpeciesInc( B, nSeq ); k++ )
			     if( ! PossMotifStartPos( Pos, loc + k * nOffset, t, B) )
			       valid = FALSE;  /* BT 10/29/2001 */			   
			 }
		     }
		 }
		   /* check if it overlaps with anybody */
	       if( valid )
		 valid = !overlaps(Pos, B->IP,loc, nMotifLen);
	       
	       it_num++;                           /* exit if a sufficient     */
	       if(it_num == 100000)
		 {
		   p_warning(error_msg);
		   free(error_msg);
		   ResetMotifCounts( B, Pos, startPos );
		   return;
		 }
	     }

	   Pos_t[loc].nMotifStartPos = TRUE;        /* mark the position in the   */
	   for(j = 0; j < nMotifLen; j++)           /* alignment         */
	     Pos_t[j + loc].nInMotif =TRUE;
	   if(i < B->IP->nNumMotifs[t][FORWARD])
	     Pos_t[loc].RevComp = FORWARD;
	   else
	     Pos_t[loc].RevComp = REVERSE;
	   nSeq = SequenceFromPosition( B, loc );
	   Pos_t[loc].nSeq = nSeq;
	   B->RP->nSites[Pos_t[loc].nSeq]++;
	   if( IsPhyloSeq( B, nSeq ) )
	     {
	       nOffset = SequenceLength( B, nSeq );
	       for( k = 1; k < SpeciesInc( B, nSeq ); k++ )
		 {		   
		   Pos_t[loc + nOffset*k].nMotifStartPos = TRUE;       
		   for(j = 0; j < nMotifLen; j++)
		     Pos_t[j + loc + nOffset*k].nInMotif = TRUE;
		   Pos_t[loc + nOffset*k].RevComp = Pos_t[loc].RevComp;
		   B->RP->nSites[Pos_t[loc + nOffset*k].nSeq]++;
		 }	       
	     }
	   i += SpeciesInc( B, nSeq );
	 }
     }
   free(error_msg);

   ResetMotifCounts( B, Pos, startPos );
}


/*********************  short overlaps *******************************/
/*                                                                   */
/*                                                                   */
/* DESCRIPTION : This function checks to see if the current position */
/*               overlaps any existing motif elements                */
/*********************************************************************/

short overlaps(PoSition **Pos,IPtype IP, int loc, int len)

{
   int j,t;
   PoSition *Pos_t; /* Pos_t=Pos[t] */

   for(t=0;t<IP->nNumMotifTypes;t++){
     Pos_t=Pos[t];
     for(j = 0; j < len; j++) {
       if(Pos_t[j + loc].nInMotif == TRUE)
	 return TRUE;
     }
   }
   return FALSE;
}


/**********************  set_indicator_vector   *************************/
/*                                                                      */
/*                                                                      */
/* DESCRIPTION : This function sets the value of the indicator vector to*/
/*               a zero (can start a new motif) or -1 (cannot --        */
/*               overlaps end of sequences.)                            */
/************************************************************************/

void set_indicator_vector(PoSition **Pos,
                          Stringstruct *Seq, Model B)
{
  int      i, j, t;
  int      pos = 0;
  PoSition *Pos_t;  /* Pos_t=Pos[t] */
  int      nSeq;		/* BT 4/21/97 */
  IPtype   IP;
  int      len;
  int      start;
  int      n;
  int      p;
  
  /* Pos{MotifType][SeqLen] */
  /* Mark all possible start position of a motif */
  /* at this point no motif segments are selected */
  
  IP = B->IP;
  
  for(t=0;t<IP->nNumMotifTypes;t++)
    {
      Pos_t=Pos[t];
      nSeq = 0;				/* BT 4/21/97 */
      for(i = 0; i < IP->nSeqLen; i++) 
	{
	  if( i >= (*Seq->nvEndLocs)[nSeq] )	/* BT 4/21/97 */
	    nSeq++;
	  Pos_t[i].nSeq = nSeq; 
	  if( ((*Seq->ProcessedSTR)[i] == 'X') || 
	      ((*Seq->ProcessedSTR)[i] == 'x') ||
	      ((*Seq->ProcessedSTR)[i] == 'u') || 
	      (IP->nAlphaLen == 4 &&
	       (((*Seq->ProcessedSTR)[i] == 'N') || 
		((*Seq->ProcessedSTR)[i] == 'n')))) /* BT 4/9/97 */
	    {
	      
	      for( j = 0; j < IP->nMotifLen[t]; j++ )		/* BT 4/21/97 */
		{
		  if( i - j >= SequenceStartPos( B, nSeq ) )
		    Pos_t[i - j].nPossStartPos = FALSE;
		}
	    }   
	  else
	    Pos_t[i].nPossStartPos = TRUE;

	  Pos_t[i].nMotifStartPos = FALSE;
	  Pos_t[i].nInMotif = FALSE;
	  Pos_t[i].bEndOfSeq = FALSE;		/* BT 4/18/97 */
	}
   }

   /* check overlap at the end */
   for(i = 0; i < IP->nNumSequences; i++) 
     {
       pos = (*Seq->nvEndLocs)[i];
       for(t=0;t<IP->nNumMotifTypes;t++)
	 {
	   Pos_t=Pos[t];
	   Pos_t[pos - 1].bEndOfSeq = TRUE;		/* BT 4/18/97 */ 
	   for(j = 1; j < (int) ((1.0 - IP->glOverlapParam) *IP->nMotifLen[t]); j++)   /* BT 3/24/97 */
	     {
	       if((pos - j) >= 0) 
		 Pos_t[pos - j].nPossStartPos =FALSE;
	     }
	 }
     }
   
   for(i = 0; i < IP->nNumSequences; i += SeqInc( B, i )) 
     {
       if( IsPhyloSeq( B, i ) )
	 {
	   for( t = 0; t < IP->nNumMotifTypes; t++)
	     {
	       Pos_t=Pos[t];
	       len = SequenceLength( B, i );
	       start = SequenceStartPos( B, i );
	       for( n = 0; n < len; n++ )
		 {
		   if( ! PossibleStartPos( Pos, start + n, t, len, i, B) )
		     {
		       for( p = 0; p < SeqInc( B, i ); p++ )
			 Pos_t[start + n + p *len].nPossStartPos = FALSE;
		     }
		 }
	     }
	 }
     }
}


void get_strings(Model B)

/**************************************************************************/
/* FUNCTION NAME : get_strings                                            */
/*                                                                        */
/* DESCRIPTION : This function reads in the sequences from the designated */
/*               input file.  The sequences are to be in fast-A format,   */
/*               meaning that there is a one line description for each    */
/*               sequence beginning with the greater than symbol (>).  It */
/*               is assumed that the elements can be stored in either     */
/*               upper case or lower case, and there can be white space   */
/*               between them.  The sequences are concatenated together   */
/*               into one string (R) and a vector is used to keep track of*/
/*               the ends of the individual sequences (nvEndLocs).        */ 
/**************************************************************************/
{
   int          i, strsize, seqsize,loc,pos, t;
   char         tmpch,*tmpmsg;
   int          nNumMotifs;
   IPtype       IP;
   Stringstruct *Seq;
   char         orgCh;
   int          headerCnt;

   IP = B->IP;
   Seq = B->Seq;

   /* initialize the buffer assuming 20 sequences, each has 50 elements */
   NEW(tmpmsg, 80, char);
   NEW((*Seq->nvEndLocs), 20, int);
   NEWP( IP->fastaHeader, 20, char );
   IP->nNumSequences = 0;
   IP->nSeqAllocCnt = 0;    /* BT 8/27/97 */
   loc = 0;
   pos=0;
   strsize = 100;
   seqsize=50;
   NEW((*Seq->R), strsize, char);
   NEW(Seq->SeqLen, 20, int);
   NEWP(Seq->Orig,20,short);
   for(i=0;i< 20;i++)
   {
       NEW(Seq->Orig[i], seqsize, short);
   }
   IP->nSeqAllocCnt += 20;    /* BT 8/27/97 */

   while((tmpch=fgetc(IP->Datafiles->fpt)) != (char) EOF)  /* BT 12/14/99 */
     {
      if(tmpch == '>'){
         /** Once we see the '>', we know that we have reached  **/
         /** another comment.  We need to mark the end location **/
         /** of the substring that came before the comment,     **/
         /** increase the number of strings, and strip out the  **/
         /** comment                                            **/

         if(loc > 0)
            (*Seq->nvEndLocs)[IP->nNumSequences - 1] = loc;

	 /* save the sequence length */
	 if(pos>0)
	   Seq->SeqLen[IP->nNumSequences-1]=pos;;

         /* reset the counters */
         IP->nNumSequences++;
	 pos=0;
	 seqsize=50;

         /* This if loop will reset the size of the End Location */
         /* Array if it is no longer large enough                */
	 if(IP->nNumSequences % 20 == 0) /* BT 12/31/97 */
	 {
            (*Seq->nvEndLocs)=(int *)realloc((*Seq->nvEndLocs), 
			      (IP->nNumSequences + 20)*sizeof(int));
            Seq->SeqLen=(int *)realloc(Seq->SeqLen, 
			      (IP->nNumSequences + 20)*sizeof(int));
            Seq->Orig=(short **)realloc(Seq->Orig, 
			      (IP->nNumSequences + 20)*sizeof(short *));
            IP->fastaHeader = (char **) realloc( IP->fastaHeader, 
			      (IP->nNumSequences + 20)*sizeof(char *));
	    IP->nSeqAllocCnt += 20;    
            for (i=IP->nNumSequences;i< IP->nNumSequences+20;i++)
	      NEW(Seq->Orig[i],seqsize,short);
         }

         /* Since this is the description line, keep on reading from */
         /* the file until we reach the end of line                  */
	 NEW( IP->fastaHeader[IP->nNumSequences-1], HEADER_ALLOC_SIZE, char );
	 headerCnt = 0;
         tmpch=fgetc(IP->Datafiles->fpt);
	 while( tmpch != '\n' && tmpch != (char) EOF )
	   {
	     if( (headerCnt > 0) && (headerCnt % HEADER_ALLOC_SIZE == 0) )
	       IP->fastaHeader[IP->nNumSequences-1] = 
		 realloc( IP->fastaHeader[IP->nNumSequences-1], 
			  (headerCnt + HEADER_ALLOC_SIZE) * sizeof( char ) ); 
	     IP->fastaHeader[IP->nNumSequences-1][headerCnt] = tmpch;
	     headerCnt++;
	     tmpch=fgetc(IP->Datafiles->fpt);
	   }
	 
	 if( tmpch == (char) EOF )
	   {
	     sprintf( tmpmsg, "Missing sequence data - sequence %d\n", 
		      IP->nNumSequences );
	     p_error( tmpmsg );
	   }

	 if( (headerCnt > 0) && (headerCnt % HEADER_ALLOC_SIZE == 0) )
	   IP->fastaHeader[IP->nNumSequences-1] = 
	     realloc( IP->fastaHeader[IP->nNumSequences-1], 
		      (headerCnt + HEADER_ALLOC_SIZE) * sizeof( char ) ); 
	 IP->fastaHeader[IP->nNumSequences-1][headerCnt] = '\0';
      }

      /** The next if statement allows us to ignore spaces, newlines, **/
      /** and tabs.  If we want to ignore anything else, add it here  **/
      orgCh = tmpch;
      if((tmpch != ' ') && (tmpch != '\n') && (tmpch != '\t')) {
	    if((int)tmpch < 97)                    /* Convert to lowercase */
	      tmpch = (char)((int)tmpch + 32);     /* letters              */
            if(IP->nAlphaLen == 4) {
               switch(tmpch) {
                  case 'a': break;
                  case 'c': break;
                  case 'n': break;               /* allow low complexity */
                  case 'x':                      /* regions              */
                     break;
                  case 't':
                     tmpch = 'b'; break;
                  case 'g':
                     tmpch = 'd'; break;
                  default:
                     sprintf(tmpmsg, "%c (ascii %d) is not a valid Nucleic Acid in sequence %d\n", 
			     orgCh, (int) ((unsigned char) orgCh + 0), IP->nNumSequences );
                     p_error(tmpmsg);
                     break;
               }
            }
            else {
               switch(tmpch) {              /* this will allow us to convert */
                  case 'v' :                /* to add into separate buckets  */
                     tmpch = 'b'; break;
                  case 'w' :
                     tmpch = 'j'; break;
                  case 'y' :
                     tmpch = 'o';  break;
	          case 'x':/* accept X as input */
		     tmpch = 'u';  break;
                  default :
                     break;
               }

               if(((int)(tmpch) < 97) || ((int)(tmpch) > 117)) {
                  sprintf(tmpmsg, "%c (ascii %d)  is not a valid amino acid abbreviation in sequence %d\n",
			  orgCh, (int) ((unsigned char) orgCh + 0), IP->nNumSequences );
                  p_error(tmpmsg);
               }

            }
         (*Seq->R)[loc] = tmpch;
	 if( IP->nNumSequences == 0 )
	   p_error( "Invalid data - probably a missing FASTA header" );
         Seq->Orig[IP->nNumSequences-1][pos] = (int)tmpch-97;
	 loc++;
	 pos++;

         /** if the string is going to be longer than the space allocated**/
         /** for it, go ahead and allocate some more space               **/

         if((loc + 1) == strsize) 
	 {
            strsize+=100;
	    (*Seq->R)=(char *)realloc((*Seq->R),(strsize)*sizeof(char));
	 }
         if((pos+2)==seqsize)
	 {
	     seqsize+=50;
	    Seq->Orig[IP->nNumSequences-1]=(short *)realloc(
					 Seq->Orig[IP->nNumSequences-1],
					 (seqsize)*sizeof(short));
         }
      }
   }
   (*Seq->nvEndLocs)[IP->nNumSequences - 1] = loc;  
   Seq->SeqLen[IP->nNumSequences - 1] = pos;

   /* Need to add in the last end*/
   /* & null terminate the string*/
   (*Seq->R)[loc] = '\0';                      
   IP->nSeqLen = strlen(*Seq->R);

   NEW(IP->nPossSites, IP->nNumMotifTypes, int);
   for(t = 0; t < IP->nNumMotifTypes; t++)
      IP->nPossSites[t] =IP->nSeqLen-(IP->nMotifLen[t]-1)*IP->nNumSequences;

   /* for site sampler the number of motif is number of sequences */
   if(IP->site_samp )     /* BT 9/14/98 */
     {
       nNumMotifs = IP->nNumSequences;

       for(t = 0; t < IP->nNumMotifTypes; t++) 
	 {
	   if(IP->RevComplement) 
	     {
	       IP->nNumMotifs[t][FORWARD] = nNumMotifs / 2 + nNumMotifs % 2;
	       IP->nNumMotifs[t][REVERSE] = nNumMotifs / 2;
	     }
	   else 
	     {
	       IP->nNumMotifs[t][FORWARD] = nNumMotifs;
	       IP->nNumMotifs[t][REVERSE] = 0;
	     }
	 }
     }

   free(tmpmsg);

   if( IP->is_defined[cl_D] && ! B->Phylo->maxPhyloSeq )
     B->Phylo->maxPhyloSeq = IP->nNumSequences - 1;

   CheckSequenceLengths( B );

   if( B->Phylo->maxPhyloSeq > IP->nNumSequences )
     p_error( "get_strings: The number of sequences given with -D exceeds the sequence count." );
}


void SetInitSequences(Model B, PoSition **Pos, int **startPos)   /* BT 8/5/98 */
{
  int    i;
  int    j;
  int    n;
  int    t;
  int    fcnt;
  int    rcnt;
  int    pos = 0;
  IPtype IP;
  int    nEnd;
  char   *msg;

  IP = B->IP;
  
  for( i = 0; i < B->IP->nNumSequences; i++ )
    B->RP->nSites[i] = 0;

  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      fcnt = 0;
      rcnt = 0;
      for( n = 0; n < IP->nSeqLen; n++ )
	{
	  if( B->InitPos[t][n].nMotifStartPos )
	    {
	      if( ! B->InitPos[t][n].RevComp )
		{
		  startPos[t][fcnt] = n;                     /* mark the starting */
		  fcnt++;
		}
	      else
		{
		  startPos[t][rcnt + IP->nNumMotifs[t][FORWARD]] = n;
		  rcnt++;
		}
	      Pos[t][n].nMotifStartPos = TRUE;        /* position in the   */
	      nEnd = MotifWidth( B, t );
	      for(j = 0; j < nEnd; j++) /* alignment         */
		Pos[t][n + j].nInMotif =TRUE;
	      Pos[t][n].RevComp = B->InitPos[t][n].RevComp;
	      Pos[t][n].nSeq = SequenceFromPosition( B, n );
	      B->RP->nSites[Pos[t][n].nSeq]++;
	      if( IP->is_defined[cl_E] && B->RP->nSites[Pos[t][n].nSeq] > B->RP->nMaxBlocks )
		{
		  NEW( msg, 1024, char);
		  sprintf( msg, "Can't initialize sequence %d  with more than %d sites",  
			   Pos[t][n].nSeq + 1, B->RP->nMaxBlocks );
		  p_error( msg );
		}
	    }
	}
      pos += fcnt + rcnt;
    }
}


void CheckSequenceLengths( Model B )
{
  int  seq;
  char msg[128];

  for( seq = 0; seq < B->IP->nNumSequences; seq++ )
    {
      if( SequenceLength( B, seq ) == 0 )
	{
	  sprintf( msg, "get_strings: Sequence %d has length = 0\n", seq+1 );
	  p_error( msg );
	}
    }
}


void ResetMotifCounts( Model B, PoSition **Pos, int **startPos )
{
  IPtype   IP;
  int      t;
  int      n;
  PoSition *pos_t;
  int      maxnum;
  int      cnt;

  IP = B->IP;

  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      pos_t = Pos[t];
      IP->nNumMotifs[t][FORWARD] = 0;
      IP->nNumMotifs[t][REVERSE] = 0;

      for( n = 0; n < IP->nSeqLen; n++ )
	{
	  if( pos_t[n].nMotifStartPos )
	    {
	      if( pos_t[n].RevComp )
		IP->nNumMotifs[t][REVERSE]++;
	      else
		IP->nNumMotifs[t][FORWARD]++;
	    }
	}
    }

  maxnum = findMaxNumMotif(IP);
  for( t = 0; t < IP->nNumMotifTypes; t++ )
    {
      pos_t = Pos[t];
      cnt = 0;
      startPos[t] = realloc( startPos[t], maxnum * sizeof( int ) );
      for( n = 0; n < IP->nSeqLen; n++ )
	{
	  if( pos_t[n].nMotifStartPos )
	    {
	      startPos[t][cnt] = n;
	      cnt++;
	    }
	}
    }
}
