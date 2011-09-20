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
/* $Id: fragadj.c,v 1.7 2007/07/20 14:19:45 Bill Exp $           */
/*                                                                        */
/* AUTHOR: Eric C. Rouchka     July 16, 1996                              */
/*         Jun Zhu   Oct. 10, 1996                                        */
/*         Bill Thompson 1/22/97                                          */
/*                         						  */
/**************************************************************************/

#include "fragadj.h"

/*------------------------------------------------------------------------*/
/*------------------------- LOCAL FUNCTION PROTOTYPES --------------------*/
/*------------------------------------------------------------------------*/



/*===============================================================*/
/* FUNCTION NAME : adjust_fragment                               */
/* DESCRIPTION : This function will shift the current alignment  */
/*               according to how many places the first position */
/*               of the motif shifts -- thus, it is not of shift */
/*               of all of the positions necessarily, just the   */
/*               first                                           */
/*===============================================================*/

void adjust_fragment(Model B, Mlist M, PoSition **Pos, ProbStruct *P)
{
  MlistEl      *curr;
  int          t, **old, **new;
  int          f1, f2, shift, revshift, l1, l2;
  int          nWidth;
  
  for(t = 0; t < B->IP->nNumMotifTypes; t++) 
    {
      old = B->F->nOldColMask;
      new = B->F->nColMask;
      f1 = FirstCol(B->F, t);
      l1 = LastCol(B->F, t);
      B->F->nColMask = old;
      f2 = FirstCol(B->F, t);
      l2 = LastCol(B->F, t);
      shift = f1 - f2;
      revshift = l1 - l2;
      curr = M[t]->Motifs;
      
      B->F->nColMask = new;			/* BT 2/7/97 */
      
      while(curr != NULL) 
	{
	  B->F->nColMask = old; 
	  nWidth = B->F->FragWidth[t];         /* BT 2/5/2001 */
	  B->F->FragWidth[t] = B->F->nOldFragWidth[t];
	  adjust_counts(B, DELETE, curr->pos,t, curr->RevComp);
	  not_in_motif(Pos,curr->pos,B,t);
	  B->F->shift[t] = shift;
	  if(!curr->RevComp) { curr->pos += shift; }
	  else               { curr->pos -= revshift; } 
	  B->F->nColMask = new;
	  B->F->FragWidth[t] = nWidth;
	  adjust_counts(B, ADD, curr->pos, t, curr->RevComp);
	  set_in_motif(Pos,curr->pos,B,t,curr->RevComp); 
	  curr = curr->next;
	}
      SetFragPos( B, t );
      
    }
  
  if( ! B->IP->is_defined[cl_sample_model] )
    update_prob(P, B, (! B->IP->is_defined[cl_b]));
}

int NextCol(Ftype F, int t, int last)

   /******************************************************************/
   /* FUNCTION NAME : NextCol                                        */
   /* DESCRIPTION   : Finds the Next Column that is turned on in the */
   /*                 Column mask vector used for fragmentation.     */
   /******************************************************************/
    
{
   register int n;
   int          maxLen;
   int          *mask;

   maxLen = F->nMaxLen[t];
   mask = F->nColMask[t];
   
   for(n = last + 1; n < maxLen; n++)
     if(mask[n] == COL_ON)
       return n;
   return -1; /* NO NEXT COLUMN */
}

int PrevCol(Ftype F, int t, int last)

   /*******************************************************************/
   /* FUNCTION NAME : PrevCol                                         */
   /* DESCRIPTION : Finds the previous column that is turned on in the*/
   /*               column mask vector used for fragmentation         */
   /*******************************************************************/

{
   register int n;
   int          *mask;

   mask = F->nColMask[t];
      
   for(n = last - 1; n >= 0; n--)
      if(mask[n] == COL_ON)
         return n;
   return -1; /* NO PREVIOUS COLUMN */
}
 
int FirstCol(Ftype F, int t)

   /*=================================================================*/
   /* FUNCTION NAME : FirstCol                                        */
   /* DESCRIPTION   : Finds the first Column in the fragmentation     */
   /*                 vector                                          */
   /*=================================================================*/
{
   int n;
   for( n = 0; n < F->nMaxLen[t]; n++)
      if(F->nColMask[t][n] == COL_ON) 
         return n;
   return -1; /* NO COLUMNS TURNED ON */
}

int LastCol(Ftype F, int t)

   /*==================================================================*/
   /* FUNCTION NAME : LastCol                                          */
   /* DESCRIPTION   : Finds the last column in the fragmentation vector*/
   /*==================================================================*/

{
   int n;
   for(n = F->nMaxLen[t] - 1; n >= 0; n--)
      if(F->nColMask[t][n] == COL_ON)
         return n;
   return -1; /* NO COLUMNS TURNED ON */
}

void printFragVec(Model B)

   /*=================================================================*/
   /* FUNCTION NAME : printFragVec                                    */
   /* DESCRIPTION   : Prints the contents of the fragmentation vector */
   /*=================================================================*/

{
   int t, n;
   
   printf("\n");
   for(t = 0; t < B->IP->nNumMotifTypes; t++) 
     {
       printf( "motif: %d ", t );
       for(n = 0; n < B->F->nMaxLen[t]; n++)
         printf("%3d", B->F->nColMask[t][n]);
       printf("\n");
     }
}

/********************   setFragVec  *********************************/
/*                                                                  */
/* DESCRIPTION   : Takes the alignment and sets the fragmenting     */
/*                 vector according to all of the possible spaces   */
/*                 given how many spaces to the left and right that */
/*                 the motif can move                               */
/*                 It is assumed that everything between the first  */
/*                 column and the last column is set properly       */
/*==================================================================*/


void setFragVec(Model B, Mlist M,  PoSition **Pos)
{
  MlistEl *curr; 
  int     j, t, diff;
  int     first, last, first_blk;
  int     max_left, max_right;
  int     currleft, currright, temp;
  int     cLeft, cRight, t2, i;
  
  for( t = 0; t < B->IP->nNumMotifTypes; t++ ) 
    {
      B->F->shift[t] = 0;
    }

  for(t = 0; t < B->IP->nNumMotifTypes; t++) 
    {
      /*      diff = B->F->FragWidth[t] - M[t]->nMotifLen; */
      diff = B->F->FragWidth[t] - 1; 	/* BT 4/2/97 */
      curr = M[t]->Motifs;
      
      if( curr != NULL) 
	{			/* BT 1/22/97 */
	  max_left = curr->left;
	  max_right = curr->right;	
	  
	  while(curr != NULL) 
	    { 
	      for(t2 = 0; t2 < B->IP->nNumMotifTypes; t2++)  /* BT 8/25/97 */
		{
		  for( cLeft = curr->left, i = curr->pos - 1; i >= curr->left; i--)   
		    {
		      if( Pos[t2][i].nInMotif )
			{
			  cLeft = i + 1;
			  i = cLeft;
			  break;
			}	     
		    }
		  /*		  for( cLeft = curr->left, i = curr->left; i < curr->pos; i++)   
		    {
		      if( Pos[t2][i].nMotifStartPos )
			{
			  cLeft = i + B->IP->nMotifLen[t2] - 1;
			  i = cLeft;
			}	     
		    } */
		  for( cRight = curr->right, 
			 i = curr->pos + B->IP->nMotifLen[t]; 
		       i < curr->right; 
		       i++)
		    {
		      if( Pos[t2][i].nMotifStartPos )
			{
			  cRight = i - 1;
			  break;
			}
		    }	       		  
		  
		  currleft = curr->pos - cLeft;
		  currright = cRight - curr->pos - diff;
		  
		  if(curr->RevComp) 
		    {
		      temp = currleft;
		      currleft = currright;
		      currright = temp;
		    } 
		  
		  if(currleft < max_left)			/* BT 1/24/97 */
		    {
		      max_left = currleft;                  /* BT 8/29/2000 */
		      max_left = max( 0, max_left );
		    }
		  if(currright < max_right)			/* BT 1/24/97 */
		    max_right = max( 0, currright);          /* BT 4/2/97 */
		}   
	      curr=curr->next;
	    }
         
	  first = FirstCol(B->F, t);
	  last = LastCol(B->F, t);
	  for(j = 0; j < first - max_left; j++)    /* Process those elements */
	    if(  B->F->nColMask[t][j] != COL_DONT_USE )
	      B->F->nColMask[t][j] = COL_BLOCKED;   /* to the left of motif   */
	  
	  for(j = max(first - max_left, 0); j < first; j++)  /* BT 1/24/97 */
	    if(  B->F->nColMask[t][j] != COL_DONT_USE )
	      B->F->nColMask[t][j] = COL_OFF;
	  
	  /* BT 3/31/97 */
	  first_blk = min( first + B->F->FragWidth[t] + max_right, B->F->nMaxLen[t]);  
	                                                       /*Process elements*/
	  for(j = last + 1; j < first_blk; j++)                /*to the right of */
	    if(  B->F->nColMask[t][j] != COL_DONT_USE )
	      B->F->nColMask[t][j] = COL_OFF;                    /*the motif       */
	  for(j = first_blk; j < B->F->nMaxLen[t]; j++)
	    if(  B->F->nColMask[t][j] != COL_DONT_USE )
	      B->F->nColMask[t][j] = COL_BLOCKED;   
	}
   }

  /* Make sure the new fragmentation does not cause overlap */
  /*
  for( t = 0; t < B->IP->nNumMotifTypes; t++ ) 
    {
      first = FirstCol(B->F, t);
      last = LastCol(B->F, t);
      max_left = 0;
      for( i = first - 1; i >= 0; i-- )
	{
	  if( B->F->nColMask[t][i] == COL_BLOCKED )
	    {
	      break;
	    }
	  else
	    {
	      max_left++;
	    }
	}

      max_right = 0;
      for( i = last + 1; i < B->F->nMaxLen[t]; i++ )
	{
	  if( B->F->nColMask[t][i] == COL_BLOCKED )
	    {
	      break;
	    }
	  else
	    {
	      max_right++;
	    }
	}

      curr = M[t]->Motifs;      
      while(curr != NULL) 
	{
	  minPos = curr->pos - max_left;
	  maxPos = curr->pos + B->F->FragWidth[t] + max_right - 1;
	  
	  for(t2 = 0; t2 < B->IP->nNumMotifTypes; t2++) 
	    {
	      first2 = FirstCol(B->F, t2);
	      last2 = LastCol(B->F, t2);
	      max_left2 = 0;
	      for( i = first2 - 1; i >= 0; i-- )
		{
		  if( B->F->nColMask[t2][i] == COL_BLOCKED )
		    {
		      break;
		    }
		  else
		    {
		      max_left2++;
		    }
		}
	      
	      max_right2 = 0;
	      for( i = last2 + 1; i <B->F->nMaxLen[t]; i++ )
		{
		  if( B->F->nColMask[t2][i] == COL_BLOCKED )
		    {
		      break;
		    }
		  else
		    {
		      max_right2++;
		    }
		}
	      
	      curr2 = M[t2]->Motifs;      
	      while(curr2 != NULL) 
		{
		  if( curr != curr2 )
		    {
		      maxPos2 =  curr2->pos + B->F->FragWidth[t2] + max_right2 - 1;
		      minPos2 = curr2->pos - max_left2;
		      if(  curr->pos >= curr2->pos && minPos <= maxPos2 && max_left > 0 )
			{
			  overlap = maxPos2 - minPos + 1;
			  i = first - max_left; 
			  for( j = 0; j < overlap; j++ )
			    {
			      if( B->F->nColMask[t][i+j] != COL_ON )
				B->F->nColMask[t][i+j] = COL_BLOCKED;
			    }
			}
		      if(  maxPos >=minPos2 && curr->pos <= curr2->pos && max_right > 0 )
			{
			  overlap = maxPos - minPos2 + 1;
			  i = last + max_right;
			  for( j = 0; j < overlap; j++ )
			    {
			      if( B->F->nColMask[t][i-j] != COL_ON )
				B->F->nColMask[t][i-j] = COL_BLOCKED;
			    }
			}
		    }
		  curr2 = curr2->next;
		}
	    }	  	  

	  curr = curr->next;
	}
      
    }
  */  

#ifdef _DEBUG_
   printFragVec(B);
#endif   
}

void setFragCnts(Model B, Mlist M)

   /*=======================================================================*/
   /* FUNCTION NAME : setFragCnts                                           */
   /* DESCRIPTION   : Sets the observed residue counts for the posible      */
   /*                 motif columns given the current alignment             */
   /*=======================================================================*/

{
   int first, last, offset;
   int i, j, t;
   MlistEl *curr;
   int  bPrint = FALSE;			/* BT 3/31/97 - Debug */
   int  nSeq;

   if(B->F->nvFragCnts == NULL) 
     {
       NEWPP(B->F->nvFragCnts, B->IP->nNumMotifTypes, double);
       for(t = 0; t < B->IP->nNumMotifTypes; t++)  
	 {
	   NEWP(B->F->nvFragCnts[t],B->F->nMaxLen[t],double);
	   for(j = 0; j < B->F->nMaxLen[t]; j++) 
	     NEW(B->F->nvFragCnts[t][j], B->IP->nAlphaLen, double);
	 }
     }
   else
     {
       for(t = 0; t < B->IP->nNumMotifTypes; t++) 
         for(j = 0; j < B->F->nMaxLen[t]; j++)
	   for(i = 0; i < B->IP->nAlphaLen; i++)
	     B->F->nvFragCnts[t][j][i] = 0.0;
     }

   for(t = 0; t < B->IP->nNumMotifTypes; t++) 
     {
       curr = M[t]->Motifs;
       first = FirstCol(B->F, t);
       last = LastCol(B->F, t);
       while(curr != NULL) 
	 {
	   offset = curr->pos - first; 
	   if(curr->RevComp) 
	     /*            offset += B->F->nMaxLen[t] - 1; */
	     offset = curr->pos + last; 		/* BT 4/7/97 */
	   for(j = 0; j < B->F->nMaxLen[t]; j++) 
	     {
	       if( ! (B->F->nColMask[t][j] == COL_BLOCKED || B->F->nColMask[t][j] == COL_DONT_USE) ) 
		 { 
		   if(!curr->RevComp && (offset + j >= 0) &&	/* BT 4/7/97 */
		      (offset + j < B->IP->nSeqLen) &&
		      ((offset + j) >= curr->left) &&	/* BT 4/2/97 */
		      ((offset + j) <= curr->right) ) 
		     {
		       i = CH2INT(offset+j, B);
		       if( i >= B->IP->nAlphaLen )
			 {
			   bPrint = TRUE;
			   printf(" printFragCnts: 1\n" );
			   printFragVec( B );
			   printf( "\n" );
			   printFragCnts( B );
			 }
		       if( B->WT == NULL )			 
			 B->F->nvFragCnts[t][j][i]++;
		       else
			 {
			   nSeq = SequenceFromPosition( B, offset + j );
			   B->F->nvFragCnts[t][j][i] += GetSeqWeight( B, nSeq, t );
			 }
		     }
		   else if(curr->RevComp && (offset - j >= 0) &&	/* BT 4/7/97 */
			   (offset - j < B->IP->nSeqLen) &&
			   ((offset - j) >= curr->left) &&	/* BT 4/2/97 */
			   ((offset - j) <= curr->right) ) 
		     {
		       i = CH2INTCOMP(offset-j, B);
		       if( i >= B->IP->nAlphaLen )
			 {
			   bPrint = TRUE;
			   printf(" printFragCnts: 2\n" );
			   printFragVec( B );
			   printf( "\n" );
			   printFragCnts( B );
			 }
		       if( B->WT == NULL )			 
			 B->F->nvFragCnts[t][j][i]++;
		       else
			 {
			   nSeq = SequenceFromPosition( B, offset + j );
			   B->F->nvFragCnts[t][j][i] += GetSeqWeight( B, nSeq, t );
			 }
		     }
		   else
		     {
		       printf(" printFragCnts: 3\n" );
		       printf( "Mtype %d\npos %d\nseq_num %d\nleft %d\nright %d\nrev %d\n",
			       curr->Mtype, curr->pos, curr->seq_num, curr->left, 
			       curr->right, curr->RevComp );
		       printFragVec( B );
		       printf( "\n" );
		       printFragCnts( B );
		     }
		 } 
	     }   
	   curr = curr->next;
	 }
     }
}   

void printFragCnts(Model B)

   /*========================================================================*/
   /* FUNCTION NAME : printFragCnts                                          */
   /* DESCRIPTION   : Prints out the counts for all of the possible motif    */
   /*                 columns given the current alignment                    */
   /*========================================================================*/
   
{
   int j, t, i;

   for(t = 0; t < B->IP->nNumMotifTypes; t++) {
      for(i = 0; i < B->IP->nAlphaLen; i++) {
         for(j = 0; j < B->F->nMaxLen[t]; j++)
            printf("%3d", (int) B->F->nvFragCnts[t][j][i]);
         printf("\n");
      }
      printf("\n");
   }
}

/**********************    PossFragStartPos  ****************************/
/*                                                                      */
/* DESCRIPTION   : returns a true or false value indicating whether or  */
/*                 not the sequence position pointed to by n is a       */
/*                 possible motif starting position given fragmentation */
/*======================================================================*/

short PossFragStartPos(PoSition **Pos, int n, int t, Model B)
{
   int      j;
   int      limit;
   PoSition *Pos_t;  /* Pos_t=Pos[t] */

   if( n < 0 )
     return FALSE;

   Pos_t=Pos[t];
   if(!B->IP->is_defined[cl_F]) 
     {     
       if( B->IP->is_defined[cl_v] ) 	/* BT 12/22/99 */   
	 limit =  B->F->FragWidth[t] - 
	   (int) ((1.0 - B->IP->glOverlapParam) * B->IP->nMotifLen[t]) + 1;
       else
	 limit = B->F->FragWidth[t] - B->IP->nMotifLen[t] + 1;
	   
       for(j = 0; j < limit; j++ )
	 {
	   if(!Pos_t[n + j].nPossStartPos || 
	      (Pos_t[n + j].nInMotif && !Pos_t[n].nMotifStartPos))
	     return FALSE;  
	 }
     }
   else
     { /* not fragmentation */
       if( !Pos_t[n].nPossStartPos || 
	   (Pos_t[n].nInMotif && !Pos_t[n].nMotifStartPos))
	 return FALSE;
     }
 
  return TRUE;
}

void StoreMask(Model B)
 
   /*==========================================================*/
   /* FUNCTION NAME : StoreMask                                */
   /* DESCRIPTION   : Stores the old fragmentation column mask */
   /*==========================================================*/
{
   int i, t;

   for(t = 0; t < B->IP->nNumMotifTypes; t++) 
     {
       B->F->nOldFragWidth[t] = B->F->FragWidth[t];
       for(i = 0; i < B->F->nMaxLen[t]; i++)
         B->F->nOldColMask[t][i] = B->F->nColMask[t][i];
     }
}


/*********************  initMask *****************************************/
/*                                                                       */
/* DESCRIPTION   : Initializes the column mask used in fragmentation to a*/
/*                 contiguous block of columns                           */
/*=======================================================================*/

void initMask(Model B)

{
   int j, t;
  
   for(t = 0; t < B->IP->nNumMotifTypes; t++) 
     {
       for( j = 0; j < B->F->nMaxLen[t]; j++ )   /* BT 07/26/2001 */
	 {
	   B->F->nColMask[t][j] = B->F->fragInitMask[t][j];
	 }
       /*      i = (B->F->nMaxLen[t] - B->IP->nMotifLen[t]) / 2;
      for(j = 0; j < i; j++)
         B->F->nColMask[t][j] = 0;
      for(j = i; j < i + B->IP->nMotifLen[t]; j++)
         B->F->nColMask[t][j] = 1;
      for(j = i + B->IP->nMotifLen[t]; j < B->F->nMaxLen[t]; j++)
      B->F->nColMask[t][j] = 0; 
      B->F->FragWidth[t] = B->IP->nMotifLen[t]; */

      B->F->FragWidth[t] = LastCol( B->F, t ) - FirstCol( B->F, t ) + 1;
      for(j = 0; j < B->IP->nMotifLen[t]; j++)
	B->F->fragPos[t][j] = 0;
      SetFragPos( B, t );
   }
   
}


/* BT 12/23/99 */
void SetFragPos( Model B, int t )
{
  Ftype  F;
  IPtype IP;
  int    first;
  int    p;
  int    i = 0;
  
  IP = B->IP;
  F = B->F;

  first = FirstCol(F, t);
  p = first;
  while( p != -1 )
    {
      F->fragPos[t][i] = p - first;
      i++;
      p = NextCol( F, t, p );
    }
}
