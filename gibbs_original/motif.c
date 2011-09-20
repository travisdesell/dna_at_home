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
/* $Id: motif.c,v 1.7 2009/04/23 18:43:54 Bill Exp $             */
/*                                                                        */
/* Author :       Eric C. Rouchka, July  3, 1996                          */
/*                Jun Zhu, October 20, 1996                               */
/*                                                                        */
/*                                                                        */
/**************************************************************************/

#include "motif.h"


/*****************  copy_motif_num  *****************************/
/*                                                              */
/* DESCRIPTION   : returns a vector containing the number of    */
/*                 motifs for each motif type in both the fwd   */
/*                 and reverse directions                       */
/*==============================================================*/

int **copy_motif_num(IPtype IP)
{
    int i, **temp;

    NEWP(temp, IP->nNumMotifTypes, int);
    for(i = 0; i < IP->nNumMotifTypes; i++) {
       NEW(temp[i], 2, int);
       temp[i][FORWARD] = IP->nNumMotifs[i][FORWARD];
       temp[i][REVERSE] = IP->nNumMotifs[i][REVERSE];
    }
    return temp;
}


/******************  zero_motifs *********************************/
/*                                                               */
/* DESCRIPTION   : Sets the number of forward and reverse motifs */
/*                 for all motif types equal to zero             */
/*===============================================================*/

void zero_motifs(IPtype IP)
{
   int i;
   for(i = 0; i < IP->nNumMotifTypes; i++) {
      IP->nNumMotifs[i][FORWARD] = 0;
      IP->nNumMotifs[i][REVERSE] = 0;
   }
}

void reset_motif_num(int **nNumMotifs, IPtype IP)

   /*================================================================*/
   /* FUNCTION NAME : reset_motif_num                                */
   /* DESCRIPTION   : Resets the number of forward and reverse motifs*/
   /*                 for all motif types                            */
   /*================================================================*/
{
   int i;

   for(i = 0; i < IP->nNumMotifTypes; i++) {
      IP->nNumMotifs[i][FORWARD] = nNumMotifs[i][FORWARD];
      IP->nNumMotifs[i][REVERSE] = nNumMotifs[i][REVERSE];
    } 
}

/* BT 4/4/97 */

short PossibleOverlapsPrevMotif( Model B, int curPos, int nMotifType,  PoSition **Pos )

   /*=================================================================*/
   /* FUNCTION NAME : OverlapsPrevMotif                               */
   /* DESCRIPTION   : Checks to see if a motif at the current         */
   /*                 position overlaps any of the previous motifs    */
   /*                 that have already been set.                     */
   /*=================================================================*/
{
     int t, i;
     int nEnd;
    
     nEnd = MotifWidth( B, nMotifType );   /* BT 12/7/2000 */
     for( t = 0; t < B->IP->nNumMotifTypes; t++ )
     {
       for( i = curPos; i < curPos+nEnd; i++ )
         {
	   if( (t != nMotifType) && (/* Pos[t][i].nMotifStartPos ||*/ Pos[t][i].nInMotif) )
	     return TRUE;    
         }
     }
     return FALSE; 
}


short OverlapsPrevMotif(int CurPos, int nNumMotif,
                        PoSition **Pos)


   /*=================================================================*/
   /* FUNCTION NAME : OverlapsPrevMotif                               */
   /* DESCRIPTION   : Checks to see if the current position overlaps  */
   /*                 any of the previous motifs that have already    */
   /*                 been set                                        */
   /*=================================================================*/

{
  int t;

  for(t = 0; t < nNumMotif; t++)
      if(Pos[t][CurPos].nInMotif && !Pos[t][CurPos].nMotifStartPos ) 
	return TRUE;
  return FALSE;
}


Mlist initializeMotifs(IPtype IP)

   /**************************************************/
   /* FUNCTION NAME : initializeMotifs               */
   /* DESCRIPTION : initializes the motifs           */
   /**************************************************/

{
   int t;
   Mlist M;

   NEWP(M, IP->nNumMotifTypes, Mlist_struct);
   for(t = 0; t < IP->nNumMotifTypes; t++) {
      NEW(M[t], 1, Mlist_struct);
      M[t]->nMotifLen = IP->nMotifLen[t];
      M[t]->nNumMotifs = 0;
      M[t]->Motifs = NULL;
   }
   return M;
}

/******************  void  set_in_motif  *************************/
/*                                                               */
/* DESCRIPTION   : Sets the position pointed to be i in motif t  */
/*                 to be a motif starting position of type t     */
/*===============================================================*/

void set_in_motif(PoSition **Pos, int i, Model B, int t, short RevComp)

{
   register int n;   /* general register    */
   PoSition *Pos_t;  /* Pos_t=Pos[t]        */
   int      nEnd;

   Pos_t=Pos[t];
   
   nEnd = MotifWidth( B, t );

   for(n=0; n < nEnd; n++)  /* BT 9/11/98 */
   {
     Pos_t[i+n].nInMotif = TRUE;
     if( Pos_t[i + n].bEndOfSeq )
        break;
   }
   Pos_t[i].nMotifStartPos = TRUE;
   Pos_t[i].RevComp = RevComp; 
}

/******************  void  not_in_motif  ***************************/
/*                                                                 */
/* DESCRIPTION   : delete the position pointed to be i in motif t  */
/*                 to be a motif starting position of type t       */
/*=================================================================*/

void not_in_motif(PoSition **Pos, int i, Model B,int t)

{
   register int n;   /* general register        */
   PoSition *Pos_t;  /* Pos_t=Pos[t]            */
   int      nEnd;

   nEnd = MotifWidth( B, t );

   Pos_t=Pos[t];
   for(n=0; n < nEnd; n++)
   {
     Pos_t[i+n].nInMotif = FALSE;
     if( Pos_t[i + n].bEndOfSeq )	/* BT 4/23/97 */
        break;
   } 
   Pos_t[i].nMotifStartPos = FALSE;
}


int overlap(int i, PoSition **Pos, int *nMotifLen, int *newpos, int t, int t1)
{
   /*===============================================*/
   /* FUNCTION NAME : overlap                       */
   /*                                               */
   /* DESCRIPTION : checks to see if the current    */
   /*               position overlaps another motif */
   /*               element                         */
   /*===============================================*/

   int flag = FALSE;
   int j;
  
   for(j = 0; ((j < nMotifLen[t]) && (!flag)); j++) {
      if(Pos[t1][i + j].nMotifStartPos) {
         flag = TRUE;
         (*newpos) = i + j;
      }
   }
   return(flag);
}


short PossibleMotifPos(int i, PoSition **Pos, IPtype IP, int t, Model B )     /* BT 7/5/97 */
{
   /*===============================================*/
   /* FUNCTION NAME : PossibleMotifPos              */
   /*                                               */
   /* DESCRIPTION : checks to see if the current    */
   /*               position is a possible motif    */
   /*               starting pos for motif type t   */
   /*===============================================*/

   short     flag = TRUE;
   int       j;
   int       t1;
   PoSition *Pos_t;  /* Pos_t=Pos[t]        */  

   for( t1 = 0; (t1 < IP->nNumMotifTypes) && flag; t1++ )
     {
       Pos_t = Pos[t1];
       for(j = 0; ((j < MotifWidth( B, t ) ) && flag); j++) 
	 {
	   if(Pos_t[i + j].nMotifStartPos || (! Pos_t[i + j].nPossStartPos)) 
	     {
	       flag = FALSE;
	       break;
	     }
	 }
     }

   return(flag);
}


/*********************** set_motif_info  ************************/
/*                                                              */
/* create a list of motif                                       */
/****************************************************************/


Mlist set_motif_info(IPtype IP, int **StartPos, Stype Seq)  
{
   Mlist M;
   int m, t;
   int *StartPos_t;  /* StartPos_t=StartPos[t] */

   NEWP(M, IP->nNumMotifTypes, Mlist_struct);
   for(t = 0; t < IP->nNumMotifTypes; t++) {
       StartPos_t=StartPos[t];
       NEW(M[t], 1, Mlist_struct);
       M[t]->nNumMotifs = 0;
       M[t]->nMotifLen = IP->nMotifLen[t];

       M[t]->Motifs = NULL;
       for(m = 0; m < IP->nNumMotifs[t][FORWARD]; m++)
          add_motif(IP, Seq, StartPos_t[m], M, t, FALSE);
       for(m = IP->nNumMotifs[t][FORWARD]; m < IP->nNumMotifs[t][REVERSE]+
                                              IP->nNumMotifs[t][FORWARD]; m++)
          add_motif(IP, Seq, StartPos_t[m], M, t, TRUE);
   }
   return M;
}



/**********************  void add_motif **********************/
/*                                                           */
/* DESCRIPTION : This function adds a motif to the alignment */
/*===========================================================*/

void add_motif(IPtype IP, Stype Seq, int pos, Mlist M, int t, short RevComp)
{
    MlistEl *curr;
    int i, FOUND;
    int tmppos;

    NEW(curr, 1, MlistEl);
    curr->pos = pos;
    curr->RevComp = RevComp;
    for(FOUND = FALSE, i=0; !FOUND; i++) 
      {                                      /** Find the current sequence **/
	if((*Seq->nvEndLocs)[i] > curr->pos) /** number according to the   **/
          FOUND = TRUE;                      /** position                  **/
      }
    curr->seq_num = i - 1;
    if(i > 1)                                  /** Set the max left shift **/
       curr->left = (*Seq->nvEndLocs)[i - 2];  /** to the left end of the **/
    else                                       /** current string         **/
       curr->left = 0;

    /*    curr->right = (*Seq->nvEndLocs)[i-1] - M[t]->nMotifLen; */ /* BT 4/2/97 */
    curr->right = (*Seq->nvEndLocs)[i-1] - 1;  
    /* for amino acid */
    
    /*    if(IP->nAlpha==20){*/
    if(TRUE){
      /* Readjust Maxleft so no processed are included */  
      for(tmppos = curr->left; tmppos < pos; tmppos++) 
	{ 
	  if( (*Seq->ProcessedSTR)[tmppos] == 'X' ||
	      (*Seq->ProcessedSTR)[tmppos] == 'x' || 
	      (*Seq->ProcessedSTR)[tmppos] == 'u') 
	    curr->left = tmppos+1;
	  else if( IP->nAlphaLen == 4 && 
		 (((*Seq->ProcessedSTR)[tmppos] == 'N') || 
		  ((*Seq->ProcessedSTR)[tmppos] == 'n')) )
          curr->left = tmppos+1; 			/* BT 3/31/96 */
	}                                                  
      /* Readjust Max right so no processed ones are included */
      for(tmppos = pos+M[t]->nMotifLen;              
	  tmppos <=curr->right; tmppos++) 
	{                
	  if((*Seq->ProcessedSTR)[tmppos] == 'X' ||
	     (*Seq->ProcessedSTR)[tmppos] == 'x' ||
	     (*Seq->ProcessedSTR)[tmppos] == 'u')  
	    {         
	      curr->right = tmppos - 1;                      
	      break;
	    }
	  else if( IP->nAlphaLen == 4 &&
		   (((*Seq->ProcessedSTR)[tmppos] == 'N') || 
		    ((*Seq->ProcessedSTR)[tmppos] == 'n')) )
	    { 
	      curr->right = tmppos - 1; 			/* BT 3/31/96 */      
	      break; 
	    }
	}
    }
    curr->next = M[t]->Motifs;
    curr->Mtype = t;
    M[t]->Motifs = curr;
    (M[t]->nNumMotifs)++;
}



/************** void delete_motif  ************************/
/*                                                        */
/* DESCRIPTION : This function deletes a motif given its  */
/*               starting position                        */
/**********************************************************/

void delete_motif(Model B, int pos, Mlist M, int t) 		/* BT 3/31/97 */

{
   MlistEl *curr, *prev;
   char    msg[32];
   int     pos2;
   int     seq;

   curr = M[t]->Motifs;
   prev = M[t]->Motifs;
   /* search through the list */
   if(curr != NULL) 
     {
       (M[t]->nNumMotifs)--;
       if(M[t]->nNumMotifs == 0)
	 {
	   free(curr);   /* BT 1/2/98 */
	   M[t]->Motifs = NULL;
	 }
       else 
	 {
	   while(curr->pos != pos) 
	     {
	       prev = curr;
	       curr = curr->next;
	       if(curr == NULL) 
		 {                                               /* The motif is not */
		   seq = SequenceFromPosition( B, pos );
		   pos2 = pos - SequenceStartPos( B, seq );
		   sprintf( msg, "Motif %d not found seq %d pos %d model %d\n", 
			    pos, seq, pos2, t);                     /* found, so exit   */
		   print_motifs(B, M, t);                        /* with an error,   */  /* BT 7/28/97 */
		   p_internal_error( msg );                               /* printing motifs  */
		 }
	     }
	   /* the current element is not the last one */
	   if(curr->next != NULL) 
	     {
	       if(curr == M[t]->Motifs) /* if it is the first one */
		 M[t]->Motifs = curr->next;
	       else 
		 {
		   prev->next = curr->next;
		   /* free(curr); */    /* BT 9/22/97 */
		 }
	     }
	   else
	     { /*the current element is the last one */ 
	       if(curr == M[t]->Motifs) /* if it is the first one */
		 M[t]->Motifs = NULL;
	       else
		 {
		   prev->next = NULL;
		   /* free(curr); */     /* BT 9/22/97 */
		 }
	     }
	   if( curr != NULL)       /* BT 9/22/97 */
	     free( curr );   
	 }
     }
}

void print_motifs(Model B, Mlist M, int t)

   /************************************************/
   /* FUNCTION NAME : print_motifs                 */
   /* DESCRIPTION : Prints out all of the motifs   */
   /************************************************/

{
  MlistEl *curr;
  int i;
  printf("MOTIF %d\n", t);
  printf("---------------\n");
  curr = M[t]->Motifs;
  printf("Num Motifs = %d\n", M[t]->nNumMotifs);
  while(curr != NULL) 
    {
      printf("PoSition %d Rev? %d Seq %d Pos %d\n", 
	     curr->pos, curr->RevComp, curr->seq_num, 
	     curr->pos - SequenceStartPos( B, curr->seq_num ));
      for(i = 0; i < B->IP->nMotifLen[t]; i++)
	printf("%c", curr_ch((*B->Seq->R)[i+ curr->pos], B->IP->nAlphaLen, TRUE ) );
      printf("\n");
      curr = curr->next;
    }
}


void free_motifs(Model B, Mlist M)
{
   MlistEl *curr, *next;
   int t;

   for(t = 0; t < B->IP->nNumMotifTypes; t++) {
      curr = M[t]->Motifs;
      while(curr != NULL) {
         next = curr->next;
         free(curr);
         curr = next;
      }
      free(M[t]);   /* BT 8/27/97 */
   }

   free(M);    /* BT 8/27/97 */
   M = NULL;
}


int AnyOverlap(Model B, int n, int t, PoSition **Pos, int *newpos, int *typ)
{
   /*===============================================*/
   /* FUNCTION NAME : overlap                       */
   /*                                               */
   /* DESCRIPTION : checks to see if the current    */
   /*               position overlaps another motif */
   /*               element                         */
   /*===============================================*/

   int j;
   int t1;
   int p;
   int seq;
   int len;
     
   seq = SequenceFromPosition( B, n );
   len = SequenceLength( B, seq );

   for(j = 0; j < MotifWidth( B, t ); j++) 
     {
       for( p = 0; p < SeqInc( B, seq ); p++ )
	 {
	   for( t1 = 0; t1 < B->IP->nNumMotifTypes; t1++ )
	     {
	       if(Pos[t1][n + (p * len) + j].nMotifStartPos) 
		 {
		   (*newpos) = n + j;
		   (*typ) = t1;
		   return( TRUE );
		 }
	     }
	 }
     }

   return( FALSE );
}
