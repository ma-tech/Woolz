#ifndef ALCPROTO_H
#define ALCPROTO_H
#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _AlcProto_h[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libAlc/AlcProto.h
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
* 
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be
* useful but WITHOUT ANY WARRANTY; without even the implied
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the GNU General Public License for more
* details.
*
* You should have received a copy of the GNU General Public
* License along with this program; if not, write to the Free
* Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
* Boston, MA  02110-1301, USA.
* \brief        Function prototypes for the Woolz type allocation 
*		library.
* \todo		-
* \bug          None known.
*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/************************************************************************
* AlcAlloc.c
************************************************************************/
extern void			*AlcCalloc(
				  size_t elCount,
				  size_t elSize);
extern void			*AlcMalloc(
				  size_t byteCount);
extern void			*AlcRealloc(
				  void *givendata,
				  size_t byteCount);
extern void			AlcFree(
				  void *data);

/************************************************************************
* AlcArray.c
************************************************************************/
extern AlcErrno			AlcBit1Calloc(
				  unsigned char **dest,
				  size_t mElem);
extern AlcErrno			AlcChar1Calloc(
				  char **dest,
				  size_t mElem);
extern AlcErrno			AlcUnchar1Calloc(
				  unsigned char **dest,
				  size_t mElem);
extern AlcErrno			AlcShort1Calloc(
				  short **dest,
				  size_t mElem);
extern AlcErrno			AlcInt1Calloc(
				  int **dest,
				  size_t mElem);
extern AlcErrno			AlcInt1Calloc(
				  int **dest,
				  size_t mElem);
extern AlcErrno			AlcFloat1Calloc(
				  float **dest,
				  size_t mElem);
extern AlcErrno			AlcDouble1Calloc(
				  double **dest,
				  size_t mElem);
extern AlcErrno			AlcBit1Malloc(
				  unsigned char **dest,
				  size_t mElem);
extern AlcErrno			AlcChar1Malloc(
				  char **dest,
				  size_t mElem);
extern AlcErrno			AlcUnchar1Malloc(
				  unsigned char **dest,
				  size_t mElem);
extern AlcErrno			AlcShort1Malloc(
				  short **dest,
				  size_t mElem);
extern AlcErrno			AlcInt1Malloc(
				  int **dest,
				  size_t mElem);
extern AlcErrno			AlcFloat1Malloc(
				  float **dest,
				  size_t mElem);
extern AlcErrno			AlcDouble1Malloc(
				  double **dest,
				  size_t mElem);
extern AlcErrno			AlcBit2Calloc(
				  unsigned char ***dest,
				  size_t mElem,
				  size_t nElem);
extern AlcErrno			AlcChar2Calloc(
				  char ***dest,
				  size_t mElem,
				  size_t nElem);
extern AlcErrno			AlcUnchar2Calloc(
				  unsigned char ***dest,
				  size_t mElem,
				  size_t nElem);
extern AlcErrno			AlcShort2Calloc(
				  short ***dest,
				  size_t mElem,
				  size_t nElem);
extern AlcErrno			AlcInt2Calloc(
				  int ***dest,
				  size_t mElem,
				  size_t nElem);
extern AlcErrno			AlcFloat2Calloc(
				  float ***,
				  size_t mElem,
				  size_t nElem);
extern AlcErrno			AlcDouble2Calloc(
				  double ***dest,
				  size_t mElem,
				  size_t nElem);
extern AlcErrno			AlcPtr2Calloc(
				  void ****dest,
				  size_t mElem,
				  size_t nElem);
extern AlcErrno			AlcBit2Malloc(
				  unsigned char ***dest,
				  size_t mElem,
				  size_t nElem);
extern AlcErrno			AlcChar2Malloc(
				  char ***dest,
				  size_t mElem,
				  size_t nElem);
extern AlcErrno			AlcUnchar2Malloc(
				  unsigned char ***dest,
				  size_t mElem,
				  size_t nElem);
extern AlcErrno			AlcShort2Malloc(
				  short ***dest,
				  size_t mElem,
				  size_t nElem);
extern AlcErrno			AlcInt2Malloc(
				  int ***dest,
				  size_t mElem,
				  size_t nElem);
extern AlcErrno			AlcFloat2Malloc(
				  float ***dest,
				  size_t mElem,
				  size_t nElem);
extern AlcErrno			AlcDouble2Malloc(
				  double ***dest,
				  size_t mElem,
				  size_t nElem);
extern AlcErrno        		AlcSymChar2Calloc(
				  char ***dest,
				  size_t nElem);
extern AlcErrno        		AlcSymUnchar2Calloc(
				  unsigned char ***dest,
				  size_t nElem);
extern AlcErrno        		AlcSymShort2Calloc(
				  short ***dest,
				  size_t nElem);
extern AlcErrno        		AlcSymInt2Calloc(
				  int ***dest,
				  size_t nElem);
extern AlcErrno        		AlcSymFloat2Calloc(
				  float ***dest,
				  size_t nElem);
extern AlcErrno        		AlcSymDouble2Calloc(
				  double ***dest,
				  size_t nElem);
extern AlcErrno        		AlcSymChar2Malloc(
				  char ***dest,
				  size_t nElem);
extern AlcErrno        		AlcSymUnchar2Malloc(
				  unsigned char ***dest,
				  size_t nElem);
extern AlcErrno        		AlcSymShort2Malloc(
				  short ***dest,
				  size_t nElem);
extern AlcErrno        		AlcSymInt2Malloc(
				  int ***dest,
				  size_t nElem);
extern AlcErrno        		AlcSymFloat2Malloc(
				  float ***dest,
				  size_t nElem);
extern AlcErrno        		AlcSymDouble2Malloc(
				  double ***dest,
				  size_t nElem);
extern AlcErrno			Alc2Free(
				  void **dat);
extern AlcErrno			AlcBit2Free(
				  unsigned char **dat);
extern AlcErrno			AlcChar2Free(
				  char **dat);
extern AlcErrno			AlcUnchar2Free(
				  unsigned char **dat);
extern AlcErrno			AlcShort2Free(
				  short **dat);
extern AlcErrno			AlcInt2Free(
				  int **dat);
extern AlcErrno			AlcFloat2Free(
				  float **dat);
extern AlcErrno			AlcDouble2Free(
				  double **dat);
extern AlcErrno			AlcBit3Calloc(
				  unsigned char ****dest,
				  size_t mElem,
				  size_t nElem,
				  size_t oElem);
extern AlcErrno			AlcChar3Calloc(
				  char ****dest,
				  size_t mElem,
				  size_t nElem,
				  size_t oElem);
extern AlcErrno			AlcUnchar3Calloc(
				  unsigned char ****dest,
				  size_t mElem,
				  size_t nElem,
				  size_t oElem);
extern AlcErrno			AlcShort3Calloc(
				  short ****dest,
				  size_t mElem,
				  size_t nElem,
				  size_t oElem);
extern AlcErrno			AlcInt3Calloc(
				  int ****dest,
				  size_t mElem,
				  size_t nElem,
				  size_t oElem);
extern AlcErrno			AlcFloat3Calloc(
				  float ****dest,
				  size_t mElem,
				  size_t nElem,
				  size_t oElem);
extern AlcErrno			AlcDouble3Calloc(
				  double ****dest,
				  size_t mElem,
				  size_t nElem,
				  size_t oElem);
extern AlcErrno			AlcBit3Malloc(
				  unsigned char ****dest,
				  size_t mElem,
				  size_t nElem,
				  size_t oElem);
extern AlcErrno			AlcChar3Malloc(
				  char ****dest,
				  size_t mElem,
				  size_t nElem,
				  size_t oElem);
extern AlcErrno			AlcUnchar3Malloc(
				  unsigned char ****dest,
				  size_t mElem,
				  size_t nElem,
				  size_t oElem);
extern AlcErrno			AlcShort3Malloc(
				  short ****dest,
				  size_t mElem,
				  size_t nElem,
				  size_t oElem);
extern AlcErrno			AlcInt3Malloc(
				  int ****dest,
				  size_t mElem,
				  size_t nElem,
				  size_t oElem);
extern AlcErrno			AlcFloat3Malloc(
				  float ****dest,
				  size_t mElem,
				  size_t nElem,
				  size_t oElem);
extern AlcErrno			AlcDouble3Malloc(
				  double ****dest,
				  size_t mElem,
				  size_t nElem,
				  size_t oElem);
extern AlcErrno			AlcPtr3Calloc(
				  void *****dest,
				  size_t mElem,
				  size_t nElem,
			          size_t oElem);
extern AlcErrno			Alc3Free(
				  void ***dat);
extern AlcErrno			AlcBit3Free(
				  unsigned char ***dat);
extern AlcErrno			AlcChar3Free(
				  char ***dat);
extern AlcErrno			AlcUnchar3Free(
				  unsigned char ***dat);
extern AlcErrno			AlcShort3Free(
				  short ***dat);
extern AlcErrno			AlcInt3Free(
				  int ***dat);
extern AlcErrno			AlcFloat3Free(
				  float ***dat);
extern AlcErrno			AlcDouble3Free(
				  double ***);
extern AlcErrno			AlcDouble1ReadAsci(
				  FILE *fP,
				  double **dstA,
				  size_t *dstNElem);
extern AlcErrno			AlcDouble2ReadAsci(
				  FILE *fP,
				  double ***dstA,
				  size_t *dstMElem,
				  size_t *dstNElem);
extern AlcErrno			AlcDouble1WriteAsci(
				  FILE *fP,
				  double *ar,
				  size_t nElem);
extern AlcErrno			AlcDouble2WriteAsci(
				  FILE *fP,
				  double **ar,
				  size_t mElem,
				  size_t nElem);

/************************************************************************
* AlcBlockStack.c
************************************************************************/
extern AlcBlockStack		*AlcBlockStackNew(
				  size_t nElem,
				  size_t elmSz,
				  AlcBlockStack *tBlk,
				  AlcErrno *dstErr);
extern AlcErrno			AlcBlockStackFree(
				  AlcBlockStack *blk);

/************************************************************************
* AlcCPQueue.c
************************************************************************/
extern AlcCPQQueue              *AlcCPQQueueNew(
                                  AlcErrno *dstErr);
extern AlcCPQItem               *AlcCPQItemNew(
                                  AlcCPQQueue *q,
                                  float priority,
                                  void *entry,
                                  AlcErrno *dstErr);
extern AlcErrno                 AlcCPQQueueFree(
                                  AlcCPQQueue *q);
extern void                     AlcCPQItemFree(
                                  AlcCPQQueue *q,
                                  AlcCPQItem *item);
extern AlcErrno                 AlcCPQEntryInsert(
                                  AlcCPQQueue *q,
                                  float priority,
                                  void *entry);
extern AlcErrno                 AlcCPQItemInsert(
                                  AlcCPQQueue *q,
                                  AlcCPQItem *item);
extern AlcCPQItem               *AlcCPQItemUnlink(
                                  AlcCPQQueue *q);

/************************************************************************
* AlcDLPList.c
************************************************************************/
extern AlcDLPList 		*AlcDLPListNew(
				  AlcErrno *dstErr);
extern AlcErrno			AlcDLPListFree(
				  AlcDLPList *list);
extern AlcDLPItem		*AlcDLPItemNew(
				  void *entry,
				  void (*freeFn)(void *),
				  AlcErrno *dstErr);
extern AlcErrno			AlcDLPListEntryAppend(
				  AlcDLPList *list,
				  AlcDLPItem *appAfter,
				  void *entry,
				  void (*freeFn)(void *));
extern AlcErrno			AlcDLPListEntryInsert(
				  AlcDLPList *list,
				  AlcDLPItem *insBefore,
				  void *entry,
				  void (*freeFn)(void *));
extern AlcDLPItem		*AlcDLPItemUnlink(
				  AlcDLPList *list,
				  AlcDLPItem *item,
				  int freeItem,
				  AlcErrno *dstErr);
extern AlcErrno			AlcDLPItemAppend(
				  AlcDLPList *list,
				  AlcDLPItem *appAfter,
				  AlcDLPItem *item);
extern AlcErrno			AlcDLPItemInsert(
				  AlcDLPList *list,
				  AlcDLPItem *insBefore,
				  AlcDLPItem *item);
extern AlcErrno			AlcDLPItemFree(
				  AlcDLPItem *item);
extern	int			AlcDLPListCount(
				  AlcDLPList *list,
				  AlcErrno *dstErr);
extern AlcDLPItem  		*AlcDLPListIterate(
				  AlcDLPList *list,
				  AlcDLPItem *item,
				  AlcDirection dir,
				  int (*iterFn)(AlcDLPList *,
				  		AlcDLPItem *,
						void *),
				  void *iterData,
				  AlcErrno *dstErr);
extern AlcErrno			AlcDLPListSort(
				  AlcDLPList *list,
				  int (*entryCompFn)(void *,
				                     void *));

/************************************************************************
* AlcFreeStack.c
************************************************************************/
extern void            		*AlcFreeStackPush(
				  void *prev,
				  void *data,
				  AlcErrno *dstErr);
extern void			*AlcFreeStackPop(
				  void *prev,
				  void **dstData,
				  AlcErrno *dstErr);
extern AlcErrno			AlcFreeStackFree(
				  void *stack);

/************************************************************************
* AlcHashTable.c
************************************************************************/
extern AlcHashTable		*AlcHashTableNew(
				  size_t tableSz,
				  int (*keyCmp)(void *, void *),
				  unsigned (*hashFn)(void *),
				  AlcErrno *dstErr);
extern AlcHashItem		*AlcHashItemNew(
				  void *entry,
				  void (*freeFn)(void *),
				  void *key,
				  AlcErrno *dstErr);
extern AlcErrno			AlcHashTableFree(
				  AlcHashTable *hTbl);
extern AlcErrno			AlcHashTableEntryInsert(
				  AlcHashTable *hTbl,
				  void *key,
				  void *entry,
				  void (*freeFn)(void *));
extern AlcErrno			AlcHashItemUnlink(
				  AlcHashTable *hTbl,
				  AlcHashItem *rItem,
				  int freeItem);
extern AlcErrno        		AlcHashTableUnlinkAll(
				  AlcHashTable *hTbl,
                                  int (*testFn)(AlcHashTable *,
                                                AlcHashItem *, void *),
                                  void *fnData,
				  int freeItems);
extern AlcErrno 		AlcHashItemInsert(
				  AlcHashTable *hTbl,
				  AlcHashItem *newItem);
extern AlcErrno			AlcHashItemFree(
				  AlcHashItem *item);
extern size_t			AlcHashTableCount(
				  AlcHashTable *hTbl,
				  AlcErrno *dstErr);
extern AlcHashItem		*AlcHashTableIterate(
				  AlcHashTable *hTbl,
				  AlcDirection dir,
				  int (*iterFn)(AlcHashTable *,
				                AlcHashItem *, void *),
				  void *iterData, AlcErrno *dstErr);
extern AlcHashItem		*AlcHashItemGet(
				  AlcHashTable *hTbl,
				  void *key,
				  AlcErrno *dstErr);
extern int			AlcHashItemOrder(
				  AlcHashTable *hTbl,
				  AlcHashItem *item0,
				  AlcHashItem *item1);
/************************************************************************
* AlcKDTree.c
************************************************************************/
extern AlcKDTTree	       *AlcKDTTreeNew(
				   AlcPointType type,
				  int dim,
				  double tol,
				  size_t nNodes,
				  AlcErrno *dstErr);
extern AlcErrno 	  	AlcKDTTreeFree(
				    AlcKDTTree *tree);
extern AlcKDTNode		*AlcKDTNodeNew(
				  AlcKDTTree *tree,
				  AlcKDTNode *parent,
				  AlcPointP key,
				  int cmp,
				  AlcErrno *dstErr);
extern void			AlcKDTNodeFree(
				  AlcKDTTree *tree,
				  AlcKDTNode *node);
extern int			AlcKDTTreeFacts(
				  AlcKDTTree *tree,
				  FILE *fP);
extern AlcKDTNode		*AlcKDTInsert(
				  AlcKDTTree *tree,
				  void *keyVal,
				  AlcKDTNode **dstFndNod,
				  AlcErrno *dstErr);
extern AlcKDTNode		*AlcKDTGetMatch(
				  AlcKDTTree *tree,
				  void *keyVal,
				  AlcErrno *dstErr);
extern AlcKDTNode	        *AlcKDTGetLeaf(
				  AlcKDTTree *tree,
				  AlcKDTNode *node,
				  AlcPointP key);
extern AlcKDTNode	       *AlcKDTGetNN(
				  AlcKDTTree *tree,
				  void *keyVal,
				  double minDist,
				  double *dstNNDist,
				  AlcErrno *dstErr);

/************************************************************************
* AlcString.c
************************************************************************/
extern char			*AlcStrDup(
				  const char *srcStr);

/************************************************************************
* AlcVector.c
************************************************************************/
extern AlcVector		*AlcVectorNew(
				  size_t elmCnt,
				  size_t elmSz,
				  size_t blkSz,
                              	  AlcErrno *dstErr);
extern AlcErrno			AlcVectorFree(
				  AlcVector *vec);
extern AlcErrno			AlcVectorExtend(
				  AlcVector *vec,
				  size_t elmCnt);
extern void			*AlcVectorItemGet(
				  AlcVector *vec,
				  size_t idx);
extern void			*AlcVectorExtendAndGet(
				  AlcVector *vec,
				  size_t idx);
extern size_t			AlcVectorCount(
				  AlcVector *vec);
extern void			AlcVectorSetArray1D(
				  AlcVector *vec,
				  size_t fIdx,
				  size_t lIdx,
				  void *aM);
extern void			*AlcVectorToArray1D(
				  AlcVector *vec,
				  size_t fIdx,
				  size_t lIdx,
				  AlcErrno *dstErr);
void				**AlcVectorToArray2D(
				  AlcVector *vec,
				  size_t fIdx,
				  size_t lIdx,
				  size_t nR,
				  size_t nC,
				  AlcErrno *dstErr);
#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* ALCPROTO_H */
