#ifndef ALCPROTO_H
#define ALCPROTO_H
#pragma ident "MRC HGU $Id$"
/*!
* \file         AlcProto.h
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \ingroup	Alc
* \brief        Function prototypes for the MRC HGU memory allocation 
*		and fundamental type library.
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
				  int elCount,
				  int elSize);
extern void			*AlcMalloc(
				  int byteCount);
extern void			*AlcRealloc(
				  void *givendata,
				  int byteCount);
extern void			AlcFree(
				  void *data);

/************************************************************************
* AlcArray.c
************************************************************************/
extern AlcErrno	AlcBit1Calloc(unsigned char **, int),
		AlcChar1Calloc(char **, int),
		AlcUnchar1Calloc(unsigned char **, int),
		AlcShort1Calloc(short **, int),
		AlcInt1Calloc(int **, int),
		AlcInt1Calloc(int **, int),
		AlcFloat1Calloc(float **, int),
		AlcDouble1Calloc(double **, int),
		AlcBit1Malloc(unsigned char **, int),
		AlcChar1Malloc(char **, int),
		AlcUnchar1Malloc(unsigned char **, int),
		AlcShort1Malloc(short **, int),
		AlcInt1Malloc(int **, int),
		AlcFloat1Malloc(float **, int),
		AlcDouble1Malloc(double **, int),
		AlcBit2Calloc(unsigned char ***, int, int),
		AlcChar2Calloc(char ***, int, int),
		AlcUnchar2Calloc(unsigned char ***, int, int),
		AlcShort2Calloc(short ***, int, int),
		AlcInt2Calloc(int ***, int, int),
		AlcFloat2Calloc(float ***, int, int),
		AlcDouble2Calloc(double ***, int, int),
		AlcBit2Malloc(unsigned char ***, int, int),
		AlcChar2Malloc(char ***, int, int),
		AlcUnchar2Malloc(unsigned char ***, int, int),
		AlcShort2Malloc(short ***, int, int),
		AlcInt2Malloc(int ***, int, int),
		AlcFloat2Malloc(float ***, int, int),
		AlcDouble2Malloc(double ***, int, int),
		Alc2Free(void **),
		AlcBit2Free(unsigned char **),
		AlcChar2Free(char **),
		AlcUnchar2Free(unsigned char **),
		AlcShort2Free(short **),
		AlcInt2Free(int **),
		AlcFloat2Free(float **),
		AlcDouble2Free(double **),
		AlcBit3Calloc(unsigned char ****, int, int, int),
		AlcChar3Calloc(char ****, int, int, int),
		AlcUnchar3Calloc(unsigned char ****, int, int, int),
		AlcShort3Calloc(short ****, int, int, int),
		AlcInt3Calloc(int ****, int, int, int),
		AlcFloat3Calloc(float ****, int, int, int),
		AlcDouble3Calloc(double ****, int, int, int),
		AlcBit3Malloc(unsigned char ****, int, int, int),
		AlcChar3Malloc(char ****, int, int, int),
		AlcUnchar3Malloc(unsigned char ****, int, int, int),
		AlcShort3Malloc(short ****, int, int, int),
		AlcInt3Malloc(int ****, int, int, int),
		AlcFloat3Malloc(float ****, int, int, int),
		AlcDouble3Malloc(double ****, int, int, int),
		Alc3Free(void ***),
		AlcBit3Free(unsigned char ***),
		AlcChar3Free(char ***),
		AlcUnchar3Free(unsigned char ***),
		AlcShort3Free(short ***),
		AlcInt3Free(int ***),
		AlcFloat3Free(float ***),
		AlcDouble3Free(double ***);
extern AlcErrno			AlcDouble2ReadAsci(
				  FILE *fP,
				  double ***dstA,
				  int *dstMElem,
				  int *dstNElem);

/************************************************************************
* AlcBlockStack.c
************************************************************************/
extern AlcBlockStack		*AlcBlockStackNew(
				  unsigned nElem,
				  unsigned elmSz,
				  AlcBlockStack *tBlk,
				  AlcErrno *dstErr);
extern AlcErrno			AlcBlockStackFree(
				  AlcBlockStack *blk);

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
				  		AlcDLPItem *, void *),
				  void *iterData,
				  AlcErrno *dstErr);
extern AlcErrno			AlcDLPListSort(
				  AlcDLPList *list,
				  int (*entryCompFn)(void *, void *));

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
				  unsigned tableSz,
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
extern int			AlcHashTableCount(
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
				  int nNodes,
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
				  unsigned int elmCnt,
				  unsigned int elmSz,
				  unsigned int blkSz,
                              	  AlcErrno *dstErr);
extern AlcErrno			AlcVectorFree(
				  AlcVector *vec);
extern AlcErrno			AlcVectorExtend(
				  AlcVector *vec,
				  unsigned int elmCnt);
extern void			*AlcVectorItemGet(
				  AlcVector *vec,
				  unsigned int idx);
extern void			*AlcVectorExtendAndGet(
				  AlcVector *vec,
				  unsigned int idx);
extern unsigned int		AlcVectorCount(
				  AlcVector *vec);
extern void			AlcVectorSetArray1D(
				  AlcVector *vec,
				  int fIdx,
				  int lIdx,
				  void *aM);
extern void			*AlcVectorToArray1D(
				  AlcVector *vec,
				  int fIdx,
				  int lIdx,
				  AlcErrno *dstErr);
void				**AlcVectorToArray2D(
				  AlcVector *vec,
				  int fIdx,
				  int lIdx,
				  int nR,
				  int nC,
				  AlcErrno *dstErr);
#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* ALCPROTO_H */
