#ifndef HGUDLPLIST_H
#define HGUDLPLIST_H
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        HGUDlpList.h
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Data structures and functions for doubly linked lists
*		of pointers.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/

#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef enum					 /* List traversal direction */
{
  HGU_DLPLIST_DIR_TOHEAD,
  HGU_DLPLIST_DIR_TOTAIL
} HGUDlpListDirection;

typedef enum				  /* State of list locking mechanism */
{
  HGU_DLPLIST_STATE_EMPTY 	= (0),
  HGU_DLPLIST_STATE_ERROR 	= (0),
  HGU_DLPLIST_STATE_UNLOCK 	= (1),
  HGU_DLPLIST_STATE_LOCK	= (1<<1),
  HGU_DLPLIST_STATE_CREATE	= (1<<2),
  HGU_DLPLIST_STATE_DESTROY	= (1<<3)
} HGUDlpListState;

#ifdef HGUDLPLIST_C
typedef struct _HGUDlpListItem			  /* Doubly linked list item */
{
  void          (*freeFn)(void *);
  void          *entry;
  struct _HGUDlpListItem *next;
  struct _HGUDlpListItem *prev;
} HGUDlpListItem;
 
typedef struct _HGUDlpList		   /* Doubly linked list of pointers */
{
  HGUDlpListState (*lockFn)(void *, HGUDlpListState);
  void          *lockData;
  int           itemCount;
  HGUDlpListItem *head;
  HGUDlpListItem *tail;
 
} HGUDlpList;
#else /* ! HGUDLPLIST_C */
typedef void HGUDlpListItem;    /* Opaque handle for doubly linked list item */
typedef void HGUDlpList; /* Opaque handle for doubly linked list of pointers */

extern int	HGUDlpListSort(HGUDlpList *list,
			       int (*entryCompFn)(void *, void *)),
		HGUDlpListItemIsHead(HGUDlpList *list, HGUDlpListItem *item),
		HGUDlpListItemIsTail(HGUDlpList *list, HGUDlpListItem *item),
		HGUDlpListCount(HGUDlpList *list),
		HGUDlpListOffset(HGUDlpList *list, HGUDlpListItem *item,
		 		 HGUDlpListDirection dir);
extern void	HGUDlpListDestroy(HGUDlpList *list);
extern void 	*HGUDlpListEntryGet(HGUDlpList *list, HGUDlpListItem *item),
		*HGUDlpListEntrySet(HGUDlpList *list, HGUDlpListItem *item,
				    void *entry);
extern HGUDlpList *HGUDlpListCreate(HGUDlpListState (*lockFn)(void *,
						              HGUDlpListState)),
		*HGUDlpListDup(HGUDlpList *list);
extern HGUDlpListItem *HGUDlpListInsert(HGUDlpList *list,
					HGUDlpListItem *before, void *entry,
					void (*freeFn)(void *)),
		*HGUDlpListAppend(HGUDlpList *list,
				  HGUDlpListItem *after, void *entry,
				  void (*freeFn)(void *)),
		*HGUDlpListExchange(HGUDlpList *list,
				    HGUDlpListItem *item0,
				    HGUDlpListItem *item1),
		*HGUDlpListDeleteAll(HGUDlpList *list),
		*HGUDlpListDelete(HGUDlpList *list, HGUDlpListItem *item),
		*HGUDlpListRemoveAll(HGUDlpList *list),
		*HGUDlpListRemove(HGUDlpList *list,
				  HGUDlpListItem *item),
		*HGUDlpListIterate(HGUDlpList *list, HGUDlpListItem *item,
				   HGUDlpListDirection dir,
				   int (*iterFn)(HGUDlpList *, HGUDlpListItem *,
						 void *),
				   void *iterData),
		*HGUDlpListHead(HGUDlpList *list),
		*HGUDlpListTail(HGUDlpList *list),
		*HGUDlpListNext(HGUDlpList *list, 
				HGUDlpListItem *item),
		*HGUDlpListPrev(HGUDlpList *list,
				HGUDlpListItem *item),
		*HGUDlpListNth(HGUDlpList *list, HGUDlpListItem *item,
			       HGUDlpListDirection dir, int num);
#endif /* HGUDLPLIST_C */

#ifdef  __cplusplus
}
#endif /* __cplusplus */

#endif /* HGUDLPLIST_H */
