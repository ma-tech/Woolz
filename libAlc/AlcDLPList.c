#pragma ident "MRC HGU $Id"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        AlcDLPList.c
* Date:         November 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      A general purpose doubly linked circular list of
*		pointers.
*		This code has been derived from the hguDlpList but
*		has been stripped down for efficiency.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <Alc.h>

static void	AlcDLPListSortPv(AlcDLPItem *low, AlcDLPItem *high,
				 int lowIdx, int highIdx,
				 int (*entryCompFn)(void *, void *));


/************************************************************************
* Function:	AlcDLPListNew
* Returns:	AlcDLPList *:		List data structure, or	NULL on
*					error.
* Purpose:	Creates a list data structure which is required by all
*		the other AlcDLPList functions.
* Global refs:	-
* Parameters:	AlcErrno *dstErr:	Destination pointer for error
*					code, may be NULL.
************************************************************************/
AlcDLPList	*AlcDLPListNew(AlcErrno *dstErr)
{
  AlcDLPList	*list = NULL;
  AlcErrno	errNum = ALC_ER_NONE;

  if((list = (AlcDLPList *)AlcCalloc(1, sizeof(AlcDLPList))) == NULL)
  {
    errNum = ALC_ER_ALLOC;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(list);
}

/************************************************************************
* Function:	AlcDLPItemNew
* Returns:	AlcDLPItem *:		List item structure, or	NULL on
*					error.
* Purpose:	Creates a list item data structure for building into
*		a AlcDLPList list.
* Global refs:	-
* Parameters:	void *entry:		New list entry.
*		void (*freeFn)(void *)): Function that will be called
*					(if not NULL) to free the entry.
*		AlcErrno *dstErr:	Destination pointer for error
*					code, may be NULL.
************************************************************************/
AlcDLPItem	*AlcDLPItemNew(void *entry, void (*freeFn)(void *),
			       AlcErrno *dstErr)
{
  AlcDLPItem 	*newItem;
  AlcErrno	errNum = ALC_ER_NONE;

  if((newItem = (AlcDLPItem *)AlcCalloc(1, sizeof(AlcDLPItem))) == NULL)
  {
    errNum = ALC_ER_ALLOC;
  }
  else
  {
    newItem->freeFn = freeFn;
    newItem->entry = entry;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(newItem);
}

/************************************************************************
* Function:	AlcDLPListFree
* Returns:	AlcErrno:		Error code.
* Purpose:	Free's the given list data structure and any list
*		items.
* Global refs:	-
* Parameters:	AlcDLPList *list:	The list data structure.
************************************************************************/
AlcErrno	AlcDLPListFree(AlcDLPList *list)
{
  AlcDLPItem 	*head,
  		*item,
		*fItem;
  AlcErrno	errNum = ALC_ER_NONE;

  if(list == NULL)
  {
    errNum =  ALC_ER_NULLPTR;
  }
  else
  {
    if((head = list->head) != NULL)
    {
      item = head;
      do
      {
        fItem = item;
	item = item->next;
	if(fItem->freeFn && fItem->entry)
	{
	  (*(fItem->freeFn))(fItem->entry);
	}
	AlcFree(fItem);
      }
      while(item != head);
    }
    AlcFree(list);
  }
  return(errNum);
}

/************************************************************************
* Function:	AlcDLPListEntryInsert
* Returns:	AlcErrno:		Error code.
* Purpose:	Inserts the given entry into the list before the given
*		item.
* Global refs:	-
* Parameters:	AlcDLPList *list:	The list data structure.
*		AlcDLPItem *insBefore:	Given item that entry is to
*					be inserted before, if NULL
*					then the item is inserted at
*					the head of the list.
*		void *entry:		New list entry.
*		void (*freeFn)(void *)): Function that will be called
*					(if not NULL) to free the entry
*					should the item be deleted.
************************************************************************/
AlcErrno	AlcDLPListEntryInsert(AlcDLPList *list, AlcDLPItem *insBefore,
				      void *entry, void (*freeFn)(void *))
{
  AlcDLPItem 	*item,
		*newItem;
  AlcErrno	errNum = ALC_ER_NONE;

  if(list == NULL)
  {
    errNum = ALC_ER_NULLPTR;
  }
  else
  {
    if((newItem = AlcDLPItemNew(entry, freeFn, &errNum)) != NULL)
    {
      if(list->head == NULL)
      {
        item = list->head = newItem;
	item->prev = item;
	item->next = item;
      }
      else
      {
        if(insBefore)
	{
	  item = insBefore;
	}
	else
	{
	  item = list->head;
	  list->head = newItem;
	}
	newItem->next = item;
	newItem->prev = item->prev;
	item->prev->next = newItem;
	item->prev = newItem;
      }
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	AlcDLPListEntryAppend
* Returns:	AlcErrno:		Error code.
* Purpose:	Appends the given entry into the list after the given
*		item.
* Global refs:	-
* Parameters:	AlcDLPList *list:	The list data structure.
*		AlcDLPItem *appAfter: 	Given item that entry is to
*					be inserted after, if NULL
*					then the item is appended at
*					the tail of the list.
*		void *entry:		New list entry.
*		void (*freeFn)(void *)): Function that will be called
*					(if not NULL) to free the entry
*					should the item be deleted.
************************************************************************/
AlcErrno	AlcDLPListEntryAppend(AlcDLPList *list, AlcDLPItem *appAfter,
				      void *entry, void (*freeFn)(void *))
{
  AlcDLPItem 	*item,
  		*newItem;
  AlcErrno	errNum = ALC_ER_NONE;

  if(list == NULL)
  {
    errNum = ALC_ER_NULLPTR;
  }
  else
  {
    if((newItem = AlcDLPItemNew(entry, freeFn, &errNum)) != NULL)
    {
      if(list->head == NULL)
      {
        item = list->head = newItem;
	item->prev = item;
	item->next = item;
      }
      else
      {
        if(appAfter)
	{
	  item = appAfter;
	}
	else
	{
	  item = list->head->prev;
	}
	newItem->prev = item;
	newItem->next = item->next;
	item->next->prev = newItem;
	item->next = newItem;
      }
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	AlcDLPItemUnlink
* Returns:	AlcDLPItem *:		Next list item after the item
*					removed, NULL on error or if
*					last item removed.
* Purpose:	Removes the item from the list, but does not free the
*		item unless the freeItem flag is set.
* Global refs:	-
* Parameters:	AlcDLPList *list:	The list data structure.
*		AlcDLPItem *item:	Item to be removed.
*		int freeItem;		Free item if non-zero.
*		AlcErrno *dstErr:	Destination pointer for error
*					code, may be NULL.
************************************************************************/
AlcDLPItem	*AlcDLPItemUnlink(AlcDLPList *list, AlcDLPItem *item,
				  int freeItem, AlcErrno *dstErr)
{
  AlcDLPItem 	*nextItem = NULL;
  AlcErrno	errNum = ALC_ER_NONE;

  if((list == NULL) || (item == NULL))
  {
    errNum = ALC_ER_NULLPTR;
  }
  else
  {
    if(item == list->head)
    {
      if(item = item->next)
      {
	nextItem = NULL;
        list->head = NULL;
      }
      else
      {
	nextItem = item->next;
        list->head = nextItem;
	item->prev->next = nextItem;
	nextItem->prev = item->prev;
      }
    }
    else
    {
      nextItem = item->next;
      item->prev->next = nextItem;
      nextItem->prev = item->prev;
    }
    if(freeItem)
    {
      AlcDLPItemFree(item);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(nextItem);
}

/************************************************************************
* Function:	AlcDLPItemInsert
* Returns:	AlcErrno:		Error code.
* Purpose:	Inserts a new item into the list before the given
*		item.
* Global refs:	-
* Parameters:	AlcDLPList *list:	The list data structure.
*		AlcDLPItem *insBefore: Given item that entry is to
*					be inserted before, if NULL
*					then the item is inserted at
*					the head of the list.
*		AlcDLPItem *newItem: New item to insert.
************************************************************************/
AlcErrno	AlcDLPItemInsert(AlcDLPList *list,
				     AlcDLPItem *insBefore,
				     AlcDLPItem *newItem)
{
  AlcDLPItem 	*item;
  AlcErrno	errNum = ALC_ER_NONE;

  if((list == NULL) || (newItem == NULL))
  {
    errNum = ALC_ER_NULLPTR;
  }
  else
  {
    if(list->head == NULL)
    {
      item = list->head = newItem;
      item->prev = item;
      item->next = item;
    }
    else
    {
      if(insBefore)
      {
	item = insBefore;
      }
      else
      {
	item = list->head;
	list->head = newItem;
      }
      newItem->next = item;
      newItem->prev = item->prev;
      item->prev->next = newItem;
      item->prev = newItem;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	AlcDLPItemAppend
* Returns:	AlcErrno:		Error code.
* Purpose:	Appends a new item into the list after the given
*		item.
* Global refs:	-
* Parameters:	AlcDLPList *list:	The list data structure.
*		AlcDLPItem *appAfter: Given item that entry is to
*					be inserted after, if NULL
*					then the item is appended at
*					the tail of the list.
*		AlcDLPItem *newItem: New item to insert.
************************************************************************/
AlcErrno	AlcDLPItemAppend(AlcDLPList *list,
				 AlcDLPItem *appAfter,
				 AlcDLPItem *newItem)
{
  AlcDLPItem 	*item;
  AlcErrno	errNum = ALC_ER_NONE;

  if((list == NULL) || (newItem == NULL))
  {
    errNum = ALC_ER_NULLPTR;
  }
  else
  {
    if(list->head == NULL)
    {
      item = list->head = newItem;
      item->prev = item;
      item->next = item;
    }
    else
    {
      if(appAfter)
      {
	item = appAfter;
      }
      else
      {
	item = list->head->prev;
      }
      newItem->prev = item;
      newItem->next = item->next;
      item->next->prev = newItem;
      item->next = newItem;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	AlcDLPItemFree
* Returns:	AlcErrno:		Error code.
* Purpose:	Free's the list item which has already been removed
*		from the list.
* Global refs:	-
* Parameters:	AlcDLPItem *item:	Item to be deleted.
************************************************************************/
AlcErrno	AlcDLPItemFree(AlcDLPItem *item)
{
  AlcErrno	errNum = ALC_ER_NONE;

  if(item == NULL)
  {
    errNum = ALC_ER_NULLPTR;
  }
  else
  {
    if(item->freeFn && item->entry)
    {
      (*(item->freeFn))(item->entry);
    }
    AlcFree(item);
  }
  return(errNum);
}

/************************************************************************
* Function:	AlcDLPListSort
* Returns:	AlcErrno:		Error code.
* Purpose:	Sorts the entire list using the given entry comparison
*		function.
* Global refs:	-
* Parameters:	AlcDLPList *list:	The list data structure.
*		int (*entryCompFn)(void *, void *): Given entry
*					comparison function which must
*					return  an integer less than,
*					equal to, or greater than zero
*					to indicate if the first entry
*					is to be  considered  less
*					than, equal to, or greater than
*					the second.
************************************************************************/
AlcErrno	AlcDLPListSort(AlcDLPList *list,
			       int (*entryCompFn)(void *, void *))
{
  int		itemCnt;
  AlcErrno	errNum = ALC_ER_NONE;

  if((list == NULL) || (entryCompFn == NULL))
  {
    errNum = ALC_ER_NULLPTR;
  }
  else if(list->head && ((itemCnt = AlcDLPListCount(list, NULL)) > 0))
  {
    AlcDLPListSortPv(list->head->prev, list->head, 0, itemCnt, entryCompFn);
  }
  return(errNum);
}

/************************************************************************
* Function:	AlcDLPListSortPv
* Returns:	AlcErrno:		Error code.
* Purpose:	Private recursive quick sort function to sorts the
*		list using the given entry comparison function.
*		The values of the indicies must be such that lowIdx is
*		less than highIdx and their difference is equal to one
*		more than the number of items between them (excluding
*		the items themselves).
* Global refs:	-
* Parameters:	AlcDLPItem *low:	Item nearest to tail in section
*					to be sorted.
*		AlcDLPItem *high:	Item nearest to head in section
*					to be sorted.
*		int lowIdx:		Index for item nearest to tail.
*		int highIdx:		Index for item nearest to head.
*		int (*entryCompFn)(void *, void *): Given entry
*					comparison function which must
*					return  an integer less than,
*					equal to, or greater than zero
*					to indicate if the first entry
*					is to be  considered  less
*					than, equal to, or greater than
*					the second.
************************************************************************/
static void	AlcDLPListSortPv(AlcDLPItem *low, AlcDLPItem *high,
				 int lowIdx, int highIdx,
				 int (*entryCompFn)(void *, void *))
{

  AlcDLPItem *pivot,
  		 *base;
  int		pivotIdx,
  		baseIdx;
  void		*tP0;

  baseIdx = lowIdx;
  base = low;    				  /* Remember the local tail */
  pivotIdx = highIdx;
  pivot = high;    				  /* Partition off the pivot */
  --highIdx;
  high = high->next;
  do
  {
    while((lowIdx < highIdx) &&
          ((*entryCompFn)(low->entry, pivot->entry) <= 0))
    {
      ++lowIdx;
      low = low->prev;
    }
    while((lowIdx < highIdx) &&
          ((*entryCompFn)(high->entry, pivot->entry) >= 0))
    {
      --highIdx;
      high = high->next;
    }
    if(lowIdx < highIdx)
    {
      tP0 = high->entry; 				 /* Exchange entries */
      high->entry = low->entry;
      low->entry = tP0;
    }
  } while(lowIdx < highIdx);
  if((lowIdx < pivotIdx) &&
     ((*entryCompFn)(low->entry, pivot->entry) > 0))
  {
    tP0 = pivot->entry;				         /* Exchange entries */
    pivot->entry = low->entry;
    low->entry = tP0;
  }
  ++lowIdx;
  low = low->prev;
  if((highIdx - baseIdx) < (pivotIdx - lowIdx))
  {
    if(lowIdx < pivotIdx)
    {
      AlcDLPListSortPv(low, pivot, lowIdx, pivotIdx, entryCompFn);
    }
    if(baseIdx < highIdx)
    {
      AlcDLPListSortPv(base, high, baseIdx, highIdx, entryCompFn);
    }
  }
  else
  {
    if(baseIdx < highIdx)
    {
      AlcDLPListSortPv(base, high, baseIdx, highIdx, entryCompFn);
    }
    if(lowIdx < pivotIdx)
    {
      AlcDLPListSortPv(low, pivot, lowIdx, pivotIdx, entryCompFn);
    }
  }
}

/************************************************************************
* Function:	AlcDLPListCount
* Returns:	int:			Number of items in list. This
*					is always >= 0.
* Purpose:	Returns the number of items in the list.
* Global refs:	-
* Parameters:	AlcDLPList *list:	The list data structure.
*		AlcErrno *dstErr:	Destination pointer for error
*					code, may be NULL.
************************************************************************/
int		AlcDLPListCount(AlcDLPList *list, AlcErrno *dstErr)
{
  int		cnt = 0;
  AlcDLPItem *item,
  		*head;
  AlcErrno	errNum = ALC_ER_NONE;

  if(list == NULL)
  {
    errNum = ALC_ER_NULLPTR;
  }
  else if((head = list->head) != NULL)
  {
    item = head;
    do
    {
      ++cnt;
      item = item->next;
    } while(item != head);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cnt);
}

/************************************************************************
* Function:	AlcDLPListIterate
* Returns:	AlcDLPItem:		Last item.
* Purpose:	Iterates the given function through the list, starting
*		with the given item. The iteration may proceed towards
*		either the head or tail of the list. The iterated
*		function must take the form:
*		  int func(AlcDLPList *, AlcDLPItem *, void *)
*		As in the example:
*		  int MyItemCount(AlcDLPList *list,
*				 AlcDLPItem *item,
*				 void *myData)
*		  {
*		    int		*count;
*
*		    if(list && item)
*		    {
*		      count = (int *)myData;
*		      ++*count;
*		    }
*		    return(list->head != item->next);
*		  }
*		Where the data parameter is the data supplied to this
*		(AlcDLPListIterate) function.
*		The iteration continues until the itterated function
*		returns zero.
* Global refs:	-
* Parameters:	AlcDLPList *list:	The list data structure.
*		AlcDLPItem *item:	First item of iteration which
*					may be NULL to indicate list
*					head (dir == _TO_TAIL) or tail
*					(dir == _TO_HEAD).
*		AlcDirection dir: 	Iteration direction, either
*					towards the head or the tail.
*		int (*iterFn)():	Function to be iterated, see
*					example above.
*		void *iterData:		Data supplied to the iterated
*					function.
*		AlcErrno *dstErr:	Destination pointer for error
*					code, may be NULL.
************************************************************************/
AlcDLPItem	*AlcDLPListIterate(AlcDLPList *list, AlcDLPItem *item,
				   AlcDirection dir,
				   int (*iterFn)(AlcDLPList *,
				   		 AlcDLPItem *, void *),
				   void *iterData, AlcErrno *dstErr)
{
  int		iterate = 1;
  AlcDLPItem 	*lastItem = NULL;
  AlcErrno	errNum = ALC_ER_NONE;

  if((list == NULL) || (iterFn == NULL))
  {
    errNum = ALC_ER_NULLPTR;
  }
  else
  {
    if((item == NULL) && (list->head))
    {
      switch(dir)
      {
        case ALC_DIRECTION_FWD:
	  item = list->head->prev;
	  break;
	case ALC_DIRECTION_REV:
	  item = list->head;
	  break;
        default:
	  errNum = ALC_ER_PARAM;
	  break;
      }
    }
    while(item && iterate)
    {
      lastItem = item;
      if(dir == ALC_DIRECTION_FWD)
      {
        item = lastItem->prev;
      }
      else
      {
        item = lastItem->next;
      }
      iterate = (*iterFn)(list, lastItem, iterData);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(lastItem);
}
