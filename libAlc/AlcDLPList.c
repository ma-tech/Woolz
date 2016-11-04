#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlcDLPList_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlc/AlcDLPList.c
* \author       Bill Hill
* \date         November 1999
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2012],
* The University Court of the University of Edinburgh,
* Old College, Edinburgh, UK.
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
* \brief        A general purpose doubly linked circular list of pointers.
* \note		This code has been derived from the hguDlpList but
*		has been stripped down for efficiency.
* \ingroup	AlcDLPList
*/
#include <stdio.h>
#include <stdlib.h>
#include <Alc.h>

static void			AlcDLPListSortPv(
				  AlcDLPItem *low,
				  AlcDLPItem *high,
				  int lowIdx,
				  int highIdx,
				  int (*entryCompFn)(void *, void *));


/*!
* \return	List data structure, or NULL on error.
* \ingroup	AlcDLPList
* \brief	Creates a list data structure which is required by all
*               the other AlcDLPList functions.
* \param	dstErr 			Destination pointer for error
*					code, may be NULL.
*/
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

/*!
* \return	List item structure, or	NULL on error.
* \ingroup	AlcDLPList
* \brief	Creates a list item data structure for building into
*		a AlcDLPList list.
* \param	entry 			New list entry.
* \param	freeFn 			Function that will be called
*					(if not NULL) to free the entry.
* \param	dstErr 			Destination pointer for error
*					code, may be NULL.
*/
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

/*!
* \return	Error code.
* \ingroup	AlcDLPList
* \brief	Free's the given list data structure and any list items.
* \param	list 		 	The list data structure.
*/
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

/* function:     AlcDLPListEntryInsert    */
/*! 
* \ingroup      AlcDLPList
* \brief        Inserts the given entry into the list before the given
*		item.
* \param	list 			The list data structure.
* \param	insBefore:		Given item that entry is to
*					be inserted before, if NULL
*					then the item is inserted at
*					the head of the list.
* \param	entry:			New list entry.
* \param	freeFn 			Function that will be called
*					(if not NULL) to free the entry
*					should the item be deleted.
*/
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

/*!
* \return	Error code.
* \ingroup	AlcDLPList
* \brief	Appends the given entry into the list after the given
*		item.
* \param	list 			The list data structure.
* \param	appAfter 		 Given item that entry is to
*					be inserted after, if NULL
*					then the item is appended at
*					the tail of the list.
* \param	entry 			New list entry.
* \param	freeFn 			 Function that will be called
*					(if not NULL) to free the entry
*					should the item be deleted.
*/
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

/*!
* \return	Next list item after the item removed, NULL on error or
*		if last item removed.
* \ingroup	AlcDLPList
* \brief	Removes the item from the list, but does not free the
*		item unless the freeItem flag is set.
* \param	list 			The list data structure.
* \param	item 			Item to be removed.
* \param	freeItem 		Free item if non-zero.
* \param	dstErr 			Destination pointer for error
*					code, may be NULL.
*/
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
      if(item == item->next)
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
      (void )AlcDLPItemFree(item);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(nextItem);
}

/*!
* \return	Error code.
* \ingroup	AlcDLPList
* \brief	Inserts a new item into the list before the given
*		item.
* \param	list 			The list data structure.
* \param	insBefore 		Given item that entry is to
*					be inserted before, if NULL
*					then the item is inserted at
*					the head of the list.
* \param	newItem 		New item to insert.
*/
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

/*!
* \return	Error code.
* \ingroup	AlcDLPList
* \brief	Appends a new item into the list after the given
*		item.
* \param	list 			The list data structure.
* \param	appAfter 		Given item that entry is to
*					be inserted after, if NULL
*					then the item is appended at
*					the tail of the list.
* \param	newItem 		New item to insert.
*/
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

/*!
* \return	Error code.
* \ingroup	AlcDLPList
* \brief	Free's the list item which has already been removed
*		from the list.
* \param	item 		 	Item to be deleted.
*/
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

/*!
* \return	Error code.
* \ingroup	AlcDLPList
* \brief	Sorts the entire list using the given entry comparison
*		function.
* \param	list:			The list data structure.
* \param	entryCompFn 		Given entry
*					comparison function which must
*					return  an integer less than,
*					equal to, or greater than zero
*					to indicate if the first entry
*					is to be  considered  less
*					than, equal to, or greater than
*					the second.
*/
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

/*!
* \return	Error code.
* \ingroup	AlcDLPList
* \brief	Private recursive quick sort function to sorts the
*		list using the given entry comparison function.
*		The values of the indicies must be such that lowIdx is
*		less than highIdx and their difference is equal to one
*		more than the number of items between them (excluding
*		the items themselves).
* \param	low 			Item nearest to tail in section
*					to be sorted.
* \param	high 			Item nearest to head in section
*					to be sorted.
* \param	lowIdx 			Index for item nearest to tail.
* \param	highIdx 		Index for item nearest to head.
* \param	entryCompFn 		Given entry
*					comparison function which must
*					return  an integer less than,
*					equal to, or greater than zero
*					to indicate if the first entry
*					is to be  considered  less
*					than, equal to, or greater than
*					the second.
*/
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

/*!
* \return	Number of items in list. This is always >= 0.
* \ingroup	AlcDLPList
* \brief	Returns the number of items in the list.
* \param	list 			The list data structure.
* \param	dstErr 			Destination pointer for error
*					code, may be NULL.
*/
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

/*!
* \return	Last item.
* \ingroup	AlcDLPList
* \brief	Iterates the given function through the list, starting
*		with the given item. The iteration may proceed towards
*		either the head or tail of the list. The iterated
*		function must take the form
* \verbatim
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
* \endverbatim
*		Where the data parameter is the data supplied to this
*		(AlcDLPListIterate) function.
*		The iteration continues until the itterated function
*		returns zero.
* \param	list 			The list data structure.
* \param	item 			First item of iteration which
*					may be NULL to indicate list
*					head (dir == _TO_TAIL) or tail
*					(dir == _TO_HEAD).
* \param	dir 		 	Iteration direction, either
*					towards the head or the tail.
* \param	iterFn 			Function to be iterated, see
*					example above.
* \param	iterData 		Data supplied to the iterated
*					function.
* \param	dstErr 			Destination pointer for error
*					code, may be NULL.
*/
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
