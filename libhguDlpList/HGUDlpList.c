#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _HGUDlpList_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libhguDlpList/HGUDlpList.c
* \author       Bill Hill
* \date         March 1999
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
* \brief	Data structures and functions for doubly linked lists
* 		of pointers.
* \ingroup	hguDlpList
*/

#define HGUDLPLIST_C

#include <stdio.h>
#include <stdlib.h>
/*#ifndef __ppc
#include <malloc.h>
#endif*/
#include <assert.h>
#include <HGUDlpList.h>

static HGUDlpListItem 		*HGUDlpListPvItemCreate(
				  void *, void (*)(void *));
static int			HGUDlpListPvUnlock(
				  HGUDlpList *);
static int			HGUDlpListPvLock(
				  HGUDlpList *);
static void			HGUDlpListPvDeleteAll(
				  HGUDlpList *);
static void			HGUDlpListPvSort(
				  HGUDlpListItem *,
				  HGUDlpListItem *,
				  int,
				  int,
				  int (*)(void *, void *));
static HGUDlpListItem 		*HGUDlpListPvDelRmItem(
				  HGUDlpList *,
				  HGUDlpListItem *,
				  int);

/*!
* \return	List data structure, or NULL on error.
* \ingroup	hguDlpList
* \brief	Creates a list data structure which is required by all
* 		the other HGUDlpList functions.
* \param	lockFn			The lock function to be used when
* 					accessing the list. If NULL
* 					then no locking is used.
*/
HGUDlpList	*HGUDlpListCreate(HGUDlpListState (*lockFn)(void *,
							    HGUDlpListState))
{
  HGUDlpList *list;

  list = (HGUDlpList *)malloc(sizeof(HGUDlpList));
  if(list)
  {
    list->lockFn = lockFn;
    list->lockData = NULL;
    list->itemCount = 0;
    list->head = NULL;
    list->tail = NULL;
    if(list->lockFn)
    {
      if((*list->lockFn)((void *)&(list->lockData),
      			 HGU_DLPLIST_STATE_CREATE | HGU_DLPLIST_STATE_LOCK) !=
         HGU_DLPLIST_STATE_LOCK)
      {
        free(list);
	list = NULL;
      }
    }
  }
  return(list);
}

/*!
* \return	List data structure, or NULL on error.
* \ingroup	hguDlpList
* \brief	Duplicates a list data structure, but NOT its items.
* 		Ie head and tail are both NULL and the item count is
* 		zero, only the lock function is common.
* \param	givenList
*/
HGUDlpList	*HGUDlpListDup(HGUDlpList *givenList)
{
  HGUDlpList *newList = NULL;

  if(givenList)
    newList = HGUDlpListCreate(givenList->lockFn);
  return(newList);
}

/*!
* \return	void
* \ingroup	hguDlpList
* \brief	Destroys the given list list data structure and any
* 		list items.
* \param	list			The list data structure.
*/
void		HGUDlpListDestroy(HGUDlpList *list)
{
  if(list)
  {
    (void )HGUDlpListPvLock(list);			      /* Lock list [ */
    (void )HGUDlpListPvDeleteAll(list);
    if(list->lockFn)
      (void )(*list->lockFn)(list->lockData,
      			     HGU_DLPLIST_STATE_DESTROY);   /* Destroy lock ] */
    free(list);
  }
}

/*!
* \return
* \ingroup	hguDlpList
* \brief	Inserts the given entry into the list before the given
*		item.
* \param	list			The list data structure.
* \param	insBefore		Given item that entry is to
* 					be inserted before, if NULL then
* 					the item is inserted at the head
* 					of the list.
* \param	entry			New list entry.
* \param	freeFn			Function that will be called
* 					(if not NULL) to free the entry
* 					should the item be deleted.
*/
HGUDlpListItem	*HGUDlpListInsert(HGUDlpList *list,
				  HGUDlpListItem *insBefore, void *entry,
				  void (*freeFn)(void *))
{
  HGUDlpListItem *tItem0,
  		*newItem = NULL;

  if(list)
  {
    if(HGUDlpListPvLock(list))				      /* Lock list [ */
    {
      if((newItem = HGUDlpListPvItemCreate(entry, freeFn)) != NULL)
      {
	if(list->itemCount == 0)     /* Item is 1st so at both head and tail */
	{
	  assert((list->head == NULL) && (list->tail == NULL));
	  list->head = newItem;
	  list->tail = newItem;
	  list->itemCount = 1;
	}
	else if((list->itemCount == 1) ||
	        (insBefore == NULL))			   /* Insert at head */
	{
	  assert((insBefore == NULL) ||
	         ((unsigned long )(list->head) == 
		  (unsigned long )(list->tail)));
	  tItem0 = list->head;
	  list->head = newItem;
	  tItem0->prev = newItem;
	  newItem->next = tItem0;
	  assert((insBefore == NULL) || (tItem0->next == NULL));
	  ++(list->itemCount);
	}
	else if(list->itemCount > 1)         /* Item anywhere except at tail */
	{
	  tItem0 = insBefore->prev;
	  newItem->prev = tItem0;
	  newItem->next = insBefore;
	  insBefore->prev = newItem;
	  if(tItem0)
	    tItem0->next = newItem;
	  else
	  {
	    assert((unsigned long )insBefore == (unsigned long)(list->head));
	    list->head = newItem;
	  }
	  ++(list->itemCount);
	}
      }
      (void )HGUDlpListPvUnlock(list);			    /* Unlock list ] */
    }
  }
  return(newItem);
}

/*!
* \return	The appended list item, or NULL on error.
* \ingroup	hguDlpList
* \brief	Appends the given entry into the list after the given
*		item.
* \param	list			The list data structure.
* \param	appAfter		Given item that entry is to
* 					be inserted after, if NULL
* 					then the item is appended at
* 					the tail of the list.
* \param	entry			New list entry.
* \param	freeFn			Function that will be called
* 					(if not NULL) to free the entry
* 					should the item be deleted.
*/
HGUDlpListItem	*HGUDlpListAppend(HGUDlpList *list,
				  HGUDlpListItem *appAfter, void *entry,
				  void (*freeFn)(void *))
{
  HGUDlpListItem *tItem0,
  		*newItem = NULL;

  if(list)
  {
    if(HGUDlpListPvLock(list))				      /* Lock list [ */
    {
      if((newItem = HGUDlpListPvItemCreate(entry, freeFn)) != NULL)
      {
	if(list->itemCount == 0)     /* Item is 1st so at both head and tail */
	{
	  assert((list->head == NULL) && (list->tail == NULL));
	  list->head = newItem;
	  list->tail = newItem;
	  list->itemCount = 1;
	}
	else if((list->itemCount == 1) ||
	        (appAfter == NULL))			   /* Append at tail */
	{
	  assert((appAfter == NULL) ||
	         ((unsigned long )(list->head) == 
		  (unsigned long )(list->tail)));
	  tItem0 = list->tail;
	  list->tail = newItem;
	  tItem0->next = newItem;
	  newItem->prev = tItem0;
	  assert((appAfter == NULL) || (tItem0->prev == NULL));
	  ++(list->itemCount);
	}
	else if(list->itemCount > 1)         /* Item anywhere except at head */
	{
	  tItem0 = appAfter->next;
	  newItem->next = tItem0;
	  newItem->prev = appAfter;
	  appAfter->next = newItem;
	  if(tItem0)
	    tItem0->prev = newItem;
	  else
	  {
	    assert((unsigned long )appAfter == (unsigned long)(list->tail));
	    list->tail = newItem;
	  }
	  ++(list->itemCount);
	}
      }
      (void )HGUDlpListPvUnlock(list);			    /* Unlock list ] */
    }
  }
  return(newItem);
}

/*!
* \return	The first of the two list items (ie  first in the
* 		parameter list) or NULL on error.
* \ingroup	hguDlpList
* \brief	Exchanges the two given list item entries and not the
*		items, so that head is still head and tail is still
*		tail.
* \param	list			The list data structure.
* \param	item0			First item (maybe before, after
* 					or even the same as the second
* 					item.
* \param	item1			The second item.
*/
HGUDlpListItem	*HGUDlpListExchange(HGUDlpList *list,
				    HGUDlpListItem *item0,
				    HGUDlpListItem *item1)
{
  void		*tEntry0;
  HGUDlpListItem *rtnItem = NULL;

  if(list && item0 && item1)
  {
    if(HGUDlpListPvLock(list))				      /* Lock list [ */
    {
      tEntry0 = item0->entry;
      item0->entry = item1->entry;
      item1->entry = tEntry0;
      rtnItem = item0;
      (void )HGUDlpListPvUnlock(list);			    /* Unlock list ] */
    }
  }
  return(rtnItem);
}

/*!
* \return	New list head which is NULL if all the items were
* 		deleted ok.
* \ingroup	hguDlpList
* \brief	Deletes all list items from the head on down to and
*		including the tail. Where delete implies both the
*		removal of items from the list and freeing the entries
*		using the item's free functions (unless either the free
*		function or the entry is NULL).
* \param	list			The list data structure.
*/
HGUDlpListItem	*HGUDlpListDeleteAll(HGUDlpList *list)
{
  if(list)
  {
    if(HGUDlpListPvLock(list))				      /* Lock list [ */
    {
      HGUDlpListPvDeleteAll(list);
      (void )HGUDlpListPvUnlock(list);			    /* Unlock list ] */
    }
  }
  return(list->head);
}

/*!
* \return	Next list item after the item deleted, may be NULL either
* 		on error or if the list item was the last in the list.
* \ingroup	hguDlpList
* \brief	Deletes the given list item from the list with the
*		given list. Where delete implies both the removal of
*		an item, the freeing of the item AND (when neither the
*		free function or the entry are NULL) it's entry too.
* \param	list			The list data structure.
* \param	item			Item to be deleted.
*/
HGUDlpListItem	*HGUDlpListDelete(HGUDlpList *list,
				  HGUDlpListItem *item)
{
  const int	delete = 1;
  HGUDlpListItem *nextItem;

  nextItem = HGUDlpListPvDelRmItem(list, item, delete);
  return(nextItem);
}

/*!
* \return	Next list item after the item removed, may be NULL
* 		either on error or if the list item was the last in the list.
* \ingroup	hguDlpList
* \brief	Removes the item from the list withe the given list.
*		Where remove implies the removal of the item from the
*		list and the freeing of the item EXCEPT for its entry.
* \param	list			The list data structure.
* \param	item			Item to be deleted.
*/
HGUDlpListItem	*HGUDlpListRemove(HGUDlpList *list,
				  HGUDlpListItem *item)
{
  const int	delete = 0;
  HGUDlpListItem *nextItem;

  nextItem = HGUDlpListPvDelRmItem(list, item, delete);
  return(nextItem);
}

/*!
* \return	Number of items in list or zero for an empty list or on
* 		error.
* \ingroup	hguDlpList
* \brief	Sorts the entire list using the given entry comparison
*		function.
* \param	list			The list data structure.
* \param	entryCompFn		Given entry comparison function which
* 					must return  an integer less than,
* 					equal to, or greater than zero
* 					to indicate if the first entry
* 					is to be  considered  less than,
* 					equal to, or greater than the second.
*/
int		HGUDlpListSort(HGUDlpList *list,
			       int (*entryCompFn)(void *, void *))
{
  int		ok = 0;

  if(list && list->itemCount && list->head && list->tail)
  {
    if(HGUDlpListPvLock(list))				      /* Lock list [ */
    {
      if(list->itemCount > 1)
	HGUDlpListPvSort(list->tail, list->head, 0, list->itemCount - 1,
			 entryCompFn);
      ok = list->itemCount;
    }
    (void )HGUDlpListPvUnlock(list);			    /* Unlock list ] */
  }
  return(ok);
}

/*!
* \return	Last item to which the iterated function was applied,
* 		may be NULL on error.
* \ingroup	hguDlpList
* \brief	Iterates the given function through the list, starting
*		with the supplied given item. The iteration may proceed
*		toward either the head or tail of the list. The
*		iterated function must take the form of:
\verbatim
                  int MyItemCount(HGUDlpList *list,
                                 HGUDlpListItem *item,
                                 void *myData)
                  {
                    int         *count;

                    if(list && item)
                    {
                      count = (int *)myData;
                      ++*count;
                    }
                    return(1);
                  }
\endverbatim
*		Where the data parameter is the data supplied to this
*		(HGUDlpListIterate) function. The iteration stops when:
*		An error is encountered, either the head or tail list
*		item has be processed, or when the supplied function
*		returns zero. So the example function would count the
*		number of items, from the given item to either the
*		list head or tail.
* \param	list			The list data structure.
* \param	item			First item of iteration which
* 					may be NULL to indicate list
* 					head (dir == HGU_DLPLIST_DIR_TO_TAIL)
* 					or tail
* 					(dir == HGU_DLPLIST_DIR_TO_HEAD).
* \param	dir			Iteration direction, either
* 					towards the head or the tail.
* \param	iterFn			Function to be iterated, see
* 					example above.
* \param	iterData		Data supplied to the iterated
* 					function.
*/
HGUDlpListItem	*HGUDlpListIterate(HGUDlpList *list, HGUDlpListItem *item,
				   HGUDlpListDirection dir,
				   int (*iterFn)(HGUDlpList *,
				   		 HGUDlpListItem *, void *),
				   void *iterData)
{
  int		iterate = 1;
  HGUDlpListItem *lastItem = NULL;

  if(list && iterFn &&
     ((dir == HGU_DLPLIST_DIR_TOHEAD) || (dir == HGU_DLPLIST_DIR_TOTAIL)))
  {
    if(HGUDlpListPvLock(list))				      /* Lock list [ */
    {
      if(item == NULL)
      {
	if(dir == HGU_DLPLIST_DIR_TOHEAD)
	  item = list->tail;
	else
	  item = list->head;
      }
      while(item && iterate)
      {
	lastItem = item;
	if(dir == HGU_DLPLIST_DIR_TOHEAD)
	  item = lastItem->prev;
	else
	  item = lastItem->next;
        iterate = (*iterFn)(list, lastItem, iterData);
      }
      (void )HGUDlpListPvUnlock(list);			    /* Unlock list ] */
    }
  }
  return(lastItem);
}

/*!
* \return	Nth item from the given item, may be NULL on error or if
* 		the nth item is beyond the list head or tail.
* \ingroup	hguDlpList
* \brief	Finds the n'th item from the given item in the list.
*		The n'th item from the head or tail can be found by
*		calling the function with item == NULL, in which case
*		the direction of approach is optimised.
* \param	list			The list data structure.
* \param	item			Given item from which to find
* 					the n'th item.
* \param	dir			The direction from given item,
* 					either towards the head or the
* 					tail.
* \param	offset			The n in n'th. If zero then
* 					returns the given item, if the
* 					direction is towards the tail
* 					and offset == 1, then this is
* 					equivalent to HGUDlpListNext().
*/
HGUDlpListItem	*HGUDlpListNth(HGUDlpList *list, HGUDlpListItem *item,
			       HGUDlpListDirection dir, int offset)
{
  HGUDlpListItem *nthItem = NULL;

  if(list &&
     ((dir == HGU_DLPLIST_DIR_TOHEAD) || (dir == HGU_DLPLIST_DIR_TOTAIL)))
  {
    if(HGUDlpListPvLock(list))				      /* Lock list [ */
    {
      if(offset < 0)	    /* Allow negative item offsets, just reverse dir */
      {
	offset = 0 - offset;
        if(dir == HGU_DLPLIST_DIR_TOHEAD)
	  dir = HGU_DLPLIST_DIR_TOTAIL;
        else
	  dir = HGU_DLPLIST_DIR_TOHEAD;
      }
      if(offset > list->itemCount)        /* Check that offset is reasonable */
	item = NULL;
      else
      {
	if(item == NULL)	  /* Offset is known to be from head or tail */
	{
	  if(offset > (list->itemCount / 2)) 	       /* Optimise direction */
	  {
	    offset = list->itemCount - (offset + 1);
	    if(dir == HGU_DLPLIST_DIR_TOHEAD)
	    {
	      dir = HGU_DLPLIST_DIR_TOTAIL;
	      item = list->head;
	    }
	    else
	    {
	      dir = HGU_DLPLIST_DIR_TOHEAD;
	      item = list->tail;
	    }
	  }
	  else
	  {
	    if(dir == HGU_DLPLIST_DIR_TOHEAD)
	      item = list->tail;
	    else
	      item = list->head;
	  }
	}
	if(dir == HGU_DLPLIST_DIR_TOHEAD)		     /* Apply offset */
	{
	  while((offset > 0)  && item)
	  {
	    item = item->prev;
	    --offset;
	  }
	}
	else
	{
	  while((offset > 0)  && item)
	  {
	    item = item->next;
	    --offset;
	  }
	}
	if(offset == 0)
	  nthItem = item;
      }
      (void )HGUDlpListPvUnlock(list);			    /* Unlock list ] */
    }
  }
  return(nthItem);
}

/*!
* \return	Number of items from the given item to the head or tail
* 		item of the list.
* \ingroup	hguDlpList
* \brief	Counts the number of items from the given item to the
*		item with a NULL next/prev item, which is at the head
*		or tail of list. The offset between an item and itself
*		is defined to be zero.
* \param	list			The list data structure.
* \param	item			Given item from which to find
* 					the offset to the head or the
* 					tail.
* \param	dir			The direction from given item,
* 					either towards the head or the
* 					tail.
*/
int		HGUDlpListOffset(HGUDlpList *list, HGUDlpListItem *item,
			         HGUDlpListDirection dir)
{
  int		count = 0;

  if(list && item &&
     ((dir == HGU_DLPLIST_DIR_TOHEAD) || (dir == HGU_DLPLIST_DIR_TOTAIL)))
  {
    if(HGUDlpListPvLock(list))				      /* Lock list [ */
    {
      (void )HGUDlpListPvUnlock(list);			    /* Unlock list ] */
      if(dir == HGU_DLPLIST_DIR_TOHEAD)
      {
        while((item = item->prev) != NULL)
	  ++count;
      }
      else
      {
        while((item = item->next) != NULL)
	  ++count;
      }
    }
  }
  return(count);
}

/*!
* \return	Non zero if the item is at the head of the list.
* \ingroup	hguDlpList
* \brief	Looks to see if the given item is at the head of the
*		given list.
* \param	list			The list data structure.
* \param	item			Is this item at the head of the
* 					list?
*/
int		HGUDlpListItemIsHead(HGUDlpList *list, HGUDlpListItem *item)
{
  int		atHead = 0;

  if(list && item)
  {
    if(HGUDlpListPvLock(list))				      /* Lock list [ */
    {
      if(item->prev == 0)
      {
	assert((unsigned long )(list->head) == (unsigned long )item);
	atHead = 1;
      }
    }
    (void )HGUDlpListPvUnlock(list);			    /* Unlock list ] */
  }
  return(atHead);
}

/*!
* \return	Non zero if the item is at the tail of the list.
* \ingroup	hguDlpList
* \brief	Looks to see if the given item is at the tail of the
*		given list.
* \param	list			The list data structure.
* \param	item			Is this item at the tail of the
* 					list?
*/
int		HGUDlpListItemIsTail(HGUDlpList *list, HGUDlpListItem *item)
{
  int		atTail = 0;

  if(list && item)
  {
    if(HGUDlpListPvLock(list))				      /* Lock list [ */
    {
      if(item->next == 0)
      {
	assert((unsigned long )(list->tail) == (unsigned long )item);
	atTail = 1;
      }
    }
    (void )HGUDlpListPvUnlock(list);			    /* Unlock list ] */
  }
  return(atTail);
}

/*!
* \return	The list item's entry, which may be NULL for an empty
* 		list or on error.
* \ingroup	hguDlpList
* \brief	Returns the list items entry
* \param	list			The list data structure.
* \param	item			Item with the entry.
*/
void		*HGUDlpListEntryGet(HGUDlpList *list,
				    HGUDlpListItem *item)
{
  void		*entry = NULL;

  if(list && item)
    entry = item->entry;
  return(entry);
}

/*!
* \return	The previous entry.
* \ingroup	hguDlpList
* \brief	Sets the given items entry and returns the previous
*		entry. Entries are NEVER freed by this function.
* \param	list			The list data structure.
* \param	item			Item for which to change entry.
* \param	newEntry		New entry to be interted into the
* 					given item.
*/
void		*HGUDlpListEntrySet(HGUDlpList *list,
				    HGUDlpListItem *item, void *newEntry)
{
  void		*oldEntry;

  if(list && item)
  {
    if(HGUDlpListPvLock(list))				      /* Lock list [ */
    {
      oldEntry = item->entry;
      item->entry = newEntry;
      (void )HGUDlpListPvUnlock(list);			    /* Unlock list ] */
    }
  }
  return(oldEntry);
}

/*!
* \return	List item at tail of list, may may be NULL either
* 		on error or if the list is empty.
* \ingroup	hguDlpList
* \brief	Returns the tail list item.
* \param	list			The list data structure.
*/
HGUDlpListItem	*HGUDlpListTail(HGUDlpList *list)
{
  HGUDlpListItem *tail = NULL;

  if(list)
    tail = list->tail;
  return(tail);
}

/*!
* \return	List item at head of list, may may be NULL either
* 		on error or if the list is empty.
* \ingroup	hguDlpList
* \brief	Returns the head list item.
* \param	list			The list data structure.
*/
HGUDlpListItem	*HGUDlpListHead(HGUDlpList *list)
{
  HGUDlpListItem *head = NULL;

  if(list)
    head = list->head;
  return(head);
}

/*!
* \return	Next list item after the given item, may be NULL on error,
* 		if the list is empty or if the given item was at the tail
* 		of the list.
* \ingroup	hguDlpList
* \brief	Returns the next list item.
* \param	list			The list data structure.
* \param	item			Item with the entry.
*/
HGUDlpListItem	*HGUDlpListNext(HGUDlpList *list, HGUDlpListItem *item)
{
  HGUDlpListItem *next = NULL;

  if(list && item)
    next = item->next;
  return(next);
}

/*!
* \return	Prev list item before the given item, may be NULL on
* 		error, if the list is empty or if the given item was
* 		at the head of the list.
* \ingroup	hguDlpList
* \brief	Returns the prev list item.
* \param	list			The list data structure.
* \param	item			Item with the entry.
*/
HGUDlpListItem	*HGUDlpListPrev(HGUDlpList *list, HGUDlpListItem *item)
{
  HGUDlpListItem *prev = NULL;

  if(list && item)
    prev = item->prev;
  return(prev);
}

/*!
* \return	Number of items in list. This is >= zero for a valid
* 		list but may be negative if the list is invalid.
* \ingroup	hguDlpList
* \brief	Returns the number of items in the list.
* \param	list			The list data structure.
*/
int		HGUDlpListCount(HGUDlpList *list)
{
  int		count = -1;
  if(list)
    count = list->itemCount;
  return(count);
}


/*!
* \return	Non zero if locked or lock is not required.
* \ingroup	hguDlpList
* \brief	Private function: Locks list data structures using the
*		installed lock function.
* \param	list			The list data structure.
*/
static int	HGUDlpListPvLock(HGUDlpList *list)
{
  int		locked = 0;

  if(list)
  {
    if((list->lockFn == NULL) ||
       ((*list->lockFn)(list->lockData,
    			HGU_DLPLIST_STATE_LOCK) == HGU_DLPLIST_STATE_LOCK))
      locked = 1;
  }
  return(locked);
}

/*!
* \return	Non zero if unlocked or lock is not required.
* \ingroup	hguDlpList
* \brief	Unlocks list data structures using the installed lock
*		function.
* \param	list			The list data structure.
*/
static int	HGUDlpListPvUnlock(HGUDlpList *list)
{
  int		unlocked = 0;

  if(list)
  {
    if((list->lockFn == NULL) ||
       ((*list->lockFn)(list->lockData, HGU_DLPLIST_STATE_UNLOCK) ==
        HGU_DLPLIST_STATE_UNLOCK))
      unlocked = 1;
  }
  return(unlocked);
}

/*!
* \return	New item with entry and free function set.
* \ingroup	hguDlpList
* \brief	Creates a new item with the entry and free function
*		set.
* \param	entry			New list entry.
* \param	freeFn			Function that will be called
* 					(if not NULL) to free the entry
* 					should the item be deleted.
*/
static HGUDlpListItem *HGUDlpListPvItemCreate(void *entry,
					      void (*freeFn)(void *))
{
  HGUDlpListItem *newItem;

  newItem = (HGUDlpListItem *)malloc(sizeof(HGUDlpListItem));
  if(newItem)
  {
    newItem->freeFn = freeFn;
    newItem->entry = entry;
    newItem->next = NULL;
    newItem->prev = NULL;
  }
  return(newItem);
}

/*!
* \return	New list head which is NULL if all the items were deleted
* 		ok.
* \ingroup	hguDlpList
* \brief	Deletes all list items from the head on down to and
*		including the tail. Where delete implies both the
*		removal of items from the list and freeing the entries
*		using the item's free functions (unless either the free
*		function or the entry is NULL).
*		Note that this function assumes that all the list data
*		structures can be accessed without acquiring a lock.
* \param	list			The list data structure.
*/
static void	HGUDlpListPvDeleteAll(HGUDlpList *list)
{
  HGUDlpListItem *tItem0,
  		*tItem1;

  if(list)
  {
    tItem0 = list->head;
    while(tItem0)
    {
      if(tItem0->freeFn && tItem0->entry)
	(*(tItem0->freeFn))(tItem0->entry);
      tItem1 = tItem0;
      tItem0 = tItem1->next;
      free(tItem1);
    }
    list->head = NULL;
    list->tail = NULL;
    list->itemCount = 0;
  }
}

/*!
* \return	Next list item after the item deleted, may be NULL either
* 		on error or if the list item was the last in the list.
* \ingroup	hguDlpList
* \brief	Deletes or removes the given list item from the list
*		with the given list. Where delete implies both the
*		removal of the item, the freeing of the item and (when
*		neither the free function or the entry are NULL) the
*		freeing of the list entry too. Remove never frees the
*		list entry.
* \param	list			The list data structure.
* \param	item			Item to be deleted or removed.
* \param	deleteFlag		Delete rather than remove item
* 					if non zero.
*/
static HGUDlpListItem *HGUDlpListPvDelRmItem(HGUDlpList *list,
				             HGUDlpListItem *item,
					     int deleteFlag)
{
  HGUDlpListItem *tItem0,
		*tItem1,
		*freeItem = NULL,
  		*nextItem = NULL;

  if(list && item)
  {
    if(HGUDlpListPvLock(list))				      /* Lock list [ */
    {
      if(list->head || list->tail)
      {
	assert(list->itemCount >= 1);
	assert((list->head != NULL) && (list->tail != NULL));
	if(item->prev == NULL)
	{
	  if(item->next == NULL)			    /* Only one item */
	  {
	    assert(list->itemCount == 1);
	    freeItem = item;
	    list->head = NULL;
	    list->tail = NULL;
	  }
	  else			    /* Item to be deleted is at head of list */
	  {
	    assert(list->itemCount >= 2);
	    assert((unsigned long )item == (unsigned long )(list->head));
	    freeItem = item;
	    tItem0 = item->next;
	    tItem0->prev = NULL;
	    list->head = tItem0;
	    nextItem = tItem0;
	  }
	}
	else if(item->next == NULL) /* Item to be deleted is at tail of list */
	{
	  assert(list->itemCount >= 2);
	  assert((unsigned long )item == (unsigned long )(list->tail));
	  freeItem = item;
	  tItem0 = item->prev;
	  tItem0->next = NULL;
	  list->tail = tItem0;
	  nextItem = NULL;
	}
	else 			     /* Item not at head or tail of the list */
	{
	  assert(list->itemCount >= 3);
	  freeItem = item;
	  tItem0 = item->prev;
	  tItem1 = item->next;
	  tItem0->next = tItem1;
	  tItem1->prev = tItem0;
	  nextItem = tItem1;
	}
	if(freeItem)
	{
	  if(deleteFlag && freeItem->freeFn && freeItem->entry)
	    (*(freeItem->freeFn))(freeItem->entry);
	  free(freeItem);
	  --(list->itemCount);
	}
      }
      (void )HGUDlpListPvUnlock(list);			    /* Unlock list ] */
    }
  }
  return(nextItem);
}

/*!
* \return	void
* \ingroup	hguDlpList
* \brief	Sorts the given section of the list using a modified
*		quick sort algorithm. This function is RECURSIVE. The
*		values of the indicies must be such that lowIdx is less
*		than highIdx and their difference is equal to one more
*		than the number of items between them (excluding the
*		items themselves that is).
* \param	low			Item nearest to tail in section
* 					to be sorted.
* \param	high			Item nearest to head in section
* 					to be sorted.
* \param	lowIdx			Index for item nearest to tail.
* \param	highIdx			Index for item nearest to head.
* \param	entryCompFn		Given entry comparison function
* 					which must return an integer less
* 					than, equal to, or greater than zero
* 					to indicate if the first entry
* 					is to be  considered less than,
* 					equal to, or greater than the second.
*/
static void	HGUDlpListPvSort(HGUDlpListItem *low, HGUDlpListItem *high,
				 int lowIdx, int highIdx,
				 int (*entryCompFn)(void *, void *))
{

  HGUDlpListItem *pivot,
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
      HGUDlpListPvSort(low, pivot, lowIdx, pivotIdx, entryCompFn);
    if(baseIdx < highIdx)
      HGUDlpListPvSort(base, high, baseIdx, highIdx, entryCompFn);
  }
  else
  {
    if(baseIdx < highIdx)
      HGUDlpListPvSort(base, high, baseIdx, highIdx, entryCompFn);
    if(lowIdx < pivotIdx)
      HGUDlpListPvSort(low, pivot, lowIdx, pivotIdx, entryCompFn);
  }
}

#undef HGUDLPLIST_C
