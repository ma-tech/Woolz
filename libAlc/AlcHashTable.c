#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _AlcHashTable_c[] = "MRC HGU $Id$";
#endif
#endif

/*!
* \file         libAlc/AlcHashTable.c
* \author       Bill Hill
* \date         November 1999
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
* \brief        A general purpose hash table.
* \ingroup	AlcHashTable
* \todo		-
* \bug          None known.
*/
#include <stdio.h>
#include <stdlib.h>
#include <Alc.h>

static AlcHashItem 		*AlcHashItemListLast(
				  AlcHashItem *item);

/*!
* \return	Hash table structure, or NULL on error.
* \ingroup	AlcHashTable
* \brief	Creates a hash table data structure which is required
*		by all the other AlcHashTable functions.
* \param	tableSz 		Hash table size.
* \param	keyCmp			Key comparison function.
* \param	hashFn			Hash function for table.
* \param	dstErr			Destination pointer for error
*					code, may be NULL.
*/
AlcHashTable	*AlcHashTableNew(size_t tableSz,
				 int (*keyCmp)(void *, void *),
				 unsigned (*hashFn)(void *),
			         AlcErrno *dstErr)
{
  AlcHashTable	*hTbl = NULL;
  AlcHashItem	**items = NULL;
  AlcErrno	errNum = ALC_ER_NONE;

  if(((hTbl = (AlcHashTable *)AlcCalloc(1, sizeof(AlcHashTable))) == NULL) ||
     ((items = (AlcHashItem **)AlcCalloc(tableSz,
     					 sizeof(AlcHashItem *))) == NULL))
  {
    if(hTbl)
    {
      AlcFree(hTbl);
    }
    if(items)
    {
      AlcFree(items);
    }
    errNum = ALC_ER_ALLOC;
  }
  else
  {
    hTbl->keyCmp = keyCmp;
    hTbl->hashFn = hashFn;
    hTbl->tableSz = tableSz;
    hTbl->table = items;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(hTbl);
}

/*!
* \return	Hash item structure, or	NULL on error.
* \ingroup	AlcHashTable
* \brief	Creates a hash item data structure for building into
*		an AlcHashTable.
* \param	entry			New table entry.
* \param	freeFn			Function that will be called
*					(if not NULL) to free the entry.
* \param	key			Hash key.
* \param	dstErr			Destination pointer for error
*					code, may be NULL.
*/
AlcHashItem	*AlcHashItemNew(void *entry, void (*freeFn)(void *),
			       void *key, AlcErrno *dstErr)
{
  AlcHashItem 	*newItem;
  AlcErrno	errNum = ALC_ER_NONE;

  if((newItem = (AlcHashItem *)AlcCalloc(1, sizeof(AlcHashItem))) == NULL)
  {
    errNum = ALC_ER_ALLOC;
  }
  else
  {
    newItem->freeFn = freeFn;
    newItem->entry = entry;
    newItem->key = key;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(newItem);
}

/*!
* \return	Error code.
* \ingroup	AlcHashTable
* \brief	Free's the given hash table data structure and any hash
*		table items.
* \param	hTbl 			The hash table data structure.
*/
AlcErrno	AlcHashTableFree(AlcHashTable *hTbl)
{
  size_t	hCnt;
  AlcHashItem	*item,
		*fItem;
  AlcHashItem 	**head;
  AlcErrno	errNum = ALC_ER_NONE;

  if(hTbl == NULL)
  {
    errNum =  ALC_ER_NULLPTR;
  }
  else
  {
    head = hTbl->table;
    hCnt = hTbl->tableSz;
    while(hCnt-- > 0)
    {
      if((item = *head++) != NULL)
      {
	while(item)
	{
	  fItem = item;
	  item = item->next;
	  if(fItem->freeFn && fItem->entry)
	  {
	    (*(fItem->freeFn))(fItem->entry);
	  }
	  AlcFree(fItem);
	}
      }
    }
    AlcFree(hTbl);
  }
  return(errNum);
}

/*!
* \return	Error code.
* \ingroup	AlcHashTable
* \brief	Inserts the given entry into the hash table.
* \param	hTbl 			The hash table data structure.
* \param	key			Hash key.
* \param	entry			New list entry.
* \param	freeFn			Function that will be called
*					(if not NULL) to free the entry
*					should the item be deleted.
*/
AlcErrno	AlcHashTableEntryInsert(AlcHashTable *hTbl, void *key,
				      void *entry, void (*freeFn)(void *))
{
  AlcHashItem	*newItem;
  AlcErrno	errNum = ALC_ER_NONE;

  if((hTbl == NULL) || (hTbl->hashFn == NULL))
  {
    errNum = ALC_ER_NULLPTR;
  }
  else
  {
    if((newItem = AlcHashItemNew(entry, freeFn, key, &errNum)) != NULL)
    {
      if((errNum = AlcHashItemInsert(hTbl, newItem)) != ALC_ER_NONE)
      {
        (void )AlcHashItemFree(newItem);
      }
    }
  }
  return(errNum);
}

/*!
* \return	Error code.
* \ingroup	AlcHashTable
* \brief	Removes the item from the hash table, but does not free
*		the item unless the freeItem flag is set.
* \param	hTbl			The hash table data structure.
* \param	rItem			Item to be removed.
* \param	freeItem		Free item if non-zero.
*/
AlcErrno	AlcHashItemUnlink(AlcHashTable *hTbl, AlcHashItem *rItem,
				  int freeItem)
{
  size_t	hIdx;
  AlcHashItem	**head;
  AlcHashItem	*pItem,
  		*nItem;
  AlcErrno	errNum = ALC_ER_NONE;

  if((hTbl == NULL) || (rItem == NULL))
  {
    errNum = ALC_ER_NULLPTR;
  }
  else
  {
    nItem = rItem->next;
    if((pItem = rItem->prev) == NULL)
    {
      hIdx = (*(hTbl->hashFn))(rItem->key) % hTbl->tableSz;
      head = hTbl->table + hIdx;
      *head = nItem;
    }
    else
    {
      pItem->next = nItem;
    }
    if(nItem)
    {
      nItem->prev = pItem;
    }
    if(freeItem)
    {
      (void )AlcHashItemFree(rItem);
    }
  }
  return(errNum);
}

/*!
* \return	Error code.
* \ingroup	AlcHashTable
* \brief	Inserts a new item into the hash table.
*		First find the table list head by generating an
*		index from the key using the hash function, then
*		insert the entry into the sorted list.
* \param	hTbl			The hash table data structure.
* \param	newItem			New item to insert.
*/
AlcErrno	AlcHashItemInsert(AlcHashTable *hTbl, AlcHashItem *newItem)
{
  int		cmp;
  size_t	hIdx;
  AlcHashItem 	**head;
  AlcHashItem	*item,
  		*pItem,
  		*nItem;
  AlcErrno	errNum = ALC_ER_NONE;

  if((hTbl == NULL) || (newItem == NULL))
  {
    errNum = ALC_ER_NULLPTR;
  }
  else
  {
    newItem->next = NULL;
    newItem->prev = NULL;
    hIdx = (*(hTbl->hashFn))(newItem->key) % hTbl->tableSz;
    head = hTbl->table + hIdx;
    if(*head == NULL)			     		 /* 0 items at index */
    {
      *head = newItem;
    }
    else if((item = *head)->next == NULL)	          /* 1 item at index */
    {
      if((cmp = (*(hTbl->keyCmp))(item->key, newItem->key)) > 0)
      {
	/* Insert new item before the list head. */
	nItem = *head;
	*head = newItem;
	newItem->next = nItem;
	nItem->prev = newItem;
      }
      else if(cmp == 0)
      {
	/* Replace list head with the new item. */
	item = *head;
	*head = newItem;
	(void )AlcHashItemFree(item);
      }
      else if(cmp < 0)
      {
	/* Insert new item after the list head. */
	pItem = *head;
	pItem->next = newItem;
	newItem->prev = pItem;
      }
    }
    else						/* > 1 item at index */
    {
      /* Search for position in ordered list. */
      pItem = NULL;
      nItem = NULL;
      while(((cmp = (*(hTbl->keyCmp))(item->key, newItem->key)) < 0) &&
	    ((nItem = item->next)  != NULL))
      {
	pItem = item;
	item = nItem;
      }
      /* Insert new item into list. */
      if(cmp > 0)
      {
	/* Insert new item before matched item */
	if(pItem == NULL)
	{
	  /* Insert new item at the head of the list. */
	  item = *head;
	  *head = newItem;
	}
	else
	{
	  pItem->next = newItem;
	  newItem->prev = pItem;
	}
	newItem->next = item;
	item->prev = newItem;
      }
      else if(cmp == 0)
      {
	/* Replace matched item with the new item */
	if(pItem == NULL)
	{
	  *head = newItem;
	}
	else
	{
	  pItem->next = newItem;
	  newItem->prev = pItem;
	}
	if((nItem = item->next) != NULL)
	{
	  nItem->prev = newItem;
	}
	newItem->next = nItem;
        (void )AlcHashItemFree(item);
      }
      else if(cmp < 0)
      {
	/* Insert new item at the end of the list. */
	item->next = newItem;
	newItem->prev = item;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Error code.
* \ingroup	AlcHashTable
* \brief	Free's the list item which has already been removed
*		from the list.
* \param	item			Item to be deleted.
*/
AlcErrno	AlcHashItemFree(AlcHashItem *item)
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
* \return	Number of items in list. This is always >= 0.
* \ingroup	AlcHashTable
* \brief	Returns the number of items in the list.
* \param	hTbl			The hash table data structure.
* \param	dstErr			Destination pointer for error
*					code, may be NULL.
*/
size_t		AlcHashTableCount(AlcHashTable *hTbl, AlcErrno *dstErr)
{
  size_t	cnt = 0,
  		hCnt;
  AlcHashItem 	*item,
  		**head;
  AlcErrno	errNum = ALC_ER_NONE;

  if(hTbl == NULL)
  {
    errNum = ALC_ER_NULLPTR;
  }
  else
  {
    hCnt = hTbl->tableSz;
    head = hTbl->table;
    while(hCnt-- > 0)
    {
      item = *head;
      while(item)
      {
        ++cnt;
	item = item->next;
      }
      ++head;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cnt);
}

/*!
* \return	Last item.
* \ingroup	AlcHashTable
* \brief	Iterates the given function through all entries of
*		the hash table, starting with either the first or
*		last item. The iteration proceeds towards either the
*		last or first item in the table.
*		The iterated function must take the form
* \verbatim
*		  int MyItemCount(AlcHashTable *hTbl,
*				  ALC_DIRECTION_FWD,
*				  void *myData)
*		  {
*		    int		*count;
*
*		    if(hTbl && item)
*		    {
*		      if(keepGoing)
*		      {
*		        count = (int *)myData;
*		        ++*count;
*		      }
*		    }
*		    return(1);
*		  }
* \endverbatim
*		Where the data parameter is the data supplied to this
*		(AlcHashTableIterate) function.
*		The iteration continues until either the iterated
*		function returns zero or the last/first item of
*		the hash table has been processed.
* \param	hTbl			The hash table data structure.
* \param	dir			Iteration direction.
* \param	iterFn			Function to be iterated, see
*					example above.
* \param	iterData		Data supplied to the iterated
*					function.
* \param	dstErr			Destination pointer for error
*					code, may be NULL.
*/
AlcHashItem	*AlcHashTableIterate(AlcHashTable *hTbl, AlcDirection dir,
				     int (*iterFn)(AlcHashTable *,
				                   AlcHashItem *, void *),
				     void *iterData, AlcErrno *dstErr)
{
  int		iterate = 1;
  AlcHashItem	**head;
  AlcHashItem 	*item,
  		*lastItem = NULL;
  AlcErrno	errNum = ALC_ER_NONE;

  if((hTbl == NULL) || (iterFn == NULL))
  {
    errNum = ALC_ER_NULLPTR;
  }
  else if(hTbl->table)
  {
    switch(dir)
    {
      case ALC_DIRECTION_FWD:
	head = hTbl->table;
	item = *head;
	break;
      case ALC_DIRECTION_REV:
	head = hTbl->table + hTbl->tableSz - 1;
	item = AlcHashItemListLast(*head);
	break;
      default:
	iterate = 0;
	errNum = ALC_ER_PARAM;
	break;
    }
    while(iterate)
    {
      if(item)
      {
        iterate = (*iterFn)(hTbl, item, iterData);
        lastItem = item;
      }
      if(dir == ALC_DIRECTION_FWD)
      { 
        if(item)
	{
	  item = item->next;
	}
	if(item == NULL)
	{
	  if(++head >= hTbl->table + hTbl->tableSz)
	  {
	    iterate = 0;
	  }
	  else
	  {
	    item = *head;
	  }
	}
      }
      else /* dir == ALC_DIRECTION_REV */
      {
        if(item)
	{
	  item = item->prev;
	}
	if(item == NULL)
	{
	  if(--head < hTbl->table)
	  {
	    iterate = 0;
	  }
	  else
	  {
	    item = AlcHashItemListLast(*head);
	  }
	}
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(lastItem);
}

/*!
* \return	Error code.
* \ingroup	AlcHashTable
* \brief	Unlinks all items which are matched by the given
*		match function. If the match pointer function is
*		NULL then all the items are unlinked.
*		The test function should return non-zero to unlink
*		the item passed to it.
* \param	hTbl			The hash table data structure.
* \param	testFn			Function to test items.
* \param	fnData			Data supplied to the test function.
* \param	freeItems		Free unlinked items if non-zero.
*/
AlcErrno	AlcHashTableUnlinkAll(AlcHashTable *hTbl,
				      int (*testFn)(AlcHashTable *,
				      		    AlcHashItem *, void *),
				      void *fnData, int freeItems)
{
  size_t	hIdx;
  AlcHashItem	**head;
  AlcHashItem 	*item,
  		*nItem;
  AlcErrno	errNum = ALC_ER_NONE;

  if(hTbl == NULL)
  {
    errNum = ALC_ER_NULLPTR;
  }
  else if(hTbl->table)
  {
    hIdx = 0;
    head = hTbl->table;
    item = *head;
    while(hIdx < hTbl->tableSz)
    {
      item = *head;
      while(item)
      {
        nItem = item->next;
        if((*testFn)(hTbl, item, fnData))
	{
	  (void )AlcHashItemUnlink(hTbl, item, freeItems);
	}
	item = nItem;
      }
      ++head;
      ++hIdx;
    }
  }
  return(errNum);
}

/*!
* \return	Matched item, NULL if no match.
* \ingroup	AlcHashTable
* \brief	Gets the item from the hash table with the matching
*		key, it's not an error if a matching item isn't found.
* \param	hTbl			The hash table data structure.
* \param	key			Given key to match.
* \param	dstErr			Destination pointer for error
*					code, may be NULL.
*/
AlcHashItem	*AlcHashItemGet(AlcHashTable *hTbl, void *key,
				AlcErrno *dstErr)
{
  int		cmp;
  size_t	hIdx;
  AlcHashItem 	**head;
  AlcHashItem	*nItem,
  		*item = NULL;
  AlcErrno	errNum = ALC_ER_NONE;

  if(hTbl == NULL)
  {
    errNum = ALC_ER_NULLPTR;
  }
  else
  {
    hIdx = (*(hTbl->hashFn))(key) % hTbl->tableSz;
    head = hTbl->table + hIdx;
    if(*head)
    {
      item = *head;
      while(((cmp = (*(hTbl->keyCmp))(item->key, key)) < 0) &&
	    ((nItem = item->next)  != NULL))
      {
        item = nItem;
      }
      if(cmp != 0)
      {
        item = NULL;
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(item);
}

/*!
* \return	>0 if item0 before item1 or
*		== 0 if item0 is item1 or
*		< 0 if item0 after item1.
* \ingroup	AlcHashTable
* \brief	Finds the order in which the given items would occur
*		in the hash table.
* \param	hTbl			The hash table.
* \param	item0			Given item0.
* \param	item1			Given item1.
*/
int		AlcHashItemOrder(AlcHashTable *hTbl,
				 AlcHashItem *item0, AlcHashItem *item1)
{
  int		order;
  size_t 	hIdx0,
  		hIdx1,
		hVal0,
  		hVal1;

  hVal0 = (*(hTbl->hashFn))(item0->key);
  hIdx0 = hVal0 % hTbl->tableSz;
  hVal1 = (*(hTbl->hashFn))(item1->key);
  hIdx1 = hVal1 % hTbl->tableSz;
  if(hIdx1 == hIdx0)
  {
    if(hVal1 == hVal0)
    {
      order = 0;
    }
    else if(hVal1 > hVal0)
    {
      order = 1;
    }
    else /* hVal1 < hVal0 */
    {
      order = -1;
    }
  }
  else if(hIdx1 > hIdx0)
  {
    order = 1;
  }
  else /* hIdx1 < hIdx0 */
  {
    order = -1;
  }
  return(order);
}

/*!
* \return	Last item (from given) in list, NULL if given NULL.
* \ingroup	AlcHashTable
* \brief	Finds the last item in the list from the given item.
* \param	item			Given item.
*/
static AlcHashItem *AlcHashItemListLast(AlcHashItem *item)
{
  AlcHashItem	*nItem;

  if(item)
  {
    while((nItem = item->next) != NULL)
    {
      item = nItem;
    }
  }
  return(item);
}
