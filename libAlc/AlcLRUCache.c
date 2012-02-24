#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlcLRUCache_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlc/AlcLRUCache.c
* \author       Bill Hill
* \date         November 2011
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
* \brief	A cache allowing rank and key access to it's entries.
* 		The cache is maintained using a maximum number of items and
* 		maximum total entry size.
* \ingroup	AlcLRUCache
*/

#include <string.h>
#include <stdio.h>
#include <Alc.h>

static AlcLRUCItem 		*AlcLRUCItemNew(
				  AlcLRUCache *cache);
static unsigned int 		AlcLRUCHash(
				  unsigned int k);
static void      		AlcLRUCItemFree(
				  AlcLRUCache *cache,
				  AlcLRUCItem *item);
static void     		AlcLRUCCheckLimits(
				  AlcLRUCache *cache,
				  size_t newEntrySz,
				  unsigned int numNewEntry);
static void     		AlcLRUCMoreItems(
				  AlcLRUCache *cache);
static void     		AlcLRUCItemRankUnlink(
				  AlcLRUCache *cache,
				  AlcLRUCItem *item);
static void     		AlcLRUCItemRankInsert(
				  AlcLRUCache *cache,
				  AlcLRUCItem *item);
static void     		AlcLRUCItemHashInsert(
				  AlcLRUCache *cache,
				  AlcLRUCItem *item);
static void     		AlcLRUCItemHashUnlink(
				  AlcLRUCache *cache,
				  AlcLRUCItem *item);
static void     		AlcLRUCRemoveItem(
				  AlcLRUCache *cache,
				  AlcLRUCItem *item);

/*!
* \return	New least recent use cache or NULL on error.
* \ingroup	AlcLRUCache
* \brief	Allocates a new least recent use removal cache.
* \param	maxItem			Maximum number of items in cache.
* \param	maxSz			Maximum total cache entry size.
* \param	keyFn			Supplied function for computing
* 					a numeric key from cache entries.
* \param	cmpFn			Supplied function for matching cache
* 					entries.
* \param	unlinkFn		Optional supplied function that is
* 					called prior to unlinking a cache
* 					item and removing it from the cache.
* \param	dstErr			Destination error pointer, may be NULL.
*/
AlcLRUCache	*AlcLRUCacheNew(unsigned int maxItem, size_t maxSz, 
                                AlcLRUCKeyFn keyFn, AlcLRUCCmpFn cmpFn,
				AlcLRUCUnlinkFn unlinkFn, AlcErrno *dstErr)
{
  AlcLRUCache	*cache = NULL;
  AlcErrno	errNum = ALC_ER_NONE;

  if((cache = (AlcLRUCache *)AlcCalloc(1, sizeof(AlcLRUCache))) == NULL)
  {
    errNum = ALC_ER_ALLOC;
  }
  else
  {
    cache->maxItem = maxItem;
    cache->maxSz = maxSz;
    cache->keyFn = keyFn;
    cache->cmpFn = cmpFn;
    cache->unlinkFn = unlinkFn;
    cache->hashTblSz = (maxItem > 256)? maxItem / 4: 64;
    cache->itemBlkSz = cache->hashTblSz;
    if((cache->hashTbl = (AlcLRUCItem **)
                         AlcCalloc(cache->hashTblSz,
    				   sizeof(AlcLRUCache *))) == NULL)
    {
      AlcFree(cache);
      cache = NULL;
      errNum = ALC_ER_ALLOC;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cache);
}

/*!
* \ingroup	AlcLRUCache
* \brief	Frees a least recent use removal cache.
* \param	cache			The cache.
* \param	unlink			Flag, which if non-zero, will remove
* 					all items from the cache before
* 					freeing the cache. If the unlink
* 					function is to be called for the
* 					item entries as they are removed from
* 					the cache then this flag must be set.
*/
void       	AlcLRUCacheFree(AlcLRUCache *cache, int unlink)
{
  if(cache)
  {
    if(unlink)
    {
      AlcLRUCEntryRemoveAll(cache);
    }
    (void )AlcBlockStackFree(cache->freeStack);
    AlcFree(cache->hashTbl);
    AlcFree(cache);
  }
}

/*!
* \return	Cache entry or NULL if not found.
* \ingroup	AlcLRUCache
* \brief	Given a cache entry (with sufficient data in the entry
* 		for the cache key generation function) this function attempts
* 		to find a matching entry in the cache. The partical entry could
* 		for example contain only an identification string if that is
* 		all that is required for key generation and entry comparison.
* 		If the entry exists in the cache then it become the most
* 		recently used entry.
* 		See AlcLRUCEntryGetWithKey().
* \param	cache			The cache.
* \param	entry			Given partial cache entry to be
* 					matched.
*/
void		*AlcLRUCEntryGet(AlcLRUCache *cache, void *entry)
{
  unsigned int	key;
  void		*found;

  key = (*(cache->keyFn))(cache, entry);
  found = AlcLRUCEntryGetWithKey(cache, key, entry);
  return(found);
}

/*!
* \return	Cache entry or NULL if not found.
* \ingroup	AlcLRUCache
* \brief	Given a cache entry with sufficient data in the entry
* 		for the cache key generation function, this function attempts
* 		to find the matching cache entry. The partical entry could
* 		for example contain only an identification string if that is
* 		all that is required for entry comparison.
* 		If The given key is used to search the cache for entries
* 		with the same key. These are then checked using the cache
* 		entry comparison function.
* 		If the entry exists in the cache then it become the most
* 		recently used entry.
* \param	cache			The cache.
* \param	key			Key generated by the cache key
* 					generation function for the entry.
* \param	entry			Given partial cache entry to be
* 					matched.
*/
void		*AlcLRUCEntryGetWithKey(AlcLRUCache *cache, unsigned int key,
				        void *entry)
{
  void		*found = NULL;
  AlcLRUCItem	*item;

  item = AlcLRUCItemFind(cache, key, entry);
  if(item)
  {
    AlcLRUCItemRankUnlink(cache, item);
    AlcLRUCItemRankInsert(cache, item);
    found = item->entry;
  }
  return(found);
}

/*!
* \return	The cache item or NULL if not in cache on return.
* \ingroup	AlcLRUCache
* \brief	Attempts to add the given entry to the cache. If the
* 		entry already exists in the cache then it become the
* 		most recently used entry, but is not modified.
* \param	cache			The cache.
* \param	entrySz			Size of cache entry for use in
* 					limiting total cache entry size.
* \param	entry			Given entry to be added to the cache.
* \param	dstNewFlg		Destination pointer set to zero or
* 					non-zero. Set to non-zero only if
* 					a new item is created. May be NULL.
*/
AlcLRUCItem	*AlcLRUCEntryAdd(AlcLRUCache *cache, size_t entrySz,
				void *entry, int *dstNewFlg)
{
  int		newFlg = 0;
  unsigned int 	key;
  AlcLRUCItem	*item = NULL;

  if(cache && ((cache->maxSz == 0) || (entrySz < cache->maxSz)))
  {
    key = (*(cache->keyFn))(cache, entry);
    item = AlcLRUCItemFind(cache, key, entry);
    if(item == NULL)
    {
      if((item = AlcLRUCItemNew(cache)) != NULL)
      {
	newFlg = 1;
	item->key = key;
	item->sz = entrySz;
	item->entry = entry;
	AlcLRUCCheckLimits(cache, entrySz, 1);
	AlcLRUCItemRankInsert(cache, item);
	AlcLRUCItemHashInsert(cache, item);
	++(cache->numItem);
	cache->curSz += entrySz;
      }
    }
  }
  if(dstNewFlg)
  {
    *dstNewFlg = newFlg;
  }
  return(item);
}

/*!
* \return	The cache item or NULL if not in cache on return.
* \ingroup	AlcLRUCache
* \brief	Attempts to add the given entry to the cache. If the
* 		entry already exists in the cache then it become the
* 		most recently used entry.
* \param	cache			The cache.
* \param	entrySz			Size of cache entry for use in
* 					limiting total cache entry size.
* \param	entry			Given entry to be added to the cache.
* \param	dstNewFlg		Destination pointer set to zero or
* 					non-zero. Set to non-zero only if
* 					a new item is created. May be NULL.
*/
AlcLRUCItem	*AlcLRUCEntryAddWithKey(AlcLRUCache *cache, size_t entrySz,
				        void *entry, unsigned int key,
					int *dstNewFlg)
{
  int		newFlg = 0;
  AlcLRUCItem	*item = NULL;

  if(cache && ((cache->maxSz == 0) || (entrySz < cache->maxSz)))
  {
    item = AlcLRUCItemFind(cache, key, entry);
    if(item == NULL)
    {
      if((item = AlcLRUCItemNew(cache)) != NULL)
      {
	newFlg = 1;
	item->key = key;
	item->sz = entrySz;
	item->entry = entry;
	AlcLRUCCheckLimits(cache, entrySz, 1);
	AlcLRUCItemRankInsert(cache, item);
	AlcLRUCItemHashInsert(cache, item);
	++(cache->numItem);
	cache->curSz += entrySz;
      }
    }
  }
  if(dstNewFlg)
  {
    *dstNewFlg = newFlg;
  }
  return(item);
}

/*!
* \ingroup	AlcLRUCache
* \brief	Removes the matching cache entry from the queue.
* 		See AlcLRUCEntryGet() for the entry use.
* \param	cache			The cache.
* \param	entry			Given partial cache entry to be
* 					matched.
*/
void		AlcLRUCEntryRemove(AlcLRUCache *cache, void *entry)
{
  unsigned int	key;

  key = (*(cache->keyFn))(cache, entry);
  AlcLRUCEntryRemoveWithKey(cache, key, entry);
}

/*!
* \ingroup	AlcLRUCache
* \brief	Removes the matching cache entry from the queue.
* 		See AlcLRUCEntryGetWithKey() for the entry use.
* \param	cache			The cache.
* \param	ley			Key generated from given entry.
* \param	entry			Given partial cache entry to be
* 					matched.
*/
void		AlcLRUCEntryRemoveWithKey(AlcLRUCache *cache, unsigned int key,
				          void *entry)
{
  AlcLRUCItem	*item;

  item = AlcLRUCItemFind(cache, key, entry);
  if(item)
  {
    if(cache->unlinkFn)
    {
      (*(cache->unlinkFn))(cache, item->entry);
    }
    AlcLRUCItemRankUnlink(cache, item);
    AlcLRUCItemHashUnlink(cache, item);
    --(cache->numItem);
    cache->curSz -= item->sz;
    AlcLRUCItemFree(cache, item);
  }
}

/*!
* \ingroup	AlcLRUCache
* \brief	Removes all cache entries.
* \param	cache			The cache.
*/
void            AlcLRUCEntryRemoveAll(AlcLRUCache *cache)
{
  unsigned int	idx;
  AlcLRUCItem	*item;

  item = cache->rankHead;
  while(item)
  {
    if(cache->unlinkFn)
    {
      (*(cache->unlinkFn))(cache, item->entry);
    }
    AlcLRUCItemFree(cache, item);
    item = item->rankNxt;
  }
  for(idx = 0; idx < cache->hashTblSz; ++idx)
  {
    cache->hashTbl[idx] = NULL;
  }
  cache->numItem = 0;
  cache->curSz = 0;
}

/*!
* \return	Number of items.
* \ingroup	AlcLRUCache
* \brief	Returns the number of items in the hash bin corresponding
* 		to the given item key value.
* 		This is probably only useful for debug and tuning.
* \param	cache			The cache.
* \param	key			Given item key value.
*/
unsigned int	AlcLRUCKeyGetNHashItem(AlcLRUCache *cache, unsigned int key)
{
  unsigned int  idx,
  		nItem = 0;
  AlcLRUCItem	*item;

  idx = AlcLRUCHash(key) % cache->hashTblSz;
  item = cache->hashTbl[idx];
  while(item)
  {
    ++nItem;
    item = item->hashNxt;
  }
  return(nItem);
}

/*!
* \ingroup	AlcLRUCache
* \brief	Set a new maximum total cache entry size.
* \param	cache			The cache.
* \param	newMaxSz		Given new maximum total cache
* 					entry size.
*/
void		AlcLRUCacheMaxSz(AlcLRUCache *cache, size_t newMaxSz)
{
  if(newMaxSz && (cache->curSz > newMaxSz))
  {
    AlcLRUCItem 	*item;

    while((cache->curSz > newMaxSz) && ((item = cache->rankTail) != NULL))
    {
      AlcLRUCRemoveItem(cache, item);
    }
  }
  cache->maxSz = newMaxSz;
}

/*!
* \ingroup	AlcLRUCache
* \brief	Prints a dump of the current cache status to a file.
* 		This is only intended for debug and tuning.
* \param	cache			The cache.
* \param	fP			File pointer opened for writing.
*/
void	AlcLRUCacheFacts(AlcLRUCache *cache, FILE *fP)
{
  int		idx;
  AlcLRUCItem	*item;

  (void )fprintf(fP,
                 "AlcLRUCache %p\n"
                 "  numItem =      %u\n"
                 "  maxItem =      %u\n"
                 "  curSz =        %lu\n"
                 "  maxSz =        %lu\n"
                 "  hashTblSz =    %u\n"
                 "  itemBlkSz =    %u\n"
		 "  keyFn =        %p\n"
		 "  cmpFn =        %p\n"
		 "  unlinkFn =     %p\n"
		 "  freeStack =    %p\n"
		 "  freeList =     %p\n"
		 "  rankHead =     %p\n"
		 "  rankTail =     %p\n"
		 "  hashTbl =      %p\n",
		 cache,
		 cache->numItem,
		 cache->maxItem,
		 cache->curSz,
		 cache->maxSz,
		 cache->hashTblSz,
		 cache->itemBlkSz,
		 cache->keyFn,
		 cache->cmpFn,
		 cache->unlinkFn,
		 cache->freeStack,
		 cache->freeList,
		 cache->rankHead,
		 cache->rankTail,
		 cache->hashTbl);
  (void )fprintf(fP, "rankList\n");
  item = cache->rankHead;
  while(item)
  {
    (void )fprintf(fP,
                   "  AlcLRUCItem %p\n"
		   "    key =        %u\n"
		   "    sz =         %lu\n"
		   "    entry =      %p\n"
		   "    rankPrv =    %p\n"
		   "    rankNxt =    %p\n"
		   "    hashNxt =    %p\n",
		   item,
		   item->key,
		   item->sz,
		   item->entry,
		   item->rankPrv,
		   item->rankNxt,
		   item->hashNxt);
    item = item->rankNxt;
  }
  (void )fprintf(fP, "hashTbl\n");
  for(idx = 0; idx < cache->hashTblSz; ++idx)
  {
    (void )fprintf(fP, "  % 8u\n", idx);
    item = cache->hashTbl[idx];
    while(item)
    {
      (void )fprintf(fP,
		     "  AlcLRUCItem %p\n"
		     "    key =        %u\n"
		     "    sz =         %lu\n"
		     "    entry =      %p\n"
		     "    rankPrv =    %p\n"
		     "    rankNxt =    %p\n"
		     "    hashNxt =    %p\n",
		     item,
		     item->key,
		     item->sz,
		     item->entry,
		     item->rankPrv,
		     item->rankNxt,
		     item->hashNxt);
      item = item->hashNxt;
    }
  }
  (void )fprintf(fP,"\n");
}

/*!
* \return	New cache item.
* \ingroup	AlcLRUCache
* \brief	Gets a new (or recycled) cache item. The returned cache item
* 		will have all fields cleared to zero.
* \param	cache			The cache.
*/
static AlcLRUCItem *AlcLRUCItemNew(AlcLRUCache *cache)
{
  AlcLRUCItem	*item = NULL;

  if(cache->freeList == NULL)
  {
    AlcLRUCMoreItems(cache);
  }
  if(cache->freeList != NULL)
  {
    item = cache->freeList;
    cache->freeList = item->hashNxt;
    memset(item, 0, sizeof(AlcLRUCItem));
  }
  return(item);
}

/*!
* \ingroup	AlcLRUCache
* \brief	Frees a cache item, which must previously have been unlinked,
* 		and adds it to the list for recycling.
* \param	cache			The cache.
* \param	item			Cache item to free.
*/
static void	 AlcLRUCItemFree(AlcLRUCache *cache, AlcLRUCItem *item)
{
  item->hashNxt = cache->freeList;
  cache->freeList = item;
}

static void	AlcLRUCCheckLimits(AlcLRUCache *cache, size_t newEntrySz,
				   unsigned int numNewEntry)
{
  AlcLRUCItem 	*item;

  while(((cache->numItem + numNewEntry > cache->maxItem) ||
         (cache->maxSz && (cache->curSz + newEntrySz > cache->maxSz))) &&
	((item = cache->rankTail) != NULL))
  {
    AlcLRUCRemoveItem(cache, item);
  }
}

/*!
* \ingroup	AlcLRUCache
* \brief	Increases the number of cache items available in it's
* 		free list.
* \param	cache			The cache.
*/
static void	AlcLRUCMoreItems(AlcLRUCache *cache)
{
  size_t	i,
  		n;
  AlcLRUCItem	*item;
  AlcBlockStack	*blk;

  if((blk = AlcBlockStackNew((size_t )cache->itemBlkSz, sizeof(AlcLRUCItem),
                             cache->freeStack, NULL)) != NULL)
  {
    cache->freeStack = blk;
    item = (AlcLRUCItem *)(cache->freeStack->elements);
    n = cache->freeStack->maxElm;
    for(i = 0; i < n; ++i)     
    {
      item->hashNxt = cache->freeList;
      cache->freeList = item;
      ++item;
    }
  }
}

/*!
* \ingroup	AlcLRUCache
* \brief	Unlinks the given item (which must be in the cache rank list)
* 		from the cache rank list.
* \param	cache			The cache.
* \param	item			Cache item to unlink from the rank
* 					list.
*/
static void	AlcLRUCItemRankUnlink(AlcLRUCache *cache, AlcLRUCItem *item)
{
  if(item == cache->rankTail)
  {
    if((cache->rankTail = item->rankPrv) == NULL)
    {
      cache->rankHead = NULL;
    }
    else
    {
      cache->rankTail->rankNxt = NULL;
    }
  }
  else if(item == cache->rankHead)
  {
    if((cache->rankHead = item->rankNxt) == NULL)
    {
      cache->rankTail = NULL;
    }
    else
    {
      cache->rankHead->rankPrv = NULL;
    }
  }
  else
  {
    item->rankPrv->rankNxt = item->rankNxt;
    item->rankNxt->rankPrv = item->rankPrv;
  }
}

/*!
* \ingroup	AlcLRUCache
* \brief	Inserts the given item at the head of the cache rank list.
* \param	cache			The cache.
* \param	item			Cache item to insert.
*/
static void	AlcLRUCItemRankInsert(AlcLRUCache *cache, AlcLRUCItem *item)
{
  if((item->rankNxt = cache->rankHead) == NULL)
  {
    cache->rankTail = item;
  }
  else
  {
    cache->rankHead->rankPrv = item;
  }
  item->rankPrv = NULL;
  cache->rankHead = item;
}

/*!
* \ingroup	AlcLRUCache
* \brief	Inserts the given item into the cache hash table.
* \param	cache			The cache.
* \param	item			Cache item to insert.
*/
static void	AlcLRUCItemHashInsert(AlcLRUCache *cache, AlcLRUCItem *item)
{
  unsigned int	idx;
  AlcLRUCItem	*cItem,
  		*pItem;

  idx = AlcLRUCHash(item->key) % cache->hashTblSz;
  if((cItem = cache->hashTbl[idx]) == NULL)
  {
    cache->hashTbl[idx] = item;
    item->hashNxt = NULL;
  }
  else if(item->key < cItem->key)
  {
    item->hashNxt = cItem;
    cache->hashTbl[idx] = item;
  }
  else
  {
    pItem = cItem;
    while(((cItem = cItem->hashNxt) != NULL) && (cItem->key < item->key))
    {
      pItem = pItem->hashNxt;
    }
    pItem->hashNxt = item;
    item->hashNxt = cItem;
  }
}

/*!
* \ingroup	AlcLRUCache
* \brief	Unlinks the given item from the cache hash table.
* \param	cache			The cache.
* \param	item			Cache item to unlink.
*/
static void	AlcLRUCItemHashUnlink(AlcLRUCache *cache, AlcLRUCItem *item)
{
  unsigned int  idx;
  AlcLRUCItem   *cItem,
  		*pItem;

  idx = AlcLRUCHash(item->key) % cache->hashTblSz;
  if((cItem = cache->hashTbl[idx]) != NULL)
  {
    if(item == cItem)
    {
      cache->hashTbl[idx] = item->hashNxt;
    }
    else
    {
      pItem = cItem;
      while(((cItem = cItem->hashNxt) != NULL) && (cItem != item))
      {
        pItem = pItem->hashNxt;
      }
      if(cItem == item)
      {
        pItem->hashNxt = cItem->hashNxt;
      }
    }
  }
}

/*!
* \return	Found cache item or NULL if no match found.
* \ingroup	AlcLRUCache
* \brief	Finds the cache item which matches the given entry.
* \param	cache			The cache
* \param	key			Key computed from the entry.
* \param	entry			Partical entry to match which must
* 					have sufficient information for the
* 					cache entry comparison function.
*/
AlcLRUCItem 	*AlcLRUCItemFind(AlcLRUCache *cache,
				 unsigned int key, void *entry)
{
  unsigned int	idx;
  AlcLRUCItem	*item;
  
  idx = AlcLRUCHash(key) % cache->hashTblSz;
  item = cache->hashTbl[idx];
  while(item && ((item->key != key) ||
                 (*(cache->cmpFn))(item->entry, entry)))
  {
    item = item->hashNxt;
  }
  return(item);
}

/*!
* \ingroup	AlcLRUCache
* \brief	Removes the last (least recently used) item.
* \param	cache			The cache.
* \param	item			Item to be removed.
*/
static void	AlcLRUCRemoveItem(AlcLRUCache *cache, AlcLRUCItem *item)
{
  if(cache->unlinkFn)
  {
    (*(cache->unlinkFn))(cache, item->entry);
  }
  AlcLRUCItemRankUnlink(cache, item);
  AlcLRUCItemHashUnlink(cache, item);
  --cache->numItem;
  cache->curSz -= item->sz;
  AlcLRUCItemFree(cache, item);
}

/*!
* \return	Hashed key.
* \ingroup	AlcLRUCache
* \brief	Given an unsigned integer key returns a hashed key.
* 		The hasing algorithm is Robert Jenkins 32 bit integer hash
* 		function: Bob Jenkins "Hash functions". Dr. Dobbs Journal,
* 		September 1997.
* \param	unsignedk
*/
static unsigned int AlcLRUCHash(unsigned int k)
{
  k = (k + 0x7ed55d16) + (k << 12);
  k = (k ^ 0xc761c23c) ^ (k >> 19);
  k = (k + 0x165667b1) + (k <<  5);
  k = (k + 0xd3a2646c) ^ (k <<  9);
  k = (k + 0xfd7046c5) + (k <<  3);
  k = (k ^ 0xb55a4f09) ^ (k >> 16);
  return(k);
}
