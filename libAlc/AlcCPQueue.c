#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _AlcCPQueue_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libAlc/AlcCPQueue.c
* \author       Bill Hill
* \date         July 2003
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
* \brief	An \f$O(1)\f$ priority queue based on algorithms from:
*		Brown R, 1988. Calendar Queues: A Fast {O}(1)Priority Queue
*		Implementation for the Simulation Event Set Problem.
*		Communications of the ACM 31(10), 1220-1227.
*		The priority queue has an \f$O(1)\f$ time for insertion
*		unlinking and holding of queue items with a favorable
*		constant when compared to tree based priority queue data
*		structures such as splay trees. Results presented by Brown
*		show that this data structure outperforms tree based
*		implementations and that simple ordered linked lists with
*		\f$O(n)\f$ operations outperform both calendar queues and
*		tree based queues when there are less than ten items in
*		the queue.
* \ingroup	AlcCPQ
* \todo         -
* \bug          None known.
*/
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <Alc.h>

static AlcErrno			AlcCPQMoreItems(
				  AlcCPQQueue *q);
static AlcErrno			AlcCPQQueueIncSize(
				  AlcCPQQueue *q);
static AlcErrno			AlcCPQQueueResize(
				  AlcCPQQueue *q,
				  int incFlag);
static AlcCPQItem 		*AlcCPQUnlinkItemInList(
				  AlcCPQItem **bucket);
static float			AlcCPQQueueNewBWidth(
				  AlcCPQQueue *q,
				  AlcErrno *errNum);
static void			AlcCPQQueueDecSize(
				  AlcCPQQueue *q);
static void			AlcCPQInsertItemInList(
				  AlcCPQQueue *q,
				  int gIdx,
				  AlcCPQItem *gItem);
#ifdef ALC_CPQ_DEBUG
static void			AlcCPQDebugOut(
				  AlcCPQQueue *q);
#endif /* ALC_CPQ_DEBUG */

/*!
* \return	New priority queue.
* \ingroup	AlcCPQ
* \brief	Creates a new priority queue.
* \param	dstErr			Destination error pointer may be NULL.
*/
AlcCPQQueue	*AlcCPQQueueNew(AlcErrno *dstErr)
{
  AlcCPQQueue	*q = NULL;
  AlcErrno	errNum = ALC_ER_NONE;

  if((q = (AlcCPQQueue *)AlcCalloc(1, sizeof(AlcCPQQueue))) == NULL)
  {
    errNum = ALC_ER_ALLOC;
  }
  else
  {
    q->nItem = 0;
    q->nBucket = 2;
    q->nItemIncThr = q->nBucket * 4;
    q->nItemDecThr = 0;
    q->lastIdx = 0;
    q->bucketMin = -1.0;
    q->itemBlockSz = 1024;
    q->bucketWidth = 1.0;
    q->maxBucket = 1024;
    q->bucketBase = 0;
    q->resizable = 1;
    q->lastPriority = 0.0;
    if((q->buckets = (AlcCPQItem **)
    		     AlcCalloc((size_t )(q->maxBucket),
		     	       sizeof(AlcCPQItem *))) == NULL)
    {
      errNum = ALC_ER_ALLOC;
      (void )AlcCPQQueueFree(q);
      q = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(q);
}

/*!
* \return	New priority queue item.
* \ingroup	AlcCPQ
* \brief	Creates a new priority queue item.
* \param	q			Given queue.
* \param	priority		Priority for the item.
* \param	entry			Entry for the item.
* \param	dstErr			Destination error pointer may be NULL.
*/
AlcCPQItem	*AlcCPQItemNew(AlcCPQQueue *q, float priority,
			      void *entry, AlcErrno *dstErr)
{
  AlcCPQItem	*item = NULL;
  AlcErrno	errNum = ALC_ER_NONE;

  if(priority < 0.0)
  {
    errNum = ALC_ER_PARAM;
  }
  else
  {
    if(q->freeItem == NULL)
    {
      errNum = AlcCPQMoreItems(q);
    }
  }
  if(errNum == ALC_ER_NONE)
  {
    item = q->freeItem;
    q->freeItem = item->next;
    if(q->freeItem)
    {
      q->freeItem->prev = NULL;
    }
    item->prev = item->next = NULL;
    item->priority = priority;
    item->entry = entry;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(item);
}

/*!
* \return	Alc error code.
* \ingroup	AlcCPQ
* \brief	Frees a priority queue. Queue item entries are not
*		freed.
* \param	q			Given queue to free.
*/
AlcErrno	AlcCPQQueueFree(AlcCPQQueue *q)
{
  AlcErrno	errNum = ALC_ER_NONE;

  if(q)
  {
    AlcFree(q->buckets);
    (void )AlcBlockStackFree(q->freeStack);
    AlcFree(q);
  }
  return(errNum);
}

/*!
* \return	Alc error code.
* \ingroup	AlcCPQ
* \brief	Frees a priority queue item by putting it onto the freeItem
*		list. The queue items entries are not freed.
* \param	q			The queue.
* \param	item			Given item to free.
*/
void		AlcCPQItemFree(AlcCPQQueue *q, AlcCPQItem *item)
{
  if(q && item)
  {
    item->next = NULL;
    if((item->prev = q->freeItem) != NULL)
    {
      item->prev->next = item;
    }
    q->freeItem = item;
  }
}

/*!
* \return	Alc error code.
* \ingroup	AlcCPQ
* \brief	Inserts an entry into the priority queue. This function
*		calls AlcCPQItemNew() followed by AlcCPQItemInsert().
* \param	q			The queue.
* \param	priority		Priority of entry.
* \param	entry			Given entry.
*/
AlcErrno	AlcCPQEntryInsert(AlcCPQQueue *q, float priority, void *entry)
{
  AlcCPQItem	*item;
  AlcErrno	errNum = ALC_ER_NONE;

  if((item = AlcCPQItemNew(q, priority, entry, &errNum)) != NULL)
  {
    errNum = AlcCPQItemInsert(q, item);
  }
  return(errNum);
}

/*!
* \return	Alc error code.
* \ingroup	AlcCPQ
* \brief	Inserts an item into the priority queue.
* \param	q			The queue.
* \param	item			Given item to insert.
*/
AlcErrno	AlcCPQItemInsert(AlcCPQQueue *q, AlcCPQItem *item)
{
  int		idx;
  AlcErrno	errNum = ALC_ER_NONE;

#ifdef ALC_CPQ_DEBUG
  (void )fprintf(stderr,
		 "AlcCPQItemInsert() 0x%lx 0x%lx %g\n",
		 (unsigned long )q, (unsigned long)item, item->priority);
  AlcCPQDebugOut(q);
#endif /* ALC_CPQ_DEBUG */
  if((q->nItem >= q->nItemIncThr) && (q->resizable != 0))
  {
    errNum = AlcCPQQueueIncSize(q);
  }
  if(errNum == ALC_ER_NONE)
  {
    idx = q->nBucket -
          (int )floor(item->priority / q->bucketWidth) % q->nBucket - 1;
    AlcCPQInsertItemInList(q, idx, item);
    ++(q->nItem);
    if(item->priority > q->lastPriority)
    {
      q->lastIdx = idx;
      q->lastPriority = item->priority;
      q->bucketMin = floor(item->priority / q->bucketWidth) * q->bucketWidth;
    }
  }
#ifdef ALC_CPQ_DEBUG
  AlcCPQDebugOut(q);
  (void )fprintf(stderr,
		 "AlcCPQItemInsert() %d\n",
		 (int )errNum);
#endif /* ALC_CPQ_DEBUG */
  return(errNum);
}

/*!
* \return	The highest priority item or NULL if the queue is empty.
* \ingroup	AlcCPQ
* \brief	Unlinks and returns the item with the highest priority from
*		the given priority queue.
* \param	q			Given queue.
*/
AlcCPQItem	*AlcCPQItemUnlink(AlcCPQQueue *q)
{
  int		idx,
  		idN,
		idM;
  AlcCPQItem	*item0 = NULL,
  		*item1;

  if(q->nItem > 0)
  {
    /* Search the buckets for the highest priority item . */
    idM = -1;
    if(q->bucketMin > 0.0)
    {
      idx = q->lastIdx;
      do
      {
	idN = (idx + 1) % q->nBucket;
	if((item0 = *(q->buckets + q->bucketBase + idx)) != NULL)
	{
	  if(item0->priority >= q->bucketMin)
	  {
	    idM = idx;
	    break;
	  }
	}
	idx = idN;
	q->bucketMin -= q->bucketWidth;
      } while((idx != q->lastIdx) && (q->bucketMin > 0.0));
    }
    if(idM < 0)
    {
      /* Do direct search of all buckets. */
      idx = 0;
      do
      {
	idM = idx;
	item0 = *(q->buckets + q->bucketBase + idx);
      } while((item0 == NULL) && (++idx < q->nBucket));
      if(item0)
      {
	while(++idx < q->nBucket) 
	{
	  if(((item1 = *(q->buckets + q->bucketBase + idx)) != NULL) &&
	     (item1->priority > item0->priority))
	  {
	    idM = idx;
	    item0 = item1;
	  }
	}
      }
    }
    /* Unlink the item. */
    item0 = (idM < q->nBucket)?
	    AlcCPQUnlinkItemInList(q->buckets + q->bucketBase + idM): NULL;
    /* Set last bucket index and lack bucket minimum. */
    if(item0)
    {
#ifdef ALC_CPQ_DEBUG
      (void )fprintf(stderr,
      		     "AlcCPQItemUnlink() %g %4d %4d %4d\n",
		     item0->priority, q->nItem, q->lastIdx, idM);
#endif /* ALC_CPQ_DEBUG */
      q->lastIdx = idM;
      q->lastPriority = item0->priority;
      q->bucketMin = floor(item0->priority / q->bucketWidth) * q->bucketWidth;
      --(q->nItem);
      if((q->nItem < q->nItemDecThr) && (q->resizable != 0))
      {
	AlcCPQQueueDecSize(q);
      }
    }
  }
  return(item0);
}

/*!
* \return	Alc error code.
* \ingroup	AlcCPQ
* \brief	Allocate more items and add them to the free item list
*		ready for use.
* \param	q			The queue.
*/
static AlcErrno	AlcCPQMoreItems(AlcCPQQueue *q)
{
  int		idx,
  		iCnt;
  AlcCPQItem	*lItem,
  		*mItem;
  AlcErrno	errNum = ALC_ER_NONE;

  if((q->freeStack = AlcBlockStackNew((size_t )q->itemBlockSz,
  				      sizeof(AlcCPQItem),
       				      q->freeStack, NULL)) == NULL)
  {
    errNum = ALC_ER_ALLOC;
  }
  else
  {
    mItem = lItem = (AlcCPQItem *)(q->freeStack->elements);
    mItem->next = q->freeItem;
    if(q->freeItem)
    {
      q->freeItem->prev = mItem;
    }
    q->freeItem = mItem;
    ++mItem;
    iCnt = (int )(q->freeStack->maxElm) - 1;
    for(idx = 1; idx < iCnt; ++idx)
    {
      mItem->prev = lItem;
      lItem->next = mItem;
      lItem = mItem;
      ++mItem;
    }
    mItem->prev = lItem;
    mItem->next = NULL;
  }
  return(errNum);
}

/*!
* \return
* \ingroup      AlcCPQ
* \brief	Insert the given item into the priority queue in the
*		given bucket.
* \param	q			The queue.
* \param	gIdx			Given bucket index.
* \param	gItem			Given item to insert.
*/
static void	AlcCPQInsertItemInList(AlcCPQQueue *q, int gIdx,
				      AlcCPQItem *gItem)
{
  AlcCPQItem	*item0,
  		*item1;
  AlcCPQItem	**bucket;

  bucket = q->buckets + q->bucketBase + gIdx;
  if(*bucket == NULL)
  {
    *bucket = gItem;
    gItem->prev = gItem->next = NULL;
  }
  else
  {
    item0 = NULL;
    item1 = *bucket;
    while(item1 && (item1->priority > gItem->priority))
    {
      item0 = item1;
      item1 = item1->next;
    }
    if(item0 == NULL)
    {
      *bucket = gItem;
      item1->prev = gItem;
      gItem->next = item1;
      gItem->prev = NULL;
    }
    else if(item1 == NULL)
    {
      item0->next = gItem;
      gItem->prev = item0;
      gItem->next = NULL;
    }
    else
    {
      item0->next = gItem;
      item1->prev = gItem;
      gItem->prev = item0;
      gItem->next = item1;
    }
  }
}

/*!
* \return	Unlinked item.
* \ingroup	AlcCPQ
* \brief	Unlink and return the highest priority item from the
*		given bucket.
*		
* \param	bucket			Bucket with item at it's head.
*/
static AlcCPQItem *AlcCPQUnlinkItemInList(AlcCPQItem **bucket)
{
  AlcCPQItem	*item;

  if((item = *bucket) != NULL)
  {
    if((*bucket = item->next) != NULL)
    {
      item->next->prev = NULL;
    }
  }
  return(item);
}

/*!
* \return	Alc error code.
* \ingroup	AlcCPQ
* \brief	Increase the queue size by a factor of two.
* \param	q			The queue.
*/
static AlcErrno	AlcCPQQueueIncSize(AlcCPQQueue *q)
{
  AlcErrno	errNum = ALC_ER_NONE;

  if(q->resizable)
  {
    q->resizable = 0;
    if(q->maxBucket < q->nBucket * 3)
    {
      q->maxBucket *= 4;
      if((q->buckets = (AlcCPQItem **)
      		       AlcRealloc(q->buckets,
		       		  q->maxBucket * sizeof(AlcCPQItem *))) == NULL)
      {
        errNum = ALC_ER_ALLOC;
      }
    }
    if(errNum == ALC_ER_NONE)
    {
      errNum = AlcCPQQueueResize(q, 1);
    }
    q->resizable = 1;
  }
  return(errNum);
}

/*!
* \return	void
* \ingroup	AlcCPQ
* \brief	Decrease the queue size by a factor of two.
* \param	q			The queue.
*/
static void	AlcCPQQueueDecSize(AlcCPQQueue *q)
{
  if(q->resizable)
  {
    q->resizable = 0;
    (void )AlcCPQQueueResize(q, 0);
    q->resizable = 1;
  }
}

/*!
* \return
* \ingroup      AlcCPQ
* \brief	Resizes the queue, either increasing or decreasing it's
*		number of buckets by a factor of two, except when the
*		number of buckets is small.
* \param	q			The queue.
* \param	incFlag			Increase the number of buckets if
*					non-zero otherwise decrease the number.
*/
static AlcErrno	AlcCPQQueueResize(AlcCPQQueue *q, int incFlag)
{
  int		idx,
  		oldNBuckets,
  		oldBucketBase;
  AlcCPQItem	*oldItem,
  		*oldNextItem;
  AlcErrno	errNum = ALC_ER_NONE;

#ifdef ALC_CPQ_DEBUG
  (void )fprintf(stderr, "AlcCPQQueueResize() FE\n");
  AlcCPQDebugOut(q);
#endif /* ALC_CPQ_DEBUG */
  oldNBuckets = q->nBucket;
  oldBucketBase = q->bucketBase;
  q->bucketWidth = AlcCPQQueueNewBWidth(q, &errNum);
  if(errNum == ALC_ER_NONE)
  {
    if(incFlag)
    {
      q->nBucket = q->nBucket * 2;
      q->nItemIncThr = q->nBucket * 2;
      if(q->nItemIncThr == 64)
      {
	q->nItemDecThr = 16;
      }
      else if(q->nItemDecThr > 0)
      {
	q->nItemDecThr *= 2;
      }
    }
    else
    {
      q->nBucket = q->nBucket / 2;
      q->nItemIncThr = q->nBucket / 2;
      if(q->nItemIncThr < 64)
      {
	q->nItemDecThr = 0;
      }
      else
      {
	q->nItemDecThr /= 2;
      }
    }
    q->nItem = 0;
    q->lastIdx = 0;
    q->lastPriority = 0.0;
    q->bucketBase = (oldBucketBase)? 0: q->maxBucket - q->nBucket;
    q->bucketMin = -1.0;
    for(idx = 0; idx < q->nBucket; ++idx)
    {
      *(q->buckets + q->bucketBase + idx) = NULL;
    }
  }
  idx = 0;
  while((errNum == ALC_ER_NONE) && (idx < oldNBuckets))
  {
    oldItem = *(q->buckets + oldBucketBase + idx);
    while((errNum == ALC_ER_NONE) && oldItem)
    {
      oldNextItem = oldItem->next;
      errNum = AlcCPQItemInsert(q, oldItem);
      oldItem = oldNextItem;
    }
    ++idx;
  }
#ifdef ALC_CPQ_DEBUG
  if(errNum == ALC_ER_NONE)
  {
    AlcCPQDebugOut(q);
  }
  (void )fprintf(stderr, "AlcCPQQueueResize() FX\n");
#endif /* ALC_CPQ_DEBUG */
  return(errNum);
}

/*!
* \return	New band width.
* \ingroup	AlcCPQ
* \brief	Compute the new priority queue band width using the algorithm
* 		of Brown:1988.
* \param	q			The queue.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static float	AlcCPQQueueNewBWidth(AlcCPQQueue *q, AlcErrno *dstErr)
{
  int		idx,
  		nSam0,
		nSam1;
  double	tD0,
  		mean,
  		sum,
		newBWidth;
  AlcErrno	errNum = ALC_ER_NONE;
  AlcCPQItem	*items[25];

  if(q->nItem < 2)
  {
    newBWidth = 1.0;
  }
  else
  {
    if(q->nItem <= 5)
    {
      nSam0 = q->nItem;
    }
    else
    {
      nSam0 = 5 + (q->nItem / 10);
    }
    if(nSam0 > 25)
    {
      nSam0 = 25;
    }
    /* Unlink nSam0 items from the queue and store them locally while
     * recording the priorities, entries and free functions. */
#ifdef ALC_CPQ_DEBUG
    (void )fprintf(stderr,
    		   "AlcCPQQueueNewBWidth() Unlinked item priorities:\n");
#endif /* ALC_CPQ_DEBUG */
    for(idx = 0; idx < nSam0; ++idx)
    {
      items[idx] = AlcCPQItemUnlink(q);
#ifdef ALC_CPQ_DEBUG
      if(items[idx])
      {
        (void )fprintf(stderr, "  %g\n", items[idx]->priority);
      }
      else
      {
        (void )fprintf(stderr, "  Item NULL\n");
      }
#endif /* ALC_CPQ_DEBUG */
    }
    /* Compute the mean priority seperation. */
    mean = (items[0]->priority - items[nSam0 - 1]->priority) / (nSam0 - 1);
    /* Recompute the mean priority seperation using only those that are
     * < twice the original mean. */
    nSam1 = 0;
    sum = 0.0;
    mean *= 2.0;
    for(idx = 1; idx < nSam0; ++idx)
    {
      tD0 = items[idx - 1]->priority - items[idx]->priority;
      if(tD0 <= mean)
      {
	++nSam1;
	sum += tD0;
      }
    }
    /* Set the new band width to 3 * this new mean. */
    newBWidth = (nSam1 > 0)? (3.0 * sum) / (nSam1 + 1): 1.0;
    /* Reinsert the items and restore the last bucket index. */
    idx = 0;
    while((errNum == ALC_ER_NONE) && (idx < nSam0))
    {
      errNum = AlcCPQItemInsert(q, items[idx]);
      ++idx;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(newBWidth);
}

#ifdef ALC_CPQ_DEBUG
static void	AlcCPQDebugOut(AlcCPQQueue *q)
{
  int		idB;
  AlcCPQItem	*item;

  (void )fprintf(stderr, "q:\n");
  (void )fprintf(stderr, "  nItem = %d\n", q->nItem);
  (void )fprintf(stderr, "  nBucket = %d\n", q->nBucket);
  (void )fprintf(stderr, "  nItemIncThr = %d\n", q->nItemIncThr);
  (void )fprintf(stderr, "  nItemDecThr = %d\n", q->nItemDecThr);
  (void )fprintf(stderr, "  lastIdx = %d\n", q->lastIdx);
  (void )fprintf(stderr, "  itemBlockSz = %d\n", q->itemBlockSz);
  (void )fprintf(stderr, "  bucketBase = %d\n", q->bucketBase);
  (void )fprintf(stderr, "  maxBucket = %d\n", q->maxBucket);
  (void )fprintf(stderr, "  resizable = %s\n", (q->resizable)? "True":"False");
  (void )fprintf(stderr, "  bucketWidth = %g\n", q->bucketWidth);
  (void )fprintf(stderr, "  bucketMin = %g\n", q->bucketMin);
  (void )fprintf(stderr, "  lastPriority = %g\n", q->lastPriority);
  for(idB = 0; idB < q->nBucket; ++idB)
  {
    (void )fprintf(stderr, "  bucket[%4d]: ", idB);
    item = *(q->buckets + q->bucketBase + idB);
    while(item)
    {
      (void )fprintf(stderr, "%g, ", item->priority);
      item = item->next;
    }
    (void )fprintf(stderr, "NULL\n");
  }
}
#endif /* ALC_CPQ_DEBUG */

#if (ALC_CPQ_TEST_MAIN == 1)
/*
* Test #1
*
* This test main builds a priority queue and then knocks it down again. 
* At each stage facts about the queue are output.
*/
int		main(int argc, char *argv[])
{
  int		idI,
  		idR,
		nRepeats = 1,
		verbose = 1;
  AlcCPQQueue	*q;
  AlcCPQItem	*item;
  float		priority[4000];
  AlcErrno	errNum = ALC_ER_NONE;
  const int	nItems = 4000;

  if(verbose)
  {
    (void )fprintf(stderr, "Building queue with priorities:\n");
  }
  for(idI = 0; idI < nItems; ++idI)
  {
    priority[idI] = ((double )rand()) / RAND_MAX;
    if(verbose)
    {
      (void )fprintf(stderr, "%g\n", priority[idI]);
    }
  }
  if((q = AlcCPQQueueNew(&errNum)) == NULL)
  {
    (void )fprintf(stderr, "%s: Failed to create a priority queue.\n", *argv);
    exit(1);
  }
  for(idR = 0; (errNum == ALC_ER_NONE) && (idR < nRepeats); ++idR)
  {
    for(idI = 0; (errNum == ALC_ER_NONE) && (idI < nItems); ++idI)
    {
      errNum = AlcCPQEntryInsert(q, priority[idI], NULL);
#ifdef ALC_CPQ_DEBUG
      if(verbose)
      {
        AlcCPQDebugOut(q);
      }
#endif /* ALC_CPQ_DEBUG */
    }
    if(errNum != ALC_ER_NONE)
    {
      (void )fprintf(stderr, "%s: Failed to enqueue item %d of %d\n",
      		     *argv, idI, nItems);
      exit(1);
    }
    if((errNum == ALC_ER_NONE) && verbose)
    {
      (void )fprintf(stderr, "Removing items from queue with priorities:\n");
    }
    idI = 0;
    do
    {
#ifdef ALC_CPQ_DEBUG
      if(verbose)
      {
        AlcCPQDebugOut(q);
      }
#endif /* ALC_CPQ_DEBUG */
      item = AlcCPQItemUnlink(q);
#ifdef ALC_CPQ_DEBUG
      if(verbose)
      {
        AlcCPQDebugOut(q);
      }
#endif /* ALC_CPQ_DEBUG */
      if(verbose)
      {
	if(item)
	{
          (void )fprintf(stderr, "%g\n", item->priority);
	}
	else
	{
	  (void )fprintf(stderr, "Item NULL\n");
	}
      }
      AlcCPQItemFree(q, item);
      ++idI;
    } while(item != NULL);
    (void )AlcCPQQueueFree(q);
    if(errNum != ALC_ER_NONE)
    {
      (void )fprintf(stderr, "%s: Failed to unlink item %d of %d repeat %d\n",
      		     *argv, idI, nItems, idR);
      exit(1);
    }
  }
  return(0);
}
#endif /* ALC_CPQ_TEST_MAIN == 1 */

#if (ALC_CPQ_TEST_MAIN == 2)

#include <sys/time.h>

/*
* Test #2
*
* This test main builds a priority queue and establishes the hold time
* at each stage of the build. At each stage the queue size and hold time
* (\f$\mu\f$s) are output.
*/
int		main(int argc, char *argv[])
{
  int		idI,
		idP,
  		idR;
  long		time0,
  		time1;
  double	timeDelta;
  AlcCPQQueue	*q;
  AlcCPQItem	*item;
  float		*priority;
  struct timeval tp;
  AlcErrno	errNum = ALC_ER_NONE;
  const int	nRepeats = 100,
  		maxItems = 100,
  		maxPriorites = 100;

  if((priority = (float *)AlcMalloc(sizeof(float) * maxPriorites)) == NULL)
  {
    (void )fprintf(stderr, "%s: Failed to allocate memory.\n", *argv);
    exit(1);
  }
  for(idP = 0; idP < maxPriorites; ++idP)
  {
    *(priority + idP) = ((double )rand()) / RAND_MAX;
  }
  if((q = AlcCPQQueueNew(&errNum)) == NULL)
  {
    (void )fprintf(stderr, "%s: Failed to create a priority queue.\n", *argv);
    exit(1);
  }
  idP = 0;
  for(idI = 0; (errNum == ALC_ER_NONE) && (idI < maxItems); ++idI)
  {
    errNum = AlcCPQEntryInsert(q, *(priority + idP), NULL);
#ifdef ALC_CPQ_DEBUG
    AlcCPQDebugOut(q);
#endif /* ALC_CPQ_DEBUG */
    (void )gettimeofday(&tp, NULL);
    for(idR = 0; (errNum == ALC_ER_NONE) && (idR < nRepeats); ++idR)
    {
      item = AlcCPQItemUnlink(q);
      idP = (idP + 1) % maxPriorites;
      item->priority = *(priority + idP);
      errNum = AlcCPQItemInsert(q, item);
    }
    time0 = (tp.tv_sec * 1000000) + tp.tv_usec;
    (void )gettimeofday(&tp, NULL);
    time1 = (tp.tv_sec * 1000000) + tp.tv_usec;
    timeDelta = (double )(time1 - time0) / nRepeats;
    (void )printf("AlcCPQueue hold time for queue size %d = %lgus\n",
                  idI, timeDelta);
  }
  if(errNum == ALC_ER_NONE)
  {
    for(idI = 0; (errNum == ALC_ER_NONE) && (idI < maxItems); ++idI)
    {
      errNum = AlcCPQEntryInsert(q, *(priority + idP), NULL);
    }
    time0 = (tp.tv_sec * 1000000) + tp.tv_usec;
    while((item = AlcCPQItemUnlink(q)) != NULL)
    {
#ifdef ALC_CPQ_DEBUG
      AlcCPQDebugOut(q);
#endif /* ALC_CPQ_DEBUG */
    }
    idI = q->nItem;
    time1 = (tp.tv_sec * 1000000) + tp.tv_usec;
    timeDelta = (double )(time1 - time0);
    (void )printf("AlcCPQueue total unlink time for queue size %d = %lgus\n",
		  idI, timeDelta);
  }
  (void )AlcCPQQueueFree(q);
  if(errNum != ALC_ER_NONE)
  {
    (void )fprintf(stderr, "%s: Test failed (%d)\n",
		   *argv, (int )errNum);
    exit(1);
  }
  return(0);
}
#endif /* ALC_CPQ_TEST_MAIN == 2 */
