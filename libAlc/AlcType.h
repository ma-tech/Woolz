#ifndef ALCTYPE_H
#define ALCTYPE_H
#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlcType_h[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlc/AlcType.h
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
* \brief        Type definitions for the Woolz type allocation
*		library.
*/

#ifndef WLZ_EXT_BIND
#ifdef __cplusplus
extern "C" {
#endif
#endif /* WLZ_EXT_BIND */

/************************************************************************
* Portability macros.                                                   *
************************************************************************/
#ifdef HAVE_STRTOK_R
#define ALC_STRTOK_R(S,D,P)     strtok_r((S),(D),(P))
#else
#define ALC_STRTOK_R(S,D,P)     strtok((S),(D))
#endif
#ifdef HAVE_TIMERSUB
#define ALC_TIMERSUB(A,B,R)	timersub((A),(B),(R))
#else
#define ALC_TIMERSUB(A,B,R) \
{ \
  (R)->tv_sec = (A)->tv_sec - (B)->tv_sec; \
  (R)->tv_usec = (A)->tv_usec - (B)->tv_usec; \
  if((R)->tv_usec < 0) \
  {  \
    (R)->tv_sec  -= 1; \
    (R)->tv_usec += 1000000; \
  } \
}
#endif


/*!
* \enum 	_AlcErrno
* \ingroup	AlcError
* \brief	Error codes returned by functions of the
*		memory allocation and fundamental type library.
*               Typedef: ::AlcErrno
*/
typedef enum _AlcErrno
{
  ALC_ER_NONE = 0,		/*!< Noe error */
  ALC_ER_ALLOC,			/*!< Memory allocation error */
  ALC_ER_NULLPTR,		/*!< Null pointer detected */
  ALC_ER_NUMELEM,		/*!< Inappropriate number of elements */
  ALC_ER_PARAM,			/*!< Inappropriate function paramerter */
  ALC_ER_READ,			/*!< File read failure */
  ALC_ER_WRITE			/*!< File write failure */
} AlcErrno;


/*!
* \struct	_AlcBlockStack
* \ingroup	AlcBlockStack
* \brief	General purpose data structure for maintaining blocks
*		of some data type. Useful for efficient memory allocation.
*		It's not a stack but a doubly linked list of blocks of
*		data which can be used as a stack, heap, list, ....
*               Typedef: ::AlcBlockStack
*/
typedef struct _AlcBlockStack
{
  size_t           elmCnt;	/*!< Number of elements used in block */
  size_t           maxElm;  	/*!< Number of elements space allocated for
  				     in block */
  void          *elements;      /*!< Block of elements */
  struct _AlcBlockStack *prev;  /*!< Previous block up in stack */
  struct _AlcBlockStack *next;  /*!< Next block down, remember it's a
  				     doubly linked list not a stack! */
} AlcBlockStack;

#ifndef WLZ_EXT_BIND
/*!
* \typedef	AlcLRUCKeyFn
* \ingroup	AlcLRUCache
* \brief	Function called to compute a ::AlcLRUCache item's numeric
* 		key given it's entry.
* 		The required function parameters are the cache and the
* 		entry.
*/
struct _AlcLRUCache;
typedef unsigned int (*AlcLRUCKeyFn)(struct _AlcLRUCache *, void *);

/*!
* \typedef	AlcLRUCCmpFn
* \ingroup	AlcLRUCache
* \brief	Function called to compare two entries in a ::AlcLRUCache.
* 		This function must return zero only iff the two entries
* 		match.
* 		The required function parameters are the two cache entry
* 		data structures with sufficient fields set to make a
* 		comparison.
*/
typedef int	(*AlcLRUCCmpFn)(const void *, const void *);

/*!
* \typedef	AlcLRUCCmpFn
* \ingroup	AlcLRUCache
* \brief	Function called when a ::AlcLRUCache item is unlinked and
* 		removed from the cache. It may be used, for example, to
* 		free the item's entry.
* 		The required function parameters are the cache and the
* 		entry that will be unlinked and removed.
*/
struct _AlcLRUCache;
struct _AlcLRUCItem;
typedef void	(*AlcLRUCUnlinkFn)(struct _AlcLRUCache *, void *);

/*!
* \struct	_AlcLRUCItem
* \ingroup	AlcLRUCache
* \brief	A cache item for a ::AlcLRUCache.
* 		Typedef: ::AlcLRUCItem
*/
typedef struct	_AlcLRUCItem
{
  unsigned int	key;		/*!< Numeric key used as input to a hash
  				     function to locate the item. */
  size_t	sz;		/*!< Size of the item's entry. */
  void          *entry;		/*!< User supplied entry. */
  struct _AlcLRUCItem *rankPrv; /*!< Previous item in priority rank table. */
  struct _AlcLRUCItem *rankNxt; /*!< Next item in priority rank table. */
  struct _AlcLRUCItem *hashNxt; /*!< Next item in hash table or free
      				     item list. */
} AlcLRUCItem;

/*!
* \struct	_AlcLRUCache
* \ingroup	AlcLRUCache
* \brief	A least recent use removal cache allowing rank and random
* 		access to it's entries. Rank access is via a doubly linked
* 		list (using rankPrv and rankNxt), while random access is
* 		via hash table (using the item's key and hashNxt).
* 		The cache will unlink items as required to maintain the
* 		set maximum number of entries and maximum total entry
* 		size.
* 		Typedef: ::AlcLRUCache
*/
typedef struct	_AlcLRUCache
{
  unsigned int	numItem;	/*!< Current number of items in the cache. */
  unsigned int	maxItem;	/*!< Maximum number of items in the cache. */
  size_t	curSz;		/*!< Current total cache size. */
  size_t	maxSz;		/*!< Maximum total cache size. */
  unsigned int	hashTblSz;	/*!< Size of the hash table. */
  unsigned int  itemBlkSz;	/*!< Number of items to allocate in a block. */
  AlcLRUCKeyFn	keyFn;		/*!< Function called to compute a numeric
  				     key for an entry. */
  AlcLRUCCmpFn	cmpFn;		/*!< Function called to compare two entries,
  				     returning zero only iff the entries
				     match. */
  AlcLRUCUnlinkFn unlinkFn;	/*!< Function called before unlinking and
  				     removing an item from the cache. */
  struct _AlcBlockStack *freeStack; /*!< Stack of blocks of items. */
  struct _AlcLRUCItem *freeList; /*!< Head of the linked list of free cache
  				      items available for use. */
  struct _AlcLRUCItem *rankHead; /*!< Head of the linked list of cache items
  				      in priority rank order. This is
				      maintained with the most recently used
				      item at the head and the least recently
				      used item at the tail. */
  struct _AlcLRUCItem *rankTail; /*!< Tail of the linked list of cache items
  				      in rank order. */
  struct _AlcLRUCItem **hashTbl; /*!< Hash table of cache items. */
} AlcLRUCache;
#endif /* WLZ_EXT_BIND */

/*!
* \struct	_AlcCPQQueue
* \ingroup	AlcCPQ
* \brief	An \f$O(1)\f$ priority queue based on the Calendar
*		Priority Queue data structure.
*		Typedef: ::AlcCPQQueue.
*/
typedef struct _AlcCPQQueue
{
  int           nItem;		/*!< Number of items in the queue. */
  int           nBucket;	/*!< Number of buckets for items. */
  int           nItemIncThr;	/*!< Threshold number of items, at
  				     which the number of buckets is
				     increased. */
  int           nItemDecThr;	/*!< Threshold number of items, at
  				     which the number of buckets is
				     decreased. */
  int           lastIdx;	/*!< The index (offset by bucketBase) of
                                     the last bucket from which an item
				     was dequeued or the index of a top
				     priority item that has been inserted.
				     A negative value is used to indicate
				     that the index is invalid. */
  int           itemBlockSz;	/*!< Number of items allocated in each
  				     block by AlcBlockStackNew(). */
  int           bucketBase;	/*!< The offset into the allocated buckets
  				     of the current buckets. */
  int           maxBucket;	/*!< The total number of buckets that have
  				     been allocated. */
  int           resizable;	/*!< Non-zero if the queue is resizable. */
  float         lastPriority;	/*!< The priority of either the last item
  				     dequeued or of the top priority item
				     that has been inserted. */
  double        bucketWidth;	/*!< The width of each bucket. */
  double        bucketMin;	/*!< The minimum priority for the bucket
  				     with (offset) index lastIdx. A negative
				     value is used to indicate that the
				     minimum priority is invalid. */
  struct _AlcCPQItem **buckets;	/*!< Contiguous array of buckets. */
  struct _AlcCPQItem *freeItem;	/*!< List of items available for use in the
  				     priority queue. */
  struct _AlcBlockStack *freeStack; /*!< Free stack used to efficiently
  				     allocate and free blocks of queue
				     items. */
} AlcCPQQueue;

/*!
* \struct	_AlcCPQItem
* \ingroup	AlcCPQ
* \brief	An item in a calendar priority queue.
*		Typedef: ::AlcCPQItem.
*/
typedef struct _AlcCPQItem
{
  float         priority;	/*!< Priority of the item, which must be
  				     \f$ \geq 0\f$. Priority is greater
				     for greater priority values. */
  void          *entry;		/*!< User supplied entry. May be NULL. */
  struct _AlcCPQItem *prev;	/*!< Previous item in bucket or free
      				     item list. */
  struct _AlcCPQItem *next;	/*!< Next item in bucket or free
      				     item list. */
} AlcCPQItem;

/*!
* \enum 	_AlcDirection
* \ingroup	AlcBlockStack
* \brief	Data structure traversal direction.
*               Typedef: ::AlcDirection
*/
typedef enum _AlcDirection
{
  ALC_DIRECTION_FWD,	        /*!< Towards end of a data structure,
  				     eg the tail of a list */
  ALC_DIRECTION_REV		/*!< Towards begining of a data structure,
  				     eg the head of a list */
} AlcDirection;


/*!
* \struct	_AlcDLPItem
* \ingroup	AlcDLPList
* \brief	A doubly linked list item.
*               Typedef: ::AlcDLPItem
*/
typedef struct _AlcDLPItem
{
#ifdef WLZ_EXT_BIND
  void          *freeFn;
#else /* WLZ_EXT_BIND */
  void          (*freeFn)(void *); /*!< Function to free list item, may be
  					NULL */
#endif /* WLZ_EXT_BIND */
  void          *entry;	    	/*!< The list item's entry */
  struct _AlcDLPItem *next;	/*!< The next item in the list, towards the
  				     tail of the list*/
  struct _AlcDLPItem *prev;     /*!< The previous item in the list, towards
  				     the head of the list */
} AlcDLPItem;

/*!
* \struct	_AlcDLPList
* \ingroup	AlcDLPList
* \brief	A doubly linked list of pointers.
*               Typedef: ::AlcDLPList
*/
typedef struct _AlcDLPList
{
  AlcDLPItem *head;	    	/*!< The head of the list */
} AlcDLPList;

#ifndef WLZ_EXT_BIND
/*!
* \struct	_AlcHashItem
* \ingroup	AlcHashTable
* \brief	A hash table item.
*               Typedef: ::AlcHashItem
*/
typedef struct _AlcHashItem
{
  void		(*freeFn)(void *); /*!< Function to free the hash item, may
  					be NULL */
  void		*entry;	    	/*!< The hash item's entry */
  void		*key;		/*!< The key which identifies the hash item */
  struct _AlcHashItem *next;	/*!< Next item in the hash tables linked list,
  				     NULL if this is the last item */
  struct _AlcHashItem *prev; 	/*!< Previous item in the hash tables linked
  				     list, NULL if this is the first item */
} AlcHashItem;

/*!
* \struct	_AlcHashTable
* \ingroup	AlcHashTable
* \brief	A hash table.
*               Typedef: ::AlcHashTable
*/
typedef struct _AlcHashTable
{
  int		(*keyCmp)(void *, void *); /*!< Key comparison function */
  size_t	(*hashFn)(void *);         /*!< The hashing function */
  size_t	tableSz;		   /*!< Number of list slots in the
  						hash table */
  AlcHashItem	**table;		   /*!< Table of lists */
} AlcHashTable;
#endif /* WLZ_EXT_BIND */


/*!
* \struct	_AlcVector
* \ingroup      AlcVector
* \brief	An extensible 1D array.
*               Typedef: ::AlcVector
*/
typedef struct _AlcVector
{
  size_t	elmSz;		/*!< Size of elements of the vector */
  size_t	blkCnt;		/*!< Number of block pointers */
  size_t	blkUse;		/*!< Number of blocks used */
  size_t	blkSz;		/*!< Number of elements in a block, must NOT be
  				     changed once vector has been created! */
  void		*freeStack;	/*!< Free stack */
  void		**blocks;	/*!< Data blocks */
} AlcVector;

/*!
* \enum		_AlcPointType
* \brief	Type of coordinate point.
*               Typedef: ::AlcPointType
*/

typedef enum _AlcPointType
{
  ALC_POINTTYPE_INT,	    	/*!< Integer point */
  ALC_POINTTYPE_DBL		/*!< Double precision floating point point */
} AlcPointType;

/*!
* \union	_AlcPointP
* \brief	Pointer to a generic coordinate.
*               Typedef: ::AlcPointP
*/
typedef union _AlcPointP
{
  void		*kV;   		/*!< Pointer to any point */
  int		*kI;		/*!< Pointer to an integer point */
  double	*kD;	     	/*!< Pointer to a double precision floating
  				     point point */
} AlcPointP;


/*!
* \struct	_AlcKDTNode
* \ingroup	AlcKDTree
* \brief	A node in a binary space partition tree (kD-tree).
*               Typedef: ::AlcKDTNode
*/
typedef struct _AlcKDTNode
{
  size_t	idx;		/*!< Index or identifier for the node */
  int		split;  	/*!< The splitting dimension */
  struct _AlcKDTNode *parent; 	/*!< The parent node, NULL if the node is
  				     the root of the tree */
  struct _AlcKDTNode *childN;  	/*!< Child node with -ve comparision result */
  struct _AlcKDTNode *childP;   /*!< Child node with +ve comparision result */
  AlcPointP	key;		/*!< The node's key value */
  AlcPointP	boundN; 	/*!< Coordinate of the minimum bounding box
  				     values */
  AlcPointP	boundP; 	/*!< Coordinate of the maximum bounding box
  				     values */
} AlcKDTNode;

/*!
* \struct	_AlcKDTTree
* \ingroup	AlcKDTree
* \brief	A binary space partition tree (kD-tree).
*               Typedef: ::AlcKDTTree
*/
typedef struct _AlcKDTTree
{
   AlcPointType	type;		/*!< The type of tree, ie the type of node
   				     key */
  int		dim;		/*!< Dimension of the tree. */
  size_t	keySz;		/*!< sizeof(key) * dimension */
  double	tol; 		/*!< Comparision tollerance for double key
  				     values */
  size_t	nNodes;  	/*!< Number of nodes in the tree */
  struct _AlcKDTNode *root;	/*!< The root node of the tree. */
  struct _AlcKDTNode *nodeStack; /*!< Stack of nodes available for use */
  size_t	nodeBlockSz; 	/*!< Number of nodes allocated in a block */
  AlcBlockStack *freeStack;	/*!< Stack of allocated node blocks */
} AlcKDTTree;

/*!
* \struct       _AlcHeapEntryCore
* \ingroup      AlcHeap
* \brief        Core heap entry data structure. All other heap entry
*               data structures must have the fields of this data
*               structure first.
*               Typedef: ::AlcHeapEntryCore
*/
typedef struct _AlcHeapEntryCore
{
  double                priority;       /*!< Entry priority highest/lowest
  					     priority at the head of the
					     queue depending on
					     AlcHeap::topPriLo. */
} AlcHeapEntryCore;

/*!
* \struct       _AlcHeap
* \ingroup      AlcHeap
* \brief        A general purpose heap data structure.
*               Typedef: ::AlcHeap
*/
typedef struct _AlcHeap
{
  int                   nEnt;           /*!< Number of entries in use:
                                             incremented when entry inserted,
                                             no change is popped,
                                             decremented when entry freed. */
  int                   maxEnt;         /*!< Number of entries allocated. */
  int                   entInc;         /*!< When allocating space for more
                                             entries, allocate space for
                                             this many at a time. */
  int                   entSz;          /*!< Size of each heap entry. */
  int			topPriLo;	/*!< If non-zero the top priority
  					     is low rather than high. This
					     should be set (default is zero)
					     before any entries are added
					     to the heap. */
  void			*data;		/*!< Application data. */
  void                  *entries;       /*!< Allocated heap entries. */
} AlcHeap;

/*!
* \struct	_AlcUFTree
* \ingroup	AlcUFTree
* \brief 	A general purpose union tree based on Robert Sedgewick's
* 		Weighted Quick Union Find.
*/
typedef struct _AlcUFTree
{
  int           *pr;                    /*!< The parent node of each node in
  					     the tree. */
  int           *sz;                    /*!< The number of nodes in the subtree
  					     of the nodes (including the node
					     itself). */
  int           nCmp;                   /*!< Number of components in the
  					     tree. */
  int           nNod;                   /*!< Number of nodes. */
  int           maxNod;                 /*!< Maximum number of nodes space
                                             allocated for. */
} AlcUFTree;

#ifndef WLZ_EXT_BIND
#ifdef __cplusplus
}					       /* Close scope of 'extern "C" */
#endif
#endif /* WLZ_EXT_BIND */

#endif /* ALCTYPE_H */
