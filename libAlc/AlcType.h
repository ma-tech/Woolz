#ifndef ALCTYPE_H
#define ALCTYPE_H
#pragma ident "MRC HGU $Id$"
/*!
* \file         AlcType.h
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
* \brief        Type definitions for the MRC HGU memory allocation
*		and fundamental type library.
* \todo		-
* \bug          None known.
*/

#ifdef __cplusplus
extern "C" {
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
  void          (*freeFn)(void *); /*!< Function to free list item, may be
  					NULL */
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
  				     list, NULL if this is thefirst item */
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
  unsigned	(*hashFn)(void *);         /*!< The hashing function */
  int		tableSz;		   /*!< Number of list slots in the
  						hash table */
  AlcHashItem	**table;		   /*!< Table of lists */
} AlcHashTable;



/*!
* \struct	_AlcVector
* \ingroup      AlcVector
* \brief	An extensible 1D array.
*               Typedef: ::AlcVector
*/
typedef struct _AlcVector
{
  unsigned int	elmSz;		/*!< Size of elements of the vector */
  unsigned int	blkCnt;		/*!< Number of block pointers */
  unsigned int	blkUse;		/*!< Number of blocks used */
  unsigned int	blkSz;		/*!< Number of elements in a block, must NOT be
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
  int		idx;		/*!< Index or identifier for the node */
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
  int		keySz;		/*!< sizeof(key) * dimension */
  double	tol; 		/*!< Comparision tollerance for double key
  				     values */
  int		nNodes;  	/*!< Number of nodes in the tree */
  struct _AlcKDTNode *root;	/*!< The root node of the tree. */
  struct _AlcKDTNode *nodeStack; /*!< Stack of nodes available for use */
  int		nodeBlockSz; 	/*!< Number of nodes allocated in a block */
  AlcBlockStack *freeStack;	/*!< Stack of allocated node blocks */
} AlcKDTTree;

#ifdef __cplusplus
}					       /* Close scope of 'extern "C" */
#endif

#endif /* ALCTYPE_H */
