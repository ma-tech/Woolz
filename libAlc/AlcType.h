#ifndef ALCTYPE_H
#define ALCTYPE_H
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        AlcType.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Header file which contains type definitions for the
*		MRC HGU memory allocation library.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 01-11-00 bill Added AlcKDTree, AlcKDTNode, AlcPointP and AlcPointType.
* 02-03-00 bill Added AlcVector.
* 01-11-99 bill Added AlcDLPList.
* 01-12-99 bill	Add AlcBlockStack and many comments.
************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

typedef enum	/* Error codes returned by MRC HGU type allocation functions */
{
  ALC_ER_NONE = 0,
  ALC_ER_ALLOC,
  ALC_ER_NULLPTR,
  ALC_ER_NUMELEM,
  ALC_ER_PARAM
} AlcErrno;

/************************************************************************
* General purpose data structure for maintaining blocks of some data
* type. Useful for efficient memory allocation. It's not a stack
* but a doubly linked list of blocks of data which can be used as a
* stack, heap, list, ....
************************************************************************/
typedef struct _AlcBlockStack
{
  int           elmCnt;                  /* Number of elements used in block */
  int           maxElm;   /* Number of elements space allocated for in block */
  void          *elements;                              /* Block of elements */
  struct _AlcBlockStack *prev;                 /* Previous block up in stack */
  struct _AlcBlockStack *next;   /* Next block down, it's a DLL not a stack! */
} AlcBlockStack;

/************************************************************************
* General purpose doubly linked list.
************************************************************************/
typedef enum					 /* List traversal direction */
{
  ALC_DIRECTION_FWD,				        /* Towards list tail */
  ALC_DIRECTION_REV					/* Towards list head */
} AlcDirection;

typedef struct _AlcDLPItem			  /* Doubly linked list item */
{
  void          (*freeFn)(void *);      /* Fn to free list item, may be NULL */
  void          *entry;				    /* The list item's entry */
  struct _AlcDLPItem *next;		 /* Next item, towards tail, in list */
  struct _AlcDLPItem *prev; 	     /* previous item, towards head, in list */
} AlcDLPItem;

typedef struct _AlcDLPList		   /* Doubly linked list of pointers */
{
  AlcDLPItem *head;				     /* The head of the list */
} AlcDLPList;

/************************************************************************
* General hash table.
************************************************************************/
typedef struct _AlcHashItem				  /* Hash table item */
{
  void		(*freeFn)(void *);      /* Fn to free hash item, may be NULL */
  void		*entry;				    /* The hash item's entry */
  void		*key;		       /* Key which identifies the hash item */
  struct _AlcHashItem *next;		     /* Next item, NULL if last item */
  struct _AlcHashItem *prev; 	        /* Previous item, NULL if first item */
} AlcHashItem;

typedef struct _AlcHashTable				       /* Hash table */
{
  int		(*keyCmp)(void *, void *); 	        /* Key comparison fn */
  unsigned	(*hashFn)(void *);		           /* The hashing fn */
  int		tableSz;		    /* Number of list slots in table */
  AlcHashItem	**table;				   /* Table of lists */
} AlcHashTable;

/************************************************************************
* General purpose vector, an extensible 1D array.
************************************************************************/

typedef struct _AlcVector
{
  unsigned int	elmSz;		/* Size of elements of the vector */
  unsigned int	blkCnt;		/* Number of block pointers */
  unsigned int	blkUse;		/* Number of blocks used */
  unsigned int	blkSz;		/* Number of elements in a block, must NOT be
  				 * changed once vector has been created! */
  void		*freeStack;	/* Free stack */
  void		**blocks;	/* Data blocks */
} AlcVector;

/*
#define ALC_VECTOR_DATA(V,I) (void*)((char*)(*((V)->blocks+(I)))+(((I)%(V)->blkSz)*(V)->elmSz))
*/

/************************************************************************
* General purpose coordinate point.
************************************************************************/
typedef enum
{
  ALC_POINTTYPE_INT,					    /* Integer point */
  ALC_POINTTYPE_DBL		     		     /* Floating point point */
} AlcPointType;

typedef union _AlcPointP
{
  void			*kV;			     /* Pointer to any point */
  int			*kI;			 /* Pointer to integer point */
  double		*kD;		  /* Pointer to floating point point */
} AlcPointP;

/************************************************************************
* General purpose binary space partition tree (kD-tree).
************************************************************************/
typedef struct _AlcKDTNode
{
  int			idx;
  int			split;			  /* The splitting dimension */
  struct _AlcKDTNode	*parent;   /* Parent node, NULL if node is tree root */
  struct _AlcKDTNode	*childN;  /* Child node with -ve comparision result  */
  struct _AlcKDTNode	*childP;  /* Child node with +ve comparision result  */
  AlcPointP		key;				 /* Node's key value */
  AlcPointP		boundN; 	/* Nodes bounding box minimum values */
  AlcPointP		boundP; 	/* Nodes bounding box maximum values */
} AlcKDTNode;

typedef struct _AlcKDTTree
{
   AlcPointType type;			/* Type of tree, ie type of node key */
  int			dim;			   /* Dimension of the tree. */
  int			keySz;			  /* sizeof(key) * dimension */
  double		tol; /* Comparision tollerance for double key values */
  int			nNodes; 	      /* Number of nodes in the tree */
  struct _AlcKDTNode	*root;			   /* Root node of the tree. */
  struct _AlcKDTNode	*nodeStack;		  /* Nodes available for use */
  int			nodeBlockSz; /* Number of nodes allocated in a block */
  AlcBlockStack 	*freeStack;	   /* Stack of allocated node blocks */
} AlcKDTTree;

#ifdef __cplusplus
}					       /* Close scope of 'extern "C" */
#endif

#endif /* ALCTYPE_H */
