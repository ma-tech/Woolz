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

typedef enum					 /* List traversal direction */
{
  ALC_DIRECTION_FWD,
  ALC_DIRECTION_REV
} AlcDirection;

typedef struct _AlcDLPItem			  /* Doubly linked list item */
{
  void          (*freeFn)(void *);
  void          *entry;
  struct _AlcDLPItem *next;
  struct _AlcDLPItem *prev;
} AlcDLPItem;

typedef struct _AlcDLPList		   /* Doubly linked list of pointers */
{
  AlcDLPItem *head;
} AlcDLPList;

typedef struct _AlcHashItem				  /* Hash table item */
{
  void		(*freeFn)(void *);
  void		*entry;
  void		*key;
  struct _AlcHashItem *next;		     /* Next item, NULL if last item */
  struct _AlcHashItem *prev; 	        /* Previous item, NULL if first item */
} AlcHashItem;

typedef struct _AlcHashTable				       /* Hash table */
{
  int		(*keyCmp)(void *, void *);
  unsigned	(*hashFn)(void *);
  int		tableSz;
  AlcHashItem	**table;
} AlcHashTable;

#ifdef __cplusplus
}					       /* Close scope of 'extern "C" */
#endif

#endif /* ALCTYPE_H */
