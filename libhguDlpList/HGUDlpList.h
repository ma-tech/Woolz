#ifndef HGUDLPLIST_H
#define HGUDLPLIST_H
#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _HGUDlpList_h[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         HGUDlpList.h
* \author       Bill Hill
* \date         March 1999
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
* \brief	Data structures and functions for doubly linked lists
* 		linked lists of pointers.
* \ingroup	hguDlpList
* \todo         -
* \bug          None known.
*/

#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */

/*!
* \enum		_HGUDlpListDirection
* \ingroup	hguDlpList
* \brief	List traversal direction;
* 		Typedef: ::HGUDlpListDirection
*/
typedef enum _HGUDlpListDirection
{
  HGU_DLPLIST_DIR_TOHEAD,
  HGU_DLPLIST_DIR_TOTAIL
} HGUDlpListDirection;

/*!
 * \enum	_HGUDlpListState
 * \ingroup	hguDlpList
 * \brief	State of list locking mechanism.
 * 		Typedef: ::HGUDlpListState
 */
typedef enum _HGUDlpListState
{
  HGU_DLPLIST_STATE_EMPTY 	= (0),
  HGU_DLPLIST_STATE_ERROR 	= (0),
  HGU_DLPLIST_STATE_UNLOCK 	= (1),
  HGU_DLPLIST_STATE_LOCK	= (1<<1),
  HGU_DLPLIST_STATE_CREATE	= (1<<2),
  HGU_DLPLIST_STATE_DESTROY	= (1<<3)
} HGUDlpListState;

#ifdef HGUDLPLIST_C

typedef struct _HGUDlpListItem			  /* Doubly linked list item */
{
  void          (*freeFn)(void *);
  void          *entry;
  struct _HGUDlpListItem *next;
  struct _HGUDlpListItem *prev;
} HGUDlpListItem;
 
typedef struct _HGUDlpList		   /* Doubly linked list of pointers */
{
  HGUDlpListState (*lockFn)(void *, HGUDlpListState);
  void          *lockData;
  int           itemCount;
  HGUDlpListItem *head;
  HGUDlpListItem *tail;
 
} HGUDlpList;

#else /* ! HGUDLPLIST_C */

/*!
 * \typedef	HGUDlpListItem
 * \brief 	Opaque handle for doubly linked list item.
 */
typedef void HGUDlpListItem;

/*!
 * \typedef	HGUDlpList
 * \brief	Opaque handle for doubly linked list of pointers.
 */
typedef void HGUDlpList;

extern int			HGUDlpListSort(
    				  HGUDlpList *list,
			          int (*entryCompFn)(void *, void *));
extern int			HGUDlpListItemIsHead(
    				  HGUDlpList *list,
				  HGUDlpListItem *item);
extern int			HGUDlpListItemIsTail(
    				  HGUDlpList *list,
				  HGUDlpListItem *item);
extern int			HGUDlpListCount(
    				  HGUDlpList *list);
extern int			HGUDlpListOffset(
    				  HGUDlpList *list,
				  HGUDlpListItem *item,
		 		  HGUDlpListDirection dir);
extern void			HGUDlpListDestroy(
    				  HGUDlpList *list);
extern void 			*HGUDlpListEntryGet(
    				  HGUDlpList *list,
				  HGUDlpListItem *item);
extern void			*HGUDlpListEntrySet(
    				  HGUDlpList *list,
				  HGUDlpListItem *item,
				  void *entry);
extern HGUDlpList 		*HGUDlpListCreate(
    				  HGUDlpListState (*lockFn)(void *,
						            HGUDlpListState)),
				  *HGUDlpListDup(HGUDlpList *list);
extern HGUDlpListItem 		*HGUDlpListInsert(HGUDlpList *list,
					HGUDlpListItem *before, void *entry,
					void (*freeFn)(void *));
extern HGUDlpListItem		*HGUDlpListAppend(
    				  HGUDlpList *list,
				  HGUDlpListItem *after,
				  void *entry,
				  void (*freeFn)(void *));
extern HGUDlpListItem		*HGUDlpListExchange(
    				  HGUDlpList *list,
				  HGUDlpListItem *item0,
				  HGUDlpListItem *item1);
extern HGUDlpListItem		*HGUDlpListDeleteAll(
    				  HGUDlpList *list);
extern HGUDlpListItem		*HGUDlpListDelete(
    				  HGUDlpList *list,
				  HGUDlpListItem *item);
extern HGUDlpListItem		*HGUDlpListRemoveAll(
    				  HGUDlpList *list);
extern HGUDlpListItem		*HGUDlpListRemove(
    				  HGUDlpList *list,
				  HGUDlpListItem *item);
extern HGUDlpListItem		*HGUDlpListIterate(
    				  HGUDlpList *list,
				  HGUDlpListItem *item,
				  HGUDlpListDirection dir,
				  int (*iterFn)(HGUDlpList *, HGUDlpListItem *,
						void *),
				  void *iterData);
extern HGUDlpListItem		*HGUDlpListHead(
    				  HGUDlpList *list);
extern HGUDlpListItem		*HGUDlpListTail(
    				  HGUDlpList *list);
extern HGUDlpListItem		*HGUDlpListNext(
    				  HGUDlpList *list, 
				  HGUDlpListItem *item);
extern HGUDlpListItem		*HGUDlpListPrev(
    				  HGUDlpList *list,
				  HGUDlpListItem *item);
extern HGUDlpListItem		*HGUDlpListNth(
    				  HGUDlpList *list,
				  HGUDlpListItem *item,
			       	  HGUDlpListDirection dir,
				  int num);
#endif /* HGUDLPLIST_C */

#ifdef  __cplusplus
}
#endif /* __cplusplus */

#endif /* HGUDLPLIST_H */
