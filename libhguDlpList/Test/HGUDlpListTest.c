#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        HGUDlpListTest.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Simple test program for libhguDlpList.a, a library for
*		handling doubly linked lists of pointers.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <HGUDlpList.h>

static void	HGUDlpListTest(void);
static int	CompareEntries(void *, void *),
		PrintItems(HGUDlpList *, HGUDlpListItem *, void *),
		FindEqualItem(HGUDlpList *, HGUDlpListItem *, void *);
static HGUDlpListState LockFunc(void *, HGUDlpListState);

int		main(int argc, char **argv)
{
  HGUDlpListTest();
  return(0);
}

static void     HGUDlpListTest(void)
{
  HGUDlpList	*list = NULL;
  HGUDlpListItem *item = NULL,
  		*head,
		*tail;
  int		count;

  (void )printf("Simple test program for libhguDlpList.a\n");
  (void )printf("=======================================\n\n");
  (void )printf("\n* Create list with dummy lock function "
  		"(lock function prints 666).\n");
  list = HGUDlpListCreate(LockFunc);
  assert(list != NULL);
  (void )printf("\n* Insert items into list.\n");
  count = 0;
  do
  {
    item = HGUDlpListInsert(list, item, (void *)count, NULL);
    ++count;
  } while(item && (count < 10));
  assert(item != NULL);
  (void )printf("\n* Find head of list.\n");
  head = HGUDlpListHead(list);
  assert(head != NULL);
  (void )printf("\n* Find tail of list.\n");
  tail = HGUDlpListTail(list);
  assert(tail != NULL);
  (void )printf("\n* Print list using iterator function.\n");
  (void )printf("\nExpected output:\n"
		"9 33\n"
		"8 33\n"
		"7 33\n"
		"6 33\n"
		"5 33\n"
		"4 33\n"
		"3 33\n"
		"2 33\n"
		"1 33\n"
		"0 33\n");
  (void )printf("\nActual output:\n");
  item = HGUDlpListIterate(list, head, HGU_DLPLIST_DIR_TOTAIL,
  			   PrintItems, (void *)33);
  assert((unsigned )item == (unsigned)tail);
  (void )printf("\n* Find entry in the list using iterator function.\n");
  item = HGUDlpListIterate(list, head, HGU_DLPLIST_DIR_TOTAIL,
  			   FindEqualItem, (void *)5);
  (void )printf("\n* Get entry for the item found.\n");
  count = (int )HGUDlpListEntryGet(list, item);
  (void )printf("\nExpected output:\n"
                "count 5\n");
  (void )printf("\nActual output:\n"
  		"count %d\n", count);
  assert(count == 5);
  (void )printf("\n* Delete this entry from list.\n");
  item = HGUDlpListDelete(list, item);
  (void )printf("\n* Print list using iterator function.\n");
  (void )printf("\nExpected output:\n"
		"9 33\n"
		"8 33\n"
		"7 33\n"
		"6 33\n"
		"4 33\n"
		"3 33\n"
		"2 33\n"
		"1 33\n"
		"0 33\n");
  (void )printf("\nActual output:\n");
  item = HGUDlpListIterate(list, head, HGU_DLPLIST_DIR_TOTAIL,
  			    PrintItems, (void *)33);
  assert((unsigned )item == (unsigned )tail);
  (void )printf("\n* Find entry in the list using iterator function.\n");
  item = HGUDlpListIterate(list, head, HGU_DLPLIST_DIR_TOTAIL,
  			   FindEqualItem, (void *)3);
  (void )printf("\n* Get entry for the item found.\n");
  count = (int )HGUDlpListEntryGet(list, item);
  (void )printf("\nExpected output:\n"
  		"count 3\n");
  (void )printf("\nActual output:\n"
  		"count %d\n", count);
  assert(count == 3);
  (void )printf("\n* Exchange entry at head of this with this item.\n");
  item = HGUDlpListExchange(list, head, item);
  (void )printf("\n* Print list using iterator function.\n");
  (void )printf("\nExpected output:\n"
  		"3 33\n"
		"8 33\n"
		"7 33\n" 
		"6 33\n"
		"4 33\n"
		"9 33\n"
		"2 33\n"
		"1 33\n"
		"0 33\n");
  (void )printf("\nActual output:\n");
  item = HGUDlpListIterate(list, head, HGU_DLPLIST_DIR_TOTAIL,
  			   PrintItems, (void *)33);
  (void )printf("\n* Find entry in the list using iterator function.\n");
  item = HGUDlpListIterate(list, head, HGU_DLPLIST_DIR_TOTAIL,
  			   FindEqualItem, (void *)7);
  (void )printf("\n* Get entry for the item found.\n");
  count = (int )HGUDlpListEntryGet(list, item);
  (void )printf("\nExpected output:\n"
  		"count 7\n");
  (void )printf("\nActual output:\n"
  		"count %d\n", count);
  assert(count == 7);
  (void )printf("\n* Set entry for the item found.\n");
  count = (int )HGUDlpListEntrySet(list, item, (void *)77);
  (void )printf("\nExpected output:\n"
  		"count 7\n");
  (void )printf("\nActual output:\n"
		"count %d\n", count);
  assert(count == 7);
  (void )printf("\n* Print list using iterator function.\n");
  (void )printf("\nExpected output:\n"
  		"3 33\n"
		"8 33\n"
		"77 33\n" 
		"6 33\n"
		"4 33\n"
		"9 33\n"
		"2 33\n"
		"1 33\n"
		"0 33\n");
  (void )printf("\nActual output:\n");
  item = HGUDlpListIterate(list, head, HGU_DLPLIST_DIR_TOTAIL,
  			   PrintItems, (void *)33);
  (void )printf("\n* Find entry in the list using iterator function.\n");
  item = HGUDlpListIterate(list, head, HGU_DLPLIST_DIR_TOTAIL,
  			   FindEqualItem, (void *)6);
  count = (int )HGUDlpListEntryGet(list, item);
  (void )printf("\nExpected output:\n"
  		"count 6\n");
  (void )printf("\nActual output:\n"
  		"count %d\n", count);
  assert(count == 6);
  (void )printf("\n* Append new item after this item.\n");
  item = HGUDlpListAppend(list, item, (void *)55, NULL);
  count = (int )HGUDlpListEntryGet(list, item);
  (void )printf("\nExpected output:\n"
                "count 55\n");
  (void )printf("\nActual output:\n"
                "count %d\n", count); 
  (void )printf("\n* Print list using iterator function.\n");
  (void )printf("\nExpected output:\n"
  		"3 33\n"
		"8 33\n"
		"77 33\n" 
		"6 33\n"
		"55 33\n"
		"4 33\n"
		"9 33\n"
		"2 33\n"
		"1 33\n"
		"0 33\n");
  (void )printf("\nActual output:\n");
  item = HGUDlpListIterate(list, head, HGU_DLPLIST_DIR_TOTAIL,
  			   PrintItems, (void *)33);
  (void )printf("\n* Sort the list entries.\n");
  count = (int )HGUDlpListSort(list, CompareEntries);
  (void )printf("\nExpected output:\n"
                "count 10\n");
  (void )printf("\nActual output:\n"
                "count %d\n", count);
  assert(count == 10);
  (void )printf("\n* Print list using iterator function.\n");
  (void )printf("\nExpected output:\n"
                "77 33\n"
                "55 33\n"
                "9 33\n"
                "8 33\n"
                "6 33\n"
                "4 33\n"
                "3 33\n"
                "2 33\n"
                "1 33\n"
                "0 33\n");
  (void )printf("\nActual output:\n");
  item = HGUDlpListIterate(list, head, HGU_DLPLIST_DIR_TOTAIL,
                           PrintItems, (void *)33);
  (void )printf("\n* Find 6th item from the head.\n");
  item = HGUDlpListNth(list, NULL, HGU_DLPLIST_DIR_TOTAIL, 6);
  count = (int )HGUDlpListEntryGet(list, item);
  (void )printf("\nExpected output:\n"
                "3\n");
  assert(count == 3);
  (void )printf("\nActual output:\n"
                "%d\n", count);
  (void )printf("\n* Find offset of this item from the head.\n");
  count = (int )HGUDlpListOffset(list, item, HGU_DLPLIST_DIR_TOHEAD);
  (void )printf("\nExpected output:\n"
  		"6\n");
  assert(count == 6);
  (void )printf("\nActual output:\n"
  		"%d\n", count);
  (void )printf("\n* Find offset of this item from the tail.\n");
  count = (int )HGUDlpListOffset(list, item, HGU_DLPLIST_DIR_TOTAIL);
  (void )printf("\nExpected output:\n"
  		"3\n");
  assert(count == 3);
  (void )printf("\nActual output:\n"
  		"%d\n", count);
  (void )printf("\n* Destroy the list and all its items.\n");
  HGUDlpListDestroy(list);
}

static int	FindEqualItem(HGUDlpList *list, HGUDlpListItem *item,
			      void *magic)
{
  int		notEqual,
  		have,
  		target;

  have = (int )HGUDlpListEntryGet(list, item);
  target = (int )magic;
  notEqual = have != target;
  return(notEqual);
}

static int	PrintItems(HGUDlpList *list, HGUDlpListItem *item, void *magic)
{
  int		number,
  		thirtyThree;
  
  number = (int )HGUDlpListEntryGet(list, item);
  thirtyThree = (int )magic;
  (void )printf("%d %d\n", number, thirtyThree);
  return(1);
}

static HGUDlpListState LockFunc(void *lockData, HGUDlpListState state)
{
  HGUDlpListState newState = HGU_DLPLIST_STATE_EMPTY;
  int		*magic = NULL;

  if(lockData && (state & (HGU_DLPLIST_STATE_CREATE|
  			   HGU_DLPLIST_STATE_LOCK|
  			   HGU_DLPLIST_STATE_UNLOCK|
			   HGU_DLPLIST_STATE_DESTROY)))
  {
    if(state & HGU_DLPLIST_STATE_CREATE)
    {
      (void )printf("Create lock\n");
      magic = malloc(sizeof(int));
      *((int **)lockData) = magic;
      if(magic)
	*magic = 666;
    }
    else
      magic = (int *)lockData;
    if(state &  HGU_DLPLIST_STATE_LOCK)
    {
      (void )printf("Lock %d\n", *magic);
      newState |= HGU_DLPLIST_STATE_LOCK;
    }
    if(state & HGU_DLPLIST_STATE_UNLOCK)
    {
      (void )printf("Unlock %d\n", *magic);
      newState |= HGU_DLPLIST_STATE_UNLOCK;
    }
    if(state & HGU_DLPLIST_STATE_DESTROY)
    {
      (void )printf("Destroy lock %d\n", *magic);
      free(lockData);
      newState |= HGU_DLPLIST_STATE_DESTROY;
    }
  }
  return(newState);
}

static int	CompareEntries(void *entry0, void *entry1)
{
  return((int )entry0 - (int )entry1);
}
