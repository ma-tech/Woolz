#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        AlcDLPListTest.c
* Date:         November 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Simple test program for AlcDLPList doubly linked lists
* 		of pointers.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <Alc.h>

static void	AlcDLPListTest(void),
		AlcPrintAndCheckErr(int errNum);
static int	CompareEntries(void *, void *),
		FindEntryEq(AlcDLPList *, AlcDLPItem *, void *),
		PrintItems(AlcDLPList *, AlcDLPItem *, void *);

int		main(int argc, char **argv)
{
  AlcDLPListTest();
  return(0);
}

static void     AlcDLPListTest(void)
{
  AlcDLPList	*list = NULL;
  AlcDLPItem 	*item0 = NULL,
		*item1;
  int		count;
  AlcErrno	errNum = ALC_ER_NONE;

  (void )printf("Simple test program for AlcDLPList functions\n");
  (void )printf("=======================================\n\n");
  (void )printf("\n* Create a new list.\n");
  list = AlcDLPListNew(&errNum);
  AlcPrintAndCheckErr(errNum);
  (void )printf("\n* Insert items into list.\n");
  count = 0;
  do
  {
    errNum = AlcDLPListEntryInsert(list, NULL, (void *)count, NULL);
    ++count;
    AlcPrintAndCheckErr(errNum);
  } while(count < 10);
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
  item0 = AlcDLPListIterate(list, list->head, ALC_DIRECTION_FWD,
  			    PrintItems, (void *)33, &errNum);
  AlcPrintAndCheckErr(errNum);
  (void )printf("\n* Sort items.\n");
  errNum = AlcDLPListSort(list, CompareEntries);
  AlcPrintAndCheckErr(errNum);
  (void )printf("\n* Find entry == 5.\n");
  item0 = AlcDLPListIterate(list, list->head, ALC_DIRECTION_FWD,
  			    FindEntryEq, (void *)5, &errNum);
  AlcPrintAndCheckErr(errNum);
  (void )printf("\n* Found entry == %d\n", (int )(item0->entry));
  (void )printf("\n* Delete the entry found.\n");
  (void )AlcDLPItemUnlink(list, item0, 1, &errNum);
  AlcPrintAndCheckErr(errNum);
  (void )printf("\n* Print list using iterator function (again).\n");
  item0 = AlcDLPListIterate(list, list->head, ALC_DIRECTION_FWD,
  			   PrintItems, (void *)33, &errNum);
  AlcPrintAndCheckErr(errNum);    		
  (void )printf("\n* Find entry == 3.\n");
  item0 = AlcDLPListIterate(list, list->head, ALC_DIRECTION_FWD,
  			    FindEntryEq, (void *)3, &errNum);
  AlcPrintAndCheckErr(errNum);
  (void )printf("\n* Unlink this item.\n");
  (void )AlcDLPItemUnlink(list, item0, 0, &errNum);
  AlcPrintAndCheckErr(errNum);
  (void )printf("\n* Find entry == 7.\n");
  item1 = AlcDLPListIterate(list, list->head, ALC_DIRECTION_FWD,
  			    FindEntryEq, (void *)7, &errNum);
  AlcPrintAndCheckErr(errNum);
  (void )printf("\n* Append the item with entry == 3 "
  	        "after item with entry == 7.\n");
  errNum = AlcDLPItemAppend(list, item1, item0);
  AlcPrintAndCheckErr(errNum);
  (void )printf("\n* Print list using iterator function (again).\n");
  item0 = AlcDLPListIterate(list, list->head, ALC_DIRECTION_FWD,
  			   PrintItems, (void *)33, &errNum);
  AlcPrintAndCheckErr(errNum);
  (void )printf("\n* Free the list.\n");
  AlcDLPListFree(list);
}

static int	PrintItems(AlcDLPList *list, AlcDLPItem *item, void *data)
{
  int		number;
  
  number = (int )(item->entry);
  (void )printf("%d %d\n", number, (int )data);
  return(item->next != list->head);
}

static int	CompareEntries(void *entry0, void *entry1)
{
  return((int )entry1 - (int )entry0);
}

static void	AlcPrintAndCheckErr(int errNum)
{
  if(errNum != ALC_ER_NONE)
  {
    printf("Error code %d\n", errNum);
    printf("Exit\n");
    exit(0);
  }
  else
  {
    printf("No error\n");
  }
}

static int	FindEntryEq(AlcDLPList *list, AlcDLPItem *item, void *data)
{
  return((item->next != list->head) &&
  	 ((int )(item->entry) != (int )data));
}
