#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        AlcHashTableTest.c
* Date:         November 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Simple test program for AlcHashTable hash table.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <Alc.h>

static void	AlcHashTableTest(void),
		AlcPrintAndCheckErr(int errNum);
static unsigned	HashFn(void *);
static int	KeyCmp(void *, void *),
		IsModZero(AlcHashTable *, AlcHashItem *, void *),
		PrintItems(AlcHashTable *, AlcHashItem *, void *);

int		main(int argc, char **argv)
{
  AlcHashTableTest();
  return(0);
}

static void     AlcHashTableTest(void)
{
  AlcHashTable	*hTbl = NULL;
  AlcHashItem 	*item0 = NULL,
		*item1;
  int		count;
  unsigned	key;
  AlcErrno	errNum = ALC_ER_NONE;

  (void )printf("Simple test program for AlcHashTable functions\n");
  (void )printf("=======================================\n\n");
  (void )printf("\n* Create a new hash table.\n");
  hTbl = AlcHashTableNew(5, KeyCmp, HashFn, &errNum);
  AlcPrintAndCheckErr(errNum);
  (void )printf("\n* Insert items into hash table.\n");
  count = 0;
  do
  {
    errNum = AlcHashTableEntryInsert(hTbl, (void *)count, (void *)count, NULL);
    ++count;
    AlcPrintAndCheckErr(errNum);
  } while(count < 20);
  (void )printf("\n* Print hash table using iterator function.\n");
  (void )printf("\nExpected output:\n"
		"15 10 5 0\n"
		"16 11 6 1\n"
		"17 12 7 2\n"
		"18 13 8 3\n"
		"19 14 9 4\n");
  (void )printf("\nActual output:\n");
  item0 = AlcHashTableIterate(hTbl, ALC_DIRECTION_FWD,
  			      PrintItems, NULL, &errNum);
  AlcPrintAndCheckErr(errNum);
  (void )printf("\n* Remove all entries divisible by 3 without remainder.\n");
  errNum = AlcHashTableUnlinkAll(hTbl, IsModZero, (void *)3, 1);
  AlcPrintAndCheckErr(errNum);
  (void )printf("\n* Print list using iterator function (again).\n");
  item0 = AlcHashTableIterate(hTbl, ALC_DIRECTION_FWD,
  			      PrintItems, NULL, &errNum);
  AlcPrintAndCheckErr(errNum);
  (void )printf("\n* Find random entries and delete them.\n");
  srandom(0);
  do
  {
    key = random() % 20;
    item0 = AlcHashItemGet(hTbl, (void *)key, &errNum);
    if(item0)
    {
      AlcPrintAndCheckErr(errNum);
      errNum = AlcHashItemUnlink(hTbl, item0, 1);
      AlcPrintAndCheckErr(errNum);
    }
    (void )printf("\n* Key == %d, found == %s %d\n",
    		  key,
		  (item0)? "": "NULL",
		  (item0)? (int )(item0->entry): 0);
    item0 = AlcHashTableIterate(hTbl, ALC_DIRECTION_FWD,
				PrintItems, NULL, &errNum);
  } while(item0 && *(hTbl->table));
  (void )printf("\n* Print list using iterator function (again).\n");
  item0 = AlcHashTableIterate(hTbl, ALC_DIRECTION_FWD,
  			      PrintItems, NULL, &errNum);
  AlcPrintAndCheckErr(errNum);
  (void )printf("\n* Free the hash table.\n");
  AlcHashTableFree(hTbl);
}

static int	IsModZero(AlcHashTable *hTbl, AlcHashItem *item, void *data)
{
  int		number,
  		mod;
  number = (int )(item->entry);
  mod = (int )data;

  return((number % mod) == 0);
}

static int	PrintItems(AlcHashTable *hTbl, AlcHashItem *item, void *data)
{
  int		number;
  
  number = (int )(item->entry);
  (void )printf("%d%s", number, (item->next)? " ": "\n");
  return(1);
}

static int	KeyCmp(void *entry0, void *entry1)
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

static unsigned	HashFn(void *data)
{
  unsigned	hashVal;

  hashVal = (unsigned )data;
  return(hashVal);
}


