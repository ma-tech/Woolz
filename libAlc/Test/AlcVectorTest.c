#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        AlcVectorTest.c
* Date:         march 2000
* Author:       Bill Hill
* Copyright:	2000 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Simple test program for AlcVector doubly linked lists
* 		of pointers.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <Alc.h>

/* #define ALC_VECTOR_DATA(V,I)  (void*)((char*)(*((V)->blocks+(I)))+(((I)%(V)->blkSz)*(V)->elmSz)) */

static void *ALC_VECTOR_DATA(AlcVector *V, unsigned int I)
{
  char *data;
  data = (char *)(*(V->blocks + (I / V->blkSz))) + ((I % V->blkSz) * V->elmSz);
  return((void *)data);
}

static void	AlcVectorTest(void),
		PrintVector(AlcVector *vec, AlcErrno errNum);

int		main(int argc, char **argv)
{
  AlcVectorTest();
  return(0);
}

static void     AlcVectorTest(void)
{
  AlcVector	*vec = NULL;
  int		idx;
  AlcErrno	errNum = ALC_ER_NONE;

  (void )printf("Simple test program for AlcVector functions\n");
  (void )printf("=======================================\n\n");
  (void )printf("\n* Create a new vector.\n");
  vec = AlcVectorNew(3, sizeof(int), 4, &errNum);
  PrintVector(vec, errNum);
  (void )printf("\n* Set values [0-2] using macro.\n");
  idx = 0;
  while(idx < 3)
  {
    *(int *)(ALC_VECTOR_DATA(vec,idx)) = idx++;
  }
  PrintVector(vec, errNum);
  (void )printf("\n* Extend to 15 values in vector.\n");
  errNum = AlcVectorExtend(vec, 15 );
  PrintVector(vec, errNum);
  (void )printf("\n* Add 100 * idx  to all values [0-14] using function.\n");
  idx = 0;
  while((idx < 15) && (errNum == ALC_ER_NONE))
  {
    *(int *)(AlcVectorItemGet(vec, idx)) += 100 * idx++;
  }
  PrintVector(vec, errNum);
  (void )printf("\n* Free vector.\n");
  errNum = AlcVectorFree(vec);
  PrintVector(NULL, errNum);
}

static void	PrintVector(AlcVector *vec, AlcErrno errNum)
{
  int		pIdx,
  		bIdx;

  (void )printf("errNum = %d\n", errNum);
  if(vec == NULL)
  {
    (void )printf("vec = NULL\n");
  }
  else
  {
    (void )printf("vec = 0x%08lx\n", (unsigned long )vec);
    (void )printf("vec->elmSz = %d\n", vec->elmSz);
    (void )printf("vec->blkCnt = %d\n", vec->blkCnt);
    (void )printf("vec->blkUse = %d\n", vec->blkUse);
    (void )printf("vec->blkSz = %d\n", vec->blkSz);
    (void )printf("vec->freeStack = 0x%lx\n", vec->freeStack);
    (void )printf("vec->blocks = 0x%lx\n", vec->blocks);
    for(pIdx = 0; pIdx < vec->blkCnt; ++pIdx)
    {
      (void )printf(" % 3d 0x%08lx | ",
                    pIdx, (unsigned long )((char *)*(vec->blocks + pIdx)));
      if((char *)*(vec->blocks + pIdx))
      {
	for(bIdx = 0; bIdx < vec->blkSz; ++bIdx)
	{
	  (void )printf("  % 3d",
			*((int *)*(vec->blocks + pIdx) + bIdx));
	}
      }
      printf("\n");
    }
  }
}
