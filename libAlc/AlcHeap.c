#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlcHeap_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlc/AlcHeap.c
* \author       Zsolt Husz, Bill Hill
* \date         February 2009
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
* \brief	A basic heap data structure.
* \ingroup	AlcHeap
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Alc.h>

/*!
* \return       New heap data structure, NULL on error.
* \ingroup      AlcHeap
* \brief        Constructs a new heap data structure.
*
* 		The following example code illustrates the use of the
* 		heap data structures and functions, although (for clarity
* 		and succinctness) without any error checking.
*   \verbatim
    typedef struct _AlcHeapEntryTest
    {
      double    priority;
      int       index; 
    } AlcHeapEntryTest;

    int             main(int argc, char *argv[])
    {
      int           idN,
                    size = 10;
      AlcHeap       *heap = NULL;
      double        *randoms = NULL;
      AlcHeapEntryTest ent;
      AlcHeapEntryTest *entP;
      AlcErrno      alcErr = ALC_ER_NONE;

      randoms = (double *)AlcMalloc(sizeof(double) * size);
      srand(0);
      for(idN = 0; idN < size; ++idN)
      {
        randoms[idN] = (double)rand() / (double )(RAND_MAX);
      }
      heap = AlcHeapNew(sizeof(AlcHeapEntryTest), (size + 1) / 2, NULL);
      for(idN = 0; idN < size; ++idN)
      {
        ent.priority = randoms[idN];
        ent.index = idN;
        AlcHeapInsertEnt(heap, &ent);
      }
      for(idN = 0; idN < size; ++idN)
      {
        entP = AlcHeapTop(heap);
        (void )printf("% 8d % 8d %1.12lf\n", idN, entP->index, entP->priority);
        AlcHeapEntFree(heap);
      }
      AlcHeapFree(heap);
      AlcFree(randoms);
    }
    \endverbatim
*
*		entSz			Size of the heap entries.
*		entInc			Number of entries to allocate at once
*					or zero for the default (1024).
* 		data			Application specific data which
* 					is ignored by heap functions,
* 					may be NULL.
*/
AlcHeap 	*AlcHeapNew(int entSz, int entInc, void *data)
{
  AlcHeap 	*heap = NULL;
  const int	minEntInc = 1024;

  if(entInc < minEntInc)
  {
    entInc = minEntInc;
  }
  if((heap = (AlcHeap *)AlcCalloc(1, sizeof(AlcHeap))) != NULL)
  {
    heap->entSz = entSz;
    heap->entInc = entInc;
    heap->data = data;
  }
#ifdef ALC_HEAP_DEBUG
  (void )fprintf(stderr, "AlcHeapNew()\n");
#endif
  return(heap);
}

/*!
* \ingroup      AlcHeap
* \brief        Frees the given heap data structure along with it's entries.
* \param        heap                   The heap to free.
*/
void     	AlcHeapFree(AlcHeap *heap)
{
  if(heap != NULL)
  {
    AlcFree(heap->entries);
    AlcFree(heap);
  }
#ifdef ALC_HEAP_DEBUG
  (void )fprintf(stderr, "AlcHeapFree()\n");
#endif
}

/*!
* \ingroup      AlcHeap
* \brief        Frees the entry at the top of the heap.
* \param        heap                   The heap.
*/
void     	AlcHeapEntFree(AlcHeap *heap)
{
  int           id0,
                id1;
  void		*ent0,
		*ent1,
		*ent2;
  double        p0,
  		p1,
		p2;

  if(heap != NULL)
  {
    if(--(heap->nEnt) >= 0)
    {
      p0 = ((AlcHeapEntryCore *)((char *)(heap->entries) +
                                 (heap->nEnt * heap->entSz)))->priority;
      id0 = 0;
      id1 = ((id0 + 1) * 2) - 1;
      while(id1 <= heap->nEnt)
      {
	ent1 = (void *)((char *)(heap->entries) + (id1 * heap->entSz));
	ent2 = (void *)((char *)ent1 + heap->entSz);
	p1 = ((AlcHeapEntryCore *)ent1)->priority;
	p2 = ((AlcHeapEntryCore *)ent2)->priority;
	if((id1 + 1 <= heap->nEnt) && (p1 > p2))
	{
	  ++id1;
	  ent1 = ent2;
	  p1 = p2;
	}
	if(p1 < p0)
	{
	  ent0 = (void *)((char *)(heap->entries) + (id0 * heap->entSz));
	  (void )memcpy(ent0, ent1, heap->entSz);
	  id0 = id1;
	  id1 = ((id0 + 1) * 2) - 1;
	}
	else
	{
	  break;
	}
      }
      ent0 = (void *)((char *)(heap->entries) + (id0 * heap->entSz));
      ent1 = (void *)((char *)(heap->entries) + (heap->nEnt * heap->entSz));
      (void )memmove(ent0, ent1, heap->entSz);
    }
  }
#ifdef ALC_HEAP_DEBUG
  (void )fprintf(stderr, "AlcHeapEntFree()\n");
  if(heap != NULL)
  {
    (void )fprintf(stderr, "heap->nEnt = %d\n", heap->nEnt);
    for(id0 = 0; id0 < heap->nEnt; ++id0)
    {
      ent0 = heap->entries + (id0 * heap->entSz);
      p0 = ((AlcHeapEntryCore *)ent0)->priority;
      (void )fprintf(stderr, "% 8d %g\n", id0, p0);
    }
    (void )fprintf(stderr, "\n");
  }
#endif
}

/*!
* \ingroup      AlcHeap
* \brief        Frees all the heap entries.
* \param        heap                   The heap.
* \param	reallyFree		If non zero the entries are really
* 					freed rather than just marked free
* 					and available for reuse.
*/
void     	AlcHeapAllEntFree(AlcHeap *heap, int reallyFree)
{
  if(heap != NULL)
  {
    heap->nEnt = 0;
    if(reallyFree != 0)
    {
      AlcFree(heap->entries);
      heap->entries = NULL;
    }
  }
#ifdef ALC_HEAP_DEBUG
  (void )fprintf(stderr, "AlcHeapAllEntFree()\n");
#endif
}

/*!
* \ingroup      AlcHeap
* \return       Error code.
* \brief        Inserts the given entry into the queue.
* \param        heap                  Given heap data structure.
* \param        ent                    Entry to insert.
*/
AlcErrno	AlcHeapInsertEnt(AlcHeap *heap, void *ent)
{
  int           id0,
                id1;
  double	p,
                p1;
  void		*ent0,
  		*ent1;
  AlcErrno   	alcErr = ALC_ER_NONE;

  if((heap != NULL) && (ent != NULL))
  {
    id0 = (heap->nEnt)++;
    if(heap->nEnt > heap->maxEnt)
    {
      heap->maxEnt += heap->entInc;
      if((heap->entries = AlcRealloc(heap->entries,
			             heap->entSz * heap->maxEnt)) == NULL)
      {
	alcErr = ALC_ER_ALLOC;
      }
    }
    if(alcErr == ALC_ER_NONE)
    {
      p = ((AlcHeapEntryCore *)ent)->priority;
      while(id0 > 0)
      {
	id1 = ((id0 + 1) / 2) - 1;
        ent1 = (void *)((char *)(heap->entries) + (id1 * heap->entSz));
	p1 = ((AlcHeapEntryCore *)ent1)->priority;
	if(p1 > p)
	{
          ent0 = (void *)((char *)(heap->entries) + (id0 * heap->entSz));
          (void )memcpy(ent0, ent1, heap->entSz);
	}
	else
	{
	  break;
	}
	id0 = id1;
      }
      ent0 = (void *)((char *)(heap->entries) + (id0 * heap->entSz));
      (void )memcpy(ent0, ent, heap->entSz);
    }
  }
#ifdef ALC_HEAP_DEBUG
  (void )fprintf(stderr, "AlcHeapInsertEnt()\n");
  if((heap != NULL) && (ent != NULL))
  {
    if(alcErr == ALC_ER_NONE)
    {
      (void )fprintf(stderr, "heap->nEnt = %d\n", heap->nEnt);
      for(id1 = 0; id1 < heap->nEnt; ++id1)
      {
	ent1 = heap->entries + (id1 * heap->entSz);
        p1 = ((AlcHeapEntryCore *)ent1)->priority;
	(void )fprintf(stderr, "% 8d %g\n", id1, p1);
      }
      (void )fprintf(stderr, "\n");
    }
  }
#endif
  return(alcErr);
}

/*!
* \return       The entry at the top of the heap.
* \ingroup      AlcHeap
* \brief        Gets the top heap entry.
* \param        heap                   Given priority queue.
*/
void		 *AlcHeapTop(AlcHeap *heap)
{
  void	 	*ent = NULL;

  if((heap != NULL) && (heap->nEnt > 0))
  {
    ent = heap->entries;
  }
#ifdef ALC_HEAP_DEBUG
  if(heap != NULL)
  {
    (void )fprintf(stderr, "AlcHeapTop()\n");
  }
#endif
  return(ent);
}
#ifdef ALC_HEAP_TEST_MAIN_1
/* Test #1
 * Very basic test - probably only useful to get profiles.
 */
#include <sys/time.h>
#include<Alg.h>

typedef struct _AlcHeapEntryTest
{
  double	priority;   /* Required as first field of any heap entry. */
  int		index;      /* Example of application data. */
} AlcHeapEntryTest;

static void	BuildAndKnockDown(AlcHeap *heap, int size, double *randoms)
{
  int		idN;
  AlcHeapEntryTest ent;
  AlcHeapEntryTest *entP;

  for(idN = 0; idN < size; ++idN)
  {
    ent.priority = randoms[idN];
    ent.index = idN;
    AlcHeapInsertEnt(heap, &ent);
  }
  for(idN = 0; idN < size; ++idN)
  {
    entP = AlcHeapTop(heap);
    AlcHeapEntFree(heap);
  }
}

int             main(int argc, char *argv[])
{
  int		idN,
                repeats = 1000,
  		size = 1000;
  AlcHeap	*heap = NULL;
  double	*randoms = NULL;
  AlcHeapEntryTest ent;
  AlcHeapEntryTest *entP;
  AlcErrno	alcErr = ALC_ER_NONE;

  randoms = (double *)AlcMalloc(sizeof(double) * size);
  AlgRandSeed(0);
  for(idN = 0; idN < size; ++idN)
  {
    randoms[idN] = AlgRandUniform();
  }
  heap = AlcHeapNew(sizeof(AlcHeapEntryTest), (size + 1) / 2, NULL);
  for(idN = 0; idN < repeats; ++idN)
  {
    BuildAndKnockDown(heap, size, randoms);
  }
  AlcHeapFree(heap);
  AlcFree(randoms);
}
#endif

#ifdef ALC_HEAP_TEST_MAIN_2
/* Test #2
 * This builds a heap and then knocks it down again.
 */

#include <sys/time.h>
#include<Alg.h>

typedef struct _AlcHeapEntryTest
{
  double	priority;   /* Required as first field of any heap entry. */
  int		index;      /* Example of application data. */
} AlcHeapEntryTest;

extern char     *optarg;
extern int      optind,
		opterr,
		optopt;

int             main(int argc, char *argv[])
{
  int		idN,
		idR,
		option,
  		repeats = 1,
  		size = 10,
		time = 0,
		usage = 0,
		verbose = 0;
  double	execTime;
  struct timeval tp;
  AlcHeap	*heap = NULL;
  double	*randoms = NULL;
  AlcHeapEntryTest ent;
  AlcHeapEntryTest *entP;
  const char   optList[] = "hr:s:tv";
  AlcErrno	alcErr = ALC_ER_NONE;

  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'r':
	if(sscanf(optarg, "%d", &repeats) != 1)
	{
	  usage = 1;
	}
        break;
      case 's':
	if(sscanf(optarg, "%d", &size) != 1)
	{
	  usage = 1;
	}
        break;
      case 't':
	time = 1;
        break;
      case 'v':
	verbose = 1;
        break;
      case 'h': /* FALLTHROUGH */
      default:
        usage = 1;
	break;
    }
  }
  if(usage != 0)
  {
    (void )fprintf(stderr,
    "Usage: %s%s",
    *argv,
    " [-h] [-r#] [-s#] [-t] [-v]\n"
    "Test for the AlcHeap data structure which builds a heap and then\n"
    "breaks it down again.\n"
    "Options:\n"
    " -h  Prints this usage information.\n"
    " -r  Number of repeats building and breaking down heap.\n"
    " -s  Size of heap (number of entries).\n"
    " -t  Output time to build and break down heap for all repeats.\n");
    (void )exit(1);
  }
  if((randoms = (double *)AlcMalloc(sizeof(double) * size)) == NULL)
  {
    (void )fprintf(stderr, "%s: Failed to allocate test data.\n", *argv);
    exit(1);
  }
  AlgRandSeed(0);
  for(idN = 0; idN < size; ++idN)
  {
    randoms[idN] = AlgRandUniform();
  }
  if((heap = AlcHeapNew(sizeof(AlcHeapEntryTest),
                        (size + 1) / 2, NULL)) == NULL)
  {
    (void )fprintf(stderr, "%s: Failed to create a heap.\n", *argv);
    exit(1);
  }
  if(time != 0)
  {
    (void )gettimeofday(&tp, NULL);
    execTime = (double )tp.tv_sec + (1.0e-6) * tp.tv_usec;
  }
  for(idR = 0; idR < repeats; ++idR)
  {
    for(idN = 0; idN < size; ++idN)
    {
      ent.priority = randoms[idN];
      ent.index = idN;
      if(AlcHeapInsertEnt(heap, &ent) != ALC_ER_NONE)
      {
	break;
      }
    }
    if(alcErr != ALC_ER_NONE)
    {
      (void )fprintf(stderr, "%s: Failed to insert entry (%d) into heap.\n",
		      *argv, idN);
      exit(1);
    }
    for(idN = 0; idN < size; ++idN)
    {
      entP = AlcHeapTop(heap);
      if(entP == NULL)
      {
	break;
      }
      if(verbose)
      {
	(void )printf("% 8d % 8d %1.12lf\n", idN, entP->index, entP->priority);
      }
      AlcHeapEntFree(heap);
    }
    if((idN != 0) && (entP == NULL))
    {
      (void )fprintf(stderr, "%s: Failed to get top entry from heap.\n",
		     *argv, idN);
      exit(1);
    }
  }
  if(time != 0)
  {
    (void )gettimeofday(&tp, NULL);
    execTime = (double )tp.tv_sec + (1.0e-6) * tp.tv_usec - execTime;
    (void )printf("%s: time for %d repeats with heap size %d = %gs.\n",
                  *argv, repeats, size, execTime);
  }
  AlcHeapFree(heap);
  AlcFree(randoms);
}
#endif
