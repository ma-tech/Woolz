#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTstObjectCache_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzTst/WlzTstObjectCache.c
* \author       Bill Hill
* \date         November 2011
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
* \brief	Test program for a least recent removal Woolz object cache.
* \ingroup 	BinWlzTst
*/

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <Wlz.h>

typedef struct _WlzTstObjCacheEntry
{
  size_t	size;
  char		*file;
  WlzObject	*obj;
} WlzTstObjCacheEntry;

static int			WlzTstObjCacheCmpFn(
				  const void *e0,
				  const void *e1);
static unsigned int 		WlzTstObjCacheKeyFn(
				  AlcLRUCache *cache,
				  const void *e);
static void			WlzTstObjCacheUnlinkFn(
				  AlcLRUCache *cache,
				  const void *e);
static WlzErrorNum 		WlzTstObjCacheAdd(
				  AlcLRUCache *cache,
				  WlzTstObjCacheEntry *ent,
				  unsigned int *hit,
				  unsigned int *colSz);

extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;

int		main(int argc, char *argv[])
{
  int		pLen,
  		option,
		debug = 0,
  		ok = 1,
		random = 0,
  		usage = 0;
  unsigned int	nAdd = 0,
		nCol = 0,
  		nHit = 0,
		objCnt = 0,
  		maxObj = 8,
  		repeats = 1,
		maxSz = 0,
		maxColSz = 0,
		maxTotCacheItems = 0;
  size_t	maxTotCacheEntSz = 0;
  char		*path;
  DIR		*dir;
  struct dirent *dP;
  struct stat stbuf;
  AlcLRUCache	*cache = NULL;
  WlzTstObjCacheEntry *objTbl = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  static char   optList[] = "dhmn:r:s:";

  opterr = 0;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'd':
	debug = 1;
	break;
      case 'm':
	random = 1;
	break;
      case 'n':
	usage = (sscanf(optarg, "%u", &maxObj) != 1);
	break;
      case 'r':
        usage = (sscanf(optarg, "%u", &repeats) != 1);
	break;
      case 's':
        usage = (sscanf(optarg, "%u", &maxSz) != 1);
	break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  if((optind + 1) == argc)
  {
    path = argv[optind];
    if((path == NULL) || (*path == '\0'))
    {
      usage = 1;
    }
  }
  else
  {
    usage = 1;
  }
  ok = !usage;
  if(ok)
  {
    if((dir = opendir(path)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open directory %s.\n",
                     *argv, path);
    }
    else
    {
      while((dP = readdir(dir)) != NULL)
      {
	int	len;

        if(((len = strlen(dP->d_name)) > 4) &&
	   (strcmp(dP->d_name + len - 4, ".wlz") == 0))
	{
	  ++objCnt;
	}
      }
      (void )closedir(dir);
    }
    if(objCnt < 1)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Directory %s contains no Woolz files.\n",
                     *argv, path);
    }
  }
  if(ok)
  {
    if((objTbl = (WlzTstObjCacheEntry *)
                 AlcCalloc(objCnt, sizeof(WlzTstObjCacheEntry))) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to allocate Woolz object table.\n",
		     *argv);
    }
  }
  if(ok)
  {
    pLen = strlen(path);
    if((dir = opendir(path)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open directory %s.\n",
                     *argv, path);
    }
    else
    {
      objCnt = 0;
      while(ok && ((dP = readdir(dir)) != NULL))
      {
	int	len;

        if(((len = strlen(dP->d_name)) > 4) &&
	   (strcmp(dP->d_name + len - 4, ".wlz") == 0))
	{
	  if((objTbl[objCnt].file = AlcStrCat3(path, "/",
	                                       dP->d_name)) == NULL)
	  {
	    ok = 0;
	    (void )fprintf(stderr,
	                   "%s: Failed to allocate file path string.\n",
			   *argv);
	  }
	  else
	  {
	    if((stat(objTbl[objCnt].file, &stbuf) == 0) &&
	       ((stbuf.st_mode & S_IFREG) != 0))
	    {
	      objTbl[objCnt].size = stbuf.st_size;
	      ++objCnt;
	    }
	  }
	}
      }
      (void )closedir(dir);
      if(objCnt < 1)
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: Directory %s contains no readable Woolz files.\n",
		       *argv, path);
      }
    }
  }
  if(ok)
  {
    if((cache = AlcLRUCacheNew(maxObj, maxSz,
    			       (AlcLRUCKeyFn )WlzTstObjCacheKeyFn,
                               (AlcLRUCCmpFn )WlzTstObjCacheCmpFn,
			       (AlcLRUCUnlinkFn )WlzTstObjCacheUnlinkFn,
			       NULL)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to create cache.\n",
		     *argv);
    }
    else if(debug)
    {
      AlcLRUCacheFacts(cache, stderr);
    }
  }
  if(ok)
  {
    unsigned int i;

    if(random)
    {
      srand(0u);
    }
    for(i = 0; i < repeats; ++i)
    {
      unsigned int j;

      for(j = 0; j < objCnt; ++j)
      {
	unsigned int k,
		     hit,
		     colSz;

	k = (random)? rand() % objCnt: j; 
	errNum = WlzTstObjCacheAdd(cache, objTbl + k, &hit, &colSz);
	if(hit)
	{
	  ++nHit;
	}
	if(colSz > 0)
	{
	  if(colSz > 1)
	  {
	    ++nCol;
	  }
	  if(colSz > maxColSz)
	  {
	    maxColSz = colSz;
	  }
	}
	if(cache->numItem > maxTotCacheItems)
	{
	  maxTotCacheItems = cache->numItem;
	}
	if(cache->curSz > maxTotCacheEntSz)
	{
	  maxTotCacheEntSz = cache->curSz;
	}
	++nAdd;
	if(debug)
	{
	  AlcLRUCacheFacts(cache, stderr);
	}
      }
    }
    (void )printf("%s:\n"
                  "  calls = %u, hits = %u\n"
		  "  max cache items = %u\n"
		  "  max cache sz = %lu\n"
		  "  hash tbl sz = %u, num collisions %u, max list len %u\n",
                  *argv,
		  nAdd, nHit,
		  maxTotCacheItems, maxTotCacheEntSz,
		  cache->hashTblSz, nCol, maxColSz);
  }
  if(cache)
  {
    AlcLRUCacheFree(cache, 1);
  }
  if(objTbl)
  {
    int		i;

    for(i = 0; i < objCnt; ++i)
    {
      AlcFree(objTbl[i].file);
    }
    AlcFree(objTbl);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-d] [-h] [-n<max obj>] [-r<repeats>] [-s<max size>] [<dir>]\n"
    "Creates a least recent use removal cache and adds all Woolz objects in\n"
    "the given directory to it. Then accesses each Woolz object through a\n"
    "cycle of repeats. The cache statistics are output.\n"
    "Options are:\n"
    "  -d  Print cache debug facts to stderr (very noisy).\n"
    "  -h  Help, prints this usage message.\n"
    "  -n  Limit on the number of objects in the cache.\n"
    "  -r  Number of repeat accesses to all objects.\n"
    "  -m  Add objects in random order (rather than sequentially).\n"
    "  -s  Total cache entry size limit.\n",
    argv[0]);

  }
  return(!ok);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTst
* \brief	Adds an object to the cache if it is not already there.
* \param	cache			The cache.
* \param	ent			Cache entry.
* \param	hit			Hit flag set on cache hit, must not
* 					be NULL.
* \param	colSz			Hash collision size, must not be NULL.
*/
static WlzErrorNum WlzTstObjCacheAdd(AlcLRUCache *cache,
				     WlzTstObjCacheEntry *ent,
				     unsigned int *hit,
				     unsigned int *colSz)
{
  FILE		*fP;
  AlcLRUCItem	*item;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  *hit = 0;
  *colSz = 0;
  if((item = AlcLRUCEntryAdd(cache, ent->size, ent, NULL)) != NULL)
  {
    *colSz = AlcLRUCKeyGetNHashItem(cache, item->key);
    if(ent->obj)
    {
      *hit = 1;
    }
    else
    {
      if((fP = fopen(ent->file, "r")) != NULL)
      {
	ent->obj = WlzReadObj(fP, &errNum);
	(void )fclose(fP);
      }
    }
  }
  return(errNum);
}

/*!
* \return	Numeric key which identifies the object for the entry.
* \ingroup	WlzTst
* \brief	Computes a hash key from the file path of the given
* 		entry.
* \param	cache			The cache (not used).
* \param	e			Cast to (WlzTstObjCacheEntry *) to
* 					get the cache entry.
*/
static unsigned int WlzTstObjCacheKeyFn(AlcLRUCache *cache, const void *e)
{
  unsigned int	key;
  WlzTstObjCacheEntry *ent;

  ent = (WlzTstObjCacheEntry *)e;
  key = AlcStrSFHash(ent->file);
  return(key);
}

/*!
* \return	Non-zero value if the cache entries are different.
* \ingroup	WlzTst
* \brief	Compares the cache entry file paths and returns zero iff
* 		they are the same. Passed cache entry pointers are cast
* 		to (WlzTstObjCacheEntry *).
* 		entry.
* \param	e0			First cache entry pointer.
* \param	e1			Second cache entry pointer.
*/
static int	WlzTstObjCacheCmpFn(const void *e0, const void *e1)
{
  int		cmp;
  WlzTstObjCacheEntry *ent0,
  		   *ent1;

  ent0 = (WlzTstObjCacheEntry *)e0;
  ent1 = (WlzTstObjCacheEntry *)e1;
  cmp = strcmp(ent0->file, ent1->file);
  return(cmp);
}

/*!
* \ingroup	WlzTst
* \brief	A function called when a cache entry is about to be removed
* 		from the cache. This function just frees the entry object
* 		but if cache entries were not kept in a table the cache entry
* 		would be freed here to.
* \param	cache			The cache.
* \param	e			Cast to (WlzTstObjCacheEntry *) to get
* 					the cache entry.
*/
static void	WlzTstObjCacheUnlinkFn(AlcLRUCache *cache, const void *e)
{
  WlzTstObjCacheEntry *ent;

  if((ent = (WlzTstObjCacheEntry *)e) != NULL)
  {
    (void )WlzFreeObj(ent->obj);
    ent->obj = NULL;
  }
}
