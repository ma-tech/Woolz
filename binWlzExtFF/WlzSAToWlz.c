#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzSAToWlz_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzExtFF/WlzSAToWlz.c
* \author       Richard Baldock
* \date         March 2005
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
* \brief        Convert the Smart Atlas polygons to a Woolz
* 		boundaries.
* \ingroup	BinWlzExtFF
*
* \par Binary
* \ref wlzsatowlz "WlzSAToWlz"
*/

/*!
\ingroup BinWlzExtFF
\defgroup wlzsatowlz WlzSAToWlz
\par Name
WlzSAToWlz - converts a Smart Atlas XML file and convert to Woolz boundaries.
\par Synopsis
\verbatim
WlzSAToWlz  [-h] [-v] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Prints usage information.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Verbose operation.</td>
  </tr>
</table>
\par Description
WlzXXX.
\par Examples
\verbatim
\endverbatim
\par File
\ref WlzFacts.c "WlzFacts.c"
\par See Also
\ref BinWlzExtFF "WlzIntro(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <Wlz.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

typedef struct _NamedBndItem {
  char		*name;
  WlzBoundList	*bnd;
} NamedBndItem;

void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s -h -v [<input file>]\n"
	  "\tRead the Smart Atlas XML file and convert to woolz\n"
	  "Version: %s\n"
	  "Options:\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -v        Verbose operation\n",
	  proc_str,
	  WlzVersion());
  return;
}

char *nextStartTag(
  FILE	*inFile)
{
  char 	*tag=NULL;
  char 	buf[256];
  char	c;
  int	i;

  while( (c = fgetc(inFile)) != EOF ){
    if( c == '<' ){
      i = 0;
      c = fgetc(inFile);
      if( c == EOF ){
	break;
      }
      else if( c != '/' ){
	buf[i++] = c;
	while( (c = fgetc(inFile)) != '>' ){
	  if( c == EOF ){
	    break;
	  }
	  buf[i++] = c;
	}
	if( c == '>' ){
	  buf[i] = '\0';
	  break;
	}
      }
    }
  }

  if( i ){
    tag = AlcStrDup(buf);
  }
  return tag;
}
  

char *getNextName(
  FILE	*inFile)
{
  char	*name=NULL;
  char	*tag;
  char	c;
  char	buffer[128];
  int	i;

  /* check file pointer */
  if( inFile == NULL ){
    return NULL;
  }

  /* search for name tag */
  while((tag = nextStartTag(inFile)) != NULL){
    if(strcmp(tag, "name") == 0){
      break;
    }
  }

  /* read until end tag */
  if( tag ){
    i = 0;
    while( (c = fgetc(inFile)) != '<' ){
      if( c == EOF ){
	break;
      }
      buffer[i++] = c;
    }
    buffer[i] = '\0';
    while( (c = fgetc(inFile)) != '>' ){
      if( c == EOF ){
	break;
      }
    }

    /* extract name */
    name = AlcStrDup(buffer);
  }

  return name;
}

WlzBoundList *getNextBoundary(
  FILE	*inFile)
{
  WlzBoundList	*bnd=NULL, *bnd1=NULL;
  WlzPolygonDomain	*poly;
  char	c, *buf, *str, *tag;
  long	start, end;
  int	i, nVtxs;
  WlzDVertex2	*vtxs;

  /* find the polyline start tag */
  while((tag = nextStartTag(inFile)) != NULL){
    if(strcmp(tag, "c") == 0){
      break;
    }
  }

  /* count characters and vertices */
  start = ftell(inFile);
  nVtxs = 0;
  while( (c = fgetc(inFile)) != '<' ){
    if( c == ',' ){
      nVtxs++;
    }
  }
  end = ftell(inFile);

  /* read vertex values and consume end tag */
  fseek(inFile, start, SEEK_SET);
  buf = (char *) AlcMalloc(sizeof(char) * (end-start+1));
  fread(buf, end-start-1, sizeof(char), inFile);
  while( (c = fgetc(inFile)) != '>' ){}

  /* now read the vertices */
  vtxs = (WlzDVertex2 *) AlcMalloc(sizeof(WlzDVertex2) * (nVtxs + 1));
  str = buf;
  for(i=0; i < (nVtxs-1); i++){
    sscanf(str, "%lf,%lf", &(vtxs[i].vtX), &(vtxs[i].vtY));
    while( *str != ' ' ){
      str++;
    }
    str++;
  }
  sscanf(str, "%lf,%lf", &(vtxs[i].vtX), &(vtxs[i].vtY));
  i++;
  /* close the boundary */
  vtxs[i] = vtxs[0];
  nVtxs++;

  /* create a boundlist and polyline structure */
  poly = WlzMakePolygonDomain(WLZ_POLYGON_DOUBLE,
			      nVtxs, (WlzIVertex2 *) vtxs,
			      nVtxs, 0, NULL);
  bnd = WlzMakeBoundList(WLZ_BOUNDLIST_PIECE, 1, poly, NULL);

  /* need to check if there is another boundary */
  start = ftell(inFile);
  fread(buf, 3, sizeof(char), inFile);
  fseek(inFile, start, SEEK_SET);
  if( strncmp(buf, "<c>", 3) == 0 ){
    if((bnd1 = getNextBoundary(inFile)) != NULL){
      /* check for down or next */
      vtxs = (WlzDVertex2 *) bnd1->poly;
      if( WlzInsidePolyEOD(*vtxs, poly, NULL) ){
	bnd1->type = WLZ_BOUNDLIST_HOLE;
	bnd->down = WlzAssignBoundList(bnd1, NULL);
      }
      else {
	bnd->next = WlzAssignBoundList(bnd1, NULL);
      }
    }
  }

  /* clear space */
  AlcFree(buf);

  return bnd;
}

void addToBndList(
  AlcDLPList	*list,
  char		*name,
  WlzBoundList	*bnd)
{
  AlcDLPItem	*bndItem;
  NamedBndItem	*namedBndItem;
  WlzBoundList	*tmpBnd;

  /* check if name exists */
  bndItem = list->head;
  namedBndItem = NULL;
  while( bndItem ){
    namedBndItem = (NamedBndItem *) bndItem->entry;
    if( strcmp(name, namedBndItem->name) == 0 ){
      break;
    }
    namedBndItem = NULL;
    bndItem = bndItem->next;
    if( bndItem == list->head ){
      break;
    }
  }

  if( namedBndItem ){
    tmpBnd = namedBndItem->bnd;
    while(tmpBnd->next){
      tmpBnd = tmpBnd->next;
    }
    tmpBnd->next = WlzAssignBoundList(bnd, NULL);
  }
  else {
    /* create a NamedBndItem */
    namedBndItem = (NamedBndItem *) AlcMalloc(sizeof(NamedBndItem));
    namedBndItem->name = name;
    namedBndItem->bnd = bnd;

    /* add to the list */
    AlcDLPListEntryAppend(list, NULL, (void *) namedBndItem, NULL);
  }

  return;
}

int		main(int argc, char *argv[])
{
  FILE		*inFile, *outFile;
  char 		optList[] = "hv";
  int		option;
  int		verboseFlg=0;
  WlzObject	*obj;
  WlzDomain	domain;
  WlzValues	values;
  WlzBoundList	*tmpBnd;
  char		*name;
  AlcDLPList	*bndDLPList;
  AlcDLPItem	*bndItem;
  NamedBndItem	*namedBndItem;
  int		noNameCntr=0;
  char		noNameBuf[32];

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'v':
      verboseFlg = 1;
      break;

    case 'h':
      usage(argv[0]);
      return WLZ_ERR_NONE;

    default:
      usage(argv[0]);
      return WLZ_ERR_UNSPECIFIED;
    }
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return( 1 );
    }
  }

  /* read XML file building boundary lists */
  /* quick fix read  - read until no new names */
  bndDLPList = AlcDLPListNew(NULL);
  name = NULL;
  while((name = getNextName(inFile)) != NULL){
    /* check for empty name */
    if( strcmp(name, " ") == 0 ){
      noNameCntr++;
      sprintf(noNameBuf, "NoName_%03d", noNameCntr);
      name = AlcStrDup(noNameBuf);
    }
    /* convert comma to underscore */
    tmpBnd = getNextBoundary(inFile);
    addToBndList(bndDLPList, name, tmpBnd);
  }

  /* write boundary lists */
  bndItem = bndDLPList->head;
  while( bndItem ){
    namedBndItem = (NamedBndItem *) bndItem->entry;
    name = (char *) AlcMalloc(sizeof(char)*(strlen(namedBndItem->name)+16));
    sprintf(name, "%s.wlz", namedBndItem->name);
    if((outFile = fopen(name, "w")) != NULL){
      domain.b = namedBndItem->bnd;
      values.core = NULL;
      obj = WlzMakeMain(WLZ_BOUNDLIST, domain, values, NULL, NULL, NULL);
      WlzWriteObj(outFile, obj);
      WlzFreeObj(obj);
      fclose(outFile);
    }
    AlcFree(name);
    bndItem = bndItem->next;
    if( bndItem == bndDLPList->head ){
      break;
    }
  }

  exit(0);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
