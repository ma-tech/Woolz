#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzBoundaryToObj.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Filter to convert a boundary list to the corresponding
*		domain object.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <Wlz.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void ShowUsage(char *arg)
{
  fprintf(stderr,
	  "Usage: %s [-f s|e|o] [-h] [<input file>]\n"
	  "Options:\n"
	  "  -f Sets the boundary fill method to be one of: simple (s),\n"
	  "     even/odd (e) or vertex (v) fill, default is even/odd.\n"
	  "  -h Help - prints this usage message\n"
	  "Converts a boundary list woolz object to a corresponding domain\n"
	  "using the given fill method.\n",
	  arg);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *nobj;
  FILE		*inFile;
  char 		optList[] = "f:h";
  int		option,
  		usage = 0;
  WlzPolyFillMode fill = WLZ_EVEN_ODD_FILL;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  if((option = getopt(argc, argv, optList)) != EOF)
  {
    switch(option)
    {
      case 'f':
	switch(*optarg)
	{
	  case 'e':
	    fill = WLZ_EVEN_ODD_FILL;
	    break;
	  case 's':
	    fill = WLZ_SIMPLE_FILL;
	    break;
	  case 'v':
	    fill = WLZ_VERTEX_FILL;
	    break;
	  default:
	    usage = 1;
	    break;
	}
	break;
      case 'h': /* FALLTHROUGH */
      default:
	usage = 1;
	break;
    }
  }

  if(usage)
  {
    ShowUsage(argv[0]);
    return 1;
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      ShowUsage(argv[0]);
      return 1;
    }
  }

  /* Read objects and convert. */
  while((obj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) != NULL) 
  {
    switch( obj->type )
    {
    case WLZ_BOUNDLIST:
      if((nobj = WlzBoundaryToObj(obj, fill, NULL)) != NULL)
      {
	(void )WlzWriteObj(stdout, nobj);
	WlzFreeObj(nobj);
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      if(obj->domain.p->type == WLZ_PLANEDOMAIN_BOUNDLIST)
      {
	if((nobj = WlzBoundaryToObj(obj, fill, NULL)) != NULL)
	{
	  (void )WlzWriteObj(stdout, nobj);
	  WlzFreeObj(nobj);
	}
      }
      else {
	(void )WlzWriteObj(stdout, obj);
      }
      break;
	
    default:
      (void )WlzWriteObj(stdout, obj);
      break;
    }

    WlzFreeObj(obj);
  }
  return WLZ_ERR_NONE;
}
