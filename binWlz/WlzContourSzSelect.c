#ifndef DOXYGEN_SHOULD_SKIP_THIS
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        WlzContourSzSelect.c
* Date:         March 2001
* Author:       Bill Hill
* Copyright:	2001 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Filters to remove small shells from contours.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <Wlz.h>
#include <string.h>


extern int	getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int             main(int argc, char *argv[])
{
  int           idx,
  		option,
  		ok = 1,
		usage = 0,
		smShellSz = -1;
  FILE		*fP = NULL;
  char		*inObjFileStr,
		*outObjFileStr;
  WlzObject	*obj;
  WlzDomain	ctrDom;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char    *errMsgStr;
  static char	optList[] = "ho:s:";
  const char	outObjFileStrDef[] = "-",
		inObjFileStrDef[] = "-";

  obj = NULL;
  ctrDom.core = NULL;
  inObjFileStr = (char *)inObjFileStrDef;
  outObjFileStr = (char *)outObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
	outObjFileStr = optarg;
	break;
      case 's':
        if(sscanf(optarg, "%d", &smShellSz) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'h':
        usage = 1;
	ok = 0;
	break;
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if(smShellSz <= 1)
  {
    ok = 0;
    usage = 1;
  }
  if(ok)
  {
    if((inObjFileStr == NULL) || (*inObjFileStr == '\0'))
    {
      ok = 0;
      usage = 1;
    }
    if(ok && (optind < argc))
    {
      if((optind + 1) != argc)
      {
        usage = 1;
	ok = 0;
      }
      else
      {
        inObjFileStr = *(argv + optind);
      }
    }
  }
  if(ok)
  {
    if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
	((fP = (strcmp(inObjFileStr, "-")?
		fopen(inObjFileStr, "r"): stdin)) == NULL) ||
	((obj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
	(errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: failed to read object from file %s\n",
		     *argv, inObjFileStr);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      (void )fclose(fP);
      fP = NULL;
    }
  }
  if(ok)
  {
    if(obj->type != WLZ_CONTOUR)
    {
      ok = 0;
      (void )fprintf(stderr, "%s Input object isn't a contour.\n",
		     argv[0]);
    }
  }
  if(ok)
  {
    errNum = WlzGMFilterRmSmShells(obj->domain.ctr->model, smShellSz);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
		     "%s: failed to filter small shells from contour (%s).\n",
		     *argv, errMsgStr);
    }
  }
  if(ok)
  {
  if((fP = (strcmp(outObjFileStr, "-")?
           fopen(outObjFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open output file %s.\n",
                     argv[0], outObjFileStr);
    }
  }
  if(ok)
  {
    errNum = WlzWriteObj(fP, obj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to write output object (%s).\n",
                     argv[0], errMsgStr);
    }
  }
  if(fP)
  {
    (void )fclose(fP);
    fP = NULL;
  }
  if(obj)
  {
    (void )WlzFreeObj(obj);
  }
  if(usage)
  {
      (void )fprintf(stderr,
      "Usage: %s %s",
      *argv,
      " [-h] -s<threshold> [<input object>]\n"
      "Options:\n"
      "  -s  Threshold shell size (edges in 2D, loops in 3D).\n"
      "  -h  Prints this usage information.\n"
      "Removes all shells in the contour which have less than the threshold\n"
      "number of elements.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
