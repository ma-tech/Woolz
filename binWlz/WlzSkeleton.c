#ifndef DOXYGEN_SHOULD_SKIP_THIS
#pragma ident "MRC HGU %W%\t%G% bill@hgu.mrc.ac.uk"
/************************************************************************
* Project:      Woolz							*
* Title:        WlzSkeleton.c			                      	*
* Date:         November 1998	                                    	*
* Author:       Bill Hill 				    		*
* Copyright:	1998 Medical Research Council, UK.			*
*		All rights reserved.					*
* Address:	MRC Human Genetics Unit,				*
*		Western General Hospital,				*
*		Edinburgh, EH4 2XU, UK.					*
* Version:	%I%							*
* Purpose:      Woolz filter which computes the skeleton of the		*
*		given objects domain.					*
* Maintenance:	Log changes below, with most recent at top of list.	*
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		option,
		cValI,
		smPass,
		ok = 1,
		usage = 0;
  WlzConnectType con;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	defCValI = 8,
  		defSmPass = 1;
  const char    *errMsg;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL,
		*outObj = NULL;
  char 		*outObjFileStr,
  		*inObjFileStr;
  static char	optList[] = "o:c:s:h",
		outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  cValI = defCValI;
  smPass = defSmPass;
  outObjFileStr = outObjFileStrDef;
  inObjFileStr = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'c':
	if(sscanf(optarg, "%d", &cValI) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 's':
	if(sscanf(optarg, "%d", &smPass) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'h':
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
     (outObjFileStr == NULL) || (*outObjFileStr == '\0'))
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
  if(ok)
  {
    switch(cValI)
    {
      case 4:
	con = WLZ_4_CONNECTED;
	break;
      case 6:
	con = WLZ_6_CONNECTED;
	break;
      case 8:
	con = WLZ_8_CONNECTED;
	break;
      case 18:
	con = WLZ_18_CONNECTED;
	break;
      case 26:
	con = WLZ_26_CONNECTED;
	break;
      default:
	usage = 1;
	ok = 0;
	break;
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_READ_EOF;
    if((inObjFileStr == NULL) ||
       (*inObjFileStr == '\0') ||
       ((fP = (strcmp(inObjFileStr, "-")?
	      fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to read object from file %s (%s).\n",
		     *argv, inObjFileStr, errMsg);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok)
  {
    if((outObj = WlzSkeleton(inObj, smPass, con, &errNum)) == NULL)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to compute skeleton (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outObjFileStr, "-")?
              fopen(outObjFileStr, "w"):
	      stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to write output object (%s).\n",
		     *argv, errMsg);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%d%s%d%sExample: %s%s",
    *argv,
    " [-o<out object>] [-c #] [-s #] [-h] [<in object>]\n"
    "Options:\n"
    "  -o  Output object file name.\n"
    "  -c  minimum connectivity (4, 6, 8, 18, 26), default ",
    defCValI,
    "\n."
    "  -s  smoothing passes, default ",
    defSmPass,
    ".\n"
    "  -h  Help, prints this usage message.\n"
    "Computes the skeleton of a 2D domain object.\n"
    "Objects are read from stdin and written to stdout unless the filenames\n"
    "are given.\n",
    *argv,
    " -o skel.wlz myobj.wlz\n"
    "The input Woolz object is read from myobj.wlz, filtered and written\n"
    "to skel.wlz.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
