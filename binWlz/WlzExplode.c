#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzExplode.c
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2002 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Woolz filter to explode 3D and compound objects into
*		2D objects.
* \todo         -
* \bug          None known.
*/
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
  int		objPlane,
  		objVecIdx,
  		objVecCount,
		useProps = 0,
  		option,
		ok = 1,
		usage = 0;
  WlzProperty	prop;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*obj,
  		*inObj = NULL;
  WlzObject	**objVec = NULL;
  WlzCompoundArray	*cmpObj;
  char 		*name = NULL,
  		*outObjFileBodyStr = NULL,
		*outObjFileExtStr = NULL,
  		*inObjFileStr;
  const char	*errMsg;
  static char	optList[] = "b:e:ph",
		outObjFileStr[FILENAME_MAX],
		outObjFileExtStrDef[] = "wlz",
		outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  inObjFileStr = inObjFileStrDef;
  outObjFileExtStr = outObjFileExtStrDef;
  strcpy(outObjFileStr, outObjFileStrDef);
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'p':
        useProps = 1;
	break;
      case 'b':
        outObjFileBodyStr = optarg;
	break;
      case 'e':
        outObjFileExtStr = optarg;
	break;
      case 'h':
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if((inObjFileStr == NULL) || (*inObjFileStr == '\0'))
  {
    ok = 0;
    usage = 1;
  }
  if(outObjFileBodyStr && outObjFileExtStr &&
     (strlen(outObjFileExtStr) +
      strlen(outObjFileExtStr)) > (FILENAME_MAX - 24))
  {
    ok = 0;
    (void )fprintf(stderr,
    		   "%s: output filename length exceeded!\n",
		   *argv);
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
    errNum = WLZ_ERR_READ_EOF;
    if((inObjFileStr == NULL) ||
       (*inObjFileStr == '\0') ||
       ((fP = (strcmp(inObjFileStr, "-")?
	      fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((inObj= WlzReadObj(fP, &errNum)) == NULL))
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
    switch( inObj->type )
    {
    case WLZ_COMPOUND_ARR_1:
    case WLZ_COMPOUND_ARR_2:
      cmpObj = (WlzCompoundArray *) inObj;
      objVec = cmpObj->o;
      objVecCount = cmpObj->n;
      break;

    case WLZ_BOUNDLIST:
      if((errNum = WlzBoundaryToPolyObjArray(inObj, &objVecCount, &objVec)) !=
	 WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to explode object (%s)\n",
		       *argv, errMsg);
      }
      break;

    default:
      if((errNum = WlzExplode3D(&objVecCount, &objVec, inObj)) !=
	 WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to explode object (%s)\n",
		       *argv, errMsg);
      }
      break;
    }
  }
  if(ok)
  {
    objVecIdx = 0;
    objPlane = 0;
    if( inObj->type == WLZ_3D_DOMAINOBJ )
    {
      objPlane = inObj->domain.p->plane1;
    }
    while((objVecIdx < objVecCount) && ok)
    {
      obj = *(objVec + objVecIdx);
      if(outObjFileBodyStr)
      {
	if(useProps && obj->plist)
	{
	  prop = WlzGetProperty(obj->plist->list, WLZ_PROPERTY_NAME, NULL);
	  if(prop.core == NULL)
	  {
	    prop = WlzGetProperty(obj->plist->list, WLZ_PROPERTY_GREY, NULL);
	  }
	  if(prop.core)
	  {
	    switch(prop.core->type)
	    {
	      case WLZ_PROPERTY_NAME:
		name = prop.name->name;
	        break;
	      case WLZ_PROPERTY_GREY:
		name = prop.greyV->name;
	        break;
	    }
	  }
	}
	if(name)
	{
	  (void )sprintf(outObjFileStr,
	  	  "%s%s.%s", outObjFileBodyStr, name, outObjFileExtStr);
	}
	else
	{
	  (void )sprintf(outObjFileStr,
		  "%s%06d.%s", outObjFileBodyStr, objPlane, outObjFileExtStr);
	}
        if((fP = fopen(outObjFileStr, "w")) == NULL)
	{
	  ok = 0;
	  (void )fprintf(stderr,
	  		 "%s: failed to open output file %s.\n",
			 *argv, outObjFileStr);
	}
      }
      else
      {
        fP = stdout;
      }
      if(ok)
      {
        if((errNum = WlzWriteObj(fP, obj)) != WLZ_ERR_NONE)
	{
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
	  		 "%s: failed to write output object (%s).\n",
			 *argv, errMsg);
	}
      }
      WlzFreeObj(obj);
      if(outObjFileBodyStr && fP)
      {
        fclose(fP);
      }
      ++objPlane;
      ++objVecIdx;
    }
    AlcFree(objVec);
  }
  if(inObj)
  {
    WlzFreeObj(inObj);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-b<file body>] [-e <file extension>] [-h] [-p] [<in object>]\n"
    "Options:\n"
    "  -b  Output object file body.\n"
    "  -e  Output object file extension (default wlz).\n"
    "  -p  Use names from the object properties instead of a numerical\n"
    "      index for the file names. If the objects don't have properties\n"
    "      with names then the numerical index will still be used.\n"
    "  -h  Help, prints this usage message.\n"
    "Explodes the input 3D Woolz domain object into 2D domain objects or\n"
    "explodes an input compound object into its constituent objects or\n"
    "explodes the input boundary object into its constituent polylines.\n"
    "An output file extension can only be given if an output file body is\n"
    "also given.\n"
    "Objects are read from stdin and written to stdout unless the filenames\n"
    "are given.\n",
    *argv,
    " -b plane myobj.wlz\n"
    "The input Woolz object is read from myobj.wlz (minimum plane coordinate\n"
    "30), exploded and written to plane000030.wlz, plane000031.wlz, etc.\n");
  }
  return(!ok);
}
