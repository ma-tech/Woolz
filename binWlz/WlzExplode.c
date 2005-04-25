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
  int		objVecIdx,
  		objVecCount,
  		objIdx = 0,
		colourOnly = 0,
		freeObjAry = 0,
		useProps = 0,
		unconcatOnly = 0,
  		option,
		ok = 1,
		usage = 0;
  WlzProperty	prop;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*inFP = NULL,
  		*outFP = NULL;
  WlzObject	*obj,
  		*inObj = NULL;
  WlzObject	**objAry = NULL;
  WlzObject	*singleObjAry[1];
  WlzCompoundArray *cmpObj = NULL;
  char 		*name = NULL,
  		*outObjFileBodyStr = NULL,
		*outObjFileExtStr = NULL,
  		*inObjFileStr;
  const char	*errMsg;
  static char	optList[] = "b:e:cpuh",
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
      case 'c':
        colourOnly = 1;
        break;
      case 'p':
        useProps = 1;
	break;
      case 'b':
        outObjFileBodyStr = optarg;
	break;
      case 'e':
        outObjFileExtStr = optarg;
	break;
      case 'u':
        unconcatOnly = 1;
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
       ((inFP = (strcmp(inObjFileStr, "-")?
	        fopen(inObjFileStr, "r"): stdin)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to open input file %s (%s).\n",
		     *argv, inObjFileStr, errMsg);
    }
  }
  while(ok)
  {
    objAry = NULL;
    freeObjAry = 0;
    if((inObj= WlzAssignObject(WlzReadObj(inFP, &errNum), NULL)) == NULL)
    {
      ok = 0;
      if(0 == objIdx)
      {
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to read object from file %s (%s).\n",
		     *argv, inObjFileStr, errMsg);
      }
    }
    if(ok)
    {
      if(colourOnly)
      {
        cmpObj = WlzRGBAToCompound(inObj, WLZ_RGBA_SPACE_RGB, &errNum);
	if((NULL == cmpObj) || (WLZ_ERR_NONE != errNum))
	{
	  ok = 0;
	}
	else
	{
	  freeObjAry = 0;
	  objAry = cmpObj->o;
	  objVecCount = cmpObj->n;
	}
      }
      else if(unconcatOnly)
      {
	freeObjAry = 0;
	objVecCount = 1;
	objAry = singleObjAry;
	objAry[0] = WlzAssignObject(inObj, NULL);
      }
      else
      {
	switch( inObj->type )
	{
	  case WLZ_2D_DOMAINOBJ:
	    freeObjAry = 0;
	    objVecCount = 1;
	    objAry = singleObjAry;
	    objAry[0] = WlzAssignObject(inObj, NULL);
	    break;
	  case WLZ_3D_DOMAINOBJ:
	    freeObjAry = 1;
	    errNum = WlzExplode3D(&objVecCount, &objAry, inObj);
	    break;
	  case WLZ_COMPOUND_ARR_1:
	  case WLZ_COMPOUND_ARR_2:
	    freeObjAry = 0;
	    cmpObj = (WlzCompoundArray *)inObj;
	    objAry = cmpObj->o;
	    objVecCount = cmpObj->n;
	    break;
	  case WLZ_BOUNDLIST:
	    freeObjAry = 1;
	    errNum = WlzBoundaryToPolyObjArray(inObj, &objVecCount, &objAry);
	    break;
	  default:
	    errNum = WLZ_ERR_OBJECT_TYPE;
	    break;
	}
      }
      if(WLZ_ERR_NONE != errNum)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: Failed to explode object (%s)\n",
		       *argv, errMsg);
      }
    }
    if(ok)
    {
      objVecIdx = 0;
      if((colourOnly == 0) && (unconcatOnly == 0) &&
         (WLZ_3D_DOMAINOBJ ==  inObj->type))
      {
	objIdx = inObj->domain.p->plane1;
      }
      while((objVecIdx < objVecCount) && ok)
      {
	obj = *(objAry + objVecIdx);
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
		    "%s%06d.%s", outObjFileBodyStr, objIdx, outObjFileExtStr);
	  }
	  if((outFP = fopen(outObjFileStr, "w")) == NULL)
	  {
	    ok = 0;
	    (void )fprintf(stderr,
			   "%s: failed to open output file %s.\n",
			   *argv, outObjFileStr);
	  }
	}
	else
	{
	  outFP = stdout;
	}
	if(ok)
	{
	  if((errNum = WlzWriteObj(outFP, obj)) != WLZ_ERR_NONE)
	  {
	    ok = 0;
	    (void )WlzStringFromErrorNum(errNum, &errMsg);
	    (void )fprintf(stderr,
			   "%s: failed to write output object (%s).\n",
			   *argv, errMsg);
	  }
	}
	WlzFreeObj(obj);
	if(outObjFileBodyStr && outFP)
	{
	  (void )fclose(outFP);
	}
	++objIdx;
	++objVecIdx;
      }
      if(objAry && freeObjAry)
      {
        AlcFree(objAry);
      }
      if(cmpObj)
      {
        (void )WlzFreeObj((WlzObject *)cmpObj);
      }
    }
    if(inObj)
    {
      (void )WlzFreeObj(inObj);
    }
  }
  if((0 == ok) && (WLZ_ERR_READ_EOF == errNum) && (0 < objIdx))
  {
    ok = 1; 	/* Recover from file read error which is probably just the end
                   of files. */
  }
  if(inFP && strcmp(inObjFileStr, "-"))
  {
    (void )fclose(inFP);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-b<file body>] [-e <file extension>] [-h] [-p] [-u] [<in object>]\n"
    "Options:\n"
    "  -b  Output object file body.\n"
    "  -c  Colour objects exploded into RGBA components only.\n"
    "  -e  Output object file extension (default wlz).\n"
    "  -p  Use names from the object properties instead of a numerical\n"
    "      index for the file names. If the objects don't have properties\n"
    "      with names then the numerical index will still be used.\n"
    "  -u  Unconcatenate objects only. Objects will not be exploded, only\n"
    "      their concatenation will be.\n"
    "  -h  Help, prints this usage message.\n"
    "Explodes the following:\n"
    "  Concatenated input objects into separate objects.\n"
    "  Input 3D Woolz domain object into 2D domain objects.\n"
    "  Input compound object into its constituent objects.\n"
    "  Input boundary object into its constituent polylines.\n"
    "An output file extension can only be given if an output file body is\n"
    "also given.\n"
    "Objects are read from stdin unless a filename is given.\n",
    *argv,
    " -b plane myobj.wlz\n"
    "The input Woolz object is read from myobj.wlz (minimum plane coordinate\n"
    "30), exploded and written to plane000030.wlz, plane000031.wlz, etc.\n");
  }
  return(!ok);
}
