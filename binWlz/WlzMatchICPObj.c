#pragma ident "MRC HGU $Id$"
/*!
* \file        	WlzMatchICPObj.c
* \author       Bill Hill
* \date         March 2002
* \version	$Id$
* \note
*		Copyright:
*		2001 Medical Research Council, UK.
*		All rights reserved.
* \par  Address:
*		MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* \brief	Object matching using ICP based registration to
*		find corresponding points.
* \ingroup	WlzTransform
* \todo         
* \bug
*/
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <Wlz.h>
#include <strings.h>

extern int	getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int             main(int argc, char *argv[])
{
  int           idx,
		brkFlg = INT_MAX,
		dbgFlg = 0,
		nMatch,
		maxItr = 200,
		minSpx = 20,
  		option,
  		ok = 1,
		usage = 0;
  double	thrMeanD = 1.0,
  		thrMaxD = 3.0;
  FILE		*fP = NULL;
  char		*inTrFileStr,
  		*outFileStr,
		*outDbgFileStr;
  char		*inObjFileStr[2];
  WlzObject	*inTrObj;
  WlzObject	*inObj[2];
  WlzVertexP	matchTP,
  		matchSP;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  static char	optList[] = "aghrb:d:D:n:o:t:i:s:x:";
  const char	outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  inTrObj = NULL;
  inObj[0] = inObj[1] = NULL;
  outFileStr = (char *)outFileStrDef;
  outDbgFileStr = NULL;
  inTrFileStr = NULL;
  inObjFileStr[0] = inObjFileStr[1] = (char *)inObjFileStrDef;
  matchTP.v = matchSP.v = NULL;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'b':
        if((sscanf(optarg, "%d", &brkFlg) != 1) || (brkFlg < 0))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'd':
	dbgFlg = 1;
        outDbgFileStr = optarg;
	break;
      case 'D':
        if((sscanf(optarg, "%d", &dbgFlg) != 1) ||
	   (dbgFlg < 0) || (dbgFlg > 2))
        {
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'n':
        if((sscanf(optarg, "%lg", &thrMeanD) != 1) || (thrMeanD <= 0.0))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 't':
        inTrFileStr = optarg;
	break;
      case 'i':
        if((sscanf(optarg, "%d", &maxItr) != 1) || (maxItr <= 0))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 's':
        if((sscanf(optarg, "%d", &minSpx) != 1) || (minSpx <= 0))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'x':
        if((sscanf(optarg, "%lg", &thrMaxD) != 1) || (thrMaxD <= 0.0))
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
  if(ok)
  {
    if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
       (outFileStr == NULL) || (*outFileStr == '\0'))
    {
      ok = 0;
      usage = 1;
    }
    if(ok && (optind < argc))
    {
      if((optind + 2) != argc)
      {
        usage = 1;
	ok = 0;
      }
      else
      {
        inObjFileStr[0] = *(argv + optind);
        inObjFileStr[1] = *(argv + optind + 1);
      }
    }
  }
  if(ok)
  {
    idx = 0;
    while(ok && (idx < 2))
    {
      if((inObjFileStr[idx] == NULL) ||
	  (*inObjFileStr[idx] == '\0') ||
	  ((fP = (strcmp(inObjFileStr[idx], "-")?
		  fopen(inObjFileStr[idx], "r"): stdin)) == NULL) ||
	  ((inObj[idx] = WlzAssignObject(WlzReadObj(fP, &errNum),
	  				 NULL)) == NULL) ||
	  (errNum != WLZ_ERR_NONE))
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to read object from file %s\n",
		       *argv, inObjFileStr[idx]);
      }
      if(fP && strcmp(inObjFileStr[idx], "-"))
      {
	fclose(fP);
        fP = NULL;
      }
      ++idx;
    }
  }
  if(ok && inTrFileStr && (*inTrFileStr != '\0'))
  {
    if(((fP = (strcmp(inTrFileStr, "-")?
    	       fopen(inTrFileStr, "r"): stdin)) == NULL) ||
       ((inTrObj = WlzAssignObject(WlzReadObj(fP, &errNum),
       				   NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s Failed to read the initial affine transform from\n"
      		     "file %s.\n",
		     *argv, inTrFileStr);
    }
    if(fP && strcmp(inTrFileStr, "-"))
    {
      fclose(fP);
      fP = NULL;
    }
    if(inTrObj &&
       ((inTrObj->type != WLZ_AFFINE_TRANS) || (inTrObj->domain.core == NULL)))
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: Initial affine transform object invalid type\n",
		     *argv);
    }
  }
  if(ok)
  {
    errNum = WlzMatchICPObjs(inObj[0], inObj[1],
    			     (inTrObj)? inTrObj->domain.t: NULL,
			     &nMatch, &matchTP, &matchSP,
			     maxItr, minSpx, brkFlg, thrMeanD, thrMaxD);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr, "%s Failed to match contours.\n",
      		     argv[0]);
    }
  }
  if(ok)
  {
    if(outDbgFileStr && (dbgFlg > 0))
    {
      errNum = WLZ_ERR_WRITE_EOF;
      if(((fP = (strcmp(outDbgFileStr, "-")?
      	        fopen(outDbgFileStr, "w"):
		stdout)) == NULL))
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to open debug object "
		       "file %s (%s).\n",
		       *argv, outDbgFileStr, errMsg);
      }
      else if((errNum = WlzWriteObj(fP, inObj[1])) != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to write debug object "
		       "to file %s (%s).\n",
		       *argv, outDbgFileStr, errMsg);
      }
      if(fP && strcmp(outDbgFileStr, "-"))
      {
        (void )fclose(fP);
      }
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outFileStr, "-")?
              fopen(outFileStr, "w"):
	      stdout)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to open correspondences "
		     "file %s (%s).\n",
		     *argv, outFileStr, errMsg);
    }
    else
    {
      /* Output the correspondces: Assumed to be WlzDVertex2. */
      for(idx = 0; idx <nMatch; ++idx)
      {
        fprintf(fP, "%g %g %g %g\n",
		(matchSP.d2 + idx)->vtX, (matchSP.d2 + idx)->vtY,
		(matchTP.d2 + idx)->vtX - (matchSP.d2 + idx)->vtX,
		(matchTP.d2 + idx)->vtY - (matchSP.d2 + idx)->vtY);
      }
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  AlcFree(matchTP.v);
  AlcFree(matchSP.v);
  (void )WlzFreeObj(inTrObj);
  (void )WlzFreeObj(inObj[0]);
  (void )WlzFreeObj(inObj[1]);
  if(usage)
  {
      (void )fprintf(stderr,
      "Usage: %s%s",
      *argv,
      " [-b#] [-d#] [-D#] [-o#] [-t#] [-i#] [-h]\n"
      "          [<input object 0>] [<input object 1>]\n"
      "Options:\n"
      "  -b  Shell breaking control value.\n"
      "  -d  Debug output file.\n"
      "  -D  Debug flag:\n"
      "        0  no debug output.\n"
      "        1  untransformed decomposed source model.\n"
      "        2  transformed decomposed source model.\n"
      "  -h  Prints this usage information.\n"
      "  -n  Threshold mean source normal vertex distance for fragment to\n"
      "      be used in correspondence computation.\n"
      "  -o  Output file name.\n"
      "  -t  Initial affine transform.\n"
      "  -i  Maximum number of iterations.\n"
      "  -s  Miimum number of simplicies.\n"
      "  -x  Threshold maximum source normal vertex distance for fragment to\n"
      "      be used in correspondence computation.\n"
      "Reads a pair of Woolz objects, which idealy are contours and then\n"
      "computes a set of correspondence points using an ICP based matching\n"
      "algorithm.\n");
  }
  return(!ok);
}
