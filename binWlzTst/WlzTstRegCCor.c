#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTstRegCCor_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzTst/WlzTstRegCCor.c
* \author       Bill Hill
* \date         October 2007
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
* \brief	Tests for registration using frequency
* 		domain cross correlation.
* \ingroup	BinWlzTst
*/

#include <float.h>
#include <limits.h>
#include <string.h>
#include <Wlz.h>

extern int	getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int             main(int argc, char *argv[])
{
  int           idx,
  		option,
		maxItr = 10,
		noise = 0,
		inv = 0,
		conv = 0,
  		ok = 1,
		usage = 0;
  double	cCor,
  		maxRot = 45 * (2.0 * WLZ_M_PI / 360.0);
  FILE		*fP = NULL;
  char		*outObjFileStr;
  char		*inObjFileStr[2];
  WlzDVertex2	maxTran;
  WlzWindowFnType winFn = WLZ_WINDOWFN_NONE;
  WlzTransformType trType = WLZ_TRANSFORM_2D_REG;
  WlzValues	nullVal;
  WlzDomain	outDom;
  WlzObject	*outObj;
  WlzObject	*inObj[2];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  static char	optList[] = "hintro:w:";
  const char	outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  maxTran.vtX = maxTran.vtY = 200.0;
  inObj[0] = inObj[1] = NULL;
  nullVal.core = NULL;
  outDom.core = NULL;
  outObj = NULL;
  outObjFileStr = (char *)outFileStrDef;
  inObjFileStr[0] = inObjFileStr[1] = (char *)inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'h':
        usage = 1;
	ok = 0;
	break;
      case 'i':
        inv = 1;
	break;
      case 'n':
        noise = 1;
	break;
      case 'r':
        trType = WLZ_TRANSFORM_2D_REG;
	break;
      case 't':
        trType = WLZ_TRANSFORM_2D_TRANS;
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'w':
	if((winFn = WlzWindowFnValue(optarg)) == WLZ_WINDOWFN_UNSPECIFIED)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if(ok)
  {
    if((inObjFileStr[0] == NULL) || (*inObjFileStr[0] == '\0') ||
       (inObjFileStr[1] == NULL) || (*inObjFileStr[1] == '\0') ||
       (outObjFileStr == NULL) || (*outObjFileStr == '\0'))
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
      }
      ++idx;
    }
  }
  if(ok)
  {
    outDom.t = WlzRegCCorObjs(inObj[0], inObj[1], NULL, trType,
			      maxTran, maxRot, maxItr,
			      winFn, noise, inv,
			      &conv, &cCor, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s Failed to register objects (%s).\n",
      		     argv[0], errMsg);
    }
  }
  if(ok)
  {
    outObj = WlzMakeMain(WLZ_AFFINE_TRANS, outDom, nullVal, NULL, NULL,
    			 &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to make affine transform object (%s).\n",
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
      		     "%s: failed to write output affine transform object "
		     "to file %s (%s).\n",
		     *argv, outObjFileStr, errMsg);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  if(ok)
  {
    (void )printf("%g\n", cCor);
  }
  WlzFreeObj(inObj[0]);
  WlzFreeObj(inObj[1]);
  if(outObj)
  {
    (void )WlzFreeObj(outObj);
  }
  else if(outDom.core)
  {
    (void )WlzFreeAffineTransform(outDom.t);
  }
  if(usage)
  {
      (void )fprintf(stderr,
      "Usage: %s%s",
      *argv,
      " [-o<output object>] [-h] [-o<output file>] [-i] [-n] [-r] [-t]\n"
      "              [-w<window fn>] [<input object 0>] [<input object 1>]\n"
      "Options:\n"
      "  -h  Prints this usage information.\n"
      "  -i  Invert both the input objects values.\n"
      "  -n  Use Gaussian noise.\n"
      "  -r  Find registration transform.\n"
      "  -t  Find tranlation only transform.\n"
      "  -w  Window function name, this must be one of: blackman,\n"
      "      hamming, hanning, none, parzen, rectangle or welch.\n"
      "  -o  Output transform object file name.\n"
      "Computes a registration transform which registers the second of the\n"
      "given objects to the first using a cross correlation algorithm.\n");
  }
  return(!ok);
}
