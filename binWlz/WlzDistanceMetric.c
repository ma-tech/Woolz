#ifndef DOXYGEN_SHOULD_SKIP_THIS
#pragma ident "MRC HGU $Id$"
/*!
* \file         binWlz/WlzDistanceMetric.c
* \author       Bill Hill
* \date         August 2003
* \version      $Id$
* \note
*               Copyright
*               2003 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Computes distance metrics between objects.
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
  int		idx,
		hFlag = 0,
		mFlag = 0,
		nFlag = 0,
		iFlag = 0,
		option,
		ok = 1,
		usage = 0;
  double	hDist,
  		mDist,
		nDist,
		iDist;
  char		hStr[64],
  		mStr[64],
		nStr[64],
		iStr[64];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*inObj[2];
  FILE		*fP = NULL;
  char 		*outFileStr;
  char  	*inObjFileStr[2];
  const char	*errMsg;
  static char	optList[] = "ho:HMNI",
		fileStrDef[] = "-";

  opterr = 0;
  inObj[0] = NULL;
  inObj[1] = NULL;
  inObjFileStr[0] = inObjFileStr[1] = outFileStr = fileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outFileStr = optarg;
	break;
      case 'H':
        hFlag = 1;
	break;
      case 'M':
        mFlag = 1;
	break;
      case 'N':
        nFlag = 1;
	break;
      case 'I':
        iFlag = 1;
	break;
      case 'h':
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if((inObjFileStr[0] == NULL) || (*inObjFileStr[0] == '\0') ||
     (inObjFileStr[1] == NULL) || (*inObjFileStr[1] == '\0') ||
     (outFileStr == NULL) || (*outFileStr == '\0'))
  {
    ok = 0;
    usage = 1;
  }
  if(ok && (optind < argc))
  {
    idx = 0;
    while((idx < 2) && (optind < argc))
    {
      inObjFileStr[idx] = *(argv + optind);
      ++optind;
      ++idx;
    }
  }
  if(ok && (optind != argc))
  {
    usage = 1;
    ok = 0;
  }
  if(ok)
  {
    /* Set default to all metrics if none have been specified. */
    if((hFlag == 0) && (mFlag == 0) && (nFlag == 0) && (iFlag == 0))
    {
      hFlag = mFlag = nFlag = iFlag = 1;
    }
    /* Read objects. */
    idx = 0;
    while((errNum == WLZ_ERR_NONE) && (idx < 2))
    {
      errNum = WLZ_ERR_READ_EOF;
      if((inObjFileStr[idx] == NULL) ||
	  (*inObjFileStr[idx] == '\0') ||
	  ((fP = (strcmp(inObjFileStr[idx], "-")?
		  fopen(inObjFileStr[idx], "r"): stdin)) == NULL) ||
	  ((inObj[idx] = WlzAssignObject(WlzReadObj(fP,
	  					    &errNum), NULL)) == NULL))
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: Failed to read object %d from file %s (%s)\n",
		       *argv, idx, inObjFileStr[idx], errMsg);
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
    /* Check object types. */
    if(inObj[0]->type != inObj[1]->type)
    {
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
    else if((inObj[0]->domain.core == NULL) ||
            (inObj[1]->domain.core == NULL))
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    else
    {
      switch(inObj[0]->type)
      {
        case WLZ_CONTOUR:
          if((inObj[0]->domain.core == NULL) ||
             (inObj[0]->domain.core == NULL) ||
             (inObj[0]->domain.ctr->model == NULL) ||
             (inObj[1]->domain.ctr->model == NULL))
          {
            errNum = WLZ_ERR_DOMAIN_NULL;
          }
          else
          {
	    errNum = WlzDistMetricGM(inObj[0]->domain.ctr->model,
	    			     inObj[1]->domain.ctr->model,
				     (hFlag)? &hDist: NULL,
				     (mFlag)? &mDist: NULL,
				     (nFlag)? &nDist: NULL,
				     (iFlag)? &iDist: NULL);
          }
          break;
        default:
          errNum = WLZ_ERR_OBJECT_TYPE;
          ok = 0;
          break;
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: Failed to compute distance metric (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    if(((fP = (strcmp(outFileStr, "-")?
              fopen(outFileStr, "w"): stdout)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed open output file (%s).\n",
		     *argv, errMsg);
    }
    else
    {
      hStr[0] = mStr[0] = nStr[0] = iStr[0] = '\0';
      if(hFlag)
      {
        (void )sprintf(hStr, " %g", hDist);
      }
      if(mFlag)
      {
        (void )sprintf(mStr, " %g", mDist);
      }
      if(nFlag)
      {
        (void )sprintf(nStr, " %g", nDist);
      }
      if(iFlag)
      {
        (void )sprintf(iStr, " %g", iDist);
      }
      (void )fprintf(fP, "%s%s%s%s\n", hStr, mStr, nStr, iStr);
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      fclose(fP);
    }
  }
  (void )WlzFreeObj(inObj[0]);
  (void )WlzFreeObj(inObj[1]);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-h] [-o <output file>] [-H] [-M] [-N] [-I]\n"
    "                         [<in obj 0>] [<in obj 1>]\n"
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -o  Output file for distance metrics.\n"
    "  -H  Hausdorff distance metric.\n"
    "  -M  Mean of nearest neighbours distance metrics.\n"
    "  -N  Median of nearest neighbours distance metrics.\n"
    "  -I  Minimum nearest neighbour distance.\n"
    "Computes distance metrics for the given pair of objects.\n"
    "If no distance metrics are specified on the command line then all\n"
    "the metrics will be be computed and output, but if any metrics are\n"
    "specified on the computed line then only those metrics will be\n"
    "computed and output.\n"
    "The selected distance metrics are written on a single line as,\n"
    "ascii floating point values, separated by white spaces characters\n"
    "and in the order: Hausdorff; mean; median and minimum distance.\n"
    "The input objects are read from stdin and the distances are written\n"
    "to stdout unless the filenames are given.\n",
    *argv,
    " -H -N -o out.num in0.wlz in1.wlz\n"
    "Two objects are read from in0.wlz and in1.wlz, the Hausdorff and\n"
    "median of nearest neighbours distance metrics are computed and\n"
    "written to the file out.num on a single line in that order.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
