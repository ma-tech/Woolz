#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzSnapFitObjs.c
* \author       Bill Hill
* \date         July 2004
* \version      $Id$
* \note
*               Copyright
*               2004 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Given an pair of objects which are in good alignment, finds
*		a set of correspondences between the two objects based only
*		on closest points.
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <float.h>
#include <Wlz.h>


extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		idN,
  		nVtx,
  		absTPMode = 0,
  		option,
		ok = 1,
		usage = 0;
  double	maxCDist = DBL_MAX,
  		minTDist = 10.0,
		minSDist = 10.0;
  WlzVertexType vType;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzVertexP	tVtxP,
  		sVtxP;
  WlzObject	*tObj = NULL,
  		*sObj = NULL,
		*trObj = NULL;
  FILE		*fP = NULL;
  char 		*tObjFileStr,
  		*sObjFileStr,
		*trObjFileStr,
		*outFileStr;
  static char	optList[] = "hAi:d:s:t:o:";

  opterr = 0;
  trObjFileStr = NULL;
  tVtxP.v = sVtxP.v = NULL;
  trObjFileStr = NULL;
  sObjFileStr = tObjFileStr = outFileStr = "-";
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'A':
        absTPMode = 1; 
	break;
      case 'i':
        trObjFileStr = optarg;
	break;
      case 'd':
        if(sscanf(optarg, "%lg", &maxCDist) != 1)
	{
	  usage = 1;
	}
	break;
      case 's':
        if(sscanf(optarg, "%lg", &minSDist) != 1)
	{
	  usage = 1;
	}
	break;
      case 't':
        if(sscanf(optarg, "%lg", &minTDist) != 1)
	{
	  usage = 1;
	}
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'h':
      default: /* FALLTHROUGH */
        usage = 1;
	break;
    }
  }
  if(ok && (optind < argc))
  {
    if((optind + 2) != argc)
    {
      usage = 1;
    }
    else
    {
      tObjFileStr = *(argv + optind);
      sObjFileStr = *(argv + optind + 1);
    }
  }
  ok = !usage;
  /* Read the target object. */
  if(ok)
  {
    if((tObjFileStr == NULL) ||
       (*tObjFileStr == '\0') ||
       ((fP = (strcmp(tObjFileStr, "-")?
	      fopen(tObjFileStr, "r"): stdin)) == NULL) ||
       ((tObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to read target object from file %s.\n",
		     *argv, tObjFileStr);
    }
    if(fP && strcmp(tObjFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  /* Read the source object. */
  if(ok)
  {
    if((sObjFileStr == NULL) ||
       (*sObjFileStr == '\0') ||
       ((fP = (strcmp(sObjFileStr, "-")?
	      fopen(sObjFileStr, "r"): stdin)) == NULL) ||
       ((sObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to read source object from file %s.\n",
		     *argv, sObjFileStr);
    }
    if(fP && strcmp(sObjFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  /* Read the optional initial transform object. */
  if(ok && trObjFileStr)
  {
    if((*trObjFileStr == '\0') ||
       ((fP = (strcmp(trObjFileStr, "-")?
	      fopen(trObjFileStr, "r"): stdin)) == NULL) ||
       ((trObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (trObj == NULL) || (trObj->type != WLZ_AFFINE_TRANS) ||
       (trObj->domain.core == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
		 "%s: Failed to read affine transform object from file %s.\n",
		     *argv, trObjFileStr);
    }
    if(fP && strcmp(trObjFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  /* Compute correspondences. */
  if(ok)
  {
    errNum = WlzSnapFit(tObj, sObj, (trObj)? trObj->domain.t: NULL,
    			&vType, &nVtx, &tVtxP, &sVtxP,
			maxCDist, minTDist, minSDist);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to compute correspondences.\n",
		     *argv);
    }
  }
  /* Output the correspondences. */
  if(ok)
  {
    if((outFileStr == NULL) ||
       (*outFileStr == '\0') ||
       ((fP = (strcmp(outFileStr, "-")?
	      fopen(outFileStr, "w"): stdout)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to write correspondences to file %s\n",
		     *argv, outFileStr);
    }
  }
  if(ok)
  {
    if(vType == WLZ_VERTEX_D2)
    {
      if(absTPMode)
      {
	for(idN = 0; idN < nVtx; ++idN)
	{
	  (void )fprintf(fP, "%g %g %g %g\n", 
	  		 (tVtxP.d2 + idN)->vtX,
			 (tVtxP.d2 + idN)->vtY,
	  		 (sVtxP.d2 + idN)->vtX,
			 (sVtxP.d2 + idN)->vtY);
	}
      }
      else
      {
	for(idN = 0; idN < nVtx; ++idN)
	{
	  (void )fprintf(fP, "%g %g %g %g\n", 
	  		 (sVtxP.d2 + idN)->vtX,
			 (sVtxP.d2 + idN)->vtY,
	  		 (tVtxP.d2 + idN)->vtX - (sVtxP.d2 + idN)->vtX,
	  		 (tVtxP.d2 + idN)->vtY - (sVtxP.d2 + idN)->vtY);
	}
      }
    }
    else /* vType == WLZ_VERTEX_D3 */
    {
      if(absTPMode)
      {
	for(idN = 0; idN < nVtx; ++idN)
	{
	  (void )fprintf(fP, "%g %g %g %g %g %g\n", 
	  		 (tVtxP.d3 + idN)->vtX,
			 (tVtxP.d3 + idN)->vtY,
			 (tVtxP.d3 + idN)->vtZ,
	  		 (sVtxP.d3 + idN)->vtX,
			 (sVtxP.d3 + idN)->vtY,
			 (sVtxP.d3 + idN)->vtZ);
	}
      }
      else
      {
	for(idN = 0; idN < nVtx; ++idN)
	{
	  (void )fprintf(fP, "%g %g %g %g %g %g\n", 
	  		 (sVtxP.d3 + idN)->vtX,
			 (sVtxP.d3 + idN)->vtY,
			 (sVtxP.d3 + idN)->vtZ,
	  		 (tVtxP.d3 + idN)->vtX - (sVtxP.d3 + idN)->vtX,
	  		 (tVtxP.d3 + idN)->vtY - (sVtxP.d3 + idN)->vtY,
	  		 (tVtxP.d3 + idN)->vtZ - (sVtxP.d3 + idN)->vtZ);
	}
      }
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      fclose(fP); fP = NULL;
    }
    
  }
  (void )WlzFreeObj(tObj);
  (void )WlzFreeObj(sObj);
  (void )WlzFreeObj(trObj);
  (void )AlcFree(tVtxP.v);
  (void )AlcFree(sVtxP.v);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-h] [-A] [-i <transform>] [-d #] [-t #] [-s #]\n"
    "        [-o<output file>] <target> <source>\n" 
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -A  Correspondences output in absolute format.\n"
    "  -i  Initial affine transform to be applied to the source object.\n"
    "  -d  Maximum distance between any target source correspondences.\n"
    "  -t  Minimum distance between target correspondence points.\n"
    "  -s  Minimum distance between source correspondence points.\n"
    "  -o  Output file for correspondences.\n"
    "Computes a set of correspondences (tie points) from the given\n"
    "source object to the given target object based on closest points.\n"
    "The correspondences are output as ascii vertices either in relative\n"
    "format:\n"
    "  <sx> <sy>[ <sz> ]<tx - sx> <ty - sy>[ <tz - sz>]\n"
    "or absolute format:\n"
    "  <tx> <ty>[ <tz> ]<sx> <sy>[ <sz>],\n"
    "where sx, sy, sz, tx, ty and tz are the source and target coordinates.\n",
    *argv,
    " -3 -d 10 -i tr.wlz -o out.num -t 20 -s 20 trg.wlz src.wlz\n"
    "Reads an initial source affine transform from tr.wlz and computes\n"
    "a set of corresponding closest points such that the minimum distance\n"
    "between and target and source pair is less than 10, no target points\n"
    "are closer than 20 and no target points have corresponding points\n"
    "closer than 30 in the source object. all distances are in pixels/\n"
    "voxels.\n");
  }
  return(!ok);
}
