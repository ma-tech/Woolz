#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzAffineTransformVertices.c
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
* \brief	Reads vertices, applies an affine transform to them
*		and then writes transformed vertices.
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
  		option,
		ok = 1,
		usage = 0,
		nV,
		nVC;
  double	**vAry = NULL;
  WlzVertex	vtx;
  WlzAffineTransform *trans = NULL;
  WlzObject	*inObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  char 		*outFileStr,
  		*inFileStr,
  		*inObjFileStr;
  static char	optList[] = "ho:t:",
		fileStrDef[] = "-";

  opterr = 0;
  inFileStr = inObjFileStr = outFileStr = fileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outFileStr = optarg;
	break;
      case 't':
        inObjFileStr = optarg;
	break;
      case 'h':
      default:
        usage = 1;
	ok = 0;
	break;
    }
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
      inFileStr = *(argv + optind);
    }
  }
  /* Read the transform. */
  if(ok)
  {
    if((inObjFileStr == NULL) ||
       (*inObjFileStr == '\0') ||
       ((fP = (strcmp(inObjFileStr, "-")?
	      fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to read transform from file %s\n",
		     *argv, inObjFileStr);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }
    if(inObj)
    {
      if((inObj->type != WLZ_AFFINE_TRANS) ||
	 (inObj->domain.core == NULL))
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: Invalid transform object read from file %s\n",
		       *argv, inObjFileStr);
      }
      else
      {
	trans = WlzAssignAffineTransform(inObj->domain.t, NULL);
      }
    }
  }
  /* Read the vertices into an array. */
  if(ok)
  {
    if((inFileStr == NULL) ||
       (*inFileStr == '\0') ||
       ((fP = (strcmp(inFileStr, "-")?
	      fopen(inFileStr, "r"): stdin)) == NULL) ||
       (AlcDouble2ReadAsci(fP, &vAry, &nV, &nVC) != ALC_ER_NONE) ||
       (nV < 1) || ((nVC != 2) && (nVC != 3)))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to read vertices from file %s\n",
		     *argv, inFileStr);
    }
    if(fP && strcmp(inFileStr, "-"))
    {
      fclose(fP);
    }
  }
  /* Transform each vertex in turn using the same array. */
  if(ok)
  {
    if(nVC == 2) /* 2D vertices */
    {
      for(idx = 0; idx < nV; ++idx)
      {
        vtx.d2.vtX = *(*(vAry + idx) + 0);
        vtx.d2.vtY = *(*(vAry + idx) + 1);
	vtx.d2 = WlzAffineTransformVertexD2(trans, vtx.d2, NULL);
	*(*(vAry + idx) + 0) = vtx.d2.vtX;
	*(*(vAry + idx) + 1) = vtx.d2.vtY;
      }
    }
    else /* 3D vertices */
    {
      for(idx = 0; idx < nV; ++idx)
      {
        vtx.d3.vtX = *(*(vAry + idx) + 0);
        vtx.d3.vtY = *(*(vAry + idx) + 1);
        vtx.d3.vtZ = *(*(vAry + idx) + 2);
	vtx.d3 = WlzAffineTransformVertexD3(trans, vtx.d3, NULL);
	*(*(vAry + idx) + 0) = vtx.d3.vtX;
	*(*(vAry + idx) + 1) = vtx.d3.vtY;
	*(*(vAry + idx) + 2) = vtx.d3.vtZ;
      }
    }
  }
  /* Write out the array. */
  if(ok)
  {
    if((outFileStr == NULL) ||
       (*outFileStr == '\0') ||
       ((fP = (strcmp(outFileStr, "-")?
	      fopen(outFileStr, "w"): stdout)) == NULL) ||
       (AlcDouble2WriteAsci(fP, vAry, nV, nVC) != ALC_ER_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to write vertices to file %s\n",
		     *argv, outFileStr);
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      fclose(fP);
    }
    
  }
  if(vAry)
  {
    (void )Alc2Free((void **)vAry);
  }
  if(inObj)
  {
    WlzFreeObj(inObj);
  }
  if(trans)
  {
    WlzFreeAffineTransform(trans);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-h] [-o<output file>] [-t <transform>] [<input file>]\n" 
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -o  Output object file name.\n"
    "  -t  Input affine transform object.\n" 
    "Reads vertices, applies an affine transform to them and then writes\n"
    "out the transformed vertices.\n"
    "The input vertices are read from stdin and the transformed vertices\n"
    "are written to stdout unless the filenames are given.\n",
    "The vertex format for both the input and output is: Space (or tab)\n"
    "separated ascii floating point, with either two or three floating point\n"
    "numbers per vertex and one vertex per line in x, y, z order. All\n"
    "vertices must be of the same type (either 2D or 3D).\n",
    *argv,
    " -t trans.wlz -o shifted.num orig.num\n"
    "Vertices are read from the file orig.num, transformed by the affine"
    "transfrom in trans.wlz and then written to shifted.num.\n");
  }
  return(!ok);
}
