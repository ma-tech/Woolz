#ifndef DOXYGEN_SHOULD_SKIP_THIS
#pragma ident "MRC HGU $Id$"
/*!
* \file         binWlz/WlzAffineTransformVertices.c
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
		disp = 0,
		ok = 1,
		usage = 0,
		nV,
		nVC;
  double	**vA0 = NULL,
  		**vA1 = NULL,
		**vA2 = NULL;
  WlzVertex	vtx;
  WlzAffineTransform *trans = NULL;
  WlzObject	*inObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  char 		*outFileStr,
  		*inFileStr,
  		*inObjFileStr;
  static char	optList[] = "dho:t:",
		fileStrDef[] = "-";

  opterr = 0;
  inFileStr = inObjFileStr = outFileStr = fileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'd':
        disp = 1;
	break;
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
       (AlcDouble2ReadAsci(fP, &vA0, (size_t *) &nV, (size_t *) &nVC) != ALC_ER_NONE) ||
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
  /* Create array for transformed vertices. */
  if(ok)
  {
    if(AlcDouble2Malloc(&vA1, nV, nVC) != ALC_ER_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: Failed to allocate memory\n",
		     *argv);
    }
  }
  /* Transform each vertex in turn using the same array. */
  if(ok)
  {
    if(nVC == 2) /* 2D vertices */
    {
      for(idx = 0; idx < nV; ++idx)
      {
        vtx.d2.vtX = vA0[idx][0];
        vtx.d2.vtY = vA0[idx][1];
	vtx.d2 = WlzAffineTransformVertexD2(trans, vtx.d2, NULL);
	vA1[idx][0] = vtx.d2.vtX;
	vA1[idx][1] = vtx.d2.vtY;
      }
    }
    else /* 3D vertices */
    {
      for(idx = 0; idx < nV; ++idx)
      {
        vtx.d3.vtX = vA0[idx][0];
        vtx.d3.vtY = vA0[idx][1];
        vtx.d3.vtZ = vA0[idx][2];
	vtx.d3 = WlzAffineTransformVertexD3(trans, vtx.d3, NULL);
	vA1[idx][0] = vtx.d3.vtX;
	vA1[idx][1] = vtx.d3.vtY;
	vA1[idx][2] = vtx.d3.vtZ;
      }
    }
  }
  /* Write out the array. */
  if(ok)
  {
    if(disp == 0)
    {
      if((outFileStr == NULL) ||
	 (*outFileStr == '\0') ||
	 ((fP = (strcmp(outFileStr, "-")?
		fopen(outFileStr, "w"): stdout)) == NULL) ||
	 (AlcDouble2WriteAsci(fP, vA1, nV, nVC) != ALC_ER_NONE))
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: Failed to write vertices to file %s\n",
		       *argv, outFileStr);
      }
    }
    else
    {
      if(AlcDouble2Malloc(&vA2, nV, nVC * 2) != ALC_ER_NONE)
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: Failed to allocate memory\n",
		       *argv);
      }
      if(ok)
      {
	if(nVC == 2) /* 2D vertices */
	{
	  for(idx = 0; idx < nV; ++idx)
	  {
	    vA2[idx][0] = vA0[idx][0];
	    vA2[idx][1] = vA0[idx][1];
	    vA2[idx][2] = vA1[idx][0] - vA0[idx][0];
	    vA2[idx][3] = vA1[idx][1] - vA0[idx][1];
	  }
	}
	else /* 3D vertices */
	{
	  for(idx = 0; idx < nV; ++idx)
	  {
	    vA2[idx][0] = vA0[idx][0];
	    vA2[idx][1] = vA0[idx][1];
	    vA2[idx][2] = vA0[idx][2];
	    vA2[idx][3] = vA1[idx][0] - vA0[idx][0];
	    vA2[idx][4] = vA1[idx][1] - vA0[idx][1];
	    vA2[idx][5] = vA1[idx][2] - vA0[idx][2];
	  }
	}
	if((outFileStr == NULL) ||
	   (*outFileStr == '\0') ||
	   ((fP = (strcmp(outFileStr, "-")?
	          fopen(outFileStr, "w"): stdout)) == NULL) ||
           (AlcDouble2WriteAsci(fP, vA2, nV, 2 * nVC) != ALC_ER_NONE))
        {
	  ok = 0;
	  (void )fprintf(stderr,
	       "%s: Failed to write vertices and displacements to file %s\n",
			 *argv, outFileStr);
	}
      }
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      fclose(fP);
    }
    
  }
  (void )Alc2Free((void **)vA0);
  (void )Alc2Free((void **)vA1);
  (void )Alc2Free((void **)vA2);
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
    " [-d] [-h] [-o<output file>] [-t <transform>] [<input file>]\n" 
    "Options:\n"
    "  -d Output the vertices in same format used by WlzAffineTransformLSq\n"
    "       <vtx x> <vtx y> <disp x> <disp y> for 2D\n"
    "      or\n"
    "       <vtx x> <vtx y> <vtx z> <disp x> <disp y> <disp z> for 3D.\n"
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
    "Vertices are read from the file orig.num, transformed by the affine\n"
    "transfrom in trans.wlz and then written to shifted.num.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
