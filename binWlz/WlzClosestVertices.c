#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzClosestVertices.c
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
* \brief	Finds the closest vertices in the given object to those in the
*		given list of vertices.
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <float.h>
#include <Wlz.h>

static int			WlzClosestTransformVtx(
				  char *tFileStr,
				  int dim,
				  int nV,
				  WlzVertexP vP);
static int			WlzClosestEnsureDVtx(
				  int nV,
				  WlzVertexP *vP,
				  WlzVertexType *vType,
				  int dim);
static void			WlzVerticesClosestKDTree(
				  AlcKDTTree *tree,
				  int nTV,
				  WlzVertexP tV,
				  int nSV,
				  WlzVertexP sV,
				  int dim);
static void			WlzVerticesClosestDumb(
				  int nTV,
				  WlzVertexP tV,
				  int nSV,
				  WlzVertexP sV,
				  int dim);

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		idN,
		nIOC,
		nIOR,
		nTV,
		nSV,
  		option,
		dim = 2,
		ok = 1,
		usage = 0;
  size_t	sz;
  WlzVertexType	tVType;
  WlzVertexP	sV,
  		tV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  int		*shfBuf = NULL;
  double	**ioA = NULL;
  AlcKDTTree 	*tree = NULL;
  WlzObject	*obj = NULL;
  FILE		*fP = NULL;
  char 		*inFileStr,
  		*inObjFileStr,
		*outFileStr,
		*sTrFileStr,
		*tTrFileStr;
  static char	optList[] = "h23s:t:u:o:";
  const int	thrForKDTree = 50; /* Use kD-tree if more than 50 source
  				    * vertices. */

  opterr = 0;
  sV.v = tV.v = NULL;
  sTrFileStr = tTrFileStr = NULL;
  inFileStr = outFileStr = inObjFileStr = "-";
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case '2':
        dim = 2;
	break;
      case '3':
        dim = 3;
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 's':
        sTrFileStr = optarg;
	break;
      case 't':
        tTrFileStr = optarg;
	break;
      case 'u':
        inFileStr = optarg;
	break;
      case 'h':
      default: /* FALLTHROUGH */
        usage = 1;
	break;
    }
  }
  if(ok && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
    }
    else
    {
      inObjFileStr = *(argv + optind);
    }
  }
  ok = !usage;
  /* Read the target object. */
  if(ok)
  {
    if((inObjFileStr == NULL) ||
       (*inObjFileStr == '\0') ||
       ((fP = (strcmp(inObjFileStr, "-")?
	      fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((obj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to read object from file %s.\n",
		     *argv, inObjFileStr);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  /* Read the source vertices. */
  if(ok)
  {
    if((inFileStr == NULL) ||
       (*inFileStr == '\0') ||
       ((fP = (strcmp(inFileStr, "-")?
	      fopen(inFileStr, "r"): stdin)) == NULL) ||
       (AlcDouble2ReadAsci(fP, &ioA, &nIOR, &nIOC) != ALC_ER_NONE) ||
       (nIOR < 1) || ((nIOC != 2) && (nIOC != 3)) ||
       ((dim == 3) && (nIOC == 2)))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to read vertices from file %s.\n",
		     *argv, inFileStr);
    }
    if(fP && strcmp(inFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok)
  {
    nSV = nIOR;
    sz = (dim == 2)? sizeof(WlzDVertex2): sizeof(WlzDVertex3);
    if((sV.v = AlcMalloc(nSV * sz)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: Failed to allocate memory.\n",
		     *argv);
    }
  }
  if(ok)
  {
    if(dim == 2)
    {
      for(idN = 0; idN < nSV; ++idN)
      {
        (sV.d2 + idN)->vtX = ioA[idN][0];
        (sV.d2 + idN)->vtY = ioA[idN][1];
      }
    }
    else
    {
      for(idN = 0; idN < nSV; ++idN)
      {
        (sV.d3 + idN)->vtX = ioA[idN][0];
        (sV.d3 + idN)->vtY = ioA[idN][1];
        (sV.d3 + idN)->vtZ = ioA[idN][2];
      }
    }
  }
  /* Get target vertices from the target object. */
  if(ok)
  {
    tV = WlzVerticesFromObj(obj, NULL, &nTV, &tVType, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to get vertices from object.\n",
		     *argv);
    }
  }
  (void )WlzFreeObj(obj); obj = NULL;
  if(ok)
  {
    ok = WlzClosestEnsureDVtx(nTV, &tV, &tVType, dim);
    if(ok == 0)
    {
      (void )fprintf(stderr,
		     "%s: Failed to convert vertices to double precision\n"
		     "file %s.\n",
		     *argv, sTrFileStr);
    }
  }
  /* Transform source vertices. */
  if(ok && sTrFileStr)
  {
    ok = WlzClosestTransformVtx(sTrFileStr, dim, nSV, sV);
    if(ok == 0)
    {
      (void )fprintf(stderr,
		     "%s: Failed to transform vertices using transform from\n"
		     "file %s.\n",
		     *argv, sTrFileStr);
    }
  }
  /* Transform target vertices. */
  if(ok && tTrFileStr)
  {
    ok = WlzClosestTransformVtx(tTrFileStr, dim, nTV, tV);
    if(ok == 0)
    {
      (void )fprintf(stderr,
		     "%s: Failed to transform vertices using transform from\n"
		     "file %s.\n",
		     *argv, tTrFileStr);
    }
  }
  /* Find the closest vertices. */
  if(ok)
  {
    if(nSV > thrForKDTree)
    {
      if(((shfBuf = (int *)AlcMalloc(nTV * sizeof(int))) == NULL) ||
         ((tree = WlzVerticesBuildTree(tVType, nTV, tV, shfBuf,
	 			       &errNum)) == NULL))
      {
        ok = 0;
	(void )fprintf(stderr,
		       "%s: Failed to compute kD-tree.\n",
		       *argv);
      }
      else
      {
        WlzVerticesClosestKDTree(tree, nTV, tV, nSV, sV, dim);
      }
    }
    else
    {
      WlzVerticesClosestDumb(nTV, tV, nSV, sV, dim);
    }
    (void )AlcFree(shfBuf); shfBuf = NULL;
    (void )AlcKDTTreeFree(tree); tree = NULL;
  }
  /* Output the closest target vertices. */
  if(ok)
  {
    if(dim == 2)
    {
      for(idN = 0; idN < nSV; ++idN)
      {
        ioA[idN][0] = (sV.d2 + idN)->vtX;
        ioA[idN][1] = (sV.d2 + idN)->vtY;
      }
    }
    else
    {
      for(idN = 0; idN < nSV; ++idN)
      {
        ioA[idN][0] = (sV.d3 + idN)->vtX;
        ioA[idN][1] = (sV.d3 + idN)->vtY;
        ioA[idN][2] = (sV.d3 + idN)->vtZ;
      }
    }
    nIOR = nSV;
    nIOC = dim;
    if((outFileStr == NULL) ||
       (*outFileStr == '\0') ||
       ((fP = (strcmp(outFileStr, "-")?
	      fopen(outFileStr, "w"): stdout)) == NULL) ||
       (AlcDouble2WriteAsci(fP, ioA, nIOR, nIOC) != ALC_ER_NONE))
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
  (void )Alc2Free((void **)ioA);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-h] [-2] [-3] [-s <transform>] [-t <transform>]\n"
    "        [-u <vertices>] [-o<output file>] [<input file>]\n" 
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -2  Vertices and object are 2D.\n"
    "  -3  Vertices and object are 3D.\n"
    "  -s  Affine transform to be applied to the given vertices.\n"
    "  -t  Affine transform o be applied to the given object.\n"
    "  -u  File containing vertices encoded as ASCII in x, y, z order.\n"
    "  -o  Output file for closest vertices to those given that have\n"
    "      been found in in the given object. These are in the same order\n"
    "      as the given vertices and have the same format.\n"
    "Reads vertices, optional affine transforms and a Woolz object which\n"
    "has a vertex representation (such as boundaries and contours) and\n"
    "finds the closest points ion the given object to those given.\n",
    *argv,
    " -3 -s tr.wlz -o out.num -u in.num surf.wlz\n"
    "Reads a target surface object from file surf.wlz and 3D vertices\n"
    "from file in.num. Transforms the vertices using an affine transform\n"
    "read from the file tr.wlz. For each transformed vertex finds the\n"
    "closest vertex in the surface object and then saves the list of\n"
    "closest surface vertices to the file out.num.\n");
  }
  return(!ok);
}

/*!
* \return	Zero on error.
* \brief	Ensures that the given vertices are double precision
*		by reallocating and converting if necessary.
* \param	nV			Number of vertices.
* \param	vP			Pointer to vertices.
* \param	vType			Type of vertices.
* \param	dim			Dimension.
*/
static int	WlzClosestEnsureDVtx(int nV, WlzVertexP *vP,
				     WlzVertexType *vType, int dim)
{
  int		ok = 1;
  WlzVertexP	tV;

  switch(*vType)
  {
    case WLZ_VERTEX_I2: /* FALLTHROUGH */
    case WLZ_VERTEX_F2:
      if(dim == 2)
      {
        if((tV.v = AlcMalloc(nV * sizeof(WlzDVertex2))) == NULL)
	{
	  ok = 0;
	}
        else
	{
	  if(*vType == WLZ_VERTEX_I2)
	  {
	    WlzValueCopyIVertexToDVertex(tV.d2, (*vP).i2, nV);
	  }
	  else /* *vType == WLZ_VERTEX_F2 */
	  {
	    WlzValueCopyFVertexToDVertex(tV.d2, (*vP).f2, nV);
	  }
	  AlcFree((*vP).v);
	  *vP = tV;
	  *vType = WLZ_VERTEX_D2;
	}
      }
      else
      {
        ok = 0;
      }
      break;
    case WLZ_VERTEX_D2:
      ok = dim == 2;
      break;
    case WLZ_VERTEX_I3: /* FALLTHROUGH */
    case WLZ_VERTEX_F3:
      if(dim == 3)
      {
        if((tV.v = AlcMalloc(nV * sizeof(WlzDVertex3))) == NULL)
	{
	  ok = 0;
	}
        else
	{
	  if(*vType == WLZ_VERTEX_I3)
	  {
	    WlzValueCopyIVertexToDVertex3(tV.d3, (*vP).i3, nV);
	  }
	  else /* *vType == WLZ_VERTEX_F3 */
	  {
	    WlzValueCopyFVertexToDVertex3(tV.d3, (*vP).f3, nV);
	  }
	  AlcFree((*vP).v);
	  *vP = tV;
	  *vType = WLZ_VERTEX_D3;
	}
      }
      else
      {
        ok = 0;
      }
      break;
    case WLZ_VERTEX_D3:
      ok = dim == 3;
      break;
    default:
      ok = 0;
      break;
  }
  return(ok);
}

/*!
* \return	Non zero value if no error.
* \brief	Reads an affine transform object from the file
*		with the given name and then transforms the given
*		vertices.
* \param	tFileStr		File name, "-" is stain.
* \param	dim			Dimension must be 2 or 3.
* \param	nV			Number of vertices.
* \param	vP			Vertex pointer with either 2D or
*					3D vertices accoring to dim.
*/
static int	WlzClosestTransformVtx(char *tFileStr, int dim,
				       int nV, WlzVertexP vP)
{
  int 		idN,
  		ok = 1;
  FILE		*fP = NULL;
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((*tFileStr == '\0') ||
     ((fP = (strcmp(tFileStr, "-")?
	    fopen(tFileStr, "r"): stdin)) == NULL) ||
     ((obj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
     (obj->type != WLZ_AFFINE_TRANS) ||
     ((dim == 2) && (obj->domain.t->type != WLZ_TRANSFORM_2D_AFFINE)) ||
     ((dim == 3) && (obj->domain.t->type != WLZ_TRANSFORM_3D_AFFINE)))
  {
    ok = 0;
  }
  if(fP && strcmp(tFileStr, "-"))
  {
    (void )fclose(fP);
  }
  if(ok)
  {
    if(dim == 2)
    {
      for(idN = 0; idN < nV; ++idN)
      {
	*(vP.d2 + idN) = WlzAffineTransformVertexD2(obj->domain.t,
						    *(vP.d2 + idN), NULL);
      }
    }
    else /* dim == 3 */
    {
      for(idN = 0; idN < nV; ++idN)
      {
	*(vP.d3 + idN) = WlzAffineTransformVertexD3(obj->domain.t,
						    *(vP.d3 + idN), NULL);
      }
    }
  }
  (void )WlzFreeObj(obj);
  return(ok);
}

/*!
* \return	void
* \brief	For each source vertex finds the closest target vertex using
*		a dumb search which is faster for small numbers of source
*		vertices.
* \param	nTV			Number of target vertices.
* \param	tV			Target vertices.
* \param	nSV			Number of source vertices.
* \param	sV			Source vertices on entry and closest
*					target vertices on return.
* \param	dim			Dimension == 2 || 3.
*/
static void	WlzVerticesClosestDumb(int nTV, WlzVertexP tV,
				       int nSV, WlzVertexP sV, int dim)
{
  int		idS,
  		idT;
  double	cDist,
  		mDist;
  WlzVertex	delta,
  		mVtx,
		sVtx,
  		tVtx;
  		

  if(dim == 2)
  {
    for(idS = 0; idS < nSV; ++idS)
    {
      sVtx.d2 = *(sV.d2 + idS);
      mVtx.d2 = *(tV.d2 + 0);
      WLZ_VTX_2_SUB(delta.d2, sVtx.d2, mVtx.d2);
      mDist = WLZ_VTX_2_SQRLEN(delta.d2);
      for(idT = 1; idT < nTV; ++idT)
      {
	tVtx.d2 = *(tV.d2 + idT);
	WLZ_VTX_2_SUB(delta.d2, sVtx.d2, tVtx.d2);
	cDist = WLZ_VTX_2_SQRLEN(delta.d2);
	if(cDist < mDist)
	{
	  mDist = cDist;
	  mVtx.d2 = tVtx.d2;
	}
      }
      *(sV.d2 + idS) = mVtx.d2;
    }
  }
  else /* dim == 3 */
  {
    for(idS = 0; idS < nSV; ++idS)
    {
      sVtx.d3 = *(sV.d3 + idS);
      mVtx.d3 = *(tV.d3 + 0);
      WLZ_VTX_3_SUB(delta.d3, sVtx.d3, mVtx.d3);
      mDist = WLZ_VTX_3_SQRLEN(delta.d3);
      for(idT = 1; idT < nTV; ++idT)
      {
	tVtx.d3 = *(tV.d3 + idT);
	WLZ_VTX_3_SUB(delta.d3, sVtx.d3, tVtx.d3);
	cDist = WLZ_VTX_3_SQRLEN(delta.d3);
	if(cDist < mDist)
	{
	  mDist = cDist;
	  mVtx.d3 = tVtx.d3;
	}
      }
      *(sV.d3 + idS) = mVtx.d3;
    }
  }
}

/*!
* \return	void
* \brief	For each of the source vertices finds the closest target
*		vertex in the kD-tree.
* \param	tree			kD-tree of target vertices.
* \param        nTV                     Number of target vertices.
* \param        tV                      Target vertices.
* \param	nSV			Number of source vertices.
* \param	sV			Source vertices on entry and closest
*					target vertices on return.
* \param	dim			Dimension == 2 || 3.
*/
static void	WlzVerticesClosestKDTree(AlcKDTTree *tree,
					 int nTV, WlzVertexP tV,
					 int nSV, WlzVertexP sV, int dim)
{
  int		idS;
  AlcKDTNode	*node;
  double	pos[3];

  if(dim == 2)
  {
    for(idS = 0; idS < nSV; ++idS)
    {
      pos[0] = (sV.d2 + idS)->vtX;
      pos[1] = (sV.d2 + idS)->vtY;
      node = AlcKDTGetNN(tree, pos, DBL_MAX, NULL, NULL);
      *(sV.d2 + idS) = *(tV.d2 + node->idx);
    }
  }
  else /* dim == 3 */
  {
    for(idS = 0; idS < nSV; ++idS)
    {
      pos[0] = (sV.d3 + idS)->vtX;
      pos[1] = (sV.d3 + idS)->vtY;
      pos[2] = (sV.d3 + idS)->vtZ;
      node = AlcKDTGetNN(tree, pos, DBL_MAX, NULL, NULL);
      *(sV.d3 + idS) = *(tV.d3 + node->idx);
    }
  }
}
