#pragma ident "MRC HGU $Id$"
/************************************************************************
* Project:      Woolz
* Title:        WlzBasisFnTransformObj.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Woolz filter which applies a basis function to a given
*		object.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <Wlz.h>

#define	IN_RECORD_MAX   (1024)

static void	WlzMeshOutputPS(WlzMeshTransform *);

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		nTiePP,
		option,
		vxCount = 0,
		vxLimit = 0,
		basisFnPolyOrder = 3,
		debugPSOutput = 0,
		outBasisTrFlag = 0,
		outMeshTrFlag = 0,
		restrictToBasisFn = 0,
		ok = 1,
		usage = 0;
  double	meshMinDist = 20.0,
  		meshMaxDist = 40.0;
  WlzDVertex2	*vx0,
  		*vx1,
		*vxVec0 = NULL,
		*vxVec1 = NULL;
  WlzObject	*inObj = NULL,
		*outObj = NULL;
  WlzMeshTransform *meshTr = NULL;
  WlzMeshGenMethod meshGenMth = WLZ_MESH_GENMETHOD_GRADIENT;
  WlzBasisFnTransform *basisTr = NULL;
  WlzBasisFnType basisFnType = WLZ_BASISFN_MQ;
  WlzInterpolationType interp = WLZ_INTERPOLATION_NEAREST;
  FILE		*fP = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char 		*rec,
  		*inObjFileStr,
		*basisFnTrFileStr = NULL,
		*tiePtFileStr = NULL,
  		*outObjFileStr;
  const char    *errMsg;
  static char	optList[] = "m:o:p:t:M:Y:cghqsyBDGLRTU",
  		inObjFileStrDef[] = "-",
		outObjFileStrDef[] = "-",
  		inRecord[IN_RECORD_MAX];

  opterr = 0;
  outObjFileStr = outObjFileStrDef;
  inObjFileStr = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'B':
        meshGenMth = WLZ_MESH_GENMETHOD_BLOCK;
	break;
      case 'D':
        debugPSOutput = 1;
	break;
      case 'c':
        basisFnType = WLZ_BASISFN_CONF_POLY;
	break;
      case 'g':
        basisFnType = WLZ_BASISFN_GAUSS;
	break;
      case 'G':
        meshGenMth = WLZ_MESH_GENMETHOD_GRADIENT;
	break;
      case 'L':
        interp = WLZ_INTERPOLATION_LINEAR;
	break;
      case 'm':
	if(sscanf(optarg, "%lg", &meshMinDist) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'M':
	if(sscanf(optarg, "%lg", &meshMaxDist) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'p':
        tiePtFileStr = optarg;
	break;
      case 'q':
        basisFnType = WLZ_BASISFN_MQ;
	break;
      case 's':
        basisFnType = WLZ_BASISFN_TPS;
	break;
      case 'y':
        basisFnType = WLZ_BASISFN_POLY;
	break;
      case 't':
        basisFnTrFileStr = optarg;
        break;
      case 'R':
        restrictToBasisFn = 1;
	break;
      case 'T':
        outBasisTrFlag = 1;
	break;
      case 'U':
        outMeshTrFlag = 1;
	break;
      case 'Y':
        if(sscanf(optarg, "%d", &basisFnPolyOrder) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'h':
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
     (((tiePtFileStr == NULL) || (*tiePtFileStr == '\0')) &&
      ((basisFnTrFileStr == NULL) || (*basisFnTrFileStr == '\0'))) ||
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
    if(basisFnTrFileStr)
    {
      errNum = WLZ_ERR_UNSPECIFIED;
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: Reading basis function transforms has not been\n"
		     "implemented yet\n",
		     *argv);
    }
    else if(tiePtFileStr)
    {
      if((fP = (strcmp(tiePtFileStr, "-")?
                fopen(tiePtFileStr, "r"): stdin)) == NULL)
      {
        ok = 0;
	errNum = WLZ_ERR_READ_EOF;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to open tie points file %s (%s).\n",
		       *argv, tiePtFileStr, errMsg);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	while((errNum == WLZ_ERR_NONE) &&
	      (fgets(inRecord, IN_RECORD_MAX - 1, fP) != NULL))
	{
	  inRecord[IN_RECORD_MAX - 1] = '\0';
	  rec = inRecord;
	  while(*rec && isspace(*rec))
	  {
	    ++rec;
	  }
	  if(*rec && (*rec != '#'))
	  {
	    if(vxCount >= vxLimit)
	    {
	      vxLimit = (vxLimit + 1024) * 2;
	      if(((vxVec0 = (WlzDVertex2 *)AlcRealloc(vxVec0,
				     vxLimit * sizeof(WlzDVertex2))) == NULL) ||
		 ((vxVec1 = (WlzDVertex2 *)AlcRealloc(vxVec1,
				     vxLimit * sizeof(WlzDVertex2))) == NULL))
	      {
		errNum = WLZ_ERR_MEM_ALLOC;
	      }
	      else
	      {
		vx0 = vxVec0 + vxCount;
		vx1 = vxVec1 + vxCount;
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      if(sscanf(rec, "%lg %lg %lg %lg", &(vx0->vtX), &(vx0->vtY),
			&(vx1->vtX), &(vx1->vtY)) != 4)
	      {
		errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      else
	      {
		vx1->vtX += vx0->vtX;
		vx1->vtY += vx0->vtY;
		++vx0;
		++vx1;
		++vxCount;
	      }
	    }
	  }
	}
	nTiePP = vxCount;
	if(errNum != WLZ_ERR_NONE)
	{
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
			 "%s: failed to read tie points file %s (%s).\n",
			 *argv, tiePtFileStr, errMsg);
	}
      }
      if(fP && strcmp(tiePtFileStr, "-"))
      {
	fclose(fP);
      }
      fP = NULL;
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_READ_EOF;
    if((inObjFileStr == NULL) ||
	(*inObjFileStr == '\0') ||
	((fP = (strcmp(inObjFileStr, "-")?
		fopen(inObjFileStr, "r"): stdin)) == NULL) ||
	((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
	(errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: failed to read object from file %s\n",
		     *argv, inObjFileStr);
    }
    if(fP)
    {
      if(strcmp(inObjFileStr, "-"))
      {
	fclose(fP);
      }
      fP = NULL;
    }
  }
  if(ok)
  {
    if(restrictToBasisFn)
    {
      basisTr = WlzBasisFnTrFromCPts(basisFnType, basisFnPolyOrder,
	  nTiePP, vxVec0,
	  nTiePP, vxVec1, &errNum);
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		     "%s: failed to compute basis function transform (%s).\n",
		     *argv, errMsg);
      }
    }
    else /* Use affine and basis functions to define mesh. */
    {
      meshTr = WlzMeshTransformFromCPts(inObj, basisFnType, basisFnPolyOrder,
      				nTiePP, vxVec0, nTiePP, vxVec1,
				meshGenMth, meshMinDist, meshMaxDist,
				&errNum);
    }
  }
  if(ok)
  {
    if(outMeshTrFlag)
    {
      errNum = WLZ_ERR_UNIMPLEMENTED;
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: Writing mesh transforms has not been\n"
		     "implemented yet\n",
		     *argv);
    }
    else if(outBasisTrFlag)
    {
      errNum = WLZ_ERR_UNIMPLEMENTED;
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: Writing basis function transforms has not been\n"
		     "implemented yet\n",
		     *argv);
    }
    else /* Transform the object and then write it out. */
    {
      if(restrictToBasisFn)
      {
	meshTr = WlzMeshFromObj(inObj, meshGenMth, meshMinDist, meshMaxDist,
	    &errNum);
	if(errNum != WLZ_ERR_NONE)
	{
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
			 "%s: failed to compute a mesh for the object (%s).\n",
			 *argv, errMsg);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzBasisFnSetMesh(meshTr, basisTr);
	  if(errNum != WLZ_ERR_NONE)
	  {
	    ok = 0;
	    (void )WlzStringFromErrorNum(errNum, &errMsg);
	    (void )fprintf(stderr,
			   "%s: failed to set the mesh displacements (%s).\n",
			   *argv, errMsg);
	  }
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(debugPSOutput)
	{
	  WlzMeshOutputPS(meshTr);
	}
	outObj = WlzMeshTransformObj(inObj, meshTr, interp, &errNum);
      }
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to transform object (%s).\n",
		       *argv, errMsg);
      }
      if(errNum == WLZ_ERR_NONE)
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
    }
  }
  (void )WlzMeshFreeTransform(meshTr);
  (void )WlzBasisFnFreeTransform(basisTr);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-o<out object>] [-p<tie points file>]\n"
    "                  [-m<min mesh dist>] [-M<max mesh dist>]\n"
    "                  [-t<basis fn transform>] [-Y<order of polynomial>]\n"
    "                  [-g] [-h] [-q] [-s] [-y] [-B] [-D] [-G] [-L] [-T]\n"
    "                  [<in object>]\n"
    "Options:\n"
    "  -B  Block mesh generation method (default).\n"
    "  -D  Output the mesh as postscript to the standard error output,\n"
    "      may be useful for debugging.\n"
    "  -G  Gradient mesh generation method.\n"
    "  -L  Use linear interpolation instead of nearest neighbour.\n"
    "  -m  Minimum mesh node separation distance (default 10.0)\n"
    "  -M  Maximum mesh node separation distance (default 100.0)\n"
    "  -R  Restrict the transformation to only a basis function transform\n"
    "      instead of the default which uses a least squares affine and\n"
    "      basis function transform.\n"
    "  -o  Output object file name.\n"
    "  -p  Tie point file.\n"
    "  -t  Basis function transform object.\n"
    "  -c  Use conformal polynomial basis function if tie points are given.\n"
    "  -g  Use Gaussian basis function if tie points are given.\n"
    "  -h  Help, prints this usage message.\n"
    "  -q  Use multi-quadric basis function if tie points are given.\n"
    "  -s  Use thin plate spline basis function (default) if tie points\n"
    "      are given.\n"
    "  -T  Output a basis function transform instead of a transformed\n"
    "      object.\n"
    "  -U  Output a mesh transform instead of a transformed object.\n"
    "  -y  Use polynomianl basis function if tie points are given.\n"
    "  -Y  Polynomial order for oolynomianl basis function (default 3).\n"
    "Computes and applies Woolz basis function transforms.\n"
    "Tie points may be read from an ascii file with the format:\n"
    "  <vertex x> <vertex y> <displacement x> <displacement y>\n"
    "Objects are read from stdin and written to stdout unless the filenames\n"
    "are given.\n",
    *argv,
    " -o tied.wlz -p points.tie myobj.wlz\n"
    "A thin plate spline basis function transform is computed from the\n"
    "tie points read from the file points.tie. This transform is then\n"
    "applied to the object is read from myobj.wlz. The resulting object\n"
    "is then written to tied.wlz.\n");
  }
  return(!ok);
}

/************************************************************************
* Function:	WlzMeshOutputPS						*
* Returns:	void							*
* Purpose:	Outputs the mesh and displace mesh of the given mesh	*
*		transform as a Postscript diagram to the standard error	*
*		output.							*
* Global refs:	-							*
* Parameters:	WlzMeshTransform *mesh:	Given mesh transform.		*
************************************************************************/
static void	WlzMeshOutputPS(WlzMeshTransform *mesh)
{
  int		elmCount,
  		nodCount;
  double	scale;
  WlzMeshNode	*nod;
  WlzMeshElem	*elm;
  WlzDVertex2	tVx,
  		offset;
  WlzDBox2	bBox;
  WlzDVertex2	elmVx[3];

  (void )fprintf(stderr, "%%!\n"
  			 "1 setlinecap\n"
			 "1 setlinejoin\n"
			 "0.01 setlinewidth\n");
  if(mesh && (mesh->type == WLZ_TRANSFORM_2D_MESH) &&
    (mesh->nNodes > 0) && (mesh->nElem > 0))
  {
    /* Compute the bounding box of the mesh transform. */
    nod = mesh->nodes;
    bBox.xMin = bBox.xMax = nod->position.vtX;
    bBox.yMin = bBox.yMax = nod->position.vtY;
    nodCount = mesh->nNodes;
    while(nodCount-- > 0)
    {
      tVx = nod->position;
      if(tVx.vtX < bBox.xMin)
      {
        bBox.xMin = tVx.vtX;
      }
      else if(tVx.vtX > bBox.xMax)
      {
        bBox.xMax = tVx.vtX;
      }
      if(tVx.vtY < bBox.yMin)
      {
        bBox.yMin = tVx.vtY;
      }
      else if(tVx.vtY > bBox.yMax)
      {
        bBox.yMax = tVx.vtY;
      }
      tVx.vtX += nod->displacement.vtX;
      tVx.vtY += nod->displacement.vtY;
      if(tVx.vtX < bBox.xMin)
      {
        bBox.xMin = tVx.vtX;
      }
      else if(tVx.vtX > bBox.xMax)
      {
        bBox.xMax = tVx.vtX;
      }
      if(tVx.vtY < bBox.yMin)
      {
        bBox.yMin = tVx.vtY;
      }
      else if(tVx.vtY > bBox.yMax)
      {
        bBox.yMax = tVx.vtY;
      }
      ++nod;
    }
    if(((tVx.vtX = fabs(bBox.xMax - bBox.xMin)) > 1.0) &&
       ((tVx.vtY = fabs(bBox.yMax - bBox.yMin)) > 1.0))
    {
      /* Compute scale and offset for an A4 page. */
      if((4.0 * tVx.vtX) > (3.0 * tVx.vtY))
      {
        scale = 500.0 / (bBox.xMax - bBox.xMin);
      }
      else
      {
        scale = 700.0 / (bBox.yMax - bBox.yMin);
      }
      offset.vtX = (bBox.xMin * scale) - 50;
      offset.vtY = (bBox.yMin * -scale) - 750;
      (void )fprintf(stderr, "255 0 0 setrgbcolor\n");
      elmCount = mesh->nElem;
      elm = mesh->elements;
      while(elmCount-- > 0)
      {
	nodCount = 3;
	while(--nodCount >= 0)
	{
	  nod = mesh->nodes + elm->nodes[nodCount];
	  elmVx[nodCount] = nod->position;
	  elmVx[nodCount].vtX = (elmVx[nodCount].vtX * scale) - offset.vtX;
	  elmVx[nodCount].vtY = (elmVx[nodCount].vtY * -scale) - offset.vtY;
	}
	(void )fprintf(stderr,
		       "%g %g moveto "
		       "%g %g lineto "
		       "%g %g lineto "
		       "%g %g lineto "
		       "stroke\n",
		       elmVx[0].vtX, elmVx[0].vtY,
		       elmVx[1].vtX, elmVx[1].vtY,
		       elmVx[2].vtX, elmVx[2].vtY,
		       elmVx[0].vtX, elmVx[0].vtY);
	++elm;
      }
      (void )fprintf(stderr, "0 0 255 setrgbcolor\n");
      elmCount = mesh->nElem;
      elm = mesh->elements;
      while(elmCount-- > 0)
      {
	nodCount = 3;
	while(--nodCount >= 0)
	{
	  nod = mesh->nodes + elm->nodes[nodCount];
	  elmVx[nodCount].vtX = nod->position.vtX + nod->displacement.vtX;
	  elmVx[nodCount].vtY = nod->position.vtY + nod->displacement.vtY;
	  elmVx[nodCount].vtX = (elmVx[nodCount].vtX * scale) - offset.vtX;
	  elmVx[nodCount].vtY = (elmVx[nodCount].vtY * -scale) - offset.vtY;
	}
	(void )fprintf(stderr,
		       "%g %g moveto "
		       "%g %g lineto "
		       "%g %g lineto "
		       "%g %g lineto "
		       "stroke\n",
		       elmVx[0].vtX, elmVx[0].vtY,
		       elmVx[1].vtX, elmVx[1].vtY,
		       elmVx[2].vtX, elmVx[2].vtY,
		       elmVx[0].vtX, elmVx[0].vtY);
	++elm;
      }
      (void )fprintf(stderr, "showpage\n");
    }
  }
}
