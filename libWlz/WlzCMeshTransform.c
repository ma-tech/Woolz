#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzCMeshTransform.c
* \author       Bill Hill
* \date         October 2004
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
* \brief	Functions for creating and applying 2D and 3D conforming
*		mesh transforms.
* \bug          None known.
*/
#include <stdio.h>
#include <float.h>
#include <Wlz.h>

typedef enum _WlzCMeshTransformType /* TODO move to WlzTransformType */
{
  WLZ_TRANSFORM_2D_CMESH,
  WLZ_TRANSFORM_3D_CMESH
} WlzCMeshTransformType;

typedef union _WlzCMeshP
{
  void		*v;
  WlzCMesh2D	*m2d;
  WlzCMesh3D	*m3d;
} WlzCMeshP;


/*!
* \struct	_WlzCMeshTransform
* \ingroup	WlzTransform
* \brief	A conforming mesh transform.
*		typedef: ::WlzCMeshTransform.
*/
typedef struct _WlzCMeshTransform
{
  int		type;			/*!< Type of transform. */
  int		linkcount;		/*!< Core. */
  void		*freeptr;		/*!< Core. */
  WlzCMeshP	mesh;			/*!< The conforming mesh. */
  AlcVector	*dspVec;		/*!< Vector for displacements. */
} WlzCMeshTransform;

static WlzErrorNum 		WlzCMeshTransMakeDispCb2D(
				  void *meshP,
				  void *nodP, 
				  void *mTrP);
static WlzErrorNum 		WlzCMeshTransMakeDisp2D(
				  WlzCMeshTransform *mTr,
				  WlzCMesh2D *mesh,
				  WlzCMeshNod2D *nod,
				  int idx);

/*!
* \return	New 2D conforming mesh transform.
* \ingroup	WlzTransform
* \brief	Creates a new conforming mesh transform.
* \param	type			Mesh type.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzCMeshTransform *WlzMakeCMeshTransform(WlzCMeshTransformType type,
					 WlzErrorNum *dstErr)
{
  WlzCMeshTransform *mTr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(type)
  {
    case WLZ_TRANSFORM_2D_CMESH:
      if((mTr = (WlzCMeshTransform *)
      	        AlcCalloc(1, sizeof(WlzCMeshTransform))) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      break;
    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mTr->type = type;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mTr);
}

/*!
* \return	New 2D conforming mesh transform.
* \ingroup	WlzTransform
* \brief	Creates a new 2D conforming mesh transform from the
*		given 2D conforming mesh.
* \param	mesh			Given 2D conforming mesh.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzCMeshTransform *WlzMakeCMeshTransform2D(WlzCMesh2D *mesh,
					WlzErrorNum *dstErr)
{
  int 		idN;
  WlzDomain	dom;
  WlzCMeshNod2D	*nod;
  WlzCMeshTransform *mTr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_TRI2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mTr = WlzMakeCMeshTransform(WLZ_TRANSFORM_2D_CMESH, &errNum);
  }
  /* Create vector for the displacements. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((mTr->dspVec = AlcVectorNew(1, sizeof(WlzDVertex2),
    				   mesh->res.nod.vec->blkSz, NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Assign the mesh to the transform. */
    dom.cm2 = mesh;
    (void )WlzAssignDomain(dom, NULL);
    mTr->mesh.m2d = mesh;
    /* Set the mesh node properties to be displacements for all existing
     * mesh nodes. */
    idN = 0;
    while((errNum == WLZ_ERR_NONE) && (idN < mesh->res.nod.maxEnt))
    {
      nod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      errNum = WlzCMeshTransMakeDisp2D(mTr, mesh, nod, idN);
      ++idN;
    }
  }
  /* Set a mesh new node callback to set the node property to be a
   * displacement for all new nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzCMeshAddNewNodCb2D(mesh, WlzCMeshTransMakeDispCb2D,
    				   mTr);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mTr);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Callback function for a 2D mesh which creates a displacement
*		property for a node.
* \param	meshP			Used to pass the 2D mesh.
* \param	nodP			Used to pass the 2D node.
* \param	mTrP			Used to pass the mesh transform.
*/
static WlzErrorNum WlzCMeshTransMakeDispCb2D(void *meshP,
					void *nodP, void *mTrP)
{
  WlzCMesh2D	*mesh;
  WlzCMeshNod2D	*nod;
  WlzCMeshTransform *mTr;
  WlzErrorNum errNum = WLZ_ERR_NONE;

  mesh = (WlzCMesh2D *)meshP;
  nod = (WlzCMeshNod2D *)nodP;
  mTr = (WlzCMeshTransform *)mTrP;
  if(mesh && nod && mTr && (nod->idx >= 0))
  {
    errNum = WlzCMeshTransMakeDisp2D(mTr, mesh, nod, nod->idx);
  }
  return(errNum);
}

/*!
* \return
* \ingroup
* \brief	Creates a displacement for the given 2D conforming mesh node.
* \param	mesh			Given 2D conforming mesh.
* \param	vec			Vector from which to allocate the
* 	 				displacement.
* \param	nod			Node to have displacement.
* \param	idx			Index of the node in it's vector,
*					must be valid.
*/
static WlzErrorNum WlzCMeshTransMakeDisp2D(WlzCMeshTransform *mTr,
					   WlzCMesh2D *mesh,
					   WlzCMeshNod2D *nod, int idx)
{
  WlzDVertex2	*dsp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dsp = (WlzDVertex2 *)AlcVectorItemGet(mTr->dspVec, idx)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    nod->prop = (void *)dsp;
    dsp->vtX = dsp->vtY = 0.0;
  }
  return(errNum);
}

#ifdef WLZ_CMESH_DEBUG_MAIN

/*!
* \return       void
* \ingroup      WlzTransform
* \brief        Debuging function for 2D mesh output in VTK format.
* \param        fP                      Given file pointer.
* \param        mesh                    Given mesh.
* \param	dspFlg			Non zero for displaced mesh output.
*/
void            WlzCMeshTransDbgOutVTK2D(FILE *fP, WlzCMesh2D *mesh,
					 int dspFlg)
{
  int           idE,
                idN,
                bCnt,
                nElm,
                nVElm,
                nNod;
  WlzDVertex2	*dsp;
  WlzCMeshElm2D  *elm;
  WlzCMeshNod2D *nod;

  if(mesh && (mesh->type == WLZ_CMESH_TRI2D) &&
    ((nNod = mesh->res.nod.maxEnt) > 0) &&
    ((nElm = mesh->res.elm.maxEnt) > 0))
  {
    (void )fprintf(fP,
                   "# vtk DataFile Version 1.0\n"
                   "WlzCMesh2D 2D\n"
                   "ASCII\n"
                   "DATASET POLYDATA\n"
                   "POINTS %d float\n",
                   nNod);
    for(idN = 0; idN < nNod; ++idN)
    {
      nod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      if(nod->idx >= 0)
      {
	if(dspFlg)
	{
	  dsp = (WlzDVertex2 *)(nod->prop);
	  (void )fprintf(fP, "%g %g 0\n",
			 nod->pos.vtX + dsp->vtX,
			 nod->pos.vtY + dsp->vtY);
	}
	else
	{
	  (void )fprintf(fP, "%g %g 0\n",
			 nod->pos.vtX, nod->pos.vtY);
        }
      }
      else
      {
        (void )fprintf(fP, "0 0 0\n");
      }
    }
    nVElm = 0;
    for(idE = 0; idE < nElm; ++idE)
    {
      elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
        ++nVElm;
      }
    }
    (void )fprintf(fP, "POLYGONS %d %d\n",
                   nVElm, nVElm * 4);
    for(idE = 0; idE < nElm; ++idE)
    {
      elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
        (void )fprintf(fP, "3 %d %d %d\n",
                       elm->edg[0].nod->idx, elm->edg[1].nod->idx,
                       elm->edg[2].nod->idx);
      }
    }
  }
}

extern char 	*optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		dspFlg = 0,
  		ok = 1,
  		option,
  		usage = 0;
  double	minElmSz = 25.0,
  		maxElmSz = 100.0;
  FILE		*fP = NULL;
  char		*inObjFileStr,
  		*outFileStr;
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*obj = NULL;
  WlzCMesh2D 	*mesh = NULL;
  WlzCMeshTransform *mTr = NULL;
  static char   optList[] = "dhm:M:o:";
  const char    inObjFileStrDef[] = "-",
  	        outFileStrDef[] = "-";

  opterr = 0;
  inObjFileStr = (char *)inObjFileStrDef;
  outFileStr = (char *)outFileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'd':
        dspFlg = 1;
	break;
      case 'm':
        if(sscanf(optarg, "%lg", &minElmSz) != 1)
	{
	  usage = 1;
	}
	break;
      case 'M':
        if(sscanf(optarg, "%lg", &maxElmSz) != 1)
	{
	  usage = 1;
	}
        break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  ok = usage == 0;
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
  }
  if(ok)
  {
    if((inObjFileStr == NULL) ||
       (*inObjFileStr == '\0') ||
       ((fP = (strcmp(inObjFileStr, "-")?
              fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((obj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to read object from file %s\n",
                     *argv, inObjFileStr);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  if(ok)
  {
    (void )WlzAssignObject(obj, NULL);
    mesh = WlzCMeshFromObj2D(obj, minElmSz, maxElmSz, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s Failed to create conforming mesh, %s.\n",
      		     argv[0],
		     errMsgStr);
    }
  }
  if(ok)
  {
    if((mTr = WlzMakeCMeshTransform2D(mesh, &errNum)) == NULL)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to make mesh transform from mesh, %s.\n",
		     argv[0]);
    }
  }
  if(ok)
  {
    if((fP = (strcmp(outFileStr, "-")?
	     fopen(outFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to open output file %s.\n",
		     argv[0], outFileStr);
    }
  }
  if(ok)
  {
    WlzCMeshTransDbgOutVTK2D(fP, mesh, dspFlg);
  }
  if(fP && strcmp(outFileStr, "-"))
  {
    (void )fclose(fP); fP = NULL;
  }
  (void )WlzFreeObj(obj);
  if(usage)
  {
    fprintf(stderr,
            "Usage: %s [-h] [-o<output file>] [-m#] [-M#] [<input object>]\n"
    	    "Computes a conforming mesh for the given input object.\n"
	    "Options are:\n"
	    "  -h  Help, prints this usage message.\n"
	    "  -o  Output file.\n"
	    "  -d  Output displaced mesh.\n"
	    "  -m  Minimum mesh element size.\n"
	    "  -M  Maximum mesh element size.\n",
	    argv[0]);

  }
  return(!ok);
}
#endif /* WLZ_CMESH_DEBUG_MAIN */
