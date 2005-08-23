#ifndef DOXYGEN_SHOULD_SKIP_THIS
#pragma ident "MRC HGU $Id$"
/*!
* \file         binWlz/WlzMarkersToDomain.c
* \author       Bill Hill
* \date         April 2003
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
* \brief	Simple binary which reads a list of 2D or 3D vertices
*		and then creates a domain with a marker located at each
*		vertex position.
* \todo         -
* \bug          None known.
*/
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <Wlz.h>

typedef enum _WlzMarkerType
{
  WLZ_MARKER_NONE,
  WLZ_MARKER_SPHERE
} WlzMarkerType;

#define WLZ_CFP_READLN_LEN	(1024)

static WlzMarkerType 		WlzMarkerTypeFromStr(
				  const char *markerStr,
				  WlzErrorNum *dstErr);
static WlzVertexP 		WlzMTDReadVtxArray(
				  FILE *fP,
				  int dim,
				  int *dstNVtx, 
				  WlzVertexType *dstVType,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzRasterizeMarkers(
				  WlzVertexType vType,
				  int nVtx,
				  WlzVertexP vtx,
				  WlzMarkerType mType,
				  int mSz,
				  WlzErrorNum *dstErr);

extern char 	*optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		idx,
  		ok = 1,
  		option,
  		usage = 0,
		dim = 2,
		nVtx = 0,
		markerSz = 1;
  FILE		*fP = NULL;
  char		*iFile,
  		*oObjFile;
  const char	*errMsg;
  WlzVertexType	vType;
  WlzVertexP	vtx;
  WlzObject	*obj = NULL;
  WlzMarkerType	markerType = WLZ_MARKER_SPHERE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  static char   optList[] = "h23o:s:t:";
  const char    iFileDef[] = "-",
  		oObjFileDef[] = "-";

  opterr = 0;
  vtx.v = NULL;
  iFile = (char *)iFileDef;
  oObjFile = (char *)oObjFileDef;
  while(ok && ((option = getopt(argc, argv, optList)) != EOF))
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
        oObjFile = optarg;
	break;
      case 's':
	if(sscanf(optarg, "%d", &markerSz) != 1)
	{
	  usage = 1;
	}
        break;
      case 't':
	markerType = WlzMarkerTypeFromStr(optarg, &errNum);
	if(errNum != WLZ_ERR_NONE)
	{
	  usage = 1;
	}
        break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  if(!usage)
  {
    if((oObjFile == NULL) || (*oObjFile == '\0') ||
       (iFile == NULL) || (*iFile == '\0'))
    {
      usage = 1;
    }
  }
  if(!usage && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
    }
    else
    {
      iFile = *(argv + optind);
    }
  }
  ok = !usage;
  if(ok)
  {
    if(strcmp(iFile, "-"))
    {
      if((fP = fopen(iFile, "r")) == NULL)
      {
        ok = 0;
	errNum = WLZ_ERR_READ_EOF;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to open vertices file %s (%s).\n",
		       *argv, iFile, errMsg);
      }
    }
    else
    {
      fP = stdin;
    }
  }
  if(ok)
  {
    vtx = WlzMTDReadVtxArray(fP, dim, &nVtx, &vType, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to read vertices from file %s (%s).\n",
		     *argv, iFile, errMsg);
    }
  }
  if(fP)
  {
    if(strcmp(iFile, "-"))
    {
      fclose(fP);
    }
    fP = NULL;
  }
  if(ok)
  {
    obj = WlzRasterizeMarkers(vType, nVtx, vtx, markerType, markerSz,
    			      &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: Failed to create domain from vertices (%s).\n",
      		     argv[0],
		     errMsg);
    }
  }
  if(ok)
  {
    if((fP = (strcmp(oObjFile, "-")? fopen(oObjFile, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to open output file %s.\n",
		     argv[0], oObjFile);
    }
  }
  if(ok)
  {
    errNum = WlzWriteObj(fP, obj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to write output object, %s.\n",
		     argv[0], errMsg);
    }
  }
  if(fP && strcmp(oObjFile, "-"))
  {
    (void )fclose(fP);
  }
  AlcFree(vtx.v);
  (void )WlzFreeObj(obj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-o<output file>] [-2] [-3] [-s #] [-t <type>]\n"
    "                         [<input file>]\n"
    "  -h  Output this usage message.\n"
    "  -o  Output file name, default is the standard output.\n"
    "  -2  Vertices and output domain are 2D.\n"
    "  -3  Vertices and output domain are 3D.\n"
    "  -s  Marker size.\n"
    "  -s  Marker type: Only valid type is 'sphere'.\n"
    "Reads a list of 2D or 3D vertices either from the given file or the\n"
    "standard input (default). These vertices are then used to create\n"
    "either a 2D or 3D domain with a marker at the position of each vertex.\n",
    argv[0]);
  }
  return(!ok);
}

/*!
* \return	New vertex array.
* \brief	Reads an array of vertices (2D or 3D) from the given file.
* \param	fP			File pointer.
* \param	dim			Dimension which MUST be either 2 or 3.
* \param	dstNVtx			Destination pointer for the number of
*					vertices read, MUST not be NULL.
* \param	dstVType		Destination pointer for the vertex
*					type, MUST not be NULL.
* \param	dstErr			Used to return Woolz error code,
*					MUST not be NULL.
*/
static WlzVertexP 	WlzMTDReadVtxArray(FILE *fP, int dim, int *dstNVtx, 
					   WlzVertexType *dstVType,
					   WlzErrorNum *dstErr)
{
  int		datCnt = 0,
  		datMax = 0;
  char		lnBuf[WLZ_CFP_READLN_LEN];
  WlzVertex	datV;
  WlzVertexP	datVP;
  WlzVertexType vType;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  datVP.v = NULL;
  while((errNum == WLZ_ERR_NONE) &&
	(fgets(lnBuf, WLZ_CFP_READLN_LEN, fP) != NULL))
  {
    lnBuf[WLZ_CFP_READLN_LEN - 1] = '\0';
    switch(dim)
    {
      case 2:
	vType = WLZ_VERTEX_I2;
        if(sscanf(lnBuf, "%d %d",
	          &(datV.i2.vtX), &(datV.i2.vtY)) != 2)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
	else
	{
	  if(datCnt >= datMax)
	  {
	    datMax = (datMax + 1024) * 2;
	    if((datVP.v = AlcRealloc(datVP.v,
	    			     sizeof(WlzIVertex2) * datMax)) == NULL)
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  *(datVP.i2 + datCnt++) = datV.i2;
	}
	break;
      case 3:
	vType = WLZ_VERTEX_I3;
        if(sscanf(lnBuf, "%d %d %d",
	          &(datV.i3.vtX), &(datV.i3.vtY), &(datV.i3.vtZ)) != 3)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
	else
	{
	  if(datCnt >= datMax)
	  {
	    datMax = (datMax + 1024) * 2;
	    if((datVP.v = AlcRealloc(datVP.v,
	    			     sizeof(WlzIVertex3) * datMax)) == NULL)
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  *(datVP.i3 + datCnt++) = datV.i3;
	}
	break;
      default:
        errNum = WLZ_ERR_PARAM_DATA;
	break;
    }
  }
  *dstErr = errNum;
  *dstVType = vType;
  *dstNVtx = datCnt;
  return(datVP);
}

/*!
* \return	Marker type.
* \brief	Gets a marker type from a string.
* \param	markerStr		Given marker string.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzMarkerType WlzMarkerTypeFromStr(const char *markerStr,
				          WlzErrorNum *dstErr)
{
  int		tI0;
  WlzMarkerType markerType = WLZ_MARKER_NONE;
  WlzErrorNum	errNum = WLZ_ERR_PARAM_TYPE;
  
  if(WlzStringMatchValue(&tI0, markerStr,
  			 "none", WLZ_MARKER_NONE,
  			 "sphere", WLZ_MARKER_SPHERE,
			 NULL))
  {
    markerType = tI0;
    errNum = WLZ_ERR_NONE;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(markerType);
}

/*!
* \return	New domain object without values.
* \brief	Constructs a domain from the union of marker domains with
*		a marker domain at each of the given vertex positions.
* \param	nVtx			Number of vertices.
* \param	vtx			Given vertices.
* \param	mType			Marker type.
* \param	mSz			Marker size.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzRasterizeMarkers(WlzVertexType vType,
				     int nVtx, WlzVertexP vtx,
				     WlzMarkerType mType, int mSz,
				     WlzErrorNum *dstErr)
{
  int		idx,
  		dim;
  WlzObjectType	oType;
  WlzObject	*mObj = NULL;
  WlzObject	*tObj[4];
  WlzIVertex3	off;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  tObj[0] = tObj[1] = tObj[2] = tObj[3] = NULL;
  if((nVtx <= 0) || (mSz <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if(vtx.v == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(vType)
    {
      case WLZ_VERTEX_I2:
	dim = 2;
	oType = WLZ_2D_DOMAINOBJ;
        break;
      case WLZ_VERTEX_I3:
        dim = 3;
	oType = WLZ_3D_DOMAINOBJ;
	break;
      default:
        errNum = WLZ_ERR_PARAM_DATA;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tObj[0] = WlzMakeSphereObject(oType, mSz, 0.0, 0.0, 0.0, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tObj[1] = WlzMakeEmpty(&errNum);
  }
  for(idx = 0; (idx < nVtx) && (errNum == WLZ_ERR_NONE); ++idx)
  {
    if(dim == 2)
    {
      off.vtX = (vtx.i2 + idx)->vtX;
      off.vtY = (vtx.i2 + idx)->vtY;
    }
    else
    {
      off = *(vtx.i3 + idx);
    }
    tObj[2] = WlzShiftObject(tObj[0], off.vtX, off.vtY, off.vtZ, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      tObj[3] = WlzUnion2(tObj[1], tObj[2], &errNum);
    }
    WlzFreeObj(tObj[1]); tObj[1] = NULL;
    WlzFreeObj(tObj[2]); tObj[2] = NULL;
    if(errNum == WLZ_ERR_NONE)
    {
      tObj[1] = tObj[3];
      tObj[3] = NULL;
    }
  }
  WlzFreeObj(tObj[0]);
  if(errNum == WLZ_ERR_NONE)
  {
    mObj = tObj[1];
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mObj);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
