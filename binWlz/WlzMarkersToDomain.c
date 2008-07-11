#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzMarkersToDomain_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzMarkersToDomain.c
* \author       Bill Hill
* \date         April 2003
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
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
* \brief	Creates a domain with a marker located at the position
* 		of each vertex read from a file.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzmarkerstodomain "WlzMarkersToDomain"
*/

/*!
\ingroup BinWlz
\defgroup wlzmarkerstodomain WlzMarkersToDomain
\par Name
WlzMarkersToDomain - creates a domain with a marker located at the position
		     of each vertex read from a file.
\par Synopsis
\verbatim
WlzMarkersToDomain [-h] [-o<output file>] [-2] [-3] [-s #] [-t <type>]
		   [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file name, default standard output.</td>
  </tr>
  <tr> 
    <td><b>-2</b></td>
    <td>Vertices and output domain are 2D.</td>
  </tr>
  <tr> 
    <td><b>-3</b></td>
    <td>Vertices and output domain are 3D.</td>
  </tr>
  <tr> 
    <td><b>-2</b></td>
    <td>marker size.</td>
  </tr>
  <tr> 
    <td><b>-t</b></td>
    <td>Marker type: Only valid type is 'sphere'.</td>
  </tr>
</table>
\par Description
Reads a list of 2D or 3D vertices either from the given file or the
standard input (default). These vertices are then used to create
either a 2D or 3D domain with a marker at the position of each vertex.
\par Examples
\verbatim
cat in.num
10 20 30
50 60 80
70 70 60

WlzMarkersToDomain -o out.wlz -3 -s 3 in.num
\endverbatim
Creates a new domain with spheres or radius 3 at coordinates (10,20,30),
(50,60,80) and (70,70,60) and writes the object to the file out.wlz.
\par File
\ref WlzMarkersToDomain.c "WlzMarkersToDomain.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <unistd.h>
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
  int		ok = 1,
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
    "  -t  Marker type: Only valid type is 'sphere'.\n"
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
        if(sscanf(lnBuf, "%lg %lg",
	          &(datV.d2.vtX), &(datV.d2.vtY)) != 2)
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
	  (datVP.i2 + datCnt)->vtX = ALG_NINT(datV.d2.vtX);
	  (datVP.i2 + datCnt)->vtY = ALG_NINT(datV.d2.vtY);
	  ++datCnt;
	}
	break;
      case 3:
	vType = WLZ_VERTEX_I3;
        if(sscanf(lnBuf, "%lg %lg %lg",
	          &(datV.d3.vtX), &(datV.d3.vtY), &(datV.d3.vtZ)) != 3)
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
	  (datVP.i3 + datCnt)->vtX = ALG_NINT(datV.d3.vtX);
	  (datVP.i3 + datCnt)->vtY = ALG_NINT(datV.d3.vtY);
	  (datVP.i3 + datCnt)->vtZ = ALG_NINT(datV.d3.vtZ);
	  ++datCnt;
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
