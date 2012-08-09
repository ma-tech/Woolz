#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzMarkersToDomain_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzMarkersToDomain.c
* \author       Bill Hill
* \date         April 2003
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
* \brief	Creates a domain with a marker located at the position
* 		of each vertex read from a file.
* \ingroup	BinWlz
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

#define WLZ_CFP_READLN_LEN	(1024)

static WlzVertexP 		WlzMTDReadVtxArray(
				  FILE *fP,
				  int dim,
				  int *dstNVtx, 
				  WlzVertexType *dstVType,
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
	markerType = WlzStringToMarkerType(optarg, &errNum);
	if(errNum != WLZ_ERR_NONE)
	{
	  usage = 1;
	}
        break;
      case 'h': /* FALLTROUGH */
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
    obj = WlzMakeMarkers(vType, nVtx, vtx, markerType, markerSz,
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
    "                          [<input file>]\n"
    "Version: %s\n"
    "Options:\n"
    "  -h  Output this usage message.\n"
    "  -o  Output file name, default is the standard output.\n"
    "  -2  Vertices and output domain are 2D.\n"
    "  -3  Vertices and output domain are 3D.\n"
    "  -s  Marker size.\n"
    "  -t  Marker type: Only valid type is 'sphere'.\n"
    "Reads a list of 2D or 3D vertices either from the given file or the\n"
    "standard input (default). These vertices are then used to create\n"
    "either a 2D or 3D domain with a marker at the position of each vertex.\n",
    argv[0],
    WlzVersion());
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

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
