#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzExtFFNodeEle_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlzExtFF/WlzExtFFNodeEle.c
* \author       Bill Hill
* \date         September 2009
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
* \brief	Functions for reading and writting Woolz objects to and from
* 		the the two part .node/.ele tetrahedral mesh file format
* 		used by Jonathan Shewchuk in Stellar.
*
*		The file formats are as follows:
*		\verbatim
 		.node Files
 		<n nodes> <dimension> <n attributes> <n boundary markers>
 		1 <x_1> <y_1> <z_1>
 		...
 		<n nodes> <x_n> <y_n> <z_n>
 
 		.ele files
		<n elements> <n nodes per element> [n attributes>
		1 <n_i> <n_j> <n_k> ... <n_n> [<a_1>  ... <a_n>]
		...
		<n elmements> <n_o> <n_p> <n_q> ... <n_n> [<a_1>  ... <a_n>]

 		\endverbatim
*		Blank lines and comments are allowed with comments being
*		the remainder of any record floowing a '#'.
* \ingroup	WlzExtFF
* \todo         -
* \bug          None known.
*/

#include <ctype.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

static char			*WlzEffReadObjNodeEleRec(
				  FILE *fP,
				  char *buf,
				  int bufLen);
/*!
* \return       Woolz error number.
* \ingroup      WlzExtFF
* \brief        Builds the node/ele file names from the given file name.
*               These strings should be free'd using AlcFree() when
*               no longer required.
* \param        fileBody       		Dest ptr for the file body.
* \param        nodeFileName		Dest ptr for the '.node' file name.
* \param        eleFileName		Dest ptr for the '.ele' file name.
* \param        gvnFileName		Given file name with .node or no
*                                       extension.
*/
WlzErrorNum     WlzEffNodeEleFileNames(char **fileBody,
                                   char **nodeFileName,
                                   char **eleFileName,
                                   const char *gvnFileName)
{
  int           tI0;
  WlzErrorNum   errFlag = WLZ_ERR_MEM_ALLOC;

  tI0 = ((int )strlen(gvnFileName) + 5) * sizeof(char);
  if(((*fileBody = (char *)AlcMalloc(tI0)) != NULL) &&
     ((*nodeFileName = (char *)AlcMalloc(tI0)) != NULL) &&
     ((*eleFileName = (char *)AlcMalloc(tI0)) != NULL))
  {
    (void )strcpy(*fileBody, gvnFileName);
    if((tI0 = (int )strlen(*fileBody) - 5) >= 0)
    {
      if((strcmp(*fileBody +  tI0, ".node") == 0) ||
         (strcmp(*fileBody +  tI0, ".ele") == 0))
      {
        *(*fileBody +  tI0) = '\0';
      }
    }
    (void )sprintf(*nodeFileName, "%s.node", *fileBody);
    (void )sprintf(*eleFileName, "%s.ele", *fileBody);
    errFlag = WLZ_ERR_NONE;
  }
  return(errFlag);
}


/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from a pair of files using the node/ele
* 		format. The given file name is used to generate the '.node'
* 		and '.ele' filenames.
* \param	gvnFileName		Given file name.
* \param	dstErr			Destination error number ptr, may be
* 					NULL.
*/
WlzObject	*WlzEffReadObjNodeEle(const char *gvnFileName,
				      WlzErrorNum *dstErr)
{
  int		cnt,
		idx,
  		idE,
  		idN,
		bnd = 0,
		dim = 0,
		nAtr = 0,
		nBnd = 0,
  		nElm = 0,
  		nNod = 0,
		nodPerElm = 0;
  FILE		*fP = NULL;
  char		*fileName = NULL,
  		*nodeFileName = NULL,
		*eleFileName = NULL,
		*str;
  int 		eBuf[4];
  WlzVertexP    vBuf;
  WlzCMeshP	mesh;
  WlzObject	*obj = NULL;
  WlzDomain	dom;
  WlzValues	val;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char		cBuf[256];
  WlzCMeshNod2D	*nBuf2[4];
  WlzCMeshNod3D	*nBuf3[4];

  vBuf.v = NULL;
  mesh.v = NULL;
  dom.core = NULL;
  val.core = NULL;
  if((gvnFileName == NULL) || (*gvnFileName == '\0'))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    errNum = WlzEffNodeEleFileNames(&fileName, &nodeFileName, &eleFileName,
                                    gvnFileName);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((nodeFileName == NULL) || (*nodeFileName == '\0') ||
       ((fP = fopen(nodeFileName, "r")) == NULL))
    {
      errNum = WLZ_ERR_READ_EOF;
    }
#ifdef _WIN32
    if(fP != NULL)
    {
      if(_setmode(_fileno(fP), 0x8000) == -1)
      {
        errNum = WLZ_ERR_READ_EOF;
      }
    }
#endif
  }
  if(fP == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    /* Optional records may be comments starting with a '#' but apart
     * from these the first line should have four integers specifying
     * the number of nodes, the dimension, the number of attributes
     * and the number of boundary markers. */
    cnt = 0;
    do
    {
      if((str = WlzEffReadObjNodeEleRec(fP, cBuf, 256)) == NULL)
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else if(*str != '#')
      {
	if((cnt = sscanf(str, "%d %d %d %d",
	                 &nNod, &dim, &nAtr, &nBnd)) != 4)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
      }
    } while((errNum == WLZ_ERR_NONE) && (cnt != 4));
  }
  /* Check for reasonable number of nodes, etc.... */
  if(errNum == WLZ_ERR_NONE)
  {
    if((nNod <= 0) || ((dim != 2) && (dim != 3)) || (nAtr != 0) || (nBnd != 1))
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    idN = 0;
    if(dim == 2)
    {
      /* Create a new 2D constrained mesh. */
      mesh.m2 = WlzCMeshNew2D(&errNum);
      /* Read in the node positions into a temporary buffer computing their
       * bounding box and then create the nodes. */
      if(errNum == WLZ_ERR_NONE)
      {
	if((vBuf.v = AlcMalloc(sizeof(WlzDVertex2) * nNod)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
      while((errNum == WLZ_ERR_NONE) && (idN < nNod))
      {
	if((str = WlzEffReadObjNodeEleRec(fP, cBuf, 256)) == NULL)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	  break;
	}
	else if(*str != '#')
	{
	  if((sscanf(str, "%d %lg %lg %d",
		     &idx,
		     &(vBuf.d2[idN].vtX), &(vBuf.d2[idN].vtY),
		     &bnd) == 4) && (idx = idN + 1))
	  {
	    if(idN == 0)
	    {
	      mesh.m2->bBox.xMin = mesh.m2->bBox.xMax = vBuf.d2[idN].vtX;
	      mesh.m2->bBox.yMin = mesh.m2->bBox.yMax = vBuf.d2[idN].vtY;
	    }
	    else
	    {
	      if(vBuf.d2[idN].vtX < mesh.m2->bBox.xMin)
	      {
		mesh.m2->bBox.xMin = vBuf.d2[idN].vtX;
	      }
	      else if(vBuf.d2[idN].vtX > mesh.m2->bBox.xMax)
	      {
		mesh.m2->bBox.xMax = vBuf.d2[idN].vtX;
	      }
	      if(vBuf.d2[idN].vtY < mesh.m2->bBox.yMin)
	      {
		mesh.m2->bBox.yMin = vBuf.d2[idN].vtY;
	      }
	      else if(vBuf.d2[idN].vtY > mesh.m2->bBox.yMax)
	      {
		mesh.m2->bBox.yMax = vBuf.d2[idN].vtY;
	      }
	    }
	    ++idN;
	  }
	  else
	  {
	    errNum = WLZ_ERR_READ_INCOMPLETE;
	    break;
	  }
	}
      }
      /* Create the nodes. */
      if(errNum == WLZ_ERR_NONE)
      {
	if(AlcVectorExtendAndGet(mesh.m2->res.nod.vec, nNod) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  errNum = WlzCMeshReassignGridCells2D(mesh.m2, nNod);
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	for(idN = 0; idN < nNod; ++idN)
	{
	   nBuf2[0] = WlzCMeshNewNod2D(mesh.m2, vBuf.d2[idN], NULL);
	   nBuf2[0]->pos = vBuf.d2[idN];
	   nBuf2[0]->flags = 0;
	}
      }
    }
    else /* dim == 3 */
    {
      /* Create a new 3D constrained mesh. */
      mesh.m3 = WlzCMeshNew3D(&errNum);
      /* Read in the node positions into a temporary buffer computing their
       * bounding box and then create the nodes. */
      if(errNum == WLZ_ERR_NONE)
      {
	if((vBuf.v = AlcMalloc(sizeof(WlzDVertex3) * nNod)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
      while((errNum == WLZ_ERR_NONE) && (idN < nNod))
      {
	if((str = WlzEffReadObjNodeEleRec(fP, cBuf, 256)) == NULL)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	  break;
	}
	else if(*str != '#')
	{
	  if((sscanf(str, "%d %lg %lg %lg %d",
		     &idx,
		     &(vBuf.d3[idN].vtX), &(vBuf.d3[idN].vtY),
		     &(vBuf.d3[idN].vtZ),
		     &bnd) == 5) && (idx = idN + 1))
	  {
	    if(idN == 0)
	    {
	      mesh.m3->bBox.xMin = mesh.m3->bBox.xMax = vBuf.d3[idN].vtX;
	      mesh.m3->bBox.yMin = mesh.m3->bBox.yMax = vBuf.d3[idN].vtY;
	      mesh.m3->bBox.zMin = mesh.m3->bBox.zMax = vBuf.d3[idN].vtZ;
	    }
	    else
	    {
	      if(vBuf.d3[idN].vtX < mesh.m3->bBox.xMin)
	      {
		mesh.m3->bBox.xMin = vBuf.d3[idN].vtX;
	      }
	      else if(vBuf.d3[idN].vtX > mesh.m3->bBox.xMax)
	      {
		mesh.m3->bBox.xMax = vBuf.d3[idN].vtX;
	      }
	      if(vBuf.d3[idN].vtY < mesh.m3->bBox.yMin)
	      {
		mesh.m3->bBox.yMin = vBuf.d3[idN].vtY;
	      }
	      else if(vBuf.d3[idN].vtY > mesh.m3->bBox.yMax)
	      {
		mesh.m3->bBox.yMax = vBuf.d3[idN].vtY;
	      }
	      if(vBuf.d3[idN].vtZ < mesh.m3->bBox.zMin)
	      {
		mesh.m3->bBox.zMin = vBuf.d3[idN].vtZ;
	      }
	      else if(vBuf.d3[idN].vtZ > mesh.m3->bBox.zMax)
	      {
		mesh.m3->bBox.zMax = vBuf.d3[idN].vtZ;
	      }
	    }
	    ++idN;
	  }
	  else
	  {
	    errNum = WLZ_ERR_READ_INCOMPLETE;
	    break;
	  }
	}
      }
      /* Create the nodes. */
      if(errNum == WLZ_ERR_NONE)
      {
	if(AlcVectorExtendAndGet(mesh.m3->res.nod.vec, nNod) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  errNum = WlzCMeshReassignGridCells3D(mesh.m3, nNod);
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	for(idN = 0; idN < nNod; ++idN)
	{
	   nBuf3[0] = WlzCMeshNewNod3D(mesh.m3, vBuf.d3[idN], NULL);
	   nBuf3[0]->pos = vBuf.d3[idN];
	   nBuf3[0]->flags = 0;
	}
      }
    }
  }
  AlcFree(vBuf.v);
  /* Close node file. */
  if(fP)
  {
    (void )fclose(fP);
    fP = NULL;
  }
  /* Read in the elements from the .ele file and create them. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((eleFileName == NULL) || (*eleFileName == '\0') ||
       ((fP = fopen(eleFileName, "r")) == NULL))
    {
      errNum = WLZ_ERR_READ_EOF;
    }
#ifdef _WIN32
    if(fP != NULL)
    {
      if(_setmode(_fileno(fP), 0x8000) == -1)
      {
        errNum = WLZ_ERR_READ_EOF;
      }
    }
#endif
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Optional records may be comments starting with a '#' but apart
     * from these the first line should have four integers specifying
     * the number of nodes, the dimension, the number of attributes
     * and the number of boundary markers. */
    cnt = 0;
    do
    {
      if((str = WlzEffReadObjNodeEleRec(fP, cBuf, 256)) == NULL)
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else if(*str != '#')
      {
	if((cnt = sscanf(str, "%d %d %d",
	                 &nElm, &nodPerElm, &nAtr)) != 3)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
      }
    } while((errNum == WLZ_ERR_NONE) && (cnt != 3));
  }
  /* Check for reasonable number of nodes, etc.... */
  if(errNum == WLZ_ERR_NONE)
  {
    if((nElm <= 0) || (nodPerElm != dim + 1) || (nAtr != 0))
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    idE = 0;
    if(dim == 2)
    {
      if(AlcVectorExtendAndGet(mesh.m2->res.elm.vec, nElm) == NULL)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
      while((errNum == WLZ_ERR_NONE) && (idE < nElm))
      {
	if((str = WlzEffReadObjNodeEleRec(fP, cBuf, 256)) == NULL)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	  break;
	}
	else if(*str != '#')
	{
	  if(sscanf(str, "%d %d %d %d",
		    &idx, eBuf + 0, eBuf + 1, eBuf + 2) == 4)
	  {
	    if((idx != idE + 1) ||
	       (eBuf[0] <= 0) || (eBuf[0] > nNod) ||
	       (eBuf[1] <= 0) || (eBuf[1] > nNod) ||
	       (eBuf[2] <= 0) || (eBuf[2] > nNod))
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	      break;
	    }
	    for(idN = 0; idN < 4; ++idN)
	    {
	      nBuf2[idN] = (WlzCMeshNod2D *)
			   AlcVectorItemGet(mesh.m2->res.nod.vec,
			                    eBuf[idN] - 1);
	    }
	    (void )WlzCMeshNewElm2D(mesh.m2, nBuf2[0], nBuf2[2], nBuf2[1],
				    1, &errNum);
	    ++idE;
	  }
	  else
	  {
	    errNum = WLZ_ERR_READ_INCOMPLETE;
	    break;
	  }
        }
      }
      if(errNum == WLZ_ERR_NONE)
      {
	WlzCMeshUpdateMaxSqEdgLen2D(mesh.m2);
	dom.cm2 = mesh.m2;
	obj = WlzMakeMain(WLZ_CMESH_2D, dom, val, NULL, NULL, &errNum);
      }
    }
    else /* dim == 3 */
    {
      if(AlcVectorExtendAndGet(mesh.m3->res.elm.vec, nElm) == NULL)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
      while((errNum == WLZ_ERR_NONE) && (idE < nElm))
      {
	if((str = WlzEffReadObjNodeEleRec(fP, cBuf, 256)) == NULL)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	  break;
	}
	else if(*str != '#')
	{
	  if(sscanf(str, "%d %d %d %d %d",
		    &idx, eBuf + 0, eBuf + 1, eBuf + 2, eBuf + 3) == 5)
	  {
	    if((eBuf[0] <= 0) || (eBuf[0] > nNod) ||
	       (eBuf[1] <= 0) || (eBuf[1] > nNod) ||
	       (eBuf[2] <= 0) || (eBuf[2] > nNod) ||
	       (eBuf[3] <= 0) || (eBuf[3] > nNod))
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	      break;
	    }
	    for(idN = 0; idN < 4; ++idN)
	    {
	      nBuf3[idN] = (WlzCMeshNod3D *)
			   AlcVectorItemGet(mesh.m3->res.nod.vec,
			                    eBuf[idN] - 1);
	    }
	    (void )WlzCMeshNewElm3D(mesh.m3,
				    nBuf3[0], nBuf3[1], nBuf3[3], nBuf3[2], 1, &errNum);
	    ++idE;
	  }
	  else
	  {
	    errNum = WLZ_ERR_READ_INCOMPLETE;
	    break;
	  }
        }
      }
      if(errNum == WLZ_ERR_NONE)
      {
	WlzCMeshUpdateMaxSqEdgLen3D(mesh.m3);
	dom.cm3 = mesh.m3;
	obj = WlzMakeMain(WLZ_CMESH_3D, dom, val, NULL, NULL, &errNum);
      }
    }
  }
  /* Close element file. */
  if(fP)
  {
    (void )fclose(fP);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(dim == 2)
    {
      (void )WlzCMeshFree2D(mesh.m2);
    }
    else if(dim == 3)
    {
      (void )WlzCMeshFree3D(mesh.m3);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	Start of string having skipped past possible white space
* 		characters.
* \ingroup	WlzExtFF
* \brief	Gets the next record from a file.
* \param	fP			Input file.
* \param	buf			Buffer for input record.
* \param	bufLen			Buffer length.
*/
static char	*WlzEffReadObjNodeEleRec(FILE *fP, char *buf, int bufLen)
{
  char		*str = NULL;

  if((str = fgets(buf, bufLen, fP)) != NULL)
  {
    buf[bufLen - 1] = '\0';
    while((*str != '\0') && isspace(*str))
    {
      ++str;
    }
  }
  return(str);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz object to a pair of files
* 		using the node/ele two part mesh file format.  The
* 		given file name is used to generate the '.node' and '.ele'
* 		filenames.
* \param	gvnFileName		Given file name with '.node' or
* 					no extension.
* \param	obj			Given woolz object.
*/
WlzErrorNum	WlzEffWriteObjNodeEle(const char *gvnFileName, WlzObject *obj)
{
  int		cnt,
  		idE,
  		idN,
		bnd = 0,
		dim = 0,
  		nElm = 0,
		nNod = 0,
		maxElm = 0,
		maxNod = 0;
  char		*fileName = NULL,
  		*nodeFileName = NULL,
		*eleFileName = NULL;
  int		*nodTbl = NULL;
  FILE		*fP = NULL;
  WlzCMeshP	mesh;
  WlzCMeshElm2D	*elm2;
  WlzCMeshElm3D	*elm3;
  WlzCMeshNod2D	*nod2[4];
  WlzCMeshNod3D	*nod3[4];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  mesh.v = NULL;
  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((gvnFileName == NULL) || (*gvnFileName == '\0'))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    errNum = WlzEffNodeEleFileNames(&fileName, &nodeFileName, &eleFileName,
                                    gvnFileName);
  }
  /* Open the node file. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((nodeFileName == NULL) || (*nodeFileName == '\0') ||
       ((fP = fopen(nodeFileName, "w")) == NULL))
    {
      errNum = WLZ_ERR_READ_EOF;
    }
#ifdef _WIN32
    if(fP != NULL)
    {
      if(_setmode(_fileno(fP), 0x8000) == -1)
      {
        errNum = WLZ_ERR_READ_EOF;
      }
    }
#endif
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(obj->type)
    {
      case WLZ_CMESH_2D:
	if((mesh.m2 = obj->domain.cm2) == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if(mesh.m2->type == WLZ_CMESH_2D)
	{
	  dim = 2;
	  nNod = mesh.m2->res.nod.numEnt;
	  nElm = mesh.m2->res.elm.numEnt;
	  maxNod = mesh.m2->res.nod.maxEnt;
	  maxElm = mesh.m2->res.elm.maxEnt;
	}
	else
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
        break;
      case WLZ_CMESH_3D:
	if((mesh.m3 = obj->domain.cm3) == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if(mesh.m3->type == WLZ_CMESH_3D)
	{
	  dim = 3;
	  nNod = mesh.m3->res.nod.numEnt;
	  nElm = mesh.m3->res.elm.numEnt;
	  maxNod = mesh.m3->res.nod.maxEnt;
	  maxElm = mesh.m3->res.elm.maxEnt;
	}
	else
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
        break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((nodTbl = (int *)AlcMalloc(maxNod * sizeof(int))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Output comment identifying the file followed by the number of nodes,
   * dimension, the number of attributes and the number of boundary markers.
   * On following records output the node psoitions. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP, "#%s\n%d %d %d %d\n", nodeFileName, nNod, dim, 0, 1) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  /* Output node positions while building a node table to avoid deleted
   * nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    cnt = 0;
    if(dim == 2)
    {
      for(idN = 0; idN < maxNod; ++idN)
      {
	nod2[0] = (WlzCMeshNod2D *)AlcVectorItemGet(mesh.m2->res.nod.vec, idN);
	if(nod2[0]->idx >= 0)
	{
	  nodTbl[idN] = cnt++;
	  if(fprintf(fP, "%d %lg %lg %d\n",
		     cnt,
		     nod2[0]->pos.vtX, nod2[0]->pos.vtY,
		     bnd) <= 0)
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    break;
	  }
	}
      }
    }
    else /* dim == 3 */
    {
      for(idN = 0; idN < maxNod; ++idN)
      {
	nod3[0] = (WlzCMeshNod3D *)AlcVectorItemGet(mesh.m3->res.nod.vec, idN);
	if(nod3[0]->idx >= 0)
	{
	  nodTbl[idN] = cnt++;
	  if(fprintf(fP, "%d %lg %lg %lg %d\n",
		     cnt,
		     nod3[0]->pos.vtX, nod3[0]->pos.vtY, nod3[0]->pos.vtZ,
		     bnd) <= 0)
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    break;
	  }
	}
      }
    }
  }
  /* Close the node file. */
  if(fP)
  {
    (void )fclose(fP);
    fP = NULL;
  }
  /* Open the element file. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((eleFileName == NULL) || (*eleFileName == '\0') ||
       ((fP = fopen(eleFileName, "w")) == NULL))
    {
      errNum = WLZ_ERR_READ_EOF;
    }
#ifdef _WIN32
    if(fP != NULL)
    {
      if(_setmode(_fileno(fP), 0x8000) == -1)
      {
        errNum = WLZ_ERR_READ_EOF;
      }
    }
#endif
  }
  /* Output comment identifying the file followed by the number of elements,
   * # the number of nodes per element and the number of attributes. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP, "#%s\n%d %d %d\n", eleFileName, nElm, dim + 1, 0) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  /* Output the node indices for each element. */
  if(errNum == WLZ_ERR_NONE)
  {
    cnt = 0;
    if(dim == 2)
    {
      for(idE = 0; idE < maxElm; ++idE)
      {
	elm2 = (WlzCMeshElm2D *)AlcVectorItemGet(mesh.m2->res.elm.vec, idE);
	if(elm2->idx >= 0)
	{
	  nod2[0] = WLZ_CMESH_ELM2D_GET_NODE_0(elm2);
	  nod2[1] = WLZ_CMESH_ELM2D_GET_NODE_1(elm2);
	  nod2[2] = WLZ_CMESH_ELM2D_GET_NODE_2(elm2);
	  if(fprintf(fP, "%d %d %d %d\n",
		     ++cnt,
		     nodTbl[nod2[0]->idx] + 1,
		     nodTbl[nod2[1]->idx] + 1,
		     nodTbl[nod2[2]->idx] + 1) <= 0)
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    break;
	  }
	}
      }
    }
    else /* dim == 3 */
    {
      for(idE = 0; idE < maxElm; ++idE)
      {
	elm3 = (WlzCMeshElm3D *)AlcVectorItemGet(mesh.m3->res.elm.vec, idE);
	if(elm3->idx >= 0)
	{
	  nod3[0] = WLZ_CMESH_ELM3D_GET_NODE_0(elm3);
	  nod3[1] = WLZ_CMESH_ELM3D_GET_NODE_1(elm3);
	  nod3[2] = WLZ_CMESH_ELM3D_GET_NODE_2(elm3);
	  nod3[3] = WLZ_CMESH_ELM3D_GET_NODE_3(elm3);
	  if(fprintf(fP, "%d %d %d %d %d\n",
		     ++cnt,
		     nodTbl[nod3[0]->idx] + 1,
		     nodTbl[nod3[1]->idx] + 1,
		     nodTbl[nod3[2]->idx] + 1,
		     nodTbl[nod3[3]->idx] + 1) <= 0)
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    break;
	  }
	}
      }
    }
  }
  /* Close the element file. */
  if(fP)
  {
    (void )fclose(fP);
    fP = NULL;
  }
  AlcFree(nodTbl);
  return(errNum);
}
