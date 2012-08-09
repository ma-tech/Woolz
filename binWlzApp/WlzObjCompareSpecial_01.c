#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzObjCompareSpecial_01_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/WlzObjCompareSpecial_01.c
* \author	Richard Baldock
* \date         June 2000
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
* \brief	Compares a pair of objects and prints some statistics.
* \ingroup	BinWlzApp
*
* \par Binary
* \ref wlzobjcomparespecial_01 "WlzObjCompareSpecial_01"
*/

/*!
\ingroup BinWlzApp
\defgroup wlzobjcomparespecial_01 WlzObjCompareSpecial_01
\par Name
WlzObjCompareSpecial_01 - compares a pair of objects and
                          prints some statistics.
\par Synopsis
\verbatim
WlzObjCompareSpecial_01  [-h]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Prints usage information.</td>
  </tr>
</table>
\par Description
WlzObjCompareSpecial_01 compares a pair of objects read from the standard
input and prints some statistics.
\par Examples
\verbatim
\endverbatim
\par File
\ref WlzObjCompareSpecial_01.c "WlzObjCompareSpecial_01.c"
\par See Also
\ref BinWlzApp "WlzIntro(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
#include <Wlz.h>

static void usage(
  char	*str)
{
  fprintf(stderr,
	  "Usage:\t%s\n"
	  "\tPrint out some stats for two input woolz objects read from\n"
	  "\tthe standard input.\n"
	  "Version: %s\n",
	  str,
	  WlzVersion());

  return;
}

double vtxDist(
  WlzDVertex3	vtx1,
  WlzDVertex3	vtx2)
{
  double dist=0;

  dist += (vtx1.vtX - vtx2.vtX) * (vtx1.vtX - vtx2.vtX);
  dist += (vtx1.vtY - vtx2.vtY) * (vtx1.vtY - vtx2.vtY);
  dist += (vtx1.vtZ - vtx2.vtZ) * (vtx1.vtZ - vtx2.vtZ);
  
  if( dist > 0.0 ){
    dist = sqrt(dist);
  }
  return dist;
}

int WlzEdgeVerticesPoly(
  WlzPolygonDomain	*poly,
  WlzDVertex3		**vtxs,
  int			numVtxs,
  WlzErrorNum		*dstErr)
{
  int		nVtxs=0;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzObject	*polyObj;
  int		i;

  /* convert to an 8-connected polyline */
  if( poly ){
    if((polyObj = WlzPolyTo8Polygon(poly, 1, &errNum)) != NULL){
      nVtxs = polyObj->domain.poly->nvertices;
      if( numVtxs > 0 ){
	*vtxs = AlcRealloc(*vtxs, sizeof(WlzDVertex3)*
			   (numVtxs + nVtxs));
      }
      else {
	numVtxs = 0;
	*vtxs = AlcMalloc(sizeof(WlzDVertex3)*nVtxs);
      }
      switch( polyObj->domain.poly->type ){
      case WLZ_POLYGON_INT:
	for(i=0; i < nVtxs; i++){
	  (*vtxs)[numVtxs+i].vtX = polyObj->domain.poly->vtx[i].vtX;
	  (*vtxs)[numVtxs+i].vtY = polyObj->domain.poly->vtx[i].vtY;
	  (*vtxs)[numVtxs+i].vtZ = 0.0;
	}
	break;

      case WLZ_POLYGON_FLOAT:
	for(i=0; i < nVtxs; i++){
	  (*vtxs)[numVtxs+i].vtX =
	    ((WlzFVertex2 *) polyObj->domain.poly->vtx)[i].vtX;
	  (*vtxs)[numVtxs+i].vtY =
	    ((WlzFVertex2 *) polyObj->domain.poly->vtx)[i].vtY;
	  (*vtxs)[numVtxs+i].vtZ = 0.0;
	}
	break;

      case WLZ_POLYGON_DOUBLE:
	for(i=0; i < nVtxs; i++){
	  (*vtxs)[numVtxs+i].vtX =
	    ((WlzDVertex2 *) polyObj->domain.poly->vtx)[i].vtX;
	  (*vtxs)[numVtxs+i].vtY =
	    ((WlzDVertex2 *) polyObj->domain.poly->vtx)[i].vtY;
	  (*vtxs)[numVtxs+i].vtZ = 0.0;
	}
	break;
      default:
	errNum = WLZ_ERR_POLYGON_TYPE;
        break;
      }
      WlzFreeObj(polyObj);
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  nVtxs += numVtxs;
  return nVtxs;
}

int WlzEdgeVerticesBound(
  WlzBoundList	*bound,
  WlzDVertex3	**vtxs,
  int		numVtxs,
  WlzErrorNum	*dstErr)
{
  int		nVtxs=numVtxs;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( bound ){
    nVtxs = WlzEdgeVerticesBound(bound->next, vtxs, nVtxs, &errNum);
    nVtxs = WlzEdgeVerticesBound(bound->down, vtxs, nVtxs, &errNum);
    nVtxs = WlzEdgeVerticesPoly(bound->poly, vtxs, nVtxs, &errNum);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return nVtxs;
}

int WlzEdgeVertices(
  WlzObject	*obj,
  WlzDVertex3	**vtxs,
  WlzErrorNum	*dstErr)
{
  int		numVtxs = 0, n;
  WlzDVertex3	*rtnVtxs=NULL, *tmpVtxs;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzObject	*obj1, *obj2;
  int		i, j;
  WlzValues	values;
  WlzPlaneDomain	*planedmn;

  /* check object */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
      if((obj1 = WlzObjToBoundary(obj, 1, &errNum)) != NULL){
	numVtxs = WlzEdgeVerticesBound(obj1->domain.b, &rtnVtxs,
				       0, &errNum);
	WlzFreeObj(obj1);
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      if( obj->domain.core == NULL){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else {
	switch( obj->domain.core->type ){
	case WLZ_PLANEDOMAIN_DOMAIN:
	  planedmn = obj->domain.p;
	  values.core = NULL;
	  for(i=planedmn->plane1; i <= planedmn->lastpl; i++){
	    if( planedmn->domains[i-planedmn->plane1].core ){
	      obj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ,
				 planedmn->domains[i-planedmn->plane1],
				 values, NULL, NULL, NULL);
	      if( (n = WlzEdgeVertices(obj2, &tmpVtxs, &errNum)) > 0){
		if( numVtxs > 0 ){
		  rtnVtxs = AlcRealloc(rtnVtxs,
				       sizeof(WlzDVertex3)*(numVtxs+n));
		}
		else {
		  rtnVtxs = AlcMalloc(sizeof(WlzDVertex3)*n);
		}
		for(j=0; j < n; j++){
		  tmpVtxs[j].vtZ = i;
		  rtnVtxs[numVtxs+j] = tmpVtxs[j];
		}
		AlcFree(tmpVtxs);
		numVtxs += n;
	      }
	      WlzFreeObj(obj2);
	    }
	  }
	  break;

	case WLZ_PLANEDOMAIN_POLYGON:
	case WLZ_PLANEDOMAIN_BOUNDLIST:
	  break;

	case WLZ_EMPTY_DOMAIN:
	  if( dstErr ){
	    *dstErr = errNum;
	  }
	  return numVtxs;

	default:
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
      }
      break;

    case WLZ_2D_POLYGON:
      numVtxs = WlzEdgeVerticesPoly(obj->domain.poly, &rtnVtxs, 0, &errNum);
      break;

    case WLZ_BOUNDLIST:
      numVtxs = WlzEdgeVerticesBound(obj->domain.b, &rtnVtxs, 0, &errNum);
      break;

    case WLZ_EMPTY_OBJ:
      if( dstErr ){
	*dstErr = errNum;
      }
      return numVtxs;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  *vtxs = rtnVtxs;
  return numVtxs;
}

int main(
  int   argc,
  char  **argv)
{
  WlzObject     	*obj1, *obj2;
  int			verboseFlg=0;
  char			*nameStr;
  double		xsize=1.0, ysize=1.0, zsize=1.0;
  int			i, j;
  double		dist, minDist;
  double		dist1, dist2, dist3, dist4;
  int			numVtxs1, numVtxs2;
  WlzDVertex3		*edgeVtxs1, *edgeVtxs2;
  WlzDVertex3		cmVtx1, cmVtx2;
  double		vol1, vol2;
  double		mass1, mass2;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check input arguments */
  nameStr = argv[0];
  argv++;
  argc--;
  while( argc > 0 )
  {
    switch( argv[0][1] )
    {
    case 'h':
      usage(nameStr);
      return 0;

    case 'v':
      verboseFlg = 1;
      break;

    default:
      usage(nameStr);
      return 1;
    }
    argc--;
    argv++;
  }

  /* read the objects */
  if((obj1 = WlzReadObj(stdin, NULL)) &&
     (obj2 = WlzReadObj(stdin, NULL)) ){

    if( obj1->type != obj2->type ){
      return 1;
    }
    if((obj1->type != WLZ_2D_DOMAINOBJ) &&
       (obj1->type != WLZ_3D_DOMAINOBJ)){
      usage(nameStr);
      return 1;
    }

    /* get edge vertices */
    numVtxs1 = WlzEdgeVertices(obj1, &edgeVtxs1, NULL);
    numVtxs2 = WlzEdgeVertices(obj2, &edgeVtxs2, NULL);

    /* renormalise with pixel sizes */
    if( obj1->type == WLZ_3D_DOMAINOBJ ){
      xsize = obj1->domain.p->voxel_size[0];
      ysize = obj1->domain.p->voxel_size[1];
      zsize = obj1->domain.p->voxel_size[2];

      for(i=0; i < numVtxs1; i++){
	edgeVtxs1[i].vtX *= xsize;
	edgeVtxs1[i].vtY *= ysize;
	edgeVtxs1[i].vtZ *= zsize;
      }
      for(i=0; i < numVtxs2; i++){
	edgeVtxs2[i].vtX *= xsize;
	edgeVtxs2[i].vtY *= ysize;
	edgeVtxs2[i].vtZ *= zsize;
      }
    }

    /* get centre of mass vertices */
    if( obj1->type == WLZ_2D_DOMAINOBJ ){
      WlzDVertex2	vtx2;

      vtx2 = WlzCentreOfMass2D(obj1, 0, &mass1, &errNum);
      cmVtx1.vtX = vtx2.vtX;
      cmVtx1.vtY = vtx2.vtY;
      cmVtx1.vtZ = 0.0;

      vtx2 = WlzCentreOfMass2D(obj2, 0, &mass2, &errNum);
      cmVtx2.vtX = vtx2.vtX;
      cmVtx2.vtY = vtx2.vtY;
      cmVtx2.vtZ = 0.0;
    }
    else {
      cmVtx1 = WlzCentreOfMass3D(obj1, 0, &mass1, &errNum);
      cmVtx2 = WlzCentreOfMass3D(obj2, 0, &mass2, &errNum);
      cmVtx1.vtX *= xsize;
      cmVtx1.vtY *= ysize;
      cmVtx1.vtZ *= zsize;
      cmVtx2.vtX *= xsize;
      cmVtx2.vtY *= ysize;
      cmVtx2.vtZ *= zsize;
    }

    /* find distance between centres of mass */
    dist1 = vtxDist(cmVtx1, cmVtx2);

    /* find cm 1 to surface 2 dist */
    if( numVtxs2 > 0 ){
      minDist = vtxDist(cmVtx1, edgeVtxs2[0]);
      for(j=0; j < numVtxs2; j++){
	dist = vtxDist(cmVtx1, edgeVtxs2[j]);
	if( dist < minDist ){
	  minDist = dist;
	}
      }
      dist2 = minDist;
      if( !WlzInsideDomain(obj2, cmVtx1.vtZ/zsize,
			   cmVtx1.vtY/ysize, cmVtx1.vtX/xsize,
			   &errNum) ){
	dist2 *= -1.0;
      }
    }
    else {
      dist2 = 0.0;
    }

    /* find surface 1 to cm 2 dist */
    if( numVtxs1 > 0 ){
      minDist = vtxDist(cmVtx2, edgeVtxs1[0]);
      for(j=0; j < numVtxs1; j++){
	dist = vtxDist(cmVtx2, edgeVtxs1[j]);
	if( dist < minDist ){
	  minDist = dist;
	}
      }
      dist3 = minDist;
      if( !WlzInsideDomain(obj1, cmVtx2.vtZ/zsize,
			   cmVtx2.vtY/ysize, cmVtx2.vtX/xsize,
			   &errNum) ){
	dist3 *= -1.0;
      }
    }
    else {
      dist3 = 0.0;
    }

    /* find min distance between surfaces */
    if( (numVtxs1 > 0) && (numVtxs2 > 0) ){
      minDist = vtxDist(edgeVtxs1[0], edgeVtxs2[0]);
      for(i=0; i < numVtxs1; i++){
	for(j=0; j < numVtxs2; j++){
	  dist = vtxDist(edgeVtxs1[i], edgeVtxs2[j]);
	  if( dist < minDist ){
	    minDist = dist;
	  }
	}
      }
      dist4 = minDist;
    }
    else {
      dist4 = 0.0;
    }

    /* get the volumes */
    if( obj1->type == WLZ_2D_DOMAINOBJ ){
      vol1 = WlzArea(obj1, &errNum);
      vol2 = WlzArea(obj2, &errNum);
    }
    else {
      vol1 = WlzVolume(obj1, &errNum);
      vol2 = WlzVolume(obj2, &errNum);
      vol1 *= xsize*ysize*zsize;
      vol2 *= xsize*ysize*zsize;
    }

    /* print it */
    fprintf(stdout, "%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n",
	    vol1, vol2, mass1, mass2, dist1, dist2, dist3, dist4);
  }

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
