#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzReadObj_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzReadObj.c
* \author       Richard Baldock, Bill Hill
* \date         March 1999
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
* \brief	Reads a Woolz object from a file stream.
* \ingroup	WlzIO
* \todo         - 
* \bug          None known.
*/

#include <stdlib.h>
#include <string.h>
#include <Wlz.h>

#ifdef HAVE_MMAP
#define WLZ_USE_MMAP
#include <errno.h>
#include <unistd.h>
#include <sys/mman.h>
#endif

/* #define WLZ_DEBUG_READOBJ */
#define WLZ_OLD_CMESH_TRANS_SUPPORT

#if defined(_WIN32) && !defined(__x86)
#define __x86
#endif

static WlzIntervalDomain 	*WlzReadIntervalDomain(
				  FILE *fp,
				  WlzErrorNum *);
static WlzPlaneDomain 		*WlzReadPlaneDomain(
				  FILE *fp,
				  WlzErrorNum *);
static WlzErrorNum		WlzReadGreyValues(
				  FILE *fp,
				  WlzObject *obj);
static WlzErrorNum		WlzReadRectVtb(
				  FILE *fp,
				  WlzObject *obj,
				  WlzObjectType type);
static WlzErrorNum 		WlzReadDomObjValues3D(
				  FILE *fP,
				  WlzObject *obj);
static WlzErrorNum 		WlzReadTiledValues(
				  FILE *fP,
				  WlzObject *obj,
				  int dim,
				  WlzObjectType type,
				  int map);
static WlzErrorNum		WlzReadVoxelValues(
				  FILE *fp,
				  WlzObject *obj);
static WlzProperty	 	WlzReadProperty(
				  FILE *fp,
				  WlzErrorNum *);
static WlzPropertyList	 	*WlzReadPropertyList(
				  FILE *fp,
				  WlzErrorNum *);
static WlzPolygonDomain 	*WlzReadPolygon(
				  FILE *fp,
				  WlzErrorNum *);
static WlzBoundList 		*WlzReadBoundList(
				  FILE *fp,
				  WlzErrorNum *);
static WlzIRect 		*WlzReadRect(
				  FILE *fp,
				  WlzErrorNum *);
static WlzHistogramDomain 	*WlzReadHistogramDomain(
				  FILE *fp,
				  WlzErrorNum *);
static WlzObject 		*WlzReadCompoundA(
				  FILE *fp,
				  WlzObjectType type,
				  WlzErrorNum *);
static WlzAffineTransform 	*WlzReadAffineTransform(
				  FILE *fp,
				  WlzErrorNum *);
static WlzWarpTrans 		*WlzReadWarpTrans(
				  FILE *fp,
				  WlzErrorNum *);
static WlzFMatchObj 		*WlzReadFMatchObj(
				  FILE *fp,
				  WlzErrorNum *);
static Wlz3DWarpTrans 		*WlzRead3DWarpTrans(
				  FILE *fp,
				  WlzErrorNum *);
static WlzGMModel 		*WlzReadGMModel(
				  FILE *fP,
				  WlzErrorNum *dstErr);
static WlzContour 		*WlzReadContour(
				  FILE *fP,
				  WlzErrorNum *dstErr);
static WlzErrorNum 		WlzReadInt(
				  FILE *fP,
				  int *iP,
				  int nI);
static WlzErrorNum 		WlzReadVertex2D(
				  FILE *fP,
				  WlzDVertex2 *vP,
				  int nV);
static WlzErrorNum 		WlzReadVertex2I(
				  FILE *fP,
				  WlzIVertex2 *vP,
				  int nV);
static WlzErrorNum 		WlzReadVertex3I(
				  FILE *fP,
				  WlzIVertex3 *vP,
				  int nV);
static WlzErrorNum 		WlzReadVertex3D(
				  FILE *fP,
				  WlzDVertex3 *vP,
				  int nV);
static WlzErrorNum 		WlzReadStr(
				  FILE *fP,
				  char **dstStr);
static WlzErrorNum 		WlzReadPixelV(FILE *fP,
				  WlzPixelV *pV,
				  int nPV);
static WlzErrorNum 		WlzReadGreyV(FILE *fP,
				  WlzGreyType type,
				  WlzGreyV *gV,
				  int nGV);
static WlzMeshTransform         *WlzReadMeshTransform2D(
				  FILE *fp,
				  WlzErrorNum *);
static WlzCMesh2D		*WlzReadCMesh2D(
				  FILE *fp,
				  WlzErrorNum *dstErr);
static WlzCMesh2D5		*WlzReadCMesh2D5(
				  FILE *fp,
				  WlzErrorNum *dstErr);
static WlzCMesh3D		*WlzReadCMesh3D(
				  FILE *fp,
				  WlzErrorNum *dstErr);
static WlzErrorNum 		WlzReadIndexedvValues(
				  FILE *fP,
				  WlzObject *obj);
#ifdef WLZ_OLD_CMESH_TRANS_SUPPORT
static WlzObject 		*WlzReadOldCMeshTransform(
				  FILE *fP,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzReadOldCMeshTransform2D(
				  FILE *fP,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzReadOldCMeshTransform3D(
				  FILE *fP,
				  WlzErrorNum *dstErr);
#endif /* WLZ_OLD_CMESH_TRANS_SUPPORT */

#ifdef WLZ_UNUSED_FUNCTIONS
static WlzErrorNum 		WlzReadBox2I(
				  FILE *fP,
				  WlzIBox2 *bP,
				  int nB);
static WlzErrorNum 		WlzReadBox2D(
				  FILE *fP,
				  WlzDBox2 *bP,
				  int nB);
static WlzErrorNum 		WlzReadBox3I(
				  FILE *fP,
				  WlzIBox3 *bP,
				  int nB);
static WlzErrorNum 		WlzReadBox3D(
				  FILE *fP,
				  WlzDBox3 *bP,
				  int nB);
#endif /* WLZ_UNUSED_FUNCTIONS */

/*!
* \return	The word value.
* \ingroup	WlzIO
* \brief	Reads the next word (int) from the input file converting
*		from DEC VAX(!) byte order.
* \param	fp			Input file.
*/

static int 	getword(FILE *fp)
{
  char cin[4], cout[4];
 fread(cin,sizeof(char),4,fp);

#if defined (__sparc) || defined (__mips) || defined (__ppc)
  cout[0] = cin[3];
  cout[1] = cin[2];
  cout[2] = cin[1];
  cout[3] = cin[0];
#endif /* __sparc || __mips */
#if defined (__x86) || defined (__alpha)
  cout[0] = cin[0];
  cout[1] = cin[1];
  cout[2] = cin[2];
  cout[3] = cin[3];
#endif /* __x86 || __alpha */
  return(*((int *) &cout[0]));
}

/*!
* \return	The short value.
* \ingroup	WlzIO
* \brief	Reads the next short from the input file converting
*		from DEC VAX(!) byte order.
* \param	fp			Input file.
*/
static int 	getshort(FILE *fp)
{
  unsigned char cin[2], cout[2];
  fread(cin,sizeof(char),2,fp);

#if defined (__sparc) || defined (__mips) || defined (__ppc)
  cout[0] = cin[1];
  cout[1] = cin[0];
#endif /* __sparc || __mips */
#if defined (__x86) || defined (__alpha)
  cout[0] = cin[0];
  cout[1] = cin[1];
#endif /* __x86 || __alpha */
  return((int) *((short *) &cout[0]));
}

/*!
* \return	The float value.
* \ingroup      WlzIO
* \brief	Reads the next float from the input file converting
*               from DEC VAX(!) byte order.
* \param	fp			Input file.
*/
static float 	getfloat(FILE *fp)
{
  char cin[4], cout[4];

fread(cin,sizeof(char),4,fp);

#if defined (__sparc) || defined (__mips) || defined (__ppc)
  cout[0] = (char )(cin[1] - 1);
  cout[1] = cin[0];
  cout[2] = cin[3];
  cout[3] = cin[2];
#endif /* __sparc || __mips */
#if defined (__x86) || defined (__alpha)
  cout[3] = (char )(cin[1] - 1);
  cout[2] = cin[0];
  cout[1] = cin[3];
  cout[0] = cin[2];
#endif /* __x86 || __alpha */
  return(*((float *) &cout[0]));
}

/*
*   Function   : getdouble						*
*   Synopsis   : get the next double from the input stream		*
*   Returns    : double:	value of next double on the input stream*
*   Parameters : FILE *fp:	input stream				*
*   Global refs: -                                                      *
*/

/*!
* \return	The  doublevalue.
* \ingroup      WlzIO
* \brief	Reads the next double from the input file converting
*               from DEC VAX(!) byte order.
* \param	fp			Input file.
*/
static double 	getdouble(FILE *fp)
{
  char cin[8], cout[8];
  fread(cin,sizeof(char),8,fp);

#if defined (__sparc) || defined (__mips) || defined (__ppc)
  cout[0] = cin[7];
  cout[1] = cin[6];
  cout[2] = cin[5];
  cout[3] = cin[4];
  cout[4] = cin[3];
  cout[5] = cin[2];
  cout[6] = cin[1];
  cout[7] = cin[0];
#endif /* __sparc || __mips */
#if defined (__x86) || defined (__alpha)
  cout[7] = cin[7];
  cout[6] = cin[6];
  cout[5] = cin[5];
  cout[4] = cin[4];
  cout[3] = cin[3];
  cout[2] = cin[2];
  cout[1] = cin[1];
  cout[0] = cin[0];
#endif /* __x86 || __alpha */
  return(*((double *) &cout[0]));
}

/*!
* \return	Woolz object type as read from file.
* \ingroup	WlzIO
* \brief	Reads the type of an object from the given input file stream.
* 		This type is returned as read from the file and may not be a
* 		valid object type if the file is not a Woolz file or it is
* 		corrupted.
* \param	fP			Input file.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObjectType	WlzReadObjType(FILE *fP, WlzErrorNum *dstErr)
{
  WlzObjectType	type;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(fP == NULL )
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
#ifdef _WIN32
  else if(_setmode(_fileno(fP), 0x8000) == -1)
  {
    errNum = WLZ_ERR_READ_EOF;
  }
#endif
  else if(feof(fP) != 0)
  {
    errNum = WLZ_ERR_READ_EOF;
  }
  else
  {
    type = (WlzObjectType )getc(fP);
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(type);
}

/*!
* \return	New Woolz object or NULL on error.
* \ingroup	WlzIO
* \brief	Reads a woolz object from the given input stream. For
*		some object types (e.g. 3D) an object may be returned with the
*		error set to WLZ_ERR_READ_INCOMPLETE. This allows partial
*		recovery of data.
* \param	fp			Input file.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 	*WlzReadObj(FILE *fp, WlzErrorNum *dstErr)
{
  WlzObjectType		type;
  WlzObject 		*obj;
  WlzDomain		domain;
  WlzValues		values;
  Wlz3DWarpTrans	*wtrans3d;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  obj = NULL;
  domain.core = NULL;
  values.core = NULL;
  type = WlzReadObjType(fp, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    switch(type)
    {

    case (WlzObjectType) EOF:
      errNum = WLZ_ERR_READ_EOF;
      break;

    case WLZ_NULL:	/* signifies no more objects - not a true error */
      errNum = WLZ_ERR_EOO;
      break;

    case WLZ_EMPTY_OBJ:
      obj =  WlzMakeMain(WLZ_EMPTY_OBJ, domain, values, NULL,
			 NULL, &errNum);
      break;

    case WLZ_2D_DOMAINOBJ:
      if(((domain.i = WlzReadIntervalDomain(fp, &errNum)) != NULL) &&
	 ((obj = WlzMakeMain(type, domain, values, NULL, NULL,
	                     &errNum)) != NULL))
      {
	if((errNum = WlzReadGreyValues(fp, obj)) == WLZ_ERR_NONE)
	{
	  obj->plist = WlzAssignPropertyList(WlzReadPropertyList(fp, NULL),
					     NULL);
	}
	/* Don't clear up on error return, instead partial object as it may be
	 * salvagable. */
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      if(((domain.p = WlzReadPlaneDomain(fp, &errNum)) != NULL) &&
	 ((obj = WlzMakeMain(type, domain, values, NULL, NULL,
	                     &errNum)) != NULL ))
      {
	if((errNum = WlzReadDomObjValues3D(fp, obj)) == WLZ_ERR_NONE)
	{
	  obj->plist = WlzAssignPropertyList(WlzReadPropertyList(fp, NULL),
					     NULL);
	}
	/* Don't clear up on error, instead return partial object as it may be
	 * salvagable. */
      }
      break;

    case WLZ_TRANS_OBJ:
      if((domain.t = WlzReadAffineTransform( fp, &errNum)) != NULL){
	if((values.obj = WlzReadObj( fp, &errNum)) != NULL){
	  if((obj = WlzMakeMain(WLZ_TRANS_OBJ, domain, values,
				NULL, NULL, &errNum)) != NULL){
	    obj->plist = WlzAssignPropertyList(WlzReadPropertyList(fp, NULL),
					       NULL);
	  } else {
	    WlzFreeAffineTransform(domain.t);
	    WlzFreeObj(values.obj);
	  }
	}
	else {
	  WlzFreeAffineTransform(domain.t);
	}
      }
      break;

    case WLZ_3D_WARP_TRANS:
      if((wtrans3d = WlzRead3DWarpTrans(fp, &errNum)) != NULL){
	wtrans3d->plist = WlzAssignPropertyList(
	  WlzReadPropertyList(fp, NULL), NULL);
      }
      obj = (WlzObject *) wtrans3d;
      break;

    case WLZ_2D_POLYGON:
      if((domain.poly = WlzReadPolygon(fp, &errNum)) != NULL){
	obj = WlzMakeMain(type, domain, values, NULL, NULL, &errNum);
      }
      break;

    case WLZ_BOUNDLIST:
      if((domain.b = WlzReadBoundList(fp, &errNum)) != NULL){
	obj = WlzMakeMain(type, domain, values, NULL, NULL, &errNum);
      }
      break;

    case WLZ_HISTOGRAM:
      if((domain.hist = WlzReadHistogramDomain(fp, &errNum)) != NULL){
	obj = WlzMakeMain(type, domain, values, NULL, NULL, &errNum);
      }
      break;

    case WLZ_CONTOUR:
      if((domain.ctr = WlzReadContour(fp, &errNum)) != NULL){
	obj = WlzMakeMain(type, domain, values, NULL, NULL, &errNum);
      }
      break;

    case WLZ_CMESH_2D:
      if(((domain.cm2 = WlzReadCMesh2D(fp, &errNum)) != NULL) &&
         ((obj = WlzMakeMain(type, domain, values, NULL, NULL,
	                     &errNum)) != NULL)){
        if((errNum = WlzReadIndexedvValues(fp, obj)) == WLZ_ERR_NONE){
	  obj->plist = WlzAssignPropertyList(
	               WlzReadPropertyList(fp, NULL), NULL);
	}
      }
      break;

    case WLZ_CMESH_2D5:
      if(((domain.cm2d5 = WlzReadCMesh2D5(fp, &errNum)) != NULL) &&
         ((obj = WlzMakeMain(type, domain, values, NULL, NULL,
	                     &errNum)) != NULL)){
        if((errNum = WlzReadIndexedvValues(fp, obj)) == WLZ_ERR_NONE){
	  obj->plist = WlzAssignPropertyList(
	               WlzReadPropertyList(fp, NULL), NULL);
	}
      }
      break;

    case WLZ_CMESH_3D:
      if(((domain.cm3 = WlzReadCMesh3D(fp, &errNum)) != NULL) &&
         ((obj = WlzMakeMain(type, domain, values, NULL, NULL,
	                     &errNum)) != NULL)){
        if((errNum = WlzReadIndexedvValues(fp, obj)) == WLZ_ERR_NONE){
	  obj->plist = WlzAssignPropertyList(
	               WlzReadPropertyList(fp, NULL), NULL);
	}
      }
      break;

    case WLZ_RECTANGLE:
      if((domain.r = WlzReadRect(fp, &errNum)) != NULL){
	obj = WlzMakeMain(type, domain, values, NULL, NULL, &errNum);
      }
      break;

    case WLZ_AFFINE_TRANS:
      if((domain.t = WlzReadAffineTransform(fp, &errNum)) != NULL){
	obj = WlzMakeMain(type, domain, values, NULL, NULL, &errNum);
      }
      break;

    case WLZ_WARP_TRANS:
      obj = (WlzObject *) WlzReadWarpTrans(fp, &errNum);
      break;

#ifdef WLZ_OLD_CMESH_TRANS_SUPPORT
    case WLZ_CMESH_TRANS:
      obj = WlzReadOldCMeshTransform(fp, &errNum);
      break;
#endif

    case WLZ_FMATCHOBJ:
      obj = (WlzObject *) WlzReadFMatchObj(fp, &errNum);
      break;

    case WLZ_COMPOUND_ARR_1:
    case WLZ_COMPOUND_ARR_2:
      obj = (WlzObject *) WlzReadCompoundA(fp, type, &errNum);
      break;

    case WLZ_PROPERTY_OBJ:
      obj = WlzMakeMain(type, domain, values,
			WlzReadPropertyList(fp, NULL), NULL, &errNum);
      break;
    case WLZ_MESH_TRANS:
      if((domain.mt = WlzReadMeshTransform2D(fp, &errNum)) != NULL){
	obj = WlzMakeMain(type, domain, values, NULL, NULL, &errNum);
      }
      break;
    /* orphans and not yet implemented object types for I/O */
    case WLZ_CONV_HULL:
    case WLZ_3D_POLYGON:
    case WLZ_CONVOLVE_INT:
    case WLZ_CONVOLVE_FLOAT:
    case WLZ_TEXT:
    case WLZ_COMPOUND_LIST_1:
    case WLZ_COMPOUND_LIST_2:
    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  if(dstErr){
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Read's 4 byte integers the given file stream into a
*               buffer of native ints (which must have room for at
*               least nI ints).
* \param	fP			Given file.
* \param	iP			Buffer for ints.
* \param	nI			Number of ints.
*/
static WlzErrorNum WlzReadInt(FILE *fP, int *iP, int nI)
{
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  while((errNum == WLZ_ERR_NONE) && (nI-- > 0))
  {
    *iP++ = getword(fP);
  }
  if(feof(fP))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzIO
* \brief	Read's 2D integer vertices from the given file
*               stream into a buffer (which must have room for
*               at least nV vertices).
* \param	fP			Given file.
* \param	vP			Buffer for 2D integer vertices.
* \param	nV			Number of vertices.
*/
static WlzErrorNum WlzReadVertex2I(FILE *fP, WlzIVertex2 *vP, int nV)
{
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  while((errNum == WLZ_ERR_NONE) && (nV-- > 0))
  {
    vP->vtY = getword(fP);
    vP->vtX = getword(fP);
    ++vP;
  }
  if(feof(fP))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzIO
* \brief	Read's 2D double vertices from the given file
*               stream into a buffer (which must have room for
*               at least nV vertices).
* \param	fP			Given file.
* \param	vP			Buffer for 2D integer vertices.
* \param	nV			Number of vertices.
*/
static WlzErrorNum WlzReadVertex2D(FILE *fP, WlzDVertex2 *vP, int nV)
{
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  while((errNum == WLZ_ERR_NONE) && (nV-- > 0))
  {
    vP->vtY = getdouble(fP);
    vP->vtX = getdouble(fP);
    ++vP;
  }
  if(feof(fP))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzIO
* \brief	Read's 3D integer vertices from the given file
*               stream into a buffer (which must have room for
*               at least nV vertices).
* \param	fP			Given file.
* \param	vP			Buffer for 3D integer vertices.
* \param	nV			Number of vertices.
*/
static WlzErrorNum WlzReadVertex3I(FILE *fP, WlzIVertex3 *vP, int nV)
{
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  while((errNum == WLZ_ERR_NONE) && (nV-- > 0))
  {
    vP->vtX = getword(fP);
    vP->vtY = getword(fP);
    vP->vtZ = getword(fP);
    ++vP;
  }
  if(feof(fP))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzIO
* \brief	Read's 3D integer vertices from the given file
*               stream into a buffer (which must have room for
*               at least nV vertices).
* \param	fP			Given file.
* \param	vP			Buffer for 3D integer vertices.
* \param	nV			Number of vertices.
*/
static WlzErrorNum WlzReadVertex3D(FILE *fP, WlzDVertex3 *vP, int nV)
{
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  while((errNum == WLZ_ERR_NONE) && (nV-- > 0))
  {
    vP->vtX = getdouble(fP);
    vP->vtY = getdouble(fP);
    vP->vtZ = getdouble(fP);
    ++vP;
  }
  if(feof(fP))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  return(errNum);
}

#ifdef WLZ_UNUSED_FUNCTIONS
/*!
* \return	Woolz error code.
* \ingroup      WlzIO
* \brief	Read's 2D integer boxes from the given file stream
*               into a buffer (which must have room for at least
*               nB bounding boxes).
* \param	fP			Given file.
* \param	bP			Ptr to 2D integer box.
* \param	nB			Number of bounding boxes.
*/
static WlzErrorNum WlzReadBox2I(FILE *fP, WlzIBox2 *bP, int nB)
{
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  while((errNum == WLZ_ERR_NONE) && (nB-- > 0))
  {
    bP->xMin = getword(fP);
    bP->yMin = getword(fP);
    bP->xMax = getword(fP);
    bP->yMax = getword(fP);
    ++bP;
  }
  if(feof(fP))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzIO
* \brief	Read's 2D double boxes from the given file stream
*               into a buffer (which must have room for at least
*               nB bounding boxes).
* \param	fP			Given file.
* \param	bP			Ptr to 2D integer box.
* \param	nB			Number of bounding boxes.
*/
static WlzErrorNum WlzReadBox2D(FILE *fP, WlzDBox2 *bP, int nB)
{
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  while((errNum == WLZ_ERR_NONE) && (nB-- > 0))
  {
    bP->xMin = getdouble(fP);
    bP->yMin = getdouble(fP);
    bP->xMax = getdouble(fP);
    bP->yMax = getdouble(fP);
    ++bP;
  }
  if(feof(fP))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzIO
* \brief	Read's 3D integer boxes from the given file stream
*               into a buffer (which must have room for at least
*               nB bounding boxes).
* \param	fP			Given file.
* \param	bP			Ptr to 3D integer box.
* \param	nB			Number of bounding boxes.
*/
static WlzErrorNum WlzReadBox3I(FILE *fP, WlzIBox3 *bP, int nB)
{
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  while((errNum == WLZ_ERR_NONE) && (nB-- > 0))
  {
    bP->xMin = getword(fP);
    bP->yMin = getword(fP);
    bP->zMin = getword(fP);
    bP->xMax = getword(fP);
    bP->yMax = getword(fP);
    bP->zMax = getword(fP);
    ++bP;
  }
  if(feof(fP))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzIO
* \brief	Read's 3D double boxes from the given file stream
*               into a buffer (which must have room for at least
*               nB bounding boxes).
* \param	fP			Given file.
* \param	bP			Ptr to 3D integer box.
* \param	nB			Number of bounding boxes.
*/
static WlzErrorNum WlzReadBox3D(FILE *fP, WlzDBox3 *bP, int nB)
{
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  while((errNum == WLZ_ERR_NONE) && (nB-- > 0))
  {
    bP->xMin = getdouble(fP);
    bP->yMin = getdouble(fP);
    bP->zMin = getdouble(fP);
    bP->xMax = getdouble(fP);
    bP->yMax = getdouble(fP);
    bP->zMax = getdouble(fP);
    ++bP;
  }
  if(feof(fP))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  return(errNum);
}
#endif /* WLZ_UNUSED_FUNCTIONS */

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Read a string as written by WlzWriteStr().
* \param	fP			Input file.
* \param	dstStr			Destination pointer for string.
*/
static WlzErrorNum WlzReadStr(FILE *fP, char **dstStr)
{
  int		len;
  char		*str = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(((len = getword(fP)) < 0) || feof(fP))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else if((str = (char *)AlcMalloc(sizeof(char) * (len + 1))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else if(fread(str, sizeof(char), len, fP) != len)
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
    AlcFree(str);
  }
  else
  {
    *(str + len) = '\0';
    *dstStr = str;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Read pixel values with type as written by WlzWritePixelV().
* \param	fP			Input file.
* \param	pV			Pixel values to set.
*/
static WlzErrorNum WlzReadPixelV(FILE *fP, WlzPixelV *pV, int nPV)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  while((errNum == WLZ_ERR_NONE) && (nPV-- > 0))
  {

    pV->type = (WlzGreyType )getc(fP);

    errNum = WlzReadGreyV(fP, pV->type, &(pV->v), 1);
    ++pV;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Read grey values as written by WlzWriteGreyV().
* \param	fP			Input file.
* \param	gType			Grey type.
* \param	gV			Grey values to set.
* \param	nGV			Number of grey values.
*/
static WlzErrorNum WlzReadGreyV(FILE *fP, WlzGreyType gType, WlzGreyV *gV,
				int nGV)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(gType)
  {
    case WLZ_GREY_INT:
      while(nGV-- > 0)
      {
	gV->inv = getword(fP);
	++gV;
      }
      break;
    case WLZ_GREY_SHORT:
      while(nGV-- > 0)
      {
	gV->shv = (short )getshort(fP);
	++gV;
      }
      break;
    case WLZ_GREY_UBYTE:
      while(nGV-- > 0)
      {
        gV->ubv = (WlzUByte )(getc(fP));
	++gV;
      }
      break;
    case WLZ_GREY_FLOAT:
      while(nGV-- > 0)
      {
	gV->flv = getfloat(fP);
	++gV;
      }
      break;
    case WLZ_GREY_DOUBLE:
      while(nGV-- > 0)
      {
	gV->dbv = getdouble(fP);
	++gV;
      }
      break;
    case WLZ_GREY_RGBA:
      while(nGV-- > 0)
      {
	gV->rgbv = (unsigned int )getword(fP);
	++gV;
      }
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
  }
  if((errNum == WLZ_ERR_NONE) && feof(fP))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  return(errNum);
}

/*!
* \return	New interval domain.
* \ingroup	WlzIO
* \brief	Reads a Woolz interval domain from the file.
* \param	fp			Given file.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzIntervalDomain *WlzReadIntervalDomain(FILE *fP,
						WlzErrorNum *dstErr)
{
  WlzObjectType		type;
  int			i, l, l1, ll, k1, kl, nints;
  WlzIntervalDomain	*idmn=NULL;
  WlzIntervalLine 	*ivln;
  WlzInterval 		*itvl0,
  			*itvl;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* read the type, currently WriteObj will write '\0' given a
     NULL pointer so we do the same here setting the error to be
     WLZ_ERR_EOO to distinguish it from an EOF error */

type = (WlzObjectType) getc(fP);

  if( type == (WlzObjectType) EOF ){
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else if( type == WLZ_NULL ){
    errNum = WLZ_ERR_EOO;
  }

  if( errNum == WLZ_ERR_NONE ){
    l1 = getword(fP);
    ll = getword(fP);
    k1 = getword(fP);
    kl = getword(fP);
    if( feof(fP) != 0 ){
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    else {
      idmn = WlzMakeIntervalDomain(type, l1, ll, k1, kl, &errNum);
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    switch (type) {

    case WLZ_INTERVALDOMAIN_INTVL:
      nints = 0;
      ivln = idmn->intvlines;
      for (l=l1; l<=ll; l++) {
	ivln->nintvs = getword(fP);
	nints += ivln->nintvs;
	ivln++;
      }
      if( feof(fP) != 0 ){
	WlzFreeIntervalDomain(idmn);
	idmn = NULL;
	errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }

      if( nints == 0 ){
	/* curious case of a no-intervals domain */
	WlzFreeIntervalDomain(idmn);
	idmn = NULL;
	errNum = WLZ_ERR_EOO;
	break;
      }
      else if( (itvl0 = (WlzInterval *)
	   AlcMalloc(nints * sizeof(WlzInterval))) == NULL){
	WlzFreeIntervalDomain(idmn);
	idmn = NULL;
	errNum = WLZ_ERR_MEM_ALLOC;
	break;
      }
      itvl = itvl0;
      ivln = idmn->intvlines;
      idmn->freeptr = AlcFreeStackPush(idmn->freeptr, (void *)itvl0, NULL);

      for (i=0; i<nints; i++,itvl++) {
	itvl->ileft = getword(fP);
	itvl->iright = getword(fP);
      }

      itvl = itvl0;

      if (feof(fP) != 0){
	WlzFreeIntervalDomain(idmn);
	idmn = NULL;
	errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }

      for (l=l1; l<=ll; l++) {
	nints = ivln->nintvs;
	errNum = WlzMakeInterval(l, idmn, nints, itvl);
	ivln++;
	itvl += nints;
      }
      errNum = WlzStandardIntervalDomain(idmn);
      break;

    case WLZ_INTERVALDOMAIN_RECT:
      break;

    default:
      /* this can't happen because the domain type has been checked
	 by WlzMakeIntervalDomain */
      break;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(idmn);
}

/*!
* \return	New plane domain.
* \ingroup	WlzIO
* \brief	Reads a Woolz plane domain from the given file.
* \param	fp			Given file.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzPlaneDomain *WlzReadPlaneDomain(FILE *fp,
					  WlzErrorNum *dstErr)
{
  WlzObjectType		type;
  WlzDomain		domain, *domains;
  WlzPlaneDomain	*planedm=NULL;
  int			i, nplanes;
  int			p1, pl, l1, ll, k1, kl;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  type = (WlzObjectType) getc(fp);


  if( type == (WlzObjectType) EOF ){
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else if( type == WLZ_NULL ){
    errNum = WLZ_ERR_EOO;
  }
  else if( type == (WlzObjectType) 2 ){
    /* some old object files have this old value - now converted
       silently but we must continue to read these files */
    type = WLZ_PLANEDOMAIN_DOMAIN;
  }

  if( errNum == WLZ_ERR_NONE ){
	p1 = getword(fp);
    pl = getword(fp);
    l1 = getword(fp);
    ll = getword(fp);
    k1 = getword(fp);
    kl = getword(fp);

    if(feof(fp) != 0){
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    else if((planedm = WlzMakePlaneDomain(type, p1, pl, l1, ll, k1, kl,
					  &errNum)) != NULL){
      domains = planedm->domains;
      (planedm->voxel_size)[0] = getfloat(fp);
      (planedm->voxel_size)[1] = getfloat(fp);
      (planedm->voxel_size)[2] = getfloat(fp);
      nplanes = pl - p1 + 1;

      /* object format includes redundant plane positions -
	 read and discard */
      for(i=0; i < nplanes; i++){
	(void) getfloat(fp);
      }
      if (feof(fp) != 0){
	WlzFreePlaneDomain(planedm);
	planedm = NULL;
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    switch (type) {

    case WLZ_PLANEDOMAIN_DOMAIN:
      for(i=0; i < nplanes; i++, domains++){
	if((domain.i = WlzReadIntervalDomain(fp, &errNum)) != NULL){
	  *domains = WlzAssignDomain(domain, NULL);
	} else if( errNum == WLZ_ERR_EOO ){
	  errNum = WLZ_ERR_NONE;
	} else {
	  break;
	}
      }
      break;

    case WLZ_PLANEDOMAIN_POLYGON:
      for(i=0; i < nplanes; i++, domains++){
	if((domain.poly = WlzReadPolygon(fp, &errNum)) != NULL){
	  *domains = WlzAssignDomain(domain, NULL);
	} else if( errNum == WLZ_ERR_EOO ){
	  errNum = WLZ_ERR_NONE;
	} else {
	  break;
	}
      }
      break;

    case WLZ_PLANEDOMAIN_BOUNDLIST:
      for(i=0; i < nplanes; i++, domains++){
	if((domain.b = WlzReadBoundList(fp, &errNum)) != NULL){
	  *domains = WlzAssignDomain(domain, NULL);
	} else if( errNum == WLZ_ERR_EOO ){
	  errNum = WLZ_ERR_NONE;
	} else {
	  break;
	}
      }
      break;

    case WLZ_PLANEDOMAIN_HISTOGRAM:
      for(i=0; i < nplanes; i++, domains++){
	if((domain.hist = WlzReadHistogramDomain(fp, &errNum)) != NULL){
	  *domains = WlzAssignDomain(domain, NULL);
	} else if( errNum == WLZ_ERR_EOO ){
	  errNum = WLZ_ERR_NONE;
	} else {
	  break;
	}
      }
      break;

    case WLZ_PLANEDOMAIN_AFFINE:
      for (i=0; i< nplanes ; i++, domains++){
	if((domain.t = WlzReadAffineTransform(fp, &errNum)) != NULL){
	  *domains = WlzAssignDomain(domain, NULL);
	} else if( errNum == WLZ_ERR_EOO ){
	  errNum = WLZ_ERR_NONE;
	} else {
	  break;
	}
      }
      break ;

    case WLZ_PLANEDOMAIN_WARP:
      for (i=0; i< nplanes ; i++, domains++){
	if((domain.wt = WlzReadWarpTrans(fp, &errNum)) != NULL){
	  *domains = WlzAssignDomain(domain, NULL);
	} else if( errNum == WLZ_ERR_EOO ){
	  errNum = WLZ_ERR_NONE;
	} else {
	  break;
	}
      }
      break ;

    default:
      /* this can't happen because the domain type has been checked
	 by WlzMakeIntervalDomain */
      break;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(planedm);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Reads a Woolz grey-value table from the given file.
* \param	fp			Input file.
* \param	obj			Object defining the domain of the
*					grey values.
*/
static WlzErrorNum WlzReadGreyValues(FILE *fp, WlzObject *obj)
{
  WlzObjectType		type;
  WlzGreyType		gtype;
  WlzIntervalWSpace 	iwsp;
  WlzValues		values;
  WlzGreyType		packing;
  int 			kstart, l1, ll, k1;
  int 			i;
  WlzPixelV 		backgrnd;
  WlzGreyP		v, g;
  size_t		table_size;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  type = (WlzObjectType) getc(fp);

  if( type == (WlzObjectType) EOF )
  {
    return WLZ_ERR_READ_INCOMPLETE;
  }
  else if( type == WLZ_NULL )
  {
    obj->values.core = NULL;
    return WLZ_ERR_NONE;
  }
  gtype = (WlzGreyType) type;

  /*
   * The "type" read from disc only coded the pixel type.
   * The shape type is determined by the object domain.
   * For the moment, choice is between standard ragged-rectangle
   * or true rectangle.
   */
  if( obj->domain.core == NULL ){
    return WLZ_ERR_DOMAIN_NULL;
  }

  switch( obj->domain.core->type ){

  case WLZ_INTERVALDOMAIN_INTVL:
    type = WlzGreyTableType(WLZ_GREY_TAB_RAGR, gtype, &errNum);
    break;

  case WLZ_INTERVALDOMAIN_RECT:
    type = WlzGreyTableType(WLZ_GREY_TAB_RECT, gtype, &errNum);
    break;

  default:
    obj->values.core = NULL;
    return( WLZ_ERR_DOMAIN_TYPE );

  }
  if( errNum != WLZ_ERR_NONE ){
    return errNum;
  }

  l1 = obj->domain.i->line1;
  ll = obj->domain.i->lastln;
  k1 = obj->domain.i->kol1;
  backgrnd.type = gtype;
  switch (type) {

  case WLZ_VALUETABLE_RAGR_INT:

    packing = (WlzGreyType) getc(fp);

    backgrnd.v.inv = getword(fp);

    /* create the value table */
    if( (values.v = WlzMakeValueTb(type, l1, ll, k1,
				   backgrnd, obj, &errNum)) == NULL ){
      return errNum;
    }
    values.v->width = obj->domain.i->lastkl - k1 + 1;
    obj->values = WlzAssignValues(values, NULL);

    /* allocate space for the pixel values, preset to background value */
    table_size = WlzLineArea(obj, NULL) * sizeof(int);
    if( (v.inp = (int *) AlcMalloc(table_size)) == NULL){
      WlzFreeValueTb(values.v);
      obj->values.v = NULL;
      return WLZ_ERR_MEM_ALLOC;
    }
    memset((void *) v.inp, backgrnd.v.inv, table_size);
    values.v->freeptr = AlcFreeStackPush(values.v->freeptr, (void *)v.inp,
    				         NULL);

    if( (errNum = WlzInitRasterScan(obj, &iwsp,
    				    WLZ_RASTERDIR_ILIC)) == WLZ_ERR_NONE ){
      while((errNum = WlzNextInterval(&iwsp)) == WLZ_ERR_NONE) {
	if (iwsp.nwlpos){
	  kstart = iwsp.lftpos;
	}
	switch (packing) {
	case WLZ_GREY_INT:
	  g.inp = v.inp+iwsp.lftpos-kstart;
	  for (i=0; i<iwsp.colrmn; i++){
	    *g.inp++ = getword(fp);
	  }
	  break;
	case WLZ_GREY_SHORT:
	  g.inp = v.inp+iwsp.lftpos-kstart;
	  for (i=0; i<iwsp.colrmn; i++){
	    *g.inp++ = getshort(fp);
	  }
	  break;
	case WLZ_GREY_UBYTE:
	  g.inp = v.inp+iwsp.lftpos-kstart;
	  for (i=0; i<iwsp.colrmn; i++){

		  *g.inp++ = getc(fp);

	  }
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
	}
	if (iwsp.intrmn == 0) {
	  (void) WlzMakeValueLine(values.v, iwsp.linpos, kstart,
				  iwsp.rgtpos, v.inp);
	  v.inp += (iwsp.rgtpos - kstart + 1);
	}
      }
      if (feof(fp) != 0){
	WlzFreeValueTb(values.v);
	obj->values.v = NULL;
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else if( errNum == WLZ_ERR_EOO ){
	errNum = WLZ_ERR_NONE;
      }
    }
    return errNum;

  case WLZ_VALUETABLE_RAGR_SHORT:

    packing = (WlzGreyType) getc(fp);

    backgrnd.v.shv = (short )getword(fp);

    /* create the value table */
    if( (values.v = WlzMakeValueTb(type, l1, ll, k1,
				   backgrnd, obj, &errNum)) == NULL ){
      return errNum;
    }
    values.v->width = obj->domain.i->lastkl - k1 + 1;
    obj->values = WlzAssignValues(values, NULL);

    /* allocate space for the pixel values, preset to background value */
    table_size = WlzLineArea(obj, NULL) * sizeof(short);
    if( (v.shp = (short *) AlcMalloc(table_size)) == NULL){
      WlzFreeValueTb(values.v);
      obj->values.v = NULL;
      return WLZ_ERR_MEM_ALLOC;
    }
    memset((void *) v.shp, (int) backgrnd.v.shv, table_size);
    values.v->freeptr = AlcFreeStackPush(values.v->freeptr, (void *)v.shp,
    					 NULL);

    if( (errNum = WlzInitRasterScan(obj, &iwsp,
    				    WLZ_RASTERDIR_ILIC)) == WLZ_ERR_NONE ){
      while((errNum = WlzNextInterval(&iwsp)) == WLZ_ERR_NONE) {
	if (iwsp.nwlpos){
	  kstart = iwsp.lftpos;
	}
	switch (packing) {
	case WLZ_GREY_SHORT:
	  g.shp = v.shp+iwsp.lftpos-kstart;
	  for (i=0; i<iwsp.colrmn; i++){
	    *g.shp++ = (short )getshort(fp);
	  }
	  break;
	case WLZ_GREY_UBYTE:
	  g.shp = v.shp+iwsp.lftpos-kstart;
	  for (i=0; i<iwsp.colrmn; i++){

		  *g.shp++ = (short )getc(fp);

	  }
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
	}
	if (iwsp.intrmn == 0) {
	  (void) WlzMakeValueLine(values.v, iwsp.linpos, kstart,
				  iwsp.rgtpos, v.inp);
	  v.shp += (iwsp.rgtpos - kstart + 1);
	}
      }
      if( feof(fp) != 0 ){
	WlzFreeValueTb(values.v);
	obj->values.v = NULL;
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else if( errNum == WLZ_ERR_EOO ){
	errNum = WLZ_ERR_NONE;
      }
    }
    return errNum;

  case WLZ_VALUETABLE_RAGR_UBYTE:

    packing = (WlzGreyType) getc(fp);

    backgrnd.v.ubv = (WlzUByte )getword(fp);

    /* create the value table */
    if( (values.v = WlzMakeValueTb(type, l1, ll, k1,
				   backgrnd, obj, &errNum)) == NULL ){
      return errNum;
    }
    values.v->width = obj->domain.i->lastkl - k1 + 1;
    obj->values = WlzAssignValues(values, NULL);

    /* allocate space for the pixel values, preset to background value */
    table_size = WlzLineArea(obj, NULL) * sizeof(WlzUByte);
    if( (v.ubp = (WlzUByte *) AlcMalloc(table_size)) == NULL){
      WlzFreeValueTb(values.v);
      obj->values.v = NULL;
      return WLZ_ERR_MEM_ALLOC;
    }
    memset((void *) v.ubp, (int) backgrnd.v.ubv, table_size);
    values.v->freeptr = AlcFreeStackPush(values.v->freeptr, (void *)v.ubp,
    					 NULL);

    if( (errNum = WlzInitRasterScan(obj, &iwsp,
    				    WLZ_RASTERDIR_ILIC)) == WLZ_ERR_NONE ){
      while((errNum = WlzNextInterval(&iwsp)) == WLZ_ERR_NONE) {
	if (iwsp.nwlpos){
	  kstart = iwsp.lftpos;
	}
	g.ubp = v.ubp+iwsp.lftpos-kstart;
	for (i=0; i<iwsp.colrmn; i++){

	  *g.ubp++ = (WlzUByte )getc(fp);

	}
	if (iwsp.intrmn == 0) {
	  (void) WlzMakeValueLine(values.v, iwsp.linpos, kstart,
				  iwsp.rgtpos, v.inp);
	  v.ubp += (iwsp.rgtpos - kstart + 1);
	}
      }
      if (feof(fp) != 0){
	WlzFreeValueTb(values.v);
	obj->values.v = NULL;
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else if( errNum == WLZ_ERR_EOO ){
	errNum = WLZ_ERR_NONE;
      }
    }
    return( WLZ_ERR_NONE );

  case WLZ_VALUETABLE_RAGR_FLOAT:

    packing = (WlzGreyType) getc(fp);

    backgrnd.v.flv = getfloat(fp);

    /* create the value table */
    if( (values.v = WlzMakeValueTb(type, l1, ll, k1,
				   backgrnd, obj, &errNum)) == NULL ){
      return errNum;
    }
    values.v->width = obj->domain.i->lastkl - k1 + 1;
    obj->values = WlzAssignValues(values, NULL);

    /* allocate space for the pixel values, preset to background */
    table_size = WlzLineArea(obj, NULL) * sizeof(float);
    if( (v.flp = (float *) AlcMalloc(table_size)) == NULL){
      WlzFreeValueTb(values.v);
      obj->values.v = NULL;
      return WLZ_ERR_MEM_ALLOC;
    }
    memset((void *) v.flp, (int) backgrnd.v.flv, table_size);
    values.v->freeptr = AlcFreeStackPush(values.v->freeptr, (void *)v.flp,
    					 NULL);

    if( (errNum = WlzInitRasterScan(obj, &iwsp,
    				    WLZ_RASTERDIR_ILIC)) == WLZ_ERR_NONE ){
      while((errNum = WlzNextInterval(&iwsp)) == WLZ_ERR_NONE) {
	if (iwsp.nwlpos){
	  kstart = iwsp.lftpos;
	}
	g.flp = v.flp+iwsp.lftpos-kstart;
	for (i=0; i<iwsp.colrmn; i++){
	  *g.flp++ = getfloat(fp);
	}
	if (iwsp.intrmn == 0) {
	  (void) WlzMakeValueLine(values.v, iwsp.linpos, kstart,
				  iwsp.rgtpos, v.inp);
	  v.flp += (iwsp.rgtpos - kstart + 1);
	}
      }
      if (feof(fp) != 0){
	WlzFreeValueTb(values.v);
	obj->values.v = NULL;
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else if( errNum == WLZ_ERR_EOO ){
	errNum = WLZ_ERR_NONE;
      }
    }
    return errNum;

  case WLZ_VALUETABLE_RAGR_DOUBLE:

    packing = (WlzGreyType) getc(fp);

    backgrnd.v.dbv = getdouble(fp);

    /* create the value table */
    if( (values.v = WlzMakeValueTb(type, l1, ll, k1,
				   backgrnd, obj, &errNum)) == NULL ){
      return errNum;
    }
    values.v->width = obj->domain.i->lastkl - k1 + 1;
    obj->values = WlzAssignValues(values, NULL);

    /* allocate space for the pixel values, preset to background */
    table_size = WlzLineArea(obj, NULL) * sizeof(double);
    if( (v.dbp = (double *) AlcMalloc(table_size)) == NULL){
      WlzFreeValueTb(values.v);
      obj->values.v = NULL;
      return WLZ_ERR_MEM_ALLOC;
    }
    memset((void *) v.dbp, (int) backgrnd.v.dbv, table_size);
    values.v->freeptr = AlcFreeStackPush(values.v->freeptr, (void *)v.dbp,
    					 NULL);

    if( (errNum = WlzInitRasterScan(obj, &iwsp,
    				    WLZ_RASTERDIR_ILIC)) == WLZ_ERR_NONE ){
      while((errNum = WlzNextInterval(&iwsp)) == WLZ_ERR_NONE) {
	if (iwsp.nwlpos){
	  kstart = iwsp.lftpos;
	}
	g.dbp = v.dbp+iwsp.lftpos-kstart;
	for (i=0; i<iwsp.colrmn; i++){
	  *g.dbp++ = getdouble(fp);
	}
	if (iwsp.intrmn == 0) {
	  (void) WlzMakeValueLine(values.v, iwsp.linpos, kstart,
				  iwsp.rgtpos, v.inp);
	  v.dbp += (iwsp.rgtpos - kstart + 1);
	}
      }
      if (feof(fp) != 0){
	WlzFreeValueTb(values.v);
	obj->values.v = NULL;
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else if( errNum == WLZ_ERR_EOO ){
	errNum = WLZ_ERR_NONE;
      }
    }
    return errNum;

  case WLZ_VALUETABLE_RAGR_RGBA:
    packing = (WlzGreyType) getc(fp);
    backgrnd.v.rgbv = getword(fp);

    /* create the value table */
    if( (values.v = WlzMakeValueTb(type, l1, ll, k1,
				   backgrnd, obj, &errNum)) == NULL ){
      return errNum;
    }
    values.v->width = obj->domain.i->lastkl - k1 + 1;
    obj->values = WlzAssignValues(values, NULL);

    /* allocate space for the pixel values, preset to background */
    table_size = WlzLineArea(obj, NULL) * sizeof(WlzUInt);
    if( (v.rgbp = (WlzUInt *) AlcMalloc(table_size)) == NULL){
      WlzFreeValueTb(values.v);
      obj->values.v = NULL;
      return WLZ_ERR_MEM_ALLOC;
    }
    memset((void *) v.rgbp, (int) backgrnd.v.rgbv, table_size);
    values.v->freeptr = AlcFreeStackPush(values.v->freeptr, (void *)v.rgbp,
    					 NULL);

    if( (errNum = WlzInitRasterScan(obj, &iwsp,
    				    WLZ_RASTERDIR_ILIC)) == WLZ_ERR_NONE ){
      while((errNum = WlzNextInterval(&iwsp)) == WLZ_ERR_NONE) {
	if (iwsp.nwlpos){
	  kstart = iwsp.lftpos;
	}
	g.rgbp = v.rgbp+iwsp.lftpos-kstart;
	for (i=0; i<iwsp.colrmn; i++){
	  *g.rgbp++ = getword(fp);
	}
	if (iwsp.intrmn == 0) {
	  (void) WlzMakeValueLine(values.v, iwsp.linpos, kstart,
				  iwsp.rgtpos, v.inp);
	  v.rgbp += (iwsp.rgtpos - kstart + 1);
	}
      }
      if (feof(fp) != 0){
	WlzFreeValueTb(values.v);
	obj->values.v = NULL;
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else if( errNum == WLZ_ERR_EOO ){
	errNum = WLZ_ERR_NONE;
      }
    }
    return errNum;

  case WLZ_VALUETABLE_RECT_INT:
  case WLZ_VALUETABLE_RECT_SHORT:
  case WLZ_VALUETABLE_RECT_UBYTE:
  case WLZ_VALUETABLE_RECT_FLOAT:
  case WLZ_VALUETABLE_RECT_DOUBLE:
  case WLZ_VALUETABLE_RECT_RGBA:
    return WlzReadRectVtb(fp, obj, type);

  default:
    /* this can't happen because the domain type has been checked
       by WlzGreyValuesTableType */
    return WLZ_ERR_VALUES_TYPE;
  }

}

/*!
* \return	Wolz error code.
* \ingroup	WlzIO
* \brief	Reads a Woolz rectangular grey table.
* \param	fp			Input file.
* \param	obj			Object defining the domain of the
*					grey values.
* \param	type			Grey table type - encodes greytype.
*/
static WlzErrorNum WlzReadRectVtb(FILE 		*fp,
				  WlzObject 	*obj,
				  WlzObjectType type)
{
  WlzGreyP		values;
  int 			i, num;
  WlzGreyType		packing;
  WlzIntervalDomain 	*idmn;
  WlzValues		vtb;
  WlzPixelV		bgd;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  packing = (WlzGreyType) getc(fp);

  if( (idmn = obj->domain.i) == NULL ){
    return( WLZ_ERR_DOMAIN_NULL );
  }

  bgd.type = WlzGreyTableTypeToGreyType( type, NULL );
  bgd.v.inv = 0;
  if((vtb.r = WlzMakeRectValueTb(type, idmn->line1, idmn->lastln,
				 idmn->kol1, idmn->lastkl-idmn->kol1 + 1,
				 bgd, NULL, &errNum)) == NULL ){
    return errNum;
  }
  num = vtb.r->width * (vtb.r->lastln - vtb.r->line1 + 1);

  /* test on pixel type to read background */
  switch( WlzGreyTableTypeToGreyType( type, NULL ) ){
  case WLZ_GREY_INT:
    vtb.r->bckgrnd.v.inv = getword(fp);
    values.inp = (int *) AlcMalloc(num*sizeof(int));
    break;
  case WLZ_GREY_SHORT:
    vtb.r->bckgrnd.v.shv = (short )getword(fp);
    values.shp = (short *) AlcMalloc(num*sizeof(short));
    break;
  case WLZ_GREY_UBYTE:
    vtb.r->bckgrnd.v.ubv = (WlzUByte )getword(fp);
    values.ubp = (WlzUByte *) AlcMalloc(num*sizeof(WlzUByte));
    break;
  case WLZ_GREY_FLOAT:
    vtb.r->bckgrnd.v.flv = getfloat(fp);
    values.flp = (float *) AlcMalloc(num*sizeof(float));
    break;
  case WLZ_GREY_DOUBLE:
    vtb.r->bckgrnd.v.dbv = getdouble(fp);
    values.dbp = (double *) AlcMalloc(num*sizeof(double));
    break;
  case WLZ_GREY_RGBA:
    vtb.r->bckgrnd.v.rgbv = getword(fp);
    values.rgbp = (WlzUInt *) AlcMalloc(num*sizeof(WlzUInt));
    break;
  default:
    return WLZ_ERR_GREY_TYPE;
    break;
  }

  if( values.inp == NULL ){
    WlzFreeValueTb(vtb.v);
    return WLZ_ERR_MEM_ALLOC;
  }
  vtb.r->freeptr = AlcFreeStackPush(vtb.r->freeptr, (void *)values.inp, NULL);
  vtb.r->values = values;
  obj->values = WlzAssignValues(vtb, NULL);

  switch( WlzGreyTableTypeToGreyType( type, NULL ) ) {

  case WLZ_GREY_INT:
    switch (packing) {
    case WLZ_GREY_INT:
      for (i=0 ; i<num ; i++)
	*values.inp++ = getword(fp);
      break;
    case WLZ_GREY_SHORT:
      for (i=0 ; i<num ; i++){
        short shv = (short )getshort(fp);
	*values.inp++ = shv;
/*fprintf(stderr, "%d\n", shv);*/}
      break;
    case WLZ_GREY_UBYTE:
		for (i=0 ; i<num ; i++){

	*values.inp++ = getc(fp);

		}
      break;
    default:
      break;
    }
    break;

  case WLZ_GREY_SHORT:
    switch (packing) {
    case WLZ_GREY_SHORT:
      for (i=0 ; i<num ; i++)
	*values.shp++ = (short )getshort(fp);
      break;
    case WLZ_GREY_UBYTE:
		for (i=0 ; i<num ; i++){

	*values.shp++ = (short )getc(fp);

		}
      break;
    default:
      break;
    }
    break;

  case WLZ_GREY_UBYTE:
    fread(values.ubp, sizeof(WlzUByte), num, fp);
    break;

  case WLZ_GREY_FLOAT:
    for (i=0 ; i<num ; i++)
      *values.flp++ = getfloat(fp);
    break;

  case WLZ_GREY_DOUBLE:
    for (i=0 ; i<num ; i++)
      *values.dbp++ = getdouble(fp);
    break;

  case WLZ_GREY_RGBA:
    for (i=0 ; i<num ; i++)
      *values.rgbp++ = getword(fp);
    break;

  default:
    return  WLZ_ERR_GREY_TYPE;
    break;
  }
  if( feof(fp) != 0 ){
    WlzFreeValueTb(vtb.v);
    obj->values.core = NULL;
    return  WLZ_ERR_READ_INCOMPLETE;
  }

  return WLZ_ERR_NONE;
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Reads a Woolz 3D value table for a 3D domain object.
* 		This function reads the grey value table type and then
* 		calls the appropriate function to complete the read.
* \param	fP			Input file.
* \param	obj			Object defining the domain of the
*					grey values. The domain is known to
*					be non NULL.
*/
static WlzErrorNum WlzReadDomObjValues3D(FILE *fP, WlzObject *obj)
{
  WlzObjectType	type;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  type = (WlzObjectType )getc(fP);
  if(type == (WlzObjectType )EOF)
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else if(type == WLZ_NULL)
  {
    obj->values.core = NULL;
  }
  else
  {
    switch(type)
    {
      case WLZ_VOXELVALUETABLE_GREY:
        errNum = WlzReadVoxelValues(fP, obj);
	break;
      case WLZ_VALUETABLE_TILED_INT:    /* FALLTHROUGH */
      case WLZ_VALUETABLE_TILED_SHORT:  /* FALLTHROUGH */
      case WLZ_VALUETABLE_TILED_UBYTE:  /* FALLTHROUGH */
      case WLZ_VALUETABLE_TILED_FLOAT:  /* FALLTHROUGH */
      case WLZ_VALUETABLE_TILED_DOUBLE: /* FALLTHROUGH */
      case WLZ_VALUETABLE_TILED_RGBA:
        errNum = WlzReadTiledValues(fP, obj, 3, type, 1);
	break;
      default:
        errNum = WLZ_ERR_VALUES_TYPE;
	break;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Reads a Woolz tiled value table from the input file.
* 		The table type has already been read and verified.
* \param	fp			Input file.
* \param	obj			Object defining the domain of the
*					grey values.
* \param	dim			The dimension of the tiled value
* 					table.
* \param	type			The grey value table type which
* 					encodes both the grey type and the
* 					value table type.
* \param	map			If non zero the tiles are memory
* 					mapped rather than read.
*/
static WlzErrorNum WlzReadTiledValues(FILE *fP, WlzObject *obj,
				      int dim, WlzObjectType type,
				      int map)
{
  WlzGreyType	gType;
  WlzTiledValues *tVal = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  gType = WlzGreyTableTypeToGreyType(type, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    int		tDim;

    if((tDim = getc(fP)) != dim)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tVal = WlzMakeTiledValues(dim, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tVal->type = type;
    tVal->dim = dim;
    tVal->kol1 = getword(fP);
    tVal->lastkl = getword(fP);
    tVal->line1 = getword(fP);
    tVal->lastln = getword(fP);
    tVal->plane1 = getword(fP);
    tVal->lastpl = getword(fP);
    errNum = WlzReadPixelV(fP, &(tVal->bckgrnd), 1);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int 	nIdx;

    tVal->tileSz = getword(fP);
    tVal->tileWidth = getword(fP);
    tVal->numTiles = getword(fP);
    tVal->nIdx[0] = getword(fP);
    tVal->nIdx[1] = getword(fP);
    tVal->nIdx[2] = getword(fP);
    nIdx = tVal->nIdx[0] * tVal->nIdx[1]  * tVal->nIdx[2];
    if((tVal->indices = (int *)AlcMalloc(nIdx * sizeof(int))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      errNum = WlzReadInt(fP, tVal->indices, nIdx);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    long long	off[2];

    off[0] = getword(fP);
    off[1] = getword(fP);
    if(feof(fP) != 0)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    else if(sizeof(long) < 8)
    {
      if(off[1] != 0)
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else
      {
        tVal->tileOffset = off[0];
      }
    }
    else
    {
      tVal->tileOffset = (long )(off[1] << 32) | (long )(off[0]);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    size_t	gSz,
      		tSz;

    gSz = WlzGreySize(gType);
    tSz = tVal->numTiles * tVal->tileSz;
    if(map == 0)
    {
      tVal->fd = -1;
      if((tVal->tiles.v = AlcMalloc(tSz * gSz)) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(fseek(fP, tVal->tileOffset, SEEK_SET) != 0)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	/* The tiles are stored using native byte ordering. */
	if(fread(tVal->tiles.v, gSz, tSz, fP) != tSz)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
      }
    }
    else
    {
#ifdef WLZ_USE_MMAP
      /* For mmap to work the file must have been opened either with
       * "rb" or "rb+" if the values in the file are to be modified. */
      if((tVal->fd = dup(fileno(fP))) < 0)
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else
      {
	tVal->tiles.v = mmap(NULL, tSz * gSz, PROT_READ | PROT_WRITE,
			     MAP_SHARED | MAP_FILE |MAP_NORESERVE,
			     tVal->fd, tVal->tileOffset);
	if(tVal->tiles.v == MAP_FAILED)
	{
	  tVal->tiles.v = mmap(NULL, tSz * gSz, PROT_READ,
			       MAP_PRIVATE | MAP_FILE |MAP_NORESERVE,
			       tVal->fd, tVal->tileOffset);
	}
	if(tVal->tiles.v == MAP_FAILED)
	{
	  tVal->fd = -1;
	  tVal->tiles.v = NULL;
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
      }
#else /* WLZ_USE_MMAP */
      tVal->tiles.v = NULL;
      errNum = WLZ_ERR_READ_INCOMPLETE;
#endif /* WLZ_USE_MMAP */
    }
  }
#ifdef WLZ_DEBUG_READOBJ
  if(tVal == NULL)
  {
    (void )fprintf(stderr, "WlzReadTiledValues() tVal == NULL\n");
  }
  else
  {
    (void )fprintf(stderr, "WlzReadTiledValues() tVal\n");
    (void )fprintf(stderr, "WlzReadTiledValues() tVal->type = %d\n",
                   (int )(tVal->type));
    (void )fprintf(stderr, "WlzReadTiledValues() tVal->linkcount = %d\n",
                   tVal->linkcount);
    (void )fprintf(stderr, "WlzReadTiledValues() tVal->dim = %d\n",
                   tVal->dim);
    (void )fprintf(stderr, "WlzReadTiledValues() tVal->kol1 = %d\n",
                   tVal->kol1);
    (void )fprintf(stderr, "WlzReadTiledValues() tVal->lastkl = %d\n",
                   tVal->lastkl);
    (void )fprintf(stderr, "WlzReadTiledValues() tVal->line1 = %d\n",
                   tVal->line1);
    (void )fprintf(stderr, "WlzReadTiledValues() tVal->lastln = %d\n",
                   tVal->lastln);
    (void )fprintf(stderr, "WlzReadTiledValues() tVal->plane1 = %d\n",
                   tVal->plane1);
    (void )fprintf(stderr, "WlzReadTiledValues() tVal->lastpl = %d\n",
                   tVal->lastpl);
    (void )fprintf(stderr, "WlzReadTiledValues() tVal->bckgrnd.type = %d\n",
                   (int )(tVal->bckgrnd.type));
    (void )fprintf(stderr, "WlzReadTiledValues() tVal->tileSz = %ld\n",
                   (long )(tVal->tileSz));
    (void )fprintf(stderr, "WlzReadTiledValues() tVal->tileWidth = %ld\n",
                   (long )tVal->tileWidth);
    (void )fprintf(stderr, "WlzReadTiledValues() tVal->numTiles = %ld\n",
                   (long )tVal->numTiles);
    (void )fprintf(stderr, "WlzReadTiledValues() tVal->nIdx = %d, %d, %d\n",
                   tVal->nIdx[0], tVal->nIdx[1],
		   (tVal->dim == 2)? 0: tVal->nIdx[2]);
    (void )fprintf(stderr, "WlzReadTiledValues() tVal->indices = 0x%lx\n",
                   (long )(tVal->indices));
    (void )fprintf(stderr, "WlzReadTiledValues() tVal->fd = %d\n",
                   tVal->fd);
    (void )fprintf(stderr, "WlzReadTiledValues() tVal->tileOffset = %ld\n",
                   tVal->tileOffset);
    (void )fprintf(stderr, "WlzReadTiledValues() tVal->tiles.v = 0x%lx\n",
                   (long )(tVal->tiles.v));
    (void )fprintf(stderr, "WlzReadTiledValues() errno = %s\n",
                   strerror(errno));
  }
#endif /* WLZ_DEBUG_READOBJ */
  if(errNum == WLZ_ERR_NONE)
  {
    obj->values.t = tVal;
    (void )WlzAssignValues(obj->values, NULL);
  }
  else
  {
    WlzFreeTiledValues(tVal);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Reads a Woolz voxel value table from the input file.
* 		The table type has already been read and verified.
* \param	fp			Input file.
* \param	obj			Object defining the domain of the
*					grey values.
*/
static WlzErrorNum WlzReadVoxelValues(FILE *fp, WlzObject *obj)
{
  int 			i, nplanes;
  WlzObject 		*tmpobj;
  WlzDomain 		*domains;
  WlzValues		*values, value;
  WlzPlaneDomain 	*planedm;
  WlzVoxelValues 	*voxtab;
  WlzPixelV		bgd;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  planedm = obj->domain.p;
  domains = planedm->domains;
  nplanes = planedm->lastpl - planedm->plane1 + 1;

  /* the background value is ill-used, currently the type is not
     written to disc therefore int assumed. The grey-types of the
     value tables are assumed to be the same therefore we reset the
     voxeltable background type and value to be that of any of the
     valuetable background values - should be changed when the file
     format is changed */
  bgd.type = WLZ_GREY_INT;
  bgd.v.inv = 0;

  if((voxtab = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
				   planedm->plane1, planedm->lastpl,
				   bgd, obj, &errNum)) != NULL){
    values = voxtab->values;
    voxtab->bckgrnd.v.inv = getword(fp);
  }
  else {
    return errNum;
  }

  for(i=0; i < nplanes; i++, values++, domains++){
    (*values).core = NULL;
    if((tmpobj = WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains, *values,
			     NULL, NULL, &errNum)) != NULL){
      if( (errNum = WlzReadGreyValues(fp, tmpobj)) == WLZ_ERR_NONE ){
	*values = WlzAssignValues(tmpobj->values, NULL);
	/* reset voxel-table background */
	if( (*values).core != NULL ){
	  switch( WlzGreyTableTypeToTableType((*values).core->type, NULL) ){
	  case WLZ_GREY_TAB_RAGR:
	    voxtab->bckgrnd = (*values).v->bckgrnd;
	    break;
	  case WLZ_GREY_TAB_RECT:
	    voxtab->bckgrnd = (*values).r->bckgrnd;
	    break;
	  case WLZ_GREY_TAB_INTL:
	    voxtab->bckgrnd = (*values).i->bckgrnd;
	    break;
	  default:
	    return WLZ_ERR_VALUES_TYPE;
	    break;
	  }
	}
      }
      else {
	(*values).core = NULL;
	WlzFreeDomain(*domains);
	(*domains).core = NULL;
      }
      WlzFreeObj( tmpobj );
    }
  }

  if( feof(fp) != 0 ){
    /* allow incomplete object - set domains of unread
       valuetables to empty */
    /* WlzFreeVoxelValueTb( voxtab );*/
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  value.vox = voxtab;
  obj->values = WlzAssignValues(value, NULL);

  return errNum;
}

/*!
* \return	New Woolz property.
* \ingroup	WlzIO
* \brief	Reads a single property from the input file.
*		Unknown properties are ignored in the hope that at least
*		the objects domain and values can be recovered.
* \param	fp			input file.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzProperty WlzReadProperty(
  FILE		*fp,
  WlzErrorNum	*dstErr)
{
  WlzObjectType	type;
  WlzProperty	rtnProp;
  int		si;
  char		*name;
  WlzPixelV	pV;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  rtnProp.core = NULL;

  type = (WlzObjectType )getc(fp);

  switch( type ){

  case (WlzObjectType) EOF:
    errNum = WLZ_ERR_READ_INCOMPLETE;
    break;

    /* a NULL property could be allowed */
  case WLZ_NULL:
    break;

  case WLZ_PROPERTY_SIMPLE:
    /* read size */
    si = getword(fp);
    if( feof(fp) != 0 || si < 0 ){
      errNum = WLZ_ERR_READ_INCOMPLETE;
      break;
    }

    /* create property list with space for the data */
    if( (rtnProp.simple = WlzMakeSimpleProperty(si, &errNum)) == NULL ){
      break;
    }

    /* The size is now correct for the amount of data */
    if( si > 0 ){
      fread(rtnProp.simple->prop, si, 1, fp);
    }
    if( feof(fp) != 0 ){
      WlzFreeSimpleProperty( rtnProp.simple );
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    break;

  case WLZ_PROPERTY_EMAP:
    /* create an empty property */
    if( (rtnProp.emap =
	 WlzMakeEMAPProperty(WLZ_EMAP_PROPERTY_GREY_MODEL,
			     NULL, NULL, NULL, NULL, NULL,
			     NULL, NULL, NULL, NULL, NULL,
			     &errNum)) == NULL ){
      break;
    }

    /* read the property values */

    rtnProp.emap->emapType = (WlzEMAPPropertyType )getc(fp);

    fread(rtnProp.emap->modelUID,
	  EMAP_PROPERTY_UID_LENGTH, 1, fp);
    fread(rtnProp.emap->anatomyUID,
	  EMAP_PROPERTY_UID_LENGTH, 1, fp);
    fread(rtnProp.emap->targetUID,
	  EMAP_PROPERTY_UID_LENGTH, 1, fp);
    fread(rtnProp.emap->targetVersion,
	  EMAP_PROPERTY_VERSION_LENGTH, 1, fp);
    fread(rtnProp.emap->stage,
	  EMAP_PROPERTY_STAGE_LENGTH, 1, fp);
    fread(rtnProp.emap->subStage,
	  EMAP_PROPERTY_STAGE_LENGTH, 1, fp);
    fread(rtnProp.emap->modelName,
	  EMAP_PROPERTY_MODELNAME_LENGTH, 1, fp);
    fread(rtnProp.emap->version,
	  EMAP_PROPERTY_VERSION_LENGTH, 1, fp);
    rtnProp.emap->creationTime = getword(fp);
    fread(rtnProp.emap->creationAuthor,
	  EMAP_PROPERTY_AUTHORNAME_LENGTH, 1, fp);
    fread(rtnProp.emap->creationMachineName,
	  EMAP_PROPERTY_MACHINENAME_LENGTH, 1, fp);
    rtnProp.emap->modificationTime = getword(fp);
    fread(rtnProp.emap->modificationAuthor,
	  EMAP_PROPERTY_AUTHORNAME_LENGTH, 1, fp);

    /* now the variable length bits */
    si = getword(fp);
    if( si > 0 ){
      rtnProp.emap->fileName = (char *) AlcMalloc(sizeof(char) *
						  (si+1));
      fread(rtnProp.emap->fileName, si, 1, fp);
      rtnProp.emap->fileName[si] = '\0';
      rtnProp.emap->freeptr = AlcFreeStackPush(rtnProp.emap->freeptr,
					       rtnProp.emap->fileName,
					       NULL);
    }
    else {
      rtnProp.emap->fileName = NULL;
    }

    si = getword(fp);
    if( si > 0 ){
      rtnProp.emap->comment = (char *) AlcMalloc(sizeof(char) *
						 (si+1));
      fread(rtnProp.emap->comment, si, 1, fp);
      rtnProp.emap->comment[si] = '\0';
      rtnProp.emap->freeptr = AlcFreeStackPush(rtnProp.emap->freeptr,
					       rtnProp.emap->comment,
					       NULL);
    }
    else {
      rtnProp.emap->comment = NULL;
    }
    break;
  case WLZ_PROPERTY_NAME:
    if((errNum = WlzReadStr(fp, &name)) == WLZ_ERR_NONE)
    {
      rtnProp.name = WlzMakeNameProperty(name, &errNum);
      AlcFree(name);
    }
    break;
  case WLZ_PROPERTY_GREY:
    if((errNum = WlzReadStr(fp, &name)) == WLZ_ERR_NONE)
    {
      if((errNum = WlzReadPixelV(fp, &pV, 1)) == WLZ_ERR_NONE)
      {
	rtnProp.greyV = WlzMakeGreyProperty(name, pV, &errNum);
      }
      AlcFree(name);
    }
    break;
  default:
    break;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnProp;
}

/*!
* \return	New property list.
* \ingroup	WlzIO
* \brief	Reads a Woolz property list from the input file.
* \param	fp			input file.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzPropertyList *WlzReadPropertyList(FILE *fp, WlzErrorNum *dstErr)
{
  int 		pSz,
  		numProps = 0;
  WlzObjectType	type;
  WlzProperty	prop;
  WlzPropertyList *pList = NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* Find number of properties */

  type = (WlzObjectType )getc(fp);

  switch(type)
  {
    case (WlzObjectType )EOF:
      errNum = WLZ_ERR_READ_INCOMPLETE;
      break;
    case WLZ_NULL:
      errNum = WLZ_ERR_EOO;
      break;
    case (WlzObjectType )1:
      /* For backward compatibility this corresponds to the old format
       * for a simple property so convert here to a property list. */
      /* Read size: Old format had a funny size definition. */
      pSz = getword(fp) - sizeof(int);
      if((feof(fp) != 0) || (pSz < 0))
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }
      /* Create property with space for the data */
      if((prop.simple = WlzMakeSimpleProperty(pSz, &errNum)) != NULL)
      {
	/* The size is now correct for the amount of data */
	if(pSz > 0)
	{
	  (void )fread(prop.simple->prop, pSz, 1, fp);
	}
	if(feof(fp) != 0)
	{
	  WlzFreeSimpleProperty(prop.simple);
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	  prop.core = NULL;
	}
      }
      /* Create a list with the simple property as it's only item. */
      if(errNum == WLZ_ERR_NONE)
      {
	(void )WlzAssignProperty(prop, NULL);
	if(((pList = WlzMakePropertyList(NULL)) == NULL) ||
	    (AlcDLPListEntryAppend(pList->list, NULL, (void *)(prop.core),
				   WlzFreePropertyListEntry) != ALC_ER_NONE))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	  (void )WlzFreePropertyList(pList);
	  pList = NULL;
	}
      }
      break;
    case (WlzObjectType )2:
      /* The new style property list. */
      numProps = getword(fp);
      if((feof(fp) != 0) || (numProps < 0))
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      /* Make a property list without any items. */
      if((errNum == WLZ_ERR_NONE) && (numProps > 0))
      {
	if((pList = WlzMakePropertyList(NULL)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
      /* Now read each property in turn and append it to the list. */
      while((errNum == WLZ_ERR_NONE) && (numProps-- > 0))
      {
	prop = WlzReadProperty(fp, &errNum);
	if(prop.core)
	{
	  (void )WlzAssignProperty(prop, NULL);
	  if(AlcDLPListEntryAppend(pList->list, NULL, (void *)prop.core,
				   WlzFreePropertyListEntry) != ALC_ER_NONE)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	    (void )WlzFreePropertyList(pList);
	    pList = NULL;
	  }
	}
      }
      break;
    default:
      errNum = WLZ_ERR_PROPERTY_TYPE;
      break;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(pList);
}

/*!
* \return	New polygon domain.
* \ingroup	WlzIO
* \brief	Reads a Woolz polygon domain from the input file.
* \param	fp			Input file.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzPolygonDomain *WlzReadPolygon(FILE *fp, WlzErrorNum *dstErr)
{
  WlzObjectType		type;
  WlzPolygonDomain	*poly=NULL;
  WlzFVertex2		*fvtx;
  WlzDVertex2		*dvtx;
  int 			nvertices, i;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  type = (WlzObjectType) getc(fp);

  if( type == (WlzObjectType) EOF ){
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  if( type == WLZ_NULL ){
    errNum = WLZ_ERR_EOO;
  }

  nvertices = getword(fp);
  if( feof(fp) != 0 ){
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }

  if((errNum == WLZ_ERR_NONE) &&
     (poly = WlzMakePolygonDomain(type, 0, NULL, nvertices, 1, &errNum)) ){
    poly->nvertices = nvertices;
    switch (type) {

    case WLZ_POLYGON_INT:
      for(i=0; i < nvertices; i++){
	((poly->vtx)+i)->vtY = getword(fp);
	((poly->vtx)+i)->vtX = getword(fp);
      }
      break;

    case WLZ_POLYGON_FLOAT:
      fvtx = (WlzFVertex2 *) poly->vtx;
      for(i=0; i < nvertices; i++){
	((fvtx)+i)->vtY = getfloat(fp);
	((fvtx)+i)->vtX = getfloat(fp);
      }
      break;

    case WLZ_POLYGON_DOUBLE:
      dvtx = (WlzDVertex2 *) poly->vtx;
      for(i=0; i < nvertices; i++){
	((dvtx)+i)->vtY = getdouble(fp);
	((dvtx)+i)->vtX = getdouble(fp);
      }
      break;
    default:
      errNum = WLZ_ERR_POLYGON_TYPE;
      break;
    }

    if( feof(fp) != 0 ){
      WlzFreePolyDmn( poly );
      poly = NULL;
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(poly);
}

/*!
* \return	New boundary list.
* \ingroup	WlzIO
* \brief	Reads a Woolz boundary list from the input file.
* \param	fp			Input file.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzBoundList *WlzReadBoundList(FILE *fp, WlzErrorNum *dstErr)
{
  WlzObjectType	type;
  WlzBoundList	*blist=NULL, *tmpblist;
  WlzPolygonDomain	*tmppoly;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  type = (WlzObjectType) getc(fp);

  switch( type ){

  case (WlzObjectType) EOF:
    errNum = WLZ_ERR_READ_INCOMPLETE;
    break;

  case WLZ_NULL:
    errNum = WLZ_ERR_EOO;
    break;

    /* dummy type written by WriteBoundList, real type read next */
  case (WlzObjectType) 1:

    type = (WlzObjectType) getc(fp);

    if( (blist = WlzMakeBoundList(type, 0, NULL, &errNum)) == NULL ){
      break;
    }

    if((tmpblist = WlzReadBoundList(fp, &errNum)) != NULL){
      blist->next = WlzAssignBoundList(tmpblist, NULL);
    }
    else if( errNum == WLZ_ERR_EOO){
      blist->next = NULL;
      errNum = WLZ_ERR_NONE;
    }
    else {
      WlzFreeBoundList(blist);
      blist = NULL;
      break;
    }

    if((tmpblist = WlzReadBoundList(fp, &errNum)) != NULL){
      blist->down = WlzAssignBoundList(tmpblist, NULL);
    }
    else if( errNum == WLZ_ERR_EOO){
      blist->down = NULL;
      errNum = WLZ_ERR_NONE;
    }
    else {
      WlzFreeBoundList(blist);
      blist = NULL;
      break;
    }

    blist->wrap = getword(fp);
    if((tmppoly = WlzReadPolygon(fp, &errNum)) != NULL){
      blist->poly = WlzAssignPolygonDomain(tmppoly, NULL);
    }
    else if( errNum == WLZ_ERR_EOO){
      blist->poly = NULL;
      errNum = WLZ_ERR_NONE;
    }
    else{
      WlzFreeBoundList(blist);
      blist = NULL;
      break;
    }

    if( feof(fp) != 0 ){
      WlzFreeBoundList(blist);
      blist = NULL;
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    break;

  default:
    break;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return( blist );
}

/*!
* \return	New Woolz rectangle.
* \ingroup	WlzIO
* \brief	Reads a Woolz rectangle from the given file.
* \param	fp			Input file.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzIRect *WlzReadRect(FILE *fp, WlzErrorNum *dstErr)
{
  WlzObjectType	type;
  WlzIRect 	*ir=NULL;
  WlzFRect 	*fr;
  int		i;
  WlzErrorNum	errNum=WLZ_ERR_NONE;


  type = (WlzObjectType) getc(fp);

  switch( type ){

  case (WlzObjectType) EOF:
    errNum = WLZ_ERR_READ_INCOMPLETE;
    break;

  case WLZ_NULL:
    errNum = WLZ_ERR_EOO;
    break;

  case WLZ_RECTANGLE_DOMAIN_INT:
    if((ir = (WlzIRect *) AlcMalloc (sizeof(WlzIRect))) == NULL){
      errNum = WLZ_ERR_MEM_ALLOC;
      break;
    }
    ir->type = type;
    ir->linkcount = 0;
    ir->freeptr = NULL;
    for(i=0; i < 4; i++){
      ir->irk[i] = getword(fp);
    }
    for(i=0; i < 4; i++){
      ir->irl[i] = getword(fp);
    }
    ir->rangle = getfloat(fp);
    break;

  case WLZ_RECTANGLE_DOMAIN_FLOAT:
    if((fr = (WlzFRect *) AlcMalloc (sizeof(WlzFRect))) == NULL){
      errNum = WLZ_ERR_MEM_ALLOC;
      break;
    }
    fr->type = type;
    fr->linkcount = 0;
    fr->freeptr = NULL;
    for(i=0; i < 4; i++){
      fr->frk[i] = getfloat(fp);
    }
    for(i=0; i < 4; i++){
      fr->frl[i] = getfloat(fp);
    }
    fr->rangle = getfloat(fp);
    ir = (WlzIRect *) fr;
    break;

  default:
    errNum = WLZ_ERR_OBJECT_TYPE;
    break;
  }

  if((errNum == WLZ_ERR_NONE) && (feof(fp) != 0) ){
    AlcFree( (void *) ir );
    ir = NULL;
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return( ir );
}

/*!
* \return	New histogram domain.
* \ingroup	WlzIO
* \brief	Reads a Woolz histogram domain.
* \param	fp			input file.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzHistogramDomain *WlzReadHistogramDomain(FILE *fp,
						  WlzErrorNum *dstErr)
{
  int		tI0,
  		numBins;
  double	origin,
  		binSize;
  int		*tIP0;
  double	*tDP0;
  WlzObjectType	type,
  		newType;
  WlzHistogramDomain *hist = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  type = (WlzObjectType )getc(fp);

  if(type == (WlzObjectType )EOF)
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else if(type == WLZ_NULL)
  {
    errNum = WLZ_ERR_EOO;
  }
  else
  {
    switch(type)
    {
      case WLZ_HISTOGRAMDOMAIN_OLD_INT:
      case WLZ_HISTOGRAMDOMAIN_OLD_FLOAT:
        (void )getword(fp);
	(void )getword(fp);
	(void )getword(fp);
	numBins = getword(fp);
	newType = (type == WLZ_HISTOGRAMDOMAIN_OLD_INT)?
		  WLZ_HISTOGRAMDOMAIN_INT:
		  WLZ_HISTOGRAMDOMAIN_FLOAT;
	if(feof(fp))
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
	else if((hist = WlzMakeHistogramDomain(newType, numBins,
					       NULL)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  hist->nBins = numBins;
	  hist->origin = 0.0;
	  hist->binSize = 1.0;
	  switch(type)
	  {
	    case WLZ_HISTOGRAMDOMAIN_OLD_INT:
	      tI0 = numBins;
	      tIP0 = hist->binValues.inp;
	      while(tI0-- > 0)
	      {
		*tIP0++ = getword(fp);
	      }
	      break;
	    case WLZ_HISTOGRAMDOMAIN_OLD_FLOAT:
	      tI0 = numBins;
	      tDP0 = hist->binValues.dbp;
	      while(tI0-- > 0)
	      {
		*tDP0++ = getfloat(fp);
	      }
	      break;
	    default:
	      break;
	  }
	}
	break;
      case WLZ_HISTOGRAMDOMAIN_INT: /* FALLTHROUGH */
      case WLZ_HISTOGRAMDOMAIN_FLOAT:
	numBins = getword(fp);
	origin = getdouble(fp);
	binSize = getdouble(fp);
	if(feof(fp))
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
	else if((hist = WlzMakeHistogramDomain(type, numBins,
					       NULL)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  switch(type)
	  {
	    case WLZ_HISTOGRAMDOMAIN_INT:
	      tI0 = numBins;
	      tIP0 = hist->binValues.inp;
	      while(tI0-- > 0)
	      {
		*tIP0++ = getword(fp);
	      }
	      break;
	    case WLZ_HISTOGRAMDOMAIN_FLOAT:
	      tI0 = numBins;
	      tDP0 = hist->binValues.dbp;
	      while(tI0-- > 0)
	      {
		*tDP0++ = getdouble(fp);
	      }
	      break;
	    default:
	      break;
	  }
	  hist->nBins = numBins;
	  hist->origin = origin;
	  hist->binSize = binSize;
	}
	break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
        break;
    }
    if((errNum == WLZ_ERR_NONE) && feof(fp))
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      if(hist)
      {
        WlzFreeHistogramDomain(hist);
	hist = NULL;
      }
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(hist);
}

/*!
* \return	New compound array object.
* \ingroup	WlzIO
* \brief	Reads a Woolz compund object.
* \param	fp			Input file.
* \param	type			Object type as read by WlzReadObj().
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzReadCompoundA(FILE			*fp,
				   WlzObjectType	type,
				   WlzErrorNum		*dstErr)
{
  WlzCompoundArray	*c=NULL;
  WlzObjectType		otype;
  int 			i, n;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  otype = (WlzObjectType) getc(fp);

  if( otype == (WlzObjectType) EOF ){
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else if( (type == WLZ_COMPOUND_ARR_1) && (otype == WLZ_NULL) ){
    errNum = WLZ_ERR_EOO;
  }
  else {
    n = getword(fp);
    if( feof(fp) != 0 ){
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }

  if((errNum == WLZ_ERR_NONE) &&
     ((c = WlzMakeCompoundArray(type, 1, n, NULL, otype, &errNum)) != NULL)){
    for(i=0; (i<n) && (errNum == WLZ_ERR_NONE); i++){
      c->o[i] = WlzAssignObject(WlzReadObj(fp, &errNum), NULL);
    }
    if( errNum == WLZ_ERR_NONE ){
      c->plist = WlzAssignPropertyList(WlzReadPropertyList(fp, NULL), NULL);
    }
    else {
      WlzFreeObj((WlzObject *) c);
      c = NULL;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return( (WlzObject *) c );
}

/*!
* \return	New affine transform.
* \ingroup	WlzIO
* \brief	Reads a woolz transform from the input file.
* \param	fp			Input file.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzAffineTransform *WlzReadAffineTransform(FILE *fp,
						  WlzErrorNum *dstErr)
{
  WlzTransformType	type;
  WlzAffineTransform	*trans=NULL;
  int			i, j;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  type = (WlzTransformType) getc( fp );

  if(type == (WlzTransformType )EOF){
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else if(type == (WlzTransformType )WLZ_NULL){
    errNum = WLZ_ERR_EOO;
  }
  else if((trans = WlzMakeAffineTransform(type, &errNum)) != NULL){

    /* set linkcount and freeptr */
    trans->linkcount = 0;
    trans->freeptr   = NULL;

    /* This code has been commented out rather than removed, just incase
     * there are any 2D affine transforms or afftine transobj's written to
     * a file somewhere using the old format.
     * dummyTx = getdouble(fp);
     * dummyTy = getdouble(fp);
     * dummyTz = getdouble(fp);
     * dummyScale = getdouble(fp);
     * dummyTheta = getdouble(fp);
     * dummyPhi = getdouble(fp);
     * dummyAlpha = getdouble(fp);
     * dummyPsi = getdouble(fp);
     * dummyXsi = getdouble(fp);
     * dummyInvert = getword(fp);
     */

    /* Read the matrix */
    for(i=0; i < 4; i++){
      for(j=0; j < 4; j++){
	trans->mat[i][j] = getdouble( fp );
      }
    }

    if( feof(fp) != 0 ){
      WlzFreeAffineTransform( trans );
      trans = NULL;
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return( trans );
}

/*!
* \return	New transform.
* \brief	Reads a Woolz FE warp transform.
* \param	fp			Input file.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzWarpTrans *WlzReadWarpTrans(FILE *fp, WlzErrorNum *dstErr)
{
  /* local variables */
  int	       	i, j;
  WlzObjectType type;
  WlzDVertex2 	*dptr;
  WlzTElement 	*eptr;
  WlzWarpTrans 	*obj=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  type = (WlzObjectType) getc(fp);

  if( type == (WlzObjectType) EOF ){
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else if( type == WLZ_NULL ){
    errNum = WLZ_ERR_EOO;
  }
  /* this is a repeat of the same char read by WlzReadObj
     and therefore is redundant, this can only be different if the file
     is corrupt */
  else if( type != WLZ_WARP_TRANS ){
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  /* make space for obj */
  else if( (obj = (WlzWarpTrans *) AlcMalloc(sizeof(WlzWarpTrans)))
	  == NULL ){
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else {
    obj->type = type;
    obj->linkcount = 0;

    obj->nelts = getword(fp);
    obj->nodes = getword(fp);
    obj->imdisp = getfloat(fp);
    obj->iterdisp = getfloat(fp);
    if( feof(fp) != 0 ){
      AlcFree( (void *) obj );
      obj = NULL;
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }

  /* read nodal coords */
  if( errNum == WLZ_ERR_NONE ){
    if( (obj->ncoords = (WlzDVertex2 *)
	 AlcMalloc(sizeof(WlzDVertex2) * obj->nodes)) == NULL ){
      AlcFree( (void *) obj );
      obj = NULL;
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else {
      dptr = obj->ncoords;
      for(i=0; i<obj->nodes; i++, dptr++){
	dptr->vtX = (double) getfloat(fp);
	dptr->vtY = (double) getfloat(fp);
      }
    }
  }

  /* read nodal displacements */
  if( errNum == WLZ_ERR_NONE ){
    if( (obj->displacements = (WlzDVertex2 *)
	 AlcMalloc(sizeof(WlzDVertex2) * obj->nodes)) == NULL ){
      AlcFree( (void *) obj->ncoords );
      AlcFree( (void *) obj );
      obj = NULL;
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else {
      dptr = obj->displacements ;
      for(i=0; i<obj->nodes; i++, dptr++){
	dptr->vtX = (double) getfloat(fp);
	dptr->vtY = (double) getfloat(fp);
      }
    }
  }

  /* read elements */
  if( errNum == WLZ_ERR_NONE ){
    if( (obj->eltlist = (WlzTElement *)
	 AlcMalloc(sizeof(WlzTElement) * obj->nelts)) == NULL ){
      AlcFree((void *)(obj->nodes));
      AlcFree((void *)(obj->ncoords));
      AlcFree((void *)obj);
      obj = NULL;
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else {
      eptr = obj->eltlist;
      for(i=0; i<obj->nelts; i++, eptr++){

	eptr->type = (WlzElementType )getc(fp);

	eptr->n = getword(fp);
	for(j=0; j<3; j++){
	  eptr->nodes[j] = getword(fp);
	}
	for (j=0; j<3; j++){
	  eptr->u[j] = getfloat(fp);
	}
	for (j=0; j<3; j++){
	  eptr->a[j] = getfloat(fp);
	}
      }
    }
  }

  /* check if EOF error has been set */
  if( (errNum == WLZ_ERR_NONE) && (feof(fp) != 0) ){
    AlcFree( (void *) obj->eltlist );
    AlcFree( (void *) obj->nodes );
    AlcFree( (void *) obj->ncoords );
    AlcFree((void *) obj);
    obj = NULL;
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(obj) ;
}

/*!
* \return	New match object.
* \ingroup	WlzIO
* \brief	Reads a FE warp match object.
* \param	fp			Input file.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzFMatchObj *WlzReadFMatchObj(FILE *fp,
				      WlzErrorNum *dstErr)
{
  int 			i, j;
  WlzFMatchPoint	*mptr;
  WlzFMatchObj		*obj=NULL;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  if( (obj = (WlzFMatchObj *) AlcMalloc(sizeof(WlzFMatchObj))) == NULL ){
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else {
    obj->type = WLZ_FMATCHOBJ;
    obj->linkcount = 0;

    obj->nopts = getword(fp);
    if( feof(fp) != 0 ){
      AlcFree((void *) obj);
      obj = NULL;
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    else {
      if( (obj->matchpts = (WlzFMatchPoint *)
	   AlcMalloc(sizeof(WlzFMatchPoint) * obj->nopts)) == NULL ){
	AlcFree((void *) obj);
	obj = NULL;
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else {
	mptr = obj->matchpts;
	for(i=0; i<obj->nopts; i++, mptr++){
	  mptr->type = (WlzMatchType )getword(fp);
	  mptr->node = getword(fp);
	  mptr->ptcoords.vtX = getfloat(fp);
	  mptr->ptcoords.vtY = getfloat(fp);

	  for(j=0; j<WLZ_MAX_NODAL_DEGREE; j++){
	    mptr->elements[j] = getword(fp);
	  }
	}
      }
    }
  }

  if( (errNum == WLZ_ERR_NONE) && (feof(fp) != 0) ){
    AlcFree((void *) obj->matchpts);
    AlcFree((void *) obj);
    obj = NULL;
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	New 3D warp transform.
* \ingroup	WlzIO
* \brief	Reads a 3D FE warp transform.
* \param	fp			Input file.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static Wlz3DWarpTrans *WlzRead3DWarpTrans(FILE *fp, WlzErrorNum *dstErr)
{
  /* local variables */
  int 			i, nplanes;
  Wlz3DWarpTrans	*obj=NULL;
  WlzFMatchObj		**intdoms;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  if( (obj = (Wlz3DWarpTrans *) AlcMalloc(sizeof(Wlz3DWarpTrans)))
     == NULL ){
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else {
    obj->type = WLZ_3D_WARP_TRANS;
    obj->linkcount = 0;
    obj->iteration = getword(fp);
    obj->currentplane = getword(fp);
    obj->maxdisp = getfloat(fp);
    if( feof(fp) != 0 ){
      AlcFree((void *) obj);
      obj = NULL;
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    if( (obj->pdom = WlzReadPlaneDomain(fp, &errNum)) == NULL ){
      AlcFree((void *) obj);
      obj = NULL;
    }
    else {
      obj->pdom->linkcount = 1;
      nplanes = obj->pdom->lastln - obj->pdom->line1 + 1;
      if( (obj->intptdoms = (WlzFMatchObj **)
	   AlcMalloc(sizeof(WlzFMatchObj *) * nplanes)) == NULL ){
	WlzFreePlaneDomain( obj->pdom );
	AlcFree((void *) obj);
	obj = NULL;
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else {
	intdoms = obj->intptdoms;
	for(i=obj->pdom->plane1; i<=obj->pdom->lastpl; i++, intdoms++){
	  if( (*intdoms = WlzReadFMatchObj(fp, &errNum)) ){
	    (*intdoms)->linkcount = 1;
	  }
	  else if( errNum != WLZ_ERR_EOO ){
	    WlzFree3DWarpTrans(obj);
	    obj = NULL;
	  }
	}
      }
    }
  }

  if( (errNum == WLZ_ERR_NONE) && (feof(fp) != 0) ){
    WlzFree3DWarpTrans(obj);
    obj = NULL;
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	New contour.
* \ingroup	WlzIO
* \brief	Reads a contour domain from the input file.
* \param	fP			Given file.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzContour *WlzReadContour(FILE *fP, WlzErrorNum *dstErr)
{
  WlzObjectType cType;
  WlzContour	*ctr = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;


  cType = (WlzObjectType )getc(fP);

  switch(cType)
  {
    case EOF:
      errNum = WLZ_ERR_READ_INCOMPLETE;
      break;
    case WLZ_NULL:
      errNum = WLZ_ERR_EOO;
      break;
    case WLZ_CONTOUR:
      if((ctr = WlzMakeContour(&errNum)) != NULL)
      {
        ctr->model = WlzAssignGMModel(WlzReadGMModel(fP, &errNum), NULL);
        if(errNum != WLZ_ERR_NONE)
	{
	  WlzFreeContour(ctr);
	  ctr = NULL;
	}
      }
      break;
    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(ctr);
}

/*!
* \return	New geometric model.
* \ingroup	WlzIO
* \brief	Read's a geometric model from the input file. For the
*		file format see WlzWriteGMModel().
* \param	fP			input file.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzGMModel *WlzReadGMModel(FILE *fP, WlzErrorNum *dstErr)
{
  int		tI0,
		idN,
		sCnt,
		encodeMtd,
  		nVertex = 0,
		nSimplex = 0,
		vgElmSz,
		vHTSz;
  int		bufI[3];
  void		*bufVG = NULL;
  WlzGMModelType mType;
  WlzGMModel	*model = NULL;
  WlzIVertex2	tIV2;
  WlzIVertex3	tIV3;
  WlzDVertex2	pos2D[2],
  		nrm2D[2];
  WlzDVertex3	pos3D[3],
  		nrm3D[3];
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  mType = (WlzGMModelType )getc(fP);

  switch(mType)
  {
    case EOF: /* FALLTHROUGH */
    case WLZ_NULL:
      errNum = WLZ_ERR_READ_INCOMPLETE;
      break;
    case WLZ_GMMOD_2I:
      vgElmSz = sizeof(WlzIVertex2);
      break;
    case WLZ_GMMOD_2D:
      vgElmSz = sizeof(WlzDVertex2);
      break;
    case WLZ_GMMOD_2N:
      vgElmSz = 2 * sizeof(WlzDVertex2);
      break;
    case WLZ_GMMOD_3I:
      vgElmSz = sizeof(WlzIVertex3);
      break;
    case WLZ_GMMOD_3D:
      vgElmSz = sizeof(WlzDVertex3);
      break;
    case WLZ_GMMOD_3N:
      vgElmSz = 2 * sizeof(WlzDVertex3);
      break;
    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Read the file encoding method, currently this is just a 'placeholder'
     * but in the future new encoding methods may be used, for example strips
     * of simplicies may be used to reduce the file size. */
    if((encodeMtd = getc(fP)) == EOF)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }

  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Read number of vertex's and simplicies (edges or loops). */
    if((errNum = WlzReadInt(fP, bufI, 2)) == WLZ_ERR_NONE)
    {
      if((nVertex = bufI[0]) < 0)
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else if(nVertex > 0)
      {
        nSimplex = bufI[1];
      }
    }
  }
  if((errNum == WLZ_ERR_NONE) && (nVertex > 0))
  {
    /* Choose an appropriate vertex (matching) hash table size and make a new
     * GM. */
    vHTSz = ((tI0 = nVertex) < 1024)? 1024: tI0;
    model = WlzGMModelNew(mType, 0, vHTSz, &errNum);
  }
  if((errNum == WLZ_ERR_NONE) && (nVertex > 0))
  {
    /* Create a vertex buffer. */
    if((bufVG = AlcMalloc(vgElmSz * nVertex)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (nVertex > 0))
  {
    /* Read the vertex's into the buffer. */
    switch(mType)
    {
      case WLZ_GMMOD_2I:
	errNum = WlzReadVertex2I(fP, (WlzIVertex2 *)bufVG, nVertex);
	break;
      case WLZ_GMMOD_2D:
	errNum = WlzReadVertex2D(fP, (WlzDVertex2 *)bufVG, nVertex);
	break;
      case WLZ_GMMOD_2N:
	errNum = WlzReadVertex2D(fP, (WlzDVertex2 *)bufVG, 2 * nVertex);
	break;
      case WLZ_GMMOD_3I:
	errNum = WlzReadVertex3I(fP, (WlzIVertex3 *)bufVG, nVertex);
	break;
      case WLZ_GMMOD_3D:
	errNum = WlzReadVertex3D(fP, (WlzDVertex3 *)bufVG, nVertex);
	break;
      case WLZ_GMMOD_3N:
	errNum = WlzReadVertex3D(fP, (WlzDVertex3 *)bufVG, 2 * nVertex);
	break;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (nVertex > 0))
  {
    /* Read the vertex indicies and build the model. */
    sCnt = 0;
    while((errNum == WLZ_ERR_NONE) && (sCnt++ < nSimplex))
    {
      switch(mType)
      {
	case WLZ_GMMOD_2I:
	  if((errNum = WlzReadInt(fP, bufI, 2)) == WLZ_ERR_NONE)
	  {
	    for(idN = 0; idN < 2; ++idN)
	    {
	      tIV2 = *(((WlzIVertex2 *)bufVG) + bufI[idN]);
	      pos2D[idN].vtX = tIV2.vtX;
	      pos2D[idN].vtY = tIV2.vtY;
	    }
	    errNum = WlzGMModelConstructSimplex2D(model, pos2D);
	  }
	  break;
	case WLZ_GMMOD_2D:
	  if((errNum = WlzReadInt(fP, bufI, 2)) == WLZ_ERR_NONE)
	  {
	    for(idN = 0; idN < 2; ++idN)
	    {
	      pos2D[idN] = *(((WlzDVertex2 *)bufVG) + bufI[idN]);
	    }
	    errNum = WlzGMModelConstructSimplex2D(model, pos2D);
	  }
	  break;
	case WLZ_GMMOD_2N:
	  if((errNum = WlzReadInt(fP, bufI, 2)) == WLZ_ERR_NONE)
	  {
	    for(idN = 0; idN < 2; ++idN)
	    {
	      pos2D[idN] = *(((WlzDVertex2 *)bufVG) + (2 * bufI[idN]) + 0);
	      nrm2D[idN] = *(((WlzDVertex2 *)bufVG) + (2 * bufI[idN]) + 1);
	    }
	    errNum = WlzGMModelConstructSimplex2N(model, pos2D, nrm2D);
	  }
	  break;
	case WLZ_GMMOD_3I:
	  if((errNum = WlzReadInt(fP, bufI, 3)) == WLZ_ERR_NONE)
	  {
	    for(idN = 0; idN < 3; ++idN)
	    {
	      tIV3 = *(((WlzIVertex3 *)bufVG) + bufI[idN]);
	      pos3D[idN].vtX = tIV3.vtX;
	      pos3D[idN].vtY = tIV3.vtY;
	      pos3D[idN].vtY = tIV3.vtZ;
	    }
	    errNum = WlzGMModelConstructSimplex3D(model, pos3D);
	  }
	  break;
	case WLZ_GMMOD_3D:
	  if((errNum = WlzReadInt(fP, bufI, 3)) == WLZ_ERR_NONE)
	  {
	    for(idN = 0; idN < 3; ++idN)
	    {
	      pos3D[idN] = *(((WlzDVertex3 *)bufVG) + bufI[idN]);
	    }
	    errNum = WlzGMModelConstructSimplex3D(model, pos3D);
	  }
	  break;
	case WLZ_GMMOD_3N:
	  if((errNum = WlzReadInt(fP, bufI, 3)) == WLZ_ERR_NONE)
	  {
	    for(idN = 0; idN < 3; ++idN)
	    {
	      pos3D[idN] = *(((WlzDVertex3 *)bufVG) + (2 * bufI[idN]) + 0);
	      nrm3D[idN] = *(((WlzDVertex3 *)bufVG) + (2 * bufI[idN]) + 1);
	    }
	    errNum = WlzGMModelConstructSimplex3N(model, pos3D, nrm3D);
	  }
	  break;
      }
    }
  }
  if(bufVG)
  {
    AlcFree(bufVG);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(model);
}


/*! 
* \ingroup      WlzIO
* \brief        reads a woolz 3D MeshTransform
* \author	J. Rao.
*
* \return       Woolz mesh transform read from the input stream.
* \param    fp	Input file pointer.
* \param    dstErr	Error return.
* \par      Source:
*                WlzReadObj.c
*/
WlzMeshTransform3D *WlzReadMeshTransform3D(FILE *fp,
				      WlzErrorNum *dstErr)
{
  /* local variables */
  int		i,
  		j;
  WlzMeshNode3D	*dptr;
  WlzMeshElem3D	*eptr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  /* WlzObjectType type; */
  WlzMeshTransform3D   	*obj=NULL;
  /*
  type = (WlzObjectType) getc(fp);
  if( type == (WlzObjectType) EOF ){
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else if( type == WLZ_NULL ){
    errNum = WLZ_ERR_EOO;
  }
  */
  /* this is a repeat of the same char read by WlzReadObj
     and therefore is redundant, this can only be different if the file
     is corrupt */
     /*
  else if( type != WLZ_WARP_TRANS ){
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  */
  /* make space for obj */
  if( (obj = (WlzMeshTransform3D *) AlcMalloc(sizeof(WlzMeshTransform3D)))
	  == NULL )
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    /*
    obj->type = type;
    obj->linkcount = 0;
    */
    obj->nElem = getword(fp);
    obj->nNodes = getword(fp);
    if( feof(fp) != 0 )
    {
      AlcFree( (void *) obj );
      obj = NULL;
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }

  /* read  nodal position and displacement */
  if( errNum == WLZ_ERR_NONE )
  {
    if( (obj->nodes = (WlzMeshNode3D *)
                      AlcMalloc(obj->nNodes * sizeof(WlzMeshNode3D))) == NULL)
    {
      AlcFree( (void *) obj->nodes );
      AlcFree( (void *) obj );
      obj = NULL;
      errNum = WLZ_ERR_MEM_ALLOC;
      printf("Failed to allocate memory!");
      exit(1);
    }
    dptr = obj->nodes;
    for(i = 0; (i < obj->nNodes) && (errNum == WLZ_ERR_NONE); i++, dptr++)
    {  /* read position */
        dptr->position.vtX = (double) getfloat(fp);
        dptr->position.vtY = (double) getfloat(fp);
        dptr->position.vtZ = (double) getfloat(fp);
       /* read displacement */
        dptr->displacement.vtX =(double) getfloat(fp);
        dptr->displacement.vtY =(double) getfloat(fp);
        dptr->displacement.vtZ =(double) getfloat(fp);
    }
  }

  /* read elements */
  if( errNum == WLZ_ERR_NONE ){
    if( (obj->elements = (WlzMeshElem3D *)AlcMalloc(  obj->nElem *
                           sizeof(WlzMeshElem3D))) == NULL )
     {
      AlcFree( (void *) obj->nodes );
      AlcFree( (void *) obj->elements );
      AlcFree( (void *) obj );
      obj = NULL;
      errNum = WLZ_ERR_MEM_ALLOC;
     }
     else
     {
     /* read elements */
     eptr = obj->elements;
       for(i = 0; (i < obj->nElem) && (errNum == WLZ_ERR_NONE); i++, eptr++)
       {/*
        if(!putc((unsigned int )eptr->type, fp)  )
         {
          errNum = WLZ_ERR_WRITE_INCOMPLETE;
         } */
         /* get the index of this element */
         eptr->idx = getword(fp);
         /* get nodes indeces */
         for(j = 0; (j < 4) && (errNum == WLZ_ERR_NONE); j++)
         {
           eptr->nodes[j] = getword(fp);
         }
         /* output its neighbours */
         /*
         for(j = 0; (j < 4) && (errNum == WLZ_ERR_NONE); j++)
         {
	  eptr->neighbours[j] = getword(fp);
         }
         */
       }
     }
   }
  /* check if EOF error has been set */
  if( (errNum == WLZ_ERR_NONE) && (feof(fp) != 0) ){
    AlcFree( (void *) obj->elements );
    AlcFree( (void *) obj->nodes );
    AlcFree((void *) obj);
    obj = NULL;
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(obj) ;
}

/*!
* \return	Mesh Transform.
* \ingroup	WlzIO
* \brief	Reads a mesh transform from the input file.
* \param	fP			Given file.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzMeshTransform *WlzReadMeshTransform2D(FILE *fp,
				      WlzErrorNum *dstErr)
{

  /* local variables */
  int		i,
  		j;
  WlzMeshNode	*dptr;
  WlzMeshElem	*eptr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzMeshTransform   	*obj=NULL;

  /* make space for obj */
  if( (obj = (WlzMeshTransform *) AlcMalloc(sizeof(WlzMeshTransform)))
	  == NULL )
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    obj->type = WLZ_TRANSFORM_2D_MESH;
    obj->linkcount = 0;

    obj->nElem = getword(fp);
    obj->nNodes = getword(fp);
    if( feof(fp) != 0 )
    {
      AlcFree( (void *) obj );
      obj = NULL;
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }

  /* read  nodal position and displacement */
  if( errNum == WLZ_ERR_NONE )
  {
    if( (obj->nodes = (WlzMeshNode *)AlcCalloc(  obj->nNodes, sizeof(WlzMeshNode))) == NULL )
    {
      AlcFree( (void *) obj->nodes );
      AlcFree( (void *) obj );
      obj = NULL;
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    if (errNum == WLZ_ERR_NONE)
    {
      dptr = obj->nodes;
      for(i = 0; (i < obj->nNodes) && (errNum == WLZ_ERR_NONE); i++, dptr++)
      {  /* read position */
        dptr->position.vtX = (double) getfloat(fp);
        dptr->position.vtY = (double) getfloat(fp);
       /* read displacement */
        dptr->displacement.vtX =(double) getfloat(fp);
        dptr->displacement.vtY =(double) getfloat(fp);
      }
    }

    /* read elements */
    if( errNum == WLZ_ERR_NONE ){
      if( (obj->elements = (WlzMeshElem *)AlcCalloc(  obj->nElem,
                           sizeof(WlzMeshElem))) == NULL )
      {
        AlcFree( (void *) obj->nodes );
        AlcFree( (void *) obj->elements );
        AlcFree( (void *) obj );
        obj = NULL;
        errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
        /* read elements */
        eptr = obj->elements;
        for(i = 0; (i < obj->nElem) && (errNum == WLZ_ERR_NONE); i++, eptr++)
        {
          /* get the index of this element */
          eptr->idx = getword(fp);
	  /* get neighbour flags */
	  eptr->flags = (unsigned int) getword(fp);
          /* get nodes indeces */
          for(j = 0; (j < 3) && (errNum == WLZ_ERR_NONE); j++)
          {
            eptr->nodes[j] = getword(fp);
          }
          /* output its neighbours */

          for(j = 0; (j < 3) && (errNum == WLZ_ERR_NONE); j++)
          {
	   eptr->neighbours[j] = getword(fp);
          }
        }
      }
    }
  }

  /* check if EOF error has been set */
  if( (errNum == WLZ_ERR_NONE) && (feof(fp) != 0) ){
    AlcFree( (void *) obj->elements );
    AlcFree( (void *) obj->nodes );
    AlcFree((void *) obj);
    obj = NULL;
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(obj) ;
}

/*!
* \return	New 2D constrained mesh.
* \ingroup	WlzIO
* \brief	reads a 2D constrained mesh using the given file pointer.
* \param	fp			Given file pointer.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzCMesh2D *WlzReadCMesh2D(FILE *fp, WlzErrorNum *dstErr)
{
  int		idE,
  		idN,
		flags,
		nElm,
		nNod;
  WlzDVertex2	pos;
  WlzCMeshElm2D	*elm;
  WlzCMesh2D	*mesh = NULL;
  int		idx[3];
  WlzCMeshNod2D	*nod[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  nNod = getword(fp);
  nElm = getword(fp);
  if(feof(fp) != 0)
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else if((nNod < 0) || (nElm < 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
#ifdef WLZ_DEBUG_READOBJ
  if(errNum == WLZ_ERR_NONE)
  {
    (void )fprintf(stderr,
    		   "WlzReadCMesh2D() "
		   "%d %d\n",
		   nNod, nElm);
  }
#endif /* WLZ_DEBUG_READOBJ */
  if(errNum == WLZ_ERR_NONE)
  {
    mesh = WlzCMeshNew2D(&errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((AlcVectorExtendAndGet(mesh->res.nod.vec, nNod) == NULL) ||
       (AlcVectorExtendAndGet(mesh->res.elm.vec, nElm) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      errNum = WlzCMeshReassignGridCells2D(mesh, nNod);
    }
  }
  /* Read nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    for(idN = 0; idN < nNod; ++idN)
    {
      flags = getword(fp);
      pos.vtX = getdouble(fp);
      pos.vtY = getdouble(fp);
#ifdef WLZ_DEBUG_READOBJ
      (void )fprintf(stderr,
                     "WlzReadCMesh2D() "
		     "% 8d 0x%08x % 8g % 8g\n",
		     idN, flags, pos.vtX, pos.vtY);
#endif /* WLZ_DEBUG_READOBJ */
      if(feof(fp) != 0)
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }
      if(errNum == WLZ_ERR_NONE)
      {
        nod[0] = WlzCMeshAllocNod2D(mesh);
	nod[0]->pos = pos;
	nod[0]->flags = flags;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Update the bounding box. */
    WlzCMeshUpdateBBox2D(mesh);
    /* Compute a new bucket grid and reassign nodes to it. */
    errNum = WlzCMeshReassignGridCells2D(mesh, nNod);
  }
  /* Read elements and create them within the mesh using the already
   * created nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    for(idE = 0; idE < nElm; ++idE)
    {
      flags = getword(fp);
      idx[0] = getword(fp);
      idx[1] = getword(fp);
      idx[2] = getword(fp);
      if(feof(fp) != 0)
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }
      if(errNum == WLZ_ERR_NONE)
      {
#ifdef WLZ_DEBUG_READOBJ
	(void )fprintf(stderr,
		       "WlzReadCMesh2D() "
		       "% 8d 0x%08x % 8d % 8d % 8d\n",
		       idE, flags, idx[0], idx[1], idx[2]);
#endif /* WLZ_DEBUG_READOBJ */
	for(idN = 0; idN < 3; ++idN)
	{
	  if((idx[idN] < 0) || (idx[idN] >= nNod))
	  {
	    errNum = WLZ_ERR_DOMAIN_DATA;
	    break;
	  }
	  nod[idN] = (WlzCMeshNod2D *)
	             AlcVectorItemGet(mesh->res.nod.vec, idx[idN]);
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
        elm = WlzCMeshNewElm2D(mesh,
	                       nod[0], nod[1], nod[2], &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        elm->flags = flags;
      }
      else
      {
        break;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Compute maximum edge length. */
    WlzCMeshUpdateMaxSqEdgLen2D(mesh);
  }
#ifdef WLZ_DEBUG_READOBJ
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzCMeshVerify2D(mesh, NULL, 1, stderr);
  }
#endif /* WLZ_DEBUG_READOBJ */
  /* Clean up if errors. */
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzCMeshFree2D(mesh);
    mesh = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mesh);
}

/*!
* \return	New 2D5 constrained mesh.
* \ingroup	WlzIO
* \brief	reads a 2D5 constrained mesh using the given file pointer.
* \param	fp			Given file pointer.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzCMesh2D5 *WlzReadCMesh2D5(FILE *fp, WlzErrorNum *dstErr)
{
  int		idE,
  		idN,
		flags,
		nElm,
		nNod;
  WlzDVertex3	pos;
  WlzCMeshElm2D5 *elm;
  WlzCMesh2D5	*mesh = NULL;
  int		idx[3];
  WlzCMeshNod2D5 *nod[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  nNod = getword(fp);
  nElm = getword(fp);
  if(feof(fp) != 0)
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else if((nNod < 0) || (nElm < 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
#ifdef WLZ_DEBUG_READOBJ
  if(errNum == WLZ_ERR_NONE)
  {
    (void )fprintf(stderr,
    		   "WlzReadCMesh2D5() "
		   "%d %d\n",
		   nNod, nElm);
  }
#endif /* WLZ_DEBUG_READOBJ */
  if(errNum == WLZ_ERR_NONE)
  {
    mesh = WlzCMeshNew2D5(&errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((AlcVectorExtendAndGet(mesh->res.nod.vec, nNod) == NULL) ||
       (AlcVectorExtendAndGet(mesh->res.elm.vec, nElm) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      errNum = WlzCMeshReassignGridCells2D5(mesh, nNod);
    }
  }
  /* Read nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    for(idN = 0; idN < nNod; ++idN)
    {
      flags = getword(fp);
      pos.vtX = getdouble(fp);
      pos.vtY = getdouble(fp);
      pos.vtZ = getdouble(fp);
#ifdef WLZ_DEBUG_READOBJ
      (void )fprintf(stderr,
                     "WlzReadCMesh2D5() "
		     "% 8d 0x%08x % 8g % 8g % 8g\n",
		     idN, flags, pos.vtX, pos.vtY, pos.vtZ);
#endif /* WLZ_DEBUG_READOBJ */
      if(feof(fp) != 0)
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }
      if(errNum == WLZ_ERR_NONE)
      {
        nod[0] = WlzCMeshAllocNod2D5(mesh);
	nod[0]->pos = pos;
	nod[0]->flags = flags;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Update the bounding box. */
    WlzCMeshUpdateBBox2D5(mesh);
    /* Compute a new bucket grid and reassign nodes to it. */
    errNum = WlzCMeshReassignGridCells2D5(mesh, nNod);
  }
  /* Read elements and create them within the mesh using the already
   * created nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    for(idE = 0; idE < nElm; ++idE)
    {
      flags = getword(fp);
      idx[0] = getword(fp);
      idx[1] = getword(fp);
      idx[2] = getword(fp);
      if(feof(fp) != 0)
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }
      if(errNum == WLZ_ERR_NONE)
      {
#ifdef WLZ_DEBUG_READOBJ
	(void )fprintf(stderr,
		       "WlzReadCMesh2D5() "
		       "% 8d 0x%08x % 8d % 8d % 8d\n",
		       idE, flags, idx[0], idx[1], idx[2]);
#endif /* WLZ_DEBUG_READOBJ */
	for(idN = 0; idN < 3; ++idN)
	{
	  if((idx[idN] < 0) || (idx[idN] >= nNod))
	  {
	    errNum = WLZ_ERR_DOMAIN_DATA;
	    break;
	  }
	  nod[idN] = (WlzCMeshNod2D5 *)
	             AlcVectorItemGet(mesh->res.nod.vec, idx[idN]);
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
        elm = WlzCMeshNewElm2D5(mesh,
	                       nod[0], nod[1], nod[2], &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        elm->flags = flags;
      }
      else
      {
        break;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Compute maximum edge length. */
    WlzCMeshUpdateMaxSqEdgLen2D5(mesh);
  }
#ifdef WLZ_DEBUG_READOBJ
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzCMeshVerify2D5(mesh, NULL, 1, stderr);
  }
#endif /* WLZ_DEBUG_READOBJ */
  /* Clean up if errors. */
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzCMeshFree2D5(mesh);
    mesh = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mesh);
}

/*!
* \return	New 3D constrained mesh.
* \ingroup	WlzIO
* \brief	reads a 3D constrained mesh using the given file pointer.
* \param	fp			Given file pointer.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzCMesh3D *WlzReadCMesh3D(FILE *fp, WlzErrorNum *dstErr)
{
  int		idE,
  		idN,
		flags,
		nElm,
		nNod;
  WlzDVertex3	pos;
  WlzDVertex3	*posBuf = NULL;
  WlzCMeshElm3D	*elm;
  WlzCMesh3D	*mesh = NULL;
  int		idx[4];
  WlzCMeshNod3D	*nod[4];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  nNod = getword(fp);
  nElm = getword(fp);
  if(feof(fp) != 0)
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else if((nNod < 0) || (nElm < 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
#ifdef WLZ_DEBUG_READOBJ
  if(errNum == WLZ_ERR_NONE)
  {
    (void )fprintf(stderr,
    		   "WlzReadCMesh3D() "
		   "%d %d\n",
		   nNod, nElm);
  }
#endif /* WLZ_DEBUG_READOBJ */
  if(errNum == WLZ_ERR_NONE)
  {
    mesh = WlzCMeshNew3D(&errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(((posBuf = AlcMalloc(sizeof(WlzDVertex3) * nNod)) == NULL) ||
       (AlcVectorExtendAndGet(mesh->res.nod.vec, nNod) == NULL) ||
       (AlcVectorExtendAndGet(mesh->res.elm.vec, nElm) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Read node positionss into a temporary buffer and compute their bounding
   * box. */
  if(errNum == WLZ_ERR_NONE)
  {
    for(idN = 0; idN < nNod; ++idN)
    {
      flags = getword(fp);
      pos.vtX = getdouble(fp);
      pos.vtY = getdouble(fp);
      pos.vtZ = getdouble(fp);
#ifdef WLZ_DEBUG_READOBJ
      (void )fprintf(stderr,
                     "WlzReadCMesh3D() "
		     "% 8d 0x%08x % 8g % 8g % 8g\n",
		     idN, flags, pos.vtX, pos.vtY, pos.vtZ);
#endif /* WLZ_DEBUG_READOBJ */
      if(feof(fp) != 0)
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }
      posBuf[idN] = pos;
      if(idN == 0)
      {
	mesh->bBox.xMin = mesh->bBox.xMax = pos.vtX;
	mesh->bBox.yMin = mesh->bBox.yMax = pos.vtY;
	mesh->bBox.zMin = mesh->bBox.zMax = pos.vtZ;
      }
      else
      {
        if(pos.vtX < mesh->bBox.xMin)
	{
	  mesh->bBox.xMin = pos.vtX;
	}
	else if(pos.vtX > mesh->bBox.xMax)
	{
	  mesh->bBox.xMax = pos.vtX;
	}
        if(pos.vtY < mesh->bBox.yMin)
	{
	  mesh->bBox.yMin = pos.vtY;
	}
	else if(pos.vtY > mesh->bBox.yMax)
	{
	  mesh->bBox.yMax = pos.vtY;
	}
        if(pos.vtZ < mesh->bBox.zMin)
	{
	  mesh->bBox.zMin = pos.vtZ;
	}
	else if(pos.vtZ > mesh->bBox.zMax)
	{
	  mesh->bBox.zMax = pos.vtZ;
	}
      }
    }
  }
  /* Setup the bucket grid using the bounding box of the node positions. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzCMeshReassignGridCells3D(mesh, nNod);
  }
  /* Create all nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    for(idN = 0; idN < nNod; ++idN)
    {
      pos = posBuf[idN];
      nod[0] = WlzCMeshNewNod3D(mesh, pos, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	nod[0]->pos = pos;
	nod[0]->flags = flags;
      }
    }
  }
  AlcFree(posBuf);
  /* Read elements and create them within the mesh using the already
   * created nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    for(idE = 0; idE < nElm; ++idE)
    {
      flags = getword(fp);
      idx[0] = getword(fp);
      idx[1] = getword(fp);
      idx[2] = getword(fp);
      idx[3] = getword(fp);
      if(feof(fp) != 0)
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }
      if(errNum == WLZ_ERR_NONE)
      {
#ifdef WLZ_DEBUG_READOBJ
	(void )fprintf(stderr,
		       "WlzReadCMesh3D() "
		       "% 8d 0x%08x % 8d % 8d % 8d % 8d\n",
		       idE, flags, idx[0], idx[1], idx[2], idx[3]);
#endif /* WLZ_DEBUG_READOBJ */
	for(idN = 0; idN < 4; ++idN)
	{
	  if((idx[idN] < 0) || (idx[idN] >= nNod))
	  {
	    errNum = WLZ_ERR_DOMAIN_DATA;
	    break;
	  }
	  nod[idN] = (WlzCMeshNod3D *)
	             AlcVectorItemGet(mesh->res.nod.vec, idx[idN]);
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
        elm = WlzCMeshNewElm3D(mesh,
	                       nod[0], nod[1], nod[2], nod[3], &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        elm->flags = flags;
      }
      else
      {
        break;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Compute maximum edge length. */
    WlzCMeshUpdateMaxSqEdgLen3D(mesh);
  }
  /* Clean up if errors. */
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzCMeshFree3D(mesh);
    mesh = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mesh);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Reads an indexed value table.
* 		Watch out that this function assumes that the given object
* 		has NULL values, which is true in it's current use.
* \param	fP			File pointer.
* \param	obj			Object that the indexed value table
* 					will be attached to.
*/
static WlzErrorNum WlzReadIndexedvValues(FILE *fP, WlzObject *obj)
{
  int		idX,
  		idV,
		rank,
  		nValues,
		vCount;
  WlzObjectType	type;
  WlzGreyType	vType;
  WlzValueAttach vAttach;
  int		*dim = NULL;
  WlzGreyP	gP;
  WlzValues 	values;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  values.core = NULL;
  type = (WlzObjectType )getc(fP);
  switch(type)
  {
    case WLZ_INDEXED_VALUES:
      break;
    case WLZ_NULL:
      errNum = WLZ_ERR_EOO;
      break;
    default:
      errNum = WLZ_ERR_READ_INCOMPLETE;
      break;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rank = getword(fP);
    if(rank < 0)
    {
      errNum = WLZ_ERR_VALUES_DATA;
    }
    else if((dim = (int *)AlcMalloc((sizeof(int) * rank))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(idX = 0; idX < rank; ++idX)
    {
      dim[idX] = getword(fP);
    }
    vType = (WlzGreyType )getc(fP);
    vAttach = (WlzValueAttach )getc(fP);
    nValues = getword(fP);
    if(feof(fP))
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    values.x = WlzMakeIndexedValues(obj, rank, dim, vType, vAttach, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(AlcVectorExtendAndGet(values.x->values, nValues) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    vCount = 1;
    if(rank > 0)
    {
      for(idX = 0; idX < rank; ++idX)
      {
        vCount *= dim[idX];
      }
    }
    for(idX = 0; idX < nValues; ++idX)
    {
      gP.v = WlzIndexedValueGet(values.x, idX);
      for(idV = 0; idV < vCount; ++idV)
      {
	switch(vType)
	{
	  case WLZ_GREY_INT:
	    gP.inp[idV] = getword(fP);
	    break;
	  case WLZ_GREY_SHORT:
	    gP.shp[idV] = getshort(fP);
	    break;
	  case WLZ_GREY_UBYTE:
	    gP.ubp[idV] = getc(fP);
	    break;
	  case WLZ_GREY_FLOAT:
	    gP.flp[idV] = getfloat(fP);
	    break;
	  case WLZ_GREY_DOUBLE:
	    gP.dbp[idV] = getdouble(fP);
	    break;
	  case WLZ_GREY_RGBA:
	    gP.rgbp[idV] = getword(fP);
	    break;
	  default:
	    errNum = WLZ_ERR_GREY_TYPE;
	    break;
	}
	if(feof(fP))
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	  break;
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    obj->values = WlzAssignValues(values, NULL);
  }
  else
  {
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;                        /* No values is allowed. */
    }
    else
    {
      (void )WlzFreeIndexedValues(values.x);
      values.x = NULL;
    }
  }
  return(errNum);
}

#ifdef WLZ_OLD_CMESH_TRANS_SUPPORT
/*!
* \return	New Woolz constrained mesh object.
* \ingroup	WlzIO
* \brief	Reads an old format constrained mesh transform into the
* 		new constrained mesh data structures. Only the object
* 		type will have already been read from the input file.
* \param	fP			File pointer for input.
* \param	dstErr			destination error pointer may be NULL.
*/
static WlzObject *WlzReadOldCMeshTransform(FILE *fP, WlzErrorNum *dstErr)
{
  WlzTransformType tType;
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  tType = (WlzTransformType )getword(fP);
  if(feof(fP) != 0)
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else
  {
    switch(tType)
    {
      case WLZ_TRANSFORM_2D_CMESH:
	obj = WlzReadOldCMeshTransform2D(fP, &errNum);
	break;
      case WLZ_TRANSFORM_3D_CMESH:
	obj = WlzReadOldCMeshTransform3D(fP, &errNum);
	break;
      default:
        errNum = WLZ_ERR_TRANSFORM_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	New Woolz constrained mesh object.
* \ingroup	WlzIO
* \brief	Reads an old format 2D constrained mesh transform into the
* 		new constrained mesh data structures. Only the object
* 		and transform types will have already been read from the
* 		input file.
* \param	fP			File pointer for input.
* \param	dstErr			destination error pointer may be NULL.
*/
static WlzObject *WlzReadOldCMeshTransform2D(FILE *fP, WlzErrorNum *dstErr)
{
  int		idX,
  		dim,
  		nDsp;
  WlzDomain	dom;
  WlzValues	val;
  WlzGreyP	gP;
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  val.core = NULL;
  dom.cm2 = WlzReadCMesh2D(fP, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    if((nDsp = dom.cm2->res.nod.maxEnt) < 0)
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    obj = WlzMakeMain(WLZ_CMESH_2D, dom, val, NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dim = 2;
    val.x = WlzMakeIndexedValues(obj, 1, &dim,
                                 WLZ_GREY_DOUBLE, WLZ_VALUE_ATTACH_NOD,
				 &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(AlcVectorExtendAndGet(val.x->values, nDsp) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(idX = 0; idX < nDsp; ++idX)
    {
      gP.v = WlzIndexedValueGet(val.x, idX);
      gP.dbp[0] = getdouble(fP);
      gP.dbp[1] = getdouble(fP);
      if(feof(fP))
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    obj->values = WlzAssignValues(val, NULL);
  }
  else
  {
    if(obj)
    {
      (void )WlzFreeObj(obj);
      obj = NULL;
    }
    else
    {
      (void )WlzFreeDomain(dom);
      if(val.core)
      {
        (void )WlzFreeIndexedValues(val.x);
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	New Woolz constrained mesh object.
* \ingroup	WlzIO
* \brief	Reads an old format 3D constrained mesh transform into the
* 		new constrained mesh data structures. Only the object
* 		and transform types will have already been read from the
* 		input file.
* \param	fP			File pointer for input.
* \param	dstErr			destination error pointer may be NULL.
*/
static WlzObject *WlzReadOldCMeshTransform3D(FILE *fP, WlzErrorNum *dstErr)
{
  int		idX,
  		dim,
  		nDsp;
  WlzDomain	dom;
  WlzValues	val;
  WlzGreyP	gP;
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  val.core = NULL;
  dom.cm3 = WlzReadCMesh3D(fP, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    if((nDsp = dom.cm3->res.nod.maxEnt) < 0)
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    obj = WlzMakeMain(WLZ_CMESH_3D, dom, val, NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dim = 3;
    val.x = WlzMakeIndexedValues(obj, 1, &dim,
                                 WLZ_GREY_DOUBLE, WLZ_VALUE_ATTACH_NOD,
				 &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(AlcVectorExtendAndGet(val.x->values, nDsp) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(idX = 0; idX < nDsp; ++idX)
    {
      gP.v = WlzIndexedValueGet(val.x, idX);
      gP.dbp[0] = getdouble(fP);
      gP.dbp[1] = getdouble(fP);
      gP.dbp[2] = getdouble(fP);
      if(feof(fP))
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    obj->values = WlzAssignValues(val, NULL);
  }
  else
  {
    if(obj)
    {
      (void )WlzFreeObj(obj);
      obj = NULL;
    }
    else
    {
      (void )WlzFreeDomain(dom);
      if(val.core)
      {
        (void )WlzFreeIndexedValues(val.x);
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}
#endif
