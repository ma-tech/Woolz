#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzWriteObj_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzWriteObj.c
* \author       Richard Baldock, Bill Hill
* \date         September 2005
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
* \brief	Functions for writing Woolz objects.
* \ingroup	WlzIO
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <Wlz.h>

#if defined(_WIN32) && !defined(__x86)
#define __x86
#endif


static WlzErrorNum	WlzWriteIntervalDomain(
			  FILE *fp,
			  WlzIntervalDomain *itvl);
static WlzErrorNum   	WlzWritePlaneDomain(
			  FILE *fp,
		    	  WlzPlaneDomain *planedm);
static WlzErrorNum	WlzWritePropertyList(
			  FILE *fp,
			  WlzPropertyList *pList);
static WlzErrorNum	WlzWriteProperty(
			  FILE *fp,
			  WlzProperty property);
static WlzErrorNum	WlzWriteValueTable(
			  FILE	*fp,
			  WlzObject *obj);
static WlzErrorNum	WlzWriteVoxelValueTable(
			  FILE *fp,
			  WlzObject *obj);
static WlzErrorNum	WlzWritePolygon(
			  FILE *fp,
			  WlzPolygonDomain *poly);
static WlzErrorNum	WlzWriteBoundList(
			  FILE *fp,
			  WlzBoundList *blist);
static WlzErrorNum      WlzWriteConvexHullValues(
                          FILE *fp,
			  WlzConvHullValues *cnvhull);
static WlzErrorNum	WlzWriteRect(
			  FILE *fp,
			  WlzIRect *rdom);
static WlzErrorNum	WlzWriteHistogramDomain(
			  FILE *fp,
			  WlzHistogramDomain *hist);
static WlzErrorNum	WlzWriteCompoundA(
			  FILE *fp,
			  WlzCompoundArray *c);
static WlzErrorNum	WlzWriteAffineTransform(
			  FILE *fp,
			  WlzAffineTransform *trans);
static WlzErrorNum	WlzWriteWarpTrans(
			  FILE *fp,
			  WlzWarpTrans *obj);
static WlzErrorNum	WlzWriteFMatchObj(
			  FILE *fp,
			  WlzFMatchObj *obj);
static WlzErrorNum	WlzWrite3DWarpTrans(
			  FILE *fp,
			  Wlz3DWarpTrans *obj);
static WlzErrorNum 	WlzWriteContour(
			  FILE *fP,
			  WlzContour *ctr);
static WlzErrorNum 	WlzWriteGMModel(
			  FILE *fP,
			  WlzGMModel *model);
static WlzErrorNum	WlzWriteInt(
			  FILE *fP,
			  int *iP,
			  int nI);
static WlzErrorNum 	WlzWriteVertex2I(
			  FILE *fP,
			  WlzIVertex2 *vP,
			  int nV);
static WlzErrorNum 	WlzWriteVertex2D(
			  FILE *fP,
			  WlzDVertex2 *vP,
			  int nV);
static WlzErrorNum 	WlzWriteVertex3I(
			  FILE *fP,
			  WlzIVertex3 *vP,
			  int nV);
static WlzErrorNum 	WlzWriteVertex3D(
			  FILE *fP,
			  WlzDVertex3 *vP,
			  int nV);
#ifdef WLZ_UNUSED_FUNCTIONS
static WlzErrorNum 	WlzWriteBox2I(
			  FILE *fP,
			  WlzIBox2 *bP,
			  int nB);
static WlzErrorNum 	WlzWriteBox2D(
			  FILE *fP,
			  WlzDBox2 *bP,
			  int nB);
static WlzErrorNum 	WlzWriteBox3I(
			  FILE *fP,
			  WlzIBox3 *bP,
			  int nB);
static WlzErrorNum 	WlzWriteBox3D(
			  FILE *fP,
			  WlzDBox3 *bP,
			  int nB);
#endif /* WLZ_UNUSED_FUNCTIONS */
static WlzErrorNum	WlzWriteStr(
			  FILE *fP,
			  char *str);
static WlzErrorNum	WlzWritePixelV(
			  FILE *fP,
			  WlzPixelV *pV,
			  int nPV);
static WlzErrorNum	WlzWriteGreyV(
			  FILE *fP,
			  WlzGreyType type,
			  WlzGreyV *gV,
			  int nGV);
static WlzErrorNum	WlzWriteMeshTransform2D(
			  FILE *fp,
			  WlzMeshTransform *mTrans);

#ifdef _OPENMP
#define putc(C,S) putc_unlocked(C,S)
#endif


/*!
* \return	Number of bytes written.
* \ingroup 	WlzIO
* \brief	Writes an integer reordered to DEC VAX(!) format.
* \param	i			Value written.
* \param	fp			Given file.
*/
static int putword(int i, FILE *fp)
{
  unsigned char *cin, cout[4];

  cin = (unsigned char *) &i;
#if defined (__sparc) || defined (__mips) || defined (__ppc)
  cout[0] = *(cin+3);
  cout[1] = *(cin+2);
  cout[2] = *(cin+1);
  cout[3] = *(cin+0);
#endif /* __sparc || __mips */
#if defined (__x86) || defined (__alpha)
  cout[0] = *(cin+0);
  cout[1] = *(cin+1);
  cout[2] = *(cin+2);
  cout[3] = *(cin+3);
#endif /* __x86 || __alpha */
  return( (int) fwrite(&cout[0], sizeof(char), 4, fp) );
}

/*!
* \return	Number of bytes written.
* \ingroup	WlzIO
* \brief	Writes a short reordered to DEC VAX(!) format.
* \param	i			Value written.
* \param	fp			Given file.
*/
static int putshort(short i, FILE *fp)
{
  unsigned char *cin, cout[2];

  cin = (unsigned char *) &i;
#if defined (__sparc) || defined (__mips) || defined (__ppc)
  cout[0] = *(cin+1);
  cout[1] = *(cin+0);
#endif /* __sparc || __mips */
#if defined (__x86) || defined (__alpha)
  cout[0] = *(cin+0);
  cout[1] = *(cin+1);
#endif /* __x86 || __alpha */
  return( (int) fwrite(&cout[0], sizeof(char), 2, fp) );
}

/*!
* \return	Number of bytes written.
* \ingroup	WlzIO
* \brief	Writes a float reordered to DEC VAX(!) format.
* \param	f			Value written.
* \param	fp			Given file.
*/
static int putfloat(float f, FILE *fp)
{
  float ff = f;
  unsigned char *cin, cout[4];

  cin = (unsigned char *) &ff;
#if defined (__sparc) || defined (__mips) || defined (__ppc)
  cout[0] = *(cin+1);
  cout[1] = *cin + 1;
  cout[2] = *(cin+3);
  cout[3] = *(cin+2);
#endif /* __sparc || __mips */
#if defined (__x86) || defined (__alpha)
  cout[3] = *(cin+1);
  cout[2] = *(cin+0);
  cout[1] = *(cin+3) + 1;
  cout[0] = *(cin+2);
#endif /* __x86 || __alpha */
  return( (int) fwrite(&cout[0], sizeof(char), 4, fp) );
}

/*!
* \return	Number of bytes written.
* \ingroup	WlzIO
* \brief	Writes a double reordered to DEC VAX(!) format.
* \param	f			Value written.
* \param	fp			Given file.
*/
static int putdouble(double d, FILE *fp)
{
  double dd = d;
  unsigned char *cin, cout[8];

  cin = (unsigned char *) &dd;
#if defined (__sparc) || defined (__mips) || defined (__ppc)
  cout[0] = *(cin+7);
  cout[1] = *(cin+6);
  cout[2] = *(cin+5);
  cout[3] = *(cin+4);
  cout[4] = *(cin+3);
  cout[5] = *(cin+2);
  cout[6] = *(cin+1);
  cout[7] = *cin;
#endif /* __sparc || __mips */
#if defined (__x86) || defined (__alpha)
  cout[7] = *(cin+7);
  cout[6] = *(cin+6);
  cout[5] = *(cin+5);
  cout[4] = *(cin+4);
  cout[3] = *(cin+3);
  cout[2] = *(cin+2);
  cout[1] = *(cin+1);
  cout[0] = *cin;
#endif /* __x86 || __alpha */
  return( (int) fwrite(&cout[0], sizeof(char), 8, fp) );
}

/*!
* \return       Woolz error number code.
* \ingroup      WlzIO
* \brief        Top-level procedure for writing an object to a file stream.
*		For historical reasons most data are written using DEC VAX
*		byte ordering.
*
* \param    	fp			File pointer for output.
* \param    	obj			Ptr to top-level object to be written.
*/
WlzErrorNum	WlzWriteObj(FILE *fp, WlzObject *obj)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

#ifdef _OPENMP
  #pragma omp critical
  {
#endif
  if(fp == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
#ifdef _WIN32
  else if(_setmode(_fileno(fp), 0x8000) == -1)
  {
    errNum = WLZ_ERR_READ_EOF;
  }
#endif
  else if(obj == NULL)
  {
    if(putc((unsigned int )obj, fp) == EOF)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
  }
  else if(putc((unsigned int )obj->type, fp) == EOF)
  {
    errNum = WLZ_ERR_WRITE_EOF;
  }
  else if(errNum == WLZ_ERR_NONE)
  {
    switch( obj->type )
    {
      case WLZ_EMPTY_OBJ:
	/* Nothing except the object type needs to be written */
	break;
      case WLZ_2D_DOMAINOBJ:
	if(((errNum = WlzWriteIntervalDomain(fp,
				obj->domain.i)) == WLZ_ERR_NONE) &&
	   ((errNum = WlzWriteValueTable(fp, obj)) == WLZ_ERR_NONE))
	{
	  errNum = WlzWritePropertyList(fp, obj->plist);
	}
	break;
      case WLZ_3D_DOMAINOBJ:
	if(((errNum = WlzWritePlaneDomain(fp,
				obj->domain.p)) == WLZ_ERR_NONE) &&
	   ((errNum = WlzWriteVoxelValueTable(fp, obj)) == WLZ_ERR_NONE))
	{
	  errNum = WlzWritePropertyList(fp, obj->plist);
	}
	break;
      case WLZ_TRANS_OBJ:
	if(((errNum = WlzWriteAffineTransform(fp,
				obj->domain.t)) == WLZ_ERR_NONE) &&
	   ((errNum = WlzWriteObj(fp, obj->values.obj)) == WLZ_ERR_NONE))
	{
	  errNum = WlzWritePropertyList(fp, obj->plist);
	}
	break;
      case WLZ_3D_WARP_TRANS:
	if((errNum = WlzWrite3DWarpTrans(fp,
				(Wlz3DWarpTrans *)obj)) == WLZ_ERR_NONE)
	{
	  errNum = WlzWritePropertyList(fp, ((Wlz3DWarpTrans *)obj)->plist);
	}
	break;
      case WLZ_2D_POLYGON:
	errNum = WlzWritePolygon(fp, obj->domain.poly);
        break;
      case WLZ_BOUNDLIST:
	errNum = WlzWriteBoundList(fp, obj->domain.b);
        break;
      case WLZ_HISTOGRAM:
	errNum = WlzWriteHistogramDomain(fp, obj->domain.hist);
	break;
      case WLZ_CONTOUR:
        errNum = WlzWriteContour(fp, obj->domain.ctr);
	break;
      case WLZ_RECTANGLE:
	errNum = WlzWriteRect(fp, obj->domain.r);
	break;
      case WLZ_AFFINE_TRANS:
	errNum = WlzWriteAffineTransform(fp, obj->domain.t);
	break;
      case WLZ_WARP_TRANS:
	errNum = WlzWriteWarpTrans(fp, (WlzWarpTrans *)obj);
	break;
      case WLZ_FMATCHOBJ:
	errNum = WlzWriteFMatchObj(fp, (WlzFMatchObj *)obj);
	break;
      case WLZ_COMPOUND_ARR_1: /* FALLTHROUGH */
      case WLZ_COMPOUND_ARR_2:
	errNum = WlzWriteCompoundA(fp, (WlzCompoundArray *)obj);
	break;
      case WLZ_PROPERTY_OBJ:
	errNum = WlzWritePropertyList(fp, obj->plist);
	break;
    case WLZ_CONV_HULL:
      if((errNum = WlzWritePolygon(fp, obj->domain.poly))
	 == WLZ_ERR_NONE){
	errNum = WlzWriteConvexHullValues(fp, obj->values.c);
      }
      break;
    case WLZ_3D_POLYGON:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    case WLZ_MESH_TRANS:
        errNum = WlzWriteMeshTransform2D(fp, obj->domain.mt);
      break;
	/* Orphans and not yet implemented object types for I/O */
      case WLZ_CONVOLVE_INT:    /* FALLTHROUGH */
      case WLZ_CONVOLVE_FLOAT:  /* FALLTHROUGH */
      case WLZ_TEXT:            /* FALLTHROUGH */
      case WLZ_COMPOUND_LIST_1: /* FALLTHROUGH */
      case WLZ_COMPOUND_LIST_2: /* FALLTHROUGH */
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
#ifdef _OPENMP
  }
#endif
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Write's the given native int values to the given file
*               as a 4 byte integer.
* \param	fP			Given file.
* \param	iP			Ptr to native ints.
* \param	nI			Number of ints.
*/
static WlzErrorNum WlzWriteInt(FILE *fP, int *iP, int nI)
{
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  while((errNum == WLZ_ERR_NONE) && (nI-- > 0))
  {
    if(!putword(*iP, fP))
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
    ++iP;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzIO
* \brief	Write's the given 2D integer verticies to the given file.
* \param	fP			Given file.
* \param	vP			Ptr to 2D integer verticies.
* \param	nV			Number of verticies.
*/
static WlzErrorNum WlzWriteVertex2I(FILE *fP, WlzIVertex2 *vP, int nV)
{
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  while((errNum == WLZ_ERR_NONE) && (nV-- > 0))
  {
    if(!putword(vP->vtY, fP) || !putword(vP->vtX, fP))
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
    ++vP;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzIO
* \brief	Write's the given 2D double verticies to the given file.
* \param	fP			Given file.
* \param	vP			Ptr to 2D double verticies.
* \param	nV			Number of verticies.
*/
static WlzErrorNum WlzWriteVertex2D(FILE *fP, WlzDVertex2 *vP, int nV)
{
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  while((errNum == WLZ_ERR_NONE) && (nV-- > 0))
  {
    if(!putdouble(vP->vtY, fP) || !putdouble(vP->vtX, fP))
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
    ++vP;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzIO
* \brief	Write's the given 3D integer verticies to the given file.
* \param	fP			Given file.
* \param	vP			Ptr to 3D integer verticies.
* \param	nV			Number of verticies.
*/
static WlzErrorNum WlzWriteVertex3I(FILE *fP, WlzIVertex3 *vP, int nV)
{
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  while((errNum == WLZ_ERR_NONE) && (nV-- > 0))
  {
    if(!putword(vP->vtX, fP) || !putword(vP->vtY, fP) ||
       !putword(vP->vtZ, fP))
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
    ++vP;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzIO
* \brief	Write's the given 3D double verticies to the given file.
* \param	fP			Given file stream.
* \param	vP			Ptr to 3D double verticies.
* \param	nV			Number of verticies.
*/
static WlzErrorNum WlzWriteVertex3D(FILE *fP, WlzDVertex3 *vP, int nV)
{
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  while((errNum == WLZ_ERR_NONE) && (nV-- > 0))
  {
    if(!putdouble(vP->vtX, fP) || !putdouble(vP->vtY, fP) ||
       !putdouble(vP->vtZ, fP))
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
    ++vP;
  }
  return(errNum);
}

#ifdef WLZ_UNUSED_FUNCTIONS
/*!
* \return	Woolz error code.
* \ingroup      WlzIO
* \brief	Write's the given 2D integer box to the given file.
* \param	fP			Given file.
* \param	bP			Ptr to 2D integer box.
* \param	nB			Number of bounding boxes.
*/
static WlzErrorNum WlzWriteBox2I(FILE *fP, WlzIBox2 *bP, int nB)
{
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  while((errNum == WLZ_ERR_NONE) && (nB-- > 0))
  {
    if(!putword(bP->xMin, fP) || !putword(bP->yMin, fP) ||
       !putword(bP->xMax, fP) || !putword(bP->yMax, fP))
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
    ++bP;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzIO
* \brief	Write's the given 2D double box to the given file.
* \param	fP			Given file.
* \param	bP			Ptr to 2D double box.
* \param	nB			Number of bounding boxes.
*/
static WlzErrorNum WlzWriteBox2D(FILE *fP, WlzDBox2 *bP, int nB)
{
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  while((errNum == WLZ_ERR_NONE) && (nB-- > 0))
  {
    if(!putdouble(bP->xMin, fP) || !putdouble(bP->yMin, fP) ||
       !putdouble(bP->xMax, fP) || !putdouble(bP->yMax, fP))
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
    ++bP;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Write's the given 3D integer box to the given file.
* \param	fP			Given file.
* \param	bP			Ptr to 3D integer box.
* \param	nB			Number of bounding boxes.
*/
static WlzErrorNum WlzWriteBox3I(FILE *fP, WlzIBox3 *bP, int nB)
{
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  while((errNum == WLZ_ERR_NONE) && (nB-- > 0))
  {
    if(!putword(bP->xMin, fP) || !putword(bP->yMin, fP) ||
       !putword(bP->zMin, fP) ||
       !putword(bP->xMax, fP) || !putword(bP->yMax, fP) ||
       !putword(bP->zMax, fP))
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
    ++bP;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzIO
* \brief	Write's the given 3D double box to the given file.
* \param	fP			Given file.
* \param	bP			Ptr to 3D double box.
* \param	nB			Number of bounding boxes.
*/
static WlzErrorNum WlzWriteBox3D(FILE *fP, WlzDBox3 *bP, int nB)
{
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  while((errNum == WLZ_ERR_NONE) && (nB-- > 0))
  {
    if(!putdouble(bP->xMin, fP) || !putdouble(bP->yMin, fP) ||
       !putdouble(bP->zMin, fP) ||
       !putdouble(bP->xMax, fP) || !putdouble(bP->yMax, fP) ||
       !putdouble(bP->zMax, fP))
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
    ++bP;
  }
  return(errNum);
}
#endif /* WLZ_UNUSED_FUNCTIONS */

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Writes an ASCII string to the given file.
*		\verbatim
		<number of characters> int, 4 bytes
		<characters>
		\endverbatim
* \param	fP			Given file.
* \param	str			The string.
*/
static WlzErrorNum	WlzWriteStr(FILE *fP, char *str)
{
  int		len;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  len = strlen(str);
  if((putword(len, fP) !=  4) || (fwrite(str, 1, len, fP) != len))
  {
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Writes the given grey values to the given file. Each
*		of the grey values may be of the same type.
* \param	fP			Given file.
* \param	type			Grey value type.
* \param	gV			The grey values.
* \param	nGV			Number of grey values to write.
*/
static WlzErrorNum	WlzWriteGreyV(FILE *fP, WlzGreyType type,
				      WlzGreyV *gV, int nGV)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(type)
  {
    case WLZ_GREY_INT:
      while((errNum == WLZ_ERR_NONE) && (nGV-- > 0))
      {
	if(putword(gV->inv, fP) != 4)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
	++gV;
      }
      break;
    case WLZ_GREY_SHORT:
      while((errNum == WLZ_ERR_NONE) && (nGV-- > 0))
      {
	if(putshort(gV->shv, fP) != 2)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
	++gV;
      }
      break;
    case WLZ_GREY_UBYTE:
      while((errNum == WLZ_ERR_NONE) && (nGV-- > 0))
      {
	if(putc(((unsigned int )(gV->ubv)), fP) == EOF)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
	++gV;
      }
      break;
    case WLZ_GREY_FLOAT:
      while((errNum == WLZ_ERR_NONE) && (nGV-- > 0))
      {
        if(putfloat(gV->flv, fP) != 4)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
	++gV;
      }
      break;
    case WLZ_GREY_DOUBLE:
      while((errNum == WLZ_ERR_NONE) && (nGV-- > 0))
      {
        if(putdouble(gV->dbv, fP) != 8)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
	++gV;
      }
      break;
    case WLZ_GREY_RGBA:
      while((errNum == WLZ_ERR_NONE) && (nGV-- > 0))
      {
	if(putword((int )(gV->rgbv), fP) != 4)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
	++gV;
      }
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Writes the given pixel values to the given file. Each
*		of the pixel values may be of a different type.
* \param	fP			Given file.
* \param	pV			The pixel values.
* \param	nGV			Number of grey values to write.
*/
static WlzErrorNum	WlzWritePixelV(FILE *fP, WlzPixelV *pV, int nPV)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  while((errNum == WLZ_ERR_NONE) && (nPV-- > 0))
  {
    if(putc((unsigned int)pV->type, fP) == EOF)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
    else
    {
      errNum = WlzWriteGreyV(fP, pV->type, &(pV->v), 1);
    }
    ++pV;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzIO
* \brief	Write's the given Woolz interval domain to the given file.
* \param	fp			Given file.
* \param	itvl			Interval domain.
*/
static WlzErrorNum WlzWriteIntervalDomain(FILE *fp, WlzIntervalDomain *itvl)
{
  int 			i,
  			j,
			nlines;
  WlzIntervalLine	*ivln;
  WlzErrorNum		errNum = WLZ_ERR_NONE;

  if(itvl == NULL)
  {
    if(putc(0,fp) == EOF)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
  }
  else
  {
    /* standardise it so no non-standard domains on disk
       - can't do this it conflicts with read-only access
       to object store */
    /*WlzStandardIntervalDomain(itvl);*/

    /* write the type and bounding box */
    if((putc((unsigned int) itvl->type, fp) == EOF) ||
       !putword(itvl->line1, fp) ||
       !putword(itvl->lastln, fp) ||
       !putword(itvl->kol1, fp) ||
       !putword(itvl->lastkl, fp) )
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
    else
    {
      switch(itvl->type)
      {
	case WLZ_INTERVALDOMAIN_INTVL:
	  nlines = itvl->lastln - itvl->line1;
	  for(i = 0; (i <= nlines) && (errNum == WLZ_ERR_NONE); i++)
	  {
	    if(!putword(itvl->intvlines[i].nintvs, fp))
	    {
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    ivln = itvl->intvlines;
	    for(i = 0; (i <= nlines) && (errNum == WLZ_ERR_NONE); i++)
	    {
	      for(j = 0; (j < ivln->nintvs) && (errNum == WLZ_ERR_NONE); j++)
	      {
		if(!putword(ivln->intvs[j].ileft, fp) ||
		   !putword(ivln->intvs[j].iright, fp))
		{
		  errNum = WLZ_ERR_WRITE_INCOMPLETE;
		}
	      }
	      ivln++;
	    }
	  }
	  break;
	case WLZ_INTERVALDOMAIN_RECT:
	  break;
	default:
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	  break;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Writes a plane domain to the given file.
* \param	fp			Given file.
* \param	planedm			Palne domain.
*/
static WlzErrorNum WlzWritePlaneDomain(FILE *fp, WlzPlaneDomain *planedm)
{
  int		i,
  		nplanes;
  float		dummy_float = 0.0;
  WlzDomain	*domains;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(planedm == NULL)
  {
    if(putc(0,fp) == EOF)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
  }
  else
  {
    if((putc((unsigned int )(planedm->type), fp) == EOF) ||
       !putword(planedm->plane1, fp) ||
       !putword(planedm->lastpl, fp) ||
       !putword(planedm->line1, fp) ||
       !putword(planedm->lastln, fp) ||
       !putword(planedm->kol1, fp) ||
       !putword(planedm->lastkl, fp) ||
       !putfloat((planedm->voxel_size)[0], fp) ||
       !putfloat((planedm->voxel_size)[1], fp) ||
       !putfloat((planedm->voxel_size)[2], fp))
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
    else
    {
      /* write dummy values of plane positions for backward
	 compatibility - should go on file-format revision */
      nplanes = planedm->lastpl - planedm->plane1 + 1;
      for(i = 0; (i < nplanes) && (errNum == WLZ_ERR_NONE); i++)
      {
	if(!putfloat(dummy_float, fp))
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      domains = planedm->domains;
      switch(planedm->type)
      {
	case WLZ_PLANEDOMAIN_DOMAIN:
	  for(i = 0; (i < nplanes) && (errNum == WLZ_ERR_NONE); i++, domains++)
	  {
	    errNum = WlzWriteIntervalDomain(fp, (*domains).i);
	  }
	  break;
	case WLZ_PLANEDOMAIN_POLYGON:
	  for(i = 0; (i < nplanes) && (errNum == WLZ_ERR_NONE); i++, domains++)
	  {
	    errNum = WlzWritePolygon(fp, (*domains).poly);
	  }
	  break;
	case WLZ_PLANEDOMAIN_BOUNDLIST:
	  for(i = 0; (i < nplanes) && (errNum == WLZ_ERR_NONE); i++, domains++)
	  {
	    errNum = WlzWriteBoundList(fp, (*domains).b);
	  }
	  break;
	case WLZ_PLANEDOMAIN_HISTOGRAM:
	  for(i = 0; (i < nplanes) && (errNum == WLZ_ERR_NONE); i++, domains++)
	  {
	    errNum = WlzWriteHistogramDomain(fp, (*domains).hist);
	  }
	  break;
	case WLZ_PLANEDOMAIN_AFFINE:
	  for(i = 0; (i < nplanes) && (errNum == WLZ_ERR_NONE); i++, domains++)
	  {
	    errNum = WlzWriteAffineTransform(fp, (*domains).t);
	  }
	  break;
	case WLZ_PLANEDOMAIN_WARP:
	  for(i = 0; (i < nplanes) && (errNum == WLZ_ERR_NONE); i++, domains++)
	  {
	    errNum = WlzWriteWarpTrans(fp, (*domains).wt);
	  }
	  break;
	default:
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	  break;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Writes a single property to the given file.
* \param	fp			Given file.
* \param	property		Property to be written.
*/
static WlzErrorNum WlzWriteProperty(FILE *fp, WlzProperty property)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(property.core == NULL)
  {
    if(putc(0,fp) == EOF)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
  }
  else
  {
    switch( property.core->type ){
    case WLZ_PROPERTY_SIMPLE:
      /* Write the size */
      if((putc((unsigned int) WLZ_PROPERTY_SIMPLE, fp) == EOF) ||
	 !putword(property.simple->size, fp))
      {
	errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
      else if(!fwrite(property.simple->prop, property.simple->size, 1, fp))
      {
	errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
      break;

    case  WLZ_PROPERTY_EMAP:
      if((putc((unsigned int) WLZ_PROPERTY_EMAP, fp) == EOF) ||
	 (putc((unsigned int) property.emap->emapType, fp) == EOF) ||
	 !fwrite(property.emap->modelUID,
		 EMAP_PROPERTY_UID_LENGTH, 1, fp) ||
	 !fwrite(property.emap->anatomyUID,
		 EMAP_PROPERTY_UID_LENGTH, 1, fp) ||
	 !fwrite(property.emap->targetUID,
		 EMAP_PROPERTY_UID_LENGTH, 1, fp) ||
	 !fwrite(property.emap->targetVersion,
		 EMAP_PROPERTY_VERSION_LENGTH, 1, fp) ||
	 !fwrite(property.emap->stage,
		 EMAP_PROPERTY_STAGE_LENGTH, 1, fp) ||
	 !fwrite(property.emap->subStage,
		 EMAP_PROPERTY_STAGE_LENGTH, 1, fp) ||
	 !fwrite(property.emap->modelName,
		 EMAP_PROPERTY_MODELNAME_LENGTH, 1, fp) ||
	 !fwrite(property.emap->version,
		 EMAP_PROPERTY_VERSION_LENGTH, 1, fp) ||
	 !putword(property.emap->creationTime, fp) ||
	 !fwrite(property.emap->creationAuthor,
		 EMAP_PROPERTY_AUTHORNAME_LENGTH, 1, fp) ||
	 !fwrite(property.emap->creationMachineName,
		 EMAP_PROPERTY_MACHINENAME_LENGTH, 1, fp) ||
	 !putword(property.emap->modificationTime, fp) ||
	 !fwrite(property.emap->modificationAuthor,
		 EMAP_PROPERTY_AUTHORNAME_LENGTH, 1, fp)){
	errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
      else {
	if(property.emap->fileName &&
	   (strlen(property.emap->fileName) > 0)){
	  if(!putword(strlen(property.emap->fileName), fp) ||
	     !fwrite(property.emap->fileName, strlen(property.emap->fileName),
		     1, fp)){
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	}
	else {
	  if(!putword(0, fp)){
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	}
	if((errNum == WLZ_ERR_NONE) && property.emap->comment &&
	   (strlen(property.emap->comment))){
	  if(!putword(strlen(property.emap->comment), fp) ||
	     !fwrite(property.emap->comment, strlen(property.emap->comment),
		     1, fp)){
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	}
	else {
	  if(!putword(0, fp)){
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	}
      }
      break;
    case WLZ_PROPERTY_NAME:
      if(putc((unsigned int) WLZ_PROPERTY_NAME, fp) == EOF)
      {
        errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
      else
      {
        errNum = WlzWriteStr(fp, property.name->name);
      }
      break;
    case WLZ_PROPERTY_GREY:
      if(putc((unsigned int) WLZ_PROPERTY_GREY, fp) == EOF)
      {
        errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
      else if((errNum = WlzWriteStr(fp, property.greyV->name)) == WLZ_ERR_NONE)
      {
        errNum = WlzWritePixelV(fp, &(property.greyV->value), 1);
      }
      break;
    default:
      errNum = WLZ_ERR_PROPERTY_TYPE;
      break;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Writes a property list to the given file.
* \param	fp			Given file.
* \param	pList			Property list.
*/
static WlzErrorNum WlzWritePropertyList(FILE *fp, WlzPropertyList *pList)
{
  AlcDLPItem	*item;
  WlzProperty	property;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((pList == NULL) || (pList->list == NULL) ||
     (AlcDLPListCount(pList->list, NULL) < 1))
   {
    if(putc(0,fp) == EOF)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
  }
  else
  {
    if((putc(2,fp) == EOF) || !putword(AlcDLPListCount(pList->list, NULL), fp))
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
    else
    {
      item = pList->list->head;
      do
      {
	if(item)
	{
	  if(item->entry != NULL)
	  {
	    property.core = (WlzCoreProperty *)(item->entry);
	    errNum = WlzWriteProperty(fp, property);
	  }
	  else
	  {
	    putc(0,fp);
	  }
	  item = item->next;
	}
      } while((errNum == WLZ_ERR_NONE) && (item != pList->list->head));
    }
  }
  return errNum;
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Writes the 2D values of a Woolz 2D domain object to the
*		given file.
* \param	fp			Given file.
* \param	obj			Object containing values that
*					are to be written to file.
*/
static WlzErrorNum WlzWriteValueTable(FILE *fp, WlzObject *obj)
{
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzGreyType		gType;
  WlzGreyP		g;
  WlzPixelV		background,
  			min,
			max;
  int 			i;
  WlzGreyType		packing;
  WlzErrorNum		errNum = WLZ_ERR_NONE;

  /* obj == NULL has been checked by WlzWriteObj() */
  if(obj->values.core == NULL)
  {
    if(putc(0,fp) == EOF)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
  }
  else
  {
    if(errNum == WLZ_ERR_NONE)
    {
      gType = WlzGreyTableTypeToGreyType(obj->values.core->type, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      /* The "type" written to disc only codes the pixel type.
         The shape type on subsequent reading is entirely
         determined by the object domain.
         For the moment, choice is between standard ragged-rectangle
         or true rectangle.  */
      if(putc((unsigned int )gType, fp) == EOF)
      {
	errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
    }
    background = WlzGetBackground(obj, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      switch(gType)
      {
	case WLZ_GREY_INT:
	  /* Calculate packing to minimise disc space */
	  if((errNum = WlzGreyRange(obj, &min, &max)) == WLZ_ERR_NONE)
	  {
	    if((min.v.inv >= 0) && (max.v.inv <= 255))
	    {
	      packing = WLZ_GREY_UBYTE;
	    }
	    else if((min.v.inv >= SHRT_MIN) && (max.v.inv <= SHRT_MAX))
	    {
	      packing = WLZ_GREY_SHORT;
	    }
	    else
	    {
	      packing = WLZ_GREY_INT;
	    }
	    if((putc((unsigned int )packing, fp) == EOF) ||
	       !putword(background.v.inv, fp))
            {
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    }
	    if((errNum = WlzInitGreyScan(obj, &iwsp, &gwsp)) == WLZ_ERR_NONE)
	    {
	      while((errNum == WLZ_ERR_NONE) &&
	            ((errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE))
	      {
		g = gwsp.u_grintptr;
		switch(packing)
		{
		  case WLZ_GREY_INT:
		    for(i = 0; (i < iwsp.colrmn) && (errNum == WLZ_ERR_NONE);
			i++, g.inp++)
		    {
		      if(!putword(*g.inp, fp))
		      {
			errNum = WLZ_ERR_WRITE_INCOMPLETE;
		      }
		    }
		    break;
		  case WLZ_GREY_SHORT:
		    for(i = 0; (i < iwsp.colrmn) && (errNum == WLZ_ERR_NONE);
			i++, g.inp++)
		    {
		      if(!putshort((short) (*g.inp), fp))
		      {
			errNum = WLZ_ERR_WRITE_INCOMPLETE;
		      }
		    }
		    break;
		  case WLZ_GREY_UBYTE:
		    for(i = 0; (i < iwsp.colrmn) && (errNum == WLZ_ERR_NONE);
		    	i++, g.inp++)
		    {
		      if(putc((unsigned int )*g.inp, fp) == EOF)
		      {
			errNum = WLZ_ERR_WRITE_INCOMPLETE;
		      }
		    }
		    break;
		  default:
		    errNum = WLZ_ERR_GREY_TYPE;
		    break;
		}
	      }
	      if(errNum == WLZ_ERR_EOO)
	      {
	        errNum = WLZ_ERR_NONE;
	      }
	    }
	  }
	  break;
	case WLZ_GREY_SHORT:
	  /* Calculate packing to minimise disc space */
	  if((errNum = WlzGreyRange(obj, &min, &max)) == WLZ_ERR_NONE)
	  {
	    if((min.v.shv >= 0) && (max.v.shv <= 255))
	    {
	      packing = WLZ_GREY_UBYTE;
	    }
	    else
	    {
	      packing = WLZ_GREY_SHORT;
	    }
	    if((putc((unsigned int )packing, fp) == EOF) ||
	       !putword(background.v.shv, fp))
	    {
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    }
	    if((errNum = WlzInitGreyScan(obj, &iwsp, &gwsp)) == WLZ_ERR_NONE)
	    {
	      while((errNum == WLZ_ERR_NONE) &&
	            ((errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE))
	      {
		g = gwsp.u_grintptr;
		switch(packing)
		{
		  case WLZ_GREY_SHORT:
		    for(i = 0; (i < iwsp.colrmn) && (errNum == WLZ_ERR_NONE);
		    	i++, g.shp++)
		    {
		      if(!putshort(*g.shp, fp))
		      {
			errNum = WLZ_ERR_WRITE_INCOMPLETE;
		      }
		    }
		    break;
		  case WLZ_GREY_UBYTE:
		    for(i = 0; (i < iwsp.colrmn) && (errNum == WLZ_ERR_NONE);
		    	i++, g.shp++)
		    {
		      if(putc((unsigned int) *g.shp, fp) == EOF)
		      {
			errNum = WLZ_ERR_WRITE_INCOMPLETE;
		      }
		    }
		    break;
		  default:
		    errNum = WLZ_ERR_GREY_TYPE;
		    break;
		}
	      }
	      if(errNum == WLZ_ERR_EOO)
	      {
	        errNum = WLZ_ERR_NONE;
	      }
	    }
	  }
	  break;
	case WLZ_GREY_UBYTE:
	  packing = WLZ_GREY_UBYTE;
	  if((putc((unsigned int )packing, fp) == EOF) ||
	     !putword(background.v.ubv, fp))
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	  else
	  {
	    errNum = WlzInitGreyScan(obj, &iwsp, &gwsp);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    while((errNum == WLZ_ERR_NONE) &&
	           ((errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE))
	    {
	      g = gwsp.u_grintptr;
	      if((int )fwrite(g.ubp, sizeof(WlzUByte),
	      		      iwsp.colrmn, fp) < iwsp.colrmn)
	      {
		errNum = WLZ_ERR_WRITE_INCOMPLETE;
	      }
	    }
	    if(errNum == WLZ_ERR_EOO)
	    {
	      errNum = WLZ_ERR_NONE;
	    }
	  }
	  break;
	case WLZ_GREY_FLOAT:
	  packing = WLZ_GREY_FLOAT;
	  if((putc((unsigned int )packing, fp) == EOF) ||
	     !putfloat(background.v.flv, fp))
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	  else
	  {
	    errNum = WlzInitGreyScan(obj, &iwsp, &gwsp);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    while((errNum == WLZ_ERR_NONE) &&
	           ((errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE))
	    {
	      g = gwsp.u_grintptr;
	      for(i = 0; (i < iwsp.colrmn) && (errNum == WLZ_ERR_NONE);
		  i++, g.flp++)
	      {
		if(!putfloat(*g.flp, fp))
		{
		  errNum = WLZ_ERR_WRITE_INCOMPLETE;
		}
	      }
	    }
	    if(errNum == WLZ_ERR_EOO)
	    {
	      errNum = WLZ_ERR_NONE;
	    }
	  }
	  break;
	case WLZ_GREY_DOUBLE:
	  packing = WLZ_GREY_DOUBLE;
	  if((putc((unsigned int )packing, fp) == EOF) ||
	     !putdouble(background.v.dbv, fp))
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	  else
	  {
	    errNum = WlzInitGreyScan(obj, &iwsp, &gwsp);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    while((errNum == WLZ_ERR_NONE) &&
	           ((errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE))
	    {
	      g = gwsp.u_grintptr;
	      for(i = 0; (i < iwsp.colrmn) && (errNum == WLZ_ERR_NONE);
	          i++, g.dbp++)
	      {
		if(!putdouble(*g.dbp, fp))
		{
		  errNum = WLZ_ERR_WRITE_INCOMPLETE;
		}
	      }
	    }
	    if(errNum == WLZ_ERR_EOO)
	    {
	      errNum = WLZ_ERR_NONE;
	    }
	  }
	  break;
	case WLZ_GREY_RGBA:
	  packing = WLZ_GREY_RGBA;
	  if((putc((unsigned int )packing, fp) == EOF) ||
	     !putword(background.v.rgbv, fp))
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	  else
	  {
	    errNum = WlzInitGreyScan(obj, &iwsp, &gwsp);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    while((errNum == WLZ_ERR_NONE) &&
	           ((errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE))
	    {
	      g = gwsp.u_grintptr;
	      for(i = 0; (i < iwsp.colrmn) && (errNum == WLZ_ERR_NONE);
	          i++, g.rgbp++)
	      {
		if(!putword(*g.rgbp, fp))
		{
		  errNum = WLZ_ERR_WRITE_INCOMPLETE;
		}
	      }
	    }
	    if(errNum == WLZ_ERR_EOO)
	    {
	      errNum = WLZ_ERR_NONE;
	    }
	  }
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Writes the voxel values of a Woolz object to the given file.
* \param	fp			Given file.
* \param	obj			Object with values.
*/
static WlzErrorNum WlzWriteVoxelValueTable(FILE *fp, WlzObject *obj)
{
  int			i, nplanes;
  WlzObject		tempobj;
  WlzDomain 		*domains;
  WlzValues		*values;
  WlzVoxelValues	*voxtab;
  WlzPlaneDomain	*planedm;
  WlzErrorNum		errNum = WLZ_ERR_NONE;

  /* check object */
  if(obj->values.core == NULL)
  {
    if(putc(0,fp) == EOF)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
  }
  else
  {
    voxtab = (WlzVoxelValues *) obj->values.vox;
    /* note here the background is written without a type and as an
       integer. On read the value is replaced by a value from one of
       the plane valuetables */
    if((putc((unsigned int) voxtab->type, fp) == EOF) ||
	!putword(voxtab->bckgrnd.v.inv, fp))
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
    else
    {
      planedm = obj->domain.p;
      switch(voxtab->type)
      {
	case WLZ_VOXELVALUETABLE_GREY:
	  nplanes = planedm->lastpl - planedm->plane1 + 1;
	  values = voxtab->values;
	  domains = planedm->domains;
	  tempobj.type = WLZ_2D_DOMAINOBJ;
	  tempobj.linkcount = 0;
	  tempobj.plist = NULL;
	  tempobj.assoc = NULL;
	  for(i=0; (i < nplanes) && (errNum == WLZ_ERR_NONE);
	      i++, domains++, values++)
	  {
	    tempobj.domain.i = (*domains).i;
	    tempobj.values.v = (*values).v;
	    errNum = WlzWriteValueTable(fp, &tempobj);
	  }
	  break;
	default:
	  errNum = WLZ_ERR_VALUES_TYPE;
	  break;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup 	WlzIO
* \brief	Writes a polygon domain to the given file.
* \param	fp			Given file.
* \param	poly			Polygon domain.
*/
static WlzErrorNum WlzWritePolygon(FILE *fp, WlzPolygonDomain *poly)
{
  int		nvertices, i;
  WlzIVertex2	*ivtx;
  WlzFVertex2	*fvtx;
  WlzDVertex2	*dvtx;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(poly == NULL)
  {
    if(putc(0,fp) == EOF)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
  }
  else
  {
    /* Write the number of vertices */
    nvertices = poly->nvertices;
    if((putc((unsigned int )poly->type, fp) == EOF) ||
       !putword(nvertices, fp))
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE && poly)
  {
    switch(poly->type)
    {
      case WLZ_POLYGON_INT:
	ivtx = poly->vtx;
	for(i = 0; (i < nvertices) && (errNum == WLZ_ERR_NONE); i++, ivtx++)
	{
	  if(!putword(ivtx->vtY, fp) || !putword(ivtx->vtX, fp))
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	}
	break;
      case WLZ_POLYGON_FLOAT:
	fvtx = (WlzFVertex2 *)poly->vtx;
	for(i = 0; (i < nvertices) && (errNum == WLZ_ERR_NONE); i++, fvtx++)
	{
	  if(!putfloat(fvtx->vtY, fp) || !putfloat(fvtx->vtX, fp))
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	}
	break;
      case WLZ_POLYGON_DOUBLE:
	dvtx = (WlzDVertex2 *)poly->vtx;
	for(i = 0; (i < nvertices) && (errNum == WLZ_ERR_NONE); i++, dvtx++)
	{
	  if(!putdouble(dvtx->vtY, fp) || !putdouble(dvtx->vtX, fp))
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	}
	break;
      default:
	errNum = WLZ_ERR_POLYGON_TYPE;
	break;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Writes a boundary list to the given file.
* \param	fp			Given file.
* \param	blist			Boundary list.
*/
static WlzErrorNum WlzWriteBoundList(FILE *fp, WlzBoundList *blist)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(blist == NULL)
  {
    if(putc(0,fp) == EOF)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
  }
  else
  {
    if((putc(1, fp) == EOF) ||
       (putc((unsigned int)blist->type, fp) == EOF))
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
    else if(((errNum = WlzWriteBoundList(fp, blist->next)) == WLZ_ERR_NONE) &&
	    ((errNum = WlzWriteBoundList(fp, blist->down)) == WLZ_ERR_NONE))
    {
      if(!putword(blist->wrap, fp))
      {
        errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
      else
      {
	errNum = WlzWritePolygon(fp, blist->poly);
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Writes a convex hull to the given file.
* \param	fp			Given file.
* \param	cnvhull			Convex hull.
*/
static WlzErrorNum      WlzWriteConvexHullValues(
  FILE *fp,
  WlzConvHullValues *cnvhull)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzChord	*ch;
  int		i;

  if(cnvhull == NULL)
  {
    if(putc(0,fp) == EOF)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
  }
  else
  {
    if((putc((unsigned int)cnvhull->type, fp) == EOF))
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
    else if(putword(cnvhull->nchords, fp) &&
	    putword(cnvhull->nsigchords, fp) &&
	    putword(cnvhull->mdlin, fp) &&
	    putword(cnvhull->mdkol, fp)){
      for(i=0, ch = cnvhull->ch; i < cnvhull->nchords; i++, ch++){
	if(!putword(ch->sig, fp) ||
	   !putword(ch->acon, fp) ||
	   !putword(ch->bcon, fp) ||
	   !putword(ch->ccon, fp) ||
	   !putdouble(ch->cl, fp) ||
	   !putword(ch->bl, fp) ||
	   !putword(ch->bk, fp) ||
	   !putword(ch->barea, fp) ||
	   !putword(ch->bd, fp)){
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  break;
	}
      }
    }
    else {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Writes an integer rectangle to the given file.
* \param	fp			Given file.
* \param	rdom			Integer rectangle.
*/
static WlzErrorNum WlzWriteRect(FILE *fp, WlzIRect *rdom)
{
  WlzFRect	*frdom;
  int		i;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(rdom == NULL)
  {
    if(putc(0,fp) == EOF)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
  }
  else
  {
    if(putc(rdom->type,fp) == EOF)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
    else
    {
      switch(rdom->type)
      {
	case WLZ_RECTANGLE_DOMAIN_INT:
	  for(i = 0; (i < 4) && (errNum == WLZ_ERR_NONE); ++i)
	  {
	    if(!putword(rdom->irk[i], fp))
	    {
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    }
	  }
	  for(i = 0; (i < 4) && (errNum == WLZ_ERR_NONE); ++i)
	  {
	    if(!putword(rdom->irl[i], fp))
	    {
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    }
	  }
	  if(!putfloat(rdom->rangle, fp))
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	  break;
	case WLZ_RECTANGLE_DOMAIN_FLOAT:
	  frdom = (WlzFRect *) rdom;
	  for(i = 0; (i < 4) && (errNum == WLZ_ERR_NONE); ++i)
	  {
	    if(!putfloat(frdom->frk[i], fp))
	    {
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    }
	  }
	  for(i = 0; (i < 4) && (errNum == WLZ_ERR_NONE); ++i)
	  {
	    if(!putfloat(frdom->frl[i], fp))
	    {
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    }
	  }
	  if(!putfloat(frdom->rangle, fp))
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	  break;
	default:
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	  break;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Writes a histogram domain to the given file.
* \param	fp			Given file.
* \param	hist			Histogram domain.
*/
static WlzErrorNum WlzWriteHistogramDomain(FILE *fp, WlzHistogramDomain *hist)
{
  int		tI0;
  int		*tIP0;
  double	*tDP0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(hist == NULL)
  {
    if(putc(0, fp) == EOF)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
  }
  else
  {
    if((putc(hist->type, fp) == EOF) ||
       (putword((hist->nBins > 0)? (hist->nBins): 0, fp) == 0) ||
       (putdouble(hist->origin, fp) == 0) ||
       (putdouble(hist->binSize, fp) == 0))
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
    else if((hist->type != WLZ_HISTOGRAMDOMAIN_INT) &&
            (hist->type != WLZ_HISTOGRAMDOMAIN_FLOAT))
    {
      errNum = WLZ_ERR_DOMAIN_TYPE;
    }
    else
    {
      if((tI0 = hist->nBins) > 0)
      {
	switch(hist->type)
	{
	  case WLZ_HISTOGRAMDOMAIN_INT:
	    tIP0 = hist->binValues.inp;
	    do
	    {
	      if(!putword(*tIP0++, fp))
	      {
	        errNum = WLZ_ERR_WRITE_INCOMPLETE;
	      }
	    } while((--tI0 > 0) && (errNum == WLZ_ERR_NONE));
	    break;
	  case WLZ_HISTOGRAMDOMAIN_FLOAT:
	    tDP0 = hist->binValues.dbp;
	    do
	    {
	      if(!putdouble(*tDP0++, fp))
	      {
	        errNum = WLZ_ERR_WRITE_INCOMPLETE;
	      }
	    } while((--tI0 > 0) && (errNum == WLZ_ERR_NONE));
	    break;
	  default:
	    break;
	}
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Writes a compound array object to the given file.
* \param	fp			Given file.
* \param	c			Compound array object.
*/
static WlzErrorNum WlzWriteCompoundA(FILE *fp, WlzCompoundArray *c)
{
  int 		i;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* No need to check for NULL because checked by WlzWriteObj() */
  if((putc((unsigned int )c->otype, fp) == EOF) || !putword(c->n, fp))
  {
    errNum = WLZ_ERR_WRITE_EOF;
  }
  else
  {
    /* Write the objects, note NULL is a legal value */
    for(i = 0; (i < c->n) && (errNum == WLZ_ERR_NONE); i++)
    {
      if(c->o[i] == NULL)
      {
	if(putc(0, fp) == EOF)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
      }
      else
      {
	errNum = WlzWriteObj(fp, c->o[i]);
      }
    }
  }
  if( errNum == WLZ_ERR_NONE )
  {
    errNum = WlzWritePropertyList(fp, c->plist);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Writes a Woolz affine transform to the given file.
* \param	fp			Given file.
* \param	trans			Affine transform.
*/
static WlzErrorNum WlzWriteAffineTransform(FILE *fp, WlzAffineTransform *trans)
{
  int		i,
  		j;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* check for NULL */
  if(trans == NULL)
  {
    if(putc((unsigned int )0, fp) == EOF)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
  }
  else
  {
    if(putc((unsigned int) trans->type, fp) == EOF)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      for(i = 0; (i < 4) && (errNum == WLZ_ERR_NONE); i++)
      {
	for(j = 0; (j < 4) && (errNum == WLZ_ERR_NONE); j++)
	{
	  if(!putdouble(trans->mat[i][j], fp))
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	}
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Writes a FE warp transform to the given file.
* \param	fp			Given file.
* \param	obj			Warp transform.
*/
static WlzErrorNum WlzWriteWarpTrans(FILE *fp, WlzWarpTrans *obj)
{
  int		i,
  		j;
  WlzDVertex2	*dptr;
  WlzTElement	*eptr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(!putc((unsigned int )obj->type, fp))
  {
    errNum = WLZ_ERR_WRITE_EOF;
  }
  else if(!putword(obj->nelts,fp) ||
	  !putword(obj->nodes,fp) ||
	  !putfloat(obj->imdisp,fp) ||
	  !putfloat(obj->iterdisp,fp))
  {
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }
  else
  {
    /* Write out nodal coords */
    dptr = obj->ncoords;
    for(i = 0; (i < obj->nodes) && (errNum == WLZ_ERR_NONE); i++, dptr++)
    {
      if(!putfloat((float) dptr->vtX, fp) ||
	 !putfloat((float) dptr->vtY, fp))
      {
	errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Write out nodal displacements */
    dptr = obj->displacements;
    for(i = 0; (i < obj->nodes) && (errNum == WLZ_ERR_NONE); i++, dptr++)
    {
      if(!putfloat((float) dptr->vtX, fp) ||
	 !putfloat((float) dptr->vtY, fp))
      {
	errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Write out elements */
    eptr = obj->eltlist;
    for(i = 0; (i < obj->nelts) && (errNum == WLZ_ERR_NONE); i++, eptr++)
    {
      if((putc((unsigned int) eptr->type, fp) == EOF) ||
         !putword(eptr->n, fp))
      {
        errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
      for(j = 0; (j < 3) && (errNum == WLZ_ERR_NONE); j++)
      {
	if (!putword(eptr->nodes[j], fp))
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
      }
      for(j = 0; (j < 3) && (errNum == WLZ_ERR_NONE); j++)
      {
	if (!putfloat(eptr->u[j], fp))
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
      }
      for(j = 0; (j < 3) && (errNum == WLZ_ERR_NONE); j++)
      {
	if (!putfloat(eptr->a[j], fp))
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Writes FE warp match features to the given file.
* \param	fp			Given file.
* \param	obj			Match features.
*/
static WlzErrorNum WlzWriteFMatchObj(FILE *fp, WlzFMatchObj *obj)
{
  int		i,
  		j;
  WlzFMatchPoint *mptr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(!putword(obj->nopts,fp))
  {
    errNum = WLZ_ERR_WRITE_EOF;
  }
  else
  {
    mptr = obj->matchpts;
    for(i = 0; (i < obj->nopts) && (errNum == WLZ_ERR_NONE); i++, mptr++)
    {
      if(!putword(mptr->type,fp) ||
	 !putword(mptr->node,fp) ||
	 !putfloat(mptr->ptcoords.vtX,fp) ||
	 !putfloat(mptr->ptcoords.vtY,fp))
      {
	errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
      for(j = 0; (j < WLZ_MAX_NODAL_DEGREE) && (errNum == WLZ_ERR_NONE); j++)
      {
	if (!putword(mptr->elements[j],fp))
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Writes a 3D FE warp transform to the given file.
* \param	fp			Given file.
* \param	obj			Warp transform.
*/
static WlzErrorNum WlzWrite3DWarpTrans(FILE *fp, Wlz3DWarpTrans *obj)
{
  int 		i;
  WlzFMatchObj	**intdoms;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(!putword(obj->iteration,fp))
  {
    errNum = WLZ_ERR_WRITE_EOF;
  }
  else if(!putword(obj->currentplane,fp) ||
           !putfloat(obj->maxdisp,fp))
  {
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }
  else if(!WlzWritePlaneDomain(fp, obj->pdom))
  {
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }
  else
  {
    intdoms = obj->intptdoms;
    for(i = obj->pdom->plane1;
        (i <= obj->pdom->lastpl) && (errNum == WLZ_ERR_NONE); i++, intdoms++)
    {
      if(!WlzWriteFMatchObj(fp, *intdoms))
      {
	errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Writes either a 2D or 3D contour to the given file.
* \param	fP			Given file.
* \param	ctr			Contour.
*/
static WlzErrorNum WlzWriteContour(FILE *fP, WlzContour *ctr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(ctr == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(ctr->type != WLZ_CONTOUR)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    if(putc((unsigned int )(ctr->type), fP) == EOF)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
    else
    {
      errNum = WlzWriteGMModel(fP, ctr->model);
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Write's a geometric model data structure to the given file.
*		The format is:
*		\verbatim
 		  model->type (byte)
 		  {
 		    case WLZ_GMMOD_2I
 		    case WLZ_GMMOD_2D
 		    case WLZ_GMMOD_2N
 		      nEdge (int)
 		    case WLZ_GMMOD_3I
 		    case WLZ_GMMOD_3D
 		    case WLZ_GMMOD_3N
 		      nLoop (int)
 		  }
 		  {
 		    case model->type == WLZ_GMMOD_2I
 		      vertexGU.vg2I->vtx (WlzIVertex2, int * 2)
 		    case model->type == WLZ_GMMOD_2D
 		      vertexGU.vg2D->vtx (WlzDVertex2, double * 2)
 		    case model->type == WLZ_GMMOD_2N
 		      vertexGU.vg2D->vtx (WlzDVertex2, double * 2)
 		      vertexGU.vg2D->nrm (WlzDVertex2, double * 2)
 		    case model->type == WLZ_GMMOD_3I
 		      vertexGU.vg3I->vtx (WlzIVertex2, int * 3)
 		    case model->type == WLZ_GMMOD_3D
 		      vertexGU.vg3D->vtx (WlzDVertex2, double * 3)
 		    case model->type == WLZ_GMMOD_3N
 		      vertexGU.vg3D->vtx (WlzDVertex2, double * 3)
 		      vertexGU.vg3D->nrm (WlzDVertex2, double * 3)
 		  } * nVertex
 		  {
 		    case WLZ_GMMOD_2I
 		    case WLZ_GMMOD_2D
 		    case WLZ_GMMOD_2N
 		    {
 		      2 vertex indicies
 		    } * nEdge
 		    case WLZ_GMMOD_3I
 		    case WLZ_GMMOD_3D
 		    case WLZ_GMMOD_3N
 		    {
 		      3 vertex indicies
 		    } * nLoop
 		  }
		\endverbatim
* \param	fP			Given file.
* \param	model			Geometric model.
*/
static WlzErrorNum WlzWriteGMModel(FILE *fP, WlzGMModel *model)
{

  int		idI,
  		iCnt,
		vCnt,
		encodeMtd = 0;
  int		bufI[3];
  AlcVector	*vec;
  WlzGMEdgeT	*tET;
  WlzGMElemP	eP;
  WlzGMResIdxTb	*resIdxTb = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    /* Check the model type. */
    switch(model->type)
    {
      case WLZ_GMMOD_2I: /* FALLTHROUGH */
      case WLZ_GMMOD_2D: /* FALLTHROUGH */
      case WLZ_GMMOD_2N: /* FALLTHROUGH */
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D: /* FALLTHROUGH */
      case WLZ_GMMOD_3N:
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Output model type and file encoding method, followed by the number of
     * verticies and the number of simplicies. Currently the file encoding
     * method is just output as '0', but in future new encoding methods
     * may be used, for example strips of simplicies to reduce the file
     * size. */
    switch(model->type)
    {
      case WLZ_GMMOD_2I: /* FALLTHROUGH */
      case WLZ_GMMOD_2D: /* FALLTHROUGH */
      case WLZ_GMMOD_2N:
	if((putc((unsigned int )(model->type), fP) == EOF) ||
	    (putc((unsigned int )encodeMtd, fP) == EOF) ||
	    !putword(model->res.vertex.numElm, fP) ||
	    !putword(model->res.edge.numElm, fP))
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
	break;
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D: /* FALLTHROUGH */
      case WLZ_GMMOD_3N:
	if((putc((unsigned int )(model->type), fP) == EOF) ||
	    (putc((unsigned int )encodeMtd, fP) == EOF) ||
	    !putword(model->res.vertex.numElm, fP) ||
	    !putword(model->res.face.numElm, fP))
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
	break;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (model->res.vertex.numElm > 0))
  {
      /* Index the verticies. */
      resIdxTb = WlzGMModelResIdx(model, WLZ_GMELMFLG_VERTEX, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      /* Output the vertex geometries. */
      idI = 0;
      vec = model->res.vertex.vec;
      iCnt = model->res.vertex.numIdx;
      vCnt = 0;
      while((errNum == WLZ_ERR_NONE) && (iCnt-- > 0))
      {
	eP.vertex = (WlzGMVertex *)AlcVectorItemGet(vec, idI++);
	if(eP.vertex->idx >= 0)
	{
	  ++vCnt;
	  switch(model->type)
	  {
	    case WLZ_GMMOD_2I:
	      errNum = WlzWriteVertex2I(fP, &(eP.vertex->geo.vg2I->vtx), 1);
	      break;
	    case WLZ_GMMOD_2D:
	      errNum = WlzWriteVertex2D(fP, &(eP.vertex->geo.vg2D->vtx), 1);
	      break;
	    case WLZ_GMMOD_2N:
	      errNum = WlzWriteVertex2D(fP, &(eP.vertex->geo.vg2N->vtx), 1);
	      if(errNum == WLZ_ERR_NONE)
	      {
	        errNum = WlzWriteVertex2D(fP, &(eP.vertex->geo.vg2N->nrm), 1);
	      }
	      break;
	    case WLZ_GMMOD_3I:
	      errNum = WlzWriteVertex3I(fP, &(eP.vertex->geo.vg3I->vtx), 1);
	      break;
	    case WLZ_GMMOD_3D:
	      errNum = WlzWriteVertex3D(fP, &(eP.vertex->geo.vg3D->vtx), 1);
	      break;
	    case WLZ_GMMOD_3N:
	      errNum = WlzWriteVertex3D(fP, &(eP.vertex->geo.vg3N->vtx), 1);
	      if(errNum == WLZ_ERR_NONE)
	      {
	        errNum = WlzWriteVertex3D(fP, &(eP.vertex->geo.vg3N->nrm), 1);
	      }
	      break;
	  }
	}
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      /* Output the vertex indicies of the simplicies. */
      idI = 0;
      vCnt = 0;
      switch(model->type)
      {
	case WLZ_GMMOD_2I:
	case WLZ_GMMOD_2D:
	case WLZ_GMMOD_2N:
	  vec = model->res.edge.vec;
	  iCnt = model->res.edge.numIdx;
	  while((errNum == WLZ_ERR_NONE) && (iCnt-- > 0))
	  {
	    ++vCnt;
	    eP.edge = (WlzGMEdge *)AlcVectorItemGet(vec, idI++);
	    if(eP.edge->idx >= 0)
	    {
	      tET = eP.edge->edgeT;
	      bufI[0] = *(resIdxTb->vertex.idxLut +
			  tET->vertexT->diskT->vertex->idx);
	      bufI[1] = *(resIdxTb->vertex.idxLut +
			  tET->opp->vertexT->diskT->vertex->idx);
	      errNum = WlzWriteInt(fP, bufI, 2);
	    }
	  }
	  break;
	case WLZ_GMMOD_3I:
	case WLZ_GMMOD_3D:
	case WLZ_GMMOD_3N:
	  vec = model->res.face.vec;
	  iCnt = model->res.face.numIdx;
	  while((errNum == WLZ_ERR_NONE) && (iCnt-- > 0))
	  {
	    eP.face = (WlzGMFace *)AlcVectorItemGet(vec, idI++);
	    if(eP.face->idx >= 0)
	    {
	      ++vCnt;
	      /* Loop IS a triangle, in 3D nothing else is allowed. */
	      tET = eP.face->loopT->edgeT;
	      bufI[0] = *(resIdxTb->vertex.idxLut +
			  tET->vertexT->diskT->vertex->idx);
	      bufI[1] = *(resIdxTb->vertex.idxLut +
			  tET->next->vertexT->diskT->vertex->idx);
	      bufI[2] = *(resIdxTb->vertex.idxLut +
			  tET->prev->vertexT->diskT->vertex->idx);
	      errNum = WlzWriteInt(fP, bufI, 3);
	    }
	  }
	  break;
      }
    }
  }
  if(resIdxTb)
  {
    WlzGMModelResIdxFree(resIdxTb);
  }
  return(errNum);
}

/* function:     WlzWriteMeshTransform3D    */
/*! 
* \ingroup      WlzIO
* \author	J. Rao.
* \brief        Write a 3D mesh transform to the given file-stream.
*
* \return       Error number.
* \param    fp	Output file-stream pointer
* \param    obj	Mesh transform to be written.
* \par      Source:
*                WlzWriteObj.c
*/
WlzErrorNum WlzWriteMeshTransform3D(
  FILE 			*fp,
  WlzMeshTransform3D 	*obj)
{
  int		i,
  		j;
  WlzMeshNode3D	*dptr;
  WlzMeshElem3D	*eptr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  /*
  if( !putc((unsigned int )obj->type, fp) )
  {
    errNum = WLZ_ERR_WRITE_EOF;
  }
  */
  if(     !putword(obj->nElem,  fp) ||
	  !putword(obj->nNodes ,fp) )
  {
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }
  else
  {
    /* Write out nodal position and displacement */
    dptr = obj->nodes;
    for(i = 0; (i < obj->nNodes) && (errNum == WLZ_ERR_NONE); i++, dptr++)
    {  /* position */
      if(!putfloat((float) dptr->position.vtX, fp) ||
         !putfloat((float) dptr->position.vtY, fp) ||
	 !putfloat((float) dptr->position.vtZ, fp))
      {
	errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
      /* displacement */
      if(!putfloat((float) dptr->displacement.vtX, fp) ||
         !putfloat((float) dptr->displacement.vtY, fp) ||
	 !putfloat((float) dptr->displacement.vtZ, fp))
      {
	errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }

    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Write out elements */
    eptr = obj->elements;
    for(i = 0; (i < obj->nElem) && (errNum == WLZ_ERR_NONE); i++, eptr++)
    {/*
      if(!putc((unsigned int )eptr->type, fp)  )
      {
        errNum = WLZ_ERR_WRITE_INCOMPLETE;
      } */
      /* output the index of this element */
      if(!putword(eptr->idx, fp))
      {
        errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
      /* output nodes indeces */
      for(j = 0; (j < 4) && (errNum == WLZ_ERR_NONE); j++)
      {
	if (!putword(eptr->nodes[j], fp))
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
      }
      /* output its neighbours */
      /*
      for(j = 0; (j < 4) && (errNum == WLZ_ERR_NONE); j++)
      {
	if (!putword(eptr->neighbours[j], fp))
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
      }
      */
    }
  }
    /* now Write out displacement */
    /*
  if(errNum == WLZ_ERR_NONE)
  {
    dptr = obj->nodes;
    for(i = 0; (i < obj->nNodes) && (errNum == WLZ_ERR_NONE); i++, dptr++)
    {
      if(!putfloat((float) dptr->displacement.vtX, fp) ||
         !putfloat((float) dptr->displacement.vtY, fp) ||
	 !putfloat((float) dptr->displacement.vtZ, fp))
      {
	errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
    }
  }
  */
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzIO
* \brief	Writes a 2D mesh transform to the given file.
* \param	fp			Given file.
* \param	obj			Mesh transform.
*/
WlzErrorNum    WlzWriteMeshTransform2D(
				  FILE *fp,
			          WlzMeshTransform *obj)
{
  int		i,
  		j;
  WlzMeshNode	*dptr;
  WlzMeshElem	*eptr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(     !putword(obj->nElem,  fp) ||
	  !putword(obj->nNodes ,fp) )
  {
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }
  else
  {
    /* Write out nodal position and displacement */
    dptr = obj->nodes;
    for(i = 0; (i < obj->nNodes) && (errNum == WLZ_ERR_NONE); i++, dptr++)
    {  /* position */
      if(!putfloat((float) dptr->position.vtX, fp) ||
         !putfloat((float) dptr->position.vtY, fp))
	{
	errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
      /* displacement */
      if(!putfloat((float) dptr->displacement.vtX, fp) ||
         !putfloat((float) dptr->displacement.vtY, fp))
	 {
	errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }

    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Write out elements */
    eptr = obj->elements;
    for(i = 0; (i < obj->nElem) && (errNum == WLZ_ERR_NONE); i++, eptr++)
    {
      /* output the index of this element */
      if(!putword(eptr->idx, fp))
      {
        errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
      /* output the neighbour flags */
      if(!putword(eptr->flags, fp))
      {
        errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
      
      /* output nodes indeces */
      for(j = 0; (j < 3) && (errNum == WLZ_ERR_NONE); j++)
      {
	if (!putword(eptr->nodes[j], fp))
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
      }
      /* output its neighbours */

      for(j = 0; (j < 3) && (errNum == WLZ_ERR_NONE); j++)
      {
	if (!putword(eptr->neighbours[j], fp))
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
      }
    }
  }
    /* now Write out displacement */

  if(errNum == WLZ_ERR_NONE)
  {
    dptr = obj->nodes;
    for(i = 0; (i < obj->nNodes) && (errNum == WLZ_ERR_NONE); i++, dptr++)
    {
      if(!putfloat((float) dptr->displacement.vtX, fp) ||
         !putfloat((float) dptr->displacement.vtY, fp))
	 {
	errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
    }
  }

  return(errNum);
}

