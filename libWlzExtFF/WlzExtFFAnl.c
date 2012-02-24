#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzExtFFAnl_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlzExtFF/WlzExtFFAnl.c
* \author       Bill Hill
* \date         February 2005
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
* \brief	Functions for reading and writting Woolz objects to
* 		the ANALYZE 7.5 format file.
* \ingroup	WlzExtFF
*/

#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

static WlzErrorNum		WlzEffReadAnlHdrKey(
				  FILE *fP,
				  WlzEffAnlDsr *dsr,
				  int *swap);
static WlzErrorNum		WlzEffReadAnlHdrImageDim(
				  FILE *fP, WlzEffAnlDsr *dsr,
				  int swap);
static WlzErrorNum		WlzEffReadAnlHdrDataHistory(
				  FILE *fP,
				WlzEffAnlDsr *dsr,
				int swap);
static WlzErrorNum		WlzEffReadAnlImgData(
				  WlzObject *obj,
  				WlzGreyType gType,
				WlzEffAnlDsr *dsr,
				FILE *fP,
				int swap);
static WlzErrorNum		WlzEffWriteAnlHdrDataHistory(
				  FILE *fP,
				  WlzEffAnlDsr *dsr);
static WlzErrorNum		WlzEffWriteAnlHdrImageDim(
				  FILE *fP,
				  WlzEffAnlDsr *dsr);
static WlzErrorNum		WlzEffWriteAnlHdrKey(
				  FILE *fP,
				  WlzEffAnlDsr *dsr);
static WlzErrorNum		WlzEffReadAnlChar(
				  FILE *fP,
				  int cnt,
				  WlzUByte *buf);
static WlzErrorNum		WlzEffReadAnlShort(
				  FILE *fP,
				  int cnt,
				  short *buf,
				  int swap);
static WlzErrorNum		WlzEffReadAnlInt(
				  FILE *fP,
				  int cnt,
				  int *buf,
				  int swap);
static WlzErrorNum		WlzEffReadAnlFloat(
				  FILE *fP,
				  int cnt,
				  float *buf,
				  int swap);
static WlzErrorNum		WlzEffReadAnlDouble(
				  FILE *fP,
				  int cnt,
				  double *buf,
				  int swap);
static WlzErrorNum		WlzEffWriteAnlChar(
				  FILE *fP,
				  int cnt,
				  WlzUByte *buf);
static WlzErrorNum		WlzEffWriteAnlShort(
				  FILE *fP,
				  int cnt,
				  short *buf);
static WlzErrorNum		WlzEffWriteAnlInt(
				  FILE *fP,
				  int cnt,
				  int *buf);
static WlzErrorNum		WlzEffWriteAnlFloat(
				  FILE *fP,
				  int cnt,
				  float *buf);
static WlzErrorNum		WlzEffWriteAnlDouble(
				  FILE *fP,
				  int cnt,
				  double *buf);
static void			WlzEffReadAnlSwap2(
				  void *buf,
				  unsigned int cnt);
static void			WlzEffReadAnlSwap4(
				  void *buf,
				  unsigned int cnt);

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Builds the ANALYZE file names from the given file name.
*		These strings should be free'd using AlcFree() when
*		no longer required.
* \param	fileBody		Dest ptr for the file body.
* \param	hdrFileName		Dest ptr for the '.hdr' file
*					name.
* \param	imgFileName		Dest ptr for the '.img' file
*					name.
* \param	gvnFileName		Given file name with '.hdr',
*					'.img' or no extension.
*/
WlzErrorNum	WlzEffAnlFileNames(char **fileBody,
				   char **hdrFileName,
				   char **imgFileName,
				   const char *gvnFileName)
{
  int		tI0;
  WlzErrorNum	errFlag = WLZ_ERR_MEM_ALLOC;

  tI0 = ((int )strlen(gvnFileName) + 5) * sizeof(char);
  if(((*fileBody = (char *)AlcMalloc(tI0)) != NULL) &&
     ((*hdrFileName = (char *)AlcMalloc(tI0)) != NULL) &&
     ((*imgFileName = (char *)AlcMalloc(tI0)) != NULL))
  {
    (void )strcpy(*fileBody, gvnFileName);
    if((tI0 = (int )strlen(*fileBody) - 4) >= 0)
    {
      if((strcmp(*fileBody +  tI0, ".hdr") == 0) ||
	 (strcmp(*fileBody +  tI0, ".img") == 0))
      {
	*(*fileBody +  tI0) = '\0';
      }
    }
    (void )sprintf(*hdrFileName, "%s.hdr", *fileBody);
    (void )sprintf(*imgFileName, "%s.img", *fileBody);
    errFlag = WLZ_ERR_NONE;
  }
  return(errFlag);
}

/*!
* \return	New Woolz object or NULL on error.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the given stream using the
*		ANALYZE 7.5 format. The given file name is used to
*		generate the '.hdr' and '.img' filenames.
* \param	gvnFileName		Given file name.
* \param	dstErr			Destination error code, may be NULL.
*/
WlzObject	*WlzEffReadObjAnl(const char *gvnFileName, WlzErrorNum *dstErr)
{
  int		idx,
  		timePoints = 0;
  FILE		*fP = NULL;
  char		*fileName = NULL,
  		*hdrFileName = NULL,
		*imgFileName = NULL;
  WlzVertex	org,
  		sz;
  WlzEffAnlDsr 	dsr;
  int 		swap;
  WlzObject	*obj = NULL;
  WlzGreyType	gType = WLZ_GREY_ERROR;
  WlzPixelV	bgdV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Compute .hdr and .img file names. */
  if((gvnFileName == NULL) || (*gvnFileName == '\0'))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    errNum = WlzEffAnlFileNames(&fileName, &hdrFileName, &imgFileName,
				gvnFileName);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((hdrFileName == NULL) || (*hdrFileName == '\0') ||
       ((fP = fopen(hdrFileName, "r")) == NULL))
    {
      errNum = WLZ_ERR_READ_EOF;
    }
#ifdef _WIN32
    else
    {
      if(fP != NULL)
      {
	if(_setmode(_fileno(fP), 0x8000) == -1)
	{
	  errNum = WLZ_ERR_READ_EOF;
	}
      }
    }
#endif
  }
  /* Read the .hdr (header) file. */
  if(errNum == WLZ_ERR_NONE)
  {
    (void )memset((void *)&dsr, 0, sizeof(WlzEffAnlDsr));
    errNum = WlzEffReadAnlHdrKey(fP, &dsr, &swap);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffReadAnlHdrImageDim(fP, &dsr, swap);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffReadAnlHdrDataHistory(fP, &dsr, swap);
  }
  if(fP)
  {
    (void )fclose(fP);
    fP = NULL;
  }
  /* Check header data is valid. */
  if(errNum == WLZ_ERR_NONE)
  {
    if('r' != dsr.hk.regular)
    {
      errNum = WLZ_ERR_OBJECT_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(dsr.dim.bitPix)
    {
      case 8:  /* FALLTHROUGH */
      case 16: /* FALLTHROUGH */
      case 32: /* FALLTHROUGH */
      case 64:
        break;
      case 1:
      default:
        errNum = WLZ_ERR_OBJECT_DATA;
	break;
    }
  }
  /* Compute the image parameters. */
  if(errNum == WLZ_ERR_NONE)
  {
    switch(dsr.dim.dataType)
    {
      case WLZEFF_ANL_DT_UNSIGNED_CHAR:
        gType = WLZ_GREY_UBYTE;
	bgdV.type = gType;
	bgdV.v.ubv = 0;
	break;
      case WLZEFF_ANL_DT_SIGNED_SHORT:
        gType = WLZ_GREY_SHORT;
	bgdV.type = gType;
	bgdV.v.shv = 0;
	break;
      case WLZEFF_ANL_DT_SIGNED_INT:
        gType = WLZ_GREY_INT;
	bgdV.type = gType;
	bgdV.v.inv = 0;
	break;
      case WLZEFF_ANL_DT_FLOAT:
        gType = WLZ_GREY_FLOAT;
	bgdV.type = gType;
	bgdV.v.flv = 0.0;
	break;
      case WLZEFF_ANL_DT_DOUBLE:
        gType = WLZ_GREY_DOUBLE;
	bgdV.type = gType;
	bgdV.v.dbv = 0.0;
	break;
      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(dsr.dim.dim[0])
    {
      case 2:
	org.i2.vtX = 0;
	org.i2.vtY = 0;
	sz.i2.vtX = dsr.dim.dim[1];
	sz.i2.vtY = dsr.dim.dim[2];
	timePoints = 1;
	if((sz.i2.vtX < 1) || (sz.i2.vtY < 1))
	{
	  errNum = WLZ_ERR_OBJECT_DATA;
	}
	break;
      case 3:
        org.i3.vtX = 0;
	org.i3.vtY = 0;
	org.i3.vtZ = 0;
	sz.i3.vtX = dsr.dim.dim[1];
	sz.i3.vtY = dsr.dim.dim[2];
	sz.i3.vtZ = dsr.dim.dim[3];
	if((sz.i3.vtX < 1) || (sz.i3.vtY < 1) || (sz.i3.vtZ < 1))
	{
	  errNum = WLZ_ERR_OBJECT_DATA;
	}
        break;
      case 4:
        org.i3.vtX = 0;
	org.i3.vtY = 0;
	org.i3.vtZ = 0;
	sz.i3.vtX = dsr.dim.dim[1];
	sz.i3.vtY = dsr.dim.dim[2];
	sz.i3.vtZ = dsr.dim.dim[3];
	timePoints = dsr.dim.dim[4];
	if((sz.i3.vtX < 1) || (sz.i3.vtY < 1) || (sz.i3.vtZ < 1) ||
	   (timePoints < 1))
	{
	  errNum = WLZ_ERR_OBJECT_DATA;
	}
        break;
      default:
	errNum = WLZ_ERR_OBJECT_DATA;
        break;
        
    }
  }
  /* Create object. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(timePoints > 1)
    {
      obj = (WlzObject *)
            WlzMakeCompoundArray(WLZ_3D_DOMAINOBJ, 1, timePoints, NULL,
	                         WLZ_3D_DOMAINOBJ, &errNum);
      if(obj)
      {
        idx = 0;
	do
	{
	  ((WlzCompoundArray *)obj)->o[idx] = WlzAssignObject(
		WlzMakeCuboid(org.i3.vtZ, org.i3.vtZ + sz.i3.vtZ - 1,
			      org.i3.vtY, org.i3.vtY + sz.i3.vtY - 1,
			      org.i3.vtX, org.i3.vtX + sz.i3.vtX - 1,
			      gType, bgdV, NULL, NULL, &errNum), NULL);
	} while((++idx < timePoints) && (WLZ_ERR_NONE == errNum));
      }
    }
    else if(3 <= dsr.dim.dim[0])
    {
      obj = WlzMakeCuboid(org.i3.vtZ, org.i3.vtZ + sz.i3.vtZ - 1,
                          org.i3.vtY, org.i3.vtY + sz.i3.vtY - 1,
			  org.i3.vtX, org.i3.vtX + sz.i3.vtX - 1,
			  gType, bgdV, NULL, NULL, &errNum);
    }
    else /* 2 == dsr.dim.dim[0] */
    {
      obj = WlzMakeRect(org.i2.vtY, org.i2.vtY + sz.i2.vtY - 1,
			org.i2.vtX, org.i2.vtX + sz.i2.vtX - 1,
			gType, NULL, bgdV, NULL, NULL, &errNum);
    }
  }
  /* Open the .img (image data) file. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((imgFileName == NULL) || (*imgFileName == '\0') ||
       ((fP = fopen(imgFileName, "r")) == NULL))
    {
      errNum = WLZ_ERR_READ_EOF;
    }
#ifdef _WIN32
    else
    {
      if(fP != NULL)
      {
	if(_setmode(_fileno(fP), 0x8000) == -1)
	{
	  errNum = WLZ_ERR_READ_EOF;
	}
      }
    }
#endif
  }
  /* Read the object values. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffReadAnlImgData(obj, gType, &dsr, fP, swap);
  }
  if(fP)
  {
    (void )fclose(fP);
  }
  (void )AlcFree(fileName);
  (void )AlcFree(hdrFileName);
  (void )AlcFree(imgFileName);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(obj);
    obj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz object to the given file(s)
*		using the ANALYZE 7.5 file format. The given file name is
*		used to generate the '.hdr' and '.img' filenames.
* \param	gvnFileName		Given file name with .hdr, .img or no
* 					extension.
* \param	obj			Given woolz object.
*/
WlzErrorNum	WlzEffWriteObjAnl(const char *gvnFileName, WlzObject *obj)
{
  int  		cnt;
  double	dMin,
  		dMax;
  FILE		*fP = NULL;
  void	 	*data = NULL;
  WlzGreyP	dataP;
  WlzEffAnlDsr 	dsr;
  WlzGreyType	gType;
  WlzVertex	org,
  		sz;
  char		*fileName = NULL,
  		*hdrFileName = NULL,
		*imgFileName = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((gvnFileName == NULL) || (*gvnFileName == '\0'))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(obj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    errNum = WlzEffAnlFileNames(&fileName, &hdrFileName, &imgFileName,
				gvnFileName);
  }
  /* Fill in the ANALYZE file header data structures. */
  if(errNum == WLZ_ERR_NONE)
  {
    (void )memset(&dsr, 0, sizeof(WlzEffAnlDsr));
    switch(obj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	org.i2.vtX = obj->domain.i->kol1;
	org.i2.vtY = obj->domain.i->line1;
	sz.i2.vtX = obj->domain.i->lastkl - org.i2.vtX + 1;
	sz.i2.vtY = obj->domain.i->lastln - org.i2.vtY + 1;
	if((sz.i2.vtX <= 0) || (sz.i2.vtY <= 0))
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	}
	else
	{
	  dsr.dim.dim[0] = 2;
	  dsr.dim.dim[1] = sz.i2.vtX;
	  dsr.dim.dim[2] = sz.i2.vtY;
	  dsr.dim.pixdim[0] = 2;
	  dsr.dim.pixdim[1] = 1.0;
	  dsr.dim.pixdim[2] = 1.0;
	}
        break;
      case WLZ_3D_DOMAINOBJ:
	org.i3.vtX = obj->domain.p->kol1;
	org.i3.vtY = obj->domain.p->line1;
	org.i3.vtZ = obj->domain.p->plane1;
	sz.i3.vtX = obj->domain.p->lastkl - org.i3.vtX + 1;
	sz.i3.vtY = obj->domain.p->lastln - org.i3.vtY + 1;
	sz.i3.vtZ = obj->domain.p->lastpl - org.i3.vtZ + 1;
	if((sz.i3.vtX <= 0) || (sz.i3.vtY <= 0) || (sz.i3.vtZ <= 0))
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	}
        else
	{
	  dsr.dim.dim[0] = 3;
	  dsr.dim.dim[1] = sz.i3.vtX;
	  dsr.dim.dim[2] = sz.i3.vtY;
	  dsr.dim.dim[3] = sz.i3.vtZ;
	  dsr.dim.pixdim[0] = 3;
	  dsr.dim.pixdim[1] = obj->domain.p->voxel_size[0];
	  dsr.dim.pixdim[2] = obj->domain.p->voxel_size[1];
	  dsr.dim.pixdim[3] = obj->domain.p->voxel_size[2];
	}
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    (void )WlzGreyStats(obj, &gType, &dMin, &dMax, NULL, NULL, NULL, NULL,
    			&errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dsr.dim.glMax = (int )ceil(dMax);
    dsr.dim.glMin = (int )floor(dMin);
    switch(gType)
    {
      case WLZ_GREY_INT:
        dsr.dim.dataType = WLZEFF_ANL_DT_SIGNED_INT;
        dsr.dim.bitPix = 32;
        break;
      case WLZ_GREY_SHORT:
        dsr.dim.dataType = WLZEFF_ANL_DT_SIGNED_SHORT;
        dsr.dim.bitPix = 16;
	break;
      case WLZ_GREY_UBYTE:
        dsr.dim.dataType = WLZEFF_ANL_DT_UNSIGNED_CHAR;
        dsr.dim.bitPix = 8;
	break;
      case WLZ_GREY_FLOAT:
        dsr.dim.dataType = WLZEFF_ANL_DT_FLOAT;
        dsr.dim.bitPix = 32;
	break;
      case WLZ_GREY_DOUBLE:
        dsr.dim.dataType = WLZEFF_ANL_DT_DOUBLE;
        dsr.dim.bitPix = 64;
	break;
      default:
	errNum = WLZ_ERR_GREY_TYPE;
        break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dsr.hk.hdrSz = 348;
    dsr.hk.regular = 'r';
    dsr.hk.extents = 16384;
    (void )strcpy(dsr.hist.descrip, "Converted from Woolz");
    (void )strcpy(dsr.hist.originator, "Woolz");
    (void )strcpy(dsr.hist.generated, "Woolz");
  }
  /* Get Woolz object's values as an array. */
  if(WLZ_ERR_NONE == errNum)
  {
    switch(obj->type)
    {
      case WLZ_2D_DOMAINOBJ:
        cnt = sz.i2.vtX * sz.i2.vtY;
	errNum = WlzToArray2D((void ***)&data, obj, sz.i2, org.i2, 0, gType);
	if(errNum == WLZ_ERR_NONE)
	{
	  switch(gType)
	  {
	    case WLZ_GREY_INT:
	      dataP.inp = *(int **)data;
	      break;
	    case WLZ_GREY_SHORT:
	      dataP.shp = *(short **)data;
	      break;
	    case WLZ_GREY_UBYTE:
	      dataP.ubp = *(WlzUByte **)data;
	      break;
	    case WLZ_GREY_FLOAT:
	      dataP.flp = *(float **)data;
	      break;
	    case WLZ_GREY_DOUBLE:
	      dataP.dbp = *(double **)data;
	      break;
	    default:
	      errNum = WLZ_ERR_GREY_TYPE;
	      break;
	  }
	}
        break;
      case WLZ_3D_DOMAINOBJ:
        cnt = sz.i3.vtX * sz.i3.vtY * sz.i3.vtZ;
	errNum = WlzToArray3D((void ****)&data, obj, sz.i3, org.i3, 0, gType);
	if(errNum == WLZ_ERR_NONE)
	{
	  switch(gType)
	  {
	    case WLZ_GREY_INT:
	      dataP.inp = **(int ***)data;
	      break;
	    case WLZ_GREY_SHORT:
	      dataP.shp = **(short ***)data;
	      break;
	    case WLZ_GREY_UBYTE:
	      dataP.ubp = **(WlzUByte ***)data;
	      break;
	    case WLZ_GREY_FLOAT:
	      dataP.flp = **(float ***)data;
	      break;
	    case WLZ_GREY_DOUBLE:
	      dataP.dbp = **(double ***)data;
	      break;
	    default:
	      errNum = WLZ_ERR_GREY_TYPE;
	      break;
	  }
	}
        break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
        break;
    }
  }
  /* Open the ANALYZE .hdr file and write the header information. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((fP = fopen(hdrFileName, "w")) == NULL)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
#ifdef _WIN32
    else if(_setmode(_fileno(fP), 0x8000) == -1)
    {
	errNum = WLZ_ERR_READ_EOF;
    }
#endif
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffWriteAnlHdrKey(fP, &dsr);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffWriteAnlHdrImageDim(fP, &dsr);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffWriteAnlHdrDataHistory(fP, &dsr);
  }
  if(fP)
  {
    (void )fclose(fP);
    fP = NULL;
  }
  /* Open the ANALYZE .img file and write the image data. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((fP = fopen(imgFileName, "w")) == NULL)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
#ifdef _WIN32
    else if(_setmode(_fileno(fP), 0x8000) == -1)
    {
	errNum = WLZ_ERR_READ_EOF;
    }
#endif
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(gType)
    {
      case WLZ_GREY_INT:
	errNum = WlzEffWriteAnlInt(fP, cnt, dataP.inp);
        break;
      case WLZ_GREY_SHORT:
	errNum = WlzEffWriteAnlShort(fP, cnt, dataP.shp);
        break;
      case WLZ_GREY_UBYTE:
	errNum = WlzEffWriteAnlChar(fP, cnt, dataP.ubp);
        break;
      case WLZ_GREY_FLOAT:
	errNum = WlzEffWriteAnlFloat(fP, cnt, dataP.flp);
        break;
      case WLZ_GREY_DOUBLE:
	errNum = WlzEffWriteAnlDouble(fP, cnt, dataP.dbp);
	break;
      default:
        break;
    }
  }
  if(fP)
  {
    (void )fclose(fP);
  }
  if(data)
  {
    switch(obj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	Alc2Free((void **)data);
	break;
      case WLZ_3D_DOMAINOBJ:
	Alc3Free((void ***)data);
	break;
      default:
	break;
    }
  }
  (void )AlcFree(fileName);
  (void )AlcFree(hdrFileName);
  (void )AlcFree(imgFileName);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads an ANALYZE header key.
* \param	fP		Open file stream.
* \param	dsr		File header data structure to be filled in.
* \param	swap		Destination pointer for the swap value,
				non-zero if byte swaping required.
*/
static WlzErrorNum	WlzEffReadAnlHdrKey(FILE *fP, WlzEffAnlDsr *dsr,
					    int *swap)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(WlzEffReadAnlInt(fP, 1, &(dsr->hk.hdrSz), 0))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else
  {
    switch(dsr->hk.hdrSz)
    {
      case 348:
        *swap = 0;
	break;
      case 1543569408:
        *swap = 1;
	break;
      default:
        errNum = WLZ_ERR_OBJECT_DATA;
	break;
    }
  }
  if(WlzEffReadAnlChar(fP, 10, (WlzUByte *)(dsr->hk.dataType)) ||
     WlzEffReadAnlChar(fP, 18, (WlzUByte *)(dsr->hk.dbName)) ||
     WlzEffReadAnlInt(fP, 1, &(dsr->hk.extents), *swap) ||
     WlzEffReadAnlShort(fP, 1, &(dsr->hk.sessionErr), *swap) ||
     WlzEffReadAnlChar(fP, 1, (WlzUByte *)&(dsr->hk.regular)) ||
     WlzEffReadAnlChar(fP, 1, (WlzUByte *)&(dsr->hk.hKeyUn0)))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads an ANALYZE header image dimensions.
* \param	fP		Open file stream.
* \param	dsr		File header data structure to be filled in.
* \param	swap		Destination pointer for the swap value.
*/
static WlzErrorNum	WlzEffReadAnlHdrImageDim(FILE *fP, WlzEffAnlDsr *dsr,
				int swap)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(WlzEffReadAnlShort(fP, 8, dsr->dim.dim, swap) ||
     WlzEffReadAnlShort(fP, 1, &(dsr->dim.unused8), swap) ||
     WlzEffReadAnlShort(fP, 1, &(dsr->dim.unused9), swap) ||
     WlzEffReadAnlShort(fP, 1, &(dsr->dim.unused10), swap) ||
     WlzEffReadAnlShort(fP, 1, &(dsr->dim.unused11), swap) ||
     WlzEffReadAnlShort(fP, 1, &(dsr->dim.unused12), swap) ||
     WlzEffReadAnlShort(fP, 1, &(dsr->dim.unused13), swap) ||
     WlzEffReadAnlShort(fP, 1, &(dsr->dim.unused14), swap) ||
     WlzEffReadAnlShort(fP, 1, &(dsr->dim.dataType), swap) ||
     WlzEffReadAnlShort(fP, 1, &(dsr->dim.bitPix), swap) ||
     WlzEffReadAnlShort(fP, 1, &(dsr->dim.dimUn0), swap) ||
     WlzEffReadAnlFloat(fP, 8, dsr->dim.pixdim, swap) ||
     WlzEffReadAnlFloat(fP, 1, &(dsr->dim.voxOffset), swap) ||
     WlzEffReadAnlFloat(fP, 1, &(dsr->dim.fUnused1), swap) ||
     WlzEffReadAnlFloat(fP, 1, &(dsr->dim.fUnused2), swap) ||
     WlzEffReadAnlFloat(fP, 1, &(dsr->dim.fUnused3), swap) ||
     WlzEffReadAnlFloat(fP, 1, &(dsr->dim.calMax), swap) ||
     WlzEffReadAnlFloat(fP, 1, &(dsr->dim.calMin), swap) ||
     WlzEffReadAnlFloat(fP, 1, &(dsr->dim.compressed), swap) ||
     WlzEffReadAnlFloat(fP, 1, &(dsr->dim.verified), swap) ||
     WlzEffReadAnlInt(fP, 1, &(dsr->dim.glMax), swap) ||
     WlzEffReadAnlInt(fP, 1, &(dsr->dim.glMin), swap))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads an ANALYZE header data history.
* \param	fP		Open file stream.
* \param	dsr		File header data structure to be filled in.
* \param	swap		Destination pointer for the swap value.
*/
static WlzErrorNum	WlzEffReadAnlHdrDataHistory(FILE *fP,
				WlzEffAnlDsr *dsr, int swap)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(WlzEffReadAnlChar(fP, 80, (WlzUByte *)(dsr->hist.descrip)) ||
     WlzEffReadAnlChar(fP, 24, (WlzUByte *)(dsr->hist.auxFile)) ||
     WlzEffReadAnlChar(fP, 1, (WlzUByte *)&(dsr->hist.orient)) ||
     WlzEffReadAnlChar(fP, 10, (WlzUByte *)(dsr->hist.originator)) ||
     WlzEffReadAnlChar(fP, 10, (WlzUByte *)(dsr->hist.generated)) ||
     WlzEffReadAnlChar(fP, 10, (WlzUByte *)(dsr->hist.scanNum)) ||
     WlzEffReadAnlChar(fP, 10, (WlzUByte *)(dsr->hist.patientId)) ||
     WlzEffReadAnlChar(fP, 10, (WlzUByte *)(dsr->hist.expDate)) ||
     WlzEffReadAnlChar(fP, 10, (WlzUByte *)(dsr->hist.expTime)) ||
     WlzEffReadAnlChar(fP, 3, (WlzUByte *)(dsr->hist.hisUn0)) ||
     WlzEffReadAnlInt(fP, 1, &(dsr->hist.views), swap) ||
     WlzEffReadAnlInt(fP, 1, &(dsr->hist.volsAdded), swap) ||
     WlzEffReadAnlInt(fP, 1, &(dsr->hist.startField), swap) ||
     WlzEffReadAnlInt(fP, 1, &(dsr->hist.fieldSkip), swap) ||
     WlzEffReadAnlInt(fP, 1, &(dsr->hist.oMax), swap) ||
     WlzEffReadAnlInt(fP, 1, &(dsr->hist.oMin), swap) ||
     WlzEffReadAnlInt(fP, 1, &(dsr->hist.sMax), swap) ||
     WlzEffReadAnlInt(fP, 1, &(dsr->hist.sMin), swap))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads the ANALYZE image data into the given object.
* \param	obj		Object with values to be set from the
*				image file.
* \param	gType		Object and image data grey type.
* \param	dsr		File header data structure to be filled in.
* \param	fP		Open file stream.
* \param 	swap		Swap bytes if non-zero.
*/
static WlzErrorNum	WlzEffReadAnlImgData(WlzObject *obj,
  				WlzGreyType gType, WlzEffAnlDsr *dsr,
				FILE *fP, int swap)
{
  int		cnt,
  		idx;
  WlzDomain	*dom2DP;
  WlzValues	*val2DP;
  WlzObject	*obj2D;
  WlzDomain	dom;
  WlzValues	val;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(obj->type)
  {
    case WLZ_2D_DOMAINOBJ:
      dom.i = obj->domain.i;
      val.r = obj->values.r;
      cnt = (dom.i->lastkl - dom.i->kol1 + 1) *
            (dom.i->lastln - dom.i->line1 + 1);
      switch(gType)
      {
        case WLZ_GREY_INT:
	  errNum = WlzEffReadAnlInt(fP, cnt, val.r->values.inp, swap);
	  break;
	case WLZ_GREY_SHORT:
	  errNum = WlzEffReadAnlShort(fP, cnt, val.r->values.shp, swap);
	  break;
	case WLZ_GREY_UBYTE:
	  errNum = WlzEffReadAnlChar(fP, cnt, val.r->values.ubp);
	  break;
	case WLZ_GREY_FLOAT:
	  errNum = WlzEffReadAnlFloat(fP, cnt, val.r->values.flp, swap);
	  break;
	case WLZ_GREY_DOUBLE:
	  errNum = WlzEffReadAnlDouble(fP, cnt, val.r->values.dbp, swap);
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
      break;
    case WLZ_3D_DOMAINOBJ:
      dom.p = obj->domain.p;
      dom2DP = dom.p->domains;
      val.vox = obj->values.vox;
      val2DP = val.vox->values;
      dom.p->voxel_size[0] = dsr->dim.pixdim[1];
      dom.p->voxel_size[1] = dsr->dim.pixdim[2];
      dom.p->voxel_size[2] = dsr->dim.pixdim[3];
      idx = dom.p->plane1;
      do
      {
        obj2D =  WlzMakeMain(WLZ_2D_DOMAINOBJ, *dom2DP++, *val2DP++,
			     NULL, NULL, &errNum);
	if(WLZ_ERR_NONE == errNum)
	{
	  errNum = WlzEffReadAnlImgData(obj2D, gType, dsr, fP, swap);
	}
        WlzFreeObj(obj2D);
      } while((++idx <= dom.p->lastpl) && (WLZ_ERR_NONE == errNum));
      break;
    case WLZ_COMPOUND_ARR_1:
      idx = 0;
      while((idx < ((WlzCompoundArray *)obj)->n) && (WLZ_ERR_NONE == errNum))
      {
	errNum = WlzEffReadAnlImgData(((WlzCompoundArray *)obj)->o[idx],
	                              gType, dsr, fP, swap);
        ++idx;
      }
      break;
    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Writes an ANALYZE header key.
* \param	fP		Open file stream.
* \param	dsr		File header data structure to output.
*/
static WlzErrorNum	WlzEffWriteAnlHdrKey(FILE *fP, WlzEffAnlDsr *dsr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(WlzEffWriteAnlInt(fP, 1, &(dsr->hk.hdrSz)) ||
     WlzEffWriteAnlChar(fP, 10, (WlzUByte *)(dsr->hk.dataType)) ||
     WlzEffWriteAnlChar(fP, 18, (WlzUByte *)(dsr->hk.dbName)) ||
     WlzEffWriteAnlInt(fP, 1, &(dsr->hk.extents)) ||
     WlzEffWriteAnlShort(fP, 1, &(dsr->hk.sessionErr)) ||
     WlzEffWriteAnlChar(fP, 1, (WlzUByte *)&(dsr->hk.regular)) ||
     WlzEffWriteAnlChar(fP, 1, (WlzUByte *)&(dsr->hk.hKeyUn0)))
  {
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Writes an ANALYZE header image dimensions.
* \param	fP		Open file stream.
* \param	dsr		File header data structure to output.
*/
static WlzErrorNum	WlzEffWriteAnlHdrImageDim(FILE *fP, WlzEffAnlDsr *dsr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(WlzEffWriteAnlShort(fP, 8, dsr->dim.dim) ||
     WlzEffWriteAnlShort(fP, 1, &(dsr->dim.unused8)) ||
     WlzEffWriteAnlShort(fP, 1, &(dsr->dim.unused9)) ||
     WlzEffWriteAnlShort(fP, 1, &(dsr->dim.unused10)) ||
     WlzEffWriteAnlShort(fP, 1, &(dsr->dim.unused11)) ||
     WlzEffWriteAnlShort(fP, 1, &(dsr->dim.unused12)) ||
     WlzEffWriteAnlShort(fP, 1, &(dsr->dim.unused13)) ||
     WlzEffWriteAnlShort(fP, 1, &(dsr->dim.unused14)) ||
     WlzEffWriteAnlShort(fP, 1, &(dsr->dim.dataType)) ||
     WlzEffWriteAnlShort(fP, 1, &(dsr->dim.bitPix)) ||
     WlzEffWriteAnlShort(fP, 1, &(dsr->dim.dimUn0)) ||
     WlzEffWriteAnlFloat(fP, 8, dsr->dim.pixdim) ||
     WlzEffWriteAnlFloat(fP, 1, &(dsr->dim.voxOffset)) ||
     WlzEffWriteAnlFloat(fP, 1, &(dsr->dim.fUnused1)) ||
     WlzEffWriteAnlFloat(fP, 1, &(dsr->dim.fUnused2)) ||
     WlzEffWriteAnlFloat(fP, 1, &(dsr->dim.fUnused3)) ||
     WlzEffWriteAnlFloat(fP, 1, &(dsr->dim.calMax)) ||
     WlzEffWriteAnlFloat(fP, 1, &(dsr->dim.calMin)) ||
     WlzEffWriteAnlFloat(fP, 1, &(dsr->dim.compressed)) ||
     WlzEffWriteAnlFloat(fP, 1, &(dsr->dim.verified)) ||
     WlzEffWriteAnlInt(fP, 1, &(dsr->dim.glMax)) ||
     WlzEffWriteAnlInt(fP, 1, &(dsr->dim.glMin)))
  {
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Writes an ANALYZE header data history.
* \param	fP		Open file stream.
* \param	dsr		File header data structure to be filled in.
*/
static WlzErrorNum	WlzEffWriteAnlHdrDataHistory(FILE *fP,
				WlzEffAnlDsr *dsr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(WlzEffWriteAnlChar(fP, 80, (WlzUByte *)(dsr->hist.descrip)) ||
     WlzEffWriteAnlChar(fP, 24, (WlzUByte *)(dsr->hist.auxFile)) ||
     WlzEffWriteAnlChar(fP, 1, (WlzUByte *)&(dsr->hist.orient)) ||
     WlzEffWriteAnlChar(fP, 10, (WlzUByte *)(dsr->hist.originator)) ||
     WlzEffWriteAnlChar(fP, 10, (WlzUByte *)(dsr->hist.generated)) ||
     WlzEffWriteAnlChar(fP, 10, (WlzUByte *)(dsr->hist.scanNum)) ||
     WlzEffWriteAnlChar(fP, 10, (WlzUByte *)(dsr->hist.patientId)) ||
     WlzEffWriteAnlChar(fP, 10, (WlzUByte *)(dsr->hist.expDate)) ||
     WlzEffWriteAnlChar(fP, 10, (WlzUByte *)(dsr->hist.expTime)) ||
     WlzEffWriteAnlChar(fP, 3, (WlzUByte *)(dsr->hist.hisUn0)) ||
     WlzEffWriteAnlInt(fP, 1, &(dsr->hist.views)) ||
     WlzEffWriteAnlInt(fP, 1, &(dsr->hist.volsAdded)) ||
     WlzEffWriteAnlInt(fP, 1, &(dsr->hist.startField)) ||
     WlzEffWriteAnlInt(fP, 1, &(dsr->hist.fieldSkip)) ||
     WlzEffWriteAnlInt(fP, 1, &(dsr->hist.oMax)) ||
     WlzEffWriteAnlInt(fP, 1, &(dsr->hist.oMin)) ||
     WlzEffWriteAnlInt(fP, 1, &(dsr->hist.sMax)) ||
     WlzEffWriteAnlInt(fP, 1, &(dsr->hist.sMin)))
  {
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads chars into the buffer from the given open stream.
* \param	fP			Given open stream.
* \param	cnt			Number of chars to read.
* \param	buf			Buffer for chars.
*/
static WlzErrorNum	WlzEffReadAnlChar(FILE *fP, int cnt, WlzUByte *buf)
{
  WlzErrorNum	errNum;
  
  errNum = (fread(buf, sizeof(char), cnt, fP) == cnt)?
           WLZ_ERR_NONE: WLZ_ERR_READ_INCOMPLETE;
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads shorts into the buffer from the given open stream.
* \param	fP			Given open stream.
* \param	cnt			Number of shorts to read.
* \param	buf			Buffer for shorts.
* \param 	swap			Swap bytes if non-zero.
*/
static WlzErrorNum	WlzEffReadAnlShort(FILE *fP, int cnt, short *buf,
					int swap)
{
  WlzErrorNum	errNum;
  
  errNum = (fread(buf, sizeof(short), cnt, fP) == cnt)?
           WLZ_ERR_NONE: WLZ_ERR_READ_INCOMPLETE;
  if(swap)
  {
    WlzEffReadAnlSwap2(buf, cnt);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads ints into the buffer from the given open stream.
* \param	fP			Given open stream.
* \param	cnt			Number of ints to read.
* \param	buf			Buffer for ints.
* \param	swap			Swap bytes if non-zero.
*/
static WlzErrorNum	WlzEffReadAnlInt(FILE *fP, int cnt, int *buf,
					int swap)
{
  WlzErrorNum	errNum;
  
  errNum = (fread(buf, sizeof(int), cnt, fP) == cnt)?
           WLZ_ERR_NONE: WLZ_ERR_READ_INCOMPLETE;
  if(swap)
  {
    WlzEffReadAnlSwap4(buf, cnt);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads floats into the buffer from the given open stream.
* \param	fP			Given open stream.
* \param	cnt			Number of floats to read.
* \param	buf			Buffer for floats.
* \param	swap			Swap bytes if non-zero.
*/
static WlzErrorNum	WlzEffReadAnlFloat(FILE *fP, int cnt, float *buf,
					int swap)
{
  WlzErrorNum	errNum;
  
  errNum = (fread(buf, sizeof(float), cnt, fP) == cnt)?
           WLZ_ERR_NONE: WLZ_ERR_READ_INCOMPLETE;
  if(swap)
  {
    WlzEffReadAnlSwap4(buf, cnt);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads doubles into the buffer from the given open stream.
* \param	fP			Given open stream.
* \param	cnt			Number of doubles to read.
* \param	buf			Buffer for doubles.
* \param	swap			Swap bytes if non-zero.
*/
static WlzErrorNum	WlzEffReadAnlDouble(FILE *fP, int cnt, double *buf,
					int swap)
{
  WlzErrorNum	errNum;
  
  errNum = (fread(buf, sizeof(double), cnt, fP) == cnt)?
           WLZ_ERR_NONE: WLZ_ERR_READ_INCOMPLETE;
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Writes chars from the buffer to the given open stream.
* \param	fP			Given open stream.
* \param	cnt			Number of chars to write.
* \param	buf			Buffer with chars.
*/
static WlzErrorNum	WlzEffWriteAnlChar(FILE *fP, int cnt, WlzUByte *buf)
{
  WlzErrorNum	errNum;
  
  errNum = (fwrite(buf, sizeof(char), cnt, fP) == cnt)?
           WLZ_ERR_NONE: WLZ_ERR_WRITE_INCOMPLETE;
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Writes shorts from the buffer to the given open stream.
* \param	fP			Given open stream.
* \param	cnt			Number of shorts to write.
* \param	buf			Buffer with shorts.
*/
static WlzErrorNum	WlzEffWriteAnlShort(FILE *fP, int cnt, short *buf)
{
  WlzErrorNum	errNum;
  
  errNum = (fwrite(buf, sizeof(short), cnt, fP) == cnt)?
           WLZ_ERR_NONE: WLZ_ERR_WRITE_INCOMPLETE;
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Writes ints from the buffer to the given open stream.
* \param	fP			Given open stream.
* \param	cnt			Number of ints to write.
* \param	buf			Buffer with ints.
*/
static WlzErrorNum	WlzEffWriteAnlInt(FILE *fP, int cnt, int *buf)
{
  WlzErrorNum	errNum;
  
  errNum = (fwrite(buf, sizeof(int), cnt, fP) == cnt)?
           WLZ_ERR_NONE: WLZ_ERR_WRITE_INCOMPLETE;
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Writes floats from the buffer to the given open stream.
* \param	fP			Given open stream.
* \param	cnt			Number of floats to write.
* \param	buf			Buffer with floats.
*/
static WlzErrorNum	WlzEffWriteAnlFloat(FILE *fP, int cnt, float *buf)
{
  WlzErrorNum	errNum;
  
  errNum = (fwrite(buf, sizeof(float), cnt, fP) == cnt)?
           WLZ_ERR_NONE: WLZ_ERR_WRITE_INCOMPLETE;
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Writes doubles from the buffer to the given open stream.
* \param	fP			Given open stream.
* \param	cnt			Number of doubles to write.
* \param	buf			Buffer with doubles.
*/
static WlzErrorNum	WlzEffWriteAnlDouble(FILE *fP, int cnt, double *buf)
{
  WlzErrorNum	errNum;
  
  errNum = (fwrite(buf, sizeof(double), cnt, fP) == cnt)?
           WLZ_ERR_NONE: WLZ_ERR_WRITE_INCOMPLETE;
  return(errNum);
}

/*!
* \return	void
* \ingroup	WlzExtFF
* \brief	Swaps bytes in data composed of two byte tuples such
*		as shorts.
* \param	buf			Buffer of data.
* \param	cnt			Number of two byte data.
*/
static void	WlzEffReadAnlSwap2(void *buf, unsigned int cnt)
{
  WlzUByte	ubv;
  WlzUByte	*ubp;

  ubp = (WlzUByte *)buf;
  while(cnt-- > 0)
  {
    ubv = *(ubp + 0);
    *(ubp + 0) = *(ubp + 1);
    *(ubp + 1) = ubv;
    ubp += 2;
  }
}

/*!
* \return	void
* \ingroup	WlzExtFF
* \brief	Swaps bytes in data composed of four byte tuples such
*		as ints and floats.
* \param	buf			Buffer of data.
* \param	cnt			Number of four byte data.
*/
static void	WlzEffReadAnlSwap4(void *buf, unsigned int cnt)
{
  WlzUByte	ubv;
  WlzUByte	*ubp;

  ubp = (WlzUByte *)buf;
  while(cnt-- > 0)
  {
    ubv = *(ubp + 0);
    *(ubp + 0) = *(ubp + 1);
    *(ubp + 1) = ubv;
    ubv = *(ubp + 2);
    *(ubp + 2) = *(ubp + 3);
    *(ubp + 3) = ubv;
    ubp += 4;
  }
}
