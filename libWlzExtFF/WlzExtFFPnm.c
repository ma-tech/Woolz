#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzExtFFPnm_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlzExtFF/WlzExtFFPnm.c
* \author       Bill Hill
* \date         March 1999
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
* \brief	Functions for reading and writting Woolz objects to and
*		from the portable anymap '.pnm' data format.
* \ingroup	WlzExtFF
*/

#include <ctype.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the given file(s) using the '.pnm'
* 		file format.
* \param	gvnFileName		Given file name.
* \param	dstErr			Destination error number ptr,
*					may be NULL.
*/
WlzObject	*WlzEffReadObjPnm(const char *gvnFileName, WlzErrorNum *dstErr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*obj = NULL;

  obj = WlzEffReadObjStack(gvnFileName, WLZEFF_FORMAT_PNM, &errNum);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz object to the given file(s)
*		using the '.pnm' file format.
* \param	gvnFileName		Given file name.
* \param	obj			Destination error number ptr,
*					may be NULL.
*/
WlzErrorNum	WlzEffWriteObjPnm(const char *gvnFileName, WlzObject *obj)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  errNum = WlzEffWriteObjStack(gvnFileName, WLZEFF_FORMAT_PNM, obj);
  return(errNum);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given 2D Woolz object to the given file
*		using the '.pnm' file format.
* \param	fNameStr		Given file name.
* \param	obj			Given woolz object.
* \param	imgSz			Required image size.
* \param	imgOrg			Required image origin.
* \param	data			Buffer of imgSz bytes.
* \param	bgd			Background value.
*/
WlzErrorNum 	WlzEffWriteObjPnm2D(const char *fNameStr, WlzObject *obj,
				    WlzIVertex2 imgSz, WlzIVertex2 imgOrg,
				    unsigned char *data,
				    unsigned char bgd)
{
  int		tI0;
  FILE		*fP = NULL;
  WlzIBox2	cutBox;
  WlzObject	*cutObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((fP = fopen(fNameStr, "w")) == NULL)
  {
    errNum = WLZ_ERR_WRITE_EOF;
  }
  else
  {
#ifdef _WIN32
  if(fP != NULL)
  {
    if(_setmode(_fileno(fP), 0x8000) == -1)
    {
      errNum = WLZ_ERR_READ_EOF;
    }
  }
#endif
    if((obj == NULL) || (obj->values.core == NULL))
    {
      (void )memset(data, (unsigned int )bgd, imgSz.vtX * imgSz.vtY);
    }
    else
    {
      cutBox.xMin = imgOrg.vtX;
      cutBox.yMin = imgOrg.vtY;
      cutBox.xMax = imgOrg.vtX + imgSz.vtX - 1;
      cutBox.yMax = imgOrg.vtY + imgSz.vtY - 1;
      cutObj = WlzCutObjToValBox2D(obj, cutBox, WLZ_GREY_UBYTE, data,
      			        0, 0.0, 0.0, &errNum);
      if(cutObj)
      {
	if(cutObj->type == WLZ_2D_DOMAINOBJ)
	{
	  /* (void )WlzPopFreePtr(cutObj->values.r->freeptr, NULL, NULL); */
	  AlcFree(cutObj->values.r->freeptr);
	  cutObj->values.r->freeptr = NULL;
	  cutObj->values.r->values.inp = NULL;
	}
	WlzFreeObj(cutObj);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tI0 = imgSz.vtX * imgSz.vtY;
    if((fprintf(fP, "%s\n%d %d\n255\n",
    		WLZEFF_PGM_MAGIC, imgSz.vtX, imgSz.vtY) < 8) ||
       (fwrite(data, sizeof(unsigned char), tI0, fP) !=  tI0))

    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  if(fP)
  {
    fclose(fP);
  }
  return(errNum);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Reads the  data from a pgm file into the data array given.
* \param	fP			Given file stream.
* \param	gvnImgSz		Dst ptr for image size which is
*					assumed valid if height and width are
*					non zero.
* \param	data			Ptr to 2D array for data.
*/
WlzErrorNum 	WlzEffReadObjPnmData2D(FILE *fP, WlzIVertex2 *gvnImgSz,
				       unsigned char ***data)
{
  int		tI0,
  		parsedFlag = 0,
  		recIdx = 0;
  char		*tCP0,
  		*rec;
  WlzIVertex2	imgSz;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char		fRecord[256];

  while((parsedFlag == 0) && (errNum == WLZ_ERR_NONE))
  {
    if(fgets(fRecord, 256, fP) == NULL)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    else
    {
      fRecord[255] = '\0';
      if((tCP0 = strrchr(fRecord,  '\n')) != NULL)
      {
	*tCP0 = '\0';
      }
      rec = fRecord;
      while(*rec && isspace(*rec))
      {
	++rec;
      }
      if(*rec && (*rec != '#'))
      {
	switch(recIdx)
	{
	  case 0:
	    if(strcmp(rec, WLZEFF_PGM_MAGIC))
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	    break;
	  case 1:
	    if(sscanf(rec, "%d %d", &(imgSz.vtX), &(imgSz.vtY)) != 2)
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	    else
	    {
	      if((gvnImgSz->vtX && (imgSz.vtX != gvnImgSz->vtX)) ||
		 (gvnImgSz->vtY && (imgSz.vtY != gvnImgSz->vtY)))
	      {
		errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      else
	      {
	        *gvnImgSz = imgSz;
	      }
	    }
	    break;
	  case 2:
	    if((sscanf(rec, "%d", &tI0) != 1) || (tI0 != 255))
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	    parsedFlag = 1;
	    break;
	}
	++recIdx;
      }
    }
  }
  if(*data == NULL)
  {
    if((AlcUnchar2Malloc(data, imgSz.vtY, imgSz.vtX) != ALC_ER_NONE) ||
       (*data == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tI0 = gvnImgSz->vtX * gvnImgSz->vtY;
    if(fread(**data, sizeof(unsigned char), tI0, fP) != tI0)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  return(errNum);
}
