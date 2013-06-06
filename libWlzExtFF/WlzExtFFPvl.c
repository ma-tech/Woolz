#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzExtFFPvl_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzExtFFPvl.c
* \author       Bill Hill
* \date         February 2013
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2013],
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
* 		the Drishti pvl.nc file format.
* \ingroup	WlzExtFF
*/

#include <errno.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Builds the drishti pvl file name from the given file name
* 		and allocates space for the data file name buffer.
*		These strings should be free'd using AlcFree() when
*		no longer required.
* \param	pvlFileName		Destination pointer for the '.pvl.nc'
* 					file name.
* \param	datFileName		Destination pointer witrh space for
* 					the '.pvl.nc.001' data file name.
* \param	gvnFileName		Given file name with or without a
* 					.pvl extension.
* \param	fCntLen			Number of additional (file number)
* 					characters to allocate space for
* 					in the data file name buffer.
*/
static WlzErrorNum WlzEffPvlFileNames(char **pvlFileName,
				      char **datFileName,
				      const char *gvnFileName,
				      int fCntLen)
{
  int		len;
  char		*fileBody = NULL;
  WlzErrorNum	errFlag = WLZ_ERR_MEM_ALLOC;


  if(++fCntLen < 4)
  {
    fCntLen = 4;
  }
  len = ((int )strlen(gvnFileName) + 7) * sizeof(char);
  if(((fileBody = (char *)AlcMalloc(len)) != NULL) &&
     ((*pvlFileName = (char *)AlcMalloc(len)) != NULL) &&
     ((*datFileName = (char *)AlcMalloc(len + fCntLen)) != NULL))
  {
    char *e;

    (void )strcpy(fileBody, gvnFileName);
    if((e = strstr(fileBody, ".pvl.nc")) != NULL)
    {
      *e = '\0';
    }
    (void )sprintf(*pvlFileName, "%s.pvl.nc", fileBody);
    errFlag = WLZ_ERR_NONE;
  }
  AlcFree(fileBody);
  return(errFlag);
}

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from Drishti pvl.nc files.
* \param	gvnFileName		Given file name.
* \param	dstErr			Destination error number ptr,
*					may be NULL.
*/
WlzObject	*WlzEffReadObjPvl(const char *gvnFileName, WlzErrorNum *dstErr)
{
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_UNIMPLEMENTED;

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz object to a pair of files
*		using the Derishti .pvl.nc file format.
* \param	gvnFileName		Given file name with .ics, .ids or
*					no extension.
* \param	obj			Given woolz object.
*/
WlzErrorNum	WlzEffWriteObjPvl(const char *gvnFileName, WlzObject *obj)
{
  FILE		*fP = NULL;
  unsigned char	 ***data = NULL;
  WlzDVertex3	voxSz;
  WlzIVertex3	objOrg,
  		objSz;
  char		*pvlFileName = NULL,
		*datFileName = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((gvnFileName == NULL) || (*gvnFileName == '\0'))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->type != WLZ_3D_DOMAINOBJ)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(obj->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(obj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(obj->values.core->type != WLZ_VOXELVALUETABLE_GREY)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    errNum = WlzEffPvlFileNames(&pvlFileName, &datFileName,
    			        gvnFileName, 3);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    objOrg.vtX = obj->domain.p->kol1;
    objOrg.vtY = obj->domain.p->line1;
    objOrg.vtZ = obj->domain.p->plane1;
    objSz.vtX = obj->domain.p->lastkl - objOrg.vtX + 1;
    objSz.vtY = obj->domain.p->lastln - objOrg.vtY + 1;
    objSz.vtZ = obj->domain.p->lastpl - objOrg.vtZ + 1;
    voxSz.vtX = obj->domain.p->voxel_size[0];
    voxSz.vtY = obj->domain.p->voxel_size[1];
    voxSz.vtZ = obj->domain.p->voxel_size[2];
    if((objSz.vtX <= 0) || (objSz.vtY <= 0) || (objSz.vtZ <= 0))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
    else
    {
      errNum = WlzToArray3D((void ****)&data, obj, objSz, objOrg,
			    0, WLZ_GREY_UBYTE);
    }
    if((errNum == WLZ_ERR_NONE) && pvlFileName)
    {
      if((fP = fopen(pvlFileName, "w")) == NULL)
      {
        errNum = WLZ_ERR_WRITE_EOF;
      }
      else
      {
	if(fprintf(fP,
                   "<!DOCTYPE Drishti_Header>\n"
		   "<PvlDotNcFileHeader>\n"
		   "  <rawfile></rawfile>\n"
		   "  <voxeltype>unsigned char</voxeltype>\n"
		   "  <pvlvoxeltype>unsigned char</pvlvoxeltype>\n"
		   "  <gridsize>%d %d %d</gridsize>\n"
		   "  <voxelunit>no units</voxelunit>\n"
		   "  <voxelsize>%g %g %g</voxelsize>\n"
		   "  <description>%s exported from Woolz.</description>\n"
		   "  <slabsize>%d</slabsize>\n"
		   "  <rawmap>0 255</rawmap>\n"
		   "  <pvlmap>0 255</pvlmap>\n"
		   "</PvlDotNcFileHeader>\n",
		   objSz.vtZ, objSz.vtY, objSz.vtX,
		   voxSz.vtX, voxSz.vtY, voxSz.vtZ,
		   pvlFileName,
		   objSz.vtX) < 1)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
	(void )fclose(fP);
      }
    }
  }
  if((errNum == WLZ_ERR_NONE) && datFileName)
  {
    size_t	dSz;

    (void )sprintf(datFileName, "%s.001", pvlFileName);
    if((fP = fopen(datFileName, "w")) == NULL)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
    else
    {
      WlzGreyP 	vp;
      char 	header[16];

#ifdef _WIN32
  if (fP != NULL){
	if(_setmode(_fileno(fP), 0x8000) == -1)
	{
		errNum = WLZ_ERR_READ_EOF;
	}
  }
#endif
      dSz = objSz.vtX * objSz.vtY * objSz.vtZ;
      header[0] = 0;
      vp.v = &(header[1]); *(vp.inp) = objSz.vtX;
      vp.v = &(header[5]); *(vp.inp) = objSz.vtY;
      vp.v = &(header[9]); *(vp.inp) = objSz.vtZ;
      if((fwrite(header, sizeof(unsigned char), 13, fP) != 13) ||
         (fwrite(**data, sizeof(unsigned char), dSz, fP) != dSz))
      {
	errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
      (void )fclose(fP);
    }
  }
  AlcFree(pvlFileName);
  AlcFree(datFileName);
  (void )AlcUnchar3Free(data);
  return(errNum);
}
