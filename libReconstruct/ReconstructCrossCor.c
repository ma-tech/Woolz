#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _ReconstructCrossCor_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libReconstruct/ReconstructCrossCor.c
* \author       Bill Hill
* \date         April 1999
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
* \brief	Provides functions for the automatic registration of a single
* 		pair of serial sections.
* \ingroup	Reconstruct
*/
#include <Reconstruct.h>
#include <unistd.h>
#include <string.h>
#include <float.h>


static RecError	RecCCorObjToFour(double **data, double *sSq, WlzObject *obj,
				 WlzIVertex2 org, WlzIVertex2 size,
				 RecPPControl *ppCtrl);

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Performs the cross correlation of the two given Woolz
*		objects and return their cross correlation together with the
*		fourier transform of the second object. The cross correlation
*		data are un-normalized, this can be done by multiplying each
*		datum by \f$\frac{1}{\sqrt{ss1 ss2}}\f$, where \f$ss1\f$ and
*		\f$ss2\f$ are the sums of squares of the two objects.
* \param	data0			Given object (obj0) data and returned
*					FFT.
* \param	data1			Given object (obj1) data and returned
*					FFT.
* \param	sSq0			Destination for sum of squares of obj0
*					after any preprocessing which may be
*					used together with sSq1 to normalize
*					the cross-correlation data.
* \param	sSq1			Destination for sum of squares of obj1
*					after any preprocessing which may be
*					used together with sSq0 to normalize
*					the cross-correlation data.
* \param	ccFlags			Crosscorrelation control flags.
* \param	obj0			First object for cross-correlation.
* \param	obj1			Second object for cross-correlation.
* \param	org			Origin of data wrt the given woolz
* 					objects, should be less than or equal
* 					to the minimum coordinate pair.
* \param	size			The size of data0 and data1.
* \param	ppCtrl			Pre-processing control data structure.
*/
RecError	RecCrossCorrelate(double **data0, double **data1,
			          double *sSq0, double *sSq1,
			          RecCcFlag ccFlags,
			          WlzObject *obj0, WlzObject *obj1,
			          WlzIVertex2 org, WlzIVertex2 size,
				  RecPPControl *ppCtrl)
{
  int		idX,
		idY;
  RecError	errFlag = REC_ERR_NONE;
  WlzIVertex2	tC2I,
		p2Size;
  double	tD1,
		tD2,
		tD3,
		tD4;
  double	*tDP1,
		*tDP2;

  REC_DBG((REC_DBG_CROSS|REC_DBG_LVL_FN|REC_DBG_LVL_1),
          ("RecCrossCorrelate FE 0x%lx 0x%lx 0x%lx 0x%lx %d 0x%lx 0x%lx"
	   " {%d %d} {%d %d}\n",
	   (unsigned long )data0, (unsigned long )data1,
	   (unsigned long )sSq0, (unsigned long )sSq1, ccFlags,
	   (unsigned long )obj0, (unsigned long )obj1,
	   org.vtX, org.vtY, size.vtX, size.vtY));
  if((data0 == NULL) || (data1 == NULL) || (obj0 == NULL) || (obj1 == NULL))
  {
     errFlag = REC_ERR_FUNC;
  }
  if(errFlag == REC_ERR_NONE)
  {
    p2Size = RecPowerOfTwoC2I(&tC2I, size);
    if((tC2I.vtX != size.vtX) || (tC2I.vtY != size.vtY) ||
       (p2Size.vtX < REC_FOUR_DIMP2_MIN) ||
       (p2Size.vtY < REC_FOUR_DIMP2_MIN) ||
       (p2Size.vtX > REC_FOUR_DIMP2_MAX) ||
       (p2Size.vtY > REC_FOUR_DIMP2_MAX))
    {
      errFlag = REC_ERR_FUNC;
    }
  }
  if((errFlag == REC_ERR_NONE) && ((ccFlags & REC_CCFLAG_DATA0VALID) == 0))
  {
    errFlag = RecCCorObjToFour(data0, sSq0, obj0, org, size, ppCtrl);
  }
  if((errFlag == REC_ERR_NONE) && ((ccFlags & REC_CCFLAG_DATA1VALID) == 0))
  {
    errFlag = RecCCorObjToFour(data1, sSq1, obj1, org, size, ppCtrl);
  }
  if(errFlag == REC_ERR_NONE)
  {
    tC2I.vtX = size.vtX / 2;
    tC2I.vtY = size.vtY / 2;
    for(idY = 0; idY < size.vtY; ++idY)
    {
      tDP1 = *(data0 + idY) + 1;
      tDP2 = *(data1 + idY) + 1;
      for(idX = 1; idX < tC2I.vtX; ++idX)
      {
	tD1 = *tDP1;
	tD2 = *(tDP1 + tC2I.vtX);
	tD3 = *tDP2;
	tD4 = -*(tDP2 + tC2I.vtX);
	*tDP1 = tD1 * tD3 - tD2 * tD4;
	*(tDP1 + tC2I.vtX) = tD1 * tD4 + tD2 * tD3;
	++tDP1;
	++tDP2;
      }
    }
    for(idX = 0; idX < size.vtX; idX += tC2I.vtX)
    {
      for(idY = 1; idY < tC2I.vtY; ++idY)
      {
	tDP1 = *(data0 + idY) + idX;
	tDP2 = *(data0 + tC2I.vtY + idY) + idX;
	tD1 = *tDP1;
	tD2 = *tDP2;
	tD3 = *(*(data1 + idY) + idX);
	tD4 = -*(*(data1 + tC2I.vtY + idY) + idX);
	*tDP1 = tD1 * tD3 - tD2 * tD4;
	*tDP2 = tD1 * tD4 + tD2 * tD3;
      }
    }
    **data0 *= **data1;
    **(data0 + tC2I.vtY) *= **(data1 + tC2I.vtY);
    *(*data0 + tC2I.vtX) *= *(*data1 + tC2I.vtX);
    *(*(data0 + tC2I.vtY) + tC2I.vtX) *= *(*(data1 + tC2I.vtY) + tC2I.vtX);
    AlgFourRealInv2D(data0, 1, size.vtX, size.vtY);
    REC_DBGW((REC_DBG_CROSS|REC_DBG_LVL_1),
	     WlzFromArray2D((void **)data0, size, org,
	     		    WLZ_GREY_UBYTE, WLZ_GREY_DOUBLE,
			    0.0, 255.0 / (1.0 + sqrt(*sSq0 * *sSq1)),
			    0, 0, NULL), 1);
  }
  REC_DBG((REC_DBG_CROSS|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecCrossCorrelate FX %d\n",
	  errFlag));
  return(errFlag);
}

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Preprocesses the given object and calculates it's (real)
*		Fourier transform.
* \param	data			Given object data and returned FFT.
* \param	sSq			Destination for sum of squares of the
*					object after any preprocessing.
* \param	obj			Given object.
* \param	org			Origin of data wrt the given woolz
* 					object.
* \param	size			The size of data.
* \param	ppCtrl			Preprocessing control data structure.
*/
static RecError	RecCCorObjToFour(double **data,
				 double *sSq, WlzObject *obj,
				 WlzIVertex2 org, WlzIVertex2 size,
				 RecPPControl *ppCtrl)
{
  RecError      errFlag = REC_ERR_NONE;
  WlzObject	*ppObj = NULL;

  REC_DBG((REC_DBG_CROSS|REC_DBG_PPROC | REC_DBG_LVL_2),
	  ("RecCCorObjToFour FE %p %p %p {%d, %d} {%d, %d} %p\n",
	   data, sSq, obj, org.vtX, org.vtY, size.vtX, size.vtY, ppCtrl));
  (void )memset(*data, 0, size.vtX * size.vtY * sizeof(double));
  ppObj = RecPreProcObj(obj, ppCtrl, &errFlag); /* Assigned by RecPreProcObj */
  REC_DBG((REC_DBG_CROSS|REC_DBG_PPROC | REC_DBG_LVL_2),
	  ("RecCCorObjToFour 01 0x%lx %d\n",
	   (unsigned long )ppObj, (ppObj)? ppObj->linkcount: 0));
  if((ppObj == NULL)  && (errFlag == REC_ERR_NONE))
  {
    errFlag = REC_ERR_WLZ;
  }
  else
  {
    REC_DBGW((REC_DBG_CROSS|REC_DBG_PPROC | REC_DBG_LVL_2), ppObj, 0);
  }
  if(errFlag == REC_ERR_NONE)
  {
    errFlag = RecErrorFromWlz(
    	      WlzToArray2D((void ***)&data, ppObj, size, org,
			   (ppCtrl->method & REC_PP_NOISE)? 1: 0,
	      	    	   WLZ_GREY_DOUBLE));
    if(errFlag == REC_ERR_NONE)
    {
      if(sSq)
      {
	(void )WlzArrayStats2D((void **)data, size, WLZ_GREY_DOUBLE,
			       NULL, NULL, NULL, sSq, NULL, NULL);
	*sSq *= size.vtX * size.vtY;
	REC_DBG((REC_DBG_CROSS|REC_DBG_PPROC | REC_DBG_LVL_2),
		("RecCCorObjToFour 02 %g\n",
		 *sSq));
      }
    }
    REC_DBGW((REC_DBG_CROSS|REC_DBG_LVL_2),
	     WlzFromArray2D((void **)data, size, org,
	     		    WLZ_GREY_UBYTE, WLZ_GREY_DOUBLE,
			    0.0, 1.0, 0, 0, NULL), 1);
    AlgFourReal2D(data, 1, size.vtX, size.vtY);
    REC_DBGW((REC_DBG_CROSS|REC_DBG_LVL_3),
	     WlzFromArray2D((void **)data, size, org,
	     		    WLZ_GREY_UBYTE, WLZ_GREY_DOUBLE,
			    0.0, 1.0, 0, 0, NULL), 1);
  }
  if(ppObj)
  {
    (void )WlzFreeObj(ppObj);
  }
  REC_DBG((REC_DBG_CROSS|REC_DBG_PPROC | REC_DBG_LVL_2),
  	  ("RecCCorObjToFour FX %d\n",
	   errFlag));
  return(errFlag);
}

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Performs a cross correlation on the two rectangular regions of
* 		interest of the given woolz objects.
* \param	newCC			Destination pointer for the normalised
* 				        [0.0 - 1.0] cross-correlation value
* 				        (may be NULL if notrequired).
* \param	data0			Data for obj0's FFT (size roiSz), if
* 					NULL the data space is freed and is
*					returned with the normalised (range
*					[0.0 - 1.0]) cross-correlation data
*					using *data0.
* \param	data1			Data for obj1's FFT as for data0,
* 					except that if not NULL the FFT  data
* 					is returned using *data1.
* \param	obj0			First object for cross-correlation.
* \param	obj1			Second object for cross-correlation.
* \param	roiCtr			The center of the region of interest
* 					wrt obj0 and obj1.
* \param	roiSz			The size of the region of interest.
* \param	ppCtrl			Preprocessing control data structure.
*/
RecError	RecCrossCorrelateROI(double *newCC,
				     double ***data0, double ***data1,
				     WlzObject *obj0, WlzObject *obj1,
				     WlzIVertex2 roiCtr, WlzIVertex2 roiSz,
				     RecPPControl *ppCtrl)
{
  RecError	errFlag = REC_ERR_NONE;
  WlzObject	*roiObj[2];
  double	**prvData[2];
  double	sSq[2];
  int		idX,
  		idY;
  double	normFac;
  WlzIVertex2	roiOrg;
  WlzIBox2	roiBox;
  WlzErrorNum	wlzErr = WLZ_ERR_NONE;

  REC_DBG((REC_DBG_CROSS|REC_DBG_LVL_FN|REC_DBG_LVL_1),
          ("RecCrossCorrelateROI FE 0x%lx 0x%lx 0x%lx 0x%lx"
	   " {%d %d} {%d %d} 0x%lx\n",
	   (unsigned long )newCC,
	   (unsigned long )data0, (unsigned long )data1,
	   (unsigned long )obj0, (unsigned long )obj1,
	   roiCtr.vtX, roiCtr.vtY, roiSz.vtX, roiSz.vtY,
	   ppCtrl));
  roiObj[0] = NULL;
  roiObj[1] = NULL;
  prvData[0] = NULL;
  prvData[1] = NULL;
  if(((obj0 == NULL) && ((data0 == NULL) || (*data0 == NULL))) ||
     ((obj1 == NULL) && ((data1 == NULL) || (*data1 == NULL))) ||
     (roiSz.vtX < REC_FOUR_DIM_MIN) || (roiSz.vtX > REC_FOUR_DIM_MAX) ||
     (roiSz.vtY < REC_FOUR_DIM_MIN) || (roiSz.vtY > REC_FOUR_DIM_MAX))
  {
    errFlag = REC_ERR_FUNC;
  }
  else if((obj0 == NULL) || (obj0->type != WLZ_2D_DOMAINOBJ) ||
          (obj0->domain.core == NULL) || (obj0->values.core == NULL) ||
          (obj1 == NULL) || (obj1->type != WLZ_2D_DOMAINOBJ) ||
          (obj1->domain.core == NULL) || (obj1->values.core == NULL))
  {
    errFlag = REC_ERR_WLZ;
  }
  if(errFlag == REC_ERR_NONE)
  {
    if(data0 && *data0)
    {
      prvData[0] = *data0;
    }
    else if(AlcDouble2Malloc(&(prvData[0]),
    			  roiSz.vtY, roiSz.vtX) != ALC_ER_NONE)
    {
      errFlag = REC_ERR_MALLOC;
    }
  }
  if(errFlag == REC_ERR_NONE)
  {
    if(data1 && *data1)
    {
      prvData[1] = *data1;
    }
    else if(AlcDouble2Malloc(&(prvData[1]),
    			  roiSz.vtY, roiSz.vtX) != ALC_ER_NONE)
    {
      errFlag = REC_ERR_MALLOC;
    }
  }
  if(errFlag == REC_ERR_NONE)
  {
    roiOrg.vtX = roiCtr.vtX - (roiSz.vtX / 2);
    roiOrg.vtY = roiCtr.vtY - (roiSz.vtY / 2);
    roiBox.xMin = roiOrg.vtX;
    roiBox.yMin = roiOrg.vtY;
    roiBox.xMax = roiBox.xMin + roiSz.vtX - 1;
    roiBox.yMax = roiBox.yMin + roiSz.vtY - 1;
    if(((roiObj[0] = WlzAssignObject(
    		     WlzClipObjToBox2D(obj0, roiBox,
		     		       &wlzErr), NULL)) == NULL) ||
       (wlzErr != WLZ_ERR_NONE) ||
       ((roiObj[1] = WlzAssignObject(
    		     WlzClipObjToBox2D(obj1, roiBox,
		     		       &wlzErr), NULL)) == NULL) ||
       (wlzErr != WLZ_ERR_NONE))
    {
      errFlag = RecErrorFromWlz(wlzErr);
    }
    else
    {
      REC_DBGW((REC_DBG_CROSS|REC_DBG_LVL_1),
	       roiObj[0], 0);
      REC_DBGW((REC_DBG_CROSS|REC_DBG_LVL_1),
	       roiObj[1], 0);
    }
  }
  if(errFlag == REC_ERR_NONE)
  {
    errFlag = RecCrossCorrelate(prvData[0], prvData[1],  &(sSq[0]), &(sSq[1]),
    				REC_CCFLAG_NONE, roiObj[0], roiObj[1],
				roiOrg, roiSz,
				ppCtrl);
  }
  if(errFlag == REC_ERR_NONE)
  {
    normFac = sqrt(sSq[0] * sSq[1]);
    normFac = (normFac > DBL_EPSILON)? (1 / normFac): 0.0;
    if(newCC)
    {
      *newCC = **(prvData[0]) * normFac;
      if(data0 && (*data0 == NULL))
      {
	*data0 = prvData[0];
	for(idY = 0; idY < roiSz.vtY; ++idY)
	{
	  for(idX = 0; idX < roiSz.vtX; ++idX)
	  {
	    *(*(*data0 + idY) + idX) *= normFac;
	  }
	}
      }
      if(data1 && (*data1 == NULL))
      {
	*data1 = prvData[1];
      }
    }
    REC_DBGW((REC_DBG_CROSS|REC_DBG_LVL_1),
	     WlzFromArray2D((void **)(prvData[0]), roiSz, roiOrg,
	     		     WLZ_GREY_UBYTE, WLZ_GREY_DOUBLE,
			     0.0, 1.0, 0, 0, NULL), 1);
  }
  if(prvData[0] && (data0 == NULL))
  {
    AlcDouble2Free(prvData[0]);
  }
  if(prvData[1] && (data1 == NULL))
  {
    AlcDouble2Free(prvData[1]);
  }
  if(roiObj[0])
  {
    (void )WlzFreeObj(roiObj[0]);
  }
  if(roiObj[1])
  {
    (void )WlzFreeObj(roiObj[1]);
  }
  REC_DBG((REC_DBG_CROSS|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecCrossCorrelateROI FX %d\n",
	   errFlag));
  return(errFlag);
}
