#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _ReconstructTranMatch_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libReconstruct/ReconstructTranMatch.c
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
* \brief	Provides functions for the automatic registration of
*		a single pair of serial sections for the Reconstruct
*		library.
* \ingroup	Reconstruct
*/
#include <Reconstruct.h>
#include <string.h>
#include <float.h>

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Uses the cross correlation of the two given objects to
*               find the shift (translation) which gives the best match
*               (cross correlation maximum) between them.
*               Data is accessed as:
* \verbatim
                  *(*(data + line) + column), x == column, y == line.
  \endverbatim
* \param	shift			Destination pointer for the
*					cross-correlation peak position.
* \param	value			Destination pointer for the
*					cross-correlation peak value.
* \param	obj0			First of two WLZ_2D_DOMAINOBJ objects.
* \param	obj1			Second of two WLZ_2D_DOMAINOBJ objects.
* \param	maxShift		Maximum shift.
* \param	ppCtrl			Pre-processing control.
*/
RecError	RecTranMatch(WlzDVertex2 *shift, double *value,
			      WlzObject *obj0, WlzObject *obj1,
			      WlzIVertex2 maxShift,
			      RecPPControl *ppCtrl)
{
  WlzIBox2	objLim;
  RecError	errFlag = REC_ERR_NONE;
  WlzIVertex2	origin,
		size;
  double	sSq0,
  		sSq1;
  double	**data0 = NULL,
  		**data1 = NULL;
  const WlzDVertex2 zeroCoord2D = {0.0, 0.0};

  REC_DBG((REC_DBG_LVL_1|REC_DBG_LVL_FN|REC_DBG_TRAN),
	  ("RecTranMatch FE 0x%lx 0x%lx 0x%lx 0x%lx {%d %d} 0x%lx\n",
	   (unsigned long )shift, (unsigned long )value,
	   (unsigned long )obj0, (unsigned long)obj1,
	   maxShift.vtX, maxShift.vtY, (unsigned long )ppCtrl));
  REC_DBGW((REC_DBG_TRAN|REC_DBG_LVL_2), obj0, 0);
  REC_DBGW((REC_DBG_TRAN|REC_DBG_LVL_2), obj1, 0);
  *value = 0.0;
  *shift = zeroCoord2D;
  objLim.xMin = WLZ_MIN(obj0->domain.i->kol1, obj1->domain.i->kol1);
  objLim.yMin = WLZ_MIN(obj0->domain.i->line1, obj1->domain.i->line1);
  objLim.xMax = WLZ_MAX(obj0->domain.i->lastkl, obj1->domain.i->lastkl);
  objLim.yMax = WLZ_MAX(obj0->domain.i->lastln, obj1->domain.i->lastln);
  origin.vtX = objLim.xMin - maxShift.vtX;
  origin.vtY = objLim.yMin - maxShift.vtY;
  size.vtX = objLim.xMax - objLim.xMin + (2 * maxShift.vtX) + 1;
  size.vtY = objLim.yMax - objLim.yMin + (2 * maxShift.vtY) + 1;
  REC_DBG((REC_DBG_TRAN|REC_DBG_LVL_2),
	  ("RecTranMatch 01 {%d %d %d %d} {%d %d} {%d %d}\n",
	   objLim.xMin, objLim.yMin, objLim.xMax, objLim.yMax,
	   origin.vtX, origin.vtY, size.vtX, size.vtY));
  (void )RecPowerOfTwoC2I(&size, size);
  REC_DBG((REC_DBG_TRAN|REC_DBG_LVL_2),
	  ("RecTranMatch 02 {%d %d}\n",
	   size.vtX, size.vtY));
  if((AlcDouble2Malloc(&data0, size.vtY, size.vtX) != ALC_ER_NONE) ||
     (data0 == NULL))
  {
    errFlag = REC_ERR_MALLOC;
  }
  if(errFlag == REC_ERR_NONE)
  {
    if((AlcDouble2Malloc(&data1, size.vtY, size.vtX) != ALC_ER_NONE) ||
       (data1 == NULL))
    {
      errFlag = REC_ERR_MALLOC;
    }
  }
  if(errFlag == REC_ERR_NONE)
  {
    errFlag = RecCrossCorrelate(data0, data1, &sSq0, &sSq1,
    			        REC_CCFLAG_NONE, obj0, obj1,
    			        origin, size, ppCtrl);
  }
  if(errFlag == REC_ERR_NONE)
  {
    REC_DBG((REC_DBG_TRAN|REC_DBG_LVL_2),
	    ("RecTranMatch 03 %f %f\n",
	     sSq0, sSq1));
    if((sSq0 <= 0.0) || (sSq1 <= 0.0))
    {
      errFlag = REC_ERR_WLZ;
    }
  }
  if(errFlag == REC_ERR_NONE)
  {
    RecTranFindPeak(shift, value, data0, size, maxShift);
    *value /= sqrt(sSq0 * sSq1);
    REC_DBG((REC_DBG_TRAN|REC_DBG_LVL_2),
	    ("RecTranMatch 04 %f {%f %f}\n",
	     *value, shift->vtX, shift->vtY));
  }
  if(data0)
  {
    (void )AlcDouble2Free(data0);
  }
  if(data1)
  {
    (void )AlcDouble2Free(data1);
  }
  REC_DBG((REC_DBG_LVL_1|REC_DBG_LVL_FN|REC_DBG_TRAN),
	  ("RecTranMatch FX %d\n",
	   errFlag));
  return(errFlag);
}

/*!
* \ingroup	Reconstruct
* \brief	Finds the maximum correlation value in the given two
*               dimensional array. Because the correlation data are
*               stored in wrap-around order and the maximum is known
*               to lie within a limited search range, only the four
*               corners of the array are searched.
* \param	shift			Destination pointer for the position
*					of the maximum value.
* \param	pValue			Destination pointer for the maximum
* 					value.
* \param	data			Data to search for maximum.
* \param	size			Size of data array.
* \param	search			Size of search range.
*/
void		RecTranFindPeak(WlzDVertex2 *shift, double *pValue,
			        double **data, WlzIVertex2 size,
			        WlzIVertex2 search)
{
  int		yIdx,
		xIdx,
		xMax,
		yMax;
  double	maxVal;
  double	*left,
		*right;

  REC_DBG((REC_DBG_LVL_1|REC_DBG_LVL_FN|REC_DBG_TRAN),
	  ("RecTranFindPeak FE 0x%lx 0x%lx 0x%lx {%d %d} {%d %d}\n",
	   (unsigned long )shift, (unsigned long )pValue, (unsigned long )data,
	   size.vtX, size.vtY, search.vtX, search.vtY));
  xMax = 0;
  yMax = 0;
  maxVal = **data;
  for(yIdx = 0; yIdx < search.vtY; ++yIdx) 	/* search two bottom corners */
  {
    left = *(data + yIdx);
    right = left + size.vtX - 1;
    for(xIdx = 0; xIdx < search.vtX; ++xIdx)
    {
      if(*left > maxVal)
      {
	maxVal = *left;
	xMax = xIdx;
	yMax = yIdx;
      }
      if(*right > maxVal)
      {
	maxVal = *right;
	xMax = -(xIdx + 1);
	yMax = yIdx;
      }
      ++left;
      --right;
    }
  }
  for(yIdx = 1; yIdx <= search.vtY; ++yIdx)	   /* search two top corners */
  {
    left = *(data + size.vtY - yIdx);
    right = left + size.vtX - 1;
    for(xIdx= 0; xIdx < search.vtX; ++xIdx)
    {
      if(*left > maxVal)
      {
	maxVal = *left;
	xMax = xIdx;
	yMax = -yIdx;
      }
      if(*right > maxVal)
      {
	maxVal = *right;
	xMax = -(xIdx + 1);
	yMax = -yIdx;
      }
      ++left;
      --right;
    }
  }
  REC_DBG((REC_DBG_LVL_1|REC_DBG_TRAN),
	  ("RecTranFindPeak 01 %d %d %f\n",
	   xMax, yMax, maxVal));
  if(pValue)
  {
    *pValue = maxVal;
  }
  if(shift)
  {
    shift->vtX = xMax;
    shift->vtY = yMax;
  }
  REC_DBG((REC_DBG_LVL_1|REC_DBG_TRAN),
	  ("RecTranFindPeak 02 {%f %f} %f\n",
	   shift? shift->vtX: HUGE_VAL, shift? shift->vtY: HUGE_VAL,
	   pValue? *pValue: HUGE_VAL));
  REC_DBG((REC_DBG_LVL_1|REC_DBG_LVL_FN|REC_DBG_TRAN),
	  ("RecTranFindPeak FX\n"));
}
