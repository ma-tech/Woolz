#pragma ident "MRC HGU $Id$"
/************************************************************************
* Project:	Mouse Atlas
* Title:        ReconstructRotMatch.c				
* Date:         April 1999
* Author:       Bill Hill                                              
* Copyright:    1999 Medical Research Council, UK.
*		All rights reserved.				
* Address:	MRC Human Genetics Unit,			
*		Western General Hospital,			
*		Edinburgh, EH4 2XU, UK.				
* Purpose:      Provides functions for the automatic registration of
*		a single pair of serial sections for the MRC Human
*		Genetics Unit reconstruction library.		
* $Revision$
* Maintenance:  Log changes below, with most recent at top of list.    
************************************************************************/
#include <Reconstruct.h>
#include <string.h>
#include <float.h>

static void	RecFindRotPeak(double *angle, double *value, double **data,
			       double angleInc, int size);

/************************************************************************
* Function:	RecRotMatch					
* Returns:	int:			Non zero if match fails.
* Purpose:	Performs polar resampling of the givn objects and then
*		uses cross correlation to find the angle of rotation
*		which gives the best match between the given objects.
*		Data are accessed as				
*		  *(*(data + line) + column), x == column, y == line.
* Global refs:	-						
* Parameters:	double *angle:		Destination for angle (radians)
*					of rotation.		
*		double *value:		Destination for cross-	
*					correlation peak value.	
*		WlzObject *obj0:	First of two type 1 objects.
*		WlzObject *obj1:	Second of two type 1 objects.
*		WlzIVertex2 cRot:	Center of rotation for objects.
*		double angleInc:	Angle increment (radians).
*		double distInc:		Distance increment.	
*		int maxRadiusFlag:	Use maximum radius for the
*					polar resampling if non zero.
*		RecPPControl *ppCtrl:	Pre-processing control.	
************************************************************************/
RecError	RecRotMatch(double *angle, double *value,
			    WlzObject *obj0, WlzObject *obj1,
			    WlzIVertex2 cRot,
			    double angleInc, double distInc, int maxRadiusFlag,
			    RecPPControl *ppCtrl)
{
  int		length,
		p2Len;
  WlzErrorNum	wlzErr = WLZ_ERR_NONE;
  RecError	errFlag = REC_ERR_NONE;
  WlzIVertex2	size,
		origin;
  WlzObject	*pObj0 = NULL,
		*pObj1 = NULL;
  double	sSq0,
		sSq1;
  double	**data0 = NULL,
  		**data1 = NULL;

  REC_DBG((REC_DBG_ROT|REC_DBG_LVL_1|REC_DBG_LVL_FN),
	  ("RecRotMatch FE 0x%x 0x%x 0x%x 0x%x {%d %d} %f %f %d 0x%x\n",
	   angle, value, obj0, obj1, cRot.vtX, cRot.vtY,
	   angleInc, distInc, maxRadiusFlag, ppCtrl));
  *value = 0.0;
  *angle = 0.0;
  length = WLZ_NINT(2.0 * WLZ_M_PI / angleInc);
  p2Len = RecPowerOfTwo(&length, length);
  REC_DBG((REC_DBG_ROT|REC_DBG_LVL_2),
	  ("RecRotMatch 01 %d %d\n",
	   length, p2Len));
  REC_DBGW((REC_DBG_ROT | REC_DBG_LVL_3), obj0, 0);
  REC_DBGW((REC_DBG_ROT | REC_DBG_LVL_3), obj1, 0);
  if(((pObj0 = WlzAssignObject(
               WlzPolarSample(obj0, cRot, angleInc, distInc,
			     length, maxRadiusFlag,
			     &wlzErr), NULL)) == NULL) ||
     (wlzErr != WLZ_ERR_NONE) ||
     ((pObj1 = WlzAssignObject(
     	       WlzPolarSample(obj1, cRot, angleInc, distInc,
			      length, maxRadiusFlag,
			      &wlzErr), NULL)) == NULL) ||
     (wlzErr != WLZ_ERR_NONE))
  {
    errFlag = REC_ERR_WLZ;
  }
  if(errFlag == REC_ERR_NONE)
  {

    REC_DBGW((REC_DBG_ROT | REC_DBG_LVL_2), pObj0, 0);
    REC_DBGW((REC_DBG_ROT | REC_DBG_LVL_2), pObj1, 0);
    origin.vtX = WLZ_MIN(pObj0->domain.i->kol1, pObj1->domain.i->kol1);
    origin.vtY = WLZ_MIN(pObj0->domain.i->line1, pObj1->domain.i->line1);
    size.vtX = WLZ_MAX(pObj0->domain.i->lastkl - pObj0->domain.i->kol1 + 1,
		       pObj1->domain.i->lastkl - pObj1->domain.i->kol1 + 1);
    size.vtY = WLZ_MAX(pObj0->domain.i->lastln - pObj0->domain.i->line1 + 1,
		       pObj1->domain.i->lastln - pObj1->domain.i->line1 + 1);
    
    (void )RecPowerOfTwoC2I(&size, size);
    REC_DBG((REC_DBG_ROT|REC_DBG_LVL_2),
	    ("RecRotMatch 02 {%d %d} {%d %d}\n",
	     origin.vtX, origin.vtY, size.vtX, size.vtY));
  }
  if(errFlag == REC_ERR_NONE)
  {
    if((AlcDouble2Malloc(&data0, size.vtY, size.vtX) != ALC_ER_NONE) ||
       (data0 == NULL))
    {
      errFlag = REC_ERR_MALLOC;
    }
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
    			        REC_CCFLAG_NONE,
    			        pObj0, pObj1, origin, size, ppCtrl);
    REC_DBG((REC_DBG_ROT|REC_DBG_LVL_2),
	    ("RecRotMatch 03 %d %g %g\n",
	     (int )errFlag, sSq0, sSq1));
    if((sSq0 <= 0.0) || (sSq1 <= 0.0))
    {
      errFlag = REC_ERR_WLZ;
    }
  }
  (void )WlzFreeObj(pObj0);
  (void )WlzFreeObj(pObj1);
  if(errFlag == REC_ERR_NONE)
  {
    RecFindRotPeak(angle, value, data0, angleInc, size.vtY);
    *value /= sqrt(sSq0 * sSq1);
    REC_DBG((REC_DBG_ROT|REC_DBG_LVL_1),
	    ("RecRotMatch 04 %g %g\n",
	     *value, *angle));
  }
  if(data0)
  {
    (void )AlcDouble2Free(data0);
  }
  if(data1)
  {
    (void )AlcDouble2Free(data1);
  }
  REC_DBG((REC_DBG_ROT|REC_DBG_LVL_1|REC_DBG_LVL_FN),
	  ("RecRotMatch FX %d\n",
	   errFlag));
  return(errFlag);
}

/************************************************************************
* Function:	RecFindRotPeak					
* Returns:	void						
* Purpose:	Finds peak rotation value in cross correlation data.
*		The search is particularly simple because only the
*		first column of the cross correlation data needs to be
*		searched for the angle of rotation.		
*		Data are in wrap around order.			
*		A least squares quadratic is fitted to the cross
*		correlation maximum.				
* Global refs:	-						
* Parameters:	double *angle:		Angle (degrees) for the cross-
*					correlation maximum.	
*		double *value:		Cross correlation maximum.
*		double **data:		Cross correlation data.	
*		double angleInc:	Angle increment used to find 
*					angle from data index.	
*		int size:		Size of data array.	
************************************************************************/
static void	RecFindRotPeak(double *angle, double *value, double **data,
			       double angleInc, int size)
{
  int		idx,
		maxIdx;
  double	maxVal;

  REC_DBG((REC_DBG_ROT|REC_DBG_LVL_1|REC_DBG_LVL_FN),
	  ("RecFindRotPeak FE 0x%lx 0x%lx 0x%lx %f %d\n",
	   (unsigned long )angle, (unsigned long )value, (unsigned long )data,
	   angleInc, size));
  maxIdx = 0;
  maxVal = (**data);
  for(idx = 1; idx < size; ++idx)
  {
    if(**(data + idx) > maxVal)
    {
      maxVal = **(data + idx);
      maxIdx = idx;
    }
  }
  if(maxIdx > (size / 2))
  {
    maxIdx -= size;
  }
  REC_DBG((REC_DBG_ROT|REC_DBG_LVL_1),
	  ("RecFindRotPeak 01 %d %f\n",
	   maxIdx, maxVal));
  if(angle)
  {
    *angle = maxIdx * angleInc;
  }
  if(value)
  {
    *value = maxVal;
  }
  REC_DBG((REC_DBG_ROT|REC_DBG_LVL_1),
	  ("RecFindRotPeak 02 %f %f\n",
	   angle? *angle: HUGE_VAL, value? *value: HUGE_VAL));
  REC_DBG((REC_DBG_ROT|REC_DBG_LVL_1|REC_DBG_LVL_FN),
	  ("RecFindRotPeak FX\n"));
}
