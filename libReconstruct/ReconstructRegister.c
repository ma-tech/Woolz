#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _ReconstructRegister_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libReconstruct/ReconstructRegister.c
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
* \brief	Functions for the automatic registration of a single
*		pair of serial sections for the Reconstruct library.
* \ingroup	Reconstruct
*/

#include <Reconstruct.h>
#include <string.h>
#include <float.h>

static RecError			RecRegPrincipal(
				  WlzAffineTransform **transf,
				  WlzObject **transfObj,
				  WlzDVertex2 cMass0,
				  WlzDVertex2 cMass1,
				  double angle,
				  WlzObject *givenObj);
static RecError			RecRegRotate(
				  WlzAffineTransform **transf,
				  WlzObject **transfObj,
				  WlzIVertex2 cMass,
				  double angle,
				  WlzObject *secObj);
static RecError			RecRegTranslate(
				  WlzAffineTransform **transf,
				  WlzObject **transfObj,
				  WlzDVertex2 trans,
				  WlzObject *secObj);

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Calculates the affine transform which when applied
*               to the second object brings it into register with the
*               first object.
* \param	dstTrans		Destination pointer for transform.
* \param	dstCrossC		Destination pointer for the
*					cross-correlation value.
* \param	dstIter			Destination pointer for the number
*                                       of iterations.
* \param	rCtrl			The registration control data
*                                       structure.
* \param	ppCtrl			Pre-processing control data
*                                       structure.
* \param	obj0			First object.
* \param	obj1			Second object.
* \param	workFn			Application supplied work
*                                       function.
* \param	workData		Application supplied data for
*                                       the work function.
* \param	eMsg			Destination pointer for messages.
*/
RecError	RecRegisterPair(WlzAffineTransform **dstTrans,
				double *dstCrossC, int *dstIter,
				RecControl *rCtrl, RecPPControl *ppCtrl,
			        WlzObject *obj0, WlzObject *obj1,
				RecWorkFunction workFn, void *workData,
				char **eMsg)
{
  int		approach = 0,
		approach0 = 0,
		approach1 = 1,
		freeObjFlag = 0;
  double	correl = 1.0,
		distInc = 1.0,
		angleInc = (2.0 * WLZ_M_PI / 512.0),
		tD0 = 0.0,
		tD1 = 0.0;
  WlzAffineTransform	*tTr;
  WlzIVertex3	samFac;
  WlzIVertex2	tIV0,
		maxShift;
  WlzDVertex2	tDV0,
  		cMass0,
		cMass1;
  WlzObject	*ppObj0 = NULL,
  		*ppObj1 = NULL,
		*trObj = NULL;
  WlzErrorNum	wlzErr = WLZ_ERR_NONE;
  RecError	errFlag = REC_ERR_NONE;
  double	theta[2];
  RecState	states[2];
  WlzAffineTransformPrim prim;
  RecPPControl	newPP;

  REC_DBG((REC_DBG_REG|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecRegisterPair FE 0x%lx 0x%lx 0x%lx 0x%lx 0x%lx\n",
	   (unsigned long )dstTrans,
	   (unsigned long )dstCrossC, (unsigned long )dstIter,
	   (unsigned long )rCtrl, (unsigned long )ppCtrl,
	   (unsigned long )obj0, (unsigned long )obj1,
	   (unsigned long )eMsg));
  (void )memset(&(states[0]), 0, sizeof(RecState));
  (void )memset(&(states[1]), 0, sizeof(RecState));
  if((rCtrl == NULL) || (ppCtrl == NULL) || (obj0 == NULL) || (obj1 == NULL))
  {
    errFlag = REC_ERR_FUNC;
  }
  else
  {
    if((obj0->type != WLZ_2D_DOMAINOBJ) || (obj1->type != WLZ_2D_DOMAINOBJ) ||
       (obj0->domain.core == NULL) || (obj1->domain.core == NULL) ||
       (obj0->values.core == NULL) || (obj1->values.core == NULL))
    {
      errFlag = REC_ERR_WLZ;
    }
  }
  if(errFlag == REC_ERR_NONE)
  {
    maxShift.vtX = rCtrl->xLim;
    maxShift.vtY = rCtrl->yLim;
    newPP = *ppCtrl;
    if(ppCtrl->method & REC_PP_SAMPLE)
    {
      maxShift.vtX = (maxShift.vtX + ppCtrl->sample.factor - 1) /
      		     ppCtrl->sample.factor;
      maxShift.vtY = (maxShift.vtY + ppCtrl->sample.factor - 1) /
      		     ppCtrl->sample.factor;
      samFac.vtX = ppCtrl->sample.factor;
      samFac.vtY = ppCtrl->sample.factor;
      samFac.vtZ = 1;
      newPP.sample.factor = 1;
      newPP.sample.function = WLZ_SAMPLEFN_NONE;
      newPP.method = (RecPPMethod )((unsigned int )newPP.method & 
      				    ~(unsigned int )(REC_PP_SAMPLE));
      obj0 = WlzAssignObject(WlzSampleObj(obj0, samFac, ppCtrl->sample.function,
      				          &wlzErr), NULL);
      if(obj0 && (wlzErr == WLZ_ERR_NONE))
      {
        obj1 = WlzAssignObject(
	       WlzSampleObj(obj1, samFac, ppCtrl->sample.function,
	       		    &wlzErr), NULL);
      }
      if(obj0 && obj1 && (wlzErr == WLZ_ERR_NONE))
      {
        freeObjFlag = 1;
      }
      errFlag = RecErrorFromWlz(wlzErr);
    }
    else
    {
      obj0 = WlzAssignObject(obj0, NULL);
      obj1 = WlzAssignObject(obj1, NULL);
      freeObjFlag = 1;
    }
  }
  if((errFlag == REC_ERR_NONE) &&
     ((rCtrl->method & REC_MTHD_PRINC) || (rCtrl->method & REC_MTHD_ROTATE) ||
      (rCtrl->method & REC_MTHD_TRANS)))
  {
    ppObj0 = RecPreProcObj(obj0, &newPP,
    			   &errFlag); 	        /* Assigned by RecPreProcObj */
    if(ppObj0 && (errFlag == REC_ERR_NONE))
    {
      REC_DBGW((REC_DBG_REG|REC_DBG_LVL_2), obj0, 0);
      cMass0 = WlzCentreOfMass2D(ppObj0, 0, NULL, &wlzErr);
      errFlag = RecErrorFromWlz(wlzErr);
      if(errFlag == REC_ERR_NONE)
      {
	if(rCtrl->method & REC_MTHD_PRINC)
	{
	  ppObj1 = RecPreProcObj(obj1, &newPP,
	  			 &errFlag); 	/* Assigned by RecPreProcObj */
	  if(ppObj1 && (errFlag == REC_ERR_NONE))
	  {
	    REC_DBGW((REC_DBG_REG|REC_DBG_LVL_2), obj1, 0);
	    cMass1 = WlzCentreOfMass2D(ppObj1, 0, NULL, &wlzErr);
	    if(wlzErr == WLZ_ERR_NONE)
	    {
	      tD0 = WlzPrincipalAngle(ppObj0, cMass0, 0, &wlzErr);
	    }
	    if(wlzErr == WLZ_ERR_NONE)
	    {
	      tD1 = WlzPrincipalAngle(ppObj1, cMass1, 0, &wlzErr);
	    }
	    errFlag = RecErrorFromWlz(wlzErr);
	    if(errFlag == REC_ERR_NONE)
	    {
	      theta[0] = tD0 - tD1;
	      theta[1] = theta[0] + WLZ_M_PI;
	      REC_DBG((REC_DBG_REG|REC_DBG_LVL_1),
		      ("RecRegisterPair 01 {%g %g} {%g %g} %f %f\n",
		       cMass0.vtX, cMass0.vtY, cMass1.vtX, cMass1.vtY,
		       tD0, tD1));
	      if(WLZ_ABS(rCtrl->rLim) < 90.0) 	 /* Use limit to reduce work */
	      {
		if(WLZ_ABS(theta[0]) < WLZ_ABS(theta[1]))
		{
		  approach0 = 0;
		}
		else
		{
		  approach0 = 1;
		}
		approach1 = approach0;
	      }
	      REC_DBG((REC_DBG_REG|REC_DBG_LVL_1),
		      ("RecRegisterPair 02 %d %d\n",
		       approach0, approach1));
	      approach = approach0;
	    }
	  }
	}
	else
	{
	  approach = 0;
	  approach0 = 0;
	  approach1 = 0;
	}
      }
    }
  }
  while((errFlag == REC_ERR_NONE) && (approach <= approach1))
  {
    states[approach].approach = approach;
    if(rCtrl->method & REC_MTHD_PRINC)
    {
      errFlag = RecRegPrincipal(&(states[approach].transform), &trObj,
		                cMass0, cMass1, theta[approach], obj1);
      if(errFlag == REC_ERR_NONE)
      {
	if(workFn)
	{
	  states[approach].lastMethod = REC_MTHD_PRINC;
	  states[approach].errFlag = errFlag;
	  (*workFn)(&(states[approach]), workData);
	  errFlag = states[approach].errFlag;
	}
      }
    }
    else if((rCtrl->method & REC_MTHD_IDENTITY) ||
	    (dstTrans == NULL) || (*dstTrans == NULL))
    {
      states[approach].transform = WlzAssignAffineTransform(
      		WlzAffineTransformFromPrimVal(WLZ_TRANSFORM_2D_AFFINE,
      				  	0.0, 0.0, 0.0,
					1.0, 0.0, 0.0,
					0.0, 0.0, 0.0,
					0, &wlzErr), NULL);
      errFlag = RecErrorFromWlz(wlzErr);
      if(errFlag == REC_ERR_NONE)
      {
        trObj = WlzAssignObject(obj1, NULL);
      }
    }
    else
    {
      if(ppCtrl->method & REC_PP_SAMPLE)
      {
	tTr = *dstTrans;
	if((wlzErr = WlzAffineTransformPrimGet(tTr, &prim)) == WLZ_ERR_NONE)
	{
	  states[approach].transform = WlzAssignAffineTransform(
		  WlzAffineTransformFromPrimVal(WLZ_TRANSFORM_2D_AFFINE,
					  prim.tx / samFac.vtX,
					  prim.ty / samFac.vtY,
					  0.0,
					  1.0, prim.theta, 0.0,
					  0.0, 0.0, 0.0,
					  0, &wlzErr), NULL);
	}
        errFlag = RecErrorFromWlz(wlzErr);
      }
      else
      {
        states[approach].transform = WlzAssignAffineTransform(*dstTrans, NULL);
      }
      trObj = WlzAssignObject(WlzAffineTransformObj(obj1,
      						    states[approach].transform,
						    WLZ_INTERPOLATION_NEAREST,
      					            &wlzErr), NULL);
      errFlag = RecErrorFromWlz(wlzErr);
    }
    if(errFlag == REC_ERR_NONE)
    {
      if(rCtrl->method & REC_MTHD_TRANS)
      {
	errFlag = RecTranMatch(&tDV0, &correl, obj0, trObj, maxShift,
			       &newPP);
	REC_DBG((REC_DBG_REG|REC_DBG_LVL_1),
		("RecRegisterPair 03 %d %f {%f %f} %d\n",
		 approach, correl, tDV0.vtX, tDV0.vtY, errFlag));
	if((errFlag == REC_ERR_NONE) && (correl > 0.3) &&
	   ((WLZ_ABS(tDV0.vtX) >= 0.5) || (WLZ_ABS(tDV0.vtY) >= 0.5)))
	{
	  errFlag = RecRegTranslate(&(states[approach].transform),
				   &trObj, tDV0, obj1);
	  REC_DBGW((REC_DBG_REG|REC_DBG_LVL_2), obj0, 0);
	  REC_DBGW((REC_DBG_REG|REC_DBG_LVL_1), trObj, 0);
	}
	if(workFn)
	{
	  states[approach].lastMethod = REC_MTHD_TRANS;
	  states[approach].errFlag = errFlag;
	  (*workFn)(&(states[approach]), workData);
	  errFlag = states[approach].errFlag;
	}
      }
      else
      {
	tDV0.vtX = 0.4;			  /* Just to get the iteration going */
	tDV0.vtY = 0.4;
      }
    }
    while((errFlag == REC_ERR_NONE) && 
	  (rCtrl->method & (REC_MTHD_TRANS|REC_MTHD_ROTATE)) &&
	  (states[approach].iteration < rCtrl->itLim) &&
	  ((WLZ_ABS(tDV0.vtX) >= 0.3) || (WLZ_ABS(tDV0.vtY) >= 0.3)))
    {
      if(rCtrl->method & REC_MTHD_ROTATE)
      {
	tIV0.vtX = WLZ_NINT(cMass0.vtX);
	tIV0.vtY = WLZ_NINT(cMass0.vtY);
	errFlag = RecRotMatch(&tD0, &correl, obj0, trObj,
			      tIV0, angleInc, distInc, 0, &newPP);
	REC_DBG((REC_DBG_REG|REC_DBG_LVL_1),
		("RecRegisterPair 04 %d %d %f %f %d\n",
		 approach, states[approach].iteration, correl,
		 tD0, errFlag));
	if((errFlag == REC_ERR_NONE) && (correl > 0.3) &&
	   (WLZ_ABS(tD0) > 0.001))
	{
	  tIV0.vtX = WLZ_NINT(cMass0.vtX);
	  tIV0.vtY = WLZ_NINT(cMass0.vtY);
	  errFlag = RecRegRotate(&(states[approach].transform),
		                 &trObj, tIV0, tD0, obj1);
	  REC_DBGW((REC_DBG_REG|REC_DBG_LVL_2), obj0, 0);
	  REC_DBGW((REC_DBG_REG|REC_DBG_LVL_1), trObj, 0);
	}
	if((errFlag == REC_ERR_NONE) && workFn)
	{
	  states[approach].lastMethod = REC_MTHD_ROTATE;
	  states[approach].errFlag = errFlag;
	  (*workFn)(&(states[approach]), workData);
	  errFlag = states[approach].errFlag;
	}
      }
      if(rCtrl->method & REC_MTHD_TRANS)
      {
	if(errFlag == REC_ERR_NONE)
	{
	  errFlag = RecTranMatch(&tDV0, &correl, obj0, trObj,
				  maxShift, &newPP);
	  REC_DBG((REC_DBG_REG|REC_DBG_LVL_1),
		  ("RecRegisterPair 03 %d %d %f {%f %f} %d\n",
		   approach, states[approach].iteration,
		   correl, tDV0.vtX, tDV0.vtY, errFlag));
	}
	if((errFlag == REC_ERR_NONE) && (correl > 0.3) &&
	   ((WLZ_ABS(tDV0.vtX) >= 0.3) || (WLZ_ABS(tDV0.vtY) >= 0.3)))
	{
	  errFlag = RecRegTranslate(&(states[approach].transform),
			            &trObj, tDV0, obj1);
	  REC_DBGW((REC_DBG_REG|REC_DBG_LVL_2), obj0, 0);
	  REC_DBGW((REC_DBG_REG|REC_DBG_LVL_1), trObj, 0);
	}
	if((errFlag == REC_ERR_NONE) && workFn)
	{
	  states[approach].lastMethod = REC_MTHD_TRANS;
	  states[approach].errFlag = errFlag;
	  (*workFn)(&(states[approach]), workData);
	  errFlag = states[approach].errFlag;
	}
      }
      if(errFlag == REC_ERR_NONE)
      {
	++(states[approach].iteration);
      }
    }
    (void )WlzFreeObj(trObj);
    trObj = NULL;
    states[approach].correl = correl;
    ++approach;
  }
  if(errFlag == REC_ERR_NONE)
  {
    if(approach0 == approach1)
    {
      approach = approach0;
    }
    else
    {
      if(states[0].correl >= states[1].correl)
      {
	approach = 0;
	(void )WlzFreeAffineTransform(states[1].transform);
	states[1].transform = NULL;
      }
      else
      {
	approach = 1;
	(void )WlzFreeAffineTransform(states[0].transform);
	states[0].transform = NULL;
      }
    }
    if(dstTrans)
    {
      (void )WlzFreeAffineTransform(*dstTrans);
      *dstTrans = NULL;
      if(ppCtrl->method & REC_PP_SAMPLE)
      {
	tTr = states[approach].transform;
	if((wlzErr = WlzAffineTransformPrimGet(tTr, &prim)) == WLZ_ERR_NONE)
	{
	  *dstTrans = WlzAssignAffineTransform(
		      WlzAffineTransformFromPrimVal(WLZ_TRANSFORM_2D_AFFINE,
					    prim.tx * samFac.vtX,
					    prim.ty * samFac.vtY,
					    0.0,
					    1.0, prim.theta, 0.0,
					    0.0, 0.0, 0.0,
					    0, &wlzErr), NULL);
	  (void )WlzFreeAffineTransform(states[approach].transform);
	}
	errFlag = RecErrorFromWlz(wlzErr);
      }
      else
      {
        *dstTrans = states[approach].transform;   /* Use existing linkcount */
      }
    }
    if(dstCrossC)
    {
      *dstCrossC = states[approach].correl;
    }
    if(dstIter)
    {
      *dstIter = states[approach].iteration;
    }
  }
  else
  {
    for(approach = 0; approach <= 1; ++approach)
    {
      if(states[approach].transform)
      {
	(void )WlzFreeAffineTransform(states[approach].transform);
      }
    }
  }
  (void )WlzFreeObj(ppObj0);
  (void )WlzFreeObj(ppObj1);
  if(freeObjFlag)
  {
    (void )WlzFreeObj(obj0);
    (void )WlzFreeObj(obj1);
  }
  REC_DBG((REC_DBG_REG|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecRegisterPair FX %d\n",
	   errFlag));
  return(errFlag);
}

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Builds a transform and then transforms the given
*               image object, using the given pair of centers of mass
*               and angle for the principal axis.
*               If *transfObj == givenObj, this function should still
*               work.
* \param	transf			Destination pointer for transform.
* \param	transfObj		Destination pointer for transformed
* 					object.
* \param	cMass0			Center of mass of 1st object.
* \param	cMass1			Center of mass of 2nd object.
* \param	angle			Principal angle.
* \param	givenObj		Given object fro transformation.
*/
static RecError	RecRegPrincipal(WlzAffineTransform **transf,
			       WlzObject **transfObj,
			       WlzDVertex2 cMass0, WlzDVertex2 cMass1,
			       double angle, WlzObject *givenObj)
{
  WlzAffineTransform	*tTr0 = NULL,
		*tTr1 = NULL,
		*newTr = NULL,
		*givenTr = NULL;
  WlzObject	*tObj = NULL;
  RecError	errFlag = REC_ERR_NONE;
  WlzErrorNum	wlzErr = WLZ_ERR_NONE;

  REC_DBG((REC_DBG_REG|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecRegPrincipal FE 0x%lx 0x%lx {%g %g} {%g %g} %g 0x%lx\n",
	   (unsigned long )transf, (unsigned long )transfObj,
	   cMass0.vtX, cMass0.vtY, cMass1.vtX, cMass1.vtY,
	   angle, (unsigned long )givenObj));
  tTr0 = WlzAffineTransformFromSpin((double )(cMass0.vtX),
  			            (double )(cMass0.vtY), angle, &wlzErr);
  if(wlzErr == WLZ_ERR_NONE)
  {
    tTr1 = WlzAffineTransformFromPrimVal(WLZ_TRANSFORM_2D_AFFINE,
    				        cMass0.vtX - cMass1.vtX,
			                cMass0.vtY - cMass1.vtY, 0.0,
				        1.0, 0.0, 0.0,
				        0.0, 0.0, 0.0,
				        0, &wlzErr);
  }
  if(wlzErr == WLZ_ERR_NONE)
  {
    givenTr = *transf;
    newTr = WlzAssignAffineTransform(
            WlzAffineTransformProduct(tTr1, tTr0, &wlzErr), NULL);
  }
  if(tTr0)
  {
    (void )WlzFreeAffineTransform(tTr0);
  }
  if(tTr1)
  {
    (void )WlzFreeAffineTransform(tTr1);
  }
  if(wlzErr == WLZ_ERR_NONE)
  {
    if(givenTr)
    {
      *transf = WlzAssignAffineTransform(
		WlzAffineTransformProduct(givenTr, newTr, &wlzErr), NULL);
      (void )WlzFreeAffineTransform(givenTr);
    }
    else
    {
      *transf = WlzAssignAffineTransform(newTr, NULL);
    }
  }
  if(newTr)
  {
    (void )WlzFreeAffineTransform(newTr);
  }
  if(wlzErr == WLZ_ERR_NONE)
  {
    tObj = *transfObj;
    *transfObj = WlzAssignObject(
    		 WlzAffineTransformObj(givenObj, *transf,
		 		       WLZ_INTERPOLATION_NEAREST,
				       &wlzErr), NULL);
    if(tObj)
    {
      (void )WlzFreeObj(tObj);
    }
  }
  errFlag = RecErrorFromWlz(wlzErr);
  REC_DBG((REC_DBG_REG|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecRegPrincipal FX %d\n",
	   (int )errFlag));
  return(errFlag);
}

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Builds a transform and then transforms the given
*               image object, using the given center of mass and
*               angle of rotation.
* \param	transf			Destination pointer for transform.
* \param	transfObj		Destination pointer for transformed
*					object.
* \param	cMass			Center of mass about which to
*					rotate the image object.
* \param	angle			Angle of roattion.
* \param	secObj			Object to be transformed.
*/
static RecError	RecRegRotate(WlzAffineTransform **transf, WlzObject **transfObj,
			     WlzIVertex2 cMass, double angle,
			     WlzObject *secObj)
{
  WlzAffineTransform	*tTransf0,
		*tTransf1;
  WlzErrorNum	wlzErr = WLZ_ERR_NONE;
  RecError	errFlag = REC_ERR_NONE;

  REC_DBG((REC_DBG_REG|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecRegRotate FE 0x%lx 0x%lx {%d %d} %f 0x%lx\n",
	   (unsigned long )transf, (unsigned long )transfObj,
	   cMass.vtX, cMass.vtY, angle, (unsigned long )secObj));
  tTransf0 = WlzAffineTransformFromSpin((double )(cMass.vtX),
  					(double )(cMass.vtY),
					angle, &wlzErr);
  if(wlzErr == WLZ_ERR_NONE)
  {
    if((tTransf1 = *transf) != NULL)
    {
      *transf = WlzAssignAffineTransform(
      	       WlzAffineTransformProduct(tTransf1, tTransf0, &wlzErr), NULL);
      if(tTransf1)
      {
        (void )WlzFreeAffineTransform(tTransf1);
      }
    }
    else
    {
      *transf = WlzAssignAffineTransform(tTransf0, NULL);
    }
  }
  if(tTransf0)
  {
    WlzFreeAffineTransform(tTransf0);
  }
  if(*transfObj)
  {
    (void )WlzFreeObj(*transfObj);
  }
  if(wlzErr == WLZ_ERR_NONE)
  {
    *transfObj = WlzAssignObject(
    		 WlzAffineTransformObj(secObj, *transf,
		 		       WLZ_INTERPOLATION_NEAREST,
				       &wlzErr), NULL);
  }
  errFlag = RecErrorFromWlz(wlzErr);
  REC_DBG((REC_DBG_REG|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecRegRotate FX %d\n",
	   (int )errFlag));
  return(errFlag);
}

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Builds a transform and then transforms the given
*               object, using the given translation.
* \param	transf			Destination pointer for transform.
* \param	transfObj		Destination pointer for transformed
* 					object.
* \param	trans			Given translation.
* \param	secObj			Oobject to be transformed.
*/
static RecError	RecRegTranslate(WlzAffineTransform **transf,
			        WlzObject **transfObj,
			        WlzDVertex2 trans,
			        WlzObject *secObj)
{
  WlzAffineTransform	*tTransf0,
		*tTransf1;
  WlzErrorNum	wlzErr = WLZ_ERR_NONE;
  RecError	errFlag = REC_ERR_NONE;

  REC_DBG((REC_DBG_REG|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecRegTranslate FE 0x%lx 0x%lx {%f %f} 0x%lx\n",
	   (unsigned long )transf, (unsigned long )transfObj,
	   trans.vtX, trans.vtY, (unsigned long )secObj));
  tTransf0 = WlzAffineTransformFromPrimVal(WLZ_TRANSFORM_2D_AFFINE,
  					   trans.vtX, trans.vtY, 0.0,
					   1.0, 0.0, 0.0,
					   0.0, 0.0, 0.0,
					   0, &wlzErr);
  if(wlzErr == WLZ_ERR_NONE)
  {
    if((tTransf1 = *transf) != NULL)
    {
      *transf = WlzAssignAffineTransform(
      		WlzAffineTransformProduct(tTransf1, tTransf0, &wlzErr), NULL);
      if(tTransf1)
      {
	(void )WlzFreeAffineTransform(tTransf1);
      }
    }
    else
    {
      *transf =  WlzAssignAffineTransform(tTransf0, NULL);
    }
  }
  if(tTransf0)
  {
    WlzFreeAffineTransform(tTransf0);
  }
  if(*transfObj)
  {
    (void )WlzFreeObj(*transfObj);
  }
  if(wlzErr == WLZ_ERR_NONE)
  {
    *transfObj = WlzAssignObject(
    		 WlzAffineTransformObj(secObj, *transf,
		 		       WLZ_INTERPOLATION_NEAREST,
				       &wlzErr), NULL);
  }
  errFlag = RecErrorFromWlz(wlzErr);
  REC_DBG((REC_DBG_REG|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecRegTranslate FX %d\n",
	   (int )errFlag));
  return(errFlag);
}

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Calculates a registration transform from the given
*               vector of tie-points. If required an error distance
*               and/or cross-correlation value are calculated.
*               The cross-correlation is the mean of cross-correlations
*               about the tie points.
* \param	dstTr			Destination pointer for the
*                                       calculated transform.
* \param	dstED			Destination pointer for the error
*                                       distance. Not calculated if
*                                       NULL.
* \param	dstCC			Destination pointer for the
*					cross-correlation value. Not
*					calculated if NULL.
* \param	tppVec			Vector of tie points.
* \param	tppCount		Number of tie point pairs.
* \param	obj0			First object (only used if
*                                       cross-correlation required).
* \param	obj1			Second object (only used if
*                                       cross-correlation required).
*/
RecError	RecRegisterTiePoints(WlzAffineTransform **dstTr,
				     double *dstED, double *dstCC,
				     RecTiePointPair *tppVec, int tppCount,
				     WlzObject *obj0, WlzObject *obj1)
{
  int		tppIdx;
  double	newCC,
  		sumCC = 0.0,
		sumED = 0.0;
  WlzAffineTransform *newTr = NULL;
  WlzDVertex2	dVtx,
  		trVtx;
  WlzDVertex2	*vtxVec0 = NULL,
  		*vtxVec1 = NULL;
  WlzObject	*trObj = NULL;
  WlzIVertex2	roiCtr,
		roiSz;
  WlzErrorNum	wlzErr = WLZ_ERR_NONE;
  RecError	errFlag = REC_ERR_NONE;
  WlzAffineTransformPrim prim;
  RecPPControl	ppCtrl;

  REC_DBG((REC_DBG_REG|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecRegisterTiePoints FE 0x%lx 0x%lx 0x%lx %d 0x%lx 0x%lx\n",
	   (unsigned long )dstTr, (unsigned long )dstED,
	   (unsigned long )tppVec, tppCount,
	   (unsigned long )obj0, (unsigned long )obj1));
  if((dstTr || dstED) && tppVec && (tppCount > 0))
  {
    if(((vtxVec0 = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) *
       					   tppCount)) == NULL) ||
       ((vtxVec1 = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) *
       					   tppCount)) == NULL))

    {
      errFlag = REC_ERR_MALLOC;
    }
    else
    {
      for(tppIdx = 0; tppIdx < tppCount; ++tppIdx)
      {
        (vtxVec0 + tppIdx)->vtX = (tppVec + tppIdx)->first.vtX;
        (vtxVec0 + tppIdx)->vtY = (tppVec + tppIdx)->first.vtY;
        (vtxVec1 + tppIdx)->vtX = (tppVec + tppIdx)->second.vtX;
        (vtxVec1 + tppIdx)->vtY = (tppVec + tppIdx)->second.vtY;
        REC_DBG((REC_DBG_REG|REC_DBG_LVL_2),
	 	("RecRegisterTiePoints 01 %d %g %g %g %g\n",
		 tppIdx,
		 (tppVec + tppIdx)->first.vtX,
		 (tppVec + tppIdx)->first.vtY,
		 (tppVec + tppIdx)->second.vtX,
		 (tppVec + tppIdx)->second.vtY));
      }
      if((newTr = WlzAffineTransformLSq2D(tppCount, vtxVec0,
					  tppCount, vtxVec1,
					  0, NULL,
					  WLZ_TRANSFORM_2D_REG,
					  &wlzErr)) == NULL)
      {
        errFlag = RecErrorFromWlz(wlzErr);
      }
      else
      {
	(void )WlzAffineTransformPrimGet(newTr, &prim);
	REC_DBG((REC_DBG_REG|REC_DBG_LVL_1),
	 	("RecRegisterTiePoints 02 %g %g %g\n",
		 prim.tx, prim.ty, prim.theta));
	if(dstED || dstCC)
	{
	  sumCC = 0.0;
	  sumED = 0.0;
	  roiSz.vtX = 16;
	  roiSz.vtY = 16;
	  if(dstCC)
	  {
	    trObj = WlzAssignObject(
	    	    WlzAffineTransformObj(obj1, newTr,
					  WLZ_INTERPOLATION_NEAREST,
					  &wlzErr), NULL);
	    if(trObj == NULL)
	    {
	      errFlag = RecErrorFromWlz(wlzErr);
	    }
	  }
	  tppIdx = 0;
	  while((tppIdx < tppCount) && (errFlag == REC_ERR_NONE))
	  {
	    trVtx = WlzAffineTransformVertexD2(newTr, (tppVec + tppIdx)->second,
					       &wlzErr);
	    if(wlzErr == WLZ_ERR_NONE)
	    {
	      if(dstED)
	      {
		dVtx.vtX = (tppVec + tppIdx)->first.vtX - trVtx.vtX;
		dVtx.vtY = (tppVec + tppIdx)->first.vtY - trVtx.vtY;
		sumED += WLZ_VTX_2_SQRLEN(dVtx);
	      }
	      if(dstCC)
	      {
		roiCtr.vtX = WLZ_NINT((((tppVec + tppIdx)->first.vtX) +
				       trVtx.vtX) / 2.0);
		roiCtr.vtY = WLZ_NINT((((tppVec + tppIdx)->first.vtY) +
				       trVtx.vtY) / 2.0);
		ppCtrl.method = REC_PP_NONE;
	        errFlag = RecCrossCorrelateROI(&newCC, NULL, NULL,
					       obj0, trObj, roiCtr, roiSz,
					       &ppCtrl);
	        if((errFlag == REC_ERR_NONE) && (newCC > DBL_EPSILON))
		{
		  sumCC += newCC;
		  REC_DBG((REC_DBG_REG|REC_DBG_LVL_2),
			  ("RecRegisterTiePoints 03 %g\n",
			   newCC));
		}
	      }
	      ++tppIdx;
	    }
	    else
	    {
	      errFlag = RecErrorFromWlz(wlzErr);
	    }
	  }
	  REC_DBG((REC_DBG_REG|REC_DBG_LVL_1),
		  ("RecRegisterTiePoints 04 %g %g\n",
		   (sumED > DBL_EPSILON)? sqrt(sumED): 0.0,
		   (tppCount > 0)? (sumCC / tppCount): 0.0));
	}
      }
    }
  }
  if(trObj)
  {
    (void )WlzFreeObj(trObj);
  }
  if(vtxVec0)
  {
    free(vtxVec0);
  }
  if(vtxVec1)
  {
    free(vtxVec1);
  }
  if(errFlag == REC_ERR_NONE)
  {
    if(dstTr)
    {
      *dstTr = newTr;
    }
    if(dstCC)
    {
      *dstCC = (tppCount > 0)? (sumCC / tppCount): 0.0;
    }
    if(dstED)
    {
      *dstED = (sumED > DBL_EPSILON)? sqrt(sumED): 0.0;
    }
  }
  else
  {
    if(newTr)
    {
      WlzFreeAffineTransform(newTr);
    }
  }
  REC_DBG((REC_DBG_REG|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecRegisterTiePoints FX %d\n",
	   (int )errFlag));
  return(errFlag);
}
