#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzImageArithmetic.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions for performing binary (ie two objects)
*		arithmetic on a pair of domain objects with grey
*		values.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
*		bill@hgu.mrc.ac.uk Added magnitude.
************************************************************************/
#include <stdarg.h>
#include <float.h>
#include <Wlz.h>

typedef	void (*WlzBinaryOperatorFn)(WlzGreyP, WlzGreyP, int);

/************************************************************************
* Function:	WlzBufXXX(I)|(D)					*
* Returns:	void							*
* Purpose:	Perform self-expanatory binary operation on pair of	*
*		buffers.						*
* Global refs:	-							*
* Parameters:	WlzGreyP gP1:		Second buffer for input, also	*
*					has the result.			*
*		WlzGreyP gP0:		First input buffer.		*
*		int count:		Number of data in buffer.	*
************************************************************************/
static void	WlzBufAddI(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  int		*iP0,
  		*iP1;
  iP0 = gP0.inp;
  iP1 = gP1.inp;
  while(count-- > 0)
  {
    *iP1 = *iP0++ + *iP1;
    ++iP1;
  }
}

static void	WlzBufAddD(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  double	*dP0,
  		*dP1;
  dP0 = gP0.dbp;
  dP1 = gP1.dbp;
  while(count-- > 0)
  {
    *dP1 = *dP1 + *dP0++;
    ++dP1;
  }
}

static void	WlzBufSubI(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  int		*iP0,
  		*iP1;
  iP0 = gP0.inp;
  iP1 = gP1.inp;
  while(count-- > 0)
  {
    *iP1 = *iP0++ - *iP1;
    ++iP1;
  }
}

static void	WlzBufSubD(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  double	*dP0,
  		*dP1;
  dP0 = gP0.dbp;
  dP1 = gP1.dbp;
  while(count-- > 0)
  {
    *dP1 = *dP0++ - *dP1;
    ++dP1;
  }
}

static void	WlzBufMulI(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  int		*iP0,
  		*iP1;
  iP0 = gP0.inp;
  iP1 = gP1.inp;
  while(count-- > 0)
  {
    *iP1 = *iP0++ * *iP1;
    ++iP1;
  }
}

static void	WlzBufMulD(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  double	*dP0,
  		*dP1;
  dP0 = gP0.dbp;
  dP1 = gP1.dbp;
  while(count-- > 0)
  {
    *dP1 = *dP0++ * *dP1;
    ++dP1;
  }
}

static void	WlzBufDivI(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  int		*iP0,
  		*iP1;
  iP0 = gP0.inp;
  iP1 = gP1.inp;
  while(count-- > 0)
  {
    *iP1 = (*iP1)? (*iP0 / *iP1): (*iP0);
    ++iP0;
    ++iP1;
  }
}

static void	WlzBufDivD(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  double	*dP0,
  		*dP1;
  dP0 = gP0.dbp;
  dP1 = gP1.dbp;
  while(count-- > 0)
  {
    *dP1 = (fabs(*dP1) > DBL_EPSILON)? (*dP0 / *dP1): (*dP0);
    ++dP0;
    ++dP1;
  }
}

static void	WlzBufModI(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  int		*iP0,
  		*iP1;
  iP0 = gP0.inp;
  iP1 = gP1.inp;
  while(count-- > 0)
  {
    *iP1 = (*iP1)? (*iP0 % *iP1): 0;
    ++iP0;
    ++iP1;
  }
}

static void	WlzBufModD(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  double	*dP0,
  		*dP1;
  dP0 = gP0.dbp;
  dP1 = gP1.dbp;
  while(count-- > 0)
  {
    *dP1 = (fabs(*dP1) > DBL_EPSILON)? (*dP0 / *dP1): (0.0);
    ++dP0;
    ++dP1;
  }
}

static void	WlzBufEQI(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  int		*iP0,
  		*iP1;
  iP0 = gP0.inp;
  iP1 = gP1.inp;
  while(count-- > 0)
  {
    *iP1 = *iP0 == *iP1;
    ++iP0;
    ++iP1;
  }
}

static void	WlzBufEQD(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  double	*dP0,
  		*dP1;

  dP0 = gP0.dbp;
  dP1 = gP1.dbp;
  while(count-- > 0)
  {
    *dP1 = (fabs(*dP0 - *dP1) > DBL_EPSILON)? (1.0): (0.0);
    ++dP0;
    ++dP1;
  }
}

static void	WlzBufNEI(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  int		*iP0,
  		*iP1;
  iP0 = gP0.inp;
  iP1 = gP1.inp;
  while(count-- > 0)
  {
    *iP1 = *iP0 != *iP1;
    ++iP0;
    ++iP1;
  }
}

static void	WlzBufNED(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  double	*dP0,
  		*dP1;

  dP0 = gP0.dbp;
  dP1 = gP1.dbp;
  while(count-- > 0)
  {
    *dP1 = (fabs(*dP0 - *dP1) < DBL_EPSILON)? (1.0): (0.0);
    ++dP0;
    ++dP1;
  }
}

static void	WlzBufGTI(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  int		*iP0,
  		*iP1;
  iP0 = gP0.inp;
  iP1 = gP1.inp;
  while(count-- > 0)
  {
    *iP1 = *iP0 > *iP1;
    ++iP0;
    ++iP1;
  }
}

static void	WlzBufGTD(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  double	*dP0,
  		*dP1;

  dP0 = gP0.dbp;
  dP1 = gP1.dbp;
  while(count-- > 0)
  {
    *dP1 = (fabs(*dP0 - *dP1) > 0.0)? (1.0): (0.0);
    ++dP0;
    ++dP1;
  }
}
static void	WlzBufGEI(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  int		*iP0,
  		*iP1;
  iP0 = gP0.inp;
  iP1 = gP1.inp;
  while(count-- > 0)
  {
    *iP1 = *iP0 >= *iP1;
    ++iP0;
    ++iP1;
  }
}
static void	WlzBufLTI(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  int		*iP0,
  		*iP1;
  iP0 = gP0.inp;
  iP1 = gP1.inp;
  while(count-- > 0)
  {
    *iP1 = *iP0 < *iP1;
    ++iP0;
    ++iP1;
  }
}
static void	WlzBufLTD(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  double	*dP0,
  		*dP1;

  dP0 = gP0.dbp;
  dP1 = gP1.dbp;
  while(count-- > 0)
  {
    *dP1 = (fabs(*dP0 - *dP1) < 0.0)? (1.0): (0.0);
    ++dP0;
    ++dP1;
  }
}

static void	WlzBufLEI(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  int		*iP0,
  		*iP1;
  iP0 = gP0.inp;
  iP1 = gP1.inp;
  while(count-- > 0)
  {
    *iP1 = *iP0 <= *iP1;
    ++iP0;
    ++iP1;
  }
}

static void	WlzBufAndI(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  int		*iP0,
  		*iP1;
  iP0 = gP0.inp;
  iP1 = gP1.inp;
  while(count-- > 0)
  {
    *iP1 = *iP0 & *iP1;
    ++iP0;
    ++iP1;
  }
}

static void	WlzBufAndD(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  double	*dP0,
  		*dP1;

  dP0 = gP0.dbp;
  dP1 = gP1.dbp;
  while(count-- > 0)
  {
    *dP1 = (unsigned int)*dP0 & (unsigned int )*dP1;
    ++dP0;
    ++dP1;
  }
}

static void	WlzBufOrI(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  int		*iP0,
  		*iP1;
  iP0 = gP0.inp;
  iP1 = gP1.inp;
  while(count-- > 0)
  {
    *iP1 = *iP0 | *iP1;
    ++iP0;
    ++iP1;
  }
}

static void	WlzBufOrD(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  double	*dP0,
  		*dP1;

  dP0 = gP0.dbp;
  dP1 = gP1.dbp;
  while(count-- > 0)
  {
    *dP1 = (unsigned int)*dP0 | (unsigned int )*dP1;
    ++dP0;
    ++dP1;
  }
}

static void	WlzBufXOrI(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  int		*iP0,
  		*iP1;
  iP0 = gP0.inp;
  iP1 = gP1.inp;
  while(count-- > 0)
  {
    *iP1 = *iP0 ^ *iP1;
    ++iP0;
    ++iP1;
  }
}

static void	WlzBufXOrD(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  double	*dP0,
  		*dP1;

  dP0 = gP0.dbp;
  dP1 = gP1.dbp;
  while(count-- > 0)
  {
    *dP1 = (unsigned int)*dP0 ^ (unsigned int )*dP1;
    ++dP0;
    ++dP1;
  }
}

static void	WlzBufMaxI(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  int		*iP0,
  		*iP1;
  iP0 = gP0.inp;
  iP1 = gP1.inp;
  while(count-- > 0)
  {
    *iP1 = (*iP0 > *iP1)? *iP0: *iP1;
    ++iP0;
    ++iP1;
  }
}
static void	WlzBufMaxD(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  double	*dP0,
  		*dP1;
  dP0 = gP0.dbp;
  dP1 = gP1.dbp;
  while(count-- > 0)
  {
    *dP1 = (*dP0 > *dP1)? *dP0: *dP1;
    ++dP0;
    ++dP1;
  }
}

static void	WlzBufMinI(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  int		*iP0,
  		*iP1;
  iP0 = gP0.inp;
  iP1 = gP1.inp;
  while(count-- > 0)
  {
    *iP1 = (*iP0 < *iP1)? *iP0: *iP1;
    ++iP0;
    ++iP1;
  }
}
static void	WlzBufMinD(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  double	*dP0,
  		*dP1;
  dP0 = gP0.dbp;
  dP1 = gP1.dbp;
  while(count-- > 0)
  {
    *dP1 = (*dP0 < *dP1)? *dP0: *dP1;
    ++dP0;
    ++dP1;
  }
}

static void	WlzBufMagI(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  double	dV0;
  int		*iP0,
  		*iP1;
  iP0 = gP0.inp;
  iP1 = gP1.inp;
  while(count-- > 0)
  {
    if(*iP0 || *iP1)
    {
      dV0 = (*iP0 * *iP0) + (*iP1 * *iP1);
      dV0 = sqrt(dV0);
      *iP1 = WLZ_NINT(dV0);
    }
    else
    {
      *iP1 = 0;
    }
    ++iP0;
    ++iP1;
  }
}
static void	WlzBufMagD(WlzGreyP gP1, WlzGreyP gP0, int count)
{
  double	*dP0,
  		*dP1;
  dP0 = gP0.dbp;
  dP1 = gP1.dbp;
  while(count-- > 0)
  {
    if((*dP1 = (*dP0 * *dP0) + (*dP1 * *dP1)) > DBL_EPSILON)
    {
      *dP1= sqrt(*dP1);
    }
    else
    {
      *dP1 = 0.0;
    }
    ++dP0;
    ++dP1;
  }
}

/************************************************************************
* Function:	WlzBinaryOperatorFnSet					*
* Returns:	WlzBinaryOperatorFn:	Pointer to appropriate function	*
*					or NULL on error.		*
* Purpose:	Finds the appropriate function for the given binary	*
*		operator and grey type.					*
* Global refs:	-							*
* Parameters:	WlzGreyType gType:	Given grey type.		*
*		WlzBinaryOperatorType op: Given operator.		*
*		WlzErrorNum *dstErr:	Destination error pointer, may	*
*					be NULL.			*
************************************************************************/
static WlzBinaryOperatorFn WlzBinaryOperatorFnSet(WlzGreyType gType,
						  WlzBinaryOperatorType op,
						  WlzErrorNum *dstErr)
{
  WlzBinaryOperatorFn binOpFn = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((gType != WLZ_GREY_INT) && (gType != WLZ_GREY_DOUBLE))
  {
    errNum = WLZ_ERR_GREY_TYPE;
  }
  else
  {
    switch(op)
    {
      case WLZ_BO_ADD:
	binOpFn = (gType == WLZ_GREY_INT)? WlzBufAddI: WlzBufAddD;
	break;
      case WLZ_BO_SUBTRACT:
	binOpFn = (gType == WLZ_GREY_INT)? WlzBufSubI: WlzBufSubD;
	break;
      case WLZ_BO_MULTIPLY:
	binOpFn = (gType == WLZ_GREY_INT)? WlzBufMulI: WlzBufMulD;
	break;
      case WLZ_BO_DIVIDE:
	binOpFn = (gType == WLZ_GREY_INT)? WlzBufDivI: WlzBufDivD;
	break;
      case WLZ_BO_MODULUS:
	binOpFn = (gType == WLZ_GREY_INT)? WlzBufModI: WlzBufModD;
	break;
      case WLZ_BO_EQ:
	binOpFn = (gType == WLZ_GREY_INT)? WlzBufEQI: WlzBufEQD;
	break;
      case WLZ_BO_NE:
	binOpFn = (gType == WLZ_GREY_INT)? WlzBufNEI: WlzBufNED;
	break;
      case WLZ_BO_GT:
	binOpFn = (gType == WLZ_GREY_INT)? WlzBufGTI: WlzBufGTD;
	break;
      case WLZ_BO_GE:
	binOpFn = (gType == WLZ_GREY_INT)? WlzBufGEI: WlzBufGTD;
	break;
      case WLZ_BO_LT:
	binOpFn = (gType == WLZ_GREY_INT)? WlzBufLTI: WlzBufLTD;
	break;
      case WLZ_BO_LE:
	binOpFn = (gType == WLZ_GREY_INT)? WlzBufLEI: WlzBufLTD;
	break;
      case WLZ_BO_AND:
	binOpFn = (gType == WLZ_GREY_INT)? WlzBufAndI: WlzBufAndD;
	break;
      case WLZ_BO_OR:
	binOpFn = (gType == WLZ_GREY_INT)? WlzBufOrI: WlzBufOrD;
	break;
      case WLZ_BO_XOR:
	binOpFn = (gType == WLZ_GREY_INT)? WlzBufXOrI: WlzBufXOrD;
	break;
      case WLZ_BO_MAX:
	binOpFn = (gType == WLZ_GREY_INT)? WlzBufMaxI: WlzBufMaxD;
	break;
      case WLZ_BO_MIN:
	binOpFn = (gType == WLZ_GREY_INT)? WlzBufMinI: WlzBufMinD;
	break;
      case WLZ_BO_MAGNITUDE:
	binOpFn = (gType == WLZ_GREY_INT)? WlzBufMagI: WlzBufMagD;
	break;
      default:
	errNum = WLZ_ERR_BINARY_OPERATOR_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(binOpFn);
}

/************************************************************************
* Function:	WlzImageArithmetic2D					*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Performs  binary (ie two objects) arithmetic on a pair	*
*		of 2D domain objects.					*
*		If the overwrite flag is set and the grey values of	*
*		the object to be overwritten are of the wrong type then	*
*		the returned object has new values as though the	*
*		overwrite flag was not set.				*
* Global refs:	-							*
* Parameters:	WlzObject *obj0:	First object.			*
*		WlzObject *obj1:	Second object.			*
*		WlzObject *obj2:	Intersection of the 1st and 2nd	*
*					objects, but values still to be	*
*					filled in.			*
*		WlzBinaryOperatorType op: Binary operator.		*
*		int overwrite:		Allow the destination object	*
*					to share values with one of	*
*					the given objects if non zero.	*
*					  0: No values shared.		*
*					  1: Values shared with obj0.	*
*					  2: Values shared with obj1.	*
*					  < 0 || > 2: Error condition.	*
************************************************************************/
static WlzErrorNum WlzImageArithmetic2D(WlzObject *obj0, WlzObject *obj1,
				        WlzObject *obj2,
				        WlzBinaryOperatorType op,
				        int overwrite)
{
  int		idx,
  		tI0;
  WlzValues	tVal;
  WlzObject	*obj[3];
  WlzGreyType	gType[4];
  WlzPixelV	bgd[3];
  WlzGreyP	buf[2];
  WlzIntervalWSpace iWsp[3];
  WlzGreyWSpace	gWsp[3];
  WlzBinaryOperatorFn binOpFn;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  obj[0] = obj[1] = NULL;
  obj[2] = obj2;
  buf[0].inp = buf[1].inp = NULL;
  obj[0] = WlzMakeMain(WLZ_2D_DOMAINOBJ, obj2->domain, obj0->values,
		       NULL, NULL, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    obj[1] = WlzMakeMain(WLZ_2D_DOMAINOBJ, obj2->domain, obj1->values,
			 NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Get background value and grey types */
    gType[0] = WlzGreyTableTypeToGreyType(obj0->values.core->type,
					  &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      gType[1] = WlzGreyTableTypeToGreyType(obj1->values.core->type,
					    &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      bgd[0] = WlzGetBackground(obj0, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      bgd[1] = WlzGetBackground(obj1, &errNum);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Promote grey type */
    if((gType[0] == WLZ_GREY_DOUBLE) || (gType[1] == WLZ_GREY_DOUBLE))
    {
      gType[2] = WLZ_GREY_DOUBLE;
      gType[3] = WLZ_GREY_DOUBLE;
    }
    else if((gType[0] == WLZ_GREY_FLOAT) || (gType[1] == WLZ_GREY_FLOAT))
    {
      gType[2] = WLZ_GREY_FLOAT;
      gType[3] = WLZ_GREY_DOUBLE;
    }
    else if((gType[0] == WLZ_GREY_INT) || (gType[1] == WLZ_GREY_INT))
    {
      gType[2] = WLZ_GREY_INT;
      gType[3] = WLZ_GREY_INT;
    }
    else if((gType[0] == WLZ_GREY_SHORT) || (gType[1] == WLZ_GREY_SHORT))
    {
      gType[2] = WLZ_GREY_SHORT;
      gType[3] = WLZ_GREY_INT;
    }
    else
    {
      gType[3] = WLZ_GREY_INT;
      switch(op)
      {
	case WLZ_BO_ADD:	/* FALLTHROGH */
	case WLZ_BO_SUBTRACT:	/* FALLTHROGH */
	case WLZ_BO_MULTIPLY:	/* FALLTHROGH */
	case WLZ_BO_DIVIDE:	/* FALLTHROGH */
	case WLZ_BO_MAGNITUDE:
          gType[2] = WLZ_GREY_SHORT;
	  break;
	case WLZ_BO_MODULUS:	/* FALLTHROGH */
	case WLZ_BO_EQ:		/* FALLTHROGH */
	case WLZ_BO_NE:		/* FALLTHROGH */
	case WLZ_BO_GT:		/* FALLTHROGH */
	case WLZ_BO_GE:		/* FALLTHROGH */
	case WLZ_BO_LT:		/* FALLTHROGH */
	case WLZ_BO_LE:		/* FALLTHROGH */
	case WLZ_BO_AND:	/* FALLTHROGH */
	case WLZ_BO_OR:		/* FALLTHROGH */
	case WLZ_BO_XOR:	/* FALLTHROGH */
	case WLZ_BO_MAX:	/* FALLTHROGH */
	case WLZ_BO_MIN:
          gType[2] = WLZ_GREY_UBYTE;
	  break;
	default:
	  errNum = WLZ_ERR_BINARY_OPERATOR_TYPE;
	  break;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    binOpFn = WlzBinaryOperatorFnSet(gType[3], op, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    (void )WlzValueConvertPixel(bgd + 0, bgd[0], gType[3]);
    (void )WlzValueConvertPixel(bgd + 1, bgd[1], gType[3]);
    bgd[2] = bgd[1];
    if(gType[3] == WLZ_GREY_INT)
    {
      buf[0].inp = &(bgd[0].v.inv);
      buf[1].inp = &(bgd[2].v.inv);
    }
    else
    {
      buf[0].dbp = &(bgd[0].v.dbv);
      buf[1].dbp = &(bgd[2].v.dbv);
    }
    binOpFn(buf[1], buf[0], 1);
    if(((overwrite == 1) && (gType[0] != gType[2])) ||
       ((overwrite == 2) && (gType[1] != gType[2])))
    {
      overwrite = 0;    /* Can only overwrite if the grey types are the same */
    }
    switch(overwrite)
    {
      case 0:
	tVal.v = WlzNewValueTb(obj2,
			       WlzGreyTableType(WLZ_GREY_TAB_RAGR,
						gType[2], NULL),
			       bgd[2], &errNum);
        break;
      case 1:
	tVal = obj[0]->values;
	break;
      case 2:
	tVal = obj[1]->values;
	break;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      obj2->values = WlzAssignValues(tVal, NULL);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Make buffers */
    tI0 = obj2->domain.i->lastkl - obj2->domain.i->kol1 + 1;
    if(gType[3] == WLZ_GREY_INT)
    {
      tI0 *= sizeof(int);
      if(((buf[0].inp = (int *)AlcMalloc(tI0)) == NULL) ||
	 ((buf[1].inp = (int *)AlcMalloc(tI0)) == NULL))
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    else
    {
      tI0 *= sizeof(double);
      if(((buf[0].dbp = (double *)AlcMalloc(tI0)) == NULL) ||
	 ((buf[1].dbp = (double *)AlcMalloc(tI0)) == NULL))
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    idx = 0;
    while((idx <= 2) && (errNum == WLZ_ERR_NONE))
    {
      errNum = WlzInitGreyScan(obj[idx], iWsp + idx, gWsp + idx);
      ++idx;
    }
  }
  while((errNum == WLZ_ERR_NONE) &&
        ((errNum = WlzNextGreyInterval(iWsp + 0)) == WLZ_ERR_NONE) &&
	((errNum = WlzNextGreyInterval(iWsp + 1)) == WLZ_ERR_NONE) &&
	((errNum = WlzNextGreyInterval(iWsp + 2)) == WLZ_ERR_NONE))
  {
    tI0 = iWsp[0].rgtpos - iWsp[0].lftpos + 1;
    WlzValueCopyGreyToGrey(buf[0], 0, gType[3],
			   gWsp[0].u_grintptr, 0, gWsp[0].pixeltype,
			   tI0);
    WlzValueCopyGreyToGrey(buf[1], 0, gType[3],
			   gWsp[1].u_grintptr, 0, gWsp[1].pixeltype,
			   tI0);
    binOpFn(buf[1], buf[0], tI0);
    WlzValueCopyGreyToGrey(gWsp[2].u_grintptr, 0, gWsp[2].pixeltype,
			   buf[1], 0, gType[3],
			   tI0);
  }
  if(errNum == WLZ_ERR_EOO)
  {
    errNum = WLZ_ERR_NONE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(overwrite != 0)
    {
      errNum = WlzSetBackground(obj2, bgd[2]);
    }
  }
  if(buf[0].inp)
  {
    AlcFree(buf[0].inp);
  }
  if(buf[1].inp)
  {
    AlcFree(buf[1].inp);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzImageArithmetic3D					*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Performs  binary (ie two objects) arithmetic on a pair	*
*		of 3D domain objects.					*
*		If the overwrite flag is set and the grey values of	*
*		the object to be overwritten are of the wrong type then	*
*		the returned object has new values as though the	*
*		overwrite flag was not set.				*
* Global refs:	-							*
* Parameters:	WlzObject *obj0:	First object.			*
*		WlzObject *obj1:	Second object.			*
*		WlzObject *obj2:	Intrsetion of the 1st and 2nd	*
*					objects, but values still to be	*
*					filled in.			*
*		WlzBinaryOperatorType op: Binary operator.		*
*		int overwrite:		Allow the destination object	*
*					to share values with one of	*
*					the given objects if non zero.	*
*					  0: No values shared.		*
*					  1: Values shared with obj0.	*
*					  2: Values shared with obj1.	*
*					  < 0 || > 2: Error condition.	*
************************************************************************/
static WlzErrorNum WlzImageArithmetic3D(WlzObject *obj0, WlzObject *obj1,
				        WlzObject *obj2,
				        WlzBinaryOperatorType op,
				        int overwrite)
{
  int		tI0,
  		tI1,
		oIdx,
		nPlanes,
		bgdFlag = 0;
  int		pIdx[3],
  		vIdx[3];
  WlzObject	*obj[3],
  		*obj2D[3];
  WlzPixelV	bgd[3];
  WlzPlaneDomain *pDom[3];
  WlzVoxelValues *vVal[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if ((obj0->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN) ||
           (obj1->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN))
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if ((obj0->values.core->type != WLZ_VOXELVALUETABLE_GREY) ||
           (obj1->values.core->type != WLZ_VOXELVALUETABLE_GREY))
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    obj[0] = obj0;
    obj[1] = obj1;
    obj[2] = WlzIntersect2(obj0, obj1, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(obj[2]->type != WLZ_EMPTY_OBJ)
    {
      pDom[2] = obj[2]->domain.p;
      vVal[2] = obj[2]->values.vox;
      pIdx[2] = 0;
      vIdx[2] = 0;
      tI0 = obj[2]->domain.p->plane1;
      tI1 = obj[2]->values.vox->plane1;
      oIdx = 2;
      while((oIdx-- > 0) && (errNum == WLZ_ERR_NONE))
      {
        pDom[oIdx] = obj[oIdx]->domain.p;
	vVal[oIdx] = obj[oIdx]->values.vox;
	pIdx[oIdx] = pDom[oIdx]->plane1 - tI0;
	vIdx[oIdx] = vVal[oIdx]->plane1 - tI1;
	bgd[oIdx] = WlzGetBackground(obj[oIdx], &errNum);
        
      }
      if(errNum == WLZ_ERR_NONE)
      {
	nPlanes = pDom[2]->lastpl - pDom[2]->plane1 + 1;
	pIdx[2] = 0;
	while((pIdx[2] < nPlanes) && (errNum == WLZ_ERR_NONE))
	{
	  oIdx = 2;
	  obj2D[2] = obj2D[1] = obj2D[0] = NULL;
	  while((oIdx-- > 0) && (errNum == WLZ_ERR_NONE))
	  {
	    obj2D[oIdx] = WlzMakeMain(WLZ_2D_DOMAINOBJ,
	    			      *(pDom[oIdx]->domains + pIdx[oIdx]),
	    			      *(vVal[oIdx]->values + vIdx[oIdx]),
				      NULL, NULL, &errNum);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzImageArithmetic2D(obj2D[0], obj2D[1], obj2D[2],
					  op, overwrite);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if((bgdFlag == 0) && (obj2D[2]->type == WLZ_2D_DOMAINOBJ))
	    {
      	      bgd[2] = WlzGetBackground(obj0, &errNum);
	      bgdFlag = 1;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    WlzFreeValues(*(vVal[2]->values + pIdx[2]));
	    *(vVal[2]->values + pIdx[2]) = WlzAssignValues(obj2D[2]->values,
	    						   NULL);
	  }
	  if(obj2D[0])
	  {
	    (void )WlzFreeObj(obj2D[0]);
	  }
	  if(obj2D[1])
	  {
	    (void )WlzFreeObj(obj2D[1]);
	  }
	  if(obj2D[2])
	  {
	    (void )WlzFreeObj(obj2D[2]);
	  }
	  ++pIdx[0];
	  ++vIdx[0];
	  ++pIdx[1];
	  ++vIdx[1];
	  ++pIdx[2];
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzSetBackground(obj[2], bgd[2]);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzImageArithmetic					*
* Returns:	WlzObject *:		New object or NULL on error.	*
* Purpose:	Performs  binary (ie two objects) arithmetic on a pair	*
*		of domain objects.					*
*		If the overwrite flag is set and the grey values of	*
*		the object to be overwritten are of the wrong type then	*
*		the returned object has new values as though the	*
*		overwrite flag was not set.				*
* Global refs:	-							*
* Parameters:	WlzObject *obj0:	First object.			*
*		WlzObject *obj1:	Second object.			*
*		WlzBinaryOperatorType op: Binary operator.		*
*		int overwrite:		Allow the destination object	*
*					to share values with one of	*
*					the given objects if non zero.	*
*					  0: No values shared.		*
*					  1: Values shared with obj0.	*
*					  2: Values shared with obj1.	*
*					  < 0 || > 2: Error condition.	*
*		WlzErrorNum *dstErr:	Destination error pointer, may	*
*					be NULL.			*
************************************************************************/
WlzObject	*WlzImageArithmetic(WlzObject *obj0, WlzObject *obj1,
				    WlzBinaryOperatorType op,
				    int overwrite,
				    WlzErrorNum *dstErr)
{
  WlzObject	*obj2 = NULL,
  		*dstObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((obj0 == NULL) || (obj1 == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;;
  }
  else if(obj0->type != obj1->type)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((overwrite < 0) || (overwrite > 2))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    switch(obj0->type)
    {
      case WLZ_EMPTY_OBJ:
        dstObj = WlzMakeEmpty(&errNum);
	break;
      case WLZ_2D_DOMAINOBJ: /* FALLTHROGH */
      case WLZ_3D_DOMAINOBJ: /* FALLTHROGH */
	if((obj0->domain.core == NULL) || (obj1->domain.core == NULL))
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if((obj0->values.core == NULL) || (obj1->values.core == NULL))
	{
	  errNum = WLZ_ERR_VALUES_NULL;
	}
	else
	{
	  obj2 = WlzIntersect2(obj0, obj1, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if(obj2->type != WLZ_EMPTY_OBJ)
	    {
	      switch(obj0->type)
	      {
		case WLZ_2D_DOMAINOBJ:
		  errNum =  WlzImageArithmetic2D(obj0, obj1, obj2, op,
		  				 overwrite);
		  break;
		case WLZ_3D_DOMAINOBJ:
		  errNum = WlzImageArithmetic3D(obj0, obj1, obj2, op,
		  				overwrite);
		  break;
	      }
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    dstObj = obj2;
	    obj2 = NULL;
	  }
	  else if(obj2)
	  {
	    (void )WlzFreeObj(obj2);
	  }
	}
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
        break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}
