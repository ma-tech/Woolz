#ifndef WLZ_ERROR_H
#define WLZ_ERROR_H
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzError.h
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Defines the Woolz error numbers.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 28-01-00 bill Add ALG_ERR_CONVERGENCE.
************************************************************************/

#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */

/************************************************************************
* Woolz error numbers: The error number for no error (WLZ_ERR_NONE)	*
* has an integer value of zero, but there are no integer values set	*
* for the other error numbers.						*
************************************************************************/
typedef enum
{
  WLZ_ERR_NONE		= 0,
  WLZ_ERR_EOO,
  WLZ_ERR_OBJECT_NULL,
  WLZ_ERR_OBJECT_TYPE,
  WLZ_ERR_OBJECT_DATA,
  WLZ_ERR_DOMAIN_NULL ,
  WLZ_ERR_DOMAIN_TYPE,
  WLZ_ERR_DOMAIN_DATA,
  WLZ_ERR_VALUES_NULL,
  WLZ_ERR_VALUES_TYPE,
  WLZ_ERR_VALUES_DATA,
  WLZ_ERR_PROPERTY_NULL,
  WLZ_ERR_PROPERTY_TYPE,
  WLZ_ERR_PARAM_NULL,
  WLZ_ERR_PARAM_TYPE,
  WLZ_ERR_PARAM_DATA,
  WLZ_ERR_INT_DATA,
  WLZ_ERR_SHORT_DATA,
  WLZ_ERR_UBYTE_DATA,
  WLZ_ERR_FLOAT_DATA,
  WLZ_ERR_DOUBLE_DATA,
  WLZ_ERR_GREY_TYPE,
  WLZ_ERR_GREY_DATA,
  WLZ_ERR_PLANEDOMAIN_TYPE,
  WLZ_ERR_PLANEDOMAIN_DATA,
  WLZ_ERR_INTERVALDOMAIN_NULL,
  WLZ_ERR_INTERVALDOMAIN_TYPE,
  WLZ_ERR_INTERVALLINE_NULL,
  WLZ_ERR_INTERVAL_NULL,
  WLZ_ERR_INTERVAL_DATA,
  WLZ_ERR_INTERVAL_ADJACENT,
  WLZ_ERR_INTERVAL_BOUND,
  WLZ_ERR_INTERVAL_NUMBER,
  WLZ_ERR_TRANSFORM_DATA,
  WLZ_ERR_TRANSFORM_TYPE,
  WLZ_ERR_VOXELVALUES_TYPE,
  WLZ_ERR_COLUMN_DATA,
  WLZ_ERR_LINE_DATA,
  WLZ_ERR_PLANE_DATA,
  WLZ_ERR_BINARY_OPERATOR_TYPE,
  WLZ_ERR_COMPTHRESH_TYPE,
  WLZ_ERR_CONNECTIVITY_TYPE,
  WLZ_ERR_INTERPOLATION_TYPE,
  WLZ_ERR_POINT_TYPE,
  WLZ_ERR_POLYGON_TYPE,
  WLZ_ERR_RASTERDIR_TYPE,
  WLZ_ERR_VECTOR_TYPE,
  WLZ_ERR_LINKCOUNT_DATA,
  WLZ_ERR_MEM_ALLOC,
  WLZ_ERR_MEM_FREE,
  WLZ_ERR_READ_EOF,
  WLZ_ERR_READ_INCOMPLETE,
  WLZ_ERR_WRITE_EOF,
  WLZ_ERR_WRITE_INCOMPLETE,
  WLZ_ERR_ALG,
  WLZ_ERR_ALG_SINGULAR,
  WLZ_ERR_ALG_HOMOGENEOUS,
  WLZ_ERR_ALG_CONVERGENCE,
  WLZ_ERR_UNSPECIFIED,
  /**********************************************************************
  * WLZ_ERR_COUNT is not an error number. It is the number of errors.	*
  * Keep it the last enumerator!					*
  **********************************************************************/
  WLZ_ERR_COUNT
} WlzErrorNum;

#ifdef  __cplusplus
}
#endif /* __cplusplus */

#endif	/* !WLZ_ERROR_H Don't put anything after this line */
