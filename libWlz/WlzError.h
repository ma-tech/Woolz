#ifndef WLZ_ERROR_H
#define WLZ_ERROR_H
#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzError.h
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Definitions of Woolz error codes.
* \ingroup	WlzError
* \todo         -
* \bug          None known.
*/

#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */

/*!
* \enum		_WlzErrorNum
* \ingroup	WlzError
* \brief	Woolz error codes (or error numbers) have an integer value
*		with WLZ_ERR_NONE, the no error code, defined to be zero.
*		Typedef: WlzErrorNum.
*/
typedef enum _WlzErrorNum
{
  WLZ_ERR_NONE		= 0,	/*!< No error, defined to be zero. */
  WLZ_ERR_EOO,			/*!< End of objects, not necessarily an
  				     error. */
  WLZ_ERR_OBJECT_NULL,		/*!< Object pointer NULL. */
  WLZ_ERR_OBJECT_TYPE,		/*!< Object type is invalid or
  				     inappropriate. */
  WLZ_ERR_OBJECT_DATA,		/*!< Object data is invalid or
  				     inappropriate. */
  WLZ_ERR_DOMAIN_NULL,		/*!< Domain pointer NULL. */
  WLZ_ERR_DOMAIN_TYPE,		/*!< Domain type is invalid or
  				     inappropriate. */
  WLZ_ERR_DOMAIN_DATA,		/*!< Domain data is invalid or
  				     inappropriate. */
  WLZ_ERR_VALUES_NULL,		/*!< Values pointer NULL. */
  WLZ_ERR_VALUES_TYPE,		/*!< Values type is invalid or
  				     inappropriate. */
  WLZ_ERR_VALUES_DATA,		/*!< Values data is invalid or
  				     inappropriate. */
  WLZ_ERR_PROPERTY_NULL,	/*!< Property pointer NULL. */
  WLZ_ERR_PROPERTY_TYPE,	/*!< Property type is invalid or
  				     inappropriate. */
  WLZ_ERR_GMELM_NULL,		/*!< Geometric model element pointer is
  				     NULL. */
  WLZ_ERR_GMELM_TYPE,		/*!< Geometric model element type is invalid
  				     or inappropriate. */
  WLZ_ERR_GMELM_DATA,		/*!< Geometric model element data is invalid
  				     or inappropriate. */
  WLZ_ERR_PARAM_NULL,		/*!< A library function parameter is NULL. */
  WLZ_ERR_PARAM_TYPE,		/*!< A library function parameter's type is
  				     invalid or inappropriate. */
  WLZ_ERR_PARAM_DATA,		/*!< A library function parameter's data is
  				     invalid or inappropriate. */
  WLZ_ERR_INT_DATA,		/*!< Data is not int. */
  WLZ_ERR_SHORT_DATA,		/*!< Data is not short. */
  WLZ_ERR_UBYTE_DATA,		/*!< Data is not unsigned byte. */
  WLZ_ERR_FLOAT_DATA,		/*!< Data is not float. */
  WLZ_ERR_DOUBLE_DATA,		/*!< Data is not double. */
  WLZ_ERR_GREY_TYPE,		/*!< Grey value type is invalid or
  				     inappropriate. */
  WLZ_ERR_GREY_DATA,		/*!< Grey value data is invalid or
  				     inappropriate. */
  WLZ_ERR_PLANEDOMAIN_TYPE,	/*!< Plane domain type is invalid or
  				     inappropriate. */
  WLZ_ERR_PLANEDOMAIN_DATA,	/*!< Plane domain data is invalid or
  				     inappropriate. */
  WLZ_ERR_INTERVALDOMAIN_NULL,	/*!< Interval domain pointer is NULL. */
  WLZ_ERR_INTERVALDOMAIN_TYPE,	/*!< Interval domain type is invalid or
  				     inappropriate. */
  WLZ_ERR_INTERVALLINE_NULL,	/*!< Interval line pointer is NULL. */
  WLZ_ERR_INTERVAL_NULL,	/*!< Interval pointer is NULL. */
  WLZ_ERR_INTERVAL_DATA,	/*!< Interval data is invalid. */
  WLZ_ERR_INTERVAL_ADJACENT,	/*!< Interval is adjacent to another
  				     interval.  */
  WLZ_ERR_INTERVAL_BOUND,	/*!< Interval bounds are invalid. */
  WLZ_ERR_INTERVAL_NUMBER,	/*!< Number of intervals is incorrect. */
  WLZ_ERR_TRANSFORM_DATA,	/*!< Transform data is invalid or
  				     inappropriate. */
  WLZ_ERR_TRANSFORM_TYPE,	/*!< Transform type is invalid or
  				     inappropriate. */
  WLZ_ERR_VOXELVALUES_TYPE,	/*!< Voxel values type is invalid or
  				     inappropriate. */
  WLZ_ERR_COLUMN_DATA,		/*!< Domain column bounds are invalid. */
  WLZ_ERR_LINE_DATA,		/*!< Domain line bounds are invalid. */
  WLZ_ERR_PLANE_DATA,		/*!< Domain plane bounds are invalid. */
  WLZ_ERR_BINARY_OPERATOR_TYPE,	/*!< Binary operator type is invalid. */
  WLZ_ERR_COMPTHRESH_TYPE,	/*!< Automatic threshold computation method
  				     is invalid. */
  WLZ_ERR_CONNECTIVITY_TYPE,	/*!< Connectivity type is invalid or
  				     inappropriate. */
  WLZ_ERR_INTERPOLATION_TYPE,	/*!< Interpolation method is invalid or
  				     inappropriate. */
  WLZ_ERR_POLYGON_TYPE,		/*!< Polygon type is invalid or
  				     inappropriate. */
  WLZ_ERR_RASTERDIR_TYPE,	/*!< Raster direction is invalid or
  				     inappropriate. */
  WLZ_ERR_LINKCOUNT_DATA,	/*!< Link count value is invalid because an
  				     object, domain or value table or some
				     other element has been freed. */
  WLZ_ERR_MEM_ALLOC,		/*!< Memory allocation failure. */
  WLZ_ERR_MEM_FREE,		/*!< Memory freeing failure. */
  WLZ_ERR_READ_EOF,		/*!< End of file encountered when reading an
  				     object, not necessarily an error. */
  WLZ_ERR_READ_INCOMPLETE,	/*!< End of file or some other error
  				     during an object being read. */
  WLZ_ERR_WRITE_EOF,		/*!< End of file encountered when writing to
  				     a file. */
  WLZ_ERR_WRITE_INCOMPLETE,	/*!< End of file or some other error during
  				     an object being written. */
  WLZ_ERR_ALG,			/*!< General algorithm failure. */
  WLZ_ERR_ALG_SINGULAR,		/*!< Algorithm failure due to a singular
  				     matrix or value. */
  WLZ_ERR_ALG_HOMOGENEOUS,	/*!< Algorithm failure due to a homogeneous
  				     matrix. */
  WLZ_ERR_ALG_CONVERGENCE,	/*!< Algorithm convergence failure. */
  WLZ_ERR_UNIMPLEMENTED,	/*!< A function has not been implemented. */
  WLZ_ERR_UNSPECIFIED,		/*!< All other errors. */
  /* Keep WLZ_ERR_COUNT the last enumerator! */
  WLZ_ERR_COUNT			/*!< Not an error but the number of errors! */
} WlzErrorNum;

#ifdef  __cplusplus
}
#endif /* __cplusplus */

#endif	/* !WLZ_ERROR_H Don't put anything after this line */
