#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzError.c
* \author       Bill Hill
* \date         November March 1999
* \version      $Id$
* \note
*               Copyright
*               2002 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Error related functions for the Woolz library.
* \ingroup	WlzError
* \todo         -
* \bug          None known.
*/
#pragma ident "MRC HGU $Id$"
#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

/*!
* \return	Woolz error code.
* \ingroup	WlzError
* \brief	Converts an Alg error code to a Woolz error code.
* \param	algErr			Given Alg error code.
*/
WlzErrorNum	WlzErrorFromAlg(AlgError algErr)
{
  WlzErrorNum	wlzErr = WLZ_ERR_ALG;

  switch(algErr)
  {
    case ALG_ERR_NONE:
      wlzErr = WLZ_ERR_NONE;
      break;
    case ALG_ERR_FUNC:
      wlzErr = WLZ_ERR_PARAM_DATA;
      break;
    case ALG_ERR_MALLOC:
      wlzErr = WLZ_ERR_MEM_ALLOC;
      break;
    case ALG_ERR_SINGULAR:
      wlzErr = WLZ_ERR_ALG_SINGULAR;
      break;
    case ALG_ERR_HOMOGENEOUS:
      wlzErr = WLZ_ERR_ALG_HOMOGENEOUS;
      break;
    case ALG_ERR_CONVERGENCE:
      wlzErr = WLZ_ERR_ALG_CONVERGENCE;
      break;
    case ALG_ERR_DIVZERO:
      wlzErr = WLZ_ERR_ALG;
      break;
  }
  return(wlzErr);
}
