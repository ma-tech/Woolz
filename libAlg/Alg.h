#ifndef ALG_H
#define ALG_H
#pragma ident "MRC HGU $Id$"
/*!
* \file         Alg.h
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \defgroup     Alg
* \brief        MRC HGU numerical algorithms library.
* \todo         -
* \bug          None known.
*/

/*!
* \mainpage     
* This library contains numerical algorithms which depend on no other
* non system library apart from the allocation and fundamental types
* library libAlc.
* 
*/

#ifndef __EXTENSIONS__
#define __EXTENSIONS__
#endif

#ifdef __cplusplus
using namespace std;
#else
#include <stdlib.h>
#endif

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <Alc.h>
#include <AlgType.h>
#include <AlgProto.h>

#ifdef __cplusplus
}
#endif

#endif /* ! ALG_H */
