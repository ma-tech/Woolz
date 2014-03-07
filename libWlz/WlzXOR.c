#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzXOR_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzXOR.c
* \author       Bill Hill
* \date         March 2014
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2014],
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
* \brief	Functions for computing the set exclusive or of
* 		objects.
* \ingroup	WlzDomainOps
*/
#include <Wlz.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/*!
* \return	Object with domain equal to the set exclusive or of
* 		the two given objects.
* \ingroup	WlzDomainOps
* \brief	Calculates the domain exclusive or of the two given
* 		objects.
* 		The exclusive or object \f$O_x\f$ is computed using
* 		\f[
		O_x = (O_0 - O_1) \cup (O_1 - O_0)
 		\f]
*		where \f$-\f$ and \f$\cup\f$ are the set difference and
*		union operators respectively. This seems to be slightly
*		quicker than the alternative:
*		\f[
		O_x = (O_0 \cup O_1) - (O_1 \cap O_0)
		\f]
* \param	o0			First object.
* \param	o1			Second object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzXORDom(WlzObject *o0, WlzObject *o1, WlzErrorNum *dstErr)
{
  WlzObject	*d0 = NULL,
  		*d1 = NULL,
		*x0 = NULL;
  WlzErrorNum	errNum0 = WLZ_ERR_NONE,
  		errNum1 = WLZ_ERR_NONE;

#ifdef _OPENMP
#pragma omp parallel shared(o0,o1)
  {
#pragma omp sections
    {
#pragma omp section
      {
#endif
	d0 = WlzDiffDomain(o0, o1, &errNum0);
#ifdef _OPENMP
      }
#pragma omp section
      {
#endif
	d1 = WlzDiffDomain(o1, o0, &errNum1);
#ifdef _OPENMP
      }
    }
  }
#endif
  if((errNum0 == WLZ_ERR_NONE) && (errNum1 != WLZ_ERR_NONE))
  {
    errNum0 = errNum1;
  }
  if(errNum0 == WLZ_ERR_NONE)
  {
    x0 = WlzUnion2(d0, d1, &errNum0);
  }

  (void )WlzFreeObj(d0);
  (void )WlzFreeObj(d1);
  if(dstErr != NULL)
  {
    *dstErr = errNum0;
  }
  return(x0);
}
