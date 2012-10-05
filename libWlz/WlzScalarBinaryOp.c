#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzScalarBinaryOp_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzScalarBinaryOp.c
* \author       Richard Baldock
* \date         February 2012
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
* \brief	Applies scalar binary operators to objects.
* \ingroup	WlzArithmetic
*/

#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

static WlzErrorNum WlzScalarBinaryOp3d(WlzObject	*o1,
				       WlzPixelV	pval,
				       WlzObject	*o3,
				       WlzBinaryOperatorType op);

static WlzErrorNum WlzBufIntIntScalarBinaryOp(
  int		*inbuf1,
  int		val,
  WlzPixelP	outPP,
  int		buflen,
  WlzBinaryOperatorType	op){

  switch( outPP.type ){
  case WLZ_GREY_INT:
    switch( op ){
    case WLZ_BO_ADD:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ + val;}
      break;
    case WLZ_BO_SUBTRACT:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ - val;}
      break;
    case WLZ_BO_MULTIPLY:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ * val;}
      break;
    case WLZ_BO_DIVIDE:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ / val;}
      break;
    case WLZ_BO_MODULUS:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ % val;}
      break;
    case WLZ_BO_EQ:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ == val;}
      break;
    case WLZ_BO_NE:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ != val;}
      break;
    case WLZ_BO_GT:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ > val;}
      break;
    case WLZ_BO_GE:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ >= val;}
      break;
    case WLZ_BO_LT:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ < val;}
      break;
    case WLZ_BO_LE:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ <= val;}
      break;
    case WLZ_BO_AND:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ & val;}
      break;
    case WLZ_BO_OR:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ | val;}
      break;
    case WLZ_BO_XOR:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ ^ val;}
      break;
    default:
      break;
    }
    break;
  case WLZ_GREY_SHORT:
    switch( op ){
    case WLZ_BO_ADD:
      while( buflen-- ){*outPP.p.shp++ = (short )(*inbuf1++ + val);}
      break;
    case WLZ_BO_SUBTRACT:
      while( buflen-- ){*outPP.p.shp++ = (short )(*inbuf1++ - val);}
      break;
    case WLZ_BO_MULTIPLY:
      while( buflen-- ){*outPP.p.shp++ = (short )(*inbuf1++ * val);}
      break;
    case WLZ_BO_DIVIDE:
      while( buflen-- ){*outPP.p.shp++ = (short )(*inbuf1++ / val);}
      break;
    case WLZ_BO_MODULUS:
      while( buflen-- ){*outPP.p.shp++ = *inbuf1++ % val;}
      break;
    case WLZ_BO_EQ:
      while( buflen-- ){*outPP.p.shp++ = *inbuf1++ == val;}
      break;
    case WLZ_BO_NE:
      while( buflen-- ){*outPP.p.shp++ = *inbuf1++ != val;}
      break;
    case WLZ_BO_GT:
      while( buflen-- ){*outPP.p.shp++ = *inbuf1++ > val;}
      break;
    case WLZ_BO_GE:
      while( buflen-- ){*outPP.p.shp++ = *inbuf1++ >= val;}
      break;
    case WLZ_BO_LT:
      while( buflen-- ){*outPP.p.shp++ = *inbuf1++ < val;}
      break;
    case WLZ_BO_LE:
      while( buflen-- ){*outPP.p.shp++ = *inbuf1++ <= val;}
      break;
    case WLZ_BO_AND:
      while( buflen-- ){*outPP.p.shp++ = *inbuf1++ & val;}
      break;
    case WLZ_BO_OR:
      while( buflen-- ){*outPP.p.shp++ = *inbuf1++ | val;}
      break;
    case WLZ_BO_XOR:
      while( buflen-- ){*outPP.p.shp++ = *inbuf1++ ^ val;}
      break;
    default:
      break;
    }
    break;
  case WLZ_GREY_UBYTE:
    switch( op ){
    case WLZ_BO_ADD:
      while( buflen-- ){*outPP.p.ubp++ = (WlzUByte )(*inbuf1++ + val);}
      break;
    case WLZ_BO_SUBTRACT:
      while( buflen-- ){*outPP.p.ubp++ = (WlzUByte )(*inbuf1++ - val);}
      break;
    case WLZ_BO_MULTIPLY:
      while( buflen-- ){*outPP.p.ubp++ = (WlzUByte )(*inbuf1++ * val);}
      break;
    case WLZ_BO_DIVIDE:
      while( buflen-- ){*outPP.p.ubp++ = (WlzUByte )(*inbuf1++ / val);}
      break;
    case WLZ_BO_MODULUS:
      while( buflen-- ){*outPP.p.ubp++ = *inbuf1++ % val;}
      break;
    case WLZ_BO_EQ:
      while( buflen-- ){*outPP.p.ubp++ = *inbuf1++ == val;}
      break;
    case WLZ_BO_NE:
      while( buflen-- ){*outPP.p.ubp++ = *inbuf1++ != val;}
      break;
    case WLZ_BO_GT:
      while( buflen-- ){*outPP.p.ubp++ = *inbuf1++ > val;}
      break;
    case WLZ_BO_GE:
      while( buflen-- ){*outPP.p.ubp++ = *inbuf1++ >= val;}
      break;
    case WLZ_BO_LT:
      while( buflen-- ){*outPP.p.ubp++ = *inbuf1++ < val;}
      break;
    case WLZ_BO_LE:
      while( buflen-- ){*outPP.p.ubp++ = *inbuf1++ <= val;}
      break;
    case WLZ_BO_AND:
      while( buflen-- ){*outPP.p.ubp++ = *inbuf1++ & val;}
      break;
    case WLZ_BO_OR:
      while( buflen-- ){*outPP.p.ubp++ = *inbuf1++ | val;}
      break;
    case WLZ_BO_XOR:
      while( buflen-- ){*outPP.p.ubp++ = *inbuf1++ ^ val;}
      break;
    default:
      break;
    }
    break;
  case WLZ_GREY_FLOAT:
    switch( op ){
    case WLZ_BO_ADD:
      while( buflen-- ){*outPP.p.flp++ = (float )(*inbuf1++ + val);}
      break;
    case WLZ_BO_SUBTRACT:
      while( buflen-- ){*outPP.p.flp++ = (float )(*inbuf1++ - val);}
      break;
    case WLZ_BO_MULTIPLY:
      while( buflen-- ){*outPP.p.flp++ = (float )(*inbuf1++ * val);}
      break;
    case WLZ_BO_DIVIDE:
      while( buflen-- ){*outPP.p.flp++ = (float )(*inbuf1++ / val);}
      break;
    case WLZ_BO_MODULUS:
      while( buflen-- ){*outPP.p.flp++ = *inbuf1++ % val;}
      break;
    case WLZ_BO_EQ:
      while( buflen-- ){*outPP.p.flp++ = *inbuf1++ == val;}
      break;
    case WLZ_BO_NE:
      while( buflen-- ){*outPP.p.flp++ = *inbuf1++ != val;}
      break;
    case WLZ_BO_GT:
      while( buflen-- ){*outPP.p.flp++ = *inbuf1++ > val;}
      break;
    case WLZ_BO_GE:
      while( buflen-- ){*outPP.p.flp++ = *inbuf1++ >= val;}
      break;
    case WLZ_BO_LT:
      while( buflen-- ){*outPP.p.flp++ = *inbuf1++ < val;}
      break;
    case WLZ_BO_LE:
      while( buflen-- ){*outPP.p.flp++ = *inbuf1++ <= val;}
      break;
    case WLZ_BO_AND:
      while( buflen-- ){*outPP.p.flp++ = *inbuf1++ & val;}
      break;
    case WLZ_BO_OR:
      while( buflen-- ){*outPP.p.flp++ = *inbuf1++ | val;}
      break;
    case WLZ_BO_XOR:
      while( buflen-- ){*outPP.p.flp++ = *inbuf1++ ^ val;}
      break;
    default:
      break;
    }
    break;
  case WLZ_GREY_DOUBLE:
    switch( op ){
    case WLZ_BO_ADD:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ + val;}
      break;
    case WLZ_BO_SUBTRACT:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ - val;}
      break;
    case WLZ_BO_MULTIPLY:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ * val;}
      break;
    case WLZ_BO_DIVIDE:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ / val;}
      break;
    case WLZ_BO_MODULUS:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ % val;}
      break;
    case WLZ_BO_EQ:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ == val;}
      break;
    case WLZ_BO_NE:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ != val;}
      break;
    case WLZ_BO_GT:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ > val;}
      break;
    case WLZ_BO_GE:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ >= val;}
      break;
    case WLZ_BO_LT:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ < val;}
      break;
    case WLZ_BO_LE:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ <= val;}
      break;
    case WLZ_BO_AND:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ & val;}
      break;
    case WLZ_BO_OR:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ | val;}
      break;
    case WLZ_BO_XOR:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ ^ val;}
      break;
    default:
      break;
    }
    break;
  default:
    return WLZ_ERR_GREY_TYPE;
  }

  return WLZ_ERR_NONE;
}

static WlzErrorNum WlzBufIntDblScalarBinaryOp(
  int		*inbuf1,
  double	val,
  WlzPixelP	outPP,
  int		buflen,
  WlzBinaryOperatorType	op){

  switch( outPP.type ){
  case WLZ_GREY_INT:
    switch( op ){
    case WLZ_BO_ADD:
      while( buflen-- ){*outPP.p.inp++ = (int )(*inbuf1++ + val);}
      break;
    case WLZ_BO_SUBTRACT:
      while( buflen-- ){*outPP.p.inp++ = (int )(*inbuf1++ - val);}
      break;
    case WLZ_BO_MULTIPLY:
      while( buflen-- ){*outPP.p.inp++ = (int )(*inbuf1++ * val);}
      break;
    case WLZ_BO_DIVIDE:
      while( buflen-- ){*outPP.p.inp++ = (int )(*inbuf1++ / val);}
      break;
    case WLZ_BO_EQ:
      while( buflen-- ){*outPP.p.inp++ = 
			  fabs(*inbuf1++ - val) < DBL_EPSILON;}
      break;
    case WLZ_BO_NE:
      while( buflen-- ){*outPP.p.inp++ = 
			  fabs(*inbuf1++ - val) > DBL_EPSILON;}
      break;
    case WLZ_BO_GT:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ > val;}
      break;
    case WLZ_BO_GE:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ >= val;}
      break;
    case WLZ_BO_LT:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ < val;}
      break;
    case WLZ_BO_LE:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ <= val;}
      break;
    case WLZ_BO_MODULUS: /* FALLTHROUGH */
    case WLZ_BO_AND: /* FALLTHROUGH */
    case WLZ_BO_OR: /* FALLTHROUGH */
    case WLZ_BO_XOR: /* FALLTHROUGH */
    default:
      return WLZ_ERR_BINARY_OPERATOR_TYPE;
    }
    break;
  case WLZ_GREY_SHORT:
    switch( op ){
    case WLZ_BO_ADD:
      while( buflen-- ){*outPP.p.shp++ = (short )(*inbuf1++ + val);}
      break;
    case WLZ_BO_SUBTRACT:
      while( buflen-- ){*outPP.p.shp++ = (short )(*inbuf1++ - val);}
      break;
    case WLZ_BO_MULTIPLY:
      while( buflen-- ){*outPP.p.shp++ = (short )(*inbuf1++ * val);}
      break;
    case WLZ_BO_DIVIDE:
      while( buflen-- ){*outPP.p.shp++ = (short )(*inbuf1++ / val);}
      break;
    case WLZ_BO_EQ:
      while( buflen-- ){*outPP.p.shp++ = 
			  fabs(*inbuf1++ - val) < DBL_EPSILON;}
      break;
    case WLZ_BO_NE:
      while( buflen-- ){*outPP.p.shp++ = 
			  fabs(*inbuf1++ - val) > DBL_EPSILON;}
      break;
    case WLZ_BO_GT:
      while( buflen-- ){*outPP.p.shp++ = *inbuf1++ > val;}
      break;
    case WLZ_BO_GE:
      while( buflen-- ){*outPP.p.shp++ = *inbuf1++ >= val;}
      break;
    case WLZ_BO_LT:
      while( buflen-- ){*outPP.p.shp++ = *inbuf1++ < val;}
      break;
    case WLZ_BO_LE:
      while( buflen-- ){*outPP.p.shp++ = *inbuf1++ <= val;}
      break;
    case WLZ_BO_MODULUS: /* FALLTHROUGH */
    case WLZ_BO_AND: /* FALLTHROUGH */
    case WLZ_BO_OR: /* FALLTHROUGH */
    case WLZ_BO_XOR: /* FALLTHROUGH */
    default:
      return WLZ_ERR_BINARY_OPERATOR_TYPE;
    }
    break;
  case WLZ_GREY_UBYTE:
    switch( op ){
    case WLZ_BO_ADD:
      while( buflen-- ){*outPP.p.ubp++ = (WlzUByte )(*inbuf1++ + val);}
      break;
    case WLZ_BO_SUBTRACT:
      while( buflen-- ){*outPP.p.ubp++ = (WlzUByte )(*inbuf1++ - val);}
      break;
    case WLZ_BO_MULTIPLY:
      while( buflen-- ){*outPP.p.ubp++ = (WlzUByte )(*inbuf1++ * val);}
      break;
    case WLZ_BO_DIVIDE:
      while( buflen-- ){*outPP.p.ubp++ = (WlzUByte )(*inbuf1++ / val);}
      break;
    case WLZ_BO_EQ:
      while( buflen-- ){*outPP.p.ubp++ = 
			  fabs(*inbuf1++ - val) < DBL_EPSILON;}
      break;
    case WLZ_BO_NE:
      while( buflen-- ){*outPP.p.ubp++ = 
			  fabs(*inbuf1++ - val) > DBL_EPSILON;}
      break;
    case WLZ_BO_GT:
      while( buflen-- ){*outPP.p.ubp++ = *inbuf1++ > val;}
      break;
    case WLZ_BO_GE:
      while( buflen-- ){*outPP.p.ubp++ = *inbuf1++ >= val;}
      break;
    case WLZ_BO_LT:
      while( buflen-- ){*outPP.p.ubp++ = *inbuf1++ < val;}
      break;
    case WLZ_BO_LE:
      while( buflen-- ){*outPP.p.ubp++ = *inbuf1++ <= val;}
      break;
    case WLZ_BO_MODULUS: /* FALLTHROUGH */
    case WLZ_BO_AND: /* FALLTHROUGH */
    case WLZ_BO_OR: /* FALLTHROUGH */
    case WLZ_BO_XOR: /* FALLTHROUGH */
    default:
      return WLZ_ERR_BINARY_OPERATOR_TYPE;
    }
    break;
  case WLZ_GREY_FLOAT:
    switch( op ){
    case WLZ_BO_ADD:
      while( buflen-- ){*outPP.p.flp++ = (float )(*inbuf1++ + val);}
      break;
    case WLZ_BO_SUBTRACT:
      while( buflen-- ){*outPP.p.flp++ = (float )(*inbuf1++ - val);}
      break;
    case WLZ_BO_MULTIPLY:
      while( buflen-- ){*outPP.p.flp++ = (float )(*inbuf1++ * val);}
      break;
    case WLZ_BO_DIVIDE:
      while( buflen-- ){*outPP.p.flp++ = (float )(*inbuf1++ / val);}
      break;
    case WLZ_BO_EQ:
      while( buflen-- ){*outPP.p.flp++ = 
			  fabs(*inbuf1++ - val) < DBL_EPSILON;}
      break;
    case WLZ_BO_NE:
      while( buflen-- ){*outPP.p.flp++ = 
			  fabs(*inbuf1++ - val) > DBL_EPSILON;}
      break;
    case WLZ_BO_GT:
      while( buflen-- ){*outPP.p.flp++ = *inbuf1++ > val;}
      break;
    case WLZ_BO_GE:
      while( buflen-- ){*outPP.p.flp++ = *inbuf1++ >= val;}
      break;
    case WLZ_BO_LT:
      while( buflen-- ){*outPP.p.flp++ = *inbuf1++ < val;}
      break;
    case WLZ_BO_LE:
      while( buflen-- ){*outPP.p.flp++ = *inbuf1++ <= val;}
      break;
    case WLZ_BO_MODULUS: /* FALLTHROUGH */
    case WLZ_BO_AND: /* FALLTHROUGH */
    case WLZ_BO_OR: /* FALLTHROUGH */
    case WLZ_BO_XOR: /* FALLTHROUGH */
    default:
      return WLZ_ERR_BINARY_OPERATOR_TYPE;
    }
    break;
  case WLZ_GREY_DOUBLE:
    switch( op ){
    case WLZ_BO_ADD:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ + val;}
      break;
    case WLZ_BO_SUBTRACT:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ - val;}
      break;
    case WLZ_BO_MULTIPLY:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ * val;}
      break;
    case WLZ_BO_DIVIDE:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ / val;}
      break;
    case WLZ_BO_EQ:
      while( buflen-- ){*outPP.p.dbp++ = 
			  fabs(*inbuf1++ - val) < DBL_EPSILON;}
      break;
    case WLZ_BO_NE:
      while( buflen-- ){*outPP.p.dbp++ = 
			  fabs(*inbuf1++ - val) > DBL_EPSILON;}
      break;
    case WLZ_BO_GT:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ > val;}
      break;
    case WLZ_BO_GE:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ >= val;}
      break;
    case WLZ_BO_LT:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ < val;}
      break;
    case WLZ_BO_LE:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ <= val;}
      break;
    case WLZ_BO_MODULUS: /* FALLTHROUGH */
    case WLZ_BO_AND: /* FALLTHROUGH */
    case WLZ_BO_OR: /* FALLTHROUGH */
    case WLZ_BO_XOR: /* FALLTHROUGH */
    default:
      return WLZ_ERR_BINARY_OPERATOR_TYPE;
    }
    break;
  default:
    return WLZ_ERR_GREY_TYPE;
  }

  return WLZ_ERR_NONE;
}

static WlzErrorNum WlzBufDblIntScalarBinaryOp(
  double	*inbuf1,
  int		val,
  WlzPixelP	outPP,
  int		buflen,
  WlzBinaryOperatorType	op){

  switch( outPP.type ){
  case WLZ_GREY_INT:
    switch( op ){
    case WLZ_BO_ADD:
      while( buflen-- ){*outPP.p.inp++ = (int )(*inbuf1++ + val);}
      break;
    case WLZ_BO_SUBTRACT:
      while( buflen-- ){*outPP.p.inp++ = (int )(*inbuf1++ - val);}
      break;
    case WLZ_BO_MULTIPLY:
      while( buflen-- ){*outPP.p.inp++ = (int )(*inbuf1++ * val);}
      break;
    case WLZ_BO_DIVIDE:
      while( buflen-- ){*outPP.p.inp++ = (int )(*inbuf1++ / val);}
      break;
    case WLZ_BO_EQ:
      while( buflen-- ){*outPP.p.inp++ = 
			  fabs(*inbuf1++ - val) < DBL_EPSILON;}
      break;
    case WLZ_BO_NE:
      while( buflen-- ){*outPP.p.inp++ = 
			  fabs(*inbuf1++ - val) > DBL_EPSILON;}
      break;
    case WLZ_BO_GT:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ > val;}
      break;
    case WLZ_BO_GE:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ >= val;}
      break;
    case WLZ_BO_LT:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ < val;}
      break;
    case WLZ_BO_LE:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ <= val;}
      break;
    case WLZ_BO_MODULUS: /* FALLTHROUGH */
    case WLZ_BO_AND: /* FALLTHROUGH */
    case WLZ_BO_OR: /* FALLTHROUGH */
    case WLZ_BO_XOR: /* FALLTHROUGH */
    default:
      return WLZ_ERR_BINARY_OPERATOR_TYPE;
    }
    break;
  case WLZ_GREY_SHORT:
    switch( op ){
    case WLZ_BO_ADD:
      while( buflen-- ){*outPP.p.shp++ = (short )(*inbuf1++ + val);}
      break;
    case WLZ_BO_SUBTRACT:
      while( buflen-- ){*outPP.p.shp++ = (short )(*inbuf1++ - val);}
      break;
    case WLZ_BO_MULTIPLY:
      while( buflen-- ){*outPP.p.shp++ = (short )(*inbuf1++ * val);}
      break;
    case WLZ_BO_DIVIDE:
      while( buflen-- ){*outPP.p.shp++ = (short )(*inbuf1++ / val);}
      break;
    case WLZ_BO_EQ:
      while( buflen-- ){*outPP.p.shp++ = 
			  fabs(*inbuf1++ - val) < DBL_EPSILON;}
      break;
    case WLZ_BO_NE:
      while( buflen-- ){*outPP.p.shp++ = 
			  fabs(*inbuf1++ - val) > DBL_EPSILON;}
      break;
    case WLZ_BO_GT:
      while( buflen-- ){*outPP.p.shp++ = *inbuf1++ > val;}
      break;
    case WLZ_BO_GE:
      while( buflen-- ){*outPP.p.shp++ = *inbuf1++ >= val;}
      break;
    case WLZ_BO_LT:
      while( buflen-- ){*outPP.p.shp++ = *inbuf1++ < val;}
      break;
    case WLZ_BO_LE:
      while( buflen-- ){*outPP.p.shp++ = *inbuf1++ <= val;}
      break;
    case WLZ_BO_MODULUS: /* FALLTHROUGH */
    case WLZ_BO_AND: /* FALLTHROUGH */
    case WLZ_BO_OR: /* FALLTHROUGH */
    case WLZ_BO_XOR: /* FALLTHROUGH */
    default:
      return WLZ_ERR_BINARY_OPERATOR_TYPE;
    }
    break;
  case WLZ_GREY_UBYTE:
    switch( op ){
    case WLZ_BO_ADD:
      while( buflen-- ){*outPP.p.ubp++ = (WlzUByte )(*inbuf1++ + val);}
      break;
    case WLZ_BO_SUBTRACT:
      while( buflen-- ){*outPP.p.ubp++ = (WlzUByte )(*inbuf1++ - val);}
      break;
    case WLZ_BO_MULTIPLY:
      while( buflen-- ){*outPP.p.ubp++ = (WlzUByte )(*inbuf1++ * val);}
      break;
    case WLZ_BO_DIVIDE:
      while( buflen-- ){*outPP.p.ubp++ = (WlzUByte )(*inbuf1++ / val);}
      break;
    case WLZ_BO_EQ:
      while( buflen-- ){*outPP.p.ubp++ = 
			  fabs(*inbuf1++ - val) < DBL_EPSILON;}
      break;
    case WLZ_BO_NE:
      while( buflen-- ){*outPP.p.ubp++ = 
			  fabs(*inbuf1++ - val) > DBL_EPSILON;}
      break;
    case WLZ_BO_GT:
      while( buflen-- ){*outPP.p.ubp++ = *inbuf1++ > val;}
      break;
    case WLZ_BO_GE:
      while( buflen-- ){*outPP.p.ubp++ = *inbuf1++ >= val;}
      break;
    case WLZ_BO_LT:
      while( buflen-- ){*outPP.p.ubp++ = *inbuf1++ < val;}
      break;
    case WLZ_BO_LE:
      while( buflen-- ){*outPP.p.ubp++ = *inbuf1++ <= val;}
      break;
    case WLZ_BO_MODULUS: /* FALLTHROUGH */
    case WLZ_BO_AND: /* FALLTHROUGH */
    case WLZ_BO_OR: /* FALLTHROUGH */
    case WLZ_BO_XOR: /* FALLTHROUGH */
    default:
      return WLZ_ERR_BINARY_OPERATOR_TYPE;
    }
    break;
  case WLZ_GREY_FLOAT:
    switch( op ){
    case WLZ_BO_ADD:
      while( buflen-- ){*outPP.p.flp++ = (float )(*inbuf1++ + val);}
      break;
    case WLZ_BO_SUBTRACT:
      while( buflen-- ){*outPP.p.flp++ = (float )(*inbuf1++ - val);}
      break;
    case WLZ_BO_MULTIPLY:
      while( buflen-- ){*outPP.p.flp++ = (float )(*inbuf1++ * val);}
      break;
    case WLZ_BO_DIVIDE:
      while( buflen-- ){*outPP.p.flp++ = (float )(*inbuf1++ / val);}
      break;
    case WLZ_BO_EQ:
      while( buflen-- ){*outPP.p.flp++ = 
			  fabs(*inbuf1++ - val) < DBL_EPSILON;}
      break;
    case WLZ_BO_NE:
      while( buflen-- ){*outPP.p.flp++ = 
			  fabs(*inbuf1++ - val) > DBL_EPSILON;}
      break;
    case WLZ_BO_GT:
      while( buflen-- ){*outPP.p.flp++ = *inbuf1++ > val;}
      break;
    case WLZ_BO_GE:
      while( buflen-- ){*outPP.p.flp++ = *inbuf1++ >= val;}
      break;
    case WLZ_BO_LT:
      while( buflen-- ){*outPP.p.flp++ = *inbuf1++ < val;}
      break;
    case WLZ_BO_LE:
      while( buflen-- ){*outPP.p.flp++ = *inbuf1++ <= val;}
      break;
    case WLZ_BO_MODULUS: /* FALLTHROUGH */
    case WLZ_BO_AND: /* FALLTHROUGH */
    case WLZ_BO_OR: /* FALLTHROUGH */
    case WLZ_BO_XOR: /* FALLTHROUGH */
    default:
      return WLZ_ERR_BINARY_OPERATOR_TYPE;
    }
    break;
  case WLZ_GREY_DOUBLE:
    switch( op ){
    case WLZ_BO_ADD:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ + val;}
      break;
    case WLZ_BO_SUBTRACT:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ - val;}
      break;
    case WLZ_BO_MULTIPLY:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ * val;}
      break;
    case WLZ_BO_DIVIDE:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ / val;}
      break;
    case WLZ_BO_EQ:
      while( buflen-- ){*outPP.p.dbp++ = 
			  fabs(*inbuf1++ - val) < DBL_EPSILON;}
      break;
    case WLZ_BO_NE:
      while( buflen-- ){*outPP.p.dbp++ = 
			  fabs(*inbuf1++ - val) > DBL_EPSILON;}
      break;
    case WLZ_BO_GT:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ > val;}
      break;
    case WLZ_BO_GE:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ >= val;}
      break;
    case WLZ_BO_LT:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ < val;}
      break;
    case WLZ_BO_LE:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ <= val;}
      break;
    case WLZ_BO_MODULUS: /* FALLTHROUGH */
    case WLZ_BO_AND: /* FALLTHROUGH */
    case WLZ_BO_OR: /* FALLTHROUGH */
    case WLZ_BO_XOR: /* FALLTHROUGH */
    default:
      return WLZ_ERR_BINARY_OPERATOR_TYPE;
    }
    break;
  default:
    return WLZ_ERR_GREY_TYPE;
  }

  return WLZ_ERR_NONE;
}

static WlzErrorNum WlzBufDblDblScalarBinaryOp(
  double	*inbuf1,
  double	val,
  WlzPixelP	outPP,
  int		buflen,
  WlzBinaryOperatorType	op){

  switch( outPP.type ){
  case WLZ_GREY_INT:
    switch( op ){
    case WLZ_BO_ADD:
      while( buflen-- ){*outPP.p.inp++ = (int )(*inbuf1++ + val);}
      break;
    case WLZ_BO_SUBTRACT:
      while( buflen-- ){*outPP.p.inp++ = (int )(*inbuf1++ - val);}
      break;
    case WLZ_BO_MULTIPLY:
      while( buflen-- ){*outPP.p.inp++ = (int )(*inbuf1++ * val);}
      break;
    case WLZ_BO_DIVIDE:
      while( buflen-- ){*outPP.p.inp++ = (int )(*inbuf1++ / val);}
      break;
    case WLZ_BO_EQ:
      while( buflen-- ){*outPP.p.inp++ = 
			  fabs(*inbuf1++ - val) < DBL_EPSILON;}
      break;
    case WLZ_BO_NE:
      while( buflen-- ){*outPP.p.inp++ = 
			  fabs(*inbuf1++ - val) > DBL_EPSILON;}
      break;
    case WLZ_BO_GT:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ > val;}
      break;
    case WLZ_BO_GE:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ >= val;}
      break;
    case WLZ_BO_LT:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ < val;}
      break;
    case WLZ_BO_LE:
      while( buflen-- ){*outPP.p.inp++ = *inbuf1++ <= val;}
      break;
    case WLZ_BO_MODULUS: /* FALLTHROUGH */
    case WLZ_BO_AND: /* FALLTHROUGH */
    case WLZ_BO_OR: /* FALLTHROUGH */
    case WLZ_BO_XOR: /* FALLTHROUGH */
    default:
      return WLZ_ERR_BINARY_OPERATOR_TYPE;
    }
    break;
  case WLZ_GREY_SHORT:
    switch( op ){
    case WLZ_BO_ADD:
      while( buflen-- ){*outPP.p.shp++ = (short )(*inbuf1++ + val);}
      break;
    case WLZ_BO_SUBTRACT:
      while( buflen-- ){*outPP.p.shp++ = (short )(*inbuf1++ - val);}
      break;
    case WLZ_BO_MULTIPLY:
      while( buflen-- ){*outPP.p.shp++ = (short )(*inbuf1++ * val);}
      break;
    case WLZ_BO_DIVIDE:
      while( buflen-- ){*outPP.p.shp++ = (short )(*inbuf1++ / val);}
      break;
    case WLZ_BO_EQ:
      while( buflen-- ){*outPP.p.shp++ = 
			  fabs(*inbuf1++ - val) < DBL_EPSILON;}
      break;
    case WLZ_BO_NE:
      while( buflen-- ){*outPP.p.shp++ = 
			  fabs(*inbuf1++ - val) > DBL_EPSILON;}
      break;
    case WLZ_BO_GT:
      while( buflen-- ){*outPP.p.shp++ = *inbuf1++ > val;}
      break;
    case WLZ_BO_GE:
      while( buflen-- ){*outPP.p.shp++ = *inbuf1++ >= val;}
      break;
    case WLZ_BO_LT:
      while( buflen-- ){*outPP.p.shp++ = *inbuf1++ < val;}
      break;
    case WLZ_BO_LE:
      while( buflen-- ){*outPP.p.shp++ = *inbuf1++ <= val;}
      break;
    case WLZ_BO_MODULUS: /* FALLTHROUGH */
    case WLZ_BO_AND: /* FALLTHROUGH */
    case WLZ_BO_OR: /* FALLTHROUGH */
    case WLZ_BO_XOR: /* FALLTHROUGH */
    default:
      return WLZ_ERR_BINARY_OPERATOR_TYPE;
    }
    break;
  case WLZ_GREY_UBYTE:
    switch( op ){
    case WLZ_BO_ADD:
      while( buflen-- ){*outPP.p.ubp++ = (WlzUByte )(*inbuf1++ + val);}
      break;
    case WLZ_BO_SUBTRACT:
      while( buflen-- ){*outPP.p.ubp++ = (WlzUByte )(*inbuf1++ - val);}
      break;
    case WLZ_BO_MULTIPLY:
      while( buflen-- ){*outPP.p.ubp++ = (WlzUByte )(*inbuf1++ * val);}
      break;
    case WLZ_BO_DIVIDE:
      while( buflen-- ){*outPP.p.ubp++ = (WlzUByte )(*inbuf1++ / val);}
      break;
    case WLZ_BO_EQ:
      while( buflen-- ){*outPP.p.ubp++ = 
			  fabs(*inbuf1++ - val) < DBL_EPSILON;}
      break;
    case WLZ_BO_NE:
      while( buflen-- ){*outPP.p.ubp++ = 
			  fabs(*inbuf1++ - val) > DBL_EPSILON;}
      break;
    case WLZ_BO_GT:
      while( buflen-- ){*outPP.p.ubp++ = *inbuf1++ > val;}
      break;
    case WLZ_BO_GE:
      while( buflen-- ){*outPP.p.ubp++ = *inbuf1++ >= val;}
      break;
    case WLZ_BO_LT:
      while( buflen-- ){*outPP.p.ubp++ = *inbuf1++ < val;}
      break;
    case WLZ_BO_LE:
      while( buflen-- ){*outPP.p.ubp++ = *inbuf1++ <= val;}
      break;
    case WLZ_BO_MODULUS: /* FALLTHROUGH */
    case WLZ_BO_AND: /* FALLTHROUGH */
    case WLZ_BO_OR: /* FALLTHROUGH */
    case WLZ_BO_XOR: /* FALLTHROUGH */
    default:
      return WLZ_ERR_BINARY_OPERATOR_TYPE;
    }
    break;
  case WLZ_GREY_FLOAT:
    switch( op ){
    case WLZ_BO_ADD:
      while( buflen-- ){*outPP.p.flp++ = (float )(*inbuf1++ + val);}
      break;
    case WLZ_BO_SUBTRACT:
      while( buflen-- ){*outPP.p.flp++ = (float )(*inbuf1++ - val);}
      break;
    case WLZ_BO_MULTIPLY:
      while( buflen-- ){*outPP.p.flp++ = (float )(*inbuf1++ * val);}
      break;
    case WLZ_BO_DIVIDE:
      while( buflen-- ){*outPP.p.flp++ = (float )(*inbuf1++ / val);}
      break;
    case WLZ_BO_EQ:
      while( buflen-- ){*outPP.p.flp++ = 
			  fabs(*inbuf1++ - val) < DBL_EPSILON;}
      break;
    case WLZ_BO_NE:
      while( buflen-- ){*outPP.p.flp++ = 
			  fabs(*inbuf1++ - val) > DBL_EPSILON;}
      break;
    case WLZ_BO_GT:
      while( buflen-- ){*outPP.p.flp++ = *inbuf1++ > val;}
      break;
    case WLZ_BO_GE:
      while( buflen-- ){*outPP.p.flp++ = *inbuf1++ >= val;}
      break;
    case WLZ_BO_LT:
      while( buflen-- ){*outPP.p.flp++ = *inbuf1++ < val;}
      break;
    case WLZ_BO_LE:
      while( buflen-- ){*outPP.p.flp++ = *inbuf1++ <= val;}
      break;
    case WLZ_BO_MODULUS: /* FALLTHROUGH */
    case WLZ_BO_AND: /* FALLTHROUGH */
    case WLZ_BO_OR: /* FALLTHROUGH */
    case WLZ_BO_XOR: /* FALLTHROUGH */
    default:
      return WLZ_ERR_BINARY_OPERATOR_TYPE;
    }
    break;
  case WLZ_GREY_DOUBLE:
    switch( op ){
    case WLZ_BO_ADD:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ + val;}
      break;
    case WLZ_BO_SUBTRACT:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ - val;}
      break;
    case WLZ_BO_MULTIPLY:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ * val;}
      break;
    case WLZ_BO_DIVIDE:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ / val;}
      break;
    case WLZ_BO_EQ:
      while( buflen-- ){*outPP.p.dbp++ = 
			  fabs(*inbuf1++ - val) < DBL_EPSILON;}
      break;
    case WLZ_BO_NE:
      while( buflen-- ){*outPP.p.dbp++ = 
			  fabs(*inbuf1++ - val) > DBL_EPSILON;}
      break;
    case WLZ_BO_GT:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ > val;}
      break;
    case WLZ_BO_GE:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ >= val;}
      break;
    case WLZ_BO_LT:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ < val;}
      break;
    case WLZ_BO_LE:
      while( buflen-- ){*outPP.p.dbp++ = *inbuf1++ <= val;}
      break;
    case WLZ_BO_MODULUS: /* FALLTHROUGH */
    case WLZ_BO_AND: /* FALLTHROUGH */
    case WLZ_BO_OR: /* FALLTHROUGH */
    case WLZ_BO_XOR: /* FALLTHROUGH */
    default:
      return WLZ_ERR_BINARY_OPERATOR_TYPE;
    }
    break;
  default:
    return WLZ_ERR_GREY_TYPE;
  }

  return WLZ_ERR_NONE;
}

/* function:     WlzScalarBinaryOp    */
/*! 
* \ingroup      WlzArithmetic
* \brief        Perform the given binary operation for each grey
value in image o1 with val putting the result in o3.
The value table of o3 can be identical to that of o1
to allow overwriting.
This function assumes that the domains of each object
are the same. Unless overwriting is required then the
functions WlzImageAdd, WlzImageSubtract etc. should be
used. These return an object which is the intersection
of the input objects.
*
* \return       Woolz error.
* \param    o1	Input object.
* \param    pval	Operand value.
* \param    o3	Oject for the return values. Setting equal to
 <tt>o1</tt> means values will be overwritten.
 * \param    op	Opertor to be applied.
* \par      Source:
*                WlzScalarBinaryOp.c
*/
WlzErrorNum WlzScalarBinaryOp(
  WlzObject		*o1,
  WlzPixelV		pval,
  WlzObject		*o3,
  WlzBinaryOperatorType op)
{
  WlzErrorNum		errNum = WLZ_ERR_NONE;

  /* check pointers */
  if( (o1 == NULL) || (o3 == NULL) ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }


  /* check object types */
  if( errNum == WLZ_ERR_NONE ){
    switch( o1->type ){

    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
    case WLZ_TRANS_OBJ:
    case WLZ_EMPTY_OBJ:
      if( (o3->type != o1->type) ){
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
      }
      break;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;

    }
  }

  /* check domains and valuetables */
  if( errNum == WLZ_ERR_NONE ){
    switch( o1->type ){

    case WLZ_2D_DOMAINOBJ:
      if( (o1->domain.core == NULL) || (o3->domain.core == NULL) ){
	errNum = WLZ_ERR_DOMAIN_NULL;
	break;
      }
      if( (o3->domain.core != o1->domain.core) ){
	errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
      }
      if( (o1->values.core == NULL) || (o3->values.core == NULL) ){
	errNum = WLZ_ERR_VALUES_NULL;
	break;
      }
      break;
      
    case WLZ_3D_DOMAINOBJ:
      if( (o1->domain.core == NULL) || (o3->domain.core == NULL) ){
	errNum = WLZ_ERR_DOMAIN_NULL;
	break;
      }
      if( (o3->domain.core != o1->domain.core) ){
	errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
      }
      if( (o1->values.core == NULL) || (o3->values.core == NULL) ){
	errNum = WLZ_ERR_VALUES_NULL;
	break;
      }
      errNum = WlzScalarBinaryOp3d(o1, pval, o3, op);
      break;
      
    case WLZ_TRANS_OBJ:
      if( (o1->values.core == NULL) || (o3->values.core == NULL) ){
	errNum = WLZ_ERR_VALUES_NULL;
	break;
      }
      errNum = WlzScalarBinaryOp(o1->values.obj, pval, o3->values.obj, op);
      break;

    case WLZ_EMPTY_OBJ:
      break;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  /* scan each object */
  if( (errNum == WLZ_ERR_NONE) && (o1->type == WLZ_2D_DOMAINOBJ) ){
    WlzIntervalWSpace	iwsp1, iwsp3;
    WlzGreyWSpace	gwsp1, gwsp3;
    WlzPixelP		o1PP, o3PP;
    int			*o1Buf=NULL;
    int			i, ival;
    double		dval;

    o1PP.p.v = o3PP.p.v = NULL;
    o1PP.type = o3PP.type = WLZ_GREY_ERROR;
    /* initialise the workspaces */
    if( (errNum = WlzInitGreyScan(o1,&iwsp1,&gwsp1)) ){
      return errNum;
    }
    if( (errNum = WlzInitGreyScan(o3,&iwsp3,&gwsp3)) ){
      return errNum;
    }

    /* initialise the buffers */
    switch( gwsp1.pixeltype ){
    case WLZ_GREY_INT:
    case WLZ_GREY_DOUBLE:
      o1PP.type = gwsp1.pixeltype;
      break;

    case WLZ_GREY_SHORT:
    case WLZ_GREY_UBYTE:
      o1PP.type = WLZ_GREY_INT;
      o1Buf = (int *) AlcMalloc(sizeof(int) *
				(iwsp1.intdmn->lastkl-iwsp1.intdmn->kol1+1));
      o1PP.p.inp = o1Buf;
      break;

    case WLZ_GREY_FLOAT:
      o1PP.type = WLZ_GREY_DOUBLE;
      o1Buf = (int *) AlcMalloc(sizeof(double) *
				(iwsp1.intdmn->lastkl-iwsp1.intdmn->kol1+1));
      o1PP.p.inp = o1Buf;
      break;
    default:
      break;
    }

    switch( pval.type ){
    case WLZ_GREY_INT:
    case WLZ_GREY_DOUBLE:
      break;

    case WLZ_GREY_SHORT:
      ival = (int) pval.v.shv;
      pval.v.inv = ival;
      pval.type = WLZ_GREY_INT;
      break;

    case WLZ_GREY_UBYTE:
      ival = (int) pval.v.ubv;
      pval.v.inv = ival;
      pval.type = WLZ_GREY_INT;
      break;

    case WLZ_GREY_FLOAT:
      dval = (double) pval.v.flv;
      pval.v.inv = (float )dval;
      pval.type = WLZ_GREY_DOUBLE;
      break;
    default:
      break;
    }

    o3PP.type = gwsp3.pixeltype;

    while ((errNum = WlzNextGreyInterval(&iwsp1)) == WLZ_ERR_NONE) {
      WlzNextGreyInterval(&iwsp3);

      /* copy values to buffers */
      switch(  gwsp1.pixeltype ){
      case WLZ_GREY_INT:
	o1PP.p.inp = gwsp1.u_grintptr.inp;
	break;

      case WLZ_GREY_SHORT:
	for(i=0; i<iwsp1.colrmn; i++){
	  o1PP.p.inp[i] = *gwsp1.u_grintptr.shp++;
	}
	break;

      case WLZ_GREY_UBYTE:
	for(i=0; i<iwsp1.colrmn; i++){
	  o1PP.p.inp[i] = (int) (*gwsp1.u_grintptr.ubp++);
	}
	break;

      case WLZ_GREY_FLOAT:
	for(i=0; i<iwsp1.colrmn; i++){
	  o1PP.p.dbp[i] = *gwsp1.u_grintptr.flp++;
	}
	break;

      case WLZ_GREY_DOUBLE:
	o1PP.p.dbp = gwsp1.u_grintptr.dbp;
	break;

      default:
	break;
      }

      switch(  gwsp3.pixeltype ){
      case WLZ_GREY_INT:
	o3PP.p.inp = gwsp3.u_grintptr.inp;
	break;

      case WLZ_GREY_SHORT:
	o3PP.p.shp = gwsp3.u_grintptr.shp;
	break;

      case WLZ_GREY_UBYTE:
	o3PP.p.ubp = gwsp3.u_grintptr.ubp;
	break;

      case WLZ_GREY_FLOAT:
	o3PP.p.flp = gwsp3.u_grintptr.flp;
	break;

      case WLZ_GREY_DOUBLE:
	o3PP.p.dbp = gwsp3.u_grintptr.dbp;
	break;

      default:
	break;
      }

      /* apply binary operation */
      switch( o1PP.type ){
      case WLZ_GREY_INT:
	switch( pval.type ){
	case WLZ_GREY_INT:
	  errNum = WlzBufIntIntScalarBinaryOp(o1PP.p.inp, pval.v.inv,
						o3PP, iwsp1.colrmn, op);
	  break;
	case WLZ_GREY_DOUBLE:
	  errNum = WlzBufIntDblScalarBinaryOp(o1PP.p.inp, pval.v.dbv,
						o3PP, iwsp1.colrmn, op);
	  break;
	default:
	  break;
	}
	break;
      case WLZ_GREY_DOUBLE:
	switch( pval.type ){
	case WLZ_GREY_INT:
	  errNum = WlzBufDblIntScalarBinaryOp(o1PP.p.dbp, pval.v.inv,
						o3PP, iwsp1.colrmn, op);
	  break;
	case WLZ_GREY_DOUBLE:
	  errNum = WlzBufDblDblScalarBinaryOp(o1PP.p.dbp, pval.v.dbv,
						o3PP, iwsp1.colrmn, op);
	  break;
	default:
	  break;
	}
	break;
      default:
	break;
      }

      if( errNum != WLZ_ERR_NONE ){
	break;
      }
    }
    if(errNum == WLZ_ERR_EOO)		/* Reset error from end of intervals */
    {
      errNum = WLZ_ERR_NONE;
    }

    /* free space */
    if( o1Buf ){
      AlcFree( o1Buf );
    }

  }

  return errNum;
}


static WlzErrorNum WlzScalarBinaryOp3d(
  WlzObject	*o1,
  WlzPixelV	pval,
  WlzObject	*o3,
  WlzBinaryOperatorType op)
{
  WlzObject		*temp1, *temp3;
  WlzPlaneDomain	*pdom;
  WlzDomain		*domains;
  WlzValues		*values1, *values3;
  int			i, nplanes;
  WlzErrorNum		errNum = WLZ_ERR_NONE;

  /* no need to check the object pointer or type because this procedure
     can only be accessed via WlzImageAdd. The domain and valuetable
     types must be checked however */

  /* check planedomain type */
  switch( o1->domain.p->type ){
  case WLZ_PLANEDOMAIN_DOMAIN:
    break;
  default:
    errNum = WLZ_ERR_PLANEDOMAIN_TYPE;
    break;
  }

  /* check voxeltable type */
  switch( o1->values.vox->type ){
  case WLZ_VOXELVALUETABLE_GREY:
    break;
  default:
    errNum = WLZ_ERR_VOXELVALUES_TYPE;
    break;
  }

  switch( o3->values.vox->type ){
  case WLZ_VOXELVALUETABLE_GREY:
    break;
  default:
    errNum = WLZ_ERR_VOXELVALUES_TYPE;
    break;
  }

  /* apply the binary operation to each plane */

  /* note a NULL domain is currently legal (to be WLZ_EMPTY_DOMAIN) and
   corresponds to a WLZ_EMPTY_OBJ */

  pdom = o1->domain.p;
  domains = pdom->domains;
  values1 = o1->values.vox->values;
  values3 = o3->values.vox->values;
  nplanes = pdom->lastpl - pdom->plane1 + 1;

  for(i=0; i < nplanes; i++, domains++, values1++, values3++){
    if( (*domains).core == NULL ){
      continue;
    }

    if((temp1 = WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains, *values1,
			    NULL, NULL, &errNum)) != NULL){
      if((temp3 = WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains, *values3,
			     NULL, NULL, &errNum)) != NULL){
	errNum = WlzScalarBinaryOp(temp1, pval, temp3, op);
	WlzFreeObj( temp3 );
      }
      WlzFreeObj( temp1 );
    }
    if( errNum != WLZ_ERR_NONE ){
      break;
    }
  }

  return errNum;
}
