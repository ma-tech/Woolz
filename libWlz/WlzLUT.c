#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzLUT_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzLUT.c
* \author       Bill Hill
* \date         November 2011
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
* \brief	Look up tables for value transformations.
* \ingroup	WlzValuesFilters
*/
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <Wlz.h>

static WlzErrorNum 		WlzLUTTransformObj2D(
				  WlzObject *rObj,
				  WlzObject *gObj,
				  WlzObject *tObj,
				  int pln,
				  int dither);
static WlzErrorNum 		WlzLUTTransformObj3D(
				  WlzObject *rObj,
				  WlzObject *gObj,
				  WlzObject *tObj,
				  int dither);

/*!
* \return	New look up table object.
* \ingroup	WlzValuesFilters
* \brief	Creates a new look up table object by calling
* 		WlzMakeLUTDomain() and WlzMakeLUTValues(). The last
* 		bin of the look up table must be greater than the first.
* \param	vType			Value type, which must be either
* 					WLZ_GREY_INT or WLZ_GREY_RGBA.
* \param	bin1			First bin of the look up table.
* \param	lastbin			Last bin of the look up table.
* \param	maxVal			Number of values to allocate (may be
* 					zero if no look up table required).
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzMakeLUTObject(WlzGreyType vType, int bin1, int lastbin,
				  WlzErrorNum *dstErr)
{
  int		nVal;
  WlzDomain	dom;
  WlzValues	val;
  WlzObject	*obj = NULL;
  WlzErrorNum  	errNum = WLZ_ERR_NONE;

  dom.core = NULL;
  val.core = NULL;
  if((nVal = lastbin - bin1 + 1) <= 0)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    dom.lut = WlzMakeLUTDomain(bin1, lastbin, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    val.lut = WlzMakeLUTValues(vType, nVal, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    obj = WlzMakeMain(WLZ_LUT, dom, val, NULL, NULL, &errNum);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeDomain(dom);
    (void )WlzFreeLUTValues(val);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	New look up table domain.
* \ingroup	WlzValuesFilters
* \brief	Makes a new look up table domain. The returned data
* 		structure should be freed using WlzFreeDomain().
* \param	bin1			First bin of look up table.
* \param	lastbin			Last bin of look up table.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzLUTDomain	*WlzMakeLUTDomain(int bin1, int lastbin, WlzErrorNum *dstErr)
{
  WlzLUTDomain  *lDom = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((lDom = (WlzLUTDomain *)
             AlcCalloc(1, sizeof(WlzLUTDomain))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    lDom->type = WLZ_LUT;
    lDom->bin1 = bin1;
    lDom->lastbin = lastbin;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(lDom);
}

/*!
* \return	New look up table values.
* \ingroup	WlzValuesFilters
* \brief	Makes a new look up table values data structure. The
* 		returned data structure should be freed using
* 		WlzFreeLUTValues().
* \param	vType			Value type, which must be either
* 					WLZ_GREY_INT or WLZ_GREY_RGBA.
* \param	maxVal			Number of values to allocate (may be
* 					zero if no look up table required).
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzLUTValues	*WlzMakeLUTValues(WlzGreyType vType, int maxVal,
				  WlzErrorNum *dstErr)
{
  size_t	gSz;
  WlzLUTValues  *lVal = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((gSz = WlzGreySize(vType)) == 0)
  {
    errNum = WLZ_ERR_GREY_TYPE;
  }
  else
  {
    switch(vType)
    {
      case WLZ_GREY_INT:  /* FALLTHROUGH */
      case WLZ_GREY_RGBA:
        break;
      default:
	errNum = WLZ_ERR_GREY_TYPE;
        break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(maxVal < 0)
    {
      maxVal = 0;
    }
    if((lVal = (WlzLUTValues *)
	       AlcCalloc(1, sizeof(WlzLUTValues))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      lVal->type = WLZ_LUT;
      lVal->vType = vType;
      lVal->maxVal = maxVal;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (maxVal > 0))
  {
    if((lVal->val.v = AlcCalloc(1, gSz * maxVal)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      lVal->freeptr = AlcFreeStackPush(lVal->freeptr, lVal->val.v, NULL);
    }
  }
  if((errNum != WLZ_ERR_NONE) && lVal)
  {
    AlcFree(lVal);
    lVal = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(lVal);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzValuesFilters
* \brief	Frees a look up table values data structure.
* \param	val			Given look up table values.
*/
WlzErrorNum	WlzFreeLUTValues(WlzValues val)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(val.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(val.core->type != WLZ_LUT)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    if(WlzUnlink(&(val.lut->linkcount), &errNum))
    {
      if(val.lut->freeptr)
      {
        (void )AlcFreeStackFree(val.lut->freeptr);
      }
      AlcFree(val.lut);
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzValuesFilters
* \brief	Creates a new look up table object, setting it's look up
* 		table values using a grey transform. See also
* 		WlzMakeLUTObject() and WlzLUTGreyTransformSet().
* \param	gTrType			Grey transform with valid transforms:
* 					WLZ_GREYTRANSFORMTYPE_IDENTITY,
* 					WLZ_GREYTRANSFORMTYPE_LINEAR,
*					WLZ_GREYTRANSFORMTYPE_GAMMA and
*					WLZ_GREYTRANSFORMTYPE_SIGMOID.
* \param	gType			Grey type for the look up table
* 					values and supplied minimum and
* 					maximum values..
* \param	il			Input minimum grey value, must be
* 					integral.
* \param	iu			Input maximum grey value, must be
* 					integral.
* \param	gol			Output minimum grey value.
* \param	gou			Output maximum grey value.
* \param	p0			First parameter, with meaning for
* 					transforms: unused, unused, gamma,
* \param	p1			Second parameter, with meaning for
* 					transforms: unused, unused, gamma,
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzLUTGreyTransformNew(WlzGreyTransformType gTrType,
					WlzGreyType gType,
					int il, int iu,
					WlzGreyV gol, WlzGreyV gou,
				        double p0, double p1,
					WlzErrorNum *dstErr)
{
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(errNum == WLZ_ERR_NONE)
  {
    obj = WlzMakeLUTObject(gType, il, iu, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzLUTGreyTransformSet(obj, gTrType, gType, il, iu, gol, gou,
                                    p0, p1);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(obj);
    obj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzValuesFilters
* \brief	Sets a given look up table using a grey transform.
* \param	obj			Given look up table object with
* 					the domain bin1 and lastbin set
* 					along with the values first and
* 					last bin values. All types, fields
* 					and arrays must be appropriate.
* \param	gTrType			Grey transform with valid transforms:
* 					WLZ_GREYTRANSFORMTYPE_IDENTITY,
* 					WLZ_GREYTRANSFORMTYPE_LINEAR,
*					WLZ_GREYTRANSFORMTYPE_GAMMA and
*					WLZ_GREYTRANSFORMTYPE_SIGMOID.
* \param	il			First integer value to set in LUT.
* \param	iu			Last integer value to set in LUT.
* \param	p0			First parameter, with meaning for
* 					transforms: unused, unused, gamma,
					mean.
* \param	p1			Second parameter, with meaning for
* 					transforms: unused, unused, unused,
					sigma (minimum clamped to 1e-6).
*/
WlzErrorNum 	WlzLUTGreyTransformSet(WlzObject *obj,
				       WlzGreyTransformType gTrType,
				       WlzGreyType gType, int il, int iu,
				       WlzGreyV gol, WlzGreyV gou,
				       double p0, double p1)
{
  int		kl,
		ku,
		ll,
  		lu,
		nk,
		nl;
  WlzLUTDomain  *lutDom;
  WlzLUTValues  *lutVal;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	eps = 1.0e-06;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->type != WLZ_LUT)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(obj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(obj->domain.lut->type != WLZ_LUT)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(obj->values.lut->type != WLZ_LUT)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    lutDom = obj->domain.lut;
    lutVal = obj->values.lut;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(gType != lutVal->vType) 
    {
      errNum = WLZ_ERR_GREY_TYPE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    ll = lutDom->bin1;
    lu = lutDom->lastbin;
    kl = ALG_MAX(ll, il) - ll;
    ku = ALG_MIN(lu, iu) - ll;
    nl = lu - ll + 1;
    nk = ku - kl + 1;
    if(nk > 0)
    {
      switch(gTrType)
      {
	case WLZ_GREYTRANSFORMTYPE_IDENTITY:
	  switch(gType)
	  {
	    case WLZ_GREY_INT:
	      {
		int	i;
		int	*p;

		p = lutVal->val.inp + kl;
		for(i = 0; i < nk; ++i)
		{
		  p[i] = kl + i;
		}
	      }
	      break;
	    case WLZ_GREY_RGBA:
	      {
		int	g,
		      	i;
		unsigned int *p;

		p = lutVal->val.rgbp + kl;;
		for(i = 0; i < nk; ++i)
		{
		  g = kl + i;
		  g = WLZ_CLAMP(g, 0, 255);
		  WLZ_RGBA_RGBA_SET(p[i], g, g, g, 255);
		}
	      }
	      break;
	    default:
	      errNum = WLZ_ERR_GREY_TYPE;
	      break;
	  }
	  break;
	case WLZ_GREYTRANSFORMTYPE_LINEAR:
	  switch(gType)
	  {
	    case WLZ_GREY_INT:
	      {
		int	i;
		double  t;
		int	*p;

		p = lutVal->val.inp + kl;
		t = (double )(gou.inv - gol.inv + 1) / (double )nk;
		for(i = 0; i < nk; ++i)
		{
		  p[i] = gol.inv + (int )((double )(kl + i) * t + eps); 
		}
	      }
	      break;
	    case WLZ_GREY_RGBA:
	      {
		int	i,
			j;
		double  w;
		int	g[4];
		double  ol[4],
			t[4];
		unsigned int *p;

		w = (double )nl;
		p = lutVal->val.rgbp + kl;
		ol[0] = WLZ_RGBA_RED_GET(gol.rgbv);
		ol[1] = WLZ_RGBA_GREEN_GET(gol.rgbv);
		ol[2] = WLZ_RGBA_BLUE_GET(gol.rgbv);
		ol[3] = WLZ_RGBA_ALPHA_GET(gol.rgbv);
		t[0] = ((double )(WLZ_RGBA_RED_GET(gou.rgbv))   - ol[0]) / w;
		t[1] = ((double )(WLZ_RGBA_GREEN_GET(gou.rgbv)) - ol[1]) / w;
		t[2] = ((double )(WLZ_RGBA_BLUE_GET(gou.rgbv))  - ol[2]) / w;
		t[3] = ((double )(WLZ_RGBA_ALPHA_GET(gou.rgbv)) - ol[3]) / w;
		for(i = 0; i < nk; ++i)
		{
		  for(j = 0; j < 4; ++j)
		  {
		    g[j] = ol[j] + (int )((double )(kl + i) * t[j] + eps);
		    g[j] = WLZ_CLAMP(g[j], 0, 255);
		  }
		  WLZ_RGBA_RGBA_SET(p[i], g[0], g[1], g[2], g[3]);
		}
	      }
	      break;
	    default:
	      errNum = WLZ_ERR_GREY_TYPE;
	      break;
	  }
	  break;
	case WLZ_GREYTRANSFORMTYPE_GAMMA:
	  switch(gType)
	  {
	    case WLZ_GREY_INT:
	      {
		int	i;
		double  ol,
			t;
		int	*p;

		p = lutVal->val.inp + kl;;
		ol = gol.inv + eps;
		t = (gou.inv - gol.inv + 1) * pow(nl, -p0);
		for(i = 0; i < nk; ++i)
		{
		  p[i] = (int )(t * pow((kl + i), p0) + ol);
		}
	      }
	      break;
	    case WLZ_GREY_RGBA:
	      {
		int	i;
		double  q;
		int	g[4];
		double  ol[4],
			t[4];
		unsigned int *p;

		p = lutVal->val.rgbp + kl;
		ol[0] = WLZ_RGBA_RED_GET(gol.rgbv)   + eps;
		ol[1] = WLZ_RGBA_GREEN_GET(gol.rgbv) + eps;
		ol[2] = WLZ_RGBA_BLUE_GET(gol.rgbv)  + eps;
		ol[3] = WLZ_RGBA_ALPHA_GET(gol.rgbv) + eps;
		q = pow(nl, -p0);
		t[0] = (WLZ_RGBA_RED_GET(gou.rgbv)   - ol[0]) * q;
		t[1] = (WLZ_RGBA_GREEN_GET(gou.rgbv) - ol[1]) * q;
		t[2] = (WLZ_RGBA_BLUE_GET(gou.rgbv)  - ol[2]) * q;
		t[3] = (WLZ_RGBA_ALPHA_GET(gou.rgbv) - ol[3]) * q;
		for(i = 0; i < nk; ++i)
		{
		  int	 j;

		  q = pow(i, -p0);
		  for(j = 0; j < 4; ++j)
		  {
		    g[j] = (int )(t[j] * q + ol[j]);
		    g[j] = WLZ_CLAMP(g[j], 0, 255);
		  }
		  WLZ_RGBA_RGBA_SET(p[i], g[0], g[1], g[2], g[3]);
		}
	      }
	      break;
	    default:
	      errNum = WLZ_ERR_GREY_TYPE;
	      break;
	  }
	  break;
	case WLZ_GREYTRANSFORMTYPE_SIGMOID:
	  {
	    double	fl,
	    		fu;

	    p1 = (p1 < DBL_EPSILON)? eps: 1.0 / (p1 * p1);  /* sigma^-2) */
	    fl =  1.0 / (1.0 + exp(-(il - p0) *  p1));
	    fu = 1.0 / ((1.0 / (1.0 + exp(-(iu - p0) * p1))) - fl);
	    switch(gType)
	    {
	      case WLZ_GREY_INT:
		{
		  int	i,
			g;
		  double t0,
			 t1;
		  int	*p;

		  p = lutVal->val.inp + kl;
		  t0 = (gou.inv - gol.inv + 1) * fu;
		  t1 = gol.inv - t0 * fl + eps;
		  for(i = 0; i < nk; ++i)
		  {
		    g = kl + i;
		    p[i] = (int )((t0 / (1.0 + exp(-(g - p0) * p1))) + t1);
		  }
		}
		break;
	      case WLZ_GREY_RGBA:
		{
		  int	i;
		  int	g[4];
		  double t0[4],
			 t1[4];
		  unsigned int *p;

		  p = lutVal->val.rgbp + kl;
		  t0[0] = (WLZ_RGBA_RED_GET(gou.rgbv) -
			   WLZ_RGBA_RED_GET(gol.rgbv)) * fu;
		  t0[1] = (WLZ_RGBA_GREEN_GET(gou.rgbv) -
			   WLZ_RGBA_GREEN_GET(gol.rgbv)) * fu;
		  t0[2] = (WLZ_RGBA_BLUE_GET(gou.rgbv) -
			   WLZ_RGBA_BLUE_GET(gol.rgbv)) * fu;
		  t0[3] = (WLZ_RGBA_ALPHA_GET(gou.rgbv) -
			   WLZ_RGBA_ALPHA_GET(gol.rgbv)) * fu;
		  t1[0] = WLZ_RGBA_RED_GET(gol.rgbv)   - t0[0] * fl + eps;
		  t1[1] = WLZ_RGBA_GREEN_GET(gol.rgbv) - t0[1] * fl + eps;
		  t1[2] = WLZ_RGBA_BLUE_GET(gol.rgbv)  - t0[2] * fl + eps;
		  t1[3] = WLZ_RGBA_ALPHA_GET(gol.rgbv) - t0[3] * fl + eps;
		  for(i = 0; i < nk; ++i)
		  {
		    int	j;
		    double q;

		    q = kl + i;
		    q = 1.0 / (1.0 + exp(-(q - p0) * p1));
		    for(j = 0; j < 4; ++j)
		    {
		      g[j] = (int )((t0[j] * q) + t1[j]);
		      g[j] = WLZ_CLAMP(g[j], 0, 255);
		    }
		    WLZ_RGBA_RGBA_SET(p[i], g[0], g[1], g[2], g[3]);
		  }
		}
		break;
	      default:
		errNum = WLZ_ERR_GREY_TYPE;
		break;
	    }
	  }
	  break;
	default:
	  errNum = WLZ_ERR_TRANSFORM_TYPE;
	  break;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Merged LUT object.
* \ingroup	WlzValuesFilters
* \brief	Merges separate grey LUTs into a single RGBA LUT. All LUTs
* 		should either be NULL or have the same domain and integer
* 		LUT values, with at least one LUT object being non NULL.
* \param	r			Red LUT, if NULL a constant value
* 					of 0 is assumed.
* \param	g			Green LUT, if NULL a constant value
* 					of 0 is assumed.
* \param	b			Blue LUT, if NULL a constant value
* 					of 0 is assumed.
* \param	a			Alpha LUT, if NULL a constant value
* 					of 255 is assumed.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject       *WlzLUTMergeToRGBA(WlzObject *r, WlzObject *g, WlzObject *b,
				   WlzObject *a, WlzErrorNum *dstErr)
{
  int		i;
  WlzObject	*mObj = NULL;
  WlzLUTDomain	*lutDom;
  unsigned int	*mP;
  int		*cP[4];
  WlzObject	*cObj[4];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  cObj[0] = r;
  cObj[1] = g;
  cObj[2] = b;
  cObj[3] = a;
  lutDom = NULL;
  /* Check the given LUT objects. */
  for(i = 0; i < 4; ++i)
  {
    if(cObj[i])
    {
      if(cObj[i]->type != WLZ_LUT)
      {
        errNum = WLZ_ERR_OBJECT_TYPE;
      }
      else if(cObj[i]->domain.core == NULL)
      {
        errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else if(cObj[i]->values.core == NULL)
      {
        errNum = WLZ_ERR_VALUES_NULL;
      }
      else if(cObj[i]->domain.core->type != WLZ_LUT)
      {
        errNum = WLZ_ERR_DOMAIN_TYPE;
      }
      else if(cObj[i]->values.core->type != WLZ_LUT)
      {
        errNum = WLZ_ERR_VALUES_TYPE;
      }
      else if(cObj[i]->values.lut->vType != WLZ_GREY_INT)
      {
        errNum = WLZ_ERR_VALUES_DATA;
      }
      else
      {
	cP[i] = cObj[i]->values.lut->val.inp;
	if(lutDom == NULL)
	{
	  lutDom = cObj[i]->domain.lut;
	}
	else if(lutDom)
	{
	  if((lutDom->bin1 != cObj[i]->domain.lut->bin1) ||
	      (lutDom->lastbin != cObj[i]->domain.lut->lastbin))
	  {
	    errNum = WLZ_ERR_DOMAIN_DATA;
	  }
	}
      }
      if(errNum != WLZ_ERR_NONE)
      {
        break;
      }
    }
  }
  if(lutDom == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  /* Create merged RGBA LUT object. */
  if(errNum == WLZ_ERR_NONE)
  {
    mObj = WlzMakeLUTObject(WLZ_GREY_RGBA, lutDom->bin1, lutDom->lastbin,
			    &errNum);
  }
  /* Merge the LUT values. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		j,
    		n;
    int		c[4];

    mP = mObj->values.lut->val.rgbp;
    n = lutDom->lastbin - lutDom->bin1 + 1;
    for(i = 0; i < n; ++i)
    {
      for(j = 0; j < 3; ++j)
      {
	c[j] = (cObj[j])? WLZ_CLAMP(cP[j][i], 0, 255): 0;
      }
      c[3] = (cObj[3])? WLZ_CLAMP(cP[3][i], 0, 255): 255;
      WLZ_RGBA_RGBA_SET(mP[i], c[0], c[1], c[2], c[3]);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mObj);
}

/*!
* \return	Value transformed object.
* \ingroup	WlzValuesFilters
* \brief	Applies a LUT to an domain object with grey values.
* \param	gObj			Given domain object with values.
* \param	tObj			Given LUT transform object.
* \param	rGType			Grey type for the return object.
* 					Ignored if the operation is done
* 					in place.
* \param	inplace			Flag, which if non-zero, specifies
* 					that the given objects grey values
* 					should be modified in place. If set
* 					then the return grey type parameter
* 					is ignored.
* \param	dither			Non-zero if value dithering is
* 					required. Dithering may make the value
* 					transformation significantly slower.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzLUTTransformObj(WlzObject *gObj, WlzObject *tObj,
				    WlzGreyType rGType,
				    int inplace, int dither,
				    WlzErrorNum *dstErr)
{
  WlzPixelV	bkg;
  WlzObject	*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((gObj == NULL) || (tObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((gObj->domain.core == NULL) || (tObj->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((gObj->values.core == NULL) || (tObj->values.core == NULL))
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if((tObj->type != WLZ_LUT) ||
          ((gObj->type != WLZ_2D_DOMAINOBJ) &&
	   (gObj->type != WLZ_3D_DOMAINOBJ)))
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bkg = WlzGetBackground(gObj, &errNum);
  }
  if((errNum == WLZ_ERR_NONE) && inplace)
  {
    rGType = WlzGreyTypeFromObj(gObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bkg = WlzLUTTransformPixelValue(tObj, bkg, rGType, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Create new grey valued return object if required. */
    if(inplace)
    {
      rObj = gObj;
    }
    else
    {
      WlzValues	rVal;
      WlzObjectType rTabType;

      rVal.core = NULL;
      rTabType = WlzGreyTableType(WLZ_GREY_TAB_RAGR, rGType, NULL);
      if(errNum == WLZ_ERR_NONE)
      {
	if(gObj->type == WLZ_2D_DOMAINOBJ)
	{
          rVal.v = WlzNewValueTb(gObj, rTabType, bkg, &errNum);
	}
	else /* gObj->type == WLZ_3D_DOMAINOBJ */
	{
	  rVal.vox = WlzNewValuesVox(gObj, rTabType, bkg, &errNum);
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
        rObj = WlzMakeMain(gObj->type, gObj->domain, rVal, NULL, NULL,
			   &errNum);
      }
      if((errNum != WLZ_ERR_NONE) && (rVal.core != NULL))
      {
        if(gObj->type == WLZ_2D_DOMAINOBJ)
	{
	  (void )WlzFreeValueTb(rVal.v);
	}
	else /* gObj->type == WLZ_3D_DOMAINOBJ */
	{
	  (void )WlzFreeVoxelValueTb(rVal.vox);
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(gObj->type == WLZ_2D_DOMAINOBJ)
    {
      errNum = WlzLUTTransformObj2D(rObj, gObj, tObj, 0, dither);
    }
    else
    {
      errNum = WlzLUTTransformObj3D(rObj, gObj, tObj, dither);
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(inplace == 0)
    {
      (void )WlzFreeObj(rObj);
    }
    rObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzValuesFilters
* \brief	Applies a grey value LUT transform to the values of the
* 		given 2D domain object setting them in the object for
* 		return. The domains of the two grey valued object must
* 		be the same (assumed, not tested).
* \param	rObj			Object with values to be set.
* \param	gObj			Given 2D domain object with values to
* 					be transformed.
* \param	tObj			LUT object.
* \param	pln			Plane value, only used with 3D value
* 					tables.
* \param	dither			Flag, dither if non-zero.
*/
static WlzErrorNum WlzLUTTransformObj2D(WlzObject *rObj, WlzObject *gObj,
					WlzObject *tObj, int pln, int dither)
{
  WlzGreyWSpace gGWSp,
  	        rGWSp;
  WlzIntervalWSpace gIWSp,
  		    rIWSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  errNum = WlzInitGreyScan(gObj, &gIWSp, &gGWSp);
  if(errNum == WLZ_ERR_NONE)
  {
    gIWSp.plnpos = pln;
    errNum = WlzInitGreyScan(rObj, &rIWSp, &rGWSp);
    if(errNum == WLZ_ERR_NONE)
    {
      rIWSp.plnpos = pln;
      while((errNum == WLZ_ERR_NONE) &&
	    ((errNum = WlzNextGreyInterval(&gIWSp)) == WLZ_ERR_NONE) &&
	    ((errNum = WlzNextGreyInterval(&rIWSp)) == WLZ_ERR_NONE))
      {
	errNum = WlzLUTTransformGreyValues(tObj,
					   rGWSp.u_grintptr, rGWSp.pixeltype,
					   gGWSp.u_grintptr, gGWSp.pixeltype,
					   gIWSp.rgtpos - gIWSp.lftpos + 1,
					   dither);
      }
      (void )WlzEndGreyScan(&rIWSp, &rGWSp);
    }
    (void )WlzEndGreyScan(&gIWSp, &gGWSp);
  }
  if(errNum == WLZ_ERR_EOO)
  {
    errNum = WLZ_ERR_NONE;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzValuesFilters
* \brief	Applies a grey value LUT transform to the values of the
* 		given 3D domain object setting them in the object for
* 		return. The domains of the two grey valued object must
* 		the same (assumed, but not tested).
* \param	rObj			Object with values to be set.
* \param	gObj			Given 3D domain object with values to
* 					be transformed.
* \param	tObj			LUT object.
* \param	dither			Flag, dither if non-zero.
*/
static WlzErrorNum WlzLUTTransformObj3D(WlzObject *rObj, WlzObject *gObj,
				        WlzObject *tObj, int dither)
{
  int		pln,
  		gTiled;
  WlzPlaneDomain *gDom,
  		 *rDom;
  WlzValues      gVal,
  		 rVal;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  gDom = gObj->domain.p;
  rDom = rObj->domain.p;
  gVal = gObj->values;
  rVal = rObj->values; 			/* Know that rObj doesn't have tiled
  				           values since it was created in
					   the calling function. */
  gTiled = WlzGreyTableIsTiled(gVal.core->type);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(pln = gDom->plane1; pln <= gDom->lastpl; ++pln)
  {
    if(errNum == WLZ_ERR_NONE)
    {
      WlzDomain *gDom2,
      		*rDom2;
      WlzValues *gVal2,
      		*rVal2;
      WlzErrorNum errNum2 = WLZ_ERR_NONE;

      gDom2 = gDom->domains + pln - gDom->plane1;
      rDom2 = rDom->domains + pln - rDom->plane1;
      if(gTiled == 0)
      {
        gVal2 = gVal.vox->values  + pln - gDom->plane1;
      }
      rVal2 = rVal.vox->values  + pln - rDom->plane1;
      if(gDom2 && rDom2 && (gTiled || gVal2) && rVal2)
      {
        WlzObject *gObj2 = NULL,
		  *rObj2 = NULL;

        if(((gObj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ, *gDom2,
	                         (gTiled)? gVal: *gVal2, NULL, NULL,
	                         &errNum2)) != NULL) &&
           ((rObj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ, *rDom2, *rVal2, NULL, NULL,
	                         &errNum2)) != NULL))
        {
          errNum2 = WlzLUTTransformObj2D(rObj2, gObj2, tObj, pln, dither);
        }
        (void )WlzFreeObj(gObj2);
        (void )WlzFreeObj(rObj2);
      }
#ifdef _OPENMP
#pragma omp critical
      {
#endif
        if(errNum2 != WLZ_ERR_NONE)
        {
          errNum = errNum2;
        }
#ifdef _OPENMP
      }
#endif
    }
  }
  return(errNum);
}

/*!
* \return	Transformed pixel value.
* \ingroup	WlzValuesFilters
* \brief	Transforms the grey value of a single pixel using the given
* 		look up table object.
* 		This function is not very efficient, use
* 		WlzLUTTransformGreyValues() (which this function calls)
* 		where possible.
* \param	tObj			Look up table transform object.
* \param	gPix			Given pixel value.
* \param	rType			Return grey type.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzPixelV	WlzLUTTransformPixelValue(WlzObject *tObj,
				WlzPixelV gPix, WlzGreyType rType,
				WlzErrorNum *dstErr)
{
  WlzGreyP 	gP,
  		rP;
  WlzPixelV	rPix;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  rPix.type = rType;
  switch(rPix.type)
  {
    case WLZ_GREY_INT:
      rP.inp = &(rPix.v.inv);
      break;
    case WLZ_GREY_SHORT:
      rP.shp = &(rPix.v.shv);
      break;
    case WLZ_GREY_UBYTE:
      rP.ubp = &(rPix.v.ubv);
      break;
    case WLZ_GREY_RGBA:
      rP.rgbp = &(rPix.v.rgbv);
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(gPix.type)
    {
      case WLZ_GREY_INT:
	gP.inp = &(gPix.v.inv);
	break;
      case WLZ_GREY_SHORT:
	gP.shp = &(gPix.v.shv);
	break;
      case WLZ_GREY_UBYTE:
	gP.ubp = &(gPix.v.ubv);
	break;
      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzLUTTransformGreyValues(tObj, rP, rType, gP, gPix.type, 1, 0);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rPix);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzValuesFilters
* \brief	Transforms the grey values in the given buffer using the
* 		look up table object, putting the values into the return
* 		buffer.
* \param	tObj			Look up table transform object.
* \param	rP			Return grey value buffer.
* \param	rType			Return grey type.
* \param	gP			Given grey value buffer.
* \param	gType			Given grey type.
* \param	nVal			Number of grey values in buffer.
* \param	dither			Flag, dither if non-zero.
*/
WlzErrorNum 	WlzLUTTransformGreyValues(WlzObject *tObj,
					  WlzGreyP rP, WlzGreyType rType,
					  WlzGreyP gP, WlzGreyType gType,
					  int nVal, int dither)
{
  int		bin1,
		lastbin,
		maxBin;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(tObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(tObj->type != WLZ_LUT)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(tObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(tObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(nVal > 0)
  {
    bin1 = tObj->domain.lut->bin1;
    lastbin = tObj->domain.lut->lastbin;
    maxBin = lastbin - bin1;      /* lastbin - bin1 *not* lastbin - bin1 + 1 */
    switch(tObj->values.lut->vType)
    {
      case WLZ_GREY_INT:
	{
	  int	*lut;

	  lut = tObj->values.lut->val.inp;
	  if(dither)
	  {
	    int		i,
	    		k;
	    int		m[3];

	    switch(rType)
	    {
	      case WLZ_GREY_INT:
		switch(gType)
		{
		  case WLZ_GREY_INT:
		    for(i = 0; i < nVal; ++i)
		    {
		      m[0] = WLZ_CLAMP(gP.inp[i] - bin1 - 1, 0, maxBin);
		      m[1] = WLZ_CLAMP(gP.inp[i] - bin1,     0, maxBin);
		      m[2] = WLZ_CLAMP(gP.inp[i] - bin1 + 1, 0, maxBin);
		      k = WlzValueDitherI(lut[m[0]], lut[m[1]], lut[m[2]]);
		      rP.inp[i] = k;
		    }
		    break;
		  case WLZ_GREY_SHORT:
		    for(i = 0; i < nVal; ++i)
		    {
		      m[0] = WLZ_CLAMP(gP.shp[i] - bin1 - 1, 0, maxBin);
		      m[1] = WLZ_CLAMP(gP.shp[i] - bin1,     0, maxBin);
		      m[2] = WLZ_CLAMP(gP.shp[i] - bin1 + 1, 0, maxBin);
		      k = WlzValueDitherI(lut[m[0]], lut[m[1]], lut[m[2]]);
		      rP.inp[i] = k;
		    }
		    break;
		  case WLZ_GREY_UBYTE:
		    for(i = 0; i < nVal; ++i)
		    {
		      m[0] = WLZ_CLAMP(gP.ubp[i] - bin1 - 1, 0, maxBin);
		      m[1] = WLZ_CLAMP(gP.ubp[i] - bin1,     0, maxBin);
		      m[2] = WLZ_CLAMP(gP.ubp[i] - bin1 + 1, 0, maxBin);
		      k = WlzValueDitherI(lut[m[0]], lut[m[1]], lut[m[2]]);
		      rP.inp[i] = k;
		    }
		    break;
		  default:
		    errNum = WLZ_ERR_GREY_TYPE;
		    break;
		}
		break;
	      case WLZ_GREY_SHORT:
		switch(gType)
		{
		  case WLZ_GREY_INT:
		    for(i = 0; i < nVal; ++i)
		    {
		      m[0] = WLZ_CLAMP(gP.inp[i] - bin1 - 1, 0, maxBin);
		      m[1] = WLZ_CLAMP(gP.inp[i] - bin1,     0, maxBin);
		      m[2] = WLZ_CLAMP(gP.inp[i] - bin1 + 1, 0, maxBin);
		      k = WlzValueDitherI(lut[m[0]], lut[m[1]], lut[m[2]]);
		      rP.shp[i] = WLZ_CLAMP(k, SHRT_MIN, SHRT_MAX);
		    }
		    break;
		  case WLZ_GREY_SHORT:
		    for(i = 0; i < nVal; ++i)
		    {
		      m[0] = WLZ_CLAMP(gP.shp[i] - bin1 - 1, 0, maxBin);
		      m[1] = WLZ_CLAMP(gP.shp[i] - bin1,     0, maxBin);
		      m[2] = WLZ_CLAMP(gP.shp[i] - bin1 + 1, 0, maxBin);
		      k = WlzValueDitherI(lut[m[0]], lut[m[1]], lut[m[2]]);
		      rP.shp[i] = WLZ_CLAMP(k, SHRT_MIN, SHRT_MAX);
		    }
		    break;
		  case WLZ_GREY_UBYTE:
		    for(i = 0; i < nVal; ++i)
		    {
		      m[0] = WLZ_CLAMP(gP.ubp[i] - bin1 - 1, 0, maxBin);
		      m[1] = WLZ_CLAMP(gP.ubp[i] - bin1,     0, maxBin);
		      m[2] = WLZ_CLAMP(gP.ubp[i] - bin1 + 1, 0, maxBin);
		      k = WlzValueDitherI(lut[m[0]], lut[m[1]], lut[m[2]]);
		      rP.shp[i] = WLZ_CLAMP(k, SHRT_MIN, SHRT_MAX);
		    }
		    break;
		  default:
		    errNum = WLZ_ERR_GREY_TYPE;
		    break;
		}
		break;
	      case WLZ_GREY_UBYTE:
		switch(gType)
		{
		  case WLZ_GREY_INT:
		    for(i = 0; i < nVal; ++i)
		    {
		      m[0] = WLZ_CLAMP(gP.inp[i] - bin1 - 1, 0, maxBin);
		      m[1] = WLZ_CLAMP(gP.inp[i] - bin1,     0, maxBin);
		      m[2] = WLZ_CLAMP(gP.inp[i] - bin1 + 1, 0, maxBin);
		      k = WlzValueDitherI(lut[m[0]], lut[m[1]], lut[m[2]]);
		      rP.ubp[i] = WLZ_CLAMP(k, 0, 255);
		    }
		    break;
		  case WLZ_GREY_SHORT:
		    for(i = 0; i < nVal; ++i)
		    {
		      m[0] = WLZ_CLAMP(gP.shp[i] - bin1 - 1, 0, maxBin);
		      m[1] = WLZ_CLAMP(gP.shp[i] - bin1,     0, maxBin);
		      m[2] = WLZ_CLAMP(gP.shp[i] - bin1 + 1, 0, maxBin);
		      k = WlzValueDitherI(lut[m[0]], lut[m[1]], lut[m[2]]);
		      rP.ubp[i] = WLZ_CLAMP(k, 0, 255);
		    }
		    break;
		  case WLZ_GREY_UBYTE:
		    for(i = 0; i < nVal; ++i)
		    {
		      m[0] = WLZ_CLAMP(gP.ubp[i] - bin1 - 1, 0, maxBin);
		      m[1] = WLZ_CLAMP(gP.ubp[i] - bin1,     0, maxBin);
		      m[2] = WLZ_CLAMP(gP.ubp[i] - bin1 + 1, 0, maxBin);
		      k = WlzValueDitherI(lut[m[0]], lut[m[1]], lut[m[2]]);
		      rP.ubp[i] = WLZ_CLAMP(k, 0, 255);
		    }
		    break;
		  default:
		    errNum = WLZ_ERR_GREY_TYPE;
		    break;
		}
		break;
	      case WLZ_GREY_RGBA:
		switch(gType)
		{
		  case WLZ_GREY_INT:
		    for(i = 0; i < nVal; ++i)
		    {
		      m[0] = WLZ_CLAMP(gP.inp[i] - bin1 - 1, 0, maxBin);
		      m[1] = WLZ_CLAMP(gP.inp[i] - bin1,     0, maxBin);
		      m[2] = WLZ_CLAMP(gP.inp[i] - bin1 + 1, 0, maxBin);
		      k = WlzValueDitherI(lut[m[0]], lut[m[1]], lut[m[2]]);
		      k = WLZ_CLAMP(k, 0, 255);
		      WLZ_RGBA_RGBA_SET(rP.rgbp[i], k, k, k, 255);

		    }
		    break;
		  case WLZ_GREY_SHORT:
		    for(i = 0; i < nVal; ++i)
		    {
		      m[0] = WLZ_CLAMP(gP.shp[i] - bin1 - 1, 0, maxBin);
		      m[1] = WLZ_CLAMP(gP.shp[i] - bin1,     0, maxBin);
		      m[2] = WLZ_CLAMP(gP.shp[i] - bin1 + 1, 0, maxBin);
		      k = WlzValueDitherI(lut[m[0]], lut[m[1]], lut[m[2]]);
		      k = WLZ_CLAMP(k, 0, 255);
		      WLZ_RGBA_RGBA_SET(rP.rgbp[i], k, k, k, 255);
		    }
		    break;
		  case WLZ_GREY_UBYTE:
		    for(i = 0; i < nVal; ++i)
		    {
		      m[0] = WLZ_CLAMP(gP.ubp[i] - bin1 - 1, 0, maxBin);
		      m[1] = WLZ_CLAMP(gP.ubp[i] - bin1,     0, maxBin);
		      m[2] = WLZ_CLAMP(gP.ubp[i] - bin1 + 1, 0, maxBin);
		      k = WlzValueDitherI(lut[m[0]], lut[m[1]], lut[m[2]]);
		      k = WLZ_CLAMP(k, 0, 255);
		      WLZ_RGBA_RGBA_SET(rP.rgbp[i], k, k, k, 255);
		    }
		    break;
		  default:
		    errNum = WLZ_ERR_GREY_TYPE;
		    break;
		}
		break;
	      default:
		errNum = WLZ_ERR_GREY_TYPE;
		break;
	    }
	  }
	  else /* dither == 0 */
	  {
	    int		i,
	    		j,
			k;

	    switch(rType)
	    {
	      case WLZ_GREY_INT:
		switch(gType)
		{
		  case WLZ_GREY_INT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.inp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      rP.inp[i] = lut[j];
		    }
		    break;
		  case WLZ_GREY_SHORT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.shp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      rP.inp[i] = lut[j];
		    }
		    break;
		  case WLZ_GREY_UBYTE:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.ubp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      rP.inp[i] = lut[j];
		    }
		    break;
		  default:
		    errNum = WLZ_ERR_GREY_TYPE;
		    break;
		}
		break;
	      case WLZ_GREY_SHORT:
		switch(gType)
		{
		  case WLZ_GREY_INT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.inp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      rP.shp[i] = WLZ_CLAMP(lut[j], SHRT_MIN, SHRT_MAX);
		    }
		    break;
		  case WLZ_GREY_SHORT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.shp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      rP.shp[i] = WLZ_CLAMP(lut[j], SHRT_MIN, SHRT_MAX);
		    }
		    break;
		  case WLZ_GREY_UBYTE:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.ubp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      rP.shp[i] = WLZ_CLAMP(lut[j], SHRT_MIN, SHRT_MAX);
		    }
		    break;
		  default:
		    errNum = WLZ_ERR_GREY_TYPE;
		    break;
		}
		break;
	      case WLZ_GREY_UBYTE:
		switch(gType)
		{
		  case WLZ_GREY_INT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.inp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      rP.ubp[i] = WLZ_CLAMP(lut[j], 0, 255);
		    }
		    break;
		  case WLZ_GREY_SHORT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.shp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      rP.ubp[i] = WLZ_CLAMP(lut[j], 0, 255);
		    }
		    break;
		  case WLZ_GREY_UBYTE:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.ubp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      rP.ubp[i] = WLZ_CLAMP(lut[j], 0, 255);
		    }
		    break;
		  default:
		    errNum = WLZ_ERR_GREY_TYPE;
		    break;
		}
		break;
	      case WLZ_GREY_RGBA:
		switch(gType)
		{
		  case WLZ_GREY_INT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.inp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      k = WLZ_CLAMP(lut[j], 0, 255);
		      WLZ_RGBA_RGBA_SET(rP.rgbp[i], k, k, k, 255);
		    }
		    break;
		  case WLZ_GREY_SHORT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.shp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      k = WLZ_CLAMP(lut[j], 0, 255);
		      WLZ_RGBA_RGBA_SET(rP.rgbp[i], k, k, k, 255);
		    }
		    break;
		  case WLZ_GREY_UBYTE:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.ubp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      k = WLZ_CLAMP(lut[j], 0, 255);
		      WLZ_RGBA_RGBA_SET(rP.rgbp[i], k, k, k, 255);
		    }
		    break;
		  default:
		    errNum = WLZ_ERR_GREY_TYPE;
		    break;
		}
		break;
	      default:
		errNum = WLZ_ERR_GREY_TYPE;
		break;
	    }
	  }
	}
	break;
      case WLZ_GREY_RGBA:
	{
	  unsigned int *lut;

	  lut = tObj->values.lut->val.rgbp;
	  if(dither)
	  {
	    int		i,
	    		j,
	    		k;
	    int		m[4];
	    unsigned int r;
	    unsigned int s[3];

	    switch(rType)
	    {
	      case WLZ_GREY_INT:
		switch(gType)
		{
		  case WLZ_GREY_INT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = WLZ_CLAMP(gP.inp[i] - bin1 - 1, 0, maxBin);
		      s[0] = lut[j];
		      j = WLZ_CLAMP(gP.inp[i] - bin1,     0, maxBin);
		      s[1] = lut[j];
		      j = WLZ_CLAMP(gP.inp[i] - bin1 + 1, 0, maxBin);
		      s[2] = lut[j];
		      m[0] = WlzValueDitherI(WLZ_RGBA_RED_GET(s[0]),
		                             WLZ_RGBA_RED_GET(s[1]),
					     WLZ_RGBA_RED_GET(s[2]));
		      m[1] = WlzValueDitherI(WLZ_RGBA_GREEN_GET(s[0]),
		                             WLZ_RGBA_GREEN_GET(s[1]),
					     WLZ_RGBA_GREEN_GET(s[2]));
		      m[2] = WlzValueDitherI(WLZ_RGBA_BLUE_GET(s[0]),
		                             WLZ_RGBA_BLUE_GET(s[1]),
					     WLZ_RGBA_BLUE_GET(s[2]));
		      m[3] = WlzValueDitherI(WLZ_RGBA_ALPHA_GET(s[0]),
		                             WLZ_RGBA_ALPHA_GET(s[1]),
					     WLZ_RGBA_ALPHA_GET(s[2]));
		      WLZ_RGBA_RGBA_SET(r, m[0], m[1], m[2], m[3]);
		      rP.inp[i] = WLZ_RGBA_MODULUS(r);
		    }
		    break;
		  case WLZ_GREY_SHORT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = WLZ_CLAMP(gP.shp[i] - bin1 - 1, 0, maxBin);
		      s[0] = lut[j];
		      j = WLZ_CLAMP(gP.shp[i] - bin1,     0, maxBin);
		      s[1] = lut[j];
		      j = WLZ_CLAMP(gP.shp[i] - bin1 + 1, 0, maxBin);
		      s[2] = lut[j];
		      m[0] = WlzValueDitherI(WLZ_RGBA_RED_GET(s[0]),
		                             WLZ_RGBA_RED_GET(s[1]),
					     WLZ_RGBA_RED_GET(s[2]));
		      m[1] = WlzValueDitherI(WLZ_RGBA_GREEN_GET(s[0]),
		                             WLZ_RGBA_GREEN_GET(s[1]),
					     WLZ_RGBA_GREEN_GET(s[2]));
		      m[2] = WlzValueDitherI(WLZ_RGBA_BLUE_GET(s[0]),
		                             WLZ_RGBA_BLUE_GET(s[1]),
					     WLZ_RGBA_BLUE_GET(s[2]));
		      m[3] = WlzValueDitherI(WLZ_RGBA_ALPHA_GET(s[0]),
		                             WLZ_RGBA_ALPHA_GET(s[1]),
					     WLZ_RGBA_ALPHA_GET(s[2]));
		      WLZ_RGBA_RGBA_SET(r, m[0], m[1], m[2], m[3]);
		      rP.inp[i] = WLZ_RGBA_MODULUS(r);
		    }
		    break;
		  case WLZ_GREY_UBYTE:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = WLZ_CLAMP(gP.ubp[i] - bin1 - 1, 0, maxBin);
		      s[0] = lut[j];
		      j = WLZ_CLAMP(gP.ubp[i] - bin1,     0, maxBin);
		      s[1] = lut[j];
		      j = WLZ_CLAMP(gP.ubp[i] - bin1 + 1, 0, maxBin);
		      s[2] = lut[j];
		      m[0] = WlzValueDitherI(WLZ_RGBA_RED_GET(s[0]),
		                             WLZ_RGBA_RED_GET(s[1]),
					     WLZ_RGBA_RED_GET(s[2]));
		      m[1] = WlzValueDitherI(WLZ_RGBA_GREEN_GET(s[0]),
		                             WLZ_RGBA_GREEN_GET(s[1]),
					     WLZ_RGBA_GREEN_GET(s[2]));
		      m[2] = WlzValueDitherI(WLZ_RGBA_BLUE_GET(s[0]),
		                             WLZ_RGBA_BLUE_GET(s[1]),
					     WLZ_RGBA_BLUE_GET(s[2]));
		      m[3] = WlzValueDitherI(WLZ_RGBA_ALPHA_GET(s[0]),
		                             WLZ_RGBA_ALPHA_GET(s[1]),
					     WLZ_RGBA_ALPHA_GET(s[2]));
		      WLZ_RGBA_RGBA_SET(r, m[0], m[1], m[2], m[3]);
		      rP.inp[i] = WLZ_RGBA_MODULUS(r);
		    }
		    break;
		  default:
		    errNum = WLZ_ERR_GREY_TYPE;
		    break;
		}
		break;
	      case WLZ_GREY_SHORT:
		switch(gType)
		{
		  case WLZ_GREY_INT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = WLZ_CLAMP(gP.inp[i] - bin1 - 1, 0, maxBin);
		      s[0] = lut[j];
		      j = WLZ_CLAMP(gP.inp[i] - bin1,     0, maxBin);
		      s[1] = lut[j];
		      j = WLZ_CLAMP(gP.inp[i] - bin1 + 1, 0, maxBin);
		      s[2] = lut[j];
		      m[0] = WlzValueDitherI(WLZ_RGBA_RED_GET(s[0]),
		                             WLZ_RGBA_RED_GET(s[1]),
					     WLZ_RGBA_RED_GET(s[2]));
		      m[1] = WlzValueDitherI(WLZ_RGBA_GREEN_GET(s[0]),
		                             WLZ_RGBA_GREEN_GET(s[1]),
					     WLZ_RGBA_GREEN_GET(s[2]));
		      m[2] = WlzValueDitherI(WLZ_RGBA_BLUE_GET(s[0]),
		                             WLZ_RGBA_BLUE_GET(s[1]),
					     WLZ_RGBA_BLUE_GET(s[2]));
		      m[3] = WlzValueDitherI(WLZ_RGBA_ALPHA_GET(s[0]),
		                             WLZ_RGBA_ALPHA_GET(s[1]),
					     WLZ_RGBA_ALPHA_GET(s[2]));
		      WLZ_RGBA_RGBA_SET(r, m[0], m[1], m[2], m[3]);
		      rP.shp[i] = WLZ_RGBA_MODULUS(r);
		    }
		    break;
		  case WLZ_GREY_SHORT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = WLZ_CLAMP(gP.shp[i] - bin1 - 1, 0, maxBin);
		      s[0] = lut[j];
		      j = WLZ_CLAMP(gP.shp[i] - bin1,     0, maxBin);
		      s[1] = lut[j];
		      j = WLZ_CLAMP(gP.shp[i] - bin1 + 1, 0, maxBin);
		      s[2] = lut[j];
		      m[0] = WlzValueDitherI(WLZ_RGBA_RED_GET(s[0]),
		                             WLZ_RGBA_RED_GET(s[1]),
					     WLZ_RGBA_RED_GET(s[2]));
		      m[1] = WlzValueDitherI(WLZ_RGBA_GREEN_GET(s[0]),
		                             WLZ_RGBA_GREEN_GET(s[1]),
					     WLZ_RGBA_GREEN_GET(s[2]));
		      m[2] = WlzValueDitherI(WLZ_RGBA_BLUE_GET(s[0]),
		                             WLZ_RGBA_BLUE_GET(s[1]),
					     WLZ_RGBA_BLUE_GET(s[2]));
		      m[3] = WlzValueDitherI(WLZ_RGBA_ALPHA_GET(s[0]),
		                             WLZ_RGBA_ALPHA_GET(s[1]),
					     WLZ_RGBA_ALPHA_GET(s[2]));
		      WLZ_RGBA_RGBA_SET(r, m[0], m[1], m[2], m[3]);
		      rP.shp[i] = WLZ_RGBA_MODULUS(r);
		    }
		    break;
		  case WLZ_GREY_UBYTE:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = WLZ_CLAMP(gP.ubp[i] - bin1 - 1, 0, maxBin);
		      s[0] = lut[j];
		      j = WLZ_CLAMP(gP.ubp[i] - bin1,     0, maxBin);
		      s[1] = lut[j];
		      j = WLZ_CLAMP(gP.ubp[i] - bin1 + 1, 0, maxBin);
		      s[2] = lut[j];
		      m[0] = WlzValueDitherI(WLZ_RGBA_RED_GET(s[0]),
		                             WLZ_RGBA_RED_GET(s[1]),
					     WLZ_RGBA_RED_GET(s[2]));
		      m[1] = WlzValueDitherI(WLZ_RGBA_GREEN_GET(s[0]),
		                             WLZ_RGBA_GREEN_GET(s[1]),
					     WLZ_RGBA_GREEN_GET(s[2]));
		      m[2] = WlzValueDitherI(WLZ_RGBA_BLUE_GET(s[0]),
		                             WLZ_RGBA_BLUE_GET(s[1]),
					     WLZ_RGBA_BLUE_GET(s[2]));
		      m[3] = WlzValueDitherI(WLZ_RGBA_ALPHA_GET(s[0]),
		                             WLZ_RGBA_ALPHA_GET(s[1]),
					     WLZ_RGBA_ALPHA_GET(s[2]));
		      WLZ_RGBA_RGBA_SET(r, m[0], m[1], m[2], m[3]);
		      rP.shp[i] = WLZ_RGBA_MODULUS(r);
		    }
		    break;
		  default:
		    errNum = WLZ_ERR_GREY_TYPE;
		    break;
		}
		break;
	      case WLZ_GREY_UBYTE:
		switch(gType)
		{
		  case WLZ_GREY_INT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = WLZ_CLAMP(gP.inp[i] - bin1 - 1, 0, maxBin);
		      s[0] = lut[j];
		      j = WLZ_CLAMP(gP.inp[i] - bin1,     0, maxBin);
		      s[1] = lut[j];
		      j = WLZ_CLAMP(gP.inp[i] - bin1 + 1, 0, maxBin);
		      s[2] = lut[j];
		      m[0] = WlzValueDitherI(WLZ_RGBA_RED_GET(s[0]),
		                             WLZ_RGBA_RED_GET(s[1]),
					     WLZ_RGBA_RED_GET(s[2]));
		      m[1] = WlzValueDitherI(WLZ_RGBA_GREEN_GET(s[0]),
		                             WLZ_RGBA_GREEN_GET(s[1]),
					     WLZ_RGBA_GREEN_GET(s[2]));
		      m[2] = WlzValueDitherI(WLZ_RGBA_BLUE_GET(s[0]),
		                             WLZ_RGBA_BLUE_GET(s[1]),
					     WLZ_RGBA_BLUE_GET(s[2]));
		      m[3] = WlzValueDitherI(WLZ_RGBA_ALPHA_GET(s[0]),
		                             WLZ_RGBA_ALPHA_GET(s[1]),
					     WLZ_RGBA_ALPHA_GET(s[2]));
		      WLZ_RGBA_RGBA_SET(r, m[0], m[1], m[2], m[3]);
		      k = WLZ_RGBA_MODULUS(r);
		      rP.ubp[i] = WLZ_CLAMP(k, 0, 255);
		    }
		    break;
		  case WLZ_GREY_SHORT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = WLZ_CLAMP(gP.shp[i] - bin1 - 1, 0, maxBin);
		      s[0] = lut[j];
		      j = WLZ_CLAMP(gP.shp[i] - bin1,     0, maxBin);
		      s[1] = lut[j];
		      j = WLZ_CLAMP(gP.shp[i] - bin1 + 1, 0, maxBin);
		      s[2] = lut[j];
		      m[0] = WlzValueDitherI(WLZ_RGBA_RED_GET(s[0]),
		                             WLZ_RGBA_RED_GET(s[1]),
					     WLZ_RGBA_RED_GET(s[2]));
		      m[1] = WlzValueDitherI(WLZ_RGBA_GREEN_GET(s[0]),
		                             WLZ_RGBA_GREEN_GET(s[1]),
					     WLZ_RGBA_GREEN_GET(s[2]));
		      m[2] = WlzValueDitherI(WLZ_RGBA_BLUE_GET(s[0]),
		                             WLZ_RGBA_BLUE_GET(s[1]),
					     WLZ_RGBA_BLUE_GET(s[2]));
		      m[3] = WlzValueDitherI(WLZ_RGBA_ALPHA_GET(s[0]),
		                             WLZ_RGBA_ALPHA_GET(s[1]),
					     WLZ_RGBA_ALPHA_GET(s[2]));
		      WLZ_RGBA_RGBA_SET(r, m[0], m[1], m[2], m[3]);
		      k = WLZ_RGBA_MODULUS(r);
		      rP.ubp[i] = WLZ_CLAMP(k, 0, 255);
		    }
		    break;
		  case WLZ_GREY_UBYTE:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = WLZ_CLAMP(gP.ubp[i] - bin1 - 1, 0, maxBin);
		      s[0] = lut[j];
		      j = WLZ_CLAMP(gP.ubp[i] - bin1,     0, maxBin);
		      s[1] = lut[j];
		      j = WLZ_CLAMP(gP.ubp[i] - bin1 + 1, 0, maxBin);
		      s[2] = lut[j];
		      m[0] = WlzValueDitherI(WLZ_RGBA_RED_GET(s[0]),
		                             WLZ_RGBA_RED_GET(s[1]),
					     WLZ_RGBA_RED_GET(s[2]));
		      m[1] = WlzValueDitherI(WLZ_RGBA_GREEN_GET(s[0]),
		                             WLZ_RGBA_GREEN_GET(s[1]),
					     WLZ_RGBA_GREEN_GET(s[2]));
		      m[2] = WlzValueDitherI(WLZ_RGBA_BLUE_GET(s[0]),
		                             WLZ_RGBA_BLUE_GET(s[1]),
					     WLZ_RGBA_BLUE_GET(s[2]));
		      m[3] = WlzValueDitherI(WLZ_RGBA_ALPHA_GET(s[0]),
		                             WLZ_RGBA_ALPHA_GET(s[1]),
					     WLZ_RGBA_ALPHA_GET(s[2]));
		      WLZ_RGBA_RGBA_SET(r, m[0], m[1], m[2], m[3]);
		      k = WLZ_RGBA_MODULUS(r);
		      rP.ubp[i] = WLZ_CLAMP(k, 0, 255);
		    }
		    break;
		  default:
		    errNum = WLZ_ERR_GREY_TYPE;
		    break;
		}
		break;
	      case WLZ_GREY_RGBA:
		switch(gType)
		{
		  case WLZ_GREY_INT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = WLZ_CLAMP(gP.inp[i] - bin1 - 1, 0, maxBin);
		      s[0] = lut[j];
		      j = WLZ_CLAMP(gP.inp[i] - bin1,     0, maxBin);
		      s[1] = lut[j];
		      j = WLZ_CLAMP(gP.inp[i] - bin1 + 1, 0, maxBin);
		      s[2] = lut[j];
		      m[0] = WlzValueDitherI(WLZ_RGBA_RED_GET(s[0]),
		                             WLZ_RGBA_RED_GET(s[1]),
					     WLZ_RGBA_RED_GET(s[2]));
		      m[1] = WlzValueDitherI(WLZ_RGBA_GREEN_GET(s[0]),
		                             WLZ_RGBA_GREEN_GET(s[1]),
					     WLZ_RGBA_GREEN_GET(s[2]));
		      m[2] = WlzValueDitherI(WLZ_RGBA_BLUE_GET(s[0]),
		                             WLZ_RGBA_BLUE_GET(s[1]),
					     WLZ_RGBA_BLUE_GET(s[2]));
		      m[3] = WlzValueDitherI(WLZ_RGBA_ALPHA_GET(s[0]),
		                             WLZ_RGBA_ALPHA_GET(s[1]),
					     WLZ_RGBA_ALPHA_GET(s[2]));
		      WLZ_RGBA_RGBA_SET(rP.rgbp[i], m[0], m[1], m[2], m[3]);
		    }
		    break;
		  case WLZ_GREY_SHORT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = WLZ_CLAMP(gP.shp[i] - bin1 - 1, 0, maxBin);
		      s[0] = lut[j];
		      j = WLZ_CLAMP(gP.shp[i] - bin1,     0, maxBin);
		      s[1] = lut[j];
		      j = WLZ_CLAMP(gP.shp[i] - bin1 + 1, 0, maxBin);
		      s[2] = lut[j];
		      m[0] = WlzValueDitherI(WLZ_RGBA_RED_GET(s[0]),
		                             WLZ_RGBA_RED_GET(s[1]),
					     WLZ_RGBA_RED_GET(s[2]));
		      m[1] = WlzValueDitherI(WLZ_RGBA_GREEN_GET(s[0]),
		                             WLZ_RGBA_GREEN_GET(s[1]),
					     WLZ_RGBA_GREEN_GET(s[2]));
		      m[2] = WlzValueDitherI(WLZ_RGBA_BLUE_GET(s[0]),
		                             WLZ_RGBA_BLUE_GET(s[1]),
					     WLZ_RGBA_BLUE_GET(s[2]));
		      m[3] = WlzValueDitherI(WLZ_RGBA_ALPHA_GET(s[0]),
		                             WLZ_RGBA_ALPHA_GET(s[1]),
					     WLZ_RGBA_ALPHA_GET(s[2]));
		      WLZ_RGBA_RGBA_SET(rP.rgbp[i], m[0], m[1], m[2], m[3]);
		    }
		    break;
		  case WLZ_GREY_UBYTE:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = WLZ_CLAMP(gP.ubp[i] - bin1 - 1, 0, maxBin);
		      s[0] = lut[j];
		      j = WLZ_CLAMP(gP.ubp[i] - bin1,     0, maxBin);
		      s[1] = lut[j];
		      j = WLZ_CLAMP(gP.ubp[i] - bin1 + 1, 0, maxBin);
		      s[2] = lut[j];
		      m[0] = WlzValueDitherI(WLZ_RGBA_RED_GET(s[0]),
		                             WLZ_RGBA_RED_GET(s[1]),
					     WLZ_RGBA_RED_GET(s[2]));
		      m[1] = WlzValueDitherI(WLZ_RGBA_GREEN_GET(s[0]),
		                             WLZ_RGBA_GREEN_GET(s[1]),
					     WLZ_RGBA_GREEN_GET(s[2]));
		      m[2] = WlzValueDitherI(WLZ_RGBA_BLUE_GET(s[0]),
		                             WLZ_RGBA_BLUE_GET(s[1]),
					     WLZ_RGBA_BLUE_GET(s[2]));
		      m[3] = WlzValueDitherI(WLZ_RGBA_ALPHA_GET(s[0]),
		                             WLZ_RGBA_ALPHA_GET(s[1]),
					     WLZ_RGBA_ALPHA_GET(s[2]));
		      WLZ_RGBA_RGBA_SET(rP.rgbp[i], m[0], m[1], m[2], m[3]);
		    }
		    break;
		  default:
		    errNum = WLZ_ERR_GREY_TYPE;
		    break;
		}
		break;
	      default:
		errNum = WLZ_ERR_GREY_TYPE;
		break;
	    }
	  }
	  else /* dither == 0 */
	  {
	    int		i,
	    		j,
	    		k;

	    switch(rType)
	    {
	      case WLZ_GREY_INT:
		switch(gType)
		{
		  case WLZ_GREY_INT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.inp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      rP.inp[i] = WLZ_RGBA_MODULUS(lut[j]);
		    }
		    break;
		  case WLZ_GREY_SHORT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.shp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      rP.inp[i] = WLZ_RGBA_MODULUS(lut[j]);
		    }
		    break;
		  case WLZ_GREY_UBYTE:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.ubp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      rP.inp[i] = WLZ_RGBA_MODULUS(lut[j]);
		    }
		    break;
		  default:
		    errNum = WLZ_ERR_GREY_TYPE;
		    break;
		}
		break;
	      case WLZ_GREY_SHORT:
		switch(gType)
		{
		  case WLZ_GREY_INT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.inp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      rP.shp[i] = WLZ_RGBA_MODULUS(lut[j]);
		    }
		    break;
		  case WLZ_GREY_SHORT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.shp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      rP.shp[i] = WLZ_RGBA_MODULUS(lut[j]);
		    }
		    break;
		  case WLZ_GREY_UBYTE:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.ubp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      rP.shp[i] = WLZ_RGBA_MODULUS(lut[j]);
		    }
		    break;
		  default:
		    errNum = WLZ_ERR_GREY_TYPE;
		    break;
		}
		break;
	      case WLZ_GREY_UBYTE:
		switch(gType)
		{
		  case WLZ_GREY_INT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.inp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      k = WLZ_RGBA_MODULUS(lut[j]);
		      rP.ubp[i] = WLZ_CLAMP(k, 0, 255);
		    }
		    break;
		  case WLZ_GREY_SHORT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.shp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      k = WLZ_RGBA_MODULUS(lut[j]);
		      rP.ubp[i] = WLZ_CLAMP(k, 0, 255);
		    }
		    break;
		  case WLZ_GREY_UBYTE:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.ubp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      k = WLZ_RGBA_MODULUS(lut[j]);
		      rP.ubp[i] = WLZ_CLAMP(k, 0, 255);
		    }
		    break;
		  default:
		    errNum = WLZ_ERR_GREY_TYPE;
		    break;
		}
		break;
	      case WLZ_GREY_RGBA:
		switch(gType)
		{
		  case WLZ_GREY_INT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.inp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      rP.rgbp[i] = lut[j];
		    }
		    break;
		  case WLZ_GREY_SHORT:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.shp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      rP.rgbp[i] = lut[j];
		    }
		    break;
		  case WLZ_GREY_UBYTE:
		    for(i = 0; i < nVal; ++i)
		    {
		      j = gP.ubp[i] - bin1;
		      j = WLZ_CLAMP(j, 0, maxBin);
		      rP.rgbp[i] = lut[j];
		    }
		    break;
		  default:
		    errNum = WLZ_ERR_GREY_TYPE;
		    break;
		}
		break;
	      default:
		errNum = WLZ_ERR_GREY_TYPE;
		break;
	    }
	  }
	}
	break;
      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  return(errNum);
}
