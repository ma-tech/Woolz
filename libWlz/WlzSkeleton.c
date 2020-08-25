#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzSkeleton_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzSkeleton.c
* \author       Jim Piper, Bill Hill, Richard Baldock
* \date         March 1997
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
* \brief	Performs a proper interval-domain skeletonisation
* 		Hilditch's  method.
* \ingroup	WlzDomainOps
*/

#include <stdio.h>
#include <Wlz.h>

/* a line of intervals plus current position information */
typedef struct _WlzSkIntvLn
{
  WlzInterval	*intp;
  int		intcount;
  int		ijntpos;
  int		ijleft;
  int		ijright;
  int		prevright;
} WlzSkIntvLn;

typedef struct _WlzExtIntv
{
  int		jjleft;
  int		jjright;
  int		nextleft;
} WlzExtIntv;

static WlzErrorNum 	WlzSkStrip8(WlzObject *, WlzDomain,
			            WlzInterval *, int, int *, int),
			WlzSkStrip4(WlzObject *, WlzObject *, WlzObject *,
			            WlzInterval *, int, int *, int);
static WlzObject 	*WlzSkeleton2D(WlzObject *, int, WlzConnectType,
			               WlzErrorNum *),
		 	*WlzSkeleton3D(WlzObject *, int, WlzConnectType,
		 		       WlzErrorNum *);

/*! 
* \return       Skeleton object, NULL on error.
* \ingroup      WlzDomainOps
* \brief        Computes the skeleton of the given object using Hilditch's
*		 method. See detail.
*
* \par 	 Skeleton algorithm
Performs a proper interval-domain skeletonisation by Hilditch's  method.
<p>This algorithm iterates a two-stage process with three interval
      domain objects 
  as follows: </p>
<p> Stage 1 - remove points according to algorithm: </p>
  <ul>
    <li> (i) is current input object.</li>
    <li> (ii) is current subset of (i) where a point deletion is
	    possible. </li>
    <li> (iii) is the set of points deleted during this pass</li>
  </ul>
<p> Stage 2 - tidy up for next pass:</p>
<ol>
  <li>check for termination - object (iii) above empty</li>
  <li> construct new, thinned input object (i) by subtracting object
	(iii) from 
    old object (i).</li>
  <li> construct new object (ii) by taking dilation of (iii) and
	intersecting 
    with new (i). </li>
</ol>
<p>The method can be further improved by reducing the number of points
      scanned 
  in the first pass to edge points. These are extracted as the
      difference between 
  the original object and its erosion. Implement Hilditch's tests 1-6
      by table 
  look up (48 words, 32 bits). Neighbour order is as follows:</p>
<table width="301" border="1" cellpadding="2">
  <tr>
    <td width="77"><table width="0%" border="0" cellpadding="2">
        <tr>
          <td>3</td>
          <td>2*</td>
          <td>1</td>
        </tr>
        <tr>
          <td>4*</td>
          <td>&nbsp;</td>
          <td>0=8</td>
        </tr>
        <tr>
          <td>5</td>
          <td>6</td>
          <td>7</td>
        </tr>
      </table></td>
    <td width="184"><p>asterisk marks neighbourswhere &quot;removed
	      this pass&quot; (RTP) 
        is relevant to algorithm </p>
      </td>
  </tr>
</table>
<p>Use neighbours 0 - 4 and RTP 2 to generate an address in range 0 -
      47 (neighbour 
  2 has three effective values: not in object previously, in object
      now, removed 
  this pass), bit values as follows:</p>
<table width="0%" border="0" cellpadding="2">
  <tr>
    <td>8</td>
    <td>16/32*</td>
    <td>4</td>
  </tr>
  <tr>
    <td>2</td>
    <td>-</td>
    <td>1</td>
  </tr>
  <tr>
    <td>-</td>
    <td>-</td>
    <td>-</td>
  </tr>
</table>
<p>or without 4 connected points: </p>
<table width="0%" border="0" cellpadding="2">
  <tr>
    <td>2</td>
    <td>16/32*</td>
    <td>1</td>
  </tr>
  <tr>
    <td>8</td>
    <td>-</td>
    <td>4</td>
  </tr>
  <tr>
    <td>-</td>
    <td>-</td>
    <td>-</td>
  </tr>
</table>
<p>Use neighbours 5-7, RTP 4, and the smoothing criterion to generate
      a 5 bit 
  bit-address, look at corresponding bit in addressed word in look-up
      table. If 
  zero, retain point. Address values as follows: </p>
<table width="0%" border="0" cellpadding="2">
  <tr>
    <td>-</td>
    <td>-</td>
    <td>-</td>
  </tr>
  <tr>
    <td>8*</td>
    <td>-</td>
    <td>-</td>
  </tr>
  <tr>
    <td>2</td>
    <td>4</td>
    <td>1</td>
  </tr>
</table>
<p>Thus lookup table is 48 words long, 32 bit unsigned values. </p>
* \par Notes
Look-up table constructed by separate program	
"newskelsetup.c".				
We must intersect potDelObj and skObj on each iteration
or alternatively explicitly look at each element of 
potDelObj and see if a member of skObj otherwise we 
may decide to delete points outside skObj and then muck
up points within skObj.
* \par      Source:
*                WlzSkeleton.c
* \param    srcObj		Input object.
* \param    smoothpasses	Number of smoothing passes to be applied.
* \param    minCon		Minimum connectivity required.
* \param    dstErr		Error return.
*/
WlzObject 	*WlzSkeleton(WlzObject *srcObj, int smoothpasses,
			     WlzConnectType minCon,
			     WlzErrorNum *dstErr)
{
  WlzObject	*skObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(srcObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(srcObj->type)
    {
      case WLZ_EMPTY_OBJ:
	skObj = WlzMakeEmpty(&errNum);
        break;
      case WLZ_2D_DOMAINOBJ:
	skObj = WlzSkeleton2D(srcObj, smoothpasses, minCon, &errNum);
        break;
      case WLZ_3D_DOMAINOBJ:
        skObj = WlzSkeleton3D(srcObj, smoothpasses, minCon, &errNum);
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
  return(skObj);
}

static WlzObject *WlzSkeleton3D(WlzObject *srcObj, int smoothpasses,
				WlzConnectType minCon,
				WlzErrorNum *dstErr)
{
  WlzObject	*skObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((minCon != WLZ_6_CONNECTED) && (minCon != WLZ_18_CONNECTED) &&
     (minCon != WLZ_26_CONNECTED))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Code not written yet! */
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(skObj);
}

static WlzObject *WlzSkeleton2D(WlzObject *srcObj, int smoothpasses,
			        WlzConnectType minCon,
			        WlzErrorNum *dstErr)
{
  int		delArea = 0,
  		maxItv,
		passCnt,
		srcItvCount;
  WlzObject	*tObj,
  		*delObj = NULL,
  		*potDelObj = NULL,
		*skObj = NULL;
  WlzDomain	delDom,
  		dumDom,
  		altDelDom;
  WlzValues	dumVal;
  WlzInterval	*altItvBase = NULL,
  		*delItvBase = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  AlcErrno	alcErr = ALC_ER_NONE;
  const int	minItv = 1024,
  		sclItv = 16;

  dumDom.core = NULL;
  delDom.core = NULL;
  altDelDom.core = NULL;
  dumVal.core = NULL;
  if((minCon != WLZ_8_CONNECTED) && (minCon != WLZ_4_CONNECTED))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((srcObj->domain.core->type != WLZ_INTERVALDOMAIN_INTVL) &&
          (srcObj->domain.core->type != WLZ_INTERVALDOMAIN_RECT))
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    srcItvCount = WlzIntervalCount(srcObj->domain.i, &errNum);
    if((srcItvCount <= 0) && (errNum == WLZ_ERR_NONE))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }

  /* check for rectangular domain which has no intervalline pointer */
  if(errNum == WLZ_ERR_NONE){
    if(srcObj->domain.core->type == WLZ_INTERVALDOMAIN_RECT){
      if((dumDom.i = WlzNewIDomain(WLZ_INTERVALDOMAIN_INTVL,
				   srcObj->domain.i, &errNum)) != NULL){
	WlzFreeDomain(srcObj->domain);
	srcObj->domain = WlzAssignDomain(dumDom, &errNum);
	dumDom.core = NULL;
      }
    }
  }

  if(errNum == WLZ_ERR_NONE)
  {
    /* Copy the input object, and make plenty of space for intermediate
     * intervals in the deleted-this-pass object. */
    maxItv = minItv + (sclItv * srcItvCount);
    if(minCon == WLZ_8_CONNECTED)
    {
      delDom.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
      				       srcObj->domain.i->line1,
				       srcObj->domain.i->lastln,
				       srcObj->domain.i->kol1,
				       srcObj->domain.i->lastkl, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
        (void )WlzAssignDomain(delDom, NULL);
        if((delItvBase = (WlzInterval *)
			 AlcMalloc(maxItv * sizeof(WlzInterval))) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  delDom.i->freeptr = AlcFreeStackPush(delDom.i->freeptr,
					       (void *)delItvBase,
					       &alcErr);
	  if(alcErr != ALC_ER_NONE)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
        altDelDom.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
					    srcObj->domain.i->line1,
					    srcObj->domain.i->lastln,
					    srcObj->domain.i->kol1,
					    srcObj->domain.i->lastkl, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        (void )WlzAssignDomain(altDelDom, NULL);
        if((altItvBase = (WlzInterval *)
			 AlcMalloc(maxItv * sizeof(WlzInterval))) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  altDelDom.i->freeptr = AlcFreeStackPush(altDelDom.i->freeptr,
						  (void *)altItvBase,
						  &alcErr);
	  if(alcErr != ALC_ER_NONE)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
      }
      skObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, dumDom, dumVal, NULL, srcObj,
	 		  &errNum);
    }
    else /* minCon == WLZ_4_CONNECTED */
    {
      delDom.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
      				       srcObj->domain.i->line1,
				       srcObj->domain.i->lastln,
				       srcObj->domain.i->kol1,
				       srcObj->domain.i->lastkl, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	(void )WlzAssignDomain(delDom, NULL);
        delObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, delDom, dumVal, NULL, srcObj,
			     &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        if((delItvBase = (WlzInterval *)
			 AlcMalloc(maxItv * sizeof(WlzInterval))) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  delDom.i->freeptr = AlcFreeStackPush(delDom.i->freeptr,
					       (void *)delItvBase,
					       &alcErr);
	  if(alcErr != ALC_ER_NONE)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(minCon == WLZ_8_CONNECTED)
    {
      /* First pass - output thinned object in delObj. */
      errNum = WlzSkStrip8(srcObj, delDom, delItvBase, maxItv, &delArea,
	  		   smoothpasses--);
      /* Subsequent passes - use output of previous pass as input,
       * swap between the two temporary objects. */
      while((errNum == WLZ_ERR_NONE) && (delArea > 0))
      {
	skObj->domain = delDom;
	delDom = altDelDom;
	altDelDom = skObj->domain;
	errNum = WlzSkStrip8(skObj, delDom, delItvBase, maxItv,
	    	             &delArea, smoothpasses--);
      }
      /* Now need to standardise the interval domain, copy for
       * output object, and free workspace. */
      
      skObj->domain.i = WlzNewIDomain(WLZ_INTERVALDOMAIN_INTVL,
	  			      delDom.i, &errNum);
      (void )WlzFreeDomain(delDom);
      (void )WlzFreeDomain(altDelDom);
    }
    else /* minCon == WLZ_4_CONNECTED */
    {
      /*
       * First pass - the potentially deletable object is input object.
       */
      potDelObj = srcObj;
      skObj = srcObj;
      errNum = WlzSkStrip4(skObj, potDelObj, delObj, delItvBase, maxItv,
	  	           &delArea, smoothpasses--);
      passCnt = 0;
      /* Subsequent passes. */
      while((errNum == WLZ_ERR_NONE) && (delArea > 0))
      {
	tObj = skObj;
	skObj = WlzDiffDomain(tObj, delObj, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  if(++passCnt > 1)
	  {
	    (void )WlzFreeObj(tObj);
	  }
	  potDelObj = skObj;
	  errNum = WlzSkStrip4(skObj, potDelObj, delObj, delItvBase, maxItv,
		               &delArea, smoothpasses--);
	}
      }
      skObj->values.core = NULL;
    }
  }
  return(skObj);
}

static WlzErrorNum WlzSkStrip4(WlzObject *skObj, WlzObject *potDelObj,
    			       WlzObject *delObj, WlzInterval *delItvBase,
			       int itvSpace, int *delArea, int smPass)
{ 
  int		b,
  		delflipflop,
  		k,
		k1,
		kl,
  		km1,
		l,
  		lastrem,
		pdl1,
		skl1,
		skll,
  		smBits,
		lintc = 0;
  unsigned int	w;
  WlzDomain	delDom,
  		skDom;
  WlzInterval	*jntp,
  		*intp = NULL;
  WlzIntervalLine *intl,
  		 *lintl;
  WlzIntervalWSpace iWSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzSkIntvLn	dm1lint = {0},
  		llint = {0},
  		lm1lint = {0},
		lp1lint = {0};
  const unsigned int wlzSkLut4[48] =
  {
000075400340, 000074600362, 036475372364, 036074170360, 
000000200000, 000074600363, 000000000000, 036074170360, 
000000200000, 000000000000, 000175200765, 000074000360, 
000000000000, 000000000000, 000000000000, 000074000360, 
000000200000, 000074600363, 000175200765, 000303601417, 
000000200001, 000074600363, 000175200765, 000303601417, 
000000200001, 000074600363, 000175200765, 000303601417, 
000000200001, 000074600363, 000175200765, 000303601417, 
000000000000, 000074600363, 000175200765, 000000000000, 
000000200001, 000074600363, 000000000000, 000000000000, 
000000200001, 000000000000, 000175200765, 000000000000, 
000000000000, 000000000000, 000000000000, 000000000000, 

/*    000075400340, 000074600362, 036475372364, 036074170360,
    000000200000, 000074600363, 000000000000, 036074170360,
    000000200000, 000000000000, 000175200765, 000074000360,
    000000000000, 000000000000, 000000000000, 000074000360,
    000000200000, 000074600363, 000175200765, 000303601417,
    000000200001, 000074600363, 000175200765, 000303601417,
    000000200001, 000074600363, 000175200765, 000303601417,
    000000200001, 000074600363, 000175200765, 000303601417,
    000000000000, 000074600363, 000175200765, 000000000000,
    000000200001, 000074600363, 000000000000, 000000000000,
    000000200001, 000000000000, 000175200765, 000000000000,
    000000000000, 000000000000, 000000000000, 000000000000*/
  };

  smBits = (smPass >=0)? 16: 0;
  *delArea = 0;
  skDom = skObj->domain;
  skl1 = skDom.i->line1;
  skll = skDom.i->lastln;
  pdl1 = potDelObj->domain.i->line1;
  /* Set up interval table for delObj. */
  delDom = delObj->domain;
  jntp = delItvBase;
  intl = delDom.i->intvlines;
  l = skll - skl1;
  for(k=0; k<=l; ++k)
  {
    (intl++)->nintvs = 0;
  }
  /* Start scanning potDel. */
  errNum = WlzInitRasterScan(potDelObj, &iWSp, WLZ_RASTERDIR_ILIC);
  while((errNum == WLZ_ERR_NONE) &&
        ((errNum = WlzNextInterval(&iWSp)) == WLZ_ERR_NONE))
  {
    /* If beginning of line, much to set up. */
    if(iWSp.nwlpos != 0)
    {
      l = iWSp.linpos - skl1;
      /* We know current line has intervals. */
      lintl = skDom.i->intvlines + l;
      if(l > 0)
      {
	intl = lintl - 1;
	if((lm1lint.intcount = (intl)->nintvs) > 0)
	{
	  lm1lint.intp = intp = intl->intvs;
	  lm1lint.ijleft = intp->ileft;
	  lm1lint.ijright = intp->iright;
	  lm1lint.ijntpos = 0;
	}
      }
      else
      {
	lm1lint.intcount = 0;
      }
      if(iWSp.linpos > pdl1)
      {
	intl = delDom.i->intvlines + l - 1;
	if((dm1lint.intcount = intl->nintvs) > 0)
	{
	  dm1lint.intp = intp = intl->intvs;
	  dm1lint.ijleft = intp->ileft;
	  dm1lint.ijright = intp->iright;
	  dm1lint.ijntpos = 0;
	}
      }
      else
      {
	dm1lint.intcount = 0;
      }
      if(iWSp.linpos < skll)
      {
	if((lp1lint.intcount = (++lintl)->nintvs) > 0)
	{
	  lp1lint.intp = intp = lintl->intvs;
	  lp1lint.ijleft = intp->ileft;
	  lp1lint.ijright = intp->iright;
	  lp1lint.ijntpos = 0;
	}
      }
      else
      {
	lp1lint.intcount = 0;
      }
      lintc = 0;
      intp = jntp;
    }
    /* Compute deletable points in this interval. */
    lastrem = 0;
    delflipflop = 0;
    k1 = iWSp.lftpos - skDom.i->kol1; 
    kl = iWSp.rgtpos - skDom.i->kol1;
    for(k=k1; k<=kl; ++k)
    {
      km1 = k-1;
      /* Here "k" refers to the left-most position in the run of three. */
      if(lm1lint.intcount <= 0)
      {
	w = 0; 
	goto exitP3a; 
      } 
      while(km1 > lm1lint.ijright)
      {
	if(++(lm1lint.ijntpos) >= lm1lint.intcount)
	{
	  lm1lint.intcount = -1; 
	  w = 0; 
	  goto exitP3a; 
	} 
	else {			  
	  /* Increment interval pointer. */
	  lm1lint.ijleft = (++lm1lint.intp)->ileft; 
	  lm1lint.ijright = lm1lint.intp->iright; 
	} 
      } 
      if(km1+2 < lm1lint.ijleft)
      {
	/* 0 0 0 */
	w = 0; 
	goto exitP3a; 
      }
      else if(km1 < lm1lint.ijleft)
      {
	if(km1+2 == lm1lint.ijleft)
	{
	  /* 0 0 1 */
	  w = 1; 
	  goto exitP3a; 
	}
	else
	{ 
	  if(km1+2 <= lm1lint.ijright)
	  {
	    /* 0 1 1 */
	    w = 5; 
	    goto exitP3a; 
	  }
	  else
	  {			  
	    /* 0 1 0 */
	    w = 4; 
	    goto exitP3a; 
	  }
	} 
      }
      else
      { 
	if(km1+2 <= lm1lint.ijright)
	{
	  /* 1 1 1 */
	  w = 7; 
	  goto exitP3a; 
	}
	else if(km1+1 == lm1lint.ijright)
	{
	  /* 1 1 0 */
	  w = 6; 
	  goto exitP3a; 
	}
	else if((lm1lint.ijntpos < lm1lint.intcount) && 
	        (km1+2 == ((WlzExtIntv *)(lm1lint.intp))->nextleft))
	{
	  /* 1 0 1 */
	  w = 3;				  
	  goto exitP3a; 
	}
	else
	{				  
	  /* 1 0 0 */
	  w = 2; 
	} 
      } 
exitP3a:
      w <<= 2;
      /* Here "km1" refers to the left-most position in the run of three. */
      if(llint.intcount <= 0)
      {
	goto exitP2; 
      } 
      while(km1 > llint.ijright)
      {
	if(++llint.ijntpos >= llint.intcount)
	{
	  llint.intcount = -1; 
	  goto exitP2; 
	} 
	else
	{			  
	  /* Increment interval pointer. */
	  llint.ijleft = (++llint.intp)->ileft; 
	  llint.ijright = llint.intp->iright; 
	} 
      } 
      if(km1+2 < llint.ijleft)
      {
	/* 0 0 0 */
	goto exitP2; 
      }
      else if(km1 < llint.ijleft)
      {
	if(km1+2 <= llint.ijright)
	{
	  /* 0 0 1, 0 1 1 */
	  w += 1; 
	  goto exitP2; 
	}
	else
	{				  
	  /* 0 1 0 */
	  goto exitP2; 
	} 
      }
      else
      { 
	if(km1+2 <= llint.ijright)
	{
	  /* 1 1 1 */
	  w += 3; 
	  goto exitP2; 
	}
	else if(km1+1 == llint.ijright)
	{
	  /* 1 1 0 */
	  w += 2; 
	  goto exitP2; 
	}
	else if((llint.ijntpos < llint.intcount) && 
	        (km1+2 == ((WlzExtIntv *)(llint.intp))->nextleft))
	{
	  /* 1 0 1 */
	  w += 3;				  
	  goto exitP2; 
	}
	else
	{
	  /* 1 0 0 */
	  w += 2; 
	}
      } 
exitP2: 
      if((dm1lint.intcount <= 0) || (k < dm1lint.ijleft))
      {
	goto exitP1; 
      }
      while(k > dm1lint.ijright)
      {
	if(++dm1lint.ijntpos >= dm1lint.intcount)
	{
	  dm1lint.intcount = -1; 
	  goto exitP1; 
	} 
	else { 
	  dm1lint.prevright = dm1lint.ijright; 
	  dm1lint.ijleft = (++dm1lint.intp)->ileft; 
	  dm1lint.ijright = dm1lint.intp->iright; 
	  if(k < dm1lint.ijleft)
	  {
	    goto exitP1; 
	  } 
	} 
      } 
      w += 16; 
exitP1: 
      w = wlzSkLut4[w];
      if(w != 0)
      {
	/* Here "k" refers to the left-most position in the run of three. */
	if(lp1lint.intcount <= 0)
	{
	  b = 0; 
	  goto exitP3b; 
	} 
	while(km1 > lp1lint.ijright)
	{
	  if(++(lp1lint.ijntpos) >= lp1lint.intcount)
	  {
	    lp1lint.intcount = -1; 
	    b = 0; 
	    goto exitP3b; 
	  } 
	  else
	  {			  
	    /* Increment interval pointer. */
	    lp1lint.ijleft = (++lp1lint.intp)->ileft; 
	    lp1lint.ijright = lp1lint.intp->iright; 
	  } 
	} 
	if(km1+2 < lp1lint.ijleft)
	{
	  /* 0 0 0 */
	  b = 0; 
	  goto exitP3b; 
	}
	else if(km1 < lp1lint.ijleft)
	{
	  if(km1+2 == lp1lint.ijleft)
	  {
	    /* 0 0 1 */
	    b = 1; 
	    goto exitP3b; 
	  }
	  else
	  { 
	    if(km1+2 <= lp1lint.ijright)
	    {
	      /* 0 1 1 */
	      b = 5; 
	      goto exitP3b; 
	    }
	    else
	    {			  
	      /* 0 1 0 */
	      b = 4; 
	      goto exitP3b; 
	    }
	  } 
	}
	else
	{ 
	  if(km1+2 <= lp1lint.ijright)
	  {
	    /* 1 1 1 */
	    b = 7; 
	    goto exitP3b; 
	  }
	  else if(km1+1 == lp1lint.ijright)
	  {
	    /* 1 1 0 */
	    b = 6; 
	    goto exitP3b; 
	  }
	  else if((lp1lint.ijntpos < lp1lint.intcount) && 
	          (km1+2 == ((WlzExtIntv *)(lp1lint.intp))->nextleft))
	  {
	    /* 1 0 1 */
	    b = 3;				  
	    goto exitP3b; 
	  }
	  else
	  {				  
	    /* 1 0 0 */
	    b = 2; 
	  } 
	} 
exitP3b:
	b += (smBits + lastrem);
	w &= (1 << b);
      }
      if(w != 0)
      {
	(*delArea)++;
	lastrem = 8;	 
	if(delflipflop == 0)
	{
	  lintc++;
	  delflipflop = 1;
	  intp->ileft = k;
	}
	intp->iright = k;
      }
      else
      {			 
	lastrem = 0;	 
	if(delflipflop == 1)
	{
	  delflipflop = 0;
	  intp++;
	}
      }
    }
    if(delflipflop == 1)
    {
      intp++;
    }
    if(iWSp.intrmn == 0)
    {
      if(itvSpace > lintc)
      {
	errNum = WlzMakeInterval(iWSp.linpos, delDom.i, lintc, jntp);
	jntp = intp;
	itvSpace -= lintc;
      }
      else
      {
        errNum = WLZ_ERR_DOMAIN_DATA;
      }
    }
  }
  if(errNum == WLZ_ERR_EOO)		/* Reset error from end of intervals */
  {
    errNum = WLZ_ERR_NONE;
  }
  return(errNum);
}

static WlzErrorNum WlzSkStrip8(WlzObject *skObj, WlzDomain delDom,
			       WlzInterval *delItvBase,
			       int itvSpace, int *delArea, int smPass)
{ 
  int		b,
		delflipflop,
    		k,
		l,
		lastrem,
		k1,
		kl,
		pdl1,
		rtp,
		smBits,
		skl1,
		skll,
		lintc = 0;
  unsigned int	w;
  WlzDomain	skDom;
  WlzInterval	*intp = NULL,
  		*jntp = NULL;
  WlzIntervalLine *intl,
  		 *lintl;
  WlzSkIntvLn	dm1lint = {0},
  		lm1lint = {0},
		lp1lint = {0};
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzIntervalWSpace iWSp;
  const unsigned int wlzSkLut8[48] =
  {
    000075400340, 000000200000, 000000200000, 000000000000, 
    000074600362, 000074600363, 000000000000, 000000000000, 
    036475372364, 000000000000, 000175200765, 000000000000, 
    036074170360, 036074170360, 000074000360, 000074000360, 
    000000200000, 000000200001, 000000200001, 000000200001, 
    000074600363, 000074600363, 000074600363, 000074600363, 
    000175200765, 000175200765, 000175200765, 000175200765, 
    000303601417, 000303601417, 000303601417, 000303601417, 
    000000000000, 000000200001, 000000200001, 000000000000, 
    000074600363, 000074600363, 000000000000, 000000000000, 
    000175200765, 000000000000, 000175200765, 000000000000, 
    000000000000, 000000000000, 000000000000, 000000000000, 
  };

  smBits = (smPass >=0)? 16: 0;
  *delArea = 0;
  skDom = skObj->domain;
  skl1 = skDom.i->line1;
  skll = skDom.i->lastln;
  pdl1 = skl1;
  jntp = delItvBase;
  intl = delDom.i->intvlines;
  l = skll - skl1;
  for(k=0; k<=l; ++k)
  {
    (intl++)->nintvs = 0;
  }

  /* Start scanning potDel. */
  errNum = WlzInitRasterScan(skObj, &iWSp, WLZ_RASTERDIR_ILIC);
  while((errNum == WLZ_ERR_NONE) &&
        ((errNum = WlzNextInterval(&iWSp)) == WLZ_ERR_NONE))
  {
    /* If beginning of line, much to set up. */
    if(iWSp.nwlpos != 0)
    {
      /* l is the number of lines above line1
	 lintl is the current line of intervals */
      l = iWSp.linpos - skl1;
      lintl = skDom.i->intvlines + l;

      /* if not line1 then set up line/interval struct for the
	 preceding line */
      if(l > 0 )
      {
	intl = lintl - 1;
	if((lm1lint.intcount = (intl)->nintvs) > 0)
	{
	  lm1lint.intp = intp = intl->intvs;
	  lm1lint.ijleft = intp->ileft;
	  lm1lint.prevright =
	    lm1lint.ijright = intp->iright;
	}
      }
      else
      {
	lm1lint.intcount = 0;
      }

      /* the preceding line in the thinned object
	 - unless pdl1 has been incremented  - it isn't */
      if(iWSp.linpos > pdl1)
      {
	intl = delDom.i->intvlines + l - 1;
	if((dm1lint.intcount = intl->nintvs) > 0)
	{
	  dm1lint.intp = intp = intl->intvs;
	  dm1lint.ijleft = intp->ileft;
	  dm1lint.prevright =
	    dm1lint.ijright = intp->iright;
	}
      }
      else
      {
	dm1lint.intcount = 0;
      }

      /* if not the last line increment lintl and set up
	 structure for the following line */
      if(iWSp.linpos < skll)
      {
	if((lp1lint.intcount = (++lintl)->nintvs) > 0)
	{
	  lp1lint.intp = intp = lintl->intvs;
	  lp1lint.ijleft = intp->ileft;
	  lp1lint.prevright =
	    lp1lint.ijright = intp->iright;
	}
      }
      else
      {
	lp1lint.intcount = 0;
      }

      /* local interval count and reset interval pointer */
      lintc = 0;
      intp = jntp;
    }

    /* Compute deletable points in this interval. */
    lastrem = 0;
    delflipflop = 0;
    k1 = iWSp.lftpos - skDom.i->kol1; 
    kl = iWSp.rgtpos - skDom.i->kol1;
    for(k=k1; k<=kl; ++k)
    {
      /* determine the 4-connected neighbours - up, left, right
	 these determine the word address or LUT entry */
      if(lm1lint.intcount <= 0 || k < lm1lint.ijleft)
      {
	w = 0; 
	goto exitP1a; 
      } 
      while(k > lm1lint.ijright)
      {
	if(--lm1lint.intcount <= 0)
	{
	  lm1lint.intcount = -1;
	  w = 0; 
	  goto exitP1a; 
	} 
	else
	{ 
	  lm1lint.prevright = lm1lint.ijright; 
	  lm1lint.ijleft = (++lm1lint.intp)->ileft; 
	  lm1lint.ijright = lm1lint.intp->iright; 
	  if(k < lm1lint.ijleft)
	  {
	    w = 0; 
	    goto exitP1a; 
	  } 
	} 
      } 
      w = 16;

      /* at this point check for 4-connected left & right */
exitP1a: 
      if(k > k1)
      {
	w += 8;
      }
      if(k < kl)
      {
	w += 4;
      }

      /* determine the bit address - depends on the next line
	 and left point */
      if((lp1lint.intcount <= 0) || (k < lp1lint.ijleft))
      {
	b = 0; 
	goto exitP1b; 
      }
      while(k > lp1lint.ijright)
      {
	if(--lp1lint.intcount <= 0)
	{
	  lp1lint.intcount = -1;
	  b = 0; 
	  goto exitP1b; 
	} 
	else
	{ 
	  lp1lint.prevright = lp1lint.ijright; 
	  lp1lint.ijleft = (++lp1lint.intp)->ileft; 
	  lp1lint.ijright = lp1lint.intp->iright; 
	  if(k < lp1lint.ijleft)
	  {
	    b = 0; 
	    goto exitP1b; 
	  } 
	} 
      } 
      b = 16; 
exitP1b:
      if(b > 0)						 /* Actually b == 16 */
      {
	/* Not a 4-connected boundary point? - quick cut out */
	if(w == 28)
	{
	  /* Retain point. */
	  lastrem = 0;	 				   /* For next point */
	  if(delflipflop == 0)
	  {
	    lintc++;
	    delflipflop = 1;
	    intp->ileft = k;
	  }
	  /* Skip over further interior points. */
	  k = WLZ_MIN((kl - 1), lm1lint.ijright);
	  k = WLZ_MIN(k, lp1lint.ijright);
	  intp->iright = k;
	  continue;	 
	}
	b = 4;
      }
      /* If neighbour 2 present, see if it still present in
       * thinned object "delObj". */
      if(w >= 16)
      {
	if((dm1lint.intcount <= 0) || (k < dm1lint.ijleft))
	{
	  rtp = 0;
	  goto exitP1c; 
	} 
	while(k > dm1lint.ijright)
	{
	  if(--dm1lint.intcount <= 0)
	  {
	    dm1lint.intcount = -1;
	    rtp = 0; 
	    goto exitP1c; 
	  } 
	  else { 
	    dm1lint.prevright = dm1lint.ijright; 
	    dm1lint.ijleft = (++dm1lint.intp)->ileft; 
	    dm1lint.ijright = dm1lint.intp->iright; 
	    if(k < dm1lint.ijleft)
	    {
	      rtp = 0; 
	      goto exitP1c; 
	    } 
	  } 
	} 
	rtp = 16; 
exitP1c:
	w += (16 - rtp);
      }
      /* Here "k" refers to the middle position in the run of three.
       * Interval stepping has already been performed.
       * First check intervals. */
      if(lm1lint.intcount != 0)
      {
	/* Past end of last interval in line? */
	if(lm1lint.intcount < 0)
	{
	  if(k-1 == lm1lint.ijright) 
	  {
	    w += 2; 
	  }
	}
	else
	{ 
	  /* Otherwise not yet past last interval, check cases. */
	  if(k > lm1lint.ijleft || k-1 == lm1lint.prevright) 
	  {
	    w += 2; 
	  }
	  if(k < lm1lint.ijright && k+1 >= lm1lint.ijleft) 
	  {
	    w++; 
	  }
	} 
      } 
      w = wlzSkLut8[w];
      if(w != 0)
      {
	{   
	  /* Here "k" refers to the middle position in the run of three.
	   * Interval checking has already been performed.
	   * First check if any intervals in the line. */
	  if(lp1lint.intcount != 0)
	  {
	    /* Past end of last interbal in line? */
	    if(lp1lint.intcount < 0)
	    {
	      if(k-1 == lp1lint.ijright) 
	      {
		b += 2; 
	      }
	    }
	    else
	    { 
	      /* Otherwise not yet past last interval, check cases. */
	      if((k > lp1lint.ijleft) || (k-1 == lp1lint.prevright) )
	      {
		b += 2; 
	      }
	      if((k < lp1lint.ijright) && (k+1 >= lp1lint.ijleft))
	      {
		b++; 
	      }
	    } 
	  } 
	}
	b += (smBits + lastrem);
	w &= (1 << b);
      }
      if(w != 0)
      {
	(*delArea)++;
	lastrem = 8;	 				   /* For next point */
	if(delflipflop == 1)
	{
	  delflipflop = 0;
	  intp++;
	}
      }
      else
      {			 
	lastrem = 0;	 				   /* For next point */
	if(delflipflop == 0)
	{
	  lintc++;
	  delflipflop = 1;
	  intp->ileft = k;
	}
	intp->iright = k;
      }
    }
    /* Clear up an open interval end. */
    if(delflipflop == 1)
    {
      intp++;
    }
    /* Tidy up at end of line. */
    if(iWSp.intrmn == 0)
    {
      if(itvSpace > lintc)
      {
	errNum = WlzMakeInterval(iWSp.linpos, delDom.i, lintc, jntp);
	jntp = intp;
	itvSpace -= lintc;
      }
      else
      {
        errNum = WLZ_ERR_DOMAIN_DATA;
      }
    }
  }
  if(errNum == WLZ_ERR_EOO)		/* Reset error from end of intervals */
  {
    errNum = WLZ_ERR_NONE;
  }
  return(errNum);
}
