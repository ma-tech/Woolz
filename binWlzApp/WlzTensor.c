#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzTensor_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzTensor.c
* \author       Angus Murray
* \date         September 2005
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
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
* \brief
* \ingroup	WlzTransform
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <float.h>
#include <Wlz.h>
#include <WlzExtFF.h>
#include <bibFile.h>

/* Externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;
#endif /* __STDC__ ] */

static void     		usage(
				  char *proc_str); 

static WlzBasisFnTransform 	*WlzAffineBasisFnTransformFromCPtsT(
				  WlzObject *obj,
				  WlzFnType basisFnType,
				  int polyOrder,
				  int nSPts, WlzDVertex2 *sPts,
				  int nDPts, WlzDVertex2 *dPts,
				  WlzErrorNum *dstErr);
static WlzCompoundArray 	*WlzBasisFnTensorTransformObjPrv(
				  WlzObject *inObj, 
				  WlzBasisFnTransform *basisTr, 
				  WlzErrorNum *dstErr);
static double 			WlzBasisFnValueMQ2DPrv(
				  WlzBasisFn *basisFn, 
				  WlzDVertex2 srcVx,
				  int component);
static double 			WlzBasisFnValueTPS2DPrv(
				  WlzBasisFn *basisFn, 
				  WlzDVertex2 srcVx,
				  int component);


int             main(int argc, char **argv)
{
  int		nTiePP,
		option,
		vxCount = 0,
		vxLimit = 0,
		basisFnPolyOrder = 3,
                index,
                relFlag,
                meshMinDist,
                meshMaxDist;
  char 		*srcFileStr,
		*bibFileStr = NULL,
                *outFileStr,
                *bibErrMsg;
  const char    *errMsg;
  FILE		*inFile = NULL,
                *outFile = NULL;
  WlzDVertex2   *sVtx2 = NULL,
                *dVtx2 = NULL,
                *srcVtx2 = NULL,
                *dstVtx2 = NULL;
  WlzDVertex3   *dVtx = NULL, 
                *sVtx = NULL,
                *srcVtx = NULL,   
                *dstVtx = NULL;
  WlzObject	*inObj = NULL;
  WlzCompoundArray  *outObj = NULL;              
  WlzBasisFnTransform  *basisTr;
  WlzTransformType     trType;
  WlzFnType basisFnType;
  WlzMeshGenMethod genMethod;
  WlzErrorNum	errNum = WLZ_ERR_NONE; 
  BibFileRecord	       *bibfileRecord;
  BibFileError         bibFileErr;
  
  /* read the argument list and check for an input file */

  static char	       optList[] = "s:b:t:h";
 
   while( (option = getopt(argc, argv, optList)) != EOF )
   {
      switch( option )
      {
        case 'h':
             usage(argv[0]);
	     return(0);
        case 'b':
	    bibFileStr = optarg;
	    break;
        case 's':
	    srcFileStr = optarg;
	    break;
        case 't':
	    outFileStr = optarg;
	    break;
        default:
              return(0);
      }
   }
   if((inFile = fopen(bibFileStr, "r")) == NULL )
   {
       printf("cannot open the input bib file.\n");
       exit(1);
   }
   /* read the bibfile until we get the WlzWarptransformParams part */

   while( !feof(inFile) ){
      bibFileErr = BibFileRecordRead(&bibfileRecord, &bibErrMsg, inFile);
      if(bibFileErr != BIBFILE_ER_EOF)
      {
	bibFileErr = WlzEffBibParseWarpTransformParamsRecord(
						bibfileRecord,
						&basisFnType,
						&trType,
						&genMethod,
						&meshMinDist,
						&meshMaxDist);
	 if (bibFileErr == BIBFILE_ER_NONE)
	   break;
      }
      else
	break;
   }

   /* read the bibfile until we get the TiePointVtxs part */
   if(  bibFileErr == BIBFILE_ER_NONE )
   {
     while( !feof(inFile) ){
       bibFileErr = BibFileRecordRead(&bibfileRecord, &bibErrMsg, inFile);
       if(bibFileErr != BIBFILE_ER_EOF)
       {
	    if(vxCount >= vxLimit)
	    {
	      vxLimit = (vxLimit + 1024) * 2;
	      if(((srcVtx = (WlzDVertex3 *)AlcRealloc(srcVtx,
			     vxLimit * sizeof(WlzDVertex3))) == NULL) ||
		 ((dstVtx = (WlzDVertex3 *)AlcRealloc(dstVtx,
			     vxLimit * sizeof(WlzDVertex3))) == NULL))
	      {
		errNum = WLZ_ERR_MEM_ALLOC;
	      }
	      else
	      {
		sVtx = srcVtx + vxCount;
		dVtx = dstVtx + vxCount;
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      bibFileErr = WlzEffBibParseTiePointVtxsRecord(
						bibfileRecord,
						&index,
						dVtx,
						sVtx,
						&relFlag);
	      if (bibFileErr == BIBFILE_ER_NONE)
	      {	      
		++sVtx;
		++dVtx;
		++vxCount;
	      }
	    }
      }
      else
	break;
     }
   }

   nTiePP = vxCount;

   BibFileRecordFree(&bibfileRecord);

   fclose(inFile);

   if(  bibFileErr == BIBFILE_ER_EOF )
   {
     errNum = WLZ_ERR_NONE;
   }
   else
   {
     printf("Read bib file error.\n");
     exit(1);
   }

   /* copy 3D verices to 2D vertices - cast won't work */
   if(((srcVtx2 = (WlzDVertex2 *)AlcRealloc(srcVtx2,
		  vxCount * sizeof(WlzDVertex2))) == NULL) ||
		  ((dstVtx2 = (WlzDVertex2 *)AlcRealloc(dstVtx2,
		  vxCount * sizeof(WlzDVertex2))) == NULL))
   {
     errNum = WLZ_ERR_MEM_ALLOC;
   }
   else
   {
     sVtx = srcVtx;
     dVtx = dstVtx;
     sVtx2 = srcVtx2;
     dVtx2 = dstVtx2;
     for (vxCount = 0; vxCount < nTiePP;  vxCount++)
     {
       sVtx2->vtX = sVtx->vtX;
       sVtx2->vtY = sVtx->vtY;
       dVtx2->vtX = dVtx->vtX;
       dVtx2->vtY = dVtx->vtY;
       ++sVtx2;
       ++dVtx2;
       ++sVtx;
       ++dVtx;
     }
   }

   errNum = WLZ_ERR_READ_EOF;
   if((srcFileStr == NULL) ||
	(*srcFileStr == '\0') ||
	((inFile = (strcmp(srcFileStr, "-")?
		fopen(srcFileStr, "r"): stdin)) == NULL) ||
       	((inObj= WlzAssignObject(WlzReadObj(inFile, &errNum), NULL)) == NULL) 
               || (errNum != WLZ_ERR_NONE))
   {
      (void )fprintf(stderr,
		     "%s: failed to read object from file %s\n",
		     *argv, srcFileStr);
   }
   if(inFile)
   {
      if(strcmp(srcFileStr, "-"))
      {
	fclose(inFile);
      }
      inFile = NULL;
   }
  if(errNum == WLZ_ERR_NONE)
  {
      basisTr = WlzAffineBasisFnTransformFromCPtsT(inObj, basisFnType, 
					    basisFnPolyOrder,
					    nTiePP, srcVtx2, nTiePP, dstVtx2, 
					    &errNum);
  }
  
  /* Calculate tensor components and then write them out. */
  if(errNum == WLZ_ERR_NONE)
  {
    outObj = WlzBasisFnTensorTransformObjPrv(inObj, basisTr, &errNum);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzStringFromErrorNum(errNum, &errMsg);
    (void )fprintf(stderr,
    	       "%s: failed to transform object (%s).\n",
       	       *argv, errMsg);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((outFile = (strcmp(outFileStr, "-")?
		  fopen(outFileStr, "w"):
		  stdout)) == NULL) ||
	          ((errNum = WlzWriteObj(outFile, (WlzObject *)outObj)) 
		   != WLZ_ERR_NONE))
    {
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
			 "%s: failed to write output compound array (%s).\n",
			 *argv, errMsg);
    }
    if(outFile && strcmp(outFileStr, "-"))
    {
      fclose(outFile);
    }
  }

  /* free memory */
  AlcFree(srcVtx);
  AlcFree(dstVtx);
  AlcFree(srcVtx2);
  AlcFree(dstVtx2);
  WlzBasisFnFreeTransform(basisTr);
  return(errNum);
}

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s [-h] \n"
	  "\tDetermine tensor object.\n"
	  "\n"
	  "\n"
	  "\tOptions are:\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -s        Source file\n"
          "\t  -b        Bib files\n"
	  "\t  -t        Output tensor file\n"
	  "\t                                                      \n"
					  "",
	  proc_str);
  return;
}

/*!
* \return	New BasisFn transform.
* \ingroup	WlzTransform
* \brief	Computes a basis function transform for the given object and a
*		set of control points using both an affine and a basis
*		function transform.
* \param	obj			Given object.
* \param	basisFnType		Required basis function type.
* \param	polyOrder		Order of polynomial, only used for
*					WLZ_FN_BASIS_2DPOLY.
* \param	nSPts			Number of source control points.
* \param	sPts			Source control points.
* \param	nDPts			Number of destination control points.
* \param	dPts			Destination control points.
* \param	dstErr			Destination error pointer, may be
*					NULL.
*/
static WlzBasisFnTransform *WlzAffineBasisFnTransformFromCPtsT(WlzObject *obj,
				WlzFnType basisFnType, int polyOrder,
				int nSPts, WlzDVertex2 *sPts,
				int nDPts, WlzDVertex2 *dPts,
				WlzErrorNum *dstErr)
{
  int		idx;
  WlzDVertex2	*dPtsT = NULL;
  WlzAffineTransform *aTr = NULL,
  		*aTrI = NULL;
  WlzBasisFnTransform *bTr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((nDPts <= 0) || (nDPts != nSPts))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((dPts == NULL) || (sPts == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  /* Compute least squares affine transform from the tie points. */
  if(errNum == WLZ_ERR_NONE)
  {
    aTr = WlzAffineTransformLSq2D(nSPts, sPts, nSPts, dPts, 0, NULL,
				  WLZ_TRANSFORM_2D_AFFINE, &errNum);
  }

  if(errNum == WLZ_ERR_NONE)
  {
    if(nSPts >= 4)
    {
      /* Create a new array of destination vertices which have the original
       * destination transformed by the inverse affine transform. */
      if((dPtsT = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) *
	      				   nSPts)) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      if(errNum == WLZ_ERR_NONE)
      {
        aTrI = WlzAffineTransformInverse(aTr, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        for(idx = 0; idx < nDPts; ++idx)
	{
          dPtsT[idx] = WlzAffineTransformVertexD2(aTrI, dPts[idx], NULL);
	}
        bTr = WlzBasisFnTrFromCPts2D(basisFnType, polyOrder,
			  	   nSPts, sPts, nSPts, dPtsT, NULL, &errNum);
      }
    }
  }
  AlcFree(dPtsT);
  (void )WlzFreeAffineTransform(aTr);
  (void )WlzFreeAffineTransform(aTrI);
  if(errNum != WLZ_ERR_NONE)
  {
      (void )WlzBasisFnFreeTransform(bTr);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(bTr);
}

/*!
* \return	Compound array.
* \ingroup	WlzTransform
* \brief	Computes a compound array of three objects, the objects 
                representing the t11, t12 qnd t22 componemts of a tensor.
* \param	inobj			Given object.
* \param        basisTensorTr           Basis Function tensor transform
* \param	dstErr			Destination error pointer, may be
*					NULL.
*/
static WlzCompoundArray *WlzBasisFnTensorTransformObjPrv(WlzObject *inObj, 
					 WlzBasisFnTransform *basisTr, 
					 WlzErrorNum *dstErr)
{
  WlzObject 		*objT11 = NULL,
                        *objT12 = NULL,
                        *objT22 = NULL;
  WlzCompoundArray	*cArray = NULL;
  WlzValues		valuesT11,
                        valuesT12,
                        valuesT22;
  WlzIntervalWSpace     iWsp,
                        iWspT11,
                        iWspT12,
                        iWspT22;
  WlzGreyWSpace	        gWspT11,
                        gWspT12,
                        gWspT22;
  WlzDVertex2           sVtx;
  WlzErrorNum		errNum=WLZ_ERR_NONE;
  WlzObjectType           vType;
  WlzGreyP	        gPixT11,
                        gPixT12,
                        gPixT22;
  WlzPixelV 		backgrnd;
  WlzFnType             basisFnType;
  double                t11Partial,
                        t12APartial,
                        t12BPartial,
                        t22Partial;
  int 		        k,
                        l;
 
  valuesT11.core = NULL;
  valuesT12.core = NULL;  
  valuesT22.core = NULL;
  basisFnType = basisTr->basisFn->type;
  /* Create value tables for the two objects */
  vType = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_FLOAT, NULL);
  backgrnd.type = WLZ_GREY_FLOAT;
  backgrnd.v.inv = 0;
  if((valuesT11.v = WlzNewValueTb(inObj, vType, backgrnd, &errNum)) != NULL)
  {
    if((valuesT12.v = WlzNewValueTb(inObj, vType, backgrnd, &errNum)) != NULL)
    {
      valuesT22.v = WlzNewValueTb(inObj, vType, backgrnd, &errNum);
    }
    else
    {
      WlzFreeValueTb(valuesT11.v);
      WlzFreeValueTb(valuesT12.v);
    }
  }
  else
  {
    WlzFreeValueTb(valuesT11.v);
  }
  /* create three Wlz objects to hold t11, t22 and t12 tensor components */
  if (errNum == WLZ_ERR_NONE)
  {
    if ((objT11 = WlzMakeMain(inObj->type, inObj->domain, valuesT11, NULL, 
			      NULL, &errNum)) != NULL) 
    {
      if ((objT12 = WlzMakeMain(inObj->type, inObj->domain, valuesT12, NULL, 
			        NULL, &errNum)) != NULL) 
      {
	if ((objT22 = WlzMakeMain(inObj->type, inObj->domain, valuesT22, 
				   NULL, NULL, &errNum)) == NULL)
	{
	  WlzFreeObj(objT11);
	  objT11 = NULL;
	  WlzFreeObj(objT12);
	  objT12 = NULL;
	}
      }
      else
      {
	WlzFreeObj(objT11);
	objT11 = NULL;
      }
    }
  }
  /* initialise workspaces */
  if (errNum == WLZ_ERR_NONE)
  {
    if ((errNum = WlzInitRasterScan(inObj, &iWsp, WLZ_RASTERDIR_ILIC)) == 
	WLZ_ERR_NONE)
    {
      if ((errNum = WlzInitGreyScan(objT11, &iWspT11, &gWspT11)) == 
	  WLZ_ERR_NONE)
      {
	if ((errNum = WlzInitGreyScan(objT12, &iWspT12, &gWspT12)) == 
	  WLZ_ERR_NONE)
	{
	  errNum = WlzInitGreyScan(objT22, &iWspT22, &gWspT22);
	}
      }
    }
  }
  /* Calculate tensor components for MQ2 only */
  if (errNum == WLZ_ERR_NONE)
  {
    switch (basisFnType)
    {
      case WLZ_FN_BASIS_2DMQ:
	while(((errNum = WlzNextInterval(&iWsp)) == WLZ_ERR_NONE) &&
	      ((errNum = WlzNextGreyInterval(&iWspT11)) == WLZ_ERR_NONE) &&
	      ((errNum = WlzNextGreyInterval(&iWspT12)) == WLZ_ERR_NONE) &&
	      ((errNum = WlzNextGreyInterval(&iWspT22)) == WLZ_ERR_NONE))
	{
	  gPixT11 = gWspT11.u_grintptr;
	  gPixT12 = gWspT12.u_grintptr;
	  gPixT22 = gWspT22.u_grintptr;
	  l = iWsp.linpos;
	  for (k = iWsp.lftpos; k <= iWsp.rgtpos; k++)
	  {
	    sVtx.vtX = k;
	    sVtx.vtY = l;
	    t11Partial = WlzBasisFnValueMQ2DPrv(basisTr->basisFn, sVtx, 0);
	    *(gPixT11.flp) = (float)t11Partial;
	    ++(gPixT11.flp);
	    t22Partial = WlzBasisFnValueMQ2DPrv(basisTr->basisFn, sVtx, 1);
	    *(gPixT22.flp) = (float)t22Partial;
	    ++(gPixT22.flp);
	    t12APartial = WlzBasisFnValueMQ2DPrv(basisTr->basisFn, sVtx, 2);
	    t12BPartial = WlzBasisFnValueMQ2DPrv(basisTr->basisFn, sVtx, 3);
	    *(gPixT12.flp) = 0.5 * ((float)t12APartial + (float)t12BPartial);
	    ++(gPixT12.flp);
	  }
	}
        if(errNum == WLZ_ERR_EOO)        
	{
	  errNum = WLZ_ERR_NONE;
	}
	break;
      case WLZ_FN_BASIS_2DTPS:
	while(((errNum = WlzNextInterval(&iWsp)) == WLZ_ERR_NONE) &&
	      ((errNum = WlzNextGreyInterval(&iWspT11)) == WLZ_ERR_NONE) &&
	      ((errNum = WlzNextGreyInterval(&iWspT12)) == WLZ_ERR_NONE) &&
	      ((errNum = WlzNextGreyInterval(&iWspT22)) == WLZ_ERR_NONE))
	{
	  gPixT11 = gWspT11.u_grintptr;
	  gPixT12 = gWspT12.u_grintptr;
	  gPixT22 = gWspT22.u_grintptr;
	  l = iWsp.linpos;
	  for (k = iWsp.lftpos; k <= iWsp.rgtpos; k++)
	  {
	    sVtx.vtX = k;
	    sVtx.vtY = l;
	    t11Partial = WlzBasisFnValueTPS2DPrv(basisTr->basisFn, sVtx, 0);
	    *(gPixT11.flp) = (float)t11Partial;
	    ++(gPixT11.flp);
	    t22Partial = WlzBasisFnValueTPS2DPrv(basisTr->basisFn, sVtx, 1);
	    *(gPixT22.flp) = (float)t22Partial;
	    ++(gPixT22.flp);
	    t12APartial = WlzBasisFnValueTPS2DPrv(basisTr->basisFn, sVtx, 2);
	    t12BPartial = WlzBasisFnValueTPS2DPrv(basisTr->basisFn, sVtx, 3);
	    *(gPixT12.flp) = 0.5 * ((float)t12APartial + (float)t12BPartial);
	    ++(gPixT12.flp);
	  }
	}
        if(errNum == WLZ_ERR_EOO)        
	{
	  errNum = WLZ_ERR_NONE;
	}
	break;
      default:
        errNum = WLZ_ERR_TRANSFORM_TYPE;
      break;
    }
  }
  /* create compound object */
  if (errNum == WLZ_ERR_NONE)
  {
    cArray = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 1, 3, NULL, 
				  objT11->type, &errNum);
    if (errNum == WLZ_ERR_NONE) 
    {
      cArray->o[0] = WlzAssignObject(objT11, NULL);
      cArray->o[1] = WlzAssignObject(objT22, NULL);
      cArray->o[2] = WlzAssignObject(objT12, NULL);
    }
  }

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cArray);
}

/*!
* \return	Partial derivative.
* \ingroup	WlzTransform
* \brief	Calculates the partial derivative for the given vertex using
*		a 2D multiquadric basis function.
* \param	basisFn			Basis function.
* \param	srcVx			Source vertex.
*/
static double 	WlzBasisFnValueMQ2DPrv(WlzBasisFn *basisFn, 
				       WlzDVertex2 srcVx, int component)
{
  int           count;
  double        tD0,
		tD1,
                tDS,
		delta;
  double        partial;
  WlzDVertex2   *basisCo,
		*cPts;
  WlzDVertex2   *polyVx;

  partial = 0.0;
  count = basisFn->nVtx;
  cPts = basisFn->vertices.d2;
  basisCo = basisFn->basis.d2;
  delta = *((double *)(basisFn->param));
  while(count-- > 0)
  {
    tD0 = srcVx.vtX - cPts->vtX;
    tD1 = srcVx.vtY - cPts->vtY;
    tDS = (tD0 * tD0) + (tD1 * tD1);
    if(tDS > DBL_EPSILON)
    {
      switch(component)
      {
        case 0:
	  tDS = tD0/sqrt(tDS + delta);
	  partial += basisCo->vtX * tDS;
	  break;
        case 1:
	  tDS = tD1/sqrt(tDS + delta);
	  partial += basisCo->vtY * tDS;
	  break;
        case 2:
	  tDS = tD1/sqrt(tDS + delta);
	  partial += basisCo->vtX * tDS;
	  break;
        case 3:
	  tDS = tD0/sqrt(tDS + delta);
	  partial += basisCo->vtY * tDS;
	  break;
      }
    }
    ++cPts;
    ++basisCo;
  }
  polyVx = basisFn->poly.d2;
  polyVx++;
  switch(component)
  {
    case 0:
      partial += polyVx->vtX;
      break;
    case 1:
      polyVx++;
      partial += polyVx->vtY;
      break;
    case 2:
      polyVx++;
      partial += polyVx->vtX;
      break;
    case 3:
      partial += polyVx->vtY;
      break;
  }

  return(partial);
}

/*!
* \return	Partial derivative.
* \ingroup	WlzTransform
* \brief	Calculates the partial derivative for the given vertex using
*		a 2D thin plate spline basis function.
* \param	basisFn			Basis function.
* \param	srcVx			Source vertex.
*/
static double 	WlzBasisFnValueTPS2DPrv(WlzBasisFn *basisFn, 
				       WlzDVertex2 srcVx, int component)
{
  int           count;
  double        tD0,
		tD1,
                tDS;
  double        partial;
  WlzDVertex2   *basisCo,
		*cPts;
  WlzDVertex2   *polyVx;

  partial = 0.0;
  count = basisFn->nVtx;
  cPts = basisFn->vertices.d2;
  basisCo = basisFn->basis.d2;
  while(count-- > 0)
  {
    tD0 = srcVx.vtX - cPts->vtX;
    tD1 = srcVx.vtY - cPts->vtY;
    tDS = (tD0 * tD0) + (tD1 * tD1);
    if(tDS > DBL_EPSILON)
    {
      switch(component)
      {
        case 0:
	  tDS = tD0 * (1 + log(tDS));
	  partial += basisCo->vtX * tDS;
	  break;
        case 1:
	  tDS = tD1 * (1 + log(tDS));
	  partial += basisCo->vtY * tDS;
	  break;
        case 2:
	  tDS = tD1 * (1 + log(tDS));
	  partial += basisCo->vtX * tDS;
	  break;
        case 3:
	  tDS = tD0 * (1 + log(tDS));
	  partial += basisCo->vtY * tDS;
	  break;
      }
    }
    ++cPts;
    ++basisCo;
  }
  polyVx = basisFn->poly.d2;
  polyVx++;
  switch(component)
  {
    case 0:
      partial += polyVx->vtX;
      break;
    case 1:
      polyVx++;
      partial += polyVx->vtY;
      break;
    case 2:
      polyVx++;
      partial += polyVx->vtX;
      break;
    case 3:
      partial += polyVx->vtY;
      break;
  }

  return(partial);
}

