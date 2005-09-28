#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzDisplacement.c
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
#include <ctype.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>
#include <bibFile.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

static void     usage(char *proc_str); 
static WlzBasisFnTransform *WlzAffineBasisFnTransformFromCPts(WlzObject *obj,
				WlzFnType basisFnType, int polyOrder,
				int nSPts, WlzDVertex2 *sPts,
				int nDPts, WlzDVertex2 *dPts,
							WlzErrorNum *dstErr);
static WlzCompoundArray *WlzBasisFnTransformObjPrv(WlzObject *inObj, 
					 WlzBasisFnTransform *basisTr, 
					 WlzErrorNum *dstErr);

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
  char 		*rec,
  		*srcFileStr,
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
  WlzBasisFnTransform *basisTr = NULL;
  WlzTransformType     trType;
  WlzFnType basisFnType;
  WlzMeshGenMethod genMethod;
  WlzInterpolationType interp = WLZ_INTERPOLATION_NEAREST;
  WlzErrorNum	errNum = WLZ_ERR_NONE; 
  BibFileRecord	       *bibfileRecord;
  BibFileError         bibFileErr;
  
  /* read the argument list and check for an input file */

  static char	       optList[] = "s:b:d:h";
 
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
        case 'd':
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
      basisTr = WlzAffineBasisFnTransformFromCPts(inObj, basisFnType, 
					    basisFnPolyOrder,
					    nTiePP, srcVtx2, nTiePP, dstVtx2, 
					    &errNum);
  }
  
  /* Calculate displacements from the source object and then write it out. */
  if(errNum == WLZ_ERR_NONE)
  {
    outObj = WlzBasisFnTransformObjPrv(inObj, basisTr, &errNum);
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
  (void )WlzBasisFnFreeTransform(basisTr);
}

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s [-h] \n"
	  "\tDetermine displacement object.\n"
	  "\n"
	  "\n"
	  "\tOptions are:\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -s        Source file\n"
          "\t  -b        Bib files\n"
	  "\t  -d        Output displacement file\n"
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
static WlzBasisFnTransform *WlzAffineBasisFnTransformFromCPts(WlzObject *obj,
				WlzFnType basisFnType, int polyOrder,
				int nSPts, WlzDVertex2 *sPts,
				int nDPts, WlzDVertex2 *dPts,
				WlzErrorNum *dstErr)
{
  int		idx;
  WlzDVertex2	tDV0;
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
    printf("nDpts = %d\n", nDPts);
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
			  	   nSPts, sPts, nSPts, dPtsT,NULL, &errNum);
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
* \brief	Computes a compound array of two objects. One with the
*               x component of displacement, the other with the y component.
* \param	inobj			Given object.
* \param        basisTr                 Basis Function transform
* \param	dstErr			Destination error pointer, may be
*					NULL.
*/
static WlzCompoundArray *WlzBasisFnTransformObjPrv(WlzObject *inObj, 
					 WlzBasisFnTransform *basisTr, 
					 WlzErrorNum *dstErr)
{
  WlzObject 		*objX = NULL,
                        *objY = NULL;
  WlzCompoundArray	*cArray = NULL;
  WlzValues		valuesX,
                        valuesY;
  WlzIntervalWSpace     iWsp,
                        iWspX,
                        iWspY;
  WlzGreyWSpace	        gWspX,
                        gWspY;
  WlzDVertex2           sVtx,
                        dVtx;
  WlzErrorNum		errNum=WLZ_ERR_NONE;
  WlzObjectType		type;
  WlzObjectType         vType;
  WlzGreyP	        gPixX,
                        gPixY;
  WlzPixelV 		backgrnd;
  WlzFnType             basisFnType;  
  int 		        k,
                        l;

  basisFnType = basisTr->basisFn->type;
  valuesX.core = NULL;
  valuesY.core = NULL;  
  /* Create value tables for the two objects */
  vType = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_FLOAT, NULL);
  backgrnd.type = WLZ_GREY_FLOAT;
  backgrnd.v.inv = 0;
  if(valuesX.v = WlzNewValueTb(inObj, vType, backgrnd, &errNum))
  {
    valuesY.v = WlzNewValueTb(inObj, vType, backgrnd, &errNum);
  }
  else
  {
    WlzFreeValueTb(valuesX.v);
  }
  /* create two Wlz objects to hold x and y displacements */
  if (errNum == WLZ_ERR_NONE)
  {
    if (objX = WlzMakeMain(inObj->type, inObj->domain, valuesX, NULL, NULL, 
			 &errNum)) 
    {
      if (!(objY = WlzMakeMain(inObj->type, inObj->domain, valuesY, NULL, 
			       NULL, &errNum))) 
      {
	WlzFreeObj(objX);
	objX = NULL;
      }
    }
  }
  /* initialise workspaces */
  if (errNum == WLZ_ERR_NONE)
  {
    if ((errNum = WlzInitRasterScan(inObj, &iWsp, WLZ_RASTERDIR_ILIC)) == 
	WLZ_ERR_NONE)
    {
      if ((errNum = WlzInitGreyScan(objX, &iWspX, &gWspX)) == WLZ_ERR_NONE)
      {
	errNum = WlzInitGreyScan(objY, &iWspY, &gWspY);
      }
    }
  }
  /* Calculate displacements */
  switch(basisFnType)
  {
    case WLZ_FN_BASIS_2DMQ:
      if (errNum == WLZ_ERR_NONE)
      {
	while(((errNum = WlzNextInterval(&iWsp)) == WLZ_ERR_NONE) &&
	     ((errNum = WlzNextGreyInterval(&iWspX)) == WLZ_ERR_NONE) &&
	     ((errNum = WlzNextGreyInterval(&iWspY)) == WLZ_ERR_NONE))
	{
	  gPixX = gWspX.u_grintptr;
	  gPixY = gWspY.u_grintptr;
	  l = iWsp.linpos;
	  for (k = iWsp.lftpos; k <= iWsp.rgtpos; k++)
	  {
	    sVtx.vtX = k;
	    sVtx.vtY = l;
	    dVtx = WlzBasisFnValueMQ2D(basisTr->basisFn, sVtx);
	    *(gPixX.flp) = (float)dVtx.vtX;
	    *(gPixY.flp) = (float)dVtx.vtY;
	    ++(gPixX.flp);
	    ++(gPixY.flp);
	  }
	}
	if(errNum == WLZ_ERR_EOO)        
	{
	  errNum = WLZ_ERR_NONE;
	}
      }
      break;
    case WLZ_FN_BASIS_2DTPS:
      if (errNum == WLZ_ERR_NONE)
      {
	while(((errNum = WlzNextInterval(&iWsp)) == WLZ_ERR_NONE) &&
	     ((errNum = WlzNextGreyInterval(&iWspX)) == WLZ_ERR_NONE) &&
	     ((errNum = WlzNextGreyInterval(&iWspY)) == WLZ_ERR_NONE))
	{
	  gPixX = gWspX.u_grintptr;
	  gPixY = gWspY.u_grintptr;
	  l = iWsp.linpos;
	  for (k = iWsp.lftpos; k <= iWsp.rgtpos; k++)
	  {
	    sVtx.vtX = k;
	    sVtx.vtY = l;
	    dVtx = WlzBasisFnValueTPS2D(basisTr->basisFn, sVtx);
	    *(gPixX.flp) = (float)dVtx.vtX;
	    *(gPixY.flp) = (float)dVtx.vtY;
	    ++(gPixX.flp);
	    ++(gPixY.flp);
	  }
	}
	if(errNum == WLZ_ERR_EOO)        
	{
	  errNum = WLZ_ERR_NONE;
	}
      }
      break;
    default:
      errNum = WLZ_ERR_TRANSFORM_TYPE;
      break;
  }
  /* create compound object */
  if (errNum == WLZ_ERR_NONE)
  {
    cArray = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 1, 2, NULL, objX->type, 
				  &errNum);
    if (errNum == WLZ_ERR_NONE) 
    {
      cArray->o[0] = WlzAssignObject(objX, NULL);
      cArray->o[1] = WlzAssignObject(objY, NULL);
    }
  }

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cArray);
}
