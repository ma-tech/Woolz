#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzGreyValueMixing_s.c
* \author       Jianguo Rao
* \date         April 2003
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
* \brief	Functions to mix the grey values of two woolz object
* 		and produce a new object using
* 		\f[
* 		o = (1 - x) o_1 + x o_2
 		\f]
* \ingroup	WlzFeatures
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <float.h>
#include <Wlz.h>
#include <math.h>

#define	IN_RECORD_MAX   (1024)

typedef struct _WlzMeshScanDElm
{
  int		valid;				        /* Non-zero if valid */
  double	xTr[3];         /* Affine transform coefficients for columns */
  double	yTr[3];           /* Affine transform coefficients for lines */
} WlzMeshScanDElm;

static void GetLookUpTableDist(int ***LUT, int mElem, int nElem, int TypeDis);

static void GetGreyValue( WlzGreyValueWSpace *gVWSp, int kz, int jy, int ix, 
                     int *intensity,    WlzErrorNum  *dstErr );
		     
static void FillGreyValue( WlzGreyValueWSpace *gVWSp, int kz, int jy, int ix, 
                     int intensity,    WlzErrorNum  *dstErr );


/*!
* \ingroup      WlzFeatures
* \return	Woolz obj contain the mixing grey value.
* \brief        calculate the distance map of a woolz obj	
* \param	sObj			Given Woolz grey value source object.
* \param	tObj			Given Woolz grey value target object.
* \param	xmiddle          	mixing parameters 
* \param	dstErr			Destination error pointer.
*/
WlzObject	*WlzGreyValueMixing_s(  WlzObject       *sObj,
                                        WlzObject       *tObj,
				        double         xmiddle,
				    WlzErrorNum        *dstErr
				 )
{
    WlzIBox3        bBox;
    WlzObject      *tempObj;
    int             M0, Mf;
    int             ix, iy, iz, xi, yi, zi, ixp, iyp, im, jm;
    int             distance;
    int             i, j, k, intensityS, intensityT, intensityW;
    int             mElem, nElem, lElem, dist, dist0;
    int          dmin, dmin0;
    AlcErrno        alcErr =   ALC_ER_NONE;
    WlzGreyValueWSpace *gVWSpS = NULL, *gVWSpT = NULL, *gVWSpW = NULL;
    WlzErrorNum  errNum = WLZ_ERR_NONE;



    /* we need bounding box transformation here which can be 
       done by using the mesh trasform */
    bBox =  WlzBoundingBox3I(sObj,  dstErr );


    
    /* Now make a new object by copy the source one */
    tempObj = WlzCopyObject(sObj, &errNum);
   if( errNum  != WLZ_ERR_NONE)
   {
        printf("Something wrong with when copy obj ");
        exit(1);
   }

   if(errNum == WLZ_ERR_NONE)
   {
      gVWSpS         = WlzGreyValueMakeWSp(sObj, &errNum);
   }
   
   if(errNum == WLZ_ERR_NONE)
   {
      gVWSpT         = WlzGreyValueMakeWSp(tObj, &errNum);
   }

   if(errNum == WLZ_ERR_NONE)
   {
      gVWSpW         = WlzGreyValueMakeWSp(tempObj, &errNum);
   }

   if(errNum != WLZ_ERR_NONE)
   {
        printf("error : %d\n", errNum);
        exit(1);
   }

   if(errNum == WLZ_ERR_NONE)
   {
      /* 3D case */
      if(   sObj->type  == WLZ_3D_DOMAINOBJ )
      {

          dmin0 = 3;
      }
      else if(sObj->type == WLZ_2D_DOMAINOBJ )  /* 2D case */
      {

          dmin0 = 2;
      }
      else
      {
          errNum = WLZ_ERR_DOMAIN_TYPE;
      }
    }
    
   if(errNum == WLZ_ERR_NONE)
   {
      /* 3D case */
      if(   tObj->type  == WLZ_3D_DOMAINOBJ )
      {

          dmin = 3;
      }
      else if(tObj->type == WLZ_2D_DOMAINOBJ )  /* 2D case */
      {
          dmin = 2;

      }
      else
      {
          errNum = WLZ_ERR_DOMAIN_TYPE;
      }
    } 

    if( dmin != dmin0 )
    {
         errNum = WLZ_ERR_DOMAIN_TYPE;
    }
   
   if(errNum == WLZ_ERR_NONE)
   {
      if( dmin0 == 3 )
      {
         mElem = bBox.xMax - bBox.xMin + 1;
         nElem = bBox.yMax - bBox.yMin + 1;
         lElem = bBox.zMax - bBox.zMin + 1;

         /* for each layer */
         for(iz = 0; iz < lElem; iz++)
         {
            zi = bBox.zMin + iz;
            for(iy = 0; iy < nElem; iy++)
            {
	        yi = bBox.yMin + iy;
                for(ix = 0; ix < mElem; ix++)
                {
		   xi = bBox.xMin + ix;
                   /* get the grey value */
                   GetGreyValue(gVWSpS, zi, yi, xi, &intensityS, &errNum);		   
                   GetGreyValue(gVWSpT, zi, yi, xi, &intensityT, &errNum);
		   /* mixing the grey value */
                   intensityW = (int) ( intensityS * (1. - xmiddle)  + intensityT * xmiddle  );
		   if(intensityW > 255)
		      intensityW = 255;
		   if(intensityW < 0)
		       intensityW = 0;
		   /* fill the grey value */
                   FillGreyValue(gVWSpW, zi, yi, xi, intensityW, &errNum);
		}

	    }

         }
      }
      else if(dmin0 == 2)
      {
         mElem = bBox.xMax - bBox.xMin + 1;
         nElem = bBox.yMax - bBox.yMin + 1;
            zi = 0;
            for(iy = 0; iy < nElem; iy++)
            {
	        yi = bBox.yMin + iy;
                for(ix = 0; ix < mElem; ix++)
                {
		   xi = bBox.xMin + ix;
                   /* get the grey value */
                   GetGreyValue(gVWSpS, zi, yi, xi, &intensityS, &errNum);		   
                   GetGreyValue(gVWSpT, zi, yi, xi, &intensityT, &errNum);
		   /* mixing the grey value */
                   intensityW = (int) ( intensityS * (1. - xmiddle)  + intensityT * xmiddle  );
		   if(intensityW > 255)
		      intensityW = 255;
		   if(intensityW < 0)
		       intensityW = 0;
		   /* fill the grey value */
                   FillGreyValue(gVWSpW, zi, yi, xi, intensityW, &errNum);
		}

	    }
      

      }
      else
      {
         errNum = WLZ_ERR_DOMAIN_TYPE;
      }

   }

   
    WlzGreyValueFreeWSp(gVWSpS);
    WlzGreyValueFreeWSp(gVWSpT);
    WlzGreyValueFreeWSp(gVWSpW);

    return tempObj;


}


static void GetGreyValue( WlzGreyValueWSpace *gVWSp, int kz, int jy, int ix, 
                     int *intensity,    WlzErrorNum  *dstErr )
{
          WlzGreyValueGetDir( gVWSp, kz, jy, ix );
         switch(gVWSp->gType)
             {
	         case WLZ_GREY_INT:
	              *intensity =  ( int ) (*(gVWSp->gVal)).inv;
	              break;
	         case WLZ_GREY_SHORT:
	              *intensity =  ( int ) (*(gVWSp->gVal)).shv;
	              break;
	         case WLZ_GREY_UBYTE:
	              *intensity =  ( int ) (*(gVWSp->gVal)).ubv;
	              break;
	         case WLZ_GREY_FLOAT:
	              *intensity =  ( int ) (*(gVWSp->gVal)).flv;
	              break;
	         case WLZ_GREY_DOUBLE:
	              *intensity =  ( int ) (*(gVWSp->gVal)).dbv;
	              break;
	         default:
	              *dstErr = WLZ_ERR_GREY_TYPE;
	              break;
              }
}


static void FillGreyValue( WlzGreyValueWSpace *gVWSp, int kz, int jy, int ix, 
                     int intensity,    WlzErrorNum  *dstErr )
{
            WlzGreyValueGetDir( gVWSp, kz, jy, ix );
            switch(gVWSp->gType)
            {
	      case WLZ_GREY_INT:
	              *(gVWSp->gPtr[0].inp) = intensity; 
	              break;
	      case WLZ_GREY_SHORT:
	              *(gVWSp->gPtr[0].shp) = intensity;
	              break;
	      case WLZ_GREY_UBYTE:
	              *(gVWSp->gPtr[0].ubp) = intensity;
	              break;
	      case WLZ_GREY_FLOAT:
	              *(gVWSp->gPtr[0].flp) = intensity;
	              break;
	      case WLZ_GREY_DOUBLE:
	              *(gVWSp->gPtr[0].dbp) = intensity;
	              break;
	      default:
	              *dstErr = WLZ_ERR_GREY_TYPE;
	              break;
            }
}
