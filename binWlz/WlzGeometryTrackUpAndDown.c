/***********************************************************************
* \Project:      Woolz
* \Title:        WlzGeometryTrackUpAndDown.c
* \Date:         28 Jan 2003 
* \Author:       J. Rao
* \Copyright:	 2003 Medical Research Council, UK.
*		 All rights reserved.
* \Address: 	 MRC Human Genetics Unit,
*		 Western General Hospital,
*		 Edinburgh, EH4 2XU, UK.
* \Purpose:      Interpolate for bib files 
*		
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <float.h>
#include <string.h>
#include <time.h>
#include <Wlz.h>

//#include <MAPaint.h>
//#include <MAWarp.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */
/*
static  WlzErrorNum WlzTetrahedronProducer(int nxmax, int nymax, int nzmax,
   double xmin, double xmax, double ymin , double ymax, double zmin, double zmax);
*/

static void GetGreyValue( WlzGreyValueWSpace *gVWSp, int kz, int jy, int ix, 
                     int *intensity,    WlzErrorNum  *dstErr );
		     
static void FillGreyValue( WlzGreyValueWSpace *gVWSp, int kz, int jy, int ix, 
                     int intensity,    WlzErrorNum  *dstErr );


static void  outPutSectionsForTest(
                  WlzObject   *WObjS, 
		  WlzObject   *WObj2D, 
		  FILE        *outFile, 
		  char        *outFileStr,
		  int          iloop,
		  WlzErrorNum *errNum);

char fileNameData[300][100];

static void usage(char *proc_str);  


int main(int	argc,
	 char	**argv)
{

  FILE	       *inFile = NULL;   /* used to read Woolz object */ 
  FILE	       *outFile = NULL;   /* used to read Woolz object */
  int           i,j,k, m, l1, label[3], individual[4];
  int           downOrUp = 0;
  int		option;
  int           numChar, numChar1, TotalN;
  int           outputAutoMeshVTK        = 0,
                outputAutoMeshWLZ        = 0,
                outputTransformedMeshWLZ = 0,
		outputCutPlaneAndCorrepSurfaceVTK        = 0,
		WARP                     = 1;
  int           basisFnPolyOrder = 3;
  int           Inumber, numP_in_Z= 10;
  int           globalCycle;
  int           binObjFlag;
  int           centorOfMass = 1;
  int           numOf2DWlzFiles = 0, numOfVtkFiles = 0;
  int           initialn0, endNum0,initialn1, endNum1;
  int           BibFileIndex[50];
  int           sectionLength_N       = 40;
  int           subSubSectionLength_L;
  int           numberOfSampleP_k     =  2;
  int           iloop=1, test=0;
  double	zConst                   = 0.;
  double        mass = 0.0;
  char          under ='_';
  char          ctemp;
  char         *inFileStr, *inFileStr3, *outFileStr;
  char          inFileStrSec[100];
  char          TwoDImageFilesNameList[700][120];
  char          TwoDImagePureFilesNameList[700][50];
  char          PureBibFilesNameList[700][50];
  char          SampleBibFilesNameList[50][120];
  char          SampleBibPureFilesName[50][50];
  char          SampleBibFilesDir[100];
  char          TwoDImageFilesDir[100];
  char         *StandContourFileName, *ToBeTrackedContourFileName;
  char         *ContourFilesNameList;
  char         *surfacePointFileName, *surfaceInPointFileName, *surfaceOutPointFileName;
  char         *Wlz2DFileName;
  char         *cstr;
  char          StringCh[100];
  const char   *errMsg;
  WlzDVertex3          *SurfacePatchPoints;
  WlzErrorNum	       errNum = WLZ_ERR_NONE;
  WlzInterpolationType interp = WLZ_INTERPOLATION_NEAREST;  /* Use the nearest neighbour */
  WlzThreeDViewStruct *wlzViewStr, *wlzViewStr1, *wlzViewStrInter;
  AlcErrno             alcErr = ALC_ER_NONE;
  WlzObject           *WObjS, *WObjTrack, *WObj2D;
  WlzObjectType        wtp = WLZ_3D_DOMAINOBJ;

  /* read the argument list and check for an input file */
  static char	optList[] = "S:s:w:t:N:n:o:O:i:l:L:j:J:h",

  opterr = 0;
 
  while( (option = getopt(argc, argv, optList)) != EOF )
  {
      switch( option )
      {
        case 'h':
             usage(argv[0]);
	     return(0);
        case 'L':
	     ContourFilesNameList       = optarg;
	     break;
        case 'S':
	     StandContourFileName       = optarg;
	     break;
        case 'N':
	     ToBeTrackedContourFileName = optarg;
	     break;
        case 'w':
	     Wlz2DFileName = optarg;
	     break;
        case 't':
	 if(sscanf(optarg, "%d", &downOrUp) != 1)
	     {
	       printf("read error");
	       exit(1);
	     }
	     break;
       case 'n':
	 if(sscanf(optarg, "%d", &numP_in_Z) != 1)
	     {
	       printf("read error");
	       exit(1);
	     }
	     break;
       case 'l':
	 if(sscanf(optarg, "%d", &iloop) != 1)
	     {
	       printf("read error");
	       exit(1);
	     }
	     break;
       case 'j':
	 if(sscanf(optarg, "%d", &numberOfSampleP_k) != 1)
	     {
	       printf("read error");
	       exit(1);
	     }
	     break;
       case 'J':
	 if(sscanf(optarg, "%d", &sectionLength_N) != 1)
	     {
	       printf("read error");
	       exit(1);
	     }
	     break;
        case 'o':
	     outFileStr = optarg;
	     break;
        case 's':
	     surfacePointFileName = optarg;
	     break;
        case 'i':
	     surfaceInPointFileName = optarg;
	     break;
        case 'O':
	     surfaceOutPointFileName = optarg;
	     break;
        default:
             return(0);
      }
  }

  
  subSubSectionLength_L = sectionLength_N/2 ;
  
  /* Read the contour Wlz file name List */
  if((inFile = fopen(ContourFilesNameList, "r")) == NULL )
  {
        printf("cannot open the 2D image files list file.\n");
        exit(1);
  }
  i=0;
  while(!feof(inFile))
  {
        fscanf(inFile, "%s", TwoDImageFilesNameList[i]);
	printf("%s\n", TwoDImageFilesNameList[i]);
	i++;
  }
  fclose(inFile); 
  inFile = NULL;
  numOf2DWlzFiles = i-1;
  printf("n of files = %d\n", numOf2DWlzFiles);
  //exit(0);



  /* Read the standard contour Wlz file */
  //printf("Read standarded contour Wlz first\n");
    i = 0;
    if((inFile = fopen(StandContourFileName, "r")) == NULL )
    {
        printf("cannot open the standard contour Wlz file .\n");
        exit(1);
    }
    
    if( !(WObjS = WlzReadObj(inFile, &errNum) ) )
    {
       printf("input Woolz Object Error.\n");
       fclose(inFile); 
       exit(1);
    }
     fclose(inFile); 
     inFile = NULL;

  /* Read the contour Wlz file to be tracked */
  //printf("Read the contour Wlz file to be tracked\n");
    if((inFile = fopen(ToBeTrackedContourFileName, "r")) == NULL )
    {
        printf("cannot open the contour Wlz file to be tracked.\n");
        exit(1);
    }
    
    if( !(WObjTrack = WlzReadObj(inFile, &errNum) ) )
    {
       printf("input Woolz Object Error.\n");
       fclose(inFile); 
       exit(1);
    }
     fclose(inFile); 
     inFile = NULL;

 
  /* Read the wlz file corresponding to standard contour Wlz */
  //printf("Read the Wlz file \n");
    if((inFile = fopen(Wlz2DFileName, "r")) == NULL )
    {
        printf("cannot open the  Wlz file \n");
        exit(1);
    }
    
    if( !(WObj2D = WlzReadObj(inFile, &errNum) ) )
    {
       printf("input Woolz Object Error.\n");
       fclose(inFile); 
       exit(1);
    }
     fclose(inFile); 
     inFile = NULL;

    /* ---- output for test ---- */
    if(test  &&  ( errNum == WLZ_ERR_NONE ) )
    {
      outPutSectionsForTest(WObjS, WObj2D, outFile, outFileStr, iloop, &errNum);
      //outPutSectionsForTest(WObjTrack, WObj2D, outFile, outFileStr, iloop, &errNum);
     }

    /*---------- extract contour and track down or up  to get the surface patch ----------- */
    if( (!test)  && ( errNum == WLZ_ERR_NONE ) )
    {
      SurfacePatchPoints = WlzGeometryTrackUpAndDown_s( WObjS,     
                                                        WObjTrack, 
							numP_in_Z,
                                                        TwoDImageFilesNameList,
                                                        numOf2DWlzFiles,
                                                        downOrUp,  
  							sectionLength_N,
  							subSubSectionLength_L,
  							numberOfSampleP_k,
							surfacePointFileName,
							surfaceInPointFileName,
							surfaceOutPointFileName, 
						       &errNum );
    }

   AlcFree(WObjS);
   AlcFree(WObj2D);
   AlcFree(WObjTrack);

  return ( 0 );
}

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s [-h] [<input file>]\n"
	  "\tAutomatically interpolate bib files.\n"
	  "\n"
	  "\tThe user needs to prepare some sample 2D Woolz File names which contains \n"
	  "\ta number of files with which the middle one will be used as a standard contour\n"
	  "\tand a start plan to track up and down.\n"
	  "\n"
	  "\tFor example:\n"
	  "\t/net/laphroaig/export/data1/rjianguo/newcastle/LizpDir/\n"
	  "\tsampleBibFilesDir/s11r2s06_template_g.wlz\n"
	  "\n"
	  "\tALSO, The base part of sample bib file's name should be the same as\n"
	  "\tthat of 2D Woolz file.\n"
	  "\t \n"
	  "\tFor example:\n"
	  "\t ABC.bib,  ABC.wlz ,\n"
	  "\twhere the base part of their name is the same\n"
	  "\n"
	  "\tPlease report bugs to RAO.JIANGUO@hgu.mrc.ac.uk\n"
	  "\t \n"
	  "\tOptions are:\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -L        input name of a file which contains the list of \n" 
          "\t                  of the 2D Woolz contour image file names  \n"
	  "\t                     \n"
	  "\t  -s        input name of a data file which used to output surface points\n"
	  "\t                     \n"
	  "\t  -w        input name of a Woolz file which will be used to output pictures\n"
	  "\t  -t        input number to decide track down or up\n"
	  "\t                     t = 0 track down(default) \n"
	  "\t                     t = 1 track up \n"
	  "\t                                                        \n"
	  "\t  -o        input name of a file to output the Woolz object to show the shells and loops\n"
	  "\t  -O        input name of a data file to output the out surface points\n"
	  "\t                        and their distance to the  the surface\n"
	  "\t  -i        input name of a data file to output the in surface points\n"
	  "\t                        and their distance to the  the surface\n"
	  "\t                                                        \n"
	  "\t  -n        input number of pixels for z-direction of vtk contour\n"
	  "\t  -j        input number of pixels to be sampled in the length of segmented line see J option\n"
	  "\t                     default is 2                                   \n"
	  "\t  -J        input length to be used to segment the line in unit of number of pixels\n"
	  "\t                     default is 40                                   \n"
	  "\t  -l        input number of lth loop for output \n"
	  "\t                  \n"
	        	  "", 
	  proc_str);
  return;
}


/*!
* - Function:   WlzEffWriteMeshTransform3DWithoutDisplacementVTK
* - Returns:    none 
* - Purpose:    output the orginal mesh.	
* - Parameters:	
*     -# *fp:               pointer pointing to a specific file.
*     -#  wmt3D:            mesh transform.
* - Author:       J. Rao, R. Baldock and B. Hill
*/
static void  outPutSectionsForTest(
                  WlzObject   *WObjS, 
		  WlzObject   *WObj2D, 
		  FILE        *outFile, 
		  char        *outFileStr, 
		  int          iloop,
		  WlzErrorNum *errNum){

  WlzGreyValueWSpace *gVWSp;
  WlzErrorNum       wErrN = WLZ_ERR_NONE;

  WlzGMModel         *gM;
  WlzGMShell         *cS, *fS;
  WlzGMLoopT         *cLT, *fLT;
  WlzGMEdgeT         *cET, *fET;
  WlzDVertex2        *v1,  *v2;
  WlzIBox2            bBoxT;
  WlzObject          *tObj;
  int ix, jy, kt, intensity;

 if( ( (WObjS->type) != WLZ_CONTOUR ) )
  {
      	wErrN = WLZ_ERR_DOMAIN_TYPE;
  }


  if( ( WObj2D->type != WLZ_2D_DOMAINOBJ )  )
  {
      	wErrN = WLZ_ERR_DOMAIN_TYPE;
  }

  if(  ( WObj2D->domain.core == NULL ) )
  {
	wErrN = WLZ_ERR_DOMAIN_NULL;
  }
  
  if(wErrN == WLZ_ERR_NONE)
  {
       /* copy it */
       tObj = WlzCopyObject(WObj2D, &wErrN);
       if( wErrN  != WLZ_ERR_NONE)
       {
          printf("Something wrong with when copy obj ");
          exit(1);
       }
   

       /* get the bounding box */
       bBoxT = WlzBoundingBox2I(tObj, &wErrN); 
       if( wErrN  != WLZ_ERR_NONE)
       {
          printf("Something wrong with getting bounding box of obj ");
          exit(1);
       }
       
       /* revise the grey value of it */   
       gVWSp        = WlzGreyValueMakeWSp(tObj, &wErrN);
       if( wErrN  != WLZ_ERR_NONE)
       {
          printf("Something wrong with when make warped obj Wsp");
          exit(1);
       }
       
       for(ix = bBoxT.xMin; ix <bBoxT.xMax; ix++)
       {
          kt = 0;
          for(jy= bBoxT.yMin; jy< bBoxT.yMax; jy++)
          {
	      // GetGreyValue(gVWSp, 0, jy, ix, &intensity, &errNum);
	      intensity = 255;
	      FillGreyValue(gVWSp, 0, jy, ix, intensity, &wErrN);
	      
          } 
       }

   }

  // visualize the data:
  if( wErrN == WLZ_ERR_NONE)
  {
    iloop = 1;

    kt = 0;
  
    gM = WObjS->domain.ctr->model;
    /* For each shell of the model. */
    cS = fS = (WlzGMShell *) gM->child;
    do
    {
       printf("Shell  %d   Shell index %d\n", kt, cS->idx);
      /* For each loop topology element of the model. */
       cLT = fLT = (WlzGMLoopT *) cS->child;
       kt++;
       do
       {
           printf("Loop index %d\n", cLT->idx);
           /* For each edge topology element of the model. */
	   cET = fET = cLT->edgeT;
	   do
	   {
              if(cET == cET->edge->edgeT) /* Print edge end points */
	      {
	         v1 = (WlzDVertex2 *) cET->vertexT->diskT->vertex->geo.vg2D;
	         v2 = (WlzDVertex2 *) cET->opp->vertexT->diskT->vertex->geo.vg2D;
		 /*
		 printf("%lg    %lg\n",cET->vertexT->diskT->vertex->geo.vg2D->vtx.vtX, \
		                       cET->vertexT->diskT->vertex->geo.vg2D->vtx.vtY );
		 printf("%lg    %lg\n",cET->opp->vertexT->diskT->vertex->geo.vg2D->vtx.vtX, \
		                       cET->opp->vertexT->diskT->vertex->geo.vg2D->vtx.vtY );
				       */
		//if(kt%3 == iloop)
		{
		 ix = (int) cET->vertexT->diskT->vertex->geo.vg2D->vtx.vtX;
		 jy = (int) cET->vertexT->diskT->vertex->geo.vg2D->vtx.vtY;
		 printf("idx %d ix  %d  iy  %d\n", cET->vertexT->diskT->vertex->idx, ix, jy);
	 	 intensity = kt * 40;
	         FillGreyValue(gVWSp, 0, jy, ix, intensity, &wErrN);
		 ix = (int) cET->opp->vertexT->diskT->vertex->geo.vg2D->vtx.vtX;
		 jy = (int) cET->opp->vertexT->diskT->vertex->geo.vg2D->vtx.vtY;
                 printf("idx %d ix  %d  iy  %d\n", cET->opp->vertexT->diskT->vertex->idx, ix, jy);
                 FillGreyValue(gVWSp, 0, jy, ix, intensity, &wErrN);
		} 
              }
	      cET = cET->next;
	   } while (cET != fET);
           cLT = cLT->next;


       } while(cLT != fLT);
       cS = cS->next;
   }  while(cS != fS );
  }

 
    if((outFile = fopen(outFileStr, "w")) == NULL )
    {
      printf("cannot open the output woolz file.\n");
      exit(1);
    }

    if( ( wErrN = WlzWriteObj(outFile, tObj) ) != WLZ_ERR_NONE )
    {
       printf("output Woolz Object Error.\n");
       fclose(outFile); 
       exit(1);
    }
    
    fclose(outFile); 
    outFile = NULL;

    WlzGreyValueFreeWSp(gVWSp);


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

