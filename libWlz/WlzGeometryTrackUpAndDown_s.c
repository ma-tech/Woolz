/* This code started on 28/01/2003 by J. Rao  
         it contains codes to find the edges 
*/
/************************************************************************
*
* Function:     WlzGeometryTrackUpAndDown_s.c	
* Returns:	WlzDVertex3
* Purpose: 	Giving the standard contour Woolz and the contour Woolz object 
*               to be tracked, The code can forme a patch surface points.
*
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <Wlz.h>

#define	NumberToTrack   (5)
#define MaxNumOfFiles   (100)
#define TOL 1.0e-5

static void testConsistencyBreakDownTheGM(  WlzGMModel *gModel[], 
                        int numOf2DWlzFiles,
			int sectionLength_N,
                        WlzErrorNum   *dstErr );

static void testBeforeBreakDownTheGM(  WlzGMModel *gModel[], 
                        int numOf2DWlzFiles,
			int sectionLength_N,
                        WlzErrorNum   *dstErr );

static void testJustAfterBreakDownTheGM(  WlzGMModel *gModel[], 
                        int numOf2DWlzFiles,
			int sectionLength_N,
                        WlzErrorNum   *dstErr );


static void svdfit(double x[], double y[], double sig[], int ndata, double a[],
            int ma, double **u, double **v, double w[], double *chisq, 
	    void (*funcs)(double, double [], int) );
	    
static void svdcmp(double **a, int m, int n, double *w, double **v);


static WlzDVertex2 *WlzDVerticesThisLoopSimpleShellGM2(WlzGMModel  *model, 
                                            int            LoopIdx, 
					    int           *numOfVertices, 
					    WlzErrorNum   *dstErr   );



static WlzVertexP WlzVerticesThisLoopOfSimpleShell2DGM(  WlzGMModel  *model, 
					    int          LoopIdx,
					    int         *numOfVertices,
					    WlzErrorNum *dstErr   );


static int *WlzGetAllHighConVertOfGM2D( WlzGMModel *sGM,
                                            int *numHV,  
					    WlzErrorNum *disErr);

static WlzErrorNum WlzRemoveVertices(
                                WlzGMModel *sGM, int *iBuf, int nD);


static void   outputSurfaceByPatch( int j_shell,
				    int i_loop,
				    int i_sec,
				    int nOfFiles,
  				    int subSubSectionLength_L,
				    int           *numOfPInS,
                                    WlzIVertex2   *sectionData[MaxNumOfFiles],
				    FILE          *surfacefp,
				    char          *surfaceStr,
				    FILE          *infp,
				    char          *inStr,
				    FILE          *outfp,
				    char          *outStr,
		                    WlzDVertex3   *SurfacePoints, 
				    int            nOfVOnTheS,
				    WlzDVertex3   *InPoints,      
				    int            nOfVWithinTheS,
				    WlzDVertex3   *OutPoints,     
				    int            nOfVOutTheS,
				    WlzErrorNum *dstErr );

static void svdfit(double x[], double y[], double sig[], int ndata, double a[],
                       int ma, double **u, double **v, double w[], double *chisq ,
		        void (*funcs)(double, double [], int) );

static void funcs(double x, double y[], int i);

double    **matrix(int nrl, int nrh, int ncl, int nch);

double     *vector(int nl, int nh);

void        free_vector(double *v, int nl, int nh);

void        free_matrix(double **m, int nrl, int nrh, int ncl, int nch);


static double polynomialF( double x, double ac[], int polyDeg);

void nrerror(char error_text[]);


static void  outputThisSection( WlzVertexP AllVerticesInThisStartLoop, 
                                int          i_s, 
                                int          subSubSectionLength_L,
                                WlzErrorNum *dstErr );

static void  outputThisSectionForView( WlzVertexP AllVerticesInThisStartLoop, 
                                int          i_s, 
                                int          subSubSectionLength_L,
                                WlzErrorNum *dstErr );


static double SixTimesOfVoulumeOfTetraHedron( WlzDVertex3 p0, 
			                      WlzDVertex3 p1,
					      WlzDVertex3 p2,
					      WlzDVertex3 p3 );

static int IsDVertexZero( WlzDVertex3  normal);

static int IsDVertexEqual( WlzDVertex3 p1, WlzDVertex3 p2);

static void RecalInOutPBy3DNorm(
                                  int            numOf2DWlzFiles,
				  double         distance,
		                  WlzDVertex3   *SurfacePoints, 
				  int            nOfVOnTheS,
				  WlzDVertex3   *InPoints,      
				  int           *nOfVWithinTheS,
				  WlzDVertex3   *OutPoints,     
				  int           *nOfVOutTheS );

static void pureThisGM(WlzGMModel *gMC, 
		       int LoopIdx, 
		       int numOfEnds, 
                       int *EndsIdx, 
		       int *bifurIdx, 
		       WlzErrorNum *dstErrgMC);

static void outputVerticesInThisLoop(WlzGMModel *gMC, int LoopIdx, WlzErrorNum *dstErr);

static void WlzGMModelPure( WlzGMModel *gMC, WlzErrorNum *dstErr);

static void GetInOutAndFromSurfacePoints( 
				   int          numberOfPixelsZ,
                                   WlzDVertex2  StandSampleP[NumberToTrack],
                                   double       distance,
				   WlzDVertex3 *InPoints,
				   int         *nOfVWithinTheS,
				   WlzDVertex3 *OutPoints,
				   int         *nOfVOutTheS,
				   WlzErrorNum *dstErr );

static void	outputSurfaceInAndOutPoints(int j_shell,
                                            int i_loop,
					    int i_sec,
					    int nOfFiles,
					    int           *numOfPInS,
                                            WlzIVertex2   *sectionData[MaxNumOfFiles],
					    FILE          *surfacefp,
					    char          *surfaceStr,
					    FILE          *infp,
					    char          *inStr,
					    FILE          *outfp,
					    char          *outStr,
		                            WlzDVertex3   *SurfacePoints, 
					    int            nOfVOnTheS,
					    WlzDVertex3   *InPoints,      
					    int            nOfVWithinTheS,
					    WlzDVertex3   *OutPoints,     
					    int            nOfVOutTheS,
					    WlzErrorNum *dstErr);
static void GetSamplePointsFromAllVertices(   WlzVertexP   sLoopVData, 
                                              int          i_s,
					      int          subSubSectionLength_L,       
                                              WlzDVertex2  StandSampleP[NumberToTrack],
                                              WlzErrorNum *dstErr );

static int NumberOfVerticesDoubleInNonCycleCaseInTheLoopOfGM( WlzGMModel  *gM, 
 						int sLoopIdx, 
						WlzErrorNum *dstErr);

static int *GetThebifurcatePointsIndexInTheNonCycleLine(  WlzGMModel *gM, 
                                           int          LoopIdx,
					   int          numOfEnds,
                                           WlzErrorNum *dstErr );


static int *GetTheEndsPointsIndexInTheNonCycleLine(  WlzGMModel *gM, 
                                           int          LoopIdx,
					   int          numOfEnds,
                                           WlzErrorNum *dstErr );

static int HowManyEndsInTheNonCycleLine(  WlzGMModel  *gM, 
                                           int          LoopIdx, 
                                           WlzErrorNum *dstErr );

static int IsACycle(  WlzGMModel  *gM, int  LoopIdx, WlzErrorNum *dstErr );

static void  WlzOutputForFormingSurface( FILE *testFile, 
					 char *Fstr, 
					 WlzVertexP sLoopVData, 
					 int nV, 
					 double zC, 
					 int sampleD, 
					 WlzErrorNum *dstErr);

static void  WlzOutputForTestGM(WlzGMModel *gM, WlzErrorNum *dstErr);

static void  WlzOutputForTestAsectionData(WlzVertexP   sLoopVData, int nV, WlzErrorNum *dstErr);

static int  nOfEdgeEndVInOneLoopOfGM(WlzGMModel *gM, int nThShell, int nThLoop);

static int  nOfVInOneLoopOfGM(WlzGMModel *gM, int nThShell, int nThLoop);

static int *WlzShellsIndexAboutGM(WlzGMModel *gM, int nS, WlzErrorNum *dstErr);

static int *WlzLoopIndexAboutGM(WlzGMModel *gM, int nThSell, int numLoops, WlzErrorNum *dstErr);

static void FactsAboutGM(WlzGMModel *gM);

static int NumberOfShellsAboutGM(WlzGMModel *gM );

static int NumberOfLoopsInThe_nThShellOfGM(WlzGMModel *gM, int nThShell );

static int NumberOfVerticesInTheLoopOfGM(WlzGMModel   *gM, 
 						int sLoopIdx, 
						WlzErrorNum *dstErr);  

static void WlzGetTheShellAndLoopNumberFromTheIndexOfTheVertex(WlzGMModel *gM,
                                           int idx, 
					   int *nShell, 
					   int *nLoop,
					   WlzErrorNum *dstErr);

static WlzVertexP WlzEdgeEndVerticesOneLoopFromGM(WlzGMModel *gM, 
                                           int   nThShell, 
					   int   nThLoop, 
					   WlzErrorNum *dstErr);

static WlzDVertex2 *WlzDEdgeEndVerticesOneLoopFromGM2(WlzGMModel *model, 
                                            int          nThShell, 
					    int          nThLoop, 
					    WlzErrorNum *dstErr   );

static WlzDVertex3 *WlzDEdgeEndVerticesOneLoopFromGM3(WlzGMModel *model, 
                                            int          nThShell, 
					    int          nThLoop, 
					    WlzErrorNum *dstErr   );

static WlzVertexP WlzVerticesOneLoopFromGM(WlzGMModel *gM, 
                                           int   nThShell, 
					   int   nThLoop, 
					   WlzErrorNum *dstErr);

					   
static WlzDVertex2 *WlzDVerticesOneLoopFromGM2(WlzGMModel *model, 
                                            int          nThShell, 
					    int          nThLoop, 
					    WlzErrorNum *dstErr   );
					    
static WlzDVertex3 *WlzDVerticesOneLoopFromGM3(WlzGMModel *model, 
                                            int          nThShell, 
					    int          nThLoop, 
					    WlzErrorNum *dstErr   );

static WlzVertexP WlzVerticesThisLoopOfGM(WlzGMModel  *model, 
                                            int            LoopIdx, 
					    int            numOfVertices, 
					    WlzErrorNum   *dstErr   );

static WlzDVertex2 *WlzDVerticesThisLoopOfGM2(WlzGMModel  *model, 
                                            int            LoopIdx, 
					    int            numOfVertices, 
					    WlzErrorNum   *dstErr   );

static WlzDVertex3 *WlzDVerticesThisLoopOfGM3(WlzGMModel  *model, 
                                            int            LoopIdx, 
					    int            numOfVertices, 
					    WlzErrorNum   *dstErr   );

static WlzDVertex2 *WlzVerticesThisIndexOfGM2(WlzGMModel  *model, 
                                            int            Idx, 
					    WlzErrorNum   *dstErr   );

static WlzDVertex2 WlzVerticesNormTriple2(WlzDVertex2 vA, WlzDVertex2 vB,
					   WlzDVertex2 vC);

static WlzDVertex2 WlzVerticesNormPair2(WlzDVertex2 v0, WlzDVertex2 v1);

static int *GroupDivid(int LoopIndex[], int *NumOfGroup, int *Continue, int nt, WlzErrorNum   *dstErr);

static int *GroupByNumber(int LoopIndex[], int *NumOfGroup, int nt,  WlzErrorNum *dstErr);

static int BigestGroup(int  GroupNum[], int nt,  WlzErrorNum *dstErr);

static void InAndOutSurfacePoints( WlzGMModel  *gM2, 
                                   int          nThLoop, 
				   int          numberOfPixelsZ,
                                   WlzDVertex2  StandSampleP[NumberToTrack],
                                   double       distance,
				   WlzDVertex3 *InPoints,
				   int         *nOfVWithinTheS,
				   WlzDVertex3 *OutPoints,
				   int         *nOfVOutTheS,
				   int         *suc,
				   WlzErrorNum *dstErr );

static void GetSamplePointsFromOneLoopOfGM( WlzGMModel  *gM, 
                                   int          sloopIdx,  
                                   WlzDVertex2  StandSampleP[NumberToTrack],
                                   WlzErrorNum *dstErr );

static void GetTrackedSamplePointsFromOneLoopOfGM( WlzGMModel  *gM2, 
                                    AlcKDTTree         *tTree,
                                   int          *LoopIdx, 
				   int           numberOfPixelsZ,
                                   double        minDis,
                                   WlzDVertex2   StandSampleP[NumberToTrack],
                                   WlzDVertex2   TrackedSamplePBlgToOneLp[NumberToTrack],
				   int          *nOfTracked,
				   int          *suc,
                                   WlzErrorNum  *dstErr );

static AlcKDTTree *GetTreeFromGM( WlzObject  *tObj,
                                  WlzErrorNum *dstErr );

static void GetTrackedPlusSomeSamplePointsFromOneLoopOfGM(
                           int          j_shell,
			   int          i_sec,
			   int          i_file,
			   int          sectionLength_L,
                           int         *numOfPInS,
                           WlzIVertex2 *sectionData,
			   WlzGMModel  *gM, 
                           int          LoopIdx, 
                           WlzDVertex2  StandSampleP[NumberToTrack],
		   	   WlzDVertex2	TrackedSamplePBlgToOneLp[NumberToTrack],
			   int		nOfTracked,
			   int         *suc,
                           WlzErrorNum *dstErr );

static void GetSamplePointsForNextTrackFromTheTrackedLoopOfGM(
                           int          i_file,
			   int          sectionLength_L,
                           int         *numOfPInS,
                           WlzIVertex2 *sectionData,
			   WlzGMModel  *gM, 
                           int          LoopIdx, 
                           WlzDVertex2  StandSampleP[NumberToTrack],
		   	   WlzDVertex2	TrackedSamplePBlgToOneLp[NumberToTrack],
			   int		nOfTracked,
			   int         *suc,
                           WlzErrorNum *dstErr );

/*!
* \ingroup	WlzFeatures:	WlzGeometryTrackUpAndDown_s.c
* \return:	WlzDVertex3:  The nearest edge point
*		
* \brief:    track down a curve lines 
* \ref:	-
* \param:	
*      \param   *sObj:	                    Woolz contour object (input)
*      \param   *tObj:	                    Woolz contour object (input)
*      \param    numberOfPixelsZ:	            the input z-direction of pixel number (input)
*      \param  **TwoDImageFilesNameList         2D image FIles name list. ( input )
*      \param    downOrUp:                      parameter to track down or up z-direction
*      \param    sectionLength_N                length used to cut a patch (unit in pixel) (input )
*      \param    subSubSectionLength_L          length smalle than the above used to sample points
*                                                  in the region
*      \param    numberOfSampleP_k              number of points will be sampled in the above section
*      \param   *surfacePointFileName:          FileNameStr to output the surface points
*      \param   *surfaceInPointFileName:        FileNameStr to output the in  surface points 
*                                                  and their distances to the surface
*      \param   *surfaceOutPointFileName:       FileNameStr to output the out surface points
*                                                  and their distances to the surface
*      \param    startShell                     the n-th Shell to begin with tracking
*      \param    endShell                       the end  Shell from where stop tracking 
*      \param    startSection                   the n-th Section to begin with tracking
*      \param    endSection                     the section to stop tracking.
*      \param   
*      \param#  *dstErr:	        Destination error pointer, may be NULL.
* \author:       J. Rao, R. Baldock and B. Hill
*/
WlzDVertex3  *WlzGeometryTrackUpAndDown_s(  
                                           int                 numberOfPixelsZ,
					   int                 startTrackingFile,
					   int                 numberOfFilesDownOrUp,
					   double	       disForInAndOutGuid, 
                                           double	       disForInAndOut, 
					   unsigned  char    **TwoDImageFilesNameList,
					   int                 numOf2DWlzFiles,
				           int                 downOrUp,
                                           int                 sectionLength_N,
                                           int                 subSubSectionLength_L,
                                           int                 numberOfSampleP_k,
					   char               *surfacePointFileName,
					   char               *surfaceInPointFileName,
					   char               *surfaceOutPointFileName,
					   int                 startShell,
					   int                 endShell,
					   int                 startSection,
					   int                 endSection,
					   double              minDis,
		                           WlzErrorNum        *dstErr
		                     )
{
  WlzErrorNum	      errNum = WLZ_ERR_NONE;
  WlzGMModel         *gM, *gMS, *gMT, *gMC;
  WlzGMModel         *gM2, *gM3,*gModel[numOf2DWlzFiles];

  int                *sNN, suc;
  AlcKDTTree         *tTree;
  AlcKDTTree         *tTreeArray[numOf2DWlzFiles];
  AlcKDTNode         *node;
  AlcErrno            errNo = ALC_ER_NONE;

  double             datD[3];
   
  WlzGMShell         *cS, *fS;
  WlzGMLoopT         *cLT, *fLT;
  WlzGMEdgeT         *cET, *fET;
  WlzDVertex2        *v1,  *v2, *vNorm=NULL;
  WlzDVertex2	      segV[3];
  WlzDVertex2         StandSampleP[NumberToTrack], StandSampleP0[NumberToTrack] ;
  WlzDVertex2         TrackedSampleP[NumberToTrack], TrackedSamplePBlgToOneLp[NumberToTrack];
  WlzDVertex2         TrackedSampleFirstP, TrackedSampleSecondP;
  WlzIVertex2         wV2, wV2I;
  WlzDVertex2         wD1, wD2;
  WlzDVertex3        *surfP, *SurfacePoints, *InPoints, *OutPoints;
  int                *numOfPInS, numInThisGM;  
  WlzIVertex2        *sectionData[MaxNumOfFiles]; 
  WlzDVertex2         fPoint, *inPoints, *outPoints;
  WlzVertexP          nData[2], vData[2];
  WlzVertexP          sLoopVData, tLoopVData, AllVerticesInThisStartLoop;
  WlzVertexP          outVerticesForGM;
  WlzVertexType       vType;
   
  WlzObject          *newSObj, *newTObj, *WObj[numOf2DWlzFiles];
  int                 indexStr[NumberToTrack], LoopIdx, ifileStar;
  int                 LoopIndex[NumberToTrack];
  int                 ShellIndex[NumberToTrack];
  int                 *GroupNum, *FixNumForEachGroup;
  int                 NumOfGroup;
  int                 sp, nbgN, trackedBigestLoopIndex;
  int                 loopDirection, ncycals;
  int                 numOfLoopsInStartShell, numOfShellsInStartGM;
  int                 numOfVerticesInTheStartLoop, numOfSections, i_s;
  int                 nOfVOnTheS, nOfVWithinTheS, nOfVOutTheS;
  int                *iBuf, numHV;
  int                 endLoop, ifile;
  int                 numberN1, numberN2;

  AlcVector          *vec;
   
  int                 idx, max, nShell, nLoop;
  int                 intensity, intensityB;
  int                 intensityT;
  int                 NumOfDiv, sampleD = 10;
  int                *sShellIdx, *tShellIdx, *sLoopIdx, *tLoopIdx;
  int                *startShellIdx, *startLoopIdx;
  int                *EndPointsIdx, *BifurcationPointsIdx, numOfEnds;
  WlzIBox2            bBoxS, bBoxT, bBoxO, bBoxWarp;
  int                 imax, nOfTracked;

  int                 i,j, k, ix, jy, it, jt, kt, test = 0, nE, i_l;
  int                 loopForFile, loopForShell, loopForLoop, loopForSection;
  int                 tempI, tempJ, tempK, tempL, tempM, tempN;
  int                 nS, nL, nV, ntS, ntL, ntV;
  int                 Continue;
  double              TempD, Tempd;
  double              di, dj, dk;
  double              Sumfr_halpha, Ialpha;
  double              HA, HTB, HATB, YATB, MI;
  double              wg[4];
  double              h,w, zC;
  double              zdvoxels;
  int                 vCnt[2], idN;
  int                 deltaX, deltaY, startI, fianlI;
  FILE                *testFile = NULL, *surfacefp, *infp= NULL, *outfp=NULL;
  char                *testFileStr = "outputForTest.dat";

   ifileStar = startTrackingFile;
   /* ifileStar = 10; */
 
   /* read the wlz object: */
   if( errNum == WLZ_ERR_NONE)
   {
       for(loopForFile=0; loopForFile<numOf2DWlzFiles; loopForFile++)
       {
	        if((testFile = fopen( (const char *) TwoDImageFilesNameList[loopForFile], "r")) == NULL )
                {
         		printf("cannot open the standard contour Wlz file .\n");
         		exit(1);
      		}

      		if( !(WObj[loopForFile] = WlzReadObj(testFile, &errNum) ) )
      		{
         		printf("input Woolz Object Error.\n");
         		fclose(testFile); 
         		exit(1);
      		}
      		fclose(testFile); 
      		testFile = NULL;
		
     		if(  ( (WObj[loopForFile]->type) != WLZ_CONTOUR ) )
      		{
       	 		errNum = WLZ_ERR_DOMAIN_TYPE;
      		}

      		if(  ( WObj[loopForFile]->domain.core == NULL ) )
      		{
	 		errNum = WLZ_ERR_DOMAIN_NULL;
      		}
       }
   }

   /*  get g-model */
   if( errNum == WLZ_ERR_NONE)
   {
        for(loopForFile=0; loopForFile<numOf2DWlzFiles; loopForFile++)
        {
      	     if( errNum == WLZ_ERR_NONE)
      	     {
                   gModel[loopForFile] = WObj[loopForFile]->domain.ctr->model;
             }
	}
   }

  /* out put for test */

   test = 1;


  if( test && (errNum == WLZ_ERR_NONE) )
  {
    testBeforeBreakDownTheGM(  gModel, numOf2DWlzFiles, sectionLength_N, &errNum );

    for(i=0; i<numOf2DWlzFiles; i++)
    {
          WlzFreeObj(WObj[i]);
    }
    exit(0);
  }


   /*  break the model into simple one: */

   if( errNum == WLZ_ERR_NONE)
   {
        for(loopForFile=0; loopForFile<numOf2DWlzFiles; loopForFile++)
        {
 		
                /* break down the shell to simple one */
      		if( errNum == WLZ_ERR_NONE)
      		{
                   iBuf  = WlzGetAllHighConVertOfGM2D( gModel[loopForFile], &numHV, &errNum);
                }
		
      		if( errNum == WLZ_ERR_NONE)
      		{
		     if(test && (loopForFile == ifileStar) )
		     {
		       /* printf("number of HV:  %d\n",numHV); */
		       /*
		       for(k=0; k<numHV; k++)
		       {
                         printf("%d\n", *(iBuf + k) );
		       }
		       */
		     }  

		    errNum =  WlzRemoveVertices( gModel[loopForFile], iBuf, numHV);
		    AlcFree(iBuf);
                }

		if(errNum != WLZ_ERR_NONE)
		{
		  printf("problems with break down  the GM into simple structure");
                  exit(1);

		}
	}
 
   }
   
  /* output the vertices */
  /*
    outVerticesForGM =	WlzVerticesFromGM( gModel[0],
				           NULL, NULL,
				  &numInThisGM, &vType,
				  &errNum);

    for(k=0; k< numInThisGM; k++)
    {
       printf("%lg  %lg  %lg\n", ( outVerticesForGM.d2 + k )->vtX,
                                 ( outVerticesForGM.d2 + k )->vtY, 0.0 );
    }

	        for(k=0; k<numOf2DWlzFiles; k++)
	        {
	            if(WObj[k])
	               WlzFreeObj(WObj[k]);
	        }	


          exit(0);


   */

   /* out for test */
   test = 1;

  if( test && (errNum == WLZ_ERR_NONE) )
  {
    testJustAfterBreakDownTheGM(  gModel, numOf2DWlzFiles, sectionLength_N, &errNum );

    for(i=0; i<numOf2DWlzFiles; i++)
    {
          WlzFreeObj(WObj[i]);
    }
    exit(0);
  }

   /*  as kd tree is very expensive to build so build once and use them afterwards */ 
   if( errNum == WLZ_ERR_NONE )
   {
        /* begin the building of kd tree: */
        for(loopForFile=0; loopForFile<numOf2DWlzFiles; loopForFile++)
        {
      		if( errNum == WLZ_ERR_NONE)
      		{
	           /* get a KD Tree for the target vertices of the Woolz contour */
       		   tTreeArray[loopForFile] = GetTreeFromGM(WObj[loopForFile], &errNum);
                }
	}

   }	 




          i=0;
          datD[0] = 639.0;
          datD[1] = -94.0;
	  datD[2] = 0.;
          /*  printf("Point %lg %lg\n",datD[0], datD[1]); */
          node    = AlcKDTGetNN(tTreeArray[i], datD, minDis, NULL, &errNo);
     	  if(  node && (errNum == ALC_ER_NONE) )
  	  {
             printf("tracked\n");
	     printf("%lg  %lg\n", *(node->key.kD), *(node->key.kD + 1) );
	     printf("index %d\n", node->idx );
  	  }
	  else
	  {
             printf(" cannot track this points!!! you have \n \
	              to give a biger mindistance \n \
	              for nearest point to track! or there \
		      are something wrong!\n");
	  }



	        for(k=0; k<numOf2DWlzFiles; k++)
	        {
	            if(WObj[k])
	               WlzFreeObj(WObj[k]);
	        }	


          exit(0);








   /* initial the outpoint  */
   if( ( surfP = (WlzDVertex3 *)AlcMalloc(sizeof(WlzDVertex3) * 100 )) == NULL ) 
   {
     errNum = WLZ_ERR_MEM_ALLOC;
   }



   /* allocate memory to store the tracked surface , in and out points. */
   if( ( SurfacePoints = ( WlzDVertex3 *)AlcMalloc(sizeof( WlzDVertex3) * ( numOf2DWlzFiles 
               * NumberToTrack  ) )) == NULL ) 
   {
      errNum = WLZ_ERR_MEM_ALLOC;
   }
   
   if( ( InPoints      = ( WlzDVertex3 *)AlcMalloc(sizeof( WlzDVertex3)
                               * (NumberToTrack * numOf2DWlzFiles ) )) == NULL ) 
   {
      errNum = WLZ_ERR_MEM_ALLOC;
   } 
   
   if( (OutPoints      = ( WlzDVertex3 *)AlcMalloc(sizeof( WlzDVertex3) 
                          * (NumberToTrack * numOf2DWlzFiles) )) == NULL ) 
   {
      errNum = WLZ_ERR_MEM_ALLOC;
   }

   if( (numOfPInS      = ( int *)AlcMalloc(sizeof( int) 
                          * numOf2DWlzFiles )) == NULL ) 
   {
      errNum = WLZ_ERR_MEM_ALLOC;
   }

   if(errNum == WLZ_ERR_NONE)
   {
       for(loopForFile=0; loopForFile<numOf2DWlzFiles; loopForFile++)
       {
          *(numOfPInS + loopForFile )  = 0;
       }
     
   }
  
   /* allocate memory to store the line sections of tracking data for each file */

   for(loopForFile=0; loopForFile<numOf2DWlzFiles; loopForFile++)
   {
     if( ( sectionData[loopForFile] = ( WlzIVertex2 *)AlcMalloc(sizeof( WlzIVertex2) *  subSubSectionLength_L  
                )) == NULL ) 
      {
         errNum = WLZ_ERR_MEM_ALLOC;
      }
      if(errNum != WLZ_ERR_NONE)
         break;
   }   

   test = 0;

  if( test && (errNum == WLZ_ERR_NONE) )
  {
    testConsistencyBreakDownTheGM(  gModel, numOf2DWlzFiles, sectionLength_N/2, &errNum );

    for(i=0; i<numOf2DWlzFiles; i++)
    {
          WlzFreeObj(WObj[i]);
    }
    exit(0);
  }

   /* initial the standard obj first: */
   if( errNum == WLZ_ERR_NONE)
   {

         if( errNum == WLZ_ERR_NONE)
         {
            /* get the number of Shells in this model */
            numOfShellsInStartGM  = NumberOfShellsAboutGM(gModel[ifileStar]);
	    printf("num of shells = %d\n", numOfShellsInStartGM);
	    if( endShell > numOfShellsInStartGM )
	    {
	        printf("not that many shells!\n");
	        for(k=0; k<numOf2DWlzFiles; k++)
	        {
	            if(WObj[k])
	               WlzFreeObj(WObj[k]);
	        }	
                exit(0);
	    }
            /* get the Index for each Shell; */
            startShellIdx = WlzShellsIndexAboutGM(gModel[ifileStar], numOfShellsInStartGM, &errNum);
         }
         zdvoxels = 0;			
   } 





   /* begin tracking */ 
   if( ( errNum == WLZ_ERR_NONE) && !test  )
   {
       if(startShell <0)
            startShell = 0;
       if(endShell < 0)
            endShell   = numOfShellsInStartGM;
	   
       /* cycle for each shell */
       for(loopForShell=startShell; loopForShell< endShell; loopForShell++)
        {
           /* get number loops in this shell */
           numOfLoopsInStartShell  =  NumberOfLoopsInThe_nThShellOfGM(gModel[ifileStar],  startShellIdx[loopForShell] );
	   printf("The  %d th Shell\n",  loopForShell );

           /* get the index for the loops */
           startLoopIdx            =  WlzLoopIndexAboutGM(gModel[ifileStar], startShellIdx[loopForShell], 
	                                      numOfLoopsInStartShell, &errNum);

	   if(errNum != WLZ_ERR_NONE)
	   {
              if( startLoopIdx  )
	            AlcFree(startLoopIdx);
	       loopForShell++;	 
	   }


           /* cycle for loops now as we only have simple shell, we will only use one loop for each shell
              for(loopForLoop=0; loopForLoop < numOfLoopsInStartShell; loopForLoop++) */
	   endLoop = numOfLoopsInStartShell;
	   for(loopForLoop=0; loopForLoop < endLoop; loopForLoop++)
	   {
              printf("The %d th loop\n",loopForLoop);

              if( errNum == WLZ_ERR_NONE)
              {
	        /* get number of vertices in this loop */
                numOfVerticesInTheStartLoop =   NumberOfVerticesInTheLoopOfGM(gModel[ifileStar], startLoopIdx[loopForLoop], &errNum);

	        if( ( numOfLoopsInStartShell == 2 ) && (loopForLoop == 0 ) )
	        {
		    if( IsACycle(gModel[ifileStar],  startLoopIdx[loopForLoop], &errNum ) )
		    {
	                     endLoop = 1;
		    }	     
		    else
		    {
                             numOfVerticesInTheStartLoop =  numOfVerticesInTheStartLoop/2 + 1;
		    }
	        }	 
	        else
	        {
                    numOfVerticesInTheStartLoop =  numOfVerticesInTheStartLoop/2 + 1;
	        }
		if(numOfVerticesInTheStartLoop < sectionLength_N/2)
		      break;
              }
      
              if( errNum == WLZ_ERR_NONE)
              {
                 /* get all the vertices in this loop (accurately it double the vertices!!) */
		 /*
                 AllVerticesInThisStartLoop  = WlzVerticesThisLoopOfGM( gModel[ifileStar], 
	                                                             startLoopIdx[loopForLoop], 
								     numOfVerticesInTheStartLoop, 
								    &errNum  );
								    */
                 AllVerticesInThisStartLoop  = WlzVerticesThisLoopOfSimpleShell2DGM( gModel[ifileStar], 
	                                                             startLoopIdx[loopForLoop], 
								    &numOfVerticesInTheStartLoop, 
								    &errNum  );
			   /* for test */					    
		           for(k=0; k<numOfVerticesInTheStartLoop; k++)
	   		   {
                			printf("%lg  %lg  %lg\n", ( AllVerticesInThisStartLoop.d2 + k )->vtX,
		                	     		          ( AllVerticesInThisStartLoop.d2 + k )->vtY ),
								  (double) 0;
	   		   }
						    
              }

              if( errNum == WLZ_ERR_NONE )
              {
                /*  sample into serveral sections */
	        if(!test)
	        {
	          numOfSections    =   numOfVerticesInTheStartLoop/sectionLength_N;
	          printf("There are %d number of vertices in This Shell\n",   numOfVerticesInTheStartLoop);
	       
	          printf("And it is divided into %d sections\n",  numOfSections );
	      
	          /*  cycle through the sections */
	          if( startSection < 0 )
	             startSection = 0;
	          if( endSection < 0 )
	             endSection = numOfSections;
	          for(loopForSection =startSection; loopForSection <endSection; loopForSection++)
	          {
	             nOfVOnTheS   = 0;
	             nOfVWithinTheS  = 0;
	             nOfVOutTheS = 0;

		     /*  we need to out put the section vertices to get a contour for later registration */
		     /*  and for check */
		     if(test)
		     {
		       printf("Start of This section: %d\n", loopForSection);
		  
		       outputThisSection(  AllVerticesInThisStartLoop, 
                                              loopForSection * sectionLength_N,
                                              subSubSectionLength_L,
                                             &errNum    );

		       printf("End of This section: %d\n", loopForSection);
		     }  
		     
                     /*  get the sampe points in the loopForSection th section for track down */
                     GetSamplePointsFromAllVertices( AllVerticesInThisStartLoop, 
                                                     loopForSection * sectionLength_N,
                                                     subSubSectionLength_L,
                                                     StandSampleP0,
                                                    &errNum );


		     /*  the again for future track up: */
                     GetSamplePointsFromAllVertices( AllVerticesInThisStartLoop, 
                                                     loopForSection * sectionLength_N,
                                                     subSubSectionLength_L,
                                                     StandSampleP,
                                                    &errNum );

		     suc = 1;	



		/* track down */
		if( (errNum == WLZ_ERR_NONE)  && suc )
		{
		       	 zdvoxels = -1;

         		 for(i= startTrackingFile-1; i >= startTrackingFile - numberOfFilesDownOrUp; i--)
          		 {     
				
	       			if( errNum == WLZ_ERR_NONE   )
       				{ 
                                    if(test)
         	        	         printf("------Target -------\n");

	  			     /*  get the tracked points in a loop and get the loop index and */
	 		 	     /*  number of sample points  */
	 		 		GetTrackedSamplePointsFromOneLoopOfGM(gModel[i], 
	                                        tTreeArray[i], 
					       &ntL, 
						numberOfPixelsZ,
	                                        minDis,
	                                        StandSampleP0,
						TrackedSamplePBlgToOneLp,
						&nOfTracked,
						&suc,
						&errNum);

              			}
                                
	      			if(!suc)
	      			{
                 			break; /*  break this section track */

	      			}

              			if(errNum == WLZ_ERR_NONE)
	      			{

	            			/*  get the sampled points and stored in the StandSampleP */
					/*  also output for contour and checking */
					if(test)
					   printf("begain target points %d th section\n", i);

					GetSamplePointsForNextTrackFromTheTrackedLoopOfGM(
					               i,
						       subSubSectionLength_L,
						       numOfPInS,
                                                       sectionData[i],
						       gModel[i],
                                                       ntL, 
                                                       StandSampleP0,
				                       TrackedSamplePBlgToOneLp,
				                       nOfTracked,
						      &suc,
                                                      &errNum );
                                        if(test)
					   printf("end target points %d th section\n", i);
      
				}		   

	              		if(errNum == WLZ_ERR_NONE)
	      			{
				        if(!suc)
	      				{
                 				   break; /*  break this section track */

	      				}

           	    			/*  store the output sample points: */
  	       	    			for(k=0; k<NumberToTrack; k++)
	            			{	
                   	  			SurfacePoints[nOfVOnTheS].vtX = StandSampleP0[k].vtX;
                   	  			SurfacePoints[nOfVOnTheS].vtY = StandSampleP0[k].vtY;
                   	  			SurfacePoints[nOfVOnTheS].vtZ = zdvoxels;
                                		nOfVOnTheS++;
	       	    			}

	   	    			/*  get the loop index we have just tracked */
           	    		        /*  calculate the in and out surface points  */

                                        GetInOutAndFromSurfacePoints( (int) zdvoxels, 
					                                StandSampleP0,
								        disForInAndOutGuid,
                                        	                        InPoints,   
								       &nOfVWithinTheS, 
					 	                        OutPoints,  
								       &nOfVOutTheS,
								       &errNum);


	      			  }
              			  /*  prepare for next section track: */
				  if(!suc)
				  {
                 				break; /*  break this section track */
				  }
				  /* if(newTObj) */
	      			  /*    WlzFreeObj(newTObj); */
                           	  zdvoxels -= (double) numberOfPixelsZ;


           		   }   /* end of track down the 2D image files array */
                }



                     /*  track up now: */

                     if( ( errNum == WLZ_ERR_NONE ) && ( suc ))
                     {
			 zdvoxels = 0;
			 /*  cycle through the 2D image files (tracking process) */


			 /*  track up first: */
         		 /* for(i=0; i<numOf2DWlzFiles; i++), track up first: */
         		 for(i= startTrackingFile; i<= startTrackingFile + numberOfFilesDownOrUp; i++)
          		 {    
	       			if( errNum == WLZ_ERR_NONE   )
       				{ 
				      
                                        if(test)
         	        		    printf("------Target -------\n");
	  				/*  get the tracked points in a loop and get the loop index and */
	 		 		/*  number of sample points  */
	 		 		GetTrackedSamplePointsFromOneLoopOfGM(gModel[i], 
	                                        tTreeArray[i], 
					       &ntL, 
						numberOfPixelsZ,
	                                        minDis,
	                                        StandSampleP,
						TrackedSamplePBlgToOneLp,
						&nOfTracked,
						&suc,
						&errNum);
	
              			}

	      			if(!suc)
	      			{
                 			break; /*  break this section track */

	      			}
              			if(errNum == WLZ_ERR_NONE)
	      			{

	            			/*  get the sampled points and stored in the StandSampleP */
					/*  also output for contour and checking */
					if(test)
					   printf("begain target points %d th section\n", i);
					/*
	     	    			GetTrackedPlusSomeSamplePointsFromOneLoopOfGM(
					               j, 
						       loopForSection, 
						       i,
						       subSubSectionLength_L,
						       numOfPInS,
                                                       sectionData[i],
						       gM2,
                                                       ntL, 
                                                       StandSampleP,
				                       TrackedSamplePBlgToOneLp,
				                       nOfTracked,
						      &suc,
                                                      &errNum );
                                        */

					GetSamplePointsForNextTrackFromTheTrackedLoopOfGM(
					               i,
						       subSubSectionLength_L,
						       numOfPInS,
                                                       sectionData[i],
						       gModel[i],
                                                       ntL, 
                                                       StandSampleP,
				                       TrackedSamplePBlgToOneLp,
				                       nOfTracked,
						      &suc,
                                                      &errNum );
                                        if(test)
					   printf("end target points %d th section\n", i);
      
				}		   

	              		if(errNum == WLZ_ERR_NONE)
	      			{

				        if(!suc)
	      				{
						/* if(newTObj) */
	         				  /*  WlzFreeObj(newTObj); */
                 				   break; /*  break this section track */

	      				}
	      

           	    			/*  store the output sample points: */
  	       	    			for(k=0; k<NumberToTrack; k++)
	            			{	
                   	  			SurfacePoints[nOfVOnTheS].vtX = StandSampleP[k].vtX;
                   	  			SurfacePoints[nOfVOnTheS].vtY = StandSampleP[k].vtY;
                   	  			SurfacePoints[nOfVOnTheS].vtZ = zdvoxels;
                                		nOfVOnTheS++;
	       	    			}

	   	    			/*  get the loop index we have just tracked */
           	    		        /*  calculate the in and out surface points  */

                                        GetInOutAndFromSurfacePoints( (int) zdvoxels, 
					                                StandSampleP,
								            disForInAndOutGuid,
                                        	                            InPoints,   
								     &nOfVWithinTheS, 
					 	                           OutPoints,  
								     &nOfVOutTheS,
								     &errNum);

	      			  }
              			  /*  prepare for next section track: */
				  if(!suc)
				  {
				                /*  if(newTObj) */
	         				 /*   WlzFreeObj(newTObj); */
                 				break; /*  break this section track */
				  }
				  /* if(newTObj) */
	      			    /*  WlzFreeObj(newTObj); */
                           	  zdvoxels += (double) numberOfPixelsZ;


           		   }   /* end of track up the 2D image files array */

       		 }


                if(suc)
	        {
                     /*  recalculate in and out points by using 3D  */
		     /*  (as they are needed in Bill's code which seems to */
		     /*   use the in and outside points to get normal for */
		     /*   the interpolate surface!!!) */

		     /* Calculate the IN and Out points by using 3D normal */
		     RecalInOutPBy3DNorm(  2*numberOfFilesDownOrUp+1,
		                           disForInAndOut,
		     			   SurfacePoints, 
					   nOfVOnTheS,
		     			   InPoints,      
				          &nOfVWithinTheS, 
		     			   OutPoints,     
				          &nOfVOutTheS
				         );
		
        	     /* output points */
		     /*
		     outputSurfaceInAndOutPoints(j, 
		                                 loopForLoop,
						 loopForSection,
                                                 2*numberOfFilesDownOrUp+1,
                             			 numOfPInS,
                             			 sectionData,
		     				 surfacefp, 
						 surfacePointFileName,
						 infp,  
						 surfaceInPointFileName,
						 outfp, 
						 surfaceOutPointFileName,
		                                 SurfacePoints, 
						 nOfVOnTheS,
						 InPoints,      
						 nOfVWithinTheS,
						 OutPoints,     
						 nOfVOutTheS,
						 &errNum
						 );
					*/	 
		    outputSurfaceByPatch(loopForShell, 
		                                 loopForLoop,
						 loopForSection,
                                                 2*numberOfFilesDownOrUp+1,
						 subSubSectionLength_L,
                             			 numOfPInS,
                             			 sectionData,
		     				 surfacefp, 
						 surfacePointFileName,
						 infp,  
						 surfaceInPointFileName,
						 outfp, 
						 surfaceOutPointFileName,
		                                 SurfacePoints, 
						 nOfVOnTheS,
						 InPoints,      
						 nOfVWithinTheS,
						 OutPoints,     
						 nOfVOutTheS,
						 &errNum
						 );



						 

	        }	 
             }/*  end of track each section */
          }
       } 
       /* prepare for the next loop tracking */
        AlcFree(AllVerticesInThisStartLoop.d2);
     } /* end of cycle for the loop */
       AlcFree(startLoopIdx);
    } /*  end of cycle for the shell */
     /* if(startLoopIdx) */

  }

 /*  if(SurfacePoints) */
     AlcFree(SurfacePoints);
  /* if(InPoints) */
     AlcFree(InPoints);
 /*  if(OutPoints) */
     AlcFree(OutPoints);

  for(i=0; i<numOf2DWlzFiles; i++)
  {
     AlcFree(tTreeArray[i]);
     AlcFree(sectionData[i]);
     WlzFreeObj(WObj[i]);
  }
  /* if(startShellIdx) */
     AlcFree(startShellIdx);

  /* if(numOfPInS) */
     AlcFree(numOfPInS);


  /* if(testFile) */
      fclose(testFile);
   exit(0);

   i = AlcKDTTreeFacts(tTree, testFile);

  return surfP; 
}



/*!
*  \return      Number of Vertices in the nThShell and nThLoop;
*  \ingroup     WlzFeatures
*  \brief       get Number of Vertices in the nThShell and nThLoop of
*               given model.
*  \param       gM                  given G model 
*  \param       nThShell            n-th shell in the model
*  \param       nThLoop             n-th loop of the n-th shell in the model
*/
static int nOfEdgeEndVInOneLoopOfGM(WlzGMModel *gM, int nThShell, int nThLoop)
{
   int nShell, nLoop, nVertices;
   WlzGMShell         *cS, *fS;
   WlzGMLoopT         *cLT, *fLT;
   WlzGMEdgeT         *cET, *fET;
  
   nShell    = 0;
   nLoop     = 0;
   nVertices = 0;
   
   /* For each shell of the model. */
    cS = fS = (WlzGMShell *) gM->child;
     do
     {
  	if ( nShell >= nThShell)
	   break;
        cS = cS->next;
        ++nShell;
     }  while( cS != fS  );
   /* For each loop topology element of the model. */
   
    cLT = fLT = (WlzGMLoopT *) cS->child;
     do
     {
	if (nLoop >= nThLoop )
	  break;
        cLT = cLT->next;
        ++nLoop;
     } while(  cLT != fLT );
    /* For each edge topology element of the model. */
    cET = fET = cLT->edgeT;
    do
    {
       /*
        if(cET == cET->edge->edgeT)
        {
	    nVertices += 2;
        }
	*/
	nVertices +=1;
	cET = cET->next;
    } while (cET != fET);

    return  nVertices;
}

/*!
*  \return      Number of unique Vertices in the nThShell and nThLoop;
*  \ingroup     WlzFeatures
*  \brief       get Number of Vertices in the nThShell and nThLoop of
*               given model.
*  \param       gM                  given G model 
*  \param       nThShell            n-th shell in the model
*  \param       nThLoop             n-th loop of the n-th shell in the model
*/
static int nOfVInOneLoopOfGM(WlzGMModel *gM, int nThShell, int nThLoop)
{
   int nShell, nLoop, nVertices;
   WlzGMShell         *cS, *fS;
   WlzGMLoopT         *cLT, *fLT;
   WlzGMEdgeT         *cET, *fET;
  
   nShell    = 0;
   nLoop     = 0;
   nVertices = 0;
   
   /* For each shell of the model. */
    cS = fS = (WlzGMShell *) gM->child;
     do
     {
  	if ( nShell >= nThShell)
	   break;
        cS = cS->next;
        ++nShell;
     }  while( cS != fS  );
   /* For each loop topology element of the model. */
   
    cLT = fLT = (WlzGMLoopT *) cS->child;
     do
     {
	if (nLoop >= nThLoop )
	  break;
        cLT = cLT->next;
        ++nLoop;
     } while(  cLT != fLT );
    /* For each edge topology element of the model. */
    cET = fET = cLT->edgeT;
    do
    {
	nVertices +=1;
	if( ( cET->vertexT->diskT->idx == cET->next->next->vertexT->diskT->idx ) )
	   break;
	cET = cET->next;
    } while (cET != fET);

    return  nVertices;
}



/*!
*  \return      none;
*  \ingroup     WlzFeatures
*  \brief       print out Number of shells, number of Loops and number of edge end vertices 
*               of the given G model 
*  \param       gM                  given G model 
*/
static void FactsAboutGM(WlzGMModel *gM)
{
   int nShell, nLoop, nVertices;
   WlzGMShell         *cS, *fS;
   WlzGMLoopT         *cLT, *fLT;
   WlzGMEdgeT         *cET, *fET;
  
   nShell    = 0;
   nLoop     = 0;
   nVertices = 0;

    /* For each shell of the model. */
    cS = fS = (WlzGMShell *) gM->child;
    do
    {
      printf("Shell:  %d  Shell idx %d\n",nShell, cS->idx);
      /* For each loop topology element of the model. */
       cLT = fLT = (WlzGMLoopT *) cS->child;
       /* printf("Next loop\n"); */
       do
       {
           printf("   Loop:  %d  Loop idx %d\n", nLoop, cLT->idx);
           /* For each edge topology element of the model. */
	   cET = fET = cLT->edgeT;
   	   do
	   {
              if(cET == cET->edge->edgeT) /*  Print edge end points  */
	      {
	         nVertices += 2;
	      }
	      cET = cET->next;
	   } while (cET != fET);
	   /*  if( nVertices > 50) */
	   printf("       number of edge end vertics in the Loop:  %d\n", nVertices);
	   nVertices = 0;
           cLT = cLT->next;
           nLoop++;
       } while(cLT != fLT);
       nLoop     = 0;
       cS = cS->next;
       ++nShell;
   }  while(cS != fS );
}

/*!
*  \return      none;
*  \ingroup     WlzFeatures
*  \brief       print out Number of shells, number of Loops and number of edge end vertices 
*               of the given G model 
*  \param       gM                  given G model 
*/
static int NumberOfShellsAboutGM(WlzGMModel *gM)
{
   int nShell;
   WlzGMShell         *cS, *fS;
  
   nShell    = 0;

    /* For each shell of the model. */
    cS = fS = (WlzGMShell *) gM->child;
    do
    {
       cS = cS->next;
       ++nShell;
    }  while(cS != fS );
   return nShell;
}

/*!
*  \return      none;
*  \ingroup     WlzFeatures
*  \brief       print out Number of shells, number of Loops and number of edge end vertices 
*               of the given G model 
*  \param       gM                  given G model 
*/
static int *WlzShellsIndexAboutGM(WlzGMModel *gM, int nS, WlzErrorNum *dstErr)
{
   int nShell;
   int *vData;
   WlzGMShell         *cS, *fS;
  
   nShell    = 0;
   
   if((vData = (int *)AlcMalloc(sizeof(int) * nS)) == NULL)
   {
     *dstErr = WLZ_ERR_MEM_ALLOC;
   }
   
   if( *dstErr == WLZ_ERR_NONE)
   {
      /* For each shell of the model. */
      cS = fS = (WlzGMShell *) gM->child;
      do
      {
         *(vData + nShell) = cS->idx;
          cS = cS->next;
        ++nShell;
      }  while(cS != fS );

    
   }
   return vData;
}


/*!
*  \return      none;
*  \ingroup     WlzFeatures
*  \brief       print out Number of shells, number of Loops and number of edge end vertices 
*               of the given G model 
*  \param       gM                  given G model 
*/
static int *WlzLoopIndexAboutGM(WlzGMModel *gM, int nThShell, int numLoops, WlzErrorNum *dstErr)
{
   int *vData;
   int nLoop;
   AlcVector          *vec;
   WlzGMShell         *cS;
   WlzGMLoopT         *cLT, *fLT;

   if(numLoops < 1)
   {
       
     if((vData = (int *)AlcMalloc(sizeof(int) * 1)) == NULL)
     {
        *dstErr = WLZ_ERR_MEM_ALLOC;
     }
     
     *vData = -19;
     return vData;
   }
  
   if((vData = (int *)AlcMalloc(sizeof(int) * numLoops)) == NULL)
   {
     *dstErr = WLZ_ERR_MEM_ALLOC;
   }
   
   if( *dstErr == WLZ_ERR_NONE)
   {

 
      vec    =  gM->res.shell.vec;
      nLoop     = 0;
      /* For each shell of the model. */
      cS = (WlzGMShell *) AlcVectorItemGet(vec, nThShell);
      cLT = fLT = (WlzGMLoopT *) cS->child;
       /* printf("Next loop\n"); */
       do
       {
           /* For each edge topology element of the model. */
           *(vData + nLoop) = cLT->idx;
            nLoop++;
            cLT = cLT->next;
       } while(cLT != fLT);
   }



   return vData;
}



/*!
* \return	Allocated edge end vertices.
* \ingroup	WlzFeatures
* \brief	Allocates a buffer which it fills with the vertices
*		from a geometric model.
* \param	model			Given geometric model.
* \param       nThShell                n-th shell in the model
* \param       nThLoop                 n-th loop of the n-th shell in the model
*					of vertices.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
static WlzVertexP WlzEdgeEndVerticesOneLoopFromGM( WlzGMModel  *model, 
                                            int          nThShell, 
					    int          nThLoop, 
					    WlzErrorNum *dstErr   )
{
 WlzVertexP vData;
 WlzErrorNum errNum = WLZ_ERR_NONE;

 vData.v = NULL;
 if( model != NULL)
 {
    switch(model->type)
    {
      case WLZ_GMMOD_2I: /* FALLTHROUGH */
      case WLZ_GMMOD_2D:
	vData.d2 = WlzDEdgeEndVerticesOneLoopFromGM2(model, nThShell, nThLoop, &errNum);
        break;
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D:
	vData.d3 = WlzDEdgeEndVerticesOneLoopFromGM3(model, nThShell, nThLoop, &errNum);
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
 }
 if(dstErr)
 {
    *dstErr = errNum;
 }
 return(vData);

}


/*!
* \return	Allocated vertices.
* \ingroup	WlzFeatures
* \brief	Allocates a buffer which it fills with the vertices
*		from a geometric model.
* \param	model			Given geometric model.
* \param       nThShell                n-th shell in the model
* \param       nThLoop                 n-th loop of the n-th shell in the model
*					of vertices.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
static WlzVertexP WlzVerticesOneLoopFromGM( WlzGMModel  *model, 
                                            int          nThShell, 
					    int          nThLoop, 
					    WlzErrorNum *dstErr   )
{
 WlzVertexP vData;
 WlzErrorNum errNum = WLZ_ERR_NONE;

 vData.v = NULL;
 if( model != NULL)
 {
    switch(model->type)
    {
      case WLZ_GMMOD_2I: /* FALLTHROUGH */
      case WLZ_GMMOD_2D:
	vData.d2 = WlzDVerticesOneLoopFromGM2(model, nThShell, nThLoop, &errNum);
        break;
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D:
	vData.d3 = WlzDVerticesOneLoopFromGM3(model, nThShell, nThLoop, &errNum);
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
 }
 if(dstErr)
 {
    *dstErr = errNum;
 }
 return(vData);

}

/*!
* \return	Allocated edge end vertices.
* \ingroup	WlzFeatures
* \brief	Allocates a buffer which it fills with the vertices
*		from a 2D GM.
* \param	model			Given model.
* \param	dstCnt			Destination ptr for the number
*					of vertices.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
static WlzDVertex2 *WlzDEdgeEndVerticesOneLoopFromGM2(WlzGMModel *model, 
                                            int            nThShell, 
					    int            nThLoop, 
					    WlzErrorNum   *dstErr   )
{
  int		idx,
  		cnt,
		vIdx;
  int      nShell, nLoop, nVertices;		
   WlzGMShell         *cS, *fS;
   WlzGMLoopT         *cLT, *fLT;
   WlzGMEdgeT         *cET, *fET;
   WlzVertexType	type;
   WlzDVertex2   *vData = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  vIdx = 0;

  cnt = nOfEdgeEndVInOneLoopOfGM(model, nThShell, nThLoop);
  if((vData = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) * cnt)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /*  go to the correction position */

 
   nShell    = 0;
   nLoop     = 0;
   nVertices = 0;
   
   /* For each shell of the model. */
    cS = fS = (WlzGMShell *) model->child;
     do
     {
	if ( nShell >= nThShell)
	   break;
        cS = cS->next;
        ++nShell;
     }  while( cS != fS  );
   /* For each loop topology element of the model. */
   
    cLT = fLT = (WlzGMLoopT *) cS->child;
     do
     {
	if (nLoop >= nThLoop )
	  break;
        cLT = cLT->next;
        ++nLoop;
     } while(  cLT != fLT );
    /* For each edge topology element of the model. */
    cET = fET = cLT->edgeT;
    do
    {
        if(cET == cET->edge->edgeT) /*  Print edge end points  */
        {

	    if(model->type == WLZ_GMMOD_2I)
	    {
	       (vData + nVertices)->vtX     = (double) cET->vertexT->diskT->vertex->geo.vg2I->vtx.vtX;
	       (vData + nVertices)->vtY     = (double) cET->vertexT->diskT->vertex->geo.vg2I->vtx.vtY;
	       (vData + nVertices + 1)->vtX = (double) cET->opp->vertexT->diskT->vertex->geo.vg2I->vtx.vtX;
	       (vData + nVertices + 1)->vtY = (double) cET->opp->vertexT->diskT->vertex->geo.vg2I->vtx.vtY;
               nVertices += 2;
	    }
            else
            {
	       (vData + nVertices)->vtX     = cET->vertexT->diskT->vertex->geo.vg2D->vtx.vtX;
	       (vData + nVertices)->vtY     = cET->vertexT->diskT->vertex->geo.vg2D->vtx.vtY;
	       (vData + nVertices + 1)->vtX = cET->opp->vertexT->diskT->vertex->geo.vg2D->vtx.vtX;
	       (vData + nVertices + 1)->vtY = cET->opp->vertexT->diskT->vertex->geo.vg2D->vtx.vtY;
               nVertices += 2;
	    }
        }
	cET = cET->next;
    } while (cET != fET);

  if(dstErr)
  {
    *dstErr = errNum;
  }
  } 
  return(vData);
}






/*!
* \return	Allocated vertices.
* \ingroup	WlzFeatures
* \brief	Allocates a buffer which it fills with the vertices
*		from a 2D GM.
* \param	model			Given model.
* \param	dstCnt			Destination ptr for the number
*					of vertices.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
static WlzDVertex2 *WlzDVerticesOneLoopFromGM2(WlzGMModel *model, 
                                            int            nThShell, 
					    int            nThLoop, 
					    WlzErrorNum   *dstErr   )
{
  int		idx,
  		cnt,
		vIdx;
  int      nShell, nLoop, nVertices;		
   WlzGMShell         *cS, *fS;
   WlzGMLoopT         *cLT, *fLT;
   WlzGMEdgeT         *cET, *fET;
   WlzVertexType	type;
   WlzDVertex2    *vData = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  vIdx = 0;

  cnt = nOfVInOneLoopOfGM(model, nThShell, nThLoop);
  if((vData = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) * cnt)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /*  go to the correction position */

 
   nShell    = 0;
   nLoop     = 0;
   nVertices = 0;
   
   /* For each shell of the model. */
    cS = fS = (WlzGMShell *) model->child;
     do
     {
	if ( nShell >= nThShell)
	   break;
        cS = cS->next;
        ++nShell;
     }  while( cS != fS  );
   /* For each loop topology element of the model. */
   
    cLT = fLT = (WlzGMLoopT *) cS->child;
     do
     {
	if (nLoop >= nThLoop )
	  break;
        cLT = cLT->next;
        ++nLoop;
     } while(  cLT != fLT );
    /* For each edge topology element of the model. */
    cET = fET = cLT->edgeT;
    do
    {
        /*  if(cET == cET->edge->edgeT)  Print end points  */
        {

	    if(model->type == WLZ_GMMOD_2I)
	    {
	       (vData + nVertices)->vtX     = (double) cET->vertexT->diskT->vertex->geo.vg2I->vtx.vtX;
	       (vData + nVertices)->vtY     = (double) cET->vertexT->diskT->vertex->geo.vg2I->vtx.vtY;
               nVertices += 1;
	    }
            else
            {
	       (vData + nVertices)->vtX     = cET->vertexT->diskT->vertex->geo.vg2D->vtx.vtX;
	       (vData + nVertices)->vtY     = cET->vertexT->diskT->vertex->geo.vg2D->vtx.vtY;
               nVertices += 1;
	    }
	   if( ( cET->vertexT->diskT->idx == cET->next->next->vertexT->diskT->idx ) )
	      break;
        }
	cET = cET->next;
    } while (cET != fET);

  if(dstErr)
  {
    *dstErr = errNum;
  }
  } 
  return(vData);
}


/*!
* \return	Allocated edge end vertices.
* \ingroup	WlzFeatures
* \brief	Allocates a buffer which it fills with the vertices
*		from a 3D GM.
* \param	model			Given model.
* \param	dstCnt			Destination ptr for the number
*					of vertices.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
static WlzDVertex3 *WlzDEdgeEndVerticesOneLoopFromGM3( WlzGMModel  *model, 
                                                int          nThShell, 
					        int          nThLoop, 
					        WlzErrorNum *dstErr   )
{
  int		idx,
  		cnt,
		vIdx;
  int      nShell, nLoop, nVertices;	
  WlzGMShell         *cS, *fS;
  WlzGMLoopT         *cLT, *fLT;
  WlzGMEdgeT         *cET, *fET;
  WlzVertexType	type;
  WlzDVertex3   *vData = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  vIdx = 0;
  cnt = nOfEdgeEndVInOneLoopOfGM(model,  nThShell, nThLoop );
  if((vData = (WlzDVertex3 *)AlcMalloc(sizeof(WlzDVertex3) * cnt)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(model->type == WLZ_GMMOD_3I)
    {
    /*  go to the correction position */

 
   nShell    = 0;
   nLoop     = 0;
   nVertices = 0;
   
   /* For each shell of the model. */
    cS = fS = (WlzGMShell *) model->child;
     do
     {
	if ( nShell >= nThShell)
	   break;
        cS = cS->next;
        ++nShell;
     }  while( cS != fS  );
   /* For each loop topology element of the model. */
   
    cLT = fLT = (WlzGMLoopT *) cS->child;
     do
     {
	if (nLoop >= nThLoop )
	  break;
        cLT = cLT->next;
        ++nLoop;
     } while(  cLT != fLT );
    /* For each edge topology element of the model. */
    cET = fET = cLT->edgeT;
    do
    {
        if(cET == cET->edge->edgeT) /*  Print edge end points  */
        {

	    if(model->type == WLZ_GMMOD_3I)
	    {
	       (vData + nVertices)->vtX     = (double)  cET->vertexT->diskT->vertex->geo.vg3I->vtx.vtX;
	       (vData + nVertices)->vtY     = (double)  cET->vertexT->diskT->vertex->geo.vg3I->vtx.vtY;
	       (vData + nVertices)->vtZ     = (double)  cET->vertexT->diskT->vertex->geo.vg3I->vtx.vtZ;
	       (vData + nVertices + 1)->vtX = (double)  cET->opp->vertexT->diskT->vertex->geo.vg3I->vtx.vtX;
	       (vData + nVertices + 1)->vtY = (double)  cET->opp->vertexT->diskT->vertex->geo.vg3I->vtx.vtY;
	       (vData + nVertices + 1)->vtZ = (double)  cET->opp->vertexT->diskT->vertex->geo.vg3I->vtx.vtZ;
                nVertices += 2;
	    }
            else
            {
	       (vData + nVertices)->vtX     = cET->vertexT->diskT->vertex->geo.vg3D->vtx.vtX;
	       (vData + nVertices)->vtY     = cET->vertexT->diskT->vertex->geo.vg3D->vtx.vtY;
	       (vData + nVertices)->vtZ     = cET->vertexT->diskT->vertex->geo.vg3D->vtx.vtZ;
	       (vData + nVertices + 1)->vtX = cET->opp->vertexT->diskT->vertex->geo.vg3D->vtx.vtX;
	       (vData + nVertices + 1)->vtY = cET->opp->vertexT->diskT->vertex->geo.vg3D->vtx.vtY;
	       (vData + nVertices + 1)->vtZ = cET->opp->vertexT->diskT->vertex->geo.vg3D->vtx.vtZ;
                nVertices += 2;
	    }
        }
	cET = cET->next;
    } while (cET != fET);

  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
 } 
  return(vData);
}




/*!
* \return	Allocated vertices.
* \ingroup	WlzFeatures
* \brief	Allocates a buffer which it fills with the vertices
*		from a 3D GM.
* \param	model			Given model.
* \param	dstCnt			Destination ptr for the number
*					of vertices.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
static WlzDVertex3 *WlzDVerticesOneLoopFromGM3( WlzGMModel  *model, 
                                                int          nThShell, 
					        int          nThLoop, 
					        WlzErrorNum *dstErr   )
{
  int		idx,
  		cnt,
		vIdx;
  int      nShell, nLoop, nVertices;	
  WlzGMShell         *cS, *fS;
  WlzGMLoopT         *cLT, *fLT;
  WlzGMEdgeT         *cET, *fET;
  WlzVertexType	type;
  WlzDVertex3  	*vData = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  vIdx = 0;
  cnt = nOfVInOneLoopOfGM(model,  nThShell, nThLoop );
  if((vData = (WlzDVertex3 *)AlcMalloc(sizeof(WlzDVertex3) * cnt)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(model->type == WLZ_GMMOD_3I)
    {
    /*  go to the correction position */

 
   nShell    = 0;
   nLoop     = 0;
   nVertices = 0;
   
   /* For each shell of the model. */
    cS = fS = (WlzGMShell *) model->child;
     do
     {
	if ( nShell >= nThShell)
	   break;
        cS = cS->next;
        ++nShell;
     }  while( cS != fS  );
   /* For each loop topology element of the model. */
   
    cLT = fLT = (WlzGMLoopT *) cS->child;
     do
     {
	if (nLoop >= nThLoop )
	  break;
        cLT = cLT->next;
        ++nLoop;
     } while(  cLT != fLT );
    /* For each edge topology element of the model. */
    cET = fET = cLT->edgeT;
    do
    {
        /* if(cET == cET->edge->edgeT)  Print edge end points  */
        {

	    if(model->type == WLZ_GMMOD_3I)
	    {
	       (vData + nVertices)->vtX     = (double)  cET->vertexT->diskT->vertex->geo.vg3I->vtx.vtX;
	       (vData + nVertices)->vtY     = (double)  cET->vertexT->diskT->vertex->geo.vg3I->vtx.vtY;
	       (vData + nVertices)->vtZ     = (double)  cET->vertexT->diskT->vertex->geo.vg3I->vtx.vtZ;
                nVertices += 1;
	    }
            else
            {
	       (vData + nVertices)->vtX     = cET->vertexT->diskT->vertex->geo.vg3D->vtx.vtX;
	       (vData + nVertices)->vtY     = cET->vertexT->diskT->vertex->geo.vg3D->vtx.vtY;
	       (vData + nVertices)->vtZ     = cET->vertexT->diskT->vertex->geo.vg3D->vtx.vtZ;
                nVertices += 1;
	    }


        }
	if( ( cET->vertexT->diskT->idx == cET->next->next->vertexT->diskT->idx ) )
	      break;

	cET = cET->next;
    } while (cET != fET);

  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
 } 
  return(vData);
}

static void WlzGetTheShellAndLoopNumberFromTheIndexOfTheVertex( WlzGMModel *gM,
                                                          int idx, 
					                  int *nShell, 
					                  int *nLoop,
					                  WlzErrorNum *dstErr)
{
   AlcVector   *vec;
   WlzGMVertex *cVT;
   WlzGMLoopT  *cLT;
   WlzGMShell  *cST;
   vec    =  gM->res.vertex.vec;
   cVT    =   (WlzGMVertex *)AlcVectorItemGet(vec, idx);
   /*  printf("here\n"); */
   cLT     =   (WlzGMLoopT *) cVT->diskT->vertexT->parent->parent;
  *nLoop   =    cLT->idx;
   cST     =    cLT->parent;
   *nShell =    cST->idx;
   /*  printf("Shell:  %d  loop:  %d\n", *nShell, *nLoop); */
}

static int NumberOfLoopsInThe_nThShellOfGM(WlzGMModel *gM, int ShellIdx )
{
   int nLoop;
   AlcVector          *vec;
   WlzGMShell         *cS;
   WlzGMLoopT         *cLT, *fLT;
   
   vec    =  gM->res.shell.vec;

   nLoop     = 0;
    /* For each shell of the model. */
    cS = (WlzGMShell *) AlcVectorItemGet(vec, ShellIdx);
    if(!cS->child)
       return nLoop;
    cLT = fLT = (WlzGMLoopT *) cS->child;
       /* printf("Next loop\n"); */
       do
       {
           /* For each edge topology element of the model. */
           nLoop++;
           cLT = cLT->next;
       } while(cLT != fLT);

       return nLoop;
}


static int NumberOfVerticesInTheLoopOfGM(WlzGMModel   *gM, 
 						int sLoopIdx, 
						WlzErrorNum *dstErr)
{
     WlzErrorNum   errNum = WLZ_ERR_NONE;
   int nVertices;
   AlcVector          *vec;
   WlzGMLoopT         *cLT;
   WlzGMEdgeT         *cET, *fET;
   

   
   vec    =  gM->res.loopT.vec;
   nVertices     = 0;
   if(sLoopIdx < 1)
   {
      return  nVertices;
   }
    /* For each LoopT of the model. */
    cLT = (WlzGMLoopT *) AlcVectorItemGet(vec, sLoopIdx);
    cET = fET = cLT->edgeT;
    do
    {

        nVertices += 1;
	cET = cET->next;
    } while (cET != fET);
   
   *dstErr = errNum;

    return nVertices;

   
}  




/*!
* \return	Allocated vertices(note it gives double the vertices number in non-cycle shell).
*                  As the start point could be anywhere in the loop, we have no other choise
*                  otherwize we will loose some of the vertices and double counted for others!
* \ingroup	WlzFeatures
* \brief	Allocates a buffer which it fills with the vertices
*		from a 3D GM.
* \param	model			Given model.
* \param	dstCnt			Destination ptr for the number
*					of vertices.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
static int NumberOfVerticesDoubleInNonCycleCaseInTheLoopOfGM(WlzGMModel   *gM, 
 						int sLoopIdx, 
						WlzErrorNum *dstErr)
{
   WlzErrorNum   errNum = WLZ_ERR_NONE;
   int nVertices;
   AlcVector          *vec;
   WlzGMLoopT         *cLT;
   WlzGMEdgeT         *cET, *fET;
   
   vec    =  gM->res.loopT.vec;
   nVertices     = 0;
    /* For each LoopT of the model. */
    cLT = (WlzGMLoopT *) AlcVectorItemGet(vec, sLoopIdx);
    cET = fET = cLT->edgeT;
    do
    {

        nVertices += 1;
	cET = cET->next;
    } while (cET != fET);
    
    return ( nVertices );

}  



/*!
* \return	Allocated vertices.
* \ingroup	WlzFeatures
* \brief	Allocates a buffer which it fills with the vertices
*		from a geometric model.
* \param	model			Given geometric model.
* \param        LoopIdx                 The loop index in the model
* \param        numOfVertices           number of vertices in this loop of the model
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
static WlzVertexP WlzVerticesThisLoopOfGM(  WlzGMModel  *model, 
					    int          LoopIdx,
					    int          numOfVertices,
					    WlzErrorNum *dstErr   )
{
 WlzVertexP vData;
 WlzErrorNum errNum = WLZ_ERR_NONE;

 vData.v = NULL;
 if( model != NULL)
 {
    switch(model->type)
    {
      case WLZ_GMMOD_2I: /* FALLTHROUGH */
      case WLZ_GMMOD_2D:
	vData.d2 = WlzDVerticesThisLoopOfGM2(model, LoopIdx, numOfVertices, &errNum);
        break;
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D:
	vData.d3 = WlzDVerticesThisLoopOfGM3(model, LoopIdx, numOfVertices, &errNum);
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
 }
 if(dstErr)
 {
    *dstErr = errNum;
 }
 return(vData);

}


/*!
* \return	Allocated vertices.
* \ingroup	WlzFeatures
* \brief	Allocates a buffer which it fills with the vertices
*		from a 2D GM.
* \param	model			Given geometric model.
* \param        LoopIdx                 The loop index in the model
* \param        numOfVertices           number of vertices in this loop of the model
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
static WlzDVertex2 *WlzDVerticesThisLoopOfGM2(WlzGMModel  *model, 
                                            int            LoopIdx, 
					    int            numOfVertices, 
					    WlzErrorNum   *dstErr   )
{
   WlzDVertex2   *vData = NULL;
   WlzErrorNum   errNum = WLZ_ERR_NONE;
   int nVertices;
   AlcVector          *vec;
   WlzGMLoopT         *cLT;
   WlzGMEdgeT         *cET, *fET;
   int test = 0;
   
  if((vData = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) * numOfVertices)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
   
   vec    =  model->res.loopT.vec;
   nVertices     = 0;
    /* get the loopT of the model. */
    cLT = (WlzGMLoopT *) AlcVectorItemGet(vec, LoopIdx);
    /* For each edge topology element of the model. */    
    cET = fET = cLT->edgeT;
    if(test)
        printf("test loopIdx %d\n ",LoopIdx );
    do
    {
        if(nVertices >= numOfVertices )
	   break;

	if(model->type == WLZ_GMMOD_2I)
	{
	    (vData + nVertices)->vtX     = (double) cET->vertexT->diskT->vertex->geo.vg2I->vtx.vtX;
	    (vData + nVertices)->vtY     = (double) cET->vertexT->diskT->vertex->geo.vg2I->vtx.vtY;
            nVertices += 1;
	}
        else
        {

	    (vData + nVertices)->vtX     = cET->vertexT->diskT->vertex->geo.vg2D->vtx.vtX;
	    (vData + nVertices)->vtY     = cET->vertexT->diskT->vertex->geo.vg2D->vtx.vtY;
	    if(test && nVertices < 200)
	    {
	      printf(" idx %d  x %lg   y %lg\n", cET->vertexT->diskT->idx,
	                                         (vData + nVertices)->vtX,
	                                         (vData + nVertices)->vtY );
	    } 
	    nVertices += 1;

	}
	/* if( ( cET->vertexT->diskT->idx == cET->next->next->vertexT->diskT->idx ) ) */
	/*       break; */
	cET = cET->next;
    } while (cET != fET  );

  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(vData);
}

/*!
* \return	Allocated vertices.
* \ingroup	WlzFeatures
* \brief	Allocates a buffer which it fills with the vertices
*		from a 2D GM.
* \param	model			Given geometric model.
* \param        LoopIdx                 The loop index in the model
* \param        numOfVertices           number of vertices in this loop of the model
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
static WlzDVertex3 *WlzDVerticesThisLoopOfGM3(WlzGMModel  *model, 
                                            int            LoopIdx, 
					    int            numOfVertices, 
					    WlzErrorNum   *dstErr   )
{
   WlzDVertex3  *vData = NULL;
   WlzErrorNum   errNum = WLZ_ERR_NONE;
   int nVertices;
   AlcVector          *vec;
   WlzGMLoopT         *cLT;
   WlzGMEdgeT         *cET, *fET;
   
  if((vData = (WlzDVertex3 *)AlcMalloc(sizeof(WlzDVertex3) * numOfVertices)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
   
   vec    =  model->res.loopT.vec;
   nVertices     = 0;
    /* get the loopT of the model. */
    cLT = (WlzGMLoopT *) AlcVectorItemGet(vec, LoopIdx);
    /* For each edge topology element of the model. */    
    cET = fET = cLT->edgeT;
    do
    {
        if(nVertices >= numOfVertices )
	   break;

	if(model->type == WLZ_GMMOD_2I)
	{
	    (vData + nVertices)->vtX   = (double) cET->vertexT->diskT->vertex->geo.vg3I->vtx.vtX;
	    (vData + nVertices)->vtY   = (double) cET->vertexT->diskT->vertex->geo.vg3I->vtx.vtY;
	    (vData + nVertices)->vtZ   = (double) cET->vertexT->diskT->vertex->geo.vg3I->vtx.vtZ;
             nVertices += 1;
	}
        else
        {
	    (vData + nVertices)->vtX   = cET->vertexT->diskT->vertex->geo.vg3D->vtx.vtX;
	    (vData + nVertices)->vtY   = cET->vertexT->diskT->vertex->geo.vg3D->vtx.vtY;
	    (vData + nVertices)->vtZ   = cET->vertexT->diskT->vertex->geo.vg3D->vtx.vtZ;
             nVertices += 1;
	}
        
	/* if( ( cET->vertexT->diskT->idx == cET->next->next->vertexT->diskT->idx ) ) */
	/*       break; */
	cET = cET->next;
    } while (cET != fET);

  }

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(vData);
}

static WlzDVertex2 *WlzVerticesThisIndexOfGM2( WlzGMModel    *model, 
                                              int            Idx, 
					      WlzErrorNum   *dstErr   )
{
   WlzDVertex2    *vData;
   WlzErrorNum   errNum = WLZ_ERR_NONE;
   AlcVector          *vec;
   vec    =  model->res.vertex.vec;
    /* get the loopT of the model. */
    vData = (WlzDVertex2 *) AlcVectorItemGet(vec, Idx);
    return vData;

}

static void  WlzOutputForTestGM(WlzGMModel *gM, WlzErrorNum   *dstErr)
{
    int          i, j;
    int          nS, nL, nV;
    int         *sShellIdx = NULL, *sLoopIdx = NULL;
    WlzVertexP   sLoopVData;
    WlzErrorNum   errNum = WLZ_ERR_NONE;
    /*  get the number of Shells in this model */
     nS        = NumberOfShellsAboutGM(gM);

    /*  get the Index for each Shell; */
     sShellIdx = WlzShellsIndexAboutGM(gM, nS, &errNum);

     printf("nS = %d\n",nS);
     for(i=0; i< nS; i++)
     {
         printf("           n-th shell = %d  Shell-index %d\n",i, sShellIdx[i]);
         nL       = NumberOfLoopsInThe_nThShellOfGM(gM,  sShellIdx[i]);
         sLoopIdx = WlzLoopIndexAboutGM(gM, sShellIdx[i], nL, &errNum);
              
         for(j=0; j<nL; j++)
         {
             printf("       n-th loop = %d  Loop-index %d\n", j,  sLoopIdx[j] );
             nV       = NumberOfVerticesInTheLoopOfGM(gM, sLoopIdx[j], &errNum); 
             printf("       num Vertices = %d\n",nV);
             sLoopVData = WlzVerticesThisLoopOfGM(gM, sLoopIdx[j], nV, &errNum  );
	     /* if(sLoopVData.v) */
                AlcFree(sLoopVData.v);					   
         }
	 /* if(sLoopIdx) */
	 {
           AlcFree(sLoopIdx);
	   sLoopIdx = NULL;
	 }  
     }
     /* if(sShellIdx) */
     {
        
        AlcFree(sShellIdx);
	sShellIdx = NULL;
     }	
     if(dstErr)
       *dstErr = errNum;
}


static void  WlzOutputForTestAsectionData(WlzVertexP   sLoopVData, int nV, WlzErrorNum *dstErr)
{
     WlzErrorNum   errNum = WLZ_ERR_NONE;
     int i;
     if( sLoopVData.d2  == NULL)
     {
         errNum = WLZ_ERR_GMELM_NULL;
	 *dstErr = errNum;
     }	 
     if(errNum == WLZ_ERR_NONE)
     {
       printf(" === nT = %d  ===\n",nV);

       for(i=0; i < nV; ++i)
       {
         printf("%d %lg %lg\n",i, (sLoopVData.d2 + i)->vtX, (sLoopVData.d2 + i)->vtY );
       }
     }


}

/*!
* \return	none.
* \ingroup	WlzFeatures
* \brief	Output the curve for surfaces
* \param	sLoopVData		Given points in a curve of a plane.
* \param        nV                      number of vertices in the curve
* \param        zC                      z-coordinate to be addied.
* \param        sampleD                 sample distance used to select points to be output.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
static void  WlzOutputForFormingSurface(FILE *testFile, char *Fstr, WlzVertexP sLoopVData, int nV, double zC, int sampleD, WlzErrorNum *dstErr)
{
     WlzErrorNum   errNum = WLZ_ERR_NONE;
     int i;
     if( sLoopVData.d2  == NULL)
     {
         errNum = WLZ_ERR_GMELM_NULL;
	 *dstErr = errNum;
     }	 
     if(errNum == WLZ_ERR_NONE)
     {
        if((testFile = fopen(Fstr, "w")) == NULL )
        {
	     errNum =  WLZ_ERR_FILE_OPEN;
        }
        if(errNum == WLZ_ERR_NONE)
	{
           for(i=0; i < nV; ++i)
           {
              if(i%sampleD == 0)
                 fprintf(testFile, "%lg %lg %lg\n", (sLoopVData.d2 + i)->vtX, (sLoopVData.d2 + i)->vtY,  zC );
           }
	}
        fclose(testFile);
	testFile = NULL;
     }
}



/*!
* \ingroup	WlzFeatures
* \return	<void>
* \brief	Computes the normal (n) to a segment (g) between the
*		given pair of vertices. There are clearly two solutions
*		to the problem of finding a normal to a line segment,
*		but this function always finds the normal vector with a
*		+ve x component.
*		If the two vertices are coincident then the normal
*		vector is set to {0, 0}.
*		With two non-coincident vertices the normal vector is
*		computed using the relationships g.n = 0 and \|n\|^2 = 1.
*		Giving:
*		\f$ nx = \frac{1}{\sqrt(1+(gx/gy)^2)} \f$ ,
*		\f$ ny = -nx\frac{gx}{gy} \f$
*		There is no need for any type checking in this function
*		because it is static and all types have been checked.
* \param	v0			First of the given pair.
* \param	v1			Second of the given pair.
*/
static WlzDVertex2 WlzVerticesNormPair2(WlzDVertex2 v0, WlzDVertex2 v1)
{
  WlzDVertex2	tV0,
  		tV1,
		nrm;

  WLZ_VTX_2_SUB(tV0, v1, v0);
  tV1.vtX = tV0.vtX * tV0.vtX;
  tV1.vtY = tV0.vtY * tV0.vtY; 
  if(tV1.vtY < DBL_EPSILON)
  {
    nrm.vtX = 0.0;
    if(tV1.vtX < DBL_EPSILON)
    {
      nrm.vtY = 0.0;
    }
    else
    {
      nrm.vtY = 1.0;
    }
  }
  else if(tV1.vtX < DBL_EPSILON)
  {
    nrm.vtX = 1.0;
    nrm.vtY = 0.0;
  }
  else
  {
    nrm.vtX = 1.0 / sqrt(1.0 + (tV1.vtX / tV1.vtY));
    nrm.vtY = -((tV0.vtX) * nrm.vtX) / tV0.vtY;
  }
  return(nrm);
}



/*!
* \ingroup      WlzFeatures
* \return				Normal vector.
* \brief	Computes the normal (n) at a vertex. This is chosen to
*		be the unit vector which bisects the angle which two
*		line segments make at the vertex.
*
*		Given two line segments specified by three vertices
*		\f$A\f$, \f$B\f$ and \f$C\f$, with a common vertex
*		\f$B\f$. Find a pair of points \f$A'\f$ and \f$C'\f$
*		on line segments \f$B \rightarrow A\f$ and
*		\f$B \rightarrow C\f$ such that they have unit distance
*		from \f$B\f$ and are in the directions of \f$A\f$ and
*		\f$C\f$.  Next find the midpoint of the two vertices
*		\f$A'\f$ and \f$C'\f$, call this point \f$D\f$.
*		Lastly find the unit vector directed from \f$B\f$ towards
*		\f$D\f$.
*		If all three vertices are coincident a zero vector is
*		retuned.
* \param	vA			First vertex.
* \param	vB			Second vertex, common to both line
*					segments.
* \param	vC			Third vertex.
*/
static WlzDVertex2 WlzVerticesNormTriple2(WlzDVertex2 vA, WlzDVertex2 vB,
					   WlzDVertex2 vC)
{
  double	tD0;
  WlzDVertex2	tV0,
		tV1,
		tV2,
		tV3,
  		vAU,
  		vCU,
		vD,
		nrm;

  WLZ_VTX_2_SUB(tV0, vA, vB);
  WLZ_VTX_2_SUB(tV1, vC, vB);
  tV2.vtX = tV0.vtX * tV0.vtX; tV2.vtY = tV0.vtY * tV0.vtY;
  tV3.vtX = tV1.vtX * tV1.vtX; tV3.vtY = tV1.vtY * tV1.vtY;
  if((tV2.vtX < DBL_EPSILON) && (tV2.vtY < DBL_EPSILON))
  {
    nrm = WlzVerticesNormPair2(vB, vC);
  }
  else if((tV3.vtX < DBL_EPSILON) && (tV3.vtY < DBL_EPSILON))
  {
    nrm = WlzVerticesNormPair2(vB, vA);
  }
  else
  {
    /* Check for colinearity and coincidence of all three vertices by
     * computing the area of the triangle ABC. */
    tD0 = WlzGeomTriangleSnArea2(vA, vB, vC);
    if((tD0 * tD0) < (DBL_EPSILON))
    {
      nrm = WlzVerticesNormPair2(vB, vC);
    }
    else
    {
      /* Compute the positions of A' and C' */
      tD0 = 1.0 / sqrt(tV2.vtX + tV2.vtY);
      WLZ_VTX_2_SCALE(vAU, tV0, tD0);
      WLZ_VTX_2_ADD(vAU, vAU, vB);
      tD0 = 1.0 / sqrt(tV3.vtX + tV3.vtY);
      WLZ_VTX_2_SCALE(vCU, tV1, tD0);
      WLZ_VTX_2_ADD(vCU, vCU, vB);
      /* Find D, the midpoint between A' and C' */
      WLZ_VTX_2_ADD(vD, vAU, vCU);
      WLZ_VTX_2_SCALE(vD, vD, 0.5);
      /* Compute the unit normal vector. */
      WLZ_VTX_2_SUB(nrm, vD, vB);
      tD0 = 1.0 / (WLZ_VTX_2_LENGTH(nrm));
      WLZ_VTX_2_SCALE(nrm, nrm, tD0);
    }
  }
  return(nrm);
}



/*!
* \ingroup      WlzFeatures
* \return			array of int	.
* \brief	Computes the number of elements in each group and its
*               continutity. 
*
* \param	LoopIndex[]		input index of loops
* \param       *NumOfGroup              output
* \param       *Continue                output
* \param        nt                      number of elements (input)
* \param   	vC			Third vertex.
*/
static int *GroupDivid(int LoopIndex[], int *NumOfGroup, int *Continue, int nt,  WlzErrorNum *dstErr)
{
  int  i, j, k, jg, IbreakPoint, FirstIofG, loopI, IthGroup;
  int *groupI, *NumOfElInGroup;
  int test = 0;

  WlzErrorNum errNum = WLZ_ERR_NONE;

  if((groupI = (int *)AlcMalloc(sizeof(int) * nt)) == NULL)
  {
     errNum = WLZ_ERR_MEM_ALLOC;
    *dstErr = errNum;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *(groupI)  = 1;
            j  = 1;
	    IbreakPoint = 1;
	    IthGroup    = 1;
    loopI = LoopIndex[0];


    for(i=1; i<nt; i++)
    {
         *(groupI+i) = -10;
    }
    do
    {  
       for(i=1; i<nt; i++)
       {
           if( LoopIndex[i] == loopI )
	   {
              *(groupI+i) = IthGroup;
	       j++;
	   }
	   else
	   {
              if(IbreakPoint == 1 && ( *(groupI + i) < 0) )
	      {
                 FirstIofG   = i;
	         IbreakPoint = 0;
	      }
	   }
       }
       if(j < nt)
       {
         loopI = LoopIndex[FirstIofG];
	 IthGroup++;
	 IbreakPoint = 1;
       }	 
    } while( j < nt);
    *NumOfGroup =  IthGroup;

     if((NumOfElInGroup = (int *)AlcMalloc(sizeof(int) * IthGroup )) == NULL)
     {
         errNum = WLZ_ERR_MEM_ALLOC;
        *dstErr = errNum;
     }
     
     /*  get the number of Element in each group */
     /*  and get whether it's continue ? */
     if(errNum == WLZ_ERR_NONE)
     {
        jg = 0;
	j  = 0;
        for(i=0; i< IthGroup; i++)
	{
          for(k=jg; k<nt; k++)
	  {
            /*   */
            if(groupI[k] == i+1)
	    {
	       j++;
	       if(j == nt)
	       {
                  NumOfElInGroup[i] = nt;
		  jg  = nt;
	       }
	       if( (i > 0) && (jg+j  == nt) )
	       {
	          NumOfElInGroup[i] = j;
		  jg += j;
	       }  
	    }   
	    else
	    {
	       NumOfElInGroup[i] = j;
	       jg  += j;
	       j   = 0;
	       break;
	    }   
	  }
	}
     }
   }  
   *Continue = 0;
   if(test) 
        printf("jg  %d  nt  %d\n",jg, nt);
   if( jg == nt  )
       *Continue = 1;
   if(*Continue == 0)
   {
       printf("non - correct:\n");
       for(i=0; i<nt; i++)
	  printf("GD  %d\n", groupI[i]);
       /* get the correct number */
       for(i=0; i< IthGroup; i++)
       {
	  NumOfElInGroup[i] = 0;
          for(k=0; k<nt; k++)
	  {
            if(groupI[k] == i+1)
               NumOfElInGroup[i] += 1;
	  }
       }
    }
    /* if(groupI) */
    {
      AlcFree(groupI);
      groupI = NULL;
    }  

    return NumOfElInGroup;
}



/*!
* \ingroup      WlzFeatures
* \return			array of int	.
* \brief	Computes the number of elements in each group and its
*               continutity. 
*
* \param	LoopIndex[]		input index of loops
* \param       *NumOfGroup              output
* \param       *Continue                output
* \param        nt                      number of elements (input)
* \param   	vC			Third vertex.
*/
static int *GroupByNumber(int LoopIndex[], int *NumOfGroup, int nt,  WlzErrorNum *dstErr)
{
  int  i, j, k, jg, IbreakPoint, FirstIofG, loopI, IthGroup;
  int *groupI;
  WlzErrorNum errNum = WLZ_ERR_NONE;

  if((groupI = (int *)AlcMalloc(sizeof(int) * nt)) == NULL)
  {
     errNum = WLZ_ERR_MEM_ALLOC;
    *dstErr = errNum;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *(groupI)  = 1;
            j  = 1;
	    IbreakPoint = 1;
	    IthGroup    = 1;
    loopI = LoopIndex[0];


    for(i=1; i<nt; i++)
    {
         *(groupI+i) = -10;
    }
    do
    {  
       for(i=1; i<nt; i++)
       {
           if( LoopIndex[i] == loopI )
	   {
              *(groupI+i) = IthGroup;
	       j++;
	   }
	   else
	   {
              if(IbreakPoint == 1 && ( *(groupI + i) < 0) )
	      {
                 FirstIofG   = i;
	         IbreakPoint = 0;
	      }
	   }
       }
       if(j < nt)
       {
         loopI = LoopIndex[FirstIofG];
	 IthGroup++;
	 IbreakPoint = 1;
       }	 
   } while( j < nt);
   *NumOfGroup =  IthGroup;

  }
    return groupI;
}

static int BigestGroup(int  GroupNum[], int nt, WlzErrorNum *dstErr)
{
  WlzErrorNum errNum = WLZ_ERR_NONE;
  int i, j, ig;
  if( GroupNum == NULL)
  {
     errNum = WLZ_ERR_PARAM_NULL;
    *dstErr = errNum;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    j  =  GroupNum[0];
    ig = 1;
    for(i=1; i<nt; i++)
    {
       if(GroupNum[i] > j)
       {
            j  = GroupNum[i];
	    ig = i+1;
       }
    }
  }
  return ig;
}

static void InAndOutSurfacePoints( WlzGMModel  *gM2, 
                                   int          nThLoop, 
				   int          numberOfPixelsZ,
                                   WlzDVertex2  StandSampleP[NumberToTrack],
                                   double       distance,
				   WlzDVertex3 *InPoints,
				   int         *nOfVWithinTheS,
				   WlzDVertex3 *OutPoints,
				   int         *nOfVOutTheS,
				   int         *suc,
				   WlzErrorNum *dstErr )
{
    int i, j, k, m,  idx, ntV, loopDirection, di;
    int inN, outN;
    WlzVertexP         tLoopVData;

    WlzErrorNum errNum = WLZ_ERR_NONE;
    WlzDVertex2         segV[3],  fPoint, bPoint, *vNorm=NULL, *inPoints=NULL, *outPoints=NULL;
    inN  = *nOfVWithinTheS;
    outN = *nOfVOutTheS;
    *suc  = 1;
 
	   /*  calculate the in and out surface points  */
        
	   /*  get the number of vertices in this loop */
           ntV  =  NumberOfVerticesInTheLoopOfGM(gM2, nThLoop, &errNum); 

           if( errNum == WLZ_ERR_NONE   )
           {
              /*  tLoopVData = WlzVerticesOneLoopFromGM(gM2,  ntS, ntL, &errNum); */
              tLoopVData = WlzVerticesThisLoopOfGM(  gM2, nThLoop, ntV, &errNum  );
	      if(errNum != WLZ_ERR_NONE)
	      {
                 *suc =0;
		 return;
	      }
           }

           if( errNum == WLZ_ERR_NONE   )
           {

	    loopDirection = 1;
	    /*  determin the direciton and goto the position: */
            /*  get the position of the index */
            for(i=0; i<ntV; i++)
            {
              if(  ( (tLoopVData.d2 + i )->vtX  == StandSampleP[0].vtX ) &&
	           ( (tLoopVData.d2 + i )->vtY  == StandSampleP[0].vtY )
	        )
		{
		  j = i;
	          break;
		}  
            }
            for(i=0; i<ntV; i++)
            {
              if(  ( (tLoopVData.d2 + i )->vtX  == StandSampleP[5].vtX ) &&
	           ( (tLoopVData.d2 + i )->vtY  == StandSampleP[5].vtY )
	        )
 	       {
	         k = i;
	         break;
	       } 
            }
	    
            for(i=0; i<ntV; i++)
            {
              if(  ( (tLoopVData.d2 + i )->vtX  == StandSampleP[NumberToTrack-1].vtX ) &&
	           ( (tLoopVData.d2 + i )->vtY  == StandSampleP[NumberToTrack-1].vtY )
	        )
 	       {
	         m = i;
	         break;
	       } 
            }

              if(k < j)
		   loopDirection = 0;

	     /*  calculate in and out points of the surface */
	    if(ntV > 60)
	    {
	      if(   (  (vNorm     = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) * NumberToTrack)) == NULL) ||
	            (  (inPoints  = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) * NumberToTrack)) == NULL) ||
		    (  (outPoints = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) * NumberToTrack)) == NULL)
		) 
              {
               	errNum = WLZ_ERR_MEM_ALLOC;
              }
	      
              if( errNum == WLZ_ERR_NONE   )
              {

                idx=0;
		if( loopDirection == 1)
		  i = k;
		else
		  i = k;
	        do
	        {
		  /*
	          segV[0].vtX = (tLoopVData.d2  + i)->vtX; 
	          segV[0].vtY = (tLoopVData.d2  + i)->vtY;
	          segV[1].vtX = (tLoopVData.d2  + i+1)->vtX; 
	          segV[1].vtY = (tLoopVData.d2  + i+1)->vtY;	      
	          segV[2].vtX = (tLoopVData.d2  + i+2)->vtX; 
	          segV[2].vtY = (tLoopVData.d2  + i+2)->vtY;
		  */
	          segV[0].vtX = StandSampleP[0].vtX; 
	          segV[0].vtY = StandSampleP[0].vtY;
	          segV[1].vtX = StandSampleP[5].vtX; 
	          segV[1].vtY = StandSampleP[5].vtY;	      
	          segV[2].vtX = StandSampleP[NumberToTrack-1].vtX; 
	          segV[2].vtY = StandSampleP[NumberToTrack-1].vtY;

		  
		  /*  printf("%lg   %lg\n", segV[1].vtX, segV[1].vtY ); */
	         *(vNorm + idx) = WlzVerticesNormTriple2(segV[0], segV[1], segV[2]);
		  /*  get the in and out points using the normal and points */
		  fPoint.vtX =  segV[1].vtX + distance * (vNorm+idx)->vtX;
		  fPoint.vtY =  segV[1].vtY + distance * (vNorm+idx)->vtY;
		  bPoint.vtX =  segV[1].vtX - distance * (vNorm+idx)->vtX;
		  bPoint.vtY =  segV[1].vtY - distance * (vNorm+idx)->vtY;
		  /*  check its out or in */
		  di = WlzGeomTriangleSnArea2(segV[0], segV[1], fPoint); 
		  if( ( ( di >  0. ) && ( loopDirection == 1) ) ||
		      ( ( di <  0. ) && ( loopDirection == 0) )
		  )
		  {
		     /*  fPoint is on the left of line segV[0]->segV[1] */
		     (inPoints +idx)->vtX  = fPoint.vtX; 
		     (inPoints +idx)->vtY  = fPoint.vtY;
		     
		     (outPoints +idx)->vtX = bPoint.vtX;
		     (outPoints +idx)->vtY = bPoint.vtY;
		      i += 4;
		      idx++;	     
  
		  }
		  else if( ( ( di <  0. ) && ( loopDirection == 1) ) ||
		           ( ( di >  0. ) && ( loopDirection == 0) )
		  )
		  {
		     (outPoints +idx)->vtX =  fPoint.vtX; 
		     (outPoints +idx)->vtY =  fPoint.vtY;
		     
		     (inPoints +idx)->vtX  =  bPoint.vtX;
		     (inPoints +idx)->vtY  =  bPoint.vtY;
                     i += 4;
		     idx++;	
		  }
		  else
		  {

                     printf("On the surface!!!!!!!!!!!!! not counted it\n");
		     break;

		  }
		

                /* } while(idx < NumberToTrack/4 ); */
                } while(idx < 1 );
			
		/* store the in points  */
                if(errNum == WLZ_ERR_NONE)
		{
		   for(i=0; i< idx; i++)
		   {

		       (InPoints + i + inN )->vtX =   (inPoints + i)->vtX;
		       (InPoints + i + inN )->vtY =   (inPoints + i)->vtY;
		       (InPoints + i + inN )->vtZ =   (double)numberOfPixelsZ;
			/*	       
                       fprintf(infp, "%6.2f   %6.2f   %6.2f   %6.2f\n", (inPoints + i)->vtX, (inPoints + i)->vtY, 
		                                                 (double)numberOfPixelsZ, distance);
			*/					 
		   }
		   inN += idx;
		}

		/* store the out points */
                if(errNum == WLZ_ERR_NONE)
		{
		   for(i=0; i< idx; i++)
		   {
		       (OutPoints + i + outN )->vtX =   (outPoints + i)->vtX;
		       (OutPoints + i + outN )->vtY =   (outPoints + i)->vtY;
		       (OutPoints + i + outN )->vtZ =   (double)numberOfPixelsZ;
		       /*
                      fprintf(outfp, "%6.2f   %6.2f   %6.2f   %6.2f\n", (outPoints + i)->vtX, (outPoints + i)->vtY, 
		                                                 (double) numberOfPixelsZ, distance);
								 */
		   }
		   outN += idx;
		}


	      }

	    }
	   }

      /* if(tLoopVData.d2) */
              AlcFree(tLoopVData.d2);
      /* if(vNorm) */
      {
              AlcFree(vNorm);
	      vNorm     = NULL;

      }	      
      /* if(inPoints) */
      {
              AlcFree(inPoints);
              inPoints  = NULL;

      }
      /* if(outPoints) */
      {
             AlcFree(outPoints);
      	     outPoints = NULL;

      }
      *nOfVOutTheS = outN;
      *nOfVWithinTheS  = inN;

}

/*!
* \ingroup      WlzFeatures
* \return			array of int	.
* \brief	Computes the standard sample points in a specified loop
*               continutity. 
*
* \param	sLoopIdx		    input index of loops
* \param       *gM                          input geometric model 
* \param        StandSampleP[NumberToTrack] output sampled points
* \param       *dstErr			    Woolz err number pointer.
*/
static void GetSamplePointsFromOneLoopOfGM(     WlzGMModel  *gM, 
                                   int          LoopIdx,
                                   WlzDVertex2  StandSampleP[NumberToTrack],
                                   WlzErrorNum *dstErr )
{
    int   i, j, k, sp, nS, nV, nL, test;
    int  *sShellIdx=NULL, *sLoopIdx=NULL;
    WlzErrorNum errNum = WLZ_ERR_NONE;
    WlzVertexP          sLoopVData;

    /*  first  */
    if(LoopIdx < 0)
    {
        /*  get the number of Shells in this model */
        nS        = NumberOfShellsAboutGM(gM);
        /*  get the Index of for each Shell; */
        sShellIdx = WlzShellsIndexAboutGM(gM, nS, &errNum);


        if( errNum == WLZ_ERR_NONE   )
        {
           test = 0;
           if(test)
           {
              WlzOutputForTestGM( gM, &errNum);
              exit(0);
           }

           /*  use the first Shell and first loop as a test */
           i        = 0;


           /*  get number loops in this shell */
           nL       =  NumberOfLoopsInThe_nThShellOfGM(gM,  sShellIdx[i] );

           /*  get the index for the loops */
           sLoopIdx =  WlzLoopIndexAboutGM(gM, sShellIdx[i], nL, &errNum);
           

	  if( errNum == WLZ_ERR_NONE )
	  {
             /*  get howManyVertices in the loop: */
             j        = 0;
	     LoopIdx  = sLoopIdx[j];
             nV       = NumberOfVerticesInTheLoopOfGM(gM, LoopIdx, &errNum); 
             if( errNum == WLZ_ERR_NONE )
	     {
	        if(nV < 60)
	        {
                  printf("section too short, get anthor one\n");
		  j++;
		  LoopIdx =  sLoopIdx[j];
	        }
		nV       = NumberOfVerticesInTheLoopOfGM(gM, LoopIdx, &errNum); 
	        if(nV < 60)
	        {
                  printf("section too short, stop\n");
		  exit(0);
	        }
                printf("num Vertices = %d\n",nV);
	     }
          }

        }
	/* if(sShellIdx) */
	{
	   AlcFree(sShellIdx);
	   sShellIdx = NULL;
	}   
	/* if(sLoopIdx) */
	{
	   AlcFree(sLoopIdx);
	   sLoopIdx = NULL;
	}   

    }

    if( errNum == WLZ_ERR_NONE   )
    {
        
        nV       = NumberOfVerticesInTheLoopOfGM(gM, LoopIdx, &errNum); 
        if( errNum == WLZ_ERR_NONE )
	{
            /* printf("num Vertices = %d\n",nV); */
             sLoopVData = WlzVerticesThisLoopOfGM( gM, LoopIdx, nV, &errNum  );
               
	     if( errNum == WLZ_ERR_NONE )
	     {
	       /*  output for forming surface */
               test = 0;
               if(test)
               {
                    WlzOutputForTestAsectionData( sLoopVData, nV, &errNum );
                    exit(0);
               } 
	     }
         }
      }
      /* get the stand sample points first */
      if(nV < 62)
      {
           printf("This section is too short, please change a section and try again!\n");
	   exit(0);
      }
      /* j = nV/NumberToTrack; */
      j = 3;

      for(i=0; i<NumberToTrack; i++)
      {
	   /*  using sample position, sp should not be too big as it will push  */
	   /*  the index out of range!!!: */
	   sp = 1;
	   k = i* j;
           StandSampleP[i].vtX =  (sLoopVData.d2 + sp + k )->vtX;
           StandSampleP[i].vtY =  (sLoopVData.d2 + sp + k )->vtY;
      }
      AlcFree(sLoopVData.d2);
      
}


/*!
* \ingroup      WlzFeatures
* \return			array of int	.
* \brief	get the tracked sample points in a tracked loop
*               continutity. 
*
* \param	tLIdx		              output index of the tracked loop
* \param       *gM2                           input geometric model 
* \param        StandSampleP[NumberToTrack]   input  sampled points to be tracked
* \param        TrackedSecondP[NumberToTrack] output sampled points
* \param        nOfTracked                    output num Of Tracked Points belong to one loop
* \param       *suc                           output output pointer to indicate whether this
*                                             track successful. 
* \param       *dstErr			      Woolz err number pointer.
*/
static void GetTrackedSamplePointsFromOneLoopOfGM( WlzGMModel  *gM2,
                                    AlcKDTTree         *tTree,
                                   int          *LoopIdx, 
				   int           numberOfPixelsZ,
				   double        minDis,
                                   WlzDVertex2   StandSampleP[NumberToTrack],
                                   WlzDVertex2   TrackedSamplePBlgToOneLp[NumberToTrack],
				   int          *nOfTracked,
				   int          *suc,
                                   WlzErrorNum *dstErr )
{
        
        int                 i, j, k, m, idx, ntL;
	int                 test = 0;
	int                 nShell, nLoop;
	int                 Continue;
        double              datD[3];
        AlcKDTNode         *node;
	AlcErrno            errNo = ALC_ER_NONE;
        int                 LoopIndex[NumberToTrack];
        int                 ShellIndex[NumberToTrack];
        int                 indexStr[NumberToTrack];
	int                 indexStrC[NumberToTrack];
        int                *GroupNum, *FixNumForEachGroup;
        int                 NumOfGroup;
        int                 sp, nbgN, trackedBigestLoopIndex;
	int                 nV;
        int                 loopDirection;
        WlzDVertex2         TrackedSampleFirstP, TrackedSampleSecondP;
	WlzDVertex2         TrackedSampleP[NumberToTrack];
        WlzVertexP          sLoopVData;

        WlzErrorNum         errNum = WLZ_ERR_NONE;
	*suc = 1;
        /*  --- find the nearst node (one or nil )---- */
	for(i=0; i<NumberToTrack; i++)
	{
          datD[0] = StandSampleP[i].vtX;
          datD[1] = StandSampleP[i].vtY;
	  datD[2] = 0.;
          /*  printf("Point %lg %lg\n",datD[0], datD[1]); */
          node    = AlcKDTGetNN(tTree, datD, minDis, NULL, &errNo);
     	  if(  node && (errNum == ALC_ER_NONE) )
  	  {
             TrackedSampleP[i].vtX = *(node->key.kD);
             TrackedSampleP[i].vtY = *(node->key.kD + 1);
	     indexStr[i]           = node->idx;
  	  }
	  else
	  {
             printf(" cannot track this points!!! you have \n \
	              to give a biger mindistance \n \
	              for nearest point to track! or there \
		      are something wrong!\n");
	     *suc = 0;
	     return;
	  }
          /*   get the loop corresponding  to the tarcked point */
          /*   it may failed as the crossed point belongs to */
	  /*   to loops !! */
	  /* printf("index = %d\n",indexStr[i] ); */
	  if(indexStr[i] > 5000000)
	  {
           printf("Something wrong\n");
	     *suc = 0;
	     return;

	  }
          WlzGetTheShellAndLoopNumberFromTheIndexOfTheVertex( gM2,
                                       indexStr[i], 
	             	              &nShell, 
			              &nLoop,
			              &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
            LoopIndex[i]  = nLoop;
	    ShellIndex[i] = nShell;
	  }

	}

	/*  check, how many loops in the target group: */
        /*  check, whether it continue ?: */
	Continue   = 0;
	/*  get the accurate number of points for each group */
	/*  and answer the continunities of the data */
	GroupNum   = GroupDivid(LoopIndex, &NumOfGroup, &Continue, NumberToTrack, &errNum);
	
        if(errNum == WLZ_ERR_NONE)
        {
	  /*  get the fixed number for each group */
          FixNumForEachGroup = GroupByNumber(LoopIndex, &NumOfGroup, NumberToTrack, &errNum);
	}
	
	if(errNum == WLZ_ERR_NONE)
        {
	  if(test)
	  {
	     printf("number of groups = %d\n",NumOfGroup );
	     printf("Cont  %d\n", Continue );
	     {
	        for(i=0; i<NumOfGroup; i++)
	        {
                    printf("Num in the %d group = %d\n",i, GroupNum[i]);
	            for(j=0; j<NumberToTrack; j++)
	            {
	                if(FixNumForEachGroup[j] == i+1)
	                {
	                     printf("idx = %d\n", j);

	                }

	            }
	        }
	     }
	  }   


          /*  we have to choose the right one */
          nbgN = BigestGroup(GroupNum, NumOfGroup, &errNum);
	  if(test)
	     printf("group  %d is the bigest group\n",nbgN);
	     
	  if(GroupNum[nbgN-1] < (NumberToTrack* 65)/100 )
	  {
              printf("track failed as no group is a dominant one\n");
	      /* if(GroupNum) */
	         AlcFree(GroupNum);
              /* if(FixNumForEachGroup)		  */
	         AlcFree(FixNumForEachGroup);
              *suc = 0;
	  }

	}

	/*  */
        if(errNum == WLZ_ERR_NONE)
        {
	     idx = 0;
  	     for(i=0; i<NumberToTrack; i++)
	     {
	         if(FixNumForEachGroup[i] == nbgN)
		 {
		   TrackedSamplePBlgToOneLp[idx].vtX =  TrackedSampleP[i].vtX;
		   TrackedSamplePBlgToOneLp[idx].vtY =  TrackedSampleP[i].vtY;
		   indexStr[idx] = indexStr[i];
                   /* fprintf(testFile, "%6.2f  %6.2f   %6.2f\n",TrackedSampleP[i].vtX,  */
	           /*                         TrackedSampleP[i].vtY, (double) numberOfPixelsZ); */
		  if(idx == 0)
		  {
		   TrackedSampleFirstP.vtX  = TrackedSampleP[i].vtX;
		   TrackedSampleFirstP.vtY  = TrackedSampleP[i].vtY;
		  }
		  else if(idx == 1)
		  {
		   TrackedSampleSecondP.vtX = TrackedSampleP[i].vtX;
		   TrackedSampleSecondP.vtY = TrackedSampleP[i].vtY;
		  } 
		  idx++;
		}			   
	     }

	}
	*nOfTracked = idx;

        if(errNum == WLZ_ERR_NONE)
        {
	    /*  get the loop index we have just tracked */
	    for(i=0; i<NumberToTrack; i++)
	    {
              if(FixNumForEachGroup[i] == nbgN )
	      {
                  /*  printf("Loop index tracked = %d\n", LoopIndex[i]); */
	   	  ntL = LoopIndex[i];
		  break;

	      }

	    }
       }
       *LoopIdx = ntL;
       if(ntL  < 0)
           *suc = 0;
    /* if(GroupNum)    */
       AlcFree(GroupNum);
    /* if(FixNumForEachGroup)    */
       AlcFree(FixNumForEachGroup);
}


/*!
* \ingroup      WlzFeatures
* \return	AlcKDTree.
* \brief	get the tracked sample points in a tracked loop
*               continutity. 
*
* \param	tLIdx		              output index of the tracked loop
* \param       *gM2                           input geometric model 
* \param        StandSampleP[NumberToTrack]   input  sampled points to be tracked
* \param        TrackedSecondP[NumberToTrack] output sampled points
* \param       *dstErr			    Woolz err number pointer.
*/
static AlcKDTTree *GetTreeFromGM( WlzObject  *tObj,
                                  WlzErrorNum *dstErr )
{
    AlcKDTTree         *tTree;
    int                 vCnt[2], idN;
    int                *sNN, test = 0;
    WlzVertexP         nData[2], vData[2];
    WlzVertexType       vType;
    WlzErrorNum         errNum = WLZ_ERR_NONE;
        /*  get the vertices from obj */

        for(idN = 0; idN < 2; ++idN)
        {
           if(errNum == WLZ_ERR_NONE)
           {
               if(idN == 0)
                 vData[idN] = WlzVerticesFromObj(tObj, nData + idN, vCnt + idN, &vType,
			       &errNum);
               else
                 vData[idN] = WlzVerticesFromObj(tObj, nData + idN, vCnt + idN, &vType,
	               	       &errNum);
           }
        }

		     

		     
        if(errNum == WLZ_ERR_NONE)
	{
            if(test)
            {
                printf("Start\n");
                idN =1;
                {

                    WlzOutputForTestAsectionData(vData[idN], vCnt[idN], &errNum);
                    if(errNum == WLZ_ERR_NONE)
                    {
                        WlzOutputForTestAsectionData(nData[idN], vCnt[idN], &errNum);
                    }	
                }
                printf("End\n");
                exit(0);
            }
	}	 
 
    if( errNum == WLZ_ERR_NONE   )
    {
	 /*  build up a KD Tree for the target vertices of the Woolz contour  */
         if(  (sNN = (int *)AlcMalloc(sizeof(int) * vCnt[1])) == NULL )
         {
              errNum = WLZ_ERR_MEM_ALLOC;
         }

         if(errNum == WLZ_ERR_NONE)
         {
              tTree = WlzVerticesBuildTree(vType, vCnt[1], vData[1], sNN, &errNum);
         }
    }
    /* if(sNN) */
       AlcFree(sNN);
    /* if(nData[0].d2)   */
       AlcFree(nData[0].d2);
    /* if(vData[0].d2)    */
       AlcFree(vData[0].d2);
    /* if(nData[1].d2)          */
       AlcFree(nData[1].d2);
    /* if(vData[1].d2)   */
       AlcFree(vData[1].d2);
      

    return tTree;

}








/*!
* \ingroup      WlzFeatures
* \return			array of int	.
* \brief	Computes the standard sample points in a specified loop
*               continutity. 
*
* \param	sLoopIdx		    input index of loops
* \param        i_sec                       input integer to indicate section 
* \param       *gM                          input geometric model 
* \param        TrackedSamplePBlgToOneLp[NumberToTrack] input tracked points
* \param        StandSampleP[NumberToTrack] output sampled points
* \param       *dstErr			    Woolz err number pointer.
*/
static void GetTrackedPlusSomeSamplePointsFromOneLoopOfGM(
                           int          j_shell,
			   int          i_sec,
			   int          i_file,
			   int          sectionLength_L,
                           int         *numOfPInS,
                           WlzIVertex2 *sectionData,
			   WlzGMModel  *gM, 
                           int          LoopIdx, 
                           WlzDVertex2  StandSampleP[NumberToTrack],
		   	   WlzDVertex2	TrackedSamplePBlgToOneLp[NumberToTrack],
			   int		nOfTracked,
			   int         *suc,
                           WlzErrorNum *dstErr )
{
    int   i, j, k, m, sp, nS, nV, nL, test, istar, iend, icount;
    int  *sShellIdx=NULL, *sLoopIdx=NULL;
    int   deltaV;
    WlzErrorNum errNum = WLZ_ERR_NONE;
    WlzVertexP          sLoopVData;
    AlcErrno            alcErr = ALC_ER_NONE;
   *suc = 1;

     /*  change to sample the bigest tracked line segment as next track starting points */
     deltaV = sectionLength_L/NumberToTrack;
    for(i=0; i<NumberToTrack; i++)
    {
          StandSampleP[i].vtX =  TrackedSamplePBlgToOneLp[i].vtX;
          StandSampleP[i].vtY =  TrackedSamplePBlgToOneLp[i].vtY;
    }

    if(nOfTracked >= NumberToTrack)  
    {
       nV       = NumberOfVerticesInTheLoopOfGM(gM, LoopIdx, &errNum); 
       if( errNum == WLZ_ERR_NONE )
       {
             sLoopVData = WlzVerticesThisLoopOfGM( gM, LoopIdx, nV, &errNum  );
       }
      
       if( errNum == WLZ_ERR_NONE )
       {

         for(i=0; i<nV; i++)
         {
           if(  ( (sLoopVData.d2 + i )->vtX  == TrackedSamplePBlgToOneLp[0].vtX ) &&
	        ( (sLoopVData.d2 + i )->vtY  == TrackedSamplePBlgToOneLp[0].vtY )
	   )
	   break;
         }
     
         for(j=0; j<nV; j++)
         {
           if(  ( (sLoopVData.d2 + j )->vtX  == TrackedSamplePBlgToOneLp[1].vtX ) &&
	        ( (sLoopVData.d2 + j )->vtY  == TrackedSamplePBlgToOneLp[1].vtY )
	   )
	   break;
         }

         for(k=0; k<nV; k++)
         {
           if(  ( (sLoopVData.d2 + k )->vtX  == TrackedSamplePBlgToOneLp[nOfTracked-1].vtX ) &&
	        ( (sLoopVData.d2 + k )->vtY  == TrackedSamplePBlgToOneLp[nOfTracked-1].vtY )
	   )
	   break;
         }
     
	 istar = i;
	 iend  = k-1;
	 
         if(k < i)
	 {
	   istar = k;
	   iend  = i-1;
	 }
	 
	 if(iend - istar > sectionLength_L)
	     iend = istar + sectionLength_L;

         icount = 0;
         for(j=istar; j<iend; j++)
         {
	     ( sectionData + icount )->vtX = 	 (int ) (sLoopVData.d2 + j )->vtX;
	     ( sectionData + icount )->vtY =     (int ) (sLoopVData.d2 + j )->vtY;
             icount++;
         }
	 
         *( numOfPInS + i_file ) = icount;
	 if(icount < 3)
	    *suc = 0;
      }
      /* if(sLoopVData.d2) */
         AlcFree(sLoopVData.d2);
      return;
    }
    else
    {
       nV       = NumberOfVerticesInTheLoopOfGM(gM, LoopIdx, &errNum); 
       if( errNum == WLZ_ERR_NONE )
       {
             sLoopVData = WlzVerticesThisLoopOfGM( gM, LoopIdx, nV, &errNum  );
       }
       /* get the stand sample points first */
       if(nV < 62)
       {
           printf("This section is too short, please change a section and try again!\n");
	  /*  if(sLoopVData.d2) */
	      AlcFree(sLoopVData.d2);
	  *suc = 0;
	   return; 
       }

       /*  get the position of the index */
       for(i=0; i<nV; i++)
       {
           if(  ( (sLoopVData.d2 + i )->vtX  == TrackedSamplePBlgToOneLp[nOfTracked-1].vtX ) &&
	        ( (sLoopVData.d2 + i )->vtY  == TrackedSamplePBlgToOneLp[nOfTracked-1].vtY )
	   )
	   break;
       }
       
       j = 2;

       for(m = nOfTracked; m<NumberToTrack; m++)
       {
	   /*  using sample position, sp should not be too big as it will push  */
	   /*  the index out of range!!!: */
	   sp = 1;
	   k  = i +  ( m - nOfTracked + 1) * j;
	   if(k < nV)
	   {
             StandSampleP[m].vtX =  (sLoopVData.d2 + k )->vtX;
             StandSampleP[m].vtY =  (sLoopVData.d2 + k )->vtY;
	   }
	   else
	   {
	    /* if(sLoopVData.d2) */
	        AlcFree(sLoopVData.d2);
                printf("not that much points to add in the loop");
		*suc = 0;
		return;
	   }
       }
    }  

    /*  store the section data for making a contour */
    if( errNum == WLZ_ERR_NONE )
    {
        for(i=0; i<nV; i++)
        {
           if(  ( (sLoopVData.d2 + i )->vtX  == TrackedSamplePBlgToOneLp[0].vtX ) &&
	        ( (sLoopVData.d2 + i )->vtY  == TrackedSamplePBlgToOneLp[0].vtY )
	   )
	   break;
        }

        /*  */
	 istar = i;
	 iend  = k-1;
	 
         if(k < i)
	 {
	   istar = k;
	   iend  = i-1;
	 }
 
         if(iend - istar > sectionLength_L)
	        iend = istar + sectionLength_L;


         icount = 0;
         for(j=istar; j<iend; j++)
         {  
	     ( sectionData + icount )->vtX = 	 (int ) (sLoopVData.d2 + j )->vtX;
	     ( sectionData + icount )->vtY =     (int ) (sLoopVData.d2 + j )->vtY;
             icount++;
         }
	 
         *( numOfPInS + i_file ) = icount;
	 if(icount < 3)
	    *suc = 0;
         /* if(sLoopVData.d2) */
             AlcFree(sLoopVData.d2);
      return;

    }

    /* if(sLoopVData.d2) */
      AlcFree(sLoopVData.d2);

}



/*!
* \ingroup      WlzFeatures
* \return			array of int	.
* \brief	Computes the standard sample points in a specified loop
*               continutity. 
*
* \param	sLoopIdx		    input index of loops
* \param        i_sec                       input integer to indicate section 
* \param       *gM                          input geometric model 
* \param        TrackedSamplePBlgToOneLp[NumberToTrack] input tracked points
* \param        StandSampleP[NumberToTrack] output sampled points
* \param       *dstErr			    Woolz err number pointer.
*/
static void GetSamplePointsForNextTrackFromTheTrackedLoopOfGM(
                           int          i_file,
			   int          sectionLength_L,
                           int         *numOfPInS,
                           WlzIVertex2 *sectionData,
			   WlzGMModel  *gM, 
                           int          LoopIdx, 
                           WlzDVertex2  StandSampleP[NumberToTrack],
		   	   WlzDVertex2	TrackedSamplePBlgToOneLp[NumberToTrack],
			   int		nOfTracked,
			   int         *suc,
                           WlzErrorNum *dstErr )
{
    int   i, j, k, m, sp, nS, nV, nL, test, istar, iend, icount;
    int  *sShellIdx=NULL, *sLoopIdx=NULL;
    int   deltaV;
    int   oldversion = 1;
    double  x0, y0, x1, y1, dis, stand;
    WlzErrorNum errNum = WLZ_ERR_NONE;
    WlzVertexP          sLoopVData;
    AlcErrno            alcErr = ALC_ER_NONE;
   *suc = 1;

     /*  change to sample the bigest tracked line segment as next track starting points */
     deltaV = sectionLength_L/NumberToTrack;

     nV       = NumberOfVerticesInTheLoopOfGM(gM, LoopIdx, &errNum); 
     if( errNum == WLZ_ERR_NONE )
     {
         sLoopVData = WlzVerticesThisLoopOfGM( gM, LoopIdx, nV, &errNum  );
     }
      
     if( errNum == WLZ_ERR_NONE )
     {
         istar = 1000000;
	 iend  = -1;
	 for(k=0; k<nOfTracked; k++)
	 {

           for(i=0; i<nV; i++)
           {
              if(  ( (sLoopVData.d2 + i )->vtX  == TrackedSamplePBlgToOneLp[k].vtX ) &&
	           ( (sLoopVData.d2 + i )->vtY  == TrackedSamplePBlgToOneLp[k].vtY )
	      )
	      {
                 if(i> iend)
		      iend  = i;
		 if(i< istar)
		      istar = i;
	      }
	   }   
         }
     
	 /*  test the validity: */
	 if( ( istar + (NumberToTrack-1) * (deltaV) ) > nV )
	 {
            *suc = 0;
	 }
         
	 /*  make sure the number tracked is a real marjority here !!! */
	 if( *suc )
	 {
            if(  (iend - istar) < ( ( sectionLength_L * 70) /100 )  )
	    {
               *suc = 0;
	    }
	 }

        /*  old version */
	 if( *suc && oldversion )
	 {
                /*  get the sample points */
	         for(j=0; j<NumberToTrack; j++)
	         {
                    m = istar+j*(deltaV);
                    StandSampleP[j].vtX =  (sLoopVData.d2 + m )->vtX;
                    StandSampleP[j].vtY =  (sLoopVData.d2 + m )->vtY;
                 }

		 
	       iend = istar + sectionLength_L;
               icount = 0;
               for(j=istar; j<iend; j++)
               {  
	          ( sectionData + icount )->vtX =  (int ) (sLoopVData.d2 + j )->vtX;
	          ( sectionData + icount )->vtY =  (int ) (sLoopVData.d2 + j )->vtY;
                  icount++;
               }
               *( numOfPInS + i_file ) = icount;

	 }


	  /*  new version */
         /*  determine the filling method */
         if( ( *suc && !oldversion ) )
	 {
	      m = istar + (NumberToTrack-1) * (deltaV);
	      x0 = (sLoopVData.d2 + m )->vtX;
	      y0 = (sLoopVData.d2 + m )->vtY;

              if(k < i )
	      {
                 x1 = TrackedSamplePBlgToOneLp[0].vtX; 
                 y1 = TrackedSamplePBlgToOneLp[0].vtY; 
	      }
	      else
	      {
                 x1 = TrackedSamplePBlgToOneLp[nOfTracked-1].vtX; 
                 y1 = TrackedSamplePBlgToOneLp[nOfTracked-1].vtY; 
	      }
              dis   = (x1-x0) * ( x1-x0 ) + (y1-y0) * (y1-y0);
	      stand = sectionLength_L * sectionLength_L;
              if(dis < stand )
	      {
                 /*  get the sample points */
	         for(j=0; j<NumberToTrack; j++)
	         {
                    m = istar+j*(deltaV);
                    StandSampleP[j].vtX =  (sLoopVData.d2 + m )->vtX;
                    StandSampleP[j].vtY =  (sLoopVData.d2 + m )->vtY;
                 }
	      }
	      else
	      {
                  for(j=0; j<NumberToTrack; j++)

            	  {
                    m = istar - j*(deltaV);
		    if(m < 0)
		    {
		        *suc = 0;
			break;
		    }	
                      StandSampleP[j].vtX =  (sLoopVData.d2 + m )->vtX;
                      StandSampleP[j].vtY =  (sLoopVData.d2 + m )->vtY;
                  }
	      }
	 }

	 if(*suc && !oldversion )
	 {
	     icount = 0;
	     if(dis < stand )
	     {
	       iend = istar + sectionLength_L;

               for(j=istar; j<iend; j++)
               {  
	          ( sectionData + icount )->vtX =  (int ) (sLoopVData.d2 + j )->vtX;
	          ( sectionData + icount )->vtY =  (int ) (sLoopVData.d2 + j )->vtY;
                  icount++;
               }
               *( numOfPInS + i_file ) = icount;
	     }  
	     else
	     {
               iend = istar - sectionLength_L;
               for(j=istar; j>iend; j--)
               {  
	          ( sectionData + icount )->vtX =  (int ) (sLoopVData.d2 + j )->vtX;
	          ( sectionData + icount )->vtY =  (int ) (sLoopVData.d2 + j )->vtY;
                  icount++;
               }
               *( numOfPInS + i_file ) = icount;
	     }

	 }

      }
      /* if(sLoopVData.d2) */
           AlcFree(sLoopVData.d2);
      return;
}














/*!
*  \return      int 0 for non-cycle (shell);
*  \ingroup     WlzFeatures
*  \brief       Determine whether the line(shell) relating with this loop is a cycle.
*               a better way to do it seems to know how many loops with a
*               shell in a given model.
*               if a shell has two loops, it's a cycle. one loops is a
*               section. This is not a correct way!!!!!
*  \param       gM                  given G model (input)
*  \param       nThLoop             n-th loop of the n-th shell in the model (input )
*/
static int IsACycle(   WlzGMModel  *gM, 
                       int          LoopIdx, 
                       WlzErrorNum *dstErr )
{
   int cycle, nLoop;
   WlzErrorNum errNum = WLZ_ERR_NONE;
   int nV1, nV2;
   int *tLoopIdx, Sidx;
   AlcVector          *vec;
   WlzGMLoopT         *cLT, *fLT;
   WlzGMShell         *cS;
  
   nLoop = 0;
   cycle = 0;
   vec    =  gM->res.loopT.vec;

   /* get the loop corresponding to the LoopIdx of the model. */
    cLT = (WlzGMLoopT *) AlcVectorItemGet(vec, LoopIdx);

   /* get the Shell it belongs to */
   if(cLT->parent == NULL)
   {
       printf("no parent shell for this loop");
       exit(0);
   }

    cS =  (WlzGMShell *) cLT->parent;
    Sidx = cS->idx;

    if(Sidx < 0)
        return(cycle);

    cLT = fLT = cS->child;
    if(!cLT)
        return(cycle);
    /* For each loop topology element of the model. */
    do
    {
        cLT = cLT->next;
        ++nLoop;
    } while(  cLT != fLT );


    if(nLoop == 2)
    {
      nV1 = -3;
      nV2 = -1;
      if(errNum == WLZ_ERR_NONE )
      {
        tLoopIdx = WlzLoopIndexAboutGM(gM, Sidx, nLoop, &errNum);
      }
      if(errNum == WLZ_ERR_NONE )
      {
         nV1 = NumberOfVerticesInTheLoopOfGM(gM, tLoopIdx[0], &errNum);
      }
      if(errNum == WLZ_ERR_NONE )
      {
         nV2 = NumberOfVerticesInTheLoopOfGM(gM, tLoopIdx[1], &errNum);
      } 
      if(tLoopIdx)
         AlcFree(tLoopIdx);
      if(nV1 == nV2)
                cycle = 1;
    }
    return(cycle);
}

/*!
*  \return      int number of ends in the non-cycle line( shell with one loop);
*  \ingroup     WlzFeatures
*  \brief       get the number of ends for a non-cycle shell.
*  \param       gM                  given G model (input) 
*  \param       nThLoop             n-th loop of the n-th shell in the model (input)
*/
static int HowManyEndsInTheNonCycleLine(   WlzGMModel  *gM, 
                                           int          LoopIdx, 
                                           WlzErrorNum *dstErr )
{
   int numberOfEnds, numBAgain;
   WlzErrorNum errNum = WLZ_ERR_NONE;
   int nVertices;
   AlcVector          *vec;
   WlzGMLoopT         *cLT;
   WlzGMShell         *cS;
   WlzGMEdgeT         *cET, *fET;
  
   numberOfEnds = 0;
   numBAgain    = 0;
   vec    =  gM->res.loopT.vec;

   /* get the loop corresponding to the LoopIdx of the model. */
    cLT = (WlzGMLoopT *) AlcVectorItemGet(vec, LoopIdx);

   /* For each edge topology element of the model. */
    cET = fET = cLT->edgeT;

   /* For each loop topology element of the model. */
    do
    {
	if( ( cET->prev->vertexT->diskT->idx == cET->next->vertexT->diskT->idx ) )
	    numberOfEnds  +=1;
	cET = cET->next;
    } while (cET != fET);

    /*  counter again */
    cET = fET = cLT->edgeT->next;

   /* For each loop topology element of the model. */
    do
    {
	if( ( cET->prev->vertexT->diskT->idx == cET->next->vertexT->diskT->idx ) )
	    numBAgain  +=1;
	cET = cET->next;
    } while (cET != fET);

   if( numBAgain > numberOfEnds )
          return(numBAgain);

   return(numberOfEnds);
}


/*!
*  \return      index of the Ends points ;
*  \ingroup     WlzFeatures
*  \brief       get the number of ends for a non-cycle shell.
*  \param       gM                  given G model (input ) 
*  \param       nThLoop             n-th loop of the n-th shell in the model (input )
*  \param       numOfEnds           number of Ends points in n-th shell in the model (input )
*/
static int *GetTheEndsPointsIndexInTheNonCycleLine(   WlzGMModel  *gM, 
                                           int          LoopIdx,
					   int          numOfEnds,
                                           WlzErrorNum *dstErr )
{
   int          numberOfEnds, numBAgain;
   int         *EndPointsIndex = NULL;
   WlzErrorNum errNum = WLZ_ERR_NONE;
   AlcVector          *vec;
   WlzGMLoopT         *cLT;
   WlzGMShell         *cS;
   WlzGMEdgeT         *cET, *fET;

  if( ( EndPointsIndex = (int *)AlcMalloc(sizeof(int) * numOfEnds )) == NULL ) 
  {
     errNum = WLZ_ERR_MEM_ALLOC;
  }
  

  if( errNum == WLZ_ERR_NONE)
  {
     numberOfEnds = 0;
     numBAgain    = 0;
     vec          =  gM->res.loopT.vec;

     /* get the loop corresponding to the LoopIdx of the model. */
      cLT = (WlzGMLoopT *) AlcVectorItemGet(vec, LoopIdx);

     /* For each edge topology element of the model. */
      cET = fET = cLT->edgeT;

     /* For each loop topology element of the model. */
     do
     {
	if( ( cET->prev->vertexT->diskT->idx == cET->next->vertexT->diskT->idx ) )
	{
	    *(EndPointsIndex + numberOfEnds ) = cET->vertexT->diskT->idx;
	    numberOfEnds  +=1;
	}    
	cET = cET->next;
     } while (cET != fET);

     if(numberOfEnds != numOfEnds)
     {

       /*  counter again */
       cET = fET = cLT->edgeT->next;

       /* For each loop topology element of the model. */
       do
       {
	 if( ( cET->prev->vertexT->diskT->idx == cET->next->vertexT->diskT->idx ) )
	 {
	    *(EndPointsIndex + numBAgain ) = cET->vertexT->diskT->idx;
	     numBAgain  +=1;
	 }    
	  cET = cET->next;
       } while (cET != fET);
     }  
  }

    *dstErr = errNum;
  return (EndPointsIndex);
}

/*!
*  \return      index of the bifurcating points ;
*  \ingroup     WlzFeatures
*  \brief       get the index of bifurcating points for a non-cycle shell.
*               The number of bifurcations equals to number of ends points
*                    minus two.
*  \param       gM                  given G model (input ) 
*  \param       nThLoop             n-th loop of the n-th shell in the model (input )
*  \param       numOfEnds           number of Ends points in n-th shell in the model (input )
*                                    > 2 (otherwise no bifurcate point)
*/
static int *GetThebifurcatePointsIndexInTheNonCycleLine(   WlzGMModel  *gM, 
                                           int          LoopIdx,
					   int          numOfEnds,
                                           WlzErrorNum *dstErr )
{
   int          numberOfBifurP, i;
   int         *BifurcatePointsIndex = NULL;
   WlzErrorNum errNum = WLZ_ERR_NONE;
   AlcVector          *vec;
   WlzGMLoopT         *cLT;
   WlzGMShell         *cS;
   WlzGMEdgeT         *cET, *fET;
   if( numOfEnds < 3 )
   {
      printf("no bifurcation in this section of line!\n");
      exit(0);
   }

   if( ( BifurcatePointsIndex = (int *)AlcMalloc(sizeof(int) * (numOfEnds-2) )) == NULL ) 
   {
      errNum = WLZ_ERR_MEM_ALLOC;
   }
  

  if( errNum == WLZ_ERR_NONE)
  {
     numberOfBifurP = 0;
     vec          =  gM->res.loopT.vec;

     /* get the loop corresponding to the LoopIdx of the model. */
      cLT = (WlzGMLoopT *) AlcVectorItemGet(vec, LoopIdx);

     /* For each edge topology element of the model. */
      cET = fET = cLT->edgeT;

     /* For each loop topology element of the model. get one of the end points as
        a start point */
     do
     {
	if( ( cET->opp->next->next->vertexT->diskT->idx != cET->prev->vertexT->diskT->idx ) )
	{
            if( numberOfBifurP == 0)
	    {
	     *(BifurcatePointsIndex + numberOfBifurP ) = cET->vertexT->diskT->idx;
	       numberOfBifurP  +=1;
	    }
	    else
	    {
	      /*  check whether it already in the bank */
	      for(i=0; i<numberOfBifurP; i++)
	      {
                  if( *(BifurcatePointsIndex + i ) == cET->vertexT->diskT->idx )
		        break;
	      }
	      if(i == numberOfBifurP)
	      {
	         /*  not in the bank */
		 *(BifurcatePointsIndex + numberOfBifurP ) = cET->vertexT->diskT->idx;
	         numberOfBifurP  +=1;
	      } 
	    }   
	}  
	printf("This Loop vertices idx  %d  ( %6.2f  %6.2f )\n", 
	                               cET->vertexT->diskT->idx,
				       cET->vertexT->diskT->vertex->geo.vg2D->vtx.vtX,
				       cET->vertexT->diskT->vertex->geo.vg2D->vtx.vtY
				       );
	cET = cET->next;
     } while (cET != fET);

  }
  if(numberOfBifurP != numOfEnds-2 )
   {
     printf("Problem with bifur  %d\n", numberOfBifurP);
     exit(0);
   }

    *dstErr = errNum;
  return (BifurcatePointsIndex);
}

/*!
* \ingroup      WlzFeatures
* \return			none
* \brief	get the standard sample points in a specified loop points
*               continutity. 
*
* \param	AllVerticesInThisStartLoop  input pointer points to the full vertices
* \param	i_s		            input index to indicate from which position we
*                                           want to extract the points (input)
* \param        StandSampleP[NumberToTrack] output sampled points
* \param       *dstErr			    Woolz err number pointer.
*/
static void GetSamplePointsFromAllVertices(   WlzVertexP   sLoopVData, 
                                              int          i_s,
					      int          subSubSectionLength_L,
                                              WlzDVertex2  StandSampleP[NumberToTrack],
                                              WlzErrorNum *dstErr )
{
    int   i, k;
    int delta;
    WlzErrorNum errNum = WLZ_ERR_NONE;

    /*  first  */

    /* get the stand sample points first */
    delta = subSubSectionLength_L/NumberToTrack;
    k = i_s;
    for(i=0; i<NumberToTrack; i++)
    {
         StandSampleP[i].vtX =  (sLoopVData.d2 + k )->vtX;
         StandSampleP[i].vtY =  (sLoopVData.d2 + k )->vtY;
	 k +=  delta;
    }

}

/*!
* \return       none.	
* \ingroup	WlzFeatures
* \brief	out put points on surface, in and out of the surface
*		from a 3D GM.
* \param       *surfacefp		Given File pointer points to the file to output points on surface.
* \param       *infp		        Given File pointer points to the file to output points in surface.
* \param       *outfp		        Given File pointer points to the file to output points out surface.
* \param       *SurfacePoints           Given pointer points to the surface points.
* \param        nOfVOnTheS              Given the number of points on the surface. 
* \param       *InPoints                Given pointer points to the surface points.
* \param        nOfVWithinTheS         Given the number of points on the surface. 
* \param       *OutPoints               Given pointer points to the points in surface.
* \param        nOfVOutTheS          Given the number of points outside  the surface. 
* \todo         
*/
static void	outputSurfaceInAndOutPoints(int j_shell,
					    int i_loop,
					    int i_sec,
					    int nOfFiles,
					    int           *numOfPInS,
                                            WlzIVertex2   *sectionData[MaxNumOfFiles],
					    FILE          *surfacefp,
					    char          *surfaceStr,
					    FILE          *infp,
					    char          *inStr,
					    FILE          *outfp,
					    char          *outStr,
		                            WlzDVertex3   *SurfacePoints, 
					    int            nOfVOnTheS,
					    WlzDVertex3   *InPoints,      
					    int            nOfVWithinTheS,
					    WlzDVertex3   *OutPoints,     
					    int            nOfVOutTheS,
					    WlzErrorNum *dstErr )
{
     WlzErrorNum errNum = WLZ_ERR_NONE;
     char    fullnameS[90];
     char    fullnameI[90];
     char    fullnameO[90];
     char    fullname_section[90];
     char    fullname_sectionR[90];
     char    fullname_sectionS[90];
     char    fullname_sectionT[90];
     char    fStr[80];
     char    tempStr[3];
     int     i, k;
     FILE   *secfp = NULL;

     strcpy( fStr, "shell_");
     sprintf(tempStr,"%d" ,j_shell);
     strcat(fStr,  tempStr  );

     strcat(fStr,  "_loop_" );
     sprintf(tempStr,"%d" ,i_loop);
     strcat(fStr,  tempStr  );
     
     strcat(fStr,  "_Section_" );
     sprintf(tempStr,"%d" ,i_sec);
     strcat(fStr,  tempStr  );
     strcat(fStr, "_");

     strcpy(fullnameS, fStr);
     strcpy(fullnameI, fStr);
     strcpy(fullnameO, fStr);
     strcpy(fullname_section,fStr);
     
     strcat(fullname_section,"_file_");
      
     strcat(fullnameS, surfaceStr);
     strcat(fullnameI, inStr);
     strcat(fullnameO, outStr);
   

     /* printf("%s\n",fullnameS); */
     /* printf("%s\n",fullnameI); */
     /* printf("%s\n",fullnameO); */

  if((infp = fopen(fullnameI, "w")) == NULL )
  {
     errNum =  WLZ_ERR_FILE_OPEN;
  }
  if((outfp = fopen(fullnameO, "w")) == NULL )
  {
     errNum =  WLZ_ERR_FILE_OPEN;
  }
  if((surfacefp = fopen(fullnameS, "w")) == NULL )
  {
     errNum =  WLZ_ERR_FILE_OPEN;
  }
  if(errNum == WLZ_ERR_NONE)
  {

       	     for(k=0; k<nOfVOnTheS; k++)
       	     {
                 fprintf(surfacefp, "%6.2f  %6.2f   %6.2f\n",
		 SurfacePoints[k].vtX, 
	         SurfacePoints[k].vtY, 
	         SurfacePoints[k].vtZ ); 

       	     }
	     
       	     for(k=0; k<nOfVWithinTheS; k++)
       	     {
                 fprintf(infp, "%6.2f  %6.2f   %6.2f\n",
		 InPoints[k].vtX, 
	         InPoints[k].vtY, 
	         InPoints[k].vtZ ); 
       	     }

       	     for(k=0; k<nOfVOutTheS; k++)
       	     {
                 fprintf(outfp, "%6.2f  %6.2f   %6.2f\n",
		 OutPoints[k].vtX, 
	         OutPoints[k].vtY, 
	         OutPoints[k].vtZ ); 
       	     }

             for(k=0; k<nOfFiles; k++)
	     {
	        strcpy(fullname_sectionR, fullname_section);
	        strcpy(fullname_sectionT, fullname_section);
	        strcpy(fullname_sectionS, fullname_section);
                sprintf(tempStr,"%d" ,k);
                strcat(fullname_sectionR,tempStr);
		strcat(fullname_sectionS,tempStr);
		strcat(fullname_sectionT,tempStr);
		strcat(fullname_sectionR,"_forJavaDraw");
		strcat(fullname_sectionS,"_forContour");
		strcat(fullname_sectionT,"_forTie");
			
                if((secfp = fopen(fullname_sectionR, "w")) == NULL )
                {
                      errNum =  WLZ_ERR_FILE_OPEN;
                }
                if(errNum == WLZ_ERR_NONE)
		{
		    for(i=0; i< *(numOfPInS + k); i++)
		    {
                       fprintf(secfp, "%d  %d\n",
		               (sectionData[k]+i)->vtX, 
		               (sectionData[k]+i)->vtY 
	                      ); 
		    }
		    if(secfp)
		    {
		       fclose(secfp);
                       secfp = NULL;
		    }
		}
		
                if((secfp = fopen(fullname_sectionS, "w")) == NULL )
                {
                      errNum =  WLZ_ERR_FILE_OPEN;
                }
                if(errNum == WLZ_ERR_NONE)
		{
		    for(i=0; i< *(numOfPInS + k) - 1; i++)
		    {
                       fprintf(secfp, "%d  %d  %d  %d\n",
		               (sectionData[k]+i)->vtX, 
		               (sectionData[k]+i)->vtY, 
			       (sectionData[k]+i+1)->vtX, 
		               (sectionData[k]+i+1)->vtY 
	                      ); 
		    }
		    if(secfp)
		    {
		       fclose(secfp);
                       secfp = NULL;
		    }
		}

                if((secfp = fopen(fullname_sectionT, "w")) == NULL )
                {
                      errNum =  WLZ_ERR_FILE_OPEN;
                }
                if(errNum == WLZ_ERR_NONE)
		{
		    for( i= ( *(numOfPInS + k) + 1)/3; i< *(numOfPInS + k); i += ( *(numOfPInS + k) + 1)/3 + 1 )
		    {
                       fprintf(secfp, "%d  %d\n",
		               (sectionData[k]+i)->vtX, 
		               (sectionData[k]+i)->vtY 
	                      ); 
		    }
		    if(secfp)
		    {
		       fclose(secfp);
                       secfp = NULL;
		    }
		}


	     }

   }


      if(surfacefp)
      {
        fclose(surfacefp);
	surfacefp = NULL;
      }
      
      if(infp)
      {
        fclose(infp);
	infp = NULL;
      }
      
      if(outfp)
      {
        fclose(outfp);
	outfp = NULL;
      }
   *dstErr = errNum;
      
}

/*!
* \return	the pure the Geometry Model to a simple case by cutting small braches.
* \ingroup	WlzFeatures
* \brief	Allocates a buffer which it fills with the vertices
*		from a 3D GM.
* \param	gMC			Given model.
* \param	LoopIdx                 Given Loop index.
* \param       *EndsIdx                 Given the indices of ends points in this Loop
* \param       *bifurIdx                Given the indices of bifurcation points in this Loop
* \param	dstErr			Destination error pointer,
*					may be NULL.
* \todo        old vertion will be delet in the future 
*/
static void GetInOutAndFromSurfacePoints( 
				   int          numberOfPixelsZ,
                                   WlzDVertex2  StandSampleP[NumberToTrack],
                                   double       distance,
				   WlzDVertex3 *InPoints,
				   int         *nOfVWithinTheS,
				   WlzDVertex3 *OutPoints,
				   int         *nOfVOutTheS,
				   WlzErrorNum *dstErr )
{
    int i, j, k, m, delta, idx, ntV, loopDirection, di, num;
    int inN, outN;

    WlzErrorNum errNum = WLZ_ERR_NONE;
    WlzDVertex2         segV[3],  fPoint, *vNorm=NULL, *inPoints=NULL, *outPoints=NULL;
    inN  = *nOfVWithinTheS;
    outN = *nOfVOutTheS;
 
     /*  calculate the in and out surface points  */
	    
     idx = 0;
     i   = 0;
     num = 2;
     delta = NumberToTrack/num;
     
     if(   (  (vNorm     = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) * NumberToTrack)) == NULL) ||
           (  (inPoints  = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) * NumberToTrack)) == NULL) ||
           (  (outPoints = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) * NumberToTrack)) == NULL)
     ) 
     {
               	errNum = WLZ_ERR_MEM_ALLOC;
     }
     if(errNum == WLZ_ERR_NONE)
     {
      do
      {
	      
	  segV[0].vtX = StandSampleP[i].vtX; 
	  segV[0].vtY = StandSampleP[i].vtY;
	  segV[1].vtX = StandSampleP[i+1].vtX; 
	  segV[1].vtY = StandSampleP[i+1].vtY;	      
	  segV[2].vtX = StandSampleP[i+2].vtX; 
	  segV[2].vtY = StandSampleP[i+2].vtY;
	  /*  printf("%lg   %lg\n", segV[1].vtX, segV[1].vtY ); */
	 *(vNorm + idx) = WlzVerticesNormTriple2(segV[0], segV[1], segV[2]);
	  /*  get the in and out points using the normal and points */
	  fPoint.vtX =  segV[1].vtX + distance * (vNorm+idx)->vtX;
	  fPoint.vtY =  segV[1].vtY + distance * (vNorm+idx)->vtY;
	  /*  check its out or in */
	  di = WlzGeomTriangleSnArea2(segV[0], segV[1], fPoint); 
	  if(   di >  0. ) 
	  {
	     /*  fPoint is on the left of line segV[0]->segV[1] */
	     (inPoints +idx)->vtX =  fPoint.vtX; 
	     (inPoints +idx)->vtY =  fPoint.vtY;
	     
	     (outPoints +idx)->vtX =   segV[1].vtX - distance * (vNorm+idx)->vtX;
	     (outPoints +idx)->vtY =   segV[1].vtY - distance * (vNorm+idx)->vtY;
	     i += delta;
	     idx++;	     
  
	  }
	  else if(  di <  0. ) 
	  {
	     (outPoints +idx)->vtX =  fPoint.vtX; 
	     (outPoints +idx)->vtY =  fPoint.vtY;
	     
	     (inPoints +idx)->vtX =   segV[1].vtX - distance * (vNorm+idx)->vtX;
	     (inPoints +idx)->vtY =   segV[1].vtY - distance * (vNorm+idx)->vtY;
             i += delta;
	     idx++;	
	  }
	  else
	  {
              printf("On the surface!!!!!!!!!!!!! not counted it\n");
	      i++;
	      if(i + 3 > NumberToTrack)
	          break;

	  }
          /* } while(idx < NumberToTrack/4 ); */
     } while(idx < 1 );
			
     /* store the in points  */
     if(errNum == WLZ_ERR_NONE)
     {
	   for(i=0; i< idx; i++)
	   {

	       (InPoints + i + inN )->vtX =   (inPoints + i)->vtX;
	       (InPoints + i + inN )->vtY =   (inPoints + i)->vtY;
	       (InPoints + i + inN )->vtZ =   (double)numberOfPixelsZ;
	   }
	   inN += idx;
     }

     /* store the out points */
     if(errNum == WLZ_ERR_NONE)
     {
	   for(i=0; i< idx; i++)
	   {
	       (OutPoints + i + outN )->vtX =   (outPoints + i)->vtX;
	       (OutPoints + i + outN )->vtY =   (outPoints + i)->vtY;
	       (OutPoints + i + outN )->vtZ =   (double)numberOfPixelsZ;
	   }    
	   outN += idx;
     }
   }  
      if(vNorm)
      {   
           /* if(vNorm) */
              AlcFree(vNorm);
	      vNorm     = NULL;

      }	      
      if(inPoints)
      {
           /* if(inPoints) */
              AlcFree(inPoints);
              inPoints  = NULL;

      }
      if(outPoints)
      {
          /* if(outPoints) */
             AlcFree(outPoints);
      	     outPoints = NULL;

      }
      *nOfVOutTheS = outN;
      *nOfVWithinTheS  = inN;

}

/*!
* \return	the pure the Geometry Model to a simple case by cutting small braches.
* \ingroup	WlzFeatures
* \brief	Allocates a buffer which it fills with the vertices
*		from a 3D GM.
* \param	gMC			Given model.
* \param	dstErr			Destination error pointer,
*					may be NULL.
* \todo         haven't finished
*/
static void WlzGMModelPure( WlzGMModel *gMC, WlzErrorNum *dstErr)
{
   int             numOfShellsInStartGM;
   int           numOfLoopsInStartShell;
   int                   *startShellIdx;
   int                    *startLoopIdx;
   int                    *EndsIdx;
   int                    *bifurIdx;
   int                              i,j;
   int                        numOfEnds;
   WlzErrorNum                   errNum = WLZ_ERR_NONE;
   
        /*  get the number of Shells in this model */
        numOfShellsInStartGM  = NumberOfShellsAboutGM(gMC);
        /*  get the Index for each Shell; */
        startShellIdx = WlzShellsIndexAboutGM(gMC, numOfShellsInStartGM, &errNum);
        
	if(errNum == WLZ_ERR_NONE)
	{
          for(i=0; i<numOfShellsInStartGM; i++)
	  {
	     /*  get number loops in this shell */
             numOfLoopsInStartShell  =  NumberOfLoopsInThe_nThShellOfGM(gMC,  startShellIdx[i] );
	     printf("The  %d th Shell\n",  i );

             /*  get the index for the loops */
             startLoopIdx            =  WlzLoopIndexAboutGM(gMC, startShellIdx[i], numOfLoopsInStartShell, &errNum);
	     
	     if(errNum == WLZ_ERR_NONE )
	     {
                for(j=0; j<numOfLoopsInStartShell; j++)
		{
		   /* outputVerticesInThisLoop(gMC, startLoopIdx[j], &errNum); */
		   if( errNum != WLZ_ERR_NONE )
		          break;
	           /*  check is it a cycle */
		   if(!IsACycle(gMC, startLoopIdx[j], &errNum))
		   {
                        if(errNum == WLZ_ERR_NONE )
			{
			   numOfEnds = HowManyEndsInTheNonCycleLine( gMC, 
                                                            startLoopIdx[j], 
                                                            &errNum );
			  if( ( numOfEnds > 2 ) && (errNum == WLZ_ERR_NONE) ) 
			  {
                             /* pure this cycle */
                               EndsIdx =   GetTheEndsPointsIndexInTheNonCycleLine( gMC, 
                                           startLoopIdx[j],
					   numOfEnds,
                                           &errNum );
					   
			       bifurIdx = GetThebifurcatePointsIndexInTheNonCycleLine(   gMC, 
                                           startLoopIdx[j],
					   numOfEnds,
                                           &errNum );
			       pureThisGM(gMC, startLoopIdx[j], numOfEnds, EndsIdx, bifurIdx, &errNum);	
			       
			       /* if(EndsIdx) */
                                  AlcFree(EndsIdx);
			      /*  if(bifurIdx)	   */
			          AlcFree(bifurIdx);

			  }   
			     
			}

		   }

		}
	     }
	     /* if(startLoopIdx) */
	        AlcFree(startLoopIdx);  

	  }

	}
        /* if(startShellIdx) */
	    AlcFree(startShellIdx);
	
	if(errNum != WLZ_ERR_NONE)
	        *dstErr = errNum;
}

/*!
* \return	none.
* \ingroup	WlzFeatures
* \brief	output vertices of a specific  loop.
*		from a 3D GM.
* \param	gMC			Given model.
* \param	LoopIdx                 Given Loop index.
* \param	dstErr			Destination error pointer,
*					may be NULL.
* \todo         
*/
static void outputVerticesInThisLoop(WlzGMModel *gMC, int LoopIdx, WlzErrorNum *dstErr)
{
   int numberOfEnds, numBAgain;
   WlzErrorNum errNum = WLZ_ERR_NONE;
   int nVertices;
   AlcVector          *vec;
   WlzGMLoopT         *cLT;
   WlzGMShell         *cS;
   WlzGMEdgeT         *cET, *fET;
  
   numberOfEnds = 0;
   numBAgain    = 0;
   vec    =  gMC->res.loopT.vec;

   /* get the loop corresponding to the LoopIdx of the model. */
    cLT = (WlzGMLoopT *) AlcVectorItemGet(vec, LoopIdx);

   /* For each edge topology element of the model. */
    cET = fET = cLT->edgeT;

   /* For each loop topology element of the model. */
    do
    {
	cET = cET->next;
	/*
	printf("%6.2f   %6.2f\n",cET->vertexT->diskT->vertex->geo.vg2D->vtx.vtX,
	                         cET->vertexT->diskT->vertex->geo.vg2D->vtx.vtY);
				 */
	printf("%d   %d\n",(int) cET->vertexT->diskT->vertex->geo.vg2D->vtx.vtX,
	                   (int) cET->vertexT->diskT->vertex->geo.vg2D->vtx.vtY);
    } while (cET != fET);


}

/*!
* \return	the pure the Geometry Model to a simple case by cutting small braches.
* \ingroup	WlzFeatures
* \brief	Allocates a buffer which it fills with the vertices
*		from a 3D GM.
* \param	gMC			Given model.
* \param	LoopIdx                 Given Loop index.
* \param       *EndsIdx                 Given the indices of ends points in this Loop
* \param       *bifurIdx                Given the indices of bifurcation points in this Loop
* \param	dstErr			Destination error pointer,
*					may be NULL.
* \todo         haven't finished
*/
static void pureThisGM(WlzGMModel *gMC, 
				int LoopIdx, 
				int numOfEnds, 
				int *EndsIdx, 
				int *bifurIdx, 
				WlzErrorNum *dstErr)
{
   int      numberOfEnds, numBAgain;
   WlzErrorNum       errNum = WLZ_ERR_NONE;
   int                 nVertices;
   int      i, j, d1, d2, d3, it;
   int      dx, dy, dz;
   int      *length;
   AlcVector          *vec;
   WlzGMLoopT         *cLT;
   WlzGMShell         *cS;
   WlzGMEdgeT         *cET, *fET, *aET, *bET, *ffET;

    vec          =  gMC->res.loopT.vec;
    
    /* get the loop corresponding to the LoopIdx of the model. */
    cLT = (WlzGMLoopT *) AlcVectorItemGet(vec, LoopIdx);

    /* For each edge topology element of the model. */
    cET = fET = cLT->edgeT;
   /*
     delete bifurcation one by one but from the shortest to the longest */
   /*  
   for(i=0; i<numOfEnds; i++)
   {
      do
      {
  	  cET = cET->next;
	  if(    cET->vertexT->diskT->vertex->idx  ==   *(bifurIdx+i)  )
	       break;
      } while (cET != fET);
      */
      /*  check its length to the nearst bifurcation points */
      /*
      aET = ffET = cET->next;
      it = 0;
      do
      {
            it++;
            aET = aET->next;
	    for(j=0; j<numOfEnds-1; j++)
	    {
                  if(    aET->vertexT->diskT->vertex->idx  ==   *(bifurIdx+j)  )
		  {
		          break;
		  }  
	    }

      } while ( aET != ffET );

   }
   


   */
   
   for(i=0; i<numOfEnds-2; i++)
   {

      do
      {
  	  cET = cET->next;
	  if(    cET->vertexT->diskT->vertex->idx  ==   *(bifurIdx+i)  )
	       break;
      } while (cET != fET);
      /*  find the shortest braches from the bifurcation points */
      it = 0;
      d1 = 5;
      aET = cET->opp->next;
      bET = ffET = cET;
      do
      {
          /*  walking along three directions */
  	  cET = cET->next;
	  dx  = cET->vertexT->diskT->vertex->idx;
	  for(j=0; j<numOfEnds; j++)
	  {
	     if(  dx  ==   *(EndsIdx+j)  )
	     {
	          d1 = 0; 
	          break;
	     }	  
	  }
	  

	  if(d1 < 5)
	  {
	       it++;
	       break;
	  }     
  	  aET = aET->next;
	  dy = aET->vertexT->diskT->vertex->idx;
	  for(j=0; j<numOfEnds; j++)
	  {
	     if(  dy  ==   *(EndsIdx+j)  )
	     {
	          d1 = 1;
	          break;
	     }	  
	  }

	  
	  if(d1 < 5)
	  {
	       it++;
	       break;
	  }     
  
  	  bET = bET->prev;
	  dz = bET->vertexT->diskT->vertex->idx;
	  for(j=0; j<numOfEnds; j++)
	  {
	     if(   dz  ==   *(EndsIdx+j)  )
	     {
	          d1 = 2;
	          break;
	     }	  
	  }

	  if(d1 < 5)
	  {    
	       it++;
	       break;
	  }     
		  
	  it++;
      } while (bET != ffET);
     

      for(j=it; j>1; j--)
      {
            cET--;
	    aET--;
	    bET++;
      }
      if( d1== 0 )
      {
             errNum = WlzGMModelDeleteV( gMC, cET->vertexT->diskT->vertex);
	     if(it > 1)
	         errNum = WlzGMModelDeleteS( gMC, cET->parent->parent);

      }
      else if( d1== 1 )
      {
             errNum = WlzGMModelDeleteV( gMC, aET->vertexT->diskT->vertex);
	     if(it > 1)
	       errNum = WlzGMModelDeleteS( gMC, aET->parent->parent);
      }
      else if( d1== 2 )
      {
             errNum = WlzGMModelDeleteV( gMC, bET->vertexT->diskT->vertex);
	     if(it > 1)
	        errNum = WlzGMModelDeleteS( gMC, bET->parent->parent);
      }



   }
  
   if(errNum != WLZ_ERR_NONE)
        *dstErr = errNum;



}

/*!
* \return       none.	
* \ingroup	WlzFeatures
* \brief	Calculate In and Out points by using 3D normal:
*		from a 3D GM.
* \param        numOf2DWlzFiles         Given number of 2D Woolz contour files .
* \param       *SurfacePoints           Given pointer points to the surface points.
* \param        nOfVOnTheS              Given the number of points on the surface. 
* \param       *InPoints                Given pointer points to the surface points.
* \param        nOfVWithinTheS          Given the number of points on the surface. 
* \param       *OutPoints               Given pointer points to the points in surface.
* \param        nOfVOutTheS          Given the number of points outside  the surface. 
* \todo         
*/
static void RecalInOutPBy3DNorm(
                                 int            numOf2DWlzFiles,
				 double         distance,
		                 WlzDVertex3   *SurfacePoints, 
				 int            nOfVOnTheS,
				 WlzDVertex3   *InPoints,      
				 int           *nOfVWithinTheS,
				 WlzDVertex3   *OutPoints,     
				 int           *nOfVOutTheS )
{
   int i, j, k;
   int nn,mm,ll;
   int l, delta, halfDelta, count, num, lastline, lines1, lines2, nW, nt;
   int jumpLines, ik, mid, last;
   WlzDVertex3 p[4], p1[3], p2[3], normal, inOrOut[3][2], InStand[100];
   WlzDVertex3 tpv1[NumberToTrack*numOf2DWlzFiles], tpv2[numOf2DWlzFiles], tpv3[numOf2DWlzFiles];
   
   double volume;
        num       =  3;
	delta     =  NumberToTrack/num;
	mid       =  NumberToTrack/2;
	last      =  NumberToTrack - 1;
	halfDelta =  delta/2;
	k         =  NumberToTrack/num;
	l         =  NumberToTrack - delta - 1;
	count     =  0;
	jumpLines =  3; /* numOf2DWlzFiles/4; */
        lastline  =  numOf2DWlzFiles - jumpLines - 1;
	lines1    =  (numOf2DWlzFiles-1)/2;

        /*  reorganize it: */
	ll = 0;
	for(nn=lines1-1; nn>=0; nn--)
	{
	  mm  =   nn * NumberToTrack;
	  
          WlzValueCopyDVertexToDVertex3( (tpv1 + (ll * NumberToTrack)),     ( SurfacePoints + mm ), NumberToTrack); 
	  ll++;
	}

 	for(nn=0; nn<lines1; nn++)
	{
	  mm  =   nn * NumberToTrack;
          WlzValueCopyDVertexToDVertex3(  ( SurfacePoints + mm ),  (tpv1 + mm), NumberToTrack); 
	}




      /*  also re Order the in and out  */
       if( *nOfVWithinTheS == numOf2DWlzFiles )
       { 
       
          ll = 0;
          for(nn=lines1-1; nn>=0; nn--)
	  {
             WlzValueCopyDVertexToDVertex3( (tpv2 + ll) ,     ( InPoints  + nn ), 1);
             WlzValueCopyDVertexToDVertex3( (tpv3 + ll) ,     ( OutPoints + nn ), 1);
	     ll++;
	  }
	  
	  for(nn=0; nn<lines1; nn++)
	  {
             WlzValueCopyDVertexToDVertex3( ( InPoints  + nn ), (tpv2 + nn) ,  1);
             WlzValueCopyDVertexToDVertex3( ( OutPoints + nn ), (tpv3 + nn) ,  1);
	  }

       }

     /*  as a reference: */
       WlzValueCopyDVertexToDVertex3(InStand,  InPoints, *nOfVWithinTheS);

       nW = *nOfVWithinTheS;


     /*  the simplest way!!! */
     /*
     for(i=jumpLines; i< numOf2DWlzFiles; i++)
     {
        lines1  =   i                *  NumberToTrack;
	lines2  = lines1 + jumpLines *  NumberToTrack;
	if(   ( i + jumpLines - 1 )  > numOf2DWlzFiles )
	    break;
        if( ( i% ( jumpLines - 1 ) == 0 ) && (i != 0) )
	    i++;
	*/    
        /*  get 3 points from each line */
	/*
          WlzValueCopyDVertexToDVertex3(p1,     ( SurfacePoints + lines1 ), 1); 
          WlzValueCopyDVertexToDVertex3(&p1[1],  ( SurfacePoints + lines1 + mid ), 1); 
          WlzValueCopyDVertexToDVertex3(&p1[2],  ( SurfacePoints + lines1 + last), 1); 
	  
          WlzValueCopyDVertexToDVertex3(p2,     ( SurfacePoints + lines2 ), 1); 
          WlzValueCopyDVertexToDVertex3(&p2[1],  ( SurfacePoints + lines2 + mid ), 1); 
          WlzValueCopyDVertexToDVertex3(&p2[2],  ( SurfacePoints + lines2 + last), 1); 
	  */
        /*    p1 ---- p1[1] ---- p1[2] */
	/*  */
	/*    p2 ---- p2[1] ---- p2[2] */

	/*  calculate the normal p1 p1[1] p1			 */
	/*
	  normal = WlzGeomTriangleNormal(p1[0], p1[1], p2[0]);
	  if( !IsDVertexZero(normal) )
	  {
	     break;
	  }   
     }
     */
     /*  get output: */
     /*
     for(i=0; i< numOf2DWlzFiles; i +=2)
     {
        
        lines1  =  i   *  NumberToTrack;
	for(j=0; j<NumberToTrack; j += 3)
	{
	   k = lines1 + j;
           (InPoints + count)->vtX  =  ( SurfacePoints + k )->vtX + distance * normal.vtX;
           (InPoints + count)->vtY  =  ( SurfacePoints + k )->vtY + distance * normal.vtY;
           (InPoints + count)->vtZ  =  ( SurfacePoints + k )->vtZ + distance * normal.vtZ;

           (OutPoints + count)->vtX  =  ( SurfacePoints + k )->vtX - distance * normal.vtX;
           (OutPoints + count)->vtY  =  ( SurfacePoints + k )->vtY - distance * normal.vtY;
           (OutPoints + count)->vtZ  =  ( SurfacePoints + k )->vtZ - distance * normal.vtZ;

	   count++;
        }
     */  
     /*  recalculate the in points and out points by scan line by line */
     for(i=0; i< numOf2DWlzFiles; i++)
     {
        lines1  =   i                *  NumberToTrack;
	lines2  = lines1 + jumpLines *  NumberToTrack;
        if( ( i% ( jumpLines - 1 ) == 0 ) && (i != 0) )
	    i++;
	if(   ( i + jumpLines - 1 )  >= numOf2DWlzFiles )
	    break;
        /*  get 3 points from each line */
          WlzValueCopyDVertexToDVertex3(p1,     ( SurfacePoints + lines1 ), 1); 
          WlzValueCopyDVertexToDVertex3(&p1[1],  ( SurfacePoints + lines1 + mid ), 1); 
          WlzValueCopyDVertexToDVertex3(&p1[2],  ( SurfacePoints + lines1 + last), 1); 
	  
          WlzValueCopyDVertexToDVertex3(p2,     ( SurfacePoints + lines2 ), 1); 
          WlzValueCopyDVertexToDVertex3(&p2[1],  ( SurfacePoints + lines2 + mid ), 1); 
          WlzValueCopyDVertexToDVertex3(&p2[2],  ( SurfacePoints + lines2 + last), 1); 


        /*    p1 ---- p1[1] ---- p1[2] */
	/*  */
	/*    p2 ---- p2[1] ---- p2[2] */

	/*  calculate the normal p1 p1[1] p1			 */
	  normal = WlzGeomTriangleNormal(p1[0], p1[1], p2[0]);
	  if( !IsDVertexZero(normal) )
	  {
	     nt = i;
	     if(nt > nW )
	         nt = nW -1;
  	     WlzValueCopyDVertexToDVertex3( &p[3], ( InStand + nt), 1 ); 
             volume = SixTimesOfVoulumeOfTetraHedron( p1[0], p1[1], p2[0], p[3]); 
             if( volume != 0.)
	     {
	       inOrOut[0][0].vtX = p1[0].vtX + distance * normal.vtX;
	       inOrOut[0][0].vtY = p1[0].vtY + distance * normal.vtY;
	       inOrOut[0][0].vtZ = p1[0].vtZ + distance * normal.vtZ;
	       
	       inOrOut[0][1].vtX = p1[0].vtX - distance * normal.vtX;
	       inOrOut[0][1].vtY = p1[0].vtY - distance * normal.vtY;
	       inOrOut[0][1].vtZ = p1[0].vtZ - distance * normal.vtZ;
	       
	       inOrOut[1][0].vtX = p1[1].vtX + distance * normal.vtX;
	       inOrOut[1][0].vtY = p1[1].vtY + distance * normal.vtY;
	       inOrOut[1][0].vtZ = p1[1].vtZ + distance * normal.vtZ;
	       
	       inOrOut[1][1].vtX = p1[1].vtX - distance * normal.vtX;
	       inOrOut[1][1].vtY = p1[1].vtY - distance * normal.vtY;
	       inOrOut[1][1].vtZ = p1[1].vtZ - distance * normal.vtZ;
	       
	       inOrOut[2][0].vtX = p2[0].vtX + distance * normal.vtX;
	       inOrOut[2][0].vtY = p2[0].vtY + distance * normal.vtY;
	       inOrOut[2][0].vtZ = p2[0].vtZ + distance * normal.vtZ;
	       
	       inOrOut[2][1].vtX = p2[0].vtX - distance * normal.vtX;
	       inOrOut[2][1].vtY = p2[0].vtY - distance * normal.vtY;
	       inOrOut[2][1].vtZ = p2[0].vtZ - distance * normal.vtZ;
				
		if(volume > 0)
		{
                       (InPoints + count)->vtX  = inOrOut[0][0].vtX;
                       (InPoints + count)->vtY  = inOrOut[0][0].vtY;
                       (InPoints + count)->vtZ  = inOrOut[0][0].vtZ;

                       (OutPoints + count)->vtX = inOrOut[0][1].vtX;
                       (OutPoints + count)->vtY = inOrOut[0][1].vtY;
                       (OutPoints + count)->vtZ = inOrOut[0][1].vtZ;

		       count++;

                       (InPoints + count)->vtX  = inOrOut[1][0].vtX;
                       (InPoints + count)->vtY  = inOrOut[1][0].vtY;
                       (InPoints + count)->vtZ  = inOrOut[1][0].vtZ;

                       (OutPoints + count)->vtX = inOrOut[1][1].vtX;
                       (OutPoints + count)->vtY = inOrOut[1][1].vtY;
                       (OutPoints + count)->vtZ = inOrOut[1][1].vtZ;
                       
		       count++;

                       (InPoints + count)->vtX  = inOrOut[2][0].vtX;
                       (InPoints + count)->vtY  = inOrOut[2][0].vtY;
                       (InPoints + count)->vtZ  = inOrOut[2][0].vtZ;

                       (OutPoints + count)->vtX = inOrOut[2][1].vtX;
                       (OutPoints + count)->vtY = inOrOut[2][1].vtY;
                       (OutPoints + count)->vtZ = inOrOut[2][1].vtZ;

		       count++;


		}
		else
		{
                       (InPoints + count)->vtX  = inOrOut[0][1].vtX;
                       (InPoints + count)->vtY  = inOrOut[0][1].vtY;
                       (InPoints + count)->vtZ  = inOrOut[0][1].vtZ;

                       (OutPoints + count)->vtX = inOrOut[0][0].vtX;
                       (OutPoints + count)->vtY = inOrOut[0][0].vtY;
                       (OutPoints + count)->vtZ = inOrOut[0][0].vtZ;
		       
		        count++;

                       (InPoints + count)->vtX  = inOrOut[1][1].vtX;
                       (InPoints + count)->vtY  = inOrOut[1][1].vtY;
                       (InPoints + count)->vtZ  = inOrOut[1][1].vtZ;

                       (OutPoints + count)->vtX = inOrOut[1][0].vtX;
                       (OutPoints + count)->vtY = inOrOut[1][0].vtY;
                       (OutPoints + count)->vtZ = inOrOut[1][0].vtZ;
		       
                        count++;

		       
                       (InPoints + count)->vtX  = inOrOut[2][1].vtX;
                       (InPoints + count)->vtY  = inOrOut[2][1].vtY;
                       (InPoints + count)->vtZ  = inOrOut[2][1].vtZ;

                       (OutPoints + count)->vtX = inOrOut[2][0].vtX;
                       (OutPoints + count)->vtY = inOrOut[2][0].vtY;
                       (OutPoints + count)->vtZ = inOrOut[2][0].vtZ;

		       count++;

		}


           }
        }



	/*  calculate the normal p1[1] p1[2] p2[2] */
	  normal = WlzGeomTriangleNormal(p1[1], p1[2], p2[2]);
	  if( !IsDVertexZero(normal) )
	  {
	     	     nt = i;
	     if(nt > nW )
	         nt = nW -1;

  	     WlzValueCopyDVertexToDVertex3( &p[3], ( InStand + nt), 1 ); 
             volume = SixTimesOfVoulumeOfTetraHedron( p1[1], p1[2], p2[2], p[3]); 
             if( volume != 0.)
	     {
	       /*
	       inOrOut[0][0].vtX = p1[1].vtX + distance * normal.vtX;
	       inOrOut[0][0].vtY = p1[1].vtY + distance * normal.vtY;
	       inOrOut[0][0].vtZ = p1[1].vtZ + distance * normal.vtZ;
	       
	       inOrOut[0][1].vtX = p1[1].vtX - distance * normal.vtX;
	       inOrOut[0][1].vtY = p1[1].vtY - distance * normal.vtY;
	       inOrOut[0][1].vtZ = p1[1].vtZ - distance * normal.vtZ;
	       */
	       inOrOut[1][0].vtX = p1[2].vtX + distance * normal.vtX;
	       inOrOut[1][0].vtY = p1[2].vtY + distance * normal.vtY;
	       inOrOut[1][0].vtZ = p1[2].vtZ + distance * normal.vtZ;
	       
	       inOrOut[1][1].vtX = p1[2].vtX - distance * normal.vtX;
	       inOrOut[1][1].vtY = p1[2].vtY - distance * normal.vtY;
	       inOrOut[1][1].vtZ = p1[2].vtZ - distance * normal.vtZ;
	       
	       inOrOut[2][0].vtX = p2[2].vtX + distance * normal.vtX;
	       inOrOut[2][0].vtY = p2[2].vtY + distance * normal.vtY;
	       inOrOut[2][0].vtZ = p2[2].vtZ + distance * normal.vtZ;
	       
	       inOrOut[2][1].vtX = p2[2].vtX - distance * normal.vtX;
	       inOrOut[2][1].vtY = p2[2].vtY - distance * normal.vtY;
	       inOrOut[2][1].vtZ = p2[2].vtZ - distance * normal.vtZ;
				
		if(volume > 0)
		{
		      /*
                       (InPoints + count)->vtX  = inOrOut[0][0].vtX;
                       (InPoints + count)->vtY  = inOrOut[0][0].vtY;
                       (InPoints + count)->vtZ  = inOrOut[0][0].vtZ;

                       (OutPoints + count)->vtX = inOrOut[0][1].vtX;
                       (OutPoints + count)->vtY = inOrOut[0][1].vtY;
                       (OutPoints + count)->vtZ = inOrOut[0][1].vtZ;

		       count++;
                      */
                       (InPoints + count)->vtX  = inOrOut[1][0].vtX;
                       (InPoints + count)->vtY  = inOrOut[1][0].vtY;
                       (InPoints + count)->vtZ  = inOrOut[1][0].vtZ;

                       (OutPoints + count)->vtX = inOrOut[1][1].vtX;
                       (OutPoints + count)->vtY = inOrOut[1][1].vtY;
                       (OutPoints + count)->vtZ = inOrOut[1][1].vtZ;
                       
		       count++;

                       (InPoints + count)->vtX  = inOrOut[2][0].vtX;
                       (InPoints + count)->vtY  = inOrOut[2][0].vtY;
                       (InPoints + count)->vtZ  = inOrOut[2][0].vtZ;

                       (OutPoints + count)->vtX = inOrOut[2][1].vtX;
                       (OutPoints + count)->vtY = inOrOut[2][1].vtY;
                       (OutPoints + count)->vtZ = inOrOut[2][1].vtZ;

		       count++;

		}
		else
		{
		    /*
                       (InPoints + count)->vtX  = inOrOut[0][1].vtX;
                       (InPoints + count)->vtY  = inOrOut[0][1].vtY;
                       (InPoints + count)->vtZ  = inOrOut[0][1].vtZ;

                       (OutPoints + count)->vtX = inOrOut[0][0].vtX;
                       (OutPoints + count)->vtY = inOrOut[0][0].vtY;
                       (OutPoints + count)->vtZ = inOrOut[0][0].vtZ;
		       
		        count++;
		     */	

                       (InPoints + count)->vtX  = inOrOut[1][1].vtX;
                       (InPoints + count)->vtY  = inOrOut[1][1].vtY;
                       (InPoints + count)->vtZ  = inOrOut[1][1].vtZ;

                       (OutPoints + count)->vtX = inOrOut[1][0].vtX;
                       (OutPoints + count)->vtY = inOrOut[1][0].vtY;
                       (OutPoints + count)->vtZ = inOrOut[1][0].vtZ;
		       
                        count++;

		       
                       (InPoints + count)->vtX  = inOrOut[2][1].vtX;
                       (InPoints + count)->vtY  = inOrOut[2][1].vtY;
                       (InPoints + count)->vtZ  = inOrOut[2][1].vtZ;

                       (OutPoints + count)->vtX = inOrOut[2][0].vtX;
                       (OutPoints + count)->vtY = inOrOut[2][0].vtY;
                       (OutPoints + count)->vtZ = inOrOut[2][0].vtZ;

		       count++;

		}


           }
        }











       /*






         */
        /*  take three points from the surface: */
/*	for(j=0; j<l; j += halfDelta)
	{  */
	   /*  get the first point  */
	 /*  WlzValueCopyDVertexToDVertex3(p,  (SurfacePoints + lines1 + j + 1 ), 1);
	   */
	   /*  get the second point */
/*	   WlzValueCopyDVertexToDVertex3(&p[1],  (SurfacePoints + lines1 + delta - 1), 1);
*/	   /*  if( !IsDVertexEqual(p[0], p[1]) ) */

           /*  get the third point */
/*	   WlzValueCopyDVertexToDVertex3(&p[2],  (SurfacePoints + lines2 + j + halfDelta ), 1);
*/
	  /*  calculate the normal			 */
/*	   normal = WlzGeomTriangleNormal(p[0], p[1], p[2]);
 */ 
/*	  if( !IsDVertexZero(normal) )
	  {
*/	     /*  detremine the in and out points: */
/*	       nt = i;
	     if(nt > nW )
	        nt = nW - 1;
	     WlzValueCopyDVertexToDVertex3( &p[3], ( InStand + nt), 1 ); 
             volume = SixTimesOfVoulumeOfTetraHedron( p[0], p[1], p[2], p[3]); 
             if( volume != 0.)
	     {
	       inOrOut[0][0].vtX = p[0].vtX + distance * normal.vtX;
	       inOrOut[0][0].vtY = p[0].vtY + distance * normal.vtY;
	       inOrOut[0][0].vtZ = p[0].vtZ + distance * normal.vtZ;
	       
	       inOrOut[0][1].vtX = p[0].vtX - distance * normal.vtX;
	       inOrOut[0][1].vtY = p[0].vtY - distance * normal.vtY;
	       inOrOut[0][1].vtZ = p[0].vtZ - distance * normal.vtZ;
	       
	       inOrOut[1][0].vtX = p[1].vtX + distance * normal.vtX;
	       inOrOut[1][0].vtY = p[1].vtY + distance * normal.vtY;
	       inOrOut[1][0].vtZ = p[1].vtZ + distance * normal.vtZ;
	       
	       inOrOut[1][1].vtX = p[1].vtX - distance * normal.vtX;
	       inOrOut[1][1].vtY = p[1].vtY - distance * normal.vtY;
	       inOrOut[1][1].vtZ = p[1].vtZ - distance * normal.vtZ;
	       
	       inOrOut[2][0].vtX = p[2].vtX + distance * normal.vtX;
	       inOrOut[2][0].vtY = p[2].vtY + distance * normal.vtY;
	       inOrOut[2][0].vtZ = p[2].vtZ + distance * normal.vtZ;
	       
	       inOrOut[2][1].vtX = p[2].vtX - distance * normal.vtX;
	       inOrOut[2][1].vtY = p[2].vtY - distance * normal.vtY;
	       inOrOut[2][1].vtZ = p[2].vtZ - distance * normal.vtZ;
				
		if(volume > 0)
		{
                       (InPoints + count)->vtX  = inOrOut[0][0].vtX;
                       (InPoints + count)->vtY  = inOrOut[0][0].vtY;
                       (InPoints + count)->vtZ  = inOrOut[0][0].vtZ;

                       (OutPoints + count)->vtX = inOrOut[0][1].vtX;
                       (OutPoints + count)->vtY = inOrOut[0][1].vtY;
                       (OutPoints + count)->vtZ = inOrOut[0][1].vtZ;

		       count++;

                       (InPoints + count)->vtX  = inOrOut[1][0].vtX;
                       (InPoints + count)->vtY  = inOrOut[1][0].vtY;
                       (InPoints + count)->vtZ  = inOrOut[1][0].vtZ;

                       (OutPoints + count)->vtX = inOrOut[1][1].vtX;
                       (OutPoints + count)->vtY = inOrOut[1][1].vtY;
                       (OutPoints + count)->vtZ = inOrOut[1][1].vtZ;
                       
		       count++;

                       (InPoints + count)->vtX  = inOrOut[2][0].vtX;
                       (InPoints + count)->vtY  = inOrOut[2][0].vtY;
                       (InPoints + count)->vtZ  = inOrOut[2][0].vtZ;

                       (OutPoints + count)->vtX = inOrOut[2][1].vtX;
                       (OutPoints + count)->vtY = inOrOut[2][1].vtY;
                       (OutPoints + count)->vtZ = inOrOut[2][1].vtZ;

		       count++;


		}
		else
		{
                       (InPoints + count)->vtX  = inOrOut[0][1].vtX;
                       (InPoints + count)->vtY  = inOrOut[0][1].vtY;
                       (InPoints + count)->vtZ  = inOrOut[0][1].vtZ;

                       (OutPoints + count)->vtX = inOrOut[0][0].vtX;
                       (OutPoints + count)->vtY = inOrOut[0][0].vtY;
                       (OutPoints + count)->vtZ = inOrOut[0][0].vtZ;
		       
		        count++;

                       (InPoints + count)->vtX  = inOrOut[1][1].vtX;
                       (InPoints + count)->vtY  = inOrOut[1][1].vtY;
                       (InPoints + count)->vtZ  = inOrOut[1][1].vtZ;

                       (OutPoints + count)->vtX = inOrOut[1][0].vtX;
                       (OutPoints + count)->vtY = inOrOut[1][0].vtY;
                       (OutPoints + count)->vtZ = inOrOut[1][0].vtZ;
		       
                        count++;

		       
                       (InPoints + count)->vtX  = inOrOut[2][1].vtX;
                       (InPoints + count)->vtY  = inOrOut[2][1].vtY;
                       (InPoints + count)->vtZ  = inOrOut[2][1].vtZ;

                       (OutPoints + count)->vtX = inOrOut[2][0].vtX;
                       (OutPoints + count)->vtY = inOrOut[2][0].vtY;
                       (OutPoints + count)->vtZ = inOrOut[2][0].vtZ;

		       count++;

		}
*/	        /*  replace the normal with new one! */
/*	      }
	   }
        


        }
	*/
     }
     *nOfVWithinTheS   = count;
     *nOfVOutTheS      = count;



}

static int IsDVertexEqual(WlzDVertex3 p1, WlzDVertex3 p2)
{
   int equal;
   equal = 0;
   if( p1.vtX == p2.vtX )
   {
      if( (p1.vtY == p2.vtY) && (p1.vtZ == p2.vtZ) )
      {
         equal = 1;
      }
   }
   return equal;
}


static int IsDVertexZero(WlzDVertex3  normal)
{
   int isZero;
   isZero = 0;
   if( normal.vtX == 0. )
   {
      if( (normal.vtY == 0. ) && (normal.vtZ == 0.) )
      {
        isZero  = 1;
      }
   }
   return isZero;


}


static double SixTimesOfVoulumeOfTetraHedron(WlzDVertex3 p0, 
			WlzDVertex3 p1,
			WlzDVertex3 p2,
			WlzDVertex3 p3
)
{
   double vol;
   double ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz;

   ax = p0.vtX;
   ay = p0.vtY;
   az = p0.vtZ;
   
   bx = p1.vtX;
   by = p1.vtY;
   bz = p1.vtZ;
   
   cx = p2.vtX;
   cy = p2.vtY;
   cz = p2.vtZ;
   
   dx = p3.vtX;
   dy = p3.vtY;
   dz = p3.vtZ;

   vol = - az * by * cx+ay * bz * cx+az * bx * cy - ax * bz * cy
         - ay * bx * cz+ax * by * cz+az * by * dx - ay * bz * dx
	 - az * cy * dx+bz * cy * dx+ay * cz * dx - by * cz * dx
	 - az * bx * dy+ax * bz * dy+az * cx * dy - bz * cx * dy
	 - ax * cz * dy+bx * cz * dy+ay * bx * dz - ax * by * dz;
	 - ay * cx * dz+by * cx * dz+ax * cy * dz - bx * cy * dz;
		 
   return vol;

}


static void  outputThisSection( WlzVertexP AllVerticesInThisStartLoop, 
                                int          i_s, 
                                int          subSubSectionLength_L,
                                WlzErrorNum *dstErr )
{
  int i, k;
  k = i_s + subSubSectionLength_L;
  for(i=i_s; i<k; i++)
  {
     printf("%d %d\n", (int ) (AllVerticesInThisStartLoop.d2 + i)->vtX, 
                       (int ) (AllVerticesInThisStartLoop.d2 + i)->vtY );
  }


}


static void  outputThisSectionForView( WlzVertexP AllVerticesInThisStartLoop, 
                                int          i_s, 
                                int          subSubSectionLength_L,
                                WlzErrorNum *dstErr )
{
  int i, k;
  k = i_s + subSubSectionLength_L;
  for(i=i_s; i<k; i++)
  {
     printf("%lg %lg %lg\n",  (AllVerticesInThisStartLoop.d2 + i)->vtX, 
                        (AllVerticesInThisStartLoop.d2 + i)->vtY, 0 );
  }


}






/*!
* \return       none.	
* \ingroup	WlzFeatures
* \brief	out put points on surface, in and out of the surface
*		from a 3D GM.
* \param       *surfacefp		Given File pointer points to the file to output points on surface.
* \param       *infp		        Given File pointer points to the file to output points in surface.
* \param       *outfp		        Given File pointer points to the file to output points out surface.
* \param       *SurfacePoints           Given pointer points to the surface points.
* \param        nOfVOnTheS              Given the number of points on the surface. 
* \param       *InPoints                Given pointer points to the surface points.
* \param        nOfVWithinTheS         Given the number of points on the surface. 
* \param       *OutPoints               Given pointer points to the points in surface.
* \param        nOfVOutTheS          Given the number of points outside  the surface. 
* \todo         
*/
static void	outputSurfaceByPatch(int j_shell,
					    int i_loop,
					    int i_sec,
					    int nOfFiles,
					    int subSubSectionLength_L,
					    int           *numOfPInS,
                                            WlzIVertex2   *sectionData[MaxNumOfFiles],
					    FILE          *surfacefp,
					    char          *surfaceStr,
					    FILE          *infp,
					    char          *inStr,
					    FILE          *outfp,
					    char          *outStr,
		                            WlzDVertex3   *SurfacePoints, 
					    int            nOfVOnTheS,
					    WlzDVertex3   *InPoints,      
					    int            nOfVWithinTheS,
					    WlzDVertex3   *OutPoints,     
					    int            nOfVOutTheS,
					    WlzErrorNum *dstErr )
{
     WlzErrorNum errNum = WLZ_ERR_NONE;
     double  *xp,  *yp, *zp,  *ac,  *sig;
     double  *xpA,  *ypA, *zpA, *sigA;
     double **u,  **v,   *w1,   chisq, **fittedx, **fittedy;
     double **uA;
     double **orderedSurfaceX, **orderedSurfaceY, *inx;
     double  *Esx, *Esy, *Efx, *Efy;
     double   x0, x1, y0, y1, X0, Y0, X1, Y1;

     char    fullnameS[90];
     char    fullnameI[90];
     char    fullnameO[90];
     char    fullname_section[90];
     char    fullname_sectionBase[90];
     char    fullname_sectionR[90];
     char    fullname_sectionS[90];
     char    fullname_sectionT[90];
     char    fullname_sectionTC[90];
     char    fullname_sectionTN[90];
     char    fStr[80];
     char    tempStr[4];
     int     i, k, j, l, m, delta, dimV;
     int     polyDeg = 3;
     FILE   *secfp = NULL;
     
     delta = subSubSectionLength_L/NumberToTrack;

     strcpy( fStr, "shell_");
     sprintf(tempStr,"%d" ,j_shell);
     strcat(fStr,  tempStr  );

     strcat(fStr,  "_loop_" );
     sprintf(tempStr,"%d" ,i_loop);
     strcat(fStr,  tempStr  );
     
     strcat(fStr,  "_Section_" );
     sprintf(tempStr,"%d" ,i_sec);
     strcat(fStr,  tempStr  );
     strcat(fStr, "_");

     strcpy(fullnameS, fStr);
     strcpy(fullnameI, fStr);
     strcpy(fullnameO, fStr);
     strcpy(fullname_section,fStr);
     
     strcat(fullname_section,"file_");
      
     strcat(fullnameS, surfaceStr);
     strcat(fullnameI, inStr);
     strcat(fullnameO, outStr);
   

     /* printf("%s\n",fullnameS); */
     /* printf("%s\n",fullnameI); */
     /* printf("%s\n",fullnameO); */

  if((infp = fopen(fullnameI, "w")) == NULL )
  {
     errNum =  WLZ_ERR_FILE_OPEN;
  }
  if((outfp = fopen(fullnameO, "w")) == NULL )
  {
     errNum =  WLZ_ERR_FILE_OPEN;
  }
  if((surfacefp = fopen(fullnameS, "w")) == NULL )
  {
     errNum =  WLZ_ERR_FILE_OPEN;
  }
  if(errNum == WLZ_ERR_NONE)
  {

       	     for(k=0; k<nOfVOnTheS; k++)
       	     {
                 fprintf(surfacefp, "%6.2f  %6.2f   %6.2f\n",
		 SurfacePoints[k].vtX, 
	         SurfacePoints[k].vtY, 
	         SurfacePoints[k].vtZ ); 

       	     }
	     /*  using the above data to get interpolate surface here: */

	     /*  modeling the data using least square method  */

             /* allocate memory for edge (Ex,Ey) */

             Efx    = vector(1, nOfFiles); /*  used to stored sampled points */
             Efy    = vector(1, nOfFiles);

             Esx    = vector(1, nOfFiles); /*  used to stored sampled points */
             Esy    = vector(1, nOfFiles);
	     inx    = vector(1, nOfFiles);

             dimV   = NumberToTrack*delta;
             fittedx   = matrix(1,nOfFiles, 1, dimV);
             fittedy   = matrix(1,nOfFiles, 1, dimV);
             u   = matrix(1,nOfFiles, 1, polyDeg);
             v   = matrix(1,polyDeg, 1, polyDeg);
             w1  = vector(1,polyDeg);
  	     xp  = vector(1,nOfFiles);
  	     yp  = vector(1,nOfFiles);
  	     zp  = vector(1,nOfFiles);
   	     ac  = vector(1,polyDeg);
  	     sig = vector(1,nOfFiles);

             uA   = matrix(1,NumberToTrack, 1, polyDeg);
  	     xpA  = vector(1,NumberToTrack);
  	     ypA  = vector(1,NumberToTrack);
  	     zpA  = vector(1,NumberToTrack);
  	     sigA = vector(1,NumberToTrack);


          /*  make sure that the surface line is in the same direction: */
          orderedSurfaceX  = matrix(1,nOfFiles, 1, NumberToTrack);
          orderedSurfaceY  = matrix(1,nOfFiles, 1, NumberToTrack);

	  for(j=1; j<=nOfFiles; j++)
	  {
	     l =  (j-1) * NumberToTrack;
	     m =  j * NumberToTrack - 1;
	     x0 =  ( SurfacePoints[m].vtX - SurfacePoints[l].vtX );
	     inx[j] = 1;
	     if( x0 > 0)
	     {
	       for(k=1; k<= NumberToTrack; k++)
	       {
                  orderedSurfaceX[j][k] = SurfacePoints[l+k-1].vtX;
                  orderedSurfaceY[j][k] = SurfacePoints[l+k-1].vtY;
	       }
	     }  
	     else
	     {
	       inx[j] =  0;
	       for(k=1; k<= NumberToTrack; k++)
	       {
                  orderedSurfaceX[j][k] = SurfacePoints[m+1-k].vtX;
                  orderedSurfaceY[j][k] = SurfacePoints[m+1-k].vtY;
	       }
	     }  
	  }   
 
         /* model the data vertically now */
         if( nOfFiles > 1 )
	 {
          for(j=1; j<=NumberToTrack; j++)
	  {

	     for(k=1; k<=nOfFiles; k++)
	     {
		xp[k]  = (double) k;
		yp[k]  = orderedSurfaceX[k][j];
		sig[k] = 1.0;
	     }
             
	     /* for x */
	     svdfit(xp, yp, sig, nOfFiles, ac, polyDeg, u, v, w1, &chisq, funcs );
	     
             for(k = 1; k <= nOfFiles; k++)
             {
                fittedx[k][j] = polynomialF( xp[k], ac, polyDeg);
             }
	     /* for y */
	     
	     for(k=1; k<=nOfFiles; k++)
	     {
		xp[k]  = (double) k;
		zp[k]  = orderedSurfaceY[k][j];
		sig[k] = 1.0;
	     }
	     
	     svdfit(xp, zp, sig, nOfFiles, ac, polyDeg, u, v, w1, &chisq, funcs );
	     
             for(k = 1; k <= nOfFiles; k++)
             {
                fittedy[k][j] = polynomialF( xp[k], ac, polyDeg);
             }

          }
	 }
	 else if( nOfFiles == 1 )
	 {
          for(j=1; j<=NumberToTrack; j++)
	  {
	     for(k=1; k<=nOfFiles; k++)
	     {
	       fittedx[k][j] = orderedSurfaceX[k][j];
	       fittedy[k][j] = orderedSurfaceY[k][j];
	     }
	  }
	 }
 
          /*  now interpolate horizontally */
          for(j=1; j<=nOfFiles; j++)
	  {
             /* for x */
	     for(k=1; k<=NumberToTrack; k++)
	     {
		xpA[k]  = (double) ( (k-1) * delta + 1);
		ypA[k]  = fittedx[j][k];
		sigA[k] = 1.0;
	     }
             
	     svdfit(xpA, ypA, sigA, NumberToTrack, ac, polyDeg, uA, v, w1, &chisq, funcs );
	     
             for(k = 1; k <= dimV ; k++)
             {
                fittedx[j][k] = polynomialF( (double) k, ac, polyDeg);
             }
	     /* for y */
	     
	     for(k=1; k<= NumberToTrack; k++)
	     {
		xpA[k]  = (double) ( (k-1)*delta + 1 );
		zpA[k]  = fittedy[j][k];
		sigA[k] = 1.0;
	     }
	     
	     svdfit(xpA, zpA, sigA, NumberToTrack, ac, polyDeg, uA, v, w1, &chisq, funcs );
	     
             for(k = 1; k <= dimV; k++)
             {
                fittedy[j][k] = polynomialF( (double) k, ac, polyDeg);
             }


          }



             /*  interpolate finished here !!! */
	     
       	     for(k=0; k<nOfVWithinTheS; k++)
       	     {
                 fprintf(infp, "%6.2f  %6.2f   %6.2f\n",
		 InPoints[k].vtX, 
	         InPoints[k].vtY, 
	         InPoints[k].vtZ ); 
       	     }

       	     for(k=0; k<nOfVOutTheS; k++)
       	     {
                 fprintf(outfp, "%6.2f  %6.2f   %6.2f\n",
		 OutPoints[k].vtX, 
	         OutPoints[k].vtY, 
	         OutPoints[k].vtZ ); 
       	     }

	     /*  output data for other purpose: */

             for(k=0; k<nOfFiles; k++)
	     {
                sprintf(tempStr,"%d" ,k);
	        strcpy(fullname_sectionBase, fullname_section);
		strcat(fullname_sectionBase, tempStr);

	        strcpy(fullname_sectionR, fullname_sectionBase);
	        strcpy(fullname_sectionT, fullname_sectionBase);
	        strcpy(fullname_sectionS, fullname_sectionBase);
	        strcpy(fullname_sectionTC, fullname_sectionBase);
	        strcpy(fullname_sectionTN, fullname_sectionBase);
			
		strcat(fullname_sectionR,"_forJavaDraw");
		strcat(fullname_sectionS,"_forContour");
		strcat(fullname_sectionT,"_forTie");
		strcat(fullname_sectionTC,"_forTContour");
		strcat(fullname_sectionTN,"_TieWithD");

                if((secfp = fopen(fullname_sectionR, "w")) == NULL )
                {
                      errNum =  WLZ_ERR_FILE_OPEN;
                }
                if(errNum == WLZ_ERR_NONE)
		{
		    for(i=0; i< *(numOfPInS + k); i++)
		    {
                       fprintf(secfp, "%d  %d\n",
		               (sectionData[k]+i)->vtX, 
		               (sectionData[k]+i)->vtY 
	                      ); 
		    }
		    if(secfp)
		    {
		       fclose(secfp);
                       secfp = NULL;
		    }
		}
		
		/*  for source contour */
                if((secfp = fopen(fullname_sectionS, "w")) == NULL )
                {
                      errNum =  WLZ_ERR_FILE_OPEN;
                }
		
                if(errNum == WLZ_ERR_NONE)
		{
		  if(inx[k] > 0)
		  {
		    for(i=0; i< *(numOfPInS + k) - 1; i++)
		    {
                       fprintf(secfp, "%d  %d  %d  %d\n",
		               (sectionData[k]+i)->vtX, 
		               (sectionData[k]+i)->vtY, 
			       (sectionData[k]+i+1)->vtX, 
		               (sectionData[k]+i+1)->vtY 
	                      ); 
		    }
		  }
		  else
		  {
		    for(i= *(numOfPInS + k)-2; i>=0; i--)
		    {
                       fprintf(secfp, "%d  %d  %d  %d\n",
		               (sectionData[k]+i)->vtX, 
		               (sectionData[k]+i)->vtY, 
			       (sectionData[k]+i+1)->vtX, 
		               (sectionData[k]+i+1)->vtY 
	                      ); 
		    }
		  }
		  
		  if(secfp)
		  {
		     fclose(secfp);
                     secfp = NULL;
		  }
		}

	       /*  for target contour	 */
                if((secfp = fopen(fullname_sectionTC, "w")) == NULL )
                {
                      errNum =  WLZ_ERR_FILE_OPEN;
                }
                if(errNum == WLZ_ERR_NONE)
		{
		    for(i=1; i<= dimV-1; i++)
		    {
                       fprintf(secfp, "%d  %d  %d  %d\n",
		              (int) fittedx[k+1][i], 
		              (int) fittedy[k+1][i], 
		              (int) fittedx[k+1][i+1], 
		              (int) fittedy[k+1][i+1] 
	                      ); 
		    }
		    if(secfp)
		    {
		       fclose(secfp);
                       secfp = NULL;
		    }
		}
                 
		/*  for tie point (old version)  */
                if((secfp = fopen(fullname_sectionT, "w")) == NULL )
                {
                      errNum =  WLZ_ERR_FILE_OPEN;
                }
                if(errNum == WLZ_ERR_NONE)
		{
		    for( i= ( *(numOfPInS + k) + 1)/3; i< *(numOfPInS + k); i += ( *(numOfPInS + k) + 1)/3 + 1 )
		    {
                       fprintf(secfp, "%d  %d\n",
		               (sectionData[k]+i)->vtX, 
		               (sectionData[k]+i)->vtY 
	                      ); 
		    }
		    if(secfp)
		    {
		       fclose(secfp);
                       secfp = NULL;
		    }
		}

               /*  for tie points with distance (new version) */
                if((secfp = fopen(fullname_sectionTN, "w")) == NULL )
                {
                      errNum =  WLZ_ERR_FILE_OPEN;
                }
                if(errNum == WLZ_ERR_NONE)
		{
		    /* for( i= 4; i<NumberToTrack ; i += NumberToTrack/3 ) */
		    {
		      /* i= 4; */
		       i= 2;
		       x1 = fittedx[k+1][(i-1)*delta+1];
		       y1 = fittedy[k+1][(i-1)*delta+1];
		       x0 = orderedSurfaceX[k+1][i];
		       y0 = orderedSurfaceY[k+1][i];
		       
                       fprintf(secfp, "%lg  %lg  %lg  %lg\n", x0, y0,  x1-x0,  y1-y0 );
		      /* 
		       i = 7;
		       X1 = fittedx[k+1][(i-1)*delta+1];
		       Y1 = fittedy[k+1][(i-1)*delta+1];
		       X0 = orderedSurfaceX[k+1][i];
		       Y0 = orderedSurfaceY[k+1][i];
		     */
		       /*  add only those points that seperated to resolve local conflicts */
		       /*
                       if(  ( abs(X0-x0) + abs(Y0-y0) > 3 ) && ( abs(Y1-y1) + abs(X1-x1) > 3 )  )
                           fprintf(secfp, "%lg  %lg  %lg  %lg\n", X0, Y0, X1-X0, Y1-Y0);
			   */
		    }
		    if(secfp)
		    {
		       fclose(secfp);
                       secfp = NULL;
		    }
		}
		

	     }

   }


      if(surfacefp)
      {
        fclose(surfacefp);
	surfacefp = NULL;
      }
      
      if(infp)
      {
        fclose(infp);
	infp = NULL;
      }
      
      if(outfp)
      {
        fclose(outfp);
	outfp = NULL;
      }


  free_matrix(u, 1,  nOfFiles, 1, polyDeg);
  free_matrix(v, 1,  polyDeg, 1, polyDeg);
  free_matrix(fittedx, 1,  nOfFiles, 1, dimV );
  free_matrix(fittedy, 1,  nOfFiles, 1, dimV );
  free_matrix(orderedSurfaceX, 1,  nOfFiles, 1, NumberToTrack  );
  free_matrix(orderedSurfaceY, 1,  nOfFiles, 1, NumberToTrack  );
      
  free_vector(w1,1,  polyDeg);
  free_vector(xp,1,  nOfFiles);
  free_vector(yp,1,  nOfFiles);
  free_vector(ac,1,  polyDeg);

  free_vector(Esy,1, nOfFiles);

  free_vector(sig,1, nOfFiles);
  free_vector(inx,1, nOfFiles);





      
   *dstErr = errNum;
      
}




/*!
* - Function:	svdfit.c ( a function from the Numeric Recipes
* - Returns:	none 
*		
* - Purpose:    modeling data by least squre fit
* - Global refs:	-
* - Parameters:	
*       -#  x[1:ndata], y[1:ndata]:  a set of data points (input) 
*       -#  sig[1...ndata]           individual standard veviations (input) default all = 1: 
*       -#  a[1:ma]                  number of coefficients y=a1+a2*x+a2*x^2+...an*x^{ma-1}
*	-#  *dstErr:	Destination error pointer, may be NULL.
* - Author:    revised by   J. Rao, R. Baldock and B. Hill
*/
static void svdfit(double x[], double y[], double sig[], int ndata, double a[],
            int ma, double **u, double **v, double w[], double *chisq, 
	    void (*funcs)(double, double [], int) )
{
        int j,i;
        double wmax,tmp,thresh,sum,*b,*afunc,*vector();
        void svdcmp(),svbksb(),free_vector();

        b=vector(1,ndata);
        afunc=vector(1,ma);
        for (i=1;i<=ndata;i++) {
                (*funcs)(x[i],afunc,ma);
                tmp=1.0/sig[i];
                for (j=1;j<=ma;j++) u[i][j]=afunc[j]*tmp;
                b[i]=y[i]*tmp;
        }
        svdcmp(u,ndata,ma,w,v);
        wmax=0.0;
        for (j=1;j<=ma;j++)
                if (w[j] > wmax) wmax=w[j];
        thresh=TOL*wmax;
        for (j=1;j<=ma;j++)
                if (w[j] < thresh) w[j]=0.0;
        svbksb(u,w,v,ndata,ma,b,a);
        *chisq=0.0;
        for (i=1;i<=ndata;i++) {
                (*funcs)(x[i],afunc,ma);
                for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
                *chisq += (tmp=(y[i]-sum)/sig[i],tmp*tmp);
        }
	/*
	for (i=1; i<ndata;i++) {
	    printf("%lg\n",afunc[i]);
	}
	*/
        free_vector(afunc,1,ma);
        free_vector(b,1,ndata);
}


static void funcs(double x, double p[], int np)  
{
/* ANSI: void (*funcs)(double,double *,int); */
  int j;
  p[1] = 1.0;
  for (j=2; j<=np; j++) p[j] = p[j-1]*x;
}


#include <math.h>

static double at,bt,ct;
#define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ? \
(ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))

static double maxarg1,maxarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static void svdcmp(double **a, int m, int n, double *w, double **v)
{
        int flag,i,its,j,jj,k,l,nm;
        double c,f,h,s,x,y,z;
        double anorm=0.0,g=0.0,scale=0.0;
        double *rv1,*vector();
        void nrerror(),free_vector();

        if (m < n) nrerror("SVDCMP: You must augment A with extra zero rows");
        rv1=vector(1,n);
        for (i=1;i<=n;i++) {
                l=i+1;
                rv1[i]=scale*g;
                g=s=scale=0.0;
                if (i <= m) {
                        for (k=i;k<=m;k++) scale += fabs(a[k][i]);
                        if (scale) {
                                for (k=i;k<=m;k++) {
                                        a[k][i] /= scale;
                                        s += a[k][i]*a[k][i];
                                }
                                f=a[i][i];
                                g = -SIGN(sqrt(s),f);
                                h=f*g-s;
                                a[i][i]=f-g;
                                if (i != n) {
                                        for (j=l;j<=n;j++) {
                                                for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
                                                f=s/h;
                                                for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
                                        }
                                }
                                for (k=i;k<=m;k++) a[k][i] *= scale;
                        }
                }
                w[i]=scale*g;
                g=s=scale=0.0;
                if (i <= m && i != n) {
                        for (k=l;k<=n;k++) scale += fabs(a[i][k]);
                        if (scale) {
                                for (k=l;k<=n;k++) {
                                        a[i][k] /= scale;
                                        s += a[i][k]*a[i][k];
                                }
                                f=a[i][l];
                                g = -SIGN(sqrt(s),f);
                                h=f*g-s;
                                a[i][l]=f-g;
                                for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
                                if (i != m) {
                                        for (j=l;j<=m;j++) {
                                                for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
                                                for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
                                        }
                                }
                                for (k=l;k<=n;k++) a[i][k] *= scale;
                        }
                }
                anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
        }
        for (i=n;i>=1;i--) {
                if (i < n) {
                        if (g) {
                                for (j=l;j<=n;j++)
                                        v[j][i]=(a[i][j]/a[i][l])/g;
                                for (j=l;j<=n;j++) {
                                        for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
                                        for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
                                }
                        }
                        for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
                }
                v[i][i]=1.0;
                g=rv1[i];
                l=i;
        }
        for (i=n;i>=1;i--) {
                l=i+1;
                g=w[i];
                if (i < n)
                        for (j=l;j<=n;j++) a[i][j]=0.0;
                if (g) {
                        g=1.0/g;
                        if (i != n) {
                                for (j=l;j<=n;j++) {
                                        for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
                                        f=(s/a[i][i])*g;
                                        for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
                                }
                        }
                        for (j=i;j<=m;j++) a[j][i] *= g;
                } else {
                        for (j=i;j<=m;j++) a[j][i]=0.0;
                }
                ++a[i][i];
        }
        for (k=n;k>=1;k--) {
                for (its=1;its<=30;its++) {
                        flag=1;
                        for (l=k;l>=1;l--) {
                                nm=l-1;
                                if (fabs(rv1[l])+anorm == anorm) {
                                        flag=0;
                                        break;
                                }
                                if (fabs(w[nm])+anorm == anorm) break;
                        }
                        if (flag) {
                                c=0.0;
                                s=1.0;
                                for (i=l;i<=k;i++) {
                                        f=s*rv1[i];
                                        if (fabs(f)+anorm != anorm) {
                                                g=w[i];
                                                h=PYTHAG(f,g);
                                                w[i]=h;
                                                h=1.0/h;
                                                c=g*h;
                                                s=(-f*h);
                                                for (j=1;j<=m;j++) {
                                                        y=a[j][nm];
                                                        z=a[j][i];
                                                        a[j][nm]=y*c+z*s;
                                                        a[j][i]=z*c-y*s;
                                                }
                                        }
                                }
                        }
                        z=w[k];
                        if (l == k) {
                                if (z < 0.0) {
                                        w[k] = -z;
                                        for (j=1;j<=n;j++) v[j][k]=(-v[j][k]);
                                }
                                break;
                        }
                        if (its == 30) nrerror("No convergence in 30 SVDCMP iterations");
                        x=w[l];
                        nm=k-1;
                        y=w[nm];
                        g=rv1[nm];
                        h=rv1[k];
                        f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
                        g=PYTHAG(f,1.0);
                        f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
                        c=s=1.0;
                        for (j=l;j<=nm;j++) {
                                i=j+1;
                                g=rv1[i];
                                y=w[i];
                                h=s*g;
                                g=c*g;
                                z=PYTHAG(f,h);
                                rv1[j]=z;
                                c=f/z;
                                s=h/z;
                                f=x*c+g*s;
                                g=g*c-x*s;
                                h=y*s;
                                y=y*c;
                                for (jj=1;jj<=n;jj++) {
                                        x=v[jj][j];
                                        z=v[jj][i];
                                        v[jj][j]=x*c+z*s;
                                        v[jj][i]=z*c-x*s;
                                }
                                z=PYTHAG(f,h);
                                w[j]=z;
                                if (z) {
                                        z=1.0/z;
                                        c=f*z;
                                        s=h*z;
                                }
                                f=(c*g)+(s*y);
                                x=(c*y)-(s*g);
                                for (jj=1;jj<=m;jj++) {
                                        y=a[jj][j];
                                        z=a[jj][i];
                                        a[jj][j]=y*c+z*s;
                                        a[jj][i]=z*c-y*s;
                                }
                        }
                        rv1[l]=0.0;
                        rv1[k]=f;
                        w[k]=x;
                }
        }
        free_vector(rv1,1,n);
}

#undef SIGN
#undef MAX
#undef PYTHAG


void free_vector(double *v, int nl, int nh)
{
        free((char*) (v+nl));
}


double *vector(nl,nh)
int nl,nh;
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl;
}

double **matrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	double **m;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

void free_matrix(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}





void nrerror(char error_text[])
{
	void exit();

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[])
{
	int jj,j,i;
	double s, *tmp, *vector();
	void free_vector();

	tmp=vector(1,n);
	for (j=1;j<=n;j++) {
		s=0.0;
		if (w[j]) {
			for (i=1;i<=m;i++) s += u[i][j]*b[i];
			s /= w[j];
		}
		tmp[j]=s;
	}
	for (j=1;j<=n;j++) {
		s=0.0;
		for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
		x[j]=s;
	}
	free_vector(tmp,1,n);
}


/*  calculate y = a[1] + a[2] x + a[3] x^2 +... + a[polyDeg x^{polyDeg-1} */

static double polynomialF( double x, double ac[], int polyDeg)
{ 
  double sum, temp;
  int i, jy;
      sum  = ac[1];
      temp = x;
      for(jy=2; jy<= polyDeg; jy++)
      {
        sum +=  ac[jy] * temp;
        temp = temp*x;
      }
   return sum;	
}


/*!
* \return       Error code.
* \ingroup      WlzTransform
* \brief        Find all vertices within all the shells of the given Geometry
*               which have high connectivity, while at the same time
*               return them in a array for later removing.
* \param        sGM                     Source model.
* \param        iBuf                    Buffer containing the indicies of
*                                       vertices to be removed from the
*                                       model.
* \param        numHV                   the number of the above vertices.
*/
static int *WlzGetAllHighConVertOfGM2D( WlzGMModel *sGM,
                                            int *numHV,  
					    WlzErrorNum *dstErr)
{
  int            nHCV, *iBuf;
  WlzGMVertex   *v0;
  WlzGMVertexT  *vT0;
  WlzGMDiskT    *dT0;
  WlzGMEdgeT    *cET,
                *fET;
  WlzGMLoopT    *cLT,
                *fLT;
  WlzErrorNum    errNum = WLZ_ERR_NONE;
  int            maxVI;
  WlzGMShell    *cS, *fS;

  /* For the model. */
  
  maxVI = sGM->res.vertex.numIdx;
  /*  alocate memory */
  if(  ( ( iBuf = (int *)AlcMalloc( maxVI * sizeof(int) ) ) == NULL)  ) 
  {
      errNum = WLZ_ERR_MEM_ALLOC;
  }
 
  /* Find all vertices with high connectivity that are in the given
   * shell. */
  if(errNum == WLZ_ERR_NONE)
  {
    cS = fS = (WlzGMShell *) sGM->child;
    nHCV = 0;
    /* For each shell of the model. */
    do
    {
      cS = cS->next;
      /* this should not happen but i have no choice, possible bug in Bill's code */
      if(!cS->child)
      {
         if(cS == fS )
	      break;
         continue; 
      }
      cLT = fLT = cS->child;
      do
      {
          cET = fET = cLT->edgeT;
          do
          {
              vT0 = cET->vertexT;
              /* Test for high connectivity */
              if((vT0 != vT0->next) && (vT0 != vT0->next->next))
              {
                 dT0 = vT0->diskT;
                 v0  = dT0->vertex;
                 /* Only delete a vertex once. */
                 if((v0->diskT == dT0) && (dT0->vertexT == vT0))
                 {
                    *(iBuf + nHCV++) = v0->idx;
                 }
              }
              cET = cET->next;
          } while(cET != fET);
          cLT = cLT->next;
      } while(cLT != fLT);

    }  while(cS != fS );
  }   
  *numHV = nHCV;
  *dstErr = errNum;
  return(iBuf);
}

/*!
* \return                               Error code.
* \ingroup      WlzTransform
* \brief        Removes all vertices in the workspace's vertex flags
*               buffer from the source model.
* \param        sGM                     Source model.
* \param        iBuf                    Buffer containing the indicies of
*                                       vertices to be removed from the
*                                       model.
* \param        nD                      Number of vertices to be removed.
*/
static WlzErrorNum WlzRemoveVertices(
                                WlzGMModel *sGM, int *iBuf, int nD)
{
  int           idD;
  AlcVector     *vec;
  WlzGMVertex   *sV;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  /* Remove all these high connectivity vertices from the model creating
   * new shells. Make sure that the transform of each of these new shells
   * is a copy of the transform of the shell with the high connectivity
   * vertex. */
  idD = 0;
  vec = sGM->res.vertex.vec;
  while((errNum == WLZ_ERR_NONE) && (idD < nD))
  {
    sV = (WlzGMVertex *)AlcVectorItemGet(vec, *(iBuf + idD));
    if(sV && (sV->idx >= 0))
    {
      errNum = WlzGMModelDeleteV(sGM, sV);
    }
    ++idD;
  }
  return(errNum);
}


/*!
* \return	Allocated vertices.
* \ingroup	WlzFeatures
* \brief	Allocates a buffer which it fills with the vertices
*		from a geometric model.
* \param	model			Given geometric model.
* \param        LoopIdx                 The loop index in the model
* \param        numOfVertices           number of vertices in this loop of the model
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
static WlzVertexP WlzVerticesThisLoopOfSimpleShell2DGM(  WlzGMModel  *model, 
					    int          LoopIdx,
					    int         *numOfVertices,
					    WlzErrorNum *dstErr   )
{
 WlzVertexP vData;
 WlzErrorNum errNum = WLZ_ERR_NONE;

 vData.v = NULL;
 if( model != NULL)
 {
    switch(model->type)
    {
      case WLZ_GMMOD_2I: /* FALLTHROUGH */
      case WLZ_GMMOD_2D:
	vData.d2 = WlzDVerticesThisLoopSimpleShellGM2(model, LoopIdx, numOfVertices, &errNum);
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
 }
 if(dstErr)
 {
    *dstErr = errNum;
 }
 return(vData);

}




/*!
* \return	Allocated vertices.
* \ingroup	WlzFeatures
* \brief	Allocates a buffer which it fills with the vertices
*		from a 2D GM.
* \param	model			Given geometric model.
* \param        LoopIdx                 The loop index in the model
* \param        numOfVertices           number of vertices in this loop of the model
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
static WlzDVertex2 *WlzDVerticesThisLoopSimpleShellGM2(WlzGMModel  *model, 
                                            int            LoopIdx, 
					    int           *numOfVertices, 
					    WlzErrorNum   *dstErr   )
{
   WlzDVertex2   *vData = NULL, *vData2 = NULL;
   WlzErrorNum   errNum = WLZ_ERR_NONE;
   int nVertices;
   AlcVector          *vec;
   WlzGMLoopT         *cLT;
   WlzGMEdgeT         *cET, *fET;
   int test = 0;
   
  if((vData = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) * (*numOfVertices) )) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
   
   vec         =  model->res.loopT.vec;
   nVertices   = 0;
    /* get the loopT of the model. */
    cLT = (WlzGMLoopT *) AlcVectorItemGet(vec, LoopIdx);
    /* For each edge topology element of the model. */    
    cET = fET = cLT->edgeT;
    if(test)
        printf("test loopIdx %d\n ",LoopIdx );

    if(IsACycle(  model, LoopIdx, &errNum) )
    {
       do
       {
           if(nVertices >= (*numOfVertices) )
	     break;

	   if(model->type == WLZ_GMMOD_2I)
	   {
	      (vData + nVertices)->vtX     = (double) cET->vertexT->diskT->vertex->geo.vg2I->vtx.vtX;
	      (vData + nVertices)->vtY     = (double) cET->vertexT->diskT->vertex->geo.vg2I->vtx.vtY;
              nVertices += 1;
	   }
           else
           {

	      (vData + nVertices)->vtX     = cET->vertexT->diskT->vertex->geo.vg2D->vtx.vtX;
	      (vData + nVertices)->vtY     = cET->vertexT->diskT->vertex->geo.vg2D->vtx.vtY;
	      nVertices += 1;

	   }
	   /* if( ( cET->vertexT->diskT->idx == cET->next->next->vertexT->diskT->idx ) ) */
	   /*       break; */
	   cET = cET->next;
       } while (cET != fET  );
   }
   else
   {
      /*  not cycle but a simple shell, should have two ends got it first: */
      /* For each loop topology element of the model. go to the one of the end points*/
      do
      {
	 if( ( cET->prev->vertexT->diskT->idx == cET->next->vertexT->diskT->idx ) )
	 {
	    break;
	 }    
	 cET = cET->next;
      } while (cET != fET);
       

      /*  get the vertices until we meet the other ends ! */
      do
      {
 	 (vData + nVertices)->vtX     = cET->vertexT->diskT->vertex->geo.vg2D->vtx.vtX;
	 (vData + nVertices)->vtY     = cET->vertexT->diskT->vertex->geo.vg2D->vtx.vtY;
         nVertices += 1;
	 if(nVertices >= (*numOfVertices)  )
	 {
            break;
	 }

	 if( ( cET->prev->vertexT->diskT->idx == cET->next->vertexT->diskT->idx ) && ( nVertices > 1) )
	 {
	    /* this part should not be happen but as bill's code to delete vertices , I have to 
	       use this part to get over it should be corrected in the future ! */
	    if( nVertices != (*numOfVertices) )
	    {
                 if((vData2 = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) * nVertices )) == NULL)
                 {
                    errNum = WLZ_ERR_MEM_ALLOC;
                 }
                 if(errNum == WLZ_ERR_NONE)
                 {
                       WlzValueCopyDVertexToDVertex( vData2,  vData, nVertices);
		       AlcFree(vData);
		       if(dstErr)
		       {
		          *dstErr = errNum;
		       }
		       *numOfVertices = nVertices;
		       return(vData2);
		 }
	    }
	    break;
	 }   
	 cET = cET->next;
      } while (cET != fET);
         
   }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(vData);
}


static void testBeforeBreakDownTheGM(  WlzGMModel *gModel[], 
                        int numOf2DWlzFiles,
			int sectionLength_N,
                        WlzErrorNum   *dstErr )
{

   int   loopForFile, loopForShell, loopForLoop;
   int   numOfShellsInStartGM, *startShellIdx, numOfLoopsInStartShell, *startLoopIdx, numOfVerticesInTheStartLoop;
   WlzVertexP   AllVerticesInThisStartLoop;
   WlzErrorNum   errNum = WLZ_ERR_NONE;

 
   /* out put for test */
   printf("BEFORE break down:::\n");
   printf("teststart\n");
   
     for(loopForFile=0; loopForFile<numOf2DWlzFiles ; loopForFile++)
     {
       printf("------ file %d test: -------\n", loopForFile); 
        /* get the number of Shells in this model */
        numOfShellsInStartGM  = NumberOfShellsAboutGM(gModel[loopForFile]);
	startShellIdx         = WlzShellsIndexAboutGM(gModel[loopForFile], numOfShellsInStartGM, &errNum);
       
       for(loopForShell=0; loopForShell< numOfShellsInStartGM; loopForShell++)
       {
	 numOfLoopsInStartShell  =  NumberOfLoopsInThe_nThShellOfGM(gModel[loopForFile],  startShellIdx[loopForShell] );
         /* get the index for the loops */
         startLoopIdx            =  WlzLoopIndexAboutGM(gModel[loopForFile], startShellIdx[loopForShell], 
	                                      numOfLoopsInStartShell, &errNum);
           printf("The  %d th Shell\n",  loopForShell); 
           printf("num of loops is    %d\n", numOfLoopsInStartShell ); 
	 for(loopForLoop=0; loopForLoop<numOfLoopsInStartShell ; loopForLoop++)
	 {
	   /* get number of vertices in this loop */
           numOfVerticesInTheStartLoop = NumberOfVerticesInTheLoopOfGM(
	 				gModel[loopForFile], 
	                                startLoopIdx[loopForLoop], &errNum);
					
	    printf("loop %d with %d number of vertices\n", loopForLoop, numOfVerticesInTheStartLoop );
	   if( numOfVerticesInTheStartLoop < sectionLength_N )
	                  break;
	   
           AllVerticesInThisStartLoop  = WlzVerticesThisLoopOfGM( gModel[loopForFile], 
	                                                          startLoopIdx[loopForLoop], 
								  numOfVerticesInTheStartLoop, 
								 &errNum  );
	   /*							    
           for(k=0; k<numOfVerticesInTheStartLoop; k++)
	   {
                printf("%lg  %lg  %lg\n", ( AllVerticesInThisStartLoop.d2 + k )->vtX,
		                          ( AllVerticesInThisStartLoop.d2 + k )->vtY, 0.0 );
	   }
	   */

	   if(AllVerticesInThisStartLoop.d2)
	       AlcFree(AllVerticesInThisStartLoop.d2);
	 }				
         if(startLoopIdx)
	   AlcFree(startLoopIdx);

       }
       
      if(startShellIdx)
       AlcFree(startShellIdx);
     } 
}


static void testJustAfterBreakDownTheGM(  WlzGMModel *gModel[], 
                        int numOf2DWlzFiles,
			int sectionLength_N,
                        WlzErrorNum   *dstErr )
{

   int   k, loopForFile, loopForShell, loopForLoop;
   int   numOfShellsInStartGM, *startShellIdx, numOfLoopsInStartShell, *startLoopIdx, numOfVerticesInTheStartLoop;
   WlzVertexP   AllVerticesInThisStartLoop;
   WlzErrorNum   errNum = WLZ_ERR_NONE;


   /* out put for test */
   printf("AFTER BREAK INTO SIMPLE GM\n");

     for(loopForFile=0; loopForFile< numOf2DWlzFiles; loopForFile++)
     {
         printf("=======file  %d test =======\n", loopForFile); 
        /* get the number of Shells in this model */
        numOfShellsInStartGM  = NumberOfShellsAboutGM(gModel[loopForFile]);
	startShellIdx         = WlzShellsIndexAboutGM(gModel[loopForFile], numOfShellsInStartGM, &errNum);

           printf("The number of Shells is  %d\n", numOfShellsInStartGM ); 


       for(loopForShell=0; loopForShell< numOfShellsInStartGM; loopForShell++)
       {
	 numOfLoopsInStartShell  =  NumberOfLoopsInThe_nThShellOfGM(gModel[loopForFile],  startShellIdx[loopForShell] );
         /* get the index for the loops */
         startLoopIdx            =  WlzLoopIndexAboutGM(gModel[loopForFile], startShellIdx[loopForShell], 
	                                      numOfLoopsInStartShell, &errNum);
           printf("The  %d th Shell\n",  loopForShell); 
          printf("num of loops is    %d\n", numOfLoopsInStartShell ); 
	 for(loopForLoop=0; loopForLoop<numOfLoopsInStartShell ; loopForLoop++)
	 {
	   if(errNum != WLZ_ERR_NONE)
	   {
              if( startLoopIdx  )
	         AlcFree(startLoopIdx);
	      break;	 
	   }
	   /* get number of vertices in this loop */
           numOfVerticesInTheStartLoop = NumberOfVerticesInTheLoopOfGM(
	 				gModel[loopForFile], 
	                                startLoopIdx[loopForLoop], &errNum);
					
	   printf("loop %d with %d number of vertices\n", loopForLoop, numOfVerticesInTheStartLoop ); 
	   /* if( numOfVerticesInTheStartLoop < 2 ) */
	  if( numOfVerticesInTheStartLoop < sectionLength_N ) 
                  break;
	   
           AllVerticesInThisStartLoop  = WlzVerticesThisLoopOfGM( gModel[loopForFile], 
	                                                             startLoopIdx[loopForLoop], 
								     numOfVerticesInTheStartLoop, 
								    &errNum  );
	   /* if( (loopForFile == 8) && (loopForShell == 62)) */
	   {
              for(k=0; k<numOfVerticesInTheStartLoop; k++)
	      {
                printf("%lg  %lg  %lg\n", ( AllVerticesInThisStartLoop.d2 + k )->vtX,
		                          ( AllVerticesInThisStartLoop.d2 + k )->vtY, (double) ( loopForFile * 5 ) );
	      }
	   } 
	   if(AllVerticesInThisStartLoop.d2)
	       AlcFree(AllVerticesInThisStartLoop.d2);
	 }				
         if(startLoopIdx)
	   AlcFree(startLoopIdx);

       }
       
      if(startShellIdx)
       AlcFree(startShellIdx);
       } 


}


static void testConsistencyBreakDownTheGM(  WlzGMModel *gModel[], 
                        int numOf2DWlzFiles,
			int sectionLength_N,
                        WlzErrorNum   *dstErr )
{

   int   k, endLoop, loopForFile, loopForShell, loopForLoop;
   int   numOfShellsInStartGM, *startShellIdx, numOfLoopsInStartShell, *startLoopIdx, numOfVerticesInTheStartLoop;
   WlzVertexP   AllVerticesInThisStartLoop;
   WlzErrorNum   errNum = WLZ_ERR_NONE;

   /* out put for test */
   printf("Consitent test with the tracking of SIMPLE GM\n");

     for(loopForFile=0; loopForFile< numOf2DWlzFiles; loopForFile++)
     {
        /* printf("=======file  %d test =======\n", loopForFile); */
        /* get the number of Shells in this model */
        numOfShellsInStartGM  = NumberOfShellsAboutGM( gModel[loopForFile] );
	startShellIdx         = WlzShellsIndexAboutGM( gModel[loopForFile], numOfShellsInStartGM, &errNum );
	/* printf("The number of Shells is: %d\n", numOfShellsInStartGM ); */
	if(errNum == WLZ_ERR_NONE)
	{

          for(loopForShell=0; loopForShell< numOfShellsInStartGM; loopForShell++)
          {
            /* get number of loops in this shell */
	    numOfLoopsInStartShell  =  NumberOfLoopsInThe_nThShellOfGM( gModel[loopForFile], startShellIdx[loopForShell] );
            /* get the index for the loops */
            startLoopIdx            =  WlzLoopIndexAboutGM(gModel[loopForFile], startShellIdx[loopForShell], 
	                                      numOfLoopsInStartShell, &errNum);
            /* printf("The  %d th Shell\n",  loopForShell); 
            printf("num of loops is    %d\n", numOfLoopsInStartShell ); */
	    endLoop = numOfLoopsInStartShell;
	    for(loopForLoop=0; loopForLoop<endLoop ; loopForLoop++)
	    {
	      if(errNum != WLZ_ERR_NONE)
	      {
                 if( startLoopIdx  )
	            AlcFree(startLoopIdx);
	         break;	 
	      }

              if( errNum == WLZ_ERR_NONE)
              {
	        /* get number of vertices in this loop */
                numOfVerticesInTheStartLoop =   NumberOfVerticesInTheLoopOfGM(gModel[loopForFile], startLoopIdx[loopForLoop], &errNum);

	        if( ( numOfLoopsInStartShell == 2 ) && (loopForLoop == 0 ) )
	        {
		    if( IsACycle(gModel[loopForFile],  startLoopIdx[loopForLoop], &errNum ) )
		    {
	                     endLoop = 1;
		    }	     
		    else
		    {
                             numOfVerticesInTheStartLoop =  numOfVerticesInTheStartLoop/2 + 1;
		    }
	        }	 
	        else
	        {
                    numOfVerticesInTheStartLoop =  numOfVerticesInTheStartLoop/2 + 1;
	        }
		if(numOfVerticesInTheStartLoop < sectionLength_N )
		      break;
              }

              /*
	      printf("loop %d with %d number of vertices\n", loopForLoop, numOfVerticesInTheStartLoop ); 
	      */
	      /* if( numOfVerticesInTheStartLoop < 2 ) */
	      if( numOfVerticesInTheStartLoop < sectionLength_N ) 
                  break;
	      /*
              AllVerticesInThisStartLoop  = WlzVerticesThisLoopOfGM( gModel[loopForFile], 
	                                                             startLoopIdx[loopForLoop], 
								     numOfVerticesInTheStartLoop, 
								    &errNum  );
								    */
	      AllVerticesInThisStartLoop  = WlzVerticesThisLoopOfSimpleShell2DGM( gModel[loopForFile], 
	                                                             startLoopIdx[loopForLoop], 
								    &numOfVerticesInTheStartLoop, 
								    &errNum  );
							    
	      /* if( (loopForFile == 8) && (loopForShell == 62)) */
	      {
                 for(k=0; k<numOfVerticesInTheStartLoop; k++)
	         {
                    printf("%lg  %lg  %lg\n", ( AllVerticesInThisStartLoop.d2 + k )->vtX,
		                          ( AllVerticesInThisStartLoop.d2 + k )->vtY, (double) ( loopForFile * 5 ) );
	         }
	      } 
	      if(AllVerticesInThisStartLoop.d2)
	          AlcFree(AllVerticesInThisStartLoop.d2);
	   }				
           if(startLoopIdx)
	      AlcFree(startLoopIdx);

         }
        } 
        if(startShellIdx)
        AlcFree(startShellIdx);
      } 
     
}
