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
#include <float.h>
#include <Wlz.h>

#define	IN_RECORD_MAX   (1024)
#define	NumberToTrack   (10)
static double SixTimesOfVoulumeOfTetraHedron(WlzDVertex3 p0, 
			WlzDVertex3 p1,
			WlzDVertex3 p2,
			WlzDVertex3 p3
);

static int IsDVertexZero(WlzDVertex3  normal);
static int IsDVertexEqual(WlzDVertex3 p1, WlzDVertex3 p2);

static void RecalculateInAndOutPointsByUsing3DNormal(
				            double         distance,
		                            WlzDVertex3   *SurfacePoints, 
					    int            nOfVInTheSFile,
					    WlzDVertex3   *InPoints,      
					    int           *nOfVInTheInFile,
					    WlzDVertex3   *OutPoints,     
					    int           *nOfVInTheOutFile );

static void pureThisGM(WlzGMModel *gMC, int LoopIdx, int numOfEnds, 
                       int *EndsIdx, int *bifurIdx, 
		       WlzErrorNum *dstErrgMC);

static void outputVerticesInThisLoop(WlzGMModel *gMC, int LoopIdx, WlzErrorNum *dstErr);

static void WlzGMModelPure( WlzGMModel *gMC, WlzErrorNum *dstErr);

static void GetInOutAndFromSurfacePoints( 
				   int          numberOfPixelsZ,
                                   WlzDVertex2  StandSampleP[NumberToTrack],
                                   double       distance,
				   WlzDVertex3 *InPoints,
				   int         *nOfVInTheInFile,
				   WlzDVertex3 *OutPoints,
				   int         *nOfVInTheOutFile,
				   WlzErrorNum *dstErr );

static void	outputSurfaceInAndOutPoints(FILE          *surfacefp, 
					    FILE          *infp,  
					    FILE          *outfp, 
		                            WlzDVertex3   *SurfacePoints, 
					    int            nOfVInTheSFile,
					    WlzDVertex3   *InPoints,      
					    int            nOfVInTheInFile,
					    WlzDVertex3   *OutPoints,     
					    int            nOfVInTheOutFile);

static void GetSamplePointsFromAllVertices(   WlzVertexP   sLoopVData, 
                                              int          i_s,
					      int          subSubSectionLength_L,       
                                              WlzDVertex2  StandSampleP[NumberToTrack],
                                              WlzErrorNum *dstErr );

static int NumberOfVerticesDoubleInNonCycleCaseInTheLoopOfGM(WlzGMModel   *gM, 
 						int sLoopIdx, 
						WlzErrorNum *dstErr);

static int *GetThebifurcatePointsIndexInTheNonCycleLine(   WlzGMModel  *gM, 
                                           int          LoopIdx,
					   int          numOfEnds,
                                           WlzErrorNum *dstErr );


static int *GetTheEndsPointsIndexInTheNonCycleLine(   WlzGMModel  *gM, 
                                           int          LoopIdx,
					   int          numOfEnds,
                                           WlzErrorNum *dstErr );

static int HowManyEndsInTheNonCycleLine(   WlzGMModel  *gM, 
                                           int          LoopIdx, 
                                           WlzErrorNum *dstErr );

static int IsACycle(   WlzGMModel  *gM, int  LoopIdx, WlzErrorNum *dstErr );

static void  WlzOutputForFormingSurface(FILE *testFile, char *Fstr, WlzVertexP sLoopVData, int nV, double zC, int sampleD, WlzErrorNum *dstErr);

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
				   int         *nOfVInTheInFile,
				   WlzDVertex3 *OutPoints,
				   int         *nOfVInTheOutFile,
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
WlzDVertex3  *WlzGeometryTrackUpAndDown_s( WlzObject          *sObj,   
                                           WlzObject          *tObj,
				           int                 numberOfPixelsZ,
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
		                           WlzErrorNum        *dstErr
		                     )
{
  WlzErrorNum	      errNum = WLZ_ERR_NONE;
  WlzGMModel         *gM, *gMS, *gMT, *gMC;
  WlzGMModel         *gM2;


  int                *sNN, suc;
  AlcKDTTree         *tTree;
  AlcKDTTree         *tTreeArray[numOf2DWlzFiles];
  AlcKDTNode         *node;
  double             datD[3];
   
  WlzGMShell         *cS, *fS;
  WlzGMLoopT         *cLT, *fLT;
  WlzGMEdgeT         *cET, *fET;
  WlzDVertex2        *v1,  *v2, *vNorm=NULL;
  WlzDVertex2	      segV[3];
  WlzDVertex2         StandSampleP[NumberToTrack];
  WlzDVertex2         TrackedSampleP[NumberToTrack], TrackedSamplePBlgToOneLp[NumberToTrack];
  WlzDVertex2         TrackedSampleFirstP, TrackedSampleSecondP;
  WlzIVertex2         wV2, wV2I;
  WlzDVertex2         wD1, wD2;
  WlzDVertex3        *surfP, *SurfacePoints, *InPoints, *OutPoints;
  WlzDVertex2         fPoint, *inPoints, *outPoints;
  WlzVertexP          nData[2], vData[2];
  WlzVertexP          sLoopVData, tLoopVData, AllVerticesInThisStartLoop;
  WlzVertexType       vType;
   
  WlzObject           *newSObj, *newTObj;
   
  int                 indexStr[NumberToTrack], LoopIdx;
  int                 LoopIndex[NumberToTrack];
  int                 ShellIndex[NumberToTrack];
  int                 *GroupNum, *FixNumForEachGroup;
  int                 NumOfGroup;
  int                 sp, nbgN, trackedBigestLoopIndex;
  int                 loopDirection, ncycals;
  int                 numOfLoopsInStartShell, numOfShellsInStartGM;
  int                 numOfVerticesInTheStartLoop, numOfSections, i_s;
  int                 nOfVInTheSFile, nOfVInTheInFile, nOfVInTheOutFile;

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

  int                 i,j, k, ix, jy, it, jt, kt, test = 0, nE;
  int                 tempI, tempJ, tempK, tempL, tempM, tempN;
  int                 nS, nL, nV, ntS, ntL, ntV;
  int                 Continue;
  double              TempD, Tempd;
  double              di, dj, dk;
  double              Sumfr_halpha, Ialpha;
  double              HA, HTB, HATB, YATB, MI;
  double              wg[4];
  double              h,w, zC;
  double              minDis = 10., distance, zdvoxels;
  int                 vCnt[2], idN;
  int                 deltaX, deltaY, startI, fianlI;
  FILE                *testFile = NULL, *surfacefp, *infp= NULL, *outfp=NULL;
  char                *testFileStr = "outputForTest.dat";
 
  /*
  for(i=0; i<numOf2DWlzFiles ; i++)
  {
	printf("%s\n", TwoDImageFilesNameList[i]);
  }

  exit(0);
  */
  distance = 10.;

  /* initial the outpoint  */
  
  if( ( surfP = (WlzDVertex3 *)AlcMalloc(sizeof(WlzDVertex3) * 100 )) == NULL ) 
  {
     errNum = WLZ_ERR_MEM_ALLOC;
  }
  if((infp = fopen(surfaceInPointFileName, "w")) == NULL )
  {
     errNum =  WLZ_ERR_FILE_OPEN;
  }
  if((outfp = fopen(surfaceOutPointFileName, "w")) == NULL )
  {
     errNum =  WLZ_ERR_FILE_OPEN;
  }
  if((surfacefp = fopen(surfacePointFileName, "w")) == NULL )
  {
     errNum =  WLZ_ERR_FILE_OPEN;
  }

  // initial the standard obj first:
  if( errNum == WLZ_ERR_NONE)
  {
      //i = numOf2DWlzFiles/2;
      i = 0;
      // read the first one as a start point
      if((testFile = fopen( (const char *)  TwoDImageFilesNameList[i], "r")) == NULL )
      {
        printf("cannot open the standard contour Wlz file .\n");
        exit(1);
      }

      if( !(newSObj = WlzReadObj(testFile, &errNum) ) )
      {
         printf("input Woolz Object Error.\n");
         fclose(testFile); 
         exit(1);
      }
      fclose(testFile);
      
      testFile = NULL;
      if( ( (newSObj->type) != WLZ_CONTOUR )  )
      {
       	 errNum = WLZ_ERR_DOMAIN_TYPE;
      }

      if(  ( newSObj->domain.core == NULL  ) )
      {
	 errNum = WLZ_ERR_DOMAIN_NULL;
      }

      if( errNum == WLZ_ERR_NONE)
      {
        // get the G Model for standard
        gM       = newSObj->domain.ctr->model;

	// copy the model
	//  gMC      = WlzGMModelCopy( gM, &errNum);

      }

      // pure the model
      if( errNum == WLZ_ERR_NONE)
      {
      
	// WlzGMModelPure(gMC, &errNum);
	//  exit(0);

        // get the number of Shells in this model
        numOfShellsInStartGM  = NumberOfShellsAboutGM(gM);
        // get the Index for each Shell;
        startShellIdx = WlzShellsIndexAboutGM(gM, numOfShellsInStartGM, &errNum);
      }

       // get the stand points from on loop of the GM
       /*
          LoopIdx = -1;
	  
          GetSamplePointsFromOneLoopOfGM( gM, 
                                    LoopIdx,
                                    StandSampleP,
                                    &errNum );
	*/			    
      zdvoxels = 0;			
   } 


   // allocate memory to store the tracked surface , in and out points.
   if( ( SurfacePoints = ( WlzDVertex3 *)AlcMalloc(sizeof( WlzDVertex3) * ( numOf2DWlzFiles * NumberToTrack  ) )) == NULL ) 
   {
      errNum = WLZ_ERR_MEM_ALLOC;
   }
   
   if( ( InPoints      = ( WlzDVertex3 *)AlcMalloc(sizeof( WlzDVertex3) * (3 * numOf2DWlzFiles ) )) == NULL ) 
   {
      errNum = WLZ_ERR_MEM_ALLOC;
   } 
   
   if( (OutPoints      = ( WlzDVertex3 *)AlcMalloc(sizeof( WlzDVertex3) * (3 * numOf2DWlzFiles) )) == NULL ) 
   {
      errNum = WLZ_ERR_MEM_ALLOC;
   }

   // as kd tree is very expensive to build so build once and use them afterwards 
   if( errNum == WLZ_ERR_NONE)
   {
        // begin the track process:
        for(i=0; i<numOf2DWlzFiles; i++)
        {
                if((testFile = fopen( (const char *) TwoDImageFilesNameList[i], "r")) == NULL )
                {
         		printf("cannot open the standard contour Wlz file .\n");
         		exit(1);
      		}

      		if( !(newTObj = WlzReadObj(testFile, &errNum) ) )
      		{
         		printf("input Woolz Object Error.\n");
         		fclose(testFile); 
         		exit(1);
      		}
      		fclose(testFile); 
      		testFile = NULL;

   
      		if(  ( (newTObj->type) != WLZ_CONTOUR ) )
      		{
       	 		errNum = WLZ_ERR_DOMAIN_TYPE;
      		}

      		if(  ( newTObj->domain.core == NULL ) )
      		{
	 		errNum = WLZ_ERR_DOMAIN_NULL;
      		}

      		if( errNum == WLZ_ERR_NONE)
      		{
	           // get a KD Tree for the target vertices of the Woolz contour
       		   tTreeArray[i] = GetTreeFromGM(newTObj, &errNum);
                }
		
      		if( errNum == WLZ_ERR_NONE)
      		{
	           // free the obj 
       		   AlcFree(newTObj);
                }
	 }

   }	 


   // cycle for each shell

   if( errNum == WLZ_ERR_NONE)
   {
       printf("number of Shells = %d\n",  numOfShellsInStartGM );
       if(startShell <0)
            startShell = 0;
       if(endShell < 0)
            endShell   = numOfShellsInStartGM;
       for(j=startShell; j< endShell; j++)
        {
           // get number loops in this shell
           numOfLoopsInStartShell  =  NumberOfLoopsInThe_nThShellOfGM(gM,  startShellIdx[j] );
	    printf("The  %d th Shell\n",  j );


           // get the index for the loops
           startLoopIdx            =  WlzLoopIndexAboutGM(gM, startShellIdx[j], numOfLoopsInStartShell, &errNum);

	   // output for test:
	   test = 0;
	   if(test)
	      outputVerticesInThisLoop( gM, startShellIdx[j], &errNum );
	   test = 0;	
	   // always use the first loop (as shell has as most as two loops and this only happens where
	   // the shell is a cycle !)
        

           if( errNum == WLZ_ERR_NONE)
           {
	      // get number of vertices in this loop
              //numOfVerticesInTheStartLoop  = NumberOfVerticesDoubleInNonCycleCaseInTheLoopOfGM(
	      //					gM, startLoopIdx[0], &errNum);
              numOfVerticesInTheStartLoop  = NumberOfVerticesInTheLoopOfGM(
						gM, startLoopIdx[0], &errNum);
           }
      
           if( errNum == WLZ_ERR_NONE)
           {
              //get all the vertices in this loop (accurately it double the vertices!!)
              AllVerticesInThisStartLoop  = WlzVerticesThisLoopOfGM( gM, 
	                                                            startLoopIdx[0], 
								    numOfVerticesInTheStartLoop, 
								    &errNum  );
           }

           if( errNum == WLZ_ERR_NONE)
           {
               // sampe into serveral sections
	       numOfSections    =   numOfVerticesInTheStartLoop/sectionLength_N;
	       printf("There are %d number of vertices in This Shell\n",   numOfVerticesInTheStartLoop);
	       
	       printf("And it is divided into %d sections\n",  numOfSections );
	      


	       // cycle through the sections
	       if( startSection < 0 )
	            startSection = 0;
	       if( endSection < 0 )
	            endSection = numOfSections;
	       for(i_s =startSection; i_s <endSection; i_s++)
	       {
	          printf("The  %d th section\n",  i_s );

	          nOfVInTheSFile   = 0;
	          nOfVInTheInFile  = 0;
	          nOfVInTheOutFile = 0;

                  // get the sampe points in the i_s th section
                  GetSamplePointsFromAllVertices( AllVerticesInThisStartLoop, 
                                              i_s * sectionLength_N,
                                              subSubSectionLength_L,
                                              StandSampleP,
                                             &errNum );
                  if( errNum == WLZ_ERR_NONE)
                  {
			zdvoxels = 0;
			// cycle through the 2D image files (tracking process)
			suc = 1;
         		for(i=0; i<numOf2DWlzFiles; i++)
         		{      
                		if((testFile = fopen( (const char *) TwoDImageFilesNameList[i], "r")) == NULL )
                		{
         				printf("cannot open the standard contour Wlz file .\n");
         				exit(1);
      				}

      				if( !(newTObj = WlzReadObj(testFile, &errNum) ) )
      				{
         				printf("input Woolz Object Error.\n");
         				fclose(testFile); 
         				exit(1);
      				}
      				fclose(testFile); 
      				testFile = NULL;

   
      				if(  ( (newTObj->type) != WLZ_CONTOUR ) )
      				{
       	 				errNum = WLZ_ERR_DOMAIN_TYPE;
      				}

      				if(  ( newTObj->domain.core == NULL ) )
      				{
	 				errNum = WLZ_ERR_DOMAIN_NULL;
      				}
				
	       			if( errNum == WLZ_ERR_NONE   )
       				{
 
         	        		printf("------Target -------\n");
          				// get the G Model for target to be tracked
         		 		gM2 = newTObj->domain.ctr->model;

	  				// get the tracked points in a loop and get the loop index and
	 		 		// number of sample points
	 		 		GetTrackedSamplePointsFromOneLoopOfGM(gM2, 
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
				        if(newTObj)
	         			    WlzFreeObj(newTObj);
                 			break; // break this section track

	      			}

              			if(errNum == WLZ_ERR_NONE)
	      			{

	            			// output is in the StandSampleP
	     	    			GetTrackedPlusSomeSamplePointsFromOneLoopOfGM( gM2,
                                                       ntL, 
                                                       StandSampleP,
				                       TrackedSamplePBlgToOneLp,
				                       nOfTracked,
						      &suc,
                                                      &errNum );
				}		      

	              		if(errNum == WLZ_ERR_NONE)
	      			{

				        if(!suc)
	      				{
	        	 			//AlcFree(tTree);
						if(newTObj)
	         				   WlzFreeObj(newTObj);
                 				   break; // break this section track

	      				}
	      

           	    			// store the output sample points:
  	       	    			for(k=0; k<NumberToTrack; k++)
	            			{	
                   	  			SurfacePoints[nOfVInTheSFile].vtX = StandSampleP[k].vtX;
                   	  			SurfacePoints[nOfVInTheSFile].vtY = StandSampleP[k].vtY;
                   	  			SurfacePoints[nOfVInTheSFile].vtZ = zdvoxels;
                                		nOfVInTheSFile++;
	       	    			}

	   	    			// get the loop index we have just tracked
           	    		         // calculate the in and out surface points 

					/*
					  bugs in this part:
					  
           	    			InAndOutSurfacePoints(gM2, ntL, (int) zdvoxels, StandSampleP, 
		                        	 distance, 
                                        	 InPoints,   &nOfVInTheInFile, 
					 	 OutPoints,  &nOfVInTheOutFile,
						&suc,
						 &errNum);
					*/	 

                                        GetInOutAndFromSurfacePoints( (int) zdvoxels, 
					                                StandSampleP,
								            distance,
                                        	                            InPoints,   
								     &nOfVInTheInFile, 
					 	                           OutPoints,  
								     &nOfVInTheOutFile,
								     &errNum);


	      			  }
              			  // prepare for next section track:
				  if(!suc)
				  {
				                if(newTObj)
	         				   WlzFreeObj(newTObj);
                 				break; // break this section track
				  }
				  if(newTObj)
	      			     WlzFreeObj(newTObj);
                           	  zdvoxels += (double) numberOfPixelsZ;


           		   }   //end of track the 2D image files array
       		 }


                if(suc)
	        {
                     // recalculate in and out points by using 3D 
		     // (as they are needed in Bill's code which seems to
		     //  use the in and outside points to get normal for
		     //  the interpolate surface!!!)

		     //Calculate the IN and Out points by using 3D normal
		     RecalculateInAndOutPointsByUsing3DNormal(  distance,
		     						SurfacePoints, nOfVInTheSFile,
		     						InPoints,      &nOfVInTheInFile, 
		     						OutPoints,     &nOfVInTheOutFile
									);
		
        	     //output points

		     outputSurfaceInAndOutPoints(surfacefp, infp,  outfp, 
		                                 SurfacePoints, nOfVInTheSFile,
						 InPoints,      nOfVInTheInFile,
						 OutPoints,     nOfVInTheOutFile);

	        }	 
             }// end of track each section

       } 


     } 
     //prepare for the next shell tracking
     AlcFree(startLoopIdx);
     AlcFree(AllVerticesInThisStartLoop.d2);


  }// end of cycle for the shell

 
  AlcFree(SurfacePoints);
  AlcFree(InPoints);
  AlcFree(OutPoints);


  for(i=0; i<numOf2DWlzFiles; i++)
  {
     AlcFree(tTreeArray[i]);
  }
  AlcFree(startShellIdx);

  // track down first
  if( errNum == WLZ_ERR_NONE)
  {
    /*
    for(j=i; j<numOf2DWlzFiles; j++)
    {
      // read the target one
      if((testFile = fopen(TwoDImageFilesNameList[j], "r")) == NULL )
      {
         printf("cannot open the standard contour Wlz file .\n");
         exit(1);
      }

      if( !(newTObj = WlzReadObj(testFile, &errNum) ) )
      {
         printf("input Woolz Object Error.\n");
         fclose(testFile); 
         exit(1);
      }
      fclose(testFile); 
      testFile = NULL;

   
      if(  ( (newTObj->type) != WLZ_CONTOUR ) )
      {
       	 errNum = WLZ_ERR_DOMAIN_TYPE;
      }

      if(  ( newTObj->domain.core == NULL ) )
      {
	 errNum = WLZ_ERR_DOMAIN_NULL;
      }




      if( errNum == WLZ_ERR_NONE)
      {


       // get a KD Tree for the target vertices of the Woolz contour
       tTree = GetTreeFromGM(newTObj, &errNum);

  
       if( errNum == WLZ_ERR_NONE   )
       {
 
          printf("------Target -------\n");
          // get the G Model for target to be tracked
          gM2 = newTObj->domain.ctr->model;

	  // get the tracked points in a loop and get the loop index and
	  // number of sample points
	  GetTrackedSamplePointsFromOneLoopOfGM(gM2, 
	                                        tTree, 
					       &ntL, 
						numberOfPixelsZ,
	                                        minDis,
	                                        StandSampleP,
						TrackedSamplePBlgToOneLp,
						&nOfTracked,
						&suc,
						&errNum);
	


        if(errNum == WLZ_ERR_NONE)
        {
             ncycals  = IsACycle(gM2, ntL, &errNum);
	     if(errNum == WLZ_ERR_NONE)
	     {
	        printf("cycle = %d\n", ncycals );
	        if(!ncycals)
	        {   
		    numOfEnds = HowManyEndsInTheNonCycleLine( gM2, 
                                                              ntL, 
                                                             &errNum ); 
                    printf("Number of Ends = %d\n", numOfEnds );
	        }			   
	        if(errNum == WLZ_ERR_NONE && (!ncycals) )
	        {
                     EndPointsIdx =  GetTheEndsPointsIndexInTheNonCycleLine( gM2, 
                                                                         ntL,
					                                 numOfEnds,
                                                                        &errNum );
		     for(nE=0; nE<numOfEnds; nE++)
		     {
                          printf("End P idx = %d\n", *(EndPointsIdx + nE) );

		     }
		     AlcFree(EndPointsIdx);
		     if( (errNum == WLZ_ERR_NONE )&& ( numOfEnds > 2) )
		     {

                         BifurcationPointsIdx = GetThebifurcatePointsIndexInTheNonCycleLine( gM2, 
                                                                       ntL,
			                                               numOfEnds,
                                                                       &errNum );
			 for(nE=0; nE<numOfEnds-2; nE++)
			 {
			    printf("Bifurcate P idx = %d\n", *(BifurcationPointsIdx + nE) );

			 }
			 AlcFree(BifurcationPointsIdx);
			 
		     }
		}


	     }
	   // output is in the StandSampleP
	     GetTrackedPlusSomeSamplePointsFromOneLoopOfGM( gM2,
                                                            ntL, 
                                                            StandSampleP,
				                            TrackedSamplePBlgToOneLp,
				   	                    nOfTracked,
                                                            &errNum );

           // output sample points from the tracked loop:
  	       for(k=0; k<NumberToTrack; k++)
	       {
                   fprintf(surfacefp, "%6.2f  %6.2f   %6.2f\n",
		   StandSampleP[k].vtX, 
	           StandSampleP[k].vtY, 
		   (double) zdvoxels  );
	       }
	   
	   // get the loop index we have just tracked
           distance = 10.;
	   // calculate the in and out surface points 
           InAndOutSurfacePoints(gM2, ntL, zdvoxels, StandSampleP, distance,infp, outfp, &errNum);

	}
        // prepare for next cycle:
	AlcFree(tTree);
	AlcFree(newTObj);
	//newTObj = NULL;
	AlcFree(gM2);
	//gM2 = NULL;

      
      }

      zdvoxels += numberOfPixelsZ;
    }
    // track up 

  }

   */

// ===================================================================
 /*
  if( ( (sObj->type) != WLZ_CONTOUR ) || ( (tObj->type) != WLZ_CONTOUR ) )
  {
      	errNum = WLZ_ERR_DOMAIN_TYPE;
  }

  if(  ( sObj->domain.core == NULL ) || ( tObj->domain.core == NULL ) )
  {
	errNum = WLZ_ERR_DOMAIN_NULL;
  }
  
  // get one loop line points from GM one
  if(errNum == WLZ_ERR_NONE)
  {

    printf("-----Standard -----\n");
    // get the G Model for standard
    gM        = sObj->domain.ctr->model;

    // get the stand points from on loop of the GM
    LoopIdx = -1;
    GetSamplePointsFromOneLoopOfGM( gM, 
                                    LoopIdx,
                                    StandSampleP,
                                    &errNum );
  }

   // get a KD Tree for the target vertices of the Woolz contour
     tTree = GetTreeFromGM(tObj, &errNum);

  
    if( errNum == WLZ_ERR_NONE   )
    {
 
          printf("------Target -------\n");
          // get the G Model for target to be tracked
          gM2 = tObj->domain.ctr->model;

	  // get the tracked points in a loop and get the loop index and
	  // number of sample points
	  GetTrackedSamplePointsFromOneLoopOfGM(gM2, 
	                                        tTree, 
					       &ntL, 
						numberOfPixelsZ,
	                                        minDis,
	                                        StandSampleP,
						TrackedSamplePBlgToOneLp,
						&nOfTracked,
						&suc,
						&errNum);
	
        if(errNum == WLZ_ERR_NONE)
        {

           // output the surface tracked:
           if((testFile = fopen(surfacePointFileName, "w")) == NULL )
           {
		     errNum =  WLZ_ERR_FILE_OPEN;
           }
           if(errNum == WLZ_ERR_NONE)
           {
	       idx = 0;
  	       for(i=0; i<nOfTracked; i++)
	       {
                   fprintf(testFile, "%6.2f  %6.2f   %6.2f\n",
		   TrackedSamplePBlgToOneLp[i].vtX, 
	           TrackedSamplePBlgToOneLp[i].vtY, 
		   (double) numberOfPixelsZ  );
	       }
	   }
	   
	   // get the loop index we have just tracked
           distance = 10.;
	   // calculate the in and out surface points 
           InAndOutSurfacePoints(gM2, ntL, numberOfPixelsZ, distance,infp, outfp, &errNum);

	}
    */
   }

   if(testFile)
      fclose(testFile);
   fclose(infp);
   fclose(outfp);
   fclose(surfacefp);
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
        if(cET == cET->edge->edgeT) // Print edge end points 
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
       //printf("Next loop\n");
       do
       {
           printf("   Loop:  %d  Loop idx %d\n", nLoop, cLT->idx);
           /* For each edge topology element of the model. */
	   cET = fET = cLT->edgeT;
   	   do
	   {
              if(cET == cET->edge->edgeT) // Print edge end points 
	      {
	         nVertices += 2;
	      }
	      cET = cET->next;
	   } while (cET != fET);
	   // if( nVertices > 50)
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
       //printf("Next loop\n");
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
   WlzDVertex2   *tDVP,
  		 *vData = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  vIdx = 0;

  cnt = nOfEdgeEndVInOneLoopOfGM(model, nThShell, nThLoop);
  if((vData = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) * cnt)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    // go to the correction position

 
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
        if(cET == cET->edge->edgeT) // Print edge end points 
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
   WlzDVertex2   *tDVP,
  		 *vData = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  vIdx = 0;

  cnt = nOfVInOneLoopOfGM(model, nThShell, nThLoop);
  if((vData = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) * cnt)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    // go to the correction position

 
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
        // if(cET == cET->edge->edgeT) // Print end points 
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
  WlzDVertex3   *tDVP,
  		*vData = NULL;
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
    // go to the correction position

 
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
        if(cET == cET->edge->edgeT) // Print edge end points 
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
  WlzDVertex3   *tDVP,
  		*vData = NULL;
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
    // go to the correction position

 
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
        //if(cET == cET->edge->edgeT) // Print edge end points 
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
   // printf("here\n");
   cLT     =   (WlzGMLoopT *) cVT->diskT->vertexT->parent->parent;
  *nLoop   =    cLT->idx;
   cST     =    cLT->parent;
   *nShell =    cST->idx;
   AlcFree(vec);
   AlcFree(cVT);
   AlcFree(cLT);
   AlcFree(cST);
   // printf("Shell:  %d  loop:  %d\n", *nShell, *nLoop);
}

static int NumberOfLoopsInThe_nThShellOfGM(WlzGMModel *gM, int nThShell )
{
   int nLoop;
   AlcVector          *vec;
   WlzGMShell         *cS;
   WlzGMLoopT         *cLT, *fLT;
   
   vec    =  gM->res.shell.vec;

   nLoop     = 0;
    /* For each shell of the model. */
    cS = (WlzGMShell *) AlcVectorItemGet(vec, nThShell);
    cLT = fLT = (WlzGMLoopT *) cS->child;
       //printf("Next loop\n");
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
    /* For each LoopT of the model. */
    cLT = (WlzGMLoopT *) AlcVectorItemGet(vec, sLoopIdx);
    cET = fET = cLT->edgeT;
    do
    {

        nVertices += 1;
	cET = cET->next;
    } while (cET != fET);
   
   *dstErr = errNum;

    if(IsACycle(gM, sLoopIdx, &errNum))
    {
               return (nVertices);
    }	   
    return ( nVertices / 2  + 1 );

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
   WlzDVertex2   *tDVP,
  		 *vData = NULL;
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
	    if(test)
	    printf(" idx %d  x %lg   y %lg\n", cET->vertexT->diskT->idx,
	                                       (vData + nVertices)->vtX,
	                                       (vData + nVertices)->vtY );
            nVertices += 1;
	}
	//if( ( cET->vertexT->diskT->idx == cET->next->next->vertexT->diskT->idx ) )
	//      break;
	cET = cET->next;
    } while (cET != fET);

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
   WlzDVertex3   *tDVP,
  		 *vData = NULL;
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
        
	//if( ( cET->vertexT->diskT->idx == cET->next->next->vertexT->diskT->idx ) )
	//      break;
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
    int         *sShellIdx, *sLoopIdx;
    WlzVertexP   sLoopVData;
    WlzErrorNum   errNum = WLZ_ERR_NONE;
    // get the number of Shells in this model
     nS        = NumberOfShellsAboutGM(gM);

    // get the Index for each Shell;
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
             AlcFree(sLoopVData.v);					   
         }
         AlcFree(sLoopIdx);
     }
     AlcFree(sShellIdx);
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
     
     // get the number of Element in each group
     // and get whether it's continue ?
     if(errNum == WLZ_ERR_NONE)
     {
        jg = 0;
	j  = 0;
        for(i=0; i< IthGroup; i++)
	{
          for(k=jg; k<nt; k++)
	  {
            // 
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
       //get the correct number
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
    AlcFree(groupI);

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
				   int         *nOfVInTheInFile,
				   WlzDVertex3 *OutPoints,
				   int         *nOfVInTheOutFile,
				   int         *suc,
				   WlzErrorNum *dstErr )
{
    int i, j, k, m,  idx, ntV, loopDirection, di;
    int inN, outN;
    WlzVertexP         tLoopVData;

    WlzErrorNum errNum = WLZ_ERR_NONE;
    WlzDVertex2         segV[3],  fPoint, bPoint, *vNorm=NULL, *inPoints=NULL, *outPoints=NULL;
    inN  = *nOfVInTheInFile;
    outN = *nOfVInTheOutFile;
    *suc  = 1;
 
	   // calculate the in and out surface points 
        
	   // get the number of vertices in this loop
           ntV  =  NumberOfVerticesInTheLoopOfGM(gM2, nThLoop, &errNum); 

           if( errNum == WLZ_ERR_NONE   )
           {
              // tLoopVData = WlzVerticesOneLoopFromGM(gM2,  ntS, ntL, &errNum);
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
	    // determin the direciton and goto the position:
            // get the position of the index
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

	     // calculate in and out points of the surface
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

		  
		  // printf("%lg   %lg\n", segV[1].vtX, segV[1].vtY );
	         *(vNorm + idx) = WlzVerticesNormTriple2(segV[0], segV[1], segV[2]);
		  // get the in and out points using the normal and points
		  fPoint.vtX =  segV[1].vtX + distance * (vNorm+idx)->vtX;
		  fPoint.vtY =  segV[1].vtY + distance * (vNorm+idx)->vtY;
		  bPoint.vtX =  segV[1].vtX - distance * (vNorm+idx)->vtX;
		  bPoint.vtY =  segV[1].vtY - distance * (vNorm+idx)->vtY;
		  // check its out or in
		  di = WlzGeomTriangleSnArea2(segV[0], segV[1], fPoint); 
		  if( ( ( di >  0. ) && ( loopDirection == 1) ) ||
		      ( ( di <  0. ) && ( loopDirection == 0) )
		  )
		  {
		     // fPoint is on the left of line segV[0]->segV[1]
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
		

                //} while(idx < NumberToTrack/4 );
                } while(idx < 1 );
			
		//store the in points 
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

		//store the out points
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

      if(tLoopVData.d2)
              AlcFree(tLoopVData.d2);
      if(vNorm)
      {
              AlcFree(vNorm);
	      vNorm     = NULL;

      }	      
      if(inPoints)
      {
              AlcFree(inPoints);
              inPoints  = NULL;

      }
      if(outPoints)
      {
             AlcFree(outPoints);
      	     outPoints = NULL;

      }
      *nOfVInTheOutFile = outN;
      *nOfVInTheInFile  = inN;

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

    // first 
    if(LoopIdx < 0)
    {
        // get the number of Shells in this model
        nS        = NumberOfShellsAboutGM(gM);
        // get the Index of for each Shell;
        sShellIdx = WlzShellsIndexAboutGM(gM, nS, &errNum);


        if( errNum == WLZ_ERR_NONE   )
        {
           test = 0;
           if(test)
           {
              WlzOutputForTestGM( gM, &errNum);
              exit(0);
           }

           // use the first Shell and first loop as a test
           i        = 0;


           // get number loops in this shell
           nL       =  NumberOfLoopsInThe_nThShellOfGM(gM,  sShellIdx[i] );

           // get the index for the loops
           sLoopIdx =  WlzLoopIndexAboutGM(gM, sShellIdx[i], nL, &errNum);
           

	  if( errNum == WLZ_ERR_NONE )
	  {
             // get howManyVertices in the loop:
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
	if(sShellIdx)
	   AlcFree(sShellIdx);
	if(sLoopIdx) 
	   AlcFree(sLoopIdx);

    }

    if( errNum == WLZ_ERR_NONE   )
    {
        
        nV       = NumberOfVerticesInTheLoopOfGM(gM, LoopIdx, &errNum); 
        if( errNum == WLZ_ERR_NONE )
	{
            //printf("num Vertices = %d\n",nV);
             sLoopVData = WlzVerticesThisLoopOfGM( gM, LoopIdx, nV, &errNum  );
               
	     if( errNum == WLZ_ERR_NONE )
	     {
	       // output for forming surface
               test = 0;
               if(test)
               {
                    WlzOutputForTestAsectionData( sLoopVData, nV, &errNum );
                    exit(0);
               } 
	     }
         }
      }
      //get the stand sample points first
      if(nV < 62)
      {
           printf("This section is too short, please change a section and try again!\n");
	   exit(0);
      }
      //j = nV/NumberToTrack;
      j = 3;

      for(i=0; i<NumberToTrack; i++)
      {
	   // using sample position, sp should not be too big as it will push 
	   // the index out of range!!!:
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
        
        int                 i,j, idx, ntL;
	int                 test = 0;
	int                 nShell, nLoop;
	int                 Continue;
        double              datD[3];
        AlcKDTNode         *node;
	AlcErrno            errNo = ALC_ER_NONE;
        int                 LoopIndex[NumberToTrack];
        int                 ShellIndex[NumberToTrack];
        int                 indexStr[NumberToTrack]; 
        int                *GroupNum, *FixNumForEachGroup;
        int                 NumOfGroup;
        int                 sp, nbgN, trackedBigestLoopIndex;
        int                 loopDirection;
        WlzDVertex2         TrackedSampleFirstP, TrackedSampleSecondP;
	WlzDVertex2         TrackedSampleP[NumberToTrack];
        WlzErrorNum         errNum = WLZ_ERR_NONE;
	*suc = 1;
        // --- find the nearst node (one or nil )----
	for(i=0; i<NumberToTrack; i++)
	{
          datD[0] = StandSampleP[i].vtX;
          datD[1] = StandSampleP[i].vtY;
	  datD[2] = 0.;
          // printf("Point %lg %lg\n",datD[0], datD[1]);
          node    = AlcKDTGetNN(tTree, datD, minDis, NULL, &errNo);
     	  if(  node && (errNum == ALC_ER_NONE) )
  	  {
             TrackedSampleP[i].vtX = *(node->key.kD);
             TrackedSampleP[i].vtY = *(node->key.kD + 1);
	     indexStr[i]           = node->idx;
	     AlcFree(node);
  	  }
	  else
	  {
             printf(" cannot track this points!!! you have \
	              to give a biger mindistance \
	              for nearest point to track!\n");
	     *suc = 0;
	     return;
	  }
          //  get the loop corresponding  to the tarcked point
          //  it may failed as the crossed point belongs to
	  //  to loops !!
	  //printf("index = %d\n",indexStr[i] );
	  if(indexStr[i] > 50000)
	  {
           printf("Something wrong\n");
	     AlcFree(node);
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

	// check, how many loops in the target group:
        // check, whether it continue ?:
	Continue = 0;
	// get the accurate number of points for each group
	// and answer the continunities of the data
	GroupNum           = GroupDivid(LoopIndex, &NumOfGroup, &Continue, NumberToTrack, &errNum);
	
        if(errNum == WLZ_ERR_NONE)
        {
	  // get the fixed number for each group
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


          // we have to choose the right one
          nbgN = BigestGroup(GroupNum, NumOfGroup, &errNum);
	  if(test)
	     printf("group  %d is the bigest group\n",nbgN);
	     
	  if(GroupNum[nbgN-1] < (NumberToTrack* 65)/100 )
	  {
              printf("track failed as no group is a dominant one\n");
	      AlcFree(GroupNum);
	      AlcFree(FixNumForEachGroup);
              *suc = 0;
	  }

	}

	//
	/*
	printf("tracked points belongs to one loop!\n");

        if((testFile = fopen(surfacePointFileName, "w")) == NULL )
        {
		     errNum =  WLZ_ERR_FILE_OPEN;
        }
	*/
        if(errNum == WLZ_ERR_NONE)
        {
	     idx = 0;
  	     for(i=0; i<NumberToTrack; i++)
	     {
	         if(FixNumForEachGroup[i] == nbgN)
		 {
		   TrackedSamplePBlgToOneLp[idx].vtX =  TrackedSampleP[i].vtX;
		   TrackedSamplePBlgToOneLp[idx].vtY =  TrackedSampleP[i].vtY;
                   //fprintf(testFile, "%6.2f  %6.2f   %6.2f\n",TrackedSampleP[i].vtX, 
	           //                        TrackedSampleP[i].vtY, (double) numberOfPixelsZ);
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
	    // get the loop index we have just tracked
	    for(i=0; i<NumberToTrack; i++)
	    {
              if(FixNumForEachGroup[i] == nbgN )
	      {
                  // printf("Loop index tracked = %d\n", LoopIndex[i]);
	   	  ntL = LoopIndex[i];
		  break;

	      }

	    }
       }
       *LoopIdx = ntL;

    AlcFree(GroupNum);
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
    int                *sNN, test;
    WlzVertexP         nData[2], vData[2];
    WlzVertexType       vType;
    WlzErrorNum         errNum = WLZ_ERR_NONE;
        // get the vertices from obj

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
            test = 0;
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
	 // build up a KD Tree for the target vertices of the Woolz contour 
         if(  (sNN = (int *)AlcMalloc(sizeof(int) * vCnt[1])) == NULL )
         {
              errNum = WLZ_ERR_MEM_ALLOC;
         }

         if(errNum == WLZ_ERR_NONE)
         {
              tTree = WlzVerticesBuildTree(vType, vCnt[1], vData[1], sNN, &errNum);
         }
    }

    AlcFree(sNN);
    AlcFree(nData[0].d2);
    AlcFree(vData[0].d2);
    AlcFree(nData[1].d2);
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
* \param       *gM                          input geometric model 
* \param        TrackedSamplePBlgToOneLp[NumberToTrack] input tracked points
* \param        StandSampleP[NumberToTrack] output sampled points
* \param       *dstErr			    Woolz err number pointer.
*/
static void GetTrackedPlusSomeSamplePointsFromOneLoopOfGM( 
				   WlzGMModel  *gM, 
                                   int          LoopIdx, 
                                   WlzDVertex2  StandSampleP[NumberToTrack],
				   WlzDVertex2	TrackedSamplePBlgToOneLp[NumberToTrack],
				   int		nOfTracked,
				   int         *suc,
                                   WlzErrorNum *dstErr )
{
    int   i, j, k, m, sp, nS, nV, nL, test;
    int  *sShellIdx=NULL, *sLoopIdx=NULL;
    WlzErrorNum errNum = WLZ_ERR_NONE;
    WlzVertexP          sLoopVData;
    *suc = 1;

    for(i=0; i<NumberToTrack; i++)
    {
          StandSampleP[i].vtX =  TrackedSamplePBlgToOneLp[i].vtX;
          StandSampleP[i].vtY =  TrackedSamplePBlgToOneLp[i].vtY;
    }
    
    if(nOfTracked >= NumberToTrack)  
    {
      return;
    }
    else
    {
       nV       = NumberOfVerticesInTheLoopOfGM(gM, LoopIdx, &errNum); 
       if( errNum == WLZ_ERR_NONE )
       {
             sLoopVData = WlzVerticesThisLoopOfGM( gM, LoopIdx, nV, &errNum  );
       }
       //get the stand sample points first
       if(nV < 62)
       {
           printf("This section is too short, please change a section and try again!\n");
	   AlcFree(sLoopVData.d2);
	  *suc = 0;
	   return; 
       }

       // get the position of the index
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
	   // using sample position, sp should not be too big as it will push 
	   // the index out of range!!!:
	   sp = 1;
	   k  = i +  ( m - nOfTracked + 1) * j;
	   if(k < nV)
	   {
             StandSampleP[m].vtX =  (sLoopVData.d2 + k )->vtX;
             StandSampleP[m].vtY =  (sLoopVData.d2 + k )->vtY;
	   }
	   else
	   {
	        AlcFree(sLoopVData.d2);
                printf("not that much points to add in the loop");
		*suc = 0;
		return;
	   }
       }
    }   
    AlcFree(sLoopVData.d2);

}

/*!
*  \return      int 0 for non-cycle (shell);
*  \ingroup     WlzFeatures
*  \brief       Determine whether the line(shell) relating with this loop is a cycle.
*               a better way to do it seems to know how many loops with a
*               shell in a given model.
*               if a shell has two loops, it's a cycle. one loops is a
*               section.
*  \param       gM                  given G model (input)
*  \param       nThLoop             n-th loop of the n-th shell in the model (input )
*/
static int IsACycle(   WlzGMModel  *gM, 
                       int          LoopIdx, 
                       WlzErrorNum *dstErr )
{
   int cycle, nLoop;
   WlzErrorNum errNum = WLZ_ERR_NONE;
   int nVertices;
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


    cLT = fLT = cS->child;
   /* For each loop topology element of the model. */
     do
     {
        cLT = cLT->next;
        ++nLoop;
     } while(  cLT != fLT );

    if(nLoop > 1)
       cycle = 1;
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

    // counter again
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

       // counter again
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
	      // check whether it already in the bank
	      for(i=0; i<numberOfBifurP; i++)
	      {
                  if( *(BifurcatePointsIndex + i ) == cET->vertexT->diskT->idx )
		        break;
	      }
	      if(i == numberOfBifurP)
	      {
	         // not in the bank
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

    // first 

    //get the stand sample points first
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
* \param        nOfVInTheSFile          Given the number of points on the surface. 
* \param       *InPoints                Given pointer points to the surface points.
* \param        nOfVInTheSFile          Given the number of points on the surface. 
* \param       *OutPoints               Given pointer points to the points in surface.
* \param        nOfVInTheSFile          Given the number of points outside  the surface. 
* \todo         
*/
static void	outputSurfaceInAndOutPoints(FILE          *surfacefp, 
					    FILE          *infp,  
					    FILE          *outfp, 
		                            WlzDVertex3   *SurfacePoints, 
					    int            nOfVInTheSFile,
					    WlzDVertex3   *InPoints,      
					    int            nOfVInTheInFile,
					    WlzDVertex3   *OutPoints,     
					    int            nOfVInTheOutFile)
{

     int k;

       		     for(k=0; k<nOfVInTheSFile; k++)
       		     {
                   	  	fprintf(surfacefp, "%6.2f  %6.2f   %6.2f\n",
		   		SurfacePoints[k].vtX, 
	           		SurfacePoints[k].vtY, 
	           		SurfacePoints[k].vtZ ); 

       		     }
       		     for(k=0; k<nOfVInTheInFile; k++)
       		     {
                   	  	fprintf(infp, "%6.2f  %6.2f   %6.2f\n",
		   		InPoints[k].vtX, 
	           		InPoints[k].vtY, 
	           		InPoints[k].vtZ ); 

       		     }
       		     for(k=0; k<nOfVInTheOutFile; k++)
       		     {
                   	  	fprintf(outfp, "%6.2f  %6.2f   %6.2f\n",
		   		OutPoints[k].vtX, 
	           		OutPoints[k].vtY, 
	           		OutPoints[k].vtZ ); 

       		     }





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
				   int         *nOfVInTheInFile,
				   WlzDVertex3 *OutPoints,
				   int         *nOfVInTheOutFile,
				   WlzErrorNum *dstErr )
{
    int i, j, k, m,  idx, ntV, loopDirection, di;
    int inN, outN;

    WlzErrorNum errNum = WLZ_ERR_NONE;
    WlzDVertex2         segV[3],  fPoint, *vNorm=NULL, *inPoints=NULL, *outPoints=NULL;
    inN  = *nOfVInTheInFile;
    outN = *nOfVInTheOutFile;
 
     // calculate the in and out surface points 
	    
     idx = 0;
     i   = NumberToTrack/4;
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
	  // printf("%lg   %lg\n", segV[1].vtX, segV[1].vtY );
	 *(vNorm + idx) = WlzVerticesNormTriple2(segV[0], segV[1], segV[2]);
	  // get the in and out points using the normal and points
	  fPoint.vtX =  segV[1].vtX + distance * (vNorm+idx)->vtX;
	  fPoint.vtY =  segV[1].vtY + distance * (vNorm+idx)->vtY;
	  // check its out or in
	  di = WlzGeomTriangleSnArea2(segV[0], segV[1], fPoint); 
	  if(   di >  0. ) 
	  {
	     // fPoint is on the left of line segV[0]->segV[1]
	     (inPoints +idx)->vtX =  fPoint.vtX; 
	     (inPoints +idx)->vtY =  fPoint.vtY;
	     
	     (outPoints +idx)->vtX =   segV[1].vtX - distance * (vNorm+idx)->vtX;
	     (outPoints +idx)->vtY =   segV[1].vtY - distance * (vNorm+idx)->vtY;
	     i += 4;
	     idx++;	     
  
	  }
	  else if(  di <  0. ) 
	  {
	     (outPoints +idx)->vtX =  fPoint.vtX; 
	     (outPoints +idx)->vtY =  fPoint.vtY;
	     
	     (inPoints +idx)->vtX =   segV[1].vtX - distance * (vNorm+idx)->vtX;
	     (inPoints +idx)->vtY =   segV[1].vtY - distance * (vNorm+idx)->vtY;
             i += 4;
	     idx++;	
	  }
	  else
	  {
              printf("On the surface!!!!!!!!!!!!! not counted it\n");
	      i++;
	      if(i + 3 > NumberToTrack)
	          break;

	  }
          //} while(idx < NumberToTrack/4 );
     } while(idx < 1 );
			
     //store the in points 
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

     //store the out points
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
              AlcFree(vNorm);
	      vNorm     = NULL;

      }	      
      if(inPoints)
      {
              AlcFree(inPoints);
              inPoints  = NULL;

      }
      if(outPoints)
      {
             AlcFree(outPoints);
      	     outPoints = NULL;

      }
      *nOfVInTheOutFile = outN;
      *nOfVInTheInFile  = inN;

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
   
        // get the number of Shells in this model
        numOfShellsInStartGM  = NumberOfShellsAboutGM(gMC);
        // get the Index for each Shell;
        startShellIdx = WlzShellsIndexAboutGM(gMC, numOfShellsInStartGM, &errNum);
        
	if(errNum == WLZ_ERR_NONE)
	{
          for(i=0; i<numOfShellsInStartGM; i++)
	  {
	     // get number loops in this shell
             numOfLoopsInStartShell  =  NumberOfLoopsInThe_nThShellOfGM(gMC,  startShellIdx[i] );
	     printf("The  %d th Shell\n",  i );

             // get the index for the loops
             startLoopIdx            =  WlzLoopIndexAboutGM(gMC, startShellIdx[i], numOfLoopsInStartShell, &errNum);
	     
	     if(errNum == WLZ_ERR_NONE )
	     {
                for(j=0; j<numOfLoopsInStartShell; j++)
		{
		   //outputVerticesInThisLoop(gMC, startLoopIdx[j], &errNum);
		   if( errNum != WLZ_ERR_NONE )
		          break;
	           // check is it a cycle
		   if(!IsACycle(gMC, startLoopIdx[j], &errNum))
		   {
                        if(errNum == WLZ_ERR_NONE )
			{
			   numOfEnds = HowManyEndsInTheNonCycleLine( gMC, 
                                                            startLoopIdx[j], 
                                                            &errNum );
			  if( ( numOfEnds > 2 ) && (errNum == WLZ_ERR_NONE) ) 
			  {
                             //pure this cycle
                               EndsIdx =   GetTheEndsPointsIndexInTheNonCycleLine( gMC, 
                                           startLoopIdx[j],
					   numOfEnds,
                                           &errNum );
					   
			       bifurIdx = GetThebifurcatePointsIndexInTheNonCycleLine(   gMC, 
                                           startLoopIdx[j],
					   numOfEnds,
                                           &errNum );
			       pureThisGM(gMC, startLoopIdx[j], numOfEnds, EndsIdx, bifurIdx, &errNum);	
			       
                               AlcFree(EndsIdx);
			       AlcFree(bifurIdx);

			  }   
			     
			}

		   }

		}
	     }
	     if(startLoopIdx)
	     AlcFree(startLoopIdx);  

	  }

	}

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
	printf("%6.2f   %6.2f\n",cET->vertexT->diskT->vertex->geo.vg2D->vtx.vtX,
	                         cET->vertexT->diskT->vertex->geo.vg2D->vtx.vtY);
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
static void pureThisGM(WlzGMModel *gMC, int LoopIdx, int numOfEnds, int *EndsIdx, int *bifurIdx, WlzErrorNum *dstErr)
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
   // delete bifurcation one by one but from the shortest to the longest
   for(i=0; i<numOfEnds; i++)
   {
      do
      {
  	  cET = cET->next;
	  if(    cET->vertexT->diskT->vertex->idx  ==   *(bifurIdx+i)  )
	       break;
      } while (cET != fET);
      // check its length to the nearst bifurcation points

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
      // find the shortest braches from the bifurcation points
      it = 0;
      d1 = 5;
      aET = cET->opp->next;
      bET = ffET = cET;
      do
      {
          // walking along three directions
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
* \param        distance                Given distance of the in or out points to the surface.
* \param       *SurfacePoints           Given pointer points to the surface points.
* \param        nOfVInTheSFile          Given the number of points on the surface. 
* \param       *InPoints                Given pointer points to the surface points.
* \param        nOfVInTheSFile          Given the number of points on the surface. 
* \param       *OutPoints               Given pointer points to the points in surface.
* \param        nOfVInTheSFile          Given the number of points outside  the surface. 
* \todo         
*/
static void RecalculateInAndOutPointsByUsing3DNormal(
				            double         distance,
		                            WlzDVertex3   *SurfacePoints, 
					    int            nOfVInTheSFile,
					    WlzDVertex3   *InPoints,      
					    int           *nOfVInTheInFile,
					    WlzDVertex3   *OutPoints,     
					    int           *nOfVInTheOutFile )
{
   int i, j, k, l, count;
   WlzDVertex3 p[4], normal, inOrOut[2];
   double volume;
	k = NumberToTrack/2;
	l = *nOfVInTheInFile - 3;

     // recalculate the in points and out points
     for(i=0; i<l; i++)
     {
        j = 0;
        // take three points from the surface:
	for(j=0; j<k; j++)
	{
	   WlzValueCopyDVertexToDVertex3(p,  (SurfacePoints + j+1), 1);
	   WlzValueCopyDVertexToDVertex3(&p[1],  (SurfacePoints + j + k - 1), 1);
	   if( !IsDVertexEqual(p[0], p[1]) )
	          break;
        }
	
	WlzValueCopyDVertexToDVertex3(&p[2],  (SurfacePoints + (i+2) * NumberToTrack +j+ 5), 1);
        // calculate the normal			
	normal = WlzGeomTriangleNormal(p[0], p[1], p[2]);
	if( !IsDVertexZero(normal) )
	{
	     // detremine the in and out points:
	     WlzValueCopyDVertexToDVertex3( &p[3], ( InPoints + j), 1 ); 
             volume = SixTimesOfVoulumeOfTetraHedron( p[0], p[1], p[2], p[3]); 
             if( volume != 0.)
	     {
	       inOrOut[0].vtX = p[2].vtX + distance * normal.vtX;
	       inOrOut[0].vtY = p[2].vtY + distance * normal.vtY;
	       inOrOut[0].vtZ = p[2].vtZ + distance * normal.vtZ;
	       
	       inOrOut[1].vtX = p[2].vtX - distance * normal.vtX;
	       inOrOut[1].vtY = p[2].vtY - distance * normal.vtY;
	       inOrOut[1].vtZ = p[2].vtZ - distance * normal.vtZ;
		
		if(volume > 0)
		{
                       (InPoints + count)->vtX  = inOrOut[0].vtX;
                       (InPoints + count)->vtY  = inOrOut[0].vtY;
                       (InPoints + count)->vtZ  = inOrOut[0].vtZ;

                       (OutPoints + count)->vtX = inOrOut[1].vtX;
                       (OutPoints + count)->vtY = inOrOut[1].vtY;
                       (OutPoints + count)->vtZ = inOrOut[1].vtZ;
		}
		else
		{
                       (InPoints + count)->vtX  = inOrOut[1].vtX;
                       (InPoints + count)->vtY  = inOrOut[1].vtY;
                       (InPoints + count)->vtZ  = inOrOut[1].vtZ;

                       (OutPoints + count)->vtX = inOrOut[0].vtX;
                       (OutPoints + count)->vtY = inOrOut[0].vtY;
                       (OutPoints + count)->vtZ = inOrOut[0].vtZ;

		}
	        // replace the normal with new one!
                count++;
	     }	
	}
        




     }

     *nOfVInTheInFile   = count;
     *nOfVInTheOutFile  = count;



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
