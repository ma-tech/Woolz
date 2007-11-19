/***********************************************************************
* \Project:      Woolz
* \Title:        WlzInterpolateBibFiles.c
* \Date:         1 Aug 2002 
* \Author:       J. Rao
* \Copyright:	 2002 Medical Research Council, UK.
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

#include <WlzExtFF.h>
#include <bibFile.h>

/* #include <MAPaint.h> */
/* #include <MAWarp.h> */

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


static void usage(char *proc_str);  

static WlzErrorNum ReadBibFile(FILE *inFile, char *inFileStr, WlzThreeDViewStruct  *wlzViewStr
);

static WlzErrorNum WriteBibFile(FILE *inFile,    char *inFileStr,  char *inFileStrw, char *inFileStrSec,
WlzThreeDViewStruct  *wlzViewStr  );


static WlzErrorNum     LinearInterpolations( WlzThreeDViewStruct *wlzViewStr, WlzThreeDViewStruct  *wlzViewStr1, 
                                       WlzThreeDViewStruct *wlzViewStrInter, int Inumber, int TotalN );

static WlzErrorNum     ChangeFixedPoint(WlzThreeDViewStruct *wlzViewStr, double x, double y, double z );

int main(int	argc,
	 char	**argv)
{

  FILE	       *inFile = NULL;   /* used to read Woolz object */ 
  FILE	       *outFile = NULL;   /* used to read Woolz object */
  int           i,j,k=0, m; /* , l1, label[3], individual[4] */
  int		option;
  int            TotalN; /* numChar, numChar1, */
  
  /* int           outputAutoMeshVTK        = 0,
                outputAutoMeshWLZ        = 0,
                outputTransformedMeshWLZ = 0,
    outputCutPlaneAndCorrepSurfaceVTK        = 0;
  */
    /*		WARP                     = 1 */
  /* int                  basisFnPolyOrder = 3; */
  int                  Inumber, num;
  int                  globalCycle;
  /* int                  binObjFlag; */
  int                  centorOfMass = 1;
  int                  numOf2DWlzFiles = 0, numOfSampleBibFiles = 0;
  /* int                  initialn0, endNum0,initialn1, endNum1; */
  int           BibFileIndex[50];
  /* double	zConst                   = 0.; */
  /* double        mass = 0.0; */
  /* char          under ='_'; */
  /* char                ctemp; */
  char                *inFileStr, *inFileStr1,  *outFileStr; /* *inFileStr3, */
  char                *inFileStrw, *outputBibFilesDir; /* , inFileStrSec[100]; */
  char                TwoDImageFilesNameList[1000][120];
  char                TwoDImagePureFilesNameList[1000][70];
  char                PureBibFilesNameList[700][70];
  char                SampleBibFilesNameList[100][120];
  char                SampleBibPureFilesName[100][70];
  char                SampleBibFilesDir[100];
  char                TwoDImageFilesDir[100];
  char                *List2DImageFilesStr, *SampleBibFilesStr;
  /* char                *cstr; */
  /* const char	      *errMsg; */
  /* WlzDVertex2          cMass, cMassS; */
  WlzErrorNum	       errNum = WLZ_ERR_NONE;
  /* WlzInterpolationType interp = WLZ_INTERPOLATION_NEAREST;  Use the nearest neighbour */
  WlzThreeDViewStruct *wlzViewStr, *wlzViewStr1, *wlzViewStrInter;
  /* AlcErrno             alcErr = ALC_ER_NONE; */
  WlzObject           *WObjS;/*, *WObj2DS, *WObj2D;*/ 
  /* WlzObjectType        wtp = WLZ_3D_DOMAINOBJ; */

  /* read the argument list and check for an input file */
  static char	optList[] = "i:I:f:w:o:d:z:M:m:n:r:R:L:B:h",

  opterr = 0;
 
   while( (option = getopt(argc, argv, optList)) != EOF )
   {
      switch( option )
      {
        case 'h':
             usage(argv[0]);
	     return(0);
        case 'i':
	    inFileStr  = optarg;
	    /*
	    outFileStr = optarg;
	    */
	    break;
        case 'I':
	    inFileStr1 = optarg;
	    break;
        case 'L':
	    List2DImageFilesStr = optarg;
	    break;
        case 'B':
	    SampleBibFilesStr = optarg;
	    break;
        case 'w':
	    inFileStrw = optarg;
	    break;
        case 'n':
	 if(sscanf(optarg, "%d", &num) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;
        case 'o':
	    outFileStr = optarg;
	    break;
        case 'd':
	    outputBibFilesDir = optarg;
	    break;
        default:
              return(0);
      }
   }


  /* Read 2D Wlz files List */
     printf("Read 2D Wlz file list\n");
    if((inFile = fopen(List2DImageFilesStr, "r")) == NULL )
    {
        printf("cannot open the 2D image files list file.\n");
        exit(1);
    }
    
    i=0;
    /* get the number of files */
    while(!feof(inFile)){

        fscanf(inFile, "%s", *(TwoDImageFilesNameList + i ));
	/*  printf("%s\n", *(TwoDImageFilesNameList + i ) ); */
	i++;
    }
      fclose(inFile); 
      inFile = NULL;

    numOf2DWlzFiles = i-1;
    printf("number of files = %d\n", numOf2DWlzFiles );

    /* extract pure wlz filename without directory */
    printf("extract pure wlz filename without directory\n");

    for(m=0; m<numOf2DWlzFiles; m++)
    {
       for(i=0; i<sizeof(TwoDImageFilesNameList[m]); i++){
          if( TwoDImageFilesNameList[m][i] == (char  )NULL )
           break; 
          if( TwoDImageFilesNameList[m][i] == '/' )
            k = i; 	   
       }
       if(k != 0)
       k++;
       for(j=k; j<i; j++)
       {
         TwoDImagePureFilesNameList[m][j-k] = TwoDImageFilesNameList[m][j];
         PureBibFilesNameList[m][j-k]       = TwoDImageFilesNameList[m][j];
       }
       TwoDImagePureFilesNameList[m][j-k] = (char ) NULL;
       PureBibFilesNameList[m][j-k]       = (char ) NULL;
       PureBibFilesNameList[m][j-k-3]     = 'b';
       PureBibFilesNameList[m][j-k-2]     = 'i';
       PureBibFilesNameList[m][j-k-1]     = 'b';
         
       /* extract directory */
       if(m == 0)
       {
          for(j=0; j<k; j++)
	  {
	     TwoDImageFilesDir[j] = TwoDImageFilesNameList[m][j];
	  }
	     TwoDImageFilesDir[j] = (char ) NULL;
	     printf("I_DIR  %s\n",TwoDImageFilesDir );
       }
       printf("%s\n", TwoDImagePureFilesNameList[m]);   
    }

    /* changeToBibFileNamse */

    /*--- Read Sample bib File name list ----*/
    if((inFile = fopen(SampleBibFilesStr, "r")) == NULL )
    {
        printf("cannot open the sample bib files list file.\n");
        exit(1);
    }
    i=0;
    while(!feof(inFile)){
        fscanf(inFile, "%s", SampleBibFilesNameList[i]);
	printf("%s\n", SampleBibFilesNameList[i]);
	i++;
    }
      fclose(inFile); 
      inFile = NULL;
    numOfSampleBibFiles = i-1;
    printf("number of files = %d\n", numOfSampleBibFiles );


     /* extract pure sample bib filename without directory */
    printf("extract pure sample bib filename without directory\n");

    for(m=0; m<numOfSampleBibFiles; m++)
    {
       for(i=0; i<sizeof(SampleBibFilesNameList[m]); i++){
          if( SampleBibFilesNameList[m][i] == (char ) NULL )
           break; 
          if( SampleBibFilesNameList[m][i] == '/' )
            k = i; 	   
       }
       if(k != 0)
       k++;
       for(j=k; j<i; j++)
       {
         SampleBibPureFilesName[m][j-k] = SampleBibFilesNameList[m][j];
       }
        SampleBibPureFilesName[m][j-k] = (char ) NULL;
       /* extract directory */
       if(m == 0)
       {
          for(j=0; j<k; j++)
	  {
	     SampleBibFilesDir[j] = SampleBibFilesNameList[m][j];
	  }
	     SampleBibFilesDir[j] = (char) NULL;
	     printf("B_DIR  %s\n",SampleBibFilesDir);
       }

       printf("%s\n", SampleBibPureFilesName[m]);   
    }


    /*---------- extract the bib File number ----------- */
    for(i=0; i<numOfSampleBibFiles; i++)
    {
       BibFileIndex[i] = -1;
       for(m=0; m<numOf2DWlzFiles; m++)
       {
          /* printf("%s  %s\n", PureBibFilesNameList[m], SampleBibPureFilesName[i] ); */
          if(strcmp(PureBibFilesNameList[m], SampleBibPureFilesName[i] ) == 0)
          {
	     BibFileIndex[i] = m;
             break;
          }
       }

       printf("m = %d\n",BibFileIndex[i]);
     }


   centorOfMass = 0;
   if(centorOfMass)
   {
      /* read Woolz object */
      if((inFile = fopen(inFileStrw, "r")) == NULL )
      {
        printf("cannot open the input woolz file.\n");
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
    }
    /* read Section viewer parameters */
    
     centorOfMass = 1;

   if(centorOfMass)
   {
      wlzViewStr      = WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, NULL);
      wlzViewStr1     = WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, NULL);
      wlzViewStrInter = WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, NULL);
   }
   /* interpolation */
 
   globalCycle = 0;
   while(globalCycle < numOfSampleBibFiles - 1)
   {
      printf("%d\n",globalCycle );
      errNum = ReadBibFile(inFile, SampleBibFilesNameList[globalCycle],   wlzViewStr  );
      errNum = ReadBibFile(inFile, SampleBibFilesNameList[globalCycle+1], wlzViewStr1 );
      errNum = ReadBibFile(inFile, SampleBibFilesNameList[globalCycle],   wlzViewStrInter  );
      TotalN = BibFileIndex[globalCycle+1] - BibFileIndex[globalCycle];


      /* cycle through between sample bib file */
      j = BibFileIndex[globalCycle] + 1;
      k = BibFileIndex[globalCycle+1];
      printf("NNNN:  %d %d\n",j,k );
      i=j;
      while ( i < k )
      {
          Inumber = i - BibFileIndex[globalCycle];
	  
	  /* linear interpolations */
	  errNum = LinearInterpolations(wlzViewStr, wlzViewStr1, wlzViewStrInter, Inumber, TotalN);

	  /* change the fixed point to (0 0 0)  */
          errNum = ChangeFixedPoint(wlzViewStrInter, 0., 0., 0.);
         
	  /* write Section viewer parameters */
	  printf("%s\n", PureBibFilesNameList[i]);
	  errNum = WriteBibFile(outFile, PureBibFilesNameList[i], inFileStrw,  TwoDImageFilesNameList[i], wlzViewStrInter);
	  i++;
      }
      globalCycle++;
   }

    /* Now cover the orginal input bib files with new  */
 
   globalCycle = 0;
   while(globalCycle < numOfSampleBibFiles)
   {
     i = BibFileIndex[globalCycle];
      errNum = ReadBibFile(inFile, SampleBibFilesNameList[globalCycle],   wlzViewStr  );
      errNum = ReadBibFile(inFile, SampleBibFilesNameList[globalCycle],   wlzViewStr1 );
      errNum = ReadBibFile(inFile, SampleBibFilesNameList[globalCycle],   wlzViewStrInter  );
      errNum = LinearInterpolations(wlzViewStr, wlzViewStr, wlzViewStrInter, 1, 1);
      /* change the fixed point to (0 0 0)  */
      errNum = ChangeFixedPoint(wlzViewStrInter, 0., 0., 0.);

      errNum = WriteBibFile(outFile, PureBibFilesNameList[i], inFileStrw,  TwoDImageFilesNameList[i], wlzViewStrInter);
      globalCycle++;
   } 

    /*  */
    centorOfMass = 0;

   if(centorOfMass)
   {
       WlzFreeObj(WObjS);
   }
   centorOfMass = 1;

   if(centorOfMass)
   {
       WlzFree3DViewStruct(wlzViewStr);
       WlzFree3DViewStruct(wlzViewStr1);
       WlzFree3DViewStruct(wlzViewStrInter);
   }  
     

  return ( 0 );
}

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s [-h] [<input file>]\n"
	  "\tAutomatically interpolate bib files.\n"
	  "\n"
	  "\tThe user needs to prepare some sample bib files by using MAPaint.\n"
	  "\tAnd the user should also prepare some file name list. Please give\n"
	  "\tfull directory including the machine name when give your file name.\n"
	  "\n"
	  "\tFor example:\n"
	  "\t/net/laphroaig/export/data1/rjianguo/newcastle/LizpDir/\n"
	  "\tsampleBibFilesDir/s11r2s06_template_g.bib\n"
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
          "\t  -w        input file name for the Woolz Object( OPT model !)\n"
	  "\t  -d        input the Directory of the output bib files\n"
	  "\t                                                    \n"
	  "\t  -L        input name of a file which contains the list\n"
	  "\t                     of the 2D raw data image file names\n"
	  "\t                                                        \n"
	  "\t  -B        input name of a file which contains the Sample\n"
	  "\t                     Bib file names \n"
	  "\t                                                        \n"
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
static WlzErrorNum ReadBibFile(FILE *inFile, char *inFileStr, WlzThreeDViewStruct  *wlzViewStr )
{
  WlzErrorNum	    errNum = WLZ_ERR_NONE;
  BibFileRecord	   *bibfileRecord;
  BibFileError      bibFileErr;
  char             *errMsg;
   if((inFile = fopen(inFileStr, "r")) == NULL )
   {
       printf("cannot open the input sample bib  file.\n");
       exit(1);
   }
   /* read the bibfile */
   bibFileErr = BibFileRecordRead(&bibfileRecord, &errMsg, inFile);
   while( bibFileErr == BIBFILE_ER_NONE ){
     if( !strncmp(bibfileRecord->name, "Wlz3DSectionViewParams", 22) ){
       /* allocate memory for wlzViewStr */

       WlzEffBibParse3DSectionViewParamsRecord(bibfileRecord, wlzViewStr);

     }
     BibFileRecordFree(&bibfileRecord);
     bibFileErr = BibFileRecordRead(&bibfileRecord, &errMsg, inFile);
   }
     BibFileRecordFree(&bibfileRecord);
     fclose(inFile);
     inFile = NULL;        
     return errNum;
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
static WlzErrorNum WriteBibFile(FILE *outFile, char *outFileStr,  char *inFileStrw, char *inFileStrSec,WlzThreeDViewStruct  *wlzViewStr)
{
   WlzErrorNum	    errNum = WLZ_ERR_NONE;
   BibFileRecord	*bibfileRecord;
   /* BibFileError      bibFileErr; */
   /* char             *errMsg; */
   char             *tmpS,
                    *dateS = NULL,
		    *hostS = NULL,
                    *userS = NULL,
                    *refFileS = NULL,
     *srcFileS; /*,     sgnlFileS; */
   time_t		tmpTime;   
   char			tmpBuf[256];
   static char		unknownS[] = "unknown";
   /* WlzBasisFnType     bt = WLZ_BASISFN_MQ; */
   WlzFnType          bt = WLZ_FN_BASIS_3DMQ;

   WlzTransformType   at = WLZ_TRANSFORM_2D_NOSHEAR;
   WlzMeshGenMethod     wmgmd  = WLZ_MESH_GENMETHOD_GRADIENT;
   WlzEffFormat    weft =  WLZEFF_FORMAT_WLZ;
   /* WlzObjectType    wtp =    WLZ_2D_DOMAINOBJ; */
   if((outFile = fopen(outFileStr, "w")) == NULL )
   {
       printf("cannot open the output bib file.\n");
       exit(1);
   }
   /* write some sort of identifier */
      bibfileRecord = 
	BibFileRecordMake("Ident", "0",
			  BibFileFieldMakeVa("Text",
					     "MAPaint 2D warp input parameters",
					     "Version",	"1",
					     NULL));
      BibFileRecordWrite(outFile, NULL, bibfileRecord);
      BibFileRecordFree(&bibfileRecord);


   /* now a comment with user, machine, date etc. */
      tmpS = getenv("USER");
      (void )sprintf(tmpBuf, "User: %s", tmpS?tmpS:unknownS);
      userS = AlcStrDup(tmpBuf);

      tmpTime = time(NULL);
      tmpS = ctime(&tmpTime);
      *(tmpS + strlen(tmpS) - 1) = '\0';
      (void )sprintf(tmpBuf, "Date: %s", tmpS?tmpS:unknownS);
      dateS = AlcStrDup(tmpBuf);

      tmpS = getenv("HOST");
      (void )sprintf(tmpBuf, "Host: %s", tmpS?tmpS:unknownS);
      hostS = AlcStrDup(tmpBuf);

      (void )sprintf(tmpBuf, "RefFile: %s", "WlzFileName");
      refFileS = AlcStrDup(tmpBuf);

      (void )sprintf(tmpBuf, "SrcFile: %s",outFileStr );
      srcFileS = AlcStrDup(tmpBuf);

      bibfileRecord = 
	BibFileRecordMake("Comment", "0",
			  BibFileFieldMakeVa("Text", userS,
					     "Text", dateS,
					     "Text", hostS,
					     "Text", refFileS,
					     "Text", srcFileS,
					     "Text", "no need",
					     NULL));
      BibFileRecordWrite(outFile, NULL, bibfileRecord);
      BibFileRecordFree(&bibfileRecord);

       /* if write a file record for the reference file */
      if( inFileStrw ){
	WlzEffBibWriteFileRecord(outFile, "MAPaintReferenceFile",
			 inFileStrw,
			  weft);
      }

      /* if defined write a file record for the source */
      if( inFileStrSec ){
	WlzEffBibWriteFileRecord(outFile, "MAPaintWarpInputSourceFile",
			  inFileStrSec,
			  weft);
      }

      /* write the section data */
      if( WlzEffBibWrite3DSectionViewParamsRecord(outFile, "Wlz3DSectionViewParams", 
					      wlzViewStr) != WLZ_ERR_NONE ){
	printf(
			"Save Warp Parameters:\n"
			"    Error in writing the bibfile\n"
			"    Please check disk space or quotas\n"
			"    Section parameters not saved"
			);
      }

      /* write the warp transform parameters */
       if( WlzEffBibWriteWarpTransformParamsRecord(outFile, "WlzWarpTransformParams",
					      bt,
					      at,
					      wmgmd,
					      20,
					      40)
     	 != WLZ_ERR_NONE ){
     	printf(
			"Save Warp Parameters:\n"
			"    Error in writing the bibfile\n"
			"    Please check disk space or quotas\n"
			"    Section parameters not saved"
			);
      }

    /* write the warp transform parameters */
   /*
   if( WlzEffBibWrite3DSectionViewParamsRecord(outFile, outFileStr, 
					      wlzViewStr) != WLZ_ERR_NONE ) {
		printf("can not output:");	
		exit(0);
					      
   }
   */
   fclose(outFile);
   outFile = NULL;
   return errNum;
}

WlzErrorNum     LinearInterpolations(  WlzThreeDViewStruct *wlzViewStr, WlzThreeDViewStruct  *wlzViewStr1, 
                                       WlzThreeDViewStruct *wlzViewStrInter, int Inumber, int TotalN )
{
  /* int i,j,k; */
   double ratio;
   WlzErrorNum	    errNum = WLZ_ERR_NONE;
   ratio = ( (double) Inumber ) /(double)(TotalN);
   /* linear interpolate  */
   wlzViewStrInter->theta =       wlzViewStr->theta 
                             +  ( wlzViewStr1->theta - wlzViewStr->theta) * ratio;
 	     
   wlzViewStrInter->phi   =       wlzViewStr->phi 
                             +  ( wlzViewStr1->phi   - wlzViewStr->phi)   * ratio;
    
   wlzViewStrInter->zeta   =      wlzViewStr->zeta 
                             +  ( wlzViewStr1->zeta   - wlzViewStr->zeta) * ratio;
   
   wlzViewStrInter->dist   =      wlzViewStr->dist 
                             +  ( wlzViewStr1->dist  - wlzViewStr->dist)  * ratio;
   
   wlzViewStrInter->scale  =      wlzViewStr->scale 
                             +  ( wlzViewStr1->scale - wlzViewStr->scale) * ratio;
   
   wlzViewStrInter->fixed.vtX  =  wlzViewStr->fixed.vtX  
                             +  ( wlzViewStr1->fixed.vtX - wlzViewStr->fixed.vtX ) * ratio;
   
   wlzViewStrInter->fixed.vtY  =  wlzViewStr->fixed.vtY  
                             +  ( wlzViewStr1->fixed.vtY - wlzViewStr->fixed.vtY ) * ratio;
    
   wlzViewStrInter->fixed.vtZ  =  wlzViewStr->fixed.vtZ  
                             +  ( wlzViewStr1->fixed.vtZ - wlzViewStr->fixed.vtZ ) * ratio;

   /* the following is not needed, if don't want to get the same fixed point */			     

   /*  change the fixed point to (0 0 0) which will change the distance parameter */
   /*
   Wlz3DViewGetPlaneEqn(wlzViewStrInter, &dstA, &dstB, &dstC, &dstD);
   */
   /* get the distance */
   /*
   wlzViewStrInter->dist = dstD/sqrt( dstA * dstA + dstB * dstB + dstC * dstC );
   wlzViewStrInter->fixed.vtX = 0.;
   wlzViewStrInter->fixed.vtY = 0.;
   wlzViewStrInter->fixed.vtZ = 0.;
   */
   return errNum;

}

WlzErrorNum     ChangeFixedPoint(WlzThreeDViewStruct *wlzViewStr, double xa, double ya, double za )
{
   WlzErrorNum	  errNum = WLZ_ERR_NONE;
   double         dstA, dstB, dstC, dstD;
   double         de, ci;

   /*  change the fixed point to (0 0 0) which will change the distance parameter */
   Wlz3DViewGetPlaneEqn(wlzViewStr, &dstA, &dstB, &dstC, &dstD);
   de =   dstA * xa + dstB * ya + dstC * za + dstD;
   ci =   sqrt( dstA * dstA + dstB * dstB + dstC * dstC );
   /* get the distance */
   wlzViewStr->dist = -de/ci;
   /*
   wlzViewStr->dist =  dstD /sqrt( dstA * dstA + dstB * dstB + dstC * dstC );
   */
   wlzViewStr->fixed.vtX  = xa;
   wlzViewStr->fixed.vtY  = ya;
   wlzViewStr->fixed.vtZ  = za;

   /* to get 3D sections you need the type  */
   wlzViewStr->type = WLZ_3D_VIEW_STRUCT;

   return errNum;
}

