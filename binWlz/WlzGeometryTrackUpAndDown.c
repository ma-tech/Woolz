#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzGeometryTrackUpAndDown_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzGeometryTrackUpAndDown.c
* \author       Jianguo Rao
* \date         January 2003
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
* \brief	Interpolates tie-point bibfiles.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzgeometrytrackupanddown "WlzGeometryTrackUpAndDown"
*/

/*!
\ingroup BinWlz
\defgroup wlzgeometrytrackupanddown WlzGeometryTrackUpAndDown
\par Name
WlzGeometryTrackUpAndDown - interpolates tie-point bibfiles.
\par Synopsis
\verbatim
WlzGeometryTrackUpAndDown [-h] [-s #] [-t #] [-n #] [-o #] [-O #] [-i #]
                          [-l #] [-j #] [-J #] [-a #] [-A #] [-b #] [-B #]
			  [-c #] [-d #] [-e #] [-E #] [-f #] [-F #] [-g #]
			  [-G #] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-L</b></td>
    <td>
    Input file which contains the list of 2D Woolz contour image
    file names.
    </td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>
    Output file for surface points.
    </td>
  </tr>
  <tr> 
    <td><b>-t</b></td>
    <td>
    Input number to decide track down or up, with 0 - track down (the default)
    and 1 - track up.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file to show the shells and loops.</td>
  </tr>
  <tr> 
    <td><b>-O</b></td>
    <td>Output file for surface points and their distance to the
    surface.</td>
  </tr>
  <tr> 
    <td><b>-n</b></td>
    <td>Number of pixels for z-direction of vtk contour.</td>
  </tr>
  <tr> 
    <td><b>-d</b></td>
    <td>Maximum distance in number of pixels allowed to be tracked,
        default 10.</td>
  </tr>
  <tr> 
    <td><b>-f</b></td>
    <td>Distance in number of pixels for in and out points, default 15.</td>
  </tr>
  <tr> 
    <td><b>-F</b></td>
    <td>Distance in number of pixels for in and out points guid,
        default is 15.</td>
  </tr>
  <tr> 
    <td><b>-g</b></td>
    <td>Number of files to tack up or down, default is 5.</td>
  </tr>
  <tr> 
    <td><b>-G</b></td>
    <td>Start tracking file number (counting from  0), default 5.</td>
  </tr>
  <tr> 
    <td><b>-j</b></td>
    <td>Length to be used to segment the line in unit of number of pixels,
        default 40.</td>
  </tr>
  <tr> 
    <td><b>-C</b></td>
    <td>Length to be used to sample the data in unit of number,
        default 30.</td>
  </tr>
  <tr> 
    <td><b>-l</b></td>
    <td>Number of lth loop for output.</td>
  </tr>
  <tr> 
    <td><b>-a</b></td>
    <td>Start shell number (not index) for tracking.</td>
  </tr>
  <tr> 
    <td><b>-A</b></td>
    <td>End shell number (not index) for tracking.</td>
  </tr>
  <tr> 
    <td><b>-b</b></td>
    <td>Start section number for tracking.</td>
  </tr>
  <tr> 
    <td><b>-B</b></td>
    <td>End section number for tracking.</td>
  </tr>
</table>
\par Description
Automatically interpolates tie point bibfiles.
The user needs to prepare some sample 2D Woolz File names which
contains a number of files with which the middle one will be used as a
standard contour and a start plan to track up and down.
The base part of sample bib file's name should be the same as
that of 2D Woolz file.
\par Examples
\verbatim
\endverbatim
\par File
\ref WlzGeometryTrackUpAndDown.c "WlzGeometryTrackUpAndDown.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <float.h>
#include <string.h>
#include <time.h>
#include <Wlz.h>

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

#ifdef WLZ_UNUSED_FUNCTIONS
static void GetGreyValue( WlzGreyValueWSpace *gVWSp, int kz, int jy, int ix, 
                     int *intensity,    WlzErrorNum  *dstErr );
		     
static void FillGreyValue( WlzGreyValueWSpace *gVWSp, int kz, int jy, int ix, 
                     int intensity,    WlzErrorNum  *dstErr );


static void  OutPutSectionsForTest(
                  WlzObject   *WObjS, 
		  WlzObject   *WObj2D, 
		  FILE        *outFile, 
		  char        *outFileStr,
		  int          iloop,
		  WlzErrorNum *errNum);
#endif /* WLZ_UNUSED_FUNCTIONS */


static void usage(char *proc_str);  


int main(int	argc,
	 char	**argv)
{

  FILE	        *inFile = NULL;   /* used to read Woolz object */ 
  /* FILE       *outFile = NULL;   used to read Woolz object */
  int           i;
  int           startShell = -10, endShell = -15, startSection = -10, endSection = -15;
  int           startLoop  = -10, endLoop  = -15;
  int           downOrUp = 0;
  int		option;
  int           numP_in_Z= 10;
  int           numOf2DWlzFiles = 0;
  int           sectionLength_N       = 40;
  int           subSubSectionLength_L = sectionLength_N /2;
  int           numberOfSampleP_k     =  2;
  int           numOfTrackUpOrDown    =  5;
  int           startTrackingFile     =  5;
  int           iloop=1, test=0;
  double        DisForInOut = 15, DisForInOutGuid = 15;
  double        minDis = 10.;
  /* char       *outFileStr; unused */
  unsigned char        **TwoDImageFilesNameList;
  char         *ContourFilesNameList = NULL;
  char         *surfacePointFileName = NULL, *surfaceInPointFileName = NULL,
               *surfaceOutPointFileName = NULL;
  /* WlzDVertex3          *SurfacePatchPoints; unused */
  WlzErrorNum	       errNum = WLZ_ERR_NONE;
  AlcErrno             alcErr = ALC_ER_NONE;

  /* read the argument list and check for an input file */
  static char	optList[] = "s:t:n:o:O:i:l:L:j:J:a:A:b:B:c:d:e:E:f:F:g:G:h";

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
        case 't':
	 if(sscanf(optarg, "%d", &downOrUp) != 1)
	     {
	       printf("read error");
	       exit(1);
	     }
	     break;
        case 'd':
	 if(sscanf(optarg, "%lg", &minDis) != 1)
	     {
	       printf("read error");
	       exit(1);
	     }
	     break;
       case 'f':
	 if(sscanf(optarg, "%lg", &DisForInOut) != 1)
	     {
	       printf("read error");
	       exit(1);
	     }
	     break;
       case 'F':
	 if(sscanf(optarg, "%lg", &DisForInOutGuid) != 1)
	     {
	       printf("read error");
	       exit(1);
	     }
	     break;
       case 'g':
	 if(sscanf(optarg, "%d", &numOfTrackUpOrDown) != 1)
	     {
	       printf("read error");
	       exit(1);
	     }
	     break;
       case 'G':
	 if(sscanf(optarg, "%d", &startTrackingFile) != 1)
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
       case 'c':
	 if(sscanf(optarg, "%d", &subSubSectionLength_L) != 1)
	     {
	       printf("read error");
	       exit(1);
	     }
	     break;
       case 'o':
	     /* outFileStr = optarg; unused */
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
       case 'a':
	 if(sscanf(optarg, "%d", &startShell) != 1)
	     {
	       printf("read error");
	       exit(1);
	     }
	     break;
       case 'A':
	 if(sscanf(optarg, "%d", &endShell) != 1)
	     {
	       printf("read error");
	       exit(1);
	     }
	     break;
       case 'b':
	 if(sscanf(optarg, "%d", &startSection) != 1)
	     {
	       printf("read error");
	       exit(1);
	     }
	     break;
       case 'B':
	 if(sscanf(optarg, "%d", &endSection) != 1)
	     {
	       printf("read error");
	       exit(1);
	     }
       case 'e':
	 if(sscanf(optarg, "%d", &startLoop) != 1)
	     {
	       printf("read error");
	       exit(1);
	     }
	     break;
       case 'E':
	 if(sscanf(optarg, "%d", &endLoop) != 1)
	     {
	       printf("read error");
	       exit(1);
	     }
	     break;
	     
        default:
             return(0);
      }
  }

  /*  allocate memory for the file name list */
  alcErr = AlcUnchar2Calloc(&TwoDImageFilesNameList, 50, 120);
  if(alcErr != ALC_ER_NONE)
  {
     printf("allocate memory for 2D array of file list");
     exit(1);

  }
     
  /* Read the contour Wlz file name List */
  if((inFile = fopen(ContourFilesNameList, "r")) == NULL )
  {
        printf("cannot open the 2D image files list file.\n");
        exit(1);
  }
  i=0;
  while(!feof(inFile))
  {
        fscanf(inFile, "%s", *(TwoDImageFilesNameList + i ));
	/*  printf("%s\n", *(TwoDImageFilesNameList + i ) ); */
	i++;
  }
  fclose(inFile); 
  inFile = NULL;
  numOf2DWlzFiles = i-1;

  if( 2 * numOfTrackUpOrDown + 1 > numOf2DWlzFiles)
  {
      printf("there are not enough 2D image files to track.\n");
      exit(0);

  }

  if( ( startTrackingFile - numOfTrackUpOrDown >0 )  ||
      ( startTrackingFile + numOfTrackUpOrDown) >  numOf2DWlzFiles )
      {
        printf("there are not enough 2D image files to track if we use the file as a start.\n");
        exit(0);
      }
    /* ---- output for test ---- */
    /*
    if(test  &&  ( errNum == WLZ_ERR_NONE ) )
    {
      OutPutSectionsForTest(WObjS, WObj2D, outFile, outFileStr, iloop, &errNum);
     }
    */
    /*---------- extract contour and track down or up  to get the surface patch ----------- */
    if( (!test)  && ( errNum == WLZ_ERR_NONE ) )
    {
      /* SurfacePatchPoints = unused */
      (void )WlzGeometryTrackUpAndDown_s(      
							numP_in_Z,
                                                        startTrackingFile,
							numOfTrackUpOrDown,
							DisForInOutGuid,
                                                        DisForInOut,
                                                        TwoDImageFilesNameList,
                                                        numOf2DWlzFiles,
                                                        downOrUp,  
  							sectionLength_N,
  							subSubSectionLength_L,
  							numberOfSampleP_k,
							surfacePointFileName,
							surfaceInPointFileName,
							surfaceOutPointFileName, 
							startShell,
							endShell,
							startSection,
							endSection,
							minDis,
						       &errNum );
    }

   Alc2Free( (void **)TwoDImageFilesNameList);

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
	  "Version: %s\n"
	  "Options:\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -L        input name of a file which contains the list of \n" 
          "\t                  of the 2D Woolz contour image file names  \n"
	  "\t                     \n"
	  "\t  -s        input name of a data file which used to output surface points\n"
	  "\t                     \n"
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
	  "\t  -d        input maximum distance in number of pixels allowed to be tracked \n"
	  "\t                     default is 10                                   \n"
	  "\t  -f        input  distance in number of pixels for in and out points \n"
	  "\t                     default is 15 \n "   
	  "\t  -F        input  distance in number of pixels for in and out points  guid\n"
	  "\t                     default is 15  \n"   
	  "\t  -g        input  number of files to tack up or down\n"
	  "\t                     default is 5  \n"   
	  "\t  -G        input  start tracking file number (counting from  0)\n"
	  "\t                     default is 5  \n"   
	  "\t  -j        input number of pixels to be sampled in the length of segmented line see J option\n"
	  "\t                     default is 2                                   \n"
	  "\t  -J        input length to be used to segment the line in unit of number of pixels\n"
	  "\t                     default is 40                                   \n"
	  "\t  -c        input length to be used to sample the data in unit of number of pixels\n"
	  "\t                     default is 30                                   \n"
	  "\t  -l        input number of lth loop for output \n"
	  "\t  -a        input the start shell number (not index ) for tracking\n"
	  "\t  -A        input the end   shell number (not index ) for tracking\n"
	  "\t  -b        input the start section number for tracking\n"
	  "\t  -B        input the end   section number for tracking\n"
	  "\t                  \n",
	  proc_str,
	  WlzVersion());
  return;
}


#ifdef WLZ_UNUSED_FUNCTIONS
/*!
* - Function:   WlzEffWriteMeshTransform3DWithoutDisplacementVTK
* - Returns:    none 
* - Purpose:    output the orginal mesh.	
* - Parameters:	
*     -# *fp:               pointer pointing to a specific file.
*     -#  wmt3D:            mesh transform.
* - Author:       J. Rao, R. Baldock and B. Hill
*/
static void  OutPutSectionsForTest(
                  WlzObject   *WObjS, 
		  WlzObject   *WObj2D, 
		  /* FILE        *outFile;  */
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
	      /*  GetGreyValue(gVWSp, 0, jy, ix, &intensity, &errNum); */
	      intensity = 255;
	      FillGreyValue(gVWSp, 0, jy, ix, intensity, &wErrN);
	      
          } 
       }

   }

  /*  visualize the data: */
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
		/* if(kt%3 == iloop) */
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
#endif /* WLZ_UNUSED_FUNCTIONS */

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
