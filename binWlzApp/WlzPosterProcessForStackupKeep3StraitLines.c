#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzPosterProcessForStackupKeep3StraitLines_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/WlzPosterProcessForStackupKeep3StraitLines.c
* \author       Jianguo Rao
* \date         April 2004
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
* \brief	Compute transforms to stack sections keeping verticals.
* \ingroup	BinWlzApp
*
* \par Binary
* \ref wlzposterprocessforstackupkeepvertical "WlzPosterProcessForStackupKeepVertical"
*/

/*!
\ingroup BinWlzApp
\defgroup wlzposterprocessforstackupkeepvertical WlzPosterProcessForStackupKeepVertical
\par Name
WlzPosterProcessForStackupKeepVertical
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>
#include <float.h>
#include <string.h>
#include <time.h>
#include <Wlz.h>

#include <WlzExtFF.h>
#include <bibFile.h>


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

static WlzErrorNum WlzTrack3DVertical( WlzVertex sip, 
                                WlzVertex *sip1p,
                                WlzThreeDViewStruct  *wlzViewStri,  
			        WlzThreeDViewStruct  *wlzViewStrip1,
				int UpOrDown);
				
static WlzErrorNum WlzGetCorssPoint( WlzVertex sip,
                                WlzVertex nStraghtline,
                                WlzVertex *sip1p,
                                WlzThreeDViewStruct  *wlzViewStri  
				);


static WlzErrorNum WlzGetDxDy( 
                                WlzVertex nStraghtline,
				double    dz, 
                                WlzVertex *sip1p
				);

static WlzErrorNum outputTheWarpedWlz( WlzVertex  vsp, 
                                       WlzVertex  vtp,
				       WlzVertex  vfp,
                                       WlzDVertex2 *vtxVec1,
                                       char *souce2DWlzStr,
				       char *out2DWlzStr
				     );

int main(int	argc,
	 char	**argv)
{

  FILE	       *inFile = NULL;   /* used to read Woolz object */ 
  int           i,j,k=0, m;
  int		option;
  int           startJ;
  int           cycleI, cycleNext;
  double           ntnm,nfnm,nsnm;
  int           numOf2DWlzFiles = 0, numOfBibFiles = 0;
  int           ok, usage1, idx;
  double        deltaZ;
  double        sx=0,sy=0,sz=0, fx=0,fy=100,fz=0, tx=100,ty=0,tz=0;
  double        Sx=0,Sy=0,Sz=60, Fx=0,Fy=100,Fz=60, Tx=100,Ty=0,Tz=60;
  double        lengthTS;
  char         *inFileStr, *inFileStr1, *outFileStr;
  char         *inFileStrw=NULL, *outputWlzFilesDir=NULL;
  char          TwoDImageFilesNameList[700][120];
  char          TwoDImagePureFilesNameList[700][120];
  char          warpedFilesNameList[700][120];
  char          BibFilesNameList[700][120];

  char          TwoDImageFilesDir[120];
  char         *List2DImageFilesStr= NULL, *ListBibFilesStr= NULL;


  WlzVertex      vs, vt, vf,  vS, vT, vF;
  WlzVertex      vs1, vt1, vf1, vsp, vtp, vfp, vST;

  WlzVertex      vsr, vtr, vfr;
  WlzVertex      ns, nf, nt, nsp, nfp, ntp, nz,nm;
  WlzDVertex2   *vtxVec1;
  WlzErrorNum	 errNum = WLZ_ERR_NONE;
  WlzThreeDViewStruct *wlzViewStr, *wlzViewStr1;
  WlzObject     *WObjS; 
  char		*cutStr[3];
  double        cutVal[3];

  /* read the argument list and check for an input file */
  static char	optList[] = "i:I:w:o:d:s:S:f:F:t:T:M:n:r:R:L:B:h",opterr = 0;
  ok = 1;
 
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
	    ListBibFilesStr = optarg;
	    break;
        case 'w':
	    inFileStrw = optarg;
	    break;
        case 'n':
	 if(sscanf(optarg, "%d", &startJ) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;
        case 's':
        case 'f':
        case 't':
        case 'S':
        case 'F':
        case 'T':
	  if(optarg)
	  {
	  	while(*optarg && isspace(*optarg))
	  	{
                    ++optarg;
		}
	 	if(*optarg == ',')
	  	{
	    		cutStr[0] = NULL;
	    		cutStr[1] = strtok(optarg, ",");
	    		cutStr[2] = strtok(optarg, ",");
	  	}
	  	else
	  	{
	    		cutStr[0] = strtok(optarg, ",");
	    		cutStr[1] = strtok(NULL, ",");
		    	cutStr[2] = strtok(NULL, ",");
	  	}
	  	if((cutStr[0] == NULL) && (cutStr[1] == NULL) && (cutStr[2] == NULL)   )
	  	{
	    		usage1 = 1;
	    		ok = 0;
	  	}
	  	else
	  	{
	    		idx = 0;
	    		while(ok && (idx < 3))
	    		{
	      			if(cutStr[idx] && (sscanf(cutStr[idx], "%lg",
					 cutVal + idx) != 1))
	      			{
					usage1 = 1;
					ok = 0;
	      			}
	      			++idx;
	    		}
	  	}
	}
	if(ok)
	{
	  switch(option)
	  {
	    case 's':
		sx = (double) cutVal[0];
		sy = (double) cutVal[1];
		sz = (double) cutVal[2];
		break;
	    case 'f':
		fx = (double) cutVal[0];
		fy = (double) cutVal[1];
		fz = (double) cutVal[2];
		break;
	    case 't':
		tx = (double) cutVal[0];
		ty = (double) cutVal[1];
		tz = (double) cutVal[2];
		break;
	    case 'S':
		Sx = (double) cutVal[0];
		Sy = (double) cutVal[1];
		Sz = (double) cutVal[2];
		break;
	    case 'F':
		Fx = (double) cutVal[0];
		Fy = (double) cutVal[1];
		Fz = (double) cutVal[2];
		break;
	    case 'T':
		Tx = (double) cutVal[0];
		Ty = (double) cutVal[1];
		Tz = (double) cutVal[2];
		break;
	  }	

	}
	break;
        case 'o':
	    outFileStr = optarg;
	    break;
        case 'd':
	    outputWlzFilesDir = optarg;
	    break;
        default:
              return(0);
      }
   }

      printf("t: %f %f %f\n",tx,ty, tz );
      printf("T: %f %f %f\n",Tx,Ty, Tz );

    /*------- allocate memory --------*/
    if (
         ((vtxVec1 = (WlzDVertex2 *)AlcCalloc(sizeof(WlzDVertex2), 3))  == NULL) 
      )
    {
       errNum = WLZ_ERR_MEM_ALLOC;
    }

    /* ------Read 2D Wlz files List------- */
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
          if( TwoDImageFilesNameList[m][i] == (char) NULL )
           break; 
          if( TwoDImageFilesNameList[m][i] == '/' )
            k = i; 	   
       }
       if(k != 0)
       k++;
       for(j=k; j<i; j++)
       {
         TwoDImagePureFilesNameList[m][j-k] = TwoDImageFilesNameList[m][j];
         warpedFilesNameList[m][j-k]        = TwoDImageFilesNameList[m][j];
       }
       TwoDImagePureFilesNameList[m][j-k] = ( char ) NULL;
       warpedFilesNameList[m][j-k]        = ( char ) NULL;
         
       /* extract directory */
       if(m == 0)
       {
          for(j=0; j<k; j++)
	  {
	     TwoDImageFilesDir[j] = TwoDImageFilesNameList[m][j];
	  }
	     TwoDImageFilesDir[j] = ( char ) NULL;
	     printf("I_DIR  %s\n",TwoDImageFilesDir );
       }
       printf("%s\n", TwoDImagePureFilesNameList[m]);   
    }

/*--- Read bib File name list ----*/
    if((inFile = fopen(ListBibFilesStr, "r")) == NULL )
    {
        printf("cannot open the 2D image files list file.\n");
        exit(1);
    }
    i=0;
    while(!feof(inFile)){
        fscanf(inFile, "%s", BibFilesNameList[i]);
	printf("%s\n", BibFilesNameList[i]);
	i++;
    }
      fclose(inFile); 
      inFile = NULL;
    numOfBibFiles = i-1;
    printf("number of files bib files = %d\n", numOfBibFiles );

/* check whether the number of bib files is the same as number of Wlz 2D sections */
    if(numOfBibFiles != numOf2DWlzFiles )
    {
         printf("Number of 2D woolz files and number of bib files are different\n");
	 exit(1);
    }

/* get the output file namelist */
   i =  (int) strlen( (const char *) outputWlzFilesDir);
   printf("%d\n", i);
   for( j = 0; j < numOfBibFiles; j++)
   {
      strcpy( warpedFilesNameList[j], (const char *) outputWlzFilesDir );
      strcat( warpedFilesNameList[j], "/");
      strcat(warpedFilesNameList[j], (const char *) TwoDImagePureFilesNameList[j]);
      printf("%s\n", warpedFilesNameList[j] );   
   }   
   
/* --------read Woolz object--------- */
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

/* ------get the fixed S , T and F ready and s, t, f  in OPT space------ */

      vS.d3.vtX = Sx; 
      vS.d3.vtY = Sy;
      vS.d3.vtZ = Sz;
      vT.d3.vtX = Tx; 
      vT.d3.vtY = Ty;
      vT.d3.vtZ = Tz;
      vF.d3.vtX = Fx; 
      vF.d3.vtY = Fy;
      vF.d3.vtZ = Fz;
      
      vs.d3.vtX = sx; 
      vs.d3.vtY = sy;
      vs.d3.vtZ = sz;
      vt.d3.vtX = tx; 
      vt.d3.vtY = ty;
      vt.d3.vtZ = tz;
      vf.d3.vtX = fx; 
      vf.d3.vtY = fy;
      vf.d3.vtZ = fz;
      
      printf("vt: %f %f %f\n",vt.d3.vtX,vt.d3.vtY, vt.d3.vtZ );
      printf("vT: %f %f %f\n",vT.d3.vtX,vT.d3.vtY, vT.d3.vtZ );
/*  -----get the three directions in OPT space------- 
    n = 1/(|r2-r1|) { (x2-x1) i  + (y2-y1) j + (z2-z1) k  }

  */
   
      WLZ_VTX_3_SUB(vST.d3,vS.d3,vs.d3);
      lengthTS  =  WLZ_VTX_3_LENGTH( vST.d3);
      ns.d3.vtX = vST.d3.vtX/ lengthTS;
      ns.d3.vtY = vST.d3.vtY/ lengthTS;
      ns.d3.vtZ = vST.d3.vtZ/ lengthTS;


      WLZ_VTX_3_SUB(vST.d3,vF.d3,vf.d3);
      lengthTS  =  WLZ_VTX_3_LENGTH( vST.d3);
      nf.d3.vtX = vST.d3.vtX/ lengthTS;
      nf.d3.vtY = vST.d3.vtY/ lengthTS;
      nf.d3.vtZ = vST.d3.vtZ/ lengthTS;

      WLZ_VTX_3_SUB(vST.d3,vT.d3,vt.d3);
      printf("vt: %f %f %f\n",vt.d3.vtX,vt.d3.vtY, vt.d3.vtZ );
      printf("vT: %f %f %f\n",vT.d3.vtX,vT.d3.vtY, vT.d3.vtZ );
      lengthTS  =  WLZ_VTX_3_LENGTH( vST.d3);
      printf("Length: %f\n",lengthTS);
      nt.d3.vtX = vST.d3.vtX/ lengthTS;
      nt.d3.vtY = vST.d3.vtY/ lengthTS;
      nt.d3.vtZ = vST.d3.vtZ/ lengthTS;

/* 
  ---- get the three points in the master plane ( startJ plane ) ----

*/

/* allocate memory and give default values for viewer strcuture */
      wlzViewStr      = WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, NULL);
      wlzViewStr1     = WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, NULL);

/* set up */
   cycleI     = startJ;
   /* read Section viewer parameters */
   errNum = ReadBibFile(inFile, (char *) BibFilesNameList[cycleI],   wlzViewStr  );
   if(errNum != WLZ_ERR_NONE)
   {
         printf("Read or parse bib file err.\n");
         exit(1);
   
   }
   inFile = NULL;

   WlzInit3DViewStruct(wlzViewStr, WObjS);


/* -------- get the corresponding three directions in Main view space -------
    (It is also the stack up space 
    nsp = Tm ns, 
    nfp = Tm nf, 
    ntp = Tm nt;
  */
  errNum = Wlz3DSectionTransformVtxR(wlzViewStr, ns.d3, &nsp.d3 );
  if(errNum != WLZ_ERR_NONE)
   {
         printf("can not get the direction in stack up space for ns err.\n");
	 exit(1);
   }
  errNum = Wlz3DSectionTransformVtxR(wlzViewStr, nf.d3, &nfp.d3 );
  if(errNum != WLZ_ERR_NONE)
   {
         printf("can not get the direction in stack up space for nf err.\n");
	 exit(1);
   }
  errNum = Wlz3DSectionTransformVtxR(wlzViewStr, nt.d3, &ntp.d3 );
  if(errNum != WLZ_ERR_NONE)
   {
         printf("can not get the direction in stack up space for np err.\n");
	 exit(1);
   }
     

  /* --- get the main plane direction in OPT space --- */
  /* it can got by inverse transform from stack up viewer space to
     OPT space */
     /* in stack up space or master plane sapce */
      nz.d3.vtX = 0.0;
      nz.d3.vtY = 0.0;
      nz.d3.vtZ = 1.0;
  
  /*  nm = Tm^-1  nz  (inverse transform used here ) nm is in OPT space */
  errNum = Wlz3DSectionTransformInvVtxR(wlzViewStr, nz.d3, &nm.d3 );
  if(errNum != WLZ_ERR_NONE)
   {
         printf("can not get the main plane direction in OPT space err.\n");
         printf("Which means your three straight line has to be rechosen\n");
	 exit(1);
   }
   /* get the nt dot nm, etc */

  ntnm  =  WLZ_VTX_3_DOT(nt.d3, nm.d3);
  nfnm  =  WLZ_VTX_3_DOT(nf.d3, nm.d3);
  nsnm  =  WLZ_VTX_3_DOT(ns.d3, nm.d3);
  ntnm = nm.d3.vtX * nt.d3.vtX +  nm.d3.vtY * nt.d3.vtY + nm.d3.vtZ * nt.d3.vtZ;


  if(  ( ( ntnm == 0 )  || ( nfnm == 0 ) ) ||( nsnm == 0 ) )
    {
      printf("At least one of your three strait lines are parallel to the main plane err.\n");
      printf("Which means your three straight line has to be rechosen or chose another plane\n");
      printf("as main plane\n");
      exit(1);
    }




  /* get crossed point in this plane p_{i} then project to viewer space */
  errNum = WlzGetCorssPoint(vs,  ns,  &vsr,  wlzViewStr);
  if(errNum != WLZ_ERR_NONE)
    {
      printf("can not find crossed point err.\n");
      exit(1);
    }
  errNum = WlzGetCorssPoint(vt,  nt,  &vtr,  wlzViewStr);
  if(errNum != WLZ_ERR_NONE)
    {
      printf("can not find crossed point err.\n");
      exit(1);
    }

  errNum = WlzGetCorssPoint(vf,  nf,  &vfr,  wlzViewStr);
  if(errNum != WLZ_ERR_NONE)
  {
         printf("can not find crossed point err.\n");
         exit(1);
   
  }

  /* use this 3 point as refernce */
  /* --------Get Tie-points -------- */


   /* get the target points for this plane in stack up space */

   (vtxVec1    )->vtX = vsr.d3.vtX;
   (vtxVec1    )->vtY = vsr.d3.vtY;
   (vtxVec1 + 1)->vtX = vtr.d3.vtX;
   (vtxVec1 + 1)->vtY = vtr.d3.vtY;
   (vtxVec1 + 2)->vtX = vfr.d3.vtX;
   (vtxVec1 + 2)->vtY = vfr.d3.vtY;

  /* get the source point */
   /* get ( sx' sy' sz' )_i  */
   printf("%f  %f\n", sx, sy);
   vsp.d3.vtX = vsr.d3.vtX;
   vsp.d3.vtY = vsr.d3.vtY;
   vsp.d3.vtZ = vsr.d3.vtZ;

   /* get ( tx' ty' tz' )_i  */
   printf("%f  %f\n", tx, ty);
   vtp.d3.vtX = vtr.d3.vtX;
   vtp.d3.vtY = vtr.d3.vtY;
   vtp.d3.vtZ = vtr.d3.vtZ;

   /* get ( fx' fy' fz' )_i  */
   printf("%f  %f\n", fx, fy);
   vfp.d3.vtX =  vfr.d3.vtX;
   vfp.d3.vtY =  vfr.d3.vtY;
   vfp.d3.vtZ =  vfr.d3.vtZ;
 
   /* ------ output the affine transformed Wlz ------*/
   errNum =        outputTheWarpedWlz( vsp, 
                                       vtp,
				       vfp,
                                       vtxVec1,
                                       TwoDImageFilesNameList[cycleI],
				       warpedFilesNameList[cycleI]
				     );

   while(cycleI < numOfBibFiles - 1)
   {

      	cycleNext = cycleI + 1;
	
   	errNum = ReadBibFile(inFile, (char *) BibFilesNameList[cycleNext],   wlzViewStr1  );
   	WlzInit3DViewStruct(wlzViewStr1, WObjS);

  	/* ------ get Tie - points -------*/

 	/* get crossed point in this plane p_{i} then project to viewer space */
   	errNum = WlzGetCorssPoint(vs,  ns,  &vs1,  wlzViewStr1);
  	if(errNum != WLZ_ERR_NONE)
   	{
        	printf("can not find crossed point err.\n");
         	exit(1);
   
  	 }
  	errNum = WlzGetCorssPoint(vt,  nt,  &vt1,  wlzViewStr1);
  	if(errNum != WLZ_ERR_NONE)
   	{
        	printf("can not find crossed point err.\n");
         	exit(1);
   
   	}

  	errNum = WlzGetCorssPoint(vf,  nf,  &vf1,  wlzViewStr1);
  	if(errNum != WLZ_ERR_NONE)
   	{
        	printf("can not find crossed point err.\n");
         	exit(1);
   
   	}


  	/*  get the source points */
   	/* get ( sx' sy' sz' )_i  */
   	/* printf("%f  %f\n", sx, sy); */
   	vsp.d3.vtX = vs1.d3.vtX;
   	vsp.d3.vtY = vs1.d3.vtY;
   	vsp.d3.vtZ = vs1.d3.vtZ;

   	/* get ( tx' ty' tz' )_i  */
   	/* printf("%f  %f\n", tx, ty); */
   	vtp.d3.vtX = vt1.d3.vtX;
   	vtp.d3.vtY = vt1.d3.vtY;
   	vtp.d3.vtZ = vt1.d3.vtZ;

   	/* get ( fx' fy' fz' )_i  */
   	/* printf("%f  %f\n", fx, fy); */
   	vfp.d3.vtX =  vf1.d3.vtX;
   	vfp.d3.vtY =  vf1.d3.vtY;
   	vfp.d3.vtZ =  vf1.d3.vtZ;

  	/*  get the target points */

   	/* find the target points by keep the same direction in
	   stack up space and the OPT space.The following are in stack up space */
        deltaZ = (double) ( cycleNext - startJ );


  	errNum = WlzGetDxDy( nsp,  deltaZ,  &vs1);
  	if(errNum != WLZ_ERR_NONE)
   	{
        	printf("can not get Dx Dy err.\n");
         	exit(1);
   
   	}
  	errNum = WlzGetDxDy( nfp,  deltaZ,  &vf1);
  	if(errNum != WLZ_ERR_NONE)
   	{
        	printf("can not get Dx Dy err.\n");
         	exit(1);
   
   	}
  	errNum = WlzGetDxDy( ntp,  deltaZ,  &vt1);
  	if(errNum != WLZ_ERR_NONE)
   	{
        	printf("can not get Dx Dy err.\n");
         	exit(1);
   
   	}

   	(vtxVec1    )->vtX = vsr.d3.vtX + vs1.d3.vtX;
   	(vtxVec1    )->vtY = vsr.d3.vtY + vs1.d3.vtY;
   	(vtxVec1 + 1)->vtX = vtr.d3.vtX + vt1.d3.vtX;
   	(vtxVec1 + 1)->vtY = vtr.d3.vtY + vt1.d3.vtY;
   	(vtxVec1 + 2)->vtX = vfr.d3.vtX + vf1.d3.vtX;
   	(vtxVec1 + 2)->vtY = vfr.d3.vtY + vf1.d3.vtY;
 
   	/* ------ output the affine transformed Wlz ------*/
	printf("----------- %d\n",cycleI);
	/* out put the source and target for checking */
   	errNum =        outputTheWarpedWlz( vsp, 
                                       vtp,
				       vfp,
                                       vtxVec1,
                                       TwoDImageFilesNameList[cycleNext],
				       warpedFilesNameList[cycleNext]
				     );

   	/* track from s_i' to  s_{i+1}' */
   	/*  scale = |T - S |/|t - s|      */
        
        cycleI++;

   }


   cycleI = startJ;
   while(  cycleI > 0  )
   {
      	cycleNext = cycleI - 1;
	
   	errNum = ReadBibFile(inFile, (char *) BibFilesNameList[cycleNext],   wlzViewStr1  );
   	WlzInit3DViewStruct(wlzViewStr1, WObjS);

  	/* ------ get Tie - points -------*/

 	/* get crossed point in this plane p_{i} then project to viewer space */
   	errNum = WlzGetCorssPoint(vs,  ns,  &vs1,  wlzViewStr1);
  	if(errNum != WLZ_ERR_NONE)
   	{
        	printf("can not find crossed point err.\n");
         	exit(1);
   
  	 }
  	errNum = WlzGetCorssPoint(vt,  nt,  &vt1,  wlzViewStr1);
  	if(errNum != WLZ_ERR_NONE)
   	{
        	printf("can not find crossed point err.\n");
         	exit(1);
   
   	}

  	errNum = WlzGetCorssPoint(vf,  nf,  &vf1,  wlzViewStr1);
  	if(errNum != WLZ_ERR_NONE)
   	{
        	printf("can not find crossed point err.\n");
         	exit(1);
   
   	}


  	/*  get the source points */
   	/* get ( sx' sy' sz' )_i  */
   	printf("%f  %f\n", sx, sy);
   	vsp.d3.vtX = vs1.d3.vtX;
   	vsp.d3.vtY = vs1.d3.vtY;
   	vsp.d3.vtZ = vs1.d3.vtZ;

   	/* get ( tx' ty' tz' )_i  */
   	printf("%f  %f\n", tx, ty);
   	vtp.d3.vtX = vt1.d3.vtX;
   	vtp.d3.vtY = vt1.d3.vtY;
   	vtp.d3.vtZ = vt1.d3.vtZ;

   	/* get ( fx' fy' fz' )_i  */
   	printf("%f  %f\n", fx, fy);
   	vfp.d3.vtX =  vf1.d3.vtX;
   	vfp.d3.vtY =  vf1.d3.vtY;
   	vfp.d3.vtZ =  vf1.d3.vtZ;

  	/*  get the target points */

   	/* find the target points by keep the same direction in
	   stack up space and the OPT space.The following are in stack up space */
        deltaZ = (double) ( cycleNext - startJ );


  	errNum = WlzGetDxDy( nsp,  deltaZ,  &vs1);
  	if(errNum != WLZ_ERR_NONE)
   	{
        	printf("can not get Dx Dy err.\n");
         	exit(1);
   
   	}
  	errNum = WlzGetDxDy( nfp,  deltaZ,  &vf1);
  	if(errNum != WLZ_ERR_NONE)
   	{
        	printf("can not get Dx Dy err.\n");
         	exit(1);
   
   	}
  	errNum = WlzGetDxDy( ntp,  deltaZ,  &vt1);
  	if(errNum != WLZ_ERR_NONE)
   	{
        	printf("can not get Dx Dy err.\n");
         	exit(1);
   
   	}


       /* target points */
   	(vtxVec1    )->vtX = vsr.d3.vtX + vs1.d3.vtX;
   	(vtxVec1    )->vtY = vsr.d3.vtY + vs1.d3.vtY;
   	(vtxVec1 + 1)->vtX = vtr.d3.vtX + vt1.d3.vtX;
   	(vtxVec1 + 1)->vtY = vtr.d3.vtY + vt1.d3.vtY;
   	(vtxVec1 + 2)->vtX = vfr.d3.vtX + vf1.d3.vtX;
   	(vtxVec1 + 2)->vtY = vfr.d3.vtY + vf1.d3.vtY;
 
   	/* ------ output the affine transformed Wlz ------*/
   	errNum =   outputTheWarpedWlz( vsp, 
                                       vtp,
				       vfp,
                                       vtxVec1,
                                       TwoDImageFilesNameList[cycleNext],
				       warpedFilesNameList[cycleNext]
				     );

   	/* track from s_i' to  s_{i+1}' */
   	/*  scale = |T - S |/|t - s|      */
        
        cycleI--;

	}



     AlcFree(vtxVec1 );

   /* */
       WlzFree3DViewStruct(wlzViewStr);
       WlzFree3DViewStruct(wlzViewStr1);

  return ( 0 );
}

static void usage(char *proc_str)
{
  (void )fprintf(stderr,
	  "Usage:\t%s [-h] [<input file>]\n"
	  "\tGet the affine transformation to transform the warped or.\n"
	  "\nGutted sections through OPT may be  by differnt angles into\n"
	  "\nThe position we think are right.\n"
	  "\tThe user needs to prepare the bib file file names list which assotiate .\n"
	  "\twith the cuting OPT. and the file name of warped section or the sections\n"
	  "\tcut through OPT. Please give\n"
	  "\tfull directory including the machine name when give your file name.\n"
	  "\tFor example:\n"
	  "\t/net/laphroaig/export/data1/rjianguo/newcastle/LizpDir/\n"
	  "\tBibFilesDir/s11r2s06_template_g.bib\n"
	  "\n"
	  "\tALSO, The base part of sample bib file's name should be the same as\n"
	  "\tthat of 2D Woolz file.\n"
	  "\tFor example:\n"
	  "\t ABC.bib,  ABC.wlz ,\n"
	  "\twhere the base part of their name is the same\n"
	  "\tThe output will be the set of files Affine transformed which will be\n"
	  "\tas  ABC_Affine.wlz\n"
	  "Version: %s\n"
	  "Options:\n"
	  "\t  -h        Help - prints this usage message\n"
          "\t  -w        input file name for the Woolz Object( OPT model !)\n"
	  "\t  -d        input the Directory of the output affineTransformed files\n"
	  "\t                                                    \n"
	  "\t  -L        input name of a file which contains the list\n"
	  "\t                     of the 2D raw data image file names\n"
	  "\t                                                        \n"
	  "\t  -B        input name of a file which contains the\n"
	  "\t                     Bib file names \n"
	  "\t  -n        input j , the j-th layer of the image as a start section\n"
	  "\t                                                    \n"
	  "\t  -s        input the sx,sy,sz      no space please!\n"
	  "\t  -f        input the fx,fy,fz      no space please!\n"
	  "\t  -t        input the tx,ty,tz      no space please!\n"
	  "\t  -S        input the Sx,Sy,Sz      no space please!\n"
	  "\t  -F        input the Fx,Fy,Fz      no space please!\n"
	  "\t  -T        input the Tx,Ty,Tz      no space please!\n",
	  proc_str,
	  WlzVersion());
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
   /* read bibfile records until the filename matches */

      bibFileErr = BibFileRecordRead(&bibfileRecord, &errMsg, inFile);
      while( bibFileErr == BIBFILE_ER_NONE ) 
      {
          if( !strncmp(bibfileRecord->name, "Wlz3DSectionViewParams", 22) ){
             /* parse the record */
             errNum = WlzEffBibParse3DSectionViewParamsRecord(bibfileRecord, wlzViewStr);
             break;
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
static WlzErrorNum WriteBibFile(FILE *outFile, char *outFileStr,  char *inFileStrw, char *inFileStrSec
,WlzThreeDViewStruct  *wlzViewStr)
{
   WlzErrorNum	    errNum = WLZ_ERR_NONE;
   BibFileRecord	   *bibfileRecord;

   char             *tmpS,
                    *dateS = NULL,
		    *hostS = NULL,
                    *userS = NULL,
                    *refFileS = NULL,
		    *srcFileS;
   time_t		tmpTime;   
   char			tmpBuf[256];
   static char		unknownS[] = "unknown";
   WlzFnType          bt = WLZ_FN_BASIS_3DMQ;

   WlzTransformType   at = WLZ_TRANSFORM_2D_NOSHEAR;
   WlzMeshGenMethod     wmgmd  = WLZ_MESH_GENMETHOD_GRADIENT;
   WlzEffFormat    weft =  WLZEFF_FORMAT_WLZ;

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

/*!
* - Function:   WlzTrack3DVertical
* - Returns:    WlzVetex 
* - Purpose:    Track vertex through neighbour layer.	
* - Parameters:	
*
*     -#   sip:              input   WlzVertex.
*     -#   sip1p:            output  WlzVertex by tracking.
*     -#   wlzViewStri:      WlzThreeDViewStruct for layer i;
*     -#   wlzViewStrip1:    WlzThreeDViewStruct for the next layer i+-1;
*     -#   UpOrDown:         > 0 Up;  <= 0 down; 
* - Author:       J. Rao, R. Baldock
*/
static WlzErrorNum WlzTrack3DVertical( WlzVertex sip, 
                                WlzVertex *sip1p,
                                WlzThreeDViewStruct  *wlzViewStri,  
			        WlzThreeDViewStruct  *wlzViewStrip1,
				int UpOrDown)
{

   WlzVertex vtemp, vtemp0, vs, vo, vso, vtx1;
   WlzVertex vs1, nI, nIP1;
   double    dtemp, dtemp1, dtemp2;
   WlzErrorNum	       errNum = WLZ_ERR_NONE;
  

   /* output for test */
   /*
   printf("Yaw:      %f\n",  wlzViewStri->theta*180.0/3.14); 
   printf("Pitch:    %f\n",  wlzViewStri->phi*180.0/3.14 );
   printf("Roll:     %f\n",  wlzViewStri->zeta*180.0/3.14 );
   printf("distance: %f\n",  wlzViewStri->dist );
   printf("scaling:  %f\n",  wlzViewStri->scale );
   printf("%f  %f  %f\n",  wlzViewStri->fixed.vtX, wlzViewStri->fixed.vtY ,wlzViewStri->fixed.vtZ  );
   */

   
   /* vtemp0.d3 = (sx', sy', sz')   */
   WlzValueCopyDVertexToDVertex3(&vtemp0.d3, &sip.d3, 1);
   vtemp0.d3.vtZ = wlzViewStri->dist;
   
   /* reverse affine transformation the vertx to get  vs.d3 = s_i here is in OPT space */
   errNum = Wlz3DSectionTransformInvVtxR(wlzViewStri, vtemp0.d3, &vs.d3 );
   
   printf("%f  %f  %f\n", vs.d3.vtX, vs.d3.vtY, vs.d3.vtZ);

   /* as the vs.d3 is in OPT space and we want vertical tracking
      we have:
      
      vs1.d3.vtX = vs.d3.vtX
      vs1.d3.vtY = vs.d3.vtY

      vs1.d3.vtZ = ? need to be find !!!!!

      vs1.d3.vtZ = vs.d3.vtZ + s

      where s  n_z dot n_{i+1} = [ O_{i+1} - r_{i}     ] dot n_{i+1}

      */
      /*
      nz.d3.vtX = 0.0;
      nz.d3.vtY = 0.0;
      nz.d3.vtZ = 1.0;
      if( UpOrDown <= 0 )
      {
   	nz.d3.vtZ = -1.0;
      }
      */



   vtx1.d3.vtX = 0.0;
   vtx1.d3.vtY = 0.0;
   vtx1.d3.vtZ = 1.0;
   /* get the normal n_i = nI.d3   of this plan i -plane do not need !!! */
   errNum = Wlz3DSectionTransformInvVtxR(wlzViewStri, vtx1.d3, &nI.d3 );
   
   /* get the normal n_i(+-1) = nIP1.d3  of the next plan i+-1 -plane */
   errNum = Wlz3DSectionTransformInvVtxR(wlzViewStrip1, vtx1.d3, &nIP1.d3 );

   /* track the s_{i+-1} vertically in OPT space */

   /* check  the n_{i+1}_z != 0 other wise throw an erro  */
   dtemp =  nIP1.d3.vtZ;
   /* dtemp = abs(dtemp); */
   if( ( dtemp <  0.00000001 ) && ( dtemp >  -0.00000001 ) )
   {
      printf("%lg\n",dtemp);
      printf(" the z-component of n_{i+1} is zero in OPT space, we can't track vertically\n ");
      printf("%lg  %lg  %lg\n", nIP1.d3.vtX, nIP1.d3.vtY, nIP1.d3.vtZ);
      exit( 1 );

   }
   else
   {
     /* not parallel to z-direction  */
     /*  s_{i+1} = S_{i} + n_z * s
                 = S_{i} + n_z * <  [ O_{i+1} - r_{i}     ] dot n_{i+1}  /   ( n_z dot n_{i+1} ) >

     */
    
     	/*  dtemp =    ( n_z dot n_{i+1} )  */
     		dtemp  =  WLZ_VTX_3_DOT(vtx1.d3, nIP1.d3);

     	/* get the  o_i+1 
   		vo.d3.vtX   = wlzViewStr1->fixed.vtX + wlzViewStr->dist * nIP1.d3.vtX;
   		vo.d3.vtY   = wlzViewStr1->fixed.vtY + wlzViewStr->dist * nIP1.d3.vtY;
   		vo.d3.vtY   = wlzViewStr1->fixed.vtZ + wlzViewStr->dist * nIP1.d3.vtZ;
    	*/
    		WLZ_VTX_3_SCALE(vtemp0.d3, nIP1.d3, wlzViewStrip1->dist);
    		WLZ_VTX_3_ADD(vo.d3, wlzViewStrip1->fixed, vtemp0.d3);

    		/* vso.d3 = vo.d3 - vs.d3 */
    		WLZ_VTX_3_SUB(vso.d3,vo.d3,vs.d3);
    
   	/* vs0.d3 dot n_i+1 */
     		dtemp1  =  WLZ_VTX_3_DOT(vso.d3, nIP1.d3);

   	/* get the s */
     		dtemp2 =  dtemp1  / dtemp;
     
      	/*  vtemp.d3   =   (  n_z * s  )  */
     		WLZ_VTX_3_SCALE(vtemp.d3, vtx1.d3,  dtemp2 );


    	/* now get the s_{i+1} */ 
     		WLZ_VTX_3_ADD(vs1.d3, vs.d3, vtemp.d3 );

   }

   /*---- get the s_{i+1}' = T_{i+1} s_{i+1}  ----*/
   errNum = Wlz3DSectionTransformVtxR(wlzViewStrip1, vs1.d3, &vtemp.d3 );
   sip1p->d3.vtX = vtemp.d3.vtX;
   sip1p->d3.vtY = vtemp.d3.vtY;
   sip1p->d3.vtZ = vtemp.d3.vtZ;



   printf("%lg  %lg  %lg\n", sip1p->d3.vtX, sip1p->d3.vtY, sip1p->d3.vtZ );


   return errNum;
}


/*!
* - Function:   outputTheWarpedWlz
* - Returns:    WlzVetex 
* - Purpose:    Track vertex through neighbour layer.	
* - Parameters:	
*
*     -#   vsp:              input   WlzVertex.
*     -#   vtp:              input   WlzVertex.
*     -#   wlzViewStri:      WlzThreeDViewStruct for layer i;
*     -#   wlzViewStrip1:    WlzThreeDViewStruct for the next layer i+-1;
*     -#   UpOrDown:         > 0 Up;  <= 0 down; 
* - Author:       J. Rao, R. Baldock
*/
static WlzErrorNum outputTheWarpedWlz( WlzVertex  vsp, 
                                       WlzVertex  vtp,
				       WlzVertex  vfp,
                                       WlzDVertex2 *vtxVec1,
                                       char *souce2DWlzStr,
				       char *out2DWlzStr
				     )
{
    WlzAffineTransform  *trans=NULL;

    WlzTransformType     tType = WLZ_TRANSFORM_2D_AFFINE;
    WlzDVertex2         *vtxVec0;
    WlzDVertex2          sv, tv;
    WlzObject           *sObj=NULL, *tObj=NULL;
    FILE                *inFile;
    WlzErrorNum	         errNum = WLZ_ERR_NONE;
    WlzInterpolationType interp = WLZ_INTERPOLATION_LINEAR;  /* Use the nearest neighbour */

   /*------- allocate memory --------*/
   if (
        ((vtxVec0 = (WlzDVertex2 *)AlcCalloc(sizeof(WlzDVertex2), 3)) == NULL) 
      )
   {
      	errNum = WLZ_ERR_MEM_ALLOC;
   }

   if( errNum == WLZ_ERR_NONE )
   {
    	(vtxVec0)->vtX     = vsp.d3.vtX;
    	(vtxVec0)->vtY     = vsp.d3.vtY;
    	(vtxVec0 + 1)->vtX = vtp.d3.vtX;
    	(vtxVec0 + 1)->vtY = vtp.d3.vtY;
    	(vtxVec0 + 2)->vtX = vfp.d3.vtX;
    	(vtxVec0 + 2)->vtY = vfp.d3.vtY;
    	trans = WlzAffineTransformLSq2D(3, vtxVec1, 3, vtxVec0, 0, NULL, tType, &errNum); 
        /* trans = WlzAffineTransformLSq2D(3, vtxVec0, 3, vtxVec1, 0, NULL, tType, &errNum); */
   }
   
   if( errNum == WLZ_ERR_NONE )
   {
      /* check whether it is right or not */
      sv.vtX = (vtxVec0)->vtX;
      sv.vtY = (vtxVec0)->vtY;
      printf("compare\n");
      
      printf("first one\n");
      printf("source: %lg   %lg\n", sv.vtX, sv.vtY );
      tv = WlzAffineTransformVertexD2(trans, sv, &errNum);
      printf("target: %lg   %lg\n", (vtxVec1)->vtX, (vtxVec1)->vtY );
      printf("source trans to target : %lg   %lg\n", tv.vtX, tv.vtY);
      
      sv.vtX = (vtxVec0 + 1)->vtX;
      sv.vtY = (vtxVec0 + 1)->vtY;
      printf("second one\n");
      printf("source: %lg   %lg\n", sv.vtX, sv.vtY );
      tv = WlzAffineTransformVertexD2(trans, sv, &errNum);
      printf("target: %lg   %lg\n", (vtxVec1 + 1)->vtX, (vtxVec1 + 1)->vtY );
      printf("source trans to target : %lg   %lg\n", tv.vtX, tv.vtY);

      sv.vtX = (vtxVec0 + 2)->vtX;
      sv.vtY = (vtxVec0 + 2)->vtY;
      printf("third one\n");
      printf("source: %lg   %lg\n", sv.vtX, sv.vtY );
      tv = WlzAffineTransformVertexD2(trans, sv, &errNum);
      printf("target: %lg   %lg\n", (vtxVec1 + 2)->vtX, (vtxVec1 + 2)->vtY );
      printf("source trans to target : %lg   %lg\n", tv.vtX, tv.vtY);
   
   }
	   
	   
   if( errNum == WLZ_ERR_NONE )
   {
    	/* read Woolz object */
        if((inFile = fopen(souce2DWlzStr, "r")) == NULL )
        {
        	printf("cannot open the input woolz file.\n");
        	exit(1);
      	}
      	if( !( sObj = WlzReadObj(inFile, &errNum) ) )
      	{
        	printf("input Woolz Object Error.\n");
        	fclose(inFile); 
        	exit(1);
      	}
      	fclose(inFile); 
      	inFile = NULL;
   }
    
   if( errNum == WLZ_ERR_NONE )
   {
   	tObj =  WlzAffineTransformObj(sObj, trans, interp, &errNum  );
   }

   if( errNum == WLZ_ERR_NONE )
   {
   	/* write Woolz object */
      	if((inFile = fopen(out2DWlzStr, "w")) == NULL )
      	{
        	printf("cannot open the output woolz file.\n");
        	exit(1);
      	}
      	errNum = WlzWriteObj(inFile, tObj); 

      	if(  errNum != WLZ_ERR_NONE )
      	{
        	printf("input Woolz Object Error.\n");
        	fclose(inFile); 
        	exit(1);
      	}
      	fclose(inFile); 
      	inFile = NULL;
	
   }
  
      WlzFreeObj(sObj);
      WlzFreeObj(tObj);
      AlcFree(vtxVec0 );
      AlcFree(trans);


   return errNum;
}

/*! input is the voxsel sip and the unit vector nStraghtline, indicating direction
    in OPT space.
    output is the *sip1p in view space!
*/
static WlzErrorNum WlzGetCorssPoint( WlzVertex sip,
                                WlzVertex nStraghtline,
                                WlzVertex *sip1p,
                                WlzThreeDViewStruct  *wlzViewStri  
				)
{
   WlzVertex    vtemp, vtemp0, vo, vso, vtx1;
   WlzVertex    vs1, nI;
   double       dtemp, dtemp1, dtemp2;
   WlzErrorNum	errNum = WLZ_ERR_NONE;
  

   /* output for test */
   /*
   printf("Yaw:      %f\n",  wlzViewStri->theta*180.0/3.14); 
   printf("Pitch:    %f\n",  wlzViewStri->phi*180.0/3.14 );
   printf("Roll:     %f\n",  wlzViewStri->zeta*180.0/3.14 );
   printf("distance: %f\n",  wlzViewStri->dist );
   printf("scaling:  %f\n",  wlzViewStri->scale );
   printf("%f  %f  %f\n",  wlzViewStri->fixed.vtX, wlzViewStri->fixed.vtY ,wlzViewStri->fixed.vtZ  );
   */

   vtx1.d3.vtX = 0.0;
   vtx1.d3.vtY = 0.0;
   vtx1.d3.vtZ = 1.0;
   /* get the normal n_i = nI.d3   of this plan i -plane do not need !!! */
   errNum = Wlz3DSectionTransformInvVtxR(wlzViewStri, vtx1.d3, &nI.d3 );

   /* in OPT space */

   /* check  the n_{i} *  nStraghtline != 0  other wise throw an erro  */
   dtemp =  WLZ_VTX_3_DOT(nStraghtline.d3, nI.d3);
   /* dtemp = abs(dtemp); */
   if( ( dtemp <  0.000001 ) && ( dtemp >  -0.000001 ) )
   {
      printf("%lg\n",dtemp);
      printf(" the  n_{i} is verticl to n_straght line in OPT space,\n");
      printf(" there is no crossed point. Please chose different straight\n");
      printf("line and  try again\n ");
      printf("%lg  %lg  %lg\n", nI.d3.vtX, nI.d3.vtY, nI.d3.vtZ);
      printf("%lg  %lg  %lg\n", nStraghtline.d3.vtX, nStraghtline.d3.vtY, nStraghtline.d3.vtZ);
      exit( 1 );
   }
   else
   {
     /* not parallel to z-direction  */
     /*  s_{i} = x_{i} + n_straight_line * s
                 = x_{i} + n_straight_line * 
		     <  [ O_{i} - x_{i}  ] dot n_{i}  /   ( n_straightline dot n_{i} ) >

     */
    
     	/*  dtemp =  ( n_straightlien dot n_{i} )  aleady got */

     	/* get the  o_i+1 
   		vo.d3.vtX   = wlzViewStr1->fixed.vtX + wlzViewStr->dist * nI.d3.vtX;
   		vo.d3.vtY   = wlzViewStr1->fixed.vtY + wlzViewStr->dist * nI.d3.vtY;
   		vo.d3.vtY   = wlzViewStr1->fixed.vtZ + wlzViewStr->dist * nI.d3.vtZ;
    	*/
    		WLZ_VTX_3_SCALE(vtemp0.d3, nI.d3, wlzViewStri->dist);
    		WLZ_VTX_3_ADD(vo.d3, wlzViewStri->fixed, vtemp0.d3);

    		/* vso.d3 = vo.d3 - vx.d3 */
    		WLZ_VTX_3_SUB(vso.d3,vo.d3,sip.d3);
    
   	/* vs0.d3 dot n_i */
     		dtemp1  =  WLZ_VTX_3_DOT(vso.d3, nI.d3);

   	/* get the s */
     		dtemp2 =  dtemp1  / dtemp;
     
      	/*  vtemp.d3   =   (  n_straightlien * s  )  */
     		WLZ_VTX_3_SCALE(vtemp.d3, nStraghtline.d3,  dtemp2 );


    	/* now get the p_{i} */ 
     		WLZ_VTX_3_ADD(vs1.d3, sip.d3, vtemp.d3 );

   }

   /*---- get the p_{i}' = T_{i} p_{i}  ----*/
   errNum = Wlz3DSectionTransformVtxR(wlzViewStri, vs1.d3, &vtemp.d3 );
   sip1p->d3.vtX = vtemp.d3.vtX;
   sip1p->d3.vtY = vtemp.d3.vtY;
   sip1p->d3.vtZ = vtemp.d3.vtZ;

   printf("%lg  %lg  %lg\n", sip1p->d3.vtX, sip1p->d3.vtY, sip1p->d3.vtZ );


   return errNum;



}

/*!
   This funciton can be used to get tartget points in Stack up 
   space.
   Where nStraight line is a unit vect indicating the directions
   in stack up space.
   dz is the length in stackup space relative to the main plane.

   sip1p will store the output results.


*/
static WlzErrorNum WlzGetDxDy(  
                                WlzVertex nStraghtline,
				double    dz, 
                                WlzVertex *sip1p
				)
{
   double       dtemp;
   WlzErrorNum	errNum = WLZ_ERR_NONE;
  
   /* in stack up space */

   /* check  the n_{z} *  nStraghtline != 0  other wise throw an erro  */
   dtemp =  nStraghtline.d3.vtZ;
   /* dtemp = abs(dtemp); */
   if( ( dtemp <  0.000001 ) && ( dtemp >  -0.000001 ) )
   {
      printf("%lg\n",dtemp);
      printf(" the  z-compenent of n_straght line in OPT space,\n");
      printf(" is zero. Please chose different straight\n");
      printf("line and  try again\n ");
      printf("%lg  %lg  %lg\n", nStraghtline.d3.vtX, nStraghtline.d3.vtY, nStraghtline.d3.vtZ);
      exit( 1 );
   }
   else
   {
     /* not crossed to z-direction  */
     /*  dx = dz (n_d dot n_x )/ (n_d dot n_z )
	 dy = dz (n_d dot n_y )/ (n_d dot n_z )
     */

       /*---- get the p_{i}' = T_{i} p_{i}  ----*/
   	sip1p->d3.vtX = dz * nStraghtline.d3.vtX/dtemp;
   	sip1p->d3.vtY = dz * nStraghtline.d3.vtY/dtemp;
   	sip1p->d3.vtZ = dz;

   }

   printf("%lg  %lg  %lg\n", sip1p->d3.vtX, sip1p->d3.vtY, sip1p->d3.vtZ );


   return errNum;

}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
