#ifndef DOXYGEN_SHOULD_SKIP_THIS
/***********************************************************************
* \Project:      Woolz
* \Title:        Wlz3DWarpMQ.c
* \Date:         18 January 2002 
* \Author:       J. Rao
* \Copyright:	 2002 Medical Research Council, UK.
*		 All rights reserved.
* \Address: 	 MRC Human Genetics Unit,
*		 Western General Hospital,
*		 Edinburgh, EH4 2XU, UK.
* \Purpose:      Woolz 3D Warping using MQ 
*		
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <Wlz.h>

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

 
int main(int	argc,
	 char	**argv)
{

  FILE	       *inFile = NULL;   /* used to read Woolz object */ 
  int		option,i;
  int           nbelow                   = 2, 
                nup                      = 0;
  int           ReadDisplacement         = 0;
  int           outputAutoMeshVTK        = 0,
                outputAutoMeshWLZ        = 0,
                outputTransformedMeshVTK = 0,
                outputTransformedMeshWLZ = 0,
		outputCutPlaneAndCorrepSurfaceVTK        = 0,
		ReadMeshTransformFromFile= 0,
		WARP                     = 1;
  int           numOfElemAlonX = 6,
                numOfElemAlonY = 6,
                numOfElemAlonZ = 6;
  double        xmin = -10.0, xmax = 10.0,
                ymin = -10.0, ymax = 10.0,
	        zmin = -10.0, zmax = 10.0,
		       zConst = 0.;
  char                *inFileStr, *outFileStr, *TiePointsFileStr;
  char                *outFileMeshTrWlzStr, *inFileMeshTrWlzStr;
  const char	      *errMsg;
  WlzErrorNum	       errNum = WLZ_ERR_NONE;
  WlzMeshTransform3D  *wmt3D = NULL;
  WlzBasisFnTransform *basisTr       = NULL;   /* the transformation functions  */
  WlzFnType            basisFnType    = WLZ_FN_BASIS_3DMQ;
  WlzInterpolationType interp = WLZ_INTERPOLATION_NEAREST;  /* Use the nearest neighbour */
  int                  basisFnPolyOrder = 3;
  WlzDVertex3 vx4,vx5;
  WlzIBox3             bBox0 = {0,0,0,-1,-1,-1}; /* store the output Bounding box */
  WlzDBox3             bBoxS = {0.,0.,0.,0.,0.,0.};
  WlzObject           *wObjS          = NULL; /* source Woolz object */
  WlzObject           *wObjW          = NULL; /* Warped Woolz object */
  WlzDVertex3         *vxVec0;   /* source points */
  WlzDVertex3         *vxVec1;   /* correponding displacements or target points */
  int                  nTiePP;   /* number of tie points */
  char                *cstr;
  WlzMeshTransform2D5 *wmtCut2D5        = NULL;    /* for plane 2D triangle mesh and 3D displacement */

  /* read the argument list and check for an input file */

  static char	optList[] = "l:m:n:x:y:z:X:Y:Z:i:j:k:O:o:p:t:M:c:C:b:u:r:R:e:E:f:F:g:G:a:A:Q:q:W:B:hGL",

  opterr = 0;
  /*
  bBox0.xMin = 0;
  bBox0.xMax = -1;
  bBox0.yMin = 0;
  bBox0.yMax = -1;
  bBox0.zMin = 0;
  bBox0.zMax = -1;
  */
 
   while( (option = getopt(argc, argv, optList)) != EOF )
   {
      switch( option )
      {
        case 'h':
             usage(argv[0]);
	     return(0);
	case 'x':
	    if(sscanf(optarg, "%lg", &bBoxS.xMin) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;

        case 'y':
	    if(sscanf(optarg, "%lg", &bBoxS.yMin) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;

        case 'z':
	    if(sscanf(optarg, "%lg", &bBoxS.zMin) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;

	case 'X':
	    if(sscanf(optarg, "%lg", &bBoxS.xMax) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;

        case 'Y':
	    if(sscanf(optarg, "%lg", &bBoxS.yMax) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;

        case 'Z':
	    if(sscanf(optarg, "%lg", &bBoxS.zMax) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;

	case 'l':
	    if(sscanf(optarg, "%d", &numOfElemAlonX) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;

        case 'm':
	    if(sscanf(optarg, "%d", &numOfElemAlonY) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;

        case 'n':
	    if(sscanf(optarg, "%d", &numOfElemAlonZ) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;

        case 'i':
	    inFileStr = optarg;
	    break;
        case 'o':
	    outFileStr = optarg;
	    break;
        case 't':
	    TiePointsFileStr = optarg;
	    break;
         case 'q':
	    outFileMeshTrWlzStr = optarg;
	    break;
        case 'B':
	    inFileMeshTrWlzStr = optarg;
	    break;
        case 'c':
	    if(sscanf(optarg, "%lg", &zConst) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;
        case 'C':
	    if(sscanf(optarg, "%d", &outputCutPlaneAndCorrepSurfaceVTK) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;
        case 'a':
	    if(sscanf(optarg, "%d", &outputAutoMeshVTK) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;
        case 'A':
	    if(sscanf(optarg, "%d", &outputAutoMeshWLZ) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;
         case 'M':
	    if(sscanf(optarg, "%d", &outputTransformedMeshVTK) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;
         case 'Q':
	    if(sscanf(optarg, "%d", &outputTransformedMeshWLZ) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;
         case 'r':
	    if(sscanf(optarg, "%d", &ReadMeshTransformFromFile) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;
         case 'b':
	    if(sscanf(optarg, "%d", &nbelow) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;
        case 'u':
	    if(sscanf(optarg, "%d", &nup) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;
        case 'R':
	    if(sscanf(optarg, "%d", &ReadDisplacement) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;
        case 'e':
	    if(sscanf(optarg, "%d", &bBox0.xMin) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;
        case 'E':
	    if(sscanf(optarg, "%d", &bBox0.xMax) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;
        case 'f':
	    if(sscanf(optarg, "%d", &bBox0.yMin) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;
        case 'F':
	    if(sscanf(optarg, "%d", &bBox0.yMax) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;
        case 'g':
	    if(sscanf(optarg, "%d", &bBox0.zMin) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;
        case 'G':
	    if(sscanf(optarg, "%d", &bBox0.zMax) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;
        case 'W':
	    if(sscanf(optarg, "%d", &WARP) != 1)
	    {
	      printf("read error");
	      exit(1);
	    }
	    break;
	case 'p':
	case 'L':
        default:
              return(0);
      }
   }
  /*

  if(errNum != WLZ_ERR_NONE) {
    (void )WlzStringFromErrorNum(errNum, &errMsg);
    (void )fprintf(stderr,
    		   "%s: failed to make boundary list (%s)\n",
		   argv[0], errMsg);
  }   */

  /* read Woolz object */
  if((inFile = fopen(inFileStr, "r")) == NULL )
  {
    printf("cannot open the input woolz file.\n");
    exit(1);
  }

  if( !(wObjS = WlzReadObj(inFile, &errNum) ) )
  {
    printf("input Woolz Object Error.\n");
    fclose(inFile); 
    exit(1);
  }
  fclose(inFile); 
  inFile = NULL;

  if(errNum == WLZ_ERR_NONE) {

    if(bBoxS.xMin == bBoxS.xMax)
    {
      bBoxS = WlzBoundingBox3D(wObjS, &errNum);
      /*
      printf("%6.2f  %6.2f %6.2f %6.2f  %6.2f %6.2f\n",
            bBoxS.xMin, bBoxS.xMax,
            bBoxS.yMin, bBoxS.yMax,
            bBoxS.zMin, bBoxS.zMax
	    );
      */	    
    }

  }

   if(errNum == WLZ_ERR_NONE) {
        /*-- get the mesh ---*/

        if(ReadMeshTransformFromFile)
        {
            if((inFile = fopen(inFileMeshTrWlzStr, "r")) == NULL )
            {
               printf("cannot open the input file of MeshTransform3D.wlz.\n");
               exit(1);
            }
            wmt3D = (WlzMeshTransform3D *) WlzReadMeshTransform3D(inFile, &errNum);
            if(errNum != WLZ_ERR_NONE)
            {
               printf("input MeshTransform3D in woolz failed");
               exit(1);
            }
            fclose(inFile);
            inFile = NULL;
   
        }
        else
        {

           if( (wmt3D = WlzTetrahedronMeshFromObj(wObjS, bBoxS,numOfElemAlonX, 
                                           numOfElemAlonY, numOfElemAlonZ, &errNum) ) != NULL);
           {
              /* output orginal mesh in VTK format */
              if(outputAutoMeshVTK)
              {
                 cstr = "orginalMeshAutoProduced.vtk";
                 if((inFile = fopen(cstr, "w")) == NULL )
                 {
                    printf("cannot open the output  data file.\n");
                    exit(1);
                 }
                 WlzEffWriteMeshTransform3DWithoutDisplacementVTK(inFile, wmt3D);
                 fclose(inFile);
                 inFile = NULL;
              }

              /* read the tie-points */
              if((inFile = fopen(TiePointsFileStr, "r")) == NULL )
              {
                     printf("cannot open the input tie point data file.\n");
                     exit(1);
              }
              errNum = read_WlzTiePoints(inFile, &nTiePP, &vxVec0, &vxVec1, ReadDisplacement);
              if(inFile)
              {
                 fclose(inFile);
                 inFile = NULL;
              }
	      
              /* get the basis transform */
              basisTr = WlzBasisFnTrFromCPts3( basisFnType, basisFnPolyOrder,
    			                nTiePP, vxVec0,
				        nTiePP, vxVec1, &errNum );
              /* get the mesh transformation  */
              { /* need consistency with 2D format */
                if( WlzGetTransformedMesh(wmt3D, basisTr) != WLZ_ERR_NONE )
                {
                     printf("Can't get the mesh transform. ");
                     exit(1);
                }
              }

              /* output transformed mesh in VTK format  */
              if(outputTransformedMeshVTK)
              {
                  cstr = "transformedMesh.vtk";
                  if((inFile = fopen(cstr, "w")) == NULL )
                  {
                    printf("cannot open the output transformedMesh.vtk file.\n");
                    exit(1);
                  }
                  WlzEffWriteMeshTransform3DWithDisplacementVTK(inFile, wmt3D);
                  fclose(inFile);
                  inFile = NULL;
               }

               /* output transformed mesh in Woolz format */
               if(outputTransformedMeshWLZ)
               {
                  if((inFile = fopen(outFileMeshTrWlzStr, "w")) == NULL )
                   {
                      printf("cannot open the output file of MeshTransform3D.wlz.\n");
                      exit(1);
                   }
                   if(WlzWriteMeshTransform3D(inFile, wmt3D) != WLZ_ERR_NONE)
                   {
                      printf("output MeshTransform3D in woolz failed");
                      exit(1);
                   }
                   fclose(inFile);
                   inFile = NULL;
	       }
	   }

	}
   }


  /*--- want a specific cut section and its origina surface ? ----*/
  if(outputCutPlaneAndCorrepSurfaceVTK)
  {  
     /* get the 2D5 meshtransform first */
     wmtCut2D5  = Wlz2D5TransformFromCut3Dmesh(zConst, wmt3D, &errNum);

     /* output cuting plane */
     if((inFile = fopen("CutPlane.vtk", "w")) == NULL )
     {
         printf("cannot open the output file for originalPlan.vtk \n");
         exit(1);
     }
     WlzEffWriteOriginalPlaneVTKByPos(inFile, wmtCut2D5);
     fclose(inFile);
     inFile = NULL;

     /* output the corresponding original surface */
     if((inFile = fopen("CorrespOriginalSurface.vtk", "w")) == NULL )
     {
         printf("cannot open the output file for CorrespOriginalSurface.vtk \n");
         exit(1);
     }
     WlzEffWriteOriginalPlaneVTKByDis(inFile, wmtCut2D5);
     fclose(inFile);
     inFile = NULL;
     AlcFree(wmtCut2D5->elements);
     AlcFree(wmtCut2D5->nodes);
     AlcFree(wmtCut2D5);
  }

  /*--- warping ----*/
  if(WARP)
  {
    wObjW = WlzMeshTransformObj_3D(wObjS,  wmt3D, interp, &errNum);
    /* output the woolz object */
    if((inFile = fopen(outFileStr, "w")) == NULL )
    {
       printf("cannot open output file.\n");
       exit(1);
    } 
  
    if( (errNum =	WlzWriteObj(inFile, wObjW) ) !=WLZ_ERR_NONE )
    {
       printf("can not out put the transformed Woolz Object Error.\n");
       fclose(inFile);
       exit(1);
    }
    fclose(inFile);
    inFile = NULL;
    AlcFree(wObjW);
  }  
  
  AlcFree(wmt3D->elements);
  AlcFree(wmt3D->nodes);
  AlcFree(wmt3D);
  AlcFree(vxVec0);
  AlcFree(vxVec1);
  AlcFree(wObjS);
  return ( 0 );
}


static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s [-h] [<input file>]\n"
	  "\tAutomatically produce a tetrahedron mesh from a woolz object. But\n"
	  "\tthe user should input some parameters to give a large cuboid which\n"
	  "\tcover the Woolz object. Sure, you can also just give a small cuboid\n"
	  "\tby transfer only part of your Woolz object. \n"
	  "\n"
	  "\tPlease report bugs to RAO.JIANGUO@hgu.mrc.ac.uk\n"
	  "\n"
	  "\tOptions are:\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -l        input the number of cuboid elements along x-direction\n"
	  "\t                                      default is  6\n"
	  "\t  -m        input the number of cuboid elements along y-direction\n"
	  "\t                                      default is  6\n"
	  "\t  -n        input the number of cuboid elements along z-direction\n"
	  "\t                                      default is  6\n"
	  "\t  -x        input the minimum x-coordinate\n"
	  "\t  -y        input the minimum y-coordinate\n"
	  "\t  -z        input the minimum z-coordinate\n"
	  "\t  -X        input the maximum X-coordinate\n"
	  "\t  -Y        input the maximum Y-coordinate\n"
	  "\t  -Z        input the maximum Z-coordinate\n"
	  "\t            default will use the input Woolz Obj boxding box\n"
	  "\t  -C        output cuting plane and corresponding\n"
	  "\t                 surface in VTK form ?\n"
	  "\t                                 1 for yes\n"
	  "\t                                 0 for  no\n"
	  "\t                            default is  no\n"
	  "\t                                          \n"
	  "\t  -c        input the constant Z-plane used to cut the transformed mesh\n"
	  "\t                                      default is    0\n"
          "\t  -i        input file name for the Woolz Object\n"
          "\t  -o        output file name for the warped Woolz Object\n"
          "\t  -t        input file name for the tiepoints\n"
          "\t  -q        input file name for output Mesh Transform in Wlz format\n"
          "\t  -B        input file name for input  Mesh Transform in Wlz Format\n"
	  "\t  -g        input the minimum z-coordinate target plane\n"
	  "\t                                      default is  0\n"
	  "\t  -G        input the maximum Z-coordinate target plane\n"
	  "\t                                      default is   -1\n"
	  "\t  -R        The tiepoint contains displacement \n"
	  "\t                                                       1 for yes\n"
	  "\t                                                       0 for  no\n"
	  "\t                                                  default is  no\n"
	  "\t                                                                \n"
	  "\t  -r        Read Mesh transform from a file\n"
	  "\t                                                       1 for yes\n"
	  "\t                                                       0 for  no\n"
	  "\t                                                  default is  no\n"
	  "\t                                                                \n"
	  "\t  -Q        output Mesh transform in Wlz format to a file\n"
	  "\t                                                       1 for yes\n"
	  "\t                                                       0 for  no\n"
	  "\t                                                  default is  no\n"
	  "\t                                                                \n"
	  "\t  -a        output the automatically produced 3D mesh VTK form ?\n"
	  "\t                                                       1 for yes\n"
	  "\t                                                       0 for  no\n"
	  "\t                                                  default is  no\n"
	  "\t                                                    \n"
	  "\t  -M        output the transformed 3D mesh VTK form?\n"
	  "\t                                           1 for yes\n"
	  "\t                                           0 for  no\n"
	  "\t                                      default is  no\n"
	  "\t                                                    \n"
	  "\t  -W        Warp the WlzObj ?\n"
	  "\t                                             1 for yes\n"
	  "\t                                             0 for  no\n"
	  "\t                                        default is  YES\n"
	  "\t                                                      \n"
					  "",
	  proc_str);
  return;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
