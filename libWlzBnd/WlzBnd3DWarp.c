    #include <WlzBnd.h>


    WlzObject  *WlzBnd3DWarpFile(char *inFileStr, char *TiePointsFileStr, char *outFileStr, WlzErrorNum *dstErr)
    {
    FILE *inFile = NULL;   /* used to read Woolz object */
    int option,i;
    int nbelow = 2, nup = 0;
    int ReadDisplacement = 0;
    int outputAutoMeshVTK = 0,
    outputAutoMeshWLZ = 0,
    outputTransformedMeshVTK = 0,
    outputTransformedMeshWLZ = 0,
    outputCutPlaneAndCorrepSurfaceVTK = 0,
    ReadMeshTransformFromFile= 0,
    WARP = 1;
    int numOfElemAlonX = 6,
    numOfElemAlonY = 6,
    numOfElemAlonZ = 6;
    double xmin = -10.0, xmax = 10.0,
    ymin = -10.0, ymax = 10.0,
    zmin = -10.0, zmax = 10.0,
    zConst = 0.;
    /*char                *inFileStr, *outFileStr, *TiePointsFileStr;*/

    const char *errMsg;
    WlzErrorNum errNum = WLZ_ERR_NONE;
    WlzMeshTransform3D *wmt3D = NULL;
    WlzBasisFnTransform *basisTr = NULL;   /* the transformation functions  */
    WlzFnType basisFnType = WLZ_FN_BASIS_3DMQ;
    WlzInterpolationType interp = WLZ_INTERPOLATION_NEAREST;  /* Use the nearest neighbour */
    int basisFnPolyOrder = 3;
    WlzDVertex3 vx4,vx5;
    WlzIBox3 bBox0 = {0,0,0,-1,-1,-1}; /* store the output Bounding box */
    WlzDBox3 bBoxS = {0.,0.,0.,0.,0.,0.};
    WlzObject *wObjS = NULL; /* source Woolz object */
    WlzObject *wObjW = NULL; /* Warped Woolz object */
    WlzDVertex3 *vxVec0;   /* source points */
    WlzDVertex3 *vxVec1;   /* correponding displacements or target points */
    int nTiePP;   /* number of tie points */
    char *cstr;
    WlzMeshTransform2D5 *wmtCut2D5 = NULL;    /* for plane 2D triangle mesh and 3D displacement */

    /* read Woolz object */
    if((inFile = fopen(inFileStr, "rb")) == NULL )
    {
    printf("cannot open the input woolz file.\n");
    errNum = WLZ_ERR_FILE_OPEN;
    }

    if( !(wObjS = WlzReadObj(inFile, &errNum) ) )
    {
    printf("input Woolz Object Error.\n");
    fclose(inFile);
    errNum = WLZ_ERR_READ_EOF;
    }
    fclose(inFile);
    inFile = NULL;

    if(errNum == WLZ_ERR_NONE) {

    if(bBoxS.xMin == bBoxS.xMax)
    {
    bBoxS = WlzBoundingBox3D(wObjS, &errNum);
    /*printf("%6.2f  %6.2f %6.2f %6.2f  %6.2f %6.2f\n",
    bBoxS.xMin, bBoxS.xMax,
    bBoxS.yMin, bBoxS.yMax,
    bBoxS.zMin, bBoxS.zMax
    );*/
    }
    }

    if(errNum == WLZ_ERR_NONE) {
    /*-- get the mesh ---*/
    if(ReadMeshTransformFromFile)
    {
    if((inFile = fopen("MeshTransform3D.wlz", "r")) == NULL )
    {
        printf("cannot open the input file of MeshTransform3D.wlz.\n");
        errNum = WLZ_ERR_FILE_OPEN;
    }
    wmt3D = (WlzMeshTransform3D *) WlzReadMeshTransform3D(inFile, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
        printf("input MeshTransform3D in woolz failed");
        errNum = WLZ_ERR_READ_INCOMPLETE;
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
        errNum = WLZ_ERR_WRITE_INCOMPLETE;
        }
        WlzEffWriteMeshTransform3DWithoutDisplacementVTK(inFile, wmt3D);
        fclose(inFile);
        inFile = NULL;
        }

        /* read the tie-points */
        if((inFile = fopen(TiePointsFileStr, "r")) == NULL )
        {
        printf("cannot open the input tie point data file.\n");
        errNum = WLZ_ERR_FILE_OPEN;
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
        errNum = WLZ_ERR_TRANSFORM_DATA;
        }
        }

        /* output transformed mesh in VTK format  */
        if(outputTransformedMeshVTK)
        {
        cstr = "transformedMesh.vtk";
        if((inFile = fopen(cstr, "w")) == NULL )
        {
        printf("cannot open the output transformedMesh.vtk file.\n");
        errNum = WLZ_ERR_WRITE_EOF;
        }
        WlzEffWriteMeshTransform3DWithDisplacementVTK(inFile, wmt3D);
        fclose(inFile);
        inFile = NULL;
        }

        /* output transformed mesh in Woolz format */
        if(outputTransformedMeshWLZ)
        {
        if((inFile = fopen("MeshTransform3D.wlz", "w")) == NULL )
        {
        printf("cannot open the output file of MeshTransform3D.wlz.\n");
        errNum = WLZ_ERR_WRITE_EOF;
        }
        if(WlzWriteMeshTransform3D(inFile, wmt3D) != WLZ_ERR_NONE)
        {
        printf("output MeshTransform3D in woolz failed");
        errNum = WLZ_ERR_WRITE_EOF;
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
    errNum = WLZ_ERR_WRITE_EOF;
    }
    WlzEffWriteOriginalPlaneVTKByPos(inFile, wmtCut2D5);
    fclose(inFile);
    inFile = NULL;

    /* output the corresponding original surface */
    if((inFile = fopen("CorrespOriginalSurface.vtk", "w")) == NULL )
    {
    printf("cannot open the output file for CorrespOriginalSurface.vtk \n");
    errNum = WLZ_ERR_WRITE_EOF;
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
    outFileStr = "outWarpWlz.wlz";
    if((inFile = fopen(outFileStr, "w")) == NULL )
    {
    printf("cannot open output file.\n");
    errNum = WLZ_ERR_WRITE_EOF;
    }

    if( (errNum = WlzWriteObj(inFile, wObjW) ) !=WLZ_ERR_NONE )
    {
    printf("can not out put the transformed Woolz Object Error.\n");
    fclose(inFile);
    errNum = WLZ_ERR_WRITE_EOF;
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

    if( dstErr )
    *dstErr = errNum;

    return wObjW;
    }


WlzObject  *WlzBnd3DWarpObj(WlzObject *ObjS, int arraySizeVec0, WlzDVertex3 *arrayVec0, int arraySizeVec1, WlzDVertex3 *arrayVec1, WlzErrorNum *dstErr)
{
  int numOfElemAlonX = 6, numOfElemAlonY = 6, numOfElemAlonZ = 6;
  WlzErrorNum errNum = WLZ_ERR_NONE;
  WlzMeshTransform3D *wmt3D = NULL;
  WlzBasisFnTransform *basisTr = NULL;
  WlzFnType basisFnType = WLZ_FN_BASIS_3DMQ;
  WlzInterpolationType interp = WLZ_INTERPOLATION_NEAREST;
  int basisFnPolyOrder = 3;
  WlzDBox3 bBoxS = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  WlzObject *ObjT = NULL;

  if (NULL == ObjS)
    errNum = WLZ_ERR_OBJECT_NULL;

  if(errNum == WLZ_ERR_NONE)
    bBoxS = WlzBoundingBox3D(ObjS, &errNum);

  if(errNum == WLZ_ERR_NONE) {
    if( (wmt3D = WlzTetrahedronMeshFromObj(ObjS, bBoxS,numOfElemAlonX,
      numOfElemAlonY, numOfElemAlonZ, &errNum) ) != NULL)
    {
      basisTr = WlzBasisFnTrFromCPts3( basisFnType, basisFnPolyOrder,
        arraySizeVec0, arrayVec0, arraySizeVec1, arrayVec1, &errNum );

    if( WlzGetTransformedMesh(wmt3D, basisTr) != WLZ_ERR_NONE )
      errNum = WLZ_ERR_TRANSFORM_DATA;
    }
  }

  ObjT = WlzMeshTransformObj_3D(ObjS,  wmt3D, interp, &errNum);

  AlcFree(wmt3D->elements);
  AlcFree(wmt3D->nodes);
  AlcFree(wmt3D);

  if( dstErr )
    *dstErr = errNum;

  return ObjT;
}


