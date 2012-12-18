#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _Wlz3DWarpMQ_S_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/Wlz3DWarpMQ_S.c
* \author       Jianguo Rao
* \date         September 2003
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
* \brief	Generates a regular tetrahedral mesh for a 3D domain.
* \ingroup	WlzTransform
*/

/* This code started at 13/09/2001 by J. Rao 13/09/2001 
         Richard and Bill have also made contributions
         it experienced a great deal changes 26/10/2001 
         Here is the first version and it can be used  
	 for 3D warping.  /18/01/2002/           */

#include <stdlib.h>
#include <ctype.h> /* this is just for isspace() function */
#include <float.h>
#include <Wlz.h>

#define	IN_RECORD_MAX   (1024)


#undef DEBUGME

typedef struct _WlzMeshScanDElm
{
  int		valid;				        /* Non-zero if valid */
  double	xTr[3];         /* Affine transform coefficients for columns */
  double	yTr[3];           /* Affine transform coefficients for lines */
} WlzMeshScanDElm;

typedef struct _WlzMeshScanItv
{
  int		elmIdx;
  int		line;
  int		lftI;
  int		rgtI;
} WlzMeshScanItv;

typedef struct _WlzMeshScanWSp2D5
{
  WlzMeshTransform2D5 *mesh;	 /* the mesh transform */
  int		nItvs; 	         /* Number of element intervals */
  WlzMeshScanItv *itvs;          /* Element intervals sorted by line then left column */
  WlzMeshScanDElm *dElm;         /* Destination mesh element data */
} WlzMeshScanWSp2D5;

static void	WlzMeshScanWSpFree(WlzMeshScanWSp2D5 *);

static int	WlzMeshItvCmp(const void *, const void *);

static int      WlzMeshScanTriElm(WlzMeshScanWSp2D5 *, int, int);

static WlzMeshScanWSp2D5 *WlzMeshScanWSpInit2D5(WlzMeshTransform2D5 *mesh,
				    	  WlzErrorNum *dstErr);

double WlzGeomTriangle2SnArea3(WlzDVertex3 *a,  int i, int j,int k);

void WlzEffWriteMeshTransform3DWithoutDisplacementVTK(FILE *fp, 
                                              WlzMeshTransform3D *wmt3D);

void WlzEffWriteMeshTransform3DWithDisplacementVTK(FILE *fp, 
                                              WlzMeshTransform3D *wmt3D);
#ifdef DEBUGME
void static WlzEffWriteCutplanTetrahedronVTK(FILE *fp, 
                                const WlzMeshTransform3D *wmt3D,  
				int *intersectIndex, 
				int nIntersect );
#endif
#ifdef DEBUGME
void static WlzEffWriteOriginalPlaneVTK(FILE *fp, 
			   const WlzDVertex3 *planepointsO, 
                           int **linkList, 
			   const int nIntersect,  
			   const int noRedundancyNumP);
#endif

void WlzEffWriteOriginalPlaneVTKByDis(FILE *fp, WlzMeshTransform2D5 *wmt2D5);

void WlzEffWriteOriginalPlaneVTKByPos(FILE *fp, WlzMeshTransform2D5 *wmt2D5);


void static Make2D5TriangleMeshFromCuttedPlane(WlzMeshTransform2D5 *wmt2D5, 
                                       WlzDVertex3 *planepoints, 
                                       int **linkList, 
				       int nIntersect,  
				       int noRedundancyNumP);
				       
void static GetTheSourceSurfaceOfTheCutPlanePoints(int nIntersect, 
                                             int *intersectIndex,
                                             const WlzMeshTransform3D  *wmt3D, 
					     WlzDVertex3 *planepoints, 
					     WlzDVertex3 *planepointsO, 
					     int **linkList
					     );

int static storedPointIndex(WlzDVertex3 tempw, WlzDVertex3 *planepoints, int icount);

WlzErrorNum  write_Wlz2D5Mesh(FILE *fp,  char *cstr, WlzMeshTransform2D5 *wmt2D5);

WlzErrorNum  Write_WlzCutScanLines(FILE *fp,  char *cstr, WlzMeshScanWSp2D5 *mSnWSp);

WlzErrorNum  WlzGetTransformedMesh(WlzMeshTransform3D *wmt3D, WlzBasisFnTransform* basisTr);

void static Make2D5MeshDisplacement( WlzMeshTransform2D5 *wmt2D5,
				     WlzDVertex3 *planepointsO );
#ifdef DEBUGME
int  static GetIntersectTetrahedronNumber(double zConst, WlzMeshTransform3D *wmt3D);
#endif


WlzObject       *WlzMeshTransformObj_3D( WlzObject            *srcObj,
				         WlzMeshTransform3D   *wmt3D,
				         WlzInterpolationType  interp,
 				         WlzErrorNum          *dstErr
                                        );



#ifdef DEBUGME
static WlzObject	*WlzMeshTransformObj3D(WlzObject            *srcObj,
				       WlzMeshTransform3D   *wmt3D,
				       WlzMeshTransform2D5  *wmt2D5,
				       WlzInterpolationType  interp,
				       WlzErrorNum          *dstErr,
                                 int         *intersectIndex,
				 WlzDVertex3 *planepoints,
				 WlzDVertex3 *planepointsO,
				 int        **linkList,
			         int          nInterSectEsimate,
                                 int          numEstTrangularsElem,
                                 int          numEstCutingPoints,
				 WlzMeshScanWSp2D5 *mSnWSp,
				 WlzIBox3        bBox0
				       );
#endif
#ifdef DEBUGME
static WlzObject *WlzMeshTransformObjPrv3D(WlzObject *srcObj,
				     WlzMeshTransform3D *mesh,
				     WlzMeshTransform2D5 *wmt2D5,
				     WlzInterpolationType interp,
				     WlzErrorNum *dstErr,
                                 int         *intersectIndex,
				 WlzDVertex3 *planepoints,
				 WlzDVertex3 *planepointsO,
				 int        **linkList,
				 int          nInterSectEsimate,
                                 int          numEstTrangularsElem,
                                 int          numEstCutingPoints,
                                 WlzMeshScanWSp2D5 *mSnWSp,
				 WlzIBox3        bBox0
				     );
#endif

WlzIBox3  WlzMQ3DTransformBBoxI3( WlzMeshTransform3D *wmt3D,  
                                  WlzErrorNum        *errNum );

WlzErrorNum static WlzMeshTransformValues3D( WlzObject           *dstObj,
				      WlzObject           *srcObj,
			 	      WlzMeshTransform3D  *mesh,
				      WlzMeshTransform2D5 *wmt2D5,
				      int                  cutPosition,
				      WlzInterpolationType interp,
                                 int         *intersectIndex,
				 WlzDVertex3 *planepoints,
				 WlzDVertex3 *planepointsO,
				 int        **linkList,
				 int          nInterSectEsimate,
                                 int          numEstTrangularsElem,
                                 int          numEstCutingPoints,
				 WlzMeshScanWSp2D5 *mSnWSp
				      );


void  static WlzGet2D5TrangularMeshFrom3DMesh(  WlzMeshTransform2D5 *wmt2D5,
                                      WlzMeshTransform3D  *mesh, 
				      int                  cutPositon, 
				      WlzErrorNum         *errNum,
				 int         *intersectIndex,
				 WlzDVertex3 *planepoints,
				 WlzDVertex3 *planepointsO,
				 int        **linkList,
				 int          nInterSectEsimate,
				 int          numEstTrangularsElem,
                                 int          numEstCutingPoints
				      );

int     WlzIsSinglePointCutMesh(const WlzMeshTransform2D5 *wmt2D5);

int     WlzNoAreaGreaterThanOne(const WlzMeshTransform2D5 *wmt2D5);


void static WlzMakeAffine2D5With3PointsTrFn( WlzDVertex2 sr0, 
                                 WlzDVertex2 sr1, 
				 WlzDVertex2 sr2,
				 WlzDVertex3 targ0,
				 WlzDVertex3 targ1, 
                                 WlzDVertex3 targ2, 
			         double **Affine2D5With3PointsTrFun );

#ifdef DEBUGME
static WlzPlaneDomain *WlzPDomFromBBox(WlzIBox3  bBox, 
				       WlzErrorNum *dstErr);
#endif

void WlzGet2D5Transform( const WlzMeshTransform3D   *wmt3D,
                      int   cutPosition,
                      WlzMeshTransform2D5        *wmt2D5,
		      WlzErrorNum                *errNum
                    );

void   WlzInitializeWmt2D5(WlzMeshTransform2D5  *wmt2D5);

#ifdef DEBUGME
void   static WlzIsoIntersectWithTetrahadronIndexAnDETC(  double                      zConst, 
                                                const WlzMeshTransform3D   *wmt3D, 
					        WlzMeshTransform2D5        *wmt2D5
					     );
#endif
#ifdef DEBUGME
void  static validityOfMinAndMax( const int   nxmax, const int   nymax, const int   nzmax, 
                           const double xmin, const double xmax, 
			   const double ymin, const double ymax,
			   const double zmin, const double zmax  );
#endif

#ifdef DEBUGME
void  static ElementsGenerator(const int nxmax, const int nymax, const int nzmax,
                        WlzMeshTransform3D   *wmt3D, 
			const int *intersectionTable);
#endif
#ifdef DEBUGME
void static GetIntersectionTable( const int nxmax, const int nymax, const int nzmax,
                             const double xmin, const double ymin, const double zmin,
                             const double xdelta, const double ydelta, const double zdelta,
			     WlzObject *wObjC,
			     int *intersectionTable 
			   );

#endif

WlzErrorNum  read_WlzTiePoints( FILE *fp, int *nTiePP, 
                     WlzDVertex3 **vxVec0, 
		     WlzDVertex3 **vxVec1, 
		     const int  ReadDisplacement  );

#ifdef DEBUGME
void static AllocateMemForElemAndNodes(WlzMeshTransform3D   *wmt3D,
                                const int nxynumber,  
			        const int nxnumber,
			        const double xmin,
			        const double ymin,
			        const double zmin,
			        const double xdelta,
			        const double ydelta,
			        const double zdelta
                                );
#endif

WlzMeshTransform3D *WlzTetrahedronMeshFromObj( WlzObject *wObjC, const WlzDBox3 bBoxS,
				       const int numOfElemAlonX,
				       const int numOfElemAlonY,
				       const int numOfElemAlonZ,
                                       WlzErrorNum *errNums
                                  );


WlzMeshTransform2D5  *Wlz2D5TransformFromCut3Dmesh(double zConst, 
                                                   WlzMeshTransform3D *wmt3D,
						   WlzErrorNum *disErr);				  
int static IsNeighbour(int *n, int *np, int *it);

/*! 
* \ingroup      WlzMesh
* \brief        Divide a cuboid into six tetrahedral elements.
*
* \return       Woolz error.
* \param    neighbourLeft	Left neighbour index
* \param    neighbourRight	right neighbour index
* \param    neighbourBack	back neighbour index
* \param    neighbourFront	fromt neighbour index
* \param    neighbourDown	down neighbour index
* \param    neighbourUp	Up neighbour index
* \param    nxmax	X size
* \param    nymax	Y size
* \param    nzmax	Z size
* \param    nx	x origin
* \param    ny	y origin
* \param    nz	z origin
* \param    indexOfNextFirst	Index for new elements
* \param    elements	Element array.
* \par      Detail
* \verbatim
                                                                      
                          @4----------------@5                        
                         /|                /|                        
                        / |               / |                         
                       /  |              /  |                         
                      /   |             /   |                         
                     /    |            /    |                         
                    /     |           /     |                         
                   /      |          /      |                         
                  /       |         /       |                         
                 @0----------------@1       |                         
                 |        |        |        |                         
                 |        @7-------|--------@6                                 
                 |       /         |       /                                 
                 |      /          |      /                                  
                 |     /           |     /                                  
                 |    /            |    /                                  
                 |   /             |   /                                  
                 |  /              |  /                                  
                 | /               | /                                  
                 |/                |/                                  
                 @3----------------@2                                 

		The tetrahedra are are assigned the following indicies:

		  tetrahedron index  	cube vertex indicies
		   0       		1, 7, 0, 3
		   1       		1, 7, 3, 2
		   2       		1, 7, 2, 6
		   3       		1, 7, 6, 5 
		   4        		1, 7, 5, 4 
		   5        		1, 7, 4, 0


            we will use the following order of surfaces for tetrahedron:
                   If the vertex indicies of a tetrahedron is in the
                   following order:
                                       a, b, c, d
                   then the surface order is:

                                       a-b-c
                                       b-c-d
                                       c-d-a
                                       d-a-b


 The relationship with orther tetrahedron is described by its neighbours
                   neighbours is defined here as a surface neighbour and
                   we assume that every surface has only on neighbour which
                   means a tetrahedron has at most four neighbours.


 Input:     the eight vertexes, 
            the indexOfTheNextFirst of the first tetrahedron of
                                    the next six new tetrahedrons. 
            the neighbours:


                   <-          neighbourLeft,  neighbourRight       ->   x,
                               neighbourUp,    neighbourDown,       ->   y
                               neighbourFront, neighbourBack,       ->   z

                                                   1 means there is a neighbour
                                                   0 no-neighbour
                
               z
             / 
            /                                  
            -------->  x
           |
           |
           
           y

\endverbatim
* \par      Source:
*                Wlz3DWarpMQ_S.c
*/
WlzErrorNum WlzTetrahedronProducerFromCube(
  int neighbourLeft,
  int neighbourRight,
  int neighbourBack,
  int neighbourFront,
  int neighbourDown,
  int neighbourUp,
  int nxmax,
  int nymax,
  int nzmax,
  int nx,
  int ny,
  int nz,
  int indexOfNextFirst,
  WlzMeshElem3D *elements)
{
  int	     /*	tI0, */
		i,
		Node0,
		Node1,
		Node2,
		Node3,
		Node4,
		Node5,
		Node6,
    Node7;
  /*	intersect; */
/* Cube's values relative to the iso-value */
  /* double 	cVal[8], */	 
/* Tetrahedron values */ 
    /* 		tVal[4];	*/		       
  /* WlzDVertex3	tPos[4]; */
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Tetrahedron to cube vertex look up table */
  /*!const int	tVxLUT[6][4] =  
   *{
   * {1, 7, 0, 3}, {1, 7, 3, 2}, {1, 7, 2, 6},
   * {1, 7, 6, 5}, {1, 7, 5, 4}, {1, 7, 4, 0}
   *};
   */
 /* Cube positions, order is {vtX, vtY, vtZ} */
  /*!
   * const WlzDVertex3 cPos[8] =	 
   *{
   *  {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
   *  {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}
   *};
   */

   
  /* Calculate the index of the in this cube  counting from 0 */ 
        
	Node0 =       (nz - 1) * (nxmax+1) * (nymax+1)     /* Number of Nodes on previous planes */
	            + (ny - 1) * (nxmax+1)                 /* Number of Nodes on previous lines in current plane */
		    + (nx - 1);
	Node1 = Node0 + 1;
	Node3 = Node1 + nxmax;
        Node2 = Node3 + 1;
	Node4 = Node0 + (nxmax+1)*(nxmax+1);
	Node5 = Node4 + 1;
	Node7 = Node5 + nxmax;
	Node6 = Node7 + 1;
  						 
  /* move to the correct postion */
  for(i=0; i<indexOfNextFirst; i++)
  {
    elements++;
  }

  /* Generate the first tetrahedron */

    elements->idx =      indexOfNextFirst;

    elements->nodes[0] = Node1;
    elements->nodes[1] = Node7;
    elements->nodes[2] = Node0;
    elements->nodes[3] = Node3;

    elements->neighbours[0] = indexOfNextFirst + 5;
    elements->neighbours[1] = -1; 
        if(neighbourLeft  == 1) elements->neighbours[1] = indexOfNextFirst - 6 + 3;
    elements->neighbours[2] = -1; 
        if(neighbourFront == 1) elements->neighbours[2] = indexOfNextFirst - 6*nx*ny + 4;
    elements->neighbours[3] = indexOfNextFirst + 1;

  /* Generate the second tetrahedron */

    elements++;   
    
    elements->idx = indexOfNextFirst + 1;

    elements->nodes[0] = Node1;
    elements->nodes[1] = Node7;
    elements->nodes[2] = Node3;
    elements->nodes[3] = Node2;

    elements->neighbours[0] = indexOfNextFirst +  0;
    elements->neighbours[1] = -1; 
        if(neighbourDown  == 1) elements->neighbours[1] = indexOfNextFirst + 6*nx + 5;
    elements->neighbours[2] = -1; 
        if(neighbourFront == 1) elements->neighbours[2] = indexOfNextFirst - 6*nx*ny + 3;
    elements->neighbours[3] = indexOfNextFirst +  2;

   /* Generate the third tetrahedron */

    elements++;   
    
    elements->idx = indexOfNextFirst + 2;

    elements->nodes[0] = Node1;
    elements->nodes[1] = Node7;
    elements->nodes[2] = Node2;
    elements->nodes[3] = Node6;

    elements->neighbours[0] = indexOfNextFirst +  1;
    elements->neighbours[1] = -1;
        if(neighbourDown  == 1) elements->neighbours[1] = indexOfNextFirst + 6*nx + 4;
    elements->neighbours[2] = -1;
        if(neighbourRight == 1) elements->neighbours[2] = indexOfNextFirst + 6 + 0;
    elements->neighbours[3] = indexOfNextFirst +  3;

   /* Generate the fourth tetrahedron */

    elements++;   
    
    elements->idx = indexOfNextFirst + 3;

    elements->nodes[0] = Node1;
    elements->nodes[1] = Node7;
    elements->nodes[2] = Node6;
    elements->nodes[3] = Node5;

    elements->neighbours[0] = indexOfNextFirst +  2;
    elements->neighbours[1] = -1;
        if(neighbourBack == 1) elements->neighbours[1] = indexOfNextFirst + 6*nx*ny + 1;
    elements->neighbours[2] = -1;
        if(neighbourFront == 1) elements->neighbours[2] = indexOfNextFirst + 6 + 5;
    elements->neighbours[3] = indexOfNextFirst +  4;

   /* Generate the fifth tetrahedron */

    elements++;   
    
    elements->idx = indexOfNextFirst + 4;

    elements->nodes[0] = Node1;
    elements->nodes[1] = Node7;
    elements->nodes[2] = Node5;
    elements->nodes[3] = Node4;

    elements->neighbours[0] =  indexOfNextFirst + 3;
    elements->neighbours[1] = -1; 
        if(neighbourBack == 1) elements->neighbours[1] = indexOfNextFirst + 6*nx*ny + 0;
    elements->neighbours[2] = -1; 
        if(neighbourUp == 1)   elements->neighbours[2] = indexOfNextFirst - 6*nx + 2;
    elements->neighbours[3] =  indexOfNextFirst + 5;

   /* Generate the last tetrahedron */

    elements++;   
    
    elements->idx = indexOfNextFirst + 5;

    elements->nodes[0] = Node1;
    elements->nodes[1] = Node7;
    elements->nodes[2] = Node4;
    elements->nodes[3] = Node0;

    elements->neighbours[0] = indexOfNextFirst + 5;
    elements->neighbours[1] = -1;
        if(neighbourLeft == 1) elements->neighbours[1] = indexOfNextFirst - 6 + 3;
    elements->neighbours[2] = -1;
        if(neighbourUp == 1)   elements->neighbours[2] = indexOfNextFirst - 6*nx + 1;
    elements->neighbours[3] = indexOfNextFirst + 0;

  return(errNum);
}

/*!
*  \return int nIntersect number of tetrahedrons whose body was cutted by the given plane.
*  \ingroup WlzMesh
*  \brief   
*          Get the index of the tetrahedrons which intersect
*          with a giving z = constant plane. 
*    
*  \param       zConst        A z = const plane used to cut the tetrahedron mesh.
*  \param       wmt3D         the tetrahedron mesh.
*   
*         the table containing the number of tetrahedron whose body intersected
*            with the give plane.
*  \param   intersectIndex    stores the intex of tetrahedron whose body been cutted by the given plane.
*  \param   linkList          linkList[i][j] stores the number of points cutted to the tetrahedra of intersectIndex[i] 
*                             in j = 0. And store the sequence number in j=0,1,...intersectIndex[i][0]-1;
*  \param   planepoints       store the cuting points sequentially according the cuting order.
*  \param   noRedundancyCutingNum store the number of cuting points.
*  \param   numTotalTrangularElem  number of total Trangular Elements.
*  \author    J. Rao, R. Baldock and B. Hill
*  \par   Detail
*  \verbatim

          x
         / \
        /   \     x
       /     \   /
      /       \ /
     x---------x        

\endverbatim 
                                          18/10/2001 J.Rao,  try to remove redundancy in the storage and
					  make it easier to connet with Bill's 2D data structure.
					  It seems OK now.
           
    

*/
int WlzIsoIntersectWithTetrahadronIndex(  double                      zConst, 
                                          const WlzMeshTransform3D   *wmt3D, 
					  int                        *intersectIndex,
                                          WlzDVertex3                *planepoints, 
					  int                       **linkList, 
					  int                        *noRedundancyCutingNum,
					  int                        *numTotalTrangularElem)
{
  int      i, j0, j1, j2, j3, k, k1, itemp;/* kmin, kmax */
  int      nIntersect = 0;
  int      nUp;
  int      iup, idown, icount, numberOfTriangular=0, isnumber=0;
  int      n0up, n1up, n2up, n3up;
  double   z0, z1, z2, z3;
  /* double   zt; */

 /* WlzDVertex3 planepoints[3000]; */  /* store the cutting points */
 /* int         linkList[3000][5];   */   /* store the linked list
                                      the first index linked with
				      the index of cutting tetrahedron 
				      The second one linked with 
				      the number of points cuttied to the tetrahedron */
  WlzDVertex3 wv3up[4];
  WlzDVertex3 wv3down[4];
  WlzDVertex3 tempw;
  /* int       iidebug = 1; */
    int       idex;
  
  icount = 0;  /* used to count how many points cutted (no redundancy counting) */

  for (i = 0; i < wmt3D->nElem; i++)
  {   /*  get the z-position of the four points relative to the zConst 
          Remember that we are cutting the transformed mesh */
      /* Also, We only cut the mesh belong to the object! */
      /* exclude the tetrahedrons not belong to the object */

         j0    =  (wmt3D->elements+i )->nodes[0];
         j1    =  (wmt3D->elements+i )->nodes[1];
         j2    =  (wmt3D->elements+i )->nodes[2];
         j3    =  (wmt3D->elements+i )->nodes[3];
    
     if( (j0 != 0) || (j1 !=0 ) || (j2 !=0) || (j3 !=0 ) )
     {

         z0    =  (wmt3D->nodes+j0)->position.vtZ + (wmt3D->nodes+j0)->displacement.vtZ   - zConst;
         z1    =  (wmt3D->nodes+j1)->position.vtZ + (wmt3D->nodes+j1)->displacement.vtZ   - zConst;
         z2    =  (wmt3D->nodes+j2)->position.vtZ + (wmt3D->nodes+j2)->displacement.vtZ   - zConst;
         z3    =  (wmt3D->nodes+j3)->position.vtZ + (wmt3D->nodes+j3)->displacement.vtZ   - zConst;
         nUp   =  0;
         n0up  =  n1up = n2up = n3up = 0;
         /* printf("%lf  %lf  %lf  %lf\n", z0, z1, z2, z3); */

         /* judge whether there is a intersection with zConst-plane 
         only 3 cases if intersected

	 1) one   point up,  three points down 
	 2) two   point up,  two   points down 
	 3) three point up,  one   points down 
	
         Not intersected only two cases:
	 4) four points all up
	 5) four points all down
         */

         if( z0 > 0. )
         {
            n0up = 1;
            nUp++;
         }
     
         if( z1 > 0.  )
         {
            n1up = 1;
             nUp++;
         }
     
         if( z2 > 0.  )
         {
            n2up = 1;
             nUp++;
         }
     
         if( z3 > 0.  )
         {
            n3up = 1;
             nUp++;
         }
      
         /* Now the cutted cases  */
         if( (nUp != 0) && (nUp != 4) )
         {
           iup = idown = 0;
           /* divide the nodes of the cutted tetrahedron into up and down nodes then store them. */
           if(n0up == 1)
	   {
            wv3up[iup].vtX = (wmt3D->nodes+j0)->position.vtX + (wmt3D->nodes+j0)->displacement.vtX; 
            wv3up[iup].vtY = (wmt3D->nodes+j0)->position.vtY + (wmt3D->nodes+j0)->displacement.vtY; 
            wv3up[iup].vtZ = (wmt3D->nodes+j0)->position.vtZ + (wmt3D->nodes+j0)->displacement.vtZ;
	    iup++;
	   }
	   else if(n0up == 0)
	   {
            wv3down[idown].vtX = (wmt3D->nodes+j0)->position.vtX + (wmt3D->nodes+j0)->displacement.vtX; 
            wv3down[idown].vtY = (wmt3D->nodes+j0)->position.vtY + (wmt3D->nodes+j0)->displacement.vtY; 
            wv3down[idown].vtZ = (wmt3D->nodes+j0)->position.vtZ + (wmt3D->nodes+j0)->displacement.vtZ; 
	    idown++;
 	   }
	  
           if(n1up == 1)
	   {
            wv3up[iup].vtX = (wmt3D->nodes+j1)->position.vtX + (wmt3D->nodes+j1)->displacement.vtX; 
            wv3up[iup].vtY = (wmt3D->nodes+j1)->position.vtY + (wmt3D->nodes+j1)->displacement.vtY; 
            wv3up[iup].vtZ = (wmt3D->nodes+j1)->position.vtZ + (wmt3D->nodes+j1)->displacement.vtZ;
	    iup++;
	   }
	   else if(n1up == 0)
	   {
            wv3down[idown].vtX = (wmt3D->nodes+j1)->position.vtX + (wmt3D->nodes+j1)->displacement.vtX; 
            wv3down[idown].vtY = (wmt3D->nodes+j1)->position.vtY + (wmt3D->nodes+j1)->displacement.vtY; 
            wv3down[idown].vtZ = (wmt3D->nodes+j1)->position.vtZ + (wmt3D->nodes+j1)->displacement.vtZ; 
	    idown++;
 	   }
		
           if(n2up == 1)
	   {
            wv3up[iup].vtX = (wmt3D->nodes+j2)->position.vtX + (wmt3D->nodes+j2)->displacement.vtX; 
            wv3up[iup].vtY = (wmt3D->nodes+j2)->position.vtY + (wmt3D->nodes+j2)->displacement.vtY; 
            wv3up[iup].vtZ = (wmt3D->nodes+j2)->position.vtZ + (wmt3D->nodes+j2)->displacement.vtZ; 
	    iup++;
	   }
	   else if(n2up == 0)
	   {
            wv3down[idown].vtX = (wmt3D->nodes+j2)->position.vtX + (wmt3D->nodes+j2)->displacement.vtX; 
            wv3down[idown].vtY = (wmt3D->nodes+j2)->position.vtY + (wmt3D->nodes+j2)->displacement.vtY; 
            wv3down[idown].vtZ = (wmt3D->nodes+j2)->position.vtZ + (wmt3D->nodes+j2)->displacement.vtZ; 
	    idown++;
 	   }
		  
           if(n3up == 1)
	   {
            wv3up[iup].vtX = (wmt3D->nodes+j3)->position.vtX + (wmt3D->nodes+j3)->displacement.vtX; 
            wv3up[iup].vtY = (wmt3D->nodes+j3)->position.vtY + (wmt3D->nodes+j3)->displacement.vtY; 
            wv3up[iup].vtZ = (wmt3D->nodes+j3)->position.vtZ + (wmt3D->nodes+j3)->displacement.vtZ; 
	    iup++;
	   }
	   else if(n3up == 0)
	   {
            wv3down[idown].vtX = (wmt3D->nodes+j3)->position.vtX + (wmt3D->nodes+j3)->displacement.vtX; 
            wv3down[idown].vtY = (wmt3D->nodes+j3)->position.vtY + (wmt3D->nodes+j3)->displacement.vtY; 
            wv3down[idown].vtZ = (wmt3D->nodes+j3)->position.vtZ + (wmt3D->nodes+j3)->displacement.vtZ; 
	    idown++;
 	   }

           /* Now get the intersect points */
	   /* for output test, we always need an sequence to make a good output 
	     for vtk visualization 
	    
	     18/10/2001 but now we need to remove the redundancy storage */
	     tempw.vtZ = 0.;
           /*  ----- Now three points cut case -----  */
	   if( nUp == 1 || nUp == 3  )
	   {
	    isnumber = 1;
            for(k = 0; k < nUp; k++) 
	    { 
 
	        for( k1=0; k1< 4 - nUp; k1++ )
		{
		   if( wv3up[k].vtZ  - wv3down[k1].vtZ  == 0. )
		     { printf("the line in the z-pline");
		       exit(1); 
		     }
                    /* Get the cutting points */
                    tempw.vtX  =  wv3down[k1].vtX +   
		                ( wv3up[k].vtX - wv3down[k1].vtX ) * ( zConst - wv3down[k1].vtZ )
		               /( wv3up[k].vtZ - wv3down[k1].vtZ ); 
			
                    tempw.vtY  =  wv3down[k1].vtY  + 
		                ( wv3up[k].vtY - wv3down[k1].vtY ) * ( zConst - wv3down[k1].vtZ )
		               /( wv3up[k].vtZ - wv3down[k1].vtZ );
			      
	            /* Ask whether this point is already stored 
		       if not return idex < 0 otherwise return the
		       index of the existing point already stored */
		    idex = -1;
                    idex = storedPointIndex(tempw, planepoints, icount);
		 
                    if(idex < 0)
		    {
		     /* store the points */
                     (planepoints+icount)->vtX  =  tempw.vtX;
                     (planepoints+icount)->vtY  =  tempw.vtY; 
                     (planepoints+icount)->vtZ  =  zConst;
		   
	             /* record the index of nodes for this cut */
               	     *( *(linkList+nIntersect) + isnumber ) = icount;
                     isnumber++;
     	             icount++;  
		    }
		    else if(idex >0)
		    { /* the point is already stored before
		       we only need to store the index */
	             /* record the index of nodes for this cut */
               	     *( *(linkList+nIntersect) + isnumber ) = idex;
                     isnumber++;
		    }
	         }
	     }	
	     /* now make sure it is counterclockwise odrdered */
	     if( WlzGeomTriangle2SnArea3(  planepoints,  *( *(linkList+nIntersect) + 1  ),
	                                      *( *(linkList+nIntersect) + 2),  
		                              *( *(linkList+nIntersect) + 3)
                             ) < 0.       )
	     { /* change the position of the last two nodes data structure
	         to make it CCW order */
               itemp =  *( *(linkList+nIntersect) + 2);
               *( *(linkList+nIntersect) + 2) = *( *(linkList+nIntersect) + 3);
	       *( *(linkList+nIntersect) + 3) = itemp;
	     }
	  }

	  /* ------ four points cutting case ------- */
	  else if(nUp == 2)
	  {
	    isnumber = 1;
	    /* k = 0 */
	    k = 0; 
	    for( k1=0; k1< 2; k1++ )
	    {
	       if( wv3up[k].vtZ  - wv3down[k1].vtZ  == 0. )
		{ printf("the line in the z-pline");
		       exit(1); 
		}
                /* Get the cutting points */
                tempw.vtX  =  wv3down[k1].vtX +   
		              ( wv3up[k].vtX - wv3down[k1].vtX ) * ( zConst - wv3down[k1].vtZ )
		             /( wv3up[k].vtZ - wv3down[k1].vtZ ); 
			
                tempw.vtY  =  wv3down[k1].vtY  + 
		              ( wv3up[k].vtY - wv3down[k1].vtY ) * ( zConst - wv3down[k1].vtZ )
		             /( wv3up[k].vtZ - wv3down[k1].vtZ ); 

	        /* Ask whether this point is already stored 
		   if not return idex < 0 otherwise return the
		   index of the existing point already stored */
		idex = -1;
                idex = storedPointIndex(tempw, planepoints, icount);
                if(idex < 0)
		{
		   /* store the points */
                   (planepoints+icount)->vtX  =  tempw.vtX;
                   (planepoints+icount)->vtY  =  tempw.vtY; 
                   (planepoints+icount)->vtZ  =  zConst;
		   
	           /* record the index of nodes for this cut */
               	   *( *(linkList+nIntersect) + isnumber ) = icount;
                   isnumber++;
     	           icount++;  
		 }
		 else if(idex >0)
		 { /* the point is already stored before
		       we only need to store the index */
	           /* record the index of nodes for this cut */
               	   *( *(linkList+nIntersect) + isnumber ) = idex;
                   isnumber++;
		 }
	      }
		
	      /* k = 1 */
	      k = 1; 
	      for( k1= 1;  k1 >= 0; k1-- )
	      {
		   if( wv3up[k].vtZ  - wv3down[k1].vtZ  == 0. )
		     { printf("the line in the z-pline");
		       exit(1); 
		     }

                     tempw.vtX  =  wv3down[k1].vtX +   
		                 ( wv3up[k].vtX - wv3down[k1].vtX ) * ( zConst - wv3down[k1].vtZ )
		                /( wv3up[k].vtZ - wv3down[k1].vtZ ); 
			
                     tempw.vtY  =  wv3down[k1].vtY  + 
		                 ( wv3up[k].vtY - wv3down[k1].vtY ) * ( zConst - wv3down[k1].vtZ )
		                /( wv3up[k].vtZ - wv3down[k1].vtZ ); 

	             /* Ask whether this point is already stored 
		        if not return idex < 0 otherwise return the
		        index of the existing point already stored */
		     idex = -1;
                     idex = storedPointIndex(tempw, planepoints, icount);
                     if(idex < 0)
		     {
		       /* store the points */
                       (planepoints+icount)->vtX  =  tempw.vtX;
                       (planepoints+icount)->vtY  =  tempw.vtY; 
                       (planepoints+icount)->vtZ  =  zConst;
		   
	               /* record the index of nodes for this cut */
               	       *( *(linkList+nIntersect) + isnumber ) = icount;
                       isnumber++;
     	               icount++;  
		      }
		      else if(idex >0)
		      { /* the point is already stored before
		           we only need to store the index */
	                /* record the index of nodes for this cut */
               	        *( *(linkList+nIntersect) + isnumber ) = idex;
                        isnumber++;
		      }
		}
	        /* now make sure that this quarilateral divided into two CCW 
	            triangles */
	        if( WlzGeomTriangle2SnArea3(  planepoints,  *( *(linkList+nIntersect) + 1 ),
	                                          *( *(linkList+nIntersect) + 2),  
		                                  *( *(linkList+nIntersect) + 3)
                                   ) < 0.       )
	        { /* change the position of the second and the last nodes */
	              itemp                          = *( *(linkList+nIntersect) + 2);
	              *( *(linkList+nIntersect) + 2) = *( *(linkList+nIntersect) + 4);
	              *( *(linkList+nIntersect) + 4) = itemp;
	        }

	    }
	 
	    /* regester the intersect tetrahedral element  */ 
	     *(intersectIndex+nIntersect) = i;

	    /* record how many points of this cut */
	    if( (nUp == 2 ) )  
	          *( *( linkList + nIntersect ) + 0 ) = 4;
	    else 
	          *( *( linkList + nIntersect ) + 0 ) = 3;

	    /* move to the next tetrahedron element */
            nIntersect++;
       }

     }
  }
  /* get the number of triangulars */
   numberOfTriangular = 0;
  for(k1=0; k1< nIntersect; k1++)
  {   if(  *(*(linkList+k1) +0)  == 3 )
         numberOfTriangular += 1;
       if( *(*(linkList+k1) +0)   == 4 )
       { /* output by divide it into two triangles */
         numberOfTriangular += 2;
       }
  }
  /* Now the number of points cutted (no redundancy counting) */
  *noRedundancyCutingNum  = icount;
  /* Now the number of triangular elments */
  *numTotalTrangularElem  = numberOfTriangular;
  /* we can just assign the value to the wmt2D5 mesh here?  */


  
  return nIntersect;  /* this is the number of cut tetrahedrons */
}

/*!
*  \return none
*  \ingroup  WlzTransform
*  \brief 
*
*    Get 3D Affine transformation from the nodes of
*    a source tetrahedron and the nodes from a target tetrahedron.
*  
*  \param   sr0    nodes of source tetrahedron
*  \param   sr1    nodes of source tetrahedron
*  \param   sr2    nodes of source tetrahedron
*  \param   sr3    nodes of source tetrahedron
*  \param   targ0  nodes of target tetrahedron
*  \param   targ1  nodes of target tetrahedron
*  \param   targ2  nodes of target tetrahedron
*  \param   targ3  nodes of target tetrahedron
*  \param   Affine3D4pointsTrFun  contains the affine transformation parameters 
*
*  \author:       J. Rao, R. Baldock and B. Hill
*  \par  Detail
*  \verbatim
  A =:
       a11  a12  a13  a14
       a21  a22  a23  a24
       a31  a32  a33  a34
        0    0    0    1
           
 \endverbatim
*/
void WlzMakeAffine3D4pointsTrFn( WlzDVertex3 sr0, 
                                 WlzDVertex3 sr1, 
				 WlzDVertex3 sr2,
                                 WlzDVertex3 sr3,  
				 WlzDVertex3 targ0,
				 WlzDVertex3 targ1, 
                                 WlzDVertex3 targ2, 
				 WlzDVertex3 targ3, 
			         double **Affine3D4pointsTrFun )
{
  int            i, j, idN;
  double	*bV = NULL,
  		*wV = NULL;
  double       **aA,
               **vA;
  double        wMax;	       
  int           nSys = 4;
  int           thresh;
  const double  tol  = 1.0e-06;
  AlgMatrix	aM,
  		vM;
  WlzErrorNum	errNum = WLZ_ERR_NONE;


  aM.core = NULL;
  vM.core = NULL;
  /*------- allocate memory --------*/
  if (
      ((wV = (double *)AlcCalloc(sizeof(double), nSys))  == NULL) ||
      ((bV = (double *)AlcMalloc(sizeof(double) * nSys)) == NULL) ||
      ((aM.rect = AlgMatrixRectNew(nSys, nSys, NULL)) == NULL) ||
      ((vM.rect = AlgMatrixRectNew(nSys, nSys, NULL)) == NULL))
  {
      errNum = WLZ_ERR_MEM_ALLOC;
  }
  
  if(errNum == WLZ_ERR_NONE)
  {
    aA = aM.rect->array;
    vA = vM.rect->array;
    /* fill the a-matrix  where a is:
   
      x0 y0 z0  1
      x1 y1 z1  1
      x2 y2 z2  1
      x3 y3 z3  1
   
    where (x1, y1, z1) is the source node one and etc.
   
    */

    *(*(aA + 0 ) + 0 ) = sr0.vtX;
    *(*(aA + 0 ) + 1 ) = sr0.vtY;
    *(*(aA + 0 ) + 2 ) = sr0.vtZ;
    *(*(aA + 0 ) + 3 ) = 1.0;

    *(*(aA + 1 ) + 0 ) = sr1.vtX;
    *(*(aA + 1 ) + 1 ) = sr1.vtY;
    *(*(aA + 1 ) + 2 ) = sr1.vtZ;
    *(*(aA + 1 ) + 3 ) = 1.0;

    *(*(aA + 2 ) + 0 ) = sr2.vtX;
    *(*(aA + 2 ) + 1 ) = sr2.vtY;
    *(*(aA + 2 ) + 2 ) = sr2.vtZ;
    *(*(aA + 2 ) + 3 ) = 1.0;

    *(*(aA + 3 ) + 0 ) = sr3.vtX;
    *(*(aA + 3 ) + 1 ) = sr3.vtY;
    *(*(aA + 3 ) + 2 ) = sr3.vtZ;
    *(*(aA + 3 ) + 3 ) = 1.0;

    /* get a11 a12 a13 and a14 */
    
    /* fill the u */
    *(bV + 0 ) = targ0.vtX - sr0.vtX;
    *(bV + 1 ) = targ1.vtX - sr1.vtX;
    *(bV + 2 ) = targ2.vtX - sr2.vtX;
    *(bV + 3 ) = targ3.vtX - sr3.vtX;

    /* */
    for(i=0;i<4; i++)
    {
      *(wV+i) = 0.0;
      for(j=0;j<4;j++)
      {
        *(*(vA+i)+j) = 0.0;
      }
    }
   
    /* Perform singular value decomposition of matrix a. */
   
    errNum = WlzErrorFromAlg(AlgMatrixSVDecomp(aM, wV, vM));
    /* */
    /*  
      for(i=0;i<4;i++)
      {
        printf("%f\n",*(wV+i));
        for(j=0;j<4;j++)
        {
          printf("%f   ",*(*(vA+i) +j) );
          printf("\n");
        }
       }
    */  
 
    if(errNum != WLZ_ERR_NONE)
    {
       printf("failed to make decompositions ");
       exit(1);
    }
    {
      /* Edit the singular values. */
      wMax = 0.0;
      for(idN = 0; idN < nSys; ++idN)
      {
        if(*(wV + idN) > wMax)
        {
	  wMax = *(wV + idN);
        }
      }
      thresh = tol * wMax;
      for(idN = 0; idN < nSys; ++idN)
      {
        if(*(wV + idN) < thresh)
        {
	  *(wV + idN) = 0.0;
        }
      }
      /* Solve for the X conponents  */
      errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
    }

    if(errNum != WLZ_ERR_NONE)
    {
       printf("failed to get Affine transform for x-component");
       exit(1);
    }
    {
      /* Store the a11, a12, a13, a14 */
       *(*(Affine3D4pointsTrFun + 0 )+ 0 ) = *(bV + 0);
       *(*(Affine3D4pointsTrFun + 0 )+ 1 ) = *(bV + 1);
       *(*(Affine3D4pointsTrFun + 0 )+ 2 ) = *(bV + 2);
       *(*(Affine3D4pointsTrFun + 0 )+ 3 ) = *(bV + 3);
     
      /* Now set for y  conponents */
       *(bV + 0 ) = targ0.vtY - sr0.vtY;
       *(bV + 1 ) = targ1.vtY - sr1.vtY;
       *(bV + 2 ) = targ2.vtY - sr2.vtY;
       *(bV + 3 ) = targ3.vtY - sr3.vtY;
   
      /* Solve for the Y conponents */
      errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
    }

    if(errNum != WLZ_ERR_NONE)
    {
       printf("failed to get Affine transform for y-component");
       exit(1);
    }
 
    {
      /* Store the a21, a22, a23, a24 */
       *(*(Affine3D4pointsTrFun + 1 )+ 0 ) = *(bV + 0);
       *(*(Affine3D4pointsTrFun + 1 )+ 1 ) = *(bV + 1);
       *(*(Affine3D4pointsTrFun + 1 )+ 2 ) = *(bV + 2);
       *(*(Affine3D4pointsTrFun + 1 )+ 3 ) = *(bV + 3);
     
      /* Now set for z  conponents */
       *(bV + 0 ) = targ0.vtZ - sr0.vtZ;
       *(bV + 1 ) = targ1.vtZ - sr1.vtZ;
       *(bV + 2 ) = targ2.vtZ - sr2.vtZ;
       *(bV + 3 ) = targ3.vtZ - sr3.vtZ;
   
      /* Solve for the Z conponents */
       errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
    }

    if(errNum != WLZ_ERR_NONE)
    {
       printf("failed to get Affine transform for z-component");
       exit(1);
    }
    {
      /* Store the a31, a32, a33, a34 */
       *(*(Affine3D4pointsTrFun + 2 )+ 0 ) = *(bV + 0);
       *(*(Affine3D4pointsTrFun + 2 )+ 1 ) = *(bV + 1);
       *(*(Affine3D4pointsTrFun + 2 )+ 2 ) = *(bV + 2);
       *(*(Affine3D4pointsTrFun + 2 )+ 3 ) = *(bV + 3);
    }
    AlcFree(bV);
    AlcFree(wV);
    AlgMatrixFree(aM);
    AlgMatrixFree(vM);
  }
 
}

/*!
* \return   None
* \ingroup  WlzGeometry
* \brief    output the orginal mesh.
* \param    fp               pointer pointing to a specific file.
* \param    wmt3D            mesh transform.
* \author:       J. Rao, R. Baldock and B. Hill
*/
void WlzEffWriteMeshTransform3DWithoutDisplacementVTK(FILE *fp,  
               WlzMeshTransform3D *wmt3D){
	       
  /* ------- output as vtk form for check ------- */
 int     i, it; 
  if(wmt3D == NULL )
  {
    printf("No MeshTransform3D.\n");
    exit(1);
  }
   
  fprintf(fp, "# vtk DataFile Version 1.0\n");
  fprintf(fp, "WlzGeoModel test output\n");
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET POLYDATA\n");
  fprintf(fp, "%s %d %s\n","POINTS",wmt3D->nNodes,"  float");
  /* 3D points */
  for(i=0; i<wmt3D->nNodes; i++)
    {
       fprintf(fp, "%lf  %lf  %f\n",
                (wmt3D->nodes+i)->position.vtX,(wmt3D->nodes+i)->position.vtY,(wmt3D->nodes+i)->position.vtZ );
    }
  /* Polygons link */

  /* only output the linked list */
  it = 0;
  /* get the number of elements for output (the tetrahedron element) */
   for(i=0; i<wmt3D->nElem; i++)
   { 
    if(    ( (wmt3D->elements+i)->nodes[0]  > 0 ) || 
           ( (wmt3D->elements+i)->nodes[1]  > 0 ) ||
           ( (wmt3D->elements+i)->nodes[2]  > 0 ) ||
	   ( (wmt3D->elements+i)->nodes[3]  > 0 )  
      )
    it++; 
   }
 
   /* there are four triangles for each tetrahedron.  */
   fprintf(fp, "%s %d %ld\n","POLYGONS", it*4, (long )it*4L*sizeof(float) ); 
  
   for(i=0; i<wmt3D->nElem; i++)
   { /*
      if(   ( (wme3D+i)->nodes[1]  != 0 ) || ( (wme3D+i)->nodes[2] !=0   ) || ( (wme3D+i)->nodes[3] != 0 )  )
      */
      if(   ( (wmt3D->elements+i)->nodes[0]  > 0  ) || 
            ( (wmt3D->elements+i)->nodes[1]  > 0  ) ||
            ( (wmt3D->elements+i)->nodes[2]  > 0  ) ||
            ( (wmt3D->elements+i)->nodes[3]  > 0  )  )
     { 

       fprintf(fp, "%d  %d  %d  %d\n", 3, (wmt3D->elements+i)->nodes[0],
	       (wmt3D->elements+i)->nodes[1], 
	       (wmt3D->elements+i)->nodes[2] );
				    
        fprintf(fp, "%d  %d  %d  %d\n", 3, (wmt3D->elements+i)->nodes[1],
                                          (wmt3D->elements+i)->nodes[2], 
				          (wmt3D->elements+i)->nodes[3] );
				    
        fprintf(fp, "%d  %d  %d  %d\n",3,(wmt3D->elements+i)->nodes[2],
                                          (wmt3D->elements+i)->nodes[3], 
				          (wmt3D->elements+i)->nodes[0] );

        fprintf(fp, "%d  %d  %d  %d\n",3, (wmt3D->elements+i)->nodes[3],
                                          (wmt3D->elements+i)->nodes[0], 
				          (wmt3D->elements+i)->nodes[1] );

     }
   }
}

/*!
* \return   None
* \ingroup  WlzGeometry
* \brief    output the transformed mesh.
* \param    fp               pointer pointing to a specific file.
* \param    wmt3D            mesh transform.
* \author:       J. Rao, R. Baldock and B. Hill
*/
void WlzEffWriteMeshTransform3DWithDisplacementVTK(FILE *fp, 
               WlzMeshTransform3D *wmt3D){
	       
  /* ------- output as vtk form for check ------- */
 int     i, it; 
  if(wmt3D == NULL )
  {
    printf("No MeshTransform3D .\n");
    exit(1);
  }
   
  fprintf(fp, "# vtk DataFile Version 1.0\n");
  fprintf(fp, "WlzGeoModel test output\n");
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET POLYDATA\n");
  fprintf(fp, "%s %d %s\n","POINTS",wmt3D->nNodes,"  float");
  /* 3D points */
  for(i=0; i<wmt3D->nNodes; i++)
    {
       fprintf(fp, "%lf  %lf  %f\n",
                (wmt3D->nodes+i)->position.vtX + (wmt3D->nodes+i)->displacement.vtX, 
                (wmt3D->nodes+i)->position.vtY + (wmt3D->nodes+i)->displacement.vtY, 
                (wmt3D->nodes+i)->position.vtZ + (wmt3D->nodes+i)->displacement.vtZ); 
    }
  /* Polygons link */

  /* only output the linked list */
  it = 0;
  /* get the number of elements for output (the tetrahedron element) */
   for(i=0; i<wmt3D->nElem; i++)
   { 
    if(    ( (wmt3D->elements+i)->nodes[0]  > 0 ) || 
           ( (wmt3D->elements+i)->nodes[1]  > 0 ) ||
           ( (wmt3D->elements+i)->nodes[2]  > 0 ) ||
	   ( (wmt3D->elements+i)->nodes[3]  > 0 )  
      )
    it++; 
   }
 
   /* there are four triangles for each tetrahedron.  */
   fprintf(fp, "%s %d %ld\n","POLYGONS",it*4, (long)it*4L*sizeof(float) ); 
  
   for(i=0; i<wmt3D->nElem; i++)
   { /*
      if(   ( (wme3D+i)->nodes[1]  != 0 ) || ( (wme3D+i)->nodes[2] !=0   ) || ( (wme3D+i)->nodes[3] != 0 )  )
      */
      if(   ( (wmt3D->elements+i)->nodes[0]  > 0 ) || 
            ( (wmt3D->elements+i)->nodes[1]  > 0 ) ||
            ( (wmt3D->elements+i)->nodes[2]  > 0 ) ||
            ( (wmt3D->elements+i)->nodes[3]  > 0 )  )
     { 

        fprintf(fp, "%d  %d  %d  %d\n", 3, (wmt3D->elements+i)->nodes[0],
                                          (wmt3D->elements+i)->nodes[1], 
				          (wmt3D->elements+i)->nodes[2] );
				    
        fprintf(fp, "%d  %d  %d  %d\n", 3, (wmt3D->elements+i)->nodes[1],
                                          (wmt3D->elements+i)->nodes[2], 
				          (wmt3D->elements+i)->nodes[3] );
				    
        fprintf(fp, "%d  %d  %d  %d\n",3,  (wmt3D->elements+i)->nodes[2],
                                          (wmt3D->elements+i)->nodes[3], 
				          (wmt3D->elements+i)->nodes[0] );

        fprintf(fp, "%d  %d  %d  %d\n",3,  (wmt3D->elements+i)->nodes[3],
                                          (wmt3D->elements+i)->nodes[0], 
				          (wmt3D->elements+i)->nodes[1] );

     }
   }
}

/*!
* \return   None
* \ingroup  WlzGeometry
* \brief    output the cut plane in VTK format to a file.
* \param    fp               pointer pointing to a specific file.
* \param    wmt3D            mesh transform.
* \param    intersectIndex   The index indicating which tetrahedron cut by the plane.
* \param    nIntersect       number of intersect with the tetrahedron.
* \author       J. Rao, R. Baldock and B. Hill
*/
#ifdef DEBUGME
void static  WlzEffWriteCutplanTetrahedronVTK(FILE *fp, 
                     const WlzMeshTransform3D *wmt3D,  int *intersectIndex, int nIntersect )
{
  int i, j, it;
  if( wmt3D == NULL )
  {
    printf("No 3D TransformMesh\n");
    exit(1);
  }
  if(intersectIndex == NULL)
  {
    printf("No IntersectIndex \n");
    exit(1);
  }
  /* ------- output as vtk form for check ------- */
  fprintf(fp, "# vtk DataFile Version 1.0\n");
  fprintf(fp, "WlzGeoModel test output\n");
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET POLYDATA\n");
  fprintf(fp, "%s %d %s\n","POINTS", wmt3D->nNodes, "  float");
  
  /* 3D points */
  
  for(i=0; i<wmt3D->nNodes; i++)
  {
       fprintf(fp, "%lf  %lf  %f\n", 
                   (wmt3D->nodes+i)->position.vtX + (wmt3D->nodes+i)->displacement.vtX, 
		   (wmt3D->nodes+i)->position.vtY + (wmt3D->nodes+i)->displacement.vtY,
		   (wmt3D->nodes+i)->position.vtZ + (wmt3D->nodes+i)->displacement.vtZ );
  }

  /* Polygons link */

  /* only output the linked list */
  it = 0;
  for(i = 0; i<nIntersect; i++)
  {
     j = *(intersectIndex+i);
    if(   (  (  wmt3D->elements+j )->nodes[0]  < 0 )    || 
             ( (wmt3D->elements+j )->nodes[1]  < 0 )    || 
             ( (wmt3D->elements+j )->nodes[2]  < 0 )    || 
	     ( (wmt3D->elements+j )->nodes[3]  < 0 )  ) 
     it++; 
  }
   
  fprintf(fp,  "%s %d %d\n","POLYGONS",it*4, it*4*sizeof(float) ); 
  
  for(i=0; i<nIntersect; i++)
  { 
      j = *(intersectIndex+i);
      if(   ( (wmt3D->elements+j )->nodes[0]  < 0 )  || 
            ( (wmt3D->elements+j )->nodes[1]  < 0 )  || 
            ( (wmt3D->elements+j )->nodes[2]  < 0 )  || 
	    ( (wmt3D->elements+j )->nodes[3]  < 0 )    )
     { 

        fprintf(fp, "%d  %d  %d  %d\n",3,(wmt3D->elements+j)->nodes[0],
                                             (wmt3D->elements+j)->nodes[1], 
				             (wmt3D->elements+j)->nodes[2] );
				    
        fprintf(fp, "%d  %d  %d  %d\n",3,(wmt3D->elements+j)->nodes[1],
                                             (wmt3D->elements+j)->nodes[2], 
				             (wmt3D->elements+j)->nodes[3] );
					     
        fprintf(fp, "%d  %d  %d  %d\n",3,(wmt3D->elements+j)->nodes[2],
                                             (wmt3D->elements+j)->nodes[3], 
				             (wmt3D->elements+j)->nodes[0] );

        fprintf(fp, "%d  %d  %d  %d\n",3,(wmt3D->elements+j)->nodes[3],
                                             (wmt3D->elements+j)->nodes[0], 
				             (wmt3D->elements+j)->nodes[1] );
     }
  }
  /* --------output for check part finished-------- */
}
#endif

/*!
* \return   None
* \ingroup  WlzGeometry
* \brief    output the original surface correspoinding to a cut plane in VTK format.
* \param    fp               pointer pointing to a specific file.
* \param    planepointsO     the  mesh array.
* \param    linkList         pointer  pointing to a array which stores the linked list.
* \param    nIntersect       number of intersect with the tetrahedron.
* \param    noRedundancyNumP number of points in the mesh nodes.
* \author:       J. Rao, R. Baldock and B. Hill
*/
#ifdef DEBUGME
void static WlzEffWriteOriginalPlaneVTK(FILE *fp, 
			   const WlzDVertex3 *planepointsO, 
                           int **linkList, 
			   const int nIntersect, 
			   const int noRedundancyNumP)
{
  int i, j, it;

  if( planepointsO == NULL )
  {
    printf("Original Plane is zero \n");
    exit(1);
  }
  /* get the number of points first */
  it = 0;
  for(i=0; i<nIntersect; i++)
  {
     for(j=0; j< *(*(linkList+i)+0); j++)
     {
        it++;
     }
  }
 
  fprintf(fp, "# vtk DataFile Version 1.0\n");
  fprintf(fp, "WlzGeoModel test output\n");
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET POLYDATA\n");
  fprintf(fp, "%s %d %s\n","POINTS", noRedundancyNumP,"  float");
  it = 0;
  for(i=0; i<noRedundancyNumP; i++)
  { 
        fprintf(fp,  "%f  %f  %f\n",(planepointsO+it)->vtX, (planepointsO+it)->vtY, (planepointsO+it)->vtZ );
         it++;
  }
  /* get intersetpointsNumber */
  it = 0;
  for(i=0; i<nIntersect; i++)
  {
       if( *(*(linkList+i) +0) == 3 )
         it +=1;
       if( *(*(linkList+i) +0) == 4 )
         it +=2; 
  }
  fprintf(fp, "%s %d %d\n","POLYGONS", it, it * sizeof(float) ); 
  for(i=0; i< nIntersect; i++)
  {
     if( *(*(linkList+i ) +0) == 3 )
        fprintf(fp, "%d  %d  %d  %d\n", 3, *(*(linkList+i ) + 1),
                                          *(*(linkList+ i ) + 2), 
				          *(*(linkList+i ) + 3) ); 
                      
        if( *(*(linkList+i ) +0) == 4 )
	{
         fprintf(fp, "%d  %d  %d   %d\n", 3, *(*(linkList+i ) + 1),
                                               *(*(linkList+i ) + 2), 
				               *(*(linkList+i ) + 3)
				                );
         fprintf(fp, "%d  %d  %d  %d\n", 3, *(*(linkList+i ) + 1),
				               *(*(linkList+i ) + 3),
				               *(*(linkList+i ) + 4) );
        }
  }
}

#endif


/*!
* \return   twice the signed area of a triangle.
* \ingroup  WlzGeometry
* \brief    Get twice of the signed area of the triangle determined by a, b, c, positive
*              if a, b, c are oriented ccw, and negative if ccw.
* \param    a 	an array of 3D vertices.
* \param    i   an index of the array indicate the first point of the triangle.
* \param    j   an index of the array indicate the second point of the triangle.
* \param    k   an index of the array indicate the third point of the triangle.
* \author:       J. Rao, R. Baldock and B. Hill
*/
double WlzGeomTriangle2SnArea3(WlzDVertex3 *a, 
                                     int i, 
				     int j, 
				     int k)
{
   return (a+i)->vtX * (a+j)->vtY - (a+i)->vtY * (a+j)->vtX +
          (a+i)->vtY * (a+k)->vtX - (a+i)->vtX * (a+k)->vtY +
          (a+j)->vtX * (a+k)->vtY - (a+k)->vtX * (a+j)->vtY; 
}

/*!
* \return   a integer  <0 or the index of the existing point
* \ingroup  WlzAccess
* \brief    output a integer  <0  or the index of the existing point 
* \param    tempw 
* \param    planepoints
* \param    icount
* \author:       J. Rao, R. Baldock and B. Hill
*/
int static storedPointIndex(WlzDVertex3 tempw, 
                            WlzDVertex3 *planepoints, 
			    int icount)
{
  int i;
  /* check if there is a element in planepoints the same as tempw */
  for(i=0; i< icount; i++)
  {
    /*
    if(  ( fabs( tempw.vtX - (planepoints + i)->vtX )  < 0.0001 )  && 
         ( fabs( tempw.vtY - (planepoints + i)->vtY )  < 0.0001 )      ) */
    if(  (  tempw.vtX == (planepoints + i)->vtX   )  && 
         (  tempw.vtY == (planepoints + i)->vtY   )      )
      {
         return i;
      }
  }
  return -1;
}


/*!
* \return   None
* \ingroup  WlzMesh
* \brief    Make a 2D5 Triangle Mesh from Cutted Plane.
* \param    wmt2D5           2D5 mesh transform.
* \param    planepoints      the  nodes of the mesh.
* \param    linkList         pointer  pointing to a array which stores the linked list.
* \param    nIntersect       number of intersect with the tetrahedron.
* \param    noRedundancyNumP number of points in the mesh nodes.
* \author:       J. Rao, R. Baldock and B. Hill
*/
void static Make2D5TriangleMeshFromCuttedPlane(WlzMeshTransform2D5 *wmt2D5, 
                                       WlzDVertex3 *planepoints, 
                                       int **linkList, 
				       int nIntersect,  
				       int noRedundancyNumP)
{
  int i, j, k, it;
  int n[3]; 
  int np[3]; 
  /* int neighbour[3]; */
   /* copy the nodes 
      correct        */
  for(i=0; i< wmt2D5->nNodes; i++)
  {
    (wmt2D5->nodes+i)->position.vtX = (planepoints+i)->vtX;
    (wmt2D5->nodes+i)->position.vtY = (planepoints+i)->vtY;
  }

  /* copy the linkList */
  it = 0;
  for(i=0; i< nIntersect; i++)
  { 
    if(  *(*(linkList+i)+0) == 3)
    {
      for(j=0; j<3; j++)
        (wmt2D5->elements+it)->nodes[j] = *(*(linkList+i)+(j+1));
      it++;	
    }
    if( *(*(linkList+i)+0) == 4)
    { /* divided into two triangles CCW order */

      /* first triangle */
      for(j=0; j<3; j++)
      {
        (wmt2D5->elements+it)->nodes[j] = *(*(linkList+i)+(j+1));
      }	
      it++;	
      /* second triangle */
        (wmt2D5->elements+it)->nodes[0] = *(*(linkList+i)+(1));
        (wmt2D5->elements+it)->nodes[1] = *(*(linkList+i)+(3));
        (wmt2D5->elements+it)->nodes[2] = *(*(linkList+i)+(4));
      it++;	
    }
     
  }

  /* get neighbours */
  /* x0, x1, x2
    neighbour 0 -> edge x1-x2's neighbour
    neighbour 1 -> edge x0-x2's neighbour
    neighbour 2 -> edge x0-x1's neighbour
  */
  /*  initialize the neighbour to negative values  */
  for(i=0; i<wmt2D5->nElem;i++)
  {  for(k=0; k<3; k++)
     (wmt2D5->elements+i)->neighbours[k]= -10; 
  }

  /* pick up one triangle */
  for(i=0; i<wmt2D5->nElem;i++)
  { /* get nodes of the triangle */
    for(k=0; k<3; k++)
       n[k] = (wmt2D5->elements+i)->nodes[k];
    /*---- pick another one ----*/   
    for(j=0; j<i; j++)
    {  /* get nodes of this triangle */ 
       for(k=0; k<3; k++)
         np[k] = (wmt2D5->elements+j)->nodes[k];
       /* check whether triangle i and j are neighbour */
       if(IsNeighbour(n, np, &it))
       {
        (wmt2D5->elements+i)->neighbours[it] = j; 
       };
    }
    
    for(j=i+1; j<wmt2D5->nElem; j++)
    {  /* get nodes of this triangle */ 
       for(k=0; k<3; k++)
         np[k] = (wmt2D5->elements+j)->nodes[k];
       /* check whether triangle i and j are neighbour */
       if(IsNeighbour(n, np, &it))
       {
        (wmt2D5->elements+i)->neighbours[it] = j; 
       };
    }
  }

}


/*!
* \return   1 for is a neighbour or 0 for not
* \ingroup  WlzAccess
* \brief    judge whether two triangles are neighbour. 
* \param    n            
* \param    np
* \param    it   contains which edge neighbour it belongs to
* \author:       J. Rao, R. Baldock and B. Hill
*/
int static IsNeighbour(int *n, int *np, int *it)
{
  int i, j,k=0;
  int itt0= 0, itt1 = 0;
  /* ---- test for the node 0 ---- */
    for(j=0; j<3; j++)
    {
     if( n[0] == np[j] )
      {
       itt0++;
       break;
      }
    }  
    
  /* ----test the node 1 ---- */
  if(itt0 == 1)
  { /* one common node already */
    if(j==0)
    {
      for(k=1;k<3; k++)
      {
         if( n[1] == np[k] )
         { /* Now, its a neighbour */
           *it = 2;  /* 2's neighbour */
            return 1; 
         }
      }
    }
    else if(j==1)
    {
      for(k=0; k<3; k +=2)
      {
         if( n[1] == np[k] )
         { /* Now, its a neighbour */
           *it = 2;  /* 2's neighbour */
            return 1; 
         }
      }
    }
    else if(j==2)
    {
      for(k=0; k<2; k++)
      {
         if( n[1] == np[k] )
         { /* Now, its a neighbour */
           *it = 2;  /* 2's neighbour */
           return 1; 
         }
      }
    }
  }
  else if(itt0 == 0)
  {
    /* no common node yet */
    for(k=0; k<3; k++)
    {
      if( n[1] == np[k] )
      {
         itt1++;
	 break;
      }

    }

  }

  /*---- test the node 2 -----*/
  
  if(itt0 == 0 && itt1 == 0)
  {/* not a neighbour */
    return 0;
  }
  else if(itt0 == 1)
  {/* we have one node common in node 0
      test the nodes 2 */
    if(j==0)
    {
      for(i=1;i<3; i++)
      {
         if( n[2] == np[i] )
         { /* Now, its a neighbour */
           *it = 1;  /* 1's neighbour */
            return 1; 
         }
      }
    }
    else if(j==1)
    {
      for(i=0; i<3; i +=2)
      {
         if( n[2] == np[i] )
         { /* Now, its a neighbour */
           *it = 1;  /* 1's neighbour */
            return 1; 
         }
      }
    }
    else if(j==2)
    {
      for(i=0; i<2; i++)
      {
         if( n[2] == np[i] )
         { /* Now, its a neighbour */
           *it = 1;  /* 1's neighbour */
           return 1; 
         }
      }
    }
  }
  else if(itt1 == 1)
  {/* we have one node common in node 0
      test the nodes 2 */
    if(k==0)
    {
      for(i=1;i<3; i++)
      {
         if( n[2] == np[i] )
         { /* Now, its a neighbour */
           *it = 0;  /* 0's neighbour */
            return 1; 
         }
      }
    }
    else if(k==1)
    {
      for(i=0; i<3; i +=2)
      {
         if( n[2] == np[i] )
         { /* Now, its a neighbour */
           *it = 0;  /* 0's neighbour */
            return 1; 
         }
      }
    }
    else if(k==2)
    {
      for(i=0; i<2; i++)
      {
         if( n[2] == np[i] )
         { /* Now, its a neighbour */
           *it = 0;  /* 0's neighbour */
           return 1; 
         }
      }
    }
  }
  return 0;
}

/* function:     write_Wlz2D5Mesh     */
/*! 
* \ingroup      WlzGeometry
* \brief        Write 2D5 mesh for inspections.
* \return       woolz error number
* \param    fp	FILE pointer opened for writing.
* \param    cstr  the out put file name. 
* \param    wmt2D5 the 2D5 mesh transform 
* \par      Source:
*                .c
*/
WlzErrorNum  write_Wlz2D5Mesh(FILE *fp, char *cstr, WlzMeshTransform2D5 *wmt2D5)
{
  int i;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  if((fp = fopen(cstr, "w")) == NULL )
  {
    printf("cannot open the output wmt2D5Mesh.dat file.\n");
    exit(1);
  }
  /* output nodes  */

  fprintf(fp, "Nodes:\n");
  fprintf(fp, "%s %d %s\n","Number of Nodes", wmt2D5->nNodes,"  Integer");
  for(i=0; i< wmt2D5->nNodes; i++)
  { 
        fprintf(fp,  "%d %f  %f\n", i, (wmt2D5->nodes+i)->position.vtX,
	                                   (wmt2D5->nodes+i)->position.vtY
	        );
  }
 
  /* output elements */
  fprintf(fp, "Nodes:\n");
  fprintf(fp, "%s %d %s\n","Number of Elements", wmt2D5->nElem,"  Integer");
  for(i=0; i< wmt2D5->nElem; i++)
  { 
        fprintf(fp,  "%d %d  %d  %d\n", i, (wmt2D5->elements+i)->nodes[0],
	                                   (wmt2D5->elements+i)->nodes[1],
	                                   (wmt2D5->elements+i)->nodes[2]
	        );
  }
  
  /* output neighbours */
  fprintf(fp, "Neighbours:\n");
  fprintf(fp, "%s %d %s\n","Number of Elments", wmt2D5->nElem,"  Integer");
  for(i=0; i< wmt2D5->nElem; i++)
  { 
        fprintf(fp,  "%d %d  %d  %d\n", i, (wmt2D5->elements+i)->neighbours[0],
	                                   (wmt2D5->elements+i)->neighbours[1],
	                                   (wmt2D5->elements+i)->neighbours[2]
	        );
  }
  fclose(fp);
  return errNum;
}

/*!
* - Function:	WlzMeshScanWSpInit2D5
* - Returns:	WlzMeshScanWSp2D5  *:	New mesh scan workspace.
* - Purpose:	Allocate and initialise a mesh scan workspace.
* - Parameters:	WlzMeshTransform2D5 *mesh:	Mesh transform.
*		WlzErrorNum *dstErr:	Destination error pointer.
* - Author:       J. Rao, R. Baldock and B. Hill
*/
static WlzMeshScanWSp2D5 *WlzMeshScanWSpInit2D5( WlzMeshTransform2D5 *mesh,
				    	         WlzErrorNum         *dstErr )
{
  int		iIdx,
  		eIdx,
		ndIdx;
  double	ndLn,
  		eLnMin,
		eLnMax;
  WlzMeshElem	*elm;
  WlzMeshScanWSp2D5 *meshSnWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
 
 

  if((meshSnWSp = (WlzMeshScanWSp2D5 *)
  		  AlcCalloc(1, sizeof(WlzMeshScanWSp2D5))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    meshSnWSp->mesh = mesh;
    /* Compute the total number of intervals in the cutting plane mesh. */
    eIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (eIdx < mesh->nElem))
    {
      elm = mesh->elements + eIdx;
      ndIdx = elm->nodes[0];
      eLnMin = eLnMax = ndLn = (mesh->nodes + ndIdx)->position.vtY;
      
      ndIdx = elm->nodes[1];
      ndLn                   = (mesh->nodes + ndIdx)->position.vtY;

      if(ndLn < eLnMin)
      {
	eLnMin = ndLn;
      }
      else if(ndLn > eLnMax)
      {
	eLnMax = ndLn;
      }
      
      ndIdx = elm->nodes[2];
      ndLn                   = (mesh->nodes + ndIdx)->position.vtY;

      if(ndLn < eLnMin)
      {
	eLnMin = ndLn;
      }
      else if(ndLn > eLnMax)
      {
	eLnMax = ndLn;
      }
      /* counter the interval here */
      meshSnWSp->nItvs += WLZ_NINT(eLnMax) - WLZ_NINT(eLnMin) + 1;
      ++eIdx;
    }

    if((meshSnWSp->itvs = (WlzMeshScanItv *)AlcMalloc(sizeof(WlzMeshScanItv) *
						    meshSnWSp->nItvs)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else if((meshSnWSp->dElm = (WlzMeshScanDElm *)
			       AlcCalloc(mesh->nElem,
					 sizeof(WlzMeshScanDElm))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Fill in the mesh scan intervals */
    eIdx = 0;
    iIdx = 0;
    while(eIdx < mesh->nElem)
    { /* pass the zero surface element */
      /*
      l        =  (meshSnWSp->mesh->elements+eIdx)->nodes[0];
      m        =  (meshSnWSp->mesh->elements+eIdx)->nodes[1];
      n        =  (meshSnWSp->mesh->elements+eIdx)->nodes[2];
     v2[0].vtX =  (meshSnWSp->mesh->nodes+l)->position.vtX;
     v2[1].vtX =  (meshSnWSp->mesh->nodes+m)->position.vtX;
     v2[2].vtX =  (meshSnWSp->mesh->nodes+n)->position.vtX;
     v2[0].vtY =  (meshSnWSp->mesh->nodes+l)->position.vtY;
     v2[1].vtY =  (meshSnWSp->mesh->nodes+m)->position.vtY;
     v2[2].vtY =  (meshSnWSp->mesh->nodes+n)->position.vtY;
     */
     /*
      if( fabs( AreaTriangle2D(v2,0,1,2) ) > 0.001  ) */
      {
        iIdx += WlzMeshScanTriElm(meshSnWSp, eIdx, iIdx);
      }
      ++eIdx;
    }
    if(iIdx == 0)
    {/* no interval */
     

    }
  }
  /* The real number of interval */

  /*
  meshSnWSp->nItvs = iIdx; 
  */
  
  if(errNum == WLZ_ERR_NONE)
  {
    /* Sort the mesh scan intervals by line and then left column */
    qsort(meshSnWSp->itvs, meshSnWSp->nItvs, sizeof(WlzMeshScanItv),
          WlzMeshItvCmp);
  }
  else
  {
    WlzMeshScanWSpFree(meshSnWSp);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(meshSnWSp);
}

/*!
* \return   None
* \ingroup  WlzMesh
* \brief    Free's a mesh scan workspace.
* \param    mSnWSp:	Mesh scan workspace.
* \author:       J. Rao, R. Baldock and B. Hill
*/
static void	WlzMeshScanWSpFree(WlzMeshScanWSp2D5 *mSnWSp)
{
  if(mSnWSp)
  {
    if(mSnWSp->itvs)
    {
      AlcFree(mSnWSp->itvs);
    }
    if(mSnWSp->dElm)
    {
      AlcFree(mSnWSp->dElm);
    }
    AlcFree(mSnWSp);
  }
}




/*!
* \return   Sorting value for qsort
* \ingroup  WlzAccess
* \brief    Callback function for qsort(3) to sort mesh element
*           intervals by line and then left left column.

* \param    cmp0	Used to pass first mesh.
* \param    cmp1	Used to pass second mesh
* \author:       J. Rao, R. Baldock and B. Hill
*/
static int	WlzMeshItvCmp(const void *cmp0, const void *cmp1)
{
  int		rtn;
  WlzMeshScanItv *itv0,
  		 *itv1;

  itv0 = (WlzMeshScanItv *)cmp0;
  itv1 = (WlzMeshScanItv *)cmp1;
  if((rtn = (itv0->line - itv1->line)) == 0)
  {
    rtn = itv0->lftI - itv1->lftI;
  }
  return(rtn);
}

/*!
* \return   Number of intervals added from the given mesh element.
* \ingroup  WlzMesh
* \brief    Scans a single triangular mesh element into mesh intervals.
*            it to handle the special situation of triangle
*             becames a line when two node comes together        10.12.2001
* \param    mSnWSp	        Mesh scan workspace.
* \param    eIdx		Element index.
* \param    iIdx		Mesh element interval index.
* \param    
* \author:       J. Rao, R. Baldock and B. Hill
*/
static int	WlzMeshScanTriElm(WlzMeshScanWSp2D5 *mSnWSp, 
                                                 int eIdx, 
						 int iIdx)
{
  int		  count,
		  kolI,
  		  lineI,
		  ndIdx0,
  		  ndIdx1,
		  ndIdx2,
		  iCnt = 0;
  double	  tD0,
  		  tD1,
		  kolD,
		  inc;
  WlzIVertex2	  dNd[3],
  		  sNd[3];
  WlzMeshNode2D5 *nod;
  WlzMeshElem    *elm;
  WlzDVertex2	  dVx0;
  /* WlzDVertex3  dVx1; unused */
  WlzMeshScanItv *itv;

  elm = mSnWSp->mesh->elements + eIdx;
  /* Compute the integer displaced nodes of the element. 
     differ from Bill's */
  for(ndIdx0 = 0; ndIdx0 < 3; ++ndIdx0)
  {
    ndIdx1 = elm->nodes[ndIdx0];
    nod = mSnWSp->mesh->nodes + ndIdx1;
    dVx0 = nod->position;
    /* dVx1 = nod->displacement; */
    tD0 = dVx0.vtX; 
    tD1 = dVx0.vtY; 
    dNd[ndIdx0].vtX = WLZ_NINT(tD0);
    dNd[ndIdx0].vtY = WLZ_NINT(tD1);
  }
  /* Sort nodes by line coordinate, min == 0, mid == 1, max == 2. */
  if(dNd[0].vtY < dNd[1].vtY)
  {
    ndIdx0 = (dNd[0].vtY < dNd[2].vtY)? 0: 2;
  }
  else
  {
    ndIdx0 = (dNd[1].vtY < dNd[2].vtY)? 1: 2;
  }
  ndIdx1 = (ndIdx0 + 1) % 3;
  ndIdx2 = (ndIdx0 + 2) % 3;
  if(dNd[ndIdx2].vtY < dNd[ndIdx1].vtY)
  {
    ndIdx1 = ndIdx2;
    ndIdx2 = (ndIdx0 + 1) % 3;
  }
  sNd[0] = dNd[ndIdx0];
  sNd[1] = dNd[ndIdx1];
  sNd[2] = dNd[ndIdx2];

  /* Compute deltas. */
  dNd[0].vtX = sNd[0].vtX - sNd[1].vtX;
  dNd[0].vtY = sNd[0].vtY - sNd[1].vtY;
  dNd[1].vtX = sNd[1].vtX - sNd[2].vtX;
  dNd[1].vtY = sNd[1].vtY - sNd[2].vtY;
  dNd[2].vtX = sNd[2].vtX - sNd[0].vtX;
  dNd[2].vtY = sNd[2].vtY - sNd[0].vtY;

  /* If the element's nodes are not coincident scan convert it. */
  /* if(dNd[2].vtY && (dNd[0].vtX || dNd[1].vtX)) */
  {
    /* Nodes: min -> max */
    itv = mSnWSp->itvs + iIdx;
    kolD = sNd[0].vtX;
    lineI = sNd[0].vtY;
    inc = (double )(dNd[2].vtX) / (double )(dNd[2].vtY);
    count = iCnt = sNd[2].vtY - sNd[0].vtY + 1;
    while(count-- > 0)
    {
      itv->elmIdx = eIdx;
      itv->line = lineI++;
      itv->lftI = itv->rgtI = WLZ_NINT(kolD);
      kolD += inc;
      ++itv;
    }
    if(dNd[0].vtY)
    {
      /* Nodes: mid -> min */
      itv = mSnWSp->itvs + iIdx;
      kolD = sNd[0].vtX;
      inc = (double )(dNd[0].vtX) / (double )(dNd[0].vtY);
      count = sNd[1].vtY - sNd[0].vtY + 1;
      while(count-- > 0)
      {
	kolI = WLZ_NINT(kolD);
        if(kolI > itv->lftI)
	{
	  itv->rgtI = kolI;
	}
	else
	{
	  itv->lftI = kolI;
	}
        kolD += inc;
	++itv;
      }
    }
    if(dNd[1].vtY)
    {
      /* Nodes: max -> mid */
      itv = mSnWSp->itvs + iIdx + iCnt - 1;
      kolD = sNd[2].vtX;
      inc = (double )(dNd[1].vtX) / (double )(dNd[1].vtY);
      count = sNd[2].vtY - sNd[1].vtY + 1;
      while(count-- > 0)
      {
	kolI = WLZ_NINT(kolD);
        if(kolI > itv->lftI)
	{
	  itv->rgtI = kolI;
	}
	else
	{
	  itv->lftI = kolI;
	}
        kolD -= inc;
	--itv;
      }
    }
  }
  return(iCnt);
}

/*! 
* \ingroup      WlzGeometry
* \brief        Write the section view parameters in the bibtex style
*               record using the bibFile library.
*
* \return       woolz error number
* \param    fp:	FILE pointer opened for writing
* \param    cstr	record name
* \param    mSnWSp 	woolz ... 
* \par		Source:
*                Wlz3DWarpMQ_S.c
*/
WlzErrorNum Write_WlzCutScanLines(
  FILE *fp,
  char *cstr,
  WlzMeshScanWSp2D5 *mSnWSp)
{

  int i;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if((fp = fopen(cstr, "w")) == NULL )
  {
    printf("cannot open the output CutScanLines.dat file.\n");
    exit(1);
  }
  /* out put the sequence of scan lines */
  fprintf(fp, "# vtk DataFile Version 1.0\n");
  fprintf(fp, "%s %d %s\n","POINTS", mSnWSp->nItvs,"  float");
   for(i=0; i< mSnWSp->nItvs; i++)
  {
     fprintf(fp,  "%d  %d  %d  %d\n", (mSnWSp->itvs + i)->elmIdx, 
                                      (mSnWSp->itvs + i)->line,
				      (mSnWSp->itvs + i)->lftI,
				      (mSnWSp->itvs + i)->rgtI    );
  }
  fclose(fp);
  return errNum;
}

/*!
* \return   Intersect Tetrahedron Number.
* \ingroup  WlzMesh
* \brief    Get Intersect Tetrahedron Number.
* \param    zConst           z= const plane.
* \param    wmt3D            mesh transform.
* \author:       J. Rao, R. Baldock and B. Hill
*/
#ifdef DEBUGME
int  static GetIntersectTetrahedronNumber(double zConst, 
                                          WlzMeshTransform3D *wmt3D)
{
  int i;
  int j0, j1, j2, j3; 
  int icount = 0;         /* used to count how many points cutted (no redundancy counting) */
  double z0, z1, z2, z3;
  int nUp, n0up, n1up, n2up, n3up;

  for (i = 0; i < wmt3D->nElem; i++)
  {   /*  get the z-position of the four points relative to the zConst 
          Remember that we are cutting the transformed mesh */
      /* Also, We only cut the mesh belong to the object! */
      /* exclude the tetrahedrons not belong to the object */

     j0    =  (wmt3D->elements+i )->nodes[0];
     j1    =  (wmt3D->elements+i )->nodes[1];
     j2    =  (wmt3D->elements+i )->nodes[2];
     j3    =  (wmt3D->elements+i )->nodes[3];
    
     if( (j0 != 0) || (j1 !=0 ) || (j2 !=0) || (j3 !=0 ) )
     {
         z0    =  (wmt3D->nodes+j0)->position.vtZ + (wmt3D->nodes+j0)->displacement.vtZ   - zConst;
         z1    =  (wmt3D->nodes+j1)->position.vtZ + (wmt3D->nodes+j1)->displacement.vtZ   - zConst;
         z2    =  (wmt3D->nodes+j2)->position.vtZ + (wmt3D->nodes+j2)->displacement.vtZ   - zConst;
         z3    =  (wmt3D->nodes+j3)->position.vtZ + (wmt3D->nodes+j3)->displacement.vtZ   - zConst;
         nUp   =  0;
         n0up  =  n1up = n2up = n3up = 0;
         /* printf("%lf  %lf  %lf  %lf\n", z0, z1, z2, z3); */

         /* judge whether there is a intersection with zConst-plane 
         only 3 cases if intersected

	 1) one   point up,  three points down 
	 2) two   point up,  two   points down 
	 3) three point up,  one   points down 
	
         Not intersected only two cases:
	 4) four points all up
	 5) four points all down
         */

         if( z0 > 0. )
         {
            n0up = 1;
            nUp++;
         }
     
         if( z1 > 0.  )
         {
            n1up = 1;
             nUp++;
         }
     
         if( z2 > 0.  )
         {
            n2up = 1;
             nUp++;
         }
     
         if( z3 > 0.  )
         {
            n3up = 1;
             nUp++;
         }
      
         /* Now the cutted cases  */
         if( (nUp != 0) && (nUp != 4) )
          icount++;
      }   

   }
   return icount;
}
#endif


/*!
* \return   None
* \ingroup  WlzMesh
* \brief    Get the source surface of the cut plane.
* \param    nIntersect:       number of intersect with the tetrahedron.
* \param    intersectIndex    the index which indicating which tetrahetron has been interseted by a give plane 
* \param    wmt3D:            3D mesh transform.
* \param    planepoints:      the  nodes of the cut mesh.
* \param    planepointsO:     the  nodes of the surface mesh corresponding to the cut mesh.
* \param    linkList:         pointer  pointing to a array which stores the linked list.
* \param    
* \param    
* \author:       J. Rao, R. Baldock and B. Hill
*/
void  static GetTheSourceSurfaceOfTheCutPlanePoints(int nIntersect, 
                                             int *intersectIndex,
                                             const WlzMeshTransform3D  *wmt3D, 
					     WlzDVertex3 *planepoints, 
					     WlzDVertex3 *planepointsO, 
					     int **linkList
					     )
{ 

  int it;
  int i, j, k;
  WlzDVertex3 sr0, sr1, sr2, sr3;
  WlzDVertex3 targ0, targ1, targ2, targ3;
  double **Affine3D4pointsTrFun;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  /* allocate memory for Affine transformation */
  if( (AlcDouble2Malloc(&Affine3D4pointsTrFun, 3, 4) !=  ALC_ER_NONE)  )
  {
      /* This error should be handled. */
      /* errNum = WLZ_ERR_MEM_ALLOC; unused */
  }

  
  it   = 0;
  for(i=0; i<nIntersect; i++)
  { 
    j = *(intersectIndex+i);                     /* here j is the intersected tetrahedron */
    /* the corresponding source nodes of the mesh */
    k = (wmt3D->elements+j )->nodes[0];
    sr0.vtX =  ( wmt3D->nodes + k )->position.vtX;
    sr0.vtY =  ( wmt3D->nodes + k )->position.vtY;
    sr0.vtZ =  ( wmt3D->nodes + k )->position.vtZ;
   
    targ0.vtX =  sr0.vtX + ( wmt3D->nodes + k )->displacement.vtX;
    targ0.vtY =  sr0.vtY + ( wmt3D->nodes + k )->displacement.vtY;
    targ0.vtZ =  sr0.vtZ + ( wmt3D->nodes + k )->displacement.vtZ;

    k =  (wmt3D->elements+j )->nodes[1]; 
    sr1.vtX =  ( wmt3D->nodes + k )->position.vtX;
    sr1.vtY =  ( wmt3D->nodes + k )->position.vtY;
    sr1.vtZ =  ( wmt3D->nodes + k )->position.vtZ;
   
    targ1.vtX =  sr1.vtX + ( wmt3D->nodes + k )->displacement.vtX;
    targ1.vtY =  sr1.vtY + ( wmt3D->nodes + k )->displacement.vtY;
    targ1.vtZ =  sr1.vtZ + ( wmt3D->nodes + k )->displacement.vtZ;

    k =  (wmt3D->elements+j )->nodes[2]; 
    sr2.vtX =  ( wmt3D->nodes + k )->position.vtX;
    sr2.vtY =  ( wmt3D->nodes + k )->position.vtY;
    sr2.vtZ =  ( wmt3D->nodes + k )->position.vtZ;
   
    targ2.vtX =  sr2.vtX + ( wmt3D->nodes + k )->displacement.vtX;
    targ2.vtY =  sr2.vtY + ( wmt3D->nodes + k )->displacement.vtY;
    targ2.vtZ =  sr2.vtZ + ( wmt3D->nodes + k )->displacement.vtZ;

    k = (wmt3D->elements+j )->nodes[3];
    sr3.vtX =  ( wmt3D->nodes + k )->position.vtX;
    sr3.vtY =  ( wmt3D->nodes + k )->position.vtY;
    sr3.vtZ =  ( wmt3D->nodes + k )->position.vtZ;
   
    targ3.vtX =  sr3.vtX + ( wmt3D->nodes + k )->displacement.vtX;
    targ3.vtY =  sr3.vtY + ( wmt3D->nodes + k )->displacement.vtY;
    targ3.vtZ =  sr3.vtZ + ( wmt3D->nodes + k )->displacement.vtZ;
 
  /*  Get 3D Affine trasformation from a tetrahedron nodes corresponding points  */
  /*
  test codes:
  input:
              source tetrahedron: 4 points (10,10,0), (10,10,10), (0,10,10), (0,0,10) 
              target tetrahedron: 4 points (0,10,0), (10,10,10), (0,10,10), (10,0,0)
  output:
               Affine3D4pointsTrFun[3][4] contains the affine transformation parameters:

  A =:
	       a11  a12  a13  a14
	       a21  a22  a23  a24
	       a31  a32  a33  a34
                0    0    0    1

  sr0.vtX = 10.;
  sr0.vtY = 10.;
  sr0.vtZ =  0.;

  sr1.vtX = 10.;
  sr1.vtY = 10.;
  sr1.vtZ = 10.;

  sr2.vtX = 0.;
  sr2.vtY = 10.;
  sr2.vtZ = 10.;

  sr3.vtX = 0.;
  sr3.vtY = 0.;
  sr3.vtZ = 10.;

  targ0.vtX =  0.;
  targ0.vtY = 10.;
  targ0.vtZ =  0.;

  targ1.vtX = 10.;
  targ1.vtY = 10.;
  targ1.vtZ = 10.;

  targ2.vtX = 0.;
  targ2.vtY = 10.;
  targ2.vtZ = 10.;

  targ3.vtX = 10.;
  targ3.vtY = 0.;
  targ3.vtZ = 0.;

  for(j=0; j< 3; j++) 
  {
                 *( *(Affine3D4pointsTrFun + j) + 0) = 0.0;
		 *( *(Affine3D4pointsTrFun + j) + 1) = 0.0;
                 *( *(Affine3D4pointsTrFun + j) + 2) = 0.0;
		 *( *(Affine3D4pointsTrFun + j) + 3) = 0.0;
  }
  */

    /* get the 3D Affine transform for this tetrahedron */ 
    /* target to source */
    /* fprintf(fp,"Alffine i = %d\n",i); */
    /* here is the target to source */
    WlzMakeAffine3D4pointsTrFn(targ0, targ1, targ2, targ3, sr0, sr1, sr2, sr3,  Affine3D4pointsTrFun); 
    /* use this transform to transfer the cutting plane back to the source space */
    /* output for check */
    /*
    for(j=0; j< 3; j++) 
    {
       fprintf(fp, "TRANS:%d %f  %f  %f  %f\n",   j, 
                 *( *(Affine3D4pointsTrFun + j) + 0),  
		 *( *(Affine3D4pointsTrFun + j) + 1),
                 *( *(Affine3D4pointsTrFun + j) + 2),  
		 *( *(Affine3D4pointsTrFun + j) + 3)     );
    }
    fprintf(fp, "\n");   */
    /* transfer the plane points back to the source space */
    for(j = 0; j< *(*(linkList+i)+0); j++)
    {
     /* 
      fprintf(fp,  "%f  %f  %f\n", (planepoints + it )->vtX,
                                   (planepoints + it )->vtY, 
                                   (planepoints + it )->vtZ  );
				   
      fprintf(fp,  "%f  %f  %f\n",
      */
      it = *(*(linkList+i)+j+1);
      (planepointsO+it)->vtX  =         (planepoints + it )->vtX                     +
                (planepoints + it )->vtX * ( *( *(Affine3D4pointsTrFun + 0) + 0) )   +
                (planepoints + it )->vtY * ( *( *(Affine3D4pointsTrFun + 0) + 1) )   +
                (planepoints + it )->vtZ * ( *( *(Affine3D4pointsTrFun + 0) + 2) )   +
	                                     *( *(Affine3D4pointsTrFun + 0) + 3);
					     
      (planepointsO+it)->vtY  =	        (planepoints + it )->vtY                     +					     
                (planepoints + it )->vtX * ( *( *(Affine3D4pointsTrFun + 1) + 0) )   +
                (planepoints + it )->vtY * ( *( *(Affine3D4pointsTrFun + 1) + 1) )   +
                (planepoints + it )->vtZ * ( *( *(Affine3D4pointsTrFun + 1) + 2) )   +
	                                     *( *(Affine3D4pointsTrFun + 1) + 3);
					     
      (planepointsO+it)->vtZ  =	        (planepoints + it )->vtZ                     +					     
                (planepoints + it )->vtX * ( *( *(Affine3D4pointsTrFun + 2) + 0) )   +
                (planepoints + it )->vtY * ( *( *(Affine3D4pointsTrFun + 2) + 1) )   +
                (planepoints + it )->vtZ * ( *( *(Affine3D4pointsTrFun + 2) + 2) )   +
			                     *( *(Affine3D4pointsTrFun + 2) + 3);
	/* it++ */;	
    } 
  }
  /*  free the memory */
  (void )AlcDouble2Free(Affine3D4pointsTrFun);
  
}

/*!
* \return   Error code. 
* \ingroup  WlzTransform
* \brief    Get a  3D  mesh transfrom from a basis transform and a 3D mesh.
* \param    wmt3D            mesh transform.
* \param    basisTr:       basis transform.
* \author:       J. Rao, R. Baldock and B. Hill
*/
WlzErrorNum  WlzGetTransformedMesh(WlzMeshTransform3D *wmt3D, 
                                 WlzBasisFnTransform* basisTr){
  WlzDVertex3 vx4;
  WlzErrorNum werro = WLZ_ERR_NONE;
  int i;

  for(i=0; i<wmt3D->nNodes; i++)
  {
    vx4 = (wmt3D->nodes + i)->position;
    
    (wmt3D->nodes + i)->displacement  = WlzBasisFnValueMQ3D(basisTr->basisFn, vx4);  

  }
  return werro;

}


/*!
* \return   None
* \ingroup  WlzTransform
* \brief    Get the mesh transform.
* \param    wmt2D5:        the 2D5 mesh transform.
* \param    planepointsO:  the nodes of the original surface.
* \author:       J. Rao, R. Baldock and B. Hill
*/
void  static Make2D5MeshDisplacement(WlzMeshTransform2D5 *wmt2D5,
				     WlzDVertex3 *planepointsO )
{
  int i;
  for(i=0; i< wmt2D5->nNodes; i++)
  { 
      (wmt2D5->nodes+i)->displacement.vtX =  (planepointsO + i)->vtX 
                                           - (wmt2D5->nodes+ i)->position.vtX;
					   
      (wmt2D5->nodes+i)->displacement.vtY =  (planepointsO + i)->vtY
                                           - (wmt2D5->nodes+ i)->position.vtY;

      (wmt2D5->nodes+i)->displacement.vtZ =  (planepointsO + i)->vtZ 
                                           -  wmt2D5->zConst ;
  }
}


/*!
* \return   None
* \ingroup  WlzGeometry
* \brief    output the orginal surface corresponding to the cut plane represented by the postion of
*           the same mesh.
* \param    fp               pointer pointing to a specific file.
* \param    wmt2D5:          the 2D5 mesh transform.
* \author:       J. Rao, R. Baldock and B. Hill
*/
void WlzEffWriteOriginalPlaneVTKByDis(FILE *fp, WlzMeshTransform2D5 *wmt2D5)
{ int i;
 if( wmt2D5 == NULL )
  {
    printf("no Mesh transform of 2D5 \n");
    exit(1);
  }
 
  fprintf(fp, "# vtk DataFile Version 1.0\n");
  fprintf(fp, "WlzGeoModel test output\n");
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET POLYDATA\n");
  fprintf(fp, "%s %d %s\n","POINTS", wmt2D5->nNodes,"  float");
  for(i=0; i< wmt2D5->nNodes; i++)
  { 
        fprintf(fp,  "%f  %f  %f\n", (wmt2D5->nodes+i)->position.vtX + (wmt2D5->nodes+i)->displacement.vtX,
                                     (wmt2D5->nodes+i)->position.vtY + (wmt2D5->nodes+i)->displacement.vtY,
                                      wmt2D5->zConst + (wmt2D5->nodes+i)->displacement.vtZ
	       );
  }
  
  fprintf(fp, "%s %d %ld\n","POLYGONS", wmt2D5->nElem,  (long)(wmt2D5->nElem) * sizeof(float) ); 
  for(i=0; i< wmt2D5->nElem ; i++)
  {
        fprintf(fp, "%d  %d  %d  %d\n", 3, (wmt2D5->elements+i)->nodes[0],
                                               (wmt2D5->elements+i)->nodes[1],
                                               (wmt2D5->elements+i)->nodes[2]
               );
  }
}

/*!
* \return   None
* \ingroup  WlzGeometry
* \brief    output the cutted plane.
* \param    fp               pointer pointing to a specific file.
* \param    wmt2D5           mesh transform.
* \author:       J. Rao, R. Baldock and B. Hill
*/
void WlzEffWriteOriginalPlaneVTKByPos(FILE *fp, WlzMeshTransform2D5 *wmt2D5)
{ int i;
 if( wmt2D5 == NULL )
  {
    printf("no Mesh transform of 2D5 \n");
    exit(1);
  }
 
  fprintf(fp, "# vtk DataFile Version 1.0\n");
  fprintf(fp, "WlzGeoModel test output\n");
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET POLYDATA\n");
  fprintf(fp, "%s %d %s\n","POINTS", wmt2D5->nNodes,"  float");
  for(i=0; i< wmt2D5->nNodes; i++)
  { 
        fprintf(fp,  "%f  %f  %f\n", (wmt2D5->nodes+i)->position.vtX,
                                     (wmt2D5->nodes+i)->position.vtY,
                                      wmt2D5->zConst
	       );
  }
  
  fprintf(fp, "%s %d %ld\n","POLYGONS", wmt2D5->nElem,  (long )(wmt2D5->nElem) * sizeof(float) ); 
  for(i=0; i< wmt2D5->nElem ; i++)
  {
        fprintf(fp, "%d  %d  %d  %d\n", 3, (wmt2D5->elements+i)->nodes[0],
                                               (wmt2D5->elements+i)->nodes[1],
                                               (wmt2D5->elements+i)->nodes[2]
               );
  }
}


/*!
* - Function:	WlzMeshTransform3DObj
* - Ingroup:    WlzTransform
* - Returns:	WlzObject *:		Transformed object, NULL on
*					error.
* - Purpose:	Transforms a woolz object using a the given mesh
*		transform.
* - Global refs:	-
* - Parameters:	
*         -#    *srcObj:	Woolz Object to be transformed.
*	  -#	*wmt3D:         Given mesh transform to apply.
*	  -#	*wmt2D5:        Given 2D5 mesh transform fram.
*	  -#	 interp:        Level of interpolation to use.
*	  -#	*dstErr:	Destination error pointer, may be NULL.
* - Author:       J. Rao, R. Baldock and B. Hill
*/
/*
*     by J. Rao, R. Baldock and B. Hill
*        30.10.2001	
*		
*/
#ifdef DEBUGME
static WlzObject	*WlzMeshTransformObj3D(WlzObject *srcObj,
				       WlzMeshTransform3D  *wmt3D,
				       WlzMeshTransform2D5 *wmt2D5,
				       WlzInterpolationType interp,
				       WlzErrorNum *dstErr,
                                 int         *intersectIndex,
				 WlzDVertex3 *planepoints,
				 WlzDVertex3 *planepointsO,
				 int        **linkList,
				 int          nInterSectEsimate,
                                 int          numEstTrangularsElem,
                                 int          numEstCutingPoints,
                                 WlzMeshScanWSp2D5 *mSnWSp,
				 WlzIBox3        bBox0
				       )
{
  WlzObject	*dstObj = NULL;
  /*
  WlzMeshTransform *aMesh = NULL;
  */
  WlzErrorNum	errNum = WLZ_ERR_NONE;


  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(wmt3D == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((wmt3D->nElem < 0) || (wmt3D->nNodes < 0))
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else
  { /*
    aMesh = WlzMeshTransformAdapt(gMesh, minElmArea, &errNum);
    */
  }
  if(errNum == WLZ_ERR_NONE)
  { /*
    dstObj = WlzMeshTransformObjPrv(srcObj, aMesh, interp, &errNum);
    */
    dstObj = WlzMeshTransformObjPrv3D(srcObj, wmt3D, wmt2D5, interp, &errNum,
                              intersectIndex,
			      planepoints,
			      planepointsO,
			      linkList,
			      nInterSectEsimate,
                              numEstTrangularsElem,
                              numEstCutingPoints,
			      mSnWSp,
			      bBox0
                                     );
  }
  /*
  if(aMesh)
  {
    (void )WlzMeshFreeTransform3D(aMesh);
  }
  */
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);

}
#endif
/*!
* - Function:	WlzMeshTransformObjPrv3D
* - Ingroup:    WlzTransform.
* - Returns:	WlzObject *:		Transformed object, NULL on
*					error.	
* - Purpose:	Private version of WlzMeshTransformObj3D() which
*		transforms a woolz object using a the given mesh
*		transform.
* - Global refs:	-
* - Parameters:	
*       -#  *srcObj:	Woolz object to be transformed.
*	-#  *mesh:      Given mesh transform to apply.
*	-#  *mesh:      Given 2D5 mesh transform frame to apply. 
*	-#   interp:    Level of interpolation touse.
*	-#  *dstErr:	Destination error pointer, may be NULL.
* - Author:       J. Rao, R. Baldock and B. Hill
*/
/*
*	        J. Rao 31.10.2001
*					 			        *
*/
#ifdef DEBUGME
static WlzObject *WlzMeshTransformObjPrv3D( WlzObject *srcObj,
				            WlzMeshTransform3D *mesh,
				            WlzMeshTransform2D5 *wmt2D5,
					    WlzInterpolationType interp,
				            WlzErrorNum *dstErr,
                                 int         *intersectIndex,
				 WlzDVertex3 *planepoints,
				 WlzDVertex3 *planepointsO,
				 int        **linkList,
				 int          nInterSectEsimate,
                                 int          numEstTrangularsElem,
                                 int          numEstCutingPoints,
                                 WlzMeshScanWSp2D5 *mSnWSp,
				 WlzIBox3        bBox0
					    )
{
  WlzDomain	  dstDom;
  WlzPlaneDomain *dstPDom   = NULL;
  WlzValues	  dstValues;
  WlzValues	  srcValues;
  WlzObject	 *tObj0,
  		 *tObj1,
		 *dstObj = NULL;
  WlzErrorNum	  errNum = WLZ_ERR_NONE;
  WlzIBox3        bBox;
  
  /*  2D object used as working variables  */
  WlzObject      *dstObj2D;
  WlzValues       dummyValues;
  WlzDomain	  dummyDom;
  int             planeIdx,
                  cutPosition,
                  planeCount;
  WlzPixelV       bckgrnd;


 
  dstDom.core      = NULL;
  srcValues.core   = NULL;
  dummyDom.core    = NULL;
  dummyValues.core = NULL;
  dstValues.core   = NULL;


  switch(srcObj->type)
  {
    case WLZ_EMPTY_OBJ:
      dstObj = WlzMakeEmpty(&errNum);
      break;
    case WLZ_3D_DOMAINOBJ:
      if(srcObj->domain.core == NULL)
      {
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else if( srcObj->domain.core->type == WLZ_PLANEDOMAIN_DOMAIN )
      { /*
        printf("bBox.xmin%d\n",bBox0.xMin);
        printf("bBox.xmax%d\n",bBox0.xMax);
        printf("bBox.ymin%d\n",bBox0.yMin);
        printf("bBox.ymax%d\n",bBox0.yMax);
        printf("bBox.zmin%d\n",bBox0.zMin);
        printf("bBox.zmax%d\n",bBox0.zMax);
  	exit(0);
	*/
	/*
        if( ( bBox0.xMax > bBox0.xMin ) &&
	     ( 
	        ( bBox0.yMax > bBox0.yMin ) &&
                ( bBox0.zMax > bBox0.zMin )
	     )	
	  )
	{
	  bBox.xMin = bBox0.xMin;
	  bBox.xMax = bBox0.xMax;
	  bBox.yMin = bBox0.yMin;
	  bBox.yMax = bBox0.yMax;
	  bBox.zMin = bBox0.zMin;
	  bBox.zMax = bBox0.zMax;
	}
	else
	*/
	{
        /* we need bounding box transformation here which can be 
	   done by using the mesh trasform */
         bBox = WlzMQ3DTransformBBoxI3(mesh, &errNum);

         if( ( bBox0.zMax > bBox0.zMin ) /*  &&
	     ( 
	        ( bBox0.xMax > bBox0.xMin ) &&
                ( bBox0.yMax > bBox0.yMin )
	     )	
	     */
	   )
	 {
	   /* bBox.xMin = bBox0.xMin;
	      bBox.xMax = bBox0.xMax;
	      bBox.yMin = bBox0.yMin;
	      bBox.yMax = bBox0.yMax; */
	   bBox.zMin = bBox0.zMin;
	   bBox.zMax = bBox0.zMax;
	 }




	 
	}
	if(errNum == WLZ_ERR_NONE)
	{
           /* Create a plane domain and a destination object. */
           dstPDom = WlzMakePlaneDomain( WLZ_PLANEDOMAIN_DOMAIN,
	                                 bBox.zMin, bBox.zMax,
	                                 bBox.yMin, bBox.yMax,
	                                 bBox.xMin, bBox.xMax,
						       &errNum );
           bckgrnd = WlzGetBackground(srcObj, &errNum);

	   dstValues.vox = WlzMakeVoxelValueTb( WLZ_VOXELVALUETABLE_GREY,
	                                        dstPDom->plane1,
					        dstPDom->lastpl,
					        bckgrnd, NULL, &errNum );
	   /* Now make a 2D object  */
	   if(errNum == WLZ_ERR_NONE)
	   {
	     dummyDom.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
	                             bBox.yMin, 
				     bBox.yMax,
	                             bBox.xMin, 
				     bBox.xMax,
				     &errNum
					     );
	   }

	   if(errNum == WLZ_ERR_NONE)
	   {
	     dstObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, dummyDom,
	                        dummyValues, NULL, NULL, &errNum);
	   }
	
           /*----------  */
	   dstDom.p = WlzPDomFromBBox(bBox,&errNum);
	   planeIdx = 0;
	   planeCount = bBox.zMax  - bBox.zMin + 1;
	   /* ------------ */
	   if(errNum == WLZ_ERR_NONE)
	   {
	     dstObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, dstDom, dstValues, NULL, NULL, &errNum);
	   }
	   /* get values */

	   
	   /* Now circle the plane one by one  */
	     cutPosition = bBox.zMin;
	   while((errNum == WLZ_ERR_NONE) && (planeCount-- > 0))
	   {
	      /* here the grey value should be added to disObj2D   */
             	      /* debug */
	      /* if(cutPosition == 20) */
	      { 
                errNum =  WlzMeshTransformValues3D(dstObj2D,
	                                         srcObj, 
					         mesh,
						 wmt2D5,
					         cutPosition,
					         interp,
                                 intersectIndex,
				 planepoints,
				 planepointsO,
				 linkList,
 				 nInterSectEsimate,
                                 numEstTrangularsElem,
                                 numEstCutingPoints,
				 mSnWSp
						 );
						 
	      }	
	      if(errNum == WLZ_ERR_NONE)
	      {
	        /*  assinDomain  */
		 *(dstObj->domain.p->domains + planeIdx) = WlzAssignDomain(dstObj2D->domain, NULL);
                /* assigin back to the 3D dstObj one layer */
	        /* printf("planeIndex: %d\n", planeIdx);  */
		*(dstObj->values.vox->values + planeIdx ) =
		  WlzAssignValues(dstObj2D->values, NULL); 
	      }
              ++planeIdx;
              cutPosition++;
	   }
	}
	dstObj2D->domain.core = NULL;
	dstObj2D->values.core = NULL;
	WlzFreeObj(dstObj2D);
 	
        /* get domain first we need a function to do this here */
	/*
	dstDom.p = ...(srcObj, trans, &errNum); */ 
	tObj1 = NULL;
	tObj0 = WlzObjToBoundary(srcObj, 1, &errNum);
      }
      break;
    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}
#endif


/*!
* - Function:     WlzMQ3DTransformBBoxI3
* - Ingroup:      WlzTransform.
* - Returns:      Wlooz Integer Bounding box. 
* - Purpose:      Get a bounding box for the transformed mesh. 
* - Global refs:  -     
* - Parameters:  
*      -#  *wmt3D:   Given mesh transform.
*      -#  *errNum:  Woolz Error number 
* - Author:       J. Rao, R. Baldock and B. Hill
*/
WlzIBox3  WlzMQ3DTransformBBoxI3(WlzMeshTransform3D *wmt3D,  WlzErrorNum *errNum)
{
  int		i;
  WlzIBox3	dstBox;
  int           tempx, tempy, tempz;
  *errNum  =    WLZ_ERR_NONE;
  if(wmt3D == NULL)
  {
    *errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  dstBox.xMin = dstBox.xMax = (wmt3D->nodes)->position.vtX + (wmt3D->nodes)->displacement.vtX;
  dstBox.yMin = dstBox.yMax = (wmt3D->nodes)->position.vtY + (wmt3D->nodes)->displacement.vtY;
  dstBox.zMin = dstBox.zMax = (wmt3D->nodes)->position.vtZ + (wmt3D->nodes)->displacement.vtZ;

  for(i=1; i<wmt3D->nNodes; i++)
  {
    tempx = (int)( (wmt3D->nodes+i)->position.vtX + (wmt3D->nodes+i)->displacement.vtX );  
    if(dstBox.xMin > tempx )
    {
       dstBox.xMin = tempx; 
    }
    if(dstBox.xMax < tempx )
    {
       dstBox.xMax = tempx; 
    }
    
    tempy = (int) ( (wmt3D->nodes+i)->position.vtY + (wmt3D->nodes+i)->displacement.vtY );  
    if(dstBox.yMin > tempy )
    {
       dstBox.yMin = tempy; 
    }
    if(dstBox.yMax < tempy )
    {
       dstBox.yMax = tempy; 
    }

    tempz = (int) ( (wmt3D->nodes+i)->position.vtZ + (wmt3D->nodes+i)->displacement.vtZ );  
    if(dstBox.zMin > tempz )
    {
       dstBox.zMin = tempz; 
    }
    if(dstBox.zMax < tempz )
    {
       dstBox.zMax = tempz; 
    }
    
  }
  return(dstBox);

}	

/*!
* - Function:     WlzMeshTransformValues3D
* - Ingroup:      WlzTransform
* - Returns:      WlzErrorNum:            Error number. 
* - Purpose:      Creates a new value table, fills in the values and 
*                 adds it to the given new object. 
* - Global refs:  -     
* - Parameters:  
*      -#  *dstObj:      Partialy transformed 2D object with a valid domain.
*      -#  *srcObj:      3D domain object which is being transformed.  
*      -#  *mesh:        Given mesh transform.
*      -#  *wmt2D5:      Given 2D5 mesh transform.
*      -#   cutPosition: The current plane considered 
*      -#   interp:      Level of interpolation to use. 
* - Author:       J. Rao, R. Baldock and B. Hill
*/
/*
*                                                                       *
*               J. Rao  31.10.2001                                      *
*                                                                       *
*/
WlzErrorNum static WlzMeshTransformValues3D(       WlzObject *dstObj,
					    WlzObject *srcObj,
					    WlzMeshTransform3D *mesh,
					    WlzMeshTransform2D5 *wmt2D5,
					    int cutPosition,
					    WlzInterpolationType interp,
                                 int         *intersectIndex,
				 WlzDVertex3 *planepoints,
				 WlzDVertex3 *planepointsO,
				 int        **linkList,
                                 int          nInterSectEsimate,
				 int          numEstTrangularsElem,
                                 int          numEstCutingPoints,
				 WlzMeshScanWSp2D5 *mSnWSp 
					    )
{
  /* int	mItvIdx; unused */
  int		i, j, k, k0, k1, k2, ip;
  int           i_lineCom, i_columnP, debugNum, i_prevC;
  int           ix,iy,iz;
  
  		
		
		
		
		
		
		
 
  double        xs, ys, zs;
  /* WlzIVertex2	dPosI; unused */
  		
  
  WlzGreyP	dGP;
  WlzGreyType	newGreyType;
  WlzPixelV	bkdV;
  WlzValues	newValues;
#ifdef WLZ_MESH_DEBUG
  WlzMeshScanItv oldItvDebug;
#endif /* WLZ_MESH_DEBUG */
  WlzMeshScanItv *mItv;

  WlzGreyValueWSpace *gVWSp = NULL;
  WlzGreyWSpace       gWSp;
  WlzIntervalWSpace   iWSp;
  WlzErrorNum         errNum = WLZ_ERR_NONE;
  WlzAffineTransform  *trans = NULL;
  
  WlzDVertex2 sr0; 
  WlzDVertex2 sr1; 
  WlzDVertex2 sr2;
  WlzDVertex3 targ0;
  WlzDVertex3 targ1; 
  WlzDVertex3 targ2;

  FILE  *fp = NULL;

  interp = WLZ_INTERPOLATION_NEAREST;  /* Use the nearest neighbour */

  /* step by step */
  /* allocate memory */
  if(errNum == WLZ_ERR_NONE)
  {
    if(((trans = (WlzAffineTransform *)
                 AlcCalloc(1, sizeof(WlzAffineTransform))) == NULL) ||
       (AlcDouble2Calloc(&trans->mat, 4, 4) != ALC_ER_NONE))
    {
      if(trans)
      {
        AlcFree(trans);
	trans = NULL;
      }
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }

  /* initialize the mat */
  for(i=0; i<4; i++)
  {
    for(j=0; j<4; j++)
    {
      trans->mat[i][j] = 0.0;
    }  
  }
  newValues.core = NULL;
  /*  this works both for 2D and 3D object  */
  bkdV = WlzGetBackground(srcObj, &errNum);

  if(errNum == WLZ_ERR_NONE)
  {  
    /*  this works both for 2D and 3D object  */
    /* newGreyType = WlzGreyTableTypeToGreyType(srcObj->values.v->type,
    					     &errNum); */
    newGreyType = WlzGreyTableTypeToGreyType(srcObj->values.vox->values->core->type,
    					     &errNum);
   }
  
  if((errNum == WLZ_ERR_NONE) && (bkdV.type != newGreyType))
  {
    errNum = WlzValueConvertPixel(&bkdV, bkdV, newGreyType);
  }
  
  if(errNum == WLZ_ERR_NONE)
  {
    /* dstObj must be WLZ_2D_DOMAINOBJ */
    newValues.v = WlzNewValueTb(dstObj,  srcObj->values.vox->values->core->type,
    				bkdV, &errNum);
  }

  if(errNum == WLZ_ERR_NONE)
  {
    dstObj->values = WlzAssignValues(newValues, &errNum);
  }

  if(errNum == WLZ_ERR_NONE)
  {
    /* mItvIdx = 0; unused */
    /* get 2D5 mesh */
    /* the following two functions can be used but it seems there are
       some differents and the first is better. The later one seems
       has some bugs */

    /* printf("Cutposition:  %d\n", cutPosition ); */
    WlzGet2D5TrangularMeshFrom3DMesh( wmt2D5, mesh, cutPosition, &errNum,
				   intersectIndex,
				   planepoints,
				   planepointsO,
				   linkList,
                                   nInterSectEsimate,
				   numEstTrangularsElem,
                                   numEstCutingPoints
                                  );
    /*
    WlzGet2D5Transform( mesh,
                     cutPosition,
                     wmt2D5,
		     &errNum
                   );
    */		   
    /* 
    WlzIsoIntersectWithTetrahadronIndexAnDETC( cutPosition, 
                                               mesh, 
					       wmt2D5
					     );

   */
  }

  /*  Here is the place to put values  */
  if(errNum == WLZ_ERR_NONE)
  {
     if(errNum == WLZ_ERR_NONE)
     {
        gVWSp = WlzGreyValueMakeWSp(srcObj, &errNum);
     }
     /*  preparations  */
     errNum   =  WlzInitGreyScan(dstObj, &iWSp, &gWSp);
     if(errNum != WLZ_ERR_NONE)
     {
        printf("error : %d\n", errNum);
        exit(1);
     }

     /*  get the mesh scanWorkSpace  */
     if( (wmt2D5->nElem > 1) && !WlzIsSinglePointCutMesh(wmt2D5) && !WlzNoAreaGreaterThanOne(wmt2D5)  )
     {
        mSnWSp = WlzMeshScanWSpInit2D5(wmt2D5, &errNum);
     }
  }
  
  if(errNum == WLZ_ERR_NONE)
  {
   /* here we should change something  */
   if(wmt2D5->nElem > 1 )
   {
     /*printf("cutPosition:  %d\n", cutPosition );
     */
     /* exclude the one points cut and the case when all cut area is less than one: */
     /* printf("SingleCut?  %d  noAreaLarger?  %d\n", WlzIsSinglePointCutMesh(wmt2D5),  WlzNoAreaGreaterThanOne(wmt2D5) ); */
     if( !WlzIsSinglePointCutMesh(wmt2D5) && !WlzNoAreaGreaterThanOne(wmt2D5) )
     {
       /* Here is the place to put values  */
       if(errNum == WLZ_ERR_NONE)
       {
 	     /* output for check  */
	     /*
	       printf("nItvs,  %d\n\n",mSnWSp->nItvs);
	     */  
             /* start from the first one  use the sCanLine to control the scan */
              mItv      =  mSnWSp->itvs;
	      /* Initialize  */
              WlzNextGreyInterval(&iWSp);
	      /* dPosI.vtX = iWSp.lftpos; */
 	      /* dPosI.vtY = iWSp.linpos; */
	      dGP       = gWSp.u_grintptr;
	      i_lineCom = iWSp.linpos;
	      i_columnP = iWSp.lftpos;
              i_prevC   =  mItv->lftI;  /* used to trace the current position to be filled with grey value */
	      if(mItv->line > 1000 || mItv->line < -1000 )
	      {
	        printf("something wrong here");
	      }
	      /* scan the interverl scan lines */
	      for( i=0; i<mSnWSp->nItvs; i++ )
	      {

		     /* output for debug */
		     /*
		     if(mItv->line == -47 )
		     {
		       printf("leftI: %d  rightI:  %d\n",mItv->lftI,mItv->rgtI );
		     }
		     */
                 
	         /* go to the correct line position */
	         while( (errNum == WLZ_ERR_NONE) && (i_lineCom < mItv->line) )
		 {
                   WlzNextGreyInterval(&iWSp);
		   /* dPosI.vtX = iWSp.lftpos; */
		   /* dPosI.vtY = iWSp.linpos; */
		   dGP       = gWSp.u_grintptr;
		   i_lineCom++;
		   i_columnP = iWSp.lftpos;

                   /* go to the correct column position */
		   if( i_lineCom == mItv->line  )
		   { 
		     /*
		     printf("\n");
		     printf("Line =:  %d\n",i_lineCom);
		     */
			     
		     i_prevC   = mItv->lftI;
		     while( (errNum == WLZ_ERR_NONE) && (i_columnP < mItv->lftI) )
		     {
		       i_columnP++;
	               switch(interp) 
	               {
		         case WLZ_INTERPOLATION_LINEAR:
		         case WLZ_INTERPOLATION_CLASSIFY_1:
		         case WLZ_INTERPOLATION_CALLBACK:



			 
	                 case WLZ_INTERPOLATION_NEAREST:
	                 case WLZ_INTERPOLATION_ORDER_2:
	                   /* get grey value   */
			   /*
			      give them a back ground value.
			    */
                           WlzGreyValueGetDir( gVWSp, 100000, 100000, 100000 );
	                   switch(gWSp.pixeltype)
	                   {
		             case WLZ_GREY_INT:
		               *(dGP.inp)++ = (*(gVWSp->gVal)).inv;
		               break;
		             case WLZ_GREY_SHORT:
		               *(dGP.shp)++ = (*(gVWSp->gVal)).shv;
		               break;
		             case WLZ_GREY_UBYTE:
		               *(dGP.ubp)++ = (*(gVWSp->gVal)).ubv;
		               break;
	  	             case WLZ_GREY_FLOAT:
		               *(dGP.flp)++ = (*(gVWSp->gVal)).flv;
		               break;
		             case WLZ_GREY_DOUBLE:
		               *(dGP.dbp)++ = (*(gVWSp->gVal)).dbv;
		               break;
		             case WLZ_GREY_RGBA:
		               *(dGP.rgbp)++ = (*(gVWSp->gVal)).rgbv;
		               break;
		             default:
		               errNum = WLZ_ERR_GREY_TYPE;
		             break;
	                   }
	                 break;
		       default:
		         errNum = WLZ_ERR_INTERPOLATION_TYPE;
			 break;
                       }
		     }
		   }
		 }
	         /* get the 3 points (xi,yi) and the 3 target points (x'i, y'i, z'i)  */
	         k         =   mItv->elmIdx;
		 if(k >= mSnWSp->nItvs)
		 {
		   printf("error elmIdx > nItvs \n");
		   exit(1);
		 }
		 /*
		    printf("elmIdx= %d\n", mItv->elmIdx);
		 */
	         k0        = ( wmt2D5->elements+k )->nodes[0];
	         k1        = ( wmt2D5->elements+k )->nodes[1];
	         k2        = ( wmt2D5->elements+k )->nodes[2];

	         sr0.vtX   = ( wmt2D5->nodes +k0 )->position.vtX ; 
	         sr0.vtY   = ( wmt2D5->nodes +k0 )->position.vtY ;
	         sr1.vtX   = ( wmt2D5->nodes +k1 )->position.vtX; 
	         sr1.vtY   = ( wmt2D5->nodes +k1 )->position.vtY;
	         sr2.vtX   = ( wmt2D5->nodes +k2 )->position.vtX; 
	         sr2.vtY   = ( wmt2D5->nodes +k2 )->position.vtY;

	         targ0.vtX = sr0.vtX + ( wmt2D5->nodes +k0 )->displacement.vtX;
	         targ0.vtY = sr0.vtY + ( wmt2D5->nodes +k0 )->displacement.vtY;
	         targ0.vtZ = (double)cutPosition + ( wmt2D5->nodes +k0 )->displacement.vtZ;
	         targ1.vtX = sr1.vtX + ( wmt2D5->nodes +k1 )->displacement.vtX;
	         targ1.vtY = sr1.vtY + ( wmt2D5->nodes +k1 )->displacement.vtY;
	         targ1.vtZ = (double)cutPosition + ( wmt2D5->nodes +k1 )->displacement.vtZ;
	         targ2.vtX = sr2.vtX + ( wmt2D5->nodes +k2 )->displacement.vtX;
	         targ2.vtY = sr2.vtY + ( wmt2D5->nodes +k2 )->displacement.vtY;
	         targ2.vtZ = (double)cutPosition + ( wmt2D5->nodes +k2 )->displacement.vtZ; 
					   
                 /* get the Affine transformation for this segment line */
                 WlzMakeAffine2D5With3PointsTrFn( sr0, sr1, sr2, targ0, targ1, targ2, 
			                    trans->mat);
		/*			    
	         if( (mItv->line == (mItv+1)->line) && ( i_prevC < mItv->rgtI ) )
		 {  */ 
		     /* move to the correct position */
		     /*
		       i_prevC =  mItv->rgtI;
		 } */
			    
	         /* within the scan line */
		 j = mItv->lftI;
		 if(j < i_prevC  )
		 {
		   j = i_prevC;
		 }
	         for(; j<= mItv->rgtI; j++)
	         {  
		   /* for debug */
		   /* printf("j=: %d\n", j); */
		    i_prevC++; 
	            /* get the source points using the Affine transform *
	               for clarity we keep it the clear form at the moment */

	            xs = (double)j          + trans->mat[0][0] * j + trans->mat[0][1] * mItv->line + trans->mat[0][2];
	            ys = (double)mItv->line + trans->mat[1][0] * j + trans->mat[1][1] * mItv->line + trans->mat[1][2];
	            zs =                      trans->mat[2][0] * j + trans->mat[2][1] * mItv->line + trans->mat[2][2];
		    /* change to the nearest integer */
		    if( xs >0 )
		    {
		       ix = (int)(xs + 0.5 );
		    }
		    else
		    {
		       ix = (int)(xs - 0.5 );
		    }
		    
		    if( ys >0 )
		    {
		       iy = (int)(ys + 0.5 );
		    }
		    else
		    {
		       iy = (int)(ys - 0.5 );
		    }
		    
		    if( zs >0 )
		    {
		       iz = (int)(zs + 0.5 );
		    }
		    else
		    {
		       iz = (int)(zs - 0.5 );
		    }
	            /* judge whether the source point is within the source object */
                    if(WlzInsideDomain(srcObj, iz, iy, ix, &errNum))
		    {
		      /* debug  */
		      /*
		      if(cutPosition == -161)
		      {
		        printf("iz iy ix = %d  %d  %d\n",iz,iy,ix);
		      }
		      */
		      WlzGreyValueGetDir( gVWSp, iz, iy, ix );
		      
		      /* debug  */
		      /*
		      if(cutPosition == -161)
		      {
		        printf("iz iy ix = %d  %d  %d\n",iz,iy,ix);
			exit(1);
		      }
		      */
		      
		    }
		    else
		    {
		     WlzGreyValueGetDir( gVWSp, 100000,100000,100000);
		    }

	            switch(interp)
	            {
		     	case WLZ_INTERPOLATION_LINEAR:
		        case WLZ_INTERPOLATION_CLASSIFY_1:
		        case WLZ_INTERPOLATION_CALLBACK:
	                case WLZ_INTERPOLATION_NEAREST:
			case WLZ_INTERPOLATION_ORDER_2:
	                  /* get grey value   */
			  /*
		          printf("xs: %f ys: %f zs %f\n",xs, ys, zs);
		          printf("xs: %d ys: %d zs %d\n",(int)xs, (int)ys, (int)zs);
			  */
                          /* WlzGreyValueGetDir( gVWSp, (int)zs, (int)ys, (int)xs ); */
		          /* assign the value to the target */
		          /* get the position to assign the grey value */
		          /* gk =  j + (mItv->line) *  (iWSp.rgtpos -  iWSp.lftpos + 1 ); */
	                  switch(gWSp.pixeltype)
	                  {
		             case WLZ_GREY_INT:
		               *(dGP.inp)++ = (*(gVWSp->gVal)).inv;
		               break;
		             case WLZ_GREY_SHORT:
		               *(dGP.shp)++ = (*(gVWSp->gVal)).shv;
		               break;
		             case WLZ_GREY_UBYTE:
		               *(dGP.ubp)++ = (*(gVWSp->gVal)).ubv;
		               break;
	  	             case WLZ_GREY_FLOAT:
		               *(dGP.flp)++ = (*(gVWSp->gVal)).flv;
		               break;
		             case WLZ_GREY_DOUBLE:
		               *(dGP.dbp)++ = (*(gVWSp->gVal)).dbv;
		               break;
		             case WLZ_GREY_RGBA:
		               *(dGP.rgbp)++ = (*(gVWSp->gVal)).rgbv;
		               break;
		             default:
		               errNum = WLZ_ERR_GREY_TYPE;
		             break;
	                  }
	                break;
		      default:
		        /* errNum = WLZ_ERR_INTERPOLATION_TYPE; unused */
			break;
                    }
 	        }
		/* ------check whether there is a gap here and filled it with back ground vlaues------- */
	        if(  (mItv->line == (mItv+1)->line)  && (mItv->rgtI < ( (mItv+1)->lftI ))  )
		{
		  /* fill the gape ! */
		  /*  debug
		  printf("pLeft=  %d   pRight=:   %d   p+1Left=: %d   p+1Right=:  %d\n",
		              mItv->lftI, mItv->rgtI, (mItv+1)->lftI, (mItv+1)->rgtI );
	          */		      
                    for(ip = i_prevC; ip <(mItv+1)->lftI; ip++)
		    {  /*  for debug */
		       /*
		       printf("ip=: %d\n", ip);
		       */
		       i_prevC++;
	               switch(interp)
	               {

		         case WLZ_INTERPOLATION_LINEAR:
		         case WLZ_INTERPOLATION_CLASSIFY_1:
		         case WLZ_INTERPOLATION_CALLBACK:
	                 case WLZ_INTERPOLATION_NEAREST:
			 case WLZ_INTERPOLATION_ORDER_2:
	                   /* get grey value   */
			   /*
			      give them a back ground value.
			    */
                           WlzGreyValueGetDir( gVWSp, 100000, 100000, 100000 ); 
	                   switch(gWSp.pixeltype)
	                   {
		             case WLZ_GREY_INT:
		               *(dGP.inp)++ = (*(gVWSp->gVal)).inv;
		               break;
		             case WLZ_GREY_SHORT:
		               *(dGP.shp)++ = (*(gVWSp->gVal)).shv;
		               break;
		             case WLZ_GREY_UBYTE:
		               *(dGP.ubp)++ = (*(gVWSp->gVal)).ubv;
		               break;
	  	             case WLZ_GREY_FLOAT:
		               *(dGP.flp)++ = (*(gVWSp->gVal)).flv;
		               break;
		             case WLZ_GREY_DOUBLE:
		               *(dGP.dbp)++ = (*(gVWSp->gVal)).dbv;
		               break;
		             case WLZ_GREY_RGBA:
		               *(dGP.rgbp)++ = (*(gVWSp->gVal)).rgbv;
		               break;
		             default:
		               errNum = WLZ_ERR_GREY_TYPE;
		             break;


			     
	                   }
	                 break;
		       default:
		         /* errNum = WLZ_ERR_INTERPOLATION_TYPE; unused */
			 break;
                       }
		    }
		}
                /* next one*/
	        ++mItv;
	      }
     if(errNum == WLZ_ERR_EOO)           /* Reset error from end of intervals */
     {
      errNum = WLZ_ERR_NONE;
     }
    }
   }
  }
 }
  /* output the 2D woolz object for debug purpose */
  debugNum = 0;
  if( debugNum == 1)
  {
    /* Write the 2D Woolz object for test */
    if((fp = fopen("Wlz2DObjout.wlz", "w")) == NULL )
    {
      printf("cannot open the output woolz file.\n");
      exit(1);
    }

    if( (errNum =	WlzWriteObj(fp, dstObj) ) != WLZ_ERR_NONE )
    {
      printf("output Woolz Object Error.\n");
      fclose(fp);
      exit(1);
    }
    fclose(fp);
    fp = NULL;
  }
  /* WlzMeshScanWSpFree(mSnWSp); */
  /*
  WlzGreyValueFreeWSp(gVWSp);
  */
  AlcFree(trans);
  return(errNum);
}



/*!
* - Function:    WlzGet2D5TrangularMeshFrom3DMesh
* - Ingroup:     WlzMesh
* - Returns:     none 
* - Purpose:     Get 2D5 triangular mesh. 
* - Global refs:  -     
* - Parameters: 
*      -#  *wmt2D5:               the 2D5 mesh transform.
*      -#  *wmt3D:                the 3D mesh transform.
*      -#   cutPosition:          plance to cut
*      -#  *errNum:               Woolz error number 
*      -#  *intersectIndex:       index which store the information which tetrahedron has been cut. 
*      -#  *planepoints:          the nodes of the cut surface. 
*      -#  *planepointsO:         the nodes of the original surface. 
*      -# **linkList:             array store the link information of the nodes. 
*      -#   nInterSectEsimate:    number of intersect guessed. 
*      -#   numEstTrangularsElem: number of triangular elements needed ( guessed). 
*      -#   numEstCutingPoints:   number of cuting points guessed. 
* - Author:       J. Rao, R. Baldock and B. Hill
*/
void  static WlzGet2D5TrangularMeshFrom3DMesh(  WlzMeshTransform2D5 *wmt2D5,
                                      WlzMeshTransform3D  *wmt3D, 
				      int cutPosition, 
				      WlzErrorNum *errNum,
				 int         *intersectIndex,
				 WlzDVertex3 *planepoints,
				 WlzDVertex3 *planepointsO,
				 int        **linkList,
                                 int          nInterSectEsimate,
				 int          numEstTrangularsElem,
				 int          numEstCutingPoints 
				      )
{
     int  nIntersect;
     int  i, j;
     int  noRedundancyCutingNum;
     int  numTotalTrangularElem;




     
     /* get the intersected plane */
     nIntersect            = 0;
     noRedundancyCutingNum = 0;
     numTotalTrangularElem = 0;
     /* Initialize the linkList */
     for(i=0; i<nInterSectEsimate; i++)
     {  for(j=0; j<5; j++)
        *(*(linkList +i )+j) = 0;
     }
     nIntersect = WlzIsoIntersectWithTetrahadronIndex(cutPosition, wmt3D, 
                                                      intersectIndex, planepoints, linkList,
                                                      &noRedundancyCutingNum, &numTotalTrangularElem ); 
     /* make the 2D5 mesh consistent with previous mesh data structure */

     /* get total number of triangle elements */
     if( numTotalTrangularElem  >  numEstTrangularsElem )
     {
      printf("memory allocated for number of Total Triangular Elements is not enough!");
      exit(1);
     }

     if( noRedundancyCutingNum  > numEstCutingPoints )
     {
      printf("memory allocated for number of Total cutting points is not enough!");
      exit(1);
     }
     
     wmt2D5->nElem   = numTotalTrangularElem;      /* total number of the tiangule elements */
     wmt2D5->nNodes  = noRedundancyCutingNum;      /* total number of 2D Nodes  */
     wmt2D5->zConst  = (double)cutPosition;

     /*
     printf("number of Elements %d  Layer  %d  nZconst %d number of points %d\n",
                     nzlayer, wmt2D5->nElem, (int)zConst, wmt2D5->nNodes);
     */

     /* Now established the 2D5 triangle mesh position without displacement */
     Make2D5TriangleMeshFromCuttedPlane(wmt2D5, planepoints, linkList, nIntersect, noRedundancyCutingNum);

     /* output for test wmt2D5 */
     /* open an file for output 
        This plane is consisted by the tetrahedrons which were cutted by the given z-const plane */
    /* transform the cutting plane points to the corresponding points in the source plane or 
        simply by get its displacements */
     GetTheSourceSurfaceOfTheCutPlanePoints( nIntersect, intersectIndex,   wmt3D, 
                                             planepoints, planepointsO,  linkList
                                           );
					   
     Make2D5MeshDisplacement(wmt2D5, planepointsO); 
     /* open an file for output */
}

/*!
* - Function:    WlzIsSinglePointCutMesh
* - Returns:     an integer >0 Yes. =0 No.
* - Ingroup:     WlzAccess
* - Purpose:     judge whether it is a single point cut. 
* - Global refs:  -     
* - Parameters: 
*      -#  *wmt2D5:               the 2D5 mesh transform.
* - Author:       J. Rao, R. Baldock and B. Hill
*/
int     WlzIsSinglePointCutMesh(const WlzMeshTransform2D5 *wmt2D5)
{ int i,k;
  WlzDVertex2 v2;
  k = 1;
  /* pick up one */
  v2.vtX = (wmt2D5->nodes)->position.vtX;
  v2.vtY = (wmt2D5->nodes)->position.vtY;
   
  for(i=1; i<wmt2D5->nElem; i++)
  {
  
    if( (  ( fabs( v2.vtX - (wmt2D5->nodes+i)->position.vtX )) > 0.01 )  || 
        (  ( fabs( v2.vtY - (wmt2D5->nodes+i)->position.vtY )) > 0.01 )     ) 
    {
      k = 0;
      break;
    }
  }
  return k;
}
/*!
* - Function:    WlzNoAreaGreaterThanOne
* - Returns:     an integer >0 Yes. =0 No.
* - Ingroup:     WlzAccess
* - Purpose:     judge whether there is no area greater than one. 
* - Global refs:  -     
* - Parameters: 
*      -#  *wmt2D5:               the 2D5 mesh transform.
* - Author:       J. Rao, R. Baldock and B. Hill
*/
int     WlzNoAreaGreaterThanOne(const WlzMeshTransform2D5 *wmt2D5)
{ int i, l,m,n, k;
  WlzDVertex2 v2[3];
  double area;
  k = 1;
    for(i=0;i<wmt2D5->nElem; i++)
    {
      l        =  (wmt2D5->elements+i)->nodes[0];
      m        =  (wmt2D5->elements+i)->nodes[1];
      n        =  (wmt2D5->elements+i)->nodes[2];
     v2[0].vtX =  (wmt2D5->nodes+l)->position.vtX;
     v2[1].vtX =  (wmt2D5->nodes+m)->position.vtX;
     v2[2].vtX =  (wmt2D5->nodes+n)->position.vtX;
     v2[0].vtY =  (wmt2D5->nodes+l)->position.vtY;
     v2[1].vtY =  (wmt2D5->nodes+m)->position.vtY;
     v2[2].vtY =  (wmt2D5->nodes+n)->position.vtY;
     area = WlzGeomTriangleSnArea2(v2[0],v2[1],v2[2]);
     /*
     printf("n1 %d n2 %d n3 %d\n", l, m ,n);
     printf("area: =  %f\n", area);
     */
     if(  fabs( area ) >2.0 )
     {
        k = 0;
        break;
     }
    }
  return k;
}

/*!
*   - Return:   None
*   - Ingroup:  WlzTransform
*   - Purpose:  Get 2D5 Affine transformation from the nodes of
*               a source trangle(in x-y plane) and the nodes from a 
*	       target trangle(in 3D plane).
*  
*   - Transform:  (x,y)  =>  (x', y', z')
*
*   - Input:
*          -#    source tetrahedron: 3 points   sr1,   sr2,   sr3
*          -#    target tetrahedron: 3 points targ1, targ2, targ3
*   - Output:
*          -#    Affine2D5With3pointsTrFun[3][3] contains the 
*	        affine transformation parameters:
*
*   - A =:
*   -     a11  a12  a13
*   -     a21  a22  a23
*   -     a31  a32  a33
*	      
*
*      - This coefficients can be used to get (x',y',z') by
*
*      -  x' = x + a11*x + a12*y + a13
*      -  y' = y + a21*x + a22*y + a23
*      -  z' =     a31*x + a22*y + a33
* - Author:       J. Rao, R. Baldock and B. Hill
*/
/*
               by   J.Rao   06/11/2001 
*/
void static WlzMakeAffine2D5With3PointsTrFn( WlzDVertex2 sr0, 
                                 WlzDVertex2 sr1, 
				 WlzDVertex2 sr2,
				 WlzDVertex3 targ0,
				 WlzDVertex3 targ1, 
                                 WlzDVertex3 targ2, 
			         double **Affine2D5With3PointsTrFun )
{
  int            i, j, idN;
  double	*bV = NULL,
  		*wV = NULL;
  double       **aA,
               **vA;
  double        wMax;	       
  int           nSys = 3;
  int           thresh;
  AlgMatrix	aM,
  		vM;
  const double  tol  = 1.0e-06;
  WlzErrorNum	errNum = WLZ_ERR_NONE;


  aM.core = NULL;
  vM.core = NULL;
  /*------- allocate memory --------*/
  if (
      ((wV = (double *)AlcCalloc(sizeof(double), nSys))  == NULL) ||
      ((bV = (double *)AlcMalloc(sizeof(double) * nSys)) == NULL) ||
      ((aM.rect = AlgMatrixRectNew(nSys, nSys, NULL)) == NULL) ||
      ((vM.rect = AlgMatrixRectNew(nSys, nSys, NULL)) == NULL))
  {
      errNum = WLZ_ERR_MEM_ALLOC;
  }
  
  if(errNum == WLZ_ERR_NONE)
  {
    aA = aM.rect->array;
    vA = vM.rect->array;
    /* fill the a-matrix  where a is:
   
      x0 y0  1
      x1 y1  1
      x2 y2  1
   
    where (x1, y1) is the source node one and etc.
   
    */

    *(*(aA + 0 ) + 0 ) = sr0.vtX;
    *(*(aA + 0 ) + 1 ) = sr0.vtY;
    *(*(aA + 0 ) + 2 ) = 1.0;

    *(*(aA + 1 ) + 0 ) = sr1.vtX;
    *(*(aA + 1 ) + 1 ) = sr1.vtY;
    *(*(aA + 1 ) + 2 ) = 1.0;

    *(*(aA + 2 ) + 0 ) = sr2.vtX;
    *(*(aA + 2 ) + 1 ) = sr2.vtY;
    *(*(aA + 2 ) + 2 ) = 1.0;

    /* get a11 a12  and a13 */
    
    /* fill the u */
    *(bV + 0 ) = targ0.vtX - sr0.vtX;
    *(bV + 1 ) = targ1.vtX - sr1.vtX;
    *(bV + 2 ) = targ2.vtX - sr2.vtX;

    /* */
    for(i=0;i<3; i++)
    {
      *(wV+i) = 0.0;
      for(j=0;j<3;j++)
      {
        *(*(vA+i)+j) = 0.0;
      }
    }
   
    /* Perform singular value decomposition of matrix a. */
   
    errNum = WlzErrorFromAlg(AlgMatrixSVDecomp(aM, wV, vM));
    /* */
    /*  
      for(i=0;i<4;i++)
      {
        printf("%f\n",*(wV+i));
        for(j=0;j<4;j++)
        {
          printf("%f   ",*(*(vA+i) +j) );
          printf("\n");
        }
       }
    */  
 
    if(errNum != WLZ_ERR_NONE)
    {
       printf("failed to make decompositions ");
       exit(1);
    }
    {
      /* Edit the singular values. */
      wMax = 0.0;
      for(idN = 0; idN < nSys; ++idN)
      {
        if(*(wV + idN) > wMax)
        {
	  wMax = *(wV + idN);
        }
      }
      thresh = tol * wMax;
      for(idN = 0; idN < nSys; ++idN)
      {
        if(*(wV + idN) < thresh)
        {
	  *(wV + idN) = 0.0;
        }
      }
      /* Solve for the X conponents  */
      errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
    }

    if(errNum != WLZ_ERR_NONE)
    {
       printf("failed to get Affine transform for x-component");
       exit(1);
    }
    {
      /* Store the a11, a12, a13 */
       *(*(Affine2D5With3PointsTrFun + 0 )+ 0 ) = *(bV + 0);
       *(*(Affine2D5With3PointsTrFun + 0 )+ 1 ) = *(bV + 1);
       *(*(Affine2D5With3PointsTrFun + 0 )+ 2 ) = *(bV + 2);
     
      /* Now set for y  conponents */
       *(bV + 0 ) = targ0.vtY - sr0.vtY;
       *(bV + 1 ) = targ1.vtY - sr1.vtY;
       *(bV + 2 ) = targ2.vtY - sr2.vtY;
   
      /* Solve for the Y conponents */
      errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
    }

    if(errNum != WLZ_ERR_NONE)
    {
       printf("failed to get Affine transform for y-component");
       exit(1);
    }
 
    {
      /* Store the a21, a22, a23 */
       *(*(Affine2D5With3PointsTrFun + 1 )+ 0 ) = *(bV + 0);
       *(*(Affine2D5With3PointsTrFun + 1 )+ 1 ) = *(bV + 1);
       *(*(Affine2D5With3PointsTrFun + 1 )+ 2 ) = *(bV + 2);
     
      /* Now set for z  conponents */
       *(bV + 0 ) = targ0.vtZ;
       *(bV + 1 ) = targ1.vtZ;
       *(bV + 2 ) = targ2.vtZ;
   
      /* Solve for the Z conponents */
       errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
    }

    if(errNum != WLZ_ERR_NONE)
    {
       printf("failed to get Affine transform for z-component");
       exit(1);
    }
    {
      /* Store the a31, a32, a33 */
       *(*(Affine2D5With3PointsTrFun + 2 )+ 0 ) = *(bV + 0);
       *(*(Affine2D5With3PointsTrFun + 2 )+ 1 ) = *(bV + 1);
       *(*(Affine2D5With3PointsTrFun + 2 )+ 2 ) = *(bV + 2);
    }
    AlcFree(bV);
    AlcFree(wV);
    AlgMatrixFree(aM);
    AlgMatrixFree(vM);
  }
 
}





/*!
* \ingroup	WlzTransform
* \return				Transformed plane domain,
*					NULL on error.
* \brief	Creates a new plane domain from a given bounding box
* \param	srcObj			Given source object.
* \param	trans			Given affine transform.
* \param	dstErr			Destination ptr for Woolz error
*					code, may be NULL.
*/
#ifdef DEBUGME
static WlzPlaneDomain *WlzPDomFromBBox(WlzIBox3  bBox, 
				       WlzErrorNum *dstErr)
{
  WlzObject      *tObj;
  /* WlzUByte	**pMsk = NULL; */
  UBYTE         **pMsk = NULL;
  WlzIVertex2	  pMskOrg,
  		  pMskSz;
  WlzPlaneDomain *dstPDom = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;


  WlzIVertex3	  idx,
  		  dPos;

           if(errNum == WLZ_ERR_NONE)
           {


            /* Create a plane domain and a destination object. */
            dstPDom = WlzMakePlaneDomain( WLZ_PLANEDOMAIN_DOMAIN,
	                                  bBox.zMin, bBox.zMax,
	                                  bBox.yMin, bBox.yMax,
	                                  bBox.xMin, bBox.xMax,
	 					       &errNum );


	      pMskOrg.vtX = bBox.xMin;
              pMskOrg.vtY = bBox.yMin;
              pMskSz.vtX = bBox.xMax - bBox.xMin + 1;
              pMskSz.vtY = bBox.yMax - bBox.yMin + 1;
	   
              if(AlcBit2Calloc(&pMsk, pMskSz.vtY+1, pMskSz.vtX) != ALC_ER_NONE)
              {
                errNum = WLZ_ERR_MEM_ALLOC;
              }
           }
	   if(errNum == WLZ_ERR_NONE)
	   {
	       idx.vtZ  = 0;
               dPos.vtZ = bBox.zMin;
               while((errNum == WLZ_ERR_NONE) && (dPos.vtZ <= bBox.zMax))
               {
                  /* Set the plane bit mask. */
                  idx.vtY  = 0;
                  dPos.vtY = bBox.yMin;
                  while(dPos.vtY <= bBox.yMax)
                  {
 	            idx.vtX  = 0;
	            dPos.vtX = bBox.xMin;
	            while(dPos.vtX <= bBox.xMax)
	            {
		    	++(idx.vtX);
	                ++(dPos.vtX);
		    }
		    ++(idx.vtY);
         	    ++(dPos.vtY);
		  }
	          /* Create a domain from the bit mask. */
                  tObj = WlzFromArray2D((void **)pMsk, pMskSz, pMskOrg, WLZ_GREY_BIT,
      			                WLZ_GREY_BIT, 0.0, 1.0, 0, 0, &errNum);
		  if(tObj)
		  {
		    *(dstPDom->domains + idx.vtZ) = WlzAssignDomain(tObj->domain, NULL);
                     (void )WlzFreeObj(tObj);
	             tObj = NULL;
		  }
		  ++(idx.vtZ);
                  ++(dPos.vtZ);
	       }
	   }

  if(pMsk)
  {
    Alc2Free((void **)pMsk);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzStandardPlaneDomain(dstPDom, NULL);
  }
  else
  {
    if(dstPDom)
    {
      (void )WlzFreePlaneDomain(dstPDom);
    }
    dstPDom = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstPDom);

}
#endif


/*!
* - Function:    WlzGet2D5Transform
* - Returns:     none 
* - InGroup      WlzMesh
* - Purpose:     Get 2D5 triangular mesh. 
* - Global refs:  -     
* - Parameters: 
*      -#  *wmt2D5:               the 2D5 mesh transform.
*      -#  *wmt3D:                the 3D mesh transform.
*      -#   cutPosition:          plance to cut
*      -#  *dsterrNum:               Woolz error number 
* - Author:       J. Rao, R. Baldock and B. Hill
*/
void WlzGet2D5Transform( const WlzMeshTransform3D   *wmt3D,
                      int   cutPosition,
                      WlzMeshTransform2D5        *wmt2D5,
                      WlzErrorNum                *dsterrNum
                    )
{ int  i,j;
  int  nInterSectEsimate;
  int  nIntersect;
  int  noRedundancyCutingNum;
  int  numTotalTrangularElem;
  int *intersectIndex   = NULL;
  int                **linkList       = NULL;   /* store the linked list
                                                   the first index linked with
		 		                   the index of cutting tetrahedron 
				                   The second one linked with 
				                   the number of points cuttied to the tetrahedron */
  WlzErrorNum errNum;
  WlzDVertex3 *planepoints  = NULL;
  WlzDVertex3 *planepointsO = NULL;


  
  errNum                = *dsterrNum;
  nIntersect            = 0;
  noRedundancyCutingNum = 0;
  numTotalTrangularElem = 0;


  nInterSectEsimate = wmt3D->nElem*6;

  if( ( intersectIndex = (int *)AlcRealloc( intersectIndex, (nInterSectEsimate) * sizeof(int) )  ) == NULL)
  {
       errNum = WLZ_ERR_MEM_ALLOC;
       printf("Failed to allocate memory for intersectIdex!");
       exit(1);
  }

  for(i=0; i<nInterSectEsimate; i++)
  {
    *(intersectIndex+i )= -10;
  }

  /* allocate memory to store the cutting points */
  if( ( planepoints = (WlzDVertex3 *)AlcRealloc(planepoints, 4 * nInterSectEsimate * sizeof(WlzDVertex3) )  ) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
    printf("Failed to allocate memory for planepoints!");
    exit(1);
  }

  /* allocate memory to store the cutting points */
  if( ( planepointsO = (WlzDVertex3 *)AlcRealloc(planepointsO, 4 * nInterSectEsimate * sizeof(WlzDVertex3) )  ) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
    printf("Failed to allocate memory for planepoints!");
    exit(1);
  }

  /* allocate memory to store the linked list of the cutting points and the index the cutted tetrahedron */
  if( (AlcInt2Malloc(&linkList, nInterSectEsimate,5) !=  ALC_ER_NONE)  )
  {
    errNum = WLZ_ERR_MEM_ALLOC;
    printf("Failed to allocate memory for link list!");
    exit(1);
  }

  /* Initialize the linkList */
  for(i=0; i<nInterSectEsimate; i++)
  {  for(j=0; j<5; j++)
        {
          *(*(linkList +i )+j) = 0;
	}  
  }
  nIntersect = WlzIsoIntersectWithTetrahadronIndex( (double)cutPosition, wmt3D, 
                                                     intersectIndex, planepoints, linkList,
                                                    &noRedundancyCutingNum, &numTotalTrangularElem ); 
     
  wmt2D5->nElem   = numTotalTrangularElem;      /* total number of the tiangule elements */
  wmt2D5->nNodes  = noRedundancyCutingNum;      /* total number of 2D Nodes  */
  wmt2D5->zConst  = (double)cutPosition;

  /* initialize the wmt2D5 */
  WlzInitializeWmt2D5(wmt2D5);

  /* Now established the 2D5 triangle mesh  */
  Make2D5TriangleMeshFromCuttedPlane(wmt2D5, planepoints, linkList, nIntersect, noRedundancyCutingNum);

  /* transform the cutting plane points to the corresponding points in the source plane or 
     simply by get its displacements */
  /*   */
  GetTheSourceSurfaceOfTheCutPlanePoints( nIntersect, intersectIndex,   wmt3D, 
                                          planepoints, planepointsO,  linkList
                                        );

  Make2D5MeshDisplacement(wmt2D5, planepointsO); 
 
  /*  free the memory */
  AlcFree(planepoints);
  AlcFree(planepointsO);
  (void)AlcInt2Free(linkList);
  AlcFree(intersectIndex);
  /*
    planepoints    = NULL;
    planepointsO   = NULL;
    linkList       = NULL;
    intersectIndex = NULL;
  */
  if(dsterrNum)
  {
    *dsterrNum = errNum;
  }
}

/*!
* \return   None
* \ingroup  WlzMesh
* \brief    initialize the 2D5 mesh.
* \param    *wmt2D5:               the 2D5 mesh transform.
* \author:       J. Rao, R. Baldock and B. Hill
*/
void   WlzInitializeWmt2D5(WlzMeshTransform2D5  *wmt2D5)
{ int i;
   
  for(i=0; i<wmt2D5->nElem; i++)
  {
    (wmt2D5->nodes+i)->position.vtX     = 0;
    (wmt2D5->nodes+i)->position.vtY     = 0;
    (wmt2D5->nodes+i)->displacement.vtX = 0;
    (wmt2D5->nodes+i)->displacement.vtY = 0;
    (wmt2D5->nodes+i)->displacement.vtZ = 0;
  }

  for(i=0; i<wmt2D5->nNodes; i++)
  {
    (wmt2D5->elements+i)->nodes[0]   = -10;
    (wmt2D5->elements+i)->nodes[1]   = -10;
    (wmt2D5->elements+i)->nodes[2]   = -10;
    (wmt2D5->elements+i)->neighbours[0] = -10;
    (wmt2D5->elements+i)->neighbours[1] = -10;
    (wmt2D5->elements+i)->neighbours[2] = -10;
  }


}

/*!
* \return   None 
* \ingroup  WlzAccess
* \brief    Get the index of the tetrahedrons which intersect
*               with a giving z = constant plane. 
* \param    zConst           z = const plane.
* \param    fp               pointer pointing to a specific file.
* \param    wmt3D            mesh transform.
* \author:       J. Rao, R. Baldock and B. Hill
*/
#ifdef DEBUGME
void  static WlzIsoIntersectWithTetrahadronIndexAnDETC(  double             zConst, 
                                                const WlzMeshTransform3D   *wmt3D, 
					        WlzMeshTransform2D5        *wmt2D5
					     )
{
  int      i, j0, j1, j2, j3, k, k1, itemp;
  int      nIntersect = 0;
  int      nUp;
  int      iup, idown, icount, numberOfTriangular=0, isnumber=0;
  int      n0up, n1up, n2up, n3up;
  double   z0, z1, z2, z3;
 

 /* WlzDVertex3 planepoints[3000]; */  /* store the cutting points */
 /* int         linkList[3000][5];   */   /* store the linked list
                                      the first index linked with
				      the index of cutting tetrahedron 
				      The second one linked with 
				      the number of points cuttied to the tetrahedron */
  WlzDVertex3 wv3up[4];
  WlzDVertex3 wv3down[4];
  WlzDVertex3 tempw;

  int       idex;

  int                  *intersectIndex;
  int                 **linkList;
  /*
  int                  *noRedundancyCutingNum;
  int                  *numTotalTrangularElem;
  */
  WlzErrorNum          errNum;
  
  int nInterSectEsimate;

  WlzDVertex3       *planepoints  = NULL;
  WlzDVertex3       *planepointsO = NULL;
  int  j;
 
/* allocate memory here */

  nInterSectEsimate = wmt3D->nElem*6;

  if( ( intersectIndex = (int *)AlcRealloc( intersectIndex, (nInterSectEsimate) * sizeof(int) )  ) == NULL)
  {
       errNum = WLZ_ERR_MEM_ALLOC;
       printf("Failed to allocate memory for intersectIdex!");
       exit(1);
  }

  for(i=0; i<nInterSectEsimate; i++)
  {
    *(intersectIndex+i )= -10;
  }

  /* allocate memory to store the cutting points */
  if( ( planepoints = (WlzDVertex3 *)AlcRealloc(planepoints, 4 * nInterSectEsimate * sizeof(WlzDVertex3) )  ) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
    printf("Failed to allocate memory for planepoints!");
    exit(1);
  }

  /* allocate memory to store the cutting points */
  if( ( planepointsO = (WlzDVertex3 *)AlcRealloc(planepointsO, 4 * nInterSectEsimate * sizeof(WlzDVertex3) )  ) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
    printf("Failed to allocate memory for planepoints!");
    exit(1);
  }

  /* allocate memory to store the linked list of the cutting points and the index the cutted tetrahedron */
  if( (AlcInt2Malloc(&linkList, nInterSectEsimate,5) !=  ALC_ER_NONE)  )
  {
    errNum = WLZ_ERR_MEM_ALLOC;
    printf("Failed to allocate memory for link list!");
    exit(1);
  }


  /* allocate memory to store the cutting points */
  if( ( planepoints = (WlzDVertex3 *)AlcRealloc(planepoints, 4 * nInterSectEsimate * sizeof(WlzDVertex3) )  ) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
    printf("Failed to allocate memory for planepoints!");
    exit(1);
  }
   
  /* allocate memory to store the cutting points in the original planes */
  if( ( planepointsO = (WlzDVertex3 *)AlcRealloc(planepointsO, 4*nInterSectEsimate*sizeof(WlzDVertex3) )  ) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
    printf("Failed to allocate memory for planepointsO!");
    exit(1);
  }

  /* Initialize the linkList */
  for(i=0; i<nInterSectEsimate; i++)
  {  for(j=0; j<5; j++)
        {
          *(*(linkList +i )+j) = 0;
	}  
  }

  icount = 0;  /* used to count how many points cutted (no redundancy counting) */

  for (i = 0; i < wmt3D->nElem; i++)
  {   /*  get the z-position of the four points relative to the zConst 
          Remember that we are cutting the transformed mesh */
      /* Also, We only cut the mesh belong to the object! */
      /* exclude the tetrahedrons not belong to the object */

         j0    =  (wmt3D->elements+i )->nodes[0];
         j1    =  (wmt3D->elements+i )->nodes[1];
         j2    =  (wmt3D->elements+i )->nodes[2];
         j3    =  (wmt3D->elements+i )->nodes[3];
    
     if( (j0 != 0) || (j1 !=0 ) || (j2 !=0) || (j3 !=0 ) )
     {

         z0    =  (wmt3D->nodes+j0)->position.vtZ + (wmt3D->nodes+j0)->displacement.vtZ   - zConst;
         z1    =  (wmt3D->nodes+j1)->position.vtZ + (wmt3D->nodes+j1)->displacement.vtZ   - zConst;
         z2    =  (wmt3D->nodes+j2)->position.vtZ + (wmt3D->nodes+j2)->displacement.vtZ   - zConst;
         z3    =  (wmt3D->nodes+j3)->position.vtZ + (wmt3D->nodes+j3)->displacement.vtZ   - zConst;
         nUp   =  0;
         n0up  =  n1up = n2up = n3up = 0;
         /* printf("%lf  %lf  %lf  %lf\n", z0, z1, z2, z3); */

         /* judge whether there is a intersection with zConst-plane 
         only 3 cases if intersected

	 1) one   point up,  three points down 
	 2) two   point up,  two   points down 
	 3) three point up,  one   points down 
	
         Not intersected only two cases:
	 4) four points all up
	 5) four points all down
         */

         if( z0 > 0. )
         {
            n0up = 1;
            nUp++;
         }
     
         if( z1 > 0.  )
         {
            n1up = 1;
             nUp++;
         }
     
         if( z2 > 0.  )
         {
            n2up = 1;
             nUp++;
         }
     
         if( z3 > 0.  )
         {
            n3up = 1;
             nUp++;
         }
      
         /* Now the cutted cases  */
         if( (nUp != 0) && (nUp != 4) )
         {
           iup = idown = 0;
           /* divide the nodes of the cutted tetrahedron into up and down nodes then store them. */
           if(n0up == 1)
	   {
            wv3up[iup].vtX = (wmt3D->nodes+j0)->position.vtX + (wmt3D->nodes+j0)->displacement.vtX; 
            wv3up[iup].vtY = (wmt3D->nodes+j0)->position.vtY + (wmt3D->nodes+j0)->displacement.vtY; 
            wv3up[iup].vtZ = (wmt3D->nodes+j0)->position.vtZ + (wmt3D->nodes+j0)->displacement.vtZ;
	    iup++;
	   }
	   else if(n0up == 0)
	   {
            wv3down[idown].vtX = (wmt3D->nodes+j0)->position.vtX + (wmt3D->nodes+j0)->displacement.vtX; 
            wv3down[idown].vtY = (wmt3D->nodes+j0)->position.vtY + (wmt3D->nodes+j0)->displacement.vtY; 
            wv3down[idown].vtZ = (wmt3D->nodes+j0)->position.vtZ + (wmt3D->nodes+j0)->displacement.vtZ; 
	    idown++;
 	   }
	  
           if(n1up == 1)
	   {
            wv3up[iup].vtX = (wmt3D->nodes+j1)->position.vtX + (wmt3D->nodes+j1)->displacement.vtX; 
            wv3up[iup].vtY = (wmt3D->nodes+j1)->position.vtY + (wmt3D->nodes+j1)->displacement.vtY; 
            wv3up[iup].vtZ = (wmt3D->nodes+j1)->position.vtZ + (wmt3D->nodes+j1)->displacement.vtZ;
	    iup++;
	   }
	   else if(n1up == 0)
	   {
            wv3down[idown].vtX = (wmt3D->nodes+j1)->position.vtX + (wmt3D->nodes+j1)->displacement.vtX; 
            wv3down[idown].vtY = (wmt3D->nodes+j1)->position.vtY + (wmt3D->nodes+j1)->displacement.vtY; 
            wv3down[idown].vtZ = (wmt3D->nodes+j1)->position.vtZ + (wmt3D->nodes+j1)->displacement.vtZ; 
	    idown++;
 	   }
		
           if(n2up == 1)
	   {
            wv3up[iup].vtX = (wmt3D->nodes+j2)->position.vtX + (wmt3D->nodes+j2)->displacement.vtX; 
            wv3up[iup].vtY = (wmt3D->nodes+j2)->position.vtY + (wmt3D->nodes+j2)->displacement.vtY; 
            wv3up[iup].vtZ = (wmt3D->nodes+j2)->position.vtZ + (wmt3D->nodes+j2)->displacement.vtZ; 
	    iup++;
	   }
	   else if(n2up == 0)
	   {
            wv3down[idown].vtX = (wmt3D->nodes+j2)->position.vtX + (wmt3D->nodes+j2)->displacement.vtX; 
            wv3down[idown].vtY = (wmt3D->nodes+j2)->position.vtY + (wmt3D->nodes+j2)->displacement.vtY; 
            wv3down[idown].vtZ = (wmt3D->nodes+j2)->position.vtZ + (wmt3D->nodes+j2)->displacement.vtZ; 
	    idown++;
 	   }
		  
           if(n3up == 1)
	   {
            wv3up[iup].vtX = (wmt3D->nodes+j3)->position.vtX + (wmt3D->nodes+j3)->displacement.vtX; 
            wv3up[iup].vtY = (wmt3D->nodes+j3)->position.vtY + (wmt3D->nodes+j3)->displacement.vtY; 
            wv3up[iup].vtZ = (wmt3D->nodes+j3)->position.vtZ + (wmt3D->nodes+j3)->displacement.vtZ; 
	    iup++;
	   }
	   else if(n3up == 0)
	   {
            wv3down[idown].vtX = (wmt3D->nodes+j3)->position.vtX + (wmt3D->nodes+j3)->displacement.vtX; 
            wv3down[idown].vtY = (wmt3D->nodes+j3)->position.vtY + (wmt3D->nodes+j3)->displacement.vtY; 
            wv3down[idown].vtZ = (wmt3D->nodes+j3)->position.vtZ + (wmt3D->nodes+j3)->displacement.vtZ; 
	    idown++;
 	   }

           /* Now get the intersect points */
	   /* for output test, we always need an sequence to make a good output 
	     for vtk visualization 
	    
	     18/10/2001 but now we need to remove the redundancy storage */
	     tempw.vtZ = 0.;
           /*  ----- Now three points cut case -----  */
	   if( nUp == 1 || nUp == 3  )
	   {
	    isnumber = 1;
            for(k = 0; k < nUp; k++) 
	    { 
 
	        for( k1=0; k1< 4 - nUp; k1++ )
		{
		   if( wv3up[k].vtZ  - wv3down[k1].vtZ  == 0. )
		     { printf("the line in the z-pline");
		       exit(1); 
		     }
                    /* Get the cutting points */
                    tempw.vtX  =  wv3down[k1].vtX +   
		                ( wv3up[k].vtX - wv3down[k1].vtX ) * ( zConst - wv3down[k1].vtZ )
		               /( wv3up[k].vtZ - wv3down[k1].vtZ ); 
			
                    tempw.vtY  =  wv3down[k1].vtY  + 
		                ( wv3up[k].vtY - wv3down[k1].vtY ) * ( zConst - wv3down[k1].vtZ )
		               /( wv3up[k].vtZ - wv3down[k1].vtZ );
			      
	            /* Ask whether this point is already stored 
		       if not return idex < 0 otherwise return the
		       index of the existing point already stored */
		    idex = -1;
                    idex = storedPointIndex(tempw, planepoints, icount);
		 
                    if(idex < 0)
		    {
		     /* store the points */
                     (planepoints+icount)->vtX  =  tempw.vtX;
                     (planepoints+icount)->vtY  =  tempw.vtY; 
                     (planepoints+icount)->vtZ  =  zConst;
		   
	             /* record the index of nodes for this cut */
               	     *( *(linkList+nIntersect) + isnumber ) = icount;
                     isnumber++;
     	             icount++;  
		    }
		    else if(idex >0)
		    { /* the point is already stored before
		       we only need to store the index */
	             /* record the index of nodes for this cut */
               	     *( *(linkList+nIntersect) + isnumber ) = idex;
                     isnumber++;
		    }
	         }
	     }	
	     /* now make sure it is counterclockwise odrdered */
	     if( WlzGeomTriangle2SnArea3(  planepoints,  *( *(linkList+nIntersect) + 1 ),
	                                       *( *(linkList+nIntersect) + 2 ),  
		                               *( *(linkList+nIntersect) + 3)
                             ) < 0.       )
	     { /* change the position of the last two nodes data structure
	         to make it CCW order */
               itemp =  *( *(linkList+nIntersect) + 2);
               *( *(linkList+nIntersect) + 2) = *( *(linkList+nIntersect) + 3);
	       *( *(linkList+nIntersect) + 3) = itemp;
	     }
	  }

	  /* ------ four points cutting case ------- */
	  else if(nUp == 2)
	  {
	    isnumber = 1;
	    /* k = 0 */
	    k = 0; 
	    for( k1=0; k1< 2; k1++ )
	    {
	       if( wv3up[k].vtZ  - wv3down[k1].vtZ  == 0. )
		{ printf("the line in the z-pline");
		       exit(1); 
		}
                /* Get the cutting points */
                tempw.vtX  =  wv3down[k1].vtX +   
		              ( wv3up[k].vtX - wv3down[k1].vtX ) * ( zConst - wv3down[k1].vtZ )
		             /( wv3up[k].vtZ - wv3down[k1].vtZ ); 
			
                tempw.vtY  =  wv3down[k1].vtY  + 
		              ( wv3up[k].vtY - wv3down[k1].vtY ) * ( zConst - wv3down[k1].vtZ )
		             /( wv3up[k].vtZ - wv3down[k1].vtZ ); 

	        /* Ask whether this point is already stored 
		   if not return idex < 0 otherwise return the
		   index of the existing point already stored */
		idex = -1;
                idex = storedPointIndex(tempw, planepoints, icount);
                if(idex < 0)
		{
		   /* store the points */
                   (planepoints+icount)->vtX  =  tempw.vtX;
                   (planepoints+icount)->vtY  =  tempw.vtY; 
                   (planepoints+icount)->vtZ  =  zConst;
		   
	           /* record the index of nodes for this cut */
               	   *( *(linkList+nIntersect) + isnumber ) = icount;
                   isnumber++;
     	           icount++;  
		 }
		 else if(idex >0)
		 { /* the point is already stored before
		       we only need to store the index */
	           /* record the index of nodes for this cut */
               	   *( *(linkList+nIntersect) + isnumber ) = idex;
                   isnumber++;
		 }
	      }
		
	      /* k = 1 */
	      k = 1; 
	      for( k1= 1;  k1 >= 0; k1-- )
	      {
		   if( wv3up[k].vtZ  - wv3down[k1].vtZ  == 0. )
		     { printf("the line in the z-pline");
		       exit(1); 
		     }

                     tempw.vtX  =  wv3down[k1].vtX +   
		                 ( wv3up[k].vtX - wv3down[k1].vtX ) * ( zConst - wv3down[k1].vtZ )
		                /( wv3up[k].vtZ - wv3down[k1].vtZ ); 
			
                     tempw.vtY  =  wv3down[k1].vtY  + 
		                 ( wv3up[k].vtY - wv3down[k1].vtY ) * ( zConst - wv3down[k1].vtZ )
		                /( wv3up[k].vtZ - wv3down[k1].vtZ ); 

	             /* Ask whether this point is already stored 
		        if not return idex < 0 otherwise return the
		        index of the existing point already stored */
		     idex = -1;
                     idex = storedPointIndex(tempw, planepoints, icount);
                     if(idex < 0)
		     {
		       /* store the points */
                       (planepoints+icount)->vtX  =  tempw.vtX;
                       (planepoints+icount)->vtY  =  tempw.vtY; 
                       (planepoints+icount)->vtZ  =  zConst;
		   
	               /* record the index of nodes for this cut */
               	       *( *(linkList+nIntersect) + isnumber ) = icount;
                       isnumber++;
     	               icount++;  
		      }
		      else if(idex >0)
		      { /* the point is already stored before
		           we only need to store the index */
	                /* record the index of nodes for this cut */
               	        *( *(linkList+nIntersect) + isnumber ) = idex;
                        isnumber++;
		      }
		}
	        /* now make sure that this quarilateral divided into two CCW 
	            triangles */
	        if( WlzGeomTriangle2SnArea3(  planepoints,  *( *(linkList+nIntersect) + 1 ),
	                                          *( *(linkList+nIntersect) + 2),  
		                                  *( *(linkList+nIntersect) + 3)
                                   ) < 0.       )
	        { /* change the position of the second and the last nodes */
	              itemp                          = *( *(linkList+nIntersect) + 2);
	              *( *(linkList+nIntersect) + 2) = *( *(linkList+nIntersect) + 4);
	              *( *(linkList+nIntersect) + 4) = itemp;
	        }

	    }
	 
	    /* regester the intersect tetrahedral element  */ 
	     *(intersectIndex+nIntersect) = i;

	    /* record how many points of this cut */
	    if( (nUp == 2 ) )  
	          *( *( linkList + nIntersect ) + 0 ) = 4;
	    else 
	          *( *( linkList + nIntersect ) + 0 ) = 3;

	    /* move to the next tetrahedron element */
            nIntersect++;
       }

     }
  }
  /* get the number of triangulars */
   numberOfTriangular = 0;
  for(k1=0; k1< nIntersect; k1++)
  {   if(  *(*(linkList+k1) +0)  == 3 )
         numberOfTriangular += 1;
       if( *(*(linkList+k1) +0)   == 4 )
       { /* output by divide it into two triangles */
         numberOfTriangular += 2;
       }
  }
  /* Now the number of points cutted (no redundancy counting) */
  /* *noRedundancyCutingNum  = icount; */
  /* Now the number of triangular elments */
  
  /* *numTotalTrangularElem  = numberOfTriangular */;
  /* we can just assign the value to the wmt2D5 mesh here?  */
     
  wmt2D5->nElem   = numberOfTriangular;      /* total number of the tiangule elements */
  wmt2D5->nNodes  = icount;                  /* total number of 2D Nodes  */
  wmt2D5->zConst  = (double)zConst;

  /* initialize the wmt2D5 */
  WlzInitializeWmt2D5(wmt2D5);

  /* Now established the 2D5 triangle mesh  */
  Make2D5TriangleMeshFromCuttedPlane(wmt2D5, planepoints, linkList, nIntersect, icount);

  /* transform the cutting plane points to the corresponding points in the source surface or 
     simply by get its displacements */
  /*   */
  GetTheSourceSurfaceOfTheCutPlanePoints( nIntersect, intersectIndex,   wmt3D, 
                                          planepoints, planepointsO,  linkList
                                        );

  Make2D5MeshDisplacement(wmt2D5, planepointsO); 
 

  /*  free the memory */
  AlcFree(planepoints);
  AlcFree(planepointsO);
  (void)AlcInt2Free(linkList);
  AlcFree(intersectIndex);
  AlcFree(planepointsO);
  AlcFree(planepoints);
    /*
    planepoints    = NULL;
    planepointsO   = NULL;
    linkList       = NULL;
    intersectIndex = NULL;
  */
  
  /* return nIntersect; */ /* this is the number of cut tetrahedrons */



}
#endif




/*!
* - Project:      Woolz
* - Return:       none 
* - Ingroup:      WlzTransform
* - Title:        validityOfMinAndMax
* - Date:         18 January 2002 
* - Author:       J. Rao, R. Baldock and B. Hill
* - Brief:        test the validity of the input 
* - Parameters: 
*     -# 	  nxmax:                   number of cuboids in x plane 
*     -# 	  nymax:                   number of cuboids in y plane 
*     -# 	  nzmax:                   number of cuboids in z plane 
*     -# 	  xmin:                    the starting positoin in x-direction.
*     -# 	  ymin:                    the starting positoin in y-direction.
*     -# 	  zmin:                    the starting positoin in z-direction.
*     -# 	  xmax:                    the maximum positoin in x-direction.
*     -# 	  ymax:                    the maximum positoin in y-direction.
*     -# 	  zmax:                    the maximum positoin in z-direction.
*/
#ifdef DEBUGME
void  static validityOfMinAndMax( const int   nxmax, const int   nymax, const int   nzmax, 
                           const double xmin, const double xmax, 
			   const double ymin, const double ymax,
			   const double zmin, const double zmax  )
{

  if(  (nxmax < 1) || (nymax < 1)  || (nzmax < 1)  )
  {
    printf("nxmax, nymax, nzmax must > 0 !");
    exit(1);
  }

  if( (xmin > xmax) || (ymin > ymax)  || (zmin > zmax) )
  {
    printf("xmin > xmax or ymin > ymax or zmin > zmax not allowed!");
    exit(1);
  }

}
#endif

/*!
* - Project:      Woolz
* - Return:       none 
* - Ingroup:      WlzTransform
* - Title:        ElementsGenerator
* - Date:         18 January 2002 
* - Author:       J. Rao, R. Baldock and B. Hill
* - Brief:        Generate wmt3D mesh element 
* - Parameters: 
*     -# 	  nxmax:                   number of cuboids in x plane 
*     -# 	  nymax:                   number of cuboids in y plane 
*     -# 	  nzmax:                   number of cuboids in z plane 
*     -# 	 *wmt3D:                   3D mesh transform.  output
*     -# 	 *intersectionTable:       intersect table. output.
*/
#ifdef DEBUGME
void  static ElementsGenerator(const int nxmax, const int nymax, const int nzmax,
                        WlzMeshTransform3D   *wmt3D, 
			const int *intersectionTable)
{
  int i;
  int neighbourLeft, neighbourRight, neighbourFront, neighbourBack,
      neighbourUp, neighbourDown, indexOfNextFirst;
  int ix, iy, iz;

  
  i                = 0;
  neighbourLeft    = 0;
  neighbourRight   = 0;
  neighbourFront   = 0;
  neighbourBack    = 0;
  neighbourUp      = 0;
  neighbourDown    = 0;
  indexOfNextFirst = 0;

  for(iz=1; iz <= nzmax; iz++)
  {
     for(iy=1; iy <= nymax; iy++)
     {
        for(ix=1; ix <= nxmax; ix++)
	{ 
	 /* only the cuboid intersect with the Woolz object will make meshes */
	 
	   if( *(intersectionTable+i) )
	   {
              if(  WlzTetrahedronProducerFromCube(
	 	        neighbourLeft, neighbourRight,
			neighbourBack, neighbourFront,
			neighbourDown, neighbourUp,
			nxmax, nymax,  nzmax,
			ix, iy, iz,
		        indexOfNextFirst, wmt3D->elements) != WLZ_ERR_NONE)
	      {
		  printf("err in produce tetrahedron. ");
	      }
	    
	    }
            indexOfNextFirst +=6;
	    i++; 
        }

     }

  }

}
#endif


/*!
* - Project:      Woolz
* - Return:       none 
* - Ingroup:      WlzTransform
* - Title:        GetIntersectionTable
* - Date:         18 January 2002 
* - Author:       J. Rao, R. Baldock and B. Hill
* - Brief:        Get intersection table
* - Parameters: 
*     -# 	  nxmax:                   number of cuboids in x plane 
*     -# 	  nymax:                   number of cuboids in y plane 
*     -# 	  nzmax:                   number of cuboids in z plane 
*     -# 	  xmin:                    the starting positoin in x-direction.
*     -# 	  ymin:                    the starting positoin in y-direction.
*     -# 	  zmin:                    the starting positoin in z-direction.
*     -# 	  xdelta:                  the size of the cuboid in x-direction.
*     -# 	  ydelta:                  the size of the cuboid in y-direction.
*     -# 	  zdelta:                  the size of the cuboid in z-direction.
*     -# 	 *wobjC:                   Woolz source obejct.
*     -# 	 *intersectionTable:       intersect table. output.
*/
#ifdef DEBUGME
void static GetIntersectionTable( const int nxmax, const int nymax, const int nzmax,
                             const double xmin, const double ymin, const double zmin,
                             const double xdelta, const double ydelta, const double zdelta,
			     WlzObject *wObjC,
			     int *intersectionTable 
			   )
{
  int i, j, k, it;
  int shiftx, shifty, shiftz;
  WlzObject         **wObjCS   = NULL;
  WlzObject          *wObj     = NULL;
  WlzObjectType       oType;
  WlzErrorNum	      errNum         = WLZ_ERR_NONE;

  /* allocate memory for objects pointer */
  if( (wObjCS = (WlzObject **)AlcRealloc(wObjCS,(nxmax*nymax*nzmax)*sizeof(WlzObject))) == NULL )
  {
      errNum = WLZ_ERR_MEM_ALLOC;
      printf("Failed to allocate memory for Woolz object!");
      exit(1);
  }

  /* generate a wlz cuboid object with xdelta, ydelta and zdelta as x-side,y-side and
    z-side  at the uper, left corner */

  oType = WLZ_3D_DOMAINOBJ;
  wObj  = WlzMakeCuboidObject( oType, xdelta/2., ydelta/2., zdelta/2., (int)xmin + xdelta/2., 
                               (int)ymin + ydelta/2., (int)zmin + zdelta/2., &errNum );
  it   = 0;
  for(k=0; k<nzmax; k++)
  {
      shiftz = k * (int) zdelta;
      for(j=0; j<nymax; j++)
      {
	  shifty = j * (int) ydelta;
	  for(i=0; i<nxmax; i++)
	  {
	    shiftx = i * (int) xdelta;
	       
	    /* shift first */
	     *(wObjCS+it) =  WlzShiftObject(wObj, shiftx, shifty, shiftz, &errNum);
	     
	    /* get intersections */
	    if(errNum != WLZ_ERR_NONE)
	    {
	       printf("Wrong with shift object");
	       printf("%d %d %d\n", i,j,k);
	       exit(1);
	    }

	    *(intersectionTable+it ) = WlzHasIntersection(*(wObjCS+it),wObjC,&errNum);
	    if(errNum != WLZ_ERR_NONE)
	    {
	       printf("Wrong with getting intersection of the shift object");
	       printf("%d %d %d\n", i,j,k);
	       exit(1);
	    }
	    it++;
	  }
      }
  }
  
  it = 0;

  for(k=0; k<nzmax; k++ )
  {
     for(j=0; j<nymax; j++)
     {
        for(i=0; i<nxmax; i++)
	{
	   WlzFreeObj( *( wObjCS + it ) );
	   it++;
	}
     }

  }


}
#endif

			   
/*!
* - Project:      Woolz
* - Return:       none 
* - Ingroup:      WlzGeometry
* - Title:        Read tiepoints
* - Date:         18 January 2002 
* - Author:       J. Rao, R. Baldock and B. Hill
* - Brief:        Read tiepoints from a file 
* - Parameters: 
*     -#         *fp:                   a file pointer pointing to the file contains tiepoints. 
*     -# 	**vxVec01:              an array pointer will contain tiepoints position.  output. 
*     -# 	**vxVec11:              an array pointer will contain tiepoints displacement or corresponding points
                                         depents on the  ReadDisplacement value. output.
*     -# 	  nTiePP:               number of tiepoints.       output.
*     -# 	  ReadDisplacement:     input to decide whether read displacement or corresponding position. 1 for yes
                                        0 for no.
*/
WlzErrorNum  read_WlzTiePoints( FILE *fp, int *nTiePP, 
                     WlzDVertex3 **vxVec01, 
		     WlzDVertex3 **vxVec11, const int  ReadDisplacement )
{
  WlzErrorNum errNum = WLZ_ERR_NONE;
  static char          inRecord[IN_RECORD_MAX];
  char                *rec;
  int                  vxCount       = 0, 
                       vxLimit       = 0;

  WlzDVertex3         *vx0=NULL, *vx1=NULL;
  WlzDVertex3         *vxVec0=NULL, *vxVec1=NULL;
  
  while( (errNum == WLZ_ERR_NONE) &&
	       (fgets(inRecord, IN_RECORD_MAX - 1, fp) != NULL) )
  {
    inRecord[IN_RECORD_MAX - 1] = '\0';
    rec = inRecord;

    while(*rec && isspace(*rec))
    {
      ++rec;
    }
    if(*rec && (*rec != '#'))
    {
       if(vxCount >= vxLimit)
       {
          vxLimit = (vxLimit + 1024) * 3;
          if(((vxVec0 = (WlzDVertex3 *)AlcRealloc(vxVec0,
			     vxLimit * sizeof(WlzDVertex3))) == NULL) ||
	     ((vxVec1 = (WlzDVertex3 *)AlcRealloc(vxVec1,
			     vxLimit * sizeof(WlzDVertex3))) == NULL))
          {
 	     errNum = WLZ_ERR_MEM_ALLOC;
          }
          else
          {
	     vx0 = vxVec0 + vxCount;
	     vx1 = vxVec1 + vxCount;
          }
       }
       if(errNum == WLZ_ERR_NONE)
       {
          if(sscanf(rec, "%lg %lg %lg %lg %lg %lg", &(vx0->vtX), &(vx0->vtY), &(vx0->vtZ),
		&(vx1->vtX), &(vx1->vtY), &(vx1->vtZ)) != 6)
          {
	    errNum = WLZ_ERR_READ_INCOMPLETE;
          }
          else
          { if(ReadDisplacement)
	    {
	      vx1->vtX += vx0->vtX;
	      vx1->vtY += vx0->vtY;
	      vx1->vtZ += vx0->vtZ;
	    }  
            ++vx0;
	    ++vx1;
	    ++vxCount;
          }
       }
    }
 }
 *nTiePP = vxCount;
    *vxVec01 =  vxVec0;
    *vxVec11 =  vxVec1;

 if(errNum != WLZ_ERR_NONE)
 {
   (void)fprintf(stderr, "failed to read tie points file.\n");
 }
 return errNum;

}

/*!
* - Project:      Woolz
* - Return:       none 
* - Ingroup:      WlzTransform
* - Title:        AllocateMemForElemAndNodes
* - Date:         18 January 2002 
* - Author:       J. Rao, R. Baldock and B. Hill
* - Brief:        Allocate memories for 3D mesh transform 
* - Parameters: 
*     -# 	  wmt3D:       a 3D mesh transform. 
*     -# 	  nxynumber:   number of cuboids in xy plane 
*     -# 	  nxnumber:    number of cuboids along the x-direction.
*     -# 	  xmin:        the starting positoin in x-direction.
*     -# 	  ymin:        the starting positoin in y-direction.
*     -# 	  zmin:        the starting positoin in z-direction.
*     -# 	  xdelta:      the size of the cuboid in x-direction.
*     -# 	  ydelta:      the size of the cuboid in y-direction.
*     -# 	  zdelta:      the size of the cuboid in z-direction.
*/
#ifdef DEBUGME
void static AllocateMemForElemAndNodes(WlzMeshTransform3D   *wmt3D, 
                                const int nxynumber,  
			        const int nxnumber,
			        const double xmin,
			        const double ymin,
			        const double zmin,
			        const double xdelta,
			        const double ydelta,
			        const double zdelta)
{ int i, ix, iy, iz; 
  if( (wmt3D->elements = (WlzMeshElem3D *)AlcRealloc(wmt3D->elements, wmt3D->nElem * sizeof(WlzMeshElem3D))) == NULL )
  {
      /* errNum = WLZ_ERR_MEM_ALLOC; */
      printf("Failed to allocate memory!");
      exit(1);
  }

  
  if( (wmt3D->nodes = (WlzMeshNode3D *)AlcRealloc(wmt3D->nodes, wmt3D->nNodes * sizeof(WlzMeshNode3D))) == NULL )
  {
      /* errNum = WLZ_ERR_MEM_ALLOC; */
      printf("Failed to allocate memory!");
      exit(1);
  }

  /*  initial the data to -10 which means we can use
       this information to know whether a tetrahedron
       belongs to the woolz object */
  for(i=0; i<wmt3D->nNodes; i++)
  {
     (wmt3D->elements+i)->nodes[0] = -10;
     (wmt3D->elements+i)->nodes[1] = -10;
     (wmt3D->elements+i)->nodes[2] = -10;
     (wmt3D->elements+i)->nodes[3] = -10;
  }

  /* generate the positions of the nodes */
  for(i = 0; i< wmt3D->nNodes; i++)
  {
      /* get ix, iy, iz first  */
      iz = i/nxynumber;
      iy = ( i % nxynumber )/nxnumber;
      ix = ( i % nxynumber ) % nxnumber; 
      
      (wmt3D->nodes+i)->position.vtX = xmin + (double) ix * xdelta; 
      (wmt3D->nodes+i)->position.vtY = ymin + (double) iy * ydelta; 
      (wmt3D->nodes+i)->position.vtZ = zmin + (double) iz * zdelta;
  }

}
#endif

/*!
* - Project:      Woolz
* - Return:       Woolz 3D MQ mesh transform
* - Ingroup:      WlzTransform
* - Title:        WlzTetrahedronMeshFromObj.c
* - Date:         18 January 2002 
* - Author:       J. Rao, R. Baldock and B. Hill
* - Brief:         Warp a Woolz object by a 3D mesh transform 
* - Parameters:   wObjC:       Woolz object to be warped.    input 
*     -# 	 wmt3D:       a 3D mesh transform.         input 
*     -# 	 bBoxS:       A bounding box which specified by the user generally 
*                             covering the input Woolz object.
*     -# 	 numOfElemAlonX:     number of cuboids along the x-direction.
*     -# 	 numOfElemAlonY:     number of cuboids along the y-direction.
*     -# 	 numOfElemAlonZ:     number of cuboids along the z-direction.
*     -#	 errNums:     Destination error pointer, may
*		              be null.                    in and output
*/
WlzMeshTransform3D *WlzTetrahedronMeshFromObj(WlzObject *wObjC, 
                                       const WlzDBox3 bBoxS,
				       const int numOfElemAlonX,
				       const int numOfElemAlonY,
				       const int numOfElemAlonZ,
				       WlzErrorNum *errNums
					 ){
  WlzMeshTransform3D   *wmt3D = NULL;
  int nxnumber, nynumber, nxynumber;
  int numOfCubes;
  int i,j,k,it, ix, iy, iz;
  double xdelta, ydelta, zdelta;
  int *intersectionTable = NULL;
  int shiftx, shifty, shiftz;
  WlzObject         **wObjCS   = NULL;
  WlzObject          *wObj     = NULL;
  WlzObjectType       oType;
  int neighbourLeft, neighbourRight, neighbourFront, neighbourBack,
      neighbourUp, neighbourDown, indexOfNextFirst;

 /* Test the validity of the input number */
  if(  ( numOfElemAlonX< 1) || ( numOfElemAlonY < 1)  || ( numOfElemAlonZ < 1)  )
  {
    printf("numOfElemAlonX, numOfElemAlonY, numOfElemAlonZ must > 0 !");
    exit(1);
  }

  if( (bBoxS.xMin > bBoxS.xMax) || (bBoxS.yMin > bBoxS.yMax)  || (bBoxS.zMin > bBoxS.zMax) )
  {
    printf("xmin > xmax or ymin > ymax or zmin > zmax not allowed!");
    exit(1);
  }

  /* get the interval of x, y and z-directions */
  xdelta = (bBoxS.xMax - bBoxS.xMin)/(double)numOfElemAlonX;
  ydelta = (bBoxS.yMax - bBoxS.yMin)/(double)numOfElemAlonY;
  zdelta = (bBoxS.zMax - bBoxS.zMin)/(double)numOfElemAlonZ;

  /* allocate memory for mesh transform */
  /* AllocateMemFor3DMesh(wmt3D); */
  
  if( ( wmt3D = (WlzMeshTransform3D *)AlcCalloc(1, sizeof(WlzMeshTransform3D) )  ) == NULL)
  {
    *errNums = WLZ_ERR_MEM_ALLOC;
    printf("Failed to allocate memory!");
    exit(1);
  }

  /* get the number of elements and number of Nodes */
  nxnumber      = numOfElemAlonX +1;               /* total number of nodes along x direction for one line */
  nynumber      = numOfElemAlonY +1;               /* total number of nodes along y direction for one line */
  nxynumber     =  nxnumber * nynumber;            /* total number of nodes in one xy planes */
  wmt3D->nElem  =  6 * numOfElemAlonX * numOfElemAlonY * numOfElemAlonZ;  /* total number of tetrahedron elements */
  wmt3D->nNodes =  nxynumber * ( numOfElemAlonZ + 1 );    /* total number of Nodes  */
  numOfCubes    =  numOfElemAlonX * numOfElemAlonY * numOfElemAlonZ;
  /* allocate memory for elements */
  if( (wmt3D->elements = (WlzMeshElem3D *)AlcRealloc(wmt3D->elements, 
                                          wmt3D->nElem * sizeof(WlzMeshElem3D))) == NULL )
  {
      /* errNum = WLZ_ERR_MEM_ALLOC; */
      printf("Failed to allocate memory!");
      exit(1);
  }

  /* allocate memory for nodes */
  if( (wmt3D->nodes = (WlzMeshNode3D *)AlcRealloc(wmt3D->nodes, wmt3D->nNodes * sizeof(WlzMeshNode3D))) == NULL )
  {
      /* errNum = WLZ_ERR_MEM_ALLOC; */
      printf("Failed to allocate memory!");
      exit(1);
  }


  /*  initial the data to -10 which means we can use
       this information to know whether a tetrahedron
       belongs to the woolz object */
  for(i=0; i<wmt3D->nNodes; i++)
  {
     (wmt3D->elements+i)->nodes[0] = -10;
     (wmt3D->elements+i)->nodes[1] = -10;
     (wmt3D->elements+i)->nodes[2] = -10;
     (wmt3D->elements+i)->nodes[3] = -10;
  }

  /* generate the positions of the nodes */
  for(i = 0; i< wmt3D->nNodes; i++)
  {
      /* get ix, iy, iz first  */
      iz = i/nxynumber;
      iy = ( i % nxynumber )/nxnumber;
      ix = ( i % nxynumber ) % nxnumber; 
      
      (wmt3D->nodes+i)->position.vtX = bBoxS.xMin + (double) ix * xdelta; 
      (wmt3D->nodes+i)->position.vtY = bBoxS.yMin + (double) iy * ydelta; 
      (wmt3D->nodes+i)->position.vtZ = bBoxS.zMin + (double) iz * zdelta;
  }
  /* allocate memory for intersection table */
  if( (intersectionTable = (int *)AlcRealloc(intersectionTable, numOfCubes *sizeof(int))) == NULL)
  {
      *errNums = WLZ_ERR_MEM_ALLOC;
      printf("Failed to allocate memory for intersectionTable!");
      exit(1);
  }
  /* *intersectionTable = 0; */
  /*-------- get the intersection table ----------*/
  /*
  GetIntersectionTable( numOfElemAlonX, numOfElemAlonY, numOfElemAlonZ,
                        bBoxS.xMin, bBoxS.yMin, bBoxS.zMin,
                         xdelta, ydelta, zdelta,
			 wObjC,
			 intersectionTable 
			);

  */
  /* allocate memory for objects pointer */
  if( (wObjCS = (WlzObject **)AlcRealloc(wObjCS,(numOfCubes)*sizeof(WlzObject))) == NULL )
  {
      *errNums = WLZ_ERR_MEM_ALLOC;
      printf("Failed to allocate memory for Woolz object!");
      exit(1);
  }

  /* generate a wlz cuboid object with xdelta, ydelta and zdelta as x-side,y-side and
    z-side  at the uper, left corner */

  oType = WLZ_3D_DOMAINOBJ;
  wObj  = WlzMakeCuboidObject( oType, xdelta/2., ydelta/2., zdelta/2., (int)bBoxS.xMin + xdelta/2., 
                               (int)bBoxS.yMin + ydelta/2., (int)bBoxS.zMin + zdelta/2., errNums );
  it   = 0;
  for(k=0; k<numOfElemAlonZ; k++)
  {
      shiftz = k * (int) zdelta;
      for(j=0; j<numOfElemAlonY; j++)
      {
	  shifty = j * (int) ydelta;
	  for(i=0; i<numOfElemAlonX; i++)
	  {
	    shiftx = i * (int) xdelta;
	       
	    /* shift first */
	     *(wObjCS+it) =  WlzShiftObject(wObj, shiftx, shifty, shiftz, errNums);
	     
	    /* get intersections */
	    if(*errNums != WLZ_ERR_NONE)
	    {
	       printf("Wrong with shift object");
	       printf("%d %d %d\n", i,j,k);
	       exit(1);
	    }

	    *(intersectionTable+it ) = WlzHasIntersection(*(wObjCS+it),wObjC,errNums);
	    if(*errNums != WLZ_ERR_NONE)
	    {
	       printf("Wrong with getting intersection of the shift object");
	       printf("%d %d %d\n", i,j,k);
	       exit(1);
	    }
	    it++;
	  }
      }
  }


  /* elements generator  */
  
  i                = 0;
  neighbourLeft    = 0;
  neighbourRight   = 0;
  neighbourFront   = 0;
  neighbourBack    = 0;
  neighbourUp      = 0;
  neighbourDown    = 0;
  indexOfNextFirst = 0;


  for(iz=1; iz <=  numOfElemAlonZ; iz++)
  {
     for(iy=1; iy <= numOfElemAlonY; iy++)
     {
        for(ix=1; ix <= numOfElemAlonX; ix++)
	{ 
	 /* only the cuboid intersect with the Woolz object will make meshes */
	 
	   if( *(intersectionTable+i) )
	   {
              if(  WlzTetrahedronProducerFromCube(
	 	        neighbourLeft, neighbourRight,
			neighbourBack, neighbourFront,
			neighbourDown, neighbourUp,
			numOfElemAlonX, numOfElemAlonY, numOfElemAlonZ,
			ix, iy, iz,
		        indexOfNextFirst, wmt3D->elements) != WLZ_ERR_NONE)
	      {
		  printf("err in produce tetrahedron. ");
	      }
	    
	    }
            indexOfNextFirst +=6;
	    i++; 
        }

     }

  }
  return(wmt3D);
}
/*!
* - Project:     Woolz
* - Return:      2D5 mesh transform, may be null
* - Ingroup:     WlzTransform
* - Title:       WlzMeshTransformObj_3D
* - Date:        18 January 2002 
* - Author:      J. Rao, R. Baldock and B. Hill
* - Brief        Warp a Woolz object by a 3D mesh transform 
* - Parameters:  srcObj:       Woolz object to be warped.  input 
*     -# 	 wmt3D:        a 3D mesh transform.        input 
*     -# 	 interp:       interpretation type.        input
*     -#	 dstErr:	      Destination error pointer, may
*		              be null.                    in and output
*/
WlzObject      *WlzMeshTransformObj_3D( WlzObject            *srcObj,
				         WlzMeshTransform3D   *wmt3D,
				         WlzInterpolationType  interp,
 				         WlzErrorNum          *dstErr
                                        ){
  WlzMeshTransform2D5 *wmt2D5        = NULL;    /* for plane 2D triangle mesh and 3D displacement */
  WlzDomain	  dstDom;
  WlzPlaneDomain *dstPDom   = NULL;
  WlzValues	  dstValues;
  WlzObject	 *dstObj = NULL;
  WlzErrorNum	  errNum = WLZ_ERR_NONE;
  WlzIBox3        bBox;
  
  /*  2D object used as working variables  */
  WlzObject      *dstObj2D = NULL;
  WlzValues       dummyValues;
  WlzDomain	  dummyDom;
  int             planeIdx,
                  cutPosition,
                  planeCount;
  WlzPixelV       bckgrnd;

   int        **linkList;
   int          nInterSectEsimate;
   int          numEstTrangularsElem;
   int          numEstCutingPoints;
   int         *intersectIndex = NULL;
   WlzDVertex3 *planepoints  = NULL;
   WlzDVertex3 *planepointsO = NULL;
   WlzMeshScanWSp2D5   *mSnWSp = NULL;

  
  dstDom.core      = NULL;
  dummyDom.core    = NULL;
  dummyValues.core = NULL;
  dstValues.core   = NULL;
  switch(srcObj->type)
  {
    case WLZ_EMPTY_OBJ:
      dstObj = WlzMakeEmpty(&errNum);
      break;
    case WLZ_3D_DOMAINOBJ:
      if(srcObj->domain.core == NULL)
      {
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else if( srcObj->domain.core->type == WLZ_PLANEDOMAIN_DOMAIN )
      { 
	{
        /* we need bounding box transformation here which can be 
	   done by using the mesh trasform */
         bBox = WlzMQ3DTransformBBoxI3(wmt3D, &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
           /* Create a plane domain and a destination object. */
           dstPDom = WlzMakePlaneDomain( WLZ_PLANEDOMAIN_DOMAIN,
	                                 bBox.zMin, bBox.zMax,
	                                 bBox.yMin, bBox.yMax,
	                                 bBox.xMin, bBox.xMax,
						       &errNum );
           bckgrnd = WlzGetBackground(srcObj, &errNum);

	   dstValues.vox = WlzMakeVoxelValueTb( WLZ_VOXELVALUETABLE_GREY,
	                                        dstPDom->plane1,
					        dstPDom->lastpl,
					        bckgrnd, NULL, &errNum );
	   /* Now make a 2D object  */
	   if(errNum == WLZ_ERR_NONE)
	   {
	     dummyDom.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
	                             bBox.yMin, 
				     bBox.yMax,
	                             bBox.xMin, 
				     bBox.xMax,
				     &errNum
					     );
	   }

	   if(errNum == WLZ_ERR_NONE)
	   {
	     dstObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, dummyDom,
	                        dummyValues, NULL, NULL, &errNum);
	   }
	
           /*----------  */
	   dstDom.p = dstPDom;
	   planeIdx = 0;
	   planeCount = bBox.zMax  - bBox.zMin + 1;
	   /* ------------ */
	   if(errNum == WLZ_ERR_NONE)
	   {
	     dstObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, dstDom, dstValues, NULL, NULL, &errNum);
	   }
	   /* get values */
	   
	   /* Now circle the plane one by one  */
	     cutPosition = bBox.zMin;


           /*----- For efficiency we allocate memory here -----*/
                  /* get the exact number first */

                  nInterSectEsimate    = wmt3D->nElem;
                  numEstTrangularsElem = 2 * wmt3D->nElem;
                  numEstCutingPoints   = 8 * wmt3D->nElem;
                  if(nInterSectEsimate <= 0 )
                   {
                     printf(" number = %d\n", nInterSectEsimate);
                     printf("Something is wrong in getting the number of intersect tetrahedron !");
                     exit(1);
                   }
  
                  if( ( intersectIndex = (int *)AlcRealloc( intersectIndex, 
		                       (nInterSectEsimate+1) * sizeof(int) )  ) == NULL)
                   {
                     errNum = WLZ_ERR_MEM_ALLOC;
                     printf("Failed to allocate memory for intersectIdex!");
                     exit(1);
                   }

                   /* allocate memory to store the cutting points */
                   if( ( planepoints = (WlzDVertex3 *)AlcRealloc(planepoints, 
		               4 * nInterSectEsimate * sizeof(WlzDVertex3) )  ) == NULL)
                   {
                     errNum = WLZ_ERR_MEM_ALLOC;
                     printf("Failed to allocate memory for planepoints!");
                     exit(1);
                   }
   
                   /* allocate memory to store the cutting points in the original planes */
                   if( ( planepointsO = (WlzDVertex3 *)AlcRealloc(planepointsO, 
		          4*nInterSectEsimate*sizeof(WlzDVertex3) )  ) == NULL)
                   {
                      errNum = WLZ_ERR_MEM_ALLOC;
                      printf("Failed to allocate memory for planepointsO!");
                      exit(1);
                   }
                   if( (AlcInt2Malloc(&linkList, nInterSectEsimate,5) !=  ALC_ER_NONE)  )
                   {
                     errNum = WLZ_ERR_MEM_ALLOC;
                     printf("Failed to allocate memory for link list!");
                     exit(1);
                   }

                   /* allocate memory for 2D5 mesh transform */
                   if( ( wmt2D5 = (WlzMeshTransform2D5 *)AlcRealloc(wmt2D5, 
		                sizeof(WlzMeshTransform2D5) )  ) == NULL)
                   {
                      errNum = WLZ_ERR_MEM_ALLOC;
                      printf("Failed to allocate memory for 2D5 mesh!");
                      exit(1);
                   }
   
                   /* It seems better to use a fixed temperary  memory for efficiency */
                   /* numEstCutingPoints = 4 * nInterSectEsimate; */
                   if( (wmt2D5->nodes = (WlzMeshNode2D5 *)AlcCalloc(sizeof(WlzMeshNode2D5), 
		                        numEstCutingPoints)) == NULL )
                   {
                      errNum = WLZ_ERR_MEM_ALLOC;
                      printf("Failed to allocate memory!");
                      exit(1);
                   }

                   /* allocate memory for elements and nodes */
                   /*  numEstTrangularsElem  = 5 * numEstCutingPoints; */
                   if( (wmt2D5->elements = (WlzMeshElem *)AlcCalloc(sizeof(WlzMeshElem), 
		                         numEstTrangularsElem)) == NULL )
                   {
                      errNum = WLZ_ERR_MEM_ALLOC;
                      printf("Failed to allocate memory!");
                      exit(1);
                   }
	   /*---------memory allocation finished here ----------*/

	   while((errNum == WLZ_ERR_NONE) && (planeCount-- > 0))
	   {
	      /* here the grey value should be added to disObj2D   */
             	      /* debug */
	      /* if(cutPosition == 20) */
	      { 
                errNum =  WlzMeshTransformValues3D(dstObj2D,
	                                         srcObj, 
					         wmt3D,
						 wmt2D5,
					         cutPosition,
					         interp,
                                 intersectIndex,
				 planepoints,
				 planepointsO,
				 linkList,
 				 nInterSectEsimate,
                                 numEstTrangularsElem,
                                 numEstCutingPoints,
				 mSnWSp
						 );
						 
	      }	
	      if(errNum == WLZ_ERR_NONE)
	      {
	        /*  assinDomain  */
		 *(dstObj->domain.p->domains + planeIdx)   = WlzAssignDomain(dstObj2D->domain, NULL);
                /* assigin back to the 3D dstObj one layer */
	        /* printf("planeIndex: %d\n", planeIdx);  */
		 *(dstObj->values.vox->values + planeIdx ) = WlzAssignValues(dstObj2D->values, NULL); 
	      }

	      
              ++planeIdx;
              cutPosition++;
	   }
	}
	 dstObj2D->domain.core = NULL; 
	 dstObj2D->values.core = NULL; 
	 WlzFreeObj(dstObj2D); 
      }
      break;
    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
  }

  
  if(dstErr)
  {
    *dstErr = errNum;
  }






  /*  free the memory */
  AlcFree(planepoints);
  AlcFree(planepointsO);
  (void)AlcInt2Free(linkList);
  AlcFree(intersectIndex);
  WlzMeshScanWSpFree(mSnWSp);
  AlcFree(wmt2D5->elements);
  AlcFree(wmt2D5->nodes);
  AlcFree(wmt2D5);
  
  return(dstObj);
}

/*!
* \return   2D5 mesh transform, may be null
* \ingroup  WlzTransform
* \brief     Get a  2D5 mesh from a transformed  3D mesh by
*                cuting it into a plane.

* \param    zConst:       the z = const cuting plane. 
* \param    wmt3D         mesh transform.
* \param    dstErr:	  Destination error pointer, may be null.
* \author:       J. Rao, R. Baldock and B. Hill
*/
WlzMeshTransform2D5  *Wlz2D5TransformFromCut3Dmesh(double zConst, 
                                                   WlzMeshTransform3D *wmt3D,
						   WlzErrorNum *dstErr)
						   {
   WlzMeshTransform2D5 *wmt2D5 = NULL;
   int  nInterSectEsimate, 
        numEstTrangularsElem, 
        numEstCutingPoints;
   int         *intersectIndex = NULL;
   int        **linkList;
   WlzDVertex3 *planepoints  = NULL;
   WlzDVertex3 *planepointsO = NULL;


                  nInterSectEsimate    = wmt3D->nElem;
                  numEstTrangularsElem = 2 * wmt3D->nElem;
                  numEstCutingPoints   = 8 * wmt3D->nElem;
                  if(nInterSectEsimate <= 0 )
                   {
                     printf(" number = %d\n", nInterSectEsimate);
                     printf("Something is wrong in getting the number of intersect tetrahedron !");
                     exit(1);
                   }
        
                   if( ( intersectIndex = (int *)AlcRealloc( intersectIndex, 
		                                  (nInterSectEsimate+1) * sizeof(int) )  ) == NULL)
                   {
                     *dstErr = WLZ_ERR_MEM_ALLOC;
                     printf("Failed to allocate memory for intersectIdex!");
                     exit(1);
                   }

                   /* allocate memory for 2D5 mesh transform */
                   if( ( wmt2D5 = (WlzMeshTransform2D5 *)AlcRealloc(wmt2D5, 
		                sizeof(WlzMeshTransform2D5) )  ) == NULL)
                   {
                      *dstErr = WLZ_ERR_MEM_ALLOC;
                      printf("Failed to allocate memory for 2D5 mesh!");
                      exit(1);
                   }
   
                   /* It seems better to use a fixed temperary  memory for efficiency */
                   /* numEstCutingPoints = 4 * nInterSectEsimate; */
                   if( (wmt2D5->nodes = (WlzMeshNode2D5 *)AlcCalloc(sizeof(WlzMeshNode2D5), 
		                        numEstCutingPoints)) == NULL )
                   {
                      *dstErr = WLZ_ERR_MEM_ALLOC;
                      printf("Failed to allocate memory!");
                      exit(1);
                   }

                   /* allocate memory for elements and nodes */
                   /*  numEstTrangularsElem  = 5 * numEstCutingPoints; */
                   if( (wmt2D5->elements = (WlzMeshElem *)AlcCalloc(sizeof(WlzMeshElem), 
		                         numEstTrangularsElem)) == NULL )
                   {
                      *dstErr = WLZ_ERR_MEM_ALLOC;
                      printf("Failed to allocate memory!");
                      exit(1);
                   }

                   /*-----------  */
                   /* allocate memory to store the cutting points */
                   if( ( planepoints = (WlzDVertex3 *)AlcRealloc(planepoints, 
		               4 * nInterSectEsimate * sizeof(WlzDVertex3) )  ) == NULL)
                   {
                     *dstErr = WLZ_ERR_MEM_ALLOC;
                     printf("Failed to allocate memory for planepoints!");
                     exit(1);
                   }
   
                   /* allocate memory to store the cutting points in the original planes */
                   if( ( planepointsO = (WlzDVertex3 *)AlcRealloc(planepointsO, 
		          4*nInterSectEsimate*sizeof(WlzDVertex3) )  ) == NULL)
                   {
                      *dstErr = WLZ_ERR_MEM_ALLOC;
                      printf("Failed to allocate memory for planepointsO!");
                      exit(1);
                   }
                   if( (AlcInt2Malloc(&linkList, nInterSectEsimate,5) !=  ALC_ER_NONE)  )
                   {
                     *dstErr = WLZ_ERR_MEM_ALLOC;
                     printf("Failed to allocate memory for link list!");
                     exit(1);
                   }

    WlzGet2D5TrangularMeshFrom3DMesh( wmt2D5, wmt3D, (int) zConst, dstErr,
				   intersectIndex,
				   planepoints,
				   planepointsO,
				   linkList,
                                   nInterSectEsimate,
				   numEstTrangularsElem,
                                   numEstCutingPoints
                                  );
  /* free memory */
  AlcFree(planepoints);
  AlcFree(planepointsO);
  (void)AlcInt2Free(linkList);
  AlcFree(intersectIndex);
  return(wmt2D5);
 }

