#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTstDistC_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzTstDistC.c
* \author       Bill Hill
* \date         September 2012
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
* \brief	Tests for distance evaluation.
* \ingroup	BinWlzTst
*/

#include <sys/time.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <stdio.h>
#include <Wlz.h>

/* Externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

typedef enum _WlzTstDistType
{
  WLZTST_DIST_TRUE = 0,                            /* True, analytic method. */
  WLZTST_DIST_MESH,                            /* Mesh fast marching method. */
  WLZTST_DIST_C4,                         /* 4 connected domain propagation. */
  WLZTST_DIST_C6,                         /* 6 connected domain propagation. */
  WLZTST_DIST_C8,                         /* 8 connected domain propagation. */
  WLZTST_DIST_C18,                       /* 18 connected domain propagation. */
  WLZTST_DIST_C26,                       /* 26 connected domain propagation. */
  WLZTST_DIST_OCT,                          /* Octagonal domain propagation. */
  WLZTST_DIST_CNT /* Keep last, not a type of distance method but the count. */
} WlzTstDistType;

typedef struct _WlzTstDistCRec
{
  int		idx;		/* Index of the node in the mesh. */
  int		inside;		/* Is mesh node inside the object's domain. */
  WlzDVertex3	pos;		/* Position of the node in the mesh. */
  double	dist[WLZTST_DIST_CNT]; /* Computed distances. */
} WlzTstDistCRec;

static WlzErrorNum		WlzTstDistSetCDist(
				  WlzObject *oC,
				  WlzDVertex3 seed,
				  double rad,
				  double width,
			          WlzDVertex3 centre);
static WlzObject 		*WlzTstDistCreateC(
				  int dim,
				  int segWidth,
				  int cRad,
				  int cGapDeg,
				  WlzDVertex3 *dstCentre,
				  WlzErrorNum *dstErr);
static double 			WlzTstDistComputeInCDist(
				  WlzDVertex2 p,
				  WlzDVertex2 s,
				  WlzDVertex2 c,
				  double rIn,
				  WlzDVertex2 qS);
static void 			WlzTstDistCTimerStart(
				  int doit,
				  struct timeval *t);
static double 			WlzTstDistCTimerStop(
				  int doit,
				  struct timeval *t,
				  double rep);

int		main(int argc, char *argv[])
{
  int		idC,
  		option,
		dim = 2,
		maxN = 0,
  		cRad = 100,
		timer = 0,
		usage = 0,
		smooth = 0,
		repeats = 1,
  		cGapDeg = 30,
		verbose = 0,
		segWidth = 100,
		skipNodes = 0,
		laplacianItr = 100,
		skipNodesOutside = 0;
  double	minElmSz = 8.0,
  		maxElmSz = 16.0,
		laplacianAlpha = 1.0;
  WlzDistanceType dstC[WLZTST_DIST_CNT]; /* Array with distance metrics. */
  unsigned int	dimC[WLZTST_DIST_CNT]; /* Array with bit masks for valid
                                          dimensions. */
  WlzDVertex2 	seed2;
  WlzDVertex3	centre,
  		seed;
  WlzObject	*oSeed = NULL;
  WlzObject	*oC[WLZTST_DIST_CNT];
  FILE		*fP = NULL;
  WlzTstDistCRec *dst = NULL;
  WlzCMeshP 	mesh;
  WlzGreyValueWSpace *gVWSpC[WLZTST_DIST_CNT];
  struct timeval times[3];
  double	distTimes[WLZTST_DIST_CNT];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char		*meshFile = NULL,
  		*outObjBuf = NULL,
  		*outStatsFile = NULL,
  		*outObjFileBase = NULL;
  const int	maxObjFileAdd = 16;
  char		*objStr[WLZTST_DIST_CNT];
  static char	optList[] = "23hkKstva:d:e:f:g:i:m:M:o:r:w:x:y:z:",
		defMeshFile[] = "null",
  		defOutObjFileBase[] =  "C",
		defOutStatsFile[] = "-";

  seed.vtX = 400.0;
  seed.vtY = 300.0;
  seed.vtZ = 50.0;
  mesh.v = NULL;
  meshFile = defMeshFile;
  for(idC = 0; idC < WLZTST_DIST_CNT; ++idC)
  {
    oC[idC] = NULL;
    gVWSpC[idC] = NULL;
    switch(idC)
    {
      case WLZTST_DIST_TRUE:
        dimC[idC] = 12;
	dstC[idC] = WLZ_EUCLIDEAN_DISTANCE;
        objStr[idC] = "dist";
	break;
      case WLZTST_DIST_MESH:
        dimC[idC] = 12;
	dstC[idC] = WLZ_EUCLIDEAN_DISTANCE;
        objStr[idC] = "mesh";
	break;
      case WLZTST_DIST_C4:
        dimC[idC] = 4;
	dstC[idC] = WLZ_4_DISTANCE;
	objStr[idC] = "c4";
	break;
      case WLZTST_DIST_C6:
	dimC[idC] = 8;
	dstC[idC] = WLZ_6_DISTANCE;
	objStr[idC] = "c6";
	break;
      case WLZTST_DIST_C8:
        dimC[idC] = 4;
	dstC[idC] = WLZ_8_DISTANCE;
	objStr[idC] = "c8";
	break;
      case WLZTST_DIST_C18:
	dimC[idC] = 8;
	dstC[idC] = WLZ_18_DISTANCE;
	objStr[idC] = "c18";
	break;
      case WLZTST_DIST_C26:
	dimC[idC] = 8;
	dstC[idC] = WLZ_26_DISTANCE;
	objStr[idC] = "c26";
	break;
      case WLZTST_DIST_OCT:
	dimC[idC] = 12;
	dstC[idC] = WLZ_OCTAGONAL_DISTANCE;
	objStr[idC] = "oct";
	break;
    }
    distTimes[idC] = 0.0;
  }
  outObjFileBase = defOutObjFileBase;
  outStatsFile = defOutStatsFile;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case '2':
        dim = 2;
	break;
      case '3':
        dim = 3;
	break;
      case 'k':
        skipNodesOutside = 1;
	break;
      case 'K':
        skipNodes = 1;
	break;
      case 's':
        smooth = 1;
	break;
      case 'v':
        verbose = 1;
	break;
      case 'e':
        meshFile = optarg;
	break;
      case 'f':
        outObjFileBase = optarg;
	break;
      case 'g':
        usage = (sscanf(optarg, "%d", &cGapDeg) != 1);
	break;
      case 'i':
        usage = (sscanf(optarg, "%d", &laplacianItr) != 1);
	break;
      case 'm':
        usage = (sscanf(optarg, "%lg", &minElmSz) != 1);
	break;
      case 'M':
        usage = (sscanf(optarg, "%lg", &maxElmSz) != 1);
	break;
      case 'o':
        outStatsFile = optarg;
	break;
      case 'd':
        usage = (sscanf(optarg, "%d", &cRad) != 1);
	break;
      case 'r':
        usage = (sscanf(optarg, "%d", &repeats) != 1);
	break;
      case 't':
        timer = 1;
	break;
      case 'w':
        usage = (sscanf(optarg, "%d", &segWidth) != 1);
	break;
      case 'x':
        usage = (sscanf(optarg, "%lg", &(seed.vtX)) != 1);
	break;
      case 'y':
        usage = (sscanf(optarg, "%lg", &(seed.vtY)) != 1);
	break;
      case 'z':
        usage = (sscanf(optarg, "%lg", &(seed.vtZ)) != 1);
	break;
      case 'h':   /* FALLTHROUGH */
      default:
        usage = 1;
	break;
    }
  }
  seed2.vtX = seed.vtX;
  seed2.vtY = seed.vtY;
  if(usage)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  if((errNum == WLZ_ERR_NONE) && strcmp(outObjFileBase, "null"))
  {
     int  len;
     char *dot;

     if((dot = strchr(outObjFileBase, '.')) != NULL)
     {
       *dot = '\0';
     }
     len = strlen(outObjFileBase);
     outObjBuf = (char *)AlcMalloc((len + maxObjFileAdd) * sizeof(char));
  }
  /* Create the C domain object. */
  if(errNum == WLZ_ERR_NONE)
  {
    oC[WLZTST_DIST_TRUE] = WlzAssignObject(
         WlzTstDistCreateC(dim, segWidth, cRad, cGapDeg, &centre,
	                   &errNum), NULL);
  }
  /* Compute true constrained distances from seed within C object. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzTstDistSetCDist(oC[WLZTST_DIST_TRUE], seed, cRad, segWidth,
                                centre);
  }
  /* Either create conforming mesh for C object's domain or read it from a
   * file. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(strcmp(meshFile, "null") == 0)
    {
      mesh = WlzCMeshFromObj(oC[WLZTST_DIST_TRUE], minElmSz, maxElmSz, NULL, 1,
			     &errNum);
      /* Smooth mesh. */
      if((errNum == WLZ_ERR_NONE) && (smooth != 0))
      {
	errNum = WlzCMeshLaplacianSmooth(mesh, laplacianItr, laplacianAlpha,
					 0, 1);
      }
      /* Squeeze out unused nodes. */
      if(errNum == WLZ_ERR_NONE)
      {
	WlzCMeshP mesh1;

	mesh1 = WlzCMeshCopy(mesh, 1, 0, NULL, NULL, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  (void )WlzCMeshFree(mesh);
	  mesh = mesh1;
	}
      }
    }
    else
    {
      WlzObject	*o = NULL;

      errNum = WLZ_ERR_READ_EOF;
      if((*meshFile != '\0') &&
	 ((fP = (strcmp(meshFile, "-")?
		fopen(meshFile, "r"): stdin)) != NULL) &&
	 ((o = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) != NULL))
      {
        if(((dim == 2) && (o->type != WLZ_CMESH_2D)) ||
	   ((dim == 3) && (o->type != WLZ_CMESH_3D)))
        {
	  errNum = WLZ_ERR_OBJECT_TYPE;
	}
	else
	{
	  mesh.v = o->domain.core;
	  mesh = WlzCMeshCopy(mesh, 1, 0, NULL, NULL, &errNum);
	  (void )WlzFreeObj(o);
	}
      }
      if(fP && strcmp(meshFile, "-"))
      {
	(void )fclose(fP);
	fP = NULL;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzDomain dom;
    WlzValues val;
    WlzObject *o;

    if(dim == 2)
    {
      dom.cm2 = mesh.m2;
      val.core = NULL;
      o = WlzAssignObject(
	  WlzMakeMain(WLZ_CMESH_2D, dom, val, NULL, NULL, &errNum), NULL);
      maxN = mesh.m2->res.nod.maxEnt;
    }
    else /* dim == 3 */
    {
      dom.cm3 = mesh.m3;
      val.core = NULL;
      o = WlzAssignObject(
	  WlzMakeMain(WLZ_CMESH_3D, dom, val, NULL, NULL, &errNum), NULL);
      maxN = mesh.m3->res.nod.maxEnt;
    }
    oC[WLZTST_DIST_MESH] = o;
  }
  /* Compute fast marching distances in mesh. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idR;
    WlzObject	*o = NULL;

    WlzTstDistCTimerStart(timer, times);
    for(idR = 0; idR < repeats; ++idR)
    {
      (void )WlzFreeObj(o);
      if(dim == 2)
      {
	o = WlzAssignObject(
	    WlzCMeshDistance2D(oC[WLZTST_DIST_MESH], WLZ_CMESH_2D, 1, &seed2,
	                       WLZ_INTERPOLATION_BARYCENTRIC, &errNum), NULL);
      }
      else /* dim == 3 */
      {
	o = WlzAssignObject(
	    WlzCMeshDistance3D(oC[WLZTST_DIST_MESH], WLZ_CMESH_3D, 1, &seed,
			       WLZ_INTERPOLATION_BARYCENTRIC, &errNum), NULL);
      }
      if(errNum != WLZ_ERR_NONE)
      {
        break;
      }
    }
    oC[WLZTST_DIST_MESH] = o;
    distTimes[WLZTST_DIST_MESH] = WlzTstDistCTimerStop(timer, times, repeats);
  }
  /* Create seed object. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(dim == 2)
    {
      oSeed = WlzAssignObject(
	      WlzMakeSinglePixelObject(WLZ_2D_DOMAINOBJ,
	                               seed.vtX, seed.vtY, 0,
				       &errNum), NULL);
    }
    else /* dim == 3 */
    {
      oSeed = WlzAssignObject(
	      WlzMakeSinglePixelObject(WLZ_3D_DOMAINOBJ,
	                               seed.vtX, seed.vtY, seed.vtZ,
                                       &errNum), NULL);
    }
  }
  /* Compute objects with required connectivity for the dimensionality. */
  if(errNum == WLZ_ERR_NONE)
  {
    for(idC = 0; idC < WLZTST_DIST_CNT; ++idC)
    {
      int		idR;

      if((idC != WLZTST_DIST_TRUE) && (idC != WLZTST_DIST_MESH) &&
	 ((dimC[idC] & (1 << dim)) != 0))
      {
	WlzTstDistCTimerStart(timer, times);
	for(idR = 0; idR < repeats; ++idR)
	{
	  (void )WlzFreeObj(oC[idC]);
	  oC[idC] = WlzAssignObject(
		    WlzDistanceTransform(oC[WLZTST_DIST_TRUE], oSeed,
					 dstC[idC], 0.0, 0.0, &errNum), NULL);
	  if(errNum != WLZ_ERR_NONE)
	  {
	    break;
	  }
	}
	distTimes[idC] = WlzTstDistCTimerStop(timer, times, repeats);
      }
      if(errNum != WLZ_ERR_NONE)
      {
	break;
      }
    }
  }
  /* Create distance evaluation records. */
  if(errNum == WLZ_ERR_NONE)
  {
    dst = (WlzTstDistCRec *)AlcCalloc(maxN, sizeof(WlzTstDistCRec));
  }
  /* Create object grey value workspaces. */
  if(errNum == WLZ_ERR_NONE)
  {
    for(idC = 0; idC < WLZTST_DIST_CNT; ++idC)
    {
      WlzObject *o;

      o = oC[idC];
      if((o != NULL) &&
         ((o->type == WLZ_2D_DOMAINOBJ) || (o->type == WLZ_3D_DOMAINOBJ)))
      {
        gVWSpC[idC] = WlzGreyValueMakeWSp(o, &errNum);
      }
    }
  }
  /* Fill in the distance evaluation records. */
  if(errNum == WLZ_ERR_NONE)
  {
    int 	idN;
    WlzCMeshP   mesh;
    WlzIndexedValues *idx;

    if(dim == 2)
    {
      mesh.m2 = oC[WLZTST_DIST_MESH]->domain.cm2;
    }
    else /* dim == 3 */
    {
      mesh.m3 = oC[WLZTST_DIST_MESH]->domain.cm3;
    }
    idx = oC[WLZTST_DIST_MESH]->values.x;
    for(idN = 0; idN < maxN; ++idN)
    {
      WlzTstDistCRec *d;
      WlzDVertex3 p;
      WlzCMeshNodP nod;

      d = dst + idN;
      if(dim == 2)
      {
        nod.n2 = AlcVectorItemGet(mesh.m2->res.nod.vec, idN);
        p.vtX = nod.n2->pos.vtX;
        p.vtY = nod.n2->pos.vtY;
        p.vtZ = 0.0;
        d->idx = nod.n2->idx;
      }
      else /* dim == 3 */
      {
        nod.n3 = AlcVectorItemGet(mesh.m3->res.nod.vec, idN);
        p = nod.n3->pos;
        d->idx = nod.n3->idx;
      }
      d->pos = p;
      d->dist[WLZTST_DIST_MESH] = *(double *)(WlzIndexedValueGet(idx, idN));
      if(WlzInsideDomain(oC[WLZTST_DIST_TRUE], p.vtZ, p.vtY, p.vtX, NULL))
      {
	d->inside = 1;
	WlzGreyValueGet(gVWSpC[WLZTST_DIST_TRUE], p.vtZ, p.vtY, p.vtX);
	WlzGreyValueGet(gVWSpC[WLZTST_DIST_OCT], p.vtZ, p.vtY, p.vtX);
        d->dist[WLZTST_DIST_TRUE] = gVWSpC[WLZTST_DIST_TRUE]->gVal[0].dbv;
	d->dist[WLZTST_DIST_OCT]  = gVWSpC[WLZTST_DIST_OCT]->gVal[0].inv;
	if(dim == 2)
	{
	  WlzGreyValueGet(gVWSpC[WLZTST_DIST_C4], p.vtZ, p.vtY, p.vtX);
	  WlzGreyValueGet(gVWSpC[WLZTST_DIST_C8], p.vtZ, p.vtY, p.vtX);
	  d->dist[WLZTST_DIST_C4] = gVWSpC[WLZTST_DIST_C4]->gVal[0].inv;
	  d->dist[WLZTST_DIST_C8] = gVWSpC[WLZTST_DIST_C8]->gVal[0].inv;
	}
	else /* dim == 3 */
	{
	  WlzGreyValueGet(gVWSpC[WLZTST_DIST_C6], p.vtZ, p.vtY, p.vtX);
	  WlzGreyValueGet(gVWSpC[WLZTST_DIST_C18], p.vtZ, p.vtY, p.vtX);
	  WlzGreyValueGet(gVWSpC[WLZTST_DIST_C26], p.vtZ, p.vtY, p.vtX);
	  d->dist[WLZTST_DIST_C6] = gVWSpC[WLZTST_DIST_C6]->gVal[0].inv;
	  d->dist[WLZTST_DIST_C18] = gVWSpC[WLZTST_DIST_C18]->gVal[0].inv;
	  d->dist[WLZTST_DIST_C26] = gVWSpC[WLZTST_DIST_C26]->gVal[0].inv;
	}
      }
    }
  }
  /* Write seed object. */
  if((errNum == WLZ_ERR_NONE) && strcmp(outObjFileBase, "null"))
  {
      if(oSeed != NULL)
      {
        errNum = WLZ_ERR_FILE_OPEN;
	(void )sprintf(outObjBuf, "%s-seed.wlz", outObjFileBase);
	if((fP = fopen(outObjBuf, "w")) != NULL)
	{
	  errNum = WlzWriteObj(fP, oSeed);
	  (void )fclose(fP);
	  fP = NULL;
	}
      }
  }
  /* Write C objects with true distances, mesh and other distances. */
  if((errNum == WLZ_ERR_NONE) && strcmp(outObjFileBase, "null"))
  {
    for(idC = 0; idC < WLZTST_DIST_CNT; ++idC)
    {
      if(oC[idC] != NULL)
      {
        errNum = WLZ_ERR_FILE_OPEN;
	(void )sprintf(outObjBuf, "%s-%s.wlz", outObjFileBase, objStr[idC]);
	if((fP = fopen(outObjBuf, "w")) != NULL)
	{
	  errNum = WlzWriteObj(fP, oC[idC]);
	  (void )fclose(fP);
	  fP = NULL;
	}
      }
      if(errNum != WLZ_ERR_NONE)
      {
	break;
      }
    }
  }
  /* Compute and output the distance evaluation statistics. */
  if((errNum == WLZ_ERR_NONE) && strcmp(outStatsFile, "null"))
  {
    errNum = WLZ_ERR_FILE_OPEN;
    fP = (strcmp(outStatsFile, "-"))? fopen(outStatsFile, "w"): stdout;
    if(fP)
    {
      int 	idN,
      		n = 0;
      double	min[WLZTST_DIST_CNT],
      		max[WLZTST_DIST_CNT],
		sum[WLZTST_DIST_CNT],
		ssq[WLZTST_DIST_CNT];
      const double eps = 0.0001;

      errNum = WLZ_ERR_NONE;
      if(verbose)
      {
        (void )fprintf(fP,
        "# Distance evaluation times using mesh and domain propogation\n"
        "# methods the domain propogation methods have connectivities of\n"
        "# 4, 6, 8, 18, 26 and octagonal distances. Symbols used are:\n"
        "# T   = true distance computed analyticaly\n"
        "# M   = mesh based distance\n"
        "# 4   = 4 connected domain propagation (2D only)\n"
        "# 6   = 6 connected domain propagation (3D only)\n"
        "# 8   = 8 connected domain propagation (2D only)\n"
        "# 18  = 18 connected domain propagation (3D only)\n"
        "# 26  = 26 connected domain propagation (3D only)\n"
        "# O   = octagonal domain propagation\n"
        "# M'  = (T - M) / T\n"
        "# 4'  = (T - 4) / T\n"
        "# 6'  = (T - 6) / T\n"
        "# 8'  = (T - 8) / T\n"
        "# 18' = (T - 18) / T\n"
        "# 26' = (T - 26) / T\n"
        "# O'  = (T - O) / T\n"
	"# also\n"
	"# I = node index\n"
	"# E = node inside domain\n"
	"# P = node position\n");
      }
      if(timer)
      {
	if(dim == 2)
	{
	  (void )fprintf(fP,
	  "# Distance evaluation times (M, 4, 8, O)\n%g %g %g %g\n\n",
	  distTimes[WLZTST_DIST_MESH],
	  distTimes[WLZTST_DIST_C4],
	  distTimes[WLZTST_DIST_C8],
	  distTimes[WLZTST_DIST_OCT]);
	}
	else /* dim == 3 */
	{
	  (void )fprintf(fP,
	  "# Distance evaluation times (M, 6, 18, 16, O)\n%g %g %g %g %g\n\n",
	  distTimes[WLZTST_DIST_MESH],
	  distTimes[WLZTST_DIST_C6],
	  distTimes[WLZTST_DIST_C18],
	  distTimes[WLZTST_DIST_C26],
	  distTimes[WLZTST_DIST_OCT]);
	}
      }
      if(skipNodes == 0)
      {
	if(dim == 2)
	{
	  (void )fprintf(fP,
	  "# Distance evaluation (I, E, P, T, M, 4, 8 O)\n");
	}
	else /* dim == 3 */
	{
	  (void )fprintf(fP,
	  "# Distance evaluation (I, E, P, T, M, 6, 18 26 O)\n");
	}
      }
      for(idN = 0; idN < maxN; ++idN)
      {
	WlzTstDistCRec *rec;

	rec = dst + idN;
	if((skipNodesOutside == 0) || (rec->inside != 0))
	{
	  if(rec->inside != 0)
	  {
	    double d0;
	    double d[WLZTST_DIST_CNT];

	    WlzValueSetDouble(d, 0.0, WLZTST_DIST_CNT);
	    d0 = rec->dist[WLZTST_DIST_TRUE];
	    if((d0 * d0) > eps)
	    {
              for(idC = 0; idC < WLZTST_DIST_CNT; ++idC)
	      {
	        if(idC == WLZTST_DIST_TRUE)
		{
		  d[WLZTST_DIST_TRUE] = d0;
		}
		else if(oC[idC] != NULL)
		{
		  d[idC] = (d0 - rec->dist[idC]) / d0;
		}
		if(n == 0)
		{
		  ssq[idC] = d[idC] * d[idC];
		  min[idC] = max[idC] = sum[idC] = d[idC];
		}
		else
		{
		  if(d[idC] < min[idC])
		  {
		    min[idC] = d[idC];
		  }
		  else if(d[idC] > max[idC])
		  {
		    max[idC] = d[idC];
		  }
		  sum[idC] += d[idC];
		  ssq[idC] += d[idC] * d[idC];
		}
	      }
	    }
	    ++n;
	  }
	  if(skipNodes == 0)
	  {
	    if(dim == 2)
	    {
	      (void )fprintf(fP,
			     "% 8d %d %g,%g %g %g %g %g %g\n",
			     rec->idx,
			     rec->inside,
			     rec->pos.vtX, rec->pos.vtY,
			     rec->dist[WLZTST_DIST_TRUE],
			     rec->dist[WLZTST_DIST_MESH],
			     rec->dist[WLZTST_DIST_C4],
			     rec->dist[WLZTST_DIST_C8],
			     rec->dist[WLZTST_DIST_OCT]);
	    }
	    else /* dim == 3 */
	    {
	      (void )fprintf(fP,
			     "% 8d %d %g,%g %g %g %g %g %g %g\n",
			     rec->idx,
			     rec->inside,
			     rec->pos.vtX, rec->pos.vtY,
			     rec->dist[WLZTST_DIST_TRUE],
			     rec->dist[WLZTST_DIST_MESH],
			     rec->dist[WLZTST_DIST_C6],
			     rec->dist[WLZTST_DIST_C18],
			     rec->dist[WLZTST_DIST_C26],
			     rec->dist[WLZTST_DIST_OCT]);
	    }
	  }
        }
      }
      if(skipNodes == 0)
      {
        (void )fprintf(fP, "\n");
      }
      if(dim == 2)
      {
	(void )fprintf(fP,
	"# Distance evaluation (T', M', 4', 8,' O'), where T^2 > %g\n"
	"n = %d\n",
	eps, n);
      }
      else /* dim == 3 */
      {
	(void )fprintf(fP,
	"# Distance evaluation (T', M', 6', 18,' 26', O'), where T^2 > %g\n"
	"n = %d\n",
	eps, n);
      }
      if(n > 0)
      {
	(void )fprintf(fP, "min     ");
	for(idC = 0; idC < WLZTST_DIST_CNT; ++idC)
	{
	  if((dimC[idC] & (1 << dim)) != 0)
	  {
	    (void )fprintf(fP, " % 8g", min[idC]);
	  }
	}
	(void )fprintf(fP, "\n");
	(void )fprintf(fP, "max     ");
	for(idC = 0; idC < WLZTST_DIST_CNT; ++idC)
	{
	  if((dimC[idC] & (1 << dim)) != 0)
	  {
	    (void )fprintf(fP, " % 8g", max[idC]);
	  }
	}
	(void )fprintf(fP, "\n");
	(void )fprintf(fP, "mean    ");
	for(idC = 0; idC < WLZTST_DIST_CNT; ++idC)
	{
	  if((dimC[idC] & (1 << dim)) != 0)
	  {
	    (void )fprintf(fP, " % 8g", sum[idC] / n);
	  }
	}
	(void )fprintf(fP, "\n");
	(void )fprintf(fP, "rms     ");
	for(idC = 0; idC < WLZTST_DIST_CNT; ++idC)
	{
	  if((dimC[idC] & (1 << dim)) != 0)
	  {
	    double x = 0.0;

	    if(n > 1)
	    {
	      x = sqrt(ssq[idC] / n);
	    }
	    (void )fprintf(fP, " % 8g", x);
	  }
	}
	(void )fprintf(fP, "\n");
	(void )fprintf(fP, "std dev ");
	for(idC = 0; idC < WLZTST_DIST_CNT; ++idC)
	{
	  if((dimC[idC] & (1 << dim)) != 0)
	  {
	    double x = 0.0;

	    if(n > 1)
	    {
	      x = sqrt((ssq[idC] - (sum[idC] * sum[idC] / n)) / (n - 1.0));
	    }
	    (void )fprintf(fP, " % 8g", x);
	  }
	}
	(void )fprintf(fP, "\n");
      }
      if(strcmp(outStatsFile, "-"))
      {
	(void )fclose(fP);
      }
    }
  }
  /* Free distance records, grey value workspace and all objects. */
  AlcFree(dst);
  AlcFree(outObjBuf);
  (void )WlzFreeObj(oSeed);
  for(idC = 0; idC < WLZTST_DIST_CNT; ++idC)
  {
    WlzGreyValueFreeWSp(gVWSpC[idC]);
    (void )WlzFreeObj(oC[idC]);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    const char	*errMsg,
    	 	*errStr;

    errStr = WlzStringFromErrorNum(errNum, &errMsg);
    (void )fprintf(stderr, "%s: Error - %s (%s)\n",
                   argv[0], errMsg, errStr);
  }
  if(usage != 0)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-2] [-3] [-m#] [-M#] [-e<file>] [-s] [-a#]\n"
    "       [-i#] [-g#] [-d#] [-w#] [-f<file base>] [-o<file>] [-k] [-t]\n"
    "       [-v] [-r#]\n"
    "Evaluates distance transforms within a 2D or 3D C domain.\n"
    "Options with current values in brackets are:\n"
    "  -h  Help, prints this usage message.\n"
    "  -2  2D objects (%s).\n"
    "  -3  3D objects (%s).\n"
    "  -m  Minimum mesh element size (%g).\n"
    "  -M  Maximum mesh element size (%g).\n"
    "  -e  Read mesh from file (%s).\n"
    "  -s  Use mesh smoothing (%s).\n"
    "  -a  Mesh smoothing alpha parameter (%g).\n"
    "  -i  Mesh smoothing itterations (%d)\n"
    "  -g  C gap angle (degrees) (%d).\n"
    "  -d  C radius (%d).\n"
    "  -w  C segment width (%d).\n"
    "  -x  Seed column coordinate (%g),\n"
    "  -y  Seed line coordinate (%g),\n"
    "  -z  Seed plane coordinate (%g),\n"
    "  -f  Output object file base, if \"null\" no objects are\n"
    "      output (\"%s\").\n"
    "  -o  Output distance statistics file, if \"null\" no objects are\n"
    "      output (\"%s\").\n"
    "  -k  Skip nodes outside domain from statistics (%s).\n"
    "  -K  Skip individual node distances from the statistics (%s).\n"
    "  -t  Include timing data in output statistics (%s).\n"
    "  -v  Verbose output (%s)\n"
    "  -r  Number of repeats for timing (%d).\n",
    argv[0],
    (dim == 2)? "true": "false",
    (dim == 3)? "true": "false",
    minElmSz,
    maxElmSz,
    meshFile,
    (smooth)? "true": "false",
    laplacianAlpha,
    laplacianItr,
    cGapDeg,
    cRad,
    segWidth,
    seed.vtX,
    seed.vtY,
    seed.vtZ,
    outObjFileBase,
    outStatsFile,
    (skipNodesOutside)? "true": "false",
    (skipNodes)? "true": "false",
    (timer)? "true": "false",
    (verbose)? "true": "false",
    repeats);
  }
  return(errNum);
}

/*!
* \return	New 2 or 3D C object with double values not yet set.
* \ingroup	BinWlzTst
* \brief	Creates a new 2 or 3D double valued C object.
* \param	dim			Required dimension, must be either 2
* 					or 3.
* \param	segWidth		Width of the C (and z thicknes if 3D).
* \param	cRad			Radius of the C.
* \param	cGapDeg			Gap angle of the C in degrees.
* \param	dstCentre		Destination pointer for the return of
* 					the centre.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzTstDistCreateC(int dim,
				    int segWidth, int cRad, int cGapDeg,
			            WlzDVertex3 *dstCentre, WlzErrorNum *dstErr)
{
  int		objSz;
  double	ang,
		ang0,
		ang1,
		rad0,
		rad1,
  		radSq,
  		radSq0,
		radSq1;
  WlzUByte	**dat = NULL;
  WlzObject	*oC = NULL,
  		*oR = NULL,
		*oCV = NULL;
  WlzPixelV	bgdV,
  		thrV;
  WlzIVertex2	pos;
  WlzDVertex2	rel,
		relSq;
  WlzDVertex3	cen;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bgdV.v.ubv = 0;
  bgdV.type = WLZ_GREY_UBYTE;
  thrV.v.ubv = 1;
  thrV.type = WLZ_GREY_UBYTE;
  ang0 = ALG_M_PI * cGapDeg / 360.0; /* Half angle of C gap in radians. */
  ang1 = (2 * ALG_M_PI) - ang0;
  objSz = (5 * (cRad + segWidth))/ 2;
  cen.vtX = cen.vtY = objSz / 2;
  cen.vtZ = (dim == 2)? 0.0: segWidth / 2.0;
  *dstCentre = cen;
  rad0 = cRad;
  rad1 = cRad + segWidth;
  radSq0 = rad0 * rad0;
  radSq1 = rad1 * rad1;
  if(AlcUnchar2Calloc(&dat, objSz, objSz) != ALC_ER_NONE)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    oR = WlzAssignObject(
         WlzMakeRect(0, objSz - 1, 0, objSz - 1,
                     WLZ_GREY_UBYTE, (int *)*dat, bgdV, NULL, NULL,
		     &errNum), NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(pos.vtY = 0; pos.vtY < objSz; ++pos.vtY)
    {
      rel.vtY = pos.vtY - cen.vtY;
      relSq.vtY = rel.vtY * rel.vtY;
      for(pos.vtX = 0; pos.vtX < objSz; ++pos.vtX)
      {
	rel.vtX = pos.vtX - cen.vtX;
	relSq.vtX = rel.vtX * rel.vtX;
	radSq = relSq.vtX + relSq.vtY;
	if((radSq >= radSq0) && (radSq <= radSq1))
	{
	  /* Pos is within the anulus with radSq >= radSq0 &&
	   * radSq <= radSq1. */
	  ang = atan2(rel.vtY, rel.vtX);
	  if(ang < 0.0)
	  {
	    ang += 2.0 * ALG_M_PI;
	  }
	  /* Pos is within the C part of the anulus with ang >= ang0 &&
	   * ang <= ang1. */
	  if((ang > ang0) && (ang < ang1))
	  {
	    dat[pos.vtY][pos.vtX] = 255;
	  }
	}
      }
    }
    oC = WlzAssignObject(
         WlzThreshold(oR, thrV, WLZ_THRESH_HIGH, &errNum), NULL);
  }
  if(oR)
  {
    AlcFree(dat);
    (void )WlzFreeObj(oR);
  }
  else if(dat)
  {
    Alc2Free((void **)dat);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzValues	val;
    WlzDomain	dom;
    WlzObjectType gTT;

    dom.core = NULL;
    val.core = NULL;
    gTT = WlzGreyValueTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_DOUBLE, NULL);
    if(dim == 2)
    {
      val.v = WlzNewValueTb(oC, gTT, bgdV, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
        oCV = WlzMakeMain(WLZ_2D_DOMAINOBJ, oC->domain, val, NULL, NULL,
			  &errNum);
        if(errNum != WLZ_ERR_NONE)
	{
	  WlzFreeValueTb(val.v);
	}
      }
    }
    else /* dim == 3 */
    {
      WlzObject	*oT = NULL;

      dom.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
		                 0, segWidth - 1,
				 oC->domain.i->line1, oC->domain.i->lastln,
				 oC->domain.i->kol1, oC->domain.i->lastkl,
				 &errNum);
      if(errNum ==  WLZ_ERR_NONE)
      {
	int	idP;

        for(idP = 0; idP < segWidth; ++idP)
	{
	  *(dom.p->domains + idP) = WlzAssignDomain(
	  			    WlzCopyDomain(WLZ_2D_DOMAINOBJ,
				    		  oC->domain, &errNum), NULL);
	  if(errNum != WLZ_ERR_NONE)
	  {
	    break;
	  }
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	oT = WlzAssignObject(
	     WlzMakeMain(WLZ_3D_DOMAINOBJ, dom, val, NULL, NULL,
	                 &errNum), NULL);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        val.vox = WlzNewValuesVox(oT, gTT, bgdV, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        oCV = WlzMakeMain(WLZ_3D_DOMAINOBJ, oT->domain, val, NULL, NULL,
	                  &errNum);
      }
      (void )WlzFreeObj(oT);
    }
  }
  (void )WlzFreeObj(oC);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(oCV);
}

/*!
* \return	Woolz error code.
* \ingroup	BinWlzTst
* \brief	Computes the distance of all pixels in the C object from the
* 		seed point.
* 		Given a circle of radius r, centre C and a point S outside
* 		the circle. A line from S touches the circle (as a tangent)
* 		at point Q.
* 		Let
* 		    u = S - C
* 		    v = Q - C
*		then
*		    u . v = |u| |v| \cos(\theta)
*		    \frac{|v|}{|u|} = \cos(\theta)
*		    |v|^2 = r^2
*		using maxima to solve for v_x and v_y gives
*		    v_x = \frac{r^2 u_x -/+ r u_y \sqrt{u^2 - r^2}}{u^2}
*		    v_y = \frac{r^2 u_y +/- r u_x \sqrt{u^2 - r^2}}{u^2}
*		where
*		    u^2 = u_x^2 + u_y^2
*		from which we get
*		    Q_x = C_x + v_x
*		    Q_y = C_y + v_y
* 		For a 3D prism formed by planar C domains then the
* 		distance is just the norm of the distance in the plane
* 		computed as above with the distance in depth.
* \param	oC			C object with double values to be set.
* \param	s			Seed position, must be withing the C
* 					object.
* \param	r			C radius.
* \param	w			Only for 3D, the c width (ie number
* 					of planes).
* \param	c			C centre.
*/
static WlzErrorNum WlzTstDistSetCDist(WlzObject *oC, WlzDVertex3 s,
			              double r, double w,
				      WlzDVertex3 c)
{
  double	t,
  		rSq,
  		uSqLn;
  int		pln,
  		pln0,
  		pln1;
  double	*gP;
  WlzDVertex2	u,
  		uSq,
  		q0,
		q1,
		qS;
  WlzObject	*oC2;
  WlzIntervalWSpace iWSp;
  WlzGreyWSpace	gWSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Compute Qs: point on inner circle of anulus which has a tangent that
   * passes through the seed point S. */
  rSq = r * r;
  WLZ_VTX_2_SUB(u, s, c);
  uSq.vtX = u.vtX * u.vtX;
  uSq.vtY = u.vtY * u.vtY;
  uSqLn = uSq.vtX + uSq.vtY;
  t = r * sqrt(uSqLn - rSq);
  q0.vtX = c.vtX + (rSq * u.vtX - u.vtY * t) / uSqLn;
  q1.vtX = c.vtX + (rSq * u.vtX + u.vtY * t) / uSqLn;
  q0.vtY = c.vtY + (rSq * u.vtY + u.vtX * t) / uSqLn;
  q1.vtY = c.vtY + (rSq * u.vtY - u.vtX * t) / uSqLn;
  /* Need to choose one of these points for Qs. This is a hack which only
   * works because I know where the seed point is! */
  qS = (q1.vtY > q0.vtY)? q1: q0;
  if(oC->type == WLZ_2D_DOMAINOBJ)
  {
    pln0 = pln1 = 0;
  }
  else /* oC->type == WLZ_3D_DOMAINOBJ */
  {
    pln0 = oC->domain.p->plane1;
    pln1 = oC->domain.p->lastpl;
  }
  for(pln = pln0; pln <= pln1; ++pln)
  {
    int		p;
    double	dzSq = 0.0;
    WlzDomain	*domains = NULL;
    WlzValues	*values = NULL;

    if(oC->type == WLZ_2D_DOMAINOBJ)
    {
      oC2 = oC;
    }
    else
    {
      p = pln - pln0;
      domains = oC->domain.p->domains;
      values = oC->values.vox->values;
      oC2 = WlzMakeMain(WLZ_2D_DOMAINOBJ, domains[p], values[p],
                        NULL, NULL, &errNum);
      dzSq = (pln - s.vtZ) * (pln - s.vtZ);
    }
    /* For each pixel in each plane compute and set it's distance from S. */
    if((oC->type == WLZ_2D_DOMAINOBJ) || (domains[p].core != NULL))
    {
      WlzDVertex2 c2,
      		  p2,
      		  s2;

      c2.vtX = c.vtX;
      c2.vtY = c.vtY;
      s2.vtX = s.vtX;
      s2.vtY = s.vtY;
      if(errNum == WLZ_ERR_NONE)
      {
	errNum = WlzInitGreyScan(oC2, &iWSp, &gWSp);
      }
      while((errNum == WLZ_ERR_NONE) &&
	    (WlzNextGreyInterval(&iWSp) == WLZ_ERR_NONE))
      {
	switch(gWSp.pixeltype)
	{
	  case WLZ_GREY_DOUBLE:
	    p2.vtY = iWSp.linpos;
	    gP = gWSp.u_grintptr.dbp;
	    for(p2.vtX = iWSp.lftpos; p2.vtX <= iWSp.rgtpos; ++p2.vtX)
	    {
	      double	d,
			dxy;

	      dxy = WlzTstDistComputeInCDist(p2, s2, c2, r, qS);
	      d = (oC->type == WLZ_2D_DOMAINOBJ)?
		  dxy: sqrt((dxy * dxy) + dzSq);
	      *gP++ = d;
	    }
	    break;
	  default:
	    errNum = WLZ_ERR_GREY_TYPE;
	    break;
	}
      }
    }
    if(oC->type != WLZ_2D_DOMAINOBJ)
    {
      (void )WlzFreeObj(oC2);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      break;
    }
  }
  return(errNum);
}

static double WlzTstDistComputeInCDist(WlzDVertex2 p, WlzDVertex2 s,
                                       WlzDVertex2 c, double r, WlzDVertex2 qS)
{
  double	d = 0.0,
  		cSQS,
		mSQS,
		rSq,
  		t,
		uSqLn;
  WlzDVertex2	u0,
		u1,
		uSq,
		q0,
		q1,
		qP;

  rSq = r * r;
  /* Compute eqn of line through s and qS with y = mSQS * x + cSQS.
   * This is also a hack since we know that this line is neither horizontal
   * or vertical. */
  mSQS = (s.vtY - qS.vtY)/(s.vtX - qS.vtX);
  cSQS = s.vtY - mSQS * s.vtX;
  /* if p is below this line or above it and to the right (but  below the
   * centre) then the distance is just the Euclidean distance p,s. */
  t = mSQS * p.vtX + cSQS;
  if(((p.vtY + DBL_EPSILON > t) || (p.vtX > qS.vtX)) && (p.vtY > c.vtY))
  {
    WLZ_VTX_2_SUB(u0, p, s);
    d = WLZ_VTX_2_LENGTH(u0);
  }
  else
  {
    /* Compute Qp: point on inner circle of anulus which has a tangent that
     * passes through the test point P. */
    WLZ_VTX_2_SUB(u0, p, c);
    uSq.vtX = u0.vtX * u0.vtX;
    uSq.vtY = u0.vtY * u0.vtY;
    uSqLn = uSq.vtX + uSq.vtY;
    t = r * sqrt(uSqLn - rSq);
    q0.vtX = c.vtX + (rSq * u0.vtX - u0.vtY * t) / uSqLn;
    q0.vtY = c.vtY + (rSq * u0.vtY + u0.vtX * t) / uSqLn;
    q1.vtX = c.vtX + (rSq * u0.vtX + u0.vtY * t) / uSqLn;
    q1.vtY = c.vtY + (rSq * u0.vtY - u0.vtX * t) / uSqLn;
    /* Need to choose one of these points for Qp. This is a hack which only
     * works because I know where the seed point is! */
    WLZ_VTX_2_SUB(u0, qS, c);
    WLZ_VTX_2_SUB(u1, qP, c);
    qP = (WlzGeomCmpAngle(u0, u1) > 0)? q0: q1;
    /* Can now compute constrained distance dPS = dPQP + dQPQS + dSQS.*/
    /* dPQP */
    WLZ_VTX_2_SUB(u0, p, qP);
    d = WLZ_VTX_2_LENGTH(u0);
    /* dSQS */
    WLZ_VTX_2_SUB(u0, s, qS);
    d += WLZ_VTX_2_LENGTH(u0);
    /* dQPQS */
    d += WlzGeomArcLength2D(qS, qP, c);
  }
  return(d);
}

static void WlzTstDistCTimerStart(int doit, struct timeval *t)
{
  if(doit)
  {
    gettimeofday(t + 0, NULL);
  }
}

static double WlzTstDistCTimerStop(int doit, struct timeval *t, double rep)
{
  double s = 0.0;

  if(doit)
  {
    gettimeofday(t + 1, NULL);
    ALC_TIMERSUB(t + 1, t + 0, t + 2);
    s = t[2].tv_sec + (0.000001 * t[2].tv_usec);
    if(rep > 0)
    {
      s /= rep;
    }
  }
  return(s);
}
