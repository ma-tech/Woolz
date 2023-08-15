#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCut2DTiles_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/WlzCut2DTiles.c
* \author       Bill Hill
* \date         July 2023
* \version      $Id$
* \par
* Address:
*               Heriot Watt
*               Edinburgh, Scotland,
*               UK EH14 4AS
* \par
* Copyright (C), [2016],
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
* \brief	Cuts 2D rectangular regions of interest from either a 3D
* 		domain object.
* \ingroup	BinWlzApp
*
* \par Binary
* \ref wlzcut2dtiles "WlzCut2DTiles"
*/

#include <stdio.h>
#include <float.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <Wlz.h>
#include <WlzExtFF.h>
#include <cJSON.h>

/*!
\ingroup BinWlz
\defgroup wlzcut2dtiles WlzCut2DTiles
\par Name
WlzCut2DTiles
\par Synopsis
\verbatim
WlzCut2DTiles [-h] [-v] [-a(n|r)]
              [-b<background tile value>] [-c<record file>]
              [-e<columns>,<lines>]
              [-E<input file format>] [-f<file base>]
              [-F<output tile format>]
              [-i]<in domain value>] [-M<max tile images>]
              [-N<max tile images per plane>] [-p(o|p|r)]
              [-r<pitch>,<yaw>,<roll>] [-R<rescale flags>
              [-s<columns>,<lines>] [-t(s|r)] [-Y] [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Prints verbose progress messages to stderr.</td>
  </tr>
  <tr> 
    <td><b>-a</b></td>
    <td>Type of within plane rotation (about the plane centre) - none
        (n),  random (r), default none.</td>
  </tr>
  <tr> 
    <td><b>-b</b></td>
    <td>Background tile value, default 0.</td>
  </tr>
  <tr> 
    <td><b>-c</b></td>
    <td>Cut all tile images using the given record file, default false.</td>
  </tr>
  <tr> 
    <td><b>-d</b></td>
    <td>Distance increment for parallel and orthogonal planes, must
        be > 0, default 1.</td>
  </tr>
  <tr> 
    <td><b>-e</b></td>
    <td>Tile overlap for raster tiling, default 0,0.</td>
  </tr>
  <tr> 
    <td><b>-E</b></td>
    <td>Input 3D image format, default Woolz (wlz).</td>
  </tr>
  <tr> 
    <td><b>-f</b></td>
    <td>Output base file name, default "".</td>
  </tr>
  <tr> 
    <td><b>-F</b></td>
    <td>Output image tile format, default Woolz (wlz).</td>
  </tr>
  <tr> 
    <td><b>-L</b></td>
    <td>Use linear interpolation instead of nearest neighbour when
        cutting planes, default false.</td>
  </tr>
  <tr> 
    <td><b>-i</b></td>
    <td>Within domain tile value, range 1-255, default 255.</td>
  </tr>
  <tr> 
    <td><b>-M</b></td>
    <td>Maximum number of tile images, default -1.</td>
  </tr>
  <tr> 
    <td><b>-N</b></td>
    <td>Maximum number of tile images for each plane cut (the same
        plane may possibly be cut more than once), default -1.</td>
  </tr>
  <tr> 
    <td><b>-p</b></td>
    <td>Type of planes to cut - orthogonal (o), parallel (p) or
      random (r). The orthogonal and parallel planes are with respect
      to the referece plane, default orthogonal.</td>
  </tr>
  <tr>
    <td><b>-r</b></td>
    <td>Reference plane angles (given using degrees), defined by pitch,
        yaw and roll (also known as phi, theta and zeta), default 0,0,0.</td>
  </tr>
  <tr>
    <td><b>-R</b></td>
    <td>Voxel rescaling flags:
        bit 1 set - use voxel-size rescaling,
        bit 2 set - enable global scaling, default 0</td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td>Size of tile images to cut, default 256,256.</td>
  </tr>
  <tr>
    <td><b>-t</b></td>
    <td>Within plane section tiling - raster (s), random (r), default
        raster.</td>
  </tr>
  <tr>
    <td><b>-Y</b></td>
    <td>Include empty tiles (entirely outside the input domain or
    of background) in output, default false.</td>
  </tr>
</table>
\par Description
Cuts 2D region of interest tiles from a 3D domain object.
\par Examples
\verbatim

  ./WlzCut2DTiles -f tiles in.wlz
\endverbatim
Cuts 2D tiles from the input file in.wlz writing the tile cut record
to tiles_cut.jsn and tile images to tiles_000000.wlz, etc....

\verbatim
  ./WlzCut2DTiles -f tiles -d prev_cut.jsn in.wlz
\endverbatim
Cuts 2D tiles from the input file in.wlz as with the previous
example except that the cutting parameters are read from the given file
prev_cut.jsn

\verbatim
  ./WlzCut2DTiles -s 64,128 -p r -a r -f tiles -M 1000 -F jpg in.wlz
\endverbatim
Cuts 2D 64x128 tiles from the input file in.wlz using random planes
and random in plane rotations, writing the tile cut record to
tiles_cut.jsn and tile images to tiles_000000.jpg, ...,
tiles_000999.jpg.
\par File
\ref WlzCut2DTiles.c "WlzCut2DTiles.c"
\par See Also
\ref BinWlzApp "WlzIntro(1)"
\ref WlzGetSubSectionFromObject "WlzGetSubSectionFromObject(3)"
\ref WlzEffReadObj "WlzEffReadObj(3)"
\ref WlzEffWriteObj "WlzEffWriteObj(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/* Externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

typedef enum _WlzC2DKeys
{
  WLZ_C2DKEY_FALSE	= 0,
  WLZ_C2DKEY_TRUE	= 1,
  WLZ_C2DKEY_NONE,
  WLZ_C2DKEY_ORTHOGONAL,
  WLZ_C2DKEY_PARALLEL,
  WLZ_C2DKEY_RANDOM,
  WLZ_C2DKEY_RASTER,
  WLZ_C2DKEY_COUNT
} WlzC2DKeys;

typedef enum _WlzC2DAngles 
{
  WLZ_C2DANG_PITCH = 0,
  WLZ_C2DANG_PHI   = 0,
  WLZ_C2DANG_YAW   = 1,
  WLZ_C2DANG_THETA = 1,
  WLZ_C2DANG_ROLL  = 2,
  WLZ_C2DANG_ZETA  = 2
} WlzC2DAngles;


static const char * const keyStrings[] =
  		{
  		  /* Must correspond to WlzC2DKeys */
		  "false",
		  "true",
		  "none",
		  "orthogonal",
		  "parallel",
		  "random",
		  "raster"
		};

static const char *prog = "WlzCut2DTiles";

static int 	ParseStrOption(const char *par);
static double	NormaliseAngle(double a, double max);

int		main(int argc, char *argv[])
{
  int		arg,
  		option,
		bgdVal = 0,
		distInc = 1,
		inDomVal = 255,
		includeEmpty = 0,
		rotType = WLZ_C2DKEY_NONE,
		plnType = WLZ_C2DKEY_ORTHOGONAL,
		plnTilType = WLZ_C2DKEY_RASTER,
		minFgSz = 8,
		maxTiles = -1,
		maxTilesPP = -1,
		maxRecTiles = 0,
		voxRescale = 0,
  		ok = 1,
  		usage = 0,
		verbose = 0;
  char		*inRecStr = NULL,
  		*baseFile = "",
	 	*inRecFile = NULL,
		*inObjFile = "-";
  char 		*filenameBuf = NULL;
  cJSON		*jRec = NULL,
  		*jRecTiles = NULL;
  FILE		*outRecFP = NULL;
  WlzObject	*inObj = NULL;
  WlzIBox3	inBBox;
  WlzEffFormat	inObjFileFmt = WLZEFF_FORMAT_WLZ,
  		tileFileFmt = WLZEFF_FORMAT_WLZ;
  WlzInterpolationType interp = WLZ_INTERPOLATION_NEAREST;
  WlzIVertex2	tileSz = {256, 256},
  		tileOverlap = {0, 0};
  WlzDVertex3	refPlnAngles = {0.0, 0.0, 0.0};
  const WlzDVertex3 fixed = {0.0, 0.0, 0.0},
		up = {0.0, 0.0, -1.0};
  WlzThreeDViewStruct *vWSp = NULL;
  const char	*errMsg;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char   optList[] = "hLRvYa:b:c:d:e:E:f:F:i:M:N:p:r:s:t:";

  opterr = 0;
  AlgRandSeed(0);
  while((usage == 0) && (errNum == WLZ_ERR_NONE) &&
        ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'v':
        verbose = 1;
	break;
      case 'L':
        interp = WLZ_INTERPOLATION_LINEAR;
	break;
      case 'b':
        if((sscanf(optarg, "%d", &(bgdVal)) != 1) ||
           (bgdVal < 0) || (bgdVal > 255))
      case 'R':
        voxRescale = 1;
	break;
      case 'a':
	arg = ParseStrOption(optarg);
        switch(arg)
	{
	  case WLZ_C2DKEY_NONE:   /* FALLTHROUGH */
	  case WLZ_C2DKEY_RANDOM:
	    rotType = arg;
	    break;
	  default:
	    usage = 1;
	    break;
	}
	break;
      case 'c':
        inRecFile = optarg;
        break;
      case 'd':
	if((sscanf(optarg, "%d", &(distInc)) != 1) || (distInc < 1))
        {
	  usage = 1;
	}
	break;
      case 'e':
	if(sscanf(optarg, "%d,%d", &(tileOverlap.vtX), &(tileOverlap.vtY)) != 2)
        {
	  usage = 1;
	}
	break;
      case 'E':
        if((inObjFileFmt = WlzEffStringExtToFormat(optarg)) == 0)
        {
	  usage = 1;
	}
	break;
      case 'f':
	baseFile = optarg;
	break;
      case 'F':
        if((tileFileFmt = WlzEffStringExtToFormat(optarg)) == 0)
        {
	  usage = 1;
	}
	break;
      case 'i':
	if((sscanf(optarg, "%d", &(inDomVal)) != 1) ||
	   (inDomVal < 1) || (inDomVal > 255))
        {
	  usage = 1;
	}
	break;
      case 'M':
        if(sscanf(optarg, "%d", &maxTiles) != 1)
	{
	  usage = 1;
	}
	break;
      case 'N':
        if(sscanf(optarg, "%d", &maxTilesPP) != 1)
	{
	  usage = 1;
	}
	break;
      case 'p':
        arg = ParseStrOption(optarg);
	switch(arg)
	{
	  case WLZ_C2DKEY_ORTHOGONAL: /* FALLTHROUGH */
	  case WLZ_C2DKEY_PARALLEL:   /* FALLTHROUGH */
	  case WLZ_C2DKEY_RANDOM:
	    plnType = arg;
	    break;
	  default:
            usage = 1;
            break;
	}
        break;
      case 'r':
        if(sscanf(optarg, "%lg,%lg,%lg",
	          &(refPlnAngles.vtX), &(refPlnAngles.vtY),
		  &(refPlnAngles.vtZ)) != 3)
        {
	  usage = 1;
	}
	/* Commandline uses degrees, internally we use radians. */
	WLZ_VTX_3_SCALE(refPlnAngles, refPlnAngles, ALG_M_PI / 180.0);
	break;
      case 's':
        if(sscanf(optarg, "%d,%d", &(tileSz.vtX), &(tileSz.vtY)) != 2)
        {
	  usage = 1;
	}
	break;
      case 't':
	arg = ParseStrOption(optarg);
	switch(arg)
	{
	  case WLZ_C2DKEY_RASTER:
	  case WLZ_C2DKEY_RANDOM:
	    plnTilType = arg;
	    break;
	  default:
            usage = 1;
	    break;
	}
	break;
      case 'Y':
        includeEmpty = 1;
	break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  if((usage == 0) && (errNum != WLZ_ERR_NONE))
  {
    ok = 0;
    (void )WlzStringFromErrorNum(errNum, &errMsg);
    (void )fprintf(stderr,
		   "%s Failed to parse command line, %s.\n",
		   argv[0],
		   errMsg);
  }
  if(ok && (usage == 0))
  {
    if((inObjFile == NULL) || (*inObjFile == '\0'))
    {
      usage = 1;
    }
    else if(optind < argc)
    {
      if((optind + 1) != argc)
      {
        usage = 1;
      }
      else
      {
        inObjFile = *(argv + optind);
      }
    }
  }
  ok = (ok != 0) && (usage == 0);
  if(ok)
  {
    /* Allocate a filename buffer in which to build filenames. */
    if((filenameBuf = (char *)AlcCalloc(strlen(baseFile) + 64,
        sizeof(char))) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to allocate filename buffer.\n", argv[0]);
    }
  }
  if(ok)
  {
    /* Open the output record file. */
    (void )sprintf(filenameBuf, "%s_cut.jsn", baseFile);
    if((outRecFP = fopen(filenameBuf, "w")) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                    "%s: failed to open output record file.\n", argv[0]);
    }
    if(ok)
    {
      int	i = 1;

      if(fprintf(outRecFP, "{\n\"commandline\": \"%s", argv[0]) < 1)
      {
        ok = 0;
      }
      while(ok && (i < argc))
      {
        ok = (fprintf(outRecFP, " %s", argv[i]) > 0);
	++i;
      }
      if(ok)
      {
        ok = (fprintf(outRecFP, "\",\n\"tiles\": [\n") > 0);
      }
      if(!ok)
      {
        (void )fprintf(stderr,
                       "%s: failed initial write to record file %s\n.",
		       argv[0], filenameBuf);
      }
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(verbose)
    {
      (void )fprintf(stderr, "%s: Reading input image from file %s\n",
	  argv[0], (strcmp(inObjFile, "-") == 0)? "stdin": inObjFile);
    }
    if(strcmp(inObjFile, "-") == 0)
    {
      if((inObj = WlzAssignObject(WlzReadObj(stdin, &errNum), NULL)) == NULL)
      {
        ok = 0;
      }
    }
    else
    {
      if((inObj = WlzAssignObject(
		  WlzEffReadObj(NULL, inObjFile, inObjFileFmt, 0, 0, 0,
			        &errNum), NULL)) == NULL)
      {
        ok = 0;
      }
    }
    if(ok == 0)
    {
      (void )fprintf(stderr,
                     "%s: failed to read object from file %s\n",
                     argv[0], inObjFile);
    }
  }
  if(ok)
  {
    /* Check the input object it needs to be a 3D domain object. */
    inBBox = WlzBoundingBox3I(inObj, &errNum);
    if(inObj->type != WLZ_3D_DOMAINOBJ)
    {
      errNum = WLZ_ERR_OBJECT_TYPE;
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Object read from %s is not a recognised 3D image\n",
		     argv[0],
		     (strcmp(inObjFile, "-") == 0)? "stdin": inObjFile);
    }
  }
  if(ok && inRecFile)
  {
    long	len = 0;
    FILE	*fP = NULL;
    const char  *errPtr = NULL;

    if(verbose)
    {
      (void )fprintf(stderr, "%s: Reading record file  %s\n",
          argv[0], (strcmp(inObjFile, "-") == 0)? "stdin": inObjFile);
    }
    errNum = WLZ_ERR_NONE;
    if(((fP = fopen(inRecFile, "r")) == NULL) ||
       (fseek(fP, 0, SEEK_END) != 0) ||
       ((len = ftell(fP)) <= 0) ||
       (fseek(fP, 0, SEEK_SET) != 0))
    {
      errNum = WLZ_ERR_READ_EOF;
    }
    else if((inRecStr = (char *)AlcCalloc(len + 1, sizeof(char))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else if(fread(inRecStr, sizeof(char), len, fP) < len)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    else
    {
      inRecStr[len] = '\0';
      if((jRec = cJSON_Parse(inRecStr)) == NULL)
      {
	errNum = WLZ_ERR_FILE_FORMAT;
	errPtr = cJSON_GetErrorPtr();
      }
      else
      {
	const cJSON *cmdln;

	/* Check that the JSON file has a matching command line and has
	 * tile data. */
        cmdln = cJSON_GetObjectItemCaseSensitive(jRec, "commandline");
	jRecTiles = cJSON_GetObjectItemCaseSensitive(jRec, "tiles");
	if((cmdln == NULL) || (jRecTiles == NULL) ||
	   (cJSON_IsString(cmdln) == 0) || (cJSON_IsArray(jRecTiles) == 0) ||
	   (cmdln->valuestring == NULL) ||
	   (strstr(cmdln->valuestring, prog) == NULL) ||
	   ((maxRecTiles = cJSON_GetArraySize(jRecTiles)) < 1))
        {
	  errNum = WLZ_ERR_FILE_FORMAT;
	}
	else
	{
	  maxRecTiles = cJSON_GetArraySize(jRecTiles);
	}
      }
    }
    if(fP)
    {
      (void )fclose(fP);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
          "%s: failed to parse cut record file %s (%s)\n",
	  argv[0], inRecFile, errMsg);
      if(errPtr)
      {
        (void )fprintf(stderr, "\t\tpossible problem before %s\n",
	    errPtr);
      }
    }
  }
  if(ok)
  {
    vWSp = WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      (void )fprintf(stderr,
                     "%s: failed to allocate view structure\n",
		     argv[0]);
      ok = 0;
    }
  }
  if(ok)
  {
    int		plCnt = 0,           /* Number of planes since inPlTlCnt == 0 */
    		tlCnt = 0,         /* Number of tiles (including empty tiles) */
		ortCnt = 0,          /* Number of orthogonal directions done. */
		inPlTlCnt = 0,           /* Number of tiles in current plane. */
		allTlCnt = 0,                        /* Total number of tiles */
		inRecTlCnt = 0,  /* Number of tile records read from cut file */
		emptyPl = 0;
    double	tlRotAngle = {0.0};
    double	plAngles[3] = {0.0, 0.0, 0.0};
    WlzIBox2	rotPlBBox,
    		tlBox;
    WlzObject	*plObj = NULL,
		*tlObj = NULL;
    

    while(errNum == WLZ_ERR_NONE)
    {
      int	tlIndex,
      		emptyTl = 0;
      cJSON	*jRecTl = NULL,
      		*jRecTlPl = NULL;
      WlzObject *clpRotPlObj = NULL;

      if(jRecTiles)
      {
        jRecTl = cJSON_GetArrayItem(jRecTiles, inRecTlCnt);
	if(((jRecTl = cJSON_GetArrayItem(
	                   jRecTiles, inRecTlCnt)) == NULL) ||
	   (cJSON_IsObject(jRecTl) == 0) ||
	   ((jRecTlPl  = cJSON_GetObjectItemCaseSensitive(
	                      jRecTl, "plane")) == NULL)  ||
	   (cJSON_IsObject(jRecTlPl) == 0))
	{
	  errNum = WLZ_ERR_FILE_FORMAT;
	}
	else
	{
	  ++inRecTlCnt;
	}
      }
      if((errNum == WLZ_ERR_NONE) && (inPlTlCnt == 0))
      {
	/* Get a new plane. */
	/* Free an existing plane. */
	(void )WlzFreeObj(plObj); plObj = NULL;
	/* Set section parameters. */
	if(jRecTlPl)
	{
	  cJSON	*jPhi = NULL,
		*jTheta = NULL,
		*jZeta = NULL;

	  if(((jPhi = cJSON_GetObjectItemCaseSensitive(
	          jRecTlPl, "phi")) == NULL) ||
	     (cJSON_IsNumber(jPhi) == 0) ||
	     ((jTheta = cJSON_GetObjectItemCaseSensitive(
	          jRecTlPl, "theta")) == NULL) ||
	     (cJSON_IsNumber(jTheta) == 0) ||
	     ((jZeta = cJSON_GetObjectItemCaseSensitive(
	          jRecTlPl, "zeta")) == NULL) ||
	     (cJSON_IsNumber(jZeta) == 0))
	  {
	    errNum = WLZ_ERR_FILE_FORMAT;
	  }
	  else
	  {
	    plAngles[WLZ_C2DANG_PHI] = jPhi->valuedouble;
	    plAngles[WLZ_C2DANG_THETA] = jTheta->valuedouble;
	    plAngles[WLZ_C2DANG_ZETA] = jZeta->valuedouble;
	  }
	}
	else
	{
	  switch(plnType)
	  {
	    case WLZ_C2DKEY_ORTHOGONAL:
	      switch(ortCnt)
	      {
	        case 0:
		  plAngles[WLZ_C2DANG_PITCH] = refPlnAngles.vtX;
                  plAngles[WLZ_C2DANG_YAW] = refPlnAngles.vtY;
                  plAngles[WLZ_C2DANG_ROLL] = refPlnAngles.vtZ;
		  break;
		case 1:
		  plAngles[WLZ_C2DANG_PITCH] = NormaliseAngle(
		      refPlnAngles.vtX + ALG_M_PI / 2.0, 2.0 * ALG_M_PI);
		  break;
		case 2:
		  plAngles[WLZ_C2DANG_YAW] = NormaliseAngle(
		      refPlnAngles.vtY + ALG_M_PI / 2.0, 2.0 * ALG_M_PI);
		  break;
		default:
		  break;
	      }
	      break;
	    case WLZ_C2DKEY_PARALLEL:
	      plAngles[WLZ_C2DANG_PITCH] = refPlnAngles.vtX;
	      plAngles[WLZ_C2DANG_YAW] = refPlnAngles.vtY;
	      plAngles[WLZ_C2DANG_ROLL] = refPlnAngles.vtZ;
	      break;
	    case WLZ_C2DKEY_RANDOM:
	      plAngles[WLZ_C2DANG_PITCH] = 2.0 * ALG_M_PI * AlgRandUniform();
	      plAngles[WLZ_C2DANG_YAW]   = ALG_M_PI * AlgRandUniform();
	      plAngles[WLZ_C2DANG_ROLL] =  2.0 * ALG_M_PI * AlgRandUniform();
	      break;
	    default:
	      errNum = WLZ_ERR_PARAM_DATA;
	      break;
	  }
	}
	/* Setup the view structure. */
	if(errNum == WLZ_ERR_NONE)
	{
	  if(inPlTlCnt == 0)
	  {
	    vWSp->phi =   plAngles[WLZ_C2DANG_PHI];
	    vWSp->theta = plAngles[WLZ_C2DANG_THETA];
	    vWSp->zeta =  plAngles[WLZ_C2DANG_ZETA];
	    vWSp->up = up;
	    vWSp->fixed = fixed;
	    /* I don't understand why this fails with dist = 0. */
	    vWSp->dist = (inBBox.xMax + inBBox.xMin +
	        inBBox.yMax + inBBox.yMin +
		inBBox.zMax + inBBox.zMin) / 6;
	    vWSp->scale = 1.0;
	    vWSp->view_mode = WLZ_ZETA_MODE;
	    vWSp->voxelRescaleFlg = voxRescale;
	  }
	}
	/* Get plane index if using cut record file. */
	if(errNum == WLZ_ERR_NONE && jRecTlPl)
	{
	  cJSON *jIndex = NULL;

	  if(((jIndex  = cJSON_GetObjectItemCaseSensitive(
	           jRecTl, "index")) == NULL) ||
             (cJSON_IsNumber(jIndex) == 0))
	  {
	    errNum = WLZ_ERR_FILE_FORMAT;
	  }
	  else
	  {
	    tlIndex = ALG_NINT(jIndex->valuedouble);
	  }
	}
	/* Set section distance. */
	if(errNum == WLZ_ERR_NONE)
	{
	  if(jRecTlPl)
	  {
	    cJSON *jDist = NULL;

	    if(((jDist = cJSON_GetObjectItemCaseSensitive(
		    jRecTlPl, "dist")) == NULL) ||
	       (cJSON_IsNumber(jDist) == 0))
	    {
	      errNum = WLZ_ERR_FILE_FORMAT;
	    }
	    else
	    {
	      vWSp->dist = jDist->valuedouble;
	    }
	  }
	  else
	  {
	    switch(plnType)
	    {
	      case WLZ_C2DKEY_ORTHOGONAL:
		switch(ortCnt)
		{
		  case 0:
		    vWSp->dist = vWSp->minvals.vtZ + plCnt;
		    break;
		  case 1:
		    vWSp->dist = vWSp->minvals.vtX  + plCnt;
		    break;
		  default:
		    vWSp->dist = vWSp->minvals.vtY  + plCnt;
		    break;
		}
		break;
	      case WLZ_C2DKEY_PARALLEL:
		vWSp->dist = vWSp->minvals.vtZ + plCnt;
		break;
	      case WLZ_C2DKEY_RANDOM:
		vWSp->dist = (vWSp->maxvals.vtZ - vWSp->minvals.vtZ) * 
			     AlgRandUniform() + vWSp->minvals.vtZ;
		break;
	    }
	  }
	  if(verbose)
	  {
	    (void )fprintf(stderr, "%s: Section distance = %g\n",
	    argv[0], vWSp->dist);
	  }
	}
	/* Cut section plane. */
	if(errNum == WLZ_ERR_NONE)
	{
	  if(verbose)
	  {
	    (void )fprintf(stderr, "%s: Initialising viewstruct with angles\n"
	    "\t\t phi = %g, theta = %g, zeta = %g (Radians),\n"
	    "\t\tfixed point %g,%g,%g, distance %g\n",
	    argv[0],
	    vWSp->phi, vWSp->theta, vWSp->zeta,
	    vWSp->fixed.vtX, vWSp->fixed.vtY, vWSp->fixed.vtZ,
	    vWSp->dist);
	  }
	  errNum = WlzInit3DViewStruct(vWSp, inObj);
	  if(errNum == WLZ_ERR_NONE)
          {
	    if(verbose)
	    {
	      (void )fprintf(stderr, "%s: Cutting plane from input object\n",
	      argv[0]);
	    }
	    plObj = WlzAssignObject(
		    WlzGetSubSectionFromObject(inObj, NULL, vWSp, interp, NULL,
					       &errNum), NULL);
	  }
          if(errNum == WLZ_ERR_NONE)
	  {
	    emptyPl = WlzIsEmpty(plObj, NULL);
	    if(!emptyPl)
	    {
	     if(!plObj || !(plObj->values.core))
	      {
		emptyPl = 1;
	      }
	      else
	      {
		WlzPixelV minVal,
			  maxVal;

		errNum = WlzGreyRange(plObj, &minVal, &maxVal);
		if(errNum == WLZ_ERR_NONE)
		{
		  double del;

		  (void )WlzValueConvertPixel(&minVal, minVal, WLZ_GREY_DOUBLE);
		  (void )WlzValueConvertPixel(&maxVal, maxVal, WLZ_GREY_DOUBLE);
		  maxVal.v.dbv -= bgdVal;
		  minVal.v.dbv -= bgdVal;
		  if((fabs(minVal.v.dbv) < 1.0 - ALG_DBL_TOLLERANCE) &&
		     (fabs(maxVal.v.dbv) < 1.0 - ALG_DBL_TOLLERANCE))
		  {
		    emptyPl = 1;
		  }
		}
	      }
	    }
	  }
	}
	/* End of getting new plane. */
      }
      /* Create tile. */
      if(errNum == WLZ_ERR_NONE)
      {
	WlzObject *rotPlObj = NULL;
	WlzAffineTransform *tlRotTr = NULL;

	/* Within plane rotation. */
	if(jRecTl)
	{
	  cJSON *jRotation;

	  if(((jRotation = cJSON_GetObjectItemCaseSensitive(
		  jRecTl, "rotation")) == NULL) ||
	     (cJSON_IsNumber(jRotation) == 0))
	  {
	    errNum = WLZ_ERR_FILE_FORMAT;
	  }
	  else
 	  {
	    tlRotAngle = jRotation->valuedouble;
	  }
	}
	else
	{
	  if(emptyPl)
	  {
	    emptyTl = 1;
	  }
	  else
	  {
	    switch(rotType)
	    {
	      case WLZ_C2DKEY_NONE:
		tlRotAngle = 0.0;
		break;
	      case WLZ_C2DKEY_RANDOM:
		/* Rotation through random angle about centre of plane */
		tlRotAngle = 2.0 * ALG_M_PI * AlgRandUniform();
		break;
	      default:
		errNum = WLZ_ERR_PARAM_DATA;
		break;
	    }
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  if(verbose)
	  {
	    (void )fprintf(stderr, "%s: Rotating cut plane by %g Radians\n",
		argv[0], tlRotAngle);
	  }
	  if(tlRotAngle * tlRotAngle > ALG_DBL_TOLLERANCE)
          {
	    (void )WlzFreeAffineTransform(tlRotTr);
	    tlRotTr = WlzAffineTransformFromRotation(
                WLZ_TRANSFORM_2D_AFFINE, 0.0, 0.0, tlRotAngle, &errNum);
            if(errNum == WLZ_ERR_NONE)
            {
	      rotPlObj = WlzAssignObject(WlzAffineTransformObj(plObj,
	          tlRotTr, interp, &errNum), NULL);
	    }
	  }
	  else
	  {
	    rotPlObj = WlzAssignObject(plObj, NULL);
	  }
	}
	/* Within plane clipping. */
	if((errNum == WLZ_ERR_NONE) && !emptyTl)
        {
	  if(jRecTlPl)
	  {
	    cJSON *jClipbox;
	    double boxVals[4];

	    if(((jClipbox = cJSON_GetObjectItemCaseSensitive(
		    jRecTl, "clipbox")) == NULL) ||
	       (cJSON_IsArray(jClipbox) == 0))
	    {
	      errNum = WLZ_ERR_FILE_FORMAT;
	    }
	    else
	    {
	      for(int i = 0; i < 4; ++i)
	      {
	        cJSON *jV;
		if(((jV = cJSON_GetArrayItem( jClipbox, i)) == NULL) ||
		   (cJSON_IsNumber(jV) == 0))
	        {
		  errNum = WLZ_ERR_FILE_FORMAT;
		  break;
		}
		boxVals[i] = jV->valuedouble;
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      tlBox.xMin = ALG_NINT(boxVals[0]);
	      tlBox.yMin = ALG_NINT(boxVals[1]);
	      tlBox.xMax = ALG_NINT(boxVals[2]);
	      tlBox.yMax = ALG_NINT(boxVals[3]);
	    }
	  }
	  else
	  {
	    rotPlBBox = WlzBoundingBox2I(rotPlObj, &errNum);
	    if(verbose)
	    {
	      (void )fprintf(stderr,
		  "%s: Rotated plane has bounding box %d,%d,%d,%d\n",
		  argv[0], rotPlBBox.xMin, rotPlBBox.yMin,
		  rotPlBBox.xMax, rotPlBBox.yMax);
	    }
	    if((rotPlBBox.xMax - rotPlBBox.xMin < minFgSz) ||
	       (rotPlBBox.yMax - rotPlBBox.yMin < minFgSz))
	    {
	      emptyTl = 1;
	    }
	    else
	    {
	      switch(plnTilType)
	      {
		case WLZ_C2DKEY_RASTER:
		  if(inPlTlCnt == 0)
		  {
		    tlBox.xMin = rotPlBBox.xMin;
		    tlBox.yMin = rotPlBBox.yMin;
		  }
		  else
		  {
		    tlBox.xMin += tileSz.vtX - tileOverlap.vtX;
		    if(tlBox.xMin > rotPlBBox.xMax)
		    {
		      tlBox.xMin = rotPlBBox.xMin;
		      tlBox.yMin  += tileSz.vtY  - tileOverlap.vtY;
		    }
		  }
		  tlBox.xMax = tlBox.xMin + tileSz.vtX - 1;
		  tlBox.yMax = tlBox.yMin + tileSz.vtY - 1;
		  break;
		case WLZ_C2DKEY_RANDOM:
		  tlBox.xMin = rotPlBBox.xMin - (tileSz.vtX / 2) + 
			   (rotPlBBox.xMax - rotPlBBox.xMin) * AlgRandUniform();
		  tlBox.yMin = rotPlBBox.yMin - (tileSz.vtY / 2) + 
			   (rotPlBBox.yMax - rotPlBBox.yMin) * AlgRandUniform();
		  tlBox.xMax = tlBox.xMin + tileSz.vtX - 1;
		  tlBox.yMax = tlBox.yMin + tileSz.vtY - 1;
		  break;
		default:
		  errNum = WLZ_ERR_PARAM_DATA;
		  break;
	      }
	    }
	  }
	}
	if((errNum == WLZ_ERR_NONE) && !emptyTl)
        {
	  if((tlBox.xMax - tlBox.xMin < minFgSz) ||
             (tlBox.yMax - tlBox.yMin < minFgSz))
          {
	    emptyTl = 1;
	  }
	  else
	  {
	    if(verbose)
	    {
	      (void )fprintf(stderr,
		  "%s: Clipping rotated plane to tile %d,%d,%d,%d\n",
		  argv[0], tlBox.xMin, tlBox.yMin, tlBox.xMax, tlBox.yMax);
	    }
	    clpRotPlObj = WlzAssignObject(
		    WlzClipObjToBox2D(rotPlObj, tlBox, &errNum), NULL);
	    if(WlzIsEmpty(clpRotPlObj, NULL))
	    {
	      emptyTl = 1;
	    }
	  }
	}
	if((errNum == WLZ_ERR_NONE) && !emptyTl &&
	   (clpRotPlObj->values.core == NULL))
        {
	  WlzPixelV	bgd,
	  		val;
	  WlzValues	newVal;

	  if(verbose)
          {
	    (void )fprintf(stderr,
                "%s: Tile has no values, setting values within domain to %d,\n"
		"\t\tbackground value 0\n",
		argv[0], inDomVal);
	  }
	  /* If the tile object has no values, create some. */
	  bgd.type = WLZ_GREY_UBYTE;
	  bgd.v.ubv = bgdVal;
	  newVal.v  = WlzNewValueTb(clpRotPlObj, WLZ_GREY_UBYTE, bgd, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    clpRotPlObj->values = WlzAssignValues(newVal, NULL);
	    val.type = WLZ_GREY_UBYTE;
	    val.v.ubv = inDomVal;
	    errNum = WlzGreySetValue(clpRotPlObj, val);
	  }
	}
	if((errNum == WLZ_ERR_NONE) && !emptyTl)
	{
	  if(verbose)
	  {
	    (void )fprintf(stderr,
		"%s: Cutting rectangular tile\n",
		argv[0]);
	  }
	  tlObj = WlzAssignObject(
		  WlzCutObjToBox2D(clpRotPlObj, tlBox, WLZ_GREY_UBYTE,
		      0, 0.0, 0.0, &errNum), NULL);
	}
        (void )WlzFreeObj(rotPlObj);
      }
      /* Set tile file name and write tile to file. */
      if(errNum == WLZ_ERR_NONE)
      {
        const char *ext;

        (void )WlzEffStringFromFormat(tileFileFmt, &ext);
	(void )sprintf(filenameBuf, "%s_%06d.%s",
	           baseFile, (jRecTiles)? tlIndex: tlCnt, ext);
	if(emptyTl)
	{
	  if(verbose)
          {
	    (void )fprintf(stderr,
	        "%s: Tile is empty, no tile image file output\n",
		argv[0]);
	  }
	}
	else
	{
	  if(verbose)
          {
	    (void )fprintf(stderr,
	        "%s: Outputting tile image to file %s\n",
		argv[0], filenameBuf);
	  }
	  errNum = WlzEffWriteObj(NULL, filenameBuf, tlObj, tileFileFmt);
	}
      }
      (void )WlzFreeObj(tlObj); tlObj = NULL;
      /* Append the cut tile record file. */
      if(errNum == WLZ_ERR_NONE)
      {
	if((!emptyTl || includeEmpty))
	{
	  if(emptyTl)
	  {
	    tlRotAngle = 0.0;
	    tlBox.xMin = tlBox.yMin = tlBox.xMax = tlBox.yMax = 0;
	  }
	  if(verbose)
	  {
	    (void )fprintf(stderr,
		"%s: Appending record to record file\n",
		argv[0]);
	  }
	  if(fprintf(outRecFP,
	      "{\n"
	      "\"file\": \"%s\",\n"
	      "\"index\": %d,\n"
	      "\"plane\":\n"
	      "{\n"
	      "\"phi\": %g,\n"
	      "\"theta\": %g,\n"
	      "\"zeta\": %g,\n"
	      "\"up\": [%g, %g, %g],\n"
	      "\"fixed\": [%g, %g, %g],\n"
	      "\"dist\": %g,\n"
	      "\"scale\": %g,\n"
	      "\"view_mode\": %d,\n"
	      "\"voxelrescale\": %d\n"
	      "},\n"
	      "\"empty\": %s,\n"
	      "\"rotation\": %g,\n"
	      "\"clipbox\": [%d, %d, %d, %d]\n"
	      "}",
	      filenameBuf,
	      (jRecTiles)? tlIndex: tlCnt,
	      vWSp->phi,
	      vWSp->theta,
	      vWSp->zeta,
	      vWSp->up.vtX, vWSp->up.vtY, vWSp->up.vtZ,
	      vWSp->fixed.vtX, vWSp->fixed.vtY, vWSp->fixed.vtZ,
	      vWSp->dist,
	      vWSp->scale,
	      vWSp->view_mode,
	      vWSp->voxelRescaleFlg,
	      keyStrings[emptyTl != 0],
	      tlRotAngle,
	      tlBox.xMin, tlBox.yMin, tlBox.xMax, tlBox.yMax) < 0)
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(jRecTiles)
	{
	  if(((maxTiles > 0) && (allTlCnt > maxTiles)) ||
	     (allTlCnt > maxRecTiles))
	  {
	    errNum = WLZ_ERR_EOO;
	  }
	}
	else
	{
	  ++tlCnt;
	  ++inPlTlCnt;
	  if(((maxTilesPP > 0) && (inPlTlCnt >= maxTilesPP)) ||
	     ((plnTilType == WLZ_C2DKEY_RASTER) &&
	      (tlBox.yMin + tileSz.vtY > rotPlBBox.yMax)))
	  {
	    inPlTlCnt = 0;
	    switch(plnType)
	    {
	      case WLZ_C2DKEY_PARALLEL:
		if(vWSp->dist > vWSp->maxvals.vtZ)
		{
		  errNum = WLZ_ERR_EOO;
		}
		break;
	      case WLZ_C2DKEY_ORTHOGONAL:
		switch(ortCnt)
		{
		  case 0:
		    if(vWSp->dist > vWSp->maxvals.vtZ)
		    {
		      plCnt = 0;
		      ++ortCnt;
		    }
		    break;
		  case 1:
		    if(vWSp->dist > vWSp->maxvals.vtX)
		    {
		      plCnt = 0;
		      ++ortCnt;
		    }
		    break;
		  default:
		    if(vWSp->dist > vWSp->maxvals.vtY)
		    {
		      errNum = WLZ_ERR_EOO;
		    }
		    break;
		}
		break;
	      default:
		break;
	    }
	  }
	  if(inPlTlCnt == 0)
	  {
	    plCnt += distInc;
	  }
	  ++tlCnt;
	  if((maxTiles > 0) && (tlCnt > maxTiles))
	  {
	    errNum = WLZ_ERR_EOO;
	  }
	  ++allTlCnt;
	  if((maxTiles > 0) && (allTlCnt > 2 * maxTiles))
	  {
	    errNum = WLZ_ERR_EOO;
	  }
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
        if(!emptyTl || includeEmpty)
	{
	  (void )fprintf(outRecFP, ",\n");
	}
      }
      else
      {
	(void )fprintf(outRecFP, "]\n}\n");
      }
      (void )WlzFreeObj(clpRotPlObj);
    }
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
    (void )WlzFreeObj(plObj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: Error cutting tiles, tile count %d "
          "(%s)\n",
          argv[0], tlCnt, errMsg);
    }
  }
  cJSON_Delete(jRec);
  (void )WlzFreeObj(inObj);
  (void )WlzFree3DViewStruct(vWSp);
  (void )AlcFree(inRecStr);
  (void )AlcFree(filenameBuf);
  if(outRecFP)
  {
    (void )fclose(outRecFP);
  }
  if(usage)
  {
    const char *inObjFileFmtStr;

    if(WlzEffStringFromFormat(inObjFileFmt, &inObjFileFmtStr) == NULL)
    {
      inObjFileFmtStr = "unknown";
    }
    (void )fprintf(stderr,
    "Usage: %s [-h] [-v] [-a(n|r)]\n"
    "\t\t[-b<background tile value>] [-c<record file>]\n"
    "\t\t[-e<columns>,<lines>]\n"
    "\t\t[-E<input file format>] [-f<file base>]\n"
    "\t\t[-F<output tile format>]\n" 
    "\t\t[-i]<in domain value>] [-M<max tile images>]\n"
    "\t\t[-N<max tile images per plane>] [-p(o|p|r)]\n"
    "\t\t[-r<pitch>,<yaw>,<roll>] [-R<rescale flags>\n"
    "\t\t[-s<columns>,<lines>] [-t(ras|ran)] [-Y] [<input object>]\n"
    "Cuts 2D region of interest tiles from a 3D domain object.\n"
    "Version: %s\n"
    "Options:\n"
    "  -h  Help, prints this usage message\n"
    "  -v  Verbose output messages, set to %s\n"
    "  -a  Type of within plane rotation (about the plane centre) - none\n"
    "      (n),  random (r), set to %s\n"
    "  -b  Background tile value, set to %d\n"
    "  -c  Cut all tile images using the given record file, set to %s\n"
    "  -d  Distance increment for parallel and orthogonal planes, must\n"
    "      be > 0, set to %d\n"
    "  -e  Tile overlap for raster tiling, set to %d,%d\n"
    "  -E  Input object file format, set to %s\n"
    "  -f  Output base file name, set to \"%s\"\n"
    "  -F  Output image tile format, set to %s\n"
    "  -L  Use linear interpolation instead of nearest neighbour when\n"
    "        cutting planes, set to %s\n"
    "  -i  Within domain tile value, range 1-255, set to %d\n"
    "  -M  Maximum number of tile images, set to %d\n"
    "  -N  Maximum number of tile images for each plane cut (the same\n"
    "      plane may possibly be cut more than once), set to %d\n"
    "  -p  Type of planes to cut - orthogonal (o), parallel (p) or\n"
    "      random (r). The orthogonal and parallel planes are with respect\n"
    "      to the referece plane, set to %s\n"
    "  -r  Reference plane angles (given using degrees), defined by pitch,\n"
    "      yaw and roll (also known as phi, theta and zeta), set to\n"
    "      %g,%g,%g\n"
    "  -R  Voxel rescaling flags:\n"
    "        bit 1 set - use voxel-size rescaling,\n"
    "        bit 2 set - enable global scaling,\n"
    "      set to %d\n"
    "  -s  Size of tile images to cut, set to %d,%d\n"
    "  -t  Within plane section tiling - raster (s), random (r), set to %s\n" 
    "  -Y  Include empty tiles (entirely outside the input domain or\n"
    "      of background) in output, set to %s\n"
    "Examples:\n"
    "  %s -f tiles in.wlz\n"
    "Cuts 2D tiles from the input file in.wlz writing the tile cut record\n"
    "to tiles_cut.jsn and tile images to tiles_000000.wlz, etc....\n"
    "  %s -f tiles -d prev_cut.jsn in.wlz\n"
    "Cuts 2D tiles from the input file in.wlz as with the previous\n"
    "example except that the cutting parameters are read from the given file\n"
    "prev_cut.jsn\n"
    "  %s -s 64,128 -p r -a r -f tiles -M 1000 -F jpg in.wlz\n"
    "Cuts 2D 64x128 tiles from the input file in.wlz using random planes\n"
    "and random in plane rotations, writing the tile cut record to\n"
    "tiles_cut.jsn and tile images to tiles_000000.jpg, ...,\n"
    "tiles_000999.jpg.\n",
    argv[0],
    WlzVersion(),
    keyStrings[verbose != 0],
    keyStrings[rotType],
    bgdVal,
    (inRecFile)? inRecFile: "false",
    distInc,
    tileOverlap.vtX, tileOverlap.vtY,
    inObjFileFmtStr,
    baseFile,
    WlzEffStringFromFormat(tileFileFmt, NULL),
    keyStrings[interp != 0],
    inDomVal,
    maxTiles,
    maxTilesPP,
    keyStrings[plnType],
    refPlnAngles.vtX * 180 / ALG_M_PI,
    refPlnAngles.vtY * 180 / ALG_M_PI,
    refPlnAngles.vtZ * 180 / ALG_M_PI,
    voxRescale,
    tileSz.vtX, tileSz.vtY,
    keyStrings[plnTilType],
    keyStrings[includeEmpty != 0],
    argv[0],
    argv[0],
    argv[0]);
    (void )fprintf(stderr,
    "Recognised file formats are:\n"
    "  Description                                       Extension\n"
    "  ***********                                       *********\n%s\n",
    WlzEffFormatTable(2, 50, 10, NULL));
  }
  return(!ok);
}

/*!
* \return	Index into keyStrings which matches the given param or
* 		-1 on no match.
* \brief	Finds first string in keyStrings which matches trhe given
* 		parameter string (up to the length of the parameter string).
* \param	par		Given parameter string.
*/
static int 	ParseStrOption(const char *par)
{
  int		n,
  		match = -1;

  if((n = strlen(par)) > 0)
  {
    for(int i = 0; i < WLZ_C2DKEY_COUNT; ++i)
    {
      if(strncmp(par, keyStrings[i], n) == 0)
      {
	match = i;
	break;
      }
    }
  }
  return(match);
}


/*!
* \return	Normalised angle.
* \brief	Normalises the given angle so that it is in the range
* 		[0.0-max].
* \param	a		Given angle.
* \param	max		Maximum angle.
*/
static double	NormaliseAngle(double a, double max)
{
  if(a > max)
  {
    double 	t;

    t = floor(a / max);
    if(t > 0)
    {
      a -= t * max;
    }
  }
  return(a);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
