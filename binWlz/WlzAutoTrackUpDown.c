#ifndef DOXYGEN_SHOULD_SKIP_THIS
#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzAutoTrackUpDown_c[] = "MRC HGU $Id$";
#endif
#endif
/***********************************************************************
* Project:      Woolz
* Title:        WlzAutoTrackUpDown.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Tracks a boundary up and down from the defined planes.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <Alg.h>
#include <Wlz.h>

typedef struct {
  double	sintheta;
  double	costheta;
  int		index;
  int		x;
  int		y;
  int		x_off;
  int		y_off;
} MatchPointStruct;

static double	nu_dist, nu_alpha, nu_kappa;
static int	spacing_PMSnake = 5;
static int	range_PMSnake = 5;
static double	nu_dist_PMSnake = 0.5;
static double	nu_alpha_PMSnake = 0.0;
static double	nu_kappa_PMSnake = 0.1;

void PMSnakeNlcSetup(
  int		spacing,
  int		range,
  double	nu_alpha,
  double	nu_kappa)
{
  spacing_PMSnake = spacing;
  range_PMSnake = range;
  nu_alpha_PMSnake = nu_alpha;
  nu_kappa_PMSnake = nu_kappa;

  return;
}

double PMSnakeNlc(
  int	i,
  int	j,
  int	jp,
  int	**op)
{
  double	cost, alpha, kappa;

  /* angle cost */
  alpha = (double) (j - jp) / spacing_PMSnake;
  cost = nu_alpha_PMSnake * alpha * alpha;

  /* curvature cost */
  kappa = 0.0;
  if( i > 0 )
  {
    kappa = (double) (j - 2*jp + op[i-1][jp]) / spacing_PMSnake /
      spacing_PMSnake;
  }
  cost += nu_kappa_PMSnake * kappa * kappa;

  return( cost );
}

void PMSnake(
  double		**local_cost,
  int			num_mpts,
  int			spacing,
  int			range,
  MatchPointStruct	*mpts,
  int			closed_path)
{
  double	**optimal_cost;
  int		**optimal_path;
  int		i, j, opt_j;
  double	dist_factor;
  double	min_cost;

  (void) AlcDouble2Calloc(&optimal_cost, num_mpts, 2*range+1);
  (void) AlcInt2Calloc(&optimal_path, num_mpts, 2*range+1);

  /* get the cost parameters from the dialog sliders */
  /* now set as globals */
  /*fprintf(stderr, "nu_dist, nu_alpha, nu_kappa = %f, %f, %f\n",
	nu_dist, nu_alpha, nu_kappa);*/

  /* add in the local distance cost - log gaussian probability
     ie quadratic */
  dist_factor = nu_dist / range / range;
  for(i=0; i < num_mpts; i++)
  {
    for(j = -range; j <= range; j++)
    {
      local_cost[i][j+range] += dist_factor * j * j;
    }
  }

  /* print out the local cost */
/*  fprintf(stderr, "num_mpts = %d, range = %d\n", num_mpts, range);
  for(i=0; i < num_mpts; i++)
  {
    for(j = -range; j <= range; j++)
    {
      fprintf(stderr, "%8.4f", local_cost[i][j+range]);
    }
    fprintf(stderr, "\n");
  }*/
  

  /* set up the parameters for the non-local cost */
  PMSnakeNlcSetup(spacing, range, nu_alpha, nu_kappa);

  /* call the dynamic programming search */
  (void) AlgDPSearch(num_mpts, 2*range+1, local_cost, optimal_cost,
		       optimal_path, PMSnakeNlc);

  /* extract the optimal path and calculate offsets */
  opt_j = 0;
  min_cost = optimal_cost[num_mpts-1][opt_j];
  for(j = -range; j <= range; j++)
  {
    if( optimal_cost[num_mpts-1][j+range] < min_cost )
    {
      min_cost = optimal_cost[num_mpts-1][j+range];
      opt_j = j+range;
    }
  }
  for(i=num_mpts-1; i >= 0; i--)
  {
    mpts[i].x_off = WLZ_NINT((opt_j - range) * mpts[i].costheta);
    mpts[i].y_off = WLZ_NINT((opt_j - range) * mpts[i].sintheta);

    /* print path and offsets */
/*    fprintf(stderr, "i = %d, opt_j = %d, (x_off, y_off) = (%d, %d)\n",
	    i, opt_j, mpts[i].x_off, mpts[i].y_off);*/

    opt_j = optimal_path[i][opt_j];
  }

  return;
}

static double edge_match_cost(
  WlzObject	*test_obj,
  WlzObject	*ref_obj,
  int		dx,
  int		dy)
{
  int		l, k;
  double	cost;
  WlzGreyP	gPix;

  /* test object pointers, if NULL return very large cost
     note 1.0 is the theoretical maximum */
  if( (test_obj == NULL) || (ref_obj == NULL) ){
    return( 10.0 );
  }

  l = ref_obj->domain.i->line1;
  k = ref_obj->domain.i->kol1;
  WlzGreyValue(test_obj, 0, l+dy, k+dx, &gPix);

  switch( WlzGreyTableTypeToGreyType(test_obj->values.core->type, NULL) )
  {
   default:
   case WLZ_GREY_INT:
     cost = *gPix.inp;
     break;
   case WLZ_GREY_SHORT:
     cost = *gPix.shp;
     break;
   case WLZ_GREY_UBYTE:
     cost = *gPix.ubp;
     break;
   case WLZ_GREY_FLOAT:
     cost = *gPix.flp;
     break;
   case WLZ_GREY_DOUBLE:
     cost = *gPix.dbp;
     break;
  }

  cost = 1.0 / (1.0 + cost);
  return( cost );
}

static double image_match_cost(
  WlzObject	*test_obj,
  WlzObject	*ref_obj,
  int		dx,
  int		dy)
{
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace	gwsp;
  WlzGreyP	grey;
  double	cost;
  int		i, a;
  WlzGreyValueWSpace *gVWSp = NULL;

  /* test object pointers, if NULL return very large cost
     note 1.0 is the theoretical maximum */
  if( (test_obj == NULL) || (ref_obj == NULL) ){
    return( 10.0 );
  }

  if( dx < -10 || dx > 10 || dy < -10 || dy > 10 ){
    fprintf(stderr, "unlikely values for offsets: (dx, dy)=(%d, %d)\n",
	    dx, dy);
  }

  /* scan through obj2 comparing values with offset obj1 */
  WlzInitGreyScan(ref_obj, &iwsp, &gwsp);
  if((gVWSp = WlzGreyValueMakeWSp(test_obj, NULL)) == NULL){
    return( 10.0 );
  }
  a = 0;
  cost = 0.0;
  while( WlzNextGreyInterval( &iwsp ) == WLZ_ERR_NONE )
  {
/*    switch( gwsp.pixeltype )
    {
     default:
     case INT_GREY:
       break;
     case SHORT_GREY:
       break;
     case UBYTE_GREY:
       break;
     case FLOAT_GREY:
       break;
     case DOUBLE_GREY:
       break;
    }
    */
    for(i=0; i < iwsp.colrmn; i++)
    {
      int	g1, g2;

      WlzGreyValueGet(gVWSp, 0, iwsp.linpos+dy, iwsp.colpos+i+dx);
      g1 = *(*gVWSp->gPtr).ubp;
      g2 = *gwsp.u_grintptr.ubp;
      cost += (double) ((g1 > g2) ? g1 - g2 : g2 - g1);
      gwsp.u_grintptr.ubp++;
      a += 1;
    }
  }
  if(gVWSp)
  {
    WlzGreyValueFreeWSp(gVWSp);
  }

  /* normalise by the object area and return */
  return( cost/a/256.0 );
}

static WlzObject *get_rect_region(
  WlzObject	*obj,
  int		line,
  int		kol,
  int		ln_size,
  int		kl_size)
{
  WlzObject	*rect_obj, *tmp_obj;
  WlzDomain	domain;
  WlzValues	values;

  domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
				   line - ln_size, line + ln_size,
				   kol - kl_size, kol + kl_size,
				   NULL);
  values.core = NULL;
  tmp_obj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
		     NULL, NULL, NULL);
  if( (rect_obj = WlzIntersect2(obj, tmp_obj, NULL)) ){
    rect_obj->values = WlzAssignValues(obj->values, NULL);
  }
  WlzFreeObj( tmp_obj );

  return( rect_obj );
}

static WlzPolygonDomain *HGU_TrackPolyline(
  WlzObject	*ref_obj,
  WlzObject	*test_obj,
  WlzPolygonDomain	*ref_polydmn,
  int		spacing,
  int		range,
  int		size,
  double	(*match_cost)(WlzObject *, WlzObject *, int, int))
{
  WlzObject		*poly_obj, *ref_region;
  WlzDomain		domain;
  WlzValues		values;
  WlzPolygonDomain	*new_polydmn;
  MatchPointStruct	*mpts;
  int			i, j, k, kmax, num_mpts;
  double		dx1, dx2, dx3, dy1, dy2, dy3, s1, s2, s3;
  int			x, y, dx, dy;
  WlzIVertex2		*vtxs;
  double		cost, **local_cost;

  /* create an 8-connected polygon */
  poly_obj = WlzPolyTo8Polygon( ref_polydmn, 1, NULL );
  new_polydmn = WlzAssignPolygonDomain(poly_obj->domain.poly, NULL);
  WlzFreeObj( poly_obj );
  new_polydmn->linkcount = 0;

  /* allocate space for the match point */
  /* check for sufficient polyline points */
  if( (num_mpts = new_polydmn->nvertices / spacing) < 5 )
  {
    return( new_polydmn );
  }
  /* make sure the last point is not matched because
     it is identical to the first */
  if( (new_polydmn->nvertices % spacing) < 2 )
  {
    num_mpts--;
  }
  mpts = (MatchPointStruct *) AlcMalloc(sizeof(MatchPointStruct) *
					(num_mpts+1));
  (void) AlcDouble2Calloc( &local_cost, num_mpts+1, 2*range+1 );

  /* scan through match points */
  for(i=0; i < num_mpts; i++)
  {
    mpts[i].index = i * spacing;
    vtxs = new_polydmn->vtx;
    x = vtxs[mpts[i].index].vtX;
    y = vtxs[mpts[i].index].vtY;
    mpts[i].x = x;
    mpts[i].y = y;

    /* simple weighted angle, note eight_poly returns the last
       vertex equal to the first when the polyline is closed */
    j = (mpts[i].index > 0) ? mpts[i].index :
      mpts[i].index + new_polydmn->nvertices - 1;

    dx1 = vtxs[mpts[i].index+1].vtX - vtxs[j-1].vtX;
    dx2 = vtxs[mpts[i].index+2].vtX - vtxs[j-2].vtX;
    dx3 = vtxs[mpts[i].index+3].vtX - vtxs[j-3].vtX;
    dy1 = vtxs[mpts[i].index+1].vtY - vtxs[j-1].vtY;
    dy2 = vtxs[mpts[i].index+2].vtY - vtxs[j-2].vtY;
    dy3 = vtxs[mpts[i].index+3].vtY - vtxs[j-3].vtY;
    s1 = sqrt(dx1*dx1 + dy1*dy1);
    s2 = sqrt(dx2*dx2 + dy2*dy2);
    s3 = sqrt(dx3*dx3 + dy3*dy3);
    s1 = (s1 < DBL_EPSILON)?1.0:s1;
    s2 = (s2 < DBL_EPSILON)?1.0:s2;
    s3 = (s3 < DBL_EPSILON)?1.0:s3;

    /* note we want the perpendicular direction */
    mpts[i].sintheta =  ( dx1/s1 + dx2/s2 + dx3/s3 ) / 3.0;
    mpts[i].costheta = -( dy1/s1 + dy2/s2 + dy3/s3 ) / 3.0;

    /* get a region from the reference image */
    ref_region = get_rect_region(ref_obj, y, x, size, size);

    /* scan through search region to generate the cost matrix */
    for(j = -range; j <= range; j++)
    {
      dx = WLZ_NINT(j * mpts[i].costheta);
      dy = WLZ_NINT(j * mpts[i].sintheta);

      /* calculate match cost */
      local_cost[i][j+range] = (*match_cost)(test_obj, ref_region, dx, dy);
    }

    if( ref_region ){
      WlzFreeObj( ref_region );
    }
  }

  /* first and last points are the same */
  mpts[num_mpts] = mpts[0];
  mpts[num_mpts].index = new_polydmn->nvertices - 1;
  for(j = -range; j <= range; j++)
  {
    local_cost[num_mpts][j+range] = local_cost[0][j+range];
  }

  /* determine offsets by minimising the global cost.
     This uses dynamic programming to include a local distance cost
     plus a non-local angle cost to mimic the snake behavior.
     Note this routine assumes that the match cost is normalised */
  PMSnake( local_cost, num_mpts+1, spacing, range, mpts, 1 );
      
  /* calculate new polyline */
  for(i=0; i < num_mpts; i++)
  {
    dx1 = mpts[i].x_off;
    dy1 = mpts[i].y_off;
    dx2 = mpts[i+1].x_off;
    dy2 = mpts[i+1].y_off;
    kmax = mpts[i+1].index - mpts[i].index;
    for(j=mpts[i].index, k=0; j < mpts[i+1].index; j++, k++)
    {
      dx = WLZ_NINT( ((kmax - k) * dx1 + k * dx2) / kmax );
      dy = WLZ_NINT( ((kmax - k) * dy1 + k * dy2) / kmax );
      vtxs[j].vtX += dx;
      vtxs[j].vtY += dy;
    }
  }
  vtxs[j] = vtxs[0];
  

  /* strip repeated vertices */
  for(i=1, j=1; i < new_polydmn->nvertices; i++)
  {
    if( (vtxs[i].vtX == vtxs[i-1].vtX) && (vtxs[i].vtY == vtxs[i-1].vtY) )
    {
      continue;
    }
    vtxs[j] = vtxs[i];
    j++;
  }
  new_polydmn->nvertices = j;

  /* clean up and return */
  AlcFree( (void *) mpts );
  return( new_polydmn );
}

static WlzBoundList *HGU_TrackBoundlist(
  WlzObject	*ref_obj,
  WlzObject	*new_obj,
  WlzBoundList	*ref_boundlist,
  int		spacing,
  int		range,
  int		size,
  double	(*match_cost)(WlzObject *, WlzObject *, int, int))
{
  WlzBoundList	*new_bndlst;

  /* check boundlist for stopping condition */
  if( ref_boundlist == NULL )
  {
    return( NULL );
  }

  /* create a new boundlist structure */
  new_bndlst = (WlzBoundList *) AlcCalloc(1, sizeof(WlzBoundList));
  new_bndlst->type = ref_boundlist->type;
  new_bndlst->linkcount = 0;
  new_bndlst->freeptr = NULL;
  new_bndlst->up = NULL;
  new_bndlst->next =
    WlzAssignBoundList(HGU_TrackBoundlist(ref_obj, new_obj,
					  ref_boundlist->next,
					  spacing, range, size,
					  match_cost), NULL);
  new_bndlst->down =
    WlzAssignBoundList(HGU_TrackBoundlist(ref_obj, new_obj,
					  ref_boundlist->down,
					  spacing, range, size,
					  match_cost), NULL);
  new_bndlst->wrap = 1;

  /* track this polyline */
  new_bndlst->poly = 
    WlzAssignPolygonDomain(HGU_TrackPolyline(ref_obj, new_obj,
					     ref_boundlist->poly,
					     spacing, range, size,
					     match_cost), NULL);

  return( new_bndlst );
}

WlzObject *HGU_TrackDomain(
  WlzObject	*ref_obj,
  WlzObject	*new_obj,
  WlzObject	*ref_domain,
  int		spacing,
  int		range,
  int		size,
  double	(*match_cost)(WlzObject *, WlzObject *, int, int))
{	
  WlzObject	*new_domain, *ref_bound;
  WlzBoundList	*new_bndlst;


  /* check objects */
  if( (ref_obj == NULL) || (new_obj == NULL) || (ref_domain == NULL) )
  {
    return( NULL );
  }

  /* get reference domain boundary */
  ref_bound = WlzObjToBoundary( ref_domain, 1, NULL );
  new_bndlst = HGU_TrackBoundlist(ref_obj, new_obj,
				 ref_bound->domain.b,
				 spacing, range, size, match_cost);

  /* create the new domain */
  new_domain = WlzBoundToObj( new_bndlst, WLZ_SIMPLE_FILL, NULL );
  WlzFreeObj( ref_bound );
  WlzFreeBoundList(new_bndlst);

  return( new_domain );
}

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s [-r#] [-s#] [-S#] [-a#] [-c#] [-d#] [-D#] [-U#]\n"
	  "\t\t[-h] [-v] [<domain file> [<grey file>]]\n"
	  "\tTrack a boundary up and down from the defined planes in\n"
	  "\tthe given domain using the poor-mans snake as in paint.\n"
	  "\tNote it is assumed that the input objects are 3D and have\n"
	  "\tappropriately matching planes etc.\n"
	  "\tOptions are:\n"
	  "\t  -r#       range - seach distance along perps. (def 5)\n"
	  "\t  -s#       spacing - spacing between perps. (def 5)\n"
	  "\t  -S#       size - side of image match square=(2S+1) (def S=4)\n"
	  "\t  -a#       angle cost parameter (def 0.0)\n"
	  "\t  -c#       curvature cost parameter (def 0.1)\n"
	  "\t  -d#       distance cost parameter (def 0.5)\n"
	  "\t  -D#       number of planes to track in down direction (def 1)\n"
	  "\t  -U#       number of planes to track in up direction (def 1)\n"
	  "\t  -h        help - prints this usage message\n"
	  "\t  -v        verbose operation\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*greyObj, *dmnObj, *newObj;
  FILE		*inFile;
  char 		optList[] = "r:s:S:a:c:d:D:U:hv";
  int		option;
  int		verboseFlg=0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  int		range, spacing, size;
  int		downTrk, upTrk;
  int		p, dp;
  WlzDomain	domain, *greyDoms, *dmnDoms, *newDoms;
  WlzObject	*refGrey, *newGrey, *refDmn, *newDmn;
  WlzValues	values, *valuess;
    
  /* set defaults, read the argument list and check for an input file */
  opterr = 0;
  range = 5;
  spacing = 5;
  size = 4;
  nu_dist = 0.5;
  nu_alpha = 0.0;
  nu_kappa = 0.1;
  downTrk = 1;
  upTrk = 1;
  
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'r':
      range = atoi(optarg);
      break;

    case 's':
      spacing = atoi(optarg);
      break;

    case 'S':
      size = atoi(optarg);
      break;

    case 'a':
      nu_alpha = atof(optarg);
      break;

    case 'c':
      nu_kappa = atof(optarg);
      break;

    case 'd':
      nu_dist = atof(optarg);
      break;

    case 'D':
      downTrk = atoi(optarg);
      break;

    case 'U':
      upTrk = atoi(optarg);
      break;

    case 'v':
      verboseFlg = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return 1;

    }
  }

  /* check for read from file, note if one file is defined then
     it is the dmnObj and the object should be on stdin */
  if( (argc - optind) >= 2 ){
    /* read the dmnObj */
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
    if( (dmnObj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) == NULL ){
      fprintf(stderr, "%s: can't read dmnObj from file %s\n", argv[0],
	      *(argv+optind));
      usage(argv[0]);
      return 1;
    }
    fclose( inFile );
    optind++;

    /* read the object */
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
    if( (greyObj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) == NULL ){
      fprintf(stderr, "%s: can't read greyObj from file %s\n", argv[0],
	      *(argv+optind));

      usage(argv[0]);
      return 1;
    }
    fclose( inFile );
  }
  else if( (argc - optind) == 1 ){
    /* read the dmnObj */
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
    if( (dmnObj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) == NULL ){
      fprintf(stderr, "%s: can't read dmnObj from file %s\n", argv[0],
	      *(argv+optind));
      usage(argv[0]);
      return 1;
    }
    fclose( inFile );
    inFile = stdin;
    if( (greyObj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) == NULL ){
      fprintf(stderr, "%s: can't read greyObj from stdin\n", argv[0]);
      usage(argv[0]);
      return 1;
    }
  }
  else {
    /* else read objects from stdin */
    if( (greyObj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) == NULL ){
      fprintf(stderr, "%s: can't read greyObj from stdin\n", argv[0]);
      usage(argv[0]);
      return 1;
    }
    if( (dmnObj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) == NULL ){
      fprintf(stderr, "%s: can't read dmnObj from stdin\n", argv[0]);
      usage(argv[0]);
      return 1;
    }
  }
    
  /* apply tracking to each non-NULL plane in the domain object
     at this point we assume that the grey and domain objects are 3D
     and have matching planes.
     */
  if((greyObj->type != WLZ_3D_DOMAINOBJ) ||
     (dmnObj->type != WLZ_3D_DOMAINOBJ) ){
    usage(argv[0]);
    return 1;
  }
  
  /* we do not track outside of the planes range of the domain object
     therefore the new domain has the same number of planes as the
     original domain object */
  domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
				dmnObj->domain.p->plane1,
				dmnObj->domain.p->lastpl,
				dmnObj->domain.p->line1,
				dmnObj->domain.p->lastln,
				dmnObj->domain.p->kol1,
				dmnObj->domain.p->lastkl,
				NULL);
  values.core = NULL;
  newObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
		       NULL, NULL, NULL);
  greyDoms = greyObj->domain.p->domains + dmnObj->domain.p->plane1 -
    greyObj->domain.p->plane1;
  valuess = greyObj->values.vox->values + dmnObj->domain.p->plane1 -
    greyObj->domain.p->plane1;
  dmnDoms = dmnObj->domain.p->domains;
  newDoms = newObj->domain.p->domains;
  for(p=domain.p->plane1; p <= domain.p->lastpl; p++,
	greyDoms++, valuess++, dmnDoms++, newDoms++){
    if( (*dmnDoms).core ){
      /* set existing planes as unchanged */
      *newDoms = WlzAssignDomain(*dmnDoms, NULL);
      

      /* track down (decreasing p) but only if pp >= plane1 */
      refGrey = WlzMakeMain(WLZ_2D_DOMAINOBJ, *greyDoms, *valuess,
			    NULL, NULL, NULL);
      refDmn = WlzMakeMain(WLZ_2D_DOMAINOBJ, *dmnDoms, values,
			   NULL, NULL, NULL);
      newGrey = WlzMakeMain(WLZ_2D_DOMAINOBJ, *greyDoms, *valuess,
			    NULL, NULL, NULL);
      newDmn = WlzMakeMain(WLZ_2D_DOMAINOBJ, *dmnDoms, values,
			   NULL, NULL, NULL);
      for(dp = -1; dp >= -downTrk; dp--){
	if( (p + dp) < domain.p->plane1 ){
	  break;
	}
	WlzFreeObj(refGrey);
	refGrey = newGrey;
	newGrey = WlzMakeMain(WLZ_2D_DOMAINOBJ, *(greyDoms+dp), *(valuess+dp),
			    NULL, NULL, NULL);
	WlzFreeObj(refDmn);
	refDmn = newDmn;
	newDmn = HGU_TrackDomain(refGrey, newGrey, refDmn,
				 spacing, range, size, image_match_cost);
	if( newDmn ){
	  *(newDoms+dp) = WlzAssignDomain(newDmn->domain, NULL);
	}
      }
      WlzFreeObj(refGrey);
      WlzFreeObj(refDmn);
      WlzFreeObj(newGrey);
      WlzFreeObj(newDmn);

      /* track up (increasing p) but only if pp <= lastpl */
      refGrey = WlzMakeMain(WLZ_2D_DOMAINOBJ, *greyDoms, *valuess,
			    NULL, NULL, NULL);
      refDmn = WlzMakeMain(WLZ_2D_DOMAINOBJ, *dmnDoms, values,
			   NULL, NULL, NULL);
      newGrey = WlzMakeMain(WLZ_2D_DOMAINOBJ, *greyDoms, *valuess,
			    NULL, NULL, NULL);
      newDmn = WlzMakeMain(WLZ_2D_DOMAINOBJ, *dmnDoms, values,
			   NULL, NULL, NULL);
      for(dp = 1; dp <= upTrk; dp++){
	if( (p + dp) > domain.p->lastpl ){
	  break;
	}
	WlzFreeObj(refGrey);
	refGrey = newGrey;
	newGrey = WlzMakeMain(WLZ_2D_DOMAINOBJ, *(greyDoms+dp), *(valuess+dp),
			    NULL, NULL, NULL);
	WlzFreeObj(refDmn);
	refDmn = newDmn;
	newDmn = HGU_TrackDomain(refGrey, newGrey, refDmn,
				 spacing, range, size, image_match_cost);
	if( newDmn ){
	  *(newDoms+dp) = WlzAssignDomain(newDmn->domain, NULL);
	}
      }
      WlzFreeObj(refGrey);
      WlzFreeObj(refDmn);
      WlzFreeObj(newGrey);
      WlzFreeObj(newDmn);
    }
  }


  WlzWriteObj(stdout, newObj);

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
