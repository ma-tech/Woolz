#pragma ident "%W%(%G%) - richard"
/************************************************************************
*   Copyright  :   1994 Medical Research Council, UK.                   *
*                  All rights reserved.                                 *
*************************************************************************
*   Address    :   MRC Human Genetics Unit,                             *
*                  Western General Hospital,                            *
*                  Edinburgh, EH4 2XU, UK.                              *
*************************************************************************
*   Project    :   Woolz Library					*
*   File       :   WlzPatchObjRegister.c				*
*************************************************************************
* This module has been copied from the original woolz library and       *
* modified for the public domain distribution. The original authors of  *
* the code and the original file headers and comments are in the        *
* HISTORY file.                                                         *
*************************************************************************
*   Author Name :  Richard Baldock					*
*   Author Login:  richard@hgu.mrc.ac.uk				*
*   Date        :  Fri Feb  6 14:42:26 1998				*
*   Synopsis    : 							*
*************************************************************************
*   Maintenance :  date - name - comments (Last changes at the top)	*
************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include <Wlz.h>
#include <Reconstruct.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static WlzIVertex2	defMaxShift;

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s [-h] [-v] [-t#[,#]] [-g] [-G] [<input file>]\n"
	  "\tWoolz in a compound woolz object assumed\n"
	  "\tto be output of patch images from xmgrab.\n"
	  "\tRegister the patches and output the single\n"
	  "\tdomain object.\n"
	  "\tOptions are:\n"
	  "\t  -t#[,#]   Set maximum translation values\n"
	  "\t            default (30,30) one input value implies\n"
	  "\t            equal maximum shift in x & y directions\n"
	  "\t  -g        Attempt to reset mean grey values to be equal\n"
	  "\t            within matched regions (on by default)\n"
	  "\t  -G        Switch off grey-level matching\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -v        Verbose operation\n"
	  "",
	  proc_str);
  return;
}

int WlzIsAdjacentPatch(
  WlzObject	*obj1,
  WlzObject	*obj2)
{
  int	adjacent=0;

  /* patches are adjacent if either the line bounds  or the column bounds
     are identical and there is some overlap */

  /* test on the columns */
  if((obj1->domain.i->kol1 == obj2->domain.i->kol1) &&
     (obj1->domain.i->lastkl == obj2->domain.i->lastkl)){
    if((obj1->domain.i->lastln < obj2->domain.i->line1) ||
       (obj1->domain.i->line1 > obj2->domain.i->lastln)){
      adjacent = 0;
    }
    else {
      adjacent = 1;
    }
  }
  else if((obj1->domain.i->line1 == obj2->domain.i->line1) &&
	  (obj1->domain.i->lastln == obj2->domain.i->lastln)){
    if((obj1->domain.i->lastkl < obj2->domain.i->kol1) ||
       (obj1->domain.i->kol1 > obj2->domain.i->lastkl)){
      adjacent = 0;
    }
    else {
      adjacent = 1;
    }
  }

  return adjacent;
}

typedef struct _patchTree{
  WlzObject	*obj;
  int		index;
  int		offsetsCalculatedFlag;
  int		depth;
  double	xOff;
  double	yOff;
  struct _patchTree	*children[4];
  int		nchildren;
  double	cost;
} WlzPatchTree;

int sortPatch(
  const void *p1,
  const void *p2)
{
  WlzPatchTree	*patch_1=*((WlzPatchTree **) p1);
  WlzPatchTree	*patch_2=*((WlzPatchTree **) p2);

  if( patch_1->cost < patch_2->cost ){
    return -1;
  }
  if( patch_1->cost > patch_2->cost ){
    return 1;
  }
  return 0;
}


WlzPatchTree *WlzMakePatchTree(
  WlzObject	*obj,
  int		depth,
  double	cost)
{
  WlzPatchTree	*patchTree;

  /* allocate space and initialise */
  patchTree = (WlzPatchTree *) AlcCalloc(1, sizeof(WlzPatchTree));
  patchTree->obj = WlzAssignObject(obj, NULL);
  patchTree->depth = depth;
  patchTree->cost = cost;

  return patchTree;
}

WlzErrorNum	WlzGetPatchTreeToDepth(
  WlzObject	**objs,
  int		nobjs,
  WlzPatchTree	*patchTree,
  int		depth)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  int		i;
  WlzGreyType	dstGType;
  double	dstMin, dstMax, dstSum, dstSumSq;
  double	dstMean, dstStdDev;

  if((patchTree == NULL) ||
     (patchTree->depth >= depth)){
    return errNum;
  }

  if( patchTree->depth == (depth - 1) ){
    for(i=0; (i < nobjs) && (patchTree->nchildren < 4); i++){
      if( objs[i] && WlzIsAdjacentPatch(patchTree->obj, objs[i]) ){
	WlzGreyStats(objs[i], &dstGType, &dstMin, &dstMax,
		     &dstSum, &dstSumSq, &dstMean, &dstStdDev,
		     NULL);
	patchTree->children[patchTree->nchildren] =
	  WlzMakePatchTree(objs[i], depth, dstMean);
	patchTree->children[patchTree->nchildren]->index = i;
	WlzFreeObj(objs[i]);
	objs[i] = NULL;
	patchTree->nchildren++;
      }
    }
    /* sort the children */
    qsort(patchTree->children, patchTree->nchildren,
	  sizeof(WlzPatchTree *), sortPatch);
  }
  else {
    for(i=0; i < patchTree->nchildren; i++){
      errNum = WlzGetPatchTreeToDepth(objs, nobjs,
				      patchTree->children[i], depth);
    }
  }

  return errNum;
}

static int numObjsLeft(
  WlzObject	**objs, 
  int		nobjs)
{
  int i, n;

  for(i=0, n=0; i < nobjs; i++){
    if( objs[i] ){
      n++;
    }
  }
  return n;
}

WlzPatchTree *WlzGetPatchTree(
  WlzObject	*obj,
  WlzObject	**objs,
  int		nobjs)
{
  int		i;
  int		nPatches;
  WlzPatchTree	*patchTree;
  WlzObject	*adjObj;

  /* allocate space and initialise */
  patchTree = (WlzPatchTree *) AlcCalloc(1, sizeof(WlzPatchTree));
  patchTree->obj = WlzAssignObject(obj, NULL);

  /* check for non-NULL adjacent objects  - depth first search */
  for(i=0, nPatches=0; (i < nobjs) && (nPatches < 4); i++){
    if( objs[i] && WlzIsAdjacentPatch(obj, objs[i]) ){
      adjObj = objs[i];
      objs[i] = NULL;
      patchTree->children[nPatches] = WlzGetPatchTree(adjObj, objs, nobjs);
      WlzFreeObj(adjObj);
      nPatches++;
    }
  }

  patchTree->offsetsCalculatedFlag = 0;
  patchTree->nchildren = nPatches;

  return patchTree;
}

WlzErrorNum WlzFreePatchTree(
  WlzPatchTree	*patchTree)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  int		i;

  /* free sub trees */
  for(i=0; (i < patchTree->nchildren) && (errNum == WLZ_ERR_NONE); i++){
    errNum = WlzFreePatchTree(patchTree->children[i]);
  }

  /* free the object */
  if( errNum == WLZ_ERR_NONE ){
    errNum = WlzFreeObj(patchTree->obj);
  }

  /* free this branch */
  AlcFree((void *) patchTree);

  return errNum;
}

double WlzMass(
  WlzObject *obj,
  WlzErrorNum *dstErr)
{
  double		mass;
  WlzIntervalWSpace 	iwsp;
  WlzGreyWSpace		gwsp;
  WlzGreyP		gptr;
  int			i;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check object */
  if( obj == (WlzObject *) NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
    mass = -1;
  }

  if( errNum == WLZ_ERR_NONE ){
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
	mass = -1;
      }
      break;
    
    case WLZ_EMPTY_OBJ:
      mass = 0;
      break;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      mass = -1;
      break;

    }
  }

  /* calculate mass */
  if( errNum == WLZ_ERR_NONE ){
    mass = 0.0;
    errNum = WlzInitGreyScan(obj, &iwsp, &gwsp);
    if( errNum == WLZ_ERR_NONE ){
      while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){
	gptr = gwsp.u_grintptr;
	switch (gwsp.pixeltype) {

	case WLZ_GREY_INT:
	  for (i=0; i<iwsp.colrmn; i++, gptr.inp++){
	    if( *gptr.inp < 0 ){
	      mass -= *gptr.inp;
	    }
	    else {
	      mass += *gptr.inp;
	    }
	  }
	  break;

	case WLZ_GREY_SHORT:
	  for (i=0; i<iwsp.colrmn; i++, gptr.shp++){
	    if( *gptr.shp < 0 ){
	      mass -= *gptr.shp;
	    }
	    else {
	      mass += *gptr.shp;
	    }
	  }
	  break;

	case WLZ_GREY_UBYTE:
	  for (i=0; i<iwsp.colrmn; i++, gptr.ubp++){
	    if( *gptr.ubp < 0 ){
	      mass -= *gptr.ubp;
	    }
	    else {
	      mass += *gptr.ubp;
	    }
	  }
	  break;

	case WLZ_GREY_FLOAT:
	  for (i=0; i<iwsp.colrmn; i++, gptr.flp++){
	    if( *gptr.flp < 0 ){
	      mass -= *gptr.flp;
	    }
	    else {
	      mass += *gptr.flp;
	    }
	  }
	  break;

	case WLZ_GREY_DOUBLE:
	  for (i=0; i<iwsp.colrmn; i++, gptr.dbp++){
	    if( *gptr.dbp < 0 ){
	      mass -= *gptr.dbp;
	    }
	    else {
	      mass += *gptr.dbp;
	    }
	  }
	  break;

	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
	}
      }
      if(errNum == WLZ_ERR_EOO)	  /* Reset error from end of intervals */ 
      {
	errNum = WLZ_ERR_NONE;
      }
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return mass;
}

WlzObject *WlzGreyShift(
  WlzObject	*obj,
  double	delta,
  WlzErrorNum	*dstErr)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzObject	*rtnObj=NULL;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzGreyP		gptr;
  int			i, idelta;
  double		r;

  /* check the object */
  /* 2D only for now to test the idea */
  if( obj ){
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
      if( obj->domain.core ){
	if(obj->domain.core->type == WLZ_EMPTY_DOMAIN){
	  rtnObj = WlzMakeEmpty(&errNum);
	}
	else {
	  if(!obj->values.core ){
	    errNum = WLZ_ERR_VALUES_NULL;
	  }
	}
      }
      else {
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      break;

    case WLZ_EMPTY_OBJ:
      rtnObj = WlzMakeEmpty(&errNum);
      break;

    case WLZ_3D_DOMAINOBJ:
    case WLZ_TRANS_OBJ:
    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }
  else {
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  /* 2D case */
  if( (errNum == WLZ_ERR_NONE) && (rtnObj == NULL) ){
    if( rtnObj = WlzCopyObject(obj, &errNum) ){
      if( (errNum = WlzInitGreyScan(rtnObj, &iwsp, &gwsp)) == WLZ_ERR_NONE ){
	while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){
	  gptr = gwsp.u_grintptr;
	  switch (gwsp.pixeltype) {
	  case WLZ_GREY_INT:
	    for(i=0; i<iwsp.colrmn; i++, gptr.inp++){
	      if( delta > 0.0 ){
		idelta = (int) (delta + ((double) (random()&0xffff))/0xffff);
	      }
	      else if( delta < 0.0 ){
		idelta = (int) (delta - ((double) (random()&0xffff))/0xffff);
	      }
	      else {
		idelta = 0;
	      }
	      *gptr.inp += idelta;
	    }
	    break;

	  case WLZ_GREY_SHORT:
	    for(i=0; i<iwsp.colrmn; i++, gptr.shp++){
	      if( delta > 0.0 ){
		idelta = (int) (delta + ((double) (random()&0xffff))/0xffff);
	      }
	      else if( delta < 0.0 ){
		idelta = (int) (delta - ((double) (random()&0xffff))/0xffff);
	      }
	      else {
		idelta = 0;
	      }
	      *gptr.shp += idelta;
	    }
	    break;

	  case WLZ_GREY_UBYTE:
	    for(i=0; i<iwsp.colrmn; i++, gptr.ubp++){
	      if( delta > 0.0 ){
		idelta = (int) (delta + ((double) (random()&0xffff))/0xffff);
	      }
	      else if( delta < 0.0 ){
		idelta = (int) (delta - ((double) (random()&0xffff))/0xffff);
	      }
	      else {
		idelta = 0;
	      }
	      idelta += *gptr.ubp;
	      *gptr.ubp = WLZ_CLAMP(idelta, 0, 255);
	    }
	    break;

	  case WLZ_GREY_FLOAT:
	    for(i=0; i<iwsp.colrmn; i++, gptr.flp++){
	      *gptr.flp += delta;
	    }
	    break;

	  case WLZ_GREY_DOUBLE:
	    for(i=0; i<iwsp.colrmn; i++, gptr.dbp++){
	      *gptr.dbp += delta;
	    }
	    break;
	  }
	}
	if( errNum == WLZ_ERR_EOO ){
	  errNum = WLZ_ERR_NONE;
	}
	if( errNum != WLZ_ERR_NONE ){
	  WlzFreeObj(rtnObj);
	  rtnObj = NULL;
	}
      }
      else {
	WlzFreeObj(rtnObj);
	rtnObj = NULL;
      }
    }
  }

  /* set return error */
  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

WlzObject *WlzGreyScale(
  WlzObject	*obj,
  double	scale,
  WlzErrorNum	*dstErr)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzObject	*rtnObj=NULL;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzGreyP		gptr;
  int			i, idelta;
  double		r;

  /* check the object */
  /* 2D only for now to test the idea */
  if( obj ){
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
      if( obj->domain.core ){
	if(obj->domain.core->type == WLZ_EMPTY_DOMAIN){
	  rtnObj = WlzMakeEmpty(&errNum);
	}
	else {
	  if(!obj->values.core ){
	    errNum = WLZ_ERR_VALUES_NULL;
	  }
	}
      }
      else {
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      break;

    case WLZ_EMPTY_OBJ:
      rtnObj = WlzMakeEmpty(&errNum);
      break;

    case WLZ_3D_DOMAINOBJ:
    case WLZ_TRANS_OBJ:
    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }
  else {
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  /* 2D case */
  if( (errNum == WLZ_ERR_NONE) && (rtnObj == NULL) ){
    if( rtnObj = WlzCopyObject(obj, &errNum) ){
      if((scale != 1.0) && 
	 (errNum = WlzInitGreyScan(rtnObj, &iwsp, &gwsp)) == WLZ_ERR_NONE ){
	while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){
	  gptr = gwsp.u_grintptr;
	  switch (gwsp.pixeltype) {
	  case WLZ_GREY_INT:
	    for(i=0; i<iwsp.colrmn; i++, gptr.inp++){
	      idelta = (int) (scale * (*gptr.inp)
			      + ((double) (random()&0xffff))/0xffff);
	      *gptr.inp = idelta;
	    }
	    break;

	  case WLZ_GREY_SHORT:
	    for(i=0; i<iwsp.colrmn; i++, gptr.shp++){
	      idelta = (int) (scale * (*gptr.shp)
			      + ((double) (random()&0xffff))/0xffff);
	      *gptr.shp = idelta;
	    }
	    break;

	  case WLZ_GREY_UBYTE:
	    for(i=0; i<iwsp.colrmn; i++, gptr.ubp++){
	      idelta = (int) (scale * (*gptr.ubp)
			      + ((double) (random()&0xffff))/0xffff);
	      *gptr.ubp = WLZ_CLAMP(idelta, 0, 255);
	    }
	    break;

	  case WLZ_GREY_FLOAT:
	    for(i=0; i<iwsp.colrmn; i++, gptr.flp++){
	      *gptr.flp *= scale;
	    }
	    break;

	  case WLZ_GREY_DOUBLE:
	    for(i=0; i<iwsp.colrmn; i++, gptr.dbp++){
	      *gptr.dbp += scale;
	    }
	    break;
	  }
	}
	if( errNum == WLZ_ERR_EOO ){
	  errNum = WLZ_ERR_NONE;
	}
	if( errNum != WLZ_ERR_NONE ){
	  WlzFreeObj(rtnObj);
	  rtnObj = NULL;
	}
      }
      else {
	WlzFreeObj(rtnObj);
	rtnObj = NULL;
      }
    }
  }

  /* set return error */
  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

double WlzGreyMeanDifference(
  WlzObject	*obj1,
  WlzObject	*obj2,
  double	samplePercent,
  WlzErrorNum	*dstErr)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  double	rtnDiff=10000.0;
  int		idiff, idiffIncr;
  WlzObject	*o1, *o2;
  int		countedArea;
  WlzIntervalWSpace 	iwsp1, iwsp2;
  WlzGreyWSpace		gwsp1, gwsp2;
  WlzGreyP		gptr1, gptr2;
  int			i;
  int			pixelIncr;
  
  /* get intersection */
  o1 = WlzIntersect2(obj1, obj2, NULL);
  o2 = WlzMakeMain(o1->type, o1->domain, o1->values, NULL, NULL, NULL);
  o1->values = WlzAssignValues(obj1->values, NULL);
  o2->values = WlzAssignValues(obj2->values, NULL);

  /* loop through grey tables checking differences */
  WlzInitGreyScan(o1, &iwsp1, &gwsp1);
  WlzInitGreyScan(o2, &iwsp2, &gwsp2);
  idiff = 0;
  countedArea = 0;
  pixelIncr = (int) (100.0 / samplePercent);
  while( (errNum = WlzNextGreyInterval(&iwsp1)) == WLZ_ERR_NONE ){
    (void) WlzNextGreyInterval(&iwsp2);
    gptr1 = gwsp1.u_grintptr;
    gptr2 = gwsp2.u_grintptr;
    switch (gwsp1.pixeltype) {
    case WLZ_GREY_INT:
      for (i=0; i<iwsp1.colrmn; i += pixelIncr,
	     gptr1.inp += pixelIncr, gptr2.inp += pixelIncr){
	idiffIncr = *gptr1.inp - *gptr2.inp;
	idiff += (idiffIncr < 0)?-idiffIncr:idiffIncr;
	countedArea++;
      }
      break;
    case WLZ_GREY_SHORT:
      for (i=0; i<iwsp1.colrmn; i += pixelIncr,
	     gptr1.shp += pixelIncr, gptr2.shp += pixelIncr){
	idiffIncr = *gptr1.shp - *gptr2.shp;
	idiff += (idiffIncr < 0)?-idiffIncr:idiffIncr;
	countedArea++;
      }
      break;
    case WLZ_GREY_UBYTE:
      for (i=0; i<iwsp1.colrmn; i += pixelIncr,
	     gptr1.ubp += pixelIncr, gptr2.ubp += pixelIncr){
	idiffIncr = *gptr1.ubp - *gptr2.ubp;
	idiff += (idiffIncr < 0)?-idiffIncr:idiffIncr;
	countedArea++;
      }
      break;
    case WLZ_GREY_FLOAT:
      for (i=0; i<iwsp1.colrmn; i += pixelIncr,
	     gptr1.flp += pixelIncr, gptr2.flp += pixelIncr){
	idiffIncr = *gptr1.flp - *gptr2.flp;
	idiff += (idiffIncr < 0)?-idiffIncr:idiffIncr;
	countedArea++;
      }
      break;
    case WLZ_GREY_DOUBLE:
      for (i=0; i<iwsp1.colrmn; i += pixelIncr,
	     gptr1.dbp += pixelIncr, gptr2.dbp += pixelIncr){
	idiffIncr = *gptr1.dbp - *gptr2.dbp;
	idiff += (idiffIncr < 0)?-idiffIncr:idiffIncr;
	countedArea++;
      }
      break;
    }
  } 
 
  rtnDiff = (double) idiff / countedArea;
  WlzFreeObj(o1);
  WlzFreeObj(o2);

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnDiff;
}

WlzErrorNum DumbRegMatch(
  WlzDVertex2	*shift,
  double	*matchVal,
  WlzObject	*obj1,
  WlzObject	*obj2,
  WlzIVertex2	maxShift)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzObject	*o2;
  double       	meanDiff, minMeanDiff;
  int		i, j;
  WlzIVertex2	rtnShift;

  /* no checks just go for it */

  minMeanDiff = 10000.0;
  rtnShift.vtX = 0;
  rtnShift.vtY = 0;
  for(i = - (int) maxShift.vtX; i <= (int) maxShift.vtX; i++){
    for(j = - (int) maxShift.vtY; j <= (int) maxShift.vtY; j++){
      o2 = WlzShiftObject(obj2, i, j, 0, NULL);
      meanDiff = WlzGreyMeanDifference(obj1, o2, 5.0, NULL);
      if( meanDiff < minMeanDiff ){
	minMeanDiff = meanDiff;
	rtnShift.vtX = i;
	rtnShift.vtY = j;
      }
      WlzFreeObj(o2);
    }
  }

  shift->vtX = rtnShift.vtX;
  shift->vtY = rtnShift.vtY;
  *matchVal = minMeanDiff;
  return errNum;
}

int WlzRegisterPatchTree(
  WlzPatchTree	*patchTree)
{
  WlzObject	*obj, *obj1, *obj2;
  int 		i;
  RecError	recErr=REC_ERR_NONE;
  WlzDVertex2	shift;
  double	ccVal;
  WlzIVertex2	maxShift;
  RecPPControl	ppCtrl;

  /* calculate the offsets for each child then
     register each child tree */
  maxShift.vtX = defMaxShift.vtX;
  maxShift.vtY = defMaxShift.vtY;
  ppCtrl.method = REC_PP_WINDOW;
  ppCtrl.window.function = WLZ_WINDOWFN_RECTANGLE;
  ppCtrl.window.size.vtX = 90;
  ppCtrl.window.size.vtY = 90;
  ppCtrl.window.offset.vtX = 0;
  ppCtrl.window.offset.vtY = 0;
  ppCtrl.sample.function = NULL;
  ppCtrl.sample.factor = 0;
  ppCtrl.erode = 0;
  for(i=0; i < patchTree->nchildren; i++){
    WlzStandardIntervalDomain(patchTree->obj->domain.i);
    WlzStandardIntervalDomain(patchTree->children[i]->obj->domain.i);

    (void) DumbRegMatch(&shift, &ccVal, patchTree->obj,
			patchTree->children[i]->obj, maxShift);
/*    recErr = RecTranMatch(&shift, &ccVal, patchTree->obj,
			  patchTree->children[i]->obj, maxShift, &ppCtrl);*/
    patchTree->children[i]->offsetsCalculatedFlag = 1;
    patchTree->children[i]->xOff = shift.vtX;
    patchTree->children[i]->yOff = shift.vtY;
    WlzRegisterPatchTree(patchTree->children[i]);
  }

  return 0;
}

WlzObject *WlzPatchTreeToObject(
  WlzPatchTree	*patchTree,
  int		alignGreysFlg)
{
  WlzObject	*rtnObj;
  WlzObject	*obj1, *objs[2];
  int		i;

  /* get the patch children adding as required */
  objs[0] = WlzAssignObject(WlzMakeMain(patchTree->obj->type,
					patchTree->obj->domain,
					patchTree->obj->values,
					NULL, NULL, NULL), NULL);
  for(i=0; i < patchTree->nchildren; i++){
    /* get object so far and remove overlap with the current object */
    obj1 = WlzPatchTreeToObject(patchTree->children[i], alignGreysFlg);
    objs[1] = WlzAssignObject(WlzDiffDomain(obj1, objs[0], NULL), NULL);

    /* reset the grey values if necessary  - approximate as a shift */
    if( alignGreysFlg ){
      WlzObject 	*obj2;
      double 		min1, max1, sum1, sumSq1, mean1, stdDev1;
      double 		min2, max2, sum2, sumSq2, mean2, stdDev2;
      WlzGreyType	gType;

      if( obj2 = WlzIntersect2(obj1, objs[0], NULL) ){
	obj2->values.core = obj1->values.core;
	WlzGreyStats(obj2, &gType, &min1, &max1,
		     &sum1, &sumSq1, &mean1, &stdDev1, NULL);
	obj2->values.core = objs[0]->values.core;
	WlzGreyStats(obj2, &gType, &min2, &max2,
		     &sum2, &sumSq2, &mean2, &stdDev2, NULL);
	obj2->values.core = NULL;
	WlzFreeObj(obj2);
	obj2 = WlzAssignObject(WlzGreyScale(objs[0], mean1 / mean2, NULL), NULL);
	WlzFreeObj(objs[0]);
	objs[0] = obj2;
      }
      
    }

    WlzFreeObj(obj1);
    obj1 = WlzUnionN(2, objs, 1, NULL);
    WlzStandardIntervalDomain(obj1->domain.i);
    WlzFreeObj(objs[0]);
    WlzFreeObj(objs[1]);
    objs[0] = obj1;
  }

  /* now shift the given object */
  rtnObj = WlzShiftObject(objs[0], (int) patchTree->xOff,
			  (int) patchTree->yOff, 0, NULL);
  WlzFreeObj(objs[0]);

  return rtnObj;
}

WlzErrorNum WlzPatchFacts(
  WlzPatchTree	*patchTree,
  FILE		*fp,
  char		**dstStr,
  int 		verbose)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  int		i;

  /* print offsets */
  fprintf(fp,
	  "WlzPatchFacts: offsets\n"
	  "==============\n"
	  " object index = %d, offsetsCalculatedFlag = %d\n"
	  " x-offset = %f, y-offset = %f\n",
	  patchTree->index,
	  patchTree->offsetsCalculatedFlag,
	  patchTree->xOff, patchTree->yOff);
  
  /* print object facts */
  if( verbose ){
    fprintf(fp,
	    "WlzPatchFacts: object detail\n"
	    "==============\n");
    WlzObjectFacts(patchTree->obj, fp, dstStr, 0);
    fprintf(fp,
	    "===================================================\n");
  }

  /* print children facts */
  for(i=0; i < patchTree->nchildren; i++){
    WlzPatchFacts(patchTree->children[i], fp, dstStr, verbose);
  }

  return errNum;
}

 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *nobj, **objs, *obj1;
  WlzPatchTree	*patchTree;
  WlzCompoundArray	*cobj;
  FILE		*inFile;
  char 		optList[] = "gGhvt:";
  int		option;
  int		verboseFlg=0;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  int		i, startIndx, depth;
  int		alignGreysFlg=1;
  WlzGreyType	dstGType;
  double	dstMin, dstMax, dstSum, dstSumSq;
  double	dstMean, dstStdDev, extMean;

  /* set the default maximum shift */
  defMaxShift.vtX = 30;
  defMaxShift.vtY = 30;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){
    case 'g':
      alignGreysFlg = 1;
      break;

    case 'G':
      alignGreysFlg = 0;
      break;

    case 't':
      switch( sscanf(optarg, "%d,%d",
		     &(defMaxShift.vtX), &(defMaxShift.vtY)) ){

      default:
      case 0:
	fprintf(stderr, "%s: no search width set\n", argv[0]);
        usage(argv[0]);
        return 1;

      case 1:
	defMaxShift.vtY = defMaxShift.vtX;
	break;

      case 2:
	break;

      }
      break;

    case 'v':
      verboseFlg = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return WLZ_ERR_UNSPECIFIED;
    }
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return WLZ_ERR_UNSPECIFIED;
    }
  }

  /* read objects and threshold if possible */
  while( (obj = WlzAssignObject(WlzReadObj(inFile, &errNum), NULL)) != NULL) 
  {
    switch( obj->type )
    {
    case WLZ_COMPOUND_ARR_1:
    case WLZ_COMPOUND_ARR_2:
      cobj = (WlzCompoundArray *) obj;
      if( cobj->n < 1 ){
	nobj = WlzMakeEmpty(&errNum);
      }
      else if( cobj->n == 1 ){
	nobj = WlzAssignObject(cobj->o[0], &errNum);
      }
      else {
	/* copy the object list so the patch tree can be built */
	objs = (WlzObject **) AlcMalloc(sizeof(WlzObject *) * cobj->n);
	for(i=0; i < cobj->n; i++){
	  objs[i] = WlzAssignObject(cobj->o[i], NULL);
	}

	/* find the most dense object - assumed to have the
	   bulk of the foreground therefore a suitable starting
	   point for matching. */
	startIndx = 0;
	WlzGreyStats(objs[0], &dstGType, &dstMin, &dstMax,
		     &dstSum, &dstSumSq, &extMean, &dstStdDev,
		     NULL);
	for(i=1; i < cobj->n; i++){
	  WlzGreyStats(objs[i], &dstGType, &dstMin, &dstMax,
		       &dstSum, &dstSumSq, &dstMean, &dstStdDev,
		       NULL);
	  if( dstMean < extMean ){
	    extMean = dstMean;
	    startIndx = i;
	  }
	}

	/* generate the patch tree - depth first search */
	/*obj1 = objs[startIndx];
	objs[startIndx] = NULL;
	patchTree = WlzGetPatchTree(obj1, objs, cobj->n);
	WlzFreeObj(obj1);*/

	/* generate the patch tree - breadth first search */
	depth = 0;
	patchTree = WlzMakePatchTree(objs[startIndx], depth, extMean);
	patchTree->index = startIndx;
	objs[startIndx] = NULL;
	while( numObjsLeft(objs, cobj->n) ){
	  depth++;
	  WlzGetPatchTreeToDepth(objs, cobj->n, patchTree, depth);
	}

	/* calculate the offsets */
	WlzRegisterPatchTree(patchTree);

	/* print out the patch data and offsets */
	if( verboseFlg ){
	  WlzPatchFacts(patchTree, stderr, NULL, 0);
	}

	/* get the patch object */
	nobj = WlzPatchTreeToObject(patchTree, alignGreysFlg);
	(void) WlzFreePatchTree(patchTree);
	(void) AlcFree((void *) objs);
      }
      if( errNum == WLZ_ERR_NONE ){
	errNum = WlzWriteObj(stdout, nobj);
      }
      if( nobj ){
	WlzFreeObj(nobj);
      }
      break;
      
    default:
      errNum = WlzWriteObj(stdout, obj);
      break;
    }

    WlzFreeObj(obj);
  }

  /* trap the WLZ_ERR_READ_EOF since this is a legal way of indicating
     the end of objects in a file */
  if( errNum = WLZ_ERR_READ_EOF ){
    errNum = WLZ_ERR_NONE;
  }

  return errNum;
}
