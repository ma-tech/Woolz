#pragma ident "MRC HGU $Id$"
/************************************************************************
*   Copyright  :   1994 Medical Research Council, UK.                   *
*                  All rights reserved.                                 *
*************************************************************************
*   Address    :   MRC Human Genetics Unit,                             *
*                  Western General Hospital,                            *
*                  Edinburgh, EH4 2XU, UK.                              *
*************************************************************************
*   Project    :   MRC HGU Image Processing Utilities			*
*   File       :   WlzFixDomains.c					*
*************************************************************************
*   Author Name :  Richard Baldock					*
*   Author Login:  richard@hgu.mrc.ac.uk				*
*   Date        :  Wed Jan 16 15:17:57 2002				*
*   $Revision$								*
*   $Name$								*
*   Synopsis    : 							*
*************************************************************************
*   Maintenance :  date - name - comments (Last changes at the top)	*
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#include <Wlz.h>
#include <HGUDlpList.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

extern WlzErrorNum WlzDomainSizeSelect(
  WlzObject		*obj,
  int			size,
  int			meshFlg,
  WlzConnectType	connectivity,
  int			*dstArraySizeObjs,
  WlzObject		***dstArrayObjs);

static WlzErrorNum WlzDomainSizeSelect3d(
  WlzObject		*obj,
  int			size,
  int			meshFlg,
  WlzConnectType	connectivity,
  int			*numObjs,
  WlzObject		***objs)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzObject	**rtnObjs;
  int		rtnNobjs, maxRtnNobjs;
  WlzObject	**objArray, *tmpObj, *obj1;
  int		i, n, p, nobjs, maxNobjs=1024;
  WlzValues	values;
  WlzDomain	domain;
  WlzPlaneDomain	*pdom;
  
  /* check parameters */
  if( size <= 0 ){
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else {
    switch( connectivity ){
    case WLZ_4_CONNECTED:
    case WLZ_8_CONNECTED:
      break;

    case WLZ_6_CONNECTED:
    case WLZ_18_CONNECTED:
    case WLZ_26_CONNECTED:
      break;

    default:
      errNum = WLZ_ERR_PARAM_DATA;
      break;
    }
  }

  /* check object */
  if( obj->domain.core == NULL ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if( obj->domain.p->type != WLZ_PLANEDOMAIN_DOMAIN ){
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }  

  /* allocate space and start running through the planes */
  if( errNum == WLZ_ERR_NONE ){
    maxRtnNobjs = 1024;
    rtnObjs = (WlzObject **) AlcCalloc(sizeof(WlzObject *), maxRtnNobjs);
    rtnNobjs = 0;
    values.core = 0;
    pdom = obj->domain.p;
    for(p=pdom->plane1; p < pdom->lastpl; p++){
      if( pdom->domains[p-pdom->plane1].core ){
	tmpObj = WlzMakeMain(WLZ_2D_DOMAINOBJ,
			     pdom->domains[p-pdom->plane1],
			     values, NULL, NULL, &errNum);
	tmpObj = WlzAssignObject(tmpObj, &errNum);
	if((errNum =
	    WlzDomainSizeSelect(tmpObj, size, meshFlg,
				connectivity,
				&nobjs, &objArray)) == WLZ_ERR_NONE ){

	  if( (rtnNobjs + nobjs) > maxRtnNobjs ){
	    maxRtnNobjs += (nobjs > 1024) ? nobjs : 1024;
	    rtnObjs = (WlzObject **)
	      AlcRealloc(rtnObjs, sizeof(WlzObject *) * maxRtnNobjs);
	  }

	  for(i=0; i < nobjs; i++){
	    domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN, p, p,
					  objArray[i]->domain.i->line1,
					  objArray[i]->domain.i->lastln,
					  objArray[i]->domain.i->kol1,
					  objArray[i]->domain.i->lastkl,
					  &errNum);
	    domain.p->domains[0] = WlzAssignDomain(objArray[i]->domain,
						   NULL);
	    obj1 = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values, NULL, NULL,
			       &errNum);
	    rtnObjs[rtnNobjs++] = WlzAssignObject(obj1, &errNum);
	    WlzFreeObj(objArray[i]);
	  }
	  AlcFree(objArray);
	}
	WlzFreeObj(tmpObj);
      }
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    *objs = rtnObjs;
    *numObjs = rtnNobjs;
  }
  else {
    *objs = NULL;
    *numObjs = 0;
  }
  return errNum;
}

WlzErrorNum WlzDomainSizeSelect(
  WlzObject		*obj,
  int			size,
  int			meshFlg,
  WlzConnectType	connectivity,
  int			*numObjs,
  WlzObject		***objs)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzObject	**objArray, *tmpObj;
  int		i, n, nobjs, maxNobjs=1024;
  WlzValues	values;

  /* check parameters */
  if( size <= 0 ){
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else {
    switch( connectivity ){
    case WLZ_4_CONNECTED:
    case WLZ_8_CONNECTED:
      break;

    case WLZ_6_CONNECTED:
    case WLZ_18_CONNECTED:
    case WLZ_26_CONNECTED:
      break;

    default:
      errNum = WLZ_ERR_PARAM_DATA;
      break;
    }
  }

  /* check object */
  if( errNum == WLZ_ERR_NONE ){
    if( obj == NULL ){
      errNum = WLZ_ERR_OBJECT_NULL;
    }
    else {
      switch( obj->type ){
      case WLZ_2D_DOMAINOBJ:
	if( obj->domain.core == NULL ){
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	break;

      case WLZ_3D_DOMAINOBJ:
	return WlzDomainSizeSelect3d(obj, size, meshFlg, connectivity,
				     numObjs, objs);

      case WLZ_EMPTY_OBJ:
	*objs = NULL;
	numObjs = 0;
	return WLZ_ERR_NONE;

      case WLZ_TRANS_OBJ:
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
      }
    }
  }

  /* 2D object, segment and select */
  if( errNum == WLZ_ERR_NONE ){
    if( (errNum = WlzLabel(obj, &nobjs, &objArray, maxNobjs, 0,
			   connectivity)) == WLZ_ERR_NONE ){
      for(i=0; i < nobjs; i++){
	if( meshFlg ){
	  if( WlzArea(objArray[i], &errNum) > size ){
	    WlzFreeObj(objArray[i]);
	    objArray[i] = NULL;
	    
	  }
	}
	else {
	  if( WlzArea(objArray[i], &errNum) <= size ){
	    WlzFreeObj(objArray[i]);
	    objArray[i] = NULL;
	  }
	}
      }

      /* collect non-NULL objects */
      for(i=0, n=0; i < nobjs; i++){
	if( objArray[i] ){
	  tmpObj = objArray[i];
	  objArray[i] = NULL;
	  objArray[n] = tmpObj;
	  n++;
	}
      }
      nobjs = n;
    }
  }
	

  if( errNum == WLZ_ERR_NONE ){
    *objs = objArray;
    *numObjs = nobjs;
  }
  else {
    *objs = NULL;
    *numObjs = 0;
  }
  return errNum;
}

extern WlzObject *WlzDomainFillHoles(
  WlzObject		*obj,
  int			size,
  WlzConnectType	connectivity,
  WlzErrorNum		*dstErr);

static WlzObject *WlzDomainFillHoles3d(
  WlzObject		*obj,
  int			size,
  WlzConnectType	connectivity,
  WlzErrorNum		*dstErr)
{
  WlzObject	*rtnObj=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check the planedomain - object and object type already checked */
  if( obj->domain.core == NULL ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if( obj->domain.p->type != WLZ_PLANEDOMAIN_DOMAIN ){
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }

  /* create a new 3D object and set domains */
  if( errNum == WLZ_ERR_NONE ){
    WlzValues	values;
    WlzDomain	domain;
    WlzObject	*obj1, *obj2;
    int		p;
    
    if( domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
				      obj->domain.p->plane1,
				      obj->domain.p->lastpl,
				      obj->domain.p->line1,
				      obj->domain.p->lastln,
				      obj->domain.p->kol1,
				      obj->domain.p->lastkl,
				      &errNum) ){
      values.core = NULL;
      rtnObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values, NULL, NULL,
			   &errNum);
      for(p=domain.p->plane1;
	  (errNum == WLZ_ERR_NONE) &&(p <= domain.p->lastpl);
	  p++){
	if( obj->domain.p->domains[p-domain.p->plane1].core ){
	  obj1 = WlzMakeMain(WLZ_2D_DOMAINOBJ,
			     obj->domain.p->domains[p-domain.p->plane1],
			     values, NULL, NULL, &errNum);
	  obj2 = WlzDomainFillHoles(obj1, size, connectivity, &errNum);
	  domain.p->domains[p-domain.p->plane1] =
	    WlzAssignDomain(obj2->domain, &errNum);
	  WlzFreeObj(obj1);
	  WlzFreeObj(obj2);
	}
      }
    }
    if( errNum != WLZ_ERR_NONE ){
      WlzFreeObj(rtnObj);
      rtnObj = NULL;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

WlzObject *WlzDomainFillHoles(
  WlzObject		*obj,
  int			size,
  WlzConnectType	connectivity,
  WlzErrorNum		*dstErr)
{
  WlzObject	*rtnObj=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzValues	values;
  WlzDomain	domain;
  WlzObject	*obj1, *backObj, *holeObj, **objs;
  int		i, numObjs;

  /* check parameters */
  if( size <= 0 ){
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else {
    switch( connectivity ){
    case WLZ_4_CONNECTED:
      break;

    case WLZ_8_CONNECTED:
    case WLZ_6_CONNECTED:
    case WLZ_18_CONNECTED:
    case WLZ_26_CONNECTED:
    default:
      errNum = WLZ_ERR_PARAM_DATA;
      break;
    }
  }

  /* check input object */
  if( errNum == WLZ_ERR_NONE ){
    if( obj == NULL ){
      errNum = WLZ_ERR_OBJECT_NULL;
    }
    else {
      switch( obj->type ){
      case WLZ_2D_DOMAINOBJ:
	if( obj->domain.core == NULL ){
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	break;

      case WLZ_3D_DOMAINOBJ:
	return WlzDomainFillHoles3d(obj, size, connectivity, dstErr);

      case WLZ_TRANS_OBJ:
	if( values.obj = WlzDomainFillHoles(obj, size, connectivity,
					    &errNum) ){
	  return WlzMakeMain(WLZ_TRANS_OBJ, obj->domain, values,
			     NULL, NULL, dstErr);
	}
	break;

      case WLZ_EMPTY_OBJ:
	return WlzMakeEmpty(dstErr);

      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
      }
    }
  }

  /* 2D object with domain */
  /* build background object and extract holes */
  if( errNum == WLZ_ERR_NONE ){
    if( domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
					 obj->domain.i->line1 - 1,
					 obj->domain.i->lastln + 1,
					 obj->domain.i->kol1 - 1,
					 obj->domain.i->lastkl + 1,
					 &errNum) ){
      values.core = NULL;
      obj1 = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
			 NULL, NULL, &errNum);
      if( backObj = WlzDiffDomain(obj1, obj, &errNum) ){
	backObj = WlzAssignObject(backObj, &errNum);
	errNum = WlzLabel(backObj, &numObjs, &objs, 1024, 0, connectivity);
	WlzFreeObj(backObj);
      }
      WlzFreeObj(obj1);
    }
  }

  /* fill holes */
  if( errNum == WLZ_ERR_NONE ){
    if( numObjs > 0 ){
      holeObj = WlzAssignObject(WlzMakeEmpty(&errNum), NULL);
      for(i=0; i < numObjs; i++){
	if( WlzArea(objs[i], &errNum) <= size ){
	  obj1 = WlzAssignObject(WlzUnion2(holeObj, objs[i], &errNum),
				 NULL);
	  WlzFreeObj(holeObj);
	  holeObj = obj1;
	}
	WlzFreeObj(objs[i]);
      }
    }
    AlcFree(objs);

    rtnObj = WlzUnion2(obj, holeObj, &errNum);
    WlzFreeObj(holeObj);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

static void usage(
  char	*str)
{
  fprintf(stderr,
	  "Usage:\n"
	  "%s -f -p -s# -h -v file1.wlz file2.wlz .... fileN.wlz\n"
	  "\tRead in all the domains and correct flying pixels\n"
	  "\tand small holes by reassigning to the closest neighbour\n"
	  "\tdomains. For the holes this is the domain itself, for\n"
	  "\tdisconnected pixels it is the domain with greatest area\n"
	  "\tof intersection with the dilated pixel. Each corrected\n"
	  "\tdomain is output to a file with extension _new, additions\n"
	  "\tand deletions are output to files with extensions _adds\n"
	  "\tan _dels respectively\n"
	  "Arguments:\n"
	  "\t-f         fill holes\n"
	  "\t-p         prune & reassign flying pixels\n"
	  "\ts#:        size parameter to determine selection, default 2\n"
	  "\t-h         print this message\n"
	  "\t-v         verbose operation\n"
	  "\n",
	  str);

  return;
}

typedef struct _FixDomainListItem {
  char		*file;
  WlzObject	*origDom;
  WlzObject	*newDom;
  WlzObject	*addsDom;
  WlzObject	*delsDom;
} FixDomainListItem;

void freeFixDomainListItem(
  void		*item)
{
  FixDomainListItem	*dmnItem = (FixDomainListItem *) item;

  AlcFree(dmnItem->file);
  if( dmnItem->origDom ){
    WlzFreeObj(dmnItem->origDom);
  }
  if( dmnItem->newDom ){
    WlzFreeObj(dmnItem->newDom);
  }
  if( dmnItem->addsDom ){
    WlzFreeObj(dmnItem->addsDom);
  }
  if( dmnItem->delsDom ){
    WlzFreeObj(dmnItem->delsDom);
  }
  AlcFree(item);

  return;
}

int main(
  int   argc,
  char  **argv)
{
  FILE		*inFile, *outFile;
  char 		optList[] = "fps:hv";
  int		option;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  char		*errMsg;
  int		verboseFlg=0;
  int		fillFlg=0, pruneFlg=0;
  int		size=2;
  HGUDlpList	*dmnList;
  HGUDlpListItem	*item, *tmpItem;
  FixDomainListItem	*newItem, *maxDomItem;
  WlzObject	*obj1, *obj2, *backObj, **objs;
  WlzDomain	domain;
  WlzValues	values;
  int		i, p, numObjs, maxDomVal;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'f':
      fillFlg = 1;
      break;

    case 'p':
      pruneFlg = 1;
      break;

    case 's':
      size = atoi(optarg);
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

  /* kludge argv[0] to reduce printout */
  if(verboseFlg){
    char	*str;
    int		c;

    c = '/';
    if( !(str = strrchr(argv[0], c)) ){
      str = argv[0];
    }
    else {
      str++;
    }
    argv[0] = str;
  }

  if( verboseFlg ){
    for(i=0; i < argc; i++){
      fprintf(stderr, "%s ", argv[i]);
    }
    fprintf(stderr, "\n");
  }

  /* read the domains - chuck away any grey-values */
  dmnList = HGUDlpListCreate(NULL);
  while( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      errNum = WLZ_ERR_FILE_OPEN;
    }
    else {
      if(verboseFlg){
	fprintf(stderr, "%s: reading %s\n", argv[0], *(argv+optind));
      }
      if( obj1 = WlzAssignObject(WlzReadObj(inFile, &errNum), NULL) ){
	values.core = NULL;
	obj2 = WlzMakeMain(obj1->type, obj1->domain, values,
			   NULL, NULL, &errNum);
	newItem = (FixDomainListItem *)
	  AlcCalloc(sizeof(FixDomainListItem), 1);
	newItem->file = strdup(*(argv+optind));
	newItem->origDom = WlzAssignObject(obj2, &errNum);
	(void) HGUDlpListInsert(dmnList, NULL, newItem,
				freeFixDomainListItem);
	WlzFreeObj(obj1);
      }
      else {
	fprintf(stderr, "%s: can't read object in file %s\n",
		argv[0], *(argv+optind));
      }
    }
    optind++;
  }

  if( errNum == WLZ_ERR_NONE ){
    if(verboseFlg){
      fprintf(stderr, "%s: read %d domains\n", argv[0],
	      HGUDlpListCount(dmnList));
    }
  }

  /* build a full volume to define a surrounding "background"  */
  if( errNum == WLZ_ERR_NONE ){
    if(verboseFlg){
      fprintf(stderr, "%s: union of each domain:", argv[0]);
    }
    objs = (WlzObject **) AlcCalloc(sizeof(WlzObject *),
				    HGUDlpListCount(dmnList));
    obj2 = NULL;
    i = 0;
    item = HGUDlpListHead(dmnList);
    while( item ){
      newItem = (FixDomainListItem *) HGUDlpListEntryGet(dmnList, item);

      objs[i++] = WlzAssignObject(newItem->origDom, &errNum);

      if(verboseFlg){
	fprintf(stderr, ".");
      }

      item = HGUDlpListNext(dmnList, item);
    }
    obj2 = WlzUnionN(i, objs, 0, &errNum);
    for(; i > 0; i--){
      WlzFreeObj(objs[i-1]);
    }
    AlcFree(objs);

    if(verboseFlg){
      fprintf(stderr, " done,\n%s: build background object: ", argv[0]);
    }

    /* now a background object */
    switch( obj2->type ){

    case WLZ_3D_DOMAINOBJ:
      domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
				    obj2->domain.p->plane1 - 1,
				    obj2->domain.p->lastpl + 1,
				    obj2->domain.p->line1 - 1,
				    obj2->domain.p->lastln + 1,
				    obj2->domain.p->kol1 - 1,
				    obj2->domain.p->lastkl + 1,
				    &errNum);
      values.core = NULL;
      obj1 = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
			 NULL, NULL, &errNum);
      for(p=obj1->domain.p->plane1; p <= obj1->domain.p->lastpl; p++){
	domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
					 obj2->domain.p->line1 - 1,
					 obj2->domain.p->lastln + 1,
					 obj2->domain.p->kol1 - 1,
					 obj2->domain.p->lastkl + 1,
					 &errNum);
	obj1->domain.p->domains[p-obj1->domain.p->plane1] = 
	  WlzAssignDomain(domain, &errNum);
      }
      break;

    case WLZ_2D_DOMAINOBJ:
      domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
					 obj2->domain.p->line1 - 1,
					 obj2->domain.p->lastln + 1,
					 obj2->domain.p->kol1 - 1,
					 obj2->domain.p->lastkl + 1,
					 &errNum);
      values.core = NULL;
      obj1 = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
			 NULL, NULL, &errNum);
      break;

    default:
      return WLZ_ERR_OBJECT_TYPE;
    }

    backObj = WlzDiffDomain(obj1, obj2, &errNum);
    backObj = WlzAssignObject(backObj, &errNum);
    WlzFreeObj(obj1);
    WlzFreeObj(obj2);

    if(verboseFlg){
      fprintf(stderr, " done.\n%s: add background to domain list.\n",
	      argv[0]);
    }

    /* add the background to the list */
    newItem = (FixDomainListItem *)
      AlcCalloc(sizeof(FixDomainListItem), 1);
    newItem->file = strdup("background");
    newItem->origDom = WlzAssignObject(backObj, &errNum);
    (void) HGUDlpListInsert(dmnList, NULL, newItem,
			    freeFixDomainListItem);
    WlzFreeObj(backObj);
  }
  

  /* now start correcting domains */
  if( errNum == WLZ_ERR_NONE ){
    if(verboseFlg){
      fprintf(stderr, "%s: start correcting domains\n", argv[0]);
    }

    /* first do the holes */
    if( fillFlg ){
      if(verboseFlg){
	fprintf(stderr, "Fill holes:\n");
      }
      item = HGUDlpListTail(dmnList);
      while( item ){
	newItem = (FixDomainListItem *) HGUDlpListEntryGet(dmnList, item);
	if( strcmp(newItem->file, "background") ){

	  if(verboseFlg){
	    fprintf(stderr, "%s", newItem->file);
	  }

	  /* fill holes */
	  if( obj1 = WlzDomainFillHoles(newItem->origDom, size,
					WLZ_4_CONNECTED, &errNum) ){
	    newItem->newDom = WlzAssignObject(obj1, &errNum);

	    /* remove intersection with other domains */
	    tmpItem = HGUDlpListTail(dmnList);
	    while( tmpItem ){
	      FixDomainListItem *domItem;
	      domItem = (FixDomainListItem *)
		HGUDlpListEntryGet(dmnList, tmpItem);
	      if( strcmp(domItem->file, newItem->file) ){
		if( domItem->newDom ){
		  if( obj2 = WlzDiffDomain(domItem->newDom,
					   obj1, &errNum) ){
		    obj2 = WlzAssignObject(obj2, NULL);
		    WlzFreeObj(domItem->newDom);
		    domItem->newDom = obj2;
		  }
		}
		else {
		  if( obj2 = WlzDiffDomain(domItem->origDom,
					   obj1, &errNum) ){
		    obj2 = WlzAssignObject(obj2, NULL);
		    WlzFreeObj(domItem->origDom);
		    domItem->origDom = obj2;
		  }
		}
	      }
	      tmpItem = HGUDlpListPrev(dmnList, tmpItem);
	    }
	  }
	  
	  if(verboseFlg){
	    fprintf(stderr, " - done\n");
	  }
	}
	item = HGUDlpListPrev(dmnList, item);
      }
    }

    /* now do the flying pixels */
    if( pruneFlg ){
      int	totalFPObjs=0;
      if(verboseFlg){
	fprintf(stderr, "Prune flying pixels:\n");
      }

      item = HGUDlpListTail(dmnList);
      while( item ){
	newItem = (FixDomainListItem *) HGUDlpListEntryGet(dmnList, item);
	if( strcmp(newItem->file, "background") ){

	  if(verboseFlg){
	    fprintf(stderr, "%s: ", newItem->file);
	  }

	  /* locate flying pixels and re-assign */
	  if( newItem->newDom ){
	    errNum = WlzDomainSizeSelect(newItem->newDom, size, 1,
					 WLZ_4_CONNECTED,
					 &numObjs, &objs);
	  }
	  else {
	    errNum = WlzDomainSizeSelect(newItem->origDom, size, 1,
					 WLZ_4_CONNECTED,
					 &numObjs, &objs);
	  }

	  if(verboseFlg){
	    fprintf(stderr, " %d flying pixel objects\n", numObjs);
	    totalFPObjs += numObjs;
	  }

	  for(i=0; (errNum == WLZ_ERR_NONE) && (i < numObjs); i++){
	    /* find maximum intersection with the dilated
	       flying pixel */
	    obj1 = WlzDilation(objs[i], WLZ_8_CONNECTED, &errNum);
	    obj1 = WlzAssignObject(obj1, NULL);
	    maxDomItem = NULL;
	    maxDomVal = 0;
	    tmpItem = HGUDlpListTail(dmnList);
	    while( tmpItem ){
	      FixDomainListItem *domItem;
	      int		val;

	      domItem = (FixDomainListItem *)
		HGUDlpListEntryGet(dmnList, tmpItem);
	      if( strcmp(domItem->file, newItem->file) ){
		if( domItem->newDom ){
		  obj2 = WlzIntersect2(obj1, domItem->newDom, &errNum);
		}
		else {
		  obj2 = WlzIntersect2(obj1, domItem->origDom, &errNum);
		}
		if( obj2 ){
		  if( obj2->type == WLZ_2D_DOMAINOBJ ){
		    val = WlzArea(obj2, &errNum);
		  }
		  else {
		    val = WlzVolume(obj2, &errNum);
		  }
		  if( val > maxDomVal ){
		    maxDomVal = val;
		    maxDomItem = domItem;
		  }
		  WlzFreeObj(obj2);
		}
		else {
		  val = 0;
		}
	      }
	      tmpItem = HGUDlpListPrev(dmnList, tmpItem);
	    }
	    WlzFreeObj(obj1);

	    /* assign to maximum overlap domain and
	       remove from current*/
	    if( maxDomItem ){
	      if( maxDomItem->newDom ){
		obj1 = WlzUnion2(maxDomItem->newDom, objs[i], &errNum);
		WlzFreeObj(maxDomItem->newDom);
		maxDomItem->newDom = WlzAssignObject(obj1, &errNum);
	      }
	      else {
		maxDomItem->newDom = WlzAssignObject(objs[i], &errNum);
	      }
	    }
	    if( newItem->newDom ){
	      obj1 = WlzDiffDomain(newItem->newDom, objs[i], &errNum);
	      WlzFreeObj(newItem->newDom);
	      newItem->newDom = WlzAssignObject(obj1, &errNum);
	    }
	    else {
	      obj1 = WlzDiffDomain(newItem->origDom, objs[i], &errNum);
	      newItem->newDom = WlzAssignObject(obj1, &errNum);
	    }
	    
	    WlzFreeObj(objs[i]);
	    if(verboseFlg){
	      fprintf(stderr, ".");
	    }
	  }
	  AlcFree(objs);
	  if(verboseFlg){
	    fprintf(stderr, " done\n", numObjs);
	  }
	}
	item = HGUDlpListPrev(dmnList, item);
      }

      if(verboseFlg){
	fprintf(stderr, "%s: total flying pixel objects fixed: %d\n",
		argv[0], totalFPObjs);
      }
    }
  }

  /* write out all changes */
  if( errNum == WLZ_ERR_NONE ){
    if(verboseFlg){
      fprintf(stderr, "%s: write out changes\n", argv[0]);
    }
    item = HGUDlpListTail(dmnList);
    while( item ){
      newItem = (FixDomainListItem *) HGUDlpListEntryGet(dmnList, item);
      if( strcmp(newItem->file, "background") ){

	/* get filename body - write to current directory */
	char *fileBody, *str, strBuf[256];
	int	c;

	if(verboseFlg){
	  fprintf(stderr, "\t%s: ", newItem->file);
	}

	c = '/';
	if( !(str = strrchr(newItem->file, c)) ){
	  str = newItem->file;
	}
	else {
	  str++;
	}
	fileBody = strdup(str);
	if( str = strstr((const char *) fileBody, (const char *) ".wlz") ){
	  str[0] = '\0';
	}

	/* write new, adds and dels */
	sprintf(strBuf, "%s_new.wlz", fileBody);
	if( outFile = fopen(strBuf, "w") ){
	  WlzWriteObj(outFile, newItem->newDom);
	  fclose(outFile);
	  if(verboseFlg){
	    fprintf(stderr, "%s ", strBuf);
	  }
	}
	else {
	  fprintf(stderr, "%s: can't open output file %s\n", argv[0], strBuf);
	}
	
	sprintf(strBuf, "%s_adds.wlz", fileBody);
	if( outFile = fopen(strBuf, "w") ){
	  if( obj1 = WlzDiffDomain(newItem->newDom, newItem->origDom,
				   &errNum) ){
	    WlzWriteObj(outFile, obj1);
	    WlzFreeObj(obj1);
	    if(verboseFlg){
	      fprintf(stderr, "%s ", strBuf);
	    }
	  }
	  fclose(outFile);
	}
	else {
	  fprintf(stderr, "%s: can't open output file %s\n", argv[0], strBuf);
	}
	
	sprintf(strBuf, "%s_dels.wlz", fileBody);
	if( outFile = fopen(strBuf, "w") ){
	  if( obj1 = WlzDiffDomain(newItem->origDom, newItem->newDom,
				   &errNum) ){
	    WlzWriteObj(outFile, obj1);
	    WlzFreeObj(obj1);
	    if(verboseFlg){
	      fprintf(stderr, "%s ", strBuf);
	    }
	  }
	  fclose(outFile);
	}
	else {
	  fprintf(stderr, "%s: can't open output file %s\n", argv[0], strBuf);
	}
	
	if(verboseFlg){
	  fprintf(stderr, "\n");
	}
      }
      item = HGUDlpListPrev(dmnList, item);
    }
  }

  return errNum;
}
