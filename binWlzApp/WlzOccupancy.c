#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzOccupancy.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Mon Jan  9 14:21:56 2006
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par Copyright:
* Copyright (C) 2005 Medical research Council, UK.
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
* \ingroup      BinWlzApp
* \brief        Calculate the occupancy of each pixel with respect
 to a series of domains.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/
 
/*!
\ingroup      BinWlzApp
\defgroup     wlzoccupancy WlzOccupancy
\par Name
WlzOccupancy - calculate pixel occupancy with respect to a series of domains.
\par Synopsis
\verbatim
WlzOccupancy [-d <domainfile>] [-m] [-n#] [-h] [-v] [<input file>]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-d</b></td>
    <td>Domain over which the occupancy will be calculated. The default
 is the union of the input domains.</td>
  </tr>
  <tr>
    <td><b>-m</b></td>
    <tdCalculate the normalised mean, i.e. multiply by (255/n).</td>
  </tr>
  <tr>
    <td><b>-n</b></td>
    <td>Maximum number of objects, default = 100.</td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Help - print help message</td>
  </tr>
  <tr>
    <td><b>-v</b></td>
    <td>Verbose operation</td>
  </tr>
</table>

\par Description
Calculate the occupancy for each pixel location within the target domain.
 The occupancy is the number of times the location is included within the
 input domain series. The target domain is an optional input or taken as
 the union of the domain series. The output object will have integer
 (WLZ_GREY_INT) grey type and the value will either be the actual
 occupancy or the normalise occupancy which will have a maximum value
 of 255.

\par Examples
\verbatim
\endverbatim

\par See Also
\par Bugs
None known
*/
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <Wlz.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s [-d <domainfile>] [-n#] [-h] [-v] [<input file>]\n"
	  "\tCalculate the occupancy given a series of domains.\n"
	  "\tOptions are:\n"
	  "\t  -d <file> Optional input domain over which occupancy\n"
	  "\t            will be calculated, default - union of input.\n"
	  "\t  -m        Calculate a normalised mean (1.e. *255/n).\n"
	  "\t  -n#       Maximum number of objects -default=100\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -v        Verbose operation\n"
	  "",
	  proc_str);
  return;
}

WlzObject *WlzAddValuesTable(
  WlzObject	*obj,
  WlzGreyType	gtype,
  WlzPixelV	bckgrnd,
  WlzErrorNum	*dstErr)
{
  WlzObject	*rtnObj=NULL;
  WlzObjectType	type;
  WlzValues	values, *valuess;
  WlzDomain	*domains;
  WlzObject	*tmpObj;
  int		p;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check inputs */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
      type = WlzGreyTableType(WLZ_GREY_TAB_RAGR, gtype, NULL);
      values.v = WlzNewValueTb(obj, type, bckgrnd, &errNum);
      rtnObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, obj->domain, values,
			   NULL, NULL, &errNum);
      break;

    case WLZ_3D_DOMAINOBJ:
      values.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
				       obj->domain.p->plane1,
				       obj->domain.p->lastpl,
				       bckgrnd, NULL, &errNum);
      rtnObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, obj->domain, values,
			   NULL, NULL, &errNum);
      domains = rtnObj->domain.p->domains;
      valuess = rtnObj->values.vox->values;
      for(p=0;
	  p < (rtnObj->domain.p->lastpl - rtnObj->domain.p->plane1 + 1);
	  p++, domains++, valuess++){
	if( (*domains).core ){
	  values.core = NULL;
	  tmpObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains, values,
			       NULL, NULL, NULL);
	  values.v = WlzNewValueTb(tmpObj, type, bckgrnd, &errNum);
	  *valuess = WlzAssignValues(values, NULL);
	  WlzFreeObj(tmpObj);
	}
      }
      break;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

WlzErrorNum	WlzGreyScalarAddValue(
  WlzObject	*obj,
  WlzPixelV	val)
{
  WlzObject		*tmpObj;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzGreyP		gptr;
  WlzPixelV		tmpVal;
  WlzDomain		*domains;
  WlzValues		*values;
  int			i, nplanes;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check object */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      if( obj->domain.i == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else if( obj->values.v == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      /* check planedomain and voxeltable */
      if( obj->domain.p == NULL ){
	errNum =  WLZ_ERR_DOMAIN_NULL;
      }
      else if( obj->domain.p->type != WLZ_PLANEDOMAIN_DOMAIN ){
	errNum = WLZ_ERR_DOMAIN_TYPE;
      }
      else if( obj->values.vox == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
      }
      else if( obj->values.vox->type != WLZ_VOXELVALUETABLE_GREY ){
	errNum = WLZ_ERR_VALUES_TYPE;
      }
      else {
	/* set range of each plane if non-empty - indicated by NULL */
	domains = obj->domain.p->domains;
	values = obj->values.vox->values;
	nplanes = obj->domain.p->lastpl - obj->domain.p->plane1 + 1;
	for(i=0; (errNum == WLZ_ERR_NONE) && (i < nplanes);
	    i++, domains++, values++){

	  if( (*domains).core == NULL || (*values).core == NULL ){
	    continue;
	  }

	  if( tmpObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains, *values,
				    NULL, NULL,	NULL) ){
	    errNum = WlzGreyScalarAddValue(tmpObj, val);
	    WlzFreeObj( tmpObj );
	  }
	}
      }
      return errNum;

    case WLZ_TRANS_OBJ:
      return WlzGreyScalarAddValue(obj->values.obj, val);

    case WLZ_EMPTY_OBJ:
      return errNum;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    errNum = WlzInitGreyScan(obj, &iwsp, &gwsp);
    WlzValueConvertPixel(&tmpVal, val, gwsp.pixeltype);
    while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){

      gptr = gwsp.u_grintptr;
      switch (gwsp.pixeltype) {

      case WLZ_GREY_INT:
	for (i=0; i<iwsp.colrmn; i++, gptr.inp++)
	  *gptr.inp += tmpVal.v.inv;
	break;

      case WLZ_GREY_SHORT:
	for (i=0; i<iwsp.colrmn; i++, gptr.shp++)
	  *gptr.shp += tmpVal.v.shv;
	break;

      case WLZ_GREY_UBYTE:
	for (i=0; i<iwsp.colrmn; i++, gptr.ubp++)
	  *gptr.ubp += tmpVal.v.ubv;
	break;

      case WLZ_GREY_FLOAT:
	for (i=0; i<iwsp.colrmn; i++, gptr.flp++)
	  *gptr.flp += tmpVal.v.flv;
	break;

      case WLZ_GREY_DOUBLE:
	for (i=0; i<iwsp.colrmn; i++, gptr.dbp++)
	  *gptr.dbp += tmpVal.v.dbv;
	break;

      case WLZ_GREY_RGBA:
	for (i=0; i<iwsp.colrmn; i++, gptr.rgbp++)
	  *gptr.rgbp += tmpVal.v.rgbv;
	break;

      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;
      }
    }
    if( errNum == WLZ_ERR_EOO ){
      errNum = WLZ_ERR_NONE;
    }
  }

  return errNum;
}
 
WlzErrorNum	WlzGreyScalarMultValue(
  WlzObject	*obj,
  WlzPixelV	val)
{
  WlzObject		*tmpObj;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzGreyP		gptr;
  WlzPixelV		tmpVal;
  WlzDomain		*domains;
  WlzValues		*values;
  int			i, nplanes;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check object */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      if( obj->domain.i == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else if( obj->values.v == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      /* check planedomain and voxeltable */
      if( obj->domain.p == NULL ){
	errNum =  WLZ_ERR_DOMAIN_NULL;
      }
      else if( obj->domain.p->type != WLZ_PLANEDOMAIN_DOMAIN ){
	errNum = WLZ_ERR_DOMAIN_TYPE;
      }
      else if( obj->values.vox == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
      }
      else if( obj->values.vox->type != WLZ_VOXELVALUETABLE_GREY ){
	errNum = WLZ_ERR_VALUES_TYPE;
      }
      else {
	/* set range of each plane if non-empty - indicated by NULL */
	domains = obj->domain.p->domains;
	values = obj->values.vox->values;
	nplanes = obj->domain.p->lastpl - obj->domain.p->plane1 + 1;
	for(i=0; (errNum == WLZ_ERR_NONE) && (i < nplanes);
	    i++, domains++, values++){

	  if( (*domains).core == NULL || (*values).core == NULL ){
	    continue;
	  }

	  if( tmpObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains, *values,
				    NULL, NULL,	NULL) ){
	    errNum = WlzGreyScalarMultValue(tmpObj, val);
	    WlzFreeObj( tmpObj );
	  }
	}
      }
      return errNum;

    case WLZ_TRANS_OBJ:
      return WlzGreyScalarMultValue(obj->values.obj, val);

    case WLZ_EMPTY_OBJ:
      return errNum;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    errNum = WlzInitGreyScan(obj, &iwsp, &gwsp);
    WlzValueConvertPixel(&tmpVal, val, gwsp.pixeltype);
    while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){

      gptr = gwsp.u_grintptr;
      switch (gwsp.pixeltype) {

      case WLZ_GREY_INT:
	for (i=0; i<iwsp.colrmn; i++, gptr.inp++)
	  *gptr.inp *= tmpVal.v.inv;
	break;

      case WLZ_GREY_SHORT:
	for (i=0; i<iwsp.colrmn; i++, gptr.shp++)
	  *gptr.shp *= tmpVal.v.shv;
	break;

      case WLZ_GREY_UBYTE:
	for (i=0; i<iwsp.colrmn; i++, gptr.ubp++)
	  *gptr.ubp *= tmpVal.v.ubv;
	break;

      case WLZ_GREY_FLOAT:
	for (i=0; i<iwsp.colrmn; i++, gptr.flp++)
	  *gptr.flp *= tmpVal.v.flv;
	break;

      case WLZ_GREY_DOUBLE:
	for (i=0; i<iwsp.colrmn; i++, gptr.dbp++)
	  *gptr.dbp *= tmpVal.v.dbv;
	break;

      case WLZ_GREY_RGBA:
	for (i=0; i<iwsp.colrmn; i++, gptr.rgbp++)
	  *gptr.rgbp *= tmpVal.v.rgbv;
	break;

      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;
      }
    }
    if( errNum == WLZ_ERR_EOO ){
      errNum = WLZ_ERR_NONE;
    }
  }

  return errNum;
}
 
WlzErrorNum	WlzGreyScalarDivValue(
  WlzObject	*obj,
  WlzPixelV	val)
{
  WlzObject		*tmpObj;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzGreyP		gptr;
  WlzPixelV		tmpVal;
  WlzDomain		*domains;
  WlzValues		*values;
  int			i, nplanes;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check object */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      if( obj->domain.i == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else if( obj->values.v == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      /* check planedomain and voxeltable */
      if( obj->domain.p == NULL ){
	errNum =  WLZ_ERR_DOMAIN_NULL;
      }
      else if( obj->domain.p->type != WLZ_PLANEDOMAIN_DOMAIN ){
	errNum = WLZ_ERR_DOMAIN_TYPE;
      }
      else if( obj->values.vox == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
      }
      else if( obj->values.vox->type != WLZ_VOXELVALUETABLE_GREY ){
	errNum = WLZ_ERR_VALUES_TYPE;
      }
      else {
	/* set range of each plane if non-empty - indicated by NULL */
	domains = obj->domain.p->domains;
	values = obj->values.vox->values;
	nplanes = obj->domain.p->lastpl - obj->domain.p->plane1 + 1;
	for(i=0; (errNum == WLZ_ERR_NONE) && (i < nplanes);
	    i++, domains++, values++){

	  if( (*domains).core == NULL || (*values).core == NULL ){
	    continue;
	  }

	  if( tmpObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains, *values,
				    NULL, NULL,	NULL) ){
	    errNum = WlzGreyScalarDivValue(tmpObj, val);
	    WlzFreeObj( tmpObj );
	  }
	}
      }
      return errNum;

    case WLZ_TRANS_OBJ:
      return WlzGreyScalarDivValue(obj->values.obj, val);

    case WLZ_EMPTY_OBJ:
      return errNum;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    errNum = WlzInitGreyScan(obj, &iwsp, &gwsp);
    WlzValueConvertPixel(&tmpVal, val, gwsp.pixeltype);
    while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){

      gptr = gwsp.u_grintptr;
      switch (gwsp.pixeltype) {

      case WLZ_GREY_INT:
	for (i=0; i<iwsp.colrmn; i++, gptr.inp++)
	  *gptr.inp /= tmpVal.v.inv;
	break;

      case WLZ_GREY_SHORT:
	for (i=0; i<iwsp.colrmn; i++, gptr.shp++)
	  *gptr.shp /= tmpVal.v.shv;
	break;

      case WLZ_GREY_UBYTE:
	for (i=0; i<iwsp.colrmn; i++, gptr.ubp++)
	  *gptr.ubp /= tmpVal.v.ubv;
	break;

      case WLZ_GREY_FLOAT:
	for (i=0; i<iwsp.colrmn; i++, gptr.flp++)
	  *gptr.flp /= tmpVal.v.flv;
	break;

      case WLZ_GREY_DOUBLE:
	for (i=0; i<iwsp.colrmn; i++, gptr.dbp++)
	  *gptr.dbp /= tmpVal.v.dbv;
	break;

      case WLZ_GREY_RGBA:
	for (i=0; i<iwsp.colrmn; i++, gptr.rgbp++)
	  *gptr.rgbp /= tmpVal.v.rgbv;
	break;

      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;
      }
    }
    if( errNum == WLZ_ERR_EOO ){
      errNum = WLZ_ERR_NONE;
    }
  }

  return errNum;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj1, *obj2, *obj, **objlist, *rtnObj;
  WlzObjectType	type = (WlzObjectType) -1;
  int 		i, n, nmax, p;
  FILE		*inFile;
  char 		optList[] = "d:mn:hv";
  int		option;
  int		meanFlg=0;
  int		verboseFlg = 0;
  WlzPixelV	bckgrnd;
  WlzValues	values;
  const char	*errMsg;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* setup the return object type and background */
  bckgrnd.type = WLZ_GREY_INT;
  bckgrnd.v.inv = 0;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  nmax = 100;
  rtnObj = NULL;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

      /* read in a target domain over which to capture the occupancy,
	 set grey type to WLZ_GREY_INT and value to zero */
    case 'd':
      if( inFile = fopen(optarg, "rb") ){
	if( obj = WlzReadObj(inFile, &errNum) ){
	  if( (rtnObj = WlzAddValuesTable(obj, WLZ_GREY_INT,
					  bckgrnd, &errNum)) == NULL ){
	    (void )WlzStringFromErrorNum(errNum, &errMsg);
	    fprintf(stderr,
		    "%s: Failed to add values table to occupancy object: %s\n",
		    argv[0], errMsg);
	    return 1;
	  }
	  errNum = WlzGreySetValue(rtnObj, bckgrnd);
	  WlzFreeObj(obj);
	}
	else {
	  fprintf(stderr, "%s: Can't read occupancy domain file\n", argv[0]);
	  usage(argv[0]);
	}
	fclose(inFile);
      }
      else {
	fprintf(stderr, "%s: Can't open occupancy domain file\n", argv[0]);
	usage(argv[0]);
      }
      break;

    case 'm':
      meanFlg = 1;
      break;

    case 'n':
      nmax = atoi(optarg);
      if( nmax < 1 ){
	fprintf(stderr, "%s: nmax = %d is invalid\n", argv[0], nmax);
	usage(argv[0]);
	return( 1 );
      }
      break;

    case 'v':
      verboseFlg = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return( 1 );

    }
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "rb")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return( 1 );
    }
  }

  /* allocate space for the object pointers */
  if( (objlist = (WlzObject **)
       AlcMalloc(sizeof(WlzObject *) * nmax)) == NULL ){
    (void )fprintf(stderr, "%s: memory allocation failed.\n",
    		   argv[0]);
    return( 1 );
  }

  /* read objects accumulating compatible types */
  n = 0;
  while( ((obj = WlzAssignObject(WlzReadObj(inFile, NULL),
  			         NULL)) != NULL) && (n < nmax) ) {

    if( type == -1 &&
	(obj->type == WLZ_2D_DOMAINOBJ || obj->type == WLZ_3D_DOMAINOBJ) ){
      type = obj->type;
    }

    if( obj->type == type ){
      objlist[n++] = WlzAssignObject(obj, NULL);
    } else {
      WlzFreeObj( obj );
    }
  }

  if( type == WLZ_EMPTY_OBJ ){
    return( 0 );
  }

  /* check for occupancy object */
  if( rtnObj ){
    /* check type against return object - must be the same */
    if( type != rtnObj->type ){
      fprintf(stderr, "%s: Occupancy domain object type does not match\n"
	      " inout domain type\n", argv[0]);
      usage(argv[0]);
      return 1;
    }
  }
  else {
    /* use the union of input domains as the occupancy object */
    if((obj1 = WlzUnionN(n, objlist, 1, &errNum)) == NULL) {
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: failed to perform union (%s).\n",
		     argv[0], errMsg);
      return(1);
    }
    rtnObj = WlzAddValuesTable(obj1, WLZ_GREY_INT, bckgrnd, &errNum);
    errNum = WlzGreySetValue(rtnObj, bckgrnd);
    WlzFreeObj(obj1);
  }

  /* now add 1 to each pixel for each domain */
  bckgrnd.v.inv = 1;
  for(i=0; i < n; i++){
    if( obj1 = WlzIntersect2(rtnObj, objlist[i], &errNum) ){
      if( obj1->type != WLZ_EMPTY_OBJ ){
	if( obj1->type == WLZ_2D_DOMAINOBJ ){
	  obj1->values = WlzAssignValues(rtnObj->values, &errNum);
	}
	else {
	  values.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
					   obj1->domain.p->plane1,
					   obj1->domain.p->lastpl,
					   bckgrnd, NULL, &errNum);
	  obj1->values = WlzAssignValues(values, &errNum);
	  for(p=obj1->domain.p->plane1; p <= obj1->domain.p->lastpl; p++){
	    values.vox->values[p-obj1->domain.p->plane1]
	      = WlzAssignValues(rtnObj->values.vox->values[p-rtnObj->domain.p->plane1],
				&errNum);
	  }
	}
	errNum = WlzGreyScalarAddValue(obj1, bckgrnd);
      }
      WlzFreeObj(obj1);
    }
  }

  /* output occupancy object */
  if( rtnObj ){
    if( meanFlg ){
      bckgrnd.v.inv = 255;
      errNum = WlzGreyScalarMultValue(rtnObj, bckgrnd);
      bckgrnd.v.inv = n;
      errNum = WlzGreyScalarDivValue(rtnObj, bckgrnd);
    }
    WlzWriteObj(stdout, rtnObj);
  }

  /* freespace so purify can check for leaks */
  WlzFreeObj(obj1);
  while( n-- ){
    WlzFreeObj(objlist[n]);
  }
  AlcFree((void *) objlist);

  return( 0 );
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
