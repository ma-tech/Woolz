#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCreateRandomDomain_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/WlzCreateRandomDomain.c
* \author       Richard Baldock
* \date         August 2006
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
* \ingroup      WlzBinApp
* \brief        Generate a set of random domains to cover a given model
* 		domain.
*               
* \par Binary
* \ref wlzcreaterandomdomain "WlzCreateRandomDomain"
*/
 
/*!
\ingroup      WlzBinApp
\defgroup     wlzcreaterandomdomain WlzCreateRandomDomain
\par Name
WlzCreateRandomDomain - 
\par Synopsis
\verbatim
WlzCreateRandomDomain

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-</b></td>
    <td> </td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file name.</td>
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

#include <stdlib.h>
#include <Wlz.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void usage(
  char	*str)
{
  fprintf(stderr,
	  "Usage:\n"
	  "%s -n <num domains> -r <radius> -t <type> -h -v  [<model-file>]\n"
	  "\tRead in model domain from file or stdin and generate a set of random\n"
	  "\tdomains within the model. <num domains> domains are created\n"
	  "\taccording to type and written to stdout.If domain sequence is\n"
	  "\tselected then a series of incrementing domains will be generated,\n"
	  "\twith each step a random selctio of pixels from the remaining set.\n"
	  "Arguments:\n"
	  "\t-m         generate a regular mesh of square domains with size\n"
	  "\t           defined by the radius and type, enough domains will\n"
	  "\t           be generated to cover the model.\n"
	  "\t-n#        number of domains to be generate (default 100)\n"
	  "\t-r#        radius parameter for random domains (default 100)\n"
	  "\t-s         create an incrementing sequence of domains\n"
	  "\t-t#        type of random domain default 1\n"
	  "\t             = 1 - circular domain (on a single plane in 3D)\n"
	  "\t             = 2 - spherical domain\n"
	  "\t-h         print this message\n"
	  "\t-v         verbose operation\n"
	  "\n",
	  str);

  return;
}

static int WlzSize(
  WlzObject	*obj,
  WlzErrorNum	*dstErr)
{
  int	size;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( obj ){
    switch(obj->type){
    case WLZ_2D_DOMAINOBJ:
      size = WlzArea(obj, &errNum);
      break;

    case WLZ_3D_DOMAINOBJ:
      size = WlzVolume(obj, &errNum);
      break;

    case WLZ_EMPTY_OBJ:
      size = 0;
      break;

    default:
      size = -1;
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }
  else {
    size = -1;
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return size;
}


int main(
  int   argc,
  char  **argv)
{
  FILE		*inFile;
  char 		optList[] = "mn:r:st:hv";
  int		option;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  int		meshFlg=0;
  int		seqFlg=0;
  int		verboseFlg=0;
  int		type=1;
  int		numDomains=100, domainCount;
  double	radius=5.0;
  int		minSize;
  WlzObject	*obj, *baseDomain, *obj1, *obj2;
  WlzDomain	domain;
  WlzValues	values;
  WlzPixelV	backgrnd;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'm':
      meshFlg = 1;
      break;

    case 'n':
      numDomains = atoi(optarg);
      if( numDomains < 1 ){
	fprintf(stderr, "%s: number of domains requested\n", argv[0]);
	return 1;
      }
      if( numDomains > 10000 ){
	fprintf(stderr, "%s: %d domains requested, are your sure? Quit with <ctrl> c\n",
		argv[0], numDomains);
      }
      break;

    case 'r':
      radius = atof(optarg);
      if( radius <= 1.0){
	radius = 5.0;
	fprintf(stderr, "%s: invalid radius, reset to 5\n", argv[0]);
      }
      break;

    case 's':
      seqFlg = 1;
      break;

    case 't':
      type = atoi(optarg);
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

  /* get the target domain */
  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "rb")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return( 1 );
    }
  }
  if( (obj = WlzReadObj(inFile, &errNum)) == NULL ){
    fprintf(stderr, "%s: failed to read a woolz object from %s\n",
	    argv[0], *(argv+optind));
    usage(argv[0]);
    return( 1 );
  }

  /* generate the base domain */
  switch( obj->type ){
  case WLZ_2D_DOMAINOBJ:
    if( meshFlg ){
      baseDomain = WlzMakeRectangleObject(radius, radius, 0.0, 0.0, &errNum);
    }
    else {
      baseDomain = WlzMakeCircleObject(radius, 0.0, 0.0, &errNum);
    }
    break;

  case WLZ_3D_DOMAINOBJ:
    switch( type ){
    case 1:
      if( meshFlg ){
	baseDomain = WlzMakeCuboidObject(WLZ_3D_DOMAINOBJ, radius, radius, 0.0,
					 0.0, 0.0, 0.0, &errNum);
      }
      else {
	obj1 = WlzMakeCircleObject(radius, 0.0, 0.0, &errNum);
	domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
				      0, 0, 0, 0, 0, 0, &errNum);
	domain.p->domains[0] = WlzAssignDomain(obj1->domain, &errNum);
	values.core = NULL;
	baseDomain = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,	
				 NULL, NULL, &errNum);
	WlzFreeObj(obj1);
      }
      break;

    case 2:
      if( meshFlg ){
	baseDomain = WlzMakeCuboidObject(WLZ_3D_DOMAINOBJ, radius, radius, radius,
					 0.0, 0.0, 0.0, &errNum);
      }
      else {
	baseDomain = WlzMakeSphereObject(WLZ_3D_DOMAINOBJ, radius,
					 0.0, 0.0, 0.0, &errNum);
      }
      break;

    default:
      fprintf(stderr, "%s: invalid random domain type: %d\n",
	      argv[0], type);
      return 1;
    }
    break;

  default:
    fprintf(stderr, "%s: invalid object type from %s, must be a 2D or 3D domain object\n",
	    argv[0], *(argv+optind));
    return 1;
  }
  minSize = WlzSize(baseDomain, &errNum) / 4;

  /* now generate random domains */
  if( seqFlg ){
    int 	pixelsPerStep, numPixels, *indices;
    int		i, j, k, index;
    long	ran;
    double	dran;
    WlzUByte	*values;

    values = AlcCalloc(sizeof(char),
		       (obj->domain.i->lastln - obj->domain.i->line1 + 1) *
		       (obj->domain.i->lastkl - obj->domain.i->kol1 + 1));
    obj1 = WlzMakeRect(obj->domain.i->line1,
		       obj->domain.i->lastln,
		       obj->domain.i->kol1,
		       obj->domain.i->lastkl,
		       WLZ_GREY_UBYTE, (int *) values, backgrnd,
		       NULL, NULL, &errNum);

    /* calc number of pixels to switch foreach step */
    numPixels = WlzArea(obj1, &errNum);
    pixelsPerStep = numPixels / numDomains;
    indices = AlcMalloc(sizeof(int) * numPixels);
    for(i=0; i < numPixels; i++){
      indices[i] = i;
    }

    /* now randomise */
    for(i=0; i < numPixels * 10; i++){
      ran = random();
      dran = ((double) (ran&0xffffff)) / 0xffffff;
      j = dran * numPixels;
      index = indices[j];
      indices[j] = indices[0];
      indices[0] = index;
    }

    /* create the domains */
    k = 0;
    WlzWriteObj(stdout, obj1);
    for(i=0; i < numDomains; i++){
      for(j=0; j < pixelsPerStep; j++, k++){
	values[indices[k]] = 255;
      }
      WlzWriteObj(stdout, obj1);
    }
    
    
  }
  else if( meshFlg ){
    int	xNum, yNum, zNum;
    int xWidth, yWidth, zWidth;
    int i, j, k;
    int	x, y, z = 0;

    domainCount = 0;
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      xWidth = baseDomain->domain.i->lastkl - baseDomain->domain.i->kol1 + 1;
      yWidth = baseDomain->domain.i->lastln - baseDomain->domain.i->line1 + 1;
      xNum = obj->domain.i->lastkl - obj->domain.i->kol1 + 1;
      xNum = (xNum % xWidth) ? xNum/xWidth + 1 : xNum/xWidth;
      yNum = obj->domain.i->lastln - obj->domain.i->line1 + 1;
      yNum = (yNum % yWidth) ? yNum/yWidth + 1 : yNum/yWidth;
      for(j=0; j < yNum; j++){
	y = obj->domain.i->line1 + j * yWidth + yWidth/2;
	for(i=0; i < xNum; i++){
	  x = obj->domain.i->kol1 + i * xWidth + xWidth/2;
	  if((obj1 = WlzShiftObject(baseDomain, x, y, z, &errNum)) != NULL){
	    if((obj2 = WlzIntersect2(obj, obj1, &errNum)) != NULL){
	      if( WlzSize(obj2, &errNum) > minSize ){
		WlzWriteObj(stdout, obj2);
		WlzFreeObj(obj2);
		domainCount++;
	      }
	    }
	    WlzFreeObj(obj1);
	  }
	}
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      xWidth = baseDomain->domain.p->lastkl - baseDomain->domain.p->kol1 + 1;
      yWidth = baseDomain->domain.p->lastln - baseDomain->domain.p->line1 + 1;
      zWidth = baseDomain->domain.p->lastpl - baseDomain->domain.p->plane1 + 1;
      xNum = obj->domain.p->lastkl - obj->domain.p->kol1 + 1;
      xNum = (xNum % xWidth) ? xNum/xWidth + 1 : xNum/xWidth;
      yNum = obj->domain.p->lastln - obj->domain.p->line1 + 1;
      yNum = (yNum % yWidth) ? yNum/yWidth + 1 : yNum/yWidth;
      zNum = obj->domain.p->lastpl - obj->domain.p->plane1 + 1;
      zNum = (zNum % zWidth) ? zNum/zWidth + 1 : zNum/zWidth;
      for(k=0; k < zNum; k++){
	z = obj->domain.p->plane1 + k * yWidth + zWidth/2;
	for(j=0; j < yNum; j++){
	  y = obj->domain.p->line1 + j * yWidth + yWidth/2;
	  for(i=0; i < xNum; i++){
	    x = obj->domain.p->kol1 + i * xWidth + xWidth/2;
	    if((obj1 = WlzShiftObject(baseDomain, x, y, z, &errNum)) != NULL){
	      if((obj2 = WlzIntersect2(obj, obj1, &errNum)) != NULL){
		if( WlzSize(obj2, &errNum) > minSize ){
		  WlzWriteObj(stdout, obj2);
		  WlzFreeObj(obj2);
		  domainCount++;
		}
	      }
	      WlzFreeObj(obj1);
	    }
	  }
	}
      }
      break;

      default:
        break;
    }
  }
  else {
    domainCount= 0;
/*  srandomdev();*/
    while( domainCount < numDomains ){
      long	xRan, yRan, zRan;
      double	xp, yp, zp;
      int		x, y, z;
      xRan = random();
      yRan = random();
      zRan = random();
      xp = ((double) (xRan&0xffffff))/0xffffff;
      yp = ((double) (yRan&0xffffff))/0xffffff;
      zp = ((double) (zRan&0xffffff))/0xffffff;

      switch( obj->type ){
      case WLZ_2D_DOMAINOBJ:
	x = obj->domain.i->kol1 + xp*(obj->domain.i->lastkl - obj->domain.i->kol1);
	y = obj->domain.i->line1 + yp*(obj->domain.i->lastln - obj->domain.i->line1);
	z = 0;
	break;

      case WLZ_3D_DOMAINOBJ:
	x = obj->domain.p->kol1 + xp*(obj->domain.p->lastkl - obj->domain.p->kol1);
	y = obj->domain.p->line1 + yp*(obj->domain.p->lastln - obj->domain.p->line1);
	z = obj->domain.p->plane1 + zp*(obj->domain.p->lastpl - obj->domain.p->plane1);
	break;
      default:
        break;
      }
      if((obj1 = WlzShiftObject(baseDomain, x, y, z, &errNum)) != NULL){
	if((obj2 = WlzIntersect2(obj, obj1, &errNum)) != NULL){
	  if( WlzSize(obj2, &errNum) > minSize ){
	    WlzWriteObj(stdout, obj2);
	    WlzFreeObj(obj2);
	    domainCount++;
	  }
	}
	WlzFreeObj(obj1);
      }
    }
  }

  return 0;
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
