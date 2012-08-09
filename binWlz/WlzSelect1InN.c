#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzSelect1InN_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzSelect1InN.c
* \author       Richard Baldock
* \date         March 1999
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
* \brief	Selects planes 1 in n from a 3D object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzselect1inn "WlzSelect1InN"
*/

/*!
\ingroup BinWlz
\defgroup wlzselect1inn WlzSelect1InN
\par Name
WlzSelect1InN - selects planes 1 in n from a 3D object.
\par Synopsis
\verbatim
WlzSelect1InN [-f#] [-l#] [-s#] [-h] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-f</b></td>
    <td>First plane, default 0.</td>
  </tr>
  <tr> 
    <td><b>-l</b></td>
    <td>Last plane. If the last plane is less than the first then all planes
        are selected. Default -1.</td>
  </tr>
  <tr> 
    <td><b>-n</b></td>
    <td>Step. Default 3, ie one in 3 planes is selected.</td>
  </tr>
</table>
\par Description
Selects planes 1 in n from a 3D object.
\par Examples
\verbatim
WlzSelect1InN -f3 -l13 -n3 in.wlz >out.wlz
\endverbatim
Selects planes 3,6,9 and 12 from in.wlz and writes the 3D object to out.wlz
with planes 1,2,3 and 4.
\par File
\ref WlzSelect1InN.c "WlzSelect1InN.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzsampleobj "WlzSampleObj(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <Wlz.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static WlzObject *WlzSelect1InN(WlzObject *obj,
				int first_pl,
				int last_pl,
				int step,
				WlzErrorNum	*dstErr);

static void usage(char *proc_str)
{
  (void )fprintf(stderr,
	  "Usage:\t%s [-f#] [-l#] [-s#] [-h] [<input file>]\n"
	  "\tSelect planes 1 in n from a 3D object\n"
	  "\tresetting voxel size to suit and\n"
	  "\twriting the new object to standard output\n"
	  "Version: %s\n"
	  "Options:\n"
	  "\t  -f#       first plane (def: 0)\n"
	  "\t  -l#       last plane (less than first => all, def: -1)\n"
	  "\t  -n#       step (def: 3)\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -v        Verbose operation\n",
	  proc_str,
	  WlzVersion());
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *nobj;
  FILE		*inFile;
  int		option;
  char 		optList[] = "f:hl:n:v";
  int		first_pl, last_pl, step;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  int		verboseFlg=0;
    
  /* set the defaults */
  first_pl = 0;
  last_pl = -1;
  step = 3;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'f':
      first_pl = atoi(optarg);
      break;

    case 'l':
      last_pl = atoi(optarg);
      break;

    case 'n':
      step = atoi(optarg);
      break;

    case 'v':
      verboseFlg = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return WLZ_ERR_PARAM_TYPE;

    }
  }
  if( verboseFlg ){
    fprintf(stderr,
	    "%s parameters set:\n"
	    "\tFirst plane, last plane, step: %d, %d, %d\n"
	    "", argv[0], first_pl, last_pl, step);
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return WLZ_ERR_PARAM_TYPE;
    }
  }

  /* read objects and select planes if possible */
  while((obj = WlzReadObj(inFile, &errNum)) != NULL) 
  {
    switch( obj->type )
    {
      case WLZ_3D_DOMAINOBJ:
	if( obj->domain.p->type == WLZ_EMPTY_OBJ ){
	  errNum = WlzWriteObj(stdout, obj);
	  break;
	}
	if( last_pl < first_pl ){
	  first_pl = obj->domain.p->plane1;
	  last_pl = obj->domain.p->lastpl;
	}
	if( (nobj = WlzSelect1InN(obj, first_pl, last_pl, step,
				  &errNum)) != NULL ){
	  errNum = WlzWriteObj(stdout, nobj);
	  WlzFreeObj(nobj);
	}
	break;

      default:
	errNum = WlzWriteObj(stdout, obj);
	break;
    }

    WlzFreeObj(obj);
  }

  return errNum;
}

static WlzObject *WlzSelect1InN(
  WlzObject	*obj,
  int 		first_pl,
  int 		last_pl,
  int 		step,
  WlzErrorNum	*dstErr)
{
  WlzObject	*nobj=NULL;
  int		p, new_p, new_plane1, new_lastpl;
  WlzDomain	domain;
  WlzValues	values;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check the object */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  if( errNum == WLZ_ERR_NONE ){
    domain.core = NULL;
    values.core = NULL;
    switch( obj->type ){

    case WLZ_3D_DOMAINOBJ:
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
	break;
      }
      if( obj->domain.core->type == WLZ_EMPTY_DOMAIN ){
	return WlzMakeEmpty(dstErr);
      }
      break;

    case WLZ_EMPTY_OBJ:
      return WlzMakeEmpty(dstErr);

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  /* check the parameters */
  if( errNum == WLZ_ERR_NONE ){
    if( first_pl > last_pl ){
      errNum = WLZ_ERR_PARAM_DATA;
    }
    else if( step <= 0 ){
      errNum = WLZ_ERR_PARAM_DATA;
    }
    /* check for empty object */
    else if((first_pl > obj->domain.p->lastpl) ||
	    (last_pl < obj->domain.p->plane1) ){
      return WlzMakeEmpty(dstErr);
    }
  }

  /* make a new planedomain and voxel table as required */
  if( errNum == WLZ_ERR_NONE ){
    for(p=first_pl; p <= last_pl; p += step){
      if( p >= obj->domain.p->plane1 ){
	new_plane1 = p / step;
	break;
      }
    }
    for(;p <= last_pl; p += step){
      if( p > obj->domain.p->lastpl ){
	break;
      }
      new_lastpl = p / step;
    }

    if((domain.p = WlzMakePlaneDomain(obj->domain.p->type,
				      new_plane1, new_lastpl,
				      obj->domain.p->line1,
				      obj->domain.p->lastln,
				      obj->domain.p->kol1,
				      obj->domain.p->lastkl,
				      &errNum)) != NULL){
      domain.p->voxel_size[0] = obj->domain.p->voxel_size[0];
      domain.p->voxel_size[1] = obj->domain.p->voxel_size[1];
      domain.p->voxel_size[2] = step * (obj->domain.p->voxel_size[2]);

      if( obj->values.core != NULL ){
	if( (values.vox =
	     WlzMakeVoxelValueTb(obj->values.vox->type,
				 new_plane1, new_lastpl,
				 obj->values.vox->bckgrnd,
				 obj, &errNum)) == NULL ){
	  (void) WlzFreePlaneDomain(domain.p);
	}
      }
    }
  }

  /* set domain and voxel values */
  if( errNum == WLZ_ERR_NONE ){
    for(p=first_pl; p <= last_pl; p += step){
      if( p < obj->domain.p->plane1 ){
	continue;
      }
      if( p > obj->domain.p->lastpl ){
	break;
      }
      new_p = (p / step) - domain.p->plane1;
      domain.p->domains[new_p] =
	WlzAssignDomain(obj->domain.p->domains[p - obj->domain.p->plane1],
			NULL);
      if( obj->values.core != NULL ){
	values.vox->values[new_p] =
	  WlzAssignValues(	
	    obj->values.vox->values[p - obj->domain.p->plane1], NULL);
      }
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    nobj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
		       NULL, NULL, &errNum);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return nobj;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
