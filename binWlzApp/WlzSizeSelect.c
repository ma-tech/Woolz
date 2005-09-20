#pragma ident "MRC HGU $Id$"
/*!
* \file         binWlzApp/WlzSizeSelect.c
* \author	Richard Baldock
* \date         July 1999
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
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
* \brief	Selects components of an object's domain on the basis
* 		of area.
* \ingroup	BinWlzApp
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzsizeselect "WlzSizeSelect"
*/

/*!
\ingroup BinWlzApp
\defgroup wlzsizeselect WlzSizeSelect
\par Name
WlzSizeSelect - selects components of an object's domain on the basis of
                area.
\par Synopsis
\verbatim
WlzSizeSelect [-h] [-v] -a<mesh_area> -c<conn> -H -m
              [-f] [-g] [-h] [-i] [-l] [-L] [-p]
              [-s #] [-m #] [-c #,#,#] [-C #,#,#] <in bib file>
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Prints usage information.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Verbose output.</td>
  </tr>
  <tr>
    <td>c</td>
    <td>Conectivity.</td>
  </tr>
  <tr> 
    <td><b>-H</b></td>
    <td>Select from the domain complement within
        the bounding box, default - use the foreground.
    </td>
  </tr>
  <tr> 
    <td><b>-m</b></td>
    <td>Keep domains <= area, the default is to keep domains > area.</td>
  </tr>
</table>
\par Description
WlzSizeSelect segments an object and select parts according to area,
either keeping all parts <= a given area or > the area.
The operation is applied by default to the object domain but can be applied to
the domain complement in which case the size filter is applied
to the holes. If the input object is 3D then the filter is
applied plane-by-plane.
\par Examples
\verbatim
WlzSizeSelect -c8 < in.wlz > out.wlz
\endverbatim
Select all foreground objects, defined by segmentation with
connectivity 8 and area > 5 pixels (i.e. discard small blobs
with area <= 5).

\verbatim
WlzSizeSelect -a10 -c4 -m < in.wlz > out.wlz
\endverbatim
Select all foreground objects, defined by segmentation with
connectivity 4 and area <= 10 pixels (i.e. collect small blobs
with area <= 10).

\verbatim
WlzSizeSelect -a10 -c4 -H -m < in.wlz > out.wlz
\endverbatim
Select all background objects, defined by segmentation with
connectivity 4 and area <= 10 pixels (i.e. collect holes
with area <= 10).
\par File
\ref WlzSizeSelect.c "WlzSizeSelect.c"
\par See Also
\ref BinWlzApp "WlzIntro(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
 
#include <Wlz.h>

static WlzObject *WlzSizeFilter(
  WlzObject	*nobj,
  int		mesh_area,
  int		mesh_flag,
  WlzConnectType connectFlg)
{
  WlzObject	**objs;
  WlzObject	*new_obj, *obj1;
  int           num_obj, i, max_i, size, max_size;

  if( (nobj == NULL) || (nobj->domain.core == NULL) )
  {
    return( NULL );
  }

  (void) WlzLabel(nobj, &num_obj, &objs, 2047, 0, connectFlg);
  if( num_obj <= 0 )
  {
    return( NULL );
  }

  /* select objects less than or equal to mesh_area */
  new_obj = NULL;
  for(i=0; i < num_obj; i++)
  {
    if( WlzArea(objs[i], NULL) <= mesh_area )
    {
      if( new_obj )
      {
	obj1 = WlzAssignObject(WlzUnion2(new_obj, objs[i], NULL), NULL);
	WlzFreeObj( new_obj );
	new_obj = WlzAssignObject(obj1, NULL);
      }
      else
      {
	new_obj = WlzAssignObject(WlzMakeMain(WLZ_2D_DOMAINOBJ,
					      objs[i]->domain, objs[i]->values,
					      NULL, NULL, NULL), NULL);
      }
    }

    if( objs[i] )
    {
      WlzFreeObj( objs[i] );
    }
  }
  AlcFree((void *) objs);

  if( !mesh_flag )
  {
    if( new_obj )
    {
      obj1 = WlzDiffDomain( nobj, new_obj, NULL );
      WlzFreeObj( new_obj );
      new_obj = obj1;
    }
    else
    {
      new_obj = WlzMakeMain(WLZ_2D_DOMAINOBJ,
			    nobj->domain, nobj->values,
			    NULL, NULL, NULL);
    }
  }

  return( new_obj );
}
 
static void usage(
  char	*str)
{
  fprintf(stderr,
	  "Usage:\n"
	  "%s -a<mesh_area> -c<conn> -H -m -h -v < infile > outfile\n"
	  "\tSegment an object and select parts according to area, either\n"
	  "\tkeeping all parts <= a given area or > the area. The operation\n"
	  "\tis applied by default to the object domain but can be applied to\n"
	  "\tthe domain complement in which case the size filter is applied\n"
	  "\tto the holes. If the input object is 3D then the filter is\n"
	  "\tapplied plane-by-plane.\n"
	  "Arguments:\n"
	  "\tmesh_area: domains <= area selected, default 5\n"
	  "\tconn:      connectivity for segmentation = 4 or 8, - default 4\n"
	  "\t-H:        select from the domain complement within\n"
	  "\t           the bounding box, default - use the foreground\n"
	  "\t-m:        keep domains <= area,\n"
	  "\t           default - keep domains > area\n"
	  "\t-h         print this message\n"
	  "\t-v         verbose operation\n"
	  "Examples:\n"
	  "\tSelect all foreground objects, defined by segmentation with\n"
	  "\tconnectivity 8 and area > 5 pixels (i.e. discard small blobs\n"
	  "\twith area <= 5):\n"
	  "\t\t%s -c8 < in.wlz > out.wlz\n\n"
	  "\tSelect all foreground objects, defined by segmentation with\n"
	  "\tconnectivity 4 and area <= 10 pixels (i.e. collect small blobs\n"
	  "\twith area <= 10):\n"
	  "\t\t%s -a10 -c4 -m < in.wlz > out.wlz\n\n"
	  "\tSelect all background objects, defined by segmentation with\n"
	  "\tconnectivity 4 and area <= 10 pixels (i.e. collect holes\n"
	  "\twith area <= 10):\n"
	  "\t\t%s -a10 -c4 -H -m < in.wlz > out.wlz\n\n"
	  "\n",
	  str, str, str, str);

  return;
}

int main(
  int   argc,
  char  **argv)
{
  WlzObject     	*obj, *nobj, *max_obj, *tmpobj, *obj1;
  WlzPlaneDomain	*planedmn, *new_planedmn;
  int           	mesh_area = 5;
  int			mesh_flag = 0;
  int			hole_flag = 0;
  int			i, plane, nplanes;
  WlzConnectType	connectFlg=WLZ_4_CONNECTED;
  int			connectivity=4;
  int			verboseFlg=0;
  char			*nameStr;
  WlzDomain		domain, *domains;
  WlzValues		values;
  
  /* check input arguments */
  nameStr = argv[0];
  argv++;
  argc--;
  while( argc > 0 )
  {
    switch( argv[0][1] )
    {
    case 'a':
      if( sscanf(*argv,"-a%d",&mesh_area) < 1 ){
	fprintf(stderr,"Setting mesh_area to 5\n");
	mesh_area = 5;
      }
      break;

    case 'c':
      if( sscanf(*argv,"-c%d",&connectivity) < 1 ){
	fprintf(stderr,"Setting connectivity to 4\n");
	connectFlg = WLZ_4_CONNECTED;
	connectivity = 4;
      }
      else if( connectivity == 4 ){
	connectFlg = WLZ_4_CONNECTED;
	connectivity = 4;
      }
      else if( connectivity == 8 ){
	connectFlg = WLZ_8_CONNECTED;
	connectivity = 8;
      }
      else {
	fprintf(stderr,"%s: Unknown connectivity - setting to 4\n");
	connectFlg = WLZ_4_CONNECTED;
	connectivity = 4;
      }
      break;

    case 'H':
      hole_flag = 1;
      break;

    case 'h':
      usage(nameStr);
      return 0;

    case 'm':
      mesh_flag = 1;
      break;

    case 'v':
      verboseFlg = 1;
      break;

    default:
      usage(nameStr);
      exit(1);
    }
    argc--;
    argv++;
  }

  if( verboseFlg ){
    fprintf(stderr,
	    "%s: mesh_area = %d,"
	    " hole_flag = %d, mesh_flag = %d, connectivity = %d\n",
	    nameStr, mesh_area, hole_flag, mesh_flag, connectivity);
  }

  while((obj = WlzReadObj(stdin, NULL)) != NULL) 
  {
    obj = WlzAssignObject(obj, NULL);
    switch( obj->type )
    {
      case WLZ_2D_DOMAINOBJ:
	if( hole_flag )
	{
	  domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
					   obj->domain.i->line1 - 2,
					   obj->domain.i->lastln + 2,
					   obj->domain.i->kol1 - 2,
					   obj->domain.i->lastkl + 2,
					   NULL);
	  values.core = NULL;

	  nobj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
			     NULL, NULL, NULL);
	  WlzAssignObject(nobj, NULL);
	  tmpobj = WlzDiffDomain(nobj, obj, NULL);
	  WlzFreeObj(obj);
	  WlzFreeObj(nobj);
	  obj = WlzAssignObject(tmpobj, NULL);
	}

        nobj = WlzSizeFilter(obj, mesh_area, mesh_flag, connectFlg);
	if( nobj )
	{
	  (void) WlzWriteObj( stdout, nobj );
	  WlzFreeObj( nobj );
	}
	break;

	/* for 3D objects ignore the grey-table for now */
      case WLZ_3D_DOMAINOBJ:
	planedmn = (WlzPlaneDomain *) obj->domain.p;
	new_planedmn = WlzMakePlaneDomain(planedmn->type,
					  planedmn->plane1, planedmn->lastpl,
					  planedmn->line1, planedmn->lastln,
					  planedmn->kol1, planedmn->lastkl,
					  NULL);
	nplanes = planedmn->lastpl - planedmn->plane1 + 1;
	for(i=0; i < 3; i++)
	{
	  new_planedmn->voxel_size[i] = planedmn->voxel_size[i];
	}

	values.core = NULL;
	for(plane=0; plane < nplanes; plane++)
	{

	  if( planedmn->domains[plane].core == NULL ){
	    new_planedmn->domains[plane].core = NULL;
	    continue;
	  }

	  tmpobj = WlzMakeMain(WLZ_2D_DOMAINOBJ, planedmn->domains[plane],
			       values, NULL, NULL, NULL);
	  tmpobj = WlzAssignObject(tmpobj, NULL);

	  if( hole_flag )
	  {
	    domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
					     tmpobj->domain.i->line1 - 2,
					     tmpobj->domain.i->lastln + 2,
					     tmpobj->domain.i->kol1 - 2,
					     tmpobj->domain.i->lastkl + 2,
					     NULL);
	    nobj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
			       NULL, NULL, NULL);
	    nobj = WlzAssignObject(nobj, NULL);
	    obj1 = WlzDiffDomain(nobj, tmpobj, NULL);
	    WlzFreeObj(tmpobj);
	    WlzFreeObj(nobj);
	    tmpobj = WlzAssignObject(obj1, NULL);
	  }

	  if( nobj = WlzSizeFilter(tmpobj, mesh_area, mesh_flag, connectFlg) )
	  {
	    new_planedmn->domains[plane] = WlzAssignDomain(nobj->domain, NULL);
	    WlzFreeObj( nobj );
	  }
	  else
	  {
	    new_planedmn->domains[plane].core = NULL;
	  }
	}

	domain.p = new_planedmn;
	nobj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
			   NULL, NULL, NULL);
	WlzWriteObj(stdout, nobj);
	WlzFreeObj( nobj );

        break;
 
      default:
        WlzWriteObj(stdout,obj);
        break;
    }
 
    WlzFreeObj(obj);
  }

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
