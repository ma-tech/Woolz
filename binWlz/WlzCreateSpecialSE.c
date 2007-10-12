#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzCreateSpecialSE_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzCreateSpecialSE.c
* \author       Richard Baldock
* \date         May 2001
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
* \brief	Creates a "special" structuring element object.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzcreatespecialse "WlzCreateSpecialSE"
*/

/*!
\ingroup BinWlz
\defgroup wlzcreatespecialse WlzCreateSpecialSE
\par Name
WlzCreateSpecialSE  -  creates a "special" structuring element object.
\par Synopsis
\verbatim
WlzCreateSpecialSE  -  [-o#] [-r#] [-s#] [-t#] [-x#] [-y#] [-z#] [-h]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Structuring element object type specified by an integer value:
    <table width="500" border="0">
    <tr> <td>0</td> <td>2D<td> </tr>
    <tr> <td>1</td> <td>3D<td> </tr>
    </table>
    </td>
  </tr>
  <tr> 
    <td><b>-r</b></td>
    <td>Structuring element radius.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Structuring element sub-type specified by an integer value:
    <table width="500" border="0">
    <tr> <td>0</td> <td>8-connected (2D)</td> </tr>
    <tr> <td>1</td> <td>4-connected (2D)</td> </tr>
    <tr> <td>2</td> <td>6-connected (3D)</td> </tr>
    <tr> <td>3</td> <td>18-connected (3D)</td> </tr>
    <tr> <td>4</td> <td>26-connected (3D)</td> </tr>
    <tr> <td>5</td> <td>octagonal (2D)</td> </tr>
    <tr> <td>6</td> <td>Euclidean (2D and 3D)</td> </tr>
    </table>
    </td>
  </tr>
  <tr>
    <td><b>-t</b></td>
    <td>Structuring element type specified by an integer value
        (Liang naming):
    <table width="500" border="0">
    <tr> <td>0</td> <td>h4 - 4-connected (default)</td> </tr>
    <tr> <td>1</td> <td>ex4 - (111)(010)(111)</td> </tr>
    <tr> <td>2</td> <td>a8 - 8-connected</td> </tr>
    <tr> <td>3</td> <td>h6 - 6 combinations of 2 missing corners</td> </tr>
    <tr> <td>4</td> <td>h5 - 4 options of 3 missing corners</td> </tr>
    <tr> <td>5</td> <td>h7 - 4 options of 1 missing corner</td> </tr>
    <tr> <td>6</td> <td>a3 - 4-connected, 4 options of 1 missing side</td> </tr>
    <tr> <td>7</td> <td>e1 - (000)(011)(000)</td> </tr>
    <tr> <td>8</td> <td>e2 - (000)(111)(000)</td> </tr>
    <tr> <td>9</td> <td>v2 - (010)(010)(010)</td> </tr>
    <tr> <td>20</td> <td>single pixel object</td> </tr>
    <tr> <td>21</td> <td>circle object</td> </tr>
    <tr> <td>22</td> <td>sphere object</td> </tr>
    <tr> <td>23</td> <td>standard structuring element - uses
                         distance type (sub-type), radius and otype</td> </tr>
    </table>
    </td>
  </tr>
</table>
\par Description
Creates a "special" structuring element object
and writes the new object to standard output.
\par Examples
\verbatim
WlzCreateSpecialSE -o1 -r32 -s4 >sphere.wlz
\endverbatim
Creates a 3D domain object in which the domain is a 26-connected sphere of
radius 32.
\par File
\ref WlzCreateSpecialSE.c "WlzCreateSpecialSE.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzdilation "WlzDilation(1)"
\ref wlzerosion "WlzErosion(1)"
\ref WlzMakeSpecialStructElement "WlzMakeSpecialStructElement(3)"
\ref WlzMakeSinglePixelObject "WlzMakeSinglePixelObject(3)"
\ref WlzMakeCircleObject "WlzMakeCircleObject(3)"
\ref WlzMakeSphereObject "WlzMakeSphereObject(3)"
\ref WlzMakeStdStructElement "WlzMakeStdStructElement(3)"
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
	  "Usage:\t%s [-o#] [-r#] [-s#] [-t#] [-x#] [-y#] [-z#] [-h] \n"
	  "\tCreate a \"special\" structuring element object\n"
	  "\tand write the new object to standard output\n"
	  "\tOptions are:\n"
	  "\t  -o#       Structuring element object type\n"
	  "\t            # =  0: 2D\n"
	  "\t            # =  1: 3D\n"
	  "\t  -r#       Structuring element radius\n"
	  "\t  -s#       Structuring element sub-type (default = 0)\n"
	  "\t            when type 23 selected this is the distance type:\n"
	  "\t            # =  0: 8-connected (2D)\n"
	  "\t            # =  1: 4-connected (2D)\n"
	  "\t            # =  2: 6-connected (2D)\n"
	  "\t            # =  3: 18-connected (3D)\n"
	  "\t            # =  4: 26-connected (3D)\n"
	  "\t            # =  5: octagonal (2D)\n"
	  "\t            # =  6: Euclidean (2D & 3D)\n"
	  "\t  -t#       Structuring element type (Liang naming):\n"
	  "\t            # =  0: h4 - 4-connected (default)\n"
	  "\t            # =  1: ex4 - (111)(010)(111)\n"
	  "\t            # =  2: a8 - 8-connected\n"
	  "\t            # =  3: h6 - 6 combinations of 2 missing corners\n"
	  "\t            # =  4: h5 - 4 options of 3 missing corners\n"
	  "\t            # =  5: h7 - 4 options of 1 missing corner\n"
	  "\t            # =  6: a3 - 4-connected, 4 options of 1 missing side\n"
	  "\t            # =  7: e1 - (000)(011)(000)\n"
	  "\t            # =  8: e2 - (000)(111)(000)\n"
	  "\t            # =  9: v2 - (010)(010)(010)\n"
	  "\t            # = 20: single pixel object\n"
	  "\t            # = 21: circle object\n"
	  "\t            # = 22: sphere object\n"
	  "\t            # = 23: standard structuring element - uses\n"
	  "\t                    distance type (sub-type), radius and otype\n"
	  "\t  -h        Help - prints this usage message\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject			*obj;
  char 				optList[] = "ho:r:s:t:x:y:z:";
  int				option;
  WlzSpecialStructElmType	seType;
  int				type, subType;
  WlzErrorNum			errNum=WLZ_ERR_NONE;
  WlzObjectType			oType;
  WlzDistanceType		dType;
  double			radius, x, y, z;

  /* set defaults */
  type = 0;
  subType = 0;
  seType = WLZ_SPEC_STRUCT_ELM_H4;
  oType = WLZ_2D_DOMAINOBJ;
  dType = WLZ_8_DISTANCE;
  radius = 1.0;
  x = y = z = 0.0;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'o':
      switch( atoi(optarg) ){
      default:
      case 0:
	oType = WLZ_2D_DOMAINOBJ;
	break;

      case 1:
	oType = WLZ_3D_DOMAINOBJ;
	break;
      }
      break;

    case 'r':
      radius = atof(optarg);
      if( radius < 1.0 ){
	fprintf(stderr, "%s: radius negative - reset to 1.0\n", argv[0]);
	radius = 1.0;
      }
      break;

    case 't':
      switch( type = atoi(optarg) ){

      case 0:
	seType = WLZ_SPEC_STRUCT_ELM_H4;
	break;
	
      case 1:
	seType = WLZ_SPEC_STRUCT_ELM_EX4;
	break;
	
      case 2:
	seType = WLZ_SPEC_STRUCT_ELM_A8;
	break;
	
      case 3:
	seType = WLZ_SPEC_STRUCT_ELM_H6;
	break;
	
      case 4:
	seType = WLZ_SPEC_STRUCT_ELM_H5;
	break;
	
      case 5:
	seType = WLZ_SPEC_STRUCT_ELM_H7;
	break;
	
      case 6:
	seType = WLZ_SPEC_STRUCT_ELM_A3;
	break;
	
      case 7:
	seType = WLZ_SPEC_STRUCT_ELM_E1;
	break;
	
      case 8:
	seType = WLZ_SPEC_STRUCT_ELM_E2;
	break;
	
      case 9:
	seType = WLZ_SPEC_STRUCT_ELM_V2;
	break;
	
      case 20:
      case 21:
      case 22:
      case 23:
	break;
	
      default:
        fprintf(stderr, "%s: structuring element type = %d is invalid\n",
		argv[0], type );
        usage(argv[0]);
        return WLZ_ERR_PARAM_DATA;

      }
      break;

    case 's':
      switch( subType = atoi(optarg) ){

      default:
      case 0:
	dType = WLZ_8_DISTANCE;
	break;

      case 1:
	dType = WLZ_4_DISTANCE;
	break;

      case 2:
	dType = WLZ_6_DISTANCE;
	break;

      case 3:
	dType = WLZ_18_DISTANCE;
	break;

      case 4:
	dType = WLZ_26_DISTANCE;
	break;

      case 5:
	dType = WLZ_OCTAGONAL_DISTANCE;
	break;

      case 6:
	dType = WLZ_EUCLIDEAN_DISTANCE;
	break;

      }
	  
      break;

    case 'h':
    default:
      usage(argv[0]);
      return WLZ_ERR_UNSPECIFIED;

    }
  }

  /* create structuring element */
  switch( type ){
  case 0:
  case 1:
  case 2:
  case 3:
  case 4:
  case 5:
  case 6:
  case 7:
  case 8:
  case 9:
    if((obj = WlzMakeSpecialStructElement(seType, subType, &errNum)) != NULL){
      WlzWriteObj(stdout, obj);
      WlzFreeObj(obj);
    }
    break;

  case 20:
    if((obj = WlzMakeSinglePixelObject(oType, (int) x, (int) y, (int) z,
				       &errNum)) != NULL){
      WlzWriteObj(stdout, obj);
      WlzFreeObj(obj);
    }
    break;

  case 21:
    if((obj = WlzMakeCircleObject(radius, x, y, &errNum)) != NULL){
      WlzWriteObj(stdout, obj);
      WlzFreeObj(obj);
    }
    break;

  case 22:
    if((obj = WlzMakeSphereObject(WLZ_3D_DOMAINOBJ, radius, x, y, z,
				  &errNum)) != NULL){
      WlzWriteObj(stdout, obj);
      WlzFreeObj(obj);
    }
    break;

  case 23:
    if((obj = WlzMakeStdStructElement(oType, dType, radius,
    				      &errNum)) != NULL){
      WlzWriteObj(stdout, obj);
      WlzFreeObj(obj);
    }
    break;

  }

  /* trap the WLZ_ERR_READ_EOF since this is a legal way of indicating
     the end of objects in a file */
  if( errNum == WLZ_ERR_READ_EOF ){
    errNum = WLZ_ERR_NONE;
  }
  return errNum;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
