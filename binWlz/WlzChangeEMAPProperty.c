#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzChangeEMAPProperty_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file		binWlz/WlzChangeEMAPProperty.c
* \author       Richard Baldock
* \date         July 2002
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
* \brief	Changes the EMAP property of a woolz object.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzchangeemapproperty "WlzChangeEMAPProperty"
*/


/*!
\ingroup      BinWlz
\defgroup     wlzchangeemapproperty WlzChangeEMAPProperty
\par Name
WlzChangeEMAPProperty - changes the EMAP property of a woolz object.
\par Synopsis
\verbatim
WlzChangeEMAPProperty [-a<anatomyUID>] [-e<modelUID>] [-E<targetUID>]
                      [-c<time>] [-C<author>] [-f<filename>] [F<target version>]
                      [-l<location>] [-m<time>] [-M<author>]
                      [-n<name>] [-r] [-s<stage>] [-S<sub-stage>]
                      [-t#] [-T<text>] [-V<version>]
                      [-h] [-v][<input file>]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-a</b></td>
    <td>anatomy component unique ID</td>
  </tr>
  <tr>
    <td><b>-e</b></td>
    <td>EMAP model unique ID</td>
  </tr>
  <tr>
    <td><b>-E</b></td>
    <td>target EMAP model unique ID for transform object</td>
  </tr>
  <tr>
    <td><b>-c</b></td>
    <td>Creation time (dd.mm.yyyy) </td>
  </tr>
  <tr>
    <td><b>-C</b></td>
    <td>Creation author </td>
  </tr>
  <tr>
    <td><b>-f</b></td>
    <td>File name associated with object </td>
  </tr>
  <tr>
    <td><b>-F</b></td>
    <td>Target model version if the object is a transform </td>
  </tr>
  <tr>
    <td><b>-l</b></td>
    <td>Machine name (or IP address) </td>
  </tr>
  <tr>
    <td><b>-m</b></td>
    <td>Modification time (dd.mm.yyyy) </td>
  </tr>
  <tr>
    <td><b>-M</b></td>
    <td>Modification author. </td>
  </tr>
  <tr>
    <td><b>-n</b></td>
    <td>EMAP Model name </td>
  </tr>
  <tr>
    <td><b>-r</b></td>
    <td>Remove the property. </td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td>Embryonic stage (default NULL). </td>
  </tr>
  <tr>
    <td><b>-S</b></td>
    <td>Embryonic sub-stage (default NULL). </td>
  </tr>
  <tr>
    <td><b>-t</b></td>
    <td>EMAP property type: </td>
  </tr>
  <tr>
    <td><b> </b></td>
    <td># = 1: grey-level EMAP-model </td>
  </tr>
  <tr>
    <td><b> </b></td>
    <td># = 2: EMAP anatomy domain </td>
  </tr>
  <tr>
    <td><b> </b></td>
    <td># = 3: other domain (default) </td>
  </tr>
  <tr>
    <td><b>-T</b></td>
    <td>Comment text </td>
  </tr>
  <tr>
    <td><b>-V</b></td>
    <td>Model version </td>
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
Change the EMAP property of a woolz object
writing the new object to standard output
Note times and authors are put in from the environment
automatically. Be sure to put in appropriate escapes
for spaces and special characters.

\par Description
The EMAP property holds information about the woolz object. 
If it is a grey-level object the the information will probably define 
the model and version of a standard EMAP atlas model and the EMAP property
type will be set to WLZ_EMAP_PROPERTY_GREY_MODEL. The property is used
track which models have been used for domains submitted to the EMAGE
database and which version of the painted domains have been used. For
domains the property type will be WLZ_EMAP_PROPERTY_DOMAIN_OTHER or
WLZ_EMAP_PROPERTY_DOMAIN_ANATOMY. Grey-level data that has been mapped or
associated with a specific EMAP model is ot type WLZ_EMAP_PROPERTY_GREY_OTHER
and a transform defining the mapping between two models is of type
WLZ_EMAP_PROPERTY_TRANSFORM.

\par Examples
\verbatim
#
# set the properties for a grey-level atlas model
# and display the result using WlzFacts
#
richard-lt.local% WlzChangeEMAPProperty -e EMAPM:00025 -C "Richard Baldock" \
-n "embryo_1_3D" -T "Program test" -V "1.0" -t "1" -f "embryo_1_3D.wlz" \
-s 14 ts14.wlz | WlzFacts
  Object type: WLZ_3D_DOMAINOBJ.
    Linkcount: 0.
    Domain type: WLZ_PLANEDOMAIN_DOMAIN.
    Linkcount: 1.
    Domain plane bounds: 0 305
    Domain line bounds: 51 223
    Domain column bounds: 60 380
    VoxelSize: 4 4 7
    Values type: WLZ_VOXELVALUETABLE_GREY.
    Linkcount: 1.
    Background type: WLZ_GREY_UBYTE.
    Background value: 255.
    Values plane bounds: 0 305
    Property list:
    Property type:WLZ_PROPERTY_EMAP
      Linkcount: 1.
      EMAP property type: WLZ_EMAP_PROPERTY_GREY_MODEL.
      Model UID: EMAPM:00025
      Anatomy UID: 
      Target UID: 
      Theiler Stage: 14
      Model name: embryo_1_3D
      Model version: 1.0
      Filename: embryo_1_3D.wlz
      Creation time: Thu Feb 23 17:36:42 2006
      Creation author: Richard Baldock
      Creation machine name: richard-lt.local
      Modification time: Thu Feb 23 17:36:42 2006
      Modification author: richard (richard)
      Comment: Program test
richard-lt.local% 
\endverbatim

\par File
\ref WlzChangeEMAPProperty.c "WlzChangeEMAPProperty.c"
\par See Also
\ref wlzfacts "WlzFacts(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
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
	  "Usage:\t%s [-a<UID>] [-e<UID>] [-E<UID>] [-F<target version>] "
	  "[-c<time>] [-C<author>] [-f<filename>] "
	  "[-l<location>] [-m<time>] [-M<author>] [-n<name>] "
	  "[-r] [-s<stage>] [-S<sub-stage>] [-t#] [-T<text>]"
	  "[-V<version>] [-h] [-v] [<input file>]\n"
	  "\tChange the EMAP property of a woolz object\n"
	  "\twriting the new object to standard output\n"
	  "\tNote times and authors are put in from the environment\n"
	  "\tautomatically. Be sure to put in appropriate escapes\n"
	  "\tfor spaces and special characters.\n"
	  "\tOptions are:\n"
	  "\t  -a<anatomy UID> anatomy UID\n"
	  "\t  -e<model UID> model EMAP UID\n"
	  "\t  -E<target UID> target model UID\n"
	  "\t  -c<time>     Creation time (dd.mm.yyyy)\n"
	  "\t  -C<author>   Creation author\n"
	  "\t  -f<filename> File name associated with object\n"
	  "\t  -F<version>  Target model version if transform object\n"
	  "\t  -l<location> Machine name (or IP address)\n"
	  "\t  -m<time>     Modification time (dd.mm.yyyy)\n"
	  "\t  -M<author>   Modification author\n"
	  "\t  -n<name>     EMAP Model name\n"
	  "\t  -r           Remove the property\n"
	  "\t  -s<stage>    Embryo stage (default NULL)\n"
	  "\t  -S<sub-stage>    Embryo sub-stage (default NULL)\n"
	  "\t  -t#          EMAP property type:\n"
	  "\t               # = %d: grey-level EMAP-model\n"
	  "\t                   %d: grey-level data other\n"
	  "\t                   %d: EMAP anatomy domain\n"
	  "\t                   %d: other domain (default)\n"
	  "\t                   %d: transform\n"
	  "\t  -T<text>     Comment text\n"
	  "\t  -V<version>  Model version\n"
	  "\t  -h           Help - prints this usage message\n"
	  "\t  -v           Verbose operation\n"
	  "",
	  proc_str, WLZ_EMAP_PROPERTY_GREY_MODEL,
	  WLZ_EMAP_PROPERTY_GREY_OTHER,
	  WLZ_EMAP_PROPERTY_DOMAIN_ANATOMY,
	  WLZ_EMAP_PROPERTY_DOMAIN_OTHER,
	  WLZ_EMAP_PROPERTY_TRANSFORM);
  return;
}

int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj;
  FILE		*inFile;
  char 		optList[] = "a:e:E:c:C:f:F:l:m:M:n:rs:S:t:T:V:hv";
  int		option;
  WlzEMAPPropertyType	emapType=WLZ_EMAP_PROPERTY_DOMAIN_OTHER;
  int		emapTypeFlg=0;
  char		*modelUID=NULL;
  char		*anatomyUID=NULL;
  char		*targetUID=NULL;
  char		*targetVersion=NULL;
  char		*stage=NULL, *subStage=NULL;
  char		*modelName=NULL;
  char		*version=NULL;
  char		*fileName=NULL;
  time_t	creationTime;
  int		cTimeFlg=0;
  char		*creationAuthor=NULL;
  char		*creationMachineName=NULL;
  time_t	modificationTime;
  int		mTimeFlg=0;
  char		*modificationAuthor=NULL;
  char		*comment=NULL;
  int		removeFlg=0;
  WlzProperty	property;
  struct tm	tm;
  WlzErrorNum	errNum;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'a':
      anatomyUID = optarg;
      break;

    case 'e':
      modelUID = optarg;
      break;

    case 'E':
      targetUID = optarg;
      break;

    case 'c':
#if defined(BSD)  || defined(_WIN32)
      if( sscanf(optarg, "%d/%d/%d", &(tm.tm_mday), &(tm.tm_mon), &(tm.tm_year)) == 3 ){
	tm.tm_sec = 0;
	tm.tm_min = 0;
	tm.tm_hour = 0;
	if( tm.tm_year > 1900 ){
	  tm.tm_year -= 1900;
	}
#else 
      if( strptime(optarg, "%d/%m/%Y", &tm) ){
#endif /* BSD */ 
	creationTime = mktime(&tm);
	cTimeFlg = 1;
      }
      else {
	usage(argv[0]);
	return( 1 );
      }
      break;

    case 'C':
      creationAuthor = optarg;
      break;

    case 'f':
      fileName = optarg;
      break;

    case 'F':
      targetVersion = optarg;
      break;

    case 'l':
      creationMachineName = optarg;
      break;

    case 'm':
#if defined(BSD)  || defined(_WIN32)
      if( sscanf(optarg, "%d/%d/%d", &(tm.tm_mday), &(tm.tm_mon), &(tm.tm_year)) == 3 ){
	tm.tm_sec = 0;
	tm.tm_min = 0;
	tm.tm_hour = 0;
	if( tm.tm_year > 1900 ){
	  tm.tm_year -= 1900;
	}
#else 
      if( strptime(optarg, "%d/%m/%Y", &tm) ){
#endif /* BSD */ 
	modificationTime = mktime(&tm);
	mTimeFlg = 1;
      }
      else {
	usage(argv[0]);
	return( 1 );
      }
      break;

    case 'M':
      modificationAuthor = optarg;
      break;

    case 'n':
      modelName = optarg;
      break;

    case 'r':
      removeFlg = 1;
      break;

    case 's':
      stage = optarg;
      break;

    case 'S':
      subStage = optarg;
      break;

    case 't':
      switch( emapType = (WlzEMAPPropertyType) atoi(optarg) ){

      case WLZ_EMAP_PROPERTY_GREY_MODEL:
      case WLZ_EMAP_PROPERTY_GREY_OTHER:
      case WLZ_EMAP_PROPERTY_DOMAIN_ANATOMY:
      case WLZ_EMAP_PROPERTY_DOMAIN_OTHER:
      case WLZ_EMAP_PROPERTY_TRANSFORM:
	break;

      default:
        fprintf(stderr, "%s: EMAP property type = %d is invalid\n",
		argv[0], emapType );
        usage(argv[0]);
        return 1;

      }
      emapTypeFlg=1;
      break;

    case 'T':
      comment = optarg;
      break;

    case 'V':
      version = optarg;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return( 1 );

    }
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
  }

  /* read objects and convert if possible */
  while( (obj = WlzReadObj(inFile, NULL)) != NULL ){
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
    case WLZ_TRANS_OBJ:
    case WLZ_3D_WARP_TRANS:
    case WLZ_PROPERTY_OBJ:
      /* get the EMAP property, if absent add one */
      if( obj->plist == NULL ){
	obj->plist = WlzMakePropertyList(NULL);
      }
      if( obj->plist ){
	property = WlzGetProperty(obj->plist->list, WLZ_PROPERTY_EMAP, NULL);
	if( property.core == NULL ){
	  property.emap = WlzMakeEMAPProperty(WLZ_EMAP_PROPERTY_GREY_MODEL,
					      NULL, NULL, NULL, NULL, NULL,
					      NULL, NULL, NULL, NULL, NULL,
					      &errNum);
	  WlzAssignProperty(property, NULL);
	  AlcDLPListEntryAppend(obj->plist->list, NULL, (void *) property.emap,
				WlzFreePropertyListEntry);
	}
      }
      if( removeFlg ){
	WlzRemoveProperty(obj->plist->list, property);
      }
      else {
	if( !emapTypeFlg ){
	  emapType = property.emap->emapType;
	}
	WlzChangeEMAPProperty(property.emap, emapType, modelUID,
			      anatomyUID, targetUID, targetVersion, stage,
			      subStage, modelName, version, fileName, comment);
	if( cTimeFlg ){
	  property.emap->creationTime = creationTime;
	}
	if( creationAuthor ){
	  strcpy(property.emap->creationAuthor, creationAuthor);
	}
	if( creationMachineName ){
	  strcpy(property.emap->creationMachineName, creationMachineName);
	}
	if( mTimeFlg ){
	  property.emap->modificationTime = modificationTime;
	}
	if( modificationAuthor ){
	  strcpy(property.emap->modificationAuthor, modificationAuthor);
	}
      }
      break;

    case WLZ_COMPOUND_ARR_1:
    case WLZ_COMPOUND_ARR_2:
      break;

    default:
      break;
    }

    (void )WlzWriteObj(stdout, obj);
    WlzFreeObj(obj);
  }

  return WLZ_ERR_NONE;
}
		

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
