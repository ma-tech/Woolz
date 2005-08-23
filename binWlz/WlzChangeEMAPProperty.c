#pragma ident "MRC HGU $Id$"
/*!
\ingroup      BinWlz
\defgroup     wlzchangeemapproperty WlzChangeEMAPProperty
\par Name
WlzChangeEMAPProperty - Change the EMAP property of a woolz object.
\par Synopsis
\verbatim
WlzChangeEMAPProperty [-c<time>] [-C<author>] [-f<filename>]
                      [-l<location>] [-m<time>] [-M<author>]
                      [-n<name>] [-r] [-s#] [-t#] [-T<text>]
                      [-V<version>]
                      [-h] [-v][<input file>]

\endverbatim
\par Options
<table width="500" border="0">
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
    <td>Theiler stage (default 1). </td>
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
domains the property type will be WLZ_EMAP_PROPERTY_OTHER_DOMAIN or
WLZ_EMAP_PROPERTY_PAINTED_DOMAIN.

\par Examples
\verbatim
#
# set the properties for a grey-level atlas model
# and display the result using WlzFacts
#
WlzChangeEMAPProperty -C "Richard Baldock" -n "embryo_2_WM_left" \
    -T "Program test" -V "1.0" -t "2" -f "embryo_2_WM_left.wlz" \
    -s 15 embryo_2_WM_left.wlz | WlzFacts

  Object type: WLZ_3D_DOMAINOBJ.
    Linkcount: 0.
    Domain type: WLZ_PLANEDOMAIN_DOMAIN.
    Linkcount: 1.
    Domain plane bounds: 0 0
    Domain line bounds: 66 316
    Domain column bounds: 44 274
    VoxelSize: 1 1 1
    Values type: WLZ_VOXELVALUETABLE_GREY.
    Linkcount: 1.
    Background type: WLZ_GREY_UBYTE.
    Background value: 255.
    Values plane bounds: 0 0
    Property list:
    Property type:WLZ_PROPERTY_EMAP
      Linkcount: 1.
      EMAP property type: WLZ_EMAP_PROPERTY_GREY_MODEL.
      Theiler Stage: 15
      Model name: embryo_2_WM_left
      Model version: 1.0
      Filename: embryo_2_WM_left.wlz
      Creation time: Fri Jul 29 17:44:54 2005
      Creation author: Richard Baldock
      Creation machine name: richard-lt.local
      Modification time: Fri Jul 29 17:44:54 2005
      Modification author: richard (richard)
      Comment: Program test

\endverbatim

\par See Also
\ref wlzfacts "WlzFacts(1)"

\par Bugs
None known
\author       richard <Richard.Baldock@hgu.mrc.ac.uk>
\date         Fri Jul 29 12:17:59 2005
\version      MRC HGU $Id$
              $Revision$
              $Name$
\par Copyright:
             1994-2003 Medical Research Council, UK.
              All rights reserved.
\par Address:
              MRC Human Genetics Unit,
              Western General Hospital,
              Edinburgh, EH4 2XU, UK.
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/************************************************************************
*   Copyright  :   1994 Medical Research Council, UK.                   *
*                  All rights reserved.                                 *
*************************************************************************
*   Address    :   MRC Human Genetics Unit,                             *
*                  Western General Hospital,                            *
*                  Edinburgh, EH4 2XU, UK.                              *
*************************************************************************
*   Project    :   MRC HGU Image Processing Utilities			*
*   File       :   WlzChangeEMAPProperty.c				*
*************************************************************************
*   Author Name :  Richard Baldock					*
*   Author Login:  richard@hgu.mrc.ac.uk				*
*   Date        :  Mon Jul 29 16:48:25 2002				*
*   $Revision$					       		*
*   $Name$								*
*   Synopsis    : 							*
*************************************************************************
*   Maintenance :  date - name - comments (Last changes at the top)	*
************************************************************************/

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
	  "Usage:\t%s [-c<time>] [-C<author>] [-f<filename>] "
	  "[-l<location>] [-m<time>] [-M<author>] [-n<name>] "
	  "[-r] [-s#] [-t#] [-T<text>] [-V<version>] [-h] [-v]"
	  "[<input file>]\n"
	  "\tChange the EMAP property of a woolz object\n"
	  "\twriting the new object to standard output\n"
	  "\tNote times and authors are put in from the environment\n"
	  "\tautomatically. Be sure to put in appropriate escapes\n"
	  "\tfor spaces and special characters.\n"
	  "\tOptions are:\n"
	  "\t  -c<time>     Creation time (dd.mm.yyyy)\n"
	  "\t  -C<author>   Creation author\n"
	  "\t  -f<filename> File name associated with object\n"
	  "\t  -l<location> Machine name (or IP address)\n"
	  "\t  -m<time>     Modification time (dd.mm.yyyy)\n"
	  "\t  -M<author>   Modification author\n"
	  "\t  -n<name>     EMAP Model name\n"
	  "\t  -r           Remove the property\n"
	  "\t  -s#          Theiler stage (default 1)\n"
	  "\t  -t#          EMAP property type:\n"
	  "\t               # = %d: grey-level EMAP-model\n"
	  "\t                   %d: EMAP anatomy domain\n"
	  "\t                   %d: other domain (default)\n"
	  "\t  -T<text>     Comment text\n"
	  "\t  -V<version>  Model version\n"
	  "\t  -h           Help - prints this usage message\n"
	  "\t  -v           Verbose operation\n"
	  "",
	  proc_str, WLZ_EMAP_PROPERTY_GREY_MODEL,
	  WLZ_EMAP_PROPERTY_ANATOMY_DOMAIN,
	  WLZ_EMAP_PROPERTY_OTHER_DOMAIN);
  return;
}

int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj;
  FILE		*inFile;
  char 		optList[] = "c:C:f:l:m:M:n:rs:t:T:V:hv";
  int		option;
  WlzEMAPPropertyType	emapType=WLZ_EMAP_PROPERTY_OTHER_DOMAIN;
  int		emapTypeFlg=0;
  int		theilerStage=1, theilerStageFlg=0;
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
    
  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'c':
#ifdef BSD 
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

    case 'l':
      creationMachineName = optarg;
      break;

    case 'm':
#ifdef BSD 
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
      theilerStage = atoi(optarg);
      theilerStageFlg = 1;
      break;

    case 't':
      switch( emapType = (WlzEMAPPropertyType) atoi(optarg) ){

      case WLZ_EMAP_PROPERTY_GREY_MODEL:
      case WLZ_EMAP_PROPERTY_ANATOMY_DOMAIN:
      case WLZ_EMAP_PROPERTY_OTHER_DOMAIN:
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
					      1, NULL, NULL, NULL, NULL, NULL);
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
	if( !theilerStageFlg ){
	  theilerStage = property.emap->theilerStage;
	}
	WlzChangeEMAPProperty(property.emap, emapType, theilerStage,
			      modelName, version, fileName, comment);
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
