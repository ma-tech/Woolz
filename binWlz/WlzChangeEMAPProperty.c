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
*   File       :   WlzChangeEMAPProperty.c				*
*************************************************************************
*   Author Name :  Richard Baldock					*
*   Author Login:  richard@hgu.mrc.ac.uk				*
*   Date        :  Mon Jul 29 16:48:25 2002				*
*   $Revision$								*
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
	obj->plist = AlcDLPListNew(NULL);
      }
      if( obj->plist ){
	property = WlzGetProperty(obj->plist, WLZ_PROPERTY_EMAP, NULL);
	if( property.core == NULL ){
	  property.emap = WlzMakeEMAPProperty(WLZ_EMAP_PROPERTY_GREY_MODEL,
					      1, NULL, NULL, NULL, NULL, NULL);
	  WlzAssignProperty(property, NULL);
	  AlcDLPListEntryAppend(obj->plist, NULL, (void *) property.emap,
				WlzFreePropertyListEntry);
	}
      }
      if( removeFlg ){
	WlzRemoveProperty(obj->plist, property);
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
		

