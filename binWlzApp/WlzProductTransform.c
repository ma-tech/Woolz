#pragma ident "MRC HGU $Id$"
/************************************************************************
*   Copyright  :   1994 Medical Research Council, UK.                   *
*                  All rights reserved.                                 *
*************************************************************************
*   Address    :   MRC Human Genetics Unit,                             *
*                  Western General Hospital,                            *
*                  Edinburgh, EH4 2XU, UK.                              *
*************************************************************************
*   Project    :   Woolz Library					*
*   File       :   WlzProductTransform.c				*
*************************************************************************
* This module has been copied from the original woolz library and       *
* modified for the public domain distribution. The original authors of  *
* the code and the original file headers and comments are in the        *
* HISTORY file.                                                         *
*************************************************************************
*   Author Name :  Richard Baldock					*
*   Author Login:  richard@hgu.mrc.ac.uk				*
*   Date        :  Mon May  3 15:24:07 1999				*
*   $Revision$								*
*   $Name$								*
*   Synopsis    : generate and write out the product transform of two	*
*		input trabsforms. If 3D go by plane number.		*
*************************************************************************
*   Maintenance :  date - name - comments (Last changes at the top)	*
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>

#include <Wlz.h>
#include <bibFile.h>
#include <Reconstruct.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s [-A] [-R] [-1] [-2] [-ffilename1] [-Ffilename2] [-3] "
	  "[-i] [-I] [-x#] [-y#] [-z#] [-s#] [-a#] [-b#] [-u#] [-v#] [-w#] "
	  "[-h] [-v] [<bibfile1>] [<bibfile2>]\n"
	  "\tFind the product transform for the given bibfiles.\n"
	  "\tBibfile data from stdin is also read and fitted in\n"
	  "\tas per the command line arguments i.e. if filename1\n"
	  "\tor bibfile1 are defined then it will be transform 2\n"
	  "\totherwise transform 1. If both transforms are defined\n"
	  "\ton the command line then stdin is ignored.\n"
	  "\tThe transform entered as arguments assume degree angles\n"
	  "\tand will be applied after transform 1 or 2.\n"
	  "\tOptions are:\n"
	  "\t  -A                calculate absolute transform(default)\n"
	  "\t  -R                calculate relative transform\n"
	  "\t  -1                transform 1 is relative (def: absolute)\n"
	  "\t  -2                transform 2 is relative (def: absolute)\n"
	  "\t  -3                command-line transform is 3D\n"
	  "\t  -a                Rotation about the z-axis.\n"
	  "\t  -b                Rotation about the y-axis.\n"
	  "\t  -i                Invert: reflect about the y-axis.\n"
	  "\t  -I                Inverse: output the inverse of the product.\n"
	  "\t  -s                Scale factor.\n"
	  "\t  -u                Shear strength.\n"
	  "\t  -v                Shear angle in x-y plane.\n"
	  "\t  -w                3D shear angle.\n"
	  "\t  -x                Column (x) translation.\n"
	  "\t  -y                Row (y) translation.\n"
	  "\t  -z                Plane (z) translation.\n"
	  "\t  -h                Help - prints this usage message\n"
	  "\t  -V                verbose operation\n"
	  "",
	  proc_str);
  return;
}

void checkTrans(
  WlzAffineTransform *trans)
{
  double	tol = 1.0e-6;

  if( trans == NULL ){
    return;
  }

  if( !finite(trans->tx) || (fabs(trans->tx) <= tol) ){
    trans->tx = 0.0;
  }
  if( !finite(trans->ty) || (fabs(trans->ty) <= tol) ){
    trans->ty = 0.0;
  }
  if( !finite(trans->tz) || (fabs(trans->tz) <= tol) ){
    trans->tz = 0.0;
  }
  if( !finite(trans->scale) || (fabs(trans->scale - 1.0) <= tol) ){
    trans->scale = 1.0;
  }
  if( !finite(trans->theta) || (fabs(trans->theta) <= tol) ){
    trans->theta = 0.0;
  }
  if( !finite(trans->phi) || (fabs(trans->phi) <= tol) ){
    trans->phi = 0.0;
  }
  if( !finite(trans->alpha) || (fabs(trans->alpha) <= tol) ){
    trans->alpha = 0.0;
  }
  if( !finite(trans->psi) || (fabs(trans->psi) <= tol) ){
    trans->psi = 0.0;
  }
  if( !finite(trans->xsi) || (fabs(trans->xsi) <= tol) ){
    trans->xsi = 0.0;
  }

  WlzAffineTransformMatrixUpdate(trans);
  return;
}


WlzObject *SecListToTransforms(
  RecSectionList	*secList,
  int		relFlg,
  WlzErrorNum	*dstErr)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzObject	*rtnObj=NULL;
  int		i, p, pMin, pMax;
  WlzDomain	domain;
  WlzValues	values;
  RecSection	*sec;
  HGUDlpListItem	*listItem;
  WlzAffineTransform	*recTrans, *tmpTrans;

  /* sort the section list and find the section range */
  RecSecListSort(secList->list, REC_SECMSK_INDEX);

  /* calculate the cumulative transforms */
  listItem = HGUDlpListTail(secList->list);
  if( relFlg ){
    RecSecCumTransfSet(secList->list, listItem);
  }

  /* define the reconstruct transform */
  recTrans = WlzAffineTransformFromPrim(WLZ_TRANSFORM_2D_AFFINE,
					0.0, 0.0, 0.0,
					secList->reconstruction.scale.vtX,
					0.0, 0.0, 0.0, 0.0,
					0.0, 0, NULL);

  /* create the transforms object */
  listItem = HGUDlpListHead(secList->list);
  sec = (RecSection *) HGUDlpListEntryGet(secList->list, listItem);
  pMin = sec->index;
  listItem = HGUDlpListTail(secList->list);
  sec = (RecSection *) HGUDlpListEntryGet(secList->list, listItem);
  pMax = sec->index;
  if( errNum == WLZ_ERR_NONE ){
    if( domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_AFFINE,
				      pMin, pMax, 0, 0, 0, 0,
				      &errNum) ){
      values.core = NULL;
      rtnObj = WlzMakeMain(WLZ_3D_DOMAINOBJ,
			   domain, values, NULL, NULL, &errNum);
    }
  }

  /* now put in the transforms */
  listItem = HGUDlpListHead(secList->list);
  while( listItem ){
    sec = (RecSection *) HGUDlpListEntryGet(secList->list, listItem);
    p = sec->index;
    i = p - rtnObj->domain.p->plane1;
    if( relFlg ){
      if( sec->cumTransform == NULL ){
	rtnObj->domain.p->domains[i].t =
	  WlzAffineTransformFromPrim(WLZ_TRANSFORM_2D_AFFINE,
				     0.0, 0.0, 0.0, 1.0,
				     0.0, 0.0, 0.0, 0.0,
				     0.0, 0, NULL);
      }
      else {
	rtnObj->domain.p->domains[i].t =
	  WlzAssignAffineTransform(sec->cumTransform, NULL);
      }
    }
    else {
      rtnObj->domain.p->domains[i].t =
	WlzAssignAffineTransform(sec->transform, NULL);
    }

    /* apply the reconstruct transform */
    tmpTrans = WlzAffineTransformProduct(rtnObj->domain.p->domains[i].t,
					 recTrans, NULL);
    WlzFreeAffineTransform(rtnObj->domain.p->domains[i].t);
    rtnObj->domain.p->domains[i].t = WlzAssignAffineTransform(tmpTrans, NULL);

    listItem = HGUDlpListNext(secList->list, listItem);
  }

  WlzFreeAffineTransform(recTrans);

  if( dstErr ){
    *dstErr = errNum;
  }
  if( (errNum != WLZ_ERR_NONE) && rtnObj ){
    WlzFreeObj(rtnObj);
    rtnObj = NULL;
  }
  return rtnObj;
}

WlzObject *ReadMAPaintRealignTransforms(
  FILE 		*fp,
  WlzErrorNum	*dstErr)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzObject	*rtnObj=NULL;
  BibFileRecord		*bibfileRecord;
  BibFileField		*bibfileField;
  BibFileError		bibFileErr=BIBFILE_ER_NONE;
  int		i, p, pMin, pMax;
  int		identFlg=0, planeFlg=0;
  WlzDomain	domain;
  WlzValues	values;
  char		*errMsg;
  int		numParsedFields=0;
  double		tx=0.0, ty=0.0, theta=0.0;
  WlzAffineTransform	*inTrans;

  /* read the file to establish the section range
     note must have an Ident record, also check for a repeat
     Ident record which indicates a following set of records
     e.g. if reading from stdin */
  bibFileErr = BibFileRecordRead(&bibfileRecord, &errMsg, fp);
  if( bibFileErr != BIBFILE_ER_NONE ){
    errNum = WLZ_ERR_UNSPECIFIED;
  }
  identFlg = 0;
  while( bibFileErr == BIBFILE_ER_NONE ) 
  {
    /* do nothing until Ident found, quit if found twice */
    if( !identFlg ){
      if( !strncmp(bibfileRecord->name, "Ident", 5) ){
	identFlg = 1;
      }
      
    }
    else {
      /* should check what type of bibfile here */
      if( !strncmp(bibfileRecord->name, "Ident", 5) ){
	BibFileRecordFree(&bibfileRecord);
	break;
      }

      if( !strncmp(bibfileRecord->name, "Section", 7) ){
	/* get the plane number */
	p = atoi(bibfileRecord->id);
	if( planeFlg ){
	  pMin = WLZ_MIN(pMin, p);
	  pMax = WLZ_MAX(pMax, p);
	}
	else {
	  pMin = p;
	  pMax = p;
	  planeFlg = 1;
	}
      }
    }
    BibFileRecordFree(&bibfileRecord);
    bibFileErr = BibFileRecordRead(&bibfileRecord, &errMsg, fp);
  }
  if( (bibFileErr != BIBFILE_ER_NONE) && (bibFileErr != BIBFILE_ER_EOF) ){
    errNum = WLZ_ERR_UNSPECIFIED;
  }
  /* this assumes that the records are at the beginning of the file */
  rewind(fp);

  /* make the tranforms object */
  if( errNum == WLZ_ERR_NONE ){
    if( domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_AFFINE,
				      pMin, pMax, 0, 0, 0, 0,
				      &errNum) ){
      values.core = NULL;
      rtnObj = WlzMakeMain(WLZ_3D_DOMAINOBJ,
			   domain, values, NULL, NULL, &errNum);
    }
  }

  /* now re-read the file, setting the transforms */
  if( errNum == WLZ_ERR_NONE ){
    bibFileErr = BibFileRecordRead(&bibfileRecord, &errMsg, fp);
    if( bibFileErr != BIBFILE_ER_NONE ){
      errNum = WLZ_ERR_UNSPECIFIED;
    }
    identFlg = 0;
    while((bibFileErr == BIBFILE_ER_NONE) &&
	  (errNum == WLZ_ERR_NONE)) 
    {
      /* do nothing until Ident found, quit if found twice */
      if( !identFlg ){
	if( !strncmp(bibfileRecord->name, "Ident", 5) ){
	  identFlg = 1;
	}
      
      }
      else {
	/* should check what type of bibfile here */
	if( !strncmp(bibfileRecord->name, "Ident", 5) ){
	  BibFileRecordFree(&bibfileRecord);
	  break;
	}

	if( !strncmp(bibfileRecord->name, "Section", 7) ){
	  /* get the plane number */
	  p = atoi(bibfileRecord->id);
	  p -= rtnObj->domain.p->plane1;

	  /* parse the record */
	  numParsedFields = BibFileFieldParseFmt
	    (bibfileRecord->field,
	     (void *) &tx, "%lg", "TransformTx",
	     (void *) &ty, "%lg", "TransformTy",
	     (void *) &theta, "%lg", "TransformTheta",
	     NULL);

	  /* make the transform for this record */
	  inTrans = WlzAffineTransformFromPrim(WLZ_TRANSFORM_2D_AFFINE,
					       tx, ty, 0.0, 1.0,
					       theta, 0.0, 0.0, 0.0,
					       0.0, 0, NULL);

	  rtnObj->domain.p->domains[p].t =
	    WlzAssignAffineTransform(inTrans, NULL);

	}
      }
      BibFileRecordFree(&bibfileRecord);
      bibFileErr = BibFileRecordRead(&bibfileRecord, &errMsg, fp);
    }
    if( (bibFileErr != BIBFILE_ER_NONE) && (bibFileErr != BIBFILE_ER_EOF) ){
      errNum = WLZ_ERR_UNSPECIFIED;
    }
  }
  

  if( dstErr ){
    *dstErr = errNum;
  }
  if( (errNum != WLZ_ERR_NONE) && rtnObj ){
    WlzFreeObj(rtnObj);
    rtnObj = NULL;
  }
  return rtnObj;
}

WlzErrorNum WriteMAPaintRealignTransforms(
  FILE		*fp,
  WlzObject	*transformsObj,
  int		relFlg)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  int			i, p;
  BibFileRecord	*record = NULL;
  BibFileField	*field0 = NULL,
		*field1 = NULL,
		*field2 = NULL;
  char		tmpBuf[256], *eMsg;
  char		*tmpS,
		*idxS = NULL,
		*dateS = NULL,
		*hostS = NULL,
		*userS = NULL,
                *fileS = NULL;
  char 		tTypeS[32],
		tTxS[32],
		tTyS[32],
		tTzS[32],
		tScaleS[32],
		tThetaS[32],
		tPhiS[32],
		tAlphaS[32],
		tPsiS[32],
		tXsiS[32],
		tInvertS[32];
  static char	unknownS[] = "unknown";
  char		*typeStr;
  time_t	tmpTime;

  /* run through the transform list writing a bibfile for each
     plane. Use correct plane numbers but do not define a
     file name (since its not known)
  */

  /* an identifier record - who, when, how */
  if( relFlg ){
    typeStr = "relative";
  }
  else {
    typeStr = "absolute";
  }
  if((((field0 = BibFileFieldMakeVa("Text",
				    "MAPaint realignment section file",
				    "Version", "1",
				    "TransformType", typeStr,
				    NULL)) == NULL) ||
      ((record = BibFileRecordMake("Ident", "0", field0)) == NULL)))
  {
    return WLZ_ERR_UNSPECIFIED;
  }
  BibFileRecordWrite(fp, &eMsg, record);
  BibFileRecordFree(&record);
  
  tmpS = getenv("USER");
  (void )sprintf(tmpBuf, "User: %s", tmpS?tmpS:unknownS);
  userS = AlcStrDup(tmpBuf);
  
  tmpTime = time(NULL);
  tmpS = ctime(&tmpTime);
  *(tmpS + strlen(tmpS) - 1) = '\0';
  (void )sprintf(tmpBuf, "Date: %s", tmpS?tmpS:unknownS);
  dateS = AlcStrDup(tmpBuf);
  
  tmpS = getenv("HOST");
  (void )sprintf(tmpBuf, "Host: %s", tmpS?tmpS:unknownS);
  hostS = AlcStrDup(tmpBuf);
  
  (void )sprintf(tmpBuf, "File: %s", unknownS);
  fileS = AlcStrDup(tmpBuf);
  
  if((((field0 = BibFileFieldMakeVa("Text", userS,
				    "Text", dateS,
				    "Text", hostS,
				    "Text", fileS,
				    NULL)) == NULL) ||
      ((record = BibFileRecordMake("Comment", "0", field0)) == NULL)))
  {
    return WLZ_ERR_UNSPECIFIED;
  }
  BibFileRecordWrite(fp, &eMsg, record);
  BibFileRecordFree(&record);
  AlcFree(dateS);
  AlcFree(hostS);
  AlcFree(userS);

  /* now the section records */
  if( transformsObj ){
    for(p=transformsObj->domain.p->plane1, i=0;
	p <= transformsObj->domain.p->lastpl; p++, i++){
      WlzAffineTransform  *transf=transformsObj->domain.p->domains[i].t;
      if( transf == NULL ){
	continue;
      }
      sprintf(tmpBuf, "%d", p);
      idxS = AlcStrDup(tmpBuf);

      (void )sprintf(tTypeS, "%d", transf->type);
      (void )sprintf(tTxS, "%g", transf->tx);
      (void )sprintf(tTyS, "%g", transf->ty);
      (void )sprintf(tTzS, "%g", transf->tz);
      (void )sprintf(tScaleS, "%g", transf->scale);
      (void )sprintf(tThetaS, "%g", transf->theta);
      (void )sprintf(tPhiS, "%g", transf->phi);
      (void )sprintf(tAlphaS, "%g", transf->alpha);
      (void )sprintf(tPsiS, "%g", transf->psi);
      (void )sprintf(tXsiS, "%g", transf->xsi);
      (void )sprintf(tInvertS, "%d", transf->invert);
      field0 = BibFileFieldMakeVa(
	"TransformType", tTypeS,
	"TransformTx", tTxS,
	"TransformTy", tTyS,
	"TransformTz", tTzS,
	"TransformScale", tScaleS,
	"TransformTheta", tThetaS,
	"TransformPhi", tPhiS,
	"TransformAlpha", tAlphaS,
	"TransformPsi", tPsiS,
	"TransformXsi", tXsiS,
	"TransformInvert", tInvertS,
	NULL);

      record = BibFileRecordMake("Section", idxS, field0);
      BibFileRecordWrite(fp, &eMsg, record);
      BibFileRecordFree(&record);
      AlcFree(idxS);
    }
  }

  return errNum;
}
 
int main(int	argc,
	 char	**argv)
{
  FILE		*inFile1=NULL, *inFile2=NULL;
  char 		optList[] = "AR123f:F:iIx:y:z:s:a:b:u:v:w:hV";
  int		option;
  int		trans1RelFlg=0, trans2RelFlg=0, trans3DFlg=0, outputRelFlg=0;
  int		invertFlg=0, inverseFlg=0, threeDFlg=0;
  WlzAffineTransform	*cmdLineTrans=NULL;
  double	tx, ty, tz, scale, theta, phi, alpha, psi, xsi;
  int		verboseFlg=0;
  RecSectionList	recSecList1, recSecList2;
  RecSectionList	*secList1=&recSecList1, *secList2=&recSecList2;
  HGUDlpListItem	*listItem, *listItem1;
  RecSection	*sec, *sec1;
  char		*errMsg = NULL;
  RecError	errFlg = REC_ERR_NONE;
  int		numSec1, numSec2;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzObject	*transformsObj1, *transformsObj2;
  int		p, i, ip;
  WlzAffineTransform	*newTransform, *tmpTrans, *tmp1Trans, *reconTrans;

  /* default transform */
  tx = 0.0;
  ty = 0.0;
  tx = 0.0;
  scale = 1.0;
  theta = 0.0;
  phi = 0.0;
  alpha = 0.0;
  psi = 0.0;
  xsi = 0.0;

  /* parse the command line */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'A':
      outputRelFlg = 0;
      break;

    case 'R':
      outputRelFlg = 1;
      break;

    case '1':
      trans1RelFlg = 1;
      break;

    case '2':
      trans2RelFlg = 1;
      break;

    case '3':
      trans3DFlg = 1;
      break;

    case 'f':
      if( (optarg != NULL) && (strlen(optarg) > 0) ){
	if( (inFile1 = fopen(optarg, "r")) == NULL ){
	  fprintf(stderr, "%s: failed to open input file %s\n",
		  argv[0], optarg);
	  usage(argv[0]);
	  return 1;
	}
      }
      else {
	usage(argv[0]);
	return 1;
      }
      break;

    case 'F':
      if( (optarg != NULL) && (strlen(optarg) > 0) ){
	if( (inFile2 = fopen(optarg, "r")) == NULL ){
	  fprintf(stderr, "%s: failed to open input file %s\n",
		  argv[0], optarg);
	  usage(argv[0]);
	  return 1;
	}
      }
      else {
	usage(argv[0]);
	return 1;
      }
      break;

    case 'i':
      invertFlg = 1;
      break;

    case 'I':
      inverseFlg = 1;
      break;

    case 'a':
      if(sscanf(optarg, "%lg", &theta) != 1)
      {
	usage(argv[0]);
	return 1;
      }
      break;

    case 'b':
      if(sscanf(optarg, "%lg", &phi) != 1)
      {
	usage(argv[0]);
	return 1;
      }
      break;

    case 's':
      if(sscanf(optarg, "%lg", &scale) != 1)
      {
	usage(argv[0]);
	return 1;
      }
      break;

    case 'u':
      if(sscanf(optarg, "%lg", &alpha) != 1)
      {
	usage(argv[0]);
	return 1;
      }
      break;

    case 'v':
      if(sscanf(optarg, "%lg", &psi) != 1)
      {
	usage(argv[0]);
	return 1;
      }
      break;

    case 'w':
      if(sscanf(optarg, "%lg", &xsi) != 1)
      {
	usage(argv[0]);
	return 1;
      }
      break;

    case 'x':
      if(sscanf(optarg, "%lg", &tx) != 1)
      {
	usage(argv[0]);
	return 1;
      }
      break;

     case 'y':
      if(sscanf(optarg, "%lg", &ty) != 1)
      {
	usage(argv[0]);
	return 1;
      }
      break;

    case 'z':
      if(sscanf(optarg, "%lg", &tz) != 1)
      {
	usage(argv[0]);
	return 1;
      }
      break;

    case 'V':
      verboseFlg = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return 1;

    }
  }

  /* create the command line transform */
  if( threeDFlg ){
    cmdLineTrans = WlzAffineTransformFromPrim(WLZ_TRANSFORM_3D_AFFINE,
					      tx, ty, tz, scale, theta, phi,
					      alpha, psi, xsi,
					      invertFlg, &errNum);
  }
  else {
     cmdLineTrans = WlzAffineTransformFromPrim(WLZ_TRANSFORM_2D_AFFINE,
					      tx, ty, tz, scale, theta, phi,
					      alpha, psi, xsi,
					      invertFlg, &errNum);
  }

  /* check for files on the command line */
  if( (inFile1 == NULL) && (optind < argc) ){
    if( (inFile1 = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
    optind++;
  }
  if( (inFile2 == NULL) && (optind < argc) ){
    if( (inFile2 = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
    optind++;
  }

  /* now read the transforms, in each case the logic is that
     transform 1 is applied then transform 2 then the command line
     transform. If transform 1 is a reconstruct bib file then the
     output will include the file information. Note in this case
     the output transform will be absolute unless specified on the
     command line */

  /* transform 1 - try to read as a Reconstruct file, on failure
     rewind and read as an MAPaint file */
  if( inFile1 != NULL ){
    errFlg = RecFileSecListRead(secList1, &numSec1, inFile1, &errMsg);
    if( errFlg == REC_ERR_NONE ){
      /* check transform attribute */
      switch( secList1->attributes.trMode ){
      case REC_TRMODE_REL:
	trans1RelFlg = 1;
	break;
      case REC_TRMODE_ABS:
	trans1RelFlg = 0;
	break;
      }
      transformsObj1 = SecListToTransforms(secList1, trans1RelFlg, &errNum);
    }
    else {
      secList1 = NULL;
      rewind(inFile1);
      transformsObj1 = ReadMAPaintRealignTransforms(inFile1, &errNum);
    }
  }
  else {
    fprintf(stderr, "%s: a file for transform 1 must be defined\n", argv[0]);
    usage(argv[0]);
    return 1;
  }

  /* transform 2 - try to read as a Reconstruct file, on failure
     rewind and read as an MAPaint file */
  if( inFile2 != NULL ){
    errFlg = RecFileSecListRead(secList2, &numSec2, inFile2, &errMsg);
    if( errFlg == REC_ERR_NONE ){
      /* check transform attribute */
      switch( secList2->attributes.trMode ){
      case REC_TRMODE_REL:
	trans2RelFlg = 1;
	break;
      case REC_TRMODE_ABS:
	trans2RelFlg = 0;
	break;
      }
      transformsObj2 = SecListToTransforms(secList2, trans2RelFlg, &errNum);
    }
    else {
      secList2 = NULL;
      rewind(inFile2);
      transformsObj2 = ReadMAPaintRealignTransforms(inFile2, &errNum);
    }
  }
  else {
    transformsObj2 = NULL;
  }

  /* form the plane by plane product of the transforms,
     note the output will be of the planes in transform 1 */
  for(p=transformsObj1->domain.p->plane1, i=0;
      p <= transformsObj1->domain.p->lastpl; p++, i++){

    /* get the existing transform */
    if( transformsObj1->domain.p->domains[i].t ){
      newTransform =
	WlzAssignAffineTransform(transformsObj1->domain.p->domains[i].t,
				 NULL);
    }
    else {
      newTransform = NULL;
    }

    /* check for transform 2 */
    if(newTransform && transformsObj2 &&
       (p >= transformsObj2->domain.p->plane1) &&
       (p <= transformsObj2->domain.p->lastpl) &&
       (transformsObj2->domain.p->
	domains[p-transformsObj2->domain.p->plane1].t)){
      tmpTrans =
	WlzAffineTransformProduct(newTransform,
				  transformsObj2->domain.p->
				  domains[p-transformsObj2->domain.p->plane1].t,
				  NULL);
      WlzFreeAffineTransform(newTransform);
      newTransform = tmpTrans;
    }

    /* apply the command line transform */
    if( newTransform && cmdLineTrans ){
      tmpTrans = WlzAffineTransformProduct(newTransform, cmdLineTrans, NULL);
      WlzFreeAffineTransform(newTransform);
      newTransform = tmpTrans;
    }

    /* put it back into transformsObj1 */
    if( transformsObj1->domain.p->domains[i].t ){
      WlzFreeAffineTransform(transformsObj1->domain.p->domains[i].t);
    }
    transformsObj1->domain.p->domains[i].t =
      WlzAssignAffineTransform(newTransform, NULL);
  }

  /* if there is a section list in Reconstruct form
     then take out the reconstruct transform */
  if( secList1 ){
    /* take out the reconstruct transform first */
    reconTrans = WlzAffineTransformFromPrim
      (WLZ_TRANSFORM_2D_AFFINE, 0.0, 0.0, 0.0,
       1.0 / secList1->reconstruction.scale.vtX,
       0.0, 0.0, 0.0, 0.0, 0.0, 0, NULL);
    listItem = HGUDlpListHead(secList1->list);
    while( listItem ){
      sec = (RecSection *) HGUDlpListEntryGet(secList1->list, listItem);
      if( sec ){
	p = sec->index;
	i = p - transformsObj1->domain.p->plane1;
	if( transformsObj1->domain.p->domains[i].t ){
	  tmpTrans = WlzAffineTransformProduct
	    (transformsObj1->domain.p->domains[i].t, reconTrans, NULL);
	  WlzFreeAffineTransform(transformsObj1->domain.p->domains[i].t);
	  transformsObj1->domain.p->domains[i].t =
	    WlzAssignAffineTransform(tmpTrans, NULL);
	}
      }
      listItem = HGUDlpListNext(secList1->list, listItem);
    }

  }

  /* check for relative flag */
  if( outputRelFlg ){
    /* use the section list if available */
    if( secList1 ){

      /* sort out the relative transform */
      listItem = HGUDlpListTail(secList1->list);
      while( listItem != HGUDlpListHead(secList1->list) ){
	sec = (RecSection *) HGUDlpListEntryGet(secList1->list, listItem);
	if( sec ){
	  p = sec->index;
	  i = p - transformsObj1->domain.p->plane1;
	  if( strncmp(sec->imageFile, "empty", 5) == 0 ){
	    WlzFreeAffineTransform(transformsObj1->domain.p->domains[i].t);
	    transformsObj1->domain.p->domains[i].t =
	      WlzAffineTransformFromPrim(WLZ_TRANSFORM_2D_AFFINE,
					 0.0, 0.0, 0.0, 1.0,
					 0.0, 0.0, 0.0, 0.0,
					 0.0, 0, NULL);
	  }
	  else {
	    /* walk backwards to previous non-empty section */
	    listItem1 = HGUDlpListPrev(secList1->list, listItem);
	    while( listItem1 != HGUDlpListHead(secList1->list) ){
	      sec1 = (RecSection *) HGUDlpListEntryGet(secList1->list, listItem1);
	      if( strncmp(sec1->imageFile, "empty", 5) ){
		break;
	      }
	      listItem1 = HGUDlpListPrev(secList1->list, listItem1);
	    }
	    ip = sec1->index - transformsObj1->domain.p->plane1;
	    tmpTrans =
	      WlzAffineTransformInverse(transformsObj1->domain.p->domains[ip].t,
					NULL);
	    tmp1Trans = WlzAffineTransformProduct(
	      transformsObj1->domain.p->domains[i].t,
	      tmpTrans, NULL);
	    checkTrans(tmp1Trans);
	    WlzFreeAffineTransform(transformsObj1->domain.p->domains[i].t);
	    WlzFreeAffineTransform(tmpTrans);
	    transformsObj1->domain.p->domains[i].t = tmp1Trans;
	  }
	}
	listItem = HGUDlpListPrev(secList1->list, listItem);
      }
    }
    else {
      i = transformsObj1->domain.p->lastpl - transformsObj1->domain.p->plane1;
      ip = i;
      for(p = transformsObj1->domain.p->lastpl;
	  p >= transformsObj1->domain.p->plane1; p--, i--){
	if( i > ip ){
	  continue;
	}
	ip--;
	if(transformsObj1->domain.p->domains[i].t){
	  /* jump over any NULL transforms */
	  while((ip >= 0) &&
		(transformsObj1->domain.p->domains[ip].t == NULL)){
	      ip--;
	  }

	  /* push the transform down the line */
	  if(ip >= 0){
	    tmpTrans =
	      WlzAffineTransformInverse(transformsObj1->domain.p->domains[ip].t,
					NULL);
	    tmp1Trans =
	      WlzAffineTransformProduct(transformsObj1->domain.p->domains[i].t,
					tmpTrans,
					NULL);
	    checkTrans(tmp1Trans);
	    WlzFreeAffineTransform(transformsObj1->domain.p->domains[i].t);
	    WlzFreeAffineTransform(tmpTrans);
	    transformsObj1->domain.p->domains[i].t =
	      WlzAssignAffineTransform(tmp1Trans, NULL);
	  }
	}
      }
    }
  }

  /* if there is a section list one then output in Reconstruct
     bibfile format else MAPaint style */
  if( secList1 ){
    /* set the relative/absolute attribute */
    if( outputRelFlg ){
      secList1->attributes.trMode = REC_TRMODE_REL;
    }
    else {
      secList1->attributes.trMode = REC_TRMODE_ABS;
    }

    /* put back the transforms */
    listItem = HGUDlpListHead(secList1->list);
    while( listItem ){
      sec = (RecSection *) HGUDlpListEntryGet(secList1->list, listItem);
      p = sec->index;
      i = p - transformsObj1->domain.p->plane1;

      /* check for NULL transform */
      if( transformsObj1->domain.p->domains[i].t == NULL ){
	transformsObj1->domain.p->domains[i].t =
	  WlzAffineTransformFromPrim(WLZ_TRANSFORM_2D_AFFINE,
				     0.0, 0.0, 0.0, 1.0,
				     0.0, 0.0, 0.0, 0.0,
				     0.0, 0, NULL);
      }

      sec->transform = 
	transformsObj1->domain.p->domains[i].t;
      
      listItem = HGUDlpListNext(secList1->list, listItem);
    }
    RecFileSecListWrite(stdout, secList1, numSec1, &errMsg);
  }
  else {
    WriteMAPaintRealignTransforms(stdout, transformsObj1, outputRelFlg);
  }

  return 0;
}
