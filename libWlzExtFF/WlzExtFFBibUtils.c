#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzExtFFBibUtils.c
* \author       Richard Baldock <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Fri Jan  4 12:31:36 2002
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Copyright:
*               1994-2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        Woolz IO bibfile utility functions to write
*		and parse records written in bibfile format.
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <WlzExtFF.h>

/*! 
* \ingroup      WlzExtFF
* \brief        Write the section view parameters in the bibtex style
*		record using the bibFile library.
* \return       woolz error number
* \param    fp	FILE pointer opened for writing
* \param    recordName	record name
* \param    wlzViewStr	woolz 3D view structure
* \par      Source:
*                WlzExtFFBibUtils.c
*/
WlzErrorNum WlzEffBibWrite3DSectionViewParamsRecord(
  FILE			*fp,
  char			*recordName,
  WlzThreeDViewStruct	*wlzViewStr)
{
  char	fixedPointStr[64];
  char	distanceStr[64];
  char	pitchStr[64];
  char	yawStr[64];
  char	rollStr[64];
  char	scaleStr[64];
  char	upVectorStr[64];
  char	viewModeStr[64];
  BibFileRecord	*bibfileRecord;
  BibFileError	bibFileErr;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* generate strings for the bibfile entry */
  sprintf(fixedPointStr, "%.1f, %.1f, %.1f", wlzViewStr->fixed.vtX,
	  wlzViewStr->fixed.vtY, wlzViewStr->fixed.vtZ);
  sprintf(distanceStr, "%.01f", wlzViewStr->dist);
  sprintf(pitchStr, "%.2f", wlzViewStr->phi * 180.0 / WLZ_M_PI);
  sprintf(yawStr, "%.2f", wlzViewStr->theta * 180.0 / WLZ_M_PI);
  sprintf(rollStr, "%.2f", wlzViewStr->zeta * 180.0 / WLZ_M_PI);
  sprintf(scaleStr, "%.2f", wlzViewStr->scale);
  sprintf(upVectorStr, "%g, %g, %g", wlzViewStr->up.vtX,
	  wlzViewStr->up.vtY, wlzViewStr->up.vtZ);

  /* view mode */
  switch( wlzViewStr->view_mode ){
  case WLZ_UP_IS_UP_MODE:
  default:
    sprintf(viewModeStr, "up-is-up");
    break;

  case WLZ_STATUE_MODE:
    sprintf(viewModeStr, "statue");
    break;

  case WLZ_ZETA_MODE:
    sprintf(viewModeStr, "absolute");
    break;
  }

  /* create the bibfile record */
  bibfileRecord = 
    BibFileRecordMake(recordName, "0",
		      BibFileFieldMakeVa("FixedPoint",	fixedPointStr,
					 "Distance",	distanceStr,
					 "Pitch",	pitchStr,
					 "Yaw",		yawStr,
					 "Roll",	rollStr,
					 "Scale",	scaleStr,
					 "UpVector",	upVectorStr,
					 "ViewMode",	viewModeStr,
					 NULL));

  /* write the bibfile and release resources */
  bibFileErr = BibFileRecordWrite(fp, NULL, bibfileRecord);
  BibFileRecordFree(&bibfileRecord);

  /* check the bibfile write error */
  if( bibFileErr != BIBFILE_ER_NONE ){
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }

  return errNum;
}

/*! 
* \ingroup      WlzExtFF
* \brief        Parse the bibfile record for a woolz 3D view structure
*		setting paramters in the given view structure.
* \return       woolz error number
* \param    bibfileRecord	bibfile record as read using libbibFile
* \param    wlzViewStr	woolz 3D view structure already allocated
* \par      Source:
*                WlzExtFFBibUtils.c
*/
WlzErrorNum WlzEffBibParse3DSectionViewParamsRecord(
  BibFileRecord		*bibfileRecord,
  WlzThreeDViewStruct	*wlzViewStr)
{
  WlzDVertex3	fixedPointVtx;
  double	distance, pitch, yaw, roll, scale;
  WlzDVertex3	upVectorVtx;
  WlzThreeDViewMode	viewMode;
  char		viewModeStr[64];
  int		numParsedFields=0;
  BibFileField	*bibfileField;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check inputs */
  if((bibfileRecord == NULL) ||
     (wlzViewStr == NULL) ){
    errNum = WLZ_ERR_PARAM_NULL;
  }

  /* parse the record */
  if( errNum == WLZ_ERR_NONE ){
    numParsedFields = BibFileFieldParseFmt
      (bibfileRecord->field,
       (void *) &(fixedPointVtx.vtX), "%lg ,%*lg ,%*lg", "FixedPoint",
       (void *) &(fixedPointVtx.vtY), "%*lg ,%lg ,%*lg", "FixedPoint",
       (void *) &(fixedPointVtx.vtZ), "%*lg ,%*lg ,%lg", "FixedPoint",
       (void *) &distance, "%lg", "Distance",
       (void *) &pitch, "%lg", "Pitch",
       (void *) &yaw, "%lg", "Yaw",
       (void *) &roll, "%lg", "Roll",
       (void *) &scale, "%lg", "Scale",
       (void *) &(upVectorVtx.vtX), "%lg ,%*lg ,%*lg", "UpVector",
       (void *) &(upVectorVtx.vtY), "%*lg ,%lg ,%*lg", "UpVector",
       (void *) &(upVectorVtx.vtZ), "%*lg ,%*lg ,%lg", "UpVector",
       (void *) &(viewModeStr[0]), "%s", "ViewMode",
       NULL);
    if( numParsedFields < 12 ){
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }

  /* doesn't read the view mode correctly - ask Bill */
  if( errNum == WLZ_ERR_NONE ){
    bibfileField = bibfileRecord->field;
    errNum = WLZ_ERR_READ_INCOMPLETE;
    while( bibfileField ){
      if( strncmp(bibfileField->name, "ViewMode", 8) == 0 ){
	strcpy(viewModeStr, bibfileField->value);
	errNum = WLZ_ERR_NONE;
	break;
      }
      bibfileField = bibfileField->next;
    }
  }

  /* set the view structure */
  if( errNum == WLZ_ERR_NONE ){
    pitch *= (WLZ_M_PI/180.0);
    yaw *= (WLZ_M_PI/180.0);
    roll *= (WLZ_M_PI/180.0);
    if( !strncmp(viewModeStr, "up-is-up", 8) ){
      viewMode = WLZ_UP_IS_UP_MODE;
    }
    else if( !strncmp(viewModeStr, "statue", 8) ){
      viewMode = WLZ_STATUE_MODE;
    }
    else {
      viewMode = WLZ_ZETA_MODE;
    }
    wlzViewStr->fixed = fixedPointVtx;
    wlzViewStr->dist = distance;
    wlzViewStr->theta = yaw;
    wlzViewStr->phi = pitch;
    wlzViewStr->zeta = roll;
    wlzViewStr->up = upVectorVtx;
    wlzViewStr->scale = scale;
    wlzViewStr->view_mode = viewMode;
  }

  return errNum;
}

/*! 
* \ingroup      WlzExtFF
* \brief        Write the details of a warp transformation parameters.
*		Note this does not write the warp transformation itself.
* \return       woolz error number
* \param    fp	FILE pointer opened for writing
* \param    recordName	record name
* \param    basisFnType	interpolation or basis function type
* \param    meshMthd	method used to define the mesh
* \param    meshMinDst	mesh minimum distance
* \param    meshMaxDst	mesh maximum distance
* \par      Source:
*                WlzExtFFBibUtils.c
*/
WlzErrorNum WlzEffBibWriteWarpTransformParamsRecord(
  FILE			*fp,
  char			*recordName,
  WlzBasisFnType	basisFnType,
  WlzMeshGenMethod	meshMthd,
  int			meshMinDst,
  int	 		meshMaxDst)
{
  char	basisFnTypeStr[32];
  char	meshMthdStr[32];
  char	meshMinDstStr[16];
  char	meshMaxDstStr[16];
  BibFileRecord	*bibfileRecord;
  BibFileError	bibFileErr;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* generate strings for the bibfile entry */
  sprintf(meshMinDstStr, "%d", meshMinDst);
  sprintf(meshMaxDstStr, "%d", meshMaxDst);

  /*  types */
  switch( basisFnType ){
  case WLZ_BASISFN_TPS:
  default:
    sprintf(basisFnTypeStr, "Thin-plate-spline");
    break;

  case WLZ_BASISFN_GAUSS:
    sprintf(basisFnTypeStr, "Gaussian");
    break;

  case WLZ_BASISFN_POLY:
    sprintf(basisFnTypeStr, "Polynomial");
    break;

  case WLZ_BASISFN_MQ:
    sprintf(basisFnTypeStr, "Multiquadric");
    break;

  case WLZ_BASISFN_CONF_POLY:
    sprintf(basisFnTypeStr, "Conformal-polynomial");
    break;
  }

  switch( meshMthd ){
  case WLZ_MESH_GENMETHOD_BLOCK:
  default:
    sprintf(meshMthdStr, "Block");
    break;

  case WLZ_MESH_GENMETHOD_GRADIENT:
    sprintf(meshMthdStr, "Gradient");
    break;
  }

  /* create the bibfile record */
  bibfileRecord = 
    BibFileRecordMake(recordName, "0",
		      BibFileFieldMakeVa("BasisFnType",		basisFnTypeStr,
					 "MeshGenMethod",	meshMthdStr,
					 "MeshMinDist",		meshMinDstStr,
					 "MeshMaxDist",		meshMaxDstStr,
					 NULL));

  /* write the bibfile and release resources */
  bibFileErr = BibFileRecordWrite(fp, NULL, bibfileRecord);
  BibFileRecordFree(&bibfileRecord);

  /* check the bibfile write error */
  if( bibFileErr != BIBFILE_ER_NONE ){
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }

  return errNum;
}

/*! 
* \ingroup      WlzExtFF
* \brief        Parse the bibfile record for the warp transformation
*		parameter settings.
* \return       Woolz error number
* \param    bibfileRecord	bibfile record as read by libbibFile
* \param    basisFnType	return for basis function type
* \param    meshMthd	return for mesh method type
* \param    meshMinDst	return for mesh minimum distance
* \param    meshMaxDst	return for mesh maximum distance
* \par      Source:
*                WlzExtFFBibUtils.c
*/
WlzErrorNum WlzEffBibParseWarpTransformParamsRecord(
  BibFileRecord		*bibfileRecord,
  WlzBasisFnType	*basisFnType,
  WlzMeshGenMethod	*meshMthd,
  int			*meshMinDst,
  int	 		*meshMaxDst)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  char		basisFnTypeStr[32];
  char		meshMthdStr[32];
  int		numParsedFields=0;
  BibFileField	*bibfileField;

  /* check inputs */
  if((bibfileRecord == NULL) ||
     (basisFnType == NULL) ||
     (meshMthd == NULL) ||
     (meshMinDst == NULL) ||
     (meshMaxDst == NULL) ){
    errNum = WLZ_ERR_PARAM_NULL;
  }

  /* parse the record */
  if( errNum == WLZ_ERR_NONE ){
    numParsedFields = BibFileFieldParseFmt
      (bibfileRecord->field,
       (void *) &(basisFnTypeStr[0]), "%s", "BasisFnType",
       (void *) &(meshMthdStr[0]), "%s", "MeshGenMethod",
       (void *) meshMinDst, "%d", "MeshMinDist",
       (void *) meshMaxDst, "%d", "MeshMaxDist",
       NULL);
    if( numParsedFields < 4 ){
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }

  /* doesn't read the strings correctly - ask Bill */
  if( errNum == WLZ_ERR_NONE ){
    bibfileField = bibfileRecord->field;
    errNum = WLZ_ERR_READ_INCOMPLETE;
    numParsedFields = 0;
    while( bibfileField ){
      if( strncmp(bibfileField->name, "BasisFnType", 11) == 0 ){
	strcpy(basisFnTypeStr, bibfileField->value);
	numParsedFields++;
      }
      if( strncmp(bibfileField->name, "MeshGenMethod", 13) == 0 ){
	strcpy(meshMthdStr, bibfileField->value);
	numParsedFields++;
      }
      bibfileField = bibfileField->next;
    }
    if( numParsedFields < 2 ){
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    else {
      errNum = WLZ_ERR_NONE;
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    if( strncmp(basisFnTypeStr, "Thin-plate-spline", 17) == 0 ){
      *basisFnType = WLZ_BASISFN_TPS;
    }
    else if( strncmp(basisFnTypeStr, "Gaussian", 8) == 0 ){
      *basisFnType = WLZ_BASISFN_GAUSS;
    }
    else if( strncmp(basisFnTypeStr, "Polynomial", 10) == 0 ){
      *basisFnType = WLZ_BASISFN_POLY;
    }
    else if( strncmp(basisFnTypeStr, "Multiquadric", 12) == 0 ){
      *basisFnType = WLZ_BASISFN_MQ;
    }
    else if( strncmp(basisFnTypeStr, "Conformal-polynomial", 20) == 0 ){
      *basisFnType = WLZ_BASISFN_CONF_POLY;
    }
    else {
      *basisFnType = WLZ_BASISFN_TPS;
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    if( strncmp(meshMthdStr, "Block", 5) == 0 ){
      *meshMthd = WLZ_MESH_GENMETHOD_BLOCK;
    }
    else if( strncmp(meshMthdStr, "Gradient", 8) == 0 ){
      *meshMthd = WLZ_MESH_GENMETHOD_GRADIENT;
    }
    else {
      *meshMthd = WLZ_MESH_GENMETHOD_BLOCK;
    }
  }

  return errNum;
}

/*! 
* \ingroup      WlzExtFF
* \brief        Write a tie-point record in bibtex style.
* \return       Woolz error number
* \param    fp	FILE pointer opened for writing
* \param    recordName	record name
* \param    index	index of the given tie points
* \param    dstVtx	destination vertex
* \param    srcVtx	source vertex
* \par      Source:
*                WlzExtFFBibUtils.c
*/
WlzErrorNum WlzEffBibWriteTiePointVtxsRecord(
  FILE		*fp,
  char		*recordName,
  int		index,
  WlzDVertex3	dstVtx,
  WlzDVertex3	srcVtx)
{
  char	indexStr[16];
  char	dstVtxStr[64];
  char	srcVtxStr[64];
  BibFileRecord	*bibfileRecord;
  BibFileError	bibFileErr;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* generate strings for the bibfile entry */
  sprintf(indexStr, "%d", index);
  sprintf(dstVtxStr, "%f, %f, %f", dstVtx.vtX,
	  dstVtx.vtY, dstVtx.vtZ);
  sprintf(srcVtxStr, "%f, %f, %f", srcVtx.vtX,
	  srcVtx.vtY, srcVtx.vtZ);

  /* create the bibfile record */
  bibfileRecord = 
    BibFileRecordMake(recordName, "0",
		      BibFileFieldMakeVa("Index",	indexStr,
					 "DstVtx",	dstVtxStr,
					 "SrcVtx",	srcVtxStr,
					 NULL));

  /* write the bibfile and release resources */
  bibFileErr = BibFileRecordWrite(fp, NULL, bibfileRecord);
  BibFileRecordFree(&bibfileRecord);

  /* check the bibfile write error */
  if( bibFileErr != BIBFILE_ER_NONE ){
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }

  return errNum;
}

/*! 
* \ingroup      WlzExtFF
* \brief        Parse the bibfile record for a tie-point pair of vertices.
*		The record should also include an index so that tie-point
*		ordering is preserved.
* \return       Woolz error number
* \param    bibfileRecord	bibfile record as read in by libbibFile
* \param    index	return for tie-point pair index within list
* \param    dstVtx	return for destination vertex
* \param    srcVtx	return for source vertex
* \par      Source:
*                WlzExtFFBibUtils.c
*/
WlzErrorNum WlzEffBibParseTiePointVtxsRecord(
  BibFileRecord		*bibfileRecord,
  int		*index,
  WlzDVertex3	*dstVtx,
  WlzDVertex3	*srcVtx)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  int		numParsedFields=0;

  /* check inputs */
  if((bibfileRecord == NULL) ||
     (index == NULL) ||
     (dstVtx == NULL) ||
     (srcVtx == NULL) ){
    errNum = WLZ_ERR_PARAM_NULL;
  }

  /* parse the record */
  if( errNum == WLZ_ERR_NONE ){
    numParsedFields = BibFileFieldParseFmt
      (bibfileRecord->field,
       (void *) index, "%d", "Index",
       (void *) &(dstVtx->vtX), "%lg ,%*lg ,%*lg", "DstVtx",
       (void *) &(dstVtx->vtY), "%*lg ,%lg ,%*lg", "DstVtx",
       (void *) &(dstVtx->vtZ), "%*lg ,%*lg ,%lg", "DstVtx",
       (void *) &(srcVtx->vtX), "%lg ,%*lg ,%*lg", "SrcVtx",
       (void *) &(srcVtx->vtY), "%*lg ,%lg ,%*lg", "SrcVtx",
       (void *) &(srcVtx->vtZ), "%*lg ,%*lg ,%lg", "SrcVtx",
       NULL);
    if( numParsedFields < 7 ){
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }

  return errNum;
}

/*! 
* \ingroup      WlzExtFF
* \brief        Write a bibfile record to record the file name and file type.
* \return       woolz error number
* \param    fp	FILE pointer
* \param    recordName	record name
* \param    fileName	file name to be recorded
* \param    fileType	file type e.g. image format
* \par      Source:
*                WlzExtFFBibUtils.c
*/
WlzErrorNum WlzEffBibWriteFileRecord(
  FILE		*fp,
  char		*recordName,
  char		*fileName,
  WlzEffFormat	fileType)
{
  BibFileRecord	*bibfileRecord;
  BibFileError	bibFileErr;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  const char	*fileTypeStr;

  if( strlen(fileName) > 511 ){
    return WLZ_ERR_PARAM_DATA;
  }

  /* create the bibfile record */
  fileTypeStr = WlzEffStringFromFormat(fileType, NULL);
  bibfileRecord = 
    BibFileRecordMake(recordName, "0",
		      BibFileFieldMakeVa("File", fileName,
					 "Type", fileTypeStr,
					 NULL));

  /* write the bibfile and release resources */
  bibFileErr = BibFileRecordWrite(fp, NULL, bibfileRecord);
  BibFileRecordFree(&bibfileRecord);

  /* check the bibfile write error */
  if( bibFileErr != BIBFILE_ER_NONE ){
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }

  return errNum;
}

/*! 
* \ingroup      WlzExtFF
* \brief        Parse a bibfile record for a filename and type.
* \return       woolz error
* \param    	bibfileRecord	bibfile record as read in by libbibFile
* \param    	dstIndex	Destination pointer for integer parsed from
*				the record ID, may be NULL.
* \param    	dstFileName	Destination pointer for file name, may be NULL.
* \param    	dstFileType	Destination pointer for the file type, may be
*				NULL.
* \par      Source:
*                WlzExtFFBibUtils.c
*/
WlzErrorNum WlzEffBibParseFileRecord(
  BibFileRecord		*bibfileRecord,
  int			*dstIndex,
  char			**dstFileName,
  WlzEffFormat		*dstFileType)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  int		numParsedFields=0;
  char		fileNameBuf[512];
  char		fileTypeBuf[32];
  BibFileField	*bibfileField;

  /* check inputs */
  if(bibfileRecord == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  /* parse the record */
  if( errNum == WLZ_ERR_NONE ){
    numParsedFields = BibFileFieldParseFmt
      (bibfileRecord->field,
       (void *) &(fileNameBuf[0]), "%s", "File",
       (void *) &(fileTypeBuf[0]), "%s", "Type",
       NULL);
    if( numParsedFields < 2 ){
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  /* doesn't read the strings correctly - ask Bill */
  if( errNum == WLZ_ERR_NONE ){
    bibfileField = bibfileRecord->field;
    errNum = WLZ_ERR_READ_INCOMPLETE;
    numParsedFields = 0;
    while( bibfileField ){
      if( strncmp(bibfileField->name, "File", 4) == 0 ){
	strcpy(fileNameBuf, bibfileField->value);
	numParsedFields++;
      }
      if( strncmp(bibfileField->name, "Type", 4) == 0 ){
	strcpy(fileTypeBuf, bibfileField->value);
	numParsedFields++;
      }
      bibfileField = bibfileField->next;
    }
    if( numParsedFields < 2 ){
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    else {
      errNum = WLZ_ERR_NONE;
    }
  }
  if(dstIndex)
  {
    *dstIndex = 0;
    sscanf(bibfileRecord->id, "%d", dstIndex);
  }
  if(dstFileName)
  {
    *dstFileName = AlcStrDup(fileNameBuf);
  }
  if(dstFileType)
  {
    *dstFileType = WlzEffStringToFormat(fileTypeBuf);
  }
  return errNum;
}
  
/*! 
* \ingroup      WlzExtFF
* \brief        Write a warp segmentation parameter record i.e. the image
*		 processing operations and paramters for extracting the
*		gene-expression domain.
* \return       woolz error
* \param    fp	FILE pointer opened for writing
* \param    recordName	bibfile record name
* \param    normFlg	normalisation flag
* \param    histoFlg	histogram equalise flag
* \param    shadeFlg	shade correct flag
* \param    gaussFlg	gauss smoothing flag
* \param    width	gauss smoothing width parameter
* \param    threshLow	ow threshold value
* \param    threshHigh	high threshold value
* \par      Source:
*                WlzExtFFBibUtils.c
*/
WlzErrorNum WlzEffBibWriteWarpInputSegmentationParamsRecord(
  FILE		*fp,
  char		*recordName,
  int		normFlg,
  int		histoFlg,
  int		shadeFlg,
  int		gaussFlg,
  double	width,
  int		threshLow,
  int	 	threshHigh)
{
  char	normStr[8];
  char	histoStr[8];
  char	shadeStr[8];
  char	gaussStr[8];
  char	widthStr[16];
  char	threshLowStr[16];
  char	threshHighStr[16];
  BibFileRecord	*bibfileRecord;
  BibFileError	bibFileErr;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* generate strings for the bibfile entry */
  sprintf(normStr, "%d", normFlg);
  sprintf(histoStr, "%d", histoFlg);
  sprintf(shadeStr, "%d", shadeFlg);
  sprintf(gaussStr, "%d", gaussFlg);
  sprintf(widthStr, "%f", width);
  sprintf(threshLowStr, "%d", threshLow);
  sprintf(threshHighStr, "%d", threshHigh);

 /* create the bibfile record */
  bibfileRecord = 
    BibFileRecordMake(recordName, "0",
		      BibFileFieldMakeVa("Normalise",	normStr,
					 "HistoEqualise",histoStr,
					 "ShadeCorrect",shadeStr,
					 "GaussSmooth",	gaussStr,
					 "GaussWidth",	widthStr,
					 "ThreshLow",	threshLowStr,
					 "ThreshHigh",	threshHighStr,
					 NULL));

  /* write the bibfile and release resources */
  bibFileErr = BibFileRecordWrite(fp, NULL, bibfileRecord);
  BibFileRecordFree(&bibfileRecord);

  /* check the bibfile write error */
  if( bibFileErr != BIBFILE_ER_NONE ){
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }

  return errNum;
}

/*! 
* \ingroup      WlzExtFF
* \brief        Parse the bibfile record for a 
*		MAPaint WarpInputSegmentationParams record.
* \return       woolz error
* \param    bibfileRecord	bibfile record as read by libbibFile
* \param    normFlg	normalisation flag
* \param    histoFlg	histogram equalise flag
* \param    shadeFlg	shade correct flag
* \param    gaussFlg	gauss smoothing flag
* \param    width	gauss smoothing width
* \param    threshLow	low threshold value
* \param    threshHigh	high threshold value
* \par      Source:
*                WlzExtFFBibUtils.c
*/
WlzErrorNum WlzEffBibParseWarpInputSegmentationParamsRecord(
  BibFileRecord		*bibfileRecord,
  int			*normFlg,
  int			*histoFlg,
  int			*shadeFlg,
  int			*gaussFlg,
  double		*width,
  int			*threshLow,
  int	 		*threshHigh)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  char		basisFnTypeStr[32];
  char		meshMthdStr[32];
  int		numParsedFields=0;

  /* check inputs */
  if((bibfileRecord == NULL) ||
     (normFlg == NULL) ||
     (histoFlg == NULL) ||
     (shadeFlg == NULL) ||
     (gaussFlg == NULL) ||
     (width == NULL) ||
     (threshLow == NULL) ||
     (threshHigh == NULL) ){
    errNum = WLZ_ERR_PARAM_NULL;
  }

  /* parse the record */
  if( errNum == WLZ_ERR_NONE ){
    numParsedFields = BibFileFieldParseFmt
      (bibfileRecord->field,
       (void *) normFlg, "%d", "Normalise",
       (void *) histoFlg, "%d", "HistoEqualise",
       (void *) shadeFlg, "%d", "ShadeCorrect",
       (void *) gaussFlg, "%d", "GaussSmooth",
       (void *) width, "%lg", "GaussWidth",
       (void *) threshLow, "%d", "ThreshLow",
       (void *) threshHigh, "%d", "ThreshHigh",
       NULL);
    if( numParsedFields < 4 ){
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }

  return errNum;
}
