

#ifndef WLZIO_H
#define WLZIO_H
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        bibFile.h
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Types and constants for the bibtex based file syntax
*		used for serial section data, ....
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/

#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>
#include <bibFile.h>
#include <Wlz.h>
#include <WlzExtFF.h>


/* From	WlzBibFileIOUtils.c */
extern WlzErrorNum write_Wlz3DSectionViewParams_Record(
  FILE			*fp,
  char			*recordName,
  WlzThreeDViewStruct	*wlzViewStr);
extern WlzErrorNum parse_Wlz3DSectionViewParams_Record(
  BibFileRecord		*bibfileRecord,
  WlzThreeDViewStruct	*wlzViewStr);
extern WlzErrorNum write_WlzWarpTransformParams_Record(
  FILE			*fp,
  char			*recordName,
  WlzBasisFnType	basisFnType,
  WlzMeshGenMethod	meshMthd,
  int			meshMinDst,
  int	 		meshMaxDst);
extern WlzErrorNum parse_WlzWarpTransformParams_Record(
  BibFileRecord		*bibfileRecord,
  WlzBasisFnType	*basisFnType,
  WlzMeshGenMethod	*meshMthd,
  int			*meshMinDst,
  int	 		*meshMaxDst);
extern WlzErrorNum write_WlzTiePointVtxs_Record(
  FILE		*fp,
  char		*recordName,
  int		index,
  WlzDVertex3	dstVtx,
  WlzDVertex3	srcVtx);
extern WlzErrorNum parse_WlzTiePointVtxs_Record(
  BibFileRecord		*bibfileRecord,
  int		*index,
  WlzDVertex3	*dstVtx,
  WlzDVertex3	*srcVtx);
extern WlzErrorNum write_File_Record(
  FILE		*fp,
  char		*recordName,
  char		*fileName,
  WlzEffFormat	fileType);
extern WlzErrorNum parse_File_Record(
  BibFileRecord		*bibfileRecord,
  int		*index,
  char		**dstFileName,
  WlzEffFormat	*dstFileType);
extern WlzErrorNum write_MAPaintWarpInputSegmentationParams_Record(
  FILE		*fp,
  char		*recordName,
  int		normFlg,
  int		histoFlg,
  int		shadeFlg,
  int		gaussFlg,
  double	width,
  int		threshLow,
  int	 	threshHigh);
extern WlzErrorNum parse_MAPaintWarpInputSegmentationParams_Record(
  BibFileRecord		*bibfileRecord,
  int			*normFlg,
  int			*histoFlg,
  int			*shadeFlg,
  int			*gaussFlg,
  double		*width,
  int			*threshLow,
  int	 		*threshHigh);

#ifdef  __cplusplus
}
#endif /* __cplusplus */

#endif /* WLZIO_H */
