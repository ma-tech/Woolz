#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzMakeProperties_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzMakeProperties.c
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
* \brief	Allocation and freeing routines of properties.
* \ingroup	WlzProperty
*/

#ifndef _WIN32
  #include <pwd.h>
#endif
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <Wlz.h>

/*!
* \return	New property list.
* \ingroup	WlzProperty
* \brief	Makes a new property list with a zero link count and a
*		linked list that has no items.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzPropertyList	*WlzMakePropertyList(WlzErrorNum *dstErr)
{
  WlzPropertyList *pList = NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if((pList = (WlzPropertyList *)
  	      AlcCalloc(1, sizeof(WlzPropertyList))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    pList->type = WLZ_PROPERTYLIST;
    if((pList->list = AlcDLPListNew(NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
      AlcFree(pList);
      pList = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(pList);
}

/*!
* \return       New simple property.
* \ingroup      WlzProperty
* \brief        Allocate space for a WlzSimpleProperty which is a
*		structure with size bytes allocated for data.
* \param    	size			Size of space allocated for the
*					property.
* \param    	dstErr			Destination error pointer, may be NULL.
*/
WlzSimpleProperty *WlzMakeSimpleProperty(int size, WlzErrorNum *dstErr)
{
  WlzSimpleProperty *p = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((p = (WlzSimpleProperty *)
       	  AlcCalloc(1,sizeof(WlzSimpleProperty))) == NULL )
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    p->type = WLZ_PROPERTY_SIMPLE;
#ifdef _OPENMP
#pragma omp critical (WlzLinkcount)
#endif
    {
      p->linkcount = 0;
    }
    p->size = size;
    if((p->prop = AlcMalloc(size)) == NULL)
    {
      AlcFree(p);
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      p->freeptr = AlcFreeStackPush(p->freeptr, p->prop, NULL);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(p);
}


/* function:     WlzFreeSimpleProperty    */
/*!
* \ingroup      WlzProperty
* \brief        Free space allocated for a WlzSimpleProperty.
*
* \return       WkzErrorNum
* \param    prop	property to be freed
* \par      Source:
*                WlzMakeProperties.c
*/
WlzErrorNum WlzFreeSimpleProperty(WlzSimpleProperty *prop)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check the object pointer and linkcount */
  if (prop == NULL){
    return( WLZ_ERR_NONE );
  }

  if( WlzUnlink(&(prop->linkcount), &errNum) ){
    if( prop->freeptr != NULL ){
      AlcFreeStackFree( prop->freeptr );
    }
    AlcFree( (void *) prop );
  }

  return errNum;
}

/*!
* \return	New name property.
* \ingroup	WlzProperty
* \brief	Makes a name property from the given name string which
*		is copied.
* \param	name			Given name.
* \param	dstErr			Destination pointer for error number,
*					may be NULL.
*/
WlzNameProperty	*WlzMakeNameProperty(char *name, WlzErrorNum *dstErr)
{
  WlzNameProperty *prop = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(((prop = (WlzNameProperty *)
  	      AlcCalloc(1, sizeof(WlzNameProperty))) == NULL) ||
     ((prop->name = AlcStrDup(name)) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
    AlcFree(prop);
  }
  else
  {
    prop->type = WLZ_PROPERTY_NAME;
    if((prop->freeptr = AlcFreeStackPush(prop->freeptr, prop->name,
    					 NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    prop = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(prop);
}

/*!
* \return	New grey value property.
* \ingroup	WlzProperty
* \brief	Makes a grey value property from the given grey value which
*		is copied.
* \param	name			Optional name string, may be NULL.
* \param	val			Given value.
* \param	dstErr			Destination pointer for error number,
*					may be NULL.
*/
WlzGreyProperty *WlzMakeGreyProperty(char *name, WlzPixelV val,
					WlzErrorNum *dstErr)
{
  WlzGreyProperty *prop = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((prop = (WlzGreyProperty *)
  	     AlcCalloc(1, sizeof(WlzGreyProperty))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    if(name)
    {
      if((prop->name = AlcStrDup(name)) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
        if((prop->freeptr = AlcFreeStackPush(prop->freeptr, prop->name,
					     NULL)) == NULL)
        {
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    prop->type = WLZ_PROPERTY_GREY;
    prop->value = val;
  }
  else if(prop)
  {
    AlcFree(prop->name);
    AlcFree(prop);
    prop = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(prop);
}

/*!
* \return	New text property.
* \ingroup	WlzProperty
* \brief	Makes a text property from the given name and text strings
* 		which are copied.
* \param	name			Optional name string, may be NULL.
* \param	text			Given text string.
* \param	dstErr			Destination pointer for error number,
*					may be NULL.
*/
WlzTextProperty *WlzMakeTextProperty(char *name, char *text,
					WlzErrorNum *dstErr)
{
  WlzTextProperty *prop = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((prop = (WlzTextProperty *)
  	     AlcCalloc(1, sizeof(WlzTextProperty))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    if(name)
    {
      if((prop->name = AlcStrDup(name)) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
        if((prop->freeptr = AlcFreeStackPush(prop->freeptr, prop->name,
					     NULL)) == NULL)
        {
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    prop->type = WLZ_PROPERTY_TEXT;
    if(text)
    {
      if((prop->text = AlcStrDup(text)) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
        if((prop->freeptr = AlcFreeStackPush(prop->freeptr, prop->text,
					     NULL)) == NULL)
        {
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
    }
  }
  if((errNum != WLZ_ERR_NONE) && (prop != NULL))
  {
    AlcFree(prop->name);
    AlcFree(prop);
    prop = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(prop);
}

/*! 
* \ingroup      WlzProperty
* \brief         Make an EMAP property structure.
*
* \return       WlzEMAPProperty
* \param    type	required EMAP property type:
 WLZ_EMAP_PROPERTY_GREY_MODEL, WLZ_EMAP_PROPERTY_ANATOMY_DOMAIN or
 WLZ_EMAP_PROPERTY_OTHER_DOMAIN.
* \param    modelUID	Unique ID for the model (e.g. EMAPM:12345)
* \param    anatomyUID	Unique ID for the anatomy component (e.g. EMAP:12345)
* \param    targetUID	Unique ID to the target model (for transform object).
* \param    targetVersion	Version of the target models
* \param    stage	Model embryonic  stage.
* \param    subStage	Model embryonic  sub-stage.
* \param    modelName	Standard name for the model.
* \param    version	Model version.
* \param    fileName	Filename when model registered.
* \param    comment	Comment string.
* \param    dstErr	Error return.
* \par      Source:
*                WlzMakeProperties.c
*/
WlzEMAPProperty *WlzMakeEMAPProperty(
  WlzEMAPPropertyType	type,
  char			*modelUID,
  char			*anatomyUID,
  char			*targetUID,
  char			*targetVersion,
  char			*stage,
  char			*subStage,
  char			*modelName,
  char			*version,
  char			*fileName,
  char			*comment,
  WlzErrorNum	*dstErr)
{
  WlzEMAPProperty	*rtnProp=NULL;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

#ifndef _WIN32
  struct passwd		*pwdStruct;

  /* check the calling parameters */
  switch( type ){
  case WLZ_EMAP_PROPERTY_GREY_MODEL:
  case WLZ_EMAP_PROPERTY_GREY_OTHER:
  case WLZ_EMAP_PROPERTY_DOMAIN_ANATOMY:
  case WLZ_EMAP_PROPERTY_DOMAIN_OTHER:
  case WLZ_EMAP_PROPERTY_TRANSFORM:
    break;

  default:
    errNum = WLZ_ERR_PARAM_DATA;
    break;
  }
  if( modelUID ){
    if( strlen(modelUID) >= EMAP_PROPERTY_UID_LENGTH ){
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  if( anatomyUID ){
    if( strlen(anatomyUID) >= EMAP_PROPERTY_UID_LENGTH ){
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  if( targetUID ){
    if( strlen(targetUID) >= EMAP_PROPERTY_UID_LENGTH ){
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  if( stage ){
    if( strlen(stage) >= EMAP_PROPERTY_STAGE_LENGTH ){
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  if( subStage ){
    if( strlen(subStage) >= EMAP_PROPERTY_STAGE_LENGTH ){
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }

  /* allocate space and copy in properties */
  if( errNum == WLZ_ERR_NONE ){
    if((rtnProp = (WlzEMAPProperty *)
                  AlcCalloc(1, sizeof(WlzEMAPProperty))) != NULL){
      rtnProp->type = WLZ_PROPERTY_EMAP;
      rtnProp->emapType = type;
      rtnProp->creationTime = time(NULL);
      rtnProp->modificationTime = rtnProp->creationTime;
      pwdStruct = getpwuid(getuid());
      if((strlen(pwdStruct->pw_name) + strlen(pwdStruct->pw_gecos)
	  + 4) < EMAP_PROPERTY_AUTHORNAME_LENGTH){
	sprintf(rtnProp->creationAuthor, "%s (%s)",
		pwdStruct->pw_name,
		pwdStruct->pw_gecos?pwdStruct->pw_gecos:"");
      }
      else {
	strncpy(rtnProp->creationAuthor, pwdStruct->pw_name,
		EMAP_PROPERTY_AUTHORNAME_LENGTH - 1);
      }
      (void )gethostname(&(rtnProp->creationMachineName[0]),
			 EMAP_PROPERTY_MACHINENAME_LENGTH);
      rtnProp->creationMachineName[EMAP_PROPERTY_MACHINENAME_LENGTH-1]
	= '\0';
      strcpy(rtnProp->modificationAuthor, rtnProp->creationAuthor);

      if( modelUID ){
	strncpy(rtnProp->modelUID, modelUID,
		EMAP_PROPERTY_UID_LENGTH - 1);
      }
      if( anatomyUID ){
	strncpy(rtnProp->anatomyUID, anatomyUID,
		EMAP_PROPERTY_UID_LENGTH - 1);
      }
      if( targetUID ){
	strncpy(rtnProp->targetUID, targetUID,
		EMAP_PROPERTY_UID_LENGTH - 1);
      }
      if( stage ){
	strncpy(rtnProp->stage, stage,
		EMAP_PROPERTY_STAGE_LENGTH - 1);
      }
      if( subStage ){
	strncpy(rtnProp->subStage, subStage,
		EMAP_PROPERTY_STAGE_LENGTH - 1);
      }
      if( modelName ){
	strncpy(rtnProp->modelName, modelName,
		EMAP_PROPERTY_MODELNAME_LENGTH - 1);
      }
      if( version ){
	strncpy(rtnProp->version, version,
		EMAP_PROPERTY_VERSION_LENGTH - 1);
      }
      if( fileName ){
	rtnProp->fileName = AlcStrDup(fileName);
	rtnProp->freeptr = AlcFreeStackPush(rtnProp->freeptr,
					    rtnProp->fileName,
					    NULL);
      }
      if( comment ){
	rtnProp->comment = AlcStrDup(comment);
	rtnProp->freeptr = AlcFreeStackPush(rtnProp->freeptr,
					    rtnProp->comment,
					    NULL);
      }
    }
    else {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
#else
  errNum = WLZ_ERR_UNIMPLEMENTED;
#endif
  return rtnProp;
}

/* function:     WlzChangeEMAPProperty    */
/*! 
* \ingroup      WlzProperty
* \brief         Change the values of an EMAP property. Each given
value will be checked against the current value to see if a change is
necessary. If the property is changed the modification time and author
will be updated. NULL for any of modelName, version, fileName, comment
implies no change. To change these please modify directly, updating the
modification time appropriately.
*
* \return       woolz error
* \param    prop	property to be changed
* \param    type	new type
* \param    modelUID	new model UID
* \param    anatomyUID	new anatomy UID
* \param    targetUID	new target UID
* \param    targetVersion	new target version
* \param    stage	new embryonic stage
* \param    subStage	new embryonic sub-stage
* \param    modelName	new model name
* \param    version	new model version
* \param    fileName	new filename
* \param    comment	new comment
* \par      Source:
*                WlzMakeProperties.c
*/
WlzErrorNum WlzChangeEMAPProperty(
  WlzEMAPProperty	*prop,
  WlzEMAPPropertyType	type,
  char			*modelUID,
  char			*anatomyUID,
  char			*targetUID,
  char			*targetVersion,
  char			*stage,
  char			*subStage,
  char			*modelName,
  char			*version,
  char			*fileName,
  char			*comment)
{
  WlzErrorNum		errNum=WLZ_ERR_NONE;

#ifndef _WIN32
   struct passwd		*pwdStruct;
  int			modifiedFlg=0;

  /* check the calling parameters */
  if( prop == NULL ){
    errNum = WLZ_ERR_PROPERTY_NULL;
  }
  else if( prop->type != WLZ_PROPERTY_EMAP ){
    errNum = WLZ_ERR_PROPERTY_TYPE;
  }
  else {
    switch( type ){
    case WLZ_EMAP_PROPERTY_GREY_MODEL:
    case WLZ_EMAP_PROPERTY_GREY_OTHER:
    case WLZ_EMAP_PROPERTY_DOMAIN_ANATOMY:
    case WLZ_EMAP_PROPERTY_DOMAIN_OTHER:
    case WLZ_EMAP_PROPERTY_TRANSFORM:
      break;

    default:
      errNum = WLZ_ERR_PARAM_DATA;
      break;
    }
  }

  /* copy in properties */
  if( errNum == WLZ_ERR_NONE ){

    if( prop->emapType != type ){
      prop->emapType = type;
      modifiedFlg = 1;
    }

    if( modelUID ){
      if( prop->modelUID ){
	if( strcmp(modelUID, prop->modelUID) ){
	  modifiedFlg = 1;
	}
      }
      else {
	modifiedFlg = 1;
      }
      strncpy(prop->modelUID, modelUID,
	      EMAP_PROPERTY_UID_LENGTH - 1);
    }

    if( anatomyUID ){
      if( prop->anatomyUID ){
	if( strcmp(anatomyUID, prop->anatomyUID) ){
	  modifiedFlg = 1;
	}
      }
      else {
	modifiedFlg = 1;
      }
      strncpy(prop->anatomyUID, anatomyUID,
	      EMAP_PROPERTY_UID_LENGTH - 1);
    }

    if( targetUID ){
      if( prop->targetUID ){
	if( strcmp(targetUID, prop->targetUID) ){
	  modifiedFlg = 1;
	}
      }
      else {
	modifiedFlg = 1;
      }
      strncpy(prop->targetUID, targetUID,
	      EMAP_PROPERTY_UID_LENGTH - 1);
    }

    if( targetVersion ){
      if( prop->targetVersion ){
	if( strcmp(targetVersion, prop->targetVersion) ){
	  modifiedFlg = 1;
	}
      }
      else {
	modifiedFlg = 1;
      }
      strncpy(prop->targetVersion, targetVersion,
	      EMAP_PROPERTY_VERSION_LENGTH - 1);
    }

    if( stage ){
      if( prop->stage ){
	if( strcmp(stage, prop->stage) ){
	  modifiedFlg = 1;
	}
      }
      else {
	modifiedFlg = 1;
      }
      strncpy(prop->stage, stage,
	      EMAP_PROPERTY_STAGE_LENGTH - 1);
    }

    if( subStage ){
      if( prop->subStage ){
	if( strcmp(subStage, prop->subStage) ){
	  modifiedFlg = 1;
	}
      }
      else {
	modifiedFlg = 1;
      }
      strncpy(prop->subStage, subStage,
	      EMAP_PROPERTY_STAGE_LENGTH - 1);
    }

    if( modelName ){
      if( prop->modelName ){
	if( strcmp(modelName, prop->modelName) ){
	  modifiedFlg = 1;
	}
      }
      else {
	modifiedFlg = 1;
      }
      strncpy(prop->modelName, modelName,
	      EMAP_PROPERTY_MODELNAME_LENGTH - 1);
    }

    if( version ){
      if( prop->version ){
	if( strcmp(version, prop->version) ){
	  modifiedFlg = 1;
	}
      }
      else {
	modifiedFlg = 1;
      }
      strncpy(prop->version, version,
	      EMAP_PROPERTY_VERSION_LENGTH - 1);
    }

    if( fileName ){
      if( prop->fileName ){
	if( strcmp(fileName, fileName) ){
	  modifiedFlg = 1;
	}
      }
      else {
	modifiedFlg = 1;
      }

      prop->fileName = AlcStrDup(fileName);
      prop->freeptr = AlcFreeStackPush(prop->freeptr,
				       prop->fileName,
				       NULL);
    }
    if( comment ){
      if( prop->comment ){
	if( strcmp(comment, comment) ){
	  modifiedFlg = 1;
	}
      }
      else {
	modifiedFlg = 1;
      }

      prop->comment = AlcStrDup(comment);
      prop->freeptr = AlcFreeStackPush(prop->freeptr,
				       prop->comment,
				       NULL);
    }

    if((errNum == WLZ_ERR_NONE) && modifiedFlg){
      prop->modificationTime = time(NULL);
      pwdStruct = getpwuid(getuid());
      sprintf(prop->modificationAuthor, "%s (%s)",
	      pwdStruct->pw_name,
	      pwdStruct->pw_gecos?pwdStruct->pw_gecos:"");
    }
  }

#else
  errNum = WLZ_ERR_UNIMPLEMENTED;
#endif
  return errNum;
}

/* function:     WlzFreeEMAPProperty    */
/*!
* \ingroup      WlzProperty
* \brief        Free an EMAP property.
*
* \return       woolz error number values: WLZ_ERR_NONE
* \param    prop	property to be freed
* \par      Source:
*                WlzMakeProperties.c
*/
WlzErrorNum WlzFreeEMAPProperty(WlzEMAPProperty *prop)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check the object pointer and linkcount */
  if (prop == NULL){
    return( WLZ_ERR_NONE );
  }

  if( WlzUnlink(&(prop->linkcount), &errNum) ){
    if( prop->freeptr != NULL ){
      AlcFreeStackFree( prop->freeptr );
    }
    AlcFree( (void *) prop );
  }

  return errNum;
}

/* function:     WlzGetProperty    */
/*!
* \return       The required property. If the core value is NULL the
*		requested property is not in the list.
* \ingroup      WlzProperty
* \brief        Get a property of a given type from a property list. It
*		is assumed that there is only one property of a given type
*		in the list.
*
* \param    plist	Property list to search
* \param    type	Property type required
* \param    dstErr	error return values: WLZ_ERR_NONE,
 WLZ_ERR_PROPERTY_NULL.
* \par      Source:
*                WlzMakeProperties.c
*/
WlzProperty WlzGetProperty(
  AlcDLPList		*plist,
  WlzObjectType		type,
  WlzErrorNum 		*dstErr)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzProperty	rtnProp;
  AlcDLPItem	*item;

  rtnProp.core = NULL;

  if( plist ){
    item = plist->head;
    do {
      if( item ){
	if(item->entry &&
	   (((WlzCoreProperty *) (item->entry))->type == type) ){
	  rtnProp.core = (WlzCoreProperty *) (item->entry);
	  break;
	}
	item = item->next;
      }
    } while( item != plist->head );
  }
  else {
    errNum = WLZ_ERR_PROPERTY_NULL;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnProp;
}

/* function:     WlzRemoveProperty    */
/*!
* \ingroup      WlzProperty
* \brief        Remove, without freeing, a property from a property list.
This is typiucally used after a call to WlzGetProperty().
*
* \return       woolz error values: WLZ_ERR_NONE
* \param    plist	Property list to search
* \param    prop	Property to be removed.
* \par      Source:
*                WlzMakeProperties.c
*/
WlzErrorNum WlzRemoveProperty(
  AlcDLPList		*plist,
  WlzProperty		prop)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  AlcDLPItem	*item;

  if( plist == NULL ){
    errNum = WLZ_ERR_PROPERTY_NULL;
  }
  else {
    item = plist->head;
    do {
      if( item ){
	if(item->entry &&
	   ((WlzCoreProperty *) (item->entry) == prop.core) ){
	  (void) AlcDLPItemUnlink(plist, item, 0, NULL);
	  AlcFree(item);
	  break;
	}
	item = item->next;
      }
    } while( item != plist->head );
  }

  return errNum;
}

/* function:     WlzFreeProperty    */
/*!
* \ingroup      WlzProperty
* \brief        Free a woolz property
*
* \return       wooz error values: WLZ_ERR_NONE, WLZ_ERR_MEM_ALLOC.
* \param    prop	property to be freed
* \par      Source:
*                WlzMakeProperties.c
*/
WlzErrorNum WlzFreeProperty(WlzProperty prop)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  AlcErrno	alcErrNum=ALC_ER_NONE;

  if(prop.core == NULL){
    return( WLZ_ERR_NONE );
  }

  if( WlzUnlink(&(prop.core->linkcount), &errNum) ){
    if( prop.core->freeptr != NULL ){
      alcErrNum = AlcFreeStackFree(prop.core->freeptr);
    }
    AlcFree( (void *) prop.core );
    if( alcErrNum != ALC_ER_NONE ){
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }

  return errNum;
}

/*!
* \return       Woolz error code.
* \ingroup      WlzProperty
* \brief        Free a complete property list (including the list
*		structure itself).
* \param    	pList			Property list to be freed.
*/
WlzErrorNum 	WlzFreePropertyList(WlzPropertyList *pList)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(pList)
  {
    if((WlzUnlink(&(pList->linkcount), &errNum)) && pList->list)
    {
      (void )AlcDLPListFree(pList->list);
      AlcFree(pList);
    }
  }
  return(errNum);
}

/*!
* \ingroup      WlzProperty
* \brief        Free a property list entry. This is a procedure that
*		can be passed e.g. to AlcDLPListEntryAppend to enable
*		automatic freeing of properties held on a property list.
* \param    	prop			Property to be freed.
*/
void 		WlzFreePropertyListEntry(void *prop)
{
  WlzCoreProperty *core = (WlzCoreProperty *)prop;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* check the object pointer and linkcount */
  if(prop)
  {
    if(WlzUnlink(&(core->linkcount), &errNum))
    {
      if(core->freeptr != NULL)
      {
	(void )AlcFreeStackFree(core->freeptr);
      }
      AlcFree((void *)prop);
    }
  }
}

/*!
* \return	Pointer to the name property in the given property list
* 		which contains matches the given name string or NULL if
* 		there is no matching name property,
* \ingroup	WlzProperty
* \brief	Finds the name property with the given name string in the
* 		given property list.
* \param	plist			Given property list.
* \param	name			Name string to be matched, or if
* 					NULL the first name property found.
*/
WlzNameProperty			*WlzPropertyListContainsName(
				  WlzPropertyList *plist,
				  char *name)
{
  AlcDLPItem *item;
  WlzNameProperty *prop = NULL;

  if(plist && plist->list && ((item = plist->list->head) != NULL))
  {
    do
    {
      WlzNameProperty *p;

      if(item->entry && ((p = (WlzNameProperty *)(item->entry)) != NULL) &&
         (p->type == WLZ_PROPERTY_NAME) && 
	 ((name == NULL) || (strcmp(p->name, name) == 0)))
      {
        prop = p;
	break;
      }
      item = item->next;
    } while(item && (item != plist->list->head));
  }
  return(prop);
}
