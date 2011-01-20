#ifndef _CIS_MODEL_PACKAGE_H_
#define _CIS_MODEL_PACKAGE_H_
#include <CIS_Model.h>
#include <CIS_3D_Model.h>
#include <CIS_2D_Model.h>

#define MAX_GLOBAL_MODEL_NUM 300

int CIS_Model_Initialize();
int CIS_Model_Release();
CIS_Model* CIS_Model_GetModel(int index);
CIS_3D_Model* CIS_3D_Model_GetModel(int index);
CIS_2D_Model* CIS_2D_Model_GetModel(int index);
int CIS_Model_GetIndex(CIS_Model *model);
int CIS_3D_Model_GetIndex(CIS_3D_Model *model);
int CIS_2D_Model_GetIndex(CIS_2D_Model *model);
int CIS_Model_GetTotalModelNum();
int CIS_3D_Model_GetTotalModelNum();
int CIS_2D_Model_GetTotalModelNum();
int CIS_Model_AddModel(CIS_Model *new_model);
int CIS_3D_Model_AddModel(CIS_3D_Model *new_model);
int CIS_2D_Model_AddModel(CIS_2D_Model *new_model);
int CIS_Model_RemoveModel(CIS_Model *model);
int CIS_3D_Model_RemoveModel(CIS_3D_Model *model);
int CIS_2D_Model_RemoveModel(CIS_2D_Model *model);
int CIS_Model_RemoveAllModel();
int CIS_Model_PrintModel(FILE *fp);
int CIS_Model_ReadModel(FILE *fp);

#endif