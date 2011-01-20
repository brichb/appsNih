#include <qapplication.h>
#include <qstatusbar.h>
#include <qmainwindow.h>
#include <qlineedit.h>
#include <qslider.h>
#include <qradiobutton.h>
#include <qcheckbox.h>
#include <qlistbox.h>
#include <qpushbutton.h>
#include <qdir.h>
#include <qlabel.h>
#include <qcursor.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qdatetime.h>
#include <qspinbox.h>
#include <qtimer.h>
#include <qimage.h>

#include <quaternion.h>

#include "XmlReader.h"

#include <NIH_OpenGL_QT_interface.h>
#include <CIS_Model_Package.h>
#include <Dicom_Reader.h>
#include "NIH_BoneSegmentation_Dlg.h"
#include <CIS_Image_Processing_Algo_3D.h>
#include <CIS_Image_Processing_Algo_2D.h>
#include <NIH_Algo_Region_Growing.h>

#include <NIH_Algo_Curve_Fitting.h>

#include <Interfile_Format.h>
#include <NIH_ITK_NIFTI_Wrapper.h>

#include "NIH_ITK_Utility.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "myQFileDialog.hpp"

#include "VertebralBodySeg.h"

#include <vector>

char global_project_filename[200], global_image_dir[200], global_ALM_batch_file[200];
QString global_ALM_multiple_fn="", global_ALM_multiple_exclude="", global_ALM_multiple_reference="";

int global_batch_start=-1, global_batch_end=-1, global_ALM_batch_id=-1;

extern int Dicom_Get_Dicom_Filename(char *dir,	char *dcm_fn);
extern void CIS_Algo_MarchingCube(int isoValue, int up_th, int low_th,
					  int s, int s_z, CIS_Array_Image3D_short *im_3d, 
					  char *surface_fn, CIS_Triangle_Mesh *surface_3d,
					  int r_x0, int r_x1, int r_y0, int r_y1, int r_z0, int r_z1);
extern void TriangleMesh_LaplacianSmooth(CIS_Triangle_Mesh *mesh, int neighRange, double maxDeform);
extern void
Curve_RigidRegistration_RotationOnly(vec3DynArray &curve1, vec3DynArray &curve2, Vec3 &cen1, Vec3 &cen2, Frame &xform);

extern int NIH_SpineSegmentation_CurvedReformation(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   VertebraStruct2D *vertebraTemplate, 
									   vec3DynArray &spinalCord3D,
									   CIS_Array_Image2D_short *img2D_sag,
									   CIS_Array_Image2D_short *img2D_cor,
									   CIS_Array_Image2D_short *img2D_sag_mask,
									   CIS_Array_Image2D_short *img2D_cor_mask,
									   bool debugMode);
extern int NIH_SpineSegmentation_SpinePartition(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   VertebraStruct2D *vertebraTemplate, 
									   vec3DynArray &spinalCord3D,
									   CIS_Array_Image2D_short *img2D_sag,
									   CIS_Array_Image2D_short *img2D_cor,
									   CIS_Array_Image2D_short *img2D_sag_mask,
									   CIS_Array_Image2D_short *img2D_cor_mask,
									   vec3DynArray &spineCenter, vec3DynArray &spineNormalX, vec3DynArray &spineNormalY, vec3DynArray &spineNormalBack,
									   doubleDynArray &spineWidthLeft, doubleDynArray &spineWidthRight, doubleDynArray &spineWidthUp, doubleDynArray &spineWidthDown, 
									   intDynArray &pedicleValley, intDynArray &backValley,
									   bool debugMode);
extern int NIH_SpineSegmentation_SpinePartition_new(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   VertebraStruct2D *vertebraTemplate, 
									   vec3DynArray &spinalCord3D,
									   CIS_Array_Image2D_short *img2D_sag,
									   CIS_Array_Image2D_short *img2D_cor,
									   CIS_Array_Image2D_short *img2D_sag_mask,
									   CIS_Array_Image2D_short *img2D_cor_mask,
									   vec3DynArray &spineCenter, vec3DynArray &spineNormalX, vec3DynArray &spineNormalY, vec3DynArray &spineNormalBack,
									   doubleDynArray &spineWidthLeft, doubleDynArray &spineWidthRight, doubleDynArray &spineWidthUp, doubleDynArray &spineWidthDown, 
									   intDynArray &pedicleValley, intDynArray &backValley,
									   bool debugMode);
extern int NIH_SpineSegmentation_ProjectColumnModel(CIS_Array_Image3D_short *img3D,
											 int projDir, int window_size, 
										   vec3DynArray &spinalCord3D,
											vec3DynArray &spineCenter, vec3DynArray &spineNormalX, vec3DynArray &spineNormalY, vec3DynArray &spineNormalBAck,
											doubleDynArray &spineWidthLeft, doubleDynArray &spineWidthRight, doubleDynArray &spineWidthUp, doubleDynArray &spineWidthDown, 
											intDynArray &pedicleValley, intDynArray &backValley,
											CIS_2D_Model_Curve *projectedCord,
											CIS_2D_Model_Curve *leftColumn, CIS_2D_Model_Curve *rightColumn,
											CIS_2D_Model_Links *pedicleModel,
											bool debugMode);


extern int NIH_BoneSegmentation_DetectMetasis(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   FeatureStructure* &detections2D, int &numDetections2D,
									   LesionStructure *lesions, int numLesions,
									   PaintingStruct *lesionVoxels,
									   bool flagApplyClassifier, bool flagApply3DClassifier, SvmCommittee *svmCommittee, float svmCutoff);
extern int NIH_BoneSegmentation_LocateDetection2D(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   FeatureStructure* &detections2D, int &numDetections2D,
									   LesionStructure *lesions, int numLesions,
									   PaintingStruct *lesionVoxels,
									   bool flagApplyClassifier, SvmCommittee *svmCommittee, float svmCutoff);
extern int NIH_BoneSegmentation_MergeMetasis(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   DetectionStructure* &detections, int &numDetections,
									   PaintingStruct * detectionVoxels);
extern int NIH_BoneSegmentation_MatchDetections(CIS_Array_Image3D_short *maskImg3D,
									   DetectionStructure *detections, int numDetections,
									   PaintingStruct *detectionVoxels,
									   LesionStructure *lesions, int numLesions,
									   PaintingStruct *lesionVoxels, 
									   int &totalMet, int &totalLytic,
									   int &foundMet, int &foundLytic,
									   int &totalDetection, int &tpDetection, int &fpDetection,
									   int &tpLyticDetection, int &fpLyticDetection,
									   int &totalMet10, int &totalLytic10,
									   int &foundMet10, int &foundLytic10,
									   float lesionSizeCutoff);
extern int NIH_BoneSegmentation_SaveSegmentation(const char *saveSegPath, const char *prefix_fn, CIS_Array_Image3D_short *maskImg3D);
extern int NIH_BoneSegmentation_LoadSegmentation(const char *loadSegPath, const char *prefix_fn, CIS_Array_Image3D_short *&maskImg3D);
extern int NIH_SpineSegmentation_DetectSpinalCord_old(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo, bool debugMode);
extern int NIH_SpineSegmentation_DetectSpinalCord(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   VertebraStruct2D *vertebraTemplate, bool debugMode);
extern int NIH_SpineSegmentation_DetectSpinalCord_new(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   VertebraStruct2D *vertebraTemplate, bool debugMode);
extern int NIH_SpineSegmentation_PreProcess(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo);
extern int NIH_SpineSegmentation_PreProcess_new(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo);
extern int SpineSegmentation_Watershed_usingITK(CIS_Array_Image2D_short *img, CIS_Array_Image2D_short *smoothedData, 
												CIS_Array_Image2D_short *waterLabel, bool debugMode);
extern int	SpineSegmentation_PostProcess_Watershed(CIS_Array_Image3D_short *maskImg3D, 
											CIS_Array_Image2D_short *subImg, 
											CIS_Array_Image2D_short *smoothedData, 
											CIS_Array_Image2D_short *waterLabel, 
											short *binArray2Dedge,
											CIS_Array_Image2D_RGB_uchar *watershedOverlay,
											//CIS_Array_Image2D_RGB_uchar *colorImage,
											int bbx0, int bby0, int curSlice,
											SpineSegmentationMaskStructure &maskStruct,
											SpineSegmentationParameter &segPara,
											SpineSegmentationInfo &segInfo,
											bool debugMode);
extern int	SpineSegmentation_PreProcess_Watershed(CIS_Array_Image3D_short *maskImg3D, 
											CIS_Array_Image2D_short *subImg, 
											int bbx0, int bby0, int curSlice,
											SpineSegmentationMaskStructure &maskStruct,
											SpineSegmentationParameter &segPara,
											SpineSegmentationInfo &segInfo);
extern int NIH_BoneSegmentation_AnalyzeMetasis2D(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   FeatureStructure* &detections2D, int &numDetections2D,
									   LesionStructure *lesions, int numLesions,
									   PaintingStruct *lesionVoxels,
									   bool flagApplyClassifier, bool flagApply3DClassifier, SvmCommittee *svmCommittee, float svmCutoff,
									   int sliceRange1=-1, int sliceRange2=-1);
extern int NIH_BoneSegmentation_AnalyzeMetasis3D(CIS_Array_Image3D_short *img3D,
										  CIS_Array_Image3D_short *maskImg3D,
										  SpineSegmentationMaskStructure &maskStruct,
										  SpineSegmentationParameter &segPara,
										  SpineSegmentationInfo &segInfo,
										  FeatureStructure* &detections3D, int &numDetections3D,
										  DetectionStructure* &detections, int &numDetections,
										  PaintingStruct *detectionVoxels,
										  LesionStructure *lesions, int numLesions,
										  PaintingStruct *lesionVoxels,
										  bool flagApply3DClassifier, SvmCommittee *svmCommittee, float svmCutoff,
										  int sliceRange1, int sliceRange2);
int NIH_SpineSegmentation_BMD_ROI(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   VertebraStruct2D *vertebraTemplate, 
									   doubleDynArray &spineWidthLeft, doubleDynArray &spineWidthRight, doubleDynArray &spineWidthUp, doubleDynArray &spineWidthDown, 
									   intDynArray &pedicleValley, intDynArray &backValley, intDynArray &roiSlices, 
									   bool debugMode, const char *export_fn);
int NIH_SpineSegmentation_BMD_ROI_3D(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   VertebraStruct2D *vertebraTemplate, 
									   vec3DynArray &spineCenter, vec3DynArray &spineNormalX, vec3DynArray &spineNormalY, vec3DynArray &spineNormalBack, 
									   doubleDynArray &spineWidthLeft, doubleDynArray &spineWidthRight, doubleDynArray &spineWidthUp, doubleDynArray &spineWidthDown, 
									   intDynArray &pedicleValley, intDynArray &backValley, intDynArray &roiSlices, 
									   bool debugMode, const char *export_fn);
int NIH_SpineSegmentation_BMD_ROI_L1L2(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   VertebraStruct2D *vertebraTemplate, 
									   doubleDynArray &spineWidthLeft, doubleDynArray &spineWidthRight, doubleDynArray &spineWidthUp, doubleDynArray &spineWidthDown, 
									   intDynArray &pedicleValley, intDynArray &backValley, intDynArray &roiSlices,
									   bool debugMode,
									   const char* export_fn,
									   int l1_slice, CIS_Array_Image3D_short *l1_img,
									   int l2_slice, CIS_Array_Image3D_short *l2_img);
int NIH_SpineSegmentation_BMD_ROI_L1L2_3D(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   VertebraStruct2D *vertebraTemplate, 
									   doubleDynArray &spineWidthLeft, doubleDynArray &spineWidthRight, doubleDynArray &spineWidthUp, doubleDynArray &spineWidthDown, 
									   intDynArray &pedicleValley, intDynArray &backValley, intDynArray &roiSlices,
									   bool debugMode,
									   const char* export_fn,
									   int l1_slice, CIS_Array_Image3D_short *l1_img,
									   int l2_slice, CIS_Array_Image3D_short *l2_img);
int NIH_SpineSegmentation_Rib_Detection(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   VertebraStruct2D *vertebraTemplate, 
									   doubleDynArray &spineWidthLeft, doubleDynArray &spineWidthRight, doubleDynArray &spineWidthUp, doubleDynArray &spineWidthDown, 
									   intDynArray &pedicleValley, intDynArray &backValley, 
									   bool debugMode);
extern int NIH_SpineSegmentation_VertebralDiskDetection(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   VertebraStruct2D *vertebraTemplate, 
									   vec3DynArray &spinalCord3D,
									   vec3DynArray &spineCenter, vec3DynArray &spineNormalX, vec3DynArray &spineNormalY, vec3DynArray &spineNormalBack, 
									   doubleDynArray &spineWidthLeft, doubleDynArray &spineWidthRight, doubleDynArray &spineWidthUp, doubleDynArray &spineWidthDown, 
									   intDynArray &pedicleValley, intDynArray &backValley, 
									   CIS_Array_Image2D_short *img2D_sag,
									   CIS_Array_Image2D_short *img2D_cor,
									   CIS_Array_Image2D_short *img2D_sag_mask,
									   CIS_Array_Image2D_short *img2D_cor_mask,
									   bool debugMode);



NIH_BoneSegmentation_Dlg *theNIH_BoneSegmentation_Dlg=NULL;

void NIH_BoneSegmentation_Dlg_LButtonDownCallBack(QMouseEvent *e, int display_id)
{
	if(theNIH_BoneSegmentation_Dlg==NULL) return;


	IntVec2 pt;
	pt.x = e->x();
	pt.y = e->y();

	if(theNIH_BoneSegmentation_Dlg->flagUseLiveWire && !theNIH_BoneSegmentation_Dlg->flagLiveWireAutoMode)
	{
		theNIH_BoneSegmentation_Dlg->StartLiveWire(pt);
		theNIH_BoneSegmentation_Dlg->liveWireInitModel->vertex_list.SetSize(0);
	}

	if(theNIH_BoneSegmentation_Dlg->flagRecordPoint && !theNIH_BoneSegmentation_Dlg->flagLiveWireAutoMode)
	{
		theNIH_BoneSegmentation_Dlg->recordedPoints.SetSize(0);
		theNIH_BoneSegmentation_Dlg->recordedPoints.Add(pt);
		theNIH_BoneSegmentation_Dlg->liveWireInitModel->vertex_list.SetSize(0);
	}

}

void NIH_BoneSegmentation_Dlg_RButtonDownCallBack(QMouseEvent *e, int display_id)
{
	if(theNIH_BoneSegmentation_Dlg==NULL) return;


	IntVec2 pt;
	pt.x = e->x();
	pt.y = e->y();
	NIH_OpenGL_TurnOnMagnifier(display_id);
	NIH_OpenGL_SetLastMousePosition(pt, display_id);
	NIH_OpenGL_Refresh_Display(display_id);

}

void NIH_BoneSegmentation_Dlg_RButtonUpCallBack(QMouseEvent *e, int display_id)
{
	if(theNIH_BoneSegmentation_Dlg==NULL) return;


	IntVec2 pt;
	pt.x = e->x();
	pt.y = e->y();

	NIH_OpenGL_TurnOffMagnifier(display_id);
	NIH_OpenGL_SetLastMousePosition(pt, display_id);
	NIH_OpenGL_Refresh_Display(display_id);
}

void NIH_BoneSegmentation_Dlg_MouseMoveCallBack(QMouseEvent *e, int display_id)
{
	if(theNIH_BoneSegmentation_Dlg==NULL) return;

	IntVec2 pt;
	pt.x = e->x();
	pt.y = e->y();

	if(e->state()==Qt::RightButton)
	{
		int dx, dy;
		IntVec2 lastMousePt=NIH_OpenGL_GetLastMousePosition(display_id);
		if(lastMousePt.x==-1) return;

		dx = pt.x-lastMousePt.x;
		dy = pt.y-lastMousePt.y;

		if(abs(dx)>5 || abs(dy)>5)
		{
			NIH_OpenGL_SetLastMousePosition(pt, display_id);
			NIH_OpenGL_Refresh_Display(display_id);
		}
	}

	if(theNIH_BoneSegmentation_Dlg->flagUseLiveWire && !theNIH_BoneSegmentation_Dlg->flagLiveWireAutoMode && 
		(e->state() & Qt::LeftButton))
	{
		theNIH_BoneSegmentation_Dlg->AddLiveWireSeed(pt);
		NIH_OpenGL_Refresh_Display(display_id);
	}

	if(theNIH_BoneSegmentation_Dlg->flagRecordPoint && !theNIH_BoneSegmentation_Dlg->flagLiveWireAutoMode &&
		(e->state() & Qt::LeftButton))
	{
		theNIH_BoneSegmentation_Dlg->recordedPoints.Add(pt);
		theNIH_BoneSegmentation_Dlg->liveWireInitModel->vertex_list.Add(pt);
		NIH_OpenGL_Refresh_Display(display_id);
	}

}

void NIH_BoneSegmentation_Dlg_KeyDownCallBack(QKeyEvent *e, int display_id)
{
	if(theNIH_BoneSegmentation_Dlg==NULL) return;
//	int curSlice,totalSlice;
	IntVec2 mousePt;
	
	switch(e->key()) 
	{
	case Qt::Key_0:
		theNIH_BoneSegmentation_Dlg->SwitchALMDisplayMode(0);
		break;
	case Qt::Key_1:
		theNIH_BoneSegmentation_Dlg->SwitchALMDisplayMode(1);
		break;
	case Qt::Key_2:
		theNIH_BoneSegmentation_Dlg->SwitchALMDisplayMode(2);
		break;
	case Qt::Key_3:
		theNIH_BoneSegmentation_Dlg->SwitchALMDisplayMode(3);
		break;
	case Qt::Key_4:
		theNIH_BoneSegmentation_Dlg->SwitchALMDisplayMode(4);
		break;
	case Qt::Key_5:
		theNIH_BoneSegmentation_Dlg->SwitchALMDisplayMode(5);
		break;
	case Qt::Key_6:
		theNIH_BoneSegmentation_Dlg->SwitchALMDisplayMode(6);
		break;
	case Qt::Key_7:
		theNIH_BoneSegmentation_Dlg->SwitchALMDisplayMode(7);
		break;
	case Qt::Key_8:
		theNIH_BoneSegmentation_Dlg->SwitchALMDisplayMode(8);
		break;
	case Qt::Key_9:
		theNIH_BoneSegmentation_Dlg->SwitchALMDisplayMode(9);
		break;
	case Qt::Key_X:
		if(theNIH_BoneSegmentation_Dlg->img3D!=NULL)
		{
			// hide or show texture plane
			if(e->state() & Qt::ControlButton) 
			{
				if(theNIH_BoneSegmentation_Dlg->tex_x->draw_mode==HIDE_MODE) theNIH_BoneSegmentation_Dlg->tex_x->draw_mode=LINE_MODE;
				else theNIH_BoneSegmentation_Dlg->tex_x->draw_mode=HIDE_MODE;
				NIH_OpenGL_Refresh_Display(theNIH_BoneSegmentation_Dlg->display_id2);
			}

			// Centralize at the texture plane
			else if(e->state() & Qt::AltButton) 
			{
				int display_id2=theNIH_BoneSegmentation_Dlg->display_id2;
				Vec3 viewAt = theNIH_BoneSegmentation_Dlg->tex_x->loc;
				double viewAngle = NIH_OpenGL_GetViewAngle(display_id2);
				double focalLength = NIH_OpenGL_GetFocalLength(display_id2);
				IntVec2 viewSize = NIH_OpenGL_GetViewSize(display_id2);
				double dist=theNIH_BoneSegmentation_Dlg->tex_x->sizeY / (tan(viewAngle*3.14159/180)) * 1.5;

				Vec3 viewDir = theNIH_BoneSegmentation_Dlg->tex_x->axis;
				Vec3 viewUp = -theNIH_BoneSegmentation_Dlg->tex_x->dirY;

				Vec3 viewFrom = viewAt-viewDir*dist;
				NIH_OpenGL_Set_ViewFrustum(display_id2, viewFrom, viewAt, viewDir, viewUp, viewSize, focalLength, viewAngle);
				NIH_OpenGL_Refresh_Display(theNIH_BoneSegmentation_Dlg->display_id2);
			}
			else {
				if(theNIH_BoneSegmentation_Dlg->viewMode==2) theNIH_BoneSegmentation_Dlg->ChangeImageSlice(theNIH_BoneSegmentation_Dlg->curSlice+1);
				if(theNIH_BoneSegmentation_Dlg->viewMode==0) theNIH_BoneSegmentation_Dlg->ChangeImageSlice_x(theNIH_BoneSegmentation_Dlg->curSlice_x+1);
				if(theNIH_BoneSegmentation_Dlg->viewMode==1) theNIH_BoneSegmentation_Dlg->ChangeImageSlice_y(theNIH_BoneSegmentation_Dlg->curSlice_y+1);
			}
		}
		break;
	case Qt::Key_Y:
		if(theNIH_BoneSegmentation_Dlg->img3D!=NULL)
		{
			// hide or show texture plane
			if(e->state() & Qt::ControlButton) 
			{
				if(theNIH_BoneSegmentation_Dlg->tex_y->draw_mode==HIDE_MODE) theNIH_BoneSegmentation_Dlg->tex_y->draw_mode=LINE_MODE;
				else theNIH_BoneSegmentation_Dlg->tex_y->draw_mode=HIDE_MODE;
				NIH_OpenGL_Refresh_Display(theNIH_BoneSegmentation_Dlg->display_id2);
			}

			// Centralize at the texture plane
			else if(e->state() & Qt::AltButton) 
			{
				int display_id2=theNIH_BoneSegmentation_Dlg->display_id2;
				Vec3 viewAt = theNIH_BoneSegmentation_Dlg->tex_y->loc;
				double viewAngle = NIH_OpenGL_GetViewAngle(display_id2);
				double focalLength = NIH_OpenGL_GetFocalLength(display_id2);
				IntVec2 viewSize = NIH_OpenGL_GetViewSize(display_id2);
				double dist=theNIH_BoneSegmentation_Dlg->tex_y->sizeY / (tan(viewAngle*3.14159/180)) * 1.5;

				Vec3 viewDir = theNIH_BoneSegmentation_Dlg->tex_y->axis;
				Vec3 viewUp = -theNIH_BoneSegmentation_Dlg->tex_y->dirY;

				Vec3 viewFrom = viewAt-viewDir*dist;
				NIH_OpenGL_Set_ViewFrustum(display_id2, viewFrom, viewAt, viewDir, viewUp, viewSize, focalLength, viewAngle);
				NIH_OpenGL_Refresh_Display(theNIH_BoneSegmentation_Dlg->display_id2);
			}
		}
		break;
	case Qt::Key_Z:
		if(theNIH_BoneSegmentation_Dlg->img3D!=NULL)
		{
			// hide or show texture plane
			if(e->state() & Qt::ControlButton) 
			{
				if(theNIH_BoneSegmentation_Dlg->tex_z->draw_mode==HIDE_MODE) theNIH_BoneSegmentation_Dlg->tex_z->draw_mode=LINE_MODE;
				else theNIH_BoneSegmentation_Dlg->tex_z->draw_mode=HIDE_MODE;
				NIH_OpenGL_Refresh_Display(theNIH_BoneSegmentation_Dlg->display_id2);
			}

			// Centralize at the texture plane
			else if(e->state() & Qt::AltButton) 
			{
				int display_id2=theNIH_BoneSegmentation_Dlg->display_id2;
				Vec3 viewAt = theNIH_BoneSegmentation_Dlg->tex_z->loc;
				double viewAngle = NIH_OpenGL_GetViewAngle(display_id2);
				double focalLength = NIH_OpenGL_GetFocalLength(display_id2);
				IntVec2 viewSize = NIH_OpenGL_GetViewSize(display_id2);
				double dist=theNIH_BoneSegmentation_Dlg->tex_z->sizeY / (tan(viewAngle*3.14159/180)) * 1.5;

				Vec3 viewDir = theNIH_BoneSegmentation_Dlg->tex_z->axis;
				Vec3 viewUp = -theNIH_BoneSegmentation_Dlg->tex_z->dirY;

				Vec3 viewFrom = viewAt-viewDir*dist;
				NIH_OpenGL_Set_ViewFrustum(display_id2, viewFrom, viewAt, viewDir, viewUp, viewSize, focalLength, viewAngle);
				NIH_OpenGL_Refresh_Display(theNIH_BoneSegmentation_Dlg->display_id2);
			}

			else {
				if(theNIH_BoneSegmentation_Dlg->viewMode==2) theNIH_BoneSegmentation_Dlg->ChangeImageSlice(theNIH_BoneSegmentation_Dlg->curSlice-1);
				if(theNIH_BoneSegmentation_Dlg->viewMode==0) theNIH_BoneSegmentation_Dlg->ChangeImageSlice_x(theNIH_BoneSegmentation_Dlg->curSlice_x-1);
				if(theNIH_BoneSegmentation_Dlg->viewMode==1) theNIH_BoneSegmentation_Dlg->ChangeImageSlice_y(theNIH_BoneSegmentation_Dlg->curSlice_y-1);
			}
		}
		break;
	case Qt::Key_V:
		theNIH_BoneSegmentation_Dlg->checkBox_overlay->setChecked(!theNIH_BoneSegmentation_Dlg->flagShowOverlay);
		break;
	case Qt::Key_S:
		mousePt = NIH_OpenGL_GetMousePosition(theNIH_BoneSegmentation_Dlg->display_id);
		theNIH_BoneSegmentation_Dlg->segInfo.bound1.x = mousePt.x;
		theNIH_BoneSegmentation_Dlg->segInfo.bound1.y = mousePt.y;
		theNIH_BoneSegmentation_Dlg->segInfo.bound1.z = theNIH_BoneSegmentation_Dlg->curSlice;;
		NIH_OpenGL_Refresh_Display(theNIH_BoneSegmentation_Dlg->display_id);
		break;
	case Qt::Key_E:
		mousePt = NIH_OpenGL_GetMousePosition(theNIH_BoneSegmentation_Dlg->display_id);
		theNIH_BoneSegmentation_Dlg->segInfo.bound2.x = mousePt.x;
		theNIH_BoneSegmentation_Dlg->segInfo.bound2.y = mousePt.y;
		theNIH_BoneSegmentation_Dlg->segInfo.bound2.z = theNIH_BoneSegmentation_Dlg->curSlice;;
		NIH_OpenGL_Refresh_Display(theNIH_BoneSegmentation_Dlg->display_id);
		break;
	case Qt::Key_F:
		mousePt = NIH_OpenGL_GetMousePosition(theNIH_BoneSegmentation_Dlg->display_id);
		if(theNIH_BoneSegmentation_Dlg->flagUseLiveWire)
		{
			theNIH_BoneSegmentation_Dlg->FreezeLiveWireSeed(mousePt);
		}
		break;
	case Qt::Key_O:
		if(theNIH_BoneSegmentation_Dlg->showOrganMask==true) theNIH_BoneSegmentation_Dlg->showOrganMask=false;
		else theNIH_BoneSegmentation_Dlg->showOrganMask=true;
		if(theNIH_BoneSegmentation_Dlg->viewMode==2) theNIH_BoneSegmentation_Dlg->ChangeImageSlice(theNIH_BoneSegmentation_Dlg->curSlice);
		if(theNIH_BoneSegmentation_Dlg->viewMode==0) theNIH_BoneSegmentation_Dlg->ChangeImageSlice_x(theNIH_BoneSegmentation_Dlg->curSlice_x);
		if(theNIH_BoneSegmentation_Dlg->viewMode==1) theNIH_BoneSegmentation_Dlg->ChangeImageSlice_y(theNIH_BoneSegmentation_Dlg->curSlice_y);

		break;
	}	// switch

}


void NIH_BoneSegmentation_Dlg::CommandLine_clicked()
{
	if(strlen(global_ALM_batch_file)>1)
	{
		ALM_Batch_file(global_ALM_batch_file);
		exit(0);
	}
	if(global_ALM_multiple_fn.length()>2)
	{
		Load_ALM_Multi_Model_File(global_ALM_multiple_fn.ascii());
		QString ALM_fn;
		ALM_fn = "C:\\Experiments\\bone_training\\Atlas\\ALM"+QString::number(totalALMModels)+"_r"+global_ALM_multiple_reference;
		if(global_ALM_multiple_exclude.length()>2) ALM_fn = ALM_fn+"_e"+global_ALM_multiple_exclude;
		ALM_fn = ALM_fn+".xml";
		printf("ALM_fn: %s", ALM_fn.ascii());
		ALM_PDM.SaveALMModel(ALM_fn.ascii());
		exit(0);
	}
	else if(strlen(global_project_filename)>1)
	{
		if(load_project(global_project_filename)==CIS_OK)
		{
			if(global_batch_start>=0 && global_batch_end>=0)
			{
				lineEdit_batch_start->setText(QString::number(global_batch_start));
				lineEdit_batch_end->setText(QString::number(global_batch_end));
				pushButton_batch_project_clicked();
				exit(0);
			}
		}
	}
	else if(strlen(global_image_dir)>1)
	{
		OpenDataSet(global_image_dir, true);
	
		// spine segmentation, partition, labeling, plus rib detection 
		printf("  1.. Preprocessing\n");
		pushButton_pre_process_clicked();
		printf("  2.. Segment spine\n");
		pushButton_spinal_cord_clicked();
		printf("  3.. Spine partition\n");
		pushButton_spine_partition_clicked();
		printf("  4.. Rib detection\n");
		pushButton_Rib_Detection_clicked();

		printf("  5.. Screenshot\n");
		flagShowBBox = true;
		checkBox_bounding_box->setChecked(flagShowBBox);
		char screenshot_root[200];
		sprintf(screenshot_root, "%s\\screenshot\\", global_image_dir);
		// create the direction if not exist
		QDir dcmdir = QString(screenshot_root);
		if(!dcmdir.exists())
		{
			dcmdir.mkdir(screenshot_root);
		}

		char ss_fn[400];
		NIH_OpenGL_Show_Background(display_id, false);
		// screen shot of sagittal view
		ChangeImageSlice_x(curSlice_x);
		sprintf(ss_fn, "%s//sag.jpg", screenshot_root);
		NIH_OpenGL_Save_Display(display_id, ss_fn);
	 	     
		// screen shot of coronal view
		ChangeImageSlice_y(curSlice_y);
		sprintf(ss_fn, "%s//cor.jpg", screenshot_root);
		NIH_OpenGL_Save_Display(display_id, ss_fn);

		cordModel->SetDrawMode(HIDE_MODE);
		smoothedCord->SetDrawMode(HIDE_MODE);

		pushButton_surface_clicked();
		sprintf(ss_fn, "%s//3d.jpg", screenshot_root);
		NIH_OpenGL_Save_Display(display_id2, ss_fn);

		exit(0);
	}
}

NIH_BoneSegmentation_Dlg::NIH_BoneSegmentation_Dlg(int _display_id, int _display_id2, QWidget* parent, const char* name, WFlags fl)
	: NIH_BoneSegmentation_Dlg_Base(parent, name, fl)
{
	Version = (float)1.0;

	display_id = _display_id;
	display_id2 = _display_id2;
	if(display_id==-1) return;

	NIH_OpenGL_Set_Display_Title(display_id, "Bone Segmentation");

	NIH_OpenGL_TurnOnLight(display_id2);
	NIH_OpenGL_Set_Display_Title(display_id2, "Surface");

	img3D = NULL;
	img2D = NULL;
	img2D_x = NULL;
	img2D_y = NULL;
	img2D_sag = NULL;
	img2D_cor = NULL;
	img2D_sag_mask = NULL;
	img2D_cor_mask = NULL;

	maskImg3D = NULL;
	organMask = NULL;
	organGTMask = NULL;
	
	curSlice = -1;
	startSlice = 0;
	endSlice = 0;

	curSlice_x = -1;
	startSlice_x = 0;
	endSlice_x = 0;

	curSlice_y = -1;
	startSlice_y = 0;
	endSlice_y = 0;

	centerx = 256;
	centery = 256;
	centerz = 0;

	strcpy(base_fn, "");
	strcpy(globalFeatureFileName, "");
	strcpy(globalDetectionFileName, "");
	strcpy(globalReportFileName, "");

	debugMode = true;
	reformationMode = false;

	// lesions and detections
	numLesions = numDetections = numDetections2D = numDetections3D = 0;
	lesions = NULL;
	detections = NULL;
	detections2D = NULL;
//	detectionVoxels = NULL;

	lesionSizeCutoff = 0.25;

	// classification
	svmCommittee = new SvmCommittee();
	flagApplyClassifier = false;
	flagApply3DClassifier = false;
	svmCutoff = 0.5;

	svmCutoffMin = 0;
	svmCutoffMax = 1;
	svmCutoffStep = 0;

	for(int i=0; i<MAX_CUTOFF_STEPS; i++)
	{
		//storeMasks[i].mask3DArray;
		storeMasks[i].curCutoff=0;
		storeMasks[i].foundMet=0;
		storeMasks[i].fpDetection=0;
		storeMasks[i].sensitivity=0;
		storeMasks[i].totalDetection=0;
		storeMasks[i].totalMet=0;
		storeMasks[i].tpDetection=0;
		for(int j=0; j<MAX_LESIONS; j++)
		{
			storeMasks[i].lesionStatus[j]=0;
		}
		for(int j=0; j<MAX_DETECTIONS; j++)
		{
			storeMasks[i].matchedLesion[j]=0;
		}
	}

	useMethod = 0;

	viewMode = 2;
	radioButton_view_axial->setChecked(true);


	zoomMode = 1;
	radioButton_zoom1->setChecked(true);

	// graph models
	mapModel1 = new CIS_2D_Model_Bitmap();
	CIS_Model_AddModel(mapModel1);
	NIH_OpenGL_AddModel(display_id, mapModel1);
	
	cmapModel1 = new CIS_2D_Model_Texture();
//	cmapModel1->colorMap = COLORMAP_JET;
	CIS_Model_AddModel(cmapModel1);
	NIH_OpenGL_AddModel(display_id, cmapModel1);

	bboxModel = new CIS_2D_Model_Rectangle();
	CIS_Model_AddModel(bboxModel);
	NIH_OpenGL_AddModel(display_id, bboxModel);

	vertebraBodyModel = new CIS_2D_Model_Rectangle();
	vertebraBodyModel->SetCurrColor(Vec3(0,1,0));
	CIS_Model_AddModel(vertebraBodyModel);
	NIH_OpenGL_AddModel(display_id, vertebraBodyModel);

	spinousProcessModel = new CIS_2D_Model_Rectangle();
	spinousProcessModel->SetCurrColor(Vec3(0,1,0));
	CIS_Model_AddModel(spinousProcessModel);
	NIH_OpenGL_AddModel(display_id, spinousProcessModel);

	paintModel = new CIS_2D_ROI_Polygon();
	paintModel->SetCurrColor(Vec3(1,0,1));
	paintModel->SetLineWidth(1.0);
	CIS_Model_AddModel(paintModel);
	NIH_OpenGL_AddModel(display_id, paintModel);

	seedModel = new CIS_2D_Model_Polygon();
	seedModel->SetCurrColor(Vec3(1,0,0));
	seedModel->SetLineWidth(1.0);
	seedModel->SetDrawMode(VERTEX_MODE);
	seedModel->vertexSize = 2;
	CIS_Model_AddModel(seedModel);
	NIH_OpenGL_AddModel(display_id, seedModel);
	seedModel->vertex_list.SetSize(1);

	metasisSurf = new CIS_3D_Model_MeshSurface();
	metasisSurf->draw_mode = SURFACE_MODE;
	metasisSurf->SetCurrAlpha(1);
	CIS_Model_AddModel(metasisSurf);
	NIH_OpenGL_AddModel(display_id2, metasisSurf);

	cordSurf = new CIS_3D_Model_MeshSurface();
	cordSurf->draw_mode = SURFACE_MODE;
	cordSurf->SetCurrAlpha(0.6);
	CIS_Model_AddModel(cordSurf);
	NIH_OpenGL_AddModel(display_id2, cordSurf);

	spineSurf = new CIS_3D_Model_MeshSurface();
	spineSurf->draw_mode = SURFACE_MODE;
	spineSurf->SetCurrAlpha(0.5);
	CIS_Model_AddModel(spineSurf);
	NIH_OpenGL_AddModel(display_id2, spineSurf);

	ribSurf = new CIS_3D_Model_MeshSurface();
	ribSurf->draw_mode = SURFACE_MODE;
	ribSurf->SetCurrAlpha(1);
	CIS_Model_AddModel(ribSurf);
	NIH_OpenGL_AddModel(display_id2, ribSurf);

	cordModel = new CIS_3D_Model_Scatter();
	CIS_Model_AddModel(cordModel);
	NIH_OpenGL_AddModel(display_id2, cordModel);

//	smoothedCord = new CIS_3D_Model_Curve_Cubic();
	smoothedCord = new CIS_3D_Model_Curve();
	CIS_Model_AddModel(smoothedCord);
	NIH_OpenGL_AddModel(display_id2, smoothedCord);

	projectedCord = new CIS_2D_Model_Curve();
	CIS_Model_AddModel(projectedCord);
	NIH_OpenGL_AddModel(display_id, projectedCord);
	projectedCord->SetDrawMode(HIDE_MODE);

	leftColumn = new CIS_2D_Model_Curve();
	CIS_Model_AddModel(leftColumn);
	NIH_OpenGL_AddModel(display_id, leftColumn);
	leftColumn->SetDrawMode(HIDE_MODE);
	leftColumn->SetCurrColor(Vec3(0,0,1));

	rightColumn = new CIS_2D_Model_Curve();
	CIS_Model_AddModel(rightColumn);
	NIH_OpenGL_AddModel(display_id, rightColumn);
	rightColumn->SetDrawMode(HIDE_MODE);
	rightColumn->SetCurrColor(Vec3(0,0,1));

	pedicleModel = new CIS_2D_Model_Links();
	CIS_Model_AddModel(pedicleModel);
	NIH_OpenGL_AddModel(display_id, pedicleModel);
	pedicleModel->SetDrawMode(HIDE_MODE);
	pedicleModel->SetCurrColor(Vec3(0,1,0));


	projectedCordSag = new CIS_2D_Model_Curve();
	CIS_Model_AddModel(projectedCordSag);
	NIH_OpenGL_AddModel(display_id, projectedCordSag);
	projectedCordSag->SetDrawMode(HIDE_MODE);

	leftColumnSag = new CIS_2D_Model_Curve();
	CIS_Model_AddModel(leftColumnSag);
	NIH_OpenGL_AddModel(display_id, leftColumnSag);
	leftColumnSag->SetDrawMode(HIDE_MODE);
	leftColumnSag->SetCurrColor(Vec3(0,0,1));

	rightColumnSag = new CIS_2D_Model_Curve();
	CIS_Model_AddModel(rightColumnSag);
	NIH_OpenGL_AddModel(display_id, rightColumnSag);
	rightColumnSag->SetDrawMode(HIDE_MODE);
	rightColumnSag->SetCurrColor(Vec3(0,0,1));

	pedicleModelSag = new CIS_2D_Model_Links();
	CIS_Model_AddModel(pedicleModelSag);
	NIH_OpenGL_AddModel(display_id, pedicleModelSag);
	pedicleModelSag->SetDrawMode(HIDE_MODE);
	pedicleModelSag->SetCurrColor(Vec3(0,1,0));


	projectedCordCor = new CIS_2D_Model_Curve();
	CIS_Model_AddModel(projectedCordCor);
	NIH_OpenGL_AddModel(display_id, projectedCordCor);
	projectedCordCor->SetDrawMode(HIDE_MODE);

	leftColumnCor = new CIS_2D_Model_Curve();
	CIS_Model_AddModel(leftColumnCor);
	NIH_OpenGL_AddModel(display_id, leftColumnCor);
	leftColumnCor->SetDrawMode(HIDE_MODE);
	leftColumnCor->SetCurrColor(Vec3(0,0,1));

	rightColumnCor = new CIS_2D_Model_Curve();
	CIS_Model_AddModel(rightColumnCor);
	NIH_OpenGL_AddModel(display_id, rightColumnCor);
	rightColumnCor->SetDrawMode(HIDE_MODE);
	rightColumnCor->SetCurrColor(Vec3(0,0,1));

	pedicleModelCor = new CIS_2D_Model_Links();
	CIS_Model_AddModel(pedicleModelCor);
	NIH_OpenGL_AddModel(display_id, pedicleModelCor);
	pedicleModelCor->SetDrawMode(HIDE_MODE);
	pedicleModelCor->SetCurrColor(Vec3(0,1,0));

	
	vertHeightRuler = new CIS_2D_Model_Links();
	CIS_Model_AddModel(vertHeightRuler);
	NIH_OpenGL_AddModel(display_id, vertHeightRuler);
	vertHeightRuler->SetDrawMode(HIDE_MODE);
	vertHeightRuler->SetCurrColor(Vec3(0,1,1));
	
	vertHeightDir = new CIS_3D_Model_Line();
	CIS_Model_AddModel(vertHeightDir);
	NIH_OpenGL_AddModel(display_id2, vertHeightDir);
	vertHeightDir->SetDrawMode(HIDE_MODE);
	vertHeightDir->SetCurrColor(Vec3(0, 1, 0));


	vertebraTemplateViewMode = 0;
	vertebraTemplate = NULL;
	templateModel = new CIS_2D_Model_Links();
	CIS_Model_AddModel(templateModel);
	NIH_OpenGL_AddModel(display_id, templateModel);
	templateModel->SetDrawMode(HIDE_MODE);
	templateModel->SetCurrColor(Vec3(0,1,0));

	// initial the 3D disk plane
	for(int p=0; p<30; p++)
	{
		diskPlaneModel[p] = new CIS_3D_Model_Plane();
		CIS_Model_AddModel(diskPlaneModel[p]);
		NIH_OpenGL_AddModel(display_id2, diskPlaneModel[p]);
		diskPlaneModel[p]->SetCurrAlpha(0.8);
		diskPlaneModel[p]->SetDrawMode(HIDE_MODE);
	}

	strcpy(patientRoot, "c:\\");
	strcpy(loadSegPath, "c:\\");
	strcpy(saveSegPath, "c:\\");

	numStudyEntries = 0;
	
	flagBatchMode = false;
	segmentationChanged = false;

	theNIH_BoneSegmentation_Dlg = this;

	batchStudyStart = batchStudyEnd=-1;

	// init the dialog status
	UpdateData(false);
	
	slider_contrast_lo->setRange(1,4000);
	slider_contrast_hi->setRange(1,4000);
	slider_contrast_hi->setValue(150);
	slider_contrast_lo->setValue(0);
	NIH_OpenGL_Set_Background_WindowLevel(display_id, 
		150, 75);

	slider_opacity->setRange(0, 10);
	slider_opacity->setValue(4);
	if(cmapModel1!=NULL) cmapModel1->curr_alpha = 0.4;

	textLabel_info->setAlignment(Qt::AlignBottom);

	NIH_OpenGL_Set_LButtonDownCallBack(display_id, NIH_BoneSegmentation_Dlg_LButtonDownCallBack);
	NIH_OpenGL_Set_RButtonDownCallBack(display_id, NIH_BoneSegmentation_Dlg_RButtonDownCallBack);
	NIH_OpenGL_Set_RButtonUpCallBack(display_id, NIH_BoneSegmentation_Dlg_RButtonUpCallBack);
	NIH_OpenGL_Set_MouseMoveCallBack(display_id, NIH_BoneSegmentation_Dlg_MouseMoveCallBack);
	NIH_OpenGL_Set_KeyDownCallBack(display_id, NIH_BoneSegmentation_Dlg_KeyDownCallBack);
	NIH_OpenGL_Set_KeyDownCallBack(display_id2, NIH_BoneSegmentation_Dlg_KeyDownCallBack);

	// set up the callback functions
	NIH_OpenGL_SetLastMousePosition(IntVec2(-1,-1), display_id);

	// set up the callback functions
	NIH_OpenGL_SetLastMousePosition(IntVec2(-1,-1), display_id2);
	NIH_OpenGL_Set_LButtonDownCallBack(display_id2, NIH_OpenGL_ViewChange_LButtonDownCallBack);
	NIH_OpenGL_Set_RButtonDownCallBack(display_id2, NIH_OpenGL_ViewChange_RButtonDownCallBack);
	NIH_OpenGL_Set_MButtonDownCallBack(display_id2, NIH_OpenGL_ViewChange_MButtonDownCallBack);
	NIH_OpenGL_Set_MouseMoveCallBack(display_id2, NIH_OpenGL_ViewChange_MouseMoveCallBack);


	flagUseLiveWire = false;
	flagRecordPoint = false;
	flagLiveWireAutoMode=false;	
	liveWireModel = new CIS_2D_Model_Polygon();
	liveWireModel->SetCurrColor(Vec3(1,0,0));
	liveWireModel->SetLineWidth(2.0);
	liveWireModel->vertexSize = 2;
	CIS_Model_AddModel(liveWireModel);
	NIH_OpenGL_AddModel(display_id, liveWireModel);
	flagImageChanged = true;
	useLiveWire->setChecked(flagUseLiveWire);

	liveWireInitModel = new CIS_2D_Model_Polygon();
	liveWireInitModel->SetCurrColor(Vec3(0,0,1));
	liveWireInitModel->SetLineWidth(2.0);
	liveWireInitModel->draw_mode = LINE_MODE;
	CIS_Model_AddModel(liveWireInitModel);
	NIH_OpenGL_AddModel(display_id, liveWireInitModel);

	groundTruthMask = NULL;

	strcpy(paintingRoot, "");


	// Active Location Model
	for(int a=0; a<MaxALMModels; a++)
	{
		singleALM[a].smoothedCord = new CIS_3D_Model_Curve();
		CIS_Model_AddModel(singleALM[a].smoothedCord);
		NIH_OpenGL_AddModel(display_id2, singleALM[a].smoothedCord);
	
		for(int v=0; v<6; v++)
		{
			singleALM[a].vertebra_3d_model[v] = new CIS_3D_Model_Cylinder();
			CIS_Model_AddModel(singleALM[a].vertebra_3d_model[v]);
			NIH_OpenGL_AddModel(display_id2, singleALM[a].vertebra_3d_model[v]);
			singleALM[a].vertebra_3d_model[v]->SetDrawMode(HIDE_MODE);
		}

		for(int m=0; m<5; m++)
		{
			singleALM[a].organSurf[m] = new CIS_3D_Model_MeshSurface();
			singleALM[a].organSurf[m]->draw_mode = HIDE_MODE;
			singleALM[a].organSurf[m]->SetCurrAlpha(0.8);
			CIS_Model_AddModel(singleALM[a].organSurf[m]);
			NIH_OpenGL_AddModel(display_id2, singleALM[a].organSurf[m]);

			singleALM[a].organLocationModel[m] = new CIS_3D_Model_Needle();
			singleALM[a].organLocationModel[m]->draw_mode = HIDE_MODE;
			CIS_Model_AddModel(singleALM[a].organLocationModel[m]);
			NIH_OpenGL_AddModel(display_id2, singleALM[a].organLocationModel[m]);
		}

		singleALM[a].organSurf[0]->SetCurrColor(Vec3(1,0.5,(float)(a%MaxALMModels)/(MaxALMModels-1)));
		singleALM[a].organLocationModel[0]->SetCurrColor(Vec3(1,0.5,(float)(a%MaxALMModels)/(MaxALMModels-1)));
		singleALM[a].organSurf[1]->SetCurrColor(Vec3(1,1,(float)(a%MaxALMModels)/(MaxALMModels-1)));
		singleALM[a].organLocationModel[1]->SetCurrColor(Vec3(1,1,(float)(a%MaxALMModels)/(MaxALMModels-1)));
		singleALM[a].organSurf[2]->SetCurrColor(Vec3((float)(a%MaxALMModels)/(MaxALMModels-1),1,1));
		singleALM[a].organLocationModel[2]->SetCurrColor(Vec3((float)(a%MaxALMModels)/(MaxALMModels-1),1,1));
		singleALM[a].organSurf[3]->SetCurrColor(Vec3((float)(a%MaxALMModels)/(MaxALMModels-1),1,0.5));
		singleALM[a].organLocationModel[3]->SetCurrColor(Vec3((float)(a%MaxALMModels)/(MaxALMModels-1),1,0.5));
		singleALM[a].organSurf[4]->SetCurrColor(Vec3((float)(a%MaxALMModels)/(MaxALMModels-1),0.5,1));
		singleALM[a].organLocationModel[4]->SetCurrColor(Vec3((float)(a%MaxALMModels)/(MaxALMModels-1),0.5,1));
	}

	// just for two models alignment
/*	for(int m=0; m<5; m++) 
	{
		singleALM[0].organSurf[m]->SetCurrColor(Vec3(1,1,0));
		singleALM[1].organSurf[m]->SetCurrColor(Vec3(0.5,0,1));
	}
*/
	totalALMModels = 0;
	curALMModel = -1;
	GTLocationModel = -1;
	refLocationModel = -1;

	showOrganMask = false;

	tex_x = new CIS_3D_Model_Texture();
	CIS_Model_AddModel(tex_x);
	NIH_OpenGL_AddModel(display_id2, tex_x);

	tex_y = new CIS_3D_Model_Texture();
	CIS_Model_AddModel(tex_y);
	NIH_OpenGL_AddModel(display_id2, tex_y);

	tex_z = new CIS_3D_Model_Texture();
	CIS_Model_AddModel(tex_z);
	NIH_OpenGL_AddModel(display_id2, tex_z);

	tex_x->SetMappingMode(MAPPING_MODE_GL);
	tex_y->SetMappingMode(MAPPING_MODE_GL);
	tex_z->SetMappingMode(MAPPING_MODE_GL);

	tex_x->SetDrawMode(HIDE_MODE);
	tex_y->SetDrawMode(HIDE_MODE);
	tex_z->SetDrawMode(HIDE_MODE);



	volumeModel = NULL;

	flagShowOverlay = true;
	checkBox_overlay->setChecked(flagShowOverlay);

	showWSOverlay = false;
	checkBox_woverlay->setChecked(showWSOverlay);

	flagShowBBox = true;
	checkBox_bounding_box->setChecked(flagShowBBox);

	flagShowPainting = true;
	checkBox_painting->setChecked(flagShowPainting);

	flagGT = false;
	checkBox_flag_GT->setChecked(flagGT);
}


NIH_BoneSegmentation_Dlg::~NIH_BoneSegmentation_Dlg()
{
	theNIH_BoneSegmentation_Dlg = NULL;
	if(img3D) delete img3D;
	if(img2D) delete img2D;
	if(img2D_x) delete img2D_x;
	if(img2D_y) delete img2D_y;
	if(img2D_sag) delete img2D_sag;
	if(img2D_cor) delete img2D_cor;
	if(img2D_sag_mask) delete img2D_sag_mask;
	if(img2D_cor_mask) delete img2D_cor_mask;

	if(maskImg3D) delete maskImg3D;
	if(organMask) delete organMask;
	if(organGTMask) delete organGTMask;

	if(lesions!=NULL) delete lesions;
	if(detections!=NULL) delete detections;
	if(detections2D!=NULL) delete detections2D;
//	if(detectionVoxels!=NULL) delete detectionVoxels;

	if(svmCommittee) delete svmCommittee;
	if(vertebraTemplate) free(vertebraTemplate);

	if(cmapModel1!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, cmapModel1);
		CIS_Model_RemoveModel(cmapModel1);
	}

	if(mapModel1!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, mapModel1);
		CIS_Model_RemoveModel(mapModel1);
	}

	if(bboxModel!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, bboxModel);
		CIS_Model_RemoveModel(bboxModel);
	}

	if(vertebraBodyModel!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, vertebraBodyModel);
		CIS_Model_RemoveModel(vertebraBodyModel);
	}

	if(spinousProcessModel!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, spinousProcessModel);
		CIS_Model_RemoveModel(spinousProcessModel);
	}

	if(paintModel!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, paintModel);
		CIS_Model_RemoveModel(paintModel);
	}
	
	if(spineSurf!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id2, spineSurf);
		CIS_Model_RemoveModel(spineSurf);
	}
	
	if(ribSurf!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id2, ribSurf);
		CIS_Model_RemoveModel(ribSurf);
	}
	
	if(cordSurf!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id2, cordSurf);
		CIS_Model_RemoveModel(cordSurf);
	}
	
	if(cordModel!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id2, cordModel);
		CIS_Model_RemoveModel(cordModel);
	}
	
	if(metasisSurf!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id2, metasisSurf);
		CIS_Model_RemoveModel(metasisSurf);
	}
	
	if(seedModel!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, seedModel);
		CIS_Model_RemoveModel(seedModel);
	}
	
	if(smoothedCord!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id2, smoothedCord);
		CIS_Model_RemoveModel(smoothedCord);
	}
	
	if(projectedCord!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, projectedCord);
		CIS_Model_RemoveModel(projectedCord);
	}
	
	if(leftColumn!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, leftColumn);
		CIS_Model_RemoveModel(leftColumn);
	}
	
	if(rightColumn!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, rightColumn);
		CIS_Model_RemoveModel(rightColumn);
	}
	
	if(pedicleModel!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, pedicleModel);
		CIS_Model_RemoveModel(pedicleModel);
	}
	
	if(projectedCordSag!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, projectedCordSag);
		CIS_Model_RemoveModel(projectedCordSag);
	}
	
	if(leftColumnSag!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, leftColumnSag);
		CIS_Model_RemoveModel(leftColumnSag);
	}
	
	if(rightColumnSag!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, rightColumnSag);
		CIS_Model_RemoveModel(rightColumnSag);
	}
	
	if(pedicleModelSag!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, pedicleModelSag);
		CIS_Model_RemoveModel(pedicleModelSag);
	}
	
	if(projectedCordCor!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, projectedCordCor);
		CIS_Model_RemoveModel(projectedCordCor);
	}
	
	if(leftColumnCor!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, leftColumnCor);
		CIS_Model_RemoveModel(leftColumnCor);
	}
	
	if(rightColumnCor!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, rightColumnCor);
		CIS_Model_RemoveModel(rightColumnCor);
	}
	
	if(pedicleModelCor!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, pedicleModelCor);
		CIS_Model_RemoveModel(pedicleModelCor);
	}
	
	if(templateModel!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, templateModel);
		CIS_Model_RemoveModel(templateModel);
	}
	
	for(int p=0; p<30; p++)
	{
		NIH_OpenGL_RemoveModel(display_id2, diskPlaneModel[p]);
		CIS_Model_RemoveModel(diskPlaneModel[p]);
	}

	for(int a=0; a<MaxALMModels; a++)
	{
		if(singleALM[a].smoothedCord!=NULL)
		{
			NIH_OpenGL_RemoveModel(display_id2, singleALM[a].smoothedCord);
			CIS_Model_RemoveModel(singleALM[a].smoothedCord);
		}
	
		for(int v=0; v<6; v++)
		{
			NIH_OpenGL_RemoveModel(display_id2, singleALM[a].vertebra_3d_model[v]);
			CIS_Model_RemoveModel(singleALM[a].vertebra_3d_model[v]);
		}

		for(int m=0; m<5; m++)
		{
			if(singleALM[a].organSurf[m]!=NULL)
			{
				NIH_OpenGL_RemoveModel(display_id, singleALM[a].organSurf[m]);
				CIS_Model_RemoveModel(singleALM[a].organSurf[m]);
			}
	
			if(singleALM[a].organLocationModel[m]!=NULL)
			{
				NIH_OpenGL_RemoveModel(display_id, singleALM[a].organLocationModel[m]);
				CIS_Model_RemoveModel(singleALM[a].organLocationModel[m]);
			}
		}
	}

	if(tex_x!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, tex_x);
		CIS_Model_RemoveModel(tex_x);
	}
	
	if(tex_y!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, tex_y);
		CIS_Model_RemoveModel(tex_y);
	}
	
	if(tex_z!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, tex_z);
		CIS_Model_RemoveModel(tex_z);
	}
	
	if(volumeModel!=NULL)
	{
		NIH_OpenGL_RemoveModel(display_id, volumeModel);
		CIS_Model_RemoveModel(volumeModel);
	}
	

	// remove the callback funcs
	NIH_OpenGL_Set_LButtonDownCallBack(display_id, NULL);
	NIH_OpenGL_Set_LButtonUpCallBack(display_id, NULL);
	NIH_OpenGL_Set_RButtonDownCallBack(display_id, NULL);
	NIH_OpenGL_Set_RButtonUpCallBack(display_id, NULL);
	NIH_OpenGL_Set_MouseMoveCallBack(display_id, NULL);
	NIH_OpenGL_Set_KeyDownCallBack(display_id, NULL);

}

void NIH_BoneSegmentation_Dlg::UpdateData(bool toEdit)
{
	if(toEdit==true)
	{
///		curSlice = lineEdit_cur_slice->text().toInt();
		batchStudyStart = lineEdit_batch_start->text().toInt();
		batchStudyEnd = lineEdit_batch_end->text().toInt();
		segPara.boneThresh = lineEdit_bone_thresh->text().toInt();
		vertebraTemplateViewMode = lineEdit_template_mode->text().toInt();
	}
	else
	{
		lineEdit_cur_slice->setText(QString::number(curSlice));
		lineEdit_batch_start->setText(QString::number(batchStudyStart));
		lineEdit_batch_end->setText(QString::number(batchStudyEnd));
		lineEdit_bone_thresh->setText(QString::number(segPara.boneThresh));
		lineEdit_template_mode->setText(QString::number(vertebraTemplateViewMode));
	}

}


void NIH_BoneSegmentation_Dlg::pushButton_load_clicked()
{
	showWSOverlay = false;
	UpdateData();

	QString fileExt;

	QString strNewFileName = QFileDialog::getOpenFileName(
		QString::null,
        "All Files (*);; Dicom File (*.dcm)",
         this,
         "Load Studies",
         "Choose a file" );
 

	if(strNewFileName!=QString::null) 
	{
	
		strcpy(base_fn,strNewFileName.ascii());
		strcpy(fn_path1,strNewFileName.ascii());
		// get the directory
		int i;
		for(i=strlen(fn_path1)-1; i>=0; i--)
		{
			if(fn_path1[i]=='\\' || fn_path1[i]=='/') break;
		}
		fn_path1[i+1]=0;

        QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );

		lineEdit_root->setText(fn_path1);

		QTime qtime;

		qtime.start();

		if(strNewFileName.find(".raw")!=-1 || strNewFileName.find(".hdr")!=-1 
			|| strNewFileName.find(".nii")!=-1 || strNewFileName.find(".gz")!=-1)
			OpenDataSet(base_fn, false);
		else OpenDataSet(fn_path1, true);

		char info[100];
		sprintf(info, "Loading time: %ds\n", qtime.elapsed()/1000);
		textLabel_info->setText(info);

		// testing nifti file
//		Write_NIFTI_File("c:\\tmp\\brain2.nii", *img3D);
//		Write_NIFTI_File("c:\\tmp\\brain2.nii.gz", *img3D);

        QApplication::restoreOverrideCursor();
	}
}


void NIH_BoneSegmentation_Dlg::pushButton_load_s_clicked()
{
	if(img3D==NULL) return;

	char fn[200];
 	
	QString strNewFileName = QFileDialog::getOpenFileName(
		QString::null,
        "Image Files (*.img)",
         this,
         "Load Segmentation",
         "Choose a file" );
 
	if(strNewFileName!=QString::null) 
	{
		strcpy(fn,strNewFileName.ascii());

        QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );
		int pos1, pos2;
		pos1 = strNewFileName.findRev('\\');
		if(pos1==-1) pos1 = strNewFileName.findRev('/');
		pos2 = strNewFileName.findRev('.');
		strNewFileName = strNewFileName.mid(pos1+1, pos2-pos1-1);

		NIH_BoneSegmentation_LoadSegmentation(loadSegPath, strNewFileName.ascii(), maskImg3D);
        QApplication::restoreOverrideCursor();
	}
}


void NIH_BoneSegmentation_Dlg::pushButton_save_s_clicked()
{
///	if(maskImg3D==NULL) return;

	char fn[200];
 	
	QString strNewFileName = QFileDialog::getSaveFileName(
		QString::null,
        "Image Files (*.img)",
         this,
         "Save Segmentation",
         "Choose a file" );
 
	if(strNewFileName!=QString::null) 
	{

		strcpy(fn,strNewFileName.ascii());

        QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );
		int pos1, pos2;
		pos1 = strNewFileName.findRev('\\');
		if(pos1==-1) pos1 = strNewFileName.findRev('/');
		pos2 = strNewFileName.findRev('.');
		strNewFileName = strNewFileName.mid(pos1+1, pos2-pos1-1);

		// BINARIZE AND save a copy of the image, for debuging
		int sizex,sizey,sizez;
		sizex = img3D->Num_Cols();
		sizey = img3D->Num_Rows();
		sizez = img3D->Num_Levels();
		short *ar;
		ar = maskImg3D->GetArray();
		int sizexyz=sizex*sizey*sizez;
		for(int k=0; k<sizexyz; k++) 
		{
			if(ar[k]==maskStruct.mask_corticalBone 
				|| ar[k]==maskStruct.mask_spongyBone
				|| ar[k]==maskStruct.mask_rib
//				|| ar[k]==maskStruct.mask_spinalCord
				|| ar[k]==maskStruct.mask_vertebra
///				|| ar[k]==maskStruct.mask_vertebralDisk
				) ar[k]=100;
			else ar[k]=0;
		}

		NIH_BoneSegmentation_SaveSegmentation(saveSegPath, strNewFileName.ascii(), maskImg3D);
        QApplication::restoreOverrideCursor();
	}
}

void NIH_BoneSegmentation_Dlg::pushButton_save_seg_clicked()
{
	if(curStudy>=0 || curStudy<numStudyEntries)
		NIH_BoneSegmentation_SaveSegmentation(saveSegPath, studyEntries[curStudy].localImagePath, maskImg3D);
}

void NIH_BoneSegmentation_Dlg::pushButton_save_feature_clicked()
{
	if(numDetections2D==0) return;

	QString strNewFileName = QFileDialog::getSaveFileName(
		QString::null,
        "csv Files (*.csv)",
         this,
         "Save Detection Features",
         "Choose a file" );
 
	if(strNewFileName!=QString::null) 
	{
		WriteOutFeatures(strNewFileName.ascii());
	}
}

int NIH_BoneSegmentation_Dlg::WriteOutFeatures( const char *feature_fn)
{
	if(numDetections2D==0) return 1;
	static bool firstTime=true;
	FILE *fp;

	if(firstTime)
	{
		fp = fopen(feature_fn, "w");
		firstTime = false;
	}
	else
	{
		fp = fopen(feature_fn, "a");
	}


	if(fp!=NULL)
	{
		/// print a header
		fprintf(fp, "CTSeriesName,");
		fprintf(fp, "centerx, centery, centerz,");
		fprintf(fp, "distToBoundary,");
		fprintf(fp, "relCoordx, relCoordy, relCoordz,");
		fprintf(fp, "area, volume, perimeter,");
		fprintf(fp, "primaryAxisLength, secondaryAxisLength,");
		fprintf(fp, "outerBorderRatio, aspectRatio,");
		fprintf(fp, "spherecity, compactness, roundness,");
		fprintf(fp, "shapeComplexity_f1, shapeComplexity_f2, shapeComplexity_f21,");
		fprintf(fp, "meanIntensity, stdevIntensity, skewnessIntensity, kurtosisIntensity,");
		fprintf(fp, "interiorIntensity, borderIntensity,");
		fprintf(fp, "outsideIntensity, outsideIntensityDev,");
		fprintf(fp, "innerOuterContrast,");
		fprintf(fp, "borderThickness,");
		fprintf(fp, "matchedLesion, matchedLesionType, matchedLesionSize, matchedOverlap\n");

		for(int i=0; i<numDetections2D; i++)
		{
			fprintf(fp, "%s,", studyEntries[curStudy].patientName.ascii());
			fprintf(fp, "%d, %d, %d,", detections2D[i].centerx,detections2D[i].centery,detections2D[i].centerz);
			fprintf(fp, "%f, ", detections2D[i].distToBoundary);
			fprintf(fp, "%f, %f, %f,", detections2D[i].relCoordx,detections2D[i].relCoordy,detections2D[i].relCoordz);
			fprintf(fp, "%f, %f, %f,", detections2D[i].area,detections2D[i].volume,detections2D[i].perimeter);
			fprintf(fp, "%f, %f,", detections2D[i].primaryAxisLength, detections2D[i].secondaryAxisLength);
			fprintf(fp, "%f, %f,", detections2D[i].outerBorderRatio, detections2D[i].aspectRatio);
			fprintf(fp, "%f, %f, %f,", detections2D[i].spherecity,detections2D[i].compactness,detections2D[i].roundness);
			fprintf(fp, "%f, %f, %f,", detections2D[i].shapeComplexity_f1,detections2D[i].shapeComplexity_f2,detections2D[i].shapeComplexity_f21);
			fprintf(fp, "%f, %f, %f, %f,", detections2D[i].meanIntensity,detections2D[i].stdevIntensity,detections2D[i].skewnessIntensity,detections2D[i].kurtosisIntensity);
			fprintf(fp, "%f, %f,", detections2D[i].interiorIntensity, detections2D[i].borderIntensity);
			fprintf(fp, "%f, %f,", detections2D[i].outsideIntensity, detections2D[i].outsideIntensityDev);
			fprintf(fp, "%f, ", detections2D[i].innerOuterContrast);
			fprintf(fp, "%f, ", detections2D[i].borderThickness);
			fprintf(fp, "%d, %d, %f, %f\n", detections2D[i].matchedLesion, detections2D[i].matchedLesionType, detections2D[i].matchedLesionSize, detections2D[i].matchedOverlap);
		}

		fclose(fp);
	}
	
	return 0;
}

int NIH_BoneSegmentation_Dlg::WriteOutFeatures3D( const char *feature_fn)
{
	if(numDetections3D==0) return 1;
	static bool firstTime=true;
	FILE *fp;

	if(firstTime)
	{
		fp = fopen(feature_fn, "w");
		firstTime = false;
	}
	else
	{
		fp = fopen(feature_fn, "a");
	}

	if(fp!=NULL)
	{
		/// print a header
		fprintf(fp, "CTSeriesName,");
		fprintf(fp, "centerx, centery, centerz,");
		fprintf(fp, "distToBoundary,");
		fprintf(fp, "relCoordx, relCoordy, relCoordz,");
		fprintf(fp, "surfaceArea, volume,");
		fprintf(fp, "primaryAxisLength, secondaryAxisLength,");
		fprintf(fp, "outerBorderRatio, aspectRatio,");
		fprintf(fp, "spherecity,");
		fprintf(fp, "shapeComplexity_f1, shapeComplexity_f2, shapeComplexity_f21,");
		fprintf(fp, "meanIntensity, stdevIntensity, skewnessIntensity, kurtosisIntensity,");
		fprintf(fp, "interiorIntensity, borderIntensity,");
		fprintf(fp, "outsideIntensity, outsideIntensityDev,");
		fprintf(fp, "innerOuterContrast,");
		fprintf(fp, "borderThickness,");
		fprintf(fp, "matchedLesion, matchedLesionType, matchedLesionSize, matchedOverlap\n");

		for(int i=0; i<numDetections3D; i++)
		{
			fprintf(fp, "%s,", studyEntries[curStudy].patientName.ascii());
			fprintf(fp, "%d, %d, %d,", detections3D[i].centerx,detections3D[i].centery,detections3D[i].centerz);
			fprintf(fp, "%f, ", detections3D[i].distToBoundary);
			fprintf(fp, "%f, %f, %f,", detections3D[i].relCoordx,detections3D[i].relCoordy,detections3D[i].relCoordz);
			fprintf(fp, "%f, %f,", detections3D[i].surfaceArea,detections3D[i].volume);
			fprintf(fp, "%f, %f,", detections3D[i].primaryAxisLength, detections3D[i].secondaryAxisLength);
			fprintf(fp, "%f, %f,", detections3D[i].outerBorderRatio, detections3D[i].aspectRatio);
			fprintf(fp, "%f,", detections3D[i].spherecity);
			fprintf(fp, "%f, %f, %f,", detections3D[i].shapeComplexity_f1,detections3D[i].shapeComplexity_f2,detections3D[i].shapeComplexity_f21);
			fprintf(fp, "%f, %f, %f, %f,", detections3D[i].meanIntensity,detections3D[i].stdevIntensity,detections3D[i].skewnessIntensity,detections3D[i].kurtosisIntensity);
			fprintf(fp, "%f, %f,", detections3D[i].interiorIntensity, detections3D[i].borderIntensity);
			fprintf(fp, "%f, %f,", detections3D[i].outsideIntensity, detections3D[i].outsideIntensityDev);
			fprintf(fp, "%f, ", detections3D[i].innerOuterContrast);
			fprintf(fp, "%f, ", detections3D[i].borderThickness);
			fprintf(fp, "%d, %d, %f, %f\n", detections3D[i].matchedLesion, detections3D[i].matchedLesionType, detections3D[i].matchedLesionSize, detections3D[i].matchedOverlap);
		}

		fclose(fp);
	}
	
	return 0;
}
int NIH_BoneSegmentation_Dlg::WriteOutDetections( const char *detection_fn)
{
	if(numDetections==0) return 1;
	static bool firstTime=true;
	FILE *fp;

	if(firstTime)
	{
		fp = fopen(detection_fn, "w");
		firstTime = false;
	}
	else
	{
		fp = fopen(detection_fn, "a");
	}


	if(fp!=NULL)
	{
		/// print a header
		fprintf(fp, "CTSeriesName,");
		fprintf(fp, "centerx, centery, centerz,");
		fprintf(fp, "lesionVol, lesionSize,");
		fprintf(fp, "matchedLesion, matchedLesionType\n");

		for(int i=0; i<numDetections; i++)
		{
			fprintf(fp, "%s,", studyEntries[curStudy].patientName.ascii());
			fprintf(fp, "%d, %d, %d,", detections[i].detectionLoc.x,detections[i].detectionLoc.y,detections[i].detectionLoc.z);
			fprintf(fp, "%f, %f,", detections[i].detectionVolume, detections[i].detectionSize);
			fprintf(fp, "%d, %d\n", detections[i].matchedLesion, detections[i].detectionType);
		}

		fclose(fp);
	}
	
	return 0;
}



void NIH_BoneSegmentation_Dlg::slider_slice_valueChanged( int )
{
	if(img3D==NULL) return;
	int newCurSlice = slider_slice->value();
	lineEdit_cur_slice->setText(QString::number(newCurSlice));
	
	if(viewMode==2) ChangeImageSlice(newCurSlice);
	if(viewMode==0) ChangeImageSlice_x(newCurSlice);
	if(viewMode==1) ChangeImageSlice_y(newCurSlice);

}

void NIH_BoneSegmentation_Dlg::slider_slice_sliderMoved( int )
{
	if(img3D==NULL) return;
	int newCurSlice = slider_slice->value();
	lineEdit_cur_slice->setText(QString::number(newCurSlice));
	
	if(viewMode==2) ChangeImageSlice(newCurSlice);
	if(viewMode==0) ChangeImageSlice_x(newCurSlice);
	if(viewMode==1) ChangeImageSlice_y(newCurSlice);

}


void NIH_BoneSegmentation_Dlg::slider_opacity_valueChanged( int )
{
	if(img3D==NULL) return;
	int newOpacity = slider_opacity->value();
	if(cmapModel1!=NULL) cmapModel1->curr_alpha = (double)newOpacity/10;
	if(volumeModel)
	{
		volumeModel->SetCurrAlpha((double)newOpacity/10);
	}

	if(viewMode==2) ChangeImageSlice(curSlice);
	if(viewMode==0) ChangeImageSlice_x(curSlice_x);
	if(viewMode==1) ChangeImageSlice_y(curSlice_y);

}

void NIH_BoneSegmentation_Dlg::slider_contrast_lo_valueChanged( int )
{
	if(img2D==NULL) return;
	int newHiLimit = slider_contrast_hi->value();
	int newLowLimit = slider_contrast_lo->value();

	NIH_OpenGL_Set_Background_WindowLevel(display_id, 
		newHiLimit-newLowLimit, (newHiLimit+newLowLimit)/2);

	if(tex_x->draw_mode!=HIDE_MODE) tex_x->SetWindowLevel(newHiLimit-newLowLimit, (newHiLimit+newLowLimit)/2);
	if(tex_y->draw_mode!=HIDE_MODE) tex_y->SetWindowLevel(newHiLimit-newLowLimit, (newHiLimit+newLowLimit)/2);
	if(tex_z->draw_mode!=HIDE_MODE) tex_z->SetWindowLevel(newHiLimit-newLowLimit, (newHiLimit+newLowLimit)/2);

	NIH_OpenGL_Refresh_Display(display_id);
	NIH_OpenGL_Refresh_Display(display_id2);
}


void NIH_BoneSegmentation_Dlg::slider_contrast_hi_valueChanged( int )
{
	if(img2D==NULL) return;
	int newHiLimit = slider_contrast_hi->value();
	int newLowLimit = slider_contrast_lo->value();

	NIH_OpenGL_Set_Background_WindowLevel(display_id, 
		newHiLimit-newLowLimit, (newHiLimit+newLowLimit)/2);

	if(tex_x->draw_mode!=HIDE_MODE) tex_x->SetWindowLevel(newHiLimit-newLowLimit, (newHiLimit+newLowLimit)/2);
	if(tex_y->draw_mode!=HIDE_MODE) tex_y->SetWindowLevel(newHiLimit-newLowLimit, (newHiLimit+newLowLimit)/2);
	if(tex_z->draw_mode!=HIDE_MODE) tex_z->SetWindowLevel(newHiLimit-newLowLimit, (newHiLimit+newLowLimit)/2);

	NIH_OpenGL_Refresh_Display(display_id);
	NIH_OpenGL_Refresh_Display(display_id2);
}

void NIH_BoneSegmentation_Dlg::spinBox_lesion_valueChanged( int )
{
	if(img3D==NULL) return;

	int curLesion = spinBox_lesion->value();

	if(curLesion<0 || curLesion>=totalPaintNum) return;
	char info[200];

	sprintf(info,"CurLesion:%d; Key:%d, Type:%d\nLoc: %d %d %d\nSize:%f; Slices:%d\nStatus:%d", 
		curLesion, lesions[curLesion].lesionKey, lesions[curLesion].lesionType, 
		lesions[curLesion].lesionLoc.x, lesions[curLesion].lesionLoc.y, lesions[curLesion].lesionLoc.z,
		lesions[curLesion].lesionSize, lesions[curLesion].lesionSlices,
		lesions[curLesion].lesionStatus);

	textLabel_info->setText(info);

	centerx = lesions[curLesion].lesionLoc.x;
	centery = lesions[curLesion].lesionLoc.y;
	centerz = lesions[curLesion].lesionLoc.z;

	if(viewMode==2)
	{
		seedModel->vertex_list[0] = Vec2(lesions[curLesion].lesionLoc.x, lesions[curLesion].lesionLoc.y);
	
		ChangeImageSlice(lesions[curLesion].lesionLoc.z);
	}
	else if(viewMode==0)
	{
		seedModel->vertex_list[0] = Vec2(lesions[curLesion].lesionLoc.y*mapModel1->size.x/img3D->Num_Rows()+mapModel1->loc.x, 
			lesions[curLesion].lesionLoc.z*mapModel1->size.y/img3D->Num_Levels()+mapModel1->loc.y);
	
		ChangeImageSlice_x(lesions[curLesion].lesionLoc.x);
	}
	else if(viewMode==1)
	{
		seedModel->vertex_list[0] = Vec2(lesions[curLesion].lesionLoc.x*mapModel1->size.x/img3D->Num_Cols()+mapModel1->loc.x, 
			lesions[curLesion].lesionLoc.z*mapModel1->size.y/img3D->Num_Levels()+mapModel1->loc.y);
	
		ChangeImageSlice_y(lesions[curLesion].lesionLoc.y);
	}
}

void NIH_BoneSegmentation_Dlg::spinBox_detection_valueChanged( int )
{
	if(img3D==NULL) return;

	int curDetection = spinBox_detection->value();

	if(curDetection<0 || curDetection>=numDetections) return;
	char info[200];

	sprintf(info,"CurDetection:%d\nLoc: %d %d %d\nVol:%.2f, Size:%.2f; z-interval: [%d, %d]\n Matched Lesion:%d", 
		curDetection, 
		detections[curDetection].detectionLoc.x, detections[curDetection].detectionLoc.y, detections[curDetection].detectionLoc.z,
		detections[curDetection].detectionVolume, detections[curDetection].detectionSize,
		detections[curDetection].minz,
		detections[curDetection].maxz,
		detections[curDetection].matchedLesion);

	textLabel_info->setText(info);

	centerx = detections[curDetection].detectionLoc.x;
	centery = detections[curDetection].detectionLoc.y;
	centerz = detections[curDetection].detectionLoc.z;


	if(viewMode==2)
	{
		seedModel->vertex_list[0] = Vec2(detections[curDetection].detectionLoc.x, detections[curDetection].detectionLoc.y);
	
		ChangeImageSlice(detections[curDetection].detectionLoc.z);
	}
	else if(viewMode==0)
	{
		seedModel->vertex_list[0] = Vec2(detections[curDetection].detectionLoc.y*mapModel1->size.x/img3D->Num_Rows()+mapModel1->loc.x, 
			detections[curDetection].detectionLoc.z*mapModel1->size.y/img3D->Num_Levels()+mapModel1->loc.y);
	
		ChangeImageSlice_x(detections[curDetection].detectionLoc.x);
	}
	else if(viewMode==1)
	{
		seedModel->vertex_list[0] = Vec2(detections[curDetection].detectionLoc.x*mapModel1->size.x/img3D->Num_Cols()+mapModel1->loc.x, 
			detections[curDetection].detectionLoc.z*mapModel1->size.y/img3D->Num_Levels()+mapModel1->loc.y);
	
		ChangeImageSlice_y(detections[curDetection].detectionLoc.y);
	}
}

void NIH_BoneSegmentation_Dlg::spinBox_detection_2D_valueChanged( int )
{
	if(img3D==NULL) return;

	int curDetection = spinBox_detection_2D->value();

	if(curDetection<0 || curDetection>=numDetections2D) return;
	char info[200];

	sprintf(info,"CurDetection:%d\nLoc: %d %d %d\nArea:%.2f, Size:%.2f; Matched Lesion:%d\nSvm:%f", 
		curDetection, 
		detections2D[curDetection].centerx, detections2D[curDetection].centery, detections2D[curDetection].centerz,
		detections2D[curDetection].area, detections2D[curDetection].primaryAxisLength,
		detections2D[curDetection].matchedLesion,
		detections2D[curDetection].svmScore);

	textLabel_info->setText(info);

	centerx = detections2D[curDetection].centerx;
	centery = detections2D[curDetection].centery;
	centerz = detections2D[curDetection].centerz;


	if(viewMode==2)
	{
		seedModel->vertex_list[0] = Vec2(detections2D[curDetection].centerx, detections2D[curDetection].centery);
	
		ChangeImageSlice(detections2D[curDetection].centerz);
	}
	else if(viewMode==0)
	{
		seedModel->vertex_list[0] = Vec2(detections2D[curDetection].centery*mapModel1->size.x/img3D->Num_Rows()+mapModel1->loc.x, 
			detections2D[curDetection].centerz*mapModel1->size.y/img3D->Num_Levels()+mapModel1->loc.y);
	
		ChangeImageSlice_x(detections2D[curDetection].centerx);
	}
	else if(viewMode==1)
	{
		seedModel->vertex_list[0] = Vec2(detections2D[curDetection].centerx*mapModel1->size.x/img3D->Num_Cols()+mapModel1->loc.x, 
			detections2D[curDetection].centerz*mapModel1->size.y/img3D->Num_Levels()+mapModel1->loc.y);
	
		ChangeImageSlice_y(detections2D[curDetection].centery);
	}

}

void NIH_BoneSegmentation_Dlg::spinBox_cutoff_valueChanged( int )
{
	if(img3D==NULL || maskImg3D==NULL) return;

	int sizex = maskImg3D->Num_Cols();
	int sizey = maskImg3D->Num_Rows();
	int sizez = maskImg3D->Num_Levels();
	int sizexy=sizex*sizey;

	short *maskA=maskImg3D->GetArray();
	int curIndex = spinBox_cutoff->value();
	int totalMet, foundMet, totalDetection, tpDetection, fpDetection, i;
	float sensitivity;

	if((!flagApplyClassifier && !flagApply3DClassifier) || svmCutoffStep<=0 || curIndex<0 || curIndex>=(svmCutoffMax-svmCutoffMin)/svmCutoffStep)
		return;

	std::cout<<"\n"<<storeMasks[0].mask3DArray[138008];

	for(i=0; i<sizexy*sizez; i++)
	{
		maskA[i]=storeMasks[curIndex].mask3DArray[i];
	}
	totalMet = storeMasks[curIndex].totalMet;
	foundMet = storeMasks[curIndex].foundMet;
	totalDetection = storeMasks[curIndex].totalDetection;
	tpDetection = storeMasks[curIndex].tpDetection;
	fpDetection = storeMasks[curIndex].fpDetection;
	sensitivity = storeMasks[curIndex].sensitivity;

	for(i=0; i<numLesions; i++)
	{
		lesions[i].lesionStatus = storeMasks[curIndex].lesionStatus[i];
	}
	for(i=0; i<numDetections; i++)
	{
		detections[i].matchedLesion = storeMasks[curIndex].matchedLesion[i];
	}

	char info[500];
	sprintf(info, "Total Mets: %d; Found: %d;\n Detection: %d; TP: %d; FP:\n %d; Cutoff: %f\n; Sensitivity: %f\n",
			totalMet, foundMet,
			totalDetection, tpDetection, fpDetection, curIndex*svmCutoffStep, sensitivity);

		textLabel_info->setText(info);

	ChangeImageSlice(slider_slice->value());
}

void NIH_BoneSegmentation_Dlg::checkBox_overlay_stateChanged( int )
{
	flagShowOverlay = checkBox_overlay->isChecked();

	if(viewMode==2) ChangeImageSlice(curSlice);
	if(viewMode==0) ChangeImageSlice_x(curSlice_x);
	if(viewMode==1) ChangeImageSlice_y(curSlice_y);

}

void NIH_BoneSegmentation_Dlg::checkBox_woverlay_stateChanged( int )
{
	showWSOverlay = checkBox_woverlay->isChecked();

	if(viewMode==2) ChangeImageSlice(curSlice);
	if(viewMode==0) ChangeImageSlice_x(curSlice_x);
	if(viewMode==1) ChangeImageSlice_y(curSlice_y);

}

void NIH_BoneSegmentation_Dlg::checkBox_bounding_box_stateChanged( int )
{
	flagShowBBox = checkBox_bounding_box->isChecked();
	
	if(flagShowBBox)
	{
		for(int p=0; p<pedicleValley.GetSize() && p<30; p++)
		{
			diskPlaneModel[p]->SetDrawMode(SURFACE_MODE);
		}
	} else
	{
		for(int p=0; p<30; p++)
		{
			diskPlaneModel[p]->SetDrawMode(HIDE_MODE);
		}
	}
	if(viewMode==2) ChangeImageSlice(curSlice);
	if(viewMode==0) ChangeImageSlice_x(curSlice_x);
	if(viewMode==1) ChangeImageSlice_y(curSlice_y);

	NIH_OpenGL_Refresh_Display(display_id);
	NIH_OpenGL_Refresh_Display(display_id2);

}

void NIH_BoneSegmentation_Dlg::checkBox_painting_stateChanged( int )
{
	flagShowPainting = checkBox_painting->isChecked();

	if(viewMode==2) ChangeImageSlice(curSlice);
	if(viewMode==0) ChangeImageSlice_x(curSlice_x);
	if(viewMode==1) ChangeImageSlice_y(curSlice_y);

}

void NIH_BoneSegmentation_Dlg::checkBox_flag_GT_stateChanged( int )
{
	flagGT = checkBox_flag_GT->isChecked();
}

void NIH_BoneSegmentation_Dlg::pushButton_spine_partition_clicked()
{
	if(maskImg3D==NULL || img3D==NULL || cordModel==NULL) return;

	reformationMode = true;
	QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );

	NIH_SpineSegmentation_CurvedReformation(img3D,
									   maskImg3D,
									   maskStruct,
									   segPara,
									   segInfo,
									   vertebraTemplate, 
									   spinalCord3D,
									   img2D_sag, img2D_cor,
									   img2D_sag_mask, img2D_cor_mask,
									   debugMode);	
	
	smoothedCord->UpdateControlList(spinalCord3D);
	
	NIH_SpineSegmentation_SpinePartition_new(img3D,
									   maskImg3D,
									   maskStruct,
									   segPara,
									   segInfo,
									   vertebraTemplate, 
									   spinalCord3D,
									   img2D_sag, img2D_cor,
									   img2D_sag_mask, img2D_cor_mask,
									   spineCenter, spineNormalX, spineNormalY, spineNormalBack,
									   spineWidthLeft, spineWidthRight, spineWidthUp, spineWidthDown, 
									   pedicleValley, backValley,
									   debugMode);

	// update the 2D superimposed models
	NIH_SpineSegmentation_ProjectColumnModel(img3D,
											 0, 512, 
											spinalCord3D,
											spineCenter, spineNormalX, spineNormalY, spineNormalBack,
											spineWidthLeft, spineWidthRight, spineWidthUp, spineWidthDown, 
											pedicleValley, backValley,
											projectedCordSag,
											leftColumnSag, rightColumnSag,
											pedicleModelSag,
											debugMode);

	NIH_SpineSegmentation_ProjectColumnModel(img3D,
											 1, 512, 
											spinalCord3D,
											spineCenter, spineNormalX, spineNormalY, spineNormalBack,
											spineWidthLeft, spineWidthRight, spineWidthUp, spineWidthDown, 
											pedicleValley, backValley,
											projectedCordCor,
											leftColumnCor, rightColumnCor,
											pedicleModelCor,
											debugMode);

	// segment the disk after the partition
	NIH_SpineSegmentation_VertebralDiskDetection(img3D,
									   maskImg3D,
									   maskStruct,
									   segPara,
									   segInfo,
									   vertebraTemplate, 
									   spinalCord3D,
									   spineCenter, spineNormalX, spineNormalY, spineNormalBack,
									   spineWidthLeft, spineWidthRight, spineWidthUp, spineWidthDown, 
									   pedicleValley, backValley,
									   img2D_sag, img2D_cor,
									   img2D_sag_mask, img2D_cor_mask,
									   debugMode);

	if(viewMode==2) ChangeImageSlice(curSlice);
	if(viewMode==0) ChangeImageSlice_x(curSlice_x);
	if(viewMode==1) ChangeImageSlice_y(curSlice_y);

	vertHeightRuler->vertex_list.SetSize(0);
	vertHeightRuler->line_list.SetSize(0);

	NIH_OpenGL_CentralizeModel(display_id2, smoothedCord);
	NIH_OpenGL_Refresh_Display(display_id2);
	QApplication::restoreOverrideCursor();
}

void NIH_BoneSegmentation_Dlg::pushButton_BMD_ROI_clicked()
{
	if(maskImg3D==NULL || img3D==NULL || cordModel==NULL) return;

	QApplication::setOverrideCursor( QCursor(Qt::WaitCursor));
	// generate the export file name
	char screenshot_root[200];
	sprintf(screenshot_root, "%s\\screenshot\\", study_root);
	// create the direction if not exist
	QDir dcmdir = QString(screenshot_root);
	if(!dcmdir.exists())
	{
		dcmdir.mkdir(screenshot_root);
	}

	char export_fn[200];
	sprintf(export_fn,"%s%s_M.csv", screenshot_root, study_name);

	NIH_SpineSegmentation_BMD_ROI(img3D,
									   maskImg3D,
									   maskStruct,
									   segPara,
									   segInfo,
									   vertebraTemplate, 
									   spineWidthLeft, spineWidthRight, spineWidthUp, spineWidthDown, 
									   pedicleValley, backValley, roiSlices,
									   debugMode, export_fn);

	if(viewMode==2) ChangeImageSlice(curSlice);
	if(viewMode==0) ChangeImageSlice_x(curSlice_x);
	if(viewMode==1) ChangeImageSlice_y(curSlice_y);
	QApplication::restoreOverrideCursor();
}

void NIH_BoneSegmentation_Dlg::pushButton_BMD_ROI_3D_clicked()
{
	if(maskImg3D==NULL || img3D==NULL || cordModel==NULL) return;

	QApplication::setOverrideCursor( QCursor(Qt::WaitCursor));
	// generate the export file name
	char screenshot_root[200];
	sprintf(screenshot_root, "%s\\screenshot\\", study_root);
	// create the direction if not exist
	QDir dcmdir = QString(screenshot_root);
	if(!dcmdir.exists())
	{
		dcmdir.mkdir(screenshot_root);
	}

	char export_fn[200];
	sprintf(export_fn,"%s%s_M3D_5.csv", screenshot_root, study_name);

	NIH_SpineSegmentation_BMD_ROI_3D(img3D,
									   maskImg3D,
									   maskStruct,
									   segPara,
									   segInfo,
									   vertebraTemplate, 
									   spineCenter, spineNormalX, spineNormalY, spineNormalBack, 
									   spineWidthLeft, spineWidthRight, spineWidthUp, spineWidthDown, 
									   pedicleValley, backValley, roiSlices,
									   debugMode, export_fn);

	if(viewMode==2) ChangeImageSlice(curSlice);
	if(viewMode==0) ChangeImageSlice_x(curSlice_x);
	if(viewMode==1) ChangeImageSlice_y(curSlice_y);
	QApplication::restoreOverrideCursor();
}


void NIH_BoneSegmentation_Dlg::pushButton_BMD_ROI_L1L2_clicked()
{
	if(maskImg3D==NULL || img3D==NULL || cordModel==NULL) return;

	QApplication::setOverrideCursor( QCursor(Qt::WaitCursor));
	// generate the export file name
	char screenshot_root[200];
	sprintf(screenshot_root, "%s\\screenshot\\", study_root);
	// create the direction if not exist
	QDir dcmdir = QString(screenshot_root);
	if(!dcmdir.exists())
	{
		dcmdir.mkdir(screenshot_root);
	}

	char export_fn[200];
	sprintf(export_fn,"%s%s_M.csv", screenshot_root, study_name);

	// read manually selected slices
	char sel_fn[200];
	sprintf(sel_fn, "%s\\L1L2_selection.txt", study_root);
	FILE *fp_sel;
	int l1_slice, l2_slice;
	char l1_fn[200], l2_fn[200];
	CIS_Array_Image3D_short *l1_img=NULL, *l2_img=NULL;

	if((fp_sel=fopen(sel_fn, "r"))!=NULL)
	{
		fscanf(fp_sel, "%d %s\n", &l1_slice, l1_fn);
		fscanf(fp_sel, "%d %s\n", &l2_slice, l2_fn);
		fclose(fp_sel);
		l1_slice-=1;
		l2_slice-=1;
		l1_img = new CIS_Array_Image3D_short(l1_fn);
		if(l1_img->Num_Cols()==0) 
		{
			printf("Can not open file %s\n", l1_fn);
			return;
		}
		l2_img = new CIS_Array_Image3D_short(l2_fn);
		if(l2_img->Num_Cols()==0) 
		{
			printf("Can not open file %s\n", l2_fn);
			return;
		}
	}
	else
	{
		printf("Can not open selection file %s\n", sel_fn);
		return;
	}
	NIH_SpineSegmentation_BMD_ROI_L1L2(img3D,
									   maskImg3D,
									   maskStruct,
									   segPara,
									   segInfo,
									   vertebraTemplate, 
									   spineWidthLeft, spineWidthRight, spineWidthUp, spineWidthDown, 
									   pedicleValley, backValley, roiSlices,
									   debugMode, export_fn, l1_slice, l1_img, l2_slice, l2_img);

	if(l1_img) delete l1_img;
	if(l2_img) delete l2_img;

	if(viewMode==2) ChangeImageSlice(curSlice);
	if(viewMode==0) ChangeImageSlice_x(curSlice_x);
	if(viewMode==1) ChangeImageSlice_y(curSlice_y);
	QApplication::restoreOverrideCursor();
}

void NIH_BoneSegmentation_Dlg::pushButton_BMD_ROI_L1L2_3D_clicked()
{
	if(maskImg3D==NULL || img3D==NULL || cordModel==NULL) return;

	QApplication::setOverrideCursor( QCursor(Qt::WaitCursor));
	// generate the export file name
	char screenshot_root[200];
	sprintf(screenshot_root, "%s\\screenshot\\", study_root);
	// create the direction if not exist
	QDir dcmdir = QString(screenshot_root);
	if(!dcmdir.exists())
	{
		dcmdir.mkdir(screenshot_root);
	}

	char export_fn[200];
	sprintf(export_fn,"%s%s_M3D_6.csv", screenshot_root, study_name);

	// read manually selected slices
	char sel_fn[200];
	sprintf(sel_fn, "%s\\L1L2_3D_selection.txt", study_root);
	FILE *fp_sel;
	int l1_slice, l2_slice;
	char l1_fn[200], l2_fn[200];
	CIS_Array_Image3D_short *l1_img=NULL, *l2_img=NULL;

	if((fp_sel=fopen(sel_fn, "r"))!=NULL)
	{
		fscanf(fp_sel, "%d %s\n", &l1_slice, l1_fn);
		fscanf(fp_sel, "%d %s\n", &l2_slice, l2_fn);
		fclose(fp_sel);
		l1_slice-=1;
		l2_slice-=1;
		l1_img = new CIS_Array_Image3D_short(l1_fn);
		if(l1_img->Num_Cols()==0) 
		{
			printf("Can not open file %s\n", l1_fn);
			return;
		}
/*		l2_img = new CIS_Array_Image3D_short(l2_fn);
		if(l2_img->Num_Cols()==0) 
		{
			printf("Can not open file %s\n", l2_fn);
			return;
		}
*/	}
	else
	{
		printf("Can not open selection file %s\n", sel_fn);
		return;
	}
	NIH_SpineSegmentation_BMD_ROI_L1L2_3D(img3D,
									   maskImg3D,
									   maskStruct,
									   segPara,
									   segInfo,
									   vertebraTemplate, 
									   spineWidthLeft, spineWidthRight, spineWidthUp, spineWidthDown, 
									   pedicleValley, backValley, roiSlices,
									   debugMode, export_fn, l1_slice, l1_img, l2_slice, l1_img);

	if(l1_img) delete l1_img;
	if(l2_img) delete l2_img;

	if(viewMode==2) ChangeImageSlice(curSlice);
	if(viewMode==0) ChangeImageSlice_x(curSlice_x);
	if(viewMode==1) ChangeImageSlice_y(curSlice_y);
	QApplication::restoreOverrideCursor();
}



void NIH_BoneSegmentation_Dlg::pushButton_Rib_Detection_clicked()
{
	if(maskImg3D==NULL || img3D==NULL || cordModel==NULL) return;

	QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );
	NIH_SpineSegmentation_Rib_Detection(img3D,
									   maskImg3D,
									   maskStruct,
									   segPara,
									   segInfo,
									   vertebraTemplate, 
									   spineWidthLeft, spineWidthRight, spineWidthUp, spineWidthDown, 
									   pedicleValley, backValley,
									   debugMode);

	if(viewMode==2) ChangeImageSlice(curSlice);
	if(viewMode==0) ChangeImageSlice_x(curSlice_x);
	if(viewMode==1) ChangeImageSlice_y(curSlice_y);
	QApplication::restoreOverrideCursor();
}

void NIH_BoneSegmentation_Dlg::pushButton_Vertebra_Seg_clicked()
{
	if(maskImg3D==NULL || img3D==NULL) return;

	// experiments to measure vertebra height
	//
	VertebraHeightMeasurementSag();

	return;

	bool hasSeed=!flagShowPainting;

	// first locate a sub image for vertebra segmentation
	// segment T1 for now, apply the segmentation after the rib is detected
	IntVec3 bb1, bb2;
	IntVec3 seed;

	if(hasSeed)
	{
		seed = IntVec3(253, 286, 98);
		bb1 = IntVec3(207, 232, 71);
		bb2 = IntVec3(308, 391, 122);
	}
	else 
	{

		if(segInfo.t12_vertebra<0)	return;	// only do this when T12 is found
		int p,z, x, y, k;
		int sizex = img3D->Num_Cols();
		int sizey = img3D->Num_Rows();
		int sizez = img3D->Num_Levels();

		int sizexy = sizex*sizey;

		int curVert= spinBox_lesion->value();

		p=segInfo.t12_vertebra+curVert;
		if(p>=pedicleValley.GetSize()-1) return;


		// find the bounding box
		bb1.z = pedicleValley[p];
		bb2.z = pedicleValley[p+1];

		bb1.x = sizex+1;
		bb2.x = -1;
		bb1.y = sizey+1;
		bb2.y = -1;

		for(z=bb1.z; z<=bb2.z; z++)
		{
/*			if(segInfo.spineBound1[z].x<bb1.x) bb1.x=segInfo.spineBound1[z].x;
			if(segInfo.spineBound1[z].y<bb1.y) bb1.y=segInfo.spineBound1[z].y;

			if(segInfo.spineBound2[z].x>bb2.x) bb2.x=segInfo.spineBound2[z].x;
			if(segInfo.spineBound2[z].y>bb2.y) bb2.y=segInfo.spineBound2[z].y;
*/
			if(segInfo.diskBound1[z].x<bb1.x) bb1.x=segInfo.diskBound1[z].x;
			if(segInfo.diskBound1[z].y<bb1.y) bb1.y=segInfo.diskBound1[z].y;

			if(segInfo.diskBound2[z].x>bb2.x) bb2.x=segInfo.diskBound2[z].x;
			if(segInfo.diskBound2[z].y>bb2.y) bb2.y=segInfo.diskBound2[z].y;

			if(segInfo.sprocessBound2[z].y>bb2.y) bb2.y= segInfo.sprocessBound2[z].y;
		}

		seed.z = (bb1.z+bb2.z)/2;
		seed.x = segInfo.diskCenter[seed.z].x;
		seed.y = segInfo.diskCenter[seed.z].y;
		seed.z -=1;
		seed.y -=5;

		// add margin to the bounding box
		bb1.x -=10;
		bb1.y -=10;
		bb1.z -=10;
		bb2.x +=10;
		bb2.y +=10;
		bb2.z += 5;
	}

	SegmentSingleVertebra(seed, bb1, bb2);
	return;
}

int NIH_BoneSegmentation_Dlg::SegmentSingleVertebra(IntVec3 seed, IntVec3 bb1, IntVec3 bb2)
{
	if(maskImg3D==NULL || img3D==NULL) return CIS_ERROR;

	QApplication::setOverrideCursor( QCursor(Qt::WaitCursor));

	// adjust the bounding box so that each dimension has even number size
	//
	if((bb2.x-bb1.x)%2==0) bb2.x +=1;
	if((bb2.y-bb1.y)%2==0) bb2.y +=1;
	if((bb2.z-bb1.z)%2==0) bb2.z +=1;

	int p,z, x, y, k;
	int sizex = img3D->Num_Cols();
	int sizey = img3D->Num_Rows();
	int sizez = img3D->Num_Levels();

	int sizexy = sizex*sizey;

	short *maskA= maskImg3D->GetArray();

	// clear the mask
//	for(k=0;k<sizexy*sizez;k++) maskA[k]=maskStruct.mask_air;

	// visualize the bounding box and seed
/*	for(z=bb1.z; z<=bb2.z; z++)
	{
		k = z*sizex*sizey;
		for(y=bb1.y; y<=bb2.y; y++)
		{
			for(x=bb1.x; x<bb2.x; x++)
			{
				if(maskA[k+y*sizex+x]!=maskStruct.mask_spinalCord)
					maskA[k+y*sizex+x] = maskStruct.mask_rib;
			}
		}
	}

	maskA[seed.z*sizex*sizey+seed.y*sizex+seed.x] = maskStruct.mask_corticalBone;
	maskA[seed.z*sizex*sizey+seed.y*sizex+seed.x-1] = maskStruct.mask_corticalBone;
	maskA[seed.z*sizex*sizey+seed.y*sizex+seed.x+1] = maskStruct.mask_corticalBone;
	maskA[seed.z*sizex*sizey+(seed.y-1)*sizex+seed.x] = maskStruct.mask_corticalBone;
	maskA[seed.z*sizex*sizey+(seed.y+1)*sizex+seed.x] = maskStruct.mask_corticalBone;
*/
	// get the sub image for vertebra segmentation
	CIS_Array_Image3D_short *subImgIn=NULL, *subImgOut=NULL;
	int sub_sizex, sub_sizey, sub_sizez;
	short *subOut;
	img3D->SubImage(subImgIn, bb1.x, bb1.y, bb1.z, bb2.x, bb2.y, bb2.z);
	img3D->SubImage(subImgOut, bb1.x, bb1.y, bb1.z, bb2.x, bb2.y, bb2.z);
	seed.x -= bb1.x;
	seed.y -= bb1.y;
	seed.z -= bb1.z;
	sub_sizex = subImgIn->Num_Cols();
	sub_sizey = subImgIn->Num_Rows();
	sub_sizez = subImgIn->Num_Levels();
	subOut = subImgOut->GetArray();

	// apply sovira's vertebra segmentation
	SegmentVertebralBody(subImgIn, subImgOut, seed.x, seed.y, seed.z);

	// put the results back to the mask
	for(z=0, k=0; z<sub_sizez; z++)
	{
		for(y=0; y<sub_sizey; y++)
		{
			for(x=0; x<sub_sizex; x++, k++)
			{
				if(subOut[k]==0) 
				{
//					if(maskA[(z+bb1.z)*sizexy+(y+bb1.y)*sizex+x+bb1.x]!=maskStruct.mask_corticalBone)
						maskA[(z+bb1.z)*sizexy+(y+bb1.y)*sizex+x+bb1.x]=maskStruct.mask_vertebra;
					subOut[k]=100;
				} else subOut[k]=0;
			}
		}
	}

	// generate surface
	spineSurf->RemoveAllMesh();
	CIS_Algo_MarchingCube(50, -1, -1, 2, 1, subImgOut, NULL, spineSurf,-1, -1, -1, -1, -1, -1);
	TriangleMesh_ApplyTranslation(spineSurf, Vec3(bb1.x*img3D->Get_Pixel_SizeX(), bb1.y*img3D->Get_Pixel_SizeY(), 0));

	NIH_OpenGL_CentralizeModel(display_id2, spineSurf);
	spineSurf->SetCurrAlpha(1.0);
	smoothedCord->SetDrawMode(HIDE_MODE);
	cordModel->SetDrawMode(HIDE_MODE);
	NIH_OpenGL_Refresh_Display(display_id2);

	delete subImgIn;
	delete subImgOut;

	QApplication::restoreOverrideCursor();
	return CIS_OK;
}


int NIH_BoneSegmentation_Dlg::VertebraHeightMeasurement()
{
	if(maskImg3D==NULL || img3D==NULL) return CIS_ERROR;
	if(segInfo.t12_vertebra<0)	return CIS_ERROR;	// only do this when T12 is found
	int p,z, x, y, k;
	int sizex = img3D->Num_Cols();
	int sizey = img3D->Num_Rows();
	int sizez = img3D->Num_Levels();
	float px = img3D->Get_Pixel_SizeX();
	float py = img3D->Get_Pixel_SizeY();
	float pz = img3D->Get_Pixel_SizeZ();

	int sizexy = sizex*sizey;

	short *maskA= maskImg3D->GetArray();

	int curVert= spinBox_lesion->value();
	p=segInfo.t12_vertebra+curVert;
	if(p>=pedicleValley.GetSize()-1) return CIS_ERROR;

	// find the bounding box and the spinal cord center
	Vec3 cordCenter1, cordCenter2;
	IntVec3 bb1, bb2;

	bb1.x = bb1.y = 512;
	bb2.x = bb2.y = 0;
	bb1.z = pedicleValley[p];
	bb2.z = pedicleValley[p+1];

	for(z=bb1.z; z<=bb2.z; z++)
	{
		if(segInfo.diskBound1[z].x<bb1.x) bb1.x=segInfo.diskBound1[z].x;
		if(segInfo.diskBound1[z].y<bb1.y) bb1.y=segInfo.diskBound1[z].y;

		if(segInfo.diskBound2[z].x>bb2.x) bb2.x=segInfo.diskBound2[z].x;
		if(segInfo.diskBound2[z].y>bb2.y) bb2.y=segInfo.diskBound2[z].y;
	}

	cordCenter1.x = segInfo.cordCenter[bb1.z].x*px;
	cordCenter1.y = segInfo.cordCenter[bb1.z].y*py;

	cordCenter2.x = segInfo.cordCenter[bb2.z].x*px;
	cordCenter2.y = segInfo.cordCenter[bb2.z].y*py;

	cordCenter1.z = img3D->GetSlicePosition(bb1.z);
	cordCenter2.z = img3D->GetSlicePosition(bb2.z);

	// 
	int start_x = (segInfo.cordCenter[bb1.z].x+segInfo.cordCenter[bb2.z].x)/2;
	int start_y = (segInfo.cordCenter[bb1.z].y+segInfo.cordCenter[bb2.z].y)/2;
	int start_z = (bb1.z+bb2.z)/2;

	Vec2 vertHeightDir2D = Vec2(cordCenter2.y-cordCenter1.y, (bb2.z-bb1.z)*pz).normalize();

	// add margin to the bounding box
	bb1.x -=10;
	bb1.y -=10;
	bb1.z -=10;
	bb2.x +=10;
	bb2.y +=10;
	bb2.z += 5;

	Vec3 stPos, curPos;
	IntVec3 endPlate1, endPlate2;
	IntVec3 curPosi;
	intVec3DynArray plateArray1, plateArray2;
	float estHeight = (cordCenter2-cordCenter1).len();

	plateArray1.SetSize(0);
	plateArray2.SetSize(0);

	// search for the vertebra height
	for(y=start_y; y>bb1.y; y--)
	{
		stPos = Vec3(start_x*px, y*py, start_z*pz);

		// search in both directions for end plates
		curPos = stPos;
		curPosi.x = (int)(curPos.x/px);
		curPosi.y = (int)(curPos.y/py);
		curPosi.z = (int)(curPos.z/pz);
		endPlate1 = curPosi;
		endPlate2 = curPosi;
		while(1)
		{
			curPosi.x = (int)(curPos.x/px);
			curPosi.y = (int)(curPos.y/py);
			curPosi.z = (int)(curPos.z/pz);

			if(maskA[curPosi.z*sizexy+curPosi.y*sizex+curPosi.x]==maskStruct.mask_spinalCord)
			{
				// hit the spinal cord, should be discarded
				break;
			}

			if(maskA[curPosi.z*sizexy+curPosi.y*sizex+curPosi.x]!=maskStruct.mask_vertebra)
			{
				endPlate1 = curPosi;
				break;
			}
			curPos.y += vertHeightDir2D.x;
			curPos.z += vertHeightDir2D.y;
		}

		// search in opposite direction
		curPos = stPos;
		while(1)
		{
			curPosi.x = (int)(curPos.x/px);
			curPosi.y = (int)(curPos.y/py);
			curPosi.z = (int)(curPos.z/pz);

			if(maskA[curPosi.z*sizexy+curPosi.y*sizex+curPosi.x]==maskStruct.mask_spinalCord)
			{
				// hit the spinal cord, should be discarded
				break;
			}

			if(maskA[curPosi.z*sizexy+curPosi.y*sizex+curPosi.x]!=maskStruct.mask_vertebra)
			{
				endPlate2 = curPosi;
				break;
			}
			curPos.y -= vertHeightDir2D.x;
			curPos.z -= vertHeightDir2D.y;
		}

		if((endPlate2-endPlate1).len()>estHeight/2)
		{
			plateArray1.Add(endPlate1);
			plateArray2.Add(endPlate2);
		}
	}

	// set up the height ruler model
	vertHeightRuler->SetDrawMode(LINE_MODE);
	vertHeightRuler->vertex_list.SetSize(0);
	vertHeightRuler->line_list.SetSize(0);
	vertHeightRuler->plot_radius = 2;

	float vsizey, vsizez;
	int isizey, isizez;
	vsizey = img3D->Get_Volume_SizeY();
	vsizez = img3D->Get_Volume_SizeZ();

	int window_size = 512;
	if(vsizey>vsizez)
	{
		isizey = window_size;
		isizez = (int)((vsizez/vsizey)*window_size);
	}
	else 
	{
		isizey = (int)((vsizey/vsizez)*window_size);
		isizez = window_size;
	}

	float projSize_x = isizey;
	float projSize_y = isizez;
	float projLoc_x = (window_size-isizey)/2;
	float projLoc_y = (window_size-isizez)/2;

	float fx = projSize_x/vsizey;
	float fy = projSize_y/vsizez;

	// Add the centerline first
	Vec2 projA1, projA2;
	int count=0;
	projA1.x = segInfo.cordCenter[pedicleValley[p]].y*py*fx+projLoc_x;
	projA1.y = pedicleValley[p]*pz*fy+projLoc_y;
	projA2.x = segInfo.cordCenter[pedicleValley[p+1]].y*py*fx+projLoc_x;
	projA2.y = pedicleValley[p+1]*pz*fy+projLoc_y;
	vertHeightRuler->vertex_list.Add(projA1);
	vertHeightRuler->vertex_list.Add(projA2);
	vertHeightRuler->line_list.Add(IntVec2(count, count+1));
	count+=2;

	int rulerCount = plateArray1.GetSize();

	if(rulerCount<10) return CIS_ERROR;


	int largestInd;
	double largestLen;
	// add first link
	// pick the largest among the first three
	//
	largestLen=0;
	for(k=0; k<3; k++)
	{
		if((plateArray1[k]-plateArray2[k]).len()>largestLen)
		{
			largestLen = (plateArray1[k]-plateArray2[k]).len();
			largestInd = k;
		}
	}
	projA1.x = plateArray1[largestInd].y*py*fx+projLoc_x;
	projA1.y = plateArray1[largestInd].z*pz*fy+projLoc_y;
	projA2.x = plateArray2[largestInd].y*py*fx+projLoc_x;
	projA2.y = plateArray2[largestInd].z*pz*fy+projLoc_y;
	vertHeightRuler->vertex_list.Add(projA1);
	vertHeightRuler->vertex_list.Add(projA2);
	vertHeightRuler->line_list.Add(IntVec2(count, count+1));
	count+=2;

	// add middle link
	// pick the largest among the middle three
	//
	largestLen=0;
	for(k=rulerCount/2-1; k<=rulerCount/2+1; k++)
	{
		if((plateArray1[k]-plateArray2[k]).len()>largestLen)
		{
			largestLen = (plateArray1[k]-plateArray2[k]).len();
			largestInd = k;
		}
	}
	projA1.x = plateArray1[largestInd].y*py*fx+projLoc_x;
	projA1.y = plateArray1[largestInd].z*pz*fy+projLoc_y;
	projA2.x = plateArray2[largestInd].y*py*fx+projLoc_x;
	projA2.y = plateArray2[largestInd].z*pz*fy+projLoc_y;
	vertHeightRuler->vertex_list.Add(projA1);
	vertHeightRuler->vertex_list.Add(projA2);
	vertHeightRuler->line_list.Add(IntVec2(count, count+1));
	count+=2;

	// add last link
	// pick the largest among the last three
	//
	largestLen=0;
	for(k=rulerCount-3; k<rulerCount; k++)
	{
		if((plateArray1[k]-plateArray2[k]).len()>largestLen)
		{
			largestLen = (plateArray1[k]-plateArray2[k]).len();
			largestInd = k;
		}
	}
	projA1.x = plateArray1[largestInd].y*py*fx+projLoc_x;
	projA1.y = plateArray1[largestInd].z*pz*fy+projLoc_y;
	projA2.x = plateArray2[largestInd].y*py*fx+projLoc_x;
	projA2.y = plateArray2[largestInd].z*pz*fy+projLoc_y;
	vertHeightRuler->vertex_list.Add(projA1);
	vertHeightRuler->vertex_list.Add(projA2);
	vertHeightRuler->line_list.Add(IntVec2(count, count+1));
	count+=2;

	// add 3D line models
	vertHeightDir->SetDrawMode(LINE_MODE);
	vertHeightDir->SetEndPoints(cordCenter1, cordCenter2);

	NIH_OpenGL_Refresh_Display(display_id2);
	return CIS_OK;
}


int NIH_BoneSegmentation_Dlg::VertebraHeightMeasurement1()
{
	if(maskImg3D==NULL || img3D==NULL) return CIS_ERROR;
	if(segInfo.t12_vertebra<0)	return CIS_ERROR;	// only do this when T12 is found
	int p,z, x, y, k;
	int sizex = img3D->Num_Cols();
	int sizey = img3D->Num_Rows();
	int sizez = img3D->Num_Levels();
	float px = img3D->Get_Pixel_SizeX();
	float py = img3D->Get_Pixel_SizeY();
	float pz = img3D->Get_Pixel_SizeZ();

	int sizexy = sizex*sizey;

	short *maskA= maskImg3D->GetArray();

	int curVert= spinBox_lesion->value();
	p=segInfo.t12_vertebra+curVert;
	if(p>=pedicleValley.GetSize()-1) return CIS_ERROR;

	// find the bounding box and the spinal cord center
	Vec3 cordCenter1, cordCenter2;
	IntVec3 bb1, bb2;

	bb1.x = bb1.y = 512;
	bb2.x = bb2.y = 0;
	bb1.z = pedicleValley[p];
	bb2.z = pedicleValley[p+1];

	for(z=bb1.z; z<=bb2.z; z++)
	{
		if(segInfo.diskBound1[z].x<bb1.x) bb1.x=segInfo.diskBound1[z].x;
		if(segInfo.diskBound1[z].y<bb1.y) bb1.y=segInfo.diskBound1[z].y;

		if(segInfo.diskBound2[z].x>bb2.x) bb2.x=segInfo.diskBound2[z].x;
		if(segInfo.diskBound2[z].y>bb2.y) bb2.y=segInfo.diskBound2[z].y;
	}

	cordCenter1.x = segInfo.cordCenter[bb1.z].x*px;
	cordCenter1.y = segInfo.cordCenter[bb1.z].y*py;

	cordCenter2.x = segInfo.cordCenter[bb2.z].x*px;
	cordCenter2.y = segInfo.cordCenter[bb2.z].y*py;

	cordCenter1.z = img3D->GetSlicePosition(bb1.z);
	cordCenter2.z = img3D->GetSlicePosition(bb2.z);

	// 
	int start_x = (segInfo.cordCenter[bb1.z].x+segInfo.cordCenter[bb2.z].x)/2;
	int start_y = (segInfo.cordCenter[bb1.z].y+segInfo.cordCenter[bb2.z].y)/2;
	int start_z = (bb1.z+bb2.z)/2;

	Vec2 vertHeightDir2D = Vec2(cordCenter2.y-cordCenter1.y, (bb2.z-bb1.z)*pz).normalize();

	// add margin to the bounding box
	bb1.x -=10;
	bb1.y -=10;
	bb1.z -=10;
	bb2.x +=10;
	bb2.y +=10;
	bb2.z += 5;

	Vec3 stPos, curPos;
	IntVec3 endPlate1, endPlate2;
	IntVec3 curPosi;
	intVec3DynArray plateArray1, plateArray2;
	float estHeight = (cordCenter2-cordCenter1).len();

	plateArray1.SetSize(0);
	plateArray2.SetSize(0);

	// search for the vertebra height
	for(y=start_y; y>bb1.y; y--)
	{
		stPos = Vec3(start_x*px, y*py, start_z*pz);

		// search in both directions for end plates
		curPos = stPos;
		curPosi.x = (int)(curPos.x/px);
		curPosi.y = (int)(curPos.y/py);
		curPosi.z = (int)(curPos.z/pz);
		if(maskA[curPosi.z*sizexy+curPosi.y*sizex+curPosi.x]!=maskStruct.mask_corticalBone &&
			maskA[curPosi.z*sizexy+curPosi.y*sizex+curPosi.x]!=maskStruct.mask_spongyBone) continue;
		
		endPlate1 = curPosi;
		endPlate2 = curPosi;
		while(1)
		{
			curPosi.x = (int)(curPos.x/px);
			curPosi.y = (int)(curPos.y/py);
			curPosi.z = (int)(curPos.z/pz);

			if(curPosi.z<bb1.z || curPosi.z>bb2.z) break;

			if(maskA[curPosi.z*sizexy+curPosi.y*sizex+curPosi.x]==maskStruct.mask_spinalCord)
			{
				// hit the spinal cord, should be discarded
				break;
			}

			if(maskA[curPosi.z*sizexy+curPosi.y*sizex+curPosi.x]==maskStruct.mask_vertebralDisk)
			{
				endPlate1 = curPosi;
				break;
			}
			curPos.y += vertHeightDir2D.x;
			curPos.z += vertHeightDir2D.y;
		}

		// search in opposite direction
		curPos = stPos;
		while(1)
		{
			curPosi.x = (int)(curPos.x/px);
			curPosi.y = (int)(curPos.y/py);
			curPosi.z = (int)(curPos.z/pz);

			if(curPosi.z<bb1.z || curPosi.z>bb2.z) break;
			if(maskA[curPosi.z*sizexy+curPosi.y*sizex+curPosi.x]==maskStruct.mask_spinalCord)
			{
				// hit the spinal cord, should be discarded
				break;
			}

			if(maskA[curPosi.z*sizexy+curPosi.y*sizex+curPosi.x]==maskStruct.mask_vertebralDisk)
			{
				endPlate2 = curPosi;
				break;
			}
			curPos.y -= vertHeightDir2D.x;
			curPos.z -= vertHeightDir2D.y;
		}

		if((endPlate2-endPlate1).len()>estHeight/2)
		{
			plateArray1.Add(endPlate1);
			plateArray2.Add(endPlate2);
		}
	}

	// set up the height ruler model
	vertHeightRuler->SetDrawMode(LINE_MODE);
//	vertHeightRuler->vertex_list.SetSize(0);
//	vertHeightRuler->line_list.SetSize(0);
	vertHeightRuler->plot_radius = 2;

	float vsizey, vsizez;
	int isizey, isizez;
	vsizey = img3D->Get_Volume_SizeY();
	vsizez = img3D->Get_Volume_SizeZ();

	int window_size = 512;
	if(vsizey>vsizez)
	{
		isizey = window_size;
		isizez = (int)((vsizez/vsizey)*window_size);
	}
	else 
	{
		isizey = (int)((vsizey/vsizez)*window_size);
		isizez = window_size;
	}

	float projSize_x = isizey;
	float projSize_y = isizez;
	float projLoc_x = (window_size-isizey)/2;
	float projLoc_y = (window_size-isizez)/2;

	float fx = projSize_x/vsizey;
	float fy = projSize_y/vsizez;

	Vec2 projA1, projA2;
	int count=vertHeightRuler->vertex_list.GetSize();

/*	// Add the centerline first
	projA1.x = segInfo.cordCenter[pedicleValley[p]].y*py*fx+projLoc_x;
	projA1.y = pedicleValley[p]*pz*fy+projLoc_y;
	projA2.x = segInfo.cordCenter[pedicleValley[p+1]].y*py*fx+projLoc_x;
	projA2.y = pedicleValley[p+1]*pz*fy+projLoc_y;
	vertHeightRuler->vertex_list.Add(projA1);
	vertHeightRuler->vertex_list.Add(projA2);
	vertHeightRuler->line_list.Add(IntVec2(count, count+1));
	count+=2;
*/
	int rulerCount = plateArray1.GetSize();

	if(rulerCount<10) return CIS_ERROR;


	int largestInd, meanInd;
	double largestLen, meanLen, closestLen;
	// add first link
	// pick the largest among the first five
	// pick the one closest to the mean
	//
	largestLen=0;
	meanLen=0;
	for(k=0; k<5; k++)
	{
		if((plateArray1[k]-plateArray2[k]).len()>largestLen)
		{
			largestLen = (plateArray1[k]-plateArray2[k]).len();
			largestInd = k;
		}
		meanLen += (plateArray1[k]-plateArray2[k]).len();
	}
	meanLen /= 5;
	closestLen = meanLen;
	for(k=0; k<5; k++)
	{
		if(fabs((plateArray1[k]-plateArray2[k]).len()-meanLen)<closestLen) 
		{
			closestLen=fabs((plateArray1[k]-plateArray2[k]).len()-meanLen);
			meanInd=k;
		}
	}

	projA1.x = plateArray1[meanInd].y*py*fx+projLoc_x;
	projA1.y = plateArray1[meanInd].z*pz*fy+projLoc_y;
	projA2.x = plateArray2[meanInd].y*py*fx+projLoc_x;
	projA2.y = plateArray2[meanInd].z*pz*fy+projLoc_y;
	vertHeightRuler->vertex_list.Add(projA1);
	vertHeightRuler->vertex_list.Add(projA2);
	vertHeightRuler->line_list.Add(IntVec2(count, count+1));
	count+=2;

	// add middle link
	// pick the largest among the middle three
	//
	largestLen=0;
	meanLen=0;
	for(k=rulerCount/2-2; k<=rulerCount/2+2; k++)
	{
		if((plateArray1[k]-plateArray2[k]).len()>largestLen)
		{
			largestLen = (plateArray1[k]-plateArray2[k]).len();
			largestInd = k;
		}
		meanLen += (plateArray1[k]-plateArray2[k]).len();
	}
	meanLen /= 5;
	closestLen = meanLen;
	for(k=rulerCount/2-2; k<=rulerCount/2+2; k++)
	{
		if(fabs((plateArray1[k]-plateArray2[k]).len()-meanLen)<closestLen) 
		{
			closestLen=fabs((plateArray1[k]-plateArray2[k]).len()-meanLen);
			meanInd=k;
		}
	}

	projA1.x = plateArray1[meanInd].y*py*fx+projLoc_x;
	projA1.y = plateArray1[meanInd].z*pz*fy+projLoc_y;
	projA2.x = plateArray2[meanInd].y*py*fx+projLoc_x;
	projA2.y = plateArray2[meanInd].z*pz*fy+projLoc_y;
	vertHeightRuler->vertex_list.Add(projA1);
	vertHeightRuler->vertex_list.Add(projA2);
	vertHeightRuler->line_list.Add(IntVec2(count, count+1));
	count+=2;


	// add last link
	// pick the largest among the last three
	//
	largestLen=0;
	meanLen=0;
	for(k=rulerCount-5; k<rulerCount; k++)
	{
		if((plateArray1[k]-plateArray2[k]).len()>largestLen)
		{
			largestLen = (plateArray1[k]-plateArray2[k]).len();
			largestInd = k;
		}
		meanLen += (plateArray1[k]-plateArray2[k]).len();
	}
	meanLen /= 5;
	closestLen = meanLen;
	for(k=rulerCount-5; k<rulerCount; k++)
	{
		if(fabs((plateArray1[k]-plateArray2[k]).len()-meanLen)<closestLen) 
		{
			closestLen=fabs((plateArray1[k]-plateArray2[k]).len()-meanLen);
			meanInd=k;
		}
	}
	projA1.x = plateArray1[meanInd].y*py*fx+projLoc_x;
	projA1.y = plateArray1[meanInd].z*pz*fy+projLoc_y;
	projA2.x = plateArray2[meanInd].y*py*fx+projLoc_x;
	projA2.y = plateArray2[meanInd].z*pz*fy+projLoc_y;
	vertHeightRuler->vertex_list.Add(projA1);
	vertHeightRuler->vertex_list.Add(projA2);
	vertHeightRuler->line_list.Add(IntVec2(count, count+1));
	count+=2;

	// add 3D line models
	vertHeightDir->SetDrawMode(LINE_MODE);
	vertHeightDir->SetEndPoints(cordCenter1, cordCenter2);

	NIH_OpenGL_Refresh_Display(display_id);
	NIH_OpenGL_Refresh_Display(display_id2);
	return CIS_OK;
}


int NIH_BoneSegmentation_Dlg::VertebraHeightMeasurementSag()
{
	if(maskImg3D==NULL || img3D==NULL || img2D_sag_mask==NULL) return CIS_ERROR;
	if(segInfo.bound1.z<0) return CIS_ERROR;
///	if(segInfo.t12_vertebra<0)	return CIS_ERROR;	// only do this when T12 is found
	if(segInfo.t12_vertebra<0)	segInfo.t12_vertebra=0;	// only do this when T12 is found

	int p,z, x, y, k;
	int sizex = img3D->Num_Cols();
	int sizey = img3D->Num_Rows();
	int sizez = img3D->Num_Levels();
	float px = img3D->Get_Pixel_SizeX();
	float py = img3D->Get_Pixel_SizeY();
	float pz = img3D->Get_Pixel_SizeZ();

	int sizexy = sizex*sizey;

	short *maskA= img2D_sag_mask->GetArray();

	// generate the export file name
	char screenshot_root[200];
	sprintf(screenshot_root, "%s\\screenshot\\", study_root);
	// create the direction if not exist
	QDir dcmdir = QString(screenshot_root);
	if(!dcmdir.exists())
	{
		dcmdir.mkdir(screenshot_root);
	}

	char export_fn[200];
	sprintf(export_fn,"%s%s_H.csv", screenshot_root, study_name);

	FILE *e_fp;
	if((e_fp = fopen(export_fn, "w"))!=NULL)
	{
		// print the header
		fprintf(e_fp,"t12_h1,t12_h2,t12_h3,l1_h1,l1_h2,l1_h3,");
		fprintf(e_fp,"l2_h1,l2_h2,l2_h3,l3_h1,l3_h2,l3_h3,");
		fprintf(e_fp,"l4_h1,l4_h2,l4_h3,l5_h1,l5_h2,l5_h3,\n");
	}


	int curVert= spinBox_lesion->value();
	for(p=segInfo.t12_vertebra; p<pedicleValley.GetSize()-1; p++)
	{

	int z1, z2;

	z1 = pedicleValley[p];
	z2 = pedicleValley[p+1];

	intDynArray filledZ;
	filledZ.SetSize(sizez);
	for(z=0; z<sizez; z++) 
	{
		filledZ[z]=-1;
	}
		
	for(k=0; k<spinalCord3D.GetSize(); k++)
	{
		z = img3D->GetSliceLevel(spinalCord3D[k].z,1);

		if(filledZ[z]==1 && z>0 && filledZ[z-1]==-1) z=z-1;
		while(z<sizez-1 && filledZ[z]!=-1) z++;
		filledZ[z] = k;
	}

	// 
	int z12 = (z1+z2)/2;
	int k1 = filledZ[z1];
	int k2 = filledZ[z2];
	int k12 = filledZ[z12];

	if(k1==-1 || k2==-1 || k12==-1)
	{
		// shouldn't happen, need double check case 010M.
		continue;
	};

	Vec2 stPos1, stPos2, stPos12;
	stPos1 = Vec2(spinalCord3D[k1].y,z1*pz);
	stPos2 = Vec2(spinalCord3D[k2].y,z2*pz);
	stPos12 = Vec2(spinalCord3D[k12].y,z12*pz);

	Vec2 normal1, normal2, normal12;
	normal1.x = spineNormalX[z1].y;
	normal1.y = spineNormalX[z1].z;
	if(normal1.x>0) normal1=-normal1;
	normal2.x = spineNormalX[z2].y;
	normal2.y = spineNormalX[z2].z;
	if(normal2.x>0) normal2=-normal2;
//	normal12.x = spineNormalX[z12].y;
//	normal12.y = spineNormalX[z12].z;
//	if(normal12.x>0) normal12=-normal12;
	normal12 = (normal1+normal2).normalize();


	Vec2 vertHeightDir2D = (stPos2-stPos1).normalize();
	float estHeight = (stPos2-stPos1).len();
	float estWidth = spineWidthUp[z12];

	Vec2 stPos, curPos, boundPos1, boundPos2;
	IntVec2 endPlate1, endPlate2;
	IntVec2 curPosi, lastStPosi, lastPosi;
	intVec2DynArray plateArray1, plateArray2;

	plateArray1.SetSize(0);
	plateArray2.SetSize(0);

	// search for the vertebra height
	float step;
	int ck;
	for(step=0; step<estWidth*1.5; step+=py/2)
	{
		stPos = stPos12+normal12*step;

		// search in both directions for end plates
		curPos = stPos;
		curPosi.x = (int)(curPos.x/py+0.5);
		curPosi.y = (int)(curPos.y/pz+0.5);
		if(curPosi==lastStPosi) continue;
		lastStPosi = curPosi;

		ck = curPosi.y*sizey+curPosi.x;
		if(curPosi.y<0 || curPosi.y>sizez-1) break;

///		if(maskA[ck]==maskStruct.mask_body) break;
		
		if(maskA[ck]!=maskStruct.mask_corticalBone &&
			maskA[ck]!=maskStruct.mask_spongyBone) continue;
		
		endPlate1 = curPosi;
		endPlate2 = curPosi;

		boundPos1 = stPos1+normal1*step;
		boundPos2 = stPos2+normal2*step;

		while(1)
		{
			lastPosi=curPosi;
			curPosi.x = (int)(curPos.x/py+0.5);
			curPosi.y = (int)(curPos.y/pz+0.5);
			ck = curPosi.y*sizey+curPosi.x;
			if(curPosi.y<0 || curPosi.y>sizez-1) break;

//			if(maskA[ck]==maskStruct.mask_spinalCord)
//			{
				// hit the spinal cord, should be discarded
//				break;
//			}

			if(curPos.y<boundPos1.y+pz/2 || curPos.y>boundPos2.y-pz/2)
			{
				endPlate1 = lastPosi;
				break;
			}
			if(maskA[ck]==maskStruct.mask_vertebralDisk ||
				maskA[ck]==maskStruct.mask_body)
			{
				endPlate1 = lastPosi;
				break;
			}
			curPos.x += vertHeightDir2D.x;
			curPos.y += vertHeightDir2D.y;
		}

		// search in opposite direction
		curPos = stPos;
		while(1)
		{
			lastPosi=curPosi;
			curPosi.x = (int)(curPos.x/py+0.5);
			curPosi.y = (int)(curPos.y/pz+0.5);
			ck = curPosi.y*sizey+curPosi.x;
			if(curPosi.y<0 || curPosi.y>sizez-1) break;

//			if(maskA[ck]==maskStruct.mask_spinalCord)
//			{
				// hit the spinal cord, should be discarded
//				break;
//			}

			if(curPos.y<boundPos1.y+pz/2 || curPos.y>boundPos2.y-pz/2)
			{
				endPlate2 = lastPosi;
				break;
			}
			if(maskA[ck]==maskStruct.mask_vertebralDisk ||
				maskA[ck]==maskStruct.mask_body)
			{
				endPlate2 = lastPosi;
				break;
			}

			curPos.x -= vertHeightDir2D.x;
			curPos.y -= vertHeightDir2D.y;
		}

		if((endPlate2-endPlate1).len()>estHeight/2)
		{
			plateArray1.Add(endPlate1);
			plateArray2.Add(endPlate2);
		}
	}

	// set up the height ruler model
	vertHeightRuler->SetDrawMode(LINE_MODE);
	vertHeightRuler->plot_radius = 2;

	float height1=0, height2=0, height3=0;

	float vsizey, vsizez;
	int isizey, isizez;
	vsizey = img3D->Get_Volume_SizeY();
	vsizez = img3D->Get_Volume_SizeZ();

	int window_size = 512;
	if(vsizey>vsizez)
	{
		isizey = window_size;
		isizez = (int)((vsizez/vsizey)*window_size);
	}
	else 
	{
		isizey = (int)((vsizey/vsizez)*window_size);
		isizez = window_size;
	}

	float projSize_x = isizey;
	float projSize_y = isizez;
	float projLoc_x = (window_size-isizey)/2;
	float projLoc_y = (window_size-isizez)/2;

	float fx = projSize_x/vsizey;
	float fy = projSize_y/vsizez;
	

	Vec2 projA1, projA2;
	int count=vertHeightRuler->vertex_list.GetSize();

/*	// Add the centerline first
	projA1.x = segInfo.cordCenter[pedicleValley[p]].y*py*fx+projLoc_x;
	projA1.y = pedicleValley[p]*pz*fy+projLoc_y;
	projA2.x = segInfo.cordCenter[pedicleValley[p+1]].y*py*fx+projLoc_x;
	projA2.y = pedicleValley[p+1]*pz*fy+projLoc_y;
	vertHeightRuler->vertex_list.Add(projA1);
	vertHeightRuler->vertex_list.Add(projA2);
	vertHeightRuler->line_list.Add(IntVec2(count, count+1));
	count+=2;
*/
	int rulerCount = plateArray1.GetSize();

	if(rulerCount<10) return CIS_ERROR;


	int largestInd, meanInd;
	double largestLen, meanLen, closestLen;
	// add first link
	// pick the largest among the first five
	// pick the one closest to the mean
	//
	largestLen=0;
	meanLen=0;
	for(k=0; k<5; k++)
	{
		if((plateArray1[k]-plateArray2[k]).len()>largestLen)
		{
			largestLen = (plateArray1[k]-plateArray2[k]).len();
			largestInd = k;
		}
		meanLen += (plateArray1[k]-plateArray2[k]).len();
	}
	meanLen /= 5;
	closestLen = meanLen;
	for(k=0; k<5; k++)
	{
		if(fabs((plateArray1[k]-plateArray2[k]).len()-meanLen)<closestLen) 
		{
			closestLen=fabs((plateArray1[k]-plateArray2[k]).len()-meanLen);
			meanInd=k;
		}
	}

	projA1.x = plateArray1[largestInd].x*py*fx+projLoc_x;
	projA1.y = plateArray1[largestInd].y*pz*fy+projLoc_y;
	projA2.x = plateArray2[largestInd].x*py*fx+projLoc_x;
	projA2.y = plateArray2[largestInd].y*pz*fy+projLoc_y;
	vertHeightRuler->vertex_list.Add(projA1);
	vertHeightRuler->vertex_list.Add(projA2);
	vertHeightRuler->line_list.Add(IntVec2(count, count+1));
	count+=2;
	height1 = largestLen;

	// add middle link
	// pick the largest among the middle three
	//
	largestLen=0;
	meanLen=0;
	for(k=rulerCount/2-2; k<=rulerCount/2+2; k++)
	{
		if((plateArray1[k]-plateArray2[k]).len()>largestLen)
		{
			largestLen = (plateArray1[k]-plateArray2[k]).len();
			largestInd = k;
		}
		meanLen += (plateArray1[k]-plateArray2[k]).len();
	}
	meanLen /= 5;
	closestLen = meanLen;
	for(k=rulerCount/2-2; k<=rulerCount/2+2; k++)
	{
		if(fabs((plateArray1[k]-plateArray2[k]).len()-meanLen)<closestLen) 
		{
			closestLen=fabs((plateArray1[k]-plateArray2[k]).len()-meanLen);
			meanInd=k;
		}
	}

	projA1.x = plateArray1[meanInd].x*py*fx+projLoc_x;
	projA1.y = plateArray1[meanInd].y*pz*fy+projLoc_y;
	projA2.x = plateArray2[meanInd].x*py*fx+projLoc_x;
	projA2.y = plateArray2[meanInd].y*pz*fy+projLoc_y;
	vertHeightRuler->vertex_list.Add(projA1);
	vertHeightRuler->vertex_list.Add(projA2);
	vertHeightRuler->line_list.Add(IntVec2(count, count+1));
	count+=2;
	height2 = (float)(plateArray1[meanInd]-plateArray2[meanInd]).len();


	// add last link
	// pick the largest among the last three
	//
	largestLen=0;
	meanLen=0;
	for(k=rulerCount-5; k<rulerCount; k++)
	{
		if((plateArray1[k]-plateArray2[k]).len()>largestLen)
		{
			largestLen = (plateArray1[k]-plateArray2[k]).len();
			largestInd = k;
		}
		meanLen += (plateArray1[k]-plateArray2[k]).len();
	}
	meanLen /= 5;
	closestLen = meanLen;
	for(k=rulerCount-5; k<rulerCount; k++)
	{
		if(fabs((plateArray1[k]-plateArray2[k]).len()-meanLen)<closestLen) 
		{
			closestLen=fabs((plateArray1[k]-plateArray2[k]).len()-meanLen);
			meanInd=k;
		}
	}
	projA1.x = plateArray1[largestInd].x*py*fx+projLoc_x;
	projA1.y = plateArray1[largestInd].y*pz*fy+projLoc_y;
	projA2.x = plateArray2[largestInd].x*py*fx+projLoc_x;
	projA2.y = plateArray2[largestInd].y*pz*fy+projLoc_y;
	vertHeightRuler->vertex_list.Add(projA1);
	vertHeightRuler->vertex_list.Add(projA2);
	vertHeightRuler->line_list.Add(IntVec2(count, count+1));
	count+=2;
	height3 = largestLen;

	fprintf(e_fp, "%.3f,%.3f,%.3f,",height1, height2, height3);

	}

	fprintf(e_fp, "\n");
	fclose(e_fp);
	// add 3D line models
//	vertHeightDir->SetDrawMode(LINE_MODE);
//	vertHeightDir->SetEndPoints(cordCenter1, cordCenter2);

	NIH_OpenGL_Refresh_Display(display_id);
	NIH_OpenGL_Refresh_Display(display_id2);
	return CIS_OK;
}



void NIH_BoneSegmentation_Dlg::pushButton_cord_clicked()
{
	if(maskImg3D==NULL || img3D==NULL || cordModel==NULL) return;

	cordModel->pt_list.SetSize(0);
	cordModel->link_list.SetSize(0);
	cordModel->link_prop.SetSize(0);
	cordModel->active_vertex_list.SetSize(0);
	cordModel->plot_radius = 1;

	int sizex, sizey, sizez;
	short *maskA;
	int x, y, z, k, k2;
	float px, py, pz;

	sizex = maskImg3D->Num_Cols();
	sizey = maskImg3D->Num_Rows();
	sizez = maskImg3D->Num_Levels();
	px = img3D->Get_Pixel_SizeX();
	py = img3D->Get_Pixel_SizeY();
	pz = img3D->Get_Pixel_SizeZ();

	maskA = maskImg3D->GetArray();

	IntVec2 center, upper, lower, left, right;
	Vec3 center3D, upper3D, lower3D, left3D, right3D;
	int count;
	vec3DynArray centerList;
	
	for(z=0; z<sizez; z++)
	{
		k= z*sizex*sizey;
		count =0;
		center = IntVec2(0,0);

		// locate the center
		for(k2=0, y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k2++)
			{
				if(maskA[k+k2]==maskStruct.mask_spinalCord)
				{
					count ++;
					center += IntVec2(x, y);
				}
			}
		}

		if(count>0)
		{
			center.x /= count;
			center.y /= count;

			center3D.x = center.x*px;
			center3D.y = center.y*py;
			center3D.z = img3D->GetSlicePosition(z);
			cordModel->pt_list.Add(center3D);

			centerList.Add(center3D);

/*			// locate the four corner of the vertebra
			upper = center;
			for(y=0, k2=0; y<center.y; y++)
			{
				for(x=0; x<sizex; x++, k2++)
				{
					if(maskA[k2+k]==maskStruct.mask_corticalBone || maskA[k2+k]==maskStruct.mask_spongyBone)
					{
						if(y<upper.y)
						{
							upper.y = y;
							upper.x = x;
						}
					}
				}
			}

			lower = center;
			for(y=center.y, k2=y*sizex; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k2++)
				{
					if(maskA[k2+k]==maskStruct.mask_corticalBone || maskA[k2+k]==maskStruct.mask_spongyBone)
					{
						if(y>lower.y)
						{
							lower.y = y;
							lower.x = x;
						}
					}
				}
			}

			right = center;
			for(y=center.y; y<sizey; y++)
			{
				k2 = y*sizex+center.x;
				for(x=center.x; x<sizex; x++, k2++)
				{
					if(maskA[k2+k]==maskStruct.mask_corticalBone || maskA[k2+k]==maskStruct.mask_spongyBone)
					{
						if(x>right.x)
						{
							right.y = y;
							right.x = x;
						}
					}
				}
			}

			left = center;
			for(y=center.y; y<sizey; y++)
			{
				k2 = y*sizex;
				for(x=0; x<center.x; x++, k2++)
				{
					if(maskA[k2+k]==maskStruct.mask_corticalBone || maskA[k2+k]==maskStruct.mask_spongyBone)
					{
						if(x<left.x)
						{
							left.y = y;
							left.x = x;
						}
					}
				}
			}

			upper3D.x = upper.x*px;
			upper3D.y = upper.y*py;
			upper3D.z = center3D.z;
			cordModel->pt_list.Add(upper3D);

			lower3D.x = lower.x*px;
			lower3D.y = lower.y*py;
			lower3D.z = center3D.z;
			cordModel->pt_list.Add(lower3D);

			right3D.x = right.x*px;
			right3D.y = right.y*py;
			right3D.z = center3D.z;
			cordModel->pt_list.Add(right3D);

			left3D.x = left.x*px;
			left3D.y = left.y*py;
			left3D.z = center3D.z;
			cordModel->pt_list.Add(left3D);
*/
		}	// if count
	}	// for z

	// Laplacian method to smooth the center line
	vec3DynArray tmp_centerList;
	tmp_centerList = centerList;

	for(k=1; k<centerList.GetSize()-1; k++)
	{
		centerList[k].x = (tmp_centerList[k-1].x+tmp_centerList[k].x+tmp_centerList[k+1].x)/3;
		centerList[k].y = (tmp_centerList[k-1].y+tmp_centerList[k].y+tmp_centerList[k+1].y)/3;
		
	}

	// smooth the centerline by fitting a b-spline curve
	//
	smoothedCord->UpdateControlList(centerList);

	// add link
/*	for(k=0; k<cordModel->pt_list.GetSize()-1; k++)
	{
		cordModel->link_list.Add(IntVec2(k, k+1));
		cordModel->link_prop.Add(0);
	}
*/
	
	// add the link
	// between plane
/*	for(k=0; k<cordModel->pt_list.GetSize()/5-1; k++)
	{
		cordModel->link_list.Add(IntVec2(k*5, k*5+5));
		cordModel->link_prop.Add(0);
	}

	
	// in plane
	for(k=0; k<cordModel->pt_list.GetSize()/5; k++)
	{
		cordModel->link_list.Add(IntVec2(k*5, k*5+1));
		cordModel->link_prop.Add(1);
		cordModel->link_list.Add(IntVec2(k*5, k*5+2));
		cordModel->link_prop.Add(1);
		cordModel->link_list.Add(IntVec2(k*5, k*5+3));
		cordModel->link_prop.Add(2);
		cordModel->link_list.Add(IntVec2(k*5, k*5+4));
		cordModel->link_prop.Add(2);

		cordModel->active_vertex_list.Add(k*5);
	}
*/
	NIH_OpenGL_CentralizeModel(display_id2, cordModel);
	NIH_OpenGL_Refresh_Display(display_id2);

	short gray;
	int *filledZ=NULL;


	// display the curved reformation
	if(viewMode==0)
	{
		spineCenter.SetSize(sizez);
		spineNormalX.SetSize(sizez);
		spineWidthUp.SetSize(sizez);
		spineWidthDown.SetSize(sizez);

		for(z=0; z<sizez; z++)
		{
			spineCenter[z] = Vec3(0,0,0);
			spineNormalX[z] = Vec3(0,0,0);
			spineWidthUp[z] = -1;
			spineWidthDown[z] = -1;
		}

		ChangeImageSlice_x(curSlice_x);

		// reformate img2D_x and cmapModel1

		filledZ = new int[sizez];
		for(z=0; z<sizez; z++) filledZ[z]=-1;
		(*img2D_x) = short(0);
		for(k=0; k<centerList.GetSize(); k++)
		{
			x = (int)(centerList[k].x/px);
			z = img3D->GetSliceLevel(centerList[k].z,1);

			// filledZ array is used to handle duplicate slices
			//
			if(filledZ[z]==1 && filledZ[z-1]==0) z=z-1;
			while(z<sizez-1 && filledZ[z]==1) z++;
			filledZ[z] = k;
		
			for(y=0; y<img2D_x->Num_Cols(); y++)
			{
				gray = img3D->FastGet(x,y,z);
				img2D_x->FastSet(y, z, gray);
			}
		}
		mapModel1->ChangeBitmap(img2D_x, mapModel1->mapping_mode);

	
		// draw the projected curve
		vec2DynArray projA, leftA, rightA, normalA;
		intDynArray leftWidth, rightWidth;
		intDynArray centerX, centerY;

		centerX.SetSize(sizez);
		centerY.SetSize(sizez);
		normalA.SetSize(sizez);

		float vsizey = img3D->Get_Volume_SizeY();
		float vsizez = img3D->Get_Volume_SizeZ();
		float fx, fy;
		int sy;
		short mask;

		fx = mapModel1->size.x/vsizey;
		fy = mapModel1->size.y/vsizez;

		projA.SetSize(centerList.GetSize());
		leftA.SetSize(centerList.GetSize());
		rightA.SetSize(centerList.GetSize());

		leftWidth.SetSize(sizez);
		rightWidth.SetSize(sizez);
		
		for(z=0; z<sizez; z++) 
		{
			filledZ[z]=-1;
			rightWidth[z] = leftWidth[z]=0;
			centerX[z] = 0;
			centerY[z] = 0;
			normalA[z] = Vec2(0,0);
		}
		
		for(k=0; k<projA.GetSize(); k++)
		{
			x = (int)(centerList[k].x/px);
			y = (int)(centerList[k].y/py);
			z = img3D->GetSliceLevel(centerList[k].z,1);

			if(filledZ[z]==1 && filledZ[z-1]==-1) z=z-1;
			while(z<sizez-1 && filledZ[z]!=-1) z++;
			filledZ[z] = k;

			spineCenter[z].x = centerList[k].x;
			spineCenter[z].y = centerList[k].y;
			spineCenter[z].z = img3D->GetSlicePosition(z);

			projA[k].x = centerList[k].y*fx+mapModel1->loc.x;
			projA[k].y = z*pz*fy+mapModel1->loc.y;

			centerX[z] = x;
			centerY[z] = y;

			// set up the spinal column
			// to the right (down)
			for(sy=y; sy<sizey; sy++)
			{
				mask = maskImg3D->FastGet(x, sy, z);
				if(mask!=maskStruct.mask_spongyBone && mask!=maskStruct.mask_corticalBone && mask!=maskStruct.mask_spinalCord)
				{
					rightWidth[z] = sy-y;
					break;
				}
			}

			// to the left (up)
			for(sy=y; sy>0; sy--)
			{
				mask = maskImg3D->FastGet(x, sy, z);
				if(mask==maskStruct.mask_body)
				{
					leftWidth[z] = y-sy;
					break;
				}
			}

		}	// for k
	
		int topSlice, bottomSlice;
		topSlice = sizez;
		bottomSlice = 0;
		for(z=0;z<sizez;z++)
		{
			if(leftWidth[z]>0)
			{
				if(z<topSlice) topSlice=z;
				if(z>bottomSlice) bottomSlice=z;
			}
		}

		int avgLeftWidth, avgRightWidth, avgLeftWidth1, avgLeftWidth2;
		int count, count1, count2;;
		avgLeftWidth = avgRightWidth = 0;
		avgLeftWidth1 = avgLeftWidth2 = 0;
		count = 0;
		count1 = count2=0;
		for(z=0; z<sizez; z++)
		{
			if(leftWidth[z]>0 && rightWidth[z]>0)
			{
				avgLeftWidth += leftWidth[z];
				avgRightWidth += rightWidth[z]; 
				count ++;

				// get the average of first half and second half
				if(z<(topSlice+bottomSlice)/2) 
				{
					avgLeftWidth1 += leftWidth[z];
					count1++;
				}
				else
				{
					avgLeftWidth2 += leftWidth[z];
					count2++;
				}
			}
		}
	 	
		if(count!=0)
		{
			avgLeftWidth /= count;
			avgRightWidth /= count;

//			avgLeftWidth += 4;
//			avgRightWidth += 4;
		}

		if(count1!=0)
		{
			avgLeftWidth1 /= count1;
		}
		
		if(count2!=0)
		{
			avgLeftWidth2 /= count2;
		}
		
		// smooth out outlier
		for(z=0; z<sizez; z++)
		{
			if(leftWidth[z]>0 && rightWidth[z]>0)
			{
				if(z<=(topSlice+bottomSlice)/2 && leftWidth[z]<(float)avgLeftWidth1*0.75)
				{
					leftWidth[z] = avgLeftWidth1;
				}
				if(z>=(topSlice+bottomSlice)/2 && leftWidth[z]<(float)avgLeftWidth2*0.75)
				{
					leftWidth[z] = avgLeftWidth2;
				}
				if(z>=(topSlice+bottomSlice)/2 && leftWidth[z]<avgLeftWidth1)
				{
					leftWidth[z] = avgLeftWidth2;
				}
			}
		}

		// smooth leftWidth and rightWidth
		intDynArray tmpWidth;
		tmpWidth = leftWidth;
		for(z=1; z<sizez-1; z++)
		{
			if(tmpWidth[z-1]==0 || tmpWidth[z+1]==0) continue;
			if(tmpWidth[z]<tmpWidth[z-1] && tmpWidth[z]<tmpWidth[z+1]) leftWidth[z] = (tmpWidth[z-1]+tmpWidth[z+1])/2;
			else leftWidth[z] = (tmpWidth[z-1]+tmpWidth[z]+tmpWidth[z+1])/3;
		}
		tmpWidth = rightWidth;
		for(z=1; z<sizez-1; z++)
		{
			if(tmpWidth[z-1]==0 || tmpWidth[z+1]==0) continue;
			rightWidth[z] = (tmpWidth[z-1]+tmpWidth[z]+tmpWidth[z+1])/3;
		}

		// assign globle variable
		for(z=0; z<sizez; z++)
		{
			spineWidthUp[z] = leftWidth[z]*py;
			spineWidthDown[z] = rightWidth[z]*py;
		}

		projectedCord->UpdateControlList(projA);
		projectedCord->SetDrawMode(LINE_MODE);

		// compute the normal from the b-spline curve
		double curveLength=projectedCord->Length();
		double stepSize=curveLength/projA.GetSize();
		double s=0;
		for(z=0;z<sizez;z++)
		{
			k = filledZ[z];
			if(k<0) continue;
			normalA[z] = projectedCord->Normal(s);
			s += stepSize;
		}

		// smooth the normal using Laplacian method
		vec2DynArray tmp_normalA;
		tmp_normalA = normalA;
		for(z=1; z<sizez-1; z++)
		{
			if(tmp_normalA[z-1].x==0 || tmp_normalA[z+1].x==0) continue;
			normalA[z] = (tmp_normalA[z-1]+tmp_normalA[z]+tmp_normalA[z+1]).normalize();
		}


		// assign globle variable
		for(z=0; z<sizez; z++)
		{
			spineNormalX[z] = Vec3(0, normalA[z].x, normalA[z].y);
		}

		// define left and right column
		z = 0;
		for(k=0; k<projA.GetSize(); k++)
		{
			for(;z<sizez;z++) if(filledZ[z]==k) break;
			leftA[k].x = projA[k].x-leftWidth[z]*px*fx;
			leftA[k].y = projA[k].y;
			rightA[k].x = projA[k].x+rightWidth[z]*px*fx;
			rightA[k].y = projA[k].y;
		}
		
		leftColumn->UpdateControlList(leftA);
		leftColumn->SetDrawMode(LINE_MODE);
		rightColumn->UpdateControlList(rightA);
		rightColumn->SetDrawMode(HIDE_MODE);

	
		// set up the pedicle model
		pedicleModel->vertex_list.SetSize(0);
		pedicleModel->line_list.SetSize(0);
		pedicleModel->SetDrawMode(LINE_MODE);
		pedicleModel->plot_radius = 0;
		count=0;
		// pedicleValley was computed from Coronal view
		for(k=0; k<pedicleValley.GetSize(); k++)
//		for(z=0; z<sizez; z++)
		{
			z = pedicleValley[k];
			if(filledZ[z]==-1) continue;

			pedicleModel->vertex_list.Add(projA[filledZ[z]]);

			pedicleModel->vertex_list.Add(projA[filledZ[z]]-normalA[z]*leftWidth[z]*px*fx);
//			pedicleModel->vertex_list.Add(projA[filledZ[z]]+normalA[z]*rightWidth[z]*px*fx);
			pedicleModel->line_list.Add(IntVec2(count, count+1));
//			pedicleModel->line_list.Add(IntVec2(count, count+2));
			count+=2;
		}
	
	}
	else if(viewMode==1)
	{
		spineCenter.SetSize(sizez);
		spineNormalY.SetSize(sizez);
		spineWidthLeft.SetSize(sizez);
		spineWidthRight.SetSize(sizez);

		for(z=0; z<sizez; z++)
		{
			spineCenter[z] = Vec3(0,0,0);
			spineNormalY[z] = Vec3(0,0,0);
			spineWidthLeft[z] = -1;
			spineWidthRight[z] = -1;
		}

		// display the curved reformation
		ChangeImageSlice_y(curSlice_y);

		filledZ = new int[sizez];

		for(z=0; z<sizez; z++) filledZ[z]=-1;
		
		(*img2D_y) = short(0);
		for(k=0; k<centerList.GetSize(); k++)
		{
			y = (int)(centerList[k].y/py);
			z = img3D->GetSliceLevel(centerList[k].z,1);

			if(filledZ[z]==1 && filledZ[z-1]==-1) z=z-1;
			while(z<sizez-1 && filledZ[z]!=-1) z++;
			filledZ[z] = k;
		
			for(x=0; x<img2D_y->Num_Cols(); x++)
			{
				gray = img3D->FastGet(x,y,z);
				img2D_y->FastSet(x, z, gray);
			}
		}
		mapModel1->ChangeBitmap(img2D_y, mapModel1->mapping_mode);

		// draw the projected curve
		vec2DynArray projA, leftA, rightA, normalA;
		intDynArray leftWidth, rightWidth;
		intDynArray centerX, centerY;

		centerX.SetSize(sizez);
		centerY.SetSize(sizez);
		normalA.SetSize(sizez);

		float vsizex = img3D->Get_Volume_SizeX();
		float vsizez = img3D->Get_Volume_SizeZ();
		float fx, fy;
		int sx;
		short mask;

		fx = mapModel1->size.x/vsizex;
		fy = mapModel1->size.y/vsizez;

		projA.SetSize(centerList.GetSize());
		leftA.SetSize(centerList.GetSize());
		rightA.SetSize(centerList.GetSize());

		leftWidth.SetSize(sizez);
		rightWidth.SetSize(sizez);
		
		for(z=0; z<sizez; z++) 
		{
			filledZ[z]=-1;
			rightWidth[z] = leftWidth[z]=0;
			centerX[z] = 0;
			centerY[z] = 0;
			normalA[z] = Vec2(0,0);
		}
		
		for(k=0; k<projA.GetSize(); k++)
		{
			x = (int)(centerList[k].x/px);
			y = (int)(centerList[k].y/py);
			z = img3D->GetSliceLevel(centerList[k].z,1);

			if(filledZ[z]==1 && filledZ[z-1]==-1) z=z-1;
			while(z<sizez-1 && filledZ[z]!=-1) z++;
			filledZ[z] = k;

			projA[k].x = centerList[k].x*fx+mapModel1->loc.x;
			projA[k].y = z*pz*fy+mapModel1->loc.y;


			spineCenter[z].x = centerList[k].x;
			spineCenter[z].y = centerList[k].y;
			spineCenter[z].z = img3D->GetSlicePosition(z);

			centerX[z] = x;
			centerY[z] = y;

			// set up the spinal column
			// to the right
			for(sx=x; sx<sizex; sx++)
			{
				mask = maskImg3D->FastGet(sx, y, z);
				if(mask!=maskStruct.mask_spongyBone && mask!=maskStruct.mask_corticalBone && mask!=maskStruct.mask_spinalCord)
				{
					rightWidth[z] = sx-x;
					break;
				}
			}
			rightA[k].x = projA[k].x+rightWidth[z]*px*fx;
			rightA[k].y = projA[k].y;

			// to the left
			for(sx=x; sx>0; sx--)
			{
				mask = maskImg3D->FastGet(sx, y, z);
				if(mask!=maskStruct.mask_spongyBone && mask!=maskStruct.mask_corticalBone && mask!=maskStruct.mask_spinalCord)
				{
					leftWidth[z] = x-sx;
					break;
				}
			}
			leftA[k].x = projA[k].x-leftWidth[z]*px*fx;
			leftA[k].y = projA[k].y;

		}	// for k
		
		int avgLeftWidth, avgRightWidth;
		int count;
		avgLeftWidth = avgRightWidth = 0;
		count = 0;
		for(z=0; z<sizez; z++)
		{
			if(leftWidth[z]>0 && rightWidth[z]>0)
			{
				avgLeftWidth += leftWidth[z];
				avgRightWidth += rightWidth[z]; 
				count ++;
			}
		}
	 	
		if(count!=0)
		{
			avgLeftWidth /= count;
			avgRightWidth /= count;

			avgLeftWidth += 4;
			avgRightWidth += 4;
		}
		
		// assign globle variable
		for(z=0; z<sizez; z++)
		{
			spineWidthLeft[z] = avgLeftWidth*px;
			spineWidthRight[z] = avgRightWidth*px;
		}

		projectedCord->UpdateControlList(projA);
		projectedCord->SetDrawMode(LINE_MODE);

		// compute the normal from the b-spline curve
		double curveLength=projectedCord->Length();
		double stepSize=curveLength/projA.GetSize();
		double s=0;
		for(z=0;z<sizez;z++)
		{
			k = filledZ[z];
			if(k<0) continue;
			normalA[z] = projectedCord->Normal(s);
			s += stepSize;
		}

		// smooth the normal using Laplacian method
		vec2DynArray tmp_normalA;
		tmp_normalA = normalA;
		for(z=1; z<sizez-1; z++)
		{
			if(tmp_normalA[z-1].x==0 || tmp_normalA[z+1].x==0) continue;
			normalA[z] = (tmp_normalA[z-1]+tmp_normalA[z]+tmp_normalA[z+1]).normalize();
		}

		// assign globle variable
		for(z=0; z<sizez; z++)
		{
			spineNormalY[z] = Vec3(normalA[z].x, 0, normalA[z].y);
		}

		// define left and right column
		for(k=0; k<projA.GetSize(); k++)
		{
			leftA[k].x = projA[k].x-avgLeftWidth*px*fx;
			leftA[k].y = projA[k].y;
			rightA[k].x = projA[k].x+avgRightWidth*px*fx;
			rightA[k].y = projA[k].y;
		}
		
		leftColumn->UpdateControlList(leftA);
		leftColumn->SetDrawMode(LINE_MODE);
		rightColumn->UpdateControlList(rightA);
		rightColumn->SetDrawMode(LINE_MODE);


		// compute the intensity integral along the normal direction
		intDynArray leftAvgI, rightAvgI, leftCount, rightCount, leftAvgI1, rightAvgI1, leftAvgI2, rightAvgI2, allCount, allAvgI;
		int count0, count1, count2;
		leftAvgI.SetSize(sizez);
		leftAvgI1.SetSize(sizez);
		leftAvgI2.SetSize(sizez);
		rightAvgI.SetSize(sizez);
		rightAvgI1.SetSize(sizez);
		rightAvgI2.SetSize(sizez);
		leftCount.SetSize(sizez);
		rightCount.SetSize(sizez);
		allCount.SetSize(sizez);
		allAvgI.SetSize(sizez);

		Vec2 startPos, curPos;
		for(z=0; z<sizez; z++)
		{
			leftAvgI[z] = leftAvgI1[z] = leftAvgI2[z] = 0;
			rightAvgI[z] = rightAvgI1[z] = rightAvgI2[z] = 0;
			leftCount[z] = rightCount[z] = 0;
			allCount[z]=0;
			allAvgI[z]=0;
			if(centerX[z]==0) continue;

			// left side
			count0 = count1 = count2 = 0;
			startPos = Vec2(centerX[z], z);
			for(k=0; k<avgLeftWidth; k++)
			{
				curPos.x = startPos.x-normalA[z].x*k;
				curPos.y = startPos.y-normalA[z].y*k*px/pz;

				mask = maskImg3D->FastGet((int)curPos.x, centerY[z], (int)curPos.y);
				gray = img3D->FastGet((int)curPos.x, centerY[z], (int)curPos.y);

				if(mask==maskStruct.mask_corticalBone || mask==maskStruct.mask_spongyBone)
				{
					count0++;
					leftAvgI[z]+=gray;
				}
				
				if(mask!=maskStruct.mask_spinalCord)
				{
					count1++;
					leftAvgI1[z]+=gray;
					allCount[z]++;
					allAvgI[z]+=gray;
				}

				count2++;
				leftAvgI2[z] += gray;
			}

			leftCount[z] = count1;
			if(count0!=0) leftAvgI[z]/= count0;
			if(count1!=0) leftAvgI1[z]/= count1;
			if(count2!=0) leftAvgI2[z]/= count2;

			// right side
			count0 = count1 = count2 = 0;
			startPos = Vec2(centerX[z], z);
			for(k=0; k<avgRightWidth; k++)
			{
//				curPos = startPos+normalA[z]*k;
				curPos.x = startPos.x-normalA[z].x*k;
				curPos.y = startPos.y-normalA[z].y*k*px/pz;

//				mask = maskImg3D->FastGet(x, centerY[z], z);
//				gray = img3D->FastGet(x, centerY[z], z);

				mask = maskImg3D->FastGet((int)curPos.x, centerY[z], (int)curPos.y);
				gray = img3D->FastGet((int)curPos.x, centerY[z], (int)curPos.y);

				if(mask==maskStruct.mask_corticalBone || mask==maskStruct.mask_spongyBone)
				{
					count0++;
					rightAvgI[z]+=gray;
				}
				
				if(mask!=maskStruct.mask_spinalCord)
				{
					count1++;
					rightAvgI1[z]+=gray;
					allCount[z]++;
					allAvgI[z]+=gray;
				}

				count2++;
				rightAvgI2[z] += gray;
			}

			rightCount[z] = count1;
			if(count0!=0) rightAvgI[z]/= count0;
			if(count1!=0) rightAvgI1[z]/= count1;
			if(count2!=0) rightAvgI2[z]/= count2;
			
			if(allCount[z]!=0) allAvgI[z] /= allCount[z];
		}	// for z

		// export to a text file
		FILE *fp;
		fp = fopen("c:\\tmp\\integral.txt", "w");
		fprintf(fp, "leftCount leftAvg leftAvg1 leftAvg2 rightCount rightAvg rightAvg1 rightAvg2 allCount allAvg\n");
		for(z=0;z<sizez; z++)
		{
			fprintf(fp, "%d %d %d %d %d %d %d %d %d %d\n", 
				leftCount[z], leftAvgI[z], leftAvgI1[z], leftAvgI2[z], 
				rightCount[z], rightAvgI[z], rightAvgI1[z], rightAvgI2[z],
				allCount[z], allAvgI[z]);
		}

		// try to locate the pedicle
		pedicleValley.SetSize(0);
		for(z=3; z<sizez-3; z++)
		{
			if(allAvgI[z]<allAvgI[z-1] && allAvgI[z]<allAvgI[z-2] && allAvgI[z]<=allAvgI[z+1] && allAvgI[z]<=allAvgI[z+2])
			{
				pedicleValley.Add(z);
			}
		}

		double avgInterval;
		avgInterval = 0;
		for(k=1; k<pedicleValley.GetSize(); k++)
		{
			avgInterval += pedicleValley[k]-pedicleValley[k-1];
		}
		avgInterval /= (double)(pedicleValley.GetSize()-1);

		// insert a valley if the interval is too big
		for(k=1; k<pedicleValley.GetSize(); k++)
		{
			if(pedicleValley[k]-pedicleValley[k-1] >= avgInterval*1.75)
			{
				int tz=(pedicleValley[k]+pedicleValley[k-1])/2;
				pedicleValley.InsertAt(k,tz);
			}
		}

		// set up the pedicle model
		pedicleModel->vertex_list.SetSize(0);
		pedicleModel->line_list.SetSize(0);
		pedicleModel->SetDrawMode(LINE_MODE);
		pedicleModel->plot_radius = 0;
		count=0;
		for(k=0; k<pedicleValley.GetSize(); k++)
		{
			z = pedicleValley[k];
			if(filledZ[z]==-1) continue;
			// report pedicle direction and normal
			fprintf(fp, "%d %.2f %.2f \n", pedicleValley[k], normalA[z].x, normalA[z].y);

			pedicleModel->vertex_list.Add(projA[filledZ[z]]);

			pedicleModel->vertex_list.Add(projA[filledZ[z]]-normalA[z]*avgLeftWidth*px*fx);
			pedicleModel->vertex_list.Add(projA[filledZ[z]]+normalA[z]*avgRightWidth*px*fx);
			pedicleModel->line_list.Add(IntVec2(count, count+1));
			pedicleModel->line_list.Add(IntVec2(count, count+2));
			count+=3;
		}

		fclose(fp);
	}	// else viewMode
	else if(viewMode==2)
	{
		// display the curved reformation
		ChangeImageSlice_y(curSlice_y);

		filledZ = new int[sizez];
		for(z=0; z<sizez; z++) filledZ[z]=0;
		(*img2D_y) = short(0);
		int integrateI, sy, mask;

		for(k=0; k<centerList.GetSize(); k++)
		{
			y = (int)(centerList[k].y/py);
			z = img3D->GetSliceLevel(centerList[k].z,1);

			if(filledZ[z]==1 && filledZ[z-1]==0) z=z-1;
			while(z<sizez-1 && filledZ[z]==1) z++;
			filledZ[z] = 1;
		
			for(x=0; x<img2D_y->Num_Cols(); x++)
			{
				integrateI = 0;
				for(sy=0; sy<y; sy++)
				{
					mask = maskImg3D->FastGet(x, sy, z);
					if(mask==maskStruct.mask_corticalBone || mask==maskStruct.mask_spongyBone)
					{
						integrateI += img3D->FastGet(x, sy, z);
					}
				}
				gray = integrateI/10;
				img2D_y->FastSet(x, z, gray);
			}
		}
		mapModel1->ChangeBitmap(img2D_y, mapModel1->mapping_mode);
	}
	
	
	NIH_OpenGL_Refresh_Display(display_id);
	projectedCord->SetDrawMode(HIDE_MODE);
	leftColumn->SetDrawMode(HIDE_MODE);
	rightColumn->SetDrawMode(HIDE_MODE);
	pedicleModel->SetDrawMode(HIDE_MODE);
	if(filledZ) delete filledZ;
}

void NIH_BoneSegmentation_Dlg::pushButton_watershed_clicked()
{
	if(img2D==NULL || img3D==NULL) return;
	if(segInfo.spineBound1.GetSize()!=img3D->Num_Levels()) return;

	QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );

	QTime qtime;

	qtime.start();

	bool debugMode_1;

	if(flagShowPainting && !flagBatchMode) debugMode_1=true;
	else debugMode_1 = false;

	int bbx0, bby0, bbx1, bby1, z;

	int stz, edz;

	//debugMode_1 =true;
	if(debugMode_1)
	{
		stz = edz = curSlice;
	}
	else
	{
		stz = 0; edz = img3D->Num_Levels()-1;
	}
	
	int sizex = maskImg3D->Num_Cols();
	int sizey = maskImg3D->Num_Rows();
	int sizez = maskImg3D->Num_Levels();
	int sizexy=sizex*sizey;
	short *maskA=maskImg3D->GetArray();
	short *imgA=img3D->GetArray();

	if(debugMode_1)
	{
		int x,y, z,k;
		
		for(z=0, k=0; z<sizez; z++)
		{
			for(y=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k++)
				{
					if(maskA[k]==maskStruct.mask_boneMetasis ||
						maskA[k]==maskStruct.mask_falseDetection) 
					{
						maskA[k]=maskStruct.mask_spongyBone;
					}
				}
			}
		}
	}

	for(z=stz; z<=edz; z++)
	{
		if(segInfo.spineBound1[z].x==-1) continue;

		ChangeImageSlice(z);

		// get a subimage around the spine on current slice
		bbx0 = segInfo.spineBound1[z].x-20;
		bby0 = segInfo.spineBound1[z].y-5;
		bbx1 = segInfo.spineBound2[z].x+20;
		bby1 = segInfo.spineBound2[z].y+5;

		if(bbx0>=bbx1 || bby0>=bby1) continue;
		if(bbx0<0 || bbx0>sizex-1 || bbx1<0 || bbx1>sizex-1 || bby0<0 || bby0>sizey-1 || bby1<0 || bby1>sizey-1) continue;

		// get subimage within the bounding box
		CIS_Array_Image2D_short *subImg, *waterLabel, *smoothedData;
		int subSizex, subSizey, subSizexy;
	
		
		subImg = img2D->SubImage(bbx0, bby0, bbx1-bbx0+1, bby1-bby0+1);
		subSizex = subImg->Num_Cols();
		subSizey = subImg->Num_Rows();
		subSizexy = subSizex*subSizey;

		waterLabel = new CIS_Array_Image2D_short(subSizex, subSizey);
		smoothedData = new CIS_Array_Image2D_short(subSizex, subSizey);

		CIS_Array_Image2D_short *binImg2Dcort;
		CIS_Array_Image2D_short *binImg2Dedge;
		short *binArray2Dcort;
		short *binArray2Dedge;

		binImg2Dcort = new CIS_Array_Image2D_short(sizex, sizey);
		binImg2Dedge = new CIS_Array_Image2D_short(sizex, sizey);
		binArray2Dcort = binImg2Dcort->GetArray();
		binArray2Dedge = binImg2Dedge->GetArray();

		int xx,yy,kk,kk2;
		for(kk2=0; kk2<sizexy; kk2++)
			binArray2Dcort[kk2] = binArray2Dedge[kk2] = 0;
		kk = z*sizexy;
		// image of cortical and spongy bone
		for(yy=0; yy<sizey; yy++)
		{
			kk2 = yy*sizex;
			for(xx=0; xx<sizex; xx++, kk2++)
			{
				if(maskA[kk+kk2]==maskStruct.mask_corticalBone || maskA[kk+kk2]==maskStruct.mask_spongyBone)
				{
					binArray2Dcort[kk2] = binArray2Dedge[kk2] = 1;
				}
			} //for x
		} //for y

		/*FILE * pFile;
		pFile = fopen("c:\\tmp\\watershedz.txt","w");
		//int xx, yy;
		for (yy=0; yy<sizey; yy++)
		{
			for (xx=0; xx<sizex; xx++)
			{
				fprintf (pFile, "%d", binImg2Dcort->FastGet(xx,yy));
			}
			fprintf( pFile, "\n");
		}
		fclose (pFile);*/
		
		CIS_IPA_Erode(binImg2Dcort, 5, false);

		//set mask to the eroded image
		kk = z*sizexy;
		for(yy=0; yy<sizey; yy++)
		{
			kk2 = yy*sizex;
			for(xx=0; xx<sizex; xx++, kk2++)
			{
				if(binArray2Dcort[kk2]==0 && (maskA[kk+kk2]==maskStruct.mask_corticalBone || maskA[kk+kk2]==maskStruct.mask_spongyBone))
				{
					maskA[kk+kk2]=0;
				}
				else
				{
					binArray2Dedge[kk2] = 0;
				}
			} //for x
		} //for y
		kk=z*sizexy;
		for(yy=0; yy<sizey; yy++)
		{
			kk2 = yy*sizex;
			for(xx=0; xx<sizex; xx++, kk2++)
			{
				if(binArray2Dedge[kk2]==1)
				{
					maskA[kk+kk2]=48;
				}
			} //for x
		}

		// Preprocess
		SpineSegmentation_PreProcess_Watershed(maskImg3D, subImg, bbx0, bby0, curSlice, 
			maskStruct, segPara, segInfo);

		// watershed on the subimage
		SpineSegmentation_Watershed_usingITK(subImg, smoothedData, waterLabel, debugMode);

		watershedOverlay = new CIS_Array_Image2D_RGB_uchar(sizex, sizey);

		// post processing the water shed result
		SpineSegmentation_PostProcess_Watershed(maskImg3D, subImg, smoothedData,  waterLabel, binArray2Dedge, watershedOverlay, bbx0, bby0, curSlice, 
			maskStruct, segPara, segInfo, debugMode);

		//to visualize the vertebra border
		int x, y, k;
		int sizexy=sizex*sizey;
		for(y=0; y<subSizey; y++)
		{
			for(x=0; x<subSizex; x++)
			{
				k =z*sizexy+(y+bby0)*sizex+(x+bbx0);
				if(maskA[k]==maskStruct.mask_spongyBone || maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_spinalCord
					|| imgA[k]>1300) subImg->FastSet(x,y,1);
				else subImg->FastSet(x,y,0);
			}
		}
		// only keep the largest blob
		intDynArray blobRank2D;
		intVec2DynArray blobCentroid2D;
		NIH_Algo_Blob_Labelling(subImg, smoothedData, blobRank2D, blobCentroid2D, true);
		int largestBlob=-1, largestSize=0;
		for(k=0; k<blobRank2D.GetSize(); k++)
		{
			if(blobRank2D[k]>largestSize)
			{
				largestSize = blobRank2D[k];
				largestBlob=k+1;
			}
		}
		for(y=0; y<subSizey; y++)
		{
			for(x=0; x<subSizex; x++)
			{
				if(smoothedData->FastGet(x,y)!=largestBlob) subImg->FastSet(x,y,0);
			}
		}

		CIS_IPA_Dilate(subImg, 2, true);
		CIS_IPA_Erode(subImg, 2, true);

		// fill holes

		for(y=0; y<subSizey; y++)
		{
			for(x=0; x<subSizex; x++)
			{
				waterLabel->FastSet(x,y, 1-subImg->FastGet(x,y));
			}
		}
		NIH_Algo_Blob_Labelling(waterLabel, smoothedData, blobRank2D, blobCentroid2D, false);
		for(y=0; y<subSizey; y++)
		{
			for(x=0; x<subSizex; x++)
			{
				if(smoothedData->FastGet(x,y)!=1) subImg->FastSet(x,y,1);
			}
		}

	/*	seedModel->vertex_list.SetSize(0);
		seedModel->draw_mode = VERTEX_MODE;
		seedModel->SetCurrColor(Vec3(0,0,1));
		// Any element of the largest blob that is adjacent to a nonelement is considered to be on the border
		for(y=1; y<subSizey-1; y++)
		{
			for(x=1; x<subSizex-1; x++)
			{
				if(subImg->FastGet(x,y)==1 && 
				(subImg->FastGet(x-1,y)==0 || subImg->FastGet(x+1,y)==0 ||subImg->FastGet(x,y-1)==0 ||subImg->FastGet(x,y+1)==0)) 
				seedModel->vertex_list.Add(Vec2(x+bbx0, y+bby0));
			}
		}*/

		delete subImg;
		delete smoothedData;
		delete waterLabel;
		delete binImg2Dcort;
		delete binImg2Dedge;
	} // for z

	// do the batchMode FROC 
	if((flagApplyClassifier || flagApply3DClassifier) && flagBatchMode && svmCutoffStep>0)
	{
		spinBox_cutoff->setRange((int) svmCutoffMin/svmCutoffStep, (int) svmCutoffMax/svmCutoffStep);
		spinBox_cutoff->setValue((int) svmCutoffMin/svmCutoffStep);
		int index = 0;

		for(float curCutoff=svmCutoffMin; curCutoff<=svmCutoffMax; curCutoff+=svmCutoffStep)
		{
			//flagApply3DClassifier = true;
			
			if(flagApply3DClassifier)
			{
				NIH_BoneSegmentation_AnalyzeMetasis3D(img3D, maskImg3D, maskStruct, segPara, segInfo, detections3D, numDetections3D, detections, numDetections, detectionVoxels,
					lesions, numLesions, paint, flagApply3DClassifier, svmCommittee, 
					curCutoff, studyEntries[curStudy].sliceRange1, studyEntries[curStudy].sliceRange2);
			}
			else
			{
				NIH_BoneSegmentation_AnalyzeMetasis2D(img3D, maskImg3D, maskStruct, segPara, segInfo, detections2D, numDetections2D, 
					lesions, numLesions, paint, flagApplyClassifier, flagApply3DClassifier, svmCommittee, 
					curCutoff, studyEntries[curStudy].sliceRange1, studyEntries[curStudy].sliceRange2);
				NIH_BoneSegmentation_MergeMetasis(img3D, maskImg3D, maskStruct, segPara, segInfo, detections, numDetections, detectionVoxels);
			}

	
			int totalMet, totalLytic, foundMet, foundLytic, totalDetection, tpDetection, fpDetection, 
				tpLyticDetection, fpLyticDetection;
			int totalMet10, totalLytic10, foundMet10, foundLytic10;

			NIH_BoneSegmentation_MatchDetections(maskImg3D, detections, numDetections, detectionVoxels, lesions, numLesions, paint, 
				totalMet, totalLytic, foundMet, foundLytic, totalDetection, tpDetection, fpDetection,
				tpLyticDetection, fpLyticDetection, totalMet10, totalLytic10, foundMet10, foundLytic10, lesionSizeCutoff);

			

			// reporting
			float sensitivity = 1, sensitivityLytic = 1;
			if(totalMet>0) 
			{
				sensitivity = (float)foundMet/(float)totalMet;
			}
			if(totalLytic>0) 
			{
				sensitivityLytic = (float)foundLytic/(float)totalLytic;
			}
			
			// report
			FILE *fp;
			if((fp = fopen(globalReportFileName, "a"))!=NULL)
			{
			//	fprintf(fp,"%s,%d,%d,%.3f,%d,%d,%d,%d,%d,%.3f,%d,%d,%d,%d,%d,%d,%d,%.2f\n",
				fprintf(fp,"%s,%d,%d,%.3f,%d,%d,%d,%.2f\n",
					studyEntries[curStudy].patientName.ascii(),
					totalMet,foundMet,sensitivity,totalDetection, tpDetection, fpDetection,
					/*totalLytic,foundLytic,sensitivityLytic,totalDetection, tpLyticDetection, fpLyticDetection, 
					totalMet10, foundMet10, totalLytic10, foundLytic10,*/
					curCutoff);
				fclose(fp);
			}

			// store the info

			//storeMasks[index].mask3DArray=maskImg3D->GetArray();
			storeMasks[index].mask3DArray.reserve(sizexy*sizez);
			for(int i=0; i<sizexy*sizez; i++)
			{
				storeMasks[index].mask3DArray.push_back(maskA[i]);
			}

			storeMasks[index].curCutoff = curCutoff;
			storeMasks[index].foundMet = foundMet;
			storeMasks[index].fpDetection = fpDetection;
			storeMasks[index].sensitivity = sensitivity;
			storeMasks[index].totalDetection = totalDetection;
			storeMasks[index].totalMet = totalMet;
			storeMasks[index].tpDetection = tpDetection;
			for(int i=0; i<numLesions; i++)
			{
				storeMasks[index].lesionStatus[i] = lesions[i].lesionStatus;
			}
			for(int i=0; i<numDetections; i++)
			{
				storeMasks[index].matchedLesion[i] = detections[i].matchedLesion;
			}
			index++;
			std::cout<<"\n"<<storeMasks[0].mask3DArray[138008];
		}
	}
	else
	{
			/*NIH_BoneSegmentation_AnalyzeMetasis2D(img3D, maskImg3D, maskStruct, segPara, segInfo, detections2D, numDetections2D, 
			lesions, numLesions, paint, flagApplyClassifier, flagApply3DClassifier, svmCommittee, svmCutoff, 
				studyEntries[curStudy].sliceRange1, studyEntries[curStudy].sliceRange2);*/
		/*if(flagBatchMode)
		{
			NIH_BoneSegmentation_AnalyzeMetasis3D(img3D, maskImg3D, maskStruct, segPara, segInfo, detections3D, numDetections3D, detections, numDetections, detectionVoxels,
					lesions, numLesions, paint, flagApply3DClassifier, svmCommittee, 
					svmCutoff, studyEntries[curStudy].sliceRange1, studyEntries[curStudy].sliceRange2);
		}
		else*/
		{
			NIH_BoneSegmentation_AnalyzeMetasis2D(img3D, maskImg3D, maskStruct, segPara, segInfo, detections2D, numDetections2D, 
			lesions, numLesions, paint, flagApplyClassifier, flagApply3DClassifier, svmCommittee, svmCutoff, -1, -1);
			NIH_BoneSegmentation_MergeMetasis(img3D, maskImg3D, maskStruct, segPara, segInfo, detections, numDetections, detectionVoxels);
		}
	
		spinBox_detection->setRange(0, numDetections-1);
		spinBox_detection->setValue(0);

		spinBox_detection_2D->setRange(0, numDetections2D-1);
		spinBox_detection_2D->setValue(0);

		int totalMet, totalLytic, foundMet, foundLytic, totalDetection, tpDetection, fpDetection, 
			tpLyticDetection, fpLyticDetection;
		int totalMet10, totalLytic10, foundMet10, foundLytic10;

		NIH_BoneSegmentation_MatchDetections(maskImg3D, detections, numDetections, detectionVoxels, lesions, numLesions, paint, 
			totalMet, totalLytic, foundMet, foundLytic, totalDetection, tpDetection, fpDetection,
			tpLyticDetection, fpLyticDetection, totalMet10, totalLytic10, foundMet10, foundLytic10, lesionSizeCutoff);

		// reporting
		float sensitivity = 1, sensitivityLytic = 1;
		if(totalMet>0) 
		{
			sensitivity = (float)foundMet/(float)totalMet;
		}
		if(totalLytic>0) 
		{
			sensitivityLytic = (float)foundLytic/(float)totalLytic;
		}
		
		char info[500];
		sprintf(info, "Total Mets: %d (%d); Found: %d (%d);\nTotal Lytic %d (%d); Found: %d (%d)\nDetection: %d; TP: %d; Lytic: %d; FP: %d\n",
			totalMet, totalMet10, foundMet, foundMet10,
			totalLytic, totalLytic10, foundLytic, foundLytic10, 
			totalDetection, tpDetection, tpLyticDetection, fpDetection);

		textLabel_info->setText(info);

		// report
		FILE *fp;
		if((fp = fopen(globalReportFileName, "a"))!=NULL)
		{
			fprintf(fp,"%s,%d,%d,%.3f,%d,%d,%d,%d,%d,%.3f,%d,%d,%d,%d,%d,%d,%d,%.0f,%.0f,%.0f\n",
				studyEntries[curStudy].patientName.ascii(),
				totalMet,foundMet,sensitivity,totalDetection, tpDetection, fpDetection,
				totalLytic,foundLytic,sensitivityLytic,totalDetection, tpLyticDetection, fpLyticDetection, 
				totalMet10, foundMet10, totalLytic10, foundLytic10,
				segInfo.avg_cortical_intensity, segInfo.avg_spongy_intensity, segInfo.avg_spinal_cord_intensity);
			fclose(fp);
		}

		// write out features
		WriteOutFeatures(globalFeatureFileName);
		//WriteOutFeatures3D(globalFeatureFileName);

		// write out 3D detections
		WriteOutDetections(globalDetectionFileName);

		sprintf(info, "Running time: %ds\n", qtime.elapsed()/1000);
		textLabel_info->setText(info);

		ChangeImageSlice(curSlice);
	}

	QApplication::restoreOverrideCursor();
}

void NIH_BoneSegmentation_Dlg::pushButton_process_clicked()
{
	processingStatus = 0;
	showWSOverlay = false;
	pushButton_pre_process_clicked();
	pushButton_spinal_cord_clicked();
    //pushButton_spine_partition_clicked();
   // pushButton_Rib_Detection_clicked();
	//pushButton_met_detection_clicked();
}

void NIH_BoneSegmentation_Dlg::pushButton_pre_process_clicked()
{
	if(img3D==NULL) return;
	UpdateData();
	
    QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );

	QTime qtime;

	qtime.start();

///	NIH_SpineSegmentation_PreProcess(img3D, maskImg3D, maskStruct, segPara, segInfo);
	NIH_SpineSegmentation_PreProcess_new(img3D, maskImg3D, maskStruct, segPara, segInfo);
	

	char info[100];
	sprintf(info, "Running time: %ds\nCortical bone: %.1f\nSpongy bone: %.1f\nAll bone: %.1f\n", 
		qtime.elapsed()/1000, segInfo.avg_cortical_intensity, segInfo.avg_spongy_intensity, segInfo.avg_bone_intensity);
	textLabel_info->setText(info);

	ChangeImageSlice(curSlice);
    QApplication::restoreOverrideCursor();
	processingStatus = 1;
}

void NIH_BoneSegmentation_Dlg::pushButton_spinal_cord_clicked()
{
	if(img3D==NULL || maskImg3D==NULL) return;
	
	UpdateData();

    QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );

	QTime qtime;

	qtime.start();

	NIH_SpineSegmentation_DetectSpinalCord_new(img3D, maskImg3D, maskStruct, segPara, segInfo, vertebraTemplate, debugMode);
///	NIH_SpineSegmentation_DetectSpinalCord(img3D, maskImg3D, maskStruct, segPara, segInfo, vertebraTemplate, debugMode);
///	NIH_SpineSegmentation_DetectSpinalCord_old(img3D, maskImg3D, maskStruct, segPara, segInfo, debugMode);
	vertebraTemplateViewMode = 1;

	ChangeImageSlice(curSlice);

	char info[100];
	sprintf(info, "Running time: %ds\n", qtime.elapsed()/1000);
	textLabel_info->setText(info);
	processingStatus = 2;

    QApplication::restoreOverrideCursor();

}

void NIH_BoneSegmentation_Dlg::pushButton_met_detection_clicked()
{
	if(img3D==NULL || maskImg3D==NULL) return;
	 
	UpdateData();

    QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );

	QTime qtime;

	qtime.start();

	NIH_BoneSegmentation_DetectMetasis(img3D, maskImg3D, maskStruct, segPara, segInfo, 
		detections2D, numDetections2D, lesions, numLesions, paint,
		flagApplyClassifier, flagApply3DClassifier, svmCommittee, svmCutoff);
	NIH_BoneSegmentation_MergeMetasis(img3D, maskImg3D, maskStruct, segPara, segInfo, detections, numDetections, detectionVoxels);
	
	spinBox_detection->setRange(0, numDetections-1);
	spinBox_detection->setValue(0);

	spinBox_detection_2D->setRange(0, numDetections2D-1);
	spinBox_detection_2D->setValue(0);

	int totalMet, totalLytic, foundMet, foundLytic, totalDetection, tpDetection, fpDetection, 
		tpLyticDetection, fpLyticDetection;
	int totalMet10, totalLytic10, foundMet10, foundLytic10;

	NIH_BoneSegmentation_MatchDetections(maskImg3D, detections, numDetections, detectionVoxels, lesions, numLesions, paint, 
		totalMet, totalLytic, foundMet, foundLytic, totalDetection, tpDetection, fpDetection,
		tpLyticDetection, fpLyticDetection, totalMet10, totalLytic10, foundMet10, foundLytic10, lesionSizeCutoff);

	// reporting
	float sensitivity = 1, sensitivityLytic = 1;
	if(totalMet>0) 
	{
		sensitivity = (float)foundMet/(float)totalMet;
	}
	if(totalLytic>0) 
	{
		sensitivityLytic = (float)foundLytic/(float)totalLytic;
	}
	
	char info[500];
	sprintf(info, "Total Mets: %d (%d); Found: %d (%d);\nTotal Lytic %d (%d); Found: %d (%d)\nDetection: %d; TP: %d; Lytic: %d; FP: %d\n",
		totalMet, totalMet10, foundMet, foundMet10,
		totalLytic, totalLytic10, foundLytic, foundLytic10, 
		totalDetection, tpDetection, tpLyticDetection, fpDetection);

	textLabel_info->setText(info);

	// report
	FILE *fp;
	if((fp = fopen(globalReportFileName, "a"))!=NULL)
	{
		fprintf(fp,"%s,%d,%d,%.3f,%d,%d,%d,%d,%d,%.3f,%d,%d,%d,%d,%d,%d,%d,%.0f,%.0f,%.0f\n",
			studyEntries[curStudy].patientName.ascii(),
			totalMet,foundMet,sensitivity,totalDetection, tpDetection, fpDetection,
			totalLytic,foundLytic,sensitivityLytic,totalDetection, tpLyticDetection, fpLyticDetection, 
			totalMet10, foundMet10, totalLytic10, foundLytic10,
			segInfo.avg_cortical_intensity, segInfo.avg_spongy_intensity, segInfo.avg_spinal_cord_intensity);
		fclose(fp);
	}

	// write out features
	WriteOutFeatures(globalFeatureFileName);

    QApplication::restoreOverrideCursor();
}

void NIH_BoneSegmentation_Dlg::pushButton_classification_clicked()
{
	return;
}


void NIH_BoneSegmentation_Dlg::pushButton_surface_clicked()
{
	if(img3D==NULL || maskImg3D==NULL || spineSurf==NULL || cordSurf==NULL) return;

	UpdateData();

	QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );

	QTime qtime;

	qtime.start();

	double px, py;
	px = img3D->Get_Pixel_SizeX();
	py = img3D->Get_Pixel_SizeY();

	// get binary image from the mask
	int sizex, sizey, sizez,  sizexyz;
	short *maskA;
	int x,y,z;

	sizex = img3D->Num_Cols();
	sizey = img3D->Num_Rows();
	sizez = img3D->Num_Levels();
	if(sizez>500) return;	// constraint the size to address memory issue;

	maskA = maskImg3D->GetArray();
	sizexyz = sizex*sizey*sizez;

	int count;
	int k;
	int lowTh=0;

	lowTh = -1;
	int stx, sty, stz, edx, edy, edz;
	int xx0=10, xx1=500, yy0=200, yy1=380, zz0=0, zz1=270;

	stx = sty = stz = edx =edy =edz = -1;


	if(segInfo.spineBound1.GetSize()>0)
	{
		stx = img3D->Num_Cols();
		sty = img3D->Num_Rows();
		stz = img3D->Num_Levels();

		edx = edy = edz = 0;

		for(z=0; z<img3D->Num_Levels(); z++)
		{
			if(segInfo.spineBound1[z].x==-1) continue;
	
			if(stx>segInfo.spineBound1[z].x) stx=segInfo.spineBound1[z].x;
			if(sty>segInfo.spineBound1[z].y) sty=segInfo.spineBound1[z].y;
			if(edx<segInfo.spineBound2[z].x) edx=segInfo.spineBound2[z].x;
			if(edy<segInfo.spineBound2[z].y) edy=segInfo.spineBound2[z].y;
			if(stz>z) stz=z;
			if(edz<z) edz=z;
		}


		stx--; sty--; stz--;
		stx--; sty--; stz--;
		edx++; edy++; edz++;
		edx++; edy++; edz++;
	}
	
	CIS_Array_Image3D_short *binImg;
	short *binArray;

	CIS_Array_Image3D_short *paintImg;
	short *paintArray;

	binImg = new CIS_Array_Image3D_short(*img3D);
	binArray = binImg->GetArray();
	paintImg = new CIS_Array_Image3D_short(*img3D);
	paintArray = paintImg->GetArray();

	int i, j, k1;

	for(i=0; i<sizexyz; i++)
		paintArray[i]=0;

	for(i=0; i<totalPaintNum; i++)
	{
		for(j=0; j<paint[i].numPaint; j++)
		{
			z = paint[i].sliceArray[j]-1;
			for(k=0; k<paint[i].paint[j].GetSize(); k++)
			{
				x = paint[i].paint[j][k].x;
				y = paint[i].paint[j][k].y;
				k1= z*sizex*sizey+y*sizex+x;
				paintArray[k1]=100;
			}	// for k
		}	// for j
	}	// for i

	stx = img3D->Num_Cols();
	sty = img3D->Num_Rows();
	stz = img3D->Num_Levels();
	edx = edy = edz = 0;
	for(z=0,k=0;z<sizez;z++) 
		for(y=0;y<sizey;y++) 
			for(x=0;x<sizex;x++, k++)
	{
		if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_spongyBone 
			|| maskA[k]==maskStruct.mask_spinalCord || maskA[k]==maskStruct.mask_boneMetasis
			|| maskA[k]==maskStruct.mask_falseDetection) 
		{
			binArray[k]=100;
			if(stx>x) stx=x;
			if(sty>y) sty=y;
			if(edx<x) edx=x;
			if(edy<y) edy=y;
			if(stz>z) stz=z;
			if(edz<z) edz=z;
		}
		else binArray[k]=0;
	}
	stx-=2; sty-=2; stz-=2;
	edx+=2; edy+=2; edz+=2;

	xx0=stx; xx1=edx; yy0=sty; yy1=edy; zz0=stz; zz1=edz;

	// compute the iso surface
	spineSurf->RemoveAllMesh();
	
	CIS_Algo_MarchingCube(50, -1, lowTh, 3, 1, binImg, NULL, spineSurf, 
		stx, edx, sty, edy, stz, edz);

//		spineSurf->colorTable.SetSize(2);
//		spineSurf->colorTable[1] = Vec3(0,1,0);
//		int numFace=spineSurf->NumFaces();
//		spineSurf->faceColorIndex.SetSize(numFace);
//		for(k=0; k<numFace; k++) spineSurf->faceColorIndex[k]=0;

	TriangleMesh_LaplacianSmooth(spineSurf, 1, 10);
	NIH_OpenGL_CentralizeModel(display_id2, spineSurf);

	spineSurf->SetCurrColor(Vec3(0.91, 0.45, 0.35));

	// surface for spinal cord

	stx = img3D->Num_Cols();
	sty = img3D->Num_Rows();
	stz = img3D->Num_Levels();
	edx = edy = edz = 0;
	for(z=0,k=0;z<sizez;z++) for(y=0;y<sizey;y++) for(x=0;x<sizex;x++, k++)
	{
		if(maskA[k]==maskStruct.mask_spinalCord) 
		{
			binArray[k]=100;
			if(stx>x) stx=x;
			if(sty>y) sty=y;
			if(edx<x) edx=x;
			if(edy<y) edy=y;
			if(stz>z) stz=z;
			if(edz<z) edz=z;
		}
		else binArray[k]=0;
	}
	stx-=2; sty-=2; stz-=2;
	edx+=2; edy+=2; edz+=2;

	// compute the iso surface
	cordSurf->RemoveAllMesh();
	
	CIS_Algo_MarchingCube(50, -1, lowTh, 3, 1, binImg, NULL, cordSurf, 
		stx, edx, sty, edy, stz, edz);

	TriangleMesh_LaplacianSmooth(cordSurf, 1, 10);

	cordSurf->SetCurrColor(Vec3(0, 0, 1));

	count=0;
	// surface for metasis
	for(k=0; k<sizexyz; k++)
	{
		if(maskA[k]==maskStruct.mask_boneMetasis) 
		{
			binArray[k]=100;
			count++;
		}
		else binArray[k]=0;
	}

	if(count>0)
	{
		// compute the iso surface
		metasisSurf->RemoveAllMesh();
	
		CIS_Algo_MarchingCube(50, -1, lowTh, 3, 1, binImg, NULL, metasisSurf, 
			stx, edx, sty, edy, stz, edz);

		TriangleMesh_LaplacianSmooth(metasisSurf, 1, 10);

		metasisSurf->SetCurrColor(Vec3(0, 1, 0));

		if(flagShowOverlay) metasisSurf->SetDrawMode(SURFACE_MODE);
		else metasisSurf->SetDrawMode(HIDE_MODE);
	} else metasisSurf->SetDrawMode(HIDE_MODE);


	stx = img3D->Num_Cols();
	sty = img3D->Num_Rows();
	stz = img3D->Num_Levels();
	edx = edy = edz = 0;
	for(z=0,k=0;z<sizez;z++) 
		for(y=0;y<sizey;y++) 
			for(x=0;x<sizex;x++, k++)
	{
		if(paintArray[k]==100) 
		{
			binArray[k]=100;
			if(stx>x) stx=x;
			if(sty>y) sty=y;
			if(edx<x) edx=x;
			if(edy<y) edy=y;
			if(stz>z) stz=z;
			if(edz<z) edz=z;
		}
		else binArray[k]=0;
	}
	stx-=2; sty-=2; stz-=2;
	edx+=2; edy+=2; edz+=2;

	// surface for painting
	paintSurf->RemoveAllMesh();

	CIS_Algo_MarchingCube(50, -1, lowTh, 3, 1, binImg, NULL, paintSurf, 
			stx, edx, sty, edy, stz, edz);

	TriangleMesh_LaplacianSmooth(paintSurf, 1, 10);

	paintSurf->SetCurrColor(Vec3(0,1,1));

	if(flagShowOverlay) paintSurf->SetDrawMode(SURFACE_MODE);
	else paintSurf->SetDrawMode(HIDE_MODE);

	// disk planes
	int p;
	Vec3 normal, dirX, dirY;
	double _sizeX, _sizeY;
	for(p=0; p<30; p++)
	{
		diskPlaneModel[p]->SetDrawMode(HIDE_MODE);
	}
	for(p=0; p<pedicleValley.GetSize() && p<30; p++)
	{
		z=pedicleValley[p];
		dirY = spineNormalX[z];
		dirX = spineNormalY[z];
		
		// image orientation is reversed
		if(img3D->GetZOrient()==1)
		{
			dirY.z = -dirY.z;
			dirX.z = -dirX.z;
		}

		normal = dirX % dirY;
		dirX = dirY % normal;
		_sizeX = (spineWidthLeft[z]+spineWidthRight[z])*1.5;
		_sizeY = (spineWidthUp[z]+spineWidthDown[z])*1.5;
		diskPlaneModel[p]->SetDimension(spineCenter[z], normal, dirX, _sizeX, _sizeY);
		if(flagShowBBox) diskPlaneModel[p]->SetDrawMode(SURFACE_MODE);
		else diskPlaneModel[p]->SetDrawMode(HIDE_MODE);

		if((normal*(spineCenter[z+1]-spineCenter[z]))<0) normal = -normal;
		frontPlane[p].Set(spineCenter[z], normal);

		dirX = spineNormalY[z];
		dirY = spineNormalBack[z];
		// image orientation is reversed
		if(img3D->GetZOrient()==1)
		{
			dirY.z = -dirY.z;
			dirX.z = -dirX.z;
		}
		normal = dirX % dirY;
		if((normal*(spineCenter[z+1]-spineCenter[z]))<0) normal = -normal;
		backPlane[p].Set(spineCenter[z], normal);
	}
	
	// show the ribs
	if(segInfo.t12_vertebra>=0)
	{
		count=0;
		stx = img3D->Num_Cols();
		sty = img3D->Num_Rows();
		stz = img3D->Num_Levels();
		edx = edy = edz = 0;
		
		// surface for ribs
		for(z=0,k=0;z<sizez;z++) for(y=0;y<sizey;y++) for(x=0;x<sizex;x++, k++)
		{
			if(maskA[k]==maskStruct.mask_rib) 
			{
				binArray[k]=100;
				if(stx>x) stx=x;
				if(sty>y) sty=y;
				if(edx<x) edx=x;
				if(edy<y) edy=y;
				if(stz>z) stz=z;
				if(edz<z) edz=z;
				count++;
			}
			else binArray[k]=0;
		}
		stx-=2; sty-=2; stz-=2;
		edx+=2; edy+=2; edz+=2;
		if(stx<xx0) xx0=stx;
		if(sty<yy0) yy0=sty;
		if(edx>xx1) xx1=edx;
		if(edy>yy1) yy1=edy;

		if(count>0)
		{
			ribSurf->RemoveAllMesh();
	
			CIS_Algo_MarchingCube(50, -1, lowTh, 3, 1, binImg, NULL, ribSurf, 
				stx, edx, sty, edy, stz, edz);

			TriangleMesh_LaplacianSmooth(ribSurf, 1, 10);

			ribSurf->SetCurrColor(Vec3(0, 1, 0));

			ribSurf->SetDrawMode(SURFACE_MODE);
		} else ribSurf->SetDrawMode(HIDE_MODE);

		// reset the color of surface based on its relative location to the partition plane
		spineSurf->colorTable.SetSize(10);
		spineSurf->colorTable[1] = Vec3(0.6,1,0.4);
		spineSurf->colorTable[2] = Vec3(0.3,1,0.7);
		spineSurf->colorTable[3] = Vec3(0.8,1,0.2);
		spineSurf->colorTable[4] = Vec3(0,1,0);
		spineSurf->colorTable[5] = Vec3(1,0,0);
		spineSurf->colorTable[6] = Vec3(1,0.8,0.2);
		spineSurf->colorTable[7] = Vec3(1,0.3,0.7);
		spineSurf->colorTable[8] = Vec3(1,0.6,0.4);
		spineSurf->colorTable[9] = Vec3(1,0.1,0.9);
		
		int vertNum=spineSurf->NumVertices();
		int faceNum=spineSurf->NumFaces();
		Vec3 v;
		double d1, d2;
		int ci;
		spineSurf->faceColorIndex.SetSize(faceNum);
		for(int i=0; i<faceNum; i++)
		{
			v = spineSurf->FaceCentroid(i);
			z = img3D->GetSliceLevel(v.z, 1);
			p=0;
			if(v.y>spineCenter[z].y) d1=backPlane[p].SignDist(v);
			else d1 = frontPlane[p].SignDist(v);

			if(d1>0)
			{
				for(p=1; p<pedicleValley.GetSize() && p<30; p++)
				{
					if(v.y>spineCenter[z].y) d2=backPlane[p].SignDist(v);
					else d2=frontPlane[p].SignDist(v);
					if(d1*d2<=0) break;
					d1=d2;
				}
			}

			ci = p-segInfo.t12_vertebra+4;

			spineSurf->faceColorIndex[i]=ci;
		}
	} else ribSurf->SetDrawMode(HIDE_MODE);


	// Setup volume rendering
	if(volumeModel==NULL)
	{
		int sizex, sizey, sizez;
		float px, py, pz;
		sizex = img3D->Num_Cols();
		sizey = img3D->Num_Rows();
		sizez = img3D->Num_Levels();
		px = img3D->Get_Pixel_SizeX();
		py = img3D->Get_Pixel_SizeY();
		pz = img3D->Get_Pixel_SizeZ();

		volumeModel = new CIS_3D_Model_Volume();
		CIS_Array_Image2D_short *imgT_2D;
		imgT_2D = new CIS_Array_Image2D_short(sizex, sizez);
		Vec3 sliceLoc;
		double dimx, dimy;
		Vec3 normal, dirX;
		dimx = sizex*px;
		dimy = sizez*pz;
		normal = Vec3(0,-1,0);
		dirX = Vec3(1,0,0);
		sliceLoc.x = px*sizex/2;
		sliceLoc.z = img3D->GetSlicePosition(sizez/2);
		volumeModel->SetTotalTextures(yy1-yy0+1);
		int x, z, k;
		short *img2DA;

		zz1 += 5;
		if(xx0<0) xx0=0; if(xx1>sizex-1) xx1=sizex-1;
		if(yy0<0) yy0=0; if(yy1>sizey-1) yy1=sizey-1;
		if(zz0<0) zz0=0; if(zz1>sizez-1) zz1=sizez-1;
		for(int y=yy1, i=0; y>=yy0; y--, i++)
		{
			img3D->GetSliceImage(*imgT_2D, y, 1, true);
			img2DA=imgT_2D->GetArray();
			sliceLoc.y=y*py;
			for(z=0, k=0; z<sizez; z++)
			{
				for(x=0; x<sizex; x++, k++)
				{
					if(x<xx0 || x>xx1 || z<zz0 || z>zz1) img2DA[k]=0;
				}
			}
			volumeModel->SetOneTexture(i, imgT_2D, sliceLoc, normal, dirX, dimx, dimy);
		}

		delete imgT_2D;
		volumeModel->SetRenderingParameter(1300, 0.1, 300, 1200);
		CIS_Model_AddModel(volumeModel);
		NIH_OpenGL_AddModel(display_id2, volumeModel);	
		volumeModel->SetDrawMode(LINE_MODE);
		volumeModel->SetCurrAlpha(0.2);
	}

	NIH_OpenGL_Refresh_Display(display_id2);

	char info[100];
	sprintf(info, "Running time: %ds\n", qtime.elapsed()/1000);
	textLabel_info->setText(info);

	delete binImg;
	QApplication::restoreOverrideCursor();
}

void NIH_BoneSegmentation_Dlg::pushButton_help_clicked()
{

	myQFileDialog *fd = new myQFileDialog(this, "File", true);
	fd->exec();
	QString allDirs = fd->selectedDirs.join("\n");
	printf("%s\n", allDirs.ascii());

 

/*	char helpText[5000];
	strcpy(helpText, "");
	strcat(helpText, "'Load' Button: Load studies\n");
	strcat(helpText, "'Process' Button: Process the image\n");
	strcat(helpText, "'Close' Button: Close the panel\n\n");
	strcat(helpText, "'Slice' Slider: Switch slices\n");
	strcat(helpText, "'Contrast' Slider: Change Window\\Level\n");
	strcat(helpText, "'Opacity' Slider: Change transparency\n\n");
	strcat(helpText, "Left Mouse + Ctrl: Define ROI\n");
	strcat(helpText, "key 'C': Clear ROI\n");
	strcat(helpText, "key 'V': equal 'Overlay' button\n");
	QMessageBox::information(NULL, "Help", helpText);
*/

}


void NIH_BoneSegmentation_Dlg::pushButton_refresh_clicked()
{
	UpdateData();
	if(viewMode==2) ChangeImageSlice(curSlice);
	if(viewMode==0) ChangeImageSlice_x(curSlice_x);
	if(viewMode==1) ChangeImageSlice_y(curSlice_y);
}


void NIH_BoneSegmentation_Dlg::pushButton_close_clicked()
{
//	if(segmentationChanged)
//	{
//		if(QMessageBox::question(NULL, "Warning", "Does you want to save the changes (Y/N)?", 
//			QMessageBox::Yes , QMessageBox::No)==QMessageBox::Yes)
//		{
//			pushButton_report_clicked();
//		}
//	}

	this->close(true);
	QApplication::exit();

}


void NIH_BoneSegmentation_Dlg::pushButton_load_project_clicked()
{
	UpdateData();

	QString filter = "Project File (*.prj)";
  
    QFileDialog* cf = new QFileDialog( this );
    cf->addFilter( filter);
	filter = "All Files (*.prj)";
    cf->addFilter( filter);
 	
	if(cf->exec()==QDialog::Accepted) 
	{
		QString strNewFileName;
	 	char proj_fn[200];
	
		strNewFileName=cf->selectedFile();

		strcpy(proj_fn,strNewFileName.ascii());

		load_project(proj_fn);
	}

	delete cf;
}


int NIH_BoneSegmentation_Dlg::load_project(const char* proj_fn)
{
 	char line[2000], value_s[2000], key_s[2000];
	QStringList studyLine;
	FILE *fp;

	numStudyEntries = 0;
	batchIndex = 1;
	if((fp=fopen(proj_fn, "r"))!=NULL)
	{
		flagApplyClassifier = false;
		flagApply3DClassifier = false;
		while(!feof(fp))
		{
			strcpy(line,"");
			fgets(line, 2000, fp);
			if(strlen(line)<2) continue;
			if(line[0]=='-' && line[1]=='1') break;

			if(strlen(line)>2 && line[0]!='#' && line[0]!=' ')
			{
				if(line[0]=='!')	// readin key-value pair
				{	
					sscanf(line, "%s %s", key_s, value_s);
					// readin PatientRootPath
					if(!strcmp(key_s, "!PathToPatient"))
					{
						strcpy(patientRoot, value_s);
						if(patientRoot[strlen(patientRoot)-1]!='\\')
							strcat(patientRoot,"\\");
					}
					if(!strcmp(key_s, "!PathToPainting"))
					{
						strcpy(paintingRoot, value_s);
						if(paintingRoot[strlen(paintingRoot)-1]!='\\')
							strcat(paintingRoot,"\\");
					}
					// readin LoadSegPath
					if(!strcmp(key_s, "!LoadSegPath"))
					{
						strcpy(loadSegPath, value_s);
						if(loadSegPath[strlen(loadSegPath)-1]!='\\')
							strcat(loadSegPath,"\\");
					}
					// readin SaveSegPath
					if(!strcmp(key_s, "!SaveSegPath"))
					{
						strcpy(saveSegPath, value_s);
						if(saveSegPath[strlen(saveSegPath)-1]!='\\')
							strcat(saveSegPath,"\\");
					}
					// readin ReportFileName
					else if(!strcmp(key_s, "!ReportFileName"))
					{
						strcpy(globalReportFileName, value_s);
						// write some header
						FILE *fp;
						if((fp = fopen(globalReportFileName, "a"))!=NULL)
						{
							fprintf(fp,"Date: %s, Time: %s\n", QDate::currentDate().toString().ascii(), QTime::currentTime().toString().ascii());
							fprintf(fp,"Model: %s, Cutoff: %f\n", svmModelFileName.ascii(), svmCutoff);
							fprintf(fp,"patientName,totalMet,foundMet,sensitivity,totalDetection,tpDetection,fpDetection,totalDetection\n");
							fclose(fp);
						}

					}
					else if(!strcmp(key_s, "!FeatureFileName"))
					{
						strcpy(globalFeatureFileName, value_s);
					}
					else if(!strcmp(key_s, "!DetectionFileName"))
					{
						strcpy(globalDetectionFileName, value_s);
					}
					else if(!strcmp(key_s, "!svmModelFileName"))
					{
						svmModelFileName = value_s;
						// read in the svm model
					    if(svmCommittee->ReadXmlModelFile((char *)svmModelFileName.ascii())==0 ) 
						{
							flagApplyClassifier = true;
						}
						else flagApplyClassifier = false;
					}
					else if(!strcmp(key_s, "!svmCutoff"))
					{
						svmCutoff = QString(value_s).toFloat();
					}
					else if(!strcmp(key_s, "!svmCutoffMin"))
					{
						svmCutoffMin = QString(value_s).toFloat();
					}
					else if(!strcmp(key_s, "!svmCutoffMax"))
					{
						svmCutoffMax = QString(value_s).toFloat();
					}
					else if(!strcmp(key_s, "!svmCutoffStep"))
					{
						svmCutoffStep = QString(value_s).toFloat();
					}
					else if(!strcmp(key_s, "!lesionSizeCutoff"))
					{
						lesionSizeCutoff = QString(value_s).toFloat();
					}
					else if(!strcmp(key_s, "!useMethod"))
					{
						useMethod = QString(value_s).toInt();
					}
					else if(!strcmp(key_s, "!batchIndex"))
					{
						batchIndex = QString(value_s).toInt();
					}
				}	// if line[0]
				else	// readin entries
				{
					// each entry should be in the form of
					// patientName, studyDate, studyId, seriesId, localImagePath
					studyLine = QStringList::split(',', line, false);
					if(studyLine.size()<5) continue;
					for(int i=0; i<studyLine.size(); i++) studyLine[i]=studyLine[i].simplifyWhiteSpace();
					studyEntries[numStudyEntries].patientName = studyLine[0];
					studyEntries[numStudyEntries].studyDate = studyLine[1];
					studyEntries[numStudyEntries].studyId = studyLine[2];
					studyEntries[numStudyEntries].seriesId = studyLine[3];
					studyEntries[numStudyEntries].localImagePath = studyLine[4];

					studyEntries[numStudyEntries].paintingFileName = "NULL";
					studyEntries[numStudyEntries].organPath = "NULL";
					if(batchIndex>5)
					{
						if(studyLine.size()>5) studyEntries[numStudyEntries].organPath = studyLine[5];
						else studyEntries[numStudyEntries].organPath = "NULL";
					}
					else
					{
						if(studyLine.size()>5) studyEntries[numStudyEntries].paintingFileName = studyLine[5];
						else studyEntries[numStudyEntries].paintingFileName = "NULL";
						if(studyLine.size()>6) studyEntries[numStudyEntries].sliceRange1 = studyLine[6].toInt()-1;
						else studyEntries[numStudyEntries].sliceRange1 = -1;
						if(studyLine.size()>7) studyEntries[numStudyEntries].sliceRange2 = studyLine[7].toInt()-1;
						else studyEntries[numStudyEntries].sliceRange2 = -1;
					}

					numStudyEntries ++;
				}	// else 
			}	// if strlen
		}	// while feof

		// fill in the list table
		char buf[2000];
		listBox_study_entry->clear();
		for(int i=0; i<numStudyEntries; i++)
		{
			sprintf(buf, "%-3d -- %15s", 
				i,
				studyEntries[i].patientName.ascii());
			listBox_study_entry->insertItem(buf);
		}

		curStudy = -1;
		batchStudyStart = 0;
		batchStudyEnd = numStudyEntries-1;
		lineEdit_batch_start->setText(QString::number(batchStudyStart));
		lineEdit_batch_end->setText(QString::number(batchStudyEnd));
		fclose(fp);

		return CIS_OK;
	}	// if fp

	return CIS_ERROR;
}
void NIH_BoneSegmentation_Dlg::pushButton_batch_project_clicked()
{
	UpdateData();

	flagBatchMode = true;

	for(curStudy=batchStudyStart; curStudy>=0 && curStudy<numStudyEntries && curStudy<=batchStudyEnd; curStudy++)
	{
		// clean out all meshes
		if(spineSurf) spineSurf->RemoveAllMesh();
		if(cordSurf) cordSurf->RemoveAllMesh();
		if(metasisSurf) metasisSurf->RemoveAllMesh();
		if(ribSurf) ribSurf->RemoveAllMesh();

		// put a report on batch processing
		FILE *fp;
		fp = fopen("c:\\tmp\\batch_progress.txt","a+");
		fprintf(fp, "%s -- Bone project: batch %d, study %d, patient %s\n", 
			QDateTime::currentDateTime().toString("yyyy-MM-dd:hh:mm:ss").ascii(),
			batchIndex, curStudy, studyEntries[curStudy].patientName.ascii());
		fclose(fp);

		segInfo.Initialize();
		LoadStudy(curStudy);

		if(batchIndex==1)
		{
			// batch 1 process the data and detect metastasis
//			NIH_BoneSegmentation_LoadSegmentation(loadSegPath, studyEntries[curStudy].localImagePath, maskImg3D);

			LoadPaintingFile(studyEntries[curStudy].paintingFileName.ascii());
	
			pushButton_pre_process_clicked();
			pushButton_spinal_cord_clicked();

			if(useMethod==0) pushButton_met_detection_clicked();
			else if(useMethod==1) pushButton_watershed_clicked();
 
			pushButton_save_seg_clicked();
		}
		else if(batchIndex==2)
		{
			// batch 2: segment and partition the spinal column, and save the result
			pushButton_pre_process_clicked();
			pushButton_spinal_cord_clicked();
			pushButton_save_seg_clicked();
			pushButton_spine_partition_clicked();

			char ss_fn[400];
			NIH_OpenGL_Show_Background(display_id, false);
			// screen shot of sagittal view
			ChangeImageSlice_x(curSlice_x);
			sprintf(ss_fn, "c:\\tmp\\ScreenShot\\%s_sag.ppm", studyEntries[curStudy].patientName.ascii());
			NIH_OpenGL_Save_Display(display_id, ss_fn);
	 	     
			// screen shot of coronal view
			ChangeImageSlice_y(curSlice_y);
			sprintf(ss_fn, "c:\\tmp\\ScreenShot\\%s_cor.ppm", studyEntries[curStudy].patientName.ascii());
			NIH_OpenGL_Save_Display(display_id, ss_fn);
			
			// may need a 3D view
		}		
		else if(batchIndex==3)
		{
			// batch 3: exam the result of segmentation, generate curved reformation and surface
			//
			NIH_BoneSegmentation_LoadSegmentation(loadSegPath, studyEntries[curStudy].localImagePath, maskImg3D);
			pushButton_spine_partition_clicked();

			char ss_fn[400];
			NIH_OpenGL_Show_Background(display_id, false);
			// screen shot of sagittal view
			ChangeImageSlice_x(curSlice_x);
			sprintf(ss_fn, "c:\\tmp\\ScreenShot\\%s_sag.ppm", studyEntries[curStudy].patientName.ascii());
			NIH_OpenGL_Save_Display(display_id, ss_fn);
	 	     
			// screen shot of coronal view
			ChangeImageSlice_y(curSlice_y);
			sprintf(ss_fn, "c:\\tmp\\ScreenShot\\%s_cor.ppm", studyEntries[curStudy].patientName.ascii());
			NIH_OpenGL_Save_Display(display_id, ss_fn);

			// may need a 3D view
		}
		else if(batchIndex==4)	// bmd and vertebra labelling
		{
			char ss_fn[400];
			char screenshot_root[200];
			sprintf(screenshot_root, "%s\\screenshot\\", study_root);

			// spine segmentation, partition, labeling, plus rib detection and ROI determination
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  1..");
			fclose(fp);			
			pushButton_pre_process_clicked();
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  2..");
			fclose(fp);			
			pushButton_spinal_cord_clicked();
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  3..");
			fclose(fp);			
			pushButton_spine_partition_clicked();
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  4..");
			fclose(fp);			
			pushButton_Rib_Detection_clicked();

			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  5..\n");
			fclose(fp);	

			flagShowBBox = true;
			checkBox_bounding_box->setChecked(flagShowBBox);

			// create the direction if not exist
			QDir dcmdir = QString(screenshot_root);
			if(!dcmdir.exists())
			{
				dcmdir.mkdir(screenshot_root);
			}


			cordModel->SetDrawMode(HIDE_MODE);
			smoothedCord->SetDrawMode(HIDE_MODE);

			pushButton_surface_clicked();
			sprintf(ss_fn, "%s%s_3d.jpg", screenshot_root,study_name);
			NIH_OpenGL_Save_Display(display_id2, ss_fn);

			cordModel->SetDrawMode(LINE_MODE);
			smoothedCord->SetDrawMode(LINE_MODE);


			NIH_OpenGL_Show_Background(display_id, false);
			// screen shot of sagittal view
			ChangeImageSlice_x(curSlice_x);
			sprintf(ss_fn, "%s%s_sag.jpg", screenshot_root,study_name);
			NIH_OpenGL_Save_Display(display_id, ss_fn);
	 	     
			// screen shot of coronal view
			ChangeImageSlice_y(curSlice_y);
			sprintf(ss_fn, "%s%s_cor.jpg", screenshot_root,study_name);
			NIH_OpenGL_Save_Display(display_id, ss_fn);

			NIH_OpenGL_Show_Background(display_id, true);
			flagShowBBox = false;
			checkBox_bounding_box->setChecked(flagShowBBox);
			pushButton_BMD_ROI_clicked();
			// screen shot of each BMD_ROI
			for(int s=0; s<roiSlices.GetSize(); s++)
			{
				ChangeImageSlice(roiSlices[s]);
				sprintf(ss_fn, "%s%s_sl%d.jpg", screenshot_root,study_name,roiSlices[s]+1);
				NIH_OpenGL_Save_Display(display_id, ss_fn);
			}

/*			pushButton_Vertebra_Seg_clicked();
			flagShowBBox = true;
			checkBox_bounding_box->setChecked(flagShowBBox);
			flagShowOverlay = true;
			checkBox_overlay->setChecked(flagShowOverlay);

			NIH_OpenGL_Show_Background(display_id, false);
			// screen shot of sagittal view
			ChangeImageSlice_x(curSlice_x);
			sprintf(ss_fn, "%s%s_vetHeight.jpg", screenshot_root,study_name);
			NIH_OpenGL_Save_Display(display_id, ss_fn);
*/		}
		else if(batchIndex==5)
		{
			// spine segmentation, partition, labeling, plus rib detection and vertebra height measurement
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  1..");
			fclose(fp);			
 			pushButton_pre_process_clicked();
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  2..");
			fclose(fp);			
			pushButton_spinal_cord_clicked();
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  3..");
			fclose(fp);			
			pushButton_spine_partition_clicked();
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  4..");
			fclose(fp);			
			pushButton_Rib_Detection_clicked();

			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  5..\n");
			fclose(fp);			
			pushButton_Vertebra_Seg_clicked();
			flagShowBBox = true;
			checkBox_bounding_box->setChecked(flagShowBBox);
			flagShowOverlay = true;
			checkBox_overlay->setChecked(flagShowOverlay);

			char screenshot_root[200];
			sprintf(screenshot_root, "%s\\screenshot\\", study_root);
			// create the directory if not exist
			QDir dcmdir = QString(screenshot_root);
			if(!dcmdir.exists())
			{
				dcmdir.mkdir(screenshot_root);
			}

			char ss_fn[400];
			NIH_OpenGL_Show_Background(display_id, false);
			// screen shot of sagittal view
			ChangeImageSlice_x(curSlice_x);
			sprintf(ss_fn, "%s%s_vetHeight.jpg", screenshot_root,study_name);
			NIH_OpenGL_Save_Display(display_id, ss_fn);
		}
		else if(batchIndex==6)
		{
			// spine segmentation, partition, labeling, plus rib detection 
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  1..");
			fclose(fp);			
			pushButton_pre_process_clicked();
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  2..");
			fclose(fp);			
			pushButton_spinal_cord_clicked();
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  3..");
			fclose(fp);			
			pushButton_spine_partition_clicked();
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  4..");
			fclose(fp);			
			pushButton_Rib_Detection_clicked();

			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  5..\n");
			fclose(fp);	

			flagShowBBox = false;
			checkBox_bounding_box->setChecked(flagShowBBox);
			char screenshot_root[200];
///			sprintf(screenshot_root, "%s\\screenshot\\", study_root);
			sprintf(screenshot_root, "c:\\tmp\\screenshot10\\", study_root);
			// create the direction if not exist
			QDir dcmdir = QString(screenshot_root);
			if(!dcmdir.exists())
			{
				dcmdir.mkdir(screenshot_root);
			}

			cordModel->SetDrawMode(HIDE_MODE);
			smoothedCord->SetDrawMode(HIDE_MODE);

			char ss_fn[400];

			pushButton_surface_clicked();
			sprintf(ss_fn, "%s%s_3d.jpg", screenshot_root,study_name);
			NIH_OpenGL_Save_Display(display_id2, ss_fn);

			NIH_OpenGL_Show_Background(display_id, false);
/*			// screen shot of coronal view
			ChangeImageSlice_y(curSlice_y);
			sprintf(ss_fn, "%s%s_cor.jpg", screenshot_root,study_name);
			NIH_OpenGL_Save_Display(display_id, ss_fn);

			// screen shot of sagittal view
			ChangeImageSlice_x(curSlice_x);
			sprintf(ss_fn, "%s%s_sag.jpg", screenshot_root,study_name);
			NIH_OpenGL_Save_Display(display_id, ss_fn);
*/
			flagShowPainting=false;
			flagShowBBox = true;
			ChangeImageSlice_x(curSlice_x);
			sprintf(ss_fn, "%s%s_sag1.jpg", screenshot_root,study_name);
			NIH_OpenGL_Save_Display(display_id, ss_fn);

		}
		else if(batchIndex==7)	// build location model
		{
			pushButton_Build_ALM_clicked();
			char fn[300];
			strcpy(fn, "c:\\Experiments\\bone_training\\Atlas\\");
			strcat(fn, studyEntries[curStudy].patientName.ascii());
			strcat(fn, ".xml");
			SaveSingleLocationModel(fn);
		}
		else if(batchIndex==9)	// output image information
		{
			FILE *fp;
			fp = fopen("c:\\tmp\\patient_stat.txt", "a+");
			fprintf(fp, "%s, %s, %s\n", studyEntries[curStudy].patientName.ascii(), patientSex.ascii(), patientAge.ascii());
			fclose(fp);
		}
		else if(batchIndex==10)	// BMD using manually identified slice
		{
			char ss_fn[400];
			char screenshot_root[200];
			sprintf(screenshot_root, "%s\\screenshot\\", study_root);

			// spine segmentation, partition, labeling, plus rib detection and ROI determination
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  1..");
			fclose(fp);			
			pushButton_pre_process_clicked();
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  2..");
			fclose(fp);			
			pushButton_spinal_cord_clicked();
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  3..");
			fclose(fp);			
/*			pushButton_spine_partition_clicked();
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  4..");
			fclose(fp);			
			pushButton_Rib_Detection_clicked();

			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  5..\n");
			fclose(fp);	
*/
			flagShowBBox = true;
			checkBox_bounding_box->setChecked(flagShowBBox);

			// create the direction if not exist
			QDir dcmdir = QString(screenshot_root);
			if(!dcmdir.exists())
			{
				dcmdir.mkdir(screenshot_root);
			}


			cordModel->SetDrawMode(HIDE_MODE);
			smoothedCord->SetDrawMode(HIDE_MODE);

///			pushButton_surface_clicked();
///			sprintf(ss_fn, "%s%s_3d.jpg", screenshot_root,study_name);
///			NIH_OpenGL_Save_Display(display_id2, ss_fn);

			cordModel->SetDrawMode(LINE_MODE);
			smoothedCord->SetDrawMode(LINE_MODE);


			NIH_OpenGL_Show_Background(display_id, false);
			// screen shot of sagittal view
			ChangeImageSlice_x(curSlice_x);
			sprintf(ss_fn, "%s%s_sag.jpg", screenshot_root,study_name);
			NIH_OpenGL_Save_Display(display_id, ss_fn);
	 	     
			// screen shot of coronal view
			ChangeImageSlice_y(curSlice_y);
			sprintf(ss_fn, "%s%s_cor.jpg", screenshot_root,study_name);
			NIH_OpenGL_Save_Display(display_id, ss_fn);

			NIH_OpenGL_Show_Background(display_id, true);
			flagShowBBox = false;
			checkBox_bounding_box->setChecked(flagShowBBox);
			pushButton_BMD_ROI_L1L2_clicked();		// for manually selected slices
///			pushButton_BMD_ROI_clicked();
			// screen shot of each BMD_ROI
			for(int s=0; s<roiSlices.GetSize(); s++)
			{
				ChangeImageSlice(roiSlices[s]);
				sprintf(ss_fn, "%s%s_sl%d.jpg", screenshot_root,study_name,roiSlices[s]+1);
				NIH_OpenGL_Save_Display(display_id, ss_fn);
			}

/*			pushButton_Vertebra_Seg_clicked();
			flagShowBBox = true;
			checkBox_bounding_box->setChecked(flagShowBBox);
			flagShowOverlay = true;
			checkBox_overlay->setChecked(flagShowOverlay);

			NIH_OpenGL_Show_Background(display_id, false);
			// screen shot of sagittal view
			ChangeImageSlice_x(curSlice_x);
			sprintf(ss_fn, "%s%s_vetHeight.jpg", screenshot_root,study_name);
			NIH_OpenGL_Save_Display(display_id, ss_fn);
*/		}
		else if(batchIndex==11)	// bmd and vertebra labelling using 3D ROI of entire vertebra body
		{
			char ss_fn[400];
			char screenshot_root[200];
			sprintf(screenshot_root, "%s\\screenshot\\", study_root);

			// spine segmentation, partition, labeling, plus rib detection and ROI determination
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  1..");
			fclose(fp);			
			pushButton_pre_process_clicked();
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  2..");
			fclose(fp);			
			pushButton_spinal_cord_clicked();
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  3..");
			fclose(fp);			
			pushButton_spine_partition_clicked();
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  4..");
			fclose(fp);			
			pushButton_Rib_Detection_clicked();

			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  5..\n");
			fclose(fp);	

			flagShowBBox = true;
			checkBox_bounding_box->setChecked(flagShowBBox);

			// create the direction if not exist
			QDir dcmdir = QString(screenshot_root);
			if(!dcmdir.exists())
			{
				dcmdir.mkdir(screenshot_root);
			}


			cordModel->SetDrawMode(HIDE_MODE);
			smoothedCord->SetDrawMode(HIDE_MODE);

			pushButton_surface_clicked();
			sprintf(ss_fn, "%s%s_3d.jpg", screenshot_root,study_name);
			NIH_OpenGL_Save_Display(display_id2, ss_fn);

			cordModel->SetDrawMode(LINE_MODE);
			smoothedCord->SetDrawMode(LINE_MODE);


			NIH_OpenGL_Show_Background(display_id, false);
			// screen shot of sagittal view
			ChangeImageSlice_x(curSlice_x);
			sprintf(ss_fn, "%s%s_sag.jpg", screenshot_root,study_name);
			NIH_OpenGL_Save_Display(display_id, ss_fn);
	 	     
			// screen shot of coronal view
			ChangeImageSlice_y(curSlice_y);
			sprintf(ss_fn, "%s%s_cor.jpg", screenshot_root,study_name);
			NIH_OpenGL_Save_Display(display_id, ss_fn);

			NIH_OpenGL_Show_Background(display_id, true);
			flagShowBBox = false;
			checkBox_bounding_box->setChecked(flagShowBBox);
			pushButton_BMD_ROI_3D_clicked();
			// screen shot of each BMD_ROI
/*			for(int s=0; s<roiSlices.GetSize(); s++)
			{
				ChangeImageSlice(roiSlices[s]);
				sprintf(ss_fn, "%s%s_sl%d.jpg", screenshot_root,study_name,roiSlices[s]+1);
				NIH_OpenGL_Save_Display(display_id, ss_fn);
			}
*/			
			// screen shot of each BMD_ROI slice between L1 and L2
			int st_p;
			if(segInfo.t12_vertebra>=0) st_p=segInfo.t12_vertebra;
			else st_p=0;
			int vertNum = pedicleValley.GetSize();
			// L1
			for(int p=pedicleValley[st_p];st_p+1<vertNum && p<pedicleValley[st_p+1]; p++)
			{
				ChangeImageSlice(p);
				sprintf(ss_fn, "%s%s_L1_%d.jpg", screenshot_root,study_name,p+1);

				NIH_OpenGL_Save_Display(display_id, ss_fn);
			}

			// L2
			for(int p=pedicleValley[st_p+1];st_p+2<vertNum && p<pedicleValley[st_p+2]; p++)
			{
				ChangeImageSlice(p);
				sprintf(ss_fn, "%s%s_L2_%d.jpg", screenshot_root,study_name,p+1);

				NIH_OpenGL_Save_Display(display_id, ss_fn);
			}

/*			pushButton_Vertebra_Seg_clicked();
			flagShowBBox = true;
			checkBox_bounding_box->setChecked(flagShowBBox);
			flagShowOverlay = true;
			checkBox_overlay->setChecked(flagShowOverlay);

			NIH_OpenGL_Show_Background(display_id, false);
			// screen shot of sagittal view
			ChangeImageSlice_x(curSlice_x);
			sprintf(ss_fn, "%s%s_vetHeight.jpg", screenshot_root,study_name);
			NIH_OpenGL_Save_Display(display_id, ss_fn);
*/		}
		else if(batchIndex==12)	// BMD using manually identified slice using 3D regions
		{
			char ss_fn[400];
			char screenshot_root[200];
			sprintf(screenshot_root, "%s\\screenshot\\", study_root);

			// spine segmentation, partition, labeling, plus rib detection and ROI determination
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  1..");
			fclose(fp);			
			pushButton_pre_process_clicked();
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  2..");
			fclose(fp);			
			pushButton_spinal_cord_clicked();
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  3..");
			fclose(fp);	

			flagShowBBox = true;
			checkBox_bounding_box->setChecked(flagShowBBox);

			// create the direction if not exist
			QDir dcmdir = QString(screenshot_root);
			if(!dcmdir.exists())
			{
				dcmdir.mkdir(screenshot_root);
			}


			cordModel->SetDrawMode(HIDE_MODE);
			smoothedCord->SetDrawMode(HIDE_MODE);

			cordModel->SetDrawMode(LINE_MODE);
			smoothedCord->SetDrawMode(LINE_MODE);

			NIH_OpenGL_Show_Background(display_id, false);
			// screen shot of sagittal view
			ChangeImageSlice_x(curSlice_x);
			sprintf(ss_fn, "%s%s_sag.jpg", screenshot_root,study_name);
			NIH_OpenGL_Save_Display(display_id, ss_fn);
	 	     
			// screen shot of coronal view
			ChangeImageSlice_y(curSlice_y);
			sprintf(ss_fn, "%s%s_cor.jpg", screenshot_root,study_name);
			NIH_OpenGL_Save_Display(display_id, ss_fn);

			NIH_OpenGL_Show_Background(display_id, true);
			flagShowBBox = false;
			checkBox_bounding_box->setChecked(flagShowBBox);
			pushButton_BMD_ROI_L1L2_3D_clicked();		// for manually selected slices
			// screen shot of each BMD_ROI
			for(int s=0; s<roiSlices.GetSize(); s++)
			{
				ChangeImageSlice(roiSlices[s]);
				sprintf(ss_fn, "%s%s_sl%d.jpg", screenshot_root,study_name,roiSlices[s]+1);
				NIH_OpenGL_Save_Display(display_id, ss_fn);
			}
		}
		else if(batchIndex==13)	// record partition location
		{
			char ss_fn[400];
			char screenshot_root[200];
			sprintf(screenshot_root, "%s\\screenshot\\", study_root);

			// spine segmentation, partition, labeling, plus rib detection and ROI determination
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  1..");
			fclose(fp);			
			pushButton_pre_process_clicked();
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  2..");
			fclose(fp);			
			pushButton_spinal_cord_clicked();
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  3..");
			fclose(fp);			
			pushButton_spine_partition_clicked();
			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  4..");
			fclose(fp);			
			pushButton_Rib_Detection_clicked();

			fp = fopen("c:\\batch_progress.txt","a+");
			fprintf(fp, "  5..\n");
			fclose(fp);	

			// record the partition location
			fp = fopen("c:\\spine_partition.txt", "a+");
			fprintf(fp, "%s,%s,%s,%s,%s,", studyEntries[curStudy].patientName.ascii(), studyEntries[curStudy].studyDate.ascii(),
				studyEntries[curStudy].studyId.ascii(), studyEntries[curStudy].seriesId.ascii(),
				studyEntries[curStudy].localImagePath.ascii());
			int st_p;
			if(segInfo.t12_vertebra>=0) st_p=segInfo.t12_vertebra;
			else st_p=0;
			int vertNum = pedicleValley.GetSize();
			for(int p=st_p+1; p<vertNum && p<st_p+6; p++)
			{
				fprintf(fp, "%d,", pedicleValley[p]);
			}
			fprintf(fp, "\n");
		}
	}

	flagBatchMode = false;

}


void NIH_BoneSegmentation_Dlg::listBox_study_entry_clicked( QListBoxItem * )
{
	UpdateData();
	int sel = listBox_study_entry->currentItem();

	if(curStudy==sel) return;

	curStudy = sel;

	QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );

	LoadStudy(curStudy);

	segInfo.Initialize();
	
	char prefix_fn[400];

	sprintf(prefix_fn, "%s_%s_%s", studyEntries[curStudy].patientName.ascii(),
		studyEntries[curStudy].studyId.ascii(),
		studyEntries[curStudy].seriesId.ascii());

	LoadPaintingFile(studyEntries[curStudy].paintingFileName.ascii());
	
	if(NIH_BoneSegmentation_LoadSegmentation(loadSegPath, studyEntries[curStudy].localImagePath, maskImg3D)==CIS_OK)
	{
		NIH_BoneSegmentation_LocateDetection2D(img3D, maskImg3D, maskStruct, segPara, segInfo, 
			detections2D, numDetections2D, lesions, numLesions, paint,
			flagApplyClassifier, svmCommittee, svmCutoff);
		
		spinBox_detection_2D->setRange(0, numDetections2D-1);
		spinBox_detection_2D->setValue(0);

		NIH_BoneSegmentation_MergeMetasis(img3D, maskImg3D, maskStruct, segPara, segInfo, detections, numDetections, detectionVoxels);
		spinBox_detection->setRange(0, numDetections-1);
		spinBox_detection->setValue(0);

		int totalMet, totalLytic, foundMet, foundLytic, totalDetection, tpDetection, fpDetection, 
			tpLyticDetection, fpLyticDetection;
		int totalMet10, totalLytic10, foundMet10, foundLytic10;

		NIH_BoneSegmentation_MatchDetections(maskImg3D, detections, numDetections, detectionVoxels, lesions, numLesions, paint, 
			totalMet, totalLytic, foundMet, foundLytic, totalDetection, tpDetection, fpDetection,
			tpLyticDetection, fpLyticDetection, totalMet10, totalLytic10, foundMet10, foundLytic10, lesionSizeCutoff);
		
		// reporting
		float sensitivity = 1, sensitivityLytic = 1;
		if(totalMet>0) 
		{
			sensitivity = (float)foundMet/(float)totalMet;
		}
		if(totalLytic>0) 
		{	
			sensitivityLytic = (float)foundLytic/(float)totalLytic;
		}
	
		char info[500];
		sprintf(info, "Total Mets: %d (%d); Found: %d (%d);\nTotal Lytic %d (%d); Found: %d (%d)\nDetection: %d; TP: %d; Lytic: %d; FP: %d\n",
			totalMet, totalMet10, foundMet, foundMet10,
			totalLytic, totalLytic10, foundLytic, foundLytic10, 
			totalDetection, tpDetection, tpLyticDetection, fpDetection);

		textLabel_info->setText(info);
	}

	// compute cortical/spongy bone ratio
/*	if(maskImg3D!=NULL)
	{
		int x, y, z, k;
		int sizex, sizey, sizez;
		short *maskA;
		int countc, counts;

		maskA = maskImg3D->GetArray();
		sizex = maskImg3D->Num_Cols();
		sizey = maskImg3D->Num_Rows();
		sizez = maskImg3D->Num_Levels();

		countc = counts = 0;
		for(z=0, k=0; z<sizez; z++)
		{
			for(y=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k++)
				{
					if(maskA[k]==maskStruct.mask_corticalBone)	countc++;
					else if(maskA[k]==maskStruct.mask_spongyBone)	counts++;

				}
			}
		}

		char info[100];
		sprintf(info, "Cortical: %d; Spongy: %d; ratio: %.2f", countc, counts, (float)countc/(countc+counts));
		textLabel_info->setText(textLabel_info->text()+info);

	}
*/

	ChangeImageSlice(curSlice);
	QApplication::restoreOverrideCursor();

}

void NIH_BoneSegmentation_Dlg::radioButton_view_toggled(bool)
{
	if(img3D==NULL) return;

	if(radioButton_view_axial->state()==QButton::On) 
	{
		viewMode=2;
		slider_slice->setValue(centerz);
		slider_slice->setRange(startSlice, endSlice);
		NIH_OpenGL_Show_Background(display_id, true);
		
		if(seedModel) seedModel->vertex_list[0] = Vec2(centerx, centery);

		ChangeImageSlice(centerz);
	}

	if(radioButton_view_sagittal->state()==QButton::On) 
	{
		viewMode=0;
		slider_slice->setValue(centerx);
		slider_slice->setRange(startSlice_x, endSlice_x);
		NIH_OpenGL_Show_Background(display_id, false);
		ChangeImageSlice_x(centerx);

		if(seedModel && mapModel1) seedModel->vertex_list[0] = Vec2(centery*mapModel1->size.x/img3D->Num_Rows()+mapModel1->loc.x, 
			centerz*mapModel1->size.y/img3D->Num_Levels()+mapModel1->loc.y);
		NIH_OpenGL_Refresh_Display(display_id);

	}

	if(radioButton_view_coronal->state()==QButton::On) 
	{
		viewMode=1;
		slider_slice->setValue(centery);
		slider_slice->setRange(startSlice_y, endSlice_y);
		NIH_OpenGL_Show_Background(display_id, false);
		ChangeImageSlice_y(centery);

		if(seedModel && mapModel1) seedModel->vertex_list[0] = Vec2(centerx*mapModel1->size.x/img3D->Num_Cols()+mapModel1->loc.x, 
			centerz*mapModel1->size.y/img3D->Num_Levels()+mapModel1->loc.y);
		NIH_OpenGL_Refresh_Display(display_id);
	}

}

void NIH_BoneSegmentation_Dlg::radioButton_zoom_toggled(bool)
{
	if(img3D==NULL) return;
	if(radioButton_zoom1->state()==QButton::On) zoomMode=1;
	if(radioButton_zoom2->state()==QButton::On) zoomMode=2;
	if(radioButton_zoom4->state()==QButton::On) zoomMode=4;

	if(viewMode==2) ChangeZoom(zoomMode);
}

void NIH_BoneSegmentation_Dlg::ChangeZoom(int newZoom) 
{
	if(img3D==NULL) return;
	if(viewMode!=2) return;
	
	int zoom_centerx, zoom_centery;
	CIS_Array_Image2D_short *sub_img2D=NULL;

	zoomMode = newZoom;

	int orgx, orgy;

	float displayPixelSize = img3D->Get_Pixel_SizeX()/(double)zoomMode;
	int sizex = img3D->Num_Cols();
	int sizey = img3D->Num_Rows();

	if(zoomMode==1)
	{
		zoom_centerx = sizex/2;
		zoom_centery = sizey/2;
	}
	else 
	{
		zoom_centerx = centerx;
		zoom_centery = centery;
	}

	switch(zoomMode)
	{
	case 2:
//		if(sub_img2D!=NULL) delete sub_img2D;
		orgx = zoom_centerx - sizex/4;
		orgy = zoom_centery - sizey/4;
		sub_img2D = img2D->SubImage(orgx,orgy,sizex/2,sizey/2);
		break;
	case 4:
//		if(sub_img2D!=NULL) delete sub_img2D;
		orgx = zoom_centerx - sizex/8;
		orgy = zoom_centery - sizey/8;
		sub_img2D = img2D->SubImage(orgx,orgy,sizex/4,sizey/4);
		break;
	}

	// shift the seed model
	if(seedModel!=NULL)
	{
		seedModel->mapZoom = zoomMode;
		seedModel->mapCenterx = zoom_centerx;
		seedModel->mapCentery = zoom_centery;
		seedModel->mapSizex = sizex;
		seedModel->mapSizey = sizey;
	}

	if(zoomMode!=1)
	{
		bboxModel->SetDrawMode(HIDE_MODE);
		vertebraBodyModel->SetDrawMode(HIDE_MODE);
		spinousProcessModel->SetDrawMode(HIDE_MODE);
		cmapModel1->SetDrawMode(HIDE_MODE);
	}

	if(zoomMode==1) NIH_OpenGL_Display_Image(display_id, img2D); //HERE is where the border is drawn
	else NIH_OpenGL_Display_Image(display_id, sub_img2D, true, true, true);


	if(sub_img2D) delete sub_img2D;
}


void NIH_BoneSegmentation_Dlg::ChangeImageSlice(int newCurSlice)
{
	if(img3D==NULL) return;
	UpdateData();

	cmapModel1->SetDrawMode(HIDE_MODE);

	NIH_OpenGL_Show_Background(display_id, true);

	projectedCordSag->SetDrawMode(HIDE_MODE);
	leftColumnSag->SetDrawMode(HIDE_MODE);
	rightColumnSag->SetDrawMode(HIDE_MODE);
	pedicleModelSag->SetDrawMode(HIDE_MODE);

	projectedCordCor->SetDrawMode(HIDE_MODE);
	leftColumnCor->SetDrawMode(HIDE_MODE);
	rightColumnCor->SetDrawMode(HIDE_MODE);
	pedicleModelCor->SetDrawMode(HIDE_MODE);

	if(newCurSlice!=-1 && newCurSlice>=startSlice && newCurSlice<=endSlice)
	{
		curSlice=newCurSlice;
		slider_slice->setValue(curSlice);

		if(img2D==NULL) img2D = img3D->GetSliceImage(curSlice-startSlice,2);
		else img3D->GetSliceImage(*img2D, curSlice-startSlice, 2);

		mapModel1->SetDrawMode(HIDE_MODE);

		char st[50];
		float curSliceZ_pos = img3D->GetSlicePosition(curSlice-startSlice);
		sprintf(st, "%.3f", curSliceZ_pos);
		textLabel_cur_location->setText(st);

		sprintf(st, "%d", curSlice);
		lineEdit_cur_slice->setText(st);

		cmapModel1->SetDrawMode(HIDE_MODE);

		

		if(flagShowOverlay && maskImg3D!=NULL && curSlice<=maskImg3D->Num_Levels())
		{
			int x, y, sizex, sizey, sizexy, k, k1, i,j,z;
			short *maskArray, *img2DA;
			sizex = maskImg3D->Num_Cols();
			sizey = maskImg3D->Num_Rows();
			sizexy = sizex*sizey;

			CIS_Array_Image2D_RGB_uchar *colorImage;
			colorImage = new CIS_Array_Image2D_RGB_uchar(sizex, sizey);
			maskArray = maskImg3D->GetArray()+(curSlice-startSlice)*sizex*sizey;
			img2DA = img2D->GetArray();


			if(showWSOverlay)
			{
				for(y=0, k=0; y<sizey; y++)
				{
					for(x=0; x<sizex; x++, k++)
					{
						colorImage->FastSet(x+segInfo.spineBound1[curSlice].x-20, y+segInfo.spineBound1[curSlice].y-5, watershedOverlay->FastGet(x,y));
					}
				}
			}

			else if(showOrganMask)
			{
				for(y=0, k=0; y<sizey; y++)
				{
					for(x=0; x<sizex; x++, k++)
					{
						if(maskArray[k]==MASK_LKIDNEY)
							colorImage->FastSet(x, y, RGBtriple_uchar(255,128,0));
						else if(maskArray[k]==MASK_RKIDNEY)
							colorImage->FastSet(x, y, RGBtriple_uchar(255,255,0));
						else if(maskArray[k]==MASK_LIVER)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,255,255));
						else if(maskArray[k]==MASK_PANCREAS)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,255,128));
						else if(maskArray[k]==MASK_SPLEEN)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,128,255));
						else 
							colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
					}
				}
			}
			
			else
			{
				for(y=0, k=0; y<sizey; y++)
				{
					for(x=0; x<sizex; x++, k++)
					{
						/*					if(segInfo.t12_vertebra>=0 && frontPlane[0].pt.x!=0)
						{
						colorImage->FastSet(x, y, ComputeVertebraLabelColor(x, y, curSlice, maskArray[k]));
						continue;
						}
						*/
						if(maskArray[k]==maskStruct.mask_air)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
					//	else if(maskArray[k]==48)
					//		colorImage->FastSet(x, y, RGBtriple_uchar(0,0,255));
						else if(maskArray[k]==maskStruct.mask_body)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
						/*else if(maskArray[k]==maskStruct.mask_corticalBone)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,0,255));
						else if(maskArray[k]==maskStruct.mask_spongyBone)
							colorImage->FastSet(x, y, RGBtriple_uchar(255,0,255));*/
						else if(maskArray[k]==maskStruct.mask_boneMetasis)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,255,0));
						else if(maskArray[k]==maskStruct.mask_falseDetection)
							//						colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
							colorImage->FastSet(x, y, RGBtriple_uchar(128,0,0));
						else if(maskArray[k]==maskStruct.mask_spinalCord)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,0,255));
						else if(maskArray[k]==maskStruct.mask_otherBone)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
					//	else if(maskArray[k]==maskStruct.mask_rib)
					//		colorImage->FastSet(x, y, RGBtriple_uchar(0,0,255));
						else if(maskArray[k]==maskStruct.mask_vertebra)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,255,255));
		///				else if(maskArray[k]==maskStruct.mask_vertebralDisk)
		///					colorImage->FastSet(x, y, RGBtriple_uchar(0,255,0));
						else 
							colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
					}	// for x
				}	// for y
			}
			// overlap the painting
			if(flagShowPainting)
			{
				for(i=0; i<totalPaintNum; i++)
				{
					for(j=0; j<paint[i].numPaint; j++)
					{
						z = paint[i].sliceArray[j]-1;
						if(z!=curSlice) continue;
						for(k=0; k<paint[i].paint[j].GetSize(); k++)
						{
							x = paint[i].paint[j][k].x;
							y = paint[i].paint[j][k].y;
							k1= y*sizex+x;
							if(maskArray[k1]==maskStruct.mask_boneMetasis)
								colorImage->FastSet(x, y, RGBtriple_uchar(0,255,255));
							else 
							{
								if(paint[i].lesionType==1) colorImage->FastSet(x, y, RGBtriple_uchar(255,255,128));
								else colorImage->FastSet(x, y, RGBtriple_uchar(128,255,255));
							}
						}	// for k
					}	// for j
				}	// for i
			}	// if flagShowPainting
			
			cmapModel1->mapping_mode=1;
			cmapModel1->size = IntVec2(sizex, sizey);
			cmapModel1->loc = Vec2(0, 0);
			cmapModel1->SetDrawMode(SURFACE_MODE);
			cmapModel1->ChangeBitmap(colorImage, cmapModel1->mapping_mode);

			delete colorImage;
		}	// flagShowOverlay


		bboxModel->SetDrawMode(HIDE_MODE);
		vertebraBodyModel->SetDrawMode(HIDE_MODE);
		spinousProcessModel->SetDrawMode(HIDE_MODE);

		if(flagShowBBox)
		{
			if(segInfo.spineBound1.GetSize()>curSlice-startSlice && segInfo.spineBound1[curSlice-startSlice].x>0) 
			{
				bboxModel->SetDrawMode(LINE_MODE);
				bboxModel->loc = (segInfo.spineBound1[curSlice-startSlice]+segInfo.spineBound2[curSlice-startSlice])/2;
				bboxModel->size = (segInfo.spineBound2[curSlice-startSlice]-segInfo.spineBound1[curSlice-startSlice]);
			}

			if(segInfo.sprocessBound1.GetSize()>curSlice-startSlice && segInfo.sprocessBound1[curSlice-startSlice].x>0) 
			{
				spinousProcessModel->SetDrawMode(LINE_MODE);
				spinousProcessModel->loc = (segInfo.sprocessBound1[curSlice-startSlice]+segInfo.sprocessBound2[curSlice-startSlice])/2;
				spinousProcessModel->size = (segInfo.sprocessBound2[curSlice-startSlice]-segInfo.sprocessBound1[curSlice-startSlice]);
			}

			if(segInfo.diskBound1.GetSize()>curSlice-startSlice && segInfo.diskBound1[curSlice-startSlice].x>0) 
			{
				vertebraBodyModel->SetDrawMode(LINE_MODE);
				vertebraBodyModel->loc = (segInfo.diskBound1[curSlice-startSlice]+segInfo.diskBound2[curSlice-startSlice])/2;
				vertebraBodyModel->size = (segInfo.diskBound2[curSlice-startSlice]-segInfo.diskBound1[curSlice-startSlice]);
			}
///			seedModel->SetDrawMode(VERTEX_MODE);
		}
///		else
///		{
///			seedModel->SetDrawMode(HIDE_MODE);
///		}

		if(vertebraTemplateViewMode==1)
		{
			templateModel->SetDrawMode(LINE_MODE);

			// form the pedicle model
			templateModel->vertex_list.SetSize(0);
			templateModel->line_list.SetSize(0);

			// add the disk
			int a, i, z;
			int counta, countp, countdl, countdr, count;
			z = curSlice;
			counta = 0;
			for(a=0; a<diskAngleCount; a++)
			{
				if(vertebraTemplate[z].diskContour[a].x!=-1)
				{
					counta++;
					templateModel->vertex_list.Add(vertebraTemplate[z].diskContour[a]);
				}
			}
			if(counta>1)
			{
				for(a=0; a<counta-1; a++) templateModel->line_list.Add(IntVec2(a,a+1));
//				templateModel->line_list.Add(IntVev2(counta-1,0));
			}

			count = counta;
			// add process
			countp = 0;
			for(i=sprocessCount-1; i>=0; i--)
			{
				if(vertebraTemplate[z].sprocessLeft[i].x!=-1)
				{
					countp++;
					templateModel->vertex_list.Add(vertebraTemplate[z].sprocessLeft[i]);
				}
			}
			if(vertebraTemplate[z].sprocessEnd.x!=-1)
			{
				countp++;
				templateModel->vertex_list.Add(vertebraTemplate[z].sprocessEnd);
			}
			for(i=0; i<sprocessCount; i++)
			{
				if(vertebraTemplate[z].sprocessRight[i].x!=-1)
				{
					countp++;
					templateModel->vertex_list.Add(vertebraTemplate[z].sprocessRight[i]);
				}
			}
			if(countp>1)
			{
				for(i=0; i<countp-1; i++) templateModel->line_list.Add(IntVec2(count+i,count+i+1));
//				templateModel->line_list.Add(IntVev2(count+countp-1,count));
			}

			count += countp;

			// pedicle
			// left
			countdl = 0;
			for(i=pedicleCount-1; i>=0; i--)
			{
				if(vertebraTemplate[z].leftPedicleUp[i].x!=-1)
				{
					countdl++;
					templateModel->vertex_list.Add(vertebraTemplate[z].leftPedicleUp[i]);
				}
			}
			if(vertebraTemplate[z].leftPedicleEnd.x!=-1)
			{
				countdl++;
				templateModel->vertex_list.Add(vertebraTemplate[z].leftPedicleEnd);
			}
			for(i=0; i<pedicleCount; i++)
			{
				if(vertebraTemplate[z].leftPedicleDown[i].x!=-1)
				{
					countdl++;
					templateModel->vertex_list.Add(vertebraTemplate[z].leftPedicleDown[i]);
				}
			}
			if(countdl>2)
			{
				for(i=0; i<countdl-1; i++) templateModel->line_list.Add(IntVec2(count+i,count+i+1));
//				templateModel->line_list.Add(IntVev2(count+countdl-1,count));
			}

			count += countdl;
			// pedicle right
			countdr = 0;
			for(i=pedicleCount-1; i>=0; i--)
			{
				if(vertebraTemplate[z].rightPedicleUp[i].x!=-1)
				{
					countdr++;
					templateModel->vertex_list.Add(vertebraTemplate[z].rightPedicleUp[i]);
				}
			}
			if(vertebraTemplate[z].rightPedicleEnd.x!=-1)
			{
				countdr++;
				templateModel->vertex_list.Add(vertebraTemplate[z].rightPedicleEnd);
			}
			for(i=0; i<pedicleCount; i++)
			{
				if(vertebraTemplate[z].rightPedicleDown[i].x!=-1)
				{
					countdr++;
					templateModel->vertex_list.Add(vertebraTemplate[z].rightPedicleDown[i]);
				}
			}
			if(countdr>2)
			{
				for(i=0; i<countdr-1; i++) templateModel->line_list.Add(IntVec2(count+i,count+i+1));
//				templateModel->line_list.Add(IntVev2(count+countdr-1,count));
			}
		}
		else if(vertebraTemplateViewMode==2)	// inital template
		{
			templateModel->SetDrawMode(LINE_MODE);

			// form the pedicle model
			templateModel->vertex_list.SetSize(0);
			templateModel->line_list.SetSize(0);

			// add the disk
			int a, i, z, angle;
			int counta, countp, countdl, countdr, count;
			Vec2 dir, v;
			z = curSlice;
			counta = 0;
			for(a=0; a<diskAngleCount; a++)
			{
				if(vertebraTemplate[z].diskCenter.x!=-1)
				{
					angle = a*diskAngleInterval;
					dir.x = sin(angle*3.14/180);
					dir.y = cos(angle*3.14/180);
					counta++;
					templateModel->vertex_list.Add(vertebraTemplate[z].diskCenter+dir*vertebraTemplate[z].diskRadius);
				}
			}
			if(counta>1)
			{
				for(a=0; a<counta-1; a++) templateModel->line_list.Add(IntVec2(a,a+1));
//				templateModel->line_list.Add(IntVev2(counta-1,0));
			}

			count = counta;
			// add process (3 points)
			countp = 0;
			if(vertebraTemplate[z].cordCenter.x!=-1)
			{
				v.x = vertebraTemplate[z].cordCenter.x-vertebraTemplate[z].cordRadius.x;
				v.y = vertebraTemplate[z].cordCenter.y+vertebraTemplate[z].cordRadius.y;
				templateModel->vertex_list.Add(v);
				v.x = vertebraTemplate[z].cordCenter.x;
				v.y = vertebraTemplate[z].cordCenter.y+vertebraTemplate[z].cordRadius.y+vertebraTemplate[z].diskRadius*2;
				templateModel->vertex_list.Add(v);
				v.x = vertebraTemplate[z].cordCenter.x+vertebraTemplate[z].cordRadius.x;
				v.y = vertebraTemplate[z].cordCenter.y+vertebraTemplate[z].cordRadius.y;
				templateModel->vertex_list.Add(v);
				countp = 3;
			}
			if(countp>1)
			{
				for(i=0; i<countp-1; i++) templateModel->line_list.Add(IntVec2(count+i,count+i+1));
//				templateModel->line_list.Add(IntVev2(count+countp-1,count));
			}

			count += countp;

			// pedicle
			// left
			countdl = 0;
			if(vertebraTemplate[z].cordCenter.x!=-1)
			{
				v.x = vertebraTemplate[z].cordCenter.x-vertebraTemplate[z].cordRadius.x;
				v.y = vertebraTemplate[z].cordCenter.y-vertebraTemplate[z].cordRadius.y;
				templateModel->vertex_list.Add(v);
				v.x = vertebraTemplate[z].cordCenter.x-vertebraTemplate[z].cordRadius.x-vertebraTemplate[z].diskRadius*2;;
				v.y = vertebraTemplate[z].cordCenter.y;
				templateModel->vertex_list.Add(v);
				v.x = vertebraTemplate[z].cordCenter.x-vertebraTemplate[z].cordRadius.x;
				v.y = vertebraTemplate[z].cordCenter.y+vertebraTemplate[z].cordRadius.y;
				templateModel->vertex_list.Add(v);
				countdl = 3;
			}

			if(countdl>2)
			{
				for(i=0; i<countdl-1; i++) templateModel->line_list.Add(IntVec2(count+i,count+i+1));
//				templateModel->line_list.Add(IntVev2(count+countdl-1,count));
			}

			count += countdl;
			// pedicle right
			countdr = 0;
			if(vertebraTemplate[z].cordCenter.x!=-1)
			{
				v.x = vertebraTemplate[z].cordCenter.x+vertebraTemplate[z].cordRadius.x;
				v.y = vertebraTemplate[z].cordCenter.y-vertebraTemplate[z].cordRadius.y;
				templateModel->vertex_list.Add(v);
				v.x = vertebraTemplate[z].cordCenter.x+vertebraTemplate[z].cordRadius.x+vertebraTemplate[z].diskRadius*2;;
				v.y = vertebraTemplate[z].cordCenter.y;
				templateModel->vertex_list.Add(v);
				v.x = vertebraTemplate[z].cordCenter.x+vertebraTemplate[z].cordRadius.x;
				v.y = vertebraTemplate[z].cordCenter.y+vertebraTemplate[z].cordRadius.y;
				templateModel->vertex_list.Add(v);
				countdr = 3;
			}
			if(countdr>2)
			{
				for(i=0; i<countdr-1; i++) templateModel->line_list.Add(IntVec2(count+i,count+i+1));
//				templateModel->line_list.Add(IntVev2(count+countdr-1,count));
			}

		}
		else
		{
			templateModel->SetDrawMode(HIDE_MODE);
		}

		tex_z->ChangeBitmap(img2D, tex_z->mapping_mode);
		tex_z->loc.z = img3D->GetSlicePosition(curSlice);
		tex_z->ComputeFrame();

		// set the cutting plane for the surface model
		for(int i=0; i<totalALMModels; i++)
		{
			for(int m=0; m<5; m++)
			{
				if(singleALM[i].organSurf[m]) singleALM[i].organSurf[m]->slicePlane.Set(tex_z->loc, tex_z->axis);
			}
		}

		NIH_OpenGL_Refresh_Display(display_id2);

	}	// if new slice

	if(viewMode==2) ChangeZoom(zoomMode); //HERE is where the border is drawn
	flagImageChanged = true;
	return;
}


void NIH_BoneSegmentation_Dlg::ChangeImageSlice_x(int newCurSlice_x)
{
	if(img3D==NULL) return;

	UpdateData();

	if(reformationMode && flagShowPainting)
	{
		bboxModel->SetDrawMode(HIDE_MODE);
		vertebraBodyModel->SetDrawMode(HIDE_MODE);
		spinousProcessModel->SetDrawMode(HIDE_MODE);

		if(flagShowBBox)
		{
			projectedCordSag->SetDrawMode(LINE_MODE);
			leftColumnSag->SetDrawMode(LINE_MODE);
			rightColumnSag->SetDrawMode(LINE_MODE);
			pedicleModelSag->SetDrawMode(LINE_MODE);
//			vertHeightRuler->SetDrawMode(LINE_MODE);
		}
		else 
		{
			projectedCordSag->SetDrawMode(HIDE_MODE);
			leftColumnSag->SetDrawMode(HIDE_MODE);
			rightColumnSag->SetDrawMode(HIDE_MODE);
			pedicleModelSag->SetDrawMode(HIDE_MODE);
//			vertHeightRuler->SetDrawMode(HIDE_MODE);
		}

		projectedCordCor->SetDrawMode(HIDE_MODE);
		leftColumnCor->SetDrawMode(HIDE_MODE);
		rightColumnCor->SetDrawMode(HIDE_MODE);
		pedicleModelCor->SetDrawMode(HIDE_MODE);

		float vsizey, vsizez;
		int isizey, isizez;
		vsizey = img3D->Get_Volume_SizeY();
		vsizez = img3D->Get_Volume_SizeZ();
		if(vsizey>vsizez)
		{
			isizey = 512;
			isizez = (int)((vsizez/vsizey)*512);
			mapModel1->size = IntVec2(isizey, isizez);
			mapModel1->loc = Vec2(256-isizey/2, 256-isizez/2);
		}
		else 
		{
			isizey = (int)((vsizey/vsizez)*512);
			isizez = 512;
			mapModel1->size = IntVec2(isizey, isizez);
			mapModel1->loc = Vec2(256-isizey/2, 256-isizez/2);
		}

		int newHiLimit = slider_contrast_hi->value();
		int newLowLimit = slider_contrast_lo->value();

		if(mapModel1) mapModel1->SetWindowLevel(newHiLimit-newLowLimit, (newHiLimit+newLowLimit)/2);
		mapModel1->mapping_mode=1;
		mapModel1->ChangeBitmap(img2D_sag, mapModel1->mapping_mode);
		mapModel1->SetDrawMode(LINE_MODE);

		// color 
		cmapModel1->SetDrawMode(HIDE_MODE);

		if(flagShowOverlay && img2D_sag_mask!=NULL)
		{
			int x, y, sizex, sizey, k;
			short *maskArray;
			sizex = img2D_sag_mask->Num_Cols();
			sizey = img2D_sag_mask->Num_Rows();

			CIS_Array_Image2D_RGB_uchar *colorImage;
			colorImage = new CIS_Array_Image2D_RGB_uchar(sizex, sizey);
			maskArray = img2D_sag_mask->GetArray();

			for(y=0, k=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k++)
				{
					if(maskArray[k]==maskStruct.mask_air)
						colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
					else if(maskArray[k]==maskStruct.mask_body)
						colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
					else if(maskArray[k]==maskStruct.mask_corticalBone)
						colorImage->FastSet(x, y, RGBtriple_uchar(255,0,0));
					else if(maskArray[k]==maskStruct.mask_spongyBone)
						colorImage->FastSet(x, y, RGBtriple_uchar(128,0,0));
					else if(maskArray[k]==maskStruct.mask_boneMetasis)
						colorImage->FastSet(x, y, RGBtriple_uchar(0,255,0));
					else if(maskArray[k]==maskStruct.mask_falseDetection)
						colorImage->FastSet(x, y, RGBtriple_uchar(128,0,0));
					else if(maskArray[k]==maskStruct.mask_spinalCord)
						colorImage->FastSet(x, y, RGBtriple_uchar(0,0,255));
					else if(maskArray[k]==maskStruct.mask_otherBone)
						colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
					else if(maskArray[k]==maskStruct.mask_rib)
						colorImage->FastSet(x, y, RGBtriple_uchar(0,0,255));
					else if(maskArray[k]==maskStruct.mask_vertebra)
						colorImage->FastSet(x, y, RGBtriple_uchar(0,255,255));
					else if(maskArray[k]==maskStruct.mask_vertebralDisk)
						colorImage->FastSet(x, y, RGBtriple_uchar(0,255,0));
					else 
						colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
				}	// for x
			}	// for y

			cmapModel1->mapping_mode=1;
			cmapModel1->size = mapModel1->size;
			cmapModel1->loc = mapModel1->loc;
			cmapModel1->SetDrawMode(SURFACE_MODE);
			cmapModel1->ChangeBitmap(colorImage, cmapModel1->mapping_mode);

			delete colorImage;
		}	// flagShowOverlay
		NIH_OpenGL_Refresh_Display(display_id);
	}
	else if(newCurSlice_x!=-1 && newCurSlice_x>=startSlice_x && newCurSlice_x<=endSlice_x)
	{
		bboxModel->SetDrawMode(HIDE_MODE);
		vertebraBodyModel->SetDrawMode(HIDE_MODE);
		spinousProcessModel->SetDrawMode(HIDE_MODE);

		projectedCordSag->SetDrawMode(HIDE_MODE);
		leftColumnSag->SetDrawMode(HIDE_MODE);
		rightColumnSag->SetDrawMode(HIDE_MODE);
		pedicleModelSag->SetDrawMode(HIDE_MODE);

		projectedCordCor->SetDrawMode(HIDE_MODE);
		leftColumnCor->SetDrawMode(HIDE_MODE);
		rightColumnCor->SetDrawMode(HIDE_MODE);
		pedicleModelCor->SetDrawMode(HIDE_MODE);

		curSlice_x=newCurSlice_x;

		char st[50];
		sprintf(st, "%d", curSlice_x);
		lineEdit_cur_slice->setText(st);

		if(img2D_x==NULL) img2D_x = img3D->GetSliceImage(curSlice_x-startSlice_x,0, true);
		else img3D->GetSliceImage(*img2D_x, curSlice_x-startSlice_x, 0, true);

		float vsizey, vsizez;
		int isizey, isizez;
		vsizey = img3D->Get_Volume_SizeY();
		vsizez = img3D->Get_Volume_SizeZ();
		if(vsizey>vsizez)
		{
			isizey = 512;
			isizez = (int)((vsizez/vsizey)*512);
			mapModel1->size = IntVec2(isizey, isizez);
			mapModel1->loc = Vec2(256-isizey/2, 256-isizez/2);
		}
		else 
		{
			isizey = (int)((vsizey/vsizez)*512);
			isizez = 512;
			mapModel1->size = IntVec2(isizey, isizez);
			mapModel1->loc = Vec2(256-isizey/2, 256-isizez/2);
		}

		int newHiLimit = slider_contrast_hi->value();
		int newLowLimit = slider_contrast_lo->value();

		if(mapModel1) mapModel1->SetWindowLevel(newHiLimit-newLowLimit, (newHiLimit+newLowLimit)/2);
		mapModel1->mapping_mode=1;
		mapModel1->ChangeBitmap(img2D_x, mapModel1->mapping_mode);
		mapModel1->SetDrawMode(LINE_MODE);


		// show the color map
		cmapModel1->SetDrawMode(HIDE_MODE);

		if(flagShowOverlay && maskImg3D!=NULL)
		{
			int x, y, sizex, sizey, sizexy, sizez, k, k1, i,j,z;
			short *maskArray, *img2DA;
			sizex = maskImg3D->Num_Cols();
			sizey = maskImg3D->Num_Rows();
			sizez = maskImg3D->Num_Levels();
			sizexy = sizex*sizey;

			CIS_Array_Image2D_RGB_uchar *colorImage;
			colorImage = new CIS_Array_Image2D_RGB_uchar(sizey, sizez);
			maskArray = maskImg3D->GetArray();
			img2DA = img2D_x->GetArray();

			if(showOrganMask)
			{
				for(y=0, k=0; y<sizez; y++)
				{
					for(x=0; x<sizey; x++)
					{
						if(maskArray[y*sizexy+x*sizex+curSlice_x]==MASK_LKIDNEY)
							colorImage->FastSet(x, y, RGBtriple_uchar(255,128,0));
						else if(maskArray[y*sizexy+x*sizex+curSlice_x]==MASK_RKIDNEY)
							colorImage->FastSet(x, y, RGBtriple_uchar(255,255,0));
						else if(maskArray[y*sizexy+x*sizex+curSlice_x]==MASK_LIVER)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,255,255));
						else if(maskArray[y*sizexy+x*sizex+curSlice_x]==MASK_PANCREAS)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,255,128));
						else if(maskArray[y*sizexy+x*sizex+curSlice_x]==MASK_SPLEEN)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,128,255));
						else 
							colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
					}
				}
			}
			else 
			{
				for(y=0; y<sizez; y++)
				{
					for(x=0; x<sizey; x++)
					{
						if(segInfo.t12_vertebra>=0 && frontPlane[0].pt.x!=0)
						{
							colorImage->FastSet(x, y, ComputeVertebraLabelColor(curSlice_x, x, y, maskArray[y*sizexy+x*sizex+curSlice_x]));
							continue;
						}

						if(maskArray[y*sizexy+x*sizex+curSlice_x]==maskStruct.mask_air)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
						else if(maskArray[y*sizexy+x*sizex+curSlice_x]==maskStruct.mask_body)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
						else if(maskArray[y*sizexy+x*sizex+curSlice_x]==maskStruct.mask_corticalBone)
							colorImage->FastSet(x, y, RGBtriple_uchar(255,0,0));
						else if(maskArray[y*sizexy+x*sizex+curSlice_x]==maskStruct.mask_spongyBone)
							colorImage->FastSet(x, y, RGBtriple_uchar(128,0,0));
						else if(maskArray[y*sizexy+x*sizex+curSlice_x]==maskStruct.mask_boneMetasis)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,255,0));
						else if(maskArray[y*sizexy+x*sizex+curSlice_x]==maskStruct.mask_falseDetection)
							colorImage->FastSet(x, y, RGBtriple_uchar(128,0,0));
						else if(maskArray[y*sizexy+x*sizex+curSlice_x]==maskStruct.mask_spinalCord)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,0,255));
						else if(maskArray[y*sizexy+x*sizex+curSlice_x]==maskStruct.mask_otherBone)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
						else if(maskArray[y*sizexy+x*sizex+curSlice_x]==maskStruct.mask_rib)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,0,255));
						else if(maskArray[y*sizexy+x*sizex+curSlice_x]==maskStruct.mask_vertebra)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,255,255));
						else if(maskArray[y*sizexy+x*sizex+curSlice_x]==maskStruct.mask_vertebralDisk)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,255,0));
						else 
							colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
					}	// for x
				}	// for y
			}

			// overlap the painting
			if(flagShowPainting)
			{
				for(i=0; i<totalPaintNum; i++)
				{
					for(j=0; j<paint[i].numPaint; j++)
					{
						z = paint[i].sliceArray[j]-1;
//						if(z!=curSlice) continue;
						for(k=0; k<paint[i].paint[j].GetSize(); k++)
						{
							x = paint[i].paint[j][k].x;
							y = paint[i].paint[j][k].y;
							if(x==curSlice_x)
							{
								k1= y*sizex+x;
								if(maskArray[z*sizexy+k1]==maskStruct.mask_boneMetasis)
									colorImage->FastSet(y, z, RGBtriple_uchar(0,255,255));
								else 
								{
									if(paint[i].lesionType==1) colorImage->FastSet(y, z, RGBtriple_uchar(255,255,128));
									else colorImage->FastSet(y, z, RGBtriple_uchar(128,255,255));
								}
							}
						}	// for k
					}	// for j
				}	// for i
			}	// if flagShowPainting
			
			cmapModel1->mapping_mode=1;
			cmapModel1->size = mapModel1->size;
			cmapModel1->loc = mapModel1->loc;
			cmapModel1->SetDrawMode(SURFACE_MODE);
			cmapModel1->ChangeBitmap(colorImage, cmapModel1->mapping_mode);

			delete colorImage;
		}	// flagShowOverlay


		if(flagShowBBox)
		{
			seedModel->SetDrawMode(LINE_MODE);
		}
		else
		{
			seedModel->SetDrawMode(HIDE_MODE);
		}

		tex_x->ChangeBitmap(img2D_x, tex_x->mapping_mode);
		tex_x->loc.x = (double)(curSlice_x)*img3D->Get_Pixel_SizeX();
		tex_x->ComputeFrame();
		tex_x->SetWindowLevel(newHiLimit-newLowLimit, (newHiLimit+newLowLimit)/2);	

		// set the cutting plane for the surface model
		for(int i=0; i<totalALMModels; i++)
		{
			for(int m=0; m<5; m++)
			{
				if(singleALM[i].organSurf[m]) singleALM[i].organSurf[m]->slicePlane.Set(tex_x->loc, tex_x->axis);
			}
		}

		NIH_OpenGL_Refresh_Display(display_id);
		NIH_OpenGL_Refresh_Display(display_id2);
	}	// if newCurSlice

	return;
}


void NIH_BoneSegmentation_Dlg::ChangeImageSlice_y(int newCurSlice_y)
{
	if(img3D==NULL) return;

	UpdateData();

	if(reformationMode && flagShowPainting)
	{
		bboxModel->SetDrawMode(HIDE_MODE);
		vertebraBodyModel->SetDrawMode(HIDE_MODE);
		spinousProcessModel->SetDrawMode(HIDE_MODE);

		projectedCordSag->SetDrawMode(HIDE_MODE);
		leftColumnSag->SetDrawMode(HIDE_MODE);
		rightColumnSag->SetDrawMode(HIDE_MODE);
		pedicleModelSag->SetDrawMode(HIDE_MODE);

		if(flagShowBBox)
		{
			projectedCordCor->SetDrawMode(LINE_MODE);
			leftColumnCor->SetDrawMode(LINE_MODE);
			rightColumnCor->SetDrawMode(LINE_MODE);
			pedicleModelCor->SetDrawMode(LINE_MODE);
		}

		float vsizex, vsizez;
		int isizex, isizez;
		vsizex = img3D->Get_Volume_SizeX();
		vsizez = img3D->Get_Volume_SizeZ();
		if(vsizex>vsizez)
		{
			isizex = 512;
			isizez = (int)((vsizez/vsizex)*512);
			mapModel1->size = IntVec2(isizex, isizez);
			mapModel1->loc = Vec2(256-isizex/2, 256-isizez/2);
		}
		else 
		{
			isizex = (int)((vsizex/vsizez)*512);
			isizez = 512;
			mapModel1->size = IntVec2(isizex, isizez);
			mapModel1->loc = Vec2(256-isizex/2, 256-isizez/2);
		}

		int newHiLimit = slider_contrast_hi->value();
		int newLowLimit = slider_contrast_lo->value();

		if(mapModel1) mapModel1->SetWindowLevel(newHiLimit-newLowLimit, (newHiLimit+newLowLimit)/2);
		mapModel1->mapping_mode=1;
		mapModel1->ChangeBitmap(img2D_cor, mapModel1->mapping_mode);
		mapModel1->SetDrawMode(LINE_MODE);

		// color 
		cmapModel1->SetDrawMode(HIDE_MODE);

		if(flagShowOverlay && img2D_cor_mask!=NULL)
		{
			int x, y, sizex, sizey, k;
			short *maskArray;
			sizex = img2D_cor_mask->Num_Cols();
			sizey = img2D_cor_mask->Num_Rows();

			CIS_Array_Image2D_RGB_uchar *colorImage;
			colorImage = new CIS_Array_Image2D_RGB_uchar(sizex, sizey);
			maskArray = img2D_cor_mask->GetArray();

			for(y=0, k=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k++)
				{
					if(maskArray[k]==maskStruct.mask_air)
						colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
					else if(maskArray[k]==maskStruct.mask_body)
						colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
					else if(maskArray[k]==maskStruct.mask_corticalBone)
						colorImage->FastSet(x, y, RGBtriple_uchar(255,0,0));
					else if(maskArray[k]==maskStruct.mask_spongyBone)
						colorImage->FastSet(x, y, RGBtriple_uchar(128,0,0));
					else if(maskArray[k]==maskStruct.mask_boneMetasis)
						colorImage->FastSet(x, y, RGBtriple_uchar(0,255,0));
					else if(maskArray[k]==maskStruct.mask_falseDetection)
						colorImage->FastSet(x, y, RGBtriple_uchar(0,0,180));
					else if(maskArray[k]==maskStruct.mask_spinalCord)
						colorImage->FastSet(x, y, RGBtriple_uchar(0,0,255));
					else if(maskArray[k]==maskStruct.mask_otherBone)
						colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
					else if(maskArray[k]==maskStruct.mask_rib)
						colorImage->FastSet(x, y, RGBtriple_uchar(0,0,255));
					else if(maskArray[k]==maskStruct.mask_vertebra)
						colorImage->FastSet(x, y, RGBtriple_uchar(0,255,255));
					else if(maskArray[k]==maskStruct.mask_vertebralDisk)
						colorImage->FastSet(x, y, RGBtriple_uchar(0,255,0));
					else 
						colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
				}	// for x
			}	// for y

			cmapModel1->mapping_mode=1;
			cmapModel1->size = mapModel1->size;
			cmapModel1->loc = mapModel1->loc;
			cmapModel1->SetDrawMode(SURFACE_MODE);
			cmapModel1->ChangeBitmap(colorImage, cmapModel1->mapping_mode);

			delete colorImage;
		}	// flagShowOverlay
		NIH_OpenGL_Refresh_Display(display_id);
	}
	else if(newCurSlice_y!=-1 && newCurSlice_y>=startSlice_y && newCurSlice_y<=endSlice_y)
	{
		bboxModel->SetDrawMode(HIDE_MODE);
		vertebraBodyModel->SetDrawMode(HIDE_MODE);
		spinousProcessModel->SetDrawMode(HIDE_MODE);

		projectedCordSag->SetDrawMode(HIDE_MODE);
		leftColumnSag->SetDrawMode(HIDE_MODE);
		rightColumnSag->SetDrawMode(HIDE_MODE);
		pedicleModelSag->SetDrawMode(HIDE_MODE);

		projectedCordCor->SetDrawMode(HIDE_MODE);
		leftColumnCor->SetDrawMode(HIDE_MODE);
		rightColumnCor->SetDrawMode(HIDE_MODE);
		pedicleModelCor->SetDrawMode(HIDE_MODE);

		curSlice_y=newCurSlice_y;

		char st[50];
		sprintf(st, "%d", curSlice_y);
		lineEdit_cur_slice->setText(st);

		if(img2D_y==NULL) img2D_y = img3D->GetSliceImage(curSlice_y-startSlice_y,1, true);
		else img3D->GetSliceImage(*img2D_y, curSlice_y-startSlice_y, 1, true);

		mapModel1->SetDrawMode(LINE_MODE);
		float vsizex, vsizez;
		int isizex, isizez;
		vsizex = img3D->Get_Volume_SizeX();
		vsizez = img3D->Get_Volume_SizeZ();
		if(vsizex>vsizez)
		{
			isizex = 512;
			isizez = (int)((vsizez/vsizex)*512);
			mapModel1->size = IntVec2(isizex, isizez);
			mapModel1->loc = Vec2(256-isizex/2, 256-isizex/2);
		}
		else 
		{
			isizex = (int)((vsizex/vsizez)*512);
			isizez = 512;
			mapModel1->size = IntVec2(isizex, isizez);
			mapModel1->loc = Vec2(256-isizex/2, 256-isizez/2);
		}

		int newHiLimit = slider_contrast_hi->value();
		int newLowLimit = slider_contrast_lo->value();

		if(mapModel1) mapModel1->SetWindowLevel(newHiLimit-newLowLimit, (newHiLimit+newLowLimit)/2);
		mapModel1->mapping_mode=1;
		mapModel1->ChangeBitmap(img2D_y, mapModel1->mapping_mode);
		mapModel1->SetDrawMode(LINE_MODE);

		// show the color map
		cmapModel1->SetDrawMode(HIDE_MODE);

		if(flagShowOverlay && maskImg3D!=NULL)
		{
			int x, y, sizex, sizey, sizexy, sizez, k, k1, i,j,z;
			short *maskArray, *img2DA;
			sizex = maskImg3D->Num_Cols();
			sizey = maskImg3D->Num_Rows();
			sizez = maskImg3D->Num_Levels();
			sizexy = sizex*sizey;

			CIS_Array_Image2D_RGB_uchar *colorImage;
			colorImage = new CIS_Array_Image2D_RGB_uchar(sizex, sizez);
			maskArray = maskImg3D->GetArray();
			img2DA = img2D_y->GetArray();

			if(showOrganMask)
			{
				for(y=0, k=0; y<sizez; y++)
				{
					for(x=0; x<sizex; x++)
					{
						if(maskArray[y*sizexy+curSlice_y*sizex+x]==MASK_LKIDNEY)
							colorImage->FastSet(x, y, RGBtriple_uchar(255,128,0));
						else if(maskArray[y*sizexy+curSlice_y*sizex+x]==MASK_RKIDNEY)
							colorImage->FastSet(x, y, RGBtriple_uchar(255,255,0));
						else if(maskArray[y*sizexy+curSlice_y*sizex+x]==MASK_LIVER)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,255,255));
						else if(maskArray[y*sizexy+curSlice_y*sizex+x]==MASK_PANCREAS)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,255,128));
						else if(maskArray[y*sizexy+curSlice_y*sizex+x]==MASK_SPLEEN)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,128,255));
						else 
							colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
					}
				}
			}
			else 
			{
				for(y=0; y<sizez; y++)
				{
					for(x=0; x<sizex; x++)
					{
						if(segInfo.t12_vertebra>=0 && frontPlane[0].pt.x!=0)
						{
							colorImage->FastSet(x, y, ComputeVertebraLabelColor(x, curSlice_y, y, maskArray[y*sizexy+curSlice_y*sizex+x]));
							continue;
						}

						if(maskArray[y*sizexy+curSlice_y*sizex+x]==maskStruct.mask_air)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
						else if(maskArray[y*sizexy+curSlice_y*sizex+x]==maskStruct.mask_body)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
						else if(maskArray[y*sizexy+curSlice_y*sizex+x]==maskStruct.mask_corticalBone)
							colorImage->FastSet(x, y, RGBtriple_uchar(255,0,0));
						else if(maskArray[y*sizexy+curSlice_y*sizex+x]==maskStruct.mask_spongyBone)
							colorImage->FastSet(x, y, RGBtriple_uchar(128,0,0));
						else if(maskArray[y*sizexy+curSlice_y*sizex+x]==maskStruct.mask_boneMetasis)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,255,0));
						else if(maskArray[y*sizexy+curSlice_y*sizex+x]==maskStruct.mask_falseDetection)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,0,180));
						else if(maskArray[y*sizexy+curSlice_y*sizex+x]==maskStruct.mask_spinalCord)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,0,255));
						else if(maskArray[y*sizexy+curSlice_y*sizex+x]==maskStruct.mask_otherBone)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
						else if(maskArray[y*sizexy+curSlice_y*sizex+x]==maskStruct.mask_rib)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,0,255));
						else if(maskArray[y*sizexy+curSlice_y*sizex+x]==maskStruct.mask_vertebra)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,255,255));
						else if(maskArray[y*sizexy+curSlice_y*sizex+x]==maskStruct.mask_vertebralDisk)
							colorImage->FastSet(x, y, RGBtriple_uchar(0,255,0));
						else 
							colorImage->FastSet(x, y, RGBtriple_uchar(0,0,0));
					}	// for x
				}	// for y
			}

			// overlap the painting
			if(flagShowPainting)
			{
				for(i=0; i<totalPaintNum; i++)
				{
					for(j=0; j<paint[i].numPaint; j++)
					{
						z = paint[i].sliceArray[j]-1;
///						if(z!=curSlice) continue;
						for(k=0; k<paint[i].paint[j].GetSize(); k++)
						{
							x = paint[i].paint[j][k].x;
							y = paint[i].paint[j][k].y;
							if(y==curSlice_y)
							{
								k1= y*sizex+x;
								if(maskArray[z*sizexy+k1]==maskStruct.mask_boneMetasis)
									colorImage->FastSet(x, z, RGBtriple_uchar(0,255,255));
								else 
								{
									if(paint[i].lesionType==1) colorImage->FastSet(x, z, RGBtriple_uchar(255,255,128));
									else colorImage->FastSet(x, z, RGBtriple_uchar(128,255,255));
								}
							}	// if y
						}	// for k
					}	// for j
				}	// for i
			}	// if flagShowPainting
			
			cmapModel1->mapping_mode=1;
			cmapModel1->size = mapModel1->size;
			cmapModel1->loc = mapModel1->loc;
			cmapModel1->SetDrawMode(SURFACE_MODE);
			cmapModel1->ChangeBitmap(colorImage, cmapModel1->mapping_mode);

			delete colorImage;
		}	// flagShowOverlay

		if(flagShowBBox)
		{
			seedModel->SetDrawMode(LINE_MODE);
		}
		else
		{
			seedModel->SetDrawMode(HIDE_MODE);
		}

		tex_y->ChangeBitmap(img2D_y, tex_y->mapping_mode);
		tex_y->loc.y = (double)(curSlice_y)*img3D->Get_Pixel_SizeY();
		tex_y->ComputeFrame();
		tex_y->SetWindowLevel(newHiLimit-newLowLimit, (newHiLimit+newLowLimit)/2);	

		// set the cutting plane for the surface model
		for(int i=0; i<totalALMModels; i++)
		{
			for(int m=0; m<5; m++)
			{
				if(singleALM[i].organSurf[m]) singleALM[i].organSurf[m]->slicePlane.Set(tex_y->loc, tex_y->axis);
			}
		}

		NIH_OpenGL_Refresh_Display(display_id);
		NIH_OpenGL_Refresh_Display(display_id2);
	}	// if newCurSlice

	return;
}

RGBtriple_uchar NIH_BoneSegmentation_Dlg::ComputeVertebraLabelColor(int x, int y, int z, short mask)
{
	RGBtriple_uchar rgb;

	rgb = RGBtriple_uchar(0,0,0);

	// the color label should be pre-computed
	// could be part of the mask, using the strategy of lung
	if(mask==maskStruct.mask_corticalBone || mask==maskStruct.mask_spongyBone || mask==maskStruct.mask_vertebralDisk)
	{
		float px, py;
		px = img3D->Get_Pixel_SizeX();
		py = img3D->Get_Pixel_SizeY();

		Vec3 v;
		v.x = px*x;
		v.y = py*y;
		v.z = img3D->GetSlicePosition(z);

		int p=0;
		double d1, d2;

		if(v.y>spineCenter[z].y) d1=backPlane[p].SignDist(v);
		else d1 = frontPlane[p].SignDist(v);

		if(d1>0)
		{
			for(p=1; p<pedicleValley.GetSize() && p<30; p++)
			{
				if(v.y>spineCenter[z].y) d2=backPlane[p].SignDist(v);
				else d2=frontPlane[p].SignDist(v);

				if(d1*d2<=0) break;
				d1=d2;
			}
		}

		if(mask==maskStruct.mask_vertebralDisk && fabs(d1)<5) return rgb;

		int ci = p-segInfo.t12_vertebra+4;

		intVec3DynArray colorTable;
		colorTable.SetSize(10);
		colorTable[0] = IntVec3(25,255,225);
		colorTable[1] = IntVec3(150,255,100);
		colorTable[2] = IntVec3(75,255,175);
		colorTable[3] = IntVec3(200,255,50);
		colorTable[4] = IntVec3(0,255,0);
		colorTable[5] = IntVec3(255,0,0);
		colorTable[6] = IntVec3(255,200,50);
		colorTable[7] = IntVec3(255,75,175);
		colorTable[8] = IntVec3(255,150,100);
		colorTable[9] = IntVec3(255,25,225);

		if(ci<0) ci=0;
		if(ci>9) ci=9;

		rgb = RGBtriple_uchar(colorTable[ci].x, colorTable[ci].y, colorTable[ci].z);
	}

	return rgb;
}

int NIH_BoneSegmentation_Dlg::OpenDataSet(char *fnDir, bool flagDcm)
{
	if(img3D!=NULL) 
	{
		delete img3D;
		img3D = NULL;
	}

	segmentationChanged = false;

	img3D = new CIS_Array_Image3D_short();

	QDir dcmdir = QString(fnDir)+"dcm\\";
	if(dcmdir.exists())
	{
		strcat(fnDir, "dcm\\");
	}

	// get study_root and study_name
	strcpy(study_root, fnDir);
	int i;
	// remove the dcm, go up the directory
	for(i=strlen(study_root)-2; i>=0; i--)
	{
		if(study_root[i]=='\\' || study_root[i]=='/') break;
	}
	study_root[i+1]=0;
	// get the study name
	for(i=strlen(study_root)-2; i>=0; i--)
	{
		if(study_root[i]=='\\' || study_root[i]=='/') break;
	}
	strcpy(study_name, study_root+(i+1));
	if(study_name[strlen(study_name)-1]=='\\' || study_name[strlen(study_name)-1]=='/') study_name[strlen(study_name)-1]=0;

	
	if(flagDcm==false) 
	{
		// read interfile
		// PET image format
		if(QString(fnDir).find(".hdr")!=-1) 
		{
			int bytes_per_pixel;
			if(Read_Interfile_File(fnDir, *img3D, bytes_per_pixel)==CIS_ERROR && bytes_per_pixel==4)
			{
				// float data type
				CIS_Array_Image3D_float *imgf;
				imgf = new CIS_Array_Image3D_float();
				if(Read_Interfile_File(fnDir, *imgf, bytes_per_pixel)==CIS_OK)
				{
					int sizex, sizey, sizez, sizexyz;
					float *arrayf;
					short *array;
					sizex = imgf->Num_Cols();
					sizey = imgf->Num_Rows();
					sizez = imgf->Num_Levels();
					sizexyz = sizex*sizey*sizez;

					img3D->SetSize(sizex, sizey, sizez);
					arrayf = imgf->GetArray();
					array = img3D->GetArray();
					for(int k=0; k<sizexyz; k++) array[k] = (short)arrayf[k];
					delete imgf;
				};
			};
		}
		else if(QString(fnDir).find(".nii")!=-1 || QString(fnDir).find(".gz")!=-1)
		{
			// nifti image
			Read_NIFTI_File(fnDir, *img3D);
		}
		else img3D->Load(fnDir);
	}
	else
	{
		// dicom images
		img3D->LoadDicomDirectory(fnDir,0);
		Dicom_Get_Dicom_Filename(fnDir, base_fn);
		char tmpfn[200];
		strcpy(tmpfn, fnDir);
		if(tmpfn[strlen(tmpfn)-1]!='\\' && tmpfn[strlen(tmpfn)-1]!='/') strcat(tmpfn, "\\");
		strcat(tmpfn, base_fn);
		strcpy(base_fn, tmpfn);
	}

	
	if(img3D->Num_Levels()<=0) 
	{
		delete img3D;
		img3D = NULL;

		return CIS_ERROR;
	}

	if(img3D!=NULL && img3D->Num_Levels()>0)
	{
		startSlice = 0;
		endSlice = startSlice+img3D->Num_Levels()-1;
		curSlice = (startSlice+endSlice)/2;

		slider_slice->setRange(startSlice, endSlice);
		slider_slice->setValue(curSlice);

		startSlice_x = 0;
		endSlice_x = startSlice_x+img3D->Num_Cols()-1;
		curSlice_x = (startSlice_x+endSlice_x)/2;

		startSlice_y = 0;
		endSlice_y = startSlice_y+img3D->Num_Rows()-1;
		curSlice_y = (startSlice_y+endSlice_y)/2;


		centerx = curSlice_x;
		centery = curSlice_y;
		centerz = curSlice;

		segInfo.bound1 = IntVec3(0,0,0);
		segInfo.bound2 = IntVec3(img3D->Num_Cols()-1, img3D->Num_Rows()-1, img3D->Num_Levels()-1);
		segInfo.spineBound1.SetSize(0);
		segInfo.spineBound2.SetSize(0);
	}

	// allocate buffers for mask
	// also allocate buffers for other segmentation
	if(img3D!=NULL)
	{
		int sizex = img3D->Num_Cols();
		int sizey = img3D->Num_Rows();
		int sizez = img3D->Num_Levels();

		if(maskImg3D!=NULL) delete maskImg3D;

		maskImg3D = new CIS_Array_Image3D_short(*img3D);
		short *maskA=maskImg3D->GetArray();
		int sizexyz=sizex*sizey*sizez;
		for(int k=0; k<sizexyz; k++) maskA[k]=0;
		
		vertebraTemplateViewMode = 0;
		if(vertebraTemplate) free(vertebraTemplate);
		vertebraTemplate = (VertebraStruct2D *)calloc(sizez, sizeof(VertebraStruct2D));

		if(img2D_x==NULL) img2D_x = new CIS_Array_Image2D_short(sizey, sizez);
		else img2D_x->SetSize(sizey, sizez);
		if(img2D_sag==NULL) img2D_sag = new CIS_Array_Image2D_short(sizey, sizez);
		else img2D_sag->SetSize(sizey, sizez);
		if(img2D_sag_mask==NULL) img2D_sag_mask = new CIS_Array_Image2D_short(sizey, sizez);
		else img2D_sag_mask->SetSize(sizey, sizez);

		if(img2D_y==NULL) img2D_y = new CIS_Array_Image2D_short(sizex, sizez);
		else img2D_y->SetSize(sizex, sizez);
		if(img2D_cor==NULL) img2D_cor = new CIS_Array_Image2D_short(sizex, sizez);
		else img2D_cor->SetSize(sizex, sizez);
		if(img2D_cor_mask==NULL) img2D_cor_mask = new CIS_Array_Image2D_short(sizex, sizez);
		else img2D_cor_mask->SetSize(sizex, sizez);

		UpdateData(FALSE);
	}

	if(strcmp(base_fn, "") && flagDcm)
	{
		DicomHeader *dh;
		dh = new DicomHeader();
		if(dh->ReadHeaderFromFile(base_fn)!=CIS_ERROR)
		{
			char *buf, ibuf[20], info[500], tmpinfo[200];
			int length;

			// readin the base value of hounsfield number
			buf = dh->GetItem(0x0028, 0x1052, length);
			if(buf!=NULL)
			{
				strncpy(ibuf, buf, length);
				ibuf[length] = 0;
				int tmpValue;
				sscanf(ibuf, "%d", &tmpValue);
				baseValue = (short)tmpValue;
			}

			// display some info
			strcpy(info,"");
			// series number
			buf = dh->GetItem(0x0008, 0x103e, length);
			if(buf!=NULL) 
			{
				sprintf(tmpinfo, "Series: %s\n", buf);
				strcat(info, tmpinfo);
			}
			// study date
			buf = dh->GetItem(0x0008, 0x0020, length);
			if(buf!=NULL) 
			{
				sprintf(tmpinfo, "Study Date: %s\n", buf);
				strcat(info, tmpinfo);
			}

			// sex
			buf = dh->GetItem(0x0010, 0x0040, length);
			if(buf!=NULL) 
			{
				patientSex = buf;
			}

			// age
			buf = dh->GetItem(0x0010, 0x1010, length);
			if(buf!=NULL) 
			{
				patientAge = buf;
			}

			// image size
			short tx, ty;
			buf = dh->GetItem(0x0028, 0x0010, length);
			memcpy((void *)&tx, buf, length);
			buf = dh->GetItem(0x0028, 0x0011, length);
			memcpy((void *)&ty, buf, length);
			sprintf(tmpinfo, "%d, %d\n", tx, ty);
			strcat(info, tmpinfo);

			textLabel_info->setText(info);

			// check the image orientation
			// if prone, then rotate the image 180 degree
			buf = dh->GetItem(0x0018, 0x5100, length);
			QString pos=QString(buf).simplifyWhiteSpace();
			if(pos=="FFP") 
			{
				int sizex = img3D->Num_Cols();
				int sizey = img3D->Num_Rows();
				int sizez = img3D->Num_Levels();
				int sizexy=sizex*sizey;
				int x, y, z, k, k1;
				short *imgA=img3D->GetArray();
				short tg;
				for(z=0; z<sizez; z++)
				{
					k1= z*sizexy;
					for(y=0,k=0; y<sizey/2; y++)
					{
						for(x=0; x<sizex; x++, k++)
						{
							tg = imgA[k+k1]; imgA[k+k1]=imgA[sizexy-1-k+k1]; imgA[sizexy-1-k+k1]=tg;
						}
					}
				}
	
			}

		}
		delete dh;
	}

	// change color map slider
	if(img3D!=NULL && img3D->Num_Levels()>0)
	{
		int grayHi, grayLow;

		grayHi = 2000;
		grayLow = 1;

		// abdomin window/level
		grayHi = 1264;
		grayLow = 904;
		
		slider_contrast_hi->setRange(grayLow,4000);
		slider_contrast_lo->setRange(grayLow,4000);
		slider_contrast_hi->setValue(grayHi);
		slider_contrast_lo->setValue(grayLow);
	
///		NIH_OpenGL_Set_Background_WindowLevel(display_id, 
///			grayHi-grayLow, (grayHi+grayLow)/2);

		// ron's window/level setting
		NIH_OpenGL_Set_Background_WindowLevel(display_id, 
			2500, 1504);
		NIH_OpenGL_Set_Background_WindowLevel(display_id, 
			360, 60+1024);

///		lineEdit_constrast_hi->setText(QString::number(grayHi));
///		lineEdit_constrast_lo->setText(QString::number(grayLow));
		
	}

	if(img3D!=NULL && img3D->Num_Levels()>0)
	{
		viewMode = 2;
		radioButton_view_axial->setChecked(true);
		ChangeImageSlice(curSlice);
	}

	// Set 3D slicer view
	if(1)
	{
		float volx, voly, volz;
		float cenx, ceny, cenz;
		volx = img3D->Get_Volume_SizeX();
		voly = img3D->Get_Volume_SizeY();
		volz = img3D->Get_Volume_SizeZ();
		cenx = volx/2;
		ceny = voly/2;
		cenz = img3D->GetSlicePosition(img3D->Num_Levels()/2);
		tex_x->SetDimension(Vec3(cenx, ceny, cenz), Vec3(-1,0,0), Vec3(0,1,0), voly, volz);
		tex_y->SetDimension(Vec3(cenx, ceny, cenz), Vec3(0,1,0), Vec3(1,0,0), volx, volz);
		tex_z->SetDimension(Vec3(cenx, ceny, cenz), Vec3(0,0,1), Vec3(1,0,0), volx, voly);
///		tex_x->SetDrawMode(LINE_MODE);
///		tex_y->SetDrawMode(LINE_MODE);
///		tex_z->SetDrawMode(LINE_MODE);
		tex_x->SetWindowLevel(360, 60+1024);
		tex_y->SetWindowLevel(360, 60+1024);
		tex_z->SetWindowLevel(360, 60+1024);
	}
/*
	// Setup volume rendering
	if(volumeModel==NULL)
	{
///		volumeModel = new CIS_3D_Model_Volume(img3D,1,-3);

		int xx0=10, xx1=500, yy0=200, yy1=380, zz0=0, zz1=270;
		int sizex, sizey, sizez;
		float px, py, pz;
		sizex = img3D->Num_Cols();
		sizey = img3D->Num_Rows();
		sizez = img3D->Num_Levels();
		px = img3D->Get_Pixel_SizeX();
		py = img3D->Get_Pixel_SizeY();
		pz = img3D->Get_Pixel_SizeZ();


		volumeModel = new CIS_3D_Model_Volume();
		CIS_Array_Image2D_short *imgT_2D;
		imgT_2D = new CIS_Array_Image2D_short(sizex, sizez);
		Vec3 sliceLoc;
		double dimx, dimy;
		Vec3 normal, dirX;
		dimx = sizex*px;
		dimy = sizez*pz;
		normal = Vec3(0,-1,0);
		dirX = Vec3(1,0,0);
		sliceLoc.x = px*sizex/2;
		sliceLoc.z = img3D->GetSlicePosition(sizez/2);
		volumeModel->SetTotalTextures(yy1-yy0+1);
		int x, z, k;
		short *img2DA;
		for(int y=yy1, i=0; y>=yy0; y--, i++)
		{
			img3D->GetSliceImage(*imgT_2D, y, 1, true);
			img2DA=imgT_2D->GetArray();
			sliceLoc.y=y*py;
			for(z=0, k=0; z<sizez; z++)
			{
				for(x=0; x<sizex; x++, k++)
				{
					if(x<xx0 || x>xx1 || z<zz0 || z>zz1) img2DA[k]=0;
				}
			}
			volumeModel->SetOneTexture(i, imgT_2D, sliceLoc, normal, dirX, dimx, dimy);
		}

		delete imgT_2D;
		volumeModel->SetRenderingParameter(1300, 0.1, 300, 1200);
		CIS_Model_AddModel(volumeModel);
		NIH_OpenGL_AddModel(display_id2, volumeModel);	
		volumeModel->SetDrawMode(LINE_MODE);
		volumeModel->SetCurrAlpha(0.5);
		NIH_OpenGL_CentralizeModel(display_id2, tex_y);
	}
*/
	// test the save template function
//	img3D->SaveDICOMviaTemplate("c:\\tmp\\3DSeries", "NULL");
/*	int newSize_x = 256;
	int newSize_y = 256;
	int newSize_z = img3D->Num_Levels();
	img3D->Resample(newSize_x, newSize_y, newSize_z);
	img3D->SaveDICOM("c:\\tmp\\3DSeries", "NULL");
*/
	return CIS_OK;
}

int NIH_BoneSegmentation_Dlg::LoadStudy(int st)
{
	if(st<0 || st>=numStudyEntries) return -1;

	// open the image
	char imgDir[400];
	sprintf(imgDir, "%s%s\\", patientRoot, studyEntries[st].localImagePath.ascii());
	OpenDataSet(imgDir, true);
	sprintf(fn_path1, "%s%s", patientRoot, studyEntries[st].localImagePath.ascii());
	if(fn_path1[strlen(fn_path1)-1]!='\\' && fn_path1[strlen(fn_path1)-1]!='/') strcat(fn_path1, "\\");

	return 0;
}


int NIH_BoneSegmentation_Dlg::LoadPaintingFile(const char *painting_fn)
{
	if(img3D==NULL) return 0;
	if(!strcmp(painting_fn,"NULL") || strlen(painting_fn)==0) return 0;

	char fn[400];
	sprintf(fn,"%s%s", paintingRoot, painting_fn);
	totalPaintNum = 0;

	if(LoadPaintingFromFile(fn, paint, totalPaintNum, true)==CIS_ERROR) return CIS_ERROR;

	int i, j, k;
	int x, y ,z;

	// fill the mask for visualization
	for(i=0; i<totalPaintNum; i++)
	{
		for(j=0; j<paint[i].numPaint; j++)
		{
			z = paint[i].sliceArray[j]-1;
			if(z==0)
			{
				int ij=0;
			}
			for(k=0; k<paint[i].paint[j].GetSize(); k++)
			{
				x = paint[i].paint[j][k].x;
				y = paint[i].paint[j][k].y;
				if(paint[i].lesionType==1) maskImg3D->FastSet(x, y, z, maskStruct.mask_paintLytic);
				else maskImg3D->FastSet(x, y, z, maskStruct.mask_paint);
			}
		}
	}

	// compute lesion statistics
	if(lesions!=NULL)
	{
		delete lesions;
		lesions=NULL;
		numLesions = 0;
	}

	numLesions = totalPaintNum;
	lesions = new LesionStructure[numLesions];

	IntVec3 lesionLoc;
	int lesionSlices, curArea2D, curArea2D_1200;
	intDynArray lesionArea_1200, lesionArea_2D, lesionArea, lesionArea_1200_2D;

	lesionArea.SetSize(numLesions);
	lesionArea_1200.SetSize(numLesions);
	lesionArea_2D.SetSize(numLesions);
	lesionArea_1200_2D.SetSize(numLesions);

	for(i=0; i<numLesions; i++)
	{
		lesionArea_1200[i] = 0;
		lesionArea_2D[i] = 0;
		lesionArea_1200_2D[i] = 0;
		lesionArea[i] = 0;

		lesions[i].lesionStatus = -1;
		lesions[i].lesionType = paint[i].lesionType;
		lesions[i].lesionSize = paint[i].lesionSize;

		lesionLoc = IntVec3(0,0,0);
		lesionSlices = 0;
		
		for(j=0; j<paint[i].numPaint; j++)
		{
			z = paint[i].sliceArray[j]-1;
			curArea2D = curArea2D_1200 = 0;
			for(k=0; k<paint[i].paint[j].GetSize(); k++)
			{
				x = paint[i].paint[j][k].x;
				y = paint[i].paint[j][k].y;
				lesionLoc += IntVec3(x,y,z);
				lesionArea[i] ++;
				curArea2D ++;

				if(img3D->FastGet(x,y,z)<1200) 
				{
					lesionArea_1200[i]++;
					curArea2D_1200 ++;
				}
			}
			lesionSlices++;

			if(curArea2D>lesionArea_2D[i]) lesionArea_2D[i]=curArea2D;
			if(curArea2D_1200>lesionArea_1200_2D[i]) lesionArea_1200_2D[i]=curArea2D_1200;
		}

		if(lesionArea[i]!=0) 
		{
			lesionLoc.x /= lesionArea[i];
			lesionLoc.y /= lesionArea[i];
			lesionLoc.z /= lesionArea[i];
		}

		lesions[i].lesionSlices = lesionSlices;
		lesions[i].lesionArea = lesionArea[i];
		lesions[i].lesionLoc = lesionLoc;
		lesions[i].lesionKey = paint[i].detectionKey;
	}

	spinBox_lesion->setRange(0, totalPaintNum-1);
	spinBox_lesion->setValue(0);
	
	char info[200];
	// count the lesionType
	int count1, count2, count3, count4;
	count1=count2=count3=count4=0;

	for(i=0; i<totalPaintNum; i++)
	{
		if(paint[i].lesionType==1) count1++;
		else if(paint[i].lesionType==2) count2++;
		else if(paint[i].lesionType==3) count3++;
		else if(paint[i].lesionType>0) count4++;
	}

	sprintf(info, "Total Mets:%d\nLytic: %d; Sclerotic: %d\nMixed: %d; Others: %d", totalPaintNum, count1, count2, count3, count4);
	textLabel_info->setText(textLabel_info->text()+info);
	
	ChangeImageSlice(curSlice);

/*	// report the lesion statistics
	float px = img3D->Get_Pixel_SizeX();
	float py = img3D->Get_Pixel_SizeY();

	FILE *fp;
	if((fp = fopen("d:\\tmp\\bm_painting_stat_test_2005_0725.csv", "a"))!=NULL)
	{
		fprintf(fp,"CTSeriesName,LesionKey,LesionCategory,LesionLocX,LesionLocY,LesionLocZ,LesionSize,LesionArea,LesionArea_1200,LesionArea_2D,LesionArea_1200_2D,NumOfSlices\n");
		for(i=0; i<totalPaintNum; i++)
		{
			fprintf(fp,"%s,%d,%d,%d,%d,%d,%f,%.3f,%.3f,%.3f,%.3f,%d\n",paint[i].seriesName,lesions[i].lesionKey,lesions[i].lesionType,
				lesions[i].lesionLoc.x,lesions[i].lesionLoc.y,lesions[i].lesionLoc.z,
				lesions[i].lesionSize,
				lesionArea[i]*px*py/100, lesionArea_1200[i]*px*py/100, 
				lesionArea_2D[i]*px*py/100, lesionArea_1200_2D[i]*px*py/100, 
				lesions[i].lesionSlices);
		}
		fclose(fp);
	}
*/
/*	FILE *fp;
	if((fp = fopen("d:\\tmp\\bm_data_stat.csv", "a"))!=NULL)
	{
		fprintf(fp,"%s,%d\n",paint[0].seriesName,img3D->Num_Levels());
		fclose(fp);
	} 
*/	return 0;
}


 
QStatusBar *globalStatusBar=NULL;

int main( int argc, char **argv )
{ 
    QApplication::setColorSpec( QApplication::CustomColor ); 
    QApplication *a;
	a = new QApplication( argc, argv );

	// parse the argument to get the project file
	strcpy(global_project_filename, "");
	strcpy(global_image_dir, "");
	strcpy(global_ALM_batch_file, "");
	global_ALM_batch_id=-1;
	global_batch_start=-1;
	global_batch_end=-1;
	for(int ar=1; ar<argc; ar+=2)
	{
		if(!strcmp(argv[ar], "-f"))
		{
			strcpy(global_project_filename, argv[ar+1]);
		}
		else if(!strcmp(argv[ar], "-s"))
		{
			sscanf(argv[ar+1], "%d", &global_batch_start);
		}
		else if(!strcmp(argv[ar], "-e"))
		{
			sscanf(argv[ar+1], "%d", &global_batch_end);
		}
		else if(!strcmp(argv[ar], "-a"))
		{
			strcpy(global_ALM_batch_file, argv[ar+1]);
		}
		else if(!strcmp(argv[ar], "-i"))
		{
			sscanf(argv[ar+1], "%d", &global_ALM_batch_id);
		}
		else if(!strcmp(argv[ar], "-m"))
		{
			global_ALM_multiple_fn = argv[ar+1];
		}
		else if(!strcmp(argv[ar], "-r"))
		{
			global_ALM_multiple_reference = argv[ar+1];
		}
		else if(!strcmp(argv[ar], "-x"))
		{
			global_ALM_multiple_exclude = argv[ar+1];
		}
	}

	if(argc==2) strcpy(global_image_dir, argv[1]);

	QT_GLView *gv=NULL;
	QMainWindow *frame;
	frame = new QMainWindow();

	NIH_OpenGL_Initialize();


	frame->setGeometry(50, 50, 1200, 1100);
	globalStatusBar = frame->statusBar();
	globalStatusBar->message("Clear");
	
    gv = new QT_GLView(frame);
	int display_id = NIH_OpenGL_Init_Display(gv);
	gv->setGeometry(5, 0, 512, 512);
///	gv->setCursor(QCursor(Qt::CrossCursor));

    QT_GLView *gv2 = new QT_GLView(frame);
	int display_id2 = NIH_OpenGL_Init_Display(gv2);
	gv2->setGeometry(5, 522, 512, 512);
	
	NIH_BoneSegmentation_Dlg *dlg = new NIH_BoneSegmentation_Dlg(display_id, display_id2, frame, "Bone Segmentation", Qt::WStyle_DialogBorder);
	dlg->setGeometry(520, 0, 700, 600);
	dlg->show();

//	QDir::setCurrent("c:\\ds\\vc\\nn1");
	QDir::setCurrent("c:\\Experiments\\bone_training\\Atlas");
    a->setMainWidget( frame );
//	frame.showMaximized();
	frame->show();

	if((global_batch_start>=0 && global_batch_end>=0) || strlen(global_image_dir)>0 || global_ALM_batch_id>=0
		|| global_ALM_multiple_fn.length()>2)
		QTimer::singleShot( 100, dlg, SLOT(CommandLine_clicked()));

/*	QKeyEvent keyE(QEvent::KeyPress, Qt::Key_B, 'B', Qt::NoButton);
	QMouseEvent mouseE(QEvent::MouseButtonPress, QPoint(0,0), Qt::LeftButton , Qt::NoButton);
	a->sendEvent(a, &keyE);
	a->sendEvent(a, &mouseE);
	a->sendEvent(dlg, &keyE);
	a->sendEvent(dlg, &mouseE);
*/
	a->exec();

	return 0;
}



void NIH_BoneSegmentation_Dlg::useLiveWire_stateChanged( int )
{
	flagUseLiveWire = useLiveWire->isChecked();
}

void NIH_BoneSegmentation_Dlg::checkBox_record_point_stateChanged( int )
{
	flagRecordPoint = checkBox_record_point->isChecked();
}

void NIH_BoneSegmentation_Dlg::checkBox_ls_show_init_model_stateChanged( int )
{
	bool show = checkBox_ls_show_init_model->isChecked();
	if(show) liveWireInitModel->draw_mode = LINE_MODE;
	else liveWireInitModel->draw_mode = HIDE_MODE;
	NIH_OpenGL_Refresh_Display(display_id);
}

void NIH_BoneSegmentation_Dlg::StartLiveWire(IntVec2 pt)
{
	if(img2D==NULL  || img2D->Num_Cols()==0) return; 

	if(flagImageChanged) 
	{
		livewire.SetImage(img2D);
		flagImageChanged = false;
	}
	
	// set parameter
	livewire.useTraining = useTraining->isChecked();
	livewire.useLevelSet = checkBox_lw_use_levelset->isChecked();
	livewire.useDirectionalSearch = checkBox_lw_directional_search->isChecked();
	livewire.searchRange = lineEdit_lw_search_range->text().toInt();
	livewire.coolLength = lineEdit_lw_cool_length->text().toInt();
	livewire.levelSetWeight = lineEdit_lw_levelset_weight->text().toFloat();
	livewire.trainingWeight = lineEdit_lw_training_weight->text().toFloat();
	livewire.trainingBin = lineEdit_lw_training_bin->text().toFloat();
	livewire.levelSetSmooth = checkBox_ls_smooth->isChecked();

	// get levelset para
	livewire.ls_para.para_diffusionConductance = lineEdit_ls_diffusion_conductance->text().toDouble();
	livewire.ls_para.para_diffusionTimeStep = lineEdit_ls_diffusion_timestep->text().toDouble();
	livewire.ls_para.para_diffusionIteration = lineEdit_ls_diffusion_iteration->text().toInt();
	livewire.ls_para.para_lap_advectionScaling = lineEdit_ls_laplacian_advection->text().toDouble();
	livewire.ls_para.para_lap_curvatureScaling = lineEdit_ls_laplacian_curvature->text().toDouble();
	livewire.ls_para.para_lap_propagationScaling = lineEdit_ls_laplacian_propagation->text().toDouble();
	livewire.ls_para.para_lap_maxError = lineEdit_ls_laplacian_max_error->text().toDouble();
	livewire.ls_para.para_lap_iteration = lineEdit_ls_laplacian_iteration->text().toInt();
	livewire.ls_para.para_gradient_sigma = lineEdit_ls_gradient_sigma->text().toDouble();
	livewire.ls_para.para_layer = lineEdit_ls_layer->text().toInt();

	if(testGradientMagCost_org->isChecked()==true) livewire.costFunction = NIH_LiveWire::COST_ORG_GRADIENT;
	else if(testGradientMagCost->isChecked()==true) livewire.costFunction = NIH_LiveWire::COST_NORMALIZED_GRADIENT;

	livewire.SetRoot(pt.x, pt.y);
	liveWireModel->vertex_list.SetSize(0);
}

void NIH_BoneSegmentation_Dlg::AddLiveWireSeed(IntVec2 pt)
{
	if(img2D==NULL || img2D->Num_Cols()==0) return; 

	livewire.AddSeed(pt.x, pt.y);
	// display path
	vec2DynArray fpath, apath;
	livewire.GetFrozenPath(fpath);
	livewire.GetActivePath(apath);
	liveWireModel->vertex_list = fpath; 
	liveWireModel->vertex_list.Append(apath);

	if(checkBox_show_active_contour->isChecked()) 
	{
		liveWireModel->color_list.SetSize(fpath.GetSize()+apath.GetSize());
		int i;
		for(i=0; i<fpath.GetSize(); i++) liveWireModel->color_list[i]=0;
		for(; i<fpath.GetSize()+apath.GetSize(); i++) liveWireModel->color_list[i]=1;
	}
	NIH_OpenGL_Refresh_Display(display_id);

}

void NIH_BoneSegmentation_Dlg::FreezeLiveWireSeed(IntVec2 pt)
{
	if(img2D==NULL) return; 

	livewire.FreezePath(pt.x, pt.y);
	// display path
	vec2DynArray fpath, apath;
	livewire.GetFrozenPath(fpath);
	livewire.GetActivePath(apath);
	liveWireModel->vertex_list = fpath; 
	liveWireModel->vertex_list.Append(apath);

	if(checkBox_show_active_contour->isChecked()) 
	{
		liveWireModel->color_list.SetSize(fpath.GetSize()+apath.GetSize());
		int i;
		for(i=0; i<fpath.GetSize(); i++) liveWireModel->color_list[i]=0;
		for(; i<fpath.GetSize()+apath.GetSize(); i++) liveWireModel->color_list[i]=1;
	}
	NIH_OpenGL_Refresh_Display(display_id);

}

int ComputeLevelSetUsingITK(CIS_Array_Image2D_short *img, 
							NIH_LevelSet_Parameter &ls_para,
						   vec2DynArray &initContour, bool closed,
						   float *levelSetDistanceArray,
						   vec2DynArray *zeroLevelSet, bool flagSmooth=false, int levelSetFunc=0);

//read points from file for livewire simulation
void NIH_BoneSegmentation_Dlg::pushButton_read_point_clicked()
{
	if(img2D==NULL) return; 

	if(checkBox_ls_show_init_model->isChecked()) liveWireInitModel->draw_mode= LINE_MODE;
	else liveWireInitModel->draw_mode= HIDE_MODE;

	QString strNewFileName;
	if(!flagBatchMode)
	{
		strNewFileName = QFileDialog::getOpenFileName(
			QString::null, "Text File (*.txt)",
	         this, "Record Points",
		     "Choose a file" );
	} else
	{
		strNewFileName = lw_recordedFileName;
	}
 
	flagLiveWireAutoMode = true;

	int levelSetFunc=0;
	if(radioButton_ls_func_laplacian->isChecked()) levelSetFunc=0;
	else if(radioButton_ls_func_geodesic->isChecked()) levelSetFunc=1;

	if(strNewFileName!=QString::null) 
	{
        QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );
		FILE *fp;
 		char buf[2000];
		QStringList line;
		recordedPoints.SetSize(0);
		bool first = true;
		int count=0;
		if((fp=fopen(strNewFileName.ascii(), "r"))!=NULL)
		{
			QTime t1;

			t1.start();
			if(checkBox_ls_just_levelset->isChecked())
			{
				vec2DynArray zeroLevelSet;
				NIH_LevelSet_Parameter ls_para;
				float *levelSetDistanceArray;
				bool levelSetSmooth;
				
				if(checkBox_ls_smooth->isChecked()) levelSetSmooth=true;
				else levelSetSmooth=false;

				int sizexy = img2D->Num_Cols()*img2D->Num_Rows();
				levelSetDistanceArray = new float[sizexy];
				for(int k=0; k<sizexy; k++) levelSetDistanceArray[k]=-sizexy;

				// get levelset para
				ls_para.para_diffusionConductance = lineEdit_ls_diffusion_conductance->text().toDouble();
				ls_para.para_diffusionTimeStep = lineEdit_ls_diffusion_timestep->text().toDouble();
				ls_para.para_diffusionIteration = lineEdit_ls_diffusion_iteration->text().toInt();
				ls_para.para_lap_advectionScaling = lineEdit_ls_laplacian_advection->text().toDouble();
				ls_para.para_lap_curvatureScaling = lineEdit_ls_laplacian_curvature->text().toDouble();
				ls_para.para_lap_propagationScaling = lineEdit_ls_laplacian_propagation->text().toDouble();
				ls_para.para_lap_maxError = lineEdit_ls_laplacian_max_error->text().toDouble();
				ls_para.para_lap_iteration = lineEdit_ls_laplacian_iteration->text().toInt();
				ls_para.para_layer = lineEdit_ls_layer->text().toInt();

				ls_para.para_gradient_sigma = lineEdit_ls_gradient_sigma->text().toDouble();

				liveWireInitModel->vertex_list.SetSize(0);
				while(!feof(fp))
				{
					strcpy(buf, "");
					fgets(buf, 2000, fp);
					if(strlen(buf)<2) continue;
					line = QStringList::split(',', buf, true);
					Vec2 pt;
					pt.x=line[0].toDouble();;
					pt.y=line[1].toDouble();;
					liveWireInitModel->vertex_list.Add(pt);
				}

				ComputeLevelSetUsingITK(img2D, ls_para, liveWireInitModel->vertex_list, true, levelSetDistanceArray, &zeroLevelSet, 
					levelSetSmooth, levelSetFunc);
				delete levelSetDistanceArray;
				liveWireModel->vertex_list = zeroLevelSet; 
				liveWireModel->draw_mode = VERTEX_MODE;
				NIH_OpenGL_Refresh_Display(display_id);
			} else
			{
				liveWireInitModel->vertex_list.SetSize(0);
				while(!feof(fp))
				{
					strcpy(buf, "");
					fgets(buf, 2000, fp);
					if(strlen(buf)<2) continue;
					line = QStringList::split(',', buf, true);
					IntVec2 pt;
					pt.x=line[0].toInt();;
					pt.y=line[1].toInt();;
					recordedPoints.Add(pt);

					if(first) StartLiveWire(pt);
					else AddLiveWireSeed(pt);

					liveWireInitModel->vertex_list.Add(pt);
					first = false;
					count++;
				}

				// evaluate the smoothness of initial contour
/*				Vec2 ti, ti_1;
				lw_smoothness = 0;
				ti = (liveWireInitModel->vertex_list[1]-liveWireInitModel->vertex_list[0]).normalize();
	
				for(int i=2; i<recordedPoints.GetSize()-2; i++)
				{
					ti_1 = (liveWireInitModel->vertex_list[i]-liveWireInitModel->vertex_list[i-1]).normalize();

					lw_smoothness += (ti*ti_1);

					ti = ti_1;
				}

				lw_smoothness /= (float)(recordedPoints.GetSize()-4);
				char info[200];
				sprintf(info,"Sm:%.3f\n", lw_smoothness);
				textLabel_info->setText(info);
*/			}
			fclose(fp);

			lw_elapsed_time = (float)t1.elapsed()/(float)1000;
			ComputeError(liveWireModel->vertex_list, groundTruth, groundTruthMask, lw_avgError, lw_maxError, lw_stdError,
				lw_sensitivity, lw_specificity, lw_dice, lw_hausdorff, lw_smoothness);
			char info[200];
			sprintf(info,"Time:%f\nAvg:%.3f, Max:%.3f, Std:%.3f, \nSen:%.3f, Spe:%.3f, Di:%.3f, Ha:%.3f\nSm:%.3f\n", lw_elapsed_time, 
				lw_avgError, lw_maxError, lw_stdError, lw_sensitivity, lw_specificity, lw_dice, lw_hausdorff, lw_smoothness);
			textLabel_info->setText(info);
		}
        QApplication::restoreOverrideCursor();
	}
	flagLiveWireAutoMode = false;
	return;
}

//write points to file for livewire simulation
void NIH_BoneSegmentation_Dlg::pushButton_write_point_clicked()
{
	QString strNewFileName = QFileDialog::getSaveFileName(
		QString::null, "Text File (*.txt)",
         this, "Record Points",
         "Choose a file" );
 
	if(strNewFileName!=QString::null) 
	{
		FILE *fp;
	
		if((fp=fopen(strNewFileName.ascii(),"w"))!=NULL)
		{
			for(int i=0; i<recordedPoints.GetSize(); i++)
			{
				fprintf(fp,"%d,%d\n", recordedPoints[i].x, recordedPoints[i].y);
			}
			fclose(fp);
		}
	}
}

void NIH_BoneSegmentation_Dlg::testImg_clicked()
{
	if(img2D==NULL) img2D = new CIS_Array_Image2D_short(512, 512);
	if(groundTruthMask==NULL) groundTruthMask = new CIS_Array_Image2D_short(512, 512);


	short bg = 50;
	short contrast=20;
	short fg = 100;
	short fg2 = 200;
	short fg4 = 200;
	short noise=20;
	float blur=3;

	short c1_x=400, c1_y=400, c1_r=50;
	short c2_x=435, c2_y=400, c2_r=25;
	short c3_x=295, c3_y=400, c3_r=50;
	short c4_x=100, c4_y=100, c4_r=50;
	short sq1_x=400, sq1_y=400, sq1_r=50;

	int c1_r2, c2_r2, c3_r2, c4_r2;
	c1_r2 = c1_r*c1_r;
	c2_r2 = c2_r*c2_r;
	c3_r2 = c3_r*c3_r;
	c4_r2 = c4_r*c4_r;

	bg = lineEdit_lw_background->text().toShort();
	contrast = lineEdit_lw_contrast->text().toShort();
	noise = lineEdit_lw_noise->text().toShort();
	blur = lineEdit_lw_blur->text().toShort();

	fg = bg +contrast;
	fg2 = bg +contrast*2;
	fg4 = bg +contrast*10;

	int x, y, k;
	int sizex, sizey;
	short *array;
	sizex = img2D->Num_Cols();
	sizey = img2D->Num_Rows();
	array = img2D->GetArray();

	// fill the background
	for(y=0, k=0; y<sizey; y++)
	{
		for(x=0; x<sizex; x++, k++) array[k]=bg;
	}

	// paint a square
	if(extra_square->isChecked())
	{
		for(y=sq1_y-sq1_r; y<=sq1_y+sq1_r; y++)
		{
			for(x=sq1_x-sq1_r; x<=sq1_x+sq1_r; x++)
			{
				array[y*sizex+x]=fg;
			}
		}
	}
	else {
		// paint first circle
		for(y=0, k=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++) 
			{
				if(((y-c1_y)*(y-c1_y)+(x-c1_x)*(x-c1_x))<c1_r2)
				{
					array[k]=fg;
				}
			}
		}

		// make a gap and a bump
		if(extra_circle->isChecked())
		{		
			int gap_x0, gap_x1, gap_y0, gap_y1;
			int gap_size;

			// bump
			gap_size =2;
			gap_x0 = c1_x+c1_r; gap_x1 = c1_x+c1_r+c1_r*1/8;
			gap_y0 = c1_y-gap_size; gap_y1 = c1_y+gap_size;
			for(y=gap_y0; y<=gap_y1; y++)
			{
				for(x=gap_x0; x<=gap_x1; x++)
				{
					array[y*sizex+x]=fg+10;
				}
			}

			// gap
			gap_size =2;
			gap_x0 = c1_x-c1_r; gap_x1 = c1_x-c1_r+c1_r*1/8;
			gap_y0 = c1_y-gap_size; gap_y1 = c1_y+gap_size;
			for(y=gap_y0; y<=gap_y1; y++)
			{
				for(x=gap_x0; x<=gap_x1; x++)
				{
					array[y*sizex+x]=bg;
				}
			}
		}	
	}

	// paint extra circle
/*	if(extra_circle->isChecked())
	{
	   for(y=0, k=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++) 
			{
				if(((y-c2_y)*(y-c2_y)+(x-c2_x)*(x-c2_x))<c2_r2)
				{
					array[k]=fg;
				}
			}
	   }	
	}
*/
	// paint third circle
	if(extra_circle_2->isChecked())
	{
		for(y=0, k=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++) 
			{
				if(((y-c3_y)*(y-c3_y)+(x-c3_x)*(x-c3_x))<c3_r2)
				{
					array[k]=fg2;
				}
			}
		}
	}

	// paint interference circle
    for(y=0, k=0; y<sizey; y++)
	{
		for(x=0; x<sizex; x++, k++) 
		{
			if(((y-c4_y)*(y-c4_y)+(x-c4_x)*(x-c4_x))<c4_r2)
			{
				array[k]=fg4;
			}
		}
	}

	if(noise>0) CIS_IPA_AddGaussianNoise(img2D, noise, 1, (short)0, (short)0);


	if(blur>0)
	{
		typedef itk::Image< float, 2 >   ImageType;
		ImageType::Pointer ITKim = ImageType::New();
		Copy_CISImage_to_ITKImage_2D(img2D,ITKim);
		
		float para_gauss_variance = blur;
		int para_gauss_kernel_width = 20;

		typedef itk::DiscreteGaussianImageFilter<itkImage_2D_float,itkImage_2D_float> GaussFilterType;
		GaussFilterType::Pointer gauss =  GaussFilterType::New();
		gauss->SetInput(ITKim);
		gauss->SetVariance( para_gauss_variance );
		gauss->SetMaximumKernelWidth( para_gauss_kernel_width);
	
		gauss->Update();

		Copy_ITKImage_to_CISImage_2D(gauss->GetOutput(), img2D);
	}

	float CNR = (float)contrast/(float)noise;
	lineEdit_lw_CNR->setText(QString::number(CNR));

	// get ground truth
	short *groundTruthA = groundTruthMask->GetArray();

	groundTruth.SetSize(0);
	for(y=0, k=0; y<sizey; y++)
	{
		for(x=0; x<sizex; x++, k++) groundTruthA[k]=0;
	}
 	if(!extra_square->isChecked())
	{
		for(y=0, k=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++) 
			{
				double dist2 = sqrt(double((y-c1_y)*(y-c1_y)+(x-c1_x)*(x-c1_x)));
				if(dist2>=c1_r-0.5 && dist2<=c1_r+0.5)
				{
					groundTruth.Add(Vec2(x,y));
				}

				if(dist2<=c1_r) groundTruthA[k]=1;
			}
		}
	} else
	{
		for(y=0, k=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++) 
			{
				if((abs(x-sq1_x)<=sq1_r+0.5 && abs(x-sq1_x)>=sq1_r-0.5 && abs(y-sq1_y)<=sq1_r+0.5) ||
					(abs(y-sq1_y)<=sq1_r+0.5 && abs(y-sq1_y)>=sq1_r-0.5 && abs(x-sq1_x)<=sq1_r+0.5))
				{
					groundTruth.Add(Vec2(x,y));
				}

				if(abs(x-sq1_x)<=sq1_r && abs(y-sq1_y)<=sq1_r) groundTruthA[k]=1;
			}
		}
	}


	QApplication::restoreOverrideCursor();
	
	flagImageChanged = true;
	NIH_OpenGL_Display_Image(display_id, img2D, true, true, true);
	slider_contrast_lo_valueChanged(0); 

}

void NIH_BoneSegmentation_Dlg::ComputeError(vec2DynArray &testContour, vec2DynArray &refContour, CIS_Array_Image2D_short *refMask, 
											float &avgError, float &maxError, float &stdError, 
											float &sensitivity, float &specificity, float &dice,
											float &hausdorff, float &smoothness)
{
	doubleDynArray errorA;
	if(refMask==NULL || refMask->Num_Cols()==0) return;
	int tsize = testContour.GetSize();
	int rsize = refContour.GetSize();
	if(tsize<=0 || rsize<=0) return;

	int i, j;
	double dist2, min_dist2;

	errorA.SetSize(tsize);

	for(i=0; i<tsize; i++)
	{
		min_dist2 = 100000;
		for(j=0; j<rsize; j++)
		{
			dist2 = (testContour[i]-refContour[j])*(testContour[i]-refContour[j]);
			if(dist2<min_dist2) min_dist2=dist2;
		}

		errorA[i] = sqrt(min_dist2);
	}

	// compute statistics
	avgError = 0;
	maxError = 0;
	for(i=0; i<tsize; i++)
	{
		avgError += errorA[i];
		if(errorA[i]>maxError) maxError=errorA[i];
	}

	avgError /= (float)tsize;

	stdError = 0;
	for(i=0; i<tsize; i++)
	{
		stdError += (errorA[i]-avgError)*(errorA[i]-avgError);
	}

	stdError /= (tsize-1);
	stdError = sqrt(stdError);

	// compute the Hausdorff distance
	float hab, hba;
	hab = 0;
	for(i=0; i<tsize; i++)
	{
		min_dist2 = 100000;
		for(j=0; j<rsize; j++)
		{
			dist2 = (testContour[i]-refContour[j])*(testContour[i]-refContour[j]);
			if(dist2<min_dist2) min_dist2=dist2;
		}
		if(hab<min_dist2) hab=min_dist2;
	}
	hab = sqrt((float)hab);

	hba = 0;
	for(i=0; i<rsize; i++)
	{
		min_dist2 = 100000;
		for(j=0; j<tsize; j++)
		{
			dist2 = (testContour[j]-refContour[i])*(testContour[j]-refContour[i]);
			if(dist2<min_dist2) min_dist2=dist2;
		}
		if(hba<min_dist2) hba=min_dist2;
	}
	hba = sqrt((float)hba);

	if(hab>hba) hausdorff=hab; 
	else hausdorff=hba;


	// Compute the sensitivity, specificity and DSC
	int NTP, NFN, NTN, NFP;

	// generate the region inside the segemented contour
	intVec2DynArray segmentRegion;
	CIS_Array_Image2D_short *segmentedMask;
	int sizex, sizey, sizexy, x, y;
	short *segMaskA, *refMaskA;
	sizex = refMask->Num_Cols();
	sizey = refMask->Num_Rows();
	sizexy = sizex*sizey;
	segmentedMask = new CIS_Array_Image2D_short(sizex, sizey);
	segMaskA = segmentedMask->GetArray();
	refMaskA = refMask->GetArray();
	CIS_Algo_Contour_GetScanWhole(testContour, segmentRegion);
	for(i=0; i<sizexy; i++) segMaskA[i]=0;
	
	for(i=0; i<segmentRegion.GetSize(); i++)
	{
		x = segmentRegion[i].x; y = segmentRegion[i].y;
		if(x>=0 && x<sizex && y>0 && y<sizey) segMaskA[x+y*sizex]=1; 
	}

	NTP = NFN = NTN = NFP = 0;

	for(i=0, y=0; y<sizey; y++)
	{
		for(x=0; x<sizex; x++, i++)
		{
			if(segMaskA[i]==1 && refMaskA[i]==1) NTP++;
			if(segMaskA[i]==0 && refMaskA[i]==0 && x>sizex/2 && y>sizey/2) NTN++;
			if(segMaskA[i]==0 && refMaskA[i]==1) NFN++;
			if(segMaskA[i]==1 && refMaskA[i]==0) NFP++;
		}
	}

	sensitivity = (float)(NTP)/(float)(NTP+NFN);
	specificity = (float)(NTN)/(float)(NTN+NFP);
	dice = (float)(2*NTP)/(float)(2*NTP+NFP+NFN);

	Vec2 ti, ti_1;
	smoothness = 0;
	ti = (testContour[1]-testContour[0]).normalize();

	for(i=2; i<tsize-2; i++)
	{
		ti_1 = (testContour[i]-testContour[i-1]).normalize();

		smoothness += (ti*ti_1);

		ti = ti_1;
	}

	smoothness /= (float)(tsize-4);

	delete segmentedMask;
	return;
}

void NIH_BoneSegmentation_Dlg::pushButton_smooth_contour_clicked()
{
	// Laplacian smoothing
	for(int i=1; i<liveWireModel->vertex_list.GetSize()-1; i++)
	{
		liveWireModel->vertex_list[i] = (liveWireModel->vertex_list[i-1]+liveWireModel->vertex_list[i+1])/2;
	}
}

void NIH_BoneSegmentation_Dlg::pushButton_lw_batch_clicked()
{
	UpdateData();

	QString filter = "Batch File (*.txt)";
  
    QFileDialog* cf = new QFileDialog( this );
    cf->addFilter( filter);
 	
	if(cf->exec()==QDialog::Accepted) 
	{
		QString strNewFileName;
	 	char proj_fn[200], line[2000];
		QString value_s, key_s;
		QStringList qline;
	
		strNewFileName=cf->selectedFile();

		strcpy(proj_fn,strNewFileName.ascii());

		FILE *fp;

		if((fp=fopen(proj_fn, "r"))!=NULL)
		{
			flagBatchMode = true;
			while(!feof(fp))
			{
				strcpy(line,"");
				fgets(line, 2000, fp);
				if(strlen(line)<2) continue;
				if(line[0]=='-' && line[1]=='1') break;

				if(strlen(line)>2 && line[0]!='#' && line[0]!=' ')
				{
					qline = QStringList::split(',', line, false);
					if(qline.size()!=2) continue;
					key_s = qline[0].simplifyWhiteSpace().lower();
					value_s = qline[1].simplifyWhiteSpace();
					
					if(key_s=="report")
					{
						lw_reportFileName = value_s;
						// print some header
						FILE *fp;
						if((fp=fopen(lw_reportFileName.ascii(), "a"))!=NULL)
						{
							fprintf(fp, "Time: %s\n", QDateTime::currentDateTime().toString("yyyy-MM-dd:hh:mm:ss").ascii());
						}
						fclose(fp);
					}
					if(key_s=="summary_file")
					{
						lw_summaryFileName = value_s;
					}
					if(key_s=="summary")
					{
						// print summary line
						FILE *fp;
						if((fp=fopen(lw_summaryFileName.ascii(), "a"))!=NULL)
						{
							fprintf(fp, "\n\nSummary: %s\n", value_s.ascii());
						}
						fclose(fp);
					}
					else if(key_s=="use_levelset")
					{
						if(value_s=="1") checkBox_lw_use_levelset->setChecked(true);
						else checkBox_lw_use_levelset->setChecked(false);
					}
					else if(key_s=="use_training")
					{
						if(value_s=="1") useTraining->setChecked(true);
						else useTraining->setChecked(false);
					}
					else if(key_s=="directional_search")
					{
						if(value_s=="1") checkBox_lw_directional_search->setChecked(true);
						else checkBox_lw_directional_search->setChecked(false);
					}
					else if(key_s=="cost_function")
					{
						if(value_s=="normalized") 
						{
							testGradientMagCost->setChecked(true);
							testGradientMagCost_org->setChecked(false);
						}
						else 
						{
							testGradientMagCost_org->setChecked(true);
							testGradientMagCost->setChecked(false);
						}
					}
					else if(key_s=="search_range")
					{
						lineEdit_lw_search_range->setText(value_s);
					}
					else if(key_s=="cool_length")
					{
						lineEdit_lw_cool_length->setText(value_s);
					}
					else if(key_s=="levelset_weight")
					{
						lineEdit_lw_levelset_weight->setText(value_s);
					}
					else if(key_s=="training_weight")
					{
						lineEdit_lw_training_weight->setText(value_s);
					}
					else if(key_s=="extra_circle")
					{
						if(value_s=="1") extra_circle->setChecked(true);
						else extra_circle->setChecked(false);
					}
					else if(key_s=="extra_circle2")
					{
						if(value_s=="1") extra_circle_2->setChecked(true);
						else extra_circle_2->setChecked(false);
					}
					else if(key_s=="extra_square")
					{
						if(value_s=="1") extra_square->setChecked(true);
						else extra_square->setChecked(false);
					}
					else if(key_s=="background")
					{
						lineEdit_lw_background->setText(value_s);
					}
					else if(key_s=="contrast")
					{
						lineEdit_lw_contrast->setText(value_s);
					}
					else if(key_s=="noise")
					{
						lineEdit_lw_noise->setText(value_s);
					}
					else if(key_s=="blur")
					{
						lineEdit_lw_blur->setText(value_s);
					}
					else if(key_s=="just_levelset")
					{
						if(value_s=="1") checkBox_ls_just_levelset->setChecked(true);
						else checkBox_ls_just_levelset->setChecked(false);
					}
					else if(key_s=="levelset_smooth")
					{
						if(value_s=="1") checkBox_ls_smooth->setChecked(true);
						else checkBox_ls_smooth->setChecked(false);
					}
					else if(key_s=="diffusion_iteration")
					{
						lineEdit_ls_diffusion_iteration->setText(value_s);
					}
					else if(key_s=="diffusion_conductance")
					{
						lineEdit_ls_diffusion_conductance->setText(value_s);
					}
					else if(key_s=="diffusion_timestep")
					{
						lineEdit_ls_diffusion_timestep->setText(value_s);
					}
					else if(key_s=="gradient_sigma")
					{
						lineEdit_ls_gradient_sigma->setText(value_s);
					}
					else if(key_s=="laplacian_iteration")
					{
						lineEdit_ls_laplacian_iteration->setText(value_s);
					}
					else if(key_s=="laplacian_curvature")
					{
						lineEdit_ls_laplacian_curvature->setText(value_s);
					}
					else if(key_s=="laplacian_propagation")
					{
						lineEdit_ls_laplacian_propagation->setText(value_s);
					}
					else if(key_s=="laplacian_max_error")
					{
						lineEdit_ls_laplacian_max_error->setText(value_s);
					}
					else if(key_s=="levelset_layer")
					{
						lineEdit_ls_layer->setText(value_s);
					}
					else if(key_s=="recorded_point")
					{
						lw_recordedFileName = value_s;

						testImg_clicked();
						livewire.SetImage(img2D);
						flagImageChanged = false;
						pushButton_read_point_clicked();
						ReportResults();
					}

				}
			}	// while feof

			flagBatchMode = false;
		}	// if fp
	}	// if cf

	return;
}

void NIH_BoneSegmentation_Dlg::ReportResults()
{
	FILE *fp;
	if((fp=fopen(lw_reportFileName.ascii(), "a"))!=NULL)
	{
		// print out the parameters
		// test image
		fprintf(fp, "BG: %s, Contrast: %s, Noise: %s, CNR:%s, Blur: %s\n", 
			lineEdit_lw_background->text().ascii(),
			lineEdit_lw_contrast->text().ascii(),
			lineEdit_lw_noise->text().ascii(),
			lineEdit_lw_CNR->text().ascii(),
			lineEdit_lw_blur->text().ascii());
		fprintf(fp, "excircle: %d, excircle2: %d, exsquare: %d\n", 
			extra_circle->isChecked(), extra_circle_2->isChecked(),
			extra_square->isChecked());
		// livewire parameters
		fprintf(fp, "LevelSet Cost: %d, Training: %d, DirectSearch: %d, LocalCost:%d\n", 
			checkBox_lw_use_levelset->isChecked(), useTraining->isChecked(),
			checkBox_lw_directional_search->isChecked(), testGradientMagCost->isChecked());
		fprintf(fp, "SearchRange: %s, CoolLength: %s, LS_W: %s, TR_W: %s\n", 
			lineEdit_lw_search_range->text().ascii(),
			lineEdit_lw_cool_length->text().ascii(),
			lineEdit_lw_levelset_weight->text().ascii(),
			lineEdit_lw_training_weight->text().ascii());

		// levelset parameters
		fprintf(fp, "Just LevelSet: %d, LS Smoothing: %d\n", 
			checkBox_ls_just_levelset->isChecked(), 
			checkBox_ls_smooth->isChecked());
		fprintf(fp, "Diffusion Iteration: %s, Conductance: %s, Timestep: %s, sigma: %s\n", 
			lineEdit_ls_diffusion_iteration->text().ascii(),
			lineEdit_ls_diffusion_conductance->text().ascii(),
			lineEdit_ls_diffusion_timestep->text().ascii(),
			lineEdit_ls_gradient_sigma->text().ascii());
		fprintf(fp, "Laplacian Iteration: %s, Curvature: %s, Propagation: %s, Error: %s\n", 
			lineEdit_ls_laplacian_iteration->text().ascii(),
			lineEdit_ls_laplacian_curvature->text().ascii(),
			lineEdit_ls_laplacian_propagation->text().ascii(),
			lineEdit_ls_laplacian_max_error->text().ascii());


		// recorded points
		fprintf(fp, "Recorded point file: %s\n", lw_recordedFileName.ascii());

		// results
		fprintf(fp, "Time:%f, Avg Error:%.3f, Max Error:%.3f, Std Error:%.3f, Sen:%.3f, Spe:%.3f, Di:%.3f, Ha:%.3f\nSm:%.3f\n\n\n", 
			lw_elapsed_time, lw_avgError, lw_maxError, lw_stdError, lw_sensitivity, lw_specificity, lw_dice, lw_hausdorff, lw_smoothness);

		fclose(fp);
	}

	if((fp=fopen(lw_summaryFileName.ascii(), "a"))!=NULL)
	{
		// results
		fprintf(fp, "Time:%f, Avg Error:%.3f, Max Error:%.3f, Std Error:%.3f, Sen:%.3f, Spe:%.3f, Di:%.3f, Ha:%.3f, Sm:%.3f\n", 
			lw_elapsed_time, lw_avgError, lw_maxError, lw_stdError, lw_sensitivity, lw_specificity, lw_dice, lw_hausdorff, lw_smoothness);

		fclose(fp);
	}
	return;
}

#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/SoQtObject.h>
#include <Inventor/SoDB.h>
#include <Inventor/SoInteraction.h>
#include <Inventor/nodekits/SoNodeKit.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoEventCallback.h>
#include <Inventor/nodes/SoText3.h>
#include <Inventor/nodes/SoNormal.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#include <Inventor/nodes/SoShapeHints.h>
#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/details/SoDetail.h>
#include <Inventor/details/SoFaceDetail.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/actions/SoSearchAction.h>
#include <Inventor/nodes/SoInfo.h>
#include <Inventor/nodes/SoTriangleStripSet.h>
#include "IvMesh.h"

// open and show an iv file in our window
void NIH_BoneSegmentation_Dlg::pushButton_surface_iv_clicked()
{
	UpdateData();

	QString strNewFileName = QFileDialog::getOpenFileName(
		QString::null,
        "IV File (*.iv)",
         this,
         "Load iv file",
         "Choose a file" );
 
	if(strNewFileName!=QString::null) 
	{
		strcpy(base_fn,strNewFileName.ascii());
	
	    SoInput input;
	    SoSeparator  * _sceneGraph;
		IvMesh * _mesh;
	    SoDB::init();

		if (input.openFile( base_fn ) )
		{
	        _sceneGraph = SoDB::readAll(&input);
		    _sceneGraph->ref();
	        _mesh = new IvMesh(*_sceneGraph);

			SoCoordinate3 &mesh_vertA = _mesh->getVertexArray();
			SoIndexedFaceSet &mesh_faceA = _mesh->getTriangleSet();
			int vertNum = mesh_vertA.point.getNum();
			int faceNum = mesh_faceA.coordIndex.getNum();
			faceNum = faceNum/4;

			vec3DynArray v_array;
			intVec3DynArray f_array;

			// convert the wing-edge structure and display the surface
/*			int k, a, b,c;

			v_array.SetSize(vertNum);
	        for (k=0;  k<vertNum;  ++k )
		    {
			    SbVec3f vec = mesh_vertA.point[k];
				v_array[k] = Vec3(vec[0], vec[1], vec[2]);
			}

			f_array.SetSize(faceNum);
	        for (k=0;  k<faceNum;  ++k )
		    {
		        a = mesh_faceA.coordIndex[k*4+0];
				b = mesh_faceA.coordIndex[k*4+1];
				c = mesh_faceA.coordIndex[k*4+2];
	           
				if ( a!=b && a!=c && b!=c )
				{
					f_array[k] = IntVec3(a,b,c);
				}
	       }

			spineSurf->BuildWingedEdgeStructure(v_array,f_array);
		
			NIH_OpenGL_CentralizeModel(display_id2, spineSurf);
			NIH_OpenGL_Refresh_Display(display_id2);
*/
			char fn_path1[200];
			strcpy(fn_path1,strNewFileName.ascii());
			// get the directory
			int i;
			for(i=strlen(fn_path1)-1; i>=0; i--)
			{
				if(fn_path1[i]=='\\' || fn_path1[i]=='/') break;
			}
			fn_path1[i+1]=0;

			QDir::setCurrent(fn_path1);

			// load sovira's fold file and convert to map file for Ananda
			int number_of_pts=0;
			float *x_array=NULL, *y_array=NULL, *z_array=NULL;

			FILE *fp;
			if((fp=fopen("foldNumber.txt","r"))!=NULL) 
			{
				char line[100];
				fgets(line, 100, fp);
				fclose(fp);
				QStringList li;
				li = QStringList::split('=', line);
				number_of_pts = li[1].toInt();
				printf("pts: %d", number_of_pts);
			} else printf("Can not open foldNumber.txt\n");

			if(number_of_pts>0) 
			{
				x_array = new float[number_of_pts];
				y_array = new float[number_of_pts];
				z_array = new float[number_of_pts];
			}
			else return;

			if((fp=fopen("FoldsX.out","rb"))!=NULL) 
			{
				fread(x_array, sizeof(float), number_of_pts, fp);
				fclose(fp);
			} else printf("Can not open FoldsX.out\n");

			if((fp=fopen("FoldsY.out","rb"))!=NULL) 
			{
				fread(y_array, sizeof(float), number_of_pts, fp);
				fclose(fp);
			} else printf("Can not open FoldsY.out\n");

			if((fp=fopen("FoldsZ.out","rb"))!=NULL) 
			{
				fread(z_array, sizeof(float), number_of_pts, fp);
				fclose(fp);
			} else printf("Can not open FoldsZ.out\n");

			// output the coordinates to a single txt file
			if((fp=fopen("Folds_coord.txt", "w"))!=NULL)
			{
				for(i=0; i<number_of_pts; i++)
				{
					fprintf(fp, "%f %f %f\n", x_array[i], y_array[i], z_array[i]);
				}
				fclose(fp);
			}
			return;

			// find the bounding box
			float bbx0, bbx1, bby0, bby1, bbz0, bbz1;
			bbx0 = 10000; bbx1=-1; bby0=10000; bby1=-1; bbz0=10000; bbz1=-1;
			for(i=0; i<vertNum; i++)
			{
				if(mesh_vertA.point[i][0]<bbx0) bbx0=mesh_vertA.point[i][0];
				if(mesh_vertA.point[i][0]>bbx1) bbx1=mesh_vertA.point[i][0];
				if(mesh_vertA.point[i][1]<bby0) bby0=mesh_vertA.point[i][1];
				if(mesh_vertA.point[i][1]>bby1) bby1=mesh_vertA.point[i][1];
				if(mesh_vertA.point[i][2]<bbz0) bbz0=mesh_vertA.point[i][2];
				if(mesh_vertA.point[i][2]>bbz1) bbz1=mesh_vertA.point[i][2];
			}
			printf("bb1: %.2f %.2f; %.2f %.2f; %.2f %.2f\n", bbx0, bbx1, bby0, bby1, bbz0, bbz1);

			bbx0 = 10000; bbx1=-1; bby0=10000; bby1=-1; bbz0=10000; bbz1=-1;
			for(i=0; i<number_of_pts; i++)
			{
				x_array[i] /=10;
				y_array[i] /=10;
				z_array[i] /=10;
				if(x_array[i]<bbx0) bbx0=x_array[i];
				if(x_array[i]>bbx1) bbx1=x_array[i];
				if(y_array[i]<bby0) bby0=y_array[i];
				if(y_array[i]>bby1) bby1=y_array[i];
				if(z_array[i]<bbz0) bbz0=z_array[i];
				if(z_array[i]>bbz1) bbz1=z_array[i];
			}
			printf("bb1: %.2f %.2f; %.2f %.2f; %.2f %.2f\n", bbx0, bbx1, bby0, bby1, bbz0, bbz1);

			// for each fold pt, search for the closest pt on the surface
			char *label;
			label = new char[vertNum];

			float tol=0.3;
			int j, minInd;
			float dist, minDist;

			for(i=0; i<vertNum; i++) label[i]=0;
			for(i=0; i<number_of_pts; i++)
			{
				minDist = 1000;
				minInd = -1;
//				x_array[i] /=10;
//				y_array[i] /=10;
//				z_array[i] /=10;
				for(j=0; j<vertNum; j++)
				{
					if(fabs(x_array[i]-mesh_vertA.point[j][0])<tol && fabs(y_array[i]-mesh_vertA.point[j][1])<tol &&
						fabs(z_array[i]-mesh_vertA.point[j][2])<tol)
					{
						dist = fabs(x_array[i]-mesh_vertA.point[j][0])+fabs(y_array[i]-mesh_vertA.point[j][1])+fabs(z_array[i]-mesh_vertA.point[j][2]);

						if(dist<minDist)
						{
							minDist =dist;
							minInd = j;
						}
					}
				}

				if(minInd>=0) 
				{
					label[minInd]=1;
///					printf("V %d: %.2f %.2f %.2f; D: %.2f; I: %d\n", i, x_array[i], y_array[i], z_array[i], minDist, minInd);
				}
				else printf("V %d: %.2f %.2f %.2f,  no match\n", i, x_array[i], y_array[i], z_array[i]);
			}

			// output the map
			if((fp=fopen("folds.map","w"))!=NULL)
			{
				for(i=0; i<vertNum; i++) 
				{
					if(label[i]>0) fprintf(fp, "%d %d\n", i, label[i]);
				}
				fclose(fp);
			}

			delete x_array;
			delete y_array;
			delete z_array;
			delete label;
		}
	
	}

}


void NIH_BoneSegmentation_Dlg::pushButton_Build_ALM_clicked()
{
	if(curStudy<0 || studyEntries[curStudy].organPath=="NULL" || img3D==NULL) return;
    QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );

	// spine not partitioned yet or t12 not found
	// segment and partition the spine
	if(pedicleValley.GetSize()<1 || segInfo.t12_vertebra<0) 
	{
		pushButton_pre_process_clicked();
		pushButton_spinal_cord_clicked();
		pushButton_spine_partition_clicked();
		pushButton_Rib_Detection_clicked();
		pushButton_surface_clicked();
		pushButton_Load_Organ_clicked();
	}

	int x, y, z;
	float px, py, pz;
	int i;

	px = img3D->Get_Pixel_SizeX();
	py = img3D->Get_Pixel_SizeY();
	pz = img3D->Get_Pixel_SizeZ();

	// first fit the center of disk into a smooth bspline curve
	vec3DynArray vlist, tangent_list, normal_list, binormal_list;
	doubleDynArray curvature_list, torsion_list;

	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(segInfo.diskCenter[z].x>0)
		{
///			vlist.Add(Vec3(segInfo.diskCenter[z].x*px, segInfo.diskCenter[z].y*py, z*pz)); 
			vlist.Add(Vec3(segInfo.diskCenter[z].x*px, segInfo.diskCenter[z].y*py, img3D->GetSlicePosition(z))); 
		}
	}

	vlist = spinalCord3D;

	NIH_Curve_Local_Frame(5, vlist, tangent_list, normal_list, binormal_list, curvature_list, torsion_list);

	singleALM[0].smoothedCord->UpdateControlList(vlist);
	singleALM[0].smoothedCord->SetCurrColor(Vec3(0,0,1));

	float closestDist, dist, centerZ;
	int closestZ, cv, centerZi;
	Vec3 vertSize;

	for(i=0; i<pedicleValley.GetSize()-1; i++)
	{
		cv = 11-segInfo.t12_vertebra+i;
		if(singleALM[0].lfw.tL[cv]==NULL) singleALM[0].lfw.tL[cv] = new VertebraLocationModel();

		centerZi = (pedicleValley[i]+pedicleValley[i+1])/2;
//		centerZ = centerZi*pz;
		centerZ = img3D->GetSlicePosition(centerZi);

		// find the closest point on the curve of spinal cord
		closestDist = 100000;
		closestZ = -1;
		for(z=0; z<vlist.GetSize(); z++)
		{
			dist = fabs((float)(vlist[z].z-centerZ));
			if(dist<closestDist)
			{
				closestDist = dist;
				closestZ = z;
			}
		}

		// set the center and orientation
		// may need to change the orientation to use the process information
		//
		singleALM[0].lfw.tL[cv]->center = vlist[closestZ];
		singleALM[0].lfw.tL[cv]->orientation.Rx = binormal_list[closestZ];
		singleALM[0].lfw.tL[cv]->orientation.Ry = normal_list[closestZ];
		singleALM[0].lfw.tL[cv]->orientation.Rz = tangent_list[closestZ];

		// use the mean of all disk size to set size
		vertSize = Vec3(0,0,0);
		for(z=pedicleValley[i]; z<=pedicleValley[i+1]; z++)
		{
	//		vertSize.x += segInfo.diskRadius[z].x*2;
			vertSize.y += segInfo.diskRadius[z].y*2;
		}
	//	vertSize.x /= (double)(pedicleValley[i+1]-pedicleValley[i]);
		vertSize.y /= (double)(pedicleValley[i+1]-pedicleValley[i]);
		vertSize.x = vertSize.y;	// force it to be same as y

		vertSize.z = (pedicleValley[i+1]-pedicleValley[i])*pz;

		singleALM[0].lfw.tL[cv]->size = vertSize;

		// set name
		if(cv<12) sprintf(singleALM[0].lfw.tL[cv]->name, "T%d",cv+1);
		else sprintf(singleALM[0].lfw.tL[cv]->name, "L%d",cv-11); 


		// set the visualization model
		if(cv>=9 && cv<=14)
		{
			singleALM[0].vertebra_3d_model[cv-9]->centroid = singleALM[0].lfw.tL[cv]->center;
			singleALM[0].vertebra_3d_model[cv-9]->loc = singleALM[0].lfw.tL[cv]->center;
			singleALM[0].vertebra_3d_model[cv-9]->dirX = singleALM[0].lfw.tL[cv]->orientation.Rx;
			singleALM[0].vertebra_3d_model[cv-9]->dirY = singleALM[0].lfw.tL[cv]->orientation.Ry;
			singleALM[0].vertebra_3d_model[cv-9]->dirZ = singleALM[0].lfw.tL[cv]->orientation.Rz;
			singleALM[0].vertebra_3d_model[cv-9]->axis = singleALM[0].lfw.tL[cv]->orientation.Rz;
			singleALM[0].vertebra_3d_model[cv-9]->sizeX = singleALM[0].lfw.tL[cv]->size.x;
			singleALM[0].vertebra_3d_model[cv-9]->sizeY = singleALM[0].lfw.tL[cv]->size.y;
			singleALM[0].vertebra_3d_model[cv-9]->sizeZ = singleALM[0].lfw.tL[cv]->size.z*0.9;
			singleALM[0].vertebra_3d_model[cv-9]->radius = (singleALM[0].lfw.tL[cv]->size.x+singleALM[0].lfw.tL[cv]->size.y)/8;
			singleALM[0].vertebra_3d_model[cv-9]->height = singleALM[0].lfw.tL[cv]->size.z*0.9;
			singleALM[0].vertebra_3d_model[cv-9]->ComputeFrame();
			singleALM[0].vertebra_3d_model[cv-9]->draw_mode = LINE_MODE;
		}
	}

	if(singleALM[0].smoothedCord->draw_mode==LINE_MODE)
	{
		singleALM[0].vertebra_3d_model[2]->SetCurrColor(Vec3(0, 1, 0));
		NIH_OpenGL_CentralizeModel(display_id2, singleALM[0].smoothedCord);
	}

	NIH_OpenGL_Refresh_Display(display_id2);

    QApplication::restoreOverrideCursor();
	return;
}

// load organ and generate surfaces
void NIH_BoneSegmentation_Dlg::pushButton_Load_Organ_clicked()
{
	if(curStudy<0 || studyEntries[curStudy].organPath=="NULL" || img3D==NULL) return;
    QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );

	QString organPath, organFn, organName, organSurfFn;
	CIS_Array_Image3D_short *organImg=NULL;
	short maxPixelValue;
	int z, i, m;
	QDir organDir;

	for(m=0; m<5; m++) 
	{
		singleALM[0].organSurf[m]->RemoveAllMesh();
		singleALM[0].organSurf[m]->SetDrawMode(HIDE_MODE);
		singleALM[0].organLocationModel[m]->SetDrawMode(HIDE_MODE);
	}

	for(m=0; m<5; m++)
	{
		if(m==0)
		{	// load left kidney
			organName = "Left Kidney";
			organPath = studyEntries[curStudy].organPath+"\\Left Kidney 100\\";
		}
		else if(m==1)
		{	// load right kidney
			organName = "Right Kidney";
			organPath = studyEntries[curStudy].organPath+"\\Right Kidney 100\\";
		}
		else if(m==2)
		{	// load Liver
			organName = "Liver";
			organPath = studyEntries[curStudy].organPath+"\\Liver 100\\";
		}
		else if(m==3)
		{	// load pancreas
			organName = "Pancreas";
			organPath = studyEntries[curStudy].organPath+"\\Pancreas 100\\";
		}
		else if(m==4)
		{	// load spleen
			organName = "Spleen";
			organPath = studyEntries[curStudy].organPath+"\\Spleen 100\\";
		}
	
		organSurfFn = organPath+organName+".off";

		organDir = organPath;
		organDir.setFilter(QDir::Files|QDir::Readable);
		organDir.setNameFilter("*.img");
		QStringList imgFiles = organDir.entryList();

		// binarized segmented data
		if(imgFiles.size()<1)
		{
			Error("organ file cannot found in ", organPath.ascii());
		}
		else 
		{
			organFn = organPath+imgFiles[0];
			if(organImg!=NULL) delete organImg;
			organImg = new CIS_Array_Image3D_short(organFn.ascii());

			maxPixelValue = organImg->Maximum_Pixel_Value();
			if(maxPixelValue==1) 
			{
				organImg->Subst_Pixel_Value(maxPixelValue, 100);
				maxPixelValue = 100;
			}

			for(z=0; z<organImg->Num_Levels() && z<img3D->Num_Levels(); z++)
				organImg->SetSlicePosition(z, img3D->GetSlicePosition(z));

			// generate surface
			CIS_Algo_MarchingCube(maxPixelValue/2, -1, -1, 4, 3, organImg, (char *)(organSurfFn.ascii()), singleALM[0].organSurf[m],-1, -1, -1, -1, -1, -1);
			TriangleMesh_LaplacianSmooth(singleALM[0].organSurf[m], 1, 10);

			singleALM[0].organSurf[m]->ComputeSize();
			singleALM[0].organSurf[m]->SetDrawMode(SURFACE_MODE);

			// set location model
			singleALM[0].organLocationModel[m]->pt_list1.SetSize(0);
			singleALM[0].organLocationModel[m]->pt_list2.SetSize(0);
			for(i=0; i<6; i++)
			{
				if(singleALM[0].vertebra_3d_model[i]!=NULL && singleALM[0].vertebra_3d_model[i]->sizeX>0) 
				{
					singleALM[0].organLocationModel[m]->pt_list1.Add(singleALM[0].vertebra_3d_model[i]->centroid);
					singleALM[0].organLocationModel[m]->pt_list2.Add(singleALM[0].organSurf[m]->centroid);

					// need the relative position and orientation
				}
			}
			singleALM[0].organLocationModel[m]->SetDrawMode(LINE_MODE);

			// set the organ location model
			if(singleALM[0].lfw.organ[m]==NULL) singleALM[0].lfw.organ[m] = new OrganLocationModel();
			strcpy(singleALM[0].lfw.organ[m]->name, organName.ascii());
			singleALM[0].lfw.organ[m]->center = singleALM[0].organSurf[m]->centroid;
			singleALM[0].lfw.organ[m]->size = Vec3(singleALM[0].organSurf[m]->sizeX, 
				singleALM[0].organSurf[m]->sizeY, singleALM[0].organSurf[m]->sizeZ);
			singleALM[0].lfw.organ[m]->orientation.Rx = singleALM[0].organSurf[m]->dirX;
			singleALM[0].lfw.organ[m]->orientation.Ry = singleALM[0].organSurf[m]->dirY;
			singleALM[0].lfw.organ[m]->orientation.Rz = singleALM[0].organSurf[m]->dirZ;
			strcpy(singleALM[0].lfw.organ[m]->surfaceFilename, organSurfFn.ascii());
			strcpy(singleALM[0].lfw.organ[m]->binaryFilename, organFn.ascii());

			// set relative location model
			for(i=0; i<6; i++)
			{
				int cv = i+9;
				if(singleALM[0].vertebra_3d_model[i]!=NULL && singleALM[0].vertebra_3d_model[i]->sizeX>0) 
				{
					singleALM[0].lfw.organ[m]->toTL[i].relativeLocation = singleALM[0].organSurf[m]->centroid-
						singleALM[0].vertebra_3d_model[i]->centroid;
					singleALM[0].lfw.organ[m]->toTL[i].relativeOrientation = singleALM[0].lfw.organ[m]->orientation*
						~singleALM[0].lfw.tL[cv]->orientation;
				}
			}

			// compute the Gaussian model, mean and variance
			//
			double meanI, varI;
			int count;
			int sizex, sizey, sizez;
			sizex = img3D->Num_Cols();
			sizey = img3D->Num_Rows();
			sizez = img3D->Num_Levels();
			int sizexo, sizeyo, sizezo;
			sizexo = organImg->Num_Cols();
			sizeyo = organImg->Num_Rows();
			sizezo = organImg->Num_Levels();

			int sizexyz=sizex*sizey*sizez;
			int sizexyzo=sizexo*sizeyo*sizezo;
			short *imgA=img3D->GetArray();
			short *organA = organImg->GetArray();

			meanI = varI = 0;
			count =0;
			for(int k=0; k<sizexyzo; k++)
			{
				if(organA[k]>0) 
				{
					if(imgA[k]>900 && imgA[k]<1200)	// don't include outliers
					{
						meanI += imgA[k];
						varI += (double)imgA[k]*(double)imgA[k];
						count++;
					}
				}
			}
			if(count>0)
			{
				meanI /= (double)count;
				varI = varI/(double)count-meanI*meanI;
				varI = sqrt(varI);
				singleALM[0].lfw.organ[m]->meanIntensity = (float)meanI;
				singleALM[0].lfw.organ[m]->varianceIntensity = (float)varI;
			}
		}	// else imgFiles
	}	// for m;

	NIH_OpenGL_Refresh_Display(display_id2);

	if(organImg) delete organImg;

    QApplication::restoreOverrideCursor();
	return;
}

// Save a single location model from one dataset using xml format
void NIH_BoneSegmentation_Dlg::pushButton_Save_ALM_Model_clicked()
{
	// get the file name
	QString strNewFileName = QFileDialog::getSaveFileName(
		QString::null,
        "Model File (*.xml)",
         this,
         "Save Model",
         "Choose a file" );
 
	if(strNewFileName==QString::null) return;

	SaveSingleLocationModel(strNewFileName.ascii());
}

void NIH_BoneSegmentation_Dlg::SaveSingleLocationModel(const char *fn)
{
	if(img3D==NULL || curStudy<0) return;

	int i;
	FILE *fp = fopen(fn,"w");
	
	if(fp == NULL)
	{
		fprintf(stderr,"Can't open xml model file %s\n",fn);
		return;
	}

	printf("here, %s\n", fn);
	fprintf(fp, "<ALMSingleModel>\n");

	// Save information about the study if available
	if(curStudy>=0 && curStudy<numStudyEntries)
	{
		fprintf(fp,"\t<studyInfo>\n");
		fprintf(fp,"\t\t<patientRoot> %s </patientRoot>\n", patientRoot);
		fprintf(fp,"\t\t<patientName> %s </patientName>\n", studyEntries[curStudy].patientName.ascii());
		fprintf(fp,"\t\t<studyDate> %s </studyDate>\n", studyEntries[curStudy].studyDate.ascii());
		fprintf(fp,"\t\t<studyId> %s </studyId>\n", studyEntries[curStudy].studyId.ascii());
		fprintf(fp,"\t\t<localImagePath> %s </localImagePath>\n", studyEntries[curStudy].localImagePath.ascii());
		fprintf(fp,"\t\t<organPath> %s </organPath>\n", studyEntries[curStudy].organPath.ascii());
		fprintf(fp, "\t</studyInfo>\n");
	}

	// Save the image information if available
	if(img3D!=NULL && img3D->Num_Levels()>0)
	{
		fprintf(fp,"\t<imageInfo>\n");
		fprintf(fp,"\t\t<imageSize> %d %d %d </imageSize>\n", img3D->Num_Cols(), img3D->Num_Rows(), img3D->Num_Levels());
		fprintf(fp,"\t\t<pixelSize> %f %f %f </pixelSize>\n", img3D->Get_Pixel_SizeX(), img3D->Get_Pixel_SizeY(), img3D->Get_Pixel_SizeZ());
		fprintf(fp,"\t\t<zOrigin> %f </zOrigin>\n", img3D->GetSlicePosition(0));
		if(img3D->GetSlicePosition(img3D->Num_Levels()-1)-img3D->GetSlicePosition(0)>0) fprintf(fp,"\t\t<zOrient> 1 </zOrient>\n");
		else fprintf(fp,"\t\t<zOrient> -1 </zOrient>\n");
		fprintf(fp,"\t</imageInfo>\n");
	}


	// Save some global spine segmentation information
	fprintf(fp,"\t<spineSegInfo>\n");
	fprintf(fp,"\t\t<bound1> %d %d %d </bound1>\n", segInfo.bound1.x, segInfo.bound1.y, segInfo.bound1.z);
	fprintf(fp,"\t\t<bound2> %d %d %d </bound2>\n", segInfo.bound2.x, segInfo.bound2.y, segInfo.bound2.z);
	fprintf(fp,"\t\t<sacrum_start> %d </sacrum_start>\n", segInfo.sacrum_start);
	fprintf(fp,"\t\t<t12_vertebra> %d </t12_vertebra>\n", segInfo.t12_vertebra);
	fprintf(fp,"\t\t<ribCageSize> %f %f %f </ribCageSize>\n", segInfo.ribCageSize.x, segInfo.ribCageSize.y, segInfo.ribCageSize.z);
	fprintf(fp,"\t</spineSegInfo>\n");

	// Save the Spine location model
	if(spinalCord3D.GetSize()>0)
	{
		fprintf(fp,"\t<spineLocation>\n");

		// save spine partition
		fprintf(fp,"\t\t<spinePartition size=\"%d\">\n", pedicleValley.GetSize());
		for(i=0; i<pedicleValley.GetSize(); i++)
		{
			fprintf(fp,"\t\t\t<v> %d </v>\n", pedicleValley[i]);
		}
		fprintf(fp,"\t\t</spinePartition>\n");

		// save spinal cord
		fprintf(fp,"\t\t<spinalCord size=\"%d\">\n", spinalCord3D.GetSize());
		for(i=0; i<spinalCord3D.GetSize(); i++)
		{
			fprintf(fp,"\t\t\t<v> %f %f %f </v>\n", spinalCord3D[i].x, spinalCord3D[i].y, spinalCord3D[i].z);
		}
		fprintf(fp,"\t\t</spinalCord>\n");

		// save the vertibra frame models
		for(i=0; i<17; i++)
		{
			if(singleALM[0].lfw.tL[i]!=NULL)
			{
				fprintf(fp,"\t\t<vertebra index=\"%d\">\n", i);
				fprintf(fp,"\t\t\t<center> %f %f %f </center>\n", singleALM[0].lfw.tL[i]->center.x, singleALM[0].lfw.tL[i]->center.y, singleALM[0].lfw.tL[i]->center.z);
				fprintf(fp,"\t\t\t<orientation> %f %f %f %f %f %f %f %f %f</orientation>\n", 
					singleALM[0].lfw.tL[i]->orientation.Rx.x, singleALM[0].lfw.tL[i]->orientation.Rx.y, singleALM[0].lfw.tL[i]->orientation.Rx.z,
					singleALM[0].lfw.tL[i]->orientation.Ry.x, singleALM[0].lfw.tL[i]->orientation.Ry.y, singleALM[0].lfw.tL[i]->orientation.Ry.z,
					singleALM[0].lfw.tL[i]->orientation.Rz.x, singleALM[0].lfw.tL[i]->orientation.Rz.y, singleALM[0].lfw.tL[i]->orientation.Rz.z
					);
				fprintf(fp,"\t\t\t<size> %f %f %f </size>\n", singleALM[0].lfw.tL[i]->size.x, singleALM[0].lfw.tL[i]->size.y, singleALM[0].lfw.tL[i]->size.z);
				fprintf(fp,"\t\t</vertebra>\n");
			}
		}

		fprintf(fp,"\t</spineLocation>\n");
	}	// if spinalCord3D

	// Save the Organ location model
	fprintf(fp,"\t<organLocation>\n");
	for(i=0; i<5; i++)
	{
		if(singleALM[0].lfw.organ[i]!=NULL)
		{
			fprintf(fp,"\t\t<organ name=\"%s\">\n", singleALM[0].lfw.organ[i]->name);
			fprintf(fp,"\t\t\t<center> %f %f %f </center>\n", singleALM[0].lfw.organ[i]->center.x, singleALM[0].lfw.organ[i]->center.y, singleALM[0].lfw.organ[i]->center.z);
			fprintf(fp,"\t\t\t<orientation> %f %f %f %f %f %f %f %f %f</orientation>\n", 
				singleALM[0].lfw.organ[i]->orientation.Rx.x, singleALM[0].lfw.organ[i]->orientation.Rx.y, singleALM[0].lfw.organ[i]->orientation.Rx.z,
				singleALM[0].lfw.organ[i]->orientation.Ry.x, singleALM[0].lfw.organ[i]->orientation.Ry.y, singleALM[0].lfw.organ[i]->orientation.Ry.z,
				singleALM[0].lfw.organ[i]->orientation.Rz.x, singleALM[0].lfw.organ[i]->orientation.Rz.y, singleALM[0].lfw.organ[i]->orientation.Rz.z
				);
			fprintf(fp,"\t\t\t<size> %f %f %f </size>\n", singleALM[0].lfw.organ[i]->size.x, singleALM[0].lfw.organ[i]->size.y, singleALM[0].lfw.organ[i]->size.z);
			fprintf(fp,"\t\t\t<Gaussian> %f %f </Gaussian>\n", singleALM[0].lfw.organ[i]->meanIntensity, singleALM[0].lfw.organ[i]->varianceIntensity);
			fprintf(fp,"\t\t\t<surfaceFilename> %s </surfaceFilename>\n", singleALM[0].lfw.organ[i]->surfaceFilename);
			fprintf(fp,"\t\t\t<binaryFilename> %s </binaryFilename>\n", singleALM[0].lfw.organ[i]->binaryFilename);

			fprintf(fp,"\t\t</organ>\n");
		}
	}
	fprintf(fp,"\t</organLocation>\n");

	fprintf(fp, "</ALMSingleModel>\n");
	fclose(fp);

	return;
}

// load a single location model given the xml file
//
int NIH_BoneSegmentation_Dlg::LoadSingleALMModel(const char *fn, int index)
{
	FILE *fp = fopen(fn,"r");
	
	if(fp == NULL)
	{
		fprintf(stderr,"Can't open xml model file %s\n",fn);
		return -1;
	}
	fclose(fp);

	strcpy(singleALM[index].modelFileName, fn);

	int i;
	QStringList valueList;
    // read the xml file
    XmlReader xmlReader(fn);
    QDomElement  ALMSingleModelElement = xmlReader.getDomElementUnique("ALMSingleModel", xmlReader.getTopNode());

	QDomElement studyInfoElement = xmlReader.getDomElementUnique("studyInfo",ALMSingleModelElement);
/*	printf("%s\n", studyInfoElement.elementsByTagName("patientRoot").item(0).toElement().text().ascii());
	printf("%s\n", studyInfoElement.elementsByTagName("patientName").item(0).toElement().text().ascii());
	printf("%s\n", studyInfoElement.elementsByTagName("studyDate").item(0).toElement().text().ascii());
	printf("%s\n", studyInfoElement.elementsByTagName("studyId").item(0).toElement().text().ascii());
	printf("%s\n", studyInfoElement.elementsByTagName("localImagePath").item(0).toElement().text().ascii());
	printf("%s\n", studyInfoElement.elementsByTagName("organPath").item(0).toElement().text().ascii());
*/
	QDomElement imageInfoElement = xmlReader.getDomElementUnique("imageInfo",ALMSingleModelElement);
	float zOrigin, zOrient;
	float pixelSizex, pixelSizey, pixelSizez;
	int isizex, isizey, isizez;

/*	printf("%s\n", imageInfoElement.elementsByTagName("imageSize").item(0).toElement().text().ascii());
	printf("%s\n", imageInfoElement.elementsByTagName("pixelSize").item(0).toElement().text().ascii());
	printf("%s\n", imageInfoElement.elementsByTagName("zOrigin").item(0).toElement().text().ascii());
	printf("%s\n", imageInfoElement.elementsByTagName("zOrient").item(0).toElement().text().ascii());
*/
	valueList = QStringList::split(' ', imageInfoElement.elementsByTagName("imageSize").item(0).toElement().text().simplifyWhiteSpace());
	isizex = valueList[0].toInt();
	isizey = valueList[1].toInt();
	isizez = valueList[2].toInt();
	valueList = QStringList::split(' ', imageInfoElement.elementsByTagName("pixelSize").item(0).toElement().text().simplifyWhiteSpace());
	pixelSizex = valueList[0].toFloat();
	pixelSizey = valueList[1].toFloat();
	pixelSizez = valueList[2].toFloat();
	valueList = QStringList::split(' ', imageInfoElement.elementsByTagName("zOrigin").item(0).toElement().text().simplifyWhiteSpace());
	zOrigin = valueList[0].toFloat();
	valueList = QStringList::split(' ', imageInfoElement.elementsByTagName("zOrient").item(0).toElement().text().simplifyWhiteSpace());
	zOrient = valueList[0].toFloat();

	QDomElement spineSegInfoElement = xmlReader.getDomElementUnique("spineSegInfo",ALMSingleModelElement);
/*	printf("%s\n", spineSegInfoElement.elementsByTagName("bound1").item(0).toElement().text().ascii());
	printf("%s\n", spineSegInfoElement.elementsByTagName("bound2").item(0).toElement().text().ascii());
	printf("%s\n", spineSegInfoElement.elementsByTagName("sacrum_start").item(0).toElement().text().ascii());
	printf("%s\n", spineSegInfoElement.elementsByTagName("t12_vertebra").item(0).toElement().text().ascii());
*/
	valueList = QStringList::split(' ', spineSegInfoElement.elementsByTagName("ribCageSize").item(0).toElement().text().simplifyWhiteSpace());
	singleALM[index].ribCageSize.x = valueList[0].toDouble();
	singleALM[index].ribCageSize.y = valueList[1].toDouble();
	singleALM[index].ribCageSize.z = valueList[2].toDouble();
	
	QDomElement spineLocationElement = xmlReader.getDomElementUnique("spineLocation",ALMSingleModelElement);

	// set the spinal cord
	doubleDynArray curvature_list, torsion_list;
	vec3DynArray vlist, tangent_list, normal_list, binormal_list;
	QDomElement spinalCordElement = xmlReader.getDomElementUnique("spinalCord",spineLocationElement);
	vlist.SetSize(spinalCordElement.attribute("size").toInt());
    QDomNodeList vElementList    = xmlReader.getDomNodeList("v", spinalCordElement);
	for(i=0; i<vElementList.count(); i++)
	{
///		printf("%s\n", vElementList.item(i).toElement().text().ascii());
		valueList = QStringList::split(' ', vElementList.item(i).toElement().text().simplifyWhiteSpace());
		vlist[i].x = valueList[0].toDouble();
		vlist[i].y = valueList[1].toDouble();
		vlist[i].z = valueList[2].toDouble();
	}
	NIH_Curve_Local_Frame(5, vlist, tangent_list, normal_list, binormal_list, curvature_list, torsion_list);
	singleALM[index].smoothedCord->UpdateControlList(vlist);
	singleALM[index].smoothedCord->loc = Vec3(0,0,0);
	singleALM[index].smoothedCord->dirX = Vec3(1,0,0);
	singleALM[index].smoothedCord->dirY = Vec3(0,1,0);
	singleALM[index].smoothedCord->dirZ = Vec3(0,0,1);
	singleALM[index].smoothedCord->SetCurrColor(Vec3(0,0,1));
	
	for(i=0; i<6; i++) singleALM[index].vertebra_3d_model[i]->draw_mode = HIDE_MODE;

	int cv;
    QDomNodeList vertElementList  = xmlReader.getDomNodeList("vertebra", spineLocationElement);
	for(i=0; i<vertElementList.count(); i++)
	{
		QDomElement vertElement = vertElementList.item(i).toElement();
///		printf("%s\n", vertElement.attribute("index").ascii());
		cv = vertElement.attribute("index").toInt();
		if(singleALM[index].lfw.tL[cv]==NULL) singleALM[index].lfw.tL[cv] = new VertebraLocationModel();

///		printf("%s\n", vertElement.elementsByTagName("center").item(0).toElement().text().ascii());
		valueList = QStringList::split(' ', vertElement.elementsByTagName("center").item(0).toElement().text().simplifyWhiteSpace());
		singleALM[index].lfw.tL[cv]->center.x = valueList[0].toDouble();
		singleALM[index].lfw.tL[cv]->center.y = valueList[1].toDouble();
		singleALM[index].lfw.tL[cv]->center.z = valueList[2].toDouble();

///		printf("%s\n", vertElement.elementsByTagName("orientation").item(0).toElement().text().ascii());
		valueList = QStringList::split(' ', vertElement.elementsByTagName("orientation").item(0).toElement().text().simplifyWhiteSpace());
		singleALM[index].lfw.tL[cv]->orientation.Rx.x = valueList[0].toDouble();
		singleALM[index].lfw.tL[cv]->orientation.Rx.y = valueList[1].toDouble();
		singleALM[index].lfw.tL[cv]->orientation.Rx.z = valueList[2].toDouble();
		singleALM[index].lfw.tL[cv]->orientation.Ry.x = valueList[3].toDouble();
		singleALM[index].lfw.tL[cv]->orientation.Ry.y = valueList[4].toDouble();
		singleALM[index].lfw.tL[cv]->orientation.Ry.z = valueList[5].toDouble();
		singleALM[index].lfw.tL[cv]->orientation.Rz.x = valueList[6].toDouble();
		singleALM[index].lfw.tL[cv]->orientation.Rz.y = valueList[7].toDouble();
		singleALM[index].lfw.tL[cv]->orientation.Rz.z = valueList[8].toDouble();

///		printf("%s\n", vertElement.elementsByTagName("size").item(0).toElement().text().ascii());
		valueList = QStringList::split(' ', vertElement.elementsByTagName("size").item(0).toElement().text().simplifyWhiteSpace());
		singleALM[index].lfw.tL[cv]->size.x = valueList[0].toDouble();
		singleALM[index].lfw.tL[cv]->size.y = valueList[1].toDouble();
		singleALM[index].lfw.tL[cv]->size.z = valueList[2].toDouble();
		if(cv<12) sprintf(singleALM[index].lfw.tL[cv]->name, "T%d",cv+1);
		else sprintf(singleALM[index].lfw.tL[cv]->name, "L%d",cv-11); 
		// set the visualization model
		if(cv>=9 && cv<=14)
		{
			singleALM[index].vertebra_3d_model[cv-9]->centroid = singleALM[index].lfw.tL[cv]->center;
			singleALM[index].vertebra_3d_model[cv-9]->loc = singleALM[index].lfw.tL[cv]->center;
			singleALM[index].vertebra_3d_model[cv-9]->dirX = singleALM[index].lfw.tL[cv]->orientation.Rx;
			singleALM[index].vertebra_3d_model[cv-9]->dirY = singleALM[index].lfw.tL[cv]->orientation.Ry;
			singleALM[index].vertebra_3d_model[cv-9]->dirZ = singleALM[index].lfw.tL[cv]->orientation.Rz;
			singleALM[index].vertebra_3d_model[cv-9]->axis = singleALM[index].lfw.tL[cv]->orientation.Rz;
			singleALM[index].vertebra_3d_model[cv-9]->sizeX = singleALM[index].lfw.tL[cv]->size.x;
			singleALM[index].vertebra_3d_model[cv-9]->sizeY = singleALM[index].lfw.tL[cv]->size.y;
			singleALM[index].vertebra_3d_model[cv-9]->sizeZ = singleALM[index].lfw.tL[cv]->size.z*0.9;
			singleALM[index].vertebra_3d_model[cv-9]->radius = (singleALM[index].lfw.tL[cv]->size.x+singleALM[index].lfw.tL[cv]->size.y)/8;
			singleALM[index].vertebra_3d_model[cv-9]->height = singleALM[index].lfw.tL[cv]->size.z*0.9;
			singleALM[index].vertebra_3d_model[cv-9]->ComputeFrame();
			singleALM[index].vertebra_3d_model[cv-9]->draw_mode = LINE_MODE;
		}
	}

	QString organPath, organFn, organName, organSurfFn;
	int m,v;
	for(m=0; m<5; m++) 
	{
		singleALM[index].organSurf[m]->RemoveAllMesh();
		singleALM[index].organSurf[m]->SetDrawMode(HIDE_MODE);
		singleALM[index].organLocationModel[m]->SetDrawMode(HIDE_MODE);
	}

	QDomElement organLocationElement = xmlReader.getDomElementUnique("organLocation",ALMSingleModelElement);
    QDomNodeList organElementList  = xmlReader.getDomNodeList("organ", organLocationElement);
	for(i=0; i<organElementList.count(); i++)
	{
		QDomElement organElement = organElementList.item(i).toElement();
///		printf("%s\n", organElement.attribute("name").ascii());
		organName = organElement.attribute("name").simplifyWhiteSpace();
		if(organName=="Left Kidney") m=0;
		else if(organName=="Right Kidney") m=1;
		else if(organName=="Liver") m=2;
		else if(organName=="Pancreas") m=3;
		else if(organName=="Spleen") m=4;

///		printf("%s\n", organElement.elementsByTagName("center").item(0).toElement().text().ascii());
///		printf("%s\n", organElement.elementsByTagName("orientation").item(0).toElement().text().ascii());
///		printf("%s\n", organElement.elementsByTagName("size").item(0).toElement().text().ascii());
///		printf("%s\n", organElement.elementsByTagName("surfaceFilename").item(0).toElement().text().ascii());

		organSurfFn = organElement.elementsByTagName("surfaceFilename").item(0).toElement().text().simplifyWhiteSpace();

		if(m>4) continue;	// only support 5 organs for now

		singleALM[index].organSurf[m]->ReadMeshFile(organSurfFn.ascii());
		TriangleMesh_LaplacianSmooth(singleALM[index].organSurf[m], 1, 10);
		singleALM[index].organSurf[m]->ComputeSize();
		singleALM[index].organSurf[m]->SetDrawMode(SURFACE_MODE);

		singleALM[index].organSurf[m]->loc = Vec3(0,0,0);
		singleALM[index].organSurf[m]->dirX = Vec3(1,0,0);
		singleALM[index].organSurf[m]->dirY = Vec3(0,1,0);
		singleALM[index].organSurf[m]->dirZ = Vec3(0,0,1);
		
		// set location model
		singleALM[index].organLocationModel[m]->pt_list1.SetSize(0);
		singleALM[index].organLocationModel[m]->pt_list2.SetSize(0);
		for(v=0; v<6; v++)
		{
			if(singleALM[index].vertebra_3d_model[v]!=NULL && singleALM[index].vertebra_3d_model[v]->sizeX>0) 
			{
				singleALM[index].organLocationModel[m]->pt_list1.Add(singleALM[index].vertebra_3d_model[v]->centroid);
				singleALM[index].organLocationModel[m]->pt_list2.Add(singleALM[index].organSurf[m]->centroid);
			}
		}
		singleALM[index].organLocationModel[m]->SetDrawMode(LINE_MODE);
		singleALM[index].organLocationModel[m]->loc = Vec3(0,0,0);
		singleALM[index].organLocationModel[m]->dirX = Vec3(1,0,0);
		singleALM[index].organLocationModel[m]->dirY = Vec3(0,1,0);
		singleALM[index].organLocationModel[m]->dirZ = Vec3(0,0,1);

		// set the organ location model
		if(singleALM[index].lfw.organ[m]==NULL) singleALM[index].lfw.organ[m] = new OrganLocationModel();
		strcpy(singleALM[index].lfw.organ[m]->name, organName.ascii());
		singleALM[index].lfw.organ[m]->center = singleALM[index].organSurf[m]->centroid;
		singleALM[index].lfw.organ[m]->size = Vec3(singleALM[index].organSurf[m]->sizeX, 
			singleALM[index].organSurf[m]->sizeY, singleALM[index].organSurf[m]->sizeZ);
		singleALM[index].lfw.organ[m]->orientation.Rx = singleALM[index].organSurf[m]->dirX;
		singleALM[index].lfw.organ[m]->orientation.Ry = singleALM[index].organSurf[m]->dirY;
		singleALM[index].lfw.organ[m]->orientation.Rz = singleALM[index].organSurf[m]->dirZ;
		strcpy(singleALM[index].lfw.organ[m]->surfaceFilename, organSurfFn.ascii());

		valueList = QStringList::split(' ', organElement.elementsByTagName("Gaussian").item(0).toElement().text().simplifyWhiteSpace());
		singleALM[index].lfw.organ[m]->meanIntensity = valueList[0].toFloat();
		singleALM[index].lfw.organ[m]->varianceIntensity = valueList[1].toFloat();

		// also load the binary file for location model 0
		// Model 0 is the reference model and it is used for transformation
		//
		if(index==0)
		{
			QString organFn = organElement.elementsByTagName("binaryFilename").item(0).toElement().text().simplifyWhiteSpace();
			CIS_Array_Image3D_short *organImg=NULL;
			organImg = new CIS_Array_Image3D_short(organFn.ascii());
			if(organMask==NULL) 
			{
				organMask = new CIS_Array_Image3D_uchar(isizex, isizey, isizez);
				organMask->Set_Pixel_SizeX(pixelSizex);
				organMask->Set_Pixel_SizeY(pixelSizey);
				organMask->Set_Pixel_SizeZ(pixelSizez);
				for(int z=0; z<isizez; z++)
				{
					organMask->SetSlicePosition(z, (float)z*zOrient+zOrigin);
				}
			}

			int sizex, sizey, sizez;
			sizex = organImg->Num_Cols();
			sizey = organImg->Num_Rows();
			sizez = organImg->Num_Levels();
			int k;
			int sizexyz=sizex*sizey*sizez;
			short *organA;
			unsigned char *maskA;
			organA = organImg->GetArray();
			maskA = organMask->GetArray();
			for(k=0; k<sizexyz; k++)
			{
				if(organA[k]>0) maskA[k] = m+1;
			}
			showOrganMask=true;
			if(organImg) delete organImg;
		}	// if index =0

		// load the ground Truth mask if the flag is set
		if(flagGT)
		{
			QString organFn = organElement.elementsByTagName("binaryFilename").item(0).toElement().text().simplifyWhiteSpace();
			CIS_Array_Image3D_short *organImg=NULL;
			organImg = new CIS_Array_Image3D_short(organFn.ascii());
			if(organGTMask==NULL) 
			{
				organGTMask = new CIS_Array_Image3D_uchar(isizex, isizey, isizez);
				organGTMask->Set_Pixel_SizeX(pixelSizex);
				organGTMask->Set_Pixel_SizeY(pixelSizey);
				organGTMask->Set_Pixel_SizeZ(pixelSizez);
				for(int z=0; z<isizez; z++)
				{
					organGTMask->SetSlicePosition(z, (float)z*zOrient+zOrigin);
				}
			}

			int sizex, sizey, sizez;
			sizex = organImg->Num_Cols();
			sizey = organImg->Num_Rows();
			sizez = organImg->Num_Levels();
			int k;
			int sizexyz=sizex*sizey*sizez;
			short *organA;
			unsigned char *maskA;
			organA = organImg->GetArray();
			maskA = organGTMask->GetArray();
			for(k=0; k<sizexyz; k++)
			{
				if(organA[k]>0) maskA[k] = m+1;
			}
			if(organImg) delete organImg;
		}	// if flagGT
	}	// if i

	// remap organMask to MaskImage
	if(organMask && maskImg3D && showOrganMask)
	{
		Frame f;
		total_probSum = total_probAvg = 0;
		for(int m=0; m<5; m++)
		{
			f=singleALM[index].organSurf[m]->GetFrame();
			if(!flagGT)
				MapOrganMask(organMask, f, m, 1, singleALM[index].lfw.organ[m]->meanIntensity, singleALM[index].lfw.organ[m]->varianceIntensity, 
					singleALM[index].lfw.organ[m]->probSum, singleALM[index].lfw.organ[m]->probAvg);
			else
			{
				MapOrganMask(organGTMask, f, m, 1, singleALM[index].lfw.organ[m]->meanIntensity, singleALM[index].lfw.organ[m]->varianceIntensity, 
					singleALM[index].lfw.organ[m]->probSum, singleALM[index].lfw.organ[m]->probAvg);
			}
			// print out for record
			printf("Organ %d: %f %f\n", m, singleALM[index].lfw.organ[m]->probSum, singleALM[index].lfw.organ[m]->probAvg);
			total_probSum += singleALM[index].lfw.organ[m]->probSum;
			total_probAvg += singleALM[index].lfw.organ[m]->probAvg;
		}
		printf("Total: %f %f\n", total_probSum, total_probAvg);
	}

	if(singleALM[index].vertebra_3d_model[2]->draw_mode==LINE_MODE)
	{
		singleALM[index].vertebra_3d_model[2]->SetCurrColor(Vec3(0, 1, 0));
	}

	// set location model for PDM model
	//
	if(ALM_PDM.pdmDegree>0)
	{
		int vecLength = 0;
		int commonOrgans=ALM_PDM.organIndexList.GetSize();
		int commonVerts=ALM_PDM.vertebraIndexList.GetSize();

		vecLength += commonOrgans*commonVerts*3;
		singleALM[index].locationVec.SetSize(vecLength);
		int i,j,m,v;
		int k=0;
		for(i=0; i<ALM_PDM.organIndexList.GetSize(); i++)
		{
			m=ALM_PDM.organIndexList[i];
			for(j=0; j<ALM_PDM.vertebraIndexList.GetSize(); j++)
			{
				v=ALM_PDM.vertebraIndexList[j];
				singleALM[index].locationVec[k] = singleALM[index].organSurf[m]->centroid.x+singleALM[index].organSurf[m]->loc.x-singleALM[index].vertebra_3d_model[v]->loc.x; k++;
				singleALM[index].locationVec[k] = singleALM[index].organSurf[m]->centroid.y+singleALM[index].organSurf[m]->loc.y-singleALM[index].vertebra_3d_model[v]->loc.y; k++;
				singleALM[index].locationVec[k] = singleALM[index].organSurf[m]->centroid.z+singleALM[index].organSurf[m]->loc.z-singleALM[index].vertebra_3d_model[v]->loc.z; k++;
			}
		}
		ALM_PDM.InversePDM(singleALM[index].locationVec, singleALM[index].pdmPara);
	}
	return 0;
}

// algin two location models based on their spine
// usually align the reference model with the ground truth model
void NIH_BoneSegmentation_Dlg::AlignTwoLocationModelsUsingSpine(int movingModel, int fixedModel)
{
	if(movingModel<0 || movingModel>=totalALMModels || fixedModel<0 || fixedModel>=totalALMModels) return;
	if(movingModel == fixedModel) return;

	// align the vertebra
	Vec3 delta_t[17], delta;
	delta = Vec3(0,0,0);
	int count=0;
	for(int v=0; v<17; v++)
	{
		delta_t[v]=Vec3(0,0,0);
		if(singleALM[movingModel].lfw.tL[v]!=NULL && singleALM[fixedModel].lfw.tL[v]!=NULL)
		{
			delta_t[v] = singleALM[fixedModel].lfw.tL[v]->center-singleALM[movingModel].lfw.tL[v]->center;

			delta += delta_t[v];
			count++;
			singleALM[movingModel].lfw.tL[v]->center=singleALM[fixedModel].lfw.tL[v]->center;
			if(v>=9 && v<=14)	// update the visual model
			{
				singleALM[movingModel].vertebra_3d_model[v-9]->centroid = singleALM[movingModel].lfw.tL[v]->center;
				singleALM[movingModel].vertebra_3d_model[v-9]->loc = singleALM[movingModel].lfw.tL[v]->center;
				singleALM[movingModel].vertebra_3d_model[v-9]->ComputeFrame();
			}
		}
	}	// for v

	if(count>0) delta/=(double)count;

	// align the spinal cord according to T12
	if(delta_t[11].len()>0)
	{
		singleALM[movingModel].smoothedCord->loc = delta_t[11];
		singleALM[movingModel].smoothedCord->ComputeFrame();
	}

	// update the organ location accordingly
	total_probSum = total_probAvg = 0;
	for(int m=0; m<5; m++)
	{
		singleALM[movingModel].organSurf[m]->loc += delta;
		singleALM[movingModel].organSurf[m]->ComputeFrame();
		for(int v=0; v<6; v++)
		{
			if(singleALM[movingModel].vertebra_3d_model[v]!=NULL)
			{
				singleALM[movingModel].organLocationModel[m]->pt_list1[v] = singleALM[movingModel].vertebra_3d_model[v]->loc;
				singleALM[movingModel].organLocationModel[m]->pt_list2[v] = singleALM[movingModel].organSurf[m]->centroid+
					singleALM[movingModel].organSurf[m]->loc;
			}
		}

		// remap to mask
		if(organMask && maskImg3D)
		{
			Frame f;
			f = singleALM[movingModel].organSurf[m]->GetFrame();
			MapOrganMask(organMask, f, m, 1, singleALM[movingModel].lfw.organ[m]->meanIntensity, 
				singleALM[movingModel].lfw.organ[m]->varianceIntensity, 
				singleALM[movingModel].lfw.organ[m]->probSum, 
				singleALM[movingModel].lfw.organ[m]->probAvg);
			// print out for record
			printf("Organ %d: %f %f\n", m, singleALM[movingModel].lfw.organ[m]->probSum, 
				singleALM[movingModel].lfw.organ[m]->probAvg);
			total_probSum += singleALM[movingModel].lfw.organ[m]->probSum;
			total_probAvg += singleALM[movingModel].lfw.organ[m]->probAvg;
		}
	}
	printf("Total: %f %f\n", total_probSum, total_probAvg);

}


void NIH_BoneSegmentation_Dlg::MapOrganMask(CIS_Array_Image3D_uchar *orgMask, Frame &f, int organIndex, int sampleRate, float probMean, float probStd, float &probSum, float &probAvg)
{
	if(orgMask==NULL || maskImg3D==NULL || img3D==NULL) return;

	int osx, osy, osz;
	float opx, opy, opz;
	unsigned char *organA;
	osx = orgMask->Num_Cols();
	osy = orgMask->Num_Rows();
	osz = orgMask->Num_Levels();
	opx = orgMask->Get_Pixel_SizeX();
	opy = orgMask->Get_Pixel_SizeY();
	opz = orgMask->Get_Pixel_SizeZ();
	organA = orgMask->GetArray();
	int osxy = osx*osy;

	int msx, msy, msz;
	float mpx, mpy, mpz;
	short *maskA;
	msx = maskImg3D->Num_Cols();
	msy = maskImg3D->Num_Rows();
	msz = maskImg3D->Num_Levels();
	mpx = img3D->Get_Pixel_SizeX();
	mpy = img3D->Get_Pixel_SizeY();
	mpz = img3D->Get_Pixel_SizeZ();
	maskA = maskImg3D->GetArray();
	int msxy=msx*msy;

	Vec3 ov, mv;
	IntVec3 iv;
	int x, y, z, k;
	// initialize the mask
	for(k=0; k<msx*msy*msz; k++) 
	{
		if(maskA[k]==organIndex+1) maskA[k]=0;
	}

	short *imgA;
	imgA = img3D->GetArray();

	probSum=0;
	float prob, sqrt_2pi, prob_std2, gray;
	sqrt_2pi = sqrt(6.283);
	prob_std2 = 2*probStd*probStd;

	int pi;
	int count=0;

	for(z=0; z<osz; z+=sampleRate)
	{
		ov.z = orgMask->GetSlicePosition(z);
		for(y=0, ov.y=0; y<osy; y+=sampleRate, ov.y+=opy*sampleRate)
		{
			for(x=0, ov.x=0; x<osx; x+=sampleRate, ov.x+=opx*sampleRate)
			{
				k = z*osxy+y*osx+x;
				if(organA[k]==organIndex+1)
				{
					mv = f*ov;
					iv.x = (int)(mv.x/mpx);
					iv.y = (int)(mv.y/mpy);
					iv.z = img3D->GetSliceLevel(mv.z, 1);
					if(iv.z>=0 && iv.z<msz && iv.x>0 && iv.x<msx && iv.y>0 && iv.y<msy)
					{
						pi = iv.z*msxy+iv.y*msx+iv.x;
						maskA[pi]=organA[k];
						gray = imgA[pi];
						prob = exp(-(gray-probMean)*(gray-probMean)/prob_std2)/(sqrt_2pi*probStd);
						probSum += prob;
						count++;
					}
				}
			}	// for x
		}	// for y
	}	// for z

	if(count>0) probAvg = probSum/(float)count;
}

// add one location model, and intanstiate against the ALM model
// confusing about single location model and ALM, need clarify it later
//
void NIH_BoneSegmentation_Dlg::pushButton_Load_ALM_Model_clicked()
{
	// get the file name
	QString strNewFileName = QFileDialog::getOpenFileName(
		QString::null,
        "Model File (*.xml)",
         this,
         "Load Model",
         "Choose a file" );
 
	if(strNewFileName==QString::null) return;

    QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );
	if(totalALMModels==MaxALMModels) totalALMModels--;

	LoadSingleALMModel(strNewFileName.ascii(),totalALMModels);

	if(singleALM[totalALMModels].smoothedCord->draw_mode==LINE_MODE)
	{
		NIH_OpenGL_CentralizeModel(display_id2, singleALM[totalALMModels].smoothedCord);
	}

	totalALMModels++;

	curALMModel = totalALMModels-1;
	spinBox_ALM_Model->setRange(0, totalALMModels-1);
	spinBox_ALM_Model->setValue(curALMModel);

	if(flagGT)	GTLocationModel = curALMModel;

	// fill the location Vector and compute inverse PDM for the single model
	//
	if(ALM_PDM.pdmRange.GetSize()>0)
	{
		int vecLength=ALM_PDM.vecLength;
		int k;
		Vec3 scale;
///		scale.x = ALM_PDM.referenceModelSize.x/singleALM[curALMModel].ribCageSize.x;
///		scale.y = ALM_PDM.referenceModelSize.y/singleALM[curALMModel].ribCageSize.y;
///		scale.z = ALM_PDM.referenceModelSize.z/singleALM[curALMModel].ribCageSize.z;
		scale.x = scale.y = scale.z = 1;

		singleALM[curALMModel].locationVec.SetSize(vecLength);

		k=0;
		for(int m=0; m<5; m++)
		{
			if(ALM_PDM.organIndexList.GetIndex(m)==-1) continue;
			for(int v=0; v<6; v++)
			{
				if(ALM_PDM.vertebraIndexList.GetIndex(v)==-1) continue;
				singleALM[curALMModel].locationVec[k] = 
					(singleALM[curALMModel].organSurf[m]->centroid.x+singleALM[curALMModel].organSurf[m]->loc.x-singleALM[curALMModel].vertebra_3d_model[v]->loc.x)*scale.x; 
				k++;
				singleALM[curALMModel].locationVec[k] = 
					(singleALM[curALMModel].organSurf[m]->centroid.y+singleALM[curALMModel].organSurf[m]->loc.y-singleALM[curALMModel].vertebra_3d_model[v]->loc.y)*scale.y; 
				k++;
				singleALM[curALMModel].locationVec[k] = 
					(singleALM[curALMModel].organSurf[m]->centroid.z+singleALM[curALMModel].organSurf[m]->loc.z-singleALM[curALMModel].vertebra_3d_model[v]->loc.z)*scale.z;
				k++;
			}
		}	// for m

		ALM_PDM.InversePDM(singleALM[curALMModel].locationVec, singleALM[curALMModel].pdmPara);

		char info[100];
		sprintf(info, "InversePDM: %f, %f, %f, %f\n", singleALM[curALMModel].pdmPara[0], singleALM[curALMModel].pdmPara[1], 
			singleALM[curALMModel].pdmPara[2], singleALM[curALMModel].pdmPara[3]);
		textLabel_info->setText(info);

		// write to a log file
		FILE *fp;
		fp = fopen("c:\\tmp\\ALM_log.txt", "a");
		fprintf(fp, "InversePDM %d: ", curALMModel);
		for(int i=0; i<singleALM[curALMModel].pdmPara.GetSize(); i++)
		{
			fprintf(fp, "%f, ", singleALM[curALMModel].pdmPara[i]);
		}
		fprintf(fp, "\n");
		if(flagGT)
		{
			fprintf(fp, "GT InversePDM %d: ", GTLocationModel);
			for(int i=0; i<singleALM[curALMModel].pdmPara.GetSize(); i++)
			{
				fprintf(fp, "%f, ", singleALM[curALMModel].pdmPara[i]);
			}
			fprintf(fp, "\n");
		}

		fclose(fp);
	}

///	if(GTLocationModel!=-1 && refLocationModel!=-1) AlignTwoLocationModelsUsingSpine(refLocationModel, GTLocationModel);

	NIH_OpenGL_Refresh_Display(display_id2);
    QApplication::restoreOverrideCursor();
}

// load multiple location model, and construct the ALM
//
void NIH_BoneSegmentation_Dlg::pushButton_Load_ALM_Multi_Model_clicked()
{
	// get the file name
	QString strNewFileName = QFileDialog::getOpenFileName(
		QString::null,
        "Txt File (*.txt)",
         this,
         "Load Multiple Models",
         "Choose a file" );
 

	if(strNewFileName==QString::null) return;

    QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );

	Load_ALM_Multi_Model_File(strNewFileName.ascii());
}

void NIH_BoneSegmentation_Dlg::Load_ALM_Multi_Model_File(const char *model_fn)
{
	FILE *fp = fopen(model_fn,"r");
	
	if(fp == NULL)
	{
		fprintf(stderr,"Can't open multiple model file %s\n",model_fn);
		return;
	}

	totalALMModels=0;
	
	int numModels, i, k;
	char line[200];
	fgets(line, 200, fp);
	sscanf(line, "%d", &numModels);

	int count=0;
	QString caseName;
	bool hasRef=false;
	for(i=0; i<numModels; i++)
	{
		fgets(line, 200, fp);
		caseName = QString(line).simplifyWhiteSpace().right(8);
		caseName=caseName.left(4);

		if(global_ALM_multiple_exclude.length()>2!=-1 && caseName==global_ALM_multiple_exclude) continue;

		if(global_ALM_multiple_reference.length()<2)
		{
			LoadSingleALMModel(QString(line).simplifyWhiteSpace(), count);
			printf("%d: %s\n", count, QString(line).simplifyWhiteSpace().ascii());
			if(count==0) 
			{
				global_ALM_multiple_reference=caseName;
				hasRef = true;
			}
		}
		else
		{
			if(caseName==global_ALM_multiple_reference) 
			{
				hasRef=true;
				LoadSingleALMModel(QString(line).simplifyWhiteSpace(), 0);
				printf("%d: %s\n", 0, QString(line).simplifyWhiteSpace().ascii());
			}
			else if(!hasRef) 
			{
				LoadSingleALMModel(QString(line).simplifyWhiteSpace(), count+1);
				printf("%d: %s\n", count+1, QString(line).simplifyWhiteSpace().ascii());
			}
			else 
			{
				LoadSingleALMModel(QString(line).simplifyWhiteSpace(), count);
				printf("%d: %s\n", count, QString(line).simplifyWhiteSpace().ascii());
			}
		}
		count++;
	}
	fclose(fp);

	totalALMModels = count;

	// align the models with model 0
	for(i=1; i<totalALMModels; i++)
	{
		AlignTwoLocationModelsUsingSpine(i, 0);
	}	// for i

	curALMModel = 0;
	spinBox_ALM_Model->setRange(0, totalALMModels-1);
	spinBox_ALM_Model->setValue(0);

	// compute statistics to build ALM
	strcpy(ALM_PDM.referenceModelFileName, singleALM[0].modelFileName);

	// first find out all the common vertebrae and commonOrgans
	int vertInModels[6];
	int commonVerts=0;
	ALM_PDM.vertebraIndexList.SetSize(0);
	for(int v=0; v<6; v++)
	{
		vertInModels[v] = 0;
		for(i=0; i<totalALMModels; i++)
		{
			if(singleALM[i].vertebra_3d_model[v]->sizeX>0) vertInModels[v]++;
		}

		if(vertInModels[v]==totalALMModels) 
		{
			ALM_PDM.vertebraIndexList.Add(v);
			commonVerts++;
		}
	}

	int organInModels[5];
	int commonOrgans=0;
	ALM_PDM.organIndexList.SetSize(0);
	for(int m=0; m<5; m++)
	{
		organInModels[m] = 0;
		for(i=0; i<totalALMModels; i++)
		{
			if(singleALM[i].organSurf[m]->sizeX>0) organInModels[m]++;
		}

		if(organInModels[m]==totalALMModels) 
		{
			ALM_PDM.organIndexList.Add(m);
			commonOrgans++;
		}
	}

	int vecLength = 0;

///	vecLength += commonVerts*10;		// each vertebra has 9 parameters, 3 locations, 4 querternion, and 3 sizes, may not need to be stored since it is extracted

///	vecLength += commonOrgans*9;	// each organ has 9 parameters, 3 locations, 3 orientation, and 3 sizes, maybe shouldn't use location

	vecLength += commonOrgans*commonVerts*3;	// relative locations between each vert and each organ, just location

	Vec3 angles;
	Quaternion q;
	for(i=0; i<totalALMModels; i++)
	{
		singleALM[i].locationVec.SetSize(vecLength);
		k=0;

/*		// vertebra location
		for(int v=0; v<6; v++)
		{
			if(vertInModels[v]==numModels)
			{
				singleALM[i].locationVec[k] = singleALM[i].vertebra_3d_model[v]->loc.x; k++;
				singleALM[i].locationVec[k] = singleALM[i].vertebra_3d_model[v]->loc.y; k++;
				singleALM[i].locationVec[k] = singleALM[i].vertebra_3d_model[v]->loc.z; k++;

				q = Quaternion(Rotation(singleALM[i].vertebra_3d_model[v]->dirX, 
					singleALM[i].vertebra_3d_model[v]->dirY, singleALM[i].vertebra_3d_model[v]->dirZ));
				singleALM[i].locationVec[k] = q.v.x; k++;
				singleALM[i].locationVec[k] = q.v.y; k++;
				singleALM[i].locationVec[k] = q.v.z; k++;
				singleALM[i].locationVec[k] = q.w; k++;

				singleALM[i].locationVec[k] = singleALM[i].vertebra_3d_model[v]->sizeX; k++;
				singleALM[i].locationVec[k] = singleALM[i].vertebra_3d_model[v]->sizeY; k++;
				singleALM[i].locationVec[k] = singleALM[i].vertebra_3d_model[v]->sizeZ; k++;
			}
		}	// for v
*/
		// organ location
/*		for(int m=0; m<5; m++)
		{
			if(organInModels[m]==numModels)
			{
				singleALM[i].locationVec[k] = singleALM[i].organSurf[m]->loc.x; k++;
				singleALM[i].locationVec[k] = singleALM[i].organSurf[m]->loc.y; k++;
				singleALM[i].locationVec[k] = singleALM[i].organSurf[m]->loc.z; k++;

				angles = Rotation(singleALM[i].organSurf[m]->dirX, 
					singleALM[i].organSurf[m]->dirY, singleALM[i].organSurf[m]->dirZ).Angles(X_Y_Z);
				singleALM[i].locationVec[k] = angles.x; k++;
				singleALM[i].locationVec[k] = angles.y; k++;
				singleALM[i].locationVec[k] = angles.z; k++;

				singleALM[i].locationVec[k] = singleALM[i].organSurf[m]->sizeX; k++;
				singleALM[i].locationVec[k] = singleALM[i].organSurf[m]->sizeY; k++;
				singleALM[i].locationVec[k] = singleALM[i].organSurf[m]->sizeZ; k++;
			}
		}
*/
		// relative vertebra-organ location
		for(int m=0; m<5; m++)
		{
			if(organInModels[m]!=totalALMModels) continue;
			for(int v=0; v<6; v++)
			{
				if(vertInModels[v]!=totalALMModels) continue;
				singleALM[i].locationVec[k] = singleALM[i].organSurf[m]->centroid.x+singleALM[i].organSurf[m]->loc.x-singleALM[i].vertebra_3d_model[v]->loc.x; k++;
				singleALM[i].locationVec[k] = singleALM[i].organSurf[m]->centroid.y+singleALM[i].organSurf[m]->loc.y-singleALM[i].vertebra_3d_model[v]->loc.y; k++;
				singleALM[i].locationVec[k] = singleALM[i].organSurf[m]->centroid.z+singleALM[i].organSurf[m]->loc.z-singleALM[i].vertebra_3d_model[v]->loc.z; k++;
			}
		}

		ALM_PDM.AddTrainingVec(singleALM[i].locationVec);
	}	// for i

	ALM_PDM.pdmDegree = totalALMModels;
	ALM_PDM.TrainPDM();

	ALM_PDM.referenceModelSize = singleALM[0].ribCageSize;

	// compute the invertPDM
	for(i=0; i<totalALMModels; i++)
	{
		ALM_PDM.InversePDM(singleALM[i].locationVec, singleALM[i].pdmPara);
	}


	if(singleALM[0].smoothedCord->draw_mode==LINE_MODE)
	{
		NIH_OpenGL_CentralizeModel(display_id2, singleALM[0].smoothedCord);
	}

	// hide some models for clarity
/*	for(i=0; i<totalALMModels; i++)
	{
		for(int m=0; m<5; m++)
		{
			singleALM[i].organLocationModel[m]->draw_mode = HIDE_MODE;
			singleALM[i].organSurf[m]->draw_mode = HIDE_MODE;
		}

		for(int v=0; v<6; v++)
		{
			singleALM[i].vertebra_3d_model[v]->draw_mode= HIDE_MODE;
		}
	}
*/

	for(i=0; i<ALM_PDM.pdmDegree; i++) ALM_PDM.pdmValue[i]=0;

	spinBox_ALM_Mode->setRange(0, ALM_PDM.pdmDegree-1);
	spinBox_ALM_Mode->setValue(0);

	textLabel1_ALM_Mode_cur->setText(QString::number(ALM_PDM.pdmValue[0]));
	textLabel1_ALM_Mode_range->setText(QString::number(ALM_PDM.pdmRange[0]));

	slider_ALM_Mode->setRange(1, 200);
	slider_ALM_Mode->setValue(100);

	NIH_OpenGL_Refresh_Display(display_id2);

	// screen capture, for documentation
	char ss_fn[200];
	sprintf(ss_fn, "c:\\tmp\\ALM.jpg");
	NIH_OpenGL_Save_Display(display_id2, ss_fn);
    QApplication::restoreOverrideCursor();

}

// assign the current location model with the mean model
//
void NIH_BoneSegmentation_Dlg::pushButton_Mean_ALM_clicked()
{
	if(totalALMModels<=0 || ALM_PDM.pdmRange.GetSize()==0 || curALMModel>=totalALMModels || curALMModel<0) return;

    QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );
	int i, v, m, k, j;

	// Hide all other models
	for(i=0; i<totalALMModels; i++)
	{
		if(i==curALMModel) continue;
		for(m=0; m<5; m++)
		{
			singleALM[i].organLocationModel[m]->draw_mode = HIDE_MODE;
			singleALM[i].organSurf[m]->draw_mode = HIDE_MODE;
		}

		for(v=0; v<6; v++)
		{
			singleALM[i].vertebra_3d_model[v]->draw_mode= HIDE_MODE;
		}

		singleALM[i].smoothedCord->draw_mode = HIDE_MODE;
	}

	// replace the model location with the mean PDM model
	singleALM[curALMModel].locationVec = ALM_PDM.meanVec;

	IntanstiateSingleALM(singleALM[curALMModel], ALM_PDM.meanVec);

	for(i=0; i<ALM_PDM.pdmDegree; i++) ALM_PDM.pdmValue[i]=0;

	spinBox_ALM_Mode->setRange(0, ALM_PDM.pdmDegree-2);
	spinBox_ALM_Mode->setValue(0);

	textLabel1_ALM_Mode_cur->setText(QString::number(ALM_PDM.pdmValue[0]));
	textLabel1_ALM_Mode_range->setText(QString::number(ALM_PDM.pdmRange[0]));

	slider_ALM_Mode->setRange(1, 200);
	slider_ALM_Mode->setValue(100);

	NIH_OpenGL_Refresh_Display(display_id2);
	pushButton_refresh_clicked();
    QApplication::restoreOverrideCursor();
	return;
}

// intansiate a model based on the model parameter
//
void NIH_BoneSegmentation_Dlg::pushButton_Intanstiate_ALM_clicked()
{
	if(totalALMModels<=0 || ALM_PDM.pdmRange.GetSize()==0 || curALMModel>=totalALMModels || curALMModel<0) return;

    QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );
	int i, v, m, k, j;

	// replace the model location with the intanstiate PDM model
	singleALM[curALMModel].locationVec = ALM_PDM.InstantiatePDM(singleALM[curALMModel].pdmPara);

	IntanstiateSingleALM(singleALM[curALMModel], singleALM[curALMModel].locationVec);

	NIH_OpenGL_Refresh_Display(display_id2);
	pushButton_refresh_clicked();
    QApplication::restoreOverrideCursor();
	return;
}


int NIH_BoneSegmentation_Dlg::IntanstiateSingleALM(NIH_ALM_Single_Model &sALM, CIS_Vector_JY_double &locationVec)
{
	int i, v, m, k, j;

	Vec3 scale;
///	scale.x = sALM.ribCageSize.x/ALM_PDM.referenceModelSize.x;
///	scale.y = sALM.ribCageSize.y/ALM_PDM.referenceModelSize.y;
///	scale.z = sALM.ribCageSize.z/ALM_PDM.referenceModelSize.z;
	scale.x = scale.y = scale.z = 1;

	k=0;

	// set vertebra location
	Rotation rot;
	Quaternion q;
/*	for(i=0; i<ALM_PDM.vertebraIndexList.GetSize(); i++)
	{
		v = ALM_PDM.vertebraIndexList[i];
		sALM.vertebra_3d_model[v]->loc.x = locationVec[k]; k++;
		sALM.vertebra_3d_model[v]->loc.y = locationVec[k]; k++;
		sALM.vertebra_3d_model[v]->loc.z = locationVec[k]; k++;

		q = Quaternion(locationVec[k], locationVec[k+1], locationVec[k+2], locationVec[k+3]);
		rot = q.Compose();
		sALM.vertebra_3d_model[v]->dirX = rot.Rx;
		sALM.vertebra_3d_model[v]->dirY = rot.Ry;
		sALM.vertebra_3d_model[v]->dirZ = rot.Rz;
		k+=4;

		sALM.vertebra_3d_model[v]->sizeX = locationVec[k]; k++;
		sALM.vertebra_3d_model[v]->sizeY = locationVec[k]; k++;
		sALM.vertebra_3d_model[v]->sizeZ = locationVec[k]; k++;

		sALM.vertebra_3d_model[v]->radius = (sALM.vertebra_3d_model[v]->sizeX+sALM.vertebra_3d_model[v]->sizeY)/8;
		sALM.vertebra_3d_model[v]->height = sALM.vertebra_3d_model[v]->sizeZ*0.9;
		sALM.vertebra_3d_model[v]->ComputeFrame();
	}
*/

	Vec3 relLoc;
	// set the relative organ location
	total_probSum = total_probAvg = 0;
	for(i=0; i<ALM_PDM.organIndexList.GetSize(); i++)
	{
		m = ALM_PDM.organIndexList[i];
		relLoc = Vec3(0,0,0);
		for(j=0; j<ALM_PDM.vertebraIndexList.GetSize(); j++)
		{
			v = ALM_PDM.vertebraIndexList[j];
			relLoc += (sALM.vertebra_3d_model[v]->loc+Vec3(locationVec[k]*scale.x, locationVec[k+1]*scale.y, locationVec[k+2]*scale.z));

			k+=3;
		}

		relLoc /= (double)ALM_PDM.vertebraIndexList.GetSize();

		sALM.organSurf[m]->loc = relLoc-sALM.organSurf[m]->centroid;
		sALM.organSurf[m]->ComputeFrame();

		sALM.organLocationModel[m]->pt_list1.SetSize(ALM_PDM.vertebraIndexList.GetSize());
		sALM.organLocationModel[m]->pt_list2.SetSize(ALM_PDM.vertebraIndexList.GetSize());

		for(j=0; j<ALM_PDM.vertebraIndexList.GetSize(); j++)
		{
			v = ALM_PDM.vertebraIndexList[j];
			sALM.organLocationModel[m]->pt_list1[j] = sALM.vertebra_3d_model[v]->loc;
			sALM.organLocationModel[m]->pt_list2[j] = sALM.organSurf[m]->centroid+sALM.organSurf[m]->loc;
		}

		// remap to mask
		if(organMask && maskImg3D)
		{
			Frame f;
			f = sALM.organSurf[m]->GetFrame();
			MapOrganMask(organMask, f, m, 1, sALM.lfw.organ[m]->meanIntensity, sALM.lfw.organ[m]->varianceIntensity, sALM.lfw.organ[m]->probSum, sALM.lfw.organ[m]->probAvg);
			// print out for record
			printf("Organ %d: %f %f\n", m, sALM.lfw.organ[m]->probSum, sALM.lfw.organ[m]->probAvg);
			total_probSum += sALM.lfw.organ[m]->probSum;
			total_probAvg += sALM.lfw.organ[m]->probAvg;
		}
	}	// for i
	printf("Total: %f %f\n", total_probSum, total_probAvg);

	return 0;
}

void NIH_BoneSegmentation_Dlg::pushButton_ALM_Refinement_clicked()
{
	if(organMask==NULL || maskImg3D==NULL) return;

    QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );

	total_probSum = total_probAvg = 0;
	for(int m=0; m<5; m++)
	{
		Vec3 originalLoc, move, maxMove;
		double dx, dy, dz;
		double range=2;
		float maxProbSum=0, maxProbAvg=0;
		Frame f;

		originalLoc = singleALM[0].organSurf[m]->loc;
		for(dz=-range; dz<=range; dz++)
		{
			for(dy=-range; dy<=range; dy++)
			{
				for(dx=-range; dx<=range; dx++)
				{
					move = Vec3(dx, dy, dz);
					singleALM[0].organSurf[m]->loc = originalLoc+move;
					singleALM[0].organSurf[m]->ComputeFrame();

					f = singleALM[0].organSurf[m]->GetFrame();
					MapOrganMask(organMask, f, m, 1, singleALM[0].lfw.organ[m]->meanIntensity, 
						singleALM[0].lfw.organ[m]->varianceIntensity, singleALM[0].lfw.organ[m]->probSum, 
						singleALM[0].lfw.organ[m]->probAvg);
					if(singleALM[0].lfw.organ[m]->probAvg>maxProbAvg)
					{
						maxProbSum = singleALM[0].lfw.organ[m]->probSum;
						maxProbAvg = singleALM[0].lfw.organ[m]->probAvg;
						maxMove = move;
					}
				}
			}
		}
		singleALM[0].organSurf[m]->loc = originalLoc+maxMove;
		singleALM[0].organSurf[m]->ComputeFrame();
		singleALM[0].lfw.organ[m]->probSum = maxProbSum;
		singleALM[0].lfw.organ[m]->probAvg = maxProbAvg;
		total_probSum += singleALM[0].lfw.organ[m]->probSum;
		total_probAvg += singleALM[0].lfw.organ[m]->probAvg;
	}	// for m

	NIH_OpenGL_Refresh_Display(display_id2);
	pushButton_refresh_clicked();
    QApplication::restoreOverrideCursor();
	return;
}


void NIH_BoneSegmentation_Dlg::slider_ALM_Mode_valueChanged( int )
{
	if(totalALMModels<=0 || ALM_PDM.pdmRange.GetSize()==0 || curALMModel>=totalALMModels || curALMModel<0) return;

	int pos = slider_ALM_Mode->value()-100;
	int curMode = spinBox_ALM_Mode->value();

    QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );
	double value = (double)pos/100.0*ALM_PDM.pdmRange[curMode]*2;
	ALM_PDM.pdmValue[curMode] = value;

	singleALM[curALMModel].locationVec = ALM_PDM.InstantiatePDM();

	IntanstiateSingleALM(singleALM[curALMModel], singleALM[curALMModel].locationVec);

	textLabel1_ALM_Mode_cur->setText(QString::number(ALM_PDM.pdmValue[curMode]));
	NIH_OpenGL_Refresh_Display(display_id2);
	pushButton_refresh_clicked();
    QApplication::restoreOverrideCursor();
	return;
}

void NIH_BoneSegmentation_Dlg::spinBox_ALM_Mode_valueChanged( int )
{
	if(totalALMModels<=0 || ALM_PDM.pdmRange.GetSize()==0) return;

	int curMode = spinBox_ALM_Mode->value();
	textLabel1_ALM_Mode_cur->setText(QString::number(ALM_PDM.pdmValue[curMode], 'g', 1));
	textLabel1_ALM_Mode_range->setText(QString::number(ALM_PDM.pdmRange[curMode], 'g', 1));

	int pos=(int)(ALM_PDM.pdmValue[curMode]/ALM_PDM.pdmRange[curMode]/2*100);

	pos += 100;

	slider_ALM_Mode->setValue(pos);

	return;
}

void NIH_BoneSegmentation_Dlg::spinBox_ALM_Model_valueChanged( int )
{
	if(totalALMModels<=0) return;

	curALMModel = spinBox_ALM_Model->value();

	return;
}

// switch the model display display mode among HIDE, SURFACE, and SLICE
//
void NIH_BoneSegmentation_Dlg::SwitchALMDisplayMode(int modelIndex)
{
	if(modelIndex<0 || modelIndex>=totalALMModels) return;

	NIH_ALM_Single_Model *sModel;
	sModel = &(singleALM[modelIndex]);
	if(sModel->organSurf[0]!=NULL)
	{
		if(sModel->organSurf[0]->draw_mode==HIDE_MODE)
		{
			for(int m=0; m<5; m++)
			{
				if(sModel->organSurf[m]) 
				{
					sModel->organSurf[m]->draw_mode = SURFACE_MODE;
					sModel->organSurf[m]->lineWidth = 1.0;
					sModel->organSurf[m]->SetCurrAlpha(0.8);
					sModel->organLocationModel[m]->draw_mode = LINE_MODE;
				}
			}	// for m
			for(int v=0; v<6; v++)
			{
				if(sModel->vertebra_3d_model[v]) sModel->vertebra_3d_model[v]->draw_mode = LINE_MODE;
			}	// for v
			sModel->smoothedCord->draw_mode = LINE_MODE;
		}	// if sModel
		else if(sModel->organSurf[0]->draw_mode==SURFACE_MODE)
		{
			for(int m=0; m<5; m++)
			{
				if(sModel->organSurf[m]) 
				{
					sModel->organSurf[m]->draw_mode = LINE_MODE | SLICE_MODE;
					sModel->organSurf[m]->lineWidth = 2.0;
					sModel->organSurf[m]->SetCurrAlpha(1);
					sModel->organLocationModel[m]->draw_mode = HIDE_MODE;
				}
			}	// for m
			for(int v=0; v<6; v++)
			{
				if(sModel->vertebra_3d_model[v]) sModel->vertebra_3d_model[v]->draw_mode = HIDE_MODE;
			}	// for v
			sModel->smoothedCord->draw_mode = HIDE_MODE;
		}	// else sModel
		else if(sModel->organSurf[0]->draw_mode==(LINE_MODE|SLICE_MODE))
		{
			for(int m=0; m<5; m++)
			{
				if(sModel->organSurf[m]) 
				{
					sModel->organSurf[m]->draw_mode = HIDE_MODE;
					sModel->organLocationModel[m]->draw_mode = HIDE_MODE;
				}
			}	// for m
			for(int v=0; v<6; v++)
			{
				if(sModel->vertebra_3d_model[v]) sModel->vertebra_3d_model[v]->draw_mode = HIDE_MODE;
			}	// for v
			sModel->smoothedCord->draw_mode = HIDE_MODE;
		}	// else sModel
	}	// if sModel

	NIH_OpenGL_Refresh_Display(display_id);
	NIH_OpenGL_Refresh_Display(display_id2);
	return;
}

// load the ALM from xml file, not just a single location model
//
void NIH_BoneSegmentation_Dlg::pushButton_Load_ALM_clicked()
{
	// get the file name
	QString fileExt;

	QString strNewFileName = QFileDialog::getOpenFileName(
		QString::null,
        "Model File (*.xml)",
         this,
         "Load Model",
         "Choose a file" );
 

    QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );
	if(strNewFileName!=QString::null) 
	{
		ALM_PDM.LoadALMModel(strNewFileName.ascii());

		LoadSingleALMModel(ALM_PDM.referenceModelFileName, 0);
		refLocationModel = 0;
		totalALMModels =1;
		curALMModel = 0;
		NIH_OpenGL_CentralizeModel(display_id2, singleALM[0].smoothedCord);
		NIH_OpenGL_Refresh_Display(display_id2);
	}
    QApplication::restoreOverrideCursor();
	return;
}

// save the ALM into a xml file
//
void NIH_BoneSegmentation_Dlg::pushButton_Save_ALM_clicked()
{
	// get the file name
	QString fileExt;

	QString strNewFileName = QFileDialog::getSaveFileName(
		QString::null,
        "Model File (*.xml)",
         this,
         "Save Model",
         "Choose a file" );
 

	if(strNewFileName!=QString::null) 
	{
		ALM_PDM.SaveALMModel(strNewFileName.ascii());
	}
	return;
}
void NIH_BoneSegmentation_Dlg::pushButton_ALM_Batch_clicked()
{
	// get the file name
	QString strNewFileName = QFileDialog::getOpenFileName(
		QString::null,
        "Txt File (*.txt)",
         this,
         "Load Batch ",
         "Choose a file" );
 
	if(strNewFileName==QString::null) return;
	ALM_Batch_file(strNewFileName.ascii());
	return;
}

void NIH_BoneSegmentation_Dlg::ALM_Batch_file(const char *batchFn)
{
	FILE *fp = fopen(batchFn,"r");
	
	if(fp == NULL)
	{
		fprintf(stderr,"Can't open batch file %s\n",batchFn);
		return;
	}

	QString projFn = "dod_atlas.prj";
	QString ALMFn = "ALM12_r132f.xml";
	int studyID=0;

	char line[200];
	QStringList list;
	int count=-1;
	while(!feof(fp))
	{
		strcpy(line,"");
		fgets(line, 200, fp);
		if(strlen(line)<2) continue;
		if(line[0]=='-' && line[1]=='1') break;

		if(strlen(line)>2 && line[0]!='#')
		{
			count++;
			if(global_ALM_batch_id!=-1 && count!= global_ALM_batch_id) continue;

			list = QStringList::split(',', line, false);
			if(list.size()!=3) continue;

			projFn = list[0].simplifyWhiteSpace();
			ALMFn = list[1].simplifyWhiteSpace();
			studyID = list[2].toInt();

			printf("count %d: %s, %s, %d\n", count, projFn.ascii(), ALMFn.ascii(), studyID);
			ALMBatchRunning(projFn, ALMFn, studyID);
		}
	}

	fclose(fp);
}

void NIH_BoneSegmentation_Dlg::ALMBatchRunning(QString &projFn, QString& ALMFn, int studyID)
{
	FILE *fp;
	fp = fopen("c:\\tmp\\ALM_log.txt", "a");

	QString GTFn;

	QString rootFn, projFn1, ALMFn1, GTFn1;
	char screenshot_root[200];
	char ss_fn[400];
	QString ALMFn2;
	ALMFn2 = ALMFn;
	ALMFn2.replace(".xml", "");

	int z_slice=50;
	int x_slice=120;
	int y_slice=260;
	

	rootFn = "C:\\Experiments\\bone_training\\Atlas\\";
	projFn1 = rootFn+projFn;
	ALMFn1 = rootFn+ALMFn;

	sprintf(screenshot_root, "%sscreenshot\\", rootFn.ascii());

	// Load project file
	load_project(projFn1.ascii());

	// Load study
	curStudy = studyID;
	LoadStudy(curStudy);
	GTFn=studyEntries[curStudy].patientName+".xml";

	totalALMModels=0;

	// log the file names
	fprintf(fp, "%s, %s, %s, %d\n", projFn.ascii(), ALMFn.ascii(), GTFn.ascii(), curStudy);

	// Load ALM
	flagGT=false;
	ALM_PDM.LoadALMModel(ALMFn1.ascii());
	LoadSingleALMModel(ALM_PDM.referenceModelFileName, 0);
	refLocationModel = 0;
	totalALMModels =1;

	// output the probSum of reference model
	fprintf(fp, "Ref Prob,");
	for(int m=0; m<5; m++)
	{
		fprintf(fp, "%f, %f, ", singleALM[0].lfw.organ[m]->probSum, singleALM[0].lfw.organ[m]->probAvg);
	}
	fprintf(fp, "%f, %f,\n", total_probSum, total_probAvg);

	// Load ground truth model (not align them)
	flagGT=true;
	GTFn1 = rootFn+GTFn;
	LoadSingleALMModel(GTFn1.ascii(),1);
	GTLocationModel = 1;
	totalALMModels=2;

	// output the probSum of GT model
	fprintf(fp, "GT Prob,");
	for(int m=0; m<5; m++)
	{
		fprintf(fp, "%f, %f, ", singleALM[1].lfw.organ[m]->probSum, singleALM[1].lfw.organ[m]->probAvg);
	}
	fprintf(fp, "%f, %f,\n", total_probSum, total_probAvg);

	flagGT=false;

	// calculate distance between reference and GT
	Vec3 delta;
	double avgDist, dist;
	fprintf(fp, "Ref-GT dist,");
	avgDist = 0;
	for(int m=0; m<5; m++)
	{
		delta = (singleALM[1].organSurf[m]->loc+singleALM[1].organSurf[m]->centroid)-
			(singleALM[0].organSurf[m]->loc+singleALM[0].organSurf[m]->centroid);
		dist = delta.len();
		avgDist += dist;
		fprintf(fp, "%f, ", dist);
	}
	avgDist /= (double)5;
	fprintf(fp, "%f,\n", avgDist);

	// load mean model
	curALMModel = refLocationModel;
	singleALM[curALMModel].locationVec = ALM_PDM.meanVec;
	IntanstiateSingleALM(singleALM[curALMModel], ALM_PDM.meanVec);

	// output the probSum of GT model
	fprintf(fp, "Mean Ref Prob,");
	for(int m=0; m<5; m++)
	{
		fprintf(fp, "%f, %f, ", singleALM[0].lfw.organ[m]->probSum, singleALM[0].lfw.organ[m]->probAvg);
	}
	fprintf(fp, "%f, %f,\n", total_probSum, total_probAvg);

	// calculate distance between reference and GT
	fprintf(fp, "Mean Ref-GT dist,");
	avgDist = 0;
	for(int m=0; m<5; m++)
	{
		delta = (singleALM[1].organSurf[m]->loc+singleALM[1].organSurf[m]->centroid)-
			(singleALM[0].organSurf[m]->loc+singleALM[0].organSurf[m]->centroid);
		dist = delta.len();
		avgDist += dist;
		fprintf(fp, "%f, ", dist);
	}
	avgDist /= (double)5;
	fprintf(fp, "%f,\n", avgDist);

	NIH_OpenGL_Show_Background(display_id, false);
	ChangeImageSlice_x(x_slice);
	sprintf(ss_fn, "%s//%s_mean_s%d.jpg", screenshot_root, ALMFn2.ascii(), x_slice);
	NIH_OpenGL_Save_Display(display_id, ss_fn);

	ChangeImageSlice_y(y_slice);
	sprintf(ss_fn, "%s//%s_mean_c%d.jpg", screenshot_root, ALMFn2.ascii(), y_slice);
	NIH_OpenGL_Save_Display(display_id, ss_fn);

	NIH_OpenGL_Show_Background(display_id, true);
	ChangeImageSlice(z_slice);
	sprintf(ss_fn, "%s//%s_mean_a%d.jpg", screenshot_root, ALMFn2.ascii(), z_slice);
	NIH_OpenGL_Save_Display(display_id, ss_fn);

	// align reference and GT
	if(GTLocationModel!=-1 && refLocationModel!=-1) AlignTwoLocationModelsUsingSpine(refLocationModel, GTLocationModel);

	// output the probSum of GT model
	fprintf(fp, "Aligned Mean Ref Prob,");
	for(int m=0; m<5; m++)
	{
		fprintf(fp, "%f, %f, ", singleALM[0].lfw.organ[m]->probSum, singleALM[0].lfw.organ[m]->probAvg);
	}
	fprintf(fp, "%f, %f,\n", total_probSum, total_probAvg);

	// calculate distance between reference and GT
	fprintf(fp, "Aligned Mean Ref-GT dist,");
	avgDist = 0;
	for(int m=0; m<5; m++)
	{
		delta = (singleALM[1].organSurf[m]->loc+singleALM[1].organSurf[m]->centroid)-
			(singleALM[0].organSurf[m]->loc+singleALM[0].organSurf[m]->centroid);
		dist = delta.len();
		avgDist += dist;
		fprintf(fp, "%f, ", dist);
	}
	avgDist /= (double)5;
	fprintf(fp, "%f,\n", avgDist);

	NIH_OpenGL_Show_Background(display_id, false);
	ChangeImageSlice_x(x_slice);
	sprintf(ss_fn, "%s//%s_align_s%d.jpg", screenshot_root, ALMFn2.ascii(), x_slice);
	NIH_OpenGL_Save_Display(display_id, ss_fn);

	ChangeImageSlice_y(y_slice);
	sprintf(ss_fn, "%s//%s_align_c%d.jpg", screenshot_root, ALMFn2.ascii(), y_slice);
	NIH_OpenGL_Save_Display(display_id, ss_fn);

	NIH_OpenGL_Show_Background(display_id, true);
	ChangeImageSlice(z_slice);
	sprintf(ss_fn, "%s//%s_align_a%d.jpg", screenshot_root, ALMFn2.ascii(), z_slice);
	NIH_OpenGL_Save_Display(display_id, ss_fn);

	// instantiate the GT location using stored PDM value
	singleALM[0].locationVec = ALM_PDM.InstantiatePDM(singleALM[GTLocationModel].pdmPara);

	IntanstiateSingleALM(singleALM[0], singleALM[0].locationVec);

	// output the probSum of GT model
	fprintf(fp, "Instantiate Ref Prob,");
	for(int m=0; m<5; m++)
	{
		fprintf(fp, "%f, %f, ", singleALM[0].lfw.organ[m]->probSum, singleALM[0].lfw.organ[m]->probAvg);
	}
	fprintf(fp, "%f, %f,\n", total_probSum, total_probAvg);

	// calculate distance between reference and GT
	fprintf(fp, "Instantiate Ref-GT dist,");
	avgDist = 0;
	for(int m=0; m<5; m++)
	{
		delta = (singleALM[1].organSurf[m]->loc+singleALM[1].organSurf[m]->centroid)-
			(singleALM[0].organSurf[m]->loc+singleALM[0].organSurf[m]->centroid);
		dist = delta.len();
		avgDist += dist;
		fprintf(fp, "%f, ", dist);
	}
	avgDist /= (double)5;
	fprintf(fp, "%f,\n", avgDist);

	NIH_OpenGL_Show_Background(display_id, false);
	ChangeImageSlice_x(x_slice);
	sprintf(ss_fn, "%s//%s_inst_s%d.jpg", screenshot_root, ALMFn2.ascii(), x_slice);
	NIH_OpenGL_Save_Display(display_id, ss_fn);

	ChangeImageSlice_y(y_slice);
	sprintf(ss_fn, "%s//%s_inst_c%d.jpg", screenshot_root, ALMFn2.ascii(), y_slice);
	NIH_OpenGL_Save_Display(display_id, ss_fn);

	NIH_OpenGL_Show_Background(display_id, true);
	ChangeImageSlice(z_slice);
	sprintf(ss_fn, "%s//%s_inst_a%d.jpg", screenshot_root, ALMFn2.ascii(), z_slice);
	NIH_OpenGL_Save_Display(display_id, ss_fn);


	// local refinement
	pushButton_ALM_Refinement_clicked();

	// output the probSum of GT model
	fprintf(fp, "Refine,");
	for(int m=0; m<5; m++)
	{
		fprintf(fp, "%f, %f, ", singleALM[0].lfw.organ[m]->probSum, singleALM[0].lfw.organ[m]->probAvg);
	}
	fprintf(fp, "%f, %f,\n", total_probSum, total_probAvg);

	// calculate distance between reference and GT
	fprintf(fp, "Refine-GT dist,");
	avgDist = 0;
	for(int m=0; m<5; m++)
	{
		delta = (singleALM[1].organSurf[m]->loc+singleALM[1].organSurf[m]->centroid)-
			(singleALM[0].organSurf[m]->loc+singleALM[0].organSurf[m]->centroid);
		dist = delta.len();
		avgDist += dist;
		fprintf(fp, "%f, ", dist);
	}
	avgDist /= (double)5;
	fprintf(fp, "%f,\n", avgDist);

	NIH_OpenGL_Show_Background(display_id, false);
	ChangeImageSlice_x(x_slice);
	sprintf(ss_fn, "%s//%s_ref_s%d.jpg", screenshot_root, ALMFn2.ascii(), x_slice);
	NIH_OpenGL_Save_Display(display_id, ss_fn);

	ChangeImageSlice_y(y_slice);
	sprintf(ss_fn, "%s//%s_ref_c%d.jpg", screenshot_root, ALMFn2.ascii(), y_slice);
	NIH_OpenGL_Save_Display(display_id, ss_fn);

	NIH_OpenGL_Show_Background(display_id, true);
	ChangeImageSlice(z_slice);
	sprintf(ss_fn, "%s//%s_ref_a%d.jpg", screenshot_root, ALMFn2.ascii(), z_slice);
	NIH_OpenGL_Save_Display(display_id, ss_fn);


	fclose(fp);
	return;
}