
//----------------------------------------------------------------
// Jianhua Yao
// DRD/CC/NIH
//
//
//----------------------------------------------------------------
// Modification History:
//
//  [<RevCode>]: mm/yy - <Description>
//
//----------------------------------------------------------------
// Description
//
// User interface for lung segmentation
// 
//----------------------------------------------------------------

#ifndef _NIH_BONE_SEGMENTATION_DIALOG_H_
#define _NIH_BONE_SEGMENTATION_DIALOG_H_

#include <CIS_2D_Model.h>
#include <CIS_Array_Image3D.h>
#include <CIS_3D_Model_Basic.h>
#include <CIS_2D_Model_Basic.h>
#include <CIS_2D_Model_Advance.h>
#include <CIS_3D_Model_Mesh.h>
#include <CIS_3D_Model_Curve.h>
#include <CIS_2D_ROI.h>
#include <qdir.h>
#include <NIH_BoneSegmentation_Dlg_Base.h>
#include <Polyp_Painting_Structure.h>
#include "NIH_SpineSegmentation_DataStructure.h"
#include <SvmCommittee.h>

#include <NIH_LiveWire.h>
#include <CIS_Algo_Contour.h>

#include "NIH_ActiveLocationModel.h"

#define MAX_STUDY_ENTRIES 2500
#define MAX_CUTOFF_STEPS 20
#define MAX_LESIONS 100
#define MAX_DETECTIONS 500

typedef struct 
{
	QString patientName, studyDate, studyId, seriesId, localImagePath, paintingFileName, organPath;
	int sliceRange1, sliceRange2;
} StudyEntryStructure;

#define MaxALMModels 12

// active location model
class NIH_ALM_Single_Model
{
public:
	LocationFramework lfw;
	CIS_Vector_JY_double locationVec;
	doubleDynArray pdmPara;

	CIS_3D_Model_Curve *smoothedCord;
	CIS_3D_Model_Cylinder *vertebra_3d_model[6];	// only show the model from t10 to L3
	CIS_3D_Model_MeshSurface *organSurf[5];		// five organs, 0: left kidney, 1: right kidney, 2: liver, 3: pancreas, 4: spleen
	CIS_3D_Model_Needle *organLocationModel[5];
	Vec3 ribCageSize;
	char modelFileName[200];

	NIH_ALM_Single_Model() {};
};

#define MASK_LKIDNEY 1
#define MASK_RKIDNEY 2
#define MASK_LIVER 3
#define MASK_PANCREAS 4
#define MASK_SPLEEN 5

class NIH_BoneSegmentation_Dlg : public NIH_BoneSegmentation_Dlg_Base
{ 
// Construction
public:
	NIH_BoneSegmentation_Dlg(int _display_id, int _display_id2, QWidget* parent = 0, const char* name = 0, WFlags fl = 0 );
	~NIH_BoneSegmentation_Dlg();

	float Version;
	
	int display_id, display_id2;

	CIS_Array_Image3D_short *img3D, *maskImg3D;
	CIS_Array_Image2D_short *img2D, *img2D_x, *img2D_y;
	char fn_path1[200], base_fn[200];
	char study_root[200], study_name[50];

	int curSlice, startSlice, endSlice;
	int curSlice_x, startSlice_x, endSlice_x;
	int curSlice_y, startSlice_y, endSlice_y;

	CIS_3D_Model_Texture *tex_x, *tex_y, *tex_z;
	CIS_3D_Model_Cube *ct_cube;

	// the view center
	int centerx, centery, centerz;
	
	// view control
	int viewMode;

	bool flagShowOverlay, segmentationChanged, flagShowPainting, flagShowBBox, showWSOverlay;

	short baseValue;

	bool debugMode;
	int processingStatus;
	
	// global variable for batch mode
	char patientRoot[300], paintingRoot[300], loadSegPath[300], saveSegPath[300], 
		globalReportFileName[300], globalFeatureFileName[300], globalDetectionFileName[300];
	StudyEntryStructure studyEntries[MAX_STUDY_ENTRIES];
	int numStudyEntries, curStudy, batchStudyStart, batchStudyEnd;
	bool flagBatchMode;
	int batchIndex;		// index to different batch processing
	QString patientAge, patientSex;

	CIS_2D_Model_Texture *cmapModel1;
	CIS_2D_Model_Bitmap *mapModel1;

	CutoffStructure storeMasks[MAX_CUTOFF_STEPS];

	// editing the segmentation
	CIS_2D_ROI_Polygon *paintModel;

	SpineSegmentationMaskStructure maskStruct;
	SpineSegmentationParameter segPara;
	SpineSegmentationInfo segInfo;

	CIS_2D_Model_Rectangle *bboxModel, *vertebraBodyModel, *spinousProcessModel;
	CIS_2D_Model_Polygon *seedModel;

	CIS_Array_Image2D_RGB_uchar *watershedOverlay;

	// local surface
	CIS_3D_Model_MeshSurface *spineSurf, *cordSurf, *metasisSurf, *ribSurf, *paintSurf;

	// disk plane
	CIS_3D_Model_Plane *diskPlaneModel[30];
	Plane3D frontPlane[30], backPlane[30];

	// variable for spine column reformation
	bool reformationMode;
	CIS_Array_Image2D_short *img2D_sag, *img2D_cor;
	CIS_Array_Image2D_short *img2D_sag_mask, *img2D_cor_mask;
	
	// spinal cord model (local coordinate system)
	vec3DynArray spinalCord3D;
	CIS_3D_Model_Scatter *cordModel;
	CIS_3D_Model_Curve *smoothedCord;

	CIS_3D_Model_Line *vertHeightDir;
	CIS_2D_Model_Links *vertHeightRuler;

	CIS_2D_Model_Curve *projectedCord;
	CIS_2D_Model_Curve *leftColumn, *rightColumn;
	CIS_2D_Model_Links *pedicleModel;

	CIS_2D_Model_Curve *projectedCordSag;
	CIS_2D_Model_Curve *leftColumnSag, *rightColumnSag;
	CIS_2D_Model_Links *pedicleModelSag;

	CIS_2D_Model_Curve *projectedCordCor;
	CIS_2D_Model_Curve *leftColumnCor, *rightColumnCor;
	CIS_2D_Model_Links *pedicleModelCor;


	CIS_3D_Model_Volume *volumeModel;

	intDynArray pedicleValley, backValley;
	vec3DynArray spineNormalX, spineNormalY, spineCenter, spineNormalBack;
	doubleDynArray spineWidthLeft, spineWidthRight, spineWidthUp, spineWidthDown;

	// vertebra template
	VertebraStruct2D *vertebraTemplate;
	CIS_2D_Model_Links *templateModel;
	int vertebraTemplateViewMode;

	// for bmd_roi
	intDynArray roiSlices;


	// painting structure
	PaintingStruct paint[MAX_SEED_NUM];
	int totalPaintNum;

	// lesions and detections
	int numLesions, numDetections, numDetections2D, numDetections3D;
	LesionStructure *lesions;
	DetectionStructure *detections;
	FeatureStructure *detections2D;
	FeatureStructure *detections3D;
	PaintingStruct detectionVoxels[MAX_SEED_NUM];

	// classification
	QString svmModelFileName;
	float svmCutoff, svmCutoffMin, svmCutoffMax, svmCutoffStep;
	SvmCommittee *svmCommittee;
	bool flagApplyClassifier;
	bool flagApply3DClassifier;

	float lesionSizeCutoff;
	int useMethod;

	int zoomMode;

	// ALM models
	int totalALMModels, curALMModel, GTLocationModel, refLocationModel;
	NIH_ALM_Single_Model singleALM[MaxALMModels];
	NIH_ALM_Single_Model meanALM;
	ActiveLocationModel_PDM ALM_PDM;
	bool showOrganMask, flagGT;
	CIS_Array_Image3D_uchar *organMask;
	CIS_Array_Image3D_uchar *organGTMask;
	float total_probSum, total_probAvg;

	void ChangeImageSlice(int newCurSlice);
	void ChangeImageSlice_x(int newCurSlice_x);
	void ChangeImageSlice_y(int newCurSlice_y);

	RGBtriple_uchar ComputeVertebraLabelColor(int x, int y, int z, short mask); 

	void ChangeZoom(int newZoom);

	int OpenDataSet(char *fnDir, bool flagDcm);
	int LoadStudy(int st);
	int LoadPaintingFile(const char *painting_fn);
	int WriteOutFeatures(const char *feature_fn);
	int WriteOutFeatures3D( const char *feature_fn);
	int WriteOutDetections(const char *detection_fn);
	void ClearROI();
	int SegmentSingleVertebra(IntVec3 seed, IntVec3 bb1, IntVec3 bb2);
	int VertebraHeightMeasurement();
	int VertebraHeightMeasurement1();
	int VertebraHeightMeasurementSag();
	int NIH_BoneSegmentation_Dlg::load_project(const char* proj_fn);

	int LoadSingleALMModel(const char *fn, int index);
	int IntanstiateSingleALM(NIH_ALM_Single_Model &sALM, CIS_Vector_JY_double &locationVec);

	void SwitchALMDisplayMode(int modelIndex);
	void SaveSingleLocationModel(const char *fn);
	void MapOrganMask(CIS_Array_Image3D_uchar *orgMask,
		Frame &f, int organIndex, int sampleRate, float probMean, float probStd, float &probSum, float &probAvg);
	
	void AlignTwoLocationModelsUsingSpine(int movingModel, int fixedModel);

	void ALMBatchRunning(QString &projFn, QString& ALMFn, int studyID);
	void ALM_Batch_file(const char *batchFn);
	void Load_ALM_Multi_Model_File(const char *model_fn);

protected:
	virtual void UpdateData(bool toEdit=true);

public slots:
    virtual void pushButton_load_s_clicked();
    virtual void pushButton_load_clicked();
    virtual void pushButton_save_s_clicked();
    virtual void slider_slice_valueChanged( int );
    virtual void slider_slice_sliderMoved( int );
    virtual void slider_contrast_lo_valueChanged( int );
    virtual void slider_contrast_hi_valueChanged( int );
    virtual void slider_opacity_valueChanged( int );
    virtual void checkBox_overlay_stateChanged( int );
	virtual void checkBox_woverlay_stateChanged( int );
    virtual void checkBox_painting_stateChanged( int );
    virtual void pushButton_help_clicked();
    virtual void pushButton_refresh_clicked();
    virtual void pushButton_close_clicked();
    virtual void pushButton_load_project_clicked();
    virtual void pushButton_batch_project_clicked();
    virtual void pushButton_pre_process_clicked();
    virtual void pushButton_surface_clicked();
    virtual void listBox_study_entry_clicked( QListBoxItem * );
    virtual void pushButton_spinal_cord_clicked();
    virtual void pushButton_met_detection_clicked();
    virtual void spinBox_lesion_valueChanged( int );
    virtual void spinBox_detection_valueChanged( int );
    virtual void spinBox_detection_2D_valueChanged( int );
	virtual void spinBox_cutoff_valueChanged( int );
    virtual void pushButton_save_seg_clicked();
    virtual void pushButton_save_feature_clicked();
    virtual void checkBox_bounding_box_stateChanged( int );
    virtual void pushButton_classification_clicked();
    virtual void pushButton_process_clicked();
    virtual void pushButton_watershed_clicked();
    virtual void pushButton_cord_clicked();
    virtual void pushButton_BMD_ROI_clicked();
    virtual void pushButton_BMD_ROI_3D_clicked();
    virtual void pushButton_BMD_ROI_L1L2_clicked();
    virtual void pushButton_BMD_ROI_L1L2_3D_clicked();
    virtual void pushButton_spine_partition_clicked();
    virtual void pushButton_Rib_Detection_clicked();
    virtual void pushButton_Vertebra_Seg_clicked();

	virtual void radioButton_view_toggled(bool );
    virtual void radioButton_zoom_toggled(bool );

	virtual void CommandLine_clicked();

	// test live wire
	NIH_LiveWire livewire;
	bool flagUseLiveWire, flagRecordPoint, flagLiveWireAutoMode;
	intVec2DynArray recordedPoints;
	bool flagImageChanged;

	vec2DynArray groundTruth;
	CIS_Array_Image2D_short *groundTruthMask;
	CIS_2D_Model_Polygon *liveWireModel, *liveWireInitModel;
    virtual void useLiveWire_stateChanged( int );
    virtual void checkBox_record_point_stateChanged( int );
    virtual void checkBox_ls_show_init_model_stateChanged( int );
    virtual void pushButton_read_point_clicked();
    virtual void pushButton_write_point_clicked();
    virtual void pushButton_smooth_contour_clicked();
    virtual void testImg_clicked();
 //   virtual void testImg_clicked1();
	void StartLiveWire(IntVec2 pt);
	void AddLiveWireSeed(IntVec2 pt);
	void FreezeLiveWireSeed(IntVec2 pt);

	void ComputeError(vec2DynArray &testContour, vec2DynArray &refContour, CIS_Array_Image2D_short *refMask,
		float &avgError, float &maxError, float &stdError, 
		float &sensitivity, float &specificity, float &dice,
		float &hausdorff, float &smoothness);

    virtual void pushButton_lw_batch_clicked();

	QString lw_reportFileName, lw_recordedFileName, lw_summaryFileName;
	float lw_elapsed_time, lw_avgError, lw_maxError, lw_stdError;
	float lw_sensitivity, lw_specificity, lw_dice, lw_hausdorff, lw_smoothness;
	void ReportResults();

	// test iv surface
    virtual void pushButton_surface_iv_clicked();

	// Active loction model
    virtual void pushButton_Build_ALM_clicked();
    virtual void pushButton_Load_Organ_clicked();
    virtual void pushButton_Save_ALM_Model_clicked();
    virtual void pushButton_Load_ALM_Model_clicked();
    virtual void pushButton_Load_ALM_Multi_Model_clicked();

	virtual void pushButton_Save_ALM_clicked();
	virtual void pushButton_Load_ALM_clicked();
	virtual void pushButton_Mean_ALM_clicked();
	virtual void pushButton_Intanstiate_ALM_clicked();
    virtual void slider_ALM_Mode_valueChanged( int );
    virtual void spinBox_ALM_Mode_valueChanged( int );
    virtual void spinBox_ALM_Model_valueChanged( int );
    virtual void checkBox_flag_GT_stateChanged( int );

	virtual void pushButton_ALM_Refinement_clicked();
	virtual void pushButton_ALM_Batch_clicked();

};

#endif
