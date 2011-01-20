#ifndef _NIH_SPINE_SEGMENTATION_DATASTRUCTURE_H_
#define _NIH_SPINE_SEGMENTATION_DATASTRUCTURE_H_

#include <VecsDynArray.h>
#include <vector>

typedef struct 
{
	int lesionKey;		// the detection key
	short lesionType;	// lesion type: lytic or scloratic
	float lesionSize;	// manually measured diameter
	int lesionArea;		// total number of pixels
	int lesionSlices; // occupied slices
	int lesionStatus; // detected or not
	IntVec3 lesionLoc;	// centroid location
	bool hit;
} LesionStructure;

typedef struct 
{
	short detectionType;
	float detectionVolume;
	float detectionSize;
	int matchedLesion; 
	int minz;
	int maxz;
	IntVec3 detectionLoc;
} DetectionStructure;

#define diskAngleInterval 10
#define diskAngleCount 36
#define sprocessCount 10
#define pedicleCount 10
#define MAX_LESIONS 100
#define MAX_DETECTIONS 500

class VertebraStruct2D
{
public:
	VertebraStruct2D();
	~VertebraStruct2D();

	// overall template (initial)
	Vec2 diskCenter, cordCenter;
	double diskRadius;
	Vec2 cordRadius;

	// cord structure
	short isKey, cordSize;
	Vec2 cordContour[diskAngleCount];
	double cordContourRadius[diskAngleCount];
	IntVec2 cord_bb0, cord_bb1;

	// disk structure
	Vec2 diskContour[diskAngleCount];
	double diskContourRadius[diskAngleCount], interpolatedRadius[diskAngleCount];
	Vec2 diskNeckLeft, diskNeckRight;
	int angleLeft, angleRight;


	// spinal process structure
	Vec2 sprocessEnd;
	Vec2 sprocessLeft[sprocessCount], sprocessRight[sprocessCount], sprocessMedial[sprocessCount];
	double sprocessLeftWidth[sprocessCount], sprocessRightWidth[sprocessCount];
	double sprocessLength, sprocessOrient;

	// pedicle structure
	Vec2 leftPedicleEnd, rightPedicleEnd;
	double leftPedicleLength, rightPedicleLength;
	Vec2 leftPedicleUp[pedicleCount], leftPedicleDown[pedicleCount], leftPedicleMedial[pedicleCount];
	Vec2 rightPedicleUp[pedicleCount], rightPedicleDown[pedicleCount], rightPedicleMedial[pedicleCount];
} ;


typedef struct 
{
	// location
	int centerx, centery, centerz;
	
	float relCoordx, relCoordy, relCoordz;
	float distToBoundary;

	// shape feature
	float area, volume, perimeter, surfaceArea;
	float primaryAxisLength, secondaryAxisLength;
	float outerBorderRatio, aspectRatio;
	float spherecity, compactness, roundness;

	// some border features
	float corticalBorderRatio, bboxBorderRatio, airBorderRatio, cordBorderRatio, ribBorderRatio;


	// radial distance measures
	float shapeComplexity_f1, shapeComplexity_f2, shapeComplexity_f21;

	// intensity features
	float meanIntensity, stdevIntensity, skewnessIntensity, kurtosisIntensity;
	float interiorIntensity, borderIntensity;
	float outsideIntensity, outsideIntensityDev;
	float innerOuterContrast;

	// neighborhood info
	float borderThickness;
	float neighborArea, neighborIntensity;



	// status
	int matchedLesion, matchedLesionType; 
	float matchedLesionSize, matchedOverlap;

	// classification
	float svmVote;
	float svmScore;
} FeatureStructure;

typedef struct
{
	// mask
	//short *mask3DArray;
	std::vector<short> mask3DArray;

	// lesion statuses
	int lesionStatus[MAX_LESIONS];

	// detections
	int matchedLesion[MAX_DETECTIONS];

	// detection stats
	int totalMet, foundMet, totalDetection, tpDetection, fpDetection;
	float sensitivity;
	float curCutoff;
} CutoffStructure;

class SpineSegmentationMaskStructure
{
public:
	SpineSegmentationMaskStructure();
	~SpineSegmentationMaskStructure();

	short mask_corticalBone, mask_spongyBone, mask_rib;
	short mask_air, mask_body;
	short mask_otherBone, mask_spinalCord, mask_boneMetasis, mask_paint, mask_paintLytic, mask_vertebralDisk;
	short mask_vertebra;
	short mask_falseDetection;
	void ReadFromFile(const char *fileName);
};


class SpineSegmentationParameter
{
public:
	SpineSegmentationParameter();
	~SpineSegmentationParameter();

	IntVec2 predefinedBound1, predefinedBound2;
		
	short boneThresh, airThresh, softThresh, maxThresh; //maxThresh to eliminate metal objects
	int closeCorticalHoleIteration;
	short cordIntensityThresh;
	int cordSizeThresh;
	int closeCordIteration;
	int metasisSizeThresh;
	int metasisIntensityThresh;
	int requiredMetContrast;

	int borderMinGrad;

	void ReadFromFile(const char *fileName);
};

class SpineSegmentationInfo
{
public:
	SpineSegmentationInfo();
	~SpineSegmentationInfo();

	// bounding box
	IntVec3 bound1, bound2;	// bounding box for all bones
	intVec2DynArray spineBound1, spineBound2;	// bounding box for spine on each slice
	intVec2DynArray diskBound1, diskBound2;		// bounding box for disk on each slice
	intVec2DynArray cordBound1, cordBound2;		// bounding box for cord on each slice
	intVec2DynArray sprocessBound1, sprocessBound2;		// bounding box for spinous process on each slice

	intVec2DynArray cordCenter, cordRadius;
	intVec2DynArray diskCenter, diskRadius;

	int sacrum_start;
	int t12_vertebra;
	Vec3 ribCageSize;
	

	// overall intensity
	float avg_cortical_intensity, avg_spongy_intensity, avg_bone_intensity, avg_spinal_cord_intensity;

	void Initialize();
	void WriteToFile(const char *fileName);
};

#endif