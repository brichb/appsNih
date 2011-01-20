#include "NIH_SpineSegmentation_DataStructure.h"

SpineSegmentationMaskStructure::SpineSegmentationMaskStructure()
{
	mask_spinalCord = 21;
	mask_boneMetasis = 22;
	mask_vertebralDisk = 23;

	mask_corticalBone = 10;
	mask_spongyBone = 11;
	mask_otherBone = 12;
	mask_falseDetection = 13;
	mask_rib = 14;
	
	mask_vertebra = 40;
	
	mask_body = 2;
	mask_air = 1;
	mask_paint = 31;
	mask_paintLytic = 32;
}

SpineSegmentationMaskStructure::~SpineSegmentationMaskStructure()
{
}

void SpineSegmentationMaskStructure::ReadFromFile(const char *fileName)
{
}


SpineSegmentationParameter::SpineSegmentationParameter()
{
	// prefined bounding box for spines
	predefinedBound1.x = 128;
	predefinedBound1.y = 256;
	predefinedBound2.x = 384;
	predefinedBound2.y = 512;

	boneThresh = 1200;			// was 1200 for bone mat
	airThresh = 300;
	closeCorticalHoleIteration = 1;
	cordIntensityThresh = 1200; // was 1100 for bone mat
	cordSizeThresh = 100;		// was 30 for bone mat
	closeCordIteration = 3;			// was 2 for bone mat
	metasisSizeThresh = 5;
	metasisIntensityThresh = 1150;
	requiredMetContrast = 10;

	maxThresh = 3000; //to eliminate metal objects

	softThresh = 1000;		// lower threshold for soft tissues
	borderMinGrad = 50;		// the minimum gradian of border for vertebra
}

SpineSegmentationParameter::~SpineSegmentationParameter()
{
}

void SpineSegmentationParameter::ReadFromFile(const char *fileName)
{
}

SpineSegmentationInfo::SpineSegmentationInfo()
{
	bound1 = bound2 = IntVec3(0,0,0);
	sacrum_start = -1;
	t12_vertebra = -1;

}

SpineSegmentationInfo::~SpineSegmentationInfo()
{
}

void SpineSegmentationInfo::Initialize()
{
	bound1 = bound2 = IntVec3(0,0,0);
	avg_cortical_intensity = avg_spongy_intensity = avg_bone_intensity = 0;
	spineBound1.SetSize(0); spineBound2.SetSize(0);
	diskBound1.SetSize(0); diskBound2.SetSize(0);
	cordBound1.SetSize(0); cordBound2.SetSize(0);
	sprocessBound1.SetSize(0); sprocessBound2.SetSize(0);
}

void SpineSegmentationInfo::WriteToFile(const char *fileName)
{
}


VertebraStruct2D::VertebraStruct2D()
{
	diskCenter = Vec2(-1,-1);
	cordCenter = Vec2(-1,-1);
	diskRadius = -1;
	cordRadius = Vec2(-1,-1);

	int a;
	for(a=0; a<diskAngleCount; a++)
	{
		diskContour[a] = Vec2(-1, -1);
		diskContourRadius[a] = -1;
		interpolatedRadius[a] = -1;
	}
	diskNeckLeft = Vec2(-1, -1);
	diskNeckRight = Vec2(-1, -1);
	angleLeft = angleRight = -1;

	sprocessEnd = Vec2(-1,-1);
	int i;
	for(i=0; i<sprocessCount; i++)
	{
		sprocessLeft[i] = sprocessRight[i] = sprocessMedial[i] = Vec2(-1, -1);
		sprocessLeftWidth[i] = sprocessRightWidth[i] = -1;
	}
	sprocessLength = sprocessOrient = -1;

	leftPedicleEnd = rightPedicleEnd = Vec2(-1, -1);
	leftPedicleLength = rightPedicleLength = -1;
	for(i=0; i<pedicleCount; i++)
	{
		leftPedicleUp[i] = leftPedicleDown[i] = rightPedicleUp[i] = rightPedicleDown[i] = Vec2(-1,-1);
		leftPedicleMedial[i] = rightPedicleMedial[i] = Vec2(-1, -1);
	}

}

VertebraStruct2D::~VertebraStruct2D()
{
}