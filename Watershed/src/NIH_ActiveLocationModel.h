#ifndef _NIH_ActiveLocationModel_H_
#define _NIH_ActiveLocationModel_H_

#include <vecs.hpp>
#include <CIS_Matrix_JY.h>

class ObjectLocationModel
{
public: 
	Vec3 center;
	Rotation orientation;
	Vec3 size;
	char name[20];
	int shape;
	int state;		// 0: not available

	float meanIntensity, varianceIntensity, probSum, probAvg;
	float volumn;

	ObjectLocationModel();
	~ObjectLocationModel();
};

class RelativeLocationModel
{
public: 
	Vec3 relativeLocation;
	Rotation relativeOrientation;
	char name[30];
	int state;		// 0: not available

	ObjectLocationModel *src, *target;

	RelativeLocationModel();
	~RelativeLocationModel();
};


class VertebraLocationModel : public ObjectLocationModel
{
public: 
	VertebraLocationModel();
	~VertebraLocationModel();
};

class OrganLocationModel : public ObjectLocationModel
{
public: 
//	RelativeLocationModel *toL1, *toL2, *toT11, *toT12;
	RelativeLocationModel toTL[6];
	char surfaceFilename[200], binaryFilename[200];

	OrganLocationModel();
	~OrganLocationModel();
};

class LocationFramework
{
public:
	VertebraLocationModel *tL[17];		// 17 vertebrae, t1 to t12, l1 to l5

	OrganLocationModel *organ[5];		// five organs, 0: left kidney, 1: right kidney, 2: liver, 3: pancreas, 4: spleen

	LocationFramework();
	~LocationFramework();
};

#define MAX_ALM_TRAINING_MODELS 100

class ActiveLocationModel_PDM
{
public:
	int numTrainingVec;
	int vecLength;
	char referenceModelFileName[200];

	Vec3 referenceModelSize;

	intDynArray vertebraIndexList, organIndexList;		// list to store the index of vertebra and organs in the PDM model

	CIS_Vector_JY_double trainingVec[MAX_ALM_TRAINING_MODELS];

	CIS_Vector_JY_double meanVec, orgVec;

	int pdmDegree;
	CIS_Matrix_JY_double pdmModel;
	doubleDynArray pdmRange, pdmValue;

	ActiveLocationModel_PDM();
	~ActiveLocationModel_PDM();

	int AddTrainingVec(CIS_Vector_JY_double &tVec);
	int TrainPDM();
	CIS_Vector_JY_double InstantiatePDM(doubleDynArray &b);
	CIS_Vector_JY_double InstantiatePDM();

	int InversePDM(CIS_Vector_JY_double &tVec, doubleDynArray &b);

	int SaveALMModel(const char *fn);
	int LoadALMModel(const char *fn);
};

#endif
