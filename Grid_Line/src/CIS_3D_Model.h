
//----------------------------------------------------------------
// Jianhua Yao
// Computer Integrated Surgery Lab
// ERC for CISST
// Computer Science Department
// The Johns Hopkins University
//
// Version 0.1 - 8/99
//
// Developmental Software - not certified for medical use
//
// (c) copyright 1998 CISST ERC
// This file contains proprietary and confidential information 
// of the Compter Integrated Surgical Systems and Technologies 
// Engineering Research Center. The contents of this file 
// may not be disclosed to third parties, translated, copied, 
// or duplicated in any form without the express written permission 
// of the director, CISST ERC.
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
// Header file for 3D model class
// 
//----------------------------------------------------------------

#ifndef _CIS_3D_MODEL_H_
#define _CIS_3D_MODEL_H_

#include <CIS_Model.h>

// The base class for all the objects (defines common fields & the basic 
// functions used to process each object.
class CIS_3D_Model : public CIS_Model
{
public:
	Vec3 loc;               // Object's location
	Vec3 axis;			  // the center axis of the object
	Vec3 dirX, dirY, dirZ;
	double sizeX, sizeY, sizeZ;
	Vec3 centroid;

	Vec3 scale;

	// material for drawing surface
//	float mat_ambient[4], mat_diffuse[4], mat_specular[4], mat_emission[4];
//	float mat_shininess;

	CIS_3D_Model();               // Constructor
	virtual ~CIS_3D_Model();

//	virtual void SetMaterial();

	virtual void ComputeFrame();
	virtual void GlobalUpdateFrame(Frame delta_f);
	virtual void UpdateFrame(Frame delta_f);
	virtual void SetFrame(Frame new_f);

	virtual void ComputeSize() {};

	virtual void OpenGLDraw() { return; };

	virtual void GenerateProjectiveContour(char *contour_fn) {return;};
	virtual void PrintToFile(FILE *fp) {return;};
	virtual void ReadFromFile(FILE *fp) {return;};
	virtual int BuildGroup() { return CIS_OK;};
	virtual int AddModel(CIS_3D_Model *newModel) { return CIS_OK;};
	virtual int RemoveModel(CIS_3D_Model *_model) { return CIS_OK;};
};

#endif