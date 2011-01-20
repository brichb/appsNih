
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
// Header file for 2D model class
// 
//----------------------------------------------------------------

#ifndef _CIS_2D_MODEL_H_
#define _CIS_2D_MODEL_H_

#include <CIS_Model.h>

// The base class for all the 2D models (defines common fields & the basic 
// functions used to process each 2D model.
class CIS_2D_Model : public CIS_Model
{
public:
	Vec2 loc;               // Model's location
	Vec2 size;
	Vec2 centroid;
	int mapping_mode;
	double zoom_x, zoom_y;

	Vec2 scale;
	double xyRatio;
	bool keepRatio;

	CIS_2D_Model();               // Constructor
	virtual ~CIS_2D_Model();

	void SetLocation(Vec2 newLoc) { loc = newLoc;};
	Vec2 GetLocation() const {return loc;};
	void SetSize(Vec2 newSize);
	Vec2 GetSize() const {return size;};
	void SetCentroid(Vec2 newCen) { centroid = newCen;};
	Vec2 GetCentroid() const {return centroid;};
	void SetScale(Vec2 newScale) { scale = newScale;};
	Vec2 GetScale() const {return scale;};

	void SetXYRatio(double newRatio) { keepRatio=true; xyRatio = newRatio;};
	double GetXYRatio() const {return xyRatio;};
	void KeepRatio(bool enable=true);

	void SetZoom(double _zoomx, double _zoomy) {zoom_x=_zoomx; zoom_y=_zoomy;};
	double GetZoomX() const {return zoom_x;};
	double GetZoomY() const {return zoom_y;};

	void SetMappingMode(int new_mode) {mapping_mode=new_mode;};
	int GetMappingMode() const {return mapping_mode;};

	virtual void Translate2D(Vec2 trans);
	virtual void Rotate2D(double angle, Vec2 rotCenter=Vec2(0,0));

	virtual void OpenGLDraw() { return; };
	virtual void ComputeSize() {return;};

	virtual void PrintToFile(FILE *fp) {return;};
	virtual void ReadFromFile(FILE *fp) {return;};
};



#endif
