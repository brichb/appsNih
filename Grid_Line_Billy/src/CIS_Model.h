
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
// Header file for model class
// 
//----------------------------------------------------------------

#ifndef _CIS_MODEL_H_
#define _CIS_MODEL_H_

#include <stdio.h>
#include <string.h>
#include "vecs.hpp"


#ifndef CIS_OK
#define CIS_OK 0
#endif

#ifndef CIS_ERROR
#define CIS_ERROR 1
#endif

#ifndef DEG2RAD
#define DEG2RAD(x)  x*3.14159265358/180
#endif

#ifndef RAD2DEG
#define RAD2DEG(x)  x*180/3.14159265358
#endif

#ifndef _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_DEPRECATE
#endif

#define MAX_MODELS_IN_GROUP	200

#define HIDE_MODE		0x00
#define LINE_MODE		0x01
#define SURFACE_MODE    0x02
#define VERTEX_MODE     0x03

#define ENTIRE_MODE     0x00
#define CONTOUR_MODE    0x10
#define SLICE_MODE      0x20
#define LAYER_MODE      0x30
#define SECTION_MODE    0x40
#define ACTIVE_LIST_MODE  0x50
#define SHELL_MODE      0x60
#define PORTION1_MODE    0x70
#define PORTION2_MODE    0x80

#define MODEL_3D_MODEL	0xf100
#define MODEL_2D_MODEL	0x1f00

#define MODEL_GROUP		0xf101
#define MODEL_GROUP_2D	0x1f01

#define MODEL_SPHERE    0xf102
#define MODEL_LINE_3D	0xf103
#define MODEL_CUBE		0xf104
#define MODEL_PATCHSURFACE	0xf105
#define MODEL_CORKSCREW		0xf106
#define MODEL_SCATTER		0xf107
#define MODEL_TEXTURE_3D	0xf108
#define MODEL_MESHSURFACE	0xf109
#define MODEL_CURVE_3D		0xf10a
#define MODEL_MESH_TETRAHEDRA	0xf10b
#define MODEL_CYLINDER		0xf10c
#define MODEL_COORDINATE_AXIS_3D	0xf10d
#define MODEL_PLANE_3D	0xf10e

#define MODEL_LINE_2D	0x1f02
#define MODEL_ELLIPSE_2D	0x1f03
#define MODEL_RECTANGLE_2D	0x1f04
#define MODEL_POLYGON_2D	0x1f05
#define MODEL_CONTOUR_2D	0x1f06
#define MODEL_CURVE_2D	0x1f07
#define MODEL_MEDIAL_AXIS_2D	0x1f08
#define MODEL_HISTOGRAM_1D	0x1f10
#define MODEL_HISTOGRAM_2D	0x1f11
#define MODEL_NEEDLE_GRAPH	0x1f12
#define MODEL_LINKS_2D	0x1f13
#define MODEL_TEXTURE_2D	0x1f14
#define MODEL_BITMAP_2D	0x1f15
#define MODEL_TEXT_2D	0x1f16
#define MODEL_COLORBAR_2D	0x1f17
#define MODEL_CROSS_2D	0x1f18
#define MODEL_ARROW_2D	0x1f19

#define PIXEL_TYPE_UCHAR	0x0001
#define PIXEL_TYPE_SHORT	0x0002
#define PIXEL_TYPE_FLOAT	0x0003
#define PIXEL_TYPE_RGB		0x0004

#define MAPPING_MODE_GL     0x0000
#define MAPPING_MODE_IMAGE     0x0001
#define MAPPING_MODE_180     0x0002
#define MAPPING_MODE_FLIPX     0x0003


// The base class for all the models (defines common fields & the basic 
// functions used to process each model.
class CIS_Model
{
public:
	unsigned int shape;     // The model's type
	int reference;			// the number of display that shows this model 
	char name[20];

	Vec3 curr_color;
	double curr_alpha;

	Vec3 orig_color;
	double orig_alpha;

	int draw_mode;
	float lineWidth;

	Frame local_global;
	double gl_matrix[16];

	CIS_Model(){};               // Constructor
	virtual ~CIS_Model(){};

	void IncreaseReference() {reference++;};
	void DecreaseReference() {reference--;};
	int GetReference() const {return reference;};

	char* GetName() { return name;  };
	void SetName(char* newname) { strncpy(name, newname, 20); };
	unsigned int GetShape() const { return shape;  };
	unsigned int GetModelCategory() const { return (shape & 0xff00);};

	virtual void SetFrame(Frame new_f)
	{
		local_global = new_f;

//		global_local = local_global.Inverse();
		gl_matrix[0] = local_global.R.Rx.x;
		gl_matrix[1] = local_global.R.Rx.y;
		gl_matrix[2] = local_global.R.Rx.z;
		gl_matrix[3] = 0;
		gl_matrix[4] = local_global.R.Ry.x;
		gl_matrix[5] = local_global.R.Ry.y;
		gl_matrix[6] = local_global.R.Ry.z;
		gl_matrix[7] = 0;
		gl_matrix[8] = local_global.R.Rz.x;
		gl_matrix[9] = local_global.R.Rz.y;
		gl_matrix[10] = local_global.R.Rz.z;
		gl_matrix[11] = 0;
		gl_matrix[12] = local_global.P.x;
		gl_matrix[13] = local_global.P.y;
		gl_matrix[14] = local_global.P.z;
		gl_matrix[15] = 1;
	};
	virtual Frame GetFrame() const {return local_global;};
	virtual void ComputeFrame() {};
	virtual void GlobalUpdateFrame(Frame delta_f) {};
	virtual void UpdateFrame(Frame delta_f) {};
	virtual void ComputeSize() {};

	virtual int GetDrawMode() const {return draw_mode;};
	virtual void SetDrawMode(int new_mode) {draw_mode=new_mode;};

	virtual void SetOrigColor(Vec3 new_color)
	{ 
		orig_color = new_color;
		curr_color = new_color;
	};
	virtual void SetCurrColor(Vec3 new_color) {curr_color = new_color;};
	virtual Vec3 GetCurrColor() const {return curr_color;};
	virtual void RestoreColor() {curr_color = orig_color;};

	virtual void SetOrigAlpha(double new_alpha)
	{ orig_alpha = new_alpha;
	  curr_alpha = new_alpha;
	};
	virtual void SetCurrAlpha(double new_alpha) {curr_alpha = new_alpha;};
	virtual double GetCurrAlpha() const {return curr_alpha;};
	virtual void RestoreAlpha() {curr_alpha = orig_alpha;};

	virtual float GetLineWidth() const {return lineWidth;};
	virtual void SetLineWidth(float new_line_width) {lineWidth=new_line_width;};
	
	virtual void PrintToFile(FILE *fp) {return;};
	virtual void ReadFromFile(FILE *fp) {return;};
	virtual void OpenGLDraw() { return; };
	virtual int BuildGroup() { return CIS_OK;};

};


#endif