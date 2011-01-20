#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <vnl/vnl_vector.h>
#include "DynArray.h"
#include "NIH_ITK_Utility.h"
//#include <NIH_OpenGL_QT_interface.h>
#include <CIS_Model_Package.h>
#include <Dicom_Reader.h>
#include <CIS_Image_Processing_Algo_3D.h>
#include <CIS_Array_Image3D.h>
#include "itkImageRegistrationMethod.h"
#include "itkVersorTransform.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkImage.h"
#include "itkResampleImageFilter.h"
#include "CIS_Algo_Contour.h"

  const    unsigned int    Dimension = 3;
  typedef  float           PixelType;
  typedef itk::Image< PixelType, Dimension >  ImageType3D;
  typedef itk::VersorTransform< double > TransformType;
  typedef itk:: NearestNeighborInterpolateImageFunction< ImageType3D, double>    InterpolatorType;
  typedef itk::ImageRegistrationMethod< ImageType3D, ImageType3D >    RegistrationType;
  typedef RegistrationType::ParametersType ParametersType;
  typedef itk::ResampleImageFilter< ImageType3D, ImageType3D >    ResampleFilterType;
  typedef TransformType::InputPointType CenterType;


class volumeforPCA{

public: 
	vec3DynArray * points;
	Vec3 center, centroid3D, primaryAxis3D, secondaryAxis3D, thirdAxis3D, axisLength3D;
	volumeforPCA(){
		points = new vec3DynArray();
	}
	void Add(Vec3 pt){
		points->Add(pt);
	}
	int GetSize(){
		return points->GetSize();
	}
	volumeforPCA & operator= (const volumeforPCA & b){
		for(int i = 0; i < b.points->GetSize(); i++)
			points->Add(b.points->GetAt(i));
		center = b.center;
		centroid3D = b.centroid3D;
		primaryAxis3D = b.primaryAxis3D;
		secondaryAxis3D = b.secondaryAxis3D;
		thirdAxis3D = b.thirdAxis3D;
		axisLength3D = b.axisLength3D;
		return *this;
	}
};

std::vector<double> Principle_Components_Finder(CIS_Array_Image3D_short * img3Dconv)
{

	//CIS_Array_Image3D_short *img3Dconv;
	//std::string Run_Name = Name;
	//img3Dconv = new CIS_Array_Image3D_short();

	//img3Dconv->Load("C:\\cvs\\CIPS\\projects\\John\\MOMA_MICCAI_2010\\PatientData\\Gahagan\\InitialGahagan_Liver.hdr");

	// Begin Centroid calculation.
	int counter = 0;
	int xsum = 0;
	int ysum = 0;
	int zsum = 0;

	for (int z = 0; z<img3Dconv->Get_SizeZ(); z++){
		for (int y = 0; y<img3Dconv->Get_SizeY(); y++){
			for (int x = 0; x<img3Dconv->Get_SizeX(); x++){
				if (img3Dconv->FastGet(x,y,z) > 0){
					xsum += x;
					ysum += y;
					zsum += z;
					counter++;
				}
			}
		}
	}

	double xcent= (double)(xsum/counter) - floor((double)(xsum/counter));
	double ycent = (double)(ysum/counter) - floor((double)(ysum/counter));
	double zcent = (double)(zsum/counter) - floor((double)(zsum/counter));
	
	int xcenter = 0, ycenter = 0, zcenter = 0;

	if(xcent < 0.5)
	{
		xcenter = floor((double)(xsum/counter));
	}
	else
	{
		xcenter = ceil((double)(xsum/counter));
	}
	if(ycent < 0.5)
	{
		ycenter = floor((double)(ysum/counter));
	}
	else
	{
		ycenter = ceil((double)(ysum/counter));
	}
	if(zcent < 0.5)
	{
		zcenter = floor((double)(zsum/counter));
	}
	else
	{
		zcenter = ceil((double)(zsum/counter));
	}

	int fixedCOM_x = xcenter;
	int fixedCOM_y = ycenter;
	int fixedCOM_z = zcenter;

	// End Centroid calculation.

	// Begin Compute Principle Components

	volumeforPCA * VolumeConv; //Convolution Image
	//volumeforPCA * VolumeAtlas; //Atlas Image

	VolumeConv = new volumeforPCA();
	//VolumeAtlas = new volumeforPCA();


	for (int z = 0; z<img3Dconv->Get_SizeZ(); z++)
		for (int y = 0; y<img3Dconv->Get_SizeY(); y++)
			for (int x = 0; x<img3Dconv->Get_SizeX(); x++)
			{
				Vec3 point;
				if (img3Dconv->FastGet(x,y,z) > 0)
				{
					point.x = x;
					point.y = y;
					point.z = z;
					VolumeConv->Add( point );
				}

			}

	


	CIS_Algo_3DPointCloud_GetSecondMoment(*VolumeConv->points, VolumeConv->centroid3D, VolumeConv->primaryAxis3D, VolumeConv->secondaryAxis3D, VolumeConv->thirdAxis3D, VolumeConv->axisLength3D);

	// JW COMMENT - Here is the first PC axis.
	Vec3 primaryAxis1 = VolumeConv->primaryAxis3D;
	Vec3 primaryAxis_conv = primaryAxis1.normalize();
	std::cout<<primaryAxis_conv.x<<std::endl;
	std::cout<<primaryAxis_conv.y<<std::endl;
	std::cout<<primaryAxis_conv.z<<std::endl;

	std::vector<double> Principle_Components;
	 Principle_Components.push_back(primaryAxis_conv.x);
	 Principle_Components.push_back(primaryAxis_conv.y);
	 Principle_Components.push_back(primaryAxis_conv.z);


	 Vec3 secondaryAxis1 = VolumeConv->secondaryAxis3D;
	Vec3 secondaryAxis_conv = secondaryAxis1.normalize();
	std::cout<<secondaryAxis_conv.x<<std::endl;
	std::cout<<secondaryAxis_conv.y<<std::endl;
	std::cout<<secondaryAxis_conv.z<<std::endl;

	 
	 Principle_Components.push_back(secondaryAxis_conv.x);
	 Principle_Components.push_back(secondaryAxis_conv.y);
	 Principle_Components.push_back(secondaryAxis_conv.z);

	 

	 Vec3 thirdAxis1 = VolumeConv->thirdAxis3D;
	Vec3 thirdAxis_conv = thirdAxis1.normalize();
	std::cout<<thirdAxis_conv.x<<std::endl;
	std::cout<<thirdAxis_conv.y<<std::endl;
	std::cout<<thirdAxis_conv.z<<std::endl;

 
	 Principle_Components.push_back(thirdAxis_conv.x);
	 Principle_Components.push_back(thirdAxis_conv.y);
	 Principle_Components.push_back(thirdAxis_conv.z);




	return Principle_Components;
 
}
