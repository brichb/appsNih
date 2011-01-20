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
#include <itkImage.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <string.h>
#include "NIH_ITK_Utility.h"
#include "NIH_Algo_Region_Growing.h"
#include "itkIndex.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkMaximumImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkAbsImageFilter.h" 
#include "itkAbsoluteValueDifferenceImageFilter.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkFastMarchingImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkPointSet.h"
#include "itkWarpImageFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkConfidenceConnectedImageFilter.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegistrationMethod.h"
#include "itkTranslationTransform.h"
#include "itkAffineTransform.h"
#include "itkVersorRigid3DTransform.h"
#include "itkDemonsRegistrationFilter.h"
#include "itkAmoebaOptimizer.h"
#include "itkEuler3DTransform.h"
#include "itkEuclideanDistancePointMetric.h"
#include "itkLevenbergMarquardtOptimizer.h"
#include "itkPointSet.h"
#include "itkPointSetToPointSetRegistrationMethod.h"
#include "itkPointSetToImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"
#include "itkGeodesicActiveContourLevelSetImageFilter.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include <itkMutualInformationImageToImageMetric.h> 
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkMutualInformationImageToImageMetric.h"
#include "itkMatchCardinalityImageToImageMetric.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include <itkNearestNeighborInterpolateImageFunction.h> 
#include "itkBSplineInterpolateImageFunction.h"
#include "itkVector.h"
#include "vnl/vnl_math.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkVersorTransformOptimizer.h"
#include "itkVersorRigid3DTransformOptimizer.h"
#include "itkResampleImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include <Dicom_Reader.h>
#include <CIS_Image_Processing_Algo_3D.h>
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImage.h"
#include "itkFixedArray.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkScalarToArrayCastImageFilter.h"
#include <CIS_Image3D.h>
#include "itkRescaleIntensityImageFilter.h"
#include "itkMRFImageFilter.h"
#include "itkDistanceToCentroidMembershipFunction.h"
#include "itkMinimumDecisionRule.h"
#include "itkImageClassifierBase.h"
#include "itkBinaryMask3DMeshSource.h"
#include "itkDeformableMesh3DFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkCovariantVector.h"
#include "itkPointSetToImageFilter.h"
#include "itkMesh.h"
#include "itkLaplacianSegmentationLevelSetImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkFastMarchingImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkZeroCrossingImageFilter.h"


/**
*Performes Active Contour
*
*@author Javed Aman, Jesse Sandberg
* NIH/CC/DRD
*/

/*
Proper usage
ActiveContour.exe <image1.hdr> <image2.hdr>
*/

CIS_Array_Image3D_short *img3D1;
CIS_Array_Image3D_short *img3Dout;

//double ComputeVolume_frust(int top_slice,int bottom_slice,double xypixel_conversion);
//double Height();
const int NumSlices = 70;
const int BufferZone = 15;
/**
* Main function
* 
*/
void main(int argc, char **argv){

	if(argc != 3){
		std::cout<<"Usage: _appsNih_Clinical_Results.exe <image1.hdr> <Starting Slice> "<<std::endl;
		exit(1);
	}
	
	img3D1 = new CIS_Array_Image3D_short();
	img3D1->Load(argv[1]);
//	img3D1->Set_Pixel_SizeZ(5);
	img3Dout = new CIS_Array_Image3D_short(img3D1->Get_SizeX(), img3D1->Get_SizeY(), img3D1->Get_SizeZ()); //NumSlices
	
	int first_slice = atoi(argv[2]);
	img3Dout->Set_SizeX(img3D1->Get_SizeX());
	img3Dout->Set_SizeY(img3D1->Get_SizeY());
	img3Dout->Set_SizeZ(NumSlices);
//	img3Dout->Set_SizeZ(img3D1->Get_SizeZ());
	img3Dout->Set_Pixel_SizeX(img3D1->Get_Pixel_SizeX());
	img3Dout->Set_Pixel_SizeY(img3D1->Get_Pixel_SizeY());
	img3Dout->Set_Pixel_SizeZ(img3D1->Get_Pixel_SizeZ());
	//Might need to set center positions

	for (int x = 0; x<img3D1->Num_Cols(); x++)
		for (int y = 0; y<img3D1->Num_Rows(); y++)
			for (int z = 0; z<BufferZone; z++)
			{
				img3Dout->FastSet(x,y,z,0);
			}

	for (int x = 0; x<img3D1->Num_Cols(); x++)
		for (int y = 0; y<img3D1->Num_Rows(); y++)
			for (int z = BufferZone; z<NumSlices; z++)
			{
				int zin = z + first_slice - BufferZone;
				if (zin<img3D1->Num_Levels())
					img3Dout->FastSet(x,y,z,img3D1->FastGet(x,y,zin));
				else
					img3Dout->FastSet(x,y,z,0);
			}
	int name_size = strlen(argv[1]);
	char *OI_name;
	OI_name	=new char[name_size+3];
//	char OI_name[20];
	for(int i = 0; i < name_size - 4; i++)
	{
		OI_name[i] = argv[1][i];
	}

	OI_name[name_size - 4] = 'O';
	OI_name[name_size - 3] = 'I';
	OI_name[name_size - 2] = '.';
	OI_name[name_size - 1] = 'h';
	OI_name[name_size + 0] = 'd';
	OI_name[name_size + 1] = 'r';
	OI_name[name_size + 2] = '\0';
	Write_Analyze_File(OI_name,*img3Dout);
//	char OI_name1 = argv[3];
//	Write_Analyze_File(OI_name1,*img3Dout);
/*
	int sizex=img3Dout->Num_Cols();
	int sizey=img3Dout->Num_Rows();
	int sizez=img3Dout->Num_Levels();
	for (int z=0; z<sizez; z++)
		for (int y=0; y<sizey; y++)
			for (int x=0; x<sizex; x++)
				if (x<sizex/2)
					img3Dout->FastSet(x,y,z,0);

	char *OCI_name;
	OCI_name=new char[name_size+4];
//	char OCI_name[20];
	for(int i = 0; i < name_size - 4; i++)
	{
		OCI_name[i] = argv[1][i];
	}

	OCI_name[name_size - 4] = 'O';
	OCI_name[name_size - 3] = 'C';
	OCI_name[name_size - 2] = 'I';
	OCI_name[name_size - 1] = '.';
	OCI_name[name_size + 0] = 'h';
	OCI_name[name_size + 1] = 'd';
	OCI_name[name_size + 2] = 'r';
	OCI_name[name_size + 3] = '\0';

	Write_Analyze_File(OCI_name,*img3Dout);
*/
	typedef itk::Image<float, 3> ImageType_3D;
	ImageType_3D::Pointer ITKimg=ImageType_3D::New();
	typedef itk::ThresholdImageFilter<ImageType_3D>FilterType;
	Copy_CISImage_to_ITKImage(img3Dout, ITKimg);
	FilterType::Pointer filter=FilterType::New();
	filter->SetInput(ITKimg);
	int Low_Thresh = 1050;
	int High_Thresh = 1200;

	filter->SetOutsideValue(0);
	filter->ThresholdOutside(Low_Thresh,High_Thresh);
	filter->Update();
	Copy_ITKImage_to_CISImage(filter->GetOutput(), img3Dout);
/*
	int name_size2 = strlen(argv[1]);
	char *TOI_name;
	TOI_name=new char[name_size2+4];
//	char TCI_name[20];
	for(int i = 0; i < name_size2 - 4; i++)
	{
		TOI_name[i] = argv[1][i];
	}
*/
	int name_size2 = strlen(argv[1]);
	char *TI_name;
	TI_name	=new char[name_size2+3];
//	char OI_name[20];
	for(int i = 0; i < name_size2 - 4; i++)
	{
		TI_name[i] = argv[1][i];
	}
	TI_name[name_size2 - 4] = 'T';
	TI_name[name_size2 - 3] = 'I';
	TI_name[name_size2 - 2] = '.';
	TI_name[name_size2 - 1] = 'h';
	TI_name[name_size2 + 0] = 'd';
	TI_name[name_size2 + 1] = 'r';
	TI_name[name_size2 + 2] = '\0';

	Write_Analyze_File(TI_name,*img3Dout);

}