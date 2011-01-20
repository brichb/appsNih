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
*Performes Pre Image Processing before Active Contour
*
Jesse Sandberg
* NIH/CC/DRD
*/

/*
Proper usage
Threshold.exe <image1.hdr>
*/

CIS_Array_Image3D_short *img3D1;

/**
* Main function
* 
*/
void main(int argc, char **argv){

	if(argc < 5){
		std::cout<<"Usage: _Threshold.exe <File.hdr> OutPutFile.hdr ThresholdValueMax ThreshMin"<<std::endl;
		exit(1);
	}
	
	img3D1 = new CIS_Array_Image3D_short();
	img3D1->Load(argv[1]);

	int sizex=img3D1->Num_Cols();
	int sizey=img3D1->Num_Rows();
	int sizez=img3D1->Num_Levels();
int max = atof(argv[3]);
int min = atof(argv[4]);

	for (int z=0; z<sizez; z++)
		for (int y=0; y<sizey; y++)
			for (int x=0; x<sizex; x++)
			{
				if (img3D1->FastGet(x,y,z) < max && img3D1->FastGet(x,y,z) > min)
					img3D1->FastSet(x,y,z,1000);
				else img3D1->FastSet(x,y,z,0);
			}
/*
	int name_size1 = strlen(argv[1]);
	char *liver_name;
	liver_name=new char[name_size1+4];
	for(int i = 0; i < name_size1 - 4; i++)
	{
		liver_name[i] = argv[1][i];
	}

	liver_name[name_size1 - 4] = 'T';
	liver_name[name_size1 - 3] = 'H';
	liver_name[name_size1 - 2] = 'R';
	liver_name[name_size1 - 1] = '.';
	liver_name[name_size1 + 0] = 'h';
	liver_name[name_size1 + 1] = 'd';
	liver_name[name_size1 + 2] = 'r';
	liver_name[name_size1 + 3] = '\0';
*/
	Write_Analyze_File(argv[2],*img3D1);

}