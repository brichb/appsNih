#include <itkImage.h>
#include <CIS_Array_Image3D.h>
#include <CIS_Array_Image2D.h>
#include <CIS_Image_Processing_Algo_2D.h>

#include "nih_itk_Utility.h"

//#include <fstream.h>
//#include <qdatetime.h>
//#include <ArithDynArray.h>
//#include <qimage.h>

#include <itkGradientAnisotropicDiffusionImageFilter.h> 
#include "itkWatershedImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkVectorCastImageFilter.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkScalarToRGBPixelFunctor.h"
#include <itkWatershedMiniPipelineProgressCommand.h> 
#include "itkRescaleIntensityImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"

#include "itkBinaryMask3DMeshSource.h"
#include "itkDeformableMesh3DFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkImage.h"
#include "itkMesh.h"
#include "itkCovariantVector.h"
#include "itkPointSetToImageFilter.h"
#include "itkImageFileReader.h"

//#include "NIH_SpineSegmentation_DataStructure.h"

#include <cis_array_image2d_rgb.h>
#include <math.h>
#include <NIH_Algo_Region_Growing.h>

extern void MyGraphCuts(unsigned short *pInImg, unsigned short* pOutImg, int **neigh, int **hist, int numRegions,
				 double InSmoothDelta, double InSmoothK, double InDataMean1, double InDataDelta1, double InDataMean2, double InDataDelta2,
				 double InDataK_sink, double InDataK_source, int nSourcePointsNum, intDynArray &pSourcePos, int nSinkPointsNum, intDynArray &pSinkPos,
				 int expansion_iteration);

int SpineSegmentation_Watershed_usingITK(CIS_Array_Image2D_short *img, CIS_Array_Image2D_short *smoothedData, 
										 CIS_Array_Image2D_short *waterLabel, bool debugMode)
{
	itkImage_2D_float::Pointer itk_img=itkImage_2D_float::New();
	typedef itk::GradientAnisotropicDiffusionImageFilter<itkImage_2D_float,
		itkImage_2D_float>  DiffusionFilterType;
	typedef itk::WatershedImageFilter<itkImage_2D_float> WatershedFilterType;
  
	typedef itk::RGBPixel<unsigned char>   RGBPixelType;
	typedef itk::Image<RGBPixelType, 2>    RGBImageType;
	typedef itk::Image<unsigned long, 2>   LabeledImageType;
	typedef itk::ImageFileWriter<RGBImageType> FileWriterType;
	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<
               itkImage_2D_float, itkImage_2D_float >  GradientMagnitudeFilterType;


	img->Save("c:\\tmp\\original.tif");
	// convert to ITK image
	Copy_CISImage_to_ITKImage_2D(img, itk_img);

	// read parameters through a file
	int diffusionIteration, useDiff;
	double conductance, timeStep;
	double gaussSigma;
	double waterLevel, waterThreshold;

	// set default values
	diffusionIteration = 3; // 3
	conductance = 3.0;
	timeStep = 0.125;
	gaussSigma = 1;
	waterLevel = 0.15; // 0.15
	waterThreshold = 0.001;// 0.001
	useDiff = 1;

	// read from a file
	FILE *fp;
	if((fp=fopen("c:\\tmp\\water_para.txt", "r"))!=NULL)
	{
		fscanf(fp, "%d %d %lf %lf %1f %lf %lf", &useDiff, &diffusionIteration, &conductance, &timeStep, &gaussSigma,
			&waterLevel, &waterThreshold);

		fclose(fp);
	}

	// anisotropic diffusion
	DiffusionFilterType::Pointer diffusion = DiffusionFilterType::New();
	diffusion->SetNumberOfIterations(diffusionIteration);
	diffusion->SetConductanceParameter( conductance );
	diffusion->SetTimeStep(timeStep);

	GradientMagnitudeFilterType::Pointer gradient = GradientMagnitudeFilterType::New();
	gradient->SetSigma(gaussSigma);

	WatershedFilterType::Pointer watershed = WatershedFilterType::New();
	watershed->SetLevel( waterLevel );
	watershed->SetThreshold( waterThreshold );

	typedef itk::Functor::ScalarToRGBPixelFunctor<unsigned long>
		ColorMapFunctorType;
	typedef itk::UnaryFunctorImageFilter<LabeledImageType,
		RGBImageType, ColorMapFunctorType> ColorMapFilterType;
	ColorMapFilterType::Pointer colormapper = ColorMapFilterType::New();

	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName("c:\\tmp\\watershed0.png");

	
	if(useDiff)
	{
		diffusion->SetInput(itk_img);
		gradient->SetInput(diffusion->GetOutput());
		watershed->SetInput(gradient->GetOutput());
//		watershed->SetInput(diffusion->GetOutput());
	}
	else 
	{
		gradient->SetInput(itk_img);
		watershed->SetInput(gradient->GetOutput());
	}

	if(debugMode)
	{
		colormapper->SetInput(watershed->GetOutput());
		writer->SetInput(colormapper->GetOutput());
	}

	typedef unsigned char WritePixelType;
	typedef itk::Image< WritePixelType, 2 > WriteImageType;
	typedef itk::RescaleIntensityImageFilter< 
               itkImage_2D_float, WriteImageType > RescaleFilterType;

	if(debugMode)
	{
		RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
/*		rescaler->SetOutputMinimum(   0 );
		rescaler->SetOutputMaximum( 255 );
		typedef itk::ImageFileWriter<WriteImageType> FileWriterType_gray;
		FileWriterType_gray::Pointer writer1 = FileWriterType_gray::New();
		writer1->SetFileName("d:\\tmp\\smoother.png");
		rescaler->SetInput( diffusion->GetOutput() );
		writer1->SetInput( rescaler->GetOutput() );
		writer1->Update();
*/
/*		RescaleFilterType::Pointer rescaler2 = RescaleFilterType::New();
		rescaler2->SetOutputMinimum(   0 );
		rescaler2->SetOutputMaximum( 255 );
		FileWriterType_gray::Pointer writer2 = FileWriterType_gray::New();
		writer2->SetFileName("d:\\tmp\\gradient.png");
		rescaler2->SetInput( gradient->GetOutput() );
		writer2->SetInput( rescaler2->GetOutput() );
		writer2->Update();
*/	}

	try 
	{
		watershed->Update();
		colormapper->Update();
		Save_ITKRGBImage(colormapper->GetOutput(), "c:\\tmp\\watershed0.jpg");
///		if(debugMode) writer->Update();
    }
	catch (itk::ExceptionObject &e)
    {
		std::cerr << e << std::endl;
    }
    
	// get some result
	Copy_ITKImage_to_CISImage_2D(diffusion->GetOutput(), smoothedData);

	// get the label
	short *labelArray;
	labelArray = waterLabel->GetArray();

	typedef itk::ImageRegionConstIterator< LabeledImageType > ConstIteratorType_label;
	ConstIteratorType_label::RegionType itk_region;
	itk_region = watershed->GetOutput()->GetLargestPossibleRegion();

	ConstIteratorType_label imgIt(watershed->GetOutput(), itk_region);

	int k;
	for(imgIt.GoToBegin(), k=0; !imgIt.IsAtEnd(); ++imgIt, ++k)
	{
		labelArray[k] = (short)imgIt.Get();
	}

	return 0;
}

// pre processing the water shed image
//


// post processing the water shed result
//
