#include "NIH_ITK_Utility.h"


template <class T1>
void Copy_CISImage_to_ITKImage(CIS_Array_Image3D<T1> *cisImage, itkImage_3D_float_pointer itk_image)
{

	itkImage_3D_float::RegionType itk_region;
	itkImage_3D_float::SizeType itk_size;

	itk_region = itk_image->GetLargestPossibleRegion();

	itk_size = itk_region.GetSize();

	int sizex, sizey, sizez, k;
	T1 *cisArray;
	sizex = cisImage->Num_Cols();
	sizey = cisImage->Num_Rows();
	sizez = cisImage->Num_Levels();
	cisArray = cisImage->GetArray();

	if(sizex!=itk_size[0] || sizey!=itk_size[1] || sizez!=itk_size[2])
	{
		itk_size[0] = sizex;
		itk_size[1] = sizey;
		itk_size[2] = sizez;
		itk_region.SetSize(itk_size);
		itk_image->SetRegions( itk_region );
		itk_image->Allocate();
	}

	IteratorType_float imgIt(itk_image, itk_region);

	for(imgIt.GoToBegin(), k=0; !imgIt.IsAtEnd(); ++imgIt, ++k)
	{
		imgIt.Set((float)cisArray[k]);
	}

	// set spacing
	double spacing[3];
	spacing[0] = cisImage->Get_Pixel_SizeX();
	spacing[1] = cisImage->Get_Pixel_SizeY();
	spacing[2] = cisImage->Get_Pixel_SizeZ();
	itk_image->SetSpacing(spacing);

	return;
}


template <class T1>
void Copy_CISImage_to_ITKImage_2D(CIS_Array_Image2D<T1> *cisImage, itkImage_2D_float_pointer itk_image)
{

	itkImage_2D_float::RegionType itk_region;
	itkImage_2D_float::SizeType itk_size;

	itk_region = itk_image->GetLargestPossibleRegion();

	itk_size = itk_region.GetSize();

	int sizex, sizey, k;
	T1 *cisArray;
	sizex = cisImage->Num_Cols();
	sizey = cisImage->Num_Rows();
	cisArray = cisImage->GetArray();

	if(sizex!=itk_size[0] || sizey!=itk_size[1])
	{
		itk_size[0] = sizex;
		itk_size[1] = sizey;
		itk_region.SetSize(itk_size);
		itk_image->SetRegions( itk_region );
		itk_image->Allocate();
	}

	IteratorType_float_2D imgIt(itk_image, itk_region);

	for(imgIt.GoToBegin(), k=0; !imgIt.IsAtEnd(); ++imgIt, ++k)
	{
		imgIt.Set((float)cisArray[k]);
	}

	// set spacing
//	double spacing[2];
///	spacing[0] = cisImage->Get_Pixel_SizeX();
///	spacing[1] = cisImage->Get_Pixel_SizeY();
///	itk_image->SetSpacing(spacing);

	return;
}


template <class T1>
void Copy_ITKImage_to_CISImage(itkImage_3D_float_pointer itk_image, CIS_Array_Image3D<T1> *cisImage)
{
	itkImage_3D_float::RegionType itk_region;
	itkImage_3D_float::SizeType itk_size;

	itk_region = itk_image->GetLargestPossibleRegion();

	itk_size = itk_region.GetSize();

	int sizex, sizey, sizez, k;
	T1 *cisArray;
	sizex = cisImage->Num_Cols();
	sizey = cisImage->Num_Rows();
	sizez = cisImage->Num_Levels();

	if(sizex!=itk_size[0] || sizey!=itk_size[1] || sizez!=itk_size[2])
	{
		sizex = itk_size[0];
		sizey = itk_size[1];
		sizez = itk_size[2];
		cisImage->SetSize(sizex, sizey, sizez);
	}
	cisArray = cisImage->GetArray();

	ConstIteratorType_float imgIt(itk_image, itk_region);

	for(imgIt.GoToBegin(), k=0; !imgIt.IsAtEnd(); ++imgIt, ++k)
	{
		cisArray[k] = (T1)imgIt.Get();
	}


	// set spacing
	itkImage_3D_float::SpacingType spacing;
	spacing = itk_image->GetSpacing();
	cisImage->Set_Pixel_SizeX(spacing[0]);
	cisImage->Set_Pixel_SizeY(spacing[1]);
	cisImage->Set_Pixel_SizeZ(spacing[2]);

	return;
}


template <class T1>
void Copy_ITKImage_to_CISImage_2D(itkImage_2D_float_pointer itk_image, CIS_Array_Image2D<T1> *cisImage)
{
	itkImage_2D_float::RegionType itk_region;
	itkImage_2D_float::SizeType itk_size;

	itk_region = itk_image->GetLargestPossibleRegion();

	itk_size = itk_region.GetSize();

	int sizex, sizey, k;
	T1 *cisArray;
	sizex = cisImage->Num_Cols();
	sizey = cisImage->Num_Rows();

	if(sizex!=itk_size[0] || sizey!=itk_size[1])
	{
		sizex = itk_size[0];
		sizey = itk_size[1];
		cisImage->SetSize(sizex, sizey);
	}
	cisArray = cisImage->GetArray();

	ConstIteratorType_float_2D imgIt(itk_image, itk_region);

	for(imgIt.GoToBegin(), k=0; !imgIt.IsAtEnd(); ++imgIt, ++k)
	{
		cisArray[k] = (T1)imgIt.Get();
	}


	// set spacing
///	itkImage_2D_float::SpacingType spacing;
///	spacing = itk_image->GetSpacing();
///	cisImage->Set_Pixel_SizeX(spacing[0]);
///	cisImage->Set_Pixel_SizeY(spacing[1]);

	return;
}



template void Copy_CISImage_to_ITKImage(CIS_Array_Image3D_float *cisImage, itkImage_3D_float_pointer itkImage);
template void Copy_CISImage_to_ITKImage(CIS_Array_Image3D_short *cisImage, itkImage_3D_float_pointer itkImage);
template void Copy_ITKImage_to_CISImage(itkImage_3D_float_pointer itkImage, CIS_Array_Image3D_short *cisImage);
template void Copy_ITKImage_to_CISImage(itkImage_3D_float_pointer itkImage, CIS_Array_Image3D_float *cisImage);

template void Copy_CISImage_to_ITKImage_2D(CIS_Array_Image2D_short *cisImage, itkImage_2D_float_pointer itkImage);
template void Copy_ITKImage_to_CISImage_2D(itkImage_2D_float_pointer itkImage, CIS_Array_Image2D_short *cisImage);
