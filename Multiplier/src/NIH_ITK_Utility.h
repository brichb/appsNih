#ifndef _NIH_ITK_Utility_h
#define _NIH_ITK_Utility_h
#include <itkImage.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkCovariantVector.h> 

#include <CIS_Array_Image3D.h>

typedef    itk::CovariantVector< float > itkVectPixelFloat;

typedef itk::Image<itkVectPixelFloat, 3> itkVectImageFloat;
typedef itkVectImageFloat::Pointer itkVectImageFloatPointer;
typedef itk::ImageRegionConstIterator< itkVectImageFloat > VectConstIteratorFloat;
//typedef itk::ImageRegionIterator< itkVectImageFloat>       IteratorType_float_2D;

typedef itk::Image<float, 3> itkImage_3D_float;
typedef itkImage_3D_float::Pointer itkImage_3D_float_pointer;
typedef itk::ImageRegionConstIterator< itkImage_3D_float > ConstIteratorType_float;
typedef itk::ImageRegionIterator< itkImage_3D_float>       IteratorType_float;

typedef itk::Image<float, 2> itkImage_2D_float;
typedef itkImage_2D_float::Pointer itkImage_2D_float_pointer;
typedef itk::ImageRegionConstIterator< itkImage_2D_float > ConstIteratorType_float_2D;
typedef itk::ImageRegionIterator< itkImage_2D_float>       IteratorType_float_2D;

template <class T1>
void Copy_CISImage_to_ITKImage(CIS_Array_Image3D<T1> *cisImage, itkImage_3D_float_pointer itk_image);
template <class T1>
void Copy_ITKImage_to_CISImage(itkImage_3D_float_pointer itk_image, CIS_Array_Image3D<T1> *cisImage);

template <class T1>
void Copy_CISImage_to_ITKImage_2D(CIS_Array_Image2D<T1> *cisImage, itkImage_2D_float_pointer itk_image);
template <class T1>
void Copy_ITKImage_to_CISImage_2D(itkImage_2D_float_pointer itk_image, CIS_Array_Image2D<T1> *cisImage);

#endif