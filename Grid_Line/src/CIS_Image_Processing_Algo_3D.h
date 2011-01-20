/* CIS_Image_Process_Algo_3D.h
 * ----------------------------------------------------------------------------
 * Image processing algorithms for 3D images
 *
 * Platforms:
 *  Win32     under MSVC++ 6.0
 *
 *
 * Author:
 *    Jianhua Yao
 *
 * Copyright (C) The Johns Hopkins University, CIS Lab, 2000
 * ----------------------------------------------------------------------------
 *
 * Modification History:
 *  [IM-AA] : 2/99  JHY Initial Creation - modified from Image.hpp in the Image Class Lib v2.0
 *  [IM-AB] : 2/00  JHY modification for CIS2 distibution
 */

#ifndef _CIS_IMAGE_PROCESS_ALGO_3D_H_
#define _CIS_IMAGE_PROCESS_ALGO_3D_H_

#include <CIS_Array_Image3D.h>
#include <CIS_Vector_Image3D.h>

template <class T> CIS_Array_Image3D<T> *
CIS_IPA_3D_PixelWiseAdd(CIS_Array_Image3D<T> *im1, CIS_Array_Image3D<T> *im2);
template <class T> CIS_Array_Image3D<T> *
CIS_IPA_3D_PixelWiseSubtract(CIS_Array_Image3D<T> *im1, CIS_Array_Image3D<T> *im2);
template <class T> CIS_Array_Image3D<T> *
CIS_IPA_3D_PixelWiseMultiply(CIS_Array_Image3D<T> *im1, CIS_Array_Image3D<T> *im2);
template <class T> CIS_Array_Image3D<T> *
CIS_IPA_3D_PixelWiseDivide(CIS_Array_Image3D<T> *im1, CIS_Array_Image3D<T> *im2);
template <class T> CIS_Array_Image3D<T> *
CIS_IPA_3D_PixelWiseAnd(CIS_Array_Image3D<T> *im1, CIS_Array_Image3D<T> *im2);
template <class T> CIS_Array_Image3D<T> *
CIS_IPA_3D_PixelWiseOr(CIS_Array_Image3D<T> *im1, CIS_Array_Image3D<T> *im2);
template <class T> CIS_Array_Image3D<T> *
CIS_IPA_3D_PixelWiseXor(CIS_Array_Image3D<T> *im1, CIS_Array_Image3D<T> *im2);
template <class T> CIS_Array_Image3D<T> *
CIS_IPA_3D_ValueAdd(CIS_Array_Image3D<T> *im1, T value);
template <class T> CIS_Array_Image3D<T> *
CIS_IPA_3D_ValueSubtract(CIS_Array_Image3D<T> *im1, T value);
template <class T> CIS_Array_Image3D<T> *
CIS_IPA_3D_ValueMultiply(CIS_Array_Image3D<T> *im1, T value);
template <class T> CIS_Array_Image3D<T> *
CIS_IPA_3D_ValueDivide(CIS_Array_Image3D<T> *im1, T value);
template <class T> CIS_Array_Image3D<T> *
CIS_IPA_3D_ValueAnd(CIS_Array_Image3D<T> *im1, T value);
template <class T> CIS_Array_Image3D<T> *
CIS_IPA_3D_ValueOr(CIS_Array_Image3D<T> *im1, T value);
template <class T> CIS_Array_Image3D<T> *
CIS_IPA_3D_ValueXor(CIS_Array_Image3D<T> *im1, T value);


template <class T> CIS_ImageErrorCode CIS_IPA_3D_Translate(CIS_Array_Image3D<T> *im, 
					 double tx, double ty, double tz, T fill=0);
template <class T> CIS_ImageErrorCode
CIS_IPA_3D_Rotate(CIS_Array_Image3D<T> *im, Vec3 rot_ang, Vec3 rot_center,T fill=0);
template <class T> CIS_ImageErrorCode
CIS_IPA_3D_Rigid_Transform(CIS_Array_Image3D<T> *im, Frame f, T fill=0);
template <class T> CIS_ImageErrorCode
CIS_IPA_3D_SubSample(CIS_Array_Image3D<T> *im_src, 
					 CIS_Array_Image3D<T> **im_dst, int xs, int ys, int zs);

template <class T> CIS_ImageErrorCode
CIS_IPA_3D_Gradient_Magnitude(CIS_Array_Image3D<T> *im, CIS_Array_Image3D<T> *im_out);
template <class T> CIS_ImageErrorCode
CIS_IPA_3D_Gradient(CIS_Array_Image3D<T> *im, CIS_Vector_Image3D_Vec3 *im_out);


template <class T> CIS_ImageErrorCode
CIS_IPA_3D_Binarize(CIS_Array_Image3D<T> *im, double th, 
				 T loValue, T hiValue);
template <class T> CIS_ImageErrorCode
CIS_IPA_3D_Gaussian_Smooth(CIS_Array_Image3D<T> *im, 
						   int ker_x, int ker_y, int ker_z, double sigma);
template <class T> CIS_ImageErrorCode
CIS_IPA_3D_Gaussian_Smooth(CIS_Array_Image3D<T> *im, 
						   double sigma);
template <class T> CIS_ImageErrorCode
CIS_IPA_3D_Normalize_Z_Interval(CIS_Array_Image3D<T> *im_src, CIS_Array_Image3D<T> **im_dst,
								double newZInterval);

template <class TPixel> CIS_ImageErrorCode
CIS_IPA_3D_Erode(CIS_Array_Image3D<TPixel> *im, int NbIteration, int connect);

template <class TPixel> CIS_ImageErrorCode
CIS_IPA_3D_Dilate(CIS_Array_Image3D<TPixel> *im, int NbIteration, int connect);


template <class _DATA_REP_IM3D_> 
CIS_ImageErrorCode
CIS_IPA_3D_AddGaussianNoise(CIS_Array_Image3D<_DATA_REP_IM3D_> *im, _DATA_REP_IM3D_ magnitude, long seed=0, 
						 _DATA_REP_IM3D_ lower_limit=0, _DATA_REP_IM3D_ upper_limit=0);

#include "CIS_Image_Processing_Algo_3D.txx"
#endif
