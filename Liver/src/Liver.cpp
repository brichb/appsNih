#include "NIH_ITK_Utility.h"
#include "NIH_Algo_Region_Growing.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkMultiplyImageFilter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"
#include "itkGeodesicActiveContourLevelSetImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkFastMarchingImageFilter.h"
#include "itkLaplacianSegmentationLevelSetImageFilter.h"
#include "itkCannySegmentationLevelSetImageFilter.h"
#include "itkFastMarchingImageFilter.h"

CIS_Array_Image3D_short *img3Dtemp;
CIS_Array_Image3D_short *img3Dp1;
CIS_Array_Image3D_short *img3Dp2;
CIS_Array_Image3D_short *img3Dorig;
CIS_Array_Image3D_short *img3Doi;
CIS_Array_Image3D_short *img3Dnr;
CIS_Array_Image3D_short *img3Dnr_1;
CIS_Array_Image3D_short *img3Dnr_2;
CIS_Array_Image3D_short *img3Dnr_3;
CIS_Array_Image3D_short *img3Dnr_4;
CIS_Array_Image3D_short *img3Dnr_5;
CIS_Array_Image3D_short *img3Dnr_6;
CIS_Array_Image3D_short *img3Dnr_7;
CIS_Array_Image3D_short *img3Dnr_8;
CIS_Array_Image3D_short *img3Doi_1;
CIS_Array_Image3D_short *img3Doi_2;
CIS_Array_Image3D_short *img3Doi_3;
CIS_Array_Image3D_short *img3Doi_4;
CIS_Array_Image3D_short *img3Doi_5;
CIS_Array_Image3D_short *img3Doi_6;
CIS_Array_Image3D_short *img3Doi_7;
CIS_Array_Image3D_short *img3Doi_8;
CIS_Array_Image3D_short *img3Dpad;

void AnisotropicDiffusion(CIS_Array_Image3D_short *img3Din, CIS_Array_Image3D_short *img3Dout);
void HomogeneityFilter(CIS_Array_Image3D_short *img3Din1, CIS_Array_Image3D_short *img3Din2);
double GetVolume(CIS_Array_Image3D_short *img3D1);
void NonrigidFix(CIS_Array_Image3D_short *img3Din2);
void CutHeart(CIS_Array_Image3D_short *img3Din, CIS_Array_Image3D_short *img3Dout);
void MaxComp(CIS_Array_Image3D_short *img3Din, CIS_Array_Image3D_short *img3Dout);
void HoleFill(CIS_Array_Image3D_short *img3Din);

struct ipixel{
	int value;
	int x;
	int y;
	int z;
};

void main(int argc, char **argv){

	if(argc != 4){
		std::cout<<"Usage: _appsNih_Liver.exe <Original_Image 1.hdr> <nonrigid result.hdr> <output file name.hdr>"<<std::endl;
		exit(1);
	}
	
	
	double nonrigid_volume = 0;
	double hom_volume_step1 = 0;
	double hom_volume_step2 = 0;
	double Thresh_rerun = 100000000;
	double volume_diff = 0;
	double volume_percent = 0;
	
	img3Doi = new CIS_Array_Image3D_short();
	img3Doi->Load(argv[1]);
	img3Dnr = new CIS_Array_Image3D_short();
	img3Dnr->Load(argv[2]);
	
	double sizex=img3Doi->Num_Cols();
	double sizey=img3Doi->Num_Rows();
	double sizez=img3Doi->Num_Levels();
	double pixelx=img3Doi->Get_Pixel_SizeX();
	double pixely=img3Doi->Get_Pixel_SizeY(); 
	double pixelz=img3Doi->Get_Pixel_SizeZ();

	int chunks = 0;
	int redo_homogen_run = -1;
	nonrigid_volume = GetVolume(img3Dnr);

	if (pixelz <= 1)
	{
		chunks = 8;
		int zlevel = sizez / chunks;
		int zremain = 7*zlevel;
		int zlevel_last = sizez - zremain;
		img3Doi_1 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
		img3Doi_2 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
		img3Doi_3 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
		img3Doi_4 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
		img3Doi_5 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
		img3Doi_6 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
		img3Doi_7 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
		img3Doi_8 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel_last);
		img3Doi_1->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_1->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_1->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		img3Doi_2->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_2->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_2->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		img3Doi_3->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_3->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_3->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		img3Doi_4->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_4->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_4->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		img3Doi_5->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_5->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_5->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		img3Doi_6->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_6->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_6->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		img3Doi_7->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_7->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_7->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		img3Doi_8->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_8->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_8->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = 0; z<zlevel; z++)
				{
					img3Doi_1->FastSet(x,y,z,img3Doi->FastGet(x,y,z));
				}
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = zlevel; z<(2*zlevel); z++)
				{
					img3Doi_2->FastSet(x,y,z-zlevel,img3Doi->FastGet(x,y,z));
				}
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = (2*zlevel); z<(3*zlevel); z++)
				{
					img3Doi_3->FastSet(x,y,z-(2*zlevel),img3Doi->FastGet(x,y,z));
				}
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = (3*zlevel); z<(4*zlevel); z++)
				{
					img3Doi_4->FastSet(x,y,z-(3*zlevel),img3Doi->FastGet(x,y,z));
				}
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = (4*zlevel); z<(5*zlevel); z++)
				{
					img3Doi_5->FastSet(x,y,z-(4*zlevel),img3Doi->FastGet(x,y,z));
				}
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = (5*zlevel); z<(6*zlevel); z++)
				{
					img3Doi_6->FastSet(x,y,z-(5*zlevel),img3Doi->FastGet(x,y,z));
				}
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = (6*zlevel); z<(7*zlevel); z++)
				{
					img3Doi_7->FastSet(x,y,z-(6*zlevel),img3Doi->FastGet(x,y,z));
				}
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = (7*zlevel); z<img3Doi->Num_Levels(); z++)
				{
					img3Doi_8->FastSet(x,y,z-(7*zlevel),img3Doi->FastGet(x,y,z));
				}
		std::cout<<"Smoothing Started 1 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_1, img3Doi_1);
		std::cout<<"Smoothing Started 2 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_2, img3Doi_2);
		std::cout<<"Smoothing Started 3 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_3, img3Doi_3);
		std::cout<<"Smoothing Started 4 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_4, img3Doi_4);
		std::cout<<"Smoothing Started 5 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_5, img3Doi_5);
		std::cout<<"Smoothing Started 6 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_6, img3Doi_6);
		std::cout<<"Smoothing Started 7 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_7, img3Doi_7);
		std::cout<<"Smoothing Started 8 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_8, img3Doi_8);



		while(redo_homogen_run != 0)
		{

			if (redo_homogen_run == 2)
			{
				delete img3Dtemp;
			}

			img3Dnr_1 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
			img3Dnr_2 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
			img3Dnr_3 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
			img3Dnr_4 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
			img3Dnr_5 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
			img3Dnr_6 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
			img3Dnr_7 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
			img3Dnr_8 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel_last);

			img3Dnr_1->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_1->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_1->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
			img3Dnr_2->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_2->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_2->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
			img3Dnr_3->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_3->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_3->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
			img3Dnr_4->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_4->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_4->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
			img3Dnr_5->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_5->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_5->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
			img3Dnr_6->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_6->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_6->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
			img3Dnr_7->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_7->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_7->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
			img3Dnr_8->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_8->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_8->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());

			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = 0; z<zlevel; z++)
					{
						img3Dnr_1->FastSet(x,y,z,img3Dnr->FastGet(x,y,z));
					}
			Write_Analyze_File("img3Dnr_1.hdr", *img3Dnr_1);
			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = zlevel; z<(2*zlevel); z++)
					{
						img3Dnr_2->FastSet(x,y,z-zlevel,img3Dnr->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = (2*zlevel); z<(3*zlevel); z++)
					{
						img3Dnr_3->FastSet(x,y,z-(2*zlevel),img3Dnr->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = (3*zlevel); z<(4*zlevel); z++)
					{
						img3Dnr_4->FastSet(x,y,z-(3*zlevel),img3Dnr->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = (4*zlevel); z<(5*zlevel); z++)
					{
						img3Dnr_5->FastSet(x,y,z-(4*zlevel),img3Dnr->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = (5*zlevel); z<(6*zlevel); z++)
					{
						img3Dnr_6->FastSet(x,y,z-(5*zlevel),img3Dnr->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = (6*zlevel); z<(7*zlevel); z++)
					{
						img3Dnr_7->FastSet(x,y,z-(6*zlevel),img3Dnr->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = (7*zlevel); z<img3Doi->Num_Levels(); z++)
					{
						img3Dnr_8->FastSet(x,y,z-(7*zlevel),img3Dnr->FastGet(x,y,z));
					}
	//		if (redo_homogen_run != 2)
	//		{
	//			std::cout<<"Attempting to get rid of the heart. . ."<<std::endl;
	//			CutHeart(img3Dnr_1,img3Dnr_1);
	//		}

			img3Dnr_1->Set_Pixel_SizeX(img3Dnr->Get_Pixel_SizeX());
			img3Dnr_1->Set_Pixel_SizeY(img3Dnr->Get_Pixel_SizeY());
			img3Dnr_1->Set_Pixel_SizeZ(img3Dnr->Get_Pixel_SizeZ());

			nonrigid_volume = GetVolume(img3Dnr_1);

			if(redo_homogen_run != 2)
			{
				std::cout<<"HomogeneityFilter Started for Part 1  ..."<<std::endl;
				HomogeneityFilter(img3Doi_1, img3Dnr_1);
				std::cout<<"HomogeneityFilter Started for Part 2  ..."<<std::endl;
				HomogeneityFilter(img3Doi_2, img3Dnr_2);
				std::cout<<"HomogeneityFilter Started for Part 3  ..."<<std::endl;
				HomogeneityFilter(img3Doi_3, img3Dnr_3);
				std::cout<<"HomogeneityFilter Started for Part 4  ..."<<std::endl;
				HomogeneityFilter(img3Doi_4, img3Dnr_4);
				std::cout<<"HomogeneityFilter Started for Part 5  ..."<<std::endl;
				HomogeneityFilter(img3Doi_5, img3Dnr_5);
				std::cout<<"HomogeneityFilter Started for Part 6  ..."<<std::endl;
				HomogeneityFilter(img3Doi_6, img3Dnr_6);
				std::cout<<"HomogeneityFilter Started for Part 7  ..."<<std::endl;
				HomogeneityFilter(img3Doi_7, img3Dnr_7);
				std::cout<<"HomogeneityFilter Started for Part 8  ..."<<std::endl;
				HomogeneityFilter(img3Doi_8, img3Dnr_8);
			}			
			if(redo_homogen_run == 2)
			{
				std::cout<<"Fixing Nonrigid Result 1  ..."<<std::endl;
				NonrigidFix(img3Dnr_1);
				std::cout<<"Fixing Nonrigid Result 2  ..."<<std::endl;
				NonrigidFix(img3Dnr_2);
				std::cout<<"Fixing Nonrigid Result 3  ..."<<std::endl;
				NonrigidFix(img3Dnr_3);
				std::cout<<"Fixing Nonrigid Result 4  ..."<<std::endl;
				NonrigidFix(img3Dnr_4);
				std::cout<<"Fixing Nonrigid Result 5  ..."<<std::endl;
				NonrigidFix(img3Dnr_5);
				std::cout<<"Fixing Nonrigid Result 6  ..."<<std::endl;
				NonrigidFix(img3Dnr_6);
				std::cout<<"Fixing Nonrigid Result 7  ..."<<std::endl;
				NonrigidFix(img3Dnr_7);
				std::cout<<"Fixing Nonrigid Result 8  ..."<<std::endl;
				NonrigidFix(img3Dnr_8);
			}
			
			img3Dtemp = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), sizez);
			img3Dtemp->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dtemp->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dtemp->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());

			for (int x = 0; x<img3Dnr_1->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_1->Num_Rows(); y++)
					for (int z = 0; z<zlevel; z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_1->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Dnr_2->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_2->Num_Rows(); y++)
					for (int z = zlevel; z<(2*zlevel); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_2->FastGet(x,y,z - zlevel));
					}
			for (int x = 0; x<img3Dnr_3->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_3->Num_Rows(); y++)
					for (int z = (2*zlevel); z<(3*zlevel); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_3->FastGet(x,y,z-(2*zlevel)));
					}
			for (int x = 0; x<img3Dnr_4->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_4->Num_Rows(); y++)
					for (int z = (3*zlevel); z<(4*zlevel); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_4->FastGet(x,y,z-(3*zlevel)));
					}
			for (int x = 0; x<img3Dnr_5->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_5->Num_Rows(); y++)
					for (int z = (4*zlevel); z<(5*zlevel); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_5->FastGet(x,y,z-(4*zlevel)));
					}
			for (int x = 0; x<img3Dnr_6->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_6->Num_Rows(); y++)
					for (int z = (5*zlevel); z<(6*zlevel); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_6->FastGet(x,y,z-(5*zlevel)));
					}
			for (int x = 0; x<img3Dnr_7->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_7->Num_Rows(); y++)
					for (int z = (6*zlevel); z<(7*zlevel); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_7->FastGet(x,y,z-(6*zlevel)));
					}
			for (int x = 0; x<img3Dnr_8->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_8->Num_Rows(); y++)
					for (int z = (7*zlevel); z<img3Dnr->Num_Levels(); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_8->FastGet(x,y,z-(7*zlevel)));
					}
			Write_Analyze_File("AfterHomogeneity.hdr",*img3Dtemp);

			delete img3Dnr_1;
			delete img3Dnr_2;
			delete img3Dnr_3;
			delete img3Dnr_4;
			delete img3Dnr_5;
			delete img3Dnr_6;
			delete img3Dnr_7;
			delete img3Dnr_8;
			delete img3Doi_1;
			delete img3Doi_2;
			delete img3Doi_3;
			delete img3Doi_4;
			delete img3Doi_5;
			delete img3Doi_6;
			delete img3Doi_7;
			delete img3Doi_8;

/*
//////////////////////////// Splitting up the MaxComp, have to or won't run, it will run out of memory
			int zlevelmax = img3Dtemp->Num_Levels() / 2;
			int zlevel_lastmax = img3Dtemp->Num_Levels() - zlevelmax;
			img3Dp1 = new CIS_Array_Image3D_short(img3Dtemp->Get_SizeX(), img3Dtemp->Get_SizeY(), zlevelmax);
			img3Dp2 = new CIS_Array_Image3D_short(img3Dtemp->Get_SizeX(), img3Dtemp->Get_SizeY(), zlevel_lastmax);

			img3Dp1->Set_Pixel_SizeX(img3Dtemp->Get_Pixel_SizeX());
			img3Dp1->Set_Pixel_SizeY(img3Dtemp->Get_Pixel_SizeY());
			img3Dp1->Set_Pixel_SizeZ(img3Dtemp->Get_Pixel_SizeZ());
			img3Dp2->Set_Pixel_SizeX(img3Dtemp->Get_Pixel_SizeX());
			img3Dp2->Set_Pixel_SizeY(img3Dtemp->Get_Pixel_SizeY());
			img3Dp2->Set_Pixel_SizeZ(img3Dtemp->Get_Pixel_SizeZ());


			for (int x = 0; x<img3Dtemp->Num_Cols(); x++)
				for (int y = 0; y<img3Dtemp->Num_Rows(); y++)
					for (int z = 0; z<zlevelmax; z++)
					{
						img3Dp1->FastSet(x,y,z,img3Dtemp->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Dtemp->Num_Cols(); x++)
				for (int y = 0; y<img3Dtemp->Num_Rows(); y++)
					for (int z = zlevelmax; z<img3Dtemp->Num_Levels(); z++)
					{
						img3Dp2->FastSet(x,y,z-zlevelmax,img3Dtemp->FastGet(x,y,z));
					}
			MaxComp(img3Dp1, img3Dp1);
			MaxComp(img3Dp2, img3Dp2);
			img3Dtemp->~CIS_Array_Image3D();
			img3Dtemp = new CIS_Array_Image3D_short(img3Dnr->Get_SizeX(), img3Dnr->Get_SizeY(), img3Dnr->Num_Levels());
			img3Dtemp->Set_Pixel_SizeX(img3Dnr->Get_Pixel_SizeX());
			img3Dtemp->Set_Pixel_SizeY(img3Dnr->Get_Pixel_SizeY());
			img3Dtemp->Set_Pixel_SizeZ(img3Dnr->Get_Pixel_SizeZ());

			for (int x = 0; x<img3Dnr->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr->Num_Rows(); y++)
					for (int z = 0; z<zlevelmax; z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dp1->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Dnr->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr->Num_Rows(); y++)
					for (int z = zlevelmax; z<img3Dnr->Num_Levels(); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dp2->FastGet(x,y,z - zlevelmax));
					}

					////////////////////////////////////
			img3Dp1->~CIS_Array_Image3D();
			img3Dp2->~CIS_Array_Image3D();
			*/
			MaxComp(img3Dtemp, img3Dtemp);
			img3Dtemp->Set_Pixel_SizeX(img3Dnr->Get_Pixel_SizeX());
			img3Dtemp->Set_Pixel_SizeY(img3Dnr->Get_Pixel_SizeY());
			img3Dtemp->Set_Pixel_SizeZ(img3Dnr->Get_Pixel_SizeZ());
			Write_Analyze_File("AfterMaxComponents.hdr",*img3Dtemp);
			HoleFill(img3Dtemp);

			if(redo_homogen_run != 2)
			{
				hom_volume_step1 = GetVolume(img3Dtemp);
				volume_diff = abs(hom_volume_step1 - nonrigid_volume);
				volume_percent = 100 * (volume_diff / nonrigid_volume);

				for (int x = 0; x<img3Dtemp->Num_Cols(); x++)
					for (int y = 0; y<img3Dtemp->Num_Rows(); y++)
						for (int z = 0; z<img3Dtemp->Num_Levels(); z++)
						{
							if (img3Dtemp->FastGet(x,y,z) > 0)
							{
								img3Dtemp->FastSet(x,y,z,1000);
							}
						}
			}


			if (redo_homogen_run == 2)
			{
				img3Dtemp->Set_Pixel_SizeX(img3Dnr->Get_Pixel_SizeX());
				img3Dtemp->Set_Pixel_SizeY(img3Dnr->Get_Pixel_SizeY());
				img3Dtemp->Set_Pixel_SizeZ(img3Dnr->Get_Pixel_SizeZ());
				Write_Analyze_File(argv[3],*img3Dtemp);
				redo_homogen_run = 0;
			}
			else if ((volume_percent < Thresh_rerun) && redo_homogen_run != 2)
			{

				std::cout<<"First Homogeneity Worked!!!!"<<std::endl;
				std::cout<<"Nonrig vol: "<<nonrigid_volume<<"cm^3"<<std::endl;
				std::cout<<"Homogen vol: "<<hom_volume_step1<<"cm^3"<<std::endl;
				std::cout<<"Difference: "<<volume_percent<<"%"<<std::endl;
				img3Dtemp->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
				img3Dtemp->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
				img3Dtemp->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
				Write_Analyze_File(argv[3],*img3Dtemp);
				redo_homogen_run = 0;
			}
			else
			{

				std::cout<<"First Homogeneity Didn't Work :( . . . Using NonRigid Result"<<std::endl;
				std::cout<<"Nonrig vol: "<<nonrigid_volume<<"cm^3"<<std::endl;
				std::cout<<"Homogen vol: "<<hom_volume_step1<<"cm^3"<<std::endl;
				std::cout<<"Difference: "<<volume_percent<<"%"<<std::endl;
				redo_homogen_run = 2;
			}
		}
	}
	else if (pixelz > 1 && pixelz < 2)
	{
		chunks = 6;
		int zlevel = sizez / chunks;
		int zremain = 5*zlevel;
		int zlevel_last = sizez - zremain;
		img3Doi_1 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
		img3Doi_2 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
		img3Doi_3 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
		img3Doi_4 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
		img3Doi_5 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
		img3Doi_6 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel_last);
		img3Doi_1->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_1->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_1->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		img3Doi_2->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_2->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_2->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		img3Doi_3->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_3->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_3->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		img3Doi_4->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_4->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_4->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		img3Doi_5->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_5->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_5->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		img3Doi_6->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_6->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_6->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = 0; z<zlevel; z++)
				{
					img3Doi_1->FastSet(x,y,z,img3Doi->FastGet(x,y,z));
				}
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = zlevel; z<(2*zlevel); z++)
				{
					img3Doi_2->FastSet(x,y,z-zlevel,img3Doi->FastGet(x,y,z));
				}
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = (2*zlevel); z<(3*zlevel); z++)
				{
					img3Doi_3->FastSet(x,y,z-(2*zlevel),img3Doi->FastGet(x,y,z));
				}
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = (3*zlevel); z<(4*zlevel); z++)
				{
					img3Doi_4->FastSet(x,y,z-(3*zlevel),img3Doi->FastGet(x,y,z));
				}
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = (4*zlevel); z<(5*zlevel); z++)
				{
					img3Doi_5->FastSet(x,y,z-(4*zlevel),img3Doi->FastGet(x,y,z));
				}
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = (5*zlevel); z<img3Doi->Num_Levels(); z++)
				{
					img3Doi_6->FastSet(x,y,z-(5*zlevel),img3Doi->FastGet(x,y,z));
				}

		std::cout<<"Smoothing Started 1 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_1, img3Doi_1);
		std::cout<<"Smoothing Started 2 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_2, img3Doi_2);
		std::cout<<"Smoothing Started 3 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_3, img3Doi_3);
		std::cout<<"Smoothing Started 4 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_4, img3Doi_4);
		std::cout<<"Smoothing Started 5 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_5, img3Doi_5);
		std::cout<<"Smoothing Started 6 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_6, img3Doi_6);

		while(redo_homogen_run != 0)
		{

			if (redo_homogen_run == 2)
			{
				delete img3Dtemp;
			}

			img3Dnr_1 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
			img3Dnr_2 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
			img3Dnr_3 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
			img3Dnr_4 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
			img3Dnr_5 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
			img3Dnr_6 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel_last);

			img3Dnr_1->Set_Pixel_SizeX(img3Dnr->Get_Pixel_SizeX());
			img3Dnr_1->Set_Pixel_SizeY(img3Dnr->Get_Pixel_SizeY());
			img3Dnr_1->Set_Pixel_SizeZ(img3Dnr->Get_Pixel_SizeZ());
			img3Dnr_2->Set_Pixel_SizeX(img3Dnr->Get_Pixel_SizeX());
			img3Dnr_2->Set_Pixel_SizeY(img3Dnr->Get_Pixel_SizeY());
			img3Dnr_2->Set_Pixel_SizeZ(img3Dnr->Get_Pixel_SizeZ());
			img3Dnr_3->Set_Pixel_SizeX(img3Dnr->Get_Pixel_SizeX());
			img3Dnr_3->Set_Pixel_SizeY(img3Dnr->Get_Pixel_SizeY());
			img3Dnr_3->Set_Pixel_SizeZ(img3Dnr->Get_Pixel_SizeZ());
			img3Dnr_4->Set_Pixel_SizeX(img3Dnr->Get_Pixel_SizeX());
			img3Dnr_4->Set_Pixel_SizeY(img3Dnr->Get_Pixel_SizeY());
			img3Dnr_4->Set_Pixel_SizeZ(img3Dnr->Get_Pixel_SizeZ());
			img3Dnr_5->Set_Pixel_SizeX(img3Dnr->Get_Pixel_SizeX());
			img3Dnr_5->Set_Pixel_SizeY(img3Dnr->Get_Pixel_SizeY());
			img3Dnr_5->Set_Pixel_SizeZ(img3Dnr->Get_Pixel_SizeZ());
			img3Dnr_6->Set_Pixel_SizeX(img3Dnr->Get_Pixel_SizeX());
			img3Dnr_6->Set_Pixel_SizeY(img3Dnr->Get_Pixel_SizeY());
			img3Dnr_6->Set_Pixel_SizeZ(img3Dnr->Get_Pixel_SizeZ());

			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = 0; z<zlevel; z++)
					{
						img3Dnr_1->FastSet(x,y,z,img3Dnr->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = zlevel; z<(2*zlevel); z++)
					{
						img3Dnr_2->FastSet(x,y,z-zlevel,img3Dnr->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = (2*zlevel); z<(3*zlevel); z++)
					{
						img3Dnr_3->FastSet(x,y,z-(2*zlevel),img3Dnr->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = (3*zlevel); z<(4*zlevel); z++)
					{
						img3Dnr_4->FastSet(x,y,z-(3*zlevel),img3Dnr->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = (4*zlevel); z<(5*zlevel); z++)
					{
						img3Dnr_5->FastSet(x,y,z-(4*zlevel),img3Dnr->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = (5*zlevel); z<img3Doi->Num_Levels(); z++)
					{
						img3Dnr_6->FastSet(x,y,z-(5*zlevel),img3Dnr->FastGet(x,y,z));
					}

	//		if (redo_homogen_run != 2)
	//		{
		//		std::cout<<"Attempting to get rid of the heart. . ."<<std::endl;
	//			CutHeart(img3Dnr_1,img3Dnr_1);
	//		}
			img3Dnr_1->Set_Pixel_SizeX(img3Dnr->Get_Pixel_SizeX());
			img3Dnr_1->Set_Pixel_SizeY(img3Dnr->Get_Pixel_SizeY());
			img3Dnr_1->Set_Pixel_SizeZ(img3Dnr->Get_Pixel_SizeZ());

			if(redo_homogen_run != 2)
			{
				std::cout<<"HomogeneityFilter Started for Part 1  ..."<<std::endl;
				HomogeneityFilter(img3Doi_1, img3Dnr_1);
				std::cout<<"HomogeneityFilter Started for Part 2  ..."<<std::endl;
				HomogeneityFilter(img3Doi_2, img3Dnr_2);
				std::cout<<"HomogeneityFilter Started for Part 3  ..."<<std::endl;
				HomogeneityFilter(img3Doi_3, img3Dnr_3);
				std::cout<<"HomogeneityFilter Started for Part 4  ..."<<std::endl;
				HomogeneityFilter(img3Doi_4, img3Dnr_4);
				std::cout<<"HomogeneityFilter Started for Part 5  ..."<<std::endl;
				HomogeneityFilter(img3Doi_5, img3Dnr_5);
				std::cout<<"HomogeneityFilter Started for Part 6  ..."<<std::endl;
				HomogeneityFilter(img3Doi_6, img3Dnr_6);
			}
			if(redo_homogen_run == 2)
			{
				std::cout<<"Fixing Nonrigid Result 1  ..."<<std::endl;
				NonrigidFix(img3Dnr_1);
				std::cout<<"Fixing Nonrigid Result 2  ..."<<std::endl;
				NonrigidFix(img3Dnr_2);
				std::cout<<"Fixing Nonrigid Result 3  ..."<<std::endl;
				NonrigidFix(img3Dnr_3);
				std::cout<<"Fixing Nonrigid Result 4  ..."<<std::endl;
				NonrigidFix(img3Dnr_4);
				std::cout<<"Fixing Nonrigid Result 5  ..."<<std::endl;
				NonrigidFix(img3Dnr_5);
				std::cout<<"Fixing Nonrigid Result 6  ..."<<std::endl;
				NonrigidFix(img3Dnr_6);
			}
			
			img3Dtemp = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), sizez);
			img3Dtemp->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dtemp->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dtemp->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());

			for (int x = 0; x<img3Dnr_1->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_1->Num_Rows(); y++)
					for (int z = 0; z<zlevel; z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_1->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Dnr_2->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_2->Num_Rows(); y++)
					for (int z = zlevel; z<(2*zlevel); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_2->FastGet(x,y,z - zlevel));
					}
			for (int x = 0; x<img3Dnr_3->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_3->Num_Rows(); y++)
					for (int z = (2*zlevel); z<(3*zlevel); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_3->FastGet(x,y,z-(2*zlevel)));
					}
			for (int x = 0; x<img3Dnr_4->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_4->Num_Rows(); y++)
					for (int z = (3*zlevel); z<(4*zlevel); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_4->FastGet(x,y,z-(3*zlevel)));
					}
			for (int x = 0; x<img3Dnr_5->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_5->Num_Rows(); y++)
					for (int z = (4*zlevel); z<(5*zlevel); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_5->FastGet(x,y,z-(4*zlevel)));
					}
			for (int x = 0; x<img3Dnr_6->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_6->Num_Rows(); y++)
					for (int z = (5*zlevel); z<img3Dnr->Num_Levels(); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_6->FastGet(x,y,z-(5*zlevel)));
					}
			Write_Analyze_File("AfterHomogeneity.hdr",*img3Dtemp);
			
			delete img3Dnr_1;
			delete img3Dnr_2;
			delete img3Dnr_3;
			delete img3Dnr_4;
			delete img3Dnr_5;
			delete img3Dnr_6;
			delete img3Doi_1;
			delete img3Doi_2;
			delete img3Doi_3;
			delete img3Doi_4;
			delete img3Doi_5;
			delete img3Doi_6;

			MaxComp(img3Dtemp, img3Dtemp);
			HoleFill(img3Dtemp);

			img3Dtemp->Set_Pixel_SizeX(img3Dnr->Get_Pixel_SizeX());
			img3Dtemp->Set_Pixel_SizeY(img3Dnr->Get_Pixel_SizeY());
			img3Dtemp->Set_Pixel_SizeZ(img3Dnr->Get_Pixel_SizeZ());

			if(redo_homogen_run != 2)
			{
				hom_volume_step1 = GetVolume(img3Dtemp);
				volume_diff = abs(hom_volume_step1 - nonrigid_volume);
				volume_percent = 100 * (volume_diff / nonrigid_volume);

				for (int x = 0; x<img3Dtemp->Num_Cols(); x++)
					for (int y = 0; y<img3Dtemp->Num_Rows(); y++)
						for (int z = 0; z<img3Dtemp->Num_Levels(); z++)
						{
							if (img3Dtemp->FastGet(x,y,z) > 0)
							{
								img3Dtemp->FastSet(x,y,z,1000);
							}
						}
			}

			if (redo_homogen_run == 2)
			{
				img3Dtemp->Set_Pixel_SizeX(img3Dnr->Get_Pixel_SizeX());
				img3Dtemp->Set_Pixel_SizeY(img3Dnr->Get_Pixel_SizeY());
				img3Dtemp->Set_Pixel_SizeZ(img3Dnr->Get_Pixel_SizeZ());
				Write_Analyze_File(argv[3],*img3Dtemp);
				redo_homogen_run = 0;
			}
			//else if ((volume_percent < Thresh_rerun) && redo_homogen_run != 2)
			else
			{

				std::cout<<"First Homogeneity Worked!!!!"<<std::endl;
				std::cout<<"Nonrig vol: "<<nonrigid_volume<<"cm^3"<<std::endl;
				std::cout<<"Homogen vol: "<<hom_volume_step1<<"cm^3"<<std::endl;
				std::cout<<"Difference: "<<volume_percent<<"%"<<std::endl;
				img3Dtemp->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
				img3Dtemp->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
				img3Dtemp->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
				Write_Analyze_File(argv[3],*img3Dtemp);
				redo_homogen_run = 0;
			}
			/*else
			{

				std::cout<<"First Homogeneity Didn't Work :( . . . Using NonRigid Result"<<std::endl;
				std::cout<<"Nonrig vol: "<<nonrigid_volume<<"cm^3"<<std::endl;
				std::cout<<"Homogen vol: "<<hom_volume_step1<<"cm^3"<<std::endl;
				std::cout<<"Difference: "<<volume_percent<<"%"<<std::endl;
				redo_homogen_run = 2;
			}*/

		}
	}
	else if (pixelz >= 2 && pixelz < 2.5)
	{
		chunks = 5;
		int zlevel = sizez / chunks;
		int zremain = 4*zlevel;
		int zlevel_last = sizez - zremain;
		img3Doi_1 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
		img3Doi_2 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
		img3Doi_3 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
		img3Doi_4 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
		img3Doi_5 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel_last);
		img3Doi_1->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_1->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_1->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		img3Doi_2->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_2->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_2->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		img3Doi_3->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_3->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_3->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		img3Doi_4->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_4->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_4->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		img3Doi_5->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_5->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_5->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = 0; z<zlevel; z++)
				{
					img3Doi_1->FastSet(x,y,z,img3Doi->FastGet(x,y,z));
				}
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = zlevel; z<(2*zlevel); z++)
				{
					img3Doi_2->FastSet(x,y,z-zlevel,img3Doi->FastGet(x,y,z));
				}
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = (2*zlevel); z<(3*zlevel); z++)
				{
					img3Doi_3->FastSet(x,y,z-(2*zlevel),img3Doi->FastGet(x,y,z));
				}
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = (3*zlevel); z<(4*zlevel); z++)
				{
					img3Doi_4->FastSet(x,y,z-(3*zlevel),img3Doi->FastGet(x,y,z));
				}
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = (4*zlevel); z<img3Doi->Num_Levels(); z++)
				{
					img3Doi_5->FastSet(x,y,z-(4*zlevel),img3Doi->FastGet(x,y,z));
				}
		std::cout<<"Smoothing Started 1 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_1, img3Doi_1);
		std::cout<<"Smoothing Started 2 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_2, img3Doi_2);
		std::cout<<"Smoothing Started 3 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_3, img3Doi_3);
		std::cout<<"Smoothing Started 4 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_4, img3Doi_4);
		std::cout<<"Smoothing Started 5 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_5, img3Doi_5);

		while(redo_homogen_run != 0)
		{

			if (redo_homogen_run == 2)
			{
				delete img3Dtemp;
			}

			img3Dnr_1 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
			img3Dnr_2 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
			img3Dnr_3 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
			img3Dnr_4 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
			img3Dnr_5 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel_last);

			img3Dnr_1->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_1->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_1->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
			img3Dnr_2->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_2->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_2->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
			img3Dnr_3->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_3->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_3->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
			img3Dnr_4->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_4->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_4->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
			img3Dnr_5->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_5->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_5->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());


			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = 0; z<zlevel; z++)
					{
						img3Dnr_1->FastSet(x,y,z,img3Dnr->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = zlevel; z<(2*zlevel); z++)
					{
						img3Dnr_2->FastSet(x,y,z-zlevel,img3Dnr->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = (2*zlevel); z<(3*zlevel); z++)
					{
						img3Dnr_3->FastSet(x,y,z-(2*zlevel),img3Dnr->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = (3*zlevel); z<(4*zlevel); z++)
					{
						img3Dnr_4->FastSet(x,y,z-(3*zlevel),img3Dnr->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = (4*zlevel); z<img3Doi->Num_Levels(); z++)
					{
						img3Dnr_5->FastSet(x,y,z-(4*zlevel),img3Dnr->FastGet(x,y,z));
					}

	//		if (redo_homogen_run != 2)
	//		{
	//			std::cout<<"Attempting to get rid of the heart. . ."<<std::endl;
	//			CutHeart(img3Dnr_1,img3Dnr_1);
	//		}

			img3Dnr_1->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_1->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_1->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());

			if(redo_homogen_run != 1 && redo_homogen_run != 2)
			{
				std::cout<<"HomogeneityFilter Started for Part 1  ..."<<std::endl;
				HomogeneityFilter(img3Doi_1, img3Dnr_1);
				std::cout<<"HomogeneityFilter Started for Part 2  ..."<<std::endl;
				HomogeneityFilter(img3Doi_2, img3Dnr_2);
				std::cout<<"HomogeneityFilter Started for Part 3  ..."<<std::endl;
				HomogeneityFilter(img3Doi_3, img3Dnr_3);
				std::cout<<"HomogeneityFilter Started for Part 4  ..."<<std::endl;
				HomogeneityFilter(img3Doi_4, img3Dnr_4);
				std::cout<<"HomogeneityFilter Started for Part 5  ..."<<std::endl;
				HomogeneityFilter(img3Doi_5, img3Dnr_5);
			}
			if(redo_homogen_run == 2)
			{
				std::cout<<"Fixing Nonrigid Result 1  ..."<<std::endl;
				NonrigidFix(img3Dnr_1);
				std::cout<<"Fixing Nonrigid Result 2  ..."<<std::endl;
				NonrigidFix(img3Dnr_2);
				std::cout<<"Fixing Nonrigid Result 3  ..."<<std::endl;
				NonrigidFix(img3Dnr_3);
				std::cout<<"Fixing Nonrigid Result 4  ..."<<std::endl;
				NonrigidFix(img3Dnr_4);
				std::cout<<"Fixing Nonrigid Result 5  ..."<<std::endl;
				NonrigidFix(img3Dnr_5);
			}
			
			img3Dtemp = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), sizez);
			img3Dtemp->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dtemp->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dtemp->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());

			for (int x = 0; x<img3Dnr_1->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_1->Num_Rows(); y++)
					for (int z = 0; z<zlevel; z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_1->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Dnr_2->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_2->Num_Rows(); y++)
					for (int z = zlevel; z<(2*zlevel); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_2->FastGet(x,y,z - zlevel));
					}
			for (int x = 0; x<img3Dnr_3->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_3->Num_Rows(); y++)
					for (int z = (2*zlevel); z<(3*zlevel); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_3->FastGet(x,y,z-(2*zlevel)));
					}
			for (int x = 0; x<img3Dnr_4->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_4->Num_Rows(); y++)
					for (int z = (3*zlevel); z<(4*zlevel); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_4->FastGet(x,y,z-(3*zlevel)));
					}
			for (int x = 0; x<img3Dnr_5->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_5->Num_Rows(); y++)
					for (int z = (4*zlevel); z<img3Dnr->Num_Levels(); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_5->FastGet(x,y,z-(4*zlevel)));
					}
			Write_Analyze_File("AfterHomogeneity.hdr",*img3Dtemp);

			delete img3Dnr_1;
			delete img3Dnr_2;
			delete img3Dnr_3;
			delete img3Dnr_4;
			delete img3Dnr_5;
			delete img3Doi_1;
			delete img3Doi_2;
			delete img3Doi_3;
			delete img3Doi_4;
			delete img3Doi_5;

			MaxComp(img3Dtemp, img3Dtemp);
			HoleFill(img3Dtemp);

			img3Dtemp->Set_Pixel_SizeX(img3Dnr->Get_Pixel_SizeX());
			img3Dtemp->Set_Pixel_SizeY(img3Dnr->Get_Pixel_SizeY());
			img3Dtemp->Set_Pixel_SizeZ(img3Dnr->Get_Pixel_SizeZ());

			if(redo_homogen_run != 2)
			{
				hom_volume_step1 = GetVolume(img3Dtemp);
				volume_diff = abs(hom_volume_step1 - nonrigid_volume);
				volume_percent = 100 * (volume_diff / nonrigid_volume);

				for (int x = 0; x<img3Dtemp->Num_Cols(); x++)
					for (int y = 0; y<img3Dtemp->Num_Rows(); y++)
						for (int z = 0; z<img3Dtemp->Num_Levels(); z++)
						{
							if (img3Dtemp->FastGet(x,y,z) > 0)
							{
								img3Dtemp->FastSet(x,y,z,1000);
							}
						}
			}

			if (redo_homogen_run == 2)
			{
				img3Dtemp->Set_Pixel_SizeX(img3Dnr->Get_Pixel_SizeX());
				img3Dtemp->Set_Pixel_SizeY(img3Dnr->Get_Pixel_SizeY());
				img3Dtemp->Set_Pixel_SizeZ(img3Dnr->Get_Pixel_SizeZ());
				Write_Analyze_File(argv[3],*img3Dtemp);
				redo_homogen_run = 0;
			}
			else if ((volume_percent < Thresh_rerun) && redo_homogen_run != 2)
			{

				std::cout<<"First Homogeneity Worked!!!!"<<std::endl;
				std::cout<<"Nonrig vol: "<<nonrigid_volume<<"cm^3"<<std::endl;
				std::cout<<"Homogen vol: "<<hom_volume_step1<<"cm^3"<<std::endl;
				std::cout<<"Difference: "<<volume_percent<<"%"<<std::endl;
				img3Dtemp->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
				img3Dtemp->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
				img3Dtemp->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
				Write_Analyze_File(argv[3],*img3Dtemp);
				redo_homogen_run = 0;
			}
			else
			{

				std::cout<<"First Homogeneity Didn't Work :( . . . Using NonRigid Result"<<std::endl;
				std::cout<<"Nonrig vol: "<<nonrigid_volume<<"cm^3"<<std::endl;
				std::cout<<"Homogen vol: "<<hom_volume_step1<<"cm^3"<<std::endl;
				std::cout<<"Difference: "<<volume_percent<<"%"<<std::endl;
				redo_homogen_run = 2;
			}

		}
	}
	else if (pixelz >= 2.5 && pixelz < 3)
	{
		chunks = 4;
		int zlevel = sizez / chunks;
		int zremain = 3*zlevel;
		int zlevel_last = sizez - zremain;
		img3Doi_1 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
		img3Doi_2 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
		img3Doi_3 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
		img3Doi_4 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel_last);
		img3Doi_1->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_1->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_1->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		img3Doi_2->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_2->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_2->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		img3Doi_3->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_3->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_3->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		img3Doi_4->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_4->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_4->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = 0; z<zlevel; z++)
				{
					img3Doi_1->FastSet(x,y,z,img3Doi->FastGet(x,y,z));
				}
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = zlevel; z<(2*zlevel); z++)
				{
					img3Doi_2->FastSet(x,y,z-zlevel,img3Doi->FastGet(x,y,z));
				}
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = (2*zlevel); z<(3*zlevel); z++)
				{
					img3Doi_3->FastSet(x,y,z-(2*zlevel),img3Doi->FastGet(x,y,z));
				}
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = (3*zlevel); z<img3Doi->Num_Levels(); z++)
				{
					img3Doi_4->FastSet(x,y,z-(3*zlevel),img3Doi->FastGet(x,y,z));
				}
		std::cout<<"Smoothing Started 1 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_1, img3Doi_1);
		std::cout<<"Smoothing Started 2 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_2, img3Doi_2);
		std::cout<<"Smoothing Started 3 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_3, img3Doi_3);
		std::cout<<"Smoothing Started 4 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_4, img3Doi_4);


		while(redo_homogen_run != 0)
		{

			if (redo_homogen_run == 2)
			{
				delete img3Dtemp;
			}

			img3Dnr_1 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
			img3Dnr_2 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
			img3Dnr_3 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
			img3Dnr_4 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel_last);

			img3Dnr_1->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_1->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_1->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
			img3Dnr_2->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_2->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_2->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
			img3Dnr_3->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_3->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_3->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
			img3Dnr_4->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_4->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_4->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());

			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = 0; z<zlevel; z++)
					{
						img3Dnr_1->FastSet(x,y,z,img3Dnr->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = zlevel; z<(2*zlevel); z++)
					{
						img3Dnr_2->FastSet(x,y,z-zlevel,img3Dnr->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = (2*zlevel); z<(3*zlevel); z++)
					{
						img3Dnr_3->FastSet(x,y,z-(2*zlevel),img3Dnr->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = (3*zlevel); z<img3Doi->Num_Levels(); z++)
					{
						img3Dnr_4->FastSet(x,y,z-(3*zlevel),img3Dnr->FastGet(x,y,z));
					}

	//		if (redo_homogen_run != 2)
	//		{
	//			std::cout<<"Attempting to get rid of the heart. . ."<<std::endl;
	//			CutHeart(img3Dnr_1,img3Dnr_1);
	//		}

			img3Dnr_1->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_1->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_1->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());

			if(redo_homogen_run != 2)
			{
				std::cout<<"HomogeneityFilter Started for Part 1  ..."<<std::endl;
				HomogeneityFilter(img3Doi_1, img3Dnr_1);
				std::cout<<"HomogeneityFilter Started for Part 2  ..."<<std::endl;
				HomogeneityFilter(img3Doi_2, img3Dnr_2);
				std::cout<<"HomogeneityFilter Started for Part 3  ..."<<std::endl;
				HomogeneityFilter(img3Doi_3, img3Dnr_3);
				std::cout<<"HomogeneityFilter Started for Part 4  ..."<<std::endl;
				HomogeneityFilter(img3Doi_4, img3Dnr_4);
			}
			if(redo_homogen_run == 2)
			{
				std::cout<<"Fixing Nonrigid Result 1  ..."<<std::endl;
				NonrigidFix(img3Dnr_1);
				std::cout<<"Fixing Nonrigid Result 2  ..."<<std::endl;
				NonrigidFix(img3Dnr_2);
				std::cout<<"Fixing Nonrigid Result 3  ..."<<std::endl;
				NonrigidFix(img3Dnr_3);
				std::cout<<"Fixing Nonrigid Result 4  ..."<<std::endl;
				NonrigidFix(img3Dnr_4);
			}
			
			img3Dtemp = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), sizez);
			img3Dtemp->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dtemp->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dtemp->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());

			for (int x = 0; x<img3Dnr_1->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_1->Num_Rows(); y++)
					for (int z = 0; z<zlevel; z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_1->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Dnr_2->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_2->Num_Rows(); y++)
					for (int z = zlevel; z<(2*zlevel); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_2->FastGet(x,y,z - zlevel));
					}
			for (int x = 0; x<img3Dnr_3->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_3->Num_Rows(); y++)
					for (int z = (2*zlevel); z<(3*zlevel); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_3->FastGet(x,y,z-(2*zlevel)));
					}
			for (int x = 0; x<img3Dnr_4->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_4->Num_Rows(); y++)
					for (int z = (3*zlevel); z<img3Dnr->Num_Levels(); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_4->FastGet(x,y,z-(3*zlevel)));
					}

			Write_Analyze_File("AfterHomogeneity.hdr",*img3Dtemp);
			
			delete img3Dnr_1;
			delete img3Dnr_2;
			delete img3Dnr_3;
			delete img3Dnr_4;
			delete img3Doi_1;
			delete img3Doi_2;
			delete img3Doi_3;
			delete img3Doi_4;

			MaxComp(img3Dtemp, img3Dtemp);
			HoleFill(img3Dtemp);

			img3Dtemp->Set_Pixel_SizeX(img3Dnr->Get_Pixel_SizeX());
			img3Dtemp->Set_Pixel_SizeY(img3Dnr->Get_Pixel_SizeY());
			img3Dtemp->Set_Pixel_SizeZ(img3Dnr->Get_Pixel_SizeZ());

			if(redo_homogen_run != 2)
			{
				hom_volume_step1 = GetVolume(img3Dtemp);
				volume_diff = abs(hom_volume_step1 - nonrigid_volume);
				volume_percent = 100 * (volume_diff / nonrigid_volume);

				for (int x = 0; x<img3Dtemp->Num_Cols(); x++)
					for (int y = 0; y<img3Dtemp->Num_Rows(); y++)
						for (int z = 0; z<img3Dtemp->Num_Levels(); z++)
						{
							if (img3Dtemp->FastGet(x,y,z) > 0)
							{
								img3Dtemp->FastSet(x,y,z,1000);
							}
						}
			}

			if (redo_homogen_run == 2)
			{
				img3Dtemp->Set_Pixel_SizeX(img3Dnr->Get_Pixel_SizeX());
				img3Dtemp->Set_Pixel_SizeY(img3Dnr->Get_Pixel_SizeY());
				img3Dtemp->Set_Pixel_SizeZ(img3Dnr->Get_Pixel_SizeZ());
				Write_Analyze_File(argv[3],*img3Dtemp);
				redo_homogen_run = 0;
			}
			else if ((volume_percent < Thresh_rerun) && redo_homogen_run != 2)
			{

				std::cout<<"First Homogeneity Worked!!!!"<<std::endl;
				std::cout<<"Nonrig vol: "<<nonrigid_volume<<"cm^3"<<std::endl;
				std::cout<<"Homogen vol: "<<hom_volume_step1<<"cm^3"<<std::endl;
				std::cout<<"Difference: "<<volume_percent<<"%"<<std::endl;
				img3Dtemp->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
				img3Dtemp->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
				img3Dtemp->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
				Write_Analyze_File(argv[3],*img3Dtemp);
				redo_homogen_run = 0;
			}
			else
			{

				std::cout<<"First Homogeneity Didn't Work :( . . . Using NonRigid Result"<<std::endl;
				std::cout<<"Nonrig vol: "<<nonrigid_volume<<"cm^3"<<std::endl;
				std::cout<<"Homogen vol: "<<hom_volume_step1<<"cm^3"<<std::endl;
				std::cout<<"Difference: "<<volume_percent<<"%"<<std::endl;
				redo_homogen_run = 2;
			}
		}
	}
	else if (pixelz >= 3 && pixelz < 5)
	{
		chunks = 3;
		int zlevel = sizez / chunks;
		int zremain = 2*zlevel;
		int zlevel_last = sizez - zremain;
		img3Doi_1 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
		img3Doi_2 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
		img3Doi_3 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel_last);
		img3Doi_1->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_1->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_1->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		img3Doi_2->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_2->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_2->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		img3Doi_3->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_3->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_3->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = 0; z<zlevel; z++)
				{
					img3Doi_1->FastSet(x,y,z,img3Doi->FastGet(x,y,z));
				}
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = zlevel; z<(2*zlevel); z++)
				{
					img3Doi_2->FastSet(x,y,z-zlevel,img3Doi->FastGet(x,y,z));
				}
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = (2*zlevel); z<img3Doi->Num_Levels(); z++)
				{
					img3Doi_3->FastSet(x,y,z-(2*zlevel),img3Doi->FastGet(x,y,z));
				}
		std::cout<<"Smoothing Started 1 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_1, img3Doi_1);
		std::cout<<"Smoothing Started 2 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_2, img3Doi_2);
		std::cout<<"Smoothing Started 3 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_3, img3Doi_3);


		while(redo_homogen_run != 0)
		{

			if (redo_homogen_run == 2)
			{
				delete img3Dtemp;
			}

			img3Dnr_1 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
			img3Dnr_2 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
			img3Dnr_3 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel_last);

			img3Dnr_1->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_1->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_1->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
			img3Dnr_2->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_2->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_2->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
			img3Dnr_3->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_3->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_3->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());

			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = 0; z<zlevel; z++)
					{
						img3Dnr_1->FastSet(x,y,z,img3Dnr->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = zlevel; z<(2*zlevel); z++)
					{
						img3Dnr_2->FastSet(x,y,z-zlevel,img3Dnr->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = (2*zlevel); z<img3Doi->Num_Levels(); z++)
					{
						img3Dnr_3->FastSet(x,y,z-(2*zlevel),img3Dnr->FastGet(x,y,z));
					}

	//		if (redo_homogen_run != 2)
	//		{
		//		std::cout<<"Attempting to get rid of the heart. . ."<<std::endl;
		//		CutHeart(img3Dnr_1,img3Dnr_1);
	//		}

			img3Dnr_1->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_1->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_1->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());


			if(redo_homogen_run != 2)
			{
				std::cout<<"HomogeneityFilter Started for Part 1  ..."<<std::endl;
				HomogeneityFilter(img3Doi_1, img3Dnr_1);
				std::cout<<"HomogeneityFilter Started for Part 2  ..."<<std::endl;
				HomogeneityFilter(img3Doi_2, img3Dnr_2);
				std::cout<<"HomogeneityFilter Started for Part 3  ..."<<std::endl;
				HomogeneityFilter(img3Doi_3, img3Dnr_3);
			}
			if(redo_homogen_run == 2)
			{
				std::cout<<"Fixing Nonrigid Result 1  ..."<<std::endl;
				NonrigidFix(img3Dnr_1);
				std::cout<<"Fixing Nonrigid Result 2  ..."<<std::endl;
				NonrigidFix(img3Dnr_2);
				std::cout<<"Fixing Nonrigid Result 3  ..."<<std::endl;
				NonrigidFix(img3Dnr_3);
			}
			
			img3Dtemp = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), sizez);
			img3Dtemp->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dtemp->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dtemp->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());

			for (int x = 0; x<img3Dnr_1->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_1->Num_Rows(); y++)
					for (int z = 0; z<zlevel; z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_1->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Dnr_2->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_2->Num_Rows(); y++)
					for (int z = zlevel; z<(2*zlevel); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_2->FastGet(x,y,z - zlevel));
					}
			for (int x = 0; x<img3Dnr_3->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_3->Num_Rows(); y++)
					for (int z = (2*zlevel); z<img3Dnr->Num_Levels(); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_3->FastGet(x,y,z-(2*zlevel)));
					}

			Write_Analyze_File("AfterHomogeneity.hdr",*img3Dtemp);

			delete img3Dnr_1;
			delete img3Dnr_2;
			delete img3Dnr_3;
			delete img3Doi_1;
			delete img3Doi_2;
			delete img3Doi_3;

			MaxComp(img3Dtemp, img3Dtemp);
			HoleFill(img3Dtemp);

			img3Dtemp->Set_Pixel_SizeX(img3Dnr->Get_Pixel_SizeX());
			img3Dtemp->Set_Pixel_SizeY(img3Dnr->Get_Pixel_SizeY());
			img3Dtemp->Set_Pixel_SizeZ(img3Dnr->Get_Pixel_SizeZ());

			if(redo_homogen_run != 2)
			{
				hom_volume_step1 = GetVolume(img3Dtemp);
				volume_diff = abs(hom_volume_step1 - nonrigid_volume);
				volume_percent = 100 * (volume_diff / nonrigid_volume);

				for (int x = 0; x<img3Dtemp->Num_Cols(); x++)
					for (int y = 0; y<img3Dtemp->Num_Rows(); y++)
						for (int z = 0; z<img3Dtemp->Num_Levels(); z++)
						{
							if (img3Dtemp->FastGet(x,y,z) > 0)
							{
								img3Dtemp->FastSet(x,y,z,1000);
							}
						}
			}


			if (redo_homogen_run == 2)
			{
				img3Dtemp->Set_Pixel_SizeX(img3Dnr->Get_Pixel_SizeX());
				img3Dtemp->Set_Pixel_SizeY(img3Dnr->Get_Pixel_SizeY());
				img3Dtemp->Set_Pixel_SizeZ(img3Dnr->Get_Pixel_SizeZ());
				Write_Analyze_File(argv[3],*img3Dtemp);
				redo_homogen_run = 0;
			}
			else if ((volume_percent < Thresh_rerun) && redo_homogen_run != 2)
			{

				std::cout<<"First Homogeneity Worked!!!!"<<std::endl;
				std::cout<<"Nonrig vol: "<<nonrigid_volume<<"cm^3"<<std::endl;
				std::cout<<"Homogen vol: "<<hom_volume_step1<<"cm^3"<<std::endl;
				std::cout<<"Difference: "<<volume_percent<<"%"<<std::endl;
				img3Dtemp->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
				img3Dtemp->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
				img3Dtemp->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
				Write_Analyze_File(argv[3],*img3Dtemp);
				redo_homogen_run = 0;
			}
			else
			{

				std::cout<<"First Homogeneity Didn't Work :( . . . Using NonRigid Result"<<std::endl;
				std::cout<<"Nonrig vol: "<<nonrigid_volume<<"cm^3"<<std::endl;
				std::cout<<"Homogen vol: "<<hom_volume_step1<<"cm^3"<<std::endl;
				std::cout<<"Difference: "<<volume_percent<<"%"<<std::endl;
				redo_homogen_run = 2;
			}


		}
	}
	else if (pixelz >= 5)
	{
		chunks = 2;
		int zlevel = sizez / chunks;
		int zremain = 1*zlevel;
		int zlevel_last = sizez - zremain;
		img3Doi_1 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
		img3Doi_2 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel_last);
		img3Doi_1->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_1->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_1->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		img3Doi_2->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
		img3Doi_2->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
		img3Doi_2->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
		
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = 0; z<zlevel; z++)
				{
					img3Doi_1->FastSet(x,y,z,img3Doi->FastGet(x,y,z));
				}
		for (int x = 0; x<img3Doi->Num_Cols(); x++)
			for (int y = 0; y<img3Doi->Num_Rows(); y++)
				for (int z = zlevel; z<img3Doi->Num_Levels(); z++)
				{
					img3Doi_2->FastSet(x,y,z-zlevel,img3Doi->FastGet(x,y,z));
				}
		std::cout<<"Smoothing Started 1 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_1, img3Doi_1);
		std::cout<<"Smoothing Started 2 ..."<<std::endl;
		AnisotropicDiffusion(img3Doi_2, img3Doi_2);

		while(redo_homogen_run != 0)
		{

			if (redo_homogen_run == 2)
			{
				delete img3Dtemp;
			}

			img3Dnr_1 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel);
			img3Dnr_2 = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), zlevel_last);
			img3Dnr_1->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_1->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_1->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
			img3Dnr_2->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_2->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_2->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());

			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = 0; z<zlevel; z++)
					{
						img3Dnr_1->FastSet(x,y,z,img3Dnr->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Doi->Num_Cols(); x++)
				for (int y = 0; y<img3Doi->Num_Rows(); y++)
					for (int z = zlevel; z<img3Doi->Num_Levels(); z++)
					{
						img3Dnr_2->FastSet(x,y,z-zlevel,img3Dnr->FastGet(x,y,z));
					}

	//		if (redo_homogen_run != 2)
	//		{
//				std::cout<<"Attempting to get rid of the heart. . ."<<std::endl;
	//			CutHeart(img3Dnr_1,img3Dnr_1);
	//		}

			img3Dnr_1->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dnr_1->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dnr_1->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());

			if(redo_homogen_run != 2)
			{
				std::cout<<"HomogeneityFilter Started for Part 1  ..."<<std::endl;
				HomogeneityFilter(img3Doi_1, img3Dnr_1);
				std::cout<<"HomogeneityFilter Started for Part 2  ..."<<std::endl;
				HomogeneityFilter(img3Doi_2, img3Dnr_2);
			}

			if(redo_homogen_run == 2)
			{
				std::cout<<"Fixing Nonrigid Result 1  ..."<<std::endl;
				NonrigidFix(img3Dnr_1);
				std::cout<<"Fixing Nonrigid Result 2  ..."<<std::endl;
				NonrigidFix(img3Dnr_2);
			}
			
			img3Dtemp = new CIS_Array_Image3D_short(img3Doi->Get_SizeX(), img3Doi->Get_SizeY(), sizez);
			img3Dtemp->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
			img3Dtemp->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
			img3Dtemp->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());

			for (int x = 0; x<img3Dnr_1->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_1->Num_Rows(); y++)
					for (int z = 0; z<zlevel; z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_1->FastGet(x,y,z));
					}
			for (int x = 0; x<img3Dnr_2->Num_Cols(); x++)
				for (int y = 0; y<img3Dnr_2->Num_Rows(); y++)
					for (int z = zlevel; z<img3Dnr->Num_Levels(); z++)
					{
						img3Dtemp->FastSet(x,y,z,img3Dnr_2->FastGet(x,y,z - zlevel));
					}
			Write_Analyze_File("AfterHomogeneity.hdr",*img3Dtemp);

			delete img3Dnr_1;
			delete img3Dnr_2;
			delete img3Doi_1;
			delete img3Doi_2;

			MaxComp(img3Dtemp, img3Dtemp);
			HoleFill(img3Dtemp);

			img3Dtemp->Set_Pixel_SizeX(img3Dnr->Get_Pixel_SizeX());
			img3Dtemp->Set_Pixel_SizeY(img3Dnr->Get_Pixel_SizeY());
			img3Dtemp->Set_Pixel_SizeZ(img3Dnr->Get_Pixel_SizeZ());

			if(redo_homogen_run != 2)
			{
				hom_volume_step1 = GetVolume(img3Dtemp);
				volume_diff = abs(hom_volume_step1 - nonrigid_volume);
				volume_percent = 100 * (volume_diff / nonrigid_volume);

				for (int x = 0; x<img3Dtemp->Num_Cols(); x++)
					for (int y = 0; y<img3Dtemp->Num_Rows(); y++)
						for (int z = 0; z<img3Dtemp->Num_Levels(); z++)
						{
							if (img3Dtemp->FastGet(x,y,z) > 0)
							{
								img3Dtemp->FastSet(x,y,z,1000);
							}
						}
			}

			if (redo_homogen_run == 2)
			{
				img3Dtemp->Set_Pixel_SizeX(img3Dnr->Get_Pixel_SizeX());
				img3Dtemp->Set_Pixel_SizeY(img3Dnr->Get_Pixel_SizeY());
				img3Dtemp->Set_Pixel_SizeZ(img3Dnr->Get_Pixel_SizeZ());
				Write_Analyze_File(argv[3],*img3Dtemp);
				redo_homogen_run = 0;
			}
			else if ((volume_percent < Thresh_rerun) && redo_homogen_run != 2)
			{

				std::cout<<"First Homogeneity Worked!!!!"<<std::endl;
				std::cout<<"Nonrig vol: "<<nonrigid_volume<<"cm^3"<<std::endl;
				std::cout<<"Homogen vol: "<<hom_volume_step1<<"cm^3"<<std::endl;
				std::cout<<"Difference: "<<volume_percent<<"%"<<std::endl;
				img3Dtemp->Set_Pixel_SizeX(img3Doi->Get_Pixel_SizeX());
				img3Dtemp->Set_Pixel_SizeY(img3Doi->Get_Pixel_SizeY());
				img3Dtemp->Set_Pixel_SizeZ(img3Doi->Get_Pixel_SizeZ());
				Write_Analyze_File(argv[3],*img3Dtemp);
				redo_homogen_run = 0;
			}
			else
			{

				std::cout<<"First Homogeneity Didn't Work :( . . . Using NonRigid Result"<<std::endl;
				std::cout<<"Nonrig vol: "<<nonrigid_volume<<"cm^3"<<std::endl;
				std::cout<<"Homogen vol: "<<hom_volume_step1<<"cm^3"<<std::endl;
				std::cout<<"Difference: "<<volume_percent<<"%"<<std::endl;
				redo_homogen_run = 2;
			}


		}
	}
	else
	{
		std::cout<<"Error with z direction pixel size: Unable to read"<<std::endl;
	}
}

double GetVolume(CIS_Array_Image3D_short *img3D1)
{
	int sizex=img3D1->Num_Cols();
	int sizey=img3D1->Num_Rows();
	int sizez=img3D1->Num_Levels();
	double result_volume = 0.0;
	double xpixel_conversion = img3D1->Get_Pixel_SizeX(); //mm
	double ypixel_conversion = img3D1->Get_Pixel_SizeY(); //mm
	double zpixel = img3D1->Get_Pixel_SizeZ(); //mm
	double xypixel_conversion = xpixel_conversion * ypixel_conversion; //area of one pixel

	for (int z=0; z<sizez; z++)
		for (int y=0; y<sizey; y++)
			for (int x=0; x<sizex; x++)
				if (img3D1->FastGet(x,y,z)>0)
				{
					result_volume += xypixel_conversion * zpixel;
				}
	result_volume = result_volume / 1000.0; // convert mm^3 to cm^3

	return result_volume;
}

void AnisotropicDiffusion(CIS_Array_Image3D_short *img3Din, CIS_Array_Image3D_short *img3Dout)
{
	int para_diffusionIteration = 5;
	double para_diffusionTimeStep = 0.0625;
		
	typedef itk::Image< float, 3 >   ImageType_3D;
	ImageType_3D::Pointer ITKim_3D = ImageType_3D::New();
	Copy_CISImage_to_ITKImage(img3Din, ITKim_3D);

	typedef itk::CurvatureAnisotropicDiffusionImageFilter< ImageType_3D,ImageType_3D> DiffusionFilterType;

	DiffusionFilterType::Pointer diffusion = DiffusionFilterType::New();

	diffusion->SetInput( ITKim_3D);
	diffusion->SetNumberOfIterations( para_diffusionIteration );
	diffusion->SetTimeStep(para_diffusionTimeStep);

	diffusion->Update();

	Copy_ITKImage_to_CISImage(diffusion->GetOutput(), img3Dout); 

	std::cout << "Smoothing Complete " << std::endl;
}

void HomogeneityFilter(CIS_Array_Image3D_short *img3Din1, CIS_Array_Image3D_short *img3Din2)
{
	typedef float PixelType;
	typedef unsigned char BinaryPixelType;
	typedef itk::Image< PixelType, 3 >   ImageType;
	ImageType::Pointer ITKim_3D1 = ImageType::New();
	ImageType::Pointer ITKim_3D2 = ImageType::New();

	Copy_CISImage_to_ITKImage(img3Din1, ITKim_3D1);
	Copy_CISImage_to_ITKImage(img3Din2, ITKim_3D2);

	typedef itk::ImageRegionIterator< ImageType > IteratorType;

	typedef itk::BinaryThresholdImageFilter< ImageType, ImageType > ThresholdingFilterType;
	ThresholdingFilterType::Pointer thresholder2 = ThresholdingFilterType::New();

	thresholder2->SetLowerThreshold( 1 );
	thresholder2->SetUpperThreshold( 101 );
	thresholder2->SetOutsideValue( 0 );
	thresholder2->SetInsideValue( 1 );
	
	thresholder2->SetInput( ITKim_3D2 );
	thresholder2->Update();

	Copy_ITKImage_to_CISImage(thresholder2->GetOutput(), img3Din2);
	Write_Analyze_File("FirstThresh.hdr", *img3Din2); 

	typedef itk::MultiplyImageFilter< ImageType, ImageType, ImageType > MultiplyFilterType;
	MultiplyFilterType::Pointer multiplier1 = MultiplyFilterType::New();
	MultiplyFilterType::Pointer multiplier2 = MultiplyFilterType::New();
	MultiplyFilterType::Pointer multiplier3 = MultiplyFilterType::New();


	int sz_x, sz_y, sz_z;
	sz_x = img3Din1->Get_SizeX();
	sz_y = img3Din1->Get_SizeY();
	sz_z = img3Din1->Get_SizeZ();

	float res_x, res_y, res_z;
	res_x = img3Din1->Get_Pixel_SizeX();
	res_y = img3Din1->Get_Pixel_SizeY();
	res_z = img3Din1->Get_Pixel_SizeZ();

////////////////////////////////////////////////
	// erode initial mask to keep just the center of liver (minimize errors after registration) to compute the min/max intensities for the convolution

	typedef itk::BinaryBallStructuringElement <BinaryPixelType, 3> StructuringElementType;
	typedef itk::GrayscaleErodeImageFilter <ImageType, ImageType, StructuringElementType> ErodeFilterType;
	typedef itk::GrayscaleDilateImageFilter <ImageType, ImageType, StructuringElementType> DilateFilterType;

	StructuringElementType structuringElement1;
	StructuringElementType structuringElement2;

	structuringElement1.SetRadius(1);  
	structuringElement1.CreateStructuringElement();
	
	structuringElement2.SetRadius(2); 
	structuringElement2.CreateStructuringElement();

	multiplier2->SetInput1( ITKim_3D1); ;
//	multiplier2->SetInput2( erode1->GetOutput() );
	multiplier2->SetInput2( thresholder2->GetOutput() );
	multiplier2->Update();

	Copy_ITKImage_to_CISImage(multiplier2->GetOutput(), img3Din2);
//	Write_Analyze_File("ErodedMultiplied.hdr", *img3Din2); 

	ThresholdingFilterType::Pointer thresholder4 = ThresholdingFilterType::New();

	thresholder4->SetLowerThreshold( 1025 ); // the liver in arterial or venous phases should not be outside of the values 1070-1190
	thresholder4->SetUpperThreshold( 1200 ); // these values are conservative
//	thresholder4->SetLowerThreshold( 1080 ); // the liver in arterial or venous phases should not be outside of the values 1070-1190
//	thresholder4->SetUpperThreshold( 1350 ); // these values are conservative
	thresholder4->SetOutsideValue( 0 );
	thresholder4->SetInsideValue( 1 );
	
	thresholder4->SetInput( multiplier2->GetOutput() );
	thresholder4->Update();

	Copy_ITKImage_to_CISImage(thresholder4->GetOutput(), img3Din2);
	Write_Analyze_File("ErodedThresh.hdr", *img3Din2); 

//	std::cout << "Erosion Done";

	// the result at this point is an eroded thersholded liver for which a histogram is computed to reject ouliers

//	Read_Analyze_File("ErodedThresh.hdr", *img3Din2);
//	Copy_CISImage_to_ITKImage(img3Din2, ITKim_3D2);

	multiplier3->SetInput1( ITKim_3D1); ;
	multiplier3->SetInput2( thresholder4->GetOutput() );
	//multiplier3->SetInput2( ITKim_3D2 );
	multiplier3->Update();

	ImageType::SizeType sizeB;
	sizeB[0] = sz_x;
	sizeB[1] = sz_y;
	sizeB[2] = sz_z;

	ImageType::IndexType startB;
	startB[0] = 0;
	startB[1] = 0;
	startB[2] = 0;

	ImageType::RegionType regionB;
	regionB.SetSize(sizeB);
	regionB.SetIndex(startB);

	int* Hist = new int[ 1201 ]; // max intensity value of 1200
//	int* Hist = new int[ 1351 ]; // max intensity value of 1200
	int tmp, i,  j, k;

	for ( i=0;  i<= 1200 ; i++ )
//	for ( i=0;  i<= 1350 ; i++ )
	{  
		Hist[i] = 0; 
	}
	
	// do the histogram
	float minB = 1025;
	float maxB = 1200;
//	float minB = 1080;
//	float maxB = 1350;
	int no = 0;
	//float sumB = 0;
	//int cptB = 0;
	IteratorType  outB( multiplier3->GetOutput() , regionB );
	for ( outB.GoToBegin(); !outB.IsAtEnd(); ++outB )
	{	
		if (outB.Get()!=0)
		{
			tmp = outB.Get();
			Hist[tmp] = Hist[tmp] + 1;
			no++;
	//		sumB = sumB + outB.Get();
	//		cptB = cptB + 1;
		}
	}
	//meanB = sumB/cptB;
	//MeanRes[0] = meanB;

	// compute the min and max of values over between lowest 15%  and under highest 10% intensities values
	float sumE = 0;
	float cptE = 0;
	float ratio = 0;
	float lim1 = 0.15; //15%
	float lim2 = 0.95; // highest 5%
	int c=0;

	for ( j=1025;  j<= 1200; j++ )
//	for ( j=1080;  j<= 1350; j++ )
	{ 
		if (ratio < lim1)
		{
			cptE = cptE + Hist[j];
		//	sumE = sumE + i*Hist[i];
			ratio = cptE/no;
		}
		else
		{
			break;
		}
	}

	/*for ( j=1080;  j<= 1200; j++ )
	{ 
		cptE =  Hist[j];
		//	sumE = sumE + i*Hist[i];
		ratio = cptE/no;
		if (ratio < lim1)
		{	break;	}
	}*/

	minB = j;

	cptE = 0;
	ratio = 0;
	/*for ( k=1200;  k>= 1060; k-- )
	{ 
		cptE = Hist[k];
		//	sumE = sumE + i*Hist[i];
			ratio = cptE/no;
		if (ratio < lim1)
		{
		break;
		}
	}*/

	for ( k=1025;  k<= 1200; k++ )
//	for ( k=1080;  k<= 1350; k++ )
	{ 
		if (ratio < lim2)
		{
			cptE = cptE + Hist[k];
		//	sumE = sumE + i*Hist[i];
			ratio = cptE/no;
					}
		else
		{
			break;
		}
	}

	maxB = k;
//	minB=k-50;
/*	if ((maxB > 1170)|(maxB < 1130)) // fully enhanced or arterial
	{ 
		minB = maxB - 35;
	}
	else if (maxB > 1150) // between arterial and venous
	{ 
		minB = maxB - 50;
	}
	else // between arterial and venous
	{
		maxB = maxB+5;
		minB= maxB - 50;
	}
*/

	std::cout << "Min Intensity Liver:";
	std::cout << minB << std::endl;
	std::cout << "Max Intensity Liver:";
	std::cout << maxB << std::endl;

//	MeanRes[1] = sumE/cptE;

/*
////////////////////////////////////////////////
	// use a gesodesic active contour to enlarge the inintial liver mask to ensure that no parts were missed after registration

	//Initiate Sigmoid Filter
	typedef itk::SigmoidImageFilter<ImageType, ImageType> SigmoidFilterType;
	SigmoidFilterType::Pointer sigmoidFilter=SigmoidFilterType::New();
	sigmoidFilter->SetOutputMinimum(0.0);
	sigmoidFilter->SetOutputMaximum(1.0);
	sigmoidFilter->SetAlpha(-0.25); //
	sigmoidFilter->SetBeta (9); //8
	
	//Initiate GeodesicActiveContourFilter
	typedef itk::GeodesicActiveContourLevelSetImageFilter<ImageType, ImageType> GeodesicActiveContourFilterType;
//	typedef itk::GeodesicActiveContourShapePriorLevelSetImageFilter<ImageType, ImageType>GeodesicActiveContourFilterType;
	GeodesicActiveContourFilterType::Pointer geodesicActiveContour=GeodesicActiveContourFilterType::New();

	geodesicActiveContour->SetPropagationScaling(0.5);
//	geodesicActiveContour->SetShapePriorScaling( 10.0 );
	geodesicActiveContour->SetCurvatureScaling(0.1);
	geodesicActiveContour->SetAdvectionScaling(0.5);
	geodesicActiveContour->SetMaximumRMSError( 0.01 );
	geodesicActiveContour->SetNumberOfIterations( 50 );
	//geodesicActiveContour->SetNumberOfLayers( 10 );

	ThresholdingFilterType::Pointer thresholder1 = ThresholdingFilterType::New();

	thresholder1->SetLowerThreshold( 0 );
	thresholder1->SetUpperThreshold( 5 );
	thresholder1->SetOutsideValue( 0 );
	thresholder1->SetInsideValue( 1 );

	//Initiate gradient magnitude Filter
	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter <ImageType, ImageType> FilterType;
	FilterType::Pointer gradientFilter= FilterType::New();

	gradientFilter->SetInput(ITKim_3D1);
	gradientFilter->SetSigma(0.7);
	gradientFilter->Update();
	Copy_ITKImage_to_CISImage(gradientFilter->GetOutput(), img3Din2);

//	Write_Analyze_File("Gradient.hdr", *img3Din2);
	std::cout<<"Gradient Filter Done!"<<std::endl;
	sigmoidFilter->SetInput(gradientFilter->GetOutput());
	sigmoidFilter->Update();

	Copy_ITKImage_to_CISImage(sigmoidFilter->GetOutput(), img3Din2);
//	Write_Analyze_File("Sigmoid.hdr", *img3Din2);
	std::cout<<"Sigmoid Filter Done!"<<std::endl;
	
	//GeodesicActiveContourPipeline
	geodesicActiveContour->SetInput(thresholder2->GetOutput()); 
	geodesicActiveContour->SetFeatureImage(sigmoidFilter->GetOutput()); 
	geodesicActiveContour->Update();

	Copy_ITKImage_to_CISImage(geodesicActiveContour->GetOutput(), img3Din2);
//	Write_Analyze_File("FirstGeodesic.hdr", *img3Din2);
	std::cout<<"Created First Geodesic!"<<std::endl;

	thresholder1->SetInput( geodesicActiveContour->GetOutput() );
	thresholder1->Update();
	
	Copy_ITKImage_to_CISImage(thresholder1->GetOutput(), img3Din2);
//	Write_Analyze_File("BigGAC.hdr", *img3Din2); 
	std::cout << "GAC Complete " << std::endl; 
*/
	multiplier1->SetInput1( ITKim_3D1); ;
//	multiplier1->SetInput2( thresholder1->GetOutput() ); //with GAC
	multiplier1->SetInput2( thresholder2->GetOutput() ); //without GAC
	multiplier1->Update();

	Copy_ITKImage_to_CISImage(multiplier1->GetOutput(), img3Din2);
//	Write_Analyze_File("Multiplied.hdr", *img3Din2); 
	Copy_CISImage_to_ITKImage(img3Din2, ITKim_3D2);


///////////////////////////////////////////////////////////////////

// convolution of 4D data

	typedef itk::ConstShapedNeighborhoodIterator< ImageType > ShapedNeighborhoodIteratorType;
	typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< ImageType > FaceCalculatorType;
	
//	bool connect26 = true;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Part I - Liver Segmentation

	ImageType::Pointer output1 = ImageType::New();
	
    Copy_CISImage_to_ITKImage(img3Din2, output1);

	// set the radius of the kernel
	ShapedNeighborhoodIteratorType::RadiusType radiusIt;
	radiusIt.Fill( 2 );
	
	FaceCalculatorType faceCalculator;
	FaceCalculatorType::FaceListType faceList;
	FaceCalculatorType::FaceListType::iterator fit;

	faceList = faceCalculator( ITKim_3D2, output1->GetRequestedRegion(), radiusIt );
	
	IteratorType out1;
  
	const float rad = static_cast<float>(2); 

	for ( fit=faceList.begin(); fit != faceList.end(); ++fit)
    {
		ShapedNeighborhoodIteratorType it( radiusIt, ITKim_3D2, *fit );
		out1 = IteratorType( output1, *fit );

		// Creates a circular structuring element by activating all the pixels less
		// than radius distance from the center of the neighborhood.
		
		for (float zz = -rad; zz <= rad; zz++)		
		{
			for (float yy = -rad; yy <= rad; yy++)
			{
				for (float xx = -rad; xx <= rad; xx++)
				{
					ShapedNeighborhoodIteratorType::OffsetType off;
					
					float dis = ::sqrt( xx*xx + yy*yy + zz*zz);
					
					if (dis <= rad)
					{   
						off[0] = static_cast<int>(xx);
						off[1] = static_cast<int>(yy);
						off[2] = static_cast<int>(zz);
						it.ActivateOffset(off);
					}
				}
			}
		}
		
		// Convolution
		for (it.GoToBegin(), out1.GoToBegin(); !it.IsAtEnd(); ++it, ++out1)
		{
			ShapedNeighborhoodIteratorType::ConstIterator ci;

			//get sum, min, max;
			int minSE = 5000;
			int maxSE = 0;
			int sumSE = 0;
			int noElemSE = 0;
			int stdSE = 0;
			int meanSE = 0;

			for (ci = it.Begin(); ci != it.End(); ci++)
			{
				//get Mean, Max, Min of the image covered by stencil
				int temp = ci.Get();

				if (temp < minSE)
				{
					minSE = temp;
				}

				if (temp > maxSE)
				{
					maxSE = temp;
				}

				sumSE +=  temp;
				noElemSE++; 
			}

			if ( ( minSE < minB ) || ( maxSE > maxB ) ) // conservative values, normally betwwen 1120(1150) and 1190
			{
			//	flag = false;
				out1.Set( 0 );
			}
			// create a  mask
			else 
			{
				out1.Set (1);		
			}
		}
    }

	std::cout << "Liver Segmentation Convolution Done" << std::endl;
	Copy_ITKImage_to_CISImage( output1, img3Din2);
	Copy_CISImage_to_ITKImage( img3Din2, ITKim_3D2 );
 ///////Skipped Max Components until later
	/*
	DilateFilterType::Pointer dilate1 =DilateFilterType::New();
	dilate1->SetKernel(structuringElement2);
	dilate1->SetInput(thresholder3->GetOutput());
	dilate1->Update();

*/
	typedef itk::VotingBinaryIterativeHoleFillingImageFilter<ImageType > FillFilterType;
	FillFilterType::Pointer fillfilter = FillFilterType::New();

	ImageType::SizeType indexRadius;
	indexRadius[0] = 4; // radius along x
	indexRadius[1] = 4; // radius along y
	indexRadius[2] = 0; // radius along z

	fillfilter->SetRadius( indexRadius );
	fillfilter->SetBackgroundValue( 0 );
	fillfilter->SetForegroundValue( 1 );
	fillfilter->SetMajorityThreshold( 1 );
	fillfilter->SetMaximumNumberOfIterations( 50);

	fillfilter->SetInput(ITKim_3D2);
	fillfilter->Update();


// closing

/*
	DilateFilterType::Pointer dilate2 =DilateFilterType::New();
	dilate2->SetKernel(structuringElement1);
	dilate2->SetInput(fillfilter->GetOutput());
	dilate2->Update();

	
	ErodeFilterType::Pointer erode1 = ErodeFilterType::New();
	erode1->SetKernel(structuringElement1);
	erode1->SetInput(dilate2->GetOutput());
	erode1->Update();
*/
	std::cout << "Liver Segmentation Done" << std::endl;
	Copy_ITKImage_to_CISImage(fillfilter->GetOutput(), img3Din2);
	
	img3Din2->Set_Pixel_SizeX(img3Din1->Get_Pixel_SizeX());
	img3Din2->Set_Pixel_SizeY(img3Din1->Get_Pixel_SizeY());
	img3Din2->Set_Pixel_SizeZ(img3Din1->Get_Pixel_SizeZ());

	for (int x = 0; x<img3Din2->Num_Cols(); x++)
		for (int y = 0; y<img3Din2->Num_Rows(); y++)
			for (int z = 0; z<img3Din2->Num_Levels(); z++)
			{
				if (img3Din2->FastGet(x,y,z) > 0)
				{
					img3Din2->FastSet(x,y,z,1000);
				}
			}

}
void NonrigidFix(CIS_Array_Image3D_short *img3Din2)
{
	typedef float PixelType;
	typedef unsigned char BinaryPixelType;
	typedef itk::Image< PixelType, 3 >   ImageType;
	ImageType::Pointer ITKim_3D2 = ImageType::New();
	Copy_CISImage_to_ITKImage(img3Din2, ITKim_3D2);
	typedef itk::ImageRegionIterator< ImageType > IteratorType;
	typedef itk::BinaryThresholdImageFilter< ImageType, ImageType > ThresholdingFilterType;
	ThresholdingFilterType::Pointer thresholder2 = ThresholdingFilterType::New();
	thresholder2->SetLowerThreshold( 1 );
	thresholder2->SetUpperThreshold( 100 );
	thresholder2->SetOutsideValue( 0 );
	thresholder2->SetInsideValue( 1 );
	thresholder2->SetInput( ITKim_3D2 );
	thresholder2->Update();
////////////////////////////////////////////////
	// erode initial mask to get rid of sharp edges then Dilate

	typedef itk::BinaryBallStructuringElement <BinaryPixelType, 3> StructuringElementType;
	typedef itk::GrayscaleErodeImageFilter <ImageType, ImageType, StructuringElementType> ErodeFilterType;
	typedef itk::GrayscaleDilateImageFilter <ImageType, ImageType, StructuringElementType> DilateFilterType;

	StructuringElementType structuringElement1;
	StructuringElementType structuringElement2;

	structuringElement1.SetRadius(2);  
	structuringElement1.CreateStructuringElement();
	
	structuringElement2.SetRadius(2); 
	structuringElement2.CreateStructuringElement();

	typedef itk::VotingBinaryIterativeHoleFillingImageFilter<ImageType > FillFilterType;
	FillFilterType::Pointer fillfilter = FillFilterType::New();

	ImageType::SizeType indexRadius;
	indexRadius[0] = 4; // radius along x
	indexRadius[1] = 4; // radius along y
	indexRadius[2] = 0; // radius along y

	fillfilter->SetRadius( indexRadius );
	fillfilter->SetBackgroundValue( 0 );
	fillfilter->SetForegroundValue( 1 );
	fillfilter->SetMajorityThreshold( 1 );
	fillfilter->SetMaximumNumberOfIterations( 50);

	fillfilter->SetInput(thresholder2->GetOutput() );
	fillfilter->Update();

	ErodeFilterType::Pointer erode1 = ErodeFilterType::New();
	erode1->SetKernel(structuringElement2);
	erode1->SetInput(fillfilter->GetOutput());
	erode1->Update();

	DilateFilterType::Pointer dilate2 =DilateFilterType::New();
	dilate2->SetKernel(structuringElement1);
	dilate2->SetInput(erode1->GetOutput());
	dilate2->Update();

	std::cout << "Non Rigid Processing Done" << std::endl;
	Copy_ITKImage_to_CISImage( dilate2->GetOutput(), img3Din2);
	
}
void CutHeart(CIS_Array_Image3D_short *img3Din, CIS_Array_Image3D_short *img3Dout)
{
	Write_Analyze_File("BeforeHeartOut.hdr", *img3Din);
	typedef float PixelType;
	typedef unsigned char BinaryPixelType;
	typedef itk::Image< PixelType, 3 >   ImageType;
	ImageType::Pointer ITKim_3Din = ImageType::New();
	
	typedef itk::BinaryBallStructuringElement <BinaryPixelType, 3> StructuringElementType;
	typedef itk::GrayscaleErodeImageFilter <ImageType, ImageType, StructuringElementType> ErodeFilterType;
	typedef itk::GrayscaleDilateImageFilter <ImageType, ImageType, StructuringElementType> DilateFilterType;

	img3Dorig = new CIS_Array_Image3D_short(img3Din->Get_SizeX(), img3Din->Get_SizeY(), img3Din->Get_SizeZ());
	img3Dorig->Set_Pixel_SizeX(img3Din->Get_Pixel_SizeX());
	img3Dorig->Set_Pixel_SizeY(img3Din->Get_Pixel_SizeY());
	img3Dorig->Set_Pixel_SizeZ(img3Din->Get_Pixel_SizeZ());
	for (int x = 0; x<img3Din->Num_Cols(); x++)
		for (int y = 0; y<img3Din->Num_Rows(); y++)
			for (int z = 0; z<img3Din->Num_Levels(); z++)
				img3Dorig->FastSet(x,y,z,img3Din->FastGet(x,y,z));
	
	for (int x = 0; x<img3Din->Num_Cols(); x++)
		for (int y = 0; y<img3Din->Num_Rows(); y++)
			for (int z = 0; z<img3Din->Num_Levels(); z++)
			{
				if (img3Din->FastGet(x,y,z) > 0)
				{
					img3Din->FastSet(x,y,z,1);
				}
			}
	int padval = 10;
	std::cout<<"1st Step: Copying Original Image Complete . . ."<<std::endl;
///// Need to Pad the end of chunk so dialation can take effect
	img3Dpad = new CIS_Array_Image3D_short(img3Din->Get_SizeX(), img3Din->Get_SizeY(), img3Din->Num_Levels() + padval);
	img3Dpad->Set_Pixel_SizeX(img3Din->Get_Pixel_SizeX());
	img3Dpad->Set_Pixel_SizeY(img3Din->Get_Pixel_SizeY());
	img3Dpad->Set_Pixel_SizeZ(img3Din->Get_Pixel_SizeZ());
	
	for (int x = 0; x<img3Dpad->Num_Cols(); x++)
		for (int y = 0; y<img3Dpad->Num_Rows(); y++)
			for (int z = 0; z<img3Dpad->Num_Levels(); z++)
			{
				if (z < img3Din->Num_Levels())
				{
					img3Dpad->FastSet(x,y,z,img3Din->FastGet(x,y,z));
				}
				else
				{
					img3Dpad->FastSet(x,y,z,0);
				}
			}
	std::cout<<"1st Step: Padding Complete . . ."<<std::endl;
	Write_Analyze_File("img3Dpad.hdr", *img3Dpad);
	Copy_CISImage_to_ITKImage(img3Dpad, ITKim_3Din);

	StructuringElementType structuringElement1;
	structuringElement1.SetRadius(4);  
	structuringElement1.CreateStructuringElement();

	StructuringElementType structuringElement2;
	structuringElement2.SetRadius(6);  
	structuringElement2.CreateStructuringElement();

	ErodeFilterType::Pointer erode1 = ErodeFilterType::New();
	erode1->SetKernel(structuringElement1);
	erode1->SetInput(ITKim_3Din);
	erode1->Update();
	std::cout<<"1st Step: Erosion Complete . . ."<<std::endl;
////////////////////////////////////////////////MAX COMP TO FIND HEART
	int sz_x, sz_y, sz_z;
	sz_x = img3Dpad->Get_SizeX();
	sz_y = img3Dpad->Get_SizeY();
	sz_z = img3Dpad->Get_SizeZ();

	float res_x, res_y, res_z;
	res_x = img3Dpad->Get_Pixel_SizeX();
	res_y = img3Dpad->Get_Pixel_SizeY();
	res_z = img3Dpad->Get_Pixel_SizeZ();
	//////////////////////////////////////////////////////////////
	//keep maximum component

	bool connect26 = true;

	float valueThresMax; //value inside the max connected component, see later;
	int minsz; 
	int maxsx;

	intDynArray blobRank;
	intVec3DynArray blobCentroid;
	
	CIS_Array_Image3D_short *bim=NULL;
	bim = new CIS_Array_Image3D_short(sz_x, sz_y, sz_z);
	CIS_Array_Image3D_short *blobImg=NULL;
	blobImg = new CIS_Array_Image3D_short(sz_x, sz_y, sz_z);

	Copy_ITKImage_to_CISImage(erode1->GetOutput(), bim);

	NIH_Algo_Blob_Labelling_3D(bim, blobImg, blobRank, blobCentroid, connect26);

	minsz = sz_z;
	maxsx = 0;
	valueThresMax = 0;

	// find Lowest Z Value blobCentroid
	for (int i=0; i<blobCentroid.GetSize(); i++)
	{
		if (blobCentroid[i].z< minsz && blobCentroid[i].x > maxsx)
		{
			minsz = blobCentroid[i].z;
			maxsx = blobCentroid[i].x;
			valueThresMax = i+1;
		}
	}
	
	Copy_CISImage_to_ITKImage(blobImg, ITKim_3Din);
	Write_Analyze_File("blobImg.hdr",*blobImg);

  	delete bim;
	delete blobImg;

	typedef itk::BinaryThresholdImageFilter< ImageType, ImageType > ThresholdingFilterType;
	ThresholdingFilterType::Pointer thresholder3 = ThresholdingFilterType::New();

	thresholder3->SetInput(ITKim_3Din);
	thresholder3->SetLowerThreshold(valueThresMax);
	thresholder3->SetUpperThreshold(valueThresMax); // Change this because of liver and spleen combo
	thresholder3->SetOutsideValue(0);
	thresholder3->SetInsideValue(1);
	thresholder3->Update();

	std::cout<<"1st Step: 'Heart' Localized (If Applicable)"<<std::endl;

	DilateFilterType::Pointer dilate3 =DilateFilterType::New();
	dilate3->SetKernel(structuringElement2);
	dilate3->SetInput(thresholder3->GetOutput());
	dilate3->Update();

	Copy_ITKImage_to_CISImage(dilate3->GetOutput(), img3Dpad);

	for (int x = 0; x<img3Din->Num_Cols(); x++)
		for (int y = 0; y<img3Din->Num_Rows(); y++)
			for (int z = 0; z<img3Din->Num_Levels(); z++)
			{
				if (img3Dpad->FastGet(x,y,z) > 0)
				{
					img3Din->FastSet(x,y,z,100);
				}
				else
				{
					img3Din->FastSet(x,y,z,0);
				}
			}
	std::cout<<"1st Step: 'Heart' Dialated . . ."<<std::endl;
	Write_Analyze_File("TheHeart_1stStep.hdr", *img3Din);
	
///// Determine if image has a problem segmenting the heart, if not just use original image
	int intersect=0;
	int union_img1=0;
	int union_img2=0;
	float overlap=0;
	for (int x = 0; x<img3Din->Num_Cols(); x++)
		for (int y = 0; y<img3Din->Num_Rows(); y++)
			for (int z = 0; z<img3Din->Num_Levels(); z++)
			{
				if (img3Din->FastGet(x,y,z)>0 && img3Dorig->FastGet(x,y,z)>0)
					intersect++;
				if (img3Din->FastGet(x,y,z)>0)
					union_img1++;
				if (img3Dorig->FastGet(x,y,z)>0)
					union_img2++;
			}
	overlap = 2.0*(float)intersect/((float)union_img1 + (float)union_img2)*100.0;
	std::cout<<"Overlap = " << overlap << "%\n";

	if (overlap >= 65.0)
	{
		std::cout<<"Heart not present in segmentation image!!!!"<<std::endl;
		std::cout<<"Returning original image with no modifications"<<std::endl;
		std::cout<<"Heart Cut Complete!!!!!"<<std::endl;
		for (int x = 0; x<img3Din->Num_Cols(); x++)
			for (int y = 0; y<img3Din->Num_Rows(); y++)
				for (int z = 0; z<img3Din->Num_Levels(); z++)
					img3Dout->FastSet(x,y,z,img3Dorig->FastGet(x,y,z));

		Write_Analyze_File("AfterHeartOut_in_before.hdr", *img3Din);
		for (int x = 0; x<img3Din->Num_Cols(); x++)
			for (int y = 0; y<img3Din->Num_Rows(); y++)
				for (int z = 0; z<img3Din->Num_Levels(); z++)
					img3Din->FastSet(x,y,z,img3Dorig->FastGet(x,y,z));

		Write_Analyze_File("AfterHeartOut_out.hdr", *img3Dout);
		Write_Analyze_File("AfterHeartOut_in_after.hdr", *img3Din);
		delete img3Dorig;
	}
	else
	{
		std::cout<<"Heart present in segmentation image :("<<std::endl;
		std::cout<<"Further Modifying . . . "<<std::endl;
		int pixel_diff = 0;

		for (int x = 0; x<img3Din->Num_Cols(); x++)
			for (int y = 0; y<img3Din->Num_Rows(); y++)
				for (int z = 0; z<img3Din->Num_Levels(); z++)
				{
					pixel_diff = img3Dorig->FastGet(x,y,z) - img3Din->FastGet(x,y,z);
					
					if (pixel_diff <= 0) 
					{
						img3Din->FastSet(x,y,z,0);
					}
					else
					{
						img3Din->FastSet(x,y,z,100);
					}
				}
		std::cout<<"1st Step: Heart Deleted From Original"<<std::endl;
		std::cout<<"1st Step: Complete"<<std::endl;
		Write_Analyze_File("LiverMinusHeart.hdr", *img3Din);

		////////////////////////////////////////////////////////////////////////////////////
		//SECOND STEP - GET RID OF RESIDUAL HEART PEICES FIRST STEP DIDN'T TAKE CARE OF
		////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"Beginning Step 2 . . ."<<std::endl;
		std::cout<<"2nd Step: Setting Image to Binary . . ."<<std::endl;
		for (int x = 0; x<img3Din->Num_Cols(); x++)
			for (int y = 0; y<img3Din->Num_Rows(); y++)
				for (int z = 0; z<img3Din->Num_Levels(); z++)
				{
					if (img3Din->FastGet(x,y,z) > 0)
					{
						img3Din->FastSet(x,y,z,1);
					}
				}
		Copy_CISImage_to_ITKImage(img3Din, ITKim_3Din);
		std::cout<<"2nd Step: Beginning Max Components . . ."<<std::endl;
		////////////////////////////////////////////////MAX COMP TO GET RID OF RESIDUALS
		int sz_x2, sz_y2, sz_z2;
		sz_x2 = img3Din->Get_SizeX();
		sz_y2 = img3Din->Get_SizeY();
		sz_z2 = img3Din->Get_SizeZ();

		float res_x2, res_y2, res_z2;
		res_x2 = img3Din->Get_Pixel_SizeX();
		res_y2 = img3Din->Get_Pixel_SizeY();
		res_z2 = img3Din->Get_Pixel_SizeZ();

		//////////////////////////////////////////////////////////////
		//keep maximum component - need help with this part

		connect26 = true;

		float valueThresMax2; //value inside the max connected component, see later;
		int maxsz; //  size of max connected component, see later;
		int minsx; //  size of min connected component
		int maxsy;

		intDynArray blobRank2;
		intVec3DynArray blobCentroid2;

		CIS_Array_Image3D_short *bim2=NULL;
		bim2 = new CIS_Array_Image3D_short(sz_x2, sz_y2, sz_z2);
		CIS_Array_Image3D_short *blobImg2=NULL;
		blobImg2 = new CIS_Array_Image3D_short(sz_x2, sz_y2, sz_z2);

		Copy_ITKImage_to_CISImage(ITKim_3Din, bim2);

		NIH_Algo_Blob_Labelling_3D(bim2, blobImg2, blobRank2, blobCentroid2, connect26);

		maxsz = 0;
		minsx = 511;
		maxsy = 0;
		valueThresMax2 = 0;
		// find Lowest Z Value blobCentroid
		for (int i=0; i<blobCentroid2.GetSize(); i++)
		{
			if (blobCentroid2[i].z> maxsz && blobCentroid2[i].x < minsx && blobCentroid2[i].y > maxsy)
			{
				maxsz = blobCentroid2[i].z;
				minsx = blobCentroid2[i].x;
				maxsy = blobCentroid2[i].y;
				valueThresMax2 = i+1;
			}
		}
		Copy_CISImage_to_ITKImage(blobImg2, ITKim_3Din);
		Write_Analyze_File("blobImg2.hdr",*blobImg2);

  		delete bim2;
		delete blobImg2;

		typedef itk::BinaryThresholdImageFilter< ImageType, ImageType > ThresholdingFilterType;
		ThresholdingFilterType::Pointer thresholder6 = ThresholdingFilterType::New();
		
		thresholder6->SetInput(ITKim_3Din);
		thresholder6->SetLowerThreshold(valueThresMax2);
		thresholder6->SetUpperThreshold(valueThresMax2); // keep only voxels in max connect component;
		thresholder6->SetOutsideValue(0);
		thresholder6->SetInsideValue(1);
		thresholder6->Update();

		Copy_ITKImage_to_CISImage(thresholder6->GetOutput(), img3Dout);
		std::cout<<"2nd Step: Max Components Complete - Liver Should Just Remain . . ."<<std::endl;

		for (int x = 0; x<img3Dout->Num_Cols(); x++)
			for (int y = 0; y<img3Dout->Num_Rows(); y++)
				for (int z = 0; z<img3Dout->Num_Levels(); z++)
				{
					if (img3Dout->FastGet(x,y,z) > 0)
					{
						img3Dout->FastSet(x,y,z,100);
					}
				}

		std::cout << "Heart Cut Complete: Hopefully Heart is out now" << std::endl;
		Write_Analyze_File("AfterHeartOut.hdr",*img3Dout);
		img3Din = img3Dout;
		delete img3Dpad;
		delete img3Dorig;
	}

}
void MaxComp(CIS_Array_Image3D_short *img3Din, CIS_Array_Image3D_short *img3Dout)
{
	for (int x = 0; x<img3Din->Num_Cols(); x++)
		for (int y = 0; y<img3Din->Num_Rows(); y++)
			for (int z = 0; z<img3Din->Num_Levels(); z++)
			{
				if (img3Din->FastGet(x,y,z) > 0)
				{
					img3Din->FastSet(x,y,z,1);
				}
			}
	int sz_x, sz_y, sz_z;
	sz_x = img3Din->Get_SizeX();
	sz_y = img3Din->Get_SizeY();
	sz_z = img3Din->Get_SizeZ();

	float res_x, res_y, res_z;
	res_x = img3Din->Get_Pixel_SizeX();
	res_y = img3Din->Get_Pixel_SizeY();
	res_z = img3Din->Get_Pixel_SizeZ();

	typedef float PixelType;
	typedef unsigned char BinaryPixelType;
	typedef itk::Image< PixelType, 3 >   ImageType;
	ImageType::Pointer ITKim_3Din = ImageType::New();

	Copy_CISImage_to_ITKImage(img3Din, ITKim_3Din);

	typedef itk::BinaryBallStructuringElement <BinaryPixelType, 3> StructuringElementType;
	typedef itk::GrayscaleErodeImageFilter <ImageType, ImageType, StructuringElementType> ErodeFilterType;
	typedef itk::GrayscaleDilateImageFilter <ImageType, ImageType, StructuringElementType> DilateFilterType;

	StructuringElementType structuringElement1;
	StructuringElementType structuringElement2;

	structuringElement1.SetRadius(1);  
	structuringElement1.CreateStructuringElement();
	
	structuringElement2.SetRadius(2); 
	structuringElement2.CreateStructuringElement();

//////////////////////////////////////////////////////////////
	//keep maximum component

	bool connect26 = true;

	float valueThresMax; //value inside the max connected component, see later;
	int maxsz; //  size of max connected component, see later;

	intDynArray blobRank;
	intVec3DynArray blobCentroid;
	
	CIS_Array_Image3D_short *bim=NULL;
	bim = new CIS_Array_Image3D_short(sz_x, sz_y, sz_z);
	CIS_Array_Image3D_short *blobImg=NULL;
	blobImg = new CIS_Array_Image3D_short(sz_x, sz_y, sz_z);

	Copy_ITKImage_to_CISImage(ITKim_3Din, bim);

	NIH_Algo_Blob_Labelling_3D(bim, blobImg, blobRank, blobCentroid, connect26);

	maxsz = 0;
	valueThresMax = 0;

	// find the maximum blobRank
	for (int i=0; i< blobRank.GetSize(); i++)
		{
		if (blobRank[i] > maxsz)
			{
			maxsz = blobRank[i];
			valueThresMax = i+1;
			}
		}

	//	Write_Analyze_File("Components2.hdr", *blobImg); 

		Copy_CISImage_to_ITKImage(blobImg, ITKim_3Din);

  	delete bim;
	delete blobImg;

	typedef itk::BinaryThresholdImageFilter< ImageType, ImageType > ThresholdingFilterType;
	ThresholdingFilterType::Pointer thresholder3 = ThresholdingFilterType::New();

	thresholder3->SetInput(ITKim_3Din);
	thresholder3->SetLowerThreshold(valueThresMax);
	thresholder3->SetUpperThreshold(valueThresMax); // keep only voxels in max connect component;
	thresholder3->SetOutsideValue(0);
	thresholder3->SetInsideValue(1);
	thresholder3->Update();
	
	DilateFilterType::Pointer dilate1 =DilateFilterType::New();
	dilate1->SetKernel(structuringElement2);
	dilate1->SetInput(thresholder3->GetOutput());
	dilate1->Update();

	std::cout << "Maximum Components Done !!!!" << std::endl;
	Copy_ITKImage_to_CISImage(dilate1->GetOutput(), img3Dout);
	for (int x = 0; x<img3Dout->Num_Cols(); x++)
		for (int y = 0; y<img3Dout->Num_Rows(); y++)
			for (int z = 0; z<img3Dout->Num_Levels(); z++)
			{
				if (img3Dout->FastGet(x,y,z) > 0)
				{
					img3Dout->FastSet(x,y,z,1000);
				}
			}
	img3Dout->Set_Pixel_SizeX(img3Din->Get_Pixel_SizeX());
	img3Dout->Set_Pixel_SizeY(img3Din->Get_Pixel_SizeY());
	img3Dout->Set_Pixel_SizeZ(img3Din->Get_Pixel_SizeZ());
	img3Din = img3Dout;
	Write_Analyze_File("AfterMaxComp.hdr",*img3Din);
}
void HoleFill(CIS_Array_Image3D_short *img3Din)
{
	std::cout<<"Starting HoleFill Step 1!!"<<std::endl;
	img3Dorig = new CIS_Array_Image3D_short(img3Din->Get_SizeX(), img3Din->Get_SizeY(), img3Din->Get_SizeZ());

	for (int x = 0; x<img3Din->Num_Cols(); x++)
		for (int y = 0; y<img3Din->Num_Rows(); y++)
			for (int z = 0; z<img3Din->Num_Levels(); z++)
			{
				if (img3Din->FastGet(x,y,z) > 0)
				{
					img3Dorig->FastSet(x,y,z,1);
					img3Din->FastSet(x,y,z,0);
				}
				else
				{
					img3Dorig->FastSet(x,y,z,0);
					img3Din->FastSet(x,y,z,1);
				}
			}

	typedef float PixelType;
	typedef unsigned char BinaryPixelType;
	typedef itk::Image< PixelType, 3 >   ImageType;
	ImageType::Pointer ITKim_3Din = ImageType::New();
	ImageType::Pointer ITKim_3Din2 = ImageType::New();

	Copy_CISImage_to_ITKImage(img3Din, ITKim_3Din);

	typedef itk::BinaryBallStructuringElement <BinaryPixelType, 3> StructuringElementType;
	typedef itk::GrayscaleErodeImageFilter <ImageType, ImageType, StructuringElementType> ErodeFilterType;
	typedef itk::GrayscaleDilateImageFilter <ImageType, ImageType, StructuringElementType> DilateFilterType;

	StructuringElementType structuringElement1;

	structuringElement1.SetRadius(1);  
	structuringElement1.CreateStructuringElement();

	StructuringElementType structuringElement2;

	structuringElement2.SetRadius(3);  
	structuringElement2.CreateStructuringElement();

	ErodeFilterType::Pointer erode1 = ErodeFilterType::New();
	erode1->SetKernel(structuringElement1);
	erode1->SetInput(ITKim_3Din);
	erode1->Update();

	Copy_ITKImage_to_CISImage(erode1->GetOutput(), img3Din);

	/*=======Flood Fill=======*/
	std::queue<ipixel> pixelq;
	ipixel pi0, piq, pie, piw;
	int pie2, piw2;
	int piy, piz;
	int replacement_color = 2;
	int target_color = 1;
	int pie_color, piw_color;
	pi0.x = 0; pi0.y = 0;pi0.z = 0; pi0.value = img3Din->FastGet(0,0,0);
	
	pixelq.push(pi0);
	while(pixelq.size() > 0){
		piq = pixelq.front();
		pixelq.pop();
		piq.value = img3Din->FastGet(piq.x, piq.y, piq.z);
		if(piq.value == target_color){
			pie = piq;
			piw = piq;
			pie_color = pie.value;
			piw_color = piw.value;

			piy = piq.y;
			piz = piq.z;
			pie2 = pie.x+1;
			while((pie_color == target_color)&&(pie2 < img3Din->Num_Cols())){
				pie_color = img3Din->FastGet(pie2, piy, piz);
				pie2++;
			}
			piw2 = piw.x-1;
			while((piw_color == target_color)&&(piw2 > 0)){
				piw_color = img3Din->FastGet(piw2, piy, piz);
				piw2--;
			}
			for(int i = pie.x; i < pie2; i++){
				img3Din->FastSet(i,piy, piz, replacement_color);
			}
			for(int i = piw.x; i > piw2; i--){
				img3Din->FastSet(i, piy, piz, replacement_color);
			}
			for(int i = piw2+1; i < pie2; i++){
				piq.x = i; piq.y = piy+1; piq.z = piz; 
				if(piq.y < img3Din->Num_Rows()){
					piq.value = img3Din->FastGet(piq.x, piq.y, piq.z);
					if(piq.value == target_color){
						pixelq.push(piq);
					}
				}
				piq.y=piy-1;
				if(piq.y > 0){
					piq.value = img3Din->FastGet(piq.x, piq.y,piq.z);
					if(piq.value == target_color){
						pixelq.push(piq);
					}
				}
				piq.y = piy;
				piq.z = piz -1;
				if(piq.z > 0){
					piq.value = img3Din->FastGet(piq.x, piq.y, piq.z);
					if(piq.value == target_color){
						pixelq.push(piq);
					}
				}
				piq.z = piz + 1;
				if(piq.z < img3Din->Num_Levels()){
					piq.value = img3Din->FastGet(piq.x, piq.y, piq.z);
					if(piq.value == target_color){
						pixelq.push(piq);
					}
				}

			}
		}
	}

	for (int x = 0; x<img3Din->Num_Cols(); x++)
		for (int y = 0; y<img3Din->Num_Rows(); y++)
			for (int z = 0; z<img3Din->Num_Levels(); z++)
			{
				if(img3Din->FastGet(x,y,z)>1)
				{
					img3Din->FastSet(x,y,z,0);
				}
			}

	Copy_CISImage_to_ITKImage(img3Din, ITKim_3Din2);

	DilateFilterType::Pointer dilate1 =DilateFilterType::New();
	dilate1->SetKernel(structuringElement2);
	dilate1->SetInput(ITKim_3Din2);
	dilate1->Update();

	Copy_ITKImage_to_CISImage(dilate1->GetOutput(), img3Din);

	for (int x = 0; x<img3Dorig->Num_Cols(); x++)
		for (int y = 0; y<img3Dorig->Num_Rows(); y++)
			for (int z = 0; z<img3Dorig->Num_Levels(); z++)
				img3Din->FastSet(x,y,z,img3Dorig->FastGet(x,y,z) + img3Din->FastGet(x,y,z));

	for (int x = 0; x<img3Din->Num_Cols(); x++)
		for (int y = 0; y<img3Din->Num_Rows(); y++)
			for (int z = 0; z<img3Din->Num_Levels(); z++)
			{
				if(img3Din->FastGet(x,y,z)>0)
				{
					img3Din->FastSet(x,y,z,1000);
				}
			}

	Write_Analyze_File("AfterFloodFill1.hdr",*img3Din);
	std::cout<<"HoleFill Step 1 Complete!!"<<std::endl;
	delete img3Dorig;

	std::cout<<"Starting HoleFill Step 2!!"<<std::endl;
	img3Dorig = new CIS_Array_Image3D_short(img3Din->Get_SizeX(), img3Din->Get_SizeY(), img3Din->Get_SizeZ());

	for (int x = 0; x<img3Din->Num_Cols(); x++)
		for (int y = 0; y<img3Din->Num_Rows(); y++)
			for (int z = 0; z<img3Din->Num_Levels(); z++)
			{
				if (img3Din->FastGet(x,y,z) > 0)
				{
					img3Dorig->FastSet(x,y,z,1);
					img3Din->FastSet(x,y,z,0);
				}
				else
				{
					img3Dorig->FastSet(x,y,z,0);
					img3Din->FastSet(x,y,z,1);
				}
			}

	/*=======Flood Fill=======*/
	pi0.x = 0; pi0.y = 0;pi0.z = 0; pi0.value = img3Din->FastGet(0,0,0);
	
	pixelq.push(pi0);
	while(pixelq.size() > 0){
		piq = pixelq.front();
		pixelq.pop();
		piq.value = img3Din->FastGet(piq.x, piq.y, piq.z);
		if(piq.value == target_color){
			pie = piq;
			piw = piq;
			pie_color = pie.value;
			piw_color = piw.value;

			piy = piq.y;
			piz = piq.z;
			pie2 = pie.x+1;
			while((pie_color == target_color)&&(pie2 < img3Din->Num_Cols())){
				pie_color = img3Din->FastGet(pie2, piy, piz);
				pie2++;
			}
			piw2 = piw.x-1;
			while((piw_color == target_color)&&(piw2 > 0)){
				piw_color = img3Din->FastGet(piw2, piy, piz);
				piw2--;
			}
			for(int i = pie.x; i < pie2; i++){
				img3Din->FastSet(i,piy, piz, replacement_color);
			}
			for(int i = piw.x; i > piw2; i--){
				img3Din->FastSet(i, piy, piz, replacement_color);
			}
			for(int i = piw2+1; i < pie2; i++){
				piq.x = i; piq.y = piy+1; piq.z = piz; 
				if(piq.y < img3Din->Num_Rows()){
					piq.value = img3Din->FastGet(piq.x, piq.y, piq.z);
					if(piq.value == target_color){
						pixelq.push(piq);
					}
				}
				piq.y=piy-1;
				if(piq.y > 0){
					piq.value = img3Din->FastGet(piq.x, piq.y,piq.z);
					if(piq.value == target_color){
						pixelq.push(piq);
					}
				}
				piq.y = piy;
				piq.z = piz -1;
				if(piq.z > 0){
					piq.value = img3Din->FastGet(piq.x, piq.y, piq.z);
					if(piq.value == target_color){
						pixelq.push(piq);
					}
				}
				piq.z = piz + 1;
				if(piq.z < img3Din->Num_Levels()){
					piq.value = img3Din->FastGet(piq.x, piq.y, piq.z);
					if(piq.value == target_color){
						pixelq.push(piq);
					}
				}

			}
		}
	}

	for (int x = 0; x<img3Din->Num_Cols(); x++)
		for (int y = 0; y<img3Din->Num_Rows(); y++)
			for (int z = 0; z<img3Din->Num_Levels(); z++)
			{
				if(img3Din->FastGet(x,y,z)>1)
				{
					img3Din->FastSet(x,y,z,0);
				}
			}
	for (int x = 0; x<img3Dorig->Num_Cols(); x++)
		for (int y = 0; y<img3Dorig->Num_Rows(); y++)
			for (int z = 0; z<img3Dorig->Num_Levels(); z++)
				img3Din->FastSet(x,y,z,img3Dorig->FastGet(x,y,z) + img3Din->FastGet(x,y,z));

	for (int x = 0; x<img3Din->Num_Cols(); x++)
		for (int y = 0; y<img3Din->Num_Rows(); y++)
			for (int z = 0; z<img3Din->Num_Levels(); z++)
			{
				if(img3Din->FastGet(x,y,z)>0)
				{
					img3Din->FastSet(x,y,z,1000);
				}
			}
	Write_Analyze_File("AfterFloodFill2.hdr",*img3Din);
	std::cout<<"HoleFill Step 2 Complete!!"<<std::endl;
	delete img3Dorig;
}