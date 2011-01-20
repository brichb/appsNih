#include <Polyp_Painting_Structure.h>
#include "NIH_SpineSegmentation_DataStructure.h"
#include <NIH_Algo_Region_Growing.h>
#include <qfile.h>
#include <qstringlist.h>
#include <CIS_Image_Processing_Algo_2D.h>
#include <CIS_Image_Processing_Algo_3D.h>
#include <CIS_Algo_ComputeFeatures.h>
#include "SvmCommittee.h"
#include <CIS_Matrix_JY.h>
#include <CIS_Vector_JY.h>
#include <numerical.h>
#include <Numerical_Recipe_Lib.h>
#include <CIS_Curve.h>
#include <CIS_Algo_Contour.h>
#include <CIS_2D_Model_Advance.h>
#include <CIS_2D_ROI.h>
#include <math.h>

const double PI = 3.1415926;

int CIS_Algo_ComputeStatisticalMoments(doubleDynArray &data, double &mean, double &variance, 
									   double &skewness, double &kurtosis, 
									   double &absoluteDeviation, double &standardDeviation);
int NIH_BoneSegmentation_ComputeMetasisFeatures2D(FeatureStructure &features, int detectionSlice, intVec2DynArray &detectionRegion, 
												  CIS_Array_Image3D_short *img3D,
													CIS_Array_Image3D_short *maskImg3D,
													SpineSegmentationMaskStructure &maskStruct,
													SpineSegmentationParameter &segPara,
													SpineSegmentationInfo &segInfo);
int NIH_BoneSegmentation_ComputeMetasisFeatures3D(FeatureStructure &features3D,
												  intVec3DynArray &detectionRegion,
												  CIS_Array_Image3D_short *img3D,
												  CIS_Array_Image3D_short *maskImg3D,
												  SpineSegmentationMaskStructure &maskStruct,
												  SpineSegmentationParameter &segPara,
												  SpineSegmentationInfo &segInfo);
int NIH_BoneSegmentation_MatchDetections2D(intVec2DynArray &detectionRegion, int detectionSlice,
									   LesionStructure *lesions, int numLesions,
									   PaintingStruct *lesionVoxels,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationParameter &segPara,
									   int &matchedLesion, int &matchedLesionType, float &matchedSize, float &matchedOverlap);
int NIH_BoneSegmentation_MatchDetections3D(intVec3DynArray &detectionRegion,
										   LesionStructure *lesions, int numLesions,
										   PaintingStruct *lesionVoxels,
										   CIS_Array_Image3D_short *maskImg3D,
										   SpineSegmentationParameter &segPara,
										   int &matchedLesion, int &matchedLesionType, float &matchedLesionSize, float &matchedOverlap);
int BernsteinSmoothing(doubleDynArray &xt, doubleDynArray &yt, doubleDynArray &yt2, int bernstein_power);
int PiecewiseBernsteinSmoothing(doubleDynArray &xt, doubleDynArray &yt, doubleDynArray &yt2, int bernstein_power, int piece_size);
int BSplineSmoothing(doubleDynArray &xt, doubleDynArray &yt, doubleDynArray &yt2, int bspline_degree);
double PedicleTemplateMatching(CIS_Array_Image3D_short *img3D, CIS_Array_Image3D_short *maskImg3D, int slice, Vec2 candidate,
							   SpineSegmentationParameter &segPara, SpineSegmentationInfo &segInfo, SpineSegmentationMaskStructure &maskStruct,
							   VertebraStruct2D *vertebraTemplate, int side, bool filled
							   );
float ComputePartitionAccumulateProfile(CIS_Array_Image3D_short *img3D,
									    CIS_Array_Image3D_short *maskImg3D,
									    SpineSegmentationInfo &segInfo, 
										SpineSegmentationMaskStructure &maskStruct,
									    CIS_Array_Image2D_short *img2D_sag,
									    CIS_Array_Image2D_short *img2D_sag_mask,
										vec3DynArray &spineCenter, doubleDynArray &spineWidthLeft, doubleDynArray &spineWidthRight,
										Vec2 &startPos, Vec2 &cNormal, Vec2 &cPerp, 
										int searchRange, float estimate_gap,
										int &profileLength);

int NIH_BoneSegmentation_SaveSegmentation(const char *saveSegPath, const char *prefix_fn, CIS_Array_Image3D_short *maskImg3D)
{
	if(maskImg3D==NULL) return CIS_ERROR;

	char img_fn[400];
	sprintf(img_fn,"%s\\%s.img", saveSegPath, prefix_fn);

	maskImg3D->Save(img_fn);

	return CIS_OK;
}

int NIH_BoneSegmentation_LoadSegmentation(const char *loadSegPath, const char *prefix_fn, CIS_Array_Image3D_short *&maskImg3D)
{
	char img_fn[400];
	sprintf(img_fn,"%s\\%s.img", loadSegPath, prefix_fn);

	QFile qf;
	if(!qf.exists(img_fn)) return CIS_ERROR;

	if(maskImg3D!=NULL) delete maskImg3D;
	maskImg3D = new CIS_Array_Image3D_short();
	maskImg3D->Load(img_fn);

	return CIS_OK;
}

// Preprocessing step for spine segmentation
// Designed for low resolution images ???
int NIH_SpineSegmentation_PreProcess(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo)
{
	int x, y, z, k, k2;
	int sizex, sizey, sizez, sizexy, sizexyz;
	short *imgA, *maskA;

	sizex = img3D->Num_Cols();
	sizey = img3D->Num_Rows();
	sizez = img3D->Num_Levels();

	imgA = img3D->GetArray();
	maskA = maskImg3D->GetArray();

	sizexyz = sizex*sizey*sizez;
	sizexy = sizex*sizey;

	// set the bounding box, could be set manually outside
	segInfo.bound1.z = 0;
	segInfo.bound2.z = sizez-1;

	int stx, sty, stz, edx, edy, edz;
//	if(segInfo.bound1.x<segInfo.bound2.x) 
//	{
//		stx = segInfo.bound1.x;	edx = segInfo.bound2.x;
//	}
//	else
//	{
//		stx = segInfo.bound2.x;	edx = segInfo.bound1.x;
//	}
//	if(segInfo.bound1.y<segInfo.bound2.y) 
//	{
//		sty = segInfo.bound1.y; edy = segInfo.bound2.y;
//	}
//	else
//	{
//		sty = segInfo.bound2.y;	edy = segInfo.bound1.y;
//	}

	if(segInfo.bound1.z<segInfo.bound2.z) 
	{
		stz = segInfo.bound1.z;	edz = segInfo.bound2.z;
	}
	else
	{
		stz = segInfo.bound2.z;	edz = segInfo.bound1.z;
	}

	// restrict the region to lower center part of the images 
	// to segment just the spine, and get ride of heart region
//	if(sty<segPara.predefinedBound1.y) sty=segPara.predefinedBound1.y;
//	if(edy>segPara.predefinedBound2.y || edy<sty) edy=segPara.predefinedBound2.y;
//	if(stx<segPara.predefinedBound1.x || edx<stx) stx=segPara.predefinedBound1.x;
//	if(edx>segPara.predefinedBound2.x || edx<stx) edx=segPara.predefinedBound2.x;

	stx=segPara.predefinedBound1.x;
	edx=segPara.predefinedBound2.x;
	sty=segPara.predefinedBound1.y;
	edy=segPara.predefinedBound2.y;


	// initialization
	// roughly categorize the pixels into three types: air, cortical bone and body based on intensity
	for(k=0; k<sizexyz; k++)
	{
		if(imgA[k]<segPara.airThresh)
		{
			maskA[k] = maskStruct.mask_air;
		}
		else if(imgA[k]>segPara.boneThresh)
		{
			maskA[k] = maskStruct.mask_corticalBone;
		}
		else
		{
			maskA[k] = maskStruct.mask_body;
		}
	}

	// close the gaps
	CIS_Array_Image3D_short *binImg;
	CIS_Array_Image3D_short *blobImg;
	short *binArray;
	short *blobArray;

	binImg = new CIS_Array_Image3D_short(sizex, sizey, sizez);
	blobImg = new CIS_Array_Image3D_short(sizex, sizey, sizez);
	binArray = binImg->GetArray();
	blobArray = blobImg->GetArray();

	// locate the bone pixels slice by slice
	//
	int sty1, ty, by, gap,ty1, by1;
	for(z=stz; z<=edz; z++)
	{
		// refine sty
		ty1 = 0;
		by1 = sizey;

		// first locate the top and bottom of body by scanning the middle line for a certain thick slab (30 pixels)
		x=sizex/2;
		gap = 0;
		k = z*sizexy+x;
		for(y=0; y<sizey; y++, k+=sizex)
		{
			if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_body)
			{
				gap++;
			}
			else gap =0;

			if(gap==30) break;
		}
		if(gap==30) ty1=y-gap;
		
		x=sizex/2;
		gap = 0;
		k = z*sizexy+x+sizex*(sizey-1);
		for(y=sizey-1; y>=0; y--, k-=sizex)
		{
			if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_body)
			{
				gap++;
			}
			else gap =0;

			if(gap==30) break;
		}
		if(gap==30) by1=y+gap;
		
		ty = sizey;
		by = 0;
		for(y=ty1; y<by1; y++)
		{
			k = z*sizexy+y*sizex;
			for(x=0; x<sizex; x++, k++)
			{
				if(maskA[k]==maskStruct.mask_corticalBone)
				{
					if(y<ty) ty=y;
					if(y>by) by=y;
				}
			}
		}	// for y

		sty1 = (by+ty1)/2;
		if(sty1>sty) sty1=sty;
///		if(sty1>ty) sty1=ty;

		for(y=sty1; y<=by1; y++)
		{
			for(x=stx; x<=edx; x++)
			{
				k=z*sizexy+y*sizex+x;
				if(maskA[k]==maskStruct.mask_corticalBone) binArray[k]=1;
				else binArray[k]=0;
			}
		}
//		if(stx<segInfo.bound1.x) segInfo.bound1.x=stx;
//		if(edx>segInfo.bound2.x) segInfo.bound2.x=stx;
//		if(sty1<segInfo.bound1.y) segInfo.bound1.y=sty1;
//		if(edy>segInfo.bound2.y) segInfo.bound2.x=edy;
	}	// for z


	// find the largest connected component of bones
	//
	intDynArray blobRank;
	intVec3DynArray blobCentroid;
	int largestBlob, largestSize;

	NIH_Algo_Blob_Labelling_3D(binImg, blobImg, blobRank, blobCentroid, false);

	largestSize = 0;
	largestBlob = -1;
	for(k=0; k<blobRank.GetSize(); k++) 
	{
		if(blobRank[k]>largestSize) 
		{
			largestSize=blobRank[k];
			largestBlob = k+1;
		}
	}

	for(k=0; k<sizexyz; k++)
	{
		if(maskA[k]==maskStruct.mask_corticalBone)
		{
			if(binArray[k]==0) maskA[k]=maskStruct.mask_otherBone;
			else if(blobArray[k]!=largestBlob)
			{
				maskA[k]=maskStruct.mask_otherBone;	// elinimate other bones
				binArray[k] = 0;
			}
		}
	}

	// do the searching slice by slice
	CIS_Array_Image2D_short *binImg2D, *blobImg2D;
	short *binArray2D, *blobArray2D;
	intDynArray blobRank2D;
	intVec2DynArray blobCentroid2D;

	binImg2D = new CIS_Array_Image2D_short(sizex, sizey);
	blobImg2D = new CIS_Array_Image2D_short(sizex, sizey);
	binArray2D = binImg2D->GetArray();
	blobArray2D = blobImg2D->GetArray();
	int count;

	// close the gap inside bone, slice by slice
	for(z=stz; z<=edz; z++)
	{
		k = z*sizexy;
		count =0;
		for(k2=0, y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++, k2++)
			{
				if(maskA[k]==maskStruct.mask_corticalBone)
				{
					binArray2D[k2] = 1;
					count++;
				}
				else binArray2D[k2]=0;
			}
		}

		if(count==0) continue;

		// close the gap on cortical bone
		CIS_IPA_Dilate(binImg2D, segPara.closeCorticalHoleIteration, true);
		CIS_IPA_Erode(binImg2D, segPara.closeCorticalHoleIteration, true);

		k = z*sizexy;
		for(k2=0, y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++, k2++)
			{
				if(binArray2D[k2]==1) 
				{
					if(maskA[k]!=maskStruct.mask_corticalBone)
						maskA[k]=maskStruct.mask_spongyBone;
				}
			}
		}
	}	// for z

	for(z=1; z<sizez-1; z++)
	{
		// restore some missing disk based on interpolation
		k=z*sizexy;
		for(y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++)
			{
				if(maskA[k]==maskStruct.mask_body && 
					(maskA[k+sizexy]==maskStruct.mask_corticalBone || maskA[k+sizexy]==maskStruct.mask_spongyBone) &&
					(maskA[k-sizexy]==maskStruct.mask_corticalBone || maskA[k-sizexy]==maskStruct.mask_spongyBone))
						maskA[k] = maskStruct.mask_spongyBone;
			}
		}
	}	// for z

/*	// remove excessive bones 

	// compute the centroid of the spine
	Vec3 centroid;
	centroid = Vec3(0,0,0);
	count =0;
	for(z=0, k=0; z<sizez; z++)
	{
		for(y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++)
			{
				if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_spongyBone)
				{
					centroid +=Vec3(x,y,z);
					count++;
				}
			}
		}
	}
	if(count!=0) centroid /= (double)count;

  
	intDynArray blobStatus;
	double largestDist, dist;

	for(z=stz; z<=edz; z++)
	{
		k = z*sizexy;
		count =0;
		for(k2=0, y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++, k2++)
			{
				if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_spongyBone)
				{
					binArray2D[k2] = 1;
					count++;
				}
				else binArray2D[k2]=0;
			}
		}

		if(count==0) continue;

		NIH_Algo_Blob_Labelling(binImg2D, blobImg2D, blobRank2D, blobCentroid2D, false);

		largestSize = 0;
		largestBlob = -1;
		for(k2=0; k2<blobRank2D.GetSize(); k2++) 
		{
			if(blobRank2D[k2]>largestSize) 
			{
				largestSize=blobRank2D[k2];
				largestBlob = k2+1;
			}
		}

		largestDist = (blobCentroid2D[largestBlob-1].x-centroid.x)*(blobCentroid2D[largestBlob-1].x-centroid.x)+
						(blobCentroid2D[largestBlob-1].y-centroid.y)*(blobCentroid2D[largestBlob-1].y-centroid.y);

		// only keep largest blob or those blob closer to the centroid of entire spine
		blobStatus.SetSize(blobRank2D.GetSize());
		for(k2=0; k2<blobStatus.GetSize(); k2++) 
		{
			dist = (blobCentroid2D[k2].x-centroid.x)*(blobCentroid2D[k2].x-centroid.x)+
						(blobCentroid2D[k2].y-centroid.y)*(blobCentroid2D[k2].y-centroid.y);
			if(dist<largestDist || k2==largestBlob-1) blobStatus[k2]=1;
			else blobStatus[k2]=0;
		}

		// remove other bones
		k = z*sizexy;
		for(k2=0, y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++, k2++)
			{
				if(binArray2D[k2]==1) 
				{
					if(blobStatus[blobArray2D[k2]-1]==0)
					{
						if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_spongyBone)
								maskA[k]=maskStruct.mask_otherBone;
					}
				}	// if binArray
			}
		}
	}	// for z
*/
	// compute the bone density

	// Compute some statistics of bones
	int count_c, count_s, count_b;
	segInfo.avg_cortical_intensity = segInfo.avg_spongy_intensity = segInfo.avg_bone_intensity = 0;
	count_c = count_s = count_b = 0;

	for(k=0; k<sizexyz; k++)
	{
		if(maskA[k]==maskStruct.mask_corticalBone)
		{
			count_c++;
			count_b++;
			segInfo.avg_cortical_intensity += imgA[k];
			segInfo.avg_bone_intensity += imgA[k];
		}
		else if(maskA[k]==maskStruct.mask_spongyBone)
		{
			count_s++;
			count_b++;
			segInfo.avg_spongy_intensity += imgA[k];
			segInfo.avg_bone_intensity += imgA[k];
		}
	}

	if(count_c!=0) segInfo.avg_cortical_intensity /= (float)count_c;
	if(count_s!=0) segInfo.avg_spongy_intensity /= (float)count_s;
	if(count_b!=0) segInfo.avg_bone_intensity /= (float)count_b;

	delete binImg;
	delete binImg2D;
	delete blobImg;
	delete blobImg2D;

	return CIS_OK;
}


// a new method for pre-processing
// Designed for high resolution data ???
// note: the algorithm only work for supine cases for now
//
int NIH_SpineSegmentation_PreProcess_new(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo)
{
	int x, y, z, k, k2;
	int sizex, sizey, sizez, sizexy, sizexyz;
	short *imgA, *maskA;
	int count;

	sizex = img3D->Num_Cols();
	sizey = img3D->Num_Rows();
	sizez = img3D->Num_Levels();

	imgA = img3D->GetArray();
	maskA = maskImg3D->GetArray();

	if(sizez>500) sizez=500;	// constraint the size to address memory issue;

	sizexyz = sizex*sizey*sizez;
	sizexy = sizex*sizey;

	// set the bounding box, could be set manually outside
	segInfo.bound1.z = 0;
	segInfo.bound2.z = sizez-1;

	int stx, sty, stz, edx, edy, edz;

	if(segInfo.bound1.z<segInfo.bound2.z) 
	{
		stz = segInfo.bound1.z;	edz = segInfo.bound2.z;
	}
	else
	{
		stz = segInfo.bound2.z;	edz = segInfo.bound1.z;
	}

	if(stz<0) stz=0;
	if(edz>sizez-1) edz=sizez-1;

	stx=segPara.predefinedBound1.x;
	edx=segPara.predefinedBound2.x;
	sty=segPara.predefinedBound1.y;
	edy=segPara.predefinedBound2.y;


	int bbx0, bbx1, bby0, bby1;
	bbx0=sizex; bbx1=0;
	bby0=sizey; bby1=0;
	// initialization using initial bone threshold
	// categorize the pixels into three classes: air, cortical bone and body
	for(z=0, k=0; z<sizez; z++)
		for(y=0; y<sizey; y++)
			for(x=0; x<sizex; x++, k++)
	{
		if(imgA[k]<segPara.airThresh)
		{
			maskA[k] = maskStruct.mask_air;
		}
		else if(imgA[k]>segPara.boneThresh)
		{
			maskA[k] = maskStruct.mask_corticalBone;
			if(x<bbx0) bbx0=x;
			if(x>bbx1) bbx1=x;
			if(y<bby0) bby0=y;
			if(y>bby1) bby1=y;
		}
		else
		{
			maskA[k] = maskStruct.mask_body;
		}
	}
	bbx0-=2; bbx1+=2; bby0-=2; bby1+=2;
	if(bbx0<0) bbx0=0;
	if(bbx1>sizex-1) bbx1=sizex-1;
	if(bby0<0) bby0=0;
	if(bby1>sizey-1) bby1=sizey-1;

	int subSizex, subSizey, sub_sizexy;
	subSizex=bbx1-bbx0+1;
	subSizey=bby1-bby0+1;
	sub_sizexy = subSizex*subSizey;

	CIS_Array_Image3D_uchar *binImg;
	CIS_Array_Image3D_short *blobImg;
	unsigned char *binArray;
	short *blobArray;

	binImg = new CIS_Array_Image3D_uchar(subSizex, subSizey, sizez);
	blobImg = new CIS_Array_Image3D_short(subSizex, subSizey, sizez);
	binArray = binImg->GetArray();
	blobArray = blobImg->GetArray();

	// do the searching slice by slice
	CIS_Array_Image2D_short *binImg2D, *blobImg2D;
	short *binArray2D, *blobArray2D;
	intDynArray blobRank2D;
	intVec2DynArray blobCentroid2D;

	binImg2D = new CIS_Array_Image2D_short(subSizex, subSizey);
	blobImg2D = new CIS_Array_Image2D_short(subSizex, subSizey);
	binArray2D = binImg2D->GetArray();
	blobArray2D = blobImg2D->GetArray();

	// close the gaps inside bone slice by slice
	int sty1, ty, by, gap,ty1, by1;
	for(z=stz; z<=edz; z++)
	{
		// refine sty
		ty1 = 0;
		by1 = sizey;

		// first locate the top and bottom of body by scanning the middle line for a certain thick slab (30 pixels)
		// work for supine cases only
		x=sizex/2;
		gap = 0;
		for(y=bby0, k=z*sizexy+y*sizex+x; y<=bby1; y++, k+=sizex)
		{
			if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_body)
			{
				gap++;
			}
			else gap =0;

			if(gap==30) break;
		}
		if(gap==30) ty1=y-gap+10;
		
		x=sizex/2;
		gap = 0;
		for(y=bby1, k=z*sizexy+y*sizex+x; y>=bby0; y--, k-=sizex)
		{
			if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_body)
			{
				gap++;
			}
			else gap =0;

			if(gap==30) break;
		}
		if(gap==30) by1=y+gap-10;
		
		// further refine by and ty
		ty = sizey;
		by = 0;
		int bcount;
		for(y=ty1; y<by1; y++)
		{
			k = z*sizexy+y*sizex+stx;
			bcount=0;
			for(x=stx; x<=edx; x++, k++)
			{
				if(maskA[k]==maskStruct.mask_corticalBone)
				{
					bcount++;
				}
				if(bcount>10)
				{
					if(y<ty) ty=y;
					if(y>by) by=y;
					break;
				}
			}
		}	// for y

//		if(z<sizez/8) sty1 = (ty+ty1)/2;		// move the region higher for thorasis spine
//		else 
			sty1 = (by+ty1)/2-10;
		if(sty1>sty) sty1=sty;
///		if(sty1>ty) sty1=ty;

		for(k2=0; k2<sub_sizexy; k2++) binArray2D[k2]=0;
		
		if(sty1<bby0) sty1=bby0;
		if(stx<bbx0) stx=bbx0;
		if(by1>bby1) by1=bby1;
		if(edx>bbx1) edx=bbx1;

		count =0;
		for(y=sty1; y<=by1; y++)
		{
			for(x=stx, k=z*sizexy+y*sizex+x, k2=(y-bby0)*subSizex+(x-bbx0); x<=edx; x++, k++, k2++)
			{
				if(maskA[k]==maskStruct.mask_corticalBone)
				{
					binArray2D[k2] = 1;
					count++;
				}
				else binArray2D[k2]=0;
			}
		}

		if(count==0) continue;
		// close the gap on cortical bone
		CIS_IPA_Dilate(binImg2D, segPara.closeCorticalHoleIteration, true);
		CIS_IPA_Erode(binImg2D, segPara.closeCorticalHoleIteration, true);

		k=z*sub_sizexy;
		for(k2=0; k2<sub_sizexy; k2++) binArray[k+k2] = binArray2D[k2];

	}	// for z


	// find the largest connected component
	intDynArray blobRank;
	intVec3DynArray blobCentroid;
	int largestBlob, largestSize;

	NIH_Algo_Blob_Labelling_3D(binImg, blobImg, blobRank, blobCentroid, false);

	largestSize = 0;
	largestBlob = -1;
	for(k=0; k<blobRank.GetSize(); k++) 
	{
		if(blobRank[k]>largestSize) 
		{
			largestSize=blobRank[k];
			largestBlob = k+1;
		}
	}

	// eliminate other bones
	for(z=0, k2=0; z<sizez; z++)
		for(y=bby0; y<=bby1; y++)
			for(x=bbx0, k=z*sizexy+y*sizex+x; x<=bbx1; x++, k++, k2++)
	{
		if(maskA[k]==maskStruct.mask_corticalBone)
		{
			if(binArray[k2]==0) maskA[k]=maskStruct.mask_otherBone;
			else if(blobArray[k2]!=largestBlob)
			{
				maskA[k]=maskStruct.mask_otherBone;	
				binArray[k2] = 0;
			}
		}
	}

	// close gap slice by slice
	for(z=stz; z<=edz; z++)
	{
		count =0;
		for(k2=0, y=bby0; y<=bby1; y++)
		{
			for(x=bbx0, k=z*sizexy+y*sizex+x; x<=bbx1; x++, k++, k2++)
			{
				if(maskA[k]==maskStruct.mask_corticalBone)
				{
					binArray2D[k2] = 1;
					count++;
				}
				else binArray2D[k2]=0;
			}
		}

		if(count==0) continue;

		// close the gap on cortical bone
		CIS_IPA_Dilate(binImg2D, segPara.closeCorticalHoleIteration, true);
		CIS_IPA_Erode(binImg2D, segPara.closeCorticalHoleIteration, true);

		for(k2=0, y=bby0; y<=bby1; y++)
		{
			for(x=bbx0, k=z*sizexy+y*sizex+x; x<=bbx1; x++, k++, k2++)
			{
				if(binArray2D[k2]==1) 
				{
					if(maskA[k]!=maskStruct.mask_corticalBone)
						maskA[k]=maskStruct.mask_spongyBone;
				}
			}
		}
	}	// for z

	for(z=1; z<sizez-1; z++)
	{
		// restore some missing disk based on interpolation
		k=z*sizexy;
		for(y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++)
			{
				if(maskA[k]==maskStruct.mask_body && 
					(maskA[k+sizexy]==maskStruct.mask_corticalBone || maskA[k+sizexy]==maskStruct.mask_spongyBone) &&
					(maskA[k-sizexy]==maskStruct.mask_corticalBone || maskA[k-sizexy]==maskStruct.mask_spongyBone))
						maskA[k] = maskStruct.mask_spongyBone;
			}
		}
	}	// for z

	// compute the bone density
	int count_c, count_s, count_b;
	segInfo.avg_cortical_intensity = segInfo.avg_spongy_intensity = segInfo.avg_bone_intensity = 0;
	count_c = count_s = count_b = 0;

	for(k=0; k<sizexyz; k++)
	{
		if(maskA[k]==maskStruct.mask_corticalBone)
		{
			count_c++;
			count_b++;
			segInfo.avg_cortical_intensity += imgA[k];
			segInfo.avg_bone_intensity += imgA[k];
		}
		else if(maskA[k]==maskStruct.mask_spongyBone)
		{
			count_s++;
			count_b++;
			segInfo.avg_spongy_intensity += imgA[k];
			segInfo.avg_bone_intensity += imgA[k];
		}
	}

	if(count_c!=0) segInfo.avg_cortical_intensity /= (float)count_c;
	if(count_s!=0) segInfo.avg_spongy_intensity /= (float)count_s;
	if(count_b!=0) segInfo.avg_bone_intensity /= (float)count_b;

	delete binImg;
	delete binImg2D;
	delete blobImg;
	delete blobImg2D;

	return CIS_OK;
}



// detect spinal cord and set up bounding box for each part of the vertebra (vetrabra body and spinal process)
// segment the vertebra structure using a template
// Designed for high res images??
int NIH_SpineSegmentation_DetectSpinalCord_new(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   VertebraStruct2D *vertebraTemplate, 
									   bool debugMode)
{
	int x, y, z, k, i;
	int sizex, sizey, sizez, sizexy, sizexyz;
	short *imgA, *maskA;
	float pixelSizez=img3D->Get_Pixel_SizeZ();

	if(pixelSizez>3) segPara.cordIntensityThresh=1100;	// tmp changed for low resolution data

	sizex = img3D->Num_Cols();
	sizey = img3D->Num_Rows();
	sizez = img3D->Num_Levels();

	imgA = img3D->GetArray();
	maskA = maskImg3D->GetArray();

	sizexyz = sizex*sizey*sizez;
	sizexy = sizex*sizey;

	// for debug purpose, reset the maskA
///	for(k=0; k<sizexyz; k++) if(maskA[k]==maskStruct.mask_spinalCord) maskA[k]=maskStruct.mask_body;

	CIS_Array_Image2D_short *binImg2D, *blobImg2D;
	short *binArray2D, *blobArray2D;
	intDynArray blobRank2D;
	intVec2DynArray blobCentroid2D;
	int largestBlob, largestSize;
	intDynArray blobStatus;

	int count, k2;
	int stx, sty, edx, edy;

	binImg2D = new CIS_Array_Image2D_short(sizex, sizey);
	blobImg2D = new CIS_Array_Image2D_short(sizex, sizey);
	binArray2D = binImg2D->GetArray();
	blobArray2D = blobImg2D->GetArray();

	for(z=0; z<sizez; z++)
	{
		vertebraTemplate[z].cordRadius = Vec2(-1,-1);
		vertebraTemplate[z].cordCenter = Vec2(-1,-1);
		vertebraTemplate[z].cord_bb0 = IntVec2(-1,-1);
		vertebraTemplate[z].cord_bb1 = IntVec2(-1,-1);
		vertebraTemplate[z].cordSize = 0;
		vertebraTemplate[z].isKey = -1;
		for(int a=0; a<diskAngleCount; a++)
		{
			vertebraTemplate[z].cordContour[a] = Vec2(-1, -1);
			vertebraTemplate[z].cordContourRadius[a] = -1;
		}
	}

	segInfo.bound1.z = -1;
	// first pass to get the candidate cord locations
	//
	for(z=0; z<sizez; z++)
	{
		// first define a bounding box
		k = z*sizexy;
		count =0;
		stx = sizex-3;
		sty = sizey-3;
		edx = 3;
		edy = 3;
		for(k2=0, y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++, k2++)
			{
				if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_spongyBone)
				{
					if(x<stx) stx=x;
					if(x>edx) edx=x;
					if(y<sty) sty=y;
					if(y>edy) edy=y;
					count++;
				}
			}
		}

		if(count<=segPara.cordSizeThresh) continue;

		// move the bounding box up if there are bone pixels above
		k = z*sizexy;
		for(y=sty; y>0 && y>sty-(edy-sty)/2; y--)
		{
			count=0;
			for(x=stx; x<edx; x++)
			{
				if(imgA[k+y*sizex+x]>segPara.boneThresh) count++;
			}
			if(count<5) break;
		}
		sty = y;

		// first pass eliminate external space
		stx-=6; sty-=6; edx+=6; edy+=6;
		if(stx<0) stx=0;
		if(sty<0) sty=0;
		if(edx>=sizex-1) edx=sizex-2;
		if(edy>=sizey-1) edy=sizey-2;

		int sub_sizex, sub_sizey;
		int sk, sx, sy, nx, ny;
		int cb;	// current blob
		CIS_Array_Image2D_short *sub_binImg2D, *sub_blobImg2D;
		short *subBinA, *subBlobA;
		sub_sizex = edx-stx+1;
		sub_sizey = edy-sty+1;

		sub_binImg2D = new CIS_Array_Image2D_short(sub_sizex, sub_sizey);
		sub_blobImg2D = new CIS_Array_Image2D_short(sub_sizex, sub_sizey);
		subBinA = sub_binImg2D->GetArray();
		subBlobA = sub_blobImg2D->GetArray();

		k = z*sizexy;
		count=0;
		for(y=sty, sk=0; y<=edy; y++)
		{
			k2 = y*sizex+stx;
			for(x=stx; x<=edx; x++, k2++, sk++)
			{
				if(maskA[k+k2]==maskStruct.mask_body && imgA[k+k2]<segPara.cordIntensityThresh)
				{
					subBinA[sk]=0;
					count ++;
				}
				else subBinA[sk]=1;
			}
		}	// for y

		if(count<segPara.cordSizeThresh) 
		{
			delete sub_binImg2D;
			delete sub_blobImg2D;
			continue;
		}

		// remove small pieces
		NIH_Algo_Blob_Labelling(sub_binImg2D, sub_blobImg2D, blobRank2D, blobCentroid2D, false);
		for(y=0, sk=0; y<sub_sizey; y++)
		{
			for(x=0; x<sub_sizex; x++, sk++)
			{
				cb = subBlobA[sk]-1;
				if(cb<=0) continue;
				if(blobRank2D[cb]<10) subBinA[sk]=0;
			}
		}
		// revert the binary image
		for(y=0, sk=0; y<sub_sizey; y++)
		{
			for(x=0; x<sub_sizex; x++, sk++)
			{
				if(subBinA[sk]==0) subBinA[sk]=1;
				else subBinA[sk]=0;
			}
		}

//		// close the hole inside
//		CIS_IPA_Dilate(sub_binImg2D, 1, false);
//		CIS_IPA_Erode(sub_binImg2D, 1, false);
		// separate cord and external space
		CIS_IPA_Erode(sub_binImg2D, segPara.closeCordIteration, true);
		CIS_IPA_Dilate(sub_binImg2D, segPara.closeCordIteration, true);

		NIH_Algo_Blob_Labelling(sub_binImg2D, sub_blobImg2D, blobRank2D, blobCentroid2D, false);

		blobStatus.SetSize(blobRank2D.GetSize());
		for(k=0; k<blobStatus.GetSize(); k++) blobStatus[k]=-1;

		intVec2DynArray blob_bb0, blob_bb1;
		intDynArray blob_perimeter, blob_perimeterIntensity;
		floatDynArray blob_aspectRatio, blob_compactness;
		blob_bb0.SetSize(blobStatus.GetSize());
		blob_bb1.SetSize(blobStatus.GetSize());
		blob_perimeter.SetSize(blobStatus.GetSize());
		blob_perimeterIntensity.SetSize(blobStatus.GetSize());
		blob_aspectRatio.SetSize(blobStatus.GetSize());
		blob_compactness.SetSize(blobStatus.GetSize());
		for(k=0; k<blobStatus.GetSize(); k++) 
		{
			blob_bb0[k]=IntVec2(edx, edy);
			blob_bb1[k]=IntVec2(stx, sty);
		}

		// Get the outer space, any blob next to border are considered outer space
		k = z*sizexy;
		count=0;
		for(y=sty, sk=0; y<=edy; y++)
		{
			k2 = y*sizex+stx;
			for(x=stx; x<=edx; x++, k2++, sk++)
			{
				if(subBlobA[sk]>=1)
				{
					cb = subBlobA[sk]-1;
					if(blobStatus[cb]==-2) continue;
					if(y<sty+segPara.closeCordIteration || y>=edy-segPara.closeCordIteration
						|| x<stx+segPara.closeCordIteration || x>=edx-segPara.closeCordIteration)
					{
						blobStatus[cb] = -2;	// outer border blob
					}
					// compute the bounding box of each blob, just y for now
					if(x<blob_bb0[cb].x) blob_bb0[cb].x=x;
					if(x>blob_bb1[cb].x) blob_bb1[cb].x=x;
					if(y<blob_bb0[cb].y) blob_bb0[cb].y=y;
					if(y>blob_bb1[cb].y) blob_bb1[cb].y=y;
				}
			}
		}	// for y

		int validBlob=0;
		for(cb=0; cb<blobStatus.GetSize(); cb++) 
		{
			if(blobStatus[cb]==-2) continue;
			
			if(blobRank2D[cb]<segPara.cordSizeThresh || blobRank2D[cb]>segPara.cordSizeThresh*25
				|| blobCentroid2D[cb].x+stx<sizex/2-sizex/8 
				|| blobCentroid2D[cb].x+stx>sizex/2+sizex/8) blobStatus[cb]=-1;		
			// too small or not in the center, this is sensitive to patient position, may not be at the center
			else
			{
				// the bounding box of cord should be in the center of the bounding box of vertebra
				// it should not near the top (the first quarter of the bounding box
				// this is sensitive to table height
				if(blob_bb0[cb].y<(sty*3+edy)/4) continue;

				blob_aspectRatio[cb] = (float)(blob_bb1[cb].x-blob_bb0[cb].x)/(float)(blob_bb1[cb].y-blob_bb0[cb].y);
				if(blob_aspectRatio[cb]<0.75 || blob_aspectRatio[cb]>2) continue;	// the aspect ratio should meet some criteria

				// do some shape and intensity analysis
				blob_perimeter[cb]=0;
				blob_perimeterIntensity[cb]=0;
				for(sy=blob_bb0[cb].y-2-sty; sy<=blob_bb1[cb].y+2-sty; sy++)
				{
					for(sx=blob_bb0[cb].x-2-stx, sk=sy*sub_sizex+sx; sx<=blob_bb1[cb].x+2-stx; sx++, sk++)
					{
						if(subBlobA[sk]==0)
						{
							k2 = (sy+sty)*sizex+sx+stx;
							if(maskA[k+k2]==maskStruct.mask_air)	// next to air pocket
							{
								blobStatus[cb]=-2;
								break;
							}
							if(subBlobA[sk-1]==cb+1 || subBlobA[sk+1]==cb+1 || subBlobA[sk-sub_sizex]==cb+1 
								|| subBlobA[sk+sub_sizex]==cb+1)
							{
								blob_perimeter[cb]++;
								blob_perimeterIntensity[cb]+=imgA[k+k2];
							}
						}
					}
				}	// for sy
				if(blobStatus[cb]==-2) continue;

				if(blob_perimeter[cb]>0) blob_perimeterIntensity[cb]/=blob_perimeter[cb];
				blob_compactness[cb] = blob_perimeter[cb]*blob_perimeter[cb]/blobRank2D[cb];
				if(blob_compactness[cb]>20) continue;
///				if(blob_perimeterIntensity[cb]<segPara.boneThresh) continue;	// not surrounded by cortical bones

				blobStatus[cb]=1;
				validBlob++;
			}
		}	// for cb
		
		if(validBlob==0) 
		{
			delete sub_binImg2D;
			delete sub_blobImg2D;
			continue;
		}

		// if more than one blob are valid, choose the one with largest perimeter intensity
		if(validBlob>1)
		{
			int largestIntensity=0;
			int largest_cb=-1;
			for(cb=0; cb<blobStatus.GetSize(); cb++) 
			{
				if(blobStatus[cb]==1 && blob_perimeterIntensity[cb]>largestIntensity)
				{
					largestIntensity = blob_perimeterIntensity[cb];
					largest_cb=cb;
				}
			}
			for(cb=0; cb<blobStatus.GetSize(); cb++) 
			{
				if(blobStatus[cb]==1 && cb!=largest_cb)
				{
					blobStatus[cb]=-2;
				}
			}
		}	// if validBlob

		// filled the template data structure
		for(cb=0; cb<blobStatus.GetSize(); cb++) 
		{
			if(blobStatus[cb]==1)
			{
				vertebraTemplate[z].cord_bb0 = blob_bb0[cb];
				vertebraTemplate[z].cord_bb1 = blob_bb1[cb];
				vertebraTemplate[z].cordRadius = (blob_bb1[cb]-blob_bb0[cb])/2;
				vertebraTemplate[z].cordCenter = Vec2(blobCentroid2D[cb].x+stx, blobCentroid2D[cb].y+sty);
				vertebraTemplate[z].isKey = 1;
				vertebraTemplate[z].cordSize = blobRank2D[cb];
			}
		}

		k = z*sizexy;

		// fill in the mask
		for(y=sty, sk=0; y<=edy; y++)
		{
			k2 = y*sizex+stx;
			for(x=stx; x<=edx; x++, k2++, sk++)
			{
				if(subBlobA[sk]>=1)
				{
					cb = subBlobA[sk]-1;
					if(blobStatus[cb]==1)
					{
						maskA[k+k2]=maskStruct.mask_spinalCord;
					}
				}
			}
		}	// for y

		delete sub_binImg2D;
		delete sub_blobImg2D;
	}	// for z

	// further check to remove outliers
	int totalKeySlices=0;

	// first use the size
	Vec2 avgCordRadius;
	double avgCordSize;
	avgCordRadius.x = avgCordRadius.y = 0;
	avgCordSize=0;
	for(z=0; z<sizez; z++)
	{
		if(vertebraTemplate[z].isKey==1)
		{
			avgCordRadius += vertebraTemplate[z].cordRadius;
			avgCordSize += vertebraTemplate[z].cordSize;
			totalKeySlices++;
		}
	}
	if(totalKeySlices>0) 
	{
		avgCordRadius /= (double)totalKeySlices;
		avgCordSize /= (double)totalKeySlices;
	}
	for(z=0; z<sizez; z++)
	{
		if(vertebraTemplate[z].isKey==1)
		{
			if(vertebraTemplate[z].cordRadius.x<avgCordRadius.x/2 || vertebraTemplate[z].cordRadius.y<avgCordRadius.y/2
				|| vertebraTemplate[z].cordRadius.x>avgCordRadius.x*2 || vertebraTemplate[z].cordRadius.y>avgCordRadius.y*2
				|| vertebraTemplate[z].cordSize<avgCordSize/2 )
			{
				vertebraTemplate[z].cord_bb0 = IntVec2(-1,-1);
				vertebraTemplate[z].cord_bb1 = IntVec2(-1,-1);
				vertebraTemplate[z].cordRadius = Vec2(-1,-1);
				vertebraTemplate[z].cordCenter = Vec2(-1,-1);
				vertebraTemplate[z].isKey = -2;
			}
 		}
	}

	// recompute average size and refiltering
	if(totalKeySlices>100)
	{
	totalKeySlices=0;
	avgCordRadius.x = avgCordRadius.y = 0;
	avgCordSize=0;
	for(z=0; z<sizez; z++)
	{
		if(vertebraTemplate[z].isKey==1)
		{
			avgCordRadius += vertebraTemplate[z].cordRadius;
			avgCordSize += vertebraTemplate[z].cordSize;
			totalKeySlices++;
		}
	}
	if(totalKeySlices>0) 
	{
		avgCordRadius /= (double)totalKeySlices;
		avgCordSize /= (double)totalKeySlices;
	}
	for(z=0; z<sizez; z++)
	{
		if(vertebraTemplate[z].isKey==1)
		{
			if(vertebraTemplate[z].cordRadius.x<avgCordRadius.x/2 || vertebraTemplate[z].cordRadius.y<avgCordRadius.y/2
				|| vertebraTemplate[z].cordRadius.x>avgCordRadius.x*2 || vertebraTemplate[z].cordRadius.y>avgCordRadius.y*2
				|| vertebraTemplate[z].cordSize<avgCordSize/2 || vertebraTemplate[z].cordSize>avgCordSize*2)
			{
				vertebraTemplate[z].cord_bb0 = IntVec2(-1,-1);
				vertebraTemplate[z].cord_bb1 = IntVec2(-1,-1);
				vertebraTemplate[z].cordRadius = Vec2(-1,-1);
				vertebraTemplate[z].cordCenter = Vec2(-1,-1);
				vertebraTemplate[z].isKey = -2;
			}
 		}
	}
	}

	// then use the overlap of bounding box as the criteria
	totalKeySlices=0;
	for(z=0; z<sizez; z++)
	{
		if(vertebraTemplate[z].isKey==1)
		{
			int totalUEval=0, totalDEval=0, totalEval=0;
			int goodEval=0;
			int z1;
			// check previous slices
			for(z1=z-1; z1>=0 && totalEval<3 && z-z1<80; z1--)
			{
				if(vertebraTemplate[z1].isKey==1)
				{
					totalEval++;
					totalUEval++;
					if(!(vertebraTemplate[z].cord_bb0.x>vertebraTemplate[z1].cord_bb1.x-vertebraTemplate[z1].cordRadius.x/3 ||
						vertebraTemplate[z1].cord_bb0.x>vertebraTemplate[z].cord_bb1.x-vertebraTemplate[z].cordRadius.x/3) &&
						!(vertebraTemplate[z].cord_bb0.y>vertebraTemplate[z1].cord_bb1.y-vertebraTemplate[z1].cordRadius.y/6 ||
						vertebraTemplate[z1].cord_bb0.y>vertebraTemplate[z].cord_bb1.y-vertebraTemplate[z].cordRadius.y/6) &&
						vertebraTemplate[z].cordRadius.x>vertebraTemplate[z1].cordRadius.x/2 &&
						vertebraTemplate[z].cordRadius.y>vertebraTemplate[z1].cordRadius.y/2)
					{
						goodEval++;
					}
				}
			}
			// check following slices
			for(z1=z+1; z1<sizez && totalEval<6 && z1-z<50; z1++)
			{
				if(vertebraTemplate[z1].isKey==1)
				{
					totalEval++;
					totalDEval++;
					if(!(vertebraTemplate[z].cord_bb0.x>vertebraTemplate[z1].cord_bb1.x-vertebraTemplate[z1].cordRadius.x/3 ||
						vertebraTemplate[z1].cord_bb0.x>vertebraTemplate[z].cord_bb1.x-vertebraTemplate[z].cordRadius.x/3) &&
						!(vertebraTemplate[z].cord_bb0.y>vertebraTemplate[z1].cord_bb1.y-vertebraTemplate[z1].cordRadius.y/6 ||
						vertebraTemplate[z1].cord_bb0.y>vertebraTemplate[z].cord_bb1.y-vertebraTemplate[z].cordRadius.y/6) &&
						vertebraTemplate[z].cordRadius.x>vertebraTemplate[z1].cordRadius.x/2 &&
						vertebraTemplate[z].cordRadius.y>vertebraTemplate[z1].cordRadius.y/2)
					{
						goodEval++;
					}
				}
			}

			if(totalEval<3 || (pixelSizez<=1.5 && totalEval-goodEval>2) 
				|| (pixelSizez<=1.5 && (float)goodEval/(float)totalEval<0.5)
				|| (pixelSizez>1.5 && totalEval-goodEval>4))	// it is a outlier
			{
				vertebraTemplate[z].isKey = -2;
			}	// for goodEval
			else totalKeySlices++;
		}	// if vertebra
	}	// for z

	// something is wrong if not enough key slices
	if(totalKeySlices<5) 
	{
		printf("Error: not enough key slices in cord localization\n");
		return CIS_ERROR;
	}

	// add a few key frame back if they satify the overlap criteria with other key frames
	for(z=0; z<sizez; z++)
	{
		if(vertebraTemplate[z].isKey==-2)
		{
			int totalEval=0;
			int goodEval=0;
			int z1;
			// check previous slices
			for(z1=z-1; z1>=0 && totalEval<5; z1--)
			{
				if(vertebraTemplate[z1].isKey==1)
				{
					totalEval++;
					if(!(vertebraTemplate[z].cord_bb0.x>vertebraTemplate[z1].cord_bb1.x-vertebraTemplate[z1].cordRadius.x/6 ||
						vertebraTemplate[z1].cord_bb0.x>vertebraTemplate[z].cord_bb1.x-vertebraTemplate[z].cordRadius.x/6) &&
						!(vertebraTemplate[z].cord_bb0.y>vertebraTemplate[z1].cord_bb1.y-vertebraTemplate[z1].cordRadius.y/6 ||
						vertebraTemplate[z1].cord_bb0.y>vertebraTemplate[z].cord_bb1.y-vertebraTemplate[z].cordRadius.y/6) &&
						vertebraTemplate[z].cordRadius.x>vertebraTemplate[z1].cordRadius.x/2 &&
						vertebraTemplate[z].cordRadius.y>vertebraTemplate[z1].cordRadius.y/2)
					{
						goodEval++;
					}
				}
			}
			// check following slices
			for(z1=z+1; z1<sizez && totalEval<5; z1++)
			{
				if(vertebraTemplate[z1].isKey==1)
				{
					totalEval++;
					if(!(vertebraTemplate[z].cord_bb0.x>vertebraTemplate[z1].cord_bb1.x-vertebraTemplate[z1].cordRadius.x/6 ||
						vertebraTemplate[z1].cord_bb0.x>vertebraTemplate[z].cord_bb1.x-vertebraTemplate[z].cordRadius.x/6) &&
						!(vertebraTemplate[z].cord_bb0.y>vertebraTemplate[z1].cord_bb1.y-vertebraTemplate[z1].cordRadius.y/6 ||
						vertebraTemplate[z1].cord_bb0.y>vertebraTemplate[z].cord_bb1.y-vertebraTemplate[z].cordRadius.y/6) &&
						vertebraTemplate[z].cordRadius.x>vertebraTemplate[z1].cordRadius.x/2 &&
						vertebraTemplate[z].cordRadius.y>vertebraTemplate[z1].cordRadius.y/2)
					{
						goodEval++;
					}
				}
			}

			if(goodEval>=4)	// it is not an outlier
			{
				vertebraTemplate[z].isKey = 1;
				totalKeySlices++;
			}	// for goodEval
		}	// if vertebra
	}	// for z

	for(z=0; z<sizez; z++)
	{
		if(vertebraTemplate[z].isKey==-2)
		{
				vertebraTemplate[z].cord_bb0 = IntVec2(-1,-1);
				vertebraTemplate[z].cord_bb1 = IntVec2(-1,-1);
				vertebraTemplate[z].cordRadius = Vec2(-1,-1);
				vertebraTemplate[z].cordCenter = Vec2(-1,-1);
		}
	}

	// interpolation to get spinal cord on non-key slices using b-spline fitting
	doubleDynArray xt, yt, yt2;
	double minInterpolate, maxInterpolate, minInterpolate2, maxInterpolate2;
	int bernstein_power = 5, piece_size=5;
	float piece_length=25;
	piece_size = (int)(piece_length/pixelSizez+0.5);

	xt.SetSize(sizez);
	yt.SetSize(sizez);
	yt2.SetSize(sizez);
	for(z=0; z<sizez; z++)
	{
		xt[z] = z;
		yt[z] = -1;
	}

	// do the interpolation twice
	for(int inn=0; inn<2; inn++)
	{
	// interpolate x coordinate of cord
	minInterpolate = -1;
	maxInterpolate = -1;
	for(z=0; z<sizez; z++)
	{
		yt[z] = -1;
		if(vertebraTemplate[z].isKey!=1) continue;

		yt[z] = vertebraTemplate[z].cordCenter.x;
		if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
		if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=0; z<sizez; z++)
	{
		if(yt2[z]==-1) continue;
		if(yt2[z]<minInterpolate) yt2[z]=minInterpolate;
		if(yt2[z]>maxInterpolate) yt2[z]=maxInterpolate;
		if(vertebraTemplate[z].isKey!=1)
			vertebraTemplate[z].cordCenter.x = yt2[z];
		else if(yt2[z]<vertebraTemplate[z].cord_bb0.x || yt2[z]>vertebraTemplate[z].cord_bb1.x)
		{
			vertebraTemplate[z].cordCenter.x = yt2[z];
			vertebraTemplate[z].isKey=-2;		// outliner
		}
	}

	// interpolate y coordinate of cord
	minInterpolate = -1;
	maxInterpolate = -1;
	for(z=0; z<sizez; z++)
	{
		yt[z] = -1;
		if(vertebraTemplate[z].isKey!=1) continue;

		yt[z] = vertebraTemplate[z].cordCenter.y;
		if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
		if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=0; z<sizez; z++)
	{
		if(yt2[z]==-1) continue;
		if(yt2[z]<minInterpolate) yt2[z]=minInterpolate;
		if(yt2[z]>maxInterpolate) yt2[z]=maxInterpolate;
		if(vertebraTemplate[z].isKey!=1)
			vertebraTemplate[z].cordCenter.y = yt2[z];
		else if(yt2[z]<vertebraTemplate[z].cord_bb0.y-vertebraTemplate[z].cordRadius.y 
			|| yt2[z]>vertebraTemplate[z].cord_bb1.y+vertebraTemplate[z].cordRadius.y)
		{
			vertebraTemplate[z].cordCenter.y = yt2[z];
			vertebraTemplate[z].isKey=-2;	// outlier
		}
	}
	}

	// interpolate x radius of cord
	minInterpolate = -1;
	maxInterpolate = -1;
	for(z=0; z<sizez; z++)
	{
		yt[z] = -1;
		if(vertebraTemplate[z].isKey!=1) continue;

		yt[z] = vertebraTemplate[z].cordRadius.x;
		if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
		if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=0; z<sizez; z++)
	{
		if(yt2[z]==-1) continue;
		if(yt2[z]<minInterpolate) yt2[z]=minInterpolate;
		if(yt2[z]>maxInterpolate) yt2[z]=maxInterpolate;
		if(vertebraTemplate[z].isKey!=1)
			vertebraTemplate[z].cordRadius.x = yt2[z];
	}

	// interpolate y radius of cord
	minInterpolate = -1;
	maxInterpolate = -1;
	for(z=0; z<sizez; z++)
	{
		yt[z] = -1;
		if(vertebraTemplate[z].isKey!=1) continue;

		yt[z] = vertebraTemplate[z].cordRadius.y;
		if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
		if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=0; z<sizez; z++)
	{
		if(yt2[z]==-1) continue;
		if(yt2[z]<minInterpolate) yt2[z]=minInterpolate;
		if(yt2[z]>maxInterpolate) yt2[z]=maxInterpolate;
		if(vertebraTemplate[z].isKey!=1)
			vertebraTemplate[z].cordRadius.y = yt2[z];
	}

	int lastKeySlice=sizez-1;
	for(z=sizez-1; z>=0; z--) 
	{
		if(vertebraTemplate[z].isKey==1)
		{
			lastKeySlice=z;
			break;
		}
	}	// for z

	// determine the slice where the sacrum starts
	// our method currently will not segment any slices after the sacrum
	// currently just use the key slice with largest spinal cord as the start of the sacrum
	//
	segInfo.sacrum_start=-1;
	int largest_cord=0;
	for(z=sizez/2; z<sizez; z++)
	{
		if(vertebraTemplate[z].isKey==1)
		{
			if(vertebraTemplate[z].cordSize>=largest_cord)
			{
				largest_cord = vertebraTemplate[z].cordSize;
				segInfo.sacrum_start = z;
			}
		}
	}	// for z
	
	if(segInfo.sacrum_start>0 && segInfo.sacrum_start>100)	// need to have enough slices
	{
		// try to extent to the end of the key slice
		for(z=segInfo.sacrum_start+1; z<sizez; z++)
		{
			if(vertebraTemplate[z].isKey!=1) break;
		}
		segInfo.sacrum_start = z-1;
		for(z=segInfo.sacrum_start+1; z<sizez; z++)
		{
			if(vertebraTemplate[z].isKey==1) vertebraTemplate[z].isKey=-2;
		}
	}
	else segInfo.sacrum_start=lastKeySlice;

	// remove the outliers in the mask
	// 
	for(z=0; z<sizez; z++)
	{
		if(vertebraTemplate[z].isKey==-2)
		{
			k = z*sizexy;

			for(y=0, k2=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k2++)
				{
					if(maskA[k+k2]==maskStruct.mask_spinalCord)
						maskA[k+k2] = maskStruct.mask_body;
				}
			}	// for y
		}
	}	// for z

	// get the first and last slice
	for(z=0; z<sizez; z++)
	{
		if(vertebraTemplate[z].cordRadius.x!=-1)
		{
			segInfo.bound1.z = z;
			break;
		}
	}
	// last key slice
/*	for(z=sizez-1; z>0; z--)
	{
		if(vertebraTemplate[z].isKey==1)
		{
			segInfo.bound2.z = z;
			break;
		}
	}
*/
	segInfo.bound2.z = segInfo.sacrum_start;

	// interpolate the non-key slices
	int tryTimes=0;
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(vertebraTemplate[z].isKey==1) continue;
		if(vertebraTemplate[z].cordRadius.x==-1) continue;

		k=z*sizexy;

		// get an initial bounding box
		IntVec2 bb0, bb1;
		bb0 = vertebraTemplate[z].cordCenter-vertebraTemplate[z].cordRadius;
		bb1 = vertebraTemplate[z].cordCenter+vertebraTemplate[z].cordRadius;

		bb0.x-=3; bb0.y-=3; bb1.x+=3; bb1.y+=3;
		if(bb0.x<0) bb0.x=0; 
		if(bb0.y<0) bb0.y=0; 
		if(bb1.x>sizex-1) bb1.x=sizex-1; 
		if(bb1.y>sizey-1) bb1.y=sizey-1; 
		int sub_sizex, sub_sizey;
		int sk, sx, sy, nx, ny;
		int cb;	// current blob
		CIS_Array_Image2D_short *sub_binImg2D, *sub_blobImg2D;
		short *subBinA, *subBlobA;
		sub_sizex = bb1.x-bb0.x+1;
		sub_sizey = bb1.y-bb0.y+1;

		sub_binImg2D = new CIS_Array_Image2D_short(sub_sizex, sub_sizey);
		sub_blobImg2D = new CIS_Array_Image2D_short(sub_sizex, sub_sizey);
		subBinA = sub_binImg2D->GetArray();
		subBlobA = sub_blobImg2D->GetArray();

		int count=0;
		// set the initial region
		for(y=bb0.y, sk=0; y<=bb1.y; y++)
		{
			for(x=bb0.x; x<=bb1.x; x++, sk++)
			{
				k2 = y*sizex+x;
				if(maskA[k+k2]==maskStruct.mask_body) 
				{
					count++;
					subBinA[sk]=0;					
				}
				else subBinA[sk]=1;
			}
		}

		if(count<segPara.cordSizeThresh) 
		{
			delete sub_binImg2D;
			delete sub_blobImg2D;
			continue;
		}

		// remove small pieces
		NIH_Algo_Blob_Labelling(sub_binImg2D, sub_blobImg2D, blobRank2D, blobCentroid2D, false);
		for(y=0, sk=0; y<sub_sizey; y++)
		{
			for(x=0; x<sub_sizex; x++, sk++)
			{
				cb = subBlobA[sk]-1;
				if(cb<=0) continue;
				if(blobRank2D[cb]<10) subBinA[sk]=0;
			}
		}
		for(y=0, sk=0; y<sub_sizey; y++)
		{
			for(x=0; x<sub_sizex; x++, sk++)
			{
				if(subBinA[sk]==0) subBinA[sk]=1;
				else subBinA[sk]=0;
			}
		}

///		// close the hole inside
///		CIS_IPA_Dilate(sub_binImg2D, 1, false);
///		CIS_IPA_Erode(sub_binImg2D, 1, false);
		// separate cord and external space
		CIS_IPA_Erode(sub_binImg2D, segPara.closeCordIteration, true);
		CIS_IPA_Dilate(sub_binImg2D, segPara.closeCordIteration, true);

		NIH_Algo_Blob_Labelling(sub_binImg2D, sub_blobImg2D, blobRank2D, blobCentroid2D, false);

		// test the validity of the center blob
		// 1. should be the largest blob
		// 2. should be in the center
		// 3. should not have too much contact to the border
		bool valid=true;
		// check center
		cb = subBlobA[sub_sizey/2*sub_sizex+sub_sizex/2]-1;
		if(cb==-1) valid=false;
		// check size
		if(valid)
		{
			for(int b=0; b<blobRank2D.GetSize(); b++)
			{
				if(blobRank2D[b]>blobRank2D[cb])
				{
					valid =false;
					break;
				}
			}
		}
		// check border
		if(valid)
		{
			int bcount_u=0;
			int bcount_d=0;
			int bcount_l=0;
			int bcount_r=0;
			for(sy=0, sk=0; sy<sub_sizey; sy++)
			{
				for(sx=0; sx<sub_sizex; sx++, sk++)
				{
					if(subBlobA[sk]==cb+1)
					{
						if(sy<=1) bcount_u++;
						if(sy>=sub_sizey-2) bcount_d++;
						if(sx<=1) bcount_l++;
						if(sx>=sub_sizex-2) bcount_r++;
					}
				}
			}
			if(bcount_u>(sub_sizex)/2) valid=false;
			if(bcount_d>(sub_sizex)/2) valid=false;
			if(bcount_l>(sub_sizey)*3/4) valid=false;
			if(bcount_r>(sub_sizey)*3/4) valid=false;
		}

		// convert the spinal cord
		// fill the mask
		if(valid)
		{
			for(sy=0, sk=0; sy<sub_sizey; sy++)
			{
				for(sx=0; sx<sub_sizex; sx++, sk++)
				{
					if(subBlobA[sk]==cb+1)
					{
						k2 = (sy+bb0.y)*sizex+(sx+bb0.x);
						if(maskA[k2+k]==maskStruct.mask_body)
							maskA[k2+k]=maskStruct.mask_spinalCord;
					}
				}
			}
			vertebraTemplate[z].isKey=0;	// interpolated
		}	// if valid

		delete sub_binImg2D;
		delete sub_blobImg2D;

		// change the centroid of cord and try again
		if(vertebraTemplate[z].isKey<0 && tryTimes==0)
		{
			tryTimes++;
			int lKey=-1,nKey=-1;
			Vec2 lCentroid, nCentroid, iCentroid;

			for(int z1=z-1; z1>=0; z1--)
			{
				if(vertebraTemplate[z1].isKey==1)
				{
					lKey = z1; 
					lCentroid = vertebraTemplate[lKey].cordCenter;
					break;
				}
			}
			for(int z1=z+1; z1<sizez; z1++)
			{
				if(vertebraTemplate[z1].isKey==1)
				{
					nKey = z1; 
					nCentroid = vertebraTemplate[nKey].cordCenter;
					break;
				}
			}
			if(lKey!=-1 && nKey!=-1)
			{
				iCentroid.x = (nCentroid.x-lCentroid.x)/(nKey-lKey)*(z-lKey)+lCentroid.x;
				iCentroid.y = (nCentroid.y-lCentroid.y)/(nKey-lKey)*(z-lKey)+lCentroid.y;
				vertebraTemplate[z].cordCenter = iCentroid;
				z--;
			}
		}	// if vertebra
		else tryTimes=0;
	}	// for z

	// redefine range
	// get the first and last slice
	for(z=segInfo.bound1.z; z<sizez; z++)
	{
		if(vertebraTemplate[z].isKey>=0)
		{
			segInfo.bound1.z = z;
			break;
		}
	}
	// last key slice
	for(z=segInfo.bound2.z; z>0; z--)
	{
		if(vertebraTemplate[z].isKey>=0)
		{
			segInfo.bound2.z = z;
			break;
		}
	}


	// refine the cord segmentation
	for(z=segInfo.bound1.z+1; z<=segInfo.bound2.z-1; z++)
	{
		// initialize the 2D image
		for(k2=0; k2<sizexy; k2++) binArray2D[k2]=0;

		IntVec2 bb0, bb1;
		bb0 = IntVec2(sizex, sizey);
		bb1 = IntVec2(0, 0);
		
		k = z*sizexy;
		count=0;
		if(vertebraTemplate[z].isKey<0)
		{
			for(k2=0, y=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k++, k2++)
				{
					if(maskA[k]==maskStruct.mask_spinalCord ||
						(maskA[k]==maskStruct.mask_body && imgA[k]<segPara.cordIntensityThresh && 
						(maskA[k-sizexy]==maskStruct.mask_spinalCord)))
					{
						binArray2D[k2] = 1;
						count++;
						if(x<bb0.x) bb0.x=x;
						if(x>bb1.x) bb1.x=x;
						if(y<bb0.y) bb0.y=y;
						if(y>bb1.y) bb1.y=y;
					}
				}
			}	// for k2
		}
		else 
		{
			for(k2=0, y=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k++, k2++)
				{
					if(maskA[k]==maskStruct.mask_spinalCord ||
						(maskA[k]==maskStruct.mask_body && imgA[k]<segPara.cordIntensityThresh && 
						(maskA[k-sizexy]==maskStruct.mask_spinalCord && maskA[k+sizexy]==maskStruct.mask_spinalCord)))
					{
						binArray2D[k2] = 1;
						count++;
						if(x<bb0.x) bb0.x=x;
						if(x>bb1.x) bb1.x=x;
						if(y<bb0.y) bb0.y=y;
						if(y>bb1.y) bb1.y=y;
					}
				}
			}	// for k2
		}

		if(count<segPara.cordSizeThresh) continue;

		bb0.x-=3; bb0.y-=3; bb1.x+=3; bb1.y+=3;
		if(bb0.x<0) bb0.x=0; 
		if(bb0.y<0) bb0.y=0; 
		if(bb1.x>sizex-1) bb1.x=sizex-1; 
		if(bb1.y>sizey-1) bb1.y=sizey-1; 

		int sub_sizex, sub_sizey;
		int sk, sx, sy, nx, ny;
		CIS_Array_Image2D_short *sub_binImg2D, *sub_blobImg2D;

		sub_sizex = bb1.x-bb0.x+1;
		sub_sizey = bb1.y-bb0.y+1;

		sub_binImg2D = new CIS_Array_Image2D_short(sub_sizex, sub_sizey);
		sub_blobImg2D = new CIS_Array_Image2D_short(sub_sizex, sub_sizey);
		short *subBinA = sub_binImg2D->GetArray();
		short *subBlobA = sub_blobImg2D->GetArray();

		for(sy=0, sk=0; sy<sub_sizey; sy++)
			for(sx=0, k2=(sy+bb0.y)*sizex+bb0.x; sx<sub_sizex; sx++, sk++, k2++)
				subBinA[sk]=binArray2D[k2];

		// close the hole inside
		CIS_IPA_Dilate(sub_binImg2D, 1, false);
		CIS_IPA_Erode(sub_binImg2D, 1, false);
		// separate small bridges
		CIS_IPA_Erode(sub_binImg2D, segPara.closeCordIteration, true);
		CIS_IPA_Dilate(sub_binImg2D, segPara.closeCordIteration, true);

		// remove bone structure
		k = z*sizexy;
		for(sy=0, sk=0; sy<sub_sizey; sy++)
		{
			for(sx=0, k2=(sy+bb0.y)*sizex+bb0.x; sx<sub_sizex; sx++, sk++, k2++)
			{
				if(subBinA[sk]==1)
				{
					if(maskA[k+k2]!=maskStruct.mask_body && maskA[k+k2]!=maskStruct.mask_spinalCord)
						subBinA[sk]=0;
				}
			}
		}


		NIH_Algo_Blob_Labelling(sub_binImg2D, sub_blobImg2D, blobRank2D, blobCentroid2D, false);

		// closest blob near the center
		int largest_blob=-1;
		int largest_size=0;
		int cb;
		for(cb=0;cb<blobRank2D.GetSize();cb++)
		{
			if(blobRank2D[cb]>largest_size)
			{
				largest_size=blobRank2D[cb];
				largest_blob=cb;
			}
		}

		int closest_blob=-1;
		int closest_dist=sizex+sizey;

		for(cb=0;cb<blobRank2D.GetSize();cb++)
		{
			int dist=abs(blobCentroid2D[cb].x-vertebraTemplate[z].cordCenter.x+bb0.x)+
				abs(blobCentroid2D[cb].y-vertebraTemplate[z].cordCenter.y+bb0.y);
			if(dist<closest_dist)
			{
				closest_dist=dist;
				closest_blob=cb;
			}
		}

		if(largest_blob==-1 || closest_blob==-1)
		{
			delete sub_binImg2D;
			delete sub_blobImg2D;
			continue;
		}

		// choose between largest_blob and closest_blob
		int selected_blob;
		if(largest_blob==closest_blob) selected_blob=largest_blob;
		else if(blobRank2D[closest_blob]>blobRank2D[largest_blob]/4) selected_blob=closest_blob;
		else selected_blob=largest_blob;

		k = z*sizexy;
		for(sy=0, sk=0; sy<sub_sizey; sy++)
		{
			for(sx=0, k2=(sy+bb0.y)*sizex+bb0.x; sx<sub_sizex; sx++, sk++, k2++)
			{
				if(subBlobA[sk]==selected_blob+1)
				{
					if(maskA[k+k2]==maskStruct.mask_body) maskA[k+k2]=maskStruct.mask_spinalCord;
				}
			}
		}
		vertebraTemplate[z].cordCenter = blobCentroid2D[selected_blob]+bb0;

		if(vertebraTemplate[z].isKey<0) vertebraTemplate[z].isKey=0;
		delete sub_binImg2D;
		delete sub_blobImg2D;
	}	// for z

	// dilate the cord to fill the gap
	IntVec2 seed;
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		k = z*sizexy;
		count = 0;
		for(k2=0; k2<sizexy; k2++, k++) 
		{
			if(maskA[k]==maskStruct.mask_spinalCord) 
			{
				count++;
				binArray2D[k2]=1;
			}
			else binArray2D[k2]=0;
		}

		if(count==0) continue;

		CIS_IPA_Dilate(binImg2D, segPara.closeCordIteration, true);

		k=z*sizexy;
		for(k2=0; k2<sizexy; k2++, k++) 
		{
			if(binArray2D[k2]==1 && maskA[k]==maskStruct.mask_body) 
				maskA[k]=maskStruct.mask_spinalCord;
		}

		// only keep the largest piece of cord
		k=z*sizexy;
		for(y=0, k2=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k2++, k++)
			{
				if(maskA[k]==maskStruct.mask_spinalCord)
				{
					binArray2D[k2]=1;
				}
				else binArray2D[k2]=0;
			}
		}
		NIH_Algo_Blob_Labelling(binImg2D, blobImg2D, blobRank2D, blobCentroid2D, false);
		if(blobRank2D.GetSize()>1)
		{
			int largestRank, largestRankI=1;
			largestRank=0;
			for(k=0; k<blobRank2D.GetSize(); k++)
			{
				if(blobRank2D[k]>largestRank)
				{
					largestRank = blobRank2D[k];
					largestRankI = k+1;
				}
			}
			vertebraTemplate[z].cordCenter = blobCentroid2D[largestRankI-1];
			k = z*sizexy;
			for(k2=0, y=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k++, k2++)
				{
					if(binArray2D[k2]==1) 
					{
						if(blobArray2D[k2]!=largestRankI)
						{
							maskA[k]=maskStruct.mask_body;
						}
					}	// if binArray
				}
			}
		}
		
		// find the bounding box of cord
		k=z*sizexy;
		vertebraTemplate[z].cord_bb0 = IntVec2(sizex, sizey);
		vertebraTemplate[z].cord_bb1 = IntVec2(0,0);
		for(y=0, k2=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k2++, k++)
			{
				if(maskA[k]==maskStruct.mask_spinalCord)
				{
					if(x<vertebraTemplate[z].cord_bb0.x) vertebraTemplate[z].cord_bb0.x = x;
					if(x>vertebraTemplate[z].cord_bb1.x) vertebraTemplate[z].cord_bb1.x = x;
					if(y<vertebraTemplate[z].cord_bb0.y) vertebraTemplate[z].cord_bb0.y = y;
					if(y>vertebraTemplate[z].cord_bb1.y) vertebraTemplate[z].cord_bb1.y = y;
				}
			}
		}
		vertebraTemplate[z].cordRadius = (vertebraTemplate[z].cord_bb1-vertebraTemplate[z].cord_bb0)/2;
	}	// for z

	// interpolate the contour in sagittal view

	// fill the interior holes with spongy bones
	for(z=0; z<sizez; z++)
	{
		// remove the slices with no spinal cord, why??
		if(vertebraTemplate[z].isKey<0)
		{
			k = z*sizexy;
			for(k2=0, y=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k++, k2++)
				{
					if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_spongyBone)
					{
						maskA[k] = maskStruct.mask_otherBone;
					}
				}
			}

			continue;	// skip to next slice
		}

		IntVec2 bb0, bb1;
		// first define a bounding box
		k = z*sizexy;
		count =0;
		bb0.x = sizex-3;
		bb0.y = sizey-3;
		bb1.x = 3;
		bb1.y = 3;
		for(k2=0, y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++, k2++)
			{
				if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_spongyBone)
				{
					if(x<bb0.x) bb0.x=x;
					if(x>bb1.x) bb1.x=x;
					if(y<bb0.y) bb0.y=y;
					if(y>bb1.y) bb1.y=y;
					count++;
				}
			}
		}
		if(count==0) continue;

		bb0.x-=6; bb0.y-=6; bb1.x+=6; bb1.y+=6;
		if(bb0.x<2) bb0.x=2;
		if(bb0.y<2) bb0.y=2;
		if(bb1.x>sizex-3) bb1.x=sizex-3;
		if(bb1.y>sizey-3) bb1.y=sizey-3;

		int sub_sizex, sub_sizey;
		int sk, sx, sy, nx, ny;
		CIS_Array_Image2D_short *sub_binImg2D, *sub_blobImg2D;

		sub_sizex = bb1.x-bb0.x+1;
		sub_sizey = bb1.y-bb0.y+1;

		sub_binImg2D = new CIS_Array_Image2D_short(sub_sizex, sub_sizey);
		sub_blobImg2D = new CIS_Array_Image2D_short(sub_sizex, sub_sizey);
		short *subBinA = sub_binImg2D->GetArray();
		short *subBlobA = sub_blobImg2D->GetArray();


		k = z*sizexy;
		count=0;
		for(y=bb0.y, sk=0; y<=bb1.y; y++)
		{
			k2 = y*sizex+bb0.x;
			for(x=bb0.x; x<=bb1.x; x++, k2++, sk++)
			{
				if(maskA[k+k2]==maskStruct.mask_body || maskA[k+k2]==maskStruct.mask_air)
				{
					subBinA[sk]=1;
					count ++;
				}
				else subBinA[sk]=0;
			}
		}	// for y

		if(count==0) continue;

		// separate interior and external space
		CIS_IPA_Erode(sub_binImg2D, segPara.closeCordIteration, true);
		CIS_IPA_Dilate(sub_binImg2D, segPara.closeCordIteration, true);

		NIH_Algo_Blob_Labelling(sub_binImg2D, sub_blobImg2D, blobRank2D, blobCentroid2D, false);

		blobStatus.SetSize(blobRank2D.GetSize());
		int cb;
		for(cb=0; cb<blobStatus.GetSize(); cb++) blobStatus[cb]=1;

		// Get the outer space, any blob next to border are considered outer space
		k = z*sizexy;
		for(y=bb0.y, sk=0; y<=bb1.y; y++)
		{
			k2 = y*sizex+bb0.x;
			for(x=bb0.x; x<=bb1.x; x++, k2++, sk++)
			{
				if(subBinA[sk]==1)
				{
					cb = subBlobA[sk]-1;
					if(cb<0 || blobStatus[cb]==-2) continue;
					if(y<bb0.y+segPara.closeCordIteration || y>=bb1.y-segPara.closeCordIteration
						|| x<bb0.x+segPara.closeCordIteration || x>=bb1.x-segPara.closeCordIteration)
					{
						blobStatus[cb] = -2;
					}
				}
			}
		}	// for y

		// assign mask
		k = z*sizexy;
		for(y=bb0.y, sk=0; y<=bb1.y; y++)
		{
			k2 = y*sizex+bb0.x;
			for(x=bb0.x; x<=bb1.x; x++, k2++, sk++)
			{
//				if(subBinA[sk]==0 && maskA[k+k2]==maskStruct.mask_body) maskA[k+k2]=maskStruct.mask_spongyBone;
//				else 
					if(subBinA[sk]==1)
				{
					cb = subBlobA[sk]-1;
					if(cb>=0 && blobStatus[cb]>=0)
					{
						maskA[k+k2]=maskStruct.mask_spongyBone;
						// also fill its neighbors
						if(maskA[k+k2+1]==maskStruct.mask_body) maskA[k+k2+1]==maskStruct.mask_spongyBone;
						if(maskA[k+k2-1]==maskStruct.mask_body) maskA[k+k2-1]==maskStruct.mask_spongyBone;
						if(maskA[k+k2+sizex]==maskStruct.mask_body) maskA[k+k2+sizex]==maskStruct.mask_spongyBone;
						if(maskA[k+k2-sizex]==maskStruct.mask_body) maskA[k+k2-sizex]==maskStruct.mask_spongyBone;
					}
				}	// if binArray
			}
		}

		delete sub_binImg2D;
		delete sub_blobImg2D;

	}	// for z


	// segment other structures


	int cordLength = segInfo.bound2.z-segInfo.bound1.z+1;

	// compute the center of the spinal cannel and radius of cannel
	// and find the bounding box for each slice
	IntVec2 c_center, c_bound1, c_bound2;
	int len_y1, len_y2, len_x1, len_x2, largest_len;
	
	segInfo.spineBound1.SetSize(sizez);
	segInfo.spineBound2.SetSize(sizez);
	segInfo.diskBound1.SetSize(sizez);
	segInfo.diskBound2.SetSize(sizez);
	segInfo.cordBound1.SetSize(sizez);
	segInfo.cordBound2.SetSize(sizez);
	segInfo.sprocessBound1.SetSize(sizez);
	segInfo.sprocessBound2.SetSize(sizez);

	segInfo.cordCenter.SetSize(sizez);
	segInfo.cordRadius.SetSize(sizez);

	segInfo.diskCenter.SetSize(sizez);
	segInfo.diskRadius.SetSize(sizez);

	for(z=0; z<sizez; z++)
	{
		segInfo.spineBound1[z] = IntVec2(-1,-1);
		segInfo.spineBound2[z] = IntVec2(-1,-1);
		segInfo.diskBound1[z] = IntVec2(-1,-1);
		segInfo.diskBound2[z] = IntVec2(-1,-1);
		segInfo.cordBound1[z] = IntVec2(-1,-1);
		segInfo.cordBound2[z] = IntVec2(-1,-1);
		segInfo.sprocessBound1[z] = IntVec2(-1,-1);
		segInfo.sprocessBound2[z] = IntVec2(-1,-1);

		segInfo.cordCenter[z] = IntVec2(-1,-1);
		segInfo.cordRadius[z] = IntVec2(-1,-1);
		segInfo.diskCenter[z] = IntVec2(-1,-1);
		segInfo.diskRadius[z] = IntVec2(-1,-1);
	}

	for(z=0; z<sizez; z++) yt[z]=-1;

	// second interpolation
	// interpolate x coordinate of cord using key slices
	minInterpolate = -1;
	maxInterpolate = -1;
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(vertebraTemplate[z].isKey==1)
		{
			yt[z] = vertebraTemplate[z].cordCenter.x;
			if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
			if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
		}
		else yt[z]=-1;
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		segInfo.cordCenter[z].x = vertebraTemplate[z].cordCenter.x;
		if(yt2[z]==-1) continue;
		if(yt2[z]<minInterpolate) yt2[z]=minInterpolate;
		if(yt2[z]>maxInterpolate) yt2[z]=maxInterpolate;
		if(vertebraTemplate[z].isKey!=1) segInfo.cordCenter[z].x = (int)yt2[z];
	}

	// interpolate y coordinate of cord
	minInterpolate = -1;
	maxInterpolate = -1;
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(vertebraTemplate[z].isKey==1)
		{
			yt[z] = vertebraTemplate[z].cordCenter.y;
			if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
			if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
		}
		else yt[z]=-1;
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		segInfo.cordCenter[z].y = vertebraTemplate[z].cordCenter.y;
		if(yt2[z]==-1) continue;
		if(yt2[z]<minInterpolate) yt2[z]=minInterpolate;
		if(yt2[z]>maxInterpolate) yt2[z]=maxInterpolate;
		if(vertebraTemplate[z].isKey!=1) segInfo.cordCenter[z].y = (int)yt2[z];
	}


	// interpolate x radius of cord
	minInterpolate = -1;
	maxInterpolate = -1;
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(vertebraTemplate[z].isKey==1)
		{
			yt[z] = vertebraTemplate[z].cordRadius.x;
			if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
			if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
		}
		else yt[z]=-1;
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		segInfo.cordRadius[z].x = vertebraTemplate[z].cordRadius.x;
		if(yt2[z]==-1) continue;
		if(yt2[z]<minInterpolate) yt2[z]=minInterpolate;
		if(yt2[z]>maxInterpolate) yt2[z]=maxInterpolate;
		if(vertebraTemplate[z].isKey!=1) segInfo.cordRadius[z].x = (int)yt2[z];
	}

	// interpolate y radius of cord
	minInterpolate = -1;
	maxInterpolate = -1;
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(vertebraTemplate[z].isKey==1)
		{
			yt[z] = vertebraTemplate[z].cordRadius.y;
			if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
			if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
		}
		else yt[z]=-1;
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		segInfo.cordRadius[z].y = vertebraTemplate[z].cordRadius.y;
		if(yt2[z]==-1) continue;
		if(yt2[z]<minInterpolate) yt2[z]=minInterpolate;
		if(yt2[z]>maxInterpolate) yt2[z]=maxInterpolate;
		if(vertebraTemplate[z].isKey!=1) segInfo.cordRadius[z].y = (int)yt2[z];
	}

	// compute the center of the disk, and radius of disk
	// interpolate the radius of the disk

	for(z=0; z<sizez; z++)
	{
		xt[z] = z;
		yt[z] = -1;
	}
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		// initial value of diskRaius
		segInfo.diskRadius[z].y = -1;

		c_center.x = segInfo.cordCenter[z].x;
		c_center.y = segInfo.cordCenter[z].y;

		// seaching from the center of cord to the top of vertebra
		int y_top, y_count;
		y_top = y_count = 0;
		k = z*sizexy;
		for(x=c_center.x-2; x<=c_center.x+2; x++)
		{
			for(y=c_center.y-segInfo.cordRadius[z].y*8; y<c_center.y-segInfo.cordRadius[z].y*3; y++)
			{
				k2 = k+y*sizex+x;
				if(maskA[k2]==maskStruct.mask_spongyBone || maskA[k2]==maskStruct.mask_corticalBone)
				{
					y_count++;
					y_top += y;
					break;
				}
			}
		}

		// estimate the diskRadius
		if(y_count>0)
		{
			y_top/=y_count;
			if(y_top<c_center.y) 
			{
				yt[z] = (c_center.y-y_top)/2;
				segInfo.diskRadius[z].y = yt[z];
			}
		}
	}	// for z

	// refine the radius of the disk, 
	// computing the radius of top half of the circle
	// then average it with the current radius
	//
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		Vec2 dir, curPos, nextPos, maxPos, maxRadius;
		int search_r1, search_r0, r, k3, grad, maxGrad, max_r;
		int angle, angle1;

		c_center.x = segInfo.cordCenter[z].x;
		c_center.y = segInfo.cordCenter[z].y-segInfo.diskRadius[z].y;

		int angleLeft, angleRight;
		angleLeft = 270;
		angleRight = 90;

		int count=0;
		int accu_radius=0;

		search_r0 = 1;
		search_r1 = segInfo.diskRadius[z].y+segInfo.cordRadius[z].y;

		for(angle=angleRight; angle<=angleLeft; angle+=diskAngleInterval)
		{
			// top 
			dir.x = sin(angle*3.14/180);
			dir.y = cos(angle*3.14/180);


			curPos.x = c_center.x;
			curPos.y = c_center.y;
			k = z*sizexy;
			maxGrad = 0;
			curPos = curPos+dir*search_r0;
			for(r=search_r0; r<=search_r1; r++)
			{
				nextPos = curPos+dir;
				k2 = k + (int)(curPos.y+0.5)*sizex + (int)(curPos.x+0.5);
				k3 = k + (int)(nextPos.y+0.5)*sizex + (int)(nextPos.x+0.5);

				if(imgA[k2]>segPara.boneThresh && imgA[k3]<segPara.boneThresh && imgA[k3]>1000)
				{
					grad = imgA[k2]-imgA[k3];
	
					if(grad>maxGrad && grad>50)
					{
						maxGrad=grad;
						max_r = r;
						maxPos = curPos;
					}
				}

				curPos += dir;
			}

			if(maxGrad>0)
			{
				accu_radius += max_r;
				count++;
			}
		}	// for angle

		if(count>0)
		{
			accu_radius /= count;
			segInfo.diskRadius[z].y = (segInfo.diskRadius[z].y+accu_radius)/2;
		}
		
	}	// for z
	
	// constrain the disk radius of the first and last few slices to avoid wild extrapolation
	// make sure it is greater than cord radius and smaller than 4 times of it
	minInterpolate = -1;
	maxInterpolate = -1;
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		yt[z] = segInfo.diskRadius[z].y;
		if(yt[z]<segInfo.cordRadius[z].y*1.5 || yt[z]>segInfo.cordRadius[z].y*6)
		{
			yt[z] = -1;
			continue;	// skip if too small or too big
		}
		if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
		if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
	}

	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(yt2[z]<minInterpolate) yt2[z]=minInterpolate;
		if(yt2[z]>maxInterpolate) yt2[z]=maxInterpolate;

		segInfo.diskCenter[z].x = segInfo.cordCenter[z].x;

		if(yt2[z]>0)
		{
			segInfo.diskRadius[z].y = (int)yt2[z];

		}
		else
		{
			// if not available, use the closest one
			if(z>1) segInfo.diskRadius[z].y=segInfo.diskRadius[z-1].y;
			else segInfo.diskRadius[z].y=2*segInfo.cordRadius[z].y;
		}
		segInfo.diskCenter[z].y = segInfo.cordCenter[z].y-(int)yt2[z];

		// visalize it 
		if(debugMode && false)
		{
			maskA[z*sizexy+(int)segInfo.cordCenter[z].y*sizex+segInfo.cordCenter[z].x] = maskStruct.mask_vertebralDisk;
			maskA[z*sizexy+(int)segInfo.diskCenter[z].y*sizex+segInfo.diskCenter[z].x] = maskStruct.mask_vertebra;
		}	
	}

	// locate and refine the vertebra template
	// it contains three consequential steps, 
	// 1. the disk
	// 2. spinal process
	// 3. left and right pedicle
	//
	
	// interpolate the contour of disk

	// based on center of the disk, search for a few points along the boundary
	// the disk is modelled as a circle
	//

	for(z=0; z<sizez; z++)
	{
		for(int a=0; a<diskAngleCount; a++)
		{
			vertebraTemplate[z].diskContour[a] = Vec2(-1, -1);
			vertebraTemplate[z].diskContourRadius[a] = -1;
			vertebraTemplate[z].interpolatedRadius[a] = -1;
		}
		vertebraTemplate[z].diskNeckLeft = Vec2(-1, -1);
		vertebraTemplate[z].diskNeckRight = Vec2(-1, -1);
	}
	
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		Vec2 dir, curPos, nextPos, maxPos, maxRadius;
		int search_r1, search_r0, r, k3, grad, maxGrad, max_r;
		int angle, angle1;

		// locate left and right border of vertebra neck
		// search from the center of spinal cord
		c_center.x = segInfo.cordCenter[z].x;
		c_center.y = segInfo.cordCenter[z].y;

		y = c_center.y;
		// left side
		k = z*sizexy;
		k2 = k+y*sizey+c_center.x-segInfo.cordRadius[z].x;
		maxGrad = 0;
		for(x=c_center.x-segInfo.cordRadius[z].x; x>=c_center.x-segInfo.cordRadius[z].x*3; x--, k2--)
		{
			if(maskA[k2]==maskStruct.mask_spinalCord)
			{
				maxGrad = 1;
				maxPos = Vec2(x,y);
			}
			else if(maskA[k2]==maskStruct.mask_corticalBone && maskA[k2-1]!=maskStruct.mask_corticalBone 
				&& imgA[k2-1]>segPara.softThresh)
			{
				grad = imgA[k2]-imgA[k2-1];

				if(grad>maxGrad && grad>segPara.borderMinGrad)	
				{
					maxGrad = grad;
					maxPos = Vec2(x,y);
				}
			}
		}
		if(maxGrad>0)
		{
			vertebraTemplate[z].diskNeckLeft = maxPos;
		}
		else
		{
			vertebraTemplate[z].diskNeckLeft = Vec2(c_center.x-segInfo.cordRadius[z].x,y);
		}


		// search the right side
		k = z*sizexy;
		k2 = k+y*sizey+c_center.x+segInfo.cordRadius[z].x;
		maxGrad = 0;
		for(x=c_center.x+segInfo.cordRadius[z].x; x<=c_center.x+segInfo.cordRadius[z].x*3; x++, k2++)
		{
			if(maskA[k2]==maskStruct.mask_spinalCord)
			{
				maxGrad = 1;
				maxPos = Vec2(x,y);
			}
			else if(maskA[k2]==maskStruct.mask_corticalBone && maskA[k2+1]!=maskStruct.mask_corticalBone 
				&& imgA[k2+1]>segPara.softThresh)
			{
				grad = imgA[k2]-imgA[k2+1];

				if(grad>maxGrad && grad>segPara.borderMinGrad)
				{
					maxGrad = grad;
					maxPos = Vec2(x,y);
				}
			}
		}
		if(maxGrad>0)
		{
			vertebraTemplate[z].diskNeckRight = maxPos;
		}
		else
		{
			vertebraTemplate[z].diskNeckRight = Vec2(c_center.x+segInfo.cordRadius[z].x,y);
		}

		// raise the neck to the top of spinal cord
///		vertebraTemplate[z].diskNeckLeft.y -= segInfo.cordRadius[z].y;
///		vertebraTemplate[z].diskNeckRight.y -= segInfo.cordRadius[z].y;

		c_center.x = segInfo.diskCenter[z].x;
		c_center.y = segInfo.diskCenter[z].y;

		// define the range of angle's of the disk
		int angleLeft, angleRight;
		angleLeft = (int)(atan2(vertebraTemplate[z].diskNeckLeft.x-c_center.x, vertebraTemplate[z].diskNeckLeft.y-c_center.y)*180/3.14)+360;
		angleRight = (int)(atan2(vertebraTemplate[z].diskNeckRight.x-c_center.x, vertebraTemplate[z].diskNeckRight.y-c_center.y)*180/3.14);

		angleLeft = (angleLeft/diskAngleInterval+1)*diskAngleInterval;
		angleRight = (angleRight/diskAngleInterval-1)*diskAngleInterval;
		vertebraTemplate[z].angleLeft = angleLeft;
		vertebraTemplate[z].angleRight = angleRight;
		
		// search for the boundary of the disk along each shooting angle
		// search along the shooting ray and get the one with maximum gradian
		//
		for(angle=angleRight; angle<=angleLeft; angle+=diskAngleInterval)
		{
			angle1 = angle/diskAngleInterval;

			// top 
			dir.x = sin(angle*3.14/180);
			dir.y = cos(angle*3.14/180);

			search_r0 = segInfo.diskRadius[z].y-segInfo.cordRadius[z].y;
			search_r1 = segInfo.diskRadius[z].y+segInfo.cordRadius[z].y;
			x = c_center.x;

			curPos.x = c_center.x;
			curPos.y = c_center.y;
			k = z*sizexy;
			maxGrad = 0;
			curPos = curPos+dir*search_r0;
			for(r=search_r0; r<=search_r1; r++)
			{
				if(curPos.y>=segInfo.cordCenter[z].y) break;

				nextPos = curPos+dir;
				k2 = k + (int)(curPos.y+0.5)*sizex + (int)(curPos.x+0.5);
				k3 = k + (int)(nextPos.y+0.5)*sizex + (int)(nextPos.x+0.5);
				if(maskA[k2]==maskStruct.mask_spinalCord)
				{
					if(maxGrad==0)
					{
						maxGrad==1;
						max_r = r;
						maxPos = curPos;
					}
					break;
				}

				if(imgA[k2]>segPara.boneThresh && imgA[k3]<segPara.boneThresh && imgA[k3]>segPara.softThresh)
				{
					grad = imgA[k2]-imgA[k3];
	
					if(grad>maxGrad && grad>segPara.borderMinGrad) 
					{
						maxGrad=grad;
						max_r = r;
						maxPos = curPos;
					}
				}

				curPos += dir;
			}

			if(maxGrad>0)
			{
				vertebraTemplate[z].diskContour[angle1] = maxPos;
				vertebraTemplate[z].diskContourRadius[angle1] = max_r;
			}
		}	// for angle
		
	}	// for z

	// interpolate every point on the contour
	// interpolate the radius
	// do the interpolation twice to refine the result
	//
	doubleDynArray yt2_x, yt2_y;
	int angle;
	yt2_x.SetSize(sizez);
	yt2_y.SetSize(sizez);
	for(angle=0; angle<360; angle+=diskAngleInterval)
	{
		int angle1 = angle/diskAngleInterval;

		int count;
		for(z=0; z<sizez; z++)
		{
			yt[z] = -1;
		}

		// interpolate the radius at each angle on the disk
		count = 0;
		minInterpolate = -1;
		maxInterpolate = -1;
		for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
		{
			if(vertebraTemplate[z].diskContourRadius[angle1]>segInfo.cordRadius[z].y)
			{
				count++;
				yt[z] = vertebraTemplate[z].diskContourRadius[angle1];
				if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
				if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
			}
		}

		if(count<cordLength/2) continue;

		PiecewiseBernsteinSmoothing(xt, yt, yt2_x, bernstein_power, piece_size);

		Vec2 dir;
		dir.x = sin(angle*3.14/180);
		dir.y = cos(angle*3.14/180);

		for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
		{
			if(yt2_x[z]==-1) continue;
			if(yt2_x[z]<minInterpolate) yt2_x[z]=minInterpolate;
			if(yt2_x[z]>maxInterpolate) yt2_x[z]=maxInterpolate;
			vertebraTemplate[z].interpolatedRadius[angle1] = yt2_x[z];
		}

	}	// for angle

	// quality assurance
	// remove outliers from each slices
	//
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		double avg_r, dev_r, r1, r2,r3;
		double count_r;
		int angle, angle1;
		Vec2 dir, curPos, nextPos, maxPos, maxRadius;
		int search_r1, search_r0, r, k3, grad, maxGrad, max_r;

		avg_r = 0;
		dev_r = 0;
		count_r = 0;

		// compute the average radius
		for(angle=0; angle<360; angle+=diskAngleInterval)
		{
			angle1 = angle/diskAngleInterval;

			r1 = -1;
			if(vertebraTemplate[z].diskContourRadius[angle1]>segInfo.cordRadius[z].y)
			{
				r1 = vertebraTemplate[z].diskContourRadius[angle1];
			}
			else if(vertebraTemplate[z].interpolatedRadius[angle1]>segInfo.cordRadius[z].y)
			{
				r1 = vertebraTemplate[z].interpolatedRadius[angle1];
			}

			if(r1>segInfo.cordRadius[z].y && r1<segInfo.diskRadius[z].y*5 && r1>segInfo.diskRadius[z].y/5)
			{
				avg_r += r1;
				count_r ++;
			}
		}

		if(count_r!=0) avg_r /= count_r;

		// compute the deviation
		for(angle=0; angle<360; angle+=diskAngleInterval)
		{
			angle1 = angle/diskAngleInterval;

			r1 = -1;
			if(vertebraTemplate[z].diskContourRadius[angle1]>segInfo.cordRadius[z].y)
			{
				r1 = vertebraTemplate[z].diskContourRadius[angle1];
			}
			else if(vertebraTemplate[z].interpolatedRadius[angle1]>segInfo.cordRadius[z].y)
			{
				r1 = vertebraTemplate[z].interpolatedRadius[angle1];
			}

			if(r1>segInfo.cordRadius[z].y && r1<segInfo.diskRadius[z].y*5 && r1>segInfo.diskRadius[z].y/5)
			{
				dev_r += (r1-avg_r)*(r1-avg_r);
				count_r ++;
			}

		}

		if(count_r>1) dev_r = sqrt(dev_r/(count_r-1));

		// pick out the outliers
		for(angle=0; angle<360; angle+=diskAngleInterval)
		{
			angle1 = angle/diskAngleInterval;
			r1 = vertebraTemplate[z].diskContourRadius[angle1];
			r2 = vertebraTemplate[z].interpolatedRadius[angle1];

			if(r1>segInfo.cordRadius[z].y)
			{
				if(fabs(r1-avg_r)<dev_r*2 && r1<segInfo.diskRadius[z].y*5 && r1>segInfo.diskRadius[z].y/5) continue;		// good cases
				else if(r2>segInfo.cordRadius[z].y && fabs(r2-avg_r)<dev_r*2 && r2<segInfo.diskRadius[z].y*5 && r2>segInfo.diskRadius[z].y/5)
				{
					vertebraTemplate[z].diskContourRadius[angle1]=r2;
				}
				else
				{
					vertebraTemplate[z].diskContourRadius[angle1]=-1;
				}
			}
			else if(r2>segInfo.cordRadius[z].y && fabs(r2-avg_r)<dev_r*2 && r2<segInfo.diskRadius[z].y*5 && r2>segInfo.diskRadius[z].y/5)
			{
				vertebraTemplate[z].diskContourRadius[angle1]=r2;				
			}
		}

		// smoothing the radius along the border on each slice using b-spline interpolation
		//
		int angleNumber=360/diskAngleInterval;
		xt.SetSize(angleNumber);
		yt.SetSize(angleNumber);
		yt2.SetSize(angleNumber);
		minInterpolate = -1;
		maxInterpolate = -1;
		for(angle1=0; angle1<angleNumber; angle1++)
		{
			xt[(int)angle1] = (int)angle1;
			yt[(int)angle1] = vertebraTemplate[z].diskContourRadius[angle1];
			if(yt[(int)angle1]==-1) continue;
			if(yt[(int)angle1]<minInterpolate || minInterpolate==-1) minInterpolate=yt[(int)angle1];
			if(yt[(int)angle1]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[(int)angle1];
		}

		BernsteinSmoothing(xt, yt, yt2, bernstein_power);
		for(angle=0; angle<360; angle+=diskAngleInterval)
		{
			angle1 = angle/diskAngleInterval;

			if(angle<vertebraTemplate[z].angleRight || angle>vertebraTemplate[z].angleLeft) continue;
			if(yt2[(int)angle1]==-1) continue;
			if(yt2[(int)angle1]<minInterpolate) yt2[(int)angle1]=minInterpolate;
			if(yt2[(int)angle1]>maxInterpolate) yt2[(int)angle1]=maxInterpolate;
			vertebraTemplate[z].diskContourRadius[angle1] = yt2[angle1];
		}

		k = z*sizexy;
		// local refinement and set contour based on radius
		for(angle=0; angle<360; angle+=diskAngleInterval)
		{
			dir.x = sin(angle*3.14/180);
			dir.y = cos(angle*3.14/180);
		
			angle1 = angle/diskAngleInterval;
			r1 = vertebraTemplate[z].diskContourRadius[angle1];
			if(r1<0 && angle>vertebraTemplate[z].angleRight && angle<vertebraTemplate[z].angleLeft) r1=segInfo.diskRadius[z].y;
			if(r1>segInfo.cordRadius[z].y)
			{
				vertebraTemplate[z].diskContour[angle1].x = -1;
				vertebraTemplate[z].diskContour[angle1].y = -1;
				// do a local search to find the local maximum gradient
				//
				search_r0 = r1-2;
				search_r1 = r1+4;
				curPos = segInfo.diskCenter[z]+dir*search_r0;
				maxGrad = 0;
				for(r=search_r0; r<=search_r1; r++)
				{
					nextPos = curPos+dir;
					if(nextPos.y<0 || curPos.y<0) break;
					k2 = k + (int)(curPos.y+0.5)*sizex + (int)(curPos.x+0.5);
					k3 = k + (int)(nextPos.y+0.5)*sizex + (int)(nextPos.x+0.5);
					if(imgA[k3]>segPara.softThresh)
					{
						grad = imgA[k2]-imgA[k3];
						if(grad>maxGrad && grad>segPara.borderMinGrad)
						{
							maxGrad=grad;
							max_r = r;
							maxPos = curPos;
						}
					}

					curPos = nextPos;
				}	// for r
		
				if(maxGrad>0)
				{
					r1 = max_r;
				}

				if(r1<segInfo.diskRadius[z].y*5 && r1>segInfo.diskRadius[z].y/5)
				{
					vertebraTemplate[z].diskContourRadius[angle1] = r1;
				}
			}	// if r1
		}	// for angle

		// set contour
		for(angle=0; angle<360; angle+=diskAngleInterval)
		{
			dir.x = sin(angle*3.14/180);
			dir.y = cos(angle*3.14/180);
			angle1 = angle/diskAngleInterval;
			r1 = vertebraTemplate[z].diskContourRadius[angle1];
			if(r1>segInfo.cordRadius[z].y)
			{
				vertebraTemplate[z].diskContour[angle1].x = segInfo.diskCenter[z].x + dir.x*r1;
				vertebraTemplate[z].diskContour[angle1].y = segInfo.diskCenter[z].y + dir.y*r1;
				
				// visualize it
				if(debugMode && false)
				{
					maskA[z*sizexy+(int)vertebraTemplate[z].diskContour[angle1].y*sizex+(int)vertebraTemplate[z].diskContour[angle1].x] = 
						maskStruct.mask_falseDetection;
				}
			}
		}
	}	// for z
	
	// convert points inside the vertebra body to different value to prevent self intersecting
	//
	vec2DynArray diskROI;
	intVec2DynArray scan;
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		k = z*sizexy;
		int angleNumber=360/diskAngleInterval;
		diskROI.SetSize(0);

		for(angle=0; angle<angleNumber; angle++)
		{
			if(vertebraTemplate[z].diskContour[angle].y>0 && vertebraTemplate[z].diskContour[angle].x>0
				&& vertebraTemplate[z].diskContour[angle].y<sizey && vertebraTemplate[z].diskContour[angle].x<sizex)
			{
				diskROI.Add(vertebraTemplate[z].diskContour[angle]);
			}
		}

		CIS_Algo_Contour_GetScanWhole(diskROI, scan);

		for(i=0; i<scan.GetSize(); i++)
		{
			x = scan[i].x;
			y = scan[i].y;
			k2 = k+y*sizex+x;

			// base 100 for the disk
			if(maskA[k2]==maskStruct.mask_spinalCord)  maskA[k2]=100+maskStruct.mask_spinalCord;
			else if(imgA[k2]>segPara.boneThresh) maskA[k2]=100+maskStruct.mask_corticalBone;
			else maskA[k2]=100+maskStruct.mask_spongyBone;
		}
	}	// for z


	// locate and interpolated spinal process
	xt.SetSize(sizez);
	yt.SetSize(sizez);
	yt2.SetSize(sizez);
	for(z=0; z<sizez; z++)
	{
		xt[z] = z;
		yt[z] = -1;
	}

	for(z=0; z<sizez; z++)
	{
		vertebraTemplate[z].sprocessEnd = Vec2(-1, -1);
		vertebraTemplate[z].sprocessLength = vertebraTemplate[z].sprocessOrient = -1;
		
		for(i=0; i<sprocessCount; i++) 
		{
			vertebraTemplate[z].sprocessMedial[i] = Vec2(-1, -1);

			vertebraTemplate[z].sprocessLeft[i] = Vec2(-1, -1);
			vertebraTemplate[z].sprocessRight[i] = Vec2(-1, -1);

			vertebraTemplate[z].sprocessLeftWidth[i] = -1;
			vertebraTemplate[z].sprocessRightWidth[i] = -1;
		}
	}
	
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{

		// search for the tip for spinal process
		Vec2 tip, tipc;
		int count;

		k = z*sizexy;
		tip.x = tip.y = -1;

		for(y=segInfo.cordCenter[z].y+segInfo.cordRadius[z].y*3; y<sizey && y<segInfo.cordCenter[z].y+segInfo.cordRadius[z].y*10; y++)
		{
			k2 = y*sizey+segInfo.cordCenter[z].x-segInfo.cordRadius[z].x;
			for(x=segInfo.cordCenter[z].x-segInfo.cordRadius[z].x; x<=segInfo.cordCenter[z].x+segInfo.cordRadius[z].x; x++, k2++)
			{
				if(maskA[k2+k]==maskStruct.mask_corticalBone)
				{
					tip.x = x;
					tip.y = y;
				}
			}
		}

		if(tip.x==-1) continue;

		// search a 2*7 neighborhood to locate the centroid
		count = 0;
		tipc = Vec2(0,0);
		for(y=tip.y-1; y<=tip.y; y++)
		{
			k2 = y*sizey+tip.x-3;
			for(x=tip.x-3; x<=tip.x+3; x++, k2++)
			{
				if(maskA[k2+k]==maskStruct.mask_corticalBone)
				{
					count++;
					tipc += Vec2(x,y);
				}
			}
		}	// for y
		
		tipc /= count;
		vertebraTemplate[z].sprocessEnd = tipc;
		
		tip = segInfo.cordCenter[z];

		vertebraTemplate[z].sprocessLength = (tipc-tip).len();

		tipc = tipc-tip;
		vertebraTemplate[z].sprocessOrient = atan2(tipc.y, tipc.x);

		// visualize
		if(debugMode && false)
		{ 
			maskA[z*sizexy+(int)vertebraTemplate[z].sprocessEnd.y*sizex+(int)vertebraTemplate[z].sprocessEnd.x] = 
				maskStruct.mask_boneMetasis;
		}
	}	// for z

	// interpolate the process tip (length and orient)
	for(z=0; z<sizez; z++)
	{
		xt[z] = z;
		yt[z] = -1;
	}
	count = 0;
	minInterpolate = -1;
	maxInterpolate = -1;
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(vertebraTemplate[z].sprocessLength!=-1)
		{
			count++;
			yt[z] = vertebraTemplate[z].sprocessLength;
			if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
			if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
		}
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2_x, bernstein_power, piece_size);

	minInterpolate2 = -1;
	maxInterpolate2 = -1;
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(vertebraTemplate[z].sprocessOrient!=-1)
		{
			yt[z] = vertebraTemplate[z].sprocessOrient;
			if(yt[z]<minInterpolate2 || minInterpolate2==-1) minInterpolate2=yt[z];
			if(yt[z]>maxInterpolate2 || maxInterpolate2==-1) maxInterpolate2=yt[z];
		}
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2_y, bernstein_power, piece_size);

	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		Vec2 dir;
		if(vertebraTemplate[z].sprocessLength==-1 || vertebraTemplate[z].sprocessLength<yt2_x[z])
		{
			if(yt2_x[z]==-1) continue;
			if(yt2_x[z]<minInterpolate) yt2_x[z]=minInterpolate;
			if(yt2_x[z]>maxInterpolate) yt2_x[z]=maxInterpolate;

			if(yt2_y[z]==-1) continue;
			if(yt2_y[z]<minInterpolate2) yt2_y[z]=minInterpolate2;
			if(yt2_y[z]>maxInterpolate2) yt2_y[z]=maxInterpolate2;

			vertebraTemplate[z].sprocessLength = yt2_x[z];
			vertebraTemplate[z].sprocessOrient = yt2_y[z];

			dir.x = cos(yt2_y[z]);
			dir.y = sin(yt2_y[z]);

			vertebraTemplate[z].sprocessEnd = segInfo.cordCenter[z]+dir*yt2_x[z];

			// visualize
			if(debugMode && false)
				maskA[z*sizexy+(int)vertebraTemplate[z].sprocessEnd.y*sizex+(int)vertebraTemplate[z].sprocessEnd.x] = 
					maskStruct.mask_falseDetection;
		}
	}


	// search for left and right border of spinal process
	//
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{

		if(vertebraTemplate[z].sprocessEnd.x==-1) continue;

		k = z*sizexy;
		Vec2 tip1, tip2, medial[sprocessCount], dir, dirx, curPos, nextPos, maxPos;
		double steps, max_r;
		int r, grad, maxGrad, k3;

		tip1 = vertebraTemplate[z].sprocessEnd;
		tip2.x = segInfo.cordCenter[z].x;
		tip2.y = segInfo.cordCenter[z].y+segInfo.cordRadius[z].y;

		dir = (tip2-tip1).normalize();
		dirx = dir.perp();
		if(dirx.x<0)
		{
			dirx.x = -dirx.x;
			dirx.y = -dirx.y;
		}

		vertebraTemplate[z].sprocessLength = (tip2-tip1).len();
		steps = vertebraTemplate[z].sprocessLength/(sprocessCount-1);

		medial[0] = tip1;
		medial[sprocessCount-1] = tip2;

		for(i=1; i<sprocessCount-1; i++) medial[i] = medial[i-1]+dir*steps;

		for(i=0; i<sprocessCount; i++) vertebraTemplate[z].sprocessMedial[i] = medial[i];

		// search left and right from medial to locate spinal process edge
		//
		for(i=0; i<sprocessCount; i++)
		{
			// left
			curPos = medial[i]-dirx;
			maxGrad = 0;
			for(r=1; r<segInfo.cordRadius[z].x; r++)
			{
				nextPos = curPos-dirx;
				k2 = k + (int)(curPos.y+0.5)*sizex + (int)(curPos.x+0.5);
				k3 = k + (int)(nextPos.y+0.5)*sizex + (int)(nextPos.x+0.5);

				if(maskA[k2]==maskStruct.mask_corticalBone && maskA[k3]==maskStruct.mask_body)
				{
					grad = imgA[k2]-imgA[k3];
					if(grad>maxGrad && grad>segPara.borderMinGrad)
					{
						maxGrad=grad;
						maxPos = curPos;
						max_r = r;
					}
				}
				curPos = nextPos;
			}	// for r

			if(maxGrad>0)
			{
				vertebraTemplate[z].sprocessLeft[i] = maxPos;
				vertebraTemplate[z].sprocessLeftWidth[i] = max_r;
			}
			

			// right
			curPos = medial[i]+dirx;
			maxGrad = 0;
			for(r=1; r<segInfo.cordRadius[z].x; r++)
			{
				nextPos = curPos+dirx;
				k2 = k + (int)(curPos.y+0.5)*sizex + (int)(curPos.x+0.5);
				k3 = k + (int)(nextPos.y+0.5)*sizex + (int)(nextPos.x+0.5);

				if(maskA[k2]==maskStruct.mask_corticalBone && maskA[k3]==maskStruct.mask_body)
				{
					grad = imgA[k2]-imgA[k3];
					if(grad>maxGrad && grad>segPara.borderMinGrad)
					{
						maxGrad=grad;
						maxPos = curPos;
						max_r = r;
					}
				}
				curPos = nextPos;
			}	// for r

			if(maxGrad>0)
			{
				vertebraTemplate[z].sprocessRight[i] = maxPos;
				vertebraTemplate[z].sprocessRightWidth[i] = max_r;
			}			
		}	// for i
	}	// for z

	// interpolate left and right border of spinal process	(by left and right width)
	//
	// left side of spinal process
	for(i=0; i<sprocessCount; i++)
	{
		Vec2 dir, dirx;

		for(z=0; z<sizez; z++)
		{
			yt[z] = yt2_x[z] = -1;
		}
	
		count=0;
		minInterpolate = -1;
		maxInterpolate = -1;
		for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
		{
			if(vertebraTemplate[z].sprocessLeftWidth[i]!=-1)
			{
				count++;
				yt[z] = vertebraTemplate[z].sprocessLeftWidth[i];
				if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
				if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
			}
		}

		if(count<cordLength/2) continue;

		PiecewiseBernsteinSmoothing(xt, yt, yt2_x, bernstein_power, piece_size);

		for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
		{
			if(vertebraTemplate[z].sprocessMedial[i].x>0 && vertebraTemplate[z].sprocessLeftWidth[i]==-1 && yt2_x[z]>1)
			{
				if(yt2_x[z]<minInterpolate) yt2_x[z]=minInterpolate;
				if(yt2_x[z]>maxInterpolate) yt2_x[z]=maxInterpolate;
		
				dir.x = cos(vertebraTemplate[z].sprocessOrient);
				dir.y = sin(vertebraTemplate[z].sprocessOrient);
				dirx = dir.perp();

				if(dirx.x<0)
				{
					dirx.x = -dirx.x;
					dirx.y = -dirx.y;
				}
				
				vertebraTemplate[z].sprocessLeftWidth[i] = yt2_x[z];

				vertebraTemplate[z].sprocessLeft[i] = vertebraTemplate[z].sprocessMedial[i]-dirx*yt2_x[z];

				// visualize
				if(debugMode && false)
					maskA[z*sizexy+(int)vertebraTemplate[z].sprocessLeft[i].y*sizex+(int)vertebraTemplate[z].sprocessLeft[i].x] = 
						maskStruct.mask_falseDetection;
			}
		}
	}	// for i
	

	// right side of spinal process
	for(i=0; i<sprocessCount; i++)
	{
		Vec2 dir, dirx;

		for(z=0; z<sizez; z++)
		{
			yt[z] = yt2_x[z] = -1;
		}
	
		count=0;
		minInterpolate = -1;
		maxInterpolate = -1;
		for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
		{
			if(vertebraTemplate[z].sprocessRightWidth[i]!=-1)
			{
				count++;
				yt[z] = vertebraTemplate[z].sprocessRightWidth[i];
				if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
				if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
			}
		}

		if(count<cordLength/2) continue;

		PiecewiseBernsteinSmoothing(xt, yt, yt2_x, bernstein_power, piece_size);

		for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
		{
			if(vertebraTemplate[z].sprocessMedial[i].x>0 && vertebraTemplate[z].sprocessRightWidth[i]==-1 && yt2_x[z]>1)
			{
				if(yt2_x[z]<minInterpolate) yt2_x[z]=minInterpolate;
				if(yt2_x[z]>maxInterpolate) yt2_x[z]=maxInterpolate;
				dir.x = cos(vertebraTemplate[z].sprocessOrient);
				dir.y = sin(vertebraTemplate[z].sprocessOrient);
				dirx = dir.perp();

				if(dirx.x<0)
				{
					dirx.x = -dirx.x;
					dirx.y = -dirx.y;
				}
				
				vertebraTemplate[z].sprocessRightWidth[i] = yt2_x[z];

				vertebraTemplate[z].sprocessRight[i] = vertebraTemplate[z].sprocessMedial[i]+dirx*yt2_x[z];

				// visualize
				if(debugMode && false)
					maskA[z*sizexy+(int)vertebraTemplate[z].sprocessRight[i].y*sizex+(int)vertebraTemplate[z].sprocessRight[i].x] = 
						maskStruct.mask_falseDetection;
			}
		}
	}	// for i
	
	
	
	// convert the process ROI
	//
	vec2DynArray processROI;

	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		k = z*sizexy;
		processROI.SetSize(0);

		for(i=sprocessCount-1; i>=0; i--)
		{
			if(vertebraTemplate[z].sprocessLeft[i].x>0 && vertebraTemplate[z].sprocessLeft[i].x<sizex 
				&& vertebraTemplate[z].sprocessLeft[i].y>0 && vertebraTemplate[z].sprocessLeft[i].y<sizey)
			{
				processROI.Add(vertebraTemplate[z].sprocessLeft[i]);
			}
		}
		if(vertebraTemplate[z].sprocessEnd.x>0 && vertebraTemplate[z].sprocessEnd.x<sizex
			&& vertebraTemplate[z].sprocessEnd.y>0 && vertebraTemplate[z].sprocessEnd.y<sizey)
		{
			processROI.Add(vertebraTemplate[z].sprocessEnd);
		}
		for(i=0; i<sprocessCount; i++)
		{
			if(vertebraTemplate[z].sprocessRight[i].x>0 && vertebraTemplate[z].sprocessRight[i].x<sizex
				&& vertebraTemplate[z].sprocessRight[i].y>0 && vertebraTemplate[z].sprocessRight[i].y<sizey)
			{
				processROI.Add(vertebraTemplate[z].sprocessRight[i]);
			}
		}

		CIS_Algo_Contour_GetScanWhole(processROI, scan);

		for(i=0; i<scan.GetSize(); i++)
		{
			x = scan[i].x;
			y = scan[i].y;
			if(x<=0 || x>=sizex || y<=0 || y>=sizey) continue;

			k2 = k+y*sizex+x;

			if(maskA[k2]==maskStruct.mask_spinalCord) continue;

			//  base 200 for spinal process
			if(imgA[k2]>segPara.boneThresh) maskA[k2]=200+maskStruct.mask_corticalBone;
			else maskA[k2]=200+maskStruct.mask_spongyBone;
		}
	}


	// locate left and right pedicle template
	//

	// left pedicle
	for(z=0; z<sizez; z++)
	{
		vertebraTemplate[z].leftPedicleEnd = Vec2(-1, -1);
		vertebraTemplate[z].leftPedicleLength = -1;
		
		for(i=0; i<sprocessCount; i++) 
		{
			vertebraTemplate[z].leftPedicleMedial[i] = Vec2(-1, -1);
			vertebraTemplate[z].leftPedicleUp[i] = Vec2(-1, -1);
			vertebraTemplate[z].leftPedicleDown[i] = Vec2(-1, -1);
		}
	}
	
	// dynamic searching for the medial model in pedicle
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(segInfo.diskRadius[z].y<0) continue; 
		if(vertebraTemplate[z].sprocessEnd.x==-1) continue;

		// search for the tip of left pedicle
		Vec2 tip, tipc, dir, dirx, curPos, nextPos, maxPos;
		double search_r0, search_r1, search_a0, search_a1;
		double a, r, step_a;
		double  max_r, curScore, maxScore;
		short grad, maxGrad;
		int k3, k1, k4;

		k = z*sizexy;
		tip.x = tip.y = -1;

		// define the search radius range
		search_r0 = segInfo.diskRadius[z].y*0.5;
		search_r1 = segInfo.diskRadius[z].y*2.5;

		// define the search angle range
		dir = (vertebraTemplate[z].sprocessEnd-segInfo.cordCenter[z]).normalize();
		dirx = dir.perp();
		if(dirx.x>0)
		{
			dirx.x = -dirx.x; dirx.y=-dirx.y;
		}
		angle =atan2(dirx.y, dirx.x)*180/3.14;
		search_a0 = angle-45;
		search_a1 = angle+5;
		step_a = 50/search_r1;
		
		tip.x = segInfo.cordCenter[z].x;
		tip.y = segInfo.cordCenter[z].y;

		maxGrad = 0;
		max_r = 0;
		maxScore = 0;
		for(a=search_a0; a<=search_a1; a+=step_a)
		{
			dir.x = cos(a*3.14/180);
			dir.y = sin(a*3.14/180);

			curPos.x = tip.x+search_r0*dir.x;
			curPos.y = tip.y+search_r0*dir.y;
			for(r=search_r0; r<=search_r1; r++)
			{
				nextPos = curPos+dir;
				k2 = k + (int)(curPos.y+0.5)*sizex + (int)(curPos.x+0.5);
				k3 = k + (int)(nextPos.y+0.5)*sizex + (int)(nextPos.x+0.5);

				if(imgA[k2]>segPara.boneThresh && imgA[k3]<segPara.boneThresh && imgA[k3]>segPara.softThresh && maskA[k2]<100 
					&& imgA[k3-segInfo.cordRadius[z].x*2]>700 && imgA[k3-segInfo.cordRadius[z].x*2*sizex]>700)
				{
					k1 = k + (int)(curPos.y-dir.y+0.5)*sizex + (int)(curPos.x-dir.x+0.5);
					k4 = k + (int)(nextPos.y+dir.y+0.5)*sizex + (int)(nextPos.x+dir.x+0.5);
					grad = (imgA[k1]+imgA[k2]-imgA[k3]-imgA[k4])/2;
	
					if(grad>segPara.borderMinGrad)
					{
						curScore = PedicleTemplateMatching(img3D, maskImg3D, z, curPos, segPara, segInfo, maskStruct, vertebraTemplate, -1, false);

						if(curScore>maxScore)
						{
							maxScore = curScore;
							maxGrad=grad;
							max_r = r;
							maxPos = curPos;
						}
					}
				}

				curPos += dir;
			}	// for r
		}	// for a


		if(maxScore>0)
		{
			vertebraTemplate[z].leftPedicleEnd = maxPos;
			vertebraTemplate[z].leftPedicleLength = max_r;
		}

		if(vertebraTemplate[z].leftPedicleEnd.y>0)
		{
			if(debugMode && false)
			{
				// visualize it 
				k2 = k + (int)(vertebraTemplate[z].leftPedicleEnd.y+0.5)*sizex + (int)(vertebraTemplate[z].leftPedicleEnd.x+0.5);
				maskA[k2] = maskStruct.mask_spinalCord;
			}

			// fill the border using the tip
			PedicleTemplateMatching(img3D, maskImg3D, z, vertebraTemplate[z].leftPedicleEnd, segPara, segInfo, maskStruct, vertebraTemplate, -1, true);
		}
	}	// for z


	// right pedicle
	for(z=0; z<sizez; z++)
	{
		vertebraTemplate[z].rightPedicleEnd = Vec2(-1, -1);
		vertebraTemplate[z].rightPedicleLength = -1;
		
		for(i=0; i<sprocessCount; i++) 
		{
			vertebraTemplate[z].rightPedicleMedial[i] = Vec2(-1, -1);
			vertebraTemplate[z].rightPedicleUp[i] = Vec2(-1, -1);
			vertebraTemplate[z].rightPedicleDown[i] = Vec2(-1, -1);
		}
	}
	
	// dynamic searching for the medial model in pedicle
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{

		if(segInfo.diskRadius[z].y<0) continue; 
		if(vertebraTemplate[z].sprocessEnd.x==-1) continue;

		// search for the tip of right pedicle
		Vec2 tip, tipc, dir, dirx, curPos, nextPos, maxPos;
		double search_r0, search_r1, search_a0, search_a1;
		double a, r, step_a;
		double  max_r, curScore, maxScore;
		short grad, maxGrad;
		int k1, k4, k3;

		k = z*sizexy;
		tip.x = tip.y = -1;

		// define the search radius range
		search_r0 = segInfo.diskRadius[z].y*0.5;
		search_r1 = segInfo.diskRadius[z].y*2.5;

		// define the search angle range
		dir = (vertebraTemplate[z].sprocessEnd-segInfo.cordCenter[z]).normalize();
		dirx = dir.perp();
		if(dirx.x<0)
		{
			dirx.x = -dirx.x; dirx.y=-dirx.y;
		}
		angle =atan2(dirx.y, dirx.x)*180/3.14;
		search_a0 = angle-5;
		search_a1 = angle+45;
		step_a = 50/search_r1;
		
		tip.x = segInfo.cordCenter[z].x;
		tip.y = segInfo.cordCenter[z].y;

		maxGrad = 0;
		max_r = 0;
		maxScore = 0;
		for(a=search_a0; a<=search_a1; a+=step_a)
		{
			dir.x = cos(a*3.14/180);
			dir.y = sin(a*3.14/180);

			curPos.x = tip.x+search_r0*dir.x;
			curPos.y = tip.y+search_r0*dir.y;
			for(r=search_r0; r<=search_r1; r++)
			{
				nextPos = curPos+dir;
				k2 = k + (int)(curPos.y+0.5)*sizex + (int)(curPos.x+0.5);
				k3 = k + (int)(nextPos.y+0.5)*sizex + (int)(nextPos.x+0.5);

				if(imgA[k2]>segPara.boneThresh && imgA[k3]<segPara.boneThresh && imgA[k3]>segPara.softThresh && maskA[k2]<100 
					&& imgA[k3+segInfo.cordRadius[z].x*2]>700 && imgA[k3-segInfo.cordRadius[z].x*2*sizex]>700)
				{
					k1 = k + (int)(curPos.y-dir.y+0.5)*sizex + (int)(curPos.x-dir.x+0.5);
					k4 = k + (int)(nextPos.y+dir.y+0.5)*sizex + (int)(nextPos.x+dir.x+0.5);
					grad = (imgA[k1]+imgA[k2]-imgA[k3]-imgA[k4])/2;
	
					if(grad>segPara.borderMinGrad)
					{
						curScore = PedicleTemplateMatching(img3D, maskImg3D, z, curPos, segPara, segInfo, maskStruct, vertebraTemplate, 1, false);

						if(curScore>maxScore)
						{
							maxScore = curScore;
							maxGrad=grad;
							max_r = r;
							maxPos = curPos;
						}
					}
				}

				curPos += dir;
			}	// for r
		}	// for a


		if(maxScore>0)
		{
			vertebraTemplate[z].rightPedicleEnd = maxPos;
			vertebraTemplate[z].rightPedicleLength = max_r;
		}

		if(vertebraTemplate[z].rightPedicleEnd.y>0)
		{
			// visualize it 
			if(debugMode && false)
			{
				k2 = k + (int)(vertebraTemplate[z].rightPedicleEnd.y+0.5)*sizex + (int)(vertebraTemplate[z].rightPedicleEnd.x+0.5);
				maskA[k2] = maskStruct.mask_spinalCord;
			}

			// fill the border using the tip
			PedicleTemplateMatching(img3D, maskImg3D, z, vertebraTemplate[z].rightPedicleEnd, segPara, segInfo, maskStruct, vertebraTemplate, 1, true);
		}
	}	// for z


	// convert the pedicle ROI
	//
	vec2DynArray pedicleROI;

	// left pedicle
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		k = z*sizexy;
		pedicleROI.SetSize(0);

		for(i=pedicleCount-1; i>=0; i--)
		{
			if(vertebraTemplate[z].leftPedicleUp[i].x>0 && vertebraTemplate[z].leftPedicleUp[i].x<sizex
				&& vertebraTemplate[z].leftPedicleUp[i].y>0 && vertebraTemplate[z].leftPedicleUp[i].y<sizey)
			{
				pedicleROI.Add(vertebraTemplate[z].leftPedicleUp[i]);
			}
		}
		if(vertebraTemplate[z].leftPedicleEnd.x>0 && vertebraTemplate[z].leftPedicleEnd.x<sizex
			&& vertebraTemplate[z].leftPedicleEnd.y>0 && vertebraTemplate[z].leftPedicleEnd.y<sizey)
		{
			pedicleROI.Add(vertebraTemplate[z].leftPedicleEnd);
		}
		for(i=0; i<pedicleCount; i++)
		{
			if(vertebraTemplate[z].leftPedicleDown[i].x>0 && vertebraTemplate[z].leftPedicleDown[i].x<sizex
				&& vertebraTemplate[z].leftPedicleDown[i].y>0 && vertebraTemplate[z].leftPedicleDown[i].y<sizey)
			{
				pedicleROI.Add(vertebraTemplate[z].leftPedicleDown[i]);
			}
		}


		CIS_Algo_Contour_GetScanWhole(pedicleROI, scan);

		for(i=0; i<scan.GetSize(); i++)
		{
			x = scan[i].x;
			y = scan[i].y;
			if(x<=0 || x>=sizex || y<=0 || y>=sizey) continue;
			k2 = k+y*sizex+x;

			if(maskA[k2]==maskStruct.mask_spinalCord) continue;

			// base 300 for left pedicle
			if(imgA[k2]>segPara.boneThresh) maskA[k2]=300+maskStruct.mask_corticalBone;
			else maskA[k2]=300+maskStruct.mask_spongyBone;
		}
	}

	// right pedicle
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		k = z*sizexy;
		pedicleROI.SetSize(0);

		for(i=pedicleCount-1; i>=0; i--)
		{
			if(vertebraTemplate[z].rightPedicleUp[i].x>0 && vertebraTemplate[z].rightPedicleUp[i].x<sizex
				&& vertebraTemplate[z].rightPedicleUp[i].y>0 && vertebraTemplate[z].rightPedicleUp[i].y<sizey)
			{
				pedicleROI.Add(vertebraTemplate[z].rightPedicleUp[i]);
			}
		}
		if(vertebraTemplate[z].rightPedicleEnd.x>0 && vertebraTemplate[z].rightPedicleEnd.x<sizex
			&& vertebraTemplate[z].rightPedicleEnd.y>0 && vertebraTemplate[z].rightPedicleEnd.y<sizey)
		{
			pedicleROI.Add(vertebraTemplate[z].rightPedicleEnd);
		}
		for(i=0; i<pedicleCount; i++)
		{
			if(vertebraTemplate[z].rightPedicleDown[i].x>0 && vertebraTemplate[z].rightPedicleDown[i].x<sizex
				&& vertebraTemplate[z].rightPedicleDown[i].y>0 && vertebraTemplate[z].rightPedicleDown[i].y<sizey)
			{
				pedicleROI.Add(vertebraTemplate[z].rightPedicleDown[i]);
			}
		}


		CIS_Algo_Contour_GetScanWhole(pedicleROI, scan);

		for(i=0; i<scan.GetSize(); i++)
		{
			x = scan[i].x;
			y = scan[i].y;
			if(x<=0 || x>=sizex || y<=0 || y>=sizey) continue;
			k2 = k+y*sizex+x;

			if(maskA[k2]==maskStruct.mask_spinalCord) continue;

			// base 400 for right pedicle
			if(imgA[k2]>segPara.boneThresh) maskA[k2]=400+maskStruct.mask_corticalBone;
			else maskA[k2]=400+maskStruct.mask_spongyBone;
		}
	}


	// extract the last portion of the template, the region around cannel
	vec2DynArray cannelROI;
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		k = z*sizexy;
		cannelROI.SetSize(0);

		bool hasUL, hasUR, hasLL, hasLR;
		hasUL=hasUR=hasLL=hasLR=false;

		// for disk
		// get the first available disk point
		for(i=0; i<diskAngleCount; i++)
		{
			if(vertebraTemplate[z].diskContour[i].x>0 && vertebraTemplate[z].diskContour[i].x<sizex
				&& vertebraTemplate[z].diskContour[i].y>0 && vertebraTemplate[z].diskContour[i].y<sizey)
			{
				cannelROI.Add(vertebraTemplate[z].diskContour[i]);
				hasUR = true;
				break;
			}
		}
		// get the last available disk point
		for(i=diskAngleCount-1; i>=0; i--)
		{
			if(vertebraTemplate[z].diskContour[i].x>0 && vertebraTemplate[z].diskContour[i].x<sizex
				&& vertebraTemplate[z].diskContour[i].y>0 && vertebraTemplate[z].diskContour[i].y<sizey)
			{
				cannelROI.Add(vertebraTemplate[z].diskContour[i]);
				hasUL = true;
				break;
			}
		}

		// left pedicle
		// last point of upper border of left pedicle
		for(i=pedicleCount-1; i>=0; i--)
		{
			if(vertebraTemplate[z].leftPedicleUp[i].x>0 && vertebraTemplate[z].leftPedicleUp[i].x<sizex
				&& vertebraTemplate[z].leftPedicleUp[i].y>0 && vertebraTemplate[z].leftPedicleUp[i].y<sizey)
			{
				cannelROI.Add(vertebraTemplate[z].leftPedicleUp[i]);
				hasUL = true;
				break;
			}
		}
		// last point of lower border of left pedicle
		for(i=pedicleCount-1; i>=0; i--)
		{
			if(vertebraTemplate[z].leftPedicleDown[i].x>0 && vertebraTemplate[z].leftPedicleDown[i].x<sizex
				&& vertebraTemplate[z].leftPedicleDown[i].y>0 && vertebraTemplate[z].leftPedicleDown[i].y<sizey)
			{
				cannelROI.Add(vertebraTemplate[z].leftPedicleDown[i]);
				hasLL = true;
				break;
			}
		}

		// spinal process
		// last point of left border of spinal process
		for(i=sprocessCount-1; i>=0; i--)
		{
			if(vertebraTemplate[z].sprocessLeft[i].x>0 && vertebraTemplate[z].sprocessLeft[i].x<sizex
				&& vertebraTemplate[z].sprocessLeft[i].y>0 && vertebraTemplate[z].sprocessLeft[i].y<sizey)
			{
				cannelROI.Add(vertebraTemplate[z].sprocessLeft[i]);
				hasLL = true;
				break;
			}
		}

		// last point of right border of spinal process
		for(i=sprocessCount-1; i>=0; i--)
		{
			if(vertebraTemplate[z].sprocessRight[i].x>0 && vertebraTemplate[z].sprocessRight[i].x<sizex
				&& vertebraTemplate[z].sprocessRight[i].y>0 && vertebraTemplate[z].sprocessRight[i].y<sizey)
			{
				cannelROI.Add(vertebraTemplate[z].sprocessRight[i]);
				hasLR = true;
				break;
			}
		}

		// right pedicle
		// last point of lower border of right pedicle
		for(i=pedicleCount-1; i>=0; i--)
		{
			if(vertebraTemplate[z].rightPedicleDown[i].x>0 && vertebraTemplate[z].rightPedicleDown[i].x<sizex
				&& vertebraTemplate[z].rightPedicleDown[i].y>0 && vertebraTemplate[z].rightPedicleDown[i].y<sizey)
			{
				cannelROI.Add(vertebraTemplate[z].rightPedicleDown[i]);
				hasLR = true;
				break;
			}
		}
		// last point of upper border of right pedicle
		for(i=pedicleCount-1; i>=0; i--)
		{
			if(vertebraTemplate[z].rightPedicleUp[i].x>0 && vertebraTemplate[z].rightPedicleUp[i].x<sizex
				&& vertebraTemplate[z].rightPedicleUp[i].y>0 && vertebraTemplate[z].rightPedicleUp[i].y<sizey)
			{
				cannelROI.Add(vertebraTemplate[z].rightPedicleUp[i]);
				hasUR = true;
				break;
			}
		}

		if(hasUR && hasUL && hasLL && hasLR)
		{
			CIS_Algo_Contour_GetScanWhole(cannelROI, scan);

			for(i=0; i<scan.GetSize(); i++)
			{
				x = scan[i].x;
				y = scan[i].y;
				if(x<=0 || x>=sizex || y<=0 || y>=sizey) continue;
				k2 = k+y*sizex+x;

				// base 500 for spinal cord region
				if(maskA[k2]==maskStruct.mask_spinalCord || maskA[k2]==maskStruct.mask_spongyBone || maskA[k2]==maskStruct.mask_corticalBone) 
					maskA[k2]+=500;
			}
		}
		else
		{
			// keep everything within the bounding box of spinal cord region
			IntVec2 bb0, bb1;
			// first define a bounding box
			k = z*sizexy;
			int count =0;
			bb0.x = sizex-3;
			bb0.y = sizey-3;
			bb1.x = 3;
			bb1.y = 3;
			for(k2=0, y=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k++, k2++)
				{
					if(maskA[k]==maskStruct.mask_spinalCord)
					{
						if(x<bb0.x) bb0.x=x;
						if(x>bb1.x) bb1.x=x;
						if(y<bb0.y) bb0.y=y;
						if(y>bb1.y) bb1.y=y;
						count++;
					}
				}
			}
			if(count==0) continue;

			bb0.x-=6; bb0.y-=6; bb1.x+=6; bb1.y+=6;
			k = z*sizexy;
			count=0;
			for(y=bb0.y; y<=bb1.y; y++)
			{
				k2 = k+y*sizex+bb0.x;
				for(x=bb0.x; x<=bb1.x; x++, k2++)
				{
					// base 500 for spinal cord region
					if(maskA[k2]==maskStruct.mask_spinalCord || maskA[k2]==maskStruct.mask_spongyBone || maskA[k2]==maskStruct.mask_corticalBone) 
						maskA[k2]+=500;
				}
			}	// for y
		}	// else hasUR
	}


	// get overall template information
	for(z=0; z<sizez; z++)
	{
		vertebraTemplate[z].diskCenter.x = segInfo.diskCenter[z].x;
		vertebraTemplate[z].diskCenter.y = segInfo.diskCenter[z].y;
		vertebraTemplate[z].diskRadius = segInfo.diskRadius[z].y;

		vertebraTemplate[z].cordCenter.x = segInfo.cordCenter[z].x;
		vertebraTemplate[z].cordCenter.y = segInfo.cordCenter[z].y;
		vertebraTemplate[z].cordRadius.x = segInfo.cordRadius[z].x;
		vertebraTemplate[z].cordRadius.y = segInfo.cordRadius[z].y;
	}

	// compute the bounding box of each part of vertebra
	for(z=0; z<sizez; z++)
	{
		if(segInfo.cordCenter[z].x==-1) continue;

		// initialize
		segInfo.spineBound1[z] = IntVec2(sizex,sizey);
		segInfo.spineBound2[z] = IntVec2(-1,-1);
		segInfo.diskBound1[z] = IntVec2(sizex,sizey);
		segInfo.diskBound2[z] = IntVec2(-1,-1);
		segInfo.cordBound1[z] = IntVec2(sizex,sizey);
		segInfo.cordBound2[z] = IntVec2(-1,-1);
		segInfo.sprocessBound1[z] = IntVec2(sizex,sizey);
		segInfo.sprocessBound2[z] = IntVec2(-1,-1);

		k2 = z*sizexy;
		for(y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k2++)
			{
				if(maskA[k2]>100)
				{
					if(x<segInfo.spineBound1[z].x) segInfo.spineBound1[z].x=x;
					if(y<segInfo.spineBound1[z].y) segInfo.spineBound1[z].y=y;
					if(x>segInfo.spineBound2[z].x) segInfo.spineBound2[z].x=x;
					if(y>segInfo.spineBound2[z].y) segInfo.spineBound2[z].y=y;

					if(maskA[k2]<200)	// disk region
					{
						if(x<segInfo.diskBound1[z].x) segInfo.diskBound1[z].x=x;
						if(y<segInfo.diskBound1[z].y) segInfo.diskBound1[z].y=y;
						if(x>segInfo.diskBound2[z].x) segInfo.diskBound2[z].x=x;
						if(y>segInfo.diskBound2[z].y) segInfo.diskBound2[z].y=y;
					}
					else if(maskA[k2]<300)	// process
					{
						if(x<segInfo.sprocessBound1[z].x) segInfo.sprocessBound1[z].x=x;
						if(y<segInfo.sprocessBound1[z].y) segInfo.sprocessBound1[z].y=y;
						if(x>segInfo.sprocessBound2[z].x) segInfo.sprocessBound2[z].x=x;
						if(y>segInfo.sprocessBound2[z].y) segInfo.sprocessBound2[z].y=y;
					}
					else
					{
						if(maskA[k2]==500+maskStruct.mask_spinalCord)
						{
							if(x<segInfo.cordBound1[z].x) segInfo.cordBound1[z].x=x;
							if(y<segInfo.cordBound1[z].y) segInfo.cordBound1[z].y=y;
							if(x>segInfo.cordBound2[z].x) segInfo.cordBound2[z].x=x;
							if(y>segInfo.cordBound2[z].y) segInfo.cordBound2[z].y=y;
						}
					}
				}
			}
		}	// for y

		// put some margin on the bound
		segInfo.spineBound1[z].x -=8;
		segInfo.spineBound1[z].y -=8;
		segInfo.spineBound2[z].x +=8;
		segInfo.spineBound2[z].y +=8;

		segInfo.diskBound1[z].x -=4;
		segInfo.diskBound1[z].y -=4;
		segInfo.diskBound2[z].x +=4;
		segInfo.diskBound2[z].y +=4;

		// leave margin to image border
		if(segInfo.spineBound1[z].x<8) segInfo.spineBound1[z].x=8;
		if(segInfo.spineBound1[z].y<8) segInfo.spineBound1[z].y=8;
		if(segInfo.spineBound2[z].x>sizex-8) segInfo.spineBound2[z].x=sizex-8;
		if(segInfo.spineBound2[z].y>sizey-8) segInfo.spineBound2[z].y=sizey-8;

		if(segInfo.diskBound1[z].x<8) segInfo.diskBound1[z].x=8;
		if(segInfo.diskBound1[z].y<8) segInfo.diskBound1[z].y=8;
		if(segInfo.diskBound2[z].x>sizex-8) segInfo.diskBound2[z].x=sizex-8;
		if(segInfo.diskBound2[z].y>sizey-8) segInfo.diskBound2[z].y=sizey-8;
		if(z==73 || z==74)
		{
			std::cout<<"\n DiskBound1["<<z<<"]= ("<<segInfo.diskBound1[z].x<<","<<segInfo.diskBound1[z].y<<")\n";
			std::cout<<"DiskBound2["<<z<<"]= ("<<segInfo.diskBound2[z].x<<","<<segInfo.diskBound2[z].y<<")";
		}

	}	// fro z

	// convert the mask back
	for(k=0; k<sizexyz; k++)
	{
		// convert those outside of template to other bone
		if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_spongyBone || maskA[k]==maskStruct.mask_spinalCord) 
			maskA[k]=maskStruct.mask_otherBone;

		if(maskA[k]==maskStruct.mask_falseDetection) maskA[k] = maskStruct.mask_spongyBone;
		if(maskA[k]==maskStruct.mask_boneMetasis) maskA[k] = maskStruct.mask_spongyBone;

		if(maskA[k]>500) maskA[k]-=500;	// spinal cord
		else if(maskA[k]>400)			// right pedicle
		{
			if(imgA[k]>1100) maskA[k]-=400;
		}
		else if(maskA[k]>300)			// left pedicle
		{
			if(imgA[k]>1100) maskA[k]-=300;
		}	
		else if(maskA[k]>200)			// spinal process
		{
			if(imgA[k]>1100) maskA[k]-=200;
		}
		else if(maskA[k]>100) maskA[k]-=100;	// disk
	}

	// refinement on each 2D slice
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		
		int count;
		k = z*sizexy;
		// clean the buffer
		for(k2=0; k2<sizexy; k2++) binArray2D[k2]=0;

		// step 1. add back corticle bone
		// for every picel > boneThresh and within diskBound, add it back
		count = 0;
		for(y=segInfo.diskBound1[z].y-4; y<=segInfo.diskBound2[z].y; y++)
		{
			for(x=segInfo.diskBound1[z].x-4,k2=y*sizex+x; x<=segInfo.diskBound2[z].x+4; x++, k2++)
			{
				if(maskA[k+k2]==maskStruct.mask_otherBone && imgA[k+k2]>segPara.boneThresh) maskA[k+k2]=maskStruct.mask_corticalBone;
			}
		}	// for y

		// step 2: do a connected component analysis
		// every spine pixel inside spineBound should be connected expect for those inside sProcessBound
		for(y=segInfo.spineBound1[z].y-4; y<=segInfo.spineBound2[z].y; y++)
		{
			for(x=segInfo.spineBound1[z].x-4,k2=y*sizex+x; x<=segInfo.spineBound2[z].x+4; x++, k2++)
			{
				if(maskA[k+k2]==maskStruct.mask_spongyBone || maskA[k+k2]==maskStruct.mask_corticalBone || maskA[k+k2]==maskStruct.mask_spinalCord) 
				{
					binArray2D[k2]=1;
					count++;
				
					if(imgA[k+k2]>segPara.maxThresh) //to eliminate metal objects
					{
						binArray2D[k2]=0;
						maskA[k+k2]=maskStruct.mask_otherBone;
						count--;
					}
				}
			}
		}	// for y
		if(count==0) continue;
		NIH_Algo_Blob_Labelling(binImg2D, blobImg2D, blobRank2D, blobCentroid2D, false);

		int largest_b, largest_rank;
		int b;
		largest_rank = 0;
		for(b=0; b<blobRank2D.GetSize(); b++)
		{
			if(blobRank2D[b]>largest_rank)
			{
				largest_rank = blobRank2D[b];
				largest_b = b+1;
			}
		}

		for(y=segInfo.spineBound1[z].y; y<=segInfo.spineBound2[z].y; y++)
		{
			for(x=segInfo.spineBound1[z].x, k2=y*sizex+x; x<=segInfo.spineBound2[z].x; x++, k2++)
			{
				if(blobArray2D[k2]>0 && blobArray2D[k2]!=largest_b) 
				{
					if(x<segInfo.sprocessBound1[z].x || x>segInfo.sprocessBound2[z].x ||
						y<segInfo.sprocessBound1[z].y || y>segInfo.sprocessBound2[z].y)
						if(maskA[k+k2]!=maskStruct.mask_spinalCord) maskA[k+k2]=maskStruct.mask_otherBone;
				}
			}
		}	// for y

		// step 3: do a close operation to close hole and small gaps
		for(y=segInfo.spineBound1[z].y; y<=segInfo.spineBound2[z].y; y++)
		{
			for(x=segInfo.spineBound1[z].x, k2=y*sizex+x; x<=segInfo.spineBound2[z].x; x++, k2++)
			{
				if(maskA[k+k2]==maskStruct.mask_otherBone) binArray2D[k2]=0;
			}
		}	// for y
		CIS_IPA_Dilate(binImg2D, segPara.closeCorticalHoleIteration, false);
		CIS_IPA_Erode(binImg2D, segPara.closeCorticalHoleIteration, false);
		for(y=segInfo.spineBound1[z].y; y<=segInfo.spineBound2[z].y; y++)
		{
			k2 = y*sizex+segInfo.spineBound1[z].x;
			for(x=segInfo.spineBound1[z].x; x<=segInfo.spineBound2[z].x; x++, k2++)
			{
				if(binArray2D[k2]==1 && 
					maskA[k+k2]!=maskStruct.mask_corticalBone
					&& maskA[k+k2]!=maskStruct.mask_spinalCord) maskA[k+k2]=maskStruct.mask_spongyBone;
			}
		}	// for y
	}	// for z

	// recompute the spineBound1 and spineBound2 after exclude otherBone
	for(z=0; z<sizez; z++)
	{
		if(segInfo.cordCenter[z].x==-1) continue;

		// initialize
		segInfo.spineBound1[z] = IntVec2(sizex,sizey);
		segInfo.spineBound2[z] = IntVec2(-1,-1);

		k2 = z*sizexy;
		for(y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k2++)
			{
				if(maskA[k2]==maskStruct.mask_spongyBone || maskA[k2]==maskStruct.mask_corticalBone)
				{
					if(x<segInfo.spineBound1[z].x) segInfo.spineBound1[z].x=x;
					if(y<segInfo.spineBound1[z].y) segInfo.spineBound1[z].y=y;
					if(x>segInfo.spineBound2[z].x) segInfo.spineBound2[z].x=x;
					if(y>segInfo.spineBound2[z].y) segInfo.spineBound2[z].y=y;
				}
			}
		}	// for y

		// put some margin on the bound
		segInfo.spineBound1[z].x -=8;
		segInfo.spineBound1[z].y -=8;
		segInfo.spineBound2[z].x +=8;
		segInfo.spineBound2[z].y +=8;

		// leave margin to image border
		if(segInfo.spineBound1[z].x<8) segInfo.spineBound1[z].x=8;
		if(segInfo.spineBound1[z].y<8) segInfo.spineBound1[z].y=8;
		if(segInfo.spineBound2[z].x>sizex-8) segInfo.spineBound2[z].x=sizex-8;
		if(segInfo.spineBound2[z].y>sizey-8) segInfo.spineBound2[z].y=sizey-8;

	}	// for z

	delete binImg2D;
	delete blobImg2D;
	return CIS_OK;
}

// detect spinal cord and set up bounding box for each part of the vertebra (vetrabra body and spinal process)
// old method for low resolution data?? Graph based method
//
int NIH_SpineSegmentation_DetectSpinalCord(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   VertebraStruct2D *vertebraTemplate, 
									   bool debugMode)
{
	int x, y, z, k;
	int sizex, sizey, sizez, sizexy, sizexyz;
	short *imgA, *maskA;

	sizex = img3D->Num_Cols();
	sizey = img3D->Num_Rows();
	sizez = img3D->Num_Levels();

	imgA = img3D->GetArray();
	maskA = maskImg3D->GetArray();

	sizexyz = sizex*sizey*sizez;
	sizexy = sizex*sizey;

	CIS_Array_Image2D_short *binImg2D, *blobImg2D;
	short *binArray2D, *blobArray2D;
	intDynArray blobRank2D;
	intVec2DynArray blobCentroid2D;
	int largestBlob, largestSize;
	intDynArray blobStatus;

	int count, k2;
	int stx, sty, edx, edy;

	binImg2D = new CIS_Array_Image2D_short(sizex, sizey);
	blobImg2D = new CIS_Array_Image2D_short(sizex, sizey);
	binArray2D = binImg2D->GetArray();
	blobArray2D = blobImg2D->GetArray();

	typedef struct 
	{
		short index;
		short type;
		short slice;
		short centerx, centery;
		short bx0, bx1, by0, by1;	// bounding box
		short size;
		short link_up, link_down;
		float overlap_up, overlap_down;
		float score;
	} CORD_NODE;


	CORD_NODE *cord;
	int max_cord, num_cord;
	short mask_tmpbase=1000;
	intDynArray blob_by0, blob_by1;

	max_cord = sizez*10;
	num_cord = 0;
	cord = (CORD_NODE *)malloc(max_cord*sizeof(CORD_NODE));

	segInfo.bound1.z = -1;
	// first pass to get the candidate cord locations
	//
	for(z=0; z<sizez; z++)
	{
		// first define a bounding box
		k = z*sizexy;
		count =0;
		stx = sizex-3;
		sty = sizey-3;
		edx = 3;
		edy = 3;
		for(k2=0, y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++, k2++)
			{
				if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_spongyBone)
				{
					if(x<stx) stx=x;
					if(x>edx) edx=x;
					if(y<sty) sty=y;
					if(y>edy) edy=y;
					count++;
				}
			}
		}

		if(count<segPara.cordSizeThresh) continue;
		if(segInfo.bound1.z==-1) segInfo.bound1.z=z;
		segInfo.bound2.z = z;

		// first pass elinimate external space
		stx-=6; sty-=6; edx+=6; edy+=6;
		if(stx<0) stx=0;
		if(sty<0) sty=0;
		if(edx>=sizex-1) edx=sizex-2;
		if(edy>=sizey-1) edy=sizey-2;

		int sub_sizex, sub_sizey;
		int sk, sx, sy, nx, ny;
		CIS_Array_Image2D_short *sub_binImg2D, *sub_blobImg2D;
		short *subBinA, *subBlobA;
		sub_sizex = edx-stx+1;
		sub_sizey = edy-sty+1;

		sub_binImg2D = new CIS_Array_Image2D_short(sub_sizex, sub_sizey);
		sub_blobImg2D = new CIS_Array_Image2D_short(sub_sizex, sub_sizey);
		subBinA = sub_binImg2D->GetArray();
		subBlobA = sub_blobImg2D->GetArray();

		k = z*sizexy;
		count=0;
		for(y=sty, sk=0; y<=edy; y++)
		{
			k2 = y*sizex+stx;
			for(x=stx; x<=edx; x++, k2++, sk++)
			{
				if(maskA[k+k2]==maskStruct.mask_body && imgA[k+k2]<segPara.cordIntensityThresh)
				{
					subBinA[sk]=1;
					count ++;
				}
				else subBinA[sk]=0;
			}
		}	// for y

		if(count<segPara.cordSizeThresh) continue;

		// separate cord and external space
		CIS_IPA_Erode(sub_binImg2D, segPara.closeCordIteration, true);
		CIS_IPA_Dilate(sub_binImg2D, segPara.closeCordIteration, true);

		// closeness test to eliminate some outer space
/*		count=0;
		int hitC, s, searchR, nk, d, delta;
		int dx[8]={0,1,1,1,0,-1,-1,-1};
		int dy[8]={-1,-1,0,1,1,1,0,-1};
		// add margin
		searchR=50;
		for(sy=0, sk=0; sy<sub_sizey; sy++)
		{
			for(sx=0; sx<sub_sizex; sx++, sk++)
			{
				// if at margin, then eliminate
				if(sx<6 || sx>sub_sizex-6 || sy<6 || sy>sub_sizey-6) 
				{
					subBinA[sk]=2;
				} else if(subBinA[sk]==0)
				{
					k2 = k+(sy+sty)*sizex+sx+stx;
					if(maskA[k2]==maskStruct.mask_air) subBinA[sk]=2;
				}
			}
		}
		for(sy=0, sk=0; sy<sub_sizey; sy++)
		{
			for(sx=0; sx<sub_sizex; sx++, sk++)
			{
				if(subBinA[sk]==1)
				{
					hitC = 0;
					// test 0 oclock direction
					ny = sy-1;
					nk = sk-sub_sizex;
					for(s=1; s<searchR; s++)
					{
						if(subBinA[nk]==2)	// hit border pixels
						{
							subBinA[sk]=2;
							break;
						}
						else if(subBinA[nk]==0)	// hit bone
						{
							hitC++;
							break;
						}
						ny--;
						nk-=sub_sizex;
					}
					if(s==searchR) subBinA[sk]=2;
					if(subBinA[sk]==2) continue;

					// search the other seven directions
					for(d=1; d<8; d++)
					{
						nx = sx+dx[d];
						ny = sy+dy[d];
						delta = dy[d]*sub_sizex+dx[d];
						nk = sk+delta;
						for(s=1; s<searchR; s++)
						{
							if(subBinA[nk]==2)	// hit border
							{
								break;
							}
							else if(subBinA[nk]==0)	// hit bone
							{
								hitC++;
								break;
							}
							nx += dx[d];
							ny += dy[d];
							nk += delta;
						}
					}
					
					if(hitC>6) count ++;
					else subBinA[sk]=2;
				}
			}
		}	// for y

		if(count<segPara.cordSizeThresh) continue;

		for(sk=0; sk<sub_sizex*sub_sizey; sk++)
			if(subBinA[sk]==2) subBinA[sk]=0;

		// separate cord and external space
		CIS_IPA_Erode(sub_binImg2D, segPara.closeCordIteration, true);
		CIS_IPA_Dilate(sub_binImg2D, segPara.closeCordIteration, true);
*/
		NIH_Algo_Blob_Labelling(sub_binImg2D, sub_blobImg2D, blobRank2D, blobCentroid2D, false);

		blobStatus.SetSize(blobRank2D.GetSize());
		for(k=0; k<blobStatus.GetSize(); k++) blobStatus[k]=-1;

		blob_by0.SetSize(blobStatus.GetSize());
		blob_by1.SetSize(blobStatus.GetSize());
		for(k=0; k<blobStatus.GetSize(); k++) 
		{
			blob_by0[k]=edy;
			blob_by1[k]=sty;
		}

		// exclude the outer space, any blob next to border is considered in outer space
		int cb;	// current blob
		k = z*sizexy;
		count=0;
		for(y=sty, sk=0; y<=edy; y++)
		{
			k2 = y*sizex+stx;
			for(x=stx; x<=edx; x++, k2++, sk++)
			{
				if(subBlobA[sk]>=1)
				{
					cb = subBlobA[sk]-1;
					if(blobStatus[cb]==-2) continue;
					if(y<sty+segPara.closeCordIteration || y>=edy-segPara.closeCordIteration
						|| x<stx+segPara.closeCordIteration || x>=edx-segPara.closeCordIteration)
					{
						blobStatus[cb] = -2;	// outer border blob
					}
					// compute the bounding box of each blob, just y for now
					if(y<blob_by0[cb]) blob_by0[cb]=y;
					if(y>blob_by1[cb]) blob_by1[cb]=y;
				}
			}
		}	// for y

		for(k2=0; k2<blobStatus.GetSize(); k2++) 
		{
			if(blobStatus[k2]==-2) continue;
			
			if(blobRank2D[k2]<segPara.cordSizeThresh || blobRank2D[k2]>segPara.cordSizeThresh*50
				|| blobCentroid2D[k2].x+stx<sizex/2-sizex/10 
				|| blobCentroid2D[k2].x+stx>sizex/2+sizex/10) blobStatus[k2]=-1;		// too small or not in the center
			else
			{
				// the bounding box of cord should be in the center of the bounding box of vertebra
				// it should not near the top (the first quarter of the bounding box
				if(blob_by0[k2]<(sty*3+edy)/4) continue;

				if(num_cord>=max_cord)
				{
					max_cord = max_cord*2;
					cord = (CORD_NODE *)realloc(cord, max_cord*sizeof(CORD_NODE));
				}

				cord[num_cord].index = mask_tmpbase+num_cord;
				cord[num_cord].type = 0;
				cord[num_cord].slice = z;
				cord[num_cord].centerx = blobCentroid2D[k2].x+stx;
				cord[num_cord].centery = blobCentroid2D[k2].y+sty;
				cord[num_cord].bx0 = blobCentroid2D[k2].x+stx;
				cord[num_cord].bx1 = blobCentroid2D[k2].x+stx;
				cord[num_cord].by0 = blobCentroid2D[k2].y+sty;
				cord[num_cord].by1 = blobCentroid2D[k2].y+sty;
				cord[num_cord].size = blobRank2D[k2];
				cord[num_cord].link_up = -1;
				cord[num_cord].link_down = -1;
				cord[num_cord].overlap_up = 0;
				cord[num_cord].overlap_down = 0;
				cord[num_cord].score = 0;

				blobStatus[k2] = mask_tmpbase+num_cord;	// blobStatus now become cord index
				num_cord++;
			}
		}	// for k2
		
		// assign a temp mask value for cord candidate pixels (different candidate has different mask value)
		k = z*sizexy;

		for(y=sty, sk=0; y<=edy; y++)
		{
			k2 = y*sizex+stx;
			for(x=stx; x<=edx; x++, k2++, sk++)
			{
				if(subBlobA[sk]>=1)
				{
					cb = subBlobA[sk]-1;
					if(blobStatus[cb]>=0)
					{
						maskA[k+k2]=blobStatus[cb];
					}
				}
			}
		}	// for y

		delete sub_binImg2D;
		delete sub_blobImg2D;
	}	// for z


//	for(k=0; k<sizexyz; k++) if(maskA[k]>=mask_tmpbase) maskA[k]=maskStruct.mask_spinalCord;
//	return CIS_ERROR;

	int org_num_cord = num_cord;

	//
	// propagate the cord candidates to fill gaps between slices, interpolation
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		// initialize the 2D image
		for(k2=0; k2<sizexy; k2++) binArray2D[k2]=0;

		k = z*sizexy;
		count =0;
		if(z==0)
		{
			for(k2=0, y=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k++, k2++)
				{
					if(maskA[k]==maskStruct.mask_body && imgA[k]<segPara.cordIntensityThresh && 
						(maskA[k+sizexy]>=mask_tmpbase && maskA[k+sizexy]<mask_tmpbase+org_num_cord) )
					{
						binArray2D[k2] = 1;
						count++;
					}
				}
			}
		}
		else if(z==sizez-1)
		{
			for(k2=0, y=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k++, k2++)
				{
					if(maskA[k]==maskStruct.mask_body && imgA[k]<segPara.cordIntensityThresh && 
						(maskA[k-sizexy]>=mask_tmpbase && maskA[k-sizexy]<mask_tmpbase+org_num_cord) )
					{
						binArray2D[k2] = 1;
						count++;
					}
				}
			}
		}
		else
		{
			for(k2=0, y=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k++, k2++)
				{
					if(maskA[k]==maskStruct.mask_body && imgA[k]<segPara.cordIntensityThresh && 
						((maskA[k-sizexy]>=mask_tmpbase && maskA[k-sizexy]<mask_tmpbase+org_num_cord) 
						|| (maskA[k+sizexy]>=mask_tmpbase && maskA[k+sizexy]<mask_tmpbase+org_num_cord)))
					{
						binArray2D[k2] = 1;
						count++;
					}
				}
			}
		}

		if(count<segPara.cordSizeThresh) continue;
		NIH_Algo_Blob_Labelling(binImg2D, blobImg2D, blobRank2D, blobCentroid2D, false);

		blobStatus.SetSize(blobRank2D.GetSize());

		// add new candidate nodes
		for(k2=0; k2<blobStatus.GetSize(); k2++) 
		{
			if(blobRank2D[k2]<segPara.cordSizeThresh) blobStatus[k2]=-1;
			else
			{
				if(num_cord>=max_cord)
				{
					max_cord = max_cord*2;
					cord = (CORD_NODE *)realloc(cord, max_cord*sizeof(CORD_NODE));
				}

				cord[num_cord].index = mask_tmpbase+num_cord;
				cord[num_cord].type = 1;
				cord[num_cord].slice = z;
				cord[num_cord].centerx = blobCentroid2D[k2].x;
				cord[num_cord].centery = blobCentroid2D[k2].y;
				cord[num_cord].bx0 = blobCentroid2D[k2].x;
				cord[num_cord].bx1 = blobCentroid2D[k2].x;
				cord[num_cord].by0 = blobCentroid2D[k2].y;
				cord[num_cord].by1 = blobCentroid2D[k2].y;
				cord[num_cord].size = blobRank2D[k2];
				cord[num_cord].link_up = -1;
				cord[num_cord].link_down = -1;
				cord[num_cord].overlap_up = 0;
				cord[num_cord].overlap_down = 0;
				cord[num_cord].score = 0;

				blobStatus[k2] = mask_tmpbase+num_cord;
				num_cord++;
			}
		}	// for k2
		
		k = z*sizexy;
		for(k2=0, y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++, k2++)
			{
				if(binArray2D[k2]==1) 
				{
					if(blobStatus[blobArray2D[k2]-1]!=-1)
					{
						maskA[k]=blobStatus[blobArray2D[k2]-1];
					}
				}	// if binArray
			}
		}
	}	// for z

//	for(k=0; k<sizexyz; k++) if(maskA[k]>=mask_tmpbase) maskA[k]=maskStruct.mask_spinalCord;
//	return CIS_ERROR;

	// compute the overlap and score of each node candidate
	intDynArray list_up, list_down, over_up, over_down;
	int n, tn, i;
	short mask_tmplast;
	float largestOverlap, c_overlap;
	short overIndex;

	mask_tmplast = mask_tmpbase+num_cord;

	for(n=0; n<num_cord; n++)
	{
		z = cord[n].slice;
		k = z*sizexy;
		list_up.SetSize(0);
		list_down.SetSize(0);
		over_up.SetSize(0);
		over_down.SetSize(0);

		for(y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++)
			{
				if(maskA[k]==mask_tmpbase+n)
				{
					// check overlap with previous slice
					if(z>0)
					{
						tn = maskA[k-sizexy]-mask_tmpbase;
						if(tn>=0 && tn<num_cord)
						{
							for(i=0; i<list_up.GetSize(); i++)
							{
								if(list_up[i]==tn)
								{
									break;
								}
							}
							if(i>=list_up.GetSize())
							{
								list_up.Add(tn);
								over_up.Add(1);
							}
							else
							{
								over_up[i]++;
							}
						}	// if tn
					}	// if z

					// check overlap with next slice
					if(z<sizez-1)
					{
						tn = maskA[k+sizexy]-mask_tmpbase;
						if(tn>=0 && tn<num_cord)
						{
							for(i=0; i<list_down.GetSize(); i++)
							{
								if(list_down[i]==tn)
								{
									break;
								}
							}
							if(i>=list_down.GetSize())
							{
								list_down.Add(tn);
								over_down.Add(1);
							}
							else
							{
								over_down[i]++;
							}
						}	// if tn
					}	// if z
				}	// if maskA
			}	// for x
		}	// for y

		// found the largest overlap node in previous and next slice
		largestOverlap = 0;
		overIndex = -1;
		for(i=0; i<list_up.GetSize(); i++)
		{
			c_overlap = (float)over_up[i]/(float)(cord[n].size+cord[list_up[i]].size-over_up[i]);
			if(c_overlap>largestOverlap && c_overlap>0.2)
			{
				largestOverlap = c_overlap;
				overIndex = list_up[i];
			}
		}
		cord[n].link_up = overIndex;
		cord[n].overlap_up = largestOverlap;

		largestOverlap = 0;
		overIndex = -1;
		for(i=0; i<list_down.GetSize(); i++)
		{
			c_overlap = (float)over_down[i]/(float)(cord[n].size+cord[list_down[i]].size-over_down[i]);
			if(c_overlap>largestOverlap && c_overlap>0.2)
			{
				largestOverlap = c_overlap;
				overIndex = list_down[i];
			}
		}
		cord[n].link_down = overIndex;
		cord[n].overlap_down = largestOverlap;

	}	// for n

	// find the longest path in the graph
	// compute the score
	int cn, nn;
	float c_score, a_score;;
	for(n=0; n<num_cord; n++)
	{
		a_score = 1;
		// going up first
		cn = n;
		c_score = 1;
		do {
			nn = cord[cn].link_up;
			if(nn==-1) break;
//			c_score = c_score*cord[cn].overlap_up;
			c_score = cord[cn].overlap_up;
			a_score += c_score;
			cn = nn;
		} while(1);

		// going down then
		cn = n;
		c_score = 1;
		do {
			nn = cord[cn].link_down;
			if(nn==-1) break;
//			c_score = c_score*cord[cn].overlap_down;
			c_score = cord[cn].overlap_down;
			a_score += c_score;
			cn = nn;
		} while(1);

		cord[n].score = a_score;
	}

	// adjust the link based on scores
	for(n=0; n<num_cord; n++)
	{
		nn = cord[n].link_down;
		if(nn!=-1 && cord[nn].link_up!=n)
		{
			if(cord[nn].link_up==-1) cord[nn].link_up=n;
			else
			{
				if(cord[n].score>cord[cord[nn].link_up].score) 
				{
					cord[nn].link_up=n;	
					cord[nn].overlap_up=cord[n].overlap_down;	
				}
			}
		}

		nn = cord[n].link_up;
		if(nn!=-1 && cord[nn].link_down!=n)
		{
			if(cord[nn].link_down==-1) cord[nn].link_down=n;
			else
			{
				if(cord[n].score>cord[cord[nn].link_down].score) 
				{
					cord[nn].link_down=n;	
					cord[nn].overlap_down=cord[n].overlap_up;	
				}
			}
		}	
	}	// for n

	// recompute the score
	for(n=0; n<num_cord; n++)
	{
		a_score = 1;
		// going up first
		cn = n;
		c_score = 1;
		do {
			nn = cord[cn].link_up;
			if(nn==-1) break;
			c_score = cord[cn].overlap_up;
			a_score += c_score;
			cn = nn;
		} while(1);

		// going down then
		cn = n;
		c_score = 1;
		do {
			nn = cord[cn].link_down;
			if(nn==-1) break;
			c_score = cord[cn].overlap_down;
			a_score += c_score;
			cn = nn;
		} while(1);

		if(a_score>=cord[n].score) cord[n].score = a_score;
		else
		{
			cn=n;	// something is wrong
		}
	}

	// only keep the cord with largest score on each slice
	intDynArray sliceCord;
	doubleDynArray sliceCordScore;
	sliceCord.SetSize(sizez);
	sliceCordScore.SetSize(sizez);
	for(z=0; z<sizez; z++)
	{
		sliceCord[z] = -1;
		sliceCordScore[z] = 0;
	}

	blobStatus.SetSize(num_cord);
	for(n=0; n<num_cord; n++)
	{
		z = cord[n].slice;
		if(cord[n].score>sliceCordScore[z])
		{
			sliceCord[z] = n;
			sliceCordScore[z] = cord[n].score;
		}

		blobStatus[n] = 0;
	}

	// correct the spinal cord if it is not a continuous chain
	int chain_st, chain_ed, all_st, all_ed;
	intDynArray segment_st, segment_ed;
	int numSeg;

	z = segInfo.bound1.z;
	for(z=segInfo.bound1.z; z<segInfo.bound2.z; z++)
	{
		if(sliceCord[z]!=-1) break;
	}
	all_st = segInfo.bound2.z;
	all_ed = segInfo.bound1.z;
	numSeg = 0;
	segment_st.SetSize(0);
	segment_ed.SetSize(0);
	while(z<segInfo.bound2.z)
	{
		chain_st = z;
		chain_ed = z;
		segment_st.Add(z);
		segment_ed.Add(z);
		while(sliceCord[z]!=-1 && z<sizez-1 && sliceCord[z+1]!=-1 && cord[sliceCord[z]].link_down==sliceCord[z+1])
		{
			if(z<all_st) all_st=z;

			z++;
			chain_ed = z;
			segment_ed[numSeg] = z;

			if(z>all_ed) all_ed=z;
		}

		numSeg++;
		z++;
	}	// while z
	
	// merge multiple segments
	intDynArray segStatus;
	int numMergedSeg, curSeg, nextSeg;
	numMergedSeg = 0;
	segStatus.SetSize(numSeg);
	for(n=0; n<numSeg; n++) segStatus[n]=-1;
	for(n=0; n<numSeg; n++)
	{
		if(segStatus[n]!=-1) continue;
		curSeg = n;
		segStatus[curSeg] = numMergedSeg;

		int cn, nn, cz, nz, ck, nk;
		int overlap;
		for(nextSeg=curSeg+1; nextSeg<numSeg; nextSeg++)
		{
			if(segStatus[nextSeg]!=-1) continue;
			cz = segment_ed[curSeg];
			nz = segment_st[nextSeg];
			// compute the overlap of cn and nn
			cn = sliceCord[cz];
			nn = sliceCord[nz];
			ck = cz*sizexy;
			nk = nz*sizexy;
			overlap = 0;
			for(y=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, ck++, nk++)
				{
					if(maskA[ck]==mask_tmpbase+cn && maskA[nk]==mask_tmpbase+nn)
					{
						overlap++;
					}	// if maskA[ck];
				}	// for x
			}	// for y

			if(overlap>segPara.cordSizeThresh/2)
			{
				segStatus[nextSeg] = numMergedSeg;
				curSeg = nextSeg;
			}
		}	// for nextSeg
		numMergedSeg++;
	}	//	for n

	// pick the longest merged seg
	int longestMergedSeg, longestLength;
	intDynArray mergedSegLength;
	mergedSegLength.SetSize(numMergedSeg);
	for(n=0; n<numMergedSeg; n++) mergedSegLength[n]=0;
	for(n=0; n<numSeg; n++)
	{
		if(segStatus[n]!=-1) mergedSegLength[segStatus[n]]+= (segment_ed[n]-segment_st[n]+1);
	}
	longestMergedSeg = 0;
	longestLength=0;
	for(n=0; n<numMergedSeg; n++) 
	{
		if(mergedSegLength[n]>longestLength) 
		{
			longestMergedSeg=n;
			longestLength = mergedSegLength[n];
		}
	}

	int lastMergedSlice=0;
	for(n=0; n<numSeg; n++)
	{
		if(segStatus[n]==longestMergedSeg)
		{
			for(z=segment_st[n]; z<=segment_ed[n]; z++)
			{
				if(sliceCord[z]!=-1) blobStatus[sliceCord[z]]=1;
			}
			if(lastMergedSlice<segment_ed[n]) lastMergedSlice=segment_ed[n];
		}
		else
		{
			for(z=segment_st[n]; z<=segment_ed[n]; z++)
			{
				sliceCord[z]=-1;
			}
		}
	}	// for n


	// trace downwards
	int pz, cz, pn, ck, pk;
	for(z=1; z<sizez; z++)
	{
		if(sliceCord[z]!=-1) continue;
		cz = z;
		pz = z-1;
		if(sliceCord[pz]==-1) continue;
		pn = sliceCord[pz];

		pk=pz*sizexy;
		ck=cz*sizexy;

		over_down.SetSize(num_cord);
		for(n=0; n<num_cord; n++) over_down[n]=0;

		for(y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, ck++, pk++)
			{
				if(maskA[pk]==mask_tmpbase+pn)
				{
					tn = maskA[ck]-mask_tmpbase;
					if(tn>=0 && tn<num_cord)
					{
						over_down[tn]++;
					}
				}	// if maskA[pk];
			}	// for x
		}	// for y

		// get the largest overlap
		largestOverlap = 0;
		overIndex = -1;
		for(n=0; n<num_cord; n++)
		{
			if(over_down[n]>largestOverlap && over_down[n]>cord[pn].size/4)
			{
				largestOverlap = over_down[n];
				overIndex = n;
			}
		}

		if(overIndex>=0)
		{
			sliceCord[cz] = overIndex;
			cord[overIndex].link_up = pn;
			cord[overIndex].overlap_up = largestOverlap;
			cord[pn].link_down = overIndex;
			cord[pn].overlap_down = largestOverlap;

			pn = overIndex;
			pz = cz;
			cz = pz+1;

			while(cord[pn].link_down>=0)
			{
				cn = cord[pn].link_down;
				sliceCord[cz] = cn;
				pn = cn;
				pz = cz;
				cz = pz+1;
			}
			if(z>lastMergedSlice) lastMergedSlice=z;
		}
		else //
		{
			// add a new node using the overlap
			largestOverlap = 0;
			k = cz*sizexy;
			for(k2=0, y=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k++, k2++)
				{
					if(maskA[k-sizexy]== mask_tmpbase+pn // && maskA[k]==maskStruct.mask_body
						&& imgA[k]<segPara.cordIntensityThresh) 
					{
						maskA[k]=mask_tmpbase+num_cord;
						largestOverlap ++;
					}
				}	// if maskA
			}	// for k2

///			if(largestOverlap<segPara.cordSizeThresh) break;

			if(num_cord>=max_cord)
			{
				max_cord = max_cord*2;
				cord = (CORD_NODE *)realloc(cord, max_cord*sizeof(CORD_NODE));
			}

			cord[num_cord].index = mask_tmpbase+num_cord;
			cord[num_cord].type = 0;
			cord[num_cord].slice = cz;
			cord[num_cord].centerx = cord[pn].centerx;
			cord[num_cord].centery = cord[pn].centery;
			cord[num_cord].bx0 = cord[pn].bx0;
			cord[num_cord].bx1 = cord[pn].bx1;
			cord[num_cord].by0 = cord[pn].by0;
			cord[num_cord].by1 = cord[pn].by1;
			cord[num_cord].size = largestOverlap;
			cord[num_cord].link_up = pn;
			cord[num_cord].link_down = -1;
			cord[num_cord].overlap_up = largestOverlap;
			cord[num_cord].overlap_down = 0;
			cord[num_cord].score = cord[pn].score;

			cord[pn].link_down = num_cord;
			cord[pn].overlap_down = largestOverlap;

			sliceCord[cz] = num_cord;

			over_down.Add(0);

			pn = num_cord;
			pz = cz;
			cz = pz+1;

			num_cord++;
		}	// else
	} // while z

	// delete all slices after lastMergedSlice
	for(z=lastMergedSlice; z<sizez; z++) sliceCord[z]=-1;

	// assign the cord
	blobStatus.SetSize(num_cord);
	for(n=0; n<num_cord; n++) blobStatus[n]=0;
	for(z=0; z<sizez; z++)
	{
		if(sliceCord[z]>=0) blobStatus[sliceCord[z]]=1;
	}

	// assign the mask
	for(k=0; k<sizexyz; k++)
	{
		if(maskA[k]>=mask_tmpbase && maskA[k]<mask_tmpbase+num_cord) 
		{
			if(blobStatus[maskA[k]-mask_tmpbase]==1) maskA[k] = maskStruct.mask_spinalCord;
///			else maskA[k] = maskStruct.mask_spongyBone;	// convert the holes
			else maskA[k] = maskStruct.mask_body;	// convert the holes
		}
	}


	// dilate the cord to fill the gap
	IntVec2 seed;
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		k = z*sizexy;
		count = 0;
		for(k2=0; k2<sizexy; k2++, k++) 
		{
			if(maskA[k]==maskStruct.mask_spinalCord) 
			{
				count++;
				binArray2D[k2]=1;
			}
			else binArray2D[k2]=0;
		}

		if(count==0) continue;

		CIS_IPA_Dilate(binImg2D, segPara.closeCordIteration, true);

		k=z*sizexy;
		for(k2=0; k2<sizexy; k2++, k++) 
		{
			if(binArray2D[k2]==1 && (maskA[k]!=maskStruct.mask_corticalBone && maskA[k]!=maskStruct.mask_spongyBone)) 
				maskA[k]=maskStruct.mask_spinalCord;
		}

		// only keep the largest piece of cord
		k=z*sizexy;
		for(y=0, k2=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k2++, k++)
			{
				if(maskA[k]==maskStruct.mask_spinalCord)
				{
					binArray2D[k2]=1;
				}
				else binArray2D[k2]=0;
			}
		}
		NIH_Algo_Blob_Labelling(binImg2D, blobImg2D, blobRank2D, blobCentroid2D, false);
		if(blobRank2D.GetSize()>1)
		{
			int largestRank, largestRankI=1;
			largestRank=0;
			for(k=0; k<blobRank2D.GetSize(); k++)
			{
				if(blobRank2D[k]>largestRank)
				{
					largestRank = blobRank2D[k];
					largestRankI = k+1;
				}
			}
			k = z*sizexy;
			for(k2=0, y=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k++, k2++)
				{
					if(binArray2D[k2]==1) 
					{
						if(blobArray2D[k2]!=largestRankI)
						{
							maskA[k]=maskStruct.mask_body;
						}
					}	// if binArray
				}
			}
		}
		
		// find the bounding box of cord
		k=z*sizexy;
		for(y=0, k2=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k2++, k++)
			{
				if(maskA[k]==maskStruct.mask_spinalCord)
				{
					if(x<cord[sliceCord[z]].bx0) cord[sliceCord[z]].bx0 = x;
					if(x>cord[sliceCord[z]].bx1) cord[sliceCord[z]].bx1 = x;
					if(y<cord[sliceCord[z]].by0) cord[sliceCord[z]].by0 = y;
					if(y>cord[sliceCord[z]].by1) cord[sliceCord[z]].by1 = y;
				}
			}
		}
	}	// for z

	int cordStartZ, cordEndZ, cordLength;
	cordStartZ=0;
	for(z=0; z<sizez; z++)
	{
		if(sliceCord[z]!=-1) 
		{
			cordStartZ = z;
			break;
		}
	}
	cordEndZ=sizez-2;
	for(; z<sizez; z++)
	{
		if(sliceCord[z]==-1) 
		{
			cordEndZ = z-1;
			break;
		}
	}
	cordLength = cordEndZ-cordStartZ+1;

	segInfo.bound1.z = cordStartZ;
	segInfo.bound2.z = cordEndZ;

	// fill the interior holes with spongy bones
	for(z=0; z<sizez; z++)
	{
		// remove the slices with no spinal cannel
		if(sliceCord[z]==-1)
		{
			k = z*sizexy;
			for(k2=0, y=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k++, k2++)
				{
					if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_spongyBone)
					{
						maskA[k] = maskStruct.mask_otherBone;
					}
				}
			}

			continue;
		}

		// first define a bounding box
		k = z*sizexy;
		count =0;
		stx = sizex-3;
		sty = sizey-3;
		edx = 3;
		edy = 3;
		for(k2=0, y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++, k2++)
			{
				if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_spongyBone)
				{
					if(x<stx) stx=x;
					if(x>edx) edx=x;
					if(y<sty) sty=y;
					if(y>sty) edy=y;
					count++;
				}
			}
		}
		stx-=6; sty-=6; edx+=6; edy+=6;
		if(stx<2) stx=2;
		if(sty<2) sty=2;
		if(edx>sizex-3) edx=sizex-3;
		if(edy>sizey-3) edy=sizey-3;

		if(count==0) continue;

		// initialize the 2D image
		for(k2=0; k2<sizexy; k2++) binArray2D[k2]=0;

		k = z*sizexy;
		count=0;
		for(y=sty; y<=edy; y++)
		{
			k2 = y*sizex+stx;
			for(x=stx; x<=edx; x++, k2++)
			{
				if(maskA[k+k2]==maskStruct.mask_body)
				{
					binArray2D[k2]=1;
					count ++;
				}
			}
		}	// for y

		if(count==0) continue;

		// separate interior and external space
		CIS_IPA_Erode(binImg2D, segPara.closeCordIteration, true);
		CIS_IPA_Dilate(binImg2D, segPara.closeCordIteration, true);

		NIH_Algo_Blob_Labelling(binImg2D, blobImg2D, blobRank2D, blobCentroid2D, false);

		blobStatus.SetSize(blobRank2D.GetSize());
		for(k=0; k<blobStatus.GetSize(); k++) blobStatus[k]=1;

		// Get the outer space, any blob next to border are considered outer space
		k = z*sizexy;
		count=0;
		for(y=sty; y<=edy; y++)
		{
			k2 = y*sizex+stx;
			for(x=stx; x<=edx; x++, k2++)
			{
				if(binArray2D[k2]==1)
				{
					if(y<sty+segPara.closeCordIteration || y>=edy-segPara.closeCordIteration
						|| x<stx+segPara.closeCordIteration || x>=edx-segPara.closeCordIteration)
					{
						blobStatus[blobArray2D[k2]-1] = -2;
					}
				}
			}
		}

		// assign mask
		k = z*sizexy;
		for(k2=0, y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++, k2++)
			{
				if(binArray2D[k2]==1) 
				{
					if(blobStatus[blobArray2D[k2]-1]>=0)
					{
						maskA[k]=maskStruct.mask_spongyBone;
						// also fill its neighbors
						if(maskA[k+1]==maskStruct.mask_body) maskA[k+1]==maskStruct.mask_spongyBone;
						if(maskA[k-1]==maskStruct.mask_body) maskA[k-1]==maskStruct.mask_spongyBone;
						if(maskA[k+sizex]==maskStruct.mask_body) maskA[k+sizex]==maskStruct.mask_spongyBone;
						if(maskA[k-sizex]==maskStruct.mask_body) maskA[k-sizex]==maskStruct.mask_spongyBone;
					}
				}	// if binArray
			}
		}

	}	// for z


	// compute the center of the spinal cannel and radius of cannel
	// and find the bounding box for each slice
	IntVec2 c_center, c_bound1, c_bound2;
	int len_y1, len_y2, len_x1, len_x2, largest_len;
	
	segInfo.spineBound1.SetSize(sizez);
	segInfo.spineBound2.SetSize(sizez);
	segInfo.diskBound1.SetSize(sizez);
	segInfo.diskBound2.SetSize(sizez);
	segInfo.cordBound1.SetSize(sizez);
	segInfo.cordBound2.SetSize(sizez);
	segInfo.sprocessBound1.SetSize(sizez);
	segInfo.sprocessBound2.SetSize(sizez);

	segInfo.cordCenter.SetSize(sizez);
	segInfo.cordRadius.SetSize(sizez);

	segInfo.diskCenter.SetSize(sizez);
	segInfo.diskRadius.SetSize(sizez);

	for(z=0; z<sizez; z++)
	{
		segInfo.spineBound1[z] = IntVec2(-1,-1);
		segInfo.spineBound2[z] = IntVec2(-1,-1);
		segInfo.diskBound1[z] = IntVec2(-1,-1);
		segInfo.diskBound2[z] = IntVec2(-1,-1);
		segInfo.cordBound1[z] = IntVec2(-1,-1);
		segInfo.cordBound2[z] = IntVec2(-1,-1);
		segInfo.sprocessBound1[z] = IntVec2(-1,-1);
		segInfo.sprocessBound2[z] = IntVec2(-1,-1);

		segInfo.cordCenter[z] = IntVec2(-1,-1);
		segInfo.cordRadius[z] = IntVec2(-1,-1);
		segInfo.diskCenter[z] = IntVec2(-1,-1);
		segInfo.diskRadius[z] = IntVec2(-1,-1);
	}

	doubleDynArray xt, yt, yt2;
	double minInterpolate, maxInterpolate, minInterpolate2, maxInterpolate2;
	int bernstein_power = 5, piece_size=5;
	float piece_length=25;
	float pixelSizez=img3D->Get_Pixel_SizeZ();
	piece_size = (int)(piece_length/pixelSizez+0.5);

	xt.SetSize(sizez);
	yt.SetSize(sizez);
	yt2.SetSize(sizez);
	for(z=0; z<sizez; z++)
	{
		xt[z] = z;
		yt[z] = -1;
	}

	// interpolate x coordinate of cord
	minInterpolate = -1;
	maxInterpolate = -1;
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;

		yt[z] = cord[sliceCord[z]].centerx;
		if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
		if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;
		if(yt2[z]==-1) continue;
		if(yt2[z]<minInterpolate) yt2[z]=minInterpolate;
		if(yt2[z]>maxInterpolate) yt2[z]=maxInterpolate;
		segInfo.cordCenter[z].x = (int)yt2[z];
	}

	// interpolate y coordinate of cord
	minInterpolate = -1;
	maxInterpolate = -1;
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;

		yt[z] = cord[sliceCord[z]].centery;
		if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
		if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;
		if(yt2[z]==-1) continue;
		if(yt2[z]<minInterpolate) yt2[z]=minInterpolate;
		if(yt2[z]>maxInterpolate) yt2[z]=maxInterpolate;
		segInfo.cordCenter[z].y = (int)yt2[z];
	}


	// interpolate x radius of cord
	minInterpolate = -1;
	maxInterpolate = -1;
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;

		yt[z] = (cord[sliceCord[z]].bx1-cord[sliceCord[z]].bx0)/2;
		if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
		if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;
		if(yt2[z]==-1) continue;
		if(yt2[z]<minInterpolate) yt2[z]=minInterpolate;
		if(yt2[z]>maxInterpolate) yt2[z]=maxInterpolate;
		segInfo.cordRadius[z].x = (int)yt2[z];
	}

	// interpolate y radius of cord
	minInterpolate = -1;
	maxInterpolate = -1;
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;

		yt[z] = (cord[sliceCord[z]].by1-cord[sliceCord[z]].by0)/2;
		if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
		if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;
		if(yt2[z]==-1) continue;
		if(yt2[z]<minInterpolate) yt2[z]=minInterpolate;
		if(yt2[z]>maxInterpolate) yt2[z]=maxInterpolate;
		segInfo.cordRadius[z].y = (int)yt2[z];
	}

	// compute the center of the disk, and radius of disk
	// interpolate the radius of the disk

	for(z=0; z<sizez; z++)
	{
		xt[z] = z;
		yt[z] = -1;
	}
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;

		// initial value of diskRaius
		segInfo.diskRadius[z].y = 2*segInfo.cordRadius[z].y;

		c_center.x = cord[sliceCord[z]].centerx;
		c_center.y = cord[sliceCord[z]].centery;

		int y_top, y_count;
		y_top = y_count = 0;
		k = z*sizexy;
		for(x=c_center.x-1; x<=c_center.x+1; x++)
		{
			for(y=c_center.y-segInfo.cordRadius[z].y*8; y<c_center.y-segInfo.cordRadius[z].y*3; y++)
			{
				k2 = k+y*sizex+x;
				if(maskA[k2]==maskStruct.mask_spongyBone || maskA[k2]==maskStruct.mask_corticalBone)
				{
					y_count++;
					y_top += y;
					break;
				}
			}
		}

		if(y_count>0)
		{
			y_top/=y_count;
			if(y_top<c_center.y) 
			{
				yt[z] = (c_center.y-y_top)/2;
				segInfo.diskRadius[z].y = yt[z];
			}
		}
	}	// for z

	// refine the radius of the disk, 
	// computing the radius of top half of the circle
	// then average it with the current radius
	//
	for(z=cordStartZ; z<=cordEndZ; z++)
	{
		if(sliceCord[z]==-1) continue;
		Vec2 dir, curPos, nextPos, maxPos, maxRadius;
		int search_r1, search_r0, r, k3, grad, maxGrad, max_r;
		int angle, angle1;

		c_center.x = segInfo.cordCenter[z].x;
		c_center.y = segInfo.cordCenter[z].y-segInfo.diskRadius[z].y;

		int angleLeft, angleRight;
		angleLeft = 270;
		angleRight = 90;

		int count=0;
		int accu_radius=0;

		search_r0 = 1;
		search_r1 = segInfo.diskRadius[z].y+segInfo.cordRadius[z].y;

		for(angle=angleRight; angle<=angleLeft; angle+=diskAngleInterval)
		{
			// top 
			dir.x = sin(angle*3.14/180);
			dir.y = cos(angle*3.14/180);


			curPos.x = c_center.x;
			curPos.y = c_center.y;
			k = z*sizexy;
			maxGrad = 0;
			curPos = curPos+dir*search_r0;
			for(r=search_r0; r<=search_r1; r++)
			{
				nextPos = curPos+dir;
				k2 = k + (int)(curPos.y+0.5)*sizex + (int)(curPos.x+0.5);
				k3 = k + (int)(nextPos.y+0.5)*sizex + (int)(nextPos.x+0.5);

				if(imgA[k2]>segPara.boneThresh && imgA[k3]<segPara.boneThresh && imgA[k3]>1000)
				{
					grad = imgA[k2]-imgA[k3];
	
					if(grad>maxGrad && grad>50)
					{
						maxGrad=grad;
						max_r = r;
						maxPos = curPos;
					}
				}

				curPos += dir;
			}

			if(maxGrad>0)
			{
				accu_radius += max_r;
				count++;
			}
		}	// for angle

		if(count>0)
		{
			accu_radius /= count;
			segInfo.diskRadius[z].y = (segInfo.diskRadius[z].y+accu_radius)/2;
		}
		
	}	// for z
	
	// constrain the disk radius of the first and last few slices to avoid wild extrapolation
	// make sure it is greater than cord radius and smaller than 4 times of it
	minInterpolate = -1;
	maxInterpolate = -1;
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1)  continue;

		yt[z] = segInfo.diskRadius[z].y;
		if(yt[z]<segInfo.cordRadius[z].y)
		{
			yt[z] = segInfo.cordRadius[z].y;
		}
		else if(yt[z]>segInfo.cordRadius[z].y*4)
		{
			yt[z] = segInfo.cordRadius[z].y*4;
		}
		if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
		if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
	}

	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;

		if(yt2[z]==-1) continue;
		if(yt2[z]<minInterpolate) yt2[z]=minInterpolate;
		if(yt2[z]>maxInterpolate) yt2[z]=maxInterpolate;

		segInfo.diskCenter[z].x = segInfo.cordCenter[z].x;

		if(yt2[z]>0)
		{
			segInfo.diskRadius[z].y = (int)yt2[z];

			segInfo.diskCenter[z].y = segInfo.cordCenter[z].y-(int)yt2[z];
		}

		// visalize it 
		if(debugMode && false)
		{
			if(yt[z]>0) maskA[z*sizexy+(int)yt[z]*sizex+cord[sliceCord[z]].centerx] = maskStruct.mask_corticalBone;
			maskA[z*sizexy+(int)yt2[z]*sizex+segInfo.diskCenter[z].x] = maskStruct.mask_boneMetasis;
		}
	
	}

	// determine the start of sacrum
	// due to the complication, the sacrum won't be used in segmentation
	//
	// the method is test the left and right side of spinal process, 
	// and check if there is sundden increase of cortical bones
	// if there are, which usually means that pelvis is very close to the spine column and was segmented as part of the spine, 
	// which also indicates that sacrum may start
	//
	segInfo.sacrum_start=-1;
	int counti, counto, last_counto, pelvis_slices=0;
	bool sacrum_exist=false;
	last_counto = 0;
	for(z=cordEndZ; z>=sizez/2; z--)		// only search the last half of the image
	{
		if(sliceCord[z]==-1) continue;
		int count=0, counta=0;

		k = z*sizexy;

		// first get the largest connected component in 2D
		//
		for(k2=0; k2<sizexy; k2++)
		{
			if(maskA[k+k2]==maskStruct.mask_spongyBone || maskA[k+k2]==maskStruct.mask_corticalBone || maskA[k+k2]==maskStruct.mask_spinalCord) 
			{
				binArray2D[k2]=1;
				count++;
			}
			else binArray2D[k2]=0;
		}

		if(count==0) continue;
		NIH_Algo_Blob_Labelling(binImg2D, blobImg2D, blobRank2D, blobCentroid2D, false);

		int largest_b, largest_rank;
		int b;
		largest_rank = 0;
		for(b=0; b<blobRank2D.GetSize(); b++)
		{
			if(blobRank2D[b]>largest_rank)
			{
				largest_rank = blobRank2D[b];
				largest_b = b+1;
			}
		}


		counti = counto = counta= 0;

		for(y=segInfo.cordCenter[z].y+segInfo.cordRadius[z].y; y<segInfo.cordCenter[z].y+segInfo.cordRadius[z].y*3; y++)
		{
			if(y>=sizey) break;
			for(x=segInfo.cordCenter[z].x-segInfo.cordRadius[z].x*8; x<segInfo.cordCenter[z].x-segInfo.cordRadius[z].x*2; x++) 
			{
				if(x<0) continue;
				k2 = y*sizex+x;
				b = blobArray2D[k2];
				if(b==largest_b && (maskA[k2+k]==maskStruct.mask_corticalBone || maskA[k2+k]==maskStruct.mask_spongyBone)) 
				{
					counto++;
				}
				else if(maskA[k2+k]==maskStruct.mask_air && y<segInfo.cordCenter[z].y+segInfo.cordRadius[z].y+5
					&& x>segInfo.cordCenter[z].x-segInfo.cordRadius[z].x*3
					&& x<segInfo.cordCenter[z].x+segInfo.cordRadius[z].x*3) counta++;
			}
			for(x=segInfo.cordCenter[z].x+segInfo.cordRadius[z].x*2; x<segInfo.cordCenter[z].x+segInfo.cordRadius[z].x*8; x++) 
			{
				if(x>=sizex) break;
				k2 = y*sizex+x;
				b = blobArray2D[k2];
				if(b==largest_b && (maskA[k2+k]==maskStruct.mask_corticalBone || maskA[k2+k]==maskStruct.mask_spongyBone)) 
				{
					counto++;
				}
				else if(maskA[k2+k]==maskStruct.mask_air && y<segInfo.cordCenter[z].y+segInfo.cordRadius[z].y+5
					&& x>segInfo.cordCenter[z].x-segInfo.cordRadius[z].x*3
					&& x<segInfo.cordCenter[z].x+segInfo.cordRadius[z].x*3) counta++;
			}
		}	// for y

		if(counto<last_counto/2 && (pelvis_slices>(cordEndZ-z)/2 || pelvis_slices>2) && counta<10)
		{
			segInfo.sacrum_start = z+1;
			break;
		}
		
		if(counto>500) pelvis_slices++;
		else pelvis_slices=0;

		// if detect pelvis, sacrum must exist
		if(pelvis_slices>10) sacrum_exist=true;

		// if sacrum exist, then the first clean vertebra is marked as sacrum_start
		if(sacrum_exist && counto<10 && last_counto<10)
		{
			segInfo.sacrum_start = z+1;
			break;
		}

		last_counto = counto;
	}

	if(segInfo.sacrum_start>0)
	{
		for(z=segInfo.sacrum_start; z<=cordEndZ; z++)
		{
			k=z*sizexy;

			// convert outer bones to other bones
			for(k2=0; k2<sizexy; k2++)
			{
				if(maskA[k2+k]==maskStruct.mask_corticalBone || maskA[k2+k]==maskStruct.mask_spongyBone 
					|| maskA[k2+k]==maskStruct.mask_spinalCord) maskA[k2+k]=maskStruct.mask_otherBone;
			}
		}

		cordEndZ = segInfo.sacrum_start-1;
		segInfo.bound2.z = cordEndZ;
	}
	cordLength = cordEndZ-cordStartZ+1;

	// locate and refine the vertebra template
	// it contains three consequential steps, 
	// 1. the disk
	// 2. spinal process
	// 3. left and right pedicle
	//
	
	// interpolate the contour of disk

	// based on center of the disk, search for a few points along the boundary
	//

	for(z=0; z<sizez; z++)
	{
		for(int a=0; a<diskAngleCount; a++)
		{
			vertebraTemplate[z].diskContour[a] = Vec2(-1, -1);
			vertebraTemplate[z].diskContourRadius[a] = -1;
			vertebraTemplate[z].interpolatedRadius[a] = -1;
		}
		vertebraTemplate[z].diskNeckLeft = Vec2(-1, -1);
		vertebraTemplate[z].diskNeckRight = Vec2(-1, -1);
	}
	
	for(z=cordStartZ; z<=cordEndZ; z++)
	{
		if(sliceCord[z]==-1) continue;
		Vec2 dir, curPos, nextPos, maxPos, maxRadius;
		int search_r1, search_r0, r, k3, grad, maxGrad, max_r;
		int angle, angle1;

		// locate left and right of vertebra neck
		c_center.x = segInfo.cordCenter[z].x;
		c_center.y = segInfo.cordCenter[z].y;

		y = c_center.y;
		// left side
		k = z*sizexy;
		k2 = k+y*sizey+c_center.x-segInfo.cordRadius[z].x;
		maxGrad = 0;
		for(x=c_center.x-segInfo.cordRadius[z].x; x>=c_center.x-segInfo.cordRadius[z].x*3; x--, k2--)
		{
			if(maskA[k2]==maskStruct.mask_spinalCord)
			{
				maxGrad = 1;
				maxPos = Vec2(x,y);
			}
			else if(maskA[k2]==maskStruct.mask_corticalBone && maskA[k2-1]!=maskStruct.mask_corticalBone && imgA[k2-1]>1000)
			{
				grad = imgA[k2]-imgA[k2-1];

				if(grad>maxGrad && grad>50)
				{
					maxGrad = grad;
					maxPos = Vec2(x,y);
				}
			}
		}
		if(maxGrad>0)
		{
			vertebraTemplate[z].diskNeckLeft = maxPos;
		}
		else
		{
			vertebraTemplate[z].diskNeckLeft = Vec2(c_center.x-segInfo.cordRadius[z].x,y);
		}


		// right side
		k = z*sizexy;
		k2 = k+y*sizey+c_center.x+segInfo.cordRadius[z].x;
		maxGrad = 0;
		for(x=c_center.x+segInfo.cordRadius[z].x; x<=c_center.x+segInfo.cordRadius[z].x*3; x++, k2++)
		{
			if(maskA[k2]==maskStruct.mask_spinalCord)
			{
				maxGrad = 1;
				maxPos = Vec2(x,y);
			}
			else if(maskA[k2]==maskStruct.mask_corticalBone && maskA[k2+1]!=maskStruct.mask_corticalBone && imgA[k2+1]>1000)
			{
				grad = imgA[k2]-imgA[k2+1];

				if(grad>maxGrad && grad>50)
				{
					maxGrad = grad;
					maxPos = Vec2(x,y);
				}
			}
		}
		if(maxGrad>0)
		{
			vertebraTemplate[z].diskNeckRight = maxPos;
		}
		else
		{
			vertebraTemplate[z].diskNeckRight = Vec2(c_center.x+segInfo.cordRadius[z].x,y);
		}

		// raise the neck to the top of spinal cord
///		vertebraTemplate[z].diskNeckLeft.y -= segInfo.cordRadius[z].y;
///		vertebraTemplate[z].diskNeckRight.y -= segInfo.cordRadius[z].y;

		c_center.x = segInfo.diskCenter[z].x;
		c_center.y = segInfo.diskCenter[z].y;

		int angleLeft, angleRight;
		angleLeft = (int)(atan2(vertebraTemplate[z].diskNeckLeft.x-c_center.x, vertebraTemplate[z].diskNeckLeft.y-c_center.y)*180/3.14)+360;
		angleRight = (int)(atan2(vertebraTemplate[z].diskNeckRight.x-c_center.x, vertebraTemplate[z].diskNeckRight.y-c_center.y)*180/3.14);

		angleLeft = (angleLeft/diskAngleInterval+1)*diskAngleInterval;
		angleRight = (angleRight/diskAngleInterval-1)*diskAngleInterval;
		vertebraTemplate[z].angleLeft = angleLeft;
		vertebraTemplate[z].angleRight = angleRight;
		
		for(angle=angleRight; angle<=angleLeft; angle+=diskAngleInterval)
		{
			angle1 = angle/diskAngleInterval;

			// top 
			dir.x = sin(angle*3.14/180);
			dir.y = cos(angle*3.14/180);

			search_r0 = segInfo.diskRadius[z].y-segInfo.cordRadius[z].y;
			search_r1 = segInfo.diskRadius[z].y+segInfo.cordRadius[z].y;
			x = c_center.x;

			curPos.x = c_center.x;
			curPos.y = c_center.y;
			k = z*sizexy;
			maxGrad = 0;
			curPos = curPos+dir*search_r0;
			for(r=search_r0; r<=search_r1; r++)
			{
				if(curPos.y>=segInfo.cordCenter[z].y) break;

				nextPos = curPos+dir;
				k2 = k + (int)(curPos.y+0.5)*sizex + (int)(curPos.x+0.5);
				k3 = k + (int)(nextPos.y+0.5)*sizex + (int)(nextPos.x+0.5);
				if(maskA[k2]==maskStruct.mask_spinalCord)
				{
					if(maxGrad==0)
					{
						maxGrad==1;
						max_r = r;
						maxPos = curPos;
					}
					break;
				}

				if(imgA[k2]>segPara.boneThresh && imgA[k3]<segPara.boneThresh && imgA[k3]>1000)
				{
					grad = imgA[k2]-imgA[k3];
	
					if(grad>maxGrad && grad>50)
					{
						maxGrad=grad;
						max_r = r;
						maxPos = curPos;
					}
				}

				curPos += dir;
			}

			if(maxGrad>0)
			{
				vertebraTemplate[z].diskContour[angle1] = maxPos;
				vertebraTemplate[z].diskContourRadius[angle1] = max_r;
			}
		}	// for angle
		
	}	// for z

	// interpolate every point on the contour
	// interpolate the radius
	// do the interpolation twice to refine the result
	//
	doubleDynArray yt2_x, yt2_y;
	int angle;
	yt2_x.SetSize(sizez);
	yt2_y.SetSize(sizez);
	for(angle=0; angle<360; angle+=diskAngleInterval)
	{
		int angle1 = angle/diskAngleInterval;

		int count;
		for(z=0; z<sizez; z++)
		{
			yt[z] = -1;
		}

		// interpolate the radius at each angle on the disk
		count = 0;
		minInterpolate = -1;
		maxInterpolate = -1;
		for(z=cordStartZ; z<=cordEndZ; z++)
		{
			if(vertebraTemplate[z].diskContourRadius[angle1]>segInfo.cordRadius[z].y)
			{
				count++;
				yt[z] = vertebraTemplate[z].diskContourRadius[angle1];
				if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
				if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
			}
		}

		if(count<cordLength/2) continue;

		PiecewiseBernsteinSmoothing(xt, yt, yt2_x, bernstein_power, piece_size);

		Vec2 dir;
		dir.x = sin(angle*3.14/180);
		dir.y = cos(angle*3.14/180);

		for(z=cordStartZ; z<=cordEndZ; z++)
		{
			if(yt2_x[z]==-1) continue;
			if(yt2_x[z]<minInterpolate) yt2_x[z]=minInterpolate;
			if(yt2_x[z]>maxInterpolate) yt2_x[z]=maxInterpolate;
			vertebraTemplate[z].interpolatedRadius[angle1] = yt2_x[z];
		}

	}	// for angle

	// quality assurance
	//
	for(z=cordStartZ; z<=cordEndZ; z++)
	{
		double avg_r, dev_r, r1, r2,r3;
		double count_r;
		int angle, angle1;
		Vec2 dir, curPos, nextPos, maxPos, maxRadius;
		int search_r1, search_r0, r, k3, grad, maxGrad, max_r;

		avg_r = 0;
		dev_r = 0;
		count_r = 0;

		// compute the average radius
		for(angle=0; angle<360; angle+=diskAngleInterval)
		{
			angle1 = angle/diskAngleInterval;

			r1 = -1;
			if(vertebraTemplate[z].diskContourRadius[angle1]>segInfo.cordRadius[z].y)
			{
				r1 = vertebraTemplate[z].diskContourRadius[angle1];
			}
			else if(vertebraTemplate[z].interpolatedRadius[angle1]>segInfo.cordRadius[z].y)
			{
				r1 = vertebraTemplate[z].interpolatedRadius[angle1];
			}

			if(r1>segInfo.cordRadius[z].y && r1<segInfo.diskRadius[z].y*5 && r1>segInfo.diskRadius[z].y/5)
			{
				avg_r += r1;
				count_r ++;
			}
		}

		if(count_r!=0) avg_r /= count_r;

		// compute the deviation
		for(angle=0; angle<360; angle+=diskAngleInterval)
		{
			angle1 = angle/diskAngleInterval;

			r1 = -1;
			if(vertebraTemplate[z].diskContourRadius[angle1]>segInfo.cordRadius[z].y)
			{
				r1 = vertebraTemplate[z].diskContourRadius[angle1];
			}
			else if(vertebraTemplate[z].interpolatedRadius[angle1]>segInfo.cordRadius[z].y)
			{
				r1 = vertebraTemplate[z].interpolatedRadius[angle1];
			}

			if(r1>segInfo.cordRadius[z].y && r1<segInfo.diskRadius[z].y*5 && r1>segInfo.diskRadius[z].y/5)
			{
				dev_r += (r1-avg_r)*(r1-avg_r);
				count_r ++;
			}

		}

		if(count_r>1) dev_r = sqrt(dev_r/(count_r-1));

		// pick out the outliers
		for(angle=0; angle<360; angle+=diskAngleInterval)
		{
			angle1 = angle/diskAngleInterval;
			r1 = vertebraTemplate[z].diskContourRadius[angle1];
			r2 = vertebraTemplate[z].interpolatedRadius[angle1];

			if(r1>segInfo.cordRadius[z].y)
			{
				if(fabs(r1-avg_r)<dev_r*2 && r1<segInfo.diskRadius[z].y*5 && r1>segInfo.diskRadius[z].y/5) continue;		// good cases
				else if(r2>segInfo.cordRadius[z].y && fabs(r2-avg_r)<dev_r*2 && r2<segInfo.diskRadius[z].y*5 && r2>segInfo.diskRadius[z].y/5)
				{
					vertebraTemplate[z].diskContourRadius[angle1]=r2;
				}
				else
				{
					vertebraTemplate[z].diskContourRadius[angle1]=-1;
				}
			}
			else if(r2>segInfo.cordRadius[z].y && fabs(r2-avg_r)<dev_r*2 && r2<segInfo.diskRadius[z].y*5 && r2>segInfo.diskRadius[z].y/5)
			{
				vertebraTemplate[z].diskContourRadius[angle1]=r2;				
			}
		}

		// smoothing on each slice
		int angleNumber=360/diskAngleInterval;
		xt.SetSize(angleNumber);
		yt.SetSize(angleNumber);
		yt2.SetSize(angleNumber);
		minInterpolate = -1;
		maxInterpolate = -1;
		for(angle1=0; angle1<angleNumber; angle1++)
		{
			xt[(int)angle1] = (int)angle1;
			yt[(int)angle1] = vertebraTemplate[z].diskContourRadius[angle1];
			if(yt[(int)angle1]==-1) continue;
			if(yt[(int)angle1]<minInterpolate || minInterpolate==-1) minInterpolate=yt[(int)angle1];
			if(yt[(int)angle1]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[(int)angle1];
		}

		BernsteinSmoothing(xt, yt, yt2, bernstein_power);
		for(angle=0; angle<360; angle+=diskAngleInterval)
		{
			angle1 = angle/diskAngleInterval;

			if(angle<vertebraTemplate[z].angleRight || angle>vertebraTemplate[z].angleLeft) continue;
			if(yt2[(int)angle1]==-1) continue;
			if(yt2[(int)angle1]<minInterpolate) yt2[(int)angle1]=minInterpolate;
			if(yt2[(int)angle1]>maxInterpolate) yt2[(int)angle1]=maxInterpolate;
			vertebraTemplate[z].diskContourRadius[angle1] = yt2[angle1];
		}

		k = z*sizexy;
		// local refinement and set contour based on radius
		for(angle=0; angle<360; angle+=diskAngleInterval)
		{
			dir.x = sin(angle*3.14/180);
			dir.y = cos(angle*3.14/180);
		
			angle1 = angle/diskAngleInterval;
			r1 = vertebraTemplate[z].diskContourRadius[angle1];
			if(r1<0 && angle>vertebraTemplate[z].angleRight && angle<vertebraTemplate[z].angleLeft) r1=segInfo.diskRadius[z].y;
			if(r1>segInfo.cordRadius[z].y)
			{
				vertebraTemplate[z].diskContour[angle1].x = -1;
				vertebraTemplate[z].diskContour[angle1].y = -1;
				// do a local search to find the local maximum gradient
				//
				search_r0 = r1-2;
				search_r1 = r1+4;
				curPos = segInfo.diskCenter[z]+dir*search_r0;
				maxGrad = 0;
				for(r=search_r0; r<=search_r1; r++)
				{
					nextPos = curPos+dir;
					k2 = k + (int)(curPos.y+0.5)*sizex + (int)(curPos.x+0.5);
					k3 = k + (int)(nextPos.y+0.5)*sizex + (int)(nextPos.x+0.5);
					if(imgA[k3]>1000)
					{
						grad = imgA[k2]-imgA[k3];
						if(grad>maxGrad && grad>50)
						{
							maxGrad=grad;
							max_r = r;
							maxPos = curPos;
						}
					}

					curPos = nextPos;
				}	// for r
		
				if(maxGrad>0)
				{
					r1 = max_r;
				}

				if(r1<segInfo.diskRadius[z].y*5 && r1>segInfo.diskRadius[z].y/5)
				{
					vertebraTemplate[z].diskContourRadius[angle1] = r1;
				}
			}	// if r1
		}	// for angle

		// set contour
		for(angle=0; angle<360; angle+=diskAngleInterval)
		{
			dir.x = sin(angle*3.14/180);
			dir.y = cos(angle*3.14/180);
			angle1 = angle/diskAngleInterval;
			r1 = vertebraTemplate[z].diskContourRadius[angle1];
			if(r1>segInfo.cordRadius[z].y)
			{
				vertebraTemplate[z].diskContour[angle1].x = segInfo.diskCenter[z].x + dir.x*r1;
				vertebraTemplate[z].diskContour[angle1].y = segInfo.diskCenter[z].y + dir.y*r1;
				
				// visualize it
				if(debugMode && false)
				{
					maskA[z*sizexy+(int)vertebraTemplate[z].diskContour[angle1].y*sizex+(int)vertebraTemplate[z].diskContour[angle1].x] = 
						maskStruct.mask_falseDetection;
				}
			}
		}
	}	// for z
	
	// convert vertex in the vertebra body to different value to prevent self intersecting
	//
	vec2DynArray diskROI;
	intVec2DynArray scan;
	for(z=cordStartZ; z<=cordEndZ; z++)
	{
		k = z*sizexy;
		int angleNumber=360/diskAngleInterval;
		diskROI.SetSize(0);

		for(angle=0; angle<angleNumber; angle++)
		{
			if(vertebraTemplate[z].diskContour[angle].y>0 && vertebraTemplate[z].diskContour[angle].x>0
				&& vertebraTemplate[z].diskContour[angle].y<sizey && vertebraTemplate[z].diskContour[angle].x<sizex)
			{
				diskROI.Add(vertebraTemplate[z].diskContour[angle]);
			}
		}

		CIS_Algo_Contour_GetScanWhole(diskROI, scan);

		for(i=0; i<scan.GetSize(); i++)
		{
			x = scan[i].x;
			y = scan[i].y;
			k2 = k+y*sizex+x;

			if(maskA[k2]==maskStruct.mask_spinalCord)  maskA[k2]=100+maskStruct.mask_spinalCord;
			else if(imgA[k2]>segPara.boneThresh) maskA[k2]=100+maskStruct.mask_corticalBone;
			else maskA[k2]=100+maskStruct.mask_spongyBone;
		}
	}

	
	// locate and interpolated spinal process
	xt.SetSize(sizez);
	yt.SetSize(sizez);
	yt2.SetSize(sizez);
	for(z=0; z<sizez; z++)
	{
		xt[z] = z;
		yt[z] = -1;
	}

	for(z=0; z<sizez; z++)
	{
		vertebraTemplate[z].sprocessEnd = Vec2(-1, -1);
		vertebraTemplate[z].sprocessLength = vertebraTemplate[z].sprocessOrient = -1;
		
		for(i=0; i<sprocessCount; i++) 
		{
			vertebraTemplate[z].sprocessMedial[i] = Vec2(-1, -1);

			vertebraTemplate[z].sprocessLeft[i] = Vec2(-1, -1);
			vertebraTemplate[z].sprocessRight[i] = Vec2(-1, -1);

			vertebraTemplate[z].sprocessLeftWidth[i] = -1;
			vertebraTemplate[z].sprocessRightWidth[i] = -1;
		}
	}
	
	for(z=cordStartZ; z<=cordEndZ; z++)
	{
		if(sliceCord[z]==-1) continue;

		// search for the tip for spinal process
		Vec2 tip, tipc;
		int count;

		k = z*sizexy;
		tip.x = tip.y = -1;

		for(y=segInfo.cordCenter[z].y+segInfo.cordRadius[z].y*3; y<sizey && y<segInfo.cordCenter[z].y+segInfo.cordRadius[z].y*10; y++)
		{
			k2 = y*sizey+segInfo.cordCenter[z].x-segInfo.cordRadius[z].x;
			for(x=segInfo.cordCenter[z].x-segInfo.cordRadius[z].x; x<=segInfo.cordCenter[z].x+segInfo.cordRadius[z].x; x++, k2++)
			{
				if(maskA[k2+k]==maskStruct.mask_corticalBone)
				{
					tip.x = x;
					tip.y = y;
				}
			}
		}

		if(tip.x==-1) continue;

		// search a 2*7 neighborhood to locate the centroid
		count = 0;
		tipc = Vec2(0,0);
		for(y=tip.y-1; y<=tip.y; y++)
		{
			k2 = y*sizey+tip.x-3;
			for(x=tip.x-3; x<=tip.x+3; x++, k2++)
			{
				if(maskA[k2+k]==maskStruct.mask_corticalBone)
				{
					count++;
					tipc += Vec2(x,y);
				}
			}
		}	// for y
		
		tipc /= count;
		vertebraTemplate[z].sprocessEnd = tipc;
		
		tip = segInfo.cordCenter[z];

		vertebraTemplate[z].sprocessLength = (tipc-tip).len();

		tipc = tipc-tip;
		vertebraTemplate[z].sprocessOrient = atan2(tipc.y, tipc.x);

		// visualize
		if(debugMode && false)
		{ 
			maskA[z*sizexy+(int)vertebraTemplate[z].sprocessEnd.y*sizex+(int)vertebraTemplate[z].sprocessEnd.x] = 
				maskStruct.mask_boneMetasis;
		}
	}	// for z

	// interpolate the process tip (length and orient)
	for(z=0; z<sizez; z++)
	{
		xt[z] = z;
		yt[z] = -1;
	}
	count = 0;
	minInterpolate = -1;
	maxInterpolate = -1;
	for(z=cordStartZ; z<=cordEndZ; z++)
	{
		if(vertebraTemplate[z].sprocessLength!=-1)
		{
			count++;
			yt[z] = vertebraTemplate[z].sprocessLength;
			if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
			if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
		}
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2_x, bernstein_power, piece_size);

	minInterpolate2 = -1;
	maxInterpolate2 = -1;
	for(z=cordStartZ; z<=cordEndZ; z++)
	{
		if(vertebraTemplate[z].sprocessOrient!=-1)
		{
			yt[z] = vertebraTemplate[z].sprocessOrient;
			if(yt[z]<minInterpolate2 || minInterpolate2==-1) minInterpolate2=yt[z];
			if(yt[z]>maxInterpolate2 || maxInterpolate2==-1) maxInterpolate2=yt[z];
		}
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2_y, bernstein_power, piece_size);

	for(z=cordStartZ; z<=cordEndZ; z++)
	{
		Vec2 dir;
		if(vertebraTemplate[z].sprocessLength==-1 || vertebraTemplate[z].sprocessLength<yt2_x[z])
		{
			if(yt2_x[z]==-1) continue;
			if(yt2_x[z]<minInterpolate) yt2_x[z]=minInterpolate;
			if(yt2_x[z]>maxInterpolate) yt2_x[z]=maxInterpolate;

			if(yt2_y[z]==-1) continue;
			if(yt2_y[z]<minInterpolate2) yt2_y[z]=minInterpolate2;
			if(yt2_y[z]>maxInterpolate2) yt2_y[z]=maxInterpolate2;

			vertebraTemplate[z].sprocessLength = yt2_x[z];
			vertebraTemplate[z].sprocessOrient = yt2_y[z];

			dir.x = cos(yt2_y[z]);
			dir.y = sin(yt2_y[z]);

			vertebraTemplate[z].sprocessEnd = segInfo.cordCenter[z]+dir*yt2_x[z];

			// visualize
			if(debugMode && false)
				maskA[z*sizexy+(int)vertebraTemplate[z].sprocessEnd.y*sizex+(int)vertebraTemplate[z].sprocessEnd.x] = 
					maskStruct.mask_falseDetection;
		}
	}


	// search for left and right border of spinal process
	//
	for(z=cordStartZ; z<=cordEndZ; z++)
	{
		if(sliceCord[z]==-1) continue;

		if(vertebraTemplate[z].sprocessEnd.x==-1) continue;

		k = z*sizexy;
		Vec2 tip1, tip2, medial[sprocessCount], dir, dirx, curPos, nextPos, maxPos;
		double steps, max_r;
		int r, grad, maxGrad, k3;

		tip1 = vertebraTemplate[z].sprocessEnd;
		tip2.x = segInfo.cordCenter[z].x;
		tip2.y = segInfo.cordCenter[z].y+segInfo.cordRadius[z].y;

		dir = (tip2-tip1).normalize();
		dirx = dir.perp();
		if(dirx.x<0)
		{
			dirx.x = -dirx.x;
			dirx.y = -dirx.y;
		}

		vertebraTemplate[z].sprocessLength = (tip2-tip1).len();
		steps = vertebraTemplate[z].sprocessLength/(sprocessCount-1);

		medial[0] = tip1;
		medial[sprocessCount-1] = tip2;

		for(i=1; i<sprocessCount-1; i++) medial[i] = medial[i-1]+dir*steps;

		for(i=0; i<sprocessCount; i++) vertebraTemplate[z].sprocessMedial[i] = medial[i];

		// search left and right from medial to locate spinal process edge
		//
		for(i=0; i<sprocessCount; i++)
		{
			// left
			curPos = medial[i]-dirx;
			maxGrad = 0;
			for(r=1; r<segInfo.cordRadius[z].x; r++)
			{
				nextPos = curPos-dirx;
				k2 = k + (int)(curPos.y+0.5)*sizex + (int)(curPos.x+0.5);
				k3 = k + (int)(nextPos.y+0.5)*sizex + (int)(nextPos.x+0.5);

				if(maskA[k2]==maskStruct.mask_corticalBone && maskA[k3]==maskStruct.mask_body)
				{
					grad = imgA[k2]-imgA[k3];
					if(grad>maxGrad && grad>50)
					{
						maxGrad=grad;
						maxPos = curPos;
						max_r = r;
					}
				}
				curPos = nextPos;
			}	// for r

			if(maxGrad>0)
			{
				vertebraTemplate[z].sprocessLeft[i] = maxPos;
				vertebraTemplate[z].sprocessLeftWidth[i] = max_r;
			}
			

			// right
			curPos = medial[i]+dirx;
			maxGrad = 0;
			for(r=1; r<segInfo.cordRadius[z].x; r++)
			{
				nextPos = curPos+dirx;
				k2 = k + (int)(curPos.y+0.5)*sizex + (int)(curPos.x+0.5);
				k3 = k + (int)(nextPos.y+0.5)*sizex + (int)(nextPos.x+0.5);

				if(maskA[k2]==maskStruct.mask_corticalBone && maskA[k3]==maskStruct.mask_body)
				{
					grad = imgA[k2]-imgA[k3];
					if(grad>maxGrad && grad>50)
					{
						maxGrad=grad;
						maxPos = curPos;
						max_r = r;
					}
				}
				curPos = nextPos;
			}	// for r

			if(maxGrad>0)
			{
				vertebraTemplate[z].sprocessRight[i] = maxPos;
				vertebraTemplate[z].sprocessRightWidth[i] = max_r;
			}			
		}	// for i
	}	// for z

	// interpolate left and right border of spinal process	(by left and right width)
	//
	// left side
	for(i=0; i<sprocessCount; i++)
	{
		Vec2 dir, dirx;

		for(z=0; z<sizez; z++)
		{
			yt[z] = yt2_x[z] = -1;
		}
	
		count=0;
		minInterpolate = -1;
		maxInterpolate = -1;
		for(z=cordStartZ; z<=cordEndZ; z++)
		{
			if(vertebraTemplate[z].sprocessLeftWidth[i]!=-1)
			{
				count++;
				yt[z] = vertebraTemplate[z].sprocessLeftWidth[i];
				if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
				if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
			}
		}

		if(count<cordLength/2) continue;

		PiecewiseBernsteinSmoothing(xt, yt, yt2_x, bernstein_power, piece_size);

		for(z=cordStartZ; z<=cordEndZ; z++)
		{
			if(vertebraTemplate[z].sprocessMedial[i].x>0 && vertebraTemplate[z].sprocessLeftWidth[i]==-1 && yt2_x[z]>1)
			{
				if(yt2_x[z]<minInterpolate) yt2_x[z]=minInterpolate;
				if(yt2_x[z]>maxInterpolate) yt2_x[z]=maxInterpolate;
		
				dir.x = cos(vertebraTemplate[z].sprocessOrient);
				dir.y = sin(vertebraTemplate[z].sprocessOrient);
				dirx = dir.perp();

				if(dirx.x<0)
				{
					dirx.x = -dirx.x;
					dirx.y = -dirx.y;
				}
				
				vertebraTemplate[z].sprocessLeftWidth[i] = yt2_x[z];

				vertebraTemplate[z].sprocessLeft[i] = vertebraTemplate[z].sprocessMedial[i]-dirx*yt2_x[z];

				// visualize
				if(debugMode && false)
					maskA[z*sizexy+(int)vertebraTemplate[z].sprocessLeft[i].y*sizex+(int)vertebraTemplate[z].sprocessLeft[i].x] = 
						maskStruct.mask_falseDetection;
			}
		}
	}	// for i
	

	// right side
	for(i=0; i<sprocessCount; i++)
	{
		Vec2 dir, dirx;

		for(z=0; z<sizez; z++)
		{
			yt[z] = yt2_x[z] = -1;
		}
	
		count=0;
		minInterpolate = -1;
		maxInterpolate = -1;
		for(z=cordStartZ; z<=cordEndZ; z++)
		{
			if(vertebraTemplate[z].sprocessRightWidth[i]!=-1)
			{
				count++;
				yt[z] = vertebraTemplate[z].sprocessRightWidth[i];
				if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
				if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
			}
		}

		if(count<cordLength/2) continue;

		PiecewiseBernsteinSmoothing(xt, yt, yt2_x, bernstein_power, piece_size);

		for(z=cordStartZ; z<=cordEndZ; z++)
		{
			if(vertebraTemplate[z].sprocessMedial[i].x>0 && vertebraTemplate[z].sprocessRightWidth[i]==-1 && yt2_x[z]>1)
			{
				if(yt2_x[z]<minInterpolate) yt2_x[z]=minInterpolate;
				if(yt2_x[z]>maxInterpolate) yt2_x[z]=maxInterpolate;
				dir.x = cos(vertebraTemplate[z].sprocessOrient);
				dir.y = sin(vertebraTemplate[z].sprocessOrient);
				dirx = dir.perp();

				if(dirx.x<0)
				{
					dirx.x = -dirx.x;
					dirx.y = -dirx.y;
				}
				
				vertebraTemplate[z].sprocessRightWidth[i] = yt2_x[z];

				vertebraTemplate[z].sprocessRight[i] = vertebraTemplate[z].sprocessMedial[i]+dirx*yt2_x[z];

				// visualize
				if(debugMode && false)
					maskA[z*sizexy+(int)vertebraTemplate[z].sprocessRight[i].y*sizex+(int)vertebraTemplate[z].sprocessRight[i].x] = 
						maskStruct.mask_falseDetection;
			}
		}
	}	// for i
	
	
	
	// convert the process ROI
	//
	vec2DynArray processROI;

	for(z=cordStartZ; z<=cordEndZ; z++)
	{
		k = z*sizexy;
		processROI.SetSize(0);

		for(i=sprocessCount-1; i>=0; i--)
		{
			if(vertebraTemplate[z].sprocessLeft[i].x>0 && vertebraTemplate[z].sprocessLeft[i].x<sizex 
				&& vertebraTemplate[z].sprocessLeft[i].y>0 && vertebraTemplate[z].sprocessLeft[i].y<sizey)
			{
				processROI.Add(vertebraTemplate[z].sprocessLeft[i]);
			}
		}
		if(vertebraTemplate[z].sprocessEnd.x>0 && vertebraTemplate[z].sprocessEnd.x<sizex
			&& vertebraTemplate[z].sprocessEnd.y>0 && vertebraTemplate[z].sprocessEnd.y<sizey)
		{
			processROI.Add(vertebraTemplate[z].sprocessEnd);
		}
		for(i=0; i<sprocessCount; i++)
		{
			if(vertebraTemplate[z].sprocessRight[i].x>0 && vertebraTemplate[z].sprocessRight[i].x<sizex
				&& vertebraTemplate[z].sprocessRight[i].y>0 && vertebraTemplate[z].sprocessRight[i].y<sizey)
			{
				processROI.Add(vertebraTemplate[z].sprocessRight[i]);
			}
		}

		CIS_Algo_Contour_GetScanWhole(processROI, scan);

		for(i=0; i<scan.GetSize(); i++)
		{
			x = scan[i].x;
			y = scan[i].y;
			if(x<=0 || x>=sizex || y<=0 || y>=sizey) continue;

			k2 = k+y*sizex+x;

			if(maskA[k2]==maskStruct.mask_spinalCord) continue;

			if(imgA[k2]>segPara.boneThresh) maskA[k2]=200+maskStruct.mask_corticalBone;
			else maskA[k2]=200+maskStruct.mask_spongyBone;
		}
	}


	// locate left and right pedicle template
	//

	// left pedicle
	for(z=0; z<sizez; z++)
	{
		vertebraTemplate[z].leftPedicleEnd = Vec2(-1, -1);
		vertebraTemplate[z].leftPedicleLength = -1;
		
		for(i=0; i<sprocessCount; i++) 
		{
			vertebraTemplate[z].leftPedicleMedial[i] = Vec2(-1, -1);
			vertebraTemplate[z].leftPedicleUp[i] = Vec2(-1, -1);
			vertebraTemplate[z].leftPedicleDown[i] = Vec2(-1, -1);
		}
	}
	
	// dynamic searching for the medial model in pedicle
	for(z=cordStartZ; z<=cordEndZ; z++)
	{
		if(sliceCord[z]==-1) continue;

		if(segInfo.diskRadius[z].y<0) continue; 
		if(vertebraTemplate[z].sprocessEnd.x==-1) continue;

		// search for the tip of left pedicle
		Vec2 tip, tipc, dir, dirx, curPos, nextPos, maxPos;
		double search_r0, search_r1, search_a0, search_a1;
		double a, r, step_a;
		double  max_r, curScore, maxScore;
		short grad, maxGrad;
		int k3, k1, k4;

		k = z*sizexy;
		tip.x = tip.y = -1;

		// define the search radius range
		search_r0 = segInfo.diskRadius[z].y*0.5;
		search_r1 = segInfo.diskRadius[z].y*2.5;

		// define the search angle range
		dir = (vertebraTemplate[z].sprocessEnd-segInfo.cordCenter[z]).normalize();
		dirx = dir.perp();
		if(dirx.x>0)
		{
			dirx.x = -dirx.x; dirx.y=-dirx.y;
		}
		angle =atan2(dirx.y, dirx.x)*180/3.14;
		search_a0 = angle-45;
		search_a1 = angle+5;
		step_a = 50/search_r1;
		
		tip.x = segInfo.cordCenter[z].x;
		tip.y = segInfo.cordCenter[z].y;

		maxGrad = 0;
		max_r = 0;
		maxScore = 0;
		for(a=search_a0; a<=search_a1; a+=step_a)
		{
			dir.x = cos(a*3.14/180);
			dir.y = sin(a*3.14/180);

			curPos.x = tip.x+search_r0*dir.x;
			curPos.y = tip.y+search_r0*dir.y;
			for(r=search_r0; r<=search_r1; r++)
			{
				nextPos = curPos+dir;
				k2 = k + (int)(curPos.y+0.5)*sizex + (int)(curPos.x+0.5);
				k3 = k + (int)(nextPos.y+0.5)*sizex + (int)(nextPos.x+0.5);

				if(imgA[k2]>segPara.boneThresh && imgA[k3]<segPara.boneThresh && imgA[k3]>1000 && maskA[k2]<100 
					&& imgA[k3-segInfo.cordRadius[z].x*2]>700 && imgA[k3-segInfo.cordRadius[z].x*2*sizex]>700)
				{
					k1 = k + (int)(curPos.y-dir.y+0.5)*sizex + (int)(curPos.x-dir.x+0.5);
					k4 = k + (int)(nextPos.y+dir.y+0.5)*sizex + (int)(nextPos.x+dir.x+0.5);
					grad = (imgA[k1]+imgA[k2]-imgA[k3]-imgA[k4])/2;
	
					if(grad>50)
					{
						curScore = PedicleTemplateMatching(img3D, maskImg3D, z, curPos, segPara, segInfo, maskStruct, vertebraTemplate, -1, false);

						if(curScore>maxScore)
						{
							maxScore = curScore;
							maxGrad=grad;
							max_r = r;
							maxPos = curPos;
						}
					}
				}

				curPos += dir;
			}	// for r
		}	// for a

/*		if(maxScore>0)
		{
			vertebraTemplate[z].leftPedicleEnd = maxPos;
			vertebraTemplate[z].leftPedicleLength = max_r;
		}
		else	// shorter the search range and redo
		{
			// define the search radius range
			search_r0 = segInfo.diskRadius[z].y/2;
			search_r1 = segInfo.diskRadius[z].y*1.5;

			// define the search angle range
			dir = (vertebraTemplate[z].sprocessEnd-segInfo.cordCenter[z]).normalize();
			dirx = dir.perp();
			if(dirx.x>0)
			{
				dirx.x = -dirx.x; dirx.y=-dirx.y;
			}
			angle =atan2(dirx.y, dirx.x)*180/3.14;
			search_a0 = angle-45;
			search_a1 = angle+5;
			step_a = 50/search_r1;
		
			tip.x = segInfo.cordCenter[z].x;
			tip.y = segInfo.cordCenter[z].y;

			maxGrad = 0;
			max_r = 0;
			maxScore = 0;
			for(a=search_a0; a<=search_a1; a+=step_a)
			{
				dir.x = cos(a*3.14/180);
				dir.y = sin(a*3.14/180);

				curPos.x = tip.x+search_r0*dir.x;
				curPos.y = tip.y+search_r0*dir.y;
				for(r=search_r0; r<=search_r1; r++)
				{
					nextPos = curPos+dir;
					k2 = k + (int)(curPos.y+0.5)*sizex + (int)(curPos.x+0.5);
					k3 = k + (int)(nextPos.y+0.5)*sizex + (int)(nextPos.x+0.5);
	
					if(imgA[k2]>segPara.boneThresh && imgA[k3]<segPara.boneThresh && imgA[k3]>1000  && maskA[k2]<100
						&& imgA[k3-segInfo.cordRadius[z].x]>700 
						&& imgA[k3-segInfo.cordRadius[z].x*sizex]>700)
					{
						k1 = k + (int)(curPos.y-dir.y+0.5)*sizex + (int)(curPos.x-dir.x+0.5);
						k4 = k + (int)(nextPos.y+dir.y+0.5)*sizex + (int)(nextPos.x+dir.x+0.5);
						grad = (imgA[k1]+imgA[k2]-imgA[k3]-imgA[k4])/2;
	
						if(grad>50)
						{
							curScore = PedicleTemplateMatching(img3D, maskImg3D, z, curPos, segPara, segInfo, maskStruct, vertebraTemplate, -1, false);

							if(curScore>maxScore)
							{
								maxScore = curScore;
								maxGrad=grad;
								max_r = r;
								maxPos = curPos;
							}
						}
					}

					curPos += dir;
				}	// for r
			}	// for a
			
			if(maxScore>0)
			{
				vertebraTemplate[z].leftPedicleEnd = maxPos;
				vertebraTemplate[z].leftPedicleLength = max_r;
			}
		}
*/
		if(maxScore>0)
		{
			vertebraTemplate[z].leftPedicleEnd = maxPos;
			vertebraTemplate[z].leftPedicleLength = max_r;
		}

		if(vertebraTemplate[z].leftPedicleEnd.y>0)
		{
			if(debugMode && false)
			{
				// visualize it 
				k2 = k + (int)(vertebraTemplate[z].leftPedicleEnd.y+0.5)*sizex + (int)(vertebraTemplate[z].leftPedicleEnd.x+0.5);
				maskA[k2] = maskStruct.mask_spinalCord;
			}

			// fill the border using the tip
			PedicleTemplateMatching(img3D, maskImg3D, z, vertebraTemplate[z].leftPedicleEnd, segPara, segInfo, maskStruct, vertebraTemplate, -1, true);
		}
	}	// for z


	// right pedicle
	for(z=0; z<sizez; z++)
	{
		vertebraTemplate[z].rightPedicleEnd = Vec2(-1, -1);
		vertebraTemplate[z].rightPedicleLength = -1;
		
		for(i=0; i<sprocessCount; i++) 
		{
			vertebraTemplate[z].rightPedicleMedial[i] = Vec2(-1, -1);
			vertebraTemplate[z].rightPedicleUp[i] = Vec2(-1, -1);
			vertebraTemplate[z].rightPedicleDown[i] = Vec2(-1, -1);
		}
	}
	
	// dynamic searching for the medial model in pedicle
	for(z=cordStartZ; z<=cordEndZ; z++)
	{
		if(sliceCord[z]==-1) continue;

		if(segInfo.diskRadius[z].y<0) continue; 
		if(vertebraTemplate[z].sprocessEnd.x==-1) continue;

		// search for the tip of right pedicle
		Vec2 tip, tipc, dir, dirx, curPos, nextPos, maxPos;
		double search_r0, search_r1, search_a0, search_a1;
		double a, r, step_a;
		double  max_r, curScore, maxScore;
		short grad, maxGrad;
		int k1, k4, k3;

		k = z*sizexy;
		tip.x = tip.y = -1;

		// define the search radius range
		search_r0 = segInfo.diskRadius[z].y*0.5;
		search_r1 = segInfo.diskRadius[z].y*2.5;

		// define the search angle range
		dir = (vertebraTemplate[z].sprocessEnd-segInfo.cordCenter[z]).normalize();
		dirx = dir.perp();
		if(dirx.x<0)
		{
			dirx.x = -dirx.x; dirx.y=-dirx.y;
		}
		angle =atan2(dirx.y, dirx.x)*180/3.14;
		search_a0 = angle-5;
		search_a1 = angle+45;
		step_a = 50/search_r1;
		
		tip.x = segInfo.cordCenter[z].x;
		tip.y = segInfo.cordCenter[z].y;

		maxGrad = 0;
		max_r = 0;
		maxScore = 0;
		for(a=search_a0; a<=search_a1; a+=step_a)
		{
			dir.x = cos(a*3.14/180);
			dir.y = sin(a*3.14/180);

			curPos.x = tip.x+search_r0*dir.x;
			curPos.y = tip.y+search_r0*dir.y;
			for(r=search_r0; r<=search_r1; r++)
			{
				nextPos = curPos+dir;
				k2 = k + (int)(curPos.y+0.5)*sizex + (int)(curPos.x+0.5);
				k3 = k + (int)(nextPos.y+0.5)*sizex + (int)(nextPos.x+0.5);

				if(imgA[k2]>segPara.boneThresh && imgA[k3]<segPara.boneThresh && imgA[k3]>1000 && maskA[k2]<100 
					&& imgA[k3+segInfo.cordRadius[z].x*2]>700 && imgA[k3-segInfo.cordRadius[z].x*2*sizex]>700)
				{
					k1 = k + (int)(curPos.y-dir.y+0.5)*sizex + (int)(curPos.x-dir.x+0.5);
					k4 = k + (int)(nextPos.y+dir.y+0.5)*sizex + (int)(nextPos.x+dir.x+0.5);
					grad = (imgA[k1]+imgA[k2]-imgA[k3]-imgA[k4])/2;
	
					if(grad>50)
					{
						curScore = PedicleTemplateMatching(img3D, maskImg3D, z, curPos, segPara, segInfo, maskStruct, vertebraTemplate, 1, false);

						if(curScore>maxScore)
						{
							maxScore = curScore;
							maxGrad=grad;
							max_r = r;
							maxPos = curPos;
						}
					}
				}

				curPos += dir;
			}	// for r
		}	// for a

/*		if(maxScore>0)
		{
			vertebraTemplate[z].rightPedicleEnd = maxPos;
			vertebraTemplate[z].rightPedicleLength = max_r;
		}
		else	// shorter the search range and redo
		{
			// define the search radius range
			search_r0 = segInfo.diskRadius[z].y/2;
			search_r1 = segInfo.diskRadius[z].y*1.5;

			// define the search angle range
			dir = (vertebraTemplate[z].sprocessEnd-segInfo.cordCenter[z]).normalize();
			dirx = dir.perp();
			if(dirx.x<0)
			{
				dirx.x = -dirx.x; dirx.y=-dirx.y;
			}
			angle =atan2(dirx.y, dirx.x)*180/3.14;
			search_a0 = angle-5;
			search_a1 = angle+45;
			step_a = 50/search_r1;
		
			tip.x = segInfo.cordCenter[z].x;
			tip.y = segInfo.cordCenter[z].y;

			maxGrad = 0;
			max_r = 0;
			maxScore = 0;
			for(a=search_a0; a<=search_a1; a+=step_a)
			{
				dir.x = cos(a*3.14/180);
				dir.y = sin(a*3.14/180);

				curPos.x = tip.x+search_r0*dir.x;
				curPos.y = tip.y+search_r0*dir.y;
				for(r=search_r0; r<=search_r1; r++)
				{
					nextPos = curPos+dir;
					k2 = k + (int)(curPos.y+0.5)*sizex + (int)(curPos.x+0.5);
					k3 = k + (int)(nextPos.y+0.5)*sizex + (int)(nextPos.x+0.5);
	
					if(imgA[k2]>segPara.boneThresh && imgA[k3]<segPara.boneThresh && imgA[k3]>1000  && maskA[k2]<100
						&& imgA[k3+segInfo.cordRadius[z].x]>700 
						&& imgA[k3-segInfo.cordRadius[z].x*sizex]>700)
					{
						k1 = k + (int)(curPos.y-dir.y+0.5)*sizex + (int)(curPos.x-dir.x+0.5);
						k4 = k + (int)(nextPos.y+dir.y+0.5)*sizex + (int)(nextPos.x+dir.x+0.5);
						grad = (imgA[k1]+imgA[k2]-imgA[k3]-imgA[k4])/2;
	
						if(grad>50)
						{
							curScore = PedicleTemplateMatching(img3D, maskImg3D, z, curPos, segPara, segInfo, maskStruct, vertebraTemplate, 1, false);

							if(curScore>maxScore)
							{
								maxScore = curScore;
								maxGrad=grad;
								max_r = r;
								maxPos = curPos;
							}
						}
					}

					curPos += dir;
				}	// for r
			}	// for a
			
			if(maxScore>0)
			{
				vertebraTemplate[z].rightPedicleEnd = maxPos;
				vertebraTemplate[z].rightPedicleLength = max_r;
			}
		}
*/
		if(maxScore>0)
		{
			vertebraTemplate[z].rightPedicleEnd = maxPos;
			vertebraTemplate[z].rightPedicleLength = max_r;
		}

		if(vertebraTemplate[z].rightPedicleEnd.y>0)
		{
			// visualize it 
			if(debugMode && false)
			{
				k2 = k + (int)(vertebraTemplate[z].rightPedicleEnd.y+0.5)*sizex + (int)(vertebraTemplate[z].rightPedicleEnd.x+0.5);
				maskA[k2] = maskStruct.mask_spinalCord;
			}

			// fill the border using the tip
			PedicleTemplateMatching(img3D, maskImg3D, z, vertebraTemplate[z].rightPedicleEnd, segPara, segInfo, maskStruct, vertebraTemplate, 1, true);
		}
	}	// for z


	// convert the pedicle ROI
	//
	vec2DynArray pedicleROI;

	// left pedicle
	for(z=cordStartZ; z<=cordEndZ; z++)
	{
		k = z*sizexy;
		pedicleROI.SetSize(0);

		for(i=pedicleCount-1; i>=0; i--)
		{
			if(vertebraTemplate[z].leftPedicleUp[i].x>0 && vertebraTemplate[z].leftPedicleUp[i].x<sizex
				&& vertebraTemplate[z].leftPedicleUp[i].y>0 && vertebraTemplate[z].leftPedicleUp[i].y<sizey)
			{
				pedicleROI.Add(vertebraTemplate[z].leftPedicleUp[i]);
			}
		}
		if(vertebraTemplate[z].leftPedicleEnd.x>0 && vertebraTemplate[z].leftPedicleEnd.x<sizex
			&& vertebraTemplate[z].leftPedicleEnd.y>0 && vertebraTemplate[z].leftPedicleEnd.y<sizey)
		{
			pedicleROI.Add(vertebraTemplate[z].leftPedicleEnd);
		}
		for(i=0; i<pedicleCount; i++)
		{
			if(vertebraTemplate[z].leftPedicleDown[i].x>0 && vertebraTemplate[z].leftPedicleDown[i].x<sizex
				&& vertebraTemplate[z].leftPedicleDown[i].y>0 && vertebraTemplate[z].leftPedicleDown[i].y<sizey)
			{
				pedicleROI.Add(vertebraTemplate[z].leftPedicleDown[i]);
			}
		}


		CIS_Algo_Contour_GetScanWhole(pedicleROI, scan);

		for(i=0; i<scan.GetSize(); i++)
		{
			x = scan[i].x;
			y = scan[i].y;
			if(x<=0 || x>=sizex || y<=0 || y>=sizey) continue;
			k2 = k+y*sizex+x;

			if(maskA[k2]==maskStruct.mask_spinalCord) continue;

			if(imgA[k2]>segPara.boneThresh) maskA[k2]=300+maskStruct.mask_corticalBone;
			else maskA[k2]=300+maskStruct.mask_spongyBone;
		}
	}

	// right pedicle
	for(z=cordStartZ; z<=cordEndZ; z++)
	{
		k = z*sizexy;
		pedicleROI.SetSize(0);

		for(i=pedicleCount-1; i>=0; i--)
		{
			if(vertebraTemplate[z].rightPedicleUp[i].x>0 && vertebraTemplate[z].rightPedicleUp[i].x<sizex
				&& vertebraTemplate[z].rightPedicleUp[i].y>0 && vertebraTemplate[z].rightPedicleUp[i].y<sizey)
			{
				pedicleROI.Add(vertebraTemplate[z].rightPedicleUp[i]);
			}
		}
		if(vertebraTemplate[z].rightPedicleEnd.x>0 && vertebraTemplate[z].rightPedicleEnd.x<sizex
			&& vertebraTemplate[z].rightPedicleEnd.y>0 && vertebraTemplate[z].rightPedicleEnd.y<sizey)
		{
			pedicleROI.Add(vertebraTemplate[z].rightPedicleEnd);
		}
		for(i=0; i<pedicleCount; i++)
		{
			if(vertebraTemplate[z].rightPedicleDown[i].x>0 && vertebraTemplate[z].rightPedicleDown[i].x<sizex
				&& vertebraTemplate[z].rightPedicleDown[i].y>0 && vertebraTemplate[z].rightPedicleDown[i].y<sizey)
			{
				pedicleROI.Add(vertebraTemplate[z].rightPedicleDown[i]);
			}
		}


		CIS_Algo_Contour_GetScanWhole(pedicleROI, scan);

		for(i=0; i<scan.GetSize(); i++)
		{
			x = scan[i].x;
			y = scan[i].y;
			if(x<=0 || x>=sizex || y<=0 || y>=sizey) continue;
			k2 = k+y*sizex+x;

			if(maskA[k2]==maskStruct.mask_spinalCord) continue;

			if(imgA[k2]>segPara.boneThresh) maskA[k2]=400+maskStruct.mask_corticalBone;
			else maskA[k2]=400+maskStruct.mask_spongyBone;
		}
	}


	// get the last portion of the template, the region around cannel
	vec2DynArray cannelROI;
	for(z=cordStartZ; z<=cordEndZ; z++)
	{
		k = z*sizexy;
		cannelROI.SetSize(0);

		// for disk
		for(i=0; i<diskAngleCount; i++)
		{
			if(vertebraTemplate[z].diskContour[i].x>0 && vertebraTemplate[z].diskContour[i].x<sizex
				&& vertebraTemplate[z].diskContour[i].y>0 && vertebraTemplate[z].diskContour[i].y<sizey)
			{
				cannelROI.Add(vertebraTemplate[z].diskContour[i]);
				break;
			}
		}
		for(i=diskAngleCount-1; i>=0; i--)
		{
			if(vertebraTemplate[z].diskContour[i].x>0 && vertebraTemplate[z].diskContour[i].x<sizex
				&& vertebraTemplate[z].diskContour[i].y>0 && vertebraTemplate[z].diskContour[i].y<sizey)
			{
				cannelROI.Add(vertebraTemplate[z].diskContour[i]);
				break;
			}
		}

		// left pedicle
		for(i=pedicleCount-1; i>=0; i--)
		{
			if(vertebraTemplate[z].leftPedicleUp[i].x>0 && vertebraTemplate[z].leftPedicleUp[i].x<sizex
				&& vertebraTemplate[z].leftPedicleUp[i].y>0 && vertebraTemplate[z].leftPedicleUp[i].y<sizey)
			{
				cannelROI.Add(vertebraTemplate[z].leftPedicleUp[i]);
				break;
			}
		}
		for(i=pedicleCount-1; i>=0; i--)
		{
			if(vertebraTemplate[z].leftPedicleDown[i].x>0 && vertebraTemplate[z].leftPedicleDown[i].x<sizex
				&& vertebraTemplate[z].leftPedicleDown[i].y>0 && vertebraTemplate[z].leftPedicleDown[i].y<sizey)
			{
				cannelROI.Add(vertebraTemplate[z].leftPedicleDown[i]);
				break;
			}
		}

		// process
		for(i=sprocessCount-1; i>=0; i--)
		{
			if(vertebraTemplate[z].sprocessLeft[i].x>0 && vertebraTemplate[z].sprocessLeft[i].x<sizex
				&& vertebraTemplate[z].sprocessLeft[i].y>0 && vertebraTemplate[z].sprocessLeft[i].y<sizey)
			{
				cannelROI.Add(vertebraTemplate[z].sprocessLeft[i]);
				break;
			}
		}
		for(i=sprocessCount-1; i>=0; i--)
		{
			if(vertebraTemplate[z].sprocessRight[i].x>0 && vertebraTemplate[z].sprocessRight[i].x<sizex
				&& vertebraTemplate[z].sprocessRight[i].y>0 && vertebraTemplate[z].sprocessRight[i].y<sizey)
			{
				cannelROI.Add(vertebraTemplate[z].sprocessRight[i]);
				break;
			}
		}

		// right pedicle
		for(i=pedicleCount-1; i>=0; i--)
		{
			if(vertebraTemplate[z].rightPedicleDown[i].x>0 && vertebraTemplate[z].rightPedicleDown[i].x<sizex
				&& vertebraTemplate[z].rightPedicleDown[i].y>0 && vertebraTemplate[z].rightPedicleDown[i].y<sizey)
			{
				cannelROI.Add(vertebraTemplate[z].rightPedicleDown[i]);
				break;
			}
		}
		for(i=pedicleCount-1; i>=0; i--)
		{
			if(vertebraTemplate[z].rightPedicleUp[i].x>0 && vertebraTemplate[z].rightPedicleUp[i].x<sizex
				&& vertebraTemplate[z].rightPedicleUp[i].y>0 && vertebraTemplate[z].rightPedicleUp[i].y<sizey)
			{
				cannelROI.Add(vertebraTemplate[z].rightPedicleUp[i]);
				break;
			}
		}

		CIS_Algo_Contour_GetScanWhole(cannelROI, scan);

		for(i=0; i<scan.GetSize(); i++)
		{
			x = scan[i].x;
			y = scan[i].y;
			if(x<=0 || x>=sizex || y<=0 || y>=sizey) continue;
			k2 = k+y*sizex+x;

			if(maskA[k2]==maskStruct.mask_spinalCord || maskA[k2]==maskStruct.mask_spongyBone || maskA[k2]==maskStruct.mask_corticalBone) 
				maskA[k2]+=500;
		}
	}


	// get overall template information
	for(z=0; z<sizez; z++)
	{
		vertebraTemplate[z].diskCenter.x = segInfo.diskCenter[z].x;
		vertebraTemplate[z].diskCenter.y = segInfo.diskCenter[z].y;
		vertebraTemplate[z].diskRadius = segInfo.diskRadius[z].y;

		vertebraTemplate[z].cordCenter.x = segInfo.cordCenter[z].x;
		vertebraTemplate[z].cordCenter.y = segInfo.cordCenter[z].y;
		vertebraTemplate[z].cordRadius.x = segInfo.cordRadius[z].x;
		vertebraTemplate[z].cordRadius.y = segInfo.cordRadius[z].y;
	}

	// define the bounding box
	for(z=0; z<sizez; z++)
	{
		if(segInfo.cordCenter[z].x==-1) continue;

		// initialize
		segInfo.spineBound1[z] = IntVec2(sizex,sizey);
		segInfo.spineBound2[z] = IntVec2(-1,-1);
		segInfo.diskBound1[z] = IntVec2(sizex,sizey);
		segInfo.diskBound2[z] = IntVec2(-1,-1);
		segInfo.cordBound1[z] = IntVec2(sizex,sizey);
		segInfo.cordBound2[z] = IntVec2(-1,-1);
		segInfo.sprocessBound1[z] = IntVec2(sizex,sizey);
		segInfo.sprocessBound2[z] = IntVec2(-1,-1);

		k2 = z*sizexy;
		for(y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k2++)
			{
				if(maskA[k2]>100)
				{
					if(x<segInfo.spineBound1[z].x) segInfo.spineBound1[z].x=x;
					if(y<segInfo.spineBound1[z].y) segInfo.spineBound1[z].y=y;
					if(x>segInfo.spineBound2[z].x) segInfo.spineBound2[z].x=x;
					if(y>segInfo.spineBound2[z].y) segInfo.spineBound2[z].y=y;

					if(maskA[k2]<200)	// disk region
					{
						if(x<segInfo.diskBound1[z].x) segInfo.diskBound1[z].x=x;
						if(y<segInfo.diskBound1[z].y) segInfo.diskBound1[z].y=y;
						if(x>segInfo.diskBound2[z].x) segInfo.diskBound2[z].x=x;
						if(y>segInfo.diskBound2[z].y) segInfo.diskBound2[z].y=y;
					}
					else if(maskA[k2]<300)	// process
					{
						if(x<segInfo.sprocessBound1[z].x) segInfo.sprocessBound1[z].x=x;
						if(y<segInfo.sprocessBound1[z].y) segInfo.sprocessBound1[z].y=y;
						if(x>segInfo.sprocessBound2[z].x) segInfo.sprocessBound2[z].x=x;
						if(y>segInfo.sprocessBound2[z].y) segInfo.sprocessBound2[z].y=y;
					}
					else
					{
						if(maskA[k2]==500+maskStruct.mask_spinalCord)
						{
							if(x<segInfo.cordBound1[z].x) segInfo.cordBound1[z].x=x;
							if(y<segInfo.cordBound1[z].y) segInfo.cordBound1[z].y=y;
							if(x>segInfo.cordBound2[z].x) segInfo.cordBound2[z].x=x;
							if(y>segInfo.cordBound2[z].y) segInfo.cordBound2[z].y=y;
						}
					}
				}
			}
		}	// for y

		// put some margin on the bound
		segInfo.spineBound1[z].x -=8;
		segInfo.spineBound1[z].y -=8;
		segInfo.spineBound2[z].x +=8;
		segInfo.spineBound2[z].y +=8;

		segInfo.diskBound1[z].x -=4;
		segInfo.diskBound1[z].y -=4;
		segInfo.diskBound2[z].x +=4;
		segInfo.diskBound2[z].y +=4;

		// leave margin to image border
		if(segInfo.spineBound1[z].x<8) segInfo.spineBound1[z].x=8;
		if(segInfo.spineBound1[z].y<8) segInfo.spineBound1[z].y=8;
		if(segInfo.spineBound2[z].x>sizex-8) segInfo.spineBound2[z].x=sizex-8;
		if(segInfo.spineBound2[z].y>sizey-8) segInfo.spineBound2[z].y=sizey-8;

		if(segInfo.diskBound1[z].x<8) segInfo.diskBound1[z].x=8;
		if(segInfo.diskBound1[z].y<8) segInfo.diskBound1[z].y=8;
		if(segInfo.diskBound2[z].x>sizex-8) segInfo.diskBound2[z].x=sizex-8;
		if(segInfo.diskBound2[z].y>sizey-8) segInfo.diskBound2[z].y=sizey-8;

	}	// for z


/*	// temp turned off 

	//trim thin and small spongybone region on the boundary
	//
	for(z=cordStartZ; z<=cordEndZ; z++)
	{
		count=0;
		k = z*sizexy;
		for(k2=0; k2<sizexy; k2++)
		{
			if(maskA[k+k2]==maskStruct.mask_spongyBone) 
			{
				binArray2D[k2]=1;
				count++;
			}
			else binArray2D[k2]=0;
		}

		if(count==0) continue;
		NIH_Algo_Blob_Labelling(binImg2D, blobImg2D, blobRank2D, blobCentroid2D, false);

		blobStatus.SetSize(blobRank2D.GetSize());

		for(i=0; i<blobRank2D.GetSize(); i++) blobStatus[i]=0;

		// get those blobs on the border
		int b;
		for(k2=0; k2<sizexy; k2++)
		{
			b = blobArray2D[k2];
			if(b!=0)
			{
				if(maskA[k+k2-1]==maskStruct.mask_body || maskA[k+k2+1]==maskStruct.mask_body ||
					maskA[k+k2-sizex]==maskStruct.mask_body || maskA[k+k2+sizex]==maskStruct.mask_body) blobStatus[b-1]=1;
				if(maskA[k+k2-1]==maskStruct.mask_otherBone || maskA[k+k2+1]==maskStruct.mask_otherBone ||
					maskA[k+k2-sizex]==maskStruct.mask_otherBone || maskA[k+k2+sizex]==maskStruct.mask_otherBone) blobStatus[b-1]=1;
			}
		}

		// open the blobs and eliminate thin blobs
		CIS_IPA_Erode(binImg2D, 1, false);
		CIS_IPA_Dilate(binImg2D, 1, false);

		for(k2=0; k2<sizexy; k2++)
		{
			b = blobArray2D[k2];
			if(b!=0 && binArray2D[k2]==0 && blobStatus[b-1]==1)
			{
				if(blobStatus[b-1]==1) maskA[k+k2]=maskStruct.mask_body;
			}
		}
		
	}
	

	// cut un-connected regions
	//
	for(z=cordStartZ; z<=cordEndZ; z++)
	{
		count=0;
		k = z*sizexy;
		for(k2=0; k2<sizexy; k2++)
		{
			if(maskA[k+k2]==maskStruct.mask_spongyBone || maskA[k+k2]==maskStruct.mask_corticalBone 
				|| maskA[k+k2]==maskStruct.mask_spinalCord) 
			{
				binArray2D[k2]=1;
				count++;
			}
			else binArray2D[k2]=0;
		}

		if(count==0) continue;

		// open the blobs and eliminate thin blobs
		CIS_IPA_Erode(binImg2D, 1, false);
		CIS_IPA_Dilate(binImg2D, 1, false);

		NIH_Algo_Blob_Labelling(binImg2D, blobImg2D, blobRank2D, blobCentroid2D, false);

		int largest_b, largest_rank;
		int b;
		largest_rank = 0;
		for(b=0; b<blobRank2D.GetSize(); b++)
		{
			if(blobRank2D[b]>largest_rank)
			{
				largest_rank = blobRank2D[b];
				largest_b = b+1;
			}
		}


		for(y=0, k2=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k2++)
			{
				if(x>=segInfo.cordBound1[z].x && x<=segInfo.cordBound2[z].x && y>segInfo.cordBound1[z].y) continue;

				if(maskA[k+k2]==maskStruct.mask_spongyBone || maskA[k+k2]==maskStruct.mask_corticalBone )
				{
					b = blobArray2D[k2];

					if(b!=largest_b) maskA[k+k2]=maskStruct.mask_otherBone;
				}
			}
		}
		
	}
	

	// recover some cortical bones
	// and close small holes inside
	for(z=cordStartZ; z<=cordEndZ; z++)
	{
		count=0;
		k = z*sizexy;
		for(k2=0; k2<sizexy; k2++)
		{
			if(maskA[k+k2]==maskStruct.mask_spongyBone || maskA[k+k2]==maskStruct.mask_corticalBone 
				|| maskA[k+k2]==maskStruct.mask_spinalCord) 
			{
				binArray2D[k2]=1;
				count++;
			}
			else binArray2D[k2]=0;
		}

		if(count==0) continue;

		// open the blobs and eliminate thin blobs
		CIS_IPA_Dilate(binImg2D, 2, false);

		for(y=0, k2=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k2++)
			{
				if(binArray2D[k2]==1 && maskA[k+k2]!=maskStruct.mask_corticalBone  && imgA[k+k2]>segPara.boneThresh)
				{
					maskA[k+k2]=maskStruct.mask_corticalBone;
				}
			}
		}
		
		CIS_IPA_Erode(binImg2D, 2, false);

		for(y=0, k2=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k2++)
			{
				if(binArray2D[k2]==1 && maskA[k+k2]==maskStruct.mask_body)
				{
					maskA[k+k2]=maskStruct.mask_spongyBone;
				}
			}
		}
		
	}
	
*/	

	// convert the mask back
	for(k=0; k<sizexyz; k++)
	{
		// convert those outside of template to other bone
		if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_spongyBone || maskA[k]==maskStruct.mask_spinalCord) 
			maskA[k]=maskStruct.mask_otherBone;

		if(maskA[k]==maskStruct.mask_falseDetection) maskA[k] = maskStruct.mask_spongyBone;
		if(maskA[k]==maskStruct.mask_boneMetasis) maskA[k] = maskStruct.mask_spongyBone;

		if(maskA[k]>500) maskA[k]-=500;	// spinal cord
		else if(maskA[k]>400)			// right pedicle
		{
			if(imgA[k]>1100) maskA[k]-=400;
		}
		else if(maskA[k]>300)			// left pedicle
		{
			if(imgA[k]>1100) maskA[k]-=300;
		}	
		else if(maskA[k]>200)			// spinal process
		{
			if(imgA[k]>1100) maskA[k]-=200;
		}
		else if(maskA[k]>100) maskA[k]-=100;	// disk
	}

	// refinement on each 2D slice
	for(z=cordStartZ; z<=cordEndZ; z++)
	{
		if(sliceCord[z]==-1) continue;
		
		int count;
		k = z*sizexy;
		// clean the buffer
		for(k2=0; k2<sizexy; k2++) binArray2D[k2]=0;

		// step 1. add back corticle bone
		// for every picel > boneThresh and within diskBound, add it back
		count = 0;
		for(y=segInfo.diskBound1[z].y-4; y<=segInfo.diskBound2[z].y; y++)
		{
			for(x=segInfo.diskBound1[z].x-4,k2=y*sizex+x; x<=segInfo.diskBound2[z].x+4; x++, k2++)
			{
				if(maskA[k+k2]==maskStruct.mask_otherBone && imgA[k+k2]>segPara.boneThresh) maskA[k+k2]=maskStruct.mask_corticalBone;
			}
		}	// for y

		// step 2: do a connected component analysis
		// every spine pixel inside spineBound should be connected expect for those inside sProcessBound
		for(y=segInfo.spineBound1[z].y-4; y<=segInfo.spineBound2[z].y; y++)
		{
			for(x=segInfo.spineBound1[z].x-4,k2=y*sizex+x; x<=segInfo.spineBound2[z].x+4; x++, k2++)
			{
				if(maskA[k+k2]==maskStruct.mask_spongyBone || maskA[k+k2]==maskStruct.mask_corticalBone || maskA[k+k2]==maskStruct.mask_spinalCord) 
				{
					binArray2D[k2]=1;
					count++;
				}
			}
		}	// for y
		if(count==0) continue;
		NIH_Algo_Blob_Labelling(binImg2D, blobImg2D, blobRank2D, blobCentroid2D, false);

		int largest_b, largest_rank;
		int b;
		largest_rank = 0;
		for(b=0; b<blobRank2D.GetSize(); b++)
		{
			if(blobRank2D[b]>largest_rank)
			{
				largest_rank = blobRank2D[b];
				largest_b = b+1;
			}
		}

		for(y=segInfo.spineBound1[z].y; y<=segInfo.spineBound2[z].y; y++)
		{
			for(x=segInfo.spineBound1[z].x, k2=y*sizex+x; x<=segInfo.spineBound2[z].x; x++, k2++)
			{
				if(blobArray2D[k2]>0 && blobArray2D[k2]!=largest_b) 
				{
					if(x<segInfo.sprocessBound1[z].x || x>segInfo.sprocessBound2[z].x ||
						y<segInfo.sprocessBound1[z].y || y>segInfo.sprocessBound2[z].y)
						maskA[k+k2]=maskStruct.mask_otherBone;
				}
			}
		}	// for y

		// step 3: do a close operation to close hole and small gaps
		for(y=segInfo.spineBound1[z].y; y<=segInfo.spineBound2[z].y; y++)
		{
			for(x=segInfo.spineBound1[z].x, k2=y*sizex+x; x<=segInfo.spineBound2[z].x; x++, k2++)
			{
				if(maskA[k+k2]==maskStruct.mask_otherBone) binArray2D[k2]=0;
			}
		}	// for y
		CIS_IPA_Dilate(binImg2D, segPara.closeCorticalHoleIteration, false);
		CIS_IPA_Erode(binImg2D, segPara.closeCorticalHoleIteration, false);
		for(y=segInfo.spineBound1[z].y; y<=segInfo.spineBound2[z].y; y++)
		{
			k2 = y*sizex+segInfo.spineBound1[z].x;
			for(x=segInfo.spineBound1[z].x; x<=segInfo.spineBound2[z].x; x++, k2++)
			{
				if(binArray2D[k2]==1 && 
					maskA[k+k2]!=maskStruct.mask_corticalBone
					&& maskA[k+k2]!=maskStruct.mask_spinalCord) maskA[k+k2]=maskStruct.mask_spongyBone;
			}
		}	// for y

	}


	delete binImg2D;
	delete blobImg2D;
	free(cord);

	return CIS_OK;
}


// detect spinal cord and set up bounding box for each part of the vertebra (vetrabra body and spinal process)
// old method, without vertebra template
//
int NIH_SpineSegmentation_DetectSpinalCord_old(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo, bool debugMode)
{
	int x, y, z, k;
	int sizex, sizey, sizez, sizexy, sizexyz;
	short *imgA, *maskA;

	sizex = img3D->Num_Cols();
	sizey = img3D->Num_Rows();
	sizez = img3D->Num_Levels();

	imgA = img3D->GetArray();
	maskA = maskImg3D->GetArray();

	sizexyz = sizex*sizey*sizez;
	sizexy = sizex*sizey;

	CIS_Array_Image2D_short *binImg2D, *blobImg2D;
	short *binArray2D, *blobArray2D;
	intDynArray blobRank2D;
	intVec2DynArray blobCentroid2D;
	int largestBlob, largestSize;
	intDynArray blobStatus;

	int count, k2;
	int stx, sty, edx, edy;

	binImg2D = new CIS_Array_Image2D_short(sizex, sizey);
	blobImg2D = new CIS_Array_Image2D_short(sizex, sizey);
	binArray2D = binImg2D->GetArray();
	blobArray2D = blobImg2D->GetArray();

	typedef struct 
	{
		short index;
		short type;
		short slice;
		short centerx, centery;
		short bx0, bx1, by0, by1;	// bounding box
		short size;
		short link_up, link_down;
		float overlap_up, overlap_down;
		float score;
	} CORD_NODE;


	CORD_NODE *cord;
	int max_cord, num_cord;
	short mask_tmpbase=1000;

	max_cord = sizez*10;
	num_cord = 0;
	cord = (CORD_NODE *)malloc(max_cord*sizeof(CORD_NODE));

	segInfo.bound1.z = 0;
	// first pass to get the candidate cord locations
	//
	for(z=0; z<sizez; z++)
	{
		// first define a bounding box
		k = z*sizexy;
		count =0;
		stx = sizex-3;
		sty = sizey-3;
		edx = 3;
		edy = 3;
		for(k2=0, y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++, k2++)
			{
				if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_spongyBone)
				{
					if(x<stx) stx=x;
					if(x>edx) edx=x;
					if(y<sty) sty=y;
					if(y>sty) edy=y;
					count++;
				}
			}
		}

		if(count==0) continue;
		if(segInfo.bound1.z==0) segInfo.bound1.z=z;
		segInfo.bound2.z = z;

		// initialize the 2D image
		for(k2=0; k2<sizexy; k2++) binArray2D[k2]=0;

		// first pass elinimate external space
		stx-=6; sty-=6; edx+=6; edy+=6;
		k = z*sizexy;
		count=0;
		for(y=sty; y<=edy; y++)
		{
			k2 = y*sizex+stx;
			for(x=stx; x<=edx; x++, k2++)
			{
				if(maskA[k+k2]!=maskStruct.mask_corticalBone && maskA[k+k2]!=maskStruct.mask_spongyBone)
				{
					binArray2D[k2]=1;
					count ++;
				}
			}
		}	// for y

		if(count==0) continue;

		// separate cord and external space
		CIS_IPA_Erode(binImg2D, segPara.closeCordIteration, false);
		CIS_IPA_Dilate(binImg2D, segPara.closeCordIteration, false);

		NIH_Algo_Blob_Labelling(binImg2D, blobImg2D, blobRank2D, blobCentroid2D, false);

		largestSize = 0;
		largestBlob = -1;
		for(k2=0; k2<blobRank2D.GetSize(); k2++) 
		{
			if(blobRank2D[k2]>largestSize) 
			{
				largestSize=blobRank2D[k2];
				largestBlob = k2+1;
			}
		}

		// remove largest blob (the external region) or high intensity region
		k = z*sizexy;
		count=0;
		for(y=sty; y<=edy; y++)
		{
			k2 = y*sizex+stx;
			for(x=stx; x<=edx; x++, k2++)
			{
				if(binArray2D[k2]==1)
				{
					if(blobArray2D[k2]==largestBlob || imgA[k+k2]>segPara.cordIntensityThresh)
					{
						binArray2D[k2]=0;

						if(blobArray2D[k2]!=largestBlob) 
						{
							maskA[k+k2]=maskStruct.mask_spongyBone;	// interior region

							// dilate a pixel to cover holes caused by erosion
							if(maskA[k+k2-1]==maskStruct.mask_body && binArray2D[k2-1]==0) maskA[k+k2-1]=maskStruct.mask_spongyBone;
							if(maskA[k+k2+1]==maskStruct.mask_body && binArray2D[k2+1]==0) maskA[k+k2+1]=maskStruct.mask_spongyBone;
							if(maskA[k+k2-sizex]==maskStruct.mask_body && binArray2D[k2-sizex]==0) maskA[k+k2-sizex]=maskStruct.mask_spongyBone;
							if(maskA[k+k2+sizex]==maskStruct.mask_body && binArray2D[k2+sizex]==0) maskA[k+k2+sizex]=maskStruct.mask_spongyBone;
						}
					}	// if blobArray
				}	// if binArray
			}	// for x
		}	// for y
		
		// second pass to get candidate cords
		NIH_Algo_Blob_Labelling(binImg2D, blobImg2D, blobRank2D, blobCentroid2D, false);
		largestBlob = -1;
		blobStatus.SetSize(blobRank2D.GetSize());

		for(k2=0; k2<blobStatus.GetSize(); k2++) 
		{
			if(k2==largestBlob-1) blobStatus[k2]=-1;	// largest blob should be the space outside the vertebra
			else
			{
				if(blobRank2D[k2]<segPara.cordSizeThresh) blobStatus[k2]=-1;
				else
				{
					if(num_cord>=max_cord)
					{
						max_cord = max_cord*2;
						cord = (CORD_NODE *)realloc(cord, max_cord*sizeof(CORD_NODE));
					}

					cord[num_cord].index = mask_tmpbase+num_cord;
					cord[num_cord].type = 0;
					cord[num_cord].slice = z;
					cord[num_cord].centerx = blobCentroid2D[k2].x;
					cord[num_cord].centery = blobCentroid2D[k2].y;
					cord[num_cord].bx0 = blobCentroid2D[k2].x;
					cord[num_cord].bx1 = blobCentroid2D[k2].x;
					cord[num_cord].by0 = blobCentroid2D[k2].y;
					cord[num_cord].by1 = blobCentroid2D[k2].y;
					cord[num_cord].size = blobRank2D[k2];
					cord[num_cord].link_up = -1;
					cord[num_cord].link_down = -1;
					cord[num_cord].overlap_up = 0;
					cord[num_cord].overlap_down = 0;
					cord[num_cord].score = 0;

					blobStatus[k2] = mask_tmpbase+num_cord;
					num_cord++;
				}
			}
		}	// for k2
		 
		// assign a temp mask value for cord candidate pixels (different candidate has different mask value)
		k = z*sizexy;
		for(k2=0, y=0; y<sizey; y++) 
		{
			for(x=0; x<sizex; x++, k++, k2++)
			{
				if(binArray2D[k2]==1)  
				{
					if(blobStatus[blobArray2D[k2]-1]!=-1)
					{
						maskA[k]=blobStatus[blobArray2D[k2]-1];
					}
					else maskA[k]=maskStruct.mask_spongyBone;
				}	// if binArray
			}
		}
	}	// for z

	int org_num_cord = num_cord;
	// dilate the cord to fill gaps between slices
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		// initialize the 2D image
		for(k2=0; k2<sizexy; k2++) binArray2D[k2]=0;

		k = z*sizexy;
		count =0;
		for(k2=0, y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++, k2++)
			{
				if(maskA[k]==maskStruct.mask_body && imgA[k]<segPara.cordIntensityThresh && 
					((maskA[k-sizexy]>=mask_tmpbase && maskA[k-sizexy]<mask_tmpbase+org_num_cord) 
					|| (maskA[k+sizexy]>=mask_tmpbase && maskA[k+sizexy]<mask_tmpbase+org_num_cord)))
				{
					binArray2D[k2] = 1;
					count++;
				}
			}
		}

		if(count==0) continue;
		NIH_Algo_Blob_Labelling(binImg2D, blobImg2D, blobRank2D, blobCentroid2D, false);

		blobStatus.SetSize(blobRank2D.GetSize());

		// add new candidate nodes
		for(k2=0; k2<blobStatus.GetSize(); k2++) 
		{
			if(blobRank2D[k2]<segPara.cordSizeThresh) blobStatus[k2]=-1;
			else
			{
				if(num_cord>=max_cord)
				{
					max_cord = max_cord*2;
					cord = (CORD_NODE *)realloc(cord, max_cord*sizeof(CORD_NODE));
				}

				cord[num_cord].index = mask_tmpbase+num_cord;
				cord[num_cord].type = 1;
				cord[num_cord].slice = z;
				cord[num_cord].centerx = blobCentroid2D[k2].x;
				cord[num_cord].centery = blobCentroid2D[k2].y;
				cord[num_cord].bx0 = blobCentroid2D[k2].x;
				cord[num_cord].bx1 = blobCentroid2D[k2].x;
				cord[num_cord].by0 = blobCentroid2D[k2].y;
				cord[num_cord].by1 = blobCentroid2D[k2].y;
				cord[num_cord].size = blobRank2D[k2];
				cord[num_cord].link_up = -1;
				cord[num_cord].link_down = -1;
				cord[num_cord].overlap_up = 0;
				cord[num_cord].overlap_down = 0;
				cord[num_cord].score = 0;

				blobStatus[k2] = mask_tmpbase+num_cord;
				num_cord++;
			}
		}	// for k2
		
		k = z*sizexy;
		for(k2=0, y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++, k2++)
			{
				if(binArray2D[k2]==1) 
				{
					if(blobStatus[blobArray2D[k2]-1]!=-1)
					{
						maskA[k]=blobStatus[blobArray2D[k2]-1];
					}
					else maskA[k]=maskStruct.mask_spongyBone;
				}	// if binArray
			}
		}
	
	}	// for z


	// compute the overlap and score of each node candidate
	intDynArray list_up, list_down, over_up, over_down;
	int n, tn, i;
	short mask_tmplast;
	float largestOverlap, c_overlap;
	short overIndex;

	mask_tmplast = mask_tmpbase+num_cord;

	for(n=0; n<num_cord; n++)
	{
		z = cord[n].slice;
		k = z*sizexy;
		list_up.SetSize(0);
		list_down.SetSize(0);
		over_up.SetSize(0);
		over_down.SetSize(0);

		for(y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++)
			{
				if(maskA[k]==mask_tmpbase+n)
				{
					// check overlap with previous slice
					tn = maskA[k-sizexy]-mask_tmpbase;
					if(tn>=0 && tn<num_cord)
					{
						for(i=0; i<list_up.GetSize(); i++)
						{
							if(list_up[i]==tn)
							{
								break;
							}
						}
						if(i>=list_up.GetSize())
						{
							list_up.Add(tn);
							over_up.Add(1);
						}
						else
						{
							over_up[i]++;
						}
					}	// if tn

					// check overlap with next slice
					tn = maskA[k+sizexy]-mask_tmpbase;
					if(tn>=0 && tn<num_cord)
					{
						for(i=0; i<list_down.GetSize(); i++)
						{
							if(list_down[i]==tn)
							{
								break;
							}
						}
						if(i>=list_down.GetSize())
						{
							list_down.Add(tn);
							over_down.Add(1);
						}
						else
						{
							over_down[i]++;
						}
					}	// if tn
				}	// if maskA
			}	// for x
		}	// for y

		// found the largest overlap node in previous and next slice
		largestOverlap = 0;
		overIndex = -1;
		for(i=0; i<list_up.GetSize(); i++)
		{
			c_overlap = (float)over_up[i]/(float)(cord[n].size+cord[list_up[i]].size-over_up[i]);
			if(c_overlap>largestOverlap)
			{
				largestOverlap = c_overlap;
				overIndex = list_up[i];
			}
		}
		cord[n].link_up = overIndex;
		cord[n].overlap_up = largestOverlap;

		largestOverlap = 0;
		overIndex = -1;
		for(i=0; i<list_down.GetSize(); i++)
		{
			c_overlap = (float)over_down[i]/(float)(cord[n].size+cord[list_down[i]].size-over_down[i]);
			if(c_overlap>largestOverlap)
			{
				largestOverlap = c_overlap;
				overIndex = list_down[i];
			}
		}
		cord[n].link_down = overIndex;
		cord[n].overlap_down = largestOverlap;

	}	// for n

	// find the longest path in the graph
	// compute the score
	int cn, nn;
	float c_score, a_score;;
	for(n=0; n<num_cord; n++)
	{
		a_score = 1;
		// going up first
		cn = n;
		c_score = 1;
		do {
			nn = cord[cn].link_up;
			if(nn==-1) break;
//			c_score = c_score*cord[cn].overlap_up;
			c_score = cord[cn].overlap_up;
			a_score += c_score;
			cn = nn;
		} while(1);

		// going down then
		cn = n;
		c_score = 1;
		do {
			nn = cord[cn].link_down;
			if(nn==-1) break;
//			c_score = c_score*cord[cn].overlap_down;
			c_score = cord[cn].overlap_down;
			a_score += c_score;
			cn = nn;
		} while(1);

		cord[n].score = a_score;
	}

	// adjust the link based on scores
	for(n=0; n<num_cord; n++)
	{
		nn = cord[n].link_down;
		if(nn!=-1 && cord[nn].link_up!=n)
		{
			if(cord[nn].link_up==-1) cord[nn].link_up=n;
			else
			{
				if(cord[n].score>cord[cord[nn].link_up].score) 
				{
					cord[nn].link_up=n;	
					cord[nn].overlap_up=cord[n].overlap_down;	
				}
			}
		}

		nn = cord[n].link_up;
		if(nn!=-1 && cord[nn].link_down!=n)
		{
			if(cord[nn].link_down==-1) cord[nn].link_down=n;
			else
			{
				if(cord[n].score>cord[cord[nn].link_down].score) 
				{
					cord[nn].link_down=n;	
					cord[nn].overlap_down=cord[n].overlap_up;	
				}
			}
		}	
	}

	// recompute the score
	for(n=0; n<num_cord; n++)
	{
		a_score = 1;
		// going up first
		cn = n;
		c_score = 1;
		do {
			nn = cord[cn].link_up;
			if(nn==-1) break;
			c_score = cord[cn].overlap_up;
			a_score += c_score;
			cn = nn;
		} while(1);

		// going down then
		cn = n;
		c_score = 1;
		do {
			nn = cord[cn].link_down;
			if(nn==-1) break;
			c_score = cord[cn].overlap_down;
			a_score += c_score;
			cn = nn;
		} while(1);

		if(a_score>=cord[n].score) cord[n].score = a_score;
		else
		{
			cn=n;	// something is wrong
		}
	}

	// only keep the cord with largest score on each slice
	intDynArray sliceCord;
	doubleDynArray sliceCordScore;
	sliceCord.SetSize(sizez);
	sliceCordScore.SetSize(sizez);
	for(z=0; z<sizez; z++)
	{
		sliceCord[z] = -1;
		sliceCordScore[z] = 0;
	}

	blobStatus.SetSize(num_cord);
	for(n=0; n<num_cord; n++)
	{
		z = cord[n].slice;
		if(cord[n].score>sliceCordScore[z])
		{
			sliceCord[z] = n;
			sliceCordScore[z] = cord[n].score;
		}

		blobStatus[n] = 0;
	}

	// correct the spinal cord if it is not a continuous chain
	int longest_st, longest_ed, chain_st, chain_ed, all_st, all_ed;
	longest_st = longest_ed = 0;
	z = segInfo.bound1.z;
	for(z=segInfo.bound1.z; z<segInfo.bound2.z; z++)
	{
		if(sliceCord[z]!=-1) break;
	}
	all_st = segInfo.bound2.z;
	all_ed = segInfo.bound1.z;
	while(z<segInfo.bound2.z)
	{
		chain_st = z;
		chain_ed = z;
		while(sliceCord[z]!=-1 && sliceCord[z+1]!=-1 && cord[sliceCord[z]].link_down==sliceCord[z+1])
		{
			if(z<all_st) all_st=z;

			z++;
			chain_ed = z;

			if(z>all_ed) all_ed=z;
		}

		if(chain_ed-chain_st>longest_ed-longest_st)
		{
			longest_st = chain_st;
			longest_ed = chain_ed;
		}

		z++;
	}	// while z
	
	// correct the chain
	int pz, cz, pn, ck, pk;
	int lastValidCord;
	// going upwards
	pz = longest_st;
	lastValidCord = longest_st;
	if(longest_st>all_st)
	{
		pn = sliceCord[longest_st];
		pz = longest_st;
		cz = longest_st-1;
		over_up.SetSize(num_cord);

		while(cz>all_st)
		{
			pk=pz*sizexy;
			ck=cz*sizexy;

			for(n=0; n<num_cord; n++)  over_up[n]=0;

			for(y=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, ck++, pk++)
				{
					if(maskA[pk]==mask_tmpbase+pn)
					{
						tn = maskA[ck]-mask_tmpbase;
						if(tn>=0 && tn<num_cord)
						{
							over_up[tn]++;
						}
					}	// if maskA[pk];
				}	// for x
			}	// for y

			// get the largest overlap
			largestOverlap = 0;
			overIndex = -1;
			for(n=0; n<num_cord; n++)
			{
				if(over_up[n]>largestOverlap)
				{
					largestOverlap = over_up[n];
					overIndex = n;
				}
			}

			if(overIndex>=0)
			{
				sliceCord[cz] = overIndex;
				cord[overIndex].link_down = pn;
				cord[overIndex].overlap_down = largestOverlap;
				cord[pn].link_up = overIndex;
				cord[pn].overlap_up = largestOverlap;	// overlap not right

				pn = overIndex;
				pz = cz;
				cz = pz-1;

				while(cord[pn].link_up>=0)
				{
					cn = cord[pn].link_up;
					sliceCord[cz] = cn;
					pn = cn;
					pz = cz;
					cz = pz-1;
				}

				lastValidCord = cz;
			}
			else //
			{
				// add a new node using the overlap
				largestOverlap = 0;
				k = cz*sizexy;
				for(k2=0, y=0; y<sizey; y++)
				{
					for(x=0; x<sizex; x++, k++, k2++)
					{
						if(maskA[k+sizexy]== mask_tmpbase+pn // && maskA[k]==maskStruct.mask_body
							&& imgA[k]<segPara.cordIntensityThresh) // also need to test the intensity
						{
							maskA[k]=mask_tmpbase+num_cord;
							largestOverlap ++;
						}
					}	// if maskA
				}	// for k2

				if(largestOverlap<segPara.cordSizeThresh) break;

				if(num_cord>=max_cord)
				{
					max_cord = max_cord*2;
					cord = (CORD_NODE *)realloc(cord, max_cord*sizeof(CORD_NODE));
				}

				cord[num_cord].index = mask_tmpbase+num_cord;
				cord[num_cord].type = 0;
				cord[num_cord].slice = cz;
				cord[num_cord].centerx = cord[pn].centerx;
				cord[num_cord].centery = cord[pn].centery;
				cord[num_cord].bx0 = cord[pn].bx0;
				cord[num_cord].bx1 = cord[pn].bx1;
				cord[num_cord].by0 = cord[pn].by0;
				cord[num_cord].by1 = cord[pn].by1;
				cord[num_cord].size = largestOverlap;
				cord[num_cord].link_up = -1;
				cord[num_cord].link_down = pn;
				cord[num_cord].overlap_up = 0;
				cord[num_cord].overlap_down = largestOverlap;	// not right
				cord[num_cord].score = cord[pn].score;

				cord[pn].link_up = num_cord;
				cord[pn].overlap_up = largestOverlap;

				sliceCord[cz] = num_cord;

				over_up.Add(0);

				pn = num_cord;
				pz = cz;
				cz = pz-1;

				num_cord++;
			}	// else
		} // while z
	}	// if longest_st

	// dis-qualify the cord from lastValidCord to all_st
	for(z=all_st; z<lastValidCord; z++) sliceCord[z] = -1;

	if(segInfo.bound1.z<pz) segInfo.bound1.z=pz;

	// trace downwards
	pz = longest_ed;
	lastValidCord = longest_ed;
	if(longest_ed<all_ed)
	{
		pn = sliceCord[longest_ed];
		pz = longest_ed;
		cz = pz+1;
		over_down.SetSize(num_cord);


		while(cz<all_ed)
		{
			pk=pz*sizexy;
			ck=cz*sizexy;

			for(n=0; n<num_cord; n++)  over_down[n]=0;

			for(y=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, ck++, pk++)
				{
					if(maskA[pk]==mask_tmpbase+pn)
					{
						tn = maskA[ck]-mask_tmpbase;
						if(tn>=0 && tn<num_cord)
						{
							over_down[tn]++;
						}
					}	// if maskA[pk];
				}	// for x
			}	// for y

			// get the largest overlap
			largestOverlap = 0;
			overIndex = -1;
			for(n=0; n<num_cord; n++)
			{
				if(over_down[n]>largestOverlap)
				{
					largestOverlap = over_down[n];
					overIndex = n;
				}
			}

			if(overIndex>=0)
			{
				sliceCord[cz] = overIndex;
				cord[overIndex].link_up = pn;
				cord[overIndex].overlap_up = largestOverlap;
				cord[pn].link_down = overIndex;
				cord[pn].overlap_down = largestOverlap;

				pn = overIndex;
				pz = cz;
				cz = pz+1;

				while(cord[pn].link_down>=0)
				{
					cn = cord[pn].link_down;
					sliceCord[cz] = cn;
					pn = cn;
					pz = cz;
					cz = pz+1;
				}

				lastValidCord = cz;
			}
			else //
			{
				// add a new node using the overlap
				largestOverlap = 0;
				k = cz*sizexy;
				for(k2=0, y=0; y<sizey; y++)
				{
					for(x=0; x<sizex; x++, k++, k2++)
					{
						if(maskA[k-sizexy]== mask_tmpbase+pn // && maskA[k]==maskStruct.mask_body
							&& imgA[k]<segPara.cordIntensityThresh) 
						{
							maskA[k]=mask_tmpbase+num_cord;
							largestOverlap ++;
						}
					}	// if maskA
				}	// for k2

				if(largestOverlap<segPara.cordSizeThresh) break;

				if(num_cord>=max_cord)
				{
					max_cord = max_cord*2;
					cord = (CORD_NODE *)realloc(cord, max_cord*sizeof(CORD_NODE));
				}

				cord[num_cord].index = mask_tmpbase+num_cord;
				cord[num_cord].type = 0;
				cord[num_cord].slice = cz;
				cord[num_cord].centerx = cord[pn].centerx;
				cord[num_cord].centery = cord[pn].centery;
				cord[num_cord].bx0 = cord[pn].bx0;
				cord[num_cord].bx1 = cord[pn].bx1;
				cord[num_cord].by0 = cord[pn].by0;
				cord[num_cord].by1 = cord[pn].by1;
				cord[num_cord].size = largestOverlap;
				cord[num_cord].link_up = pn;
				cord[num_cord].link_down = -1;
				cord[num_cord].overlap_up = largestOverlap;
				cord[num_cord].overlap_down = 0;
				cord[num_cord].score = cord[pn].score;

				cord[pn].link_down = num_cord;
				cord[pn].overlap_down = largestOverlap;

				sliceCord[cz] = num_cord;

				over_down.Add(0);

				pn = num_cord;
				pz = cz;
				cz = pz+1;

				num_cord++;
			}	// else
		} // while z
	}	// if longest_ed

	if(segInfo.bound2.z>pz) segInfo.bound2.z=pz;

	// dis-qualify the cord from lastValidCord to all_ed
	for(z=lastValidCord; z<all_ed; z++) sliceCord[z] = -1;

	// assign the cord
	blobStatus.SetSize(num_cord);
	for(n=0; n<num_cord; n++) blobStatus[n]=0;
	for(z=0; z<sizez; z++)
	{
		if(sliceCord[z]>=0) blobStatus[sliceCord[z]]=1;
	}

	// assign the mask
	for(k=0; k<sizexyz; k++)
	{
		if(maskA[k]>=mask_tmpbase && maskA[k]<mask_tmpbase+num_cord) 
		{
			if(blobStatus[maskA[k]-mask_tmpbase]==1) maskA[k] = maskStruct.mask_spinalCord;
			else maskA[k] = maskStruct.mask_spongyBone;	// convert the holes
		}
	}

	// dilate the cord to fill the gap
	IntVec2 seed;
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		k = z*sizexy;
		count = 0;
		for(k2=0; k2<sizexy; k2++, k++) 
		{
			if(maskA[k]==maskStruct.mask_spinalCord) 
			{
				count++;
				binArray2D[k2]=1;
			}
			else binArray2D[k2]=0;
		}

		if(count==0) continue;

		CIS_IPA_Dilate(binImg2D, segPara.closeCordIteration, true);

		seed.x = cord[sliceCord[z]].centerx;
		seed.y = cord[sliceCord[z]].centery;
		NIH_Algo_Region_Growth_2D(binImg2D, seed, 4, blobImg2D, -1, -1, -1, -1);

		k=z*sizexy;
		for(k2=0; k2<sizexy; k2++, k++) 
		{
			if(blobArray2D[k2]==1 && (maskA[k]==maskStruct.mask_spongyBone || maskA[k]==maskStruct.mask_body)) 
				maskA[k]=maskStruct.mask_spinalCord;
		}
		
		// find the bounding box
		k=z*sizexy;
		for(y=0, k2=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k2++, k++)
			{
				if(maskA[k]==maskStruct.mask_spinalCord)
				{
					if(x<cord[sliceCord[z]].bx0) cord[sliceCord[z]].bx0 = x;
					if(x>cord[sliceCord[z]].bx1) cord[sliceCord[z]].bx1 = x;
					if(y<cord[sliceCord[z]].by0) cord[sliceCord[z]].by0 = y;
					if(y>cord[sliceCord[z]].by1) cord[sliceCord[z]].by1 = y;
				}
			}
		}
		
	}	// for z

	// remove excessive bones on each 2D slices
	double dist, minDistUp, minDistDown;
	int blobUp, blobDown;

	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;

		k = z*sizexy;
		count =0;
		for(k2=0, y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++, k2++)
			{
				if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_spongyBone)
				{
					binArray2D[k2] = 1;
					count++;
				}
				else binArray2D[k2]=0;
			}
		}

		if(count==0) continue;

		NIH_Algo_Blob_Labelling(binImg2D, blobImg2D, blobRank2D, blobCentroid2D, false);

		blobStatus.SetSize(blobRank2D.GetSize());
		blobUp = blobDown = -1;
		minDistUp = minDistDown = (sizex+sizey)*(sizex+sizey);
		for(k2=0; k2<blobRank2D.GetSize(); k2++) 
		{
			blobStatus[k2]=0;

			if(blobRank2D[k2]<segPara.cordSizeThresh*2) 
			{
				blobStatus[k2]=1;
				continue;
			}
			
			dist = (blobCentroid2D[k2].x-cord[sliceCord[z]].centerx)*(blobCentroid2D[k2].x-cord[sliceCord[z]].centerx)+
					(blobCentroid2D[k2].y-cord[sliceCord[z]].centery)*(blobCentroid2D[k2].y-cord[sliceCord[z]].centery);

			if(blobCentroid2D[k2].y<=cord[sliceCord[z]].centery && dist<minDistUp && blobRank2D[k2]>segPara.cordSizeThresh*10)
			{
				minDistUp = dist;
				blobUp = k2;
			}
			else if(blobCentroid2D[k2].y>cord[sliceCord[z]].centery && dist<minDistDown)
			{
				minDistDown = dist;
				blobDown = k2;
			}
		}

		if(blobUp>=0) 
		{
			blobStatus[blobUp]=1;
			// every one below it should also be included
			for(k2=0; k2<blobRank2D.GetSize(); k2++) 
			{
				if(blobStatus[k2]==0 && blobCentroid2D[k2].y<=cord[sliceCord[z]].centery)
				{
					dist = (blobCentroid2D[k2].x-cord[sliceCord[z]].centerx)*(blobCentroid2D[k2].x-cord[sliceCord[z]].centerx)+
							(blobCentroid2D[k2].y-cord[sliceCord[z]].centery)*(blobCentroid2D[k2].y-cord[sliceCord[z]].centery);
					if(dist<minDistUp) blobStatus[k2]=1;
				}
			}

		}
		else
		{
			// just chose the closest one
			for(k2=0; k2<blobRank2D.GetSize(); k2++) 
			{
				if(blobStatus[k2]==1) continue;
				dist = (blobCentroid2D[k2].x-cord[sliceCord[z]].centerx)*(blobCentroid2D[k2].x-cord[sliceCord[z]].centerx)+
					(blobCentroid2D[k2].y-cord[sliceCord[z]].centery)*(blobCentroid2D[k2].y-cord[sliceCord[z]].centery);

				if(blobCentroid2D[k2].y<=cord[sliceCord[z]].centery && dist<minDistUp && blobRank2D[k2]>segPara.cordSizeThresh*10)
				{
					minDistUp = dist;
					blobUp = k2;
				}
			}
			if(blobUp>=0) blobStatus[blobUp]=1;
		}

		if(blobDown>=0) blobStatus[blobDown]=1;

		// remove other bones
		k = z*sizexy;
		for(k2=0, y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++, k2++)
			{
				if(binArray2D[k2]==1) 
				{
					if(blobStatus[blobArray2D[k2]-1]==0)
					{
						maskA[k]=maskStruct.mask_otherBone;
					}
				}	// if binArray
			}
		}
	}	// for z

	// find the bounding box for each slice
	IntVec2 c_center, c_bound1, c_bound2;
	int len_y1, len_y2, len_x1, len_x2, largest_len;
	
	segInfo.spineBound1.SetSize(sizez);
	segInfo.spineBound2.SetSize(sizez);
	segInfo.diskBound1.SetSize(sizez);
	segInfo.diskBound2.SetSize(sizez);
	segInfo.cordBound1.SetSize(sizez);
	segInfo.cordBound2.SetSize(sizez);
	segInfo.sprocessBound1.SetSize(sizez);
	segInfo.sprocessBound2.SetSize(sizez);

	for(z=0; z<sizez; z++)
	{
		segInfo.spineBound1[z] = IntVec2(-1,-1);
		segInfo.spineBound2[z] = IntVec2(-1,-1);
		segInfo.diskBound1[z] = IntVec2(-1,-1);
		segInfo.diskBound2[z] = IntVec2(-1,-1);
		segInfo.cordBound1[z] = IntVec2(-1,-1);
		segInfo.cordBound2[z] = IntVec2(-1,-1);
		segInfo.sprocessBound1[z] = IntVec2(-1,-1);
		segInfo.sprocessBound2[z] = IntVec2(-1,-1);
	}

	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;

		segInfo.cordBound1[z].x = cord[sliceCord[z]].bx0;
		segInfo.cordBound1[z].y = cord[sliceCord[z]].by0;
		segInfo.cordBound2[z].x = cord[sliceCord[z]].bx1;
		segInfo.cordBound2[z].y = cord[sliceCord[z]].by1;

		c_center.x = cord[sliceCord[z]].centerx;
		c_center.y = cord[sliceCord[z]].centery;

		y = c_center.y;
		k = z*sizexy+y*sizex+c_center.x;
		for(x=c_center.x; x<sizex; x++, k++)
		{
			if(maskA[k]==maskStruct.mask_body) break;
		}
		c_bound2.x = x;

		k = z*sizexy+y*sizex+c_center.x;
		for(x=c_center.x; x>0; x--, k--)
		{
			if(maskA[k]==maskStruct.mask_body) break;
		}
		c_bound1.x = x;

		// locate the upper bound
		k = z*sizexy;
		c_bound1.y = sizey-1;
		int gap;
		for(x=c_bound1.x; x<c_bound2.x; x++)
		{
			gap = 0;	// keep a 10-pixel gaps
			for(y=c_center.y; y>10 && gap<10; y--)
			{
				if(maskA[k+y*sizex+x]!=maskStruct.mask_corticalBone && maskA[k+y*sizex+x]!=maskStruct.mask_spongyBone) 
				{
					gap ++;
				}
				else
				{
					gap = 0;
				}
			}

			if(y<c_bound1.y) c_bound1.y = y;
		}

		// locate the lower bound
		k = z*sizexy;
		c_bound2.y = 0;
		for(x=c_bound1.x; x<c_bound2.x; x++)
		{
			for(y=c_center.y; y<sizey; y++)
			{
				if(maskA[k+y*sizex+x]==maskStruct.mask_body) break;
			}

			if(y>c_bound2.y) c_bound2.y = y;
		}

		len_x1 = c_center.x-c_bound1.x;
		len_x2 = c_bound2.x-c_center.x;
		len_y1 = c_center.y-c_bound1.y;
		len_y2 = c_bound2.y-c_center.y;

		// adjust the len, using the longest of vertebra body and process
		if(len_y1>len_y2) largest_len=len_y1;
		else largest_len=len_y2;

		largest_len = largest_len*3/2;	// add some margins

		// make the spine bbox to be a square centered at the spinal cord
		//
		c_bound1.x = c_center.x-largest_len;
		c_bound1.y = c_center.y-largest_len;
		c_bound2.x = c_center.x+largest_len;
		c_bound2.y = c_center.y+largest_len;

		// adjust the bbox if on the border of the image
		if(c_bound1.x<0) c_bound1.x=0;
		if(c_bound1.y<0) c_bound1.y=0;
		if(c_bound2.x>sizex-1) c_bound2.x=sizex-1;
		if(c_bound2.y>sizey-1) c_bound2.y=sizey-1;

		segInfo.spineBound1[z] = c_bound1;
		segInfo.spineBound2[z] = c_bound2;

		// get the initial bbox of disk, assume it is four times as wide as the spinal canel
		segInfo.diskBound1[z].x = c_center.x - (segInfo.cordBound2[z].x-segInfo.cordBound1[z].x)*2;
		segInfo.diskBound2[z].x = c_center.x + (segInfo.cordBound2[z].x-segInfo.cordBound1[z].x)*2;
		segInfo.diskBound1[z].y = c_center.y - len_y1;
		segInfo.diskBound2[z].y = c_center.y;
		if(segInfo.diskBound1[z].x<0) segInfo.diskBound1[z].x=0;
		if(segInfo.diskBound2[z].x>sizex-1) segInfo.diskBound2[z].x=sizex-1;

	}

	// cut all the bones outside the bounding box 
	// and recover some bones inside the box
	//
	int k1;
	for(z=0; z<sizez; z++)
	{
		// if no bounding box, which means spine is not located, then all bones in the slices should be other bones
		if(segInfo.spineBound1[z].y==-1) 
		{
			k=z*sizexy;
			for(k2=0; k2<sizexy; k2++, k++)
			{
				if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_spongyBone || 
					maskA[k]==maskStruct.mask_spinalCord) maskA[k]=maskStruct.mask_otherBone;
			}
			continue;
		}

		k=z*sizexy;

		for(y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++)
			{
				// convert bone outside bounding box to other bone
				if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_spongyBone)
				{
					if(y<segInfo.spineBound1[z].y || y>segInfo.spineBound2[z].y || x<segInfo.spineBound1[z].x || x>segInfo.spineBound2[z].x)
						maskA[k] = maskStruct.mask_otherBone;
				}
/*				// convert bones inside bounding box to spongy bone
				else if(maskA[k]==maskStruct.mask_otherBone) 
				{
					if(y>=segInfo.spineBound1[z].y && y<=segInfo.spineBound2[z].y && x>=segInfo.spineBound1[z].x && x<=segInfo.spineBound2[z].x)
						maskA[k] = maskStruct.mask_spongyBone;					
				}
*/
			}
		}	// for y

		// restore some missing disk based on interpolation
		k=z*sizexy;
		for(y=0; y<=segInfo.diskBound2[z].y; y++)
		{
			k1 = k+y*sizex+segInfo.diskBound1[z].x;
			for(x=segInfo.diskBound1[z].x; x<=segInfo.diskBound2[z].x; x++, k1++)
			{
				if(maskA[k1]==maskStruct.mask_body && 
					(maskA[k1+sizexy]!=maskStruct.mask_body && maskA[k1+sizexy]!=maskStruct.mask_otherBone) &&
					(maskA[k1-sizexy]!=maskStruct.mask_body && maskA[k1-sizexy]!=maskStruct.mask_otherBone))
					maskA[k1] = maskStruct.mask_spongyBone;
			}
		}
	}	// for z

	// eliminate ribs
	//
	int count1;
	IntVec2 cordCenter, leftSeed, rightSeed;
	bool foundLeft, foundRight, flagStop;
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;

		// first find out if lung region is within the bounding box, 
		// if yes, usually require to remove ribs
		count = count1 = 0;
		k=z*sizexy;
		for(y=segInfo.spineBound1[z].y; y<=segInfo.spineBound2[z].y; y++)
		{
			k1 = k+y*sizex+segInfo.spineBound1[z].x;
			for(x=segInfo.spineBound1[z].x; x<=segInfo.spineBound2[z].x; x++, k1++)
			{
				if(maskA[k1]==maskStruct.mask_air) count++;
				if(imgA[k1]>1300)	//maskA[k1]==maskStruct.mask_corticalBone && 
				{
					count1 ++;
					binArray2D[k1-k] = 1;
				}
				else binArray2D[k1-k] = 0;
			}
		}

		if(count==0 || count1==0) continue;

//		CIS_IPA_Erode(binImg2D, 1, true);
//		CIS_IPA_Dilate(binImg2D, 1, true);
		
		cordCenter = (segInfo.cordBound1[z]+segInfo.cordBound2[z])/2;
		// locate the seed point for left and right ribs

		// left side
		y = cordCenter.y;
		foundLeft = false;
		for(y=segInfo.cordBound1[z].y-5; y<=segInfo.cordBound2[z].y && !foundLeft;y++)
		{
			for(x=segInfo.diskBound1[z].x; x>segInfo.spineBound1[z].x; x--)
			{
				if(binArray2D[y*sizex+x]==1)
				{
					foundLeft = true;
					leftSeed = IntVec2(x, y);
				}
				if(maskA[k+y*sizex+x]==maskStruct.mask_air) break;
			}
		}

		if(foundLeft)
		{
			NIH_Algo_Region_Growth_2D(binImg2D, leftSeed, 8, blobImg2D, -1, -1, -1, -1);

			CIS_IPA_Dilate(blobImg2D, 1, true);

			for(y=0, k1=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k1++)
				{
					if(blobArray2D[k1]==1 && x<segInfo.cordBound1[z].x-5)
					{
						maskA[k1+k] = maskStruct.mask_rib;
					}
				}
			}

			// eliminate spongy bone inside and between ribs
			for(y=segInfo.spineBound1[z].y; y<=segInfo.spineBound2[z].y; y++)
			{
				k1 = k+y*sizex+segInfo.spineBound1[z].x;
				flagStop=false;
				for(x=segInfo.spineBound1[z].x; x<=segInfo.cordBound1[z].x-5 && !flagStop; x++, k1++)
				{
					if(maskA[k1]==maskStruct.mask_rib)
					{
						for(;x<=segInfo.cordBound1[z].x-5; x++, k1++)
						{
							if(binArray2D[k1-k]==1 && maskA[k1]!=maskStruct.mask_rib)
							{
								flagStop = true;
								break;
							}

							if(maskA[k1]==maskStruct.mask_spongyBone || maskA[k1]==maskStruct.mask_corticalBone) 
								maskA[k1]=maskStruct.mask_rib;
						}
						flagStop = true;
					}	// if maskA
				}	// for x
			}	// for y

		}	// if foundLeft

		// right side
		y = cordCenter.y;
		foundRight = false;
		for(y=segInfo.cordBound1[z].y-5; y<=segInfo.cordBound2[z].y && !foundRight;y++)
		{
			for(x=segInfo.diskBound2[z].x; x<segInfo.spineBound2[z].x; x++)
			{
				if(binArray2D[y*sizex+x]==1)
				{
					foundRight = true;
					rightSeed = IntVec2(x, y);
				}
				if(maskA[k+y*sizex+x]==maskStruct.mask_air) break;
			}
			if(binArray2D[y*sizey+x]==1)
			{
				foundRight = true;
				rightSeed = IntVec2(x, y);
				break;
			}
		}

		if(foundRight)
		{
			NIH_Algo_Region_Growth_2D(binImg2D, rightSeed, 8, blobImg2D, -1, -1, -1, -1);
			CIS_IPA_Dilate(blobImg2D, 1, true);

			for(y=0, k1=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k1++)
				{
					if(blobArray2D[k1]==1 && x>segInfo.cordBound2[z].x+5)
					{
						maskA[k1+k] = maskStruct.mask_rib;
					}
				}
			}

			// eliminate spongy bone inside and between ribs
			for(y=segInfo.spineBound1[z].y; y<=segInfo.spineBound2[z].y; y++)
			{
				k1 = k+y*sizex+segInfo.spineBound2[z].x;
				flagStop=false;
				for(x=segInfo.spineBound2[z].x; x>=segInfo.cordBound2[z].x+5 && !flagStop; x--, k1--)
				{
					if(maskA[k1]==maskStruct.mask_rib)
					{
						for(;x>=segInfo.cordBound2[z].x+5; x--, k1--)
						{
							if(binArray2D[k1-k]==1 && maskA[k1]!=maskStruct.mask_rib)
							{
								flagStop = true;
								break;
							}

							if(maskA[k1]==maskStruct.mask_spongyBone || maskA[k1]==maskStruct.mask_corticalBone) 
								maskA[k1]=maskStruct.mask_rib;
						}
						flagStop = true;
					}	// if maskA
				}	// for x
			}	// for y
		
		}	// if foundRight

	}

	// 3D region growing to get rid of disconnected bones
	//
	CIS_Array_Image3D_short *binImg;
	CIS_Array_Image3D_short *blobImg;
	short *binArray;
	short *blobArray;

	binImg = new CIS_Array_Image3D_short(sizex, sizey, sizez);
	blobImg = new CIS_Array_Image3D_short(sizex, sizey, sizez);
	binArray = binImg->GetArray();
	blobArray = blobImg->GetArray();

	for(k=0; k<sizexyz; k++)
	{
		if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_spongyBone || maskA[k]==maskStruct.mask_spinalCord) binArray[k]=1;
		else binArray[k]=0;
	}

	// find the largest connected component
	intDynArray blobRank;
	intVec3DynArray blobCentroid;

	NIH_Algo_Blob_Labelling_3D(binImg, blobImg, blobRank, blobCentroid, false);

	largestSize = 0;
	largestBlob = -1;
	for(k=0; k<blobRank.GetSize(); k++) 
	{
		if(blobRank[k]>largestSize) 
		{
			largestSize=blobRank[k];
			largestBlob = k+1;
		}
	}

	for(k=0; k<sizexyz; k++)
	{
		if(binArray[k]==1)
		{
			if(blobArray[k]!=largestBlob)
			{
				maskA[k]=maskStruct.mask_otherBone;
				binArray[k] = 0;
			}
		}
	}

	// close the interior holes if the slice was previously cut-off by the arbitrary threshold 
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(segInfo.spineBound1[z].y==-1) continue;

		if(segInfo.diskBound1[z].y<segPara.predefinedBound1.y)		// only when the spineBound is small than half of the image size
		{
			k = z*sizexy;
			count =0;
			for(k2=0, y=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k++, k2++)
				{
					if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_spongyBone)
					{
						binArray2D[k2] = 1;
						count++;
					}
					else binArray2D[k2]=0;
				}
			}

			if(count==0) continue;

			CIS_IPA_Dilate(binImg2D, segPara.closeCorticalHoleIteration, true);
			CIS_IPA_Erode(binImg2D, segPara.closeCorticalHoleIteration, true);

			k = z*sizexy;
			for(k2=0, y=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k++, k2++)
				{
					if(binArray2D[k2]==1) 
					{
						if(maskA[k]!=maskStruct.mask_corticalBone && maskA[k]!=maskStruct.mask_spongyBone
							&& maskA[k]!=maskStruct.mask_spinalCord)
							maskA[k]=maskStruct.mask_spongyBone;
					}
				}
			}
		}	// if spineBound1[z]
	}	// for z

	// adjust outlines of spinal cord bounding box
	
	// adjust the disk bounding box, and spinous process bounding box
	int gap;
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;

		// vertebra body
		segInfo.diskBound1[z].x=segInfo.cordBound1[z].x;
		segInfo.diskBound2[z].x=segInfo.cordBound2[z].x;
		segInfo.diskBound1[z].y=segInfo.cordBound1[z].y;
		segInfo.diskBound2[z].y=(segInfo.cordBound1[z].y+segInfo.cordBound2[z].y)/2;
		// find the tip of vertebra body
		k = z*sizexy;
		for(x=segInfo.cordBound1[z].x; x<=segInfo.cordBound2[z].x; x++)
		{
			for(y=segInfo.spineBound1[z].y; y<segInfo.cordBound1[z].y; y++)
			{
				if(maskA[k+y*sizex+x]==maskStruct.mask_corticalBone || maskA[k+y*sizex+x]==maskStruct.mask_spongyBone) break;
			}
			if(y<segInfo.diskBound1[z].y) segInfo.diskBound1[z].y=y;
		}
		for(y=segInfo.diskBound1[z].y; y<segInfo.cordBound1[z].y; y++)
		{
			// find the left end
			gap = 0;
			for(x=segInfo.cordBound1[z].x; x>=segInfo.spineBound1[z].x && gap<5; x--)
			{
				if(maskA[k+y*sizex+x]!=maskStruct.mask_corticalBone && maskA[k+y*sizex+x]!=maskStruct.mask_spongyBone 
					&& maskA[k+y*sizex+x]!=maskStruct.mask_spinalCord) 
				{
					gap ++;
				}
				else gap=0;
			}
			if(x+gap<segInfo.diskBound1[z].x) segInfo.diskBound1[z].x=x+gap;

			// find the right end
			gap = 0;
			for(x=segInfo.cordBound2[z].x; x<segInfo.spineBound2[z].x && gap<5; x++)
			{
				if(maskA[k+y*sizex+x]!=maskStruct.mask_corticalBone && maskA[k+y*sizex+x]!=maskStruct.mask_spongyBone
					&& maskA[k+y*sizex+x]!=maskStruct.mask_spinalCord) 
				{
					gap ++;
				}
				else gap=0;
			}
			if(x-gap>segInfo.diskBound2[z].x) segInfo.diskBound2[z].x=x-gap;
		}

		// add margin
		segInfo.diskBound1[z].x -=2;
		segInfo.diskBound2[z].x +=2;
		segInfo.diskBound1[z].y -=2;

		// spinous process
		segInfo.sprocessBound1[z].x=(segInfo.cordBound1[z].x+segInfo.diskBound1[z].x*2)/3;
		segInfo.sprocessBound2[z].x=(segInfo.cordBound2[z].x+segInfo.diskBound2[z].x*2)/3;
		segInfo.sprocessBound1[z].y=(segInfo.cordBound1[z].y+segInfo.cordBound2[z].y)/2;;
		segInfo.sprocessBound2[z].y=segInfo.cordBound2[z].y;
		for(x=segInfo.cordBound1[z].x; x<=segInfo.cordBound2[z].x; x++)
		{
			for(y=segInfo.cordBound1[z].y; y<=segInfo.spineBound2[z].y; y++)
			{
				if(maskA[k+y*sizex+x]==maskStruct.mask_corticalBone || maskA[k+y*sizex+x]==maskStruct.mask_spongyBone)
				{
					if(y>segInfo.sprocessBound2[z].y) segInfo.sprocessBound2[z].y=y;
				}
			}
		}
		segInfo.sprocessBound2[z].y += 2;
	}	// for z
	

	// adjust outliers of the bounding box of disk
	// adjust height first, then adjust width (left and right separately)
	// -- can't be out of 2sd of all bbox
	// -- width should be smaller than height
	// -- the disk can't be three times larger than the cord,
	// -- can't be less than 1.5 times of cord
	doubleDynArray heightArray, leftWidthArray, rightWidthArray;
	double mean, variance, skewness, kurtosis, absoluteDeviation, standardDeviation;

	heightArray.SetSize(0);
	leftWidthArray.SetSize(0);
	rightWidthArray.SetSize(0);
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;
		heightArray.Add(segInfo.diskBound2[z].y-segInfo.diskBound1[z].y);
		rightWidthArray.Add(segInfo.diskBound2[z].x-segInfo.cordBound2[z].x);
		leftWidthArray.Add(segInfo.cordBound1[z].x-segInfo.diskBound1[z].x);
	}
	// adjust height
	CIS_Algo_ComputeStatisticalMoments(heightArray, mean, variance, skewness, kurtosis, absoluteDeviation, standardDeviation);
	for(z=segInfo.bound1.z, k=0; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;
		if(fabs(heightArray[k]-mean)>standardDeviation*2)
		{
			// force it to be mean
			segInfo.diskBound1[z].y = segInfo.diskBound2[z].y-mean;
			heightArray[k] = mean;
		}

		k++;
	}

	// adjust width
	for(z=segInfo.bound1.z, k=0; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;
		if(leftWidthArray[k]>heightArray[k]/2)
		{
			// force it to be mean
			segInfo.diskBound1[z].x = segInfo.cordBound1[z].x-heightArray[k]/2;
			leftWidthArray[k] = heightArray[k]/2;
		}
		k++;
	}
	CIS_Algo_ComputeStatisticalMoments(leftWidthArray, mean, variance, skewness, kurtosis, absoluteDeviation, standardDeviation);
	for(z=segInfo.bound1.z, k=0; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;
		if(fabs(leftWidthArray[k]-mean)>standardDeviation*2)
		{
			// force it to be mean
			segInfo.diskBound1[z].x = segInfo.cordBound1[z].x-mean;
			leftWidthArray[k] = mean;
		}
		k++;
	}


	for(z=segInfo.bound1.z, k=0; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;
		if(rightWidthArray[k]>heightArray[k]/2)
		{
			// force it to be mean
			segInfo.diskBound2[z].x = segInfo.cordBound2[z].x+heightArray[k]/2;
			rightWidthArray[k] = heightArray[k]/2;
		}
		k++;
	}
	CIS_Algo_ComputeStatisticalMoments(rightWidthArray, mean, variance, skewness, kurtosis, absoluteDeviation, standardDeviation);
	for(z=segInfo.bound1.z, k=0; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;
		if(fabs(rightWidthArray[k]-mean)>standardDeviation*2)
		{
			// force it to be mean
			segInfo.diskBound2[z].x = segInfo.cordBound2[z].x+mean;
			rightWidthArray[k] = mean;
		}
		k++;
	}

	// also adjust the bound of sprocess
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;
		segInfo.sprocessBound1[z].x=(segInfo.cordBound1[z].x+segInfo.diskBound1[z].x*2)/3;
		segInfo.sprocessBound2[z].x=(segInfo.cordBound2[z].x+segInfo.diskBound2[z].x*2)/3;
	}

	// fit a ellipse inside the bounding box to fill a hole on the boundary 
	// make the ellipse a little smaller than the bounding box
	// current method:
	//		shooting ball to fill the eclipes
	Vec2 *rDir, curP, bp, *rDirInc;
	double rayStep;
	float ballDiameter;
	int j, incStep;
	float px=img3D->Get_Pixel_SizeX();
	int shootingRayNum=24;

	rayStep = 2*3.14159/(double)shootingRayNum;
//	ballDiameter = 2*ballRadius/px;

	incStep = 2;
	rDir = new Vec2[shootingRayNum];
	rDirInc = new Vec2[shootingRayNum];

	for(i=0; i<shootingRayNum; i++)
	{
		rDir[i].x = cos(i*rayStep);
		rDir[i].y = sin(i*rayStep);

		rDirInc[i].x = rDir[i].x*incStep;
		rDirInc[i].y = rDir[i].y*incStep;
	}

	bool flagHit, flagBorder;
	int curk;
	int bbx0, bbx1, bby0, bby1;
	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;

		bbx0 = segInfo.diskBound1[z].x;
		bbx1 = segInfo.diskBound2[z].x;
		bby0 = segInfo.diskBound1[z].y;
		bby1 = segInfo.diskBound2[z].y;

		// first get the border of the disk
		(*binImg2D) = (short)0;
		for(y=bby0; y<=bby1; y++)
		{
			k = z*sizexy+y*sizex+bbx0;
			k1 = y*sizex+bbx0;
			for(x=bbx0; x<=bbx1; x++, k++, k1++)
			{
				if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_spongyBone)
				{
					if(maskA[k+1]==maskStruct.mask_body || maskA[k-1]==maskStruct.mask_body 
						|| maskA[k+sizex]==maskStruct.mask_body || maskA[k-sizex]==maskStruct.mask_body)
					{
						binArray2D[k1]=2;
					};
				}
			}
		}	// for y


		// shooting ray to fill gaps along the border of disk
		ballDiameter =  segInfo.diskBound2[z].y-segInfo.diskBound1[z].y;
		k = z*sizexy;
		for(y=bby0; y<=bby1; y++)
		{
			k1 = y*sizex+bbx0;
			for(x=bbx0; x<=bbx1; x++, k1++)
			{
				if(binArray2D[k1]==2)
				{
					bp = Vec2(x,y);
					for(i=0; i<shootingRayNum; i++)
					{
						flagHit = false;
						flagBorder = false;
						curP = bp;
						for(j=1; j<ballDiameter; j+=incStep)
						{
							curP += rDirInc[i];
							curk = k+(int)curP.x+(int)curP.y*sizex;

							if(curP.x<=bbx0 || curP.x>=bbx1 || curP.y<=bby0 || curP.y>=bby1)
							{
								flagBorder = true;
								break;
							}

							if(binArray2D[curk-k]!=2 && (maskA[curk]==maskStruct.mask_corticalBone || maskA[curk]==maskStruct.mask_spongyBone ))
							{
								flagHit = true;
								break;
							}
						}	// for j

						if(flagHit && !flagBorder && j>1)
						{
							curP = bp;
							for(n=1; n<j; n++)
							{
								curP += rDir[i];
								curk = (int)curP.x+(int)curP.y*sizex;
								if(binArray2D[curk]==0) binArray2D[curk] = 1;
							}
						}	// if flagHit
						
					}	// for i

				}	// if maskArray
			}	// for x
		}	// for y

		// close gap to fill holes missed by shooting balls (2d based)
		// open operation to remove outliers (2d based)
		for(y=bby0; y<=bby1; y++)
		{
			k1 = y*sizex+bbx0;
			for(x=bbx0; x<=bbx1; x++, k1++)
			{
				if(binArray2D[k1]==2) binArray2D[k1]=0;
			}	// for x
		}	// y
		CIS_IPA_Dilate(binImg2D, 2, true);
		CIS_IPA_Erode(binImg2D, 3, true);
		CIS_IPA_Dilate(binImg2D, 1, true);

		k = z*sizexy;
		for(y=bby0; y<=bby1; y++)
		{
			k1 = y*sizex+bbx0;
			for(x=bbx0; x<=bbx1; x++, k1++)
			{
				if(binArray2D[k1]==1 && maskA[k+k1]==maskStruct.mask_body) maskA[k+k1]=maskStruct.mask_spongyBone;
			}	// for x
		}	// y		
	}	// for z


	// re-adjust spinal bbox
	for(z=segInfo.bound1.z, k=0; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;
		c_center.x = cord[sliceCord[z]].centerx;
		c_center.y = cord[sliceCord[z]].centery;

		if(fabs((double)(segInfo.spineBound1[z].x-c_center.x))>fabs((double)(segInfo.diskBound1[z].x-c_center.x))*1.5)
			segInfo.spineBound1[z].x = c_center.x-fabs((double)(segInfo.diskBound1[z].x-c_center.x))*1.5;
		if(fabs((double)(segInfo.spineBound2[z].x-c_center.x))>fabs((double)(segInfo.diskBound2[z].x-c_center.x))*1.5)
			segInfo.spineBound2[z].x = c_center.x+fabs((double)(segInfo.diskBound2[z].x-c_center.x))*1.5;

		if(fabs((double)(segInfo.spineBound1[z].y-c_center.y))>fabs((double)(segInfo.diskBound1[z].y-c_center.y))*1.5)
			segInfo.spineBound1[z].y = c_center.y-fabs((double)(segInfo.diskBound1[z].y-c_center.y))*1.5;
		if(fabs((double)(segInfo.spineBound2[z].y-c_center.y))>fabs((double)(segInfo.sprocessBound2[z].y-c_center.y))*1.5)
			segInfo.spineBound2[z].y = c_center.y+fabs((double)(segInfo.sprocessBound2[z].y-c_center.y))*1.5;
	}


	// cut all the bones outside the bounding box 
	//
	for(z=segInfo.bound1.z, k=0; z<=segInfo.bound2.z; z++)
	{
		if(sliceCord[z]==-1) continue;

		k=z*sizexy;

		for(y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++)
			{
				// convert bone outside bounding box to other bone
				if(maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_spongyBone)
				{
					if(y<segInfo.spineBound1[z].y || y>segInfo.spineBound2[z].y || x<segInfo.spineBound1[z].x || x>segInfo.spineBound2[z].x)
						maskA[k] = maskStruct.mask_otherBone;
				}
			}
		}	// for y
	}	// for z

	
	
	// recompute the bone intensity
	int count_cortical, count_spongy, count_cord;

	segInfo.avg_spinal_cord_intensity = segInfo.avg_spongy_intensity = segInfo.avg_cortical_intensity = 0;
	count_cortical = count_spongy = count_cord = 0;

	for(k=0; k<sizexyz; k++)
	{
		if(maskA[k]==maskStruct.mask_spinalCord)
		{
			count_cord++;
			segInfo.avg_spinal_cord_intensity += imgA[k];
		}
		else if(maskA[k]==maskStruct.mask_spongyBone)
		{
			count_spongy++;
			segInfo.avg_spongy_intensity += imgA[k];
		}
		else if(maskA[k]==maskStruct.mask_corticalBone && imgA[k]>segPara.boneThresh)
		{
			count_cortical++;
			segInfo.avg_cortical_intensity += imgA[k];
		}
	}

	if(count_cord!=0) segInfo.avg_spinal_cord_intensity /= (float)count_cord;
	if(count_spongy!=0) segInfo.avg_spongy_intensity /= (float)count_spongy;
	if(count_cortical!=0) segInfo.avg_cortical_intensity /= (float)count_cortical;

	delete binImg;
	delete blobImg;
	delete binImg2D;
	delete blobImg2D;
	free(cord);

	delete rDir;
	delete rDirInc;

	return CIS_OK;
}



// analyze the bone metastsis in 2d, extract 2D features
int NIH_BoneSegmentation_AnalyzeMetasis2D(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   FeatureStructure* &detections2D, int &numDetections2D,
									   LesionStructure *lesions, int numLesions,
									   PaintingStruct *lesionVoxels,
									   bool flagApplyClassifier, bool flagApply3DClassifier, SvmCommittee *svmCommittee, float svmCutoff,
									   int sliceRange1, int sliceRange2)
{
	int x, y, z, k;
	int sizex, sizey, sizez, sizexy, sizexyz;
	short *imgA, *maskA;
	QStringList featureNameList, featureValueList;
	double svmVote, svmScore, sumRes;

	sizex = img3D->Num_Cols();
	sizey = img3D->Num_Rows();
	sizez = img3D->Num_Levels();

	imgA = img3D->GetArray();
	maskA = maskImg3D->GetArray();

	sizexyz = sizex*sizey*sizez;
	sizexy = sizex*sizey;

	CIS_Array_Image2D_short *binImg2D, *blobImg2D;
	short *binArray2D, *blobArray2D;
	intDynArray blobRank2D;
	intVec2DynArray blobCentroid2D;
	intDynArray blobStatus;

	int count, k2, b;

	binImg2D = new CIS_Array_Image2D_short(sizex, sizey);
	blobImg2D = new CIS_Array_Image2D_short(sizex, sizey);
	binArray2D = binImg2D->GetArray();
	blobArray2D = blobImg2D->GetArray();
	
	int maxDetectionNumber;
	intVec2DynArray detectionRegion;
	int matchedLesion, matchedLesionType;
	float matchedOverlap, matchedLesionSize;

	FeatureStructure curFeature;
	numDetections2D = 0;
	maxDetectionNumber = sizez*10;	// assume 10 detection per slice

	if(detections2D!=NULL) free(detections2D);
	detections2D = (FeatureStructure *)malloc(maxDetectionNumber*sizeof(FeatureStructure));

	// assemble the feature Name list
	if(flagApplyClassifier)
	{
		featureNameList.push_back("relCoordx");
		featureNameList.push_back("relCoordy");
		featureNameList.push_back("relCoordz");
		featureNameList.push_back("distToBoundary");
		featureNameList.push_back("area");
		featureNameList.push_back("volume");
		featureNameList.push_back("perimeter");
		featureNameList.push_back("primaryAxisLength");
		featureNameList.push_back("secondaryAxisLength");
		featureNameList.push_back("outerBorderRatio");
		featureNameList.push_back("aspectRatio");
		featureNameList.push_back("spherecity");
		featureNameList.push_back("compactness");
		featureNameList.push_back("roundness");
		featureNameList.push_back("shapeComplexity_f1");
		featureNameList.push_back("shapeComplexity_f2");
		featureNameList.push_back("shapeComplexity_f21");
		featureNameList.push_back("meanIntensity");
		featureNameList.push_back("stdevIntensity");
		featureNameList.push_back("skewnessIntensity");
		featureNameList.push_back("kurtosisIntensity");
		featureNameList.push_back("interiorIntensity");
		featureNameList.push_back("borderIntensity");
		featureNameList.push_back("outsideIntensity");
		featureNameList.push_back("outsideIntensityDev");
		featureNameList.push_back("innerOuterContrast");
		featureNameList.push_back("borderThickness");
	}


	int errorindex [22];
			for(int i=0; i<22; i++) 
				errorindex[i]=0;

	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		//if spinal cord not present in current slice
		if(segInfo.diskBound1[z]==-1) continue;

		//make binArray2D for boneMetasis in spinal cord and count detection points
		for(k2=0; k2<sizexy; k2++) 
			binArray2D[k2]=0;

		k = z*sizexy;
		count =0;
		for(y=segInfo.diskBound1[z].y; y<=segInfo.diskBound2[z].y; y++)
		{
			k2 = y*sizex+segInfo.diskBound1[z].x;
			for(x=segInfo.diskBound1[z].x; x<=segInfo.diskBound2[z].x; x++, k2++)
			{
				if(maskA[k+k2]==maskStruct.mask_boneMetasis)
				{
					binArray2D[k2] = 1;
					count++;
				}
			}
		}
		//in spinal process
		for(y=segInfo.sprocessBound1[z].y; y<=segInfo.sprocessBound2[z].y; y++)
		{
			k2 = y*sizex+segInfo.sprocessBound1[z].x;
			for(x=segInfo.sprocessBound1[z].x; x<=segInfo.sprocessBound2[z].x; x++, k2++)
			{
				if(maskA[k+k2]==maskStruct.mask_boneMetasis)
				{
					binArray2D[k2] = 1;
					count++;
				}
			}
		}

		if(count==0) continue;
	
		// smooth the detection region
		CIS_IPA_Erode(binImg2D, 1, false);
		CIS_IPA_Dilate(binImg2D, 1, false);
///		CIS_IPA_Dilate(binImg2D, 1, false);
		//CIS_IPA_Erode(binImg2D, 1, true);

		
		
		NIH_Algo_Blob_Labelling(binImg2D, blobImg2D, blobRank2D, blobCentroid2D, false);

		// analyze the shape blob by blob

		blobStatus.SetSize(blobRank2D.GetSize());
		//for(int i=0; i<blobRank2D.GetSize(); i++) std::cout<<"\n Blob Rank: "<<blobRank2D[i]<<"\n";

		for(b=0; b<blobStatus.GetSize(); b++) 
		{
			if(blobRank2D[b]<segPara.metasisSizeThresh) blobStatus[b]=-1;	// too small
//			else if(blobCentroid2D[b].y>c_center.y) blobStatus[b]=-2;	// lower than spinal cord
			else blobStatus[b] = 1;
		}	// for b
		
		// compute detection features, and export to a csv file
		//
		for(b=0; b<blobStatus.GetSize(); b++)
		{
			if(sliceRange1>0 && z<sliceRange1) blobStatus[b]=-1;
			if(sliceRange2>0 && z>=sliceRange2) blobStatus[b]=-1;

			if(blobStatus[b]==-1) continue;

			// fill the detection region for a single blob
			detectionRegion.SetSize(0);
			for(k2=0, y=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k++, k2++)
				{
					if(blobArray2D[k2]==b+1) 
					{
						detectionRegion.Add(IntVec2(x,y));
					}
				}
			}

			// compute features
			NIH_BoneSegmentation_ComputeMetasisFeatures2D(curFeature, z, detectionRegion, img3D, maskImg3D, maskStruct, segPara, segInfo);

			// filtering to remove obvious false detections

			
			
			// the metasis can't appear at bottom of the vetebra body and on the side of spinal canal
			// usually are artifact from space between vertebra
			if(curFeature.relCoordy<0 && curFeature.outerBorderRatio>0.2 &&
				((curFeature.relCoordy>-0.2 && fabs(curFeature.relCoordx)>0.1) )) 
			{
				// in disk
				blobStatus[b]=-1;
				std::cout<<"Slice "<<z<<" Area "<<curFeature.area<<" IntI "<<curFeature.interiorIntensity<<" OutI "<<curFeature.outsideIntensity<<" NeighI "<<curFeature.neighborIntensity<<" Error 1 \n";
				errorindex[0]++;
				continue;
			}
			else if(curFeature.relCoordy>0  && fabs(curFeature.relCoordx)<0.5 && curFeature.area < 0.5)
			{
				// in spinal process (relative center)
				blobStatus[b]=-1;
				std::cout<<"Slice "<<z<<" Area "<<curFeature.area<<" IntI "<<curFeature.interiorIntensity<<" OutI "<<curFeature.outsideIntensity<<" NeighI "<<curFeature.neighborIntensity<<" Error 2 \n";
				errorindex[1]++;
				continue;
			}
//			else if(curFeature.relCoordy<0 && curFeature.relCoordy>-0.5 && curFeature.outerBorderRatio>0.4)
//			{
//				// in disk
//				blobStatus[b]=-1;
//				continue;
//			}
			else if(curFeature.relCoordy<0 && fabs(curFeature.relCoordx)>0.3 && curFeature.outerBorderRatio>0.4)
			{
				// in disk, close to boundary and outer border too big
				blobStatus[b]=-1;
				std::cout<<"Slice "<<z<<" Area "<<curFeature.area<<" IntI "<<curFeature.interiorIntensity<<" OutI "<<curFeature.outsideIntensity<<" NeighI "<<curFeature.neighborIntensity<<" Error 3 \n";
				errorindex[2]++;
				continue;
			}
			else if(curFeature.cordBorderRatio>0.33)
			{
				// next to cord
				blobStatus[b]=-1;
				std::cout<<"Slice "<<z<<" Area "<<curFeature.area<<" IntI "<<curFeature.interiorIntensity<<" OutI "<<curFeature.outsideIntensity<<" NeighI "<<curFeature.neighborIntensity<<" Error 4 \n";
				errorindex[3]++;
				continue;
			}
			else if(curFeature.relCoordx<0.1 && curFeature.area<0.25 && curFeature.cordBorderRatio>0.1)
			{
				// next to cord
				blobStatus[b]=-1;
				std::cout<<"Slice "<<z<<" Area "<<curFeature.area<<" IntI "<<curFeature.interiorIntensity<<" OutI "<<curFeature.outsideIntensity<<" NeighI "<<curFeature.neighborIntensity<<" Error 5 \n";
				errorindex[4]++;
				continue;
			}
			else if(curFeature.relCoordy>0 && curFeature.outerBorderRatio>0.2)
			{
				// in process and outer border too big
				blobStatus[b]=-1;
				std::cout<<"Slice "<<z<<" Area "<<curFeature.area<<" IntI "<<curFeature.interiorIntensity<<" OutI "<<curFeature.outsideIntensity<<" NeighI "<<curFeature.neighborIntensity<<" Error 6 \n";
				errorindex[5]++;
				continue;
			}
			else if(curFeature.outerBorderRatio>0.3)
			{
				// outer border too big
				blobStatus[b]=-1;
				std::cout<<"Slice "<<z<<" Area "<<curFeature.area<<" IntI "<<curFeature.interiorIntensity<<" OutI "<<curFeature.outsideIntensity<<" NeighI "<<curFeature.neighborIntensity<<" Error 7 \n";
				errorindex[6]++;
				continue;
			}
			else if(curFeature.aspectRatio>2 && curFeature.outerBorderRatio>0.3)
			{
				// too long and on the border
				blobStatus[b]=-1;
				std::cout<<"Slice "<<z<<" Area "<<curFeature.area<<" IntI "<<curFeature.interiorIntensity<<" OutI "<<curFeature.outsideIntensity<<" NeighI "<<curFeature.neighborIntensity<<" Error 8 \n";
				errorindex[7]++;
				continue;
			}
			else if(curFeature.aspectRatio>3 && curFeature.outerBorderRatio+curFeature.corticalBorderRatio>0.95 && curFeature.outerBorderRatio>0.1)
			{
				// too long and on the gap between rib and vertebra
				blobStatus[b]=-1;
				std::cout<<"Slice "<<z<<" Area "<<curFeature.area<<" IntI "<<curFeature.interiorIntensity<<" OutI "<<curFeature.outsideIntensity<<" NeighI "<<curFeature.neighborIntensity<<" Error 9 \n";
				errorindex[8]++;
				continue;
			}
			
		/*	if(curFeature.area <50 && curFeature.corticalBorderRatio>0.5)
			{
				blobStatus[b]=-1;
				std::cout<<"Slice "<<z<<" Area "<<curFeature.area<<" CBR "<<curFeature.corticalBorderRatio<<" Error 10 \n";
				errorindex[9]++;
				continue;
			}*/

			// next to bbounding box and small
			if(curFeature.area <0.1 &&  curFeature.bboxBorderRatio>0.05 && curFeature.relCoordy<0)
			{
				blobStatus[b]=-1;
				std::cout<<"Slice "<<z<<" Area "<<curFeature.area<<" IntI "<<curFeature.interiorIntensity<<" OutI "<<curFeature.outsideIntensity<<" NeighI "<<curFeature.neighborIntensity<<" Error 11 \n";
				errorindex[10]++;
				continue;
			}

			if(curFeature.bboxBorderRatio>0.2 && curFeature.aspectRatio>3)
			{
				blobStatus[b]=-1;
				std::cout<<"Slice "<<z<<" Area "<<curFeature.area<<" IntI "<<curFeature.interiorIntensity<<" OutI "<<curFeature.outsideIntensity<<" NeighI "<<curFeature.neighborIntensity<<" Error 12 \n";
				errorindex[11]++;
				continue;
			}

			// next to bbounding box and small
			if(curFeature.area <0.25 &&  curFeature.bboxBorderRatio>0.3 && curFeature.relCoordy>0)
			{
				blobStatus[b]=-1;
				std::cout<<"Slice "<<z<<" Area "<<curFeature.area<<" IntI "<<curFeature.interiorIntensity<<" OutI "<<curFeature.outsideIntensity<<" NeighI "<<curFeature.neighborIntensity<<" Error 13 \n";
				errorindex[12]++;
				continue;
			}

			// next to air
			if(curFeature.airBorderRatio>0.01)
			{
				blobStatus[b]=-1;
				std::cout<<"Slice "<<z<<" Area "<<curFeature.area<<" IntI "<<curFeature.interiorIntensity<<" OutI "<<curFeature.outsideIntensity<<" NeighI "<<curFeature.neighborIntensity<<" Error 14 \n";
				errorindex[13]++;
				continue;
			}

			// parallel to cord
/*			if(curFeature.centery>segInfo.cordBound1[z].y+2 && curFeature.centery<segInfo.cordBound2[z].y-2)
			{
				blobStatus[b]=-1;
				continue;
			}
*/
			// next to cord
			if(curFeature.cordBorderRatio>0.1 && curFeature.outerBorderRatio>0.2)
			{
				blobStatus[b]=-1;
				std::cout<<"Slice "<<z<<" Area "<<curFeature.area<<" IntI "<<curFeature.interiorIntensity<<" OutI "<<curFeature.outsideIntensity<<" NeighI "<<curFeature.neighborIntensity<<" Error 15 \n";
				errorindex[14]++;
				continue;
			}

			// next to rib
			if(curFeature.ribBorderRatio>0.01)
			{
				blobStatus[b]=-1;
				std::cout<<"Slice "<<z<<" Area "<<curFeature.area<<" IntI "<<curFeature.interiorIntensity<<" OutI "<<curFeature.outsideIntensity<<" NeighI "<<curFeature.neighborIntensity<<" Error 16 \n";
				errorindex[15]++;
				continue;
			}

			// long object
			if(curFeature.aspectRatio>3.5)
			{
				blobStatus[b]=-1;
				std::cout<<"Slice "<<z<<" Area "<<curFeature.area<<" IntI "<<curFeature.interiorIntensity<<" OutI "<<curFeature.outsideIntensity<<" NeighI "<<curFeature.neighborIntensity<<" Error 17 \n";
				errorindex[16]++;
				continue;
			}

			// intensity too high or too low
			if(curFeature.interiorIntensity>2300 || curFeature.interiorIntensity<1000)
			{
				blobStatus[b]=-1;
				std::cout<<"Slice "<<z<<" Area "<<curFeature.area<<" IntI "<<curFeature.interiorIntensity<<" OutI "<<curFeature.outsideIntensity<<" NeighI "<<curFeature.neighborIntensity<<" Error 18 \n";
				errorindex[17]++;
				continue;
			}

			// not enough contrast
			if(fabs(curFeature.neighborIntensity-curFeature.interiorIntensity)<150 && fabs(curFeature.outsideIntensity-curFeature.interiorIntensity)<100)
			{
				blobStatus[b]=-1;
				std::cout<<"Slice "<<z<<" Area "<<curFeature.area<<" IntI "<<curFeature.interiorIntensity<<" OutI "<<curFeature.outsideIntensity<<" NeighI "<<curFeature.neighborIntensity<<" Error 19 \n";
				errorindex[18]++;
				continue;
			}

			// not enough contrast
			if(fabs(curFeature.outsideIntensity-curFeature.interiorIntensity)<100)
			{
				blobStatus[b]=-1;
				std::cout<<"Slice "<<z<<" Area "<<curFeature.area<<" IntI "<<curFeature.interiorIntensity<<" OutI "<<curFeature.outsideIntensity<<" NeighI "<<curFeature.neighborIntensity<<" Error 20 \n";
				errorindex[19]++;
				continue;
			}

			// not enough contrast for small detection
/*			if(curFeature.area<0.2 && (curFeature.outsideIntensity-curFeature.meanIntensity<50 || curFeature.neighborIntensity-curFeature.meanIntensity<100))
			{
				blobStatus[b]=-1;
				continue;
			}
*/
			// not enough contrast for large detections
///			if(curFeature.area>curFeature.neighborArea/5 && curFeature.neighborIntensity-curFeature.interiorIntensity<150 
///				&& curFeature.outsideIntensity-curFeature.interiorIntensity<150)
///			{
///				blobStatus[b]=-1;
///				continue;
///			}
			if(curFeature.area>curFeature.neighborArea/5 && fabs(curFeature.outsideIntensity-curFeature.interiorIntensity)<150)
			{
				blobStatus[b]=-1;
				std::cout<<"Slice "<<z<<" Area "<<curFeature.area<<" IntI "<<curFeature.interiorIntensity<<" OutI "<<curFeature.outsideIntensity<<" NeighI "<<curFeature.neighborIntensity<<" Error 21 \n";
				errorindex[20]++;
				continue;
			}
			if(curFeature.area>0.5 && fabs(curFeature.outsideIntensity-curFeature.interiorIntensity)<100)
			{
				blobStatus[b]=-1;
				std::cout<<"Slice "<<z<<" Area "<<curFeature.area<<" IntI "<<curFeature.interiorIntensity<<" OutI "<<curFeature.outsideIntensity<<" NeighI "<<curFeature.neighborIntensity<<" Error 22 \n";
				errorindex[21]++;
				continue;
			}

			// close to boundary and small
/*			if(curFeature.area<0.2 && curFeature.distToBoundary<0.3)
			{
				blobStatus[b]=-1;
				continue;
			}
*/

			
			// pass the detections to a classifier
			if(flagApplyClassifier)
			{
				featureValueList.clear();
				featureValueList.push_back(QString::number(curFeature.relCoordx));
				featureValueList.push_back(QString::number(curFeature.relCoordy));
				featureValueList.push_back(QString::number(curFeature.relCoordz));
				featureValueList.push_back(QString::number(curFeature.distToBoundary));
				featureValueList.push_back(QString::number(curFeature.area));
				featureValueList.push_back(QString::number(curFeature.volume));
				featureValueList.push_back(QString::number(curFeature.perimeter));
				featureValueList.push_back(QString::number(curFeature.primaryAxisLength));
				featureValueList.push_back(QString::number(curFeature.secondaryAxisLength));
				featureValueList.push_back(QString::number(curFeature.outerBorderRatio));
				featureValueList.push_back(QString::number(curFeature.aspectRatio));
				featureValueList.push_back(QString::number(curFeature.spherecity));
				featureValueList.push_back(QString::number(curFeature.compactness));
				featureValueList.push_back(QString::number(curFeature.roundness));
				featureValueList.push_back(QString::number(curFeature.shapeComplexity_f1));
				featureValueList.push_back(QString::number(curFeature.shapeComplexity_f2));
				featureValueList.push_back(QString::number(curFeature.shapeComplexity_f21));
				featureValueList.push_back(QString::number(curFeature.meanIntensity));
				featureValueList.push_back(QString::number(curFeature.stdevIntensity));
				featureValueList.push_back(QString::number(curFeature.skewnessIntensity));
				featureValueList.push_back(QString::number(curFeature.kurtosisIntensity));
				featureValueList.push_back(QString::number(curFeature.interiorIntensity));
				featureValueList.push_back(QString::number(curFeature.borderIntensity));
				featureValueList.push_back(QString::number(curFeature.outsideIntensity));
				featureValueList.push_back(QString::number(curFeature.outsideIntensityDev));
				featureValueList.push_back(QString::number(curFeature.innerOuterContrast));
				featureValueList.push_back(QString::number(curFeature.borderThickness));

//				svmCommittee->Predict(featureNameList, featureValueList, svmVote, svmScore);
				svmCommittee->PredictWithCombinationSumRes(featureNameList, featureValueList, svmVote, svmScore, sumRes);

				curFeature.svmVote = svmVote;
				curFeature.svmScore = sumRes;

				if(curFeature.svmScore<svmCutoff) blobStatus[b]=-2;
			}	// if flagApplyClassifier
			
			if(blobStatus[b]==-2) continue;
			
			// test if the detection overlaped with real lesions
			NIH_BoneSegmentation_MatchDetections2D(detectionRegion, z, lesions, numLesions, lesionVoxels, maskImg3D, segPara, 
				matchedLesion, matchedLesionType, matchedLesionSize, matchedOverlap);
			curFeature.matchedLesion = matchedLesion;
			curFeature.matchedLesionType = matchedLesionType;
			curFeature.matchedLesionSize = matchedLesionSize;
			curFeature.matchedOverlap = matchedOverlap;

			if(numDetections2D>=maxDetectionNumber)
			{
				maxDetectionNumber *=2;
				detections2D = (FeatureStructure *) realloc(detections2D, maxDetectionNumber*sizeof(FeatureStructure));
			}
			detections2D[numDetections2D] = curFeature;
			numDetections2D ++;
		}	// for b

		// fill the mask
		k = z*sizexy;
		for(k2=0, y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++, k2++)
			{
				if(binArray2D[k2]==1) 
				{
					if(blobStatus[blobArray2D[k2]-1]==1)
					{
						maskA[k]=maskStruct.mask_boneMetasis;
					}
					else if(blobStatus[blobArray2D[k2]-1]==-2  && maskA[k]==maskStruct.mask_boneMetasis)
					{
						maskA[k]=maskStruct.mask_falseDetection;
					}
					else if(blobStatus[blobArray2D[k2]-1]==-1 && maskA[k]==maskStruct.mask_boneMetasis)
					{
						maskA[k]=maskStruct.mask_falseDetection;
///						maskA[k]=maskStruct.mask_spongyBone;
					}
				}	// if binArray
				else if(maskA[k]==maskStruct.mask_boneMetasis)
					maskA[k]=maskStruct.mask_spongyBone;
			}
		}	// for k2
	}	// for z

	for(int i=0; i<22; i++)
			{
				std::cout<<"Error "<<i+1<<" frequency: "<<errorindex[i]<<"\n";
			}


	delete binImg2D;
	delete blobImg2D;
	return CIS_OK;
}


// analyze the bone metastasis in 3d, extract 3D features
int NIH_BoneSegmentation_AnalyzeMetasis3D(CIS_Array_Image3D_short *img3D,
										  CIS_Array_Image3D_short *maskImg3D,
										  SpineSegmentationMaskStructure &maskStruct,
										  SpineSegmentationParameter &segPara,
										  SpineSegmentationInfo &segInfo,
										  FeatureStructure* &detections3D, int &numDetections3D,
										  DetectionStructure* &detections, int &numDetections,
										  PaintingStruct *detectionVoxels,
										  LesionStructure *lesions, int numLesions,
										  PaintingStruct *lesionVoxels,
										  bool flagApply3DClassifier, SvmCommittee *svmCommittee, float svmCutoff,
										  int sliceRange1, int sliceRange2)
{
	int x, y, z, k;
	int sizex, sizey, sizez, sizexy, sizexyz;
	short *imgA, *maskA;
	QStringList featureNameList, featureValueList;
	double svmVote, svmScore, sumRes;

	sizex = img3D->Num_Cols();
	sizey = img3D->Num_Rows();
	sizez = img3D->Num_Levels();

	imgA = img3D->GetArray();
	maskA = maskImg3D->GetArray();

	sizexyz = sizex*sizey*sizez;
	sizexy = sizex*sizey;

	CIS_Array_Image3D_short *binImg3D;
	CIS_Array_Image3D_short *blobImg3D;
	short *binArray3D;
	short *blobArray3D;
	intDynArray blobRank3D;
	intVec3DynArray blobCentroid3D;
	intDynArray blobStatus;

	int count, k2, b;

	binImg3D = new CIS_Array_Image3D_short(sizex, sizey, sizez);
	blobImg3D = new CIS_Array_Image3D_short(sizex, sizey, sizez);
	binArray3D = binImg3D->GetArray();
	blobArray3D = blobImg3D->GetArray();

	int maxDetectionNumber;
	intVec3DynArray detectionRegion;
	int matchedLesion, matchedLesionType;
	float matchedOverlap, matchedLesionSize;

	FeatureStructure curFeature;
	DetectionStructure curDetection;
	numDetections3D = 0;
	numDetections = 0;
	maxDetectionNumber = sizez*10;	// assume 10 detection per slice

	if(detections3D!=NULL)
	detections3D = (FeatureStructure *)malloc(maxDetectionNumber*sizeof(FeatureStructure));
	detections = (DetectionStructure *)malloc(maxDetectionNumber*sizeof(DetectionStructure));

	// assemble the feature name list
	if(flagApply3DClassifier)
	{
		featureNameList.push_back("relCoordx");
		featureNameList.push_back("relCoordy");
		featureNameList.push_back("relCoordz");
		featureNameList.push_back("distToBoundary");
		featureNameList.push_back("surfaceArea");
		featureNameList.push_back("volume");
		featureNameList.push_back("primaryAxisLength");
		featureNameList.push_back("secondaryAxisLength");
		featureNameList.push_back("outerBorderRatio");
		featureNameList.push_back("aspectRatio");
		featureNameList.push_back("spherecity");
		featureNameList.push_back("shapeComplexity_f1");
		featureNameList.push_back("shapeComplexity_f2");
		featureNameList.push_back("shapeComplexity_f21");
		featureNameList.push_back("meanIntensity");
		featureNameList.push_back("stdevIntensity");
		featureNameList.push_back("skewnessIntensity");
		featureNameList.push_back("kurtosisIntensity");
		featureNameList.push_back("interiorIntensity");
		featureNameList.push_back("borderIntensity");
		featureNameList.push_back("outsideIntensity");
		featureNameList.push_back("outsideIntensityDev");
		featureNameList.push_back("innerOuterContrast");
	}

	//initialization
	count = 0;
	for(k=0; k<sizexyz; k++)
	{
		if(maskA[k]==maskStruct.mask_boneMetasis) 
		{
			count++;
			binArray3D[k]=1;
		}
		else binArray3D[k]=0;
	}

	// smooth the detection region
	CIS_IPA_3D_Erode(binImg3D, 1, false);
	CIS_IPA_3D_Dilate(binImg3D, 1, false);

	NIH_Algo_Blob_Labelling_3D(binImg3D, blobImg3D, blobRank3D, blobCentroid3D, false);

	// analyze the shape blob by blob

	blobStatus.SetSize(blobRank3D.GetSize());

	for(b=0; b<blobRank3D.GetSize(); b++)
	{
		if(blobRank3D[b]<segPara.metasisSizeThresh) blobStatus[b]=-1;
		else blobStatus[b]=1;
	}

	// compute detection features, and export to a csv file
	//
	int minz, maxz;
	for(b=0; b<blobRank3D.GetSize(); b++)
	{
		if(blobStatus[b]==-1) continue;
		maxz = -1;
		minz = sizez;
		// fill the detection region for a single blob
		detectionRegion.SetSize(0);
		for(k=0, z=0; z<sizez; z++)
		{
			for(y=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k++)
				{
					if(blobArray3D[k]==b+1)
					{
						detectionRegion.Add(IntVec3(x,y,z));
						if(z>maxz)
						{
							maxz = z;
						}
						if(z<minz)
						{
							minz = z;
						}
					}
				}
			}
		}
		curDetection.minz = minz;
		curDetection.maxz = maxz;

		// compute features
		NIH_BoneSegmentation_ComputeMetasisFeatures3D(curFeature, detectionRegion, img3D, maskImg3D, maskStruct, segPara, segInfo);
		curDetection.detectionLoc = IntVec3(curFeature.centerx, curFeature.centery, curFeature.centerz);
		curDetection.detectionVolume = curFeature.volume;
		curDetection.detectionType = -1;
		curDetection.matchedLesion = -1;
		curDetection.detectionSize = curFeature.primaryAxisLength;

		// filtering to remove obvious false detections (omitted for now)

		if(flagApply3DClassifier)
		{
				featureValueList.clear();
				featureValueList.push_back(QString::number(curFeature.relCoordx));
				featureValueList.push_back(QString::number(curFeature.relCoordy));
				featureValueList.push_back(QString::number(curFeature.relCoordz));
				featureValueList.push_back(QString::number(curFeature.distToBoundary));
				featureValueList.push_back(QString::number(curFeature.surfaceArea));
				featureValueList.push_back(QString::number(curFeature.volume));
				featureValueList.push_back(QString::number(curFeature.perimeter));
				featureValueList.push_back(QString::number(curFeature.primaryAxisLength));
				featureValueList.push_back(QString::number(curFeature.secondaryAxisLength));
				featureValueList.push_back(QString::number(curFeature.outerBorderRatio));
				featureValueList.push_back(QString::number(curFeature.aspectRatio));
				featureValueList.push_back(QString::number(curFeature.spherecity));
				featureValueList.push_back(QString::number(curFeature.shapeComplexity_f1));
				featureValueList.push_back(QString::number(curFeature.shapeComplexity_f2));
				featureValueList.push_back(QString::number(curFeature.shapeComplexity_f21));
				featureValueList.push_back(QString::number(curFeature.meanIntensity));
				featureValueList.push_back(QString::number(curFeature.stdevIntensity));
				featureValueList.push_back(QString::number(curFeature.skewnessIntensity));
				featureValueList.push_back(QString::number(curFeature.kurtosisIntensity));
				featureValueList.push_back(QString::number(curFeature.interiorIntensity));
				featureValueList.push_back(QString::number(curFeature.borderIntensity));
				featureValueList.push_back(QString::number(curFeature.outsideIntensity));
				featureValueList.push_back(QString::number(curFeature.outsideIntensityDev));
				featureValueList.push_back(QString::number(curFeature.innerOuterContrast));

				svmCommittee->PredictWithCombinationSumRes(featureNameList, featureValueList, svmVote, svmScore, sumRes);

				curFeature.svmVote = svmVote;
				curFeature.svmScore = sumRes;

				if(curFeature.svmScore<svmCutoff) blobStatus[b] =-2;
		} // if flagApplyClassifier

		// test if the detection overlapped with real lesions
		NIH_BoneSegmentation_MatchDetections3D(detectionRegion, lesions, numLesions, lesionVoxels, maskImg3D, segPara, matchedLesion, matchedLesionType, matchedLesionSize, matchedOverlap);
		curFeature.matchedLesion = matchedLesion;
		curFeature.matchedLesionType = matchedLesionType;
		curFeature.matchedLesionSize = matchedLesionSize;
		curFeature.matchedOverlap = matchedOverlap;

		if(numDetections3D>=maxDetectionNumber)
		{
			maxDetectionNumber *=2;
			detections3D = (FeatureStructure *) realloc(detections3D, maxDetectionNumber*sizeof(FeatureStructure));
		}
		if(numDetections>=maxDetectionNumber)
		{
			maxDetectionNumber *=2;
			detections = (DetectionStructure *) realloc(detections, maxDetectionNumber*sizeof(DetectionStructure));
		}
		detections3D[numDetections3D] = curFeature;
		numDetections3D++;
		detections[numDetections] = curDetection;
		numDetections++;
	} // for b

	// fill the mask
	for(k=0, z=0; z<sizez; z++)
	{
		for(y=0; y<sizex; y++)
		{
			for(x=0; x<sizex; x++, k++)
			{
				if(binArray3D[k]==1)
				{
					if(blobStatus[blobArray3D[k]-1]==1)
					{
						maskA[k]=maskStruct.mask_boneMetasis;
					}
					else if(blobStatus[blobArray3D[k]-1]==-2 && maskA[k]==maskStruct.mask_boneMetasis)
					{
						maskA[k]=maskStruct.mask_falseDetection;
					}
					else if(blobStatus[blobArray3D[k]-1]==-1 && maskA[k]==maskStruct.mask_boneMetasis)
					{
						maskA[k]=maskStruct.mask_falseDetection;
					}
				} // if binArray
				else if(maskA[k]==maskStruct.mask_boneMetasis)
					maskA[k]=maskStruct.mask_spongyBone;
			} // for x
		} // for y
	} // for z

	// fill in the detection voxel structure
//	if(detectionVoxels!=NULL) delete detectionVoxels;
//	detectionVoxels = new PaintingStruct[numDetections];

	for(k=0; k<blobRank3D.GetSize(); k++)
	{
		detectionVoxels[k].numPaint = 0;
	}
	
	int isBound, curPaint, curBlob;
	for(z=0, k=0; z<sizez; z++)
	{
		for(y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++)
			{
				curBlob = blobArray3D[k]-1;

				if(curBlob>=0)
				{
					if(blobStatus[curBlob]==-1) continue;
					// check if it on boundary
					if(blobArray3D[k+1]==0 || blobArray3D[k-1]==0 || blobArray3D[k+sizex]==0 || blobArray3D[k-sizex]==0)
					{
						isBound = 1;
					}
					else isBound = 0;
					
					// check if new slice
					curPaint = detectionVoxels[curBlob].numPaint-1;
					if(curPaint>=99) continue;		// too many slices, shouldn't be right!!!
					if(detectionVoxels[curBlob].numPaint==0 || 
						detectionVoxels[curBlob].sliceArray[curPaint]!=z+1)
					{
						curPaint = detectionVoxels[curBlob].numPaint;
						detectionVoxels[curBlob].numPaint++;
						detectionVoxels[curBlob].sliceArray[curPaint] = z+1;
						detectionVoxels[curBlob].paint[curPaint].SetSize(0);
						detectionVoxels[curBlob].paintS[curPaint].SetSize(0);
					}
					detectionVoxels[curBlob].paint[curPaint].Add(IntVec2(x,y));
					detectionVoxels[curBlob].paintS[curPaint].Add(isBound);
				}	// if curBlob
			}	// for x
		}	// for y
	}
	

	delete binImg3D;
	delete blobImg3D;
	return CIS_OK;
}
int NIH_BoneSegmentation_DetectMetasis(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   FeatureStructure* &detections2D, int &numDetections2D,
									   LesionStructure *lesions, int numLesions,
									   PaintingStruct *lesionVoxels,
									   bool flagApplyClassifier, bool flagApply3DClassifier, SvmCommittee *svmCommittee, float svmCutoff)
{
	int x, y, z, k;
	int sizex, sizey, sizez, sizexy, sizexyz;
	short *imgA, *maskA;

	sizex = img3D->Num_Cols();
	sizey = img3D->Num_Rows();
	sizez = img3D->Num_Levels();

	imgA = img3D->GetArray();
	maskA = maskImg3D->GetArray();

	sizexyz = sizex*sizey*sizez;
	sizexy = sizex*sizey;

	int count, k2;

	// locate potential bone met detections
	// should move down to metasis detection program
	// using a fix threshold, should use adaptive method or watershed method
	segPara.metasisIntensityThresh = segInfo.avg_cortical_intensity-200;
//	segPara.metasisIntensityThresh = segInfo.avg_spongy_intensity;

	int maxDetectionNumber;
	intVec2DynArray detectionRegion;
	int matchedLesion, matchedLesionType;
	float matchedOverlap, matchedLesionSize;


	int meanIntensity, meanIntensity_90, meanIntensity_10, metaThresh_2D, minIntensity, maxIntensity, intensity;
	int count10, count90;
	int hist[5000];

/*CIS_Array_Image2D_short *binImg2Dcort;
short *binArray2Dcort;

binImg2Dcort = new CIS_Array_Image2D_short(sizex, sizey);
binArray2Dcort = binImg2Dcort->GetArray();*/


	for(z=segInfo.bound1.z; z<=segInfo.bound2.z; z++)
	{
		if(segInfo.diskBound1[z]==-1) continue;

		// only detect the region within vertebra body and spinus process

		// first compute the mean intensity in the vertebra body region
		meanIntensity = meanIntensity_90 =0;
		minIntensity = 5000;
		maxIntensity = 0;
		for(k2=0; k2<5000; k2++) hist[k2]=0;

		/*for(k=0; k<sizexy ;k++) binArray2Dcort[k] = 0;*/
		k = z*sizexy;
	/*	for(y=segInfo.diskBound1[z].y; y<=segInfo.diskBound2[z].y; y++)
		{
			k2 = y*sizex+segInfo.diskBound1[z].x;
			for(x=segInfo.diskBound1[z].x; x<=segInfo.diskBound2[z].x; x++, k2++)
			{
				if(maskA[k+k2]==maskStruct.mask_spongyBone || maskA[k+k2]==maskStruct.mask_corticalBone)
				{
					binArray2Dcort[k2] = 1;
				}
			}
		}
		CIS_IPA_Erode(binImg2Dcort, 2, true);
		k = z*sizexy;
		for(y=segInfo.diskBound1[z].y; y<=segInfo.diskBound2[z].y; y++)
		{
			k2 = y*sizex+segInfo.diskBound1[z].x;
			for(x=segInfo.diskBound1[z].x; x<=segInfo.diskBound2[z].x; x++, k2++)
			{
				if(binArray2Dcort[k2]==0)
				{
					maskA[k+k2]==0;
				}
			}
		}*/

		

		k = z*sizexy;
		count =0;
		for(y=segInfo.diskBound1[z].y; y<=segInfo.diskBound2[z].y; y++)
		{
			k2 = y*sizex+segInfo.diskBound1[z].x;
			for(x=segInfo.diskBound1[z].x; x<=segInfo.diskBound2[z].x; x++, k2++)
			{
				if(maskA[k+k2]==maskStruct.mask_spongyBone || maskA[k+k2]==maskStruct.mask_corticalBone)
				{
					intensity = imgA[k+k2];
					meanIntensity += intensity;
					hist[intensity] ++;
					if(intensity<minIntensity) minIntensity=intensity;
					if(intensity>maxIntensity) maxIntensity=intensity;

					count++;
				}
			}
		}	// for y

		if(count==0) continue;
		meanIntensity /= count;

		count10 = 0;
		meanIntensity_10 = 0;
		for(k2=minIntensity; k2<=maxIntensity; k2++)
		{
			if(count10<count/5)
			{
				count10 += hist[k2];
				meanIntensity_10 += hist[k2]*k2;
			}
			else break;
		}
		if(count10==0) continue;
		meanIntensity_10 /= count10;
		
		count90 = 0;
		meanIntensity_90 = 0;
		for(k2=maxIntensity; k2>=minIntensity; k2--)
		{
			if(count90<count/10)
			{
				count90 += hist[k2];
				meanIntensity_90 += hist[k2]*k2;
			}
			else break;
		}
		if(count90==0) continue;
		meanIntensity_90 /= count90;
		
		if(meanIntensity-100>meanIntensity_10) 
			metaThresh_2D = meanIntensity-100;
		else metaThresh_2D = meanIntensity_10;
		
		//Label as possible metasis if intensity smaller/larger than metaThresh_2D (lytic/sclerotic)
		k = z*sizexy;
		count =0;
		for(y=segInfo.diskBound1[z].y; y<=segInfo.diskBound2[z].y; y++)
		{
			k2 = y*sizex+segInfo.diskBound1[z].x;
			for(x=segInfo.diskBound1[z].x; x<=segInfo.diskBound2[z].x; x++, k2++)
			{
				if((maskA[k+k2]==maskStruct.mask_spongyBone || maskA[k+k2]==maskStruct.mask_corticalBone) 
					&& imgA[k+k2]>metaThresh_2D) //right here
				{
					maskA[k+k2] = maskStruct.mask_boneMetasis;
					count++;
				}
			}
		}
		for(y=segInfo.sprocessBound1[z].y; y<=segInfo.sprocessBound2[z].y; y++)
		{
			k2 = y*sizex+segInfo.sprocessBound1[z].x;
			for(x=segInfo.sprocessBound1[z].x; x<=segInfo.sprocessBound2[z].x; x++, k2++)
			{
				if((maskA[k+k2]==maskStruct.mask_spongyBone || maskA[k+k2]==maskStruct.mask_corticalBone) 
					&& imgA[k+k2]>segPara.metasisIntensityThresh) //right here
				{
					maskA[k+k2] = maskStruct.mask_boneMetasis;
					count++;
				}
			}
		}
if(z==70)
		{
		FILE * pFile;
		pFile = fopen("c:\\tmp\\iosfjiodfj.txt","w");
		int xx, yy;
		for (yy=0; yy<sizey; yy++)
		{
			for (xx=0; xx<sizex; xx++)
			{
				fprintf (pFile, "%d", maskImg3D->FastGet(xx,yy,70));
			}
			fprintf( pFile, "\n");
		}
		fclose (pFile);
		}
		if(count==0) continue;
		
	}

	NIH_BoneSegmentation_AnalyzeMetasis2D(img3D, maskImg3D, maskStruct, segPara, segInfo, detections2D, numDetections2D, 
		lesions, numLesions, lesionVoxels, flagApplyClassifier, flagApply3DClassifier, svmCommittee, svmCutoff, -1, -1);

	return CIS_OK;
}

int NIH_BoneSegmentation_LocateDetection2D(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   FeatureStructure* &detections2D, int &numDetections2D,
									   LesionStructure *lesions, int numLesions,
									   PaintingStruct *lesionVoxels,
									   bool flagApplyClassifier, SvmCommittee *svmCommittee, float svmCutoff)
{
	int x, y, z, k;
	int sizex, sizey, sizez, sizexy, sizexyz;
	short *imgA, *maskA;
	QStringList featureNameList, featureValueList;
	double svmVote, svmScore, sumRes;

	sizex = img3D->Num_Cols();
	sizey = img3D->Num_Rows();
	sizez = img3D->Num_Levels();

	imgA = img3D->GetArray();
	maskA = maskImg3D->GetArray();

	sizexyz = sizex*sizey*sizez;
	sizexy = sizex*sizey;

	CIS_Array_Image2D_short *binImg2D, *blobImg2D;
	short *binArray2D, *blobArray2D;
	intDynArray blobRank2D;
	intVec2DynArray blobCentroid2D;

	int count, k2, b;

	binImg2D = new CIS_Array_Image2D_short(sizex, sizey);
	blobImg2D = new CIS_Array_Image2D_short(sizex, sizey);
	binArray2D = binImg2D->GetArray();
	blobArray2D = blobImg2D->GetArray();
	
	int maxDetectionNumber;
	intVec2DynArray detectionRegion;
	int matchedLesion, matchedLesionType;
	float matchedOverlap, matchedLesionSize;

	FeatureStructure curFeature;
	numDetections2D = 0;
	maxDetectionNumber = sizez*10;	// assume 10 detection per slice

	if(detections2D!=NULL) free(detections2D);
	detections2D = (FeatureStructure *)malloc(maxDetectionNumber*sizeof(FeatureStructure));

	// assemble the feature Name list
	if(flagApplyClassifier)
	{
		featureNameList.push_back("relCoordx");
		featureNameList.push_back("relCoordy");
		featureNameList.push_back("relCoordz");
		featureNameList.push_back("distToBoundary");
		featureNameList.push_back("area");
		featureNameList.push_back("volume");
		featureNameList.push_back("perimeter");
		featureNameList.push_back("primaryAxisLength");
		featureNameList.push_back("secondaryAxisLength");
		featureNameList.push_back("outerBorderRatio");
		featureNameList.push_back("aspectRatio");
		featureNameList.push_back("spherecity");
		featureNameList.push_back("compactness");
		featureNameList.push_back("roundness");
		featureNameList.push_back("shapeComplexity_f1");
		featureNameList.push_back("shapeComplexity_f2");
		featureNameList.push_back("shapeComplexity_f21");
		featureNameList.push_back("meanIntensity");
		featureNameList.push_back("stdevIntensity");
		featureNameList.push_back("skewnessIntensity");
		featureNameList.push_back("kurtosisIntensity");
		featureNameList.push_back("interiorIntensity");
		featureNameList.push_back("borderIntensity");
		featureNameList.push_back("outsideIntensity");
		featureNameList.push_back("outsideIntensityDev");
		featureNameList.push_back("innerOuterContrast");
		featureNameList.push_back("borderThickness");
	}

	for(z=0; z<sizez; z++)
	{
		for(k2=0; k2<sizexy; k2++) binArray2D[k2]=0;

		k = z*sizexy;
		count =0;
		for(y=0, k2=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k2++)
			{
				if(maskA[k+k2]==maskStruct.mask_boneMetasis) 
				{
					binArray2D[k2] = 1;
					count++;
				}
			}
		}

		if(count==0) continue;
		
		NIH_Algo_Blob_Labelling(binImg2D, blobImg2D, blobRank2D, blobCentroid2D, false);

		// compute detection features, and export to a csv file
		//
		for(b=0; b<blobRank2D.GetSize(); b++)
		{
			// fill the detection region
			detectionRegion.SetSize(0);
			for(k2=0, y=0; y<sizey; y++)
			{
				for(x=0; x<sizex; x++, k++, k2++)
				{
					if(blobArray2D[k2]==b+1) 
					{
						detectionRegion.Add(IntVec2(x,y));
					}
				}
			}

			// computer features
			NIH_BoneSegmentation_ComputeMetasisFeatures2D(curFeature, z, detectionRegion, img3D, maskImg3D, maskStruct, segPara, segInfo);

			// pass the detections to a classifier
			if(flagApplyClassifier)
			{
				featureValueList.clear();
				featureValueList.push_back(QString::number(curFeature.relCoordx));
				featureValueList.push_back(QString::number(curFeature.relCoordy));
				featureValueList.push_back(QString::number(curFeature.relCoordz));
				featureValueList.push_back(QString::number(curFeature.distToBoundary));
				featureValueList.push_back(QString::number(curFeature.area));
				featureValueList.push_back(QString::number(curFeature.volume));
				featureValueList.push_back(QString::number(curFeature.perimeter));
				featureValueList.push_back(QString::number(curFeature.primaryAxisLength));
				featureValueList.push_back(QString::number(curFeature.secondaryAxisLength));
				featureValueList.push_back(QString::number(curFeature.outerBorderRatio));
				featureValueList.push_back(QString::number(curFeature.aspectRatio));
				featureValueList.push_back(QString::number(curFeature.spherecity));
				featureValueList.push_back(QString::number(curFeature.compactness));
				featureValueList.push_back(QString::number(curFeature.roundness));
				featureValueList.push_back(QString::number(curFeature.shapeComplexity_f1));
				featureValueList.push_back(QString::number(curFeature.shapeComplexity_f2));
				featureValueList.push_back(QString::number(curFeature.shapeComplexity_f21));
				featureValueList.push_back(QString::number(curFeature.meanIntensity));
				featureValueList.push_back(QString::number(curFeature.stdevIntensity));
				featureValueList.push_back(QString::number(curFeature.skewnessIntensity));
				featureValueList.push_back(QString::number(curFeature.kurtosisIntensity));
				featureValueList.push_back(QString::number(curFeature.interiorIntensity));
				featureValueList.push_back(QString::number(curFeature.borderIntensity));
				featureValueList.push_back(QString::number(curFeature.outsideIntensity));
				featureValueList.push_back(QString::number(curFeature.outsideIntensityDev));
				featureValueList.push_back(QString::number(curFeature.innerOuterContrast));
				featureValueList.push_back(QString::number(curFeature.borderThickness));

//				svmCommittee->Predict(featureNameList, featureValueList, svmVote, svmScore);

				svmCommittee->PredictWithCombinationSumRes(featureNameList, featureValueList, svmVote, svmScore, sumRes);

				curFeature.svmVote = svmVote;
				curFeature.svmScore = sumRes;

			}
			
			// test if the detection overlaped with real lesions
			NIH_BoneSegmentation_MatchDetections2D(detectionRegion, z, lesions, numLesions, lesionVoxels, maskImg3D, segPara, 
				matchedLesion, matchedLesionType, matchedLesionSize, matchedOverlap);
			curFeature.matchedLesion = matchedLesion;
			curFeature.matchedLesionType = matchedLesionType;
			curFeature.matchedLesionSize = matchedLesionSize;
			curFeature.matchedOverlap = matchedOverlap;

			if(numDetections2D>=maxDetectionNumber)
			{
				maxDetectionNumber *=2;
				detections2D = (FeatureStructure *) realloc(detections2D, maxDetectionNumber*sizeof(FeatureStructure));
			}
			detections2D[numDetections2D] = curFeature;
			numDetections2D ++;
		}	// for b

	}	// for z

	delete binImg2D;
	delete blobImg2D;
	return CIS_OK;
}




int NIH_BoneSegmentation_ComputeMetasisFeatures2D(FeatureStructure &features, int detectionSlice, intVec2DynArray &detectionRegion, 
												  CIS_Array_Image3D_short *img3D,
													CIS_Array_Image3D_short *maskImg3D,
													SpineSegmentationMaskStructure &maskStruct,
													SpineSegmentationParameter &segPara,
													SpineSegmentationInfo &segInfo)
{
	if(maskImg3D==NULL || img3D==NULL || detectionRegion.GetSize()==0) return CIS_ERROR;

	// initialize the feature values;
	features.centerx = features.centery = features.centerz = 0;
	features.distToBoundary=0;
	features.relCoordx = features.relCoordy = features.relCoordz=-1;

	features.area = features.volume = features.perimeter = 0;

	features.outerBorderRatio = features.aspectRatio = -1;
	features.spherecity = features.compactness = 0;

	features.roundness = 0;
	features.borderThickness=0;

	features.svmScore = 0;
	features.svmVote = 0;

	// float radial distance measures
	features.shapeComplexity_f1 = 0;
	features.shapeComplexity_f2 = 0;
	features.shapeComplexity_f21 = 0;

	// intensity features
	features.meanIntensity = features.stdevIntensity = features.skewnessIntensity = features.kurtosisIntensity = 0;
	features.interiorIntensity = features.borderIntensity = 0;
	features.outsideIntensity = features.outsideIntensityDev = 0;

	// status
	features.matchedLesion = features.matchedLesionType = -1;
	features.matchedOverlap = 0;
	features.matchedLesionSize = -1;
	
	// get some info about the images
	int x, y, k, k2,i;
	int sizex, sizey, sizez, sizexy, sizexyz;
	float pixelSizex, pixelSizey, pixelSizez;
	short *imgA, *maskA;

	sizex = img3D->Num_Cols();
	sizey = img3D->Num_Rows();
	sizez = img3D->Num_Levels();

	pixelSizex = img3D->Get_Pixel_SizeX()/10;
	pixelSizey = img3D->Get_Pixel_SizeY()/10;
	pixelSizez = img3D->Get_Pixel_SizeZ()/10;

	imgA = img3D->GetArray();
	maskA = maskImg3D->GetArray();

	sizexyz = sizex*sizey*sizez;
	sizexy = sizex*sizey;

	// fill a 2D image with the detection region
	CIS_Array_Image2D_short *filledRegion;
	short *filledArray;
	int bbx0, bbx1, bby0, bby1;		// bounding box for the region
	int regionSize = detectionRegion.GetSize();

	filledRegion = new CIS_Array_Image2D_short(sizex, sizey);
	filledArray = filledRegion->GetArray();

	// fill the array
	bbx0 = sizex;
	bbx1 = 0;
	bby0 = sizey;
	bby1 = 0;
	for(k=0; k<regionSize; k++)
	{
		x = detectionRegion[k].x;
		y = detectionRegion[k].y;
		
		filledArray[x+y*sizex] = 1;
		
		// update the bounding box
		if(x<bbx0) bbx0=x;
		if(x>bbx1) bbx1=x;
		if(y<bby0) bby0=y;
		if(y>bby1) bby1=y;

		// compute centroid
		features.centerx +=x;
		features.centery +=y;
	}

	features.centerx /= (float)regionSize;
	features.centery /= (float)regionSize;
	features.centerz = detectionSlice;

	// compute the relative location to the "origin" of the disk
	if(segInfo.diskBound1.GetSize()>detectionSlice && segInfo.diskBound1[detectionSlice].x!=-1)
	{
		float diskSizex, diskSizey;
		float diskOrgx, diskOrgy; // org= origin! i.e. midpoint

		diskSizex = segInfo.diskBound2[detectionSlice].x-segInfo.diskBound1[detectionSlice].x+1;
		diskSizey = segInfo.diskBound2[detectionSlice].y-segInfo.diskBound1[detectionSlice].y+1;
		diskOrgx = (segInfo.diskBound2[detectionSlice].x+segInfo.diskBound1[detectionSlice].x)/2;
		diskOrgy = segInfo.diskBound2[detectionSlice].y; //so the top of the disk?

		features.relCoordx = (features.centerx-diskOrgx)/diskSizex;
		features.relCoordy = (features.centery-diskOrgy)/diskSizey;

		//relative to the origin of the spinal process
		if(features.relCoordy>0)	
		{
			// appears in the spinal process, use another bounding box to compute coordinates
			diskSizex = segInfo.sprocessBound2[detectionSlice].x-segInfo.sprocessBound1[detectionSlice].x;
			diskSizey = segInfo.sprocessBound2[detectionSlice].y-segInfo.sprocessBound1[detectionSlice].y;
			diskOrgx = (segInfo.sprocessBound2[detectionSlice].x+segInfo.sprocessBound1[detectionSlice].x)/2;
			diskOrgy = segInfo.sprocessBound1[detectionSlice].y;

			features.relCoordx = (features.centerx-diskOrgx)/diskSizex;
			features.relCoordy = (features.centery-diskOrgy)/diskSizey;
		}
	}

	// locate the border 
	int innerBorder, outerBorder, corticalBorder, bboxBorder, airBorder, cordBorder, ribBorder;
	intVec2DynArray borderRegion;

	borderRegion.SetSize(0);
	innerBorder = outerBorder = 0;
	corticalBorder = bboxBorder = airBorder = cordBorder = ribBorder = 0;
	k = detectionSlice*sizexy;
	for(y=bby0; y<=bby1; y++)
	{
		k2 = y*sizex+bbx0;
		for(x=bbx0; x<=bbx1; x++, k2++)
		{
			if(filledArray[k2]==1)
			{
				if(filledArray[k2-sizex]==0 || filledArray[k2+sizex]==0 || filledArray[k2-1]==0 || filledArray[k2+1]==0)
				{
					borderRegion.Add(IntVec2(x, y));
					filledArray[k2]=2;

					features.perimeter += 1;

//					if(x<=bbx0 || x>=bbx1 || y<=bby0 || y>=bby1) outerBorder++;
//					else 
					if((maskA[k+k2-sizex]==maskStruct.mask_body || maskA[k+k2-sizex]==maskStruct.mask_otherBone || maskA[k+k2-sizex]==maskStruct.mask_rib) ||
						(maskA[k+k2+sizex]==maskStruct.mask_body || maskA[k+k2+sizex]==maskStruct.mask_otherBone || maskA[k+k2+sizex]==maskStruct.mask_rib) ||
						(maskA[k+k2-1]==maskStruct.mask_body || maskA[k+k2-1]==maskStruct.mask_otherBone || maskA[k+k2-1]==maskStruct.mask_rib) ||
						(maskA[k+k2+1]==maskStruct.mask_body || maskA[k+k2+1]==maskStruct.mask_otherBone || maskA[k+k2+1]==maskStruct.mask_rib)) outerBorder++;
					else innerBorder++;

					if(segInfo.diskBound1.GetSize()>detectionSlice)
					{
						if((y<segInfo.diskBound2[detectionSlice].y  && 
							(x<=segInfo.diskBound1[detectionSlice].x+1 || x>=segInfo.diskBound2[detectionSlice].x-2 
							|| y<=segInfo.diskBound1[detectionSlice].y)) || y==segInfo.diskBound2[detectionSlice].y ||
							(y>segInfo.diskBound2[detectionSlice].y  && 
							(x<=segInfo.sprocessBound1[detectionSlice].x+1 || x>=segInfo.sprocessBound2[detectionSlice].x-2)))
							bboxBorder++;
					}
					
					if(maskA[k+k2-sizex]==maskStruct.mask_corticalBone ||
						maskA[k+k2+sizex]==maskStruct.mask_corticalBone ||
						maskA[k+k2-1]==maskStruct.mask_corticalBone ||
						maskA[k+k2+1]==maskStruct.mask_corticalBone) corticalBorder++;

					if(maskA[k+k2-sizex]==maskStruct.mask_air ||
						maskA[k+k2+sizex]==maskStruct.mask_air ||
						maskA[k+k2-1]==maskStruct.mask_air ||
						maskA[k+k2+1]==maskStruct.mask_air) airBorder++;

					if(maskA[k+k2-sizex]==maskStruct.mask_rib ||
						maskA[k+k2+sizex]==maskStruct.mask_rib ||
						maskA[k+k2-1]==maskStruct.mask_rib ||
						maskA[k+k2+1]==maskStruct.mask_rib) ribBorder++;

					if(maskA[k+k2-sizex]==maskStruct.mask_spinalCord ||
						maskA[k+k2+sizex]==maskStruct.mask_spinalCord ||
						maskA[k+k2-1]==maskStruct.mask_spinalCord ||
						maskA[k+k2+1]==maskStruct.mask_spinalCord) cordBorder++;
				}
			}
		}	// for x
	}	// for y

	features.perimeter *= pixelSizex;

	features.area = regionSize*pixelSizex*pixelSizey;
	features.volume = features.area*pixelSizez;
	features.outerBorderRatio = (float)outerBorder/(float)(innerBorder+outerBorder);
	features.corticalBorderRatio = (float)corticalBorder/(float)(innerBorder+outerBorder);
	features.bboxBorderRatio = (float)bboxBorder/(float)(innerBorder+outerBorder);
	features.airBorderRatio = (float)airBorder/(float)(innerBorder+outerBorder);
	features.cordBorderRatio = (float)cordBorder/(float)(innerBorder+outerBorder);
	features.ribBorderRatio = (float)ribBorder/(float)(innerBorder+outerBorder);

	// compute aspect ratio
	Vec2 centroid, primaryAxis, secondaryAxis, axisLength;
	CIS_Algo_Contour_GetSecondMoment(borderRegion, centroid, primaryAxis, secondaryAxis, axisLength);
	features.aspectRatio = axisLength.x/axisLength.y;
	features.primaryAxisLength = axisLength.x*pixelSizex;
	features.secondaryAxisLength = axisLength.y*pixelSizex;

	features.compactness = features.perimeter*features.perimeter/features.area;
	features.roundness = (4*features.area)/(PI*features.primaryAxisLength*features.primaryAxisLength);
	features.spherecity = (pow(PI, 1/3)*pow(6*features.volume, 2/3))/(2*features.area+features.perimeter*pixelSizez);
	// sphericity is the ratio of the surface area of a sphere with the same volume as the feature, to the surface area of the feature


	// compute the radial distance measures
	double shapeComplexity_f1, shapeComplexity_f2, shapeComplexity_f21;
	CIS_Algo_ComputeRadialDistanceMeasures(borderRegion, shapeComplexity_f1, shapeComplexity_f2, shapeComplexity_f21);
	features.shapeComplexity_f1 = shapeComplexity_f1;
	features.shapeComplexity_f2 = shapeComplexity_f2;
	features.shapeComplexity_f21 = shapeComplexity_f21;

	
	// locate a ring of outside region (ringThick=2) of the detection (regular spine)
	//
	CIS_Array_Image2D_short *binImg;
	short *binArray;
	binImg = new CIS_Array_Image2D_short(sizex, sizey);
	binArray = binImg->GetArray();
	
	for(k2=0; k2<sizexy; k2++) 
	{
		if(filledArray[k2]!=0) binArray[k2]=1;
		else binArray[k2]=0;
	}

	// dilate 2 to get the outside
	CIS_IPA_Dilate(binImg, 3, false);

	k = detectionSlice*sizexy;
	for(k2=0; k2<sizexy; k2++) 
	{
		if(binArray[k2]==1 && filledArray[k2]==0)
		{
			if(maskA[k+k2]==maskStruct.mask_corticalBone || maskA[k+k2]==maskStruct.mask_spongyBone)
				filledArray[k2] = 3;
		}
	}

	
	doubleDynArray intensityArray;
	double mean, variance, skewness, kurtosis, absoluteDeviation, standardDeviation;

	// intensity of entire region
	intensityArray.SetSize(regionSize);

	k=detectionSlice*sizexy;
	for(i=0; i<regionSize; i++)
	{
		intensityArray[i] = imgA[k+detectionRegion[i].y*sizex+detectionRegion[i].x];
	}

	CIS_Algo_ComputeStatisticalMoments(intensityArray, mean, variance, skewness, kurtosis, absoluteDeviation, standardDeviation);

	features.meanIntensity = mean;
	features.stdevIntensity = standardDeviation;
	features.skewnessIntensity = skewness;
	features.kurtosisIntensity = kurtosis;

	// get the interior and border intensity and outside intensity
	doubleDynArray interiorIntensity, borderIntensity, outsideIntensity;

	k=detectionSlice*sizexy;
	for(k2=0; k2<sizexy; k2++) 
	{
		if(filledArray[k2]==1) interiorIntensity.Add(imgA[k+k2]);
		else if(filledArray[k2]==2) borderIntensity.Add(imgA[k+k2]);
		else if(filledArray[k2]==3) outsideIntensity.Add(imgA[k+k2]);
	}

	CIS_Algo_ComputeStatisticalMoments(interiorIntensity, mean, variance, skewness, kurtosis, absoluteDeviation, standardDeviation);
	features.interiorIntensity = mean;
	CIS_Algo_ComputeStatisticalMoments(borderIntensity, mean, variance, skewness, kurtosis, absoluteDeviation, standardDeviation);
	features.borderIntensity = mean;
	CIS_Algo_ComputeStatisticalMoments(outsideIntensity, mean, variance, skewness, kurtosis, absoluteDeviation, standardDeviation);
	features.outsideIntensity = mean;
	features.outsideIntensityDev = standardDeviation;

	if(features.outsideIntensity!=0) features.innerOuterContrast = features.interiorIntensity/features.outsideIntensity;
	else features.innerOuterContrast=0;
	

	// compute the distance to border
	int minDist, dist;

	k2 = detectionSlice*sizexy;
	minDist = sizex;
	// to right
	for(x=features.centerx, k=features.centery*sizex+features.centerx; x<sizex; x++, k++)
	{
		if(maskA[k2+k]==maskStruct.mask_body || maskA[k2+k]==maskStruct.mask_otherBone || maskA[k2+k]==maskStruct.mask_rib || maskA[k2+k]==maskStruct.mask_air)
			break;
	}
	minDist = x-features.centerx;
	// to left
	if(minDist>features.centerx) minDist=features.centerx;
	k = features.centery*sizex+features.centerx-minDist;
	if(maskA[k2+k]==maskStruct.mask_body || maskA[k2+k]==maskStruct.mask_otherBone || maskA[k2+k]==maskStruct.mask_rib || maskA[k2+k]==maskStruct.mask_air)
	{
		for(x=features.centerx-minDist; x<features.centerx; x++, k++)
		{
			if(maskA[k2+k]==maskStruct.mask_spongyBone || maskA[k2+k]==maskStruct.mask_corticalBone || maskA[k2+k]==maskStruct.mask_boneMetasis)
				break;
		}
		minDist = features.centerx-x;
	}
	// to up
	if(minDist>features.centery) minDist=features.centery;
	k = (features.centery-minDist)*sizex+features.centerx;
	if(maskA[k2+k]==maskStruct.mask_body || maskA[k2+k]==maskStruct.mask_otherBone || maskA[k2+k]==maskStruct.mask_rib || maskA[k2+k]==maskStruct.mask_air)
	{
		for(y=features.centery-minDist; y<features.centery; y++, k+=sizex)
		{
			if(maskA[k2+k]==maskStruct.mask_spongyBone || maskA[k2+k]==maskStruct.mask_corticalBone || maskA[k2+k]==maskStruct.mask_boneMetasis)
				break;
		}
		minDist = features.centery-y;
	}
	// to down
	if(minDist>sizey-1-features.centery) minDist=sizey-1-features.centery;
	k = (features.centery+minDist)*sizex+features.centerx;
	if(maskA[k2+k]==maskStruct.mask_body || maskA[k2+k]==maskStruct.mask_otherBone || maskA[k2+k]==maskStruct.mask_rib || maskA[k2+k]==maskStruct.mask_air)
	{
		for(y=features.centery+minDist; y>features.centery; y--, k-=sizex)
		{
			if(maskA[k2+k]==maskStruct.mask_spongyBone || maskA[k2+k]==maskStruct.mask_corticalBone || maskA[k2+k]==maskStruct.mask_boneMetasis)
				break;
		}
		minDist = y-features.centery;
	}
	features.distToBoundary = pixelSizex*minDist;

	// compute the neighborhood info
	int neighborI, neighborA;
	neighborA = neighborI = 0;
	if(segInfo.diskBound1.GetSize()>detectionSlice && segInfo.diskBound1[detectionSlice].x!=-1)
	{
		if(features.relCoordy<0)
		{
			for(y=segInfo.diskBound1[detectionSlice].y; y<=segInfo.diskBound2[detectionSlice].y; y++)
			{
				k=y*sizex+segInfo.diskBound1[detectionSlice].x;
				for(x=segInfo.diskBound1[detectionSlice].x; x<=segInfo.diskBound2[detectionSlice].x; x++, k++)
				{
					if(maskA[k2+k]==maskStruct.mask_spongyBone || maskA[k2+k]==maskStruct.mask_corticalBone) // || maskA[k2+k]==maskStruct.mask_boneMetasis)
					{
						neighborI += imgA[k2+k];
						neighborA ++;
					}
				}
			}
		}
		else
		{
			for(y=segInfo.sprocessBound1[detectionSlice].y; y<=segInfo.sprocessBound2[detectionSlice].y; y++)
			{
				k=y*sizex+segInfo.sprocessBound1[detectionSlice].x;
				for(x=segInfo.sprocessBound1[detectionSlice].x; x<segInfo.sprocessBound2[detectionSlice].x; x++, k++)
				{
					if(maskA[k2+k]==maskStruct.mask_spongyBone || maskA[k2+k]==maskStruct.mask_corticalBone) //  || maskA[k2+k]==maskStruct.mask_boneMetasis)
					{
						neighborI += imgA[k2+k];
						neighborA ++;
					}
				}
			}
		}	// else features
	}	// if segInfo

	if(neighborA!=0)
	{
		features.neighborArea = neighborA*pixelSizex*pixelSizey;
		features.neighborIntensity = neighborI/neighborA;
	}
	else
	{
		features.neighborArea = 0;
		features.neighborIntensity = 0;
	}

	delete filledRegion;
	delete binImg;

	return CIS_OK;
}


int NIH_BoneSegmentation_ComputeMetasisFeatures3D(FeatureStructure &features3D,
												  intVec3DynArray &detectionRegion,
												  CIS_Array_Image3D_short *img3D,
												  CIS_Array_Image3D_short *maskImg3D,
												  SpineSegmentationMaskStructure &maskStruct,
												  SpineSegmentationParameter &segPara,
												  SpineSegmentationInfo &segInfo)
{
	if(maskImg3D==NULL || img3D==NULL || detectionRegion.GetSize()==0) return CIS_ERROR;

	// initialize the feature values
	features3D.centerx = features3D.centery = features3D.centerz = 0;
	features3D.distToBoundary = 0;
	features3D.relCoordx = features3D.relCoordy = features3D.relCoordz = -1;

	features3D.volume = features3D.surfaceArea = 0;

	features3D.outerBorderRatio = features3D.aspectRatio = -1;
	features3D.spherecity = features3D.compactness = 0;
	
	features3D.borderThickness = 0;

	features3D.svmScore = 0;
	features3D.svmVote = 0;

	// float radial distance measures
	features3D.shapeComplexity_f1 = 0;
	features3D.shapeComplexity_f2 = 0;
	features3D.shapeComplexity_f21 = 0;

	//intensity features
	features3D.meanIntensity = features3D.stdevIntensity = features3D.skewnessIntensity = features3D.kurtosisIntensity = 0;
	features3D.interiorIntensity = features3D.borderIntensity = 0;
	features3D.outsideIntensity = features3D.outsideIntensityDev = 0;

	//status
	features3D.matchedLesion = features3D.matchedLesionType = -1;
	features3D.matchedOverlap = 0;
	features3D.matchedLesionSize = -1;

	// get some info about the images
	int x, y, z, k, k2, i, j;
	int sizex, sizey, sizez, sizexy, sizexyz;
	short *imgA, *maskA;
	float pixelSizex, pixelSizey, pixelSizez;

	pixelSizex = img3D->Get_Pixel_SizeX()/10;
	pixelSizey = img3D->Get_Pixel_SizeY()/10;
	pixelSizez = img3D->Get_Pixel_SizeZ()/10;

	sizex = img3D->Num_Cols();
	sizey = img3D->Num_Rows();
	sizez = img3D->Num_Levels();

	imgA = img3D->GetArray();
	maskA = maskImg3D->GetArray();

	sizexyz = sizex*sizey*sizez;
	sizexy = sizex*sizey;

	double cov_mat[3][3];
	IntVec3 centroid, v;
	CIS_Matrix_JY_double cov_mat1(3,3), e_vector(3,3);
	CIS_Vector_JY_double e_value(3);
	int nrot;
	float aspectRatio10, aspectRatio20, aspectRatio21;

	// fill a 3D image with the detection region
	CIS_Array_Image3D_short *filledRegion;
	short *filledArray;
	int bbx0, bbx1, bby0, bby1, bbz0, bbz1; // bounding box for the region
	int regionSize = detectionRegion.GetSize();

	filledRegion = new CIS_Array_Image3D_short(sizex, sizey, sizez);
	filledArray = filledRegion->GetArray();

	// fill the array
	bbx0 = sizex;
	bbx1 = 0;
	bby0 = sizey;
	bby1 = 0;
	bbz0 = sizez;
	bbz1 = 0;
	for(k=0; k<regionSize; k++)
	{
		x = detectionRegion[k].x;
		y = detectionRegion[k].y;
		z = detectionRegion[k].z;

		filledArray[x+y*sizex+z*sizexy] = 1;

		// update the bounding box
		if(x<bbx0) bbx0=x;
		if(x>bbx1) bbx1=x;
		if(y<bby0) bby0=y;
		if(y>bby1) bby1=y;
		if(z<bbz0) bbz0=z;
		if(z>bbz1) bbz1=z;

		// compute centroid
		features3D.centerx +=x;
		features3D.centery +=y;
		features3D.centerz +=z;
	}

	features3D.centerx /= (float)regionSize;
	features3D.centery /= (float)regionSize;
	features3D.centerz /= (float)regionSize;
	centroid.x = features3D.centerx;
	centroid.y = features3D.centery;
	centroid.z = features3D.centerz;

	// compute the relative location (relCoordx, relCoordy) to the "origin" of the disk at centerz
	if(segInfo.diskBound1.GetSize()>features3D.centerz && segInfo.diskBound1[features3D.centerz].x!=-1)
	{
		float diskSizex, diskSizey;
		float diskOrgx, diskOrgy;

		diskSizex = segInfo.diskBound2[features3D.centerz].x-segInfo.diskBound1[features3D.centerz].x+1;
		diskSizey = segInfo.diskBound2[features3D.centerz].y-segInfo.diskBound1[features3D.centerz].y+1;
		diskOrgx = (segInfo.diskBound2[features3D.centerz].x+segInfo.diskBound1[features3D.centerz].x)/2;
		diskOrgy = segInfo.diskBound2[features3D.centerz].y;

		features3D.relCoordx = (features3D.centerx-diskOrgx)/diskSizex;
		features3D.relCoordy = (features3D.centery-diskOrgy)/diskSizey;
	
		//relative to the origin of the spinal process
		if(features3D.relCoordy>0)
		{
			//appears in the spinal process, use another bounding box to compute coordinates
			diskSizex = segInfo.sprocessBound2[features3D.centerz].x-segInfo.sprocessBound1[features3D.centerz].x;
			diskSizey = segInfo.sprocessBound2[features3D.centerz].y-segInfo.sprocessBound1[features3D.centerz].y;
			diskOrgx = (segInfo.sprocessBound2[features3D.centerz].x+segInfo.sprocessBound1[features3D.centerz].x)/2;
			diskOrgy = segInfo.sprocessBound1[features3D.centerz].y;

			features3D.relCoordx = (features3D.centerx-diskOrgx)/diskSizex;
			features3D.relCoordy = (features3D.centery-diskOrgy)/diskSizey;
		}
	}

	// locate the border
	int innerBorder, outerBorder, corticalBorder, bboxBorder, airBorder, cordBorder, ribBorder;
	intVec3DynArray borderRegion;

	borderRegion.SetSize(0);
	innerBorder = outerBorder = 0;
	corticalBorder = bboxBorder = airBorder = cordBorder = ribBorder = 0;

	// to avoid problems at the z-edges
	if(bbz0==0) bbz0=1;
	if(bbz1==sizez) bbz1=sizez-1;

	for(z=bbz0; z<=bbz1; z++)
	{
		k = z*sizexy;
		for(y=bby0; y<=bby1; y++)
		{
			k2 = y*sizex+bbx0;
			for(x=bbx0; x<bbx1; x++, k2++)
			{
				if(filledArray[k+k2]==1)
				{
					features3D.volume +=1;

					if(filledArray[k+k2-sizex]==0 || filledArray[k+k2+sizex]==0 ||
						filledArray[k+k2-1]==0 || filledArray[k+k2+1]==0 ||
						filledArray[k+k2-sizexy]==0 || filledArray[k+k2+sizexy]==0)
					{
						borderRegion.Add(IntVec3(x, y, z));
						filledArray[k+k2]=2;

						features3D.surfaceArea += 1;

						if(x<=bbx0 || x>=bbx1 || y<=bby0 || y>=bby1) outerBorder++;
						else if((maskA[k+k2-sizex]==maskStruct.mask_body || maskA[k+k2-sizex]==maskStruct.mask_otherBone || maskA[k+k2-sizex]==maskStruct.mask_rib) ||
							(maskA[k+k2+sizex]==maskStruct.mask_body || maskA[k+k2+sizex]==maskStruct.mask_otherBone || maskA[k+k2+sizex]==maskStruct.mask_rib) ||
							(maskA[k+k2-1]==maskStruct.mask_body || maskA[k+k2-1]==maskStruct.mask_otherBone || maskA[k+k2-1]==maskStruct.mask_rib) ||
							(maskA[k+k2+1]==maskStruct.mask_body || maskA[k+k2+1]==maskStruct.mask_otherBone || maskA[k+k2+1]==maskStruct.mask_rib) ||
							(maskA[k+k2-sizexy]==maskStruct.mask_body || maskA[k+k2-sizexy]==maskStruct.mask_otherBone || maskA[k+k2-sizexy]==maskStruct.mask_rib) ||
							(maskA[k+k2+sizexy]==maskStruct.mask_body || maskA[k+k2+sizexy]==maskStruct.mask_otherBone || maskA[k+k2+sizexy]==maskStruct.mask_rib))
							outerBorder++;
						else innerBorder++;

						if(segInfo.diskBound1.GetSize()>z)
						{
							if((y<segInfo.diskBound2[z].y && 
								(x<=segInfo.diskBound1[z].x+1 || x>=segInfo.diskBound2[z].x-2
								|| y<=segInfo.diskBound1[z].y)) || y==segInfo.diskBound2[z].y ||
								(y>segInfo.diskBound2[z].y &&
								(x<=segInfo.sprocessBound1[z].x+1 || x>=segInfo.sprocessBound2[z].x-2)))
								bboxBorder++;
						}

						if(maskA[k+k2-sizex]==maskStruct.mask_corticalBone ||
							maskA[k+k2+sizex]==maskStruct.mask_corticalBone ||
							maskA[k+k2-1]==maskStruct.mask_corticalBone ||
							maskA[k+k2+1]==maskStruct.mask_corticalBone ||
							maskA[k+k2-sizexy]==maskStruct.mask_corticalBone ||
							maskA[k+k2+sizexy]==maskStruct.mask_corticalBone) corticalBorder++;

						if(maskA[k+k2-sizex]==maskStruct.mask_air ||
							maskA[k+k2+sizex]==maskStruct.mask_air ||
							maskA[k+k2-1]==maskStruct.mask_air ||
							maskA[k+k2+1]==maskStruct.mask_air ||
							maskA[k+k2-sizexy]==maskStruct.mask_air ||
							maskA[k+k2+sizexy]==maskStruct.mask_air) airBorder++;

						if(maskA[k+k2-sizex]==maskStruct.mask_rib ||
							maskA[k+k2+sizex]==maskStruct.mask_rib ||
							maskA[k+k2-1]==maskStruct.mask_rib ||
							maskA[k+k2+1]==maskStruct.mask_rib ||
							maskA[k+k2-sizexy]==maskStruct.mask_rib ||
							maskA[k+k2+sizexy]==maskStruct.mask_rib) ribBorder++;

						if(maskA[k+k2-sizex]==maskStruct.mask_spinalCord ||
							maskA[k+k2+sizex]==maskStruct.mask_spinalCord ||
							maskA[k+k2-1]==maskStruct.mask_spinalCord ||
							maskA[k+k2+1]==maskStruct.mask_spinalCord ||
							maskA[k+k2-sizexy]==maskStruct.mask_spinalCord ||
							maskA[k+k2+sizexy]==maskStruct.mask_spinalCord) cordBorder++;
					} // if on border
				} // if filled
			} // for x
		} // for y
	} // for z

	features3D.surfaceArea *= pixelSizex*pixelSizey;
	features3D.volume *= pixelSizex*pixelSizey*pixelSizez;

	features3D.outerBorderRatio = (float)outerBorder/(float)(innerBorder+outerBorder);
	features3D.corticalBorderRatio = (float)corticalBorder/(float)(innerBorder+outerBorder);
	features3D.bboxBorderRatio = (float)bboxBorder/(float)(innerBorder+outerBorder);
	features3D.airBorderRatio = (float)airBorder/(float)(innerBorder+outerBorder);
	features3D.cordBorderRatio = (float)cordBorder/(float)(innerBorder+outerBorder);
	features3D.ribBorderRatio = (float)ribBorder/(float)(innerBorder+outerBorder);

	for(i=0; i<3; i++)
		for(j=0; j<3;j++) cov_mat[i][j] = 0;

	for(i=0; i<detectionRegion.GetSize(); i++)
	{
		v = detectionRegion[i];
		cov_mat[0][0] += (v.x-centroid.x)*(v.x-centroid.x);
		cov_mat[0][1] += (v.x-centroid.x)*(v.y-centroid.y);
		cov_mat[0][2] += (v.x-centroid.x)*(v.z-centroid.z);

		cov_mat[1][0] += (v.y-centroid.y)*(v.x-centroid.x);
		cov_mat[1][1] += (v.y-centroid.y)*(v.y-centroid.y);
		cov_mat[1][2] += (v.y-centroid.y)*(v.z-centroid.z);

		cov_mat[2][0] += (v.z-centroid.z)*(v.x-centroid.x);
		cov_mat[2][1] += (v.z-centroid.z)*(v.y-centroid.y);
		cov_mat[2][2] += (v.z-centroid.z)*(v.z-centroid.z);
	}

	for(i=0;i<3;i++) 
	{
		for(j=0;j<3;j++) 
		{
			cov_mat1(i,j) = cov_mat[i][j]/detectionRegion.GetSize();
		}
	}

	NR_Eigen_jacobi(cov_mat1, 3, e_value, e_vector, &nrot);
	for(i=0; i<3; i++) e_value(i) = sqrt(e_value(i))*3;

	aspectRatio10 = e_value(1)/e_value(0);
	aspectRatio20 = e_value(2)/e_value(0);
	aspectRatio21 = e_value(2)/e_value(1);

	features3D.aspectRatio = aspectRatio20;
	features3D.primaryAxisLength = e_value(0)*pixelSizex;
	features3D.secondaryAxisLength = e_value(2)*pixelSizex;

	features3D.spherecity = (pow(PI, 1/3)*pow(6*features3D.volume, 2/3))/(features3D.surfaceArea);

	// compute the radial distance measures
	double shapeComplexity_f1, shapeComplexity_f2, shapeComplexity_f21;
	CIS_Algo_ComputeRadialDistanceMeasures3D(borderRegion, shapeComplexity_f1, shapeComplexity_f2, shapeComplexity_f21);
	features3D.shapeComplexity_f1 = shapeComplexity_f1;
	features3D.shapeComplexity_f2 = shapeComplexity_f2;
	features3D.shapeComplexity_f21 = shapeComplexity_f21;

	// locate a shell of outside region (shellThick = 2) of the detection
	//
	CIS_Array_Image3D_short *binImg;
	short *binArray;
	binImg = new CIS_Array_Image3D_short(sizex, sizey, sizez);
	binArray = binImg->GetArray();

	for(k2=0; k2<sizexyz; k2++)
	{
		if(filledArray[k2]!=0) binArray[k2]=1;
		else binArray[k2]=0;
	}

	// dilate 2 to get the outside
	CIS_IPA_3D_Dilate(binImg, 3, false);

	for(k2=0; k2<sizexyz; k2++)
	{
		if(binArray[k2]==1 && filledArray[k2]==0)
		{
			if(maskA[k2]==maskStruct.mask_corticalBone || maskA[k2]==maskStruct.mask_spongyBone)
				filledArray[k2] = 3;
		}
	}

	doubleDynArray intensityArray;
	double mean, variance, skewness, kurtosis, absoluteDeviation, standardDeviation;

	// intensity of entire region
	intensityArray.SetSize(regionSize);

	for(i=0; i<regionSize; i++)
	{
		intensityArray[i] = imgA[detectionRegion[i].z*sizexy+detectionRegion[i].y*sizex+detectionRegion[i].x];
	}

	CIS_Algo_ComputeStatisticalMoments(intensityArray, mean, variance, skewness, kurtosis, absoluteDeviation, standardDeviation);

	features3D.meanIntensity = mean;
	features3D.stdevIntensity = standardDeviation;
	features3D.skewnessIntensity = skewness;
	features3D.kurtosisIntensity = kurtosis;

	// get the interior and border intensity and outside intensity
	doubleDynArray interiorIntensity, borderIntensity, outsideIntensity;

	for(k2=0; k2<sizexyz; k2++) 
	{
		if(filledArray[k2]==1) interiorIntensity.Add(imgA[k2]);
		else if(filledArray[k2]==2) borderIntensity.Add(imgA[k2]);
		else if(filledArray[k2]==3) outsideIntensity.Add(imgA[k2]);
	}

	CIS_Algo_ComputeStatisticalMoments(interiorIntensity, mean, variance, skewness, kurtosis, absoluteDeviation, standardDeviation);
	features3D.interiorIntensity = mean;
	CIS_Algo_ComputeStatisticalMoments(borderIntensity, mean, variance, skewness, kurtosis, absoluteDeviation, standardDeviation);
	features3D.borderIntensity = mean;
	CIS_Algo_ComputeStatisticalMoments(outsideIntensity, mean, variance, skewness, kurtosis, absoluteDeviation, standardDeviation);
	features3D.outsideIntensity = mean;
	features3D.outsideIntensityDev = standardDeviation;

	if(features3D.outsideIntensity!=0) features3D.innerOuterContrast = features3D.interiorIntensity/features3D.outsideIntensity;
	else features3D.innerOuterContrast=0;

	// compute the distance to border
	int minDist, dist;

	minDist = sizex;
	// to right
	for(x=features3D.centerx, k=features3D.centerz*sizexy+features3D.centery*sizex+features3D.centerx; x<sizex; x++, k++)
	{
		if(maskA[k]==maskStruct.mask_body || maskA[k]==maskStruct.mask_otherBone || maskA[k]==maskStruct.mask_rib || maskA[k]==maskStruct.mask_air)
			break;
	}
	minDist = x-features3D.centerx;
	// to left
	if(minDist>features3D.centerx) minDist=features3D.centerx;
	k = features3D.centerz*sizexy+features3D.centery*sizex+features3D.centerx-minDist;
	if(maskA[k]==maskStruct.mask_body || maskA[k]==maskStruct.mask_otherBone || maskA[k]==maskStruct.mask_rib || maskA[k]==maskStruct.mask_air)
	{
		for(x=features3D.centerx-minDist; x<features3D.centerx; x++, k++)
		{
			if(maskA[k]==maskStruct.mask_spongyBone || maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_boneMetasis)
				break;
		}
		minDist = features3D.centerx-x;
	}
	// to up
	if(minDist>features3D.centery) minDist=features3D.centery;
	k = features3D.centerz*sizexy+(features3D.centery-minDist)*sizex+features3D.centerx;
	if(maskA[k]==maskStruct.mask_body || maskA[k]==maskStruct.mask_otherBone || maskA[k]==maskStruct.mask_rib || maskA[k]==maskStruct.mask_air)
	{
		for(y=features3D.centery-minDist; y<features3D.centery; y++, k+=sizex)
		{
			if(maskA[k]==maskStruct.mask_spongyBone || maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_boneMetasis)
				break;
		}
		minDist = features3D.centery-y;
	}
	// to down
	if(minDist>sizey-1-features3D.centery) minDist=sizey-1-features3D.centery;
	k = features3D.centerz*sizexy+(features3D.centery+minDist)*sizex+features3D.centerx;
	if(maskA[k]==maskStruct.mask_body || maskA[k]==maskStruct.mask_otherBone || maskA[k]==maskStruct.mask_rib || maskA[k]==maskStruct.mask_air)
	{
		for(y=features3D.centery+minDist; y>features3D.centery; y--, k-=sizex)
		{
			if(maskA[k]==maskStruct.mask_spongyBone || maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_boneMetasis)
				break;
		}
		minDist = y-features3D.centery;
	}
	// to top
	if(minDist>sizez-1-features3D.centerz) minDist=sizez-1-features3D.centerz;
	k = (features3D.centerz+minDist)*sizexy+features3D.centery*sizex+features3D.centerx;
	if(maskA[k]==maskStruct.mask_body || maskA[k]==maskStruct.mask_otherBone || maskA[k]==maskStruct.mask_rib || maskA[k]==maskStruct.mask_air)
	{
		for(z=features3D.centerz+minDist; z>features3D.centerz; z--, k-=sizexy)
		{
			if(maskA[k]==maskStruct.mask_spongyBone || maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_boneMetasis)
				break;
		}
		minDist = z-features3D.centerz;
	}
	// to bottom
	if(minDist>features3D.centerz) minDist=features3D.centerz;
	k = (features3D.centerz-minDist)*sizexy+features3D.centery*sizex+features3D.centerx;
	if(maskA[k]==maskStruct.mask_body || maskA[k]==maskStruct.mask_otherBone || maskA[k]==maskStruct.mask_rib || maskA[k]==maskStruct.mask_air)
	{
		for(z=features3D.centerz-minDist; z<features3D.centerz; z++, k+=sizexy)
		{
			if(maskA[k]==maskStruct.mask_spongyBone || maskA[k]==maskStruct.mask_corticalBone || maskA[k]==maskStruct.mask_boneMetasis)
				break;
		}
		minDist = features3D.centerz-z;
	}
	features3D.distToBoundary = pixelSizex*minDist;

	// compute the neighborhood info
	int neighborI, neighborA;
	neighborA = neighborI = 0;
	if(segInfo.diskBound1.GetSize()>features3D.centerz && segInfo.diskBound1[features3D.centerz].x!=-1)
	{
		for(z=bbz0; z<bbz1; z++)
		{
			if(features3D.relCoordy<0)
			{
				for(y=segInfo.diskBound1[z].y; y<=segInfo.diskBound2[z].y; y++)
				{
					k=z*sizexy+y*sizex+segInfo.diskBound1[z].x;
					for(x=segInfo.diskBound1[z].x; x<=segInfo.diskBound2[z].x; x++, k++)
					{
						if(maskA[k]==maskStruct.mask_spongyBone || maskA[k]==maskStruct.mask_corticalBone) // || maskA[k]==maskStruct.mask_boneMetasis)
						{
							neighborI += imgA[k];
							neighborA ++;
						}
					}
				}
			}
			else
			{
				for(y=segInfo.sprocessBound1[z].y; y<=segInfo.sprocessBound2[z].y; y++)
				{
					k=z*sizexy+y*sizex+segInfo.sprocessBound1[z].x;
					for(x=segInfo.sprocessBound1[z].x; x<segInfo.sprocessBound2[z].x; x++, k++)
					{
						if(maskA[k]==maskStruct.mask_spongyBone || maskA[k]==maskStruct.mask_corticalBone) //  || maskA[k]==maskStruct.mask_boneMetasis)
						{
							neighborI += imgA[k];
							neighborA ++;
						}
					}
				}
			}	// else features
		} // for z
	}	// if segInfo

	if(neighborA!=0)
	{
		features3D.neighborArea = neighborA*pixelSizex*pixelSizey*pixelSizez; // actually volume for 3D case
		features3D.neighborIntensity = neighborI/neighborA;
	}
	else
	{
		features3D.neighborArea = 0;
		features3D.neighborIntensity = 0;
	}

	delete filledRegion;
	delete binImg;

	return CIS_OK;
}
int NIH_BoneSegmentation_MergeMetasis(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   DetectionStructure* &detections, int &numDetections,
									   PaintingStruct* detectionVoxels)
{
	if(maskImg3D==NULL || img3D==NULL) return CIS_ERROR;

	int x, y, z, k;
	int sizex, sizey, sizez, sizexy, sizexyz;
	short *imgA, *maskA;
	float px, py, pz;

	px = img3D->Get_Pixel_SizeX()/10;
	py = img3D->Get_Pixel_SizeY()/10;
	pz = img3D->Get_Pixel_SizeZ()/10;

	sizex = img3D->Num_Cols();
	sizey = img3D->Num_Rows();
	sizez = img3D->Num_Levels();

	imgA = img3D->GetArray();
	maskA = maskImg3D->GetArray();

	sizexyz = sizex*sizey*sizez;
	sizexy = sizex*sizey;

	CIS_Array_Image3D_short *binImg;
	CIS_Array_Image3D_short *blobImg;
	short *binArray;
	short *blobArray;

	binImg = new CIS_Array_Image3D_short(sizex, sizey, sizez);
	blobImg = new CIS_Array_Image3D_short(sizex, sizey, sizez);
	binArray = binImg->GetArray();
	blobArray = blobImg->GetArray();

	// initialization
	for(k=0; k<sizexyz; k++)
	{
		if(maskA[k]==maskStruct.mask_boneMetasis) binArray[k]=1;
		else binArray[k]=0;
	}

	intDynArray blobRank;
	intVec3DynArray blobCentroid;

	NIH_Algo_Blob_Labelling_3D(binImg, blobImg, blobRank, blobCentroid, false);

	if(detections!=NULL)
	{
		delete detections;
		detections = NULL;
		numDetections = 0;
	}
	numDetections = blobRank.GetSize();
	detections = new DetectionStructure[numDetections];

	double cov_mat[3][3];
	Vec3 centroid, v;
	vec3DynArray pixelA, detectionRegion;
	int curBlob, i, j;
	float fx, fy, fz;
	CIS_Matrix_JY_double cov_mat1(3,3), e_vector(3,3);
	CIS_Vector_JY_double e_value(3);
	int nrot;
	float aspectRatio10, aspectRatio20, aspectRatio21;

	// tmp commented out to speed up
	for(int b=0; b<numDetections; b++)
	{
		detections[b].detectionLoc = blobCentroid[b];
		detections[b].detectionVolume = blobRank[b]*px*py*pz;
		detections[b].detectionType = -1;
		detections[b].matchedLesion = -1;

		// calculate the detection size
		curBlob = b+1;
		pixelA.SetSize(0);
		detectionRegion.SetSize(0);
		for(z=0, k=0, fz=0; z<sizez; z++, fz+=pz)
		{
			for(y=0, fy=0; y<sizey; y++, fy+=py)
			{
				for(x=0, fx=0; x<sizex; x++, k++, fx+=px)
				{
					if(blobArray[k]==curBlob)
					{
						pixelA.Add(Vec3(fx, fy, fz));
						detectionRegion.Add(Vec3(x, y, z));
					}
				}
			}
		}	// for z

	
		for(i=0;i<3;i++) 
			for(j=0;j<3;j++) cov_mat[i][j] = 0;

		centroid = Vec3(0,0,0);
		for(i=0;i<pixelA.GetSize();i++)
		{
			centroid += pixelA[i];
		}
		if(pixelA.GetSize()>0) centroid /= (double)(pixelA.GetSize());

		for(i=0;i<pixelA.GetSize();i++)
		{
			v = pixelA[i];
			cov_mat[0][0] += (v.x-centroid.x)*(v.x-centroid.x);
			cov_mat[0][1] += (v.x-centroid.x)*(v.y-centroid.y);
			cov_mat[0][2] += (v.x-centroid.x)*(v.z-centroid.z);

			cov_mat[1][0] += (v.y-centroid.y)*(v.x-centroid.x);
			cov_mat[1][1] += (v.y-centroid.y)*(v.y-centroid.y);
			cov_mat[1][2] += (v.y-centroid.y)*(v.z-centroid.z);

			cov_mat[2][0] += (v.z-centroid.z)*(v.x-centroid.x);
			cov_mat[2][1] += (v.z-centroid.z)*(v.y-centroid.y);
			cov_mat[2][2] += (v.z-centroid.z)*(v.z-centroid.z);
		}
	
		for(i=0;i<3;i++) 
		{
			for(j=0;j<3;j++) 
			{
				cov_mat1(i,j) = cov_mat[i][j]/pixelA.GetSize();
			}
		}

		NR_Eigen_jacobi(cov_mat1, 3, e_value, e_vector, &nrot);
		for(i=0; i<3; i++) e_value(i) = sqrt(e_value(i))*3;

		aspectRatio10 = e_value(1)/e_value(0);
		aspectRatio20 = e_value(2)/e_value(0);
		aspectRatio21 = e_value(2)/e_value(1);

		detections[b].detectionSize = e_value(0);
	} // for b



	// fill in the detection voxel structure
//	if(detectionVoxels!=NULL) delete detectionVoxels;
//	detectionVoxels = new PaintingStruct[numDetections];

	for(k=0; k<numDetections; k++)
	{
		detectionVoxels[k].numPaint = 0;
	}
	
	int isBound, curPaint;
	for(z=0, k=0; z<sizez; z++)
	{
		for(y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++)
			{
				curBlob = blobArray[k]-1;

				if(curBlob>=0)
				{
					// check if it on boundary
					if(blobArray[k+1]==0 || blobArray[k-1]==0 || blobArray[k+sizex]==0 || blobArray[k-sizex]==0)
					{
						isBound = 1;
					}
					else isBound = 0;
					
					// check if new slice
					curPaint = detectionVoxels[curBlob].numPaint-1;
					if(curPaint>=99) continue;		// too many slices, shouldn't be right!!!
					if(detectionVoxels[curBlob].numPaint==0 || 
						detectionVoxels[curBlob].sliceArray[curPaint]!=z+1)
					{
						curPaint = detectionVoxels[curBlob].numPaint;
						detectionVoxels[curBlob].numPaint++;
						detectionVoxels[curBlob].sliceArray[curPaint] = z+1;
						detectionVoxels[curBlob].paint[curPaint].SetSize(0);
						detectionVoxels[curBlob].paintS[curPaint].SetSize(0);
					}
					detectionVoxels[curBlob].paint[curPaint].Add(IntVec2(x,y));
					detectionVoxels[curBlob].paintS[curPaint].Add(isBound);
				}	// if curBlob
			}	// for x
		}	// for y
	}	// for z

	delete blobImg;
	delete binImg;

	return CIS_OK;
}


int NIH_BoneSegmentation_MatchDetections(CIS_Array_Image3D_short *maskImg3D,
									   DetectionStructure *detections, int numDetections,
									   PaintingStruct *detectionVoxels,
									   LesionStructure *lesions, int numLesions,
									   PaintingStruct *lesionVoxels, 
									   int &totalMet, int &totalLytic,
									   int &foundMet, int &foundLytic,
									   int &totalDetection, int &tpDetection, int &fpDetection,
									   int &tpLyticDetection, int &fpLyticDetection,
									   int &totalMet10, int &totalLytic10,
									   int &foundMet10, int &foundLytic10,
									   float lesionSizeCutoff)
{
	if(maskImg3D==NULL) return CIS_ERROR;
	
	// first assemble a blob image
	
	CIS_Array_Image3D_short *detectionBlobMap;
	int x, y, z, k;
	int sizex, sizey, sizez, sizexy, sizexyz;
	int i, j, n;

	sizex = maskImg3D->Num_Cols();
	sizey = maskImg3D->Num_Rows();
	sizez = maskImg3D->Num_Levels();

	sizexyz = sizex*sizey*sizez;
	sizexy = sizex*sizey;

	detectionBlobMap = new CIS_Array_Image3D_short(sizex, sizey, sizez);
	short *blobArray = detectionBlobMap->GetArray();

	// fill in the detectionBlobMap
	for(i=0; i<numDetections; i++)
	{
		for(j=0; j<detectionVoxels[i].numPaint; j++)
		{
			z = detectionVoxels[i].sliceArray[j]-1;
			for(n=0; n<detectionVoxels[i].paint[j].GetSize(); n++)
			{
				x = detectionVoxels[i].paint[j][n].x;
				y = detectionVoxels[i].paint[j][n].y;
				k = z*sizexy+y*sizex+x;
				blobArray[k] = i+1;
			}
		}
	}	// for i

	foundMet=0; foundLytic=0;
	totalMet=0; totalLytic=0;

	foundMet10=0; foundLytic10=0;
	totalMet10=0; totalLytic10=0;

	bool found;
	for(k=0; k<numDetections; k++) 
	{
		detections[k].matchedLesion = -1;
		detections[k].detectionType = -1;
	}

	for(i=0; i<numLesions; i++)
	{
		found=false;
		for(j=0; j<lesionVoxels[i].numPaint; j++)
		{
			z = lesionVoxels[i].sliceArray[j]-1;
			for(n=0; n<lesionVoxels[i].paint[j].GetSize(); n++)
			{
				x = lesionVoxels[i].paint[j][n].x;
				y = lesionVoxels[i].paint[j][n].y;
				k = z*sizexy+y*sizex+x;
				// check if lesion overlaps any blobs-- may overlap more than one!
				if(blobArray[k]>0) 
				{
					detections[blobArray[k]-1].matchedLesion = lesions[i].lesionKey;	// may match to multiple lesions, only catch the last one 
					detections[blobArray[k]-1].detectionType = lesions[i].lesionType;	// may match to multiple lesions, only catch the last one 
					found = true;
					std::cout<<"\n Lesion: "<<lesions[i].lesionKey<<" found by: "<<blobArray[k]-1;
				}
			}
		}	// for j
		if(found) 
		{
			foundMet++;
			if(lesions[i].lesionSize>=lesionSizeCutoff) foundMet10++;
			if(lesions[i].lesionType==1) 
			{
				foundLytic++;
				if(lesions[i].lesionSize>=lesionSizeCutoff) foundLytic10++;
			}

			lesions[i].lesionStatus = 1;
		}
		else
		{
			lesions[i].lesionStatus = 0;
		}

		totalMet++;
		if(lesions[i].lesionSize>=lesionSizeCutoff) totalMet10++;
		
		if(lesions[i].lesionType==1) 
		{
			totalLytic++;
			if(lesions[i].lesionSize>=lesionSizeCutoff) totalLytic10++;
		}
		lesions[i].hit=false; //initialize so we can keep track of which lesions have already been detected, to avoid double-counting
	}	// for i

	tpDetection=0; tpLyticDetection=0;
	totalDetection=0; fpDetection=0; fpLyticDetection=0;
	totalDetection = numDetections;
	int key = 0;

	for(k=0; k<numDetections; k++) 
	{
		if(detections[k].matchedLesion>0) 
		{
			std::cout<<"\n Lesion: "<<detections[k].matchedLesion<<" detected by: "<<k;
			key = detections[k].matchedLesion;
			for(i=0; i<numLesions; i++)
			{
				if(lesions[i].lesionKey==key) break;
			}
			if(lesions[i].hit==false)
			{
				tpDetection++;
				lesions[i].hit=true;
				//if(detections[k].detectionType==1) tpLyticDetection++;
			}
			else
			{
				fpDetection--;
			}
			//tpDetection++;
		}
	}

	fpDetection = fpDetection+totalDetection-tpDetection;
	fpLyticDetection = totalDetection-tpLyticDetection;

	delete detectionBlobMap;
	return CIS_OK;
}

									   

int NIH_BoneSegmentation_MatchDetections2D(intVec2DynArray &detectionRegion, int detectionSlice,
									   LesionStructure *lesions, int numLesions,
									   PaintingStruct *lesionVoxels,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationParameter &segPara,
									   int &matchedLesion, int &matchedLesionType, float &matchedLesionSize, float &matchedOverlap)
{
	matchedLesion = matchedLesionType = -1;
	matchedOverlap = 0;
	matchedLesionSize = -1;
	
	if(numLesions==0 || detectionRegion.GetSize()==0) return CIS_ERROR;

	CIS_Array_Image2D_short *detectionBlobMap;
	int x, y, z, k;
	int sizex, sizey, sizez, sizexy;
	int i, j, n;

	sizex = maskImg3D->Num_Cols();
	sizey = maskImg3D->Num_Rows();
	sizez = maskImg3D->Num_Levels();

	sizexy = sizex*sizey;

	detectionBlobMap = new CIS_Array_Image2D_short(sizex, sizey);
	short *blobArray = detectionBlobMap->GetArray();

	for(n=0; n<detectionRegion.GetSize(); n++)
	{
		x = detectionRegion[n].x;
		y = detectionRegion[n].y;
		detectionBlobMap->FastSet(x, y, 1);
	}
	
	bool found;
	int overlap, maxOverlap;
	int curMatched, maxMatched;

	maxOverlap = overlap = 0;
	curMatched = maxMatched = -1;
	for(i=0; i<numLesions; i++)
	{
		found=false;
		overlap = 0;
		for(j=0; j<lesionVoxels[i].numPaint; j++)
		{
			z = lesionVoxels[i].sliceArray[j]-1;
			if(z!=detectionSlice) continue;

			for(n=0; n<lesionVoxels[i].paint[j].GetSize(); n++)
			{
				x = lesionVoxels[i].paint[j][n].x;
				y = lesionVoxels[i].paint[j][n].y;
				k = x+y*sizex;

				if(blobArray[k]>0) 
				{
					found = true;
					overlap ++;
					curMatched = i;
				}
			}
		}	// for j


		if(found && overlap>maxOverlap) 
		{
			maxOverlap = overlap;
			maxMatched = curMatched;
		}
	}	// for i

	if(maxMatched>=0)
	{
		matchedLesion = lesions[maxMatched].lesionKey;
		if(lesions[maxMatched].lesionArea>0) matchedOverlap = (float)maxOverlap/(float)lesions[maxMatched].lesionArea;
		else matchedOverlap=0;
		matchedLesionType = lesions[maxMatched].lesionType;
		matchedLesionSize = lesions[maxMatched].lesionSize;
	}

	delete detectionBlobMap;
	return CIS_OK;
}


int NIH_BoneSegmentation_MatchDetections3D(intVec3DynArray &detectionRegion,
										   LesionStructure *lesions, int numLesions,
										   PaintingStruct *lesionVoxels,
										   CIS_Array_Image3D_short *maskImg3D,
										   SpineSegmentationParameter &segPara,
										   int &matchedLesion, int &matchedLesionType, float &matchedLesionSize, float &matchedOverlap)
{
	matchedLesion = matchedLesionType = -1;
	matchedOverlap = 0;
	matchedLesionSize = -1;

	if(numLesions==0 || detectionRegion.GetSize()==0) return CIS_ERROR;

	CIS_Array_Image3D_short *detectionBlobMap;
	int x, y, z, k;
	int sizex, sizey, sizez, sizexy;
	int i, j, n;

	sizex = maskImg3D->Num_Cols();
	sizey = maskImg3D->Num_Rows();
	sizez = maskImg3D->Num_Levels();

	sizexy = sizex*sizey;

	detectionBlobMap = new CIS_Array_Image3D_short(sizex, sizey, sizez);
	short *blobArray = detectionBlobMap->GetArray();

	for(n=0; n<detectionRegion.GetSize(); n++)
	{
		x = detectionRegion[n].x;
		y = detectionRegion[n].y;
		z = detectionRegion[n].z;
		detectionBlobMap->FastSet(x,y,z,1);
	}

	bool found;
	int overlap, maxOverlap;
	int curMatched, maxMatched;

	maxOverlap = overlap = 0;
	curMatched = maxMatched = -1;
	for(i=0; i<numLesions; i++)
	{
		found = false;
		overlap = 0;
		for(j=0; j<lesionVoxels[i].numPaint; j++)
		{
			z = lesionVoxels[i].sliceArray[j]-1;
			for(n=0; n<lesionVoxels[i].paint[j].GetSize(); n++)
			{
				x = lesionVoxels[i].paint[j][n].x;
				y = lesionVoxels[i].paint[j][n].y;
				k = z*sizexy+y*sizex+x;
				if(blobArray[k]>0)
				{
					found = true;
					overlap++;
					curMatched = i;
				}
			}
		} // for j

		if(found && overlap>maxOverlap)
		{
			maxOverlap = overlap;
			maxMatched = curMatched;
		}
	} // for i
	if(maxMatched>=0)
	{
		matchedLesion = lesions[maxMatched].lesionKey;
		if(lesions[maxMatched].lesionArea>0) matchedOverlap = (float)maxOverlap/(float)lesions[maxMatched].lesionArea;
		else matchedOverlap = 0;
		matchedLesionType = lesions[maxMatched].lesionType;
		matchedLesionSize = lesions[maxMatched].lesionSize;
	}
	
	delete detectionBlobMap;
	return CIS_OK;
}

int BernsteinSmoothing(doubleDynArray &xt, doubleDynArray &yt, doubleDynArray &yt2, int bernstein_power)
{
	int pt_num = xt.GetSize();

	if(pt_num==0 || yt.GetSize()!=pt_num || yt2.GetSize()!=pt_num) return -1;

	double *rxt, *ryt, *coef;
	int r_num, i;

	r_num=0;
	rxt = new double[pt_num];
	ryt = new double[pt_num];
	coef = new double[bernstein_power];

	// initialize 
	for(i=0; i<pt_num; i++) yt2[i]=-1;

	for(i=0; i<pt_num; i++)
	{
		if(yt[i]!=-1)
		{
			rxt[r_num] = xt[i];
			ryt[r_num] = yt[i];
			r_num++;
		}
	}

	if(r_num>bernstein_power)
	{
		Fit_Bernstein(bernstein_power, rxt, ryt, r_num, coef);

		double t_len, t, t0, t1, y1;

		t0 = rxt[0];
		t1 = rxt[r_num-1];
		t_len = t1-t0;

		for(i=0; i<pt_num; i++)
		{
			t = (double)(xt[i]-t0)/t_len;
			y1 = Eval_Bernstein(bernstein_power,coef,t);
			yt2[i] = y1;
		}

	}


	delete rxt;
	delete ryt;
	delete coef;
	
	return 0;
}


int PiecewiseBernsteinSmoothing(doubleDynArray &xt, doubleDynArray &yt, doubleDynArray &yt2, int bernstein_power, int piece_size)
{
	int pt_num = xt.GetSize();

	if(pt_num==0 || yt.GetSize()!=pt_num || yt2.GetSize()!=pt_num) return -1;

	double *rxt, *ryt, *coef, *mat;
	int r_num, i;

	r_num=0;
	rxt = new double[pt_num];
	ryt = new double[pt_num];
	mat = new double[pt_num];
	coef = new double[bernstein_power];

	// initialize 
	for(i=0; i<pt_num; i++) yt2[i]=-1;

	for(i=0; i<pt_num; i++)
	{
		if(yt[i]!=-1)
		{
			rxt[r_num] = xt[i];
			ryt[r_num] = yt[i];
			mat[r_num] = i;
			r_num++;
		}
	}

	if(r_num>bernstein_power)
	{
		double t_len, t, t0, t1, y1;

		int *fit_st, *fit_ed, *eval_st, *eval_ed;
		int piece_num, p;

		piece_num = (r_num-1)/piece_size;
		if(piece_num==0) piece_num=1;

		fit_st = new int[piece_num];
		fit_ed = new int[piece_num];
		eval_st = new int[piece_num];
		eval_ed = new int[piece_num];

		// special case for first piece and last piece
		// first piece;
		fit_st[0] = 0;
		eval_st[0] = 0;
		if(r_num>piece_size*2) 
		{
			fit_ed[0]=piece_size*2;
			eval_ed[0]=mat[piece_size];
		}
		else 
		{
			fit_ed[0]=r_num-1;
			eval_ed[0]=pt_num-1;
		}

		for(p=1; p<piece_num; p++)
		{
			fit_st[p] = (p-1)*piece_size;
			eval_st[p] = eval_ed[p-1];

			if((p+2)*piece_size<r_num)
			{
				fit_ed[p] = (p+2)*piece_size;
				eval_ed[p] = mat[(p+1)*piece_size];
			}
			else
			{
				fit_ed[p] = r_num-1;
				eval_ed[p] = pt_num-1;
			}
		}


		for(p=0; p<piece_num; p++)
		{
			Fit_Bernstein(bernstein_power, &rxt[fit_st[p]], &ryt[fit_st[p]], fit_ed[p]-fit_st[p]+1, coef);

			t0 = rxt[fit_st[p]];
			t1 = rxt[fit_ed[p]];
			t_len = t1-t0;

			for(i=eval_st[p]; i<=eval_ed[p]; i++)
			{
				t = (double)(xt[i]-t0)/t_len;
				if(t>1.1 || t<-0.1) continue;		// avoid excessive extrapolation
				if(t==1) t=0.99;
				y1 = Eval_Bernstein(bernstein_power,coef,t);
				yt2[i] = y1;
			}
		}

		delete fit_st;
		delete fit_ed;
		delete eval_st;
		delete eval_ed;
	}


	delete rxt;
	delete ryt;
	delete mat;
	delete coef;
	
	return 0;
}



int BSplineSmoothing(doubleDynArray &xt, doubleDynArray &yt, doubleDynArray &yt2, int bspline_degree)
{
	int pt_num = xt.GetSize();

	if(pt_num==0 || yt.GetSize()!=pt_num) return -1;

	vec2DynArray new_list;
	int r_num, i;

	r_num=0;

	for(i=0; i<pt_num; i++)
	{
		if(yt[i]!=-1)
		{
			new_list.Add(Vec2(xt[i], yt[i]));
			r_num++;
		}
	}

	if(r_num>bspline_degree)
	{
		CIS_BSpline_Curve_2D curve(bspline_degree, r_num, new_list);
	
		// re-sample the curve
		double length= curve.Length();
		double s, s_step=length/(double)(pt_num*2);
		Vec2 pt;
	

		for(i=0; i<pt_num; i++) yt2[i]=-1;
		
		for(s=0, i=0; i<pt_num && s<=length; s+=s_step)
		{
			pt = curve.Position(s);

			if(fabs(pt.x-xt[i])<=0.5 || pt.x>xt[i]) 
			{
				yt2[i]=pt.y;
				i++;
			}
		}

		int closest;
		for(i=0; i<pt_num; i++) 
		{
			if(yt2[i]==-1)		// closet neighbor to interpolate
			{
				closest=i-1;
				while(closest>=0 && yt2[closest]==-1) closest--;

				if(closest>=0 && yt2[closest]!=-1) yt2[i]=yt2[closest];
				else
				{
					closest=i+1;
					while(closest<pt_num && yt2[closest]==-1) closest++;
					if(closest<pt_num && yt2[closest]!=-1) yt2[i]=yt2[closest];
				}				
			}
		}

	}	// if r_num


	return 0;
}

double PedicleTemplateMatching(CIS_Array_Image3D_short *img3D, CIS_Array_Image3D_short *maskImg3D, int slice, Vec2 candidate,
							   SpineSegmentationParameter &segPara, SpineSegmentationInfo &segInfo, SpineSegmentationMaskStructure &maskStruct,
							   VertebraStruct2D *vertebraTemplate, int side, bool filled	// side: -1 left, 1 right
							   )
{
	int i, k, k2, x, y, z, sizex, sizey, sizez, sizexy;
	short *imgA, *maskA;
	double score = 0;

	z = slice;
	imgA = img3D->GetArray();
	maskA = maskImg3D->GetArray();
	sizex = img3D->Num_Cols();
	sizey = img3D->Num_Rows();
	sizez = img3D->Num_Levels();

	sizexy = sizex*sizey;

	k = z*sizexy;
	Vec2 tip1, tip2, medial[20], org_medial[20], dir, dirx, curPos, nextPos, maxPos;
	double steps;
	int r, grad, maxGrad, k3;
	double max_ru[20], max_rd[20], medial_shift[20], medial_shift_interpolated[20];

	tip1 = candidate;
	tip2.x = segInfo.cordCenter[z].x;
	tip2.y = segInfo.cordCenter[z].y;

	dir = (tip2-tip1).normalize();
	tip2 = tip2-dir*segInfo.cordRadius[z].x;

	dirx = dir.perp();
	if(dirx.y<0)
	{
		dirx.x = -dirx.x;
		dirx.y = -dirx.y;
	}

	double pedicleLength;
	pedicleLength = (tip2-tip1).len();
	steps = pedicleLength/(pedicleCount-1);

	medial[0] = tip1;
	medial[pedicleCount-1] = tip2;

	for(i=1; i<pedicleCount-1; i++) medial[i] = medial[i-1]+dir*steps;

	for(i=0; i<pedicleCount; i++) medial_shift[i] = medial_shift_interpolated[i]=0;

	for(i=0; i<pedicleCount; i++) org_medial[i] = medial[i];

	// two iteration to smoothout the medial axis

	for(int it=0; it<2; it++)
	{
		// search up and down from medial to locate pedicle border
		//
		for(i=0; i<pedicleCount; i++)
		{
			max_ru[i] = max_rd[i] = 0;
	
			// up
			curPos = medial[i]-dirx;
			maxGrad = 0;
			for(r=1; r<segInfo.cordRadius[z].y*1.5; r++)
			{
				nextPos = curPos-dirx;
				k2 = k + (int)(curPos.y+0.5)*sizex + (int)(curPos.x+0.5);
				k3 = k + (int)(nextPos.y+0.5)*sizex + (int)(nextPos.x+0.5);

				if(maskA[k2]>100) break;		// exit if already touch other structure

				if(imgA[k2]>segPara.boneThresh && imgA[k3]<segPara.boneThresh && 
					imgA[k3+side*segInfo.cordRadius[z].x]>700 && imgA[k3-segInfo.cordRadius[z].x*sizex]>700)
				{
					grad = imgA[k2]-imgA[k3];
					if(grad>maxGrad && grad>50)
					{
						maxGrad=grad;
						maxPos = curPos;
						max_ru[i] = r;
					}
				}
				curPos = nextPos;
			}	// for r

			// find the intensity match if can't get the gradient match
			if(max_ru[i]==0)
			{
				curPos = medial[i];
				for(r=0; r<segInfo.cordRadius[z].y*1.5; r++)
				{
					k2 = k + (int)(curPos.y+0.5)*sizex + (int)(curPos.x+0.5);
					if(imgA[k2]<segPara.boneThresh || maskA[k2]>100) break;
					curPos -= dirx;
				}
				max_ru[i]=r;
			}

			// down
			curPos = medial[i]+dirx;
			maxGrad = 0;
			for(r=1; r<segInfo.cordRadius[z].y*1.5; r++)
			{
				nextPos = curPos+dirx;
				k2 = k + (int)(curPos.y+0.5)*sizex + (int)(curPos.x+0.5);
				k3 = k + (int)(nextPos.y+0.5)*sizex + (int)(nextPos.x+0.5);
	
				if(maskA[k2]>100) break;		// exit if already touch other structure

				if(imgA[k2]>segPara.boneThresh && imgA[k3]<segPara.boneThresh && 
					imgA[k3+side*segInfo.cordRadius[z].x]>700 && imgA[k3-segInfo.cordRadius[z].x*sizex]>700)
				{
					grad = imgA[k2]-imgA[k3];
					if(grad>maxGrad && grad>50)
					{
						maxGrad=grad;
						maxPos = curPos;
						max_rd[i] = r;
					}
				}
				curPos = nextPos;
			}	// for r

			// find the intensity match if can't get the gradient match
			if(max_rd[i]==0)
			{
				curPos = medial[i];
				for(r=0; r<segInfo.cordRadius[z].y*1.5; r++)
				{
					k2 = k + (int)(curPos.y+0.5)*sizex + (int)(curPos.x+0.5);
					if(imgA[k2]<segPara.boneThresh || maskA[k2]>100) break;
					curPos += dirx;
				}
				max_rd[i]=r;
			}
		}	// for i

		// adjust medial axis
		for(i=1; i<pedicleCount-1; i++)
		{
			if(max_rd[i]>0 && max_ru[i]==0)
			{
				medial[i] += dirx*(max_rd[i]/2);
				max_rd[i] /=2;
				medial_shift[i] += max_rd[i]/2;
///				medial_shift_interpolated[i] = medial_shift[i] + max_rd[i]/2;
			}
			else if(max_rd[i]==0 && max_ru[i]>0)
			{
				medial[i] -= dirx*(max_ru[i]/2);
				max_ru[i] /=2;
				medial_shift[i] -= max_ru[i]/2;
///				medial_shift_interpolated[i] = medial_shift[i] - max_ru[i]/2;
			}
			else
			{
				medial[i] += dirx*((max_rd[i]-max_ru[i])/2);
				medial_shift[i] += (max_rd[i]-max_ru[i])/2;
				max_rd[i] -= (max_rd[i]-max_ru[i])/2;
				max_ru[i] += (max_rd[i]-max_ru[i])/2;
///				medial_shift_interpolated[i] = medial_shift[i] + (max_rd[i]-max_ru[i])/2;
			}
		}

/*		doubleDynArray xt, yt, yt2;
		int bernstein_power=2;

		xt.SetSize(pedicleCount);
		yt.SetSize(pedicleCount);
		yt2.SetSize(pedicleCount);

		for(i=0; i<pedicleCount; i++)
		{
			xt[i] = i;
			yt[i] = yt2[i] = -1;
		}

		// smooth the medial axis
		for(i=0; i<pedicleCount; i++) yt[i]=medial_shift_interpolated[i];
		BernsteinSmoothing(xt, yt, yt2, bernstein_power);
		for(i=1; i<pedicleCount-1; i++) 
		{
			medial_shift_interpolated[i]=yt2[i];

			medial[i] = org_medial[i]+dirx*medial_shift_interpolated[i];

			if(max_rd[i]>0) max_rd[i] += (medial_shift_interpolated[i]-medial_shift[i]);
			if(max_ru[i]>0) max_ru[i] -= (medial_shift_interpolated[i]-medial_shift[i]);
			medial_shift[i] = medial_shift_interpolated[i];
		}
*/	
	}

	// compute the score
	
/*	// use the weighed perimeter as the score
	score = 0;
	score += pedicleLength*2;	// put more weights on the pedicle length
	for(i=0; i<pedicleCount*1/3; i++)
	{
		score += max_ru[i]*1.5;
		score += max_rd[i]*1.5;
	}
	for(i=pedicleCount*1/3; i<pedicleCount*2/3; i++)
	{
		score += max_ru[i];
		score += max_rd[i];
	}
*/
	// use the covered area as the score
	double last_ru=0, last_rd=0;
	score = 0;
	k = z*sizexy;
	for(i=0; i<pedicleCount; i++)
	{
		score += max_ru[i]*steps;
		score += max_rd[i]*steps;

		// penalize the gap
		if(max_ru[i]==0 && max_rd[i]==0) score -= (last_ru+last_rd)/2*steps;
		else if(max_ru[i]>0 && max_rd[i]==0) score -= max_ru[i]/2*steps;
		else if(max_ru[i]==0 && max_rd[i]>0) score -= max_rd[i]/2*steps;

		k2 = k+medial[i].y*sizex+medial[i].x;
		if(imgA[k2]<segPara.boneThresh-100) score -= (max_ru[i]+max_rd[i])*steps;

		if(max_ru[i]>0) last_ru=max_ru[i];
		if(max_rd[i]>0) last_rd=max_rd[i];
	}

/*	// count located border
	int borderCount = 0, mid_borderCount=0, tip_borderCount=0;
	for(i=0; i<pedicleCount; i++)
	{
		if(max_rd[i]>0) borderCount++;
		if(max_ru[i]>0) borderCount++;
	}
	for(i=0; i<pedicleCount/2; i++)
	{
		if(max_rd[i]>0) tip_borderCount++;
		if(max_ru[i]>0) tip_borderCount++;
	}
	// just count the mid section
	for(i=pedicleCount*1/3; i<pedicleCount*2/3; i++)
	{
		if(max_rd[i]>0) mid_borderCount++;
		if(max_ru[i]>0) mid_borderCount++;
	}

	if(borderCount<pedicleCount) score = 0;	// forced 0
	if(mid_borderCount<pedicleCount/3) score = 0;	// forced 0
	if(tip_borderCount<pedicleCount/2) score = 0;	// forced 0
*/

	// fill the data structure
	if(filled)
	{
		doubleDynArray xt, yt, yt2;
		int bernstein_power=2;

		xt.SetSize(pedicleCount);
		yt.SetSize(pedicleCount);
		yt2.SetSize(pedicleCount);

		for(i=0; i<pedicleCount; i++)
		{
			xt[i] = i;
			yt[i] = yt2[i] = -1;
		}

		// smooth the medial axis
		for(i=0; i<pedicleCount; i++) yt[i]=medial_shift[i];
		BernsteinSmoothing(xt, yt, yt2, bernstein_power);
		for(i=0; i<pedicleCount; i++) 
		{
			medial_shift_interpolated[i]=yt2[i];

//			medial[i] += dirx*(medial_shift_interpolated[i]-medial_shift[i]);
			medial[i] = org_medial[i]+dirx*medial_shift_interpolated[i];

			if(max_rd[i]>0) max_rd[i] += medial_shift[i];
			if(max_ru[i]>0) max_ru[i] -= medial_shift[i];
//			if(max_rd[i]>0) max_rd[i] -= medial_shift_interpolated[i]-medial_shift[i];
//			if(max_ru[i]>0) max_ru[i] += medial_shift_interpolated[i]-medial_shift[i];
		}
		
		// interpolate the border

		bernstein_power=3;
		// down
		for(i=0; i<pedicleCount; i++)
		{
			if(max_rd[i]>0) yt[i] = max_rd[i];
			else yt[i]=-1;
		}
		// curb the end point
		if(yt[0]==-1) yt[0] = 2;
		if(yt[pedicleCount-1]==-1) yt[pedicleCount-1] = segInfo.cordRadius[z].y;
		BernsteinSmoothing(xt, yt, yt2, bernstein_power);
		for(i=0; i<pedicleCount; i++) 
		{
			if(max_rd[i]<=0 && yt2[i]>segInfo.cordRadius[z].y*2) max_rd[i]=segInfo.cordRadius[z].y*2;	// curb the interpolation
			else max_rd[i] = yt2[i];
		}

		// up
		for(i=0; i<pedicleCount; i++)
		{
			if(max_ru[i]>0) yt[i] = max_ru[i];
			else yt[i]=-1;
		}
		if(yt[0]==-1) yt[0] = 2;
		if(yt[pedicleCount-1]==-1) yt[pedicleCount-1] = segInfo.cordRadius[z].y;
		BernsteinSmoothing(xt, yt, yt2, bernstein_power);
		for(i=0; i<pedicleCount; i++) 
		{
			if(max_ru[i]<=0 && yt2[i]>segInfo.cordRadius[z].y*2) max_ru[i]=segInfo.cordRadius[z].y*2;
			else max_ru[i] = yt2[i];
		}

		// need some smoothing before the output
		//

		// local refinement
		// up
		for(i=0; i<pedicleCount; i++)
		{
			curPos = org_medial[i]-dirx*max_ru[i];
			for(r=0;r<segInfo.cordRadius[z].y; r++)
			{
				k2 = k+curPos.y*sizex+curPos.x;
				if(imgA[k2]<segPara.boneThresh || maskA[k2]>100) break;

				max_ru[i] ++;
				curPos -= dirx;
			}
		}

		// down
		for(i=0; i<pedicleCount; i++)
		{
			curPos = org_medial[i]+dirx*max_rd[i];
			for(r=0;r<segInfo.cordRadius[z].y; r++)
			{
				k2 = k+curPos.y*sizex+curPos.x;
				if(imgA[k2]<segPara.boneThresh || maskA[k2]>100) break;

				max_rd[i] ++;
				curPos += dirx;
			}
		}

		// smooth it a second time
		// up
		for(i=0; i<pedicleCount; i++)
		{
			if(max_ru[i]>0) yt[i] = max_ru[i];
			else yt[i]=-1;
		}
		if(yt[0]==-1) yt[0] = 2;
		if(yt[pedicleCount-1]==-1) yt[pedicleCount-1] = segInfo.cordRadius[z].y;
		BernsteinSmoothing(xt, yt, yt2, bernstein_power);
		for(i=0; i<pedicleCount; i++) 
		{
			if(yt2[i]>segInfo.cordRadius[z].y*2) max_ru[i]=segInfo.cordRadius[z].y*2;
			else max_ru[i] = yt2[i];
		}
		// down
		for(i=0; i<pedicleCount; i++)
		{
			if(max_rd[i]>0) yt[i] = max_rd[i];
			else yt[i]=-1;
		}
		// curb the end point
		if(yt[0]==-1) yt[0] = 2;
		if(yt[pedicleCount-1]==-1) yt[pedicleCount-1] = segInfo.cordRadius[z].y;
		BernsteinSmoothing(xt, yt, yt2, bernstein_power);
		for(i=0; i<pedicleCount; i++) 
		{
			if(yt2[i]>segInfo.cordRadius[z].y*2) max_rd[i]=segInfo.cordRadius[z].y*2;	// curb the interpolation
			else max_rd[i] = yt2[i];
		}

		if(side==-1)
		{
			for(i=0; i<pedicleCount; i++)
			{
///				vertebraTemplate[z].leftPedicleMedial[i] = medial[i];
///				if(max_ru[i]>0) vertebraTemplate[z].leftPedicleUp[i] = medial[i]-dirx*max_ru[i];
///				else vertebraTemplate[z].leftPedicleUp[i] = Vec2(-1,-1);

				vertebraTemplate[z].leftPedicleMedial[i] = medial[i];
				if(max_ru[i]>0) vertebraTemplate[z].leftPedicleUp[i] = org_medial[i]-dirx*max_ru[i];
				else vertebraTemplate[z].leftPedicleUp[i] = Vec2(-1,-1);

				if(max_rd[i]>0) vertebraTemplate[z].leftPedicleDown[i] = org_medial[i]+dirx*max_rd[i];
				else vertebraTemplate[z].leftPedicleDown[i] = Vec2(-1,-1);

				// visualize
				k2 = k + (int)(vertebraTemplate[z].leftPedicleMedial[i].y+0.5)*sizex + (int)(vertebraTemplate[z].leftPedicleMedial[i].x+0.5);
				maskA[k2] = maskStruct.mask_spinalCord;
				if(max_ru[i]>0)
				{
					k2 = k + (int)(vertebraTemplate[z].leftPedicleUp[i].y+0.5)*sizex + (int)(vertebraTemplate[z].leftPedicleUp[i].x+0.5);
					maskA[k2] = maskStruct.mask_falseDetection;
				}
				if(max_rd[i]>0)
				{
					k2 = k + (int)(vertebraTemplate[z].leftPedicleDown[i].y+0.5)*sizex + (int)(vertebraTemplate[z].leftPedicleDown[i].x+0.5);
					maskA[k2] = maskStruct.mask_falseDetection;
				}

			}
		}
		else if(side==1)
		{
			for(i=0; i<pedicleCount; i++)
			{
///				vertebraTemplate[z].rightPedicleMedial[i] = medial[i];
///				if(max_ru[i]>0) vertebraTemplate[z].rightPedicleUp[i] = medial[i]-dirx*max_ru[i];
///				if(max_rd[i]>0) vertebraTemplate[z].rightPedicleDown[i] = medial[i]+dirx*max_rd[i];

				vertebraTemplate[z].rightPedicleMedial[i] = medial[i];
				if(max_ru[i]>0) vertebraTemplate[z].rightPedicleUp[i] = org_medial[i]-dirx*max_ru[i];
				if(max_rd[i]>0) vertebraTemplate[z].rightPedicleDown[i] = org_medial[i]+dirx*max_rd[i];

				// visualize
				k2 = k + (int)(vertebraTemplate[z].rightPedicleMedial[i].y+0.5)*sizex + (int)(vertebraTemplate[z].rightPedicleMedial[i].x+0.5);
				maskA[k2] = maskStruct.mask_spinalCord;
				if(max_ru[i]>0)
				{
					k2 = k + (int)(vertebraTemplate[z].rightPedicleUp[i].y+0.5)*sizex + (int)(vertebraTemplate[z].rightPedicleUp[i].x+0.5);
					maskA[k2] = maskStruct.mask_falseDetection;
				}
				if(max_rd[i]>0)
				{
					k2 = k + (int)(vertebraTemplate[z].rightPedicleDown[i].y+0.5)*sizex + (int)(vertebraTemplate[z].rightPedicleDown[i].x+0.5);
					maskA[k2] = maskStruct.mask_falseDetection;
				}
			}
		}
	}

	return score;
}

// reformat the data along the spinal cord in sagittal and coronal direction
//
int NIH_SpineSegmentation_CurvedReformation(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   VertebraStruct2D *vertebraTemplate, 
									   vec3DynArray &spinalCord3D,
									   CIS_Array_Image2D_short *img2D_sag,
									   CIS_Array_Image2D_short *img2D_cor,
									   CIS_Array_Image2D_short *img2D_sag_mask,
									   CIS_Array_Image2D_short *img2D_cor_mask,
									   bool debugMode)
{
	int sizex, sizey, sizez, sizexy;
	short *maskA, *imgA, *sagA, *corA, *sagMaskA, *corMaskA;
	int x, y, z, k, k2;
	float px, py, pz;

	sizex = maskImg3D->Num_Cols();
	sizey = maskImg3D->Num_Rows();
	sizez = maskImg3D->Num_Levels();
	px = img3D->Get_Pixel_SizeX();
	py = img3D->Get_Pixel_SizeY();
	pz = img3D->Get_Pixel_SizeZ();
	sizexy = sizex*sizey;

	maskA = maskImg3D->GetArray();
	imgA = img3D->GetArray();

	Vec3 center3D;
	IntVec2 center;
	int count;
	
	spinalCord3D.SetSize(0);

	// use the center of mass of the spinal cord
	for(z=0; z<sizez; z++)
	{
		k= z*sizex*sizey;
		count =0;
		center = IntVec2(0,0);

		// locate the center
		for(k2=0, y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k2++)
			{
				if(maskA[k+k2]==maskStruct.mask_spinalCord)
				{
					count ++;
					center += IntVec2(x, y);
				}
			}
		}

		if(count>0)
		{
			center.x /= count;
			center.y /= count;

			center3D.x = center.x*px;
			center3D.y = center.y*py;
			center3D.z = img3D->GetSlicePosition(z);

			spinalCord3D.Add(center3D);
		}
	}

	// smooth the curve
	//
	doubleDynArray xt, yt, yt2;
	int bernstein_power = 5, piece_size;

	float piece_length=25;
	piece_size = (int)(piece_length/pz+0.5);

	xt.SetSize(sizez);
	yt.SetSize(sizez);
	yt2.SetSize(sizez);
	for(z=0; z<sizez; z++)
	{
		xt[z] = z;
		yt[z] = -1;
	}

	// x coordinate of cord
	for(z=0; z<spinalCord3D.GetSize(); z++)
	{
		yt[z] = spinalCord3D[z].x;
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=0; z<spinalCord3D.GetSize(); z++)
	{
		 spinalCord3D[z].x = yt2[z];
	}

	// y coordinate of cord
	for(z=0; z<spinalCord3D.GetSize(); z++)
	{
		yt[z] = spinalCord3D[z].y;
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=0; z<spinalCord3D.GetSize(); z++)
	{
		 spinalCord3D[z].y = yt2[z];
	}


	short gray;
	intDynArray filledZ;

	filledZ.SetSize(sizez);
	for(z=0; z<sizez; z++) filledZ[z]=-1;

	img2D_sag->SetSize(sizey, sizez);
	img2D_cor->SetSize(sizex, sizez);
	img2D_sag_mask->SetSize(sizey, sizez);
	img2D_cor_mask->SetSize(sizex, sizez);
	sagA = img2D_sag->GetArray();
	sagMaskA = img2D_sag_mask->GetArray();
	corA = img2D_cor->GetArray();
	corMaskA = img2D_cor_mask->GetArray();

	int k3;
	for(k=0; k<spinalCord3D.GetSize(); k++)
	{
		// sagittal reformation
		x = (int)(spinalCord3D[k].x/px);
		z = img3D->GetSliceLevel(spinalCord3D[k].z,1);

		// filledZ array is used to handle duplicate slices
		//
		if(filledZ[z]==1 && filledZ[z-1]==0) z=z-1;
		while(z<sizez-1 && filledZ[z]==1) z++;
		filledZ[z] = k;
		
		// saggital reformation
		k2 = z*sizexy+x;
		k3 = z*sizey;
		for(y=0; y<img2D_sag->Num_Cols(); y++, k2+=sizex, k3++)
		{
			gray = imgA[k2];
			sagA[k3] = gray;
			gray = maskA[k2];
			sagMaskA[k3] = gray;
		}

		// coronal reformation
		y = (int)(spinalCord3D[k].y/py);

		k2 = z*sizexy+y*sizex;
		k3 = z*sizey;
		for(x=0; x<img2D_cor->Num_Cols(); x++, k2++, k3++)
		{
			gray = imgA[k2];
			corA[k3] = gray;
			gray = maskA[k2];
			corMaskA[k3] = gray;
		}
	}

	// interpolate the missing lines
	// find the top and bottom
	int tz=0;
	while(filledZ[tz]==-1) tz++;
	int bz=sizez-1;
	while(filledZ[bz]==-1) bz--;
	for(z=tz; z<=bz; z++)
	{
		if(filledZ[z]!=-1) continue;

		// find the previous and next line
		int pz, nz;
		pz=z-1;
		while(pz>=0 && filledZ[pz]==-1) pz--;
		nz=z+1;
		while(nz<sizez && filledZ[nz]==-1) nz++;
		
		// sagital image
		short gray, count;
		for(y=0; y<img2D_sag->Num_Cols(); y++)
		{
			gray=0;
			count=0;
			if(pz>=0)
			{
				count++;
				gray+=img2D_sag->FastGet(y, pz);
			}
			if(nz<sizez)
			{
				count++;
				gray+=img2D_sag->FastGet(y, nz);
			}
			if(count>0) gray/=count;
			img2D_sag->FastSet(y, z,gray); 
		}

		// coronal image
		for(y=0; y<img2D_cor->Num_Cols(); y++)
		{
			gray=0;
			count=0;
			if(pz>=0)
			{
				count++;
				gray+=img2D_cor->FastGet(y, pz);
			}
			if(nz<sizez)
			{
				count++;
				gray+=img2D_cor->FastGet(y, nz);
			}
			if(count>0) gray/=count;
			img2D_cor->FastSet(y, z,gray); 
		}
	}

	// interpolate the upper border of spinal cord to make it smooth
	for(z=0; z<sizez; z++)
	{
		xt[z] = z;
		yt[z] = -1;
	}

	// x coordinate of cord
	for(k=0; k<spinalCord3D.GetSize(); k++)
	{
		// sagittal reformation
		y = (int)(spinalCord3D[k].y/py);
		z = img3D->GetSliceLevel(spinalCord3D[k].z,1);
		k3 = z*sizey+y;
		for(;y>0;y--,k3--)
		{
			if(sagMaskA[k3]!=maskStruct.mask_spinalCord)
			{
				yt[z] = y;
				break;
			}
		}
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	
	for(z=tz; z<=bz; z++)
	{
		if(filledZ[z]==-1 || yt2[z]==-1) continue;

		x = (int)(spinalCord3D[filledZ[z]].x/px);

		// find the previous and next line
		for(y=yt2[z], k3 = z*sizey+y; y>0; y--, k3--)
		{
			if(sagMaskA[k3]==maskStruct.mask_spinalCord)
			{
				sagMaskA[k3]=maskStruct.mask_spongyBone;
				maskA[x+y*sizex+z*sizexy]=maskStruct.mask_spongyBone;
			}
			else break;
		}
	}

	return CIS_OK;
}



int NIH_SpineSegmentation_SpinePartition(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   VertebraStruct2D *vertebraTemplate, 
									   vec3DynArray &spinalCord3D,
									   CIS_Array_Image2D_short *img2D_sag,
									   CIS_Array_Image2D_short *img2D_cor,
									   CIS_Array_Image2D_short *img2D_sag_mask,
									   CIS_Array_Image2D_short *img2D_cor_mask,
									   vec3DynArray &spineCenter, vec3DynArray &spineNormalX, vec3DynArray &spineNormalY, vec3DynArray &spineNormalBack, 
									   doubleDynArray &spineWidthLeft, doubleDynArray &spineWidthRight, doubleDynArray &spineWidthUp, doubleDynArray &spineWidthDown, 
									   intDynArray &pedicleValley, intDynArray &backValley, 
									   bool debugMode)
{
	if(spinalCord3D.GetSize()<10) return CIS_ERROR;
	int sizex, sizey, sizez, sizexy;
	short *maskA, *imgA;
	int x, y, z, k;
	float px, py, pz;

	sizex = maskImg3D->Num_Cols();
	sizey = maskImg3D->Num_Rows();
	sizez = maskImg3D->Num_Levels();
	px = img3D->Get_Pixel_SizeX();
	py = img3D->Get_Pixel_SizeY();
	pz = img3D->Get_Pixel_SizeZ();
	sizexy = sizex*sizey;

	maskA = maskImg3D->GetArray();
	imgA = img3D->GetArray();

	float vsizex = img3D->Get_Volume_SizeX();
	float vsizey = img3D->Get_Volume_SizeY();
	float vsizez = img3D->Get_Volume_SizeZ();

	// compute the projected curve on sagittal plane (x plane)
	intDynArray leftWidth, rightWidth, upWidth, downWidth;
	intDynArray filledZ;
	vec2DynArray projSag, projCor;
	vec2DynArray normalSag, normalCor;

	filledZ.SetSize(sizez);

	int sx, sy;
	short mask;

	upWidth.SetSize(sizez);
	downWidth.SetSize(sizez);
	leftWidth.SetSize(sizez);
	rightWidth.SetSize(sizez);
		
	spineCenter.SetSize(sizez);
	spineNormalX.SetSize(sizez);
	spineNormalY.SetSize(sizez);
	spineWidthUp.SetSize(sizez);
	spineWidthDown.SetSize(sizez);
	spineWidthLeft.SetSize(sizez);
	spineWidthRight.SetSize(sizez);

	for(z=0; z<sizez; z++)
	{
		spineCenter[z] = Vec3(0,0,0);
		spineNormalX[z] = Vec3(0,0,0);
		spineNormalY[z] = Vec3(0,0,0);
		spineWidthUp[z] = -1;
		spineWidthDown[z] = -1;
		spineWidthLeft[z] = -1;
		spineWidthRight[z] = -1;
	}

	for(z=0; z<sizez; z++) 
	{
		filledZ[z]=-1;
		rightWidth[z] = leftWidth[z]=-1;
		upWidth[z] = downWidth[z]=-1;
	}
		
	projSag.SetSize(spinalCord3D.GetSize());
	projCor.SetSize(spinalCord3D.GetSize());

	for(k=0; k<spinalCord3D.GetSize(); k++)
	{
		z = img3D->GetSliceLevel(spinalCord3D[k].z,1);

		projSag[k].x = spinalCord3D[k].y;
		projSag[k].y = z*pz;
		projCor[k].x = spinalCord3D[k].x;
		projCor[k].y = z*pz;

		if(filledZ[z]==1 && filledZ[z-1]==-1) z=z-1;
		while(z<sizez-1 && filledZ[z]!=-1) z++;
		filledZ[z] = k;

		spineCenter[z].x = spinalCord3D[k].x;
		spineCenter[z].y = spinalCord3D[k].y;
		spineCenter[z].z = img3D->GetSlicePosition(z);

		// set up the spinal column
		x = (int)(spinalCord3D[k].x/px+0.5);
		y = (int)(spinalCord3D[k].y/py+0.5);
		// to the right
		int count = 0;
		for(sx=x; sx<sizex; sx++)
		{
			mask = maskImg3D->FastGet(sx, y, z);
			if(mask!=maskStruct.mask_spongyBone && mask!=maskStruct.mask_corticalBone && mask!=maskStruct.mask_spinalCord)
			{
				count++;
			}
			else count=0;

			if(count==3)	// search for the border, break when meet 3 non-vertebra pixels
			{
				rightWidth[z] = sx-x;
				break;
			}
		}

		// to the left
		count =0;
		for(sx=x; sx>0; sx--)
		{
			mask = maskImg3D->FastGet(sx, y, z);
			if(mask!=maskStruct.mask_spongyBone && mask!=maskStruct.mask_corticalBone && mask!=maskStruct.mask_spinalCord)
			{
				count++;
			}
			else count=0;

			if(count==3)	// search for the border, break when meet 3 non-vertebra pixels
			{
				leftWidth[z] = x-sx;
				break;
			}
		}

		// set up the spinal column
		x = (int)(spinalCord3D[k].x/px+0.5);
		y = (int)(spinalCord3D[k].y/py+0.5);
		// upward
		for(sy=y; sy>0; sy--)
		{
			mask = maskImg3D->FastGet(x, sy, z);
			if(mask==maskStruct.mask_body)		
			{
				upWidth[z] = y-sy;	// stop at the first non-vertebra pixel
				break;
			}
		}

		// downward
		for(sy=y; sy<=y+upWidth[z]*2; sy++)
		{
			for(sx=x-leftWidth[z]; sx<=x+rightWidth[z]; sx++)
			{
				mask = maskImg3D->FastGet(sx, sy, z);
				if(mask==maskStruct.mask_spongyBone || mask==maskStruct.mask_corticalBone || mask==maskStruct.mask_spinalCord)
				{
					downWidth[z] = sy-y;	// stop at the last vertebra pixel
				}
			}
		}
	}	// for k

	// smooth upWidth and downWidth using piecewise b-spine
	//
	doubleDynArray xt, yt, yt2;
	int bernstein_power = 5, piece_size;
	float piece_length=25;
	piece_size = (int)(piece_length/pz+0.5);

	xt.SetSize(sizez);
	yt.SetSize(sizez);
	yt2.SetSize(sizez);
	for(z=0; z<sizez; z++)
	{
		xt[z] = z;
		yt[z] = -1;
	}
	// upWidth
	for(z=0; z<sizez; z++)
	{
		yt[z] = upWidth[z];
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=0; z<sizez; z++)
	{
		if(upWidth[z]!=-1) upWidth[z] = yt2[z];
	}

	// DownWidth
	for(z=0; z<sizez; z++)
	{
		yt[z] = downWidth[z];
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=0; z<sizez; z++)
	{
		if(downWidth[z]!=-1) downWidth[z] = yt2[z];
	}

	// smooth upWidth and downWidth
	// laplacian smoothing
	intDynArray tmpWidth;
	tmpWidth = upWidth;
	for(z=1; z<sizez-1; z++)
	{
		if(tmpWidth[z-1]==0 || tmpWidth[z+1]==0) continue;
		if(tmpWidth[z]<tmpWidth[z-1] && tmpWidth[z]<tmpWidth[z+1]) upWidth[z] = (tmpWidth[z-1]+tmpWidth[z+1])/2;
		else upWidth[z] = (tmpWidth[z-1]+tmpWidth[z]+tmpWidth[z+1])/3;
	}
	tmpWidth = downWidth;
	for(z=1; z<sizez-1; z++)
	{
		if(tmpWidth[z-1]==0 || tmpWidth[z+1]==0) continue;
		downWidth[z] = (tmpWidth[z-1]+tmpWidth[z]+tmpWidth[z+1])/3;
	}

	// assign globle variable, return with the function
	for(z=0; z<sizez; z++)
	{
		spineWidthUp[z] = upWidth[z]*py;
		spineWidthDown[z] = downWidth[z]*py;
	}

	// compute the left and right width
	int avgLeftWidth, avgRightWidth;
	int count, count1, count2;;
	avgLeftWidth = avgRightWidth = 0;
	count = 0;
	for(z=0; z<sizez; z++)
	{
		if(leftWidth[z]>0 && rightWidth[z]>0)
		{
			avgLeftWidth += leftWidth[z];
			avgRightWidth += rightWidth[z]; 
			count ++;
		}
	}
	 	
	if(count!=0)
	{
		avgLeftWidth /= count;
		avgRightWidth /= count;

		avgLeftWidth += 4;
		avgRightWidth += 4;
	}
		
	// assign globle variables
	for(z=0; z<sizez; z++)
	{
		spineWidthLeft[z] = avgLeftWidth*px;
		spineWidthRight[z] = avgRightWidth*px;
//		spineWidthLeft[z] = leftWidth[z]*px;
//		spineWidthRight[z] = rightWidth[z]*px;
	}

	// compute the normal in sag and cor plane separately
	//
	CIS_2D_Model_Curve *projectedCord;
	projectedCord = new CIS_2D_Model_Curve();

	normalSag.SetSize(sizez);
	normalCor.SetSize(sizez);
	for(z=0; z<sizez; z++)
	{
		normalSag[z] = Vec2(0,0);
		normalCor[z] = Vec2(0,0);
	}

	// sagittal plane
	projectedCord->UpdateControlList(projSag);

	// compute the normal from the b-spline curve
	double curveLength=projectedCord->Length();
	double stepSize=curveLength/projSag.GetSize();
	double s=0;
	for(z=0;z<sizez;z++)
	{
		k = filledZ[z];
		if(k<0) continue;
		normalSag[z] = projectedCord->Normal(s);
		s += stepSize;
	}

	// smooth the normal using Laplacian method
/*	vec2DynArray tmp_normalA;
	tmp_normalA = normalSag;
	for(z=1; z<sizez-1; z++)
	{
		if((tmp_normalA[z-1].x==0 && tmp_normalA[z-1].y==0) || (tmp_normalA[z+1].x==0 && tmp_normalA[z+1].y==0)) continue;
		normalSag[z] = (tmp_normalA[z-1]+tmp_normalA[z]+tmp_normalA[z+1]).normalize();
	}
*/
	// smooth the normal using piecewise b-spline
	for(z=0; z<sizez; z++)
	{
		xt[z] = z;
		yt[z] = -1;
	}

	for(z=0; z<sizez; z++)
	{
		if(normalSag[z].x!=0 || normalSag[z].y!=0)
		{
			yt[z] = atan2(normalSag[z].y, normalSag[z].x);
		}
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=0; z<spinalCord3D.GetSize(); z++)
	{
		if(normalSag[z].x!=0 || normalSag[z].y!=0)
		{
			normalSag[z].x = cos(yt2[z]);
			normalSag[z].y = sin(yt2[z]);
		}
	}

	// assign globle variable
	for(z=0; z<sizez; z++)
	{
		if(normalSag[z].x<0) 
		{
			normalSag[z].x = -normalSag[z].x;
			normalSag[z].y = -normalSag[z].y;
		}
		spineNormalX[z] = Vec3(0, normalSag[z].x, normalSag[z].y);
	}

	// coronal plane
	projectedCord->UpdateControlList(projCor);

	// compute the normal from the b-spline curve
	curveLength=projectedCord->Length();
	stepSize=curveLength/projCor.GetSize();
	s=0;
	for(z=0;z<sizez;z++)
	{
		k = filledZ[z];
		if(k<0) continue;
		normalCor[z] = projectedCord->Normal(s);
		s += stepSize;
	}

	// smooth the normal using Laplacian method
/*	tmp_normalA = normalCor;
	for(z=1; z<sizez-1; z++)
	{
		if((tmp_normalA[z-1].x==0 && tmp_normalA[z-1].y==0) || (tmp_normalA[z+1].x==0 && tmp_normalA[z+1].y==0)) continue;
		normalCor[z] = (tmp_normalA[z-1]+tmp_normalA[z]+tmp_normalA[z+1]).normalize();
	}
*/
	// smooth the normal using piecewise b-spline
	for(z=0; z<sizez; z++)
	{
		xt[z] = z;
		yt[z] = -1;
	}

	for(z=0; z<sizez; z++)
	{
		if(normalCor[z].x!=0 || normalCor[z].y!=0)
		{
			yt[z] = atan2(normalCor[z].y, normalCor[z].x);
		}
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=0; z<spinalCord3D.GetSize(); z++)
	{
		if(normalCor[z].x!=0 || normalCor[z].y!=0)
		{
			normalCor[z].x = cos(yt2[z]);
			normalCor[z].y = sin(yt2[z]);
		}
	}

	// assign globle variable
	for(z=0; z<sizez; z++)
	{
		if(normalCor[z].x<0) 
		{
			normalCor[z].x = -normalCor[z].x;
			normalCor[z].y = -normalCor[z].y;
		}
		spineNormalY[z] = Vec3(normalCor[z].x, 0, normalCor[z].y);
	}


	delete projectedCord;

	// compute the intensity integral along the normal direction on coronal view
	intDynArray leftAvgI, rightAvgI, leftCount, rightCount, leftAvgI1, rightAvgI1, leftAvgI2, rightAvgI2, allCount, allAvgI, allAvgI2;
	int count0;
	short gray;
	leftAvgI.SetSize(sizez);
	leftAvgI1.SetSize(sizez);
	leftAvgI2.SetSize(sizez);
	rightAvgI.SetSize(sizez);
	rightAvgI1.SetSize(sizez);
	rightAvgI2.SetSize(sizez);
	leftCount.SetSize(sizez);
	rightCount.SetSize(sizez);
	allCount.SetSize(sizez);
	allAvgI.SetSize(sizez);
	allAvgI2.SetSize(sizez);

	Vec2 startPos, curPos;
	float posx, posy;
	for(z=0; z<sizez; z++)
	{
		leftAvgI[z] = leftAvgI1[z] = leftAvgI2[z] = 0;
		rightAvgI[z] = rightAvgI1[z] = rightAvgI2[z] = 0;
		leftCount[z] = rightCount[z] = 0;
		allCount[z]=0;
		allAvgI[z]=0;
		allAvgI2[z]=0;
		if(spineCenter[z].x==0) continue;

		// left side
		count0 = count1 = count2 = 0;
		startPos.x = spineCenter[z].x;
		startPos.y = z*pz;

		curPos = startPos;

		for(k=0; k<avgLeftWidth; k++)
		{
			curPos -= normalCor[z];
			posx = curPos.x/px;
			posy = curPos.y/pz;
			x = (int)(posx+0.5);
			y = (int)(posy+0.5);

			gray = (short)img2D_cor->GetPixelAtPosition(posx, posy);
			mask = img2D_cor_mask->FastGet(x, y);

			if(gray<=0) continue;

			if(mask==maskStruct.mask_corticalBone || mask==maskStruct.mask_spongyBone)
			{
				count0++;
				leftAvgI[z]+=gray;
			}
				
			if(mask!=maskStruct.mask_spinalCord)
			{
				count1++;
				leftAvgI1[z]+=gray;
				allCount[z]++;
				allAvgI[z]+=gray;
			}

			count2++;
			leftAvgI2[z] += gray;
		}

		leftCount[z] = count1;
		if(count0!=0) leftAvgI[z]/= count0;
		if(count1!=0) leftAvgI1[z]/= count1;
		if(count2!=0) leftAvgI2[z]/= count2;

		// right side
		count0 = count1 = count2 = 0;

		curPos = startPos;
		for(k=0; k<avgRightWidth; k++)
		{
			curPos += normalCor[z];
			posx = curPos.x/px;
			posy = curPos.y/pz;
			x = (int)(posx+0.5);
			y = (int)(posy+0.5);

			gray = (short)img2D_cor->GetPixelAtPosition(posx, posy);
			mask = img2D_cor_mask->FastGet(x, y);

			if(gray<=0) continue;

			if(mask==maskStruct.mask_corticalBone || mask==maskStruct.mask_spongyBone)
			{
				count0++;
				rightAvgI[z]+=gray;
			}
				
			if(mask!=maskStruct.mask_spinalCord)
			{
				count1++;
				rightAvgI1[z]+=gray;
				allCount[z]++;
				allAvgI[z]+=gray;
			}

			count2++;
			rightAvgI2[z] += gray;
		}

		rightCount[z] = count1;
		if(count0!=0) rightAvgI[z]/= count0;
		if(count1!=0) rightAvgI1[z]/= count1;
		if(count2!=0) rightAvgI2[z]/= count2;
			
	}	// for z

	// compute the average of neighboring three slices
	//
	if(allCount[0]+allCount[1]>0)
		allAvgI2[0] = (allAvgI[0]+allAvgI[1])/(allCount[0]+allCount[1]);
	for(z=1; z<sizez-1; z++)
	{
		if(allCount[z-1]+allCount[z]+allCount[z+1]>0)
			allAvgI2[z] = (allAvgI[z-1]+allAvgI[z]+allAvgI[z+1])/(allCount[z-1]+allCount[z]+allCount[z+1]);
	}
	if(allCount[sizez-1]+allCount[sizez-2]>0)
		allAvgI2[sizez-1] = (allAvgI[sizez-1]+allAvgI[sizez-2])/(allCount[sizez-1]+allCount[sizez-2]);

	for(z=0; z<sizez; z++)
	{
		if(allCount[z]!=0) allAvgI[z] /= allCount[z];
	}

	// compute the intensity integral along the normal direction on sagital view
	// also adjust the normal to get best cut on the sagital view
	intDynArray upAvgI, upCount, upAvgI2, upAvgIm;
	upAvgI.SetSize(sizez);
	upAvgI2.SetSize(sizez);
	upAvgIm.SetSize(sizez);
	upCount.SetSize(sizez);

	float angleSearchRange = 20*3.14159/180;
	float angleStep = 1*3.14159/180;
	for(z=0; z<sizez; z++)
	{
		upAvgI[z] = upAvgI2[z] = 0;
		upAvgIm[z]=0;
		upCount[z] = 0;
		if(spineCenter[z].x==0) continue;

		// up side
		count0 = count1 = count2 = 0;
		startPos.x = spineCenter[z].y;
		startPos.y = z*pz;

		curPos = startPos;

		int searchRange=(int)(upWidth[z]*1.5);
		// adjust the angle within 15 degree with 2 degree step
		float midAngle = atan2(-normalSag[z].y, -normalSag[z].x);
		float angleRange1 = midAngle-angleSearchRange;
		float angleRange2 = midAngle+angleSearchRange;

		float cAngle, mAngle;
		Vec2 cNormal;
		float mAccu=1e10, cAccu, cAccum, mAccum=0;
		int mCount=0;

		for(cAngle=angleRange1; cAngle<=angleRange2; cAngle+=angleStep)
		{
			cNormal = Vec2(cos(cAngle), sin(cAngle));
			cAccu = 0;
			cAccum = 0;
			count0 = 0;
			count1 = 0;
			curPos = startPos;
			for(k=0; k<searchRange; k++)
			{
				curPos += cNormal;
				posx = curPos.x/py;
				posy = curPos.y/pz;
				x = (int)(posx+0.5);
				y = (int)(posy+0.5);

				mask = img2D_sag_mask->FastGet(x, y);
		
				if(mask==maskStruct.mask_corticalBone || mask==maskStruct.mask_spongyBone)
				{
					gray = (short)img2D_sag->GetPixelAtPosition(posx, posy);
					if(gray<=0) break;
					count1++;
					cAccum += gray;
					if(count1<=7) continue;	// skip the first 7 pixels, which somtimes could be big for disk
					count0++;
					cAccu += gray;
				}
			}
			if(count1>15)
			{
				cAccu /= count0;
				cAccum /= count1;
				if(cAccu<mAccu)
				{
					mAccu = cAccu;
					mCount = count0;
					mAngle= cAngle;
				}
				if(cAccum>mAccum)
				{
					mAccum = cAccum;
				}
			}
		}	// for cAngle

		if(mCount>0)
		{
			upCount[z] = mCount;
			upAvgI[z] = mAccu;
			normalSag[z].x = -cos(mAngle);
			normalSag[z].y = -sin(mAngle);

			upAvgIm[z] = mAccum;	// count all pixels
		}
	}	// for z

	// compute the average of neighboring three slices
	//
	for(z=0; z<sizez; z++) upAvgI2[z]=upAvgI[z];
	if(upCount[0]>0 && upCount[1]>0)
		upAvgI2[0] = (upAvgI[0]+upAvgI[1])/2;
	for(z=1; z<sizez-1; z++)
	{
		if(upCount[z-1]>0 && upCount[z+1]>0 && upCount[z]>0) 
			upAvgI2[z] = (upAvgI[z-1]+upAvgI[z]+upAvgI[z+1])/3;
		else upAvgI2[z] = upAvgI[z];
	}
	if(upCount[sizez-1]>0 && upCount[sizez-2]>0)
		upAvgI2[sizez-1] = (upAvgI[sizez-1]+upAvgI[sizez-2])/2;

	// assign globle variable
	for(z=0; z<sizez; z++)
	{
		if(normalSag[z].x<0) 
		{
			normalSag[z].x = -normalSag[z].x;
			normalSag[z].y = -normalSag[z].y;
		}
		spineNormalX[z] = Vec3(0, normalSag[z].x, normalSag[z].y);
	}

	float estimate_gap = 20;
	int estimate_int = (int)(estimate_gap/pz+0.5);
	int minP, minN, zz;


	// two methods to partition
	//	for low resolution data (pz>2), use the coronal view
	//	for high resolution data (pz<2), use the sagital view
	//
	pedicleValley.SetSize(0);
	if(pz>2)	//	low resolution
	{
		estimate_int-=1;
		// try to locate the valleys on coronal view
		for(z=1; z<sizez-1; z++)
		{
			minP=allAvgI[z-1];
			for(zz=z-1; zz>=0 && zz>=z-estimate_int; zz--) if(allAvgI[zz]<minP) minP=allAvgI[zz];
			minN=allAvgI[z+1];
			for(zz=z+1; zz<sizez && zz<=z+estimate_int; zz++) if(allAvgI[zz]<minN) minN=allAvgI[zz];

			// locate the valley
			if(allAvgI[z]<minP && allAvgI[z]<=minN)
			{
				pedicleValley.Add(z);
			}
		}
	}
	else {	// high resolution data

		// try to locate the disk on sag view

		for(z=1; z<sizez-1; z++)
		{
			if(upAvgI2[z]==0) continue;
			minP=upAvgI2[z-1];
			for(zz=z-1; zz>=0 && zz>=z-estimate_int; zz--) if(upAvgI2[zz]>0 && upAvgI2[zz]<minP) minP=upAvgI2[zz];
			minN=upAvgI2[z+1];
			for(zz=z+1; zz<sizez && zz<=z+estimate_int; zz++) if(upAvgI2[zz]>0 && upAvgI2[zz]<minN) minN=upAvgI2[zz];

			// locate the valley
			if(upAvgI2[z]<minP && upAvgI2[z]<=minN)
			{
				pedicleValley.Add(z);
			}
		}
	}

	double avgInterval;
	avgInterval = 0;
	for(k=1; k<pedicleValley.GetSize(); k++)
	{
		avgInterval += pedicleValley[k]-pedicleValley[k-1];
	}
	avgInterval /= (double)(pedicleValley.GetSize()-1);
	avgInterval *= pz;

	// insert a valley if the interval is too big
	for(k=1; k<pedicleValley.GetSize(); k++)
	{
		float gap = (pedicleValley[k]-pedicleValley[k-1])*pz;
		if( gap>= avgInterval*1.75 ||
			(avgInterval<estimate_gap*2.5 && gap>estimate_gap*2.5))
		{
			int tz=(pedicleValley[k]+pedicleValley[k-1])/2;
			for(z=tz-estimate_int/2; z<tz+estimate_int/2; z++)
			{
				if(z<0 || z>=sizez || upAvgI2[z]==0) continue;
				minP=upAvgI2[z-1];
				for(zz=z-1; zz>=0 && zz>=z-estimate_int/2; zz--) if(upAvgI2[zz]>0 && upAvgI2[zz]<minP) minP=upAvgI2[zz];
				minN=upAvgI2[z+1];
				for(zz=z+1; zz<sizez && zz<=z+estimate_int/2; zz++) if(upAvgI2[zz]>0 && upAvgI2[zz]<minN) minN=upAvgI2[zz];

				// locate the valley
				if(upAvgI2[z]<minP && upAvgI2[z]<=minN)
				{
					pedicleValley.InsertAt(k,z);
					break;
				}
			}
		}
	}

	// check if I need to insert one valley at the very beginning
	if(pedicleValley.GetSize()>3)
	{
		// average interval of first two vertebra
		int avI = (pedicleValley[2]-pedicleValley[0])/2;
		if(pedicleValley[0]>avI+2)
		{
			int tz = pedicleValley[0]-avI;
			pedicleValley.InsertAt(0,tz);
		}
	}

	if(debugMode)
	{
		FILE *fp;
		fp = fopen("d:\\tmp\\integral.txt", "w");
		fprintf(fp, "leftCount leftAvg leftAvg1 leftAvg2 rightCount rightAvg rightAvg1 rightAvg2 allCount allAvg allAvg2 upCount upAvg upAvg2 upAvgm\n");
		for(z=0;z<sizez; z++)
		{
			fprintf(fp, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", 
				leftCount[z], leftAvgI[z], leftAvgI1[z], leftAvgI2[z], 
				rightCount[z], rightAvgI[z], rightAvgI1[z], rightAvgI2[z],
				allCount[z], allAvgI[z], allAvgI2[z], upCount[z], upAvgI[z], upAvgI2[z], upAvgIm[z]);
		}
		// print out the pedicles
		for(z=0; z<pedicleValley.GetSize(); z++)
		{
			fprintf(fp, "%d\n", pedicleValley[z]);
		}
		fclose(fp);

	}

	// in low resolution case
	// smooth the normals at pedicle using piecewise b-spline
	if(pz>2)
	{
		for(z=0; z<sizez; z++)
		{
			xt[z] = z;
			yt[z] = -1;
		}

		for(k=0; k<pedicleValley.GetSize(); k++)
		{
			z = pedicleValley[k];
			if(normalSag[z].y!=0 || normalSag[z].x!=0)
				yt[z] = atan2(normalSag[z].y, normalSag[z].x);
			else yt[z]=-1;
		}
		PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
		for(k=0; k<pedicleValley.GetSize(); k++)
		{
			z = pedicleValley[k];
			normalSag[z].x = cos(yt2[z]);
			normalSag[z].y = sin(yt2[z]);
			spineNormalX[z].y = normalSag[z].x;
			spineNormalX[z].z = normalSag[z].y;
		}
	}

	// compute the partition of spinous process (back of vertebra)
	// only for high resolution cases

	if(pz<2)
	{
		// compute the intensity integral along the normal direction on sagital view
		// also adjust the normal to get best cut on the sagital view
		vec2DynArray normalSagBack;
		normalSagBack.SetSize(sizez);
		intDynArray backAvgI, backCount, backAvgI2;
		backAvgI.SetSize(sizez);
		backAvgI2.SetSize(sizez);
		backCount.SetSize(sizez);

		for(z=0; z<sizez; z++)
		{
			normalSagBack[z] = normalSag[z];
			backAvgI[z] = backAvgI2[z] = 0;
			backCount[z] = 0;
			if(spineCenter[z].x==0) continue;

			// back side
			count0 = count1 = count2 = 0;
			startPos.x = spineCenter[z].y;
			startPos.y = z*pz;

			curPos = startPos;

			int searchRange=(int)(upWidth[z]);
			// adjust the angle within 15 degree with 2 degree step
			float midAngle = atan2(normalSag[z].y, normalSag[z].x)+3.14159/4;
			float angleSearchRange = 15*3.14159/180;
			float angleStep = 1*3.14159/180;
			float angleRange1 = midAngle-angleSearchRange;
			float angleRange2 = midAngle+angleSearchRange;

			float cAngle, mAngle;
			Vec2 cNormal;
			float mAccu=1e10, cAccu;
			int mCount=0;

			for(cAngle=angleRange1; cAngle<=angleRange2; cAngle+=angleStep)
			{
				cNormal = Vec2(cos(cAngle), sin(cAngle));
				cAccu = 0;
				count0 = 0;
				count1 = 0;
				curPos = startPos;
				for(k=0; k<searchRange; k++)
				{
					curPos += cNormal;
					posx = curPos.x/py;
					posy = curPos.y/pz;
					x = (int)(posx+0.5);
					y = (int)(posy+0.5);
	
					mask = img2D_sag_mask->FastGet(x, y);
		
					if(mask==maskStruct.mask_body || mask==maskStruct.mask_corticalBone || mask==maskStruct.mask_spongyBone)
					{
						gray = (short)img2D_sag->GetPixelAtPosition(posx, posy);
						if(gray<=0) break;
						count0++;
						cAccu += gray;
					}
				}
				if(count0>20)
				{
					cAccu /= count0;
					if(cAccu<mAccu)
					{
						mAccu = cAccu;
						mCount = count0;
						mAngle= cAngle;
					}
				}
			}	// for cAngle

			if(mCount>0)
			{
				backCount[z] = mCount;
				backAvgI[z] = mAccu;
				normalSagBack[z].x = cos(mAngle);
				normalSagBack[z].y = sin(mAngle);
			}
		}	// for z

		// assign globle variable
		spineNormalBack.SetSize(sizez);
		for(z=0; z<sizez; z++)
		{
			spineNormalBack[z] = Vec3(0, normalSagBack[z].x, normalSagBack[z].y);
		}

		int stz, edz, minz, lastz=-1;
		backValley.SetSize(0);
		for(zz=-1; zz<pedicleValley.GetSize()-1; zz++)
		{
			if(zz==-1) 
			{
				stz=0;
				edz=pedicleValley[0];
			}
			else 
			{
				stz=pedicleValley[zz];
				edz = pedicleValley[zz+1];

				if(lastz>0 && stz<lastz+(pedicleValley[zz+1]-pedicleValley[zz])/2) stz=lastz+(pedicleValley[zz+1]-pedicleValley[zz])/2;
				if(lastz>0 && edz>lastz+(pedicleValley[zz+1]-pedicleValley[zz])*3/2) edz=lastz+(pedicleValley[zz+1]-pedicleValley[zz])*3/2;
			}

			minP = 10000;
			minz = -1;
			for(z=stz; z<edz; z++)
			{
				if(backAvgI[z]>0 && backAvgI[z]<minP)
				{
					minP = backAvgI[z];
					minz = z;
				}
			}
			if(minz!=-1) 
			{
				backValley.Add(minz);
				lastz = minz;
			}
		}	//for zz
	}	// if pz

	return CIS_OK;
}



// new algorithm for partition
//
int NIH_SpineSegmentation_SpinePartition_new(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   VertebraStruct2D *vertebraTemplate, 
									   vec3DynArray &spinalCord3D,
									   CIS_Array_Image2D_short *img2D_sag,
									   CIS_Array_Image2D_short *img2D_cor,
									   CIS_Array_Image2D_short *img2D_sag_mask,
									   CIS_Array_Image2D_short *img2D_cor_mask,
									   vec3DynArray &spineCenter, vec3DynArray &spineNormalX, vec3DynArray &spineNormalY, vec3DynArray &spineNormalBack, 
									   doubleDynArray &spineWidthLeft, doubleDynArray &spineWidthRight, doubleDynArray &spineWidthUp, doubleDynArray &spineWidthDown, 
									   intDynArray &pedicleValley, intDynArray &backValley, 
									   bool debugMode)
{
	// spinalCord not extracted
	if(spinalCord3D.GetSize()<10) return CIS_ERROR;

	int sizex, sizey, sizez, sizexy;
	short *maskA, *imgA;
	int x, y, z, k;
	float px, py, pz;

	sizex = maskImg3D->Num_Cols();
	sizey = maskImg3D->Num_Rows();
	sizez = maskImg3D->Num_Levels();
	px = img3D->Get_Pixel_SizeX();
	py = img3D->Get_Pixel_SizeY();
	pz = img3D->Get_Pixel_SizeZ();
	sizexy = sizex*sizey;

	maskA = maskImg3D->GetArray();
	imgA = img3D->GetArray();

	float vsizex = img3D->Get_Volume_SizeX();
	float vsizey = img3D->Get_Volume_SizeY();
	float vsizez = img3D->Get_Volume_SizeZ();

	// compute the projected curve on sagittal plane (x plane)
	intDynArray leftWidth, rightWidth, upWidth, downWidth;
	intDynArray filledZ;
	vec2DynArray projSag, projCor;
	vec2DynArray normalSag, normalCor;

	filledZ.SetSize(sizez);

	int sx, sy;
	short mask;

	upWidth.SetSize(sizez);
	downWidth.SetSize(sizez);
	leftWidth.SetSize(sizez);
	rightWidth.SetSize(sizez);
		
	spineCenter.SetSize(sizez);
	spineNormalX.SetSize(sizez);
	spineNormalY.SetSize(sizez);
	spineWidthUp.SetSize(sizez);
	spineWidthDown.SetSize(sizez);
	spineWidthLeft.SetSize(sizez);
	spineWidthRight.SetSize(sizez);

	for(z=0; z<sizez; z++)
	{
		spineCenter[z] = Vec3(0,0,0);
		spineNormalX[z] = Vec3(0,0,0);
		spineNormalY[z] = Vec3(0,0,0);
		spineWidthUp[z] = -1;
		spineWidthDown[z] = -1;
		spineWidthLeft[z] = -1;
		spineWidthRight[z] = -1;
	}

	for(z=0; z<sizez; z++) 
	{
		filledZ[z]=-1;
		rightWidth[z] = leftWidth[z]=-1;
		upWidth[z] = downWidth[z]=-1;
	}
		
	projSag.SetSize(spinalCord3D.GetSize());
	projCor.SetSize(spinalCord3D.GetSize());

	for(k=0; k<spinalCord3D.GetSize(); k++)
	{
		z = img3D->GetSliceLevel(spinalCord3D[k].z,1);

		projSag[k].x = spinalCord3D[k].y;
		projSag[k].y = z*pz;
		projCor[k].x = spinalCord3D[k].x;
		projCor[k].y = z*pz;

		if(filledZ[z]==1 && filledZ[z-1]==-1) z=z-1;
		while(z<sizez-1 && filledZ[z]!=-1) z++;
		filledZ[z] = k;

		spineCenter[z].x = spinalCord3D[k].x;
		spineCenter[z].y = spinalCord3D[k].y;
		spineCenter[z].z = img3D->GetSlicePosition(z);

		// set up the spinal column
		x = (int)(spinalCord3D[k].x/px+0.5);
		y = (int)(spinalCord3D[k].y/py+0.5);
		// to the right
		int count = 0;
		for(sx=x; sx<sizex; sx++)
		{
			mask = maskImg3D->FastGet(sx, y, z);
			if(mask!=maskStruct.mask_spongyBone && mask!=maskStruct.mask_corticalBone && mask!=maskStruct.mask_spinalCord)
			{
				count++;
			}
			else count=0;

			if(count==3)	// search for the border, break when meet 3 non-vertebra pixels
			{
				rightWidth[z] = sx-x;
				break;
			}
		}

		// to the left
		count =0;
		for(sx=x; sx>0; sx--)
		{
			mask = maskImg3D->FastGet(sx, y, z);
			if(mask!=maskStruct.mask_spongyBone && mask!=maskStruct.mask_corticalBone && mask!=maskStruct.mask_spinalCord)
			{
				count++;
			}
			else count=0;

			if(count==3)	// search for the border, break when meet 3 non-vertebra pixels
			{
				leftWidth[z] = x-sx;
				break;
			}
		}

		// set up the spinal column
		x = (int)(spinalCord3D[k].x/px+0.5);
		y = (int)(spinalCord3D[k].y/py+0.5);
		// upward
		for(sy=y; sy>0; sy--)
		{
			mask = maskImg3D->FastGet(x, sy, z);
			if(mask==maskStruct.mask_body)		
			{
				upWidth[z] = y-sy;	// stop at the first non-vertebra pixel
				break;
			}
		}

		// downward
		for(sy=y; sy<=y+upWidth[z]*2; sy++)
		{
			for(sx=x-leftWidth[z]; sx<=x+rightWidth[z]; sx++)
			{
				mask = maskImg3D->FastGet(sx, sy, z);
				if(mask==maskStruct.mask_spongyBone || mask==maskStruct.mask_corticalBone || mask==maskStruct.mask_spinalCord)
				{
					downWidth[z] = sy-y;	// stop at the last vertebra pixel
				}
			}
		}
	}	// for k

	// smooth upWidth and downWidth using piecewise b-spine
	//
	doubleDynArray xt, yt, yt2;
	int bernstein_power = 5, piece_size;
	double minInterpolate, maxInterpolate, minInterpolate2, maxInterpolate2;
	float piece_length=25;
	piece_size = (int)(piece_length/pz+0.5);

	xt.SetSize(sizez);
	yt.SetSize(sizez);
	yt2.SetSize(sizez);
	for(z=0; z<sizez; z++)
	{
		xt[z] = z;
		yt[z] = -1;
	}
	// upWidth
	minInterpolate = -1;
	maxInterpolate = -1;
	for(z=0; z<sizez; z++)
	{
		if(upWidth[z]>(segInfo.diskBound2[z].y-segInfo.diskBound1[z].y)/2) yt[z] = upWidth[z];
		if(yt[z]==-1) continue;
		if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
		if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=0; z<sizez; z++)
	{
		upWidth[z]=segInfo.diskBound2[z].y-segInfo.diskBound1[z].y;
		if(yt2[z]==-1) continue;
		if(yt2[z]<minInterpolate) yt2[z]=minInterpolate;
		if(yt2[z]>maxInterpolate) yt2[z]=maxInterpolate;
		upWidth[z] = yt2[z];
	}

	// DownWidth
	minInterpolate = -1;
	maxInterpolate = -1;
	for(z=0; z<sizez; z++)
	{
		yt[z] = downWidth[z];
		if(yt[z]==-1) continue;
		if(yt[z]<minInterpolate || minInterpolate==-1) minInterpolate=yt[z];
		if(yt[z]>maxInterpolate || maxInterpolate==-1) maxInterpolate=yt[z];
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=0; z<sizez; z++)
	{
		if(yt2[z]==-1) continue;
		if(yt2[z]<minInterpolate) yt2[z]=minInterpolate;
		if(yt2[z]>maxInterpolate) yt2[z]=maxInterpolate;
		downWidth[z] = yt2[z];
	}

	// smooth upWidth and downWidth
	// laplacian smoothing
	intDynArray tmpWidth;
	tmpWidth = upWidth;
	for(z=1; z<sizez-1; z++)
	{
		if(tmpWidth[z-1]==0 || tmpWidth[z+1]==0) continue;
		if(tmpWidth[z]<tmpWidth[z-1] && tmpWidth[z]<tmpWidth[z+1]) upWidth[z] = (tmpWidth[z-1]+tmpWidth[z+1])/2;
		else upWidth[z] = (tmpWidth[z-1]+tmpWidth[z]+tmpWidth[z+1])/3;
	}
	tmpWidth = downWidth;
	for(z=1; z<sizez-1; z++)
	{
		if(tmpWidth[z-1]==0 || tmpWidth[z+1]==0) continue;
		downWidth[z] = (tmpWidth[z-1]+tmpWidth[z]+tmpWidth[z+1])/3;
	}

	// assign globle variable, return with the function
	for(z=0; z<sizez; z++)
	{
		spineWidthUp[z] = upWidth[z]*py;
		spineWidthDown[z] = downWidth[z]*py;
	}

	// compute the left and right width
	int avgLeftWidth, avgRightWidth;
	int count, count1, count2;;
	avgLeftWidth = avgRightWidth = 0;
	count = 0;
	for(z=0; z<sizez; z++)
	{
		if(leftWidth[z]>0 && rightWidth[z]>0)
		{
			avgLeftWidth += leftWidth[z];
			avgRightWidth += rightWidth[z]; 
			count ++;
		}
	}
	 	
	if(count!=0)
	{
		avgLeftWidth /= count;
		avgRightWidth /= count;

		avgLeftWidth += 4;
		avgRightWidth += 4;
	}
		
	// assign globle variables
	for(z=0; z<sizez; z++)
	{
		spineWidthLeft[z] = avgLeftWidth*px;
		spineWidthRight[z] = avgRightWidth*px;
//		spineWidthLeft[z] = leftWidth[z]*px;
//		spineWidthRight[z] = rightWidth[z]*px;
	}

	// compute the normal in sag and cor plane separately
	//
	CIS_2D_Model_Curve *projectedCord;
	projectedCord = new CIS_2D_Model_Curve();

	normalSag.SetSize(sizez);
	normalCor.SetSize(sizez);
	for(z=0; z<sizez; z++)
	{
		normalSag[z] = Vec2(0,0);
		normalCor[z] = Vec2(0,0);
	}

	// sagittal plane
	projectedCord->UpdateControlList(projSag);

	// compute the normal from the b-spline curve
	double curveLength=projectedCord->Length();
	double stepSize=curveLength/projSag.GetSize();
	double s=0;
	for(z=0;z<sizez;z++)
	{
		k = filledZ[z];
		if(k<0) continue;
		normalSag[z] = projectedCord->Normal(s);
		s += stepSize;
	}

	// smooth the normal using Laplacian method
/*	vec2DynArray tmp_normalA;
	tmp_normalA = normalSag;
	for(z=1; z<sizez-1; z++)
	{
		if((tmp_normalA[z-1].x==0 && tmp_normalA[z-1].y==0) || (tmp_normalA[z+1].x==0 && tmp_normalA[z+1].y==0)) continue;
		normalSag[z] = (tmp_normalA[z-1]+tmp_normalA[z]+tmp_normalA[z+1]).normalize();
	}
*/
	// smooth the normal using piecewise b-spline
	for(z=0; z<sizez; z++)
	{
		xt[z] = z;
		yt[z] = -1;
	}

	for(z=0; z<sizez; z++)
	{
		if(normalSag[z].x!=0 || normalSag[z].y!=0)
		{
			yt[z] = atan2(normalSag[z].y, normalSag[z].x);
		}
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=0; z<spinalCord3D.GetSize(); z++)
	{
		if(normalSag[z].x!=0 || normalSag[z].y!=0)
		{
			normalSag[z].x = cos(yt2[z]);
			normalSag[z].y = sin(yt2[z]);
		}
	}

	// assign globle variable
	for(z=0; z<sizez; z++)
	{
		if(normalSag[z].x<0) 
		{
			normalSag[z].x = -normalSag[z].x;
			normalSag[z].y = -normalSag[z].y;
		}
		spineNormalX[z] = Vec3(0, normalSag[z].x, normalSag[z].y);
	}

	// coronal plane
	projectedCord->UpdateControlList(projCor);

	// compute the normal from the b-spline curve
	curveLength=projectedCord->Length();
	stepSize=curveLength/projCor.GetSize();
	s=0;
	for(z=0;z<sizez;z++)
	{
		k = filledZ[z];
		if(k<0) continue;
		normalCor[z] = projectedCord->Normal(s);
		s += stepSize;
	}

	// smooth the normal using Laplacian method
/*	tmp_normalA = normalCor;
	for(z=1; z<sizez-1; z++)
	{
		if((tmp_normalA[z-1].x==0 && tmp_normalA[z-1].y==0) || (tmp_normalA[z+1].x==0 && tmp_normalA[z+1].y==0)) continue;
		normalCor[z] = (tmp_normalA[z-1]+tmp_normalA[z]+tmp_normalA[z+1]).normalize();
	}
*/
	// smooth the normal using piecewise b-spline
	for(z=0; z<sizez; z++)
	{
		xt[z] = z;
		yt[z] = -1;
	}

	for(z=0; z<sizez; z++)
	{
		if(normalCor[z].x!=0 || normalCor[z].y!=0)
		{
			yt[z] = atan2(normalCor[z].y, normalCor[z].x);
		}
	}
	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=0; z<spinalCord3D.GetSize(); z++)
	{
		if(normalCor[z].x!=0 || normalCor[z].y!=0)
		{
			normalCor[z].x = cos(yt2[z]);
			normalCor[z].y = sin(yt2[z]);
		}
	}

	// assign globle variable
	for(z=0; z<sizez; z++)
	{
		if(normalCor[z].x<0) 
		{
			normalCor[z].x = -normalCor[z].x;
			normalCor[z].y = -normalCor[z].y;
		}
		spineNormalY[z] = Vec3(normalCor[z].x, 0, normalCor[z].y);
	}


	delete projectedCord;

	// compute the intensity integral along the normal direction on coronal view
	intDynArray leftAvgI, rightAvgI, leftCount, rightCount, leftAvgI1, rightAvgI1, leftAvgI2, rightAvgI2, allCount, allAvgI, allAvgI2;
	int count0;
	short gray;
	leftAvgI.SetSize(sizez);
	leftAvgI1.SetSize(sizez);
	leftAvgI2.SetSize(sizez);
	rightAvgI.SetSize(sizez);
	rightAvgI1.SetSize(sizez);
	rightAvgI2.SetSize(sizez);
	leftCount.SetSize(sizez);
	rightCount.SetSize(sizez);
	allCount.SetSize(sizez);
	allAvgI.SetSize(sizez);
	allAvgI2.SetSize(sizez);

	Vec2 startPos, curPos;
	float posx, posy;
	for(z=0; z<sizez; z++)
	{
		leftAvgI[z] = leftAvgI1[z] = leftAvgI2[z] = 0;
		rightAvgI[z] = rightAvgI1[z] = rightAvgI2[z] = 0;
		leftCount[z] = rightCount[z] = 0;
		allCount[z]=0;
		allAvgI[z]=0;
		allAvgI2[z]=0;
		if(spineCenter[z].x==0) continue;

		// left side
		count0 = count1 = count2 = 0;
		startPos.x = spineCenter[z].x;
		startPos.y = z*pz;

		curPos = startPos;

		for(k=0; k<avgLeftWidth; k++)
		{
			curPos -= normalCor[z];
			posx = curPos.x/px;
			posy = curPos.y/pz;
			x = (int)(posx+0.5);
			y = (int)(posy+0.5);

			gray = (short)img2D_cor->GetPixelAtPosition(posx, posy);
			mask = img2D_cor_mask->FastGet(x, y);

			if(gray<=0) continue;

			if(mask==maskStruct.mask_corticalBone || mask==maskStruct.mask_spongyBone)
			{
				count0++;
				leftAvgI[z]+=gray;
			}
				
			if(mask!=maskStruct.mask_spinalCord)
			{
				count1++;
				leftAvgI1[z]+=gray;
				allCount[z]++;
				allAvgI[z]+=gray;
			}

			count2++;
			leftAvgI2[z] += gray;
		}

		leftCount[z] = count1;
		if(count0!=0) leftAvgI[z]/= count0;
		if(count1!=0) leftAvgI1[z]/= count1;
		if(count2!=0) leftAvgI2[z]/= count2;

		// right side
		count0 = count1 = count2 = 0;

		curPos = startPos;
		for(k=0; k<avgRightWidth; k++)
		{
			curPos += normalCor[z];
			posx = curPos.x/px;
			posy = curPos.y/pz;
			x = (int)(posx+0.5);
			y = (int)(posy+0.5);

			gray = (short)img2D_cor->GetPixelAtPosition(posx, posy);
			mask = img2D_cor_mask->FastGet(x, y);

			if(gray<=0) continue;

			if(mask==maskStruct.mask_corticalBone || mask==maskStruct.mask_spongyBone)
			{
				count0++;
				rightAvgI[z]+=gray;
			}
				
			if(mask!=maskStruct.mask_spinalCord)
			{
				count1++;
				rightAvgI1[z]+=gray;
				allCount[z]++;
				allAvgI[z]+=gray;
			}

			count2++;
			rightAvgI2[z] += gray;
		}

		rightCount[z] = count1;
		if(count0!=0) rightAvgI[z]/= count0;
		if(count1!=0) rightAvgI1[z]/= count1;
		if(count2!=0) rightAvgI2[z]/= count2;
			
	}	// for z

	// compute the average of neighboring three slices
	//
	if(allCount[0]+allCount[1]>0)
		allAvgI2[0] = (allAvgI[0]+allAvgI[1])/(allCount[0]+allCount[1]);
	for(z=1; z<sizez-1; z++)
	{
		if(allCount[z-1]+allCount[z]+allCount[z+1]>0)
			allAvgI2[z] = (allAvgI[z-1]+allAvgI[z]+allAvgI[z+1])/(allCount[z-1]+allCount[z]+allCount[z+1]);
	}
	if(allCount[sizez-1]+allCount[sizez-2]>0)
		allAvgI2[sizez-1] = (allAvgI[sizez-1]+allAvgI[sizez-2])/(allCount[sizez-1]+allCount[sizez-2]);

	for(z=0; z<sizez; z++)
	{
		if(allCount[z]!=0) allAvgI[z] /= allCount[z];
	}

	float estimate_gap = 20;
	int estimate_int = (int)(estimate_gap/pz+0.5);

	// do a preliminary testing to mark potential vertebral disk
	// the test is done on curve reformated saggittal view
	// those on disk map will be treated differently when partition the spine
	short *maskA_sag;
	short *imgA_sag;
	short *diskMap;
	int k4;
	maskA_sag = img2D_sag_mask->GetArray();
	imgA_sag = img2D_sag->GetArray();
/*	diskMap = new short[sizey*sizez];
	for(k=0; k<sizey*sizez; k++) diskMap[k]=0;

	int diskT=7, corticalT=2;

	for(z=segInfo.bound1.z+10; z<segInfo.bound2.z-10; z++)
	{
///		Vec2 normal = Vec2(-normalSag[z].y, -normalSag[z].x);
///		Vec2 cp, np, np2;
		x = (int)(spineCenter[z].x/px+0.5);

		for(y=segInfo.diskBound1[z].y+2, k=z*sizey+y; y<segInfo.diskBound2[z].y-2; y++, k++)
		{
			if(maskA_sag[k]==maskStruct.mask_spinalCord) break;

			k4 = x+y*sizex+z*sizexy;

///			startPos.x = y*py;
///			startPos.y = z*pz;
///			cp = startPos;
		
			if(maskA_sag[k]==maskStruct.mask_spongyBone)
			{
				// test in four directions
				bool t1, t2, t3, t4, t5, t6;
				int k2, k3;
				t1=t2=t3=t4=t5=t6=true;
				// test left
///				np = cp-normal; np2 = cp-normal*2;
				if(maskA_sag[k-1]!=maskStruct.mask_spongyBone || maskA_sag[k-2]!=maskStruct.mask_spongyBone) 
				{ t1=false; continue;};
				// test right
				if(maskA_sag[k+1]!=maskStruct.mask_spongyBone || maskA_sag[k+2]!=maskStruct.mask_spongyBone) 
				{ t2=false; continue;};

				// test inside
				if(maskA[k4-1]!=maskStruct.mask_spongyBone || maskA[k4-2]!=maskStruct.mask_spongyBone) 
				{ t5=false; continue;};
				// test inside
				if(maskA[k4+1]!=maskStruct.mask_spongyBone || maskA[k4+2]!=maskStruct.mask_spongyBone) 
				{ t6=false; continue;};

				// test up
				for(k2=1; k2<=diskT; k2++)
				{
					if(maskA_sag[k-k2*sizey]==maskStruct.mask_corticalBone) break;
				}
				if(k2<=diskT)
				{
					for(k3=0; k3<corticalT; k3++)
					{
						if(maskA_sag[k-(k3+k2)*sizey]!=maskStruct.mask_corticalBone) break;
					}
					if(k3<corticalT) {t3=false; continue;}
				}
				else {t3=false; continue;};

				// test down
				for(k2=1; k2<=diskT; k2++)
				{
					if(maskA_sag[k+k2*sizey]==maskStruct.mask_corticalBone) break;
//					if(imgA_sag[k+k2*sizey]>1300) break;
				}
				if(k2<=diskT)
				{
					for(k3=0; k3<corticalT; k3++)
					{
						if(maskA_sag[k+(k3+k2)*sizey]!=maskStruct.mask_corticalBone) break;
//						if(imgA_sag[k+(k3+k2)*sizey]<1300) break;
					}
					if(k3<corticalT) {t4=false; continue;}
				}
				else {t4=false; continue;};

				// pass all the test, mark as potential disk
				if(t1 && t2 && t3 && t4 && t5 && t6) diskMap[k]=1;
			}	// if maskA_sag
		}	// for y
	}	// for z
*/
/*	// visualize it for debugging purpose
	for(k=0; k<sizey*sizez; k++) if(diskMap[k]==1) maskA_sag[k]=maskStruct.mask_vertebralDisk;

	return CIS_OK;
*/
	// compute the intensity integral along the normal direction on sagital view
	// also adjust the normal to get best cut on the sagital view
	intDynArray upAvgI, upCount, upAvgI2, upAvgIm;
	upAvgI.SetSize(sizez);
	upAvgI2.SetSize(sizez);
	upAvgIm.SetSize(sizez);
	upCount.SetSize(sizez);

	int x3;
	float angleSearchRange = 15*3.14159/180;
	float angleStep = 1*3.14159/180;
	int k2;
	
/*	for(z=0; z<sizez; z++)
	{
		upAvgI[z] = upAvgI2[z] = 0;
		upAvgIm[z]=0;
		upCount[z] = 0;
		if(spineCenter[z].x==0) continue;

		// up side
		count0 = count1 = count2 = 0;
		startPos.x = spineCenter[z].y;
		startPos.y = z*pz;

		curPos = startPos;

		int searchRange=(int)(upWidth[z]*1.5);
		// adjust the angle within 15 degree with 2 degree step
		float midAngle = atan2(-normalSag[z].y, -normalSag[z].x);
		float angleRange1 = midAngle-angleSearchRange;
		float angleRange2 = midAngle+angleSearchRange;

		float cAngle, mAngle;
		Vec2 cNormal;
		float mAccu=1e10, cAccu, cAccum, mAccum=0;
		int mCount=0;

		for(cAngle=angleRange1; cAngle<=angleRange2; cAngle+=angleStep)
		{
			cNormal = Vec2(cos(cAngle), sin(cAngle));
			cAccu = 0;
			cAccum = 0;
			count0 = 0;
			count1 = 0;
			curPos = startPos;
			for(k=0; k<searchRange; k++)
			{
				curPos += cNormal;
				posx = curPos.x/py;
				posy = curPos.y/pz;
				x = (int)(posx+0.5);
				y = (int)(posy+0.5);
				if(y<0 || y>segInfo.bound2.z-2) break;
				x3 = (int)(spineCenter[y].x/px+0.5);
				if(x3<=0) break;

				k2 = y*sizey+x;
				k4 = x3+y*sizexy+x*sizex;
				mask = maskA_sag[k2];
		
				if(mask==maskStruct.mask_corticalBone || mask==maskStruct.mask_spongyBone)
				{
					int gray1;
//					gray = (short)img2D_sag->GetPixelAtPosition(posx, posy);
					gray1 = imgA[k4];
					if(gray1<=0) break;
					int count=0;
					int k5;
					for(k5=1;k5<spineWidthLeft[y]/px/2; k5++)
					{
						if(maskA[k4-k5]==maskStruct.mask_corticalBone || maskA[k4-k5]==maskStruct.mask_spongyBone)
						{
							gray1 += imgA[k4-k5];
							count++;
						}
					}
					for(k5=1;k5<spineWidthRight[y]/px/2; k5++)
					{
						if(maskA[k4+k5]==maskStruct.mask_corticalBone || maskA[k4+k5]==maskStruct.mask_spongyBone)
						{
							gray1 += imgA[k4+k5];
							count++;
						}
					}

					gray1 /= (count+1);

					if(diskMap[k2]==1) gray1=gray1*3/4;	// compensate for the disk region, only use 3/4 of intensity

					gray = gray1;
					count1++;
					cAccum += gray;
					if(count1<=7) continue;	// skip the first 7 pixels, which somtimes could be big for disk
					count0++;
					cAccu += gray;
				}
			}
			if(count1>15)
			{
				cAccu /= count0;
				cAccum /= count1;
				if(cAccu<mAccu)
				{
					mAccu = cAccu;
					mCount = count0;
					mAngle= cAngle;
				}
				if(cAccum>mAccum)
				{
					mAccum = cAccum;
				}
			}
		}	// for cAngle

		if(mCount>0)
		{
			upCount[z] = mCount;
			upAvgI[z] = mAccu;
			normalSag[z].x = -cos(mAngle);
			normalSag[z].y = -sin(mAngle);

			upAvgIm[z] = mAccum;	// count all pixels
		}
	}	// for z
*/
	// compute the intensity integral along the normal direction on sagital view
	// also adjust the normal to get best cut on the sagital view
	// use information on both the valley and 
	for(z=0; z<sizez; z++)
	{
		upAvgI[z] = upAvgI2[z] = 0;
		upAvgIm[z]=0;
		upCount[z] = 0;
		if(spineCenter[z].x==0) continue;

		// up side
		count0 = count1 = count2 = 0;
		startPos.x = spineCenter[z].y;
		startPos.y = z*pz;

		int searchRange=(int)(upWidth[z]*1.5);
		// adjust the angle within 15 degree with 2 degree step
		float midAngle = atan2(-normalSag[z].y, -normalSag[z].x);
		float angleRange1 = midAngle-angleSearchRange;
		float angleRange2 = midAngle+angleSearchRange;

		float cAngle, mAngle;
		Vec2 cNormal, cPerp;
		float mAccu=1e10;
		int mCount=0;

		for(cAngle=angleRange1; cAngle<=angleRange2; cAngle+=angleStep)
		{
			cNormal = Vec2(cos(cAngle), sin(cAngle));
			cPerp = cNormal.perp().normalize();
			int profileLength;
			
			float score = ComputePartitionAccumulateProfile(img3D, maskImg3D, segInfo, maskStruct, img2D_sag, img2D_sag_mask,
				spineCenter, spineWidthLeft, spineWidthRight, startPos, cNormal, cPerp, searchRange, estimate_gap,
				profileLength);

			if(score<mAccu && profileLength>0)
			{
				mAccu = score;
				mCount = profileLength;
				mAngle= cAngle;
			}
		}	// for cAngle

		if(mCount>0)
		{
			upCount[z] = mCount;
			upAvgI[z] = mAccu;
			normalSag[z].x = -cos(mAngle);
			normalSag[z].y = -sin(mAngle);

///			upAvgIm[z] = mAccum;	// count all pixels
			upAvgIm[z] = mAccu;	// count all pixels
		}
	}	// for z


///	delete diskMap;

	// compute the average of neighboring three slices
	//
	for(z=0; z<sizez; z++) upAvgI2[z]=upAvgI[z];
	if(upCount[0]>0 && upCount[1]>0)
		upAvgI2[0] = (upAvgI[0]+upAvgI[1])/2;
	else upAvgI2[0] = 0;

	for(z=1; z<sizez-1; z++)
	{
		if(upCount[z-1]>0 && upCount[z+1]>0 && upCount[z]>0) 
			upAvgI2[z] = (upAvgI[z-1]+upAvgI[z]+upAvgI[z+1])/3;
		else if(upCount[z-1]==0 && upCount[z+1]==0) upAvgI2[z]=0;
		else upAvgI2[z] = upAvgI[z];
	}
	if(upCount[sizez-1]>0 && upCount[sizez-2]>0)
		upAvgI2[sizez-1] = (upAvgI[sizez-1]+upAvgI[sizez-2])/2;
	else upAvgI2[sizez-1]=0;

	// assign globle variable
	for(z=0; z<sizez; z++)
	{
		if(normalSag[z].x<0) 
		{
			normalSag[z].x = -normalSag[z].x;
			normalSag[z].y = -normalSag[z].y;
		}
		spineNormalX[z] = Vec3(0, normalSag[z].x, normalSag[z].y);
	}

	int minP, minPi, minN, minNi, zz, closestMinP, closestMinN, downStreakP, downStreakN, downCount;

	// two methods to partition
	//	for low resolution data (pz>2), use the coronal view
	//	for high resolution data (pz<2), use the sagital view
	//
	pedicleValley.SetSize(0);
	if(pz>2)	//	low resolution
	{
		estimate_int-=1;
		// try to locate the valleys on coronal view
		for(z=1; z<sizez-1; z++)
		{
			minP=allAvgI[z-1];
			for(zz=z-1; zz>=0 && zz>=z-estimate_int; zz--) if(allAvgI[zz]<minP) minP=allAvgI[zz];
			minN=allAvgI[z+1];
			for(zz=z+1; zz<sizez && zz<=z+estimate_int; zz++) if(allAvgI[zz]<minN) minN=allAvgI[zz];

			// locate the valley
			if(allAvgI[z]<minP && allAvgI[z]<=minN)
			{
				pedicleValley.Add(z);
			}
		}
	}
	else {	// high resolution data

		// try to locate the disk on sag view

		for(z=1; z<sizez-1; z++)
		{
			if(upAvgI2[z]==0) continue;
			minP=upAvgI2[z-1];
			minPi = 1;
			downStreakP=0;
			closestMinP=0;
			downCount = 0;
			for(zz=z-1; zz>=0 && zz>=z-estimate_int; zz--) 
			{
				// count the streak of downward
				if(upAvgI2[zz]<=upAvgI2[zz+1]) downCount++;
				else downCount=0;

				if(upAvgI2[zz]!=0 && upAvgI2[zz]<minP) 
				{
					minP=upAvgI2[zz];
					minPi = z-zz;
					downStreakP = downCount;
				}

				if(closestMinP==0 && upAvgI2[zz]<upAvgI2[z]) closestMinP = z-zz;
			}

			minN=upAvgI2[z+1];
			minNi = 1;
			downStreakN=0;
			closestMinN=0;
			downCount = 0;
			for(zz=z+1; zz<sizez && zz<=z+estimate_int; zz++) 
			{
				// count the streak of downward
				if(upAvgI2[zz]<=upAvgI2[zz-1]) downCount++;
				else downCount=0;

				if(upAvgI2[zz]!=0 && upAvgI2[zz]<minN) 
				{
					minN=upAvgI2[zz];
					minNi = zz-z;
					downStreakN = downCount;
				}

				if(closestMinN==0 && upAvgI2[zz]<upAvgI2[z]) closestMinN = zz-z;
			}

			// locate the valley
			if((upAvgI2[z]<minP && upAvgI2[z]<=minN)
				|| 
			   ((upAvgI2[z]<minP || (closestMinP>=estimate_int*3/4 && minPi-downStreakP<=closestMinP)) && 
			   (upAvgI2[z]<=minN ||(closestMinN>=estimate_int*3/4 && minNi-downStreakN<=closestMinN))))
			{
				if(pedicleValley.GetSize()==0) pedicleValley.Add(z);
				else if(z-pedicleValley[pedicleValley.GetSize()-1]>estimate_int*3/4) pedicleValley.Add(z);
				else if(z-pedicleValley[pedicleValley.GetSize()-1]>estimate_int/2 &&  upAvgI2[z]<minP && upAvgI2[z]<=minN) 
					pedicleValley.Add(z);
				else if(upAvgI2[z]<upAvgI2[pedicleValley[pedicleValley.GetSize()-1]])
				{
					pedicleValley.RemoveAt(pedicleValley.GetSize()-1);
					pedicleValley.Add(z);
				}
			}
		}
	}

	double avgInterval, interval;
	double largestInterval, smallestInterval;
	avgInterval = 0;
	for(k=1; k<pedicleValley.GetSize(); k++)
	{
		avgInterval += pedicleValley[k]-pedicleValley[k-1];
	}
	// get rid of the largest and smallest interval
	largestInterval=0;
	smallestInterval=sizez;
	for(k=1; k<pedicleValley.GetSize(); k++)
	{
		interval = pedicleValley[k]-pedicleValley[k-1];
		if(interval>largestInterval)
			largestInterval=interval;
		if(interval<smallestInterval)
			smallestInterval=interval;
	}
	if(pedicleValley.GetSize()>4) avgInterval = (avgInterval-largestInterval-smallestInterval)/(double)(pedicleValley.GetSize()-3);
	else if(pedicleValley.GetSize()>1) avgInterval /= (double)(pedicleValley.GetSize()-1);
	avgInterval *= pz;

	// insert or adjust a valley if the interval is too big
	for(k=1; k<pedicleValley.GetSize(); k++)
	{
		float gap = (pedicleValley[k]-pedicleValley[k-1])*pz;
		if( (k>pedicleValley.GetSize()/2 && gap>= avgInterval*1.75) ||		// for L verte
			(k<=pedicleValley.GetSize()/2 && gap>= avgInterval*1.66) ||		// for T verte
			(avgInterval<estimate_gap*2.4 && gap>estimate_gap*2.4))
		{
			int adjustLevel=0;
			int upInterval=0, downInterval=0;
			if(k>1) upInterval=(pedicleValley[k]-pedicleValley[k-2])*pz;
			if(k<pedicleValley.GetSize()-1) downInterval=(pedicleValley[k+1]-pedicleValley[k-1])*pz;

			int tz;
			if((upInterval==0 || upInterval>avgInterval*2.5)
				&& (downInterval==0 || downInterval>avgInterval*2.5))
			{
				adjustLevel=0;
				tz =(pedicleValley[k]+pedicleValley[k-1])/2;
			}
			else if((upInterval>downInterval || upInterval==0) && (k>=1 && k<pedicleValley.GetSize()-1)) 
			{
				tz=(pedicleValley[k+1]+pedicleValley[k-1])/2;
				adjustLevel=1;
			}
			else 
			{
				tz=(pedicleValley[k]+pedicleValley[k-2])/2;
				adjustLevel=-1;
			}

			Vec2 cNormal, cPerp;
			int minZ=tz;
			// use the average for insertion, re-evaluate the score
			cNormal = normalSag[pedicleValley[k]]+normalSag[pedicleValley[k-1]];
			cNormal = -cNormal.normalize();
			cPerp = cNormal.perp().normalize();

			float mAccu=1e10;
			int mCount=0;
			for(z=tz-estimate_int/4; z<tz+estimate_int/4; z++)
			{
				if(spineCenter[z].x==0) continue;
				startPos.x = spineCenter[z].y;
				startPos.y = z*pz;
				int searchRange=(int)(upWidth[z]*1.5);
				int profileLength;
			
				float score = ComputePartitionAccumulateProfile(img3D, maskImg3D, segInfo, maskStruct, img2D_sag, img2D_sag_mask,
					spineCenter, spineWidthLeft, spineWidthRight, startPos, cNormal, cPerp, searchRange, estimate_gap,
					profileLength);

				if(score<mAccu && profileLength>0)
				{
					mAccu = score;
					mCount = profileLength;
					minZ = z;
				}
			}	// for z

		
			if(minZ>0)
			{
				if(adjustLevel==0) pedicleValley.InsertAt(k,minZ);
				else if(adjustLevel==1) pedicleValley[k]=minZ;
				else pedicleValley[k-1]=minZ;

				normalSag[minZ] = -cNormal;
			}
		}	// if gap
	}	// for k

	// check if I need to insert one valley at the very beginning
	if(pedicleValley.GetSize()>3)
	{
		if(pedicleValley[0]>avgInterval/pz+estimate_int/5)
		{
			int tz = pedicleValley[0]-(int)(avgInterval/pz);
			int minZ=tz;

			Vec2 cNormal, cPerp;
			// use the average for insertion, re-evaluate the score
			cNormal = normalSag[pedicleValley[0]];
			cNormal = -cNormal.normalize();
			cPerp = cNormal.perp().normalize();

			float mAccu=1e10;
			int mCount=0;
			for(z=tz-estimate_int/4; z<tz+estimate_int/4; z++)
			{
				if(spineCenter[z].x==0) continue;
				startPos.x = spineCenter[z].y;
				startPos.y = z*pz;
				int searchRange=(int)(upWidth[z]*1.5);
				int profileLength;
			
				float score = ComputePartitionAccumulateProfile(img3D, maskImg3D, segInfo, maskStruct, img2D_sag, img2D_sag_mask,
					spineCenter, spineWidthLeft, spineWidthRight, startPos, cNormal, cPerp, searchRange, estimate_gap,
					profileLength);

				if(score<mAccu && profileLength>0)
				{
					mAccu = score;
					mCount = profileLength;
					minZ = z;
				}
			}	// for z

/*			minP=upAvgI2[minZ];
			for(z=tz-estimate_int/4; z<tz+estimate_int/4; z++)
			{
				if(z<0 || z>=sizez || upAvgI2[z]==0) continue;

				if(upAvgI2[z]<minP)
				{
					minP=upAvgI2[z];
					minZ = z;
				}
			}	// for z
*/
			if(minZ>0)
			{
				pedicleValley.InsertAt(0,minZ);
				normalSag[minZ] = -cNormal;
			}
		}
	}	// pedicle

	// check if I need to insert one valley at the end
	if(pedicleValley.GetSize()>4)
	{
		int ps=pedicleValley.GetSize();
		if(pedicleValley[ps-1]<segInfo.bound2.z-avgInterval/pz-estimate_int/3)
		{
			int tz = pedicleValley[ps-1]+avgInterval/pz;
			int minZ=tz;
			
			Vec2 cNormal, cPerp;
			// use the average for insertion, re-evaluate the score
			cNormal = normalSag[pedicleValley[ps-1]];
			cNormal = -cNormal.normalize();
			cPerp = cNormal.perp().normalize();

			float mAccu=1e10;
			int mCount=0;
			for(z=tz-estimate_int/4; z<tz+estimate_int/4; z++)
			{
				if(spineCenter[z].x==0) continue;
				startPos.x = spineCenter[z].y;
				startPos.y = z*pz;
				int searchRange=(int)(upWidth[z]*1.5);
				int profileLength;
			
				float score = ComputePartitionAccumulateProfile(img3D, maskImg3D, segInfo, maskStruct, img2D_sag, img2D_sag_mask,
					spineCenter, spineWidthLeft, spineWidthRight, startPos, cNormal, cPerp, searchRange, estimate_gap,
					profileLength);

				if(score<mAccu && profileLength>0)
				{
					mAccu = score;
					mCount = profileLength;
					minZ = z;
				}
			}	// for z

/*			
			minP=upAvgI2[minZ];
			for(z=tz-estimate_int/4; z<tz+estimate_int/4; z++)
			{
				if(z<0 || z>=sizez || upAvgI2[z]==0) continue;

				if(upAvgI2[z]<minP)
				{
					minP=upAvgI2[z];
					minZ = z;
				}
			}	// for z
*/
			if(minZ>0)
			{
				pedicleValley.Add(minZ);
				normalSag[minZ] = -cNormal;
			}
		}
	}	// if pedicleValley

	// remove or adjust a valley if two valleys are too close

	// recompute the averageInterval
	avgInterval = 0;
	for(k=1; k<pedicleValley.GetSize(); k++)
	{
		avgInterval += pedicleValley[k]-pedicleValley[k-1];
	}
	// get rid of the largest and smallest interval
	largestInterval=0;
	smallestInterval=sizez;
	for(k=1; k<pedicleValley.GetSize(); k++)
	{
		interval = pedicleValley[k]-pedicleValley[k-1];
		if(interval>largestInterval)
			largestInterval=interval;
		if(interval<smallestInterval)
			smallestInterval=interval;
	}
	if(pedicleValley.GetSize()>4) avgInterval = (avgInterval-largestInterval-smallestInterval)/(double)(pedicleValley.GetSize()-3);
	else if(pedicleValley.GetSize()>1) avgInterval /= (double)(pedicleValley.GetSize()-1);
	avgInterval *= pz;

	// remove or adjust a valley if two valleys are too close
	for(k=1; k<pedicleValley.GetSize() && pedicleValley.GetSize()>4; k++)
	{
		float gap = (pedicleValley[k]-pedicleValley[k-1])*pz;
		bool tooClose=false;
		if(gap<avgInterval/2) 
			tooClose=true;	// the gap is too small
		else if(k>=2 && gap<(pedicleValley[k-1]-pedicleValley[k-2])*pz*0.66) 
			tooClose=true;		// the gap is less than half of its previous one
		else if(k<pedicleValley.GetSize()-1 && gap<(pedicleValley[k+1]-pedicleValley[k])*pz*0.66) 
			tooClose=true;	// the gap is less than half of the next one
		else if(k<pedicleValley.GetSize()-1 && fabs((pedicleValley[k+1]-pedicleValley[k-1])*pz-avgInterval)<fabs(gap-avgInterval)
			&& fabs((pedicleValley[k+1]-pedicleValley[k-1])*pz-avgInterval)<fabs((pedicleValley[k+1]-pedicleValley[k])*pz-avgInterval))
			tooClose=true;	// the sum of two adjacent interval is closer to avgInterval than each individual one

		if( tooClose)
		{
			if(k==1)
			{
				if((pedicleValley[k+1]-pedicleValley[k-1])*pz<avgInterval*1.5)
				{
					// keep the one with more evenly distribution
					if(fabs((pedicleValley[k+1]-pedicleValley[k])*pz-avgInterval)>fabs((pedicleValley[k+1]-pedicleValley[k-1])*pz-avgInterval)) 
						pedicleValley.RemoveAt(k);
					else pedicleValley.RemoveAt(k-1);

					// remove the higher one
///					if(upAvgI2[pedicleValley[k]]>upAvgI2[pedicleValley[k-1]]) pedicleValley.RemoveAt(k);
///					else pedicleValley.RemoveAt(k-1);
					k--;
				}
				else 
				{
					// adjust one of valleies (this case should be k)
					int tz=(pedicleValley[k+1]+pedicleValley[k-1])/2;
					int minZ=tz;
					minP=upAvgI2[minZ];
					for(z=tz-estimate_int/3; z<tz+estimate_int/3; z++)
					{
						if(z<0 || z>=sizez || upAvgI2[z]==0) continue;
						if(upAvgI2[z]<minP)
						{
							minP=upAvgI2[z];
							minZ = z;
						}
					}	// for z
					if(minZ>0)
					{
						pedicleValley[k]=minZ;
					}
				}
			}	// if k
			else if(k==pedicleValley.GetSize()-1)
			{
				if((pedicleValley[k]-pedicleValley[k-2])*pz<avgInterval*1.5)
				{
					// keep the one with more evenly distribution
					if(fabs((pedicleValley[k]-pedicleValley[k-2])*pz-avgInterval)>fabs((pedicleValley[k-1]-pedicleValley[k-2])*pz-avgInterval)) 
						pedicleValley.RemoveAt(k);
					else pedicleValley.RemoveAt(k-1);

					// remove the higher one
///					if(upAvgI2[pedicleValley[k]]>upAvgI2[pedicleValley[k-1]]) pedicleValley.RemoveAt(k);
///					else pedicleValley.RemoveAt(k-1);
					k--;
				}
				else 
				{
					// adjust one of valleies (this case should be k-1)
					int tz=(pedicleValley[k]+pedicleValley[k-2])/2;
					int minZ=tz;
					minP=upAvgI2[minZ];
					for(z=tz-estimate_int/3; z<tz+estimate_int/3; z++)
					{
						if(z<0 || z>=sizez || upAvgI2[z]==0) continue;
						if(upAvgI2[z]<minP)
						{
							minP=upAvgI2[z];
							minZ = z;
						}
					}	// for z
					if(minZ>0)
					{
						pedicleValley[k-1]=minZ;
					}
				}
			}	// else if k
			else // all other cases
			{
				if((pedicleValley[k]-pedicleValley[k-2])*pz<avgInterval*1.5 && (pedicleValley[k+1]-pedicleValley[k-1])*pz<avgInterval*1.5)
				{
					// keep the one with more evenly distribution
					if((pedicleValley[k+1]-pedicleValley[k])<(pedicleValley[k-1]-pedicleValley[k-2])) pedicleValley.RemoveAt(k);
					else pedicleValley.RemoveAt(k-1);
					k--;
				}
				else 
				{
					if((pedicleValley[k-1]-pedicleValley[k-2])>pedicleValley[k+1]-pedicleValley[k])
					{
						// adjust one of valleies (this case should be k-1)
						int tz=(pedicleValley[k]+pedicleValley[k-2])/2;
						int minZ=tz;
						minP=upAvgI2[minZ];
						for(z=tz-estimate_int/3; z<tz+estimate_int/3; z++)
						{
							if(z<0 || z>=sizez || upAvgI2[z]==0) continue;
							if(upAvgI2[z]<minP)
							{
								minP=upAvgI2[z];
								minZ = z;
							}
						}	// for z
						if(minZ>0)
						{
							pedicleValley[k-1]=minZ;
						}
					}
					else
///					if((pedicleValley[k+1]-pedicleValley[k-1])*pz>avgInterval*2)
					{
						// adjust one of valleies (this case should be k)
						int tz=(pedicleValley[k+1]+pedicleValley[k-1])/2;
						int minZ=tz;
						minP=upAvgI2[minZ];
						for(z=tz-estimate_int/3; z<tz+estimate_int/3; z++)
						{
							if(z<0 || z>=sizez || upAvgI2[z]==0) continue;
							if(upAvgI2[z]<minP)
							{
								minP=upAvgI2[z];
								minZ = z;
							}
						}	// for z
						if(minZ>0)
						{
							pedicleValley[k]=minZ;
						}
					}

				}
			}	// else if k

		}	// if gap
	}	// for k

	// recompute the averageInterval
	avgInterval = 0;
	for(k=1; k<pedicleValley.GetSize(); k++)
	{
		avgInterval += pedicleValley[k]-pedicleValley[k-1];
	}
	// get rid of the largest and smallest interval
	largestInterval=0;
	smallestInterval=sizez;
	for(k=1; k<pedicleValley.GetSize(); k++)
	{
		interval = pedicleValley[k]-pedicleValley[k-1];
		if(interval>largestInterval)
			largestInterval=interval;
		if(interval<smallestInterval)
			smallestInterval=interval;
	}
	if(pedicleValley.GetSize()>4) avgInterval = (avgInterval-largestInterval-smallestInterval)/(double)(pedicleValley.GetSize()-3);
	else if(pedicleValley.GetSize()>1) avgInterval /= (double)(pedicleValley.GetSize()-1);
	avgInterval *= pz;

	// check if first and last valley are off mark (angles are too different)
	if(pedicleValley.GetSize()>3 && (pedicleValley[1]-segInfo.bound1.z)*pz<avgInterval*1.2 && (pedicleValley[1]-pedicleValley[0])*pz<avgInterval)
	{
		double a0 = atan2(normalSag[pedicleValley[0]].y, normalSag[pedicleValley[0]].x);
		double a1 = atan2(normalSag[pedicleValley[1]].y, normalSag[pedicleValley[1]].x);
		double a2 = atan2(normalSag[pedicleValley[2]].y, normalSag[pedicleValley[2]].x);
		if(a2<a1)
		{
			if(a0-a1>15*3.14159/180 || a0-a1<-5*3.14159/180) pedicleValley.RemoveAt(0);
		}
		else
		{
			if(a0-a1<-15*3.14159/180 || a0-a1>5*3.14159/180) pedicleValley.RemoveAt(0);
		}
	}

	int psize = pedicleValley.GetSize();
	if(psize>3 && (segInfo.bound2.z-pedicleValley[psize-2])*pz<avgInterval*1.2 && (pedicleValley[psize-1]-pedicleValley[psize-2])*pz<avgInterval)
	{
		double a0 = atan2(normalSag[pedicleValley[psize-1]].y, normalSag[pedicleValley[psize-1]].x);
		double a1 = atan2(normalSag[pedicleValley[psize-2]].y, normalSag[pedicleValley[psize-2]].x);
		double a2 = atan2(normalSag[pedicleValley[psize-3]].y, normalSag[pedicleValley[psize-3]].x);
		if(a2<a1)
		{
			if(a0-a1>15*3.14159/180 || a0-a1<-5*3.14159/180) pedicleValley.RemoveAt(psize-1);
		}
		else
		{
			if(a0-a1<-15*3.14159/180 || a0-a1>5*3.14159/180) pedicleValley.RemoveAt(psize-1);
		}
	}

	// re-interpolate the normalSag using piecewise b-spline
	intDynArray isKey;
	isKey.SetSize(sizez);

	for(z=0; z<sizez; z++)
	{
		isKey[z]=0;
		xt[z] = z;
		yt[z] = -1;
	}

	for(k=0; k<=pedicleValley.GetSize(); k++)
	{
		if(k<pedicleValley.GetSize()) isKey[pedicleValley[k]]=1;

		double midAngle, angleRange1, angleRange2, angle;
		int z1, z2;

		if(k==0)
		{
			z1 = 0;
			z2 = pedicleValley[k];
			midAngle = atan2(normalSag[z2].y, normalSag[z2].x);
			angleRange2 = atan2(normalSag[pedicleValley[1]].y, normalSag[pedicleValley[1]].x);
			angleRange1 = 2*midAngle-angleRange2;
			if(angleRange1>angleRange2)
			{
				angle = angleRange1; angleRange1=angleRange2; angleRange2=angle;
			}
		}
		else if(k==pedicleValley.GetSize())
		{
			z1 = pedicleValley[k-1];
			z2 = sizez-1;
			midAngle = atan2(normalSag[z1].y, normalSag[z1].x);
			angleRange2 = atan2(normalSag[pedicleValley[k-2]].y, normalSag[pedicleValley[k-2]].x);
			angleRange1 = 2*midAngle-angleRange2;
			if(angleRange1>angleRange2)
			{
				angle = angleRange1; angleRange1=angleRange2; angleRange2=angle;
			}
		}
		else 
		{
			z1 = pedicleValley[k-1];
			z2 = pedicleValley[k];
			angleRange2 = atan2(normalSag[z2].y, normalSag[z2].x);
			angleRange1 = atan2(normalSag[z1].y, normalSag[z1].x);
			if(angleRange1>angleRange2)
			{
				angle = angleRange1; angleRange1=angleRange2; angleRange2=angle;
			}
		}

		for(z=z1; z<z2; z++)
		{
			if(normalSag[z].x!=0 || normalSag[z].y!=0)
			{
				angle = atan2(normalSag[z].y, normalSag[z].x);
				if(angle>=angleRange1 && angle<=angleRange2) yt[z]=angle;
			}
		}
	}	// for k

	PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
	for(z=0; z<spinalCord3D.GetSize(); z++)
	{
		if((normalSag[z].x!=0 || normalSag[z].y!=0) && isKey[z]==0)
		{
			normalSag[z].x = cos(yt2[z]);
			normalSag[z].y = sin(yt2[z]);
		}
	}

	for(z=0; z<sizez; z++)
	{
		if(normalSag[z].x<0) 
		{
			normalSag[z].x = -normalSag[z].x;
			normalSag[z].y = -normalSag[z].y;
		}
		spineNormalX[z] = Vec3(0, normalSag[z].x, normalSag[z].y);
	}


	// export the profile for debugging purpose
	if(debugMode)
	{
		FILE *fp;
		fp = fopen("c:\\tmp\\integral.txt", "w");
		fprintf(fp, "leftCount leftAvg leftAvg1 leftAvg2 rightCount rightAvg rightAvg1 rightAvg2 allCount allAvg allAvg2 upCount upAvg upAvg2 upAvgm\n");
		for(z=0;z<sizez; z++)
		{
			fprintf(fp, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", 
				leftCount[z], leftAvgI[z], leftAvgI1[z], leftAvgI2[z], 
				rightCount[z], rightAvgI[z], rightAvgI1[z], rightAvgI2[z],
				allCount[z], allAvgI[z], allAvgI2[z], upCount[z], upAvgI[z], upAvgI2[z], upAvgIm[z]);
		}
		// print out the pedicles
		for(z=0; z<pedicleValley.GetSize(); z++)
		{
			fprintf(fp, "%d\n", pedicleValley[z]);
		}
		fclose(fp);

	}

	// in low resolution case
	// smooth the normals at pedicle using piecewise b-spline
	if(pz>2)
	{
		for(z=0; z<sizez; z++)
		{
			xt[z] = z;
			yt[z] = -1;
		}

		for(k=0; k<pedicleValley.GetSize(); k++)
		{
			z = pedicleValley[k];
			if(normalSag[z].y!=0 || normalSag[z].x!=0)
				yt[z] = atan2(normalSag[z].y, normalSag[z].x);
			else yt[z]=-1;
		}
		PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
		for(k=0; k<pedicleValley.GetSize(); k++)
		{
			z = pedicleValley[k];
			normalSag[z].x = cos(yt2[z]);
			normalSag[z].y = sin(yt2[z]);
			spineNormalX[z].y = normalSag[z].x;
			spineNormalX[z].z = normalSag[z].y;
		}
	}

	// compute the partition of spinous process (back of vertebra)
	// only for high resolution cases

	if(pz<2)
	{
		// compute the intensity integral along the normal direction on sagital view
		// also adjust the normal to get best cut on the sagital view
		vec2DynArray normalSagBack;
		normalSagBack.SetSize(sizez);
		intDynArray backAvgI, backCount, backAvgI2;
		backAvgI.SetSize(sizez);
		backAvgI2.SetSize(sizez);
		backCount.SetSize(sizez);

		for(z=0; z<sizez; z++)
		{
			normalSagBack[z] = normalSag[z];
			backAvgI[z] = backAvgI2[z] = 0;
			backCount[z] = 0;
			if(spineCenter[z].x==0) continue;

			// back side
			count0 = count1 = count2 = 0;
			startPos.x = spineCenter[z].y;
			startPos.y = z*pz;

			curPos = startPos;

			int searchRange=(int)(upWidth[z]);
			// adjust the angle within 15 degree with 2 degree step
			float midAngle = atan2(normalSag[z].y, normalSag[z].x)+3.14159/4;
			float angleSearchRange = 15*3.14159/180;
			float angleStep = 1*3.14159/180;
			float angleRange1 = midAngle-angleSearchRange;
			float angleRange2 = midAngle+angleSearchRange;

			float cAngle, mAngle;
			Vec2 cNormal;
			float mAccu=1e10, cAccu;
			int mCount=0;

			for(cAngle=angleRange1; cAngle<=angleRange2; cAngle+=angleStep)
			{
				cNormal = Vec2(cos(cAngle), sin(cAngle));
				cAccu = 0;
				count0 = 0;
				count1 = 0;
				curPos = startPos;
				for(k=0; k<searchRange; k++)
				{
					curPos += cNormal;
					posx = curPos.x/py;
					posy = curPos.y/pz;
					x = (int)(posx+0.5);
					y = (int)(posy+0.5);

					mask = img2D_sag_mask->FastGet(x, y);
			
					if(mask==maskStruct.mask_body || mask==maskStruct.mask_corticalBone || mask==maskStruct.mask_spongyBone)
					{
						gray = (short)img2D_sag->GetPixelAtPosition(posx, posy);
						if(gray<=0) break;
						count0++;
						cAccu += gray;
					}
				}
				if(count0>20)
				{
					cAccu /= count0;
					if(cAccu<mAccu)
					{
						mAccu = cAccu;
						mCount = count0;
						mAngle= cAngle;
					}
				}
			}	// for cAngle

			if(mCount>0)
			{
				backCount[z] = mCount;
				backAvgI[z] = mAccu;
				normalSagBack[z].x = cos(mAngle);
				normalSagBack[z].y = sin(mAngle);
			}
		}	// for z

		// assign globle variable
		spineNormalBack.SetSize(sizez);
		for(z=0; z<sizez; z++)
		{
			spineNormalBack[z] = Vec3(0, normalSagBack[z].x, normalSagBack[z].y);
		}

		int stz, edz, minz, lastz=-1;
		backValley.SetSize(0);
		for(zz=-1; zz<pedicleValley.GetSize()-1; zz++)
		{
			if(zz==-1) 
			{
				stz=0;
				edz=pedicleValley[0];
			}
			else 
			{
				stz=pedicleValley[zz];
				edz = pedicleValley[zz+1];

				if(lastz>0 && stz<lastz+(pedicleValley[zz+1]-pedicleValley[zz])/2) stz=lastz+(pedicleValley[zz+1]-pedicleValley[zz])/2;
				if(lastz>0 && edz>lastz+(pedicleValley[zz+1]-pedicleValley[zz])*3/2) edz=lastz+(pedicleValley[zz+1]-pedicleValley[zz])*3/2;
			}

			minP = 10000;
			minz = -1;
			for(z=stz; z<edz; z++)
			{
				if(backAvgI[z]>0 && backAvgI[z]<minP)
				{
					minP = backAvgI[z];
					minz = z;
				}
			}
			if(minz!=-1) 
			{
				backValley.Add(minz);
				lastz = minz;
			}
		}	// for zz
	}	// if pz

	return CIS_OK;
}

float ComputePartitionAccumulateProfile(CIS_Array_Image3D_short *img3D,
									    CIS_Array_Image3D_short *maskImg3D,
									    SpineSegmentationInfo &segInfo, 
										SpineSegmentationMaskStructure &maskStruct,
									    CIS_Array_Image2D_short *img2D_sag,
									    CIS_Array_Image2D_short *img2D_sag_mask,
										vec3DynArray &spineCenter, doubleDynArray &spineWidthLeft, doubleDynArray &spineWidthRight,
										Vec2 &startPos, Vec2 &cNormal, Vec2 &cPerp, 
										int searchRange, float estimate_gap,
										int &profileLength)
{
	float score = 100000;
	profileLength = 0;

	float cAccu = 0;
	float cAccum = 0;
	int count0 = 0;
	int count1 = 0;
	Vec2 curPos = startPos;
	float posx, posy;
	int x, y, x3, k2, k4;
	float px, py, pz;
	short mask, gray;
	int gray1;
	short *maskA, *imgA, *maskA_sag;
	imgA = img3D->GetArray();
	maskA = maskImg3D->GetArray();
	maskA_sag = img2D_sag_mask->GetArray();
	px = img3D->Get_Pixel_SizeX();
	py = img3D->Get_Pixel_SizeY();
	pz = img3D->Get_Pixel_SizeZ();
	int sizex = img3D->Num_Cols();
	int sizey = img3D->Num_Rows();
	int sizez = img3D->Num_Levels();
	int sizexy = sizex*sizey;

	for(int j=0; j<searchRange; j++)
	{
		curPos += cNormal;
		posx = curPos.x/py;
		posy = curPos.y/pz;
		x = (int)(posx+0.5);
		y = (int)(posy+0.5);
		if(y<0 || y>segInfo.bound2.z-2) break;
		x3 = (int)(spineCenter[y].x/px+0.5);
		if(x3<=0) break;

		k2 = y*sizey+x;
		k4 = x3+y*sizexy+x*sizex;
		mask = maskA_sag[k2];
		
		if(mask==maskStruct.mask_corticalBone || mask==maskStruct.mask_spongyBone)
		{
			gray = (short)img2D_sag->GetPixelAtPosition(posx, posy);
			if(gray<=100) break;

			gray1 = imgA[k4];
			if(gray1<=0) break;
			int count=0;
			int k5;
			for(k5=1;k5<spineWidthLeft[y]/px/2; k5++)
			{
				if(maskA[k4-k5]==maskStruct.mask_corticalBone || maskA[k4-k5]==maskStruct.mask_spongyBone)
				{
					gray1 += imgA[k4-k5];
					count++;
				}
			}
			for(k5=1;k5<spineWidthRight[y]/px/2; k5++)
			{
				if(maskA[k4+k5]==maskStruct.mask_corticalBone || maskA[k4+k5]==maskStruct.mask_spongyBone)
				{
					gray1 += imgA[k4+k5];
					count++;
				}
			}

			gray1 /= (count+1);

			gray = gray1;
			count1++;
			cAccum += gray;
			if(count1<=7) continue;	// skip the first 7 pixels, which somtimes could be big for disk
			count0++;
			cAccu += gray;
		}
	}	// for j

	if(count0<=15) return score;

	cAccu /= count0;
	cAccum /= count1;

	// check both sides to locate the border of vert body
	float d, best_ud, best_dd;
	float best_aAccu_u, best_aAccu_d;
	// find the maximum in upward direction
	best_ud = -1;
	best_aAccu_u=0;
	for(d=pz; d<estimate_gap/3; d+=pz)
	{
		curPos = startPos+cPerp*d;
		float cAccu_u;
		int count0_u;
		cAccu_u = 0;
		count0_u = 0;
		for(int j=0; j<searchRange; j++)
		{
			curPos += cNormal;
			posx = curPos.x/py;
			posy = curPos.y/pz;
			x = (int)(posx+0.5);
			y = (int)(posy+0.5);
			if(y<0 || y>segInfo.bound2.z-2) break;

			k2 = y*sizey+x;
			mask = maskA_sag[k2];
		
			if(mask==maskStruct.mask_corticalBone || mask==maskStruct.mask_spongyBone)
			{
				gray = (short)img2D_sag->GetPixelAtPosition(posx, posy);
				if(gray<=100) break;
				count0_u++;
				cAccu_u += gray;
			}
		}

		if(count0_u>22)
		{
			cAccu_u /= count0_u;
			if(cAccu_u>best_aAccu_u)
			{
				best_aAccu_u=cAccu_u;
				best_ud = d;
			}
		}
	}	// for d

	if(best_ud==-1 || best_aAccu_u<cAccu+10) return score;

	// find the maximum in downward direction
	best_dd = -1;
	best_aAccu_d=0;
	for(d=pz; d<estimate_gap/3; d+=pz)
	{
		curPos = startPos-cPerp*d;
		float cAccu_d;
		int count0_d;
		cAccu_d = 0;
		count0_d = 0;
		for(int j=0; j<searchRange; j++)
		{
			curPos += cNormal;
			posx = curPos.x/py;
			posy = curPos.y/pz;
			x = (int)(posx+0.5);
			y = (int)(posy+0.5);
			if(y<0 || y>segInfo.bound2.z-2) break;
	
			k2 = y*sizey+x;
			mask = maskA_sag[k2];
		
			if(mask==maskStruct.mask_corticalBone || mask==maskStruct.mask_spongyBone)
			{
				gray = (short)img2D_sag->GetPixelAtPosition(posx, posy);
				if(gray<=100) break;
				count0_d++;
				cAccu_d += gray;
			}
		}

		if(count0_d>22)
		{
			cAccu_d /= count0_d;
			if(cAccu_d>best_aAccu_d)
			{
				best_aAccu_d=cAccu_d;
				best_dd = d;
			}
		}
	}	// for d

	if(best_dd==-1 || best_aAccu_d<cAccu+10) return score;

	score = cAccu-(best_aAccu_u+best_aAccu_d)/2+fabs(best_aAccu_u-best_aAccu_d)/4;
	profileLength = count0;

	return score;
}

// segment the vertebral disk based on the spine partition
//
int NIH_SpineSegmentation_VertebralDiskDetection(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   VertebraStruct2D *vertebraTemplate, 
									   vec3DynArray &spinalCord3D,
									   vec3DynArray &spineCenter, vec3DynArray &spineNormalX, vec3DynArray &spineNormalY, vec3DynArray &spineNormalBack, 
									   doubleDynArray &spineWidthLeft, doubleDynArray &spineWidthRight, doubleDynArray &spineWidthUp, doubleDynArray &spineWidthDown, 
									   intDynArray &pedicleValley, intDynArray &backValley, 
									   CIS_Array_Image2D_short *img2D_sag,
									   CIS_Array_Image2D_short *img2D_cor,
									   CIS_Array_Image2D_short *img2D_sag_mask,
									   CIS_Array_Image2D_short *img2D_cor_mask,
									   bool debugMode)
{
	if(spinalCord3D.GetSize()<10) return CIS_ERROR;
	if(pedicleValley.GetSize()<1) return CIS_ERROR;

	int sizex, sizey, sizez, sizexy;
	short *maskA, *imgA;
	int x, y, z, k, k2;
	float px, py, pz;

	sizex = maskImg3D->Num_Cols();
	sizey = maskImg3D->Num_Rows();
	sizez = maskImg3D->Num_Levels();
	px = img3D->Get_Pixel_SizeX();
	py = img3D->Get_Pixel_SizeY();
	pz = img3D->Get_Pixel_SizeZ();
	sizexy = sizex*sizey;

	maskA = maskImg3D->GetArray();
	imgA = img3D->GetArray();

	double diskThick = 10;

	int p, ck, cck;
	Vec3 startP, scanV, cV, ccV, dirX, dirY, normal;
	double stepSize, sx, sy, _sizeX, _sizeY, r;
	IntVec3 cVi, ccVi, bb1, bb2;
	CIS_Array_Image3D_short *binImg=NULL, *blobImg=NULL;

	for(p=0; p<pedicleValley.GetSize(); p++)
	{
		z = pedicleValley[p];
		dirY = spineNormalX[z];
		dirX = spineNormalY[z];
		
		if(dirY.len()==0 || dirX.len()==0) continue;

		normal = dirX % dirY;
		dirX = dirY % normal;
		dirY = normal % dirX;
		_sizeX = (spineWidthLeft[z]+spineWidthRight[z])*1.5;
		_sizeY = (spineWidthUp[z])*1.5;

		if(dirY.y>0)
		{
			dirY = -dirY;
		}
		startP = Vec3(segInfo.cordCenter[z].x*px, segInfo.cordCenter[z].y*py, z*pz);
		stepSize = px/2;
		bb1 = IntVec3(sizex, sizey, sizez);
		bb2 = IntVec3(0,0,0);

		scanV = startP;
		for(sy=0; sy<_sizeY; sy+=stepSize, scanV += (dirY*stepSize))
		{
			for(int dir=-1; dir<=1; dir+=2)
			{
				for(sx=0; sx<_sizeX/2; sx+=stepSize)
				{
					cV = scanV+dirX*(sx*dir);
					cVi = IntVec3(cV.x/px, cV.y/py, cV.z/pz);
					ck = cVi.x+cVi.y*sizex+cVi.z*sizexy;
					if(cVi.z<1) continue;

					if(maskA[ck]==maskStruct.mask_corticalBone || maskA[ck]==maskStruct.mask_spongyBone)
					{
						for(r=-diskThick; r<=diskThick; r+=stepSize)
						{
							ccV = cV+normal*r;
							ccVi = IntVec3(ccV.x/px, ccV.y/py, ccV.z/pz);
							if(ccVi.z<1) continue;
							cck = ccVi.x+ccVi.y*sizex+ccVi.z*sizexy;
							if(maskA[cck]==maskStruct.mask_spongyBone)
							{
								maskA[cck]=maskStruct.mask_vertebralDisk;

								if(ccVi.x<bb1.x) bb1.x=ccVi.x;
								if(ccVi.y<bb1.y) bb1.y=ccVi.y;
								if(ccVi.z<bb1.z) bb1.z=ccVi.z;
								if(ccVi.x>bb2.x) bb2.x=ccVi.x;
								if(ccVi.y>bb2.y) bb2.y=ccVi.y;
								if(ccVi.z>bb2.z) bb2.z=ccVi.z;
							}
						}	// for r
					}	// if maskA
				}	// for sx
			}	// for dir
		}

		if(bb1.z>bb2.z) continue;

		// perform blob analysis to remove artefacts
		//
		bb1.x-=1; bb1.y-=1; bb1.z-=1;
		bb2.x+=1; bb2.y+=1; bb2.z+=1;
		maskImg3D->SubImage(binImg, bb1.x, bb1.y, bb1.z, bb2.x, bb2.y, bb2.z);
		maskImg3D->SubImage(blobImg, bb1.x, bb1.y, bb1.z, bb2.x, bb2.y, bb2.z);
		int sub_sizex, sub_sizey, sub_sizez, sub_sizexyz, sub_sizexy;
		short *blobA, *binA;
		sub_sizex = binImg->Num_Cols();
		sub_sizey = binImg->Num_Rows();
		sub_sizez = binImg->Num_Levels();
		sub_sizexy = sub_sizex*sub_sizey;
		sub_sizexyz= sub_sizex*sub_sizey*sub_sizez;
		blobA = blobImg->GetArray();
		binA = binImg->GetArray();
		for(k=0; k<sub_sizexyz; k++)
		{
			if(blobA[k]==maskStruct.mask_vertebralDisk) binA[k]=1;
			else binA[k]=0;
		}

		intDynArray blobRank;
		intVec3DynArray blobCentroid;

		NIH_Algo_Blob_Labelling_3D(binImg, blobImg, blobRank, blobCentroid, false);

		int largestBlob=-1, largestSize=0;
		for(k=0; k<blobRank.GetSize(); k++)
		{
			if(blobRank[k]>largestSize)
			{
				largestSize = blobRank[k];
				largestBlob=k+1;
			}
		}

		for(z=0, k=0; z<sub_sizez; z++) for(y=0; y<sub_sizey; y++) for(x=0; x<sub_sizex; x++, k++)
		{
			if(binA[k]==1 && blobA[k]!=largestBlob) 
			{
				maskA[(z+bb1.z)*sizexy+(y+bb1.y)*sizex+x+bb1.x]=maskStruct.mask_spongyBone;
				binA[k]=0;
			}
		}

		// 2D analysis on sagital view to eliminate false disk
		//
		CIS_Array_Image2D_short *binImg2D=NULL, *blobImg2D=NULL;
		binImg2D = new CIS_Array_Image2D_short(sub_sizey, sub_sizez);
		blobImg2D = new CIS_Array_Image2D_short(sub_sizey, sub_sizez);
		intDynArray blobRank2D;
		intVec2DynArray blobCentroid2D;
		short *binA2D, *blobA2D;
		int count;
		binA2D = binImg2D->GetArray();
		blobA2D = blobImg2D->GetArray();

		for(x=0; x<sub_sizex; x++)
		{
			binImg->GetSliceImage(*binImg2D, x, 0, true);
			count=0;
			for(k=0; k<sub_sizey*sub_sizez; k++)
				if(binA2D[k]==1) count++;
			if(count==0) continue;
			NIH_Algo_Blob_Labelling(binImg2D, blobImg2D, blobRank2D, blobCentroid2D, false);
			largestBlob=-1, largestSize=0;
			for(k=0; k<blobRank2D.GetSize(); k++)
			{
				if(blobRank2D[k]>largestSize)
				{
					largestSize = blobRank2D[k];
					largestBlob=k+1;
				}
			}

			for(z=0, k=0; z<sub_sizez; z++) for(y=0; y<sub_sizey; y++, k++) 
			{
				if(binA2D[k]==1 && blobA2D[k]!=largestBlob) 
				{
					maskA[(z+bb1.z)*sizexy+(y+bb1.y)*sizex+x+bb1.x]=maskStruct.mask_spongyBone;
					binA2D[k]=0;
					binA[z*sub_sizexy+y*sub_sizex+x]=0;
				}
			}
		}

		// further scan to get rid of outliers
		scanV = startP;
		double cstart, cend, lstart, lend;
		int numSeg;
		for(sy=0; sy<_sizeY; sy+=stepSize, scanV += (dirY*stepSize))
		{
			for(int dir=-1; dir<=1; dir+=2)
			{
				for(sx=0; sx<_sizeX/2; sx+=stepSize)
				{
					cV = scanV+dirX*(sx*dir);
					cVi = IntVec3(cV.x/px, cV.y/py, cV.z/pz);
					ck = cVi.x+cVi.y*sizex+cVi.z*sizexy;
					if(cVi.z<1) continue;

					if(maskA[ck]==maskStruct.mask_corticalBone || maskA[ck]==maskStruct.mask_spongyBone || maskA[ck]==maskStruct.mask_vertebralDisk)
					{
						cstart = cend = -1000;
						lstart = lend = -1000;
						numSeg = 0;
						for(r=-diskThick; r<=diskThick; r+=stepSize)
						{
							ccV = cV+normal*r;
							ccVi = IntVec3(ccV.x/px, ccV.y/py, ccV.z/pz);
							cck = ccVi.x+ccVi.y*sizex+ccVi.z*sizexy;
							if(ccVi.z<1) continue;
							if(maskA[cck]==maskStruct.mask_vertebralDisk)
							{
								if(cstart>-1000) cend=r;
								else {cstart=r; cend=r; numSeg++;};
							} else {
								if(cend-cstart>lend-lstart)
								{
									lstart = cstart;
									lend = cend;
								}
								cstart = cend = -1000;
							}
						}	// for r

						// only keep the longest segment
						if(numSeg>1 && lstart>-1000)
						{
							for(r=-diskThick; r<=diskThick; r+=stepSize)
							{
								if(r>=lstart && r<=lend) continue;
								ccV = cV+normal*r;
								ccVi = IntVec3(ccV.x/px, ccV.y/py, ccV.z/pz);
								cck = ccVi.x+ccVi.y*sizex+ccVi.z*sizexy;
								if(ccVi.z<bb1.z) continue;
								if(maskA[cck]==maskStruct.mask_vertebralDisk)
								{
									binA[ccVi.x-bb1.x+(ccVi.y-bb1.y)*sub_sizex+(ccVi.z-bb1.z)*sub_sizex*sub_sizey]=0;
									maskA[cck]=maskStruct.mask_spongyBone;
								} 
							}	// for r
						}
					}	// if maskA
				}	// for sx
			}	// for dir
		}

		// do a dilation to fill gaps
		CIS_IPA_3D_Dilate(binImg, 2, 10);
		for(z=0, k=0; z<sub_sizez; z++) for(y=0; y<sub_sizey; y++) for(x=0; x<sub_sizex; x++, k++)
		{
			if(binA[k]==1) 
			{
				ck = (z+bb1.z)*sizexy+(y+bb1.y)*sizex+x+bb1.x;
				if(maskA[ck]==maskStruct.mask_spongyBone) maskA[ck]=maskStruct.mask_vertebralDisk;
			}
		}

		CIS_IPA_3D_Dilate(binImg, 1, 10);
		CIS_IPA_3D_Erode(binImg, 1, 10);
		for(z=0, k=0; z<sub_sizez; z++) for(y=0; y<sub_sizey; y++) for(x=0; x<sub_sizex; x++, k++)
		{
			if(binA[k]==1) 
			{
				ck = (z+bb1.z)*sizexy+(y+bb1.y)*sizex+x+bb1.x;
				if(maskA[ck]==maskStruct.mask_corticalBone) maskA[ck]=maskStruct.mask_vertebralDisk;
			}
		}

		delete binImg; binImg=NULL;
		delete blobImg; blobImg=NULL;

		delete binImg2D; binImg2D=NULL;
		delete blobImg2D; blobImg2D=NULL;

	}	// for p

	intDynArray filledZ;
	filledZ.SetSize(sizez);
	for(z=0; z<sizez; z++) 
	{
		filledZ[z]=-1;
	}
		
	for(k=0; k<spinalCord3D.GetSize(); k++)
	{
		z = img3D->GetSliceLevel(spinalCord3D[k].z,1);

		if(filledZ[z]==1 && filledZ[z-1]==-1) z=z-1;
		while(z<sizez-1 && filledZ[z]!=-1) z++;
		filledZ[z] = k;
	}

	int k3;
	short gray;
	short *sagMaskA, *corMaskA, *sagA;
	sagMaskA = img2D_sag_mask->GetArray();
	corMaskA = img2D_cor_mask->GetArray();
	sagA = img2D_sag->GetArray();

	// reset the curved reformation
	for(k=0; k<spinalCord3D.GetSize(); k++)
	{
		// sagittal reformation
		x = (int)(spinalCord3D[k].x/px);
		z = img3D->GetSliceLevel(spinalCord3D[k].z,1);

		if(z<0 || z>sizez-1) continue;
		// filledZ array is used to handle duplicate slices
		//
//		if(filledZ[z]==1 && filledZ[z-1]==0) z=z-1;
//		while(z<sizez-1 && filledZ[z]==1) z++;
//		filledZ[z] = k;
		
		// saggital reformation
		k2 = z*sizexy+x;
		k3 = z*sizey;
		for(y=0; y<img2D_sag->Num_Cols(); y++, k2+=sizex, k3++)
		{
			gray = maskA[k2];
			if(gray==maskStruct.mask_vertebralDisk) sagMaskA[k3] = gray;
		}

		// coronal reformation
		y = (int)(spinalCord3D[k].y/py);

		k2 = z*sizexy+y*sizex;
		k3 = z*sizey;
		for(x=0; x<img2D_cor->Num_Cols(); x++, k2++, k3++)
		{
			gray = maskA[k2];
			if(gray==maskStruct.mask_vertebralDisk) corMaskA[k3] = gray;
		}
	}

	// refine the disk detection on curved sag reformation
	Vec2 normalX;
	for(p=0; p<pedicleValley.GetSize(); p++)
	{
		z = pedicleValley[p];
		k = filledZ[z];
		if(k==-1) continue;

		y = (int)(spinalCord3D[k].y/py);
		normalX.x = spineNormalX[z].y;
		normalX.y = spineNormalX[z].z;
		if(normalX.x>0) normalX=-normalX;
		_sizeY = (spineWidthUp[z])*1.5;

		Vec2 end1, end2, bb1, bb2;
		IntVec2 bbi1, bbi2;

		end1 = Vec2(y*py, z*pz);
		end2 = end1+normalX*_sizeY;
		bb1.x=end2.x; bb1.y=end2.y-diskThick*1.5;
		// if end1.y < end2.y
		if(end1.y-diskThick*1.5<bb1.y) bb1.y=end1.y-diskThick*1.5;
		bb2.x=end1.x; bb2.y=end1.y+diskThick*1.5;
		// if end2.y > end1.y
		if(end2.y+diskThick*1.5>bb2.y) bb2.y=end2.y+diskThick*1.5;
		
		if(bb1.y<1) bb1.y=1;
		if(bb2.y>sizez-1) bb2.y=sizez-1;
		if(bb1.y>sizez-1) bb1.y=sizez-1; //is this OK?
		
		bbi1 = IntVec2((int)(bb1.x/py), (int)(bb1.y/pz));
		bbi2 = IntVec2((int)(bb2.x/py), (int)(bb2.y/pz));
		int sub_sizex = bbi2.x-bbi1.x+1;
		int sub_sizey = bbi2.y-bbi1.y+1;
	//	if(sub_sizey < 0) sub_sizey-= 2 * sub_sizey;

		CIS_Array_Image2D_short *binImg2D=NULL, *blobImg2D=NULL;
		binImg2D = new CIS_Array_Image2D_short(sub_sizex, sub_sizey);
		blobImg2D = new CIS_Array_Image2D_short(sub_sizex, sub_sizey);
		intDynArray blobRank2D;
		intVec2DynArray blobCentroid2D;
		short *binA2D, *blobA2D;
		int count=0;
		binA2D = binImg2D->GetArray();
		blobA2D = blobImg2D->GetArray();

		for(y=0,k=0; y<sub_sizey; y++)
		{
			for(x=0, k2=(y+bbi1.y)*sizey+bbi1.x; x<sub_sizex; x++, k++, k2++)
			{
				if(sagMaskA[k2]==maskStruct.mask_vertebralDisk) 
				{
					binA2D[k]=1; count++;
				}
				else binA2D[k]=0;
			}
		}
		if(count==0) continue;

		CIS_IPA_Erode(binImg2D, 1, true);
		CIS_IPA_Dilate(binImg2D, 1, true);
		NIH_Algo_Blob_Labelling(binImg2D, blobImg2D, blobRank2D, blobCentroid2D, false);
		int largestBlob=-1, largestSize=0;
		for(k=0; k<blobRank2D.GetSize(); k++)
		{
			if(blobRank2D[k]>largestSize)
			{
				largestSize = blobRank2D[k];
				largestBlob=k+1;
			}
		}

		for(y=0,k=0; y<sub_sizey; y++)
		{
			for(x=0, k2=(y+bbi1.y)*sizey+bbi1.x; x<sub_sizex; x++, k++, k2++)
			{
				if(binA2D[k]==1 && blobA2D[k]!=largestBlob) 
				{
					sagMaskA[k2] = maskStruct.mask_spongyBone;
					binA2D[k]=0;
				}
			}
		}

		// histogram analysis to use adaptive threshold
/*		int meanI, stdI, meanI2;
		meanI = stdI = meanI2=0;
		count = 0;
		for(y=0,k=0; y<sub_sizey; y++)
		{
			for(x=0, k2=(y+bbi1.y)*sizey+bbi1.x; x<sub_sizex; x++, k++, k2++)
			{
				if(sagMaskA[k2]==maskStruct.mask_vertebralDisk) 
				{
					count++;
					meanI += sagA[k2];
					meanI2 += sagA[k2]*sagA[k2];
				}
			}
		}
		if(count>0)
		{
			meanI /= count;
			stdI = (int)sqrt((double)(meanI2/count-meanI*meanI));
		}

		for(y=0,k=0; y<sub_sizey; y++)
		{
			for(x=0, k2=(y+bbi1.y)*sizey+bbi1.x; x<sub_sizex; x++, k++, k2++)
			{
				if(sagMaskA[k2]==maskStruct.mask_vertebralDisk && sagA[k2]>meanI+stdI) 
				{
					sagMaskA[k2] = maskStruct.mask_spongyBone;
					binA2D[k]=0;
				}
			}
		}
*/
		// searching for the real boundary of vertebral disk
		// search along the normalX

		Vec2 curPos, normalY, cPos, nPos;
		IntVec2 lastPosi, curPosi, cPosi, nPosi;
		short minInt, maxGrad, gray, ngray, grad;
		int nk, countD=0;

		double s, s1, plate1=0, plate2=0, valley=0;
		curPos = end1;
		normalY = normalX.perp();

		for(s=0; s<_sizeY; s+=py/2)
		{
			curPos = end1+normalX*s;
			curPosi = IntVec2((int)(curPos.x/py+0.5), (int)(curPos.y/pz+0.5));
			if(lastPosi==curPosi) continue;
			lastPosi = curPosi;
			count =0;
			minInt = 10000;
			// search the minimum intensity along normalY as the valley
			// one direction
			for(s1=0; s1<diskThick/2; s1+=py)
			{
				cPos = curPos+normalY*s1;
				cPosi = IntVec2((int)(cPos.x/py+0.5), (int)(cPos.y/pz+0.5));
				if(cPosi.y<1 || cPosi.y>sizez-2) continue;
				k = cPosi.x+cPosi.y*sizey;
				if(sagMaskA[k]!=maskStruct.mask_vertebralDisk) continue;
				count++;
				gray = (sagA[k]+sagA[k-sizey]+sagA[k+sizey]+sagA[k-1]+sagA[k+1])/5;
				if(gray<minInt)
				{
					minInt=gray;
					valley = s1;
				}
			}
			// opposite direction
			for(s1=0; s1<diskThick/2; s1+=py)
			{
				cPos = curPos-normalY*s1;
				cPosi = IntVec2((int)(cPos.x/py+0.5), (int)(cPos.y/pz+0.5));
				if(cPosi.y<1 || cPosi.y>sizez-2) continue;
				k = cPosi.x+cPosi.y*sizey;
				if(sagMaskA[k]!=maskStruct.mask_vertebralDisk) continue;
				count++;
				gray = (sagA[k]+sagA[k-sizey]+sagA[k+sizey]+sagA[k-1]+sagA[k+1])/5;
				if(gray<minInt)
				{
					minInt=gray;
					valley = -s1;
				}
			}	// for s1

			if(count<1) continue;
			countD++;

			// start from valley, search for plate1 and plate2 location based on gradient
			maxGrad=0;
			for(s1=py*2; s1<diskThick*2; s1+=py)
			{
				cPos = curPos+normalY*(valley+s1);
				cPosi = IntVec2((int)(cPos.x/py+0.5), (int)(cPos.y/pz+0.5));
				if(cPosi.y<bbi1.y || cPosi.y>bbi2.y) break; 
				if(cPosi.y<2) continue;
				k = cPosi.x+cPosi.y*sizey;
				if(sagMaskA[k]!=maskStruct.mask_vertebralDisk) continue;
				nPos = cPos+normalY;
				nPosi = IntVec2((int)(nPos.x/py+0.5), (int)(nPos.y/pz+0.5));
				nk = nPosi.x+nPosi.y*sizey;
				gray = (sagA[k]+sagA[k-sizey]+sagA[k+sizey]+sagA[k-1]+sagA[k+1])/5;
				ngray = (sagA[nk]+sagA[nk-sizey]+sagA[nk+sizey]+sagA[nk-1]+sagA[nk+1])/5;
				grad = ngray-gray;
				if(grad>maxGrad)
				{
					maxGrad=grad;
					plate1 = valley+s1;
				}
			}	// for s1
			if(countD<diskThick && plate1>diskThick/2) plate1=diskThick/2; // limit the range at the begining
			// erase everything after plate1
			if(maxGrad>0)
			{
				for(s1=plate1+py; s1<diskThick*2; s1+=py/2)
				{
					cPos = curPos+normalY*s1;
					cPosi = IntVec2((int)(cPos.x/py+0.5), (int)(cPos.y/pz+0.5));
					if(cPosi.y<bbi1.y || cPosi.y>bbi2.y) break; 
					if(cPosi.y<1) continue;
					k = cPosi.x+cPosi.y*sizey;
					if(sagMaskA[k]==maskStruct.mask_vertebralDisk) sagMaskA[k]=maskStruct.mask_spongyBone;
				}
			}

			// work on the opposite direction for erase
			maxGrad=0;
			for(s1=py*2; s1<diskThick*2; s1+=py)
			{
				cPos = curPos+normalY*(valley-s1);
				cPosi = IntVec2((int)(cPos.x/py+0.5), (int)(cPos.y/pz+0.5));
				if(cPosi.y<2) continue;
				k = cPosi.x+cPosi.y*sizey;
				if(cPosi.y<bbi1.y || cPosi.y>bbi2.y) break; 
				if(sagMaskA[k]!=maskStruct.mask_vertebralDisk) continue;
				nPos = cPos-normalY;
				nPosi = IntVec2((int)(nPos.x/py+0.5), (int)(nPos.y/pz+0.5));
				nk = nPosi.x+nPosi.y*sizey;
				gray = (sagA[k]+sagA[k-sizey]+sagA[k+sizey]+sagA[k-1]+sagA[k+1])/5;
				ngray = (sagA[nk]+sagA[nk-sizey]+sagA[nk+sizey]+sagA[nk-1]+sagA[nk+1])/5;
				grad = ngray-gray;
				if(grad>maxGrad)
				{
					maxGrad=grad;
					plate2 = s1-valley;
				}
			}	// for s1
			if(countD<diskThick && plate2>diskThick/2) plate2=diskThick/2;	// limit the range at the begining
			// erase everything after plate2
			if(maxGrad>0)
			{
				for(s1=plate2+py; s1<diskThick*2; s1+=py/2)
				{
					cPos = curPos-normalY*s1;
					cPosi = IntVec2((int)(cPos.x/py+0.5), (int)(cPos.y/pz+0.5));
					if(cPosi.y<bbi1.y || cPosi.y>bbi2.y) break; 
					if(cPosi.y<1) continue;
					k = cPosi.x+cPosi.y*sizey;
					if(sagMaskA[k]==maskStruct.mask_vertebralDisk) sagMaskA[k]=maskStruct.mask_spongyBone;
				}
			}
		}	// for s;

		// fill gaps
		count=0;
		for(y=0,k=0; y<sub_sizey; y++)
		{
			for(x=0, k2=(y+bbi1.y)*sizey+bbi1.x; x<sub_sizex; x++, k++, k2++)
			{
				if(sagMaskA[k2]==maskStruct.mask_vertebralDisk) 
				{
					binA2D[k]=1; count++;
				}
				else binA2D[k]=0;
			}
		}
		if(count==0) continue;

		CIS_IPA_Dilate(binImg2D, 1, true);
		CIS_IPA_Erode(binImg2D, 1, true);
		NIH_Algo_Blob_Labelling(binImg2D, blobImg2D, blobRank2D, blobCentroid2D, false);
		largestBlob=-1; largestSize=0;
		for(k=0; k<blobRank2D.GetSize(); k++)
		{
			if(blobRank2D[k]>largestSize)
			{
				largestSize = blobRank2D[k];
				largestBlob=k+1;
			}
		}

		for(y=0,k=0; y<sub_sizey; y++)
		{
			for(x=0, k2=(y+bbi1.y)*sizey+bbi1.x; x<sub_sizex; x++, k++, k2++)
			{
				if(binA2D[k]==1) 
				{
					if(blobA2D[k]!=largestBlob)
					{
						sagMaskA[k2] = maskStruct.mask_spongyBone;
						binA2D[k]=0;
					}
					else sagMaskA[k2] = maskStruct.mask_vertebralDisk;
				}
			}
		}

		// remove disjointed part

		// fit smooth boundaries

		delete binImg2D;
		delete blobImg2D;
	}	// for p

	return CIS_OK;
}


// compute graphical models for visualization purpose
//
int NIH_SpineSegmentation_ProjectColumnModel(CIS_Array_Image3D_short *img3D,
											 int projDir, int window_size, 
										   vec3DynArray &spinalCord3D,
											vec3DynArray &spineCenter, vec3DynArray &spineNormalX, vec3DynArray &spineNormalY, vec3DynArray &spineNormalBack,
											doubleDynArray &spineWidthLeft, doubleDynArray &spineWidthRight, doubleDynArray &spineWidthUp, doubleDynArray &spineWidthDown, 
											intDynArray &pedicleValley, intDynArray &backValley,
											CIS_2D_Model_Curve *projectedCord,
											CIS_2D_Model_Curve *leftColumn, CIS_2D_Model_Curve *rightColumn,
											CIS_2D_Model_Links *pedicleModel,
											bool debugMode)
{
	if(spinalCord3D.GetSize()<10) return CIS_ERROR;

	float vsizex = img3D->Get_Volume_SizeX();
	float vsizey = img3D->Get_Volume_SizeY();
	float vsizez = img3D->Get_Volume_SizeZ();
	float px = img3D->Get_Pixel_SizeX();
	float py = img3D->Get_Pixel_SizeY();
	float pz = img3D->Get_Pixel_SizeZ();
	int sizez = img3D->Num_Levels();

	int x, y, z, k;
	float fx, fy;
	vec2DynArray projA, leftA, rightA, normalA;
	float isizey, isizez;
	if(vsizey>vsizez)
	{
		isizey = window_size;
		isizez = (int)((vsizez/vsizey)*window_size);
	}
	else 
	{
		isizey = (int)((vsizey/vsizez)*window_size);
		isizez = window_size;
	}

	float projSize_x = isizey;
	float projSize_y = isizez;
	float projLoc_x = (window_size-isizey)/2;
	float projLoc_y = (window_size-isizez)/2;

	if(projDir==0)	// saggital
	{
		fx = projSize_x/vsizey;
		fy = projSize_y/vsizez;
	}
	else if(projDir==1)	// cornoral
	{
		fx = projSize_x/vsizex;
		fy = projSize_y/vsizez;
	}

	intDynArray filledZ;
	filledZ.SetSize(sizez);
	for(z=0; z<sizez; z++) 
	{
		filledZ[z]=-1;
	}
		
	for(k=0; k<spinalCord3D.GetSize(); k++)
	{
		z = img3D->GetSliceLevel(spinalCord3D[k].z,1);

		if(filledZ[z]==1 && filledZ[z-1]==-1) z=z-1;
		while(z<sizez-1 && filledZ[z]!=-1) z++;
		filledZ[z] = k;
	}

	projA.SetSize(spinalCord3D.GetSize());
	leftA.SetSize(spinalCord3D.GetSize());
	rightA.SetSize(spinalCord3D.GetSize());

	for(k=0; k<spinalCord3D.GetSize(); k++)
	{
		z = img3D->GetSliceLevel(spinalCord3D[k].z,1);

		if(projDir==0)
		{
			projA[k].x = spinalCord3D[k].y*fx+projLoc_x;
			projA[k].y = z*pz*fy+projLoc_y;
		}
		else if(projDir==1)
		{
			projA[k].x = spinalCord3D[k].x*fx+projLoc_x;
			projA[k].y = z*pz*fy+projLoc_y;
		}
	}

	projectedCord->UpdateControlList(projA);

	// define left and right column
	z = 0;
	for(k=0; k<projA.GetSize(); k++)
	{
		z = img3D->GetSliceLevel(spinalCord3D[k].z,1);

		if(projDir==0)
		{
			leftA[k].x = projA[k].x-spineWidthUp[z]*fx;
			rightA[k].x = projA[k].x+spineWidthDown[z]*fx;
		}
		else if(projDir==1)
		{
			leftA[k].x = projA[k].x-spineWidthLeft[z]*fx;
			rightA[k].x = projA[k].x+spineWidthRight[z]*fx;
		}

		leftA[k].y = projA[k].y;
		rightA[k].y = projA[k].y;
	}
		
	leftColumn->UpdateControlList(leftA);
	rightColumn->UpdateControlList(rightA);


	// set up the pedicle model
	Vec2 normal;
	int count;
	if(projDir==0)
	{
		pedicleModel->vertex_list.SetSize(0);
		pedicleModel->line_list.SetSize(0);
		pedicleModel->plot_radius = 0;
		count=0;
		// pedicleValley was computed from Coronal (sag) view
		for(k=0; k<pedicleValley.GetSize(); k++)
		{
			z = pedicleValley[k];
			if(filledZ[z]==-1) continue;

			normal.x = spineNormalX[z].y;
			normal.y = spineNormalX[z].z;
			pedicleModel->vertex_list.Add(projA[filledZ[z]]);
			pedicleModel->vertex_list.Add(projA[filledZ[z]]-normal*spineWidthUp[z]*fx);
			pedicleModel->line_list.Add(IntVec2(count, count+1));
			count+=2;
		}
		// add the back of vertebra
		for(k=0; k<backValley.GetSize(); k++)
		{
			z = backValley[k];
			if(filledZ[z]==-1) continue;

			normal.x = spineNormalBack[z].y;
			normal.y = spineNormalBack[z].z;
			pedicleModel->vertex_list.Add(projA[filledZ[z]]);
			pedicleModel->vertex_list.Add(projA[filledZ[z]]+normal*spineWidthDown[z]*fx);
			pedicleModel->line_list.Add(IntVec2(count, count+1));
			count+=2;
		}
	}
	else if(projDir==1)
	{
		// set up the pedicle model
		pedicleModel->vertex_list.SetSize(0);
		pedicleModel->line_list.SetSize(0);
		pedicleModel->plot_radius = 0;
		count=0;
		for(k=0; k<pedicleValley.GetSize(); k++)
		{
			z = pedicleValley[k];
			if(filledZ[z]==-1) continue;

			normal.x = spineNormalY[z].x;
			normal.y = spineNormalY[z].z;

			pedicleModel->vertex_list.Add(projA[filledZ[z]]);

			pedicleModel->vertex_list.Add(projA[filledZ[z]]-normal*spineWidthLeft[z]*fx);
			pedicleModel->vertex_list.Add(projA[filledZ[z]]+normal*spineWidthRight[z]*fx);
			pedicleModel->line_list.Add(IntVec2(count, count+1));
			pedicleModel->line_list.Add(IntVec2(count, count+2));
			count+=3;
		}
	}
	

	return CIS_OK;
}

// automatically trace the ROI for BMD calculate
// for each vertebra, need one ellipse ROI in the center of the vertebra body
// and one polygon around the spinal process on the same slice
//
int NIH_SpineSegmentation_BMD_ROI(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   VertebraStruct2D *vertebraTemplate, 
									   doubleDynArray &spineWidthLeft, doubleDynArray &spineWidthRight, doubleDynArray &spineWidthUp, doubleDynArray &spineWidthDown, 
									   intDynArray &pedicleValley, intDynArray &backValley, intDynArray &roiSlices,
									   bool debugMode,
									   const char* export_fn)
{
	int sizex, sizey, sizez, sizexy, sizex2;
	short *maskA, *imgA;
	int x, y, z, k, k2;
	float px, py, pz;

	sizex = maskImg3D->Num_Cols();
	sizey = maskImg3D->Num_Rows();
	sizez = maskImg3D->Num_Levels();
	px = img3D->Get_Pixel_SizeX();
	py = img3D->Get_Pixel_SizeY();
	pz = img3D->Get_Pixel_SizeZ();
	sizexy = sizex*sizey;
	sizex2 =sizex*2;

	maskA = maskImg3D->GetArray();
	imgA = img3D->GetArray();

	int vertNum = pedicleValley.GetSize();
	// variables to store Nicoli's histogram
	int binNumber, upperEnd_tr, lowerEnd_tr, upperEnd_fat, lowerEnd_fat;

	float **trab_ct, **fat_ct;
	int **trab_feq, **fat_feq;
	int p;

	binNumber = 100;
	upperEnd_tr = 500;
	lowerEnd_tr = 0;
	upperEnd_fat = 300;
	lowerEnd_fat = -200;

	trab_ct = new float*[vertNum];
	for(p=0; p<vertNum; p++) trab_ct[p] = new float[binNumber];
	fat_ct = new float*[vertNum];
	for(p=0; p<vertNum; p++) fat_ct[p] = new float[binNumber];
	trab_feq = new int*[vertNum];
	for(p=0; p<vertNum; p++) trab_feq[p] = new int[binNumber];
	fat_feq = new int*[vertNum];
	for(p=0; p<vertNum; p++) fat_feq[p] = new int[binNumber];

	roiSlices.SetSize(0);
	for(p=-1; p<pedicleValley.GetSize()-1; p++)
	{
		Vec2 diskCenter, diskRadius;
		int a, count;

		// get the middle slice
		if(p==-1) z = pedicleValley[0]/2;
		else z = (pedicleValley[p]+pedicleValley[p+1])/2;

		if(segInfo.diskBound1[z].x==-1) continue;						// bounding box not right
		if(segInfo.diskBound1[z].x>segInfo.diskBound2[z].x-2) continue;	// bounding box not right
		if(segInfo.diskBound1[z].y>segInfo.diskBound2[z].y-2) continue;	// bounding box not right
		if(segInfo.sprocessBound1[z].y>segInfo.sprocessBound2[z].y-2) continue;	// bounding box not right

		k = z*sizexy;
		roiSlices.Add(z);

		// init center
		diskCenter = (segInfo.diskBound1[z]+segInfo.diskBound2[z])/2;

		// adjust the y coordinate of diskCenter
		int uy, ly;
		// check the upper border
		uy=0;
		count=0;
		for(x=(int)diskCenter.x-2; x<=(int)diskCenter.x+2; x++)
		{
			for(y=segInfo.spineBound1[z].y; y<diskCenter.y; y++)
			{
				k2 = k+y*sizex+x;
				if((imgA[k2]>segPara.boneThresh || maskA[k2]==maskStruct.mask_corticalBone || maskA[k2]==maskStruct.mask_spongyBone)
					&& (imgA[k2+sizex]>segPara.boneThresh || maskA[k2+sizex]==maskStruct.mask_corticalBone || maskA[k2+sizex]==maskStruct.mask_spongyBone)
					&& (imgA[k2+sizex2]>segPara.boneThresh || maskA[k2+sizex2]==maskStruct.mask_corticalBone || maskA[k2+sizex2]==maskStruct.mask_spongyBone))
				{
					count++;
					uy += y;
					break;
				}
			}
		}
		if(count>0) uy/=count;
		else uy = segInfo.diskBound1[z].y;

		// check the lower border
		ly=0;
		count=0;
		for(x=(int)diskCenter.x-2; x<=(int)diskCenter.x+2; x++)
		{
			for(y=(diskCenter.y+segInfo.diskBound2[z].y)/2; y<segInfo.diskBound2[z].y+3; y++)
			{
				k2 = k+y*sizex+x;
				if(maskA[k2]==maskStruct.mask_spinalCord)
				{
					count++;
					ly += y;
					break;
				}
			}
		}
		if(count>0) ly/=count;
		else ly = segInfo.diskBound2[z].y;
		diskCenter.y = (uy+ly)/2;
		diskRadius.y = (ly-uy)/2;


		// adjust the x coordinate of diskCenter
		int rx, lx;
		// check the right border
		rx=0;
		count=0;
		for(y=(int)diskCenter.y-2; y<=(int)diskCenter.y+2; y++)
		{
			for(x=segInfo.spineBound2[z].x; x>diskCenter.x; x--)
			{
				k2 = k+y*sizex+x;
				if((imgA[k2]>segPara.boneThresh || maskA[k2]==maskStruct.mask_corticalBone || maskA[k2]==maskStruct.mask_spongyBone) 
					&& (imgA[k2-1]>segPara.boneThresh || maskA[k2-1]==maskStruct.mask_corticalBone || maskA[k2-1]==maskStruct.mask_spongyBone) 
					&& (imgA[k2-2]>segPara.boneThresh || maskA[k2-2]==maskStruct.mask_corticalBone || maskA[k2-2]==maskStruct.mask_spongyBone) )
				{
					count++;
					rx += x;
					break;
				}
			}
		}
		if(count>0) rx/=count;
		else rx = segInfo.diskBound2[z].x;

		// check the left border
		lx=0;
		count=0;
		for(y=(int)diskCenter.y-2; y<=(int)diskCenter.y+2; y++)
		{
			for(x=segInfo.spineBound1[z].x; x<diskCenter.x; x++)
			{
				k2 = k+y*sizex+x;
				if((imgA[k2]>segPara.boneThresh || maskA[k2]==maskStruct.mask_corticalBone || maskA[k2]==maskStruct.mask_spongyBone)
					&& (imgA[k2+1]>segPara.boneThresh || maskA[k2+1]==maskStruct.mask_corticalBone || maskA[k2+1]==maskStruct.mask_spongyBone)
					&& (imgA[k2+2]>segPara.boneThresh || maskA[k2+2]==maskStruct.mask_corticalBone || maskA[k2+2]==maskStruct.mask_spongyBone))
				{
					count++;
					lx += x;
					break;
				}
			}
		}
		if(count>0) lx/=count;
		else lx = segInfo.diskBound1[z].x;

		diskCenter.x = (rx+lx)/2;
		diskRadius.x = (rx-lx)/2;


		// get the lower muscle and fat ROI
		Vec2 mroiCenter, mroiLeft, mroiRight;
		vec2DynArray lowerBorder;
		mroiCenter = (segInfo.sprocessBound1[z]+segInfo.sprocessBound2[z]*2)/3;
///		mroiRight = mroiCenter+Vec2(diskRadius.x*2.5,0);
///		mroiLeft = mroiCenter+Vec2(-diskRadius.x*2.5,0);
		mroiRight = mroiCenter+Vec2(diskRadius.x*2,0);	// new smaller rois
		mroiLeft = mroiCenter+Vec2(-diskRadius.x*2,0);

		for(x=(int)mroiLeft.x; x<=(int)mroiRight.x; x++)
		{
			for(y=segInfo.sprocessBound2[z].y; y<sizey; y++)
			{
				if(imgA[k+y*sizex+x]<segPara.airThresh)
				{
					lowerBorder.Add(Vec2(x,y-10));
					break;
				}
			}
		}

		// fit smooth spline curves
		int bernstein_power = 5, piece_size=10;
		doubleDynArray xt, yt, yt2;

		xt.SetSize(lowerBorder.GetSize());
		yt.SetSize(lowerBorder.GetSize());
		yt2.SetSize(lowerBorder.GetSize());

		for(k2=0; k2<lowerBorder.GetSize(); k2++)
		{
			xt[k2] = lowerBorder[k2].x;
			yt[k2] = lowerBorder[k2].y;
		}
		PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
		for(k2=0; k2<lowerBorder.GetSize(); k2++)
		{
			lowerBorder[k2].y = yt2[k2];
		}

		// erase all mask 
		k=z*sizexy;
		for(k2=0; k2<sizexy; k2++)
			maskA[k+k2] = maskStruct.mask_body;

		// diskROI
		// only use 2/3 of the disk as the roi
///		diskRadius *= 0.66;
		diskRadius *= 0.5;		// a new radius
		diskRadius *= 2;	// need diameter in the model

		// compute histogram assigned the mask
		float minG, maxG, stepG;
		int ind;
		short hu;
		vec2DynArray scan;
		CIS_2D_ROI_Ellipse *ell;
		ell = new CIS_2D_ROI_Ellipse(diskCenter, diskRadius);
		ell->GetEntireScan(scan);
		delete ell;

		minG = 100000;
		maxG = 0;
		for(k2=0; k2<scan.GetSize(); k2++)
		{
			ind = k+(int)scan[k2].x+(int)scan[k2].y*sizex;
			if(imgA[ind]<minG) minG=imgA[ind];
			if(imgA[ind]>maxG) maxG=imgA[ind];
			maskA[ind] = maskStruct.mask_boneMetasis;
		}
///		stepG = (maxG-minG)/255;
///		trab_ct[p+1][0] = minG-1024;	// convert to HU
///		trab_feq[p+1][0]=0;
		stepG = (upperEnd_tr-lowerEnd_tr)/binNumber;
		trab_ct[p+1][0] = lowerEnd_tr;	// convert to HU
		trab_feq[p+1][0]=0;
		for(ind=1; ind<binNumber; ind++) 
		{
			trab_ct[p+1][ind] = trab_ct[p+1][ind-1]+stepG;
			trab_feq[p+1][ind]=0;
		}
		for(k2=0; k2<scan.GetSize(); k2++)
		{
			ind = k+(int)scan[k2].x+(int)scan[k2].y*sizex;
			hu = imgA[ind]-1024;
			if(hu>lowerEnd_tr && hu<upperEnd_tr)
			{
				trab_feq[p+1][(int)((hu-lowerEnd_tr)/stepG)]++;
			}

			// one slice up
			ind = k+sizexy+(int)scan[k2].x+(int)scan[k2].y*sizex;
			hu = imgA[ind]-1024;
			if(hu>lowerEnd_tr && hu<upperEnd_tr) 
			{
				trab_feq[p+1][(int)((hu-lowerEnd_tr)/stepG)]++;
			}

			// one slice down
			ind = k-sizexy+(int)scan[k2].x+(int)scan[k2].y*sizex;
			hu = imgA[ind]-1024;
			if(hu>lowerEnd_tr && hu<upperEnd_tr) 
			{
				trab_feq[p+1][(int)((hu-lowerEnd_tr)/stepG)]++;
			}

			// two slice up
			ind = k+sizexy*2+(int)scan[k2].x+(int)scan[k2].y*sizex;
			hu = imgA[ind]-1024;
			if(hu>lowerEnd_tr && hu<upperEnd_tr) 
			{
				trab_feq[p+1][(int)((hu-lowerEnd_tr)/stepG)]++;
			}

			// two slice down
			ind = k-sizexy*2+(int)scan[k2].x+(int)scan[k2].y*sizex;
			hu = imgA[ind]-1024;
			if(hu>lowerEnd_tr && hu<upperEnd_tr) 
			{
				trab_feq[p+1][(int)((hu-lowerEnd_tr)/stepG)]++;
			}

		}

		// muscle ROI
		lowerBorder.Add(mroiRight);
		lowerBorder.Add(mroiLeft);

		CIS_2D_ROI_Polygon *polyROI;
		polyROI = new CIS_2D_ROI_Polygon(lowerBorder);
		polyROI->GetEntireScan(scan);
		delete polyROI;

		minG = 100000;
		maxG = 0;
		k=z*sizexy;
		for(k2=0; k2<scan.GetSize(); k2++)
		{
			ind = k+(int)scan[k2].x+(int)scan[k2].y*sizex;
			if(imgA[ind]<minG) minG=imgA[ind];
			if(imgA[ind]>maxG) maxG=imgA[ind];
			maskA[ind] = maskStruct.mask_boneMetasis;
		}
///		stepG = (maxG-minG)/255;
///		fat_ct[p+1][0] = minG-1024;		// convert to HU
///		fat_feq[p+1][0]=0;
		stepG = (upperEnd_fat-lowerEnd_fat)/binNumber;
		fat_ct[p+1][0] = lowerEnd_fat;	// convert to HU
		fat_feq[p+1][0]=0;

		for(ind=1; ind<binNumber; ind++) 
		{
			fat_ct[p+1][ind] = fat_ct[p+1][ind-1]+stepG;
			fat_feq[p+1][ind]=0;
		}
		for(k2=0; k2<scan.GetSize(); k2++)
		{
			ind = k+(int)scan[k2].x+(int)scan[k2].y*sizex;
			hu = imgA[ind]-1024;
			if(hu>lowerEnd_fat && hu<upperEnd_fat)
			{
				fat_feq[p+1][(int)((hu-lowerEnd_fat)/stepG)]++;
			}

			// one slice up
			ind = k+sizexy+(int)scan[k2].x+(int)scan[k2].y*sizex;
			hu = imgA[ind]-1024;
			if(hu>lowerEnd_fat && hu<upperEnd_fat)
			{
				fat_feq[p+1][(int)((hu-lowerEnd_fat)/stepG)]++;
			}

			// one slice down
			ind = k-sizexy+(int)scan[k2].x+(int)scan[k2].y*sizex;
			hu = imgA[ind]-1024;
			if(hu>lowerEnd_fat && hu<upperEnd_fat)
			{
				fat_feq[p+1][(int)((hu-lowerEnd_fat)/stepG)]++;
			}
			// two slice up
			ind = k+sizexy+(int)scan[k2].x+(int)scan[k2].y*sizex;
			hu = imgA[ind]-1024;
			if(hu>lowerEnd_fat && hu<upperEnd_fat)
			{
				fat_feq[p+1][(int)((hu-lowerEnd_fat)/stepG)]++;
			}

			// two slice down
			ind = k-sizexy+(int)scan[k2].x+(int)scan[k2].y*sizex;
			hu = imgA[ind]-1024;
			if(hu>lowerEnd_fat && hu<upperEnd_fat)
			{
				fat_feq[p+1][(int)((hu-lowerEnd_fat)/stepG)]++;
			}
		}
	}	// for p

	// output Nicoli's histogram into a file
	FILE *e_fp;
	if((e_fp = fopen(export_fn, "w"))!=NULL)
	{
		// print the header
		fprintf(e_fp,"t12_trabctnum,t12_trabfreq,t12_fmctnum,t12_fmfreq,l1_trabctnum,l1_trabfreq,l1_fmctnum,l1_fmfreq,");
		fprintf(e_fp,"l2_trabctnum,l2_trabfreq,l2_fmctnum,l2_fmfreq,l3_trabctnum,l3_trabfreq,l3_fmctnum,l3_fmfreq,");
		fprintf(e_fp,"l4_trabctnum,l4_trabfreq,l4_fmctnum,l4_fmfreq,l5_trabctnum,l5_trabfreq,l5_fmctnum,l5_fmfreq,\n");

		int st_p;
		if(segInfo.t12_vertebra>=0) st_p=segInfo.t12_vertebra;
		else st_p=0;

		for(int i=0; i<binNumber; i++)
		{
			for(p=st_p; p<vertNum && p<st_p+6; p++)
			{
				fprintf(e_fp, "%.3f,%d,%.3f,%d,",trab_ct[p][i], trab_feq[p][i], fat_ct[p][i], fat_feq[p][i]);
			}
			fprintf(e_fp, "\n");
		}
		fclose(e_fp);
	}


	// release the memory
	for(p=0; p<vertNum; p++) 
	{
		delete trab_ct[p];
		delete fat_ct[p];
		delete trab_feq[p];
		delete fat_feq[p];
	}
	delete trab_ct;
	delete trab_feq;
	delete fat_ct;
	delete fat_feq;

	return CIS_OK;
}

// automatically trace the ROI for BMD calculate, 3D technique
// for each vertebra, Segment the entire vertebra
//
// and one polygon around the spinal process on the same slice
//
int NIH_SpineSegmentation_BMD_ROI_3D(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   VertebraStruct2D *vertebraTemplate, 
									   vec3DynArray &spineCenter, vec3DynArray &spineNormalX, vec3DynArray &spineNormalY, vec3DynArray &spineNormalBack, 
									   doubleDynArray &spineWidthLeft, doubleDynArray &spineWidthRight, doubleDynArray &spineWidthUp, doubleDynArray &spineWidthDown, 
									   intDynArray &pedicleValley, intDynArray &backValley, intDynArray &roiSlices,
									   bool debugMode,
									   const char* export_fn)
{
	int sizex, sizey, sizez, sizexy, sizex2, sizexyz;
	short *maskA, *imgA;
	int x, y, z, k, k2, mz;
	float px, py, pz;

	sizex = maskImg3D->Num_Cols();
	sizey = maskImg3D->Num_Rows();
	sizez = maskImg3D->Num_Levels();
	px = img3D->Get_Pixel_SizeX();
	py = img3D->Get_Pixel_SizeY();
	pz = img3D->Get_Pixel_SizeZ();
	sizexy = sizex*sizey;
	sizexyz = sizexy*sizez;
	sizex2 =sizex*2;

	maskA = maskImg3D->GetArray();
	imgA = img3D->GetArray();

	CIS_Array_Image3D_uchar *binImg3D;
	unsigned char *binA3D;

	binImg3D = new CIS_Array_Image3D_uchar(sizex, sizey, sizez);
	binA3D = binImg3D->GetArray();

	float erode_mm = 5;
	int erode_xy, erode_z;
	erode_xy = (int)(erode_mm/px);
	erode_z = (int)(erode_mm/pz)-1;

	int vertNum = pedicleValley.GetSize();
	// variables to store Nicoli's histogram
	int binNumber, upperEnd_tr, lowerEnd_tr, upperEnd_fat, lowerEnd_fat;

	float **trab_ct, **fat_ct;
	int **trab_feq, **fat_feq;
	int p;

///	binNumber = 100;
	binNumber = 500;
	upperEnd_tr = 500;
	lowerEnd_tr = 0;
	upperEnd_fat = 300;
	lowerEnd_fat = -200;

	trab_ct = new float*[vertNum];
	for(p=0; p<vertNum; p++) trab_ct[p] = new float[binNumber];
	fat_ct = new float*[vertNum];
	for(p=0; p<vertNum; p++) fat_ct[p] = new float[binNumber];
	trab_feq = new int*[vertNum];
	for(p=0; p<vertNum; p++) trab_feq[p] = new int[binNumber];
	fat_feq = new int*[vertNum];
	for(p=0; p<vertNum; p++) fat_feq[p] = new int[binNumber];

	roiSlices.SetSize(0);
	for(p=-1; p<pedicleValley.GetSize()-1; p++)
	{
		Vec2 diskCenter, diskRadius;
		int a, count;
		int sz;

		// get the middle slice
		if(p==-1)
		{
			mz = pedicleValley[0]/2;
			sz = 1;
		}
		else 
		{
			sz = pedicleValley[p];
			mz = (pedicleValley[p]+pedicleValley[p+1])/2;
		}

		if(segInfo.diskBound1[mz].x==-1) continue;						// bounding box not right
		if(segInfo.diskBound1[mz].x>segInfo.diskBound2[mz].x-2) continue;	// bounding box not right
		if(segInfo.diskBound1[mz].y>segInfo.diskBound2[mz].y-2) continue;	// bounding box not right
		if(segInfo.sprocessBound1[mz].y>segInfo.sprocessBound2[mz].y-2) continue;	// bounding box not right

		roiSlices.Add(mz);

		// get the mean direction and width
		Vec3 startP, scanV, cV, ccV, dirX, dirY, mdirX, mdirY, normal;
		double stepSize, sx, sy, _sizeX, _sizeY;
		IntVec3 cVi, ccVi;
		count=0;
		mdirX = mdirY = Vec3(0,0,0);
		_sizeX = _sizeY = 0;
		for(z=sz; z<pedicleValley[p+1]; z++)
		{
			if(z<=0) continue;
			dirY = spineNormalX[z];
			dirX = spineNormalY[z];
			if(dirY.len()==0 || dirX.len()==0) continue;
			normal = dirX % dirY;
			dirX = dirY % normal;
			dirY = normal % dirX;
			if(dirY.y>0)
			{
				dirY = -dirY;
			}
			dirX = dirY % normal;

			count++;
			mdirX += dirX;
			mdirY += dirY;
			_sizeX += (spineWidthLeft[z]+spineWidthRight[z])*1.5;
			_sizeY += (spineWidthUp[z])*1.5;
		}

		if(count>0)
		{
			_sizeX /= (double)count;
			_sizeY /= (double)count;
			mdirX = mdirX.normalize();
			mdirY = mdirY.normalize();
			normal = mdirX % mdirY;
			mdirX = mdirY % normal;
		}

		float stepG;
		int ind, ck;
		short hu;

		// initialize the histogram
		stepG = (upperEnd_tr-lowerEnd_tr)/binNumber;
		trab_ct[p+1][0] = lowerEnd_tr;	// convert to HU
		trab_feq[p+1][0]=0;
		for(ind=1; ind<binNumber; ind++) 
		{
			trab_ct[p+1][ind] = trab_ct[p+1][ind-1]+stepG;
			trab_feq[p+1][ind]=0;
		}

		stepG = (upperEnd_fat-lowerEnd_fat)/binNumber;
		fat_ct[p+1][0] = lowerEnd_fat;	// convert to HU
		fat_feq[p+1][0]=0;
		for(ind=1; ind<binNumber; ind++) 
		{
			fat_ct[p+1][ind] = fat_ct[p+1][ind-1]+stepG;
			fat_feq[p+1][ind]=0;
		}

		// initialize the erosion mask
		for(ck=0; ck<sizexyz; ck++) binA3D[ck]=0;

		// go through all slices occupied by this vertebra
		//
		for(z=sz; z<pedicleValley[p+1]; z++)
		{
			stepG = (upperEnd_tr-lowerEnd_tr)/binNumber;
			if(segInfo.cordCenter[z].x<0) continue;
			startP = Vec3(segInfo.cordCenter[z].x*px, segInfo.cordCenter[z].y*py, z*pz);
			stepSize = px/2;
			scanV = startP;
			for(sy=0; sy<_sizeY; sy+=stepSize, scanV += (dirY*stepSize))
			{
				for(int dir=-1; dir<=1; dir+=2)
				{
					for(sx=0; sx<_sizeX/2; sx+=stepSize)
					{
						cV = scanV+dirX*(sx*dir);
						cVi = IntVec3(cV.x/px, cV.y/py, cV.z/pz);
						ck = cVi.x+cVi.y*sizex+cVi.z*sizexy;
						if(cVi.z<1) continue;

						if(maskA[ck]==maskStruct.mask_corticalBone || maskA[ck]==maskStruct.mask_spongyBone)
						{
/*							hu = imgA[ck]-1024;
							if(hu>lowerEnd_tr && hu<upperEnd_tr)
							{
								trab_feq[p+1][(int)((hu-lowerEnd_tr)/stepG)]++;
							}
							// include these pixels in 3D ROI of the vertebra
							maskA[ck] = maskStruct.mask_boneMetasis;
*/
							binA3D[ck]=1;
						}	// if maskA
						if(maskA[ck-1]==maskStruct.mask_corticalBone || maskA[ck-1]==maskStruct.mask_spongyBone)
						{
							binA3D[ck-1]=1;
						}
						if(maskA[ck+1]==maskStruct.mask_corticalBone || maskA[ck+1]==maskStruct.mask_spongyBone)
						{
							binA3D[ck+1]=1;
						}
						if(maskA[ck-sizex]==maskStruct.mask_corticalBone || maskA[ck-sizex]==maskStruct.mask_spongyBone)
						{
							binA3D[ck-sizex]=1;
						}
						if(maskA[ck+sizex]==maskStruct.mask_corticalBone || maskA[ck+sizex]==maskStruct.mask_spongyBone)
						{
							binA3D[ck+sizex]=1;
						}
						if(maskA[ck-sizexy]==maskStruct.mask_corticalBone || maskA[ck-sizexy]==maskStruct.mask_spongyBone)
						{
							binA3D[ck-sizexy]=1;
						}
						if(maskA[ck+sizexy]==maskStruct.mask_corticalBone || maskA[ck+sizexy]==maskStruct.mask_spongyBone)
						{
							binA3D[ck+sizexy]=1;
						}
					}	// for sx
				}	// for dir
			}	//for sy

			// Get muscle-fat region
			//
			// init center
			diskCenter = (segInfo.diskBound1[z]+segInfo.diskBound2[z])/2;

			k = z*sizexy;
			// adjust the y coordinate of diskCenter
			int uy, ly;
			// check the upper border
			uy=0;
			count=0;
			for(x=(int)diskCenter.x-2; x<=(int)diskCenter.x+2; x++)
			{
				for(y=segInfo.spineBound1[z].y; y<diskCenter.y; y++)
				{
					k2 = k+y*sizex+x;
					if((imgA[k2]>segPara.boneThresh || maskA[k2]==maskStruct.mask_corticalBone || maskA[k2]==maskStruct.mask_spongyBone)
						&& (imgA[k2+sizex]>segPara.boneThresh || maskA[k2+sizex]==maskStruct.mask_corticalBone || maskA[k2+sizex]==maskStruct.mask_spongyBone)
						&& (imgA[k2+sizex2]>segPara.boneThresh || maskA[k2+sizex2]==maskStruct.mask_corticalBone || maskA[k2+sizex2]==maskStruct.mask_spongyBone))
					{
						count++;
						uy += y;
						break;
					}
				}
			}
			if(count>0) uy/=count;
			else uy = segInfo.diskBound1[z].y;

			// check the lower border
			ly=0;
			count=0;
			for(x=(int)diskCenter.x-2; x<=(int)diskCenter.x+2; x++)
			{
				for(y=(diskCenter.y+segInfo.diskBound2[z].y)/2; y<segInfo.diskBound2[z].y+3; y++)
				{
					k2 = k+y*sizex+x;
					if(maskA[k2]==maskStruct.mask_spinalCord)
					{
						count++;
						ly += y;
						break;
					}
				}
			}
			if(count>0) ly/=count;
			else ly = segInfo.diskBound2[z].y;
			diskCenter.y = (uy+ly)/2;
			diskRadius.y = (ly-uy)/2;


			// adjust the x coordinate of diskCenter
			int rx, lx;
			// check the right border
			rx=0;
			count=0;
			for(y=(int)diskCenter.y-2; y<=(int)diskCenter.y+2; y++)
			{
				for(x=segInfo.spineBound2[z].x; x>diskCenter.x; x--)
				{
					k2 = k+y*sizex+x;
					if((imgA[k2]>segPara.boneThresh || maskA[k2]==maskStruct.mask_corticalBone || maskA[k2]==maskStruct.mask_spongyBone) 
						&& (imgA[k2-1]>segPara.boneThresh || maskA[k2-1]==maskStruct.mask_corticalBone || maskA[k2-1]==maskStruct.mask_spongyBone) 
						&& (imgA[k2-2]>segPara.boneThresh || maskA[k2-2]==maskStruct.mask_corticalBone || maskA[k2-2]==maskStruct.mask_spongyBone) )
					{
						count++;
						rx += x;
						break;
					}
				}
			}
			if(count>0) rx/=count;
			else rx = segInfo.diskBound2[z].x;

			// check the left border
			lx=0;
			count=0;
			for(y=(int)diskCenter.y-2; y<=(int)diskCenter.y+2; y++)
			{
				for(x=segInfo.spineBound1[z].x; x<diskCenter.x; x++)
				{
					k2 = k+y*sizex+x;
					if((imgA[k2]>segPara.boneThresh || maskA[k2]==maskStruct.mask_corticalBone || maskA[k2]==maskStruct.mask_spongyBone)
						&& (imgA[k2+1]>segPara.boneThresh || maskA[k2+1]==maskStruct.mask_corticalBone || maskA[k2+1]==maskStruct.mask_spongyBone)
						&& (imgA[k2+2]>segPara.boneThresh || maskA[k2+2]==maskStruct.mask_corticalBone || maskA[k2+2]==maskStruct.mask_spongyBone))
					{
						count++;
						lx += x;
						break;
					}
				}
			}
			if(count>0) lx/=count;
			else lx = segInfo.diskBound1[z].x;

			diskCenter.x = (rx+lx)/2;
			diskRadius.x = (rx-lx)/2;


			// get the lower muscle and fat ROI
			Vec2 mroiCenter, mroiLeft, mroiRight;
			vec2DynArray lowerBorder;
			mroiCenter = (segInfo.sprocessBound1[z]+segInfo.sprocessBound2[z]*2)/3;
			mroiRight = mroiCenter+Vec2(diskRadius.x*2,0);	// new smaller rois
			mroiLeft = mroiCenter+Vec2(-diskRadius.x*2,0);

			for(x=(int)mroiLeft.x; x<=(int)mroiRight.x; x++)
			{
				for(y=segInfo.sprocessBound2[z].y; y<sizey; y++)
				{
					if(imgA[k+y*sizex+x]<segPara.airThresh)
					{
						lowerBorder.Add(Vec2(x,y-10));
						break;
					}
				}
			}

			// fit smooth spline curves dor muscle border
			int bernstein_power = 5, piece_size=10;
			doubleDynArray xt, yt, yt2;

			xt.SetSize(lowerBorder.GetSize());
			yt.SetSize(lowerBorder.GetSize());
			yt2.SetSize(lowerBorder.GetSize());

			for(k2=0; k2<lowerBorder.GetSize(); k2++)
			{
				xt[k2] = lowerBorder[k2].x;
				yt[k2] = lowerBorder[k2].y;
			}
			PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
			for(k2=0; k2<lowerBorder.GetSize(); k2++)
			{
				lowerBorder[k2].y = yt2[k2];
			}

			// compute histogram assigned the mask
			vec2DynArray scan;

			// muscle ROI
			lowerBorder.Add(mroiRight);
			lowerBorder.Add(mroiLeft);

			CIS_2D_ROI_Polygon *polyROI;
			polyROI = new CIS_2D_ROI_Polygon(lowerBorder);
			polyROI->GetEntireScan(scan);
			delete polyROI;

			k=z*sizexy;
			stepG = (upperEnd_fat-lowerEnd_fat)/binNumber;
			for(k2=0; k2<scan.GetSize(); k2++)
			{
				ind = k+(int)scan[k2].x+(int)scan[k2].y*sizex;
				if(maskA[ind]==maskStruct.mask_corticalBone || maskA[ind]==maskStruct.mask_spongyBone) continue;
///				maskA[ind] = maskStruct.mask_boneMetasis;
				hu = imgA[ind]-1024;
				if(hu>lowerEnd_fat && hu<upperEnd_fat)
				{
					fat_feq[p+1][(int)((hu-lowerEnd_fat)/stepG)]++;
				}
			}
		}	// for z

		// get the eroded mask
		CIS_IPA_3D_Dilate(binImg3D, 2, true);
		CIS_IPA_3D_Erode(binImg3D, 2, true);

		int erode_stepz, erode_stepx, erode_stepy;
		erode_stepx = erode_xy;
		erode_stepy = erode_xy*sizex;
		erode_stepz = erode_z*sizexy;
		for(z=erode_z; z<sizez-erode_z-1; z++)
		{
			ck=z*sizexy;
			for(y=0; y<sizey; y++)
			{
				for(x=0;x<sizex; x++, ck++)
				{
					if(binA3D[ck]==1)
					{
						if(binA3D[ck-erode_stepz]==1 && binA3D[ck+erode_stepz]==1 &&
							binA3D[ck-erode_stepy]==1 && binA3D[ck+erode_stepy]==1 &&
							binA3D[ck-erode_stepx]==1 && binA3D[ck+erode_stepx]==1)
						{
							hu = imgA[ck]-1024;
							if(hu>lowerEnd_tr && hu<upperEnd_tr)
							{
								trab_feq[p+1][(int)((hu-lowerEnd_tr)/stepG)]++;
							}
							// include these pixels in 3D ROI of the vertebra
							maskA[ck] = maskStruct.mask_boneMetasis;						}
					}
				}
			}
		}	// for z

	}	// for p

	// output Nicoli's histogram into a file
	FILE *e_fp;
	if((e_fp = fopen(export_fn, "w"))!=NULL)
	{
		// print the header
		fprintf(e_fp,"t12_trabctnum,t12_trabfreq,t12_fmctnum,t12_fmfreq,l1_trabctnum,l1_trabfreq,l1_fmctnum,l1_fmfreq,");
		fprintf(e_fp,"l2_trabctnum,l2_trabfreq,l2_fmctnum,l2_fmfreq,l3_trabctnum,l3_trabfreq,l3_fmctnum,l3_fmfreq,");
		fprintf(e_fp,"l4_trabctnum,l4_trabfreq,l4_fmctnum,l4_fmfreq,l5_trabctnum,l5_trabfreq,l5_fmctnum,l5_fmfreq,\n");

		int st_p;
		if(segInfo.t12_vertebra>=0) st_p=segInfo.t12_vertebra;
		else st_p=0;

		for(int i=0; i<binNumber; i++)
		{
			for(p=st_p; p<vertNum && p<st_p+6; p++)
			{
				fprintf(e_fp, "%.3f,%d,%.3f,%d,",trab_ct[p][i], trab_feq[p][i], fat_ct[p][i], fat_feq[p][i]);
			}
			fprintf(e_fp, "\n");
		}

		// compute mean and standard deviation of the bone ROI using histogram
		for(p=st_p; p<vertNum && p<st_p+6; p++)
		{
			long sum=0;
			long sum2=0;
			long count=0;
			for(int i=0; i<binNumber; i++)
			{
				count += trab_feq[p][i];
				sum += trab_feq[p][i]*trab_ct[p][i];
				sum2 += trab_feq[p][i]*(trab_ct[p][i]*trab_ct[p][i]);
			}

			long dev=0;
			if(count>0) 
			{
				sum/=count;
				dev = (long)sqrt((double)(sum2/count-sum*sum));
			}

			// export
			fprintf(e_fp, "Vert T%d, %d, %d\n",p+1, sum, dev);
		}


		fclose(e_fp);
	}


	// release the memory
	for(p=0; p<vertNum; p++) 
	{
		delete trab_ct[p];
		delete fat_ct[p];
		delete trab_feq[p];
		delete fat_feq[p];
	}
	delete trab_ct;
	delete trab_feq;
	delete fat_ct;
	delete fat_feq;
	delete binImg3D;

	return CIS_OK;
}





int NIH_SpineSegmentation_BMD_ROI_L1L2(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   VertebraStruct2D *vertebraTemplate, 
									   doubleDynArray &spineWidthLeft, doubleDynArray &spineWidthRight, doubleDynArray &spineWidthUp, doubleDynArray &spineWidthDown, 
									   intDynArray &pedicleValley, intDynArray &backValley, intDynArray &roiSlices,
									   bool debugMode,
									   const char* export_fn,
									   int l1_slice, CIS_Array_Image3D_short *l1_img,
									   int l2_slice, CIS_Array_Image3D_short *l2_img)
{
	int sizex, sizey, sizez, sizexy, sizex2;
	short *maskA, *imgA, *imgAL;
	int x, y, z, k, k2;
	float px, py, pz;

	sizex = maskImg3D->Num_Cols();
	sizey = maskImg3D->Num_Rows();
	sizez = maskImg3D->Num_Levels();
	px = img3D->Get_Pixel_SizeX();
	py = img3D->Get_Pixel_SizeY();
	pz = img3D->Get_Pixel_SizeZ();
	sizexy = sizex*sizey;
	sizex2 =sizex*2;

	maskA = maskImg3D->GetArray();
	imgA = img3D->GetArray();

	int vertNum = 2;		// only l1 and l2
	// variables to store Nicoli's histogram
	int binNumber, upperEnd_tr, lowerEnd_tr, upperEnd_fat, lowerEnd_fat;

	float **trab_ct, **fat_ct;
	int **trab_feq, **fat_feq;
	int p;

	binNumber = 100;
	upperEnd_tr = 500;
	lowerEnd_tr = 0;
	upperEnd_fat = 300;
	lowerEnd_fat = -200;

	trab_ct = new float*[vertNum];
	for(p=0; p<vertNum; p++) trab_ct[p] = new float[binNumber];
	fat_ct = new float*[vertNum];
	for(p=0; p<vertNum; p++) fat_ct[p] = new float[binNumber];
	trab_feq = new int*[vertNum];
	for(p=0; p<vertNum; p++) trab_feq[p] = new int[binNumber];
	fat_feq = new int*[vertNum];
	for(p=0; p<vertNum; p++) fat_feq[p] = new int[binNumber];

	imgAL = l1_img->GetArray();
	for(k2=0; k2<sizexy; k2++) imgAL[k2]+=1024;
	for(y=0; y<sizey/2; y++)
		for(x=0; x<sizex; x++)
		{
			k2= y*sizex+x;
			k = (sizey-1-y)*sizex+x;
			p=imgAL[k2]; imgAL[k2]=imgAL[k]; imgAL[k]=p;
		}
	k=l1_slice*sizexy;
	for(k2=0; k2<sizexy; k2++) imgA[k+k2]=imgAL[k2];
	imgAL = l2_img->GetArray();
	for(k2=0; k2<sizexy; k2++) imgAL[k2]+=1024;
	for(y=0; y<sizey/2; y++)
		for(x=0; x<sizex; x++)
		{
			k2= y*sizex+x;
			k = (sizey-1-y)*sizex+x;
			p=imgAL[k2]; imgAL[k2]=imgAL[k]; imgAL[k]=p;
		}
	k=l2_slice*sizexy;
	for(k2=0; k2<sizexy; k2++) imgA[k+k2]=imgAL[k2];

	roiSlices.SetSize(0);
	for(p=-1; p<1; p++)
	{
		Vec2 diskCenter, diskRadius;
		int a, count;

		// get the middle slice
		if(p==-1) 
		{
			z = l1_slice;
			imgAL = l1_img->GetArray();
		}
		else 
		{
			z = l2_slice;
			imgAL = l2_img->GetArray();
		}

		if(segInfo.diskBound1[z].x==-1 || 
			segInfo.diskBound1[z].x>segInfo.diskBound2[z].x-2 || 
			segInfo.diskBound1[z].y>segInfo.diskBound2[z].y-2 ||
			segInfo.sprocessBound1[z].y>segInfo.sprocessBound2[z].y-2)
		{
			printf("Unable to generate ROI on slice %d\n", z);
			continue;	// bounding box not right
		}

		k = z*sizexy;
		roiSlices.Add(z);

		// init center
		diskCenter = (segInfo.diskBound1[z]+segInfo.diskBound2[z])/2;

		// adjust the y coordinate of diskCenter
		int uy, ly;
		// check the upper border
		uy=0;
		count=0;
		for(x=(int)diskCenter.x-2; x<=(int)diskCenter.x+2; x++)
		{
			for(y=segInfo.spineBound1[z].y; y<diskCenter.y; y++)
			{
				k2 = k+y*sizex+x;
				if((imgA[k2]>segPara.boneThresh || maskA[k2]==maskStruct.mask_corticalBone || maskA[k2]==maskStruct.mask_spongyBone)
					&& (imgA[k2+sizex]>segPara.boneThresh || maskA[k2+sizex]==maskStruct.mask_corticalBone || maskA[k2+sizex]==maskStruct.mask_spongyBone)
					&& (imgA[k2+sizex2]>segPara.boneThresh || maskA[k2+sizex2]==maskStruct.mask_corticalBone || maskA[k2+sizex2]==maskStruct.mask_spongyBone))
				{
					count++;
					uy += y;
					break;
				}
			}
		}
		if(count>0) uy/=count;
		else uy = segInfo.diskBound1[z].y;

		// check the lower border
		ly=0;
		count=0;
		for(x=(int)diskCenter.x-2; x<=(int)diskCenter.x+2; x++)
		{
			for(y=(diskCenter.y+segInfo.diskBound2[z].y)/2; y<segInfo.diskBound2[z].y+3; y++)
			{
				k2 = k+y*sizex+x;
				if(maskA[k2]==maskStruct.mask_spinalCord)
				{
					count++;
					ly += y;
					break;
				}
			}
		}
		if(count>0) ly/=count;
		else ly = segInfo.diskBound2[z].y;
		diskCenter.y = (uy+ly)/2;
		diskRadius.y = (ly-uy)/2;


		// adjust the x coordinate of diskCenter
		int rx, lx;
		// check the right border
		rx=0;
		count=0;
		for(y=(int)diskCenter.y-2; y<=(int)diskCenter.y+2; y++)
		{
			for(x=segInfo.spineBound2[z].x; x>diskCenter.x; x--)
			{
				k2 = k+y*sizex+x;
				if((imgA[k2]>segPara.boneThresh || maskA[k2]==maskStruct.mask_corticalBone || maskA[k2]==maskStruct.mask_spongyBone) 
					&& (imgA[k2-1]>segPara.boneThresh || maskA[k2-1]==maskStruct.mask_corticalBone || maskA[k2-1]==maskStruct.mask_spongyBone) 
					&& (imgA[k2-2]>segPara.boneThresh || maskA[k2-2]==maskStruct.mask_corticalBone || maskA[k2-2]==maskStruct.mask_spongyBone) )
				{
					count++;
					rx += x;
					break;
				}
			}
		}
		if(count>0) rx/=count;
		else rx = segInfo.diskBound2[z].x;

		// check the left border
		lx=0;
		count=0;
		for(y=(int)diskCenter.y-2; y<=(int)diskCenter.y+2; y++)
		{
			for(x=segInfo.spineBound1[z].x; x<diskCenter.x; x++)
			{
				k2 = k+y*sizex+x;
				if((imgA[k2]>segPara.boneThresh || maskA[k2]==maskStruct.mask_corticalBone || maskA[k2]==maskStruct.mask_spongyBone)
					&& (imgA[k2+1]>segPara.boneThresh || maskA[k2+1]==maskStruct.mask_corticalBone || maskA[k2+1]==maskStruct.mask_spongyBone)
					&& (imgA[k2+2]>segPara.boneThresh || maskA[k2+2]==maskStruct.mask_corticalBone || maskA[k2+2]==maskStruct.mask_spongyBone))
				{
					count++;
					lx += x;
					break;
				}
			}
		}
		if(count>0) lx/=count;
		else lx = segInfo.diskBound1[z].x;

		diskCenter.x = (rx+lx)/2;
		diskRadius.x = (rx-lx)/2;


		// get the lower muscle and fat ROI
		Vec2 mroiCenter, mroiLeft, mroiRight;
		vec2DynArray lowerBorder;
		mroiCenter = (segInfo.sprocessBound1[z]+segInfo.sprocessBound2[z]*2)/3;
///		mroiRight = mroiCenter+Vec2(diskRadius.x*2.5,0);
///		mroiLeft = mroiCenter+Vec2(-diskRadius.x*2.5,0);
		mroiRight = mroiCenter+Vec2(diskRadius.x*2,0);	// new smaller rois
		mroiLeft = mroiCenter+Vec2(-diskRadius.x*2,0);

		for(x=(int)mroiLeft.x; x<=(int)mroiRight.x; x++)
		{
			for(y=segInfo.sprocessBound2[z].y; y<sizey; y++)
			{
				if(imgAL[y*sizex+x]==1024)
				{
					lowerBorder.Add(Vec2(x,y-10));
					break;
				}
			}
		}

		// fit smooth spline curves
		int bernstein_power = 5, piece_size=10;
		doubleDynArray xt, yt, yt2;

		xt.SetSize(lowerBorder.GetSize());
		yt.SetSize(lowerBorder.GetSize());
		yt2.SetSize(lowerBorder.GetSize());

		for(k2=0; k2<lowerBorder.GetSize(); k2++)
		{
			xt[k2] = lowerBorder[k2].x;
			yt[k2] = lowerBorder[k2].y;
		}
		PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
		for(k2=0; k2<lowerBorder.GetSize(); k2++)
		{
			lowerBorder[k2].y = yt2[k2];
		}

		// erase all mask 
		k=z*sizexy;
		for(k2=0; k2<sizexy; k2++)
			maskA[k+k2] = maskStruct.mask_body;

		// diskROI
		// only use 2/3 of the disk as the roi
///		diskRadius *= 0.66;
		diskRadius *= 0.5;		// a new radius
		diskRadius *= 2;	// need diameter in the model

		// compute histogram assigned the mask
		float minG, maxG, stepG;
		int ind;
		short hu;
		vec2DynArray scan;
		CIS_2D_ROI_Ellipse *ell;
		ell = new CIS_2D_ROI_Ellipse(diskCenter, diskRadius);
		ell->GetEntireScan(scan);
		delete ell;

		minG = 100000;
		maxG = 0;
		for(k2=0; k2<scan.GetSize(); k2++)
		{
			ind = k+(int)scan[k2].x+(int)scan[k2].y*sizex;
			if(imgA[ind]<minG) minG=imgA[ind];
			if(imgA[ind]>maxG) maxG=imgA[ind];
			maskA[ind] = maskStruct.mask_boneMetasis;
		}
///		stepG = (maxG-minG)/255;
///		trab_ct[p+1][0] = minG-1024;	// convert to HU
///		trab_feq[p+1][0]=0;
		stepG = (upperEnd_tr-lowerEnd_tr)/binNumber;
		trab_ct[p+1][0] = lowerEnd_tr;	// convert to HU
		trab_feq[p+1][0]=0;
		for(ind=1; ind<binNumber; ind++) 
		{
			trab_ct[p+1][ind] = trab_ct[p+1][ind-1]+stepG;
			trab_feq[p+1][ind]=0;
		}
		for(k2=0; k2<scan.GetSize(); k2++)
		{
			ind = k+(int)scan[k2].x+(int)scan[k2].y*sizex;
			hu = imgA[ind]-1024;
			if(hu>lowerEnd_tr && hu<upperEnd_tr) 
			{
				trab_feq[p+1][(int)((hu-lowerEnd_tr)/stepG)]++;
			}

			// one slice up
			ind = k+sizexy+(int)scan[k2].x+(int)scan[k2].y*sizex;
			hu = imgA[ind]-1024;
			if(hu>lowerEnd_tr && hu<upperEnd_tr) 
			{
				trab_feq[p+1][(int)((hu-lowerEnd_tr)/stepG)]++;
			}

			// one slice down
			ind = k-sizexy+(int)scan[k2].x+(int)scan[k2].y*sizex;
			hu = imgA[ind]-1024;
			if(hu>lowerEnd_tr && hu<upperEnd_tr) 
			{
				trab_feq[p+1][(int)((hu-lowerEnd_tr)/stepG)]++;
			}

		}

		// muscle ROI
		lowerBorder.Add(mroiRight);
		lowerBorder.Add(mroiLeft);

		CIS_2D_ROI_Polygon *polyROI;
		polyROI = new CIS_2D_ROI_Polygon(lowerBorder);
		polyROI->GetEntireScan(scan);
		delete polyROI;

		minG = 100000;
		maxG = 0;
		k=z*sizexy;
		for(k2=0; k2<scan.GetSize(); k2++)
		{
			ind = k+(int)scan[k2].x+(int)scan[k2].y*sizex;
			if(imgA[ind]<minG) minG=imgA[ind];
			if(imgA[ind]>maxG) maxG=imgA[ind];
			maskA[ind] = maskStruct.mask_boneMetasis;
		}
///		stepG = (maxG-minG)/255;
///		fat_ct[p+1][0] = minG-1024;		// convert to HU
///		fat_feq[p+1][0]=0;
		stepG = (upperEnd_fat-lowerEnd_fat)/binNumber;
		fat_ct[p+1][0] = lowerEnd_fat;	// convert to HU
		fat_feq[p+1][0]=0;

		for(ind=1; ind<binNumber; ind++) 
		{
			fat_ct[p+1][ind] = fat_ct[p+1][ind-1]+stepG;
			fat_feq[p+1][ind]=0;
		}
		for(k2=0; k2<scan.GetSize(); k2++)
		{
			ind = k+(int)scan[k2].x+(int)scan[k2].y*sizex;
			hu = imgA[ind]-1024;
			if(hu>lowerEnd_fat && hu<upperEnd_fat)
			{
				fat_feq[p+1][(int)((hu-lowerEnd_fat)/stepG)]++;
			}

			// one slice up
			ind = k+sizexy+(int)scan[k2].x+(int)scan[k2].y*sizex;
			hu = imgA[ind]-1024;
			if(hu>lowerEnd_fat && hu<upperEnd_fat)
			{
				fat_feq[p+1][(int)((hu-lowerEnd_fat)/stepG)]++;
			}

			// one slice down
			ind = k-sizexy+(int)scan[k2].x+(int)scan[k2].y*sizex;
			hu = imgA[ind]-1024;
			if(hu>lowerEnd_fat && hu<upperEnd_fat)
			{
				fat_feq[p+1][(int)((hu-lowerEnd_fat)/stepG)]++;
			}
		}
	}	// for p

	// output Nicoli's histogram into a file
	FILE *e_fp;
	if((e_fp = fopen(export_fn, "w"))!=NULL)
	{
		// print the header
		fprintf(e_fp,"l1_trabctnum,l1_trabfreq,l1_fmctnum,l1_fmfreq,l2_trabctnum,l2_trabfreq,l2_fmctnum,l2_fmfreq,\n");

		for(int i=0; i<binNumber; i++)
		{
			for(p=0; p<vertNum; p++)
			{
				fprintf(e_fp, "%.3f,%d,%.3f,%d,",trab_ct[p][i], trab_feq[p][i], fat_ct[p][i], fat_feq[p][i]);
			}
			fprintf(e_fp, "\n");
		}
		fclose(e_fp);
	}


	// release the memory
	for(p=0; p<vertNum; p++) 
	{
		delete trab_ct[p];
		delete fat_ct[p];
		delete trab_feq[p];
		delete fat_feq[p];
	}
	delete trab_ct;
	delete trab_feq;
	delete fat_ct;
	delete fat_feq;

	return CIS_OK;
}



int NIH_SpineSegmentation_BMD_ROI_L1L2_3D(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   VertebraStruct2D *vertebraTemplate, 
									   doubleDynArray &spineWidthLeft, doubleDynArray &spineWidthRight, doubleDynArray &spineWidthUp, doubleDynArray &spineWidthDown, 
									   intDynArray &pedicleValley, intDynArray &backValley, intDynArray &roiSlices,
									   bool debugMode,
									   const char* export_fn,
									   int l1_slice, CIS_Array_Image3D_short *l1_img,
									   int l2_slice, CIS_Array_Image3D_short *l2_img)
{
	int sizex, sizey, sizez, sizexy, sizex2;
	short *maskA, *imgA, *imgAL;
	int x, y, z, k, k2;
	float px, py, pz;

	sizex = maskImg3D->Num_Cols();
	sizey = maskImg3D->Num_Rows();
	sizez = maskImg3D->Num_Levels();
	px = img3D->Get_Pixel_SizeX();
	py = img3D->Get_Pixel_SizeY();
	pz = img3D->Get_Pixel_SizeZ();
	sizexy = sizex*sizey;
	sizex2 =sizex*2;

	maskA = maskImg3D->GetArray();
	imgA = img3D->GetArray();

	CIS_Array_Image2D_short *binImg2D;
	short *binA2D;

	binImg2D = new CIS_Array_Image2D_short(sizex, sizey);
	binA2D = binImg2D->GetArray();

	int vertNum = 2;		// only l1 and l2
	// variables to store Nicoli's histogram
	int binNumber, upperEnd_tr, lowerEnd_tr, upperEnd_fat, lowerEnd_fat;

	float **trab_ct, **fat_ct;
	int **trab_feq, **fat_feq;
	int p;

///	binNumber = 100;
	binNumber = 500;
	upperEnd_tr = 500;
	lowerEnd_tr = 0;
	upperEnd_fat = 300;
	lowerEnd_fat = -200;

	trab_ct = new float*[vertNum];
	for(p=0; p<vertNum; p++) trab_ct[p] = new float[binNumber];
	fat_ct = new float*[vertNum];
	for(p=0; p<vertNum; p++) fat_ct[p] = new float[binNumber];
	trab_feq = new int*[vertNum];
	for(p=0; p<vertNum; p++) trab_feq[p] = new int[binNumber];
	fat_feq = new int*[vertNum];
	for(p=0; p<vertNum; p++) fat_feq[p] = new int[binNumber];

	imgAL = l1_img->GetArray();
///	for(k2=0; k2<sizexy; k2++) imgAL[k2]+=1024;
	int k3;
	for(z=0; z<sizez; z++)
	{
		k3=z*sizexy;
		for(y=0; y<sizey/2; y++)
			for(x=0; x<sizex; x++)
			{
				k2= y*sizex+x;
				k = (sizey-1-y)*sizex+x;
				p=imgAL[k3+k2]; imgAL[k3+k2]=imgAL[k3+k]; imgAL[k3+k]=p;
			}
	}

/*	k=l1_slice*sizexy;
	for(k2=0; k2<sizexy; k2++) imgA[k+k2]=imgAL[k2];
	imgAL = l2_img->GetArray();
	for(k2=0; k2<sizexy; k2++) imgAL[k2]+=1024;
	for(y=0; y<sizey/2; y++)
		for(x=0; x<sizex; x++)
		{
			k2= y*sizex+x;
			k = (sizey-1-y)*sizex+x;
			p=imgAL[k2]; imgAL[k2]=imgAL[k]; imgAL[k]=p;
		}
	k=l2_slice*sizexy;
	for(k2=0; k2<sizexy; k2++) imgA[k+k2]=imgAL[k2];
*/

	roiSlices.SetSize(0);
	for(p=-1; p<1; p++)
	{
		Vec2 diskCenter, diskRadius;
		int a, count;
		int sz, mz, ez;

		// get the middle slice
		if(p==-1) 
		{
			mz = l1_slice;
			imgAL = l1_img->GetArray();
		}
		else 
		{
			mz = l2_slice;
			imgAL = l2_img->GetArray();
		}

		if(segInfo.diskBound1[mz].x==-1 || 
			segInfo.diskBound1[mz].x>segInfo.diskBound2[mz].x-2 || 
			segInfo.diskBound1[mz].y>segInfo.diskBound2[mz].y-2 ||
			segInfo.sprocessBound1[mz].y>segInfo.sprocessBound2[mz].y-2)
		{
			printf("Unable to generate ROI on slice %d\n", mz);
			continue;	// bounding box not right
		}

		k = mz*sizexy;
		roiSlices.Add(mz);

		// init center
		diskCenter = (segInfo.diskBound1[mz]+segInfo.diskBound2[mz])/2;
		x = (int)diskCenter.x;
		y = (int)diskCenter.y;
		sz=mz;
		for(z=mz-1; z>=0; z--)
		{
			if(imgAL[z*sizexy+y*sizey+x]==-32768) break;
		}
		sz = z+1;

		ez=mz;
		for(z=mz+1; z<sizez-1; z++)
		{
			if(imgAL[z*sizexy+y*sizey+x]==-32768) break;
		}
		ez = z-1;

		// initialize histogram
		int ind, ck;
		short hu;
		float stepG = (upperEnd_tr-lowerEnd_tr)/binNumber;
		trab_ct[p+1][0] = lowerEnd_tr;	// convert to HU
		trab_feq[p+1][0]=0;
		for(ind=1; ind<binNumber; ind++) 
		{
			trab_ct[p+1][ind] = trab_ct[p+1][ind-1]+stepG;
			trab_feq[p+1][ind]=0;
		}

		stepG = (upperEnd_fat-lowerEnd_fat)/binNumber;
		fat_ct[p+1][0] = lowerEnd_fat;	// convert to HU
		fat_feq[p+1][0]=0;
		for(ind=1; ind<binNumber; ind++) 
		{
			fat_ct[p+1][ind] = fat_ct[p+1][ind-1]+stepG;
			fat_feq[p+1][ind]=0;
		}

		for(z=sz+2; z<=ez-2; z++)
		{
			if(z<=0) continue;

			diskCenter = (segInfo.diskBound1[z]+segInfo.diskBound2[z])/2;
			// adjust the y coordinate of diskCenter
			int uy, ly;
			// check the upper border
			uy=0;
			count=0;
			for(x=(int)diskCenter.x-2; x<=(int)diskCenter.x+2; x++)
			{
				for(y=segInfo.spineBound1[z].y; y<diskCenter.y; y++)
				{
					k2 = k+y*sizex+x;
					if((imgA[k2]>segPara.boneThresh || maskA[k2]==maskStruct.mask_corticalBone || maskA[k2]==maskStruct.mask_spongyBone)
						&& (imgA[k2+sizex]>segPara.boneThresh || maskA[k2+sizex]==maskStruct.mask_corticalBone || maskA[k2+sizex]==maskStruct.mask_spongyBone)
						&& (imgA[k2+sizex2]>segPara.boneThresh || maskA[k2+sizex2]==maskStruct.mask_corticalBone || maskA[k2+sizex2]==maskStruct.mask_spongyBone))
					{
						count++;
						uy += y;
						break;
					}
				}
			}
			if(count>0) uy/=count;
			else uy = segInfo.diskBound1[z].y;

			// check the lower border
			ly=0;
			count=0;
			for(x=(int)diskCenter.x-2; x<=(int)diskCenter.x+2; x++)
			{
				for(y=(diskCenter.y+segInfo.diskBound2[z].y)/2; y<segInfo.diskBound2[z].y+3; y++)
				{
					k2 = k+y*sizex+x;
					if(maskA[k2]==maskStruct.mask_spinalCord)
					{
						count++;
						ly += y;
						break;
					}
				}
			}
			if(count>0) ly/=count;
			else ly = segInfo.diskBound2[z].y;
			diskCenter.y = (uy+ly)/2;
			diskRadius.y = (ly-uy)/2;


			// adjust the x coordinate of diskCenter
			int rx, lx;
			// check the right border
			rx=0;
			count=0;
			for(y=(int)diskCenter.y-2; y<=(int)diskCenter.y+2; y++)
			{
				for(x=segInfo.spineBound2[z].x; x>diskCenter.x; x--)
				{
					k2 = k+y*sizex+x;
					if((imgA[k2]>segPara.boneThresh || maskA[k2]==maskStruct.mask_corticalBone || maskA[k2]==maskStruct.mask_spongyBone) 
						&& (imgA[k2-1]>segPara.boneThresh || maskA[k2-1]==maskStruct.mask_corticalBone || maskA[k2-1]==maskStruct.mask_spongyBone) 
						&& (imgA[k2-2]>segPara.boneThresh || maskA[k2-2]==maskStruct.mask_corticalBone || maskA[k2-2]==maskStruct.mask_spongyBone) )
					{
						count++;
						rx += x;
						break;
					}
				}
			}
			if(count>0) rx/=count;
			else rx = segInfo.diskBound2[z].x;

			// check the left border
			lx=0;
			count=0;
			for(y=(int)diskCenter.y-2; y<=(int)diskCenter.y+2; y++)
			{
				for(x=segInfo.spineBound1[z].x; x<diskCenter.x; x++)
				{
					k2 = k+y*sizex+x;
					if((imgA[k2]>segPara.boneThresh || maskA[k2]==maskStruct.mask_corticalBone || maskA[k2]==maskStruct.mask_spongyBone)
						&& (imgA[k2+1]>segPara.boneThresh || maskA[k2+1]==maskStruct.mask_corticalBone || maskA[k2+1]==maskStruct.mask_spongyBone)
						&& (imgA[k2+2]>segPara.boneThresh || maskA[k2+2]==maskStruct.mask_corticalBone || maskA[k2+2]==maskStruct.mask_spongyBone))
					{
						count++;
						lx += x;
						break;
					}
				}
			}
			if(count>0) lx/=count;
			else lx = segInfo.diskBound1[z].x;

			diskCenter.x = (rx+lx)/2;
			diskRadius.x = (rx-lx)/2;

			// get the lower muscle and fat ROI
			Vec2 mroiCenter, mroiLeft, mroiRight;
			vec2DynArray lowerBorder;
			mroiCenter = (segInfo.sprocessBound1[z]+segInfo.sprocessBound2[z]*2)/3;
			mroiRight = mroiCenter+Vec2(diskRadius.x*2,0);	// new smaller rois
			mroiLeft = mroiCenter+Vec2(-diskRadius.x*2,0);

			k=z*sizexy;
			for(x=(int)mroiLeft.x; x<=(int)mroiRight.x; x++)
			{
				for(y=segInfo.sprocessBound2[z].y; y<sizey; y++)
				{
					if(imgAL[k+y*sizex+x]==-32768)
					{
						lowerBorder.Add(Vec2(x,y-10));
						break;
					}
				}
			}

			// fit smooth spline curves
			int bernstein_power = 5, piece_size=10;
			doubleDynArray xt, yt, yt2;

			xt.SetSize(lowerBorder.GetSize());
			yt.SetSize(lowerBorder.GetSize());
			yt2.SetSize(lowerBorder.GetSize());

			for(k2=0; k2<lowerBorder.GetSize(); k2++)
			{
				xt[k2] = lowerBorder[k2].x;
				yt[k2] = lowerBorder[k2].y;
			}
			PiecewiseBernsteinSmoothing(xt, yt, yt2, bernstein_power, piece_size);
			for(k2=0; k2<lowerBorder.GetSize(); k2++)
			{
				lowerBorder[k2].y = yt2[k2];
			}

			// get bone ROI
			k = z*sizexy;
			stepG = (upperEnd_tr-lowerEnd_tr)/binNumber;

			// erosion to get trabecular bone
			//
			int erosion_it=6;
			for(ck=0; ck<sizexy; ck++) binA2D[ck]=0;
			for(y=segInfo.diskBound1[z].y; y<=segInfo.diskBound2[z].y; y++)
			{
				for(x=segInfo.diskBound1[z].x; x<segInfo.diskBound2[z].x; x++)
				{
					ck = k+ y*sizey+x;
					if(maskA[ck]==maskStruct.mask_corticalBone || maskA[ck]==maskStruct.mask_spongyBone)
					{
						binA2D[y*sizey+x]=1;
					}
				}
			}

			CIS_IPA_Erode(binImg2D, erosion_it, false);

			for(y=segInfo.diskBound1[z].y; y<=segInfo.diskBound2[z].y; y++)
			{
				for(x=segInfo.diskBound1[z].x; x<segInfo.diskBound2[z].x; x++)
				{
					ck = k+ y*sizey+x;
///					if(maskA[ck]==maskStruct.mask_corticalBone || maskA[ck]==maskStruct.mask_spongyBone)
					if(binA2D[y*sizey+x]==1)
					{
						hu = imgA[ck]-1024;
						if(hu>lowerEnd_tr && hu<upperEnd_tr)
						{
							trab_feq[p+1][(int)((hu-lowerEnd_tr)/stepG)]++;
						}
						// include these pixels in 3D ROI of the vertebra
						maskA[ck] = maskStruct.mask_boneMetasis;
					}	// if maskA
				}
			}

			// compute histogram assigned the mask
			float minG, maxG, stepG;
			int ind;
			short hu;
			vec2DynArray scan;

			// muscle ROI
			lowerBorder.Add(mroiRight);
			lowerBorder.Add(mroiLeft);

			CIS_2D_ROI_Polygon *polyROI;
			polyROI = new CIS_2D_ROI_Polygon(lowerBorder);
			polyROI->GetEntireScan(scan);
			delete polyROI;

			k=z*sizexy;
			stepG = (upperEnd_fat-lowerEnd_fat)/binNumber;

			for(k2=0; k2<scan.GetSize(); k2++)
			{
				ind = k+(int)scan[k2].x+(int)scan[k2].y*sizex;
				if(maskA[ind]==maskStruct.mask_corticalBone || maskA[ind]==maskStruct.mask_spongyBone) continue;
				maskA[ind] = maskStruct.mask_boneMetasis;
				hu = imgA[ind]-1024;
				if(hu>lowerEnd_fat && hu<upperEnd_fat)
				{
					fat_feq[p+1][(int)((hu-lowerEnd_fat)/stepG)]++;
				}
			}
		}	// for z

	}	// for p

	// output Nicoli's histogram into a file
	FILE *e_fp;
	if((e_fp = fopen(export_fn, "w"))!=NULL)
	{
		// print the header
		fprintf(e_fp,"l1_trabctnum,l1_trabfreq,l1_fmctnum,l1_fmfreq,l2_trabctnum,l2_trabfreq,l2_fmctnum,l2_fmfreq,\n");

		for(int i=0; i<binNumber; i++)
		{
			for(p=0; p<vertNum; p++)
			{
				fprintf(e_fp, "%.3f,%d,%.3f,%d,",trab_ct[p][i], trab_feq[p][i], fat_ct[p][i], fat_feq[p][i]);
			}
			fprintf(e_fp, "\n");
		}

		// compute mean and standard deviation of the bone ROI using histogram
		for(p=0; p<vertNum; p++)
		{
			long sum=0;
			long sum2=0;
			long count=0;
			for(int i=0; i<binNumber; i++)
			{
				count += trab_feq[p][i];
				sum += trab_feq[p][i]*trab_ct[p][i];
				sum2 += trab_feq[p][i]*(trab_ct[p][i]*trab_ct[p][i]);
			}

			if(count>0) sum/=count;
			long dev = (long)sqrt((double)(sum2/count-sum*sum));

			// export
			fprintf(e_fp, "Vert T%d, %d, %d\n",p+1, sum, dev);
		}

		fclose(e_fp);
	}


	// release the memory
	for(p=0; p<vertNum; p++) 
	{
		delete trab_ct[p];
		delete fat_ct[p];
		delete trab_feq[p];
		delete fat_feq[p];
	}
	delete trab_ct;
	delete trab_feq;
	delete fat_ct;
	delete fat_feq;

	delete binImg2D;
	return CIS_OK;
}






//
// automatically detect the rib cage based on the extracted and partitioned spinal column
//
int NIH_SpineSegmentation_Rib_Detection(CIS_Array_Image3D_short *img3D,
									   CIS_Array_Image3D_short *maskImg3D,
									   SpineSegmentationMaskStructure &maskStruct,
									   SpineSegmentationParameter &segPara,
									   SpineSegmentationInfo &segInfo,
									   VertebraStruct2D *vertebraTemplate, 
									   doubleDynArray &spineWidthLeft, doubleDynArray &spineWidthRight, doubleDynArray &spineWidthUp, doubleDynArray &spineWidthDown, 
									   intDynArray &pedicleValley, intDynArray &backValley, 
									   bool debugMode)
{
	if(segInfo.bound1.z<0) return CIS_ERROR;
	int sizex, sizey, sizez, sizexy, sizex2;
	short *maskA, *imgA;
	int x, y, z, k, k2, z1, z2;
	float px, py, pz;

	sizex = maskImg3D->Num_Cols();
	sizey = maskImg3D->Num_Rows();
	sizez = maskImg3D->Num_Levels();

	if(sizez>500) sizez=500;		// memory issue;

	px = img3D->Get_Pixel_SizeX();
	py = img3D->Get_Pixel_SizeY();
	pz = img3D->Get_Pixel_SizeZ();
	sizexy = sizex*sizey;
	sizex2 =sizex*2;

	maskA = maskImg3D->GetArray();
	imgA = img3D->GetArray();

	z1 = segInfo.bound1.z;
	if(segInfo.sacrum_start>0) 
	{
		z2=(int)(segInfo.sacrum_start-(30/pz)*3);	// the range for rib searching should be well above sacrum
		if(segInfo.sacrum_start+30<sizez) sizez = segInfo.sacrum_start+30;
	}
	else 
	{
		z2=segInfo.bound2.z-1;
		if(segInfo.bound2.z<sizez) sizez = segInfo.bound2.z;
	}

	// first creat a binary mask for all rib candidate pixels
	CIS_Array_Image3D_uchar *binImg;
	CIS_Array_Image3D_short *blobImg;
	unsigned char *binArray;
	short *blobArray;

	binImg = new CIS_Array_Image3D_uchar(sizex, sizey, sizez);
	blobImg = new CIS_Array_Image3D_short(sizex, sizey, sizez);
	binArray = binImg->GetArray();
	blobArray = blobImg->GetArray();

	CIS_Array_Image2D_short *binImg2D, *blobImg2D;
	short *binArray2D, *blobArray2D;
	intDynArray blobRank2D;
	intVec2DynArray blobCentroid2D;

	binImg2D = new CIS_Array_Image2D_short(sizex, sizey);
	blobImg2D = new CIS_Array_Image2D_short(sizex, sizey);
	binArray2D = binImg2D->GetArray();
	blobArray2D = blobImg2D->GetArray();

	for(z=z1; z<=z2; z++)
	{
		// invalid bounding boxes for vertebra
		if(segInfo.spineBound1[z].x<0 || segInfo.spineBound2[z].x<0 || 
			segInfo.spineBound1[z].y<0 || segInfo.spineBound2[z].y<0) continue;
		if(segInfo.diskBound1[z].x<0 || segInfo.diskBound2[z].x<0 || 
			segInfo.diskBound1[z].y<0 || segInfo.diskBound2[z].y<0) continue;
///		if(segInfo.sprocessBound1[z].x<0 || segInfo.sprocessBound2[z].x<0 || 
///			segInfo.sprocessBound1[z].y<0 || segInfo.sprocessBound2[z].y<0) continue;
		if(segInfo.spineBound1[z].x>=segInfo.spineBound2[z].x || 
			segInfo.spineBound1[z].y>=segInfo.spineBound2[z].y) continue;
		if(segInfo.diskBound1[z].x>=segInfo.diskBound2[z].x || 
			segInfo.diskBound1[z].y>=segInfo.diskBound2[z].y) continue;
///		if(segInfo.sprocessBound1[z].x>=segInfo.sprocessBound2[z].x || 
///			segInfo.sprocessBound1[z].y>=segInfo.sprocessBound2[z].y) continue;

		// get any candidate pixels
		k = z*sizexy;
		for(y=0; y>=0 && y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++)
			{
				if(imgA[k]>segPara.boneThresh && maskA[k]!=maskStruct.mask_corticalBone && maskA[k]!=maskStruct.mask_spongyBone 
					&& maskA[k]!=maskStruct.mask_vertebralDisk && maskA[k]!=maskStruct.mask_spinalCord)
				{
					binArray[k] = 1;
//					maskA[k] = maskStruct.mask_rib;
				}
				else binArray[k] = 0;
			}
		}
		// wipe off anything right above the vertebra body
		for(y=0; y>=0 && y<segInfo.diskBound1[z].y; y++)
		{
			for(x=segInfo.spineBound1[z].x, k=z*sizexy+y*sizex+x; x<=segInfo.spineBound2[z].x; x++, k++)
			{
				binArray[k] = 0;
			}
		}

		// wipe off anything below the bounding box, with a margin for error
		for(y=(segInfo.spineBound2[z].y*5+sizey)/6; y>=0 && y<sizey; y++)
		{
			for(x=0, k=z*sizexy+y*sizex+x; x<sizex; x++, k++)
			{
				binArray[k] = 0;
			}
		}

		// wipe off the region below the spinal process
/*		if(segInfo.sprocessBound2[z].y>0)
		{
			for(y=segInfo.sprocessBound2[z].y; y<sizey; y++)
			{
				for(x=0, k=z*sizexy+y*sizex; x<sizex; x++, k++)
				{
					binArray[k] = 0;
				}
			}
		}
*/
		// wipe off the region below the spinal cord
		for(y=segInfo.cordBound2[z].y; y>=0 && y<sizey; y++)
		{
			for(x=(segInfo.diskBound1[z].x*3+segInfo.diskBound2[z].x)/4, k=z*sizexy+y*sizex+x; 
				x<=(segInfo.diskBound1[z].x+segInfo.diskBound2[z].x*3)/4; x++, k++)
			{
				binArray[k] = 0;
			}
		}

		k = z*sizexy;
		for(y=0, k2=0; y<sizey; y++)
			for(x=0; x<sizex; x++, k++, k2++)
				binArray2D[k2]=binArray[k];

		// erosion in 2D to remove small blobs, which are more like noise
		CIS_IPA_Erode(binImg2D, 1, true);
		CIS_IPA_Dilate(binImg2D, 1, true);

		// 2D blob analysis to remove big blob, which is more likely fluid pockets
		NIH_Algo_Blob_Labelling(binImg2D, blobImg2D, blobRank2D, blobCentroid2D, false);

		k = z*sizexy;
		short cb;
		for(y=0, k2=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++, k2++)
			{
				cb = blobArray2D[k2];
				if(cb>0 && blobRank2D[cb-1]<2000 && blobRank2D[cb-1]>9) binArray[k]=1;
				else binArray[k]=0;
			}
		}
	}	// for z

	// may need expansion to recover some regions misclassified as vertebra
	// recover misclassified vertebra in the region between spineBound and diskBound
	for(z=z1+1; z<=z2; z++)
	{
		// invalid bounding boxes for vertebra
		if(segInfo.spineBound1[z].x<0 || segInfo.spineBound2[z].x<0 || 
			segInfo.spineBound1[z].y<0 || segInfo.spineBound2[z].y<0) continue;
		if(segInfo.diskBound1[z].x<0 || segInfo.diskBound2[z].x<0 || 
			segInfo.diskBound1[z].y<0 || segInfo.diskBound2[z].y<0) continue;
		if(segInfo.spineBound1[z].x>=segInfo.spineBound2[z].x || 
			segInfo.spineBound1[z].y>=segInfo.spineBound2[z].y) continue;
		if(segInfo.diskBound1[z].x>=segInfo.diskBound2[z].x || 
			segInfo.diskBound1[z].y>=segInfo.diskBound2[z].y) continue;

		// get any candidate pixels
		k = z*sizexy;
		for(y=segInfo.diskBound2[z].y; y<segInfo.spineBound2[z].y; y++)
		{
			for(x=segInfo.spineBound1[z].x,k2=k+y*sizex+x; x<=(segInfo.diskBound1[z].x*7+segInfo.diskBound2[z].x)/8; x++, k2++)
			{
				if(imgA[k2]>segPara.boneThresh && binArray[k2]==0)
				{
					if(binArray[k2-sizexy]==1 && maskA[k2-sizexy]!=maskStruct.mask_corticalBone
						&& maskA[k2-sizexy]!=maskStruct.mask_vertebralDisk
						&& maskA[k2-sizexy]!=maskStruct.mask_spongyBone) binArray[k2] = 1;
				}
			}
		}
		for(y=segInfo.spineBound1[z].y; y<=segInfo.diskBound2[z].y; y++)
		{
			for(x=segInfo.spineBound1[z].x,k2=k+y*sizex+x; x<=segInfo.diskBound1[z].x; x++, k2++)
			{
				if(imgA[k2]>segPara.boneThresh && binArray[k2]==0)
				{
					if(binArray[k2-sizexy]==1 && maskA[k2-sizexy]!=maskStruct.mask_corticalBone
						&& maskA[k2-sizexy]!=maskStruct.mask_vertebralDisk
						&& maskA[k2-sizexy]!=maskStruct.mask_spongyBone) binArray[k2] = 1;
				}
			}
		}
		for(y=segInfo.diskBound2[z].y; y<segInfo.spineBound2[z].y; y++)
		{
			for(x=(segInfo.diskBound1[z].x+segInfo.diskBound2[z].x*7)/8,k2=k+y*sizex+x; x<=segInfo.spineBound2[z].x; x++, k2++)
			{
				if(imgA[k2]>segPara.boneThresh && binArray[k2]==0)
				{
					if(binArray[k2-sizexy]==1 && maskA[k2-sizexy]!=maskStruct.mask_corticalBone
						&& maskA[k2-sizexy]!=maskStruct.mask_vertebralDisk
						&& maskA[k2-sizexy]!=maskStruct.mask_spongyBone) binArray[k2] = 1;
				}
			}
		}
		for(y=segInfo.spineBound1[z].y; y<=segInfo.diskBound2[z].y; y++)
		{
			for(x=segInfo.diskBound2[z].x,k2=k+y*sizex+x; x<=segInfo.spineBound2[z].x; x++, k2++)
			{
				if(imgA[k2]>segPara.boneThresh && binArray[k2]==0)
				{
					if(binArray[k2-sizexy]==1 && maskA[k2-sizexy]!=maskStruct.mask_corticalBone
						&& maskA[k2-sizexy]!=maskStruct.mask_vertebralDisk
						&& maskA[k2-sizexy]!=maskStruct.mask_spongyBone) binArray[k2] = 1;
				}
			}
		}
	}	// for z

	// reverse direction
	for(z=z2-1; z>=z1; z--)
	{
		// invalid bounding boxes for vertebra
		if(segInfo.spineBound1[z].x<0 || segInfo.spineBound2[z].x<0 || 
			segInfo.spineBound1[z].y<0 || segInfo.spineBound2[z].y<0) continue;
		if(segInfo.diskBound1[z].x<0 || segInfo.diskBound2[z].x<0 || 
			segInfo.diskBound1[z].y<0 || segInfo.diskBound2[z].y<0) continue;
		if(segInfo.spineBound1[z].x>=segInfo.spineBound2[z].x || 
			segInfo.spineBound1[z].y>=segInfo.spineBound2[z].y) continue;
		if(segInfo.diskBound1[z].x>=segInfo.diskBound2[z].x || 
			segInfo.diskBound1[z].y>=segInfo.diskBound2[z].y) continue;

		// get any candidate pixels
		k = z*sizexy;
		for(y=segInfo.diskBound2[z].y; y<segInfo.spineBound2[z].y; y++)
		{
			for(x=segInfo.spineBound1[z].x,k2=k+y*sizex+x; x<=(segInfo.diskBound1[z].x*7+segInfo.diskBound2[z].x)/8; x++, k2++)
			{
				if(imgA[k2]>segPara.boneThresh && binArray[k2]==0)
				{
					if(binArray[k2+sizexy]==1 && maskA[k2+sizexy]!=maskStruct.mask_corticalBone
						&& maskA[k2+sizexy]!=maskStruct.mask_vertebralDisk
						&& maskA[k2+sizexy]!=maskStruct.mask_spongyBone) binArray[k2] = 1;
				}
			}
		}
		for(y=segInfo.spineBound1[z].y; y<=segInfo.diskBound2[z].y; y++)
		{
			for(x=segInfo.spineBound1[z].x,k2=k+y*sizex+x; x<=segInfo.diskBound1[z].x; x++, k2++)
			{
				if(imgA[k2]>segPara.boneThresh && binArray[k2]==0)
				{
					if(binArray[k2+sizexy]==1 && maskA[k2+sizexy]!=maskStruct.mask_corticalBone
						&& maskA[k2+sizexy]!=maskStruct.mask_vertebralDisk
						&& maskA[k2+sizexy]!=maskStruct.mask_spongyBone) binArray[k2] = 1;
				}
			}
		}
		for(y=segInfo.diskBound2[z].y; y<segInfo.spineBound2[z].y; y++)
		{
			for(x=(segInfo.diskBound1[z].x+segInfo.diskBound2[z].x*7)/8,k2=k+y*sizex+x; x<=segInfo.spineBound2[z].x; x++, k2++)
			{
				if(imgA[k2]>segPara.boneThresh && binArray[k2]==0)
				{
					if(binArray[k2+sizexy]==1 && maskA[k2+sizexy]!=maskStruct.mask_corticalBone
						&& maskA[k2+sizexy]!=maskStruct.mask_vertebralDisk
						&& maskA[k2+sizexy]!=maskStruct.mask_spongyBone) binArray[k2] = 1;
				}
			}
		}
		for(y=segInfo.spineBound1[z].y; y<=segInfo.diskBound2[z].y; y++)
		{
			for(x=segInfo.diskBound2[z].x,k2=k+y*sizex+x; x<=segInfo.spineBound2[z].x; x++,k2++)
			{
				if(imgA[k2]>segPara.boneThresh && binArray[k2]==0)
				{
					if(binArray[k2+sizexy]==1 && maskA[k2+sizexy]!=maskStruct.mask_corticalBone
						&& maskA[k2+sizexy]!=maskStruct.mask_vertebralDisk
						&& maskA[k2+sizexy]!=maskStruct.mask_spongyBone) binArray[k2] = 1;
				}
			}
		}
	}	// for z



	// 3D blob analysis
	intDynArray blobRank, blobStatus;
	intVec3DynArray blobCentroid;
	NIH_Algo_Blob_Labelling_3D(binImg, blobImg, blobRank, blobCentroid, true);

	blobStatus.SetSize(blobRank.GetSize());

	int ribBlobSizeThresh=900;

	// mark small blobs
	for(k=0; k<blobStatus.GetSize(); k++) 
	{
		if(blobRank[k]<ribBlobSizeThresh) blobStatus[k]=0;
		else blobStatus[k]=1;
	}
	

	// mark all blobs that in the vicinity of spine as potential ribs
	int b;
	for(z=z1; z<z2; z++)
	{
		k=z*sizexy;
		for(y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++)
			{
				b = blobArray[k]-1;
				if(b>=0 && blobStatus[b]==1)
				{
					if((((x>=segInfo.spineBound1[z].x-3 
						|| (z>0 && x>=segInfo.spineBound1[z-1].x-3 && segInfo.spineBound1[z-1].x>0)
						|| (x>=segInfo.spineBound1[z+1].x-3 && segInfo.spineBound1[z+1].x>0))
						&& x<(segInfo.diskBound1[z].x*3+segInfo.diskBound2[z].x)/4)
						||
						(x>(segInfo.diskBound1[z].x+segInfo.diskBound2[z].x*3)/4
						&& (x<=segInfo.spineBound2[z].x+3 || 
						(z>0 && x<=segInfo.spineBound2[z-1].x+3)
						|| x<=segInfo.spineBound2[z+1].x+3)))
						&& y>=(segInfo.diskBound1[z].y+segInfo.diskBound2[z].y*3)/4 
						&& y<=segInfo.spineBound2[z].y) blobStatus[b]=2;
				}
			}
		}
	}	// for z

	// wipe out non-ribs blobs 
	for(z=z1; z<z2; z++)
	{
		k=z*sizexy;
		for(y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++)
			{
				b = blobArray[k]-1;
				if(b>=0 && blobStatus[b]==2)
				{
					binArray[k] = 1;
					maskA[k] = maskStruct.mask_rib;
				}
				else binArray[k] = 0;
			}
		}
	}

	// count all the rib candidates
	int ribCandidateNum=0;		// number of rib candidates
	intDynArray ribBlobCorr;
	for(b=0; b<blobStatus.GetSize(); b++)
	{
		if(blobStatus[b]==2)
		{
			blobStatus[b] = ribCandidateNum;
			ribBlobCorr.Add(b);
			ribCandidateNum++;
		} else blobStatus[b]=-1;
	}

	// match z slice with vertebra num;
	intDynArray zVert;
	int p;
	zVert.SetSize(sizez);
	for(z=0; z<sizez; z++) zVert[z]=-1;
	for(z=0; z<pedicleValley[0]; z++) zVert[z]=0;	// the first vert
	// all other verts
	for(p=0; p<pedicleValley.GetSize()-1; p++)
	{
		for(z=pedicleValley[p]; z<pedicleValley[p+1]; z++) zVert[z]=p+1;	
	}

	// pick the largest component of rib candidates within the vicinity of each vertebra as the real rib
	intDynArray vertToRib_left, vertToRib_right;		// vertebra to rib correspondense
	intDynArray ribToVert;		// rib to vertebra correspondense
	int *corr_mat;				// correspondence matrix
	int vertNum = pedicleValley.GetSize();	// number of vertebra
	corr_mat = new int[ribCandidateNum*vertNum];
	for(k=0; k<ribCandidateNum*vertNum; k++) corr_mat[k]=0;
	// fill the correspondence matrix
	int cv, cb;
	for(z=z1+1; z<z2; z++)
	{
		k=z*sizexy;
		for(y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++)
			{
				b = blobArray[k]-1;
				if(b>=0 && blobStatus[b]>=0)	
				{

					if((((x>=segInfo.spineBound1[z].x-1 || x>=segInfo.spineBound1[z-1].x-1 && segInfo.spineBound1[z-1].x>0)
						&& x<(segInfo.diskBound1[z].x*3+segInfo.diskBound2[z].x)/4)
						||
						(x>(segInfo.diskBound1[z].x+segInfo.diskBound2[z].x*3)/4
						&& (x<=segInfo.spineBound2[z].x+1 || x<=segInfo.spineBound2[z-1].x+1)))&&
						y>=(segInfo.diskBound1[z].y+segInfo.diskBound2[z].y*3)/4 
						&& (y<=segInfo.spineBound2[z].y||y<=segInfo.sprocessBound1[z].y))					
///						if(x>=segInfo.spineBound1[z].x && x<=segInfo.spineBound2[z].x &&
///						y>=(segInfo.diskBound1[z].y+segInfo.diskBound2[z].y*3)/4 && y<=segInfo.spineBound2[z].y) 
					{
						cb = blobStatus[b];	// curblob
						cv = zVert[z];		// curVert
						if(cb>=0 && cv>=0) corr_mat[cb*vertNum+cv]++;
					}
				}
			}
		}
	}	// z

	vertToRib_left.SetSize(vertNum);
	for(k=0; k<vertToRib_left.GetSize(); k++) vertToRib_left[k]=-1;
	vertToRib_right.SetSize(vertNum);
	for(k=0; k<vertToRib_right.GetSize(); k++) vertToRib_right[k]=-1;
	ribToVert.SetSize(ribCandidateNum);
	for(k=0; k<ribToVert.GetSize(); k++) ribToVert[k]=-1;

	// get the rib candidate to vertebra correpondence
	for(cb=0; cb<ribToVert.GetSize(); cb++)
	{
		int mv=-1, mp=0;
		for(cv=0; cv<vertNum; cv++)
		{
			if(corr_mat[cb*vertNum+cv]>mp)
			{
				mp = corr_mat[cb*vertNum+cv];
				mv = cv;
			}
		}
		ribToVert[cb] = mv;
	}

	// get the vertebra to rib correspondence
	for(cb=0; cb<ribCandidateNum; cb++)
	{
		cv = ribToVert[cb];
		if(cv<0) continue;
		b = ribBlobCorr[cb];
		if(blobRank[b]<ribBlobSizeThresh) continue;	// can't be too small

		// centroid need to be outside of spineBound, with certain margin
		if((blobCentroid[b].x>segInfo.spineBound1[blobCentroid[b].z].x+1 ||
			(blobCentroid[b].z!=0 && blobCentroid[b].x>segInfo.spineBound1[blobCentroid[b].z-1].x+1) ||
			blobCentroid[b].x>segInfo.spineBound1[blobCentroid[b].z+1].x+1)
			&&
		   (blobCentroid[b].x<segInfo.spineBound2[blobCentroid[b].z].x-1 ||
			(blobCentroid[b].z!=0 && blobCentroid[b].x<segInfo.spineBound2[blobCentroid[b].z-1].x-1) ||
			blobCentroid[b].x<segInfo.spineBound2[blobCentroid[b].z+1].x-1)) continue;

		if(blobCentroid[b].x<segInfo.spineBound1[blobCentroid[b].z].x)
		{
			if(vertToRib_left[cv]<0 || blobRank[b]>blobRank[vertToRib_left[cv]])
			{
				vertToRib_left[cv] = b;
			}
		} else if(blobCentroid[b].x>segInfo.spineBound2[blobCentroid[b].z].x)
		{
			if(vertToRib_right[cv]<0 || blobRank[b]>blobRank[vertToRib_right[cv]])
			{
				vertToRib_right[cv] = b;
			}
		} else if(blobCentroid[b].x<segInfo.cordCenter[blobCentroid[b].z].x)
		{
			if(vertToRib_left[cv]<0 || blobRank[b]>blobRank[vertToRib_left[cv]])
			{
				vertToRib_left[cv] = b;
			}
		} else {
			if(vertToRib_right[cv]<0 || blobRank[b]>blobRank[vertToRib_right[cv]])
			{
				vertToRib_right[cv] = b;
			}
		}
	}	// for cb

	for(b=0; b<blobStatus.GetSize(); b++) if(blobStatus[b]>=0) blobStatus[b]=0;

	segInfo.t12_vertebra = -1;

	int last_cv=-1;
	for(cv=0; cv<vertNum; cv++)
	{
		if(vertToRib_left[cv]>=0 && (last_cv==-1 || cv-last_cv<2 || (cv-last_cv<=2 &&vertToRib_right[cv]>=0))) 
		{
			segInfo.t12_vertebra = cv;
			blobStatus[vertToRib_left[cv]]=2;
			last_cv = cv;
		}
		if(vertToRib_right[cv]>=0 && (last_cv==-1 || cv-last_cv<2 || (cv-last_cv<=2 &&vertToRib_left[cv]>=0))) 
		{
			segInfo.t12_vertebra = cv;
			blobStatus[vertToRib_right[cv]]=2;
			last_cv = cv;
		}
	}

	// check if there are extra T vert (>6)
	if(pedicleValley.GetSize()-segInfo.t12_vertebra>6) segInfo.t12_vertebra++;

	// convert non-rib to vertebra structure
	for(z=z1; z<z2; z++)
	{
		k=z*sizexy;
		for(y=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++)
			{
				b = blobArray[k]-1;
				if(b>=0 && blobStatus[b]>=0)
				{
					if(blobStatus[b]==0)
					{
						binArray[k] = 0;
						if(x>=segInfo.spineBound1[z].x && x<=segInfo.spineBound2[z].x &&
							y>=(segInfo.diskBound1[z].y+segInfo.diskBound2[z].y)/2 && y<=segInfo.spineBound2[z].y) 
								maskA[k] = maskStruct.mask_corticalBone;
						else 	maskA[k] = maskStruct.mask_otherBone;
					} else maskA[k] = maskStruct.mask_rib;
				}
			}
		}
	}	// for z

	// refine the rib segmentation

	//  fill the holes
	// a little awkward
	for(z=z1; z<=z2; z++)
	{
		int count=0;
		k = z*sizexy;
		for(y=0, k2=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++, k2++)
			{
				if(maskA[k]==maskStruct.mask_rib)
				{
					binArray2D[k2] = 1;
					if(imgA[k+sizexy]>segPara.boneThresh) maskA[k+sizexy]=maskStruct.mask_rib;	// fill next slice, why not previous slice? since we do it the filling slice by slice
					count++;
				}
				else binArray2D[k2] = 0;
			}
		}

		if(count==0) continue;

		CIS_IPA_Dilate(binImg2D, 1, true);
		k = z*sizexy;
		for(y=0, k2=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++, k2++)
			{
				if(binArray2D[k2]==1 && maskA[k]!=maskStruct.mask_rib && imgA[k]>segPara.boneThresh)
				{
					maskA[k]=maskStruct.mask_rib;
				}
			}
		}
		CIS_IPA_Erode(binImg2D, 1, true);
		k = z*sizexy;
		for(y=0, k2=0; y<sizey; y++)
		{
			for(x=0; x<sizex; x++, k++, k2++)
			{
				if(binArray2D[k2]==1 && maskA[k]!=maskStruct.mask_rib)
				{
					maskA[k]=maskStruct.mask_rib;		// will this ever reached? will, no intensity filter
				}
			}
		}
	}	// for z


	// refinement on each 2D slice to recover the missed part of vertebra body
	for(z=segInfo.bound1.z; z<segInfo.bound2.z; z++)
	{
		int count;
		k = z*sizexy;
		// clean the buffer
		for(k2=0; k2<sizexy; k2++) binArray2D[k2]=0;

		// step 1. add back corticle bone
		// for every pixel > boneThresh and within diskBound, add it back
		count = 0;
		for(y=segInfo.spineBound1[z].y; y<=segInfo.spineBound2[z].y; y++)
		{
			for(x=segInfo.spineBound1[z].x-10,k2=y*sizex+x; x<=segInfo.spineBound2[z].x+10; x++, k2++)
			{
				if(maskA[k+k2]==maskStruct.mask_otherBone && imgA[k+k2]>segPara.boneThresh+100) maskA[k+k2]=maskStruct.mask_corticalBone;
			}
		}	// for y

		// step 2: do a connected component analysis
		// every spine pixel inside spineBound should be connected expect for those inside sProcessBound
		for(y=segInfo.spineBound1[z].y; y<=segInfo.spineBound2[z].y; y++)
		{
			for(x=segInfo.spineBound1[z].x-10,k2=y*sizex+x; x<=segInfo.spineBound2[z].x+10; x++, k2++)
			{
				if(maskA[k+k2]==maskStruct.mask_spongyBone || maskA[k+k2]==maskStruct.mask_corticalBone 
					|| maskA[k+k2]==maskStruct.mask_spinalCord || maskA[k+k2]==maskStruct.mask_vertebralDisk) 
				{
					binArray2D[k2]=1;
					count++;
				}
			}
		}	// for y
		if(count==0) continue;
		NIH_Algo_Blob_Labelling(binImg2D, blobImg2D, blobRank2D, blobCentroid2D, false);

		int largest_b, largest_rank;
		int b;
		largest_rank = 0;
		for(b=0; b<blobRank2D.GetSize(); b++)
		{
			if(blobRank2D[b]>largest_rank)
			{
				largest_rank = blobRank2D[b];
				largest_b = b+1;
			}
		}

		for(y=segInfo.spineBound1[z].y; y<=segInfo.spineBound2[z].y; y++)
		{
			for(x=segInfo.spineBound1[z].x-10, k2=y*sizex+x; x<=segInfo.spineBound2[z].x+10; x++, k2++)
			{
				if(blobArray2D[k2]>0 && blobArray2D[k2]!=largest_b) 
				{
					if(x<segInfo.sprocessBound1[z].x || x>segInfo.sprocessBound2[z].x ||
						y<segInfo.sprocessBound1[z].y || y>segInfo.sprocessBound2[z].y)
						maskA[k+k2]=maskStruct.mask_otherBone;
				}
			}
		}	// for y

		// step 3: do a close operation to close hole and small gaps
		for(y=segInfo.spineBound1[z].y; y<=segInfo.spineBound2[z].y; y++)
		{
			for(x=segInfo.spineBound1[z].x-10, k2=y*sizex+x; x<=segInfo.spineBound2[z].x+10; x++, k2++)
			{
				if(maskA[k+k2]==maskStruct.mask_otherBone) binArray2D[k2]=0;
			}
		}	// for y
		CIS_IPA_Dilate(binImg2D, segPara.closeCorticalHoleIteration, false);
		CIS_IPA_Erode(binImg2D, segPara.closeCorticalHoleIteration, false);
		for(y=segInfo.spineBound1[z].y; y<=segInfo.spineBound2[z].y; y++)
		{
			for(x=segInfo.spineBound1[z].x-10, k2=y*sizex+x; x<=segInfo.spineBound2[z].x+10; x++, k2++)
			{
				if(binArray2D[k2]==1 && 
					maskA[k+k2]!=maskStruct.mask_corticalBone
					&& maskA[k+k2]!=maskStruct.mask_spinalCord && maskA[k+k2]!=maskStruct.mask_vertebralDisk 
					&& maskA[k+k2]!=maskStruct.mask_rib) maskA[k+k2]=maskStruct.mask_spongyBone;
			}
		}	// for y

	}	// for z

	// compute the rib cage size
	if(segInfo.t12_vertebra>=0)
	{
		segInfo.ribCageSize = Vec3(0,0,0);
		int maxx, maxy, x0, x1, y0, y1;
		maxx = maxy= 0;
		for(z=segInfo.bound1.x; z<=pedicleValley[segInfo.t12_vertebra]; z++)
		{
			x0=sizex; x1=0; y0=sizey; y1=0;
			for(y=1, k=z*sizexy+sizex; y<sizey-1; y++)
			{
				for(x=0; x<sizex; x++, k++)
				{
					if(imgA[k]>segPara.boneThresh)
					{
						if(x<x0 || x>x1 || y<y0 ||y>y1)
						{
							// check all neighbors to be more robust
							if(imgA[k-1]>segPara.boneThresh && imgA[k+1]>segPara.boneThresh &&
								imgA[k-sizex]>segPara.boneThresh && imgA[k+sizex]>segPara.boneThresh)
							{
								if(x<x0) x0=x;
								if(x>x1) x1=x;
								if(y<y0) y0=y;
								if(y>y1) y1=y;
							}
						}
					}
				}
			}	// for y

			if(x1-x0>maxx) maxx=x1-x0;
			if(y1-y0>maxy) maxy=y1-y0;
		}

		// determine the veterbra length of T11, T12, L1, L2
		int lt11,lt12,ll1,ll2;
		int pSize = pedicleValley.GetSize();
		if(segInfo.t12_vertebra<pSize-1) ll1=pedicleValley[segInfo.t12_vertebra+1]-pedicleValley[segInfo.t12_vertebra];
		else ll1=0;
		if(segInfo.t12_vertebra<pSize-2) ll2=pedicleValley[segInfo.t12_vertebra+2]-pedicleValley[segInfo.t12_vertebra+1];
		else ll2=ll1;
		if(segInfo.t12_vertebra>0) lt12=pedicleValley[segInfo.t12_vertebra]-pedicleValley[segInfo.t12_vertebra-1];
		else lt12=0;
		if(segInfo.t12_vertebra>1) lt11=pedicleValley[segInfo.t12_vertebra-1]-pedicleValley[segInfo.t12_vertebra-2];
		else lt11=lt12;
		if(ll1==0) ll1= ll2 = lt12;
		if(lt12==0) lt12= lt11 = ll1;

		segInfo.ribCageSize.x = (float)maxx*px;
		segInfo.ribCageSize.y = (float)maxy*py;
		segInfo.ribCageSize.z = (float)(lt11+lt12+ll1+ll2)*pz;
	}	// if segInfo

	delete binImg;
	delete blobImg;
	delete binImg2D;
	delete blobImg2D;
	delete corr_mat;

	return CIS_OK;
}