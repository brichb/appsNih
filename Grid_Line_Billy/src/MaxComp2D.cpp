

#include "connected.h"
#include <functional>
#include <stdlib.h>
#include <stdio.h>

#include "NIH_ITK_Utility.h"
//#include <NIH_OpenGL_QT_interface.h>
 #include <CIS_Model_Package.h>
//#include <Dicom_Reader.h>
#include <CIS_Image_Processing_Algo_3D.h>
#include <CIS_Array_Image3D.h>

struct histo2D
{
int intensity;
int count;
};


//void MaxComp2D(CIS_Array_Image3D_short *Left_Overs,int maxslice, int lothresh, int hithresh)
int MaxComp2D(CIS_Array_Image3D_short *Left_Overs, CIS_Array_Image3D_short *Max_Component, int maxslice, int x_low, int x_high, int y_low, int y_high, int Survey)
{

	CIS_Array_Image2D_short *img2Din, *img2Dout;
	img2Din = new CIS_Array_Image2D_short();
	img2Dout = new CIS_Array_Image2D_short();

	img2Din = Left_Overs->GetSliceImage(maxslice,2,false);
	img2Dout = Left_Overs->GetSliceImage(maxslice,2,false);

	int width = img2Din->Num_Cols();
	int height = img2Din->Num_Rows();

	//create array

	int *img2Din_array; int *img2Dout_array;
	img2Din_array = new int[width * height];
	img2Dout_array = new int[width * height];

	int counter = 0;
	for (int y = 0; y<img2Din->Num_Rows(); y++){
		for (int x = 0; x<img2Din->Num_Cols(); x++){
			img2Din_array[counter] = img2Din->FastGet(x,y);
			counter++;
		}
	}
	
	// do connected comp - 8 connectivity
	ConnectedComponents cc(30);
    cc.connected(img2Din_array, img2Dout_array, width, height,
		 std::equal_to<int>(),
		 constant<bool,true>());

	// write output to original image

	counter = 0;
	for (int y = 0; y<img2Din->Num_Rows(); y++){
		for (int x = 0; x<img2Din->Num_Cols(); x++){
			Max_Component->FastSet(x,y,maxslice,img2Dout_array[counter]);
			Left_Overs->FastSet(x,y,maxslice,img2Dout_array[counter]);
			counter++;
		}
	}
// RUN ABOVE HERE
	counter = 0;
	for (int y = 0; y<img2Din->Num_Rows(); y++){
		for (int x = 0; x<img2Din->Num_Cols(); x++){
			img2Dout->FastSet(x,y,img2Dout_array[counter]);
			counter++;
		}
	}
	

	// Fit histogram to blobs and reassign labels based on size

	int data_size_xy = img2Dout->Num_Cols()*img2Dout->Num_Rows();

	int value;
	std::vector<int> intensity;
	intensity.reserve(data_size_xy);
	
	for (int y = 0; y < img2Dout->Num_Rows(); y++){
		for (int x = 0; x < img2Dout->Num_Cols(); x++){
			value = (int)img2Dout->FastGet(x,y);
			intensity.push_back(value);
		}
	}

	int max_intensity = std::numeric_limits<int>::min();
	int min_intensity = std::numeric_limits<int>::max();

	for(int i=0; i<intensity.size(); i++)
	{
		if(intensity[i] > max_intensity)
			max_intensity = intensity[i];

		if(intensity[i] < min_intensity)
			min_intensity = intensity[i];
	}	
	
	int	numbins = max_intensity - min_intensity+1;
	std::vector<int> count;
	count.reserve(numbins);

	for(int i=0; i <numbins; i++)
		count.push_back(0);

	if(numbins == 1)
	{
		int poo = 0;
	}

	int bin_width = (max_intensity - min_intensity)/(numbins-1);
	
	int bin = 0;
	
		for(int i=0; i<data_size_xy; i++)
		{
			bin = floor((intensity[i] - min_intensity)/((float) bin_width)+0.5);
			if(bin < 0)
				bin = 0;
			count[bin] ++;
		}

	std::vector<histo2D> blobs;
	histo2D blob_temp;

	for(int i = 0; i<numbins; i++)
	{
		blob_temp.intensity = min_intensity + i*bin_width;
		blob_temp.count = count[i];
		blobs.push_back(blob_temp);
	}

	// JW- The next few lines I added in to take only the biggest blob from the image.
	int background = 0;
	int big = blobs[0].count;
	for(int n = 0; n < numbins; n++)
	{
		if(blobs[n].count > big)
		{
			big = blobs[n].count;
			background = n;
		}
	}
	
	blobs[background].count = 0; // Erases count of the background

	int max_component = 0;
	 big = blobs[0].count;
	for(int n = 0; n < numbins; n++)
	{
		if(blobs[n].count > big)
		{
			big = blobs[n].count;
			max_component = n;
		}
	}



 //RESETS IMAGES, TURNED OFF WHEN SEEKING ONLY BEST SET OF COLLECTOR  PLANES
	if(Survey == 1)
	{
		for (int y = 0; y<= height; y++)
			for (int x = 0; x<= width; x++)
				{
					if(Left_Overs->FastGet(x,y,maxslice) == blobs[max_component].intensity)
					{
						Left_Overs->FastSet(x,y,maxslice,0); // Will erase the max_component from each slice.
						Max_Component->FastSet(x,y,maxslice,1);
					}
					else if(Left_Overs->FastGet(x,y,maxslice) == blobs[background].intensity)
					{
						Left_Overs->FastSet(x,y,maxslice,0);
						Max_Component->FastSet(x,y,maxslice,0);
					}
					else
					{
						Max_Component->FastSet(x,y,maxslice,0); // Will erase the max_component from each slice.
						Left_Overs->FastSet(x,y,maxslice,1);
					}
				}

	 
		blobs[max_component].count = 0; // Erases the main component count.
	 
	}






	

	int poop = 0;
	for(int i = 0; i < numbins; i++)
	{
		 big = blobs[i].count;
		 if(big > 100) //JW - 200 is arbitrary, you just want to skip over small error "spikes" in the image.
		 {
			poop++;
		 }
	}
	return poop;
//	std::sort(blobs.begin(),blobs.end(), blobcount);
/*
	// relabel blobs (remove first element - and keep next 3 biggest elements)

	// create duplicate image
	

		for (int y = 0; y<Left_Overs->Num_Rows(); y++){
			for (int x = 0; x<Left_Overs->Num_Cols(); x++){
					if (Left_Overs->FastGet(x,y,maxslice) == blobs[1].intensity)
						Left_Overs->FastSet(x,y,maxslice,blobs.size()+4);
					else if (Left_Overs->FastGet(x,y,maxslice) == blobs[2].intensity)
						Left_Overs->FastSet(x,y,maxslice,blobs.size()+3);
					else if (Left_Overs->FastGet(x,y,maxslice) == blobs[3].intensity)
						Left_Overs->FastSet(x,y,maxslice,blobs.size()+2);
					else if (Left_Overs->FastGet(x,y,maxslice) == blobs[4].intensity)
						Left_Overs->FastSet(x,y,maxslice,blobs.size()+1);
					else if (Left_Overs->FastGet(x,y,maxslice) < blobs.size())
						Left_Overs->FastSet(x,y,maxslice,0);
			}
		}


	// Write_Analyze_File("C:\\cvs\\CIPS\\projects\\FinishedRegistrations\\3726745\\3726745_short_mean_0000_slicemaxslice_conncomp_sorted.hdr", *Left_Overs);
	

	std::cout<<"Stop"<<std::endl;
*/
}
