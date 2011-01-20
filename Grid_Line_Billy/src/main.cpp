#include "Grid_Line.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
extern void Label_Maker(std::string * , int);

void main(int argc, char **argv)
{
	srand(time(NULL));
	CIS_Array_Image3D_short *img3Din_1, *img3Din_2, *img3DBIG_copy_1, *img3DBIG_copy_2, *img3Din_turn, *img3Dcopy, *img3Dcopy_2, *img3Dwhat, *img3Dwhat_2, *img3Dwhat_3, *img3Dwhat_4,*img3Dskin, *img3Dskin_2, * img3DSpecs, * img3DSpecs_2, * img3DSpecs_3, *img3DAll, *img3DAll_2;

    string STRING;
	ifstream infile;
	infile.open ("SegmentLibrary.txt");
    std::vector<string> files;  
	while(!infile.eof()) // To get you all the lines.
        {
	        getline(infile,STRING); // Saves the line in STRING.
	        //cout<<STRING; // Prints our STRING.
				
			files.push_back(STRING);
		}
	
	
	cout<<files[2]; //Prints out VECTOR
	infile.close();
	//system ("pause");
	{
	srand(time(NULL));

	// Input image name to propigate in output images.
	std::string Run_Name = argv[1];
	int Run_Name_Length = Run_Name.length();
	Label_Maker(&Run_Name, Run_Name_Length);

	std::string Run_Name_2 = argv[2];
	int Run_Name_Length_2 = Run_Name_2.length();
	Label_Maker(&Run_Name_2, Run_Name_Length_2);
	// Import the input image.
	img3Din_1 = new CIS_Array_Image3D_short();
	img3Din_1->Load(argv[1]);

	img3Din_2 = new CIS_Array_Image3D_short();
	img3Din_2->Load(argv[2]);

	

	//Input: 
	//1.  The two images to compare
	//2.  A user defined DOOO THSIISS

	



	double slice_gap = atof(argv[3]);
	double line_gap_1 = atof(argv[4]);
	double line_gap_2 = atof(argv[5]);
	double User_Point_x = atof(argv[6]);
	double User_Point_y = atof(argv[7]);
	double User_Point_z = atof(argv[8]);
	double User_or_PCA = atof(argv[9]);

	std::string name_1 = argv[3];
	std::string name_2 = argv[4];
	std::string name_3 = argv[5];
	std::string name_4 = argv[6];
	std::string name_5 = argv[7];
	std::string name_6 = argv[8];
	std::string name_7 = argv[9];


	Run_Name.append("_");
	Run_Name.append(name_1);
	Run_Name.append("_");
	Run_Name.append(name_2);
	Run_Name.append("_");
	Run_Name.append(name_3);
	Run_Name.append("_");
	Run_Name.append(name_4);
	Run_Name.append("_");
	Run_Name.append(name_5);
	Run_Name.append("_");
	Run_Name.append(name_6);
	Run_Name.append("_");
	Run_Name.append(name_7);

	Run_Name_2.append("_");
	Run_Name_2.append(name_1);
	Run_Name_2.append("_");
	Run_Name_2.append(name_2);
	Run_Name_2.append("_");
	Run_Name_2.append(name_3);
	Run_Name_2.append("_");
	Run_Name_2.append(name_4);
	Run_Name_2.append("_");
	Run_Name_2.append(name_5);
	Run_Name_2.append("_");
	Run_Name_2.append(name_6);
	Run_Name_2.append("_");
	Run_Name_2.append(name_7);

 
 
	// Define variables denoting maximum number of x,y, and z slices in the image.
	double xmax=img3Din_1->Num_Cols(), ymax=img3Din_1->Num_Rows(), zmax=img3Din_1->Num_Levels();


	double xmax_1=img3Din_1->Num_Cols(), ymax_1=img3Din_1->Num_Rows(), zmax_1=img3Din_1->Num_Levels();
	double xmax_2=img3Din_2->Num_Cols(), ymax_2=img3Din_2->Num_Rows(), zmax_2=img3Din_2->Num_Levels();

	// The next variables represent pixel sizes of the input image.
	double xpixel_1 = img3Din_1->Get_Pixel_SizeX(), ypixel_1 = img3Din_1->Get_Pixel_SizeY(), zpixel_1 = img3Din_1->Get_Pixel_SizeZ();
	double xpixel_2 = img3Din_2->Get_Pixel_SizeX(), ypixel_2 = img3Din_2->Get_Pixel_SizeY(), zpixel_2 = img3Din_2->Get_Pixel_SizeZ();

	int z_add = 0;
	 if(zmax_1 < 200)
	{
		z_add = (int)( (double)(200 - zmax_1)/2.0);
		zmax_1 = zmax_1 + 2*z_add;
	}

	xmax = xmax_1;
	ymax = ymax_1; 
	zmax = zmax_1;


	// Reset the image to be (0,1) instead of (0,1000)
	for(int z = 0; z < zmax; z++)
		for(int y = 0; y < ymax; y++)
			for(int x = 0; x < xmax; x++)
			{
				if(img3Din_1->FastGet(x,y,z) != 0)
				{
					img3Din_1->FastSet(x,y,z,1);
				}
			}


	img3DBIG_copy_1 = new CIS_Array_Image3D_short(xmax, ymax, zmax);
	img3DBIG_copy_1->Set_Pixel_SizeX(img3Din_1->Get_Pixel_SizeX());
	img3DBIG_copy_1->Set_Pixel_SizeY(img3Din_1->Get_Pixel_SizeY());
	img3DBIG_copy_1->Set_Pixel_SizeZ(img3Din_1->Get_Pixel_SizeZ());

	for (int x = 0 ; x < xmax  ; x++)
		for (int y = 0; y < ymax ; y++)
			for (int z = 0 + z_add; z < zmax - z_add  ; z++)
			{
				img3DBIG_copy_1->FastSet(x,y,z,img3Din_1->FastGet(x ,y ,z - z_add ));
 			}

 

	img3Dskin = new CIS_Array_Image3D_short(xmax, ymax, zmax);
	img3Dskin->Set_Pixel_SizeX(img3Din_1->Get_Pixel_SizeX());
	img3Dskin->Set_Pixel_SizeY(img3Din_1->Get_Pixel_SizeY());
	img3Dskin->Set_Pixel_SizeZ(img3Din_1->Get_Pixel_SizeZ());
	for (int x = 0; x < xmax ; x++)
		for (int y = 0; y < ymax; y++)
			for (int z = 0 + z_add; z < zmax - z_add ; z++)
			{
				img3Dskin->FastSet(x,y,z,0);
			}


	img3DSpecs = new CIS_Array_Image3D_short(xmax, ymax, zmax);
	img3DSpecs->Set_Pixel_SizeX(img3Din_1->Get_Pixel_SizeX());
	img3DSpecs->Set_Pixel_SizeY(img3Din_1->Get_Pixel_SizeY());
	img3DSpecs->Set_Pixel_SizeZ(img3Din_1->Get_Pixel_SizeZ());
	for (int x = 0; x < xmax ; x++)
		for (int y = 0; y < ymax; y++)
			for (int z = 0 + z_add; z < zmax - z_add; z++)
			{
				img3DSpecs->FastSet(x,y,z,0);
			}
	// Generates two images to be used in the Skin_it/Randomizer functions.
 	img3Dwhat = new CIS_Array_Image3D_short(xmax, ymax, zmax);
	img3Dwhat->Set_Pixel_SizeX(img3Din_1->Get_Pixel_SizeX());
	img3Dwhat->Set_Pixel_SizeY(img3Din_1->Get_Pixel_SizeY());
	img3Dwhat->Set_Pixel_SizeZ(img3Din_1->Get_Pixel_SizeZ());
	for (int x = 0; x < xmax ; x++)
		for (int y = 0; y < ymax; y++)
			for (int z = 0 + z_add; z < zmax - z_add ; z++)
			{
				img3Dwhat->FastSet(x,y,z,0);
			}

	img3Dwhat_2 = new CIS_Array_Image3D_short(xmax, ymax, zmax);
	img3Dwhat_2->Set_Pixel_SizeX(img3Din_1->Get_Pixel_SizeX());
	img3Dwhat_2->Set_Pixel_SizeY(img3Din_1->Get_Pixel_SizeY());
	img3Dwhat_2->Set_Pixel_SizeZ(img3Din_1->Get_Pixel_SizeZ());
	for (int x = 0; x < xmax ; x++)
		for (int y = 0; y < ymax; y++)
			for (int z = 0 + z_add; z < zmax - z_add ; z++)
			{
				img3Dwhat_2->FastSet(x,y,z,0);
			}

	img3DAll = new CIS_Array_Image3D_short(xmax, ymax, zmax);
	img3DAll->Set_Pixel_SizeX(img3Din_1->Get_Pixel_SizeX());
	img3DAll->Set_Pixel_SizeY(img3Din_1->Get_Pixel_SizeY());
	img3DAll->Set_Pixel_SizeZ(img3Din_1->Get_Pixel_SizeZ());
	for (int x = 0; x < xmax ; x++)
		for (int y = 0; y < ymax; y++)
			for (int z = 0 + z_add; z < zmax - z_add ; z++)
			{
				img3DAll->FastSet(x,y,z,0);
			}

	 z_add = 0;
	 if(zmax_2 < 200)
	{
		z_add = (int)( (double)(200 - zmax_2)/2.0);
		zmax_2 = zmax_2 + 2*z_add;
	}

 	xmax = xmax_2;
	ymax = ymax_2; 
	zmax = zmax_2;
	// Reset the image to be (0,1) instead of (0,1000)
	for(int z = 0; z < zmax; z++)
		for(int y = 0; y < ymax; y++)
			for(int x = 0; x < xmax; x++)
			{
				if(img3Din_2->FastGet(x,y,z) != 0)
				{
					img3Din_2->FastSet(x,y,z,1);
				}
			}

 

	img3DBIG_copy_2 = new CIS_Array_Image3D_short(xmax, ymax, zmax);
	img3DBIG_copy_2->Set_Pixel_SizeX(img3Din_2->Get_Pixel_SizeX());
	img3DBIG_copy_2->Set_Pixel_SizeY(img3Din_2->Get_Pixel_SizeY());
	img3DBIG_copy_2->Set_Pixel_SizeZ(img3Din_2->Get_Pixel_SizeZ());

	for (int x = 0; x < xmax  ; x++)
		for (int y = 0; y < ymax ; y++)
			for (int z = 0 + z_add; z < zmax - z_add  ; z++)
			{
				img3DBIG_copy_2->FastSet(x,y,z,img3Din_2->FastGet(x ,y ,z - z_add ));
 			}

	img3Dcopy_2 = new CIS_Array_Image3D_short(zmax,ymax, xmax);
	img3Dcopy_2->Set_Pixel_SizeX(img3Din_2->Get_Pixel_SizeX());
	img3Dcopy_2->Set_Pixel_SizeY(img3Din_2->Get_Pixel_SizeY());
	img3Dcopy_2->Set_Pixel_SizeZ(img3Din_2->Get_Pixel_SizeZ());
	for (int x = 0; x < xmax  ; x++)
		for (int y = 0; y < ymax ; y++)
			for (int z = 0 + z_add; z < zmax - z_add ; z++)
			{
				img3Dcopy_2->FastSet(x,y,z,img3Din_2->FastGet(x ,y ,z - z_add));
			}


	img3Dskin_2 = new CIS_Array_Image3D_short(xmax, ymax, zmax);
	img3Dskin_2->Set_Pixel_SizeX(img3Din_2->Get_Pixel_SizeX());
	img3Dskin_2->Set_Pixel_SizeY(img3Din_2->Get_Pixel_SizeY());
	img3Dskin_2->Set_Pixel_SizeZ(img3Din_2->Get_Pixel_SizeZ());
	for (int x = 0; x < xmax ; x++)
		for (int y = 0; y < ymax; y++)
			for (int z = 0 + z_add; z < zmax - z_add ; z++)
			{
				img3Dskin_2->FastSet(x,y,z,0);
			}


	img3DAll_2 = new CIS_Array_Image3D_short(xmax, ymax, zmax);
	img3DAll_2->Set_Pixel_SizeX(img3Din_2->Get_Pixel_SizeX());
	img3DAll_2->Set_Pixel_SizeY(img3Din_2->Get_Pixel_SizeY());
	img3DAll_2->Set_Pixel_SizeZ(img3Din_2->Get_Pixel_SizeZ());
	for (int x = 0; x < xmax ; x++)
		for (int y = 0; y < ymax; y++)
			for (int z = 0 + z_add; z < zmax - z_add; z++)
			{
				img3DAll_2->FastSet(x,y,z,0);
			}


	img3Dwhat_3 = new CIS_Array_Image3D_short(xmax, ymax, zmax);
	img3Dwhat_3->Set_Pixel_SizeX(img3Din_2->Get_Pixel_SizeX());
	img3Dwhat_3->Set_Pixel_SizeY(img3Din_2->Get_Pixel_SizeY());
	img3Dwhat_3->Set_Pixel_SizeZ(img3Din_2->Get_Pixel_SizeZ());
	for (int x = 0; x < xmax ; x++)
		for (int y = 0; y < ymax; y++)
			for (int z = 0 + z_add; z < zmax - z_add; z++)
			{
				img3Dwhat_3->FastSet(x,y,z,0);
			}

	img3Dwhat_4 = new CIS_Array_Image3D_short(xmax, ymax, zmax);
	img3Dwhat_4->Set_Pixel_SizeX(img3Din_2->Get_Pixel_SizeX());
	img3Dwhat_4->Set_Pixel_SizeY(img3Din_2->Get_Pixel_SizeY());
	img3Dwhat_4->Set_Pixel_SizeZ(img3Din_2->Get_Pixel_SizeZ());
	for (int x = 0; x < xmax ; x++)
		for (int y = 0; y < ymax; y++)
			for (int z = 0 + z_add; z < zmax - z_add ; z++)
			{
				img3Dwhat_4->FastSet(x,y,z,0);
			}

	img3DSpecs_2 = new CIS_Array_Image3D_short(xmax, ymax, zmax);
	img3DSpecs_2->Set_Pixel_SizeX(img3Din_2->Get_Pixel_SizeX());
	img3DSpecs_2->Set_Pixel_SizeY(img3Din_2->Get_Pixel_SizeY());
	img3DSpecs_2->Set_Pixel_SizeZ(img3Din_2->Get_Pixel_SizeZ());
	for (int x = 0; x < xmax ; x++)
		for (int y = 0; y < ymax; y++)
			for (int z = 0 + z_add; z < zmax - z_add ; z++)
			{
				img3DSpecs_2->FastSet(x,y,z,0);
			}


	Grid_Line * test;
	test = new Grid_Line(img3Din_1, img3Din_2, img3DBIG_copy_1, img3DBIG_copy_2, img3Dskin, img3Dskin_2, img3Dwhat, img3Dwhat_2,img3Dwhat_3,img3Dwhat_4, img3DAll, img3DAll_2, img3DSpecs, img3DSpecs_2, Run_Name, Run_Name_2,  slice_gap, line_gap_1 , line_gap_2 ,User_Point_x, User_Point_y, User_Point_z ,User_or_PCA);



	int poopy = 0;

	std::cout << "Skinning Object 1 ..." << "\n";
	test->Skin_it(1, 0, 1, img3DBIG_copy_1, img3Dskin);
 	Write_Analyze_File((char *)(new std::string(Run_Name + (*(new std::string("_SKIN.hdr")))))->c_str(), *img3Dskin);
	std::cout << "Skinning Object 2 ..." << "\n";
	test->Skin_it(1, 0, 2, img3DBIG_copy_2, img3Dskin_2);
 	Write_Analyze_File((char *)(new std::string(Run_Name_2 + (*(new std::string("_SKIN.hdr")))))->c_str(), *img3Dskin_2);

	std::cout << "Beginning search over Connection Planes for Object 1..." << "\n";
	test->Plane_Maker(img3Dskin, img3DBIG_copy_1, img3Dwhat, img3Dwhat_2, img3Dwhat_3, img3DSpecs_2, img3DSpecs, img3DAll, 1 );
	std::cout << "Beginning search over Connection Planes for Object 2 and Parameterization proccess..." << "\n";
	test->Plane_Maker(img3Dskin_2, img3DBIG_copy_2, img3Dwhat_3, img3Dwhat_4, img3Dskin, img3DBIG_copy_1 , img3DSpecs_2, img3DAll_2, 2 );


}

}
