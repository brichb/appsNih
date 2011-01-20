// Grid_Line made by Jeremy Watt jermwatt@gmail.com
// Last updated on 5/25/10.

#include <CIS_Array_Image3D.h>
#include <vector>
#include <time.h>
#include <fstream>



#ifndef ___GRID_LINE__
#define ___GRID_LINE__
class Grid_Line {
public:
	Grid_Line( CIS_Array_Image3D_short * ,CIS_Array_Image3D_short * ,  CIS_Array_Image3D_short * , CIS_Array_Image3D_short * , CIS_Array_Image3D_short * , CIS_Array_Image3D_short * , CIS_Array_Image3D_short * , CIS_Array_Image3D_short * , CIS_Array_Image3D_short * , CIS_Array_Image3D_short * , CIS_Array_Image3D_short * , CIS_Array_Image3D_short * , CIS_Array_Image3D_short *, CIS_Array_Image3D_short * , std::string Run_Name, std::string Run_Name_2,  double slice_gap, double line_gap_1 , double line_gap_2 ,double Input_Point_x, double Input_Point_y, double Input_Point_z ,double User_or_PCA ) ; 
 
	int Feature(int, int, int, CIS_Array_Image3D_short *);
	void Skin_it(int skin_1, int skin_2,int max,  CIS_Array_Image3D_short*  , CIS_Array_Image3D_short* );
	
    void Plane_Maker( CIS_Array_Image3D_short * ,CIS_Array_Image3D_short *, CIS_Array_Image3D_short * , CIS_Array_Image3D_short * , CIS_Array_Image3D_short * , CIS_Array_Image3D_short * , CIS_Array_Image3D_short * , CIS_Array_Image3D_short *, int Ford); 

	void Comparison(int,  std::vector<int>, std::vector<std::vector<int>>, std::vector<std::vector<std::vector<int>>>, std::vector<std::vector<std::vector<std::vector<int>>>>, CIS_Array_Image3D_short * , CIS_Array_Image3D_short * , CIS_Array_Image3D_short *, CIS_Array_Image3D_short *, CIS_Array_Image3D_short *);
	void centroid(int z_slice, int dim, int type, CIS_Array_Image3D_short *imgIN);

	CIS_Array_Image3D_short *img3Din_1, *img3Din_2, *img3DBIG_copy_1, *img3DBIG_copy_2, *img3Din_turn, *img3Dcopy, *img3Dcopy_2, *img3Dwhat, *img3Dwhat_2, *img3Dwhat_3, *img3Dwhat_4,*img3Dskin, *img3Dskin_2, * img3DSpecs, * img3DSpecs_2, * img3DSpecs_3, *img3DAll, *img3DAll_2;

private:
	int xmax, ymax, zmax;
	int xcenter, ycenter, zcenter, Final_Radius;
	int  zmax_original_1, zmax_original_2;
 	int xmax_1, ymax_1, zmax_1;
	int xmax_2, ymax_2, zmax_2;
	int Final_1, Final_2, Final_3, Final_1_1, Final_2_1, Final_3_1, Final_1_2, Final_2_2, Final_3_2;
	int mop;
	int Final_Radius_1, Final_Radius_2;

	double xpixel_1, ypixel_1, zpixel_1;
	double xpixel_2, ypixel_2, zpixel_2;
	double PI, xpixel, ypixel, zpixel;
    double slice_gap, line_gap_1 , line_gap_2 ,User_Point_x, User_Point_y, User_Point_z ,User_or_PCA, edge_close_3, cut_off_portion;
	double theta_rad_rotate, phi_rad_rotate , theta_rad, phi_rad;

	std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>>> MIN_MAX;
	std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> ALL_POINTS, ALL_POINTS_1, ALL_POINTS_2;
	std::vector<std::vector<std::vector<std::vector<int>>>> VERT, Faces_1, Faces_2;
	std::vector<std::vector<std::vector<int>>> HORIZ, Normals_1, Normals_2;
	std::vector<std::vector<int>>  Dual, Highs_Lows_1, Highs_Lows_2;
	std::vector<int> Point, Projections_1, Projections_2;

	std::vector<std::vector<std::vector<std::vector<double>>>> ALL_MEAS, ALL_MEAS_1, ALL_MEAS_2;
 	std::vector<std::vector<std::vector<double>>> DOUBLE_VERT;
	std::vector<std::vector<double>> DOUBLE_HORIZ;
	std::vector<double> DOUBLE_Dual,  ALL_ROTATED_VECTORS_1, ALL_ROTATED_VECTORS_2, Maxs_1, Maxs_2;

	std::string Run_Name, Run_Name_2;

};

//#include "Grid_Line.cpp"
#endif ___GRID_LINE__