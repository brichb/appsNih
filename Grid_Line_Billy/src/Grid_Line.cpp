#include "Grid_Line.h"

extern int MaxComp2D(CIS_Array_Image3D_short *, CIS_Array_Image3D_short *, int slice, int xl, int xh, int yl, int yh, int s);
extern void Label_Maker(std::string * , int);
extern std::vector<double> Principle_Components_Finder(CIS_Array_Image3D_short * img3Dconv);

Grid_Line::Grid_Line(CIS_Array_Image3D_short * image_in, CIS_Array_Image3D_short * image_in_2, CIS_Array_Image3D_short * image_copy, CIS_Array_Image3D_short * image_copy_2,  CIS_Array_Image3D_short * image_skin, CIS_Array_Image3D_short * image_skin_2, CIS_Array_Image3D_short * image_temp,CIS_Array_Image3D_short * image_temp_2, CIS_Array_Image3D_short * image_temp_3, CIS_Array_Image3D_short * image_temp_4, CIS_Array_Image3D_short * image_all, CIS_Array_Image3D_short * image_all_2,  CIS_Array_Image3D_short * image_spec,  CIS_Array_Image3D_short * image_spec_2,   std::string Run, std::string Run_2,  double slice_mes, double line_1 , double line_2, double ball_divide, double point_divide, double clip_1 ,double clip_2)
{

	////// The Images ///////
	img3Din_1 = image_in;
	img3Din_2 = image_in_2;
	img3DBIG_copy_1 = image_copy;
	img3DBIG_copy_2 = image_copy_2;
	img3Dwhat = image_temp;
	img3Dwhat_2 = image_temp_2;
	img3Dwhat_3 = image_temp_3;
	img3Dwhat_4 = image_temp_4;
	img3DAll = image_all;
	img3DAll_2 = image_all_2;
	img3Dskin = image_skin;
	img3Dskin_2 = image_skin_2;
	img3DSpecs = image_spec;
 	img3DSpecs_2 = image_spec_2;

	////// The ints initialized here /////// 

	mop = 1;
	theta_rad_rotate = phi_rad_rotate = theta_rad = phi_rad = 0;
	Final_Radius = 0, Final_Radius_1 = 0, Final_Radius_2 = 0;
	xcenter = ycenter = zcenter = 0;
	xmax = ymax = zmax = 0;
	xmax_1=img3DBIG_copy_1->Num_Cols(), ymax_1=img3DBIG_copy_1->Num_Rows(), zmax_1=img3DBIG_copy_1->Num_Levels();
	xmax_2=img3DBIG_copy_2->Num_Cols(), ymax_2=img3DBIG_copy_2->Num_Rows(), zmax_2=img3DBIG_copy_2->Num_Levels();

	Final_1 = Final_2 = Final_3 = 0;
	Final_1_1 = Final_2_1 = Final_3_1 = Final_1_2 = Final_2_2 = Final_3_2;
	// The next variables represent pixel sizes of the input image.
	xpixel_1 = img3Din_1->Get_Pixel_SizeX(), ypixel_1 = img3Din_1->Get_Pixel_SizeY(), zpixel_1 = img3Din_1->Get_Pixel_SizeZ();
	xpixel_2 = img3Din_2->Get_Pixel_SizeX(), ypixel_2 = img3Din_1->Get_Pixel_SizeY(), zpixel_2 = img3Din_2->Get_Pixel_SizeZ();

	////// The doubles initialized here /////// 
	PI = 3.141592653589793238462643383279502884197;

	slice_gap = slice_mes;
	line_gap_1 = line_1 ; 
	line_gap_2 = line_2;
	User_Point_x = ball_divide;
	User_Point_y = point_divide;
	User_Point_z = clip_1;
	User_or_PCA = clip_2;

 	////// The vectors initialized here /////// 
	for(int i = 0; i < 3; i++)
	{
		Point.push_back(0);
	}
	
	for(int i = 0; i < 2; i++)
	{
		Dual.push_back(Point);
		DOUBLE_Dual.push_back(0);
	}
	for(int i = 0; i < 30; i++)
	{
		HORIZ.push_back(Dual);
		DOUBLE_HORIZ.push_back(DOUBLE_Dual);
	}
	for(int i = 0; i < 30; i++)
	{
		VERT.push_back(HORIZ);
		DOUBLE_VERT.push_back(DOUBLE_HORIZ);
	}

	for(int i = 0; i <= 1; i++)
	{
		ALL_POINTS.push_back(VERT);
		ALL_POINTS_1.push_back(VERT);
		ALL_POINTS_2.push_back(VERT);
		ALL_MEAS.push_back(DOUBLE_VERT);
		ALL_MEAS_1.push_back(DOUBLE_VERT);
		ALL_MEAS_2.push_back(DOUBLE_VERT);

	}

	for(int i = 0; i <= 1; i++)
	{
		MIN_MAX.push_back(ALL_POINTS);
	}


	////// The strings initialized here /////// 
	Run_Name_2 = Run_2;
	Run_Name = Run ;

}

void Grid_Line::centroid(int z_slice, int dim, int type, CIS_Array_Image3D_short *imgIN)
{
	// Centroid/Center of Mass is the pixel (xcenter,ycenter,zcenter) where xcenter is the average of all x-coordinates
	// of pixels in the object, and likewise with ycenter and zcenter.
	int i = 0, j = 0, k = 0, n = 0;
	 xcenter = 0, ycenter = 0, zcenter = 0;
	if(dim == 3)
	{
		for(int z = 0; z < zmax; z++)
		{
			for(int y = 0; y < ymax; y++)
			{
				for(int x = 0; x < xmax; x++)
				{
					if(imgIN->FastGet(x,y,z) != 0) 
					{
						i = i + x;
						j = j + y;
						k = k + z;
						n = n + 1;
					}
				}
			}
		}
		if(n > 0)
		{
			// floor( input + 0.5) rounds a double to the nearest integer.
 			xcenter = (int)floor((double)(i/n) + 0.5);
			ycenter = (int)floor((double)(j/n) + 0.5);
			zcenter = (int)floor((double)(k/n) + 0.5);
 
		}
			// If the centroid of your object is (0,0,0) you know you've loaded in a blank image.
			std::cout << "The centroid of this object is (x,y,z) = (" << xcenter << "," << ycenter << "," << zcenter << ")" << "\n";
	}

 
}



int Grid_Line::Feature(int x, int y, int z, CIS_Array_Image3D_short * Skinned_Image)
{
	// "Feature" was made a seperate function so that one can easily 'plug in' any surface-feature one wishes into
	// the Program.  It is currently set to capture the volume of intersection between a paralellapiped (e.g. a 3D
	// Rectangle) and the interior of the Object.  You can easily exchange the paralellapiped with an Ellipsoid by
	// adjusting the limits of search and uncommenting the lines below.

	int feature_counter = 0;

	for(int i = x - Final_1; i <= x + Final_1; i++)
		for(int j = y - Final_2; j <= y + Final_2; j++)
			for(int k = z - Final_3; k <= z + Final_3; k++)
			{
				//double dist = (i - x)*(i - x)  + (j - y)*(j - y) + (double)(zmax)/(double)(xmax)*(k - z)*(k - z);
				
				//if( Skinned_Image->FastGet(i,j,k) != 0 && dist < 50  )
				//{
					feature_counter++;
				//}
			}

	return feature_counter;

}


void Grid_Line::Plane_Maker( CIS_Array_Image3D_short * Skinned_Image, CIS_Array_Image3D_short * UnSkinned_Image, CIS_Array_Image3D_short * Ray_Image, CIS_Array_Image3D_short * Major_Axis_Plane_Projection, CIS_Array_Image3D_short * Second_Skinned_Image, CIS_Array_Image3D_short * Second_UnSkinned_Image, CIS_Array_Image3D_short * Intersection_Points, CIS_Array_Image3D_short * Intersection_Planes, int Ford)
{
	////////////////////////////////////////////////////////////////////////////////////////////////////
	// FUNCTION DESCRIPTION:
	// 
	/* 

	SECTION 1:


	i)	 Using the Input Points PCAs (Principle Components) of the input object in Skinned_Image, or User defined points, 
		 we rotate a set of symmetric points from one hemisphere of verticies of a dodecahedron to be centered about the
		 Main Point (PCA or user defined point).  These points are all loaded into ALL_TEST_VECTORS and then, those 
		 that need to be rotated, are converted below.
	
	ii)  For each of these points we then find Planes that: 
		 a)  Have this point as a normal vector and
		 b)  Stem from the object's center of mass (xcenter, ycenter, zcenter)

 	iii) Because we use MaxComp2D.cpp which takes as input a z-Plane, and because we wish to visualize the results 
		(e.g.  points collected in each Mother Plane) we make each Plane thick initially so that an injective map 
		 to the zcenter-plane  may be calculated.  For each input vector n its corresponding Plane and injective
		 map to the zcenter plane is found as follows

		a)  Each Plane is foound by searching over the whole image for pixels (x,y,z) satisfying
		|n1(x - xcenter) + n2(y - ycenter) + n3(z - zcenter)| < epsilon where n = (n1, n2, n3) and epsilon is 
		chosen large enough to leave the resulting Plane, which we call the 'Mother Plane' for n, a few pixels 'thick'.
		b) Using the ratios of the largest x/y/z to smallest x/y/z pixel in the thick plane we determine which primary
		plane (a x-plane, y-plane, or z-plane) maps injectively to the Mother Plane.  
		c)  We 'thin' the thick Mother Plane and find the injective map by projecting the proper primary plane
		pixel-by-pixel onto our Mother Plane, keeping only the first pixel of the Plane intersected by each pixel
		projection from the proper primary plane.  Below is a simple 2-D illustration:

		a) The thick Mother Plane is generated

							\\\                    
							 \\\
							  \\\
							   \\\
								\\\
								 \\\
						 Thick Mother Plane

		b) The proper plane is found for the injective map to the zcenter plane.

								|			\\\                    
								|			 \\\
								|			  \\\
								|			   \\\
								|               \\\
								|                \\\
							Primary plane     Thick Mother Plane


		
		c) The proper primary plane is projected onto the thick MotherPlane. Keeping only the first pixel of intersection 
		of the projection, and storing the coordinates of said pixel, gives the map and the thin  Mother Plane.
		 
								  Projection
								|   ----->	\\\                    
								| 	----->	 \\\
								|	----->	  \\\
								|	----->	   \\\
								|   ----->      \\\
								|   ----->       \\\
							Primary plane    Thick Mother Plane

							  <map>			\                   
							  <map>			 \
							  <map>			  \
							  <map>		       \
							  <map>             \
							  <map>              \
						   Injective Map   Thinned Mother Plane
							(a matrix)

		The Map can then be treated as one mapping to the zcenter-plane (by searching it properly).

							 <map> <map> <map> <map> <map> <map>


	SECTION 2:
	In this section we
		i)   Using the Inputed Points we construct a Normal passing through the Object's centroid 
			 (xcenter,ycenter,zcenter) for each Mother Plane found in Section 1. 
		ii)  Each Mother Plane is then 'slid' along this normal in both directions to find a tight top and bottom 
			 translated Mother Plane, we call a 'Daughter Plane', which contains the Object.  In other words:
			 we search parallel versions of each Mother Plane by adding/subtracting components of the Mother Plane's
			 Normal in order to find two parallel Planes 'above' and 'below' the Mother Plane which tightly contain
			 our Object.  We find the distance (via the Mother Plane Normal) between this 'top' and 'bottom'
			 (call this d).

								   ____
								  /    \
								 /   .  \
								/________\
						The Object with '.' centroid

							   Mother Plane
								   __|_
								  /  | \
							<--------|--------> (Normal)
								/____|___\
									 |
				The Object with Plane through centroid, and Plane Normal

 							   |   ____   |
							   |  /    \  |
							<--|----------|--> 
							   |/________\|
							   |  ---d--- |
				The Object with tight 'top' (left) and 'bottom' (right) Daughter Planes


		iii) We then interect the Object between the 'top' and 'bottom' Planes a user-defined number m of times by 
			 translating the 'bottom' Plane d/m along the Plane Normal.  We call each parallel Plane intersected with
			 the Object an Intersecting Plane.

 							   |   |___|  |
							   |  /|   |\ |
							   |-->|-->|->|  
							   |/__|___|_\|
							   |   |   |  |
				The Object intersected with user-defined number of paraellel Planes

		iv)  Each Intersecting Plane is translated to the zcenter-plane, which then fed into MaxComp2D, returning 
		the number of 2D Connected Components it contains (called the Connect # of the particular Intersecting Plane).

			|
			|
			| ------------>    ==========   ------------->  MaxComp2D( ========== ) ----------> Connect #
			|
			|
		Intersecting Plane -> Translated to zcenter-plane ->   fed into MaxComp2D -> returns num of 2D Connected
																					 Components in Intersectig Plane

		v)  The average Connect # is calculated for each set of Intersecting Planes. 
		
		
	SECTION 3:
	In this section we
		i)   Find the minimum sum average Connect #'s among the two sets of corresponding Intersecting Planes for our
			 two input Objects.  We call these the Winning Connection Planes for each Object.  If several sets of Intersection
			 Planes have the minimum sum of Connect #'s we pick the pair most similiar to the Mother Plane whose normal is (0,1,0), 
			 since from visual inspection Connection Planes parallel to this Mother Plane tend to cut the gut-organs well e.g.
			 the slices have a ver low Connect #, and  each 2-D curve is well behaved.
		ii)  For each Winning Plane we find the bounding square containing all (hopefully one) components where the 
			 sides of said square are ortogonal to the Normal Vector for the Plane examined.


 							   |   |___|  |
							   |  /|   |\ |
							   |-->|-->|->|  
							   |/__|___|_\|
							   |   |   |  |
				The Object intersected with user-defined number of Winning Planes

										   ______________________
					|				      |					     |
					|				      |				 ___	 |
					|				      |		 __		/   \    |
					|				      |	    /  \___/ 	|	 |
					|  -------------->    |	   /		   /     |
					|				      |	  |		      |	     |
					|				      |	   \		 /	     |
					|			          |		\_______/	     |
					|				      |	    				 |
					|				      |______________________|
				 Intersection         Intersection Plane Rotated parallel 
				   Plane          to this window,with 2D Component of Object



		iii)   The tightest square containing all (hopefully only one) Object components with sides orthogonal
				to the Plane Normal is found.

					 ______________________
					|					   |				  ________________
					|				 ___   |			  	 |           ___  |
					|		 __		/   \  |				 |	 __		/   \ |
					|	    /  \___/ 	|  |				 |  /  \___/ 	| |
					|	   /		   /   |	--------->	 | /		   /  |
					|	  |		      |	   |				 | |		  |   |	  
					|	   \		 /	   |				 | \		 /	  |
					|		\_______/	   |				 |  \_______/	  |
					|______________________|				  ----------------
				  Tightest square for Intersection Plane Object components if found

		iv)    Sides of the bounding rectangle are partitioned a user-defined number of times, and points at 
			   these partitions on the rectangle are projected onto the Object component(s), which we refer to
			   as Intersection Lines. The first point of intersection with the Object component(s) from each 
			   projection is stored, we call these points Intersection Points. 

		
							|		|
						 ___|_______|____
						|   |       |__  |                              o
						|	|_	   /   \ |                      o 
					____|__/  \___/    |_|___                  o           o
						| /		      /  |   -------------->
					____|_|		     |___|___	              o          o
						| \		    /	 |                             o
						|  \_______/|	 |                       o
						 ---|-------|-----
							|       |
						Intersection Lines				Points on the Object component	
	

		 v)    We calculate the Objects' Surface Feature at each Intersection Point.  The current feature, we call
		       the Manay Feature, is the volume of intersection of a sphereoid/parallelepiped (or any object really),
			   called the Seed here, intersected with the interior of Object.  
			   a)  We define the axes of Seed by d/m, the distance via Plane Normal between two neighboring Intersector
			       Planes.
			   b)  We normalize the volume of the Seed by dividing by the Seed's full volume.

		 
		 vi)   We compare the Manay Feature of each point on corresponding Connection Planes, and map the difference
			   back onto each Object's surface.

	      v)   Having found all of the Intersection Points and calculated the Feature differences at each of these points, we now 
			   calculate the global difference between the two surfaces and map local differences simultaneously back to both Objects' surfaces.

			   The local differences are mapped simultaneously back onto the two Objects' surfaces.  Differences are mapped back 
			   from largest to smallest, with larger differences taking precedence on the Objects' surfaces, and by color-rank.

	*/


	

	// Here we set the correct maximum slice numbers depending on which Object is being proccessed.
	if(Ford == 1)
	{
		xmax = xmax_1;
		ymax = ymax_1;
		zmax = zmax_1;
	}
	else // Ford == 2
	{
		xmax = xmax_2;
		ymax = ymax_2;
		zmax = zmax_2;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////

	////// BEGIN SECTION 1 ///////////////////////////////////////////////////////////////////////////
 
	// We'll have to first center the image in a larger space, of course.
	// We'll port these into the PC so they are not calculated again.
 
	// This finds the centroid of the Object.
	std::cout << "Finding the centroid of Object " << Ford << "." << "\n";
	centroid(0,3,0,UnSkinned_Image);

	std::vector<std::vector<double>> Rotation, Rotation_2;
	std::vector<double> ALL_TEST_VECTORS;
	std::vector<double> ALL_Vectors, Line, ALL_ROTATED_VECTORS;

	if(User_or_PCA == 1)
	{

		std::cout << "Finding the PCAs of Object " << Ford << "." << "\n";
		std::vector<double> Principles = Principle_Components_Finder(UnSkinned_Image);

		for(int i = 0; i < Principles.size(); i++)
		{
			ALL_TEST_VECTORS.push_back(Principles[i]);
		}

		std::cout << "Your Main Point is the first PCA: (" << Principles[0] << "," << Principles[1] << "," << Principles[2] << ")" << "\n";


		ALL_TEST_VECTORS.push_back(0 );
		ALL_TEST_VECTORS.push_back(0.5257311121191336 );
		ALL_TEST_VECTORS.push_back(0.85065080835204);
		
		ALL_TEST_VECTORS.push_back(0 );
		ALL_TEST_VECTORS.push_back(-0.5257311121191336 );
		ALL_TEST_VECTORS.push_back(0.85065080835204);

		ALL_TEST_VECTORS.push_back(0.85065080835204 );
		ALL_TEST_VECTORS.push_back(0 );
		ALL_TEST_VECTORS.push_back(0.5257311121191336);

	}
	else
	{

		std::cout << "Your Main Point was inputed by the User: (" << User_Point_x << "," << User_Point_y << "," << User_Point_z << ")" << "\n";
		ALL_TEST_VECTORS.push_back(User_Point_x);
		ALL_TEST_VECTORS.push_back(User_Point_y);
		ALL_TEST_VECTORS.push_back(User_Point_z);

		ALL_TEST_VECTORS.push_back(0 );
		ALL_TEST_VECTORS.push_back(0.5257311121191336 );
		ALL_TEST_VECTORS.push_back(0.85065080835204);
		
		ALL_TEST_VECTORS.push_back(0 );
		ALL_TEST_VECTORS.push_back(-0.5257311121191336 );
		ALL_TEST_VECTORS.push_back(0.85065080835204);

		ALL_TEST_VECTORS.push_back(0.85065080835204 );
		ALL_TEST_VECTORS.push_back(0 );
		ALL_TEST_VECTORS.push_back(0.5257311121191336);

		ALL_TEST_VECTORS.push_back(-0.85065080835204 );
		ALL_TEST_VECTORS.push_back(0 );
		ALL_TEST_VECTORS.push_back(0.5257311121191336);

		ALL_TEST_VECTORS.push_back(0 );
		ALL_TEST_VECTORS.push_back(0.9341723589627157 );
		ALL_TEST_VECTORS.push_back(0.3568220897730899);
 

	}
	/*

	i)	 Using the Inputed Points PCAs (Principle Components) of the input object in Skinned_Image, or User defined points, 
		 we rotate a set of symmetric points from one hemisphere of verticies of a dodecahedron to be centered about the
		 Main Point (PCA or user defined point).  These points are all loaded into ALL_TEST_VECTORS and then, those 
		 that need to be rotated, are converted below.

		 Note: the Main Point is the FIRST point inputed into ALL_TEST_VECTORS below.
	 */


	// Note: As of 7/1/10 due to memory issues I had only been able to load 6 - 9 points at a run.

	// These are other points from the upper hemisphere of the 32-vertice dodecahedron.
	/*
	ALL_TEST_VECTORS.push_back(0 );
	ALL_TEST_VECTORS.push_back(0.5257311121191336 );
	ALL_TEST_VECTORS.push_back(0.85065080835204);
	
	ALL_TEST_VECTORS.push_back(0 );
	ALL_TEST_VECTORS.push_back(-0.5257311121191336 );
	ALL_TEST_VECTORS.push_back(0.85065080835204);

	ALL_TEST_VECTORS.push_back(0.85065080835204 );
	ALL_TEST_VECTORS.push_back(0 );
	ALL_TEST_VECTORS.push_back(0.5257311121191336);

	ALL_TEST_VECTORS.push_back(-0.85065080835204 );
	ALL_TEST_VECTORS.push_back(0 );
	ALL_TEST_VECTORS.push_back(0.5257311121191336);

	ALL_TEST_VECTORS.push_back(0 );
	ALL_TEST_VECTORS.push_back(0.9341723589627157 );
	ALL_TEST_VECTORS.push_back(0.3568220897730899);

	ALL_TEST_VECTORS.push_back(0 );
	ALL_TEST_VECTORS.push_back(-0.9341723589627157 );
	ALL_TEST_VECTORS.push_back(0.3568220897730899);

	ALL_TEST_VECTORS.push_back(0.3568220897730899 );
	ALL_TEST_VECTORS.push_back(0 );
	ALL_TEST_VECTORS.push_back(0.9341723589627157);

	ALL_TEST_VECTORS.push_back(-0.3568220897730899 );
	ALL_TEST_VECTORS.push_back(0 );
	ALL_TEST_VECTORS.push_back(0.9341723589627157);
	
	ALL_TEST_VECTORS.push_back(0.5773502691896257 );
	ALL_TEST_VECTORS.push_back(0.5773502691896257 );
	ALL_TEST_VECTORS.push_back(0.5773502691896257);

	ALL_TEST_VECTORS.push_back(-0.5773502691896257 );
	ALL_TEST_VECTORS.push_back(0.5773502691896257 );
	ALL_TEST_VECTORS.push_back(0.5773502691896257);

	ALL_TEST_VECTORS.push_back(0.5773502691896257 );
	ALL_TEST_VECTORS.push_back(-0.5773502691896257 );
	ALL_TEST_VECTORS.push_back(0.5773502691896257);

	ALL_TEST_VECTORS.push_back(-0.5773502691896257 );
	ALL_TEST_VECTORS.push_back(-0.5773502691896257 );
	ALL_TEST_VECTORS.push_back(0.5773502691896257);
	*/
	

	// The lines below rotate all of the Input Points to be centered at the Main Point.
	
		std::cout << "Centering Input Points about the Main Point for Object " << Ford << "." << "\n";

		double x = xcenter, y = ycenter, z = zcenter;

		int TEST_LENGTH =  (int)(((double)ALL_TEST_VECTORS.size())/3.0);
		
		// The following lines compute rotation matrix 'Rotation' which, multiplied by the other
		// Input Points rotates them about a vector orthogonal to the Main Point and the z-axis
		double phi = acos(ALL_TEST_VECTORS[2]);
		double s = sin(phi);
		double c = cos(phi);

		// (cross_x, cross_y, cross_z) = cross-product of the Main Point and the z-axis (0,0,1),
		// e.g. it is a vector orthogonal to the plane spanned by the Main Point and (0,0,1).  We
		// will rotate all the other Input Points about this axis towards the Main Point. (see acom-
		// -panying Manuscript for picture of this)
		double cross_x = ALL_TEST_VECTORS[1];
		double cross_y = - ALL_TEST_VECTORS[0];
		double cross_z = 0;
		
		// 'Rotation' is the matrix of rotation about the cross product towards the Main Point.  For
		// equations for the cross product and rotation matrix search Wikipedia for 'cross-product'
		// and 'Rotation Matrix -> Rotation Matrix given an Axis and Angle'.

		Line.push_back(cross_x*cross_x  + (1 - cross_x*cross_x)*c);
		Line.push_back(cross_x*cross_y*(1 - c) - 0*s);
		Line.push_back(cross_x*0*(1 - c) + cross_y*s);
		Rotation.push_back(Line);
		Line.clear();

		Line.push_back(cross_x*cross_y*(1 - c) + 0*s);
		Line.push_back(cross_y*cross_y + (1 - cross_y*cross_y)*c);
		Line.push_back(cross_y*0*(1 - c) -cross_x*s);
		Rotation.push_back(Line);
		Line.clear();

		Line.push_back(cross_x*0*(1 - c) - cross_y*s);
		Line.push_back(cross_y*0*(1 - c) + cross_x*s);
		Line.push_back(0*0 + (1 - 0*0)*c);
		Rotation.push_back(Line);
		Line.clear();

		// We load the Main Point into ALL_ROTATED_VECTORS_1/2 which will be accessed below when we
		// test Intersecting Planes, and when we finally sample points (we store in two seperate
		// vectors as well because each set, while corresponding, will be different for the two
		// Objects).
		// ALL_ROTATED_VECTORS is a general storage vector to carry each Object's points through 
		// Section 1.
		for(int i = 0; i < 1; i++)
		{
			ALL_ROTATED_VECTORS.push_back(ALL_TEST_VECTORS[3*i]);
			ALL_ROTATED_VECTORS.push_back(ALL_TEST_VECTORS[3*i + 1]);
			ALL_ROTATED_VECTORS.push_back(ALL_TEST_VECTORS[3*i + 2]);

			if(Ford == 1) // The Main Point for Object 1
			{
				ALL_ROTATED_VECTORS_1.push_back(ALL_TEST_VECTORS[3*i]);
				ALL_ROTATED_VECTORS_1.push_back(ALL_TEST_VECTORS[3*i + 1]);
				ALL_ROTATED_VECTORS_1.push_back(ALL_TEST_VECTORS[3*i + 2]);
			}
			if(Ford == 2)// The Main Point for Object 2
			{
				ALL_ROTATED_VECTORS_2.push_back(ALL_TEST_VECTORS[3*i]);
				ALL_ROTATED_VECTORS_2.push_back(ALL_TEST_VECTORS[3*i + 1]);
				ALL_ROTATED_VECTORS_2.push_back(ALL_TEST_VECTORS[3*i + 2]);
			}
		}
		
		// Now we load the remainder of the points, rotated, into the corresponding vectors.
		if(ALL_ROTATED_VECTORS.size() > 3)
		{
			for(int i = 1; i < TEST_LENGTH - 1; i++)
			{	
				double x = ALL_TEST_VECTORS[3*i]*Rotation[0][0] + ALL_TEST_VECTORS[3*i + 1]*Rotation[0][1] + ALL_TEST_VECTORS[3*i + 2]*Rotation[0][2];
				double y = ALL_TEST_VECTORS[3*i]*Rotation[1][0] + ALL_TEST_VECTORS[3*i + 1]*Rotation[1][1] + ALL_TEST_VECTORS[3*i + 2]*Rotation[1][2];
				double z = ALL_TEST_VECTORS[3*i]*Rotation[2][0] + ALL_TEST_VECTORS[3*i + 1]*Rotation[2][1] + ALL_TEST_VECTORS[3*i + 2]*Rotation[2][2];
				
				ALL_ROTATED_VECTORS.push_back(x);
				ALL_ROTATED_VECTORS.push_back(y);
				ALL_ROTATED_VECTORS.push_back(z);

				if(Ford == 1) // The rotated Input Points for Object 1
				{
					ALL_ROTATED_VECTORS_1.push_back(ALL_TEST_VECTORS[3*i]);
					ALL_ROTATED_VECTORS_1.push_back(ALL_TEST_VECTORS[3*i + 1]);
					ALL_ROTATED_VECTORS_1.push_back(ALL_TEST_VECTORS[3*i + 2]);
				}
				if(Ford == 2) // The rotated Input Points for Object 2
				{
					ALL_ROTATED_VECTORS_2.push_back(ALL_TEST_VECTORS[3*i]);
					ALL_ROTATED_VECTORS_2.push_back(ALL_TEST_VECTORS[3*i + 1]);
					ALL_ROTATED_VECTORS_2.push_back(ALL_TEST_VECTORS[3*i + 2]);
				}

			}
		}
		ALL_TEST_VECTORS.clear();
	

 
	// 'Face_Planes' will carry all Mother Planes generated by the Input Points.
	std::vector<std::vector<std::vector<std::vector<int>>>>	Face_Planes; 
	// 'Projection_Ordering' will contain a marker designating which primary plane the injective
	//  mapping is to.
	std::vector<int> Projection_Ordering; // 

	int ALL_LENGTH = (int)( (double)(ALL_ROTATED_VECTORS.size())/3.0);
	 
	for(int PCA_index = 0; PCA_index < ALL_LENGTH; PCA_index++)
	{

		/*

		ii)  For each of the Input Points we then find Planes that: 
			 a)  Have this point as a normal vector and
			 b)  Stem from the object's center of mass (xcenter, ycenter, zcenter)

 		iii) Because we use MaxComp2D.cpp which takes as input a z-Plane, and because we wish to visualize the results 
			(e.g.  points collected in each Mother Plane) we make each Plane thick initially so that an injective map 
			 to the zcenter-plane  may be calculated.  For each input vector n its corresponding Plane and injective
			 map to the zcenter plane is found as follows

			a)  Each Plane is foound by searching over the whole image for pixels (x,y,z) satisfying
			|n1(x - xcenter) + n2(y - ycenter) + n3(z - zcenter)| < epsilon where n = (n1, n2, n3) and epsilon is 
			chosen large enough to leave the resulting Plane, which we call the 'Mother Plane' for n, a few pixels 'thick'.
			b) Using the ratios of the largest x/y/z to smallest x/y/z pixel in the thick plane we determine which primary
			plane (a x-plane, y-plane, or z-plane) maps injectively to the Mother Plane.  
			c)  We 'thin' the thick Mother Plane and find the injective map by projecting the proper primary plane
			pixel-by-pixel onto our Mother Plane, keeping only the first pixel of the Plane intersected by each pixel
			projection from the proper primary plane.  Below is a simple 2-D illustration:

			a) The thick Mother Plane is generated

								\\\                    
								 \\\
								  \\\
								   \\\
									\\\
									 \\\
							 Thick Mother Plane

		*/

		// test_* is the * component of the PCA'th Input Point.
		double test_x = ALL_ROTATED_VECTORS[3*PCA_index];
		double test_y = ALL_ROTATED_VECTORS[3*PCA_index + 1]; 
		double test_z = ALL_ROTATED_VECTORS[3*PCA_index + 2];

		std::cout << "Generating the thick Mother Plane for the " << PCA_index + 1 << "th Input Point for Object " << Ford << "." << "\n";

		// low_*, high_* are the lowest/highest * components of the Plane.
		int low_x = xmax, high_x = 0, low_y = ymax, high_y = 0, low_z = zmax, high_z = 0;
		for(int z = 0; z < zmax; z++)
			for(int y = 0; y < ymax; y++)
				for(int x = 0; x < xmax; x++)
				{
					if( fabs( (double)test_x*(double)(x - xcenter )/(double)xmax  + test_y*(double)(y - ycenter)/(double)ymax + test_z*(double)(z - zcenter)/(double)zmax) < .009)
					{
						if( x < low_x)
						{
							low_x = x;
						}
						if(x > high_x)
						{
							high_x = x;
						}
						if( y < low_y)
						{
							low_y = y;
						}
						if( y > high_y)
						{
							high_y = y;
						}
						if( z < low_z )
						{
							low_z = z;
						}
						if( z > high_z )
						{
							high_z = z;
						}

						Ray_Image->FastSet(x,y,z,1000);
					}
				}

 

	/*

		

	iii)

		b) The proper plane is found for the injective map to the zcenter plane.

								|			\\\                    
								|			 \\\
								|			  \\\
								|			   \\\
								|               \\\
								|                \\\
							Primary plane     Thick Mother Plane

		c) The proper primary plane is projected onto the thick MotherPlane. Keeping only the first pixel of intersection 
		of the projection, and storing the coordinates of said pixel, gives the map and the thin  Mother Plane.
		 
								  Projection
								|   ----->	\\\                    
								| 	----->	 \\\
								|	----->	  \\\
								|	----->	   \\\
								|   ----->      \\\
								|   ----->       \\\
							Primary plane    Thick Mother Plane

							  <map>			\                   
							  <map>			 \
							  <map>			  \
							  <map>		       \
							  <map>             \
							  <map>              \
						   Injective Map   Thinned Mother Plane
							(a matrix)

		The Map can then be treated as one mapping to the zcenter-plane (by searching it properly).

							 <map> <map> <map> <map> <map> <map>
	*/


		std::cout << "Thinning the Mother Plane for the " << PCA_index + 1 << "th Input Point and creating " << "\n"; 
		std::cout << "this Mother Plane's Injective Map to the zcenter-plane for Object " << Ford << "." << "\n";

		// 'Mother_Plane' is a matrix which will contain a single Mother Plane, and is input to Face_Planes.
		std::vector<std::vector<std::vector<int>>> Mother_Plane; //Projection_Map;
		std::vector<std::vector<int>> Mother_Plane_Row;  
		std::vector<int> Mother_Plane_Point;  

		int no_run_twice = 0; // This is to assure we don't project two primary planes onto the Mother Plane.
	
		// What if high_x = xmax and high_y = ymax && high_z = zmax?  We revert to "map_z" case below, put
		// first pruposefully because in most of our cases xmax = ymax > zmax.  If working with different 
		// sized images the cases map_z/x/y will need to be re-ordered (or a general case made).
		//if( high_x == xmax -1 && high_y == ymax -1)
		if(fabs(test_z) >= fabs(test_y) && fabs(test_z) >= fabs(test_x))
		{
			Projection_Ordering.push_back(0);
			no_run_twice++;
			for(int y = 0; y < ymax; y++)
				for(int x = 0; x < xmax; x++)
				{
					for(int z = 0; z < zmax; z++)
					{
							if(Ray_Image->FastGet(x,y,z) != 0)
							{
								if( Major_Axis_Plane_Projection->FastGet(x,y,zcenter) != 0 )
								{
									Ray_Image->FastSet(x,y,z,0);
								}
								else
								{

									Major_Axis_Plane_Projection->FastSet(x,y,zcenter,Ray_Image->FastGet(x,y,z)); //  Projection of x-y top plane onto our Mother Plane.
									
									Mother_Plane_Point.push_back(x);
									Mother_Plane_Point.push_back(y);
									Mother_Plane_Point.push_back(z);

									Mother_Plane_Row.push_back(Mother_Plane_Point);

									Mother_Plane_Point.clear();
								
								}
							}

					}
					if(x == xmax - 1) 
					{
						Mother_Plane.push_back(Mother_Plane_Row);

						Mother_Plane_Row.clear();
					}
				}
		}

	

	
		if(fabs(test_y) >= fabs(test_x) && fabs(test_y) >= fabs(test_z) && no_run_twice == 0)	
		{
			Projection_Ordering.push_back(1);
			no_run_twice++;
			for(int z = 0; z < zmax; z++)
				for(int x = 0; x < xmax; x++)
				{
					for(int y = 0; y < ymax; y++)
					{
						if(Ray_Image->FastGet(x,y,z) != 0)
							{
								if( Major_Axis_Plane_Projection->FastGet(x,ycenter,z) != 0 )
								{
									Ray_Image->FastSet(x,y,z,0);
								}
								else
								{

									Major_Axis_Plane_Projection->FastSet(x,ycenter,z,Ray_Image->FastGet(x,y,z));  //  Projection of xz top plane onto our Mother Plane.
									
									Mother_Plane_Point.push_back(x);
									Mother_Plane_Point.push_back(y);
									Mother_Plane_Point.push_back(z);

									Mother_Plane_Row.push_back(Mother_Plane_Point);

									Mother_Plane_Point.clear();

								}
						}

					}
					if(x == xmax - 1) 

					{
						Mother_Plane.push_back(Mother_Plane_Row);

						Mother_Plane_Row.clear();

					}		
							
						
				}
		}


		if( fabs(test_x) >= fabs(test_y) && fabs(test_x) >= fabs(test_z) &&  no_run_twice == 0)
		{
			Projection_Ordering.push_back(2);
			no_run_twice++;
			for(int z = 0; z < zmax; z++)
				for(int y = 0; y < ymax; y++)
				{
					for(int x = 0; x < xmax; x++)
					{
						if(Ray_Image->FastGet(x,y,z) != 0)
							{
									if( Major_Axis_Plane_Projection->FastGet(xcenter,y,z) != 0 )
									{
										Ray_Image->FastSet(x,y,z,0);
									}
									else
									{

										Major_Axis_Plane_Projection->FastSet(xcenter,y,z,Ray_Image->FastGet(x,y,z));  // Projection of yz top plane onto our Mother Plane.
										
										Mother_Plane_Point.push_back(x);
										Mother_Plane_Point.push_back(y);
										Mother_Plane_Point.push_back(z);

										Mother_Plane_Row.push_back(Mother_Plane_Point);

										Mother_Plane_Point.clear();

									}
							}

					}

						if(y == ymax - 1) 
						{
							Mother_Plane.push_back(Mother_Plane_Row);

							Mother_Plane_Row.clear();

						}

				}
		}
		

		// Because each Mother Plane was thick to begin wtih, and was then thinned, no point of the thin Mother Plane
		// lies on the Object's centroid (which we want).  To correct for this error we simply drag the entire Mother
		// Plane to the Object's centroid by the very middle pixel of the Mother Plane.

		// row/column give us the very middle of the Mother Plane.
		int row = (int)((double)(Mother_Plane.size())/2.0);
		int column = (int)((double)(Mother_Plane[row].size())/2.0);

		int xt = Mother_Plane[row][column][0];
		int yt = Mother_Plane[row][column][1];
		int zt = Mother_Plane[row][column][2];
		
		// This nested for-loop drags the Mother Plane by its center pixel to the Object's centroid.
		for(int i = 0; i < Mother_Plane.size(); i++)
			for(int j = 0; j < Mother_Plane[i].size(); j++)
			{
				Mother_Plane[i][j][0] = Mother_Plane[i][j][0] - xt + xcenter;
				Mother_Plane[i][j][1] = Mother_Plane[i][j][1] - yt + ycenter;
				Mother_Plane[i][j][2] = Mother_Plane[i][j][2] - zt + zcenter;
			} 


		Face_Planes.push_back(Mother_Plane);

		// Reset the 
		Mother_Plane.clear();
 
		// Reset images used
		for(int z = 0; z < zmax; z++)
			for(int y = 0; y < ymax; y++)
				for(int x = 0; x < xmax; x++)
				{	
					Ray_Image->FastSet(x,y,z,0);
					Major_Axis_Plane_Projection->FastSet(x,y,z,0);
				}

	}


	////// END SECTION 1 ///////////////////////////////////////////////////////////////////////////

	///////BEGIN SECTION 2//////////////////////////////////////////////////////////////////////////
 
	 
	std::cout << "Beginning Section 2 for Object " << Ford << "." << "\n";

		
	// Highs_Lows will contain the parameter of the Mother Plane Normal necessary to add/subtract from the each
	// Mother Plane in order to find the two Daughter Planes bounding the object.
	std::vector<std::vector<int>> Highs_Lows;
	std::vector<int> HLs; // HLs will be a single set of parameters for one Mother Plane, loaded into Highs_Lows
	std::vector<double> Max_Comp_Average; // Will contain the average Connect # for each set of Intersector Planes
	std::vector<std::vector<std::vector<int>>> Normals; // Will contain Mother Plane Normals
 
		for(int PCA_index = 0; PCA_index < ALL_LENGTH; PCA_index++)
		{

			std::vector<int> Temp_Point, Max_Comp_Scores; // Max_Comp_Scores carries the Connect # for each Intersector Plane
			std::vector<std::vector<int>> R_Vector; // R_Vector is the temporary vector containing the Mother Plane Normal
			
		/*

		i)   Using the Inputed Points we construct a Normal passing through the Object's centroid 
			 (xcenter,ycenter,zcenter) for each Mother Plane found in Section 1. 

		 */
			std::cout << "Constructing Normal for " << PCA_index + 1 << "th Input Point for Objectt " << Ford << ".\n";

			double test_x = ALL_ROTATED_VECTORS[3*PCA_index];
			double test_y = ALL_ROTATED_VECTORS[3*PCA_index + 1]; 
			double test_z = ALL_ROTATED_VECTORS[3*PCA_index + 2];
			int morp = 1;

			// This while-loop finds the very last pixel of this Mother Plane Normal contained in the image space.
			while((double)morp*(test_x) + (double)xcenter < xmax && (double)morp*(test_x)+ (double)xcenter > 0 && (double)morp*(test_y) + (double)ycenter < ymax && (double)morp*(test_y) + (double)ycenter > 0  && (double)morp*(test_z) + (double)zcenter < zmax && (double)morp*(test_z) + (double)zcenter > 0 )
			{
				int x = (int)floor( (double)morp*test_x + 0.5) + xcenter;
				int y = (int)floor( (double)morp*test_y + 0.5) + ycenter;
				int z = (int)floor( (double)morp*test_z + 0.5) + zcenter;
				morp++;
			}
			morp--;
			int up_morp = morp;

			morp = 0;

			// This while-loop finds the very first pixel of this Mother Plane Normal contained in the image space.
			while((double)morp*(test_x) + (double)xcenter < xmax && (double)morp*(test_x)+ (double)xcenter > 0 && (double)morp*(test_y) + (double)ycenter < ymax && (double)morp*(test_y) + (double)ycenter > 0  && (double)morp*(test_z) + (double)zcenter < zmax && (double)morp*(test_z) + (double)zcenter > 0 )
			{
				int x = (int)floor( (double)morp*test_x + 0.5) + xcenter;
				int y = (int)floor( (double)morp*test_y + 0.5) + ycenter;
				int z = (int)floor( (double)morp*test_z + 0.5) + zcenter;
				morp--;
			}
			morp++;

			int down_morp = morp;
			morp = down_morp;

			// This while-loop generates the Mother Plane Normal from first to last pixel in the image space.

			while(morp < up_morp)
			{
				int x = (int)floor( (double)morp*test_x + 0.5) + xcenter;
				int y = (int)floor( (double)morp*test_y + 0.5) + ycenter;
				int z = (int)floor( (double)morp*test_z + 0.5) + zcenter;

				Temp_Point.push_back(x);
				Temp_Point.push_back(y);
				Temp_Point.push_back(z);

				R_Vector.push_back(Temp_Point);
				Temp_Point.clear();

				morp++;
			}
			
			down_morp = -down_morp; // This is the indicie such that Normal[down_morp] = (xcenter,ycenter,zcenter)

 
		

			/*

			ii)  Each Mother Plane is then 'slid' along this normal in both directions to find a tight top and bottom 
				 translated Mother Plane, we call a 'Daughter Plane', which contains the Object.  In other words:
				 we search parallel versions of each Mother Plane by adding/subtracting components of the Mother Plane's
				 Normal in order to find two parallel Planes 'above' and 'below' the Mother Plane which tightly contain
				 our Object.  We find the distance (via the Mother Plane Normal) between this 'top' and 'bottom'
				 (call this d).

									   ____
									  /    \
									 /   .  \
									/________\
							The Object with '.' centroid

								   Mother Plane
									   __|_
									  /  | \
								<--------|--------> (Normal)
									/____|___\
										 |
					The Object with Plane through centroid, and Plane Normal

 								   |   ____   |
								   |  /    \  |
								<--|----------|--> 
								   |/________\|
								   |  ---d--- |
					The Object with tight 'top' (left) and 'bottom' (right) Daughter Planes

			*/
			std::cout << "Finding Daughter Planes for " << PCA_index + 1 << "th Input Point for Objectt " << Ford << ".\n";

			int high = 0; // The parameter of the Mother Plane Normal corresponding to the 'top' Daughter Plane
			int low = 0; // The parameter of the Mother Plane Normal corresponding to the 'bottom' Daughter Plane
	
			int p = 0; // 'p' parameterizes the Mother Plane Normal, we add points on the Normal to the Mother Plane to generate the Daughter Planes we test
			int sum = 0; // sum is a counter which counts how many pixels of the Object lie in the current Daughter Plane
			int check = 0; // We shift the Mother Plane 'up' until its Daughter Plane contains none of the Object

			//  In the following while loop the TOP Daughter Plane is found.
			while(check == 0 )
			{
				for(int j = 0; j < Face_Planes[PCA_index].size(); j++)
					for(int k = 0; k < Face_Planes[PCA_index][j].size(); k++)
					{
						sum = sum + Skinned_Image->FastGet(Face_Planes[PCA_index][j][k][0] + R_Vector[p + down_morp][0] - xcenter , Face_Planes[PCA_index][j][k][1] + R_Vector[p + down_morp][1] - ycenter, Face_Planes[PCA_index][j][k][2] + R_Vector[p + down_morp][2] - zcenter);
					}

				int poopy = 0;
				if(sum > 0)
				{
					if(p + down_morp == R_Vector.size() - 1)
					{
						check = 1;
						high = p + down_morp;
						sum = 0;
					}
					p++;
					sum = 0;
				}
				else
				{
					check = 1;
					high = p + down_morp;
					sum = 0;
				}
			}

			//  In the following while loop the BOTTOM Daughter Plane is found.
	 
			check = 0; // We shift the Mother Plane 'down' until its Daughter Plane contains none of the Object
			p = 0;
			while(check == 0 )
			{
				for(int j = 0; j < Face_Planes[PCA_index].size(); j++)
					for(int k = 0; k < Face_Planes[PCA_index][j].size(); k++)
					{
						sum = sum + Skinned_Image->FastGet(Face_Planes[PCA_index][j][k][0] + R_Vector[down_morp - p][0] - xcenter, Face_Planes[PCA_index][j][k][1] + R_Vector[down_morp - p][1] - ycenter, Face_Planes[PCA_index][j][k][2] + R_Vector[down_morp - p][2] - zcenter);
					}
				if(sum > 0)
				{
					if(down_morp - p == 0)
					{
						check = 1;
						low = 0;
						sum = 0;
					}
					p++;
					sum = 0;
				}
				else
				{
					check = 1;
					low = down_morp - p;
					sum = 0;
				}
			}

 
		// All the vectors are updated/reset.
		Normals.push_back(R_Vector);
		HLs.push_back(high);
		HLs.push_back(low);
		Highs_Lows.push_back(HLs);
		HLs.clear();




	/*

			iii) We then interect the Object between the 'top' and 'bottom' Planes a user-defined number m of times by 
				 translating the 'bottom' Plane d/m along the Plane Normal.  We call each parallel Plane intersected with
				 the Object an Intersecting Plane.

 								   |   |___|  |
								   |  /|   |\ |
								   |-->|-->|->|  
								   |/__|___|_\|
								   |   |   |  |
					The Object intersected with user-defined number of paraellel Planes

			iv)  Each Intersecting Plane is translated to the zcenter-plane, which then fed into MaxComp2D, returning 
			the number of 2D Connected Components it contains (called the Connect # of the particular Intersecting Plane).

				|
				|
				| ------------>    ==========   ------------->  MaxComp2D( ========== ) ----------> Connect #
				|
				|
			Intersecting Plane -> Translated to zcenter-plane ->   fed into MaxComp2D -> returns num of 2D Connected
																						 Components in Intersectig Plane

	*/

		std::cout << "Construcing Connection Planes finding average Connect # Planes for " << "\n";
		std::cout << PCA_index + 1 << "th Input Point for Objectt " << Ford << ".\n";

		// 'top_bot_partition' is the distance 'd' between the two Daughter Planes
	double top_bot_partition = (double)(high - low)/slice_gap;

		// Remember we want (ideally) a point-to-point comparison across two Objects, but what if any of the
		// two Daughter Planes are too close together this indicates the Object is  an innappropriate candidate 
		// for our parameterization method, thus the program ends if any Daughter Planes are found to be too close.

		//Billy testing more planes
	/*if(top_bot_partition < slice_gap + 1)
		//{
			//std::cout << "The " << PCA_index + 1 << "'th set of Daughter Planes were too close together for your given slice gap of " << slice_gap << ".  You'll need to rerun with a smaller slie gap." << "\n";
			//exit(0);
		}*/


		int partition_num = 0; // 'partition_num' indexes the Connection Planes
		// To avoid Connection Planes near the tips of the Object, where there are not many points, we use 
		// 'extra' and 'extra_2' to scoot us in from the absolute tip of the Object (where the Daughter Planes lie).
		double extra = top_bot_partition/(0.5*slice_gap);
		double extra_2 = extra;
		if(extra < 2)
		{
			extra_2 = 2;
		}

		top_bot_partition = ((double)high- extra_2 - (double)low - extra_2)/slice_gap;

		if(top_bot_partition < slice_gap + 1)
		{
			top_bot_partition = 1;
			std::cout << "The " << PCA_index + 1 << "'th set of Daughter Planes were too close together for your given slice gap of " << slice_gap << ".  You'll need to rerun with a smaller slie gap." << "\n";
			exit(0);
		}

		double JUMP = (double)low + extra_2;
		
		std::cout << "Construcing Connection Planes finding average Connect # Planes for " << "\n";
		std::cout << PCA_index + 1 << "th Input Point for Objectt " << Ford << ".\n";

		while(JUMP < (double)high ) 
		{
			int part = (int)floor(JUMP + 0.5);

			for(int j = 0; j < Face_Planes[PCA_index].size(); j++)
				for(int k = 0; k < Face_Planes[PCA_index][j].size(); k++)
				{
					// 'Major_Axis_Plane_Projection' is the Intersection Plane projected onto the zcenter-slice.  
					// The 'Projection_Ordering' determines the correct shape-preserving map to the zcenter-slice.
					if(Projection_Ordering[PCA_index] == 0 && UnSkinned_Image->FastGet(Face_Planes[PCA_index][j][k][0] + R_Vector[part][0] - xcenter ,Face_Planes[PCA_index][j][k][1] +  R_Vector[part][1] - ycenter ,Face_Planes[PCA_index][j][k][2] +  R_Vector[part][2] - zcenter) != 0)
					{
						Major_Axis_Plane_Projection->FastSet(Face_Planes[PCA_index][j][k][0] + R_Vector[part][0] - xcenter, Face_Planes[PCA_index][j][k][1] + R_Vector[part][1] - ycenter , zcenter, 1000);
					}
					if(Projection_Ordering[PCA_index] == 1 && UnSkinned_Image->FastGet(Face_Planes[PCA_index][j][k][0] + R_Vector[part][0] - xcenter ,Face_Planes[PCA_index][j][k][1] +  R_Vector[part][1] - ycenter ,Face_Planes[PCA_index][j][k][2] +  R_Vector[part][2] - zcenter) != 0)
					{
						Major_Axis_Plane_Projection->FastSet(Face_Planes[PCA_index][j][k][0] + R_Vector[part][0] - xcenter , Face_Planes[PCA_index][j][k][2] + R_Vector[part][2] - zcenter , zcenter, 1000);
					}
					if(Projection_Ordering[PCA_index] == 2 && UnSkinned_Image->FastGet(Face_Planes[PCA_index][j][k][0] + R_Vector[part][0] - xcenter ,Face_Planes[PCA_index][j][k][1] +  R_Vector[part][1] - ycenter ,Face_Planes[PCA_index][j][k][2] +  R_Vector[part][2] - zcenter) != 0)
					{
						Major_Axis_Plane_Projection->FastSet(Face_Planes[PCA_index][j][k][1] + R_Vector[part][1] - ycenter , Face_Planes[PCA_index][j][k][2] + R_Vector[part][2] - zcenter , zcenter, 1000);
					}
				}
 
			
				// Here we compute the Connect # of this Intersection Plane
				int components = MaxComp2D(Major_Axis_Plane_Projection, Major_Axis_Plane_Projection, zcenter, 0, xmax, 0, zmax, 0);
				
				// In case we select a Daughter Plane as one Intersection Plane, which I believe is impossible.
				if( components == 0)
				{
					partition_num--;
				}
				//Write_Analyze_File("Plane_Check.hdr", *Major_Axis_Plane_Projection); 

				// Reset the MaxComp2D image.
				for(int x = 0; x < xmax; x++)
					for(int y = 0;  y <  ymax; y++)
					{
						Major_Axis_Plane_Projection->FastSet(x,y,zcenter,0);
					}
				Max_Comp_Scores.push_back(components);

			partition_num++;
			JUMP = JUMP + top_bot_partition;
		}
		
		// Here we take calculate the average Connect # for this set of Interesction Planes
		double average_score = 0;
		for(int i = 0; i < partition_num; i++)
		{
			average_score = average_score + Max_Comp_Scores[i];
		}
		
		average_score = average_score/(double)(partition_num );

		if(Ford == 1)
		{
			Maxs_1.push_back(average_score);
		}
		else // Ford == 2
		{
			Maxs_2.push_back(average_score);
		}
		Max_Comp_Scores.clear();
		R_Vector.clear();
 

 
	}


	
	if(Ford == 1)
	{
 
		Highs_Lows_1 = Highs_Lows;
		Normals_1 = Normals;
		Projections_1 = Projection_Ordering;

		Highs_Lows.clear();
		Face_Planes.clear();
		Normals.clear();
		Projection_Ordering.clear();

	}

	if(Ford  == 2)
	{
 
		Highs_Lows_2 = Highs_Lows;
		Normals_2 = Normals;
		Projections_2 = Projection_Ordering;

		Highs_Lows.clear();
		Face_Planes.clear();
		Normals.clear();
		Projection_Ordering.clear();


		////// END SECTION 2 ///////////////////////////////////////////////////////////////////////////

		///////BEGIN SECTION 3//////////////////////////////////////////////////////////////////////////
		/*SECTION 3:
		In this section we
			i)   Find the minimum sum of average Connect #'s among the two sets of corresponding Intersecting Planes for our
				 two input Objects.  We call these the Winning Connection Planes for each Object.  If several sets of Intersection
				 Planes have the minimum sum of Connect #'s we pick the pair most similiar to the Mother Plane whose normal is (0,1,0), 
				 since from visual inspection Connection Planes parallel to this Mother Plane tend to cut the gut-organs well e.g.
				 the slices have a ver low Connect #, and  each 2-D curve is well behaved.
		*/



		
		double Winning_Average = Maxs_1[0] + Maxs_2[0];
		int Winning_Index = 0, Winning_index_2;
		double dist_1 = 0, dist_2 = 5;
		/*

		This for loop finds the minimum sum of average Connect #'s of the two objects; in ties preference is given to 
		averages with an input number in the y-direction (0,1,0)

		*/
		for(int i = 1; i < Maxs_1.size(); i++)
		{
			if(Maxs_1[i] + Maxs_2[i] <= Winning_Average) // && Maxs_1[i] >= 1 && Maxs_2[i] >= 1)
			{
				dist_1 = sqrt((double)(ALL_ROTATED_VECTORS_1[0]*ALL_ROTATED_VECTORS_1[0] + (ALL_ROTATED_VECTORS_1[1] - 1)*(ALL_ROTATED_VECTORS_1[1] - 1) + ALL_ROTATED_VECTORS_1[2]*ALL_ROTATED_VECTORS_1[2] ));
				dist_1 = dist_1 + sqrt((double)(ALL_ROTATED_VECTORS_2[0]*ALL_ROTATED_VECTORS[0] + (ALL_ROTATED_VECTORS_2[1] - 1)*(ALL_ROTATED_VECTORS_2[1] - 1) + ALL_ROTATED_VECTORS_2[2]*ALL_ROTATED_VECTORS_2[2] ));

				if(dist_1 < dist_2)
				{
					Winning_Average = Maxs_1[i] + Maxs_2[i];
					Winning_Index = i;
					dist_2 = dist_1;
				}
			}
		}


		std::ofstream myfile("Comparisons.txt",std::ios::app);
		myfile << "\n" <<  Run_Name << " average Slice-Connectedness " << Maxs_1[Winning_Index] << " and " << "\n" <<  Run_Name_2 << "average Slice-Connectedness  " <<  Maxs_2[Winning_Index] << "\n";
		std::cout << Run_Name << " average Slice-Connectedness " << Maxs_1[Winning_Index] << " and " << "\n" <<  Run_Name_2 << "average Slice-Connectedness  " <<  Maxs_2[Winning_Index] << "\n";



		myfile << " Average Maxs_1 " ;
		for(int i = 0; i < Maxs_1.size(); i++)
		{
			myfile << Maxs_1[i] << " " ;
		}
		myfile << "\n";
		myfile << " Average Maxs_2 " ;
		for(int i = 0; i < Maxs_2.size(); i++)
		{
			myfile << Maxs_2[i] << " " ;
		}
		myfile << "\n";
	
	


	std::vector<std::vector<std::vector<std::vector<int>>>>	Faces_1, Faces_2;

	// Now that the Winning Index is found we reconstruct the Mother Planes associated to it for each Object.

	std::cout << "Reconstructing the Winning Mother Planes for both Objects. " << "\n";
	for(int PCA_index = 0; PCA_index <= 1; PCA_index++)
	{
		std::vector<int> Projection_Ordering; 
		/*

		ii)  For each of the Input Points we then find Planes that: 
			 a)  Have this point as a normal vector and
			 b)  Stem from the object's center of mass (xcenter, ycenter, zcenter)

 		iii) Because we use MaxComp2D.cpp which takes as input a z-Plane, and because we wish to visualize the results 
			(e.g.  points collected in each Mother Plane) we make each Plane thick initially so that an injective map 
			 to the zcenter-plane  may be calculated.  For each input vector n its corresponding Plane and injective
			 map to the zcenter plane is found as follows

			a)  Each Plane is foound by searching over the whole image for pixels (x,y,z) satisfying
			|n1(x - xcenter) + n2(y - ycenter) + n3(z - zcenter)| < epsilon where n = (n1, n2, n3) and epsilon is 
			chosen large enough to leave the resulting Plane, which we call the 'Mother Plane' for n, a few pixels 'thick'.
			b) Using the ratios of the largest x/y/z to smallest x/y/z pixel in the thick plane we determine which primary
			plane (a x-plane, y-plane, or z-plane) maps injectively to the Mother Plane.  
			c)  We 'thin' the thick Mother Plane and find the injective map by projecting the proper primary plane
			pixel-by-pixel onto our Mother Plane, keeping only the first pixel of the Plane intersected by each pixel
			projection from the proper primary plane.  Below is a simple 2-D illustration:

			a) The thick Mother Plane is generated

								\\\                    
								 \\\
								  \\\
								   \\\
									\\\
									 \\\
							 Thick Mother Plane

		*/

		// test_* is the * component of the PCA'th Input Point.
		double test_x = 0, test_y = 0, test_z = 0;
		if(PCA_index == 0)
		{
		  test_x = ALL_ROTATED_VECTORS_1[3*Winning_Index];
		  test_y = ALL_ROTATED_VECTORS_1[3*Winning_Index + 1]; 
		  test_z = ALL_ROTATED_VECTORS_1[3*Winning_Index + 2];
		}
		else // PCA_index = 1
		{
		  test_x = ALL_ROTATED_VECTORS_1[3*Winning_Index];
		  test_y = ALL_ROTATED_VECTORS_1[3*Winning_Index + 1]; 
		  test_z = ALL_ROTATED_VECTORS_1[3*Winning_Index + 2];
		}


		std::cout << "Generating the Winning thick Mother Plane for the " << PCA_index + 1 << "th Input Point for Object " << Ford << "." << "\n";

		// low_*, high_* are the lowest/highest * components of the Plane.
		int low_x = xmax, high_x = 0, low_y = ymax, high_y = 0, low_z = zmax, high_z = 0;
		for(int z = 0; z < zmax; z++)
			for(int y = 0; y < ymax; y++)
				for(int x = 0; x < xmax; x++)
				{
					if( fabs( (double)test_x*(double)(x - xcenter )/(double)xmax  + test_y*(double)(y - ycenter)/(double)ymax + test_z*(double)(z - zcenter)/(double)zmax) < .009)
					{
						if( x < low_x)
						{
							low_x = x;
						}
						if(x > high_x)
						{
							high_x = x;
						}
						if( y < low_y)
						{
							low_y = y;
						}
						if( y > high_y)
						{
							high_y = y;
						}
						if( z < low_z )
						{
							low_z = z;
						}
						if( z > high_z )
						{
							high_z = z;
						}

						Ray_Image->FastSet(x,y,z,1000);
					}
				}

 

	/*

		

	iii)

		b) The proper plane is found for the injective map to the zcenter plane.

								|			\\\                    
								|			 \\\
								|			  \\\
								|			   \\\
								|               \\\
								|                \\\
							Primary plane     Thick Mother Plane

		c) The proper primary plane is projected onto the thick MotherPlane. Keeping only the first pixel of intersection 
		of the projection, and storing the coordinates of said pixel, gives the map and the thin  Mother Plane.
		 
								  Projection
								|   ----->	\\\                    
								| 	----->	 \\\
								|	----->	  \\\
								|	----->	   \\\
								|   ----->      \\\
								|   ----->       \\\
							Primary plane    Thick Mother Plane

							  <map>			\                   
							  <map>			 \
							  <map>			  \
							  <map>		       \
							  <map>             \
							  <map>              \
						   Injective Map   Thinned Mother Plane
							(a matrix)

		The Map can then be treated as one mapping to the zcenter-plane (by searching it properly).

							 <map> <map> <map> <map> <map> <map>
	*/


		std::cout << "Thinning the Winning Mother Plane for the " << PCA_index + 1 << "\n";
		std::cout <<"th Input Point and creating this Mother Plane's Injective Map to the zcenter-plane for Object " << Ford << "." << "\n";

		// 'Mother_Plane' is a matrix which will contain a single Mother Plane, and is input to Face_Planes.
		std::vector<std::vector<std::vector<int>>> Mother_Plane; //Projection_Map;
		std::vector<std::vector<int>> Mother_Plane_Row;  
		std::vector<int> Mother_Plane_Point;  

		int no_run_twice = 0; // This is to assure we don't project two primary planes onto the Mother Plane.
	
		// What if high_x = xmax and high_y = ymax && high_z = zmax?  We revert to "map_z" case below, put
		// first pruposefully because in most of our cases xmax = ymax > zmax.  If working with different 
		// sized images the cases map_z/x/y will need to be re-ordered (or a general case made).
		//if( high_x == xmax -1 && high_y == ymax -1)
		if(fabs(test_z) >= fabs(test_y) && fabs(test_z) >= fabs(test_x))
		{
			Projection_Ordering.push_back(0);
			no_run_twice++;
			for(int y = 0; y < ymax; y++)
				for(int x = 0; x < xmax; x++)
				{
					for(int z = 0; z < zmax; z++)
					{
							if(Ray_Image->FastGet(x,y,z) != 0)
							{
								if( Major_Axis_Plane_Projection->FastGet(x,y,zcenter) != 0 )
								{
									Ray_Image->FastSet(x,y,z,0);
								}
								else
								{

									Major_Axis_Plane_Projection->FastSet(x,y,zcenter,Ray_Image->FastGet(x,y,z)); //  Projection of x-y top plane onto our Mother Plane.
									
									Mother_Plane_Point.push_back(x);
									Mother_Plane_Point.push_back(y);
									Mother_Plane_Point.push_back(z);

									Mother_Plane_Row.push_back(Mother_Plane_Point);

									Mother_Plane_Point.clear();
								
								}
							}

					}
					if(x == xmax - 1) 
					{
						Mother_Plane.push_back(Mother_Plane_Row);

						Mother_Plane_Row.clear();
					}
				}
		}

	

	
		if(fabs(test_y) >= fabs(test_x) && fabs(test_y) >= fabs(test_z) && no_run_twice == 0)	
		{
			Projection_Ordering.push_back(1);
			no_run_twice++;
			for(int z = 0; z < zmax; z++)
				for(int x = 0; x < xmax; x++)
				{
					for(int y = 0; y < ymax; y++)
					{
						if(Ray_Image->FastGet(x,y,z) != 0)
							{
								if( Major_Axis_Plane_Projection->FastGet(x,ycenter,z) != 0 )
								{
									Ray_Image->FastSet(x,y,z,0);
								}
								else
								{

									Major_Axis_Plane_Projection->FastSet(x,ycenter,z,Ray_Image->FastGet(x,y,z));  //  Projection of xz top plane onto our Mother Plane.
									
									Mother_Plane_Point.push_back(x);
									Mother_Plane_Point.push_back(y);
									Mother_Plane_Point.push_back(z);

									Mother_Plane_Row.push_back(Mother_Plane_Point);

									Mother_Plane_Point.clear();

								}
						}

					}
					if(x == xmax - 1) 

					{
						Mother_Plane.push_back(Mother_Plane_Row);

						Mother_Plane_Row.clear();

					}		
							
						
				}
		}


		if( fabs(test_x) >= fabs(test_y) && fabs(test_x) >= fabs(test_z) &&  no_run_twice == 0)
		{
			Projection_Ordering.push_back(2);
			no_run_twice++;
			for(int z = 0; z < zmax; z++)
				for(int y = 0; y < ymax; y++)
				{
					for(int x = 0; x < xmax; x++)
					{
						if(Ray_Image->FastGet(x,y,z) != 0)
							{
									if( Major_Axis_Plane_Projection->FastGet(xcenter,y,z) != 0 )
									{
										Ray_Image->FastSet(x,y,z,0);
									}
									else
									{

										Major_Axis_Plane_Projection->FastSet(xcenter,y,z,Ray_Image->FastGet(x,y,z));  // Projection of yz top plane onto our Mother Plane.
										
										Mother_Plane_Point.push_back(x);
										Mother_Plane_Point.push_back(y);
										Mother_Plane_Point.push_back(z);

										Mother_Plane_Row.push_back(Mother_Plane_Point);

										Mother_Plane_Point.clear();

									}
							}

					}

						if(y == ymax - 1) 
						{
							Mother_Plane.push_back(Mother_Plane_Row);

							Mother_Plane_Row.clear();

						}

				}
		}
		

		// Because each Mother Plane was thick to begin wtih, and was then thinned, no point of the thin Mother Plane
		// lies on the Object's centroid (which we want).  To correct for this error we simply drag the entire Mother
		// Plane to the Object's centroid by the very middle pixel of the Mother Plane.

		// row/column give us the very middle of the Mother Plane.
		int row = (int)((double)(Mother_Plane.size())/2.0);
		int column = (int)((double)(Mother_Plane[row].size())/2.0);

		int xt = Mother_Plane[row][column][0];
		int yt = Mother_Plane[row][column][1];
		int zt = Mother_Plane[row][column][2];
		
		// This nested for-loop drags the Mother Plane by its center pixel to the Object's centroid.
		for(int i = 0; i < Mother_Plane.size(); i++)
			for(int j = 0; j < Mother_Plane[i].size(); j++)
			{
				Mother_Plane[i][j][0] = Mother_Plane[i][j][0] - xt + xcenter;
				Mother_Plane[i][j][1] = Mother_Plane[i][j][1] - yt + ycenter;
				Mother_Plane[i][j][2] = Mother_Plane[i][j][2] - zt + zcenter;
			} 


		if(PCA_index == 0)
		{
			Faces_1.push_back(Mother_Plane);
		}
		else // PCA_index = 1
		{
			Faces_2.push_back(Mother_Plane);
		}

		// Reset the Mother Plane Matrix.
		Mother_Plane.clear();
 
		// Reset images used
		for(int z = 0; z < zmax; z++)
			for(int y = 0; y < ymax; y++)
				for(int x = 0; x < xmax; x++)
				{	
					Ray_Image->FastSet(x,y,z,0);
					Major_Axis_Plane_Projection->FastSet(x,y,z,0);
				}


		}
	
	ALL_ROTATED_VECTORS_1.clear();
	ALL_ROTATED_VECTORS_2.clear();


	// Here we pick the Winning Daughter Planes, Projection Orderings, etc. 
	// Now we go through our vectors containing the Mother Planes, Highs and Lows, and Projection Orderings and choose
	// the 'Winning Index'th entry in each, and copy them into new vectors used below.
	std::vector<int> Temp_Proj_1, Temp_Proj_2;  // Will contain the projections (e.g. which face the Mother Plane is projected to)
	std::vector<std::vector<int>> Temp_HL_1, Temp_HL_2; // Will contain the Highs and Lows, e.g. the translation factor giving the tightest top and bottom Mother Plane for each object
	std::vector<std::vector<std::vector<int>>> Temp_Norm_1, Temp_Norm_2; // Will contain the 'Winning Index'th Mother Plane Normals

	std::vector<int> Point;
	std::vector<std::vector<int>> Line;


	Temp_Proj_1.push_back(Projections_1[Winning_Index]);
	Temp_Proj_2.push_back(Projections_2[Winning_Index]);

	Projections_1.clear();
	Projections_2.clear();

	Point.push_back(Highs_Lows_1[Winning_Index][0]);
	Point.push_back(Highs_Lows_1[Winning_Index][1]);
	Temp_HL_1.push_back(Point);
	Point.clear();
	Highs_Lows_1.clear();

	Point.push_back(Highs_Lows_2[Winning_Index][0]);
	Point.push_back(Highs_Lows_2[Winning_Index][1]);
	Temp_HL_2.push_back(Point);
	Point.clear();
	Highs_Lows_2.clear();

	for(int p = 0; p < Normals_1[Winning_Index].size(); p++)
	{
		Point.push_back(Normals_1[Winning_Index][p][0]);
		Point.push_back(Normals_1[Winning_Index][p][1]);
		Point.push_back(Normals_1[Winning_Index][p][2]);

		Line.push_back(Point);

		Point.clear();
		 
	}
		
		Temp_Norm_1.push_back(Line);
		Line.clear();
		Normals_1.clear();

	for(int p = 0; p < Normals_2[Winning_Index].size(); p++)
	{
		Point.push_back(Normals_2[Winning_Index][p][0]);
		Point.push_back(Normals_2[Winning_Index][p][1]);
		Point.push_back(Normals_2[Winning_Index][p][2]);

		Line.push_back(Point);

		Point.clear();
		
	}

	
	// Setting to 1000 for easier visulaization in Volview
	Temp_Norm_2.push_back(Line);
		Normals_2.clear();
		Line.clear();

		for(int x = 0; x < xmax_1; x++)
			for(int y = 0; y < ymax_1; y++)
				for(int z = 0; z < zmax_1; z++)
				{
					if(Second_Skinned_Image->FastGet(x,y,z) != 0)
					{
						Second_Skinned_Image->FastSet(x,y,z,1000);
					}
				}
		for(int x = 0; x < xmax_2; x++)
			for(int y = 0; y < ymax_2; y++)
				for(int z = 0; z < zmax_2; z++)
				{
					if(Skinned_Image->FastGet(x,y,z) != 0)
					{
						Skinned_Image->FastSet(x,y,z,1000);
					}
				}

		Winning_Index = 0;

		// ANNOTATE ALL INPUTS HERE
		// 'Comparison' function parameterizes each Intersection Plane
		std::cout << "Beginning parameterization of Connection Planes for Object " << mop << "\n";
		Comparison(Winning_Index, Temp_Proj_1, Temp_HL_1, Temp_Norm_1,  Faces_1,Second_Skinned_Image ,Second_UnSkinned_Image , Ray_Image, Intersection_Planes, Major_Axis_Plane_Projection);
		std::cout << "Beginning parameterization of Connection Planes for Object " << mop << "\n";	
		Comparison(Winning_Index, Temp_Proj_2, Temp_HL_2, Temp_Norm_2, Faces_2, Skinned_Image, UnSkinned_Image, Ray_Image, Intersection_Planes, Major_Axis_Plane_Projection);

		
		/*

		Having found all of the Intersection Points and calculated the Feature differences at each of these points, we now 
		calculate the global difference between the two surfaces and map local differences simultaneously back to both Objects' surfaces.

		*/

		// ALL_POINTS_1/2 contains all Intersection Points from the parameterization in 'Comparison' for Object 1/2
		// ALL_MEAS_1/2 contains all of the Seed measurements for all Intersection Points from the parameterization in 'Comparison' for Object 1/2
		// ALL_MEAS contains |ALL_MEAS_1 - ALL_MEAS_2|
		// We first calculate all of the Local Differences.

		int meas_counter = 0; // The counter for the total number of Intersection Points on the two Objects.
		double total_diff = 0;
		for(int i = 0; i < 30; i++)  // Slice intersected with Plane 
		{
			
			for(int h = 0; h <= 1; h++) {// Planes
				
				for(int j = 0; j < 30; j++) { // Points in each of those slices
					for(int k = 0; k <= 1; k++) // Min/Max points
					{
						if( (ALL_POINTS_1[h][i][j][k][0] != 0 && ALL_POINTS_1[h][i][j][k][1] != 0 && ALL_POINTS_1[h][i][j][k][2] != 0) || ( ALL_POINTS_2[h][i][j][k][0] != 0 && ALL_POINTS_2[h][i][j][k][1] != 0 && ALL_POINTS_2[h][i][j][k][2] != 0 ) ) 
						{
							ALL_MEAS[h][i][j][k] =   fabs( ALL_MEAS_1[h][i][j][k] - ALL_MEAS_2[h][i][j][k] );
							
							total_diff = total_diff + ALL_MEAS[h][i][j][k];

							meas_counter++;
							
						}
						i = i;
					}

				//std::cout<<j<<". Plane_diff: "<<total_diff<<"\n";		
				}
			}
			
		}
		
	 //std::cout<<"total_diff: "<<total_diff<<"\n";

		// Now we find the max and min local Feature differences between the two Objects' surfaces.
		double max = 0, min = 1;
		for(int i = 0; i < 30; i++)  // Slice intersected with Plane 
		{
			for(int h = 0; h <= 1; h++) // Planes
				for(int j = 0; j < 30; j++)  // Points in each of those slices
					for(int k = 0; k <= 1; k++) // Min/Max points
					{
						if( (ALL_POINTS_1[h][i][j][k][0] != 0 && ALL_POINTS_1[h][i][j][k][1] != 0 && ALL_POINTS_1[h][i][j][k][2] != 0) || ( ALL_POINTS_2[h][i][j][k][0] != 0 && ALL_POINTS_2[h][i][j][k][1] != 0 && ALL_POINTS_2[h][i][j][k][2] != 0 ) ) 
						{
							ALL_MEAS[h][i][j][k] =   fabs( ALL_MEAS_1[h][i][j][k] - ALL_MEAS_2[h][i][j][k] );
						
							if(max < ALL_MEAS[h][i][j][k])
							{
								max = ALL_MEAS[h][i][j][k];
							}
							if(min > ALL_MEAS[h][i][j][k])
							{
								min = ALL_MEAS[h][i][j][k];
							}
						}
					}
		}


		double cutoff = (max - min);
		double partition = cutoff/15.0; // This dtermines the number of color-ranks local differences are fitted into below.



		total_diff = sqrt(total_diff)/(double)(meas_counter);
		myfile << "Total Diff " << total_diff <<  " Max point diff " << max << " and Total Points " << meas_counter << "\n" << "\n" ;
		myfile << " Radi_1 are " << Final_1_1 << " " << Final_2_1 << " " << Final_3_1 << "\n";
		myfile << " Radi_2 are " << Final_1_2 << " " << Final_2_2 << " " << Final_3_2 << "\n";
		myfile.close();


		std::cout << "Mapping local differences back to both Objects' surfaces. " << "\n";


		/*

		 Now we map local differences simultaneously back onto the two Objects' surfaces.  Differences are mapped back 
		 from largest to smallest, with larger differences taking precedence on the Objects' surfaces, and by color-rank.
		
		 */
		int h_hold = 0, i_hold = 0, j_hold = 0, k_hold = 0;
		double max_now = 0;
		int all_meas = 0;
		if(cutoff != 0)
		{
			while(all_meas <= (double)(meas_counter) ) // you can control what portion of the local differences you wish to map onto both surfaces by increasing/decreasing cut_off_portion (currently set to 1; we map back all differences).
			{
				// This nested for loop finds the next biggest local difference.
				for(int i = 0; i < 30; i++)  // Slice intersected with Plane 
				{
					for(int h = 0; h <= 1; h++) // Planes
						for(int j = 0; j < 30; j++)  // Points in each of those slices
							for(int k = 0; k <= 1; k++) // Min/Max points
							{
								if(ALL_MEAS[h][i][j][k] > max_now)
								{
									max_now = ALL_MEAS[h][i][j][k];
									h_hold = h;
									i_hold = i;
									j_hold = j;
									k_hold = k;
								}
							}
				}

				//  This is the point on the Second Object corresponding to this max difference.
				int x = ALL_POINTS_2[h_hold][i_hold][j_hold][k_hold][0];		
				int y = ALL_POINTS_2[h_hold][i_hold][j_hold][k_hold][1];
				int z = ALL_POINTS_2[h_hold][i_hold][j_hold][k_hold][2];
				
				// This next while loop fits our chosen difference into the right color-rank  (of 15 at present).
				double crawl = 0;
				  int crawl_count = -1;
				  while(max_now >= min + crawl)
				  {
					  crawl_count++;
					  crawl = crawl + partition;
				  }

			   // This nested for loop maps the chosen difference back onto the first Object's surface.
				for(int i = x - Final_1_2; i <= x + Final_1_2; i++)
					for(int j = y - Final_2_2; j <= y + Final_2_2; j++)
						for(int k = z - Final_3_2; k <= z + Final_3_2; k++)
						{
							if( Second_Skinned_Image->FastGet(i,j,k) == 1000)
							{
								Second_Skinned_Image->FastSet(i,j,k,5000 + crawl_count*1000);
							}
						}
						  
				  x = ALL_POINTS_1[h_hold][i_hold][j_hold][k_hold][0];		
				  y = ALL_POINTS_1[h_hold][i_hold][j_hold][k_hold][1];
				  z = ALL_POINTS_1[h_hold][i_hold][j_hold][k_hold][2];

				// This nested for loop maps the chosen difference back onto Object 2's surface.
				for(int i = x - Final_1_1; i <= x + Final_1_1; i++)
					for(int j = y - Final_2_1; j <= y + Final_2_1; j++)
						for(int k = z - Final_3_1; k <= z + Final_3_1; k++)
						{
							if(  Skinned_Image->FastGet(i,j,k) == 1000)
							{
								 Skinned_Image->FastSet(i,j,k, 5000 + crawl_count*1000);
							}
	 
						}

					max_now = 0;
					ALL_MEAS[h_hold][i_hold][j_hold][k_hold] = 0; // We erase the chosen difference, allowing us to find the next biggest one on the following iteration.
					all_meas++;	
				
			}
		}

		// Finally, we print out all local differences mapped back onto the two surfaces.

		Write_Analyze_File((char *)(new std::string(Run_Name + (*(new std::string("_HIGHS_MARKED.hdr")))))->c_str(), *Second_Skinned_Image);
		Write_Analyze_File((char *)(new std::string(Run_Name_2 + (*(new std::string("_HIGHS_MARKED.hdr")))))->c_str(), *Skinned_Image);

		int pooooooo = 0;

	}

 
}



void Grid_Line::Comparison(int Winning_Index, std::vector<int> Projection_Ordering, std::vector<std::vector<int>> Highs_Lows, std::vector<std::vector<std::vector<int>>> Normals, std::vector<std::vector<std::vector<std::vector<int>>>> Face_Planes,  CIS_Array_Image3D_short * Skinned_Image, CIS_Array_Image3D_short * UnSkinned_Image, CIS_Array_Image3D_short * Intersection_Points, CIS_Array_Image3D_short * Intersection_Planes,CIS_Array_Image3D_short * Major_Axis_Plane_Projection )
	{


	// Sets image max variables for eeach Object
	if(mop == 1)
	{
		xmax = xmax_1;
		ymax = ymax_1;
		zmax = zmax_1;
	 
	}
	if(mop == 2)
	{
		xmax = xmax_2;
		ymax = ymax_2;
		zmax = zmax_2;
 
	}

	/*

		 v)    We calculate the Objects' Surface Feature at each Intersection Point.  The current feature, we call
		       the Manay Feature, is the volume of intersection of a sphereoid/parallelepiped (or any object really),
			   called the Seed here, intersected with the interior of Object.  
			   a)  We define the axes of Seed by d/m, the distance via Plane Normal between two neighboring Intersector
			       Planes.
			   b)  We normalize the volume of the Seed by dividing by the Seed's full volume.

		 
		 vi)   We compare the Manay Feature of each point on corresponding Connection Planes, and map the difference
			   back onto each Object's surface.

		*/

	int z_top = 0, z_bot = 0, x_top = 0, x_bot = 0, y_top = 0, y_bot = 0;

	std::vector<std::vector<std::vector<int>>> Projected_Normals; // used as above, only for each bounding rectangle
	std::vector<std::vector<int>> Temp_Normal_1, Temp_Highs_Lows; // Will contain the two lines used in constructing side 1 and side 2 of each bounding rectangle
	std::vector<int> Points;

	/*
		ii)  For each Winning Plane we find the bounding square containing all (hopefully one) components where the 
			 sides of said square are ortogonal to the Normal Vector for the Plane examined.


 							   |   |___|  |
							   |  /|   |\ |
							   |-->|-->|->|  
							   |/__|___|_\|
							   |   |   |  |
				The Object intersected with user-defined number of Winning Planes

										   ______________________
					|				      |					     |
					|				      |				 ___	 |
					|				      |		 __		/   \    |
					|				      |	    /  \___/ 	|	 |
					|  -------------->    |	   /		   /     |
					|				      |	  |		      |	     |
					|				      |	   \		 /	     |
					|			          |		\_______/	     |
					|				      |	    				 |
					|				      |______________________|
				 Intersection         Intersection Plane Rotated parallel 
				   Plane          to this window,with 2D Component of Object

	*/

	/*
	 Here we strip off the two orthogonal vectors passing through (xcenter,ycenter,zcenter) from the Mother Plane to 
	 be used in parameterizing each of its Connection Planes.  These are found by simply by taking the middle row and column
	 of 'Face_Planes'; the Matrix containing the Mother Plane.  With these two vectors we make the bounding rectangle
	 containing all of the Object in each Intersection Plane.
	*/

	// The sides of the bounding rectangles are taken from the middle row and middle column of the matrix representing the
	// Mother Plane.  They are then translated to produce each bounding rectangle on each Connection Plane.
	int rows = (int)((double)(Face_Planes[Winning_Index].size())/2.0);
	int columns = (int)((double)(Face_Planes[Winning_Index][rows].size())/2.0);

	// The line taken from the middle row of the Mother Plane matrix (Face Planes)
  	for(int j = 0; j < Face_Planes[Winning_Index][rows].size(); j++)
	{
		Points.push_back(Face_Planes[Winning_Index][rows][j][0]);
		Points.push_back(Face_Planes[Winning_Index][rows][j][1]);
		Points.push_back(Face_Planes[Winning_Index][rows][j][2]);

		Temp_Normal_1.push_back(Points);

		Points.clear();

	}

	int start = 0;
	while(  Face_Planes[Winning_Index][start][columns].empty() )
	{	
		start++;
	}

	int end = Face_Planes[Winning_Index].size() - 1;
	while(  Face_Planes[Winning_Index][end][columns].empty() )
	{	
		end--;
	}
 
	std::vector<std::vector<int>> Temp_Normal_2;

	// The line taken from the middle column of the Mother Plane matrix (Face Planes)
 	for(int j = start; j <= end; j++)
	{
		Points.push_back(Face_Planes[Winning_Index][j][columns][0]);
		Points.push_back(Face_Planes[Winning_Index][j][columns][1]);
		Points.push_back(Face_Planes[Winning_Index][j][columns][2]);

		Temp_Normal_2.push_back(Points);

		Points.clear();

	}

	Projected_Normals.push_back(Temp_Normal_1);
	Projected_Normals.push_back(Temp_Normal_2);


	int top_size_1 = Temp_Normal_1.size() - 1;
	int bot_size_1 = 0;
	int top_size_2 = Temp_Normal_2.size() - 1;
	int bot_size_2 =  0;

	// We begin a while loop generating the Winning Connection Planes, inside it will also parameterize each slice.
	// But notice getting the Connection Planes themselves is exactly as done above.  
	double top_bot_partition = (double)(Highs_Lows[Winning_Index][0] - Highs_Lows[Winning_Index][1])/slice_gap;
	int partition_num = 0;
 
	double extra = top_bot_partition/(0.5*slice_gap);
	double extra_2 = extra;
	if(extra < 2)
	{
		extra_2 = 2;
	}

	double JUMP = (double)Highs_Lows[Winning_Index][1] + extra_2;

	Final_Radius = (int)(floor(sqrt( (double) ( ( Normals[Winning_Index][(int)floor(JUMP + 0.5)][0]- 
		Normals[Winning_Index][(int)floor(JUMP + top_bot_partition + 0.5)][0])*(Normals[Winning_Index][(int)floor(JUMP + 0.5)][0]- 
		Normals[Winning_Index][(int)floor(JUMP + top_bot_partition + 0.5)][0] ) + ( Normals[Winning_Index][(int)floor(JUMP + 0.5)][1] -
		Normals[Winning_Index][(int)floor(JUMP + top_bot_partition + 0.5)][1])*(Normals[Winning_Index][(int)floor(JUMP + 0.5)][1] - 
		Normals[Winning_Index][(int)floor(JUMP + top_bot_partition + 0.5)][1] ) + ( Normals[Winning_Index][(int)floor(JUMP + 0.5)][2] - 
		Normals[Winning_Index][(int)floor(JUMP + top_bot_partition + 0.5)][2])*(Normals[Winning_Index][(int)floor(JUMP + 0.5)][2] - 
		Normals[Winning_Index][(int)floor(JUMP + top_bot_partition + 0.5)][2] ) ) ) + 0.5 ) );
	Final_Radius = Final_Radius;

	int level = 0;

	if(Projection_Ordering[Winning_Index] == 0)
	{
		level = zcenter;
	}
	if(Projection_Ordering[Winning_Index] == 1)
	{
		level = ycenter;
	}			
	if(Projection_Ordering[Winning_Index] == 2)
	{				
		level = xcenter;
	}
	// These next if loops record the Seed Radi for each Object.
	if(mop == 1)
	{
		Final_Radius_1 = Final_Radius;
	}
	else // mop == 2
	{
		Final_Radius_2 = Final_Radius;
	}

	double max_pix = xpixel;
	if( max_pix < ypixel)
	{
		max_pix = ypixel;
	}
	if(max_pix < zpixel)
	{
		max_pix = zpixel;
	}

	if(mop == 1)
	{

  		Final_1_1 = (int)floor( xpixel/max_pix*( (double)(Final_Radius)/2.0) + 0.5);
		Final_2_1 = (int)floor( ypixel/max_pix* ( (double)(Final_Radius)/2.0) + 0.5);
		Final_3_1 = (int)floor( zpixel/max_pix*( (double)(Final_Radius)/2.0) + 0.5);

		if(Final_1_1 < 1)
		{
			Final_1_1 = 1;
		}
		if(Final_2_1 < 1)
		{
			Final_2_1 = 1;
		}	
		if(Final_3_1 < 1)
		{
			Final_3_1 = 1;
		}
 


		Final_1 = Final_1_1;
		Final_2 = Final_2_1;
		Final_3 = Final_3_1;
	}

	if(mop == 2)
	{
 		Final_1_2 = (int)floor( xpixel/max_pix*( (double)(Final_Radius)/2.0) + 0.5);
		Final_2_2 = (int)floor( ypixel/max_pix* ( (double)(Final_Radius)/2.0) + 0.5);
		Final_3_2 = (int)floor( zpixel/max_pix*( (double)(Final_Radius)/2.0) + 0.5);



		if(Final_1_2 < 1)
		{
			Final_1_2 = 1;
		}
		if(Final_2_2 < 1)
		{
			Final_2_2 = 1;
		}	
		if(Final_3_2 < 1)
		{
			Final_3_2 = 1;
		}

		Final_1 = Final_1_2;
		Final_2 = Final_2_2;
		Final_3 = Final_3_2;
	}

	// The full volume of the Seed is calculated (in this case the Seed is a parallelapiped.
	int ball_volume = 0;
	for(int i = xcenter - Final_1; i <= xcenter + Final_1; i++)
		for(int j = ycenter - Final_2; j <= ycenter + Final_2; j++)
			for(int k = zcenter - Final_3; k <= zcenter + Final_3; k++)
			{
				//Billy attempting to use sphere seed
				//double dist = (i - xcenter)*(i - xcenter)  + (j - ycenter)*(j - ycenter) + (k - zcenter)*(k - zcenter);
				//if(dist < Final_Radius*Final_Radius)
				{
					ball_volume++;
				}
			}




	top_bot_partition = (double)( (Highs_Lows[Winning_Index][0] - extra_2 - (Highs_Lows[Winning_Index][1] + extra_2) ) )/slice_gap;
	
	if(top_bot_partition < 1)
	{
		top_bot_partition = 1;
	}


		// For visualization purposes we copy the Object Skin to Intersection_Planes_Image.
		for(int x = 0; x < xmax; x++)
			for(int y = 0; y < ymax; y++)
				for(int z = 0; z < zmax; z++)
				{
					Intersection_Planes->FastSet(x,y,z,Skinned_Image->FastGet(x,y,z));
				}
			
	while(JUMP < Highs_Lows[Winning_Index][0] )// + extra_2/2.0 )
	{
		std::cout << "Generating and Parameterizing " << partition_num + 1 << "the Connection Planes for Object " << mop << ".\n";

		int part = (int)floor(JUMP + 0.5);
		int loooook = 0;
		// This prints the Chosen Collector Planes
		for(int j = 0; j < Face_Planes[Winning_Index].size(); j++)
			for(int k = 0; k < Face_Planes[Winning_Index][j].size(); k++)
			{
				 
				if( UnSkinned_Image->FastGet(Face_Planes[Winning_Index][j][k][0] + Normals[Winning_Index][part][0] - xcenter , Face_Planes[Winning_Index][j][k][1] + Normals[Winning_Index][part][1] - ycenter, Face_Planes[Winning_Index][j][k][2] +  Normals[Winning_Index][part][2] - zcenter) != 0)
				{
					Intersection_Planes->FastSet(Face_Planes[Winning_Index][j][k][0] + Normals[Winning_Index][part][0] - xcenter , Face_Planes[Winning_Index][j][k][1] + Normals[Winning_Index][part][1] - ycenter, Face_Planes[Winning_Index][j][k][2] +  Normals[Winning_Index][part][2] - zcenter, (partition_num + 2)*1000);
 						
					if(Projection_Ordering[Winning_Index] == 0)
					{
						Major_Axis_Plane_Projection->FastSet(Face_Planes[Winning_Index][j][k][0] + Normals[Winning_Index][part][0] - xcenter ,Face_Planes[Winning_Index][j][k][1] + Normals[Winning_Index][part][1] - ycenter , level, 500);
					}
					if(Projection_Ordering[Winning_Index] == 1)
					{
						Major_Axis_Plane_Projection->FastSet(Face_Planes[Winning_Index][j][k][0] + Normals[Winning_Index][part][0] - xcenter , level, Face_Planes[Winning_Index][j][k][2] + Normals[Winning_Index][part][2] - zcenter, 500);
					}
					if(Projection_Ordering[Winning_Index] == 2)
					{
						Major_Axis_Plane_Projection->FastSet(level, Face_Planes[Winning_Index][j][k][1] + Normals[Winning_Index][part][1] - ycenter, Face_Planes[Winning_Index][j][k][2] + Normals[Winning_Index][part][2] - zcenter, 500);
					}

				
				}
			}


			std::vector<std::vector<int>> New_Highs_Lows;
			std::vector<int> Temp;


			 top_size_1 = Temp_Normal_1.size() - 1;
			 bot_size_1 = 0;
			 top_size_2 = Temp_Normal_2.size() - 1;
			 bot_size_2 =  0;



			Temp.push_back(0);
			Temp.push_back(0);
			New_Highs_Lows.push_back(Temp);
			New_Highs_Lows.push_back(Temp);
			New_Highs_Lows.push_back(Temp);
			Temp.clear();
	

		/*
			iii)   The tightest square containing all (hopefully only one) Object components with sides orthogonal
					to the Plane Normal is found.  These orthogonal vectors are taken from each Mother Plane.

						 ______________________  
						|					   |				  ________________
						|				 ___   |			  	 |           ___  |
						|		 __		/   \  |				 |	 __		/   \ |
						|	    /  \___/ 	|  |				 |  /  \___/ 	| |
						|	   /		   /   |	--------->	 | /		   /  |
						|	  |		      |	   |				 | |		  |   |	  
						|	   \		 /	   |				 | \		 /	  |
						|		\_______/	   |				 |  \_______/	  |
						|______________________|				  ----------------
					  Tightest square for Intersection Plane Object components if found

		*/
		/*
		 
		
		'top_size_1' and 'bot_size_1' are parameters of Normal_1 such that Normal_2 + Normal_1(top_size_1), which is Normal_2
		shifted 'up' by Normal_1(top_size_1), is the closest translated Normal_2 to the Object without intersecting it on one side
		of the Intersection Plane.  
		Likewise for Normal_2 + Normal_1(bot_size_1) etc.  
		A drawing may help exlpain this (I hope):
		*/
		/*
					
								Normal_1 + Normal_2(begin)
								|
								|					   |Normal_1 + Normal_2(end)
 Normal_2 + Normal_1(begin  )___| _____________________|					
								|					   |					  ________________________ Normal_2 +  Normal_1(top_size_1)
								|				 ___   |					 |           ___  |
								|		 __		/   \  |					 |	 __		/   \ |
								|	    /  \___/ 	|  |					 |  /  \___/ 	| |
								|	   /		   /   |		--------->	 | /		   /  |
								|	  |		      |	   |					 | |		  |   |	  
								|	   \		 /	   |					 | \		 /	  |
								|		\_______/	   |					 |  \_______/	  |
Normal_2 + Normal_1(end  )   ___|______________________|					 |----------------|------- Normal_2 +  Normal_1(bot_size_1)
																			 |                |
																			 |                |
																			 |			Normal_1 + Normal_2(bot_size_2)
																			 |Normal_1 + Normal_2(top_size_2))


		*/
			

			
			int xp = 0, yp = 0, zp = 0;
			int lounter = bot_size_2;  
			int sum = 0;
			int side = top_size_1 ; 
			int lintersection = 0;

			// This 'while' finds the 'top' of the Bounding Rectangle made by Normal_2 translated by Normal_1
			while( sum == 0)
			{
				while(lintersection == 0 && lounter < top_size_2)
				{
					xp = (Projected_Normals[1][lounter][0]) + (Normals[Winning_Index][part][0] - xcenter)  +  (Projected_Normals[0][side][0] - xcenter);
					yp = (Projected_Normals[1][lounter][1]) + (Normals[Winning_Index][part][1] - ycenter)  +  (Projected_Normals[0][side][1]  - ycenter);
					zp = (Projected_Normals[1][lounter][2]) + (Normals[Winning_Index][part][2] - zcenter)  +  (Projected_Normals[0][side][2]  - zcenter);

					lounter++;
					lintersection = lintersection + UnSkinned_Image->FastGet(xp,yp,zp);
 				}
				sum = lintersection;
				lintersection = 0;
				side--;
				lounter = bot_size_2;

				if(side == bot_size_1)
				{
					sum = 1;
				}

			}
			side++;
			sum = 0;

			Temp.push_back(side);
			

			lounter = bot_size_2;
			side = bot_size_1;
			 

			// This 'while' finds the 'bottom' of the Bounding Rectangle made by Normal_2 translated by Normal_1
			while( sum == 0)
			{
				while(lintersection == 0 && lounter < top_size_2)
				{
					xp = (Projected_Normals[1][lounter][0]) + (Normals[Winning_Index][part][0] - xcenter)  +  (Projected_Normals[0][side][0]  - xcenter);
					yp = (Projected_Normals[1][lounter][1]) + (Normals[Winning_Index][part][1] - ycenter)  +  (Projected_Normals[0][side][1]  - ycenter);
					zp = (Projected_Normals[1][lounter][2]) + (Normals[Winning_Index][part][2] - zcenter)  +  (Projected_Normals[0][side][2] - zcenter);

					lounter++;
					lintersection = lintersection + UnSkinned_Image->FastGet(xp,yp,zp);

				}
				sum = lintersection;
				lintersection = 0;
				side++;
				lounter = bot_size_2;

				if(side == top_size_1)
				{
					sum = 1;
				}

			}

			side--;

			Temp.push_back(side);

			if(Temp[0] > Temp[1])
			{
				top_size_1 = Temp[0];
				bot_size_1 = Temp[1];
			}
			else
			{
				top_size_1 = Temp[1];
				bot_size_1 = Temp[0];
			}

			Temp.clear();




		   sum = 0;
		   side = top_size_2;
		   lounter = bot_size_1;
		   lintersection = 0;
			// This 'while' finds the 'top' of the Bounding Rectangle made by Normal_1 translated by Normal_2
			while( sum == 0)
			{
				while(lintersection == 0 && lounter < top_size_1)
				{
					xp = (Projected_Normals[0][lounter][0]) + (Normals[Winning_Index][part][0] - xcenter)  +  (Projected_Normals[1][side][0]  - xcenter);
					yp = (Projected_Normals[0][lounter][1]) + (Normals[Winning_Index][part][1] - ycenter)  +  (Projected_Normals[1][side][1]  - ycenter);
					zp = (Projected_Normals[0][lounter][2]) + (Normals[Winning_Index][part][2] - zcenter)  +  (Projected_Normals[1][side][2]   - zcenter);

					lounter++;
					lintersection = lintersection + UnSkinned_Image->FastGet(xp,yp,zp);
 

				}
				sum = lintersection;
				lintersection = 0;
				side--;
				lounter = bot_size_1;

				if(side == bot_size_2)
				{
					sum = 1;
				}
			}
			side++;
			sum = 0;

			Temp.push_back(side);


			side = bot_size_2;
		    lounter = bot_size_1;
			// This 'while' finds the 'bottom' of the Bounding Rectangle made by Normal_1 translated by Normal_2
			while( sum == 0)
			{
				while(lintersection == 0 && lounter < top_size_1)
				{
					xp = (Projected_Normals[0][lounter][0]) + (Normals[Winning_Index][part][0] - xcenter)  +  (Projected_Normals[1][side][0]  - xcenter);
					yp = (Projected_Normals[0][lounter][1]) + (Normals[Winning_Index][part][1] - ycenter)  +  (Projected_Normals[1][side][1]  - ycenter);
					zp = (Projected_Normals[0][lounter][2]) + (Normals[Winning_Index][part][2] - zcenter)  +  (Projected_Normals[1][side][2]   - zcenter);

					lounter++;
					lintersection = lintersection + UnSkinned_Image->FastGet(xp,yp,zp);
 

				}
				sum = lintersection;
				lintersection = 0;
				side++;

				lounter = bot_size_1;

				if(side == top_size_2)
				{
					sum = 1;
				}

			}

			side--;

			Temp.push_back(side);

			if(Temp[0] > Temp[1])
			{
				top_size_2 = Temp[0];
				bot_size_2 = Temp[1];
			}
			else
			{
				top_size_2 = Temp[1];
				bot_size_2 = Temp[0];
			}
			Temp.clear();

 

			/*

			iv)    Sides of the bounding rectangle are partitioned a user-defined number of times, and points at 
				   these partitions on the rectangle are projected onto the Object component(s) (we call these projections
				   Intersection  Lines). The first point of intersection with the Object component(s) from each projection
				   is stored, we call these points Intersection Points. 

			
								|		|
							 ___|_______|____
							|   |       |__  |                              o
							|	|_	   /   \ |                      o 
						____|__/  \___/    |_|___                  o           o
							| /		      /  |   -------------->
						____|_|		     |___|___	              o          o
							| \		    /	 |                             o
							|  \_______/|	 |                       o
							 ---|-------|-----
								|       |
							Intersecting Lines				Points on the Object component	
		
			*/



			// top_bot_partition_3 is the distance between two opposite sides of the Bounding Rectangle.
			double top_bot_partition_3 = (double)(top_size_2 -  bot_size_2  )/line_gap_1;
			int partition_num_3 = 0;
		 

			if(top_bot_partition_3 < 2)
			{
				top_bot_partition_3 = 2;
			}
		 
			// We now produce Intersection Lines/Intersection Points across each of the four sides of the Bounding Rectangle, one
			// side at a time. 

			// SIDE 1.
			// As we did when selecting the Connection Planes, when finding the Intersecting Lines we first shift in slightly
			// by bonus/bonus_3 from each side.
			double bonus = top_bot_partition_3/(0.5*line_gap_1);
			double bonus_3 = bonus;
			if(bonus < 2)
			{
				bonus_3 = 2;
			}

			top_bot_partition_3 = ((double)top_size_2 - bonus_3 - (double)bot_size_2 - bonus_3)/line_gap_1;

			if(top_bot_partition_3 < 1)
			{
				top_bot_partition_3 = 1;
			}

			double JUMP_3 = bot_size_2  +  bonus_3;

			while(JUMP_3 < top_size_2)// - bonus_3)
			{
				int part_3 = (int)floor(JUMP_3 + 0.5);
				
				int x = 0, y = 0, z = 0;
				int intersection = 0;
				int counter = bot_size_1;

				// Here we project from the 'left' side to the 'right' side.
				while(intersection == 0 && counter < top_size_1)
				{
					x = (Projected_Normals[0][counter][0]) + (Normals[Winning_Index][part][0] - xcenter)  +  (Projected_Normals[1][part_3][0] - xcenter);
					y = (Projected_Normals[0][counter][1]) + (Normals[Winning_Index][part][1] - ycenter)  +  (Projected_Normals[1][part_3][1] - ycenter);
					z = (Projected_Normals[0][counter][2]) + (Normals[Winning_Index][part][2] - zcenter)  +  (Projected_Normals[1][part_3][2] - zcenter);

					counter++;
					intersection = UnSkinned_Image->FastGet(x,y,z);

				}



				if(intersection != 0)
				{
								// At the first intersection of the projection we print out the point to the vector, image, and projected image.

					Intersection_Points->FastSet(x,y,z,2000);
					
					// Load the Intersection Point into ALL_POINTS, our vector of Intersection Points.
					ALL_POINTS[0][partition_num][partition_num_3][0][0] = x;
					ALL_POINTS[0][partition_num][partition_num_3][0][1] = y;
					ALL_POINTS[0][partition_num][partition_num_3][0][2] = z;

					// Also, this point is projected onto an approprite primary plane for visualization purposes.
					if(Projection_Ordering[Winning_Index] == 0)
					{
						Major_Axis_Plane_Projection->FastSet(x ,y , level, 1000);
					}
					if(Projection_Ordering[Winning_Index] == 1)
					{
						Major_Axis_Plane_Projection->FastSet(x,level, z, 1000);
					}
					if(Projection_Ordering[Winning_Index] == 2)
					{
						Major_Axis_Plane_Projection->FastSet(level, y, z, 1000);
					}

				}


				// SIDE 2 (opposite of SIDE 1).
				x = 0, y = 0, z = 0;
				intersection = 0;
				counter = top_size_1;
				// Here we project from the 'right' side to the 'left' side
				while(intersection == 0 && counter > bot_size_1)
				{
					x = (Projected_Normals[0][counter][0]) + (Normals[Winning_Index][part][0] - xcenter)  +  (Projected_Normals[1][part_3][0] - xcenter);
					y = (Projected_Normals[0][counter][1]) + (Normals[Winning_Index][part][1] - ycenter)  +  (Projected_Normals[1][part_3][1] - ycenter);
					z = (Projected_Normals[0][counter][2]) + (Normals[Winning_Index][part][2] - zcenter)  +  (Projected_Normals[1][part_3][2] - zcenter);
			
					counter--;
					intersection = UnSkinned_Image->FastGet(x,y,z);
				}

				if(intersection != 0)
				{
					// At the first intersection of the projection we print out the point to the vector, image, and projected image.

					Intersection_Points->FastSet(x,y,z,2000);

					ALL_POINTS[0][partition_num][partition_num_3][1][0] = x;
					ALL_POINTS[0][partition_num][partition_num_3][1][1] = y;
					ALL_POINTS[0][partition_num][partition_num_3][1][2] = z;

					if(Projection_Ordering[Winning_Index] == 0)
					{
						Major_Axis_Plane_Projection->FastSet(x ,y , level, 1000);
					}
					if(Projection_Ordering[Winning_Index] == 1)
					{
						Major_Axis_Plane_Projection->FastSet(x,level, z, 1000);
					}
					if(Projection_Ordering[Winning_Index] == 2)
					{
						Major_Axis_Plane_Projection->FastSet(level, y, z, 1000);
					}

				}
				

				partition_num_3++;
				JUMP_3 = JUMP_3 + top_bot_partition_3;
			}



			// Now we find Intersection Lines/Intersection Points produced along the other two parallel sides of the Bounding Rectangle.
			top_bot_partition_3 = (double)( top_size_1 -  bot_size_1  )/line_gap_2;
			 partition_num_3 = 0;
		
			if(top_bot_partition_3 < 2)
			{
				top_bot_partition_3 = 2;
			}
		 
		 
			 bonus = top_bot_partition_3/(0.5*line_gap_2);
			 bonus_3 = bonus;
			if(bonus < 2)
			{
				bonus_3 = 2;
			}

			top_bot_partition_3 = ((double) top_size_1 - bonus_3 - (double) bot_size_1 - bonus_3)/line_gap_2;

			if(top_bot_partition_3 < 1)
			{
				top_bot_partition_3 = 1;
			}

			 JUMP_3 =  bot_size_1  +  bonus_3;

			while(JUMP_3 <  top_size_1)// - bonus_3)
			{
				int part_3 = (int)floor(JUMP_3 + 0.5);
 
		
				int x = 0, y = 0, z = 0;
				int intersection = 0;
				int counter =  bot_size_2;

				// In this while loop we project from the 'top' side to 'bottom' sdie
				while(intersection == 0 && counter <  top_size_2)
				{
					x = (Projected_Normals[1][counter][0]) + (Normals[Winning_Index][part][0] - xcenter)  +  (Projected_Normals[0][part_3][0] - xcenter);
					y = (Projected_Normals[1][counter][1]) + (Normals[Winning_Index][part][1] - ycenter)  +  (Projected_Normals[0][part_3][1] - ycenter);
					z = (Projected_Normals[1][counter][2]) + (Normals[Winning_Index][part][2] - zcenter)  +  (Projected_Normals[0][part_3][2] - zcenter);

					counter++;
					intersection = UnSkinned_Image->FastGet(x,y,z);
				}



				if(intersection != 0)
				{
					// At the first intersection of the projection we print out the point to the vector, image, and projected image.
					Intersection_Points->FastSet(x,y,z,2000);
				

					ALL_POINTS[1][partition_num][partition_num_3][0][0] = x;
					ALL_POINTS[1][partition_num][partition_num_3][0][1] = y;
					ALL_POINTS[1][partition_num][partition_num_3][0][2] = z;

					if(Projection_Ordering[Winning_Index] == 0)
					{
						Major_Axis_Plane_Projection->FastSet(x ,y , level, 1000);
					}
					if(Projection_Ordering[Winning_Index] == 1)
					{
						Major_Axis_Plane_Projection->FastSet(x,level, z, 1000);
					}
					if(Projection_Ordering[Winning_Index] == 2)
					{
						Major_Axis_Plane_Projection->FastSet(level, y, z, 1000);
					}
				
				}

				
				x = 0, y = 0, z = 0;
				intersection = 0;

				counter = top_size_2;
				// In this while loop we project from the 'bottom' side to 'top' sdie
				while(intersection == 0 && counter > bot_size_2)
				{
					x = (Projected_Normals[1][counter][0]) + (Normals[Winning_Index][part][0] - xcenter)  +  (Projected_Normals[0][part_3][0] - xcenter);
					y = (Projected_Normals[1][counter][1]) + (Normals[Winning_Index][part][1] - ycenter)  +  (Projected_Normals[0][part_3][1] - ycenter);
					z = (Projected_Normals[1][counter][2]) + (Normals[Winning_Index][part][2] - zcenter)  +  (Projected_Normals[0][part_3][2] - zcenter);
		
					counter--;
					intersection = UnSkinned_Image->FastGet(x,y,z);

				}

				if(intersection != 0)
				{
					Intersection_Points->FastSet(x,y,z,2000);

					
					ALL_POINTS[1][partition_num][partition_num_3][1][0] = x;
					ALL_POINTS[1][partition_num][partition_num_3][1][1] = y;
					ALL_POINTS[1][partition_num][partition_num_3][1][2] = z;


					if(Projection_Ordering[Winning_Index] == 0)
					{
						Major_Axis_Plane_Projection->FastSet(x ,y , level, 1000);
					}
					if(Projection_Ordering[Winning_Index] == 1)
					{
						Major_Axis_Plane_Projection->FastSet(x,level, z, 1000);
					}
					if(Projection_Ordering[Winning_Index] == 2)
					{
						Major_Axis_Plane_Projection->FastSet(level, y, z, 1000);
					}

				}
				

				partition_num_3++;
				JUMP_3 = JUMP_3 + top_bot_partition_3;
			}

			New_Highs_Lows.clear();


		partition_num++;
		JUMP = JUMP + top_bot_partition;

		level++;

	}
	

	std::cout << "Parameterization complete, calculating feature at Intersection Points..." << "\n";

	/*
	// We now erase any tangencies, e.g. where we collected only one point with an intersection, where MIN = MAX.
	for(int h = 0; h <= 1; h++) // Planes
			for(int j = 0; j < 30; j++)  // Grid
				for(int k = 0; k < 30; k++)  // Grid
					{
						if( ALL_POINTS[h][j][k][0][0] - ALL_POINTS[h][j][k][1][0] == 0 &&  ALL_POINTS[h][j][k][0][1] - ALL_POINTS[h][j][k][1][1] == 0 &&  ALL_POINTS[h][j][k][0][2] - ALL_POINTS[h][j][k][1][2] == 0)
						{
							ALL_POINTS[h][j][k][0][0] = 0;
							ALL_POINTS[h][j][k][1][0] = 0;
							ALL_POINTS[h][j][k][0][1] = 0;
							ALL_POINTS[h][j][k][1][1] = 0;
							ALL_POINTS[h][j][k][0][2] = 0;
							ALL_POINTS[h][j][k][1][2] = 0;
						}

					}

	// Unlick this portion to find and replace any points within a given radius of one another on each Intersection Plane.
	
	// Now here is where we eliminate close points.
	for(int h = 0; h <= 1; h++) // Planes
		for(int i = 0; i < 30; i++)  // Slice intersected with Plane 
			for(int j = 0; j < 30; j++)  // Points in each of those slices
				for(int k = 0; k <= 1; k++) // Min/Max points
				{
					for(int m = 0; m <= 1; m++)// Planes
						for(int n = 0; n < 30; n++)// Points in each of those slices
							for(int p = 0; p <= 1; p++)// Min/Max points
							{
								if(  (ALL_POINTS[h][i][j][k][0] != 0 && ALL_POINTS[h][i][j][k][1] != 0 && ALL_POINTS[h][i][j][k][2] != 0 ) &&  (ALL_POINTS[m][i][n][p][0] != 0 && ALL_POINTS[m][i][n][p][1] != 0 && ALL_POINTS[m][i][n][p][2] != 0 ) )
								{
									
									double dist = 0;
									if(mop == 1)
									{
										dist = (double)(ALL_POINTS[h][i][j][k][0] - ALL_POINTS[m][i][n][p][0])*(ALL_POINTS[h][i][j][k][0] - ALL_POINTS[m][i][n][p][0]) +  (double)(ALL_POINTS[h][i][j][k][1] - ALL_POINTS[m][i][n][p][1])*(ALL_POINTS[h][i][j][k][1] - ALL_POINTS[m][i][n][p][1]) +  (double)zmax_original_1/(double)(xmax)*(double)(ALL_POINTS[h][i][j][k][2] - ALL_POINTS[m][i][n][p][2])*(ALL_POINTS[h][i][j][k][2] - ALL_POINTS[m][i][n][p][2]);
									}

									if(mop == 2)
									{
										dist = (double)(ALL_POINTS[h][i][j][k][0] - ALL_POINTS[m][i][n][p][0])*(ALL_POINTS[h][i][j][k][0] - ALL_POINTS[m][i][n][p][0]) +  (double)(ALL_POINTS[h][i][j][k][1] - ALL_POINTS[m][i][n][p][1])*(ALL_POINTS[h][i][j][k][1] - ALL_POINTS[m][i][n][p][1]) +  (double)zmax_original_2/(double)(xmax)*(double)(ALL_POINTS[h][i][j][k][2] - ALL_POINTS[m][i][n][p][2])*(ALL_POINTS[h][i][j][k][2] - ALL_POINTS[m][i][n][p][2]);
									}



									if(dist < (double)(Final_1)/2.0)
									{
										ALL_POINTS[m][i][n][p][0] =	ALL_POINTS[h][i][j][k][0];
										ALL_POINTS[m][i][n][p][1] =	ALL_POINTS[h][i][j][k][1];
										ALL_POINTS[m][i][n][p][2] =	ALL_POINTS[h][i][j][k][2];
									}
								}
							}
				}
	*/
	

	for(int z = 0; z < zmax; z++)
		for(int y = 0; y < ymax; y++)
			for(int x = 0; x < xmax; x++)
			{
				Major_Axis_Plane_Projection->FastSet(x,y,z,0);
				Intersection_Points->FastSet(x,y,z,0);
			}

	for(int z = 0; z < xmax; z++)
		for(int x = 0; x < xmax; x++)
			for(int y = 0; y < ymax; y++)
			{
				Intersection_Points->FastSet(x,y,z,Skinned_Image->FastGet(x,y,z));
			}

		int nevel = 0;

		// These variables are used in projecting the Intersection Points onto consecutive primary planes, whichi are printed out.
		if(Projection_Ordering[Winning_Index] == 0)
		{
			nevel = zcenter;
		}
		if(Projection_Ordering[Winning_Index] == 1)
		{
			nevel = ycenter;
		}			
		if(Projection_Ordering[Winning_Index] == 2)
		{				
			nevel = xcenter;
		}

		if(mop == 1)
		{
			ALL_POINTS_1 = ALL_POINTS;
		}
		else // if mop == 2
		{
			ALL_POINTS_2 = ALL_POINTS;
		}
		 
		
		/*

		 v)    We calculate the Objects' Surface Feature at each Intersection Point.  The current feature, we call
		       the Manay Feature, is the volume of intersection of a sphereoid/parallelepiped (or any object really),
			   called the Seed here, intersected with the interior of Object.  
			   a)  We define the axes of Seed by d/m, the distance via Plane Normal between two neighboring Intersector
			       Planes.
			   b)  We normalize the volume of the Seed by dividing by the Seed's full volume.

		 
		 vi)   We compare the Manay Feature of each point on corresponding Connection Planes, and map the difference
			   back onto each Object's surface.

		*/
//Billy was here, figure out what this is, added int partial ball volume
		int partial_ball_volume;	
		for(int i = 0; i < 30; i++)  // Slice intersected with Plane 
			{
				for(int h = 0; h <= 1; h++) // Planes
					for(int j = 0; j < 30; j++)  // Points in each of those slices
						for(int k = 0; k <= 1; k++) // Min/Max points
						{
							if( ALL_POINTS[h][i][j][k][0] != 0 && ALL_POINTS[h][i][j][k][1] != 0 && ALL_POINTS[h][i][j][k][2] != 0 ) 
							{
								if(mop == 1)
								{
									//Billy took out int partial ball volume
									partial_ball_volume = Feature(ALL_POINTS[h][i][j][k][0],ALL_POINTS[h][i][j][k][1],ALL_POINTS[h][i][j][k][2], UnSkinned_Image);
									ALL_MEAS_1[h][i][j][k] = (double)(partial_ball_volume)/(double)(ball_volume);
									Intersection_Points->FastSet(ALL_POINTS[h][i][j][k][0],ALL_POINTS[h][i][j][k][1],ALL_POINTS[h][i][j][k][2],5000);
								
								}
								
								if(mop == 2)
								{
									//Billy
									partial_ball_volume = Feature(ALL_POINTS[h][i][j][k][0],ALL_POINTS[h][i][j][k][1],ALL_POINTS[h][i][j][k][2], UnSkinned_Image);
									ALL_MEAS_2[h][i][j][k] = (double)(partial_ball_volume)/(double)(ball_volume);
									Intersection_Points->FastSet(ALL_POINTS[h][i][j][k][0],ALL_POINTS[h][i][j][k][1],ALL_POINTS[h][i][j][k][2],5000);
								}

								// This can be commented out, is done above.
								if(Projection_Ordering[Winning_Index] == 0)
								{
									Major_Axis_Plane_Projection->FastSet(ALL_POINTS[h][i][j][k][0] ,ALL_POINTS[h][i][j][k][1] , nevel, 200);
								}
								if(Projection_Ordering[Winning_Index] == 1)
								{
									Major_Axis_Plane_Projection->FastSet(ALL_POINTS[h][i][j][k][0],nevel, ALL_POINTS[h][i][j][k][2], 2000);
								}
								if(Projection_Ordering[Winning_Index] == 2)
								{
									Major_Axis_Plane_Projection->FastSet(nevel, ALL_POINTS[h][i][j][k][1], ALL_POINTS[h][i][j][k][2], 2000);
								}
								

 								ALL_POINTS[h][i][j][k][0] = 0;
								ALL_POINTS[h][i][j][k][1] = 0;
								ALL_POINTS[h][i][j][k][2] = 0;
							}
						}
					nevel++;
			}




	

	for(int x = 0; x < xmax; x++)
		for(int y = 0; y < ymax; y++)
			for(int z = 0; z < zmax; z++)
			{
				if(UnSkinned_Image->FastGet(x,y,z) != 0)
				{
					UnSkinned_Image->FastSet(x,y,z,500);
				}
				if(Intersection_Planes->FastGet(x,y,z) == 0)
				{
					Intersection_Planes->FastSet(x,y,z,UnSkinned_Image->FastGet(x,y,z));
				}
			}



	// At last, we print out the set of Connection Planes, Intersection Points, and Projected points.
	if(mop == 1)
	{
		Write_Analyze_File((char *)(new std::string(Run_Name + (*(new std::string("_INTERSECTION_PLANES.hdr")))))->c_str(), *Intersection_Planes);
		Write_Analyze_File((char *)(new std::string(Run_Name + (*(new std::string("_INTERSECTION_POINTS.hdr")))))->c_str(), *Intersection_Points);
		Write_Analyze_File((char *)(new std::string(Run_Name + (*(new std::string("_PROJECTED_POINTS.hdr")))))->c_str(), *Major_Axis_Plane_Projection);
	}
	if(mop == 2)
	{
		Write_Analyze_File((char *)(new std::string(Run_Name_2 + (*(new std::string("_INTERSECTION_PLANES.hdr")))))->c_str(), *Intersection_Planes);
		Write_Analyze_File((char *)(new std::string(Run_Name_2 + (*(new std::string("_INTERSECTION_POINTS.hdr")))))->c_str(), *Intersection_Points);
		Write_Analyze_File((char *)(new std::string(Run_Name_2 + (*(new std::string("_PROJECTED_POINTS.hdr")))))->c_str(), *Major_Axis_Plane_Projection);
	}

	Normals.clear();
	Face_Planes.clear();
	Projection_Ordering.clear();

	// We reset all of the manipulated images for Object 2's parameterization.
	for(int x = 0; x < xmax; x++)
		for(int y = 0; y < ymax; y++)
			for(int z = 0; z < zmax; z++)
			{
				Intersection_Points->FastSet(x,y,z,0);
 				Intersection_Planes->FastSet(x,y,z,0);
				Major_Axis_Plane_Projection->FastSet(x,y,z,0);
			}

	mop++;
}


void Grid_Line::Skin_it(int skin_pick, int skin_2, int max, CIS_Array_Image3D_short* INimage , CIS_Array_Image3D_short* OutImage)
{
	/*

	Removes all material except "skin" of binary object.  This includes "tops", which are 
	 two-dimensional skin-pieces occuring whenever the object's outline in two dimensional (x-y)
	 (e.g. a two-dimensional closed curve) does not change continuously through the z-slices.

	 This searches through all pixels of non-zero intensity and keeps only those with less than six neighboring pixels
	 with non-zero intensity.

	 'INimage' is the input image with full Object, 'OutImage' is the output skinned image.

	 If, instead, you want a slice-by-slice 2D-skin (e.g. no tops) adjust where marked in the code.

	*/
	if(max == 1)
	{
		xmax = xmax_1;
		ymax = ymax_1;
		zmax = zmax_1;
	}
	else
	{
		xmax = xmax_2;
		ymax = ymax_2;
		zmax = zmax_2;
	}

	int j = 0, k = 0, l = 0, jwt = 0, kwt = 0, lwt = 0, jarea = 0, karea = 0, larea = 0;
	// The next int will be the surface area of O.
	double surface_area = 0;
	int sum = 0;
	int skin_num = 0;
	int up = 0, down = 0, forward = 0, backward = 0, left = 0, right = 0, temp_up = 0, temp_down = 0, temp_left = 0, temp_right = 0, temp_forward = 0, temp_backward = 0;
	double uparea = 0, downarea = 0, forwardarea = 0, backwardarea = 0, leftarea = 0, rightarea = 0;

	for (int xs = 0; xs < xmax; xs++)
		for (int ys = 0; ys < ymax; ys++)
			for (int zs = 0; zs < zmax; zs++)
				if (INimage->FastGet(xs,ys,zs) != 0)
				{
					if (xs == 0 || ys == 0 || zs == 0)
					{
						if (xs == 0)
						{
							j = INimage->FastGet(xs+1,ys,zs);
							jwt = 1;
							jarea = (1-j)*ypixel*zpixel;
						}
						else
						{
							j = INimage->FastGet(xs+1,ys,zs) + INimage->FastGet(xs-1,ys,zs);
								jarea = (2-j)*ypixel*zpixel;
						}
						if (ys == 0)
						{
							k = INimage->FastGet(xs,ys+1,zs);
							kwt = 1;
							karea = (1-k)*xpixel*zpixel;
						}
						else
						{
							k = INimage->FastGet(xs,ys+1,zs) + INimage->FastGet(xs,ys-1,zs);
								karea = (2-k)*xpixel*zpixel;
						}
					   if (zs == 0)
						{
							l = INimage->FastGet(xs,ys,zs+1);
								lwt = 1;
							larea = (1-k)*xpixel*ypixel;
						}
						else
						{  
							l = INimage->FastGet(xs,ys,zs+1)+ INimage->FastGet(xs,ys,zs-1);
							sum = j + k + l;
							larea = (2-k)*xpixel*ypixel;
						}
						if (sum < (6 - (jwt + kwt + lwt)))
						{
							OutImage->FastSet(xs,ys,zs,skin_pick);
						}
						else
						{
							OutImage->FastSet(xs,ys,zs,skin_2);
						}
					}
					else
					{
						// All of the chosen pixel's neighbors.
						up = INimage->FastGet(xs,ys,zs+1);
						down = INimage->FastGet(xs,ys,zs-1);
						forward = INimage->FastGet(xs,ys+1,zs);
						backward = INimage->FastGet(xs,ys-1,zs);
						right = INimage->FastGet(xs+1,ys,zs);
						left = INimage->FastGet(xs-1,ys,zs);

						if(skin_pick == 1)
						{
							temp_up = up;
							
							up = INimage->FastGet(xs,ys,zs+1);

							if(temp_up == 10)
							{
								temp_up = 1;
							}

							down = INimage->FastGet(xs,ys,zs-1);
							temp_down = down;
							if(temp_down == 10)
							{
								temp_down = 1;
							}

							forward = INimage->FastGet(xs,ys+1,zs);
							temp_forward = forward;
							if(temp_forward == 10)
							{
								temp_forward = 1;
							}

							backward = INimage->FastGet(xs,ys-1,zs);
							temp_backward = backward;
							if(temp_backward == 10)
							{
								temp_backward = 1;
							}

							right = INimage->FastGet(xs+1,ys,zs);
							temp_right = right;
							if(temp_right == 10)
							{
								temp_right = 1;
							}

							left = INimage->FastGet(xs-1,ys,zs);
							temp_left = left;
							if(temp_left == 10)
							{
								temp_left = 1;
							}
							sum = temp_up + temp_down + temp_forward + temp_backward + temp_left + temp_right;

							int mopey = INimage->FastGet(xs,ys,zs);
							uparea = (1 - temp_up)*zpixel*zpixel;
							downarea = (1 - temp_down)*zpixel*zpixel;
							forwardarea = (1 - temp_forward)*ypixel*ypixel;
							backwardarea = (1 - temp_backward)*ypixel*ypixel;
							rightarea = (1 - temp_right)*xpixel*xpixel;
							leftarea = (1 - temp_left)*xpixel*xpixel;
							skin_num = 6;
						}
						
						if(skin_pick == 2)
						{
							sum = left + right + forward + backward;
							skin_num = 4;
						}

						if (sum < skin_num) // sum < 6 for 3-d skinning, sum < 4 for 2d skinning
						{
							OutImage->FastSet(xs,ys,zs,1);
							//surface_area = surface_area + (uparea + downarea + leftarea + rightarea + forwardarea + backwardarea);
						}
						else
						{
							OutImage->FastSet(xs,ys,zs,0); // change '0' to '10' and undo tops above to include tops
						}
				}
			}
}


/*void main(int argc, char **argv)
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

	

	Input: 
	1.  The two images to compare
	2.  A user defined DOOO THSIISS

	



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


}*/
