#ifndef __auxiliarFunctionsMeshes_cpp
#define __auxiliarFunctionsMeshes_cpp

#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPoints.h"

#include "vtkDoubleArray.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"

extern void getPoints_compressedBreast( vtkSmartPointer<vtkUnstructuredGrid> niftysim_mesh,
									std::vector<double> &f_points)
{
	vtkPoints* initial_points = niftysim_mesh->GetPoints();
	vtkDataArray* disp = niftysim_mesh->GetPointData()->GetArray(0);

	double pt[3] = {0.0, 0.0, 0.0};
	double d[3] = {0.0, 0.0, 0.0};
	double temp[3] = {0.0, 0.0, 0.0};

	for(int i=0; i<initial_points->GetNumberOfPoints(); i++)
	{
		initial_points->GetPoint(i, pt);
		disp->GetTuple(i, d);

		temp[0] = pt[0] + d[0];
		temp[1] = pt[1] + d[1];
		temp[2] = pt[2] + d[2];

		f_points[ 3*i ] = temp[0];
		f_points[ 3*i +1 ] = temp[1];
		f_points[ 3*i +2 ] = temp[2];
	}
}

extern void createUnstructuredGrid(std::vector<double> point_vector, 
							std::vector<int> cells_vector,
							vtkSmartPointer<vtkUnstructuredGrid> &grid)
{
	int numberofpoints = point_vector.size()/3;
	int numberofcells = cells_vector.size()/4;

	// ========================================================
	vtkPoints *points = vtkPoints::New();
	vtkCellArray* cells = vtkCellArray::New();

	points->SetNumberOfPoints( numberofpoints );
	double auxPoint[3] = {0.0, 0.0, 0.0};

	for(int pointId = 0; pointId < numberofpoints; pointId++)
	{
		auxPoint[0] = point_vector[ pointId*3 ];
		auxPoint[1] = point_vector[ pointId*3 +1];
		auxPoint[2] = point_vector[ pointId*3 +2];

		points->SetPoint(pointId, auxPoint);
	}

	vtkSmartPointer<vtkIdTypeArray> idCells = vtkSmartPointer<vtkIdTypeArray>::New();

	idCells->SetNumberOfComponents(5);
	idCells->SetNumberOfTuples( numberofcells );

	vtkIdType* tuple = new vtkIdType[4];

	for(int cellId=0; cellId < numberofcells; cellId++ )
	{
		tuple[0] = 4; 
		tuple[1] = cells_vector[ cellId*4 ];
		tuple[2] = cells_vector[ cellId*4 +1];
		tuple[3] = cells_vector[ cellId*4 +2];
		tuple[4] = cells_vector[ cellId*4 +3];

		idCells-> SetTupleValue( cellId, tuple);
	}

	cells->SetCells( numberofcells, idCells);

	grid->SetPoints(points);
	grid->SetCells(VTK_TETRA, cells);
}


/*
void get_NiftySim_deformedMesh(vtkSmartPointer<vtkUnstructuredGrid> initial_mesh,
							   vtkSmartPointer<vtkUnstructuredGrid> niftysim_mesh,
							   vtkSmartPointer<vtkUnstructuredGrid> &grid)
{
	// ===========================================================
	vtkPoints* initial_points = initial_mesh->GetPoints();
	// ===========================================================
	vtkDoubleArray* points_data = vtkDoubleArray::New();

	vtkCellArray *cells = niftysim_mesh->GetCells();
	vtkDataArray* disp = niftysim_mesh->GetPointData()->GetArray(0);
				
	vtkSmartPointer<vtkPoints> n_points = vtkSmartPointer<vtkPoints>::New();

	//n_mesh = mesh;
	vtkPoints* points = niftysim_mesh->GetPoints();
	
	points_data->SetNumberOfComponents(3);
	//points_data->SetNumberOfTuples( points->GetNumberOfPoints());
	
	// ================================================
	double pt_i[3] = {0.0, 0.0, 0.0};
	// ================================================
	double pt[3] = {0.0, 0.0, 0.0};
	double d[3] = {0.0, 0.0, 0.0};
	double temp[3] = {0.0, 0.0, 0.0};
	double temp_inv[3] = {0.0, 0.0, 0.0};
	
	for(int i=0; i<points->GetNumberOfPoints(); i++)
	{
		// ============================================
		initial_points->GetPoint(i, pt_i);
		// ============================================
		points->GetPoint(i, pt);
		disp->GetTuple(i, d); // Eliminar esta linea para hacer pruebas, reanudar cuando est'e niftysim
					//std::cout <<"[" << d[0] << ", "<< d[1] << ", "<< d[2] << "]" <<std::endl;
		temp[0] = pt[0] + d[0];
		temp[1] = pt[1] + d[1];
		temp[2] = pt[2] + d[2];

		temp_inv[0] = pt_i[0]-pt[0];
		temp_inv[1] = pt_i[1]-pt[1];
		temp_inv[2] = pt_i[2]-pt[2];

		n_points->InsertNextPoint(temp);
		points_data->InsertNextTuple(&temp_inv[0]);
	}
		
	grid->SetPoints(n_points);
	grid->SetCells( VTK_TETRA,cells);
	grid->GetPointData()->SetVectors(points_data);
	//grid->Update();

	double bound[6];
	n_points->GetBounds(bound);
//	std::cout << "Bound Nifty Sim mesh: [" <<bound[0] << ", " <<bound[1] << ", " <<bound[2] << ", " <<bound[3] << ", " <<bound[4] << ", " <<bound[5] << "] " <<std::endl;
}
*/
/*
void getNewMeshFromTwoSamples(vtkSmartPointer<vtkUnstructuredGrid> initial_mesh,
							  vtkSmartPointer<vtkUnstructuredGrid> final_mesh,
							  vtkSmartPointer<vtkUnstructuredGrid> &disp_grid)
{
	// ===========================================================
	vtkPoints* initial_points = initial_mesh->GetPoints();
	// ===========================================================
	vtkDoubleArray* points_data = vtkDoubleArray::New();
		points_data->SetNumberOfComponents(3);
	vtkCellArray *cells = final_mesh->GetCells();
	vtkPoints* points = final_mesh->GetPoints();
	// ===========================================================
	vtkSmartPointer<vtkPoints> n_points = vtkSmartPointer<vtkPoints>::New();
	// ================================================
		double pt_i[3] = {0.0, 0.0, 0.0};
		double pt[3] = {0.0, 0.0, 0.0};

		double temp_inv[3] = {0.0, 0.0, 0.0};
			
		for(int i=0; i<points->GetNumberOfPoints(); i++)	{
				// ============================================
				initial_points->GetPoint(i, pt_i);
				points->GetPoint(i, pt);
			
				temp_inv[0] = pt_i[0]-pt[0];
				temp_inv[1] = pt_i[1]-pt[1];
				temp_inv[2] = pt_i[2]-pt[2];

				n_points->InsertNextPoint(pt);
				points_data->InsertNextTuple(&temp_inv[0]);
		}
		
		disp_grid->SetPoints(n_points);
		disp_grid->SetCells( VTK_TETRA,cells);
		disp_grid->GetPointData()->SetVectors(points_data);
		//grid->Update();
}
*/

/*
// From mesh to image vector type
void newImage(ImageVector3Type::Pointer image, float* bounds, float* spacing)
{
	ImageVector3Type::PixelType base;
		base[0] =0;
		base[1] =0;
		base[2] =0;
	
	ImageVector3Type::IndexType start;
		start.Fill(0);

	ImageVector3Type::PointType origin;
		origin[0] = bounds[0];
		origin[1] = bounds[2];
		origin[2] = bounds[4];

	ImageVector3Type::SizeType size;
		size[0] = (bounds[1] - bounds[0]) / spacing[0];
		size[1] = (bounds[3] - bounds[2]) / spacing[1];
		size[2] = (bounds[5] - bounds[4]) / spacing[2];

//	ImageType::DirectionType direction;
//		direction[0][0] =1;  direction[0][1] = 0;   direction[0][2] =0;
//		direction[1][0] =0;  direction[1][1] = 0;   direction[1][2] =1;
//		direction[2][0] =0;  direction[2][1] = -1;   direction[2][2] =0;
		
	ImageVector3Type::RegionType region;
		region.SetSize( size );
		region.SetIndex( start );

	image->SetOrigin(origin);
	image->SetSpacing(spacing);
	//image->SetDirection( direction );
	image->SetRegions( region );
	image->Allocate();
	image->FillBuffer(base);
}
*/
/*
void get_deformationField(vtkSmartPointer<vtkUnstructuredGrid> tetrahedral_mesh,
						  ImageVector3Type::Pointer &deformationfield_image)
{
	vtkPoints* points = tetrahedral_mesh->GetPoints(); //->points pasa a ser n_points
		
	//vtkDataArray* disp = grid->GetPointData()->GetArray(0); ->disp pasa a ser points_data
		
		float bounds[6] =  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		//for(int i=0; i<6; i++) { bounds[i] = points->GetBounds()[i]; std::cout<<  n_points->GetBounds()[i] << std::endl; }
			bounds[0] = points->GetBounds()[0]-10;	bounds[1] = points->GetBounds()[1]+10; 
			bounds[2] = points->GetBounds()[2]-10;	bounds[3] = points->GetBounds()[3]+10; 
			bounds[4] = points->GetBounds()[4]-10;	bounds[5] = points->GetBounds()[5]+10; 
		
		float spacing[3] = {3, 3, 3};
		int size[3] = {0,0,0};

		size[0] =(int) ((bounds[1] - bounds[0]) / spacing[0]);
		size[1] =(int) ((bounds[3] - bounds[2]) / spacing[1]);
		size[2] =(int) ((bounds[5] - bounds[4]) / spacing[2]);

		vtkSmartPointer<vtkImageGridSource> image_grid = vtkSmartPointer<vtkImageGridSource>::New();
			image_grid->SetGridSpacing(spacing[0], spacing[1], spacing[2]);
			image_grid->SetGridOrigin(bounds[0], bounds[2], bounds[4]);
			
			image_grid->SetDataExtent(0, size[0]-1, 0, size[1]-1, 0, size[2]-1);
			image_grid->SetDataOrigin(bounds[0], bounds[2], bounds[4]);
			image_grid->SetDataSpacing(spacing[0], spacing[1], spacing[2]);

			std::cout << "IMAGE GRID................" << std::endl;

		vtkSmartPointer<vtkProbeFilter> probe = vtkSmartPointer<vtkProbeFilter>::New();
			probe->SetSourceData( tetrahedral_mesh );
			//probe->SetSource( grid  );
			probe->SetInputConnection ( image_grid->GetOutputPort() );
			probe->SpatialMatchOn();
		try{	
			probe->Update();
		} catch(int e) {
			std::cout << "Ya esta el probe jodiendo! " << e << std::endl;
		}

		//	std::cout << " probe update: " <<std::endl;
			
	// =============== Image - Deformation field ==================
		int count_tuple = 0;

		ImageVector3Type::Pointer image = ImageVector3Type::New(); 
		newImage(image, bounds, spacing);

		itk::Point<float,3> punt;
		double tuple_disp[3] = {0.0, 0.0, 0.0};
		ImageVector3Type::PixelType pixel;
		
		// Iterator
		itk::ImageRegionIterator<ImageVector3Type> it= itk::ImageRegionIterator<ImageVector3Type>(image, image->GetLargestPossibleRegion() );
		for(it.GoToBegin(); !it.IsAtEnd(); ++it)
		{
			image->TransformIndexToPhysicalPoint( it.GetIndex() , punt);
			//std::cout << "Pixel : [" <<punt[0] << ", "<<punt[1]<< ", " << punt[2] << "] ; y Punto : [" << probe->GetOutput()->GetPoint(count_tuple)[0] << ", " << probe->GetOutput()->GetPoint(count_tuple)[1] << ", " << probe->GetOutput()->GetPoint(count_tuple)[2] << "]."<<std::endl; 
			
			probe->GetOutput()->GetPointData()->GetArray(0)->GetTuple(count_tuple, tuple_disp);
			
			// if(probe->GetOutput()->GetPoint(count_tuple)[0]==punt[0] &&
			//	probe->GetOutput()->GetPoint(count_tuple)[1]==punt[1] &&
			//	probe->GetOutput()->GetPoint(count_tuple)[2]==punt[2] ){ 

				pixel[0]=tuple_disp[0];
				pixel[1]=tuple_disp[1];
				pixel[2]=tuple_disp[2];

				it.Set(pixel);
			// }else{ 
			//	pixel[0]=0;
			//	pixel[1]=0;
			//	pixel[2]=0;
			//	it.Set(pixel);
		//	}
			count_tuple+=1;
		}

	ImageVector3Type::DirectionType direction;
		direction[0][0] =1;  direction[0][1] = 0;   direction[0][2] =0;
		direction[1][0] =0;  direction[1][1] = 1;   direction[1][2] =0;
		direction[2][0] =0;  direction[2][1] = 0;   direction[2][2] =1;

		image->SetDirection(direction);

	ImageVector3Type::SpacingType spacingOutput;
		spacingOutput[0] = 0.5;
		spacingOutput[1] = 0.5;
		spacingOutput[2] = 0.5;

	ImageVector3Type::PointType originOutput = image->GetOrigin();

	ImageVector3Type::SizeType   sizeOutput;
		sizeOutput[0] = (int) ceil(size[0]*spacing[0]/spacingOutput[0]);
		sizeOutput[1] = (int) ceil(size[1]*spacing[1]/spacingOutput[1]);
		sizeOutput[2] = (int) ceil(size[2]*spacing[2]/spacingOutput[2]);

	LinearInterpolateVectorType::Pointer interpolator = LinearInterpolateVectorType::New();

	ResampleVectorType::Pointer resample = ResampleVectorType::New();
		resample->SetInput( image );
		resample->SetInterpolator( interpolator );

		resample->SetOutputSpacing( spacingOutput);
		resample->SetOutputOrigin( originOutput);
		resample->SetOutputDirection( direction );
		resample->SetSize(sizeOutput);

		resample->Update();
	
	deformationfield_image= resample->GetOutput();
}
*/
/*
void warpImage(Image3floatType::Pointer image, 
			   ImageVector3Type::Pointer deformation_field, 
			   Image3floatType::Pointer &warpimage)
{
	splineInterpolationMeshType::Pointer interpolator = splineInterpolationMeshType::New();
		interpolator->SetSplineOrder(3);
		
	//nearestInterpolationType::Pointer interpolator = nearestInterpolationType::New();

	warpMeshType::Pointer warpFilter = warpMeshType::New();
		warpFilter->SetInterpolator(interpolator);
		warpFilter->SetOutputOrigin( deformation_field->GetOrigin() );
		warpFilter->SetOutputSize( deformation_field->GetLargestPossibleRegion().GetSize() );
		warpFilter->SetOutputSpacing( deformation_field->GetSpacing() );
		warpFilter->SetOutputDirection( deformation_field->GetDirection() );

		warpFilter-> SetDisplacementField( deformation_field );
		warpFilter-> SetInput( image );

		warpFilter-> Update();

	warpimage = warpFilter->GetOutput();
}

*/


/* Feature-based registration */
/*
void getCompressedImage(std::vector<double> initial_points,
				   std::vector <double> deformed_points,
						std::vector<int> elements,
					Image3floatType::Pointer image_new_segmentation,  // esta no la necesito, .... creo !
					ImageVector3Type::Pointer &deformation_field,
					Image3floatType::Pointer &compressed_image,
					std::string output_dir)
{
	vtkSmartPointer<vtkUnstructuredGrid> mesh_grid_ini = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkUnstructuredGrid> mesh_grid_medium = vtkSmartPointer<vtkUnstructuredGrid>::New();
	//vtkSmartPointer<vtkUnstructuredGrid> mesh_grid_final = vtkSmartPointer<vtkUnstructuredGrid>::New();
	
	// new_mesh.set_translation( deformed_points, translation, final_points);

		createUnstructuredGrid( initial_points, elements, mesh_grid_ini); 
		createUnstructuredGrid( deformed_points, elements, mesh_grid_medium);
	//new_mesh.createUnstructuredGrid( final_points, elements, mesh_grid_final);

	vtkSmartPointer<vtkUnstructuredGrid> mesh_grid_ultimate = vtkSmartPointer<vtkUnstructuredGrid>::New();
	
		getNewMeshFromTwoSamples( mesh_grid_ini, mesh_grid_medium, mesh_grid_ultimate  );
	//	getNewMeshFromTwoSamples( mesh_grid_ini, mesh_grid_final, mesh_grid_ultimate );

	// ====================== Deformation Field =============================================
		get_deformationField(mesh_grid_ultimate, deformation_field);

	// 2. Write Deformation Field
	std::string defField_name = output_dir + "\\defField_image_2.mhd";

	WriterVector3Type::Pointer writer_deformationField = WriterVector3Type::New();
		writer_deformationField->SetInput( deformation_field  );
		writer_deformationField->SetFileName( defField_name );
		writer_deformationField->SetUseCompression( true );
		try{
			writer_deformationField->Update();
		} catch(itk::ExceptionObject & err) {
			std::cout << std::endl;
			std::cout << "Writer exception! " << std::endl;
			std::cout << err << std::endl;
			// return EXIT_FAILURE;
		}

	// Image3floatType::Pointer warpimage = Image3floatType::New();
	//	warpImage( image_new_segmentation, deformation_field, warpimage);
//		warpImage( image_new_segmentation, deformation_field, compressed_image);
	// compressed_image = warpimage;

//	std::cout << " ========================== WarpImageInformaci'on ================ " <<std::endl;
/*	Image3floatType::PointType warp_origin;
		warp_origin[0] = warpimage->GetOrigin()[0]; 
		warp_origin[1] = warpimage->GetOrigin()[1]; 
		warp_origin[2] = warpimage->GetOrigin()[2]; 
	
	Image3floatType::SizeType warp_size;
		warp_size[0] = warpimage->GetLargestPossibleRegion().GetSize()[0]; 
		warp_size[1] = warpimage->GetLargestPossibleRegion().GetSize()[1];
		warp_size[2] = warpimage->GetLargestPossibleRegion().GetSize()[2];

	Image3floatType::SpacingType warp_spacing;
		warp_spacing[0] = warpimage->GetSpacing()[0] ;
		warp_spacing[1] = warpimage->GetSpacing()[1] ;
		warp_spacing[2] = warpimage->GetSpacing()[2] ;

	Image3floatType::DirectionType warp_direction = warpimage->GetDirection();
*/
//}

#endif