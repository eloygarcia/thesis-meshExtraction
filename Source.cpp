#include "auxiliarHeaders.h"

// ========================== MAIN ==============================
int main( int argc, char* argv[])
{
	std::cout << std::endl;

	//======================================================
	std::cout << std::endl;
	std::cout << "Mammogram file name:  \t" << argv[1] << std::endl;
	std::cout << std::endl;
	std::cout << "MRI file name:  \t" << argv[3] << std::endl;
	std::cout << std::endl;

// ================== Read Data ========================
	// mammogram
	std::string inputMammogramfilename = argv[1];
	
	// mri
	std::string inputSegmentationfilename = argv[2];
	std::string inputBreastMaskfilename = argv[3];
		// New Segmentation.
	std::string newSegmentation_filename = argv[4];
	
	// output 
	std::string output_dir = argv[5];

// ================= Mammogram actions. =================
	mammogram * input_mammogram = new mammogram;
		input_mammogram->read_mammogram_with_metadata( inputMammogramfilename );

	std::cout << std::endl;
	std::cout << "Mammogram Information: " << std::endl;
		// ==================== Mammmogram information =====================
		std::string temp_side = input_mammogram->getImageLaterality(true);
		if ((temp_side!="R ") && (temp_side!="L ")) {
			std::cout << "Check side breast !" << std::endl; 
			return EXIT_FAILURE; 
		}
		
// =================== MRI actions ======================
	MRimage * input_mri = new MRimage;
	
	// 1. Necesitas la segmentación completa de Albert para la localización 
	// autómatica de las condiciones de Contorno (BC)
	input_mri->read_segmentation(inputSegmentationfilename ); 
	// 2. La máscara para la extración de la malla superficial.
	input_mri->read_breastmask( inputBreastMaskfilename );
	// 3. Segmentación probabilistica del tejido glandular !!
	input_mri->read_probabilistic_segmentation( newSegmentation_filename );
	
		input_mri->set_approxDegree( 2 );  // 1 or 2 !!
		input_mri->SetIsoVoxel( 5. );
		input_mri->SetSideOfinterest( temp_side.c_str() );

	vtkSmartPointer<vtkPolyData> total_mesh = vtkSmartPointer<vtkPolyData>::New();
		input_mri->getmodel( total_mesh );

	std::string ts_ply = output_dir + "\\tissues.ply";
	vtkSmartPointer<vtkPLYWriter> writeMeshPly = vtkSmartPointer<vtkPLYWriter>::New();
		writeMeshPly->SetInputData( total_mesh );
		writeMeshPly->SetFileName( ts_ply.c_str() );
		writeMeshPly->SetFileTypeToASCII();
		writeMeshPly->Update();
	
	std::cout << std::endl;
		std::cout << "MRI Information: " << std::endl;

/* Necesitamos que guarde la malla de Tetrahedros inicial para poder hallar la correlación inicial-final !! */
	std::string initialgridname = output_dir + "\\initialGrid.vtk";
	vtkSmartPointer<vtkUnstructuredGridWriter> initialGridWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
		initialGridWriter->SetInputData( input_mri->get_tetrahedral_mesh() );
		initialGridWriter->SetFileName( initialgridname.c_str() );
		initialGridWriter->Update();
		initialGridWriter->Write();

	return EXIT_SUCCESS;


}
