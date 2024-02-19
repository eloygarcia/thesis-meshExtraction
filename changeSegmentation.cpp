#include "changeSegmentation.h"

changeSegmentation::changeSegmentation(void)
{
	//this->m_newSegmentation = ImageSegmentation3floatType::New();
}

changeSegmentation::~changeSegmentation(void)
{
}

// public
 void changeSegmentation::SetNewSegmentation( std::string newSegmentation_filename )
 {
	 ReaderSegmentationType::Pointer reader_new_segmentation = ReaderSegmentationType::New();
		reader_new_segmentation->SetFileName( newSegmentation_filename );
		try{
			reader_new_segmentation->Update();
		} catch (itk::ExceptionObject & excp ) {
			std::cout << " " << std::endl;
			std::cout << " New Segmentation Fails!" << std::endl;
			std::cout << excp << std::endl;
			std::cout << " " << std::endl;
			//return EXIT_FAILURE;
		}

	//ImageSegmentation3floatType::Pointer image_new_segmentation = reader_new_segmentation->GetOutput();

	ImageSegmentation3floatType::DirectionType direction_n;
		direction_n[0][0] =1;  direction_n[0][1] = 0;   direction_n[0][2] =0;
		direction_n[1][0] =0;  direction_n[1][1] = 1;   direction_n[1][2] =0;
		direction_n[2][0] =0;  direction_n[2][1] = 0;   direction_n[2][2] =1;

		//image_new_segmentation->SetDirection( direction_n );

	this->m_newSegmentation = reader_new_segmentation->GetOutput();
	this->m_newSegmentation->SetDirection( direction_n );
 }

 void changeSegmentation::SetMaskSegmentation( std::string maskSegmentation_filename )
 {
	 ReaderSegmentationType::Pointer reader_mask_segmentation = ReaderSegmentationType::New();
		reader_mask_segmentation->SetFileName( maskSegmentation_filename );
		try{
			reader_mask_segmentation->Update();
		} catch (itk::ExceptionObject & excp ) {
			std::cout << " " << std::endl;
			std::cout << " Mask Segmentation Fails!" << std::endl;
			std::cout << excp << std::endl;
			std::cout << " " << std::endl;
			//return EXIT_FAILURE;
		}

	//ImageSegmentation3floatType::Pointer image_new_segmentation = reader_new_segmentation->GetOutput();

	ImageSegmentation3floatType::DirectionType direction_n;
		direction_n[0][0] =1;  direction_n[0][1] = 0;   direction_n[0][2] =0;
		direction_n[1][0] =0;  direction_n[1][1] = 1;   direction_n[1][2] =0;
		direction_n[2][0] =0;  direction_n[2][1] = 0;   direction_n[2][2] =1;

		//image_new_segmentation->SetDirection( direction_n );

	this->m_maskSegmentation = reader_mask_segmentation->GetOutput();
	this->m_maskSegmentation->SetDirection( direction_n );
 }

 // private
 void changeSegmentation::getROI(ImageSegmentation3floatType::Pointer imagen, const char* side, ImageSegmentation3floatType::Pointer & temp_imagen)
 {
	  // Eliminar la parte central
	ImageSegmentation3floatType::SpacingType spacing =	imagen->GetSpacing();
	ImageSegmentation3floatType::SizeType size = imagen->GetLargestPossibleRegion().GetSize();

	int x_offset = 17 / spacing[0];
	int x_slice = size[0]/2-1;

	ImageSegmentation3floatType::IndexType idx;
	itk::ImageRegionIterator<ImageSegmentation3floatType> it = 
		itk::ImageRegionIterator<ImageSegmentation3floatType>(imagen, imagen->GetLargestPossibleRegion() );
	for(it.GoToBegin(); !it.IsAtEnd(); ++it) {
		idx = it.GetIndex();
		if(idx[0]>=(x_slice - x_offset) && idx[0]<=(x_slice+x_offset)) it.Set(0);
	}
	
	// Definir el lado del pecho que queremos.
	ImageSegmentation3floatType::IndexType startPoint;

		if(strcmp(side, "R ")==0){
			//std::cout << " Side R " << std::endl;// Pecho derecho: 0;  Rigth Breast: startPoint[0]=0;
			startPoint[0] = 0;
			startPoint[1] = 0;
			startPoint[2] = 0; }
		else if (strcmp(side,"L ")==0) {// Pecho Izquierdo: size[0]/2-1; Left Breast: startPoint[0]=size[0]/2;
			//std::cout << " Side L " << std::endl;
			startPoint[0]= size[0]/2 -1;
			startPoint[1] = 0;
			startPoint[2] = 0;}
		else { 
			std::cout << "Check side breast !" << std::endl; 
			//return EXIT_FAILURE; 
		}

	ImageSegmentation3floatType::SizeType sizeRoi; 
		sizeRoi[0]=size[0]/2 -1; 
		sizeRoi[1]=size[1];
		sizeRoi[2]=size[2];
	ImageSegmentation3floatType::RegionType desireRoi;
		desireRoi.SetIndex(startPoint);
		desireRoi.SetSize(sizeRoi);
		
	RoiSegmentationType::Pointer Roi = RoiSegmentationType::New();
		Roi-> SetRegionOfInterest(desireRoi);
		Roi-> SetInput( imagen );
		Roi-> Update();

	temp_imagen = Roi->GetOutput();
 }

 void changeSegmentation::ExtractGlandularTissue()
 {
	 std::cout << "Comienza la extracción de tejido Denso !! " << std::endl;
	 
	 ImageSegmentation3floatType::SpacingType spaz = m_maskSegmentation->GetSpacing();

	 // Dilating.........!
	KernelType ball1;
	KernelType::SizeType radious;
		radious[0] = ceil(10/spaz[0]);
		radious[1] = ceil(10/spaz[1]);
		radious[2] = ceil(10/spaz[2]);

		ball1.SetRadius( radious );
		ball1.CreateStructuringElement();
 
	DilateSegmentationType::Pointer dilating = DilateSegmentationType::New();
		dilating->SetInput( this->m_maskSegmentation );
		dilating->SetKernel( ball1 );
		dilating->SetForegroundValue( 0 );
		dilating->Update();
	
	 this->m_newSegmentation->SetOrigin( this->m_maskSegmentation->GetOrigin() );
	 
	 MaskSegmentationType::Pointer mask = MaskSegmentationType::New();
		mask->SetInput( this->m_newSegmentation );
		//mask->SetMaskImage( this->m_maskSegmentation );
		mask->SetMaskImage( dilating->GetOutput() );
		mask->SetCoordinateTolerance( 0.0001 );
		mask->SetMaskingValue( 0 );
//		mask->SetOutsideValue( 0 );
		try{
			mask->Update();
		} catch(itk::ExceptionObject & excp) {
			std::cout << "Masking dense tissue Exception" << std::endl;
			std::cout << excp << std::endl;
		}

	this->m_newSegmentation = mask->GetOutput();

	ImageSegmentation3floatType::Pointer newSeg =ImageSegmentation3floatType::New();
	 getROI(mask->GetOutput(), this->side, newSeg);
	 //getRegionOfInterest(this->m_maskSegmentation, this->side, this->m_maskSegmentation);

	MedianFilterType::Pointer median = MedianFilterType::New();
	MedianFilterType::InputSizeType radius;
		radius.Fill(2);

		median->SetInput( newSeg );
		median->SetRadius(radius);
		
	 ThresholdSegmentationType::Pointer thresholding = ThresholdSegmentationType::New();
		thresholding->SetInput( median->GetOutput() );
		thresholding->SetLowerThreshold( 0.5 );
		thresholding->SetUpperThreshold( 1. );
		thresholding->SetInsideValue( 1 );
		thresholding->SetOutsideValue( 0);
		thresholding->Update();

	// Apply surface mask

	// ? Resampling.
	
	ImageSegmentation3floatType::SizeType size = thresholding->GetOutput()->GetLargestPossibleRegion().GetSize();
	ImageSegmentation3floatType::SpacingType spacing = thresholding->GetOutput()->GetSpacing();
	ImageSegmentation3floatType::PointType originRoi = thresholding->GetOutput()->GetOrigin();
	ImageSegmentation3floatType::DirectionType directionRoi = thresholding->GetOutput()->GetDirection();

	ImageSegmentation3floatType::SpacingType spacingOutput;
		spacingOutput[0] = this->iso_vox_size;
		spacingOutput[1] = this->iso_vox_size;
		spacingOutput[2] = this->iso_vox_size;
	ImageSegmentation3floatType::SizeType sizeOutput;
		sizeOutput[0] = size[0]*spacing[0]/spacingOutput[0] ;
		sizeOutput[1] = size[1]*spacing[1]/spacingOutput[1] ;
		sizeOutput[2] = size[2]*spacing[2]/spacingOutput[2] ;

	AffTransType::Pointer transf = AffTransType::New();
	InterpolatorSegmentationType::Pointer interpolatorSeg = InterpolatorSegmentationType::New();
	ResampleSegmentationType::Pointer resampSeg = ResampleSegmentationType::New();
		resampSeg-> SetInterpolator(interpolatorSeg);
		resampSeg-> SetTransform(transf);

		resampSeg-> SetDefaultPixelValue(0);
		resampSeg-> SetOutputOrigin(originRoi);
		resampSeg-> SetOutputDirection(directionRoi);

		resampSeg-> SetOutputSpacing(spacingOutput);
		resampSeg-> SetSize(sizeOutput);

		resampSeg-> SetInput( thresholding->GetOutput() );
		resampSeg-> Update();
	
		
	// Dilating.........!
	FlatType::RadiusType radius1;
		radius1.Fill(1);
	FlatType box1 = FlatType::Box(radius1);

	FlatType::RadiusType radius2;
		radius2.Fill(2);
	FlatType box2 = FlatType::Box(radius2);
/* 
	KernelType ball2;
		ball2.SetRadius(2);
		ball2.CreateStructuringElement();
/*
	ErodeType::Pointer erode2 = ErodeType::New();
		erode2->SetInput( resampSeg->GetOutput() );
		erode2->SetKernel( ball2 );
		erode2->SetForegroundValue( 1 );
		erode2->Update();
/*
/*
	DilateSegmentationType::Pointer dilating2 = DilateSegmentationType::New();
		dilating2->SetInput( resampSeg->GetOutput() );
		dilating2->SetKernel( ball2 );
		dilating2->SetForegroundValue( 1 );
		dilating2->Update();
/*
/*
	Closing3floatType::Pointer closing = Closing3floatType::New();
		closing->SetInput( resampSeg->GetOutput() );
		closing->SetKernel( ball2 );
		closing->SetForegroundValue( 1 );
		closing->Update();
		/**/

/*
	ErodeFlatType::Pointer erode2 = ErodeFlatType::New();
		erode2->SetInput( resampSeg->GetOutput() );
		erode2->SetKernel( box1 );
		erode2->SetForegroundValue( 1 );
		erode2->Update();
*/
	DilateFlatType::Pointer dilating2 = DilateFlatType::New();
		dilating2->SetInput( resampSeg->GetOutput() );
		dilating2->SetKernel( box1 );
		dilating2->SetForegroundValue( 1 );
		dilating2->Update();
/*
	ClosingFlatType::Pointer closing = ClosingFlatType::New();
		closing->SetInput( resampSeg->GetOutput() );
		closing->SetKernel( box );
		closing->SetForegroundValue( 1 );
		closing->Update();
*/
	// From itk to vtk
	ConnectTissueType::Pointer connecting = ConnectTissueType::New();
		//connecting->SetInput( this->m_maskSegmentation );
		//connecting->SetInput( closing->GetOutput() );
		 connecting->SetInput( dilating2->GetOutput() );
		// connecting->SetInput( erode2->GetOutput() );
		//connecting->SetInput( resampSeg->GetOutput() );
		try {
			connecting-> Update();
		} catch (itk::ExceptionObject & err)  {
			std::cerr << " " << std::endl;
			std::cerr << "Conection Fail! " << std::endl; 
			std::cerr << err << std::endl;
		}

	// vtk.
	vtkSmartPointer<vtkMarchingCubes> marchingCubes = vtkSmartPointer<vtkMarchingCubes>::New();
		marchingCubes-> SetInputData(connecting-> GetOutput());
		marchingCubes-> ComputeNormalsOn();
		marchingCubes-> ComputeScalarsOff();
		marchingCubes-> SetValue(0,1);
		marchingCubes-> Update();

	vtkSmartPointer<vtkSmoothPolyDataFilter> smoothy = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
		smoothy->SetRelaxationFactor(0.02);
		smoothy->SetNumberOfIterations(100);
		smoothy->SetInputData( marchingCubes->GetOutput() );
		smoothy->Update();

	vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
		connectivityFilter->SetInputConnection( smoothy->GetOutputPort() );
		connectivityFilter->SetExtractionModeToLargestRegion(); 
		connectivityFilter->Update();

	//this->m_GlandularMesh = smoothy->GetOutput();
	this->m_GlandularMesh = connectivityFilter->GetOutput();

//	std::cout << "Comienza la visualización del tejido denso!" << std::endl;
/* Visualización */
/* 
	double a[3] = {1.0, 0.0, 0.0};
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		//mapper->SetInputData( smoothy->GetOutput() );
		mapper->SetInputData( connectivityFilter->GetOutput() );

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper( mapper );
		actor->GetProperty()->GetColor( a );
		//actor->GetProperty()->SetRepresentationToWireframe();

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
		renderer->AddActor( actor );
		renderer->ResetCamera();
		renderer->SetBackground( 1.0,  1.0, 1.0 );


	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
		renderWindow->AddRenderer( renderer );

	vtkSmartPointer<vtkRenderWindowInteractor> Riterator = vtkSmartPointer<vtkRenderWindowInteractor>::New();
		Riterator->SetRenderWindow( renderWindow );

	renderWindow->Render();
	Riterator->Initialize();
	Riterator->Start();
/**/
 }

 void changeSegmentation::approxTissue()
 {
	 std::cout << "Comeanzando a aproximar tejido" << std::endl;

	 vtkPoints* allPoints = this->m_GlandularMesh->GetPoints();

	double temp_body[3]= {0.0, 0.0, 0.0};
	double temp_breast[3]= {0.0, 0.0, 0.0};
	
	std::vector<double> points;
	std::vector<int> id;


	for(vtkIdType id_breast = 0; id_breast<allPoints->GetNumberOfPoints(); id_breast++){
		allPoints->GetPoint(id_breast, temp_breast );

		//points_surface.push_back(1);
		points.push_back(temp_breast[0]);
		points.push_back(temp_breast[1]);
		points.push_back(temp_breast[2]);
		
		id.push_back( id_breast+1 );
		// id.push_back( id_breast );
	}
	
	// // //

	const vtkIdType numCells = this->m_GlandularMesh->GetNumberOfCells();
	std::vector<int> connect_surface;
	
	for(vtkIdType cellId=0 ; cellId<numCells ; ++cellId)
	{		
		// Get all cells and its points.
		vtkCell* cell = this->m_GlandularMesh-> GetCell(cellId);

		if( cell!=NULL ) { 
			vtkIdList* list;
			list = cell->GetPointIds();
			connect_surface.push_back( list->GetId( 0 ) );
			connect_surface.push_back( list->GetId( 1 ) );
			connect_surface.push_back( list->GetId( 2 ) );
			//std::cout<< list->GetId(0)<< ", " << list->GetId(1)<< ", " << list->GetId(2)<< std::endl;
		}
	}

 // // // 
	
	double control;

	 // // // 
	double pt = 0.0;
	float temp_point[3] = {0.0,0.0,0.0};

	if(this->m_degree==1){
		for(int l =0; l<(points.size()/3); l){ 
			control = (this->m_solution[0]-1) +
				points[ 3*l ] * this->m_solution[1] + points[ 3*l+1 ] * this->m_solution[2];
			
			if((points[ 3*l +2 ] - control) > 0){
					pt = (this->m_solution[0]-1) + 
						points[3*l] * this->m_solution[1] + points[3*l+1] * this->m_solution[2];
					//points_surface[4*last+1] = points_interface[4*ll+1];
					//points_surface[4*last+2] = points_interface[4*ll+2];
					points[3*l+2] = pt;
			}
		}
	} else if(this->m_degree==2){
		for(int l =0; l<(points.size()/3); l++){ 
			control = (this->m_solution[0]-1) +
				points[ 3*l ]*this->m_solution[1] +	points[ 3*l+1 ]*this->m_solution[2] + pow(points[3*l],2)*this->m_solution[3] +
					pow(points[3*l+1],2)*this->m_solution[4] + (points[3*l]*points[3*l+1])*this->m_solution[5];
			
			if((points[ 3*l+2 ] - control) > 0){
				pt = (this->m_solution[0]-1) +
					points[3*l] * this->m_solution[1] + points[3*l+1] * this->m_solution[2] + pow(points[3*l],2) * this->m_solution[3] +
					pow(points[3*l+1],2) * this->m_solution[4] + (points[3*l]*points[3*l+1]) * this->m_solution[5];
					//points_surface[4*last+1] = points_interface[4*ll+1];
					//points_surface[4*last+2] = points_interface[4*ll+2];
					points[3*l+2] = pt;
				}
			}
		}

	std::list<int> borderlisto;
	cleaningMesh cleo;
	// AQUI HEMOS CAMBIADO EL CODIGO
	bool repeat = true;
	//while(repeat) {
/*	std::cout << "comienza limpieza de la malla" << std::endl;
		cleo.removeDuplicatePoints(points, connect_surface);
		 cleo.removeDuplicatedCells(connect_surface);
	
		cleo.foundNonManifold(connect_surface, borderlisto);
		cleo.clearMesh(connect_surface, borderlisto, repeat);
		cleo.performNewMesh(points, connect_surface);
		borderlisto.clear();
*/	//}

	/// 
	std::cout << "comienza construcción de la malla" << std::endl;
	vtkSmartPointer<vtkPoints> mesh_points = vtkSmartPointer<vtkPoints>::New();

	for(int last_count = 0; last_count<points.size()/3; ++last_count){ // recordar cambiar
		temp_point[0]=points[3*last_count];
		temp_point[1]=points[3*last_count+1];
		temp_point[2]=points[3*last_count+2];

		mesh_points ->InsertNextPoint ( temp_point );
	}

	vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
		for( int rat =0 ; rat <(connect_surface.size())/3 ; rat++){
			polys -> InsertNextCell( 3 );
			polys -> InsertCellPoint( connect_surface[3*rat] );
			polys -> InsertCellPoint( connect_surface[3*rat+1] );
			polys -> InsertCellPoint( connect_surface[3*rat+2] );
		}

	vtkSmartPointer<vtkPolyData> polygonPolyData = vtkSmartPointer<vtkPolyData>::New();
		polygonPolyData->SetPoints(mesh_points);
		polygonPolyData->SetPolys(polys);
		polygonPolyData->BuildCells();  // PolyData Mesh !!

	this->m_GlandularMesh->Reset();
	this->m_GlandularMesh = polygonPolyData;

/* Visualización */
/*
	double a[3] = {1.0, 0.0, 0.0};
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputData( m_GlandularMesh );

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper( mapper );
		actor->GetProperty()->GetColor( a );
		//actor->GetProperty()->SetRepresentationToWireframe();

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
		renderer->AddActor( actor );
		renderer->ResetCamera();
		renderer->SetBackground( 1.0,  1.0, 1.0 );


	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
		renderWindow->AddRenderer( renderer );

	vtkSmartPointer<vtkRenderWindowInteractor> Riterator = vtkSmartPointer<vtkRenderWindowInteractor>::New();
		Riterator->SetRenderWindow( renderWindow );

	renderWindow->Render();
	Riterator->Initialize();
	Riterator->Start();
/**/

 }
 
 // public
 void changeSegmentation::JoinMeshes(vtkSmartPointer<vtkPolyData> surf_mesh, vtkSmartPointer<vtkPolyData> &total_mesh)
 {
	 std::cout << "Comenzando a unir mallas " << std::endl;
	 vtkSmartPointer<vtkAppendPolyData> appfake = vtkSmartPointer<vtkAppendPolyData>::New();
		appfake->SetInputData( surf_mesh );
		appfake->Update();
		
	vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
		connectivityFilter->SetInputConnection( appfake->GetOutputPort() );
		connectivityFilter->SetExtractionModeToLargestRegion(); 
		connectivityFilter->Update();

	 ExtractGlandularTissue();
	 approxTissue();

	 vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
		appendFilter->AddInputData( this->m_GlandularMesh);
		appendFilter->AddInputData( connectivityFilter->GetOutput() );
		appendFilter->Update();

/*	vtkSmartPointer<vtkCleanPolyData> clean = vtkSmartPointer<vtkCleanPolyData>::New();
		clean->SetInputData( appendFilter->GetOutput() );
		clean->ConvertStripsToPolysOn();
		clean->ConvertPolysToLinesOn();
		clean->ConvertLinesToPointsOn();
		clean->PointMergingOn();
		//clean->SetTolerance(0.01);
		clean->Update();
		
	total_mesh = clean->GetOutput();
	*/
		total_mesh = appendFilter->GetOutput();
/*
	std::string ts_ply = "C:\\Users\\Eloy\\Desktop\\tissues_20_10.ply";
	vtkSmartPointer<vtkPLYWriter> writeMeshPly = vtkSmartPointer<vtkPLYWriter>::New();
		writeMeshPly->SetInputData( total_mesh );
		writeMeshPly->SetFileName( ts_ply.c_str() );
		writeMeshPly->SetFileTypeToASCII();
		writeMeshPly->Update();


		/* Visualización */
/*
	double a[3] = {1.0, 0.0, 0.0};
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputData( appendFilter->GetOutput() );

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper( mapper );
		actor->GetProperty()->GetColor( a );
		actor->GetProperty()->SetRepresentationToWireframe();

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
		renderer->AddActor( actor );
		renderer->ResetCamera();
		renderer->SetBackground( 1.0,  1.0, 1.0 );

	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
		renderWindow->AddRenderer( renderer );

	vtkSmartPointer<vtkRenderWindowInteractor> Riterator = vtkSmartPointer<vtkRenderWindowInteractor>::New();
		Riterator->SetRenderWindow( renderWindow );

	renderWindow->Render();
	Riterator->Initialize();
	Riterator->Start();
/**/
 }