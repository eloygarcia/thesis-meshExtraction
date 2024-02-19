#ifndef __MRimage_cpp
#define __MRimage_cpp

#include "MRimage.h"

MRimage::MRimage(void)
{
	
	this->direction[0][0] = 1;	this->direction[0][1] = 0;	this->direction[0][2] = 0;
	this->direction[0][1] = 0;	this->direction[1][1] = 1;	this->direction[1][2] = 0;
	this->direction[0][2] = 0;	this->direction[1][2] = 0;	this->direction[2][2] = 1;
	// breast mask
	this->m_breastmask = ImageMriType::New();
	this->has_breastmask = false;
	
	// albert segmentation
	this->m_albertsegmentation = ImageMriType::New();
	this->has_segmentation = false;

	// mesh
	this->m_model = vtkSmartPointer<vtkUnstructuredGrid>::New();

	// New Segmentation !
	new_segmentation = new changeSegmentation;
	// New Mesh !
	new_mesh = new mesh;

	istheremodel = false;
}

MRimage::~MRimage(void)
{
	
}


void MRimage::closeBrestMask(  ImageMriType::Pointer mask, ImageMriType::Pointer &closingMask )
{
	ImageMriType::SpacingType spacing= mask->GetSpacing();
	ImageMriType::SizeType size = mask->GetLargestPossibleRegion().GetSize();

	int thres = ceil( 10/spacing[0] ); // 5 mm en X
	
	std::cout << std::endl;
	std::cout << "El size es: [" << size[0] << ", " << size[1] << ", " << size[2] << "]"<< std::endl;
	std::cout << "El thres es : " << thres << std::endl;
	std::cout << std::endl;

	LabelType::Pointer stats = LabelType::New();
		stats->SetLabelInput(mask);
		stats->SetInput(mask);
		stats->Update();

	itk::ImageRegionIterator<ImageMriType> it = itk::ImageRegionIterator<ImageMriType>(mask, mask->GetLargestPossibleRegion() );
	for(it.Begin(); !it.IsAtEnd(); it++) {
		ImageMriType::IndexType idx = it.GetIndex();
		if( idx[0]< thres || idx[0]> size[0]-thres ) it.Set( 0 );
	}

	closingMask = mask;
}

// Read data.
// --------------------------------------------------------
void MRimage::read_breastmask(std::string inputfilename)
{
	
	ReaderMriType::Pointer reader_mask= ReaderMriType::New();
		reader_mask->SetFileName( inputfilename);
		try{
			reader_mask->Update();
		}catch(itk::ExceptionObject & err){
			std::cout <<" " << std::endl;
			std::cout <<"Image File Reader Fails! " << std::endl;
			std::cout << err << std::endl;
			//return EXIT_FAILURE;
		}

	this->breastmaskfilename = inputfilename ;
	// this->m_breastmask = reader_mask->GetOutput();
		closeBrestMask(  reader_mask->GetOutput(), this->m_breastmask);
	this->m_breastmask->SetDirection( this->direction );
	this->has_breastmask  = true;
}
// --------------------------------------------------------
void MRimage::read_segmentation(std::string inputfilename)
{
	ReaderMriType::Pointer reader_seg= ReaderMriType::New();
		reader_seg->SetFileName( inputfilename);
		try{
			reader_seg->Update();
		}catch(itk::ExceptionObject & err){
			std::cout <<" " << std::endl;
			std::cout <<"Image File Reader Fails! " << std::endl;
			std::cout << err << std::endl;
			//return EXIT_FAILURE;
		}
		
	this->albertsegmentationfilename = inputfilename;
	this->m_albertsegmentation = reader_seg->GetOutput();
	this->m_albertsegmentation->SetDirection( this->direction );
	this->has_segmentation  = true;
}
// --------------------------------------------------------
void MRimage::read_probabilistic_segmentation( std::string inputfilename )
{
	new_segmentation->SetNewSegmentation( inputfilename );
	new_segmentation->SetMaskSegmentation( this->breastmaskfilename );
}

// aux. Functions
// -------------------------------------------------------------
void MRimage::getRegionOfInterest(ImageMriType::Pointer image, const char* side, ImageMriType::Pointer &new_image)
{
	// From Albert Gubern paper
	ImageMriType::SpacingType spacing = image->GetSpacing();
	ImageMriType::SizeType size = image->GetLargestPossibleRegion().GetSize();
	//std::cout << size[0] << "," << size[1] << "," <<size[2] << std::endl;
	// 17 mm, previous 25
	int x_offset = 17 / spacing[0]; //pixels
	
	// 7 mm, previous 10
	//int z_offset = 7 / spacing[2]; //pixels
	//int st = sternumPoint + 20/spacing[2];
	
    //int z_slice = size[2]/2-1;
	int x_slice = size[0]/2-1;

	ImageMriType::IndexType idx;
	itk::ImageRegionIterator<ImageMriType> it = itk::ImageRegionIterator<ImageMriType>(image, image->GetLargestPossibleRegion());
	for(it.GoToBegin(); !it.IsAtEnd(); ++it){
		idx = it.GetIndex();
		//std::cout << idx << std::endl;
		if(idx[0]>=(x_slice - x_offset) && idx[0]<=(x_slice + x_offset)){
			it.Set(0);	
		}
	}
	
	// =============== Region of interest ======================
	ImageMriType::IndexType startPoint;

		if (strcmp(side,"R ")==0) {
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
		
	ImageMriType::SizeType sizeRoi; 
		sizeRoi[0]=size[0]/2 -1; 
		sizeRoi[1]=size[1];
		sizeRoi[2]=size[2];
	ImageMriType::RegionType desireRoi;
		desireRoi.SetIndex(startPoint);
		desireRoi.SetSize(sizeRoi);
		
	RoiMriType::Pointer Roi = RoiMriType::New();
		Roi-> SetRegionOfInterest(desireRoi);
		Roi-> SetInput( image );
		Roi-> Update();
	
	//	std::cout << " " << std::endl;
	//	std::cout << "Roi ok!" << std::endl;
		//std::cout << "Start Point: [" << startPoint[0] << ", " << startPoint[1] << ", " << startPoint[2] << "]" << std::endl;
	//ImageType::Pointer segmentation = Roi_segmentation->GetOutput();

	new_image = Roi->GetOutput();
}

// -------------------------------------------------------------
void MRimage::resampleFunction (ImageMriType::Pointer image, float iso_vox_size, ImageMriType::Pointer &outputimage)
{
	ImageMriType::SizeType size = image->GetLargestPossibleRegion().GetSize();
	ImageMriType::SpacingType spacing = image->GetSpacing();
	ImageMriType::PointType origin = image->GetOrigin();
	ImageMriType::DirectionType direction = image->GetDirection();
	
	ImageMriType::SpacingType spacingOutput;
		spacingOutput[0] = iso_vox_size;
		spacingOutput[1] = iso_vox_size;
		spacingOutput[2] = iso_vox_size;

	ImageMriType::SizeType sizeOutput;
		sizeOutput[0] = size[0]*spacing[0]/spacingOutput[0] ;
		sizeOutput[1] = size[1]*spacing[1]/spacingOutput[1] ;
		sizeOutput[2] = size[2]*spacing[2]/spacingOutput[2] ;

	InterpolatorMriType::Pointer interpolator = InterpolatorMriType::New();
	
	AffineTransformType::Pointer transform = AffineTransformType::New();		

	ResampleMriType::Pointer resample = ResampleMriType::New();
		resample-> SetInterpolator(interpolator);
		resample-> SetTransform(transform);

		resample-> SetDefaultPixelValue(0);
		resample-> SetOutputOrigin(origin);
		resample-> SetOutputDirection(direction);

		resample-> SetOutputSpacing(spacingOutput);
		resample-> SetSize(sizeOutput);

		resample-> SetInput(image);
		resample-> Update();

	outputimage = resample->GetOutput();
}

// -------------------------------------------------------------
void MRimage::closingFunction (ImageMriType::Pointer image, BallType ball, ImageMriType::Pointer &outputImage, int r)
{
	ClosingMriType::Pointer closing = ClosingMriType::New();
		closing->SetInput( image );
		closing->SetKernel( ball );
		closing->SetForegroundValue(r);
		closing->Update();

	outputImage = closing->GetOutput();
}

// -------------------------------------------------------------
void MRimage::visualizationFunction (ImageMriType::Pointer image, vtkSmartPointer<vtkMarchingCubes> & marchingOutput) //,  vtkSmartPointer<vtkActor> &actorOutput)
{	
	ConnectorType::Pointer connector = ConnectorType::New();
		connector-> SetInput(image );
		try {
			connector-> Update();  }
		catch (itk::ExceptionObject & err)  {
			std::cerr << " " << std::endl;
			std::cerr << "Conection Fail! " << std::endl; 
			std::cerr << err << std::endl;
		//	return EXIT_FAILURE;
		}

	/* VTK */
	vtkSmartPointer<vtkMarchingCubes> marchingCubes = vtkSmartPointer<vtkMarchingCubes>::New();
		marchingCubes-> SetInputData(connector-> GetOutput());
		marchingCubes-> ComputeNormalsOn();
		marchingCubes-> ComputeScalarsOff();
		marchingCubes-> SetValue(0,1);
		marchingCubes-> Update();

	marchingOutput = marchingCubes;
/*
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		//mapper->SetInputConnection(smooth->GetOutputPort());
		mapper->SetInputConnection( marchingCubes-> GetOutputPort());
		mapper->ScalarVisibilityOff();

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		//actor->GetProperty()->SetRepresentationToWireframe();

		actorOutput = actor;
*/
}

void MRimage::Gaussian_solver(std::vector<double> points_interface, float* &solution, int degree)
{	
	if(degree==1)
	{
			float b_new[3] = {0.0, 0.0, 0.0};
			for(int count=0; count<(points_interface.size()/3);count++){
				b_new[0] += points_interface[(3*count)+2];
				b_new[1] += points_interface[(3*count)]   * points_interface[(3*count)+2];
				b_new[2] += points_interface[(3*count)+1] * points_interface[(3*count)+2];
			}

			float A_new[3][3] = { {0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0} };
			for(int bis=0; bis<(points_interface.size()/3); bis++){
				A_new[0][0] += 1;
						A_new[0][1] += points_interface[3*bis];
								A_new[0][2] += points_interface[3*bis+1];
				
				A_new[1][0] = A_new[0][1];
						A_new[1][1] += points_interface[3*bis]*points_interface[3*bis];
								A_new[1][2] += points_interface[3*bis]*points_interface[3*bis+1];
				
				A_new[2][0] = A_new[0][2];
						A_new[2][1] = A_new[1][2];
								A_new[2][2] += points_interface[3*bis+1]*points_interface[3*bis+1];
			}

			int n = 3;
			solution = new float[3];
			//float x[3] = {0.0, 0.0, 0.0};
			float factor = 0;
			for( int k=0; k<n-1; k++){
				for(int i=k+1; i<n; i++){
					factor = A_new[i][k]/A_new[k][k];
					for( int j=k+1; j<n; j++){
						A_new[i][j] = A_new[i][j] - factor * A_new[k][j];
					}
					A_new[i][k] = 0.0;
					b_new[i] = b_new[i] - factor * b_new[k];
				}
			}
/*
			std::cout << " " << std::endl;
			std::cout << "b: " << b_new[0] << ", " << b_new[1]<< ", " << b_new[2]<< std::endl;
			std::cout << "" << std::endl;
			std::cout << "A: " << A_new[0][0] << ", " << A_new[0][1]<< ", " << A_new[0][2]<< std::endl;
			std::cout << "A: " << A_new[1][0] << ", " << A_new[1][1]<< ", " << A_new[1][2]<< std::endl;
			std::cout << "A: " << A_new[2][0] << ", " << A_new[2][1]<< ", " << A_new[2][2]<< std::endl;
			std::cout << " "<< std::endl;
/**/
			//Sleep(3000);

			solution[n-1]=b_new[n-1]/A_new[n-1][n-1];
			float sum = 0.0;
			for(int i=n-2; i>-1;  i--){
				sum = 0.0;
				for(int j=i+1; j<n; j++){
					sum += A_new[i][j]*solution[j];
				}
				solution[i] = (b_new[i]-sum)/A_new[i][i];
			}

			// Return point & solution !
			this->solution[0] = solution[0];
			this->solution[1] = solution[1];
			this->solution[2] = solution[2];
	}
	else if (degree==2)
	{
		// Orden 6. 1, x, y, xy, x^2, y^2
			float b_new[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			for(int count=0; count<(points_interface.size()/3);count++){
				b_new[0] += points_interface[(3*count)+2];
				b_new[1] += points_interface[(3*count)]*points_interface[(3*count)+2];
				b_new[2] += points_interface[(3*count)+1]*points_interface[(3*count)+2];

				b_new[3] += pow((points_interface[(3*count)]),2) *points_interface[(3*count)+2]; // x^2
				b_new[4] += pow((points_interface[(3*count)+1]),2) *points_interface[(3*count)+2]; // y^2
				b_new[5] += (points_interface[(3*count)] * points_interface[(3*count)+1]) *points_interface[(3*count)+2]; // xy
			}

			float A_new[6][6] = { {0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0} };
			for(int bis=0; bis<(points_interface.size()/3); bis++){
				A_new[0][0] += 1;
						A_new[0][1] += points_interface[3*bis];
								A_new[0][2] += points_interface[3*bis+1];
									A_new[0][3] += pow(points_interface[3*bis],2); // x^2
										A_new[0][4] += pow(points_interface[3*bis+1],2); // y^2
											A_new[0][5] += (points_interface[3*bis] * points_interface[3*bis+1]); // xy
				
				A_new[1][0] = A_new[0][1];
						A_new[1][1] += points_interface[3*bis]*points_interface[3*bis];
								A_new[1][2] += points_interface[3*bis]*points_interface[3*bis+1];
									A_new[1][3] += points_interface[3*bis]* pow(points_interface[3*bis],2); // x^2
										A_new[1][4] += points_interface[3*bis]* pow(points_interface[3*bis+1],2); // y^2
											A_new[1][5] += points_interface[3*bis]* (points_interface[3*bis] * points_interface[3*bis+1]); // xy
				
				A_new[2][0] = A_new[0][2];
						A_new[2][1] = A_new[1][2];
								A_new[2][2] += points_interface[3*bis+1]*points_interface[3*bis+1];
									A_new[2][3] += points_interface[3*bis+1]* pow(points_interface[3*bis],2); // x^2
										A_new[2][4] += points_interface[3*bis+1]* pow(points_interface[3*bis+1],2); // y^2
											A_new[2][5] += points_interface[3*bis+1]* (points_interface[3*bis] * points_interface[3*bis+1]); // xy

				A_new[3][0] = A_new[0][3];
						A_new[3][1] = A_new[1][3];
								A_new[3][2] = A_new[2][3];
									A_new[3][3] += pow(points_interface[3*bis],2) * pow(points_interface[3*bis],2); // x^2
										A_new[3][4] += pow(points_interface[3*bis],2) * pow(points_interface[3*bis+1],2); // y^2
											A_new[3][5] += pow(points_interface[3*bis],2) * (points_interface[3*bis] * points_interface[3*bis+1]); // xy

				A_new[4][0] = A_new[0][4];
						A_new[4][1] = A_new[1][4];
								A_new[4][2] = A_new[2][4];
									A_new[4][3] = A_new[3][4];
										A_new[4][4] += pow(points_interface[3*bis+1],2) * pow(points_interface[3*bis+1],2); // y^2
										A_new[4][5] += pow(points_interface[3*bis+1],2) * (points_interface[3*bis] * points_interface[3*bis+1]); // xy

				A_new[5][0] = A_new[0][5];
						A_new[5][1] = A_new[1][5];
								A_new[5][2] = A_new[2][5];
									A_new[5][3] = A_new[3][5];
										A_new[5][4] = A_new[4][5];
											A_new[5][5] += pow((points_interface[3*bis] * points_interface[3*bis+1]), 2); // xy
			}

			int n = 6;
			solution = new float[6];
			//float x[3] = {0.0, 0.0, 0.0};
			float factor = 0;
			for( int k=0; k<n-1; k++){
				for(int i=k+1; i<n; i++){
					factor = A_new[i][k]/A_new[k][k];
					for( int j=k+1; j<n; j++){
						A_new[i][j] = A_new[i][j] - factor * A_new[k][j];
					}
					A_new[i][k] = 0.0;
					b_new[i] = b_new[i] - factor * b_new[k];
				}
			}
/*
			std::cout << " " << std::endl;
			std::cout << "b: " << b_new[0] << ", " << b_new[1]<< ", " << b_new[2]<< std::endl;
			std::cout << "" << std::endl;
			std::cout << "A: " << A_new[0][0] << ", " << A_new[0][1]<< ", " << A_new[0][2]<< std::endl;
			std::cout << "A: " << A_new[1][0] << ", " << A_new[1][1]<< ", " << A_new[1][2]<< std::endl;
			std::cout << "A: " << A_new[2][0] << ", " << A_new[2][1]<< ", " << A_new[2][2]<< std::endl;
			std::cout << " "<< std::endl;
/**/
			//Sleep(3000);

			solution[n-1]=b_new[n-1]/A_new[n-1][n-1];
			float sum = 0.0;
			for(int i=n-2; i>-1;  i--){
				sum = 0.0;
				for(int j=i+1; j<n; j++){
					sum += A_new[i][j]*solution[j];
				}
				solution[i] = (b_new[i]-sum)/A_new[i][i];
			}

			// Return point & solution !
			this->solution[0] = solution[0];
			this->solution[1] = solution[1];
			this->solution[2] = solution[2];
			this->solution[3] = solution[3];
			this->solution[4] = solution[4];
			this->solution[5] = solution[5];
	}
}

void MRimage::controlFunction(std::vector<double> &points_interface, std::vector<int> &id_interface,
		std::vector<double> points_surface, std::vector<int> id_surface, float* solution, bool &eval, int degree)
{
	//std::cout <<"S_0 size  " << S_0.size() << std::endl;
	double control;
	eval = false;
	int count =0;
	int count2 =0;
	bool is_inside = false;
	
	for(int not_rep=0; not_rep<id_surface.size(); not_rep++){
		is_inside =  false;
		for(int not_rep2=0; not_rep2<id_interface.size(); not_rep2++){
			if( id_interface[not_rep2] == id_surface[not_rep]) is_inside = true;
		}
		count+=1;
		if(!is_inside){
			if(degree==1) {
				control = solution[0] + points_surface[ 3*not_rep ] * solution[1] + points_surface[ 3*not_rep+1 ] *solution[2];
				if((points_surface[ 3*not_rep+2 ] - control) > 0)
				{  // antes >0
					//points_interface.push_back(1);
					points_interface.push_back( points_surface[ 3*not_rep ] );
					points_interface.push_back( points_surface[ 3*not_rep+1 ] );
					points_interface.push_back( points_surface[ 3*not_rep+2 ] );
					id_interface.push_back( id_surface[not_rep]  );

					eval = true;
					count2+=1; 
				}
			}else if(degree==2)	{
				control = solution[0] + points_surface[ 3*not_rep ]*solution[1] +
					points_surface[ 3*not_rep+1 ]*solution[2] + pow(points_surface[3*not_rep],2)*solution[3] +
					pow(points_surface[3*not_rep+1],2)*solution[4] + (points_surface[3*not_rep]*points_surface[3*not_rep+1])*solution[3];

				if((points_surface[ 3*not_rep+2 ] - control) > 0)
				{  // antes >0
					//points_interface.push_back(1);
					points_interface.push_back( points_surface[ 3*not_rep ] );
					points_interface.push_back( points_surface[ 3*not_rep+1 ] );
					points_interface.push_back( points_surface[ 3*not_rep+2 ] );
					id_interface.push_back( id_surface[not_rep]  );

					eval = true;
					count2+=1;
				}
			}	
		}		
	}
/*
	std::cout << "Puntos "<< count << " de " << connect_surface.size() <<" Con " << count2 << " nuevos  Puntos introducidos" << std::endl;
	std::cout << "Size A_0 = " << connect_interface.size() << std::endl;
/**/
}

void MRimage::upAndDownCut(ImageMriType::Pointer image, ImageMriType::Pointer mask, ImageMriType::Pointer &outputimage)
{
	ImageMriType::Pointer newImage = ImageMriType::New();
		newImage->SetOrigin( image->GetOrigin() );
		newImage->SetRegions( image->GetLargestPossibleRegion() );
		newImage->SetSpacing( image->GetSpacing() );
		newImage->SetDirection( image->GetDirection() );
		newImage->Allocate();
		newImage->FillBuffer(1);

	LabelType::Pointer stats = LabelType::New();
		stats->SetLabelInput(mask);
		stats->SetInput(image);
		stats->Update();

	LabelType::BoundingBoxType bb = stats->GetBoundingBox( 1 );

	itk::ImageRegionIterator<ImageMriType> it = itk::ImageRegionIterator<ImageMriType>(newImage, newImage->GetLargestPossibleRegion() );
	for(it.Begin(); !it.IsAtEnd(); it++) {
		ImageMriType::IndexType idx = it.GetIndex();
		if( idx[0]<bb[0] || idx[0]>bb[1] || idx[1]<bb[2] || idx[1]>bb[3] ) it.Set( 0 );
	}

	outputimage = newImage;
}

// -------------------------------------------------------------
void MRimage::getBreastSurfaceModel( const char* side, vtkSmartPointer<vtkPolyData> & surf_mesh)
{
	ImageMriType::Pointer temp_segmentation = ImageMriType::New();
	ImageMriType::Pointer temp_mask = ImageMriType::New();

	getRegionOfInterest(this->m_albertsegmentation, side, temp_segmentation);
	getRegionOfInterest(this->m_breastmask, side, temp_mask);

	/* ************************************************************************************************ */	

	ThresholdMriType::Pointer threshold = ThresholdMriType::New();
		threshold -> SetInput( temp_segmentation );
		threshold -> SetInsideValue( 1 );
		threshold -> SetOutsideValue( 0 );
		// Hay que usar las segmentaciones del cuerpo completo de Albert!!
		threshold -> SetLowerThreshold( 2 ); // Con las segmentaciones de Albert: Threshold = 2;
		threshold -> Update();
			
	/* ************************************************************************************************ */	
	
	ImageMriType::Pointer newImage = ImageMriType::New();
		newImage->SetOrigin( threshold->GetOutput()->GetOrigin() );
		newImage->SetRegions( threshold->GetOutput()->GetLargestPossibleRegion() );
		newImage->SetSpacing( threshold->GetOutput()->GetSpacing() );
		newImage->SetDirection( threshold->GetOutput()->GetDirection() );
		newImage->Allocate();
		newImage->FillBuffer(0);

	upAndDownCut( threshold->GetOutput(), temp_mask, newImage);

	AndMriFilter::Pointer andA = AndMriFilter::New();
		andA->SetInput( 0, threshold->GetOutput() );
		andA->SetInput( 1, newImage );
		andA->Update();

	/* ************************************************************************************************ */	
	ImageMriType::Pointer &resamp_breast = ImageMriType::New();
	ImageMriType::Pointer &resamp_body = ImageMriType::New() ;

	//resampleFunction(threshold->GetOutput(), this->iso_vox_size, resamp_body);
	resampleFunction(andA->GetOutput(), this->iso_vox_size, resamp_body);
	resampleFunction(temp_mask, this->iso_vox_size, resamp_breast);


	/* ************************************************************************************************ */	
	/* Guarda este código por si acaso, puede ser interesante comprobar el resto de las imágenes */
	/*
	typedef itk::AddImageFilter<ImageType, ImageType> AddType;

	AddType::Pointer addfilter = AddType::New();
		addfilter->SetInput1( segmentation_2 );
		addfilter->SetInput2( reader_mask->GetOutput() );
		addfilter->Update();
		*/

	/* ************************************************************************************************ */	

	SubstractMriType::Pointer substractfilter = SubstractMriType::New();
		substractfilter -> SetInput1( resamp_body );
		substractfilter -> SetInput2( resamp_breast );
		substractfilter -> Update();

	ThresholdMriType::Pointer threshold2 = ThresholdMriType::New();
		threshold2 -> SetInput( substractfilter->GetOutput() );
		threshold2 -> SetInsideValue( 1 );
		threshold2 -> SetOutsideValue( 0 );
		threshold2 -> SetLowerThreshold( 1 );
		threshold2 -> Update();
	
	/* ************************************************************************************************ */	
	// closing
	ImageMriType::Pointer &closing_breast = ImageMriType::New();
	ImageMriType::Pointer &closing_body = ImageMriType::New();

	BallType ballKernel;
		ballKernel.SetRadius( 2 );
		ballKernel.CreateStructuringElement();
//		std::cout << ballKernel.GetRadius() << std::endl;

	closingFunction(resamp_breast, ballKernel, closing_breast, 1);
	closingFunction(threshold2->GetOutput(), ballKernel, closing_body,0);

	// dilation
	BallType ballKernel2;
		ballKernel2.SetRadius( 1 );
		ballKernel2.CreateStructuringElement();

	DilateMriType::Pointer dilateFilter = DilateMriType::New();
		dilateFilter-> SetInput( closing_body );
		dilateFilter-> SetKernel( ballKernel2 );
		dilateFilter-> SetDilateValue( 1 );
		dilateFilter-> Update();  
	/* ************************************************************************************************ */	
	AndMriFilter::Pointer andfilter = AndMriFilter::New();
		andfilter->SetInput(0, dilateFilter->GetOutput() );
		andfilter->SetInput(1, closing_breast );
		andfilter->Update();
		
	/* ************************************************************************************************ */	
	/*
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(andfilter->GetOutput() );
		//writer -> SetInput(andfilter->GetOutput() );
		writer -> SetFileName( outputfilename );
		try {
			writer->Update();   }
		catch (itk::ExceptionObject & err)   {
			std::cerr << " " << std::cerr;
			std::cerr << "Writer Fails!" << std::cerr;
			std::cerr << err << std::cerr;
			// return EXIT_FAILURE;
		}
	*/

/* 
*  Ahora hay que pasar de itk a vtk, hacer un dataset de nodos y optimizar por mínimos cuadrados.
*/
	vtkSmartPointer<vtkMarchingCubes> &marching_body = vtkSmartPointer<vtkMarchingCubes>::New();
	vtkSmartPointer<vtkMarchingCubes> &marching_breast = vtkSmartPointer<vtkMarchingCubes>::New();

	visualizationFunction( andfilter->GetOutput(),  marching_body); 
	visualizationFunction( closing_breast, marching_breast); 

	/* *** Point SET *** */
  
	vtkPoints* allPoints_body = marching_body->GetOutput()->GetPoints();
	vtkPoints* allPoints_breast = marching_breast->GetOutput()->GetPoints();
  
/*

*/	
	double * temp_body = new double[3]; // = {0.0, 0.0, 0.0};
	double * temp_breast = new double[3]; //= {0.0, 0.0, 0.0};
	
	std::vector<double> points_interface;
	std::vector<int> id_interface;

	std::vector<double> points_surface;
	std::vector<int> id_surface;

	for(vtkIdType id_breast = 0; id_breast<allPoints_breast->GetNumberOfPoints(); id_breast++){
		allPoints_breast->GetPoint(id_breast, temp_breast );

		//points_surface.push_back(1);
		points_surface.push_back(temp_breast[0]);
		points_surface.push_back(temp_breast[1]);
		points_surface.push_back(temp_breast[2]);
		
		id_surface.push_back( id_breast+1 );
		
		for(vtkIdType id_body = 0; id_body<allPoints_body->GetNumberOfPoints(); id_body++){
			allPoints_body-> GetPoint( id_body, temp_body );
			if( temp_body[0]==temp_breast[0] && temp_body[1]==temp_breast[1]  && temp_body[2]==temp_breast[2] ){
				//points_interface.push_back(1);
				points_interface.push_back(temp_body[0]);
				points_interface.push_back(temp_body[1]);
				points_interface.push_back(temp_body[2]);
				
				id_interface.push_back( id_breast+1 );
			}
		 // Hasta aqui BIEN!! ******
		}
	}

	const vtkIdType numCells = marching_breast->GetOutput()->GetNumberOfCells();
	std::vector<int> connect_surface;
	
	for(vtkIdType cellId=0 ; cellId<numCells ; ++cellId)
	{		
		// Get all cells and its points.
		vtkCell* cell = marching_breast-> GetOutput()-> GetCell(cellId);

		if( cell!=NULL ) { 
			vtkIdList* list;
			list = cell->GetPointIds();
			connect_surface.push_back( list->GetId( 0 ) );
			connect_surface.push_back( list->GetId( 1 ) );
			connect_surface.push_back( list->GetId( 2 ) );
			//std::cout<< list->GetId(0)<< ", " << list->GetId(1)<< ", " << list->GetId(2)<< std::endl;
		}
	}

	float *x;
	// int degree=1; // <------------- HAY QUE CAMBIARLO MAS ADELANTE!!!
	int degree = this->m_degree;

	if(degree==1) { std::cout<<std::endl; std::cout << "Degreeeee = 1" << std::endl; this->solution = new float[3]; }
	if(degree==2) { std::cout<<std::endl; std::cout << "Degreeeee = 2" << std::endl; this->solution = new float[6]; }

	bool eval = true;
	int count = 0;
	//Gaussian_solver(A_0,A_x,A_y, b_z, x);
	Gaussian_solver( points_interface, x, degree);
	while( eval && count<4 ){
		controlFunction(points_interface, id_interface, points_surface, id_surface, x, eval, degree);
		Gaussian_solver(points_interface, x, degree);
		count++;
	}

	std::cout << "Gaussian solution: [" <<  x[0] << ", " <<  x[1] << ", " << x[2] << std::endl;
/*
	 Recomponer los puntos:
	 b = [A_0, A_x, A_y] * x;
*/
	double pt = 0.0;
	float * temp_point = new float[3];// = {0.0,0.0,0.0};

	if(degree==1){
		for(int last =0; last<(points_surface.size()/3); last++){ 
			for(int ll =0; ll<(points_interface.size()/3); ll++){
				if(id_surface[last] == id_interface[ll]){
					pt = x[0] + points_interface[3*ll] * x[1] + points_interface[3*ll+1] * x[2];
					//points_surface[4*last+1] = points_interface[4*ll+1];
					//points_surface[4*last+2] = points_interface[4*ll+2];
					points_surface[3*last+2] = pt;
				}
			}
		}
	} else if(degree==2){
		for(int last =0; last<(points_surface.size()/3); last++){ 
			for(int ll =0; ll<(points_interface.size()/3); ll++){
				if(id_surface[last] == id_interface[ll]){
					pt = x[0] + points_interface[3*ll] * x[1] + points_interface[3*ll+1] * x[2] +
						pow(points_interface[3*ll],2) * x[3] + pow(points_interface[3*ll+1],2) * x[4] + (points_interface[3*ll]*points_interface[3*ll+1]) * x[5];
					//points_surface[4*last+1] = points_interface[4*ll+1];
					//points_surface[4*last+2] = points_interface[4*ll+2];
					points_surface[3*last+2] = pt;
				}
			}
		}
	}

// ======================= New PolyData Mesh ======================
	vtkSmartPointer<vtkPoints> mesh_points = vtkSmartPointer<vtkPoints>::New();

	for(int last_count = 0; last_count<id_surface.size(); ++last_count){ // recordar cambiar
		temp_point[0]=points_surface[3*last_count];
		temp_point[1]=points_surface[3*last_count+1];
		temp_point[2]=points_surface[3*last_count+2];

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
		

// ================  return surface mesh ===============
	surf_mesh = polygonPolyData;

		// == Destroying temporal arrays !
		delete [] temp_body; 
		delete [] temp_breast;
	
		points_interface.~vector();
		id_interface.~vector();

		points_surface.~vector();
		id_surface.~vector();

		connect_surface.~vector();

		delete [] temp_point;
}

// Meshing
// --------------------------------------------------------
void MRimage::getmodel(vtkSmartPointer<vtkPolyData> &total_mesh )
{
	vtkSmartPointer<vtkPolyData> surf_mesh = vtkSmartPointer<vtkPolyData>::New();
	getBreastSurfaceModel( m_side, surf_mesh);

		new_mesh->set_surfacemesh( surf_mesh );
		new_mesh->improve_quality( surf_mesh );
	
/*
		std::cout << "mapeando la superficie" << std::endl;
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputData(surf_mesh);

	double b[3] = {1.0, 1.0, 1.0};
  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor( b );
		//actor->GetProperty()->SetRepresentationToWireframe();

	vtkSmartPointer<vtkPolyDataMapper> mapper_2 = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper_2->SetInputData(surf_mesh);
	double a[3] = {1.0, 0.0, 0.0};
	vtkSmartPointer<vtkActor> actor_2 = vtkSmartPointer<vtkActor>::New();
		actor_2->SetMapper(mapper_2);
		actor_2->GetProperty()->SetColor( a );
		actor_2->GetProperty()->SetRepresentationToWireframe();
 
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
		renderer->AddActor(actor);
		renderer->AddActor(actor_2);
		renderer->ResetCamera();
		renderer->SetBackground(0.2, 0.4, 0.5);
// vtkSmartPointer<vtkRenderer> renderer2 = vtkSmartPointer<vtkRenderer>::New();
//		renderer2->AddActor(actor_2);
//		renderer2->ResetCamera();
//		renderer2->SetBackground(0.2, 0.4, 0.5);
  
	vtkSmartPointer<vtkRenderWindow> renderWindow2 = vtkSmartPointer<vtkRenderWindow>::New();
		renderWindow2->AddRenderer(renderer);
		//renderWindow2->AddRenderer( renderer );
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor2 = vtkSmartPointer<vtkRenderWindowInteractor>::New();
		renderWindowInteractor2->SetRenderWindow(renderWindow2);

  renderWindow2->Render();
  renderWindowInteractor2->Initialize();
  renderWindowInteractor2->Start(); 
/* */

	new_segmentation->SetSide( m_side );
	new_segmentation->SetDegree( m_degree );
	new_segmentation->SetFinalSpacing( iso_vox_size );
	new_segmentation->SetSolution( solution );

	vtkSmartPointer<vtkPolyData> t_mesh = vtkSmartPointer<vtkPolyData>::New();
		new_segmentation->JoinMeshes(surf_mesh, t_mesh);
	
/* Visualización */
/* 
		std::cout << "Visualización de unión de la malla" << std::endl;
	double a2[3] = {0.0, 0.0, 0.0};
	vtkSmartPointer<vtkPolyDataMapper> mappera = vtkSmartPointer<vtkPolyDataMapper>::New();
		//mappera->SetInputData( surf_mesh );
		mappera->SetInputData( t_mesh );

	vtkSmartPointer<vtkActor> actora = vtkSmartPointer<vtkActor>::New();
		actora->SetMapper( mappera );
		actora->GetProperty()->SetColor( a2 );
		actora->GetProperty()->SetRepresentationToWireframe();

	vtkSmartPointer<vtkRenderer> renderera = vtkSmartPointer<vtkRenderer>::New();
		renderera->AddActor( actora );
		renderera->ResetCamera();
		renderera->SetBackground( 1.0, 1.0, 1.0 );

	vtkSmartPointer<vtkRenderWindow> renderWindowa = vtkSmartPointer<vtkRenderWindow>::New();
		renderWindowa->AddRenderer( renderera );

	vtkSmartPointer<vtkRenderWindowInteractor> Riteratora = vtkSmartPointer<vtkRenderWindowInteractor>::New();
		Riteratora->SetRenderWindow( renderWindowa );

	renderWindowa->Render();
	Riteratora->Initialize();
	Riteratora->Start();
/**/

		//new_mesh->improve_quality(t_mesh);
		//new_mesh.tetgen_ext(surf_mesh);
		new_mesh->set_boundaryConditions( solution );
		new_mesh->set_degree( m_degree );
		
		new_mesh->tetgen_ext(t_mesh); // Esta es la ultima, tiene que llamarse Update()

	total_mesh = t_mesh;
	istheremodel = true;
}

std::vector<double> MRimage::getPoints_model()
{
	if( istheremodel ) 
	{
		return new_mesh->m_initial_points;
	}
	else
	{
		std::cout << "No hay modelo! Revisa el código ! " << std::endl;
		std::vector<double> kk;
		return kk;
	}

}

std::vector<int> MRimage::getElements_model()
{
	if( istheremodel ) 
	{
		return new_mesh->m_elements;
	}
	else
	{
		std::cout << "No hay modelo! Revisa el código ! " << std::endl;
		std::vector<int> kk;
		return kk;
	}

}

std::vector<int>  MRimage::getTissueElements_model()
{
	if( istheremodel )
	{
		return new_mesh->m_tetrahedron_tissue;
	}
	else
	{
		std::cout << "No hay modelo! Revisa el código ! " << std::endl;
		std::vector<int> kk;
		return kk;
	}
}

std::vector<int> MRimage::getBoundaryConditions_model()
{
	if( istheremodel )
	{
		return new_mesh->m_boundaryConditions;
	}
	else
	{
		std::cout << "No hay modelo! Revisa el código ! " << std::endl;
		std::vector<int> kk;
		return kk;
	}
}

ImageSegmentation3floatType::Pointer MRimage::getImage_NewSegmentation()
{
	return new_segmentation->GetImagen();
}

#endif