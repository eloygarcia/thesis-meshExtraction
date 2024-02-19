#ifndef __mesh_cpp
#define __mesh_cpp

#include "mesh.h"

mesh::mesh(void)
{
}

mesh::~mesh(void)
{
//	m_points_surface.~vector();
//	m_cell_vector.~vector();
//	m_cell_element_list.~list();

//	m_initial_points.~vector();
//	m_final_points.~vector();
//	m_elements.~vector();
//	m_cell_element_list.~list();

//	m_tetrahedron_tissue.~vector();
//	m_intermedial_compression_points.~vector();

//	m_boundaryConditions.~vector();
}

// ===================== Surface mesh ======================
void mesh::set_surfacemesh(vtkSmartPointer<vtkPolyData> surface_mesh)
{
	std::cout << "set_surfaceMesh" << std::endl;
	// include surface mesh 
	this->m_surface_mesh = surface_mesh;
/**/
	vtkSmartPointer<vtkSmoothPolyDataFilter> smooth = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
		smooth->SetRelaxationFactor(0.01);
		smooth->SetNumberOfIterations(100);
		smooth->SetInputData( this->m_surface_mesh );
		smooth->Update();

	// ==================== Return Surface mesh ========================
	surface_mesh = smooth->GetOutput();
/**/
	// ========================= Points ============================
	vtkSmartPointer<vtkPoints> points = surface_mesh->GetPoints();
	double pt[3] = {0.0,0.0,0.0};
	std::vector<double> points_vector;
	for(vtkIdType id_p =0; id_p<points->GetNumberOfPoints(); id_p++)
	{
		 points->GetPoint(id_p,pt);
		 //poly_points->InsertNextPoint(pt);
		 points_vector.push_back( pt[0] );
		 points_vector.push_back( pt[1] );
		 points_vector.push_back( pt[2] );
	}
	this->m_points_surface = points_vector;
	//std::cout << "Points vector intial size: " << points_vector.size() << std::endl;
	
	// ========================= connectivity ========================
	std::vector<int> cell_vector;
	for(vtkIdType id_c=0; id_c<surface_mesh->GetNumberOfCells(); id_c++)
	{
		vtkCell* cell = surface_mesh->GetCell(id_c);
		vtkIdList* list = cell->GetPointIds();
		cell_vector.push_back( list->GetId(0) );
		cell_vector.push_back( list->GetId(1) );
		cell_vector.push_back( list->GetId(2) );
		// -------------------------------------
		this->m_cell_unique_list.push_back( list->GetId(0) );
		this->m_cell_unique_list.push_back( list->GetId(1) );
		this->m_cell_unique_list.push_back( list->GetId(2) );
	}
	this->m_cell_vector = cell_vector;
	
	this->m_cell_unique_list.sort();
	this->m_cell_unique_list.unique();
	//std::cout << "Cell vector initial size: " << cell_vector.size() << std::endl;

	// ==== Delete ====
//	points_vector.~vector();
//	cell_vector.~vector();
}

// =================== quality =======================
void mesh::improve_quality(vtkSmartPointer<vtkPolyData> &surface_mesh)
{
	std::cout << "improveQuality" << std::endl;
	std::list<int> borderlist;
	cleaningMesh * cle = new cleaningMesh();
	// AQUI HEMOS CAMBIADO EL CODIGO
	bool repeat = true;
	while(repeat) {
		cle->removeDuplicatePoints(this->m_points_surface, this->m_cell_vector);
		cle->removeDuplicatedCells(this->m_cell_vector);
	
		cle->foundNonManifold(this->m_cell_vector, borderlist);
		cle->clearMesh(this->m_cell_vector, borderlist, repeat);
		cle->performNewMesh(this->m_points_surface, this->m_cell_vector);
		borderlist.clear();
		//borderlist.~list();

	// ============================ Hasta Aqui ================

	double pt_m[3];// = {0.0, 0.0, 0.0};
	vtkSmartPointer<vtkPoints> poly_points = vtkSmartPointer<vtkPoints>::New();
	for(int id_p=0; id_p<this->m_points_surface.size()/3; id_p++)
	{
		pt_m[0] = this->m_points_surface[id_p*3];
		pt_m[1] = this->m_points_surface[id_p*3+1];
		pt_m[2] = this->m_points_surface[id_p*3+2];
		poly_points->InsertNextPoint(pt_m);
	}
	
//	std::cout << "Points vector final size: " << poly_points->GetNumberOfPoints() << std::endl;

	vtkSmartPointer<vtkCellArray> cell_array = vtkSmartPointer<vtkCellArray>::New();
	for(int c=0; c<this->m_cell_vector.size()/3; c++)
	{
		cell_array->InsertNextCell( 3 );
		cell_array->InsertCellPoint( this->m_cell_vector[3*c] );
		cell_array->InsertCellPoint( this->m_cell_vector[3*c+1] );
		cell_array->InsertCellPoint( this->m_cell_vector[3*c+2] );
	}
	
//	std::cout << "Cell vector final size: " << cell_array->GetNumberOfCells() << std::endl;
	
	// ======================= VISUALIZACION DE PUNTOS ==============	

	vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
		polyData->SetPoints( poly_points );
		polyData->SetPolys( cell_array );
		polyData->BuildCells();
/**/

	vtkSmartPointer<vtkSmoothPolyDataFilter> smooth = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
		smooth->SetRelaxationFactor(0.02); // antes 0.01	
		smooth->SetNumberOfIterations(150); // 100
		smooth->SetInputData( polyData );
		smooth->Update();

	// ==================== Return Surface mesh ========================
	surface_mesh = smooth->GetOutput();

/*
	vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivity_filter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
		connectivity_filter->SetInputData( smooth->GetOutput() );
		connectivity_filter->SetExtractionModeToLargestRegion();
		connectivity_filter->Update();

	vtkSmartPointer<vtkCleanPolyData> clean = vtkSmartPointer<vtkCleanPolyData>::New();
		clean->SetInputData( connectivity_filter->GetOutput() );
		clean->Update();

/*

	vtkSmartPointer<vtkPolyDataWriter> writer_grid =vtkSmartPointer<vtkPolyDataWriter>::New();
		writer_grid->SetInputData( smooth->GetOutput() );
		writer_grid->SetFileName( "C:\\Users\\Eloy\\Desktop\\prueba_20_10.vtk");
		writer_grid->Update();
	/*

	vtkSmartPointer<vtkSmoothPolyDataFilter> smooth = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
		smooth->SetRelaxationFactor(0.01);
		smooth->SetNumberOfIterations(100);
		smooth->SetInputData( connectivity_filter->GetOutput() );
		smooth->Update();

	// ==================== Return Surface mesh ========================
	surface_mesh = smooth->GetOutput();

/*
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputData( smooth->GetOutput() );
		
		double color[3] = {1, 0, 0};
	vtkSmartPointer<vtkProperty> prop = vtkSmartPointer<vtkProperty>::New();
		prop->SetRepresentationToWireframe();
		prop->SetColor(color);

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper( mapper );
		actor->SetProperty(prop);

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
		renderer->AddActor( actor );
		renderer->ResetCamera();
		renderer->SetBackground(0.2, 0.4, 0.5);
  
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
		renderWindow->AddRenderer(renderer);
	
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
		renderWindowInteractor->SetRenderWindow(renderWindow);

	renderWindow->Render();
	renderWindowInteractor->Initialize();
	renderWindowInteractor->Start(); 
/**/
	 }

	 // ==== Delete ====
//	 delete cle;
//	 borderlist.~list();

}

// =================== Create unstructured Grid =======================
void mesh::createUnstructuredGrid( std::vector<double> point_vector, 
								  std::vector<int> cells_vector,
								  std::vector<int> bC_list,
								  std::vector<int> tissue_list,
								  vtkSmartPointer<vtkUnstructuredGrid> &grid)
{
	std::cout << "createUnstructuredGrid" << std::endl;
	int numberofpoints = point_vector.size()/3;
	int numberofcells = cells_vector.size()/4;

	// ========================================================
	vtkPoints * points = vtkPoints::New();
	vtkCellArray * cells = vtkCellArray::New();

	points->SetNumberOfPoints( numberofpoints );
	double * auxPoint = new double[3];
		auxPoint[0] = 0.0;
		auxPoint[1] = 0.0;
		auxPoint[2] = 0.0;

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

	vtkIntArray * data = vtkIntArray::New();
		data->SetName("PectoralBC");
	bool is_inside = false;
	for(int bc = 0; bc< numberofpoints; bc++ )	{
		is_inside = false;
		for(int j=0; j<bC_list.size(); j++){
			if (bc==bC_list[j]) 
			{
					is_inside = true;
					std::cout << "Is INSIDE " << bc << std::endl;
			}
		}

		if(is_inside) data->InsertNextValue(1);
		else data->InsertNextValue(0);
	}
	
	//vtkIdType * id_tissue;
	vtkIntArray * tissues = vtkIntArray::New();
		tissues->SetName("Tissue");
	for( int t=0; t<tissue_list.size(); t++){
		tissues->InsertNextValue(  tissue_list[t] );
	}

	grid->SetPoints(points); // Points
	grid->SetCells(VTK_TETRA, cells); // Cells

	grid->GetPointData()->SetScalars( data ); // BoundaryConditions
	grid->GetCellData()->SetScalars( tissues ); // Tissues
}

// ==================== tetgen mesh generator =======================
void mesh::tetgen_ext(vtkSmartPointer<vtkPolyData> surface_mesh)
{
	std::cout << std::endl;
	std::cout << " ======================== TETGEN ======================== " << std::endl;
	std::cout << std::endl;
	// Tetgen

	tetgenio in, out;
	tetgenio::facet *f;
	tetgenio::polygon *p;
	//int i;

	// All indices start from 1.
	in.firstnumber=1;

	//vtkPoints* allPoints = polyData->GetPoints();
	vtkPoints* allPoints = surface_mesh->GetPoints();

	const vtkIdType numPolyDataPoints = allPoints->GetNumberOfPoints();
	
	// Inicializar lista de puntos.
	//in.numberofpoints = numPolyDataPoints;
	//in.pointlist = new REAL[in.numberofpoints*3];

	const vtkIdType numCells = surface_mesh->GetNumberOfCells();

	in.numberofpoints = numCells*3;
	in.pointlist = new REAL[in.numberofpoints*3];
//	std::cout <<" " << std::endl;
//	std::cout <<"Number of points: " << in.numberofpoints << std::endl;

	//in.numberoffacets = 1;
	in.numberoffacets = numCells;
//	std::cout <<"Number of facets: " << in.numberoffacets << std::endl;
	in.facetlist = new tetgenio::facet[ in.numberoffacets ];
	in.facetmarkerlist = new int [ in.numberoffacets ];

	//in.numberofregions(2);

	/* ESTA PARTE TAL VEZ DEBA SER CAMBIADA */

			 // Facet (cellId)
			//f= &in.facetlist[ 0 ];
			//f-> numberofpolygons = numCells;
			//f-> polygonlist = new tetgenio::polygon[ f-> numberofpolygons ];
			// Holes.
			//f-> numberofholes = 0;
			//f-> holelist = NULL;
			
			//std::cout << " " << std::endl;
			//std::cout << "Antes de introducir los valores en las celdas" << std::endl;

	int count=0;
	double pt[3]={0.0, 0.0, 0.0};

	// Introducir celdas y puntos.
	for(vtkIdType cellId=0 ; cellId<numCells ; ++cellId)
	{		
		// Get all cells and its points.
		vtkCell* cell = surface_mesh-> GetCell(cellId);

		if( cell!=NULL ) {

			/// LAS LINEAS QUE FALLAN SON ESTAS. PONIENDO LAS LINEAS DE ARRIBA LEE TODOS LOS PUNTOS Y LOS INTRODUCE EN LA LISTA. ES AL CREAR LAS FACETS.

			// Facet (cellId)
			f= &in.facetlist[ cellId ];
		//std::cout << cellId << std::endl;
			f-> numberofpolygons = 1;
			f-> polygonlist = new tetgenio::polygon[ f->numberofpolygons ];
			// Holes.
			f-> numberofholes = 0;
			f-> holelist = NULL;
			vtkPoints* points = cell-> GetPoints();
			p = &f -> polygonlist[ 0 ];  // Puta tonteria. Era este punto.  Puse cellId en lugar de 0.!!! IDIoTA.
			
			// Get all points of the cell.
			const vtkIdType numCellPoints = points->GetNumberOfPoints();
			p-> numberofvertices = numCellPoints ;
			p-> vertexlist = new int[p->numberofvertices];
			
			if(numCellPoints==3) {
				//std::cout << cellId << std::endl;
			for(vtkIdType pointId = 0 ; pointId<numCellPoints ; ++pointId)
			{
				points-> GetPoint(pointId, pt);

				// Conseguir la actual id del punto y escribirlo en su posicion.
				vtkIdType actualId = cell -> GetPointId( pointId );
				
					//in.pointlist[actualId*3] = pt[0];
					//in.pointlist[actualId*3+1] = pt[1];
					//in.pointlist[actualId*3+2] = pt[2];
				
					in.pointlist[count*3] = pt[0];
					in.pointlist[count*3+1] = pt[1];
					in.pointlist[count*3+2] = pt[2];

				//std::cout << " " << std::endl;
				
				p-> vertexlist[pointId] = count+1; 
				//p-> vertexlist[pointId] = actualId;

				//std::cout << actualId << std::endl;	
				//std::cout << "[ " << in.pointlist[actualId*3] << ", " << in.pointlist[actualId*3+1] << ", " <<in.pointlist[actualId*3+2] << "] " << std::endl;
				//std::cout << " " << std::endl;

				++count;
				//std::cout<<count<<std::endl;
			}
			
			in.facetmarkerlist[cellId] = 0;
			}
		}
	}

	//std::cout << "Total Points: " << count << std::endl;

	//in.facetmarkerlist[0] = 1;

	// std::cout << "Tetgen Correcto!!" << std::endl;

	// Comprobacion del codigo
	//in.save_nodes("barin");
	//in.save_poly("barin");
	//in.save_faces("barin");
	//std::cout << "Saved In" << std::endl;

	// Tetrahedralize
	//tetrahedralize("pAYVq1.5a12g", &in, &out);  // Tengo que buscar el punto de los 100.000 elementos
	tetrahedralize("pAYVq1.5a10g", &in, &out);  // -V -> verbose;
	//tetrahedralize("pAYq1.5a8g", &in, &out);  
	//tetrahedralize("p", &in, &out);

	// Output mesh.
	//out.save_nodes("barout");
	//out.save_elements("barout");
	//out.save_faces("barout");
	//std::cout << "Saved Out" << std::endl;
	std::cout << std::endl;

	// WRITE MODEL

	std::vector<double> nodesVect;
	std::vector<int> elementsVect;
	std::vector<int> constraintVect;

	for(int i=0; i<out.numberofpoints*3; ++i){
		nodesVect.push_back( out.pointlist[ i ] );
	}

	for(int i=0; i<out.numberoftetrahedra*4; ++i){
		//std::cout << "Tetra = " << out.tetrahedronlist[ i] -1 << std::endl;
		elementsVect.push_back( out.tetrahedronlist[ i ]-1 );
		this->m_cell_element_list.push_back( out.tetrahedronlist[ i ]-1 );
	}

	std::vector<int> elementattribute;

	for(int i=0; i<out.numberoftetrahedra; i++) {
		// std::cout << out.tetrahedronattributelist[ i ] << std::endl;
		elementattribute.push_back( out.tetrahedronattributelist[i] );
	}
	
	//std::cout << "Numero de Regiones ! : " << out.numberofregions << std::endl;

	// ========================== Return nodes & elements ===================
//	std::cout << "Number of Points: " << nodesVect.size()/4 << std::endl;
//	std::cout << "Number of Tetrahedral Elements: " << elementsVect.size()/4 << std::endl;

	this->m_initial_points = nodesVect;
	this->m_final_points = nodesVect;
	
	this->m_elements = elementsVect; 
	this->m_tetrahedron_tissue = elementattribute; // Tipo de tejido para el tetrahedro[i]

	this->m_cell_element_list.sort();
	this->m_cell_element_list.unique();

	// out.save_elements( "C:\\Users\\Eloy\\Desktop\\elements.ele");
	// out.save_nodes( "C:\\Users\\Eloy\\Desktop\\nodes.node");

	std::cout << std::endl;
	std::cout << "Boundary Conditions Size : " << m_boundaryConditions.size() << std::endl;
	get_boundaryConditions();
	std::cout << "Boundary Conditions Size : " << m_boundaryConditions.size() << std::endl;
	std::cout << std::endl;

	/* CONSTRUYENDO LA MALLA DE TETRAHEDROS PARA DEVOLVER */
	vtkSmartPointer<vtkUnstructuredGrid> unstrGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	createUnstructuredGrid( nodesVect, 
								elementsVect,
								m_boundaryConditions,
								m_tetrahedron_tissue,
								unstrGrid);


	this->m_tetrahedral_mesh = unstrGrid;
	/**/

	// ==== Delete ====
//	nodesVect.~vector();
//	elementsVect.~vector();
//	constraintVect.~vector();
//
//	elementattribute.~vector();

}

// ==================== chossing boundary conditions =======================
// ------------------------ Pectoral muscle nodes -----------------------
void mesh::get_boundaryConditions()
{
	std::cout << "Get_boundaryConditions" << std::endl;
//	std::cout << "Solution: " << solution[0] << ", " << solution[1] << ", " << solution[2] << std::endl;

	double * temp_point = new double[3];
		temp_point[0] = 0.0;
		temp_point[1] = 0.0;
		temp_point[2] = 0.0;

	float dist = 0.0;
	float min = 9999;
	int count = 0;

	for(int i=0; i<(this->m_initial_points.size()/3); i++)
	{
		temp_point[0] = this->m_initial_points[3*i];
		temp_point[1] = this->m_initial_points[3*i+1];
		temp_point[2] = this->m_initial_points[3*i+2];
		if (m_degree==1) {
			//float temp = (solution[0] + solution[1] * temp_point[0] + solution[2] * temp_point[1]);
			// Distancia de un punto a un plano.
			dist = temp_point[2] - (m_solution[0] + m_solution[1] * temp_point[0] + m_solution[2] * temp_point[1]);

			if(-dist<min) min=-dist;
			//std::cout<< "Temp Point 2: " << temp_point[2] << "Temp: " << temp <<  "Dist: " << dist << std::endl;
			if(dist>=-1.0) {
				this->m_boundaryConditions.push_back( i );
				//std::cout << "Boundary Condition: " << i << std::endl;
				count+=1;
			}
		}else if(m_degree=2){
			//float temp = solution[0] + temp_point[0] * solution[1] + temp_point[1] * solution[2] +
			//			pow(temp_point[0],2) * solution[3] + pow(temp_point[1],2) * solution[4] + (temp_point[0]*temp_point[1]) * solution[5];
			// Distancia de un punto a un plano.
			dist = temp_point[2] - ( m_solution[0] + (temp_point[0]*m_solution[1]) + (temp_point[1]*m_solution[2]) +
						pow(temp_point[0],2)*m_solution[3] + pow(temp_point[1],2)*m_solution[4] + ((temp_point[0]*temp_point[1])*m_solution[5]) );

			if(-dist<min) min=-dist;
			//std::cout<< "Temp Point 2: " << temp_point[2] << "Temp: " << temp <<  "Dist: " << dist << std::endl;
			if(dist>=-1.0) {
				this->m_boundaryConditions.push_back( i );
				//std::cout << "Boundary Condition: " << i << std::endl;
				count+=1;
			}
		}
	}
//	std::cout << "Min: " << min << std::endl;
//	std::cout << "Count: " << count << std::endl;
}


#endif