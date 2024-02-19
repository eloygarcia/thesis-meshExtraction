#pragma once

#ifndef __mesh_h
#define __mesh_h

#define _USE_MATH_DEFINES // for C++
#include <cmath>

#include "cleaningMesh.h"
#include "tetgen.h"

#include <math.h>

// Vtk
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCell.h"
#include "vtkIdList.h"

#include "vtkUnstructuredGrid.h"
#include "vtkTetra.h"

#include "vtkDataArray.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"

#include "vtkCellData.h"

#include "vtkSmoothPolyDataFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"

#include "vtkImageGridSource.h"
#include "vtkProbeFilter.h"

#include "vtkPolyDataConnectivityFilter.h"
#include "vtkPolyDataWriter.h"
#include "vtkCleanPolyData.h"


class mesh
{
public:
	mesh(void);
	~mesh(void);

public:
	// PolyData mesh.
	vtkSmartPointer<vtkPolyData> m_surface_mesh;
	std::vector<double> m_points_surface;
	std::vector<int> m_cell_vector;
	std::list<int> m_cell_unique_list;

	// Points & elements !
	vtkSmartPointer<vtkUnstructuredGrid> m_tetrahedral_mesh;

	std::vector<double> m_initial_points;
	std::vector<double> m_final_points;
	std::vector<int> m_elements;
	std::list<int> m_cell_element_list;
	
	std::vector<int> m_tetrahedron_tissue;
	std::vector<double> m_intermedial_compression_points;

	// Boundary Conditions !!
	std::vector<int> m_boundaryConditions;
	int m_degree;
	float * m_solution;

public:
	// load mesh
	void set_surfacemesh(vtkSmartPointer<vtkPolyData> surface_mesh);
	void set_boundaryConditions(float* solution){m_solution = solution; };

	vtkSmartPointer<vtkUnstructuredGrid>  get_tetrahedralMesh() {return this->m_tetrahedral_mesh;};
	void set_degree( int degree){ m_degree = degree; };
	void get_boundaryConditions();

public:
	// aux. functions
	void improve_quality(vtkSmartPointer<vtkPolyData> &surface_mesh);
	void tetgen_ext(vtkSmartPointer<vtkPolyData> surface_mesh);

	void createUnstructuredGrid( std::vector<double> point_vector,
									std::vector<int> cells_vector,	
									std::vector<int> bC_list,
									std::vector<int> tissue_list,
									vtkSmartPointer<vtkUnstructuredGrid> &grid);

};

#include "mesh.cpp"

#endif

