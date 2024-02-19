#pragma once

#ifndef __cleaningMesh_h
#define __cleaningMesh_h

// Include
#include <iostream>
#include <vector>
#include <list>

// Vtk
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"

#include "vtkPoints.h"
#include "vtkCell.h"
#include "vtkCellArray.h"

// -- CLASS CLEANING MESH -- 
class cleaningMesh
{
public:
	cleaningMesh(void);
	~cleaningMesh(void);

private:

public:
	
	void removeDuplicatePoints(std::vector<double> &points_vector, std::vector<int> &cell_vector);
	void removeDuplicatedCells(std::vector<int> &cell_vector);

	void performNewMesh( std::vector<double> &points_vector, std::vector<int> &cell_vector);

	void searchCount( int k, int a, int b, std::vector<int> cell_vector,  int &ni , int &nd);
	void foundNonManifold(std::vector<int> cell_vector, std::list<int> &borderlist);
	void clearMesh(std::vector<int> &cell_vector, std::list<int> borderlist, bool &repeat);
};

#endif