#pragma once

#define _USE_MATH_DEFINES // for C++
#include <cmath>

#include <math.h>

#include "cleaningMesh.h"

// Itk
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageRegionIterator.h"
#include "itkRegionOfInterestImageFilter.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkFlatStructuringElement.h"
                                                
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"

#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkAffineTransform.h"
#include "itkResampleImageFilter.h"

#include "itkMaskImageFilter.h"

#include "itkMedianImageFilter.h"

// From itk to vtk
#include "itkImageToVTKImageFilter.h"

// Vtk
#include "vtkSmartPointer.h"
#include "vtkMarchingCubes.h"

#include "vtkPolyData.h"

#include "vtkAppendPolyData.h"

#include "vtkPoints.h"
#include "vtkCell.h"
#include "vtkIdList.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataArray.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"

#include "vtkSmoothPolyDataFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"

#include "vtkPolyDataConnectivityFilter.h"
#include "vtkPolyDataWriter.h"
#include "vtkCleanPolyData.h"

#include "vtkCellArray.h"

#include "vtkPLYWriter.h"
// model Writer!
//#include "xmlModelWriter.h"

// typedef itk
typedef itk::Image<float,3> ImageSegmentation3floatType;
typedef itk::ImageFileReader<ImageSegmentation3floatType> ReaderSegmentationType;
typedef itk::ImageFileWriter<ImageSegmentation3floatType> WriterSegmentationType;

typedef itk::BinaryThresholdImageFilter<ImageSegmentation3floatType,ImageSegmentation3floatType> ThresholdSegmentationType;
typedef itk::RegionOfInterestImageFilter<ImageSegmentation3floatType, ImageSegmentation3floatType> RoiSegmentationType;

typedef itk::FlatStructuringElement<3> FlatType;
typedef itk::BinaryBallStructuringElement<float,3> KernelType;

typedef itk::BinaryMorphologicalClosingImageFilter<ImageSegmentation3floatType, ImageSegmentation3floatType,KernelType> Closing3floatType;
typedef itk::BinaryMorphologicalClosingImageFilter<ImageSegmentation3floatType, ImageSegmentation3floatType,FlatType> ClosingFlatType;

typedef itk::BinaryErodeImageFilter<ImageSegmentation3floatType, ImageSegmentation3floatType,KernelType> ErodeType;
typedef itk::BinaryErodeImageFilter<ImageSegmentation3floatType, ImageSegmentation3floatType,FlatType> ErodeFlatType;

typedef itk::NearestNeighborInterpolateImageFunction<ImageSegmentation3floatType, double> InterpolatorSegmentationType;
typedef itk::AffineTransform<double,3> AffTransType;
typedef itk::ResampleImageFilter<ImageSegmentation3floatType, ImageSegmentation3floatType> ResampleSegmentationType;

typedef itk::MaskImageFilter<ImageSegmentation3floatType, ImageSegmentation3floatType, ImageSegmentation3floatType> MaskSegmentationType;

typedef itk::BinaryDilateImageFilter<ImageSegmentation3floatType,ImageSegmentation3floatType,KernelType> DilateSegmentationType;
typedef itk::BinaryDilateImageFilter<ImageSegmentation3floatType,ImageSegmentation3floatType,FlatType> DilateFlatType;

typedef itk::MedianImageFilter<ImageSegmentation3floatType,ImageSegmentation3floatType> MedianFilterType;

// typedef from itk to vtk
typedef itk::ImageToVTKImageFilter<ImageSegmentation3floatType> ConnectTissueType;

class changeSegmentation
{
public:
	changeSegmentation(void);
	~changeSegmentation(void);

/* */
// public:
	
private:

	const char* side;
	int iso_vox_size;

	ImageSegmentation3floatType::Pointer m_newSegmentation;
	ImageSegmentation3floatType::Pointer m_maskSegmentation;
	ImageSegmentation3floatType::Pointer m_petiteRoi;
	
	vtkSmartPointer<vtkPolyData> m_GlandularMesh;
	
	float* m_solution;
	int m_degree;


/* */
public:
	// Set new segmentation:
	void SetNewSegmentation( std::string newSegmentation_filename);
	void SetSide( const char* side){this->side = side;};
	void SetFinalSpacing( int vox_size ){ this->iso_vox_size = vox_size;};
	void SetMaskSegmentation( std::string maskSegmentation_filename);
	
	void SetSolution(float* x_sol ){this->m_solution =x_sol;};
	void SetDegree(int degree){this->m_degree = degree;};

	// joining surfaces & tissue meshes
	void JoinMeshes( vtkSmartPointer<vtkPolyData> surface_mesh, vtkSmartPointer<vtkPolyData> & total_mesh);

	// Get new segmentation,
	ImageSegmentation3floatType::Pointer GetImagen() {return this->m_newSegmentation;};

private:
	void getROI(ImageSegmentation3floatType::Pointer imagen, const char* side, ImageSegmentation3floatType::Pointer & temp_imagen);
	void ExtractGlandularTissue();

	void approxTissue();

};

