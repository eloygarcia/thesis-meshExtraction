#pragma once

// itk
#include "itkImage.h"
#include "itkImageFileReader.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"

#include "itkGDCMImageIO.h"
#include "itkMetaDataDictionary.h"

#include "itkRegionOfInterestImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageRegionIterator.h"

#include "itkAddImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkBinaryFillholeImageFilter.h"

#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkAffineTransform.h"

#include "itkBinaryMorphologicalClosingImageFilter.h"

#include "itkBinaryContourImageFilter.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkAndImageFilter.h"

#include "itkImageFileWriter.h"

// temporal vtk
#include "vtkVertexGlyphFilter.h"
#include "vtkPolyData.h"
#include "vtkIdList.h"
#include "vtkPolygon.h"
#include "vtkCellArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSetMapper.h"
#include "vtkTriangle.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkFeatureEdges.h"

#include "vtkSmoothPolyDataFilter.h"

#include "vtkSTLWriter.h"
#include "vtkUnstructuredGridWriter.h"
// Itk.

// From Itk to Vtk
#include "itkImageToVTKImageFilter.h"

// VTK
#include "vtkCell.h"

#include "vtkSmartPointer.h"
#include "vtkMarchingCubes.h"

#include "vtkSmoothPolyDataFilter.h"
#include "vtkCleanPolyData.h"
#include "vtkProperty.h"

#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
//#include "vtkPLYWriter.h"
//#include "vtkSTLWriter.h"

#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"

// typedef.
typedef itk::Image<unsigned short, 3> ImageMriType;
typedef itk::ImageFileReader<ImageMriType> ReaderMriType;

typedef itk::RegionOfInterestImageFilter<ImageMriType, ImageMriType> RoiMriType;
typedef itk::BinaryThresholdImageFilter<ImageMriType, ImageMriType> ThresholdMriType;

typedef itk::SubtractImageFilter<ImageMriType, ImageMriType> SubstractMriType;

typedef itk::ResampleImageFilter<ImageMriType, ImageMriType> ResampleMriType;
typedef itk::NearestNeighborInterpolateImageFunction<ImageMriType, double> InterpolatorMriType;
typedef itk::AffineTransform<double, 3> AffineTransformType;

typedef itk::BinaryContourImageFilter<ImageMriType,ImageMriType> ContourMriType;

typedef itk::BinaryBallStructuringElement< ImageMriType::PixelType, 3> BallType;
typedef itk::BinaryMorphologicalClosingImageFilter<ImageMriType, ImageMriType,BallType> ClosingMriType;

typedef itk::BinaryErodeImageFilter<ImageMriType, ImageMriType, BallType> ErodeMriType;
typedef itk::BinaryDilateImageFilter<ImageMriType, ImageMriType, BallType> DilateMriType;
typedef itk::AndImageFilter<ImageMriType> AndMriFilter;

typedef itk::ImageFileWriter<ImageMriType> Writer3UsType;

// typedef -- from itk to vtk
typedef itk::ImageToVTKImageFilter<ImageMriType> ConnectorType;

// typedef
typedef itk::GDCMImageIO  IOType;
typedef itk::MetaDataDictionary MetaDataType;

typedef itk::ShapeLabelObject< unsigned short, 3 > ShapeLabelObjectMriType;
typedef itk::LabelMap< ShapeLabelObjectMriType > LabelMapMriType;
typedef itk::LabelImageToShapeLabelMapFilter< ImageMriType, LabelMapMriType> I2LMriType;

// oriented 
#include "itkOrientImageFilter.h"
typedef itk::OrientImageFilter<ImageMriType,ImageMriType> OrientedType;

// statisttics 
#include "itkLabelStatisticsImageFilter.h"
typedef itk::LabelStatisticsImageFilter<ImageMriType, ImageMriType> LabelType;

class MRSurface
{
public:
	MRSurface(void);
	~MRSurface(void);

// variables
public:

	// breast mask
	ImageMriType::Pointer m_breastmask;
	bool has_breastmask;

	// albert segmentation
	ImageMriType::Pointer m_albertsegmentation;
	bool has_segmentation;

	ImageMriType::DirectionType direction;

public:

	// mesh
	vtkSmartPointer<vtkUnstructuredGrid> m_model;
	std::vector<int> m_boundaryConditions;

private:
	float iso_vox_size;
	int m_degree;
	char * m_side;
	changeSegmentation * new_segmentation;
	mesh * new_mesh;

	std::string breastmaskfilename;
	std::string albertsegmentationfilename;
	std::string probabilisticsegmentationfilename;

// metodos
public:
	
	// read data
	void read_breastmask(std::string inputfilename);
	void read_segmentation(std::string inputfilename);
	void read_probabilistic_segmentation( std::string inputfilename );

	void SetIsoVoxel( float iso_vox ){ this->iso_vox_size = iso_vox;};
	void SetSideOfinterest( const char* side ){this->m_side = (char*) side;};

	// meshing
	void getmodel(vtkSmartPointer<vtkPolyData> &surf_mesh);
	
	float* solution;
	//void write_niftySim_model();

	void set_approxDegree(int degree) {this->m_degree = degree;};

protected:
	// aux. functions
	void getRegionOfInterest(ImageMriType::Pointer image,const char* side, ImageMriType::Pointer &new_image);
	void resampleFunction(ImageMriType::Pointer image, float iso_vox_size, ImageMriType::Pointer &outputimage);
	void closingFunction (ImageMriType::Pointer image, BallType ball, ImageMriType::Pointer &outputImage, int r);
	void visualizationFunction (ImageMriType::Pointer image,  vtkSmartPointer<vtkMarchingCubes> & marchingOutput); //, vtkSmartPointer<vtkActor> &actorOutput);

	void Gaussian_solver(std::vector<double> points_interface, float* &solution, int degree);
	void controlFunction(std::vector<double> &points_interface, std::vector<int> &connect_interface,
		std::vector<double> points_surface, std::vector<int> connect_surface, float* solution, bool &eval, int degree);

	void upAndDownCut(ImageMriType::Pointer image, ImageMriType::Pointer mask, ImageMriType::Pointer &outputimage);

public:
	void getBreastSurfaceModel( const char* side, vtkSmartPointer<vtkPolyData> &surf_mesh);
};


