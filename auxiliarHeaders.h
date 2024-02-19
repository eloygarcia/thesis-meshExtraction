#ifndef __auxiliarHeaders_h
#define __auxiliarHeaders_h

// #include "vld.h"

#define _USE_MATH_DEFINES// for C++

#include <cmath>
#include <math.h>
#include <iostream>

#include <string>
#include <string.h>

#include <vector>

// classes
#include "mammogram.h"
#include "MRimage.h"

// #include "IntensityBased.h"

// #include "HillClimbingOpt.h"
// #include "StochasticHillClimbing.h"
// #include "SimulatedAnnealing.h"

// #include "auxiliarDefinitions.h"

// #include "MechanicalProperties.h"

// Projections !
// #include "projection.h"
// #include "classicProjection.h"

// #include "NewProjection.h"

// #ifdef _GPU_
//	 #include "CudaProjection.h"
// #endif


// #include "changeSegmentation.h"

// #include "xmlModelWriter.h"

// #include "auxiliarDefinitions.h" // Se encuentra en las classes de registro:
									//			IntensityBased.h
									//			FeaturedBased.h
/*
#include "auxiliarFunctions_Transformations.h"

// Itk
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkVector.h"
#include "itkImageFileWriter.h"
*/

// Vtk
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkDataSetMapper.h"

#include "vtkPointData.h"

#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkPLYWriter.h"

#include "vtkPlaneSource.h"
#include "vtkSphereSource.h"
#include "vtkLineSource.h"

/*
// typdef
typedef itk::Image<float,3> Image3floatType;
typedef itk::Image<float,2> Image2floatType;
typedef itk::ImageFileReader<Image3floatType> Reader3floatType;

typedef itk::Vector<float,3> vectorPixelType;
typedef itk::Image<vectorPixelType,3> ImageVector3Type;
typedef itk::ImageFileWriter<ImageVector3Type> WriterVector3Type;

typedef itk::ImageFileWriter<Image3floatType> Writer3floatType;
typedef itk::ImageFileWriter<Image2floatType> Writer2floatType;
typedef itk::ImageFileWriter<ImageMriType> WriterUShortType;

// Image metric:
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkTranslationTransform.h"

typedef itk::MeanSquaresImageToImageMetric<Image2floatType, Image2floatType> MSQmetricType;
typedef itk::LinearInterpolateImageFunction<Image2floatType, double> InterpolatorType;
typedef itk::TranslationTransform<double,2> TranslationTransformType;

typedef itk::ResampleImageFilter<ImageVector3Type, ImageVector3Type> ResampleVectorType;
typedef itk::LinearInterpolateImageFunction<ImageVector3Type, double> LinearInterpolateVectorType;

// HACE FALTA ALGO DE ESTO DE AQUI EN ADELANTE?? //

// DRR
#include "itkResampleImageFilter.h"
#include "itkTranslationTransform.h"
#include "itkRayCastInterpolateImageFunction.h"
typedef itk::ResampleImageFilter< Image3floatType, Image3floatType> FilterType;
typedef itk::TranslationTransform<double,2> Transform2Type;
typedef itk::TranslationTransform<double,3> Transform3Type;
typedef itk::RayCastInterpolateImageFunction < Image3floatType, double > RayCastType;
typedef itk::BSplineInterpolateImageFunction<Image3floatType> SplineInterpolatorType;
typedef itk::BSplineInterpolateImageFunction<Image2floatType> spline2dInterpolatorType;
typedef itk::ExtractImageFilter<Image3floatType, Image2floatType> ExtractFilterType;

/*
// Dice Index
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkLabelOverlapMeasuresImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"

//typedef itk::Image<unsigned short, 2> UShort2DImage;
//typedef itk::BinaryThresholdImageFilter<Image2floatType, UShort2DImage> ThresholdDCType;
typedef itk::Image<unsigned short,2> LabelAType;
typedef itk::BinaryThresholdImageFilter<Image2floatType,LabelAType> ThresholdDCType;
typedef itk::BinaryImageToLabelMapFilter<LabelAType> ImageToLabelType;
typedef itk::LabelOverlapMeasuresImageFilter<LabelAType> OverlapMeasuresType;
typedef itk::BinaryBallStructuringElement<unsigned short, 2> StructuringElementType;  
typedef itk::BinaryMorphologicalClosingImageFilter<LabelAType, LabelAType, StructuringElementType> BinaryMorphologicalClosingImageFilterType; // Esto esta´a mal
*/

/*
// include para cuda proyecction

#include "itkImportImageFilter.h"
typedef itk::ImportImageFilter<float, 3> ImportFloatImageFilterType;

// Deformación de la imagen 
#include "itkBSplineInterpolateImageFunction.h"
#include "itkWarpImageFilter.h"

typedef itk::BSplineInterpolateImageFunction<Image3floatType> splineInterpolationMeshType;
typedef itk::WarpImageFilter<Image3floatType, Image3floatType, ImageVector3Type> warpMeshType;


*/
#endif