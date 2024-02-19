cmake_minimum_required(VERSION 2.8)
project(GettingMeshes)

find_package(ITK)
include(${ITK_USE_FILE})

find_package(ItkVtkGlue)
include(${ItkVtkGlue_USE_FILE})

find_package(VTK)
include(${VTK_USE_FILE})

add_library(	xRayMammograms
					mammogram.cpp	)

add_library(	MagneticResonanceImages
					MRimage.cpp
					changeSegmentation.cpp
					mesh.cpp
					cleaningMesh.cpp
					tetgen.cxx
					predicates.cxx	)

add_executable(	GettingMeshes
						Source.cpp	)

target_link_libraries(	GettingMeshes
						${ITK_LIBRARIES} 
						${VTK_LIBRARIES} 
						${Glue}
						xRayMammograms
						MagneticResonanceImages
						)