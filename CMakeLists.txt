option (GRAPHICAL "Graphical interface. Default value ON. For setting OFF: -DGRAPHICAL=OFF" ON)
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(TPMS)



if(GRAPHICAL)
  find_package(VTK COMPONENTS

    #DesignTPMS.cpp
    CommonCore			#vtkNew.h, vtkVersion.h
    CommonDataModel		#vtkImageData.h
    FiltersCore			#vtkFlyingEdges3D.h, vtkMarchingCubes.h
    FiltersSources  #vtkDecimate
    FiltersExtraction 
    FiltersGeometry
    FiltersModeling 
    FiltersFlowPaths
    IOGeometry			#vtkSTLWriter.h
    #FiltersModeling	#vtkFillHolesFilter.h
    
    #Utils.h for GRAPHICAL
    RenderingCore		#vtkRenderer.h, vtkRenderWindow.h, vtkRenderWindowInteractor.h,
						        #vtkPolyDataMapper.h, vtkActor.h, vtkProperty.h
    CommonColor			#vtkNamedColors.h
	  InteractionStyle	#For moving the object in the OpenGL windows
	  RenderingOpenGL2	#For opening the window in OpenGL
  
  )
else()
  find_package(VTK COMPONENTS

    #DesignTPMS.cpp
    CommonCore			#vtkNew.h, vtkVersion.h
    CommonDataModel		#vtkImageData.h
    FiltersCore			#vtkFlyingEdges3D.h, vtkMarchingCubes.h
    IOGeometry			#vtkSTLWriter.h
    FiltersSources  #vtkDecimate
    FiltersExtraction 
    FiltersGeometry
    FiltersModeling 
    FiltersFlowPaths
  
  )
endif()


if (NOT VTK_FOUND)
  message(FATAL_ERROR "mainTPMS: Unable to find the VTK build folder.")
endif()

# Prevent a "command line is too long" failure in Windows.
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")

add_executable(mainTPMS mainTPMS.cpp Utils.cpp TPMSsolidGenerator.cpp TPMSsheetGenerator.cpp Tpms.cpp)

target_link_libraries(mainTPMS PRIVATE ${VTK_LIBRARIES}
)
# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS mainTPMS
  MODULES ${VTK_LIBRARIES}
)
