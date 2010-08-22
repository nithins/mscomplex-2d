#Define the project we are working on
project(mscomplex2d)

# Look for various software we need
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(Boost COMPONENTS program_options thread serialization REQUIRED)

find_package(OpenCL REQUIRED)

set(dmsc_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_dataset.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_dataset.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_mscomplex.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_mscomplex_ensure.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_mscomplex.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_datamanager.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_datamanager.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/bitonic_sort.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/bitonic_sort.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/prefix_scan.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/prefix_scan.cpp
  )

option(NO_GUI "build without qt gui" OFF)

if(NOT NO_GUI)
  SET(NO_OPENGL OFF CACHE bool "no opengl" FORCE)
else(NOT NO_GUI)
  SET(NO_OPENGL ON CACHE bool "no opengl" FORCE)
endif(NOT NO_GUI)

add_subdirectory(utls)

find_package(PythonInterp REQUIRED)

set(dmsc_ocl_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/ocl/assigngradient.cl
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/ocl/bfs_watershed.cl
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/ocl/bitonic_sort.cl
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/ocl/collate_critpts.cl
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/ocl/common_funcs.cl
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/ocl/scan_kernel.cl
  )

set(dmsc_ocl_SRCS
  ${CMAKE_CURRENT_BINARY_DIR}/dmsc/oclsources.h
  ${CMAKE_CURRENT_BINARY_DIR}/dmsc/oclsources.cpp
  )

add_custom_command(
  OUTPUT ${dmsc_ocl_SRCS}
  COMMAND ${PYTHON_EXECUTABLE} ${utls_SOURCE_DIR}/copy_files_to_c_header.py
  ${CMAKE_CURRENT_BINARY_DIR}/dmsc/oclsources ${dmsc_ocl_FILES}
  DEPENDS ${dmsc_ocl_FILES} ${utls_SOURCE_DIR}/copy_files_to_c_header.py
  )

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}/dmsc/
  ${utls_SOURCE_DIR}/include
  ${OPENCL_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  )

set(${PROJECT_NAME}_SRCS
  ${${PROJECT_NAME}_SRCS}
  ${dmsc_SRCS}
  ${dmsc_ocl_SRCS}
  ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
  )

set(${PROJECT_NAME}_LINKLIBS
  ${${PROJECT_NAME}_LINKLIBS}
  utls
  ${Boost_LIBRARIES}
  ${OPENCL_LIBRARIES}
  )


if(NOT NO_GUI)
  set(QT_MIN_VERSION "4.5.0")
  set(QT_USE_QTGUI TRUE)
  set(QT_USE_QTOPENGL TRUE)
  set(QT_USE_QTXML TRUE)

  find_package(Qt4 REQUIRED)

  include(${QT_USE_FILE})

  find_package(QGLViewer REQUIRED)

  SET(QXT_USE_QXTCORE TRUE)
  SET(QXT_USE_QXTGUI TRUE)

  find_package(Qxt REQUIRED)

  set(viewer_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/viewer/grid_viewer.h
    ${CMAKE_CURRENT_SOURCE_DIR}/viewer/grid_viewer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/viewer/grid_viewer_mainwindow.h
    ${CMAKE_CURRENT_SOURCE_DIR}/viewer/grid_viewer_mainwindow.cpp
    )

  QT4_WRAP_CPP(viewer_MOC_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/viewer/grid_viewer_mainwindow.h
    )

  set(viewer_UIS
    ${CMAKE_CURRENT_SOURCE_DIR}/viewer/grid_viewer_mainwindow.ui
    )

  QT4_WRAP_UI(viewer_UI_SRCS ${viewer_UIS})

  set(viewer_shader_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/viewer/shaders/grad_vert.glsl
    ${CMAKE_CURRENT_SOURCE_DIR}/viewer/shaders/grad_geom.glsl
    ${CMAKE_CURRENT_SOURCE_DIR}/viewer/shaders/sphere_vert.glsl
    ${CMAKE_CURRENT_SOURCE_DIR}/viewer/shaders/sphere_geom.glsl
    ${CMAKE_CURRENT_SOURCE_DIR}/viewer/shaders/sphere_frag.glsl
    ${CMAKE_CURRENT_SOURCE_DIR}/viewer/shaders/cylinder_vert.glsl
    ${CMAKE_CURRENT_SOURCE_DIR}/viewer/shaders/cylinder_geom.glsl
    ${CMAKE_CURRENT_SOURCE_DIR}/viewer/shaders/cylinder_frag.glsl
    )

  set(viewer_shader_SRCS
    ${CMAKE_CURRENT_BINARY_DIR}/viewer/shadersources.h
    ${CMAKE_CURRENT_BINARY_DIR}/viewer/shadersources.cpp
    )

  add_custom_command(
    OUTPUT ${viewer_shader_SRCS}
    COMMAND ${PYTHON_EXECUTABLE} ${utls_SOURCE_DIR}/copy_files_to_c_header.py
    ${CMAKE_CURRENT_BINARY_DIR}/viewer/shadersources ${viewer_shader_FILES}
    DEPENDS ${viewer_shader_FILES} ${utls_SOURCE_DIR}/copy_files_to_c_header.py
    )

  include_directories(
    ${QT_ADDITIONAL_INCLUDE_PATHS}
    ${QGLVIEWER_INCLUDE_DIR}
    ${QXT_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/viewer/
    ${CMAKE_CURRENT_BINARY_DIR}/viewer/
    )

  set(${PROJECT_NAME}_SRCS
    ${${PROJECT_NAME}_SRCS}
    ${viewer_SRCS}
    ${viewer_MOC_SRCS}
    ${viewer_UI_SRCS}
    ${viewer_shader_SRCS}
  )

  set(${PROJECT_NAME}_LINKLIBS
    ${${PROJECT_NAME}_LINKLIBS}
    ${QT_LIBRARIES}
    ${QT_ADDITIONAL_LIBRARIES}
    ${QGLVIEWER_LIBRARY}
    ${QXT_LIBRARIES}
    )
endif(NOT NO_GUI)

configure_file(${PROJECT_SOURCE_DIR}/config.h.in ${PROJECT_BINARY_DIR}/config.h)

add_executable(${PROJECT_NAME} ${${PROJECT_NAME}_SRCS})

target_link_libraries(${PROJECT_NAME} ${${PROJECT_NAME}_LINKLIBS})
