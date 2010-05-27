#Define the project we are working on
project(mscomplex-2d)

# Look for various software we need
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(QT_MIN_VERSION "4.5.0")
find_package(Qt4 REQUIRED)

if(NOT QT4_FOUND)
  message(FATAL_ERROR "Failed to find Qt 4.5 or greater.")
endif(NOT QT4_FOUND)

set(QT_USE_QTGUI TRUE)
set(QT_USE_QTOPENGL TRUE)
set(QT_USE_QTXML TRUE)

include(${QT_USE_FILE})

find_package(QGLViewer)

if(NOT QGLVIEWER_FOUND)
  message(FATAL_ERROR "Failed to find QGLViewer .. quitting")
endif(NOT QGLVIEWER_FOUND)

find_package(Boost COMPONENTS program_options thread serialization)

if(NOT Boost_FOUND)
  message(FATAL_ERROR "Cannot find boost")
endif(NOT Boost_FOUND)

find_package(OpenCL)

if(NOT OPENCL_FOUND)
  message(FATAL_ERROR "could not find openCL ..")
endif(NOT OPENCL_FOUND)

add_subdirectory(utls)

set(dmsc_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/main.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid.h

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/prefix_scan.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/prefix_scan.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_dataset.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_dataset.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_mscomplex.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_mscomplex.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_datamanager.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_datamanager.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/bitonic_sort.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/bitonic_sort.cpp
  )

QT4_ADD_RESOURCES(dmsc_RC_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/ocl/mscomplex.qrc
  )

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

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/
  ${CMAKE_CURRENT_SOURCE_DIR}/viewer/
  ${utls_SOURCE_DIR}/include
  ${CMAKE_CURRENT_BINARY_DIR}
  ${OPENCL_INCLUDE_DIRS}
  ${QT_ADDITIONAL_INCLUDE_PATHS}
  ${QGLVIEWER_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
  )

add_executable(${PROJECT_NAME}
  ${dmsc_SRCS}
  ${dmsc_RC_SRCS}
  ${viewer_SRCS}
  ${viewer_MOC_SRCS}
  ${viewer_UI_SRCS}
  )

target_link_libraries(
  ${PROJECT_NAME}
  utls
  ${Boost_LIBRARIES}
  ${OPENCL_LIBRARIES}
  ${QT_LIBRARIES}
  ${QT_ADDITIONAL_LIBRARIES}
  ${QGLVIEWER_LIBRARY}
)
