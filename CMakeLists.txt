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
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/src/main.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/include/discreteMorseDS.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/include/discreteMorseAlgorithm.h

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/include/rectangle_complex.h

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/include/prefix_scan.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/src/prefix_scan.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/include/grid_dataset.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/src/grid_dataset.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/include/grid_mscomplex.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/src/grid_mscomplex.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/include/grid_datamanager.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/src/grid_datamanager.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/include/grid_viewer.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/src/grid_viewer.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/include/grid_viewer_mainwindow.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/src/grid_viewer_mainwindow.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/include/bitonic_sort.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/src/bitonic_sort.cpp
  )

QT4_WRAP_CPP(dmsc_MOC_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/include/grid_viewer_mainwindow.h
  )

set(dmsc_UIS
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/forms/grid_viewer_mainwindow.ui
  )

QT4_ADD_RESOURCES(dmsc_RC_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/resources/mscomplex.qrc
  )

QT4_WRAP_UI(dmsc_UI_SRCS ${dmsc_UIS})

include_directories(
  ${OPENCL_INCLUDE_DIRS}
  ${QT_ADDITIONAL_INCLUDE_PATHS}
  ${QGLVIEWER_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
  ${utls_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/include
  ${CMAKE_CURRENT_BINARY_DIR}
  )

add_executable(${PROJECT_NAME}
  ${dmsc_SRCS}
  ${dmsc_MOC_SRCS}
  ${dmsc_UI_SRCS}
  ${dmsc_RC_SRCS}
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
