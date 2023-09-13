cmake_minimum_required(VERSION 3.5)

project(guarding)



set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(Qt5_core "/opt/Qt5.7.0/5.7/gcc_64/lib/cmake/Qt5_core")
set(Qt5_DIR "/opt/Qt5.7.0/5.7/gcc_64/lib/cmake/Qt5")
set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(include)


# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )
if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()


# Boost and its components
find_package( Boost REQUIRED )
if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# Qt5 and its components
find_package(Qt5 COMPONENTS Widgets REQUIRED)
if ( NOT Qt5_FOUND )

  message(STATUS "This project requires Qt5, and will not be compiled.")

  return()  

endif()

set_target_properties(Qt5::Widgets PROPERTIES MAP_IMPORTED_CONFIG_COVERAGE "RELEASE")

    qt5_wrap_ui( DT_UI_FILES mainwindow.ui)
    
    add_executable(guarding
        src/glwidget.cpp
        src/bezier.cpp
        src/helper.cpp
        src/main.cpp
        src/mainwindow.cpp
        src/refinement.cpp
        src/variational.cpp
        
        include/bezier.h
        include/glwidget.h
        include/helper.h
        include/mainwindow.h
        include/refinement.h
        include/variational.h
        mainwindow.ui

        ${DT_UI_FILES}
    )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS guarding )

target_link_libraries(guarding PRIVATE CGAL::CGAL Qt5::Widgets)
