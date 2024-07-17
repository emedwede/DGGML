#always start with the minimum required version
cmake_minimum_required(VERSION 3.16)

#set the project name
project(DGGML VERSION 1.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#set the C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17) 

#add_compile_options(-O3)

# Add the serial host mode by default
# Currently turning it of would break the code
option(ENABLE_SERIAL_HOST "Enable Serial Host Mode" ON)
if(ENABLE_SERIAL_HOST)
	add_definitions(-DDGGML_ENABLE_HOST_SERIAL)
endif()

# Add debugging option
option(USE_DEBUG "Enter debug mode" ON)
if (USE_DEBUG)
	#add_compile_definitions(DGGML_DEBUG)
endif()

# In the build script simpy do this:
#
# Set the variable SUNDIALS_DIR to the SUNDIALS install directory.
# When using the cmake CLI command, this can be done like so:
#   cmake -D SUNDIALS_DIR=/path/to/sundials/installation

#try to find the sundials package for ODE solving
find_package(SUNDIALS 6.7.0 QUIET)

if(NOT SUNDIALS_FOUND)
	message(STATUS "SUNDIALS 6.7.0 not found, fetching with FetchContent...")
	set(FETCHCONTENT_QUIET FALSE)
	# Include FetchContent module
	include(FetchContent)

	message("Downloading SUNDIALS, this may take a few minutes...")
	# Specify Sundials version and URL
	FetchContent_Declare(
			sundials
			GIT_REPOSITORY https://github.com/LLNL/sundials.git
			GIT_TAG v6.7.0
			GIT_PROGRESS TRUE
	)

	# Fetch and make available
	FetchContent_MakeAvailable(sundials)
endif()

# TODO: see if there is a better way than include_directories command
include_directories(dependencies/YAGL/include)
include_directories(thirdparty/simdjson)
include_directories(src)
include_directories(src/Grids)
include_directories(src/Utlities)
include_directories(src/Grammar)
include_directories(src/CellComplexes)
include_directories(src/CellComplexes/Utils)
include_directories(src/Solvers)
add_subdirectory(dependencies/YAGL)
add_subdirectory(thirdparty/vtu11)
add_subdirectory(src)
#add_subdirectory(tests)
add_subdirectory(tools)
add_subdirectory(examples)

