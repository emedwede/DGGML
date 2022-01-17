#always start with the minimum required version
cmake_minimum_required(VERSION 3.16)

#set the project name
project(CajeteProto VERSION 0.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#set the C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17) 

# Add the serial host mode by default
# Currently turning it of would break the code
option(ENABLE_SERIAL_HOST "Enable Serial Host Mode" ON)
if(ENABLE_SERIAL_HOST)
	add_definitions(-DCAJETE_ENABLE_HOST_SERIAL)
endif()

# Add debugging option
option(USE_DEBUG "Enter debug mode" ON)
if (USE_DEBUG)
	add_definitions(-DCAJETE_DEBUG)
endif()

include_directories(dependencies/YAGL/include)
include_directories(thirdparty/simdjson)
include_directories(src)
include_directories(src/PlantModel)
include_directories(src/Grids)
include_directories(src/CellComplexes)
include_directories(src/CellComplexes/Utils)
add_subdirectory(dependencies/YAGL)
add_subdirectory(thirdparty/vtu11)
add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(tools)
add_subdirectory(simulators)

