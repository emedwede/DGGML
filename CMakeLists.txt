#always start with the minimum required version
cmake_minimum_required(VERSION 3.16)

#set the project name
project(CajeteProto VERSION 0.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#set the C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17) 

add_subdirectory(src)
add_subdirectory(tests)
