#add the configuration file, can be kept up to date between commits
configure_file(CajeteConfig.hpp.in CajeteConfig.hpp)

#set(SOURCE_FILES cell_list_simple.cpp)
#set(EXEC cell_list_simple)
file(GLOB HEADER_FILES "*.hpp")
file(GLOB SOURCE_FILES "*.cpp")

#set(LIB_NAME cajeteproto)
#add_library(cajeteproto ${SOURCE_FILES} ${HEADER_FILES})
