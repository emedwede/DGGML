#TODO: Upgrade testing workflow
add_executable(memory-manager-test tests_main.cpp memory-manager-test.cpp)
#add_executable(graph-test tests_main.cpp graph-test.cpp)
add_executable(factory-test tests_main.cpp factory-test.cpp)
add_executable(file-writer-test tests_main.cpp file-writer-test.cpp)

target_link_libraries(file-writer-test PRIVATE vtu11::vtu11)
target_link_libraries(factory-test PRIVATE vtu11::vtu11)
