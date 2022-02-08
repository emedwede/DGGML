configure_file("settings.json" "settings.json" COPYONLY)

#TODO: Upgrade testing workflow
add_executable(memory-manager-test tests_main.cpp memory-manager-test.cpp)
#add_executable(graph-test tests_main.cpp graph-test.cpp)
add_executable(factory-test tests_main.cpp factory-test.cpp  ../thirdparty/simdjson/simdjson.cpp)
add_executable(file-writer-test tests_main.cpp file-writer-test.cpp)
add_executable(grid-test tests_main.cpp grid-test.cpp)
add_executable(cplex-test tests_main.cpp cplex-test.cpp)
add_executable(plant-matcher-test tests_main.cpp plant-matcher-test.cpp)
add_executable(binning-and-sort-test tests_main.cpp binning-and-sort-test.cpp)
add_executable(math-test tests_main.cpp math-test.cpp)

target_link_libraries(file-writer-test PRIVATE vtu11::vtu11)
target_link_libraries(cplex-test PRIVATE vtu11::vtu11)
target_link_libraries(factory-test PRIVATE vtu11::vtu11)
target_link_libraries(plant-matcher-test PRIVATE vtu11::vtu11)
target_link_libraries(binning-and-sort-test PRIVATE vtu11::vtu11)
