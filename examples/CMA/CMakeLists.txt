configure_file("settings.json" "settings.json" COPYONLY)
configure_file("cross_compare.json" "cross_compare.json" COPYONLY)
include_directories(PlantModel)
add_executable(mt_dgg_simulator mt_dgg_simulator.cpp ../../thirdparty/simdjson/simdjson.cpp)
add_executable(configuration_file_generator configuration_file_generator.cpp ../../thirdparty/simdjson/simdjson.cpp)

target_link_libraries(mt_dgg_simulator PRIVATE vtu11::vtu11)
target_link_libraries(mt_dgg_simulator PUBLIC SUNDIALS::arkode SUNDIALS::nvecserial)
