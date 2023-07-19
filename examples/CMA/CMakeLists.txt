file(MAKE_DIRECTORY "results")

configure_file("settings.json" "settings.json" COPYONLY)

add_executable(mt_dgg_simulator mt_dgg_simulator.cpp ../../thirdparty/simdjson/simdjson.cpp)

target_link_libraries(mt_dgg_simulator PRIVATE vtu11::vtu11)
target_link_libraries(mt_dgg_simulator PUBLIC SUNDIALS::arkode SUNDIALS::nvecserial)
