file(MAKE_DIRECTORY "results")

configure_file("particle_settings.json" "particle_settings.json" COPYONLY)
include_directories(ParticleModel)
add_executable(particle_simulator particle_simulator.cpp ../../thirdparty/simdjson/simdjson.cpp ParticleModel/particleModel.h)

target_link_libraries(particle_simulator PRIVATE vtu11::vtu11)
target_link_libraries(particle_simulator PUBLIC SUNDIALS::arkode SUNDIALS::nvecserial)
