//#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_RUNNER //used for our custom main

#include <iostream>

#include "catch.hpp"


#include "DGGMLConfig.hpp"

int main(int argc, char *argv[]) {
    //report the version
   std::cout << "\n\n";
    std::cout << argv[0] << " Version " << DGGML_VERSION_MAJOR << "."
        << DGGML_VERSION_MINOR << " ... Welcome to DGGML, the Graph Grammar Simulator Prototype!\n";
    std::cout << "\n\nUsage: " << argv[0] << "\n\n";
 
   
   int result = Catch::Session().run(argc, argv);
   
   return result;
}

