#ifndef VTK_FILE_WRITER_HPP
#define VTK_FILE_WRITER_HPP

#include <iostream>
#include <string>

#include "FileWriter.hpp"

namespace Cajete 
{

//Overides the file_writer interface, but preserves the algorithm
class VtkFileWriter : public FileWriter {
    protected:
        void create_file() const override {
            std::cout << "Creating the vtk file\n";
        }

        void open_file() const override {
            std::cout << "Opening the vtk file\n";
        }

        void format_data() const override {
            std::cout << "Formating the data to be vtk compatibile\n";
        }

        void write_file() const override {
            std::cout << "Writing data to file in serial\n";
        }

        void close_file() const override {
            std::cout << "Closing the file\n";
        }

    private:
        std::string file;
};

} //end namespace Cajete

#endif 
