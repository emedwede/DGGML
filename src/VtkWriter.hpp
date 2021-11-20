#ifndef VTK_FILE_WRITER_HPP
#define VTK_FILE_WRITER_HPP

#include <iostream>
#include <string>

#include "FileWriter.hpp"

#include "vtu11/vtu11.hpp"

namespace Cajete 
{

namespace VTK
{
    
struct MeshData
{
    std::vector<double> points;
    std::vector<vtu11::VtkIndexType> connectivity;
    std::vector<vtu11::VtkIndexType> offsets;
    std::vector<vtu11::VtkCellType> types;

    vtu11::Vtu11UnstructuredMesh mesh( )
    {
        return { points, connectivity, offsets, types };
    }
};

}
//Overides the file_writer interface, but preserves the algorithm
template <typename DataType>
class VtkFileWriter : public FileWriter<DataType> {
    using data_type = typename FileWriter<DataType>::data_type;
    protected:
        void create_file() const override {
            std::cout << "File is created when written\n";
        }

        void open_file() const override {
            std::cout << "Current file writing does not support appending\n";
        }

        void format_data(data_type data) const override {
            std::cout << "Formating the data to be vtk compatibile\n";
        }

        void write_file(std::string name) const override {
            name += ".vtu";
            std::cout << "Writing data to " << name << "\n";
        }

        void close_file() const override {
            std::cout << "Closing the file\n";
        }

    private:
        std::string filename;

        VTK::MeshData meshData;

        std::vector<double> pointData;
        std::vector<double> cellData;
        std::vector<vtu11::DataSetInfo> dataSetInfo;
};

} //end namespace Cajete

#endif 
