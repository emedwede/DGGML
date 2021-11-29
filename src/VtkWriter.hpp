#ifndef VTK_FILE_WRITER_HPP
#define VTK_FILE_WRITER_HPP

#include <iostream>
#include <string>

#include "FileWriter.hpp"

#include "vtu11/vtu11.hpp"

#include <vector>
#include <unordered_map>

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

        void format_data(data_type data) override {
            std::cout << "Formating the data to be vtk compatibile\n";
            std::cout << data;
            
            //first rekey the points for vtk format
            vtu11::VtkIndexType count = 0;
            for(auto i = data.node_list_begin(); i != data.node_list_end(); i++)
            {
                rekey.insert({i->first, count});
                count++;
            }
            for(auto i = data.node_list_begin(); i != data.node_list_end(); i++) 
            {
                points.push_back(i->second.getData());
                points.push_back(0.0);
                points.push_back(0.0);
                keyData.push_back(i->first);
                types.push_back(4); // vtk object type
                if(offsets.size() == 0)
                    offsets.push_back(data.out_neighbors(i->second).size());
                else
                {
                    auto last = offsets.back();
                    auto size = data.out_neighbors(i->second).size();
                    offsets.push_back(last+size);
                }
                for(auto j = data.out_neighbors_begin(i->second); j != data.out_neighbors_end(i->second); j++) 
                {
                    connectivity.push_back(rekey[*j]);
                }
            }
        }

        void write_file(std::string name) override {
            name += ".vtu";
            std::cout << "Writing data to " << name << "\n";
            
            std::vector<vtu11::DataSetInfo> dataSetInfo
            {
                {"Unique ID", vtu11::DataSetType::PointData, 1}
            };
            vtu11::Vtu11UnstructuredMesh mesh {points, connectivity, offsets, types};
            vtu11::writeVtu(name, mesh, dataSetInfo, {keyData}, "Ascii");
        }

        void close_file() override {
            std::cout << "Closing the file\n";
            rekey.clear();
            points.clear();
            keyData.clear();
            types.clear();
            offsets.clear();
            connectivity.clear();
            cellData.clear();
        }

    private:
        std::string filename;

        VTK::MeshData meshData;
        
        std::unordered_map<typename DataType::key_type, vtu11::VtkIndexType> rekey;        
        std::vector<double> points;
        std::vector<double> keyData;
        std::vector<vtu11::VtkCellType> types;
        std::vector<vtu11::VtkIndexType> offsets;
        std::vector<vtu11::VtkIndexType> connectivity;
        std::vector<double> cellData;
        //std::vector<vtu11::DataSetInfo> dataSetInfo;
};

} //end namespace Cajete

#endif 
