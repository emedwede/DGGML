#ifndef VTK_FILE_WRITER_HPP
#define VTK_FILE_WRITER_HPP

#include <iostream>
#include <string>

#include "FileWriter.hpp"

#include "CartesianGrid2D.hpp"

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
                points.push_back(i->second.getData().position[0]);
                points.push_back(i->second.getData().position[1]);
                points.push_back(i->second.getData().position[2]);

                comp_type.push_back(i->second.getData().type);

                keyData.push_back(i->first);
                //types.push_back(4); // vtk object type
                //if(offsets.size() == 0)
                //    offsets.push_back(0); //offsets.push_back(data.out_neighbors(i->second).size());
                //else
                //{
                    //auto last = offsets.back();
                    //auto size = data.out_neighbors(i->second).size();
                    //offsets.push_back(last+size);
                //}
                for(auto j = data.out_neighbors_begin(i->second); j != data.out_neighbors_end(i->second); j++) 
                {
                    std::size_t last;
                    if( offsets.size() == 0 )
                        last = 0;
                    else 
                        last = offsets.back();

                    types.push_back(4);
                    offsets.push_back(last+2);
                    connectivity.push_back(rekey[i->first]);
                    connectivity.push_back(rekey[*j]);
                }
            }
        }

        void write_file(std::string name) override {
            name += ".vtu";
            std::cout << "Writing data to " << name << "\n";
            
            std::vector<vtu11::DataSetInfo> dataSetInfo
            {
                {"Unique ID", vtu11::DataSetType::PointData, 1},
                {"Node Type", vtu11::DataSetType::PointData, 1}   
            };
            vtu11::Vtu11UnstructuredMesh mesh {points, connectivity, offsets, types};
            vtu11::writeVtu(name, mesh, dataSetInfo, {keyData, comp_type}, "Ascii");
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
            comp_type.clear();
        }

    private:
        std::string filename;

        VTK::MeshData meshData;
        
        std::unordered_map<typename DataType::key_type, vtu11::VtkIndexType> rekey;        
        std::vector<double> points;
        std::vector<double> keyData;
        std::vector<vtu11::VtkCellType> types;
        std::vector<double> comp_type;
        std::vector<vtu11::VtkIndexType> offsets;
        std::vector<vtu11::VtkIndexType> connectivity;
        std::vector<double> cellData;
        //std::vector<vtu11::DataSetInfo> dataSetInfo;
};

//Overides the file_writer interface, but preserves the algorithm
class GridFileWriter : public FileWriter<std::pair<CartesianGrid2D, std::vector<int>>> {
    using data_type = typename FileWriter<std::pair<CartesianGrid2D, std::vector<int>>>::data_type;
    protected:
        void create_file() const override {
            std::cout << "File is created when written\n";
        }

        void open_file() const override {
            std::cout << "Current file writing does not support appending\n";
        }

        void format_data(data_type _data) override {
            std::cout << "Formating the data to be vtk compatibile\n";
            auto data = _data.first;
            auto labels = _data.second;
            auto n = data._nx;
            auto m = data._ny;
            auto dx = data._dx;
            auto dy = data._dy;
            auto px = data._px;
            auto py = data._py;
 
            for(auto i = 0; i < data.totalNumPoints(); i++)
            {
                double xp, yp, zp = 0.0;
                data.cardinalLatticeToPoint(xp, yp, i);
                points.push_back(xp);
                points.push_back(yp);
                points.push_back(zp);    
            }

            for(auto i = 0; i < n; i++)
            {
                for(auto j = 0; j < m; j++)
                {
                    std::size_t cardinal;
                    cardinal = data.cardinalLatticeIndex(i, j);
                    connectivity.push_back(cardinal);
                    cardinal = data.cardinalLatticeIndex(i+1, j);
                    connectivity.push_back(cardinal);
                    cardinal = data.cardinalLatticeIndex(i+1, j+1);
                    connectivity.push_back(cardinal);
                    cardinal = data.cardinalLatticeIndex(i, j+1);
                    connectivity.push_back(cardinal);
                }
            }

            for(auto i = 0; i < data.totalNumCells(); i++)
            {
                offsets.push_back((i+1)*4);
                types.push_back(7);
                if(labels.size() == data.totalNumCells())
                    cellData.push_back(labels[i]);
                else
                    cellData.push_back(1);
            }
        }

        void write_file(std::string name) override {
            name += ".vtu";
            std::cout << "Writing data to " << name << "\n";
            
            std::vector<vtu11::DataSetInfo> dataSetInfo
            {
                {"Geocell", vtu11::DataSetType::CellData, 1},
            };
            vtu11::Vtu11UnstructuredMesh mesh {points, connectivity, offsets, types};
            vtu11::writeVtu(name, mesh, dataSetInfo, {cellData}, "Ascii");
        }

        void close_file() override {
            std::cout << "Closing the file\n";
            points.clear();
            types.clear();
            offsets.clear();
            connectivity.clear();
            cellData.clear();
        }

    private:
        std::string filename;

        VTK::MeshData meshData;
        
        std::vector<double> points;
        std::vector<vtu11::VtkCellType> types;
        std::vector<vtu11::VtkIndexType> offsets;
        std::vector<vtu11::VtkIndexType> connectivity;
        std::vector<double> cellData;
        //std::vector<vtu11::DataSetInfo> dataSetInfo;
};

} //end namespace Cajete

#endif 
