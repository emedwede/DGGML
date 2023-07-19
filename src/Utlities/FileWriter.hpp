#ifndef FILE_WRITER_HPP
#define FILE_WRITER_HPP

#include <string>

namespace DGGML
{
    //Template Desing Pattern for the file writer workflow

template <typename DataType>
class FileWriter {
    public:
        
        using data_type = DataType;
        //save graph data as filename 
        void save(DataType data, std::string name) {
            create_file();
            
            open_file();

            format_data(data);

            write_file(name);

            close_file();
        }

    protected:

        virtual void create_file() const = 0;

        virtual void open_file() const = 0;

        virtual void format_data(DataType data) = 0;

        virtual void write_file(std::string name) = 0;

        virtual void close_file() = 0;
};

} //end namespace DGGML

#endif
