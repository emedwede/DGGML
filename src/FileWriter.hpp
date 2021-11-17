#ifndef FILE_WRITER_HPP
#define FILE_WRITER_HPP

namespace Cajete 
{
    //Template Desing Pattern for the file writer workflow

class FileWriter {
    public:

        void save() {
            create_file();
            
            open_file();

            format_data();

            write_file();

            close_file();
        }

    protected:

        virtual void create_file() const = 0;

        virtual void open_file() const = 0;

        virtual void format_data() const = 0;

        virtual void write_file() const = 0;

        virtual void close_file() const = 0;
};

} //end namespace cajete

#endif
