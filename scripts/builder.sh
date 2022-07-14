#!/bin/bash

arch=$(uname -m)
compiler=g++
compiler_version=$(g++ --version)
env_prefix=/home/erock/Documents/Code/CajeteRevamped
project_name=CajeteProto
build_name=build_basic
source_prefix=${env_prefix}/${project_name}
build_base_dir=${env_prefix}/${build_name}
build_prefix=${build_base_dir}/build/${arch}-${compiler}
install_prefix=${build_base_dir}/install/${arch}-${compiler}
build_dir=${build_prefix}/${project_name}
install_dir=${install_prefix}/${project_name}

sundials_dir=/home/erock/Extras/sundials-build

source_dir=${source_prefix}

#compiler_path=

if [[ "$1" == "-c" ]]; then
    echo Removing ${build_base_dir} ...
    rm -rf ${build_base_dir}
fi

echo 'Build in ' ${build_dir}
echo 'Install in ' ${install_dir}

mkdir -p ${build_dir} && cd ${build_dir}

echo `pwd`
cmake   -D CMAKE_CXX_EXTENSIONS=Off \
        -D CMAKE_INSTALL_PREFIX=${install_path} \
        -D SUNDIALS_DIR=${sundials_dir} \
        ${source_dir}
