#!/bin/sh

cd build_directory
cmake ..
cmake --build ${PWD} --parallel $(nproc)
cp compile_commands.json ..

time ./simulator



