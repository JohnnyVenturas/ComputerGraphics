#!/bin/sh

cd build_directory

cmake --build ${PWD} --parallel $(nproc)
cp compile_commands.json ..

time ./raytracer

firefox image.png

