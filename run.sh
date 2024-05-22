#!/bin/sh

cmake --build ${PWD} --parallel $(nproc)

time ./raytracer

xdg-open image.png

