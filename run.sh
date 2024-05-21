#!/bin/sh

cmake --build ${PWD} --parallel $(nproc)

./raytracer

firefox image.png

