@echo off
chcp 65001
rm -rf build
mkdir build
cd build
cmake ..
cmake --build . --config Release
cd Release

cd ..
cd ..
