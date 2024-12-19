#!/bin/bash

# 定义构建目录
BUILD_DIR=build

# 检查构建目录是否存在，如果不存在则创建
if [ ! -d "$BUILD_DIR" ]; then
  mkdir -p "$BUILD_DIR"
fi

# 进入构建目录
cd "$BUILD_DIR"

# 检查操作系统类型并执行对应的构建命令
if [ "$(uname)" == "Darwin" ]; then
  # macOS
  cmake ..
  make
  ./MyExecutable
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
  # Linux
  cmake ..
  make
  # ./MyExecutable
elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ] || [ "$(expr substr $(uname -s) 1 10)" == "MINGW64_NT" ]; then
  # Windows (MinGW)
  cmake -G "MinGW Makefiles" ..
  mingw32-make
  # ./MyExecutable.exe
  # ./MyExecutable1.exe
else
  echo "Unsupported operating system"
fi
