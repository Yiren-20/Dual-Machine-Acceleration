cmake_minimum_required(VERSION 3.10)
project(MyProject)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(include)

# 设置编译器选项
if (MSVC)

    set(CMAKE_CXX_FLAGS /source-charset:gb2312)
    find_package(OpenMP REQUIRED)
    if (OPENMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
    set(TESTWINDOWS src/test_windows.cpp)
    set(WINCLIENT src/windows_client.cpp)
    set(WINSERVER src/windows_server.cpp)
    add_executable(test_windows ${TESTWINDOWS})
    add_executable(windows_server ${WINSERVER})
    add_executable(windows_client ${WINCLIENT})
else ()

    set(CMAKE_CXX_FLAGS "-O2 -fopenmp -g -msse")
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
    set(SOURCE_FILES src/test_linux.cpp)
    set(SOURCE_FILES1 src/linux_server.cpp)
    set(SOURCE_FILES2 src/linux_client.cpp)
    add_executable(test_linux ${SOURCE_FILES})
    add_executable(linux_server ${SOURCE_FILES1})
    add_executable(linux_client ${SOURCE_FILES2})

endif()


