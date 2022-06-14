cmake_minimum_required(VERSION 3.10.2)

project(PACE22_FVSP_HUST_SCP)

#set(CMAKE_EXE_LINKER_FLAGS  -static)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O1")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} - O3")

set(CMAKE_CXX_STANDARD 17)

link_directories(${PROJECT_SOURCE_DIR}/lib)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(PACE22_FVSP_HUST_SCP src/myList.cpp src/reduction.cpp src/Topo.cpp)

target_link_libraries(${PROJECT_NAME} libBC.a libNuMVC.a)


















