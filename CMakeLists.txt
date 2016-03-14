cmake_minimum_required(VERSION 3.3)
project(lab1)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")

set(SOURCE_FILES lab1.c)
add_executable(lab1 ${SOURCE_FILES})
