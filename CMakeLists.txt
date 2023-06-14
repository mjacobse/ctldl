cmake_minimum_required(VERSION 3.1.0)
project(ctldl)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(ctldl_warnings INTERFACE)
if(MSVC)
	target_compile_options(ctldl_warnings INTERFACE /W4)
else()
	target_compile_options(ctldl_warnings INTERFACE
		-Wall -Wextra -Wpedantic -pedantic-errors -Wconversion -Wsign-conversion
		-fdiagnostics-color
	)
endif()

add_library(ctldl_includes INTERFACE)
target_include_directories(ctldl_includes INTERFACE .)

add_subdirectory(ctldl)
add_subdirectory(examples)
add_subdirectory(external)
add_subdirectory(tests)
