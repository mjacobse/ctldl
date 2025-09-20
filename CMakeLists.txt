cmake_minimum_required(VERSION 3.1.0)
project(ctldl)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

option(
	CTLDL_BUILD_BENCHMARKS
	"Enable building benchmarks. Default: ${PROJECT_IS_TOP_LEVEL}. Values: { ON, OFF }."
	${PROJECT_IS_TOP_LEVEL}
)
if(CTLDL_BUILD_BENCHMARKS)
	find_package(benchmark REQUIRED)
endif()

add_library(ctldl_includes INTERFACE)
target_include_directories(ctldl_includes INTERFACE .)

include(CTest)
if(BUILD_TESTING)
	find_package(Boost COMPONENTS unit_test_framework REQUIRED)
	add_subdirectory(tests)
endif()

add_subdirectory(ctldl)
add_subdirectory(examples)
