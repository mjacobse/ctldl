add_library(ctldl INTERFACE)

target_link_libraries(ctldl INTERFACE
    ctldl_includes
)

add_subdirectory(factor_data)
add_subdirectory(fileio)
add_subdirectory(symbolic)
