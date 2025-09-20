add_library(ctldl_fileio_mtx)
target_sources(ctldl_fileio_mtx PRIVATE
    mtx_check.cpp
    mtx_foreach_entry.cpp
    mtx_read_header.cpp
)

target_link_libraries(ctldl_fileio_mtx PRIVATE
    ctldl_includes
)
