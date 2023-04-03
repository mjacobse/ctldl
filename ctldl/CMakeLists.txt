add_library(ctldl INTERFACE)
target_sources(ctldl
    INTERFACE
        ${CMAKE_CURRENT_LIST_DIR}/factorization_repeating_block_tridiagonal.hpp
)

target_link_libraries(ctldl INTERFACE
    ctldl_includes
)
