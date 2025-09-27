add_library(ctldl)
add_library(ctldl::ctldl ALIAS ctldl)

target_sources(ctldl PRIVATE
    utility/contracts.cpp
)

target_link_libraries(ctldl PUBLIC
    ctldl_includes
)

add_subdirectory(factor_data)
add_subdirectory(fileio)
add_subdirectory(symbolic)
