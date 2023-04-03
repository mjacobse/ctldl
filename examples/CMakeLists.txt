add_executable(example_tridiagonal
    example_tridiagonal.cpp
)
add_executable(example_nos4
    example_nos4.cpp
)
target_link_libraries(example_tridiagonal PRIVATE
    ctldl
    ctldl_warnings
)
target_link_libraries(example_nos4 PRIVATE
    ctldl
    ctldl_warnings
)
