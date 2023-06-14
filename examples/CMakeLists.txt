function(add_example_target TARGET_NAME TARGET_SOURCE)
    add_executable(${TARGET_NAME} ${TARGET_SOURCE})
    target_link_libraries(${TARGET_NAME} PRIVATE
        ctldl
        ctldl_warnings
    )
endfunction()

function(add_example EXAMPLE_NAME)
    add_example_target(example_${EXAMPLE_NAME} example_${EXAMPLE_NAME}.cpp)
endfunction()

function(add_example_mtx EXAMPLE_NAME)
    add_example_target(example_${EXAMPLE_NAME} example_mtx.cpp)
    target_include_directories(example_${EXAMPLE_NAME} PRIVATE
        mtx_includes/${EXAMPLE_NAME}
    )
    target_link_libraries(example_${EXAMPLE_NAME} PRIVATE
        ctldl_fileio_mtx
    )
endfunction()

add_example(nos2)
add_example(nos4)
add_example(tridiagonal)

add_example_mtx(mhdb416)

add_example_mtx(laufkatze101)
add_example_mtx(opal11)
add_example_mtx(spline51)
