function(add_example_target TARGET_NAME TARGET_SOURCE)
    add_executable(${TARGET_NAME} ${TARGET_SOURCE})
    target_link_libraries(${TARGET_NAME} PRIVATE
        ctldl::ctldl
    )
endfunction()

function(add_example EXAMPLE_NAME)
    add_example_target(example_${EXAMPLE_NAME} example_${EXAMPLE_NAME}.cpp)
endfunction()

function(add_example_target_mtx PREFIX EXAMPLE_NAME)
    set(TARGET_NAME ${PREFIX}_mtx_${EXAMPLE_NAME})
    add_example_target(${TARGET_NAME} ${PREFIX}_mtx.cpp)
    string(REPLACE _ / MTX_INCLUDE_DIRECTORY ${EXAMPLE_NAME})
    target_include_directories(${TARGET_NAME} PRIVATE
        mtx_includes/${MTX_INCLUDE_DIRECTORY}
    )
    target_link_libraries(${TARGET_NAME} PRIVATE
        ctldl::ctldl_fileio_mtx
    )
endfunction()

function(add_example_mtx EXAMPLE_NAME)
    add_example_target_mtx(example ${EXAMPLE_NAME})
    if(CTLDL_BUILD_BENCHMARKS)
        add_example_target_mtx(benchmark ${EXAMPLE_NAME})
        target_link_libraries(benchmark_mtx_${EXAMPLE_NAME} PRIVATE
            benchmark::benchmark
        )
    endif()
endfunction()
