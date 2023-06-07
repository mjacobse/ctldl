function(add_example EXAMPLE_NAME)
    set(TARGET_NAME example_${EXAMPLE_NAME})
    add_executable(${TARGET_NAME}
        ${TARGET_NAME}.cpp
    )
    target_link_libraries(${TARGET_NAME} PRIVATE
        ctldl
        ctldl_warnings
    )
endfunction()

add_example(nos2)
add_example(nos4)
add_example(tridiagonal)
