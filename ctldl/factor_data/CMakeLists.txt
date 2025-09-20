if(BUILD_TESTING)
    add_executable(tests_unit_factor_data
        runner_tests_unit_factor_data.cpp
        factorization_repeating_block_tridiagonal_arrowhead_linked.test.cpp
    )
    target_link_libraries(tests_unit_factor_data PRIVATE
        ctldl
        Boost::unit_test_framework
    )

    add_test(NAME unit_factor_data COMMAND tests_unit_factor_data)
endif()
