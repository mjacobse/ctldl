if(BUILD_TESTING)
    add_executable(tests_unit_symbolic
        runner_tests_unit_symbolic.cpp
        compute_elimination_tree_blocked.test.cpp
        filled_in_sparsity_blocked.test.cpp
        filled_in_sparsity_repeating_arrowhead.test.cpp
    )
    target_link_libraries(tests_unit_symbolic PRIVATE
        ctldl
        ctldl_warnings
        Boost::unit_test_framework
    )

    add_test(NAME unit_symbolic COMMAND tests_unit_symbolic)
endif()
