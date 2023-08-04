find_package(Boost COMPONENTS unit_test_framework REQUIRED)

add_executable(tests_integration
    runner_tests_integration.cpp
    test_multiply_factorize_solve_correct_repeating.cpp
    test_multiply_factorize_solve_correct_single.cpp
    utility/solution_generator.cpp
)
target_link_libraries(tests_integration PRIVATE
    ctldl
    ctldl_warnings
    Boost::unit_test_framework
)

add_test(NAME integration COMMAND tests_integration)
