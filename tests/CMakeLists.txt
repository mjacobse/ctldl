add_subdirectory(utility)

add_executable(tests_integration
    multiply_factorize_solve_correct/repeating/larger_examples_good_permutation.cpp
    multiply_factorize_solve_correct/repeating/random_examples.cpp
    multiply_factorize_solve_correct/repeating/small_examples_all_permutations.cpp
    multiply_factorize_solve_correct/repeating_arrowhead_linked/random_examples.cpp
    multiply_factorize_solve_correct/single/larger_examples_good_permutation.cpp
    multiply_factorize_solve_correct/single/small_examples_all_permutations.cpp
    runner_tests_integration.cpp
    utility/solution_generator.cpp
)
target_link_libraries(tests_integration PRIVATE
    ctldl
    ctldl_warnings
    Boost::unit_test_framework
)

add_test(NAME integration COMMAND tests_integration)
