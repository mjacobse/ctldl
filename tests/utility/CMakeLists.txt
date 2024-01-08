add_executable(tests_unit_testutility
    random/procedural_test_permutation.test.cpp
    runner_tests_unit_testutility.cpp
)
target_link_libraries(tests_unit_testutility PRIVATE
    ctldl
    ctldl_warnings
    Boost::unit_test_framework
)

add_test(NAME unit_testutility COMMAND tests_unit_testutility)
