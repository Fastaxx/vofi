add_executable(gen_golden gen_golden.c)
target_include_directories(gen_golden PRIVATE "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(gen_golden vofi m)

add_executable(test_regression test_regression.c)
target_include_directories(test_regression PRIVATE "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(test_regression vofi m)
add_test(NAME regression_2d3d
         COMMAND test_regression "${CMAKE_CURRENT_SOURCE_DIR}/golden.dat")

add_executable(stress_4d stress_4d.c)
target_include_directories(stress_4d PRIVATE "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(stress_4d vofi m)
add_test(NAME stress_4d COMMAND stress_4d 250 1)

add_executable(test_4d test_4d.c)
target_include_directories(test_4d PRIVATE "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(test_4d vofi m)
add_test(NAME hyper_4d COMMAND test_4d)
