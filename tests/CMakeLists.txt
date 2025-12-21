add_executable(test_basis test_basis.cpp)
target_link_libraries(test_basis PRIVATE nurbs)
add_test(NAME test_basis COMMAND test_basis)

add_executable(test_curve test_curve.cpp)
target_link_libraries(test_curve PRIVATE nurbs)
add_test(NAME test_curve COMMAND test_curve)

add_executable(test_surface test_surface.cpp)
target_link_libraries(test_surface PRIVATE nurbs)
add_test(NAME test_surface COMMAND test_surface)
