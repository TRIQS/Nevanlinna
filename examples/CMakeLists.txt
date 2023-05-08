add_executable(ndarray_example ndarray_example.cpp)
target_compile_definitions(ndarray_example PUBLIC -DDATA_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data")
target_link_libraries(ndarray_example ${PROJECT_NAME}::${PROJECT_NAME}_c)