add_executable(ndarray_example ndarray_example.cpp)
target_compile_definitions(ndarray_example PUBLIC -DDATA_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data")
target_link_libraries(ndarray_example ${PROJECT_NAME}::${PROJECT_NAME}_c)

add_executable(ndarray_caratheodory_example ndarray_caratheodory_example.cpp)
target_compile_definitions(ndarray_caratheodory_example PUBLIC -DDATA_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data")
target_link_libraries(ndarray_caratheodory_example ${PROJECT_NAME}::${PROJECT_NAME}_c)