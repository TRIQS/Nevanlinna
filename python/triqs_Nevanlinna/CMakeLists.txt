# Configure the version
configure_file(version.py.in version.py)

# All Python files. Copy them in the build dir to have a complete package for the tests.
file(GLOB_RECURSE python_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
file(GLOB_RECURSE wrap_generators RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *_desc.py)
list(REMOVE_ITEM python_sources "${wrap_generators}")

if(NOT TRIQS_Support)
  list(REMOVE_ITEM python_sources "solver.py")
  list(REMOVE_ITEM wrap_generators "solver_core_desc.py")
endif()

foreach(file ${python_sources})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# Install python files to proper location
set(PYTHON_LIB_DEST ${CPP2PY_PYTHON_LIB_DEST_ROOT}/${PROJECT_NAME})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${PYTHON_LIB_DEST})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION ${CPP2PY_PYTHON_LIB_DEST_ROOT} FILES_MATCHING PATTERN "*.py" PATTERN "*_desc.py" EXCLUDE)

# Build and install any python modules
foreach(gen ${wrap_generators})
  string(REPLACE "_desc.py" "" gen ${gen})
  get_filename_component(module_name ${gen} NAME_WE)
  get_filename_component(module_dir ${gen} DIRECTORY)
  add_cpp2py_module(NAME ${module_name} DIRECTORY ${module_dir})
  add_library(${PROJECT_NAME}::${module_name} ALIAS ${module_name})
  target_link_libraries(${module_name} ${PROJECT_NAME}_c nda::nda_py $<$<BOOL:${TRIQS_Support}>:triqs_py>)
  install(TARGETS ${module_name} DESTINATION ${PYTHON_LIB_DEST}/${module_dir})
endforeach()
