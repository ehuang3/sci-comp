project(example_test)

file(GLOB src "*.cpp")
file(GLOB hdr "*.h")

add_executable(example ${src} ${hdr})
target_link_libraries(example sci_comp)
set_target_properties(example PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
