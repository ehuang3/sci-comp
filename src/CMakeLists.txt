project(sci-comp)

# Enable multi-threaded compilation.
# We do this here and not in the root folder since the example apps
# do not have enough source files to benefit from this.
# if(MSVC)
# 	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
# endif()

# Setup directory structure
set(SUBDIRS utils)
foreach(SUBDIR ${SUBDIRS})
	file(GLOB_RECURSE dir_files ${SUBDIR}/*.h ${SUBDIR}/*.hpp ${SUBDIR}/*.cpp)
	source_group(${SUBDIR} FILES ${dir_files})
endforeach()

# Search all header and source files
file(GLOB_RECURSE sci_comp_srcs *.cpp)
file(GLOB_RECURSE sci_comp_hdrs *.h *.hpp)

# Library
add_library(sci_comp ${sci_comp_srcs} ${sci_comp_hdrs})
# target_link_libraries(sci_comp ${DART_LIBRARIES} ${DARTExt_LIBRARIES})
set_target_properties(sci_comp PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib" ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
