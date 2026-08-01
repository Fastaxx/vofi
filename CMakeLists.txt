cmake_minimum_required(VERSION 3.10)

project(vofi_project
	LANGUAGES C
)

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

add_subdirectory(src)

option(BUILD_TESTING "Build the vofi test suite" ON)
if(BUILD_TESTING)
	enable_testing()
	add_subdirectory(tests)
endif()

