cmake_minimum_required(VERSION 2.8)

project(vofi_project
	LANGUAGES C
)

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

add_subdirectory(src)

