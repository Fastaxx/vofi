cmake_minimum_required(VERSION 2.8)

project(vofi_project
	LANGUAGES C
)

add_library(vofi SHARED
	src/checkboundary.c
	src/checkconsistency.c
	src/getarclength.c
	src/getcc.c
	src/getintersections.c
	src/getlimits.c
	src/getmin.c
	src/gettype.c
	src/getzero.c
	src/integrate.c
	src/orderdirs.c
	src/tecplot.c
	src/triangulate.c
)
target_include_directories(vofi PUBLIC include)

install(TARGETS vofi DESTINATION lib)
install(FILES include/vofi.h DESTINATION include)

