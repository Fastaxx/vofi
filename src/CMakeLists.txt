add_library(vofi
	checkboundary.c
	checkconsistency.c
	getarclength.c
	getcc.c
	getintersections.c
	getlimits.c
	getmin.c
	gettype.c
	checkboundary4d.c
	getintersections4d.c
	getlimits4d.c
	getmin4d.c
	getzero.c
	integrate.c
	integrate4d.c
	orderdirs.c
	orderdirs4d.c
	tecplot.c
	triangulate.c
)

target_include_directories(vofi PRIVATE "${PROJECT_SOURCE_DIR}/include")

set_target_properties(vofi PROPERTIES
	POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS})

install(TARGETS vofi DESTINATION lib)
install(FILES "${PROJECT_SOURCE_DIR}/include/vofi.h" DESTINATION include)

