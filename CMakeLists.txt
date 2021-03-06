cmake_minimum_required(VERSION 3.10.0)
project(probe)

if (WIN32)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif ()

set(probe_SOURCES
	abin.c
	atomprops.c
	autobondrot.c
	dots.c
	geom3d.c
	hybrid_36_c.c
	readPDBrecs.c
	parse.c
	select.c
	stdconntable.c
	utility.c
	probe.c
)

set(probe_HEADERS
	abin.h
	atomprops.h
	autobondrot.h
	dots.h
	geom3d.h
	hybrid_36_c.h
	parse.h
	readPDBrecs.h
	select.h
	stdconntable.h
	utility.h
)

add_executable(probe ${probe_SOURCES} ${probe_HEADERS})
if (NOT WIN32)
  target_link_libraries(probe m)
endif()
