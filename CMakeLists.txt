cmake_minimum_required(VERSION 2.6)
project(probe C)

if(UNIX)
  find_library(MATH_LIBRARY m)
endif()

add_executable(probe
  probe.c dots.c abin.c readPDBrecs.c geom3d.c utility.c select.c
  parse.c atomprops.c stdconntable.c autobondrot.c hybrid_36_c.c)
target_link_libraries(probe ${MATH_LIBRARY})

install(TARGETS probe DESTINATION bin)
