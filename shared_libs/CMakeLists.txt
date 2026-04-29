cmake_minimum_required(VERSION 3.13)
project(my_aspect_plugins LANGUAGES CXX)

# 👇 直接写死 ASPECT 路径
set(ASPECT_DIR "/home/bbkhangq/softwares/aspect_test/build_self_gravitation")

find_package(Aspect REQUIRED HINTS ${ASPECT_DIR})

DEAL_II_INITIALIZE_CACHED_VARIABLES()

set(TARGET geoid_self_gravitation)
add_library(${TARGET} SHARED
  geoid_self_gravitation.cc
  geoid_self_gravitation.h
)

ASPECT_SETUP_PLUGIN(${TARGET})

set_target_properties(${TARGET} PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
)