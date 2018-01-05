PROJECT(FMIX)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)

set(SRC ${CMAKE_SOURCE_DIR}/src)
set(SOURCES
  ${SRC}/forces_ext.cpp
  ${SRC}/geometry.cpp
  ${SRC}/solide.cpp
  ${SRC}/main.cpp
  )

add_executable(Mka3d ${SOURCES})
