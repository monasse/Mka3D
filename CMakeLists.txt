PROJECT(FMIX)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)

set(SRC ${CMAKE_SOURCE_DIR}/src)
set(SOURCES
  ${SRC}/Particule.cpp
  ${SRC}/lien.cpp
  ${SRC}/system.cpp
  ${SRC}/Arbre.cpp
  ${SRC}/main.cpp
  )

add_executable(asynchrone ${SOURCES})
