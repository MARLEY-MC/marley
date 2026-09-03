# CMakeLists.txt for the marg4 (MARLEY + Geant4) example
# Geant4 was found by the parent CMakeLists.txt; its variables
# (Geant4_LIBRARIES, Geant4_INCLUDE_DIRS) are available here.

add_executable(marg4
  src/marg4.cc
  src/DetectorConstruction.cc
  src/EventAction.cc
  src/PhysicsList.cc
  src/NeutronSource.cc
)

target_link_libraries(marg4 PRIVATE marley_lib ${Geant4_LIBRARIES})

target_include_directories(marg4 PRIVATE
  ${Geant4_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Suppress -Wdeprecated-copy from Geant4's bundled CLHEP headers
# (GCC 15+ makes this warning default-on)
target_compile_options(marg4 PRIVATE -Wno-deprecated-copy -Wno-unused-private-field)

set_target_properties(marg4 PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
