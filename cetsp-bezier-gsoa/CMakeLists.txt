cmake_minimum_required(VERSION 2.8)
project(cetsp-bezier-gsoa)

set(SOURCE_FILES
   src/velocity_profile.cc
   src/bezier_maneuver.cc
   src/ring.cc
   src/gsoa_cetsp_bezier.cc
   )

add_executable(tcetsp-bezier-gsoa src/tcetsp-bezier-gsoa.cc ${SOURCE_FILES})
target_include_directories(tcetsp-bezier-gsoa PUBLIC ${CMAKE_BINARY_DIR}/include) 
target_link_libraries(tcetsp-bezier-gsoa crl-algorithm crl-gui)

install(TARGETS tcetsp-bezier-gsoa DESTINATION .)
install(FILES tcetsp-bezier-gsoa.cfg DESTINATION .)
install(DIRECTORY etc DESTINATION .)
