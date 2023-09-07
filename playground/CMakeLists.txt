add_library(utility OBJECT
	Utils.cpp
)
target_link_libraries(utility PUBLIC SeQuant::SeQuant)
target_include_directories(utility PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(playground playground.cpp)
target_link_libraries(playground PRIVATE SeQuant::SeQuant utility)

add_executable(icMRCC icMRCC.cpp)
target_link_libraries(icMRCC PRIVATE SeQuant::SeQuant utility)

add_executable(repro repro.cpp)
target_link_libraries(repro PRIVATE SeQuant::SeQuant utility)

add_executable(playground_ccd ccd.cpp)
target_link_libraries(playground_ccd PRIVATE SeQuant::SeQuant utility)
set_target_properties(playground_ccd PROPERTIES OUTPUT_NAME "ccd")

add_executable(cc_sequant cc_sequant.cpp)
target_link_libraries(cc_sequant PRIVATE SeQuant::SeQuant)

add_executable(compare_spintrace compare.cpp)
target_link_libraries(compare_spintrace PRIVATE SeQuant::SeQuant)
