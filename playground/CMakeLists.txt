add_executable(playground playground.cpp)
target_link_libraries(playground PRIVATE SeQuant::SeQuant)

add_executable(icMRCC icMRCC.cpp)
target_link_libraries(icMRCC PRIVATE SeQuant::SeQuant)
