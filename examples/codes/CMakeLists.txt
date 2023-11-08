set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

add_library(utility OBJECT
	Utils.cpp
)
target_link_libraries(utility PUBLIC SeQuant::SeQuant)
target_include_directories(utility PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(cc_sequant cc_sequant.cpp)
target_link_libraries(cc_sequant PRIVATE SeQuant::SeQuant)

