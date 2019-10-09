cmake_minimum_required(VERSION 3.7)
project(sci_comp_exercises)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/esl")

file(GLOB_RECURSE CPP_FILES "src/seq/[^-]*.cpp")
foreach(CPP_FILE ${CPP_FILES})
	get_filename_component(F_NAME ${CPP_FILE} NAME_WE)
	add_executable(${F_NAME} ${CPP_FILE})
	target_include_directories(${F_NAME} PUBLIC "include")
	target_compile_features(${F_NAME} PUBLIC cxx_std_17)
	target_compile_definitions(${F_NAME} PUBLIC _USE_MATH_DEFINES)
	target_compile_options(${F_NAME} PUBLIC -Wall -Wpedantic -Wextra -Werror=return-type -m64 -march=native 
		$<$<CONFIG:DEBUG>:-O0 -g> $<$<CONFIG:RELEASE>:-Wno-unused-parameter -Wno-deprecated-declarations -O3 -DNDEBUG>)
	target_link_libraries(${F_NAME} esl)
endforeach()

find_package(MPI)

if(MPI_CXX_FOUND)
	file(GLOB CPP_FILES "src/mpi/[^-]*.cpp")
		foreach(CPP_FILE ${CPP_FILES})
			get_filename_component(F_NAME ${CPP_FILE} NAME_WE)
			add_executable(${F_NAME} ${CPP_FILE})
			target_include_directories(${F_NAME} PUBLIC "include" ${MPI_INCLUDE_PATH})
			target_link_libraries(${F_NAME} ${MPI_LIBRARIES})
			if(MPI_COMPILE_FLAGS)
				set_target_properties(${F_NAME} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
			endif()
			if(MPI_LINK_FLAGS)
				set_target_properties(${F_NAME} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
			endif()
			target_compile_features(${F_NAME} PUBLIC cxx_std_17)
			target_compile_definitions(${F_NAME} PUBLIC _USE_MATH_DEFINES)
			target_compile_options(${F_NAME} PUBLIC -Wall -Wpedantic -Wextra -Werror=return-type -m64 -march=native 
				$<$<CONFIG:DEBUG>:-O0 -g> $<$<CONFIG:RELEASE>:-Wno-unused-parameter -Wno-deprecated-declarations -O3 -DNDEBUG>)
			target_link_libraries(${F_NAME} esl)
	endforeach()
endif()

