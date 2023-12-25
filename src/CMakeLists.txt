set(TARGET_NAME "seidel_main")

file(GLOB_RECURSE SRCS "*.[ch]pp")

add_subdirectory(linalg)

add_executable(${TARGET_NAME} ${SRCS})
target_link_libraries(${TARGET_NAME} PRIVATE linalg)
