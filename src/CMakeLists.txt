set(TARGET_NAME "seidel_main")

file(GLOB_RECURSE SRCS "*.[ch]pp")

add_executable(${TARGET_NAME} ${SRCS})