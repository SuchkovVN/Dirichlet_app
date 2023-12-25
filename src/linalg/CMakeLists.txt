set(TARGET_NAME "linalg")

file(GLOB_RECURSE SRCS "*.[ch]pp")
add_library(${TARGET_NAME} STATIC ${SRCS})