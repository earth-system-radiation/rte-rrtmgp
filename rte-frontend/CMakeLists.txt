file(GLOB sourcefiles
    "*.F90")

add_library(rtef STATIC ${sourcefiles})
