file(GLOB sourcefiles
    "../gas-optics/*.F90"
    "*.F90")

add_library(rrtmgpf STATIC ${sourcefiles})
