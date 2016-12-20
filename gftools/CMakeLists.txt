include_directories(${GFTools_SOURCE_DIR}/gftools)
include_directories(${GFTools_SOURCE_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${ALPSCore_INCLUDES})

set (base_src
    defaults.hpp
    tools.hpp
    num_io.hpp
    grid_base.hpp 
    matsubara_grid.hpp 
    real_grid.hpp 
    kmesh.hpp 
    nkmesh.hpp
    index_grid.hpp
    hdf5.hpp
    container.hpp container.hxx
    grid_object.hpp grid_object.hxx
    )

add_library(libgftools OBJECT ${base_src})

install(DIRECTORY . DESTINATION include/gftools
          FILES_MATCHING PATTERN "*.hpp" PATTERN "*.hxx"
        )
