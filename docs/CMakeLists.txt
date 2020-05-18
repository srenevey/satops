# first we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)

# check if Doxygen is installed
if (BUILD_DOC)
    find_package(Doxygen)
    if (DOXYGEN_FOUND)

        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

        # note the option ALL which allows to build the docs together with the application
        add_custom_target( doc ALL
                COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_CONFIG}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMENT "Generating API documentation with Doxygen"
                VERBATIM )
    else (DOXYGEN_FOUND)
        message("Doxygen need to be installed to generate the doxygen documentation")
    endif (DOXYGEN_FOUND)
endif()