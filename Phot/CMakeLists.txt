CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

#===============================================================================
# Load elements_subdir macro here 
# Examples:
#   For declaring a project module:
#         elements_subdir(ElementsExamples)
#===============================================================================
elements_subdir(Phot)

#===============================================================================
# Load elements_depends_on_subdirs macro here 
#   For creating a dependency onto an other accessible module
#         elements_depends_on_subdirs(ElementsKernel)
#===============================================================================
elements_depends_on_subdirs(ElementsKernel)

#===============================================================================
# Add the find_package macro (a pure CMake command) here to locate the 
# libraries.
# Examples:
#          find_package(CppUnit)
#===============================================================================

#===============================================================================
# Declare the library dependencies here
# Example:
#         elements_add_library(ElementsExamples src/Lib/*.cpp
#                     LINK_LIBRARIES Boost ElementsKernel
#                     INCLUDE_DIRS Boost ElementsKernel
#                     PUBLIC_HEADERS ElementsExamples)
#===============================================================================

#===============================================================================
# Declare the executables here
# Example:
# elements_add_executable(ElementsProgramExample src/Program/ProgramExample.cpp
#                        LINK_LIBRARIES Boost ElementsExamples
#                        INCLUDE_DIRS Boost ElementsExamples)
#===============================================================================

#===============================================================================
# Declare the Boost tests here
# Example:
# elements_add_unit_test(BoostClassExample tests/src/Boost/ClassExample_test.cpp
#                       EXECUTABLE BoostClassExample_test
#                       INCLUDE_DIRS ElementsExamples
#                       LINK_LIBRARIES ElementsExamples TYPE Boost)
#===============================================================================

#===============================================================================
# Declare the Python programs here
# Examples :
# elements_add_python_program(PythonProgramExample 
#                             ElementsExamples.PythonProgramExample)
#===============================================================================
elements_add_python_program(sextractor Phot.sextractor)
elements_add_python_program(mergecatalogs Phot.mergecatalogs)
elements_add_python_program(makeregionfile Phot.makeregionfile)
elements_add_python_program(truewcs Phot.truewcs)
elements_add_python_program(generate_weight_map Phot.generate_weight_map)

#===============================================================================
# Use the following macro for python modules, scripts and aux files:
#  elements_install_python_modules()
#  elements_install_scripts()
#===============================================================================
#elements_install_scripts()
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR})
file(COPY tests/config DESTINATION ${CMAKE_BINARY_DIR})

elements_install_scripts()
elements_install_python_modules()

#===============================================================================
# Add the elements_install_conf_files macro 
# Examples:
#          elements_install_conf_files()
#===============================================================================
elements_install_conf_files()
elements_install_aux_files()
