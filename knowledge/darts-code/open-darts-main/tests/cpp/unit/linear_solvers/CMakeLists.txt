
# Adds all tests. To do that, it is assumed that each .cpp file is a test and the
# test's name is the filename without the extension.

# Collect all the names of the cpp files
file(GLOB CPPFILES  "*.cpp")

# For every cpp found
foreach(CPPFILE ${CPPFILES})
	# Extract the actual file name
	get_filename_component(FILENAME ${CPPFILE} NAME)

	# Extract the filename minus the cpp. This will be the name of target
	get_filename_component(TARGETNAME_BASE ${CPPFILE} NAME_WE)
	
	# To distinguish tests from different components we add a prefix 
	set(TARGETNAME "solvers_${TARGETNAME_BASE}")

	# Make the target
	add_executable(${TARGETNAME} ${FILENAME})

	target_link_libraries(${TARGETNAME}
		opendarts_compiler_flags
		linear_solvers_test_common
    linear_solvers
  )

  # Place the target test in the tests folder and rename the source folder to src
  set_property(TARGET ${TARGETNAME} PROPERTY FOLDER "tests//unit//solvers")
  source_group(src FILES ${FILENAME})

	# Add the test
	add_test(NAME ${TARGETNAME} COMMAND ${TARGETNAME})

	# Install the test executable
	# install(TARGETS ${TARGETNAME} DESTINATION bin/openDARTS/tests/unit/linear_solvers)
endforeach()
