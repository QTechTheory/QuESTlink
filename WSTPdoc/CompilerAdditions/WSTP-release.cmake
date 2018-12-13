#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "wstp" for configuration "Release"
set_property(TARGET wstp APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(wstp PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/CompilerAdditions/wstp.framework/Versions/4.36/wstp"
  IMPORTED_SONAME_RELEASE "@executable_path/../Frameworks/wstp.framework/Versions/4.36/wstp"
  )

list(APPEND _IMPORT_CHECK_TARGETS wstp )
list(APPEND _IMPORT_CHECK_FILES_FOR_wstp "${_IMPORT_PREFIX}/CompilerAdditions/wstp.framework/Versions/4.36/wstp" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
