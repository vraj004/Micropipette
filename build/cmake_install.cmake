# Install script for directory: /Users/vrajagopal/Documents/cell_migration/sims/opencmiss/Micropipette

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "DEBUG")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/." TYPE EXECUTABLE FILES "/Users/vrajagopal/Documents/cell_migration/sims/opencmiss/Micropipette/build/Micropipette")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/./Micropipette" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/./Micropipette")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/vrajagopal/opencmiss/install/x86_64_darwin/clang-7.3-F6.1/openmpi_debug/debug/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/./Micropipette")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/vrajagopal/opencmiss/install/x86_64_darwin/clang-7.3-F6.1/no_mpi/debug/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/./Micropipette")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/vrajagopal/opencmiss/install/x86_64_darwin/clang-7.3-F6.1/openmpi_debug/debug/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/./Micropipette")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/vrajagopal/opencmiss/install/x86_64_darwin/clang-7.3-F6.1/mpi/openmpi/debug/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/./Micropipette")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/./Micropipette")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/Users/vrajagopal/Documents/cell_migration/sims/opencmiss/Micropipette/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
