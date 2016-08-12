
                cmake_minimum_required(VERSION 3.0)
                project(verify_mpi_compat VERSION 1.0 LANGUAGES CXX)
                file(WRITE compiler_info.cmake "
                    set(MPI_CXX_COMPILER_VERSION ${CMAKE_CXX_COMPILER_VERSION})
                    set(MPI_CXX_COMPILER_ID ${CMAKE_CXX_COMPILER_ID})
                    ")
            