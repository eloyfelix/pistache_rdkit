cmake_minimum_required(VERSION 3.6)
project(PistacheRDKit)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(PkgConfig)
pkg_check_modules(Pistache REQUIRED IMPORTED_TARGET libpistache)

include_directories("/usr/include/rdkit")
set(RDKit_LIBS RDKitSmilesParse RDKitFileParsers RDKitGraphMol RDKitFingerprints RDKitInchi RDKitRDInchiLib RDKitFMCS RDKitFilterCatalog RDKitRDGeneral RDKitCatalogs RDKitDescriptors)

add_executable(${PROJECT_NAME} src/server.cc)
target_link_libraries(${PROJECT_NAME} PkgConfig::Pistache pthread ${RDKit_LIBS})
