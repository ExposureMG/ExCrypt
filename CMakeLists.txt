cmake_minimum_required(VERSION 3.20)
project(ExCrypt CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(EXCRYPT_SOURCES
    "src/excrypt_aes.c"
    "src/excrypt_bn_key.cpp"
    "src/excrypt_bn_mod.cpp"
    "src/excrypt_bn_pkcs1.cpp"
    "src/excrypt_bn_rsa.cpp"
    "src/excrypt_bn_sig.c"
    "src/excrypt_bn.c"
    "src/excrypt_des.c"
    "src/excrypt_ecc.c"
    "src/excrypt_md5.c"
    "src/excrypt_mem.cpp"
    "src/excrypt_parve.c"
    "src/excrypt_rc4.c"
    "src/excrypt_rotsum.c"
    "src/excrypt_sha.c"
    "src/excrypt_sha2.c"
    "src/exkeys.cpp"
    "src/rijndael.c"
)

set(EXCRYPT_HEADERS
    "src/excrypt_aes.h"
    "src/excrypt_bn.h"
    "src/excrypt_des_data.h"
    "src/excrypt_des.h"
    "src/excrypt_ecc.h"
    "src/excrypt_md5.h"
    "src/excrypt_mem.h"
    "src/excrypt_parve.h"
    "src/excrypt_rc4.h"
    "src/excrypt_rotsum.h"
    "src/excrypt_sha.h"
    "src/excrypt_sha2.h"
    "src/excrypt.h"
    "src/exkeys.h"
    "src/fp61.hpp"
    "src/rijndael.h"
)

add_library(ExCrypt STATIC 
    ${EXCRYPT_SOURCES}
    ${EXCRYPT_HEADERS}
)

target_include_directories(ExCrypt PUBLIC "src")

if(WIN32)
    target_link_libraries(ExCrypt PRIVATE bcrypt)
endif()