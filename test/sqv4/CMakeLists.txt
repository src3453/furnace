cmake_minimum_required(VERSION 3.16)
project(sqv4_check LANGUAGES C CXX)
set(ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../..")
add_executable(sqv4_check codec.cpp
  ${ROOT}/src/engine/sample.cpp ${ROOT}/src/engine/filter.cpp
  ${ROOT}/src/engine/safeReader.cpp ${ROOT}/src/engine/safeWriter.cpp
  ${ROOT}/src/engine/brrUtils.c
  ${ROOT}/src/engine/platform/sound/3ws8pn/sound.cpp
  ${ROOT}/extern/adpcm/bs_codec.c ${ROOT}/extern/adpcm/oki_codec.c
  ${ROOT}/extern/adpcm/yma_codec.c ${ROOT}/extern/adpcm/ymb_codec.c
  ${ROOT}/extern/adpcm/ymz_codec.c ${ROOT}/extern/adpcm-xq-s/adpcm-lib.c)
target_compile_features(sqv4_check PRIVATE cxx_std_17)
target_compile_definitions(sqv4_check PRIVATE FMT_HEADER_ONLY)
target_include_directories(sqv4_check PRIVATE ${ROOT}/extern/fmt/include)
if(MSVC)
  target_compile_options(sqv4_check PRIVATE /UNDEBUG /utf-8)
else()
  target_compile_options(sqv4_check PRIVATE -UNDEBUG)
endif()
enable_testing()
add_test(NAME sqv4_codec_sample_core COMMAND sqv4_check)

add_executable(s3w2_revision_check revision.cpp
  ${ROOT}/src/engine/platform/sound/3ws8pn/sound.cpp)
target_compile_features(s3w2_revision_check PRIVATE cxx_std_17)
if(MSVC)
  target_compile_options(s3w2_revision_check PRIVATE /UNDEBUG /utf-8)
else()
  target_compile_options(s3w2_revision_check PRIVATE -UNDEBUG)
endif()
add_test(NAME s3w2_revisions COMMAND s3w2_revision_check)
