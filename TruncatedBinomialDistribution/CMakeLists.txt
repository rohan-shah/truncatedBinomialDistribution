cmake_minimum_required(VERSION 3.1)
set(CppFiles TruncatedBinomialDistribution.cpp)
set(HeaderFiles includeMPFRTruncatedBinomial.h serializeGMPTruncatedBiniomal.hpp TruncatedBinomialDistribution.h)
add_library(truncatedBinomialDistribution ${CppFiles} ${HeaderFiles})
target_include_directories(truncatedBinomialDistribution PRIVATE ${Boost_INCLUDE_DIRS})
target_include_directories(truncatedBinomialDistribution PUBLIC ${PROJECT_SOURCE_DIR}/TruncatedBinomialDistribution)
target_compile_definitions(truncatedBinomialDistribution PUBLIC BOOST_ALL_NO_LIB)
if(WIN32)
	target_link_libraries(truncatedBinomialDistribution PUBLIC mpir mpfr)
else()
	target_compile_features(truncatedBinomialDistribution PRIVATE cxx_rvalue_references)
	target_compile_options(truncatedBinomialDistribution PUBLIC "-fPIC")
endif()
