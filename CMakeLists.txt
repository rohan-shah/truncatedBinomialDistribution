cmake_minimum_required(VERSION 3.1)
project(truncatedBinomialDistribution)
set(Boost_USE_STATIC_LIBS       OFF)
set(Boost_USE_MULTITHREADED     ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package(Boost 1.53.0 REQUIRED)
if(WIN32)
	find_package(mpir)
	find_package(mpfr)
endif()
add_subdirectory(TruncatedBinomialDistribution)
add_subdirectory(TruncatedBinomialTest)
export(TARGETS truncatedBinomialDistribution FILE truncatedBinomialDistributionConfig.cmake)