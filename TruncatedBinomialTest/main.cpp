#include "TruncatedBinomialDistribution.h"
#include <cstddef>
#include <boost/random/mersenne_twister.hpp>
#include <iostream>
#include <boost/math/distributions/binomial.hpp>
namespace networkReliability
{
	void runTest(std::size_t binomialSize, float probability, std::size_t testSize, std::size_t truncationLevel, std::size_t lastAllowedValue, boost::mt19937& randomSource)
	{
		::TruncatedBinomialDistribution::TruncatedBinomialDistribution dist(binomialSize, truncationLevel, lastAllowedValue, probability);
		boost::scoped_array<std::size_t> counts(new std::size_t[binomialSize+1]);
		memset(counts.get(), 0, sizeof(std::size_t)*(binomialSize+1));
		for(std::size_t i = 0; i < testSize; i++)
		{
			counts[dist(randomSource)]++;
		}

		boost::math::binomial targetDistribution((double)binomialSize, probability);
		float absDifference = 0;
		float truncationFactor = (float)(boost::math::cdf(targetDistribution, lastAllowedValue) - boost::math::cdf(targetDistribution, truncationLevel-1));
		for(std::size_t i = truncationLevel; i <= lastAllowedValue; i++)
		{
			absDifference += (float)fabs(boost::math::pdf(targetDistribution, i)/truncationFactor -  (float)(counts[i]) / (float)(testSize));
		}
		std::cout << "Absolute difference from runTest(" << binomialSize << ", " << probability << ", " << testSize << ", " << truncationLevel << ") was " << absDifference << std::endl;
		if(absDifference > 0.0025f)
		{
			throw std::runtime_error("Estimated absolute difference of target distribution to empirical distribution was > 0.0025f");
		}
	}
	void runTest(std::size_t binomialSize, float probability, std::size_t testSize, std::size_t lastAllowedValue, boost::mt19937& randomSource)
	{
		for(std::size_t truncationLevel = 1; truncationLevel < 4; truncationLevel++)
		{
			runTest(binomialSize, probability, testSize, truncationLevel, lastAllowedValue, randomSource);
		}
	}
	int main(int argc, char **argv)
	{
		const std::size_t testSize = 10000000;
		boost::mt19937 randomSource;
		randomSource.seed(1);
		std::size_t binomialSize[] = {::TruncatedBinomialDistribution::TruncatedBinomialDistribution::acceptanceRejectionThreshold, ::TruncatedBinomialDistribution::TruncatedBinomialDistribution::acceptanceRejectionThreshold-1};
		float probability[] = {0.01f, 0.1f, 0.5f};
		std::ptrdiff_t lastAllowedValue[] = {0L, -1};
		for(int i = 0; i < 2; i++)
		{
			for(int j = 0; j < 3; j++)
			{
				for(int k = 0; k < 2; k++)
				{
					runTest(binomialSize[i], probability[j], testSize, binomialSize[i] + lastAllowedValue[k], randomSource);
				}
			}
		}
		std::cout << "All tests passed" << std::endl;
		return 0;
	}
}
int main(int argc, char **argv)
{
	return networkReliability::main(argc, argv);
}	
