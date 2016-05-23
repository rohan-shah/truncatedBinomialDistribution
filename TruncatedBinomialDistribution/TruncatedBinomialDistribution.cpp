#include "TruncatedBinomialDistribution.h"
#include <boost/lexical_cast.hpp>
#include "includeMPFRTruncatedBinomial.h"
#include <boost/random/binomial_distribution.hpp>
#include <boost/math/distributions/binomial.hpp>
namespace TruncatedBinomialDistribution
{
	TruncatedBinomialDistribution::TruncatedBinomialDistribution(std::ptrdiff_t n, std::ptrdiff_t firstAllowedValue, std::ptrdiff_t lastAllowedValue, mpfr_class probability)
		:firstAllowedValue(firstAllowedValue), lastAllowedValue(lastAllowedValue), n(n), probability(probability), cumulativeProbabilities(new cumulativeProbabilityType[lastAllowedValue - firstAllowedValue + 1])
	{
		if (lastAllowedValue > n) throw std::runtime_error("Input lastAllowedValue must be smaller than or equal to input n");
		
		boost::math::binomial_distribution<mpfr_class> preciseDistribution(n, probability);

		lowerLimitProb = 0;
		mpfr_class upperLimitProb = 1;
		if (firstAllowedValue > 0) lowerLimitProb = boost::math::cdf(preciseDistribution, firstAllowedValue - 1);
		if(lastAllowedValue < n) upperLimitProb = boost::math::cdf(preciseDistribution, lastAllowedValue);

		mpfr_class difference = upperLimitProb - lowerLimitProb;
		//Hopefully, if this is difference is equal to zero it's because both the upper and lower limits is equal to one. In this case, hopefully calculating using the complements of the cdfs will fix things.
		useComplement = difference == 0;
		if(useComplement)
		{
			lowerLimitProb = 1;
			upperLimitProb = 0;
			if (firstAllowedValue > 0) lowerLimitProb = boost::math::cdf(boost::math::complement(preciseDistribution, firstAllowedValue - 1));
			if(lastAllowedValue < n) upperLimitProb = boost::math::cdf(boost::math::complement(preciseDistribution, lastAllowedValue));
			difference = lowerLimitProb - upperLimitProb;
		}
		if(difference == 0)
		{
			throw std::runtime_error("Attempted to use complements of cdfs, but still couldn't compute a non-zero binomial probability");
		}
		differenceInverse = 1/difference;
		isComputed.resize(lastAllowedValue - firstAllowedValue + 1, false);
	}
	int TruncatedBinomialDistribution::operator()(boost::mt19937& randomSource) const
	{
		boost::math::binomial_distribution<mpfr_class> preciseDistribution(n, probability);
		boost::random::uniform_01<> probabilityDist;
		double randUniform = probabilityDist(randomSource);
		for(std::ptrdiff_t i = 0; i < lastAllowedValue-firstAllowedValue+1; i++)
		{
			if(!isComputed[i])
			{
				if(useComplement) cumulativeProbabilities[i] = (lowerLimitProb - boost::math::cdf(boost::math::complement(preciseDistribution, i + firstAllowedValue)))*differenceInverse;
				else cumulativeProbabilities[i] = (boost::math::cdf(preciseDistribution, i+firstAllowedValue) - lowerLimitProb)*differenceInverse;
				isComputed[i] = true;
			}
			if(randUniform <= cumulativeProbabilities[i]) return (int)(i+firstAllowedValue);
		}
		std::stringstream ss;
		ss << "Internal error: " << firstAllowedValue << " " << lastAllowedValue << std::endl;;
		ss << "Cumulative probabilities:";
		for(std::ptrdiff_t i = 0; i < lastAllowedValue-firstAllowedValue+1; i++)
		{
			ss << " " << cumulativeProbabilities[i].convert_to<double>();
		}
		throw std::runtime_error(ss.str());
	}
	TruncatedBinomialDistribution::TruncatedBinomialDistribution(TruncatedBinomialDistribution&& other)
		: firstAllowedValue(other.firstAllowedValue), lastAllowedValue(other.lastAllowedValue), n(other.n), probability(other.probability)
	{
		cumulativeProbabilities.swap(other.cumulativeProbabilities);
		isComputed.swap(other.isComputed);
		lowerLimitProb = other.lowerLimitProb;
		differenceInverse = other.differenceInverse;
		useComplement = other.useComplement; 
	}
	TruncatedBinomialDistribution::TruncatedBinomialDistribution()
	{}
	const TruncatedBinomialDistribution::cumulativeProbabilityType& TruncatedBinomialDistribution::getCumulativeProbability(int value) const
	{
		if(!isComputed[value])
		{
			boost::math::binomial_distribution<mpfr_class> preciseDistribution(n, probability);
			if(useComplement) cumulativeProbabilities[value] = (lowerLimitProb - boost::math::cdf(boost::math::complement(preciseDistribution, value + firstAllowedValue)))*differenceInverse;
			else cumulativeProbabilities[value] = (boost::math::cdf(preciseDistribution, value+firstAllowedValue) - lowerLimitProb)*differenceInverse;
			isComputed[value] = true;
		}
		return cumulativeProbabilities[value];
	}
	const mpfr_class& TruncatedBinomialDistribution::getProbability() const
	{
		return probability;
	}
	std::ptrdiff_t TruncatedBinomialDistribution::getFirstAllowedValue() const
	{
		return firstAllowedValue;
	}
	std::ptrdiff_t TruncatedBinomialDistribution::getN() const
	{
		return n;
	}
}
