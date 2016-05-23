#ifndef TRUNCATED_BINOMIAL_DISTRIBUTION
#define TRUNCATED_BINOMIAL_DISTRIBUTION
#include <boost/scoped_array.hpp>
#include <boost/noncopyable.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <map>
#include "includeMPFRTruncatedBinomial.h"
namespace TruncatedBinomialDistribution
{
	class TruncatedBinomialDistribution : private boost::noncopyable
	{
	public:
		typedef mpfr_class cumulativeProbabilityType;
		friend struct TruncatedBinomialDistributionCollection;
		/*
		 * Truncation level is the highest level that's NOT allowed to occur. So a truncationLevel of 0 means we have a 
		 * Zero-truncated distribution
		*/
		TruncatedBinomialDistribution(std::ptrdiff_t n, std::ptrdiff_t firstAllowedValue, std::ptrdiff_t lastAllowedValue, mpfr_class probability);
		TruncatedBinomialDistribution(TruncatedBinomialDistribution&& other);
		TruncatedBinomialDistribution& operator=(TruncatedBinomialDistribution&& other);
		void swap(TruncatedBinomialDistribution&& other);
		int operator()(boost::mt19937& randomSource) const;
		const cumulativeProbabilityType& getCumulativeProbability(int value) const;
		const mpfr_class& getProbability() const;
		std::ptrdiff_t getFirstAllowedValue() const;
		std::ptrdiff_t getLastAllowedValue() const;
		std::ptrdiff_t getN() const;
		struct key
		{
			std::ptrdiff_t firstAllowedValue, lastAllowedValue;
			std::ptrdiff_t n;
		};
		key getKey()
		{
			key ret;
			ret.firstAllowedValue = firstAllowedValue;
			ret.lastAllowedValue = lastAllowedValue;
			ret.n = n;
			return ret;
		};
		struct sorter
		{
			bool operator()(const key& first, const key& second)
			{
				if(first.firstAllowedValue == second.firstAllowedValue)
				{
					if(first.lastAllowedValue == second.lastAllowedValue)
					{
						return first.n < second.n;
					}
					return first.lastAllowedValue < second.lastAllowedValue;
				}
				return first.firstAllowedValue < second.firstAllowedValue;
			}
		};
		const static std::size_t acceptanceRejectionThreshold = 10;
		//Required for some std collections. Shouldn't really be used. 
		//TruncatedBinomialDistribution(const TruncatedBinomialDistribution& other);
	private:
		//only for use by serialization code
		TruncatedBinomialDistribution();
		//The first allowed value for the binomial. So for a zero-truncated binomial it would be 1.
		std::ptrdiff_t firstAllowedValue, lastAllowedValue;
		//The number of Bernoulli trials.
		std::ptrdiff_t n;
		mpfr_class probability;
		//Should we use the complements of the CDF?
		bool useComplement;
		//The inverse of the probability that we're conditioning by
		mpfr_class differenceInverse;
		//The probability of the lower limit
		mpfr_class lowerLimitProb;
		//bit-set telling us whether we've computed that CDF value
		mutable std::vector<bool> isComputed;

		//For the exact method: These two go together...
		boost::scoped_array<cumulativeProbabilityType> cumulativeProbabilities;
	};
	struct TruncatedBinomialDistributionCollection
	{
		typedef std::map<TruncatedBinomialDistribution::key, TruncatedBinomialDistribution, TruncatedBinomialDistribution::sorter> mapType;
		mapType data;
		void swap(TruncatedBinomialDistributionCollection& other)
		{
			data.swap(other.data);
		}
	};
}
#endif
