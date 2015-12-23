#pragma once
#include <map>
#include <assert.h>
//#include <hash_map>
//#include <unordered_map>

template<typename Scalar, typename IndexType, int DimsCount>
struct Polynomial
{
  typedef Polynomial<Scalar, IndexType, DimsCount> PolynomialType;
  struct Pows
  {
    bool operator <(const Pows& other) const
    {
      int diffIndex = 0;
      for(;((pows[diffIndex] == other.pows[diffIndex]) && (diffIndex < DimsCount - 1)); diffIndex++);
      return pows[diffIndex] < other.pows[diffIndex];
    }
    bool operator == (const Pows& other) const
    {
      for(IndexType i = 0; i < DimsCount; i++)
      {
        if(pows[i] != other.pows[i]) 
        {
          return false;
        }
      }
      return true;
    }
    IndexType pows[DimsCount];
  };

  struct PowsHasher
  {
    const static size_t bucket_size = 10;
    const static size_t min_buckets = (1 << 3);
    PowsHasher()
    {
    }
    size_t operator()(const Pows &key) const
    {
      size_t hash_value = 1;
      for(IndexType i = 0; i < DimsCount; i++)
      {
        hash_value *= 10;
        hash_value += key.pows[i];
      }
      return hash_value;
    }

    bool operator()(const Pows &left, const Pows &right) const
    {
      return left < right;
    }    
  };

  Polynomial(Scalar value = 0)
  {
    Pows zeroPows;
    for(IndexType coordIndex = 0; coordIndex < DimsCount; coordIndex++)
    {
      zeroPows.pows[coordIndex] = 0;
    }
    AddTerm(zeroPows, value);
  }

  Polynomial &operator *=(Scalar mult)
  {
    for(typename TermsMap::iterator it = terms.begin(); it != terms.end(); it++)
    {
      it->second *= mult;
    }
    return *this;
  }

  Polynomial &operator +=(const Polynomial &other)
  {
    for(typename TermsMap::const_iterator it = other.terms.begin(); it != other.terms.end(); it++)
    {
      terms[it->first] += it->second;
    }
    return *this;
  }

  Polynomial &operator -=(Polynomial &other)
  {
    for(typename TermsMap::const_iterator it = other.terms.begin(); it != other.terms.end(); it++)
    {
      terms[it->first] -= it->second;
    }
    return *this;
  }

  Polynomial operator *(const Polynomial &other) const
  {
    PolynomialType res;
    for(typename TermsMap::const_iterator it0 = terms.begin(); it0 != terms.end(); it0++)
    {
      Pows pows0 = it0->first;
      for(typename TermsMap::const_iterator it1 = other.terms.begin(); it1 != other.terms.end(); it1++)
      {
        Pows pows1 = it1->first;
        Pows resPows;

        for(IndexType coordIndex = 0; coordIndex < DimsCount; coordIndex++)
        {
          resPows.pows[coordIndex] = pows0.pows[coordIndex] + pows1.pows[coordIndex];
        }

        if(it0->second * it1->second == 0) continue;

        if(res.terms.find(resPows) == res.terms.end())
        {
          res.terms[resPows] = it0->second * it1->second;
        }else
        {
          res.terms[resPows] += it0->second * it1->second;
        }
      }
    }
    return res;
  }

  template<typename OtherPolynomial>
  OtherPolynomial Substitute   (const OtherPolynomial *coords) const
  {
    OtherPolynomial res(0);
    for(typename TermsMap::const_iterator it = terms.begin(); it != terms.end(); it++)
    {
      if(it->second == 0) continue;

      OtherPolynomial termRes(Scalar(it->second));
      for(IndexType coordIndex = 0; coordIndex < DimsCount; coordIndex++)
      {
        for(IndexType powValue = 0; powValue < it->first.pows[coordIndex]; powValue++)
        {
          termRes = termRes * coords[coordIndex];
        }
      }

      res += termRes;
    }

    return res;
  }
  PolynomialType Differentiate(IndexType *derivativeDegrees) const
  {
    PolynomialType res;

    for(typename TermsMap::const_iterator it = terms.begin(); it != terms.end(); it++)
    {
      bool ok = 1;
      Pows termPows;
      Scalar coef = it->second;
      for(IndexType coordIndex = 0; coordIndex < DimsCount; coordIndex++)
      {
        if(derivativeDegrees[coordIndex] == 0)
        {
          termPows.pows[coordIndex] = it->first.pows[coordIndex];
        }else
        if(derivativeDegrees[coordIndex] == 1)
        {
          coef *= it->first.pows[coordIndex];
          termPows.pows[coordIndex] = it->first.pows[coordIndex] - 1;
        }else
        {
          assert(0);
        }
        if(termPows.pows[coordIndex] < 0) ok = 0;
      }
      if (ok && coef != 0)
      {
        res.AddTerm(termPows, coef);
      }
    }
    return res;
  }
  Scalar GetValue(Scalar *coords) const
  {
    Scalar res = 0;
    for(typename TermsMap::const_iterator it = terms.begin(); it != terms.end(); it++)
    {
      Scalar term(it->second);
      for(IndexType coordIndex = 0; coordIndex < DimsCount; coordIndex++)
      {
        term *= Pow(coords[coordIndex], it->first.pows[coordIndex]);
      }
      res += term;
    }
    return res;
  }

  void AddTerm(Pows pows, Scalar coef)
  {
    terms[pows] = coef;
  }

  //computes volume of a unit tetrahedra with dimsCount dimentions
  Scalar ComputeSubspaceIntegral(IndexType dimsCount) const
  {
    assert(dimsCount <= DimsCount);
    Scalar res = 0;
    for(typename TermsMap::const_iterator it = terms.begin(); it != terms.end(); it++)
    {
      Scalar nom(1.0);
      IndexType denom = dimsCount;

      for(IndexType coordIndex = 0; coordIndex < dimsCount; coordIndex++)
      {
        nom *= Factorial(it->first.pows[coordIndex]);
        denom += it->first.pows[coordIndex];
      }
      res += it->second * Scalar(nom) / Scalar(Factorial(denom));
    }
    return res;
  }
private:

  Scalar Pow(Scalar a, Scalar b) const
  {
    if(fabs(b) < Scalar(1e-5)) return Scalar(1.0);
    if(fabs(a) < Scalar(1e-5)) return 0;
    return exp(b * log(a));
  }

  Scalar Factorial(IndexType num) const
  {
    Scalar res(1.0);
    for(IndexType i = 1; i <= num; i++)
    {
      res *= Scalar(i);
    }
    return res;
  }

  typedef std::map<Pows, Scalar> TermsMap;
//  typedef stdext::hash_map<Pows, Scalar, PowsHasher> TermsMap;
//  typedef std::unordered_map<Pows, Scalar, PowsHasher> TermsMap;
  TermsMap terms;
};
