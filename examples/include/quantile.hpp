#ifndef _QUANTILE_H_
#define _QUANTILE_H_

#include <algorithm>
#include <cmath>
#include <vector>

template<typename T>
static inline double Lerp(T v0, T v1, T t)
{
    return (1 - t)*v0 + t*v1;
}

void Quantile_print(const std::vector<double>& inData, double scale=1, bool inv=false)
{
    printf("[");
    for(size_t i=0; i<inData.size()-1; ++i)
    {
        double val=inData[i];
        if(inv)
        {
            val=1/val;
        }
        printf("%f, ", scale*val);
    }
    double val=inData[inData.size()-1];
    if(inv)
    {
        val=1/val;
    }
    printf("%f]", scale*val);
}

template<typename T>
static inline std::vector<T> Quantile(std::vector<T>& inData, const std::vector<T>& probs)
{
    if (inData.empty())
    {
        return std::vector<T>();
    }

    if (1 == inData.size())
    {
        return std::vector<T>(1, inData[0]);
    }
    //duplicate elements, else it would seg fault
    while(inData.size() < 5)
    {
        std::vector<T> dupData = inData;
        inData.insert(inData.end(), dupData.begin(), dupData.end());
    }
    std::vector<T> data = inData;
    std::sort(data.begin(), data.end());
    std::vector<T> quantiles;

    for (size_t i = 0; i < probs.size(); ++i)
    {
        T poi = Lerp<T>(-0.5, data.size() - 0.5, probs[i]);

        size_t left = std::max(int64_t(std::floor(poi)), int64_t(0));
        size_t right = std::min(int64_t(std::ceil(poi)), int64_t(data.size() - 1));

        T datLeft = data.at(left);
        T datRight = data.at(right);

        T quantile = Lerp<T>(datLeft, datRight, poi - left);

        quantiles.push_back(quantile);
    }

    return quantiles;
}

#endif
