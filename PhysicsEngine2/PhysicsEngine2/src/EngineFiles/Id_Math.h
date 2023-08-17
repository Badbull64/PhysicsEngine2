#pragma once
#ifndef ID_MAHT_H
#include <iostream>
#include <cmath>
const float tempVal = std::sqrtf(2);
const float PI_MT = 3.1415926f;
float myErfInv2(float x) {
    float tt1, tt2, lnx, sgn;
    sgn = (x < 0) ? -1.0f : 1.0f;

    x = (1 - x) * (1 + x);        // x = 1 - x*x;
    lnx = logf(x);

    tt1 = 2 / (PI_MT * 0.147) + 0.5f * lnx;
    tt2 = 1 / (0.147) * lnx;

    return(sgn * sqrtf(-tt1 + sqrtf(tt1 * tt1 - tt2)));
}

float calc_Var(const std::vector<float>& data) {
    if (data.empty()) {
        std::cerr << "Error: Empty data vector in calculateVarianceOnline()" << std::endl;
        return 0.0f;  // Return 0 as a default value for an empty dataset (you may handle this differently)
    }

    float mean = data[0];  // Initial estimate of the mean
    float sumSquaredDifferences = 0.0f;
    for (size_t i = 1; i < data.size(); ++i) {
        float diff = data[i] - mean;
        mean += diff / static_cast<float>(i + 1);
        sumSquaredDifferences += diff * (data[i] - mean);
    }

    float variance = sumSquaredDifferences / static_cast<float>(data.size() - 1);
    return variance;
}



float inv_normal(float p) {
    return tempVal * myErfInv2(2.0f*p-1.0f);
}
#endif // !ID_MAHT_H

