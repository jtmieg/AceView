#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include "sa.common.h"

template<typename T>
void saGPUSort(T* cp, long int number_of_records, int (*cmp)(const void *va, const void *vb))
{
/*
    thrust::host_vector<T> h_vec;
    h_vec.reserve(number_of_records);
    for (unsigned int i=0;i < number_of_records;i++) {
        h_vec.push_back(cp[i]);
    }

    thrust::device_vector<T> d_vec = h_vec;

    thrust::sort(d_vec.begin(), d_vec.end(), [&](const T& a, const T& b) {
        return cmp(&a, &b) > 0;});

    thrust::copy(d_vec.begin(), d_vec.end(), h_vec.begin());

    for (unsigned int i=0;i < h_vec.size();i++) {
        cp[i] = h_vec[i];
    }
*/

    thrust::host_vector<int> h_vec;
    h_vec.reserve(number_of_records);
    for (unsigned int i=0;i < number_of_records;i++) {
        h_vec.push_back((int)cp[i]);
    }

    thrust::device_vector<int> d_vec = h_vec;

    thrust::sort(d_vec.begin(), d_vec.end(), [&](const T& a, const T& b) {
        return cmp(&a, &b) > 0;});

    thrust::copy(d_vec.begin(), d_vec.end(), h_vec.begin());

    for (unsigned int i=0;i < h_vec.size();i++) {
        cp[i] = h_vec[i];
    }

    
}


