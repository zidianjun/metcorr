#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// before using
// gcc -fPIC -shared lib.c -o lib.so

float* group_by(float* arr, int N, int len, int report,
                float bin_size, float max_sep) {
    /*
    arr is [f] + [x] + [y], each of them with length N.
    return [bin_inds] + [bin_stat].
    */
    int res_len = len * 2;
    float* res = (float*) malloc(sizeof(float) * res_len);
    memset(res, 0, sizeof(float) * res_len);
    float d;
    for (int p1 = 0; p1 < N; p1++) {
        for (int p2 = 0; p2 < N; p2++) {
            d = sqrt((arr[p1+1*N] - arr[p2+1*N]) * (arr[p1+1*N] - arr[p2+1*N]) +
                     (arr[p1+2*N] - arr[p2+2*N]) * (arr[p1+2*N] - arr[p2+2*N]));
            int idx;
            if (d < max_sep) {  // For acceleration.
                if (d <= bin_size && d > 1e-6) {idx = 1;}
                else {idx = (int)(d / bin_size);}
                res[idx+len] = ((arr[p1] * arr[p2] + res[idx+len] * res[idx]) /
                                (res[idx] + 1));
                res[idx] += 1;
            }
        }
        if (report == 1){
            printf("[Two point correlation matching process] %d/%d \r", p1, N);
        }
    }

    return res;
}

