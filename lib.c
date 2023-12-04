#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// before using
// gcc -fPIC -shared lib.c -o lib.so

float* group_by(float* arr, int N, int len, int report,
                float min_sep, float max_sep, float bin_size) {
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
            float dx = arr[p1+1*N] - arr[p2+1*N];
            float dy = arr[p1+2*N] - arr[p2+2*N];
            float d = sqrt(dx * dx + dy * dy);
            if ((d >= min_sep) && (d < max_sep)) {
                int idx = (int)((d - min_sep) / bin_size);
                res[idx + len] = (
                    (arr[p1] * arr[p2] + res[idx + len] * res[idx]) / (res[idx] + 1));
                res[idx] += 1;
            }
        }
        if (report == 1){
            printf("[Two point correlation matching process] %d/%d \r", p1, N);
        }
    }

    return res;
}



float* group(float* arr, int N, int len_rad, int len_azi, int report,
             float min_sep, float max_sep, float bin_size,
             float min_pa,  float max_pa,  float azi_size) {
    /*
    arr is [f] + [x] + [y], each of them with length N.
    return [bin_inds] + [bin_stat].
    */
    int res_len = len_rad * len_azi * 2;
    float* res = (float*) malloc(sizeof(float) * res_len);
    memset(res, 0, sizeof(float) * res_len);
    for (int p1 = 0; p1 < N; p1++) {
        for (int p2 = 0; p2 < N; p2++) {
            float dx = arr[p1+1*N] - arr[p2+1*N];
            float dy = arr[p1+2*N] - arr[p2+2*N];
            float d = sqrt(dx * dx + dy * dy);
            float pa;
            if (fabsf(dx) < 1e-6) {pa = 0;}
            else {pa = (atan(dy / dx) + asin(1)) / asin(1) * 90;}
            if ((d >= min_sep) && (d < max_sep)) {
                int idx_rad = (int)((d - min_sep) / bin_size);
                int idx_azi = (int)((pa - min_pa) / azi_size);
                int idx = idx_azi * len_rad + idx_rad;
                res[idx + len_rad*len_azi] = (
                    (arr[p1] * arr[p2] + res[idx + len_rad * len_azi] * res[idx]) /
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




