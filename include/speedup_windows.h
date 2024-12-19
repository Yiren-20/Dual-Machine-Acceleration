#define _WINSOCK_DEPRECATED_NO_WARNINGS
#include <WinSock2.h>
#include <omp.h>
#include <immintrin.h>
#include <assert.h>
#include <iostream> 
#include <iomanip>
#include <windows.h>
#include <intrin.h>    //(include immintrin.h)
#include <string>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <thread>
#define MAX_THREADS 64
#undef max
#undef min
#define MyAdd "127.0.0.1" // 这里应该是服务器的ip地址
#define MyPort 1111 // 这里应该是服务器的端口号
using namespace std;

// 快排
void kuaipai(float* r, const int& l, const int& h)//r为线性表，l是表的下界，h是表的上界
{
	float x;
	int i, j;
	x = r[l - 1]; i = l; j = h;
	while (i < j)
	{
		while ((i < j) && (r[j - 1] >= x))
			j--;
		r[i - 1] = r[j - 1];
		while (i < j && r[i - 1] <= x)
			i++;
		r[j - 1] = r[i - 1];
	}
	r[i - 1] = x;
	int m1 = i - 1, m2 = i + 1;
	if (l <= m1 && m1 <= h)kuaipai(r, l, m1);
	if (l <= m2 && m2 <= h)kuaipai(r, m2, h);
}

// 朴素快排
void sort(const float data[], const int& len, float result[])
{
	for (int i = 0; i < len; i++)
	{
		result[i] = log10(sqrt(data[i]));
	}
	kuaipai(result, 1, len+1);

}

// 检查是否正确排序
int check(float result[], const int& len)
{
	float sign;
	float new_sign;

	sign = result[1] - result[0];
	for (int i = 1; i < len - 1; i++)
	{
		new_sign = result[i + 1] - result[i];
		if (new_sign * sign < 0)
		{
			cout << "sorting error！i=" << i
				<< "时，result[i]=" << result[i]
				<< ",result[i + 1]=" << result[i + 1]
				<< endl;
			return 0;
		}
	}
    cout<<"sorting successfully!"<<endl;
	return 1;
};

// 朴素求sum
float sum(const float data[], const int& len) {
    double sum = 0;
    for (int i = 0; i < len; i++) {
        sum += log(sqrt(data[i]));
    }
    return sum;
}

// 朴素求max
float max(const float data[], const int& len) {
    float max = 0;
    float temp;
    for (int i = 0; i < len; i++) {
        temp = log(sqrt(data[i]));
        if (max < temp) {
            max = temp;
        }
    }
    return max;
}

// 求sse256寄存器中8个数的和
float horizontal_sum(__m256 sum_sse) {
    float sum8[8];
    float sum1 = 0;
    _mm256_store_ps(sum8, sum_sse);
    for (int i = 0; i < 8; i++) sum1 += sum8[i];
    return sum1;
}

// 求sse256寄存器中8个数的最大值
float horizontal_max(__m256 max_sse) {
    float max8[8];
    _mm256_store_ps(max8, max_sse);
    float max1 = max8[0];
    for (int i = 1; i < 8; i++) {
        if (max1 < max8[i]) {
            max1 = max8[i];
        }
    }
    return max1;
}

// 加速sum
float sumSpeedUpOMP(const float data[], const int& len) {
    assert(len % 8 == 0 && "要求N为8的倍数");
    int nb_iters = len / 8;
    __m256 global_sum = _mm256_setzero_ps();
    __m256* ptr = (__m256*)data;
#pragma omp parallel
    {
        __m256 sum = _mm256_setzero_ps();
#pragma omp for
        for (int i = 0; i < nb_iters; ++i) {
            __m256 right = _mm256_log_ps(_mm256_sqrt_ps(ptr[i]));
            sum = _mm256_add_ps(sum, right);
        }
#pragma omp critical
        global_sum = _mm256_add_ps(sum, global_sum);
    }

    return horizontal_sum(global_sum);
}

// 加速maxOMP
float maxSpeedUpOMP(const float data[], const int& len) {
    assert(len % 8 == 0 && "要求N为8的倍数");
    int nb_iters = len / 8;
    __m256 global_max = _mm256_setzero_ps();
    __m256* ptr = (__m256*)data;
#pragma omp parallel
    {
        __m256 max = _mm256_setzero_ps();
#pragma omp for
        for (int i = 0; i < nb_iters; ++i) {
            __m256 right = _mm256_log_ps(_mm256_sqrt_ps(ptr[i]));
            max = _mm256_max_ps(max, right);
        }
#pragma omp critical
        global_max = _mm256_max_ps(max, global_max);
    }

    return horizontal_max(global_max);
}

// 加速max
float maxSpeedUp(const float data[], const int& len) {
    assert(len % 8 == 0 && "要求N为8的倍数");
    int nb_iters = len / 8;
    __m256 max = _mm256_setzero_ps();
    __m256* ptr = (__m256*)data;
    for (int i = 0; i < nb_iters; ++i, ++ptr) {
        __m256 right = _mm256_log_ps(_mm256_sqrt_ps(*ptr));
        max = _mm256_max_ps(max, right);
    }

    return horizontal_max(max);
}

// 合并两个排好序的数组
void merge(const float* src1, const int& len1, const float* src2, const int& len2,
    float* dst) {
    const float* end1 = src1 + len1;
    const float* end2 = src2 + len2;

    while (true) {
        if (*src1 < *src2) {
            *dst++ = *src1;
            if (++src1 >= end1) break;
        }
        else {
            *dst++ = *src2;
            if (++src2 >= end2) break;
        }
    }
    for (; src1 < end1; ++src1, ++dst) {
        *dst = *src1;
    }
    for (; src2 < end2; ++src2, ++dst) {
        *dst = *src2;
    }
}

// 加速sort（已废弃）
void sortSpeedUpOMP(const float data[], const int& len, float result[],
    float temparr[]) {
#pragma omp parallel shared(result, temparr)
        {
            int thread_num = omp_get_num_threads();
            int thread_id = omp_get_thread_num();

#pragma omp master
            {
                DWORD pos;
                _BitScanForward(&pos, thread_num);
                // 2, 8, 32线程时归约奇数次，用temparr做缓存
                // 1, 4, 16线程时归约偶数次，直接写到result
                if (pos % 2 != 1) {
                    std::swap(result, temparr);
                }
            }
#pragma omp flush(result, temparr)

#pragma omp for nowait
            for (int i = 0; i < len; i += 8) {
                // result[i] = log(sqrt(data[i]));
                __m256 temp = _mm256_log_ps(_mm256_sqrt_ps(_mm256_load_ps(&data[i])));
                _mm256_store_ps(&temparr[i], temp);
            }

            kuaipai(temparr, (len * thread_id / thread_num),
                (len * (thread_id + 1) / thread_num));

#pragma omp barrier
            int num_groups = thread_num / 2;
            while (num_groups != 0) {
                if (thread_id < num_groups) {
                    int section_len = len / num_groups / 2;
                    int section1 = len * thread_id / num_groups;
                    int section2 = section1 + section_len;
                    merge(temparr + section1, section_len, temparr + section2, section_len,
                        result + section1);
                }

#pragma omp barrier
                num_groups /= 2;

                if (thread_id < num_groups) {
                    int section_len = len / num_groups / 2;
                    int section1 = len * thread_id / num_groups;
                    int section2 = section1 + section_len;
                    merge(result + section1, section_len, result + section2, section_len,
                        temparr + section1);
                }

#pragma omp barrier
                num_groups /= 2;
            }
        }
}

// 加速sort
void sort_omp(const float data[], const int& len, float result[])
{
	int num = len / 64;
	float* sortResults = new float[len];
#pragma omp parallel for
	for (int i = 0; i < 64; i++)
	{
		memcpy(&result[num * i], &data[num * i], num * sizeof(float));
		kuaipai(result, num * i, num * (i + 1) - 1);
	}
	memcpy(sortResults, result, len * sizeof(float));
	for (int t = 64 / 2; t > 0; t = t / 2)
	{
		int n = len / t / 2;
#pragma omp parallel for
		for (int i = 0; i < t; i++)
			//Merge(&sortResults[n * 2 * i], &sortResults[n * (2 * i + 1)], n, &result[2 * n * i]);
            merge(&sortResults[n * 2 * i], n,&sortResults[n * (2 * i + 1)], n, &result[2 * n * i]);
	}
}

// 回收堆上空间
void deleteraw(float* rawFloatData)
{
	delete[]rawFloatData;
	rawFloatData = NULL; 
}

