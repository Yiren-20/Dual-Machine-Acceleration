#include <iostream> 
#include <iomanip>
#include <immintrin.h>
#include <omp.h>
#include <stdio.h>
#include <cstring>
#include <math.h>
#include <cfloat>
#include <chrono>
#include <stdint.h>
#include <assert.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <cstdlib>
#include <ctime>
#include <time.h>
#include<sys/types.h>
#include<stdlib.h>
#include<sys/time.h>
#undef max

using namespace std;
#define __SSE2_MATH__ 1
#define __SSE2__ 1
#define MAX_THREADS 64 

#define MyAdd "127.0.0.1" // 这里应该是服务器的ip地址
#define MyPort 1111 // 这里应该是服务器的端口号

inline __m128 _mm_log_ps(__m128 x)
{
	static const int _ps_min_norm_pos[4] __attribute__((aligned(16))) = { 0x00800000, 0x00800000, 0x00800000, 0x00800000 };
	static const int _ps_inv_mant_mask[4] __attribute__((aligned(16))) = { ~0x7f800000, ~0x7f800000, ~0x7f800000, ~0x7f800000 };
	static const int _pi32_0x7f[4] __attribute__((aligned(16))) = { 0x7f, 0x7f, 0x7f, 0x7f };
	static const float _ps_1[4] __attribute__((aligned(16))) = { 1.0f, 1.0f, 1.0f, 1.0f };
	static const float _ps_0p5[4] __attribute__((aligned(16))) = { 0.5f, 0.5f, 0.5f, 0.5f };
	static const float _ps_sqrthf[4] __attribute__((aligned(16))) = { 0.707106781186547524f, 0.707106781186547524f, 0.707106781186547524f, 0.707106781186547524f };
	static const float _ps_log_p0[4] __attribute__((aligned(16))) = { 7.0376836292E-2f, 7.0376836292E-2f, 7.0376836292E-2f, 7.0376836292E-2f };
	static const float _ps_log_p1[4] __attribute__((aligned(16))) = { -1.1514610310E-1f, -1.1514610310E-1f, -1.1514610310E-1f, -1.1514610310E-1f };
	static const float _ps_log_p2[4] __attribute__((aligned(16))) = { 1.1676998740E-1f, 1.1676998740E-1f, 1.1676998740E-1f, 1.1676998740E-1f };
	static const float _ps_log_p3[4] __attribute__((aligned(16))) = { -1.2420140846E-1f, -1.2420140846E-1f, -1.2420140846E-1f, -1.2420140846E-1f };
	static const float _ps_log_p4[4] __attribute__((aligned(16))) = { 1.4249322787E-1f, 1.4249322787E-1f, 1.4249322787E-1f, 1.4249322787E-1f };
	static const float _ps_log_p5[4] __attribute__((aligned(16))) = { -1.6668057665E-1f, -1.6668057665E-1f, -1.6668057665E-1f, -1.6668057665E-1f };
	static const float _ps_log_p6[4] __attribute__((aligned(16))) = { 2.0000714765E-1f, 2.0000714765E-1f, 2.0000714765E-1f, 2.0000714765E-1f };
	static const float _ps_log_p7[4] __attribute__((aligned(16))) = { -2.4999993993E-1f, -2.4999993993E-1f, -2.4999993993E-1f, -2.4999993993E-1f };
	static const float _ps_log_p8[4] __attribute__((aligned(16))) = { 3.3333331174E-1f, 3.3333331174E-1f, 3.3333331174E-1f, 3.3333331174E-1f };
	static const float _ps_log_q1[4] __attribute__((aligned(16))) = { -2.12194440e-4f, -2.12194440e-4f, -2.12194440e-4f, -2.12194440e-4f };
	static const float _ps_log_q2[4] __attribute__((aligned(16))) = { 0.693359375f, 0.693359375f, 0.693359375f, 0.693359375f };

	__m128 one = *(__m128*)_ps_1;
	__m128 invalid_mask = _mm_cmple_ps(x, _mm_setzero_ps());
	/* cut off denormalized stuff */
	x = _mm_max_ps(x, *(__m128*)_ps_min_norm_pos);
	__m128i emm0 = _mm_srli_epi32(_mm_castps_si128(x), 23);

	/* keep only the fractional part */
	x = _mm_and_ps(x, *(__m128*)_ps_inv_mant_mask);
	x = _mm_or_ps(x, _mm_set1_ps(0.5f));

	emm0 = _mm_sub_epi32(emm0, *(__m128i *)_pi32_0x7f);
	__m128 e = _mm_cvtepi32_ps(emm0);
	e = _mm_add_ps(e, one);

	__m128 mask = _mm_cmplt_ps(x, *(__m128*)_ps_sqrthf);
	__m128 tmp = _mm_and_ps(x, mask);
	x = _mm_sub_ps(x, one);
	e = _mm_sub_ps(e, _mm_and_ps(one, mask));
	x = _mm_add_ps(x, tmp);

	__m128 z = _mm_mul_ps(x, x);
	__m128 y = *(__m128*)_ps_log_p0;
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, *(__m128*)_ps_log_p1);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, *(__m128*)_ps_log_p2);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, *(__m128*)_ps_log_p3);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, *(__m128*)_ps_log_p4);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, *(__m128*)_ps_log_p5);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, *(__m128*)_ps_log_p6);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, *(__m128*)_ps_log_p7);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, *(__m128*)_ps_log_p8);
	y = _mm_mul_ps(y, x);

	y = _mm_mul_ps(y, z);
	tmp = _mm_mul_ps(e, *(__m128*)_ps_log_q1);
	y = _mm_add_ps(y, tmp);
	tmp = _mm_mul_ps(z, *(__m128*)_ps_0p5);
	y = _mm_sub_ps(y, tmp);
	tmp = _mm_mul_ps(e, *(__m128*)_ps_log_q2);
	x = _mm_add_ps(x, y);
	x = _mm_add_ps(x, tmp);
	x = _mm_or_ps(x, invalid_mask); // negative arg will be NAN

	return x;
}

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
			cout << "排序有误！i=" << i
				<< "时，result[i]=" << result[i]
				<< ",result[i + 1]=" << result[i + 1]
				<< endl;
			return 0;
		}
	}
	return 1;
};

// 朴素求num
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

// 求sse128寄存器中4个数的和
float horizontal_sum(__m128 sum_sse) {
    float sum4[4];
    float sum1 = 0;
    _mm_store_ps(sum4, sum_sse);
    for (int i = 0; i < 4; i++) sum1 += sum4[i];
    return sum1;
}

// 求sse128寄存器中4个数的最大值
float horizontal_max(__m128 max_sse) {
    float max4[4];
    _mm_store_ps(max4, max_sse);
    float max1 = max4[0];
    for (int i = 1; i < 4; i++) {
        if (max1 < max4[i]) {
            max1 = max4[i];
        }
    }
    return max1;
}

// 加速sum
float sumSpeedUpOMP(const float data[], const int& len) {
    assert(len % 4 == 0 && "要求N为8的倍数");
    int nb_iters = len / 4;
    __m128 global_sum = _mm_setzero_ps();
    __m128* ptr = (__m128*)data;
#pragma omp parallel
    {
        __m128 sum = _mm_setzero_ps();
#pragma omp for
        for (int i = 0; i < nb_iters; ++i) {
            __m128 right = _mm_log_ps(_mm_sqrt_ps(ptr[i]));
            sum = _mm_add_ps(sum, right);
        }
#pragma omp critical
        global_sum = _mm_add_ps(sum, global_sum);
    }

    return horizontal_sum(global_sum);
}

float maxSpeedUp(const float data[], const int& len) {
    assert(len % 4 == 0 && "要求N为8的倍数");
    int nb_iters = len / 4;
    __m128 max = _mm_setzero_ps();
    __m128* ptr = (__m128*)data;
    for (int i = 0; i < nb_iters; ++i, ++ptr) {
        __m128 right = _mm_log_ps(_mm_sqrt_ps(*ptr));
        max = _mm_max_ps(max, right);
    }

    return horizontal_max(max);
}

// sse加速max
float maxSpeedUpOMP(const float data[], const int& len)
{
	assert(len % 4 == 0 && "要求N为8的倍数");
    int nb_iters = len / 4;
    __m128 global_max = _mm_setzero_ps();
    __m128* ptr = (__m128*)data;
#pragma omp parallel
    {
        __m128 max = _mm_setzero_ps();
#pragma omp for
        for (int i = 0; i < nb_iters; ++i) {
            __m128 right = _mm_log_ps(_mm_sqrt_ps(ptr[i]));
            max = _mm_max_ps(max, right);
        }
#pragma omp critical
        global_max = _mm_max_ps(max, global_max);
    }

    return horizontal_max(global_max);
}

// merge
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
void sortSpeedUpOMP(const float data[], const int & len, float result[],
    float temparr[]) {
#pragma omp parallel shared(result, temparr)
        {
            int thread_num = omp_get_num_threads();
            int thread_id = omp_get_thread_num();

#pragma omp master
            {
                uint32_t pos;
                pos = __builtin_ctz(thread_num);
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
                __m128 temp = _mm_log_ps(_mm_sqrt_ps(_mm_load_ps(&data[i])));
                _mm_store_ps(&temparr[i], temp);
            }

            kuaipai(temparr, len * thread_id / thread_num+1,
                len * (thread_id + 1) / thread_num+1);

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

// 回收堆上空间
void deleteraw(float* rawFloatData)
{
	delete[]rawFloatData;
	rawFloatData = NULL;
}

// merge
void merge_array(float arr1[], float arr2[], const int& len, float result[])
{
	int i, j, k;
	i = j = k = 0;
	while (i < len && j < len)
	{
		if (arr1[i] < arr2[j])
			result[k++] = arr1[i++];
		else
			result[k++] = arr1[j++];
	}
	while (i < len)
		result[k++] = arr1[i++];
	while (j < len)
		result[k++] = arr2[j++];
}

// 加速sort
void sort_omp(const float data[], const int& len, float result[])
{
	int num = len / 64;
	float* sortResults = new float[len];
	auto p = sortResults, q = result;
#pragma omp parallel for
	for (int i = 0; i < 64; i++)
	{
		memcpy(&result[num * i], &data[num * i], num * sizeof(float));
		kuaipai(result, num * i, num * (i + 1) - 1);
	}
	// memcpy(sortResults, result, len * sizeof(float));
	for (int t = 64 / 2; t > 0; t = t / 2)
	{
		int n = len / t / 2;
#pragma omp parallel for
		for (int i = 0; i < t; i++)
			// Merge(&sortResults[n * 2 * i], &sortResults[n * (2 * i + 1)], n, &result[2 * n * i]);
			merge_array(&p[n * 2 * i], &p[n * (2 * i + 1)], n, &q[2 * n * i]);
			//merge(&p[n * 2 * i],n,&p[n * (2 * i + 1)],n,&q[2 * n * i]);
		//memcpy(sortResults, result, len * sizeof(float));
		auto tmp = p ; p = q ; q = tmp;
	}
	memcpy(sortResults, result, len * sizeof(float));
	delete[] sortResults;
}

// merge
void merge_one(float data[], const int& left, const int& mid, const int& right,float result[])
{
	int i = left, j = mid + 1, k = left;
	while (i <= mid && j <= right)
	{
		if (data[i] > data[j])
			result[k++] = data[i++];
		else
			result[k++] = data[j++];
	}
	while (i <= mid)
		result[k++] = data[i++];
	while (j <= right)
		result[k++] = data[j++];

	for (int i = left; i <= right; i++)
		data[i] = result[i];
}

// merge
void merge_big(float data[], const int& len,const int& sortedsize, float result[])
{
	int left = 0, mid = 0, right = 0;
	// 开始合并子序列,子序列长度为1,2,4,8....成倍增长(初始子序列长度必须为1,思考为什么).
	for (int i = sortedsize; i < len; i *= 2)
	{
		for (int j = 0; j < len; j += 2 * i)
		{// 划分子序列.
			left = j;
			mid = (left + i - 1) < (len) ? (left + i - 1) : len - 1;
			right = (mid + i) < (len) ? (mid + i) : len - 1;
			merge_one(data, left, mid, right, result);
		}
	}
}