#include "speedup_windows.h"
#define SUBDATANUM 2000000
#define DATANUM (SUBDATANUM * MAX_THREADS)

void initial(float* rawFloatData)
{
	for (size_t i = 0; i < DATANUM; i++)//数据初始化
	{
		rawFloatData[i] = float(i + 1);
	}
}

void initial_rand(float* rawFloatData)
{
	for (size_t i = 0; i < DATANUM; i++)//数据初始化
	{
		rawFloatData[i] = rand() % DATANUM + 2; //数值不大于DATANUM的随机数,+1为了防止0出现
	}
}

int main()
{
    LARGE_INTEGER start;
    LARGE_INTEGER end;
    float summary = 0;     // Result of summation
    float maxi = 0;        // Result of finding the maximum value
    float sum_time = 0.0f; // Time taken for summation without acceleration
    float max_time = 0.0f; // Time taken for finding the maximum without acceleration
    float sort_time = 0.0f; // Time taken for quicksort without acceleration
    int rightOrNot1 = 0;   // Check for sorting correctness

    float* rawFloatData = new float[DATANUM];     // Allocate space for the data array
    float* finalFloatData = new float[DATANUM];   // Allocate space for the sorted result array
    float* tempArr = new float[DATANUM];

    initial(rawFloatData);   // Initialize the data

    std::cout << "Summation part: " << std::endl;
    QueryPerformanceCounter(&start); // Start
    summary = sum(rawFloatData, DATANUM);
    QueryPerformanceCounter(&end); // End
    std::cout << std::setprecision(4) << "Without acceleration version, output sum result: " << summary << std::endl;
    std::cout << "Time Consumed: " << (end.QuadPart - start.QuadPart) << std::endl;
    sum_time = end.QuadPart - start.QuadPart;

    summary = 0;
    QueryPerformanceCounter(&start); // Start
#pragma omp parallel for reduction(+:summary)
    for (long i = 0; i < MAX_THREADS; i++)
        summary = summary + sum(rawFloatData + i * SUBDATANUM, SUBDATANUM);
    QueryPerformanceCounter(&end); // End
    std::cout << std::setprecision(4) << "OpenMP accelerated version, output sum result: " << summary << std::endl;
    std::cout << "Time Consumed: " << (end.QuadPart - start.QuadPart) << "     " << "Acceleration Ratio: " << sum_time / (end.QuadPart - start.QuadPart) << std::endl;

    QueryPerformanceCounter(&start); // Start
    summary = 0;
#pragma omp parallel for reduction(+:summary)
    for (long i = 0; i < MAX_THREADS; i++)
        summary = summary + sumSpeedUpOMP(rawFloatData + i * SUBDATANUM, SUBDATANUM);
    QueryPerformanceCounter(&end); // End
    std::cout << std::setprecision(4) << "Multiple SSE + OpenMP accelerated version, output sum result: " << summary << std::endl;
    std::cout << "Time Consumed: " << (end.QuadPart - start.QuadPart) << "     " << "Acceleration Ratio: " << sum_time / (end.QuadPart - start.QuadPart) << std::endl;

    std::cout << "-----------------------------------------------------------------------------------" << std::endl;
    std::cout << "Maximum value part: " << std::endl;

    QueryPerformanceCounter(&start); // Start
    maxi = max(rawFloatData, DATANUM);
    std::cout << "Output maximum value result (without acceleration): " << maxi << std::endl;
    QueryPerformanceCounter(&end); // End
    std::cout << "Time Consumed: " << (end.QuadPart - start.QuadPart) << std::endl;
    max_time = (end.QuadPart - start.QuadPart);

    QueryPerformanceCounter(&start); // Start
    maxi = maxSpeedUp(rawFloatData, DATANUM);
    QueryPerformanceCounter(&end); // End
    std::cout << "Output maximum value result (SSE accelerated version): " << maxi << std::endl;
    std::cout << "Time Consumed: " << (end.QuadPart - start.QuadPart) << "     " << "Acceleration Ratio: " << max_time / (end.QuadPart - start.QuadPart) << std::endl;

    QueryPerformanceCounter(&start); // Start
    maxi = 0;
#pragma omp parallel for 
    for (long i = 0; i < MAX_THREADS; i++)
    {
        float tmp = maxSpeedUpOMP(rawFloatData + i * SUBDATANUM, SUBDATANUM);
        maxi = (maxi > tmp ? maxi : tmp);
    }
    QueryPerformanceCounter(&end); // End
    std::cout << "Output maximum value result (multiple SSE + OpenMP accelerated version): " << maxi << std::endl;
    std::cout << "Time Consumed: " << (end.QuadPart - start.QuadPart) << "     " << "Acceleration Ratio: " << max_time / (end.QuadPart - start.QuadPart) << std::endl;

    initial_rand(rawFloatData);
    std::cout << "-----------------------------------------------------------------------------------" << std::endl;
    std::cout << "Sorting part: " << std::endl;

    // Quick Sort
    QueryPerformanceCounter(&start); // Start
    sort(rawFloatData, DATANUM, finalFloatData);
    QueryPerformanceCounter(&end); // End
    rightOrNot1 = check(finalFloatData, DATANUM);
    std::cout << "No speedup sort: " << std::endl;
    std::cout << "Time Consumed: " << (end.QuadPart - start.QuadPart) << std::endl;
    sort_time = (end.QuadPart - start.QuadPart);

    memset(finalFloatData, -1, sizeof(finalFloatData));
    memset(tempArr, 0, sizeof(tempArr));
    QueryPerformanceCounter(&start); // Start
    //sortSpeedUpOMP(rawFloatData, DATANUM, finalFloatData, tempArr);
    sort_omp(rawFloatData, DATANUM, finalFloatData);
    QueryPerformanceCounter(&end); // End
    rightOrNot1 = check(finalFloatData, DATANUM);
    std::cout << "Speedup sort: " << std::endl;
    std::cout << "Time Consumed: " << (end.QuadPart - start.QuadPart) << "     " << "Acceleration Ratio: " << sort_time / (end.QuadPart - start.QuadPart) << std::endl;
    system("pause");
    delete[] rawFloatData;
    delete[] finalFloatData;
    delete[] tempArr;
}

