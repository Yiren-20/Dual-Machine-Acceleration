#include "speedup_linux.h"
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
	float summary = 0;  //求和的结果
	float maxi = 0;     //求最大值的结果
	double sum_time = 0.0f;    //不加速求和用的时间
	double max_time = 0.0f;    //不加速求最大值用的时间
	double sort_time = 0.0f;   //不加速快速排序的时间
	float origin_max = 0.0f;
	int rightOrNot1 = 0;//检查排序

	float* rawFloatData = new float[DATANUM];  //开辟数组空间
	float* finalFloatData = new float[DATANUM];//开辟排序结果数组空间
	float* tempArr = new float[DATANUM];
	initial(rawFloatData);  	//数据初始化

	std::cout << "Summation part: " << std::endl;
	timeval start, end;
	gettimeofday(&start, NULL);//sum start 
	summary = sum(rawFloatData, DATANUM);
	gettimeofday(&end, NULL);
	std::cout << std::setprecision(4) << "Without acceleration version, output sum result: " << summary << std::endl;
	sum_time = 1e3 * (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1e3;
	cout << "Sum Time Consumed:" << sum_time << "ms" << endl;

	summary = 0;
	gettimeofday(&start, NULL);//sum start //start  
#pragma omp parallel for reduction(+:summary)
	for (long i = 0; i < MAX_THREADS; i++)
		summary = summary + sum(rawFloatData + i * SUBDATANUM, SUBDATANUM);

	gettimeofday(&end, NULL);
	std::cout << std::setprecision(4) << "OpenMP accelerated version, output sum result: " << summary << std::endl;
	cout << "Sum Time Consumed:" << 1e3 * (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1e3 << "ms" << endl;


	gettimeofday(&start, NULL);//sum start //start
	summary = 0;
#pragma omp parallel for reduction(+:summary)
	for (long i = 0; i < MAX_THREADS; i++)
		summary = summary + sumSpeedUpOMP(rawFloatData + i * SUBDATANUM, SUBDATANUM);
	gettimeofday(&end, NULL);
	std::cout << std::setprecision(4) << "Multiple SSE + OpenMP accelerated version, output sum result: " << summary << std::endl;
	cout << "Sum Time Consumed:" << 1e3 * (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1e3 << "ms. Acceleration Ratio: " <<sum_time/(1e3 * (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1e3)<< endl;
	// std::cout << "Time Ratio:" << sumtime1/sumtime2 << std::endl;

	cout << "-----------------------------------------------------------------------------------" << endl;
	std::cout << "Maximum value part: " << std::endl;

	gettimeofday(&start, NULL);//sum start //start  
	maxi = max(rawFloatData, DATANUM);
	std::cout << "Output maximum value result (without acceleration): " << maxi << std::endl;
	gettimeofday(&end, NULL);
	max_time = 1e3 * (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1e3;
	cout << "Max Time Consumed:" << 1e3 * (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1e3 << "ms"<< endl;

	gettimeofday(&start, NULL);//sum start //start  
	maxi = maxSpeedUp(rawFloatData, DATANUM);
	gettimeofday(&end, NULL);
	std::cout << "Output maximum value result (SSE accelerated version): " << maxi << std::endl;
	cout << "Max Time Consumed:" << 1e3 * (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1e3 << "ms" << endl;

	gettimeofday(&start, NULL);//sum start //start
	maxi = 0;
#pragma omp parallel for 
	for (long i = 0; i < MAX_THREADS; i++)
	{
		float tmp = maxSpeedUpOMP(rawFloatData + i * SUBDATANUM, SUBDATANUM);
		maxi = (maxi > tmp ? maxi : tmp);
	}
	gettimeofday(&end, NULL);
	std::cout << "Output maximum value result (multiple SSE + OpenMP accelerated version): " << maxi << std::endl;
	cout << "Max Time Consumed:" << 1e3 * (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1e3 << "ms. Acceleration Ratio: " <<max_time/(1e3 * (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1e3)<< endl;
	
	cout << "-----------------------------------------------------------------------------------" << endl;
	std::cout << "Sorting part: " << std::endl;
	initial_rand(rawFloatData);
	//快排
	gettimeofday(&start, NULL);//sum start //start  
	sort(rawFloatData, DATANUM, finalFloatData);
	gettimeofday(&end, NULL);
	rightOrNot1 = check(finalFloatData, DATANUM);
	std::cout << "No speedup sort: " << std::endl;
	sort_time = 1e3 * (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1e3;
	cout << "Sort Time Consumed:" << 1e3 * (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1e3 << "ms" << endl;

	gettimeofday(&start, NULL);//sum start //start 
	// sortSpeedUpOMP(rawFloatData, DATANUM, finalFloatData, tempArr);
	sort_omp(rawFloatData, DATANUM, finalFloatData);
	gettimeofday(&end, NULL);
	rightOrNot1 = check(finalFloatData, DATANUM);
	std::cout << "Speedup sort: " << std::endl;
	cout << "Sort Time Consumed:" << 1e3 * (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1e3 << "ms. Acceleration Ratio: "<<sort_time/(1e3 * (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1e3) << endl;

	deleteraw(rawFloatData);
	deleteraw(finalFloatData);
	deleteraw(tempArr);

}