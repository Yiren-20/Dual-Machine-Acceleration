#pragma comment(lib,"ws2_32.lib")

#include "speedup_windows.h"
#define SUBDATANUM 1000000
#define DATANUM  (SUBDATANUM * MAX_THREADS)   /*这个数值是总数据量*/

void initial(float* rawFloatData)
{
	for (size_t i = 0; i < DATANUM; i++)//数据初始化
	{
		rawFloatData[i] = float(i + 1 + DATANUM);
	}
}

void initial_rand(float* rawFloatData)
{
	for (size_t i = 0; i < DATANUM; i++)//数据初始化
	{
		rawFloatData[i] = rand() % DATANUM + 1; //数值不大于DATANUM的随机数,+1为了防止0出现
	}
}

//两边比较大小
void Compare(float result_whx[], float result_lyr[], float result_final[], int len_lyr)
{
	float* head_main = new float;
	float* head_client = new float;
	head_main = result_lyr;
	head_client = result_whx;
	int count_whx = 0;
	int count_lyr = 0;
	int lyr_or_whx = 0;//谁的数组里的数更小？；lyr-0；whx-1；
	for (int i = 0; i < 2 * DATANUM; i++)
	{
		float minFloat = FLT_MAX;
		if (*head_main < *head_client)
		{
			minFloat = *head_main;
			lyr_or_whx = 0;
		}
		else
		{
			minFloat = *head_client;
			lyr_or_whx = 1;
		}
		if (count_whx < DATANUM - 1 && lyr_or_whx == 1)
		{
			head_client++;
			count_whx++;
		}
		if (count_lyr < len_lyr - 1 && lyr_or_whx == 0)
		{
			head_main++;
			count_lyr++;
		}
		result_final[i] = minFloat;
	}
}

int main()
{
    LARGE_INTEGER start;
    LARGE_INTEGER end;
    LARGE_INTEGER start1; // Communication time
    LARGE_INTEGER end1;

    cout << "---------------server here----------------" << endl;
    WSAData wsaData;
    WORD DllVersion = MAKEWORD(2, 1);
    if (WSAStartup(DllVersion, &wsaData) != 0)
    {
        MessageBoxA(NULL, "WinSock startup error", "Error", MB_OK | MB_ICONERROR);
        return 0;
    }
    cout << "Waiting for connection......" << endl;
    SOCKADDR_IN addr;
    int addrlen = sizeof(addr);
    addr.sin_addr.s_addr = inet_addr(MyAdd); //target PC (replace "MyAdd" with your actual address)
    addr.sin_port = htons(MyPort); // server Port
    addr.sin_family = AF_INET; // IPv4 Socket

    SOCKET sListen = socket(AF_INET, SOCK_STREAM, NULL);
    bind(sListen, (SOCKADDR*)&addr, sizeof(addr));
    listen(sListen, SOMAXCONN);

    float* rawFloatData_main = new float[DATANUM];
    float* rawFloatData_me_final = new float[DATANUM];
    float* rawFloatData_you = new float[DATANUM];
    float* rawFloatData_final = new float[2 * DATANUM];
    // Data initialization
    for (size_t i = 0; i < DATANUM; i++)
    {
        rawFloatData_main[i] = float(i + 1 + DATANUM);
    }
    float summary_main = 0.0f;
    float summary_client = 0.0f;
    float max_main = 0.0f;
    float max_client = 0.0f;
    int receivesuccess = 0;
    int receiveonce = 1000;
    int received = 0;

    SOCKET newConnection; // Build a new socket for a new connection. The sListen is just listening and not used to exchange data.
    newConnection = accept(sListen, (SOCKADDR*)&addr, &addrlen); // newConnection is used to exchange data with the client
    if (newConnection == 0)
    {
        std::cout << "Bad connection." << std::endl;
    }
    else
    {
        std::cout << "Good connection." << std::endl;
        // Calculate sum
        QueryPerformanceCounter(&start);
        summary_main = 0;
#pragma omp parallel for
        for (long i = 0; i < MAX_THREADS; i++)
            summary_main = summary_main + sumSpeedUpOMP(rawFloatData_main + i * SUBDATANUM, SUBDATANUM);
        //Sleep(1);
        recv(newConnection, (char*)&summary_client, sizeof(float), NULL);
        summary_main = summary_main + summary_client;
        QueryPerformanceCounter(&end);
        cout << "Sum size=";
        cout << summary_main << std::endl;
        cout << "Time spent=";
        cout << (end.QuadPart - start.QuadPart) << std::endl;

        // Compare max
        QueryPerformanceCounter(&start);
        max_main = 0;
#pragma omp parallel for 
        for (long i = 0; i < MAX_THREADS; i++)
        {
            float tmp = maxSpeedUpOMP(rawFloatData_main + i * SUBDATANUM, SUBDATANUM);
            max_main = (max_main > tmp ? max_main : tmp);
        }
        recv(newConnection, (char*)&max_client, sizeof(float), NULL);
        if (max_client > max_main)
        {
            max_main = max_client;
        }
        QueryPerformanceCounter(&end);
        cout << "Maximum value size=";
        cout << max_main << std::endl;
        cout << "Time spent=";
        cout << (end.QuadPart - start.QuadPart) << std::endl;

        // Perform an array comparison
        // Compare your data first
        initial_rand(rawFloatData_main);// The data is about to be sorted, so initialize rawFloatData as a random number array
        QueryPerformanceCounter(&start);
        // float* tempArr = new float[DATANUM];
        // sortSpeedUpOMP(rawFloatData_main, DATANUM, rawFloatData_me_final, tempArr);
        sort_omp(rawFloatData_main, DATANUM, rawFloatData_me_final);
        int rightOrNot1 = check(rawFloatData_me_final, DATANUM);
        cout << "Is server's sorting successful: ";
        if (rightOrNot1 == 1)
            cout << "Sorted successfully!" << endl;
        // Data reception part
        QueryPerformanceCounter(&start1);
        while (1)
		{
			receivesuccess = recv(newConnection, (char*)&rawFloatData_you[received], receiveonce * sizeof(float), NULL);
			if (receivesuccess == -1)
			{
				std::cout << "error";
			}
			else
			{
				received = received + receivesuccess / sizeof(float);
			}
			//cout << received << endl;//输出收到了多少数据，这一句加上后调试的时候可以更直观看到有没有出现掉包，以及传输的进度
			if (received >= DATANUM - 50)
			{
				break;
			}
		}
        QueryPerformanceCounter(&end1);
        cout << "Communication time: " << (end1.QuadPart - start1.QuadPart) << std::endl;
        Compare(rawFloatData_main, rawFloatData_you, rawFloatData_final, received);
        rightOrNot1 = check(rawFloatData_final, received + DATANUM);
        cout << "Is client's sorting successful: ";
        if (rightOrNot1 == 1)
            cout << "Sorted successfully!" << endl;
        QueryPerformanceCounter(&end);
        cout << "Total time spent: " << (end.QuadPart - start.QuadPart) << "     " << std::endl;
    }
    closesocket(newConnection);
    WSACleanup();
    system("pause");
    return 0;
}