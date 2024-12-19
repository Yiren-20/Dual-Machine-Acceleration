#pragma comment(lib,"ws2_32.lib")

#include "speedup_windows.h"
#define SUBDATANUM 1000000
#define DATANUM (SUBDATANUM * MAX_THREADS)  /*这个数值是总数据量*/
#define SENDONCE 1000

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

int main()
{
    cout << "----------------client here----------------" << endl;
    WSAData wsaData;
    WORD DllVersion = MAKEWORD(2, 1);
    if (WSAStartup(DllVersion, &wsaData) != 0)
    {
        MessageBoxA(NULL, "Winsock startup error", "Error", MB_OK | MB_ICONERROR);
        exit(1);
    }

    SOCKADDR_IN addr; // Address assigned to the Connection socket
    int sizeofaddr = sizeof(addr);
    addr.sin_addr.s_addr = inet_addr(MyAdd); // Address = localhost (replace "MyAdd" with your actual address)
    addr.sin_port = htons(MyPort); // Port = 1111 (replace with your actual port)
    addr.sin_family = AF_INET; // IPv4 Socket

    SOCKET Connection = socket(AF_INET, SOCK_STREAM, NULL);
    if (connect(Connection, (SOCKADDR*)&addr, sizeofaddr) != 0) // Connection
    {
        MessageBoxA(NULL, "Bad Connection", "Error", MB_OK | MB_ICONERROR);
        return 0;
    }

    srand(time(NULL));

    LARGE_INTEGER start;
    LARGE_INTEGER end;

    float* rawFloatData_main = new float[DATANUM];
    float* rawFloatDataResult = new float[DATANUM];

    initial(rawFloatData_main);

    // Calculate sum
    float summary_main = 0;
#pragma omp parallel for
    for (long i = 0; i < MAX_THREADS; i++)
        summary_main = summary_main + sumSpeedUpOMP(rawFloatData_main + i * SUBDATANUM, SUBDATANUM);
    send(Connection, (char*)&summary_main, sizeof(float), NULL);
    std::cout << "Sum result, summary_client = " << summary_main << std::endl;
    std::cout << "Sent. Waiting......" << std::endl;

    // Compare max
    float max_main = 0;
#pragma omp parallel for 
    for (long i = 0; i < MAX_THREADS; i++)
    {
        float tmp = maxSpeedUpOMP(rawFloatData_main + i * SUBDATANUM, SUBDATANUM);
        max_main = (max_main > tmp ? max_main : tmp);
    }
    send(Connection, (char*)&max_main, sizeof(float), NULL);
    std::cout << "Maximum value is: max_client = " << max_main << std::endl;
    std::cout << "Sent. Waiting......" << std::endl;

    // Sorting
    initial_rand(rawFloatData_main);
    // float* tempArr = new float[DATANUM];
    // sortSpeedUpOMP(rawFloatData_main, DATANUM, rawFloatDataResult, tempArr);
    sort_omp(rawFloatData_main, DATANUM, rawFloatDataResult);
    int rightOrNot = 0;
    rightOrNot = check(rawFloatDataResult, DATANUM);
    std::cout << "Is sorting successful: ";
    if (rightOrNot == 1)
        std::cout << "Sorted successfully. Waiting for sending......" << std::endl;
    int sent = 0;
    for (int i = 0; i < DATANUM / SENDONCE; i++)
    {
        sent = send(Connection, (char*)&rawFloatDataResult[i * SENDONCE], SENDONCE * sizeof(float), NULL);
    }
    cout << "Sent. Waiting......" << endl;
    // Cleanup
    deleteraw(rawFloatData_main);
    deleteraw(rawFloatDataResult);
    closesocket(Connection);
    WSACleanup();
    system("pause");
    return 0;
}