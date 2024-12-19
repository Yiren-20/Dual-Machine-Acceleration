#include "speedup_linux.h"
#define SUBDATANUM 1000000//125000  //2000000
#define CPU_THREADS 64
#define DATANUM (SUBDATANUM * CPU_THREADS) 
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
    struct sockaddr_in addr; // Address assigned to the Connection socket
    socklen_t sizeofaddr = sizeof(addr);
    addr.sin_family = AF_INET;
    addr.sin_port = htons(MyPort);                        //端口号
    addr.sin_addr.s_addr = inet_addr(MyAdd);
    int Connection = socket(AF_INET, SOCK_STREAM, NULL);
    unsigned char  service_type = 0xe0;
    if(setsockopt(Connection, SOL_IP/*IPPROTO_IP*/, IP_TOS, (void *)&service_type, sizeof(service_type)) < 0) {
        printf("setsockopt(IP_TOS) failed:\n");
    }
    if (connect(Connection, (struct sockaddr*)&addr, sizeofaddr) != 0) // Connection
    {
        cout<<"Bad Connection."<<endl;
        return 0;
    }
    float* rawFloatData_main = new float[DATANUM];
    float* rawFloatDataResult = new float[DATANUM];
    initial(rawFloatData_main);
    float summary_main = 0;
#pragma omp parallel for
    for (long i = 0; i < MAX_THREADS; i++)
        summary_main = summary_main + sumSpeedUpOMP(rawFloatData_main + i * SUBDATANUM, SUBDATANUM);
    send(Connection, (char*)&summary_main, sizeof(float), NULL);
    std::cout << "Sum result, summary_client = " << summary_main << std::endl;
    std::cout << "Sent. Waiting......" << std::endl;
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

    initial_rand(rawFloatData_main);
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
    
    deleteraw(rawFloatData_main);
    deleteraw(rawFloatDataResult);
    close(Connection);
    cout<<"well done!"<<endl;
    return 0;

}