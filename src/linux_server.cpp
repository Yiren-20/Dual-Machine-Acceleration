#include "speedup_linux.h"
#define SUBDATANUM 1000000//125000  //2000000
#define CPU_THREADS 64
#define DATANUM (SUBDATANUM * CPU_THREADS) 

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

	cout << "---------------server here----------------" << endl;
	cout << "Waiting for connection......" << endl;
	struct sockaddr_in addr;
	socklen_t addrlen = sizeof(addr);
	addr.sin_family = AF_INET;
    addr.sin_port = htons(MyPort);                        //端口号
    addr.sin_addr.s_addr = inet_addr(MyAdd);   // 服务端的IP地址  192.168.253.129/127.0.0.1
	int sListen = socket(AF_INET, SOCK_STREAM,0);

	unsigned char  service_type = 0xe0;
	if(setsockopt(sListen, SOL_IP/*IPPROTO_IP*/, IP_TOS, (void *)&service_type, sizeof(service_type)) < 0) {
		printf("setsockopt(IP_TOS) failed:\n");
    }

	bind(sListen, (struct sockaddr*)&addr, sizeof(addr));
	listen(sListen, SOMAXCONN);
	float* rawFloatData_main = new float[DATANUM];
    float* rawFloatData_me_final = new float[DATANUM];
    float* rawFloatData_you = new float[DATANUM];
    float* rawFloatData_final = new float[2 * DATANUM];
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
	int newConnection;
	timeval start, end;
	newConnection = accept(sListen, (struct sockaddr*)&addr, &addrlen);
	if(newConnection == 0)
	{
		std::cout << "Bad connection." << std::endl;
	}
	else{
		std::cout << "Good connection." << std::endl;
		gettimeofday(&start, NULL);
		summary_main = 0;
#pragma omp parallel for
        for (long i = 0; i < MAX_THREADS; i++)
            summary_main = summary_main + sumSpeedUpOMP(rawFloatData_main + i * SUBDATANUM, SUBDATANUM);
        //Sleep(1);
        recv(newConnection, (char*)&summary_client, sizeof(float), NULL);
        summary_main = summary_main + summary_client;
		gettimeofday(&end, NULL);
		cout << "Sum size=";
        cout << summary_main << std::endl;
        cout << "Time spent=";
        cout << 1e3 * (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1e3 << "ms" <<  std::endl;
		gettimeofday(&start, NULL);
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
		gettimeofday(&end, NULL);
		cout << "Maximum value size=";
        cout << max_main << std::endl;
        cout << "Time spent=";
        cout << 1e3 * (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1e3 << "ms" <<  std::endl;
	
		initial_rand(rawFloatData_main);
		gettimeofday(&start, NULL);
		        sort_omp(rawFloatData_main, DATANUM, rawFloatData_me_final);
        int rightOrNot1 = check(rawFloatData_me_final, DATANUM);
        cout << "Is server's sorting successful: ";
        if (rightOrNot1 == 1)
            cout << "Sorted successfully!" << endl;
        // Data reception part
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
        Compare(rawFloatData_main, rawFloatData_you, rawFloatData_final, received);
        rightOrNot1 = check(rawFloatData_final, received + DATANUM);
        cout << "Is client's sorting successful: ";
        if (rightOrNot1 == 1)
            cout << "Sorted successfully!" << endl;
        gettimeofday(&end, NULL);
        cout << 1e3 * (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1e3 << "ms" <<  std::endl;
	
	}
	close(newConnection);
	cout<<"well done!"<<endl;
	return 0;
}