#pragma once
#include<chrono>

using namespace std;

struct TransDof {
	float rx;
	float ry;
	float rz;
	float x;
	float y;
	float z;
};

struct OutTransYale {
	float ma00, ma01, ma02, ma03;
	float ma10, ma11, ma12, ma13;
	float ma20, ma21, ma22, ma23;
};

//extern "C" void fun(void);




struct FrameTransYale {
	unsigned int index;
	TransDof transDofs;
	OutTransYale transMat;
	//chrono::system_clock::time_point timePoint;
	__int64 timePoint;
};




struct FrameTransM {
	OutTransYale trans;
	chrono::system_clock::time_point timePoint;
};


