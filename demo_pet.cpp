////
//
////using namespace std;
////queue<FrameTransYale> transQueue;
////queue<FrameTransYale> baseQueue;
////extern ProjectConfig projectConfig;
////
////extern "C" __declspec(dllexport) void start_Acq(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax);
////extern "C" __declspec(dllexport) int Get_data(FrameTransYale& pdata);
////extern "C" __declspec(dllexport) void end_Acq();
////extern "C" __declspec(dllexport) void CapturePhantom();
////extern "C" __declspec(dllexport) int CaptureCalibration(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax, char *dicompath);
////HeadMotionTrack* Xdemodemo;
////
////
////BOOL APIENTRY DllMain(HMODULE hModule, DWORD  ul_reason_for_call, LPVOID lpReserved)
////{
////	switch (ul_reason_for_call)
////	{
////		case DLL_PROCESS_ATTACH:
////		case DLL_THREAD_ATTACH:
////		case DLL_THREAD_DETACH:
////		case DLL_PROCESS_DETACH:
////		break;
////	}
////	return TRUE;
////}
////
////void start_Acq(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax)
////{
////
////
////	Xdemodemo = new HeadMotionTrack();
////	Xdemodemo->initialize(8192, 8192, baseQueue, transQueue);
////
////	projectConfig.algoConfig.xMin = xmin;
////	projectConfig.algoConfig.xMax = xmax;
////	projectConfig.algoConfig.yMin = ymin;
////	projectConfig.algoConfig.yMax = ymax;
////	projectConfig.algoConfig.zMin = zmin;
////	projectConfig.algoConfig.zMax = zmax;
////
////	std::cout << "Enter start_Acq()\n"
////		<< "xmin: " << projectConfig.algoConfig.xMin
////		<< " xmax: " << projectConfig.algoConfig.xMax
////		<< " ymin: " << projectConfig.algoConfig.yMin
////		<< " ymax: " << projectConfig.algoConfig.yMax
////		<< " zmin: " << projectConfig.algoConfig.zMin
////		<< " zmax: " << projectConfig.algoConfig.zMax
////		<< " dicompath: " << projectConfig.directoryConfig.petDicom
////		<< std::endl;
////
////	Xdemodemo->startMotionTrack();
////	Xdemodemo->startTrack_cam();
////	std::cout << "enter into start acq" << endl;
////}
////
////int Get_data(FrameTransYale& pdata)
////{
////	//std::cout << "enter into Get_data acq 2019" << endl;
////	if (!transQueue.empty())
////	{
////		//FrameTrans pdata;
////		pdata = transQueue.front();
////		std::cout << "********" << pdata.transDofs.rx << std::endl;
////		transQueue.pop();
////		return 0;
////	}
////	else
////		return 1;
////	//std::cout << "transQueue is empty" << endl;
////
////}
////
////void end_Acq()
////{
////	std::cout << "enter into end_Acq" << endl;
////	int is_stop = Xdemodemo->endTrack_cam();
////	std::cout << "enter into end_Acq" << endl;
////	delete Xdemodemo;
////	Xdemodemo = NULL;
////}
////
////
////void CapturePhantom()
////{
////	capturePhantom(NULL);
////}
////
////
////int CaptureCalibration(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax, char* dicompath)
////{
////	projectConfig.algoConfig.xMin = xmin;
////	projectConfig.algoConfig.xMax = xmax;
////	projectConfig.algoConfig.yMin = ymin;
////	projectConfig.algoConfig.yMax = ymax;
////	projectConfig.algoConfig.zMin = zmin;
////	projectConfig.algoConfig.zMax = zmax;
////	projectConfig.directoryConfig.petDicom = dicompath;
////
////	std::cout << "Enter CaptureCalibration()\n"
////		<< "xmin: " << projectConfig.algoConfig.xMin
////		<< " xmax: " << projectConfig.algoConfig.xMax
////		<< " ymin: " << projectConfig.algoConfig.yMin
////		<< " ymax: " << projectConfig.algoConfig.yMax
////		<< " zmin: " << projectConfig.algoConfig.zMin
////		<< " zmax: " << projectConfig.algoConfig.zMax
////		<< " dicompath: " << projectConfig.directoryConfig.petDicom
////		<< std::endl;
////
////	Xdemodemo = new HeadMotionTrack();
////	int res = Xdemodemo->alignment();
////	return res;
////}
////
////
////
////
//
//
//#include "motion_correction.h"
//#include "dataFormat.h"
//#include "project_config.h"
//using namespace std;
//
//queue<FrameTransYale> transQueue;
//queue<FrameTransYale> baseQueue;
////extern ProjectConfig projectConfig;
//
//extern "C" __declspec(dllexport) void start_Acq(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax);
//extern "C" __declspec(dllexport) int Get_data(FrameTransYale & pdata);
//extern "C" __declspec(dllexport) int end_Acq();
//extern "C" __declspec(dllexport) int CapturePhantom(char* pointCloudModelOffline);
//extern "C" __declspec(dllexport) int OfflineAlignment(char* dicomFolder, char* pointCloudFile, char* matrixPl2plSavePath, char* matrixPl2dcmSavePath);
//extern "C" __declspec(dllexport) void ResetProjectConfig(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax);
//extern "C" __declspec(dllexport) int loadMatrixFile(char* matrixPl2plPath, char* matrixPl2dcmPath);
//extern "C" __declspec(dllexport) int setLaserCamer(int laserIntensity);
//extern "C" __declspec(dllexport) int CalculateMatrix(string dicomPath);
//extern "C" __declspec(dllexport) int OpenGUIWindow();
//extern "C" __declspec(dllexport) int OpenCamera();
//extern "C" __declspec(dllexport) int StartSavePointCloud(char* pointCloudFolder);
////extern "C" __declspec(dllexport) int SaveMotionCapture();
//extern "C" __declspec(dllexport) int ReadConfigFile(char* configpath);
//extern "C" __declspec(dllexport) int StopMotionCapture();
//extern "C" __declspec(dllexport) int StartMotionCapture(char* sixDofResultSavePath, char* yaleDataResultSavePath);
//extern "C" __declspec(dllexport) int Initialize();
//extern "C" __declspec(dllexport) int CloseCamera();
//extern "C" __declspec(dllexport) int StopSavePointCloud();
//extern "C" __declspec(dllexport) void CorrectTimeDelay(int delayTime);
//extern "C" __declspec(dllexport) int SaveConfigFile(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax,
//	int delayTime,
//	int laserPower,
//	char* Matrix_pl2pl, char* Matrix_pl2dcm,
//	char* pointCloudModelOffline,
//	char* sixDofResult,
//	char* yaleDataResult,
//	char* offlineFolder);
//
//HeadMotionTrack* Xdemodemo;
//
//
//
//BOOL APIENTRY DllMain(HMODULE hModule, DWORD  ul_reason_for_call, LPVOID lpReserved)
//{
//	switch (ul_reason_for_call)
//	{
//	case DLL_PROCESS_ATTACH:
//	case DLL_THREAD_ATTACH:
//	case DLL_THREAD_DETACH:
//	case DLL_PROCESS_DETACH:
//		break;
//	}
//	return TRUE;
//}
//
//
//int Initialize()
//{
//	Xdemodemo = new HeadMotionTrack();
//	Xdemodemo->initialize();
//	//Xdemodemo->startSaveOutputMotion();
//	Xdemodemo->changeGuiShow();
//	return 0;
//}
//
//void start_Acq(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax)
//{
//	//int is_init = Xdemodemo->initialize();
//	//// 初始化
//	////Xdemodemo = new HeadMotionTrack();
//	//Xdemodemo->setTransQueue(&baseQueue, &transQueue);
//	//// 參數設置
//	//Xdemodemo->projectConfig.algoConfig.xMin = xmin;
//	//Xdemodemo->projectConfig.algoConfig.xMax = xmax;
//	//Xdemodemo->projectConfig.algoConfig.yMin = ymin;
//	//Xdemodemo->projectConfig.algoConfig.yMax = ymax;
//	//Xdemodemo->projectConfig.algoConfig.zMin = zmin;
//	//Xdemodemo->projectConfig.algoConfig.zMax = zmax;
//	//std::cout << "Enter start_Acq()\n"
//	//	<< "xmin: " << Xdemodemo->projectConfig.algoConfig.xMin
//	//	<< " xmax: " << Xdemodemo->projectConfig.algoConfig.xMax
//	//	<< " ymin: " << Xdemodemo->projectConfig.algoConfig.yMin
//	//	<< " ymax: " << Xdemodemo->projectConfig.algoConfig.yMax
//	//	<< " zmin: " << Xdemodemo->projectConfig.algoConfig.zMin
//	//	<< " zmax: " << Xdemodemo->projectConfig.algoConfig.zMax
//	//	<< " dicompath: " << Xdemodemo->projectConfig.directoryConfig.petDicom
//	//	<< std::endl;
//
//	// 開始采集
//	Xdemodemo->startMotionTrack();
//	//Sleep(1000000);
//	std::cout << "enter into start acq" << endl;
//}
//
//int Get_data(FrameTransYale& pdata)
//{
//	std::cout << "enter into Get_data" << endl;
//	if (!transQueue.empty())
//	{
//		//FrameTrans pdata;
//		pdata = transQueue.front();
//		std::cout << "********" << pdata.transDofs.rx << std::endl;
//		transQueue.pop();
//		return 0;
//	}
//
//	else
//	{
//		std::cout << "transQueue is empty" << std::endl;
//		return 1;
//	}
//
//}
//
//int end_Acq()
//{
//	std::cout << "enter into end_Acq" << endl;
//	return Xdemodemo->stopMotionTrack();
//	///*delete Xdemodemo;
//	//Xdem*/odemo = NULL;
//}
//
//int CapturePhantom(char* pointCloudModelOffline)
//{
//	Xdemodemo->projectConfig.directoryConfig.pointCloudModelOffline = pointCloudModelOffline;
//	std::cout << "pointCloudModelOffline = " << Xdemodemo->projectConfig.directoryConfig.pointCloudModelOffline << std::endl;
//	return Xdemodemo->capturePhantom();
//}
//
//int OfflineAlignment(char* dicomFolder, char* pointCloudFile, char* matrixPl2plSavePath, char* matrixPl2dcmSavePath)
//{
//	Xdemodemo->projectConfig.directoryConfig.petDicom = dicomFolder;
//	Xdemodemo->projectConfig.directoryConfig.pointCloudModelOffline = pointCloudFile;
//	Xdemodemo->projectConfig.directoryConfig.Matrix_pl2pl = matrixPl2plSavePath;
//	Xdemodemo->projectConfig.directoryConfig.Matrix_pl2dcm = matrixPl2dcmSavePath;
//
//	std::cout << "Enter CaptureCalibration()\n"
//		<< "xmin: " << Xdemodemo->projectConfig.algoConfig.xMin
//		<< " xmax: " << Xdemodemo->projectConfig.algoConfig.xMax
//		<< " ymin: " << Xdemodemo->projectConfig.algoConfig.yMin
//		<< " ymax: " << Xdemodemo->projectConfig.algoConfig.yMax
//		<< " zmin: " << Xdemodemo->projectConfig.algoConfig.zMin
//		<< " zmax: " << Xdemodemo->projectConfig.algoConfig.zMax
//		<< " dicompath: " << Xdemodemo->projectConfig.directoryConfig.petDicom
//		<< " pointCloudFile: " << Xdemodemo->projectConfig.directoryConfig.pointCloudModelOffline
//		<< " Matrix_pl2pl: " << Xdemodemo->projectConfig.directoryConfig.Matrix_pl2pl
//		<< " Matrix_pl2dcm: " << Xdemodemo->projectConfig.directoryConfig.Matrix_pl2dcm
//		<< std::endl;
//
//	int res = Xdemodemo->offlineAlignment();
//	return res;
//}
//
//void ResetProjectConfig(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax)
//{
//	Xdemodemo->projectConfig.algoConfig.xMin = xmin;
//	Xdemodemo->projectConfig.algoConfig.xMax = xmax;
//	Xdemodemo->projectConfig.algoConfig.yMin = ymin;
//	Xdemodemo->projectConfig.algoConfig.yMax = ymax;
//	Xdemodemo->projectConfig.algoConfig.zMin = zmin;
//	Xdemodemo->projectConfig.algoConfig.zMax = zmax;
//
//	std::cout
//		<< "xmin: " << Xdemodemo->projectConfig.algoConfig.xMin
//		<< " xmax: " << Xdemodemo->projectConfig.algoConfig.xMax
//		<< " ymin: " << Xdemodemo->projectConfig.algoConfig.yMin
//		<< " ymax: " << Xdemodemo->projectConfig.algoConfig.yMax
//		<< " zmin: " << Xdemodemo->projectConfig.algoConfig.zMin
//		<< " zmax: " << Xdemodemo->projectConfig.algoConfig.zMax
//		<< " dicompath: " << Xdemodemo->projectConfig.directoryConfig.petDicom
//		<< std::endl;
//}
//
//int loadMatrixFile(char* matrixPl2plPath, char* matrixPl2dcmPath)
//{
//	std::cout << matrixPl2plPath << std::endl;
//	std::cout << matrixPl2dcmPath << std::endl;
//	Xdemodemo->projectConfig.directoryConfig.Matrix_pl2pl = matrixPl2plPath;
//	Xdemodemo->projectConfig.directoryConfig.Matrix_pl2dcm = matrixPl2dcmPath;
//
//	std::cout << "Matrix_pl2pl: " << Xdemodemo->projectConfig.directoryConfig.Matrix_pl2pl << endl;
//	std::cout << "Matrix_pl2dcm: " << Xdemodemo->projectConfig.directoryConfig.Matrix_pl2dcm << endl;
//	return Xdemodemo->readMatrixFile();
//}
//
//int setLaserCamer(int laserIntensity)
//{
//	std::cout << "laser power: " << laserIntensity << endl;
//	try
//	{
//		Xdemodemo->setCamera(laserIntensity);
//	}
//	catch (...)
//	{
//		return 1;
//	}
//	return 0;
//}
//
//int CalculateMatrix(string dicomPath = "")
//{
//	if (dicomPath != "")
//	{
//		Xdemodemo->projectConfig.directoryConfig.petDicom = dicomPath;
//	}
//	int align_result = Xdemodemo->offlineAlignment();
//	return align_result;
//}
//
//int OpenGUIWindow()
//{
//	int is_open = Xdemodemo->openGuiWindow();
//	return is_open;
//}
//
//int OpenCamera()
//{
//	int is_open = Xdemodemo->openCamera();
//	return is_open;
//}
//
//int StartSavePointCloud(char* pointCloudFolder)
//{
//	if (pointCloudFolder != NULL)
//	{
//		Xdemodemo->projectConfig.directoryConfig.offlineFolder = pointCloudFolder;
//	}
//	std::cout << "PointCloud Folder: " << pointCloudFolder << std::endl;
//	return Xdemodemo->startSavePointCloud();
//}
//
//
//
//int ReadConfigFile(char* configpath)
//{
//	string configFilePath = configpath;
//	std::cout << "configFilePath = " << configFilePath << std::endl;
//	//return 0;
//	return Xdemodemo->readConfigFile(configpath);
//}
//
//
//int StopMotionCapture()
//{
//	return Xdemodemo->stopMotionTrack();
//}
//
//
//int StartMotionCapture(char* sixDofResultSavePath, char* yaleDataResultSavePath)
//{
//	std::cout << "enter StartMotionCapture" << std::endl;
//
//	if (sixDofResultSavePath != NULL)
//	{
//		Xdemodemo->projectConfig.directoryConfig.sixDofResult = sixDofResultSavePath;
//	}
//
//	if (yaleDataResultSavePath != NULL)
//	{
//		Xdemodemo->projectConfig.directoryConfig.yaleDataResult = yaleDataResultSavePath;
//	}
//	std::cout << "sixDofResult = " << Xdemodemo->projectConfig.directoryConfig.sixDofResult << std::endl;
//	std::cout << "yaleDataResult = " << Xdemodemo->projectConfig.directoryConfig.yaleDataResult << std::endl;
//
//	return Xdemodemo->startMotionTrack();
//}
//
//int CloseCamera()
//{
//	return Xdemodemo->closeCamera();
//}
//
//int StopSavePointCloud()
//{
//	return Xdemodemo->stopSavePointCloud();
//}
//
//void CorrectTimeDelay(int delayTime)
//{
//	Xdemodemo->projectConfig.camConfig.timeDelay = delayTime;
//}
//
//int SaveConfigFile(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax,
//	int delayTime,
//	int laserPower,
//	char* Matrix_pl2pl, char* Matrix_pl2dcm,
//	char* pointCloudModelOffline,
//	char* sixDofResult,
//	char* yaleDataResult,
//	char* offlineFolder)
//{
//	Xdemodemo->projectConfig.algoConfig.xMin = xmin;
//	Xdemodemo->projectConfig.algoConfig.xMax = xmax;
//	Xdemodemo->projectConfig.algoConfig.yMin = ymin;
//	Xdemodemo->projectConfig.algoConfig.yMax = ymax;
//	Xdemodemo->projectConfig.algoConfig.zMin = zmin;
//	Xdemodemo->projectConfig.algoConfig.zMax = zmax;
//	Xdemodemo->projectConfig.camConfig.timeDelay = delayTime;
//	Xdemodemo->projectConfig.camConfig.laserPower = laserPower;
//	Xdemodemo->projectConfig.directoryConfig.Matrix_pl2pl = Matrix_pl2pl;
//	Xdemodemo->projectConfig.directoryConfig.Matrix_pl2dcm = Matrix_pl2dcm;
//	Xdemodemo->projectConfig.directoryConfig.pointCloudModelOffline = pointCloudModelOffline;
//	Xdemodemo->projectConfig.directoryConfig.sixDofResult = sixDofResult;
//	Xdemodemo->projectConfig.directoryConfig.yaleDataResult = yaleDataResult;
//	Xdemodemo->projectConfig.directoryConfig.offlineFolder = offlineFolder;
//
//	std::cout
//		<< "xmin: " << Xdemodemo->projectConfig.algoConfig.xMin
//		<< " xmax: " << Xdemodemo->projectConfig.algoConfig.xMax
//		<< " ymin: " << Xdemodemo->projectConfig.algoConfig.yMin
//		<< " ymax: " << Xdemodemo->projectConfig.algoConfig.yMax
//		<< " zmin: " << Xdemodemo->projectConfig.algoConfig.zMin
//		<< " zmax: " << Xdemodemo->projectConfig.algoConfig.zMax
//		<< " delaytime: " << Xdemodemo->projectConfig.camConfig.timeDelay
//		<< " laserpower: " << Xdemodemo->projectConfig.camConfig.laserPower
//		<< " Matrix_pl2pl: " << Xdemodemo->projectConfig.directoryConfig.Matrix_pl2pl
//		<< " Matrix_pl2dcm: " << Xdemodemo->projectConfig.directoryConfig.Matrix_pl2dcm
//		<< " pointCloudModelOffline: " << Xdemodemo->projectConfig.directoryConfig.pointCloudModelOffline
//		<< " sixDofResult: " << Xdemodemo->projectConfig.directoryConfig.sixDofResult
//		<< " yaleDataResult: " << Xdemodemo->projectConfig.directoryConfig.yaleDataResult
//		<< " offlineFolder: " << Xdemodemo->projectConfig.directoryConfig.offlineFolder
//		<< std::endl;
//
//	return Xdemodemo->saveConfigFile();
//}