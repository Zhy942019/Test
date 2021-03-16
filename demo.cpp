#include "motion_correction.h"
#include <queue>
#include <chrono>
#include <functional>
#include <thread>
#include "dataFormat.h"
#include "project_config.h"
#include <iostream>


using namespace std;


queue<FrameTransYale> transQueue;
queue<FrameTransYale> baseQueue;
queue<string> statusQueue;
extern ProjectConfig projectConfig;


int main()
{


	{


// 新建对象：
		HeadMotionTrack demo; 
		demo.initialize(); // 生成movetrack.txt
		demo.setTransQueue(&baseQueue, &transQueue);//有用吗？？？并发???
		//demo.setStatusQueue(&statusQueue);
		//demo.saveConfigFile(); 
		demo.readConfigFile();
		// std::thread thr1(std::mem_fn(&HeadMotionTrack::openGui), &demo);





		//// 保存数据-跑起来
		//// Zhaohui：最开始
		//demo.startSavePointCloud();
		//demo.changeGuiShow();
		//demo.startMotionTrack();
		//demo.openCamera();
		//cout << "test1" << endl;
		//cout << "test2" << endl;




		//// 离线track：trackFile()-跑起来
		//// Zhaohui：
		//// 修改数据: zmin = 0+1=1，zmax =58+1=59
		//demo.trackFile(demo.projectConfig.directoryConfig.offlineFolder, baseQueue, transQueue); //有用吗？？？并发???
		//cout << "test1" << endl;
		//cout << "test2" << endl;




		//// 在线track：motionTrack()-跑起来
		//// Zhaohui：
		//// demo.startSavePointCloud();
		//demo.changeGuiShow();
		//demo.openGuiWindow();
		//demo.startMotionTrack(); // motionTrackFlag = 1; needRefFlag = 1; needTrackBaseFlag = 1;
		//demo.openCamera(); // 相机开，while，相机关
		////cout << "test1" << endl;
		////cout << "test2" << endl;




		//// 保护readMatrixFile()
		//// Zhaohui 2021.03.04
		//demo.readMatrixFile();
		//cout << "test1" << endl;
		//cout << "test2" << endl;




		//// 保护readConfigFile()
		//// Zhaohui 2021.03.05
		//demo.readConfigFile();
		//cout << "test1" << endl;
		//cout << "test2" << endl;





		//// Alignment：离线校准-跑起来
		//// Zhaohui 2021.03.05-2021.03.09
		//// 外部配置文件 必须手动修改：camPose = 1， MR因为小所以不进行mean操作
		//// head_motion_track.h 必须手动修改：alignOnlineFlag = 0 或者 offlineAlignment()也可以alignOnlineFlag = 0
		//// 2.2：mr线
		//demo.projectConfig.directoryConfig.petDicom = "D:\\vs_projekt\\0301motion_correction\\calibration\\202012223dcm\\t1_quick3d_tra_fs_601"; 
		//// 1.1：表面cad点云：Template 
		//demo.projectConfig.directoryConfig.pointCloudTemplate = "D:\\vs_projekt\\0301motion_correction\\calibration\\alignment_model2.pcd"; 
		//// 1.2: 表面cam点云：cam  // xyz+color：nx6 或 xyz： nx3 都行
		//demo.projectConfig.directoryConfig.pointCloudModelOffline = "D:\\vs_projekt\\0301motion_correction\\calibration\\20201223_901_pcd\\points-7-seg.xyz"; 
		//demo.offlineAlignment();
		//demo.readMatrixFile();
		//cout << "test1" << endl;
		//cout << "test2" << endl;





		//// Alignment：在线校准-跑起来：无相机
		//// Zhaohui 2021.03.09-2021.03.10
		//// head_motion_track.h 必须手动修改： alignOnlineFlag = 1 !!!!!!!!!!!!!!!!
		//// 没有变
		//demo.projectConfig.directoryConfig.petDicom = "D:\\vs_projekt\\0301motion_correction\\calibration\\202012223dcm\\t1_quick3d_tra_fs_601"; 
		//// 没有变
		//demo.projectConfig.directoryConfig.pointCloudTemplate = "D:\\vs_projekt\\0301motion_correction\\calibration\\alignment_model2.pcd"; 
		//// 可实时存一个: 所以可以像家煦一样把这项注掉
		//// pointCloudModelOnline -> Q
		//demo.startAlignment(); // needAlignFlag = 1; needRefFlag = 1; // 为什么不加上：alignOnlineFlag = 1
		//demo.changeGuiShow();
		//demo.openGuiWindow();
		//demo.openCamera(); // 需要相机  // alignOnlineFlag = 1, 直接使用Q
		//cout << "test1" << endl;
		//cout << "test2" << endl;





		//// Alignment：离线校准 面 + pclAllInOne()
		//// Zhaohui 2021.03.10-2021.03.10
		//// 1：10个原始分割
		//// 2：10个分割点云生成一个点云
		//// 3：融合到系统 2021.03.12
		//// 4: 把code分割过的，人工去除噪点，再运行看结果 2021.03.15
		//// 5: 或者人直接分割原始10个点云再运行 2021.03.15
		//// demo.offlineAlignment() // alignment() // pclToPclAlign()：增加顺序功能pclAllInOne()
		//demo.projectConfig.directoryConfig.pointCloudTemplate = "D:\\vs_projekt\\0301motion_correction\\calibration\\alignment_model2.pcd"; 
		//demo.projectConfig.directoryConfig.pointCloudModelOffline = "D:\\vs_projekt\\0301motion_correction\\calibration\\20201223_901_pcd"; 
		//demo.offlineAlignment();
		//// 把分割过的，人工去除噪点，再运行看结果
		//cout << "test1" << endl;
		//cout << "test2" << endl;

		



		// Alignment：离线校准 线 + dicom正确排序
		// Zhaohui 2021.03.11-2021.03.11
		// 1: dicom正确排序
		// 2: 把dicom融合 2021.03.15
		// demo.offlineAlignment() // alignment() // pclToDcmAlign()中：dicom正确排序
		demo.projectConfig.directoryConfig.petDicom= "D:\\vs_projekt\\0301motion_correction\\VTC"; // HC VTC
		demo.offlineAlignment();

		cout << "test1" << endl;

		cout << "test2" << endl;





		//// xml: allInOneNum
		//// Zhaohui 2021.03.12-2021.03.12
		//// 往xml加入新目录
		//// 从xml新目录读取
		//demo.saveConfigFile(); 
		//demo.readConfigFile();
		//cout << demo.projectConfig.algoConfig.allInOneNum << endl;
		//cout << demo.projectConfig.algoConfig.gridSize2 << endl;
		//cout << "test1" << endl;
		//cout << "test2" << endl;




















/////////////////////////////////////////////////////////////////////////////////////原始:全部注销


////ThreadMes mes;
////mes.pBaseQueue = &baseQueue;
////mes.pTransQueue = &transQueue;
////CreateThread(NULL, 0, Threadproc, &mes, 0, NULL);
//
//
////// 新建对象：
////		HeadMotionTrack demo;
////		demo.initialize(); // 生成movetrack.txt
////		demo.setTransQueue(&baseQueue, &transQueue);//有用吗？？？
////		//demo.setStatusQueue(&statusQueue);
//
//
//
//// config.xml
//		// 生成
//		//demo.saveConfigFile(); // 里面定义了createProjectConfig(configFilePath.c_str(), projectConfig);，输入了projectConfig数据
//		//
//		// read
//		demo.readConfigFile();
//		
//// Alignment
//		//std::thread thr1(std::mem_fn(&HeadMotionTrack::openGui), &demo);
//		//demo.changeStatusLog();
//		//demo.openGuiWindow();
//		//demo.startAlignment();
//		//demo.projectConfig.directoryConfig.pointCloudModelOffline = "../temp/captureForAlignmentOffline.xyz";
//		//Sleep(5000);
//		//demo.setCamAxis();
//		//demo.startSavePointCloud();
//
//
//
//// 显示窗口：分割后的实时点云
//		demo.changeGuiShow(); 
//
//
//		// 显示：分割后-参考帧？？？
//		//demo.changeGuiShow();
//		//
//		//demo.projectConfig.algoConfig.xMax = 100;
//		//demo.projectConfig.algoConfig.yMax = 100;
//		//demo.projectConfig.algoConfig.zMax = 600;
//		//demo.projectConfig.algoConfig.xMin = -100;
//		//demo.projectConfig.algoConfig.yMin = -100;
//		//demo.projectConfig.algoConfig.zMin = 150;
//		//
//		//demo.projectConfig.camConfig.laserPower = 66;
//		//demo.projectConfig.algoConfig.gridSize2 = 3.0;
//		
//
//// 打开窗口
//		demo.openGuiWindow(); 
//
//
//// 保存离线数据: 分割后的点云
//		//demo.startSavePointCloud();
//
//
//// 在线
//		demo.startMotionTrack(); // motionTrackFlag = 1; needRefFlag = 1; needTrackBaseFlag = 1;
//
//
//// 在线 
//		demo.openCamera(); // 相机开，while，相机关
//
//
//
//
//// 离线
//		// 读取数据数量: zmin，zmax
//		//demo.readConfigFile();
//		// 跑离线，跑测试
//		// demo.trackFile(demo.projectConfig.directoryConfig.offlineFolder, baseQueue, transQueue);
//
//		
//
//
//// 离线校准
//		//demo.projectConfig.directoryConfig.pointCloudTemplate = 
//		//demo.projectConfig.directoryConfig.petDicom = "E:\\data\\CT11";
//		//demo.projectConfig.directoryConfig.pointCloudModelOffline = "E:\\data\\problemData\\pet2pcd\\captureForAlignmentOffline.xyz"; // xyz+color：nx6 // 不行：只需要xyz： nx3
//		//demo.offlineAlignment();




		Sleep(1000000000000);
		//Sleep(100000);
		//demo.offlineAlignment();
		//projectConfig.algoConfig.xMax = 123;

		//demo.startTrack_file("E:\\data\\20201205data\\movBedYZ", baseQueue, transQueue); // no need
		//cout << "waiting" << endl;
		//Sleep(100000);
		//cout << "waiting end" << endl;
		//demo.endTrack_cam();


	}


	system("pause");
	return 0;


}



