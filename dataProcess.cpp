#include "point_cloud.hpp"
#include "pointCloud_io.hpp"
#include "matrix.h"
#include "dataFormat.h"
#include "dataProcess.h"
#include "matrix.h"

#include <Eigen/Dense>

#define M_PI 3.1415926535897932384626433832795

int preProcess_rotate(PointCloud<PointXYZ>& pointCloud) {

	PointCloud<PointXYZ> pointCloudT;
	for (int i = 0; i < pointCloud.size(); i = i + 4) {
		pointCloudT.points.push_back(pointCloud.points[i]);
	}

	float RR[3][3];
	RR[0][0] = 1.0000000000000000; RR[0][1] = 0.00000000000000000; RR[0][2] = 0.0000000000000000;
	RR[1][0] = 0.0000000000000000; RR[1][1] = -0.45399049973954675; RR[1][2] = 0.89100652418836790;
	RR[2][0] = 0.0000000000000000; RR[2][1] = -0.89100652418836790; RR[2][2] = -0.45399049973954675;

	for (int i = 0; i < pointCloudT.size(); i++) {
		PointXYZ temp;
		temp.x = pointCloudT.points[i].x * RR[0][0] + pointCloudT.points[i].y * RR[0][1] + pointCloudT.points[i].z * RR[0][2];
		temp.y = pointCloudT.points[i].x * RR[1][0] + pointCloudT.points[i].y * RR[1][1] + pointCloudT.points[i].z * RR[1][2];
		temp.z = pointCloudT.points[i].x * RR[2][0] + pointCloudT.points[i].y * RR[2][1] + pointCloudT.points[i].z * RR[2][2];
		pointCloudT.points[i] = temp;
	}
	float minZ = pointCloudT.points[0].z;
	for (int i = 0; i < pointCloudT.size(); i++) {
		if (pointCloudT.points[i].z < 50) continue;
		if (minZ > pointCloudT.points[i].z) {
			minZ = pointCloudT.points[i].z;
		}
	}
	pointCloud.clear();
	for (int i = 0; i < pointCloudT.size(); i++) {
		if (pointCloudT.points[i].z < minZ + 180 && pointCloudT.points[i].z>50 && pointCloudT.points[i].x > -150 && pointCloudT.points[i].x < 150) {
			PointXYZ temp;
			temp.x = pointCloudT.points[i].x;
			temp.y = pointCloudT.points[i].y;
			temp.z = pointCloudT.points[i].z;
			pointCloud.points.push_back(temp);
		}
	}

	return 0;
}

int preProcess_rectangle(PointCloud<PointXYZ>& pointCloud) {

	PointCloud<PointXYZ> pointCloudT;
	for (int i = 0; i < pointCloud.size(); i = i + 1) {
		pointCloudT.points.push_back(pointCloud.points[i]);
	}
	pointCloud.clear();
	for (int i = 0; i < pointCloudT.size(); i++) {
		PointXYZ temp = pointCloudT.points[i];
		if (temp.x<80 && temp.x>-80 && temp.y > 50 && temp.y < 180 && temp.z < 550 && temp.z > 220) {
			//temp.x = -temp.x;
			//temp.y = -temp.y;
			pointCloud.points.push_back(temp);
		}
	}
	return 0;
}

int preProcess_rectangle(PointCloud<PointXYZ>& pointCloud, const ProjectConfig& projectConfig) {

	PointCloud<PointXYZ> pointCloudT;
	for (int i = 0; i < pointCloud.size(); i = i + 1) {
		pointCloudT.points.push_back(pointCloud.points[i]);
	}
	pointCloud.clear();
	for (int i = 0; i < pointCloudT.size(); i++) {
		PointXYZ temp = pointCloudT.points[i];
		if (temp.x<projectConfig.algoConfig.xMax && temp.x>projectConfig.algoConfig.xMin && temp.y > projectConfig.algoConfig.yMin && 
			temp.y < projectConfig.algoConfig.yMax && temp.z < projectConfig.algoConfig.zMax && temp.z > projectConfig.algoConfig.zMin) {
			//temp.x = -temp.x;
			//temp.y = -temp.y;
			pointCloud.points.push_back(temp);
		}
	}
	return 0;
}



// 
void matrix2angle(Matrix &result_trans, TransDof &result_angle)
{
	float ax, ay, az;
	if (result_trans.val[2][0] == 1 || result_trans.val[2][0] == -1)
	{
		az = 0;
		float dlta;
		dlta = atan2(result_trans.val[0][1], result_trans.val[0][2]);
		if (result_trans.val[2][0] == -1)
		{
			ay = M_PI / 2;
			ax = az + dlta;
		}
		else
		{
			ay = -M_PI / 2;
			ax = -az + dlta;
		}
	}
	else
	{
		ay = -asin(result_trans.val[2][0]);
		ax = atan2(result_trans.val[2][1] / cos(ay), result_trans.val[2][2] / cos(ay));
		az = atan2(result_trans.val[1][0] / cos(ay), result_trans.val[0][0] / cos(ay));
	}
	result_angle.rx = ax;
	result_angle.ry = ay;
	result_angle.rz = az;
	result_angle.x = result_trans.val[0][3];
	result_angle.y = result_trans.val[1][3];
	result_angle.z = result_trans.val[2][3];
}



//
void angle2matrix(Matrix& result_trans, TransDof& result_angle)
{
	::Eigen::Vector3d ea0(result_angle.rz, result_angle.ry, result_angle.rx);
	::Eigen::Matrix3d R;
	R = ::Eigen::AngleAxisd(ea0[0], ::Eigen::Vector3d::UnitZ()) * 
		::Eigen::AngleAxisd(ea0[1], ::Eigen::Vector3d::UnitY()) * 
		::Eigen::AngleAxisd(ea0[2], ::Eigen::Vector3d::UnitX());
	Eigen::Matrix4f Ti = Eigen::Matrix4f::Identity();

	result_trans.val[0][0] = R(0, 0); result_trans.val[0][1] = R(0, 1); result_trans.val[0][2] = R(0, 2); result_trans.val[0][3] = result_angle.x;
	result_trans.val[1][0] = R(1, 0); result_trans.val[1][1] = R(1, 1); result_trans.val[1][2] = R(1, 2); result_trans.val[1][3] = result_angle.y;
	result_trans.val[2][0] = R(2, 0); result_trans.val[2][1] = R(2, 1); result_trans.val[2][2] = R(2, 2); result_trans.val[2][3] = result_angle.z;
}




int preProcess_xyRotate(PointCloud<PointXYZ>& pointCloud) {

	for (int i = 0; i < pointCloud.size(); i = i + 1) {
		pointCloud.points[i].x = -pointCloud.points[i].x;
		pointCloud.points[i].y = -pointCloud.points[i].y;
	}

	return 0;
}

int plTransform(PointCloud<PointXYZ>& pointCloud, Matrix& result_trans) {
	for (int i = 0; i < pointCloud.size(); i++) {
		PointXYZ temp = pointCloud.points[i];
		pointCloud.points[i].x = result_trans.val[0][0] * temp.x + result_trans.val[0][1] * temp.y + result_trans.val[0][2] * temp.z + result_trans.val[0][3];
		pointCloud.points[i].y = result_trans.val[1][0] * temp.x + result_trans.val[1][1] * temp.y + result_trans.val[1][2] * temp.z + result_trans.val[1][3];
		pointCloud.points[i].z = result_trans.val[2][0] * temp.x + result_trans.val[2][1] * temp.y + result_trans.val[2][2] * temp.z + result_trans.val[2][3];
	}
	return 0;
}

int plTransform(vector<PointCloud<PointXYZ>>& pointCloudRef, Matrix& result_trans) {
	for (int j = 0; j < pointCloudRef.size(); j++) {
		PointCloud<PointXYZ>& pointCloud = pointCloudRef[j];
		for (int i = 0; i < pointCloud.size(); i++) {
			PointXYZ temp = pointCloud.points[i];
			pointCloud.points[i].x = result_trans.val[0][0] * temp.x + result_trans.val[0][1] * temp.y + result_trans.val[0][2] * temp.z + result_trans.val[0][3];
			pointCloud.points[i].y = result_trans.val[1][0] * temp.x + result_trans.val[1][1] * temp.y + result_trans.val[1][2] * temp.z + result_trans.val[1][3];
			pointCloud.points[i].z = result_trans.val[2][0] * temp.x + result_trans.val[2][1] * temp.y + result_trans.val[2][2] * temp.z + result_trans.val[2][3];
		}
	}
	return 0;
}

