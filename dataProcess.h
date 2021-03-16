#pragma once

#include "point_cloud.hpp"
#include "pointCloud_io.hpp"
#include "matrix.h"
#include "dataFormat.h"
#include "project_config.h"

//int postProcess_toCPU(float*& d_outputX, float*& d_outputY, float*& d_outputZ, unsigned int size, PointCloud<PointXYZ>& pointCloud);
//int preProcess_toGPU(float*& d_outputX, float*& d_outputY, float*& d_outputZ, float* h_input, unsigned int size);
int preProcess_rotate(PointCloud<PointXYZ>& pointCloud);
int preProcess_rectangle(PointCloud<PointXYZ>& pointCloud, const ProjectConfig& projectConfig);
int preProcess_rectangle(PointCloud<PointXYZ>& pointCloud);
void matrix2angle(Matrix& result_trans, TransDof& result_angle);
void angle2matrix(Matrix& result_trans, TransDof& result_angle);
int preProcess_xyRotate(PointCloud<PointXYZ>& pointCloud);
int plTransform(PointCloud<PointXYZ>& pointCloud, Matrix& result_trans);
int plTransform(vector<PointCloud<PointXYZ>>& pointCloudRef, Matrix& result_trans);

