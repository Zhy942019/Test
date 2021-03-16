#pragma once
#include "matrix.h"
#include "point_cloud.hpp"

class IcpCu
{

public:


	IcpCu();


	int initialize(int pointFixSize, int pointMoveSize);
	int deInitialize();
	unsigned int pointFixSpace, pointMoveSpace;
	unsigned int pointFixSize, pointMovSize;
	unsigned int combineSize;

	float* d_pointCloudFX, *d_pointCloudFY, *d_pointCloudFZ;
	float* d_pointCloudMX, *d_pointCloudMY, *d_pointCloudMZ;
	float* d_combineCloudX, *d_combineCloudY, *d_combineCloudZ;

	float* d_dist2Matrix;
	float* minDist;
	float* d_pointClNormalX, *d_pointClNormalY, *d_pointClNormalZ;

	unsigned int kNearNum;
	unsigned int *d_index;
	unsigned int *d_nearestIndex;

	int getFix(float* h_input, unsigned int size);
	int getMove(float* h_input, unsigned int size);
	int prepareFix();
	int computeFFDistance();
	int sortDistance();
	int computeNormal();
	int findMoveMin();
	int ComputeAb();
	int icpRegistraion(); // 分析过了：没问题
	int transform(float* &d_outputX, float* &d_outputY, float* &d_outputZ, unsigned int size, float* R, float* b);
	int transform(float* &d_outputX, float* &d_outputY, float* &d_outputZ, unsigned int size,
		float R11, float R12, float R13, float R21, float R22, float R23, float R31, float R32, float R33,
		float b1, float b2, float b3);
	int IcpCu::transform(float*& d_outputX, float*& d_outputY, float*& d_outputZ, unsigned int size,
		Matrix& trans);
	int verifyRb();
	int getTransformMatrix();
	int combineMove();
	int refreshFix();
	int set0(float* &d_outputX, float* &d_outputY, float* &d_outputZ, unsigned int size);
	int ComputeConfidence(float voxelSize, float& confidence);
	int gpuToHost(float*& d_outputX, float*& d_outputY, float*& d_outputZ, unsigned int size, PointCloud<PointXYZ>& pointCloud);
	int ComputeExpression(float voxelSize, float& expression);
	float alpha;
	float beta;
	float* A, *b;
	float* A_, *b_;
	int *d_Ipiv; /* pivoting sequence */
	int *d_info; /* error info */
	int lwork; /* size of workspac e */
	float *d_work; /* device workspace for getrf */
	float* UVW;
	float* R_;
	float distThres;

	float* RT;
	float* bT;

	Matrix transMovToFix;

	~IcpCu();


};

