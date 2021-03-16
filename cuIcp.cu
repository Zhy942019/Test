#include "cuIcp.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cublas_v2.h> 
#include "svd3_cuda.hpp"
#include <cusolverDn.h>
#include "matrix.h"
#include <vector>
#include <iostream>
#include <algorithm>



static cusolverStatus_t status;
static cublasHandle_t handle;
static cusolverDnHandle_t cusolverH;

__global__ void cuComputeDistance3FF(float* pointCRefX, float* pointCRefY, float* pointCRefZ, unsigned int sizeRef,
	unsigned int space, float*QR) {
	unsigned int bidX = blockIdx.x;
	unsigned int bidY = blockIdx.y;
	unsigned int tid = threadIdx.x;

	__shared__ float shared_pointQueryX[32];
	__shared__ float shared_pointQueryY[32];
	__shared__ float shared_pointQueryZ[32];

	float pointRefX, pointRefY, pointRefZ;

	if ((bidX * blockDim.x + tid) < sizeRef) {
		pointRefX = pointCRefX[bidX * blockDim.x + tid];
		pointRefY = pointCRefY[bidX * blockDim.x + tid];
		pointRefZ = pointCRefZ[bidX * blockDim.x + tid];
	}
	if ((bidY * blockDim.x + tid) < sizeRef) {
		shared_pointQueryX[tid] = pointCRefX[bidY * blockDim.x + tid];
		shared_pointQueryY[tid] = pointCRefY[bidY * blockDim.x + tid];
		shared_pointQueryZ[tid] = pointCRefZ[bidY * blockDim.x + tid];
	}
	__syncthreads();


	for (int i = 0;i < blockDim.x;i++) {

		if ((bidY * blockDim.x + i) < sizeRef || (bidX * blockDim.x + tid) < sizeRef) {
			QR[(bidY * blockDim.x + i)*space + bidX * blockDim.x + tid] = (shared_pointQueryX[i] - pointRefX)*(shared_pointQueryX[i] - pointRefX) +
				(shared_pointQueryY[i] - pointRefY)*(shared_pointQueryY[i] - pointRefY) +
				(shared_pointQueryZ[i] - pointRefZ)*(shared_pointQueryZ[i] - pointRefZ);
		}
	}

}

__global__ void cuSortDistance3(unsigned int* kIndex, float*QR, unsigned int kNear, unsigned int sizeRef, unsigned int space, float* maxIst) {
	unsigned int tid = threadIdx.x;
	unsigned int bid = blockIdx.x;
	__shared__ volatile unsigned int nearPoint[16];
	__shared__ volatile float nearDist[16];
	float maxDist;
	__shared__ float tempDist[32];

	if (tid < kNear) {
		nearDist[tid] = 8e8;
	}
	if (tid == 0) {
		maxDist = 8e8;
	}

	for (int i = 0; i < sizeRef; i = i + blockDim.x) {
		if (i + 32 > sizeRef) {
			unsigned int temp = sizeRef - i;
			if (tid >= temp) continue;
			tempDist[tid] = QR[bid * space + i + tid];
			__syncthreads();
			if (tid != 0) continue;
			for (int ii = 0; ii < temp; ii++) {
				if (tempDist[ii] < maxDist) {
					for (int j = 0; j < kNear; j++) {
						if (tempDist[ii] < nearDist[j]) {
							for (int k = kNear - 1; k > j; k--) {
								nearPoint[k] = nearPoint[k - 1];
								nearDist[k] = nearDist[k - 1];
							}
							nearPoint[j] = i + ii;
							nearDist[j] = tempDist[ii];
							maxDist = nearDist[kNear - 1];
							break;
						}
					}
				}
			}

		}
		else {
			tempDist[tid] = QR[bid * space + i + tid];
			__syncthreads();
			if (tid != 0) continue;
			for (int ii = 0; ii < 32; ii++) {
				if (tempDist[ii] < maxDist) {
					for (int j = 0; j < kNear; j++) {
						if (tempDist[ii] < nearDist[j]) {
							for (int k = kNear - 1; k > j; k--) {
								nearPoint[k] = nearPoint[k - 1];
								nearDist[k] = nearDist[k - 1];
							}
							nearPoint[j] = i + ii;
							nearDist[j] = tempDist[ii];
							maxDist = nearDist[kNear - 1];
							break;
						}
					}
				}
			}
		}

	}
	if (tid < kNear) {
		kIndex[bid*kNear + tid] = nearPoint[tid];
	}
	if (tid == 0) {
		maxIst[bid] = maxDist;
	}
}

__global__ void computeNormal2(float* pointCX, float* pointCY, float* pointCZ, unsigned int size, float* pointCNormalX, float* pointCNormalY, float* pointCNormalZ, unsigned int* kIndex, unsigned int kNear) {
	unsigned int xindex = blockIdx.x*blockDim.x + threadIdx.x;
	float3 point[16] = {};
	float3 mu = {};
	for (int i = 0; i < kNear; i++) {
		point[i].x = pointCX[kIndex[xindex*kNear + i]];
		point[i].y = pointCY[kIndex[xindex*kNear + i]];
		point[i].z = pointCZ[kIndex[xindex*kNear + i]];
		mu.x += point[i].x;
		mu.y += point[i].y;
		mu.z += point[i].z;
	}

	mu.x = mu.x / 16.0;
	mu.y = mu.y / 16.0;
	mu.z = mu.z / 16.0;

	float QMatrix[16][3];
	for (int i = 0; i < 16; i++) {
		QMatrix[i][0] = point[i].x - mu.x;
		QMatrix[i][1] = point[i].y - mu.y;
		QMatrix[i][2] = point[i].z - mu.z;
	}

	float QMatrixT[3][16];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 16; j++) {
			QMatrixT[i][j] = QMatrix[j][i];
		}
	}
	float HMatrix[3][3] = {};
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 16; k++)
				HMatrix[i][j] += QMatrixT[i][k] * QMatrix[k][j];

	float U[3][3], S[3], V[3][3];
	svd(HMatrix[0][0], HMatrix[0][1], HMatrix[0][2],
		HMatrix[1][0], HMatrix[1][1], HMatrix[1][2],
		HMatrix[2][0], HMatrix[2][1], HMatrix[2][2],

		U[0][0], U[0][1], U[0][2],
		U[1][0], U[1][1], U[1][2],
		U[2][0], U[2][1], U[2][2],

		S[0], S[1], S[2],

		V[0][0], V[0][1], V[0][2],
		V[1][0], V[1][1], V[1][2],
		V[2][0], V[2][1], V[2][2]
		);

	// normal
	pointCNormalX[xindex] = U[0][2];
	pointCNormalY[xindex] = U[1][2];
	pointCNormalZ[xindex] = U[2][2];
}

__global__ void cuFindQMin5(float* pointCRefX, float* pointCRefY, float* pointCRefZ, unsigned int sizeRef,
	float* pointCQueryX, float* pointCQueryY, float* pointCQueryZ, unsigned int sizeQuery, unsigned int*index, float* minDist) {
	unsigned int bid = blockIdx.x;
	unsigned int tid = threadIdx.x;
	const unsigned int shared_num = 256;
	float pointQX = pointCQueryX[bid];
	float pointQY = pointCQueryY[bid];
	float pointQZ = pointCQueryZ[bid];

	__shared__ volatile float    shared_dist[shared_num + 1];
	__shared__ float shared_index[shared_num + 1];

	unsigned int shared_off = shared_num;
	float shared_PointRX0 = pointCRefX[tid];
	float shared_PointRY0 = pointCRefY[tid];
	float shared_PointRZ0 = pointCRefZ[tid];
	float shared_dist0 = (pointQX - shared_PointRX0) * (pointQX - shared_PointRX0) + (pointQY - shared_PointRY0) * (pointQY - shared_PointRY0) + (pointQZ - shared_PointRZ0) * (pointQZ - shared_PointRZ0);
	float shared_index0 = tid;
	while (tid + shared_off < sizeRef) {
		if (tid + shared_off >= sizeRef) break;
		float tempPointRX = pointCRefX[tid + shared_off];
		float tempPointRY = pointCRefY[tid + shared_off];
		float tempPointRZ = pointCRefZ[tid + shared_off];
		float tempDist = (pointQX - tempPointRX) * (pointQX - tempPointRX) + (pointQY - tempPointRY) * (pointQY - tempPointRY) + (pointQZ - tempPointRZ) * (pointQZ - tempPointRZ);
		float tempIndex = tid + shared_off;
		if (tempDist < shared_dist0) {
			shared_index0 = tempIndex;
			shared_dist0 = tempDist;
		}
		shared_off += shared_num;
	}
	shared_dist[tid] = shared_dist0;
	shared_index[tid] = shared_index0;
	__syncthreads();
	int i = blockDim.x >> 1;
	while (i != 0) {
		if (tid < i) {
			if (shared_dist[tid] > shared_dist[tid + i]) {
				shared_dist[tid] = shared_dist[tid + i];
				shared_index[tid] = shared_index[tid + i];
			}
		}
		i >>= 1;
		__syncthreads();
	}
	__syncthreads();
	if (tid == 0) {
		index[bid] = shared_index[0];
		minDist[bid] = shared_dist[0];
	}

}

__global__ void cuComputeAb3(float*A, float*b, float* pointCRefX, float* pointCRefY, float* pointCRefZ, unsigned int sizeRef, float* pointCQueryX, 
	float* pointCQueryY, float* pointCQueryZ, unsigned int sizeQuery, float* pointClNormalX, float* pointClNormalY, float* pointClNormalZ, unsigned int*index, float* minDist, float distThres) {
	unsigned int xIndex = blockIdx.x*blockDim.x + threadIdx.x;
	float queryPointX = pointCQueryX[xIndex];
	float queryPointY = pointCQueryY[xIndex];
	float queryPointZ = pointCQueryZ[xIndex];
	float refPointX = pointCRefX[index[xIndex]];
	float refPointY = pointCRefY[index[xIndex]];
	float refPointZ = pointCRefZ[index[xIndex]];

	float normalX;
	float normalY;
	float normalZ;

	if (minDist[xIndex] < distThres) {
		normalX = pointClNormalX[index[xIndex]];
		normalY = pointClNormalY[index[xIndex]];
		normalZ = pointClNormalZ[index[xIndex]];
	}
	else {
		normalX = 0;
		normalY = 0;
		normalZ = 0; 
	}

	A[xIndex * 6 + 0] = normalZ*queryPointY - normalY*queryPointZ;
	A[xIndex * 6 + 1] = normalX*queryPointZ - normalZ*queryPointX;
	A[xIndex * 6 + 2] = normalY*queryPointX - normalX*queryPointY;
	A[xIndex * 6 + 3] = normalX;
	A[xIndex * 6 + 4] = normalY;
	A[xIndex * 6 + 5] = normalZ;
	b[xIndex] = normalX*refPointX + normalY*refPointY + normalZ*refPointZ - normalX*queryPointX - normalY*queryPointY - normalZ*queryPointZ;
}

__global__ void svd3_b_(float* input, float* ouputdata)
{
	svd(
		1, -input[2], input[1],
		input[2], 1, -input[0],
		-input[1], input[0], 1,

		ouputdata[0], ouputdata[1], ouputdata[2],
		ouputdata[3], ouputdata[4], ouputdata[5],
		ouputdata[6], ouputdata[7], ouputdata[8],

		ouputdata[9], ouputdata[10], ouputdata[11],

		ouputdata[12], ouputdata[13], ouputdata[14],
		ouputdata[15], ouputdata[16], ouputdata[17],
		ouputdata[18], ouputdata[19], ouputdata[20]
		);
}

__global__ void cuTransform(float* pointClIX, float* pointClIY, float* pointClIZ, unsigned int sizeI, float*rotateMatrix, float*transMatrix, float* pointClOX, float* pointClOY, float* pointClOZ) {
	unsigned int xIndex = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int tid = threadIdx.x;
	__shared__ float rMatrix[9];
	__shared__ float tMatrix[3];
	if (tid < 9) {
		rMatrix[tid] = rotateMatrix[tid];
	}
	if (tid < 3) {
		tMatrix[tid] = transMatrix[tid];
	}
	__syncthreads();
	if (xIndex < sizeI) {
		pointClOX[xIndex] = rMatrix[0] * pointClIX[xIndex] + rMatrix[3] * pointClIY[xIndex] + rMatrix[6] * pointClIZ[xIndex] + tMatrix[0];
		pointClOY[xIndex] = rMatrix[1] * pointClIX[xIndex] + rMatrix[4] * pointClIY[xIndex] + rMatrix[7] * pointClIZ[xIndex] + tMatrix[1];
		pointClOZ[xIndex] = rMatrix[2] * pointClIX[xIndex] + rMatrix[5] * pointClIY[xIndex] + rMatrix[8] * pointClIZ[xIndex] + tMatrix[2];
	}
}

__global__ void cuTransform2(float* pointClIX, float* pointClIY, float* pointClIZ, unsigned int sizeI,
	float R11, float R12, float R13, float R21, float R22, float R23, float R31, float R32, float R33,
	float b1, float b2, float b3, float* pointClOX, float* pointClOY, float* pointClOZ) {
	unsigned int xIndex = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int tid = threadIdx.x;

	if (xIndex < sizeI) {
		float rMatrix[9];
		float tMatrix[3];
		rMatrix[0] = R11; rMatrix[1] = R21; rMatrix[2] = R31;
		rMatrix[3] = R12; rMatrix[4] = R22; rMatrix[5] = R32;
		rMatrix[6] = R13; rMatrix[7] = R23; rMatrix[8] = R33;
		tMatrix[0] = b1; tMatrix[1] = b2; tMatrix[2] = b3;

		pointClOX[xIndex] = rMatrix[0] * pointClIX[xIndex] + rMatrix[3] * pointClIY[xIndex] + rMatrix[6] * pointClIZ[xIndex] + tMatrix[0];
		pointClOY[xIndex] = rMatrix[1] * pointClIX[xIndex] + rMatrix[4] * pointClIY[xIndex] + rMatrix[7] * pointClIZ[xIndex] + tMatrix[1];
		pointClOZ[xIndex] = rMatrix[2] * pointClIX[xIndex] + rMatrix[5] * pointClIY[xIndex] + rMatrix[8] * pointClIZ[xIndex] + tMatrix[2];

	}
}

//int cuIcp::getBaseTrans(float* h_input, unsigned int size) {
//	float meanX = 0;
//	float meanY = 0;
//	float meanZ = 0;
//
//	for (int i = 0; i < size; i++) {
//		meanX += h_input[i * 3 + 0];
//		meanY += h_input[i * 3 + 1];
//		meanZ += h_input[i * 3 + 2];
//	}
//	this->baseTrans[0] = meanX /= size;
//	this->baseTrans[1] = meanY /= size;
//	this->baseTrans[2] = meanZ /= size;
//	return 0;
//}

int IcpCu::getFix(float* h_input, unsigned int size) {

	float* hX, * hY, * hZ;
	hX = new float[size];
	hY = new float[size];
	hZ = new float[size];
	//if (withTrans) {
	//	for (int i = 0; i < size; i++) {
	//		hX[i] = h_input[i * 3 + 0] - baseTrans[0];
	//		hY[i] = h_input[i * 3 + 1] - baseTrans[1];
	//		hZ[i] = h_input[i * 3 + 2] - baseTrans[2];
	//	}
	//}
	//else {
		for (int i = 0; i < size; i++) {
			hX[i] = h_input[i * 3 + 0];
			hY[i] = h_input[i * 3 + 1];
			hZ[i] = h_input[i * 3 + 2];
		}
	//}
	if (pointFixSpace < size) {
		cudaMemcpy(d_pointCloudFX, hX, sizeof(float) * pointFixSpace, cudaMemcpyHostToDevice);
		cudaMemcpy(d_pointCloudFY, hY, sizeof(float) * pointFixSpace, cudaMemcpyHostToDevice);
		cudaMemcpy(d_pointCloudFZ, hZ, sizeof(float) * pointFixSpace, cudaMemcpyHostToDevice);
		this->pointFixSize = pointFixSpace;
	}
	else {
		cudaMemcpy(d_pointCloudFX, hX, sizeof(float) * size, cudaMemcpyHostToDevice);
		cudaMemcpy(d_pointCloudFY, hY, sizeof(float) * size, cudaMemcpyHostToDevice);
		cudaMemcpy(d_pointCloudFZ, hZ, sizeof(float) * size, cudaMemcpyHostToDevice);
		this->pointFixSize = size;
	}

	delete[] hX;
	delete[] hY;
	delete[] hZ;
	
	return 0;
}




//  pointCloudMov：j，分，down，face，gpu
int IcpCu::getMove(float* h_input, unsigned int size) {


	float* hX, * hY, * hZ;
	hX = new float[size];
	hY = new float[size];
	hZ = new float[size];


	//if (withTrans) {
	//	for (int i = 0; i < size; i++) {
	//		hX[i] = h_input[i * 3 + 0] - baseTrans[0];
	//		hY[i] = h_input[i * 3 + 1] - baseTrans[1];
	//		hZ[i] = h_input[i * 3 + 2] - baseTrans[2];
	//	}
	//}
	//else {
		for (int i = 0; i < size; i++) {
			hX[i] = h_input[i * 3 + 0];
			hY[i] = h_input[i * 3 + 1];
			hZ[i] = h_input[i * 3 + 2];
		}
	//}


	if (pointMoveSpace < size) {
		cudaMemcpy(d_pointCloudMX, hX, sizeof(float) * pointMoveSpace, cudaMemcpyHostToDevice);
		cudaMemcpy(d_pointCloudMY, hY, sizeof(float) * pointMoveSpace, cudaMemcpyHostToDevice);
		cudaMemcpy(d_pointCloudMZ, hZ, sizeof(float) * pointMoveSpace, cudaMemcpyHostToDevice);
		this->pointMovSize = pointMoveSpace;
	}
	else {
		cudaMemcpy(d_pointCloudMX, hX, sizeof(float) * size, cudaMemcpyHostToDevice);
		cudaMemcpy(d_pointCloudMY, hY, sizeof(float) * size, cudaMemcpyHostToDevice);
		cudaMemcpy(d_pointCloudMZ, hZ, sizeof(float) * size, cudaMemcpyHostToDevice);
		this->pointMovSize = size;
	}

	
	delete[] hX;
	delete[] hY;
	delete[] hZ;


	return 0;
}




// 一堆函数的合并
int IcpCu::prepareFix() {
	computeFFDistance();
	sortDistance();
	computeNormal();
	return 0;
}




int IcpCu::computeFFDistance() {

	dim3 GridDim((this->pointFixSize + 31) / 32, (this->pointFixSize + 31) / 32, 1);
	dim3 BlockDim(32, 1, 1);

	cuComputeDistance3FF << <GridDim, BlockDim >> >(d_pointCloudFX, d_pointCloudFY, d_pointCloudFZ, pointFixSize, pointFixSpace, d_dist2Matrix);

	return 0;
}
int IcpCu::sortDistance() {

	float* maxDist;
	cudaMalloc(&maxDist, sizeof(float) * pointFixSize);

	cuSortDistance3 << <this->pointFixSize, 32 >> > (d_index, d_dist2Matrix, kNearNum, pointFixSize, pointFixSpace, maxDist);

	return 0;
}
int IcpCu::computeNormal() {

	dim3 sortGridDim((this->pointFixSize + 31) / 32);
	dim3 sortBlockDim(32, 1, 1);

	computeNormal2 << <sortGridDim, sortBlockDim >> > (d_pointCloudFX, d_pointCloudFY, d_pointCloudFZ, pointFixSize, d_pointClNormalX, d_pointClNormalY, d_pointClNormalZ, d_index, kNearNum);

	return 0;
}

int IcpCu::findMoveMin() {

	//float* ttt4 = new float[pointMoveSize];
	//cudaMemcpy(ttt4, d_pointCloudMX, pointMoveSize * sizeof(float), cudaMemcpyDeviceToHost);
	//float* ttt5 = new float[pointFixSize];
	//cudaMemcpy(ttt5, d_pointCloudFX, pointFixSize * sizeof(float), cudaMemcpyDeviceToHost);
	//float* ttt6 = new float[pointMoveSize];
	//cudaMemcpy(ttt6, d_pointCloudMZ, pointMoveSize * sizeof(float), cudaMemcpyDeviceToHost);

	cuFindQMin5 << < this->pointMovSize, 256 >> >(d_pointCloudFX, d_pointCloudFY, d_pointCloudFZ, pointFixSize, d_pointCloudMX, d_pointCloudMY, d_pointCloudMZ, pointMovSize, d_nearestIndex, minDist);

	//float* ttt4 = new float[pointMoveSize];
	//cudaMemcpy(ttt4, minDist, pointMoveSize * sizeof(float), cudaMemcpyDeviceToHost);

	//cudaMemcpy(ttt4, d_pointCloudMX, pointMoveSize * sizeof(float), cudaMemcpyDeviceToHost);
	//cudaMemcpy(ttt5, d_pointCloudFX, pointFixSize * sizeof(float), cudaMemcpyDeviceToHost);
	//cudaMemcpy(ttt6, d_pointCloudMZ, pointMoveSize * sizeof(float), cudaMemcpyDeviceToHost);

	return 0;
}

int IcpCu::ComputeAb() {

	dim3 sortGridDim((this->pointMovSize + 31) / 32);
	dim3 sortBlockDim(32, 1, 1);

	cuComputeAb3 << <sortGridDim, sortBlockDim >> > (A, b, d_pointCloudFX, d_pointCloudFY, d_pointCloudFZ, pointFixSize,
		d_pointCloudMX, d_pointCloudMY, d_pointCloudMZ, pointMovSize, d_pointClNormalX, d_pointClNormalY, d_pointClNormalZ, d_nearestIndex, minDist, distThres);

	return 0;
}




// RT, bT
int IcpCu::icpRegistraion() {

	cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, 6, 6, pointMovSize, &alpha, A, 6, A, 6, &beta, A_, 6);
	cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, 6, 1, pointMovSize, &alpha, A, 6, b, 1, &beta, b_, 6);

	cusolverDnSgetrf(cusolverH, 6, 6, A_, 6, d_work, d_Ipiv, d_info);

	cusolverDnSgetrs(cusolverH, CUBLAS_OP_N, 6, 1, A_, 6, d_Ipiv, b_, 6, d_info);

	svd3_b_ << <1, 1 >> > (b_, UVW);

	cublasSgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 3, 3, 3, &alpha, &UVW[0], 3, &UVW[12], 3, &beta, R_, 3);

	cudaMemcpy(RT, R_, 9 * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(bT, b_, 6 * sizeof(float), cudaMemcpyDeviceToHost); 


	return 0;
}




int IcpCu::transform(float* &d_outputX, float* &d_outputY, float* &d_outputZ, unsigned int size, float *R, float* b) {

	dim3 sortGridDim((this->pointMovSize + 31) / 32);
	dim3 sortBlockDim(32, 1, 1);

	cuTransform << <sortGridDim, sortBlockDim >> >(d_outputX, d_outputY, d_outputZ, size, R, b, d_combineCloudX, d_combineCloudY, d_combineCloudZ);
	cudaMemcpy(d_outputX, d_combineCloudX, size * sizeof(float), cudaMemcpyDeviceToDevice);
	cudaMemcpy(d_outputY, d_combineCloudY, size * sizeof(float), cudaMemcpyDeviceToDevice);
	cudaMemcpy(d_outputZ, d_combineCloudZ, size * sizeof(float), cudaMemcpyDeviceToDevice);

	return 0;
}

int IcpCu::transform(float* &d_outputX, float* &d_outputY, float* &d_outputZ, unsigned int size, 
	float R11, float R12, float R13, float R21, float R22, float R23, float R31, float R32, float R33, 
	float b1, float b2, float b3) {

	dim3 sortGridDim((size + 31) / 32);
	dim3 sortBlockDim(32, 1, 1);

	cuTransform2 << <sortGridDim, sortBlockDim >> > (d_outputX, d_outputY, d_outputZ, size,
		R11,  R12,  R13,  R21,  R22,  R23,  R31,  R32,  R33,
		 b1,  b2,  b3, d_combineCloudX, d_combineCloudY, d_combineCloudZ);
	cudaMemcpy(d_outputX, d_combineCloudX, size * sizeof(float), cudaMemcpyDeviceToDevice);
	cudaMemcpy(d_outputY, d_combineCloudY, size * sizeof(float), cudaMemcpyDeviceToDevice);
	cudaMemcpy(d_outputZ, d_combineCloudZ, size * sizeof(float), cudaMemcpyDeviceToDevice);
	return 0;
}

int IcpCu::transform(float*& d_outputX, float*& d_outputY, float*& d_outputZ, unsigned int size,
	Matrix & trans) {
	transform(d_outputX, d_outputY, d_outputZ, size,
		trans.val[0][0], trans.val[0][1], trans.val[0][2], trans.val[1][0], trans.val[1][1], trans.val[1][2],
		trans.val[2][0], trans.val[2][1], trans.val[2][2], trans.val[0][3], trans.val[1][3], trans.val[2][3]);
	return 0;
}

int IcpCu::verifyRb() {

	Matrix R(3, 3);
	R.val[0][0] = RT[0]; R.val[0][1] = RT[3]; R.val[0][2] = RT[6];
	R.val[1][0] = RT[1]; R.val[1][1] = RT[4]; R.val[1][2] = RT[7];
	R.val[2][0] = RT[2]; R.val[2][1] = RT[5]; R.val[2][2] = RT[8];

	if (R.det() < 0) {
		std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		Matrix B = Matrix::eye(3);
		B.val[2][2] = R.det();
		Matrix U(3, 3);
		Matrix V(3, 3);
		cudaMemcpy(U.val[0], &this->UVW[0], 9 * sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(V.val[0], &this->UVW[12], 9 * sizeof(float), cudaMemcpyDeviceToHost);
		R = V * B*~U;
		R=~R;
		memcpy(RT, R.val[0], 9 * sizeof(float));
		cudaMemcpy(R_, R.val[0], 9 * sizeof(float), cudaMemcpyHostToDevice);
		std::cout << "verify" << std::endl;
	}
	
	return 0;
}




// transMovToFix
int IcpCu::getTransformMatrix() {


	transMovToFix.val[0][0] = RT[0]; transMovToFix.val[0][1] = RT[3]; transMovToFix.val[0][2] = RT[6]; transMovToFix.val[0][3] = bT[3];
	transMovToFix.val[1][0] = RT[1]; transMovToFix.val[1][1] = RT[4]; transMovToFix.val[1][2] = RT[7]; transMovToFix.val[1][3] = bT[4];
	transMovToFix.val[2][0] = RT[2]; transMovToFix.val[2][1] = RT[5]; transMovToFix.val[2][2] = RT[8]; transMovToFix.val[2][3] = bT[5];
	transMovToFix.val[3][0] = 	  0; transMovToFix.val[3][1] =     0; transMovToFix.val[3][2] =     0; transMovToFix.val[3][3] =     1;

	return 0;
}




int IcpCu::combineMove() {

	cudaMemcpy(d_combineCloudX, d_pointCloudFX, pointFixSize * sizeof(float), cudaMemcpyDeviceToDevice);
	cudaMemcpy(d_combineCloudY, d_pointCloudFY, pointFixSize * sizeof(float), cudaMemcpyDeviceToDevice);
	cudaMemcpy(d_combineCloudZ, d_pointCloudFZ, pointFixSize * sizeof(float), cudaMemcpyDeviceToDevice);

	cudaMemcpy(d_combineCloudX + pointFixSize, d_pointCloudMX, pointMovSize * sizeof(float), cudaMemcpyDeviceToDevice);
	cudaMemcpy(d_combineCloudY + pointFixSize, d_pointCloudMY, pointMovSize * sizeof(float), cudaMemcpyDeviceToDevice);
	cudaMemcpy(d_combineCloudZ + pointFixSize, d_pointCloudMZ, pointMovSize * sizeof(float), cudaMemcpyDeviceToDevice);

	combineSize = pointFixSize + pointMovSize;

	return 0;
}

int IcpCu::refreshFix() {

	if (pointFixSpace < combineSize) {
		cudaMemcpy(d_pointCloudFX, d_combineCloudX, pointFixSpace * sizeof(float), cudaMemcpyDeviceToDevice);
		cudaMemcpy(d_pointCloudFY, d_combineCloudY, pointFixSpace * sizeof(float), cudaMemcpyDeviceToDevice);
		cudaMemcpy(d_pointCloudFZ, d_combineCloudZ, pointFixSpace * sizeof(float), cudaMemcpyDeviceToDevice);
	}
	else {
		cudaMemcpy(d_pointCloudFX, d_combineCloudX, combineSize * sizeof(float), cudaMemcpyDeviceToDevice);
		cudaMemcpy(d_pointCloudFY, d_combineCloudY, combineSize * sizeof(float), cudaMemcpyDeviceToDevice);
		cudaMemcpy(d_pointCloudFZ, d_combineCloudZ, combineSize * sizeof(float), cudaMemcpyDeviceToDevice);
	}

	pointFixSize = combineSize;

	return 0;
}

int IcpCu::set0(float* &d_outputX, float* &d_outputY, float* &d_outputZ, unsigned int size) {

	cudaMemset(d_outputX, 0, size * sizeof(float));
	cudaMemset(d_outputY, 0, size * sizeof(float));
	cudaMemset(d_outputZ, 0, size * sizeof(float));

	return 0;
}




// confidence
int IcpCu::ComputeConfidence(float voxelSize, float& confidence) {

	// 堆-指针-地址：float数组
	float* Dist = new float[pointMovSize];

	cudaMemcpy(Dist, minDist, pointMovSize * sizeof(float), cudaMemcpyDeviceToHost);

	int inlierNum = 0;
	for (int i = 0; i < pointMovSize; i++) {
		if (Dist[i] < (voxelSize*0.8) * (voxelSize * 0.8)) {
			inlierNum++;
		}
	}

	confidence = float(inlierNum) / pointMovSize;

	// 堆地址收回
	delete[] Dist;
	return 0;
}


//#include <sstream>
//#include<fstream>
//using namespace std;
//int writeRawFile(std::string fileName, PointCloud<PointXYZ> pointCloud) {
//	std::ofstream outDataFile1;
//	outDataFile1.open(fileName, std::ios::out | std::ios::trunc);
//	//if (outDataFile1 == NULL)
//	//{
//	//	printf("Cannot?open?point?cloud?file.\n");
//	//	return -1;
//	//}
//
//	std::string buf;
//	for (int i = 0; i < pointCloud.size(); i++) {
//		outDataFile1 << pointCloud.points[i].x << " " << pointCloud.points[i].y << " " << pointCloud.points[i].z << endl;
//
//	}
//	return 0;
//}
static int tempIndex = 1;
int IcpCu::ComputeExpression(float voxelSize, float& expression) {

	float* Dist = new float[pointMovSize];
	cudaMemcpy(Dist, minDist, pointMovSize * sizeof(float), cudaMemcpyDeviceToHost);

	PointCloud<PointXYZ> pclTemp;
	gpuToHost(d_pointCloudMX, d_pointCloudMY, d_pointCloudMZ, pointMovSize, pclTemp); 
	PointCloud<PointXYZ> pclTemp2;
	int inlierNum = 0;
	for (int i = 0; i < pointMovSize; i++) {
		if (Dist[i] < (voxelSize + voxelSize) * (voxelSize + voxelSize) && Dist[i] > (voxelSize * 0.8) * (voxelSize * 0.8)) {
			inlierNum++;
			pclTemp2.points.push_back(pclTemp.points[i]);
		}
	}
	//std::stringstream ss4;
	//ss4 << "..\\track_" << tempIndex << "ex.xyz";
	//std::string ss(tempIndex);
	//writeRawFile(ss4.str(), pclTemp2);
	tempIndex++;

	expression = float(inlierNum) /pointMovSize;
	//std::sort(Dist, Dist + pointMovSize);
	//expression = Dist[pointMovSize/2];
	delete[] Dist;
	return 0;
}

int IcpCu::gpuToHost(float*& d_outputX, float*& d_outputY, float*& d_outputZ, unsigned int size, PointCloud<PointXYZ>& pointCloud) {

	pointCloud.points.reserve(4000);
	pointCloud.points.clear();

	float* hX, * hY, * hZ;
	hX = new float[size];
	hY = new float[size];
	hZ = new float[size];

	cudaMemcpy(hX, d_outputX, sizeof(float) * size, cudaMemcpyDeviceToHost);
	cudaMemcpy(hY, d_outputY, sizeof(float) * size, cudaMemcpyDeviceToHost);
	cudaMemcpy(hZ, d_outputZ, sizeof(float) * size, cudaMemcpyDeviceToHost);

	for (int i = 0; i < size; i++) {
		PointXYZ temp;
		temp.x = hX[i];
		temp.y = hY[i];
		temp.z = hZ[i];
		pointCloud.points.push_back(temp);
	}

	delete[] hX;
	delete[] hY;
	delete[] hZ;

	return 0;
}

int IcpCu::initialize(int pointFixSpace, int pointMoveSpace) {
	this->pointFixSpace = pointFixSpace;
	this->pointMoveSpace = pointMoveSpace;

	cudaMalloc(&d_pointCloudFX, sizeof(float) * pointFixSpace);
	cudaMalloc(&d_pointCloudFY, sizeof(float) * pointFixSpace);
	cudaMalloc(&d_pointCloudFZ, sizeof(float) * pointFixSpace);

	cudaMalloc(&d_pointCloudMX, sizeof(float) * pointMoveSpace);
	cudaMalloc(&d_pointCloudMY, sizeof(float) * pointMoveSpace);
	cudaMalloc(&d_pointCloudMZ, sizeof(float) * pointMoveSpace);

	cudaMalloc(&d_combineCloudX, sizeof(float) * (pointMoveSpace + pointFixSpace));
	cudaMalloc(&d_combineCloudY, sizeof(float) * (pointMoveSpace + pointFixSpace));
	cudaMalloc(&d_combineCloudZ, sizeof(float) * (pointMoveSpace + pointFixSpace));

	cudaMalloc(&d_dist2Matrix, this->pointFixSpace * this->pointFixSpace * sizeof(float));

	kNearNum = 16;
	cudaMalloc(&d_index, kNearNum * this->pointMoveSpace * sizeof(unsigned int));
	cudaMalloc(&d_nearestIndex, this->pointMoveSpace * sizeof(unsigned int));

	cudaMalloc(&minDist, sizeof(float) * this->pointMoveSpace);

	cudaMalloc(&d_pointClNormalX, sizeof(float) * this->pointFixSpace);
	cudaMalloc(&d_pointClNormalY, sizeof(float) * this->pointFixSpace);
	cudaMalloc(&d_pointClNormalZ, sizeof(float) * this->pointFixSpace);

	status = CUSOLVER_STATUS_SUCCESS;
	alpha = 1.0;
	beta = 0.0;
	cublasCreate(&handle);
	cusolverH = NULL;
	status = cusolverDnCreate(&cusolverH);

	cudaMalloc(&A, 6 * pointFixSpace * sizeof(float));
	cudaMalloc(&b, pointMoveSpace * sizeof(float));

	cudaMalloc(&A_, 6 * 6 * sizeof(float));
	cudaMalloc(&b_, 6 * sizeof(float));

	lwork = 0; /* size of workspace */
	cudaMalloc((void**)&d_Ipiv, sizeof(int) * 6);
	cudaMalloc((void**)&d_info, sizeof(int));

	cudaMalloc(&UVW, 21 * sizeof(float));
	cudaMalloc(&R_, 9 * sizeof(float));

	status = cusolverDnSgetrf_bufferSize(
		cusolverH,
		6,
		6,
		A_,
		6,
		&lwork);
	cudaMalloc((void**)&d_work, sizeof(float) * 27);

	RT = new float[9]{};
	bT = new float[6]{};

	return 0;
}

int IcpCu::deInitialize() {
	cudaFree(d_pointCloudFX);
	cudaFree(d_pointCloudFY);
	cudaFree(d_pointCloudFZ);

	cudaFree(d_pointCloudMX);
	cudaFree(d_pointCloudMY);
	cudaFree(d_pointCloudMZ);

	cudaFree(d_combineCloudX);
	cudaFree(d_combineCloudY);
	cudaFree(d_combineCloudZ);

	cudaFree(d_pointClNormalX);
	cudaFree(d_pointClNormalY);
	cudaFree(d_pointClNormalZ);

	cudaFree(d_dist2Matrix);
	cudaFree(d_index);
	cudaFree(d_nearestIndex);
	cudaFree(minDist);

	cudaFree(A);
	cudaFree(b);
	cudaFree(A_);
	cudaFree(b_);
	cudaFree(d_Ipiv);
	cudaFree(d_info);
	cudaFree(UVW);
	cudaFree(R_);
	cudaFree(d_work);

	cublasDestroy(handle);
	cusolverDnDestroy(cusolverH);

	delete[] RT;
	delete[] bT;

	return 0;
}

IcpCu::IcpCu()
{
	transMovToFix = Matrix::eye(4);
}

IcpCu::~IcpCu()
{

}
