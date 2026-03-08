#pragma once


#ifdef __cplusplus
extern "C"
{
#endif

typedef struct DRAFrameInfo
{
	int nChannels;
	int nSampleRate;
	int nFrameSize;
}DRAFrameInfo;

void* DRADecCreate(void);
void DRADecDestroy(void* pDecoder);
int DRADecGetFrameInfo(void* pDecoder, DRAFrameInfo* pDRAFrameInfo);
int DRADecSendData(void* pDecoder, unsigned char* pData, int nLen);
int DRADecRecvFrame(void* pDecoder, unsigned char** ppPCMData, int* nPCMLen);

#ifdef __cplusplus
}
#endif