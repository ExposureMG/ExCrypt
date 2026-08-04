#pragma once

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

uint32_t ExCryptEccGetCurveParameters(uint32_t curve_num, uint32_t** curve_param);

void XeCryptUidEccEncode(uint8_t cpu_key[16]);
int XeCryptUidEccDecode(uint8_t cpu_key[16]);

#ifdef __cplusplus
}
#endif
