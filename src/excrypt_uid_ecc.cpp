#include "excrypt.h"

#include <cstdint>
#include <span>

namespace {

// Galois Field GF(2^7) multiplication modulo reduction polynomial x^7 + x + 1 (0x83)
[[nodiscard]] constexpr uint8_t gf_mul(uint8_t a, uint8_t b) noexcept {
    a &= 0x7F;
    b &= 0x7F;
    uint32_t p = 0;
    for (int i = 0; i < 7; ++i) {
        if ((b >> i) & 1) {
            p ^= static_cast<uint32_t>(a) << i;
        }
    }
    for (int i = 12; i >= 7; --i) {
        if ((p >> i) & 1) {
            p ^= 0x83u << (i - 7);
        }
    }
    return static_cast<uint8_t>(p);
}

constexpr uint32_t ECC_LFSR_POLY = 0x360325u;
constexpr uint8_t ALPHA_POS1 = 0x02;
constexpr uint8_t ALPHA_NEG1 = 0x41;
constexpr uint8_t ALPHA_NEG2 = 0x61;

} // namespace

extern "C" {

void XeCryptUidEccEncode(uint8_t cpu_key[16]) {
    if (!cpu_key) {
        return;
    }

    uint32_t ecc_state = 0;
    uint8_t overall_parity = 0;

    for (int k = 0; k < 106; ++k) {
        const uint8_t bit = (cpu_key[k / 8] >> (k % 8)) & 1;
        overall_parity ^= bit;
        ecc_state ^= bit;
        if (ecc_state & 1) {
            ecc_state ^= ECC_LFSR_POLY;
        }
        ecc_state >>= 1;
    }

    for (int k = 106; k < 127; ++k) {
        const uint8_t ecc_bit = ecc_state & 1;
        overall_parity ^= ecc_bit;
        cpu_key[k / 8] = static_cast<uint8_t>((cpu_key[k / 8] & ~(1 << (k % 8))) | (ecc_bit << (k % 8)));
        ecc_state >>= 1;
    }

    cpu_key[15] = static_cast<uint8_t>((cpu_key[15] & 0x7F) | (overall_parity << 7));
}

int XeCryptUidEccDecode(uint8_t cpu_key[16]) {
    if (!cpu_key) {
        return -1;
    }

    uint8_t s1 = 0;
    uint8_t s3 = 0;
    uint8_t s5 = 0;
    uint8_t parity_error = (cpu_key[15] >> 7) & 1;

    for (int k = 0; k < 127; ++k) {
        const uint8_t bit = (cpu_key[k / 8] >> (k % 8)) & 1;
        parity_error ^= bit;
        s1 = gf_mul(s1 ^ bit, ALPHA_POS1);
        s3 = gf_mul(s3 ^ bit, 0x08);
        s5 = gf_mul(s5 ^ bit, 0x20);
    }

    if (parity_error == 0 && s1 == 0 && s3 == 0 && s5 == 0) {
        return 0;
    }

    const uint8_t s1_sq = gf_mul(s1, s1);          // S1^2
    const uint8_t t1 = gf_mul(s1, s1_sq) ^ s3;     // S1^3 + S3
    const uint8_t t2 = gf_mul(s3, s1_sq) ^ s5;     // S3*S1^2 + S5

    if (t1 == 0 && t2 != 0) {
        return -1;
    }

    uint8_t target = 0;
    uint8_t coef_a = 0;
    uint8_t coef_b = 0;
    uint8_t coef_c = 0;
    uint8_t rem_errors = 0;

    if (t1 == 0 && t2 == 0) {
        coef_a = 1;
        coef_b = 0;
        coef_c = 0;
        target = s1;
        rem_errors = (s1 == 0) ? 0 : 1;
    } else {
        coef_a = t1;
        coef_b = t2;
        coef_c = gf_mul(t1, t1) ^ gf_mul(s1, t2); // T1^2 + S1*T2 in GF(2^7)
        target = gf_mul(s1, t1);
        rem_errors = (coef_c == 0) ? 2 : 3;
    }

    uint8_t total_errors = rem_errors;

    if (rem_errors % 2 != parity_error) {
        if (rem_errors == 3) {
            return -1;
        }
        total_errors += 1;
        cpu_key[15] ^= 0x80;
    }

    for (int i = 126; i >= 0; --i) {
        coef_a = gf_mul(coef_a, ALPHA_POS1); // * alpha^+1 (0x02)
        coef_b = gf_mul(coef_b, ALPHA_NEG1); // * alpha^-1 (0x41)
        coef_c = gf_mul(coef_c, ALPHA_NEG2); // * alpha^-2 (0x61)

        if ((coef_a ^ coef_b ^ coef_c) == target) {
            cpu_key[i / 8] ^= static_cast<uint8_t>(1 << (i % 8));
            rem_errors -= 1;
        }
    }

    if (rem_errors == 0) {
        return total_errors;
    }

    return -1;
}

} // extern "C"
