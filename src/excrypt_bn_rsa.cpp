#include "excrypt.h"
#include "fp61.hpp"

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <bcrypt.h>
#endif

#include <algorithm>
#include <cstring>
#include <memory>

namespace {

    int bn_cmp(const uint64_t* a, const uint64_t* b, size_t n) {
        for (size_t i = n; i-- > 0;) {
            if (a[i] < b[i]) {
                return -1;
            }
            if (a[i] > b[i]) {
                return 1;
            }
        }
        return 0;
    }

    bool bn_is_zero(const uint64_t* a, size_t n) {
        for (size_t i = 0; i < n; ++i) {
            if (a[i] != 0) {
                return false;
            }
        }
        return true;
    }

    bool bn_is_one(const uint64_t* a, size_t n) {
        if (n == 0 || a[0] != 1) {
            return false;
        }
        for (size_t i = 1; i < n; ++i) {
            if (a[i] != 0) {
                return false;
            }
        }
        return true;
    }

    bool bn_gt_one(const uint64_t* a, size_t n) {
        for (size_t i = 1; i < n; ++i) {
            if (a[i] != 0) {
                return true;
            }
        }
        return n > 0 && a[0] > 1;
    }

    uint64_t bn_add_inplace(uint64_t* a, const uint64_t* b, size_t n) {
        uint64_t carry = 0;
        for (size_t i = 0; i < n; ++i) {
            const uint64_t ai = a[i];
            const uint64_t sum1 = ai + b[i];
            const uint64_t carry1 = sum1 < ai ? 1ULL : 0ULL;
            const uint64_t sum2 = sum1 + carry;
            const uint64_t carry2 = sum2 < sum1 ? 1ULL : 0ULL;
            a[i] = sum2;
            carry = carry1 + carry2;
        }
        return carry;
    }

    uint64_t bn_sub_inplace(uint64_t* a, const uint64_t* b, size_t n) {
        uint64_t borrow = 0;
        for (size_t i = 0; i < n; ++i) {
            const uint64_t ai = a[i];
            const uint64_t bi = b[i];
            const uint64_t diff1 = ai - bi;
            const uint64_t borrow1 = ai < bi ? 1ULL : 0ULL;
            const uint64_t diff2 = diff1 - borrow;
            const uint64_t borrow2 = diff1 < borrow ? 1ULL : 0ULL;
            a[i] = diff2;
            borrow = borrow1 | borrow2;
        }
        return borrow;
    }

    void bn_add_mod(const uint64_t* a, const uint64_t* b, uint64_t* out, const uint64_t* mod,
                    size_t n) {
        uint64_t tmp[32]{};
        std::memcpy(tmp, a, n * sizeof(uint64_t));
        const uint64_t carry = bn_add_inplace(tmp, b, n);
        if (carry != 0 || bn_cmp(tmp, mod, n) >= 0) {
            bn_sub_inplace(tmp, mod, n);
        }
        std::memcpy(out, tmp, n * sizeof(uint64_t));
    }

    void bn_sub_mod(const uint64_t* a, const uint64_t* b, uint64_t* out, const uint64_t* mod,
                    size_t n) {
        uint64_t tmp[32]{};
        if (bn_cmp(a, b, n) >= 0) {
            std::memcpy(tmp, a, n * sizeof(uint64_t));
            bn_sub_inplace(tmp, b, n);
        } else {
            std::memcpy(tmp, mod, n * sizeof(uint64_t));
            bn_sub_inplace(tmp, b, n);
            bn_add_inplace(tmp, a, n);
        }
        std::memcpy(out, tmp, n * sizeof(uint64_t));
    }

    void bn_mul(const uint64_t* a, size_t an, const uint64_t* b, size_t bn, uint64_t* out) {
        std::memset(out, 0, (an + bn) * sizeof(uint64_t));
        for (size_t i = 0; i < an; ++i) {
            uint64_t carry = 0;
            for (size_t j = 0; j < bn; ++j) {
                uint64_t hi = 0;
                uint64_t lo = 0;
                CAT_MUL128(hi, lo, a[i], b[j]);
                const uint64_t prev = out[i + j];
                lo += prev;
                if (lo < prev) {
                    ++hi;
                }
                lo += carry;
                if (lo < carry) {
                    ++hi;
                }
                out[i + j] = lo;
                carry = hi;
            }
            out[i + bn] = carry;
        }
    }

    void bn_mul_scalar(const uint64_t* a, size_t an, uint64_t s, uint64_t* out) {
        uint64_t carry = 0;
        for (size_t i = 0; i < an; ++i) {
            uint64_t hi = 0;
            uint64_t lo = 0;
            CAT_MUL128(hi, lo, a[i], s);
            lo += carry;
            if (lo < carry) {
                ++hi;
            }
            out[i] = lo;
            carry = hi;
        }
        out[an] = carry;
    }

    void bn_mod(const uint64_t* a, size_t an, const uint64_t* mod, size_t mn, uint64_t* rem) {
        std::memset(rem, 0, mn * sizeof(uint64_t));
        if (bn_is_zero(mod, mn)) {
            return;
        }
        for (size_t i = an; i-- > 0;) {
            for (int bit = 63; bit >= 0; --bit) {
                uint64_t carry = (a[i] >> bit) & 1ULL;
                for (size_t k = 0; k < mn; ++k) {
                    const uint64_t next_carry = rem[k] >> 63;
                    rem[k] = (rem[k] << 1) | carry;
                    carry = next_carry;
                }
                if (carry != 0 || bn_cmp(rem, mod, mn) >= 0) {
                    bn_sub_inplace(rem, mod, mn);
                }
            }
        }
    }

    void mont_mul(const uint64_t* a, const uint64_t* b, uint64_t* out, uint64_t inv,
                  const uint64_t* mod, uint32_t n) {
        if (bn_is_zero(a, n) || bn_is_zero(b, n)) {
            std::memset(out, 0, n * sizeof(uint64_t));
            return;
        }
        ExCryptBnQwNeModMul(a, b, out, inv, mod, n);
        if (bn_cmp(out, mod, n) >= 0) {
            bn_sub_inplace(out, mod, n);
        }
    }

    void mont_r_and_r2(const uint64_t* mod, uint32_t n, uint64_t* r_mod, uint64_t* r2_mod) {
        std::memset(r_mod, 0, n * sizeof(uint64_t));
        r_mod[0] = 1;
        const size_t bits = static_cast<size_t>(n) * 64;
        for (size_t i = 0; i < bits; ++i) {
            bn_add_mod(r_mod, r_mod, r_mod, mod, n);
        }
        std::memcpy(r2_mod, r_mod, n * sizeof(uint64_t));
        for (size_t i = 0; i < bits; ++i) {
            bn_add_mod(r2_mod, r2_mod, r2_mod, mod, n);
        }
    }

    void mont_pow(const uint64_t* base_bar, const uint64_t* exp, size_t exp_limbs,
                  const uint64_t* r_mod, uint64_t inv, const uint64_t* mod, uint32_t n,
                  uint64_t* out_bar) {
        uint64_t acc[32]{};
        std::memcpy(acc, r_mod, n * sizeof(uint64_t));

        size_t top_limb = exp_limbs;
        while (top_limb > 0 && exp[top_limb - 1] == 0) {
            --top_limb;
        }
        if (top_limb == 0) {
            std::memcpy(out_bar, acc, n * sizeof(uint64_t));
            return;
        }

        int top_bit = 63;
        while (top_bit > 0 && ((exp[top_limb - 1] >> top_bit) & 1ULL) == 0) {
            --top_bit;
        }

        for (size_t i = top_limb; i-- > 0;) {
            const int start_bit = (i == top_limb - 1) ? top_bit : 63;
            for (int bit = start_bit; bit >= 0; --bit) {
                mont_mul(acc, acc, acc, inv, mod, n);
                if (((exp[i] >> bit) & 1ULL) != 0) {
                    mont_mul(acc, base_bar, acc, inv, mod, n);
                }
            }
        }
        std::memcpy(out_bar, acc, n * sizeof(uint64_t));
    }

    void bn_mod_pow(const uint64_t* base, const uint64_t* exp, size_t exp_limbs,
                    const uint64_t* mod, uint32_t n, uint64_t* out) {
        uint64_t r_mod[32]{};
        uint64_t r2_mod[32]{};
        mont_r_and_r2(mod, n, r_mod, r2_mod);
        const uint64_t inv = ExCryptBnQwNeModInv(mod[0]);

        uint64_t base_bar[32]{};
        mont_mul(base, r2_mod, base_bar, inv, mod, n);

        uint64_t res_bar[32]{};
        mont_pow(base_bar, exp, exp_limbs, r_mod, inv, mod, n, res_bar);

        uint64_t one[32]{};
        one[0] = 1;
        mont_mul(res_bar, one, out, inv, mod, n);
    }

    struct RsaPrv2048Ne {
        uint32_t exp{};
        uint64_t n[32]{};
        uint64_t p[16]{};
        uint64_t q[16]{};
        uint64_t dp[16]{};
        uint64_t dq[16]{};
        uint64_t u[16]{};
    };

    bool load_and_validate_rsa2048_prv(const EXCRYPT_RSA* privkey, RsaPrv2048Ne& out) {
        if (!privkey) {
            return false;
        }
        const uint32_t num_digits = _byteswap_ulong(privkey->num_digits);
        const uint32_t exp = _byteswap_ulong(privkey->pub_exponent);
        if (num_digits != 0x20 || (exp != 3 && exp != 0x10001)) {
            return false;
        }
        const auto* key = reinterpret_cast<const EXCRYPT_RSAPRV_2048*>(privkey);
        out.exp = exp;
        ExCryptBnQw_SwapDwQwLeBe(key->modulus, out.n, 32);
        ExCryptBnQw_SwapDwQwLeBe(key->prime1, out.p, 16);
        ExCryptBnQw_SwapDwQwLeBe(key->prime2, out.q, 16);
        ExCryptBnQw_SwapDwQwLeBe(key->exponent1, out.dp, 16);
        ExCryptBnQw_SwapDwQwLeBe(key->exponent2, out.dq, 16);
        ExCryptBnQw_SwapDwQwLeBe(key->coefficient, out.u, 16);

        if ((out.n[31] >> 56) == 0 || !bn_gt_one(out.p, 16) || !bn_gt_one(out.q, 16) ||
            (out.p[0] & 1ULL) == 0 || (out.q[0] & 1ULL) == 0) {
            return false;
        }

        uint64_t pq[32]{};
        bn_mul(out.p, 16, out.q, 16, pq);
        if (bn_cmp(pq, out.n, 32) != 0) {
            return false;
        }

        uint64_t one[16]{};
        one[0] = 1;

        uint64_t p_minus_1[16]{};
        std::memcpy(p_minus_1, out.p, sizeof(p_minus_1));
        bn_sub_inplace(p_minus_1, one, 16);
        uint64_t edp[17]{};
        bn_mul_scalar(out.dp, 16, exp, edp);
        uint64_t rem[16]{};
        bn_mod(edp, 17, p_minus_1, 16, rem);
        if (!bn_is_one(rem, 16)) {
            return false;
        }

        uint64_t q_minus_1[16]{};
        std::memcpy(q_minus_1, out.q, sizeof(q_minus_1));
        bn_sub_inplace(q_minus_1, one, 16);
        uint64_t edq[17]{};
        bn_mul_scalar(out.dq, 16, exp, edq);
        bn_mod(edq, 17, q_minus_1, 16, rem);
        if (!bn_is_one(rem, 16)) {
            return false;
        }

        uint64_t uq[32]{};
        bn_mul(out.u, 16, out.q, 16, uq);
        bn_mod(uq, 32, out.p, 16, rem);
        if (!bn_is_one(rem, 16)) {
            return false;
        }

        return true;
    }

} // namespace

extern "C" {

uint8_t kStaticPrivateExponent1024[] = {
    0x51, 0xEC, 0x1F, 0x9D, 0x56, 0x26, 0xC2, 0xFC, 0x10, 0xA6, 0x67, 0x64, 0xCB, 0x3A, 0x6D, 0x4D,
    0xA1, 0xE7, 0x4E, 0xA8, 0x42, 0xF0, 0xF4, 0xFD, 0xFA, 0x66, 0xEF, 0xC7, 0x8E, 0x10, 0x2F, 0xE4,
    0x1C, 0xA3, 0x1D, 0xD0, 0xCE, 0x39, 0x2E, 0xC3, 0x19, 0x2D, 0xD0, 0x58, 0x74, 0x79, 0xAC, 0x08,
    0xE7, 0x90, 0xC1, 0xAC, 0x2D, 0xC6, 0xEB, 0x47, 0xE8, 0x3D, 0xCF, 0x4C, 0x6D, 0xFF, 0x51, 0x65,
    0xD4, 0x6E, 0xBD, 0x0F, 0x15, 0x79, 0x37, 0x95, 0xC4, 0xAF, 0x90, 0x9E, 0x2B, 0x50, 0x8A, 0x0A,
    0x22, 0x4A, 0xB3, 0x41, 0xE5, 0x89, 0x80, 0x73, 0xCD, 0xFA, 0x21, 0x02, 0xF5, 0xDD, 0x30, 0xDD,
    0x07, 0x2A, 0x6F, 0x34, 0x07, 0x81, 0x97, 0x7E, 0xB2, 0xFB, 0x72, 0xE9, 0xEA, 0xC1, 0x88, 0x39,
    0xAC, 0x48, 0x2B, 0xA8, 0x4D, 0xFC, 0xD7, 0xED, 0x9B, 0xF9, 0xDE, 0xC2, 0x45, 0x93, 0x4C, 0x4C};

BOOL ExCryptBnQwNeRsaPrvCrypt(const uint64_t* input, uint64_t* output, const EXCRYPT_RSA* key) {
#ifndef _WIN32
    return false;
#else
    if (!input || !output || !key)
        return false;

    uint32_t key_digits = _byteswap_ulong(key->num_digits);
    if (key_digits <= 0 || key_digits > 0x40)
        return false;

    // TODO: currently only works with 1024-bit/16-digit keys, due to kStaticPrivateExponent1024
    // above Need to replace that with some way of calculating it from the rest of the privkey
    // instead
    if (key_digits != 16)
        return false;

    // TODO: check this
    // if (ExCryptBnQwNeCompare(input, xecrypt_modulus, key_digits) >= 0)
    //  return false;

    uint32_t modulus_size = key_digits * 8; // = 16
    uint32_t prime_count = key_digits / 2;
    uint32_t prime_size = prime_count * 8;

    // Convert XECRYPT blob into BCrypt format
    ULONG key_size = sizeof(BCRYPT_RSAKEY_BLOB) + sizeof(uint32_t) + // exponent
                     modulus_size +                                  // modulus
                     prime_size +                                    // prime1
                     prime_size +                                    // prime2
                     prime_size +                                    // exponent1
                     prime_size +                                    // exponent2
                     prime_size +                                    // coefficient
                     modulus_size;                                   // private exponent
    auto key_buf = std::make_unique<uint8_t[]>(key_size);
    auto* key_header = reinterpret_cast<BCRYPT_RSAKEY_BLOB*>(key_buf.get());

    key_header->Magic = BCRYPT_RSAFULLPRIVATE_MAGIC;
    key_header->BitLength = modulus_size * 8;
    key_header->cbPublicExp = sizeof(uint32_t);
    key_header->cbModulus = modulus_size;
    key_header->cbPrime1 = key_header->cbPrime2 = prime_size;

    // Copy in exponent/modulus, luckily these are BE inside BCrypt blob
    uint32_t* key_exponent = reinterpret_cast<uint32_t*>(&key_header[1]);
    *key_exponent = key->pub_exponent;

    // ...except other fields need to be reversed in 64-bit chunks for BCrypt to
    // make use of them properly for some reason
    uint64_t* key_modulus = reinterpret_cast<uint64_t*>(&key_exponent[1]);
    auto* xecrypt_modulus = reinterpret_cast<const uint64_t*>(&key[1]);
    std::reverse_copy(xecrypt_modulus, xecrypt_modulus + key_digits, key_modulus);

    uint64_t* key_prime1 = reinterpret_cast<uint64_t*>(&key_modulus[key_digits]);
    auto* xecrypt_prime1 = reinterpret_cast<const uint64_t*>(&xecrypt_modulus[key_digits]);

    std::reverse_copy(xecrypt_prime1, xecrypt_prime1 + (prime_count), key_prime1);

    uint64_t* key_prime2 = reinterpret_cast<uint64_t*>(&key_prime1[prime_count]);
    auto* xecrypt_prime2 = reinterpret_cast<const uint64_t*>(&xecrypt_prime1[prime_count]);

    std::reverse_copy(xecrypt_prime2, xecrypt_prime2 + prime_count, key_prime2);

    uint64_t* key_exponent1 = reinterpret_cast<uint64_t*>(&key_prime2[prime_count]);
    auto* xecrypt_exponent1 = reinterpret_cast<const uint64_t*>(&xecrypt_prime2[prime_count]);

    std::reverse_copy(xecrypt_exponent1, xecrypt_exponent1 + prime_count, key_exponent1);

    uint64_t* key_exponent2 = reinterpret_cast<uint64_t*>(&key_exponent1[prime_count]);
    auto* xecrypt_exponent2 = reinterpret_cast<const uint64_t*>(&xecrypt_exponent1[prime_count]);

    std::reverse_copy(xecrypt_exponent2, xecrypt_exponent2 + prime_count, key_exponent2);

    uint64_t* key_coefficient = reinterpret_cast<uint64_t*>(&key_exponent2[prime_count]);
    auto* xecrypt_coefficient = reinterpret_cast<const uint64_t*>(&xecrypt_exponent2[prime_count]);

    std::reverse_copy(xecrypt_coefficient, xecrypt_coefficient + prime_count, key_coefficient);

    uint64_t* key_privexponent = reinterpret_cast<uint64_t*>(&key_coefficient[prime_count]);

    // X360 uses a static private exponent / "D" value
    std::memcpy(key_privexponent, kStaticPrivateExponent1024, 0x80);

    BCRYPT_ALG_HANDLE hAlgorithm = NULL;
    NTSTATUS status =
        BCryptOpenAlgorithmProvider(&hAlgorithm, BCRYPT_RSA_ALGORITHM, MS_PRIMITIVE_PROVIDER, 0);

    if (!BCRYPT_SUCCESS(status))
        return false;

    auto* buf = key_buf.get();
    BCRYPT_KEY_HANDLE hKey = NULL;
    status =
        BCryptImportKeyPair(hAlgorithm, NULL, BCRYPT_RSAFULLPRIVATE_BLOB, &hKey, buf, key_size, 0);

    if (!BCRYPT_SUCCESS(status)) {
        if (hAlgorithm)
            BCryptCloseAlgorithmProvider(hAlgorithm, 0);

        return false;
    }

    // Byteswap & reverse the input into output, as BCrypt wants MSB first
    uint8_t* output_bytes = reinterpret_cast<uint8_t*>(output);
    ExCryptBnQw_SwapDwQwLeBe(input, output, key_digits);
    std::reverse(output_bytes, output_bytes + modulus_size);

    // BCryptDecrypt only works with private keys, fortunately BCryptEncrypt
    // performs the right actions needed for us to decrypt the input
    ULONG result_size = 0;

    status = BCryptDecrypt(hKey, output_bytes, modulus_size, nullptr, nullptr, 0, output_bytes,
                           modulus_size, &result_size, BCRYPT_PAD_NONE);

    // assert(result_size == modulus_size);

    if (BCRYPT_SUCCESS(status)) {
        // Reverse data & byteswap again so data is as game expects
        std::reverse(output_bytes, output_bytes + modulus_size);
        ExCryptBnQw_SwapDwQwLeBe(output, output, key_digits);
    }

    if (hKey)
        BCryptDestroyKey(hKey);

    if (hAlgorithm)
        BCryptCloseAlgorithmProvider(hAlgorithm, 0);

    return BCRYPT_SUCCESS(status);
#endif
}

BOOL ExCryptBnQwNeRsaPubCrypt(const uint64_t* input, uint64_t* output, const EXCRYPT_RSA* key) {
#ifndef _WIN32
    return false;
#else
    if (!input || !output || !key)
        return false;

    uint32_t key_digits = _byteswap_ulong(key->num_digits);
    if (key_digits <= 0 || key_digits > 0x40)
        return false;

    auto* xecrypt_modulus = reinterpret_cast<const uint64_t*>(&key[1]);

    // TODO: check this
    // if (ExCryptBnQwNeCompare(input, xecrypt_modulus, key_digits) >= 0)
    //  return false;

    uint32_t modulus_size = key_digits * 8;

    // Convert XECRYPT blob into BCrypt format
    ULONG key_size = sizeof(BCRYPT_RSAKEY_BLOB) + sizeof(uint32_t) + modulus_size;
    auto key_buf = std::make_unique<uint8_t[]>(key_size);

    auto* key_header = reinterpret_cast<BCRYPT_RSAKEY_BLOB*>(key_buf.get());
    key_header->Magic = BCRYPT_RSAPUBLIC_MAGIC;
    key_header->BitLength = modulus_size * 8;
    key_header->cbPublicExp = sizeof(uint32_t);
    key_header->cbModulus = modulus_size;
    key_header->cbPrime1 = key_header->cbPrime2 = 0;

    // Copy in exponent/modulus, luckily these are BE inside BCrypt blob
    uint32_t* key_exponent = reinterpret_cast<uint32_t*>(&key_header[1]);
    *key_exponent = key->pub_exponent;

    // ...except modulus needs to be reversed in 64-bit chunks for BCrypt to make
    // use of it properly for some reason
    auto* key_modulus = reinterpret_cast<uint64_t*>(&key_exponent[1]);
    std::reverse_copy(xecrypt_modulus, xecrypt_modulus + key_digits, key_modulus);

    BCRYPT_ALG_HANDLE hAlgorithm = NULL;
    NTSTATUS status =
        BCryptOpenAlgorithmProvider(&hAlgorithm, BCRYPT_RSA_ALGORITHM, MS_PRIMITIVE_PROVIDER, 0);

    if (!BCRYPT_SUCCESS(status))
        return false;

    BCRYPT_KEY_HANDLE hKey = NULL;
    status = BCryptImportKeyPair(hAlgorithm, NULL, BCRYPT_RSAPUBLIC_BLOB, &hKey, key_buf.get(),
                                 key_size, 0);

    if (!BCRYPT_SUCCESS(status)) {
        if (hAlgorithm)
            BCryptCloseAlgorithmProvider(hAlgorithm, 0);

        return false;
    }

    // Byteswap & reverse the input into output, as BCrypt wants MSB first
    uint8_t* output_bytes = reinterpret_cast<uint8_t*>(output);
    ExCryptBnQw_SwapDwQwLeBe(input, output, key_digits);
    std::reverse(output_bytes, output_bytes + modulus_size);

    // BCryptDecrypt only works with private keys, fortunately BCryptEncrypt
    // performs the right actions needed for us to decrypt the input
    ULONG result_size = 0;
    status = BCryptEncrypt(hKey, output_bytes, modulus_size, nullptr, nullptr, 0, output_bytes,
                           modulus_size, &result_size, BCRYPT_PAD_NONE);

    // assert(result_size == modulus_size);

    if (BCRYPT_SUCCESS(status)) {
        // Reverse data & byteswap again so data is as game expects
        std::reverse(output_bytes, output_bytes + modulus_size);
        ExCryptBnQw_SwapDwQwLeBe(output, output, key_digits);
    }

    if (hKey)
        BCryptDestroyKey(hKey);

    if (hAlgorithm)
        BCryptCloseAlgorithmProvider(hAlgorithm, 0);

    return BCRYPT_SUCCESS(status);
#endif
}

BOOL ExCryptBnQwBeRsaPrvValidate(const EXCRYPT_RSA* privkey) {
    RsaPrv2048Ne key{};
    return load_and_validate_rsa2048_prv(privkey, key) ? 1 : 0;
}

BOOL ExCryptBnQwBeSigCreate(EXCRYPT_SIG* sig, const uint8_t* hash, const uint8_t* salt,
                            const EXCRYPT_RSA* privkey) {
    if (!sig || !hash || !salt) {
        return 0;
    }
    RsaPrv2048Ne key{};
    if (!load_and_validate_rsa2048_prv(privkey, key)) {
        return 0;
    }

    EXCRYPT_SIG formatted{};
    ExCryptBnQwBeSigFormat(&formatted, hash, salt);

    uint64_t msg[32]{};
    ExCryptBnQw_SwapDwQwLeBe(reinterpret_cast<const uint64_t*>(&formatted), msg, 32);
    if (bn_cmp(msg, key.n, 32) >= 0) {
        return 0;
    }

    // Pre-scale by R^(e - 1) mod N to cancel the Montgomery factor divided out by ExCryptBnQwBeSigDifference.
    uint64_t r_mod_n[32]{};
    uint64_t r2_mod_n[32]{};
    mont_r_and_r2(key.n, 32, r_mod_n, r2_mod_n);
    const uint64_t inv_n = ExCryptBnQwNeModInv(key.n[0]);
    const uint64_t scale_exp = key.exp - 1;
    uint64_t scale_bar[32]{};
    mont_pow(r2_mod_n, &scale_exp, 1, r_mod_n, inv_n, key.n, 32, scale_bar);

    uint64_t scaled[32]{};
    mont_mul(msg, scale_bar, scaled, inv_n, key.n, 32);

    // CRT exponentiation mod p and mod q.
    uint64_t s_p[16]{};
    uint64_t s_q[16]{};
    bn_mod(scaled, 32, key.p, 16, s_p);
    bn_mod(scaled, 32, key.q, 16, s_q);

    uint64_t m1[16]{};
    uint64_t m2[16]{};
    bn_mod_pow(s_p, key.dp, 16, key.p, 16, m1);
    bn_mod_pow(s_q, key.dq, 16, key.q, 16, m2);

    uint64_t m2_p[16]{};
    bn_mod(m2, 16, key.p, 16, m2_p);

    uint64_t diff[16]{};
    bn_sub_mod(m1, m2_p, diff, key.p, 16);

    uint64_t u_diff[32]{};
    bn_mul(key.u, 16, diff, 16, u_diff);

    uint64_t h[16]{};
    bn_mod(u_diff, 32, key.p, 16, h);

    uint64_t sig_ne[32]{};
    bn_mul(h, 16, key.q, 16, sig_ne);

    uint64_t m2_pad[32]{};
    std::memcpy(m2_pad, m2, sizeof(m2));
    if (bn_add_inplace(sig_ne, m2_pad, 32) != 0 || bn_cmp(sig_ne, key.n, 32) >= 0) {
        return 0;
    }

    ExCryptBnQw_SwapDwQwLeBe(sig_ne, reinterpret_cast<uint64_t*>(sig), 32);
    return 1;
}
};
