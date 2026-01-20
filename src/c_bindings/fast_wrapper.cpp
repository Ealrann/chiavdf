#include "fast_wrapper.h"

#include <atomic>
#include <mutex>
#include <vector>

#include "../vdf.h"
#include "../create_discriminant.h"

// Runtime configuration knobs required by `parameters.h`.
// These are `extern` variables there, but each binary defines them explicitly.
bool use_divide_table = false;
int gcd_base_bits = 50;
int gcd_128_max_iter = 3;
std::string asmprefix = "cel_";
bool enable_all_instructions = false;

namespace {
std::once_flag init_once;

void init_chiavdf_fast() {
    init_gmp();
    set_rounding_mode();

    // Match the vdf_client runtime selection for AVX2.
    if (hasAVX2()) {
        gcd_base_bits = 63;
        gcd_128_max_iter = 2;
    } else {
        gcd_base_bits = 50;
        gcd_128_max_iter = 3;
    }

    // Ensure we run the one-wesolowski path by default.
    fast_algorithm = false;
    two_weso = false;
    quiet_mode = true;
}

ChiavdfByteArray empty_result() { return ChiavdfByteArray{nullptr, 0}; }

uint64_t get_block(uint64_t i, uint64_t k, uint64_t T, integer& B) {
    integer res = FastPow(2, T - k * (i + 1), B);
    mpz_mul_2exp(res.impl, res.impl, k);
    res = res / B;
    auto res_vector = res.to_vector();
    return res_vector.empty() ? 0 : res_vector[0];
}

class ProgressOneWesolowskiCallback final : public OneWesolowskiCallback {
  public:
    ProgressOneWesolowskiCallback(
        integer& D,
        form& f,
        uint64_t wanted_iter,
        uint64_t progress_interval,
        ChiavdfProgressCallback progress_cb,
        void* progress_user_data)
        : OneWesolowskiCallback(D, f, wanted_iter),
          progress_interval(progress_interval),
          progress_cb(progress_cb),
          progress_user_data(progress_user_data),
          next_progress(progress_interval) {}

    void OnIteration(int type, void* data, uint64_t iteration) override {
        OneWesolowskiCallback::OnIteration(type, data, iteration);

        if (progress_cb == nullptr || progress_interval == 0) {
            return;
        }

        uint64_t done = iteration + 1;
        if (done > wanted_iter) {
            return;
        }

        if (done >= next_progress) {
            progress_cb(next_progress, progress_user_data);
            next_progress += progress_interval;
        }
    }

  private:
    uint64_t progress_interval;
    ChiavdfProgressCallback progress_cb;
    void* progress_user_data;
    uint64_t next_progress;
};

class StreamingOneWesolowskiCallback final : public WesolowskiCallback {
  public:
    StreamingOneWesolowskiCallback(
        integer& D,
        uint64_t wanted_iter,
        uint32_t k,
        uint32_t l,
        uint64_t limit,
        integer& B,
        uint64_t progress_interval,
        ChiavdfProgressCallback progress_cb,
        void* progress_user_data)
        : WesolowskiCallback(D),
          wanted_iter(wanted_iter),
          k(k),
          l(l),
          kl(static_cast<uint64_t>(k) * static_cast<uint64_t>(l)),
          limit(limit),
          B(B),
          progress_interval(progress_interval),
          progress_cb(progress_cb),
          progress_user_data(progress_user_data),
          next_progress(progress_interval) {
        form id = form::identity(D);
        buckets.resize(static_cast<size_t>(l) * (1ULL << k), id);
    }

    void OnIteration(int type, void* data, uint64_t iteration) override {
        iteration++;
        if (iteration > wanted_iter) {
            return;
        }

        if (progress_cb != nullptr && progress_interval != 0 && iteration >= next_progress) {
            progress_cb(next_progress, progress_user_data);
            next_progress += progress_interval;
        }

        if (iteration % kl == 0) {
            uint64_t pos = iteration / kl;
            if (pos < limit) {
                form checkpoint;
                SetForm(type, data, &checkpoint);
                process_checkpoint(pos, checkpoint);
            }
        }

        if (iteration == wanted_iter) {
            SetForm(type, data, &result);
            has_result = true;
        }
    }

    void process_checkpoint(uint64_t i, const form& checkpoint) {
        for (uint32_t j = 0; j < l; j++) {
            uint64_t p = i * static_cast<uint64_t>(l) + static_cast<uint64_t>(j);
            uint64_t needed = static_cast<uint64_t>(k) * (p + 1);
            if (wanted_iter < needed) {
                continue;
            }
            uint64_t b = get_block(p, k, wanted_iter, B);
            nucomp_form(bucket(j, b), bucket(j, b), checkpoint, D, L);
        }
    }

    bool ok() const { return has_result; }

    const form& y() const { return result; }

    form finalize_proof() {
        PulmarkReducer reducer;
        form id = form::identity(D);

        uint64_t k1 = k / 2;
        uint64_t k0 = k - k1;
        form x = id;

        for (int64_t j = static_cast<int64_t>(l) - 1; j >= 0; j--) {
            x = FastPowFormNucomp(x, D, integer(static_cast<uint64_t>(1) << k), L, reducer);

            for (uint64_t b1 = 0; b1 < (1ULL << k1); b1++) {
                form z = id;
                for (uint64_t b0 = 0; b0 < (1ULL << k0); b0++) {
                    nucomp_form(z, z, bucket(static_cast<uint32_t>(j), b1 * (1ULL << k0) + b0), D, L);
                }
                z = FastPowFormNucomp(
                    z,
                    D,
                    integer(static_cast<uint64_t>(b1 * (1ULL << k0))),
                    L,
                    reducer);
                nucomp_form(x, x, z, D, L);
            }

            for (uint64_t b0 = 0; b0 < (1ULL << k0); b0++) {
                form z = id;
                for (uint64_t b1 = 0; b1 < (1ULL << k1); b1++) {
                    nucomp_form(z, z, bucket(static_cast<uint32_t>(j), b1 * (1ULL << k0) + b0), D, L);
                }
                z = FastPowFormNucomp(z, D, integer(b0), L, reducer);
                nucomp_form(x, x, z, D, L);
            }
        }

        reducer.reduce(x);
        return x;
    }

  private:
    form& bucket(uint32_t j, uint64_t b) {
        size_t idx = static_cast<size_t>(j) * (1ULL << k) + static_cast<size_t>(b);
        return buckets[idx];
    }

    const form& bucket(uint32_t j, uint64_t b) const {
        size_t idx = static_cast<size_t>(j) * (1ULL << k) + static_cast<size_t>(b);
        return buckets[idx];
    }

    uint64_t wanted_iter;
    uint32_t k;
    uint32_t l;
    uint64_t kl;
    uint64_t limit;
    integer B;
    uint64_t progress_interval;
    ChiavdfProgressCallback progress_cb;
    void* progress_user_data;
    uint64_t next_progress;

    std::vector<form> buckets;
    form result;
    bool has_result = false;
};
} // namespace

extern "C" ChiavdfByteArray chiavdf_prove_one_weso_fast(
    const uint8_t* challenge_hash,
    size_t challenge_size,
    const uint8_t* x_s,
    size_t x_s_size,
    size_t discriminant_size_bits,
    uint64_t num_iterations) {
    return chiavdf_prove_one_weso_fast_with_progress(
        challenge_hash,
        challenge_size,
        x_s,
        x_s_size,
        discriminant_size_bits,
        num_iterations,
        /*progress_interval=*/0,
        /*progress_cb=*/nullptr,
        /*progress_user_data=*/nullptr);
}

extern "C" ChiavdfByteArray chiavdf_prove_one_weso_fast_with_progress(
    const uint8_t* challenge_hash,
    size_t challenge_size,
    const uint8_t* x_s,
    size_t x_s_size,
    size_t discriminant_size_bits,
    uint64_t num_iterations,
    uint64_t progress_interval,
    ChiavdfProgressCallback progress_cb,
    void* progress_user_data) {
    try {
        std::call_once(init_once, init_chiavdf_fast);

        if (challenge_hash == nullptr || challenge_size == 0 || x_s == nullptr || x_s_size == 0) {
            return empty_result();
        }
        if (num_iterations == 0) {
            return empty_result();
        }

        std::vector<uint8_t> challenge_hash_bytes(challenge_hash, challenge_hash + challenge_size);
        integer D = CreateDiscriminant(challenge_hash_bytes, static_cast<int>(discriminant_size_bits));
        integer L = root(-D, 4);

        form x = DeserializeForm(D, x_s, x_s_size);

        std::atomic<bool> stopped(false);
        ProgressOneWesolowskiCallback weso(
            D,
            x,
            num_iterations,
            progress_interval,
            progress_cb,
            progress_user_data);

        // Run the fast repeated-squaring engine to `num_iterations`.
        // The callback stores all intermediates needed for the proof.
        FastStorage* fast_storage = nullptr;
        repeated_square(num_iterations, x, D, L, &weso, fast_storage, stopped);

        // Now generate the compact proof from the stored intermediates.
        Proof proof = ProveOneWesolowski(num_iterations, D, x, &weso, stopped);
        if (proof.y.empty() || proof.proof.empty()) {
            return empty_result();
        }

        const size_t total = proof.y.size() + proof.proof.size();
        uint8_t* out = new uint8_t[total];
        std::copy(proof.y.begin(), proof.y.end(), out);
        std::copy(proof.proof.begin(), proof.proof.end(), out + proof.y.size());
        return ChiavdfByteArray{out, total};
    } catch (...) {
        return empty_result();
    }
}

extern "C" ChiavdfByteArray chiavdf_prove_one_weso_fast_streaming(
    const uint8_t* challenge_hash,
    size_t challenge_size,
    const uint8_t* x_s,
    size_t x_s_size,
    const uint8_t* y_ref_s,
    size_t y_ref_s_size,
    size_t discriminant_size_bits,
    uint64_t num_iterations) {
    return chiavdf_prove_one_weso_fast_streaming_with_progress(
        challenge_hash,
        challenge_size,
        x_s,
        x_s_size,
        y_ref_s,
        y_ref_s_size,
        discriminant_size_bits,
        num_iterations,
        /*progress_interval=*/0,
        /*progress_cb=*/nullptr,
        /*progress_user_data=*/nullptr);
}

extern "C" ChiavdfByteArray chiavdf_prove_one_weso_fast_streaming_with_progress(
    const uint8_t* challenge_hash,
    size_t challenge_size,
    const uint8_t* x_s,
    size_t x_s_size,
    const uint8_t* y_ref_s,
    size_t y_ref_s_size,
    size_t discriminant_size_bits,
    uint64_t num_iterations,
    uint64_t progress_interval,
    ChiavdfProgressCallback progress_cb,
    void* progress_user_data) {
    try {
        std::call_once(init_once, init_chiavdf_fast);

        if (challenge_hash == nullptr || challenge_size == 0 || x_s == nullptr || x_s_size == 0 ||
            y_ref_s == nullptr || y_ref_s_size == 0) {
            return empty_result();
        }
        if (num_iterations == 0) {
            return empty_result();
        }

        std::vector<uint8_t> challenge_hash_bytes(challenge_hash, challenge_hash + challenge_size);
        integer D = CreateDiscriminant(challenge_hash_bytes, static_cast<int>(discriminant_size_bits));
        integer L = root(-D, 4);

        form x = DeserializeForm(D, x_s, x_s_size);
        form y_ref = DeserializeForm(D, y_ref_s, y_ref_s_size);

        uint32_t k;
        uint32_t l;
        if (num_iterations >= (1 << 16)) {
            ApproximateParameters(num_iterations, l, k);
        } else {
            k = 10;
            l = 1;
        }
        if (k == 0) {
            k = 1;
        }
        if (l == 0) {
            l = 1;
        }

        uint64_t kl = static_cast<uint64_t>(k) * static_cast<uint64_t>(l);
        uint64_t limit = num_iterations / kl;
        if (num_iterations % kl) {
            limit++;
        }

        integer B = GetB(D, x, y_ref);

        std::atomic<bool> stopped(false);
        StreamingOneWesolowskiCallback weso(
            D,
            num_iterations,
            k,
            l,
            limit,
            B,
            progress_interval,
            progress_cb,
            progress_user_data);

        weso.process_checkpoint(/*i=*/0, x);

        FastStorage* fast_storage = nullptr;
        repeated_square(num_iterations, x, D, L, &weso, fast_storage, stopped);

        if (!weso.ok()) {
            return empty_result();
        }
        if (!(weso.y() == y_ref)) {
            return empty_result();
        }

        form proof_form = weso.finalize_proof();

        int d_bits = D.num_bits();
        std::vector<unsigned char> y_serialized = SerializeForm(y_ref, d_bits);
        std::vector<unsigned char> proof_serialized = SerializeForm(proof_form, d_bits);

        if (y_serialized.empty() || proof_serialized.empty()) {
            return empty_result();
        }

        const size_t total = y_serialized.size() + proof_serialized.size();
        uint8_t* out = new uint8_t[total];
        std::copy(y_serialized.begin(), y_serialized.end(), out);
        std::copy(proof_serialized.begin(), proof_serialized.end(), out + y_serialized.size());
        return ChiavdfByteArray{out, total};
    } catch (...) {
        return empty_result();
    }
}

extern "C" void chiavdf_free_byte_array(ChiavdfByteArray array) { delete[] array.data; }
