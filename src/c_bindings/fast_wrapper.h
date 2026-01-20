#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint8_t* data;
    size_t length;
} ChiavdfByteArray;

typedef void (*ChiavdfProgressCallback)(uint64_t iters_done, void* user_data);

// Computes a compact (witness_type=0) Wesolowski proof using the fast engine.
//
// On success, returns `y || proof` where:
// - `y` is the serialized output form (typically 100 bytes for 1024-bit discriminants)
// - `proof` is the serialized witness form (same size as `y`)
//
// On failure, returns `{NULL, 0}`.
ChiavdfByteArray chiavdf_prove_one_weso_fast(
    const uint8_t* challenge_hash,
    size_t challenge_size,
    const uint8_t* x_s,
    size_t x_s_size,
    size_t discriminant_size_bits,
    uint64_t num_iterations);

// Same as `chiavdf_prove_one_weso_fast`, but optionally invokes `progress_cb` from
// the proving thread every `progress_interval` iterations completed.
//
// If `progress_cb` is NULL or `progress_interval` is 0, no progress is reported.
ChiavdfByteArray chiavdf_prove_one_weso_fast_with_progress(
    const uint8_t* challenge_hash,
    size_t challenge_size,
    const uint8_t* x_s,
    size_t x_s_size,
    size_t discriminant_size_bits,
    uint64_t num_iterations,
    uint64_t progress_interval,
    ChiavdfProgressCallback progress_cb,
    void* progress_user_data);

// Computes a compact (witness_type=0) Wesolowski proof using the "streaming"
// bucket-accumulation algorithm (Trick 1), which requires the expected output
// `y_ref` up front (as used by bluebox compaction jobs).
//
// On success, returns `y || proof` (same format as `chiavdf_prove_one_weso_fast`).
ChiavdfByteArray chiavdf_prove_one_weso_fast_streaming(
    const uint8_t* challenge_hash,
    size_t challenge_size,
    const uint8_t* x_s,
    size_t x_s_size,
    const uint8_t* y_ref_s,
    size_t y_ref_s_size,
    size_t discriminant_size_bits,
    uint64_t num_iterations);

// Same as `chiavdf_prove_one_weso_fast_streaming`, but optionally invokes
// `progress_cb` from the proving thread every `progress_interval` iterations.
ChiavdfByteArray chiavdf_prove_one_weso_fast_streaming_with_progress(
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
    void* progress_user_data);

void chiavdf_free_byte_array(ChiavdfByteArray array);

#ifdef __cplusplus
}
#endif
