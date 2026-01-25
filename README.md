# chiavdf (WesoForge fork)

This repository is a fork of the official Chia Network `chiavdf` implementation:
https://github.com/Chia-Network/chiavdf

`chiavdf` is Chia’s class-group VDF engine (evaluation + Wesolowski proving) used by timelords and full nodes.

This fork exists to support **Bluebox proof compaction** at scale by providing **lower-level Rust bindings** than the upstream `rust_bindings/` layer, with a focus on running many independent workers efficiently (multi-process / multi-core).

WesoForge (our compaction client, consuming these bindings):
https://github.com/Ealrann/WesoForge

This README is intentionally high-level and assumes you already know chiavdf’s primitives and the one‑Wesolowski (“compact witness”) workflow.

## Streaming One‑Wesolowski (known `y_ref`) — ~3× lower memory

For Bluebox compaction, each job already includes the expected output `y_ref` (the `VDFInfo.output` from the block). That means we can compute the Wesolowski prime `B = GetB(D, x0, y_ref)` **before** starting the squaring loop.

With `B` known up front, we can switch from the upstream “store intermediates then scan them” approach to a **single-pass streaming prover**:

- During squaring, whenever we hit a checkpoint `f(i·k·l)`, we immediately multiply it into the correct `ys[j][b]` buckets.
- We never materialize the full `O(ceil(T/(k·l)))` array of intermediate forms.
- Finalization (“folding” buckets into the proof form) is unchanged.

In practice, this reduces the memory footprint significantly (we typically see ~3× less RAM for the compaction workload), which matters when you run many workers in parallel.

## GetBlock optimization — incremental mapping, no lookup table

Streaming shifts more work into the “checkpoint update” path, where we repeatedly compute:

- `b = GetBlock(p, k, T, B)` to map each `(p, k, T, B)` to its bucket index.

The naïve implementation relies on modular exponentiation per call, and a straightforward optimization is to precompute a `GetBlock` lookup table — but that costs additional RAM (and allocation/initialization time).

This fork implements an **in-flight / incremental** `GetBlock` mapping:

- We keep a small rolling state derived from `B`, `k`, and `T`.
- When calls follow the expected sequential pattern (`p = 0, 1, 2, …`), each next `GetBlock(p)` becomes a handful of big‑integer ops (multiply/mod/div) rather than a full exponentiation.
- No large table is allocated.

The mapping is mathematically identical; this is purely an implementation-level optimization.

## (k, l) tuner — memory-budgeted parameter selection

The one‑Wesolowski proof parameters `(k, l)` trade off:

- checkpoint frequency (`~T/(k·l)` checkpoints),
- bucket state size (`~l·2^k` forms),
- and finalization work (`~l·2^k` NUCOMP operations + folding).

Upstream’s `ApproximateParameters()` is a solid heuristic, but for compaction we want an explicit knob to bias choices based on how much memory a worker is allowed to use.

This fork adds a lightweight **(k, l) tuner** that searches a small grid of candidates under a configurable memory budget and picks parameters that should be more throughput‑friendly.

Today, the gains are modest (and sometimes negligible), but it’s a useful foundation for future optimizations and for workloads where `(k, l)` becomes a real bottleneck.
