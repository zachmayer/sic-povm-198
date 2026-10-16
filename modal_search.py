"""Run search.py on many cloud cores with Modal (https://modal.com). Needs a Modal account.

  uvx modal run modal_search.py --d 198 --containers 25 --minutes 60

Each container runs one single-threaded search per core, each with its own seed and stopping at
its first SIC, which is written here as fiducial_<d>_seed<seed>_trial<trial>.txt. The search that
found fiducial_198.txt was this size: 25 containers of 32 cores for an hour.
"""

import multiprocessing
import os
from pathlib import Path

import modal

CORES = 32

image = (
    modal.Image.debian_slim(python_version="3.12").uv_pip_install("numpy", "scipy").add_local_python_source("search")
)
app = modal.App("sic-search", image=image)

with image.imports():
    import search


@app.function(cpu=CORES, memory=512 * CORES, timeout=24 * 3600)
def run(d: int, seeds: list[int], minutes: float) -> list[tuple[int, int, tuple[int, str] | None]]:
    # Modal sets OPENBLAS_NUM_THREADS to the core count. Each search needs one thread, and 32
    # searches with 32 threads each thrash; spawned workers load numpy with this setting instead.
    os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    with multiprocessing.get_context("spawn").Pool(len(seeds)) as pool:
        results = pool.starmap(search.search, [(d, seed, minutes) for seed in seeds])
    return [(seed, count, found) for seed, (count, found) in zip(seeds, results, strict=True)]


@app.local_entrypoint()
def main(d: int = 198, containers: int = 25, minutes: float = 60, seed: int = 0) -> None:
    seeds = [list(range(seed + c * CORES, seed + (c + 1) * CORES)) for c in range(containers)]
    trials = sics = 0
    for results in run.map([d] * containers, seeds, [minutes] * containers):
        for s, count, found in results:
            trials += count
            if found:
                k, text = found
                path = Path(f"fiducial_{d}_seed{s}_trial{k}.txt")
                path.write_text(text)
                sics += 1
                print(f"SIC found: {path}", flush=True)
    print(f"d = {d}: {trials} trials, {sics} SICs")
