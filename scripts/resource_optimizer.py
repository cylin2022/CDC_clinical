#!/usr/bin/env python3
"""System Resource Optimizer for CPU, RAM, and I/O allocation."""

import os
import shutil
import threading
import time

_CACHE_LOCK = threading.Lock()
_CACHE_STATUS = {
    "warmed": False,
    "warming": False,
    "last_warm_sec": 0.0,
    "bytes_warmed": 0,
    "error": "",
}


def read_meminfo_gb():
    info = {}
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2:
                    key = parts[0].rstrip(":")
                    info[key] = int(parts[1]) / (1024.0 * 1024.0)  # KB -> GB
    except Exception:
        pass
    return {
        "total_gb": round(info.get("MemTotal", 0.0), 1),
        "available_gb": round(info.get("MemAvailable", 0.0), 1),
        "cached_gb": round(info.get("Cached", 0.0), 1),
    }


def read_loadavg():
    try:
        with open("/proc/loadavg") as f:
            p = f.read().split()
            return {"load_1m": float(p[0]), "load_5m": float(p[1]), "load_15m": float(p[2])}
    except Exception:
        return {"load_1m": 0.0, "load_5m": 0.0, "load_15m": 0.0}


def get_shm_avail_gb():
    try:
        st = shutil.disk_usage("/dev/shm")
        return round(st.free / (1024.0 ** 3), 1)
    except Exception:
        return 0.0


def get_system_resources():
    cpus = os.cpu_count() or 4
    mem = read_meminfo_gb()
    load = read_loadavg()
    shm_free_gb = get_shm_avail_gb()
    return {
        "cpu_cores": cpus,
        "mem_total_gb": mem["total_gb"],
        "mem_available_gb": mem["available_gb"],
        "mem_cached_gb": mem["cached_gb"],
        "shm_free_gb": shm_free_gb,
        **load,
        "cache_status": dict(_CACHE_STATUS),
    }


def compute_execution_plan(queued_count=1, active_count=0, requested_threads=0):
    """Compute optimal concurrent workers and threads per job from live CPU/RAM."""
    res = get_system_resources()
    cpus = res["cpu_cores"]
    avail_gb = res["mem_available_gb"]
    load_1m = res["load_1m"]

    # Reserve ~10% CPUs (min 2) for OS / web server / I/O
    reserve_cpus = max(2, cpus // 10)
    usable_cpus = max(1, cpus - reserve_cpus)

    # Each concurrent job uses ~25 GB working RAM (Kraken2 DB is shared via mmap page cache)
    # and ~12 GB /dev/shm when building temporary competitive minimap2 index.
    ram_budget_jobs = max(1, int(max(avail_gb - 48.0, 25.0) // 25.0))
    cpu_budget_jobs = max(1, usable_cpus // 8)
    shm_budget_jobs = max(1, int(max(res["shm_free_gb"] - 64.0, 12.0) // 12.0))
    max_workers = min(12, ram_budget_jobs, cpu_budget_jobs, shm_budget_jobs)

    # Scale threads based on current + queued demand
    demand = max(1, active_count + queued_count)
    target_workers = min(max_workers, demand)

    # Subtract approximate external CPU usage beyond our active jobs if host is busy
    effective_cpus = max(4, int(usable_cpus - max(0.0, load_1m - active_count * 12) * 0.5))
    auto_threads = max(8, min(64, effective_cpus // target_workers))

    # If user explicitly requested high threads and only 1 job is running, honor up to usable_cpus
    if requested_threads and requested_threads > 0 and target_workers == 1:
        threads_per_job = min(usable_cpus, max(requested_threads, auto_threads))
    else:
        threads_per_job = auto_threads

    use_mmap_kraken = avail_gb >= 140.0
    use_shm_tmp = res["shm_free_gb"] >= 32.0

    return {
        "max_workers": max_workers,
        "target_workers": target_workers,
        "threads_per_job": threads_per_job,
        "use_mmap_kraken": use_mmap_kraken,
        "use_shm_tmp": use_shm_tmp,
        "resources": res,
    }


def warm_kraken2_cache(db_dir, async_run=True):
    """Warm Kraken2 *.k2d files into OS page cache so --memory-mapping runs at RAM speed."""
    k2_dir = os.path.join(db_dir, "kraken2", "k2_pluspf_20261022")
    if not os.path.isdir(k2_dir):
        return

    def _warm():
        with _CACHE_LOCK:
            if _CACHE_STATUS["warming"] or _CACHE_STATUS["warmed"]:
                return
            _CACHE_STATUS["warming"] = True
        t0 = time.time()
        total_bytes = 0
        try:
            mem = read_meminfo_gb()
            if mem["available_gb"] < 150.0:
                _CACHE_STATUS["error"] = f"Skipped cache warm (avail {mem['available_gb']} GB < 150 GB)"
                return
            # If already mostly cached (>140 GB page cache), skip heavy read
            if mem["cached_gb"] >= 140.0:
                _CACHE_STATUS["warmed"] = True
                _CACHE_STATUS["last_warm_sec"] = round(time.time() - t0, 2)
                return
            for fn in ("taxo.k2d", "opts.k2d", "hash.k2d"):
                fp = os.path.join(k2_dir, fn)
                if not os.path.isfile(fp):
                    continue
                with open(fp, "rb") as f:
                    while True:
                        buf = f.read(16 * 1024 * 1024)
                        if not buf:
                            break
                        total_bytes += len(buf)
            _CACHE_STATUS["warmed"] = True
            _CACHE_STATUS["bytes_warmed"] = total_bytes
            _CACHE_STATUS["last_warm_sec"] = round(time.time() - t0, 2)
        except Exception as exc:
            _CACHE_STATUS["error"] = str(exc)
        finally:
            _CACHE_STATUS["warming"] = False

    if async_run:
        threading.Thread(target=_warm, daemon=True).start()
    else:
        _warm()

