"""
Script to run the PML and ROI Transition Factorial Experiment on MYO001.
Configurations:
- Config B: PML 35 cells, sharp ROI boundary (_pml35_sharp)
- Config C: PML 15 cells, tapered ROI boundary (_pml15_taper)
- Config D: PML 35 cells, tapered ROI boundary (_pml35_taper)

Config A (PML 15, sharp ROI) was already completed during the verification step (_lu_verify).
"""
import subprocess
import sys
import time
import os

configs = [
    {
        "name": "Config B (PML 35, sharp ROI)",
        "suffix": "_pml35_sharp",
        "npml": 35,
        "roi_taper": False,
        "taper_cells": 5
    },
    {
        "name": "Config C (PML 15, tapered ROI)",
        "suffix": "_pml15_taper",
        "npml": 15,
        "roi_taper": True,
        "taper_cells": 10
    },
    {
        "name": "Config D (PML 35, tapered ROI)",
        "suffix": "_pml35_taper",
        "npml": 35,
        "roi_taper": True,
        "taper_cells": 10
    }
]

dataset_path = "data_cubdl_prueba_1/MYO001.hdf5"

if not os.path.exists(dataset_path):
    print(f"Error: Dataset {dataset_path} not found.")
    sys.exit(1)

print("=" * 60)
print("STARTING PML & ROI TAPER FACTORIAL EXPERIMENT ON MYO001")
print("=" * 60)

total_start = time.time()

for config in configs:
    print(f"\n>>> Running: {config['name']}")
    print(f"Parameters: npml={config['npml']}, roi_taper={config['roi_taper']}, taper_cells={config['taper_cells']}")
    
    cmd = [
        "python", "-u", "main_5MHz_fwi_real_data.py",
        dataset_path, "42", config["suffix"],
        "--npml", str(config["npml"]),
        "--taper-cells", str(config["taper_cells"])
    ]
    if config["roi_taper"]:
        cmd.append("--roi-taper")
        
    start_time = time.time()
    
    # We pipe stdout and stderr to a log file so we have records of each run
    log_filename = f"articulo/fwi_run{config['suffix']}.log"
    print(f"Logging to {log_filename}...")
    
    with open(log_filename, "w", encoding="utf-8") as log_file:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1
        )
        
        # Read output in real-time and print to console while writing to log file
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log_file.write(line)
            log_file.flush()
            
        process.wait()
        
    elapsed = time.time() - start_time
    print(f"\n>>> Finished: {config['name']} in {elapsed:.1f}s ({elapsed/60:.1f} min) with return code {process.returncode}")
    print("-" * 60)

total_elapsed = time.time() - total_start
print(f"\n==================================================")
print(f"FACTORIAL EXPERIMENT COMPLETED")
print(f"Total time: {total_elapsed/60:.1f} minutes")
print(f"==================================================")
