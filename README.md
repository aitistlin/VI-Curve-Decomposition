# VI-Curve-Decomposition

This repository contains Python code and experimental data for analyzing voltage–current (V–I) curves in arc high-impedance fault (AHIF) scenarios. The core function is to calculate active and reactive power and to decompose the V–I trajectory into capacitive and inductive areas, which may assist in understanding the nature of the fault.

## Contents

- `ahif_vi_analysis.py`: Main analysis script for computing power and area features from experimental AHIF data.
- `demo_data/`: Contains example CSV files captured directly from an oscilloscope during AHIF simulations.

## Data Format

Each `.csv` file in `demo_data/` contains three columns:

- `"in s"` – Time in seconds  
- `"C1 in V"` – Voltage measurement (V)  
- `"C2 in A"` – Current measurement (A)  

> ⚠️ **Note:** Due to a setup oversight during oscilloscope scaling, the current values need to be multiplied by 3 to reflect the correct primary-side current.

## Usage

Run the script with your `.csv` file as input:

```bash
python ahif_vi_analysis.py demo_data/sample1.csv
