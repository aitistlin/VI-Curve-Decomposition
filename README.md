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

During the use of the Python script, the following parameters need to be modified in line 283 of the code:
num: Indicates the starting index for data analysis.

N: Specifies the number of fundamental frequency cycles to be analyzed.

base_fs: Represents the fundamental frequency (in Hz), typically 50 Hz for power systems.

filename: The name of the CSV file containing the data.

file_path: Automatically constructs the full path to the data file based on the current script location and the demo_data folder.

These parameters allow you to control the data source and the segment of data used for analysis.

```bash
  num = 32850
  N = 400
  filename = '131.CSV'
  base_fs = 50
  file_path = os.path.join(os.path.dirname(__file__), 'demo_data', filename)
  data = pd.read_csv(file_path)
