# VSF: Variably Saturated Zone Water Flux Estimation Based on Temperature Data

## Introduction

VSF is a Python package for estimating water and vapor fluxes in variably saturated zone based on temperature data. Detailed information can be found in the following paper:

> Meng Su, Rui Ma, Yunquan Wang, Ziyong Sun, Maosheng Yin, Mengyan Ge. Tracing the soil water flux in variably saturated zones using temperature data. (Under review in Water Resources Research)

## Requirements

Recent versions of the following packages for Python 3.8+ are required:

- NumPy 1.20+
- SciPy 1.10+
- Pandas 1.2+
- optimparallel 0.1.2
- pyswarms 1.3.0

## Usage

### Quickstart

1. 	Prepare the csv data file. Four data series are required for VSF. The first through third columns are temperatures monitored by the sensor (numbered 1) at the depth of interest and two sensors (numbered 0 and 2) above and below the depth of interest, the fourth column is measured volumetric water content (VWC, $\theta$) which depth can be specified with `--depth_vwc` option.
2. Run the following command in terminal:

```bash
python3 main.py <path to csv data file> <depth of sensor 1 and 2 relative to sensor 0 in meters> <time duration between two data records in seconds> --cp <heat capacity of moist soil>
```

3. The output file is a compressed `.npz` archive, which contains the following optimized parameters:
- `lambdas`: $\lambda$ at each discretized nodes, in $Js^{-1}m^{-1} {}^\circ C^{-1}$
- `flux`: Liquid water flux, in $m/s$
- `vflux`: Vapor flux, in $m/s$

### Command Line Interface Options

| Option | Description | Default | Unit |
| --- | --- | --- | --- |
| `data` | Path to csv data file | None | -- |
| `depth` | Depth of sensor 1 and 2 relative to sensor 0 | None | $m$ |
| `dt` | Time interval between two adjacent records | 1800 | $s$ |
| `--output` | Path to output file | out.npz | -- |
| `--freq` | Number of data records used in every stress period, can be omitted if using all data as a single stress period | 0 | -- |
| `--depth_vwc` | Depth of VWC sensor relative to sensor 0, can be omitted if temperature and VWC sensors are at the same depth | -1 | $m$ |
| `--dx` | Spatial step size | 0.05 | $m$ |
| `--cw` | Heat capacity of water | 4.2e6 | $Jm^{-3} {}^\circ C^{-1}$ |
| `--cv` | Heat capacity of vapor | 100 | $Jm^{-3} {}^\circ C^{-1}$ |
| `--cs` | Heat capacity of solid | 3.6e6 | $Jm^{-3} {}^\circ C^{-1}$ |
| `--solid_frac` | Volumetric fraction of solid phase to calculate heat capacity of moist soil | -1 | -- |
| `--cp` | Heat capacity of moist soil, can be omitted if using computed result of Cw, Cv and Cs | -1 | $Js^{-1}m^{-3} {}^\circ C^{-1}$ |
| `--theta_coef_preset` | Using preset regression coefficients of Chung and Horton's empirical eauation (1987) (Available: 0~2) | 0 | -- |
| `--theta_coef_custom` | Using custom regression coefficients of Chung and Horton's empirical eauation (1987) | None | -- |
| `--optim` | Optimization algorithm (Available: "bfgs", "pso", "nelder-mead", "tnc", "slsqp", "powell") | bfgs | -- |
| `--n_processes` | Number of processes used in parallel, only meaningful to L-BFGS-G and PSO algorithm | 6 | -- |
| `--lambda_lb` | Lower bound of thermal conductivity | 0.5 | $Js^{-1}m^{-1} {}^\circ C^{-1}$ |
| `--lambda_ub` | Upper bound of thermal conductivity | 3 | $Js^{-1}m^{-1} {}^\circ C^{-1}$ |
| `--flux_lb` |  Lower bound of water flux | -2e-06 | $m/s$ |
| `--flux_ub` | Upper bound of water flux | 2e-06 | $m/s$ |
| `--loss_coef` | Weighting coefficients of loss function | [0.7, 0.2, 0.1] | -- |
| `--seed` | Random seed | 42 | -- |
| `--filter_outputs` | Passing optimization outputs through a Butterworth filter to remove outliers | False | -- |

> Preset regression coefficients of Chung and Horton's empirical eauation (1987)

| No. | `b1` | `b2` | `b3` | Notes |
| --- | --- | --- | --- | --- |
| 0 | -0.197 | -0.962 | 2.521 | Clay |
| 1 | 0.243 | 0.393 | 1.534 | Loam |
| 2 | 0.228 | -2.406 | 4.909 | Sand |

### Examples

In this section, three examples are presented to help users be familiar with VSF package.

Examples 1, 2 and 3 shared the same following conditions in the data csv file named as `data.csv`: the distance of sensor 1 (at the depth of interest) and sensor 2 (below the depth of interest) relative to sensor 0 (above the depth of interest) are $0.25$ and $0.5m$, respectively; and the time interval between two adjacent data records is 1800s.

> Example 1: Using all data to form a single stress period with predefined heat capacity of moist soil ($C_p$, `-–cp` in options) $2.5×10^6 Jm^{-3} ℃^{-1}$. Note that if the length of the data records is very large, the data records length could be divided into several stress periods. Parameter estimation would be conducted for each stress period, and the estimated results from the former stress period would be used as the initial parameter values for the next stress period. Example 1 is the simplest application of VSF package with the following command:

```bash
python3 main.py data.csv 0.25 0.5 1800 --cp 2.5e6
```

> Example 2: Considering that the data records period may be very long that experiences different flow conditions in some applications, the whole period can be divided into several stress periods. Here, each stress period has a length of 1 day, which covers $48$ data records (`--freq 48`); the total simulation time can cover arbitrarily days and will be automatically divided into several stress periods by VSF package. In this example, we illustrated the alternative method for calculating the heat capacity of moist soil ($C_p=C_n θ_n+C_w θ + C_v (1-θ_n-θ)$), based on the default heat capacity of solid ($C_n$), water ($C_w$) and vapor ($C_v$), and customized ($0.55$) volumetric fraction of solid phase ($θ_n$, `--solid_frac` in options), instead of directly assigning the value. Then the command for Example 2 is as follows,

```bash
python3 main.py data.csv 0.25 0.5 1800 --freq 48 --solid_frac 0.55
```

> Example 3: Compared with other two examples, Example 3 provide more customized options.
> * Setting spatial step size (`--dx`) to $0.01m$, heat capacity of water (`--cw`) to $3.6 \times 10^6~Jm^{-3} {}^\circ C^{-1}$, heat capacity of vapor (`--cv`) to $90~Jm^{-3} {}^\circ C^{-1}$, heat capacity of moist soil to $2.8 \times 10^6~Jm^{-3} {}^\circ C^{-1}$, regression coefficients of Chung and Horton's empirical eauation (1987) `b1` to $-0.24$, `b2` to $-0.4$, `b3` to $1.5$.
> * Using truncated Newton algorithm (TNC) as optimization method (`--optim tnc`), adjusting weighting coefficients of loss function to $0.9, 0.05, 0.05$.
> * Setting `random` module and NumPy random seed to 9000, using Butterworth filter to remove optimized flux outliers

```bash
python3 main.py data.csv 0.25 0.5 1800 --freq 48 --dx 0.01 --cw 3.6e6 --cv 90 --cp 2.8e6 --theta_coef_custom -0.24 0.4 1.5 --optim tnc --loss_coef 0.9 0.05 0.05 --seed 9000 --filter_outputs
```

