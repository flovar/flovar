import argparse
from typing import Tuple

import numpy as np

from utils import EXIT_CODE, THETA_COEF

class Params:
    data: str
    output: str
    freq: int
    depth: Tuple[float, float]
    mid: int
    total: int
    depth_vwc: int
    bnds: Tuple[Tuple[float, float], ...]
    inits: np.ndarray
    cw: float
    cv: float
    cs: float
    solid_frac: float
    cp: float
    dx: float
    dt: float
    theta_coef: Tuple[float, float, float]
    theta_coef_preset: int
    theta_coef_custom: Tuple[float, float, float]
    optim: str
    n_process: int
    lambda_ub: float
    lambda_lb: float
    flux_lb: float
    flux_ub: float
    loss_coef: Tuple[float, float, float]
    seed: int
    filter_outputs: bool
    def __init__(self, args: argparse.Namespace):
        super().__init__()
        args_dict = vars(args)
        for key in args_dict.keys():
            setattr(self, key, args_dict[key])


def get_args() -> Params:

    parser = argparse.ArgumentParser(
        description="VSF: Variablely Saturation Zone Flux Estimation from Temperature",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("data",
                        default=None,
                        type=str,
                        help="Path to csv data file")
    parser.add_argument("depth",
                        default=None,
                        nargs=2,
                        type=float,
                        help="Depth of sensor 1 and 2 relative to sensor 0")
    parser.add_argument("dt",
                        default=1800,
                        nargs="?",
                        type=float,
                        help="Time duration between two data records")
    parser.add_argument("--output",
                        default="out.npz",
                        type=str,
                        help="Path to output file")
    parser.add_argument("--freq",
                        default=0,
                        type=int,
                        help="Number of data records used in every stress period, can be omitted if using all data as a single stress period")
    parser.add_argument("--depth_vwc",
                        default=-1,
                        type=float,
                        help="Depth of VWC sensor relative to sensor 0, can be omitted if temperature and VWC sensors are at the same depth")
    parser.add_argument("--dx",
                        default=0.05,
                        type=float,
                        help="Spatial step size")
    parser.add_argument("--cw",
                        default=4.2e6,
                        type=float,
                        help="Heat capacity of water")
    parser.add_argument("--cv",
                        default=100,
                        type=float,
                        help="Heat capacity of vapor")
    parser.add_argument("--cs",
                        default=3.6e6,
                        type=float,
                        help="Heat capacity of solid")
    parser.add_argument("--solid_frac",
                        default=-1,
                        type=float,
                        help="Volumetric fraction of solid phase to calculate heat capacity of moist soil")
    parser.add_argument("--cp",
                        default=-1,
                        type=float,
                        help="Heat capacity of moist soil, can be omitted if using computed result of Cw, Cv and Cs")
    parser.add_argument("--theta_coef_preset",
                        default=0,
                        type=int,
                        choices=[0, 1, 2],
                        help="Using preset regression coefficients of Chung and Horton's empirical eauation (1987)")
    parser.add_argument("--theta_coef_custom",
                        default=None,
                        nargs=3,
                        type=int,
                        help="Using custom regression coefficients of Chung and Horton's empirical eauation (1987)")

    parser.add_argument("--optim",
                        default="bfgs",
                        type=str,
                        choices=["bfgs", "pso", "nelder-mead", "tnc", "slsqp", "powell"],
                        help="Optimization algorithm"
                        )
    parser.add_argument("--n_processes",
                        default=6,
                        type=int,
                        help="Number of processes used in parallel, only meaningful to L-BFGS-G and PSO algorithm")
    parser.add_argument("--lambda_lb",
                        default=0.5,
                        type=float,
                        help="Lower bound of thermal conductivity")
    parser.add_argument("--lambda_ub",
                        default=3,
                        type=float,
                        help="Upper bound of thermal conductivity")
    parser.add_argument("--flux_lb",
                        default=-2e-6,
                        type=float,
                        help="Lower bound of water flux")
    parser.add_argument("--flux_ub",
                        default=2e-6,
                        type=float,
                        help="Upper bound of water flux")
    parser.add_argument("--loss_coef",
                        default=[0.7, 0.2, 0.1],
                        nargs=3,
                        type=float,
                        help="Weighting coefficients of loss function")
    parser.add_argument("--seed",
                        default=42,
                        type=int,
                        help="Random seed")
    parser.add_argument("--filter_outputs",
                        default=False,
                        action="store_true",
                        help="Passing optimization outputs through a Butterworth filter to remove outliers")
    
    args = Params(parser.parse_args())
    return preprocess_args(args)

def preprocess_args(args: Params) -> Params:
    if args.solid_frac == -1 and args.cp == -1:
        print("[FATAL] Solid fraction or heat capacity of moist soil must be specified")
        exit(EXIT_CODE["HP_CP_MISSING"])
    args.mid       = round(args.depth[0] / args.dx)
    args.total     = round(args.depth[1] / args.dx)
    args.depth_vwc = round(args.depth_vwc / args.dx) if args.depth_vwc != -1 else args.mid
    args.bnds      = ((args.lambda_lb, args.lambda_ub),) * args.total + ((args.flux_lb, args.flux_ub), (args.flux_lb * 0.1, args.flux_ub * 0.1))
    args.inits     = np.full((args.total+2), args.lambda_lb)
    args.inits[-2] = (args.flux_lb + args.flux_ub) / 2
    args.inits[-1] = args.inits[-2] * 0.1

    if args.theta_coef_custom is not None:
        args.theta_coef = args.theta_coef_custom
    else:
        args.theta_coef = THETA_COEF[args.theta_coef_preset]
    return args
