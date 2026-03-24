from pathlib import Path
from typing import Tuple
from time import time
import random

from optimparallel import minimize_parallel
from scipy import optimize
from pyswarms.single import GlobalBestPSO
import pandas as pd
import numpy as np

from functions import upsample, forward, loss, butter_bandpass_filter
from settings import get_args, Params
from utils import EXIT_CODE


def read_data(data: str) -> np.ndarray:
    data_pth = Path(data)
    if not data_pth.exists() or data_pth.is_dir():
        print(f"[FATAL] Data file {data_pth.name} cannot be found")
        exit(EXIT_CODE["FILE_NOT_FOUND"])
    try:
        data_array = pd.read_csv(data_pth, header=None).to_numpy()
        return data_array
    except:
        print(f"[FATAL] Data file {data_pth.name} cannot be read")
        exit(EXIT_CODE["FILE_NOT_READABLE"])


def seed_everything(seed: int):
    random.seed(seed)
    np.random.seed(seed)


def fitness_fcn(x: np.ndarray,
                init_data: np.ndarray,
                raw_data: np.ndarray,
                delta_x: float,
                delta_time:float,
                m1: int,
                m_vwc:int,
                theta: float,
                cp: float,
                cw: float,
                cv: float,
                theta_coef: Tuple[float, float, float],
                alpha: float,
                beta: float,
                gamma: float) -> float:
    predict_y = forward(x, init_data, raw_data, delta_x, 
                        delta_time, m1, cp, cw, cv)
    return loss(predict_y[:-1], raw_data[:-1, 1], x[m_vwc], theta, theta_coef,
                x[:-2], alpha, beta, gamma)


def fitness_fcn_pso(xs: np.ndarray, **kwargs):
    losses = []
    for x in xs:
        losses.append(fitness_fcn(x, **kwargs))
    return losses


def pso_optimize(args: Params,
                 init_data: np.ndarray,
                 raw_data: np.ndarray,
                 theta: float,
                 cp: float) -> np.ndarray:
    options   = {'c1': 0.5, 'c2': 0.3, 'w':0.9}
    optimizer = GlobalBestPSO(n_particles=10,
                              dimensions=len(args.bnds), 
                              options=options, 
                              bounds=np.array(args.bnds).transpose(1, 0))
    kwargs    = {
        "init_data":   init_data,
        "raw_data":    raw_data,
        "delta_x":     args.dx,
        "delta_time":  args.dt,
        "m1":          args.mid,
        "m_vwc":       args.depth_vwc,
        "theta":       theta,
        "cp":          cp,
        "cw":          args.cw,
        "cv":          args.cv,
        "theata_coef": args.theta_coef,
        "alpha":       args.loss_coef[0],
        "beta":        args.loss_coef[1],
        "gamma":       args.loss_coef[2]
    }
    start            = time()
    min_loss, res    = optimizer.optimize(fitness_fcn_pso, 80, 
                                         n_processes=args.n_process,
                                         verbose=False, **kwargs)
    print(f"[DBG] Time used {(time() - start)*1000}ms, loss: {min_loss}")
    return res
            

def common_optimize(args: Params,
                    init_data: np.ndarray,
                    raw_data: np.ndarray,
                    theta: float,
                    cp: float,
                    optim: str,
                    inits: np.ndarray) -> np.ndarray:
    forward_args = (init_data,
                    raw_data,
                    args.dx,
                    args.dt,
                    args.mid,
                    args.depth_vwc,
                    theta,
                    cp,
                    args.cw,
                    args.cv,
                    args.theta_coef,
                    args.loss_coef[0],
                    args.loss_coef[1],
                    args.loss_coef[2],
    )
    if optim == "bfgs":
        start = time()
        res = minimize_parallel(fitness_fcn,
                                inits,
                                args=forward_args,
                                bounds=args.bnds)
        print(f"[DBG] Time used {(time() - start)*1000}ms, loss: {res.fun}")
        return res.x
    else:
        if optim == "nelder-mead":
            optim = "Nelder-Mead"
        if optim == "tnc":
            optim = "TNC"
        if optim == "slsqp":
            optim = "SLSQP"
        if optim == "powell":
            optim = "Powell"
        start = time()
        res = optimize.minimize(fitness_fcn,
                                inits,
                                args=forward_args,
                                bounds=args.bnds,
                                method=optim)
        print(f"[DBG] Time used {(time() - start)*1000}ms, loss: {res.fun}")
        return res.x

def optimize_single(data: np.ndarray, args: Params, inits: np.ndarray) -> np.ndarray:
    init_data = upsample(np.array([0, args.depth[0], args.depth[1]]), 
                         data[0, :-1],
                         np.arange(0, args.depth[1] + args.dx, args.dx))
    raw_data  = data[:, :-1]
    theta     = float(np.mean(data[:, -1]))
    cp        = args.cs * args.solid_frac + args.cw * theta + args.cv * \
                (1 - args.solid_frac - theta) if args.cp == -1 else args.cp
    if args.optim == "pso":
        return pso_optimize(args, init_data, raw_data, theta, cp)
    elif args.optim in ["bfgs", "nelder-mead", "tnc", "slsqp", "powell"]:
        return common_optimize(args,
                               init_data,
                               raw_data,
                               theta,
                               cp,
                               args.optim,
                               inits)
    else:
        print(f"[FATAL] Optimizer {args.optim} is not supported")
        exit(EXIT_CODE["OPTIMIZER_NOT_SUPPORTED"])

def main():
    args = get_args()
    data = read_data(args.data)
    seed_everything(args.seed)
    if args.freq > 0:
        n_optimization = len(data) // args.freq
        inits = args.inits
        optis = []
        for idx, i in enumerate(range(0, len(data), args.freq)):
            print(f"[INFO] Optimizing {idx+1} / {n_optimization}")
            res = optimize_single(data[i:i+args.freq], args, inits)
            optis.append(res)
        optis = np.array(optis)
        if args.filter_outputs:
            optis[:, -2] = butter_bandpass_filter(optis[:, -2], 1, 30, 300, 4)
            optis[:, -1] = butter_bandpass_filter(optis[:, -1], 1, 30, 300, 4)
    else:
        print(f"[INFO] Optimizing single period")
        optis = optimize_single(data, args, args.inits).reshape(1, -1)
    np.savez_compressed(args.output,
                        lambdas=optis[:, :-2],
                        flux=optis[:, -2],
                        vflux=optis[:, -1])
    print("[INFO] Optimization results:")
    print("[INFO] water flux (cm/h): ", end="")
    if args.freq > 0:
        for opti in optis:
            print(f"{opti[-2]*3600*100:.4f}", end=" ")
        print()
    else:
        print(f"{optis[0, -2]*3600*100:.4f}")
    print("[INFO] vapor flux (cm/h): ", end="")
    if args.freq > 0:
        for opti in optis:
            print(f"{opti[-1]*3600*100:.4f}", end=" ")
        print()
    else:
        print(f"{optis[0, -1]*3600*100:.4f}")
    print("[INFO] Done")


if __name__ == "__main__":
    main()