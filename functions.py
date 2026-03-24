from typing import Tuple

from scipy.signal import butter, lfilter
from numpy import linalg as LA
from scipy import interpolate
import numpy as np


def butter_bandpass_filter(data: np.ndarray,
                           lowcut: float,
                           highcut: float,
                           fs: int,
                           order=5) -> np.ndarray:
    nyq  = 0.5 * fs
    low  = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    y    = lfilter(b, a, data)
    return np.array(y)

def upsample(depth: np.ndarray, temp: np.ndarray, space: np.ndarray) -> np.ndarray:
    spline = interpolate.UnivariateSpline(depth, temp, k=2)
    return np.array(spline(space))


def calc_theta(target_lamd: float, theta_coef: Tuple[float, float, float]) -> float:
    param_a       = theta_coef[1]**2
    param_b       = -(2*theta_coef[1]*target_lamd-2*theta_coef[0]*theta_coef[1]+theta_coef[2]**2)
    param_c       = (target_lamd-theta_coef[0])**2
    theta_dual    = np.roots([param_a, param_b, param_c])
    predict_theta = min(abs(np.real(theta_dual[0])), abs(np.real(theta_dual[1])))
    return predict_theta


def loss(predict_y: np.ndarray,
         y: np.ndarray,
         target_lamd: float,
         theta: float,
         theta_coef: Tuple[float, float, float],
         lamd=None,
         alpha=0.7,
         beta=0.2,
         gamma=0.1):
    temp_rmse     = np.sqrt(LA.norm(predict_y - y)**2/len(y))
    vwc_mae       = abs(calc_theta(target_lamd, theta_coef) - theta)
    # minimize distance between lamds
    lamd_residual = 0 if lamd is None else np.sqrt(LA.norm(lamd[1:] - lamd[:-1])**2/(len(lamd) - 1))
    loss          = alpha * temp_rmse + beta * 10 * vwc_mae + gamma * lamd_residual
    return loss


def forward(x: np.ndarray,
            init_data: np.ndarray,
            raw_data: np.ndarray,
            delta_x: float,
            delta_time: float,
            mid: int,
            cp: float,
            cw: float,
            cv: float):
    data      = init_data.copy()
    lamd      = x[:-2]
    qL        = x[-2]
    qV        = x[-1]
    llength   = len(lamd)
    dlength   = raw_data.shape[0]
    a         = np.zeros((llength - 1, llength - 1))
    b         = np.zeros((llength - 1,))
    predict_y = np.zeros((dlength, ))

    flat_a = a.ravel()
    flat_a[::llength]            = np.hstack((1+(lamd[2:]+2*lamd[1:-1]+lamd[:-2])*delta_time/(4*cp*delta_x**2), -1)) # type: ignore
    flat_a[1::llength]           = -(lamd[1:-1]+lamd[2:])*delta_time/(4*cp*delta_x**2)+(cw*qL+cv*qV)*delta_time/(4*delta_x*cp)  # type: ignore
    flat_a[llength - 1::llength] = np.hstack((-(lamd[2:-1]+lamd[1:-2])*delta_time/(4*cp*delta_x**2)-(cw*qL+cv*qV)*delta_time/(4*delta_x*cp), 1))  # type: ignore
    
    for i in range(dlength - 1):
        b[:-1]        = delta_time/(4*cp*delta_x**2)*((lamd[2:]+lamd[1:-1])*(data[2:-1] - data[1:-2])-(lamd[1:-1]+lamd[:-2])*(data[1:-2]-data[:-3]))+data[1:-2]-(cw*qL+cv*qV)*delta_time/(4*delta_x*cp)*(data[2:-1]-data[:-3])  # type: ignore
        b[0]         += ((lamd[1] + lamd[0]) * delta_time / (4 * cp * delta_x ** 2) + (cw * qL + cv * qV) * delta_time / (4 * delta_x * cp)) * raw_data[i+1, 0]
        temp_i        = LA.lstsq(a, b, rcond=None)[0].reshape((-1))
        data[1 : -1]  = temp_i
        data[0]       = raw_data[i + 1, 0]
        data[llength] = raw_data[i + 1, 2]
        predict_y[i]  = temp_i[mid - 1]
    
    predict_y[-1] = predict_y[-2]
    
    return predict_y
