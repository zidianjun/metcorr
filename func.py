
from ctypes import *
import numpy as np
from scipy.stats import binned_statistic
import os
import platform

suffix = 'linux' if platform.platform().split('-')[0] == 'Linux' else 'macos'
abs_path = os.path.dirname(os.path.abspath(__file__))
lib = CDLL(abs_path + '/lib_' + suffix + '.so')



def _bin_array(r, f, bins): # phase=0
    stat = binned_statistic(r, f, bins=bins)
    mask = ~np.isnan(stat.statistic)
    bin_r = stat.bin_edges[:-1][mask]
    bin_f = stat.statistic[mask]
    return bin_r, bin_f

def _step(rad, met, bin_rad, bin_met):
    rad_matrix = abs(np.subtract.outer(rad, bin_rad))
    min_value = np.expand_dims(np.min(rad_matrix, axis=1), axis=1)
    # np.where will return ALL the min values!
    x, y = np.where(rad_matrix == min_value)
    # If having recurring items, remove them.
    y = y[np.insert(np.diff(x) != 0, 0, True)]
    step_func = bin_met[y]
    fluc = met - step_func
    return fluc

def _bin_stat(f, x, y, min_sep, max_sep, bin_size, report=False):
    three_arrays = np.concatenate([f, x, y], axis=0)
    lib.group_by.argtypes = [POINTER(c_float)] + [c_int] * 3 + [c_float] * 3
    lib.group_by.restype = POINTER(c_float)
    c_array = np.ctypeslib.as_ctypes(three_arrays.astype(np.float32))
    length = int((max_sep - min_sep) / bin_size)
    c_res = lib.group_by(c_array, len(f), length, int(report),
                         min_sep, max_sep, bin_size)
    py_res = cast(c_res, POINTER(c_float * (length * 2))).contents
    res = np.array(list(py_res), dtype=float).reshape(2, -1)
    return res[1]

def _whole_stat(f, x, y, min_sep, max_sep, bin_size, min_pa, max_pa, azi_size,
                report=False):
    three_arrays = np.concatenate([f, x, y], axis=0)
    lib.group.argtypes = [POINTER(c_float)] + [c_int] * 4 + [c_float] * 6
    lib.group.restype = POINTER(c_float)
    c_array = np.ctypeslib.as_ctypes(three_arrays.astype(np.float32))
    len1 = int((max_sep - min_sep) / bin_size)
    len2 = int((max_pa - min_pa) / azi_size)
    c_res = lib.group(c_array, len(f), len1, len2, int(report),
                      min_sep, max_sep, bin_size, min_pa, max_pa, azi_size)
    py_res = cast(c_res, POINTER(c_float * (len1 * len2 * 2))).contents
    res = np.array(list(py_res), dtype=float).reshape(2, -1)
    return res[1].reshape((len2, len1))

def _tpcf(f, x, y, rad_bin, azi_bin=None, report=False):
    mean2, sigma2 = np.mean(f) ** 2, np.std(f) ** 2  # mean is 0
    min_sep, max_sep, bin_size = rad_bin
    max_sep = np.arange(*rad_bin)[-1] + bin_size
    if azi_bin is None:
        bin_scorr = _bin_stat(f, x, y, report=report,
                              min_sep=min_sep, max_sep=max_sep, bin_size=bin_size)
        return (bin_scorr - mean2) / sigma2
    else:
        min_pa, max_pa, azi_size = azi_bin
        max_pa = np.arange(*azi_bin)[-1] + azi_size
        bin_scorr = _whole_stat(f, x, y, report=report,
                                min_sep=min_sep, max_sep=max_sep, bin_size=bin_size,
                                min_pa=min_pa, max_pa=max_pa, azi_size=azi_size)
        res = (bin_scorr - mean2) / sigma2
        res[:, 0] = 1.
        return res

def _incl_to_cosi(incl, q0=0.):
    cos = np.cos(incl * np.pi / 180)
    return np.sqrt((cos ** 2 - q0 ** 2) / (1 - q0 ** 2)) if cos > q0 else 0.




def inv_where(val_arr, bool_arr, padding=np.nan):
    """
    Recover the val_arr 
    Inverse function of np.where().
    Calculate x = np.where(bool_arr, val_arr, padding),
        when val_arr and bool_arr are known.
    Parameters:
        val_arr: 1D np.array
            The array that has valid values.

        bool_arr: 1D np.array
            The array that has boolean values showing where
                val_arr is valid

        padding: float.
            The value to be padded into the places where
                bool_arr is False
            Defalted to be np.nan.
    
    returns:
        1D np.array
        Padding val_arr using padding following the indication of bool_arr.
    """
    if int(np.sum(bool_arr)) != len(val_arr):
        raise ValueError("The number of 'True' in bool_arr should be equal to" +
                         "the length of val_arr!")
    res = np.ones(len(bool_arr)) * padding
    flag = 0
    for i in range(len(bool_arr)):
        if bool_arr[i]:
            res[i] = val_arr[flag]
            flag += 1
    return res


def deproject(shape, cen_coord, PA, incl, q0=0.):
    """
    Deproject the galaxy coordinates using rotation matrix.
    Parameters:
        shape: tuple
            The shape of the original metallicity map from an IFU.

        cen_coord: 2-element tuple
            The coordinates of the galaxy center, shaped as (center_x, center_y)

        PA: float (in unit of degree)
            The position angle. PA = 0 means that the semi long axis is
                aligned to the north.

        q0: float
            A factor related with the intrinsic galaxy disk thickness.
            Defaulted to be 0, meaning that the disk is infinitely thin.
    
    returns:
        A tuple of (X, Y)
    """
    height, width = shape
    cx, cy = cen_coord
    cosi = _incl_to_cosi(incl, q0=q0)
    theta = (PA + 90) * np.pi / 180  # Aligned to x axis
    dep_mat = np.array([[np.cos(theta), np.sin(theta)],
                        [-np.sin(theta) / cosi, np.cos(theta) / cosi]])
   
    x0, y0 = np.meshgrid(range(width), range(height))
    x0, y0 = x0.reshape(-1), y0.reshape(-1)
    xy_mat = np.stack([x0 - cx, y0 - cy], axis=0)
    X, Y = np.dot(dep_mat, xy_mat)
    return X, Y


def fluc_map(x, y, z, bins):
    """
    Compute the (metallicity) fluctuation map after removing the radial profile.
    Parameters:
        x, y, and z: 1D np.array
            x and y are the coordinates and z is the map (generally metallicities).
            Their sizes must be identical.

        bins: 1D np.array
            Radial bins
    
    returns:
        A 1D np.array containing fluctuations with the same size as z has.
    """
    if x.size != y.size or y.size != z.size:
        raise ValueError("The sizes of the three arrays must be equal!")
    x_arr, y_arr, z_arr = x.reshape(-1), y.reshape(-1), z.reshape(-1)
    good = ~np.isnan(x_arr) & ~np.isnan(y_arr) & ~np.isnan(z_arr)
    x_arr, y_arr, z_arr = x_arr[good], y_arr[good], z_arr[good]
    r_arr = np.sqrt(x_arr ** 2 + y_arr ** 2)
    bin_r, bin_z = _bin_array(r_arr, z_arr, bins=bins)
    return _step(r_arr, z_arr, bin_r, bin_z)


def corr_func(x, y, z, rad_bin, azi_bin=None, report=False):
    """
    Compute the two-point correlation function of a deprojected galaxy.
    Parameters:
        x, y, and z: 1D or 2D array
            x and y are the coordinates and z is generally the metallicity fluctuations.
            Their sizes must be identical.

        rad_bin: list
            The bins for separations, in the same unit as x_arr.
            Should be as [min separation, max separation, radial bin width]

        azi_bin: 3-element list or int, in the unit of degree.
            The bins for azimuthal expansion.
            Defaulted to be None.
            Following astronomical convention, PA = 0 means aligned to the x axis 
                (the same as mathematical convention).
            Should be as [min pitch angle, max pitch angle, azimuthal bin width]
            If float, indicate the bin width only, as [0, 180, azimuthal bin width]

        report: bool, optional
            If True, then print the procedures.
    
    returns:
        If azi_bin is None, return a 1D array, having the same shape as
            np.arange(min separation, max separation, radial bin width)
        Else return a 2D array, having the shape of 
            np.arange(min separation, max separation, radial bin width) and
            np.arange(min pitch angle, max pitch angle, azimuthal bin width)
    """
    if x.size != y.size or y.size != z.size:
        raise ValueError("The sizes of the three arrays must be equal!")
    x_arr, y_arr, z_arr = x.reshape(-1), y.reshape(-1), z.reshape(-1)
    good = ~np.isnan(x_arr) & ~np.isnan(y_arr) & ~np.isnan(z_arr)
    x_arr, y_arr, z_arr = x_arr[good], y_arr[good], z_arr[good]
    azi_bins = [0., 180., 180. / azi_bin] if type(azi_bin) == int else azi_bin
    return _tpcf(z_arr, x_arr, y_arr, rad_bin=rad_bin, azi_bin=azi_bins, report=report)




