# metcorr

Metallicity two-point correlation of galaxies (please cite Li et al. [2023](https://ui.adsabs.harvard.edu/abs/2023MNRAS.518..286L/abstract)). Deprojection code can also be found here.

Demo code:

    from metcorr import deproject, fluc_map, corr_func
    X, Y = deproject(met_map, cen_coord, PA, incl)
    x, y = X * pixel_size, Y * pixel_size
    met_fluc = fluc_map(x, y, metallicity, bins=np.arange(*rad_bin))
    corr_func(x, y, met_fluc, rad_bin=rad_bin)

Note that corr_func() is designed for read-world scales. Thus, the parameters x and y should be converted to kpc from pixel coordinates.

    deproject(shape, cen_coord=(0, 0), PA=0., incl=0., q0=0.):
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


    fluc_map(x, y, z, bins):
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


    corr_func(x, y, z, rad_bin, azi_bin=None, report=False):
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
            Following astronomical convention, PA = 0 mean the north
                (aligned to the y axis).
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
