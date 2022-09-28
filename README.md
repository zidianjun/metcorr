# metcorr
metallicity two-point correlation


from metcorr import deproject, corr_func

X, Y = deproject(met_map, cen_coord, PA, b2a)

x, y = X * pixel_size, Y * pixel_size

corr_func(x, y, metallicity)



help(deproject)

deproject(image, cen_coord=(0, 0), PA=0., b2a=1., q0=0.):
    """
    Deproject the galaxy coordinates using rotation matrix.
    Parameters:
        image: 2D array
            In general the original metallicity map from an IFU.

        cen_coord: 2-element tuple
            The coordinates of the galaxy center, shaped as (center_x, center_y)

        PA: float (in unit of degree)
            The position angle. PA = 0 means that the semi long axis is aligned to x axis.

        q0: float
            A factor related with the intrinsic galaxy disk thickness.
            q0 = 0 means that the disk is infinitely thin.
    
    returns:
        A tuple of (X, Y)
    """


help(corr_func)

corr_func(x, y, met, bin_size=.2, max_sep=5., report=False, adp=False):
              
    """
    Parameters:
        x, y, and met: 1D array.
            x and y are the coordinates and met is the metallicities.
        bin_size: float (in unit of kpc)
            It could be as small as the physical spatial resolution.
        max_sep: float (in unit of kpc)
            It is a maximum separation. It is for the sake of saving time
            since a two-point correlation is very close to zero and
            has no information at very large separation.
        report: bool, optional
            If True, then print the procedures and
            how long the two-point correlation function takes.
        adp: bool, optional
            If True, then removing the radial metallicity gradient will be
            processed in adaptive bins.
    
    returns:
        sep and ksi: 1D array
        sep is the separation distance, simply an array as
            [0, 1*bin_size, 2*bin_size, ...].
        ksi is the two-point correlation, as
            [1, 0.X, 0.Y, ...].
    """
