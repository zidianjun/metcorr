# metcorr
metallicity two-point correlation


    from metcorr import deproject, corr_func

    X, Y = deproject(met_map, cen_coord, PA, b2a)

    x, y = X * pixel_size, Y * pixel_size

    corr_func(x, y, metallicity)

