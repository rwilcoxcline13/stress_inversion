def spherical2cartesian(spherical_coords):

    if spherical_coords.shape == (3,):

        spherical_coords = spherical_coords.reshape((3,1))

    cartesian_coords = np.zeros(spherical_coords.shape)

    for i in range(len(cartesian_coords)):

        r = spherical_coords[i, 0]
        angle1 = spherical_coords[i, 1]
        angle2 = spherical_coords[i, 2]
        angle3 = spherical_coords[i, 3]

        cartesian_coords[i, 0] = r*np.cos(angle1)
        cartesian_coords[i, 1] = r*np.sin(angle1)*np.cos(angle2)
        cartesian_coords[i, 2] = r*np.sin(angle1)*np.sin(angle2)*np.cos(angle3)
        cartesian_coords[i, 3] = r*np.sin(angle1)*np.sin(angle2)*np.sin(angle3)

    return cartesian_coords
