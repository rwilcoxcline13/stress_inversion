def cartesian2spherical(cartesian_coords):

    spherical_coords = np.zeros(cartesian_coords.shape)

    for i in range(len(cartesian_coords)):

        #radius

        spherical_coords[i,0] = np.sqrt(np.sum(cartesian_coords[i, :]**2))

        #angle 1

        spherical_coords[i, 1] = np.arccos(cartesian_coords[i, 0]/spherical_coords[i, 0])

        #angle 2

        spherical_coords[i, 2] = np.arccos(cartesian_coords[i, 1]/(np.sqrt(cartesian_coords[i, 3]**2+cartesian_coords[i, 2]**2+cartesian_coords[i, 1]**2)))

        #angle 3

        if cartesian_coords[i, 3] >= 0:

            spherical_coords[i, 3] = np.arccos(cartesian_coords[i, 2]/(np.sqrt(cartesian_coords[i,3]**2+cartesian_coords[i, 2]**2)))

        else:

            spherical_coords[i, 3] = 2*np.pi-np.arccos(cartesian_coords[i, 2]/(np.sqrt(cartesian_coords[i,3]**2+cartesian_coords[i, 2]**2)))


        return spherical_coords
