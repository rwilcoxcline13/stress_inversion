def rotationmat3D(rotation_angle, rotation_axis):

    norm = np.linalg.norm(rotation_axis)
    unitvector = rotation_axis/norm
    u = unitvector[0]
    v = unitvector[1]
    w = unitvector[2]
    c = np.cos(rotation_angle)
    s = np.sin(rotation_angle)
    R = np.zeros((3, 3))
    R[0, 0] = u**2 + (v**2+w**2)*c
    R[0, 1] = u*v*(1-c) - w*s
    R[0, 2] = u*w*(1-c) + v*s
    R[1, 0] = u*v*(1-c) + w*s
    R[1, 1] = v**2 + (u**2 + w**2)*c
    R[1, 2] = v*w*(1-c) - u*s
    R[2, 0] = u*w*(1-c) - v*s
    R[2, 1] = v*w*(1-c) + u*s
    R[2, 2] = w**2 + (u**2 + v**2)*c
    Rmat = np.asmatrix(R)

    return Rmat
        
