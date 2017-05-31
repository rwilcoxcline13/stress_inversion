function rotation_matrix(axis, rotation_angle)

  #=
  Creates a rotation matrix for a rotation around an arbitrary axis as
  described by the Euler-Rodrigues formula
  =#

    axis_normalized = axis/sqrt(vecdot(axis, axis))
    a = cos(rotation_angle/2)
    b = -axis_normalized[1]*sin(rotation_angle/2)
    c = -axis_normalized[2]*sin(rotation_angle/2)
    d = -axis_normalized[3]*sin(rotation_angle/2)
    rot_mat = [a^2+b^2-c^2-d^2 2*(b*c-a*d) 2*(b*d+a*c);
               2*(b*c+a*d) a^2+c^2-b^2-d^2 2*(c*d-a*b);
               2*(b*d-a*c) 2*(c*d+a*b) (a^2+d^2-b^2-c^2)]

    return rot_mat

function stress_rotation(phi_MCS, theta_MCS, rho = 0)

    #= phi_MCS is the trend of the most compressive stress measured
        off of east toward north (positive)

       theta_MCS is the plunge of the most compressive stress measured
        downward off the of the horizontal (positive)

       rho is a rotation about the MCS measured clockwise (positive)

       all angles are measured in radians from [-pi, pi]
    =#

    R_strike = rotation_matrix([0, 0, -1], phi_MCS)
    R_dip = rotation_matrix([1, 0, 0], theta_MCS)
    R_inter = dot(R_dip, R_strike)
    rho_axis = R_inter[:, 1]
    R_rho = rotation_matrix(-rho_axis, rho)
    R_stress = dot(R_rho, R_inter)

    return R_stress
