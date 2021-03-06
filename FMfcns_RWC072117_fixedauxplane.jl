function OrderedEigenSystem(M)
    #
    #
    # compute the eigensystem with ordering such that smallest eigenvalue is in first vector, etc
    #
    #
    # writted by Larissa Lu (UM/UROP) summer 2016
    #

    # find the P,N, and T axes
    (D, V) = eig(M);

    #concatenate [1;2;3] onto D to make a matrix with 2 columns
    unsortedValues = hcat(D, [1; 2; 3]);
    #sort the eigenvalues in the matrix while reordering pos accordingly
    sortedValues = sortrows(unsortedValues);

    #store the sorted eigenvalues in a temporary vector
    temp = vec(sortedValues[:, 1]);
    #store the position values as integers
    spos = round(Int64, sortedValues[:, 2]);
    #turn sorted values into a diagonal matrix
    D = diagm(temp);

    #make a new empty temp matrix to store our reordered V
    tempV = zeros(3,3);
    #sort eigenvectors that correlate to our sorted eigenvalues
    for i = 1:3
        tempV[:, i] = V[:, spos[i]];
    end
    V = tempV;

    return D, V

end

function rotation_matrix(axis, rotation_angle)

    axis_normalized = axis/sqrt(vecdot(axis, axis))
    a = cos(rotation_angle/2)
    b = -axis_normalized[1]*sin(rotation_angle/2)
    c = -axis_normalized[2]*sin(rotation_angle/2)
    d = -axis_normalized[3]*sin(rotation_angle/2)
    rot_mat = [a^2+b^2-c^2-d^2 2*(b*c-a*d) 2*(b*d+a*c);
               2*(b*c+a*d) a^2+c^2-b^2-d^2 2*(c*d-a*b);
               2*(b*d-a*c) 2*(c*d+a*b) (a^2+d^2-b^2-c^2)]

    return rot_mat
end

function moment_tensor(phi, delta, lamda)

#Calculates moment tensor using Aki-Richards (seismo) convention

  Mxx = -(sin(delta)*cos(lamda)*sin(2*phi)+sin(2*delta)*sin(lamda)*sin(phi)^2)
  Mxy = (sin(delta)*cos(lamda)*cos(2*phi)+1/2*sin(2*delta)*sin(lamda)*sin(2*phi))
  Mxz = -(cos(delta)*cos(lamda)*cos(phi)+cos(2*delta)*sin(lamda)*sin(phi))
  Myy = (sin(delta)*cos(lamda)*sin(2*phi)-sin(2*delta)*sin(lamda)*cos(phi)^2)
  Myz = -(cos(delta)*cos(lamda)*sin(phi)-cos(2*delta)*sin(lamda)*cos(phi))
  Mzz = sin(2*delta)*sin(lamda)

  M = [Mxx Mxy Mxz; Mxy Myy Myz; Mxz Myz Mzz]

  return M
end

function moment_tensor_2_PTN(M)

#Calculates the P, T, and N axes

    D, V = OrderedEigenSystem(M)

    P = V[:, 1]
    N = V[:, 2]
    T = V[:, 3]

    return P, N, T
  end

function PTN_2_fault_geometry(P, T, N)

  #= Calculates the normal and slip vectors of fault and auxilary plane given
  P, T, N axes =#

        normal = 1/sqrt(2)*(T+P)
        slip = 1/sqrt(2)*(T-P)

        normal_aux = 1/sqrt(2)*(T-P)
        slip_aux = 1/sqrt(2)*(T+P)

        return normal, normal_aux, slip, slip_aux
end

function fault_rotation(phi, delta)

  N = [-1, 0, 0]
  E = [0, 1, 0]
  D = [0, 0, -1]

  R_strike = rotation_matrix(D, -phi)
  Np = R_strike*N
  Ep = R_strike*E
  Dp = R_strike*D

  R_dip = rotation_matrix(Np, -delta)
  Npp = R_dip*Np
  Epp = R_dip*Ep
  Dpp = R_dip*Dp
  Rf = R_dip*R_strike

  return Rf

end

function strike_dip(normal_vec)


  strike_calc = atan2(normal_vec[1], normal_vec[2])


  D = [0, 0, -1]
  dip_calc = acos(dot(D, normal_vec))

  if dip_calc > pi/2
    dip_calc = pi-dip_calc
    strike_calc = pi+strike_calc
  end

    strike_calc = 2*pi-strike_calc

    if strike_calc >=2*pi
        strike_calc = strike_calc - 2*pi
    end




  return [strike_calc dip_calc]
end

function calc_rake(normal, slip)

    if normal[3] < 0

        sl = -rotation_matrix([0, 0, 1], pi/2)*slip
        n = -rotation_matrix([0, 0, 1], -pi/2)*normal
    else
        sl = rotation_matrix([0, 0, 1], pi/2)*slip
        n = rotation_matrix([0, 0, 1], -pi/2)*normal
    end

    h = [-sl[2], sl[1], 0]

    #forced rounding to prevent rounding error divergence 08/24/17
    rake = acos(round(dot(h, n)/norm(h), 8))

    if sl[3] > 0
        rake = rake
    else
        rake = -rake
    end


    return rake
end

function fault_parameters_2_nodal_planes(strike, dip, rake)

    phi = strike*pi/180
    delta = dip*pi/180
    lamda = rake*pi/180

    M = moment_tensor(phi, delta, lamda)

    #Calculate P, N, and T axes in standard convention and rotate to seismo

    P, N, T = moment_tensor_2_PTN(M)



    #Calculate and Rotate normal and slip vectors of fault and auxilary plane.

    normal, normal_aux, slip, slip_aux = PTN_2_fault_geometry(P, T, N)

    #Calculate strike and dip

    np1_sd = strike_dip(normal)
    np2_sd = strike_dip(normal_aux)



    calcrake1 = calc_rake(normal_aux, slip_aux)
    calcrake2 = calc_rake(normal, slip)


    np1 = [np1_sd[1], np1_sd[2], calcrake1]
    np2 = [np2_sd[1], np2_sd[2], calcrake2]


    return [np1, np2]

end
