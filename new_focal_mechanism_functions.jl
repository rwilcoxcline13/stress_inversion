function OrderedEigenSystem(M)
    #
    #
    # compute the eigensystem with ordering such that smallest eigenvalue is in first vector, etc
    #
    #
    #
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
end

function moment_tensor(phi, delta, lamda)

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

    D, V = OrderedEigenSystem(M)

    P = V[:, 1]
    N = V[:, 2]
    T = V[:, 3]

    return P, N, T
  end

function PTN_2_fault_geometry(P, T, N)

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


function fault_parameters_2_nodal_planes(strike, dip, rake)

  phi = strike*pi/180
  delta = dip*pi/180
  lamda = rake*pi/180


  R1 = rotation_matrix([1, 0, 0], -pi)
  R2 = rotation_matrix([0, 0, -1], -pi)
  R = R2*R1
  Rinv = R'

  M = moment_tensor(phi, delta, lamda)
  M = Rinv*M*Rinv'

  P, N, T = moment_tensor_2_PTN(M)

  #Convert to seismo
  P = R*P
  N = R*N
  T = R*T

  if rake >= 90

    P = -P
    T = -T
  end

  Rcheck = rotation_matrix([0, 0, 1], -pi)

  #rotates eigenvectors by 180

  T = Rcheck*T
  P = Rcheck*P


  normal, normal_aux, slip, slip_aux = PTN_2_fault_geometry(P, T, N)


  normal_seismo = R*normal
  normal_aux_seismo = R*normal_aux
  slip_seismo = R*slip
  slip_aux_seismo = R*slip_aux

  x = [1, 0, 0]
  y = [0, 1, 0]
  z = [0, 0, 1]

  N = R*x
  E = R*y
  D = R*z

  strikecalc2 = -(atan2(normal_aux_seismo[2], normal_aux_seismo[1])-pi/2)
  dipcalc2 = (acos(normal_aux[3]))
  strike_direction2 = rotation_matrix(D, -strikecalc2)*N



  strikecalc1 = -(atan2(normal_seismo[2], normal_seismo[1])-pi/2)
  dipcalc1 = (acos(normal[3]))

  strikes = [strikecalc2, strikecalc1]
  dips = [dipcalc2, dipcalc1]
  rakes = zeros(2)

  for i = 1:2

    if dips[i] > pi/2
      dips[i] = pi-dips[i]
      strikes[i] = pi+strikes[i]
    end
  end



  slips = [slip_aux_seismo, slip_seismo]


  for i = 1:2
    Rf = fault_rotation(strikes[i], dips[i])
    azimuthal = Rf*N
    rakes[i] = pi-(acos(vecdot(azimuthal, slips[i]))/(norm(azimuthal)*norm(slips[i])))
  end


  if rake > 90

    rakes = pi - rakes
    strikes = 2*pi + strikes

    for i = 1:2

      if strikes[i] > 2*pi

        strikes[i] = strikes[i] - 2*pi
      end

    end


  elseif rake < 0

    for i = 1:2

      if strikes[i] > 2*pi

        strikes[i] = strikes[i] - 2*pi
        rakes[i] = -rakes[i]
        rakes[1:end .!=i] = rakes[1:end .!=i] - pi
      end

      if strikes[i] < 0

        strikes[i] = 2*pi+strikes[i]
        rakes[i] = rakes[i] - pi
        rakes[1:end .!=i] = -rakes[1:end .!=i]
      end

    end
  end



  NP1 = [strikes[1]*180/pi dips[1]*180/pi rakes[1]*180/pi]
  NP2 = [strikes[2]*180/pi dips[2]*180/pi rakes[2]*180/pi]

  return NP1, NP2

end
