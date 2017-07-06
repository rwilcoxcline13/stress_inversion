

#=
function OrderedEigenSystem(M)
    # KEEP ORDERED EIGENSYSTEM FOR CONSISTENCY

    # More efficient sorting of eigevectors RWC 061417
    # compute the eigensystem with ordering such that smallest eigenvalue is in first vector, etc
    #
    #
    #
    #

    # find the P,N, and T axes
    (D, V) = eig(M);

    sorted_idx = sortperm(D)
    D_sorted = D[sorted_idx]
    V_sorted = V[:, sorted_idx]
    return D, V

end
=#

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

  #Rotation matrix for a fault plane in seismo convention

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

function strike_dip(normal_vec)|

  strike_calc = atan2(normal_vec[2], normal_vec[1])

  strike_calc = strike_calc - pi/2

  if strike_calc >= 2*pi
    strike_calc = strike_calc - 2*pi
  end

  if strike_calc < 0
    strike_calc = strike_calc + 2*pi
  end

  D = [0, 0, -1]
  dip_calc = acos(dot(D, normal_vec))

  if dip_calc > pi/2
    dip_calc = pi-dip_calc
    strike_calc = pi+strike_calc
  end


  return [strike_calc dip_calc]
end

function fault_parameters_2_nodal_planes(strike, dip, rake)

  phi = strike*pi/180
  delta = dip*pi/180
  lamda = rake*pi/180

  #Define rotation matrices to transform between seismo and standard conventions
  #R rotates standard convention to seismo

  R1 = rotation_matrix([1, 0, 0], -pi)
  R2 = rotation_matrix([0, 0, -1], -pi)
  R = R2*R1
  Rinv = R'

  #Calculate and rotate moment tensor to standard convention.

  M = moment_tensor(phi, delta, lamda)
  M = Rinv*M*Rinv'

  #Calculate P, N, and T axes in standard convention and rotate to seismo

  P, N, T = moment_tensor_2_PTN(M)

    #Convert to seismo
  P = R*P
  N = R*N
  T = R*T

  Rcheck = rotation_matrix([0, 0, 1], -pi)

  #rotates eigenvectors by 180

  #Ensures fault normal vectors will be pointing upwards


  if dip == 90
    if rake >= 90

      P = -P
      T = -T
      end
    end


  #Calculate and Rotate normal and slip vectors of fault and auxilary plane.

  normal, normal_aux, slip, slip_aux = PTN_2_fault_geometry(P, T, N)

    #= Rotates standard basis vectors to seismo basis vectors. This is down to
    properly calculate the rake =#

  x = [1, 0, 0]
  y = [0, 1, 0]
  z = [0, 0, 1]

  N = R*x
  E = R*y
  D = R*z


  np1_sd = strike_dip(normal)
  np2_sd = strike_dip(normal_aux)

  strikes = [np1_sd[1], np2_sd[1]]
  dips = [np1_sd[2], np2_sd[2]]
  slips = [slip, slip_aux]
  rakes = zeros(2)

  rakes[1] = acos(sin(dips[2])*sin(strikes[1]-strikes[2]))
  rakes[2] = acos(sin(dips[1]*sin(strikes[2]-strikes[1])))

  for i = 1:2
    if strikes[i] >= 2*pi
      strikes[i] = strikes[i] - 2*pi
    end

    if rakes[i] == 0
      strikes[i] = strikes[i] - pi
    end
  end



  np1 = [strikes[1], dips[1], rakes[1]]
  np2 = [strikes[2], dips[2], rakes[2]]

  return [np1, np2]

end
