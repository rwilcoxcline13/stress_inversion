include("FMfcns_RWC072117_fixedauxplane.jl")
include("StressInvUtils_RWC_08282017.jl")
include("old.jl")
include("rotationmatrix_to_euler.jl")

function rotation_to_euler(rotm)

    #= Function to calculate the approximate orientation angles from P,N,T axes
    Input: Rotation Matrix (Matrix of Eigenvectors)
    Output: angles = [rotation_about_x; rotation_about_y; rotation_about_z]
    =#

    R31 = rotm[3, 1]

    if abs(R31) != 1

        theta1 = -asin(R31)
        theta2 = pi-theta1

        psi1 = atan2(rotm[3,2]/cos(theta1), rotm[3, 3]/cos(theta1))
        psi2 = atan2(rotm[3,2]/cos(theta2), rotm[3, 3]/cos(theta2))

        phi1 = atan2(rotm[2, 1]/cos(theta1), rotm[1,1]/cos(theta1))
        phi2 = atan2(rotm[2, 1]/cos(theta2), rotm[1,1]/cos(theta2))

        angles = [psi1 psi2; theta1 theta2; phi1 phi2]

    else

        phi = 0

        if R31 == -1

            theta = pi/2
            psi = phi+atan2(rotm[1,2], rotm[1, 3])

            angles = [psi, theta, phi]

        else

            theta = -pi/2
            psi = -phi+atan2(-rotm[1, 2], -rotm[1, 3])

            angles = [psi, theta, phi]

        end

    end

    return angles

end

function standard_rotation_matrix(psi, theta, phi)

    Rx = [1 0 0; 0 cos(psi) -sin(psi); 0 sin(psi) cos(psi)]
    Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)]
    Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1]

    R = Rz*Ry*Rx

    return R

end


function  starting_direction(M)

    P, N, T = moment_tensor_2_PTN(M)
    rotm = [P N T]
    angles = rotation_to_euler(rotm)

    sigma_principal = [-1 0 0; 0 -1/2 0; 0 0 -1/4]

    if size(angles, 2) == 2

    for i=1:2

        R = standard_rotation_matrix(angles[1, i], angles[2, i], angles[3, i])
        sigma_geo = R*sigma_principal*R'
        Pcheck, Ncheck, Tcheck = moment_tensor_2_PTN(sigma_geo)
        calculated_axes = round([Pcheck'; Ncheck'; Tcheck'],3)

        if calculated_axes == round(rotm, 3)


            starting_angles = angles[:, i]


        end
    end
end


    return starting_angles

end
