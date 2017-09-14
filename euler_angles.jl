include("FMfcns_RWC072117_fixedauxplane.jl")
include("StressInvUtils_RWC_08282017.jl")
include("old.jl")
include("rotationmatrix_to_euler.jl")
include("starting_angles.jl")



strike = 0
dip = pi/4
rake = pi/2


M = moment_tensor(strike, dip, rake)

P, N, T = moment_tensor_2_PTN(M)




rotm = [P N T]


angles = rotation_to_euler(rotm)

size(angles,2)

sigma_tensor = [-1 0 0; 0 -1/2 0; 0 0 -1/4]

Rx1 = [1 0 0; 0 cos(angles[1, 1]) -sin(angles[1, 1]); 0 sin(angles[1, 1]) cos(angles[1, 1])]
Ry1 = [cos(angles[2, 1]) 0 sin(angles[2, 1]); 0 1 0; -sin(angles[2, 1]) 0 cos(angles[2, 1])]
Rz1 = [cos(angles[3, 1]) -sin(angles[3, 1]) 0; sin(angles[3, 1]) cos(angles[3, 1]) 0; 0 0 1]

R1 = Rz1*Ry1*Rx1

sigma_geo1 = R1*sigma_tensor*R1'

eig(sigma_geo1)






Rx2 = [1 0 0; 0 cos(angles[1, 2]) -sin(angles[1, 2]); 0 sin(angles[1, 2]) cos(angles[1, 2])]
Ry2 = [cos(angles[2, 2]) 0 sin(angles[2, 2]); 0 1 0; -sin(angles[2, 2]) 0 cos(angles[2, 2])]
Rz2 = [cos(angles[3, 2]) -sin(angles[3, 2]) 0; sin(angles[3, 2]) cos(angles[3, 2]) 0; 0 0 1]

R2 = Rz2*Ry2*Rx2

M



sigma_geo1


sigma_geo2 = R2*sigma_tensor*R2'

check = round(rotm, 3)

Pcheck2, Ncheck2, Tcheck2 = moment_tensor_2_PTN(sigma_geo2)
calculated_axes2 = round([Pcheck2'; Ncheck2'; Tcheck2'],3)






Pcheck1, Ncheck1, Tcheck1 = moment_tensor_2_PTN(sigma_geo1)
calculated_axes2 = round([Pcheck1'; Ncheck1'; Tcheck1'], 3)




round(rotm, 3)





test = starting_direction(M)
