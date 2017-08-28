#
# included functions: primary author, date added to FMfcns.jl (see below for written/modification dates)
#
# Flatten_and_Stack: RWC, August 2017
# GetParam: EA Hetland, April 2017
# writebin: EA Hetland, April 2017
# EstimateKappa: EA Hetland, April 2017
# rotationmat3D: Bileschi/Lu, April 2017
# FaultOrientation_Matrix: Medina-Luna/Hetland, April 2017
# StressOrientation_Matrix: Medina-Luna/Hetland, April 2017
# Resolve_ShearStressDir: Medina-Luna, April 2017
# Resolve_StressTensorOnFault: Medina-Luna, April 2017
# SphereStep: Hetland/Lu, April 2017
# StressStep: Lu/Wilcox-Cline/Hetland, April 2017
# SphericalCoord:: Hetland, April 2017
# Compute_StressPtsAngles: Medina-Luna, April 2017
# Compute_StressPtsAng_PStress: Medina-Luna, April 2017
# LambertProjectS2: Hetland, April 2017
#
# use and modification history in each function
#

function flatten_and_stack(data, number_groups)
    
    data_1st_slice = data[:, :, 1]
    group = ones((size(data_1st_slice, 1), 1))
    data_stack = hcat(data_1st_slice, group)
    
    for id = 2:num_groups
    
        new_level = data[:, :, id]
        group = id*ones((size(new_level, 1), 1))
        data_group = hcat(new_level, group)
        data_stack = vcat(data_stack, data_group)
    
    end
    
    return data_stack
    
end

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

function GetParam(input, param, flag, default)
#
# gets a parameter value from the input control file
#
# EA Hetland, Univ. Michigan, April 2017
#
    lne = filter(x->startswith(lowercase(x),lowercase(param)),input)
    if length(lne)==0
        val = default
    else
        val = split(lne[1],"::")[2]
        if flag=="I"
            val = parse(Int,val)
        end
        if flag=="F"
            val = parse(Float64,val)
        end
    end
    println("\t",param," = ",val)
    return val
end

function writebin(name,ar);
#
# writes array into a flat binary file with first line the array size
#
# EA Hetland, Univ. Michigan, circa 2016
#
	 fid=open(name,"w");
	 write(fid,hcat([Float64(size(ar)[1]) Float64(size(ar)[2])],reshape(transpose(ar),1,length(ar))));
	 close(fid);
	 return 1
end
function EstimateKappa(x)
#
# estimates a kappa in the Von Mises PDF such that the standard
# deviation of the Von Mises PDF is equal to x
#
# polynomial fit using Mathematica (KappaVsStandardDeviationInversion.nb)
#
# x = desired circular standard deviation
#
# EA Hetland, Univ. Michigan, April 2017
#

    x = x*180/pi
    kappa = 3.03538 - 0.239742*x + 0.0120115*x^2 - 0.000374394*x^3 +
        6.73866e-6*x^4 - 6.37376e-8*x^5 + 2.43639e-10*x^6
    kappa = 10^kappa

    return kappa
end



function rotation_matrix_stress_utils(axis, rotation_angle)

  #=
   =======================================================================
  Creates a rotation matrix for a rotation around an arbitrary axis as
  described by the Euler-Rodrigues formula

  NOTE: Use this matrix to explicitly rotate about an arbitrary axis. This
  matrix is more easily used to transform  between Aki-Richards convention
  and standard RH coordinate systems than the other rotationmat3D function
   =======================================================================

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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# NOTE: OLD ROTATION MATRICES

#function rotationmat3D(radians, axis)
#
#function R= rotationmat3D(radians,axis)
#
#creates a rotation matrix such that R * x
#operates on x by rotating x around the origin r radians around line
#connecting the origin to the point "Axis"
#
#Example:
#rotate around a random direction a random amount and then back
#the result should be an Identity matrix

#r = rand(4,1);
#rotationmat3D(r(1),[r(2),r(3),r(4)]) * rotationmat3D(-r(1),[r(2),r(3),r(4)])

#example2:
#rotate around z axis 45 degrees
#Rtest = rotationmat3D(pi/4,[0 0 1])#

#ported from Matlab code by Bileschi 2009 (from MatlabCentral) to Julia by Larrissa Lu (UM/UROP) Summer 2016


    #NOTE: Will not be converting this part to julia
    #if nargin == 1
    #    if (length(rotX) == 3)
    #        rotY = rotX(2);
    #        rotX = rotZ(3);
    #        rotX = rotX(1);
    #    end
    #end

    #axisLength = vecnorm(axis);

    #if (axisLength < eps())
        #if L is less than eps(Float64), then send an error message
    #    println("error: axis direction must be a non-zero vector")

    #if there is no error, then run the rest of the function
    #else
        #create a unit vector and assign variables
    #    unitAxis = axis / axisLength;
    #    u = unitAxis[1];
    #    v = unitAxis[2];
    #    w = unitAxis[3];
    #    c = cos(radians);
    #    s = sin(radians);
#
#        #create a 3x3 matrix of zero
#        R = zeros(3,3);
#        R[1,1] = u^2 + (v^2 + w^2) * c;
#        R[1,2] = u * v * (1 - c) - w * s;
#        R[1,3] = u * w * (1 - c) + v * s;
#        R[2,1] = u * v * (1 - c) + w * s;
#        R[2,2] = v^2 + (u^2 + w^2) * c;
#        R[2,3] = v * w * (1 - c) - u * s;
#        R[3,1] = u * w * (1 - c) - v * s;
#        R[3,2] = v * w * (1 - c) + u * s;
#        R[3,3] = w^2 + (u^2 + v ^2) * c;
#    end
#    return R
#end

#function FaultOrientation_Matrix(phiF, thetaF, lambda = 0)
#FOrientation = FaultOrientation(strike,dip)
#FOrientation = FaultOrientation(strike,dip,rake)

#Computes a rotation matrix for a fault orientation of strike,
#dip, and optionally rake. All angles should in radians.

#INPUT: (in radians)
#phiF: strike of fault CCW from North, positive rotation.
#Negative value yields clockwise rotation from North.
#    i.e. rotation about the z-axis of phi, phi > 0 measured from the Y-axis

#thetaF: fault dip
#	i.e. Clockwise rotation about the y'-axis, theta > 0 measured from the
#    x'-axis

#lambda: fault rake
#	i.e. CCW rotation about the z''-axis of lambda,
#    lambda > 0 measured from y''-axis

#e.x., phi=30, theta=30 is a fault plane striking 30 CCW from N, and
#dipping 30 to the NE

# INPUT:
#   strike = strike of fault CCW from east
#            i.e. rotation about the z-axis of phi, phi>0 measured from the x-axis
# # dip = dip of fault
#       i.e. rotation about the y'-axis, dip>0 plunges from the x'-axis
#   rake = rake of slip vector
#          i.e. rotation about the z''-axis, rake>0 measured from y''-axis
#
# e.x., phi=45, theta=30 is a fault plane striking 45 from E, and
# dipping 30 to the NE
#
# OUTPUT:
#   R = rotation matrix that transforms a horizontal flaut with trace
#       along x-axis and extending to y>0 to the plane given by
#       (strike,dip), initial fault normal = [0 0 1], initial fault
#       slip vector = [1 0 0]

#Matlab code written by Lorena Medina-Luna & EA Hetland
#ported to Julia by Larissa Lu (UM/UROP) Summer 2017



#    Rstrike_Fault = [cos(phiF) -sin(phiF) 0 ;
#                     sin(phiF)  cos(phiF) 0 ;
#                     0          0         1];
#
#    Rdip_Fault    =[1 0 0;
#                    0 cos(thetaF)  sin(thetaF) ;
#                    0 -sin(thetaF) cos(thetaF)];
#
    #creates the axis vector to find rake matrix
#    rakeAxis = (Rstrike_Fault) * (Rdip_Fault) * [0;0;1];
#    Rrake_Fault = rotationmat3D(lambda, rakeAxis);
#
#    FaultOrientation = (Rrake_Fault) * (Rstrike_Fault) * (Rdip_Fault);
#
#    return FaultOrientation
#end

#function StressOrientation_Matrix(phi_MCS, theta_MCS, rho_MCS = 0)

#SOrientation = StressOrientation(strike,dipMCS)
#SOrientation = StressOrientation(strike,dipMCS,rhoMCS)

#Computes a rotation matrix for stress orientation of trend,
#plunge, and rotation about MCS. All angles should in radians.

#INPUT: (in radians)
#phi_MCS: trend of MCS - CCW rotation about the z-axis of phi,
#    phi >0 CCW rotation measured from the x-axis (E-W) towards North
#    phi <0 CW rotation measured from x-axis (E-W) towards South
#theta_MCS: plunge of MCS- CW rotation about the y'-axis
#    theta >0 measured from the x'-axis
# varagin: rho_MCS- rotation about the x''-axis, phi_MCS
#	rho_MCS >0 CCW rotation about phiMCS
#	rho_MCS <0 CW rotation about phiMCS

#e.x., phi_MCS=30, theta_MCS=30 is a stress rotated 30 deg from E toward the
#north, so that the MCS is trending N60E and plunging 30NE

#Matlab code written by Lorena Medina-Luna & EA Hetland
#ported to Julia by Larissa Lu (UM/UROP) Summer 2017



#    Rstrike_MCS = [cos(phi_MCS) -sin(phi_MCS) 0 ;
#                   sin(phi_MCS)  cos(phi_MCS) 0 ;
#                   0          0         1      ];

#    Rdip_MCS    = [cos(theta_MCS) 0 sin(theta_MCS) ;
#                         0        1        0       ;
#                  -sin(theta_MCS) 0 cos(theta_MCS)];

#    rhoAxis = (Rstrike_MCS) * (Rdip_MCS) * [1;0;0];
#    Rrho_MCS = rotationmat3D(rho_MCS, rhoAxis);

#    StressOrientation = (Rrho_MCS) * (Rstrike_MCS) * (Rdip_MCS);

#    return StressOrientation


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function FaultOrientation_Matrix(strike, dip)

    Rstrike = rotation_matrix_stress_utils([0, 0, 1], strike)
    yp = Rstrike*[0, 1, 0]
    Rdip = rotation_matrix_stress_utils(yp, dip)
    Rfault = Rdip*Rstrike

    return Rfault

end



function StressOrientation_Matrix(phi, theta, rho)

    Rphi = rotation_matrix_stress_utils([0, 0, 1], phi)
    yp = Rphi*[0, -1, 0]
    Rtheta = rotation_matrix_stress_utils(yp, theta)
    Ri = Rtheta*Rphi
    xpp = Ri*[1, 0, 0]
    Rrho = rotation_matrix_stress_utils(xpp, rho)
    Rstress = Rrho*Ri

    return Rstress

end

function Resolve_ShearStressDir(pstress, sangle, fangle)
#=
Calculate value of FaultRake = arctan(TauS/TauD)
Calculate Tau_rake_parallel (shear slip in direction of fault rake)

NOTE: in this function, TauD and TauS seem like the equivalent of TRParallel and
TRPerpendicular in Resolve_StressTensorOnFault

Matlab code written by Lorena Medina
ported to Julia by Larissa Lu (UM/UROP) Summer 2017


=#

    #ensures input fault rake = 0
    #I changed this part a bit from the MATLAB script because you can't access
    #or modify anything out of the bounds of an array
    if length(fangle) == 3
        if fangle[3] != 0
            println("changing non-zero rake to zero")
            fangle[3] = 0;
        end
    #if there are only two elements, add 0 as the third element
    elseif length(fangle) == 2
        push!(fangle, 0);
    end


    #using function "Stress_onFault" to determine TauD and TauS
    TauDip, TauStrike =
        Resolve_StressTensorOnFault(pstress, sangle, fangle);


    RakeCalc = atan2(TauStrike, TauDip);


    return RakeCalc
end


function Resolve_StressTensorOnFault(pstress, sangle, fangle)

#=
=======================================================================
  Resolve the shear stress acting on a fault plane <----Modified
                (old description:  Function to determine sense of
                 slip on an oriented fault and oriented stress tensor)

Tau_rake_parallel = TauD  when  Fault_rake = 0
Tau_rake_perpendicular = TauS when Fault_rake = 0

#NOTE UPDATED RWC 06/14/17: See Comments


======================================================================
=#
    #set the values
    S1 = pstress[1];
    S2 = pstress[2];
    S3 = pstress[3];

    #fangle = [Fault_strike, Fault_dip, Fault_rake]
    phiF   = fangle[1]
    thetaF = fangle[2];
    rakeF  = fangle[3]

    if size(pstress) == (3,3)
        stressDirection = pstress;

    else
        #sangle = [Stress_strike_MCS, Stress_dip_MCS, Stress_dip_ICS]
        phi_MCS   = sangle[1];
        theta_MCS = sangle[2];
        theta_ICS = sangle[3];

        #define the principal stress tensor
        Sprincipal = [S1 0.0 0.0;
                      0.0 S2 0.0;
                      0.0 0.0 S3];

        #----Stress Geographic
        #rotation matrix for stress tensor
        stressRotation  = StressOrientation_Matrix(phi_MCS, theta_MCS, theta_ICS);
        #direction of stress (w/strke & dip)
        stressgeo = stressRotation * Sprincipal * stressRotation';
        stressDirection = (stressgeo'+stressgeo)./2
    end


        # RWC Eigenvalues were never used in this part of the code so commented out

        #A = (stressDirection' + stressDirection) ./ 2;

        #diagonalized tensor is too different from original
        #if maximum(abs(A[:] - stressDirection[:])) > (1e3 * eps())
        #    SigmaN = 0;
        #    Tau_rake_parallel = 0;
        #    Tau_rake_perpendicular = 0;
        #    D = eye(3);
        #    V = eye(3);
        #
        #diagonalized tensor is negligibly different from original
        #else
            #stressDirection = A;
            #compute the principal directions
            #(D, V) = eig(stressDirection);

            #concatenate [1;2;3] onto D to make a matrix with 2 columns
            #unsortedValues = hcat(D, [1; 2; 3]);
            #sort the eigenvalues in the matrix while reordering pos accordingly
            #sortedValues = sortrows(unsortedValues);

            #store the sorted eigenvalues in a temporary vector
            #temp = vec(sortedValues[:, 1]);
            #store the position values as integers
            #spos = round(Int64, sortedValues[:, 2]);
            #turn sorted values into a diagonal matrix
            #D = diagm(temp);

            #make a new empty temp matrix to store our reordered V
            #tempV = zeros(3,3);
            #sort eigenvectors that correlate to our sorted eigenvalues
            #for i = 1:3
            #    tempV[:, i] = V[:, spos[i]];
            #end
            #V = tempV;

        faultRotation = FaultOrientation_Matrix(phiF, thetaF)

        #rotate stress onto fault plane
        stressOnFault = faultRotation' * stressDirection * faultRotation;

        #RWC Changed to work with new rotation matrices
          #  Tau_rake_parallel = stressOnFault[1,3];
          #  Tau_rake_perpendicular = stressOnFault[2,3];

        Tau_rake_parallel = stressOnFault[2,3];
        Tau_rake_perpendicular = stressOnFault[1,3];

    #end


    return Tau_rake_parallel, Tau_rake_perpendicular
end

function SphereStep(alpha, theta, delta)
# α = [lon; lat] of current step; theta = azimuth of step; delta = angle of step size
#
# Julia code written by Larissa Lu (UM/UROP) Summer 2017 based on
# formulae worked out by EA Hetland in Mathematica
# (RandomWalkInSpherical.nb)
#


    lon = alpha[1]
    lat = alpha[2]

    arcTanLonY = cos(lat)cos(delta)sin(lon) + sin(delta) * (cos(theta)sin(lat)sin(lon) +
                                                    cos(lon)sin(theta));
    arcTanLonX = cos(lat)cos(lon)cos(delta) + cos(lon)cos(theta)sin(lat)sin(delta) -
                    sin(lon)sin(delta)sin(theta) +
                    sqrt(cos(lat)^2 * cos(delta)^2 +
                        2 * cos(lat)cos(delta)cos(theta)sin(lat)sin(delta) +
                        sin(delta)^2 * (cos(theta)^2 * sin(lat)^2 + sin(theta)^2)
                        );
    arcTanLatY = sqrt(cos(lat)^2 * cos(delta)^2 +
                        2 * cos(lat)cos(delta)cos(theta)sin(lat)sin(delta) +
                        sin(delta)^2 * (cos(theta)^2 * sin(lat)^2 + sin(theta)^2)
                        );
    arcTanLatX = -cos(delta)sin(lat) + cos(lat)cos(theta)sin(delta);

    #arctan in mathematica is (x,y) but atan2 is (y,x)
    LON = 2 * atan2(arcTanLonY, arcTanLonX);
    LAT = -(pi / 2) + atan2(arcTanLatY, arcTanLatX);

    return LON, LAT
end

function StressStep(S3coords,RScoords,S3stepSize,RSstepSize)
    #
    # code written by Larissa Lu (UM/UROP), Russel Wilcox-Cline, and EA Hetland
    #
    # first written Summer 2017, modified Fall 2016-Winter 2017,
    # expanded and ported to separate code by EA Hetland April 207


    # sample rho from a uniform distrubtion EAH: included the step
    # size into rho step instead of being hardwired EAH: changed
    # index as current step = k, not k-1, to reflect new code
    # structure rhostep = rand() * 2 * S3stepSize - S3stepSize
    RHOstepSize=S3stepSize;
    rhostep = (rand() * 2   - 1)*(RHOstepSize*pi/180)
    #rhostep = 0.0
    trial_rho = mod(S3coords[3] + rhostep, pi)

    # sample values for phi and theta
    # EAH: changed variable from phi to step_bearing to avoid confusion with use of phi as trend of MCS
    step_bearing = rand() * 2pi - pi



    #trial_theta_phi = collect(SphereStep(S3coords[k - 1, :], step_bearing, S3stepSize));
    # EAH: modified the step size in the theta,phi step to account for rhostep
    trial_theta_phi = collect(SphereStep(S3coords[1:2], step_bearing, S3stepSize*cos(pi*rhostep/RHOstepSize/2)));
    #trial_theta_phi = collect(SphereStep(S3coords[k,1:2], step_bearing, S3stepSize));


    # concatenate into a single vector. The vector S3angles is what will be
    # used as an input into Resolve_ShearStressDir
    S3angles = vec(trial_theta_phi)
    push!(S3angles, trial_rho)

    # do uniform step in random direction in relative stress plane
    ratio_step = [-10.0 -10.0]
    while  RScoords[1]+ratio_step[1] < 0.0 ||
        RScoords[1]+ratio_step[1] > 1.0 ||
        RScoords[2]+ratio_step[2] < 0.21 ||
        RScoords[2]+ratio_step[2] > 1.00

        step_bearing = rand() * 2pi - pi
        ratio_step = RSstepSize*[cos(step_bearing) sin(step_bearing)]

    end


    RSratios = RScoords + vec(ratio_step)

    return S3angles, RSratios

end


function SphericalCoords(A)
    #=
    %
    % X = SphericalCoords(A)
    %
    % if A = [strike dip], X=[x,y,z]
    % if A = [x,y,z], X=[strike,dip]
    %
    % each row of A is a separate sample
    %
    % all points returned are on the unit sphere, and all angles are in
    % radians
    %
    % dip = 0 is horizontal
    % dip = +90 is vertical down
    % dip = -90 is vertical up
    %

    original Matlab written by EA Hetland, circa 2012, translated to Julia April 2017

=#

    if size(A,2)==3
        X = zeros(size(A,1),2);
        r = sqrt(sum(A.^2,2));
        X[:,1] = atan2(A[:,2],A[:,1]);
        X[:,2] = acos(A[:,3]./r) - pi/2;
    elseif size(A,1)==2
        X = zeros(size(A,1)-1,3);
        r = 1;
        #% change dip convention from dip=(0,90,180) from (vertical
        #% up,horizontal,vertical down) to (horizontal,vertical
        #% down,horizontal)
        A2 = pi/2 + A[2];
        X[1] = r*cos(A[1])*sin(A2);
        X[2] = r*sin(A[1]).*sin(A2);
        X[3] = r*cos(A2);
    end

    return X

end


function Compute_StressPtsAngles(V);
    #=
    % Compute the principal directions and orient points
    % lower hemisphere projection using the eigenvectors
    % of a stress tensor
    % ex:
    % [V D] = eig(StressTensor);
    %     [tmp,spos]=sort(diag(D));
    %     D = diag(tmp);
    %     V = V(1:3,spos);
    %

    original matlab code written by Lorena Medina-Luna (UM) circa 2012 or 2013
    translated to Julia by EA Hetland, April 2017

    =#


    PtsMCS = V[:,1]';
    PtsICS = V[:,2]';
    PtsLCS = V[:,3]';

    # convert to spherical coordinates (angles)
    AngMCS = SphericalCoords(PtsMCS);
    AngICS = SphericalCoords(PtsICS);
    AngLCS = SphericalCoords(PtsLCS);


    pos = find(x->(x >=  pi/2), AngMCS[:,1]);
    neg = find(x->(x <= -pi/2), AngMCS[:,1]);
    AngMCS[pos,1] = AngMCS[pos,1] - pi;
    AngMCS[neg,1] = AngMCS[neg,1] + pi;
    AngMCS[pos,2] = -AngMCS[pos,2];
    AngMCS[neg,2] = -AngMCS[neg,2];

    neg = find(x->(x <= 0.0), AngICS[:,1]);
    AngICS[neg,1] =  AngICS[neg,1] + pi;
    AngICS[neg,2] = -AngICS[neg,2];


    neg = intersect(
                   find(x->(x <= 0.0), AngLCS[:,1]),
                   find(x->(x <= 0.0), AngLCS[:,2])
                   )
    pos = intersect(
                   find(x->(x >= 0.0), AngLCS[:,1]),
                   find(x->(x <= 0.0), AngLCS[:,2])
                   )
    AngLCS[neg,1] =  AngLCS[neg,1] + pi;
    AngLCS[pos,1] =  AngLCS[pos,1] - pi;
    AngLCS[neg,2] = -AngLCS[neg,2];
    AngLCS[pos,2] = -AngLCS[pos,2];


    #% Plot lower hemisphere
    ps = find(x->(x > 0.0), PtsMCS[:,3]);
    PtsMCS[ps,1:3] = -PtsMCS[ps,1:3];

    ps = find(x->(x > 0.0), PtsICS[:,3]);
    PtsICS[ps,1:3] = -PtsICS[ps,1:3];

    ps = find(x->(x > 0.0), PtsLCS[:,3]);
    PtsLCS[ps,1:3] = -PtsLCS[ps,1:3];

    return (AngMCS, AngICS, AngLCS, PtsMCS, PtsICS, PtsLCS)

end

function  Compute_StressPtsAng_PStress(S1, Sratio, Sangle);
    #=
    %  function to calculate for points and angles given

    % --------
    % INPUT:
    % S1 = initial S1 value as a vector, same length as TestNum
    % Sratio = [Delta R3]
    %    R3 = linear values of the ratio of LCS to MCS (S3/S1)
    %    Delta = Ratio values of intermediate Mohr circle to large Mohr circle: (S2 - S3) / (S1 - S3)
    %Sangle = Principal stress orientations, [Phi_MCS, Theta_MCS, Rho_MCS]
    % TestNum = length of models
    %
    % OUTPUT:
    % M/I/LCS_Pts = X, Y, Z points in a Lambert Lower Hemisphere projection
    % M/I/LCS_Ang = strike, dip angles corresponding to the lower hemisphere points
    %  ------------

    original matlab code written by Lorena Medina-Luna (UM) circa 2012 or 2013
    translated to Julia, and modified to return Lambert (x,y) points, by EA Hetland, April 2017

    =#

    TestNum = size(Sangle,1)

    Delta = Sratio[:,1];
    R3 = Sratio[:,2];

    MCS_Ang = zeros(TestNum,2);
    ICS_Ang = zeros(TestNum,2);
    LCS_Ang = zeros(TestNum,2);

    MCS_Pts = zeros(TestNum,3);
    ICS_Pts = zeros(TestNum,3);
    LCS_Pts = zeros(TestNum,3);


    S2 = S1.*(Delta + R3 - (Delta .*R3));
    S3 = S1.*R3;

    phi_MCS = Sangle[:,1];
    theta_MCS = Sangle[:,2];
    rho_MCS = Sangle[:,3];

    for ii = 1:TestNum;
        Sprincipal = [
                      S1 0.0 0.0;
                      0.0 S2[ii] 0.0;
                      0.0 0.0 S3[ii]
                      ];

        #% Rotation matrix for Stress Tensor
        StressRotation = StressOrientation_Matrix( phi_MCS[ii], theta_MCS[ii], rho_MCS[ii] );
        #% Stress Geographic
        #% Direction of Stress (w/strike & dip)
        Stress_Direction = StressRotation*Sprincipal*StressRotation';

        Stress_Direction = (Stress_Direction'+Stress_Direction)./2;


        #% compute the principal directions
        (D, V) = OrderedEigenSystem(Stress_Direction);

        # this is alternative to computing the eigensystem, but commented out; there are non-zeros nums close to zero
        #V = zeros(3,3)
        #V[:,1] = StressRotation*[1.0;0.0;0.0];
        #V[:,2] = StressRotation*[0.0;1.0;0.0];
        #V[:,3] = StressRotation*[0.0;0.0;1.0];
        #for vi=1:3
        #    V[:,vi] = V[:,vi]/sqrt(sum(V[:,vi].^2));
        #end
        #println(Ve)
        #println(V)

        (AngMCS, AngICS, AngLCS, PtsMCS, PtsICS, PtsLCS) = Compute_StressPtsAngles(V);

        MCS_Pts[ii,:] = PtsMCS;
        ICS_Pts[ii,:] = PtsICS;
        LCS_Pts[ii,:] = PtsLCS;
        MCS_Ang[ii,:] = AngMCS;
        ICS_Ang[ii,:] = AngICS;
        LCS_Ang[ii,:] = AngLCS;

    end

    MCS_Pts = LambertProjectS2(MCS_Pts);
    ICS_Pts = LambertProjectS2(ICS_Pts);
    LCS_Pts = LambertProjectS2(LCS_Pts);

    return (MCS_Pts, ICS_Pts, LCS_Pts, MCS_Ang, ICS_Ang, LCS_Ang)

end

function LambertProjectS2(X)
    #
    # project lower hemisphere points with Lambert equal area projections
    #
    # EA Hetland, Univ. Michigan, Apr 2017

    x = repmat(sqrt(2./(1-X[:,3])),1,2).*X[:,1:2]

    return x

end
