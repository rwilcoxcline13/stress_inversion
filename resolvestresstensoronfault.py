def Resolve_StressTensorOnFault(pstress, sangle, fangle):


    S1 = pstress[0]
    S2 = pstress[1]
    S3 = pstress[2]


    # fangle = [Fault_Strike, Fault_dip, Fault_rake]
    # phiF = Fault_Strike, thetaF = Fault_dip, rakeF = Fault_rake
    # if dealing with principal stresses

    phiF = fangle[0]
    thetaF = fangle[1]
    rakeF = fangle[2]


    if pstress.shape == (3,3):

        stressDirection = pstress

    else:

        phi_MCS = sangle[0]
        theta_MCS = sangle[1]
        theta_ICS = sangle[2]

        #define the principal stress tensor
        Sprincipal = np.array(([S1, 0, 0,], [0, S2, 0], [0, 0, S3]))

        #------Stress Geographic
        stressRotation = stressorientation_matrix(phi_MCS, theta_MCS, theta_ICS)

        stressDirection = stressRotation*Sprincipal*stressRotation.T

        A = (stressDirection.T + stressDirection) / 2

        #diagonlized tensor is too differnt from origianl

    if np.amax(np.abs(A[:]-stressDirection[:]))> (1e3*np.spacing(1)):
        SigmaN = 0
        Tau_rake_parallel = 0
        Tau_rake_perpendicular = 0
        D = np.identity(3)
        V = np.identity(3)
        #diagnolized tensor is negligibly different from the original

    else:
        stressDirection = A
        stressDirection = np.transpose(stressDirection)
        #compute the principal directions
        D, V = np.linalg.eig(stressDirection)

    #-----Stress on a fault

    # rotation matrix for a fault plane
        faultRotation = fault_orientation_matrix(phiF, thetaF, rakeF)

        #Rotate stress onto fault plane
        stressOnFault = np.transpose(faultRotation)*stressDirection*faultRotation
        SigmaN = stressOnFault[2, 2]
        Tau_rake_parallel = stressOnFault[0,2]
        Tau_rake_perpendicular = stressOnFault[1, 2]

        return SigmaN, Tau_rake_parallel, Tau_rake_perpendicular, D, V
