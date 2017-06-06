def stressorientation_matrix(phi_MCS, theta_MCS, rho_MCS=0):

    Rstrike_MCS = np.zeros((3,3))
    Rstrike_MCS[0, 0] = np.cos(phi_MCS)
    Rstrike_MCS[0, 1] = -(np.sin(phi_MCS))
    Rstrike_MCS[1, 0] = np.sin(phi_MCS)
    Rstrike_MCS[1, 1] = np.cos(phi_MCS)
    Rstrike_MCS[2, 2] = 1
    Rstrike_MCS=np.asmatrix(Rstrike_MCS)

    Rdip_MCS = np.zeros((3,3))
    Rdip_MCS[0, 0] = np.cos(theta_MCS)
    Rdip_MCS[0, 2] = np.sin(theta_MCS)
    Rdip_MCS[1, 1] = 1
    Rdip_MCS[2, 0] = -(np.sin(theta_MCS))
    Rdip_MCS[2, 2] = np.cos(theta_MCS)
    Rstrike_MCS=np.asmatrix(Rstrike_MCS)

    rhoaxis = (Rstrike_MCS)*(Rdip_MCS)*np.matrix('1;0;0')
    Rrho_MCS = rotationmat3D(rho_MCS, rhoaxis)

    StressOrientation = (Rrho_MCS)*(Rstrike_MCS)*(Rdip_MCS)

    return StressOrientation
