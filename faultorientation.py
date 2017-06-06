
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

#lambda: fault rake,
#	i.e. CCW rotation about the z''-axis of lambda,
#    lambda > 0 measured from y''-axis

#e.x., phi=30, theta=30 is a fault plane striking 30 CCW from N, and
#dipping 30 to the NE




def fault_orientation_matrix(phiF, thetaF, lamda):

    Rstrike_fault = np.zeros((3,3))
    Rstrike_fault[0, 0] = np.cos(phiF)
    Rstrike_fault[0, 1] = -np.sin(phiF)
    Rstrike_fault[1, 0] = np.sin(phiF)
    Rstrike_fault[1, 1] = np.cos(phiF)
    Rstrike_fault[2, 2] = 1
    Rstrike_fault = np.asmatrix(Rstrike_fault)

    Rdip_fault = np.zeros((3, 3))

    Rdip_fault[0 , 0] = 1
    Rdip_fault[1, 1] = np.cos(thetaF)
    Rdip_fault[1, 2] = np.sin(thetaF)
    Rdip_fault[2, 1] = -(np.sin(thetaF))
    Rdip_fault[2, 2] = np.cos(thetaF)
    Rdip_fault = np.asmatrix(Rdip_fault)

    rakeaxis = (Rstrike_fault)*(Rdip_fault)*np.matrix('0; 0 ; 1')

    rakeaxis = np.array(rakeaxis).flatten()

    Rrake_fault = rotationmat3D(lamda, rakeaxis)

    Faultorientation = (Rrake_fault)*(Rstrike_fault)*(Rdip_fault)


    return Faultorientation
