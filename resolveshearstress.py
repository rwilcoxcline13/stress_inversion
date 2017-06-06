def Resolve_ShearStressDir(pstress, sangle, fangle):
    #Notes:
    #Calculate value of FaultRake = arctan(Tau_S/Tau_D)
    #Calculate Tau_rake_parallel (shear slip in direction of fault rake)
    # In this function Tau_D and Tau_S seem like the equivalent of TRParallel
    # and TRPerpendicular in Resolve_StressTensorOnFault
    #End Notes

    #Ensures input fault rake = 0

    if len(fangle) == 3:

        if fangle[2] != 0:

            print('Changing Non-Zero Rake to Zero')
            fangle[2] = 0

    elif len(fangle) == 2: #If there are only two elements, add a third element that is 0

        np.append(fangle, 0)

    #Using function "Stress_onFault" to determine Tau_{d} and Tau_{s}
    SigmaN, TauDip, TauStrike, D, V = Resolve_StressTensorOnFault(pstress, sangle, fangle)

    #Calculate the fault rake

    RakeCalc = np.arctan2(TauStrike, TauDip) * (180/np.pi)

    #Calculate Chear Slip in the Rake Direction

    Tau_rake_parallel = TauDip/np.cos(np.deg2rad(RakeCalc))

    return RakeCalc, SigmaN, Tau_rake_parallel, D, V, TauDip, TauStrike
