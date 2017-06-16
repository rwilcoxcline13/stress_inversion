include("FMfcns_new.jl")
include("StressInvUtils2_CORRECT_902pm.jl")

if length(ARGS)==0
    error("\n\tplease supply an input file\n")
end

ControlFile = ARGS[1]
println("\nreading input file: ",ControlFile,"\n")

fio = open(ControlFile)
input = readlines(fio)
close(fio)

# remove all blank lines
input = input[find(x->(x>1),map(length,input))]

# remove comment lines
input = input[find(map(x->(!ismatch(r"^#",x)),input))]

# strip off newline characters
input=map(chomp,input)

directory   = GetParam(input,"OutputDirectory","",homedir())
if directory[end]!='!'
    directory = string(directory,'/');
end
outfile     = GetParam(input,"OutputFileBase","","DefaultStressPosterior")
fmfile      = GetParam(input,"FocalMechanisms","","FocalMechanisms.csv")
fmgroup      = GetParam(input,"FocalMechanismGroups","","all")
nodalplanes = GetParam(input,"NodalPlanes","I",0)

bufferSize    = GetParam(input,"BufferSize","I",10000)
numIterations = GetParam(input,"NumberIterations","I",8)
S3stepSize    = GetParam(input,"S3StepSize","F",0.314)
RSstepSize    = GetParam(input,"RSStepSize","F",0.05)

S3init = [GetParam(input,"PhiInitial","F",0.0)
          GetParam(input,"ThetaInitial","F",0.0)
          GetParam(input,"RhoInitial","F",0.0)]
RSinit = [GetParam(input,"DeltaInitial","F",0.5)
          GetParam(input,"R3Initial","F",0.25)]

println("\n")

# read in focal mechanism file
# see format for focal mechanisms file, note the first line is a header line, no comments allowed
focalmechs = readcsv(fmfile)[2:end,:];

# window out desired focal mechanism group(s)
if fmgroup!="all"
    ps = []
    fmvec = split(fmgroup,"/")
    for fmi = 1:length(fmvec)
        append!(ps,find(x->(ismatch(Regex(fmvec[fmi]),x)),focalmechs[:,end]))
    end
    focalmechs = focalmechs[ps,:]
end



# FM is the geometry of the one of the focal mechanism nodal planes
FM = focalmechs[:,4:6];


FMn = size(FM,1);
RakeErrors = repmat(focalmechs[:,10],1,2).*pi./180.;

println("read in ",FMn," focal mechanisms with grouping(s) ",fmgroup)

# assemble the fault plane information and rakes
FaultPlanes = zeros(FMn,4)
RakesObserved = zeros(FMn,2);





# calculate auxiliary plane
for fi = 1:FMn
    fault_parameters = FM[fi, :];
    NP1, NP2 = fault_parameters_2_nodal_planes(fault_parameters[1], fault_parameters[2], fault_parameters[3]);
    FaultPlanes[fi,:] = [NP1[1] NP1[2] NP2[1] NP2[2]];
    RakesObserved[fi,:] = [NP1[3] NP2[3]];


    if float(fault_parameters[1]) != round(NP1[1],2 )

      FaultPlanes[fi, :] = [NP2[1] NP2[2] NP1[1] NP1[2]]
      RakesObserved[fi,:] = [NP2[3] NP1[3]]

    end

    RakeErrors[fi,:] = [10.0 10.0]*pi/180.0;
end

println(RakeErrors)

# estimake kappa in a Von Mises PDF which s.t.d. given by the RakeErrors
kappa = zeros(FMn,2)

for j = 1:FMn
    for n = 1:2
        kappa[j,n] = EstimateKappa(RakeErrors[j,n])
    end
end


# change to operate on one nodal plane at a time
if nodalplanes == 0
    NPn = 2;
    #FaultPlaneKappa = 18.0;
else
    NPn = 1
    FaultPlanes = FaultPlanes[:,(nodalplanes-1)*2+1:nodalplanes*2]
    RakesObserved = RakesObserved[:,nodalplanes];
    RakeErrors = RakeErrors[:,nodalplanes];
    kappa = kappa[:,nodalplanes];
end
FaultPlaneKappa = 68.0;
FPstepSize = 2*pi/180.0;


println(string("writing out nodal plane geometries to file base ",outfile))
writebin(string(directory,outfile,"-NodalPlanes",".bin"), hcat(FaultPlanes,RakesObserved,kappa))

# re-order for ease in below computations
kappa = kappa'
RakesObserved = RakesObserved'

#Convert to radians RWC 06142017



# inversion is insensitive to S1 w/o other stress constraints, keep dummy S1 here
S1 = -1;

#Define S3coords as the steps through an S3 model space.
#  S3coords = [phi, theta, rho]
#
#Define RScoords as the steps through the relative stress model space.
#  RScoords = [Delta, R3]
#              Delta = (S2-S1)/(S1-S3)
#              R3 = S3/S1

S3coords = zeros(bufferSize, 3)
RScoords = zeros(bufferSize, 2)

#Create empty array for rakes calculated from the model parameters

Calculated_Rakes = zeros(bufferSize,NPn,FMn)
Optimal_Nodal_Plane = zeros(bufferSize, FMn)

#Define coordinates for the initial step
S3coords[1, :] = S3init
RScoords[1,:] = RSinit

#Define empty arrays for MCMC accept/reject
# EAH: removed misfit_angle as a separate variable below, do not need to store this information
# EAH: misfit_angle = zeros(numSteps, number_of_faults)
PredictedRakes = zeros(NPn,FMn)
fault_residuals = zeros(FMn)

#
# start the MCMC interation
#
S3angles = vec(S3coords[1,:])
RSratios = vec(RScoords[1,:])



#
# below only calculate the likelihood of the fault plane with maximum rake likelihood, so just compute one NP here
#



last_fault_planes_likelihood_exponent = 0.0
for j = 1:FMn

    last_fault_planes_likelihood_exponent += FaultPlaneKappa*vecdot(SphericalCoords(FaultPlanes[j,1:2]),SphericalCoords(FaultPlanes[j,1:2]));
end
NewFaultPlanes = zeros(size(FaultPlanes));
NewFaultPlanes[:,:] = FaultPlanes[:,:];

println(NewFaultPlanes)



# set the initial likelihood to 1, first trial always saved in below selection
last_likelihood_exponent = 0.0;

for iteration = 1:numIterations
    for k = 1:bufferSize

        # calculate predicted rakes for all fault planes
        pstress = [S1; S1*(RSratios[2] + RSratios[1] - RSratios[2]*RSratios[1]); S1*RSratios[2]];
        for j = 1:FMn
            for n = 1:NPn

                PredictedRakes[n,j] = Resolve_ShearStressDir(pstress, S3angles, vec(NewFaultPlanes[j, (n-1)*2+1:n*2]));

            end
        end

        #Calculate likelihoods using VM-Fisher Distribution
        for j = 1:FMn
            (fault_residuals[j],Optimal_Nodal_Plane[k,j]) = findmax(kappa[:,j].*cos(RakesObserved[:,j]-PredictedRakes[:,j]))
        end
        trial_likelihood_exponent = sum(fault_residuals)


#=
        # calculate the likelihood of the fault planes from the prior
        new_fault_planes_likelihood_exponent = 0.0

        for j = 1:FMn
            np = round(Int,Optimal_Nodal_Plane[k,j]);
            new_fault_planes_likelihood_exponent += FaultPlaneKappa*vecdot(SphericalCoords(   FaultPlanes[j,(np-1)*2+1:np*2]),
                                                                           SphericalCoords(NewFaultPlanes[j,(np-1)*2+1:np*2]));
        end
=#


        # use Gaussian for testing
        #AngDiff = abs(mod(Estimated_Rakes + pi - Rake_Observed,2*pi)-pi);
        #trial_likelihood_exponent = sum(-0.5.*AngDiff.^2./(Rake_Error.^2))

        # find ratio of likelhood between this trial and last
        LR = exp(trial_likelihood_exponent - last_likelihood_exponent)


        # select or reject
        if LR>=rand() || k==1

            println("accept")
            # accept new step
            S3coords[k, :] = S3angles
            RScoords[k, :] = RSratios
            Calculated_Rakes[k,:,:] = PredictedRakes
            last_likelihood_exponent = trial_likelihood_exponent
            #last_fault_planes_likelihood_exponent = new_fault_planes_likelihood_exponent
        else

            # reject and stay on previous step
            S3coords[k, :] = S3coords[k-1, :]
            RScoords[k, :] = RScoords[k-1, :]
            Calculated_Rakes[k, :, :] = Calculated_Rakes[k-1, :, :]


        end




        # select a new stress sample for next step
        # this is not needed if on k==numSteps, but just avoid the if structure to evaluate every loop
        (S3angles,RSratios) = StressStep(vec(S3coords[k,:]),vec(RScoords[k,:]),S3stepSize,RSstepSize);



#=
        # take a new step for the fault planes
        step_bearing = rand() * 2pi - pi
        for j = 1:FMn
            for i = 1:2:2*NPn
                NewFaultPlanes[j,i:i+1] = collect(SphereStep(NewFaultPlanes[j,i:i+1], step_bearing, FPstepSize));
            end
        end
=#
    end



    println(iteration)
    (MCS_Pts, ICS_Pts, LCS_Pts, MCS_Ang, ICS_Ang, LCS_Ang) = Compute_StressPtsAng_PStress(-1,RScoords,S3coords);
    stress_output = hcat(RScoords,S3coords, MCS_Pts, ICS_Pts, LCS_Pts, MCS_Ang, ICS_Ang, LCS_Ang)

    println(string("writing out interation ",iteration," of ",numIterations," (",bufferSize," samples) to file base ",outfile))
    writebin(string(directory,outfile,"-",iteration,".bin"), hcat(RScoords,S3coords,ICS_Ang,LCS_Ang,MCS_Pts,ICS_Pts,LCS_Pts))
    writedlm("/home/rwcline/Desktop/stress_params.txt", stress_output, ", ")


    for j = 1:FMn
        writebin(string(directory,outfile,"-",iteration,"-PredictedRakes_FM-",j,".bin"),
                 hcat(Calculated_Rakes[:,:,j],Optimal_Nodal_Plane[:,j]))
        rake_output = hcat(Calculated_Rakes[:,:,j],Optimal_Nodal_Plane[:,j])
        writedlm(string("/home/rwcline/Desktop/rakes_", j, ".txt"), rake_output, ", ")
    end


end
