tic()

include("FMfcns_RWC072117_fixedauxplane.jl");
include("StressInvUtils_RWC_09142017.jl");
using Distributions;


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
ManualOrientation = GetParam(input, "ManualOrientation", "I", 0)
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

FM = Array(focalmechs[:,4:6]);
FMG = Array(focalmechs[:, end]);

#Sort focal mechanism by groups and determine the indicies of the focal mechanisms in each group


sorted_idx  = sortperm(FMG);

FMG = FMG[sorted_idx];
FM = FM[sorted_idx, :];

uni = unique(FMG);
uni_idx = indexin(uni, FMG);
num_groups = length(uni_idx);

group_idx = Array(Vector{Int}, num_groups)




for gr = 1:num_groups

    if gr == 1

        group_idx[gr] = 1:uni_idx[gr]

    else

        group_idx[gr] = uni_idx[gr-1]+1:uni_idx[gr]

    end

end



mean_orientation = zeros(3, 1, num_groups)



for group = 1:num_groups


    idx = group_idx[group]


    str = FM[idx, 1]*pi/180

    di = FM[idx, 2]*pi/180

    la = FM[idx, 3]*pi/180


    num_focs_group = length(str)

    PTN_orientation = zeros(3, num_focs_group)


    for fm = 1:num_focs_group

        M = moment_tensor(str[fm], di[fm], la[fm])
        P, N, T = moment_tensor_2_PTN(M)
        rotm = [P N T]
        angles = rotation_to_euler(rotm)
        sigma_principal = [-1 0 0; 0 -1/2 0; 0 0 -1/4]
        PTN_orientation[:, fm] = (angles[:, 1])

    end

    mean_orientation[:, 1, group] = mean(PTN_orientation, 2)

end

FMn = size(FM,1);
FaultPlanes = zeros(FMn,4);
RakesObserved = zeros(FMn,2);
RakeErrors = repmat(focalmechs[:,10],1,2).*pi./180.;
DipErrors = repmat(focalmechs[:,8],1,2).*pi./180.;
MLE_estimates = zeros(num_groups, 5);

for fi = 1:FMn
    fault_parameters = FM[fi, :];


    NP1, NP2 = fault_parameters_2_nodal_planes(fault_parameters[1], fault_parameters[2], fault_parameters[3]);

    FaultPlanes[fi,:] = [NP1[1] NP1[2] NP2[1] NP2[2]];
    RakesObserved[fi,:] = [NP1[3] NP2[3]];



    if float(fault_parameters[1]) != round(NP1[1],2 )

      FaultPlanes[fi, :] = [NP2[1] NP2[2] NP1[1] NP1[2]]
      RakesObserved[fi,:] = [NP2[3] NP1[3]]

    end

end


kappa_rake = zeros(FMn,2)
println(size(kappa_rake))
kappa_dip = zeros(FMn,2)


for j = 1:FMn
    for n = 1:2
        kappa_rake[j,n] = EstimateKappa(RakeErrors[j,n])
        kappa_dip[j,n] = EstimateKappa(DipErrors[j,n])
    end
end

ndips = 10^4
sampled_dips = zeros(FMn, ndips, 2);

for j = 1:FMn

    for p = 2:2:4
        v_dist = VonMises(FaultPlanes[j, p],kappa_dip[j, Int(p/2)])
        sampled_dips[j, :, Int(p/2)] = rand(v_dist, ndips);
    end

end



# change to operate on one nodal plane at a time
if nodalplanes == 0
    NPn = 2;

else
    NPn = 1;
    FaultPlanes = FaultPlanes[:,(nodalplanes-1)*2+1:nodalplanes*2];
    RakesObserved = RakesObserved[:,nodalplanes];
    RakeErrors = RakeErrors[:,nodalplanes];
    DipErrors = DipErrors[:, nodalplanes];
    kappa_rake = kappa_rake[:, nodalplanes];
    kappa_dip = kappa_dip[:,nodalplanes];
end

# re-order for ease in below computations
kappa_dip = kappa_dip';
RakesObserved = RakesObserved';


S3coords = zeros(bufferSize, 3, length(uni_idx));
RScoords = zeros(bufferSize, 2, length(uni_idx));


if ManualOrientation == 0

    S3coords[1, 1, :] = mean_orientation[3];
    S3coords[1, 2, :] = mean_orientation[2];
    S3coords[1, 3, :] = mean_orientation[1];

else



    S3coords[1, 1, :] = S3init[1];
    S3coords[1, 2, :] = S3init[2];
    S3coords[1, 3, :] = S3init[3];

end

RScoords[1,1,:] = RSinit[1];
RScoords[1,2,:] = RSinit[2];

S1 = -1;


Calculated_Rakes = zeros(bufferSize,NPn,FMn);
Optimal_Nodal_Plane = zeros(bufferSize, FMn);

#Start of MCMC
#Start the iteration over the total number of focal mechanism groups

last_likelihood_exponent = 0.0;

for id = 1:num_groups


    S3angles = S3coords[1,:, id]
    RSratios = RScoords[1,:, id]


    if id == 1
        grouped_focs = 1:uni_idx[id]


    else
        grouped_focs = uni_idx[id-1]+1:uni_idx[id]

    end



    for iteration = 1:numIterations

        for k = 1:bufferSize

            pstress = [S1; S1*(RSratios[2] + RSratios[1] - RSratios[2]*RSratios[1]); S1*RSratios[2]]

            # Number of focal mechanisms in a group
            foc_idx = group_idx[id]

            #Iterate through each focal mechanism in a group



            PredictedRakes = zeros(NPn, length(foc_idx))
            fault_residuals = zeros(length(foc_idx))



            for j = foc_idx

                #This loop should ouput two predicted rakes (1 for each nodal plane)



                #Each nodal plane wiill have a distribution of dips
                rake_distribution_from_dip = zeros(ndips, 2)

                #Iterate through each nodal plane


                for n = 1:NPn

                    #Iterate through each dip
                    for d = 1:size(sampled_dips, 2)

                        fp = [FaultPlanes[j, (n-1)*2+1] sampled_dips[j, d, n]]
                        rake_distribution_from_dip[d, n] = Resolve_ShearStressDir(pstress, S3angles, vec(fp))


                    end


                    PredictedRakes[n,j] = mean(rake_distribution_from_dip[:, n])


                end


                (fault_residuals[j],Optimal_Nodal_Plane[k,j]) = findmax(kappa_rake[:, j].*cos(PredictedRakes[:, j]- RakesObserved[:, j]))

                #println(findmax(kappa_rake[:, j].*cos(PredictedRakes- RakesObserved)))
                #println("FUCK")

            end




        #(fault_residuals[k, :, j],Optimal_Nodal_Plane[k,j]) = findmax(kappa_rake[:, j].*cos(PredictedRakes- RakesObserved))

            trial_likelihood_exponent = sum(fault_residuals)


            LR = exp(trial_likelihood_exponent - last_likelihood_exponent)




            if LR>=rand() || k==1





            # accept new step
                    S3coords[k, :, id] = S3angles
                    RScoords[k, :, id] = RSratios
                    Calculated_Rakes[k,:,grouped_focs] = PredictedRakes[:, grouped_focs]

                    last_likelihood_exponent = trial_likelihood_exponent

                else

            # reject and stay on previous step
                    S3coords[k, :, id] = S3coords[k-1, :, id]
                    RScoords[k, :, id] = RScoords[k-1, :, id]
                    Calculated_Rakes[k, :, grouped_focs] = Calculated_Rakes[k-1, :, grouped_focs]



                end

                # select a new stress sample for next step
                # this is not needed if on k==numSteps, but just avoid the if structure to evaluate every loop

                (S3angles,RSratios) = StressStep(vec(S3coords[k,:, id]),vec(RScoords[k,:, id]),S3stepSize,RSstepSize);

                step_bearing = rand() * 2pi - pi


        end

    end

    (MCS_Pts, ICS_Pts, LCS_Pts, MCS_Ang, ICS_Ang, LCS_Ang) = Compute_StressPtsAng_PStress(-1,RScoords[:, :, id],S3coords[:, :, id]);
    stress_output = hcat(RScoords,S3coords, MCS_Pts, ICS_Pts, LCS_Pts, MCS_Ang, ICS_Ang, LCS_Ang)
    writedlm("/home/rwcline/Desktop/stress_params.txt", stress_output[:, :, id], ",")


    for j = 1:FMn
    rake_output = hcat(Calculated_Rakes[:,:,j],Optimal_Nodal_Plane[:,j])
    writedlm(string("/home/rwcline/Desktop/rakes_", j, ".txt"), rake_output, ", ")
    end
end
toc()
