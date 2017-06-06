include("Sp2CartN.jl")
include("SphereStep.jl")
include("Resolve_ShearStressDir.jl")
include("FaultOrientation_Matrix.jl")
include("StressOrientation_Matrix.jl")
include("Resolve_StressTensorOnFault.jl")

tic()
numSteps = 10^6;
stepSize = pi/100;
focus = [pi/4 pi/2];
kappa = 100;
mu = Sp2CartN(focus);

#preallocate for the for-loop instead of using a while loop
coords = zeros(numSteps, 2);
coords[1, :] = [0 -pi/2];
xo = Sp2CartN(coords[1, :]);

Rake_Observed = 14.4;
Calculated_Rakes = zeros(numSteps, 1);

faultangle = vec([150*pi/180, 87*pi/180]);
S1 = -1;
S2 = -2;
S3 = -1/4;

pstress = [S1;S2;S3];
misfit_angle = zeros(numSteps);
likelihood = zeros(numSteps);

for k = 1:numSteps

  if k == 1

    modelparams = vec(coords[k, :])
    push!(modelparams, pi/5)
    Estimated_Rake = Resolve_ShearStressDir(pstress, modelparams, faultangle);
    Calculated_Rakes[k] = Estimated_Rake
    misfit_angle[k] = (Rake_Observed-Calculated_Rakes[k])*pi/180
    likelihood_step = exp(cos(misfit_angle[k]))
    likelihood[k] = likelihood_step

  else

    #randomly generate some float between -pi and pi
    phi = rand() * 2pi - pi;
    #collect the tuple as an array
    trial = collect(SphereStep(coords[k - 1, :], phi, stepSize));
    modelparams = vec(trial)
    push!(modelparams, pi/5)

    Estimated_Rake = Resolve_ShearStressDir(pstress, modelparams, faultangle);
    Calculated_Rakes[k] = Estimated_Rake
    misfit_angle[k] = (Rake_Observed-Calculated_Rakes[k])*pi/180
    likelihood_step = exp(cos(misfit_angle[k]))
    likelihood[k] = likelihood_step

    LR = likelihood[k]/likelihood[k-1]
    random_number = rand()

    if LR >= rand()

      coords[k, :] = trial[1:2, 1]

    else

      coords[k, :] = coords[k-1, :]
    end
end
end

X = zeros(numSteps, 3)

for i = 1:numSteps

  X[i, :] = Sp2CartN(coords[i, :])
end

writedlm("/home/rwcline/Desktop//model.txt", X, ", ")
writedlm("/home/rwcline/Desktop//model2.txt", coords, ", ")
