
include("SphereStep.jl")
include("Resolve_ShearStressDir.jl")
include("FaultOrientation_Matrix.jl")
include("StressOrientation_Matrix.jl")
include("Resolve_StressTensorOnFault.jl")

tic()
numSteps = 100;
stepSize = pi/10;
S1 = -1;
S2 = -1/2;
S3 = -1/4;
kappa = 500;

pstress = [S1;S2;S3];
number_of_faults = 1


faultangle = [pi/3 pi/3]



coords = zeros(numSteps, 2);
coords[1, :] = [pi/4 pi/4];

Rake_Observed = [20.4]
Calculated_Rakes = zeros(numSteps, number_of_faults)

misfit_angle = zeros(numSteps, number_of_faults);
residual = zeros(numSteps, number_of_faults);
likelihood = zeros(numSteps);




for k = 1:numSteps
  if k == 1

    for j = 1:number_of_faults

      modelparams = vec(coords[k, :])
      push!(modelparams, 0)
      Estimated_Rake = Resolve_ShearStressDir(pstress, modelparams, vec(faultangle[j, :]));
      Calculated_Rakes[k, j] = Estimated_Rake
      misfit_angle[k, j] = (Rake_Observed[j]-Calculated_Rakes[k, j])*pi/180
      residual[k, j] = cos(misfit_angle[k,j])
    end

    likelihood_step = exp(kappa*(sum(residual[k, :])));
    likelihood[k] = likelihood_step;

  else


    #randomly generate some float between -pi and pi
    phi = rand() * 2pi - pi;
    #collect the tuple as an array
    trial = collect(SphereStep(coords[k - 1, :], phi, stepSize));
  for j = 1:number_of_faults
    modelparams = vec(trial)
    push!(modelparams, 0)
    Estimated_Rake = Resolve_ShearStressDir(pstress, modelparams, vec(faultangle[j, :]));
    Calculated_Rakes[k, j] = Estimated_Rake
    misfit_angle[k, j] = (Rake_Observed[j]-Calculated_Rakes[k, j])*pi/180
    residual[k, j] = cos(misfit_angle[k,j])
  end


  likelihood_step = exp(kappa*(sum(residual[k, :])));
  likelihood[k] = likelihood_step;
  LR = likelihood[k]/likelihood[k-1]

  if LR >= rand()

    coords[k, :] = trial[1:2, 1]

  else

    coords[k, :] = coords[k-1, :]
    Calculated_Rakes[k, :] = Calculated_Rakes[k-1, :]
  end

  end
end

coords

writedlm("/home/rwcline/Desktop/Research/Stress_Inversion/stress_inversion_results/model2.txt", coords, ", ")
writedlm("/home/rwcline/Desktop/Research/Stress_Inversion/stress_inversion_results/model3.txt", Calculated_Rakes, ", ")

toc()
