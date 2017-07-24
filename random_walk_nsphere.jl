include("Distributions")

function weight(kappa, dim)

    #Follows algorithm from Wood 1994

    S = dim - 1
    b = (-2*kappa+(4*kappa^2+(S)^2)^(1/2))/(S)
    x = (1 - b) / (1 + b)
    c = kappa*x+S*log(1-x^2)

    while true

    beta = Beta(S/2, S/2)
    z = rand(beta, 1)
    w = (1-(1 + b)*z)/(1-(1-b)*z)
    w = w[1]
    u = rand()
    logu = log(u)
    test_val = kappa*w+S*log(1-x*w)-c
    if test_val > log(u)
        return w
    end
  end
  return w
end

function orthonormal_projection(mu)

    v = randn(size(mu))
    projection_mu_v = mu*(vecdot(mu,v))/norm(mu)
    ortho = v - projection_mu_v

    return ortho/norm(ortho)
end

function generate_sample(mu, kappa, nsteps)

  #Generates a random vector on an n-sphere

  dim = length(mu)
  random_vector = zeros((nsteps, dim))

  for i = 1:nsteps

    #generate sample off of mean direction
    w = weight(kappa, dim)

    #sample a point v on the unit sphere that's orthogonal to mu
    v = orthonormal_projection(mu)

    # compute new point
    random_vector[i, :] = v*sqrt(1-w^2)+w*mu
    end


    return random_vector
  end
