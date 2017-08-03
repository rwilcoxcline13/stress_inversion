using DataFrames
using PyPlot

function tril_min(matrix)

  idxs = []
  vals = []

  for row = 2:size(matrix, 1)

    col = range(1, 1, row-1)

    for c in col

      idx = [row, c]
      push!(idxs, [row, c])
      push!(vals, matrix[row, c])

    end

  end

  idx_min = indmin(vals)
  loc_mat_min = idxs[idx_min]

  return loc_mat_min
end



#load focal mechanism Data

FM_data = readtable("NAF_FocalMechs.csv")
lats, lons = FM_data[1:5, 1], FM_data[1:5, 2]
num_events = range(1, 1, length(lats))

#Create distance matrix
distance_matrix = sqrt((broadcast(-, reshape(lats, (length(lats), 1)), reshape(lats, (1, length(lats))))).^2+broadcast(-, reshape(lons, (length(lons), 1)), reshape(lons, (1, length(lons)))).^2)
distance_matrix = reshape(vec(distance_matrix), (5, 5))

min_loc = tril_min(distance_matrix)
distance_matrix[min_loc[1], min_loc[2]]


plot(lats, lons, marker = "o", linestyle = "none", color = "blue")
plot(lats[min_loc[1]], lons[min_loc[1]], marker = "o", color = "red", linestyle = "none")
plot(lats[min_loc[2]], lons[min_loc[2]], marker = "o", color = "red", linestyle = "none")

min_loc[1]

min_loc[2]


#Single Link






distance_matrix[:, min_loc[1]] = 0.5*distance_matrix[:, min_loc[1]]+0.5*distance_matrix[:, min_loc[2]]+0*distance_matrix[min_loc[1], min_loc[2]]+-0.5*abs(distance_matrix[:, min_loc[1]]-distance_matrix[:, min_loc[2]])
distance_matrix[min_loc[1], :] = 0.5*distance_matrix[:, min_loc[1]]+0.5*distance_matrix[:, min_loc[2]]+0*distance_matrix[min_loc[1], min_loc[2]]+-0.5*abs(distance_matrix[:, min_loc[1]]-distance_matrix[:, min_loc[2]])
distance_matrix = distance_matrix[1: 1:size(distance_matrix, 2).!=min_loc[2],1:size(distance_matrix, 2).!=min_loc[2]]

min_loc = tril_min(distance_matrix)
distance_matrix[min_loc[1], min_loc[2]]

distance_matrix[:, min_loc[1]] = 0.5*distance_matrix[:, min_loc[1]]+0.5*distance_matrix[:, min_loc[2]]+0*distance_matrix[min_loc[1], min_loc[2]]+-0.5*abs(distance_matrix[:, min_loc[1]]-distance_matrix[:, min_loc[2]])
distance_matrix[min_loc[1], :] = 0.5*distance_matrix[:, min_loc[1]]+0.5*distance_matrix[:, min_loc[2]]+0*distance_matrix[min_loc[1], min_loc[2]]+-0.5*abs(distance_matrix[:, min_loc[1]]-distance_matrix[:, min_loc[2]])
distance_matrix = distance_matrix[1: 1:size(distance_matrix, 2).!=min_loc[2],1:size(distance_matrix, 2).!=min_loc[2]]

min_loc = tril_min(distance_matrix)
distance_matrix[min_loc[1], min_loc[2]]

distance_matrix[:, min_loc[1]] = 0.5*distance_matrix[:, min_loc[1]]+0.5*distance_matrix[:, min_loc[2]]+0*distance_matrix[min_loc[1], min_loc[2]]+-0.5*abs(distance_matrix[:, min_loc[1]]-distance_matrix[:, min_loc[2]])
distance_matrix[min_loc[1], :] = 0.5*distance_matrix[:, min_loc[1]]+0.5*distance_matrix[:, min_loc[2]]+0*distance_matrix[min_loc[1], min_loc[2]]+-0.5*abs(distance_matrix[:, min_loc[1]]-distance_matrix[:, min_loc[2]])
distance_matrix = distance_matrix[1: 1:size(distance_matrix, 2).!=min_loc[2],1:size(distance_matrix, 2).!=min_loc[2]]

min_loc = tril_min(distance_matrix)
distance_matrix[min_loc[1], min_loc[2]]

distance_matrix[:, min_loc[1]] = 0.5*distance_matrix[:, min_loc[1]]+0.5*distance_matrix[:, min_loc[2]]+0*distance_matrix[min_loc[1], min_loc[2]]+-0.5*abs(distance_matrix[:, min_loc[1]]-distance_matrix[:, min_loc[2]])
distance_matrix[min_loc[1], :] = 0.5*distance_matrix[:, min_loc[1]]+0.5*distance_matrix[:, min_loc[2]]+0*distance_matrix[min_loc[1], min_loc[2]]+-0.5*abs(distance_matrix[:, min_loc[1]]-distance_matrix[:, min_loc[2]])
distance_matrix = distance_matrix[1: 1:size(distance_matrix, 2).!=min_loc[2],1:size(distance_matrix, 2).!=min_loc[2]]





















#distance_matrix = distance_matrix[1: 1:size(distance_matrix, 2).!=min_loc[1],1:size(distance_matrix, 2).!=min_loc[1]]
#distance_matrix
#distance_matrix[:, 2] = 0
#distance_matrix[2, :] = 0

#distance_matrix
#distance_matrix = distance_matrix[1: 1:size(distance_matrix, 2).!=min_loc[2],1:size(distance_matrix, 2).!=min_loc[2]]
