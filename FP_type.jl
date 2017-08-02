function FP_type(rake)

  #=
   fpt = FaultPlaneType(rake)
 [fpt,fptstring] = FaultPlaneType(rake)

 computes the a fault plane type parameter in range [-1,1] based on
 algorithm in Shearer et al. (2006), but modified to use the mean
 rake to better represent the total fault plane type for focal
 mechanisms that are highly oblique

 INPUT:
   rake = [rake1 rake2] each event on one line
          rake should be in radians

 OUTPUT:
   fpt = fault plane type parameter
         fpt = -1 normal fault
         fpt =  0 strike-slip
         fpt = +1 thrust
   fptstring = {'normal','strike-slip','thrust'}
              ftpstring{round(ftp)+2} is the closest fault type

EA Hetland, Univ. Michigan, July 2013
   ehetland@alum.mit.edu

Ported to Julia and updated to agree with Shearer: Russel Wilcox Cline August 2017

=#

  fault_types_names = ["Normal", "Strike-Slip", "Thrust"]
  rake = pi/180*rake
  println(rake)
  fpt = zeros(size(rake, 1))
  idx = find(rake .> pi/2)

  rake[idx] = (pi - abs(rake[idx])).*(rake[idx]./abs(rake[idx]))

  #Changed from fpt = -1.*mean(rake,2)*2/pi; in order to agree with Shearer
  fpt = mean(rake, 2)*2/pi

  fpt_idx = trunc(Int, round(fpt)+2)
  fault_types = fault_types_names[fpt_idx]

  return fpt, fault_types
end
