# THIS CODE PERFORMS A VARIABLE NEIGHBORHOOD SEARCH OF A MIP
# GIVEN AN INCUMBENT, AN LP RELAXATION SOLUTION, 
# A SET OF VARIABLES TO SEARCH AND A MIP CORRESPONDING TO IT

#using JuMP

function getBounds(master_var)

   prev_l = Any[]
   prev_u = Any[]
   for i = 1:length(master_var)
      prev_l_sub = Any[]
      prev_u_sub = Any[]
      for j = 1:length(master_var[i])
         push!(prev_l_sub,getLower(master_var[i][j]))
         push!(prev_u_sub,getUpper(master_var[i][j]))
      end
      push!(prev_l, prev_l_sub)
      push!(prev_u, prev_u_sub)
   end
   return prev_l, prev_u
end

function fixLocalSolutionNeighborhood(incumbent, delta, master_var, k_neighborhood)

   for i in 1:Int(k_neighborhood)
      ind1 = delta[i].ind1
      ind2 = delta[i].ind2
      setLower(master_var[ind1][ind2], round(incumbent[ind1][ind2]))
      setUpper(master_var[ind1][ind2], round(incumbent[ind1][ind2]))
   end

end

function setWarmStart(master_var, incumbent)
   for i = 1:length(master_var)
      for j = 1:length(master_var[i])
         setValue(master_var[i][j], round(incumbent[i][j]))
      end
   end
end

function updateIncumbent(incumbent, master_var)
   for i = 1:length(master_var)
      for j = 1:length(master_var[i])
         incumbent[i][j] = getValue(master_var[i][j])
      end
   end
end

function revertBounds(master_var, prev_l, prev_u)
   for i = 1:length(master_var)
      for j = 1:length(master_var[i])
         setLower(master_var[i][j], prev_l[i][j])
         setUpper(master_var[i][j], prev_u[i][j])
      end
   end
end

type sortVector
    value::Float64
    ind1::Int64
    ind2::Int64
    sortVector() = new(0.0, 0, 0)
    sortVector(x,y,z) = new(x, y, z)
end

function myComparator(a::sortVector, b::sortVector)
    if a.value < b.value
        return true
    else
        return false
    end
end

# FOLLOWING CODE PERFORMS A VNS OF mip_model STARTING FROM INCUMBENT
# EMPLOYS A LOCAL DESCENT STEP AND REINITILIZATION MECHANIC
function VariableNeighborhoodSearch(mip_model, incumbent, lp_relaxation, master_var, supp_incumbent, supp_var, MAX_CPUTIME)

   OPTIMALITY_TOL = 1e-3

   N = 0
   delta = Any[]
   for i in 1:length(master_var)
      N += length(master_var[i])
      for j = 1:length(master_var[i])
         push!(delta, sortVector())
      end
   end

   (prev_l, prev_u) = getBounds(master_var)

   prevObjValue = 1e+10
   reinitialize = false
   noimprovement = 0
   consecutiveReinitializationCounter = 0
   MAX_NO_IMPROVEMENT = 5
   MAX_REINITIALIZATION = 4

   start = time()
   iter = 0
   # MAIN VNS LOOP
   while (time() - start < MAX_CPUTIME)
      iter += 1
      # IF A NUMBER OF CONSECUTIVE REINITIALIZATIONS MADE, STOP, RETURN THE CURRENT INCUMBENT
      if consecutiveReinitializationCounter >= MAX_REINITIALIZATION
         break
      end

      # COMPUTE DEVIATION FROM LP RELAXATION, DELTA
      k = 1
      for i = 1:length(master_var)
         for j = 1:length(master_var[i])
            delta[k].value = incumbent[i][j] - lp_relaxation[i][j]
            delta[k].ind1 = i
            delta[k].ind2 = j
            k += 1
         end
      end

      # ORDER DELTA DEPENDING ON REINITILIZATION FLAG
      noimprovement = 0
      if reinitialize
	     println("REINITIALIZING...")
         reinitialize = false
         d = 0.5
         shuffle!(delta)  
         consecutiveReinitializationCounter += 1
      else
         d = 2.0
         sort!(delta, lt=myComparator)
      end

      # COMPUTE NONZERO COMPONENTS OF DELTA
      NumNonZero = 0
      for j = 1:N
         if (delta[j].value > 1e-4)
            NumNonZero += 1
         end
      end
      # COMPUTE NEIGHBORHOOD SIZE
      k_step = round(NumNonZero / d)
      k_neighborhood = N - k_step
      @show d, NumNonZero, k_step, N, k_neighborhood
      # @show k_neighborhood
      # LOCAL DESCENT LOOP
      while k_neighborhood >= 0
	     println("NEIGHBORHOOD SIZE: $(k_neighborhood) FREED NEIGHBORHOOD: $(k_step) CPUTIME: $(time()-start)")
   	     setSolver(mip_model, CplexSolver(CPX_PARAM_MIPDISPLAY=0,CPX_PARAM_TILIM=150,CPX_PARAM_CLOCKTYPE=2))
         fixLocalSolutionNeighborhood(incumbent, delta, master_var, k_neighborhood)
   	     # WARM START FROM INCUMBENT
     	 setWarmStart(master_var, incumbent)
         
         status = solve(mip_model)
         if ((status == :Optimal || status == :Feasible || status == :CPXMIP_TIME_LIM_FEAS) && getObjectiveValue(mip_model) < prevObjValue - OPTIMALITY_TOL)
            println("BETTER OBJECTIVE FOUND: $(getObjectiveValue(mip_model))")
            consecutiveReinitializationCounter = 0
            prevObjValue = getObjectiveValue(mip_model)

            updateIncumbent(incumbent, master_var)
            updateIncumbent(supp_incumbent, supp_var)
            revertBounds(master_var, prev_l, prev_u)
	        break
         else
            revertBounds(master_var, prev_l, prev_u)
            noimprovement += 1
            if noimprovement > MAX_NO_IMPROVEMENT
               reinitialize = true
               break
            end
            k_step = (floor(k_step / 2.0) > 1 ? floor(k_step / 2.0) : 1);
            k_neighborhood = k_neighborhood - k_step
         end
      end
      if k_neighborhood < 0
         reinitialize = true
      end 
   end
   return prevObjValue
end
