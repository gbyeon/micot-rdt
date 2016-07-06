using ArgParse

include("ORDGDP_pmap.jl")
include("solver.jl")

# TODO... need to specify the solver and the algorithm at the command line
# TOOO... need to make this a Julia module... will make things easier in the long run



# general function for parsing the command line and putting in some default values and parameters
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--file", "-f"
            help = "an argument for the reslient design json file"
            default = "data/34Bus_Ice_Harden_Damageable_Rural_70Percent.json"
        "--mip_solver", "-m"
          help = "an argument for the mixed integer programming solver json file." 
    end

    return parse_args(s)
end


# The main entry point into running a power flow solver
function main(parsed_args) 
  rdt_file = parsed_args["file"]  
  mip_solver_file = parsed_args["mip_solver"]  
  if mip_solver_file == nothing
    mip_solver = build_solver(CBC_SOLVER)  
  else
    mip_solver = build_solver_file(mip_solver_file)
  end  
          
  solveORDGDP(rdt_file, mip_solver)      
 end

# this prevents main from being run automatically, unless you are running from the commandline
if isinteractive() == false
  main(parse_commandline())
end