using ArgParse

include("ORDGDP_pmap.jl")

# TODO... need to specify the solver and the algorithm at the command line
# TOOO... need to make this a Julia module... will make things easier in the long run



# general function for parsing the command line and putting in some default values and parameters
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--file", "-f"
            help = "an argument for the reslient design json file"
            default = "data/34Bus_Ice_Harden_Damageable_Rural_70Percent.json"
    end

    return parse_args(s)
end


# The main entry point into running a power flow solver
function main(parsed_args) 

  rdt_file = parsed_args["file"]
    
  solveORDGDP(rdt_file)  
    
 end

# this prevents main from being run automatically, unless you are running from the commandline
if isinteractive() == false
  main(parse_commandline())
end