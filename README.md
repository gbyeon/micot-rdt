# micot-rdt
This the repository for MICOT's RDT capability. The main point of entry is the cli.jl file.

# Installation

1. Install Julia http://julialang.org/downloads/
2. Launch Julia
3. Install the following packages from within Julia
	A. Pkg.add("Optim")
	B. Pkg.add("JuMP")
	C. Pkg.add("Cbc")
	D. Pkg.add("Clp")
	E. Pkg.add("SCS")
	F. Pkg.add("Ipopt")
	G. Pkg.add("ECOS")
	H. Pkg.add("GLPK")
	I. Pkg.add("JSON")
	J. Pkg.add("CoinOptServices")
	K. Pkg.add("AmplNLWriter")
	L. Pkg.add("ArgParse")
	M. Pkg.add("FactCheck")
	N. Optionally install the bridges to the commerical packages if you have them (CPLEX, Gurobi, etc.--see http://jump.readthedocs.io/en/latest/installation.html)
4. Add the julia binary executible to your path variable
5. Test running the code from the command line by calling "julia cli.jl --file  data/34Bus_Ice_Harden_Damageable_Rural_70Percent.json (note that this will fail at the moment as the capability is being rewritten in Julia to support non commerical optimization packages)

# Input API
	{
	"$schema": "http://json-schema.org/draft-04/schema#",
	"type": "object",
	"title": "Resilient Design Tool (RDT) JSON schema",
	"description": "These schema describes the fields for the JSON schema used by the RDT Tool",
	"properties": {
	"critical_load_met": {
			"type": "number",
			"description": "This is a number between 0 and 1 indicates the percentage of critical load that must be met in each damage scenario",
			"default": 0.98
		},
		"total_load_met": {
			"type": "number",
			"description": "This is a number between 0 and 1 indicates the percentage of non-critical load that must be met in each damage scenario",
			"default": 0.5
		},
		"chance_constraint": {
			"type": "number",
			"description": "This is a number between 0 and 1 indicates percentage of damage scenarios where the critical and total load met constraints must be met",
			"default": 1.0
		},
		"phase_variation": {
			"type": "number",
			"description": "This is a number that controls the amount of phase imbalance may occur at transformers in the system.",
			"default": 0.15
		},
	"buses": {
			"type": "array",
			"description": "This block describes the necessary bus (node) data for running an instance of RDT.",
			"items": {
				"type": "object",
				"description": "Each entry here contains information about a single bus",
				"properties": {
					"id": {
						"type": "string",
						"description": "Unique identifier for the bus",
					},
					"min_voltage": {
						"type": "number",
						"description": "Minimum voltage level for the bus in p.u.",
						"default": 0.8
					},
					"max_voltage": {
						"type": "number",
						"description": "Maximum voltage level for the bus in p.u.",
						"default": 1.2
					},
					"ref_voltage": {
						"type": "number",
						"description": "Reference voltage for the bus",
						"default": 1.0
					},
					"has_phase": {
						"type": "array",
						"description": "Array indicating whether or not a bus has a phase",
						"items": {
							"type": "boolean",
							"description": "Entries of the array",
							"default": true
						}
					},
					"has_generator": {
						"type": "boolean",
						"description": "Flag indicating whether or not the generator has a bus.",
						"default": false
					}
				},
				"required": [
					"id",
				]
			}
		},

	"loads": {
			"type": "array",
			"description": "This block describes the loads of the power system.",
			"items": {
				"type": "object",
				"description": "Each entry of the array provides information on one load.  More than one load per bus is allowed.",
				"properties": {
					"id": {
						"type": "string",
						"description": "Unique identifier for the load",
					},
					"node_id": {
						"type": "string",
						"description": "Identifier of the bus the load is attached to",
					},
					"has_phase": {
						"type": "array",
						"description": "Array that indicates whether this load serves a phase or not",
						"items": {
							"type": "boolean",
							"description": "There should be three entries in this array",
						}
					},
					"max_real_phase": {
						"type": "array",
						"description": "Maximum or desired real load",
						"items": {
							"type": "number",
							"description": "A 3 entry array indicating the phase real load",
						}
					},
					"max_reactive_phase": {
						"type": "array",
						"description": "Maximum or desired reactive load",
						"items": {
							"type": "number",
							"description": "A 3 entry array indicating the phase reactive load.",
						}
					},
					"is_critical": {
						"type": "boolean",
						"description": "a flag indicating whether or not a load is critical or not.",
						"default": false
					}
				},
				"required": [
					"id",
					"node_id",
					"has_phase",
					"max_real_phase",
					"max_reactive_phase",
				]
			}
		},

	"lines": {
			"type": "array",
			"description": "This is an array of the lines in the power system model.",
			"items": {
				"type": "object",
				"description": "Each entry contains a block containing information about a single line",
				"properties": {
					"id": {
						"type": "string",
						"description": "Unique identifier for the line",
					},
					"node1_id": {
						"type": "string",
						"description": "Identifier of the first node the line is connected to",
					},
					"node2_id": {
						"type": "string",
						"description": "Identifier of the second the line is connected to",
					},
					"has_phase": {
						"type": "array",
						"description": "Vector of booleans corresponding to whether of not the line carries a phase",
						"items": {
							"type": "boolean",
							"description": "A vector of 3 boolean values for the phases",
						}
					},
					"capacity": {
						"type": "number",
						"description": "MVA capacity of the line",
						"default": 1e30
					},
					"length": {
						"type": "number",
						"description": "The length of the line. Units need to be consistent with the line code entries",
						"default": 1.0
					},
					"num_phases": {
						"type": "integer",
						"description": "The number of phases the line carries.",
						"default": 3
					},
					"is_transformer": {
						"type": "boolean",
						"description": "Flag indicating whether or not the line is a transformer",
						"default": false
					},
					"line_code": {
						"type": "string",
						"description": "Line code identifier to get additional information about the line",
					},
					"construction_cost": {
						"type": "number",
						"description": "The cost of constructing this line.",
						"default": 1e30
					},
					"harden_cost": {
						"type": "number",
						"description": "The cost of hardening this line.",
						"default": 1e30
					},
					"switch_cost": {
						"type": "number",
						"description": "The cost of building a switch.",
						"default": 1e30
					},
					"is_new": {
						"type": "boolean",
						"description": "Flag for whether or not this line is a new construction option.",
						"default": false
					},
					"can_harden": {
						"type": "boolean",
						"description": "Flag for whether or not this line can be hardened.",
						"default": false
					},
					"can_add_switch": {
						"type": "boolean",
						"description": "Flag for whether or not a switch can be built here",
						"default": false
					},
					"has_switch": {
						"type": "boolean",
						"description": "Flag for whether or not a switch already exists at the line",
						"default": false
					}
				},
				"required": [
					"id",
					"node1_id",
					"node2_id",
					"has_phase",
					"length",
					"num_phases",
					"line_code",
				]
			}
		}


	"lines_codes": {
			"type": "array",
			"description": "This is an array of the line codes of the model.  This used to compactly model aspects of a line that are common across lines.  For example, impedance values.",
			"items": {
				"type": "object",
				"description": "Each entry contains a block containing information about a single line code",
				"properties": {
					"line code: {
						"type": "string",
						"description": "Unique identifier for the line code",
					},
					"num_phases": {
						"type": "integer",
						"description": "Number of phases for the line code",
					},
					"rmatrix": {
						"type": "array",
						"description": "An array of the resistance terms for the line (per length). The entries correspond to the rows of the rmatrix",
						"items": {
							"type": "array",
							"description": "The first entry is phase A. The second entry is phase B. The third entry is phase C",
				"items" : {
				 "type": "number",
								 "description": "Column entries for the rmatrix"
				}
						}
					},
					"xmatrix": {
						"type": "array",
						"description": "An array of the reactance terms for the line (per length). The entries correspond to the rows of the rmatrix",
						"items": {
							"type": "array",
							"description": "The first entry is phase A. The second entry is phase B. The third entry is phase C",
				"items" : {
				 "type": "number",
								 "description": "Column entries for the rmatrix"

				}
						}
					}
				},
				"required": [
					"line_code",
					"num_phases",
					"rmatrix",
					"xmatrix"
				]
			}
		},


	"generators": {
			"type": "array",
			"description": "This block describes the generators in the model.",
			"items": {
				"type": "object",
				"description": "Each entry of this array contains information about one generator in the model.  More than one generator can exist at a bus",
				"properties": {
					"id": {
						"type": "string",
						"description": "Unique identifier for the generator",
					},
					"node_id": {
						"type": "string",
						"description": "identifier of the bus the generator is connected to",
					},
					"has_phase": {
						"type": "array",
						"description": "Array of flags indicating whether or not a generator has a phase or not",
						"items": {
							"type": "boolean",
							"default": true
						}
					},
					"max_real_phase": {
						"type": "array",
						"description": "Existing maximum real power output of the generator for each phase",
						"items": {
							"type": "number",
							"default": 1.7976931348623e+308
						}
					},
					"max_reactive_phase": {
						"type": "array",
						"description": "Existing maximum raective power output of the generator for each phase",
						"items": {
							"type": "number",
							"default": 1.7976931348623e+308
						}
					},
					"microgrid_cost": {
						"type": "number",
						"description": "Cost per MW capacity of building (additional) distributed generation, i.e. expanded capacity",
						"default": 1e30
					},
					"microgrid_fixed_cost": {
						"type": "number",
						"description": "One-time fixed cost for building (additional) distributed generation, i.e. expanded capacity",
						"default": 0
					},
					"max_microgrid": {
						"type": "number",
						"description": "Maximum additional capacity that can be built at this generator.",
						"default": 0
					},
					"is_new": {
						"type": "boolean",
						"description": "Flag indicating whether or not new generation can be built",
						"default": false
					}
				},
				"required": [
					"id",
					"node_id",
					"has_phase",
					"max_real_phase",
					"max_reactive_phase"
				]
			}
		}
		
		"scenarios": {
      "type": "array",
      "description": "This block contains information about the damage scenarios",
      "items": {
        "type": "object",
        "description": "Each entry in the array contains information about a single damage scenario",
        "properties": {
          "id": {
            "type": "string",
            "description": "A unique identifier for the scenario",
          },
          "hardened_disabled_lines": {
            "type": "array",
            "description": "A list of identifiers for lines that are damaged even after being hardened",
            "items": {}
          },
          "disabled_lines": {
            "type": "array",
            "description": "A list of identifiers for lines that are damaged if they are not hardened",
            "items": {}
          }
        },
        "required": [
          "id",
        ]
      }
    } 
		
		
		
	},
	"required": [
	"buses",
	"loads",
	"lines",
    "line_codes",
	"generators",
	"scenarios"
	]
	}

# MIP Solver Settings API

The RDT tool requires a mixed integer programming solver to perform some of the optimization steps.  This can be provided with a simple json input file that contains the name of the solver and any optional parameters available to the solver (see https://jump.readthedocs.io/en/latest/ for a discussion on possible parameters)

For example, consider this json string

{
	"mip_solver": "IPOPT",
	"tol": 1e-6,
	"print_level": 0
}

The mip_solver is defined as IPOPT.  The other two arguments are solver parameters specific to ipopt.  Note that Ipopt cannot actaully solve problems that have discrete problems, so this may not be the best example.

Valid MIP solvers include "CPLEX", "MOSEK", "BONMIN", "GUROBI", and "CBC"

This optional parameter is provided with the --mip_solver flag at the command lilne.  If no parameters are provided, CBC is used as the MIP solver.

# Algorithm Settings API

TBD

# Output API

TBD






