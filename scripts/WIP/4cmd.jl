using DrWatson
@quickactivate :SCI499
using ArgParse

#create parse command line function
function parse_commandline()
# Create a new ArgParseSettings object to take in command line arguments
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--mixing_matrix_parameter"
            arg_type = Array
            default = 1.0
        "--beta"
            arg_type = Float64
            default = 0.5
        "--gamma"
            arg_type = Float64
            default = 0.1
        "--pop_info_file"
            arg_type = String
            default = "data/vic_pop2021.csv"
    end
     # Parse the command line arguments

    return parse_args(s)
end

function main()
    # Parse the command line arguments
    args = parse_commandline()
    # Print the arguments
    println(args)

end

main()
    