using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "SEED"
            help = "Seed for randomness"
            arg_type = Int
            default = 1
        "T"
            help = "Number of time periods"
            arg_type = Int
            default = 100
        "τ"
            arg_type = Int
            default = 10
        "m"
            arg_type = Float64
            default = 100.0
        "m"
            arg_type = Int
            default = 101
        "boltz_const"
            arg_type = Float64
            default = 1.0
        "sample_std"
            arg_type = Float64
            default = 4.0
        "D"
            arg_type = Float64
            default = 5.0
        "nu"
            arg_type = Float64
            default = 0.001
        "α"
            arg_type = Float64
            default = 1.0
        "λ"
            arg_type = Float64
            default = 1.0
        "μ"
            arg_type = Float64
            default = 0.5
    end

    return parse_args(s)
end
