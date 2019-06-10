using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "SEED"
            help = "Seed for randomness"
            required = true
            arg_type = Int
            default = 1
        "T"
            help = "Number of time periods"
            required = true
            arg_type = Int
            default = 2300
        "τ"
            required = true
            arg_type = Int
            default = 10
        "initial_mid_price"
            required = true
            arg_type = Float64
            default = 238.745
        "n_spatial_points"
            required = true
            arg_type = Int
            default = 501
        "boltz_const"
            required = true
            arg_type = Float64
            default = 2
        "sample_std"
            required = true
            arg_type = Float64
            default = 7.415
        "σ"
            required = true
            arg_type = Float64
            default = 0.001
        "D"
            required = true
            arg_type = Float64
            default = 5.0
        "η"
            required = true
            arg_type = Float64
            default = 0.001
        "λ"
            required = true
            arg_type = Float6
            default = 1.0
        "μ"
            required = true
            arg_type = Float64
            default = 0.5
    end

    return parse_args(s)
end
