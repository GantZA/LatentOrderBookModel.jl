include("summary_stats.jl")


# Method of Simulated Moments Basic
# Used to construct an objective function that measures errors relating to
# mean, kurtosis and Kolmogorov-Smirnov (KS) test and generalized Hurst Exponent
# when comparing a log price time series measure from the data and a simualted
# log price time series.

# Minimize the Objective Function using Nelder-Mead

function weight_matrix(log_prices, b, num_bootstap)
    # Step 1: Apply a Moving Block Bootstrap to the Measured Series
    n = size(log_prices, 1)

    b_samples = zeros(n, num_bootstap)
    block_ind = 1:n-b+1

    for i = 1:num_bootstap
        rand_blocks = sample(block_ind, Int(n/b))
        sample_ind = transpose(repeat(rand_blocks, 1, b))
        sample_ind = sample_ind[:]
        addition_vec = repeat(0:b-1,Int(n/b))
        sample_ind = sample_ind + addition_vec

        b_samples[:,i] = log_prices[sample_ind]
    end
    # Step 2: Calculate Distributions for Each Moment and Test Statistic
    dist = get_summary_stats(b_samples, log_prices)
    W = inv(cov(dist))
    return W
end



function select_moments(log_prices)
    mean_log = mean(log_prices)
    std_log = std(log_prices)
    kurt_log = normal_kurtosis(log_prices)
    ks_stat = 0.0
    hurst_log = generalized_hurts_exp(log_prices)


    return [mean_log, std_log, kurt_log, ks_stat, hurst_log]
end

function objective_value(simualted_log_prices, log_prices, b, num_bootstap)

    log_price_stats = select_moments(log_prices)

    num_replications = size(simualted_log_prices, 2)
    G = zeros(num_replications, 5)

    summary_stats_mat = get_summary_stats(simualted_log_prices, log_prices)
    print(size(summary_stats_mat))
    G = summary_stats_mat - transpose(repeat(log_price_stats,1,num_replications))
    G = transpose(mean(G, dims = 1))
    W = weight_matrix(log_prices, b, num_bootstap)
    return (transpose(G) * W * G)[1,1]
end
