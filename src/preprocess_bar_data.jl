using CSV
using DataFrames

function process_bar_data(file_name)
    bar_data = CSV.read(file_name,
        header = [:price])

    bar_data.log_price = log.(bar_data.price)
    return bar_data
end
