# This script converts input MAT files into a JSON equivalent that we can import
# into CasADi.

import JSON
import PowerModels

function _convert_to_json(filename::String)
    data = PowerModels.parse_file(filename)
    PowerModels.standardize_cost_terms!(data, order = 2)
    PowerModels.calc_thermal_limits!(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    for branch in values(ref[:branch])
        branch["g"], branch["b"] = PowerModels.calc_branch_y(branch)
        branch["tr"], branch["ti"] = PowerModels.calc_branch_t(branch)
    end
    open(replace(filename, ".m" => ".json"), "w") do io
        write(io, JSON.json(ref))
    end
    return
end

data_directory = joinpath(dirname(@__DIR__), "data")
for file in filter(f -> endswith(f, ".m"), readdir(data_directory))
    _convert_to_json(joinpath(data_directory, file))
end
