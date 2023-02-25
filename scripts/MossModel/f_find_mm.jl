using Match

function find_mm(path, subdir::Bool)::Vector{String}
    let 
        files = []
        for entry in readdir(path, join = true)
            let f = entry
                if f|>isfile
                    ext = split(f, ".")[2]
                    @match ext begin
                        "csv" => push!(files, f)
                    end
                elseif f|>isdir && subdir
                    for sub_entry in readdir(f, join = true)
                        let sub_f = sub_entry
                            if sub_f|>isfile
                                ext = split(sub_f, ".")[2]
                                @match ext begin
                                    "csv" => push!(files, sub_f)
                                end
                            end
                        end
                    end
                end
            end
        end
        files|>sort!
        files
    end
end



find_mm(datadir("Moss2018"), true)