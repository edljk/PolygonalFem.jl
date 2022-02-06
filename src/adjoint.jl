"""
    fill_vertices!(pvector, parray)
"""
function fill_vertices!(pvector, parray)
    dim = size(parray[1], 1)
    # check dimensions 
    if sum(length.(parray)) != length(pvector)
        eeror("pb of dimensions in fill_vertices $(sum(length.(parray))) != $(length(pvector)) ")
    end
    c = 1
    for k = 1:length(parray)
        for l = 1:dim # fill by column 
            for m = 1:size(parray[k], 1)
                parray[k][m, l] = pvector[c]
                c += 1
            end
        end
    end
    return nothing 
end
#-------------------------------------------------------------------------------
function prodadjoint()
end