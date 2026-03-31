using BenchmarkTools

resol = 10_000
nfields = 1

a_long_first = rand(resol, nfields)
r_long_first = range(1.0, 10.0, resol)
b_long_first = @benchmark abs2.($a_long_first) ./ $r_long_first

a_long_second = rand(nfields, resol)
r_long_second = reshape(r_long_first, 1, resol)
b_long_second = @benchmark abs2.($a_long_second) ./ $r_long_second

@show judge(minimum(b_long_first), minimum(b_long_second))
