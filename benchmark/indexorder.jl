using BenchmarkTools

suite = BenchmarkGroup()

suite["spherical"] = BenchmarkGroup()

begin
    resol_r = 10_000
    nfields = 1

    a_long_first = rand(resol_r, nfields)
    r_long_first = range(1.0, 10.0, resol_r)
    suite["spherical"]["r first"] = @benchmarkable abs2.($a_long_first) ./ $r_long_first

    a_long_second = rand(nfields, resol_r)
    r_long_second = reshape(r_long_first, 1, resol_r)
    suite["spherical"]["field first"] = @benchmarkable abs2.($a_long_second) ./
        $r_long_second
end

results = run(suite)
@show judge(minimum(results["spherical"]["r first"]), minimum(results["spherical"]["field first"]))
