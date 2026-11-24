using Test
using WittBifrostEOS
const W = WittBifrostEOS

@testset "WittBifrostEOS" begin
    @testset "shim" begin
        # the committed shim is what deps/witt/build.jl generates from the committed source
        dir = mktempdir()
        withenv("WITTBIFROSTEOS_BUILD_DIR" => dir) do
            include(joinpath(pkgdir(W), "deps", "witt", "build.jl")); Base.invokelatest(getfield(Main, :make_shim))
        end
        @test read(joinpath(dir, "witt_shim.cpp"), String) == read(joinpath(pkgdir(W), "deps", "witt", "witt_shim.cpp"), String)
        @test isfile(libwitt_path())
        @test use_libwitt!(libwitt_path()) == libwitt_path()
        @test witt_fix_gasc_h2()
        lib = libwitt_handle()
        for s in (:witt_create_backend, :witt_continuum_state, :witt_kurucz_populations, :witt_kurucz_populations_given_ne,
                  :witt_given_state_parallel, :witt_thermodynamics_from_pgas_parallel)           # FFNOInversion's entry points
            @test W.Libdl.dlsym(lib, s; throw_error = false) !== nothing
        end
        other = joinpath(mktempdir(), "libwitt.so"); cp(libwitt_path(), other)
        @test_throws InvalidInput use_libwitt!(other)                                       # one library per process
        @test use_libwitt!(other; strict = false) == libwitt_path() != other
    end

    m = EOSModel(default_pf_file(); prec = 1e-10, maxit = 5000)
    nlte = W.with_nlte(W.with_switch(m, 15019.647166203165, "robust_witt"), load_nlte_correction(default_nlte_correction()))
    lte = W.with_switch(m, 15019.647166203165, "robust_witt")
    @testset "model" begin
        s, code = witt_state(m.ctxs[1], 5778.0, 1.0e4)
        @test code == 0 && 1e-5 < s.rho < 1e-3 && s.ne > 0
        # the corona is LTE in both (the correction is 0 above 10^4.6 K); the chromosphere departs
        for (T, P) in ((3.0e6, 0.1), (1.0e5, 1.0))
            a, _, oka = model_fields(lte, log(T), log(P)); b, _, okb = model_fields(nlte, log(T), log(P))
            @test oka && okb && all(collect(a) .≈ collect(b))
        end
        a, _, _ = model_fields(lte, log(6500.0), log(0.1)); b, fl, _ = model_fields(nlte, log(6500.0), log(0.1))
        @test fl & EOS_NLTE != 0 && abs(a[2] - b[2]) > 0.05
    end

    @testset "table" begin
        env = Envelope(4000.0, 2.0e4, -1.0, 4.0, "test")
        tab = build_table(nlte, env; dlog10T = 0.02, dlog10P = 0.2, refine_max_levels = 0, log = devnull)
        path = joinpath(mktempdir(), "t.h5")
        meta = Dict{String,Any}("eos_model" => "nlte", "pf_file" => default_pf_file(), "pf_sha256" => file_sha256(default_pf_file()),
            "prec" => 1e-10, "maxit" => 5000, "abundances" => "wittmann_default", "T_sw" => 15019.647166203165, "hot_branch" => "robust_witt",
            "nlte_correction" => default_nlte_correction(), "nlte_correction_sha256" => file_sha256(default_nlte_correction()))
        save_table(path, tab, meta)
        t = load_table(path; model = "nlte")
        @test eos_model(t) == "nlte" && t.X == tab.X
        @test_throws InvalidInput load_table(path; model = "lte")
        X, Y = log(6000.0), log(10.0)
        r, n, c, st = state(t, X, Y); v, _, _ = model_fields(nlte, X, Y)
        @test isapprox(r, v[1]; atol = 1e-3) && isapprox(n, v[2]; atol = 1e-3)
        @test all(isfinite, hydrogen(t, X, Y)[1:2])
        # outside the table: the model it was built from
        e = ExtendedEOS(t)
        @test collect(extended_fields(e, X, Y))[1:2] == [r, n]
        X2 = log(1.0e6)
        @test !in_table(t, X2, Y)
        @test collect(extended_fields(e, X2, Y)) == collect(model_fields(nlte, X2, Y)[1])

        # state_d: the values of state and the X, Y derivatives of the same interpolant (cell interiors and edges)
        worst = 0.0
        for X in range(log(4200.0), log(1.9e4); length = 23), Y in range(-0.9, 3.9; length = 17)
            d = state_d(t, X, Y)
            @test collect(d[1:3]) == collect(state(t, X, Y)[1:3]) && d[10] == state(t, X, Y)[4]
            h = 1e-6
            for f in 1:3
                fdX = (state(t, X + h, Y)[f] - state(t, X - h, Y)[f]) / 2h
                fdY = (state(t, X, Y + h)[f] - state(t, X, Y - h)[f]) / 2h
                worst = max(worst, abs(fdX - d[2+2f]) / max(abs(fdX), 1e-2), abs(fdY - d[3+2f]) / max(abs(fdY), 1e-2))
            end
        end
        @test worst < 1e-6
        @test all(isnan, state_d(t, X2, Y)[1:9]) && state_d(t, X2, Y)[10] == EOS_OUT_OF_ENVELOPE
        X3, Y3 = tab.X[7], tab.Y[5]                          # on a node: one-sided derivative of the cell found
        @test all(isfinite, state_d(t, X3, Y3)[1:9])
    end
end
