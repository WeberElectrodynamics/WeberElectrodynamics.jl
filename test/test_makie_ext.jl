using CairoMakie

# animate_weber ends with display(fig). CairoMakie's default display on macOS
# opens a Preview window; in headless CI it may fail. Push a no-op display on top
# of the stack so display(fig) dispatches to it and does nothing. We still exercise
# the full construction path (buffer, observables, figure build).
struct _NullDisplay <: AbstractDisplay end
Base.display(::_NullDisplay, ::Any) = nothing

function _with_null_display(f)
    Base.Multimedia.pushdisplay(_NullDisplay())
    try
        return f()
    finally
        Base.Multimedia.popdisplay()
    end
end

function _try_animate(f)
    try
        _with_null_display(f)
        return true
    catch e
        msg = sprint(showerror, e)
        if occursin("display", lowercase(msg)) || occursin("screen", lowercase(msg))
            return true
        end
        rethrow()
    end
end

@testset "Makie extension" begin
    @testset "animate_weber streaming mode" begin
        prob = make_coulomb_like_problem(tspan = (0.0, 0.1), dt = 0.01)
        @test _try_animate() do
            WeberElectrodynamics.animate_weber(
                prob;
                buffer_size = 50,
                tail_length = 20,
                compute_batch = 1,
            )
        end
    end

    @testset "animate_weber replay mode" begin
        prob = make_coulomb_like_problem(tspan = (0.0, 0.1), dt = 0.01)
        sol = solve(prob)
        @test _try_animate() do
            WeberElectrodynamics.animate_weber(
                sol;
                buffer_size = 50,
                tail_length = 20,
                compute_batch = 1,
            )
        end
    end

    @testset "animate_weber argument validation" begin
        prob = make_coulomb_like_problem(tspan = (0.0, 0.1), dt = 0.01)
        sol = solve(prob)

        @test_throws AssertionError WeberElectrodynamics.animate_weber(prob; buffer_size = 0)
        @test_throws AssertionError WeberElectrodynamics.animate_weber(prob; tail_length = 0)
        @test_throws AssertionError WeberElectrodynamics.animate_weber(prob; compute_batch = 0)
        @test_throws AssertionError WeberElectrodynamics.animate_weber(
            prob; buffer_size = 10, tail_length = 20,
        )
        @test_throws AssertionError WeberElectrodynamics.animate_weber(sol; stride = 0)

        # 1D problem should reject animation
        sys1d = WeberSystem(2, 1)
        prob1d = WeberProblem(
            sys1d, (0.0, 0.1),
            [1.0, -1.0], [0.1, -0.1];
            masses = [1.0, 1.0], charges = [1.0, -1.0],
            c = 1e10, dt = 0.01,
        )
        @test_throws AssertionError WeberElectrodynamics.animate_weber(prob1d)
    end
end
