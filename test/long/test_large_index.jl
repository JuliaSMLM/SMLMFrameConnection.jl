using Test, SMLMFrameConnection, SMLMData

# Above 2^24 = 16,777,216 localizations, a Float32 index column rounded odd indices to
# their even neighbors: track labels were written to the wrong localizations and some
# localizations were never addressed. 16,779,216 localizations in pairs that share a
# track_id, through the grouping step (about 90 s and 6 GB).
@testset "cluster indices stay exact above 2^24 localizations" begin
    n = 2^24 + 2000
    em = Vector{Emitter2DFit{Float64}}(undef, n)
    for i in 1:n
        t = (i + 1) ÷ 2
        em[i] = Emitter2DFit{Float64}(1.0 + 1e-6t, 1.0, 1000.0, 10.0, 0.01, 0.01, 0.0,
                                      30.0, 3.0, isodd(i) ? 1 : 2, 1, t, i)
    end
    smld = BasicSMLD(em, IdealCamera(1:512, 1:512, 0.1), 2, 1, Dict{String,Any}())
    clusters = SMLMFrameConnection.organizeclusters(smld)
    @test length(clusters) == n ÷ 2
    addressed = zeros(Int, n)
    wrong = 0
    for c in clusters
        tid = Int(c[1, 8])
        for j in Int.(c[:, 9])
            addressed[j] += 1
            em[j].track_id == tid || (wrong += 1)
        end
    end
    @test wrong == 0
    @test all(==(1), addressed)
end
