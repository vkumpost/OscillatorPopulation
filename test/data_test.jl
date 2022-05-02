@testset "load_biolum" begin

    # Load the available plates
    df = load_biolum("Plate U1")
    @test names(df) == ["Time", "Light", "WT A (1)", "WT A (2)", "WT A (3)",
        "WT A (4)", "WT A (5)", "WT A (6)", "WT A (7)", "WT A (8)", "WT B (1)",
        "WT B (2)", "WT B (3)", "WT B (4)", "WT B (5)", "WT B (6)", "WT B (7)",
        "WT B (8)"
    ]
    @test df[1, "Time"] ≈ 6.16
    @test df[end, "Time"] ≈ 499.89
    @test length(df[!, "Time"]) == 510

    df = load_biolum("Plate U2")
    @test names(df) == ["Time", "Light", "DAP 49 (1)", "DAP 49 (2)",
        "DAP 49 (3)", "DAP 49 (4)"
    ]
    @test df[1, "Time"] ≈ 8.28
    @test df[end, "Time"] ≈ 371.87
    @test length(df[!, "Time"]) == 835

    df = load_biolum("Plate U3")
    @test names(df) == ["Time", "Light", "Pool 1 (1)", "Pool 1 (2)",
        "Pool 1 (3)", "Pool 1 (4)", "Pool 1 (5)", "Pool 1 (6)", "Pool 1 (7)",
        "Pool 1 (8)", "Pool 2 (1)", "Pool 2 (2)", "Pool 2 (3)", "Pool 2 (4)",
        "Pool 2 (5)", "Pool 2 (6)", "Pool 2 (7)", "Pool 2 (8)", "Pool 3 (1)",
        "Pool 3 (2)", "Pool 3 (3)", "Pool 3 (4)", "Pool 3 (5)", "Pool 3 (6)",
        "Pool 3 (7)", "Pool 3 (8)", "Pool 4 (1)", "Pool 4 (2)", "Pool 4 (3)",
        "Pool 4 (4)", "Pool 4 (5)", "Pool 4 (6)", "Pool 4 (7)", "Pool 4 (8)"
    ]
    @test df[1, "Time"] ≈ 5.8
    @test df[end, "Time"] ≈ 201.3
    @test length(df[!, "Time"]) == 171

    df = load_biolum("Plate U4")
    @test names(df) == ["Time", "Light", "Pool 1 (1)", "Pool 1 (2)",
        "Pool 1 (3)", "Pool 1 (4)", "Pool 1 (5)", "Pool 1 (6)", "Pool 1 (7)",
        "Pool 1 (8)", "Pool 2 (1)", "Pool 2 (2)", "Pool 2 (3)", "Pool 2 (4)",
        "Pool 2 (5)", "Pool 2 (6)", "Pool 2 (7)", "Pool 2 (8)"
    ]
    @test df[1, "Time"] ≈ 4.92
    @test df[end, "Time"] ≈ 146.17
    @test length(df[!, "Time"]) == 126

    df = load_biolum("Plate D1 A")
    @test names(df) == ["Time", "Light", "Forskolin 5 uM (1)",
        "Forskolin 5 uM (2)", "Forskolin 5 uM (3)", "Forskolin 5 uM (4)",
        "Forskolin 10 uM (1)", "Forskolin 10 uM (2)", "Forskolin 10 uM (3)",
        "Forskolin 10 uM (4)", "Forskolin 15 uM (1)", "Forskolin 15 uM (2)",
        "Forskolin 15 uM (3)", "Forskolin 15 uM (4)", "DBC 0.5 mM (1)",
        "DBC 0.5 mM (2)", "DBC 0.5 mM (3)", "DBC 0.5 mM (4)", "DBC 1 mM (1)",
        "DBC 1 mM (2)", "DBC 1 mM (3)", "DBC 1 mM (4)", "DBC 3 mM (1)",
        "DBC 3 mM (2)", "DBC 3 mM (3)", "DBC 3 mM (4)", "U0126 10 uM (1)",
        "U0126 10 uM (2)", "U0126 10 uM (3)", "U0126 10 uM (4)",
        "U0126 20 uM (1)", "U0126 20 uM (2)", "U0126 20 uM (3)",
        "U0126 20 uM (4)", "U0126 40 uM (1)", "U0126 40 uM (2)",
        "U0126 40 uM (3)", "U0126 40 uM (4)", "Control DMSO (1)",
        "Control DMSO (2)", "Control DMSO (3)", "Control DMSO (4)"
    ]
    @test df[1, "Time"] ≈ 12.18
    @test df[end, "Time"] ≈ 105.91
    @test length(df[!, "Time"]) == 92

    df = load_biolum("Plate D2 A")
    @test names(df) == ["Time", "Light", "Forskolin 5 uM (1)",
        "Forskolin 5 uM (2)", "Forskolin 5 uM (3)", "Forskolin 5 uM (4)",
        "Forskolin 10 uM (1)", "Forskolin 10 uM (2)", "Forskolin 10 uM (3)",
        "Forskolin 10 uM (4)", "Forskolin 15 uM (1)", "Forskolin 15 uM (2)",
        "Forskolin 15 uM (3)", "Forskolin 15 uM (4)", "DBC 0.5 mM (1)",
        "DBC 0.5 mM (2)", "DBC 0.5 mM (3)", "DBC 0.5 mM (4)", "DBC 1 mM (1)",
        "DBC 1 mM (2)", "DBC 1 mM (3)", "DBC 1 mM (4)", "DBC 3 mM (1)",
        "DBC 3 mM (2)", "DBC 3 mM (3)", "DBC 3 mM (4)", "U0126 10 uM (1)",
        "U0126 10 uM (2)", "U0126 10 uM (3)", "U0126 10 uM (4)",
        "U0126 20 uM (1)", "U0126 20 uM (2)", "U0126 20 uM (3)",
        "U0126 20 uM (4)", "U0126 40 uM (1)", "U0126 40 uM (2)",
        "U0126 40 uM (3)", "U0126 40 uM (4)", "Control DMSO (1)",
        "Control DMSO (2)", "Control DMSO (3)", "Control DMSO (4)"
    ]
    @test df[1, "Time"] ≈ 11.15
    @test df[end, "Time"] ≈ 173.07
    @test length(df[!, "Time"]) == 174

    df = load_biolum("Plate D1 B")
    @test names(df) == ["Time", "Light", "PMA 0.5 uM (1)", "PMA 0.5 uM (2)", 
        "PMA 0.5 uM (3)", "PMA 0.5 uM (4)", "PMA 1 uM (1)", "PMA 1 uM (2)",
        "PMA 1 uM (3)", "PMA 1 uM (4)", "PMA 3 uM (1)", "PMA 3 uM (2)",
        "PMA 3 uM (3)", "PMA 3 uM (4)", "EGF 30 ng/ml (1)", "EGF 30 ng/ml (2)",
        "EGF 30 ng/ml (3)", "EGF 30 ng/ml (4)", "EGF 50 ng/ml (1)",
        "EGF 50 ng/ml (2)", "EGF 50 ng/ml (3)", "EGF 50 ng/ml (4)",
        "EGF 80 ng/ml (1)", "EGF 80 ng/ml (2)", "EGF 80 ng/ml (3)",
        "EGF 80 ng/ml (4)", "Ro-318220 2 uM (1)", "Ro-318220 2 uM (2)",
        "Ro-318220 2 uM (3)", "Ro-318220 2 uM (4)", "Ro-318220 5 uM (1)",
        "Ro-318220 5 uM (2)", "Ro-318220 5 uM (3)", "Ro-318220 5 uM (4)",
        "Ro-318220 8 uM (1)", "Ro-318220 8 uM (2)", "Ro-318220 8 uM (3)",
        "Ro-318220 8 uM (4)", "Control DMSO (1)", "Control DMSO (2)",
        "Control DMSO (3)", "Control DMSO (4)"
    ]
    @test df[1, "Time"] ≈ 11.33
    @test df[end, "Time"] ≈ 97.85
    @test length(df[!, "Time"]) == 85

    df = load_biolum("Plate D2 B")
    @test names(df) == ["Time", "Light", "PMA 0.5 uM (1)", "PMA 0.5 uM (2)", 
        "PMA 0.5 uM (3)", "PMA 0.5 uM (4)", "PMA 1 uM (1)", "PMA 1 uM (2)",
        "PMA 1 uM (3)", "PMA 1 uM (4)", "PMA 3 uM (1)", "PMA 3 uM (2)",
        "PMA 3 uM (3)", "PMA 3 uM (4)", "EGF 30 ng/ml (1)", "EGF 30 ng/ml (2)",
        "EGF 30 ng/ml (3)", "EGF 30 ng/ml (4)", "EGF 50 ng/ml (1)",
        "EGF 50 ng/ml (2)", "EGF 50 ng/ml (3)", "EGF 50 ng/ml (4)",
        "EGF 80 ng/ml (1)", "EGF 80 ng/ml (2)", "EGF 80 ng/ml (3)",
        "EGF 80 ng/ml (4)", "Ro-318220 2 uM (1)", "Ro-318220 2 uM (2)",
        "Ro-318220 2 uM (3)", "Ro-318220 2 uM (4)", "Ro-318220 5 uM (1)",
        "Ro-318220 5 uM (2)", "Ro-318220 5 uM (3)", "Ro-318220 5 uM (4)",
        "Ro-318220 8 uM (1)", "Ro-318220 8 uM (2)", "Ro-318220 8 uM (3)",
        "Ro-318220 8 uM (4)", "Control DMSO (1)", "Control DMSO (2)",
        "Control DMSO (3)", "Control DMSO (4)"
    ]
    @test df[1, "Time"] ≈ 10.05
    @test df[end, "Time"] ≈ 151.5
    @test length(df[!, "Time"]) == 202

    # Select specific event-bound times from the plates
    df = load_biolum("Plate U2"; first_event=2)
    @test names(df) == ["Time", "Light", "DAP 49 (1)", "DAP 49 (2)",
        "DAP 49 (3)", "DAP 49 (4)"
    ]
    t_start = (23.33 + 23.76) / 2
    @test all(Vector(df[1, :]) .≈ [23.76 - t_start, 1, 8838, 8635, 7313, 8540])
    @test all(Vector(df[end, :]) .≈ [371.87 - t_start, 1, 5291, 5178, 977, 5291])
    @test length(df[!, "Time"]) == 799

    df = load_biolum("Plate U2"; last_event=8)
    @test names(df) == ["Time", "Light", "DAP 49 (1)", "DAP 49 (2)",
        "DAP 49 (3)", "DAP 49 (4)"
    ]
    @test all(Vector(df[1, :]) .≈ [8.28, 1, 10400, 11070, 6800, 9030])
    @test all(Vector(df[end, :]) .≈ [266.52, 0, 5061, 5097, 2614, 5248])
    @test length(df[!, "Time"]) == 590

    df = load_biolum("Plate U2"; first_event=2, last_event=8)
    @test names(df) == ["Time", "Light", "DAP 49 (1)", "DAP 49 (2)",
        "DAP 49 (3)", "DAP 49 (4)"
    ]
    t_start = (23.33 + 23.76) / 2
    @test all(Vector(df[1, :]) .≈ [23.76 - t_start, 1, 8838, 8635, 7313, 8540])
    @test all(Vector(df[end, :]) .≈ [266.52 - t_start, 0, 5061, 5097, 2614, 5248])
    @test length(df[!, "Time"]) == 554

    df = load_biolum("Plate U2"; max_time=254)
    @test names(df) == ["Time", "Light", "DAP 49 (1)", "DAP 49 (2)",
        "DAP 49 (3)", "DAP 49 (4)"
    ]
    @test all(Vector(df[1, :]) .≈ [8.28, 1, 10400, 11070, 6800, 9030])
    @test all(Vector(df[end, :]) .≈ [253.62, 1, 4442, 4378, 2387, 4652])
    @test length(df[!, "Time"]) == 560

    df = load_biolum("Plate U2"; first_event=2, max_time=254)
    @test names(df) == ["Time", "Light", "DAP 49 (1)", "DAP 49 (2)",
        "DAP 49 (3)", "DAP 49 (4)"
    ]
    @test all(Vector(df[1, :]) .≈ [23.76 - t_start, 1, 8838, 8635, 7313, 8540])  # row 38
    @test all(Vector(df[end, :]) .≈ [277.27 - t_start, 1, 4477, 4552, 2443, 4701])  # row 616
    @test length(df[!, "Time"]) == 579

end
