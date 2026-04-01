const build_dir = @__DIR__
PackageCompiler.create_sysimage(
    ["DDTOToolbox"];
    sysimage_path=joinpath(build_dir, "ddtotoolbox.so"),
    precompile_execution_file=[joinpath(build_dir, "precompile_execution.jl")],
    # incremental=false,
    # filter_stdlibs=false
)
