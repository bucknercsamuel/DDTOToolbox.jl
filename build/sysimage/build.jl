const build_dir = @__DIR__
PackageCompiler.create_sysimage(
    ["DDTOSCP"];
    sysimage_path=joinpath(build_dir, "ddtoscp.so"),
    precompile_execution_file=[joinpath(build_dir, "precompile_execution.jl")],
    # incremental=false,
    # filter_stdlibs=false
)