const sysimage_path = "build/sysimage/ddtoscp.so"
const build_dir = @__DIR__
println("Creating sysimage in $target_name")

PackageCompiler.create_sysimage(
    ["DDTOSCP"];
    sysimage_path=sysimage_path,
    precompile_execution_file=[joinpath(build_dir, "precompile_execution.jl")],
    precompile_statements_file=[joinpath(build_dir, "precompile_statements.jl")],
    incremental=false,
    filter_stdlibs=false
)