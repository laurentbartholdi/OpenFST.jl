using CxxWrap, Libdl

openfstver = "openfst-1.8.2"
openfstdir = "$(@__DIR__)/$openfstver/"

run(`tar xfz $(@__DIR__)/$openfstver.tar.gz -C $(@__DIR__)`)
run(pipeline(openfstver*".patch",`patch -p0`))
cd(openfstdir) do
    run(`./configure --prefix=$(@__DIR__)`)
    run(`make -j install`)
end
# delete!(ENV,"LDFLAGS")
# delete!(ENV,"CPPFLAGS")

run(`cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$(CxxWrap.prefix_path()) -DJulia_EXECUTABLE=$(joinpath(Sys.BINDIR,"julia")) .`)
run(`cmake --build . --config Release`)
