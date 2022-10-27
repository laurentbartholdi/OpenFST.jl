using CxxWrap, Libdl

openfstver = "openfst-1.8.2"
openfstdir = "$(@__DIR__)/$openfstver/"

# run(`tar xfz $(@__DIR__)/$openfstver.tar.gz -C $(@__DIR__)`)
cd(()->run(`./configure --prefix=$(@__DIR__)`), openfstdir)
cd(()->run(`make -j install`), openfstdir)

# mkdir build && cd build
run(`cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$(CxxWrap.prefix_path()) .`)
run(`cmake --build . --config Release`)
