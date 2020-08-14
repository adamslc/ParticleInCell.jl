using Coverage
using Printf

cd(joinpath(@__DIR__, "..")) do
    clean_folder("src")
    clean_folder("test")
end
