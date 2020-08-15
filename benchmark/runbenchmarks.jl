using BenchmarkCI

cd(joinpath(@__DIR__, "..")) do
    BenchmarkCI.judge()
    BenchmarkCI.postjudge()
end
