using Coverage
using Printf

cd(joinpath(@__DIR__, "..")) do
    covered_lines, total_lines = get_summary(process_folder())
    percentage = covered_lines / total_lines * 100
    @printf "Coverage: %.1f%% covered\n" percentage
end
