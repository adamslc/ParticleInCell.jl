@info "Loading packages"
using HDF5
using Makie
using MakieLayout

function create_pts(i, xlabel, ylabel, prefix)
    @debug "Reading data for dump $i"
    f = h5open("$(prefix)_$i.h5", "r")
    xvals = f[xlabel][:]
    yvals = f[ylabel][:]
    close(f)

    return [Point(x, y) for (x, y) in zip(xvals, yvals)]
end

function create_display(prefix)
    directory = prefix[1:findlast('/', prefix)-1]
    files = readdir(directory)
    num_dumps = 0
    for file in files
        if startswith(file, prefix[findlast('/', prefix)+1:end] * "_")
            num_dumps += 1
        end
    end

    @info "Creating scene" num_dumps
    scene, layout = layoutscene(10)
    display(scene)

    @info "Populating scene"
    ax1 = layout[1,1] = LAxis(scene, xgridvisible=false, ygridvisible=false)
    ax2 = layout[2,1] = LAxis(scene, xgridvisible=false, ygridvisible=false)
    ax3 = layout[3,1] = LAxis(scene, xgridvisible=false, ygridvisible=false)
    linkxaxes!(ax1, ax2, ax3)
    xlims!(ax1, 0, 1)

    slide = LSlider(scene, range=0:num_dumps-1, startvalue=0)
    slider_text = lift(i -> "Dump=$i", slide.value)
    text_obj = LText(scene, slider_text)
    controls = layout[end+1, :] = GridLayout(tellwidth=false)
    controls[:h] = [slide, text_obj]

    dist_pts = lift(i -> create_pts(i, "positions", "velocities", prefix), slide.value)
    electrons = scatter!(ax1, dist_pts, color=:blue, markersize=2px)

    f = h5open("$(prefix)_0.h5", "r")
    pos0 = f["positions"][:]
    vel0 = f["velocities"][:]
    close(f)
    dist0 = lines!(ax1, pos0, vel0, color=:black, linewidth=0.5, linestyle=:dash)

    layout[1, 2] = LLegend(scene, [electrons, dist0], ["Electrons", "Initial distribution"])

    rho_pts = lift(i -> create_pts(i, "grid_pos", "rho", prefix), slide.value)
    rp = lines!(ax2, rho_pts, color=:red, xautolimits=false)

    rho_hd_pts = lift(i -> create_pts(i, "grid_pos_hd", "rho_hd", prefix), slide.value)
    rp_hd = lines!(ax2, rho_hd_pts, color=:orange, xautolimits=false)

    layout[2, 2] = LLegend(scene, [rp, rp_hd], ["Rho", "HD Rho"])

    elec_pts = lift(i -> create_pts(i, "grid_pos", "elec", prefix), slide.value)
    ep = lines!(ax3, elec_pts, color=:red, xautolimits=false)

    layout[3, 2] = LLegend(scene, [ep], ["E"])

    f = h5open("$(prefix)_$(slide.value[]).h5", "r")
    grid_pos = f["grid_pos"][:]
    close(f)
    lims1 = lift(rect -> [rect.origin[2], rect.origin[2] + rect.widths[2]], ax1.limits)
    lims2 = lift(rect -> [rect.origin[2], rect.origin[2] + rect.widths[2]], ax2.limits)
    lims3 = lift(rect -> [rect.origin[2], rect.origin[2] + rect.widths[2]], ax3.limits)
    for g in grid_pos
        lines!(ax1, [g, g], lims1, color=:black, linewidth=0.5, linestyle=:dot, xautolimits=false, yautolimits=false)
        lines!(ax2, [g, g], lims2, color=:black, linewidth=0.5, linestyle=:dot, xautolimits=false, yautolimits=false)
        lines!(ax3, [g, g], lims3, color=:black, linewidth=0.5, linestyle=:dot, xautolimits=false, yautolimits=false)
    end

    on(scene.events.keyboardbuttons) do but
        ispressed(but, Keyboard.left)  && set_close_to!(slide, slide.value[]-1)
        ispressed(but, Keyboard.right) && set_close_to!(slide, slide.value[]+1)
    end

    return scene, layout
end

function ui(default_str)
    old_str = default_str
    scene = nothing
    layout = nothing
    while true
        print("> ")
        str = readline()

        if str == "exit" || str == "quit"
            println("Goodbye")
            exit()
        elseif str == "save"
            print("Enter filename: ")
            plot_str = readline()
            Makie.save(plot_str, scene)
        elseif str == ""
            try
                scene, layout = create_display(old_str)
            catch
                println("ERROR")
            end
        else
            old_str = str
            try
                scene, layout = create_display(str)
            catch
                println("ERROR")
            end
        end
    end

    return
end

ui("data/dump")
