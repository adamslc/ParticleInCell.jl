@info "Loading packages"
using HDF5
using Makie
using MakieLayout

function update_plot(scene, layout, i)
    @info "Reading data for dump $i"
    f = h5open("data/dump_$i.h5", "r")
    positions = f["positions"][:]
    velocities = f["velocities"][:]

    grid_pos = f["grid_pos"][:]
    rho = f["rho"][:]
    phi = f["phi"][:]
    elec = f["elec"][:]

    grid_pos_hd = f["grid_pos_hd"][:]
    rho_hd = f["rho_hd"][:]

    @info "Updating plots"
    ax1 = layout[1, 1] = LAxis(scene)
    xlims!(ax1, 0, 1)
    ylims!(ax1, -1.2e5, 1.2e5)
    electrons = scatter!(ax1, positions, velocities, color=:blue, markersize=2px)

    layout[1, 2] = LLegend(scene, [electrons], ["Electrons"])

    ax2 = layout[2, 1] = LAxis(scene)
    xlims!(ax2, 0, 1)
    ylims!(ax1, -0.001, 0.001)
    rp = plot!(ax2, grid_pos, rho, color=:red)
    rp_hd = plot!(ax2, grid_pos_hd, rho_hd, color=:orange)

    layout[2, 2] = LLegend(scene, [rp, rp_hd], ["Rho", "HD Rho"])

    return
end

function create_pts(i, xlabel, ylabel)
    @debug "Reading data for dump $i"
    f = h5open("data/dump_$i.h5", "r")
    xvals = f[xlabel][:]
    yvals = f[ylabel][:]
    close(f)

    return [Point(x, y) for (x, y) in zip(xvals, yvals)]
end

@info "Creating scene"
scene, layout = layoutscene(10, resolution=(1200, 700))
display(scene)

@info "Populating scene"
ax1 = layout[1,1] = LAxis(scene)
ax2 = layout[2,1] = LAxis(scene)
linkxaxes!(ax1, ax2)
xlims!(ax1, 0, 1)

slide = LSlider(scene, range=0:100, startvalue=0)

slider_text = lift(i -> "Dump=$i", slide.value)
text_obj = LText(scene, slider_text)

controls = layout[end+1, :] = GridLayout(tellwidth=false)
controls[:h] = [slide, text_obj]

dist_pts = lift(i -> create_pts(i, "positions", "velocities"), slide.value)
electrons = scatter!(ax1, dist_pts, color=:blue, markersize=2px)

layout[1, 2] = LLegend(scene, [electrons], ["Electrons"])

rho_pts = lift(i -> create_pts(i, "grid_pos", "rho"), slide.value)
rp = lines!(ax2, rho_pts, color=:red)

rho_hd_pts = lift(i -> create_pts(i, "grid_pos_hd", "rho_hd"), slide.value)
rp_hd = lines!(ax2, rho_hd_pts, color=:orange)

layout[2, 2] = LLegend(scene, [rp, rp_hd], ["Rho", "HD Rho"])

on(scene.events.keyboardbuttons) do but
    ispressed(but, Keyboard.left)  && set_close_to!(slide, slide.value[]-1)
    ispressed(but, Keyboard.right) && set_close_to!(slide, slide.value[]+1)
end

@info "Ready. Press [ENTER] to exit."
readline()
