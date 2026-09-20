
module QuantumSymbiosis

using GLMakie

export run_quantum_lever_simulation

function run_quantum_lever_simulation(z_max::Float64=11.0)
    # 1. Khởi tạo cửa sổ nền đen tuyệt đối
    fig = Figure(size = (1200, 800), backgroundcolor = :black, focus_on_show = true)
    
    # 2. Định hình lưới 3D và nhuộm đen toàn bộ khung xương viền
    ax = Axis3(fig, 
               title = "QUANTUM SYMBIOSIS: FRACTAL HYPERGRAPH (1+1=1)", 
               titlecolor = :white, titlesize = 20,
               xlabel = "X", ylabel = "Y", zlabel = "Z",
               xlabelcolor = :white, ylabelcolor = :white, zlabelcolor = :white,
               backgroundcolor = :black,
               xspinecolor_1 = :black, xspinecolor_2 = :black,
               yspinecolor_1 = :black, yspinecolor_2 = :black,
               zspinecolor_1 = :black, zspinecolor_2 = :black)
    
    z_centers = [1.5, 3.5, 5.5, 7.5, 9.5]
    t_axis = 0.0:0.04:8.0
    t_fusion = 4.0

    for z_c in z_centers
        z_c > z_max && break
        x_yang, y_yang, z_yang = Float64[], Float64[], Float64[]
        x_yin, y_yin, z_yin = Float64[], Float64[], Float64[]
        
        for t in t_axis
            radius = max(0.02, 2.0 - t * 0.45)
            theta = 8.0 * t
            if t < t_fusion
                push!(x_yang, radius * sin(theta)); push!(y_yang, radius * cos(theta)); push!(z_yang, z_c + 0.5)
                push!(x_yin, -radius * sin(theta)); push!(y_yin, -radius * cos(theta)); push!(z_yin, z_c - 0.5)
            else
                push!(x_yang, 0.0); push!(y_yang, 0.0); push!(z_yang, z_c)
                push!(x_yin, 0.0); push!(y_yin, 0.0); push!(z_yin, z_c)
            end
        end
        lines!(ax, x_yang, y_yang, z_yang, color = :orange, linewidth = 2.5)
        lines!(ax, x_yin, y_yin, z_yin, color = :dodgerblue, linewidth = 2.5)
        meshscatter!(ax, [0.0], [0.0], [z_c], color = :magenta, markersize = 0.18)
        for angle in 0:(2*pi/8):(2*pi - 0.1)
            meshscatter!(ax, [0.5 * sin(angle)], [0.5 * cos(angle)], [z_c], color = :cyan, markersize = 0.06)
        end
    end
    
    xlims!(ax, -3.0, 3.0)
    ylims!(ax, -3.0, 3.0)
    zlims!(ax, 0.0, z_max)
    
    # 3. Ép toàn không gian bao quanh hòa vào màn đêm
    fig.scene.backgroundcolor = :black
    
    display(fig)
    return fig
end

end # module
