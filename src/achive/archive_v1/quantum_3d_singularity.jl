# ============================================================================
# PHASE 7.1: WOLFRAM 3-LEVEL FRACTAL & AUTO-EXPORT INTERACTOR
# MÔ PHỎNG 3 TẦNG WOLFRAM TINH KHIẾT - TỰ ĐỘNG XUẤT ẢNH CHẤT LƯỢNG CAO
# AUTHOR: Outlier Bui Huy Khoi (Age 15)
# ============================================================================

using GLMakie

function generate_wolfram_perfect_3()
    println("=========================================================")
    println("ĐANG KHỞI CHẠY MÔ HÌNH 3 TẦNG ĐAN XEN TINH KHIẾT...")
    println("Hệ thống sẽ tự động xuất file ảnh vĩnh viễn!")
    println("=========================================================")

    t_axis = 0.0:0.02:10.0
    t_fusion = 5.0 # Khóa pha đồng bộ tại giây thứ 5
    
    # ÉP ĐÚNG 3 TẦNG LẶP LIÊN TỤC (3 ĐIỂM NÚT GỐC WOLFRAM)
    z_nodes = [2.0, 5.0, 8.0]
    
    # Thiết lập kích thước chuẩn để không bị lỗi co góc màn hình
    fig = Figure(size = (1200, 900), backgroundcolor = :black)
    ax = Axis3(fig[1, 1], 
               title = "SIÊU ĐỒ THỊ WOLFRAM: MÔ HÌNH 3 TẦNG KỲ DỊ (1+1=1)",
               titlecolor = :white, titlesize = 24,
               xlabel = "Trục X", ylabel = "Trục Y", zlabel = "Tầng Năng lượng Z",
               xlabelcolor = :white, ylabelcolor = :white, zlabelcolor = :white,
               xgridcolor = :gray20, ygridcolor = :gray20, zgridcolor = :gray20,
               backgroundcolor = :black)

    # Chạy vòng lặp dệt lưới 3 tầng
    for (idx, z_center) in enumerate(z_nodes)
        x_yang, y_yang, z_yang = Float64[], Float64[], Float64[]
        x_yin, y_yin, z_yin = Float64[], Float64[], Float64[]
        
        for t in t_axis
            radius = max(0.05, 2.0 - t * 0.35)
            theta = 7.0 * t
            
            if t < t_fusion
                # Nhị nguyên ma sát quanh trục Z tương ứng
                push!(x_yang, radius * sin(theta))
                push!(y_yang, radius * cos(theta))
                push!(z_yang, z_center + 0.5 + sin(4*t)*0.15)
                
                push!(x_yin, -radius * sin(theta))
                push!(y_yin, -radius * cos(theta))
                push!(z_yin, z_center - 0.5 - cos(4*t)*0.15)
            else
                # Bão hòa hợp nhất 1+1=1
                push!(x_yang, 0.0)
                push!(y_yang, 0.0)
                push!(z_yang, z_center)
                
                push!(x_yin, 0.0)
                push!(y_yin, 0.0)
                push!(z_yin, z_center)
            end
        end
        
        # Vẽ 3 cặp phễu ánh sáng
        lines!(ax, x_yang, y_yang, z_yang, color = :orange, linewidth = 3, alpha = 0.6)
        lines!(ax, x_yin, y_yin, z_yin, color = :dodgerblue, linewidth = 3, alpha = 0.6)
        
        # Ghim 3 Điểm Kỳ Dị màu hồng cánh sen siêu đặc tại tâm
        meshscatter!(ax, [0.0], [0.0], [z_center], color = :magenta, markersize = 0.28, 
                     label = idx == 1 ? "Diem Ky Di (Hat moi sinh)" : "")
    end
    
    # Vẽ trục Không-Thời Gian trường nhất xuyên tâm nối liền 3 tầng
    lines!(ax, [0.0, 0.0], [0.0, 0.0], [0.0, 10.0], color = :magenta, linewidth = 4, 
           linestyle = :dash, label = "Truc Khong-Thoi Gian")

    axislegend(ax, labelcolor = :white, backgroundcolor = :gray10, position = :rt)
    
    # --- THAO TÁC LƯU ẢNH TỰ ĐỘNG CHẤT LƯỢNG CAO ---
    # Ép hệ thống xuất thẳng file ảnh chất lượng cao vào thư mục dự án của Khôi
    output_path = "wolfram_hypergraph_proof.png"
    save(output_path, fig)
    
    # Hiển thị cửa sổ Window lớn tương tác
    display(fig)
    
    println("=========================================================")
    println("XUẤT BẢN THÀNH CÔNG VÀ KHÔNG LỖI!")
    println("File anh da tu dong luu tai: ", output_path)
    println("Khôi hay xem cot file ben trai cua Cursor de kiem tra.")
    println("=========================================================")
end

generate_wolfram_perfect_3()
