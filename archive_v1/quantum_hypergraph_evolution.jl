# ============================================================================
# PHASE 7.6: PERFECT FULL-WINDOW RECURSIVE HYPERGRAPH EVOLUTION
# MÔ HÌNH NÂNG CẤP ĐỘNG: SỬA LỖI TRÀN MÀN HÌNH - TỰ ĐỘNG DỆT LƯỚI WOLFRAM
# AUTHOR: Outlier Bui Huy Khoi (Age 15)
# ============================================================================

using GLMakie

function run_hypergraph_evolution_perfect()
    println("=========================================================")
    println("ĐANG CẤU HÌNH LẠI Ô LƯỚI ĐỒ HỌA OPENGL TRÀN MÀN HÌNH...")
    println("Hệ thống đang đồng bộ tần số quét với card màn hình...")
    println("=========================================================")

    # 1. THIẾT LẬP KÍCH THƯỚC VÀ ÉP GRID TRÀN TOÀN CỬA SỔ
    fig = Figure(size = (1200, 900), backgroundcolor = :black)
    
    # Ép trục Axis3 ghim chặt vào ô lưới trung tâm [1, 1] của Fig để không bị co nhỏ
    ax = Axis3(fig[1, 1], 
               title = "SIÊU ĐỒ THỊ TIẾN HÓA FRACTAL: CHUỖI KỲ DỊ LIÊN TỤC (1+1=1)",
               titlecolor = :white, titlesize = 22, titlegap = 15,
               xlabel = "Trục X", ylabel = "Trục Y", zlabel = "Tầng Tiến Hóa Z",
               xlabelcolor = :white, ylabelcolor = :white, zlabelcolor = :white,
               xgridcolor = :gray15, ygridcolor = :gray15, zgridcolor = :gray15,
               backgroundcolor = :black)

    # ĐỊNH NGHĨA CHUỖI 5 TÂM ĐIỂM KỲ DỊ GỐC WOLFRAM
    z_centers = [1.5, 3.5, 5.5, 7.5, 9.5]
    t_axis = 0.0:0.04:8.0
    t_fusion = 4.0

    # 2. VÒNG LẶP DỆT MẠNG LƯỚI FRACTAL ĐA TẦNG VÀ TỰ SINH HẠT NHÂN
    for (node_idx, z_c) in enumerate(z_centers)
        x_yang, y_yang, z_yang = Float64[], Float64[], Float64[]
        x_yin, y_yin, z_yin = Float64[], Float64[], Float64[]
        
        for t in t_axis
            # Quy luật thu hẹp bán kính phễu hướng tâm phi tuyến tính
            radius = max(0.02, 2.0 - t * 0.45)
            theta = 8.0 * t  # Tốc độ xoáy cuộn tăng cao để bện xoắn mịn hơn
            
            if t < t_fusion
                # GIAI ĐOẠN MA SÁT NHỊ NGUYÊN (1+1=2)
                push!(x_yang, radius * sin(theta))
                push!(y_yang, radius * cos(theta))
                push!(z_yang, z_c + 0.5 + sin(5*t)*0.1)
                
                push!(x_yin, -radius * sin(theta))
                push!(y_yin, -radius * cos(theta))
                push!(z_yin, z_c - 0.5 - cos(5*t)*0.1)
            else
                # GIAI ĐOẠN KHÓA PHA HOÀN TOÀN (1+1=1)
                push!(x_yang, 0.0)
                push!(y_yang, 0.0)
                push!(z_yang, z_c)
                
                push!(x_yin, 0.0)
                push!(y_yin, 0.0)
                push!(z_yin, z_c)
            end
        end
        
        # Vẽ các luồng sáng Âm-Dương bện xoắn mịn màng
        lines!(ax, x_yang, y_yang, z_yang, color = :orange, linewidth = 2.5, alpha = 0.5)
        lines!(ax, x_yin, y_yin, z_yin, color = :dodgerblue, linewidth = 2.5, alpha = 0.5)
        
        # ĐÓNG DẤU ĐIỂM KỲ DỊ SIÊU ĐẶC MÀU HỒNG CÁNH SEN ATOM
        meshscatter!(ax, [0.0], [0.0], [z_c], color = :magenta, markersize = 0.22,
                     label = node_idx == 1 ? "Điểm Kỳ Dị Hạt Nhân" : "")
        
        # --- NÂNG CẤP ĐỘNG: TỰ SINH MẠNG LƯỚI NÚT XUNG QUANH TÂM (WOLFRAM BUBBLES) ---
        # Tạo thêm 8 hạt vệ tinh nhỏ tự động mọc ra xung quanh mỗi Điểm Kỳ Dị để dệt lưới Không-Thời gian
        for angle in 0:(2*pi/8):(2*pi - 0.1)
            scatter_r = 0.4  # Bán kính mọc lưới Fractal
            meshscatter!(ax, [scatter_r * sin(angle)], [scatter_r * cos(angle)], [z_c], 
                         color = :cyan, markersize = 0.06, alpha = 0.6)
        end
    end

    # 3. TRỤC TRƯỜNG NHẤT KHÔNG-THỜI GIAN KẾT NỐI TOÀN HỆ THỐNG
    lines!(ax, [0.0, 0.0], [0.0, 0.0], [0.0, 11.0], color = :magenta, linewidth = 3, 
           linestyle = :dash, label = "Trục Không-Thời Gian")

    axislegend(ax, labelcolor = :white, backgroundcolor = :gray10, position = :rt)
    
    # 4. TỰ ĐỘNG XUẤT FILE ẢNH CHẤT LƯỢNG CAO PHIÊN BẢN NÂNG CẤP LỚN
    output_path = "wolfram_hypergraph_advanced.png"
    save(output_path, fig)
    
    # Đẩy giao diện hiển thị trực tiếp
    display(fig)
    println("=========================================================")
    println("XỬ LÝ TRÀN MÀN HÌNH THÀNH CÔNG RỰC RỠ!")
    println("File ảnh cấu trúc nâng cấp đã lưu vĩnh viễn: ", output_path)
    println("=========================================================")
end

run_hypergraph_evolution_perfect()
