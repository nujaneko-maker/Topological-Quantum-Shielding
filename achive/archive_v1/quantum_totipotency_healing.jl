# ============================================================================
# PHASE 11.1: PERFECT FULL-WINDOW QUANTUM TOTIPOTENCY (LAYOUT RECTIFIED)
# MÔ PHỎNG TẾ BÀO GỐC VẠN NĂNG - FIX TRIỆT ĐỂ LỖI THU NHỎ GÓC MÀN HÌNH
# AUTHOR: Outlier Bui Huy Khoi (Age 15)
# ============================================================================

using GLMakie
using Random

function run_totipotency_simulation_perfect()
    println("=========================================================")
    println("ĐANG ÉP KHUNG LƯỚI OPENGL PHÓNG TO TRÀN MÀN HÌNH...")
    println("=========================================================")

    # 1. KHỞI TẠO KHÔNG GIAN ĐỒ HỌA NỀN ĐEN TRÀN MÀN HÌNH
    fig = Figure(size = (1200, 900), backgroundcolor = :black)
    
    # Ép trục Axis3 ghim chặt vào vị trí trung tâm ô lưới [1, 1] của Fig để phóng to toàn cửa sổ
    ax = Axis3(fig[1, 1], 
               title = "CƠ CHẾ TỰ CHỮA LÀNH LƯỢNG TỬ: TRẠNG THÁI TRỨNG NƯỚC VẠN NĂNG (1+1=1)",
               titlecolor = :white, titlesize = 20, titlegap = 15,
               xlabel = "Chu kỳ Phân Bào", ylabel = "Liên Kết Hydro", zlabel = "Mức Độ Tổn Thương (Entropy)",
               xlabelcolor = :white, ylabelcolor = :white, zlabelcolor = :white,
               xgridcolor = :gray15, ygridcolor = :gray15, zgridcolor = :gray15,
               backgroundcolor = :black)

    # Khởi tạo 300 phần tử tế bào biểu diễn luồng sinh học
    Random.seed!(777)
    num_cells = 300
    
    angles = rand(num_cells) .* 2 * pi
    radii = 1.5 .+ rand(num_cells) .* 2.5
    
    x_init = radii .* sin.(angles)
    y_init = radii .* cos.(angles)
    z_init = 45.0 .+ rand(num_cells) .* 15.0 

    t_axis = 0.0:0.04:10.0
    t_fusion = 5.0 

    # 2. VÒNG LẶP ÉP QUỸ ĐẠO TẾ BÀO SỤP ĐỔ VÀ TÁI TẠO
    for i in 1:num_cells
        x_path, y_path, z_path = Float64[], Float64[], Float64[]
        
        current_x = x_init[i]
        current_y = y_init[i]
        current_z = z_init[i]
        
        for t in t_axis
            push!(x_path, current_x)
            push!(y_path, current_y)
            push!(z_path, current_z)
            
            if t < t_fusion
                # GIAI ĐOẠN LÃO HÓA / TỔN THƯƠNG (0 - 5 GIÂY)
                freq = 6.0 * t
                current_x = radii[i] * sin(angles[i] + freq) + rand(-0.15:0.01:0.15)
                current_y = radii[i] * cos(angles[i] + freq) + rand(-0.15:0.01:0.15)
                current_z = z_init[i] + rand(-2.0:0.1:2.0) 
            else
                # GIAI ĐOẠN KHÓA PHA VẠN NĂNG (5 - 10 GIÂY)
                decay = exp(-2.5 * (t - t_fusion))
                
                freq_lock = 6.0 * t_fusion
                current_x = radii[i] * sin(angles[i] + freq_lock) * decay
                current_y = radii[i] * cos(angles[i] + freq_lock) * decay
                current_z = z_init[i] * decay 
            end
        end
        
        cell_color = i % 2 == 0 ? :green : :lime
        lines!(ax, x_path, y_path, z_path, color = cell_color, linewidth = 1.3, alpha = 0.3)
    end

    # 3. ĐÓNG DẤU LÕI PHÔI THAI TRỨNG NƯỚC TỐI CAO (MÀU HỒNG CÁNH SEN)
    meshscatter!(ax, [0.0], [0.0], [0.0], color = :magenta, markersize = 0.38, 
                 label = "Lõi Phôi Thai Totipotency (0 Entropy)")

    # Ép camera ghim chặt tự động căng biên độ
    autolimits!(ax)

    # Nét đứt ranh giới kích hoạt dòng năng lượng tự chữa lành tại giây thứ 5.0
    lines!(ax, [5.0, 5.0], [-4.0, 4.0], [0.0, 60.0], color = :white, linewidth = 2, linestyle = :dash)

    axislegend(ax, labelcolor = :white, backgroundcolor = :gray10, position = :rt)

    # 4. TỰ ĐỘNG XUẤT BẢN FILE ẢNH CHỨNG MINH SINH HỌC LƯỢNG TỬ
    output_image = "quantum_totipotency_proof.png"
    save(output_image, fig)

    # Bật mở cửa sổ 3D OpenGL trực tiếp
    display(fig)
    println("=========================================================")
    println("XỬ LÝ PHÓNG TO TRÀN MÀN HÌNH TẾ BÀO GỐC THÀNH CÔNG RỰC RỠ!")
    println("File ảnh tự chữa lành lượng tử đã lưu tại: ", output_image)
    println("=========================================================")
end

run_totipotency_simulation_perfect()
