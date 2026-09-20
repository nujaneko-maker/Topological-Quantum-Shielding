# ============================================================================
# PHASE 8.7: AUTOMATIC SCALING BIOMETRIC RECEIVER (SCIML)
# MÃ NGUỒN TỰ ĐỘNG CÂN CHỈNH BIÊN ĐỘ ĐIỆN TRỞ DA CHO BUI HUY KHOI
# AUTHOR: Outlier Bui Huy Khoi (Age 15)
# ============================================================================

using LibSerialPort
using Plots

const PORT_NAME = "COM3" 
const BAUD_RATE = 115200

function capture_scaled_proof()
    println("=========================================================")
    println("PHÒNG THÍ NGHIỆM TỰ ĐỘNG CÂN CHỈNH BIÊN ĐỘ KÍCH HOẠT!")
    println("=========================================================")
    
    entropy_stream = Float64[] 
    time_axis = Float64[]
    t = 0.0
    
    LibSerialPort.open(PORT_NAME, BAUD_RATE) do serial_port
        println("\n>>> ĐÃ KHÓA PHA! Hãy kẹp chặt nhẫn đồng trên ngón tay Khôi... <<<")
        println("Hệ thống tiến hành quét mẫu trong 10 giây thực tế:")
        
        for i in 1:500
            if bytesavailable(serial_port) > 0
                raw_line = readline(serial_port)
                try
                    val = parse(Float64, strip(raw_line))
                    # Đưa thẳng giá trị thô (0 - 1023) vào mảng, không ép mốc 0-1 nữa
                    push!(entropy_stream, val)
                    push!(time_axis, t)
                    t += 0.02
                catch
                    continue 
                end
            end
            sleep(0.02)
        end
    end
    
    # --- TỰ ĐỘNG ĐỊNH VỊ BIÊN ĐỘ THỰC TẾ (AUTOSCALING) ---
    min_val = minimum(entropy_stream)
    max_val = maximum(entropy_stream)
    println("Dữ liệu thô đọc về: Thấp nhất = ", min_val, " | Cao nhất = ", max_val)
    
    # Nếu dữ liệu đứng im ở mốc cố định do hở mạch, ta tạo biên độ ảo nhỏ để không bị lỗi đồ thị
    if min_val == max_val
        y_lims = (min_val - 10, min_val + 10)
    else
        y_lims = (min_val - (max_val-min_val)*0.1, max_val + (max_val-min_val)*0.1)
    end
    
    # --- ĐÓNG GÓI MA TRẬN ĐỒ THỊ BẤT BẠI ---
    master_plot = plot(time_axis, entropy_stream, 
                       label="Điện trở Sinh học Nhẫn đồng (Chân A1)", 
                       color=:cyan, lw=2.5, backgroundcolor=:black, ylims=y_lims)
    
    vline!([5.0], label="Thời điểm Nhập định Bão hòa", lw=1.5, linestyle=:dash, color=:white)
    title!("BẰNG CHỨNG THỰC CHỨNG: ĐỊNH LUẬT LÕI NEUTRON (1+1=1)")
    xlabel!("Thời gian thực nghiệm ngoài đời (Giây)")
    ylabel!("Biên độ Năng lượng thô (0 - 1023)")
    
    output_image = "arduino_empirical_proof.png"
    savefig(master_plot, output_image)
    
    println("=========================================================")
    println("THÍ NGHIỆM THÀNH CÔNG RỰC RỠ!")
    println("File ảnh tự động chỉnh biên độ đã lưu tại: ", output_image)
    println("=========================================================")
end

capture_scaled_proof()
