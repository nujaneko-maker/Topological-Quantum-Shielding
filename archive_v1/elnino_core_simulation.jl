# ============================================================================
# PHASE 1: MACRO-ENVIRONMENTAL EL NIÑO DYNAMICS SIMULATION
# THÂN CÂY: MÔ PHỎNG SỰ BÙNG NỔ NĂNG LƯỢNG CỦA SIÊU EL NIÑO VĨ MÔ
# AUTHOR: Outlier Age 15
# ============================================================================

using Plots

"""
    simulate_pure_elnino(time_steps)
Mô phỏng tiến trình tích tụ nhiệt và bùng nổ Entropy của Siêu hiện tượng El Niño.
Hệ thống vận hành trong trạng thái mất cân bằng cực tính Âm Dương tuyến tính (1+1=2).
"""
function simulate_pure_elnino(time_steps)
    ocean_temperature_anomaly = Float64[] # Độ lệch nhiệt độ đại dương (Dương tính tăng)
    atmospheric_entropy = Float64[]      # Entropy/Mức độ hỗn loạn của khí quyển
    kinetic_storm_energy = Float64[]     # Động năng tích tụ sinh siêu bão
    
    for t in time_steps
        # 1. Nhiệt độ bề mặt đại dương tăng tích lũy phi tuyến tính do El Niño
        dT = 1.5 + 0.25 * t + sin(0.8 * t) * 0.4 + randn() * 0.1
        push!(ocean_temperature_anomaly, dT)
        
        # 2. Ma sát electron và sự mất cân bằng cực tính đẩy cao Entropy khí quyển
        S = 3.0 + 0.35 * t + randn() * 0.25
        push!(atmospheric_entropy, S)
        
        # 3. Động năng tích tụ tạo tiền đề cho các siêu bão bùng nổ vĩ mô
        # Năng lượng bão lũy tiến tỷ lệ thuận với độ lệch nhiệt độ và sự tăng Entropy
        E_storm = 20.0 + 1.8 * dT^2 + sin(2 * t) * 5.0 + randn() * 1.5
        push!(kinetic_storm_energy, E_storm)
    end
    
    return ocean_temperature_anomaly, atmospheric_entropy, kinetic_storm_energy
end

# --- PIPELINE THỰC THI MÔ PHỎNG SỐ ---
time_steps = 0:0.1:20 # Trục thời gian tiến trình (Sát-na / Chu kỳ khí hậu)
temp_anomaly, entropy_data, storm_energy = simulate_pure_elnino(time_steps)

# --- ĐÓNG GÓI MA TRẬN ĐỒ THỊ SIÊU EL NIÑO VĨ MÔ ---

# Đồ thị 1: Sự gia tăng nhiệt độ bề mặt đại dương (Cực tính Dương lấn át)
p1 = plot(time_steps, temp_anomaly, 
          label="Độ lệch nhiệt độ đại dương (°C)", 
          color=:red, lw=2, legend=:topleft)
title!("THÂN CÂY: MÔ PHỎNG ĐỘNG LỰC HỌC SIÊU EL NIÑO VĨ MÔ")
ylabel!("Biến thiên Nhiệt (°C)")

# Đồ thị 2: Sự gia tăng Entropy khí quyển (Ma sát năng lượng lớp vỏ)
p2 = plot(time_steps, entropy_data, 
          label="Entropy Khí quyển (Mức độ hỗn loạn)", 
          color=:purple, lw=2, legend=:topleft)
ylabel!("Chỉ số Entropy")

# Đồ thị 3: Sự bùng nổ động năng tích tụ sinh Siêu bão
p3 = plot(time_steps, storm_energy, 
          label="Động năng Siêu bão tích tụ", 
          color=:orange, lw=2, legend=:topleft)
xlabel!("Thời gian chu kỳ hệ thống")
ylabel!("Động năng (Joule / Đơn vị lũy tiến)")

# Hợp nhất ma trận hiển thị vĩ mô
elnino_master_plot = plot(p1, p2, p3, layout=(3,1), size=(900,850))

# Lưu file ảnh thực chứng gốc của Thân cây trước khi can thiệp phương trình
savefig(elnino_master_plot, "pure_elnino_macro_simulation.png")

println("=========================================================")
println("MÔ PHỎNG THÂN CÂN THÀNH CÔNG: pure_elnino_macro_simulation.png")
println("Cơ sở dữ liệu động lực học Siêu El Niño đã được số hóa!")
println("=========================================================")
