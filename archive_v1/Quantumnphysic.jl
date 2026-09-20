using DifferentialEquations
using Plots

# ==============================================================================
# 🏛️ QUANTUM ZENO DYNAMIC LOOP: CO-EVOLUTION OF BLACK HOLE & WHITE HOLE
# Thuật toán: "Dừng đúng lúc 20/80" tích hợp Toán tử Quan sát lượng tử O(t)
# Tác giả: Bui Huy Khoi (nujaneko-maker)
# ==============================================================================

"""
    quantum_zeno_energy!(du, u, p, t)

Hàm vi phân hệ động lực lượng tử phi tuyến mô tả sự co giãn của 2 cực đối ngẫu.
u[1] = Y1(t) : Trạng thái Lỗ Đen (Cực Dương)
u[2] = Y2(t) : Trạng thái Lỗ Trắng (Cực Âm)
u[3] = E(t)  : Công năng tích lũy lũy tiến trong Siêu tụ 5V
"""
function quantum_zeno_energy!(du, u, p, t)
    Y1, Y2, E = u
    alpha, beta, epsilon = p # Tham số đòn bẩy và biên phanh epsilon
    
    # 1. Đo khoảng cách nhị nguyên (Hỗn loạn Entropy hiện hành)
    delta = abs(Y1 - Y2)
    
    # 2. TOÁN TỬ QUAN SÁT RỜI RẠC O(t) - Tần số quét cao tương đương 50Hz của Arduino
    # Tự động chuyển pha đóng ngắt: 20% Quan sát giữ cấu trúc, 80% Thả lỏng sinh công
    observer_state = sin(2 * pi * 5 * t) > 0.0 ? 1.0 : 0.0
    psi_lever = alpha * (1.0 - exp(-beta * delta))
    
    if observer_state == 1.0
        # 🟣 20% PHA QUAN SÁT (LỖ ĐEN): Ép hệ thống dừng đúng lúc tại biên bão hòa để khóa pha
        du[1] = -5.0 * (Y1 - Y2)
        du[2] = 5.0 * (Y1 - Y2)
        
        # Tích lũy năng lượng bản thể lõi bão hòa tiệm cận 20%
        du[3] = psi_lever * (0.2 * epsilon)
    else
        # 🟠 80% PHA THẢ LỎNG (LỖ TRẮNG): Phun trào dòng vi phân biến thiên tự do để sinh công
        du[1] = 8.0 * sin(5*t)
        du[2] = -8.0 * sin(5*t)
        
        # Tận dụng tối đa dòng biến thiên vi phân động lực 80%
        actual_derivative = abs(du[1] - du[2])
        du[3] = psi_lever * (0.8 * actual_derivative)
    end
end

# 🛠️ THIẾT LẬP THAM SỐ THỰC NGHIỆM GỐC QUAN SÁT
# [Trạng thái ban đầu Lỗ Đen, Trạng thái ban đầu Lỗ Trắng, Năng lượng gốc 0J]
u0 = [50.0, 10.0, 0.0]
tspan = (0.0, 5.0)          # Quét thang thời gian 5 giây
params = [1.5, 2.5, 1.0]     # [alpha, beta, biên bão hòa tiệm cận epsilon = 1.0]

# Định nghĩa bài toán vi phân
prob = ODEProblem(quantum_zeno_energy!, u0, tspan, params)

# SỬA ĐỔI CỐT LÕI: Tắt adaptive, khóa chặt bước nhảy dt=0.002 để ép thực thể quan sát hoạt động
sol = solve(prob, Tsit5(), adaptive=false, dt=0.002, reltol=1e-6, abstol=1e-6)

# ==============================================================================
# 📊 BÓC TÁCH DỮ LIỆU PHẲNG VÀ TRỰC QUAN HÓA AN TOÀN TRÊN REPL GKS GQT
# ==============================================================================
thang_thoi_gian = sol.t
trang_thai_Y1 = [u[1] for u in sol.u] # Mảng phẳng Lỗ Đen
trang_thai_Y2 = [u[2] for u in sol.u] # Mảng phẳng Lỗ Trắng
nang_luong_E  = [u[3] for u in sol.u] # Mảng phẳng Năng lượng tích lũy

# Đồ thị 1: Vũ điệu co giãn răng cưa
p1 = plot(thang_thoi_gian, trang_thai_Y1, label="Y1(t) (Trạng thái Lỗ Đen)", color=:purple, lw=2)
plot!(thang_thoi_gian, trang_thai_Y2, label="Y2(t) (Trạng thái Lỗ Trắng)", color=:orange, lw=2)
hline!([30.0], label="Điểm bão hòa 1+1=1", color=:black, linestyle=:dot, alpha=0.5)
ylabel!("Biên độ Trạng thái")
title!("Vũ điệu Co giãn Vô hạn dưới Toán tử Quan sát O(t)")
grid!(true, alpha=0.3)

# Đồ thị 2: Dòng năng lượng tích lũy Zeno vọt dốc đi lên
p2 = plot(thang_thoi_gian, nang_luong_E, label="E_output(t) (Năng lượng tích lũy Zeno)", color=:green, lw=2.5)
xlabel!("Thời gian tau (Giây)")
ylabel!("Công năng tích lũy (J)")
grid!(true, alpha=0.3)

# Xuất hiển thị bảng biểu đồ đôi hoàn chỉnh
plot(p1, p2, layout=(2,1), size=(800,600))
