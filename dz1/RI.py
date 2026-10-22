import numpy as np
import matplotlib.pyplot as plt

# ========================= ПАРАМЕТРЫ =========================
c = 3e8                         # скорость света, м/с

# Частоты
f_start = 10e9                  # 10 ГГц
f_stop = 16e9                   # 16 ГГц
num_freq = 16                   # число частот
freqs = np.linspace(f_start, f_stop, num_freq)


# Передающая АР (квадратная, шаг 0.25 м)
Nt = 8                          # 8x8 = 64 элемента
step_t = 0.25
x_t = np.linspace(-(Nt-1)*step_t/2, (Nt-1)*step_t/2, Nt)
y_t = x_t.copy()
X_t, Y_t = np.meshgrid(x_t, y_t)
tx_coords = np.stack([X_t.ravel(), Y_t.ravel(), np.zeros_like(X_t.ravel())], axis=1)

# Приёмная АР (квадратная, шаг 0.288 м)
Nr = 8
step_r = 0.288
x_r = np.linspace(-(Nr-1)*step_r/2, (Nr-1)*step_r/2, Nr)
y_r = x_r.copy()
X_r, Y_r = np.meshgrid(x_r, y_r)
rx_coords = np.stack([X_r.ravel(), Y_r.ravel(), np.zeros_like(X_r.ravel())], axis=1)

# Объект: крест из 9 точек в плоскости z = z0
z0 = 1.0                        # расстояние до объекта, м
arm_len = 0.05                  # длина плеча, м
step_obj = 0.025                # шаг между точками

obj_points = []
# Центр
obj_points.append([0, 0, z0])
# Горизонтальное плечо
for x in np.arange(-arm_len, arm_len+1e-9, step_obj):
    if abs(x) > 1e-9:
        obj_points.append([x, 0, z0])
# Вертикальное плечо
for y in np.arange(-arm_len, arm_len+1e-9, step_obj):
    if abs(y) > 1e-9:
        obj_points.append([0, y, z0])
obj_points = np.array(obj_points)
num_points = obj_points.shape[0]
print(f"Объект: {num_points} точечных отражателей")

# Весовые коэффициенты (окно Хэмминга, α=0.54)
alpha = 0.54
m = np.arange(num_freq)
w_m = alpha - (1 - alpha) * np.cos(2 * np.pi * m / (num_freq - 1))

# ==================== ФОРМИРОВАНИЕ СИГНАЛА ====================
# Для простоты интенсивность всех отражателей a_p = 1
# Комплексные амплитуды V_m(n,k) - сумма по отражателям
def compute_Vm(tx, rx, obj_pts, freqs):
    # tx: (N_tx, 3), rx: (N_rx, 3), obj_pts: (P, 3)
    N_t = tx.shape[0]
    N_r = rx.shape[0]
    P = obj_pts.shape[0]
    M = len(freqs)
    V = np.zeros((N_t, N_r, M), dtype=complex)
    for p in range(P):
        rp = obj_pts[p]
        dist_t = np.linalg.norm(tx[:, None, :] - rp[None, None, :], axis=2)   # (N_t, 1)
        dist_r = np.linalg.norm(rx[None, :, :] - rp[None, None, :], axis=2)   # (1, N_r)
        tau = (dist_t + dist_r) / c                     # (N_t, N_r)
        phase = -np.outer(tau.ravel(), 2*np.pi*freqs)   # (N_t*N_r, M)
        V += np.exp(1j * phase).reshape(N_t, N_r, M)
    return V   # (N_t, N_r, M)

V_m = compute_Vm(tx_coords, rx_coords, obj_points, freqs)
print("Массив V_m сформирован:", V_m.shape)

# ==================== ФОКУСИРОВКА ИЗОБРАЖЕНИЯ ====================
def focus_image(tx, rx, Vm, freqs, z_plane, x_range, y_range, res, w):
    """
    Обратные проекции для восстановления РИ в плоскости z = z_plane.
    tx, rx: координаты передатчиков и приёмников (N_tx,3), (N_rx,3)
    Vm: массив (N_tx, N_rx, M) комплексных амплитуд
    freqs: (M,)
    z_plane: координата z плоскости изображения
    x_range, y_range: (xmin, xmax, ymin, ymax)
    res: разрешение сетки (количество точек на метр)
    w: весовые коэффициенты (M,)
    """
    x = np.linspace(x_range[0], x_range[1], int((x_range[1]-x_range[0])*res)+1)
    y = np.linspace(y_range[0], y_range[1], int((y_range[1]-y_range[0])*res)+1)
    X, Y = np.meshgrid(x, y)
    Z = z_plane * np.ones_like(X)
    points = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)  # (N_pix, 3)
    
    Nt = tx.shape[0]
    Nr = rx.shape[0]
    M = len(freqs)
    image = np.zeros(points.shape[0], dtype=complex)
    
    for n in range(Nt):
        for k in range(Nr):
            # Бистатические задержки для всех точек изображения
            d_t = np.linalg.norm(points - tx[n], axis=1)
            d_r = np.linalg.norm(points - rx[k], axis=1)
            tau = (d_t + d_r) / c                  # (N_pix,)
            # Сумма по частотам
            phase = np.outer(tau, 2*np.pi*freqs)   # (N_pix, M)
            S = np.sum(w * Vm[n, k, :] * np.exp(1j * phase), axis=1)
            image += S
    return X, Y, image.reshape(X.shape)

# Область восстановления (плоскость z0)
xlim = (-0.2, 0.2)
ylim = (-0.2, 0.2)
resolution = 80                 # точек/метр -> примерно 33x33 пикселя
X_img, Y_img, img = focus_image(tx_coords, rx_coords, V_m, freqs, z0,
                                xlim, ylim, resolution, w_m)

# Визуализация
plt.figure(figsize=(7,6))
plt.imshow(np.abs(img), extent=[xlim[0], xlim[1], ylim[0], ylim[1]],
           origin='lower', cmap='hot', aspect='equal')
plt.colorbar(label='Амплитуда')
plt.scatter(obj_points[:,0], obj_points[:,1], c='cyan', s=20, marker='+', label='Истинные точки')
plt.title(f'Радиоизображение креста, z0 = {z0} м, полоса {f_start/1e9}-{f_stop/1e9} ГГц')
plt.xlabel('x, м')
plt.ylabel('y, м')
plt.legend()
plt.tight_layout()
plt.show()