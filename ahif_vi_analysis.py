import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import pearsonr
from tqdm import tqdm
import os

green_color_code = "\033[38;5;2m"
# 设置全局字体为新罗马，默认大小为32
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 32
# 设置全局线条宽度为2，包括坐标轴线条和刻度线条
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.width'] = 2


# 检查两条线段是否相交
def check_intersection(p1, p2, q1, q2):
    def ccw(a, b, c):
        return (c[1] - a[1]) * (b[0] - a[0]) > (b[1] - a[1]) * (c[0] - a[0])

    return (ccw(p1, q1, q2) != ccw(p2, q1, q2)) and (ccw(p1, p2, q1) != ccw(p1, p2, q2))


# 计算交点
def compute_intersection(p1, p2, q1, q2):
    A1 = p2[1] - p1[1]
    B1 = p1[0] - p2[0]
    C1 = A1 * p1[0] + B1 * p1[1]

    A2 = q2[1] - q1[1]
    B2 = q1[0] - q2[0]
    C2 = A2 * q1[0] + B2 * q1[1]

    determinant = A1 * B2 - A2 * B1
    if determinant == 0:
        return None  # 平行，无交点

    x = (B2 * C1 - B1 * C2) / determinant
    y = (A1 * C2 - A2 * C1) / determinant
    return (x, y)


# 检查自交并计算所有交点
def find_intersections(x, y):
    intersections = []
    indices = []  # 存储交点索引的列表
    num_points = len(x)

    for i in range(num_points):
        for j in range(i + 2, num_points - (i == 0)):  # 避免相邻边
            p1 = (x[i], y[i])
            p2 = (x[(i + 1) % num_points], y[(i + 1) % num_points])
            q1 = (x[j], y[j])
            q2 = (x[(j + 1) % num_points], y[(j + 1) % num_points])

            if check_intersection(p1, p2, q1, q2):
                intersection = compute_intersection(p1, p2, q1, q2)
                if intersection:
                    intersections.append(intersection)
                    indices.append((i, (i + 1) % num_points, j, (j + 1) % num_points))  # 记录交点所在的索引
    return intersections, indices  # 返回交点和索引


# 分割单周期数据
def resort_vichar(x, y):
    # 找到 y 中的最大值的索引
    max_index = np.argmax(x)

    # 重新排序 x 和 y
    sorted_x = np.concatenate((x[max_index:], x[:max_index]))
    sorted_y = np.concatenate((y[max_index:], y[:max_index]))
    return sorted_x, sorted_y


# 分割单周期数据
def seg_vichar(x, y, intersections, indices):
    split_points = []
    for idx_group in indices:
        for i in range(0, len(idx_group) - 1, 2):
            mid_point = (idx_group[i] + idx_group[i + 1]) / 2
            split_points.append(mid_point)
    # print(split_points)
    split_points.sort()
    split_points = [i + 1 for i in split_points]
    split_points = [0] + split_points + [len(x)]

    x_segd = []
    for i in range(len(split_points) - 1):
        x_segd.append(list(x[int(split_points[i]):int(split_points[i + 1])]))
    for i in range(len(intersections)):
        x_segd[i].append(intersections[i][0])
        x_segd[i + 1].insert(0, intersections[i][0])  # 在列表左边插入元素
        x_segd[-i - 1].insert(0, intersections[i][0])  # 在列表左边插入元素
        x_segd[-i - 2].append(intersections[i][0])
    # 输出结果
    for i, sublist in enumerate(x_segd):
        # print(f"x_List {i + 1}: {sublist}")
        pass

    y_segd = []
    for i in range(len(split_points) - 1):
        y_segd.append(list(y[int(split_points[i]):int(split_points[i + 1])]))
    for i in range(len(intersections)):
        y_segd[i].append(intersections[i][1])
        y_segd[i + 1].insert(0, intersections[i][1])  # 在列表左边插入元素
        y_segd[-i - 1].insert(0, intersections[i][1])  # 在列表左边插入元素
        y_segd[-i - 2].append(intersections[i][1])
    # 输出结果
    for i, sublist in enumerate(y_segd):
        # print(f"y_List {i + 1}: {sublist}")
        pass
    x_segd_contact = []
    y_segd_contact = []
    for i in range(int((len(x_segd) - 1) / 2)):
        x_segd_contact.append(x_segd[i] + x_segd[-i - 1])
        y_segd_contact.append(y_segd[i] + y_segd[-i - 1])
    x_segd_contact.append(x_segd[int((len(x_segd) - 1) / 2)])
    y_segd_contact.append(y_segd[int((len(y_segd) - 1) / 2)])

    return x_segd_contact, y_segd_contact


# 鞋带公式算面积
def polygon_area(vertices):
    n = len(vertices)  # Number of vertices
    area = 0
    for i in range(n):
        x1, y1 = vertices[i]
        x2, y2 = vertices[(i + 1) % n]  # Ensures the last point connects to the first
        area += x1 * y2 - x2 * y1
    return area / 2


def getcv(data, num, normalize=False):
    N = 1
    start_point = num
    end_point = num + 10000
    fs = 8000 / (data['in s'].iloc[8000] - data['in s'].iloc[0])
    timestamp = data['in s'].iloc[start_point:end_point]
    T = 1 / fs

    timestamp = data['in s'].iloc[start_point:start_point + int((fs / 50) * N)]
    current = data['C2 in A'].iloc[start_point:start_point + int((fs / 50) * N)]

    current = current - data['C2 in A'].mean()
    timestamp = 10 * (timestamp - timestamp.iloc[0]).reset_index(drop=True)

    voltage = data['C1 in V'].iloc[start_point:start_point + int((fs / 50) * N)].tolist()
    correlation, _ = pearsonr(current, pd.Series(voltage))
    if correlation < 0:
        voltage = [-i for i in voltage]

    voltage = pd.Series(voltage)
    # 将电流和电压数据归一化到-1到1范围,这段代码可选，注意修改标幺值
    max_current = max(current)
    min_current = min(current)
    max_voltage = max(voltage)
    min_voltage = min(voltage)
    max_timestamp = max(timestamp)
    min_timestamp = min(timestamp)

    if normalize:
        current = 2 * (current - min_current) / (max_current - min_current) - 1
        voltage = 2 * (voltage - min_voltage) / (max_voltage - min_voltage) - 1
    timestamp = 20 * (timestamp - min_timestamp) / (max_timestamp - min_timestamp)

    return current, voltage, timestamp, int(fs / 50)


def get_area(x, y):
    x, y = resort_vichar(x, y)

    intersections, indices = find_intersections(x, y)
    if intersections != [] and len(intersections) < 5:
        x_segd, y_segd = seg_vichar(x, y, intersections, indices)

        parts = []
        parts_area = []
        for i in range(len(x_segd)):
            parts.append([(x_point, y_point) for x_point, y_point in zip(x_segd[i], y_segd[i])])
        for i in range(len(x_segd)):
            parts_area.append(polygon_area(parts[i]))
        parts_area_L = [i for i in parts_area if i < 0]
        parts_area_C = [i for i in parts_area if i > 0]
    else:
        area = polygon_area([(i, j) for i, j in zip(x, y)])
        if area > 0:
            parts_area_L = [0]
            parts_area_C = [area]
        else:
            parts_area_L = [area]
            parts_area_C = [0]

    return abs(sum(parts_area_L)), abs(sum(parts_area_C))


# 给电流可视化
def show_current(data, num, N):
    start_point = num
    end_point = num + 10000
    fs = 8000 / (data['in s'].iloc[8000] - data['in s'].iloc[0])
    timestamp = data['in s'].iloc[start_point:end_point]
    T = 1 / fs
    timestamp = data['in s'].iloc[start_point:start_point + int((fs / 50) * N)]
    current = data['C2 in A'].iloc[start_point:start_point + int((fs / 50) * N)]
    current = current - data['C2 in A'].mean()
    timestamp = (timestamp - timestamp.iloc[0]).reset_index(drop=True)
    # 绘制 current
    plt.figure(figsize=(16, 6))
    plt.plot(timestamp, current, label='Current')
    plt.xlim([0, 0.02 * N])
    plt.xlabel("Time(s)")
    plt.ylabel("Current (A)")

    plt.legend(loc='upper left', bbox_to_anchor=(0, 1.25), frameon=False)
    plt.tight_layout()
    plt.show()


from scipy.fft import fft


def calculate_SPQD(voltage, current, sampling_rate):
    # 确保电压和电流信号是 NumPy 数组
    voltage = np.asarray(voltage)
    current = np.asarray(current)

    # 获取信号长度
    N = len(voltage)

    # 进行快速傅里叶变换 (FFT)
    V_fft = fft(voltage)
    I_fft = fft(current)

    # 计算频率分量对应的频率
    freqs = np.fft.fftfreq(N, d=1 / sampling_rate)

    # 计算电压和电流各次谐波的幅值
    V_magnitude = np.abs(V_fft) / N
    I_magnitude = np.abs(I_fft) / N

    # 计算电压和电流各次谐波的相位
    V_phase = np.angle(V_fft)
    I_phase = np.angle(I_fft)

    # 初始化总功率
    S_total = 0.0
    Q_total = 0.0
    P_total = 0.0

    # 遍历每个谐波频率分量（忽略直流分量 freq[0]）
    for n in range(1, N // 2):
        # 计算相位差
        phase_diff = V_phase[n] - I_phase[n]
        # 计算该谐波的有功功率分量
        P_n = V_magnitude[n] * I_magnitude[n] * np.cos(phase_diff)
        # 计算该谐波的无功功率
        Q_n = V_magnitude[n] * I_magnitude[n] * np.sin(phase_diff)
        # 累加至总功率
        P_total += P_n
        Q_total += Q_n

    U_eff = np.sqrt(sum([(V_magnitude[n]) ** 2 for n in range(1, N // 2)]))
    I_eff = np.sqrt(sum([(I_magnitude[n]) ** 2 for n in range(1, N // 2)]))

    S_total = U_eff * I_eff
    D_total = np.sqrt(S_total ** 2 - P_total ** 2 - Q_total ** 2)
    # 这里仅仅计算了正频率部分，考虑负频率应该乘2
    return 2 * S_total, 2 * P_total, 2 * Q_total, 2 * D_total


if __name__ == '__main__':

    # 代码开始点
    num = 32850
    N = 400
    filename='131.CSV'
    base_fs=50
    file_path = os.path.join(os.path.dirname(__file__), 'demo_data', filename)
    data = pd.read_csv(file_path)


    # 是否归一化
    normalize_label = False
    print(f"First Point:{num}")
    if normalize_label:
        print("The area has been normalized!")



    show_current(data, num, N)
    current, voltage, timestamp, Nf = getcv(data, num, normalize=normalize_label)
    fs = 50 * Nf
    print(f'Sampling rate{fs}')

    L_area = []
    C_area = []
    S = []
    P = []
    Q = []
    D = []
    Q_byArea = []

    # 假设 L_area, C_area, P, Q 已经初始化
    for i in tqdm(range(N), desc=f"{green_color_code}process", ncols=100):
        current, voltage, timestamp, Nf = getcv(data, num, normalize=normalize_label)

        S_i, P_i, Q_i, D_i = calculate_SPQD(voltage, current, fs)

        num = num + Nf
        x = current.to_numpy()
        y = voltage.to_numpy()

        L_area_i, C_area_i = get_area(x, y)

        # 容性无功计算方法
        Q_byArea_i = -sum([((C_area_i - L_area_i) / (2 * math.pi))])

        L_area.append(L_area_i)
        C_area.append(C_area_i)
        Q_byArea.append(Q_byArea_i)

        S.append(S_i)
        P.append(P_i)
        Q.append(Q_i)
        D.append(D_i)

    if base_fs==50:
        t_values = [x * 0.02 for x in range(len(L_area))]  # 横轴重新赋时间
    else:
        t_values = [x *(1/fs) for x in range(len(L_area))]

    # 绘制感容面积
    plt.figure(figsize=(16, 6))
    plt.plot(t_values, [3 * x for x in L_area], label='L_Area', color='blue', linestyle='-', linewidth=3)
    plt.plot(t_values, [3 * x for x in C_area], label='C_Area', color='red', linestyle='-', linewidth=3)
    plt.xlabel('Time(s)')
    plt.ylabel('Area')
    plt.legend(loc='upper left', bbox_to_anchor=(0, 1.25), frameon=False, ncol=2)
    plt.xlim([0, 0.02 * N])
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(16, 6))
    # plt.plot(t_values, P, label='P', marker='s', color='orange', linestyle='-')
    plt.plot(t_values, Q, label='Q', marker='s', color='green', linestyle='-')
    # plt.plot(t_values, D, label='D', marker='s', color='blue', linestyle='-')
    plt.plot(t_values, Q_byArea, label='Q_byArea', marker='s', color='red', linestyle='-')

    plt.legend(loc='upper left', bbox_to_anchor=(0, 1.25), frameon=False, ncol=3)

    plt.ylabel('P or Q')
    # 计算功率因数和角
    PF_tan = [q / p for p, q in zip(P, Q)]
    # 计算功率因数角（弧度）
    power_factor_angle = [math.degrees(math.atan(pf)) for pf in PF_tan]

    # 创建右侧坐标轴
    ax2 = plt.gca().twinx()
    # 绘制功率因数角
    plt.plot(t_values, power_factor_angle, label='Power Factor Angle (deg)', marker='o', color='black', linestyle='-')
    plt.legend(loc='upper right',bbox_to_anchor=(1, 1.25), frameon=False)

    plt.xlabel('Time (s)')
    plt.ylabel('P & -Q')
    ax2.set_ylabel('PF_Angle (deg)')
    plt.xlim([0, 0.02 * N])
    plt.tight_layout()
    plt.show()
