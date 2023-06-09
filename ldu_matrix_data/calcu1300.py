import matplotlib.pyplot as plt
import array
import numpy as np
idx0 = []

with open('uidx.txt', 'r') as file:
    x0 = [float(x) for x in file.readlines()]
    for one in x0:
        idx0.append(one / 1381)

idx1 = []
with open('uidx.txt', 'r') as file: 
    x1 = [88360 - float(y) for y in file.readlines()]
    for one in x1:
        idx1.append(one / 1381)
# bucket_slice = bucket[:100]
bucket = array.array('i', [0] * 64)

for one in idx0:
    bucket[int(one)] += 1

for one in idx1:
    bucket[int(one)] += 1

x = np.arange(64)
bucket_slice = bucket[:64]
# 绘制柱状图
plt.bar(x, bucket_slice)

# 设置x轴的刻度标签为区间的值
plt.xticks(x, x)

# 设置标题和轴标签
plt.title("Bucket Histogram")
plt.xlabel("Value")
plt.ylabel("Frequency")

# 显示图形
plt.show()
