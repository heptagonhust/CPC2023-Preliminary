import matplotlib.pyplot as plt

# 读取横坐标文件
with open('uidx.txt', 'r') as file:
    x_coordinates = [float(x) for x in file.readlines()]

# 读取纵坐标文件
with open('lidx.txt', 'r') as file:
    y_coordinates = [0 - float(y) for y in file.readlines()]

with open('lidx.txt', 'r') as file:
    x_coordinates0 = [float(x) for x in file.readlines()]

with open('uidx.txt', 'r') as file:
    y_coordinates0 = [0 - float(x) for x in file.readlines()]



# 绘制二维坐标图
plt.scatter(x_coordinates, y_coordinates, marker='.', s=0.1)
plt.scatter(x_coordinates0, y_coordinates0, marker='.',s=0.1)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Coordinate Visualization')
plt.show()
