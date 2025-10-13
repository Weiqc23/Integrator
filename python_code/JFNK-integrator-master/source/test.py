import numpy as np
from points import gauss_lobatto  # 假设 points.py 提供该函数
n = 5
nodes_py = gauss_lobatto(n)  # 应返回 [-1, ..., 1] 区间的 5 个节点
print("Python nodes:", nodes_py)