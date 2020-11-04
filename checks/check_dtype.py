import numpy as np
dt = [('ID_SRC', '<i8'), ('ID_CLUSTER', '<i8'), ('detUID', '2a')]
x = [(1,1,"a"),(2,2,"b")]
print(x)
print(np.shape(x))
ar = np.array(x[0], dtype=dt)
print(ar)

y = [(1,2),(1,2), ("a","b")]
ar = np.array(x[0], dtype=dt, order="F")
print(ar)
