import numpy as np

dt = np.dtype([('time', np.float64), ('id', np.int64), ('x', np.float64), \
('y', np.float64), ('z', np.float64)])

data = np.fromfile("test.bin", dt)

print(data)
