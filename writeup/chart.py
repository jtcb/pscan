#!/usr/bin/env python3
# plot

import matplotlib.pyplot as plt
import numpy as np

counts = [50, 100, 150, 200, 250]

mem_par = [2.012, 5.450, 11.728, 20.827, 32.622]
raw_par = [2.366, 9.311, 22.410, 39.260, 61.583]

mem_ser = [3.508, 14.220, 33.009, 58.370, 92.271]
raw_ser = [6.836, 28.076, 65.354, 116.102, 183.114]

n_groups = 5

fig, ax = plt.subplots()

index = np.arange(n_groups)
bar_width = 0.16

opacity = 0.4

rects1 = plt.bar(index + -1 * bar_width, mem_par, bar_width,
                 alpha=opacity,
                 color='b',
                 label='Parallel (memoized)')

rects2 = plt.bar(index + 0 * bar_width, raw_par, bar_width,
                 alpha=opacity,
                 color='g',
                 label='Parallel')

rects3 = plt.bar(index + 1 * bar_width, mem_ser, bar_width,
                 alpha=opacity,
                 color='r',
                 label='Serial (memoized)')

rects4 = plt.bar(index + 2 * bar_width, raw_ser, bar_width,
                 alpha=opacity,
                 color="k",
                 label='Serial')

plt.xlabel('Sequence Count')
plt.ylabel('Time (seconds)')
plt.title('Clustering Performance')
plt.xticks(index + bar_width, ('50', '100', '150', '200', '250'))
plt.legend(loc="upper left")

plt.tight_layout()
plt.show()

