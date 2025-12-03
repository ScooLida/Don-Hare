import statistics

data = [10, 12, 23, 23, 16, 23, 21, 16]

# 1. Среднее (Mean)
mean_val = statistics.mean(data)

# 2. Разброс: Стандартное отклонение (Standard Deviation)
# Показывает, насколько сильно данные отклоняются от среднего
std_dev = statistics.stdev(data)
