import json
import os
import matplotlib.pyplot as plt

# Определяем путь к файлам
script_dir = os.path.dirname(os.path.abspath(__file__))
config_file_path = os.path.join(script_dir, 'config.json')
plots_file_path = os.path.join(script_dir, 'plots.json')

# Загрузка конфигурации
with open(config_file_path, 'r') as config_file:
    config = json.load(config_file)

# Загрузка данных из plots.json
with open(plots_file_path, 'r') as data_file:
    data = json.load(data_file)

# Настройки графиков
plt.style.use('ggplot')
colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']

h_params = config.get("H_PARAMS", [0.2, 0.1, 0.05])

# График отклонения первой производной
plt.figure(figsize=(10, 6))
for idx, h_param in enumerate(h_params):
    h_param_str = str(h_param)
    if h_param_str in data['First derivative']:
        deviations = data['First derivative'][h_param_str]
        x_values = [item['x'] for item in deviations]
        y_values = [item['deviation'] for item in deviations]
        plt.plot(x_values, y_values, label=f'h = {h_param_str}', color=colors[idx % len(colors)])

plt.title('First Derivative Deviation vs X')
plt.xlabel('X')
plt.ylabel('Deviation')
plt.legend()
plt.grid(True)
plt.savefig('first_derivative_deviation.png')
plt.show()

# Дополнительный график для первой производной с отсечением первой и последней точек
plt.figure(figsize=(10, 6))
for idx, h_param in enumerate(h_params):
    h_param_str = str(h_param)
    if h_param_str in data['First derivative']:
        deviations = data['First derivative'][h_param_str]
        x_values = [item['x'] for item in deviations]
        y_values = [item['deviation'] for item in deviations]

        # Отсекаем первую и последнюю точки
        if len(x_values) > 2 and len(y_values) > 2:
            x_values = x_values[1:-1]
            y_values = y_values[1:-1]

        plt.plot(x_values, y_values, label=f'h = {h_param_str}', color=colors[idx % len(colors)])

plt.title('First Derivative Deviation vs X (Without Edge Points)')
plt.xlabel('X')
plt.ylabel('Deviation')
plt.legend()
plt.grid(True)
plt.savefig('first_derivative_deviation_no_edges.png')
plt.show()

# График отклонения второй производной с отсечением первой и последней точек
plt.figure(figsize=(10, 6))
for idx, h_param in enumerate(h_params):
    h_param_str = str(h_param)
    if h_param_str in data['Second derivative']:
        deviations = data['Second derivative'][h_param_str]
        x_values = [item['x'] for item in deviations]
        y_values = [item['deviation'] for item in deviations]

        # Отсекаем первую и последнюю точки
        if len(x_values) > 2 and len(y_values) > 2:
            x_values = x_values[1:-1]
            y_values = y_values[1:-1]

        plt.plot(x_values, y_values, label=f'h = {h_param_str}', color=colors[idx % len(colors)])

plt.title('Second Derivative Deviation vs X (Without Edge Points)')
plt.xlabel('X')
plt.ylabel('Deviation')
plt.legend()
plt.grid(True)
plt.savefig('second_derivative_deviation_no_edges.png')
plt.show()

# График количества шагов интегрирования в зависимости от точности
if "Integral" in data:
    precisions = [item['precision'] for item in data['Integral']]
    steps = [item['steps'] for item in data['Integral']]

    plt.figure(figsize=(10, 6))
    plt.plot(precisions, steps, marker='o', linestyle='-')
    plt.title('Integral Steps vs Precision')
    plt.xlabel('Precision')
    plt.ylabel('Number of Steps')
    plt.xscale('log')
    plt.gca().invert_xaxis()
    plt.grid(True)
    plt.savefig('integral_steps_vs_precision.png')
    plt.show()
else:
    print("Key 'Integral' not found in data.")

print("Графики успешно сохранены.")
