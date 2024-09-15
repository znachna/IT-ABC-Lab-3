import os
import json
import matplotlib.pyplot as plt

script_dir = os.path.dirname(os.path.abspath(__file__))
json_file_path = os.path.join(script_dir, "plots.json")

with open(json_file_path, 'r') as file:
    data = json.load(file)

def plot_integral_deviation(data):
    if "Integral deviation" in data:
        integral_data = data["Integral deviation"]
        precisions = [item['precision'] for item in integral_data]
        steps = [item['steps'] for item in integral_data]

        plt.figure(figsize=(6, 4))
        plt.plot(precisions, steps, marker='o')
        plt.xlabel('Precision')
        plt.ylabel('Steps')
        plt.title('Integral Deviation')
        plt.grid(True)
        plt.show()

def plot_derivatives(data, derivative_type, h_params, remove_edges=False):
    plt.figure(figsize=(10, 6))

    for h_param in h_params:
        h_param_str = f"{h_param:.6f}"
        if h_param_str in data[derivative_type]:
            x_vals = [item['x'] for item in data[derivative_type][h_param_str]]
            deviations = [item['deviation'] for item in data[derivative_type][h_param_str]]

            if remove_edges and len(x_vals) > 2:
                x_vals = x_vals[1:-1]
                deviations = deviations[1:-1]

            plt.plot(x_vals, deviations, label=f'H_PARAM: {h_param_str}', marker='o')

    plt.xlabel('x')
    plt.ylabel('Deviation')
    plt.title(f'{derivative_type} for Different H_PARAM')
    plt.legend()
    plt.grid(True)
    plt.show()

plot_integral_deviation(data)

h_params = [0.05, 0.1, 0.2]

if "First derivative" in data:
    plot_derivatives(data, "First derivative", h_params)
    plot_derivatives(data, "First derivative", h_params, remove_edges=True)

if "Second derivative" in data:
    plot_derivatives(data, "Second derivative", h_params, remove_edges=True)
