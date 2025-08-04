import numpy as np
import h5py as h5
import os
import builtins 
import re
import matplotlib.pyplot as plt
import plotly.express as px
import pandas as pd
import plotly.graph_objects as go

def plot_3d_conductances(
    DB_data,
    PB_data,
    indices,
    group_labels=('DB', 'PB'),
    conductance_labels=('gX', 'gY', 'gZ'),
    highlight_point=None,  # A tuple: (x_val, y_val, z_val)
    highlight_label='highlight',
    axis_limits=None  # A dict: {'x': (min, max), 'y': (min, max), 'z': (min, max)}
):
    x_idx, y_idx, z_idx = indices
    x_label, y_label, z_label = conductance_labels

    df_DB = pd.DataFrame({
        x_label: DB_data[:, x_idx],
        y_label: DB_data[:, y_idx],
        z_label: DB_data[:, z_idx],
        'group': [group_labels[0]] * len(DB_data)
    })

    df_PB = pd.DataFrame({
        x_label: PB_data[:, x_idx],
        y_label: PB_data[:, y_idx],
        z_label: PB_data[:, z_idx],
        'group': [group_labels[1]] * len(PB_data)
    })

    df = pd.concat([df_DB, df_PB], ignore_index=True)

    fig = px.scatter_3d(df, x=x_label, y=y_label, z=z_label, color='group')

    if highlight_point:
        x_val, y_val, z_val = highlight_point
        fig.add_trace(go.Scatter3d(
            x=[x_val],
            y=[y_val],
            z=[z_val],
            mode='markers+text',
            marker=dict(size=6, symbol='diamond', color='purple'),
            text=[highlight_label],
            textposition='top center',
            name='highlight'
        ))

    # Apply axis limits if provided
    if axis_limits:
        fig.update_layout(
            scene=dict(
                xaxis=dict(range=axis_limits.get('x')),
                yaxis=dict(range=axis_limits.get('y')),
                zaxis=dict(range=axis_limits.get('z'))
            )
        )

    fig.show()

def extract_last_values(filename):
    last_values = []
    with open(filename, 'r') as file:
        for line in file:
            if "Last Value:" in line:
                # Extract values between brackets
                match = re.search(r'\[([^\]]+)\]', line)
                if match:
                    values_str = match.group(1)
                    # Convert to list of floats
                    values = [float(val) for val in values_str.split()]
                    last_values.append(values)
    return last_values

def extract_conductances(file_path):
    matrix = []
    with open(file_path, 'r') as file:
        for line in file:
            # Look for lines that start with "Final Conductances"
            if "Final Conductances|" in line:
                # Extract the five values using regex
                match = re.findall(r"g\w+=([\d\.eE+-]+)", line)
                if match and len(match) == 5:
                    # Convert to float and append to matrix
                    matrix.append([float(val) for val in match])
    
    # Convert to numpy array for matrix operations
    conductance_matrix = np.array(matrix)
    
    # Check dimensions
    if conductance_matrix.shape != (29, 5):
        raise ValueError(f"Expected a 29x5 matrix but got {conductance_matrix.shape}")
    
    return conductance_matrix


import matplotlib.pyplot as plt

def plot_error_comparison(vec1, vec2, labels=None, reference_lines=(0.05, 0.07, 0.3)):
    """
    Plots two subplots side by side showing error values against reference lines.

    Parameters:
    - vec1, vec2: Arrays or lists of 3 floats each (e.g., DB_E and PB_E).
    - labels: Tuple of labels for the subplot titles. Default is ("Vector 1", "Vector 2").
    - reference_lines: Tuple of the reference values for F, S, D (default: (0.05, 0.07, 0.3)).
    """

    if labels is None:
        labels = ("Vector 1", "Vector 2")

    colors = ["limegreen", "darkorange", "mediumorchid"]
    names = ["F", "S", "D"]
    tex_labels = [r'$\bar{F}$', r'$\bar{S}$', r'$\bar{D}$']

    plt.rcParams['text.usetex'] = True
    fig, axs = plt.subplots(1, 2, figsize=(8, 2), sharey=True)

    for i, (ax, vec, label) in enumerate(zip(axs, [vec1, vec2], labels)):
        for j in range(3):
            ax.axhline(reference_lines[j] - vec[j], color=colors[j], label=names[j])
            ax.axhline(y=reference_lines[j], linestyle='--', color=colors[j], label=tex_labels[j])
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 0.35)
        ax.set_title(label, fontsize=10)
        ax.legend(loc='lower right', fontsize=8)

    plt.tight_layout()
    plt.show()

import matplotlib.pyplot as plt
import numpy as np

def plot_conductance_comparison(ini_g, conductance_matrix, row1, row2, labels=None, titles=None, figsize=(8, 2.5), ylim=(0, 100)):
    """
    Plots side-by-side bar charts comparing initial and final conductances for two rows.

    Parameters:
    - ini_g: np.ndarray, initial conductances (shape Nx5)
    - conductance_matrix: np.ndarray, final conductances (shape Nx5)
    - row1, row2: int, indices of the rows to compare
    - labels: list of str, x-axis labels for conductances
    - titles: tuple of str, subplot titles for row1 and row2
    - figsize: tuple, size of the figure
    - ylim: tuple, y-axis limits
    """
    if labels is None:
        labels = ['gNap', 'gNa', 'gK', 'gCa', 'gCan']
    if titles is None:
        titles = (f'Row {row1}', f'Row {row2}')
    
    x = np.arange(len(labels))
    width = 0.35

    f2_g_1 = ini_g[row1, :]
    f2_final_g_1 = conductance_matrix[row1, :]
    f2_g_2 = ini_g[row2, :]
    f2_final_g_2 = conductance_matrix[row2, :]

    fig, axs = plt.subplots(1, 2, figsize=figsize, sharey=True)

    # First subplot
    axs[0].bar(x - width/2, f2_g_1, width, label='Initial g')
    axs[0].bar(x + width/2, f2_final_g_1, width, label='Final g')
    axs[0].set_title(titles[0])
    axs[0].set_xticks(x)
    axs[0].set_xticklabels(labels)
    axs[0].set_yscale('log')
    axs[0].set_ylim(*ylim)
    axs[0].set_ylabel('Conductance Values')
    axs[0].legend()

    # Second subplot
    axs[1].bar(x - width/2, f2_g_2, width, label='Initial g')
    axs[1].bar(x + width/2, f2_final_g_2, width, label='Final g')
    axs[1].set_title(titles[1])
    axs[1].set_xticks(x)
    axs[1].set_xticklabels(labels)
    axs[1].set_yscale('log')
    axs[1].set_ylim(*ylim)
    axs[1].legend()
    axs[0].grid(axis='y', linestyle='--', alpha=0.7)
    axs[1].grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()

import matplotlib.pyplot as plt
import numpy as np

def plot_voltage_traces(file_handles, index1, index2, duration=120000, figsize=(8, 2),
                        title="Final Voltage", ylim=None):
    """
    Plots voltage traces from two file handles side by side.

    Parameters:
    - file_handles: list of dicts with 'V' and 'tme' keys
    - index1, index2: indices of file_handles to compare
    - duration: number of samples to plot (default: 120000)
    - figsize: tuple for figure size
    - title: figure title
    - ylim: tuple (ymin, ymax) for voltage axis; if None, auto-scale
    """
    f1 = file_handles[index1]
    f2 = file_handles[index2]

    V1 = np.array(f1['V'])
    t1 = np.array(f1['tme']) / 60000  # ms to minutes

    V2 = np.array(f2['V'])
    t2 = np.array(f2['tme']) / 60000

    plt.rcParams['text.usetex'] = True

    fig, axes = plt.subplots(1, 2, figsize=figsize, sharex=False, sharey=False)

    axes[0].plot(t1[:duration], V1[:duration], 'k')
    axes[0].set_xlabel('Time (min)')
    axes[0].set_ylabel('Voltage (mV)')
    if ylim:
        axes[0].set_ylim(*ylim)

    axes[1].plot(t2[:duration], V2[:duration], 'k')
    axes[1].set_xlabel('Time (min)')
    axes[1].set_ylabel('Voltage (mV)')
    if ylim:
        axes[1].set_ylim(*ylim)

    plt.figtext(0.5, 1, title, size="large", ha='center')
    plt.tight_layout()
    # plt.subplots_adjust(hspace=0.6, wspace=0.4)
    plt.show()


import numpy as np

def process_simulation_data(file1, file2, F, S, D, PB_list_1, DB_list_1, PB_list_2, DB_list_2, extract_last_values_func):
    """
    Processes simulation data from two files, adjusts values based on F, S, D, and separates PB and DB entries.

    Parameters:
        file1 (str): Path to the first .txt file
        file2 (str): Path to the second .txt file
        F (float): Desired F value to adjust
        S (float): Desired S value to adjust
        D (float): Desired D value to adjust
        PB_list_1 (list of int): PB indices for file1
        DB_list_1 (list of int): DB indices for file1
        PB_list_2 (list of int): PB indices for file2
        DB_list_2 (list of int): DB indices for file2
        extract_last_values_func (function): Function to extract last values from a file

    Returns:
        tuple: (PB, DB) arrays, each with shape (n, 3) corresponding to adjusted [F, S, D] differences
    """
    lv1 = np.array(extract_last_values_func(file1))
    lv1[:, 16] = F - lv1[:, 16]
    lv1[:, 17] = S - lv1[:, 17]
    lv1[:, 18] = D - lv1[:, 18]

    lv2 = np.array(extract_last_values_func(file2))
    lv2[:, 16] = F - lv2[:, 16]
    lv2[:, 17] = S - lv2[:, 17]
    lv2[:, 18] = D - lv2[:, 18]

    DB_1 = lv1[DB_list_1, 16:19]
    PB_1 = lv1[PB_list_1, 16:19]
    DB_2 = lv2[DB_list_2, 16:19]
    PB_2 = lv2[PB_list_2, 16:19]

    DB = np.vstack((DB_1, DB_2))
    PB = np.vstack((PB_1, PB_2))

    return PB, DB



def plot_3d_target(
    DB_data,
    PB_data,
    indices,
    group_labels=('DB', 'PB'),
    conductance_labels=('gX', 'gY', 'gZ'),
    highlight_points=None,  # List of tuples: [(x, y, z, label, color), ...]
    axis_limits=None  # Dict: {'x': (min, max), 'y': (min, max), 'z': (min, max)}
):
    x_idx, y_idx, z_idx = indices
    x_label, y_label, z_label = conductance_labels

    df_DB = pd.DataFrame({
        x_label: DB_data[:, x_idx],
        y_label: DB_data[:, y_idx],
        z_label: DB_data[:, z_idx],
        'group': [group_labels[0]] * len(DB_data)
    })

    df_PB = pd.DataFrame({
        x_label: PB_data[:, x_idx],
        y_label: PB_data[:, y_idx],
        z_label: PB_data[:, z_idx],
        'group': [group_labels[1]] * len(PB_data)
    })

    df = pd.concat([df_DB, df_PB], ignore_index=True)

    fig = px.scatter_3d(df, x=x_label, y=y_label, z=z_label, color='group')

    # Add multiple highlight points with individual colors
    if highlight_points:
        for point in highlight_points:
            if len(point) == 5:
                x_val, y_val, z_val, label, color = point
            else:
                raise ValueError("Each highlight point must be a tuple of (x, y, z, label, color)")
            fig.add_trace(go.Scatter3d(
                x=[x_val],
                y=[y_val],
                z=[z_val],
                mode='markers+text',
                marker=dict(size=6, symbol='diamond', color=color),
                text=[label],
                textposition='top center',
                name=label,
                showlegend=True
            ))

    # Apply axis limits if provided
    if axis_limits:
        fig.update_layout(
            scene=dict(
                xaxis=dict(range=axis_limits.get('x')),
                yaxis=dict(range=axis_limits.get('y')),
                zaxis=dict(range=axis_limits.get('z'))
            )
        )

    fig.show()

