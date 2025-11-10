import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from scipy.interpolate import interp1d
from tmm import coh_tmm  # For the TMM method
from berreman import StackModel  # For the Berreman 4x4 method
import os

###############################################################################
# 1. PEDOT Refractive Index Data Handler
###############################################################################
class PEDOTRefractiveIndex:
    """Class to handle PEDOT:PSS refractive index data from Excel file with voltage dependency"""
    
    def __init__(self, excel_file='300rpm extinction coefficient.xlsx'):
        self.wavelengths = None
        self.voltages = None
        self.k_matrix = None
        self.voltage_interpolators = {}
        self.load_data(excel_file)
    
    def load_data(self, excel_file):
        """Load voltage-dependent extinction coefficient data from Excel file"""
        try:
            excel_data = pd.read_excel(excel_file, sheet_name=0, header=None)
            
            print("Raw data shape:", excel_data.shape)
            
            # Extract voltage headers from row 0
            voltage_row = excel_data.iloc[0, 1:].values
            
            # Extract wavelengths from column A, starting from row 2
            self.wavelengths = pd.to_numeric(excel_data.iloc[2:, 0], errors='coerce').values
            valid_wave_indices = ~np.isnan(self.wavelengths)
            self.wavelengths = self.wavelengths[valid_wave_indices]
            
            # Process voltage columns
            self.voltages = []
            k_data = []
            
            for i, voltage_header in enumerate(voltage_row):
                if pd.isna(voltage_header):
                    continue
                    
                voltage_str = str(voltage_header).replace('mV', '').replace('mv', '').strip()
                
                try:
                    voltage = float(voltage_str)
                    self.voltages.append(voltage)
                    
                    k_values = pd.to_numeric(excel_data.iloc[2:, i+1], errors='coerce').values
                    k_values = k_values[valid_wave_indices]
                    k_data.append(k_values)
                    
                except ValueError:
                    continue
            
            if len(self.voltages) == 0:
                print("No voltage columns found! Using default values.")
                self.voltages = np.array([0.0])
                k_data = [np.full_like(self.wavelengths, 0.05)]
            
            self.voltages = np.array(self.voltages)
            self.k_matrix = np.column_stack(k_data)
            
            # Create interpolation functions for each voltage
            for i, voltage in enumerate(self.voltages):
                k_values = self.k_matrix[:, i]
                valid_k_indices = ~np.isnan(k_values)
                if np.sum(valid_k_indices) > 1:
                    self.voltage_interpolators[voltage] = interp1d(
                        self.wavelengths[valid_k_indices], 
                        k_values[valid_k_indices], 
                        kind='linear',
                        bounds_error=False,
                        fill_value='extrapolate'
                    )
            
            print(f"Loaded PEDOT:PSS extinction data: {len(self.wavelengths)} wavelength points")
            print(f"Wavelength range: {self.wavelengths.min():.1f} - {self.wavelengths.max():.1f} nm")
            print(f"Available voltages: {sorted(self.voltages)} mV")
            
        except Exception as e:
            print(f"Error loading Excel file: {e}")
            # Initialize with default values
            self.wavelengths = np.linspace(400, 800, 100)
            self.voltages = np.array([0.0])
            self.k_matrix = np.full((len(self.wavelengths), 1), 0.05)
            self.voltage_interpolators = {}
    
    def get_k(self, wavelength, voltage):
        """Get extinction coefficient at specific wavelength and voltage"""
        if voltage in self.voltage_interpolators:
            return float(self.voltage_interpolators[voltage](wavelength))
        else:
            if len(self.voltages) < 2:
                return 0.05
            
            voltage_diff = np.abs(self.voltages - voltage)
            closest_idx = np.argsort(voltage_diff)
            
            if voltage_diff[closest_idx[0]] == 0:
                return float(self.voltage_interpolators[self.voltages[closest_idx[0]]](wavelength))
            
            v1, v2 = self.voltages[closest_idx[0]], self.voltages[closest_idx[1]]
            k1 = self.voltage_interpolators[v1](wavelength)
            k2 = self.voltage_interpolators[v2](wavelength)
            
            k_interp = k1 + (k2 - k1) * (voltage - v1) / (v2 - v1)
            return float(k_interp)
    
    def get_available_voltages(self):
        """Get list of available voltages"""
        return sorted(self.voltages)

# Initialize PEDOT refractive index data
pedot_data = PEDOTRefractiveIndex()

###############################################################################
# 2. Reflectivity Calculation Functions
###############################################################################
def calculate_reflectivity_vs_angle_tmm(wavelength, voltage, n_real, layer_thicknesses, 
                                        polarization='s', delta_k=0, manual_k=None):
    """
    Calculate reflectivity using TMM (isotropic).
    """
    degree = np.pi / 180
    
    n_bk = 1.518
    n_ito = 1.8529 + 0.00316*1j
    n_water = 1.333
    
    if manual_k is not None:
        k_pedot = manual_k + delta_k
    else:
        k_pedot = pedot_data.get_k(wavelength, voltage) + delta_k
    n_pedot = n_real + 1j * k_pedot
    
    n_list = [n_bk, n_ito, n_pedot, n_water]
    d_list = [np.inf, layer_thicknesses['ito'], layer_thicknesses['pedot'], np.inf]
    
    angles = np.linspace(20, 89.9, 500)
    reflectivities = []
    
    for angle in angles:
        result = coh_tmm(polarization, n_list, d_list, angle * degree, wavelength)
        reflectivities.append(result['R'])
    
    return angles, np.array(reflectivities)

def calculate_reflectivity_vs_angle_berreman(wavelength, voltage, n_real_xy, k_xy, n_real_z, k_z, 
                                             layer_thicknesses, polarization='s', delta_k=0):
    """
    Calculate reflectivity using Berreman 4x4 matrix method (anisotropic).
    
    Parameters:
    -----------
    n_real_xy, k_xy : float
        In-plane (ordinary, TE) refractive index: n_xy = n_real_xy + k_xy*j
    n_real_z, k_z : float
        Out-of-plane (extraordinary, TM) refractive index: n_z = n_real_z + k_z*j
    """
    degree = np.pi / 180
    
    n_bk = 1.518
    n_ito = 1.8529 + 0.00316*1j
    n_water = 1.333
    
    # Apply delta_k
    k_xy_modified = k_xy + delta_k
    k_z_modified = k_z + delta_k
    
    # Calculate permittivity from n and k: ε = (n² - k²) + 2nk*j
    eps_xy_real = n_real_xy**2 - k_xy_modified**2
    eps_xy_imag = 2 * n_real_xy * k_xy_modified
    eps_xy = eps_xy_real + eps_xy_imag * 1j
    
    eps_z_real = n_real_z**2 - k_z_modified**2
    eps_z_imag = 2 * n_real_z * k_z_modified
    eps_z = eps_z_real + eps_z_imag * 1j
    
    angles = np.linspace(20, 89.9, 500)
    reflectivities = []
    
    for angle in angles:
        angle_rad = angle * degree
        
        # Extract ITO refractive index components
        n_ito_real = 1.8529
        k_ito = 0.00316
        
        # Calculate ITO permittivity explicitly as: ε = (n² - k²) + 2nk*j
        eps_ito_real = n_ito_real**2 - k_ito**2
        eps_ito_imag = 2 * n_ito_real * k_ito
        eps_ito_value = eps_ito_real + eps_ito_imag * 1j
        
        # Create ITO permittivity tensor
        eps_ito = np.array([
            [eps_ito_value, 0, 0],
            [0, eps_ito_value, 0],
            [0, 0, eps_ito_value]
        ])
        
        # Anisotropic PEDOT:PSS
        eps_pedot = np.array([
            [eps_xy, 0, 0],
            [0, eps_xy, 0],
            [0, 0, eps_z]
        ])
        
        # Create stack model
        model = StackModel(
            eps_list=[eps_ito, eps_pedot],
            thickness_nm_list=[layer_thicknesses['ito'], layer_thicknesses['pedot']],
            n_entry=n_bk,
            n_exit=n_water,
            wl_nm=wavelength,
            theta_in_rad=angle_rad
        )
        
        try:
            reflectance, _ = model.get_refl_trans(method="SM")
            
            if polarization == 'p':
                R = reflectance[0, 0]  # R_p_to_p (TM)
            else:  # 's'
                R = reflectance[1, 1]  # R_s_to_s (TE)
                
            reflectivities.append(R)
        except Exception as e:
            print(f"Error at angle {angle}°: {str(e)}")
            reflectivities.append(0)
    
    return angles, np.array(reflectivities)

###############################################################################
# 3. Plotting Functions
###############################################################################
def plot_combined_te_tm_tmm(wavelength, voltage, n_real_te, k_te, n_real_tm, k_tm, 
                           layer_thicknesses, shift_te=0, shift_tm=0):
    """
    Create combined TE/TM plot using TMM method.
    """
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    # Calculate TE mode (s-polarization)
    angles_te, R_te = calculate_reflectivity_vs_angle_tmm(
        wavelength, voltage, n_real_te, layer_thicknesses, 
        polarization='s', manual_k=k_te
    )
    
    # Calculate TM mode (p-polarization)
    angles_tm, R_tm = calculate_reflectivity_vs_angle_tmm(
        wavelength, voltage, n_real_tm, layer_thicknesses, 
        polarization='p', manual_k=k_tm
    )
    
    # Plot simulations
    ax.plot(angles_te, R_te, 'b-', linewidth=2, 
            label=f'TE mode (TMM: n={n_real_te:.3f}+{k_te:.4f}i)')
    ax.plot(angles_tm, R_tm, 'r-', linewidth=2, 
            label=f'TM mode (TMM: n={n_real_tm:.3f}+{k_tm:.4f}i)')
    
    # Plot TE experimental data
    try:
        if os.path.exists('eploy_22s_Spol.xlsx'):
            df_te = pd.read_excel('eploy_22s_Spol.xlsx')
            if len(df_te.columns) >= 2:
                external_angles_te = df_te.iloc[:, 0].values
                external_R_te = df_te.iloc[:, 1].values
                mapped_angles_te = -1 * external_angles_te + 273 + shift_te
                sort_idx = np.argsort(mapped_angles_te)
                mapped_angles_te = mapped_angles_te[sort_idx]
                external_R_te = external_R_te[sort_idx]
                
                shift_info = f" (shift: {shift_te:.1f}°)" if shift_te != 0 else ""
                ax.plot(mapped_angles_te, external_R_te, 'go', linestyle='-', linewidth=2, 
                       markersize=4, label=f'TE exp{shift_info}')
    except Exception as e:
        print(f"Error reading TE experimental data: {e}")
    
    # Plot TM experimental data
    try:
        if os.path.exists('eploy_22s_Ppol.xlsx'):
            df_tm = pd.read_excel('eploy_22s_Ppol.xlsx')
            if len(df_tm.columns) >= 2:
                external_angles_tm = df_tm.iloc[:, 0].values
                external_R_tm = df_tm.iloc[:, 1].values
                mapped_angles_tm = -1 * external_angles_tm + 273 + shift_tm
                sort_idx = np.argsort(mapped_angles_tm)
                mapped_angles_tm = mapped_angles_tm[sort_idx]
                external_R_tm = external_R_tm[sort_idx]
                
                shift_info = f" (shift: {shift_tm:.1f}°)" if shift_tm != 0 else ""
                ax.plot(mapped_angles_tm, external_R_tm, 'ms', linestyle='-', linewidth=2, 
                       markersize=4, label=f'TM exp{shift_info}')
    except Exception as e:
        print(f"Error reading TM experimental data: {e}")
    
    ito_thickness = layer_thicknesses.get('ito', 0)
    pedot_thickness = layer_thicknesses.get('pedot', 0)
    
    ax.set_xlabel('Incident Angle (degrees)', fontsize=12)
    ax.set_ylabel('Reflectivity', fontsize=12)
    
    title = f'TMM (Isotropic): Combined TE/TM Reflectivity\n'
    title += f'λ = {wavelength} nm, V = {voltage} mV\n'
    title += f'ITO = {ito_thickness} nm, PEDOT = {pedot_thickness} nm\n'
    title += f'TE: n = {n_real_te:.3f}+{k_te:.4f}i, TM: n = {n_real_tm:.3f}+{k_tm:.4f}i'
    
    ax.set_title(title, fontsize=14)
    ax.set_xlim(20, 90)
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3)
    ax.set_xticks(np.arange(20, 91, 10))
    ax.legend(loc='best', framealpha=0.9, fontsize=9)
    
    plt.tight_layout()
    return fig

def plot_combined_te_tm_berreman(wavelength, voltage, n_te, k_te, n_tm, k_tm, 
                                 layer_thicknesses, shift_te=0, shift_tm=0):
    """
    Create combined TE/TM plot using Berreman method.
    """
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    # Calculate TE mode (s-polarization)
    angles_te, R_te = calculate_reflectivity_vs_angle_berreman(
        wavelength, voltage, n_te, k_te, n_tm, k_tm,
        layer_thicknesses, polarization='s', delta_k=0
    )
    
    # Calculate TM mode (p-polarization)
    angles_tm, R_tm = calculate_reflectivity_vs_angle_berreman(
        wavelength, voltage, n_te, k_te, n_tm, k_tm,
        layer_thicknesses, polarization='p', delta_k=0
    )
    
    # Plot simulations
    ax.plot(angles_te, R_te, 'b-', linewidth=2, 
            label=f'TE mode (Berreman: n_xy={n_te:.3f}+{k_te:.4f}i)')
    ax.plot(angles_tm, R_tm, 'r-', linewidth=2, 
            label=f'TM mode (Berreman: n_z={n_tm:.3f}+{k_tm:.4f}i)')
    
    # Plot TE experimental data
    try:
        if os.path.exists('eploy_22s_Spol.xlsx'):
            df_te = pd.read_excel('eploy_22s_Spol.xlsx')
            if len(df_te.columns) >= 2:
                external_angles_te = df_te.iloc[:, 0].values
                external_R_te = df_te.iloc[:, 1].values
                mapped_angles_te = -1 * external_angles_te + 273 + shift_te
                sort_idx = np.argsort(mapped_angles_te)
                mapped_angles_te = mapped_angles_te[sort_idx]
                external_R_te = external_R_te[sort_idx]
                
                shift_info = f" (shift: {shift_te:.1f}°)" if shift_te != 0 else ""
                ax.plot(mapped_angles_te, external_R_te, 'go', linestyle='-', linewidth=2, 
                       markersize=4, label=f'TE exp{shift_info}')
    except Exception as e:
        print(f"Error reading TE experimental data: {e}")
    
    # Plot TM experimental data
    try:
        if os.path.exists('eploy_22s_Ppol.xlsx'):
            df_tm = pd.read_excel('eploy_22s_Ppol.xlsx')
            if len(df_tm.columns) >= 2:
                external_angles_tm = df_tm.iloc[:, 0].values
                external_R_tm = df_tm.iloc[:, 1].values
                mapped_angles_tm = -1 * external_angles_tm + 273 + shift_tm
                sort_idx = np.argsort(mapped_angles_tm)
                mapped_angles_tm = mapped_angles_tm[sort_idx]
                external_R_tm = external_R_tm[sort_idx]
                
                shift_info = f" (shift: {shift_tm:.1f}°)" if shift_tm != 0 else ""
                ax.plot(mapped_angles_tm, external_R_tm, 'ms', linestyle='-', linewidth=2, 
                       markersize=4, label=f'TM exp{shift_info}')
    except Exception as e:
        print(f"Error reading TM experimental data: {e}")
    
    ito_thickness = layer_thicknesses.get('ito', 0)
    pedot_thickness = layer_thicknesses.get('pedot', 0)
    
    ax.set_xlabel('Incident Angle (degrees)', fontsize=12)
    ax.set_ylabel('Reflectivity', fontsize=12)
    
    title = f'Berreman 4x4 (Anisotropic): Combined TE/TM Reflectivity\n'
    title += f'λ = {wavelength} nm, V = {voltage} mV\n'
    title += f'ITO = {ito_thickness} nm, PEDOT = {pedot_thickness} nm\n'
    title += f'TE: n_xy = {n_te:.3f}+{k_te:.4f}i, TM: n_z = {n_tm:.3f}+{k_tm:.4f}i'
    
    ax.set_title(title, fontsize=14)
    ax.set_xlim(20, 90)
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3)
    ax.set_xticks(np.arange(20, 91, 10))
    ax.legend(loc='best', framealpha=0.9, fontsize=9)
    
    plt.tight_layout()
    return fig

def plot_all_four_combined(wavelength, voltage, 
                          n_tmm_te, k_tmm_te, n_tmm_tm, k_tmm_tm,
                          n_te, k_te, n_tm, k_tm,
                          layer_thicknesses, 
                          shift_tmm_te=0, shift_tmm_tm=0,
                          shift_te=0, shift_tm=0):
    """
    Create combined plot showing all 4 curves: TE-TMM, TM-TMM, TE-Berreman, TM-Berreman.
    Similar to 8screens4_2.py combined mode.
    """
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    
    # Calculate TMM - TE mode (s-polarization)
    angles_tmm_te, R_tmm_te = calculate_reflectivity_vs_angle_tmm(
        wavelength, voltage, n_tmm_te, layer_thicknesses, 
        polarization='s', manual_k=k_tmm_te
    )
    
    # Calculate TMM - TM mode (p-polarization)
    angles_tmm_tm, R_tmm_tm = calculate_reflectivity_vs_angle_tmm(
        wavelength, voltage, n_tmm_tm, layer_thicknesses, 
        polarization='p', manual_k=k_tmm_tm
    )
    
    # Calculate Berreman - TE mode (s-polarization)
    angles_brm_te, R_brm_te = calculate_reflectivity_vs_angle_berreman(
        wavelength, voltage, n_te, k_te, n_tm, k_tm,
        layer_thicknesses, polarization='s', delta_k=0
    )
    
    # Calculate Berreman - TM mode (p-polarization)
    angles_brm_tm, R_brm_tm = calculate_reflectivity_vs_angle_berreman(
        wavelength, voltage, n_te, k_te, n_tm, k_tm,
        layer_thicknesses, polarization='p', delta_k=0
    )
    
    # Plot all 4 simulation curves
    ax.plot(angles_tmm_te, R_tmm_te, 'r-o', markevery=20, markersize=4, linewidth=1.5, 
            label=f'TE-TMM (n={n_tmm_te:.3f}+{k_tmm_te:.4f}i)', zorder=10)
    ax.plot(angles_brm_te, R_brm_te, 'b-s', markevery=20, markersize=4, linewidth=1.5, 
            label=f'TE-Berreman (n_xy={n_te:.3f}+{k_te:.4f}i)')
    ax.plot(angles_brm_tm, R_brm_tm, 'g-^', markevery=20, markersize=4, linewidth=1.5, 
            label=f'TM-Berreman (n_z={n_tm:.3f}+{k_tm:.4f}i)')
    ax.plot(angles_tmm_tm, R_tmm_tm, 'm-d', markevery=20, markersize=4, linewidth=1.5, 
            label=f'TM-TMM (n={n_tmm_tm:.3f}+{k_tmm_tm:.4f}i)', zorder=5)
    
    # Plot TE experimental data
    try:
        if os.path.exists('eploy_22s_Spol.xlsx'):
            df_te = pd.read_excel('eploy_22s_Spol.xlsx')
            if len(df_te.columns) >= 2:
                external_angles_te = df_te.iloc[:, 0].values
                external_R_te = df_te.iloc[:, 1].values
                # Use average shift for TE experimental data
                shift_te_avg = (shift_tmm_te + shift_te) / 2
                mapped_angles_te = -1 * external_angles_te + 273 + shift_te_avg
                sort_idx = np.argsort(mapped_angles_te)
                mapped_angles_te = mapped_angles_te[sort_idx]
                external_R_te = external_R_te[sort_idx]
                
                shift_info = f" (shift: {shift_te_avg:.1f}°)" if shift_te_avg != 0 else ""
                ax.plot(mapped_angles_te, external_R_te, 'k*', markersize=6, 
                       label=f'TE exp{shift_info}', zorder=15)
    except Exception as e:
        print(f"Error reading TE experimental data: {e}")
    
    # Plot TM experimental data
    try:
        if os.path.exists('eploy_22s_Ppol.xlsx'):
            df_tm = pd.read_excel('eploy_22s_Ppol.xlsx')
            if len(df_tm.columns) >= 2:
                external_angles_tm = df_tm.iloc[:, 0].values
                external_R_tm = df_tm.iloc[:, 1].values
                # Use average shift for TM experimental data
                shift_tm_avg = (shift_tmm_tm + shift_tm) / 2
                mapped_angles_tm = -1 * external_angles_tm + 273 + shift_tm_avg
                sort_idx = np.argsort(mapped_angles_tm)
                mapped_angles_tm = mapped_angles_tm[sort_idx]
                external_R_tm = external_R_tm[sort_idx]
                
                shift_info = f" (shift: {shift_tm_avg:.1f}°)" if shift_tm_avg != 0 else ""
                ax.plot(mapped_angles_tm, external_R_tm, 'kx', markersize=6, 
                       label=f'TM exp{shift_info}', zorder=15)
    except Exception as e:
        print(f"Error reading TM experimental data: {e}")
    
    ito_thickness = layer_thicknesses.get('ito', 0)
    pedot_thickness = layer_thicknesses.get('pedot', 0)
    
    ax.set_xlabel('Incident Angle (degrees)', fontsize=12)
    ax.set_ylabel('Reflectivity', fontsize=12)
    
    title = f'Combined: All 4 Polarizations (TE-TMM, TE-Berreman, TM-Berreman, TM-TMM)\n'
    title += f'λ = {wavelength} nm, V = {voltage} mV\n'
    title += f'ITO = {ito_thickness} nm, PEDOT = {pedot_thickness} nm'
    
    ax.set_title(title, fontsize=14)
    ax.set_xlim(20, 90)
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3)
    ax.set_xticks(np.arange(20, 91, 10))
    ax.legend(loc='best', framealpha=0.9, fontsize=9)
    
    plt.tight_layout()
    return fig

###############################################################################
# 4. Tkinter GUI
###############################################################################
class ReflectivityGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("TMM & Berreman Reflectivity Calculator")
        self.root.geometry("600x1100")
        
        self.default_thicknesses = {
            'ito': 15,
            'pedot': 87
        }
        
        # Create main canvas and scrollbar
        canvas = tk.Canvas(root)
        scrollbar = ttk.Scrollbar(root, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        canvas.bind_all("<MouseWheel>", _on_mousewheel)
        
        main_frame = ttk.Frame(scrollable_frame, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Title
        title_label = ttk.Label(main_frame, text="TMM & Berreman 4x4 Reflectivity Calculator", 
                               font=('Arial', 14, 'bold'))
        title_label.pack(pady=10)
        
        # Layer thickness frame
        thickness_frame = ttk.LabelFrame(main_frame, text="Layer Thicknesses (nm)", padding="10")
        thickness_frame.pack(fill=tk.X, padx=10, pady=5)
        
        self.thickness_vars = {}
        layers = [
            ('ITO', 'ito', self.default_thicknesses['ito']),
            ('PEDOT:PSS', 'pedot', self.default_thicknesses['pedot'])
        ]
        
        for i, (label, key, default) in enumerate(layers):
            ttk.Label(thickness_frame, text=f"{label}:").grid(row=i, column=0, sticky='w', padx=5, pady=3)
            self.thickness_vars[key] = tk.DoubleVar(value=default)
            entry = ttk.Entry(thickness_frame, textvariable=self.thickness_vars[key], width=10)
            entry.grid(row=i, column=1, padx=5, pady=3)
        
        structure_label = ttk.Label(thickness_frame, 
                                  text="Structure: Glass | ITO | PEDOT:PSS | Water",
                                  font=('Arial', 9, 'italic'))
        structure_label.grid(row=len(layers), column=0, columnspan=2, pady=5)
        
        reset_btn = ttk.Button(thickness_frame, text="Reset to Defaults", 
                              command=self.reset_thicknesses)
        reset_btn.grid(row=len(layers)+1, column=0, columnspan=2, pady=10)
        
        # Wavelength settings
        wavelength_frame = ttk.LabelFrame(main_frame, text="Wavelength Settings", padding="10")
        wavelength_frame.pack(fill=tk.X, padx=10, pady=5)
        
        ttk.Label(wavelength_frame, text="Wavelength (nm):").grid(row=0, column=0, sticky='w', padx=5, pady=3)
        self.wavelength_var = tk.DoubleVar(value=561.0)
        ttk.Entry(wavelength_frame, textvariable=self.wavelength_var, width=10).grid(row=0, column=1, padx=5, pady=3)
        
        # Voltage settings
        voltage_frame = ttk.LabelFrame(main_frame, text="Voltage Settings", padding="10")
        voltage_frame.pack(fill=tk.X, padx=10, pady=5)
        
        ttk.Label(voltage_frame, text="Voltage (mV):").grid(row=0, column=0, sticky='w', padx=5, pady=3)
        self.voltage_var = tk.DoubleVar(value=0.0)
        ttk.Entry(voltage_frame, textvariable=self.voltage_var, width=10).grid(row=0, column=1, padx=5, pady=3)
        
        available_voltages = pedot_data.get_available_voltages()
        if available_voltages:
            voltage_info = f"Available: {min(available_voltages):.0f} to {max(available_voltages):.0f} mV"
        else:
            voltage_info = "No voltage data loaded"
        ttk.Label(voltage_frame, text=voltage_info, font=('Arial', 8, 'italic')).grid(row=1, column=0, columnspan=2, pady=2)
        
        # ==================== TMM SECTION ====================
        tmm_frame = ttk.LabelFrame(main_frame, text="TMM (Isotropic) Refractive Indices", padding="10")
        tmm_frame.pack(fill=tk.X, padx=10, pady=5)
        
        ttk.Label(tmm_frame, text="TMM uses single n and k for each mode (isotropic)", 
                 font=('Arial', 9, 'italic'), foreground='darkblue').grid(row=0, column=0, columnspan=4, pady=5)
        
        # TE Mode for TMM
        ttk.Label(tmm_frame, text="TE Mode:", font=('Arial', 9, 'bold')).grid(row=1, column=0, columnspan=2, sticky='w', pady=5)
        
        ttk.Label(tmm_frame, text="n (TE):").grid(row=2, column=0, sticky='w', padx=5, pady=2)
        self.n_tmm_te_var = tk.DoubleVar(value=1.41)
        ttk.Entry(tmm_frame, textvariable=self.n_tmm_te_var, width=10).grid(row=2, column=1, padx=5, pady=2)
        
        ttk.Label(tmm_frame, text="k (TE):").grid(row=3, column=0, sticky='w', padx=5, pady=2)
        self.k_tmm_te_var = tk.DoubleVar(value=0.05)
        ttk.Entry(tmm_frame, textvariable=self.k_tmm_te_var, width=10).grid(row=3, column=1, padx=5, pady=2)
        
        ttk.Label(tmm_frame, text="Shift (TE):").grid(row=4, column=0, sticky='w', padx=5, pady=2)
        self.shift_tmm_te_var = tk.DoubleVar(value=0.0)
        ttk.Entry(tmm_frame, textvariable=self.shift_tmm_te_var, width=10).grid(row=4, column=1, padx=5, pady=2)
        
        # TM Mode for TMM
        ttk.Label(tmm_frame, text="TM Mode:", font=('Arial', 9, 'bold')).grid(row=1, column=2, columnspan=2, sticky='w', pady=5)
        
        ttk.Label(tmm_frame, text="n (TM):").grid(row=2, column=2, sticky='w', padx=5, pady=2)
        self.n_tmm_tm_var = tk.DoubleVar(value=1.41)
        ttk.Entry(tmm_frame, textvariable=self.n_tmm_tm_var, width=10).grid(row=2, column=3, padx=5, pady=2)
        
        ttk.Label(tmm_frame, text="k (TM):").grid(row=3, column=2, sticky='w', padx=5, pady=2)
        self.k_tmm_tm_var = tk.DoubleVar(value=0.05)
        ttk.Entry(tmm_frame, textvariable=self.k_tmm_tm_var, width=10).grid(row=3, column=3, padx=5, pady=2)
        
        ttk.Label(tmm_frame, text="Shift (TM):").grid(row=4, column=2, sticky='w', padx=5, pady=2)
        self.shift_tmm_tm_var = tk.DoubleVar(value=0.0)
        ttk.Entry(tmm_frame, textvariable=self.shift_tmm_tm_var, width=10).grid(row=4, column=3, padx=5, pady=2)
        
        # ==================== BERREMAN SECTION ====================
        berreman_frame = ttk.LabelFrame(main_frame, text="Berreman 4x4 (Anisotropic) Refractive Indices", padding="10")
        berreman_frame.pack(fill=tk.X, padx=10, pady=5)
        
        ttk.Label(berreman_frame, text="Berreman uses separate n_xy, k_xy (in-plane) and n_z, k_z (out-of-plane)", 
                 font=('Arial', 9, 'italic'), foreground='darkgreen').grid(row=0, column=0, columnspan=4, pady=5)
        
        # TE values (in-plane, ordinary)
        ttk.Label(berreman_frame, text="In-plane (TE/ordinary):", font=('Arial', 9, 'bold')).grid(row=1, column=0, columnspan=2, sticky='w', pady=5)
        
        ttk.Label(berreman_frame, text="n_te (xy):").grid(row=2, column=0, sticky='w', padx=5, pady=2)
        self.n_te_var = tk.DoubleVar(value=1.41)
        ttk.Entry(berreman_frame, textvariable=self.n_te_var, width=10).grid(row=2, column=1, padx=5, pady=2)
        
        ttk.Label(berreman_frame, text="k_te (xy):").grid(row=3, column=0, sticky='w', padx=5, pady=2)
        self.k_te_var = tk.DoubleVar(value=0.05)
        ttk.Entry(berreman_frame, textvariable=self.k_te_var, width=10).grid(row=3, column=1, padx=5, pady=2)
        
        ttk.Label(berreman_frame, text="Shift (TE):").grid(row=4, column=0, sticky='w', padx=5, pady=2)
        self.shift_te_var = tk.DoubleVar(value=0.0)
        ttk.Entry(berreman_frame, textvariable=self.shift_te_var, width=10).grid(row=4, column=1, padx=5, pady=2)
        
        # TM values (out-of-plane, extraordinary)
        ttk.Label(berreman_frame, text="Out-of-plane (TM/extraordinary):", font=('Arial', 9, 'bold')).grid(row=1, column=2, columnspan=2, sticky='w', pady=5)
        
        ttk.Label(berreman_frame, text="n_tm (z):").grid(row=2, column=2, sticky='w', padx=5, pady=2)
        self.n_tm_var = tk.DoubleVar(value=1.269)
        ttk.Entry(berreman_frame, textvariable=self.n_tm_var, width=10).grid(row=2, column=3, padx=5, pady=2)
        
        ttk.Label(berreman_frame, text="k_tm (z):").grid(row=3, column=2, sticky='w', padx=5, pady=2)
        self.k_tm_var = tk.DoubleVar(value=0.05)
        ttk.Entry(berreman_frame, textvariable=self.k_tm_var, width=10).grid(row=3, column=3, padx=5, pady=2)
        
        ttk.Label(berreman_frame, text="Shift (TM):").grid(row=4, column=2, sticky='w', padx=5, pady=2)
        self.shift_tm_var = tk.DoubleVar(value=0.0)
        ttk.Entry(berreman_frame, textvariable=self.shift_tm_var, width=10).grid(row=4, column=3, padx=5, pady=2)
        
        # Experimental data checkbox
        exp_frame = ttk.LabelFrame(main_frame, text="Experimental Data", padding="10")
        exp_frame.pack(fill=tk.X, padx=10, pady=5)
        
        self.show_exp_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(exp_frame, text="Show experimental data (eploy_22s_Spol.xlsx, eploy_22s_Ppol.xlsx)", 
                       variable=self.show_exp_var).pack()
        
        # Buttons
        buttons_frame = ttk.Frame(main_frame)
        buttons_frame.pack(pady=10)
        
        tmm_btn = ttk.Button(buttons_frame, text="Plot TMM", 
                            command=self.plot_tmm, style='Accent.TButton')
        tmm_btn.grid(row=0, column=0, padx=5, pady=5)
        
        berreman_btn = ttk.Button(buttons_frame, text="Plot Berreman", 
                                 command=self.plot_berreman, style='Accent.TButton')
        berreman_btn.grid(row=0, column=1, padx=5, pady=5)
        
        combined_btn = ttk.Button(buttons_frame, text="📊 Plot All 4 Combined", 
                                 command=self.plot_all_combined, style='Accent.TButton')
        combined_btn.grid(row=0, column=2, padx=5, pady=5)
        
        # Status label
        self.status_label = ttk.Label(main_frame, text="Ready", foreground="green")
        self.status_label.pack()
        
        # Info section
        info_frame = ttk.LabelFrame(main_frame, text="Information", padding="10")
        info_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        info_text = ttk.Label(info_frame, text=
                             "Two Calculation Methods Available:\n\n"
                             "1. TMM (Transfer Matrix Method)\n"
                             "   - Treats PEDOT:PSS as isotropic material\n"
                             "   - Uses single n and k for each mode\n"
                             "   - Faster calculation\n\n"
                             "2. Berreman 4×4 Matrix Method\n"
                             "   - Treats PEDOT:PSS as anisotropic material\n"
                             "   - Uses n_te, k_te (in-plane/ordinary)\n"
                             "   - Uses n_tm, k_tm (out-of-plane/extraordinary)\n"
                             "   - More accurate for anisotropic materials\n\n"
                             "3. Combined Plot (All 4)\n"
                             "   - Shows TE-TMM, TM-TMM, TE-Berreman, TM-Berreman\n"
                             "   - Direct comparison of all methods on one graph\n"
                             "   - Similar to 8screens4_2.py combined mode\n\n"
                             "Layer structure: Glass | ITO | PEDOT | Water\n"
                             "Both methods use layer thickness for interference calculations.",
                            justify=tk.LEFT)
        info_text.pack(padx=5, pady=5)
    
    def reset_thicknesses(self):
        """Reset all thickness values to defaults"""
        for key, var in self.thickness_vars.items():
            var.set(self.default_thicknesses[key])
    
    def get_layer_thicknesses(self):
        """Get current layer thickness values"""
        return {key: var.get() for key, var in self.thickness_vars.items()}
    
    def plot_tmm(self):
        """Plot TMM reflectivity"""
        try:
            self.status_label.config(text="Calculating TMM reflectivity...", foreground="blue")
            self.root.update()
            
            # Get parameters
            wavelength = self.wavelength_var.get()
            voltage = self.voltage_var.get()
            n_te = self.n_tmm_te_var.get()
            k_te = self.k_tmm_te_var.get()
            n_tm = self.n_tmm_tm_var.get()
            k_tm = self.k_tmm_tm_var.get()
            shift_te = self.shift_tmm_te_var.get()
            shift_tm = self.shift_tmm_tm_var.get()
            layer_thicknesses = self.get_layer_thicknesses()
            
            # Create combined TE/TM plot
            fig = plot_combined_te_tm_tmm(
                wavelength, voltage, 
                n_te, k_te, n_tm, k_tm,
                layer_thicknesses,
                shift_te=shift_te,
                shift_tm=shift_tm
            )
            
            self.status_label.config(text="Ready", foreground="green")
            plt.show()
            
        except Exception as e:
            self.status_label.config(text=f"Error: {str(e)}", foreground="red")
            messagebox.showerror("Error", f"An error occurred: {str(e)}\n\n{e.__class__.__name__}")
    
    def plot_berreman(self):
        """Plot Berreman reflectivity"""
        try:
            self.status_label.config(text="Calculating Berreman reflectivity...", foreground="blue")
            self.root.update()
            
            # Get parameters
            wavelength = self.wavelength_var.get()
            voltage = self.voltage_var.get()
            n_te = self.n_te_var.get()
            k_te = self.k_te_var.get()
            n_tm = self.n_tm_var.get()
            k_tm = self.k_tm_var.get()
            shift_te = self.shift_te_var.get()
            shift_tm = self.shift_tm_var.get()
            layer_thicknesses = self.get_layer_thicknesses()
            
            # Create combined TE/TM plot
            fig = plot_combined_te_tm_berreman(
                wavelength, voltage, 
                n_te, k_te, n_tm, k_tm,
                layer_thicknesses,
                shift_te=shift_te,
                shift_tm=shift_tm
            )
            
            self.status_label.config(text="Ready", foreground="green")
            plt.show()
            
        except Exception as e:
            self.status_label.config(text=f"Error: {str(e)}", foreground="red")
            messagebox.showerror("Error", f"An error occurred: {str(e)}\n\n{e.__class__.__name__}")
    
    def plot_all_combined(self):
        """Plot all 4 combined: TE-TMM, TM-TMM, TE-Berreman, TM-Berreman"""
        try:
            self.status_label.config(text="Calculating all 4 methods...", foreground="blue")
            self.root.update()
            
            # Get parameters
            wavelength = self.wavelength_var.get()
            voltage = self.voltage_var.get()
            layer_thicknesses = self.get_layer_thicknesses()
            
            # TMM parameters
            n_tmm_te = self.n_tmm_te_var.get()
            k_tmm_te = self.k_tmm_te_var.get()
            n_tmm_tm = self.n_tmm_tm_var.get()
            k_tmm_tm = self.k_tmm_tm_var.get()
            shift_tmm_te = self.shift_tmm_te_var.get()
            shift_tmm_tm = self.shift_tmm_tm_var.get()
            
            # Berreman parameters
            n_te = self.n_te_var.get()
            k_te = self.k_te_var.get()
            n_tm = self.n_tm_var.get()
            k_tm = self.k_tm_var.get()
            shift_te = self.shift_te_var.get()
            shift_tm = self.shift_tm_var.get()
            
            # Create combined plot with all 4 curves
            fig = plot_all_four_combined(
                wavelength, voltage,
                n_tmm_te, k_tmm_te, n_tmm_tm, k_tmm_tm,
                n_te, k_te, n_tm, k_tm,
                layer_thicknesses,
                shift_tmm_te=shift_tmm_te,
                shift_tmm_tm=shift_tmm_tm,
                shift_te=shift_te,
                shift_tm=shift_tm
            )
            
            self.status_label.config(text="Ready", foreground="green")
            plt.show()
            
        except Exception as e:
            self.status_label.config(text=f"Error: {str(e)}", foreground="red")
            messagebox.showerror("Error", f"An error occurred: {str(e)}\n\n{e.__class__.__name__}")

# Run the GUI
if __name__ == "__main__":
    root = tk.Tk()
    app = ReflectivityGUI(root)
    root.mainloop()
