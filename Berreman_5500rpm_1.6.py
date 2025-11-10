import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk, messagebox
from scipy.interpolate import interp1d
from tmm import coh_tmm  # For the TMM method
from berreman import StackModel  # For the Berreman 4x4 method
import os

# Manual angle shift to add to all experimental data plots (in degrees)
MANUAL_ANGLE_SHIFT = 0  # Adjust if a global offset is required

###############################################################################
# 1. PEDOT Refractive Index Data Handler (5500rpm dataset)
###############################################################################
class PEDOTRefractiveIndex:
    """Class to handle PEDOT:PSS refractive index data from Excel file with voltage dependency"""
    
    def __init__(self, excel_file='300rpm extinction coefficient.xlsx'):
        self.wavelengths = None
        self.voltages = None
        self.k_matrix = None  # k values as function of wavelength and voltage
        self.voltage_interpolators = {}  # Dictionary of interpolators for each voltage
        self.load_data(excel_file)
    
    def load_data(self, excel_file):
        """Load voltage-dependent extinction coefficient data from Excel file"""
        try:
            excel_data = pd.read_excel(excel_file, sheet_name=0, header=None)
            
            print("Raw data shape:", excel_data.shape)
            print("First few rows:\n", excel_data.iloc[:3, :5])
            
            # Extract voltage headers from row 0 (first row)
            voltage_row = excel_data.iloc[0, 1:].values  # Skip first column (wavelength label)
            
            # Extract wavelengths from column A, starting from row 2 (index 2)
            self.wavelengths = pd.to_numeric(excel_data.iloc[2:, 0], errors='coerce').values
            valid_wave_indices = ~np.isnan(self.wavelengths)
            self.wavelengths = self.wavelengths[valid_wave_indices]
            
            # Process voltage columns
            self.voltages = []
            k_data = []
            
            print("Processing voltage columns:")
            for i, voltage_header in enumerate(voltage_row):
                if pd.isna(voltage_header):
                    continue
                
                voltage_str = str(voltage_header).replace('mV', '').replace('mv', '').strip()
                print(f"Header: {voltage_header} -> Voltage string: {voltage_str}")
                
                try:
                    voltage = float(voltage_str)
                    self.voltages.append(voltage)
                    
                    k_values = pd.to_numeric(excel_data.iloc[2:, i + 1], errors='coerce').values
                    k_values = k_values[valid_wave_indices]
                    k_data.append(k_values)
                    
                    print(
                        f"Processed {voltage} mV, k range: "
                        f"{np.nanmin(k_values):.6f} to {np.nanmax(k_values):.6f}"
                    )
                except ValueError:
                    print(f"Could not parse voltage from header: {voltage_header}")
                    continue
            
            if len(self.voltages) == 0:
                print("No voltage columns found! Using default values.")
                self.voltages = np.array([0.0])
                k_data = [np.full_like(self.wavelengths, 0.2)]  # Default k value
            
            self.voltages = np.array(self.voltages)
            self.k_matrix = np.column_stack(k_data)
            
            print(f"Loaded PEDOT:PSS extinction data: {len(self.wavelengths)} wavelength points")
            print(f"Wavelength range: {self.wavelengths.min():.1f} - {self.wavelengths.max():.1f} nm")
            print(f"Available voltages: {sorted(self.voltages)} mV")
            
            # Create interpolation functions for each voltage
            for i, voltage in enumerate(self.voltages):
                k_values = self.k_matrix[:, i]
                valid_k_indices = ~np.isnan(k_values)
                if np.sum(valid_k_indices) > 1:  # Need at least 2 points for interpolation
                    self.voltage_interpolators[voltage] = interp1d(
                        self.wavelengths[valid_k_indices],
                        k_values[valid_k_indices],
                        kind='linear',
                        bounds_error=False,
                        fill_value='extrapolate'
                    )
            
        except Exception as e:
            print(f"Error loading Excel file: {e}")
            # Initialize with default values
            self.wavelengths = np.linspace(400, 800, 100)
            self.voltages = np.array([0.0])
            self.k_matrix = np.full((len(self.wavelengths), 1), 0.2)
            self.voltage_interpolators = {}
    
    def get_k(self, wavelength, voltage):
        """Get extinction coefficient at specific wavelength and voltage"""
        if voltage in self.voltage_interpolators:
            return float(self.voltage_interpolators[voltage](wavelength))
        else:
            if len(self.voltages) < 2:
                return 0.2  # Default fallback
            
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
    """Calculate reflectivity using TMM (isotropic)."""
    degree = np.pi / 180
    
    n_bk = 1.518
    n_ito = 1.8529 + 0.00316 * 1j
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
    """Calculate reflectivity using Berreman 4x4 matrix method (anisotropic)."""
    degree = np.pi / 180
    
    n_bk = 1.518
    n_ito = 1.8529 + 0.00316 * 1j
    n_water = 1.333
    
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
        
        # ITO permittivity tensor (isotropic)
        n_ito_real = 1.8529
        k_ito = 0.00316
        eps_ito_real = n_ito_real**2 - k_ito**2
        eps_ito_imag = 2 * n_ito_real * k_ito
        eps_ito_value = eps_ito_real + eps_ito_imag * 1j
        eps_ito = np.diag([eps_ito_value, eps_ito_value, eps_ito_value])
        
        # PEDOT:PSS permittivity tensor (anisotropic)
        eps_pedot = np.array([
            [eps_xy, 0, 0],
            [0, eps_xy, 0],
            [0, 0, eps_z]
        ])
        
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
                R = reflectance[0, 0]  # TM
            else:
                R = reflectance[1, 1]  # TE
            reflectivities.append(R)
        except Exception as e:
            print(f"Berreman error at {angle:.2f}°: {e}")
            reflectivities.append(0)
    
    return angles, np.array(reflectivities)

###############################################################################
# 3. Plotting Functions
###############################################################################
def plot_combined_te_tm_tmm(wavelength, voltage, n_real_te, k_te, n_real_tm, k_tm,
                            layer_thicknesses, shift_te=0, shift_tm=0):
    """Create combined TE/TM plot using TMM method."""
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    angles_te, R_te = calculate_reflectivity_vs_angle_tmm(
        wavelength, voltage, n_real_te, layer_thicknesses,
        polarization='s', manual_k=k_te
    )
    angles_tm, R_tm = calculate_reflectivity_vs_angle_tmm(
        wavelength, voltage, n_real_tm, layer_thicknesses,
        polarization='p', manual_k=k_tm
    )
    
    ax.plot(angles_te, R_te, 'b-', linewidth=2,
            label=f'TE mode (TMM: n={n_real_te:.3f}+{k_te:.4f}i)')
    ax.plot(angles_tm, R_tm, 'r-', linewidth=2,
            label=f'TM mode (TMM: n={n_real_tm:.3f}+{k_tm:.4f}i)')
    
    try:
        if os.path.exists('5500rpm_1.6_s.xlsx'):
            df_te = pd.read_excel('5500rpm_1.6_s.xlsx')
            if len(df_te.columns) >= 2:
                external_angles_te = df_te.iloc[:, 0].values
                external_R_te = df_te.iloc[:, 1].values
                mapped_angles_te = -external_angles_te + 273 + shift_te + MANUAL_ANGLE_SHIFT
                sort_idx = np.argsort(mapped_angles_te)
                ax.plot(mapped_angles_te[sort_idx], external_R_te[sort_idx], 'go-', linewidth=2,
                        markersize=4, label=f'TE exp' + (f" (shift: {shift_te:.1f}°)" if shift_te else ""))
    except Exception as e:
        print(f"Error reading TE experimental data: {e}")
    
    try:
        if os.path.exists('5500rpm_1.6_p.xlsx'):
            df_tm = pd.read_excel('5500rpm_1.6_p.xlsx')
            if len(df_tm.columns) >= 2:
                external_angles_tm = df_tm.iloc[:, 0].values
                external_R_tm = df_tm.iloc[:, 1].values
                mapped_angles_tm = -external_angles_tm + 273 + shift_tm + MANUAL_ANGLE_SHIFT
                sort_idx = np.argsort(mapped_angles_tm)
                ax.plot(mapped_angles_tm[sort_idx], external_R_tm[sort_idx], 'ms-', linewidth=2,
                        markersize=4, label=f'TM exp' + (f" (shift: {shift_tm:.1f}°)" if shift_tm else ""))
    except Exception as e:
        print(f"Error reading TM experimental data: {e}")
    
    ito_thickness = layer_thicknesses.get('ito', 0)
    pedot_thickness = layer_thicknesses.get('pedot', 0)
    
    ax.set_xlabel('Incident Angle (degrees)')
    ax.set_ylabel('Reflectivity')
    ax.set_title(
        'TMM (Isotropic): Combined TE/TM Reflectivity\n'
        f'λ = {wavelength} nm, V = {voltage} mV\n'
        f'ITO = {ito_thickness} nm, PEDOT = {pedot_thickness} nm\n'
        f'TE: n = {n_real_te:.3f}+{k_te:.4f}i, TM: n = {n_real_tm:.3f}+{k_tm:.4f}i'
    )
    ax.set_xlim(20, 90)
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3)
    ax.set_xticks(np.arange(20, 91, 10))
    ax.legend(loc='best', framealpha=0.9, fontsize=9)
    
    plt.tight_layout()
    return fig

def plot_combined_te_tm_berreman(wavelength, voltage, n_te, k_te, n_tm, k_tm,
                                 layer_thicknesses, shift_te=0, shift_tm=0):
    """Create combined TE/TM plot using Berreman method."""
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    angles_te, R_te = calculate_reflectivity_vs_angle_berreman(
        wavelength, voltage, n_te, k_te, n_tm, k_tm,
        layer_thicknesses, polarization='s', delta_k=0
    )
    angles_tm, R_tm = calculate_reflectivity_vs_angle_berreman(
        wavelength, voltage, n_te, k_te, n_tm, k_tm,
        layer_thicknesses, polarization='p', delta_k=0
    )
    
    ax.plot(angles_te, R_te, 'b-', linewidth=2,
            label=f'TE mode (Berreman: n_xy={n_te:.3f}+{k_te:.4f}i)')
    ax.plot(angles_tm, R_tm, 'r-', linewidth=2,
            label=f'TM mode (Berreman: n_z={n_tm:.3f}+{k_tm:.4f}i)')
    
    try:
        if os.path.exists('5500rpm_1.6_s.xlsx'):
            df_te = pd.read_excel('5500rpm_1.6_s.xlsx')
            if len(df_te.columns) >= 2:
                external_angles_te = df_te.iloc[:, 0].values
                external_R_te = df_te.iloc[:, 1].values
                mapped_angles_te = -external_angles_te + 273 + shift_te + MANUAL_ANGLE_SHIFT
                sort_idx = np.argsort(mapped_angles_te)
                ax.plot(mapped_angles_te[sort_idx], external_R_te[sort_idx], 'go-', linewidth=2,
                        markersize=4, label=f'TE exp' + (f" (shift: {shift_te:.1f}°)" if shift_te else ""))
    except Exception as e:
        print(f"Error reading TE experimental data: {e}")
    
    try:
        if os.path.exists('5500rpm_1.6_p.xlsx'):
            df_tm = pd.read_excel('5500rpm_1.6_p.xlsx')
            if len(df_tm.columns) >= 2:
                external_angles_tm = df_tm.iloc[:, 0].values
                external_R_tm = df_tm.iloc[:, 1].values
                mapped_angles_tm = -external_angles_tm + 273 + shift_tm + MANUAL_ANGLE_SHIFT
                sort_idx = np.argsort(mapped_angles_tm)
                ax.plot(mapped_angles_tm[sort_idx], external_R_tm[sort_idx], 'ms-', linewidth=2,
                        markersize=4, label=f'TM exp' + (f" (shift: {shift_tm:.1f}°)" if shift_tm else ""))
    except Exception as e:
        print(f"Error reading TM experimental data: {e}")
    
    ito_thickness = layer_thicknesses.get('ito', 0)
    pedot_thickness = layer_thicknesses.get('pedot', 0)
    
    ax.set_xlabel('Incident Angle (degrees)')
    ax.set_ylabel('Reflectivity')
    ax.set_title(
        'Berreman 4x4 (Anisotropic): Combined TE/TM Reflectivity\n'
        f'λ = {wavelength} nm, V = {voltage} mV\n'
        f'ITO = {ito_thickness} nm, PEDOT = {pedot_thickness} nm\n'
        f'TE: n_xy = {n_te:.3f}+{k_te:.4f}i, TM: n_z = {n_tm:.3f}+{k_tm:.4f}i'
    )
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
    """Create combined plot showing all 4 curves: TE-TMM, TM-TMM, TE-Berreman, TM-Berreman."""
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    
    angles_tmm_te, R_tmm_te = calculate_reflectivity_vs_angle_tmm(
        wavelength, voltage, n_tmm_te, layer_thicknesses,
        polarization='s', manual_k=k_tmm_te
    )
    angles_tmm_tm, R_tmm_tm = calculate_reflectivity_vs_angle_tmm(
        wavelength, voltage, n_tmm_tm, layer_thicknesses,
        polarization='p', manual_k=k_tmm_tm
    )
    angles_brm_te, R_brm_te = calculate_reflectivity_vs_angle_berreman(
        wavelength, voltage, n_te, k_te, n_tm, k_tm,
        layer_thicknesses, polarization='s', delta_k=0
    )
    angles_brm_tm, R_brm_tm = calculate_reflectivity_vs_angle_berreman(
        wavelength, voltage, n_te, k_te, n_tm, k_tm,
        layer_thicknesses, polarization='p', delta_k=0
    )
    
    ax.plot(angles_tmm_te, R_tmm_te, 'r-o', markevery=20, markersize=4, linewidth=1.5,
            label=f'TE-TMM (n={n_tmm_te:.3f}+{k_tmm_te:.4f}i)', zorder=10)
    ax.plot(angles_brm_te, R_brm_te, 'b-s', markevery=20, markersize=4, linewidth=1.5,
            label=f'TE-Berreman (n_xy={n_te:.3f}+{k_te:.4f}i)')
    ax.plot(angles_brm_tm, R_brm_tm, 'g-^', markevery=20, markersize=4, linewidth=1.5,
            label=f'TM-Berreman (n_z={n_tm:.3f}+{k_tm:.4f}i)')
    ax.plot(angles_tmm_tm, R_tmm_tm, 'm-d', markevery=20, markersize=4, linewidth=1.5,
            label=f'TM-TMM (n={n_tmm_tm:.3f}+{k_tmm_tm:.4f}i)', zorder=5)
    
    try:
        if os.path.exists('5500rpm_1.6_s.xlsx'):
            df_te = pd.read_excel('5500rpm_1.6_s.xlsx')
            if len(df_te.columns) >= 2:
                external_angles_te = df_te.iloc[:, 0].values
                external_R_te = df_te.iloc[:, 1].values
                shift_te_avg = (shift_tmm_te + shift_te) / 2
                mapped_angles_te = -external_angles_te + 273 + shift_te_avg + MANUAL_ANGLE_SHIFT
                sort_idx = np.argsort(mapped_angles_te)
                ax.plot(mapped_angles_te[sort_idx], external_R_te[sort_idx], 'k*', markersize=6,
                        label=f'TE exp' + (f" (shift: {shift_te_avg:.1f}°)" if shift_te_avg else ""), zorder=15)
    except Exception as e:
        print(f"Error reading TE experimental data: {e}")
    
    try:
        if os.path.exists('5500rpm_1.6_p.xlsx'):
            df_tm = pd.read_excel('5500rpm_1.6_p.xlsx')
            if len(df_tm.columns) >= 2:
                external_angles_tm = df_tm.iloc[:, 0].values
                external_R_tm = df_tm.iloc[:, 1].values
                shift_tm_avg = (shift_tmm_tm + shift_tm) / 2
                mapped_angles_tm = -external_angles_tm + 273 + shift_tm_avg + MANUAL_ANGLE_SHIFT
                sort_idx = np.argsort(mapped_angles_tm)
                ax.plot(mapped_angles_tm[sort_idx], external_R_tm[sort_idx], 'kx', markersize=6,
                        label=f'TM exp' + (f" (shift: {shift_tm_avg:.1f}°)" if shift_tm_avg else ""), zorder=15)
    except Exception as e:
        print(f"Error reading TM experimental data: {e}")
    
    ito_thickness = layer_thicknesses.get('ito', 0)
    pedot_thickness = layer_thicknesses.get('pedot', 0)
    
    ax.set_xlabel('Incident Angle (degrees)')
    ax.set_ylabel('Reflectivity')
    ax.set_title(
        'Combined: TE-TMM, TE-Berreman, TM-Berreman, TM-TMM\n'
        f'λ = {wavelength} nm, V = {voltage} mV\n'
        f'ITO = {ito_thickness} nm, PEDOT = {pedot_thickness} nm'
    )
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
        self.root.title("TMM & Berreman Reflectivity Calculator (5500rpm 1.6)")
        self.root.geometry("600x1100")
        
        self.default_thicknesses = {
            'ito': 15,
            'pedot': 30.9
        }
        
        # Create main canvas and scrollbar for long layout
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
            canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
        canvas.bind_all("<MouseWheel>", _on_mousewheel)
        
        main_frame = ttk.Frame(scrollable_frame, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        ttk.Label(main_frame,
                  text="TMM & Berreman 4x4 Reflectivity Calculator",
                  font=('Arial', 14, 'bold')).pack(pady=10)
        
        # Layer thickness frame
        thickness_frame = ttk.LabelFrame(main_frame, text="Layer Thicknesses (nm)", padding="10")
        thickness_frame.pack(fill=tk.X, padx=10, pady=5)
        
        self.thickness_vars = {}
        for i, (label, key, default) in enumerate([
            ('ITO', 'ito', self.default_thicknesses['ito']),
            ('PEDOT:PSS', 'pedot', self.default_thicknesses['pedot'])
        ]):
            ttk.Label(thickness_frame, text=f"{label}:").grid(row=i, column=0, sticky='w', padx=5, pady=3)
            self.thickness_vars[key] = tk.DoubleVar(value=default)
            ttk.Entry(thickness_frame, textvariable=self.thickness_vars[key], width=10).grid(row=i, column=1, padx=5, pady=3)
        
        ttk.Label(thickness_frame,
                  text="Structure: Glass | ITO | PEDOT:PSS | Water",
                  font=('Arial', 9, 'italic')).grid(row=2, column=0, columnspan=2, pady=5)
        
        ttk.Button(thickness_frame, text="Reset to Defaults",
                   command=self.reset_thicknesses).grid(row=3, column=0, columnspan=2, pady=10)
        
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
        voltage_info = (f"Available: {min(available_voltages):.0f} to {max(available_voltages):.0f} mV"
                        if available_voltages else "No voltage data loaded")
        ttk.Label(voltage_frame, text=voltage_info, font=('Arial', 8, 'italic')).grid(row=1, column=0, columnspan=2, pady=2)
        
        # TMM section
        tmm_frame = ttk.LabelFrame(main_frame, text="TMM (Isotropic) Refractive Indices", padding="10")
        tmm_frame.pack(fill=tk.X, padx=10, pady=5)
        ttk.Label(tmm_frame,
                  text="Uses single n and k for each mode (isotropic)",
                  font=('Arial', 9, 'italic'), foreground='darkblue').grid(row=0, column=0, columnspan=4, pady=5)
        
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
        
        # Berreman section
        berreman_frame = ttk.LabelFrame(main_frame, text="Berreman 4x4 (Anisotropic) Refractive Indices", padding="10")
        berreman_frame.pack(fill=tk.X, padx=10, pady=5)
        ttk.Label(berreman_frame,
                  text="Uses n_xy, k_xy (in-plane) and n_z, k_z (out-of-plane)",
                  font=('Arial', 9, 'italic'), foreground='darkgreen').grid(row=0, column=0, columnspan=4, pady=5)
        
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
        
        # Experimental data visibility
        exp_frame = ttk.LabelFrame(main_frame, text="Experimental Data", padding="10")
        exp_frame.pack(fill=tk.X, padx=10, pady=5)
        self.show_exp_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            exp_frame,
            text="Show experimental data (5500rpm_1.6_s.xlsx / 5500rpm_1.6_p.xlsx)",
            variable=self.show_exp_var
        ).pack()
        
        # Control buttons
        buttons_frame = ttk.Frame(main_frame)
        buttons_frame.pack(pady=10)
        ttk.Button(buttons_frame, text="Plot TMM", command=self.plot_tmm, style='Accent.TButton'
                   ).grid(row=0, column=0, padx=5, pady=5)
        ttk.Button(buttons_frame, text="Plot Berreman", command=self.plot_berreman, style='Accent.TButton'
                   ).grid(row=0, column=1, padx=5, pady=5)
        ttk.Button(buttons_frame, text="📊 Plot All 4 Combined", command=self.plot_all_combined, style='Accent.TButton'
                   ).grid(row=0, column=2, padx=5, pady=5)
        
        self.status_label = ttk.Label(main_frame, text="Ready", foreground="green")
        self.status_label.pack()
        
        info_frame = ttk.LabelFrame(main_frame, text="Information", padding="10")
        info_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        ttk.Label(
            info_frame,
            text=(
                "Two calculation methods are available:\n\n"
                "1. TMM (Transfer Matrix Method)\n"
                "   - Treats PEDOT:PSS as isotropic\n"
                "   - Uses single n and k per mode\n"
                "   - Faster calculation\n\n"
                "2. Berreman 4×4 Matrix Method\n"
                "   - Treats PEDOT:PSS as anisotropic\n"
                "   - Uses n_xy / k_xy (in-plane) and n_z / k_z (out-of-plane)\n"
                "   - More accurate for anisotropic films\n\n"
                "3. Combined plot shows TE/TM from both methods alongside experimental data."
            ),
            justify=tk.LEFT
        ).pack(padx=5, pady=5)
    
    def reset_thicknesses(self):
        for key, value in self.default_thicknesses.items():
            self.thickness_vars[key].set(value)
    
    def get_layer_thicknesses(self):
        return {key: var.get() for key, var in self.thickness_vars.items()}
    
    def plot_tmm(self):
        try:
            self.status_label.config(text="Calculating TMM reflectivity...", foreground="blue")
            self.root.update()
            
            wavelength = self.wavelength_var.get()
            voltage = self.voltage_var.get()
            layer_thicknesses = self.get_layer_thicknesses()
            
            fig = plot_combined_te_tm_tmm(
                wavelength, voltage,
                self.n_tmm_te_var.get(), self.k_tmm_te_var.get(),
                self.n_tmm_tm_var.get(), self.k_tmm_tm_var.get(),
                layer_thicknesses,
                shift_te=self.shift_tmm_te_var.get(),
                shift_tm=self.shift_tmm_tm_var.get()
            )
            
            self.status_label.config(text="Ready", foreground="green")
            plt.show()
        except Exception as e:
            self.status_label.config(text=f"Error: {e}", foreground="red")
            messagebox.showerror("Error", f"TMM plot failed:\n{e}")
    
    def plot_berreman(self):
        try:
            self.status_label.config(text="Calculating Berreman reflectivity...", foreground="blue")
            self.root.update()
            
            wavelength = self.wavelength_var.get()
            voltage = self.voltage_var.get()
            layer_thicknesses = self.get_layer_thicknesses()
            
            fig = plot_combined_te_tm_berreman(
                wavelength, voltage,
                self.n_te_var.get(), self.k_te_var.get(),
                self.n_tm_var.get(), self.k_tm_var.get(),
                layer_thicknesses,
                shift_te=self.shift_te_var.get(),
                shift_tm=self.shift_tm_var.get()
            )
            
            self.status_label.config(text="Ready", foreground="green")
            plt.show()
        except Exception as e:
            self.status_label.config(text=f"Error: {e}", foreground="red")
            messagebox.showerror("Error", f"Berreman plot failed:\n{e}")
    
    def plot_all_combined(self):
        try:
            self.status_label.config(text="Calculating combined plot...", foreground="blue")
            self.root.update()
            
            wavelength = self.wavelength_var.get()
            voltage = self.voltage_var.get()
            layer_thicknesses = self.get_layer_thicknesses()
            
            fig = plot_all_four_combined(
                wavelength, voltage,
                self.n_tmm_te_var.get(), self.k_tmm_te_var.get(),
                self.n_tmm_tm_var.get(), self.k_tmm_tm_var.get(),
                self.n_te_var.get(), self.k_te_var.get(),
                self.n_tm_var.get(), self.k_tm_var.get(),
                layer_thicknesses,
                shift_tmm_te=self.shift_tmm_te_var.get(),
                shift_tmm_tm=self.shift_tmm_tm_var.get(),
                shift_te=self.shift_te_var.get(),
                shift_tm=self.shift_tm_var.get()
            )
            
            self.status_label.config(text="Ready", foreground="green")
            plt.show()
        except Exception as e:
            self.status_label.config(text=f"Error: {e}", foreground="red")
            messagebox.showerror("Error", f"Combined plot failed:\n{e}")

# Run the GUI
if __name__ == "__main__":
    root = tk.Tk()
    app = ReflectivityGUI(root)
    root.mainloop()


