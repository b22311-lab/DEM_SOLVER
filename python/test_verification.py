#!/usr/bin/env python3
"""
Verification Tests for DEM Solver
Compares numerical results against analytical solutions for:
- Test 1: Free Fall
- Test 2: Constant Velocity  
- Test 3: Particle Bounce
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import glob
import re
import argparse
from pathlib import Path

# Configuration for publication-quality plots
rcParams['font.size'] = 11
rcParams['font.family'] = 'serif'
rcParams['axes.labelsize'] = 12
rcParams['axes.titlesize'] = 13
rcParams['xtick.labelsize'] = 10
rcParams['ytick.labelsize'] = 10
rcParams['legend.fontsize'] = 10
rcParams['figure.figsize'] = (10, 6)
rcParams['figure.dpi'] = 150

# For the tuned bounce case, 8 s captures about five rebounds.
TEST3_ANALYSIS_END_TIME = 8.0
TEST3_TARGET_BOUNCES = 5

# Try to import vtk reader
try:
    import pyvista as pv
    USE_PYVISTA = True
except ImportError:
    USE_PYVISTA = False
    print("WARNING: PyVista not found. Installing: pip install pyvista")
    exit(1)


class VTKReader:
    """Read particle data from VTK files"""
    
    def __init__(self, vtk_file):
        self.file = vtk_file
        self.data = {}
        self.load_data()
    
    def load_data(self):
        """Load VTK file and extract data"""
        try:
            mesh = pv.read(self.file)
            
            # Extract points (positions in 3D)
            points_3d = mesh.points
            
            # Convert 3D to 2D (drop Z)
            if points_3d.shape[1] >= 2:
                self.data['position'] = points_3d[:, :2].copy()
            
            # Extract velocity field
            if 'velocity' in mesh.array_names:
                vel_3d = mesh['velocity']
                if vel_3d.shape[1] >= 2:
                    self.data['velocity'] = vel_3d[:, :2].copy()
            
            # Extract scalars
            if 'height' in mesh.array_names:
                self.data['height'] = mesh['height'].copy()
            
            if 'kinetic_energy' in mesh.array_names:
                self.data['ke'] = mesh['kinetic_energy'].copy()
            
            if 'radius' in mesh.array_names:
                self.data['radius'] = mesh['radius'].copy()
        
        except Exception as e:
            print(f"ERROR reading {self.file}: {e}")


def get_vtk_files(pattern='output_*.vtk'):
    """Get sorted list of VTK output files"""
    files = sorted(glob.glob(pattern), 
                   key=lambda x: int(re.search(r'_(\d+)\.vtk', x).group(1)))
    return files


def extract_time_from_vtk_header(vtk_file):
    """Read physical time from VTK header line written by the solver.

    Expected second line format:
        Particle data at time <value>
    """
    try:
        with open(vtk_file, 'r') as f:
            _ = f.readline()
            header_line = f.readline().strip()

        match = re.search(r'time\s+([\d.eE+-]+)', header_line)
        if match:
            return float(match.group(1))
    except Exception:
        pass

    return None

def extract_time_series(vtk_files, dt=0.001, write_interval=10):
    """Extract particle trajectories from all VTK files
    
    Args:
        dt: simulation timestep (0.001 s from config)
        write_interval: how many timesteps between writes (10 from config)
    """
    data_list = []
    
    iter_nums = []
    vtk_header_times = []

    for vtk_file in vtk_files:
        reader = VTKReader(vtk_file)
        if reader.data and 'position' in reader.data:
            data_list.append(reader.data)
            match = re.search(r'output_(\d+)\.vtk', vtk_file)
            if match:
                iter_nums.append(int(match.group(1)))
            else:
                iter_nums.append(len(iter_nums) + 1)
            vtk_header_times.append(extract_time_from_vtk_header(vtk_file))
    
    if not data_list:
        print("ERROR: No valid data extracted from VTK files")
        return None
    
    # Prefer physical time from VTK header when present; fallback to iter * dt.
    if vtk_header_times and all(t is not None for t in vtk_header_times):
        times = np.array(vtk_header_times)
        time_source = "vtk_header"
    else:
        times = np.array(iter_nums) * dt
        time_source = "filename_iter_x_dt"
    
    # Extract field data
    n_steps = len(data_list)
    n_particles = data_list[0]['position'].shape[0]
    
    result = {
        'time': times,
        'iter': np.array(iter_nums),
        'position': np.array([d.get('position', np.zeros((n_particles, 2))) for d in data_list]),
        'velocity': np.array([d.get('velocity', np.zeros((n_particles, 2))) for d in data_list]),
        'height': np.array([d.get('height', np.zeros(n_particles)) for d in data_list]),
        'ke': np.array([d.get('ke', np.zeros(n_particles)) for d in data_list]),
        'radius': np.array([d.get('radius', np.zeros(n_particles)) for d in data_list]),
    }
    
    print(f"\nExtracted data:")
    print(f"  Timesteps: {n_steps}")
    print(f"  Particles: {n_particles}")
    print(f"  Time range: {times[0]:.4f} - {times[-1]:.4f} s")
    print(f"  dt: {dt:.4f} s")
    print(f"  Time source: {time_source}")
    
    return result

def test1_freefall(data):
    """Test 1: Free Fall Verification"""
    print("\n" + "="*60)
    print("TEST 1: FREE FALL VERIFICATION")
    print("="*60)
    
    # For single particle
    t = data['time']
    pos = data['position'][:, 0, :]  # Select particle 0
    y_num = pos[:, 1]  # Y position
    
    vel = data['velocity'][:, 0, :]  # particle 0 velocity
    vy_num = vel[:, 1]
    
    ke_num = data['ke'][:, 0]  # Numerical KE from VTK
    
    # Initial conditions (from config)
    y0 = 4.5
    g = 9.81
    
    # Particle properties from config
    density = 1e-09
    rad_init = 0.05
    mass = density * (4.0/3.0 * np.pi * rad_init**3)
    
    # Analytical solutions
    y_analytical = y0 - 0.5 * g * t**2
    vy_analytical = -g * t
    ke_analytical = 0.5 * mass * vy_analytical**2
    
    # Compute errors (skip first few steps for stability)
    valid_idx = (t > 0.01) & (y_num > 0.1)  # Avoid ground contact
    if np.sum(valid_idx) > 5:
        y_err = np.abs(y_num[valid_idx] - y_analytical[valid_idx])
        rmse_y = np.sqrt(np.mean(y_err**2))
        max_error_y = np.max(y_err)
        
        vy_err = np.abs(vy_num[valid_idx] - vy_analytical[valid_idx])
        rmse_vy = np.sqrt(np.mean(vy_err**2))
        max_error_vy = np.max(vy_err)
        
        ke_err = np.abs(ke_num[valid_idx] - ke_analytical[valid_idx])
        rmse_ke = np.sqrt(np.mean(ke_err**2))

        
        
        print(f"Initial height y₀:     {y0:.6f} m")
        print(f"Gravity g:             {g:.2f} m/s²")
        print(f"Mass:                  {mass:.3e} kg")
        print(f"Simulation time:       {t[-1]:.3f} s")
        print(f"Timestep dt:           0.001 s")
        print(f"Number of steps:       {len(t)}")
        print(f"\nRMS Errors:")
        print(f"  Position:  {rmse_y:.2e} m ({100*rmse_y/y0:.3f}% of initial height)")
        print(f"  Velocity:  {rmse_vy:.2e} m/s ({100*rmse_vy/(-g*t[-1]):.3f}% of final velocity)")
        print(f"  KE:        {rmse_ke:.2e} J ({100*rmse_ke/ke_analytical[valid_idx][-1]:.3f}% of final KE)")
        print(f"\nMax Errors:")
        print(f"  Position:  {max_error_y:.2e} m")
        print(f"  Velocity:  {max_error_vy:.2e} m/s")
    
    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    
    # Position vs time
    ax = axes[0, 0]
    ax.plot(t, y_num, 'b.-', label='Numerical', linewidth=2, markersize=4)
    ax.plot(t, y_analytical, 'r--', label='Analytical', linewidth=2)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Height y [m]')
    ax.set_title('Position vs Time - Free Fall')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Velocity vs time
    ax = axes[0, 1]
    ax.plot(t, vy_num, 'b.-', label='Numerical', linewidth=2, markersize=4)
    ax.plot(t, vy_analytical, 'r--', label='Analytical', linewidth=2)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Velocity vᵧ [m/s]')
    ax.set_title('Velocity vs Time - Free Fall')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Position error
    ax = axes[1, 0]
    error_y = np.abs(y_num - y_analytical)
    ax.semilogy(t, error_y + 1e-12, 'g.-', linewidth=2, markersize=4)  # Add 1e-12 to avoid log(0)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Absolute Error [m]')
    ax.set_title('Position Error vs Time')
    ax.grid(True, alpha=0.3, which='both')
    
    # Kinetic energy (FIXED: use correct mass)
    ax = axes[1, 1]
    ax.plot(t, ke_num, 'b.-', label='Numerical', linewidth=2, markersize=4)
    ax.plot(t, ke_analytical, 'r--', label='Analytical', linewidth=2)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Kinetic Energy KE [J]')
    ax.set_title('Kinetic Energy vs Time')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('test1_freefall_verification.png', dpi=150, bbox_inches='tight')
    print("\nPlot saved: test1_freefall_verification.png")
    plt.close()


def test2_constant_velocity(data):
    """Test 2: Constant Velocity Verification"""
    print("\n" + "="*60)
    print("TEST 2: CONSTANT VELOCITY VERIFICATION")
    print("="*60)
    
    t = data['time']
    pos = data['position'][:, 0, :]  # particle 0
    x_num = pos[:, 0]
    y_num = pos[:, 1]
    
    vel = data['velocity'][:, 0, :]  # particle 0 velocity
    vx_num = vel[:, 0]
    vy_num = vel[:, 1]
    
    # Initial conditions
    x0 = 1.0
    y0 = 5.0
    vx = 1.0
    vy = 0.0
    
    # Analytical (constant velocity)
    x_analytical = x0 + vx * t
    y_analytical = y0 + vy * t
    vx_analytical = np.ones_like(t) * vx
    vy_analytical = np.zeros_like(t) * vy
    
    print(f"Initial position:      ({x0:.3f}, {y0:.3f}) m")
    print(f"Initial velocity:      ({vx:.3f}, {vy:.3f}) m/s")
    print(f"Gravity:               DISABLED (0.0 m/s²)")
    print(f"Simulation time:       {t[-1]:.3f} s")
    
    # Error analysis
    x_error = np.abs(x_num - x_analytical)
    y_error = np.abs(y_num - y_analytical)
    vx_error = np.abs(vx_num - vx_analytical)
    vy_error = np.abs(vy_num - vy_analytical)
    
    rmse_x = np.sqrt(np.mean(x_error**2))
    rmse_y = np.sqrt(np.mean(y_error**2))
    rmse_vx = np.sqrt(np.mean(vx_error**2))
    rmse_vy = np.sqrt(np.mean(vy_error**2))
    
    print(f"\nRMS Errors:")
    print(f"  Position x: {rmse_x:.2e} m")
    print(f"  Position y: {rmse_y:.2e} m")
    print(f"  Velocity vx: {rmse_vx:.2e} m/s")
    print(f"  Velocity vy: {rmse_vy:.2e} m/s")
    print(f"\nDisplacements:")
    print(f"  Actual x displacement:   {x_num[-1] - x_num[0]:.6f} m")
    print(f"  Expected x displacement: {vx * t[-1]:.6f} m")
    print(f"  Actual y displacement:   {y_num[-1] - y_num[0]:.6f} m")
    print(f"  Expected y displacement: {vy * t[-1]:.6f} m")
    
    # Check if gravity accidentally enabled
    if rmse_y > 0.01:
        print("\nWARNING: Large Y error suggests gravity might be enabled!")
        avg_ay = np.mean(np.diff(vy_num) / np.diff(t))
        print(f"  Average Y acceleration: {avg_ay:.3f} m/s²")
    
    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    
    # X position with units
    ax = axes[0, 0]
    ax.plot(t, x_num, 'b.-', label='Numerical', linewidth=2, markersize=4)
    ax.plot(t, x_analytical, 'r--', label='Analytical', linewidth=2)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Position x [m]')
    ax.set_title('X Position vs Time - Constant Velocity')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Y position with units
    ax = axes[0, 1]
    ax.plot(t, y_num, 'b.-', label='Numerical', linewidth=2, markersize=4)
    ax.axhline(y0, color='r', linestyle='--', label='Analytical (y=const)', linewidth=2)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Position y [m]')
    ax.set_title('Y Position vs Time - Should be Constant')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Velocity with units
    ax = axes[1, 0]
    ax.plot(t, vx_num, 'b.-', label='vₓ (numerical)', linewidth=2, markersize=3)
    ax.axhline(vx, color='b', linestyle='--', label='vₓ (analytical)', linewidth=2)
    ax.plot(t, vy_num, 'r.-', label='vᵧ (numerical)', linewidth=2, markersize=3)
    ax.axhline(vy, color='r', linestyle='--', label='vᵧ (analytical)', linewidth=2)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Velocity [m/s]')
    ax.set_title('Velocity Components vs Time - Should be Constant')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # Errors with units
    ax = axes[1, 1]
    ax.semilogy(t, x_error + 1e-12, 'b.-', label='Error in x', linewidth=2, markersize=4)
    ax.semilogy(t, y_error + 1e-12, 'r.-', label='Error in y', linewidth=2, markersize=4)
    ax.semilogy(t, vx_error + 1e-12, 'g.-', label='Error in vₓ', linewidth=2, markersize=4)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Absolute Error [m or m/s]')
    ax.set_title('Error vs Time')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    plt.savefig('test2_constant_velocity_verification.png', dpi=150, bbox_inches='tight')
    print("\nPlot saved: test2_constant_velocity_verification.png")
    plt.close()

def test3_particle_bounce(data):
    """Test 3: Particle Bounce Verification"""
    print("\n" + "="*60)
    print("TEST 3: PARTICLE BOUNCE VERIFICATION")
    print("="*60)
    
    t_all = data['time']
    mask = t_all <= (TEST3_ANALYSIS_END_TIME + 1.0e-12)
    t = t_all[mask]
    pos = data['position'][mask, 0, :]  # particle 0
    y_num = pos[:, 1]
    
    vel = data['velocity'][mask, 0, :]
    vy_num = vel[:, 1]
    
    ke_num = data['ke'][mask, 0]
    
    # Initial conditions
    y0 = 4.5
    g = 9.81
    
    print(f"Initial height y0:     {y0:.6f} m")
    print(f"Gravity g:             {g:.2f} m/s²")
    print(f"Simulation time:       {t[-1]:.3f} s")
    print(f"Analysis end time:     {TEST3_ANALYSIS_END_TIME:.3f} s")
    print(f"Final height:          {y_num[-1]:.6f} m")
    
    # Detect impacts near the ground using velocity sign changes and local minima.
    impacts_candidates = []
    if len(t) >= 3:
        dt_local = np.median(np.diff(t))
    else:
        dt_local = 0.01

    radius_guess = 0.05
    if 'radius' in data:
        radius_local = data['radius'][mask, 0]
        if np.size(radius_local) > 0:
            radius_guess = float(np.nanmedian(radius_local))

    ground_level = 0.0
    contact_band = max(0.005, 0.25 * radius_guess)

    for i in range(1, len(vy_num)-1):
        near_ground = y_num[i] <= (ground_level + radius_guess + contact_band)
        sign_change = (vy_num[i-1] < 0.0 and vy_num[i] >= 0.0) or (vy_num[i-1] < 0.0 and vy_num[i+1] > 0.0)
        local_min = (y_num[i] <= y_num[i-1]) and (y_num[i] <= y_num[i+1])
        if near_ground and (sign_change or local_min):
            impacts_candidates.append(i)

    impacts = []
    min_sep_steps = max(1, int(0.05 / max(dt_local, 1.0e-12)))
    for idx in impacts_candidates:
        if (not impacts) or (idx - impacts[-1] >= min_sep_steps):
            impacts.append(idx)
        else:
            if y_num[idx] < y_num[impacts[-1]]:
                impacts[-1] = idx
    
    print(f"Number of detected impacts: {len(impacts)}")
    print(f"Target impacts:             {TEST3_TARGET_BOUNCES}")

    if len(impacts) < TEST3_TARGET_BOUNCES:
        print("WARNING: Fewer than target impacts detected; increase bounce test end_time.")
    
    if len(impacts) > 1:
        print(f"Impact times: {t[impacts][:5]}")  # First 5 impacts
        print(f"Heights at impact: {y_num[impacts][:5]}")
    
    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    
    # Height vs time
    ax = axes[0, 0]
    ax.plot(t, y_num, 'b-', linewidth=2)
    ax.axhline(0, color='k', linestyle='--', alpha=0.5, label='Ground')
    if impacts:
        ax.plot(t[impacts], y_num[impacts], 'ro', markersize=6, label='Impacts')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Height y (m)')
    ax.set_title('Particle Height vs Time')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Velocity vs time
    ax = axes[0, 1]
    ax.plot(t, vy_num, 'g-', linewidth=2)
    ax.axhline(0, color='k', linestyle='--', alpha=0.5)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Vertical Velocity $v_y$ (m/s)')
    ax.set_title('Velocity vs Time')
    ax.grid(True, alpha=0.3)
    
    # Energy dissipation
    ax = axes[1, 0]
    ke_initial = ke_num[0]
    if ke_initial > 0:
        ax.semilogy(t, ke_num / ke_initial, 'purple', linewidth=2)
        ax.set_ylabel('Normalized KE')
    else:
        ax.semilogy(t, ke_num, 'purple', linewidth=2)
        ax.set_ylabel('Kinetic Energy')
    ax.set_xlabel('Time (s)')
    ax.set_title('Energy Dissipation')
    ax.grid(True, alpha=0.3, which='both')
    
    # Phase space
    ax = axes[1, 1]
    ax.plot(y_num, vy_num, 'b-', linewidth=1, alpha=0.6)
    ax.plot(y_num[0], vy_num[0], 'go', markersize=10, label='Start')
    ax.plot(y_num[-1], vy_num[-1], 'r^', markersize=10, label='End')
    if impacts:
        ax.plot(y_num[impacts], vy_num[impacts], 'ko', markersize=5, alpha=0.5, label='Impacts')
    ax.set_xlabel('Height y (m)')
    ax.set_ylabel('Velocity $v_y$ (m/s)')
    ax.set_title('Phase Space (Height vs Velocity)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('test3_particle_bounce_verification.png', dpi=150, bbox_inches='tight')
    print("\nPlot saved: test3_particle_bounce_verification.png")
    plt.close()


def parse_config_file(config_file='inputs/input_config.in'):
    """Parse Fortran namelist config file and extract time_step
    
    Args:
        config_file: path to the .in config file
        
    Returns:
        time_step (float) or None if not found
    """
    try:
        with open(config_file, 'r') as f:
            content = f.read()
        
        # Look for time_step in simulation_control namelist
        match = re.search(r'time_step\s*=\s*([\d.eE+-]+)', content)
        if match:
            dt = float(match.group(1))
            print(f"  Detected time_step from config: {dt} s")
            return dt
        else:
            print(f"  WARNING: Could not find time_step in {config_file}, using default 0.001")
            return 0.001
    except Exception as e:
        print(f"  ERROR reading config: {e}, using default 0.001")
        return 0.001


def main():
    """Main verification routine"""
    parser = argparse.ArgumentParser(description="DEM verification plot generator")
    parser.add_argument(
        "--test",
        choices=["all", "test1", "test2", "test3"],
        default="all",
        help="Run all verification plots or only one test",
    )
    parser.add_argument(
        "--config",
        default="inputs/input_config.in",
        help="Configuration file used for dt parsing",
    )
    args = parser.parse_args()

    print("\n" + "="*60)
    print("DEM SOLVER VERIFICATION TESTS")
    print("="*60)
    
    # Find VTK files
    vtk_files = get_vtk_files()
    if not vtk_files:
        print("ERROR: No output_*.vtk files found. Run simulation first.")
        return
    
    print(f"\nFound {len(vtk_files)} VTK files")
    
    # Auto-detect timestep from config file
    dt = parse_config_file(args.config)
    
    # Extract data with correct dt
    data = extract_time_series(vtk_files, dt=dt, write_interval=10)
    
    if data is None or data['time'] is None:
        print("ERROR: Could not extract time series")
        return
    
    if args.test in ("all", "test1"):
        try:
            test1_freefall(data)
        except Exception as e:
            print(f"Test 1 failed with error: {e}")
            import traceback
            traceback.print_exc()

    if args.test in ("all", "test2"):
        try:
            test2_constant_velocity(data)
        except Exception as e:
            print(f"Test 2 failed with error: {e}")
            import traceback
            traceback.print_exc()

    if args.test in ("all", "test3"):
        try:
            test3_particle_bounce(data)
        except Exception as e:
            print(f"Test 3 failed with error: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "="*60)
    print("Verification complete!")
    print("="*60)


if __name__ == '__main__':
    main()