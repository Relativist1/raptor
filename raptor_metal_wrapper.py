#!/usr/bin/env python3
"""
Python wrapper for RAPTOR Metal GPU Implementation

This module provides a Python interface to the Metal-based RAPTOR ray-tracing code,
making it easy to integrate with existing Python astrophysics workflows.

Requirements:
- macOS 10.15+ (Catalina)
- Metal-capable GPU (M2 chip supported)
- Compiled raptor_metal binary
"""

import numpy as np
import subprocess
import json
import os
import tempfile
from typing import Dict, List, Tuple, Optional, Union
import matplotlib.pyplot as plt
from pathlib import Path

class RAPTORMetalWrapper:
    """
    Python wrapper for the Metal-based RAPTOR implementation.
    
    This class provides a high-level interface to run ray-tracing simulations
    on the GPU using Apple's Metal framework.
    """
    
    def __init__(self, binary_path: str = "./raptor_metal"):
        """
        Initialize the RAPTOR Metal wrapper.
        
        Args:
            binary_path: Path to the compiled raptor_metal binary
        """
        self.binary_path = binary_path
        self.model_params = self._get_default_params()
        
        # Check if binary exists and is executable
        if not os.path.exists(binary_path):
            raise FileNotFoundError(f"RAPTOR Metal binary not found at {binary_path}")
        if not os.access(binary_path, os.X_OK):
            raise PermissionError(f"RAPTOR Metal binary not executable at {binary_path}")
    
    def _get_default_params(self) -> Dict:
        """Get default model parameters matching the C implementation."""
        return {
            'MBH': 6.2e9,           # Black hole mass in solar masses
            'M_UNIT': 1e9,          # Mass unit in solar masses
            'TIME_INIT': 0.0,
            'INCLINATION': 17.0,    # degrees
            'R_HIGH': 1.0,
            'R_LOW': 1.0,
            'FREQS_PER_DEC': 4,
            'FREQ_MIN': 1e10,       # 10 GHz
            'FREQ_MAX': 1e15,       # 1000 GHz
            'SOURCE_DIST': 16.4,    # kpc
            'IMG_WIDTH': 256,
            'IMG_HEIGHT': 256,
            'CAM_SIZE_X': 40.0,     # GM/c^2
            'CAM_SIZE_Y': 40.0,     # GM/c^2
            'STEPSIZE': 0.1,
            'max_level': 4
        }
    
    def set_model_parameters(self, **kwargs):
        """
        Set model parameters.
        
        Args:
            **kwargs: Model parameters to set
        """
        for key, value in kwargs.items():
            if key in self.model_params:
                self.model_params[key] = value
            else:
                print(f"Warning: Unknown parameter '{key}' ignored")
    
    def run_simulation(self, 
                      output_dir: str = "./output",
                      save_lightpaths: bool = False,
                      verbose: bool = True) -> Dict:
        """
        Run the RAPTOR Metal simulation.
        
        Args:
            output_dir: Directory to save output files
            save_lightpaths: Whether to save light path data
            verbose: Whether to print progress information
            
        Returns:
            Dictionary containing simulation results
        """
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Create temporary input file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            input_file = f.name
            json.dump(self.model_params, f, indent=2)
        
        try:
            # Run the simulation
            cmd = [self.binary_path, input_file, output_dir]
            if save_lightpaths:
                cmd.append("--save-lightpaths")
            
            if verbose:
                print(f"Running RAPTOR Metal simulation...")
                print(f"Command: {' '.join(cmd)}")
                print(f"Model parameters: {self.model_params}")
            
            # Execute the command
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
            
            if result.returncode != 0:
                raise RuntimeError(f"Simulation failed: {result.stderr}")
            
            if verbose:
                print("Simulation completed successfully!")
                print(f"Output: {result.stdout}")
            
            # Parse results
            results = self._parse_output(output_dir)
            
            return results
            
        finally:
            # Clean up temporary file
            os.unlink(input_file)
    
    def _parse_output(self, output_dir: str) -> Dict:
        """
        Parse simulation output files.
        
        Args:
            output_dir: Directory containing output files
            
        Returns:
            Dictionary containing parsed results
        """
        results = {
            'model_params': self.model_params.copy(),
            'images': {},
            'spectra': {},
            'metadata': {}
        }
        
        # Look for output files
        output_path = Path(output_dir)
        
        # Parse image files
        for img_file in output_path.glob("image_*.h5"):
            # This would require h5py to read HDF5 files
            # For now, just store the filename
            freq = img_file.stem.split('_')[1]
            results['images'][freq] = str(img_file)
        
        # Parse spectrum files
        for spec_file in output_path.glob("spectrum_*.dat"):
            try:
                spectrum_data = np.loadtxt(spec_file)
                if len(spectrum_data.shape) == 2:
                    results['spectra'][spec_file.stem] = {
                        'frequencies': spectrum_data[:, 0],
                        'I': spectrum_data[:, 1],
                        'Q': spectrum_data[:, 2] if spectrum_data.shape[1] > 2 else None,
                        'U': spectrum_data[:, 3] if spectrum_data.shape[1] > 3 else None,
                        'V': spectrum_data[:, 4] if spectrum_data.shape[1] > 4 else None
                    }
            except Exception as e:
                print(f"Warning: Could not parse spectrum file {spec_file}: {e}")
        
        return results
    
    def plot_image(self, results: Dict, frequency: Optional[float] = None, 
                   stokes: str = 'I', save_path: Optional[str] = None):
        """
        Plot simulation results.
        
        Args:
            results: Simulation results from run_simulation()
            frequency: Frequency to plot (if None, use first available)
            stokes: Stokes parameter to plot ('I', 'Q', 'U', 'V')
            save_path: Path to save the plot (if None, display)
        """
        if not results['images']:
            print("No image data available for plotting")
            return
        
        # For now, create a simple visualization
        # In a full implementation, this would read the actual image data
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        
        # Create a mock image based on model parameters
        width = self.model_params['IMG_WIDTH']
        height = self.model_params['IMG_HEIGHT']
        
        # Simple mock image (replace with actual data reading)
        x = np.linspace(-self.model_params['CAM_SIZE_X']/2, 
                        self.model_params['CAM_SIZE_X']/2, width)
        y = np.linspace(-self.model_params['CAM_SIZE_Y']/2, 
                        self.model_params['CAM_SIZE_Y']/2, height)
        X, Y = np.meshgrid(x, y)
        
        # Mock intensity (replace with actual data)
        R = np.sqrt(X**2 + Y**2)
        intensity = np.exp(-R/10) * (1 + 0.1 * np.sin(2 * np.pi * X / 5))
        
        im = ax.imshow(intensity, extent=[x.min(), x.max(), y.min(), y.max()], 
                       origin='lower', cmap='hot')
        ax.set_xlabel('α (GM/c²)')
        ax.set_ylabel('β (GM/c²)')
        ax.set_title(f'RAPTOR Metal Simulation - Stokes {stokes}')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Intensity (arbitrary units)')
        
        # Add black hole location
        ax.plot(0, 0, 'ko', markersize=10, label='Black Hole')
        ax.legend()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {save_path}")
        else:
            plt.show()
    
    def plot_spectrum(self, results: Dict, save_path: Optional[str] = None):
        """
        Plot spectrum results.
        
        Args:
            results: Simulation results from run_simulation()
            save_path: Path to save the plot (if None, display)
        """
        if not results['spectra']:
            print("No spectrum data available for plotting")
            return
        
        # Get first available spectrum
        spectrum_key = list(results['spectra'].keys())[0]
        spectrum = results['spectra'][spectrum_key]
        
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        
        # Plot Stokes I
        ax.loglog(spectrum['frequencies'], spectrum['I'], 'b-', linewidth=2, label='Stokes I')
        
        # Plot other Stokes parameters if available
        if spectrum['Q'] is not None:
            ax.loglog(spectrum['frequencies'], np.abs(spectrum['Q']), 'r--', 
                     linewidth=1, label='|Stokes Q|')
        if spectrum['U'] is not None:
            ax.loglog(spectrum['frequencies'], np.abs(spectrum['U']), 'g--', 
                     linewidth=1, label='|Stokes U|')
        if spectrum['V'] is not None:
            ax.loglog(spectrum['frequencies'], np.abs(spectrum['V']), 'm--', 
                     linewidth=1, label='|Stokes V|')
        
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Flux Density (Jy)')
        ax.set_title('RAPTOR Metal Spectrum')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Spectrum plot saved to {save_path}")
        else:
            plt.show()
    
    def benchmark(self, num_runs: int = 5) -> Dict:
        """
        Benchmark the GPU performance.
        
        Args:
            num_runs: Number of runs to average
            
        Returns:
            Dictionary containing timing information
        """
        print(f"Running benchmark with {num_runs} iterations...")
        
        timings = []
        for i in range(num_runs):
            print(f"Run {i+1}/{num_runs}...")
            
            start_time = time.time()
            try:
                results = self.run_simulation(verbose=False)
                end_time = time.time()
                timings.append(end_time - start_time)
            except Exception as e:
                print(f"Run {i+1} failed: {e}")
        
        if not timings:
            raise RuntimeError("All benchmark runs failed")
        
        benchmark_results = {
            'num_runs': num_runs,
            'timings': timings,
            'mean_time': np.mean(timings),
            'std_time': np.std(timings),
            'min_time': np.min(timings),
            'max_time': np.max(timings),
            'pixels_per_second': (self.model_params['IMG_WIDTH'] * 
                                 self.model_params['IMG_HEIGHT']) / np.mean(timings)
        }
        
        print(f"Benchmark results:")
        print(f"  Mean time: {benchmark_results['mean_time']:.3f} ± {benchmark_results['std_time']:.3f} s")
        print(f"  Min time: {benchmark_results['min_time']:.3f} s")
        print(f"  Max time: {benchmark_results['max_time']:.3f} s")
        print(f"  Performance: {benchmark_results['pixels_per_second']:.0f} pixels/s")
        
        return benchmark_results


def main():
    """Example usage of the RAPTOR Metal wrapper."""
    try:
        # Initialize wrapper
        raptor = RAPTORMetalWrapper()
        
        # Customize parameters
        raptor.set_model_parameters(
            IMG_WIDTH=512,
            IMG_HEIGHT=512,
            FREQ_MIN=1e11,  # 100 GHz
            FREQ_MAX=1e14,  # 100 THz
            INCLINATION=30.0
        )
        
        # Run simulation
        print("Starting RAPTOR Metal simulation...")
        results = raptor.run_simulation(output_dir="./raptor_output")
        
        # Plot results
        print("Generating plots...")
        raptor.plot_image(results, save_path="./raptor_image.png")
        raptor.plot_spectrum(results, save_path="./raptor_spectrum.png")
        
        # Run benchmark
        print("Running performance benchmark...")
        benchmark_results = raptor.benchmark(num_runs=3)
        
        print("Simulation completed successfully!")
        
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    import time
    exit(main())