import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import correlate2d
from scipy.ndimage import map_coordinates
import os
import argparse

def radial_profile(data, center=None):
    """Computes the radial average of 2D data."""
    y, x = np.indices(data.shape)
    if center is None:
        center = np.array([(x.max() - x.min()) / 2.0, (y.max() - y.min()) / 2.0])

    r = np.hypot(x - center[0], y - center[1])

    # Integer part of radii for binning
    r_int = r.astype(int)

    # Sum of data for each radius
    tbin = np.bincount(r_int.ravel(), data.ravel())
    # Number of pixels for each radius
    nr = np.bincount(r_int.ravel())

    # Avoid division by zero for radii with no pixels
    radialprofile = np.zeros_like(tbin, dtype=float)
    non_zero_nr = nr > 0
    radialprofile[non_zero_nr] = tbin[non_zero_nr] / nr[non_zero_nr]
    return radialprofile

def analyze_snapshot(npy_file_path, output_dir):
    """Performs FFT and ACF analysis on a single snapshot."""
    if not os.path.exists(npy_file_path):
        print(f"Error: Snapshot file not found at {npy_file_path}")
        return

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    base_filename = os.path.splitext(os.path.basename(npy_file_path))[0]

    data = np.load(npy_file_path)
    # Ensure data is float for FFT/ACF
    data = data.astype(np.float64)
    # Detrend data (subtract mean) - important for FFT/ACF
    data -= np.mean(data)

    # --- 2D FFT Analysis ---
    fft_data = np.fft.fft2(data)
    fft_shifted = np.fft.fftshift(fft_data)
    power_spectrum_2d = np.abs(fft_shifted)**2

    plt.figure(figsize=(7, 6))
    plt.imshow(np.log10(power_spectrum_2d + 1e-9), cmap='viridis') # Add small epsilon to avoid log(0)
    plt.colorbar(label='Log10(Power)')
    plt.title(f'2D Power Spectrum - {base_filename}')
    plt.xticks([])
    plt.yticks([])
    fft_plot_path = os.path.join(output_dir, f'{base_filename}_fft_2d.png')
    plt.savefig(fft_plot_path)
    plt.close()
    print(f"Saved 2D FFT plot to: {fft_plot_path}")

    # Radially Averaged FFT Power Spectrum
    radial_fft = radial_profile(power_spectrum_2d)
    # Determine appropriate k-axis (spatial frequency)
    # k = 2*pi / L, where L is wavelength. Max k is related to Nyquist frequency.
    # For plotting, use pixel units for k for now.
    k_values = np.arange(len(radial_fft))

    plt.figure(figsize=(8, 5))
    plt.plot(k_values[1:len(radial_fft)//2], radial_fft[1:len(radial_fft)//2]) # Avoid DC (k=0) and high frequencies if noisy
    plt.xlabel('Spatial Frequency k (pixels^-1)')
    plt.ylabel('Radially Averaged Power')
    plt.title(f'Radially Averaged Power Spectrum - {base_filename}')
    plt.grid(True)
    radial_fft_plot_path = os.path.join(output_dir, f'{base_filename}_fft_radial.png')
    plt.savefig(radial_fft_plot_path)
    plt.close()
    print(f"Saved Radial FFT plot to: {radial_fft_plot_path}")

    # --- 2D ACF Analysis ---
    # Normalize data for ACF (often helps interpretation, std=1)
    # data_norm = (data - np.mean(data)) / np.std(data)
    # Using original detrended data as per common practice for physical patterns
    acf_data = correlate2d(data, data, mode='same', boundary='wrap') # 'wrap' for periodic-like boundary
    # Shift ACF so that zero lag is at the center
    # For correlate2d 'same', max is already near center, but let's ensure it for plotting
    # This manual shift might not be perfect if correlate2d doesn't center perfectly for all shapes.
    # A more robust way is to use fftshift on an ACF computed via FFTs, but correlate2d is simpler.
    # For now, assume correlate2d with 'same' centers it well enough for visual inspection.

    plt.figure(figsize=(7, 6))
    plt.imshow(acf_data, cmap='viridis')
    plt.colorbar(label='Autocorrelation')
    plt.title(f'2D Autocorrelation Function - {base_filename}')
    plt.xticks([])
    plt.yticks([])
    acf_plot_path = os.path.join(output_dir, f'{base_filename}_acf_2d.png')
    plt.savefig(acf_plot_path)
    plt.close()
    print(f"Saved 2D ACF plot to: {acf_plot_path}")

    # Radially Averaged ACF
    radial_acf = radial_profile(acf_data)
    lag_values = np.arange(len(radial_acf))

    plt.figure(figsize=(8, 5))
    plt.plot(lag_values[0:len(radial_acf)//2], radial_acf[0:len(radial_acf)//2]) # Plot up to half the domain size
    plt.xlabel('Lag Distance (pixels)')
    plt.ylabel('Radially Averaged Autocorrelation')
    plt.title(f'Radially Averaged ACF - {base_filename}')
    plt.grid(True)
    radial_acf_plot_path = os.path.join(output_dir, f'{base_filename}_acf_radial.png')
    plt.savefig(radial_acf_plot_path)
    plt.close()
    print(f"Saved Radial ACF plot to: {radial_acf_plot_path}")

    print(f"Analysis complete for {base_filename}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform FFT and ACF analysis on a single .npy snapshot.')
    parser.add_argument('snapshot_file', type=str, help='Path to the .npy snapshot file.')
    parser.add_argument('output_directory', type=str, help='Directory to save the output plots.')
    
    args = parser.parse_args()
    
    analyze_snapshot(args.snapshot_file, args.output_directory)
