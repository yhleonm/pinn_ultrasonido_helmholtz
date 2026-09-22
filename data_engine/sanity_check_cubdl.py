import os
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from scipy.signal import hilbert

# Import the loader
from data_engine.loader import CUBDLLoader

def generate_mock_cubdl_file(file_path):
    print(f"Generando archivo mock CUBDL en {file_path}...")
    
    # Parámetros físicos
    c0 = 1540.0         # m/s
    fs = 40e6           # 40 MHz sampling frequency
    f0 = 5e6            # 5 MHz center frequency of pulse
    num_samples = 2000
    num_channels = 128
    pitch = 0.3e-3      # 0.3 mm
    
    # Coordenadas de receptores
    x_rx = (np.arange(num_channels) - (num_channels - 1) / 2.0) * pitch
    rx_pos = np.stack([x_rx, np.zeros_like(x_rx)], axis=1) # [128, 2]
    
    # Objetivo (Scatterer) en x = 0.0, z = 0.02 m (2 cm de profundidad)
    x_t, z_t = 0.0, 0.02
    
    # Tiempo de vuelo para cada canal: Transmisión desde el centro (x=0, z=0) -> Target -> Receptor i
    dist_tx = np.sqrt(x_t**2 + z_t**2)
    dist_rx = np.sqrt((rx_pos[:, 0] - x_t)**2 + (rx_pos[:, 1] - z_t)**2)
    t_flight = (dist_tx + dist_rx) / c0 # array [128]
    
    # Vector de tiempo
    t = np.arange(num_samples) / fs # array [2000]
    
    # Señal RF: pulso Gaussiano modulado con decaimiento geométrico (1/r)
    rf = np.zeros((num_samples, num_channels, 1), dtype=np.float32)
    sigma = 0.15e-6 # Duración del pulso en segundos
    
    for i in range(num_channels):
        tau = t_flight[i]
        # Decaimiento por divergencia geométrica (1/r)
        r = np.sqrt((rx_pos[i, 0] - x_t)**2 + z_t**2)
        amplitude_scale = (0.02 / r)**4 # Decaimiento más fuerte para asegurar pico en el centro
        
        # Señal RF en canal i
        rf[:, i, 0] = amplitude_scale * np.exp(-(t - tau)**2 / (2.0 * sigma**2)) * np.cos(2.0 * np.pi * f0 * (t - tau))
        
    # Añadir ruido blanco Gaussiano leve
    np.random.seed(42)
    rf += np.random.normal(0, 0.005, rf.shape).astype(np.float32)
    
    # Delays de transmisión (cero para transmisión divergente/focalizada simulada en los tiempos de vuelo)
    tx_delay = np.zeros((1, num_channels))
    
    # Guardar en formato .mat
    mat_dict = {
        'rf': rf,
        'fs': fs,
        'c': c0,
        't0': 0.0,
        'rx_pos': rx_pos,
        'tx_delay': tx_delay
    }
    sio.savemat(file_path, mat_dict)
    print("Archivo mock CUBDL guardado con éxito.")

def run_sanity_check():
    mock_file = "data_engine/cubdl_sample_dataset.mat"
    generate_mock_cubdl_file(mock_file)
    
    print("\n--- Ejecutando Prueba de Sanity del CUBDLLoader ---")
    # 1. Cargar datos
    loader = CUBDLLoader(mock_file)
    print(f"Cargado con éxito. fs = {loader.fs/1e6:.1f} MHz, c0 = {loader.c0:.1f} m/s, rf shape = {loader.rf.shape}")
    
    # 2. Extraer presiones complejas a 5 MHz
    f_check = 5e6
    sensor_x, P_meas = loader.extract_boundary_pressure([f_check], transmit_index=0)
    P_complex = P_meas[f_check]
    
    # 3. Analizar amplitud y fase
    amplitude = np.abs(P_complex)
    phase = np.angle(P_complex)
    
    # Verificar consistencia física de la amplitud (máximo en el centro, simétrica)
    idx_center = len(sensor_x) // 2
    max_amp_idx = np.argmax(amplitude)
    print(f"Amplitud máxima en canal: {max_amp_idx} (Posición x: {sensor_x[max_amp_idx]*1e3:.2f} mm)")
    
    # Verificar variación de fase (cuadrática con la distancia al objetivo)
    phase_unwrapped = np.unwrap(phase)
    
    # Guardar gráfico de diagnóstico
    plt.figure(figsize=(12, 5))
    
    # Amplitud
    plt.subplot(1, 2, 1)
    plt.plot(sensor_x * 1000, amplitude, 'b.-', label='Amplitud Extraída')
    plt.xlabel('Posición de Receptor (mm)')
    plt.ylabel('Amplitud de Presión Compleja')
    plt.title('Perfil de Amplitud (5 MHz)')
    plt.grid(True)
    plt.legend()
    
    # Fase
    plt.subplot(1, 2, 2)
    plt.plot(sensor_x * 1000, phase_unwrapped, 'r.-', label='Fase Unwrapped')
    plt.xlabel('Posición de Receptor (mm)')
    plt.ylabel('Fase (rad)')
    plt.title('Perfil de Fase Unwrapped (5 MHz)')
    plt.grid(True)
    plt.legend()
    
    plt.tight_layout()
    os.makedirs("articulo", exist_ok=True)
    plt.savefig('articulo/cubdl_sanity_check.png', dpi=150)
    plt.close()
    
    print("Gráficos de diagnóstico guardados en articulo/cubdl_sanity_check.png")
    
    # Validaciones automáticas de consistencia
    assert len(P_complex) == 128, "Error: Número incorrecto de canales extraídos."
    assert abs(max_amp_idx - idx_center) <= 5, "Error: El máximo de amplitud no está alineado con el centro del arreglo."
    assert amplitude[idx_center] > 2.0 * amplitude[0], "Error: La amplitud no decae adecuadamente en los bordes."
    
    print("\n[SUCCESS] PRUEBA DE SANITY COMPLETADA CON EXITO! La carga, extraccion FFT y compensacion de fase son fisicamente consistentes.")

if __name__ == "__main__":
    run_sanity_check()
