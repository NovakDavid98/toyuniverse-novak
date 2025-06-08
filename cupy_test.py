import cupy as cp
import numpy as np

print(f"CuPy version: {cp.__version__}")
print(f"NumPy version: {np.__version__}")

try:
    print(f"Number of available GPUs: {cp.cuda.runtime.getDeviceCount()}")
    if cp.cuda.runtime.getDeviceCount() > 0:
        device_id = cp.cuda.runtime.getDevice()
        print(f"Current GPU Device ID: {device_id}")
        props = cp.cuda.runtime.getDeviceProperties(device_id)
        print(f"GPU Name: {props['name'].decode()}") # props['name'] is bytes
        print(f"Compute Capability: {props['major']}.{props['minor']}")
        print(f"Total Global Memory: {props['totalGlobalMem'] / (1024**2):.2f} MB")

        print("\n--- Simple CuPy Test ---")
        # Create a CuPy array
        x_gpu = cp.arange(12, dtype=cp.float32).reshape(3, 4)
        print("x_gpu (CuPy array):\n", x_gpu)
        print("type(x_gpu):", type(x_gpu))

        # Perform an operation on the GPU
        y_gpu = cp.sin(x_gpu) * 2
        print("\ny_gpu = cp.sin(x_gpu) * 2 (CuPy array):\n", y_gpu)

        # Transfer to CPU
        y_cpu = cp.asnumpy(y_gpu)
        print("\ny_cpu (NumPy array after cp.asnumpy(y_gpu)):\n", y_cpu)
        print("type(y_cpu):", type(y_cpu))

        # Create a NumPy array and transfer to GPU
        a_cpu = np.array([[1,2],[3,4]], dtype=np.float32)
        print("\na_cpu (NumPy array):\n", a_cpu)
        a_gpu = cp.asarray(a_cpu)
        print("a_gpu (CuPy array after cp.asarray(a_cpu)):\n", a_gpu)
        print("type(a_gpu):", type(a_gpu))
        
        print("\nCuPy test completed successfully.")
    else:
        print("No GPU devices found by CuPy.")

except cp.cuda.runtime.CUDARuntimeError as e:
    print(f"CuPy CUDA Runtime Error: {e}")
except Exception as e:
    print(f"An unexpected error occurred: {e}")
