from openmm import Platform

print("Available OpenMM platforms:")
for i in range(Platform.getNumPlatforms()):
    p = Platform.getPlatform(i)
    print(f"  {i}: {p.getName()}")
    
# Try CUDA
try:
    cuda = Platform.getPlatformByName("CUDA")
    print(f"\nCUDA available! Speed: {cuda.getSpeed()}")
except Exception as e:
    print(f"\nCUDA NOT available: {e}")
