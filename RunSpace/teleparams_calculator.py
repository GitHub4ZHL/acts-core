#!/usr/bin/env python3
import math
from scipy.stats import norm

# Pixel
pixel_size = 55
pixel_number_x = 256
pixel_number_y = 512
print(f"pixel size: {pixel_size} um")
print(f"pixel number: {pixel_number_x} for x and {pixel_number_y} for y")

# Plane
plane_size_x = pixel_size/1000 * pixel_number_x
plane_size_y = pixel_size/1000 * pixel_number_y
plane_spacing = 30
print(f"plane size: {plane_size_x} mm for x and {plane_size_y} mm for y")
print(f"plane spacing: {plane_spacing} mm")

# Beam
mu = 0                    # Vertex position mean value set to (0, 0, 0, 0)
x = 5                     # like the diameter of collimator
t = 2.5                   # ns
percentile = 0.95         # Ratio in diameter range >= 0.95
z = norm.ppf(percentile)  # z=\frac{x-\mu}{\sigma}
sigma_p = (x-mu) / z      # standard deviation for position
sigma_t = (t-mu) / z      # standard deviation for time
print(f"standard deviation of position gaussian distribution: {sigma_p} mm")
print(f"standard deviation of time gaussian distribution: {sigma_t} ns")
# Phi
# 
phi_min = 0
tan_phi_max = ((plane_size_x/2)-x) / (plane_spacing*5 + 30) 
phi_max = math.degrees(math.atan(tan_phi_max))
print(f"phi range: from {phi_min} to {phi_max} degree")

# Eta
theta_min = 90 - phi_max
theta_max = 90 + phi_max
eta_min = -math.log(math.tan(math.radians(theta_min) / 2))
eta_max = -math.log(math.tan(math.radians(theta_max) / 2))
print(f"eta range: from {eta_min} to {eta_max}")

