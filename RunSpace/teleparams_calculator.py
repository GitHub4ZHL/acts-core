#!/usr/bin/env python3
import math
from scipy.stats import norm
from scipy.stats import multivariate_normal
import numpy as np
import matplotlib.pyplot as plt


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

# Draw
n_points = 100
radius = 6
# position
x_p = np.linspace(-radius, radius, n_points)
y_p = np.linspace(-radius, radius, n_points)
X_p, Y_p = np.meshgrid(x_p, y_p)
points = np.column_stack((X_p.ravel(), Y_p.ravel()))

mean = np.array([0, 0])
covariance_p = np.array([[sigma_p**2, 0], [0, sigma_p**2]])

pdf_values_p = multivariate_normal.pdf(points, mean=mean, cov=covariance_p)
pdf_values_p = pdf_values_p.reshape(n_points, n_points)

plt.figure(figsize=(8, 8))
plt.imshow(pdf_values_p, extent=[-radius, radius, -radius, radius], origin='lower', cmap='coolwarm')
plt.title('Position Distribution')
plt.xlabel('mm')
plt.ylabel('mm')
plt.grid(False)
#plt.savefig("position-distribution.png", dpi=300)
plt.savefig("position-distribution.pdf", dpi=300)

# time
x_t = np.linspace(-radius, radius, n_points)
y_t = np.linspace(-radius, radius, n_points)
X_t, Y_t = np.meshgrid(x_t, y_t)
points = np.column_stack((X_t.ravel(), Y_t.ravel()))
covariance_t = np.array([[sigma_t**2, 0], [0, sigma_t**2]])

pdf_values_t = multivariate_normal.pdf(points, mean=mean, cov=covariance_t)
pdf_values_t = pdf_values_t.reshape(n_points, n_points)

plt.figure(figsize=(8, 8))
plt.imshow(pdf_values_t, extent=[-radius, radius, -radius, radius], origin='lower', cmap='coolwarm')
plt.title('Time Distribution')
plt.xlabel('ns')
plt.ylabel('ns')
plt.grid(False)
#plt.savefig("time-distribution.png", dpi=300)
plt.savefig("time-distribution.pdf", dpi=300)
