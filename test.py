import rftpy as rft

radius = rft.Earth.radius
mu = rft.Earth.mu

v = rft.vis_viva(mu, radius, a=8000)

print(f"{v: 1.5f}")
