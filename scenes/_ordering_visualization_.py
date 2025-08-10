import numpy as np
import matplotlib.pyplot as plt

# Parameters
r_max = 2
gradient = np.array([3.0, 1.0])
gradient /= np.linalg.norm(gradient)

# Generate neighbor offsets
neighbors = []
for dx in range(-r_max, r_max + 1):
    for dy in range(-r_max, r_max + 1):
        if dx == 0 and dy == 0:
            continue
        neighbors.append(np.array([dx, dy]))

#cosine similarity
def similarity(o):
    if o[0] == round(gradient[0]) and o[1] == round(gradient[1]):
        return 2
    return np.dot(o, gradient) / np.linalg.norm(o)

# Sort: first by descending similarity, then by ascending distance
neighbors = sorted(neighbors, key=lambda o: (-similarity(o), np.linalg.norm(o)))

# Plot
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_aspect('equal', adjustable='box')

# Draw grid lines
for x in range(-r_max, r_max + 2):
    ax.axvline(x - 0.5, color='lightgray', linewidth=1)
for y in range(-r_max, r_max + 2):
    ax.axhline(y - 0.5, color='lightgray', linewidth=1)

# Plot cells with center positions
order_labels = {}
for order, offset in enumerate(neighbors, start=1):
    cx, cy = offset  # cell center coords
    ax.scatter(cx, cy, c='lightblue', s=500, edgecolor='k')
    order_labels[(cx, cy)] = str(order)

# Plot origin cell
ax.scatter(0, 0, c='orange', s=500, edgecolor='k', zorder=3)
order_labels[(0, 0)] = "Start"

# Add labels
for (cx, cy), label in order_labels.items():
    ax.text(cx, cy, label, ha='center', va='center', fontsize=10, fontweight='bold')

# Draw gradient vector
ax.arrow(0, 0, gradient[0], gradient[1], head_width=0.15, head_length=0.2, fc='r', ec='r', linewidth=2)
ax.text(gradient[0] * 1.1, gradient[1] * 1.1, r"$\hat{\mathbf{g}}$", color='r', fontsize=14)

ax.set_xlim(-r_max - 1, r_max + 1)
ax.set_ylim(-r_max - 1, r_max + 1)
ax.axis('off')
plt.show()


