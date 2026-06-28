# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 09:53:11 2026

@author: Achim
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# -----------------------------
# Parameter
# -----------------------------
L = 10.0                 # Seitenlänge der Fläche
n_steps = 500            # Anzahl Zeitschritte
human_step = 0.2         # Schrittlänge des Menschen
dog_speed = 0.04         # Geschwindigkeit des Hundes

# -----------------------------
# Anfangspositionen
# -----------------------------
human = np.array([L/2, L/2])
dog = np.array([2.0, 2.0])

human_path = [human.copy()]
dog_path = [dog.copy()]

# -----------------------------
# Simulation
# -----------------------------
for _ in range(n_steps):

    # Random Walk des Menschen
    angle = np.random.uniform(0, 2*np.pi)
    step = human_step * np.array([np.cos(angle), np.sin(angle)])
    human = (human + step) % L

    # Kürzester Verbindungsvektor wegen periodischer Randbedingungen
    delta = human - dog
    delta -= L * np.round(delta / L)

    dist = np.linalg.norm(delta)
    if dist > 1e-10:
        dog += dog_speed * delta / dist

    dog %= L

    human_path.append(human.copy())
    dog_path.append(dog.copy())

human_path = np.array(human_path)
dog_path = np.array(dog_path)

# -----------------------------
# Animation
# -----------------------------
fig, ax = plt.subplots(figsize=(6,6))
ax.set_xlim(0, L)
ax.set_ylim(0, L)
ax.set_aspect("equal")
ax.axis('off')

human_line, = ax.plot([], [], 'b-', lw=1, alpha=0.4)
dog_line, = ax.plot([], [], 'r-', lw=1, alpha=0.4)
connection_line, = ax.plot([], [], color='black', linestyle='dashed', lw=0.5, alpha=0.4)

human_dot, = ax.plot([], [], 'bo', label="human")
dog_dot, = ax.plot([], [], 'ro', label="dog")

ax.legend()


def init():
    human_line.set_data([], [])
    dog_line.set_data([], [])
    human_dot.set_data([], [])
    dog_dot.set_data([], [])
    return human_line, dog_line, human_dot, dog_dot


def update(frame):
    human_line.set_data(human_path[:frame+1,0], human_path[:frame+1,1])
    dog_line.set_data(dog_path[:frame+1,0], dog_path[:frame+1,1])
    connection_line.set_data([human_path[frame+1,0],dog_path[frame+1,0]],
                             [human_path[frame+1,1],dog_path[frame+1,1]])

    human_dot.set_data([human_path[frame,0]], [human_path[frame,1]])
    dog_dot.set_data([dog_path[frame,0]], [dog_path[frame,1]])

    ax.set_title("Pursuit Curve")

    return human_line, dog_line, human_dot, dog_dot


ani = FuncAnimation(
    fig,
    update,
    frames=n_steps,
    init_func=init,
    interval=30,
    blit=True
)

ani.save('Hunde5.gif',fps=10,dpi=300)

plt.show()