import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, 10, 100)
good = np.sin(x) + 1
bad  = np.sin(x) + 1
bad[50:] += 20   # outlier part

threshold = 3

plt.plot(x, good, label="Good Algo")

# bad algo below threshold = solid
mask = bad <= threshold
plt.plot(x[mask], bad[mask], label="Bad Algo (normal)", color="red")

# bad algo above threshold = dotted
mask = bad > threshold
plt.plot(x[mask], bad[mask], linestyle=":", color="red", label="Bad Algo (outlier)")

plt.legend()
plt.show()


fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6,4),
                               gridspec_kw={'height_ratios': [1, 3]})

# Top axis (zoomed in)
ax2.plot(x, good, label="Good Algo")
ax2.plot(x, bad, label="Bad Algo")
ax2.set_ylim(0, 3.5)

# Bottom axis (full range)
ax1.plot(x, good, label="Good Algo")
ax1.plot(x, bad, label="Bad Algo")
ax1.set_ylim(15, 25)

# Hide spines between axes
ax1.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax1.tick_params(labeltop=False)  # no tick labels on top subplot
ax2.xaxis.tick_bottom()

plt.legend()
plt.show()

fig, ax1 = plt.subplots()

ax2 = ax1.twinx()  # second y-axis

ax1.plot(x, good, label="Good Algo", color="blue")
ax2.plot(x, bad, label="Bad Algo", color="red", linestyle=":")

ax1.set_ylabel("Good Scale")
ax2.set_ylabel("Bad Scale")
