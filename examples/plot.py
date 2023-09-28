import numpy as np
import adaptiveDesignProcedure as adp
import matplotlib.pyplot as plt

path = "plams_workdir/adp.results/"

x_CO_model = np.linspace(0.001,0.999,201)
TOF_CO2_model, = adp.predict( x_CO_model.reshape(-1,1), path ).T

ax = plt.axes()
ax.set_xlabel('Molar Fraction CO', fontsize=14)
ax.set_ylabel("TOF (mol/s/site)", fontsize=14)
ax.plot(x_CO_model, TOF_CO2_model, color='red', linestyle='-', lw=2, zorder=0)
plt.tight_layout()
plt.savefig("apd.png")

