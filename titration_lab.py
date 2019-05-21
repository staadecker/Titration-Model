import math
import seaborn as sns
import pandas
import matplotlib.pyplot as plt
import numpy as np

MAX_ITERATIONS = 1000

CHANGE_CUTOFF = 0.0005

PH_CURVE_RESOLUTION = 0.01  # Greater number lower resolution
CONCENTRATION_CURVE_RESOLUTION = 1  # Greater number lower resolution

CH3COO = "Acetate"
CH3COOH = "Acetic acid"
OH = "Hydroxide"
H3O = "Hydronium"

KW = 10.0 ** -14
KA = 1.76 * (10.0 ** -5)
KB = KW / KA


def get_ph(initial_oh):
    iterations = []
    concentrations = []
    molecules = []

    def add_row(index, value, molecule):
        iterations.append(index)
        concentrations.append(value)
        molecules.append(molecule)

    oh = initial_oh + 10.0 ** (-7)
    ch3coo = 0.
    ch3cooh = 1.
    h3o = 10.0 ** (-7)

    add_row(0, ch3coo, CH3COO)
    add_row(0, ch3cooh, CH3COOH)
    add_row(0, oh, OH)
    add_row(0, h3o, H3O)

    next_entry = 1

    for i in range(1, MAX_ITERATIONS):
        # WATER AUTO-IONIZATION
        change_w = 0.5 * (- oh - h3o + math.sqrt((oh + h3o) ** 2 - 4 * oh * h3o + 4 * KW))

        if change_w < 0:
            change_w = -min(oh, h3o, -change_w)

        oh += change_w
        h3o += change_w

        # ACETIC ACID DISSOCIATION
        change_a = (- h3o - ch3coo - KA + math.sqrt(((h3o + ch3coo + KA) ** 2) - 4 * (h3o * ch3coo - KA * ch3cooh))) / 2

        if change_a > 0:
            change_a = min(change_a, ch3cooh)
        else:
            change_a = - min(-change_a, ch3coo, h3o)

        ch3cooh -= change_a
        ch3coo += change_a
        h3o += change_a

        #  ACETATE BASE
        change_b = (- oh - ch3cooh - KB + math.sqrt(((oh + ch3cooh + KB) ** 2) - 4 * (oh * ch3cooh - KB * ch3coo))) / 2

        if change_b > 0:
            change_b = min(change_b, ch3coo)
        else:
            change_b = - min(-change_b, ch3cooh, oh)

        ch3coo -= change_b
        oh += change_b
        ch3cooh += change_b

        if i % CONCENTRATION_CURVE_RESOLUTION == 0:
            add_row(next_entry, ch3coo, CH3COO)
            add_row(next_entry, ch3cooh, CH3COOH)
            add_row(next_entry, oh, OH)
            add_row(next_entry, h3o, H3O)
            next_entry += 1

        if change_b < CHANGE_CUTOFF and change_a < CHANGE_CUTOFF and change_w < CHANGE_CUTOFF:
            break

    df = pandas.DataFrame({"Iteration": iterations, "Concentration": concentrations, "Type": molecules})

    # sns.set()
    # sns.lineplot(x='Iteration', y='Concentration', data=df, hue="Type")

    h3o_concentration = df.loc[(df['Iteration'] == (next_entry - 1)) & (df['Type'] == H3O), "Concentration"]

    return -math.log10(h3o_concentration)


x_values = np.arange(0, 2, PH_CURVE_RESOLUTION)
ph_values = []

for c in x_values:
    ph_values.append(get_ph(c))

ph_df = pandas.DataFrame({"pH": ph_values, "Initial concentration of OH (mol/L)": x_values})

sns.set()
plot = sns.lineplot(x="Initial concentration of OH (mol/L)", y="pH", data=ph_df)
plot.axes.set_ylim(1, 14)

plt.show()
