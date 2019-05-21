import math
import seaborn as sns
import pandas
import matplotlib.pyplot as plt
import numpy as np

LIMIT_TO_CONSIDER_DONE = 0.0005

RESOLUTION = 1
MAX_CONCENTRATION = 1

CH3COO = "Acetate"
CH3COOH = "Acetic acid"
OH = "Hydroxide"
H3O = "Hydronium"

KW = 10.0 ** -14
KA = 1.76 * (10.0 ** -5)
KB = 5.68 * (10.0 ** -10)


def get_pH(initial_oh):
    iterations = []
    concentrations = []
    molecules = []

    def add_row(index, value, molecule):
        iterations.append(index)
        concentrations.append(value)
        molecules.append(molecule)

    oh = initial_oh
    ch3coo = 0.
    ch3cooh = 1.
    h3o = 10.0 ** (-7)
    add_row(0, ch3coo, CH3COO)
    add_row(0, ch3cooh, CH3COOH)
    add_row(0, oh, OH)
    add_row(0, h3o, H3O)
    i = 1
    while i < 1000:

        initial_ch3cooh = ch3cooh

        for _ in range(RESOLUTION):
            # q_w = h3o * oh
            # if ch3cooh == 0:
            #     q_a = 10 ** 10
            # else:
            #     q_a = h3o * ch3coo / ch3cooh
            # if ch3coo == 0:
            #     q_b = 10 ** 10
            # else:
            #     q_b = ch3cooh * oh / ch3coo

            change_w = 0.5 * (- oh - h3o + math.sqrt((oh + h3o) ** 2 - 4 * oh * h3o + 4 * KW))

            if change_w < 0:
                change_w = -min(oh, h3o, -change_w)

            # print("change w" + str(change_w))
            oh += change_w
            h3o += change_w

            change_a = 0.5 * (- h3o - ch3coo - KA + math.sqrt(
                ((h3o + ch3coo + KA) ** 2) - 4 * h3o * ch3coo + 4 * KA * ch3cooh))

            if change_a > 0:
                change_a = min(change_a, ch3cooh)
            else:
                change_a = - min(-change_a, ch3coo, h3o)

            # print("change a" + str(change_a))
            ch3cooh -= change_a
            ch3coo += change_a
            h3o += change_a

            change_b = 0.5 * (- oh - ch3cooh - KB + math.sqrt(
                ((oh + ch3cooh + KB) ** 2) - 4 * oh * ch3cooh + 4 * KB * ch3coo))

            if change_b > 0:
                change_b = min(change_b, ch3coo)
            else:
                change_b = - min(-change_b, ch3cooh, oh)

            # print("change b" + str(change_b))

            ch3coo -= change_b
            oh += change_b
            ch3cooh += change_b

        add_row(i, ch3coo, CH3COO)
        add_row(i, ch3cooh, CH3COOH)
        add_row(i, oh, OH)
        add_row(i, h3o, H3O)

        if change_b < LIMIT_TO_CONSIDER_DONE and change_a < LIMIT_TO_CONSIDER_DONE and change_w < LIMIT_TO_CONSIDER_DONE:
            break

        i += 1
    df = pandas.DataFrame({"Iteration": iterations, "Concentration": concentrations, "Type": molecules})
    # print(df)
    # sns.set()
    # sns.lineplot(x='Iteration', y='Concentration', data=df, hue="Type")

    # print()
    h3o_concentration = df.loc[(df['Iteration'] == i) & (df['Type'] == H3O), "Concentration"]
    if h3o_concentration.item() == 0:
        # print("pH: infinity")
        return -1
    else:
        pH = -math.log10(h3o_concentration)
        # print(f"ph: {pH}")
        return pH


x_values = np.arange(0, 2, 0.01)
y_values = []

for c in x_values:
    y_values.append(get_pH(c))

sns.set()
sns.lineplot(x=x_values, y=y_values)

plt.show()
