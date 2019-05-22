import math
import seaborn as sns
import pandas
import matplotlib.pyplot as plt
import numpy as np

MAX_ITERATIONS = 1000

CHANGE_CUTOFF = 1 * (10 ** -20)

PH_CURVE_RESOLUTION = 0.01  # Greater number lower resolution
CONCENTRATION_CURVE_RESOLUTION = 1  # Greater number lower resolution

NAOH = "Sodium hydroxide"
CH3COO = "Acetate"
CH3COOH = "Acetic acid"
OH = "Hydroxide"
H3O = "Hydronium"
NA = "Sodium ions"

ACID = "Acetic acid dissociation"
BASE = "Base dissociation"
ACID_BASE = "Acid Base Reaction"
WATER = "Water auto-ionization"

KW = 10.0 ** -14
K_ACID = 1.76 * (10.0 ** -5)
K_ACID_BASE = KW / K_ACID
K_BASE = 0.63


class Database:
    def __init__(self, x_axis_name, y_axis_name):
        self.x_axis_name = x_axis_name
        self.y_axis_name = y_axis_name
        self.x_values = []
        self.values = []
        self.series = []
        self.df = None

    def add_row(self, x_value, value, series):
        self.x_values.append(x_value)
        self.values.append(value)
        self.series.append(series)

    def get_pandas_database(self):
        self.df = pandas.DataFrame(
            {self.x_axis_name: self.x_values, self.y_axis_name: self.values, "Type": self.series})

    def get_final_value_of_series(self, last_x_value, series):
        if self.df is None:
            self.get_pandas_database()

        return self.df.loc[(self.df[self.x_axis_name] == last_x_value) & (
                self.df['Type'] == series), self.y_axis_name].values[0]

    def graph(self, logarithmic=False):
        if self.df is None:
            self.get_pandas_database()

        sns.lineplot(x=self.x_axis_name, y=self.y_axis_name, hue="Type", data=self.df)

        if logarithmic:
            plt.yscale('symlog', linthreshy=10 ** -16)  # Small linthreshy required to make function work

    def add_point_in_data(self, x_value, new_data):
        for row in new_data.keys():
            self.add_row(x_value, new_data[row], row)


def calculate_ph(hydronium_concentration):
    return -math.log10(hydronium_concentration)


def get_ph(initial_naoh):
    database = Database("Iterations", "Concentration (mol/L)")
    change_database = Database("Iterations", "Change")

    chemicals = {
        OH: 10.0 ** (-7),
        CH3COO: 0.,
        CH3COOH: 1.,
        H3O: 10.0 ** (-7),
        NAOH: initial_naoh,
        NA: 0.
    }

    equations = {
        WATER: [lambda: (- chemicals[OH] - chemicals[H3O] + math.sqrt(
            (chemicals[OH] + chemicals[H3O]) ** 2 - 4 * chemicals[OH] * chemicals[H3O] + 4 * KW)) / 2, 0., [(OH, True),
                                                                                                            (H3O,
                                                                                                             True)]],
        ACID: [lambda: (- chemicals[H3O] - chemicals[CH3COO] - K_ACID + math.sqrt(
            ((chemicals[H3O] + chemicals[CH3COO] + K_ACID) ** 2) - 4 * (
                    chemicals[H3O] * chemicals[CH3COO] - K_ACID * chemicals[CH3COOH]))) / 2, 0., [(CH3COOH, False),
                                                                                                  (H3O, True),
                                                                                                  (CH3COO, True)]],
        BASE: [lambda: (- chemicals[OH] - chemicals[NA] - K_BASE + math.sqrt(
            ((chemicals[OH] + chemicals[NA] + K_BASE) ** 2) - 4 * (
                    chemicals[NA] * chemicals[OH] - K_BASE * chemicals[NAOH]))) / 2, 0., [(NAOH, False),
                                                                                          (OH, True),
                                                                                          (NA, True)]],
        ACID_BASE: [lambda: (- chemicals[OH] - chemicals[CH3COOH] - K_ACID_BASE + math.sqrt(
            ((chemicals[OH] + chemicals[CH3COOH] + K_ACID_BASE) ** 2) - 4 * (
                    chemicals[OH] * chemicals[CH3COOH] - K_ACID_BASE * chemicals[CH3COO]))) / 2, 0., [(CH3COO, False),
                                                                                                      (CH3COOH, True),
                                                                                                      (OH, True)]]
    }

    database.add_point_in_data(0, chemicals)

    next_entry = 1
    for i in range(1, MAX_ITERATIONS):
        print()
        is_done = True
        for equation in equations.keys():
            required_change = equations[equation][0]()

            # Limit the change if there is not enough reactant
            if required_change != 0:
                print(f"{equation}   {required_change}")
                for (chemical, positive) in equations[equation][2]:
                    if positive:
                        if required_change < 0:
                            required_change = - min(-required_change, chemicals[chemical])
                    else:
                        if required_change > 0:
                            required_change = min(required_change, chemicals[chemical])

                # Apply the required change
                for (chemical, positive) in equations[equation][2]:
                    if positive:
                        chemicals[chemical] += required_change
                    else:
                        chemicals[chemical] -= required_change

                equations[equation][1] += required_change

                if required_change > CHANGE_CUTOFF:
                    is_done = False

            change_database.add_row(i, required_change, equation)

        if i % CONCENTRATION_CURVE_RESOLUTION == 0:
            database.add_point_in_data(next_entry, chemicals)
            next_entry += 1

        if is_done:
            print(i)
            break

    ph = calculate_ph(database.get_final_value_of_series(next_entry - 1, H3O))

    total_changes = {}

    for equation_type in equations.keys():
        total_changes[equation_type] = equations[equation_type][1]

    return ph, total_changes


def main():
    naoh_values = np.arange(0, 2, PH_CURVE_RESOLUTION)

    ph_database = Database("Volume of NaOH added (mol/L)", "pH")
    change_database = Database("Volume of NaOH added (mol/L)", "Change from equations (mol/L)")

    for c in naoh_values:
        (ph, total_changes) = get_ph(c)
        ph_database.add_point_in_data(c, {"pH": ph})
        change_database.add_point_in_data(c, total_changes)

    sns.set()

    ph_database.graph()
    plt.ylim(0, 14)
    plt.show()

    change_database.graph(logarithmic=True)
    plt.show()


if __name__ == "__main__":
    main()
