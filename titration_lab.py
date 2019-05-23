import math
import seaborn as sns
import pandas
import matplotlib.pyplot as plt
import numpy as np

MAX_ITERATIONS = 1000
CHANGE_CUTOFF = 1 * (10 ** -20)
PH_CURVE_RESOLUTION = 0.01  # Greater number lower resolution

# Chemicals
NAOH = "Sodium hydroxide"
CH3COO = "Acetate"
CH3COOH = "Acetic acid"
OH = "Hydroxide"
H3O = "Hydronium"
NA = "Sodium ions"

# Reaction names
ACID = "Acetic acid dissociation"
BASE = "Base dissociation"
ACID_BASE = "Acid Base Reaction"
WATER = "Water auto-ionization"

# Equilibrium constants
KW = 10.0 ** -14
K_ACID = 1.76 * (10.0 ** -5)
K_ACID_BASE = K_ACID / KW  # Because the equations can subtract
K_BASE = 0.63

# Reactions
# First parameter is the rearranged Keq formula to isolate for the change variable (ICE table)
# Second parameter is what chemicals form the reaction. True is a product. False is a reactant.
REACTIONS = {
    ACID: (lambda chemicals: (- chemicals[H3O] - chemicals[CH3COO] - K_ACID + math.sqrt(
        ((chemicals[H3O] + chemicals[CH3COO] + K_ACID) ** 2) - 4 * (
                chemicals[H3O] * chemicals[CH3COO] - K_ACID * chemicals[CH3COOH]))) / 2,
           [(CH3COOH, False), (H3O, True), (CH3COO, True)]),
    BASE: (lambda chemicals: (- chemicals[OH] - chemicals[NA] - K_BASE + math.sqrt(
        ((chemicals[OH] + chemicals[NA] + K_BASE) ** 2) - 4 * (
                chemicals[NA] * chemicals[OH] - K_BASE * chemicals[NAOH]))) / 2,
           [(NAOH, False), (OH, True), (NA, True)]),
    WATER: (lambda chemicals: (- chemicals[OH] - chemicals[H3O] + math.sqrt(
        (chemicals[OH] + chemicals[H3O]) ** 2 - 4 * chemicals[OH] * chemicals[H3O] + 4 * KW)) / 2,
            [(OH, True), (H3O, True)]),
    ACID_BASE: (lambda chemicals: 2 / (- chemicals[OH] - chemicals[CH3COOH] - K_ACID_BASE + math.sqrt(
        ((chemicals[OH] + chemicals[CH3COOH] + K_ACID_BASE) ** 2) - 4 * (
                chemicals[OH] * chemicals[CH3COOH] - K_ACID_BASE * chemicals[CH3COO]))),
                [(CH3COOH, False), (OH, False), (CH3COO, True)])
}

# NaOH is not included because it is the variable (not constant)
INITIAL_CONCENTRATIONS = {
    OH: 10.0 ** (-7),
    CH3COO: 0,
    CH3COOH: 1,
    H3O: 10.0 ** (-7),
    NA: 0
}

sns.set()  # Set the style of the graphs


class LinearGraph:
    def __init__(self, x_axis_name, y_axis_name, series_name):
        self.x_axis_name = x_axis_name
        self.y_axis_name = y_axis_name
        self.series_name = series_name

        self.x_values = []
        self.y_values = []
        self.series = []

        self.database = None

    def add_data_point(self, x_value, y_value, series):
        self.x_values.append(x_value)
        self.y_values.append(y_value)
        self.series.append(series)

    def add_data_for_each_series(self, x_value, series_with_data):
        for (series, value) in series_with_data.items():
            self.add_data_point(x_value, value, series)

    def _create_database(self):
        self.database = pandas.DataFrame(
            {self.x_axis_name: self.x_values, self.y_axis_name: self.y_values, self.series_name: self.series})

    def find_y_value(self, x_value, series):
        if self.database is None:
            self._create_database()

        return self.database.loc[(self.database[self.x_axis_name] == x_value) & (
                self.database[self.series_name] == series), self.y_axis_name].values[0]

    def graph(self, logarithmic_y_axis=False, title=None):
        if self.database is None:
            self._create_database()

        plot = sns.lineplot(x=self.x_axis_name, y=self.y_axis_name, hue=self.series_name, data=self.database)

        if logarithmic_y_axis:
            plt.yscale('symlog', linthreshy=10 ** -16)  # Small linthreshy required to make function work

        if title is not None:
            plot.set_title(title)

        plt.show()


def calculate_ph(hydronium_concentration):
    return -math.log10(hydronium_concentration)


def calculate_equilibrium(initial_naoh_concentration):
    concentrations = INITIAL_CONCENTRATIONS.copy()
    concentrations[NAOH] = initial_naoh_concentration

    total_changes = {}

    for name in REACTIONS.keys():
        total_changes[name] = 0

    # concentration_graph = Graph("Iterations", "Concentration (mol/L)", "Molecule")
    # change_in_concentration_graph = Graph("Iterations", "Change in concentration", "Reaction")
    # concentration_graph.add_data_for_each_series(0, chemicals)

    for i in range(1, MAX_ITERATIONS):  # Starts at 1 because when we use graphs there's already data at 0
        change_is_minimal = True  # Default True. Set to False on a big change

        for (name, (equation, changes)) in REACTIONS.items():
            required_change = equation(concentrations)

            # Limit the change if there is not enough reactant or product
            if required_change != 0:
                for (chemical, positive) in changes:
                    if positive:
                        if required_change < 0:
                            required_change = - min(-required_change, concentrations[chemical])
                    else:
                        if required_change > 0:
                            required_change = min(required_change, concentrations[chemical])

                required_change /= 2  # So that all reactions have a chance to contribute to the equilibrium

                # Apply the required change
                for (chemical, positive) in changes:
                    if positive:
                        concentrations[chemical] += required_change
                    else:
                        concentrations[chemical] -= required_change

                total_changes[name] += required_change

                if required_change > CHANGE_CUTOFF:
                    change_is_minimal = False

            # change_in_concentration_graph.add_data_point(i, required_change, equation)
        # concentration_graph.add_data_for_each_series(i, chemicals)

        if change_is_minimal:
            break

    return concentrations, total_changes


def main():
    ph_graph = LinearGraph("Initial concentration of NaOH (mol/L)", "pH", "")
    # change_graph = LinearGraph("Volume of NaOH added (mol/L)", "Change from equations (mol/L)", "Reaction")

    for c in np.arange(0, 2, PH_CURVE_RESOLUTION):
        (chemicals, total_changes) = calculate_equilibrium(c)
        ph_graph.add_data_for_each_series(c, {"pH": calculate_ph(chemicals[H3O])})
        # change_graph.add_data_for_each_series(c, total_changes)

    plt.ylim(0, 14)
    ph_graph.graph(title="pH as sodium hydroxide is added to 1M acetic acid")

    # change_graph.graph(logarithmic_y_axis=True)
    # change_graph.graph(logarithmic_y_axis=False)


if __name__ == "__main__":
    main()
