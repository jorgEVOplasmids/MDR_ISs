import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import seaborn as sns
from scipy.optimize import minimize
from pprint import pprint
from matplotlib.patches import Wedge
import time
import pickle
import os
from scipy.stats import entropy
import networkx as nx
import matplotlib.colors as mcolors
import matplotlib.lines as mlines


def load_strain_parameters(strainIDs, pathPARAMS):
    """
    Load plasmid-free and plasmid-bearing strain parameters from .pkl files.

    Parameters:
    strainIDs (list): List of strain identifiers (e.g., ["K253", "K168"]).
    pathPARAMS (str): Path to the folder containing parameter .pkl files.

    Returns:
    dict: Dictionary of strain parameters with structure:
          strains_species[species_idx]["0"] → plasmid-free
          strains_species[species_idx]["p"] → plasmid-bearing
    """
    strains_species = {}
    for species_idx, strainID in enumerate(strainIDs):
        file_0 = f'params_{strainID}_0.pkl'
        file_p = f'params_{strainID}_p.pkl'

        with open(os.path.join(pathPARAMS, file_0), 'rb') as f:
            strains_0 = pickle.load(f)
        print(f"Plasmid-free parameters for strain {strainID} loaded successfully from:\n{os.path.join(pathPARAMS, file_0)}\n")

        with open(os.path.join(pathPARAMS, file_p), 'rb') as f:
            strains_p = pickle.load(f)
        print(f"Plasmid-bearing parameters for strain {strainID} loaded successfully from:\n{os.path.join(pathPARAMS, file_p)}\n")

        strains_species[species_idx] = {
            "0": strains_0,
            "p": strains_p
        }

    return strains_species


def gillespie_bacterial_growth_multi_species(
    strains, plasmid_matrix, populations_0, populations_p, antibiotic_concentration,
    initial_resource, simulation_time, dt=0.1, verbose=False
):
    """
    Simulates bacterial growth for multiple species with plasmid-free and plasmid-bearing populations.
    Includes resource competition and inter-species plasmid conjugation.

    Parameters:
    strains (dict): Dictionary with parameters for plasmid-free ("0") and plasmid-bearing ("p") strains for each species.
    populations_0 (np.array): 3D array of plasmid-free populations [species][SNP][IS].
    populations_p (np.array): 3D array of plasmid-bearing populations [species][SNP][IS].
    antibiotic_concentration (float): Antibiotic concentration in the environment.
    initial_resource (float): Initial shared resource concentration.
    simulation_time (float): Total simulation time.
    dt (float): Time step for the simulation.

    Returns:
    tuple: Time values, population values, and resource values during the simulation.
    """
    t = 0
    R = initial_resource  # Shared resource pool
    num_species, num_mutationsSNP, num_mutationsIS = populations_0.shape

    # Tracking results
    t_values = [t]
    population_values = [(populations_0.copy(), populations_p.copy())]
    resource_values = [R]

    while t < simulation_time and R > 0 and (
        np.sum(populations_0) > 0 or np.sum(populations_p) > 0
    ):
        if verbose:
            print(f"--- Time {t:.2f} ---")
            print(f"Resource level: {R:.2f}")
            print(f"Total plasmid-free population: {np.sum(populations_0)}")
            print(f"Total plasmid-bearing population: {np.sum(populations_p)}")

        for species in range(num_species):
            if verbose:
                print(f"\nProcessing species {species}...")

            # Initialize counters for events
            births, deaths, mutations, transpositions, segregations, conjugations = 0, 0, 0, 0, 0, 0

            # Process plasmid-free populations
            for k, strain in enumerate(strains[species]["0"]):
                i, j = divmod(k, num_mutationsIS)
                half_saturation_resource = strain['half_saturation_resource']
                half_saturation_antibiotic = strain['half_saturation_antibiotic']

                # Effective birth rate
                birth_rate_eff = strain['birth_rate'] * (R / (R + half_saturation_resource))
                birth_rate_eff = max(birth_rate_eff, 0)

                # Expected births
                expected_births = (
                    max(0, np.random.poisson(birth_rate_eff * populations_0[species, i, j] * dt))
                    if populations_0[species, i, j] > 0 else 0
                )
                births += expected_births

                # Effective death rate
                death_rate_eff = strain['death_rate'] * (
                    antibiotic_concentration / (antibiotic_concentration + half_saturation_antibiotic)
                )
                death_rate_eff = max(death_rate_eff, 0)

                # Expected deaths
                expected_deaths = (
                    max(0, np.random.poisson(death_rate_eff * populations_0[species, i, j] * dt))
                    if populations_0[species, i, j] > 0 else 0
                )
                deaths += expected_deaths

                # Update population
                populations_0[species, i, j] = max(
                    0, populations_0[species, i, j] + expected_births - expected_deaths
                )
                populations_0[species, i, j] = (
                    populations_0[species, i, j] if populations_0[species, i, j] >= 1 else 0
                )

                # Resource consumption
                R = max(0, R - expected_births * strain['consumption_rate'])

                # Mutations and transpositions
                if strain['mutation_rate'] > 0 and i < num_mutationsSNP - 1:
                    mutation_events = np.random.poisson(strain['mutation_rate'] * populations_0[species, i, j])
                    populations_0[species, i, j] = max(0, populations_0[species, i, j] - mutation_events)
                    populations_0[species, i + 1, j] += mutation_events
                    mutations += mutation_events

                if strain['transposition_rate'] > 0 and j < num_mutationsIS - 1:
                    transposition_events = np.random.poisson(strain['transposition_rate'] * populations_0[species, i, j])
                    populations_0[species, i, j] = max(0, populations_0[species, i, j] - transposition_events)
                    populations_0[species, i, j + 1] += transposition_events
                    transpositions += transposition_events

            # Process plasmid-bearing populations
            for k, strain in enumerate(strains[species]["p"]):
                i, j = divmod(k, num_mutationsIS)
                half_saturation_resource = strain['half_saturation_resource']
                half_saturation_antibiotic = strain['half_saturation_antibiotic']

                # Effective birth rate
                birth_rate_eff = strain['birth_rate'] * (R / (R + half_saturation_resource))
                birth_rate_eff = max(birth_rate_eff, 0)

                # Expected births
                expected_births = (
                    max(0, np.random.poisson(birth_rate_eff * populations_p[species, i, j] * dt))
                    if populations_p[species, i, j] > 0 else 0
                )
                births += expected_births

                # Effective death rate
                death_rate_eff = strain['death_rate'] * (
                    antibiotic_concentration / (antibiotic_concentration + half_saturation_antibiotic)
                )
                death_rate_eff = max(death_rate_eff, 0)

                # Expected deaths
                expected_deaths = (
                    max(0, np.random.poisson(death_rate_eff * populations_p[species, i, j] * dt))
                    if populations_p[species, i, j] > 0 else 0
                )
                deaths += expected_deaths

                # Update population
                populations_p[species, i, j] = max(
                    0, populations_p[species, i, j] + expected_births - expected_deaths
                )
                populations_p[species, i, j] = (
                    populations_p[species, i, j] if populations_p[species, i, j] >= 1 else 0
                )

                # Resource consumption
                R = max(0, R - expected_births * strain['consumption_rate'])

                # Mutations, transpositions, segregations, and conjugations
                if strain['mutation_rate'] > 0 and i < num_mutationsSNP - 1:
                    mutation_events = np.random.poisson(strain['mutation_rate'] * populations_p[species, i, j])
                    populations_p[species, i, j] = max(0, populations_p[species, i, j] - mutation_events)
                    populations_p[species, i + 1, j] += mutation_events
                    mutations += mutation_events

                if strain['transposition_rate'] > 0 and j < num_mutationsIS - 1:
                    transposition_events = np.random.poisson(strain['transposition_rate'] * populations_p[species, i, j])
                    populations_p[species, i, j] = max(0, populations_p[species, i, j] - transposition_events)
                    populations_p[species, i, j + 1] += transposition_events
                    transpositions += transposition_events

                # Segregations
                segregation_events = np.random.poisson(strain['segregation_rate'] * expected_births)
                populations_p[species, i, j] = max(0, populations_p[species, i, j] - segregation_events)
                populations_0[species, i, j] += segregation_events
                segregations += segregation_events

                # Inter-species conjugations using plasmid transmission matrix
                for target_species in range(num_species):
                    if plasmid_matrix[species, target_species] > 0:  # Ensure conjugation is allowed
                        for m in range(num_mutationsSNP):
                            for n in range(num_mutationsIS):
                                #print(f"Species {species} → Target Species {target_species}, Conjugation Rate: {plasmid_matrix[species, target_species]}")

                                # Use the transmission matrix value instead of fixed conjugation rate
                                conjugation_events = np.random.poisson(
                                    plasmid_matrix[species, target_species] * 
                                    populations_p[species, i, j] * populations_0[target_species, m, n] * dt
                                )

                                # Update populations
                                populations_0[target_species, m, n] = max(0, populations_0[target_species, m, n] - conjugation_events)
                                populations_p[target_species, m, n] += conjugation_events
                                conjugations += conjugation_events


            # Print event statistics for the species
            if verbose:
                print(f"Species {species} summary:")
                print(f"  Births: {births}")
                print(f"  Deaths: {deaths}")
                print(f"  Mutations: {mutations}")
                print(f"  Transpositions: {transpositions}")
                print(f"  Segregations: {segregations}")
                print(f"  Conjugations: {conjugations}")

        # Update time and store results
        t += dt
        t_values.append(t)
        population_values.append((populations_0.copy(), populations_p.copy()))
        resource_values.append(R)

    return t_values, population_values, resource_values



def runSimulationIS_multi_species(strains,plasmid_matrix, initial_populations_0, initial_populations_p, num_days, antibiotic_concentration, initial_resource, simulation_time, dilution=0.1, verbose=False):
    results = []

    # Initialize populations
    populations_0 = np.array(initial_populations_0, dtype=float).copy()
    populations_p = np.array(initial_populations_p, dtype=float).copy()

    for day in range(num_days):
        # Simulate one day
        time_points, population_values, resource_values = gillespie_bacterial_growth_multi_species(
            strains, plasmid_matrix, populations_0, populations_p, antibiotic_concentration, initial_resource, simulation_time, 0.1, verbose
        )

        # Save results
        day_results = {
            'day': day + 1,
            'time_points': time_points,
            'population_values': population_values,
            'resource_values': resource_values,
            'final_populations_0': population_values[-1][0],  # Last state of plasmid-free populations
            'final_populations_p': population_values[-1][1],  # Last state of plasmid-bearing populations
            'resource': resource_values[-1]
        }
        results.append(day_results)

        # Apply dilution for the next day
        populations_0 *= dilution
        populations_p *= dilution

        # Stop if all populations are extinct
        if np.all(populations_0 <= 0) and np.all(populations_p <= 0):
            print(f"All populations went to zero on day {day + 1}. Stopping simulation.")
            break

    return results

def initialize_populations(num_species, num_mutationsSNP, num_mutationsIS, strainIDs, initial_values_by_strain=None, default_initial_values={"0": 1e6, "p": 0}):
    """
    Initialize populations for multiple species based on strain-specific or default initial values.

    Parameters:
    - num_species (int): Number of species.
    - num_mutationsSNP (int): Number of SNP mutation levels.
    - num_mutationsIS (int): Number of IS transposition levels.
    - strainIDs (list): List of strain IDs.
    - initial_values_by_strain (dict): Initial values for each strain (optional).
    - default_initial_values (dict): Default initial values (optional).

    Returns:
    - populations_0 (np.array): Initial plasmid-free populations.
    - populations_p (np.array): Initial plasmid-bearing populations.
    """

    # Use default initial values if none are provided
    if default_initial_values is None:
        default_initial_values = {"0": 1e6, "p": 100}

    if initial_values_by_strain is None:
        initial_values_by_strain = {}

    # Initialize populations
    populations_0 = np.zeros((num_species, num_mutationsSNP, num_mutationsIS))  # Plasmid-free
    populations_p = np.zeros((num_species, num_mutationsSNP, num_mutationsIS))  # Plasmid-bearing

    # Assign initial populations
    for species_idx, strainID in enumerate(strainIDs):
        # Get initial values for the current strain, falling back to defaults
        strain_initial_values = initial_values_by_strain.get(strainID, default_initial_values)

        # Assign populations to the first mutation and transposition level
        populations_0[species_idx, 0, 0] = strain_initial_values.get("0", 0)
        populations_p[species_idx, 0, 0] = strain_initial_values.get("p", 0)

    return populations_0, populations_p


def runManySimulations(plasmid_matrix, num_reps, runSimulationIS_multi_species, strains, initial_populations_0, initial_populations_p, num_days, antibiotic_concentration, initial_resource, simulation_time, dilution_factor):
    """
    Runs the simulation multiple times and computes the mean of final populations and resource values across repetitions.

    Parameters:
    num_reps (int): Number of repetitions of the simulation.
    runSimulationIS_multi_species (function): The simulation function to run.
    strains (list): List of strain parameters.
    initial_populations_0 (list): Initial plasmid-free populations.
    initial_populations_p (list): Initial plasmid-bearing populations.
    num_days (int): Number of simulation days.
    antibiotic_concentration (float): Antibiotic concentration used in the simulation.
    initial_resource (float): Initial resource concentration.
    simulation_time (float): Total simulation time.
    dilution_factor (float): Dilution factor applied daily.

    Returns:
    list: Averaged results for each day, containing mean `final_populations_0`, `final_populations_p`, and `resource`.
    """
    # List to store results from all repetitions
    all_results = []

    # Run simulations and collect results
    for rep in range(num_reps):
        #print("Simulation %s"%(rep+1))
        print(".", end='')
        results = runSimulationIS_multi_species(
            strains=strains,
            plasmid_matrix=plasmid_matrix,
            initial_populations_0=initial_populations_0,
            initial_populations_p=initial_populations_p,
            num_days=num_days,
            antibiotic_concentration=antibiotic_concentration,
            initial_resource=initial_resource,
            simulation_time=simulation_time,
            dilution=dilution_factor,
            verbose=False
        )
        all_results.append(results)

    # Initialize mean results
    num_days = len(all_results[0])
    mean_results = []
    for day_idx in range(num_days):
        mean_results.append({
            'day': day_idx + 1,
            'final_populations_0': None,
            'final_populations_p': None,
            'resource': None
        })

    # Compute means for each day
    for day_idx in range(num_days):
        # Initialize accumulators
        sum_pop_0 = None
        sum_pop_p = None
        sum_resource = 0

        for rep in range(num_reps):
            day_result = all_results[rep][day_idx]

            if sum_pop_0 is None:
                sum_pop_0 = np.zeros_like(day_result['final_populations_0'])
            if sum_pop_p is None:
                sum_pop_p = np.zeros_like(day_result['final_populations_p'])

            sum_pop_0 += day_result['final_populations_0']
            sum_pop_p += day_result['final_populations_p']
            sum_resource += day_result['resource']

        # Compute averages
        mean_results[day_idx]['final_populations_0'] = sum_pop_0 / num_reps
        mean_results[day_idx]['final_populations_p'] = sum_pop_p / num_reps
        mean_results[day_idx]['resource'] = sum_resource / num_reps

    return mean_results





def extractSubpopulationDensitiesMultiSpecies(results, day, Nmuts, Nins, plasmid, num_species):
    """
    Extracts the densities of subpopulations for all species on the specified day.

    Parameters:
    results (list): Simulation results containing populations and resources over time.
    day (int): Day for which to extract densities.
    Nmuts (int): Number of mutation levels.
    Nins (int): Number of transposition levels.
    plasmid (str): Type of plasmid ('free' or 'bearing').
    num_species (int): Number of species.

    Returns:
    numpy array: 3D array of population densities [species][SNP][IS].
    """
    if plasmid == 'free':
        populations = results[day - 1]['final_populations_0']
    elif plasmid == 'bearing':
        populations = results[day - 1]['final_populations_p']
    else:
        raise ValueError("Plasmid type must be 'free' or 'bearing'.")

    if populations.shape != (num_species, Nmuts, Nins):
        raise ValueError(f"Expected populations of size {(num_species, Nmuts, Nins)}, but got {populations.shape}.")
    return populations

def plotSubpopulationsGridCombinedPie(populations_free, populations_bearing, day, Nmuts, Nins, species_colors, species_labels):
    """
    Plots a grid of subpopulations, with pie charts at each point showing populations of plasmid-free
    and plasmid-bearing strains for all species, and includes a legend.

    Parameters:
    populations_free (numpy array): 3D array of plasmid-free population sizes [species][SNP][IS].
    populations_bearing (numpy array): 3D array of plasmid-bearing population sizes [species][SNP][IS].
    day (int): Current day for the title.
    Nmuts (int): Number of mutation levels.
    Nins (int): Number of transposition levels.
    species_colors (list): List of colors for each species.
    species_labels (list): List of labels for each species.
    """
    fig, ax = plt.subplots(figsize=(8, 8))

    # Compute maximum population across all grid points
    total_pops = populations_free.sum(axis=0) + populations_bearing.sum(axis=0)
    max_pop = np.max(total_pops) if np.max(total_pops) > 0 else 1  # Avoid division by zero

    # Iterate over grid points
    for i in range(Nmuts):  # SNP mutations (rows)
        for j in range(Nins):  # IS transpositions (columns)
            # Compute total population at this grid point
            total_pop = populations_free[:, i, j].sum() + populations_bearing[:, i, j].sum()
            if total_pop > 0:
                # Calculate the radius to scale area proportionally to population size
                radius = np.sqrt(total_pop / max_pop) * 0.4  # Adjust scaling factor for display

                start_angle = 0
                for species_idx, species_color in enumerate(species_colors):
                    # Free population slice
                    species_pop_free = populations_free[species_idx, i, j]
                    angle_free = (species_pop_free / total_pop) * 360 if total_pop > 0 else 0
                    wedge_free = Wedge((j, i), radius, start_angle, start_angle + angle_free,
                                       facecolor=species_color, edgecolor='black', alpha=0.4, hatch='//')
                    ax.add_patch(wedge_free)
                    start_angle += angle_free

                    # Bearing population slice
                    species_pop_bearing = populations_bearing[species_idx, i, j]
                    angle_bearing = (species_pop_bearing / total_pop) * 360 if total_pop > 0 else 0
                    wedge_bearing = Wedge((j, i), radius, start_angle, start_angle + angle_bearing,
                                          facecolor=species_color, edgecolor='black', alpha=0.8)
                    ax.add_patch(wedge_bearing)
                    start_angle += angle_bearing

    # Add legend
    legend_handles = []
    for species_idx, species_label in enumerate(species_labels):
        # Create patches for legend
        legend_handles.append(
            plt.Line2D(
                [0], [0],
                color=species_colors[species_idx],
                lw=10,
                label=f"{species_label}"
            )
        )

    ax.legend(
        handles=legend_handles,
        loc='upper right',
        bbox_to_anchor=(1.25, 1),  # Position legend outside the top right of the plot
        fontsize=12,
        title="Strains",
        title_fontsize=14
    )

    # Set labels, title, and grid layout
    ax.set_xlabel('# accumulated IS1', fontsize=14)
    ax.set_ylabel('# accumulated SNP', fontsize=14)
    ax.set_title(f'Day {day}', fontsize=16)
    ax.set_xlim(-0.5, Nins - 0.5)
    ax.set_ylim(-0.5, Nmuts - 0.5)
    ax.set_xticks(np.arange(Nins))
    ax.set_yticks(np.arange(Nmuts))
    ax.grid(True)

    # Set axis to square for equal scaling
    ax.set_aspect('equal', 'box')

    plt.tight_layout()



def plotSubpopulationsGridDaysCombinedPie(results, days, Nmuts, Nins, initial_populations_0, initial_populations_p, num_species, species_colors, species_labels, outPath=''):
    """
    Plots subpopulation grids for all species combined across multiple days.

    Parameters:
    results (list): Simulation results containing populations and resources over time.
    days (list): Days to plot.
    Nmuts (int): Number of mutation levels.
    Nins (int): Number of transposition levels.
    initial_populations_0 (numpy array): Initial plasmid-free populations [species][SNP][IS].
    initial_populations_p (numpy array): Initial plasmid-bearing populations [species][SNP][IS].
    num_species (int): Number of species.
    species_colors (list): List of colors for each species.
    outPath (str): Output path to save plots.
    """
    for day in days:
        if day == 0:
            populations_day_free = initial_populations_0
            populations_day_bearing = initial_populations_p
        else:
            populations_day_free = np.array([
                extractSubpopulationDensitiesMultiSpecies(results, day=day, Nmuts=Nmuts, Nins=Nins, plasmid='free', num_species=num_species)[species_idx]
                for species_idx in range(num_species)
            ])
            populations_day_bearing = np.array([
                extractSubpopulationDensitiesMultiSpecies(results, day=day, Nmuts=Nmuts, Nins=Nins, plasmid='bearing', num_species=num_species)[species_idx]
                for species_idx in range(num_species)
            ])

        # Create the plot for all species combined
        plotSubpopulationsGridCombinedPie(populations_day_free, populations_day_bearing, day, Nmuts, Nins, species_colors, species_labels)

        if outPath != '':
            # Define the file path with day appended
            filename = f"{outPath}_day{day}.png"

            # Save the plot to the specified path
            plt.savefig(filename, format='png')
            print("Exporting %s" % filename)
            plt.close()  # Close the plot to free memory after saving


def plotStrainDensitiesOverTime(
    results, num_days, num_species, species_colors, species_labels,
    initial_populations_0, initial_populations_p, title='', outPath=''
):
    """
    Plots the density of each strain (y-axis) against time (x-axis) for plasmid-free and plasmid-bearing subpopulations.

    Parameters:
    results (list): Simulation results containing populations over time.
    num_days (int): Total number of days in the simulation.
    num_species (int): Number of species in the simulation.
    species_colors (list): List of colors for each species.
    species_labels (list): List of labels for each species.
    initial_populations_0 (numpy array): Initial plasmid-free populations [species][SNP][IS].
    initial_populations_p (numpy array): Initial plasmid-bearing populations [species][SNP][IS].
    outPath (str): Path to save the plot (optional).
    """
    days = np.arange(0, num_days + 1)

    # Initialize figure
    fig, ax = plt.subplots(figsize=(8, 4))

    # Iterate through species
    for species_idx in range(num_species):
        # Collect time-series data for plasmid-free and plasmid-bearing populations
        densities_free = []
        densities_bearing = []

        for day in days:
            if day == 0:  # Initial populations
                densities_free.append(initial_populations_0[species_idx].sum())
                densities_bearing.append(initial_populations_p[species_idx].sum())
            else:  # Populations at subsequent days
                densities_free.append(results[day - 1]['final_populations_0'][species_idx].sum())
                densities_bearing.append(results[day - 1]['final_populations_p'][species_idx].sum())

        # Plot densities for plasmid-free (dotted line) and plasmid-bearing (solid line)
        ax.plot(
            days, densities_free, linestyle='dotted', color=species_colors[species_idx],
            label=f"{species_labels[species_idx]} (plasmid-free)"
        )
        ax.plot(
            days, densities_bearing, linestyle='solid', color=species_colors[species_idx],
            label=f"{species_labels[species_idx]} (plasmid-bearing)"
        )

    # Add labels, title, legend, and grid
    ax.set_xlabel('Time (days)', fontsize=16)
    ax.set_ylabel('Density', fontsize=16)
    ax.set_title(title, fontsize=16)
    ax.legend(fontsize=12, loc='upper right', bbox_to_anchor=(1.4, 1), title="Strains", title_fontsize=14)
    ax.set_yscale('log')
    ax.grid(False)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)

    # Adjust layout
    plt.tight_layout()

    # Save or show plot
    if outPath:
        filename = f"{outPath}_strain_densities_over_time.png"
        plt.savefig(filename, format='png', dpi=300)
        print(f"Exporting {filename}")
        plt.close()
    else:
        plt.show()


def calculateFrequenciesPerStrain(
    results, num_days, num_species, initial_populations_0, initial_populations_p
):
    """
    Calculates the frequencies per strain per day, normalized by the total population on each day.

    Parameters:
    results (list): Simulation results containing populations over time.
    num_days (int): Total number of days in the simulation.
    num_species (int): Number of species in the simulation.
    initial_populations_0 (numpy array): Initial plasmid-free populations [species][SNP][IS].
    initial_populations_p (numpy array): Initial plasmid-bearing populations [species][SNP][IS].

    Returns:
    dict: A dictionary containing frequencies for plasmid-free and plasmid-bearing strains per species.
    """
    # Create arrays to store the frequencies of plasmid-free and plasmid-bearing populations
    frequencies_free = np.zeros((num_species, num_days + 1))
    frequencies_bearing = np.zeros((num_species, num_days + 1))

    # Calculate total populations per day
    total_populations = []
    for day in range(num_days + 1):
        if day == 0:  # Initial populations (day 0)
            total_population = initial_populations_0.sum() + initial_populations_p.sum()
        else:  # Subsequent days
            total_population = (
                results[day - 1]['final_populations_0'].sum() +
                results[day - 1]['final_populations_p'].sum()
            )
        total_populations.append(total_population)

    total_populations = np.array(total_populations)

    # Iterate over species to calculate frequencies for plasmid-free and plasmid-bearing
    for species_idx in range(num_species):
        for day in range(num_days + 1):
            if day == 0:  # Initial populations (day 0)
                free_pop = initial_populations_0[species_idx].sum()
                bearing_pop = initial_populations_p[species_idx].sum()
            else:  # Subsequent days
                free_pop = results[day - 1]['final_populations_0'][species_idx].sum()
                bearing_pop = results[day - 1]['final_populations_p'][species_idx].sum()

            total_population = total_populations[day]
            # Calculate frequencies for plasmid-free and plasmid-bearing
            frequencies_free[species_idx, day] = free_pop / total_population if total_population > 0 else 0
            frequencies_bearing[species_idx, day] = bearing_pop / total_population if total_population > 0 else 0

    return frequencies_free, frequencies_bearing

def plotStrainFrequenciesStackedArea(
    frequencies_free, frequencies_bearing, num_days, num_species,
    species_colors, species_labels, title='', outPath=''
):
    """
    Plots the frequency of each strain (y-axis) against time (x-axis) using a single stacked area chart,
    combining plasmid-free and plasmid-bearing subpopulations.

    Parameters:
    frequencies_free (numpy array): Array of frequencies for plasmid-free populations [species][days].
    frequencies_bearing (numpy array): Array of frequencies for plasmid-bearing populations [species][days].
    num_days (int): Total number of days in the simulation.
    num_species (int): Number of species in the simulation.
    species_colors (list): List of colors for each species.
    species_labels (list): List of labels for each species.
    outPath (str): Path to save the plot (optional).
    """
    days = np.arange(0, num_days + 1)  # Includes day 0

    # Initialize figure
    fig, ax = plt.subplots(figsize=(8, 4))

    # Storage for cumulative frequencies for stacking
    cumulative_frequencies = np.zeros(len(days), dtype=float)

    for species_idx in range(num_species):
        # Add plasmid-free subpopulation
        prev_cumulative = cumulative_frequencies.copy()
        cumulative_frequencies += frequencies_free[species_idx, :]

        ax.fill_between(
            days, 100 * prev_cumulative, 100 * cumulative_frequencies,
            color=species_colors[species_idx], alpha=0.4, hatch='//'
        )

        # Add plasmid-bearing subpopulation
        prev_cumulative = cumulative_frequencies.copy()
        cumulative_frequencies += frequencies_bearing[species_idx, :]

        ax.fill_between(
            days, 100 * prev_cumulative, 100 * cumulative_frequencies,
            color=species_colors[species_idx], alpha=0.8,
            label=f"{species_labels[species_idx]} "
        )

    # Add labels, title, legend, and grid
    ax.set_xlabel('Time (days)', fontsize=16)
    ax.set_ylabel('Frequency (%)', fontsize=16)
    ax.set_title(title, fontsize=16)
    ax.legend(fontsize=12, loc='upper right', bbox_to_anchor=(1.2, 1), title="Strains", title_fontsize=14)
    ax.set_ylim(0, 100)
    ax.set_xlim(0, num_days)
    ax.grid(False)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)

    # Adjust layout
    plt.tight_layout()

    # Save or show plot
    if outPath:
        filename = f"{outPath}_strain_frequencies_stacked_area.png"
        plt.savefig(filename, format='png', dpi=300)
        print(f"Exporting {filename}")
        plt.close()
    else:
        plt.show()



def plotSubpopulationDynamics(
    results, days, Nmuts, Nins, initial_populations_0, initial_populations_p, species_colors, species_labels, outPath=''
):
    """
    Plots the densities of subpopulations with varying IS and SNP levels over time.
    Line width represents IS or SNP levels, and colors represent strains.

    Parameters:
    results (list of dict): Simulation results with population densities.
    days (list of int): List of days to include in the plot.
    Nmuts (int): Number of mutation levels (SNPs).
    Nins (int): Number of transposition levels (ISs).
    initial_populations_0 (np.array): Initial population of plasmid-free cells [species][SNP][IS].
    initial_populations_p (np.array): Initial population of plasmid-bearing cells [species][SNP][IS].
    species_colors (list): List of colors for each species.
    species_labels (list): Labels for each species.
    outPath (str): Path to save the plot (optional).
    """
    time_points = [0]  # Day 0

    # Initialize data storage for each species
    subpop_densities_B0 = [[[np.sum(initial_populations_0[species, i, j]) for j in range(Nins)] for i in range(Nmuts)] for species in range(len(species_colors))]
    subpop_densities_Bp = [[[np.sum(initial_populations_p[species, i, j]) for j in range(Nins)] for i in range(Nmuts)] for species in range(len(species_colors))]

    # Extract subpopulation data from results
    max_day = len(results)
    for day in days:
        if day < max_day:
            time_points.append(day)
            for species_idx in range(len(species_colors)):
                day_densities_B0 = results[day - 1]['final_populations_0'][species_idx]
                day_densities_Bp = results[day - 1]['final_populations_p'][species_idx]
                for i in range(Nmuts):
                    for j in range(Nins):
                        subpop_densities_B0[species_idx][i][j] = np.append(
                            subpop_densities_B0[species_idx][i][j], np.sum(day_densities_B0[i, j])
                        )
                        subpop_densities_Bp[species_idx][i][j] = np.append(
                            subpop_densities_Bp[species_idx][i][j], np.sum(day_densities_Bp[i, j])
                        )

    # Plot subpopulation dynamics
    fig, ax = plt.subplots(figsize=(8, 4))
    for species_idx, species_color in enumerate(species_colors):
        for i in range(Nmuts):
            for j in range(Nins):
                line_width = 0.5 + i * 2 + j * 2  # Line width increases with SNP and IS levels

                # Plot plasmid-free (dotted lines)
                ax.plot(
                    time_points,
                    subpop_densities_B0[species_idx][i][j],
                    linestyle='--',
                    linewidth=line_width,
                    color=species_color,
                    #label=f"{species_labels[species_idx]} (SNP {i}, IS {j})" if time_points.index(0) == 0 else None,
                    alpha=0.8,
                )

                # Plot plasmid-bearing (solid lines)
                ax.plot(
                    time_points,
                    subpop_densities_Bp[species_idx][i][j],
                    linestyle='-',
                    linewidth=line_width,
                    color=species_color,
                    #label=f"{species_labels[species_idx]} (SNP {i}, IS {j})" if time_points.index(0) == 0 else None,
                    alpha=0.8,
                )

    # Add labels and legend
    ax.set_xlabel('Time (days)', fontsize=14)
    ax.set_ylabel('Density', fontsize=14)
    #ax.set_yscale('log')
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    #ax.set_ylim([1e1, 1e9])
    #ax.legend(fontsize=10, loc='upper left', title='Subpopulations', title_fontsize=12)

    # Save or show the plot
    if outPath:
        filename = f"{outPath}_subpopulation_dynamics.png"
        plt.savefig(filename, format='png', dpi=300)
        print(f"Exporting {filename}")
        plt.close()
    else:
        plt.show()



def plotMutationAndISTranspositionHistogramsFromResults(
    results, Nmuts, Nins, num_species, species_colors, species_labels, outPath=''
):
    """
    Plots stacked histograms showing the distribution of cell populations across mutation levels and IS transpositions
    for multiple species from the simulation `results`.

    Parameters:
    results (list of dict): Simulation results containing final populations over days.
    Nmuts (int): Number of mutation levels.
    Nins (int): Number of transposition levels.
    num_species (int): Number of species.
    species_colors (list): List of colors for each species.
    species_labels (list): Labels for each species.
    outPath (str): Path to save the plot (optional).
    """
    # Initialize arrays to store cumulative populations across days
    total_populations_mutations = np.zeros((num_species, Nmuts))
    total_populations_transpositions = np.zeros((num_species, Nins))

    # Iterate over results to accumulate final populations
    for day_result in [results[-1]]: #final day
        final_populations_0 = day_result["final_populations_0"]
        final_populations_p = day_result["final_populations_p"]

        # Sum populations by SNP (axis 2) and IS (axis 1)
        total_populations_mutations += final_populations_0.sum(axis=2) + final_populations_p.sum(axis=2)
        total_populations_transpositions += final_populations_0.sum(axis=1) + final_populations_p.sum(axis=1)

    # Total population across all SNP and IS levels
    total_population = total_populations_mutations.sum()

    # Variables to store cumulative frequencies excluding Level 0
    snp_cumulative_freq = 0
    is_cumulative_freq = 0

    # Debugging print statements to verify calculations
    print("Mutation Totals by SNP Level and Frequencies (stacked by Species):")
    for level in range(Nmuts):
        level_total = total_populations_mutations[:, level].sum()
        frequency = (level_total / total_population) * 100 if total_population > 0 else 0
        print(f"  SNP Level {level}: Total Population = {level_total}, Frequency = {frequency:.2f}%")
        if level > 0:  # Exclude Level 0
          snp_cumulative_freq += frequency

    print("\nTransposition Totals by IS Level and Frequencies (stacked by Species):")
    for level in range(Nins):
        level_total = total_populations_transpositions[:, level].sum()
        frequency = (level_total / total_population) * 100 if total_population > 0 else 0
        print(f"  IS Level {level}: Total Population = {level_total}, Frequency = {frequency:.2f}%")
        if level > 0:  # Exclude Level 0
          is_cumulative_freq += frequency

    # Display cumulative frequencies
    print("\nCumulative Frequencies Excluding Level 0:")
    print(f"  SNP Levels > 0: {snp_cumulative_freq:.2f}%")
    print(f"  IS Levels > 0: {is_cumulative_freq:.2f}%")

    # Plot histograms
    fig, axes = plt.subplots(1, 2, figsize=(8, 4), sharey=True)

    # Stacked SNP-level histogram
    bottom = np.zeros(Nmuts)
    for species_idx in range(num_species):
        axes[0].bar(
            range(Nmuts),
            total_populations_mutations[species_idx],
            bottom=bottom,
            color=species_colors[species_idx],
            edgecolor='black',
            label=species_labels[species_idx],
            alpha=0.8,
        )
        bottom += total_populations_mutations[species_idx]
    axes[0].set_xlabel('# SNP', fontsize=14)
    axes[0].set_ylabel('Total Cell Population', fontsize=14)
    axes[0].set_xticks(range(Nmuts))

    #axes[0].set_yscale('log')
    axes[0].set_ylim([1, 1.1*total_population])
    #axes[0].legend(fontsize=10, title='Species', loc='upper left')

    # Stacked IS-level histogram
    bottom = np.zeros(Nins)
    for species_idx in range(num_species):
        axes[1].bar(
            range(Nins),
            total_populations_transpositions[species_idx],
            bottom=bottom,
            color=species_colors[species_idx],
            edgecolor='black',
            label=species_labels[species_idx],
            alpha=0.8,
        )
        bottom += total_populations_transpositions[species_idx]
    axes[1].set_xlabel('# IS1', fontsize=14)
    #axes[1].set_yscale('log')
    axes[1].set_ylim([1, 1.1*total_population])
    axes[1].legend(fontsize=10, title='Species', loc='upper right')
    axes[1].set_xticks(range(Nins))


    # Adjust layout
    plt.tight_layout()

    if outPath:
        filename = f"{outPath}_mutation_transposition_histograms_stacked.pdf"
        plt.savefig(filename, format='pdf', bbox_inches='tight')
        print(f"Exporting {filename}")
        plt.close()
    else:
        plt.show()

    return snp_cumulative_freq, is_cumulative_freq

def plotCumulativeFrequenciesSNPandIS(snp_cumulative_freq, is_cumulative_freq, title='', outPath=''):
    """
    Plots cumulative frequencies as two bars (SNP and IS).

    Parameters:
    snp_cumulative_freq (float): Cumulative frequency for SNP levels > 0.
    is_cumulative_freq (float): Cumulative frequency for IS levels > 0.
    title (str): Plot title.
    outPath (str): Path to save the plot (optional).
    """
    # Data for the two bars
    categories = ['SNP', 'IS1']
    frequencies = [snp_cumulative_freq, is_cumulative_freq]

    # Create the plot
    plt.figure(figsize=(3, 6))
    bars = plt.bar(categories, frequencies, color=['skyblue', 'salmon'], edgecolor='black')

    # Add values on top of the bars
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2, height + 1, f'{height:.2f}%',
                 ha='center', va='bottom', fontsize=12)

    # Add labels and formatting
    plt.ylabel('Cumulative Frequency (%)', fontsize=14)
    plt.title(title, fontsize=16)
    plt.ylim(0, 105)  # The maximum value is 100% since it's a cumulative frequency

    # Save or show the plot
    plt.tight_layout()
    if outPath:
        filename = f"{outPath}_cumulative_frequencies.pdf"
        plt.savefig(filename, format='pdf', bbox_inches='tight')
        print(f"Exporting {filename}")
        plt.close()
    else:
        plt.show()



def computeAccumulatedMutations(results, Nmuts, Nins, num_species):
    """
    Computes the total number of accumulated mutations (SNPs + IS1s) for each species.

    Parameters:
    results (dict): Simulation results for a single case.
    Nmuts (int): Number of mutation levels.
    Nins (int): Number of transposition levels.
    num_species (int): Number of species/strains.

    Returns:
    np.array: A 1D array (species) containing the total number of accumulated mutations per species.
    """
    accumulated_mutations = np.zeros(num_species)

    # Get the final day's population data
    final_day = results[-1]  # Get the last recorded time point
    final_populations_0 = final_day["final_populations_0"]  # Plasmid-free
    final_populations_p = final_day["final_populations_p"]  # Plasmid-bearing

    # Compute total accumulated mutations per species
    for species_idx in range(num_species):
        # SNP Mutations: Weight by level
        for level in range(1, Nmuts):  # Skip level 0 (non-mutant)
            accumulated_mutations[species_idx] += level * np.sum(final_populations_0[species_idx, level, :])
            accumulated_mutations[species_idx] += level * np.sum(final_populations_p[species_idx, level, :])

        # IS Transpositions: Weight by level
        for level in range(1, Nins):  # Skip level 0 (non-mutant)
            accumulated_mutations[species_idx] += level * np.sum(final_populations_0[species_idx, :, level])
            accumulated_mutations[species_idx] += level * np.sum(final_populations_p[species_idx, :, level])

    return accumulated_mutations


def computeMutationAndISTranspositionStatistics(results, Nmuts, Nins, num_species):
    """
    Computes statistics on the distribution of cell populations across mutation (SNP) levels
    and IS transpositions for multiple species from the simulation `results`.

    Parameters:
    results (list of dict): Simulation results containing final populations over days.
    Nmuts (int): Number of mutation levels.
    Nins (int): Number of transposition levels.
    num_species (int): Number of species.

    Returns:
    dict: A dictionary containing the total populations, frequencies, and cumulative frequencies.
    """
    # Initialize arrays to store cumulative populations across days
    total_populations_mutations = np.zeros((num_species, Nmuts))
    total_populations_transpositions = np.zeros((num_species, Nins))

    # Extract final day's population data
    final_day = results[-1]  # Get results from the last simulation day
    final_populations_0 = final_day["final_populations_0"]
    final_populations_p = final_day["final_populations_p"]

    # Sum populations by SNP (axis 2) and IS (axis 1)
    total_populations_mutations += final_populations_0.sum(axis=2) + final_populations_p.sum(axis=2)
    total_populations_transpositions += final_populations_0.sum(axis=1) + final_populations_p.sum(axis=1)

    # Total population across all SNP and IS levels
    total_population = total_populations_mutations.sum()

    # Compute frequencies and cumulative frequencies
    snp_frequencies = []
    is_frequencies = []
    snp_cumulative_freq = 0
    is_cumulative_freq = 0

    for level in range(Nmuts):
        level_total = total_populations_mutations[:, level].sum()
        frequency = (level_total / total_population) * 100 if total_population > 0 else 0
        snp_frequencies.append({"SNP_Level": level, "Total_Population": level_total, "Frequency": frequency})
        if level > 0:  # Exclude Level 0
            snp_cumulative_freq += frequency

    for level in range(Nins):
        level_total = total_populations_transpositions[:, level].sum()
        frequency = (level_total / total_population) * 100 if total_population > 0 else 0
        is_frequencies.append({"IS_Level": level, "Total_Population": level_total, "Frequency": frequency})
        if level > 0:  # Exclude Level 0
            is_cumulative_freq += frequency

    # Compute the total number of Non-Mutant cells (cells that have neither SNPs nor IS transpositions)
    #total_non_mutants = np.sum(final_populations_0[:, 0, 0] + final_populations_p[:, 0, 0])

    total_mutant_cells = computeAccumulatedMutations(results, Nmuts, Nins, num_species)


    # Return the statistics in a dictionary
    return {
        "SNP_Frequencies": snp_frequencies,
        "IS_Frequencies": is_frequencies,
        "Total_Population": total_population,
        "Total_Mutant_Cells": total_mutant_cells,
        "SNP_Cumulative_Frequency": snp_cumulative_freq,
        "IS_Cumulative_Frequency": is_cumulative_freq
    }



def computeAccumulatedMutations(results, Nmuts, Nins, num_species):
    """
    Computes the total number of accumulated mutations (SNPs + IS1s) for each species.

    Parameters:
    results (dict): Simulation results for a single case.
    Nmuts (int): Number of mutation levels.
    Nins (int): Number of transposition levels.
    num_species (int): Number of species/strains.

    Returns:
    np.array: A 1D array (species) containing the total number of accumulated mutations per species.
    """
    accumulated_mutations = np.zeros(num_species)

    # Get the final day's population data
    final_day = results[-1]  # Get the last recorded time point
    final_populations_0 = final_day["final_populations_0"]  # Plasmid-free
    final_populations_p = final_day["final_populations_p"]  # Plasmid-bearing

    # Compute total accumulated mutations per species
    for species_idx in range(num_species):
        # SNP Mutations: Weight by level
        for level in range(1, Nmuts):  # Skip level 0 (non-mutant)
            accumulated_mutations[species_idx] += level * np.sum(final_populations_0[species_idx, level, :])
            accumulated_mutations[species_idx] += level * np.sum(final_populations_p[species_idx, level, :])

        # IS Transpositions: Weight by level
        for level in range(1, Nins):  # Skip level 0 (non-mutant)
            accumulated_mutations[species_idx] += level * np.sum(final_populations_0[species_idx, :, level])
            accumulated_mutations[species_idx] += level * np.sum(final_populations_p[species_idx, :, level])

    return accumulated_mutations



def printMutationAndISTranspositionStatistics(stats):
    """
    Prints mutation (SNP) and IS transposition statistics in a human-readable format.

    Parameters:
    stats (dict): A dictionary containing SNP and IS frequencies, total population, and cumulative frequencies.
    """

    # Convert NumPy array to scalar (sum across species if needed)
    total_mutant_cells = stats['Total_Mutant_Cells']
    if isinstance(total_mutant_cells, np.ndarray):
        total_mutant_cells = total_mutant_cells.sum()  # Sum across all species

    # Print total population and non-mutant cell count
    print("\n=== Population Overview ===")
    print(f"  Total Population: {stats['Total_Population']:.4g}")
    print(f"  Total Accumulated Mutations: {total_mutant_cells:.4g}")  # Now correctly formatted

    print("\n=== Mutation spectrum ===\n")

    for entry in stats["SNP_Frequencies"]:
        print(f"  SNP Level {entry['SNP_Level']}: Total Population = {entry['Total_Population']:.4g}, Frequency = {entry['Frequency']:.2f}%")

    for entry in stats["IS_Frequencies"]:
        print(f"  IS Level {entry['IS_Level']}: Total Population = {entry['Total_Population']:.4g}, Frequency = {entry['Frequency']:.2f}%")

    print("\nFrequency of Mutants:")
    print(f"  SNP Levels > 0: {stats['SNP_Cumulative_Frequency']:.2f}%")
    print(f"  IS Levels > 0: {stats['IS_Cumulative_Frequency']:.2f}%")



def computeCumulativeStrainDensities(results, num_species, num_days):
    """
    Computes cumulative strain densities over time by summing all IS and SNP subpopulations
    for both plasmid-free and plasmid-bearing populations.

    Parameters:
    results (list of dict): Simulation results containing final populations over multiple days.
    num_species (int): Number of species in the simulation.
    num_days (int): Total number of days in the simulation.

    Returns:
    tuple:
        - cumulative_densities_free (numpy array): Cumulative densities of plasmid-free strains [species][days].
        - cumulative_densities_bearing (numpy array): Cumulative densities of plasmid-bearing strains [species][days].
    """
    # Initialize arrays to store cumulative strain densities for each species over time
    cumulative_densities_free = np.zeros((num_species, num_days + 1))  # +1 to include day 0
    cumulative_densities_bearing = np.zeros((num_species, num_days + 1))

    # Extract initial populations (Day 0)
    initial_day = results[0]  # First entry in results
    for species_idx in range(num_species):
        cumulative_densities_free[species_idx, 0] = np.sum(initial_day["final_populations_0"][species_idx])
        cumulative_densities_bearing[species_idx, 0] = np.sum(initial_day["final_populations_p"][species_idx])

    # Loop over simulation days and sum up densities
    for day in range(1, num_days + 1):
        day_results = results[day - 1]  # Note: results are 0-indexed
        for species_idx in range(num_species):
            cumulative_densities_free[species_idx, day] = np.sum(day_results["final_populations_0"][species_idx])
            cumulative_densities_bearing[species_idx, day] = np.sum(day_results["final_populations_p"][species_idx])

    return cumulative_densities_free, cumulative_densities_bearing

def printCumulativeStrainDensities(cumulative_densities_free, cumulative_densities_bearing, num_species, species_labels):
    """
    Prints the final cumulative strain densities in a human-readable format,
    marking the strain with the highest total density with an asterisk (*).

    Parameters:
    cumulative_densities_free (numpy array): Cumulative densities of plasmid-free strains [species][days].
    cumulative_densities_bearing (numpy array): Cumulative densities of plasmid-bearing strains [species][days].
    num_species (int): Number of species in the simulation.
    species_labels (list): List of species names corresponding to the indices.

    Returns:
    None
    """
    print("\n=== Final Cumulative Strain Densities ===\n")

    # Compute total densities for all species
    total_densities = cumulative_densities_free[:, -1] + cumulative_densities_bearing[:, -1]

    # Determine the index of the species with the highest total density
    max_idx = np.argmax(total_densities)

    # Print headers
    print(f"{'Species':<10}{'Plasmid-Free':>20}{'Plasmid-Bearing':>20}{'Total':>20}")
    print("-" * 70)

    for species_idx in range(num_species):
        free_density = cumulative_densities_free[species_idx, -1]  # Last time point
        bearing_density = cumulative_densities_bearing[species_idx, -1]  # Last time point
        total_density = free_density + bearing_density

        # Add an asterisk (*) to the species with the highest total density
        label = f"{species_labels[species_idx]} *" if species_idx == max_idx else species_labels[species_idx]

        print(f"{label:<10}{free_density:>20.2e}{bearing_density:>20.2e}{total_density:>20.2e}")



def computeStrainDiversity(results):
    """
    Computes Shannon entropy and Simpson's diversity index at the **strain level**,
    summing over all SNP and IS subpopulations.

    Parameters:
    results (list of dict): Simulation results with population densities.

    Returns:
    dict: Dictionary containing Shannon entropy, Simpson's index, and total number of surviving strains.
    """
    final_day = results[-1]  # Extract final day's results

    # Sum over SNP and IS levels -> get total population per strain
    total_pop_per_strain_0 = np.sum(final_day['final_populations_0'], axis=(1, 2))  # Sum SNP and IS levels
    total_pop_per_strain_p = np.sum(final_day['final_populations_p'], axis=(1, 2))  # Sum SNP and IS levels

    # Combine plasmid-free and plasmid-bearing populations
    total_pop_per_strain = total_pop_per_strain_0 + total_pop_per_strain_p

    # Remove zero-density strains
    nonzero_strain_pops = total_pop_per_strain[total_pop_per_strain > 0]

    # Compute diversity only if there are surviving strains
    total_population = np.sum(nonzero_strain_pops)
    if total_population == 0:
        return {"Shannon_entropy": 0, "Simpson_index": 0, "Total_strains": 0}

    # Compute relative abundances per strain
    proportions = nonzero_strain_pops / total_population

    # Shannon entropy (H')
    shannon_entropy = entropy(proportions, base=2)

    # Simpson's diversity index (D)
    simpson_index = 1 - np.sum(proportions**2)

    return {
        "Shannon_entropy": shannon_entropy,
        "Simpson_index": simpson_index,
        "Total_strains": len(nonzero_strain_pops)
    }


def printStrainDiversity(diversity_stats):
    """
    Prints strain diversity metrics in a human-readable format.

    Parameters:
    diversity_stats (dict): Dictionary containing Shannon entropy, Simpson's index, and strain counts.

    Returns:
    None
    """
    print("\n=== Strain Diversity Metrics ===\n")
    print(f"{'Metric':<20}{'Value'}")
    print("-" * 35)
    print(f"{'Shannon Entropy':<20}{diversity_stats['Shannon_entropy']:.4f}")
    print(f"{'Simpson Index':<20}{diversity_stats['Simpson_index']:.4f}")
    print(f"{'Total Surviving Strains':<20}{diversity_stats['Total_strains']}")



def computePlasmidFraction(results):
    """
    Computes cumulative densities of plasmid-bearing and plasmid-free cells,
    returning their total counts and plasmid fraction.

    Parameters:
    results (list of dict): Simulation results containing final populations over multiple days.

    Returns:
    dict: Dictionary containing final cumulative densities and plasmid fraction.
    """
    final_day = results[-1]  # Extract final day's results

    # Sum over all species, SNP, and IS levels
    total_population_free = np.sum(final_day['final_populations_0'])  # Plasmid-free
    total_population_bearing = np.sum(final_day['final_populations_p'])  # Plasmid-bearing

    total_population = total_population_free + total_population_bearing

    # Compute plasmid fraction
    plasmid_fraction = total_population_bearing / total_population if total_population > 0 else 0

    return {
        "Plasmid-Free Density": total_population_free,
        "Plasmid-Bearing Density": total_population_bearing,
        "Total Density": total_population,
        "Plasmid Fraction": plasmid_fraction
    }


def printPlasmidFraction(plasmid_stats):
    """
    Prints the plasmid-bearing and plasmid-free cell densities,
    along with the plasmid fraction, in a human-readable format.

    Parameters:
    plasmid_stats (dict): Dictionary containing cumulative densities and plasmid fraction.

    Returns:
    None
    """
    print("\n=== Final Plasmid Fraction Summary ===\n")
    print(f"{'Category':<25}{'Density'}")
    print("-" * 45)
    print(f"{'Plasmid-Free Cells':<25}{plasmid_stats['Plasmid-Free Density']:.2e}")
    print(f"{'Plasmid-Bearing Cells':<25}{plasmid_stats['Plasmid-Bearing Density']:.2e}")
    print(f"{'Total Cells':<25}{plasmid_stats['Total Density']:.2e}")
    print(f"{'Plasmid Fraction':<25}{plasmid_stats['Plasmid Fraction']:.4f}")


def plotWeightedCumulativeMutationsByTypePerStrain(
    results, Nmuts, Nins, num_species, species_colors, species_labels, ymax=5e7, title='', outPath=''
):
    """
    Plots cumulative weighted SNP and IS1 mutations (levels multiplied by populations) as two columns,
    stacking contributions from different strains using distinct colors.

    Parameters:
    results (list of dict): Simulation results containing final populations over days.
    Nmuts (int): Number of mutation levels.
    Nins (int): Number of transposition levels.
    num_species (int): Number of species/strains.
    species_colors (list): List of colors for each strain.
    species_labels (list): Labels for each species/strain.
    ymax (float): Maximum y-axis value.
    title (str): Plot title.
    outPath (str): Path to save the plot (optional).
    """
    # Initialize cumulative weighted mutations for each species
    weighted_snps = np.zeros(num_species)
    weighted_is = np.zeros(num_species)

    # Iterate over results to calculate weighted mutations (final day only)
    for species_idx in range(num_species):
        final_populations_0 = results[-1]["final_populations_0"][species_idx]
        final_populations_p = results[-1]["final_populations_p"][species_idx]

        # Sum weighted SNPs: Levels * Populations (skip Level 0)
        for level in range(1, Nmuts):
            weighted_snps[species_idx] += level * np.sum(final_populations_0[level, :])
            weighted_snps[species_idx] += level * np.sum(final_populations_p[level, :])

        # Sum weighted IS1s: Levels * Populations (skip Level 0)
        for level in range(1, Nins):
            weighted_is[species_idx] += level * np.sum(final_populations_0[:, level])
            weighted_is[species_idx] += level * np.sum(final_populations_p[:, level])

    # Plot two stacked bars (SNPs and IS1s) with contributions from all strains
    x_positions = [0, 1]  # Positions for SNP and IS1 bars
    bottom_snp = np.zeros(1)  # Start bottom for SNP stacking
    bottom_is = np.zeros(1)   # Start bottom for IS1 stacking

    plt.figure(figsize=(4, 6))

    for species_idx in range(num_species):
        # Stack contributions for SNP
        plt.bar(
            x_positions[0], weighted_snps[species_idx],
            bottom=bottom_snp,
            color=species_colors[species_idx],
            edgecolor='black'
        )
        bottom_snp += weighted_snps[species_idx]

        # Stack contributions for IS1
        plt.bar(
            x_positions[1], weighted_is[species_idx],
            bottom=bottom_is,
            color=species_colors[species_idx],
            edgecolor='black'
        )
        bottom_is += weighted_is[species_idx]

    # Add labels and formatting
    plt.xticks(x_positions, ['SNP', 'IS1'], fontsize=14)
    plt.ylabel('# Mutant Cells', fontsize=16)
    plt.title(title, fontsize=16)
    plt.ylim([0, ymax])
    plt.yticks(fontsize=14)

    # Create a legend with one entry per species
    handles = [plt.Rectangle((0, 0), 1, 1, color=species_colors[i], edgecolor='black') for i in range(num_species)]
    plt.legend(handles, species_labels, fontsize=12, loc='upper left', bbox_to_anchor=(1.05, 1), title="Strains")

    # Save or show plot
    plt.tight_layout()
    if outPath:
        filename = f"{outPath}_cumulative_mutations_stacked.pdf"
        plt.savefig(filename, format='pdf', bbox_inches='tight')
        print(f"Exporting {filename}")
        plt.close()
    else:
        plt.show()


# Function to plot the HGT network
def plot_HGT_network(matrix, strainIDs, species_colors, outPath='', node_alpha=0.4):
    G = nx.DiGraph()  # Directed graph for HGT

    # Add nodes
    for i, strain in enumerate(strainIDs):
        G.add_node(strain)

    # Add edges with weights
    edges = []
    weights = []
    for i, from_strain in enumerate(strainIDs):
        for j, to_strain in enumerate(strainIDs):
            if matrix[i, j] > 0:  # Only add edges if there is transmission
                G.add_edge(from_strain, to_strain, weight=matrix[i, j])
                edges.append((from_strain, to_strain))
                weights.append(matrix[i, j])  # Store edge weight

    # Define positions for better layout
    pos = nx.spring_layout(G, seed=42)  # Spring layout for clear visualization

    fig, ax = plt.subplots(figsize=(6, 6))

    # Ensure edges appear on top of nodes
    ax.set_axisbelow(False)

    # Draw edges with adjusted arrow positioning
    edge_widths = np.array(weights) * 15000  # Scale edge thickness
    nx.draw_networkx_edges(
        G, pos, edgelist=edges, width=edge_widths, edge_color='gray', arrows=True, alpha=0.8,
        connectionstyle="arc3,rad=0.1",  # Slight curvature to avoid node overlap
        arrowsize=15  # Increase arrowhead size for better visibility
    )

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_size=1500, node_color=species_colors, alpha=node_alpha)

    # Draw labels
    nx.draw_networkx_labels(G, pos, labels={strain: strain for strain in strainIDs}, font_size=12, font_weight='bold', font_color='black')

    if outPath:
        filename = f"{outPath}_network.pdf"
        plt.savefig(filename, format='pdf', bbox_inches='tight')
        print(f"Exporting {filename}")
        plt.close()
    else:
        plt.show()



def plotResults(matrix, populations_0, populations_p, results, num_mutationsSNP, num_mutationsIS, num_days, num_species, species_colors, strainIDs, pathFIGURES, expe_path, antibiotic_concentration):
    """
    Generates multiple plots based on the simulation results.

    Parameters:
    results (list of dict): Simulation results containing final populations over multiple days.
    num_mutationsSNP (int): Number of SNP mutation levels.
    num_mutationsIS (int): Number of IS transposition levels.
    num_days (int): Total number of days in the simulation.
    num_species (int): Number of species in the simulation.
    species_colors (list): List of colors for each species.
    strainIDs (list): Labels for each species.
    pathFIGURES (str): Base directory where figures will be saved.
    expe_path (str): Subdirectory for experiment-specific figures.
    antibiotic_concentration (float or str): Antibiotic concentration label for title purposes.

    Returns:
    None
    """
    
    # 0. Plot network
    plot_HGT_network(matrix, strainIDs, species_colors, outPath=expe_path, node_alpha=0.45)


    # 1. Plot subpopulation grid as pie charts
    #plotSubpopulationsGridDaysCombinedPie(
    #    results=results,
    #    days=list(range(num_days + 1)),  # Assuming we want to plot all days
    #    Nmuts=num_mutationsSNP,
    #    Nins=num_mutationsIS,
    #    initial_populations_0=populations_0,
    #    initial_populations_p=populations_p,
    #    num_species=num_species,
    #    species_colors=species_colors,
    #    species_labels=strainIDs,
    #    outPath=outPath
    #)

    # 2. Plot strain densities over time
    plotStrainDensitiesOverTime(
        results=results,
        num_days=num_days,
        num_species=num_species,
        species_colors=species_colors,
        species_labels=strainIDs,
        initial_populations_0=populations_0,
        initial_populations_p=populations_p,
        title=f'A={antibiotic_concentration}',
        outPath=expe_path
    )

    # 3. Compute and plot SNP and IS transposition histograms
    snp_cumulative_freq, is_cumulative_freq = plotMutationAndISTranspositionHistogramsFromResults(
        results=results,
        Nmuts=num_mutationsSNP,
        Nins=num_mutationsIS,
        num_species=num_species,
        species_colors=species_colors,
        species_labels=strainIDs,
        outPath=expe_path
    )

    # 4. Calculate frequencies per day per strain
    frequencies_free, frequencies_bearing = calculateFrequenciesPerStrain(
        results=results,
        num_days=num_days,
        num_species=num_species,
        initial_populations_0=populations_0,
        initial_populations_p=populations_p
    )

    # 5. Plot stacked area chart of strain frequencies
    plotStrainFrequenciesStackedArea(
        frequencies_free, frequencies_bearing, num_days, num_species, species_colors, strainIDs, '', expe_path
    )

    # 6. Plot subpopulation dynamics
    plotSubpopulationDynamics(
        results, list(range(num_days + 1)), num_mutationsSNP, num_mutationsIS, populations_0, populations_p, species_colors, strainIDs, expe_path
    )

    # 7. Plot SNP and IS cumulative frequencies
    plotCumulativeFrequenciesSNPandIS(
        snp_cumulative_freq=snp_cumulative_freq,
        is_cumulative_freq=is_cumulative_freq,
        title="",
        outPath=expe_path
    )


    # 8. Cumulative Mutations per strain
    ymax=3.5e8
    plotWeightedCumulativeMutationsByTypePerStrain(
        results=results,
        Nmuts=num_mutationsSNP,
        Nins=num_mutationsIS,
        num_species=num_species,
        species_colors=species_colors,
        species_labels=strainIDs,
        ymax=ymax,
        title='',
        outPath=expe_path
    )


# Function to plot scatter plots
def plot_parameter_scatter(strainIDs, strain_parameters, strain_color_map):
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    for species_idx, conditions in strain_parameters.items():
        strainID = strainIDs[species_idx]  # Convert index to actual strain name
        color = strain_color_map.get(strainID, "black")  # Default to black if strain not in dict

        # Iterate over plasmid-free ("0") and plasmid-bearing ("p") conditions
        for condition, levels in conditions.items():
            strain_label = strainID if condition == "0" else f"{strainID}p"  # Plasmid-bearing gets 'K253p'

            for level_idx, params in enumerate(levels):
                alpha = 0.3 + 0.7 * (1 - level_idx / len(levels))  # Higher level -> lower transparency

                # Extract parameters
                death_rate = params['death_rate']
                half_saturation_antibiotic = params['half_saturation_antibiotic']
                birth_rate = params['birth_rate']
                consumption_rate = params['consumption_rate']
                mutation_rate = params['mutation_rate']
                transposition_rate = params['transposition_rate']

                # Subplot 1: Death rate vs. Half-saturation antibiotic
                # Use the same circular marker for both conditions
                marker = "o"
                is_plasmid_bearing = (condition == "p")

                # Choose fill style
                facecolor = color if is_plasmid_bearing else 'none'
                edgecolor = color

                # Subplot 1
                axes[0].scatter(
                    death_rate, half_saturation_antibiotic,
                    facecolors=facecolor,
                    edgecolors=edgecolor,
                    alpha=alpha,
                    marker=marker,
                    label=strain_label if level_idx == 0 else "",
                    s=100
                )

                # Subplot 2
                axes[1].scatter(
                    birth_rate, consumption_rate,
                    facecolors=facecolor,
                    edgecolors=edgecolor,
                    alpha=alpha,
                    marker=marker,
                    s=100
                )

                # Subplot 3
                axes[2].scatter(
                    mutation_rate, transposition_rate,
                    facecolors=facecolor,
                    edgecolors=edgecolor,
                    alpha=alpha,
                    marker=marker,
                    s=100
                )


    # Set labels and titles
    axes[0].set_xlabel("Death Rate", fontsize=20)
    axes[0].set_ylabel("Half-Saturation Antibiotic", fontsize=20)
    axes[0].set_title("Death Rate vs. Half-Saturation Antibiotic", fontsize=20)
    axes[0].tick_params(axis='both', labelsize=18)

    axes[1].set_xlabel("Birth Rate", fontsize=20)
    axes[1].set_ylabel("Consumption Rate", fontsize=20)
    axes[1].set_title("Birth Rate vs. Consumption Rate", fontsize=20)
    axes[1].tick_params(axis='both', labelsize=18)

    axes[2].set_xlabel("SNP Rate ", fontsize=20)
    axes[2].set_ylabel("IS Rate", fontsize=20)
    axes[2].set_title("Mutation Rate vs. Transposition Rate", fontsize=20)
    axes[2].tick_params(axis='both', labelsize=18)

    # Adjust ticks
    for ax in axes:
        ax.grid(True)

    # Legend positioned outside to the right
    axes[0].legend(fontsize=12, loc='upper left', bbox_to_anchor=(1, 1))

    plt.tight_layout()
    plt.show()


def plotTimeDependentStrainDensities(results_list, num_species, species_colors, species_labels, simulation_time, title='', outPath=''):
    """
    Plots the time-dependent densities of plasmid-free and plasmid-bearing populations
    for each strain in isolation, using horizontal subplots.

    Parameters:
    results_list (list): List of simulation results (one dictionary per strain).
    num_species (int): Number of species in the simulation.
    species_colors (list): List of colors for each species.
    species_labels (list): List of labels for each species.
    title (str): Plot title.
    outPath (str): Path to save the plot (optional).
    """
    ymax=1.5e9
    fig, axes = plt.subplots(nrows=1, ncols=num_species, figsize=(3 * num_species, 4), sharey=True)
    #fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 4), sharey=True)

    for species_idx in range(num_species):
        ax = axes[species_idx]


        # Ensure results_list has enough data
        if species_idx >= len(results_list):
            print(f"Skipping species {species_idx} (no data)")
            continue

        species_results = results_list[species_idx]  # Dictionary with keys ['strain', 'results']
        if 'results' not in species_results or len(species_results['results']) == 0:
            print(f"Skipping species {species_idx} due to missing data.")
            continue

        # Extract time series data for the FIRST simulated day
        first_day_results = species_results['results'][0]
        #time_points = first_day_results['time_points']


        population_values = first_day_results['population_values']
        time_points = np.linspace(0, simulation_time, len(population_values))


        if len(time_points) == 0 or len(population_values) == 0:
            print(f"Skipping species {species_idx} due to missing time-dependent data.")
            continue

        # Extract plasmid-free and plasmid-bearing densities over time
        densities_free = [pop[0][species_idx].sum() for pop in population_values]  # Summing over mutations
        densities_bearing = [pop[1][species_idx].sum() for pop in population_values]  # Summing over mutations

        # Plot time-dependent densities
        ax.plot(time_points, densities_free, linestyle='dotted', color=species_colors[species_idx], lw=3, label="Plasmid-Free")
        ax.plot(time_points, densities_bearing, linestyle='solid', color=species_colors[species_idx], lw=3, label="Plasmid-Bearing")

        # Set labels and formatting
        ax.set_xlabel('Time', fontsize=16)
        ax.set_title(species_labels[species_idx], fontsize=16)
        #ax.set_yscale('log')
        ax.grid(False)
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        ax.set_ylim([0, ymax])
        ax.set_xticks([0, 6, 12, 18, 24])

        if species_idx == 0:  # Only show ylabel on first subplot
            ax.set_ylabel('Density', fontsize=16)

    # Add a single legend for all plots
    handles, labels = ax.get_legend_handles_labels()
    #fig.legend(handles, labels, fontsize=14, loc='upper right', bbox_to_anchor=(1.17, 0.92))

    black_free = mlines.Line2D([], [], color='black', linestyle='dotted', lw=3, label="Plasmid-Free")
    black_bearing = mlines.Line2D([], [], color='black', linestyle='solid', lw=3, label="Plasmid-Bearing")

    # Add the legend with black lines
    fig.legend(handles=[black_free, black_bearing], fontsize=14, loc='upper right', bbox_to_anchor=(1.17, 0.92))
    # Adjust layout
    plt.tight_layout()

    # Save or show the plot
    if outPath:
        filename = f"{outPath}_time_dependent_strain_densities.png"
        plt.savefig(filename, format='png', dpi=300)
        print(f"Exporting {filename}")
        plt.close()
    else:
        plt.show()


def binary_filename_to_matrix(binary_string, shape):
    """
    Converts a single-line binary string back into a numpy matrix.

    Parameters:
    binary_string (str): Binary matrix representation.
    shape (tuple): Original shape of the matrix (rows, columns).

    Returns:
    numpy array: Reconstructed binary matrix.
    """
    binary_list = np.array(list(map(int, binary_string))).reshape(shape)
    return binary_list



def matrix_to_binary_filename(matrix):
    """
    Converts a numerical matrix to a single-line binary string where values >0 are '1' and 0 remains '0'.
    The output can be used as part of a filename.

    Parameters:
    matrix (numpy array): Input matrix.

    Returns:
    str: Single-line binary string representation of the matrix.
    """
    binary_matrix = (matrix > 0).astype(int)  # Convert to 0s and 1s
    binary_string = "".join(map(str, binary_matrix.flatten()))  # Flatten and convert to a single string
    return binary_string


def save_simulation_results(plasmid_matrix, num_reps, initial_values_by_strain, num_mutationsSNP, num_mutationsIS,
                            strainIDs, results, this_diversity_metrics, this_stats,
                            this_cumulative_densities_free, this_cumulative_densities_bearing,
                            plasmid_fraction_stats, save_path=""):
    """
    Saves simulation results and metadata into a pickle file.

    Parameters:
    plasmid_matrix (numpy array): The plasmid transmission matrix.
    num_reps: (int): Number of simulation repetitions.
    initial_values_by_strain (dict): Initial strain population values.
    num_mutationsSNP (int): Number of SNP mutation levels.
    num_mutationsIS (int): Number of IS transposition levels.
    strainIDs (list): List of strain identifiers.
    results (dict): Simulation results.
    this_diversity_metrics (dict): Strain diversity metrics.
    this_stats (dict): Mutation and IS transposition statistics.
    this_cumulative_densities_free (numpy array): Cumulative plasmid-free densities.
    this_cumulative_densities_bearing (numpy array): Cumulative plasmid-bearing densities.
    plasmid_fraction_stats (dict): Plasmid fraction statistics.
    save_path (str): Directory to save the file (default: current directory).

    Returns:
    str: The filename of the saved pickle file.
    """
    # Convert matrix to a binary string for filename
    binary_filename = matrix_to_binary_filename(plasmid_matrix)
    filename = f"sim_{binary_filename}.pkl"
    full_path = os.path.join(save_path, filename)

    # Save results as a dictionary
    data_to_save = {
        "plasmid_matrix": plasmid_matrix,
        "num_reps": num_reps,
        "initial_values_by_strain": initial_values_by_strain,
        "num_mutationsSNP": num_mutationsSNP,
        "num_mutationsIS": num_mutationsIS,
        "strainIDs": strainIDs,
        "results": results,
        "diversity_metrics": this_diversity_metrics,
        "mutation_statistics": this_stats,
        "cumulative_densities_free": this_cumulative_densities_free,
        "cumulative_densities_bearing": this_cumulative_densities_bearing,
        "plasmid_fraction_stats": plasmid_fraction_stats
    }

    # Save as pickle
    with open(full_path, 'wb') as file:
        pickle.dump(data_to_save, file)

    print(f"\nSimulation results saved to: {full_path}")
    return full_path

def load_simulation_results(file_path):
    """
    Loads simulation results and metadata from a pickle file.

    Parameters:
    file_path (str): Path to the pickle file.

    Returns:
    dict: Dictionary containing all stored simulation data.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    with open(file_path, 'rb') as file:
        data = pickle.load(file)

    print(f"Simulation results loaded from: {file_path}")
    return data


def plot_dose_response_curve(dose_response_results, strainIDs, species_colors, outPath=''):
    """
    Plots dose-response curves for each strain, showing the final cell density at different drug concentrations.

    Parameters:
    dose_response_results (dict): Dictionary with dose-response data per strain.
    strainIDs (list): List of strain names.
    species_colors (list): List of colors for plotting.
    """
    plt.figure(figsize=(6, 4))

    for idx, strain in enumerate(strainIDs):
        data = dose_response_results[strain]
        concentrations = [d["antibiotic_concentration"] for d in data]
        final_populations_0 = [d["final_population_0"] for d in data]
        final_populations_p = [d["final_population_p"] for d in data]

        # Plot plasmid-free as dotted line
        plt.plot(concentrations, final_populations_0, '--', color=species_colors[idx], label=f"{strain} (Plasmid-Free)")

        # Plot plasmid-bearing as solid line
        plt.plot(concentrations, final_populations_p, '-', color=species_colors[idx], label=f"{strain} (Plasmid-Bearing)")

    # Formatting
    plt.xlabel("Antibiotic Concentration", fontsize=16)
    plt.ylabel("Final Density (cells/mL)", fontsize=16)
    plt.xscale("log")  # Log scale for better visualization
    plt.ylim([1e5, 1.5e9])  # Adjust y-axis range if needed
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=12, loc='upper left', bbox_to_anchor=(1, 1.025))

    plt.grid(True)


    # Save or show the plot
    if outPath:
        filename = f"{outPath}doseResponse.pdf"
        plt.savefig(filename, format='pdf', bbox_inches='tight', pad_inches=0.2)
        print(f"Exporting {filename}")
        plt.show()
        plt.close()
    else:
        plt.show()


# Convert plasmid matrix to NetworkX graph
def matrix_to_graph(matrix):
    G = nx.from_numpy_array(np.array(matrix), create_using=nx.DiGraph)
    return G

# Compute connectivity metrics
def compute_network_measures(matrix):
    G = matrix_to_graph(matrix)

    density = nx.density(G)
    mean_degree = np.mean([d for n, d in G.degree()])

    # Get the largest strongly connected component
    giant_component_size = len(max(nx.strongly_connected_components(G), key=len))

    # Compute average path length only if the graph is strongly connected
    avg_path_length = nx.average_shortest_path_length(G) if nx.is_strongly_connected(G) else None

    # Clustering coefficient for undirected version of the graph
    clustering = nx.average_clustering(G.to_undirected())

    return {
        "Density": density,
        "Mean Degree": mean_degree,
        "Giant Component Size": giant_component_size,
        "Avg Path Length": avg_path_length,
        "Clustering Coefficient": clustering
    }


