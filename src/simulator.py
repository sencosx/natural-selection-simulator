# Copyright (c) 2025 Francesco Giuseppe Gillio

'''
Computer-Aided Simulation Framework
to Investigate Natural Selection and Species Extinction Dynamics
'''

import argparse
import copy
import itertools
import math
import os
import queue
import random
import warnings

import matplotlib.pyplot as plt
import numpy
import scipy

from utils import *


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--selections', default=100, type=int,
                        help='number of selection simulations')
    parser.add_argument('--timeframe', default=100, type=int,
                        help='simulation timeframe')
    parser.add_argument('--rate', default=0.15, type=float,
                        help='reproduction rate')
    parser.add_argument('--alpha', default=0.25, type=float,
                        help='improvement probability parameter')
    parser.add_argument('--improve', default=0.5, type=float,
                        help='improvement rate')
    parser.add_argument('--lifetime', default=15, type=int,
                        help='average lifetime')
    parser.add_argument('--directory', default='results', type=str,
                        help='path to the output location for checkpoint storage')
    return parser.parse_args()


def simulation(homo,
               clock,
               rate,
               alpha,
               improve,
               lifetime,
               districts):
    # initialize the simulation time [years]
    time = 0
    # the colony register
    world = grid(districts)
    # the clan register
    nature = dict()
    # the species register
    history = dict()
    # the exrtinction register
    extinction = dict()
    # randomly choose the origin sites for each clan
    origins = random.sample(
        [colony for colony in world.keys()], len(homo.keys()))
    for species in ['sapiens', 'neanderthalensis', 'erectus']:
        history[species] = dict()
        history[species][time] = 0
    # for each clan
    for clan, colony in zip(homo.keys(), origins):
        (species, population) = homo[clan]
        # initialize the dict data type to archive the clan individuals
        individuals = {code: Individual(code,
                                        time,
                                        clan,
                                        colony,
                                        None,
                                        species,
                                        numpy.random.uniform(low=0,
                                                             high=lifetime,
                                                             size=None)) for code in range(population)}
        # insert the clan into the clan register
        nature[clan] = Clan(clan,
                            species,
                            individuals,
                            colony,
                            population,
                            0)
        # insert the colony into the colony register
        world[colony].colonization(clan, species)
        # insert the species into the species register
        history[species][time] += population
        extinction[species] = 0
    # initialize the FES
    FES = queue.PriorityQueue()
    # for each clan
    for clan in nature.keys():
        # the clan species
        species = nature[clan].species
        # the clan colony
        colony = nature[clan].colonies[-1]
        # randomly choose the next migration time
        m = numpy.random.uniform(low=0, high=5, size=None)
        # for each individual within the clan
        for individual in nature[clan].individuals.values():
            # insert the death of the individual into the FES
            FES.put((individual.death(), (individual.code,
                                          individual.clan,
                                          individual.species,
                                          individual.colony,
                                          'death')))
            # insert the next birth by the individual into the FES
            if individual.death() - time > 0:
                x = individual.begetting(rate)
                if time + x < individual.death():
                    FES.put((time + x, (individual.code,
                                        individual.clan,
                                        individual.species,
                                        individual.colony,
                                        'birth')))
        # insert the next clan migration into the FES
        FES.put((time + m, (00,
                            clan,
                            species,
                            colony,
                            'departure')))
    # event loop
    while bool(FES.queue):
        # get the future event from the FES
        (time, (code, clan, species, colony, event)) = FES.get()
        print(end='\r|%-100s|' % ('-' * (100 * (int(time)) // clock)))
        # event loop termination criteria
        if time > clock:
            break
        if event == 'birth':
            birth(time,
                  code,
                  clan,
                  species,
                  colony,
                  rate,
                  alpha,
                  improve,
                  nature,
                  history,
                  world,
                  FES)
        if event == 'death':
            FES = death(time,
                        code,
                        clan,
                        species,
                        colony,
                        rate,
                        alpha,
                        improve,
                        nature,
                        history,
                        world,
                        FES)
        if event == 'departure':
            FES = departure(time,
                            code,
                            clan,
                            species,
                            colony,
                            rate,
                            alpha,
                            improve,
                            nature,
                            history,
                            world,
                            FES)
        if event == 'arrival':
            FES = arrival(time,
                          code,
                          clan,
                          species,
                          colony,
                          rate,
                          alpha,
                          improve,
                          nature,
                          history,
                          world,
                          FES)
        if event == 'clash':
            FES = clash(time,
                        code,
                        clan,
                        species,
                        colony,
                        rate,
                        alpha,
                        improve,
                        nature,
                        history,
                        world,
                        FES)
        # update the history on the counts
        history[species][time] = sum(nature[clan].history(
        ) for clan in nature.keys() if not nature[clan].species != species)
        for species in history.keys():
            if history[species][list(history[species])[-1]] < 1:
                extinction[species] = 1
        if sum([value for value in extinction.values()]) == len(homo.keys()):
            break
    return world, history, nature, extinction


def output(x, cofidence):
    if numpy.array(x).size < 2:
        return {'mean': numpy.mean(x),
                'deviation': 0,
                'accuracy': 0,
                'interval': 0}
    degrees = 1
    mean = numpy.mean(x)
    deviation = math.sqrt(numpy.var(x, ddof=degrees))
    sem = deviation / math.sqrt(numpy.array(x).size)
    interval = scipy.stats.t.interval(cofidence, degrees, mean, sem)
    error = ((max(interval) - min(interval)) / 2) / mean
    accuracy = 1 - error if 1 - error > 0 else 0
    return {'mean': mean,
            'deviation': deviation,
            'accuracy': accuracy,
            'interval': interval}


plt.rcParams['font.family'] = 'serif'

plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

font = {'fontsize': '7.5', 'color': 'black'}

homo = {'alpha': ('sapiens', 25),
        'beta': ('sapiens', 25),
        'gamma': ('sapiens', 25),
        'delta': ('sapiens', 25),
        'epsilon': ('sapiens', 25),

        'zeta': ('neanderthalensis', 25),
        'eta': ('neanderthalensis', 25),
        'theta': ('neanderthalensis', 25),
        'iota': ('neanderthalensis', 25),
        'kappa': ('neanderthalensis', 25),

        'lambda': ('erectus', 25),
        'omicron': ('erectus', 25),
        'rho': ('erectus', 25),
        'sigma': ('erectus', 25),
        'tau': ('erectus', 25)}

results = {'sapiens': list(),
           'neanderthalensis': list(),
           'erectus': list()}


def main():
    args = parse_args()
    warnings.filterwarnings("ignore")

    if not os.path.exists(args.directory):
        os.mkdir(args.directory)

    selections = args.selections
    confidence = 0.95
    districts = 25

    iter = 0
    rate = args.rate
    alpha = args.alpha
    clock = args.timeframe
    improve = args.improve
    lifetime = args.lifetime

    for selection in range(selections):
        iter += 1
        world, history, nature, extinction = simulation(homo,
                                                        clock,
                                                        rate,
                                                        alpha,
                                                        improve,
                                                        lifetime,
                                                        districts)
        colors = {'sapiens': 'darkblue',
                  'neanderthalensis': 'darkorange',
                  'erectus': 'darkgreen'}

        for species in history.keys():
            plt.plot(history[species].keys(),
                     history[species].values(),
                     label=f'homo {species}',
                     color=colors[species],
                     linestyle='--',
                     marker='.',
                     markersize=1.5,
                     linewidth=0.5)
            results[species].append(extinction[species])
            print(
                f'\nextinction rate at iteration {iter}: {species} - {output(results[species], confidence)["mean"]:.3f} with accuracy {output(results[species], confidence)["accuracy"]:.3f}')
        plt.tick_params(axis='both',
                        labelsize='7.5',
                        color='black')
        plt.legend(prop={'size': '10'})
        plt.savefig(f'{args.directory}/history.png', dpi=300)
        plt.close()

        for species in history.keys():
            X = [colony.x for colony in world.values()]
            Y = [colony.y for colony in world.values()]
            plt.scatter(X, Y, s=25, c='black', alpha=0.25)
            for colony in world.keys():
                for neighbor in world[colony].neighbors:
                    plt.plot([world[colony].x,
                              world[neighbor].x],
                             [world[colony].y,
                              world[neighbor].y],
                             color='black',
                             linewidth=0.5,
                             alpha=0.25)
            X, Y, colors = [], [], {'sapiens': 'darkblue',
                                    'neanderthalensis': 'darkorange',
                                    'erectus': 'darkgreen'}
            for clan in [clan for clan in nature.keys()
                         if nature[clan].species in [species]]:
                for colony in nature[clan].colonies:
                    X.append(world[colony].x)
                    Y.append(world[colony].y)
            plt.scatter(X, Y,
                        s=75,
                        c=colors[species],
                        label=f'homo {species}')
            plt.legend(prop={'size': '10'})
            plt.tick_params(axis='both',
                            labelsize='7.5',
                            color='black')
            plt.savefig(f'{args.directory}/{species}.png', dpi=300)
            plt.close()


if __name__ == '__main__':
    main()