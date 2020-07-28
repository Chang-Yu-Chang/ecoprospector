#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 27 2019
@author: changyuchang
"""
import numpy as np
import scipy as sp


# Design migration_factor (a sequence of binary factors)
def no_migration(community_function):
    """
    No migration
    """
    # Number of wells
    n_wells = len(community_function)

    # No migration
    migration_factor = np.zeros(n_wells)

    return migration_factor

def parent_migration(community_function):
    """
    Parent migration, migrate into all wells
    """
    # Number of wells
    n_wells = len(community_function)

    # All migration
    migration_factor = np.ones(n_wells)

    #dont migrate into winner
    winner_index = np.where(community_function >= np.max(community_function))[0][::-1] # Reverse the list so the higher
    migration_factor[winner_index] = 0
    return migration_factor

def directed_selection_migrate(community_function):
    """
    Sample new communities from species pool, coalesce the migrant communities to the species pools
    """
    # Number of wells
    n_wells = len(community_function)

    # Compute the cutoff based on the number of wells
    cut_off_percent = (np.sqrt(n_wells))/n_wells

    # Sort the community function in this transfer
    sorted_community_function = np.sort(community_function)

    # Community function value cutoff for selecting communities
    cut_off = sorted_community_function[int(np.round(len(community_function)*(1-cut_off_percent)))]

    # Winner wells
    winner_index = np.where(community_function >= cut_off)[0][::-1]

    # Migration factor. A list of whether to migrate the community or not
    migration_factor = np.ones(n_wells) # Migrate all the wells except for the new wells that contain the winner replicate
    migration_factor[range(len(winner_index))] = 0 # Don't migrate to the winner wells

    return migration_factor

def migrate_half(community_function):
    # Number of wells
    n_wells = len(community_function)

    # Migration
    migration_factor = [1, 0] * int(n_wells/2)

    return migration_factor


def migrate_random(community_function):
    # Number of wells
    n_wells = len(community_function)

    # Migration
    migration_factor = np.random.binomial(1, 0.5, size = n_wells)

    return migration_factor
