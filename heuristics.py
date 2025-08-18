import math
from re import S
import gurobipy as gb
import numpy as np
import data_inizialization as di
import networkx as nx
from data import data
from instance import OptimizationInstance



# HEURISTIC HR

def HR_changes(istance: OptimizationInstance):
    """
    Heuristic HR for the optimization instance.
    This heuristic is designed to provide a feasible solution quickly.
    """
    R = istance.R

    RHR = set(R)

    stop_condition = False

    HR_istance = 0
    
    return HR_istance











# HEURISTIC HRBC

