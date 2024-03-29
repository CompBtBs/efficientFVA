# -*- coding: utf-8 -*-
"""
Created on Sun May 15 22:25:01 2022

@author: Tyrion
"""
import cobra as cb
import pandas as pd
from class_fastfva import eFVA
from cobrapy_bigg_client import client

cb.Configuration.solver="gurobi"

#download specific model from BiGG database
model=client.download_model("e_coli_core,save=False,file_format="json")
#inizialize efficient FVA object and compute connected sets
fva_object=eFVA(model)

#compute FVA and save the results in pandas DataFrame
df=fva_object.fastFVA() 
