# -*- coding: utf-8 -*-

import cobra as cb
import pandas as pd
import numpy as np
from scipy.sparse.csgraph import connected_components
from cobra.util.array import create_stoichiometric_matrix

class eFVA:
    def __init__(self,model):
        self.model=model
   
        reactions=[reaction.id for reaction in model.reactions]    
        Df_system=pd.DataFrame(data=create_stoichiometric_matrix(model) ,index=range(len(model.metabolites)),columns=reactions)

        valori_nonzero=Df_system.astype(bool).sum(axis=1)
        valori_nonzero=valori_nonzero[valori_nonzero==2]

        Df_system=Df_system.loc[valori_nonzero.index]   
        Df_system=Df_system.loc[:,(Df_system!=0).any()] 
    
        lista_reazioni=list(Df_system.columns)

        Df2=pd.DataFrame(data=0,index=lista_reazioni,columns=lista_reazioni)
        for index, row in Df_system.iterrows(): 
            lista=row[row!=0].index
            Df2.loc[lista[0],lista[1]]=1

        self.Df_connected=pd.DataFrame(index=lista_reazioni,columns=["labels"])
         
        n_components, labels =connected_components(Df2.values, directed=False, return_labels=True)
     
 
        self.Df_connected["labels"]=labels
        self.Df_connected.reset_index(inplace=True)
        self.Df_connected=self.Df_connected.groupby("labels").agg(lambda x:list(x))
       
        self.reaction_computed=list(self.Df_connected["index"].apply(lambda x:x[0]).values)
        self.reaction_not_computed=[reaction for reaction in lista_reazioni if reaction not in self.reaction_computed]
        self.reaction_fva=self.reaction_computed.copy()
        self.reaction_fva.extend([reaction for reaction in reactions if reaction not in lista_reazioni])
 

        rows=pd.DataFrame(data=0,index=self.reaction_computed,columns=lista_reazioni)
      
        for el in self.reaction_computed:
            rows.loc[el,el]=1 

        Df_system=Df_system.append(rows,ignore_index=True)

        self._b_system=np.zeros(len(valori_nonzero.index))
        self._Df_system=np.linalg.pinv(Df_system.to_numpy(dtype=float))
        self._df=pd.DataFrame(index=Df_system.columns,columns=["minimum","maximum"])
 
    def fastFVA(
            self,
            processes=1,
            fraction_of_optimum=0,
            round_elementd=6,
            ):

        df2=cb.flux_analysis.flux_variability_analysis(self.model,
                                                       self.reaction_fva,
                                                       processes=processes,
                                                       fraction_of_optimum=fraction_of_optimum)
        
        b_system_max=np.append(self._b_system,df2.loc[self.reaction_computed,"maximum"].values,axis=None)
        b_system_min=np.append(self._b_system,df2.loc[self.reaction_computed,"minimum"].values,axis=None)

        sol_max=self._Df_system.dot(b_system_max)
        sol_min=self._Df_system.dot(b_system_min)

        for min_val,max_val,i in zip(sol_min,sol_max,range(len(sol_min))):
            if min_val>max_val:
                sol_min[i]=max_val
                sol_max[i]=min_val

        self._df["minimum"]=sol_min
        self._df["maximum"]=sol_max

        df2=df2.append(self._df.loc[self.reaction_not_computed,:]).round(round_elementd)

        return df2

    def fastBlocked(self,processes=1):
        
        reactions=[self.model.reactions.get_by_id(reaction) for reaction in self.reaction_fva]
        
        lista_blocked=cb.flux_analysis.find_blocked_reactions(self.model,
                                                      reactions,
                                                      processes=processes)
        lista_blocked_total=list()
  
        for index, reaction_computed in zip(self.Df_connected.index,self.reaction_computed):
            row=self.Df_connected.loc[index,"index"]
            if reaction_computed in lista_blocked:

                lista_blocked_total.extend(row)
     
        lista_blocked_total.extend([reaction for reaction in lista_blocked if reaction not in lista_blocked_total]) 
     
        return lista_blocked_total
