'''
This file <WSC.py> ---> Weights Store to .csv file
Aims at transforming the .pt file to .csv file 

Authored by R. Ma @ Harbin Institute of Technology
Updated April 2022

'''
import torch
import csv
import codecs
from multiprocessing.connection import deliver_challenge
from NNCT import model # From <NNCT.py> import the trained FFNN model
import numpy as np
import pandas as pd



def main1():
    print(model.state_dict().keys())
    
    print("Model's state_dict:")
    for param_tensor in model.state_dict():
        print(param_tensor,"\t",model.state_dict()[param_tensor].size())
    
    print(list(model.linear_tanh_stack.parameters()))
    
  #  loaded_paras = torch.load('F:\\Huawei\\Code23\\save_weight0.pt')
  #  print(list(loaded_paras.linear_tanh_stack.parameters()))
 #   for param_tensor in loaded_paras:
 #     print(param_tensor, "\t", loaded_paras[param_tensor].size())
    
    
if __name__=='__main__':
    main1()
 
    