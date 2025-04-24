'''
This file <NNES.py> 
Aims at Training DNN for NN-based estimator for learning uncertainty in CBF

Authored by R. Ma @ Harbin Institute of Technology
Updated June 2022

'''

# import modules 
from tkinter import Y
import torch
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from torchvision import datasets
from torchvision.transforms import ToTensor
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd   
from sklearn import preprocessing
import os
import time


# Get cpu or gpu device for training. --> CUDA
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Using {device} device")


# Define feature dataset for training
class FeatureDataset(Dataset):  # troch.utils.data --> generate torch tensors and have torch compatible dataset
    def __init__(self,file_name):
        
        # read csv file and load row data into variables
        '''
        We do not seperate the dataset into train and test set. Instead, we consider complete file as training set. 
        This is because at training the torch model can separate data accordingly for given 'train:test' ratio.
        '''
        file_out = pd.read_csv(file_name) 
        xx = file_out.iloc[:,0:2].values   
        uu = file_out.iloc[:, 2].values   
        '''
        scale our dataset using StandardScalar object to standardized the features
        '''
        sc = preprocessing.StandardScaler()  # scikit-learn --> data splitting and pre-processing
        x_train = sc.fit_transform(xx)      
        u_train = uu 
        
        # converting to troch tensors
        '''
        cast data into float32 may be essential as default dtype returned after scaling is of double and not accepted by torch model.
        '''
        self.X_train = torch.tensor(x_train, dtype=torch.float32)
        self.U_train = torch.tensor(u_train, dtype=torch.float32)
        print (self.X_train.shape,self.U_train.shape)
        
        
    def __len__(self):
        return len(self.U_train)
    
    def __getitem__(self,idx):
        return self.X_train[idx], self.U_train[idx]
        

# Define FFNN model
class FFNNController(torch.nn.Module):
    def __init__(self):
        super(FFNNController,self).__init__()
        self.flatten = torch.nn.Flatten() 
        self.linear_tanh_stack = torch.nn.Sequential(
            torch.nn.Linear(2,32,bias = False),  
            torch.nn.Tanh(),
            torch.nn.Linear(32,16,bias = False),
            torch.nn.Tanh(),
            torch.nn.Linear(16,1,bias = False)   
        ) 
        
    def forward(self, x):
       # x = self.flatten(x) # check whether flatten is necessary?
        logits = self.linear_tanh_stack(x)
        return logits
    
model = FFNNController().to(device)
print(model)

batch_size = 32

# Load the training data   
#feature_data = FeatureDataset('F:\\Huawei\\Code23\\exp_dataCBF01.csv')  # Load the stored data file
feature_data = FeatureDataset('F:\\Huawei\\Code23\\exp_dataCBF02.csv')  # Load the stored data file
len_feature = len(feature_data) # Examine the length of feature_data (the number of state-control pairs)
print(len_feature)

train_size = int(0.8 * len_feature)
test_size = len_feature - train_size
train_dataset, test_dataset = torch.utils.data.random_split(feature_data,[train_size,test_size])
train_dataloader = DataLoader(train_dataset,batch_size,shuffle=True,num_workers=1)
test_dataloader = DataLoader(test_dataset)
len_train = len(train_dataset)  # Calculate the length of train data --> len(train_dataloader)*batch_size=len(train_dataset)
len_test = len(test_dataset) # Calculate the length of test data
print(len_train,len_test)

for xx, uu in test_dataloader:
    print(f"Shape of xx [N, C, H, W]: {xx.shape} {xx.dtype}")
    print(f"Shape of uu: {uu.shape} {uu.dtype}")
    break

# Test the training data has been loaded
hh = 1 # flag with no meanings
print(hh)


# Define loss function and optimizer
loss_ffnn = torch.nn.L1Loss()
optimizer = torch.optim.Adam(model.parameters(),lr=1e-3)

'''
In a single training loop, the model makes predicitons on the training dataset (fed to it in batches)
and propagates the prediction error to adjust the model parameters
'''
def train(dataloader,model,loss_ffnn,optimizer):
    size = len(dataloader.dataset)
    model.train()
    for batch, (xx,uu) in enumerate (dataloader):
        xx, uu = xx.to(device), uu.to(device)
        
        # Compute prediction error
        pred = model(xx)
        loss = loss_ffnn(pred,uu)
        
        # Backpropagitation
        optimizer.zero_grad() # set gradients to 0 for each batch
        loss.backward()   # backward propagation
        optimizer.step() # update weights
        
        if batch % 100 == 0:
            loss, current =loss.item(), batch * len(xx)
            print(f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")

'''
Check the model's performance against the test dataset to ensure it is learning
'''

def test(dataloader,model,loss_ffnn):
    size = len(dataloader.dataset)
    num_batches = len(dataloader)
    model.eval()
    test_loss, correct =0, 0
    with torch.no_grad():
        for xx , uu in dataloader:
            xx , uu = xx.to(device), uu.to(device)
            pred = model(xx)
            test_loss += loss_ffnn(pred,uu).item()
            correct += (pred.argmax(1)==uu).type(torch.float).sum().item()
    test_loss /= num_batches
    correct /= size
    print(f"Test Error: \n Accuracy: {(100*correct):>0.1f}%, Avg loss: {test_loss:>8f} \n")

'''
The training process is conducted over several epoches (iterations), during each of which the model learns parameters to make 
better predictions;
Print the model's accuracy and loss at each epoch
'''  


epochs = 6
def main():
    for t in range (epochs):
        print(f"Epoch {t+1}\n-------------------------------")
        train(train_dataloader,model,loss_ffnn,optimizer)
        test(test_dataloader,model,loss_ffnn)
    print("Done")
    
    # print the state_dict for models --> https://blog.csdn.net/weixin_40522801/article/details/106563354
    print("Model's state_dict:")
    for param_tensor in model.state_dict():
        print(param_tensor,"\t",model.state_dict()[param_tensor].size())
        
    # print the state_dict for optimizer
    print("Optimizer's state_dict:")
    for var_name in optimizer.state_dict():
        print(var_name,"\t",optimizer.state_dict()[var_name])
    
    #torch.save(model.state_dict(),'F:\\Huawei\\Code23\\save_weightE0.pt') # Save weights eposidic learning 1
    torch.save(model.state_dict(),'F:\\Huawei\\Code23\\save_weightE1.pt') # Save weights eposidic learning 2
        

if __name__=='__main__':
    main()