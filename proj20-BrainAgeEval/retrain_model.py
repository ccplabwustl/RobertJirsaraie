#!/usr/bin/python
#################

import tensorflow as tf ; import keras
from keras import layers ; from keras import optimizers
from keras.preprocessing import image ; from keras import models
from keras.applications.resnet50 import preprocess_input, decode_predictions
from sklearn.model_selection import cross_val_predict, cross_val_score, KFold
from keras.preprocessing.image import ImageDataGenerator ; import math
from keras import regularizers ; from time import time
from tensorflow.python.keras.callbacks import CSVLogger
from sklearn.model_selection import train_test_split
import numpy as np ; import pandas as pd ; import cv2
import h5py ; from keras.models import load_model
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import OneHotEncoder
from statistics import median, median_grouped
import os ; from os import listdir, fsencode 
import matplotlib.pyplot as plt ; import sys
from sklearn.svm import LinearSVR, SVR

DataFreeze=sys.argv[1] ; ImgDir=sys.argv[2] ; ModelDir=sys.argv[3] ; SPLIT=sys.argv[4] ; SIZE=sys.argv[5]

######
### Aggregate Datasets so that Each File Path has an Associated Age
######

master = pd.read_csv('{}'.format(DataFreeze)) ; df = pd.read_csv('{}/Retrain.csv'.format(ImgDir))
if SIZE == "half":
	master = master[master.split==1]
	master = master.reset_index() ; del master['index']

df["age"] = ""
for index, row in master.iterrows():
	INDICES=df.index[df['sub'] == row['sub']].tolist()
	df.age[INDICES]=row['age']

master['Split'] = 0 ; INDEX=0
kf = KFold(n_splits=5,shuffle = True,random_state=42) ; kf.get_n_splits(master['sub']) 
for train_index, test_index in kf.split(master['sub']):
	INDEX=INDEX+1 ; master.loc[test_index,'Split']=INDEX

validate = df[df['sub'].isin(master[master.Split==int(SPLIT)].iloc[:,0])]
validate["age"]=pd.to_numeric(validate["age"]) ; validate = validate.reset_index()
train = df[df['sub'].isin(master[master.Split!=int(SPLIT)].iloc[:,0])]
train["age"]=pd.to_numeric(train["age"]) ; train = train.reset_index() 

######
### Read and Prepare Data For Deep Learning
######

datagenT=ImageDataGenerator(rescale=1. / 255, horizontal_flip= True, vertical_flip = True,featurewise_center = False, featurewise_std_normalization = False)
train_generator=datagenT.flow_from_dataframe(
dataframe=train,
directory=ImgDir,
x_col="path",
y_col="age",
has_ext=True,
batch_size=64,
seed=42,
shuffle=True,
class_mode="other")

datagenV=ImageDataGenerator(rescale=1. / 255,horizontal_flip= False, vertical_flip = False,featurewise_center = False, featurewise_std_normalization = False)
valid_generator=datagenV.flow_from_dataframe(
dataframe=validate,
directory=ImgDir,
x_col="path",
y_col="age",
has_ext=True,
batch_size=64,
seed=42,
shuffle=True,
class_mode="other")

######
### Start Retraining Model With 
######

model = load_model('{}/DBN_Pediatric.h5'.format(ModelDir)) ; model.summary()

layers_to_train = 4
for layer in model.layers[:(-1*(layers_to_train))]:
    layer.trainable = False

for layer in model.layers:
    print(layer.name,layer.trainable)

######
### Retrain Model with First Epoch Freezing the First Layer
######

model.compile(loss=keras.losses.mean_squared_error,
              optimizer=optimizers.Adam(lr =7e-5,decay = 3e-7),
              metrics=['mae'])

print("Now Training Epoch: 1")
STEP_SIZE_TRAIN=train_generator.n//train_generator.batch_size
STEP_SIZE_VALID=valid_generator.n//valid_generator.batch_size
model.fit_generator(generator=train_generator,
	steps_per_epoch=STEP_SIZE_TRAIN,
	validation_data=valid_generator,
	validation_steps=STEP_SIZE_VALID,
	epochs=1,workers = 8)
model.save('{}/DBN_tPediatric_{}_epoch-1_split-{}.h5'.format(ModelDir,SIZE,SPLIT))

######
### Unlock All Layers To Continue Retraining Until Convergance
######

for layer in model.layers[:]:
	layer.trainable = True

for ITERATION in range(1,29):
	print("Now Training Epoch: {}".format(ITERATION+1))
	STEP_SIZE_TRAIN=train_generator.n//train_generator.batch_size
	STEP_SIZE_VALID=valid_generator.n//valid_generator.batch_size
	model.fit_generator(generator=train_generator,
		steps_per_epoch=STEP_SIZE_TRAIN,
		validation_data=valid_generator,
		validation_steps=STEP_SIZE_VALID,
		epochs=1,workers = 8)
	model.save('{}/DBN_tPediatric_{}_epoch-{}_split-{}.h5'.format(ModelDir,SIZE,ITERATION+1,SPLIT))

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
