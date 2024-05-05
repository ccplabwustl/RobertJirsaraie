#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#######################

import pandas as pd, sys, joblib, shap, numpy as np, seaborn as sns, matplotlib.pyplot as plt, os
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV, cross_val_score, KFold
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from sklearn.ensemble import GradientBoostingRegressor,AdaBoostRegressor
from sklearn.inspection import permutation_importance
from sklearn.preprocessing import MinMaxScaler
from sklearn.neighbors import KernelDensity
from sklearn.utils import resample
from scipy.special import softmax
import scipy.stats as stats
from datetime import date

######
### Define the Inner & Outer Loops For The Cross Validation
######

TODAY=sys.argv[2]

def OuterCVLoop(df, anticlust, colname):
	test=df.loc[df['anticlusters'] == anticlust] 
	train=df.loc[df['anticlusters'] != anticlust] 
	train_labels = train[colname].values ; test_labels = test[colname].values
	train_features = train.iloc[:,22:train.shape[1]] ; test_features = test.iloc[:,22:test.shape[1]]
	if MODE == 'FEAT':
		if SUBSET == 'Multimod':
			feature_subset = featmap['features'].tolist()
		else:		
			feature_subset = featmap[featmap['feature_type'] == SUBSET]['features'].tolist()
	elif MODE == 'NET':
		if SUBSET == 'Global':
			feature_subset = featmap[featmap['network_label'].notna()]['features'].tolist()
		else:
			feature_subset = featmap[featmap['network_label'] == SUBSET]['features'].tolist()
	test_features = test_features[feature_subset]
	train_features = train_features[feature_subset]
	return [train_labels, train_features, test_labels, test_features]

def InnerAdaBoost(train_labels, train_features, test_labels, test_features, anticlust, colname):
	parameters = {
			"n_estimators": [5, 15, 25, 35, 50, 75, 100, 150, 200, 250, 300, 400, 550, 750, 1000], 
			"learning_rate": [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.75, 1.0]
			}	
	grid = GridSearchCV(estimator=AdaBoostRegressor(),
			cv=KFold(n_splits=10),
			param_grid=parameters,
			n_jobs=-1,
			scoring="neg_mean_squared_error")
	model = grid.fit(X = train_features, y = train_labels)
	hyperparameters = model.best_estimator_.get_params()
	predictions = model.predict(test_features)
	errors = pd.DataFrame({
					'dimension':[colname],
					'anticlust':[anticlust],
					'mse':[mean_squared_error(test_labels,predictions)],
	            	'mae':[mean_absolute_error(test_labels,predictions)],
	            	'cod':[r2_score(test_labels,predictions)],
	            	'rate':[hyperparameters['learning_rate']],
	            	'estimators':[hyperparameters['n_estimators']]
	            })
	importance = model.best_estimator_.feature_importances_
	importance = pd.DataFrame(importance.reshape(1,-1), columns=train_features.columns, index=[f'{colname}X{anticlust}'])
	predictions = pd.DataFrame({'sub': predict.iloc[test_features.index,]['sub'].values, f'{colname}': predictions})
	return [predictions, errors, importance, model]

def InnerGTBoost(train_labels, train_features, test_labels, test_features, anticlust, colname, weight_obs):
	parameters = {"n_estimators": list(range(100, 1200, 100)),"max_depth": list(range(3, 8, 1))}
	grid = GridSearchCV(estimator=GradientBoostingRegressor(),
				cv=KFold(n_splits=10),
				param_grid=parameters,
				n_jobs=-1,
				scoring="neg_mean_squared_error")
	if (weight_obs == False):
		model = grid.fit(X = train_features, y = train_labels)
	if (weight_obs == True):
		kde_grid = KernelDensity(kernel="gaussian",bandwidth=0.25)
		kde_grid.fit(train_labels.reshape(-1,1))
		weights = 1-np.exp(kde_grid.score_samples(train_labels.reshape(-1,1)))
		model = grid.fit(X = train_features, y = train_labels, sample_weight = weights)
	hyperparameters = model.best_estimator_.get_params()
	predictions = model.predict(test_features)
	errors = pd.DataFrame({
					'dimension':[colname],
					'anticlust':[anticlust],
					'mse':[mean_squared_error(test_labels,predictions)],
	            	'mae':[mean_absolute_error(test_labels,predictions)],
	            	'cod':[r2_score(test_labels,predictions)],
	            	'depth':[hyperparameters['max_depth']],
	            	'estimators':[hyperparameters['n_estimators']]
	            })
	importance = model.best_estimator_.feature_importances_
	importance = pd.DataFrame(importance.reshape(1,-1), columns=train_features.columns, index=[f'{colname}X{anticlust}'])
	predictions = pd.DataFrame({'sub': predict.iloc[test_features.index,]['sub'].values, f'{colname}': predictions})
	return [predictions, errors, importance, model]

def FindModels(root_path):
	MODELS = []
	for root, dirs, files in os.walk(root_path):
		for file in files:
			if file.endswith(".dat"):
				file_path = os.path.join(root, file)
				MODELS.append(file_path)
	return MODELS

######
### Main Analysis With All Neuroimaging Features 
######

#Read in Datasets
MODE=sys.argv[1].split('-')[0]
SUBSET=sys.argv[1].split('-')[1]
root_dir = "/scratch/jirsaraie/proj22-NeuroMap"
work_dir=f"{root_dir}/analysis/{TODAY}_ADA/{MODE}-{SUBSET}"
df = pd.read_csv(f'{root_dir}/data/freeze_sens_755x2340_20230601.csv')
featmap = pd.read_csv(f'{root_dir}/data/featmap_2319x5_20230531.csv')
os.makedirs(work_dir,exist_ok=True)

#Store Results For All Models
predict=df['sub'].to_frame()
perform = pd.DataFrame(columns=['mse','mae','cod','estimators','rate'])
explain = pd.DataFrame(index=range(0),columns=range(df.shape[1]-22)) ; explain.columns = df.columns[22:df.shape[1]]

#Train the Psychopathology Models Within Each AntiCluster
for anticlust in range(1,6):
	for colname in [col for col in df.columns if '_score' in col]:
		train_labels, train_features, test_labels, test_features = OuterCVLoop(df, anticlust, colname) 
		predictions, errors, importance, model = InnerAdaBoost(train_labels, train_features, test_labels, test_features, anticlust, colname)
		predict = pd.merge(predict, predictions, left_on='sub', right_on='sub', how='outer')
		perform = pd.concat([perform, errors], sort=False)
		try:
			explain = pd.concat([explain, importance],sort=False)
		except TypeError:
   			 explain = importance
		joblib.dump(model,f"{work_dir}/model_{colname}_{anticlust}.dat")
		os.chmod(f"{work_dir}/model_{colname}_{anticlust}.dat",777)

#Save Model Predictions, Performance, and Explainability
PREDICT=predict ; predict=pd.DataFrame(predict['sub']) 
for dimension in ['pfact','inter','exter']:
	if PREDICT.filter(regex=dimension).shape[1] != 5:
		print(f"! Critical Error Detected !") ; #sys.exit()
		with open(f"{work_dir}/ERROR_DETECTED.txt", 'w') as file:
			pass
	else:
		predict[f'PREDICT-{dimension}_{os.path.basename(work_dir)}'] = PREDICT.filter(regex=dimension).mean(axis=1)

predict.to_csv(f'{work_dir}/frozen-predict.csv',index=False) 
perform=perform[['dimension','anticlust','mse','mae','cod','rate','estimators']]
perform.to_csv(f'{work_dir}/frozen-perform.csv',index=False) 
explain = explain.dropna(axis=1, how='all')
explain['dimension'] = explain.index.str.split("X").str.get(0)
explain['anticlust'] = explain.index.str.split("X").str.get(1)
order = ['dimension', 'anticlust'] + [col for col in explain.columns if col not in ['dimension', 'anticlust']]
explain.reindex(columns=order).to_csv(f'{work_dir}/frozen-explain.csv',index=False) 

######
### Perform Pseudo-Bootstrapping with the Resample Package
######

bootstrap_statistics = pd.DataFrame() ; STRATIFIED="FALSE"
for MODEL in FindModels(work_dir):
	print(MODEL) 
	model=joblib.load(MODEL)
	MODE=MODEL.split('/')[-2].split('-')[0]
	SUBSET=MODEL.split('/')[-2].split('-')[1]
	colname=os.path.basename(MODEL).split('_')[1] + "_score"
	anticlust=int(os.path.basename(MODEL).split('_')[3].replace(".dat", ""))
	train_labels, train_features, test_labels, test_features = OuterCVLoop(df, anticlust, colname)
	for i in range(0,2000):
		if STRATIFIED == "TRUE":
			resampled_features, resampled_labels = resample(test_features, test_labels, replace=True)
		else:
			sample_indices = np.random.choice(len(test_features), size=len(test_features), replace=True)
			resampled_features = test_features.values[sample_indices]
			resampled_labels = test_labels[sample_indices]
		predictions = model.predict(resampled_features)
		errors = pd.DataFrame({
					'dimension':[colname],
					'anticluster':[anticlust],
					'mse':[mean_squared_error(resampled_labels,predictions)],
	            	'mae':[mean_absolute_error(resampled_labels,predictions)],
	            	'cod':[r2_score(resampled_labels,predictions)]
				}); bootstrap_statistics = pd.concat([bootstrap_statistics, errors])
bootstrap_statistics.to_csv(os.path.dirname(MODEL)+'/frozen-resampl.csv',index=False) 

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######