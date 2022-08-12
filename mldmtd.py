import numpy as np
import pandas as pd

import gsw
import netCDF4

from scipy.optimize import curve_fit
from scipy.stats import linregress
from sklearn.metrics import r2_score
from netCDF4 import Dataset

########################
# Get data from netcdf
########################
def getProfileData(file):
	try:
		nc = Dataset(file)
		JULD = nc.variables['JULD']
		JULD = netCDF4.num2date(JULD[0],JULD.units)
		LATITUDE = nc.variables['LATITUDE'][0]
		LONGITUDE = nc.variables['LONGITUDE'][0]
		PRES_ADJUSTED = nc.variables['PRES_ADJUSTED'][0]
		TEMP_ADJUSTED = nc.variables['TEMP_ADJUSTED'][0]
		PSAL_ADJUSTED = nc.variables['PSAL_ADJUSTED'][0]
		nc.close()
	except:
		return pd.DataFrame([])

	n = len(PRES_ADJUSTED)
	df = pd.DataFrame({'date': np.repeat(JULD, n), 'latitude': np.repeat(LATITUDE, n), 'longitude': np.repeat(LONGITUDE, n), 'pres': PRES_ADJUSTED, 'psal': PSAL_ADJUSTED, 'temp': TEMP_ADJUSTED})
	df = df.dropna()

	########################
	# Calculate absolute salinity and conservative temperature
	########################
	df['asal'] = gsw.SA_from_SP(df.psal, df.pres, df.longitude, df.latitude)
	df['ctemp'] = gsw.CT_from_t(df.asal, df.temp, df.pres)
	return df

########################
# Sigmoid function
########################
def fsigmoid(x, a, b):
    return 1.0 / (1.0 + np.exp(-a*(x-b)))

########################
# Normalization
########################
def norm(y, y_min, y_max):
	return (y - y_min)/(y_max-y_min)

########################
# Denormalization
########################
def unnorm(y, y_min, y_max):
	return y*(y_max-y_min)+y_min

########################
# Calculate thermocline
########################
def thermocline(df, m_precision=0.01, threshold=0.2):
	########################
	# parameters:
	# 	df: Dataset, columns=['pres','asal','ctemp']
	# 	m_precision: Accuracy in meters
	#	threshold: Threshold for determining MLD and MTD
	########################
	x_true = df.pres.to_numpy()
	y_true = df.ctemp.to_numpy()

	########################
	# Calculates the buoyancy frequency squared (N2)
	########################
	N2, p_mid = gsw.Nsquared(df.asal, df.ctemp, df.pres)
	N2 = np.vstack((N2, p_mid)).T

	########################
	# Obtain the pressure range where N2 is greater
	########################
	try:
		index = np.argmax(N2[:,0])
	except:
		return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

	max_pres = x_true[index]*2
	index = x_true < max_pres
	x_true = x_true[index] 
	y_true = y_true[index]
	N2 = N2[index[:-1]]

	########################
	# Normalize data
	########################
	sign = 1
	if (np.mean(y_true[:2]) > np.mean(y_true[-2:])):
		y_true = -y_true # It is inverted to resemble the function
		sign = -1

	y_min = np.min(y_true) # Min before normalizing
	y_max = np.max(y_true) # Max before normalizing
	y_true = norm(y_true, y_min, y_max) # Normalization

	bounds = ([np.min(y_true), np.min(x_true)], # Lower and upper limits of the parameters
			  [np.max(y_true), np.max(x_true)])

	########################
	# Fit data to sigmoid
	########################
	try:
		popt, pcov = curve_fit(fsigmoid, x_true, y_true, method='dogbox', bounds=bounds)
	except:
		return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

	x_pred = np.linspace(np.min(x_true), np.max(x_true), int(np.max(x_true)-np.min(x_true)))
	y_pred = fsigmoid(x_pred, *popt)

	########################
	# Calculate coefficient of determination (R2)
	########################	
	r2 = r2_score(y_true, fsigmoid(x_true, *popt))

	########################
	# Denormalization data
	########################
	y_true = unnorm(y_true, y_min, y_max) # Denormalize real data
	y_pred = unnorm(y_pred, y_min, y_max) # Denormalize predicted data

	########################
	#  Calculate MLD
	######################### 
	pres_10m = 10
	temp_10m = unnorm(fsigmoid(pres_10m, *popt), y_min, y_max)

	pres_mld = np.nan
	temp_mld = np.nan
	for k in np.arange(pres_10m + m_precision, x_true[-1], m_precision): # check every 10 cm
		temp_mld = unnorm(fsigmoid(k, *popt), y_min, y_max)
		if (temp_mld - temp_10m) >= threshold:
			pres_mld = k
			break

	########################
	#  Calculate MTD
	######################### 
	pres_deep = x_true[-1]
	temp_deep = unnorm(fsigmoid(pres_deep, *popt), y_min, y_max)
	
	pres_thermo = np.nan
	temp_thermo = np.nan
	if ~np.isnan(pres_mld):
		for k in np.arange(pres_deep + m_precision, pres_mld, -m_precision): # verificación cada 10 cm
			temp_thermo = unnorm(fsigmoid(k, *popt), y_min, y_max)
			if (temp_deep - temp_thermo) >= threshold:
				pres_thermo = k
				break

	if pres_mld > pres_thermo:
		return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

	return pres_thermo, temp_thermo*sign, pres_mld, temp_mld*sign, r2, N2, x_pred, y_pred*sign
