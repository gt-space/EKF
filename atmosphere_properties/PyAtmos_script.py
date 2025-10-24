import numpy as np
import csv


from pyatmos import download_sw_nrlmsise00, read_sw_nrlmsise00
swfile = download_sw_nrlmsise00() 
swdata = read_sw_nrlmsise00(swfile) 


from pyatmos import nrlmsise00


def density(alt, t, lat, lon, swdata):
    nr_output = nrlmsise00(t, (lat, lon, alt), swdata)
    density = nr_output.rho
    return density

def temperature(alt, t, lat, lon, swdata):
    nr_output = nrlmsise00(t, (lat, lon, alt), swdata)
    temperature = nr_output.T
    return temperature

def molar_mass(alt, t, lat, lon, swdata):
    nr_output = nrlmsise00(t, (lat, lon, alt), swdata)
    composition = nr_output.nd.keys()
    abundances = nr_output.nd
    total = 0.0
    molar_mass = 0.0
    # weights from Pubchem
    weights = {'He': 0.004002602, 'O': 0.0159994, 'N2': 0.028014, 'O2': 0.031999, 'Ar': 0.039948, 'H': 0.00100794, 'N': 0.01400674, 'ANM O': 0.0159994}
    for species in composition:
        total += abundances[species]
    for species in composition:
        molar_mass += weights[species] * abundances[species] / total
    return molar_mass


# Set a specific time and location (time in year/month/date 24 hour time)
t = '2025-01-12 14:00:00' # time(UTC) 
lat, lon = 35.3478813171387, -117.806831359863 # latitude, longitude in [degrees]
alt = 0 # initial altitude, in km
n = 0 # initial step n
alt_step = 2 # the altitude, in km, between each step
n_max = 100 # number of steps to take (the array size will be n_max + 1, and the max altitude will be alt_step * n_max)
altitude_list = []
density_list = []
temperature_list = []
weight_list = []


while True:
    altitude_list += [alt]
    density_list += [density(alt, t, lat, lon, swdata)]
    temperature_list += [temperature(alt, t, lat, lon, swdata)]
    weight_list += [molar_mass(alt, t, lat, lon, swdata)]
    alt += alt_step
    n += 1
    if n > n_max:
        break


with open('pyatmos_output.csv', 'w') as f:
    write = csv.writer(f)
    write.writerow(altitude_list)
    write.writerow(altitude_list)
    write.writerow(density_list)
    write.writerow(temperature_list)
    write.writerow(weight_list)


print(len(altitude_list)) # prints the lenngth of the array size for confirmation
print(alt-alt_step) # prints the max altitude for confirmation


