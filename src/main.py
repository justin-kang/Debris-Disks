import numpy as np
import matplotlib.pyplot as plt
from util import *
from star import Star
from dust import Dust
from scipy.integrate import simps, quad

# choose whether to display results or not
SHOW = 1

# the two distances to check
DIST_1 = 10 * AU
DIST_2 = 130 * AU
# the frequency range to go over for the stellar spectra
STAR_FREQ = 3 * np.logspace(13, 15, num=10000000)

# PART 1
fomalhaut = Star('Fomalhaut',1.842*R_SUN,8590)
spectrum_1 = fomalhaut.spectrum(DIST_1)
spectrum_2 = fomalhaut.spectrum(DIST_2)
# plot the resulting spectra
if SHOW:
    plt.figure()
    plt.plot(STAR_FREQ, spectrum_1(STAR_FREQ))
    plt.title('Spectrum of ' + fomalhaut.name() 
        + ' at ' + str(ftoi(DIST_1/AU)) + ' AU')
    plt.xlabel('Frequency (Hz)')
    plt.tight_layout()
    plt.figure()
    plt.plot(STAR_FREQ, spectrum_2(STAR_FREQ))
    plt.title('Spectrum of ' + fomalhaut.name()
        + ' at ' + str(ftoi(DIST_2/AU)) + ' AU')
    plt.xlabel('Frequency (Hz)')
    plt.tight_layout()
    plt.show()

# PART 2
print('PART 2\n')
dusts = []
dust_1 = Dust(1000,DIST_1)
dusts_1 = [Dust(0.1,DIST_1), Dust(1,DIST_1), Dust(10,DIST_1)]
dusts.append((dust_1, dust_1.power(fomalhaut)))
for dust in dusts_1:
    dusts.append((dust, dust.power(fomalhaut,False)))
dust_2 = Dust(1000,DIST_2)
dusts_2 = [Dust(0.1,DIST_2), Dust(1,DIST_2), Dust(10,DIST_2)]
dusts.append((dust_2, dust_2.power(fomalhaut)))
for dust in dusts_2:
    dusts.append((dust, dust.power(fomalhaut,False)))
if SHOW:
    for dust in dusts:
        a = ftoi(dust[0].semimajor_axis() / AU)
        radius = ftoi(dust[0].radius() * 10**6)
        print('Power absorbed by', radius, 
            'µ dust at distance', a, 'AU:', sig_figs(dust[1],4), 'W')

# PART 3
print('PART 3\n')
# the frequency range to go over for the dust spectra
DUST_FREQ = 3 * np.logspace(15, 17, num=10000000)
temp_dusts = []
for dust in dusts:
    temp = dust[0].temp(dust[1],fomalhaut)
    spectrum = dust[0].spectrum(dust[1],fomalhaut)
    temp_dusts.append((dust[0], dust[1], temp, spectrum))
dusts = temp_dusts
if SHOW:
    for dust in dusts:
        a = ftoi(dust[0].semimajor_axis() / AU)
        radius = ftoi(dust[0].radius() * 10**6)
        print('Equilibrium temperature of', radius, 
            'µ dust at distance', a, 'AU:', sig_figs(dust[2],4), 'K')
    '''
    for dust in dusts:
        a = ftoi(dust[0].semimajor_axis() / AU)
        radius = ftoi(dust[0].radius() * 10**6)
        plt.figure()
        plt.plot(DUST_FREQ, dust[3](DUST_FREQ))
        plt.title('Spectrum of ' + str(radius) + 'µ dust at ' + str(a) + 'AU')
        plt.xlabel('Frequency (Hz)')
        plt.tight_layout()
    plt.show()
    '''

# PART 4
'''
print('PART 4\n')
Go dig up an infrared+millimeter SED of Fomalhaut’s two (warm+cold) debris 
belts out of the literature, and compare it to the SEDs computed in Part 3. 
For each dust grain size/type, how many dust grains are needed to roughly 
match the observed SED of each component? What is the corresponding dust mass, 
if you assume the Fomalhaut debris ring is made out of grains of that size and 
that dust grains have a density somewhere around 2 g/cc? Which dust grain size 
produces the best fit to the SED shape?
'''


# PART 5
'''
print('PART 5\n')
In class we talked about radiation pressure and Poynting-Robertson drag in 
terms of ideal blackbodies. However, for previous problems you’ve already 
calculated the actual energy absorption and emission of real particles, 
which should allow you to calculate better versions of each force. For each 
particle size, calculate the real radiation pressure and Poynting-Robertson 
drag at each radius. Estimate a characteristic timescale for dust grains of 
each size to be removed from the system, assuming they start at 10 or 130 AU.
'''
