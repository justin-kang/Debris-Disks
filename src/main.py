import numpy as np
import matplotlib.pyplot as plt
from util import *
from star import Star
from dust import Dust

# choose whether to display results of the 5 parts, then graphs
SHOW = [1,1,1,1,1,1]
# the two distances to check
DIST_1 = 10.0 * AU
DIST_2 = 130.0 * AU

# PART 1
fomalhaut = Star('Fomalhaut',1.92*M_SUN,1.842*R_SUN,8590.0)
# distance from us to fomalhaut
distance = 7.70 * PC
# plot Fomalhaut's spectra
if SHOW[0]:
    # the frequency range to go over for the stellar spectra
    STAR_FREQ = 3 * np.logspace(13, 15, num=10000000)
    plt.figure()
    plt.loglog(STAR_FREQ, fomalhaut.spectrum(DIST_1)(STAR_FREQ))
    plt.loglog(STAR_FREQ, fomalhaut.spectrum(DIST_2)(STAR_FREQ))
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Flux Density (Jy)')
    plt.legend(('10 AU','130 AU'),loc='best',prop={'size':10})
    plt.tight_layout()

# PART 2
dusts = []
dusts_1 = [Dust(0.1,fomalhaut,DIST_1), Dust(1,fomalhaut,DIST_1), \
    Dust(10,fomalhaut,DIST_1)]
for dust in dusts_1:
    dusts.append((dust, dust.power()))
dust_1 = Dust(1000,fomalhaut,DIST_1)
dusts.append((dust_1, dust_1.power(perfect=True)))
dusts_2 = [Dust(0.1,fomalhaut,DIST_2), Dust(1,fomalhaut,DIST_2), \
    Dust(10,fomalhaut,DIST_2)]
for dust in dusts_2:
    dusts.append((dust, dust.power()))
dust_2 = Dust(1000,fomalhaut,DIST_2)
dusts.append((dust_2, dust_2.power(perfect=True)))
if SHOW[1]:
    print('PART 2')
    for dust in dusts:
        a = ftoi(dust[0].semimajor_axis() / AU)
        radius = ftoi(dust[0].radius() * 1e6)
        print('Power absorbed by', radius, 
            'µ dust at distance', a, 'AU:', sig_figs(dust[1]), 'W')

# PART 3
# the frequency range to go over for the dust spectra
DUST_FREQ = dusts[0][0].freq()
if SHOW[2]:
    print('\nPART 3')
    for dust in dusts:
        a = ftoi(dust[0].semimajor_axis() / AU)
        radius = ftoi(dust[0].radius() * 1e6)
        print('Equilibrium temperature of', radius, 
            'µ dust at distance', a, 'AU:', sig_figs(dust[0].temp()), 'K')
    # plot the emission spectra of the dust
    freq = DUST_FREQ[:round(len(DUST_FREQ)/2)]
    plt.figure()
    for dust in dusts[0:4]:
        lumi = dust[0].spectrum(distance)(DUST_FREQ)
        lumi = lumi[:round(len(lumi)/2)]
        plt.loglog(freq, lumi)
    plt.xlabel('Wavelength (µm)')
    plt.ylabel('Flux (Jy)')
    plt.legend(('0.1 Micron','1.0 Micron','10.0 Micron','1000 Micron'),
        loc='best',prop={'size':10})
    plt.tight_layout()
    plt.figure()
    for dust in dusts[4:]:
        lumi = dust[0].spectrum(distance)(DUST_FREQ)
        lumi = lumi[:round(len(lumi)/2)]
        plt.loglog(freq, lumi)
    plt.xlabel('Wavelength (µm)')
    plt.ylabel('Flux (Jy)')
    plt.legend(('0.1 Micron','1.0 Micron','10.0 Micron','1000 Micron'),
        loc='best',prop={'size':10})
    plt.tight_layout()

# PART 4
if SHOW[3]:
    print('\nPART 4')
    ref_peaks = np.array([10,10,10,10,10**(-.5),10**(-.5),10**(-.5),10**(-.5)])
    mass = []
    spectra = []
    for dust in dusts:
        mass.append(dust[0].mass())
        spectra.append(dust[0].spectrum(distance)(DUST_FREQ))
    exp_peaks = np.amax(spectra,axis=1)
    num = np.divide(ref_peaks, exp_peaks)
    masses = np.multiply(mass, num)
    for i in range(len(dusts)):
        a = ftoi(dusts[i][0].semimajor_axis() / AU)
        radius = ftoi(dusts[i][0].radius() * 1e6)
        print('Needed number of', radius, 'µ dust at distance', a, 'AU:', 
            sig_figs(num[i]))
        print('Needed mass of', radius, 'µ dust at distance', a, 'AU:',
            sig_figs(masses[i]), 'kg')

# PART 5
if SHOW[4]:
    print('\nPART 5')
    mstar = fomalhaut.mass()
    YR = 60 * 60 * 24 * 365
    for dust in dusts:
        pr = dust[0].pr_drag()
        rp = dust[0].radiation_pressure()
        time = dust[0].removal_time() / YR
        a = ftoi(dust[0].semimajor_axis() / AU)
        radius = ftoi(dust[0].radius() * 1e6)
        print('Poynting-Robertson drag on', radius, 'µ dust at', a, 'AU:', 
            sig_figs(pr), 'N')
        print('Radiation pressure on', radius, 'µ dust at', a, 'AU:',
            sig_figs(rp), 'N')
        print('Time to remove', radius, 'µ dust at semimajor axis', 
            a, 'AU:', sci_not(time), 'yr')

# GRAPHS
if SHOW[5]:
    plt.show()