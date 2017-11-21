import numpy as np
import matplotlib.pyplot as plt
from util import *
from star import Star
from dust import Dust
from scipy.integrate import simps, quad

# choose whether to display results or not
SHOW = [0,0,0,1,0]

# the two distances to check
DIST_1 = 10.0 * AU
DIST_2 = 130.0 * AU
# the frequency range to go over for the stellar spectra
STAR_FREQ = 3 * np.logspace(13, 15, num=10000000)

# PART 1
fomalhaut = Star('Fomalhaut',1.92*M_SUN,1.842*R_SUN,8590.0)
spectrum_1 = fomalhaut.spectrum(DIST_1)
spectrum_2 = fomalhaut.spectrum(DIST_2)
# plot the resulting spectra
if SHOW[0]:
    plt.figure()
    plt.semilogx(STAR_FREQ, spectrum_1(STAR_FREQ))
    plt.title('Spectrum of ' + fomalhaut.name() 
        + ' at ' + str(ftoi(DIST_1/AU)) + ' AU')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Flux Density (Jy)')
    plt.tight_layout()
    plt.figure()
    plt.semilogx(STAR_FREQ, spectrum_2(STAR_FREQ))
    plt.title('Spectrum of ' + fomalhaut.name()
        + ' at ' + str(ftoi(DIST_2/AU)) + ' AU')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Flux Density (Jy)')
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
        a = ftoi(dust[0].axis() / AU)
        radius = ftoi(dust[0].radius() * 10**6)
        print('Power absorbed by', radius, 
            'µ dust at distance', a, 'AU:', sig_figs(dust[1],4), 'W')

# PART 3
DUST_FREQ = dusts[0][0].freq()
if SHOW[2]:
    print('\nPART 3')
    for dust in dusts:
        a = ftoi(dust[0].axis() / AU)
        radius = ftoi(dust[0].radius() * 10**6)
        print('Equilibrium temperature of', radius, 
            'µ dust at distance', a, 'AU:', sig_figs(dust[0].temp(),4), 'K')
    freq = DUST_FREQ[:round(len(DUST_FREQ)/2)]
    plt.figure()
    for dust in dusts[0:4]:
        a = ftoi(dust[0].axis() / AU)
        radius = ftoi(dust[0].radius() * 1e6)
        lumi = dust[0].spectrum()(DUST_FREQ)
        lumi = lumi[:round(len(lumi)/2)]
        plt.semilogx(freq, lumi)
    plt.title('Spectrum of '+str(radius)+' µ dust at '+str(a)+' AU')
    plt.xlabel('Wavelength (µm)')
    plt.ylabel('Flux (Jy)')
    plt.legend(('0.1 Micron','1.0 Micron','10.0 Micron','1000 Micron'),
        loc='best',prop={'size':10})
    plt.tight_layout()
    plt.figure()
    for dust in dusts[4:]:
        a = ftoi(dust[0].axis() / AU)
        radius = ftoi(dust[0].radius() * 1e6)
        lumi = dust[0].spectrum()(DUST_FREQ)
        lumi = lumi[:round(len(lumi)/2)]
        plt.semilogx(freq, lumi)
    plt.title('Spectrum of '+str(radius)+' µ dust at '+str(a)+' AU')
    plt.xlabel('Wavelength (µm)')
    plt.ylabel('Flux (Jy)')
    plt.legend(('0.1 Micron','1.0 Micron','10.0 Micron','1000 Micron'),
        loc='best',prop={'size':10})
    plt.tight_layout()

# PART 4
if SHOW[3]:
    print('\nPART 4')
    # distance from us to fomalhaut
    distance = 7.70 * PC
    # TODO: right numbers
    num = np.array([1e35, 1e32, 1e030, 1e26, 1e38, 1e35, 1e32, 1e28])
    spectra = []
    for i, dust in enumerate(dusts):
        factor = (dust[0].radius()/distance)**2
        spectra.append(num[i] * factor * dust[0].spectrum()(DUST_FREQ))
    total = np.sum(spectra,axis=0)
    # TODO: is this the right conversion from freq to wl?
    WL = C/DUST_FREQ
    def half(arr):
        return arr[:round(len(arr)/2.5)]
    plt.figure()
    plt.semilogx(half(C/DUST_FREQ)*1e6, half(total))
    for spectrum in spectra:
        plt.semilogx(half(C/DUST_FREQ)*1e6, half(spectrum))
    plt.legend(('Total','0.1µ, 10 AU','1µ, 10 AU','10µ, 10 AU','1mm 10 AU',
        '0.1µ, 130 AU','1µ, 130 AU','10µ, 130 AU','1mm, 130 AU'),
        loc='best',prop={'size':10})
    for i, dust in enumerate(dusts):
        a = ftoi(dust[0].axis() / AU)
        radius = ftoi(dust[0].radius() * 10**6)
        print('Mass of', radius, 'µ dust at distance', a, 'AU:', 
            sig_figs(num[i]*dust[0].mass()), 'kg')

# PART 5
dust_times = []
for dust in dusts:
    pr = dust[0].pr_drag()
    # TODO: does RP need gravity adjustment?
    rp = dust[0].radiation_pressure()
    force = rp - pr
    # PR > RP -> the dust will be pulled in to the star
    if force < 0:
        radius = fomalhaut.radius()
        pull = dust[0].pr_drag(radius) - dust[0].radiation_pressure(radius)
        # the dust stabilizes before reaching the star
        if pull < 0:
            dust_times.append((dust[0], 'inf'))
        # TODO: the dust gets pulled into the star
        else:
            dust_times.append((dust[0], 0))
    # PR = RP -> the dust is stable where it is
    elif force == 0:
        dust_times.append((dust[0], 'inf'))
    # TODO PR < RP -> the dust will be blown out
    else:
        dust_times.append((dust[0], 0))
if SHOW[4]:
    print('\nPART 5')
    for dust_time in dust_times:
        dust = dust_time[0]
        time = dust_time[1]
        a = ftoi(dust.axis() / AU)
        radius = ftoi(dust.radius() * 10**6)
        print('Timescale to remove', radius, 'µ dust at',
            'initial semimajor axis', a, 'AU:', time, 's')

plt.show()