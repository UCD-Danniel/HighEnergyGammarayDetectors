import argparse as ap
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import glob


def ParseArgs():
    parser = ap.ArgumentParser()
    parser.add_argument('folder', help='Type in the folder containing spectra.')
    return parser.parse_args()

## Fit plots
def gaussprofile(x, mu, sig, a):
    return a * np.exp(-((x - mu) ** 2) / (2 * sig ** 2))

def lineplot(x, m, c):
    return m * x + c

def respoly(x, a, b, c):
    return a + (b * x) + (c * x**2)

def efflog(x, a, b, c):
    ln_x = np.log(x)
    return a + b * ln_x + c * (ln_x)**2

## Plot functions
def channelenergyplot(x, y, err):

    popt, pcov = curve_fit(lineplot, x, y)
    m, c = popt
    m_err = np.sqrt(pcov[0,0])
    print(f'Channel-energy relation: ({m:.2e} +/- {m_err:.2e}) keV/channel')

    plt.plot(x, lineplot(x, *popt), label='Best Fit Line')
    plt.errorbar(x, y, xerr=err, fmt='o', color='r', label=f'{detector} Data')  
    plt.xlabel("Channel")
    plt.ylabel("Energy (keV)")
    plt.title(f"{detector} Channel-energy Plot")
    plt.minorticks_on()
    plt.grid()
    plt.legend()
    plt.savefig(f'{pfolder}/{detector}_channelenergy.png') 
    plt.show()

    
def resolutionplot(x, y, err):
    popt, pcov = curve_fit(respoly, x, y)
    a, b, c = popt

    x_fit = np.linspace(min(x), max(x), 200)  
    y_fit = respoly(x_fit, *popt)  
    
    plt.plot(x_fit, y_fit, label='Best Fit Curve')
    plt.errorbar(x, y, yerr=err, fmt='o', color='orange', label=f'{detector} Data')  
    plt.xlabel("Energy (keV)")
    plt.ylabel("Percent Resolution (%)")
    plt.title(f"{detector} Resolution-energy plot")
    plt.grid()
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.savefig(f'{pfolder}/{detector}_resolutionenergy.png')
    plt.show()


def efficiencyplot(x, y, err):  
    
    popt, pcov = curve_fit(efflog, x, y)
    a, b, c = popt

    x_fit = np.linspace(min(x), max(x), 200)  
    y_fit = efflog(x_fit, *popt)  
    
    plt.plot(x_fit, y_fit, label='Best Fit Curve')
    plt.errorbar(x, y, yerr=err, fmt='o', color='g', label=f'{detector} Data')  
    plt.xlabel("Energy (keV)")
    plt.ylabel("Efficiency")
    plt.title(f"{detector} Efficiency-energy Plot")
    plt.minorticks_on()
    plt.xscale('log')
    #plt.yscale('log')
    plt.grid()
    plt.legend()
    plt.savefig(f'{pfolder}/{detector}_efficiencyenergy.png')
    plt.show()


def effangleplot(x, y, n):
    
    abseffangle = []
    abseffangle_err = []
    onabseffangle = 0
    onabseffangle_err = 0
    
    for i, val in enumerate(y):
        en = np.arange(len(val))
        initial_guess = [np.mean(en), np.std(en), max(val)] 
        bounds = [[0, 0, 0],  # Lower bounds
                 [np.inf, np.inf, np.inf]]  # Upper bounds 
            
        mask = (en >= masklower[n]) & (en <= maskhigher[n])  
            
        popt, pcov = curve_fit(gaussprofile, en[mask], val[mask], p0=initial_guess, bounds=bounds)     
        
        _, sig, a = popt

        sig_err = np.sqrt(pcov[1,1])
        a_err = np.sqrt(pcov[2,2])            
        
        # Absolute Efficiency                                 
        area = a*sig*np.sqrt(2*np.pi)   
        count_rate = area / time[n]
        total_ph = activity[n] * phyield[n] 
        
        abs_eff = count_rate / total_ph
        abs_eff_err = abs_eff * (a_err / a)
        abseffangle.append(abs_eff)
        abseffangle_err.append(abs_eff_err)
        if i == 0:
            onabseffangle = abs_eff
            onabseffangle_err = abs_eff_err
            
    y = np.array(abseffangle) / onabseffangle
    err = np.array(abseffangle_err)

    popt, pcov = curve_fit(lineplot, x, y)
    m, c = popt
    m_err = np.sqrt(pcov[0,0])
    print(f'Efficiency-angle relation: ({m:.2e} +/- {m_err:.2e}) efficiency/angle')

    plt.plot(x, lineplot(x, *popt), color='cornflowerblue', label='Best Fit Line')    
    plt.errorbar(x, y, yerr=err, fmt='o', color='orange', label=f'{source[n]} Data')
                
    plt.xlabel("Angle (degrees)")
    plt.ylabel("Efficiency ratio (Off/On)")
    plt.title(f"{detector} Efficiency-Angle Plot")
    plt.minorticks_on()
    plt.grid()
    plt.legend()
    plt.savefig(f'{pfolder}/{detector}_effangle.png') 
    plt.show()

def spectraplot(vals):
    """ Spectral data Plot Function""" 
    channels = []
    channels_err = []
    fwhmlist = []
    fwhmlist_err = []
    abseff = []
    abseff_err = []
    for i, val in enumerate(vals):
        x = np.arange(len(val))
        plt.figure(figsize = (9,5)) 
        plt.plot(x, val, label="Spectrum")
        plt.title(f'Spectrum of {source[i]} using {detector} detector with gaussian profile')
        plt.xlabel('Energy (keV)')
        plt.ylabel('Counts')
        plt.grid()
    
        initial_guess = [np.mean(x), np.std(x), max(val)] 
        bounds = [[0, 0, 0],  # Lower bounds
                  [np.inf, np.inf, np.inf]]  # Upper bounds 
        
        mask = (x  >= masklower[i]) & (x <= maskhigher[i])  
        
        popt, pcov = curve_fit(gaussprofile, x[mask], val[mask], p0=initial_guess, bounds=bounds)     
              
        # Assigning popts to variables
        mu, sig, a = popt
        # Assigning pcovs to variables
        mu_err = np.sqrt(pcov[0,0])
        sig_err = np.sqrt(pcov[1,1])
        a_err = np.sqrt(pcov[2,2])
        channels.append(mu)
        channels_err.append(mu_err)
            
        # FWHM
        fwhm = 2 * sig * np.sqrt(2*np.log(2))
        fwhm_err = fwhm * (sig_err / sig)
        fwhmlist.append(fwhm)
        fwhmlist_err.append(fwhm_err)
        
        # Absolute Efficiency                                 
        area = a*sig*np.sqrt(2*np.pi)   
        count_rate = area / time[i]
        total_ph = activity[i] * phyield[i] 
        
        abs_eff = count_rate / total_ph
        abs_eff_err = abs_eff * (a_err / a)
        abseff.append(abs_eff)
        abseff_err.append(abs_eff_err)
        
        plt.plot(x, gaussprofile(x, *popt), c='r', label='Gauss Profile')
        plt.legend()
        plt.show()

    # Resolution
    fwhmlist = np.array(fwhmlist)
    fwhmlist_err = np.array(fwhmlist_err)
    
    reslist = 100*(fwhmlist / energy)
    reslist_err = reslist * (fwhmlist_err / fwhmlist)
    
    # Intrinsic Efficiency
    d = 15 
    solidangle = 2*np.pi* (1- d/np.sqrt(d**2 + ar**2))
    abseff = np.array(abseff)   
    abseff_err = np.array(abseff_err)
    
    inteff = abseff * (4*np.pi / solidangle)
    inteff_err = inteff * (abseff_err / abseff)

    return channels, channels_err, reslist, reslist_err, abseff, abseff_err, inteff, inteff_err


def main():
    global detector
    global source 
    global energy 
    global sfile  
    global phyield 
    global activity 
    global time 
    global masklower 
    global maskhigher 
    global ar
    global ameffanglei
    global cseffanglei    
    global pfolder
    
    args = ParseArgs()
    
    # File path
    folder = args.folder
    pfolder = f'{folder}/PlotFolder'
    detector = str(input("Please input detector [bgo, nai, cdte]: "))
    
    # Get list of files
    files = glob.glob(f'{folder}/{detector}_*')
    bg = glob.glob(f'{folder}/Background_{detector}*')
    bg = bg[0]
        
    ### Background ###
    bg_info = []
    bg_data = []
        
    with open(bg, 'r') as file:
        indicator = False
        for line in file:
            if not indicator:
                bg_info.append(line)
            if '0 1023' in line or '<<DATA>>' in line:
                indicator = True
                continue
        
            if indicator:
                if '$ROI:' in line or '<<END>>' in line:
                    break
                else:
                    bg_data.append(int(line))
            
    bg_data = np.asarray(bg_data)
        
        
    ### Source Spectra ###
    angles = np.array([0, 105, 120, 15, 30, 45, 60, 75, 90])
    # Different spectra sorted by source   
    csfiles = []   
    amfiles = []
    bafiles = []
    cofiles = []
        
    for i, file_path in enumerate(files):
        info = [] # header info
        data = [] # data block
    
        parts = file_path.split('_')
        detec, source, angle = parts
        
        with open(file_path, 'r') as file:
             indicator = False
             for line in file:
                if not indicator:
                    info.append(line)
        
                if '0 1023' in line or '<<DATA>>' in line:
                    indicator = True
                    continue
        
                if indicator:
                    if '$ROI:' in line or '<<END>>' in line:
                        break
                    else:
                        data.append(int(line)) 
        data = np.asarray(data)
    
        if source == 'CS':
            csfiles.append(data)
        elif source == 'AM':
            amfiles.append(data)
        elif source == 'BA':
            bafiles.append(data)
        elif source == 'CO':
            cofiles.append(data) 
    
    # Photopeaks, source info order
    if detector == 'bgo':
        detector = 'BGO'
        source = ['Am', 'Ba', 'Ba', 'Cs', 'Co', 'Co']
        energy = np.array([59.5409, 80.9979, 356.0129, 661.657, 1173.228, 1332.492])
        sfile = [amfiles[0], bafiles[0], bafiles[0], csfiles[0], cofiles[0], cofiles[0]]
        phyield = [0.3578, 0.329, 0.6205, 0.8499, 0.9985, 0.999826]
        activity = [412550, 20350, 20350, 162430, 1147, 1147]
        time = [150, 300, 300, 150, 300, 300]
        masklower = [20, 30, 170, 320, 560, 650]
        maskhigher = [30, 40, 200, 350, 600, 680]
        ar = 2.94      # detector face radius
        ameffanglei = 0
        cseffanglei = 3
    
    elif detector == 'nai':
        detector = 'NaI(Ti)'
        source = ['Am', 'Am', 'Ba', 'Ba', 'Cs']
        energy = np.array([26.3446, 59.541, 80.997, 356.013, 661.657])
        sfile = [amfiles[0], amfiles[0], bafiles[0], bafiles[0], csfiles[0]]
        phyield = [0.024, 0.3578, 0.329, 0.6205, 0.8499,]
        activity = [412550, 20350, 20350, 20350, 162430]
        time = [150, 150, 300, 300, 150]        
        masklower = [10, 25, 45, 180, 325]
        maskhigher = [25, 40, 55, 210, 400]
        ar = 2.94
        ameffanglei = 1
        cseffanglei = 4
    
    elif detector == 'cdte':
        detector = 'CdTe'
        source = ['Am', 'Am', 'Am', 'Ba']
        energy = np.array([26.3446, 33.541, 59.541, 80.997])
        sfile = [amfiles[0], amfiles[0], amfiles[0], bafiles[0]]
        phyield = [0.024, 0.00121, 0.3578, 0.329]
        activity = [412550, 412550, 412550, 20350]          
        time = [150, 150, 150, 300]        
        masklower = [220, 310, 1150, 1550]
        maskhigher = [300, 370, 1250, 1650]
        ar = 0.5        
        ameffanglei = 2
        cseffanglei = None
    
    ## Plots
    channels, channels_err, reslist, reslist_err, abseff, abseff_err, inteff, inteff_err = spectraplot(sfile)
    channelenergyplot(np.array(channels), energy, np.array(channels_err))
    resolutionplot(energy, reslist, reslist_err)
    efficiencyplot(energy, inteff, inteff_err)
    effangleplot(angles, amfiles, ameffanglei)
    if ameffanglei < 2:
        effangleplot(angles, csfiles, cseffanglei)

if __name__ == "__main__":
    main()