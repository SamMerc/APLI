#Importing libraries
import astropy.io.fits as pf
import numpy as np
import os
import math as mt
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
import scipy.stats as st
from stingray import Lightcurve
import statistics
from scipy.stats import chi2
from stingray.pulse import epoch_folding_search
from stingray.pulse import get_orbital_correction_from_ephemeris_file


def get_GTIs(data_file):
    '''
    Function to get the Good Time Intervals (GTIs).
    Parameters
    ----------
    :param data_file: file name, containing the source's GTI data file
    Returns
    ----------
    :param new_gtis: tuples list, containing tuples of start and stop time of the GTIs
    '''
    #Retrieving the GTIs
    GTI_phase_file = data_file
    #Re-formatting the GTIs so that we can use them in the Stingray package
    new_gtis=[]
    for i in GTI_phase_file.data:
        new_gtis.append((i['START'], i['STOP']))
    return new_gtis


def get_orbital_correction(data_file, param_file_name):
    '''
    Function to get the event arrival times with orbital correction.
    Parameters
    ----------
    :param data_file: file name, containing the source's event arrival times and 
    other useful data (base time)
    :param param_file_name: string, containing the name for the parameter file used
    for the orbital correction
    Returns
    ----------
    :param correct_orbit_time: array, containing event arrival times with orbital correction
    '''
    #Reference time in MJD
    base_time = data_file.header['MJDREFI']+data_file.header['MJDREFF']
    
    #Converting start time of observations in seconds
    start_time = data_file.header['TSTART']/(24*3600)
    
    #Converting end time of observations in seconds
    end_time = data_file.header['TSTOP']/(24*3600)
    
    #Getting the correct_time function from the parameters defined above 
    #and the orbit_t2.txt file containing the orbital parameters of the source
    correct_time = get_orbital_correction_from_ephemeris_file(base_time+start_time, base_time+end_time, parfile=param_file_name)[0]
    
    #Applying the correct_time function to our event arrival times
    correct_orbit_time = correct_time(data_file.data['TIME'], base_time)
    
    return correct_orbit_time

def get_pulse_freq(f_min, f_max, N, correct_time_array, GTI, seg_size):
    '''
    Function to get the pulse frequency using the event arrival times with orbital correction.
    Parameters
    ----------
    :param f_min: float, the low end of range of frequencies we want to test.
    :param f_max: float, the high end of range of frequencies we want to test.
    :param N: int, the number of frequencies to test between f_min and f_max.
    :param correct_time_array: array, containing the source's event arrival times 
    with orbital correction.
    :param GTI: array, containing the GTI obtained from get_GTIs function.
    :param seg_size: int, length of segments to to be averaged in periodogram.
    Returns
    ----------
    :param correct_L: tuple list, containing a list of frequencies and a list of corresponding power.
    :param best_freq: float, frequency corresponding the the maximum power. 
    '''
    
    #Defining frequencies to try - required for epoch_folding_search
    #Can either take a large range but then the resolution is low
    #Or small range with high resolution -> could maybe do this iteratively
    trial_freqs = np.linspace(f_min, f_max, N)
    
    #Getting the power as a function of frequency from 1. the event arrival times and 
    #2. the event arrival times with orbital correction to show the difference
    correct_L = epoch_folding_search(correct_time_array, trial_freqs, 
                                     segment_size=seg_size, gti=GTI)
    
    best_freq = correct_L[0][np.argmax(correct_L[1])]
    
    return correct_L, best_freq


def segment_time(PI_data, correct_time, N):
    '''
    Function to make N segments from the event arrival times and PI energies.
    Parameters
    ----------
    :param PI_data: array, containing the energy (PI) data for each event arrival time.
    :param correct_time: array, containing the event arrival times with orbital correction.
    :param N: int, the number of segments to make.
    Returns
    ----------
    :param Time_segments: nested lists, containing N segments of event arrival times.
    :param PI_segments: nested lists, containing N segments of energy corresponding to the event arrival time. 
    '''
    #We get the size of bins to use
    seg_size = int(len(correct_time)/N)
    
    #Making lists that will contain the N segments of time and PI data.
    #We extract the PI data to get the periodogram for each segment.
    Time_segments = []
    PI_segments = []

    #Populating the time and PI segment lists
    for i in range(N):
        Time_segments.append(correct_time[seg_size*i:seg_size*(i+1)])
        PI_segments.append(PI_data[seg_size*i:seg_size*(i+1)])
    return Time_segments, PI_segments


def get_power_and_freq(time_segments, energy_segments, n):
    '''
    Function to get the power and frequency from the Lomb-Scargle Periodogram
    of the event arrival times array and corresponding energy array.
    Parameters
    ----------
    :param time_segments: array_like, event arrival times
    :param energy_segments: array_like, energy associated to event arrival times (in PI)
    :param n: int, number of segments we consider (useful for the pulse profile matrix)
    Returns
    ----------
    :param freq_segments: array, containing frequencies of the Lomb-Scargle Periodogram
    :param Power_segments: array, containing powers associated to each frequency of the Lomb-Scargle Periodogram
    '''
    
    #Defining the lists that will contain frequency and power values
    freq_segments=[]
    Power_segments=[]
    
    #Populating the list
    for i in range(n):
        freq_temp, power_temp = LombScargle(time_segments[i], energy_segments[i]).autopower()
        freq_segments.append(freq_temp)
        Power_segments.append(power_temp)   
    return freq_segments, Power_segments


def Harmonic_funk(order_fit, limit, ref_time, time_segments, freq_segments, power_segments, guess_frequency):
    '''
    Function to get the first harmonic frequencies for each segment.
    We will then fit a polynomial to these frequencies to make an accurate epoch_folding function.
    Parameters
    ----------
    :param order_fit: int, the order of the poynomial fit we want to do
    :param limit: float, the cut-off we apply to remove the red noise from the Periodogram
    :param ref_time: float, a reference time used to calculate the elements in bin list (further detailed below)
    :param time_segments: array_like, event arrival times
    :param freq_segments: array_like, frequency values from the Periodogram computed above
    :param power_segments: array_like, power values from the Periodogram computed above
    Returns
    ----------
    :param reg: list, the coefficients of the polynomial fit to the harmonic frequencies. Can be used to find the 
    derivative of frequency with respect to segments.
    :param bins: list, containing the average time per segment since the start of observations.
    '''
    
    #Making the harmonic frequency list and a bin list that we 
    #will use to plot the frequencies and fit a linear regression. 
    #We fill the bin list with the average time per segment SINCE the 
    #start of observations. This is so that the regression 
    #can be used later for folding our time lists.
    Harmonic=[]
    bins = []
    
    #Populating the Harmonic and bins lists 
    for i in range(len(time_segments)):
        bins.append(float(np.mean(time_segments[i]-ref_time)))
    
    #Cutting off frequencies (and power associated) 
    #that are due to red noise (low frequencies)
    #Obtained by plotting the periodogram and checking when the red noise decreases -> change it?
        temp_1 = freq_segments[i] > limit
    
    #Getting the cut-off frequencies and powers
        temp_2 = freq_segments[i][temp_1]
        temp_3 = power_segments[i][temp_1]
    
    #Getting the harmonic frequency i.e. frequency of max power
        temp_freq_max = temp_2[np.argmax(temp_3)]
        Harmonic.append(float(temp_freq_max))
        
    #Get regression coefficients.
    reg = np.polyfit(bins, Harmonic, order_fit)
    
    #We plot the Harmonic frequencies as a function of the average time
    #of segments with respect to start time of observations.
    plt.plot(bins, Harmonic, '.', label='Measured Harmonic freq.')
    
    #Getting modelled harmonic values
    new_Harmonics=[]
    predict=np.poly1d(reg)
    for i in bins:
        new_Harmonics.append(predict(i))
    
    #Plot modelled harmonics on top of the observed ones.
    plt.plot(bins, new_Harmonics, label='Expected Harmonic freq. (Lin.Reg.)')
    plt.axhline(y=guess_frequency, color='green')
    plt.xlabel('Segment number')
    plt.ylabel('Main Hamornic Frequency (Hz)')
    plt.legend()
    
    #Return the polynomial coefficients and the bins array (for further plotting purposes)
    return reg, bins


def phase_fold(ref_time, time_array, coeffs):
    '''
    Epoch folding function used to plot the phase folded
    pulse profiles.
    Parameters
    ----------
    :param ref_time: float, a reference time used to calculate the phase
    :param time_array: array_like, event arrival times
    :param coeffs: array_like, the polynomial fitting coefficients gotten from Harmonic_funk
    special case, if coeffs is one value put it in a list! Can't give value by itself.
    Returns
    ----------
    :param phasefold_time: list, the epoch folded event arrival times with respect to the reference time
    '''
    #We adapt our epoch folding to take into account the frequencies of harmonics, assumed to 
    #correspond to the "pulse period".
    #We differentiate the case where we do not fit a polynomial and use instead an average value 
    if len(coeffs)==1:
        A = (time_array - ref_time)*coeffs[0]
        return A%1
    #In the other case we get the epoch-folded phase using Taylor expansion
    else:
        A = (time_array - ref_time)*coeffs[-1]
        for i in range(2, len(coeffs)+1):
            A+=(1/mt.factorial(i))*coeffs[-i]*(time_array - ref_time)**i
        phasefold_time = A%1
    return phasefold_time


def chisquared(model, expec, uncertainty):
    '''
    Function to calculate the chi-squared statistic given data points
    with their associated uncertainties and theoretical data points.
    Parameters
    ----------
    :param model: array_like, model data 
    :param expec: array_like, observational data 
    :param uncertainty: array_like, uncertainties on the observational data 
    Returns
    ----------
    :param Xi: float, the (unreduced) chisquared value
    '''

    #Set chi-squared to 0 and progressivelly add on
    Xi=0
    for i in range(len(model)):
        Xi+=((model[i]-expec[i])/uncertainty[i])**2
    return Xi


def model_pulse(order, time, counts):
    '''
    Function to model the pulse period
    with an n-order sinusoidal function.
    Parameters
    ----------
    :param order: int, order of the sinusoidal function
    :param time: array_like, epoch-folded lightcurve times used to compute sinusoidal function
    :param counts: array_like, lightcurve counts we want to model using the n-order sinusoidal function
    Returns
    ----------
    :param A: list, the theoretical counts obtained from the n-order sinusoidal function
    :param phases: list, all the phases obtained from the Fourier transform our the lightcurve counts
    '''
    
    #Fourier transform the pulse profiles
    counts_fft = np.fft.rfft(counts)
    
    #Get the phase and amplitudes from Fourier transformed pulse profiles
    phases = np.arctan2(counts_fft.imag, counts_fft.real)
    amplitudes = np.sqrt(counts_fft.imag**2 + counts_fft.real**2)
    
    
    #Getting the sinusoidal model
    base_phi = 2*np.pi*time
    A = 0.5*counts_fft.real[0]
    for i in range(order):
        A+=amplitudes[i+1]*np.cos((i+1)*base_phi + phases[i+1])
        
    #Normalizing it
    A = A/len(counts_fft)
    
    #Return model and all the phases (useful for later)
    return A, phases


def pulse_profile_matrix(segments, ref_time, reg_coeffs, title_str, cutoff, dt=0.01, 
                         confi_level=0.1, error=True, save_img=False, save_img_path='~/'):
    '''
    Function to make and display the pulse profile matrix
    Parameters
    ----------
    :param segments: array_like, even arrival time segments used to make the lightcurve.
    :param ref_time: float, reference time used in the epoch folding function.
    :param reg_coeffs: array_like, coefficients used in the epoch folding function.
    obtained from the regression in Harmonic_funk.
    :param title_str: string, to be used as the main title for the figure.
    :param cutoff: int, the number of segments we want to use.
    :param dt: float, value used in the making of a Lightcurve with Stingray.
    :param confi_level: float, value used in find_optimal function as the cutoff for the
    survival function.
    :param error: bool, True if we want to plot the error bars on each pulse profile.
    False if not.
    :param save_img: bool, True if we want to save the pulse profile matrix to a pdf.
    :param save_img_path: string, corresponding to the path of the pdf we want to save 
    the pulse profile matrix. Only used if save_img=True.
    Returns
    ----------
    orders: list, containing the orders of the sinusoidal functions used
        to model each segment's pulse profile
    phases: nested lists, containing all the phases for each segment's pulse profile
    folded_counts: nested lists, containing all the lightcurve counts for each segment's pulse profile
    folded_phases: nested lists, containing all the epoch-folded
        lightcurve times for each segment's pulse profile
    '''
    
    #Making a function to calculate the order of sinusoidal function to use
    #for the model to overplot on pulse profile of segments. We do this using chisquared
    #TO DO update documentation

    def find_optimal_order(counts, time, conf_level = 0.1):
        '''
        Function to calculate the optimal order of sinusoidal function to use
        for a given pulse profile. The optimal function is found with comparison
        of chisquared values.
        Parameters
        ----------
        :param counts: array_like, epoch-folded lightcurve times used to compute sinusoidal function
        :param time: array_like, lightcurve counts we want to model using the n-order sinusoidal function
        Returns
        :param conf_level: float, confidence level used for the survival function test.
        ----------
        :param n: int, the order of the sinusoidal function to use to model the given pulse profile.
        '''
        n=0
        SFs=[]
        std = [np.sqrt(j) for j in counts]
        #Go through possible values of order and compute the chisquared for each of them
        #Stop once the chisquared has reached <0.001.
        for i in range(1,int(len(counts)/2)):
            model = model_pulse(i, time, counts)[0]
            #We use the chisquared calculator defined above
            #We take into account the degrees of freedom when calculating the
            #survival function with chi2.sf -> chi2.sf gives us our goodness of fit value
            free_deg = len(counts)-(1 + 2*i)
            chisq=chisquared(model, counts, std)
            #devising a failsafe in case we get n=0 (need to re-work this)
            SFs.append(chi2.sf(chisq, free_deg))
            #print('for degree ', i, ' we have SF ', chi2.sf(chisq, free_deg))
            n=i
            if chi2.sf(chisq, free_deg) > conf_level:
                break
        #print(np.diff(SFs))
        if n==int(len(counts)/2)-1 or n>10:
            for i in range(len(np.diff(SFs))):
                if np.diff(SFs)[i]<0:
                    n=i
                    break
            print('Warning: Did not find an order that satisfies the confidence level, '
                  'failsafe system implements:', n, 'with SF=', SFs[n])
        if n==1:
            n=2        
        return n
    
    #Create array that will contain the phase value of 
    #1st harmonic for each segment 
    phases = []
    
    #Create figure variables in the case where we make the normal pulse profile matrix
    if cutoff <= 4:
        fig, axs = plt.subplots(1, 4, figsize=(10,10))
    else:
        fig, axs = plt.subplots(mt.ceil(cutoff/4), 4, sharex=True, figsize=(15, 15))
    fig.subplots_adjust(hspace=.05, wspace=.1)
    axs = axs.ravel()
            
    #Making a list that will contain the orders of sinusoidal model used for each segment 
    orders=[]
    
    #Making lists to store the folded counts and phases for bootstrapping method 
    folded_counts = []
    folded_phases = []
    
    #Making the pulse profile for each segment 
    #and overplotting a first-order sinusoidal model
    #using the phase of first harmonics of each segment.
    for i in range(cutoff):
        #Computing the phase folded time
        test_phase_fold = phase_fold(ref_time, segments[i], reg_coeffs)
    #Making the lightcurve using the phase folded time
        test_lc = Lightcurve.make_lightcurve(test_phase_fold, dt)
        folded_counts.append(test_lc.counts)
        folded_phases.append(test_lc.time)
            
    #Finding order of model to use
        order=find_optimal_order(test_lc.counts, test_lc.time, confi_level)
        orders.append(order)
        print(order)
            
    #Create n-order sinusoidal model using the phase of first harmonics
        model, model_phase = model_pulse(order, test_lc.time, test_lc.counts)

    #Make a list containing the phases for future use 
        phases.append(model_phase)

    #Plot the pulse profile and the cosine obtained from phase of first "order" harmonic frequencies
        if error:
            axs[i].errorbar(test_lc.time, test_lc.counts, yerr = np.sqrt(test_lc.counts), fmt='.')
        else:
            axs[i].plot(test_lc.time, test_lc.counts, '.')
        axs[i].plot(test_lc.time, model)
        if i%4 != 0:
            axs[i].set_yticklabels([])
        fig.suptitle(title_str)
        if save_img:
            plt.savefig(save_img_path+'PulseProfileMatrix.pdf')
    
    #Returning the order of sinusoidal model used, all the phases for each segment, 
    #the folded counts and times for bootstrapping'''
    return orders, phases, folded_counts, folded_phases


def Plot_phases(order, xdata, phase_list):
    '''
    Function to plot the phases for each segment. We differentiate 
    the different order harmonics by different colors. 
    Parameters
    ----------
    :param order: int, order of the sinusoidal model used corresponding to the number of harmonics used
    :param xdata: array_like, data used for plotting. By convention we use the average time 
    per segment with respect to the start of observations. See output of Harmonic_funk.
    :param phase_list: array_like, list of all the phases for every segment's pulse profile. 
    See output of pulse_profile_matrix.
    Returns
    ----------
    No returns. Outputs a plot of the phases of harmonics used by each segment,
    distinguishing the different order of harmonics used.
    '''
    
    #Making lists containing the phases we only care about i.e. the ones used in the sinusoidal model
    relevant_phases=[]
    temp=[]
    
    #We make a small double nested for loop to extract the phases
    #we actually used to make our model
    
    for i in range(order):
        for j in range(len(phase_list)):
        #Don't forget to increment by one the phase we extract
        #because the first value in each list of phases is always 0
            temp.append(phase_list[j][i+1])
        relevant_phases.append(temp)
        temp=[]
    
    #Make n plots for the n phases used in our n-order sinusoidal model 
    for i in range(order):
        plt.plot(xdata, relevant_phases[i], '.', label=str(i)+'th phase')
        plt.ylabel('Phase')
        plt.xlabel('Average time per segment')
        plt.title('Plot of the first '+str(order)+' harmonics\' phases')
        plt.legend()
    plt.show()
    return 


def segment_timewise(time, n):
    '''
    Function to segment the event arrival times chronologically into n segments.
    Parameters
    ----------
    :param time: array_like, containing event arrival times (with orbital correction).
    :param nbin: int, number of segments we want to make.
    Returns
    time_segments: nested lists, containing the chronological times segments.
    ----------
    '''
    #Getting the size of the time array
    size = len(time)
    
    #Getting the size of each segment
    seg_size = int(size/n)
    
    #Creating a list that will contain the time array segments
    time_segments = []
    
    #Populating time_segments
    for i in range(seg_size):
        time_segment = np.array(time[i*seg_size:(i+1)*seg_size])
        time_segments.append(time_segment)

    return time_segments


def segment_energywise(time_data, energy_data, nmin, nmax, nbin, SNR=False, SNR_freq=0, SNR_dt=0.01, target_SNR=50):
    '''
    Function to segment the event arrival times depending on their energy levels, 
    considering events between a certain range of energies.
    Parameters
    ----------
    :param time_data: array_like, containing event arrival times (with orbital correction).
    :param energy_data: array_like, containing event associated to event arrival times.
    :param nmin: float, minimum energy to consider.
    :param nmax: float, maximum energy to consider.
    :param nbin: int, number of segments we want to make. If nbin<0 we use a logarithmic segmentation.
    If nbin>0 we use a linear segmentation.
    :param SNR: bool, telling us whether we want an adaptive binning depending on the SNR or not.
    :param SNR_freq: float, frequency to use for the phase folding in the SNR calculator.   
    :param SNR_dt: float, value used in the making of a Lightcurve with Stingray. This is used in the
    SNR calculator.
    :param target_SNR: float, SNR value we want each bin to have if we use an SNR adaptive binning.
    Returns
    bins_time: nested lists in the SNR case, containing the times associated to the energies of the segments.
    bins: nest lists in the SNR case, containing the energies for the segments.
    energy_time_segments: nested lists, containing the times associated to
    the energies for nbin segments.
    energy_segments: nested lists, containing the energies for nbin segments in increasing energy order.
    ----------
    '''
    
    def SNR_calculator(time_data, energy_data, freq, dt):
        '''
        Function to segment the event arrival times depending on their energy levels, 
        considering events between a certain range of energies.
        Parameters
        ----------
        :param time_data: array_like, containing event arrival times (with orbital correction).
        :param energy_data: array_like, containing event associated to event arrival times.
        :param dt: float, value used in the making of a Lightcurve with Stingray.
        :param freq: float, frequency to use for the phase folding.
        Returns
        SNR: Signal-to-Noise ratio.
        ----------
        '''
        mphase_fold = phase_fold(time_data[0], time_data, freq)
        lc = Lightcurve.make_lightcurve(mphase_fold, dt)
        S = np.sum(lc.counts)
        N = np.sqrt(S)
        return S/N
    
    #Differentiating the cases depending on sign of nbin.
    #For both cases we create a list of boundaries we will use to create the energy
    #segments.
    if nbin<0:
        bins = 10**(np.linspace(np.log10(nmin), np.log10(nmax), np.abs(nbin)+1))
    else:
        bins = np.linspace(nmin, nmax, nbin+1)

    if not SNR:
    #Initiating time and energy segments' lists to populate them in the next step.
        energy_time_segments = []
        energy_time_segment = []
        energy_segments = []
        energy_segment = []
        if nbin<0:
            for i in range(np.abs(nbin)):
                for j in range(len(energy_data)):
                    if bins[i] < energy_data[j] < bins[i+1]:
                        energy_segment.append(energy_data[j]*4e-2)
                        energy_time_segment.append(time_data[j])
                energy_segments.append(np.array(energy_segment))
                energy_time_segments.append(np.array(energy_time_segment))
                energy_time_segment=[]
                energy_segment=[]
        else: 
            for i in range(nbin):
                for j in range(len(energy_data)):
                    if bins[i] < energy_data[j] < bins[i+1]:
                        energy_segment.append(energy_data[j]*4e-2)
                        energy_time_segment.append(time_data[j])
                energy_segments.append(np.array(energy_segment))
                energy_time_segments.append(np.array(energy_time_segment))
                energy_time_segment=[]
                energy_segment=[]
           
    if SNR:
        bins_energy = []
        bins_time = []
        SNRs=[]
        important_ind=[0]
        for i in range(10, len(time_data), 500):
            SNRs.append(SNR_calculator(time_data[:i], energy_data[:i], SNR_freq, SNR_dt))
            print(i, SNR_calculator(time_data[:i], energy_data[:i], SNR_freq, SNR_dt))
            
        for i in range(len(SNRs)):
            if SNRs[i]%50<1 or -SNRs[i]%50<1:
                print(SNRs[i])
                important_ind.append(10+500*i)
                print('The index in the actual array is:', 10+500*i)
        for i in range(len(important_ind)):
            bins_energy.append(energy_data[important_ind[i]:important_ind[i+1]])
            bins_time.append(time_data[important_ind[i]:important_ind[i+1]])      
    if SNR:
        return bins_time, bins
    else:
        return energy_time_segments, energy_segments

def bootstrap_generate(counts, k):
    '''
    Function to make k realizations of a given pulse profile using a Poisson distribution
    to describe the number of counts.
    Parameters
    ----------
    :param counts: array_like containing the lightcurve counts for a given pulse profile.
    :param k: int, the number of realizations of the pulse profile to make
    Returns
    ----------
    :param fake_profiles: nested lists, the simulated k realizations of the given pulse profile. 
    Can be used for plotting to make sure output of function is plausible.
    '''
    #Setting up lists for the fake pulse profile values
    #and matrix multiplication step
    fake_profiles = []
    fake_profile = []
    row=[]
    
    #Making k realizations of our pulse profile using a Poisson distibution centered
    #on the number of photon counts of each point in the pulse profile.
    for i in range(len(counts)):
        fake_profile.append(st.poisson.rvs(mu=counts[i], size=k))
    
    #CF : here maybe you could just transpose the matrix
    #Re-making the k realizations so that instead of having an n x k matrix 
    #We have a k x n matrix -> will make it easier to plot the k realizations of the 
    #pulse profile.
    for j in range(k):
        for i in range(len(fake_profile)):
            row.append(fake_profile[i][j])
        fake_profiles.append(row)
        row=[]
    return fake_profiles


def bootstrap_total(counts, k, func, *args):
    '''
    Function to calculate the standard deviation of a specific quantity obtained from
    the function func.
    Parameters
    ----------
    :param counts: array_like containing the lightcurve counts for a given pulse profile
    :param k: int, the number of realizations of the pulse profile to make
    :param func: function, to calculate a certain quantity/statistic of the pulse profile
    :param *args: any additional arguments required for func
    Returns
    ----------
    :param std: float, standard deviation on the quantity/statistic for the given given pulse profile.
    '''
    #Get the k realizations of the pulse profile from the bootstrap_generate function
    profiles = bootstrap_generate(counts, k)
    #Prepare an array that will contain the quantity we want for every realization
    values = []
    #Populate the values list with the quantity using the input func function and its arguments
    for i in range(len(profiles)):
        values.append(func(profiles[i], *args))
    
    #Calculating the standard deviation of the quantity of interest using the k realizations of it
    std = np.sqrt(np.sum((np.array(values)-np.mean(values))**2)/(k-1.0))
    
    #Compared to the plot in your Elaboration.ipynb, the RMS errors bars in my RMSvsEnergy plot 
    # at the end of Day3.ipynb are smaller.
    return std

def get_first_harmonic_phase(counts):
    '''
    Function to calculate the phase of the first harmonic frequency for a given pulse profile.
    Parameters
    ----------
    :param counts: array_like, containing the lightcurve counts for a given pulse profile
    Returns
    ----------
    :param phase_l: float, phase of the first harmonic frequency for the given pulse profile.
    '''
    #Calculating the phases for the lightcurve counts using the Fourier transform
    #and taking the first element of the list of phases (technically the second element
    #because the very first phase is always 0).
    fft = np.fft.rfft(counts)
    phase_l = np.arctan2(fft.imag, fft.real)[1]
    return phase_l

def RMS_calculator(counts, k):
    '''
    Function to calculate the Root-Mean Squared (RMS) for a given pulse profile.
    Parameters
    ----------
    :param counts: array_like, containing the lightcurve counts for a given pulse profile
    :param k: int, number of amplitudes/phases used for the sinusoidal model of the pulse
    profile we consider. Corresponds to the order of the model.
    Returns
    ----------
    :param RMS: float, RMS for the given pulse profile.
    '''
    #We calculate the amplitudes associated to each harmonics' phase using a 
    #Fourier transform of the pulse profile. We then perform a RMS method on these
    #amplitudes and output the result.
    counts_fft = np.fft.rfft(counts)
    
    #Hi Carlo, the error I was getting with the RMS was because I was using these equations: 
    #amps = np.sqrt(counts_fft.imag**2 + counts_fft.real**2)[1:k]
    #RMS = np.sqrt(np.sum(amps**2))/amps[0]
    
    #Instead of these equations 
    amps = np.sqrt(counts_fft.imag**2 + counts_fft.real**2)
    RMS = np.sqrt(np.sum(amps[1:k]**2))/amps[0]
    
    #As you can see this is simply an indexing error but I am still a little confused as to
    #why it is happening.
    
    #CF you were summing all harmonics, before, not just the first "k"
    return RMS


#Note: the below function is used with Renkulab and as a result works with model 
#instances from xspec and pyxmmas!
def BIC_calc(chi2, n_spectra):
    '''
    Function to calculate the Bayesian Imformation Criterion (BIC)
    given the chi-squared for our model. 
    Parameters
    ----------
    :param chi2: float, chi-squared value to calcule the BIC for.
    Returns
    ----------
    BIC: float, BIC value.
    '''
 
    #Initialize the number of data points used
    n=0
    #Initialize the number of free parameters used 
    k=0
    
    #Summing the number of data points of each 
    #individual spectra
    for i in range(1, n_spectra+1):
        n += len(xspec.AllData(i).noticed)

    #Get the models for each spectra
    models = [xspec.AllModels(j) for j in range(1, n_spectra + 1)]

    #Make a function to get the number of free parameters
    def get_frozen(nn):
        '''
        Function to calculate the number of free parameters. 
        Parameters
        ----------
        :param nn: object, model for a spectra
        Returns
        ----------
        k: int, number of free parameters
        '''
        #Initializing k
        k=0
        #Iterating over the components of the model
        for i in nn.componentNames:
            #Iterating over the parameters of each component
            for j in getattr(nn,i).parameterNames:
                #If the parameters are frozen, the output is True
                #So we add to the sum when the output is False
                if not getattr(getattr(nn,i), j).frozen:
                    k+=1
        return k

    #Getting the total number of free parameters for all the spectra used
    for i in range(n_spectra):
        k+=get_frozen(models[i])
        
    #Return the BIC
    return chi2 + k*np.log(n)