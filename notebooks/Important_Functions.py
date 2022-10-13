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


def Harmonic_funk(order_fit, limit, ref_time, time_segments, freq_segments, power_segments):
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
    for i in range(20):
        bins.append(float(np.mean(time_segments[i]-ref_time)))
    
    #Cutting off frequencies (and power associated) 
    #that are due to red noise (low frequencies)
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
    plt.plot(bins, new_Harmonics, '.', label='Expected Harmonic freq. (Lin.Reg.)')
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


def pulse_profile_matrix(segments, ref_time, reg_coeffs, special_case=False, special_case_segments=2, dt=0.01):
    '''
    Function to make and display the pulse profile matrix
    Parameters
    ----------
    :param segments: int, time segments to use
    :param ref_time: float, reference time used in the epoch folding function
    :param reg_coeffs: array_like, coefficients used in the epoch folding function 
    obtained from the regression in Harmonic_funk
    :param special_case: boolean, to say whether we are getting a "normal" 
    pulse profile matrix or a special one that
    :param special_case_segments: int, value pertaining to the number of segments we have in the special case
    :param dt: float, value used in the making of a Lightcurve with Stingray
    Returns
    ----------
    :param orders: list, containing the orders of the sinusoidal functions used 
    to model each segment's pulse profile
    :param phases: nested lists, containing all the phases for each segment's pulse profile
    :param folded_counts: nested lists, containing all the lightcurve counts for each segment's pulse profile
    :param folded_phases: nested lists, containing all the epoch-folded 
    lightcurve times for each segment's pulse profile
    '''
    
    #Making a function to calculate the order of sinusoidal function to use
    #for the model to overplot on pulse profile of segments. We do this using chisquared
    def find_optimal_order(counts, time):
        '''
        Function to calculate the optimal order of sinusoidal function to use
        for a given pulse profile. The optimal function is found with comparison
        of chisquared values.
        Parameters
        ----------
        :param counts: array_like, epoch-folded lightcurve times used to compute sinusoidal function
        :param time: array_like, lightcurve counts we want to model using the n-order sinusoidal function
        Returns
        ----------
        :param n: int, the order of the sinusoidal function to use to model the given pulse profile.
        '''
        n=0
        failsafe=[]
        std = [np.sqrt(j) for j in counts]
        #Go through possible values of order and compute the chisquared for each of them
        #Stop once the chisquared has reached <0.001.
        for i in range(1,int(len(counts)/2)):
            model = model_pulse(i, time, counts)[0]
            #We use the chisquared calculator defined above
            #We take into account the degrees of freedom when calculating the
            #survival function with chi2.sf -> chi2.sf gives us our goodness of fit value
            free_deg = 1 + 2*i
            chisq=chisquared(model, counts, std)
            #devising a failsafe in case we get n=0 (need to re-work this)
            failsafe.append(chi2.sf(chisq, free_deg))
            if chi2.sf(chisq, free_deg) > .1:
                n=i
                break
        if n==0:
            n=failsafe[2]
        return n
    
    #Create array that will contain the phase value of 
    #1st harmonic for each segment 
    phases = []

    #Create figure variables in the case where we make the normal pulse profile matrix
    if not special_case:
        f, ((ax1, ax2, ax3, ax4, ax5), 
    (ax6, ax7, ax8, ax9, ax10), 
    (ax11, ax12, ax13, ax14, ax15), 
    (ax16, ax17, ax18, ax19, ax20)) = plt.subplots(4, 5, sharex=True, figsize=(10, 10), sharey=True)
        axs = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,
       ax11,ax12,ax13,ax14,ax15,ax16,ax17,ax18,ax19,ax20] 
    
    #Making a list that will contain the orders of sinusoidal model used for each segment 
    orders=[]
    
    #Making lists to store the folded counts and phases for bootstrapping method 
    folded_counts = []
    folded_phases = []
    
    #Making the pulse profile for each segment 
    #and overplotting a first-order sinusoidal model
    #using the phase of first harmonics of each segment.
    if not special_case:
        for i in range(20):
        
        #Computing the phase folded time
            test_phase_fold = phase_fold(ref_time, segments[i], reg_coeffs)

    #Making the lightcurve using the phase folded time
            test_lc = Lightcurve.make_lightcurve(test_phase_fold, dt)
            folded_counts.append(test_lc.counts)
            folded_phases.append(test_lc.time)
            
    #Finding order of model to use
            order=find_optimal_order(test_lc.counts, test_lc.time)
            orders.append(order)
            print(order)
            
    #Create n-order sinusoidal model using the phase of first harmonics
            model, model_phase = model_pulse(order, test_lc.time, test_lc.counts)

    #Make a list containing the phases for future use 
            phases.append(model_phase)

    #Plot the pulse profile and the cosine obtained from phase of first "order" harmonic frequencies
            axs[i].plot(test_lc.time, test_lc.counts, '.')
            axs[i].plot(test_lc.time, model, '.')
        axs[2].set_title('Pulse Profile Matrix')
        axs[5].set_ylabel('Photon count')
        axs[17].set_xlabel('Phase')
    else:
        for i in range(special_case_segments):
        
        #Computing the phase folded time
            test_phase_fold = phase_fold(ref_time, segments[i], reg_coeffs)

    #Making the lightcurve using the phase folded time
            test_lc = Lightcurve.make_lightcurve(test_phase_fold, dt)
            folded_counts.append(test_lc.counts)
            folded_phases.append(test_lc.time)
    #Finding order of model to use
    
            order=find_optimal_order(test_lc.counts, test_lc.time)
            orders.append(order)
            print(order)
    #Create n-order sinusoidal model using the phase of first harmonics
            model, model_phase = model_pulse(order, test_lc.time, test_lc.counts)

    #Make a list containing the phases for future use 
            phases.append(model_phase)

    #Plot the pulse profile and the cosine obtained from phase of first harmonic frequency
            plt.plot(test_lc.time, test_lc.counts, '.')
            plt.plot(test_lc.time, model, '.')
            plt.show()
    #Returning the order of sinusoidal model used, all the phases for each segment, 
    #the folded counts and times for bootstrapping
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


def segment_energywise(time_data, energy_data, nmin, nmax, nbin):
    '''
    Function to segment the event arrival times depending on their energy levels, 
    considering events between a certain range of energies.
    Parameters
    ----------
    :param time_data: array_like, containing event arrival times.
    :param energy_data: array_like, containing event associated to event arrival times.
    :param nmin: float, minimum energy to consider.
    :param nmax: float, maximum energy to consider.
    :param nbin: int, number of segments we want to make. If nbin<0 we use a logarithmic segmentation.
    If nbin>0 we use a linear segmentation.
    Returns
    :param energy_time_segments: nested lists, containing the times associated to 
    the energies for nbin segments.
    :param energy_segments: nested lists, containing the energies for nbin segments in increasing energy order.
    ----------
    '''
    if nbin < 0:
        size_bin = mt.floor(np.log(nmax-nmin)/nbin)
    else:
        size_bin = mt.floor((nmax-nmin)/nbin)
    
    #Defining lists we want to populate with the event arrival times 
    #and energy segments
    energy_time_segments = []
    energy_time_segment = []
    energy_segments = []
    energy_segment = []
    
    #Distinguishing the cases and populating
    if nbin<0:
        for i in range(nbin):
            for j in range(len(energy_data)):
                if np.exp(size_bin)*i*40 < Time_phase_file[1].data['PI'][j] < np.exp(size_bin)*(i+1)*40:
                    energy_segment.append(Time_phase_file[1].data['PI'][j]*40)
                    energy_time_segment.append(correct_orbit_time[j])
            energy_segments.append(energy_segment)
            energy_time_segments.append(energy_time_segment)
            energy_time_segment=[]
            energy_segment=[]
    else: 
        for i in range(nbin):
            for j in range(len(energy_data)):
                if size_bin*i < Time_phase_file[1].data['PI'][j] < size_bin*(i+1):
                    energy_segment.append(Time_phase_file[1].data['PI'][j]*40)
                    energy_time_segment.append(correct_orbit_time[j])
            energy_segments.append(energy_segment)
            energy_time_segments.append(energy_time_segment)
            energy_time_segment=[]
            energy_segment=[]
            
    #Returning the event arrival times in terms of their energy
    return energy_time_segments, energy_segments


def bootstrap(pulse_profile, k):
    '''
    Function to get the uncertainty of the first harmonic's phase for a given pulse profile.
    This is done by making k realizations of the given pulse profile using a Poisson distribution
    to describe the number of counts.
    Parameters
    ----------
    :param pulse_profile: tuple, containing the epoch-folded lightcurve times and lightcurve counts for
    a given pulse profile.
    :param k: int, the number of realizations of the pulse profile to make
    Returns
    ----------
    :param fake_profiles: nested lists, the simulated k realizations of the given pulse profile. 
    Can be used for plotting to make sure output of function is plausible.
    :param std: float, standard deviation on the phase of first harmonic for a given pulse profile.
    :param true_phase: float, the phase of first harmonic for the given pulse profile we consider.
    '''
    #Getting the time and counts data
    time, counts = pulse_profile
    
    #Setting up lists for the fake pulse profile values
    #and matrix multiplication step as well as the k realizations 
    #of the first harmonic phase
    fake_profiles = []
    fake_profile = []
    row=[]
    phases_k=[]
        
    #Getting the actual phase of first harmonic for the pulse profile we are considering
    #to make sure our results are sensible 
    true_phase = np.arctan2(np.fft.rfft(counts).imag, np.fft.rfft(counts).real)[1]
        
    #Making k realizations of our pulse profile using a Poisson distibution centered
    #on the number of photon counts of each point in the pulse profile.
    for i in range(len(counts)):
        fake_profile.append(st.poisson.rvs(mu=counts[i], size=k))
        
    #Re-making the k realizations so that instead of having an n x k matrix 
    #We have a k x n matrix -> will make it easier to plot the k realizations of the 
    #pulse profile.
    for j in range(k):
        for i in range(len(fake_profile)):
            row.append(fake_profile[i][j])
        fake_profiles.append(row)
        row=[]
    
    #Calculating the phase of first harmonic for the k realizations 
    #and populating the phases_k list. 
    for l in range(len(fake_profiles)):
        fft = np.fft.rfft(fake_profiles[l])
        phase_l = np.arctan2(fft.imag, fft.real)
        phases_k.append(phase_l[1])
        
    #Getting the standard deviation of our phase measurement using the 
    #bootstrapping method)
    std = (1/k)*np.sum((np.array(phases_k)-np.mean(phases_k))**2)
    
     
    return fake_profiles, std, true_phase


