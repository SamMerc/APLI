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
import lmfit as lf

def get_GTIs(data_file):
    '''
    Function to get the Good Time Intervals (GTIs).
    Parameters
    ----------
    :param data_file: file name, containing the source's GTI data file.
    Returns
    ----------
    :param new_gtis: tuples list, containing tuples of start and stop time of the GTIs.
    '''
    #Retrieving the GTIs.
    GTI_phase_file = data_file
    #Re-formatting the GTIs so that we can use them in the Stingray package.
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
    other useful data (start and end times).
    :param param_file_name: string, containing the name for the parameter file used
    for the orbital correction.
    Returns
    ----------
    :param correct_orbit_time: array, containing event arrival times with orbital correction.
    '''
    #Reference time in MJD.
    base_time = data_file.header['MJDREFI']+data_file.header['MJDREFF']
    
    #Converting start time of observations in seconds.
    start_time = data_file.header['TSTART']/(24*3600)
    
    #Converting end time of observations in seconds.
    end_time = data_file.header['TSTOP']/(24*3600)
    
    #Getting the correct_time function from the parameters defined above 
    #and the file containing the orbital parameters of the source.
    correct_time = get_orbital_correction_from_ephemeris_file(base_time+start_time, base_time+end_time, parfile=param_file_name)[0]
    
    #Applying the correct_time function to our event arrival times.
    correct_orbit_time = correct_time(data_file.data['TIME'], base_time)
    
    return correct_orbit_time

def get_pulse_freq(N, correct_time_array, use_single, time_array, gtis, distance, seg_size, cutoff, dt=0.05):
    '''
    Function to get the pulse frequency using the event arrival times with orbital correction.
    Parameters
    ----------
    :param N: int, the number of frequencies to test between f_min and f_max.
    :param correct_time_array: array, containing the source's concatenated event arrival times 
    with orbital correction from both detectors.
    :param use_single: bool, whether to use the concatenated data from both detectors (False) 
    or data from a single detector (True). 
    This is particularly useful for large datasets to reduce periodogram computation time.
    :param time_array: array, containing the source's event arrival times 
     with orbital correction from a single detector.
    :param gtis: array, containing the GTI pairs.
    :param distance: float, used to define the range of frequencies around the periodogram best 
    frequency we want to probe with the epoch folding search algorithm.
    :param seg_size: int, length of segments to to be averaged in periodogram.
    :param cutoff: float, frequency delimiting the low frequency tail of the periodogram for
    which the power is too high. We are interested in the high frequency regime where the harmonics
    can be found.
    :param dt: float, value used in the making of a Lightcurve with Stingray.
    Returns
    ----------
    :param correct_L: tuple list, containing a list of frequencies and a list of corresponding power.
    :param best_freq: float, pulse frequency obtained from the periodogram and epoch folding
    search algorithms. 
    :param periodogram_freq: float, pulse frequency obtained from the periodogram. 
    '''
    
    #Differentiating the cases where we use data from a single detector or both
    #detectors.
    if use_single:
        #Separating the data even further, by using every n data points, where
        #n is a random number between 20 and 50.
        
        #Creating a new array of times.
        new_times = [time_array[0]]
        
        #Creating an index to figure out where we are in the original data array.
        index=0
        
        #Iterating over the original data array.
        while index<len(time_array)-50:
            #Choosing a random number and a corresponding datum.
            x=np.random.randint(20, 50)
            index+=x
            #Adding the data point to the new data array.
            new_times.append(time_array[index])

        #Making a lightcurve from the new data array we just created.
        new_l = Lightcurve.make_lightcurve(new_times, dt, gti=gtis)
        new_l.apply_gtis()
        
        #Getting a periodogram from the lightcurve.
        print('Starting periodogram')
        f_temp, power_temp = LombScargle(new_l.time, new_l.counts).autopower()
    else:
        #Separating the data even further, by using every n data points, where
        #n is a random number between 20 and 50.
        
        #Creating a new array of times.
        new_times = [correct_time_array[0]]
        
        #Creating an index to figure out where we are in the original data array.
        index=0
        
        #Iterating over the original data array.
        while index<len(correct_time_array):
            #Choosing a random number and a corresponding datum.
            x=np.random.randint(1, 10)
            index+=x
            #Adding the data point to the new data array.
            new_times.append(correct_time_array[index])
        
        #Making a lightcurve from the new data array we just created.
        new_l = Lightcurve.make_lightcurve(correct_time_array, dt=0.01, gti=gtis)
        new_l.apply_gtis()
        
        #Getting a periodogram from the lightcurve.
        print('Starting periodogram')
        f_temp, power_temp = LombScargle(new_l.time, new_l.counts).autopower()
    
    #Plotting the periodogram (for debugging purposes and sanity checks).
    plt.loglog(f_temp, power_temp)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power')
    plt.title('Periodogram')
    plt.show()
    
    #Considering only the frequencies after the cutoff.
    power_temp = power_temp[f_temp>cutoff]
    f_temp = f_temp[f_temp>cutoff]
    
    #We also cutoff the high frequencies at the Nyquist frequency.
    #High-frequencies are dominated by dead time intervals which we want to avoid.
    
    #We retrieve the frequency corresponding to the maximum power in our cutoff
    #interval of frequencies.
    periodogram_freq = f_temp[np.argmax(power_temp[max(f_temp)/2>f_temp])]
    
    #Creating a range of frequencies we want to probe with the epoch folding search.
    f_min =  periodogram_freq - distance
    f_max = periodogram_freq + distance
    
    #Making an array of frequencies for the epoch folding search algorithm (Stingray implementation).
    trial_freqs = np.linspace(f_min, f_max, N)
    
    print('Finished periodogram step, found frequency:', periodogram_freq)
    print('Now using Epoch-Folding Search around frequencies [', f_min, ',', f_max,']')
    
    #Performing the epoch folding search algorithm.
    #We only perform epoch folding search on the concatenated data as 
    #it is a fast algorithm, and we have already narrowed down the range of
    #frequencies to study with our periodogram step.
    correct_L = epoch_folding_search(correct_time_array, trial_freqs, 
                                     segment_size=seg_size, gti=gtis)
    
    #We retrieve the frequency corresponding to the maximum power in our epoch folding search algorithm.
    best_freq = correct_L[0][np.argmax(correct_L[1])]
    print('Finished epoch folding step, found frequency:', best_freq)
    print('Refinement on the pulse frequency of:', np.abs(periodogram_freq - best_freq), 'Hz')
    return correct_L, best_freq, periodogram_freq

def phase_fold(ref_time, time_array, coeffs):
    '''
    Epoch folding function used to plot the phase folded.
    pulse profiles.
    Parameters
    ----------
    :param ref_time: float, a reference time used to calculate the phase.
    :param time_array: array_like, event arrival times.
    :param coeffs: array_like, polynomial fitting coefficients obtained by fitting the Time-Phase plot.
    special case, if coeffs is one value put it in a list! Can't give value by itself.
    Returns
    ----------
    :param phasefold_time: list, the epoch folded event arrival times with respect to the reference time.
    '''
        
    #Case where we use a constant pulse frequency.  
    if len(coeffs)==1:
        A = (time_array - ref_time)*coeffs[0]
        return A%1
    #Case where we get a pulse frequency derivative from fitting the Time-phase plot.
    else:
        #Initializing the epoch folded data.
        A = (time_array - ref_time)*coeffs[0]
        
        #Adding the terms from a Taylor expansion - detailed in our report.
        for i in range(1, len(coeffs)):
            A+=(1/mt.factorial(i+1))*coeffs[i]*(time_array - ref_time)**(i+1)
        phasefold_time = A%1
    return phasefold_time


def segment_timewise(time, n):
    '''
    Function to segment the event arrival times chronologically into n bins.
    Parameters
    ----------
    :param time: array_like, containing event arrival times, with or without orbital correction.
    :param nbin: int, number of segments we want to make.
    Returns
    ----------
    time_segments: nested lists, containing the chronological times segments.
    '''
    #Getting the size of the time array.
    size = len(time)
    
    #Getting the size of each segment.
    seg_size = int(size/n)
    
    #Creating a list that will contain the time array segments.
    time_segments = np.ones((n, seg_size))
    
    #Populating time_segments.
    for i in range(n):
        time_segments[i] = np.array(time[i*seg_size:(i+1)*seg_size])

    return time_segments


def segment_energywise(time_data, energy_data, nmin, nmax, nbin, threshold, load_bins, concatenation=False, load=False):
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
    :param threshold: int, reference length of the segments that we want. Used for logarithmic binning.
    :param load_bins: array, manually built bins that can be used to segment the data.
    :param concatenation: bool, True if you want to concatenate the segments such that they all have
    a minimum number of points, obtained from threshold parameter.
    :param load: bool, whether to use the loaded bins or not.
    Returns
    ----------
    energy_time_segments: nested lists, containing the times associated to
    the energies for nbin segments.
    energy_segments: nested lists, containing the energies for nbin segments in increasing energy order.
    '''
    
    #Differentiating the cases depending on sign of nbin, and on whether or not 
    #the bins are loaded manually.
    #For both cases we create a list of boundaries we will use to create the energy
    #segments.
    if not load:
        if nbin<0:
            bins = 10**(np.linspace(np.log10(nmin), np.log10(nmax), np.abs(nbin)+1))
        else:
            bins = np.linspace(nmin, nmax, nbin+1)

        #Initializing time and energy segments lists to populate them in the next step.
        
        #The energy_time_segments and energy_segments will contain energy_time_segment and
        #energy_segment lists.
        energy_time_segments = []
        energy_time_segment = []
        energy_segments = []
        energy_segment = []
    
        #If we have logarithmic bins.
        if nbin<0:
            #Loop over the bins.
            for i in range(np.abs(nbin)):
                #Loop over the data.
                for j in range(len(energy_data)):
                    #Check if the data falls in the bins - if so, add to the subsegments.
                    if bins[i] < energy_data[j] < bins[i+1]:
                        energy_segment.append(energy_data[j])
                        energy_time_segment.append(time_data[j])
                #Append the subsegments to our master list.
                energy_segments.append(np.array(energy_segment))
                energy_time_segments.append(np.array(energy_time_segment))
                energy_time_segment=[]
                energy_segment=[]
            
            #We implement a routine to concatenate the lists that contain 
            # a low number of events - this is done to improve the building 
            #of pulse profiles at high-energies as the number of photon counts
            #decreases exponentially in those regimes.
            if concatenation:
                #Initialize a new list of segments and sublists for each segment.
                new_segments = []
                new_segments_t = []
                #Initializing the index to go through the master list of energy segments made prior.
                i=0
                #Iterating over all the segments - if the size of the sub-segments are below a 
                #reference value, we concatenate them together.
                while i<len(energy_segments):
                
                    new_segment = energy_segments[i]
                    new_segment_t = energy_time_segments[i]
                    j=i
            
                    #Concatenating until the size of the subsegment reaches our reference value.
                    while len(new_segment)<threshold and j<len(energy_segments)-1:
                        j+=1
                        new_segment=np.concatenate((new_segment, energy_segments[j]))
                        new_segment_t=np.concatenate((new_segment_t, energy_time_segments[j]))
    
                    i=j+1
                    #Adding our subsegments to the corresponding master list.
                    new_segments.append(new_segment)
                    new_segments_t.append(new_segment_t)

                #The last segment is skipped in our iteration and is generally the one
                #that suffers the most from decrease in photon counts. 
                #To double-check this, we check the size of the last segment. If it is indeed
                #lower than our reference value, we concatenate it to the second-to-last segment.
                if len(new_segments[-1])<threshold:
                    new_segments[-2]=np.concatenate((new_segments[-2], new_segments[-1]))
                    new_segments = new_segments[:-1]

                    new_segments_t[-2]=np.concatenate((new_segments_t[-2], new_segments_t[-1]))
                    new_segments_t = new_segments_t[:-1]
        
            #Updating the master lists.
                energy_time_segments = new_segments_t
                energy_segments = new_segments
        
        #If we have linear bins.
        else: 
            #Loop over the bins.
            for i in range(nbin):
                #Loop ove the data.
                for j in range(len(energy_data)):
                    #Check if the data falls in the bins - if so, add to the subsegments.
                    if bins[i] < energy_data[j] < bins[i+1]:
                        energy_segment.append(energy_data[j])
                        energy_time_segment.append(time_data[j])
                #Append the subsegments to our master list.
                energy_segments.append(np.array(energy_segment))
                energy_time_segments.append(np.array(energy_time_segment))
                energy_time_segment=[]
                energy_segment=[]
                
    #Adding a case if we want to load the bins manually
    if load:
        #Initializing time and energy segments lists to populate them in the next step.
        
        #The energy_time_segments and energy_segments will contain energy_time_segment and
        #energy_segment lists.
        energy_time_segments = []
        energy_time_segment = []
        energy_segments = []
        energy_segment = []

        #Loop over the bins.
        for i in range(len(load_bins)-1):
           #Loop ove the data.
            for j in range(len(energy_data)):
                    #Check if the data falls in the bins - if so, add to the subsegments.
                if load_bins[i] < energy_data[j] < load_bins[i+1]:
                    energy_segment.append(energy_data[j])
                    energy_time_segment.append(time_data[j])
                #Append the subsegments to our master list.
            energy_segments.append(np.array(energy_segment))
            energy_time_segments.append(np.array(energy_time_segment))
            energy_time_segment=[]
            energy_segment=[]
        
     
    return energy_time_segments, energy_segments

    
def chisquared(model, expec, uncertainty):
    '''
    Function to calculate the chi-squared statistic given data points
    with their associated uncertainties and theoretical data points.
    Parameters
    ----------
    :param model: array_like, model data.
    :param expec: array_like, observational data.
    :param uncertainty: array_like, uncertainties on the observational data.
    Returns
    ----------
    :param Xi: float, the (unreduced) chisquared value.
    '''

    #Initialize chi-squared to 0.
    Xi=0
    #Progressively add to chi-squared.
    for i in range(len(model)):
        Xi+=((model[i]-expec[i])/uncertainty[i])**2
    return Xi


def model_pulse(order, time, counts):
    '''
    Function to describe the pulse period
    with an n-order sinusoidal function.
    Parameters
    ----------
    :param order: int, order of the sinusoidal function, corresponding to the number of harmonics used.
    :param time: array_like, epoch-folded lightcurve times used to compute sinusoidal function.
    :param counts: array_like, lightcurve counts we want to model using the n-order sinusoidal function.
    Returns
    ----------
    :param A: list, the theoretical counts obtained from the n-order sinusoidal function.
    :param phases: list, phases of all the harmonics of the lightcurve counts for a given pulse profile.
    '''
    
    #Fourier transform the pulse profile.
    counts_fft = np.fft.rfft(counts)

    #Get the phase and amplitudes from Fourier transformed pulse profile.
    phases = np.arctan2(counts_fft.imag, counts_fft.real)
    amplitudes = np.sqrt(counts_fft.imag**2 + counts_fft.real**2)
    

    #Creating the sinusoidal model.
    base_phi = 2*np.pi*time
    A = 0.5*counts_fft.real[0]
    for i in range(order):
        A+=amplitudes[i+1]*np.cos((i+1)*base_phi + phases[i+1])
        
    #Normalizing the model.
    A = A/len(counts_fft)
    
    return A, phases


def pulse_profile_matrix(segments, ref_time, reg_coeffs, title_str, dt=0.01, 
                         confi_level=0.1, error=False, save_img=False, save_img_path='~/'):
    '''
    Function to make and display the pulse profile matrix.
    Parameters
    ----------
    :param segments: array_like, even arrival time segments from which we will make the
    pulse profiles in the matrix.
    :param ref_time: float, reference time used in the phase folding function.
    :param reg_coeffs: array_like, coefficients used in the phase folding function.
    :param title_str: string, used as the main title for the figure.
    :param dt: float, bin size for the Stingray light curve.
    :param confi_level: float, value used in find_optimal function as the cutoff for the
    survival function.
    :param error: bool, True if we want to plot the error bars on each pulse profile. False if not.
    :param save_img: bool, True if we want to save the pulse profile matrix to a pdf.
    :param save_img_path: string, corresponding to the path of the pdf where we want to save 
    the pulse profile matrix. Only used if save_img=True.
    Returns
    ----------
    orders: list, containing the orders of the sinusoidal functions used to model 
    each segment's pulse profile.
    phases: nested lists, containing the phases of harmonics for each pulse profile, 
    obtained from the model_pulse function. 
    folded_counts: nested lists, containing the lightcurve counts for each pulse profile.
    folded_times: nested lists, containing the epoch-folded lightcurve times for each pulse profile.
    '''
    
    def find_optimal_order(counts, time, conf_level = 0.1):
        '''
        Function to calculate the optimal order of sinusoidal function to use
        for a given pulse profile. The optimal function is found with comparison
        of chisquared values.
        Parameters
        ----------
        :param counts: array_like, epoch-folded lightcurve times used to compute sinusoidal function.
        :param time: array_like, lightcurve counts we want to model using the n-order sinusoidal function.
        Returns
        :param conf_level: float, confidence level used for the survival function test.
        ----------
        :param n: int, the order of the sinusoidal function used to describe the given pulse profile.
        '''
        #Initializing the order and the values of the survival function
        #We make a list of survival function values for our failsafe system.
        n=0
        SFs=[]
        
        #Getting the error values on the counts - used for calculating chi-squared.
        std = [np.sqrt(j) for j in counts]
        
        #Go through possible values of order and compute the chisquared for each of them
        #Stop once the survival function has reached a specified confidence level.
        for i in range(1,int(len(counts)/2)):
            
            #Obtain the model counts.
            model = model_pulse(i, time, counts)[0]
            
            #We use the chisquared calculator defined above
            #We take into account the degrees of freedom when calculating the
            #survival function. 
            
            #The number of degrees of freedom is the number of points minus the 
            #number of parameters in the model. Each iteration of the model has a constant term for
            #the amplitude of the 0-th harmonics as well as the phase and amplitude for each harmonic.
            #This is why we use (1+2i).
            free_deg = len(counts)-(1 + 2*i)
            chisq=chisquared(model, counts, std)
            
            #We update our list of survival function values.
            SFs.append(chi2.sf(chisq, free_deg))
            
            #We update the order to the value used in the loop, 
            #and we stop the loop is the value of the survival function reaches
            #our indicated confidence level.
            n=i
            if chi2.sf(chisq, free_deg) > conf_level:
                break
        print('Found an order of ', n,' with our chi-squared minimization routine.')

        #We instore a failsafe system in case the survival function never reaches 
        #the confidence level provided OR if the order is over 10, which generally 
        #leads to overfitting. 
        if n==int(len(counts)/2)-1 or n>10:
            #We iterate over the change in survival function. Once we reach a negative
            #value for the change in survival function AND the order is below 10, we have 
            #reached the overfitting regime and we use the corresponding value for our fit.
            for i in range(len(np.diff(SFs))):
                if np.diff(SFs)[i]<0 and i<10:
                    n=i
                    break
                    
            #We add a condition in the case where the negative change in survival function occurs at orders
            #above 10. In that case the data is still being overfit and we default to an order of 9.
                elif np.diff(SFs)[i]<0 and i>=10:
                    n=9
                    break
                else:
                    #We add a condition in the case where the change in survival function never reaches 
                    #a negative values. In that case, we use the minimum value of the survival function,
                    # skipping the first order as it will always yield a better fit.
                    n = np.diff(SFs).tolist().index(min(np.diff(SFs)[1:])) 
            print('Warning: This is an overfitting regime, '
                  'so our failsafe system implements:', n)
            
        #If our failsafe system and normal optimization system yield an order of 0 or 1, we 
        #default the value to 2. Most pulse profiles are poorly described by a single 
        #sinusoid but a 2nd order sinusoidal model is much better.
        if n<=1:
            n=2
            print('Warning: The order found was too low so we default to:', n)
        print('\n')
        return n
    
    #Create a list that will contain the harmonic phases for all the pulse profiles.
    phases = []
    
    #Create figure.
    fig, axs = plt.subplots(mt.ceil(len(segments)/4), 4, figsize=(15, 15))
    fig.subplots_adjust(hspace=.05, wspace=.1)
    axs = axs.ravel()
    
    for ax in axs:
        ax.set_axis_off()
    
    #Making a list that will contain the orders of sinusoidal model used for each segment (debugging 
    #purposes and sanity check).
    orders=[]
    
    #Making lists to store the folded counts and phases of each pulse profile - used by bootstrapping method.
    folded_counts = []
    folded_times = []
    
    #Making the pulse profile for each segment 
    #and overplotting a first-order sinusoidal model
    #using the phase of first harmonics of each segment.
    for i in range(len(segments)):
        
        #Computing the phase folded time.
        test_phase_fold = phase_fold(ref_time, segments[i], reg_coeffs)
        
        #Making the lightcurve using the phase folded time.
        test_lc = Lightcurve.make_lightcurve(test_phase_fold, dt)
        
        #Adding data to the corresponding lists.
        folded_counts.append(test_lc.counts)
        folded_times.append(test_lc.time)
            
        #Finding the optimal number of harmonics to use.
        order=find_optimal_order(test_lc.counts, test_lc.time, confi_level)
        
        #Adding this value to the corresponding list.
        orders.append(order)
   
        #Creating nth-order sinusoidal model from the order and lightcurve computed above.
        model, model_phase = model_pulse(order, test_lc.time, test_lc.counts)

        #Populating the list containing the phases of each pulse profile for future use. 
        phases.append(model_phase)

        #Plot the pulse profile and the descriptive model.
        axs[i].set_axis_on()
        #Plotting the errorbars depending on the value of the input boolean 'error'.
        if error:
            #Since the counts are drawn from a Poisson distribution, their error is simply the square root.
            axs[i].errorbar(test_lc.time, test_lc.counts, yerr = np.sqrt(test_lc.counts), fmt='.')
        else:
            axs[i].plot(test_lc.time, test_lc.counts, '.')
        axs[i].plot(test_lc.time, model)
        
        #Doing some additional formatting of the pulse profile matrix.
        if i%4 != 0:
            axs[i].set_yticklabels([])
        if i<len(segments)-4:
            axs[i].set_xticklabels([])
        axs[i].tick_params(axis='x', labelsize=13)
        axs[i].tick_params(axis='y', labelsize=13)
        fig.suptitle(title_str, y=0.92)
    fig.supxlabel('Pulse Phase', y=0.07, fontsize=16)
    fig.supylabel('Photon counts', x=0.07, fontsize=16)
    
    #Saving the pulse profile matrix to a pdf.
    if save_img:
            plt.savefig(save_img_path)
    
    return orders, phases, folded_counts, folded_times

def bootstrap_generate(counts, k):
    '''
    Function to make k realizations of a given pulse profile using a Poisson distribution
    to describe the number of counts.
    Parameters
    ----------
    :param counts: array_like, containing the counts for a given pulse profile.
    :param k: int, the number of realizations to make of the pulse profile.
    Returns
    ----------
    :param fake_profiles: nested lists, the simulated k realizations of the given pulse profile. 
    Can be used for plotting for sanity checks.
    '''
    #Setting up lists for the fake pulse profile values.
    fake_profiles = []
    
    #Making k realizations of our pulse profile using a Poisson distibution centered
    #on the number of photon counts of each point in the pulse profile.
    for i in range(len(counts)):
        fake_profiles.append(st.poisson.rvs(mu=counts[i], size=k))
    
    #Transposing the matrix of fake profiles so that we have one row for each realization
    #rather than one column - easier to read out.
    fake_profiles = np.array(fake_profiles).T
    
    return fake_profiles


def bootstrap_total(counts, k, func, *args):
    '''
    Function to calculate the standard deviation of a specific quantity obtained from
    the function func using Poisson realizations of a pulse profile.
    Parameters
    ----------
    :param counts: arraylike, containing the lightcurve counts for a given pulse profile.
    :param k: int, the number of realizations to make of the pulse profile.
    :param func: function, to calculate a certain quantity/statistic of the pulse profile.
    :param *args: additional input parameters required for func.
    Returns
    ----------
    :param std: float, standard deviation on the quantity/statistic for the given given pulse profile.
    We will use this value as the error on the quantity considered.
    '''
    
    #Get the k realizations of the pulse profile from the bootstrap_generate function.
    profiles = bootstrap_generate(counts, k)
    
    #Prepare an array that will contain the quantity we want for every realization.
    values = []
    
    #Populate the values list with the quantity using the input func function and its arguments.
    for i in range(len(profiles)):
        values.append(func(profiles[i], *args))
    
    #Calculating the standard deviation of the quantity of interest using the k realizations of it.
    std = np.sqrt(np.sum((np.array(values)-np.mean(values))**2)/(k-1.0))

    return std

def RMS_calculator(counts, k):
    '''
    Function to calculate the Root-Mean Squared (RMS) for a given pulse profile.
    Parameters
    ----------
    :param counts: array_like, containing the lightcurve counts for a given pulse profile.
    :param k: int, number of amplitudes/phases used for the sinusoidal model of the pulse
    profile we consider. Corresponds to the order of the model.
    Returns
    ----------
    :param RMS: float, RMS for the given pulse profile.
    '''
    #Fourier transform the counts.
    counts_fft = np.fft.rfft(counts)
    
    #Get the amplitudes of each harmonic.
    amps = np.sqrt(counts_fft.imag**2 + counts_fft.real**2)
    
    #Getting the RMS value from Parseval's theorem.
    RMS = np.sqrt(np.sum(amps[1:k]**2))/amps[0]
    
    return RMS

def plot_time_vs_phase(time_array, counts_array, model_phases, save_img, save_img_path, num_realizations=1000):
    '''
    Function to plot the first and second harmonic phase for each pulse profile
    in the pulse profile matrix, as a function of time.
    Parameters
    ----------
    :param time_array: array, event arrival times used to make each pulse profile in 
    the matrix.
    :param model_phases: array, all harmonic phases of each pulse profile in the matrix.
    :param save_img: bool, True if we want to save the plot to a pdf.
    :param save_img_path: string, corresponding to the path of the pdf we want to save 
    the plot. Only used if save_img=True.
    Returns
    ----------
    :param first_bestfit_params: array, Best-fit parameters from our 2nd order polynomial fit of the 
    first harmonic phases as a function of time.
    :param second_bestfit_params: array, Best-fit parameters from our 2nd order polynomial fit of the 
    second harmonic phases as a function of time.
    Plot of Phase against Time.
    '''
    
    #Creating a 2nd order polynomial to fit to the Time-phase plot.
    def quadratic(x, a, b, c):
        return a*x**2 + b*x + c
    
    #Building a model.
    mod = lf.Model(quadratic)
    
    #Extracting the first and second harmonics.
    first_harmonics = [model_phases[i][1]/(2*np.pi) for i in range(len(model_phases))]
    second_harmonics = [model_phases[i][2]/(2*np.pi) for i in range(len(model_phases))]
    
    #Extracting the average time in each bin used to make the pulse profiles.
    times = [np.mean(time_array[i] - time_array[0][0]) for i in range(len(model_phases))]
        
    #Initializing the arrays.
    first_harmonics_error = np.ones(len(first_harmonics))
    second_harmonics_error = np.ones(len(second_harmonics))
    #Getting the error on the first and second harmonic phases, with our bootstrapping method.
    for i in range(len(counts_array)):
        first_harmonics_error[i] = bootstrap_total(counts_array[i], num_realizations, harmonic_calc, 1)
        second_harmonics_error[i] = bootstrap_total(counts_array[i], num_realizations, harmonic_calc, 2)
    
    #Fitting our 2nd-order polynomial to the pulsed fraction plot.
    results_1 = mod.fit(np.array(first_harmonics), x=times, a=1, b=1, c=1, weights = 1/first_harmonics_error)
    results_2 = mod.fit(np.array(second_harmonics), x=times,  a=1, b=1, c=1, weights = 1/second_harmonics_error)
    
    #Creating a figure.
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4), sharey=True)
    fig.tight_layout()
    
    #Plotting the phase of first and second harmonics as a function of time.
    ax1.errorbar(times, first_harmonics, yerr = first_harmonics_error, fmt = '.', label='Observations')
    ax1.plot(times, results_1.best_fit, 'r', label='$2^{\mathrm{nd}}$ order polynomial fit')
    ax2.errorbar(times, second_harmonics, yerr = second_harmonics_error, fmt = '.', label='Observations')
    ax2.plot(times, results_2.best_fit, 'r', label='$2^{\mathrm{nd}}$ order polynomial fit')
    ax1.set_xlabel('Time (s)')
    ax2.set_xlabel('Time (s)')
    ax1.set_ylabel('Phase')
    ax1.set_title('First Harmonic Phase')
    ax2.set_title('Second Harmonic Phase')
    ax1.legend()
    ax2.legend()
    #Saving the plot.
    if save_img:
            plt.savefig(save_img_path, bbox_inches='tight')
    plt.show()
     
    #Extracting the best-fit parameters.
    first_bestfit_params = [results_1.params['a'].value, results_1.params['b'].value, results_1.params['c'].value]  
    second_bestfit_params = [results_2.params['a'].value, results_2.params['b'].value, results_2.params['c'].value]  

    return first_bestfit_params, second_bestfit_params

def harmonic_calc(counts, k):
    '''
    Function to calculate the phase of the k-th harmonic for
    a given pulse profile. Primarly used in our bootstrapping method.
    Parameters
    ----------
    :param counts: array, containing the counts for a given pulse profile.
    :param k: int, the order of harmonic we want to retrieve.
    Returns
    ----------
    harmonic_ph: float, phase of the k-th harmonic.
    '''
    #Getting the Fourier transform of the counts.
    counts_fft = np.fft.rfft(counts)
    
    #Getting the phases of all the harmonics. 
    phases = np.arctan2(counts_fft.imag, counts_fft.real)
    
    #Getting the phase of the k-th harmonic.
    harmonic_ph = phases[k]/(2*np.pi)
    
    return harmonic_ph


def BIC_calc(chi2, n_spectra):
    '''
    Function to calculate the Bayesian Imformation Criterion (BIC)
    given the chi-squared for our model. Important: this function 
    is used with Renkulab and as a result works with model 
    instances from xspec and pyxmmas!
    Parameters
    ----------
    :param chi2: float, chi-squared value to calcule the BIC for.
    :param n_spectra: int, the number of spectra fit.
    Returns
    ----------
    BIC: float, BIC value.
    '''
 
    #Initialize the number of data points used.
    n=0
    #Initialize the number of free parameters used.
    k=0
    
    #Summing the number of data points of each 
    #individual spectra.
    for i in range(1, n_spectra+1):
        n += len(xspec.AllData(i).noticed)

    #Get the models for each spectra.
    models = [xspec.AllModels(j) for j in range(1, n_spectra + 1)]

    #Make a function to get the number of free parameters.
    def get_frozen(nn):
        '''
        Function to calculate the number of free parameters. 
        Parameters
        ----------
        :param nn: object, model for a spectra.
        Returns
        ----------
        k: int, number of free parameters.
        '''
        #Initializing k.
        k=0
        #Iterating over the components of the model.
        for i in nn.componentNames:
            #Iterating over the parameters of each component.
            for j in getattr(nn,i).parameterNames:
                #If the parameters are frozen, the output is True
                #So we add to the sum when the output is False.
                if not getattr(getattr(nn,i), j).frozen:
                    k+=1
        return k

    #Getting the total number of free parameters for all the spectra used.
    for i in range(n_spectra):
        k+=get_frozen(models[i])
        
    #Return the BIC.
    BIC = chi2 + k*np.log(n)
    
    return BIC