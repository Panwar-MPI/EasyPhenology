#EasyPhenology

'''
This script contains functions for calculating Phnological transition dates (PTDS) that are Start of season (SOS), peak of season (POS),
end of season (EOS) and growing season length (GSL) using the daily timeseries of Gross primary
production (GPP) in eddy covariance data.

It can be applied to any other variable that is in daily time step, example GCC, LAI etc

The script contains the following primary functions:

     1. integral_smoothing(df,knots):  Smooths the daily value of GPP using our integral smoothing method.
        This function requires dataframe with columns 'time', 'year', 'doy' and 'Var', in the given order.
        These values are in daily time scale. integral_smoothing function smooths variable Var, in this case GPP.
        First step is to interpolate the missing value. Second step is to calculate the cumulative of the variable
        for each calender year. The cumulative of the variable is then smoothed. Taking the differentiation of smoothed 
        cumulative value the smoothed value of Var is then obtained. Spline function is used for smoothing. 
        spline is a piecewise regression, thereby it is sensitive to the number of knots, provided by the user. 
        Usually, knots for an annual time series can vary from 8 to 15.


     2. direct_smoothing(df,knots): Smooths the daily value of GPP using the traditional smoothing method.
        This function requires dataframe with columns 'time', 'year', 'doy' and 'Var', in the given order.
        This functions smooths variable Var, in this case GPP. The first step is to interpolate the missing value.
        Then the smoothing of daily GPP signal is done for each calender month. Spline function is used for smoothing.
         


     3. EasyPhenology(df, Threshold_value, Smoothing,knots): It calculates PTDs using threshold and first derivative method. 
       The input is a dataframe with columns in order 'time', 'year', 'doy' and 'Var'. Here Var is GPP. The Threshold_value 
       is the fixed percentage of the annual maximum GPP and can vary from 0 to 1. If smoothing='True', integral smoothing is used,
       if smoothing='False', then the direct smoothing method is used. 

       The outputs are two dataframes df_pheno_out and df_smooth.
       The first dataframe contains columns "Year", "SOS", "POS", "EOS", "GSL", "SOS_der", "EOS_der" and "GSL_der".
       Here 'der' refers to PTDs calculated using derivative method, The second dataframe output is the smoothed value
       of dataframe using integral smoothing. Less than 0 values indicate days in the previous years and more than 365 value
       indicates days in the next year. See the function for more information.

Links:
Gitlab: https://git.bgc-jena.mpg.de/apanwar/phenofeedbacks.git


Contact:
Annu Panwar: apanwar@bgc-jena.mpg.de


'''


#import required libraries
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import linregress
from scipy import stats
from scipy import signal
from scipy.signal import argrelextrema
import warnings
pd.options.mode.chained_assignment = None  # default='warn'
warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=RuntimeWarning)

from tsmoothie.smoother import *


def integral_smoothing(df, knots):

    '''

    Smooths the daily value of GPP using integral smoothing

    This function requires dataframe with columns 'time', 'year', 'doy' and 'Var', in the given order.
    These values should be in daily time scale. This functions smooths variable Var. First step is to interpolate
    the missing value. Second step is to calculatee cumulative of the variable Var for each year. The
    cumulative Var is then smoothed. Taking the differentiation of smoothed cumulative value the smoothed value of
    Var is obtained.

    Here Var is GPP.

    Parameters
    ----------
    df : a dataframe with columns : time, year, doy, Var. Var is the variable of interest

    Returns
    -------
    df_smooth : the dataframe with columns time,year, doy and Var, where Var is the smoothed values

    '''

    # Create empty list for Var_smooth
    Var_smooth=[]

    # Create empty dataframe for smoothed values of Variable Var
    df_smooth=pd.DataFrame()

    # Total number of years in the site
    Var_Years=np.unique(df['year'])

    for j, element in enumerate(Var_Years):  # Loop for each year
        df_jyear=df.loc[df['year'] ==Var_Years[j]] #the dataframe for year j

        #Loop only if year has nan values for less than 50 days
        if df_jyear['Var'].isnull().sum()<50:

            # If it is not the first year (j>0) then also get the points where the threshold line crosses the previous year
            if j>0:
                df_jyear_p=df.loc[df['year'] ==Var_Years[j-1]] #the dataframe for year j-1

            # If it is not the last year (j<len(Var_Years)-1) then also get the points where the threshold line crosses the next year
            if j<len(Var_Years)-1:
                df_jyear_n=df.loc[df['year'] ==Var_Years[j+1]] #the dataframe for year j+1

            #For smoothing add 20 days before and after the year from the same year end points
            if j>0:
                df_jyear_b= df_jyear_p.tail(20)
            else:
                df_jyear_b= df_jyear.head(20)

            if j<len(Var_Years)-1:
                df_jyear_a= df_jyear_n.head(20)
            else:
                df_jyear_a= df_jyear.tail(20)


            df_jyear=pd.concat([df_jyear_b, df_jyear, df_jyear_a ])

            #interpolate if there are some nan
            df_jyear['Var']=df_jyear['Var'].interpolate(method='linear')

            # Calculate Cumulative of Variable from raw Var values
            Var_cum = df_jyear['Var'].cumsum()

            #Smooth the cumulative signal with savgol
            #Cumulative curve is already quite smooth using different smoothing method can be tested
            #Var_cum_smooth= signal.savgol_filter(Var_cum, window_length=201, polyorder=3, mode="nearest")
            #Var_smooth=np.gradient(Var_cum_smooth,1)


            #smooth the signal with spline
            smoother = SplineSmoother(n_knots=knots, spline_type='natural_cubic_spline')
            smoother.smooth(Var_cum)
            Var_cum_smooth=smoother.smooth_data[0]


            #Take the differentiation of the smoothed cumulative values
            Var_smooth=np.gradient(Var_cum_smooth,1)
            df_jyear['Var']=Var_smooth

            #Remove extra 10 days before and after the yearly values
            df_jyear=df_jyear.iloc[20:385,:]

            #Store the smoothed values for each year
            df_smooth = pd.concat([df_smooth, df_jyear])

        else:

            # if year has nan values for more than 50 days return all nan
            df_jyear["Var"]=np.nan

            #Store nan values for years with more than 50 days of nan values
            df_smooth = pd.concat([df_smooth, df_jyear])

    return df_smooth






def direct_smoothing(df,knots):

    '''

        Smooths the daily value of GPP using the traditional smoothing method.
        This function requires dataframe with columns 'time', 'year', 'doy' and 'Var', in the given order.
        This functions smooths variable Var, in this case GPP. The first step is to interpolate the missing value.
        Then the smoothing of daily GPP signal is done for each calender month. Spline function is used for smoothing. 


     Parameters
     ----------
        df : a dataframe with columns : time, year, doy, Var. Var is the variable of interest

     Returns
        -------
        df_smooth : the dataframe with columns time,year, doy and Var, where Var is the smoothed values

    '''

    # Create empty list for Var_smooth
    Var_smooth=[]

    # Create empty dataframe for smoothed values of Variable Var
    df_smooth=pd.DataFrame()

    # Total number of years in the site
    Var_Years=np.unique(df['year'])

    for j, element in enumerate(Var_Years):  # Loop for each year
        df_jyear=df.loc[df['year'] ==Var_Years[j]] #the dataframe for year j

        #Loop only if year has nan values for less than 50 days
        if df_jyear['Var'].isnull().sum()<50:

            # If it is not the first year (j>0) then also get the points where the threshold line crosses the previous year
            if j>0:
                df_jyear_p=df.loc[df['year'] ==Var_Years[j-1]] #the dataframe for year j-1

            # If it is not the last year (j<len(Var_Years)-1) then also get the points where the threshold line crosses the next year
            if j<len(Var_Years)-1:
                df_jyear_n=df.loc[df['year'] ==Var_Years[j+1]] #the dataframe for year j+1

            #For smoothing add 20 days before and after the year from the same year end points
            if j>0:
                df_jyear_b= df_jyear_p.tail(20)
            else:
                df_jyear_b= df_jyear.head(20)

            if j<len(Var_Years)-1:
                df_jyear_a= df_jyear_n.head(20)
            else:
                df_jyear_a= df_jyear.tail(20)


            df_jyear=pd.concat([df_jyear_b, df_jyear, df_jyear_a ])

            #interpolate if there are some nan
            df_jyear['Var']=df_jyear['Var'].interpolate(method='linear')

            # Calculate Cumulative of Variable from raw Var values
            Var_cum = df_jyear['Var'].cumsum()

            #Smooth the cumulative signal
            #Cumulative curve is already quite smooth using different smoothing method can be tested
            #Var_cum_smooth= signal.savgol_filter(Var_cum, window_length=201, polyorder=3, mode="nearest")
            #Var_smooth=signal.savgol_filter(df_jyear['Var'], window_length=201, polyorder=3, mode="nearest")

            #smooth the signal with spline
            smoother = SplineSmoother(n_knots=knots, spline_type='natural_cubic_spline')
            smoother.smooth(df_jyear['Var'])
            Var_smooth=smoother.smooth_data[0]

            #Take the differentiation of the smoothed cumulative values
            #Var_smooth=np.gradient(Var_cum_smooth,1)
            df_jyear['Var']=Var_smooth

            #Remove extra 10 days before and after the yearly values
            df_jyear=df_jyear.iloc[20:385,:]

            #Store the smoothed values for each year
            df_smooth = pd.concat([df_smooth, df_jyear])

        else:

            # if year has nan values for more than 50 days return all nan
            df_jyear["Var"]=np.nan

            #Store nan values for years with more than 50 days of nan values
            df_smooth = pd.concat([df_smooth, df_jyear])

    return df_smooth








def EasyPhenology(df, Threshold_value, Smoothing,knots):

    '''

    Produces phenological transition dates (PTDs) for the given dataframe

    This function requires dataframe with columns 'time', 'year', 'doy' and 'Var'. These values
    should be in daily time scale. This function first smooths Var using integral smoothing and
    then the PTDs are produced.
    The Threshold_value is the fixed percentage of the annual maximum GPP and can vary from 0 to 1. 
    If smoothing='True', integral smoothing is used,
    if smoothing='False', then the direct smoothing method is used. 

    The outputs are two dataframes df_pheno_out and df_smooth.
    The first dataframe contains columns "Year", "SOS", "POS", "EOS", "GSL", "SOS_der", "EOS_der" and "GSL_der".
    Here 'der' refers to PTDs calculated using derivative method, The second dataframe output is the smoothed value
    of dataframe using integral smoothing. Less than 0 values indicate days in the previous years and more than 365 value
    indicates days in the next year. See the function for more information.

    Returns
    -------
    df_pheno_out: the dataframe with columns "Year", "SOS", "POS", "EOS", "GSL", "SOS_der", "EOS_der","GSL_der"
    df_smooth : the dataframe with columns time,year, doy and Var, where Var is the smoothed values
    '''


    # Define rows in the output dataframe, Pheno_out
    rows = []

    #number of total years
    Var_Years=np.unique(df['year'])

    #Smooth Var using integral_smoothing
    if  Smoothing == "True":
        df=integral_smoothing(df, knots)

    #Take the raw Var
    if  Smoothing == "False":
        df=direct_smoothing(df, knots)


    for j, element in enumerate(Var_Years): #Loop over years

        df_jyear=df.loc[df['year'] ==Var_Years[j]]  #the dataframe for year j

        # peak of season, POS
        pos=df.iloc[df_jyear.Var.argmax(), 2] # it is when the maximum of var (GPP here) occurs


        #SOS, EOS, GSL from threshold method

        #Normalize the year such that all values lie between 0 to 1, 1 being value at pos
        Var_n=(df_jyear.Var -min(df_jyear.Var ))/(max(df_jyear.Var )-min(df_jyear.Var )  )

        #Points threshold line crosses
        line_cross=np.argwhere(np.diff(np.sign(Var_n - Threshold_value))).flatten()

        # If it is not the first year (j>0) then also get the points where the threshold line crosses the previous year
        if j>0:
            df_jyear_p=df.loc[df['year'] ==Var_Years[j-1]] #the dataframe for year j-1

            #Normalize the previous year such that all values lie between 0 to 1, 1 being value at pos in year j
            Var_n_p=(df_jyear_p.Var -min(df_jyear.Var ))/(max(df_jyear.Var )-min(df_jyear.Var )  )

            #The points threshold line crosses in previous year
            line_cross_p=np.argwhere(np.diff(np.sign(Var_n_p - Threshold_value))).flatten()

        # If it is not the last year (j<len(Var_Years)-1) then also get the points where the threshold line crosses the next year
        if j<len(Var_Years)-1:
            df_jyear_n=df.loc[df['year'] ==Var_Years[j+1]] #the dataframe for year j+1

            #Normalize the next year such that all values lie between 0 to 1, 1 being value at pos in year j
            Var_n_n=(df_jyear_n.Var -min(df_jyear.Var ))/(max(df_jyear.Var )-min(df_jyear.Var )  )

            #The points threshold line crosses in next year
            line_cross_n=np.argwhere(np.diff(np.sign(Var_n_n - Threshold_value))).flatten()

        # If number of nan in the year is less than 50 then only start the search for PTDs
        if df_jyear['Var'].isnull().sum()<50:
            eos=9999 #false value
            sos=9999 #false value
            gsl=9999

             #CASE 1: If threshold line crosses 2 points, normal case, sos<pos and eos>pos
            if len(line_cross)==2 and line_cross[0]<pos and line_cross[1]>pos: # order sos-pos-eos
                sos=line_cross[0]+1
                eos=line_cross[1]+1


            #CASE 2: If threshold line crosses >=1 point before pos
            if len(line_cross)>=1 and line_cross[-1]<pos and j<len(Var_Years): #order sos-pos-eos(next year)
                sos=line_cross[-1] +1

                if len(line_cross_n)>=1: # If the threshold line passes through next year
                    eos=line_cross_n[0] +1 #eos calculated from next year, , not possible for last year (j<len(Var_Years))

            #CASE 3: If threshold line crosses >=1 point after pos,
            if len(line_cross)>=1 and line_cross[0]>pos and j>0: #order sos(previous year)-pos-eos
                eos=line_cross[0] +1
                if len(line_cross_p)>=1: # If the threshold line passes through previous year
                    sos=line_cross_p[-1] +1 #sos calculated from previous year, not possible for first year (j>0)


            #CASE 4: If threshold line crosses >= 2 points, and pos is in between
            if len(line_cross)>=2 and line_cross[0]<pos and line_cross[-1]>pos: # order sos-pos-eos
                sos= line_cross[line_cross < pos].min() +1
                eos= line_cross[line_cross > pos].min() +1





            # growing season length
            #For CASE 1 and CASE 4
            if eos<9999 and sos<9999 and sos<pos and eos>pos:
                gsl=eos-sos

            #For CASE 2
            if eos<9999 and sos<9999 and sos<pos and eos<pos:
                gsl=365- sos + eos

            #FOR CASE 3
            if eos<9999 and sos<9999 and sos>pos and eos>pos:
                gsl= eos + 365-sos



            #SOS, EOS, GSL from first derivative method

            eos_der=9999
            sos_der=9999


            #CASE 1: If threshold line crosses 2 points, normal case, sos<pos and eos>pos
            if sos<pos and eos>pos: # order sos-pos-eos
                Var_1der=np.gradient(df_jyear.Var,1)
                Var_1der= signal.savgol_filter(Var_1der, window_length=101, polyorder=3, mode="nearest")
                sos_der=Var_1der[0:pos].argmax() +1
                eos_der=Var_1der[pos:-1].argmin() +1
                

            #CASE 2: If threshold line crosses >=1 point before pos
            if len(line_cross)>=1 and line_cross[-1]<pos and j<len(Var_Years): #order sos-pos-eos(next year)
                Var_1der=np.gradient(df_jyear.Var,1)
                Var_1der= signal.savgol_filter(Var_1der, window_length=101, polyorder=3, mode="nearest")
                sos_der=Var_1der[0:pos].argmax() +1

                if len(line_cross_n)>=1: # If the threshold line passes through next year
                    frames = [df_jyear[pos:-1], df_jyear_n]
                    result = pd.concat(frames)
                    result=result.reset_index(drop=True)
                    Var_1der=np.gradient(result.Var,1)
                    Var_1der= signal.savgol_filter(Var_1der, window_length=101, polyorder=3, mode="nearest")
                    eos_der=Var_1der[pos:-1].argmin() + 1


            #CASE 3: If threshold line crosses >=1 point after pos,
            if len(line_cross)>=1 and line_cross[0]>pos and j>0: #order sos(previous year)-pos-eos
                Var_1der=np.gradient(df_jyear.Var,1)
                Var_1der= signal.savgol_filter(Var_1der, window_length=101, polyorder=3, mode="nearest")
                eos_der=Var_1der.argmin() + 1

                if len(line_cross_p)>=1: # If the threshold line passes through previous year
                    sos=line_cross_p[-1] +1 #sos calculated from previous year, not possible for first year (j>0)
                    frames = [df_jyear_p[eos_der:-1], df_jyear[0:pos]] #eos of previous year to pos of present year
                    result = pd.concat(frames)
                    result=result.reset_index(drop=True)
                    Var_1der=np.gradient(result.Var,1)
                    Var_1der= signal.savgol_filter(Var_1der, window_length=101, polyorder=3, mode="nearest")
                    sos_der=Var_1der.argmax() #position of sos_der
                    #correct for previouse or same year
                    pyl=len(df_jyear_p[eos_der:-1]) #length of previous year from eos to end of season
                    sos_der=sos_der-len(df_jyear_p[eos_der:-1]) 

            #CASE 4: If threshold line crosses >= 2 points, and pos is in between

            if len(line_cross)>=2 and line_cross[0]<pos and line_cross[-1]>pos: # order sos-pos-eos
                Var_1der=np.gradient(df_jyear.Var,1)
                Var_1der= signal.savgol_filter(Var_1der, window_length=201, polyorder=3, mode="nearest")
                sos_der=Var_1der[0:pos].argmax() +1
                eos_der=Var_1der.argmin() +1

            # growing season length
            gsl_der=eos_der -sos_der

            if gsl_der<0:
                gsl_der=gsl_der+365


            rows.append([Var_Years[j], sos, pos, eos, gsl, sos_der, eos_der, gsl_der]) #Raw sos, pos, eos [0 to 365]

        # If number of nan in the year is greater than 50 then PTDs are nan
        else:
            rows.append([Var_Years[j], np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])

        df_pheno = pd.DataFrame(rows, columns=["Year", "SOS", "POS", "EOS", "GSL", "SOS_der", "EOS_der","GSL_der"])
        df_pheno_out = pd.DataFrame(rows, columns=["Year", "SOS", "POS", "EOS", "GSL", "SOS_der", "EOS_der","GSL_der"])

        #correct for order sos,pos,eos

        #eos<pos ( eos in next year), CASE 2
        df_pheno_out['EOS'] = df_pheno_out.apply(lambda x: x['EOS'] +365 if x['EOS'] <= x['POS'] and x['EOS']< 9999 else x['EOS'], axis=1)

        #sos>pos ( sos in previous year), CASE 3
        df_pheno_out['SOS'] = df_pheno_out.apply(lambda x: x['SOS'] -365 if x['SOS'] >= x['POS'] and x['SOS']< 9999 else x['SOS'], axis=1)

        #Replace default value of 9999 to nan
        #The dataframe with smoothed Var [ output 1]
        df_pheno_out=df_pheno_out.replace({9999:np.nan})

        #The dataframe with smoothed Var [ output 2]
        df_smooth=df

    return df_pheno_out, df_smooth




#Check the pylint score of the code
#import pylint.lint
#pylint_opts = ['--disable=line-too-long', 'EasyPhenology.py']
#pylint.lint.Run(pylint_opts)
