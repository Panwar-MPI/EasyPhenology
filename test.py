from turtle import onclick
import pandas as pd  
import matplotlib.pyplot as plt                                                                                                              
                                                                 
from EasyPhenology import integral_smoothing                                     
from EasyPhenology import EasyPhenology                                    
                                                               
if __name__ == '__main__':                                                       
    df = pd.read_pickle('Data/df_input.pkl')  # df_input is taken from site DE-THa Germany. df_input has columns in order time, year, doy, Var
    df_smooth=  integral_smoothing(df)        #use function integral_smoothing(df) to get df_smooth with column Var smoothed values 
                                              # Function EasyPhenology(df, Threshold_value) also give df_smooth, see below
    df_pheno, df_smooth=  EasyPhenology(df, 0.5)    #Use EasyPhenology(df, Threshold_value) to get phenophases and smoothed Var value     
                                                    #Here Threhold_value=0.5                         

    #Print the smoothed Var value                           
    print(integral_smoothing(df))  
    #Print the phenophaes and smoothed Var value                           
    print(EasyPhenology(df, 0.5))


    #Plot smoothed values and phenophases                          
    fig = plt.figure(figsize=(14,10))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    #Plot smoothed values and raw values of Var
    ax1.plot(df.time,df['Var'], color='grey', alpha=0.6)
    ax1.plot(df_smooth.time,df_smooth['Var'],color='black')
    ax1.set_title('Integral Smoothing')
    ax1.set_ylabel('GPP (mgC $m^{-2}$)')
    ax1.text(df_smooth.time.mean(), 15, 'Raw values', style='italic', color='grey')
    ax1.text(df_smooth.time.mean(), 12, 'Smoothed values', style='italic', color='black')

    #Plot Phenophases
    ax2.plot(df_pheno.Year, df_pheno.SOS, 'b')
    ax2.plot(df_pheno.Year, df_pheno.POS, 'k')
    ax2.plot(df_pheno.Year, df_pheno.EOS, 'r')
    ax2.plot(df_pheno.Year, df_pheno.SOS_der, 'b--')
    ax2.plot(df_pheno.Year, df_pheno.EOS_der, 'r--')
    ax2.set_title('PTDs')
    ax2.set_xlabel('Year')  
    ax2.set_ylabel('Day of year')
    ax2.text(2000, 280, 'EOS', style='italic', color='red')
    ax2.text(2000, 50, 'SOS', style='italic', color='blue')
    ax2.text(2000, 200, 'POS', style='italic', color='black')
    ax2.text(2002, 350, 'Dashed lines= Derivative method', style='italic', color='black')
    ax2.text(2002, 400, 'Solid lines= Threshold method', style='italic', color='black')
    ax2.set_ylim(-100, 450)
    fig.tight_layout()
