import numpy as np
import pandas as pd

CO2_density = 1964 # mg/l
O2_density = 1430 # mg/l
CO2_solubility = 1449 # mg/l/1atm
O2_solubility = 40 # mg/l/1atm

CO2_ratio = CO2_density / CO2_solubility
O2_ratio = O2_density / O2_solubility

CO2_MW = 44 # g/mol
O2_MW = 32

DCO2 = 1.91 # x10-5 cm2/s 25 °C
DO2 = 2.42 # x10-5 cm2/s 25 °C
Dr=DO2/DCO2

class Utils:
    
    @staticmethod
    def M_to_ppm(M, mw):
        return M * mw * 1000 
    
    @staticmethod
    def first_order_CO2(time, D, CO2_eq, CO2_start): 
        return (CO2_start-CO2_eq) * np.exp(-D*time) + CO2_eq
    
    @staticmethod
    def first_order_H(time, D, H_eq, H_start): 
        return (H_start-H_eq) * np.exp(-D*time) + H_eq
    
    @staticmethod
    def pH_to_CO2(pH, KH, pK1=6.52, pK2=10.3):
        HCO3 = Utils.KH_to_HCO3(KH)
        HCO3_corrected = Utils.get_HCO3_from_alk(HCO3, pH, pK2)
        log10_CO2 = pK1 + np.log10(HCO3_corrected) - pH
        CO2 = np.power(10, log10_CO2)
        return CO2 * 44 * 1e3
    
    @staticmethod
    def CO2_to_pH(CO2, KH, pK1=6.52, pK2=10.3):
        HCO3 = Utils.KH_to_HCO3(KH)
        # HCO3_corrected = Utils.get_HCO3_from_alk(HCO3, pH, pK2)
        CO2 /= 44 * 1e3
        pH = pK1 + np.log10(HCO3/CO2)
        return pH
    
    @staticmethod
    def KH_to_HCO3(KH):
        ppm = 21.8 * KH
        molar = Utils.ppm_to_molar(ppm, 61.0168)
        return molar
    
    @staticmethod
    def ppm_to_molar(ppm, mw):
        mass = ppm*1e-3
        return mass / mw
    
    @staticmethod
    def get_HCO3_from_alk(alk, pH, pK2=10.3):
        # alk as [HCO3]- + [CO3--]
        pH_diff = pK2 - pH
        HCO3_CO3_ratio = np.power(10, pH_diff)
        HCO3 = (alk * HCO3_CO3_ratio)/(1+HCO3_CO3_ratio)
        return HCO3
    
class Tank:
    
    def __init__(self, volume, headspace, KH=13, room_CO2_equilibrium=0.58, room_O2_equilibrium = 8.24):
        """
        Class to describe a model aquarium 
        
        Args:
            volume: water volume in liters
            headspace: headspace volume in liters
            KH: carbonate hardness in dKH units
        """
        
        self.volume = volume
        self.headspace= headspace
        self.KH = KH
        
        self.CO2_mass_ratio = (CO2_density * self.headspace) / (CO2_solubility * self.volume)
        self.O2_mass_ratio = (O2_density * self.headspace) / (O2_solubility * self.volume)
        
        self.CO2_water_concentration = room_CO2_equilibrium
        self.get_CO2_values_from_conc()

        
        self.O2_water_concentration = room_O2_equilibrium
        self.get_O2_values_from_conc()
        
        self.data_to_df()
        
    def get_CO2_values_from_conc(self):
        """
        Calculates all CO2 equilibrium values from the water CO2 concentration
        """
        
        self.CO2_water_amount = self.CO2_water_concentration * self.volume
        self.CO2_headspace_amount =(CO2_ratio * self.CO2_water_concentration) * self.headspace
        self.CO2_total_amount = self.CO2_water_amount + self.CO2_headspace_amount
        self.CO2_total_concentration = self.CO2_total_amount / (self.volume + self.headspace)
        if self.headspace > 0:
            self.CO2_headspace_concentration = self.CO2_headspace_amount / self.headspace
        else:
            self.CO2_headspace_concentration = np.nan
            
        
     
    def get_O2_values_from_conc(self):
        """
        Calculates all O2 equilibrium values from the water O2 concentration
        """
        
        self.O2_water_amount = self.O2_water_concentration * self.volume
        self.O2_headspace_amount =(O2_ratio * self.O2_water_concentration) * self.headspace
        self.O2_total_amount = self.O2_water_amount + self.O2_headspace_amount
        self.O2_total_concentration = self.O2_total_amount / (self.volume + self.headspace)
        if self.headspace > 0:
            self.O2_headspace_concentration = self.O2_headspace_amount / self.headspace
        else:
            self.O2_headspace_concentration = np.nan
            
        self.O2_amount_ratio = self.O2_headspace_amount / self.O2_water_amount
        
    def data_to_df(self):
        """
        Collects all data into a dataframe
        """
        
        cols = ["CO2 (mg/l)", "CO2 (mg)", "O2 (mg/l)", "O2 (mg)", "pH"]
        index = ["headspace", "water", "total"]
        
        headspace = [self.CO2_headspace_concentration, self.CO2_headspace_amount,
                    self.O2_headspace_concentration, self.O2_headspace_amount, np.nan]
        
        water = [self.CO2_water_concentration, self.CO2_water_amount,
                    self.O2_water_concentration, self.O2_water_amount, 
                 Utils.CO2_to_pH(self.CO2_water_concentration, self.KH)]
        
        total = [self.CO2_total_concentration, self.CO2_total_amount,
                    self.O2_total_concentration, self.O2_total_amount, np.nan]
        
        self.data = pd.DataFrame([headspace, water,total], index=index, columns=cols)
        
        self.data = self.data.round(3)
        
        
    def __repr__(self):
        
        self.data_to_df()
        
        return repr(self.data)
                
    
    def use_O2(self, O2_used, RQ=1):
        """
        Calculates the new O2 and CO2 values after consumption of an amount of O2
        
        Args:
            O2_used: the amount of O2 used by aquarium biomass in miligrams
            RQ: respiratory quotient
        """
        
        self.O2_total_amount -= O2_used
        
        if self.O2_total_amount < 0:
            self.O2_total_amount += O2_used
            
        else:  
            CO2_produced = (Utils.M_to_ppm(RQ, CO2_MW)/Utils.M_to_ppm(1, O2_MW)) * O2_used
            self.CO2_total_amount += CO2_produced
            
            self.O2_water_amount = self.O2_total_amount / (1 + self.O2_mass_ratio)
            self.O2_water_concentration = self.O2_water_amount / self.volume
            self.get_O2_values_from_conc()

            self.CO2_water_amount = self.CO2_total_amount / (1 + self.CO2_mass_ratio)
            self.CO2_water_concentration = self.CO2_water_amount / self.volume
            self.get_CO2_values_from_conc()
