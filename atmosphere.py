#========================
# Author: Kaique Fernandes
# Library: Atmosphere and Anenometry for FlightDynamics
#========================
import numpy as np
import math as m


class atmosphere:
    '''
    atmosphere library defined after International Standard Atmosphere
    '''
    def __init__(self, altitude: float):
        # constant definition
        self.R = 287
        self.k = -0.0065 # K/m 
        self.g = 9.80665 # m/s2
        self.temperature_0_K = 288.15 # 
        self.gama = 1.4 
        self.pho_0 = 1.225 # kg/m at sea level
        self.pressure_0 = 101325 # Pa
        self.altitude = altitude
        #self.true_airspeed = 100
        self.alt_tropo = 11000.0

        self.temperature_k = self.temperature()
        self.pressure_pa = self.pressure()
        self.air_density_kgpm3 = self.air_density()
        self.sound_speed_mps = self.sound_speed()

    def mach_number(self, true_airspeed: float) -> float:
        return (true_airspeed/self.sound_speed_mps)

    def true_airspeed_2_calibrated(self, true_airspeed: float) -> float:
        self.dynamic_pressure = self.pressure_pa * (((self.gama - 1)/2 * (true_airspeed/self.sound_speed_mps)**2 + 1)**(self.gama/(self.gama - 1)) - 1)

        return  (m.sqrt(2*340.262**2/(self.gama - 1) 
                *((self.dynamic_pressure/self.pressure_0 + 1)**((self.gama - 1)/self.gama)
                 - 1)))

    def temperature(self) -> float:
        if self.altitude <= self.alt_tropo:
            return self.temperature_0_K + self.k * self.altitude
        else:
            return self.temperature_0_K + self.k * self.alt_tropo

    def pressure(self) -> float:
        if self.altitude <= self.alt_tropo:
            return self.pressure_0 * (self.temperature_k/self.temperature_0_K) ** (-self.g / (self.R * self.k))
        else:
            pressure_11 = atmosphere(self.alt_tropo).pressure_pa
            return pressure_11 * m.exp(-self.g / (self.R * self.temperature_k) * (self.altitude - self.alt_tropo))

    def sound_speed(self) -> float:
        return (self.temperature_k * self.R * self.gama)**(1/2)
    
    def air_density(self) -> float:
        return self.pressure_pa / (self.R * self.temperature_k)


#if __name__ == "__main__":
#    print(atmosphere(0))