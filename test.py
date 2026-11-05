#========================
# Author: Kaique Fernandes
# Library: Equation of Motion for FlightDynamics
#========================
import numpy as np
import math as m
type Vector = list[float]

class EqM:
    '''
    Equation of Motions Set
    '''
    def __init__(self, 
                 simulation_time,
                 ground_forces_and_moments_N_Nm: Vector, 
                 aero_forces_and_moments_N_Nm: Vector, 
                 thrust_force_and_moments_N_Nm: Vector):
        
        # constant definition
        self.simulation_time_s = simulation_time
        self.time_step_s = 0.5
        self.time_s = 0
        self.g_mps2 = 9.80665
        self.mass_kg = 8000
        self.weight_N = np.array([0, 0, self.mass_kg * self.g_mps2]).transpose()
        self.inertia_tensor_kgm2 = np.array([[]])

        #self.LBE = np.array[LBE]
        self.ground_forces_and_moments_N_Nm = np.array([ground_forces_and_moments_N_Nm]).transpose()
        self.aero_forces_and_moments_N_Nm = np.array([aero_forces_and_moments_N_Nm]).transpose()
        self.thrust_force_and_moments_N_Nm = np.array([thrust_force_and_moments_N_Nm]).transpose()

        self.p_radps = 0
        self.q_radps = 0
        self.r_radps = 0
        self.body_rates_radps = np.array([[self.p_radps],
                                          [self.q_radps],
                                          [self.r_radps]])
        
        self.u_radps = 10
        self.v_radps = 100
        self.w_radps = 0
        self.body_velocities_mps = np.array([[self.u_radps],
                                            [self.v_radps],
                                            [self.w_radps]])

        self.theta_deg = 0
        self.psi_deg = 0
        self.phi_deg = 0
        self.euler_angles_deg = np.array([[np.degrees(self.theta_deg)],
                                            [np.degrees(self.psi_deg)],
                                            [np.degrees(self.phi_deg)]])
        
        self.test_var = 1

        for _ in range(int(self.simulation_time_s/self.time_step_s)):
            self.test_var = self.Integrator(self.test_var, self.test_func)

    def test_func(self):

        return 1

    def Integrator(self, variable, function):

        var_in_time = function()*self.time_step_s
        variable += var_in_time 

        return variable

if __name__ == '__main__':
    print(EqM(10,[0,0,0],[0,0,0],[0,0,0]))
