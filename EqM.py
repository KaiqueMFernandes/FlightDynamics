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
                 simulation_time: float,
                 mass_kg: float,
                 inertia_tensor_kgm2: Vector,
                 ground_forces_and_moments_N_Nm: Vector, 
                 aero_forces_and_moments_N_Nm: Vector, 
                 thrust_force_and_moments_N_Nm: Vector):
        
        # constant definition
        self.simulation_time_s = simulation_time
        self.time_step_s = 0.1
        self.time_s = 0
        self.g_mps2 = 9.80665
        self.mass_kg = mass_kg
        self.weight_N = np.array([0, 0, self.mass_kg * self.g_mps2]).transpose()
        self.inertia_tensor_kgm2 = np.array([inertia_tensor_kgm2])

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
        
        self.u_radps = 0
        self.v_radps = 0
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
        
        self.x_m = 0
        self.y_m = 0
        self.z_m = 0
        self.position_deg = np.array([[np.degrees(self.x_m)],
                                            [np.degrees(self.y_m)],
                                            [np.degrees(self.z_m)]])
        
        for i in range(int(self.simulation_time_s/self.time_step_s)):
            self.body_rates_radps = self.Body_Rates()
            self.body_velocities_mps = self.Body_Velocities()


    def External_forces_and_moments(self) -> Vector:
        self.weight_N_body = self.Body_2_Earth(self.weight_N)

        self.external_forces_N_body = (self.weight_N_body + self.ground_forces_and_moments_N_Nm[:3] +
                                                    self.aero_forces_and_moments_N_Nm[:3] +
                                                    self.thrust_force_and_moments_N_Nm[:3])
        
        self.external_moments_Nm_body = (self.ground_forces_and_moments_N_Nm[3:6] +
                                self.aero_forces_and_moments_N_Nm[3:6] +
                                self.thrust_force_and_moments_N_Nm[3:6])
        
        def load_factor(self) -> Vector:
            x = self.external_forces_N_body[0] / self.weight_N_body
            y = self.external_forces_N_body[1] / self.weight_N_body
            z = - self.external_forces_N_body[2] / self.weight_N_body
            load_factors = [x, y, z]
            return np.array([load_factors]).transpose    
        
        return self.external_forces_N_body, self.external_moments_Nm_body, load_factor()

    def Body_Rates(self) -> Vector:
        moments = self.External_forces_and_moments()[1]
        
        product = np.cross(self.body_rates_radps, 
                            np.matmul(self.inertia_tensor_kgm2, self.body_rates_radps))
        
        body_rates_radps = np.matmul(np.invert(self.inertia_tensor_kgm2), 
                                    (moments - product))

        return body_rates_radps
    
    def Body_Velocities(self) -> Vector:
        forces = self.External_forces_and_moments()[0]     
    
        product = np.cross(self.body_rates_radps, self.body_velocities_mps)
        
        body_velocities_mps = (forces/self.mass_kg - product)

        return body_velocities_mps

    def Inertial_Data(self) -> Vector:
        self.inertial_velocity_mps = self.Body_2_Earth(self.Body_Velocities())
        
        x_m = self.Integrator(self.x_m, self.inertial_velocity_mps)[0]
        y_m = self.Integrator(self.y_m, self.inertial_velocity_mps)[1]
        z_m = self.Integrator(self.z_m, self.inertial_velocity_mps)[2]

        iterated_position = np.array([[x_m],
                                    [y_m],
                                    [z_m]])

        return iterated_position

    def Trajectory_Data(self):
        self.ground_speed_mps = (self.inertial_velocity_mps[1]** 2 + self.inertial_velocity_mps[0]**2)**0.5
        self.gamma_deg = m.degrees(m.atan(-self.inertial_velocity_mps[2] / self.ground_speed_mps))

        self.track_angle_deg = m.degrees(m.atan(-self.inertial_velocity_mps[1] / self.inertial_velocity_mps[0]))

    def Euler_Angles(self):
        p = self.p_radps
        q = self.q_radps
        r = self.r_radps
        phi = m.radians(self.phi_deg)
        psi = m.radians(self.psi_deg)
        theta = m.radians(self.theta_deg)

        phi_rad = self.Integrator(phi, Euler_Rates)[0]
        psi_rad = self.Integrator(psi, Euler_Rates)[1]
        theta_rad = self.Integrator(theta, Euler_Rates)[2]
        
        def Euler_Rates(self):
            if m.cos(theta) == 0:
                phi_dot = (p + 
                         q * m.sin(phi) * m.tan(theta) + 
                         r * m.cos(phi) * m.tan(theta))

                theta_dot = (q * m.cos(phi) - 
                            r * m.sin(phi))

                psi_dot = ((q * m.sin(phi) + 
                            r * m.cos(phi)) / 
                            0.001)

            else:
                phi_dot = (p + 
                         q * m.sin(phi) * m.tan(theta) + 
                         r * m.cos(phi) * m.tan(theta))

                theta_dot = (q * m.cos(phi) - 
                            r * m.sin(phi))

                psi_dot = ((q * m.sin(phi) + 
                            r * m.cos(phi)) / 
                            m.cos(theta))

            euler_rates_radps = np.array([[phi_dot],
                                        [theta_dot],
                                        [psi_dot]])
            return euler_rates_radps  
        
        euler_angles_rad = np.array([[phi_rad],
                                    [theta_rad],
                                    [psi_rad]])

        return euler_angles_rad

    def Body_2_Earth(self, matrix: Vector):
        phi = m.radians(self.phi_deg)
        psi = m.radians(self.psi_deg)
        theta = m.radians(self.theta_deg)

        c = m.cos
        s = m.sin

        L_E_B = np.array([[c(theta)*c(psi), s(phi)*s(theta)*c(psi) - c(phi)*s(psi), c(phi)*s(theta)*c(psi) + s(phi)*s(psi)],
                          [c(theta)*s(psi), c(phi)*c(psi) + s(phi)*s(theta)*s(psi), -s(phi)*c(psi) + c(phi)*s(theta)*s(psi)],
                          [-s(theta),       s(phi)*c(theta),                        c(phi)*c(theta)]])

        return np.matmul(L_E_B, matrix)

    def Integrator(self, variable, function):
        """
        Numerical Integrator - sums currently value to iteration value
        
        Parameters
        ----------

        variable : variable to be added to.
        function : function to be iterated.

        """
        delta = function()*self.time_step_s
        variable += delta 

        return variable

    def Integrator_mid(self, function: function, a: float, b: float):
        """
        Numerical Integrator using middle point method
        
        Parameters
        ----------

        function : function to be integrated.
        a : lower limit of integration.
        b : upper limit of integration.

        """
        n = 100
        delta_x = (b-a)/n
        sum = 0

        for i in range(n):
            sum += function(i*delta_x+a) * delta_x
        
        return sum

