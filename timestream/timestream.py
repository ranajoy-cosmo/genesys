import numpy as np
from genesys import Genesys_Class

"""
CLASS FOR HANDLING TIMESTREAM OBJECTS
"""

class TStream(Genesys_Class):
    def __init__(self):
        pass

    def gen_t_steps(self, segment_length, segment, sampling_rate):
        """
        Convention: segment 1 is the first and starts at t=0
        Caveat: The first time sample will always start at the beginning of a segment. 
        This ensures equal sample numbers per segment
        """
        delta_t = 1.0 / sampling_rate
        n_samples = int(sampling_rate * segment_length)
        t_start = (segment - 1) * segment_length
        delta_t = 1.0 / sampling_rate
        self.t_steps = t_start + delta_t * np.arange(n_samples)
#
    #  def generate_scan_tstream(self):
        #  theta, phi, psi = self.detector.pointing.get_pointing_and_pol_angle(self.segment, self.sim_config['segment_length'], self.detector.params['sampling_rate'], self.n_samples, self.sim_config['coordinate_system'])
        #  hit_pix = self.detector.pointing.get_hit_pix(theta, phi, self.detector.sky_map.nside)
        #  sim_pol_type = self.sim_config['sim_pol_type']
        #  if sim_pol_type == 'noise':
            #  self.ts['signal'] = np.zeros(self.n_samples)
        #  if sim_pol_type == 'I':
            #  self.ts['signal'] = self.detector.sky_map.sky_map[hit_pix]
        #  if sim_pol_type == 'QU':
            #  self.ts['signal'] = self.detector.sky_map.sky_map[0][hit_pix]*np.cos(psi) + self.detector.sky_map.sky_map[1][hit_pix]*np.sin(psi)
        #  if sim_pol_type == '_QU':
            #  self.ts['signal'] = self.detector.sky_map.sky_map[1][hit_pix]*np.cos(psi) + self.detector.sky_map.sky_map[2][hit_pix]*np.sin(psi)
        #  if sim_pol_type == 'IQU':
            #  self.ts['signal'] = self.detector.sky_map.sky_map[0][hit_pix] + self.detector.sky_map.sky_map[1][hit_pix]*np.cos(psi) + self.detector.sky_map.sky_map[2][hit_pix]*np.sin(psi)
#
        #  self.ts['signal'] += self.noise.simulate_ts_noise(self.n_samples)
#
        #  self.ts['theta'] = theta
        #  self.ts['phi'] = phi
        #  self.ts['psi'] = psi
#
