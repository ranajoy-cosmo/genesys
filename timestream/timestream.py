import numpy as np
import copy
from genesys import Genesys_Class
import math
from genesys.noise.noise import Noise

"""
Class for handling timestream objects
"""

class TStream(Genesys_Class):
    def __init__(self, sim_config, detector, segment):
        self.sim_config = sim_config
        self.detector = detector
        self.segment = segment
        self.ts = {}
        self.set_n_samples()
        self.noise = Noise(detector.params['noise'])

    def generate_scan_tstream(self):
        theta, phi, psi = self.detector.pointing.get_pointing_and_pol_angle(self.segment, self.sim_config['segment_length'], self.detector.params['sampling_rate'], self.n_samples, self.sim_config['coordinate_system'])
        hit_pix = self.detector.pointing.get_hit_pix(theta, phi, self.detector.sky_map.nside)
        sim_pol_type = self.sim_config['sim_pol_type']
        if sim_pol_type == 'noise':
            self.ts['signal'] = np.zeros(self.n_samples)
        if sim_pol_type == 'I':
            self.ts['signal'] = self.detector.sky_map.sky_map[hit_pix]
        if sim_pol_type == 'QU':
            self.ts['signal'] = self.detector.sky_map.sky_map[0][hit_pix]*np.cos(psi) + self.detector.sky_map.sky_map[1][hit_pix]*np.sin(psi)
        if sim_pol_type == '_QU':
            self.ts['signal'] = self.detector.sky_map.sky_map[1][hit_pix]*np.cos(psi) + self.detector.sky_map.sky_map[2][hit_pix]*np.sin(psi)
        if sim_pol_type == 'IQU':
            self.ts['signal'] = self.detector.sky_map.sky_map[0][hit_pix] + self.detector.sky_map.sky_map[1][hit_pix]*np.cos(psi) + self.detector.sky_map.sky_map[2][hit_pix]*np.sin(psi)

        self.ts['signal'] += self.noise.simulate_ts_noise(self.n_samples)

        self.ts['theta'] = theta
        self.ts['phi'] = phi
        self.ts['psi'] = psi

    def set_n_samples(self):
        delta_t = 1.0 / self.detector.params['sampling_rate']
        self.n_samples = int(math.ceil(self.sim_config['segment_length'] / delta_t))

