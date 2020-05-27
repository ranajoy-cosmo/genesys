import os
import shutil
import argparse
from genesys import global_paths, load_param_file

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', type=str, help='Name of config file under relative path config_files/')
    parser.add_argument('out_file', type=str, help='Name of output file in the current directory')
    in_args = parser.parse_args()

    current_dir = os.path.dirname(__file__)
    config_file = os.path.join(current_dir,'config_files', in_args.config_file)
    config = load_param_file(file_path=config_file)

    tod_dir = os.path.join(global_paths['output_dir'], config['sim_tag'], "tod")

    instrument_dir = os.path.join(global_paths['base_dir'], 'instruments', config['instrument_dir'])
    
    shutil.move(in_args.out_file, os.path.join(tod_dir, "out.txt"))
    shutil.copy(config_file, os.path.join(tod_dir, "config.yaml"))
    shutil.copytree(instrument_dir, os.path.join(tod_dir, "instrument_params"))
