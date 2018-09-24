#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# This routine loads the config and checks for any inconsistencies or errors before the simulation starts
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
def pre_sim_general_config_check(config):
    check_flag = 1

    # Checking main config variables
    main_var_list = ['sim_tag', 'scan_tag', 'overwrite', 'sim_pol_type', 'coordinate_system', 'tod_type', 'oversampling_rate', 'notes', 'focal_plane_config', 'noise_type', 't_segment', 'detector_segment_dict', 'ts_data_products', 'output_dir', 'beam_type', 'write_beam']
    main_var_exists = list(map(lambda var_name: var_name in dir(config), main_var_list))
    if False in main_var_exists:
        check_flag = 0 
        if rank == 0:
            prompt(colored("WARNING: ", color="red") + "{} is not defined.\n".format(list(compress(main_var_list, np.logical_not(main_var_exists)))))

    sim_pol_type_list = ['IQU', 'QU', 'I', '_QU', 'noise_only']
    if config.sim_pol_type not in sim_pol_type_list:
        check_flag = 0
        if rank == 0:
            prompt(colored("WARNING: ", color="red") + "sim_pol_ype is defines as \"{}\". It should be one of {}.\n".format(config.sim_pol_type, sim_pol_type_list))
    if config.sim_pol_type != "noise_only" and "nside_in" not in dir(config):
        check_flag = 0
        if rank == 0:
            prompt(colored("WARNING: ", color="red") + "nside_in is not defined\n")

    if config.tod_type not in ['signal', 'gradient']:
        ckeck_flag = 0 
        if rank == 0:
            prompt(colored("WARNING: ", color="red") + "tod_type is defined as \"{}\". It should be either \"signal\" or \"gradient\".\n".format(config.sim_pol_type))
    if config.tod_type == "gradient":
        grad_type_list = ['zero_order', 'grad_co', 'grad_cross', 'grad_coXgrad_co', 'grad_crossXgrad_cross', 'grad_coXgrad_cross']
        if not set(config.grad_type).issubset(set(grad_type_list)):
            check_flag = 0
            if rank == 0:
                prompt(colored("WARNING: ", color="red") + "grad_type is set to {}. Please check as it is not one of the defined types\n".format(config.grad_type))

    if config.coordinate_system not in ['ecliptic', 'galactic']:
        check_flag = 0
        if rank == 0:
            prompt(colored("WARNING: ", color="red") + "coordinate_system is defined as \"{}\". It should be either \"ecliptic\" or \"galactic\".\n".format(config.coordinate_system))
    #-----------------------------------------------------------

    # Checking that the scan strategy variables are in place
    if config.scan_strategy_module != None:
        scan_strategy = importlib.import_module("genesys.scan_strategy." + config.scan_strategy_module).scan_strategy
        config.__dict__.update(scan_strategy.__dict__)

    scan_strategy_var_list = ['t_year', 't_prec', 't_spin', 'sampling_rate', 'alpha', 'beta', 'scan_strategy_name', 'scan_strategy_note']
    scan_strategy_var_exists = list(map(lambda var_name: var_name in dir(config), scan_strategy_var_list))
    if False in scan_strategy_var_exists:
        check_flag = 0
        if rank == 0:
            prompt(colored("WARNING: ", color="red") + "{} is not defined for the scan strategy.\n".format(list(compress(scan_strategy_var_list, np.logical_not(scan_strategy_var_exists)))))
    #-----------------------------------------------------------

    # Checking that the I/O variables are in place
    ts_data_product_list = ['signal', 'pointing_vec', 'pol_ang', 'noise', 'mask']
    if not set(config.ts_data_products).issubset(set(ts_data_product_list)):
        check_flag = 0
        if rank == 0:
            prompt(colored("WARNING: ", color="red") + "ts_data_products are {}. Please check as it is not one of the defined types\n".format(config.ts_data_products))

    if config.output_dir == None:
        config.output_dir = global_paths.data_dir
    #-----------------------------------------------------------

    # Noise and beam settings
    beam_type_list = ['pencil', 'full_simulated', 'from_file']
    if not config.beam_type in beam_type_list:
        check_flag = 0
        if rank == 0:
            prompt(colored("WARNING: ", color="red") + "beam_type is {}. Please check as it is not one of the defined types\n".format(config.beam_type))
    
    if config.beam_type == "full_simulated":
        if "beam_cutoff" not in dir(config):
            check_flag = 0
            if rank == 0:
                prompt(colored("WARNING: ", color="red") + "beam_cutoff not defined")

    noise_type_list = ['white', '1_over_f', 'none']
    if config.noise_type not in noise_type_list:
        check_flag = 0
        if rank == 0:
            prompt(colored("WARNING: ", color="red") + "noise_type is {}. Please check as it is not one of the defined types\n".format(config.noise_type))
    #-----------------------------------------------------------

    return check_flag
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def pre_sim_detector_config_check(config):
    return
