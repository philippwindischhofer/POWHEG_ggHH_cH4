# pwhg_run.py

import yaml
import os
from argparse import ArgumentParser

def parse_template_file(src,dest,replacements):
    # read template file
    with open(src, 'r') as src_file:
        string = src_file.read()
    # apply replacements
    string = string % replacements
    # write parsed file
    with open(dest, 'w') as dest_file:
        dest_file.write(string)

if __name__ == '__main__':

    parser = ArgumentParser(description='starts parallel POWHEG-BOX jobs')

    parser.add_argument('-c', '--configuration', required=True, type=str, dest="CONFIGURATION",
                        metavar='Configuration YAML File', help="YAML file containing nodes/tasks per node/launch script etc.")
    parser.add_argument('-s', '--startseed', required=False, type=int, dest="STARTSEED",
                        default=1, metavar='Start seed', help="Start seed for this batch of runs [1].")
    parser.add_argument('-p', '--parallelstage', required=True, type=int, dest="PARALLELSTAGE",
                        metavar='parallelstage (for powheg.input)', help="The parallelstage to set in the powheg.input.")
    parser.add_argument('-x', '--xgriditeration', required=True, type=int, dest="XGRIDITERATION",
                        metavar='xgriditeration (for powheg.input)', help="The xgriditeration to be set in the powheg.input.")
    parser.add_argument('-d', '--dependency', required=False, type=str, dest="DEPENDENCY",
                        default='', metavar='dependency (slurm)', help="The dependency of the slurm job e.g... afterok:JOBID.")
    # Parse arguments
    args = parser.parse_args()
    configuration_file = args.CONFIGURATION
    startseed = args.STARTSEED
    parallelstage = args.PARALLELSTAGE
    xgriditeration = args.XGRIDITERATION
    dependency = args.DEPENDENCY

    # Parse configuration YAML
    with open(configuration_file, 'r') as f:
        try:
            configuration = yaml.load(f)
        except yaml.YAMLError as exc:
            print(exc)

    # Add passed parameters to the configuration dictionary
    configuration["startseed"] = startseed
    configuration["parallelstage"] = parallelstage
    configuration["xgriditeration"] = xgriditeration
    configuration["dependency"] = dependency

    # Parse file paths/names
    configuration_filename_base = os.path.splitext(os.path.basename(str(configuration_file)))[0]
    script_filename_base = os.path.splitext(str(configuration["script_template"]))[0]

    # Create code to be used to identify log files
    log_file_code = configuration_filename_base+'-s'+str(configuration["startseed"])+'-p'+\
                    str(configuration["parallelstage"])+'-x'+str(configuration["xgriditeration"])

    configuration["jobname"] = log_file_code
    configuration["launch_log_filename"] = log_file_code+".log"
    configuration["powheg_log_filename"] = "run-"+log_file_code
    configuration["script_filename"] = script_filename_base+'-'+log_file_code+".shx"
    configuration["powheg_backup_filename"] = 'powheg-'+log_file_code+'.input.bak'
    configuration["submit_command"] = configuration["batch_command"]+' '+configuration["script_filename"]

    # Create batch script
    parse_template_file(configuration["script_template"], configuration["script_filename"], configuration)

    # Create powheg.input (and backup)
    parse_template_file(configuration["powheg_input_template"],configuration["powheg_backup_filename"],configuration)
    parse_template_file(configuration["powheg_input_template"],"powheg.input",configuration)

    print('-- configuration --')
    print(configuration)
    with open(configuration["launch_log_filename"], 'w') as f:
        f.write(str(configuration))

    print('-- launching --')
    os.system(configuration["submit_command"])
