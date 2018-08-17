#!/usr/bin/env python

import argparse
import os
import multiprocessing
import subprocess
import shutil
import glob
import json
from utils import logger


""" Script configuration """
tmp_folder_name = '.tmp'
uploads_folder_name = 'uploads'
json_results_file_name = '.results.json'
ini_file_script = 'prepare_ini_file.py'
pydock_bin = 'pydock3'
sampling_script = 'run_ftdock.sh'
scoring_script = 'parallel_scoring.py'
models_dest_folder = 'models'
models_prefix = 'mug_'
results_csv_file = 'result.csv'
mock_folder_dna = "/home/user/bin/mug/mock/3mfk"
mock_folder_protein = "/home/user/bin/mug/mock/3mfk_monomers"
top_models = 10
""" End of configuration """


class CommandLineParser(object):
    """Parses command line"""
    @staticmethod
    def valid_file(file_name):
        if not os.path.exists(file_name):
            raise argparse.ArgumentTypeError("The file does not exist")
        return file_name
    
    @staticmethod
    def valid_integer_number(ivalue):
        try:
            ivalue = int(ivalue)
        except:
            raise argparse.ArgumentTypeError("%s is an invalid value" % ivalue)
        if ivalue <= 0:
            raise argparse.ArgumentTypeError("%s is an invalid value" % ivalue)
        return ivalue


class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, new_path):
        self.new_path = os.path.expanduser(new_path)

    def __enter__(self):
        self.saved_path = os.getcwd()
        os.chdir(self.new_path)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.saved_path)


def read_config(config_json_file):
    data = None
    receptor_id = None
    ligand_id = None
    project_path = None
    project_name = None
    models = None
    scoring = None
    try:
        with open(config_json_file) as data_file:
            data = json.load(data_file)
            if data['input_files'][0]['name'] == 'ligand':
                ligand_id = data['input_files'][0]['value']
                receptor_id = data['input_files'][1]['value']
            else:
                ligand_id = data['input_files'][1]['value']
                receptor_id = data['input_files'][0]['value']
            for argument in data['arguments']:
                if argument['name'] == 'execution':
                    project_path = argument['value']
                    project_name = os.path.basename(project_path)
                if argument['name'] == 'models':
                    models = argument['value']
                if argument['name'] == 'scoring':
                    scoring = argument['value']
    except Exception, e:
        logger.error('Error reading config JSON: %s' % str(e))
    return receptor_id, ligand_id, project_path, project_name, int(models), scoring


def read_metadata(metadata_json_file):
    metadata = {}
    try:
        with open(metadata_json_file) as data_file:
            data = json.load(data_file)
            for argument in data:
                metadata[argument['_id']] = argument
    except Exception, e:
        logger.error('Error reading metadata JSON: %s' % str(e))
    return metadata


def get_top_from_ene(ene_file, top=10):
    """Parses the top models conformations from a .ene file"""
    top_list = []
    line_count = 0
    with open(ene_file) as input_file:
        for line in input_file:
            line_count += 1
            if line_count > 2:
                fields = line.split()
                conf = fields[0]
                top_list.append(conf)
                if line_count == top+2:
                    return top_list
    return top_list


def ene_to_csv(ene_file, csv_file, top=100, has_header=True):
    """Energy file to CSV file format"""
    with open(ene_file) as input_file:
        with open(csv_file, 'w') as output_file:
            lines_to_write=top
            if has_header:
                lines_to_write += 1
            line_count = 0
            for line in input_file:
                if line and line[0] not in ['#', '-'] and line_count <= lines_to_write:
                    output_file.write((','.join([field.strip() for field in line.split()])) + os.linesep)
                line_count += 1


def prepare_workspace(project_path, log_file):
    """Prepares the workspace"""
    logger.progress("Preparing workspace", status="RUNNING")
    # Create project path if required
    if not os.path.exists(project_path):
        os.makedirs(project_path)

    # Create temporal working path if required
    tmp_path = os.path.join(project_path, tmp_folder_name)
    if not os.path.exists(tmp_path):
        os.makedirs(tmp_path)

    # Calculate the uploads path
    source_data_path = os.path.join(project_path, uploads_folder_name)

    # Calculate the results path
    results_path = project_path

    logger.progress("Preparing workspace", status="DONE")
    return source_data_path, tmp_path, results_path


def setup_molecules(working_path, receptor_pdb, ligand_pdb, project_name):
    with cd(working_path):
        logger.progress("Setup", status="RUNNING")
        command = "%s %s %s %s" % (ini_file_script, project_name, receptor_pdb, ligand_pdb)
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        #print 'Waiting for pid %s' % str(process.pid)
        process.wait()
        command = "%s %s setup > setup.log" % (pydock_bin, project_name)
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        #print 'Waiting for pid %s' % str(process.pid)
        process.wait()
        logger.progress("Setup", status="DONE")
        return os.path.join(working_path, "%s_rec.pdb" % project_name), os.path.join(working_path, "%s_lig.pdb" % project_name)


def sampling(working_path, receptor_pdb, ligand_pdb, project_name, num_cores, mock=False, mock_folder=""):
    with cd(working_path):
        logger.progress("Sampling", status="RUNNING")
        if not mock:
            command = "%s %s %s %s %s" % (sampling_script, project_name, receptor_pdb, ligand_pdb, str(num_cores))
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
            #print 'Waiting for pid %s' % str(process.pid)
            process.wait()
            command = "%s %s rotftdock" % (pydock_bin, project_name)
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
            #print 'Waiting for pid %s' % str(process.pid)
            process.wait()
        else:
            shutil.copy2(os.path.join(mock_folder, '3mfk.ftdock'), os.path.join(working_path, "%s.ftdock" % project_name))
            shutil.copy2(os.path.join(mock_folder, '3mfk.rot'), os.path.join(working_path, "%s.rot" % project_name))
        logger.progress("Sampling", status="DONE")
        return os.path.join(working_path, "%s.ftdock" % project_name)


def scoring(working_path, project_name, num_cores, scoring_module="dockser", mock=False, mock_folder=""):
    with cd(working_path):
        logger.progress("Scoring", status="RUNNING")
        if not mock:
            command = "%s %s %s %s" % (scoring_script, project_name, str(num_cores), scoring_module)
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
            #print 'Waiting for pid %s' % str(process.pid)
            process.wait()
        else:
            shutil.copy2(os.path.join(mock_folder, '3mfk.ene'), os.path.join(working_path, "%s.ene" % project_name))
        logger.progress("Scoring", status="DONE")
        return os.path.join(working_path, "%s.ene" % project_name)


def create_top_structures(working_path, models_refix, project_name, top, file_name):
    with cd(working_path):
        with open(file_name, 'w') as output:
            num_model = 1
            for conf in top:
                try:
                    pdb_file_name = "%s%s_%s.pdb" % (models_prefix, project_name, conf)
                    with open(pdb_file_name) as input_pdb:
                        output.write('MODEL %d\n' % num_model)
                        for line in input_pdb:
                            output.write(line)
                        output.write('ENDMDL\n')
                    num_model += 1
                except IOError:
                    pass


def generate_models(working_path, project_name, num_models):
    with cd(working_path):
        logger.progress("Generating models", status="RUNNING")
        command = "%s %s makePDB 1 %s %s.ene %s" % (pydock_bin, project_name, str(num_models), project_name, models_prefix)
	#print 'Generating structures: %s' % str(command)
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        #print 'Waiting for pid %s' % str(process.pid)
        process.wait()
        # Keep the top
        top = get_top_from_ene("%s.ene" % project_name, top=top_models)
	#print 'Top structures are: %s' % str(top)
        create_top_structures(working_path, models_prefix, project_name, top, 'top_structures.pdb')
        # Create top 10
        for i in range(top_models):
            shutil.copy2("%s%s_%s.pdb" % (models_prefix, project_name, top[0]) , 'top_%d.pdb' % (i+1))
        models_path = os.path.join(working_path, models_dest_folder)
        if not os.path.exists(models_path):
            os.makedirs(models_path)
        else:
            shutil.rmtree(models_path)
        models = glob.glob("%s*.pdb" % models_prefix)
        for model in models:
            shutil.move(model, models_path)
        logger.progress("Generating models", status="DONE")
        return True


def clean_workspace(working_path, project_name):
    """Cleans the workspace from temporal folder and scratch files"""
    with cd(working_path):
        logger.progress("Cleaning", status="RUNNING")
        # Remove scoring temporal folders
        temp_folders = glob.glob('tmp_pyDock*')
        for folder in temp_folders:
            try:
                shutil.rmtree(folder)
            except:
                pass
        # Remove scratch sampling files
        scratch_files = glob.glob('scratch*')
        for scratch_file in scratch_files:
            try:
                os.remove(scratch_file)
            except:
                pass
        # Remove specific files
        try:
            os.remove('%s.ftdock.log' % project_name)
        except:
            pass
        logger.progress("Cleaning", status="DONE")


def create_compress_results(working_path, project_name):
    with cd(working_path):
        to_move = glob.glob('*')
        if not os.path.exists(project_name):
            os.makedirs(project_name)
        for thing in to_move:
            try:
                shutil.move(thing, project_name)
            except:
                pass
        command = "tar zcf %s.tgz %s" % (project_name, project_name)
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        process.wait()
        return "%s.tgz" % project_name


def prepare_results(working_path, results_path, project_name, num_models):
    with cd(working_path):
        # Clean workspace from temporal results
        clean_workspace(working_path, project_name)
        # Move top PDB to results folder
        try:
            shutil.move('top_structures.pdb', results_path)
        except:
            pass
        for i in range(top_models):
            try:
                shutil.move('top_%d.pdb' % (i+1), results_path)
            except:
                pass
        # Create CSV file
        ene_file = "%s.ene" % project_name
        csv_file = os.path.join(results_path, results_csv_file)
        ene_to_csv(ene_file, csv_file, top=num_models)
        # Create compress file
        tgz_file = create_compress_results(working_path, project_name)
        try:
            shutil.move(tgz_file, results_path)
        except:
            pass


def mark_as_complete(results_path, project_name):
    json_file_name = os.path.join(results_path, json_results_file_name)
    with open(json_file_name, 'w') as output:
        content = """
{
"output_files": [
        {
            "name": "top_structures",
            "source_id": [
                ""
            ],
            "taxon_id": "",
            "meta_data": {
            },
            "file_path": "%s/top_structures.pdb"
        },
        {
            "name": "results",
            "source_id": [
                ""
            ],
            "taxon_id": "",
            "meta_data": {
            },
            "file_path": "%s/%s.tgz"
        },
        {
            "name": "energy_table",
            "source_id": [
                ""
            ],
            "taxon_id": "",
            "meta_data": {
            },
            "file_path": "%s/%s"
        },
        {
            "name": "top10",
            "source_id": [
                ""
            ],
            "taxon_id": "",
            "meta_data": {
            },
            "file_path": "%s/top_1.pdb"
        },
        {
            "name": "top10",
            "source_id": [
                ""
            ],
            "taxon_id": "",
            "meta_data": {
            },
            "file_path": "%s/top_2.pdb"
        },
        {
            "name": "top10",
            "source_id": [
                ""
            ],
            "taxon_id": "",
            "meta_data": {
            },
            "file_path": "%s/top_3.pdb"
        },
        {
            "name": "top10",
            "source_id": [
                ""
            ],
            "taxon_id": "",
            "meta_data": {
            },
            "file_path": "%s/top_4.pdb"
        },
        {
            "name": "top10",
            "source_id": [
                ""
            ],
            "taxon_id": "",
            "meta_data": {
            },
            "file_path": "%s/top_5.pdb"
        },
        {
            "name": "top10",
            "source_id": [
                ""
            ],
            "taxon_id": "",
            "meta_data": {
            },
            "file_path": "%s/top_6.pdb"
        },
        {
            "name": "top10",
            "source_id": [
                ""
            ],
            "taxon_id": "",
            "meta_data": {
            },
            "file_path": "%s/top_7.pdb"
        },
        {
            "name": "top10",
            "source_id": [
                ""
            ],
            "taxon_id": "",
            "meta_data": {
            },
            "file_path": "%s/top_8.pdb"
        },
        {
            "name": "top10",
            "source_id": [
                ""
            ],
            "taxon_id": "",
            "meta_data": {
            },
            "file_path": "%s/top_9.pdb"
        },
        {
            "name": "top10",
            "source_id": [
                ""
            ],
            "taxon_id": "",
            "meta_data": {
            },
            "file_path": "%s/top_10.pdb"
        }
        ]
}
""" % (results_path, results_path, project_name, results_path, results_csv_file, results_path, results_path, results_path, results_path, results_path, results_path, results_path, results_path, results_path, results_path)
        output.write(content)

    return json_file_name


def check_output(file_name):
    """Check if file exists and contains actual data"""
    try:
        if os.stat(file_name).st_size > 0:
            return True
        else:
            # Empty file
            return False
    except OSError:
        # No file
        return False


def run_pipeline(args, num_cores):
    # Prepare all required parameters to run the pipeline
    receptor_id, ligand_id, project_path, project_name, num_models, scoring_function = read_config(args.config)
    
    metadata = read_metadata(args.in_metadata)

    receptor_pdb_file = metadata[receptor_id]['file_path']
    ligand_pdb_file = metadata[ligand_id]['file_path']

    # Log file
    log_file = args.log_file

    # Prepare workspace and get the relevant paths for the pipeline
    source_data_path, working_path, results_path = prepare_workspace(project_path, log_file)

    # Setup molecules
    receptor_pdb, ligand_pdb = setup_molecules(working_path, receptor_pdb_file, ligand_pdb_file, project_name)
   
    # Activate mocks if required
    mocking = False
    mock_folder = ""
    rec_file = os.path.basename(metadata[receptor_id]['file_path'])
    lig_file = os.path.basename(metadata[ligand_id]['file_path'])
    if rec_file == '3mfk_homodimer.pdb' and lig_file == '3mfk_dna.pdb':
        mocking = True
        mock_folder = mock_folder_dna
        logger.info("Mocking Protein-DNA: %s" % mock_folder)
    if (rec_file == '3mfk_monomer2.pdb' and lig_file == '3mfk_monomer1.pdb') or (rec_file == '3mfk_monomer1.pdb' and lig_file == '3mfk_monomer2.pdb'):
        mocking = True
        mock_folder = mock_folder_protein
        logger.info("Mocking Protein-Protein: %s" % mock_folder)
 
    # Sampling step
    sampling_output_file = sampling(working_path, receptor_pdb, ligand_pdb, project_name, num_cores, mocking, mock_folder)
    if not check_output(sampling_output_file):
        logger.error('Sampling process, FTDock output file not found')
        raise SystemExit

    # Energetic scoring step
    scoring_output_file = scoring(working_path, project_name, num_cores, scoring_function, mocking, mock_folder)
    if not check_output(scoring_output_file):
        logger.error('Scoring process, energy table file not found')
        raise SystemExit

    # Generating models step 
    generate_models(working_path, project_name, num_models)

    # Prepare results and cleaning step
    prepare_results(working_path, results_path, project_name, num_models)

    # Finishing the pipeline
    json_file_name = mark_as_complete(results_path, project_name)


if __name__ == "__main__":

    # Parse command line
    parser = argparse.ArgumentParser(prog="docking_dna")
        
    # Config file
    parser.add_argument("--config", help="Configuration JSON file", 
                        type=CommandLineParser.valid_file, metavar="config", required=True)
    # Metadata
    parser.add_argument("--in_metadata", help="Project metadata", metavar="in_metadata", required=True)
    # Output metadata
    parser.add_argument("--out_metadata", help="Output metadata", metavar="output_metadata", required=True)

    # Log file
    parser.add_argument("--log_file", help="Log file", metavar="log_file", required=True)

    args = parser.parse_args()

    # Number of cores available
    num_cores = multiprocessing.cpu_count()

    # Protein-DNA docking pipeline
    run_pipeline(args, num_cores)

