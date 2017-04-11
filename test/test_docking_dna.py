"""
Testing module for docking_dna
"""
import os
import shutil
import filecmp
from nose import with_setup
from ..docking_dna import mark_as_complete


test_scratch_folder = 'scratch'
test_project_folder = 'test_project'
test_project_name = 'test_project'


class RegressionTest(object):
    def __init__(self):
        self.path = ''
        self.test_path = ''
        self.golden_data_path = ''

    def ini_test_path(self):
        try:
            if self.test_path:
                shutil.rmtree(self.test_path)
        except OSError:
            pass
        os.makedirs(self.test_path)

    def clean_test_path(self):
        try:
            shutil.rmtree(self.test_path)
        except OSError:
            pass


class TestDockingDNA(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = os.path.join(self.path, test_scratch_folder)

        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/'

    def teardown(self):
        self.clean_test_path()

    def test_mark_as_complete(self):
        os.chdir(self.test_path)
        os.makedirs(test_project_folder)

        json_file_name = mark_as_complete(test_project_folder, test_project_name)

        assert filecmp.cmp(os.path.join(self.golden_data_path, 'results.json'),
                            os.path.join(self.test_path, test_project_folder, '.results.json'))
