#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   rfd.py
## @brief  main test files for RFdiffusion
## @author Sergey Lyskov


import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py')  # A bit of Python magic here, what we trying to say is this: from __init__ import *, but init is calculated from file location

_api_version_ = '1.0'

import os, tempfile
import urllib.request


_models_urls_ = '''
http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt
http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt
http://files.ipd.uw.edu/pub/RFdiffusion/60f09a193fb5e5ccdc4980417708dbab/Complex_Fold_base_ckpt.pt
http://files.ipd.uw.edu/pub/RFdiffusion/74f51cfb8b440f50d70878e05361d8f0/InpaintSeq_ckpt.pt
http://files.ipd.uw.edu/pub/RFdiffusion/76d00716416567174cdb7ca96e208296/InpaintSeq_Fold_ckpt.pt
http://files.ipd.uw.edu/pub/RFdiffusion/5532d2e1f3a4738decd58b19d633b3c3/ActiveSite_ckpt.pt
http://files.ipd.uw.edu/pub/RFdiffusion/12fc204edeae5b57713c5ad7dcb97d39/Base_epoch8_ckpt.pt
http://files.ipd.uw.edu/pub/RFdiffusion/f572d396fae9206628714fb2ce00f72e/Complex_beta_ckpt.pt
http://files.ipd.uw.edu/pub/RFdiffusion/1befcb9b28e2f778f53d47f18b7597fa/RF_structure_prediction_weights.pt
'''.split()

def run_main_test_suite(repository_root, working_dir, platform, config):
    full_log = ''

    python_environment = local_python_install(platform, config)

    models_dir = repository_root + '/models'
    if not os.path.isdir(models_dir): os.makedirs(models_dir)

    for url in _models_urls_:
        file_name = models_dir + '/' + url.split('/')[-1]
        tmp_file_name = file_name + '.tmp'
        if not os.path.isfile(file_name):
            print(f'downloading {url}...')
            full_log += f'downloading {url}...\n'
            urllib.request.urlretrieve(url, tmp_file_name)
            os.rename(tmp_file_name, file_name)

    #with tempfile.TemporaryDirectory(dir=working_dir) as tmpdirname:
    tmpdirname = working_dir+'/.ve'
    if True:

        #ve = setup_persistent_python_virtual_environment(python_environment, packages='numpy torch omegaconf scipy opt_einsum dgl')
        #ve = setup_python_virtual_environment(working_dir+'/.ve', python_environment, packages='numpy torch omegaconf scipy opt_einsum dgl e3nn icecream pyrsistent wandb pynvml decorator jedi hydra-core')
        ve = setup_python_virtual_environment(tmpdirname, python_environment, packages='numpy torch omegaconf scipy opt_einsum dgl e3nn icecream pyrsistent wandb pynvml decorator jedi hydra-core')

        execute('Installing local se3-transformer package...', f'cd {repository_root}/env/SE3Transformer && {ve.bin}/pip3 install --editable .')
        execute('Installing RFdiffusion package...', f'cd {repository_root} && {ve.bin}/pip3 install --editable .')

        res, output = execute('running unit tests...', f'{ve.activate} && cd {repository_root} && python -m unittest', return_='tuple', add_message_and_command_line_to_output=True)
        #res, output = execute('running unit tests...', f'cd {repository_root} && {ve.bin}/pytest', return_='tuple')

        results = {
            _StateKey_ : _S_failed_ if res else _S_passed_,
            _LogKey_ : full_log + '\n' + output,
            _ResultsKey_ : {},
        }

    return results



def run(test, repository_root, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    if test == '': return run_main_test_suite(repository_root, working_dir, platform, config)
    else: raise BenchmarkError('Unknow scripts test: {}!'.format(test))