# -*- coding: utf-8 -*-

# GmshCFD
# Copyright (C) 2023 Adrien Crovato
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

def _parse_args():
    """Parse command line arguments
    """
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--clean', help='clean workspace', action='store_true')
    parser.add_argument('file', help='python file')
    args = parser.parse_args()
    return args

def _setup_workdir(fname, clean, verb=True):
    """Create a single directory for the given test
    """
    import os, os.path
    # Build the name of the workspace folder
    ori_dir = os.path.abspath(os.getcwd()) + os.path.sep
    if verb: print(f'Setting run directory for "{fname}"')
    common = os.path.commonprefix((fname, ori_dir))
    resdir = fname[len(common):].replace(os.path.sep, '_')
    resdir = os.path.splitext(resdir)[0] # remove ".py"
    new_dir = os.path.join('workspace', resdir)

    # Create the directory and clean it
    if not os.path.isdir(new_dir):
        if verb: print('- creating', new_dir)
        os.makedirs(new_dir)
    elif os.path.isdir(new_dir) and clean:
        if verb: print('- cleaning', new_dir)
        import shutil
        for f in os.listdir(new_dir):
            fpth = os.path.join(new_dir, f)
            if os.path.isfile(fpth):
                os.remove(fpth)
            elif os.path.isdir(fpth):
                shutil.rmtree(fpth)

    # Change directory
    if verb: print('- changing to', new_dir)
    os.chdir(new_dir)

def main():
    """Entry point for gmshcfd-run"""
    import os, time, socket

    # Get file
    args = _parse_args()
    fname = os.path.abspath(args.file)
    if not os.path.isfile(fname):
        raise FileNotFoundError(f'File not found: {fname}')

    # Create workspace and change to it
    _setup_workdir(fname, args.clean)

    # Run file
    print(f'[{socket.gethostname()} - {time.strftime("%c")}] Running:', fname)
    script = open(fname, 'r', encoding='utf-8').read()
    exec(compile(script, fname, 'exec'), {'__file__': fname, '__name__':'__main__'})
