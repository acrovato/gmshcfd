#!/usr/bin/env python3
# -*- coding: utf8 -*-
# test encoding: à-é-è-ô-ï-€

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

# Run a script as if the software was installed

def parse_args():
    """Parse command line arguments
    """
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--clean', help='clean workspace', action='store_true')
    parser.add_argument('file', help='python files')
    args = parser.parse_args()
    return args

def setup_workdir(testname, clean, verb=True):
    """Create a single directory for the given test
    """
    import os, os.path

    # build the name of the workspace folder
    dir1 = os.path.dirname(os.path.abspath(os.path.dirname(__file__))) + os.sep
    if verb: print(f'Setting workspace for "{testname}"')
    common = os.path.commonprefix((testname, dir1))
    resdir = testname[len(common):].replace(os.sep, "_")
    resdir = os.path.splitext(resdir)[0] # remove ".py"
    wdir = os.path.join('workspace', resdir)

    # create the directory
    if not os.path.isdir(wdir):
        if verb: print('- creating', wdir)
        os.makedirs(wdir)
    elif os.path.isdir(wdir) and clean:
        if verb: print('- cleaning', wdir)
        import shutil
        for f in os.listdir(wdir):
            fpth = os.path.join(wdir, f)
            if os.path.isfile(fpth):
                os.remove(fpth)
            elif os.path.isdir(fpth):
                shutil.rmtree(fpth)
    # change dir
    if verb: print('- changing to', wdir)
    os.chdir(wdir)

def main():
    import os, time, socket
    import gmshcfd

    # prepare run
    args = parse_args()
    testname = os.path.abspath(args.file)
    if not os.path.isfile(testname):
        raise Exception(f'File not found: {testname}')
    setup_workdir(testname, args.clean)

    # Run
    print('Starting test', testname)
    print('Time:', time.strftime('%c'))
    print('Hostname:', socket.gethostname())
    script = open(testname, 'r', encoding='utf-8').read()
    exec(compile(script, testname, 'exec'), {'__file__': testname, '__name__':'__main__'})

if __name__ == '__main__':
    main()
