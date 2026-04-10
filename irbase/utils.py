import gzip
import bz2
import itertools
import sys
from pathlib import PurePath, Path

def open_compressed(filename, mode='rb'):
    if filename == '-':
        if 'r' in mode:
            if 'b' in mode:
                return sys.stdin.buffer
            else:
                return sys.stdin
        elif 'w' in mode:
            if 'b' in mode:
                return sys.stdout.buffer
            else:
                return sys.stdout
    else:
        if not isinstance(filename, PurePath):
            filename = Path(filename)
        
        if filename.suffix == '.gz':
            return gzip.open(filename, mode)
        elif filename.suffix == '.bz2':
            return bz2.open(filename, mode)     
        else:
            return open(filename, mode=mode)

def slice_from_range(range_):
    return slice(range_['start'], range_['stop'])

def make_range(start, stop):
    result = {'start': start, 'stop': stop}
    return result

def make_named_range(name, start, stop):
    result = {'name': name, 'start': start, 'stop': stop}
    return result

def make_tail_range(start, stop, length):
    stop = stop - length
    if stop == 0:
        stop = None
    return make_range(start, stop)
