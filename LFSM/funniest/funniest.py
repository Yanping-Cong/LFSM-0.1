from markdown import markdown
import h5py
import os
class funnist(object):
    def __init__(self):
        _path = os.path.dirname(os.path.abspath(__file__))
        _path = os.path.split(_path)[0]
        self.file_dir = _path
    def read_file(self):
        with h5py.File(self.file_dir+'/input/sky_map_10-2300MHz/drao_10MHz.hdf5','r') as f:
            print (f.keys())
    def read_file1(self):
        with h5py.File('../input/sky_map_10-2300MHz/drao_10MHz.hdf5','r') as f:
            print (f.keys())
    def read_file2(self):
        with h5py.File('sky_map_10-2300MHz/drao_10MHz.hdf5','r') as f:
            print (f.keys())

    def joke(self):
        return markdown(u'Wenn ist das Nunst\u00fcck git und Slotermeyer?'
                    u'Ja! ... **Beiherhund** das Oder die Flipperwaldt '
                    u'gersput.')
