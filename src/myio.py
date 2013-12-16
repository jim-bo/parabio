'''
I/O related code.
'''
import sys
import os
import numpy as np
import scipy.stats as stats
import subprocess
import logging

### classes ###

class Data(object):
    
    def __init__(self, *args):
        ''' constructor '''
        
        # return if no args.
        if len(args) == 0:
            return
            
        # test type of argument.
        if os.path.isfile(args[0]) == True:
            try:
                self._load_txt(args[0])
            except:
                logging.error('error in loading file')
                sys.exit(1)
        elif os.path.isdir(args[0]) == True:
            try:
                self._load_npy(args[0])
            except:
                logging.error('error in loading dir')
                sys.exit(1)        
        else:
            logging.error('unknown input')
            sys.exit(1)
            
        
    def _load_txt(self, file_path):
        """ loads data from txt file """
        
        # load the lines.
        fin = open(file_path, "rb")
        lines = fin.readlines()
        fin.close()
        
        # grab header info.
        tmp = lines[0].strip().split()
        self.nclass = int(tmp[0])
        self.nsamps = int(tmp[1])
        self.nfeats = int(tmp[2])
        
        # create annotation arrays.
        self.class_ids = np.array(lines[1].strip().split())
        self.samps_ids = np.array(lines[2].strip().split())
        self.feats_ids = np.array(lines[3].strip().split())
        
        # create the important arrays.
        self.targets = np.zeros(self.nsamps, dtype=np.int)
        self.expr = np.zeros((self.nsamps, self.nfeats), dtype=np.float)
        self.digi = np.zeros((self.nsamps, self.nfeats), dtype=np.int)
        
        # populate the targets.
        tmp = lines[4].strip().split()
        for i in range(self.nsamps):
            self.targets[i] = int(tmp[i])
        
        # populate the expr matrix.
        i = 0
        for ridx in range(5, 5+self.nsamps):
            tmp = lines[ridx].strip().split()
            for j in range(self.nfeats):
                self.expr[i,j] = float(tmp[j])
            i += 1
            
    def _load_npy(self, in_dir):
        """ loads data from dir """
        
        # load them.
        self.class_ids = np.load("%s/class_ids.npy" % in_dir)
        self.samps_ids = np.load("%s/samps_ids.npy" % in_dir)
        self.feats_ids = np.load("%s/feats_ids.npy" % in_dir)
        self.targets = np.load("%s/targets.npy" % in_dir)
        self.expr = np.load("%s/expr.npy" % in_dir)
        self.digi = np.load("%s/digi.npy" % in_dir)
        
        # compute summary stats.
        self.nclass = self.class_ids.shape[0]
        self.nsamps = self.samps_ids.shape[0]
        self.nfeats = self.feats_ids.shape[0]

    def save_npy(self, out_dir):
        """ save npy dir """
        if os.path.isdir(out_dir) == False:
            subprocess.call(['mkdir','-p',out_dir])
        np.save("%s/class_ids.npy" % out_dir, self.class_ids)
        np.save("%s/samps_ids.npy" % out_dir, self.samps_ids)
        np.save("%s/feats_ids.npy" % out_dir, self.feats_ids)
        np.save("%s/targets.npy" % out_dir, self.targets)
        np.save("%s/expr.npy" % out_dir, self.expr)
        np.save("%s/digi.npy" % out_dir, self.digi)

    def save_txt(self, out_file):
        """ save txt file """
        fout = open(out_file, "wb")
        fout.write("%i %i %i\n" % (self.nclass, self.nsamps, self.nfeats))
        fout.write("%s\n" % " ".join(self.class_ids))
        fout.write("%s\n" % " ".join(self.samps_ids))
        fout.write("%s\n" % " ".join(self.feats_ids))
        fout.write("%s\n" % " ".join([str(x) for x in self.targets]))
        for i in range(self.expr.shape[0]):
            fout.write("%s\n" % " ".join([str(x) for x in self.expr[i]]))
        for i in range(self.digi.shape[0]):
            fout.write("%s\n" % " ".join([str(x) for x in self.digi[i]]))
        fout.close()



### functions ###

def load_redis(ro, data_dir):
    ''' loads data into redis '''
    
    # load the data.
    do = Data(data_dir)
    
    # save the arrays into redis.
    ro.set_numpy('class_ids', do.class_ids)
    ro.set_numpy('samps_ids', do.samps_ids)
    ro.set_numpy('feats_ids', do.feats_ids)
    ro.set_numpy('targets', do.targets)
    ro.set_numpy('expr', do.expr)
    ro.set_numpy('digi', do.digi)

def load_data(fpath):
    ''' loads data into numpy array'''

    # read in data.
    fin = open(fpath)
    lines = fin.readlines()
    fin.close()
    lines.reverse()
    
    # get feature annotation.
    tokens = lines.pop().strip().split()
    features = np.zeros(len(tokens)-3, dtype="S255")    
    for i in range(1,len(tokens)-2):
        features[i-1] = tokens[i]
        
    # reverse back.
    lines.reverse()
        
    # get data and class annotation.
    expr = np.zeros((len(lines),features.shape[0]), dtype=np.float)
    categories = np.zeros(len(lines), dtype="S255")
    quality = np.zeros(len(lines), dtype=np.float)
    for i in range(len(lines)):
        tokens = lines[i].strip().split()
        
        # fill data for this row.
        for j in range(1, len(tokens) - 2):
            expr[i][j-1] = float(tokens[j])

        # fill quality.
        if tokens[-2] == "Y":
            quality[i] = 0.0
        else:
            quality[i] = 1.0
            
        # fill annotation.
        categories[i] = tokens[-1]
        
    # make object and return it.
    return Data(features, categories, quality, expr)

def write_matrix(matrix, file_path):
    ''' writes matrix to file'''
    
    with open(file_path, "wb") as fout:
        for i in range(matrix.shape[0]):
            fout.write('\t'.join([str(x) for x in matrix[i]]) + "\n")
            
def write_vector(vector, file_path):
    ''' writes matrix to file'''
    
    with open(file_path, "wb") as fout:
        fout.write('\t'.join([str(x) for x in vector]) + "\n")

def exec_cmd(cmd):
    """ executes command """
    cmd = [str(x) for x in cmd]
    with open('/dev/null', 'wb') as fout:
        if subprocess.call(cmd, stdout=fout, stderr=fout) != 0:
            logging.error('problem running code')
            logging.error(' '.join(cmd))
            sys.exit(1)     
