#!/usr/bin/python
'''
executes nucmer in parallel using multiprocessing
(future will support MPI and distributed task-queues)
'''
### imports ###

# system
import subprocess
import warnings
import argparse
import logging
import time
import sys
import os
import warnings
import itertools
import multiprocessing
import shutil
import traceback

# special
from pyfasta import Fasta, MemoryRecord
import pyfasta.split_fasta

# logging level.
#logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(message)s', )
logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s', )

# local
#import src.run

# hack to silence argparser.
warnings.filterwarnings('ignore', category=DeprecationWarning)

### configuration ###

### definitions ###

### classes ###      

class Consumer(multiprocessing.Process):
    
    def __init__(self, task_queue, result_queue):
        ''' create task and result ququeue'''
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue

    def run(self):
        ''' run each task '''
        
        # loop until poison.
        proc_name = self.name
        while True:
            
            # get next task.
            next_task = self.task_queue.get()
            
            # poison
            if next_task is None:
                self.task_queue.task_done()
                break
                
            # run the task.
            answer = next_task()
            self.task_queue.task_done()
            
            # return result.
            self.result_queue.put(answer)
    
class TaskNucmer(object):
    def __init__(self, cwd, prefix, ref, qry, nuc_opts, del_opts, sho_opts):
        self.cwd = cwd
        self.prefix = prefix
        self.ref = ref
        self.qry = qry
        self.nuc_opts = nuc_opts
        self.del_opts = del_opts
        self.sho_opts = sho_opts
        
    def __call__(self):
        logging.debug("starting nucmer: %s" % self.cwd)
        
        # open log and error files.
        log_file = "%s/paranuc.out" % self.cwd
        
        # create nucmer commands.
        cmd1, cmd2, cmd3 = _commands(self.ref, self.qry, self.nuc_opts, self.del_opts, self.sho_opts)
        
        # fill a script with them.
        txt = "#!/bin/bash\n"
        txt += ' '.join([str(x) for x in cmd1]) + ' || exit 1\n'
        if self.del_opts != "":
            txt += ' '.join([str(x) for x in cmd2]) + ' || exit 1\n'
        if self.sho_opts != "":
            txt += ' '.join([str(x) for x in cmd3]) + ' || exit 1\n'

        # try to write script. 
        try:
            script_file = "%s/nucmer.sh" % self.cwd
            with open(script_file, "wb") as fout:
                fout.write(txt)
            subprocess.call(["chmod", "a+x", script_file])
            
            # try to execute.
            try:
                with open(log_file, 'w') as log:
                    ret = subprocess.call([script_file],cwd=self.cwd,stdout=log,stderr=log)
            except:
                msg = traceback.format_exc()
                ret = 1                 
                return 1, msg
        except:
            msg = traceback.format_exc()
            ret = 1          
            return 1, msg  
        
        # return error code and log-file.
        return ret, log_file
        
    def __str__(self):
        return str(prefix)
    
class NucObj(object):
    ''' abstract Nucmer alignment experiment '''
    
    # required files
    out_dir = ""
    ref_fasta = ""
    qry_fasta = ""
    dlt_file = ""
    ref_names = []
    ref_dirs = []
    
    def __init__(self, ref_fasta, qry_fasta, out_dir):
        ''' files for nucmer alignment '''

        # sanity input.
        if out_dir == None or ref_fasta == None or qry_fasta == None:
            logging.error("missing required input")
            sys.exit(1)
        
        for x in [ref_fasta, qry_fasta]:
            if os.path.isfile(x) == False:
                logging.error("bad file: %s" % x)
                sys.exit(1)
                
        # prepare output base.
        for x in [out_dir]:
            if os.path.isdir(x) == False:
                logging.info("creating output directory: %s" % x)
                os.makedirs(x)
       
        # update hard-coded links.
        self.out_dir = os.path.abspath(out_dir)
        self.ref_fasta = os.path.abspath(ref_fasta)
        self.qry_fasta = os.path.abspath(qry_fasta)
        self.dlt_file = '%s/nucmer.delta' % self.out_dir
       
    def _no_empty(self, lista, listb):
        ''' removes empty entries '''
        
        # check for empty fasta.
        tmpa = list()
        tmpb = list()
        for i in range(len(listb)):
            
            # open it.
            try:
                z = Fasta(listb[i], record_class=MemoryRecord)
            
                # check for empty.
                if len(z.keys()) == 0:
                    continue

                # add to temp.
                tmpa.append(lista[i])
                tmpb.append(listb[i])

            except:
                logging.warning("bad fasta file")
            
        # sort back.
        return tmpa, tmpb

       
    def split_seqs(self, num_jobs, max_ref=5, max_qry=20):
        ''' splits reference and query into appropriate number of splits '''
        
        # load data into memory.
        r = Fasta(self.ref_fasta, record_class=MemoryRecord)
        q = Fasta(self.qry_fasta, record_class=MemoryRecord)
        
        ## reference ##
        # split according to criteria.
        if len(r) < max_ref:
            max_ref = len(r)
            
        if max_ref > num_jobs:
            max_ref = 1
        
        if len(q) < max_qry:
            max_qry = len(q)

        if num_jobs < max_qry:
            max_qry = num_jobs

        if (max_ref * max_qry) > num_jobs:
            max_qry = int(float(num_jobs) / float(max_ref))
        
        # count number of seqs.
        sc = len(r.keys())
        
        # create split info.
        self.ref_names = ["ref_%i" % x for x in range(max_ref)]
        self.ref_files = ["%s/%s.fasta" % (self.out_dir, x) for x in self.ref_names]
        
        # split according to rules.
        pyfasta.split_fasta.without_kmers(r, self.ref_files)
        self.ref_names, self.ref_files = self._no_empty(self.ref_names, self.ref_files)
        
        ## query ##
        # create split info.
        self.qry_names = ["qry_%i" % x for x in range(max_qry)]
        self.qry_files = ["%s/%s.fasta" % (self.out_dir, x) for x in self.qry_names]
        
        # split according to rules.
        pyfasta.split_fasta.without_kmers(q, self.qry_files)
        self.qry_names, self.qry_files = self._no_empty(self.qry_names, self.qry_files)
                
### private functions ###
def _setup_multi(num_cpu):
    
    # Establish communication queues
    tasks = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()

    # Start consumers
    logging.info('creating %d consumers' % num_cpu)
    consumers = [ Consumer(tasks, results) for i in xrange(num_cpu) ]
    for w in consumers:
        w.start()

    return tasks, results

def _kill_multi(tasks, num_cpu):
    for i in xrange(num_cpu):
        logging.debug("waiting on: %i" % i)
        tasks.put(None)
    tasks.join()    
    
def _combine_delta(suffix, out_file, flist, nucobj, args):
    
    # open all delta files.
    fout = open(out_file, "wb")
    fout.write("%s %s\n" % (nucobj.ref_fasta, nucobj.qry_fasta))
    fout.write("NUCMER\n")
    for rn, rf, qn, qf, cwd, prefix in flist:
                
        # verify non-empty file.
        dlt_tmp = '%s/out.%s' % (cwd, suffix)
        if os.path.isfile(dlt_tmp) == False or os.path.getsize(dlt_tmp) == 0:
            logging.error("error in nucmer: %s" % dlt_tmp)
            sys.exit()
        
        # open file and skip first two lines.
        with open(dlt_tmp, "rb") as fin:
            
            # skip.
            fin.readline()
            fin.readline()
            
            # copy remainder.
            for line in fin:
                fout.write(line)
                
    # close delta file.
    fout.close()
    
def _combine_coords(out_file, flist, nucobj, args):
    
    # open all delta files.
    tmp_file = out_file + ".tmp"
    fout = open(tmp_file, "wb")
    fout.write("%s %s\n" % (nucobj.ref_fasta, nucobj.qry_fasta))
    fout.write("NUCMER\n")
    for rn, rf, qn, qf, cwd, prefix in flist:
                
        # verify non-empty file.
        dlt_tmp = '%s/out.coords' % (cwd)
        if os.path.isfile(dlt_tmp) == False or os.path.getsize(dlt_tmp) == 0:
            logging.error("error in nucmer: %s" % dlt_tmp)
            sys.exit()
        
        # open file and skip first two lines.
        with open(dlt_tmp, "rb") as fin:
            
            # skip.
            fin.readline()
            fin.readline()
            
            # copy remainder.
            for line in fin:
                fout.write(line)
                
    # close delta file.
    fout.close()
    
    # sort results.
    script_file = "%s.sh" % out_file
    cmd = ['cat', tmp_file, '|', 'sort', '-k13', '-k1n', '-k2n', '>',out_file]
    with open(script_file, "wb") as fout:
        fout.write("#!/bin/bash\n")
        fout.write(' '.join(cmd) + '\n')
        fout.write('rm -f %s %s\n' % (tmp_file, script_file))
    subprocess.call(['chmod','u+x',script_file])
    subprocess.call([script_file])
    
def _commands(ref, qry, nuc_opts, del_opts, sho_opts):
    cmd1 = ['nucmer'] + [x for x in nuc_opts.split()] + [ref, qry]
    if del_opts != "":
        cmd2 = ['delta-filter'] + [x for x in del_opts.split()] + ["out.delta", ">", "out.fdelta"]
    if sho_opts != "":
        cmd3 = ['show-coords'] + [x for x in sho_opts.split()] + ['out.fdelta', ">", "out.coords"]
    return cmd1, cmd2, cmd3
    
### callable functions ###
    
def smp(args):
    ''' runs nucmer '''
    
    # create experiment object.
    nucobj = NucObj(args.ref_fasta, args.qry_fasta, args.tmp_dir)
    
    # create the appropriate split.
    nucobj.split_seqs(args.num_cpu, max_ref=args.max_ref, max_qry=args.max_qry)
    
    # setup processing.
    tasks, results = _setup_multi(args.num_cpu)
    
    # run over each combo and run nucmer.
    flist = list()
    a = zip(nucobj.ref_names, nucobj.ref_files)
    b = zip(nucobj.qry_names, nucobj.qry_files)
    for x, y in itertools.product(a, b):
        
        # simplify.
        rn, rf = x
        qn, qf = y
        
        # create directory.
        cwd = "%s/%s_%s" % (nucobj.out_dir, rn, qn)
        prefix = "%s/out" % cwd
        if os.path.isdir(cwd) == False:
            logging.debug("creating internal directory: %s" % cwd)
            os.makedirs(cwd)     
               
        # enqueue task.
        tasks.put(TaskNucmer(cwd, prefix, rf, qf,\
            args.nuc_opts.strip('"'),\
            args.del_opts.strip('"'),\
            args.sho_opts.strip('"'))\
        )
        
        # save for later use.
        flist.append((rn,rf,qn,qf,cwd,prefix))
            
    # wait for all tasks to finish and kill workers.
    logging.info("waiting for success")    
    _kill_multi(tasks, args.num_cpu)
    
    # collecting tasks.
    error = False
    for i in range(len(flist)):
        r, msg = results.get()
        if r != 0:
            logging.error('non-zero exit: %s' % msg)
            error = True
    if error == True:
        sys.exit(1)
    
    # combine results.
    logging.info("combining results")
    res1 = "%s/out.delta" % args.nuc_dir
    res2 = None
    res3 = None
    
    _combine_delta("delta", res1, flist, nucobj, args)
    if args.del_opts.strip('"') != "":
        res2 = "%s/out.fdelta" % args.nuc_dir
        _combine_delta("fdelta", res2, flist, nucobj, args)
    if args.sho_opts.strip('"') != "":
        res3 = "%s/out.coords" % args.nuc_dir
        _combine_coords(res3, flist, nucobj, args)

    # finish.
    logging.info("completed succesfully")
    return res1, res2, res3
 

def cluster(args):
    ''' runs nucmer on cluster'''
    ### SAME ###
    # create experiment object.
    nucobj = NucObj(args.ref_fasta, args.qry_fasta, args.tmp_dir)
    
    # create the appropriate split.
    nucobj.split_seqs(args.num_cpu, max_ref=args.max_ref, max_qry=args.max_qry)
    ### SAME ###
    
    # cluster details.
    name = "%s_%s" % (args.ref_fasta.split("/")[-1], args.qry_fasta.split("/")[-1])
    bdir = '/home/CAM/jlindsay/projects/nucmer/%s' % args.prj_name
    acmd = 'jlindsay@sig1-submit-ext.cam.uchc.edu'
    
    # create master script.
    llog_file = '%s/paranuc.log' % nucobj.out_dir
    lmaster = '%s/master.sh' % nucobj.out_dir
    rmaster = '%s/master.sh' % bdir
    
    # open log.
    llog = open(llog_file, "wb")
    
    # write master header.
    mout = open(lmaster, 'wb')
    mout.write("#!/bin/bash\n")
    
    # clear cluster and make dir.
    logging.info("making remote directories")
    #subprocess.call(["ssh %s 'mkdir -p %s'" % (acmd, bdir)], shell=True, stdout=llog, stderr=llog)
    subprocess.call(["ssh %s 'rm -rf %s'" % (acmd, bdir)], shell=True, stdout=llog, stderr=llog)
    #subprocess.call(["ssh %s 'mkdir -p %s'" % (acmd, bdir)], shell=True, stdout=llog, stderr=llog)
    
    # build script.
    ulist = [lmaster]
    flist = list()
    a = zip(nucobj.ref_names, nucobj.ref_files)
    b = zip(nucobj.qry_names, nucobj.qry_files)
    todo = 0
    for x, y in itertools.product(a, b):
    
        # simplify.
        rn, lrf = x
        qn, lqf = y
        lcwd = "%s/%s_%s" % (nucobj.out_dir, rn, qn)
        lprefix = "%s/nucmer" % lcwd
        lscript = '%s.sh' % (lprefix)
        
        rrf = '%s/%s' % (bdir, lrf.split("/")[-1])
        rqf = '%s/%s' % (bdir, lqf.split("/")[-1])
        rcwd = "%s/%s_%s" % (bdir, rn, qn)
        rprefix = "%s/nucmer" % rcwd
        rscript = '%s.sh' % (rprefix)
        
        # save for later use.
        flist.append((rn,lrf,qn,lqf,lcwd,lprefix))
        ulist.append(lcwd)
        ulist.append(lrf)
        ulist.append(lqf)
        ulist.append(lscript)
        
        # create local directory.
        if os.path.isdir(lcwd) == False:
            os.makedirs(lcwd)     
        
        # create command.
        cmd1, cmd2, cmd3 = _commands(rrf, rqf, args.nuc_opts.strip('"'), args.del_opts.strip('"'), args.sho_opts.strip('"'))
        cmd1[0] = cmd1[0].replace("nucmer","/home/CAM/jlindsay/.local/bin/nucmer")
        cmd2[0] = cmd2[0].replace("delta-filter","/home/CAM/jlindsay/.local/bin/delta-filter")
        cmd3[0] = cmd3[0].replace("show-coords","/home/CAM/jlindsay/.local/bin/show-coords")
        
        txt = ""
        txt += ' '.join([str(x) for x in cmd1]) + ' || exit 1\n'
        if args.del_opts != "":
            txt += ' '.join([str(x) for x in cmd2]) + ' || exit 1\n'
        if args.sho_opts != "":
            txt += ' '.join([str(x) for x in cmd3]) + ' || exit 1\n'   
            
        # write to script.
        with open(lscript, "wb") as fout:
            fout.write("#!/bin/bash\n")
            fout.write("#PBS -o %s.qsub\n" % rprefix)
            fout.write("#PBS -j oe\n")
            fout.write("#PBS -l pmem=10gb\n")
            fout.write("cd %s\n" % rcwd)
            fout.write(txt)
            
        # add line to master script.
        mout.write("qsub %s\n" % (rscript))
        
    # wait for jerbs to finish then zip up results.
    mout.write("/home/CAM/jlindsay/clearq.sh 30\n")
    mout.write("tar cvfz results.tar.gz ref_*/out.*\n")
    mout.close()
    
    # create file list relative tmp directory.
    file_list = "%s/file_list.txt" % args.tmp_dir
    with open(file_list, "wb") as fout:
        for f in ulist:
            f = os.path.relpath(f).replace(args.tmp_dir,"")
            fout.write(f + '\n')
    
    # copy using rsync.
    logging.info("uploading data to cluster")
    cmd = ['rsync','-avz', '--files-from=file_list.txt', '-e','ssh','./','%s:%s' % (acmd, bdir)]
    subprocess.call(cmd, cwd=args.tmp_dir)
    
    # write execution script locally.
    escript = '%s/cmds.sh' % args.tmp_dir
    with open(escript, 'wb') as fout:
        fout.write('cd %s;\n' % bdir)
        fout.write('sh master.sh\n')
        
    # run it remotely.
    logging.info("running on cluster")
    subprocess.call(["cat %s | ssh %s" % (escript,acmd)], shell=True, stdout=llog, stderr=llog)
    
    # retrieve results.
    logging.info("retrieving results")
    subprocess.call(["scp %s:%s/results.tar.gz %s/results.tar.gz" % (acmd, bdir, args.tmp_dir)], shell=True)
    subprocess.call(["cd %s; tar xvfz results.tar.gz" % (args.tmp_dir)], shell=True, stdout=llog, stderr=llog)
    
    # close log.
    llog.close()
    
    ### SAME ###
    # combine results.
    logging.info("combining results")
    res1 = "%s/out.delta" % args.nuc_dir
    res2 = None
    res3 = None
    
    _combine_delta("delta", res1, flist, nucobj, args)
    if args.del_opts.strip('"') != "":
        res2 = "%s/out.fdelta" % args.nuc_dir
        _combine_delta("fdelta", res2, flist, nucobj, args)
    if args.sho_opts.strip('"') != "":
        res3 = "%s/out.coords" % args.nuc_dir
        _combine_coords(res3, flist, nucobj, args)

    # finish.
    logging.info("completed succesfully")
    return res1, res2, res3
    ### SAME ###
        
### script ###

def _add_options(subpp):
    subpp.add_argument('-nucmer', dest='nuc_opts', type=str, required=True, help='nucmer: use quotes to include all nucmer options excluding reference, query and output dir')
    subpp.add_argument('-delta', dest='del_opts', type=str, default="", help='delta filter: use quotes to include all nucmer options excluding reference, query and output dir')
    subpp.add_argument('-show', dest='sho_opts', type=str, default="", help='show-coords: use quotes to include all nucmer options excluding reference, query and output dir')
    subpp.add_argument('-p', dest='num_cpu', type=int, required=True, help='number of jobs to use')
    subpp.add_argument('-maxref', dest='max_ref', type=int, default=10, help='max number of ref splits')
    subpp.add_argument('-maxqry', dest='max_qry', type=int, default=20, help='max number of qry splits')
    subpp.add_argument('-r', dest='ref_fasta', type=str, required=True, help='reference fasta')
    subpp.add_argument('-q', dest='qry_fasta', type=str, required=True, help='query fasta')
    subpp.add_argument('-d', dest='nuc_dir', type=str, required=True, help='nucmer dir')
    subpp.add_argument('-t', dest='tmp_dir', type=str, required=True, help='output dir')

if __name__ == '__main__':

    # mode parser.
    main_p = argparse.ArgumentParser()
    subp = main_p.add_subparsers(help='sub-command help')

    # SMP computer.
    subpp = subp.add_parser('smp', help='runs nucmer using multiprocess')
    _add_options(subpp)
    subpp.set_defaults(func=smp)

    # cluster
    subpp = subp.add_parser('cluster', help='runs nucmer using multiprocess')
    _add_options(subpp)
    subpp.add_argument('-n', dest='prj_name', type=str, required=True, help='unique project name on cluster')
    subpp.set_defaults(func=cluster)
 
    # parse args.
    args = main_p.parse_args()
    args.func(args)
