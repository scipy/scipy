""" Cow --- Cluster Of Workstations.
    These modules facilitate the use of a cluster
    of workstations as a single computational engine.
    It will hopefully simplify parallelizing 
    "embarrasingly parallel" task.  If your problem 
    requires closely coupled parallelism, take a 
    look at pyMPI.
    
    Note: we might want to move the server list out
    of the module where it can be passed as an argument.
"""

import sync_cluster, socket
import os # for getuid()

ClusterError = 'ClusterError'
TimeoutError = 'TimeoutError'
      
class machine_cluster:
    def __init__(self,server_list):        
        self.workers=[]
        self.worker_by_name={}
        worker_id = 1
        for host,port in server_list:
            #add the uid of the person starting the job to the port to
            #hopefully avoid port conflicts.  This is only guaranteed to
            #work on a cluster that shares users.
            port = port #+ 760 #+ os.getuid() #304 (balaji's id) 
            new_worker = sync_cluster.standard_sync_client(host,port,worker_id)
            self.workers.append(new_worker)
            self.worker_by_name[host] = new_worker
            worker_id = worker_id + 1
                                   
    def start(self,force_restart=0,timeout=60):
        """ Start all the worker processes on remote machines.
            The timeout value is specified in seconds and
            defaults to 60.  The timeout starts counting down
            only after rsh has tried to start all the remote
            processes.  This means the actual time in the
            function could be much longer than 60 seconds - depending
            on how long rsh takes.
            
            Its possible the 60 second time out will be
            to short for large clusters - but I hope not!
        """    
        if not force_restart and self.is_running():
            return
        # start the worker processes.
        for worker in self.workers:
            worker.start_server()
        if not self.is_running(timeout):
            raise TimeoutError  
            
    def stop(self):
        for worker in self.workers:
            import sys; sys.stdout.flush()
            try: worker.exec_code('import sys;sys.exit(0)')
            except:
                #should really do something here to
                # trap non-SystemExit errors.  
                pass
    def restart(self):
        self.stop()
        self.start()

    def is_running(self,timeout=0):
        # wait for them to start           
        import time
        st = time.time()            
        still_waiting = 1
        while still_waiting: 
            try:
                # Send a simple command to all workers
                # and wait till they handle it successfully
                self.exec_code("1==1") 
            except ClusterError:
                still_waiting = 1
                elapsed = time.time() - st
                if elapsed > timeout:
                    # We've run out of time.
                    return 0
            else:
                still_waiting = 0    
            wait_time = time.time() - st                                
        # should we somehow dessiminate worker topology (ids)
        # to all machines here?
        return 1
        
    def _send(self,package,addendum=None):
        """addendum is either None, or a list
           of addendums <= in length to
           the number of workers
        """
        
        if addendum: 
            N = len(addendum)
            assert(N <= len(self.workers))
        else:
            N = len(self.workers)            
        self.send_exc = {}
        self.had_send_error = []
        for i in range(N):
            try:
                if not addendum:
                    self.workers[i].send(package)
                else:
                    self.workers[i].send(package,addendum[i])   
            except socket.error, msg:    
                import sys
                err_type, err_msg = sys.exc_info()[1]
                self.had_send_error.append(self.workers[i])
                try: self.send_exc[(err_type,err_msg)].append(self.workers[i].id)
                except: self.send_exc[(err_type,err_msg)] = [self.workers[i].id]
            # else - handle other errors?               
        self.Nsent = N            
    def _recv(self):    
        self.had_recv_error = []
        self.recv_exc = {}        
        results = []
        import sys; 
        #only listen on workers involved in calculation.
        for worker in self.workers[:self.Nsent]:
            if worker in self.had_send_error:
                results.append(None)
            else:    
                try:
                    print worker.id,
                    sys.stdout.flush()
                    results.append(worker.recv())
                                        
                except sync_cluster.RemoteError:    
                    import sys
                    err = sys.exc_info()[1]
                    # Force the err msg (err[1]) to be a string.
                    # This dimishes info content, but makes sure
                    # that the sames erros are hashed correctly
                    # in the dictionary.  
                    err_type,err_msg, err_traceback = err
                    err = err_type,str(err_msg), err_traceback
                    self.had_recv_error.append(worker)                    
                    try: self.recv_exc[err].append(worker.id)
                    except: self.recv_exc[err] = [worker.id]
                    results.append(None)
                except sync_cluster.RemoteCrashError:
                    # Gotta be more intelligent here...
                    msg =  'Error! Remote worker %d appears to have crashed.' \
                            % worker.id
                    raise sync_cluster.RemoteCrashError,msg           
                # else handle other errors    
        print                
        return tuple(results)                
    
    def _send_recv(self,package,addendum=None):
        self._send(package,addendum)
        self.last_results = self._recv()
        if(len(self.send_exc) or len(self.recv_exc)):
            self.handle_error()
        return self.last_results    

    def handle_error(self):
        # perhaps do some nifty stuff here to 
        # mark bad workers, try to restart, etc.
        msg = ''
        Nworkers = len(self.workers)
        Nsend_errors = len(self.had_send_error)
        Nsend_error_types = len(self.send_exc.keys())
        Nrecv_errors = len(self.had_recv_error)
        Nrecv_error_types = len(self.recv_exc.keys())        
        if (Nsend_errors == Nworkers and
            Nsend_error_types == 1):    
            sock_err_type,err_msg = self.send_exc.keys()[0]
            if sock_err_type == 111:
                # An attempt at helpful info for a common problem.
                msg =       '\n\nConnection refused on all workers.\n'
                msg = msg + '    Perhaps restarting the cluster would help.\n'
                msg = msg + '    Use Your_Clusters_Name_Here.restart()'
            else:
                msg =       'A Socket error occured sending to all workers.\n\t'
                msg = msg + str(sock_err_type) + ': ' + str(err_msg)                
        elif Nsend_errors:
            msg = '\n\nThe following errors occured when sending data:\n\t'
            for err,guilty_workers in self.send_exc.items():                
                msg = msg + str(err) + '\n\t'
                msg = msg + 'Guilty workers: ' + str(guilty_workers) + '\n'    
        
        if (Nrecv_errors == Nworkers and
              Nrecv_error_types == 1):
                err,dummy = self.recv_exc.items()[0]                    
                err_type, err_msg, err_traceback = err
                msg =       '\n\nThe same error occured on all workers:\n\t'
                msg = msg + str(err_type) + ': ' + str(err_msg)                
                msg = msg + err_traceback
        elif Nrecv_errors:
            msg = '\n\nThe following errors occured on workers:\n\t'
            for err,guilty_workers in self.recv_exc.items():                
                err_type, err_msg, err_traceback = err
                msg = msg + str(err_type) + ': ' + str(err_msg) + '\n'              
                msg = msg + 'Guilty workers: ' + str(guilty_workers) + '\n'    
                msg = msg + err_traceback
                

        raise ClusterError, msg
        
    def load(self):        
        import scipy.common.proc
        res = self.apply(scipy.common.proc.load_avg,())
        return res

    def info(self):
        import scipy.common.proc
        res = self.apply(scipy.common.proc.machine_info,())
        return res

    def load_summary(self):        
        import string
        import scipy.common.proc
        results = self.load()
        for i in range(len(self.workers)):            
            name = string.split(self.workers[i].host,'.')[0]
            res = results[i]
            s = "%6s: %1.2f," % (name[-6:], res['load_1'])
            print s,
            if not ((i+1) % 5):
                print
    
    def info_summary(self):
        import string
        results = self.info()
        labels = "%-8s  %-9s  %-4s  %-8s  %-8s  %-4s" % \
                 ('MACHINE','CPU','GHZ','MB TOTAL', 
                  'MB FREE','LOAD')
        print labels          
        for i in range(len(self.workers)):            
            name = string.split(self.workers[i].host,'.')[0]
            res = results[i]
            s = "%-8s %2dx%-6s  %4.1f  %8.1f  %8.1f   %4.2f" %  \
                (name[-8:], res['cpu_count'],res['cpu_type'][-6:], \
                 res['cpu_speed'],res['mem_total'],res['mem_free'],\
                 res['load_1'])
            print s
    
    def ps_list(self,sort_by='cpu',**filters):
        import operator
        import scipy.common.proc
        res = self.apply(scipy.common.proc.ps_list,())
        psl = reduce(operator.add,res)
        psl = scipy.common.proc.ps_sort(psl,sort_by,**filters)        
        return psl
 
    def ps(self,sort_by='cpu',**filters):
        psl = self.ps_list(sort_by,**filters)
        if len(psl):
            print psl[0].labels_with_name()
        for i in psl: print i.str_with_name()

    def nice(self,increment=10):
        """* increment the interpreter daemon on all remote machines.
             hmmm. this doesn't seem to work. 
         see os.nice()
        *"""
        res = self.apply(os.nice,(increment,))
        return res

    def renice(self,process_list,level):
        """* change the nice level of multiple remote processes.
             
             Once niced down, a process cannot be reniced back up.
             This is a Linux issue.
        *"""    
        res = []
        pids = {}
        for process in process_list:
            if hasattr(process,'machine'):
                try:
                    worker = self.worker_by_name[process.machine] 
                except KeyError:
                    worker = self.worker_by_name[process.long_machine] 
                pid = process.pid
            else:
                worker = self.workers[process[0]] 
                pid = process[1]
            try:
                pids[worker] = pids[worker] + ' ' + str(pid)
            except:
                pids[worker] = str(pid)
        for worker,value in pids.items():
            arg = 'renice %d -p %s' % (level,value)
            res.append(worker.apply(os.system,(arg,)))
        return res

    def kill(self,process_list,signal = 'TERM'):
        """* change the nice level of multiple remote processes.
             
             Once niced down, a process cannot be reniced back up.
             This is a Linux issue.
        *"""    
        res = []
        pids = {}
        for process in process_list:
            if hasattr(process,'machine'):
                try:
                    worker = self.worker_by_name[process.machine] 
                except KeyError:
                    worker = self.worker_by_name[process.long_machine] 
                pid = process.pid
            else:
                worker = self.workers[process[0]] 
                pid = process[1]
            try:
                pids[worker] = pids[worker] + ' ' + str(pid)
            except:
                pids[worker] = str(pid)
        for worker,value in pids.items():
            arg = 'kill -s ' + signal + ' %s' % (level,value)
            res.append(worker.apply(os.system,(arg,)))
        return res
    
    #def system(self,cmd):
    #    return self.workers.apply(os.system,(cmd,))
    
    def system(self,cmd):
        return self.exec_code(('import os;f=os.popen("%s");res = f.read(-1);f.close();' \
               % cmd),returns=['res'])
    def reload(self,module):
        try: 
            cmd = 'import %s; reload(%s)' % ((module.__name__,)*2)
        except AttributeError:
            cmd = 'import %s; reload(%s)' % ((module,)*2)
        self.workers.exec_code('cmd')    
        
# mirror all of sync_client functions
    # They assumes all clients have the same packing procedures.
    def exec_code(self,code,inputs=None,returns=None,global_vars=None):
        #use the first worker as a server
        package = self.workers[0].exec_code_pack(code,inputs,returns,global_vars)
        return self._send_recv(package)

    def apply(self,function,args=(),keywords=None):        
        package = self.workers[0].apply_pack(function,args,keywords)
        return self._send_recv(package)

    def loop_apply(self,function,loop_var,args,keywords=None):
        #----------------------------------------------------
        # Prepare main package for sending
        # almost verbatim from loop_apply_pack in sync_cluster
        #----------------------------------------------------
        #if type(loop_var) == type(1):
        #    loop_var = function.func_code.co_varnames[loop_var]
        all_keywords = {}
        if keywords: all_keywords.update(keywords)
        #more_keywords = sync_cluster.args_to_keywords(function,args)
        #sync_cluster.catch_keyword_conflicts(more_keywords,all_keywords)
        #all_keywords.update(more_keywords)

        # pull out the loop variable.
        if type(loop_var) != type(''):
            loop_var = int(loop_var)
            loop_data = args[loop_var]
            # no need to pack and send since it'll be in the "addendum"
            args = list(args)
            args[loop_var] = None
            args = tuple(args)
        else:
            loop_data = all_keywords[loop_var]    
            # no need to pack and send since it'll be in the "addendum"
            del all_keywords[loop_var]             
        contents={'_command':sync_cluster.loop_func,'function':function,
                  'args':args,'keywords':all_keywords,'loop_var':loop_var}
        package = self.workers[0].packer.pack(contents) 
        return self.loop_send_recv(package,loop_data,loop_var)

    def loop_code(self,code,loop_var,inputs=None,returns=None,global_vars=None):
        the_inputs = {}
        the_inputs.update(inputs)    
        loop_data = the_inputs[loop_var]
        the_inputs[loop_var] = None #make it small for packing
        package = self.workers[0].loop_code_pack(code,loop_var,
                        the_inputs,returns,global_vars)    
        return self.loop_send_recv(package,loop_data,loop_var)
        
    def loop_send_recv(self,package,loop_data,loop_var):        
        #----------------------------------------------------
        # Now split the loop data evenly among the workers,
        # pack them up as addendums to the original package,
        # and send them off for processing.
        #----------------------------------------------------
        job_groups = equal_balance(loop_data,len(self.workers))
        addendums = []
        for grp in job_groups:
            addendums.append({loop_var:grp})
        results = self._send_recv(package,addendums)
        # Nothing done here to figure out the output format.
        # It is always returned as a tuple
        # Probably will be handier to have it as a 
        # Numeric array sometimes.
        out = ()
        for result_group in results:
            out = out + result_group
        return out

    def __getitem__(self, key):
        # currently allowing tuples also!
        #assert(type(key) is type(''))
        package = self.workers[0].get_pack(key)    
        return self._send_recv(package)

    def __setitem__(self, key, item):
        assert(type(key) is type(''))
        package = self.workers[0].update_pack({key:item})    
        return self._send_recv(package)

    def __delitem__(self, key):
        # currently allowing tuples also!
        # assert(type(key) is type(''))
        package = self.workers[0].del_pack(key)    
        return self._send_recv(package)
        
    def update(self, dict):
        package = self.workers[0].update_pack(dict)    
        return self._send_recv(package)

def equal_balance(jobs,Nworkers):
    """ Split jobs into Nworkers equal groups.    
        When an equal split is not possible,
        the larger groups occur at the front
        of the list.
    """
    
    #no jobs to do - return empty group list.
    if not len(jobs): return ()
    Ntotal_jobs = len(jobs)
    
    # find the number of jobs each wroker must do
    # for everyone to have equal work loads
    group_size = Ntotal_jobs / Nworkers
    
    # if there are jobs left over, some of the workers
    # will need to do 1 extra job.
    if Ntotal_jobs % Nworkers:
        group_size = group_size + 1

    # after some algebra, we can solve for the
    # number, a, of workers that need to do extra work
    a = Ntotal_jobs + Nworkers - Nworkers*group_size                       
    if a*group_size < Ntotal_jobs:
        b = Nworkers - a
    else:
        b = 0        
    
    # a workers do an extra job, b workers do standard
    # number of jobs.        
    group_sizes = a*[group_size] + b*[group_size-1]
    
    # now split the jobs up into groups for each of
    # the workers.
    last = 0
    job_groups = []
    for size in group_sizes:
        next = last+size
        job_groups.append(jobs[last:next])
        last = next
#    sum = 0        
#    for grp in job_groups:
#        sum = sum + len(grp)        
#    assert(sum,Ntotal_jobs)           
    return tuple(job_groups)

def t_func(i):
    return i    


if __name__ == '__main__':
    borg = cluster(server_list)    
    borg.start()
