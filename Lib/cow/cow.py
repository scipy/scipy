""" Cow --- Cluster Of Workstations.
    Description

        cow is a package that includes tools for interacting with a cluster
        of workstations as a single computational engine.  It is well suited
        for "embarassingly parallel" problems that fit comfortably into a 
        master-slave paradigm.  For such problems, it is often possible to
        convert them from parallel to serial execution by changing a few lines
        of code.  Parallel Genetic Algorithms or Monte Carlo simulations are 
        two classes of problems that can benefit significantly from using
        cow.  If your problem requires closely coupled parallelism, take a 
        look at pyMPI.

        cow also includes a reasonably complete set of tools administering a 
        cluster of machines from Python such as process control, gathering
        system information and executing shell commands on remote processes
        simultaneously.
            
        cow should graze happily in heterogeneous cluster environments, though
        it hasn't tested extensively there.
        
        See scipy.cow.machine_cluster for more information and examples.
"""

import sync_cluster, socket
import time
import os # for getuid()

ClusterError = 'ClusterError'
TimeoutError = 'TimeoutError'
      
class machine_cluster:
    """ Treats a cluster of workstations as a single computational engine.
    
        Description
        
            machine_cluster simplifies interacting with a cluster of 
            workstations.  machine_cluster starts a group of slave
            interpreters.  Setting and retreiving global variables on these
            machines is accomplished through its dictionary interface.  
            You can also call functions or execute code blocks on the remote
            machines.  A couple of routines provide automatic parallelization
            of looping constructs.  It also provides routines for cluster
            wide process control similar to ps and cluster wide machine
            information on load, CPU, and memory information.
            
            
        Caveats
        
            Some functionality is only available on Linux. Unix is better
            supported in general than MSWindows.  See user manual for details.
        
            Cow assumes you have ssh access to all slave machines.  If not,
            the start() method will fail and you'll need to start all slaves
            by hand using::
            
                python sync_cluster.py server <port>
    
        Example Usage::
        
            # start two slave interpreters on the local machine at the
            # specified ports
            >>> slave_list = [ ('127.0.0.1',10000), ('127.0.0.1',10001)]
            >>> cluster = scipy.cow.machine_cluster(slave_list)
            >>> cluster.start()
            # set and retreive a global variable on the slave interpreters
            >>> cluster['a'] = 1 
            >>> cluster['a']
            (1, 1)
            >>> import string
            # process a list of strings in parallel converting them to
            # upper case (illustrative example only -- dismal performance)
            >>> string_list = ['aaa','bbb','ccc','ddd']
            >>> cluster.loop_apply(string.upper,0,(string_list,))
            ['AAA','BBB','CCC','DDD']
            >>> cluster.info()
            MACHINE   CPU        GHZ   MB TOTAL  MB FREE   LOAD
            bull      2xP3       0.5     960.0     930.0   0.00
            bull      2xP3       0.5     960.0     930.0   0.00
    """
    def __init__(self,server_list):        
        """ machine_cluster(slave_list) --> cluster_object
        
            Description
            
                Create a cluster object from a list of slave machines.
                slave_list is a list of 2-tuples of the form (address, port)
                that specify the network address and port where the slave
                interpreters will live and listen.  The address is always
                a string and can be either the machine name or its IP address.
                The port should be an unused port on the slave machine.
                Always use a port number higher than 1024.
            
            Example::
                
                # example 1 using IP addresses
                >>> slave_list = [ ('127.0.0.1',10000), ('127.0.0.1',10001)]
                >>> cluster = scipy.cow.machine_cluster(slave_list)                
                # example 2 using machine names                
                >>> slave_list = [ ('node0',11500), ('node1',11500)]
                >>> cluster = scipy.cow.machine_cluster(slave_list)
        """
        self.workers=[]
        self.worker_by_name={}
        worker_id = 1
        for host,port in server_list:
            # Add the uid here can help with port conflicts, but only works
            # on Unix clusters.  We really need to work out a daemon service
            # model that makes the port mess transparent.
            port = port #+ os.getuid()
            new_worker = sync_cluster.standard_sync_client(host,port,worker_id)
            self.workers.append(new_worker)
            self.worker_by_name[host] = new_worker
            worker_id = worker_id + 1
                                   
    def start(self,force_restart=0,timeout=60):
        """ start(force_restart=0, timeout=60) -- Start remote slave processes.
        
            Description
            
                Start the remote slave interpreters in the cluster.  The timeout
                value is specified in seconds and defaults to 60.  The timeout
                starts counting down only after ssh/rsh has tried to start
                all the remote processes.  This means the actual time in the
                function could be much longer than 60 seconds - depending
                on how long rsh/ssh takes.
                
                Its possible the 60 second time out will be to short for large
                clusters - but I hope not!
                
            Caveats                
                
                start() is not supported on MSWindows because of the lack of
                standard/robust support for remote startup and background 
                processing in the CMD shell.
        """    
        if not force_restart and self.is_running():
            return
        # start the worker processes.
        for worker in self.workers:
            worker.start_server()

        if not self.is_running():
            print '                      Starting Servers'
            print ' |----|----|----15---|----|----30---|----|----45---' \
                   '|----|----60'
            print '0.',
            stop_watch = timer()
            stop_watch.start()
            minute = 0
        import sys    
        while not self.is_running():
            if stop_watch.current_lap() > 1:
                sys.stdout.write('.')
                stop_watch.mark_lap()               
            elapsed = stop_watch.elapsed()
            if (elapsed - minute * 60) > 60:
                minute = minute + 1
                print
                print minute,
            if elapsed > timeout:
                raise TimeoutError  
        print 'servers running!'
            
    def stop(self):
        """ stop() -- Tell all remote slaves to terminate.
        
            Description
            
                stop calls sys.exit(0) on all the slave processes so that 
                they will terminate gracefully.  Note that if, for some
                reason, you are unable to connect to a remote processes due
                to some socket error, you'll have to kill the slave process
                by hand.
        """
        for worker in self.workers:
            import sys; sys.stdout.flush()
            try: worker.exec_code('import sys;sys.exit(0)')
            except:
                #should really do something here to
                # trap non-SystemExit errors.  
                pass
    def restart(self):
        """ restart() -- terminate all remote slaves and restart them.
        
            restart is useful when you would like to reset all slave
            interpreters to a known state.
        """
        self.stop()
        self.start()

    def is_running(self,timeout=0):
        """ is_running(timeout=0) --> 0 or 1
        
            check all the slave processes in the cluster are up and running.
            if timeout is specified, is_running will continually check if the 
            cluster is_running until it either gets a positive result or gives
            up and returns 0 after the specified number of seconds.
        """
        
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
        """ _send(package, addendum=None) -- send a package to all slaves.
        
            Description
            
                _send takes a package packed up by a packer object (see
                sync_cluster) and sends it to each of the slave processes.
                addendum is either None or a list with the same length as
                there are slave processes.  Each entry is a small package
                of additional information that is to be sent to a specific
                slave process.  It contains data that is only needed by that
                process.
                
            Implementation Notes
            
                The send is done synchronously to each worker in turn.  The
                entire package is sent to slave0 before moving on and sending
                the message to slave1.
                
                If a socket error occurs while trying to send data to a
                given slave, the offending worker is pushed into the 
                self.has_send_error list.  Also, self.send_exc is a dictionary
                that stores the (err_type,err_msg) as the key and the offending
                worker as the value.  This information is used in recv to 
                skip receiving from slaves who failed on send and also for 
                error reporting.
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
                err_type, err_msg = str,sys.exc_info()[:2]
                self.had_send_error.append(self.workers[i])
                key = (err_type,err_msg)
                try:
                    self.send_exc[key].append(self.workers[i].id)
                except: 
                    self.send_exc[key] = [self.workers[i].id]
            # else - handle other errors?               
        self.Nsent = N
                    
    def _recv(self):
        """ _recv() -- retreive results from slave processes.
        
            Description
            
                Retreive results from all slave processes that were 
                successfully sent a package.  If an error occurs while
                receiving from one of the slaves, the error is noted and
                the results from the other slaves are retreived.  A tuple
                the results from all workers is returned by _recv.  An
                entry of None is placed in the tuple for any worker that
                had an error.      
                
                The recv is done synchronously, waiting for slave0 to return
                its results before moving on to slave1 to recv its results.          
        """

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
                    sys.stdout.flush()
                    results.append(worker.recv())
                except sync_cluster.RemoteError:    
                    import sys
                    err = sys.exc_info()[1]
                    # Force the err msg (err[1]) to be a string.
                    # This dimishes info content, but makes sure
                    # that the sames errors are hashed correctly
                    # in the dictionary. (does it?)   
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
        #print                
        return tuple(results)                
    
    def _send_recv(self,package,addendum=None):
        """ _send_recv(package,addendum=None) --> results
        
            Description
            
                send a message to each worker in turn and then immediately
                began listening for the results.  All sends are done before
                listening for results from any of the slave processes.  See
                _send and _recv for more information.
                
                If an error occurs during either the send or recv phases,
                the handle_error() method is called.  If know errors are found,
                a tuple containing the results from each slave is returned.
                
                If an error does occur and an exception is raised, it is still
                possible to retreive the set of results that executed correctly
                from the last_results attribute.
        """        
        self._send(package,addendum)
        self.last_results = self._recv()
        if(len(self.send_exc) or len(self.recv_exc)):
            self.handle_error()
        return self.last_results    

    def handle_error(self):
        """ handle_error() -- make sense of send and recv errors
        
            Description
            
                Error handling attempts to examine the errors that occuer
                during remote execution and report them in the least verbose
                manner.  If the same error occurs on all slaves, it tries to
                only report it once.  Otherwise it reports all the errors
                that occur on slaves and prints the slaves traceback.
                
                Currently error handling is pretty simplistic.  It'd be nice
                if socket errors were viewed as severe and the slave either
                restarted or marked as dead and its work distributed among
                the other workers.
                
        """
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
        
    ##############################################################
    # slave processor info
    ##############################################################
    def load(self):
        """ load() -- print human readable load information for slave hosts 
            
            Description
            
                The load value printed is the 1 minute load average that
                is commonly printed by uptime on Unix machines.
                
                load depends on the implementation of scipy_proc on each
                slave's host OS. It will not work for Windows slave processes.
                However, if you are using a Windows master to control a Linux
                cluster of slaves, it should work fine.                

            Example::

                >>> slave_list = [('n0',10000), ('n1', 10001)]
                >>> cluster = scipy.cow.machine_cluster(slave_list)
                >>> cluster.start()
                >>> cluster.load()
                    n0: 0.00, n1: 0.00
        """        
        import string
        import scipy_distutils.proc as scipy_proc
        results = self.load_list()
        for i in range(len(self.workers)):            
            name = string.split(self.workers[i].host,'.')[0]
            res = results[i]
            s = "%6s: %1.2f," % (name[-6:], res['load_1'])
            print s,
            if not ((i+1) % 5):
                print
    
    def info(self):
        """ info() -- print human readable info about the slave hosts 
            
            Description
            
                Print out the each slave interpreters host name, number
                and type of processors, memory usage, and current load 
                information in human readable form.
                
                info depends on the implementation of scipy_proc on each
                slave's host OS. It will not work for Windows slave processes.
                However, if you are using a Windows master to control a Linux
                cluster of slaves, it should work fine.                

            Example::

                >>> slave_list = [('n0',10000), ('n1', 10001)]
                >>> cluster = scipy.cow.machine_cluster(slave_list)
                >>> cluster.start()
                >>> cluster.info()
                MACHINE   CPU        GHZ   MB TOTAL  MB FREE   LOAD
                n0        2xP3       0.4     192.0      13.0   0.00
                n1        2xP3       0.4     192.0      22.0   0.00
        """                
        import string
        results = self.info_list()
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

    def load_list(self):        
        """ load_list() -- Return a list of slave load information dictionaries
            
            Description
            
                Retreive a dictionary with information about the load on each
                host processor.  The dictionaries have three keys, load1, 
                load5, and load15 indicating the 1, 5, and 15 minute load 
                averages for the processor.  These could be useful for (as 
                yet unimplemented) load balancing schemes.

                load_list depends on the implementation of scipy_proc on each
                slave's host OS. It will not work for Windows slave processes.
                However, if you are using a Windows master to control a Linux
                cluster of slaves, it should work fine.                

            Example::

                >>> slave_list = [('n0',10000), ('n1', 10001)]
                >>> cluster = scipy.cow.machine_cluster(slave_list)
                >>> cluster.start()
                >>> cluster.load_info()
                ({'load_5': 0.0, 'load_1': 0.0, 'load_15': 0.0},
                 {'load_5': 0.0, 'load_1': 0.0, 'load_15': 0.0})
        """        
        import scipy_distutils.proc as scipy_proc
        res = self.apply(scipy_proc.load_avg,())
        return res

    def info_list(self):
        """ info() -- print human readable info about the slave hosts 
            
            Description
            
                Print out the each slave interpreters host name, number
                and type of processors, memory usage, and current load 
                information in human readable form.
                
                info depends on the implementation of scipy_proc on each
                slave's host OS. It will not work for Windows slave processes.
                However, if you are using a Windows master to control a Linux
                cluster of slaves, it should work fine.                

            Example::

                >>> slave_list = [('n0',10000), ('n1', 10001)]
                >>> cluster = scipy.cow.machine_cluster(slave_list)
                >>> cluster.start()
                >>> cluster.info()
                MACHINE   CPU        GHZ   MB TOTAL  MB FREE   LOAD
                n0        2xP3       0.4     192.0      13.0   0.00
                n1        2xP3       0.4     192.0      22.0   0.00
    
        """                
        import scipy_distutils.proc as scipy_proc
        res = self.apply(scipy_proc.machine_info,())
        return res

    ##############################################################
    # slave process information and control
    ##############################################################
    def ps(self,sort_by='cpu',**filters):
        """ ps(sort_by='cpu',**filters) -- list processes on slave machines.
        
            Description
                
                List all the processes on all remote slave machines.  This
                is like a cluster-wide Unix ps command and is output in a
                similar human readable form.  The sort_by argument allows
                you to sore the process list by various fields including,
                pid, cpu, user, machine, memory, state and command.  keyword
                arguments are used as filters to limit the number of processes
                displayed.  For example, the keyword, user='ej' will only list
                processes for user ej and cpu='>10' will only list processes 
                using more th 50% of the cpu cycles.
                
                ps depends on the implementation of scipy_proc on each
                slave's host OS. It will not work for Windows slave processes.
                However, if you are using a Windows master to control a Linux
                cluster of slaves, it should work fine.                

            Example::

                >>> slave_list = [('n0',10000), ('n1', 10001)]
                >>> cluster = scipy.cow.machine_cluster(slave_list)
                >>> cluster.start()
                >>> cluster.ps(user='ej')
                MACHINE USER       PID  %CPU  %MEM TOTAL MB   ...
                n0     ej        22915  99.9   2.1    4.027   ...
                n0     ej        22916  99.9   2.1    4.055   ...
                n1     ej        22915  99.9   2.1    4.027   ...
                n1     ej        22916  99.9   2.1    4.055   ...
                ...
                
        """
        psl = self.ps_list(sort_by,**filters)
        if len(psl):
            print psl[0].labels_with_name()
        for i in psl: print i.str_with_name()
    
    def ps_list(self,sort_by='cpu',**filters):
        """ ps_list(self,sort_by='cpu',**filters) -- get cluster processes
        
            Description
            
                Return a list containing one scipy_proc.process objects for 
                each process running on the cluster host machines.  process 
                objects contain a ton of information about cpu, memory, etc.
                used by the process.
                
                See ps for more information.
            
            Example:
                        
                >>> slave_list = [('n0',10000), ('n1', 10001)]
                >>> cluster = scipy.cow.machine_cluster(slave_list)
                >>> cluster.start()
                >>> p = cluster.ps_list()
                >>> for i in p: print p.pid
                ...
                22890 22889 22889 22890 1071 1071 ...
        """
        import operator
        import scipy_distutils.proc as scipy_proc
        res = self.apply(scipy_proc.ps_list,())
        psl = reduce(operator.add,res)
        psl = scipy_proc.ps_sort(psl,sort_by,**filters)        
        return psl
 
    def nice(self,increment=10):
        """ nice(increment=10) --> success_list
        
            increment all slave interpreter's nice value by increment.            
            hmmm. this doesn't seem to work. see os.nice()
        """
        res = self.apply(os.nice,(increment,))
        return res

    def renice(self,process_list,level):
        """ renice(process_list,level) -- set nice value multiple processes 

            Description
            
                Change the nice level of multiple remote processes. 
                process_list is a list of scipy_proc.process objects. 
                level is the new nice value for the listed processes.
             
            Caveats

                Once niced down, a process cannot be reniced back up.
                This is a Linux issue.
        """    
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
        """ kill(self,process_list,signal = 'TERM') -- Signal process list.

            Description
            
                Send a signal to all of the scipy_proc.process objects in
                the process_list.  This is usually used to kill the processes.
                The signal may be given as a signal name or number.
        """    
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
        
    def system(self,cmd):
        """ system(cmd) -- execute cmd on all remote machines

            A list of all the remote responses is returned.  Unlike
            os.system which returns the exit value of the cmd string,
            this function returns the text output by the command.
        """            
        code = 'import os;f=os.popen("%s");res = f.read(-1);f.close();' % cmd
        return self.exec_code(code,returns=['res'])

    def reload(self,module):
        """ reload(module) -- reload module on all remote interpreters
        
            module can either be the name of a module or the actual
            module object.
        """        
        try: 
            code = 'import %s; reload(%s)' % ((module.__name__,)*2)
        except AttributeError:
            code = 'import %s; reload(%s)' % ((module,)*2)
        self.workers.exec_code(code)    
        
    ##############################################################
    # remote code and function execution
    ##############################################################
    
    # mirror all of sync_client functions
    # They assumes all clients have the same packing procedures.
    def exec_code(self,code,inputs=None,returns=None):
        """ exec_code(code,inputs=None,returns=None)
        
            Similar to Python's exec statement. Execute the same code fragment
            on all remote interpreter. inputs is a dictionary of variable 
            values to use when executing the code. returns is a list of
            variable names that should be returned after executing the code.  
            If one value is specified, the value for that variable is returned.
            If multiple values are specified, a tuple is returned.
            
            exec_code returns a list of the values requested variables,
            one entry for each slave.
        """
        #use the first worker to package up the cmd.
        package = self.workers[0].exec_code_pack(code,inputs,returns)
        return self._send_recv(package)

    def apply(self,function,args=(),keywords=None):        
        """ apply(function,args=(),keywords=None)
        
            Similar to Python's builtin apply method.  Execute the given
            function with the argument list, args, and keyword arguments,
            keywords, on each of the slave processes.
            
            apply returns a list of the results from calling function,
            one result for each slave.
        """        
        package = self.workers[0].apply_pack(function,args,keywords)
        return self._send_recv(package)

    def loop_apply(self,function,loop_var,args=(),keywords=None):
        """ loop_apply(function,loop_var, args=(),keywords=None)
        
            Description

                Call function with the given args and keywords.  One of the
                arguments or keywords is actually a sequence of arguments.  
                This sequence is looped over, calling function once for each
                value in the sequence. loop_var indicates which variable to
                loop over.  If an integer, loop_var indexes the args list.
                If a string, it specifies a keyword variable.  The loop sequence
                is divided as evenly as possible between the worker nodes and
                executed in parallel.

            Example::
                
                >>> slave_list = [('n0',10000), ('n1', 10001)]
                >>> cluster = scipy.cow.machine_cluster(slave_list)
                >>> cluster.start()
                >>> import string
                >>> string_list = ['aaa','bbb','ccc','ddd']
                >>> cluster.loop_apply(string.upper,0,(string_list,))
                ['AAA','BBB','CCC','DDD']
            
        """                
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

    def loop_code(self,code,loop_var,inputs=None,returns=None):
        """ loop_code(code,loop_var,inputs=None,returns=None)
        
            Description
        
                Similar to exec_code and loop_apply.  Here loop_var indicates 
                the variable name in the inputs dictionary that is looped over.
            
            Example::
                
                >>> slave_list = [('n0',10000), ('n1', 10001)]
                >>> cluster = scipy.cow.machine_cluster(slave_list)
                >>> cluster.start()
                >>> import string
                >>> a = [1, 2, 3, 4]
                >>> cluster.loop_code("b=a*2",'a',{'a':a},('b',))
                (2, 4, 6, 8)
        """
        the_inputs = {}
        the_inputs.update(inputs)    
        loop_data = the_inputs[loop_var]
        the_inputs[loop_var] = None #make it small for packing
        package = self.workers[0].loop_code_pack(code,loop_var,
                                                 the_inputs,returns)
        return self.loop_send_recv(package,loop_data,loop_var)
    
    # array specific routines
    def row_split(self,name,sequence):
        """experimental"""
        import scipy
        q=scipy.split(sequence,len(self.workers))
        self.loop_code(name+'=_q_','_q_',inputs={'_q_':q},returns=(),
                        global_vars=(name,))
    def row_gather(self,name):
        """experimental"""
        from Numeric import concatenate
        return concatenate(self[name])
        
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

    ##############################################################
    # dictionary interface
    ##############################################################

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

import operator
class timer:

    def __init__(self):
        self.reset()
    def reset(self):
        self.total_t = 0
        self.lap_start = 0
        self.running = 0
        self.lap_list = [0]        
    def start(self):
        if not self.running:
            self.running = 1
            self.start_t = time.time()                
        else:
            print 'already running: use reset() to start the timer back at 0.'
                
    def stop(self):
        self.running = 0
        elapsed_t = time.time() - self.start_t
        self.total_t = self.total_t + elapsed_t
        self.lap_list[-1] = self.lap_list[-1] + elapsed_t
        return self.total_t
    def elapsed(self):
        if self.running:
            return self.total_t + (time.time() - self.start_t)
        else:
            return self.total_t
    def mark_lap(self):
        self.lap_start = self.elapsed()
        self.lap_list.append(self.lap_start-self.lap_list[-1])
    def get_laps(self):
        return self.lap_list[1:] +[self.elapsed()-self.lap_list[-1]]               
    def current_lap(self):
            return self.elapsed()-self.lap_start
        
                

if __name__ == '__main__':
    borg = cluster(server_list)    
    borg.start()
