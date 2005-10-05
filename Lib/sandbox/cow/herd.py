"""
    Main functions on a cluster:

    cluster.start()      -- start the cluster
    cluster['a'] = 1     -- assignment "a=1" on all machines
    all_a = cluster['a'] -- retreive a list of 'a' variable from all machines
    del cluster['a']     -- delete 'a' from all machines
    
    cluster.exec_code(code,input_vars,returns,global_vars)
        Execute a code fragment on a remote machine
            code        -- the code fragment to execute
            input_vars  -- dictionary of variables to use
            returns     -- list of variables to return from call
            global_vars -- variables to use from global namespace

    cluster.apply(function,args,keyword_args)
        Execute a function on a remote machine
            function     -- function to execute on remote machine
            args         -- argument list for function
            keyword_args -- dictionary of keyword arguments

    cluster.loop_apply(function,loop_var,args,keyword_args)
        Execute a function on a remote machine
            function     -- function to execute on remote machine
            loop_var     -- the variable to loop over
            args         -- argument list for function
            keyword_args -- dictionary of keyword arguments
            
    cluster.loop_code(code,loop_var,input_vars,returns,global_vars)
        Execute a code fragment on a remote machine
            code        -- the code fragment to execute
            loop_var    -- the variable to loop over
            input_vars  -- dictionary of variables to use
            returns     -- list of variables to return from call
            global_vars -- variables to use from global namespace
            
    >>> import herd
    >>> cluster = herd.cluster()
    >>> cluster.start()
    >>> cluster['a'] = 1    #set a=1 on all machines
    >>> print cluster['a']  #retreive a from all machines
    (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    >>> cluster.exec_code('import mom')
"""
Mrange = range(1,17)    
import cow # cow stands for "Cluster Of Workstations"
server_list = []
for i in Mrange:
    server_list.append(('cow%d.ee.duke.edu' % i,10000))
#for i in Mrange:
#    server_list.append(('cow%d.ee.duke.edu' % i,10001))

cluster = cow.machine_cluster(server_list)


