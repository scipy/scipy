#!/usr/bin/env python
import cPickle, cStringIO
import time
pickle = cPickle
#import zlib
import socket
import sys

import sync_cluster #yes I'm importing the current module

bufsize = 1<<12 #4K buffer, for Linux
shell = "ssh -X"

RemoteError = 'RemoteError'
RemoteCrashError = 'RemoteError'
PackError = 'PackError'
UnpackError = 'UnpackError'
NotImplemented = 'NotImplemented'


class pickle_packer:
    """* Pickle and unpickle an object for transfer over a socket.

        Description:
        This class takes arbitrary python objects and turns them
        into a character string suitable for transfer over a socket.
        It is currently a thin wrapper around cPickle routines, and
        thus can only pack "picklable" objects (see the python manual).
        This is not much of a limitation.

        Why do this instead of using pickle directly:
        The class is abstracted in the hopes that future packers
        might pack items much more quickly and "in place" instead of
        creating intermediate string buffers etc.  I now think a
        slightly different architecture is necessary for this to
        be possible.  The packer will need to be something like
        a "producer" in medusa.
    *"""
    def pack(self,item):
        """* Pack item and return resulting string. *"""
        try:
            packed  = pickle.dumps(item,1)
        except pickle.PicklingError, why:
            raise PackError, why
        return packed
    def unpack(self, x):
        """* Unpack string x and return the resulting object."""
        # If we get an empty string back, more than likely
        # the remote worker crashed. - They are always supposed
        # to send back something.
        #if x == '':
        #    raise RemoteCrashError
        try:
            res = pickle.loads(x)
        except pickle.UnpicklingError, why:
            raise UnpackError, why
        return res

    def fpack(self,item,f):
        """* Pack item into the file object f. Nothing is returned.
            f can also be a cStringIO object. fpack can be
            called multiple times with the same file object
            to pack several items into f.
        *"""
        try:
            packed  = pickle.dumps(item,1)
        except pickle.PicklingError, why:
            raise PackError, why

    def funpack(self, f):
        """* Unpack an object from thr file object f.  The
            resulting object is returned. f can also be
            a cStringIO object.  if multiple items are
            packed within f, funpack can be called multiple
            times to retreive the objects.
        *"""
        try:
            res = pickle.load(f)
        except pickle.UnpicklingError, why:
            raise UnpackError, why
        return res

class blocking_channel:
    """* Manages sending and receiving data through a blocking socket.

        Description:
        The channel is meant for a single exchange of data between
        two programs.  Generally, data is sent down the write channel
        (wfile) using write().  The write channel is then closed
        with close_write() signally to the other end that all
        data has been sent.  The program on the other end then reads
        the data, processes it, and sends a response back down the
        read channel (rfile).  You can read the the results using recv().

        Many programs, like chat, keep a socket open for many exchanges.
        This is possible when sending text data back and forth because
        you can indicate the end of a sent message with a delimiter
        character.  When sending binary data, however, it is not possible
        to break up messages with a delimiter because the delimiter
        might occur randomly within the binary data indicating a false
        message break.  Closing the write channel is an unambiguous
        end of message indicator.

        This channel is used by the standard_client to send snippets
        of python code along with pickled objects to a remote machine.
        The remote machine then executes the code and returns the results
        to the original machine.

        Because the underlying sockets are blocking, this method may not
        be the most efficient way to communicate with many different
        clients.  Non-blocking sockets, as are used in medusa, may
        be more efficient, but the programming model is slightly more
        confusing.  See ***.
    *"""
    def __init__(self,host,port,log = None):
        """* Open a channel (socket) to a remote machine.  The
             socket attempts to connect to the host on the
             indicated port.
        *"""
        self.sock = socket.socket( socket.AF_INET, socket.SOCK_STREAM)
        #print 'connecting: ', host, port
        self.sock.connect((host,port))
        self.host,self.port = host,port
        #print 'creating files'
        self.rfile = self.sock.makefile('rb', bufsize)
        self.wfile = self.sock.makefile('wb', bufsize)
        import sys
        self.log = log
        #self.log = sys.stdout
    def write(self,x):
        """* Send data to the remote machine.  This function
            can be called multiple times.  The net affect, though
            is that all data sent is concatenated together to
            form a single message read on the other end.

            The remote machine will not read any of the data sent
            until close_write() is called to signal that all data
            is sent.  Calls to write() after close_write() has
            been called will result in an error.
        *"""
        self.log_msg( ("writing %d bytes" % len(x)) )
        self.wfile.write(x)
    def close_write(self):
        """* Close the write channel.  This signals that all data
             has been sent to the remote end of the connection.
             Further calls to write() will result in an error.
        *"""
        self.log_msg("closing write channel")
        self.wfile.flush()
        self.wfile.close()
        self.sock.shutdown(1) #tells other side transfer is over
    def read(self):
        """* Return a string of data returned by the remote connection.
            read() should not be called until after close_write(). It
            is a blocking call and will wait indefinitely if the
            remote connection does not send back data and close the
            connection.

            read() stops reading data when the rfile is closed
            by remote connection signally the end of data.  It then
            closes the socket.  read() can only be called once.
            Further calls will result in an error.
        *"""
        x = self.rfile.read()
        self.log_msg( ("read %d bytes" % len(x)) )
        self.rfile.close()
        self.sock.close()
        return x
    def write_read(self,data):
        """* Send data to the remote connection, wait for it to
            process the data, and return a response.  The
            response is returned as a string.  After calling this
            function the socket is closed and further calls
            to write(), read() or wrtie_read() will result in
            errors.
        *"""
        sent = self.write(data)
        self.close_write()
        return self.read()
    def log_msg(self,msg):
        if self.log:
            pre_msg = 'blocking_channel on %s, %d:' % (self.host,self.port)
            self.log.write(pre_msg+msg+'\n')
            self.log.flush()

class standard_sync_client:
    """* Client side of a remote Python interpreter.

        This class allows you to connect to remote machines
        and send them arbitrary python code and objects.
        The remote machine executes the code and sends it
        back to you.  It is intended to work with
        standard_sync_handler on the server side.  Together,
        they make a powerful building block for "embarassingly
        parallel" distributed computing.  Note, however,
        that this class used blocking sockets (hence the
        "sync" in the name).  A non-blocking client
        may lead to more efficient use of available bandwidth.

        Note that, when working with a cluster of computers,
        using wscluster is much easier than trying to manage
        the machine client on your own.  It handles splitting
        up the processing across a cluster of workstations much
        more transparently.

        There are several ways to used standard_sync_client:

            1) As a dictionary for setting and getting
               global variables on a remote interpreter.
            2) To execute a function in a similar fashion
               to calling apply()
            3) To execute a code fragment similar to
               calling exec
            4) loop() calls the same function with multiple
               different inputs. (SIMD)

        This functionality is accessible in an two different
        formats, "easy to use" or "most efficient".  ** something
        about load balancing statistics here **

        The client and server machines are assumed to have
        identical copies of all module files.  Also, any
        object, function or class sent to remote machines
        must be defined in a module.  Objects are "pickled"
        before being sent and pickling objects
        from __main__ causes problems when unpickling
        them on the remote end.

        apply() example:
        Suppose you have a function process_image() that
        processes an image and returns some result.

            # in module radar.py
            # ** make this a real function
            # that does something time consuming **
            def process_image(image):
                H = fft2(image)
                # some other expensive operations
                return result

       Here's one way to execute this on remote machines.
       This assumes a server is up and running on the machines.

            # 1. Create an 256x256 image of random pixels
            #    between 0 and 1.
            import scipy.stats 
            image1 = scipy.stats.uniform(0.,1.,(256,256))
            image2 = scipy.stats.uniform(0.,1.,(256,256))
            # 2. Create connections to remote machines
            #    named cow3 and cow4 on port 10000
            remote1 = standard_sync_client('cow3',10000)
            remote2 = standard_sync_client('cow4',10000)
            # 3. Call apply() on the remote machine.
            result1 = remote1.apply(radar.process_image,(image1,))
            result2 = remote2.apply(radar.process_image,(image2,))

        The remote apply() has identical calling semantics
        to the standard one.  See the primer (below) on the standard
        Pyton apply for more insight into how this function works.

        if you time the previous calls, you'll note that they don't
        occur twice as fast simply because they are processed
        on two different machines.  Calls to apply()
        are blocking.  Image1 is processed completely before
        image2 is processed here.  More than likely you'd rather
        the work happen in parallel.  This requires you to
        separate the read and write processes.

            # 3. Call apply() on the remote machine.
            #   a) call apply pack to prepare the messages.
            #   b) send to each machine.
            #   c) read from each machine.
            task1 = remote1.apply_pack(radar.process_image,(image1,))
            task2 = remote2.apply_pack(radar.process_image,(image2,))
            remote1.send(task1)
            remote2.send(task1)
            result1 = remote1.recv(task1)
            result2 = remote1.recv(task2)

        The send and recv are still blocking, but if the processing
        time for the images are much larger than the packing and
        transmit time, you'll get your images processed in about
        half the time.

        exec() example:
        Executing code fragments proceeds in a similar fashion.
        It is modeled after the exec statement, but takes
        more parameters as a consequence of being remote.
        Additionly, exec is a function here instead of a statement.
        # **more here**

        apply() primer:
        The apply() function in Python is used to "apply"
        a set of arguments to a function takes three arguements.
        The first is a function, the second is a tuple
        of arguments, and the third is a dictionary of
        keyword arguements. for example if the function is:

            def silly_func(a,b,c,d=2):
                return a+b+c+d

       then apply can be called by passing all arguments
       in through the argument list

           >> r = apply(silly_func,(1,2,3,4))
           r = 10

        or some can be passed in using the argument list
        and others using keywords

           >> r = apply(silly_func,(1,),{'b':2,'c':3,'d':4})
           r = 10

        Note the "," after 1 here. The second
        argument must be tuple.  Python
        interprets (1) as the number 1 and (1,) as
        a tuple with 1 as its only entry.  Such
        silliness goes away if the function takes
        more than one arguement.

        See the standard Python tutorial for more
        info on calling functions with apply().
    *"""
    def __init__(self,host,port,my_id=-1,log=None):
        self.packer = pickle_packer()
        self.host,self.port = host,port
        self.id = my_id
        self.log = log
    def start_server(self):
        """* expands to something like:
        rsh -n 127.0.0.1 "python ~/wrk/sync_cluster.py server 10000 >&/dev/null </dev/null &"
        *"""
        import sync_cluster, os
        module_name = os.path.abspath(sync_cluster.__file__)
        #cmd = 'rsh -n %s "python %s server %d >&/dev/null </dev/null &"' % \
        #      (self.host, module_name, self.port)
        cmd = '%s %s "python %s server %d >&/tmp/crud%d </dev/null &"' % \
              (shell,self.host, module_name, self.port, self.port)
        self.log_msg(cmd)
        #print cmd
        os.system(cmd)

    def exec_code_pack(self,code,inputs=None,returns=None,global_vars=None):
        contents={'_command':sync_cluster.exec_code,'code':code,
                  'inputs':inputs,'returns':returns,
                  'global_vars':global_vars }
        return self.packer.pack(contents)

    def loop_code_pack(self,code,loop_var,inputs=None,returns=None,global_vars=None):
        contents={'_command':sync_cluster.loop_code,'code':code,
                  'loop_var':loop_var,'inputs':inputs,
                  'returns':returns,'global_vars':global_vars }
        return self.packer.pack(contents)
    def apply_pack(self,function,args,keywords=None):
        if not keywords: keywords = {}
        #more_keywords = args_to_keywords(function,args)
        #catch_keyword_conflicts(more_keywords,keywords)
        #keywords.update(more_keywords)
        contents={'_command':sync_cluster.apply_func,'function':function,
                  'args':args,'keywords':keywords}
        return self.packer.pack(contents)
    def loop_apply_pack(self,function,loop_var,args,keywords=None):
        # NOT USED -- see cow.machine_cluster.loop_apply.
        all_keywords = {}
        if keywords: all_keywords.update(keywords)
        more_keywords = args_to_keywords(function,args)
        catch_keyword_conflicts(more_keywords,all_keywords)
        all_keywords.update(more_keywords)
        contents={'_command':sync_cluster.loop_func,'function':function,
                  'keywords':all_keywords,'loop_var':loop_var}
        return self.packer.pack(contents)

    def apply(self,function,args,keywords=None):
        package = self.apply_pack(function,args,keywords)
        self.send(package)
        return self.recv()
        
    def exec_code(self,code,inputs=None,returns=None,global_vars=None):
        package = self.exec_code_pack(code,inputs,returns,global_vars)
        self.send(package)
        return self.recv()

    def loop_apply(self,function,loop_var,args,keywords=None):
        # NOT USED -- see cow.machine_cluster.loop_apply.
        package = self.loop_apply_pack(function,loop_var,args,keywords)
        self.send(package)
        return self.recv()

    def loop_code(self,code,loop_var,inputs=None,returns=None,global_vars=None):
        # NOT USED -- see cow.machine_cluster.loop_code.
        package = self.loop_code_pack(code,loop_var,inputs,returns,global_vars)
        self.send(package)
        return self.recv()
    def send(self,package,addendum=None):
        """* addendum is a dictionary.  It is not packed...*"""
        # build a socket
        self.channel = blocking_channel(self.host,self.port)
        self.channel.write(package)
        if addendum:
            extra_package = self.packer.pack(addendum)
            self.channel.write(extra_package)
        self.channel.close_write()
    def recv(self):
        package = self.channel.read()
        contents = self.packer.unpack(package)
        #print 'self.id, ',self.id, contents
        self.catch_exception(contents)
        #print 'after exception'
        #try:    print  contents['_exec_time']
        #except: pass    
        return contents['result']
    def get_load_info(contents):
        #use this to read execution time info from the package...
        pass
    def log_msg(self,msg):
        if self.log:
            pre_msg = 'client %d:' % self.id
            self.log.write(pre_msg+msg+'\n')
            self.log.flush()

    #-------------------------------------------
    # Partial mimicking of a Dictionary.
    # Need some thought on how to deal with keys(),
    # items(), len(), clear(), and values()
    #
    # Only strings can be used as keys.
    #-------------------------------------------

    def get_pack(self,keys):
        contents={'_command':sync_cluster.get_keys,
                  'keys':keys}
        return self.packer.pack(contents)
    def set_pack(self,key,item):
        # relies on update_pack.
        contents={'_command':sync_cluster.update,
                  'global_dict':{key:item}}
        return self.packer.pack(contents)
    def update_pack(self,global_dict):
        contents={'_command':sync_cluster.update,
                  'global_dict':global_dict}
        return self.packer.pack(contents)
    def del_pack(self,keys):
        contents={'_command':sync_cluster.del_keys,
                  'keys':keys}
        return self.packer.pack(contents)

    def __getitem__(self, key):
        # currently allowing tuples also!
        #assert(type(key) is type(''))
        package = self.get_pack(key)
        self.send(package)
        return self.recv()

    def __setitem__(self, key, item):
        assert(type(key) is type(''))
        package = self.update_pack({key:item})
        self.send(package)
        self.recv() #must receive to catch exceptions

    def __delitem__(self, key):
        # currently allowing tuples also!
        # assert(type(key) is type(''))
        package = self.del_pack(key)
        self.send(package)
        self.recv() #must receive to catch exceptions

    def update(self, dict):
        package = self.update_pack(dict)
        self.send(package)
        self.recv() #must receive to catch exceptions

    def clear(self):raise NotImplemented
    def keys(self): raise NotImplemented
    def items(self): raise NotImplemented
    def values(self): raise NotImplemented
    def has_key(self, key): raise NotImplemented
    def get(self, key, failobj=None): raise NotImplemented

    #-------------------------------------------
    # End of dictionary.
    #-------------------------------------------

    def catch_exception(self,x):
        """* exception_traceback is really a string representation
            of the error.  The actual traceback can't be sent
            because of pickleing issues.
        *"""
        try:
            err = x['type'],x['msg'],x['exception_traceback']
            raise RemoteError, err
        except (TypeError, KeyError, AttributeError): pass

def args_to_keywords(function,args):
    arg_keywords = {}
    N = len(args)
    for i in range(N):
        arg_name = function.func_code.co_varnames[i]
        arg_keywords[arg_name] = args[i]
    return arg_keywords

def catch_keyword_conflicts(kw1, kw2):
    keys1 = kw1.keys()
    keys2 = kw2.keys()
    for key in keys1:
        if key in keys2:
            raise TypeError, ('keyword parameter "%s" redefined' % key)

def read_log():
    """ Read the results of the log file.  
    
        This assumes that stdout and the log file are the same.    
    """
    #global log_file
    #sys.stdout.close()
    #f = open(log_file,'r')
    #results = f.read(-1)
    #f.close()
    #sys.stdout = open(log_file,'a+')
    #return results
    pass
    
import SocketServer
class standard_sync_handler(SocketServer.StreamRequestHandler):
    verbose = 1
    packer = pickle_packer()
    def setup(self):
        """* This function is to fix a bug in the python library.
             The newest release of Python should have that bug fixed,
             and this function can be removed
        *"""
        self.connection = self.request
        print 'here:',self.connection
        self.rfile = self.connection.makefile('rb', bufsize)
        self.wfile = self.connection.makefile('wb', bufsize)
        
        #SocketServer.StreamRequestHandler.setup(self)
        #import tempfile,os
        #global log_file
        #dr = tempfile.gettempdir()
        #log_file = os.path.join(dr,'cow.log')
        #f = open(log_file,'w')
        #sys.stdout = f
        self.log = sys.stdout

    def handle(self):
        import time
        kill_when_finished = 0
        try:
            recv_msg = self.recv()
            #print 'received msg'
            send_msg = self.process(recv_msg)
            #print 'processed msg'
            packed_msg = self.packer.pack(send_msg)
            if self.verbose: print '\tserver: sending %d bytes' % len(packed_msg)
            #print 'send msg:', send_msg
        except SystemExit:
            import sys
            err = sys.exc_info()
            packed_msg = self.pack_exception(err)
            del err
            kill_when_finished = 1
            #print 'should exit'
        except:
            import sys
            err = sys.exc_info()
            packed_msg = self.pack_exception(err)
            #print 'packed exception:', err[0],err[1]
            del err

        self.wfile.write(packed_msg)
        #print 'sent msg'
        #print '\tserver: sent %d bytes' % len(packed_msg)
        if self.verbose: print '\tserver: sent %d bytes' % len(packed_msg)
        if kill_when_finished:
            print 'killing'
            self.finish() # make sure this socke this closed
            #signal main server to die
            import sys,os
            #ppid = os.getppid ()
            #print ppid
            force_kill()
            sys.exit(0)
        #print 'kill_when_finished', kill_when_finished

    def finish(self):
        """* This function shouldn't have to be defined, but, on occasion,
            when an Exception is thrown in the code processed by loop_code
            ( and maybe others ), the socket (request) doesn't die, even
            when all the files are closed.  As a result, the other end
            blocks forever.  We close the socket explicitly here, until
            we find the real problem.  The exceptions causing the problems
            are occuring in Fortran routines, but the execptions seem
            to be caught properly on this end and sent down the channel.
        *"""
        self.wfile.flush()
        self.wfile.close()
        self.rfile.close()
        self.request.close() # hmm.  Shouldn't need this, but occasionaly seems necessary.
        #print 'should be closed'

    def recv(self):
        if self.verbose: print '\tserver: connection established'
        #need to do something here to handle extra messages...
        packed_msg =self.rfile.read()
        if self.verbose: print '\tserver: received %d bytes' % len(packed_msg)

        packed_msg = cStringIO.StringIO(packed_msg)
        task = self.packer.funpack(packed_msg)
        current_location = packed_msg.tell()
        packed_msg.seek(0,2)
        end_location = packed_msg.tell()
        if(end_location != current_location):
            packed_msg.seek(current_location)
            addendum = self.packer.funpack(packed_msg)
        else:
            addendum = None
        task['addendum'] = addendum
        return task

    def process(self,task):
        # find out the command type (this is really a function)
        command = task['_command']
        # remove the command from the task.  The task now
        # represents
        del task['_command']
        t1 = time.time()
        result= apply(command,(),task)
        t2 = time.time()
        result_dict = {'result':result,'_exec_time':t2-t1}
        #add load info here
        return result_dict
    def pack_exception(self,err):
        import traceback,string
        msg = {'exception_traceback':
               '########### Traceback text from remote machine ###########\n' \
               + string.join(traceback.format_exception(err[0],err[1],err[2]),'') \
               + '################# End remote traceback ##################',
               'type': err[0],
               'msg': err[1]}
        return self.packer.pack(msg)

    def log_msg(self,msg):
        if self.log_name:
            pre_msg = 'client %d:' % self.id
            self.log.write(pre_msg+msg+'\n')

#there may be some issues with globals scope in these
def apply_func(function, args, keywords=None,addendum=None):
    if not keywords: keywords = {}
    if addendum: keywords.update(addendum)
    return apply(function,args ,keywords)

def exec_code(code,inputs,returns,global_vars,addendum=None):
    if addendum: inputs.update(addendum)
    if not returns: returns = ()
    if type(returns) == type(''):
        raise TypeError, 'returns must be a sequence object - not a string'
    exec_code = build_globals(global_vars)
    exec_code = exec_code + build_inputs(inputs)
    exec_code = exec_code + code
    globals()['_inputs'] = inputs
    exec exec_code in globals(), globals()
    #perhaps do something here to catch errors
    if len(returns) == 1:
        results = eval(returns[0])
    elif len(returns) > 1:
        results = []
        for i in returns:
            results.append(eval(i))
        results = tuple(results)
    else:
        results = None
    return results

def loop_func(function,loop_var,args,keywords,addendum=None):
    if not keywords: keywords = {}
    if addendum: 
        keywords.update(addendum)
    result = []
    _loop_data = keywords[loop_var]
    del keywords[loop_var] #not strictly necessary
    if type(loop_var) == type(''):                                
        for _i in _loop_data:
            keywords[loop_var] = _i
            result.append(apply(function,args ,keywords))
    elif type(loop_var) == type(1):
        args_list = list(args)
        for _i in _loop_data:
            args_list[loop_var] = _i
            args = tuple(args_list)
            result.append(apply(function,args,keywords))
    return tuple(result)

def loop_code(code,loop_var,inputs,returns,global_vars,addendum=None):
    if type(returns) == type(''):
        raise TypeError, 'returns must be a sequence object - not a string'
    if addendum: inputs.update(addendum)
    globals()['_loop_data'] = inputs[loop_var]
    globals()['_returns'] = returns 
    #added to set all inputs in the global namespace
    globals().update(inputs)
    del inputs[loop_var] #not strictly necessary
    exec_code = build_loop_code(code,loop_var,inputs,returns,global_vars)
    exec exec_code in globals(), globals()
    return _all_results

#------------------------------------------------
# dictionary like routines for global variables.
# All ignore addendum.
#------------------------------------------------
def set_key(key,item,addendum=None):
    globals()[key] = item

def update(global_dict,addendum=None):
    globals().update(global_dict)

def get_keys(keys,addendum=None):
    if type(keys) == type(''):
        return globals()[keys]
    else:
        results = []
        for i in keys:
            results.append(globals()[i])
        return tuple(results)

def del_keys(keys,addendum=None):
    if type(keys) == type(''):
        del globals()[keys]
    else:
        for i in keys:
            del globals()[i]

#------------------------------------------------
# utilities and template strings used for building
# loop_code and exec_code strings.
#------------------------------------------------

import string

def build_globals(global_vars):
    # if we have global variables specified, we
    # need to append there definition:
    #    global a, b, c, etc.
    # at the beginning of the code to be executed
    # This ensure that global variables set by the
    # code will actually be set globally instead
    # of locally
    if global_vars:
        return 'global %s;' % string.join(global_vars,',')
    else:
        return ''

def build_inputs(inputs):
    if inputs:
        asgn = map(lambda x: '%s = _inputs["%s"];' % (x,x), inputs.keys())
        return string.join(asgn,'')
    else:
        return ''

def indent(code):
    import string
    code = string.replace(code,'\n','\n    ')
    return code

loop_code_template = """
%(global_code)s
%(input_code)s
_all_results = []
for %(loop_var)s in _loop_data:
    %(code)s
    _result = []
    if len(_returns) == 1:
        _result = eval(_returns[0])
    elif len(_returns) > 1:
        _result = []
        for _j in _returns:
            _result.append(eval(_j))
        _result = tuple(_result)
    else:
        _result = None
    _all_results.append(_result)
_all_results = tuple(_all_results)
"""

def build_loop_code(code,loop_var,inputs,returns,global_vars):
    code_entries ={}
    code_entries['global_code'] = '' # build_globals(global_vars)
    code_entries['input_code'] = '' # build_inputs(inputs)
    code_entries['loop_var'] = loop_var
    code_entries['code'] = indent(code)
    exec_code = loop_code_template % code_entries
    return exec_code
#------------------------------------------------
server_pid = 0
def force_kill():
    """* this is silliness, but I couldn't figure out
        a more elegant approach to killing the server
        process. calling os.exit() in the handler
        only seemed to kill the handler thread -
        not the listening thread
    *"""
    global server_pid
    import os
    print server_pid
    if os.name == 'nt':
        import win32api
        ph = win32api.OpenProcess(1,0,server_pid)
        win32api.TerminateProcess(ph,0)
    else:
        os.kill(server_pid,15) # 15 = TERM 9 = ABRT

class MyThreadingTCPServer(SocketServer.ThreadingTCPServer):
    """ Threaded Server for handling request for python commands
    
        This class was added as of Python2.1 in response
        to the bug reported at:
           
           http://sourceforge.net/tracker/index.php?
              func=detail&aid=417845&group_id=5470&atid=105470
        
        Hopefully the issue will get fixed in the standard library
        by 2.2, and this can go away.
    """
    def handle_request(self):
        """Handle one request, possibly blocking."""
        try:
            request, client_address = self.get_request()
        except socket.error:
            return
        if self.verify_request(request, client_address):
            try:
                self.process_request(request, client_address)
            except:
                self.handle_error(request, client_address)
                self.close_request(request)

    def finish_request(self, request, client_address):
        """Finish one request by instantiating RequestHandlerClass."""
        print 'finish request:', request
        self.RequestHandlerClass(request, client_address, self)
        self.close_request(request)

default_host = socket.gethostname()
def server(host=default_host,port=10000):
    import os
    global server_pid
    server_pid = os.getpid()
    sync_cluster.server_pid = server_pid
    print "starting server on %s:%s" % (host,port)
    print server_pid
    #the_server=SocketServer.TCPServer( (host, port), standard_sync_handler)
    #the_server=SocketServer.ForkingTCPServer( (host, port), standard_sync_handler)
    the_server=MyThreadingTCPServer( (host, port), standard_sync_handler)
    __name__ = '__main__'
    the_server.serve_forever()

def client(host,port,length=1e5):
    import time
    t1 = time.time()
    machine = standard_sync_client(host,port)
    print 'sock open time: ', time.time() - t1

    func = sync_cluster.test_function
    args = (1,2)
    t1 = time.time()
    addendum = {'b':4}
    task = machine.apply_pack(func,args)
    machine.send(task,addendum)
    results = machine.recv()
    t2 = time.time()

#   print 'sent: (beginning)', msg[:5]
#   print 'recv: (beginning)',response[:5]
    print 'roundtrip time:', t2-t1, 'sec'
    #print t2 - t1, (8*2*float(len(msg))/(t2-t1)    / 1e6)

usage = """
sync_cluster mode port host length
mode: server, daemon, client, test
port: server or daemon - port to start server on
      client or test - port to connect to
host - client or test - machine to connect to
length = client - length of message to send
"""

host = default_host
port = 10000

if __name__ == '__main__':
    import sys
    import string
    host = default_host
    if len(sys.argv) >= 5: length = string.atoi(sys.argv[4])
    if len(sys.argv) >= 4: host = sys.argv[3]
    if len(sys.argv) >= 3: port = string.atoi(sys.argv[2])
    if len(sys.argv) >= 2:
        if sys.argv[1] == 'server':
            server(host,port)
        elif sys.argv[1] == 'daemon':
            import os, sync_cluster
            pid = os.fork()
            if pid==0:
                print sync_cluster.__file__, port
                # Child process.
                # This must never return, hence os._exit()!
                try:
                    os.setsid()
                    server(host,port)
                    print "yikes - shouldn't get here"
                finally:
                    os._exit(1)
            else:
                #parent process
                print 'parent'
                os._exit(0)

        elif sys.argv[1] == 'client':
            client(host,port)
        elif sys.argv[1] == 'test':
            import sync_cluster_test
            sync_cluster_test.test(host,port)
        else: print usage
    else:
        print usage
    print 'done'
