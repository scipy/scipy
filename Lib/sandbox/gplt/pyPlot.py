import os, sys,time, tempfile
from types import *
from Numeric import *
from scipy_base.fastumath import *

try:
    import win32pipe
    plot_popen = win32pipe.popen
    default_terminal = 'windows'
except:
    import os
    plot_popen = os.popen
    default_terminal = 'x11'
    #default_terminal = 'png color'

DEBUG = 0

def is1D(a):
    as = shape(a)
    if(len(as) == 1):
        return 1
    if(as[0] == 1 or as[1]==1):
        return 1
    return 0
    
def row(a):
    return reshape(asarray(a),[1,-1])
def col(a):
    return reshape(asarray(a),[-1,1])

SizeMismatch = 'SizeMismatch'
SizeError = 'SizeError'
NotImplemented = 'NotImplemented'

#b = array((0.,0.,0.))
#b1= array((0.,1.,-1))
#bad = b1/b     
#INF = bad[1]
#NEG_INF = bad[2]
b = array((0.,0.))
b1= array((1.,-1))
bad = b1/b      
INF = bad[0]
NEG_INF = bad[1]
BIG = 1e20
SMALL = 1e-20
tmp_x=0


class Plot:
    m_rmax = SMALL  # used to set up polar plots
    m_rmin = BIG    # used to set up polar plots
    
    m_texthandle = 1
    # I have a little C utility program that starts wgnuplot on win32 machines
    # and redirects stdin to the command window of the program.  While the path 
    # of the program can change, I think it will only work if the name of the file
    # is wgnuplot because of how the program goes about finding the command window.
    # Use 'gnuplot' to specify the path of your wgnuplot program.  Also, to help
    # with debugging, you can specify the '-v' arg so that the wgnuplot window is 
    # visible. Otherwise the window is hidden. 
    if sys.platform == 'win32':
        d,f = os.path.split(__file__)
        gnuhelper = os.path.join(d,'gnuplot_helper.exe')
        gnuplot = os.path.join(d,'wgnuplot.exe')
    else:
        gnuhelper = ''
        gnuplot = 'gnuplot'    
###Eric's configurations
###GUI gnuplot
#   gnuhelper = 'c:\programs\gnuplot\gnuplot.exe'
#   gnuplot = 'c:\programs\gnuplot\wgnuplot.exe'
###console gnuplot  
#   gnuhelper = ''
#   gnuplot = 'c:\programs\gnuplot\gnuplot2.exe'
#unix???
#localhost' # doesn't need a hlper function
#   gnuplot = 'gnuplot'
###David's configurations
# I think you have the console version
### GUI gnuplot
#   gnuhelper = 'c:\gnuplot\gnuplot.exe'
#   gnuplot = 'c:\gnuplot\wgnuplot.exe'
### console gnuplot
#   gnuhelper = ''
#   gnuplot = 'c:\gnuplot\gnuplot.exe'

    def __init__(self, *cmd):
        visible = ''
        if (len(cmd) == 1 and type(cmd[0]) == StringType):
            if(cmd[0][0] == 'v' or cmd[0][0] == 'V' or DEBUG == 1):
                visible = ' -v'
        helper = ''
        if self.gnuhelper != '':
            helper = self.gnuhelper + visible + ' ' 
#       self.g = plot_popen(helper + self.gnuplot, 'w')
#unix???    borrowed from Konrad Hinsen's gnuplot interface
#        self.g = plot_popen(helper + self.gnuplot  + ' -persist 1> /dev/null 2>&1', 'w')
        if sys.platform == 'win32':
            self.g = plot_popen(helper + ' ' + self.gnuplot , 'w')
        else: 
            self.g = plot_popen(helper + self.gnuplot  + ' -persist 1> /dev/null 2>&1', 'w')
        self._defaults()
        time.sleep(.2) # put a pause here because in windows, wgnuplot doesn't seem to get the first plot message
        if (len(cmd) > 0 and type(cmd[0]) != StringType):
            self.plot(cmd)
        self.m_tmpfiles = []    
        self.m_label_coord = 'first'
        self.m_angle = 'radians'
        self.m_mapping = 'cartesian'
        self.m_hold = 0
    def valid(self): 
        try:
            # no idea why we need to send so much garbage. 
            self.g.write('  \n') #test the pipe
            self.g.flush()
            self._send('   \n')
            self._send('   \n')
        except (IOError, ValueError): 
            return 0 
        return 1
    def plot_func(self,func):
        self._init_plot()
        if(type(func) == StringType):
            self._send('plot '+ func)
        else:
            print 'error: requires a function as a string'
        
            
#same as plot, but makes a polar plot 
    def plot(self,*data):
        self._init_plot()
        self._plot(data)    

    def polar(self,*data):
        self._init_plot()
        self._plot(data)    
        #m_rmax is calcalculated in _prepare_plot
        #so this has to be done after it is called
        #!requires a replot - might restructure
        # _plot code here so we don't have to plot
        # twice here
        self._init_polar()

    def _init_polar(self):
        cmd = 'set size square;'
        cmd = cmd + 'set label "90" at graph .5, .98 center;'
        cmd = cmd + 'set label "270" at graph .5, .02 center;'
        cmd = cmd + 'set label "0" at graph .98, .5 center;'
        cmd = cmd + 'set label "180" at graph 0.02, 0.5 center;'
        cmd = cmd + 'set polar;'
        if(self.m_angle == 'radians'):
            ttic = pi/180 * 15
        else:
            ttic = 15       
        cmd = cmd + 'set grid polar %f;' % ttic
        cmd = cmd + 'set noborder;'     
        #need more itellignet choice for tic marks
        #shouldn't gnuplot handle this a bit better??
        mn,mx,inc = self._figure_rtics(self.m_rmin,self.m_rmax,5)
        psize = mx * 1.1
        cmd = cmd + 'set xrange [%f:%f];' % (-psize,psize)
        cmd = cmd + 'set yrange [%f:%f];' % (-psize,psize)
        cmd = cmd + 'set rrange [%f:%f];' % (mn,mx)
        cmd = cmd + 'set xtics axis nomirror;'
        cmd = cmd + 'set xtics %f,%f,%f;' % (-mx,inc,mx)
        cmd = cmd + 'set noytics;'
#       cmd = cmd + 'set ytics axis nomirror;'
#       cmd = cmd + 'set ytics %f,%f,%f;' % (mn,inc,mx)
        cmd = cmd + 'set clip;'
#       self._send(cmd)
        self._replot(cmd)

    def _plot(self,data):
        remains = data
        counter = 0
        plotcmd = ''
        tmpf = []
        while len(remains) > 0:
            group, remains = self._get_plot_group(remains)
            if(len(group) > 0):
                try:
                    cmd, file = self._prepare_plot(group)
                except SizeMismatch, data:                  
                    print "Warning: ignoring plot group ", counter
                    print "The arrays should have the same number of " + data[0]
                    print 'array sizes: ', data[1],'and', data[2]
                else:
                    plotcmd = plotcmd + cmd
                    tmpf.append(file)                   
            counter = counter + 1
            
        # if the plot command was valid, and we aren't
        # adding to the current plot, delete all the old plot files 
        if(len(plotcmd) > 0):
            if(not self.m_hold):            
                self._delete_files()
            plotcmd = plotcmd[0:-2] # get rid of trailing comma
            if DEBUG == 1:
                print plotcmd
            if(self.m_hold == 1 and len(self.m_tmpfiles) > 0):
                self._send('replot ' + plotcmd)         
            else:               
                self._send('plot ' + plotcmd)           
            self.m_tmpfiles = self.m_tmpfiles + tmpf
        
    def _prepare_plot(self, group):
        plotinfo = ''
        x = group[0]
        ar_num = 1
        if(len(group) > 1):
            if( type(group[1]) == StringType):
                plotinfo = group[1]
            else:
                ar_num = 2
                y = group[1]                
        if(len(group) == 3):
            plotinfo = group[2]
        
#force 1D arrays to 2D columns
        if(is1D(x) ):
            x = col(x)
        if(ar_num == 2 and is1D(y) ):
            y = col(y)

        xs = shape(x)
        if(ar_num == 2):
            ys = shape(y)   
#test that x and y have compatible shapes
        if(ar_num == 2):
            #check that each array has the same number of rows
            if(xs[0] != ys[0] ):
                raise SizeMismatch, ('rows', xs, ys)
            #check that x.cols = y.cols
            #no error x has 1 column
            if(xs[1] > 1 and xs[1] != ys[1]):
                raise SizeMismatch, ('cols', xs, ys)
#now write the data to a temp file and then plot it         
        filename = tempfile.mktemp()
        f = open(filename, 'w')
        using = ()
        #plot x against index
        if(ar_num == 1):
            mn,mx = self._find_min_max(x)
            if(mn < self.m_rmin):
                self.m_rmin = mn
            if(mx > self.m_rmax):
                self.m_rmax = mx        
            for i in range(0,xs[0]):
                for j in range(0,xs[1] ):
                    f.write( self._format( (x[i][j],) ))
                f.write('\n')           
            for j in range(0,xs[1] ):                       
                using = using + (' using 0:%d ' % (j+1),)
        #plot x vs y                    
        elif(ar_num ==2):
            mn,mx = self._find_min_max(y)
            if(mn < self.m_rmin):
                self.m_rmin = mn
            if(mx > self.m_rmax):
                self.m_rmax = mx
            
            #x is effectively 1D
            if(xs[1] == 1):
                for i in range(0,xs[0] ):
#                   f.write( "%g " % x[i][0])
                    f.write( self._format( (x[i][0],) ))
                    for j in range(0,ys[1] ):
#                       f.write( "%g " % y[i][j])
                        f.write( self._format( (y[i][j],) ))
                    f.write('\n')           
                for j in range(0,ys[1] ):
                    using = using + (' using 1:%d ' % (j+2),)
            #x is 2D                    
            else:
                for i in range(0,xs[0] ):
                    for j in range(0,xs[1] ):
#                       f.write( "%g %g " % (x[i][j], y[i][j]))
                        f.write( self._format( (x[i][j],y[i][j]) ) )
                    f.write('\n')           
                for j in range(0,xs[1] ):
                    using = using + (' using %d:%d ' % (j*2+1,j*2+2),)
        if plotinfo == '':
            #defualt to printing with lines 
            plotinfo = " notitle w l "
        cmd = ''
        # gnuplot doesn't seem to like backslashes
        fn =    string.replace(filename,'\\','/')
        for us in using:
            cmd = cmd + ' "' + fn + '" ' + us + plotinfo + ', '
        return cmd, filename    
        
    def _get_plot_group(self,data):
        group = ()
        remains = data
        state = 0
        finished = 0
        while(len(remains) > 0 and not finished):
            el = remains[0]
            if(state == 0):
                el = asarray(el)
                state = 1 
            elif(state == 1):
                if(type(el) == StringType):
                    finished = 1
                else:
                    el = asarray(el)
                state = 2
            elif(state == 2):
                finished = 1
                if(type(el) != StringType):
                    break                   
            group = group + (el,)
            remains = remains[1:]       
        return group, remains           

    def autoscale(self):
        self._replot('set autoscale xy')

    def _log(self,axis,st='on'):
        from string import upper
        s = upper(st)
        if(s == 'ON'):
           self._replot('set logscale %s' % axis)
        elif(s == 'OFF'):
            self._replot('set nologscale %s' % axis)
        else:
            print 'Command ' + s + 'is invalid.\n use "on" or "off."'           
            return

    def logx(self,st='on'):
        self._log('x',st)

    def logy(self,st='on'):
        self._log('y',st)
         
    def _set_range(self,axis,rng):
        if(rng == 'auto'):
            cmd = 'set autoscale ' + axis       
        else:
            cmd = 'set ' + axis + 'range ' + ' [%g:%g] ' % (rng[0],rng[1])      
        self._replot(cmd)   
    def xaxis(self,rng):
        self._set_range('x',rng)    
    def yaxis(self,rng):
        self._set_range('y',rng)    
    def zaxis(self,rng):
        self._set_range('z',rng)    #not sure that auto will work for z-range
    def raxis(self,rng):
        self._set_range('r',rng)    
    def taxis(self,rng):
        self._set_range('t',rng)    

    def grid(self,st):
        from string import upper
        s = upper(st)
        if(s == 'ON'):
            cmd = 'set grid;'
        elif(s == 'OFF'):
            cmd = 'set nogrid;'
        else:
            print 'Command ' + s + 'is invalid.\n use "on" or "off."'           
            return
        self._replot(cmd+'show grid')
    def hold(self,st):
        from string import upper
        s = upper(st)
        if(s == 'ON'):
            self.m_hold = 1
        elif(s == 'OFF'):
            self.m_hold = 0
        else:
            print 'Command ' + s + 'is invalid.\n use "on" or "off."'           
            return
    def title(self,t):
        cmd = 'set title "' + t + '"'
        self._replot(cmd)

    def xtitle(self,t):
        cmd = 'set xlabel "' + t + '"'
        self._replot(cmd)

    def ytitle(self,t):
        cmd = 'set ylabel "' + t + '"'
        self._replot(cmd)
    def ztitle(self,t):
        cmd = 'set zlabel "' + t + '"'
        self._replot(cmd)
        
#   text( (x,y,z) ,text, *justification, *symbol)   
#later maybe let this handle arrays of text
    def text(self,xyz,*t):
        if(len(xyz) ==2 ):
            point = '%g, %g' % xyz
        elif(len(xyz) ==3 ):    
            point = '%g, %g, %g' % xyz
        else:   
            #should raise exception
            print "invalid size for location vector"
            return
            
        just = ''
        font = ''
        txt=t[0]    
        if(len(t) > 1):
            options = ['center','left','right']
            if self._valid_option(t[1],options):
                just = ' ' + t[1] + ' '
        if(len(t) > 2):
            #no check for correctness of font
            font = ' font "' + t[2] + '"'
        cmd = 'set label ' + ('%d' % self.m_texthandle) + ' "' + txt + '" at ' + self.m_label_coord + ' ' + point + ' ' + just + ' ' + font
        self.m_texthandle = self.m_texthandle + 1
        self._replot(cmd)
        return self.m_texthandle - 1
        
    def label_coord(self,c):
        valid = ['first','second','graph','screen']
        if self._valid_option(c,valid):
            self.m_label_coord = c

    def legend(self,cmd=''):
        if cmd == 'hide':
            self._replot('set nokey')
        else:
            self._replot('set key ' + cmd)

    def close(self):
        self._delete_files()
        if not self.g.closed:
            self._send('exit')
            self.g.close()  

#3d plotting stuff      
#   surf(x,y,z) or surf(z)
    def surf(self,*data):
        self._init_plot()
        cmd, file= self._prepare_surf(data)                 
        self._complete_surf(cmd,file)

#   plot3d(x,y,z)   or plot3d(z)
    def plot3d(self,*data):
        self._init_plot()
        cmd, file= self._prepare_plot3d(data)                   
        self._complete_surf(cmd,file)

        
    def _complete_surf(self,plotcmd,file):
        if(len(plotcmd) > 0):
            if(not self.m_hold):            
                self._delete_files()
                                        
            plotcmd = plotcmd[0:-2] # get rid of trailing comma
            if DEBUG == 1:
                print plotcmd
            if(self.m_hold == 1 and len(self.m_tmpfiles) > 0):
                plotcmd = 'replot ' + plotcmd
            else:       
                plotcmd = 'splot ' + plotcmd
                #force view to default view     
                plotcmd = 'set view 60,30,1,1; ' + plotcmd
            self._send(plotcmd)         
            self.m_tmpfiles.append(file)

            
    def _prepare_surf(self, group):
        plotinfo = ''
        if(len(group) == 1):
            z, = group
            ar_num = 1
        elif(len(group) == 2):
            z,plotinfo = group
            ar_num = 1
        elif(len(group) == 3):
            ar_num = 3
            x,y,z = group
        elif(len(group) == 4):
            ar_num = 3
            x,y,z,plotinfo = group
        else:
            print 'error in _prepare_surf: group can only 4 terms'
            return
        z = asarray(z)
        if(ar_num == 3):
            x = asarray(x)
            y = asarray(y)

#test the shapes of the arrays
        if(ar_num == 1 and is1D(z)):
            z = row(z); #force 1D to 2D 
        zs = shape(z)
        if(len(zs) != 2):
            raise SizeError, ('z is not 2 dimensional', zs)     
        if(ar_num == 3):
            xs = shape(x)
            ys = shape(y)
            if(xs == ys == zs):
                case = 2
            else:   
                case = 1
                if(not is1D(x)):
                    raise SizeError, ('x is not 1 dimensional', xs)     
                if(not is1D(y)):
                    raise SizeError, ('y is not 1 dimensional', ys)
                #force x and y arrays to 1 dimension                
                x = ravel(x)
                y = ravel(y)            
                xs = shape(x)[0]
                ys = shape(y)[0]
    
                #check that each arrays have compatible shapes
                #note: this is sorta backward (put also more visually
                #intuitive) x matches with columns, and y matches with rows
                if(xs != zs[1]):
                    raise SizeMismatch, ('x axis', xs, zs)
                if(ys != zs[0]):
                    raise SizeMismatch, ('y axis', ys, zs)
                                
#now write the data to a temp file and then plot it         
        filename = tempfile.mktemp()
        f = open(filename, 'w')
        # gnuplot doesn't seem to like backslashes
        fn =    string.replace(filename,'\\','/')   
        #plot x against index
        if(ar_num == 1):
            fileformat = ' matrix '             
            for i in range(0,zs[0]):
                for j in range(0,zs[1] ):
                    f.write( "%g " % z[i][j])
                f.write('\n')
        elif(ar_num == 3):
            fileformat = ''             
            if case == 1:
                    for i in range(0,xs):
                        for j in range(0,ys):
                            f.write( self._format( (x[i],y[j],z[j][i]) ) )
                            f.write('\n')
                        f.write('\n')                                       
            else:                   
                for i in range(0,xs[0]):
                    for j in range(0,xs[1]):
#                       f.write( "%g %g %g\n" % (x[i][j], y[i][j], z[i][j]) )               
                        f.write( self._format( (x[i][j], y[i][j], z[i][j]) ))
                        f.write('\n')
                    f.write('\n')
        if plotinfo == '':
            #defualt to printing with lines 
            plotinfo = " notitle w l "                      
        using = ''
        if(self.m_mapping == 'spherical' and ar_num == 3):
            using = 'using ($3*sin($1)*cos($2)):($3*sin($1)*sin($2)):($3*cos($1))'
        elif(self.m_mapping == 'cylindrical' and ar_num == 3):
            using = 'using ($3*cos($1)):($3*sin($1)):2'
        cmd = ' "' +  fn + '" ' + fileformat + ' ' + using + ' ' + plotinfo + ', '
        return cmd, filename    

    def _prepare_plot3d(self, group):
        plotinfo = ''
        data = group[0]
        if type(data) == TupleType or type(data) == ListType:
            data = map(lambda x:asarray(x),data) # convert to list
        elif not (type(data) == ListType):
            data = [asarray(data)]          
        if(len(group) == 2):
            plotinfo = group[1]

#test the shapes of all the arrays
        for datum in data:
            if(shape(datum)[1] != 3):
                raise SizeError, ('the array must have 3 columns', xs)      
#sort the arrays in desending order based on how many entries they have
#this lets us put all the data in a single file
        data.sort(lambda x,y: cmp(-shape(x)[0],-shape(y)[0]) )
        max_entries = shape(data[0])[0];
#open the file
        filename = tempfile.mktemp()
        f = open(filename, 'w')
        # gnuplot doesn't seem to like backslashes
        fn =    string.replace(filename,'\\','/')   
#now write the data to a temp file and then plot it         
        for j in range(0,max_entries):
            for i in range(0,len(data)):
                try:
                    d = data[i][j]
#                   f.write( "%g %g %g " % (d[0],d[1],d[2]) )
                    f.write( self._format( (d[0],d[1],d[2]) ))
                except IndexError:
                    #ignore
                    i=i
            f.write('\n')
        using = ()      
        if DEBUG == 1:
            print self.m_mapping    
        for i in range(0,len(data)):
            us = 'using '
            c = (i*3+1,i*3+2,i*3+3)
            if(self.m_mapping == 'spherical' ):
                us = us + ( '($%d*sin($%d)*cos($%d)):($%d*sin($%d)*sin($%d)):($%d*cos($%d)) ' % 
                         (c[2],   c[0],    c[1],  c[2],   c[0],    c[1],  c[2],   c[0]) ) 
            elif(self.m_mapping == 'cylindrical'):
                using = us + ( '($%d*cos($%d)):($%d*sin($%d)):%d ' %
                           (c[2],  c[0], c[2],  c[0], c[1]) )
            else:               
                us = us + (' %d:%d:%d ' % c)            
            using = using + (us,)
        if plotinfo == '':
            #defualt to printing with lines 
            plotinfo = " notitle w l "  
        cmd = ''                        
        for us in using:
            cmd = cmd + ' "' + fn + '" ' + us + plotinfo + ', '                     
        return cmd, filename    

    def mapping(self,val):
        options = ['cartesian', 'spherical', 'cylindrical']
        if self._valid_option(val,options):
            self.m_mapping = val
#       We'll do the mapping via a ''using' statement
#       because I'd rather use polar-spherical coordinates      
#           self._replot('set mapping ' + val)
            
    def angles(self,val):
        options = ['degrees',  'radians']

        if self._valid_option(val,options):
            self._replot('set angles ' + val)
            self.m_angle = val 
        #set the view of 3d plots val = (xrot, *zrot, *scale, *zscale)
    def view(self,val):
        if(val == 'default'):
            self._replot('set view 60,30,1,1')
            return
            
        if (0 <= val[0] and val[0] <= 180):
            cmd = 'set view ' + ("%f " % val[0] )           
        else:
            print 'error: xrot must be between 0 and 180'           
            return
        if(len(val) > 1):
            zrot = val[1]
            if (0 <= zrot and zrot <= 360):
                cmd = cmd + (",%f " % val[1] )          
            else:
                print 'error: zrot must be between 0 and 360'           
                return
        if(len(val) > 2):
            cmd = cmd + (",%f " % val[2] )          
        if(len(val) > 3):
            cmd = cmd + (",%f " % val[3] )          
        self._replot(cmd)
    # location as percentage of z axis above x,y plane where
    # plotting starts       
    def ticlevel(self,val):
        cmd = 'set ticslevel ' + "%f" % val + ';show tics'
        self._replot(cmd)

    def hidden(self,val):
        options = ['remove',  'show']

        if self._valid_option(val,options):
            no=' o'
            if(val == 'remove'):
                no = '';
            self._replot('set ' + no + 'hidden3d')
    
    def output(self,filename,type,options = ''):
        #changed 7/19/99 ej: always switches to current working directory
        import os, time
        self._send("cd '%s'" % os.path.abspath(os.curdir))
        self._replot('set terminal ' + type + ' ' + options + ";set output '" + filename + "'")
        self._replot('set terminal ' + default_terminal + ';set output')
    def size(self,size=(1.,1.)):
        cmd = 'set size %f,%f' % size
        self._replot(cmd)
    def png(self,filename):
        self.output(filename,'png color')
#contour plotting stuff
#   def contour(self,val):      
                
#utilities
    def _defaults(self):
        #this changes things for the windows terminals
        #just re-arrange some of the colors
        #??? how do I force the defaults to these values
        #??? i.e. can I tell the terminal what I want the 
        #??? default line colors to be?
        cmd ='set linestyle 1 lt 1;'    #bright blue
        cmd = cmd +'set linestyle 2 lt 10;' #dark green
        cmd = cmd +'set linestyle 3 lt 3;'  #red
        cmd = cmd +'set linestyle 4 lt 5;'  #magenta
        cmd = cmd +'set linestyle 5 lt 12;'     #purple
        cmd = cmd +'set linestyle 6 lt 14;'     #cyan
        cmd = cmd +'set linestyle 7 lt 15;'     #yellow
        cmd = cmd +'set linestyle 8 lt 5;'  #dark blue
        cmd = cmd +'set grid;'
        cmd = cmd +'show grid\n'     
        self._send(cmd)
    def _init_plot(self):
        self.m_rmin = BIG
        self.m_rmax = SMALL

        self._send('reset')
        self.grid('on')
        self.angles(self.m_angle)
        
    def _valid_option(self,opt,optlist):
        try:
            optlist.index(opt)
        except ValueError:
            print opt, ' is an invalid choice.'
            print 'Valid choices are ', optlist
            return 0
        else:
            return 1

#used to format data sent to gnuplot files.
#would not be necessary if INF, -INF, and IND
#didn't exist
    def _format(self,x):
        #############################
        ## complications added to check for INF, -INF
        #############################
        #note: does not check for IND at the moment
        str = ''
        for val in x:
            if(val == INF or val == NEG_INF):
                str = str + ' - '
            else:
                str = '%s %g' % (str,val)
        return str
    def _find_min_max(self,x):
        global tmp_x
        #find min...
        temp_min = BIG 
        mn = minimum.reduce(minimum.reduce(x))
        if ( mn < temp_min):
            temp_min = mn

        temp_max = SMALL    
        mx = maximum.reduce(maximum.reduce(x))
        if ( mx > temp_max):
            temp_max = mx
        #############################
        ## complications added to check for INF, -INF
        #############################
        #the following garbage is to check for an infinite or
        #undefined value... there has got to be a smarter way
        # don't compare to 1#IND yet cause it doesn't compare correctly
        if(temp_min == NEG_INF):
            mask = equal(x,NEG_INF)
            tmp_x = choose(mask,(x,BIG))
            mn = minimum.reduce(minimum.reduce(tmp_x))
            if ( mn < temp_min):
                temp_min = mn
        if(temp_max == INF):
            mask = equal(x,INF)
            tmp_x = choose(mask,(x,SMALL))  
            mx = maximum.reduce(maximum.reduce(tmp_x))
            if ( mx < temp_max):
                temp_max = mx
        return (temp_min, temp_max)
    
    #translated from gnuplot code
    def _figure_rtics(self,minval,maxval,guess=20):
        xr = abs(maxval-minval)
        l10 = log10(xr)
        fl = floor(l10)
        xnorm = 10**(l10-fl)
        #you can change the value of 5
        # to something if you want more tics...
        posns = guess / xnorm
        
        if (posns > 40):
            tics = 0.05     # eg 0, .05, .10, ... 
        if (posns > 20):
            tics = 0.1      # eg 0, .1, .2, ... 
        elif (posns > 10):
            tics = 0.2      #/* eg 0,0.2,0.4,... */
        elif (posns > 4):
            tics = 0.5      #/* 0,0.5,1, */
        elif (posns > 1):
            tics = 1        #/* 0,1,2,.... */
        elif (posns > 0.5):
            tics = 2        #/* 0, 2, 4, 6 */
        else:
            # getting desperate... the ceil is to make sure we
            # go over rather than under - eg plot [-10:10] x*x
            # gives a range of about 99.999 - tics=xnorm gives
            # tics at 0, 99.99 and 109.98  - BAD !
            # This way, inaccuracy the other way will round
            # up (eg 0->100.0001 => tics at 0 and 101
            # I think latter is better than former   
            tics = ceil(xnorm);

        tic = tics * 10.0**fl;
        results = [0,0,0]
        results[0] = tic * floor(minval / tic); 
        results[1] = tic * ceil(maxval / tic);  
        results[2] = tic
        return results
        
    def _delete_files(self):
        for file in self.m_tmpfiles:
            try:
                os.unlink(file)
            except:
                #do nothing
                file = file
            else:
                self.m_tmpfiles.remove(file)
    
    def _replot(self,cmd):
        self._send(cmd + ';replot\n')

#this is a compromise: for small commands, send directly to gnuplot
# for large commands, put them in a file and load the file from gnuplot
    def _send(self,cmd):
        if(len(cmd) < 200):
            self.g.write(cmd + '\n')
            self.g.flush()
        else:           
            filename = tempfile.mktemp()
            f = open(filename, 'w')
            f.write(cmd + '\n')
            f.close();
            fn =    string.replace(filename,'\\','/')
            self.g.write('load "' + fn +  '"\n')
            self.g.flush()  
            self.m_tmpfiles.append(filename)
        time.sleep(.15)             
#I'd rather use some function similar to the following to send data to gnuplot
#but I'm having problems with getting all the data transferred to gnuplot in 
#large(huge) plot commands - it hangs after 2042 characters.  I don't know
#if the problem is with the pipe, with gnuplots command line buffer size, or
#my helper program (probable??)     
#NOTE: Using spy++ to watch WM_CHAR messages, I was able to confirm that 
# all the characters are indeed sent to wgnuplot by the gnuhelper program.
# It looks as if wgnuplot message pump just can't handle it if a barrage 
# of events are sent its way.  If I put in huge delays, the problem goes away,
# but for now, I think it is more efficient just to use files.
    def _send_try2_but_still_doesnt_work(self,cmd):
        bufsize = 50
        idx = range(0,len(cmd),bufsize)
        idx.append(len(cmd))
        for i in range(len(idx)-1):         
#           print i, idx[i]+1, idx[i+1]
            print cmd[idx[i]:idx[i+1]]
            self.g.write(cmd[idx[i]:idx[i+1]])
            self.g.flush()  
            time.sleep(.05)         
        self.g.write('\n')
        self.g.flush()      
        #without the sleep command, it appears things get lost in the pipe!
        #I don't understand this
        time.sleep(.15)   
#and really it should be this easy...
    def _send_best_but_doesnt_work(self,cmd):
        self.g.write(cmd + '\n')
        self.g.flush()  
#destructor
    def __del__(self):
        try: self.close()
        except (ValueError, IOError, AttributeError): pass

#this version of prepare surf does all the coordinate mapping in python.
#Initially I thought this was a good idea, but now I have no idea why I thought
#that.  I'll keep it until I'm sure it isn't a good idea
    def _prepare_surf_slower(self, group):
        plotinfo = ''
        if(len(group) == 1):
            z, = group
            ar_num = 1
        elif(len(group) == 2):
            z,plotinfo = group
            ar_num = 1
        elif(len(group) == 3):
            ar_num = 3
            x,y,z = group
        elif(len(group) == 4):
            ar_num = 3
            x,y,z,plotinfo = group
        else:
            print 'error in _prepare_surf: group can only 4 terms'
            return
        z = asarray(z)
        if(ar_num == 3):
            x = asarray(x)
            y = asarray(y)

#test the shapes of the arrays
        if(ar_num == 1 and is1D(z)):
            z = row(z); #force 1D to 2D 
        zs = shape(z)
        if(len(zs) != 2):
            raise SizeError, ('z is not 2 dimensional', zs)     
        if(ar_num == 3):
            xs = shape(x)
            ys = shape(y)
            if(xs == ys == zs):
                case = 2
            else:   
                case = 1
                if(not is1D(x)):
                    raise SizeError, ('x is not 1 dimensional', xs)     
                if(not is1D(y)):
                    raise SizeError, ('y is not 1 dimensional', ys)
                #force x and y arrays to 1 dimension                
                x = ravel(x)
                y = ravel(y)            
                xs = shape(x)[0]
                ys = shape(y)[0]
    
                #check that each arrays have compatible shapes
                #note: this is sorta backward (put also more visually
                #intuitive) x matches with columns, and y matches with rows
                if(xs != zs[1]):
                    raise SizeMismatch, ('x axis', xs, zs)
                if(ys != zs[0]):
                    raise SizeMismatch, ('y axis', ys, zs)
                                
#now write the data to a temp file and then plot it         
        filename = tempfile.mktemp()
        f = open(filename, 'w')
        # gnuplot doesn't seem to like backslashes
        fn =    string.replace(filename,'\\','/')   
        #plot x against index
        if(ar_num == 1):
            fileformat = ' matrix '             
            for i in range(0,zs[0]):
                for j in range(0,zs[1] ):
                    f.write( "%g " % z[i][j])
                f.write('\n')
        elif(ar_num == 3):
            fileformat = ''             
            if case == 1:
                if(self.m_mapping == 'spherical'):
                    sinx = sin(x)
                    cosx = cos(x)       
                    siny = sin(y)
                    cosy = cos(y)
                    for i in range(0,xs):
                        for j in range(0,ys):
                            v1 = z[j][i] * sinx[i]*cosy[j]
                            v2 = z[j][i] * sinx[i]*siny[j]
                            v3 = z[j][i] * cosx[i]
                            f.write( self._format( (v1,v2,v3) ) )
                            f.write('\n')
                        f.write('\n')                       
                elif(self.m_mapping == 'cylindrical'):
                    sinx = sin(x)
                    cosx = cos(x)       
                    for i in range(0,xs):
                        for j in range(0,ys):
                            v1 = z[j][i] * cosx[i]
                            v2 = z[j][i] * sinx[i]
                            v3 = y[j]
                            f.write( self._format( (v1,v2,v3) ) )
                            f.write('\n')
                        f.write('\n')                       
                else:
                    for i in range(0,xs):
                        for j in range(0,ys):
                            v1 = x[i]
                            v2 = y[j]
                            v3 = z[j][i]
                            f.write( self._format( (v1,v2,v3) ) )
                            f.write('\n')
                        f.write('\n')                       
                        
            else:                   
                for i in range(0,xs[0]):
                    for j in range(0,xs[1]):
#                       f.write( "%g %g %g\n" % (x[i][j], y[i][j], z[i][j]) )               
                        f.write( self._format( (x[i][j], y[i][j], z[i][j]) ))
                        f.write('\n')
                    f.write('\n')
        if plotinfo == '':
            #defualt to printing with lines 
            plotinfo = " notitle w l "                      
        using = ''
#       if(self.m_mapping == 'spherical' and ar_num == 3):
#           using = 'using ($3*sin($1)*cos($2)):($3*sin($1)*sin($2)):($3*cos($1))'
#       elif(self.m_mapping == 'cylindrical' and ar_num == 3):
#           using = 'using ($3*cos($1)):($3*sin($1)):2'
        cmd = ' "' +  fn + '" ' + fileformat + ' ' + using + ' ' + plotinfo + ', '
        return cmd, filename    
