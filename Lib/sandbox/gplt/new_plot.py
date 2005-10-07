#this just adds errorbar plotting capabilities
from Numeric import *
from scipy_base.scimath import *
from types import StringType
import pyPlot
from pyPlot import *

def squeeze(a):
	s = shape(a)
	dims = len(s) - sum(equal(1,s))
	new_s = compress(not_equal(1,s),s)
	return reshape(a,new_s)


class Plot(pyPlot.Plot):
	errorbar_lt = 1
	def errorbar(self,*data):
		self._init_plot()
		self._errorbar_plot(data)	

	def _errorbar_plot(self,data):
		remains = data
		counter = 0
		plotcmd = ''
		tmpf = []
		if not self.m_hold: self.errorbar_lt = 1
		while len(remains) > 0:
			group, remains = self._errorbar_get_plot_group(remains)
			if(len(group) > 0):
				try:
					cmd, file = self._errorbar_prepare_plot(group)
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
		
	def _errorbar_prepare_plot(self, group):
		plotinfo = ''
		x = group[0]
		y = group[1]				
		e = group[2]
		ar_num = 3
		if(len(group) > 3):
			plotinfo = group[3]
		
#force 1D arrays to 2D columns
		if(is1D(x) ): x = col(x)
		if(is1D(y) ): y = col(y)
		if(is1D(e) ): e = col(e)
		xs = shape(x)
		ys = shape(y)	
		es = shape(e)
#test that x and y have compatible shapes
		#check that each array has the same number of rows
		if(xs[0] != ys[0] or xs[0] != es[0]):
			raise SizeMismatch, ('rows', xs, ys, es)
		#check that x.cols = y.cols
		#no error x has 1 column
		if(xs[1] > 1 and (xs[1] != ys[1] or xs[1] !=es[1])):
				raise SizeMismatch, ('cols', xs, ys, es)
		if(ys[1] != es[1]):
			raise SizeMismatch, ('cols', xs, ys, es)		
#now write the data to a temp file and then plot it			
		filename = tempfile.mktemp()
		f = open(filename, 'w')
		using = ()

		mn,dummy = self._find_min_max(y-e)
		dummy,mx = self._find_min_max(y+e)
		if(mn < self.m_rmin):
			self.m_rmin = mn
		if(mx > self.m_rmax):
			self.m_rmax = mx
		
		#x is effectively 1D
		if(xs[1] == 1):
			for i in range(0,xs[0] ):
				f.write( self._format( (x[i][0],) ))
				for j in range(0,ys[1] ):
					f.write( self._format( (y[i][j],e[i][j]) ))
				f.write('\n')			
			for j in range(0,ys[1] ):
				using = using + (' using 1:%d:%d' % (j+2,j+3),)
			#x is 2D					
		else:
			for i in range(0,xs[0] ):
				for j in range(0,xs[1] ):
					f.write( self._format( (x[i][j],y[i][j],e[i][j]) ) )
				f.write('\n')			
			for j in range(0,xs[1] ):
				using = using + (' using %d:%d:%d' % (j*2+1,j*2+2,j*2+3),)
		if plotinfo == '':
			#defualt to printing with lines	
			plotinfo = " notitle w lp  lt " +`self.errorbar_lt` + ' '
		cmd = ''
		# gnuplot doesn't seem to like backslashes
		fn = 	string.replace(filename,'\\','/')
		for us in using:
			cmd = cmd + ' "' + fn + '" ' + us + 'notitle w errorbars lt ' +`self.errorbar_lt` + ', '
			cmd = cmd + ' "' + fn + '" ' + us + plotinfo + ', '
		self.errorbar_lt = (self.errorbar_lt + 1) % 15
		return cmd, filename	
		
	def _errorbar_get_plot_group(self,data):
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
				el = asarray(el)
				state = 2
			elif(state == 2):
				if(type(el) == StringType):
					finished = 1
				else:
					el = asarray(el)
					state = 3
			elif(state == 3):
				finished = 1
				if(type(el) != StringType):
					break					
			group = group + (el,)
			remains = remains[1:]		
		return group, remains	
